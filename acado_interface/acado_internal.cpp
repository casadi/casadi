/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A minimalistic computer algebra system with automatic differentiation 
 *    and framework for dynamic optimization.
 *    Copyright (C) 2010 by Joel Andersson, Moritz Diehl et al., K.U.Leuven. All rights reserved.
 *
 *    CasADi is free software; you can redistribute it and/or
 *    modify it under the terms of the GNU Lesser General Public
 *    License as published by the Free Software Foundation; either
 *    version 3 of the License, or (at your option) any later version.
 *
 *    CasADi is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *    Lesser General Public License for more details.
 *
 *    You should have received a copy of the GNU Lesser General Public
 *    License along with CasADi; if not, write to the Free Software
 *    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 *
 */

#include "acado_internal.hpp"
#include "acado_integrator_backend.hpp"

#include <acado_optimal_control.hpp>
#include <casadi/stl_vector_tools.hpp>

#include <cassert>
#include <limits>

using namespace std;

namespace CasADi{

AcadoInternal::AcadoInternal(const FX& ffcn, const FX& mfcn, const FX& cfcn, const FX& rfcn) : ffcn_(ffcn), mfcn_(mfcn), cfcn_(cfcn), rfcn_(rfcn){
  // Options
  addOption("integrator_tolerance",     OT_REAL);
  addOption("absolute_tolerance",       OT_REAL);
  addOption("kkt_tolerance",            OT_REAL);
  addOption("number_of_shooting_nodes", OT_INTEGER,  20);
  addOption("max_num_iterations",       OT_INTEGER);
  addOption("hessian_approximation",    OT_STRING);
  addOption("dynamic_sensitivity",      OT_STRING); // forward_sensitivities or backward_sensitivities
  addOption("exact_jacobian",           OT_BOOLEAN,  true);
  addOption("start_time",               OT_REAL,  0.0);
  addOption("final_time",               OT_REAL,  1.0);
  addOption("print_level",              OT_STRING,  "low"); // "none", "low", "medium", "high", "debug"
  addOption("auto_init",                OT_BOOLEAN,  false); // initialize differential and angebraic states by a forward integration
  addOption("max_num_integrator_steps", OT_INTEGER);
  addOption("relaxation_parameter",     OT_REAL);
  addOption("periodic_bounds",          OT_INTEGERVECTOR);
  addOption("integrator",               OT_STRING);
  
  
// Set pointers to null
  t_=0;
  xd_=0;
  xa_=0;
  u_=0;
  p_=0;
  xdot_=0;
  f_=0;
  arg_ = 0;
  ocp_ = 0;
  algorithm_ = 0;
  
  // The following will make sure that all counters are set to zero
  ACADO::AlgebraicState().clearStaticCounters();
  ACADO::Control().clearStaticCounters();
  ACADO::DifferentialState().clearStaticCounters();
  ACADO::DifferentialStateDerivative().clearStaticCounters();
  ACADO::Disturbance().clearStaticCounters();
  ACADO::IntegerControl().clearStaticCounters();
  ACADO::IntegerParameter().clearStaticCounters();
  ACADO::IntermediateState().clearStaticCounters();
  ACADO::Parameter().clearStaticCounters();
  
#ifdef ACADO_HAS_USERDEF_INTEGRATOR
  // Point the pointer in the ACADO integrator class to the AcadoIntegratorBackend create function
  if(ACADO::Integrator::integrator_creator_ || ACADO::Integrator::integrator_user_data_)
    throw CasadiException("AcadoInternal::AcadoInternal: An instance already exists");
#endif
  
}

AcadoInternal::~AcadoInternal(){

  // Free memory
  if(t_) delete t_;
  if(xd_) delete[] xd_;
  if(xa_) delete[] xa_;
  if(u_) delete[] u_;
  if(p_) delete[] p_;
  if(xdot_) delete[] xdot_;
  
  if(f_) delete f_;
  if(ocp_) delete ocp_;
  if(algorithm_) delete algorithm_;
  if(arg_) delete arg_;

  ACADO::Integrator::integrator_creator_ = 0;
  ACADO::Integrator::integrator_user_data_ = 0;
}


void AcadoInternal::init(){
  // Initialize the functions and get dimensions
  // DAE
  ffcn_.init();
  
  // Get dimensions
  nt_ = ffcn_.f_.input(ACADO_FCN_T).numel();
  nxd_ = ffcn_.f_.input(ACADO_FCN_XD).numel();
  nxa_ = ffcn_.f_.input(ACADO_FCN_XA).numel();
  nx_ = nxd_ + nxa_;
  nu_  = ffcn_.f_.input(ACADO_FCN_U).numel();
  np_  = ffcn_.f_.input(ACADO_FCN_P).numel();
  nxdot_  = ffcn_.f_.input(ACADO_FCN_XDOT).numel();

  // Objective
  mfcn_.init();

  // Path constraints
  if(!cfcn_.f_.isNull()){
    cfcn_.init();
    nc_ = cfcn_.f_.output().numel();
  } else {
    nc_ = 0;
  }
  
  // Initial constraint
  if(!rfcn_.f_.isNull()){
    rfcn_.init();
    nr_ = rfcn_.f_.output().numel();
  } else {
    nr_ = 0;
  }
  
  // Print:
  cout << "nt = " << nt_ << endl;
  cout << "nxd = " << nxd_ << endl;
  cout << "nxa = " << nxa_ << endl;
  cout << "nxa = " << nxa_ << endl;
  cout << "nu = " << nu_ << endl;
  cout << "np = " << np_ << endl;
  cout << "nxdot = " << nxdot_ << endl;
  cout << "nr = " << nr_ << endl;
  cout << "nc = " << nc_ << endl;

  // Number of shooting nodes
  n_nodes_ = getOption("number_of_shooting_nodes").toInt();

  // Input dimensions
  input_.resize(ACADO_NUM_IN);
  input(ACADO_X_GUESS).setSize(nx_,n_nodes_+1);
  input(ACADO_U_GUESS).setSize(nu_,n_nodes_+1);
  input(ACADO_P_GUESS).setSize(np_);

  input(ACADO_LBX).setSize(nx_);
  input(ACADO_UBX).setSize(nx_);
  input(ACADO_LBX0).setSize(nx_);
  input(ACADO_UBX0).setSize(nx_);
  input(ACADO_LBXF).setSize(nx_);
  input(ACADO_UBXF).setSize(nx_);
  
  input(ACADO_LBU).setSize(nu_);
  input(ACADO_UBU).setSize(nu_);
  
  input(ACADO_LBP).setSize(np_);
  input(ACADO_UBP).setSize(np_);
  
  input(ACADO_LBC).setSize(nc_);
  input(ACADO_UBC).setSize(nc_);
  
  input(ACADO_LBR).setSize(nr_);
  input(ACADO_UBR).setSize(nr_);
  
  // Output dimensions
  output_.resize(ACADO_NUM_OUT);
  output(ACADO_X_OPT).setSize(nx_,n_nodes_+1);
  output(ACADO_U_OPT).setSize(nu_,n_nodes_+1);
  output(ACADO_P_OPT).setSize(np_);
  output(ACADO_COST).setSize(1);

  // Initialize
  FXNode::init();

  // Set all bounds except initial constraints to +- infinity by default
  for(int i=ACADO_LBX; i<ACADO_UBC; i = i+2){
    vector<double> &lb = input(i).data();
    fill(lb.begin(),lb.end(),-numeric_limits<double>::infinity());
    
    vector<double> &ub = input(i+1).data();
    fill(ub.begin(),ub.end(),numeric_limits<double>::infinity());
  }
  
  // variables
  t_ = new ACADO::TIME();
  xd_ = new ACADO::DifferentialState[nxd_];
  xa_ = new ACADO::AlgebraicState[nxa_];
  u_ = new ACADO::Control[nu_];
  p_ = new ACADO::Parameter[np_];
  xdot_ = new ACADO::DifferentialStateDerivative[nxdot_];
  f_ = new ACADO::DifferentialEquation();

  // Augmented state vector
  arg_ = new ACADO::IntermediateState(nt_+nxd_+nxa_+nu_+np_+nxdot_);
  int ind=0;
  for(int i=0; i<nt_; ++i)    (*arg_)(ind++) = t_[i];
  for(int i=0; i<nxd_; ++i)   (*arg_)(ind++) = xd_[i];
  for(int i=0; i<nxa_; ++i)   (*arg_)(ind++) = xa_[i];
  for(int i=0; i<nu_; ++i)    (*arg_)(ind++) = u_[i];
  for(int i=0; i<np_; ++i)    (*arg_)(ind++) = p_[i];
  for(int i=0; i<nxdot_; ++i) (*arg_)(ind++) = xdot_[i];

  // Create an ocp object
  double t_start = getOption("start_time").toDouble();
  double t_final = getOption("final_time").toDouble();
  ocp_ = new ACADO::OCP(t_start, t_final, n_nodes_);

  // Pass objective function
  ocp_->minimizeMayerTerm( (*mfcn_.fcn_)(*arg_) );
  
  // Pass dynamic equation
  ocp_->subjectTo( *f_ << (*ffcn_.fcn_)(*arg_) );
  
}

void AcadoInternal::evaluate(int fsens_order, int asens_order){ 
  // Initial constraint function
  if(!rfcn_.f_.isNull()){
    const vector<double>& lbr = input(ACADO_LBR).data();
    ACADO::Vector lb(lbr.size(),&lbr[0]);
    const vector<double>& ubr = input(ACADO_UBR).data();
    ACADO::Vector ub(ubr.size(),&ubr[0]);
    ocp_->subjectTo( ACADO::AT_START, lb <= (*rfcn_.fcn_)(*arg_) <= ub);
  }

  // Path constraint function
  if(!cfcn_.f_.isNull()){
    const vector<double>& lbc = input(ACADO_LBC).data();
    ACADO::Vector lb(lbc.size(),&lbc[0]);
    const vector<double>& ubc = input(ACADO_UBC).data();
    ACADO::Vector ub(ubc.size(),&ubc[0]);
    ocp_->subjectTo( lb <= (*cfcn_.fcn_)(*arg_) <=  ub );
  }

  // State bounds
  vector<double> &lbx = input(ACADO_LBX).data();
  vector<double> &ubx = input(ACADO_UBX).data();
  for(int i=0; i<nxd_; ++i)
    ocp_->subjectTo( lbx[i] <= xd_[i] <=  ubx[i] );
  for(int i=nxd_; i<nx_; ++i)
    ocp_->subjectTo( lbx[i] <= xa_[i-nxd_] <=  ubx[i] );

  // Pass bounds on state at initial time
  vector<double> &lbx0 = input(ACADO_LBX0).data();
  vector<double> &ubx0 = input(ACADO_UBX0).data();
  for(int i=0; i<nxd_; ++i)
    ocp_->subjectTo( ACADO::AT_START, lbx0[i] <= xd_[i] <=  ubx0[i] );
  for(int i=nxd_; i<nx_; ++i)
    ocp_->subjectTo( ACADO::AT_START, lbx0[i] <= xa_[i-nxd_] <=  ubx0[i] );

//     ocp_->subjectTo( AT_END  , xd_[1] ==  0.0 );
//     ocp_->subjectTo( AT_END  , xd_[2] ==  0.0 );

  // Control bounds
  vector<double> &lbu = input(ACADO_LBU).data();
  vector<double> &ubu = input(ACADO_UBU).data();
  for(int i=0; i<nu_; ++i)
    ocp_->subjectTo( lbu[i] <= u_[i] <= ubu[i] );

  // Parameter bounds
  vector<double> &lbp = input(ACADO_LBP).data();
  vector<double> &ubp = input(ACADO_UBP).data();
  for(int i=0; i<np_; ++i)
    ocp_->subjectTo( lbp[i] <= p_[i] <= ubp[i] );

  // Periodic boundary condition
  if(hasSetOption("periodic_bounds")){
    const vector<int>& periodic = getOption("periodic_bounds").toIntVector();
    if(periodic.size()!=nx_) throw CasadiException("wrong dimension for periodic_bounds");
    for(int i=0; i<nxd_; ++i)
      if(periodic[i])
        ocp_->subjectTo( 0.0, xd_[i], -xd_[i], 0.0);

    for(int i=nxd_; i<nx_; ++i)
      if(periodic[i])
        ocp_->subjectTo( 0.0, xa_[i-nxd_], -xa_[i-nxd_], 0.0);
  }
  
  algorithm_ = new ACADO::OptimizationAlgorithm(*ocp_);
  
  // set print level
  ACADO::PrintLevel printlevel;
  if(getOption("print_level")=="none")        printlevel = ACADO::NONE;
  else if(getOption("print_level")=="low")    printlevel = ACADO::LOW;
  else if(getOption("print_level")=="medium") printlevel = ACADO::MEDIUM;
  else if(getOption("print_level")=="high")   printlevel = ACADO::HIGH;
  else if(getOption("print_level")=="debug")  printlevel = ACADO::DEBUG;
  else throw CasadiException("Illegal print level. Allowed are \"none\", \"low\", \"medium\", \"high\", \"debug\"");
  algorithm_->set(ACADO::INTEGRATOR_PRINTLEVEL, printlevel );

  // Set integrator
  if(hasSetOption("integrator")){
    Option integ = getOption("integrator");
    ACADO::IntegratorType itype;
    if(integ=="rk4")           itype=ACADO::INT_RK4;
    else if(integ=="rk12")     itype=ACADO::INT_RK12;
    else if(integ=="rk23")     itype=ACADO::INT_RK23;
    else if(integ=="rk45")     itype=ACADO::INT_RK45;
    else if(integ=="rk78")     itype=ACADO::INT_RK78;
    else if(integ=="bdf")      itype=ACADO::INT_BDF;
    else if(integ=="discrete") itype=ACADO::INT_DISCRETE;
    else if(integ=="unknown")  itype=ACADO::INT_UNKNOWN;
    #ifdef ACADO_HAS_USERDEF_INTEGRATOR
    else if(integ=="casadi"){
      ACADO::Integrator::integrator_creator_ = &ACADO::AcadoIntegratorBackend::create;
      ACADO::Integrator::integrator_user_data_ = this;
      itype=ACADO::INT_UNKNOWN;
    }
    #endif
  else throw CasadiException("AcadoInternal::evaluate: no such integrator: " + integ.toString());
    algorithm_->set(ACADO::INTEGRATOR_TYPE, itype);
  };
  
  // Set integrator tolerance
  if(hasSetOption("integrator_tolerance")) algorithm_->set( ACADO::INTEGRATOR_TOLERANCE, getOption("integrator_tolerance").toDouble());
  if(hasSetOption("absolute_tolerance")) algorithm_->set( ACADO::ABSOLUTE_TOLERANCE, getOption("absolute_tolerance").toDouble());
  if(hasSetOption("kkt_tolerance")) algorithm_->set( ACADO::KKT_TOLERANCE, getOption("kkt_tolerance").toDouble());
  if(hasSetOption("max_num_iterations")) algorithm_->set( ACADO::MAX_NUM_ITERATIONS, getOption("max_num_iterations").toInt() );
  if(hasSetOption("max_num_integrator_steps")) algorithm_->set( ACADO::MAX_NUM_INTEGRATOR_STEPS, getOption("max_num_integrator_steps").toInt() );
  if(hasSetOption("relaxation_parameter")) algorithm_->set( ACADO::RELAXATION_PARAMETER, getOption("relaxation_parameter").toDouble());
  
  if(hasSetOption("dynamic_sensitivity"))
    if(getOption("dynamic_sensitivity") == "forward_sensitivities")
      algorithm_->set( ACADO::DYNAMIC_SENSITIVITY,  ACADO::FORWARD_SENSITIVITY );
    else if(getOption("dynamic_sensitivity") == "backward_sensitivities")
      algorithm_->set( ACADO::DYNAMIC_SENSITIVITY,  ACADO::BACKWARD_SENSITIVITY );
    else 
      throw CasadiException("illegal dynamic_sensitivity");

  if(hasSetOption("hessian_approximation")){
    int hess;
    Option op = getOption("hessian_approximation");
    if(op=="exact_hessian")                 hess = ACADO::EXACT_HESSIAN;
    else if(op == "constant_hessian")       hess = ACADO::CONSTANT_HESSIAN;
    else if(op == "full_bfgs_update")       hess = ACADO::FULL_BFGS_UPDATE;
    else if(op == "block_bfgs_update")      hess = ACADO::BLOCK_BFGS_UPDATE;
    else if(op == "gauss_newton")           hess = ACADO::GAUSS_NEWTON;
    else if(op == "gauss_newton_with_block_bfgs") hess = ACADO::GAUSS_NEWTON_WITH_BLOCK_BFGS;
    else throw CasadiException("illegal hessian approximation");
    
    algorithm_->set( ACADO::HESSIAN_APPROXIMATION,  hess);
  }

  // should the states be initialized by a forward integration?
  bool auto_init = getOption("auto_init").toInt();

  // Initialize differential states
  if(nxd_>0){
    // Initial guess
    vector<double> &x0 = input(ACADO_X_GUESS).data();
    
    // Assemble the variables grid
    ACADO::VariablesGrid xd(nxd_, n_nodes_+1);
    for(int i=0; i<n_nodes_+1; ++i){
      ACADO::Vector v(nxd_,&x0[i*nx_]);
      xd.setVector(i,v);
    }
    
    // Pass to acado
    algorithm_->initializeDifferentialStates(xd,auto_init ? ACADO::BT_TRUE : ACADO::BT_FALSE);
  }
    
  // Initialize algebraic states
  if(nxa_>0){
    // Initial guess
    vector<double> &x0 = input(ACADO_X_GUESS).data();
    
    // Assemble the variables grid
    ACADO::VariablesGrid xa(nxa_, n_nodes_+1);
    for(int i=0; i<n_nodes_+1; ++i){
      ACADO::Vector v(nxa_,&x0[i*nx_+nxd_]);
      xa.setVector(i,v);
    }
    
    // Pass to acado
    algorithm_->initializeAlgebraicStates(xa,auto_init ? ACADO::BT_TRUE : ACADO::BT_FALSE);
  }
    
  // Initialize controls
  if(nu_>0){
    // Initial guess
    vector<double> &u0 = input(ACADO_U_GUESS).data();
    
    // Assemble the variables grid
    ACADO::VariablesGrid u(nu_, n_nodes_+1);
    for(int i=0; i<n_nodes_+1; ++i){
      ACADO::Vector v(nu_,&u0[i*nu_]);
      u.setVector(i,v);
    }
    
    // Pass to acado
    algorithm_->initializeControls(u);
  }
    
  // Initialize parameters
  if(np_>0){
    // Initial guess
    vector<double> &p0 = input(ACADO_P_GUESS).data();
    
    // Assemble the variables grid
    ACADO::VariablesGrid p(np_, n_nodes_+1);
    for(int i=0; i<n_nodes_+1; ++i){
      ACADO::Vector v(np_,&p0[0]); // NB!
      p.setVector(i,v);
    }
    
    // Pass to acado
    algorithm_->initializeParameters(p);
  }

  // Solve
  algorithm_->solve();

  // Get the optimal state trajectory
  if(nxd_>0){
    vector<double> &xopt = output(ACADO_X_OPT).data();
    ACADO::VariablesGrid xd;
    algorithm_->getDifferentialStates(xd);
    assert(xd.getNumPoints()==n_nodes_+1);
    for(int i=0; i<n_nodes_+1; ++i){
      // Copy to result
      ACADO::Vector v = xd.getVector(i);
      &xopt[i*nx_] << v;
    }
  }
  if(nxa_>0){
    vector<double> &xopt = output(ACADO_X_OPT).data();
    ACADO::VariablesGrid xa;
    algorithm_->getAlgebraicStates(xa);
    assert(xa.getNumPoints()==n_nodes_+1);
    for(int i=0; i<n_nodes_+1; ++i){
      // Copy to result
      ACADO::Vector v = xa.getVector(i);
      &xopt[i*nx_ + nxd_] << v;
    }
  }

  // Get the optimal control trajectory
  if(nu_>0){
    vector<double> &uopt = output(ACADO_U_OPT).data();
    ACADO::VariablesGrid u;
    algorithm_->getControls(u);
    assert(u.getNumPoints()==n_nodes_+1);
    for(int i=0; i<n_nodes_+1; ++i){
      // Copy to result
      ACADO::Vector v = u.getVector(i);
      &uopt[i*nu_] << v;
    }
  }

  // Get the optimal parameters
  if(np_>0){
    vector<double> &popt = output(ACADO_P_OPT).data();
    ACADO::Vector p;
    algorithm_->getParameters(p);
    &popt[0] << p;
  }
  
  // Get the optimal cost
  double cost = algorithm_->getObjectiveValue();
  output(ACADO_COST).set(cost);
}



} // namespace CasADi


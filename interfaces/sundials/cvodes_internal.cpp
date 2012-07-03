/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A symbolic framework for dynamic optimization.
 *    Copyright (C) 2010 by Joel Andersson, Moritz Diehl, K.U.Leuven. All rights reserved.
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

#include "cvodes_internal.hpp"
#include "casadi/stl_vector_tools.hpp"
#include "casadi/sx/sx_tools.hpp"
#include "casadi/fx/linear_solver_internal.hpp"
#include "casadi/fx/mx_function.hpp"

using namespace std;
namespace CasADi{
namespace Sundials{

CVodesInternal* CVodesInternal::clone() const{
  // Return a deep copy
  CVodesInternal* node = new CVodesInternal(f_,g_);
  node->setOption(dictionary());
  node->jac_f_ = jac_f_;
  node->jac_ = jac_;
  node->linsol_ = linsol_;
  return node;
}

CVodesInternal::CVodesInternal(const FX& f, const FX& g) : SundialsInternal(f,g){
  addOption("linear_multistep_method",     OT_STRING,  "bdf","bdf|adams");
  addOption("nonlinear_solver_iteration",  OT_STRING,  "newton","","newton|functional");
  addOption("fsens_all_at_once",           OT_BOOLEAN,true); // calculate all right hand sides of the sensitivity equations at once
  addOption("disable_internal_warnings",   OT_BOOLEAN,false, "Disable CVodes internal warning messages");
  addOption("monitor",      OT_STRINGVECTOR, GenericType(),  "", "res|resB|resQB|reset", true);
    
  mem_ = 0;

  x0_ = x_ = 0;
  q_ = 0;

  is_init = false;
  isInitAdj_ = false;
  disable_internal_warnings_ = false;
}

void CVodesInternal::freeCVodes(){
  if(mem_) CVodeFree(&mem_);

  // ODE integration
  if(x0_) N_VDestroy_Serial(x0_);
  if(x_) N_VDestroy_Serial(x_);
  if(q_) N_VDestroy_Serial(q_);
  
  // Forward problem
  for(vector<N_Vector>::iterator it=xF0_.begin(); it != xF0_.end(); ++it)   if(*it) N_VDestroy_Serial(*it);
  for(vector<N_Vector>::iterator it=xF_.begin(); it != xF_.end(); ++it)     if(*it) N_VDestroy_Serial(*it);
  for(vector<N_Vector>::iterator it=qF_.begin(); it != qF_.end(); ++it)   if(*it) N_VDestroy_Serial(*it);
  
  // Adjoint problem
  for(vector<N_Vector>::iterator it=xA0_.begin(); it != xA0_.end(); ++it)   if(*it) N_VDestroy_Serial(*it);
  for(vector<N_Vector>::iterator it=xA_.begin(); it != xA_.end(); ++it)     if(*it) N_VDestroy_Serial(*it);
  for(vector<N_Vector>::iterator it=qA_.begin(); it != qA_.end(); ++it)   if(*it) N_VDestroy_Serial(*it);
}

CVodesInternal::~CVodesInternal(){
  freeCVodes();
}

void CVodesInternal::updateNumSens(bool recursive){
  // Not supported re-initalization needed
  init();
}

void CVodesInternal::init(){
  // Initialize the base classes
  SundialsInternal::init();

  // Free memory if already initialized
  if(is_init) freeCVodes();
  
  // Read options
  monitor_rhsB_  = monitored("resB");
  monitor_rhs_   = monitored("res");
  monitor_rhsQB_ = monitored("resQB");
  
  // Try to generate a jacobian if none provided
  if(!linsol_.isNull() && jac_.isNull()){
    SXFunction f = shared_cast<SXFunction>(f_);
    if(!f.isNull()){
      // Get the Jacobian in the Newton iteration
      SX gamma("gamma");
      SXMatrix jac = SXMatrix::eye(nx_) - gamma * f.jac(DAE_X,DAE_ODE);
      
      // Jacobian function
      vector<SXMatrix> jac_in(Sundials::M_NUM_IN);
      jac_in[M_T] = f.inputSX(DAE_T);
      jac_in[M_Y] = f.inputSX(DAE_X);
      jac_in[M_P] = f.inputSX(DAE_P);
      jac_in[M_GAMMA] = gamma;
      SXFunction M(jac_in,jac);
      
      // Pass sparsity to linear solver
      linsol_.setSparsity(jac.sparsity());
      
      // Save function
      jac_ = M;
    }
  }
  
  if(!jac_.isNull()) jac_.init();
  if(!linsol_.isNull()) linsol_.init();

  // Get the number of forward and adjoint directions
  nfdir_f_ = f_.getOption("number_of_fwd_dir");
  nadir_f_ = f_.getOption("number_of_adj_dir");

  // Set state derivative and its derivatives to zero (explicit integrator)
  f_.input(DAE_XDOT).setAll(0);
  for(int i=0; i<nfdir_f_; ++i) f_.fwdSeed(DAE_XDOT,i).setAll(0);
  for(int i=0; i<nadir_f_; ++i) f_.adjSens(DAE_XDOT,i).setAll(0);
  
  // Sundials return flag
  int flag;

  if(getOption("linear_multistep_method")=="adams")  lmm_ = CV_ADAMS;
  else if(getOption("linear_multistep_method")=="bdf") lmm_ = CV_BDF;
  else throw CasadiException("Unknown linear multistep method");

  if(getOption("nonlinear_solver_iteration")=="newton") iter_ = CV_NEWTON;
  else if(getOption("nonlinear_solver_iteration")=="functional") iter_ = CV_FUNCTIONAL;
  else throw CasadiException("Unknown nonlinear solver iteration");

  // Create CVodes memory block
  mem_ = CVodeCreate(lmm_,iter_);
  if(mem_==0) throw CasadiException("CVodeCreate: Creation failed");

  // Allocate n-vectors for ivp
  x0_ = N_VMake_Serial(nx_,&input(INTEGRATOR_X0).front());
  x_ = N_VMake_Serial(nx_,&output(INTEGRATOR_XF).front());

  // Disable internal warning messages?
  disable_internal_warnings_ = getOption("disable_internal_warnings");

  // Set error handler function
  flag = CVodeSetErrHandlerFn(mem_, ehfun_wrapper, this);
  if(flag != CV_SUCCESS) cvodes_error("CVodeSetErrHandlerFn",flag);

  // Initialize CVodes
  double t0 = 0;
  flag = CVodeInit(mem_, rhs_wrapper, t0, x0_);
  if(flag!=CV_SUCCESS) cvodes_error("CVodeInit",flag);

  // Set tolerances
  flag = CVodeSStolerances(mem_, reltol_, abstol_);
  if(flag!=CV_SUCCESS) cvodes_error("CVodeInit",flag);

  // Maximum number of steps
  CVodeSetMaxNumSteps(mem_, getOption("max_num_steps").toInt());
  if(flag != CV_SUCCESS) cvodes_error("CVodeSetMaxNumSteps",flag);
  
  // attach a linear solver
  if(getOption("linear_solver")=="dense"){
    // Dense jacobian
    flag = CVDense(mem_, nx_);
    if(flag!=CV_SUCCESS) cvodes_error("CVDense",flag);
    if(exact_jacobian_){
      // Create jacobian if it does not exist
      if(jac_f_.isNull()) jac_f_ = f_.jacobian(DAE_X,DAE_ODE);
      jac_f_.init();
      
      // Pass to CVodes
      flag = CVDlsSetDenseJacFn(mem_, djac_wrapper);
      if(flag!=CV_SUCCESS) cvodes_error("CVDlsSetDenseJacFn",flag);
      
    }
  } else if(getOption("linear_solver")=="banded") {
    // Banded jacobian
    flag = CVBand(mem_, nx_, getOption("upper_bandwidth").toInt(), getOption("lower_bandwidth").toInt());
    if(flag!=CV_SUCCESS) cvodes_error("CVBand",flag);
    if(exact_jacobian_){
      flag = CVDlsSetBandJacFn(mem_, bjac_wrapper);
      if(flag!=CV_SUCCESS) cvodes_error("CVDlsSetBandJacFn",flag);
    }
  } else if(getOption("linear_solver")=="iterative") {
    // Sparse (iterative) solver  

    // Preconditioning type
    int pretype;
    if(getOption("pretype")=="none")               pretype = PREC_NONE;
    else if(getOption("pretype")=="left")          pretype = PREC_LEFT;
    else if(getOption("pretype")=="right")         pretype = PREC_RIGHT;
    else if(getOption("pretype")=="both")          pretype = PREC_BOTH;
    else                                           throw CasadiException("Unknown preconditioning type");

    // Max dimension of the Krylov space
    int maxl = getOption("max_krylov");

    // Attach the sparse solver  
    if(getOption("iterative_solver")=="gmres"){
      flag = CVSpgmr(mem_, pretype, maxl);
      if(flag!=CV_SUCCESS) cvodes_error("CVBand",flag);      
    } else if(getOption("iterative_solver")=="bcgstab") {
      flag = CVSpbcg(mem_, pretype, maxl);
      if(flag!=CV_SUCCESS) cvodes_error("CVSpbcg",flag);
    } else if(getOption("iterative_solver")=="tfqmr") {
      flag = CVSptfqmr(mem_, pretype, maxl);
      if(flag!=CV_SUCCESS) cvodes_error("CVSptfqmr",flag);
    } else throw CasadiException("CVODES: Unknown sparse solver");
      
    // Attach functions for jacobian information
    if(exact_jacobian_){
      flag = CVSpilsSetJacTimesVecFn(mem_, jtimes_wrapper);
      if(flag!=CV_SUCCESS) cvodes_error("CVSpilsSetJacTimesVecFn",flag);      
    }
    
    // Add a preconditioner
    if(bool(getOption("use_preconditioner"))){
      // Make sure that a Jacobian has been provided
      if(jac_.isNull()) throw CasadiException("CVodesInternal::init(): No Jacobian has been provided.");

      // Make sure that a linear solver has been providided
      if(linsol_.isNull()) throw CasadiException("CVodesInternal::init(): No user defined linear solver has been provided.");

      // Pass to IDA
      flag = CVSpilsSetPreconditioner(mem_, psetup_wrapper, psolve_wrapper);
      if(flag != CV_SUCCESS) cvodes_error("CVSpilsSetPreconditioner",flag);
    }    
  } else if(getOption("linear_solver")=="user_defined") {
    initUserDefinedLinearSolver();
  } else throw CasadiException("Unknown linear solver ");

  // Set user data
  flag = CVodeSetUserData(mem_,this);
  if(flag!=CV_SUCCESS) cvodes_error("CVodeSetUserData",flag);

  // Quadrature equations
  if(nq_>0){
    // Allocate n-vectors for quadratures
    q_ = N_VMake_Serial(nq_,&output(INTEGRATOR_QF).front());

    // Initialize quadratures in CVodes
    N_VConst(0.0, q_);
    flag = CVodeQuadInit(mem_, rhsQ_wrapper, q_);
    if(flag != CV_SUCCESS) cvodes_error("CVodeQuadInit",flag);
    
    // Should the quadrature errors be used for step size control?
    if(getOption("quad_err_con").toInt()){
      flag = CVodeSetQuadErrCon(mem_, true);
      if(flag != CV_SUCCESS) cvodes_error("IDASetQuadErrCon",flag);
      
      // Quadrature error tolerances
      flag = CVodeQuadSStolerances(mem_, reltol_, abstol_); // TODO: vector absolute tolerances
      if(flag != CV_SUCCESS) cvodes_error("CVodeQuadSStolerances",flag);
    }
  }
  
    // Forward sensitivity problem
    if(nfdir_>0){
      // Allocate n-vectors
      xF0_.resize(nfdir_,0);
      xF_.resize(nfdir_,0);
      for(int i=0; i<nfdir_; ++i){
        xF0_[i] = N_VMake_Serial(nx_,&fwdSeed(INTEGRATOR_X0,i).data()[0]);
        xF_[i] = N_VMake_Serial(nx_,&fwdSens(INTEGRATOR_XF,i).data()[0]);
      }

      // Allocate n-vectors for quadratures
      if(nq_>0){
        qF_.resize(nfdir_,0);
        for(int i=0; i<nfdir_; ++i){
          qF_[i] = N_VMake_Serial(nq_,&fwdSens(INTEGRATOR_QF,i).front());
        }
      }
      
    // Calculate all forward sensitivity right hand sides at once?
    bool all_at_once = getOption("fsens_all_at_once");
      
    // Get the sensitivity method
    if(getOption("sensitivity_method")=="simultaneous") ism_ = CV_SIMULTANEOUS;
    else if(getOption("sensitivity_method")=="staggered") ism_ = all_at_once ? CV_STAGGERED : CV_STAGGERED1;
    else throw CasadiException("CVodes: Unknown sensitivity method");
  
    // Initialize forward sensitivities
    if(finite_difference_fsens_){
      // Use finite differences to calculate the residual in the forward sensitivity equations
      if(all_at_once){
        flag = CVodeSensInit(mem_,nfdir_,ism_,0,getPtr(xF0_));
        if(flag != CV_SUCCESS) cvodes_error("CVodeSensInit",flag);
      } else {
        flag = CVodeSensInit1(mem_,nfdir_,ism_,0,getPtr(xF0_));
        if(flag != CV_SUCCESS) cvodes_error("CVodeSensInit1",flag);
      }
      
      // Pass pointer to parameters
      flag = CVodeSetSensParams(mem_,&input(INTEGRATOR_P).data().front(),0,0);
      if(flag != CV_SUCCESS) cvodes_error("CVodeSetSensParams",flag);

      //  CVodeSetSensDQMethod

    } else {
      if(all_at_once){
        // Use AD to calculate the residual in the forward sensitivity equations
        flag = CVodeSensInit(mem_,nfdir_,ism_,rhsS_wrapper,getPtr(xF0_));
        if(flag != CV_SUCCESS) cvodes_error("CVodeSensInit",flag);
      } else {
        flag = CVodeSensInit1(mem_,nfdir_,ism_,rhsS1_wrapper,getPtr(xF0_));
        if(flag != CV_SUCCESS) cvodes_error("CVodeSensInit",flag);
      }
    }
    
    // Set tolerances
    vector<double> fsens_abstol(nfdir_,fsens_abstol_);
    
    flag = CVodeSensSStolerances(mem_,fsens_reltol_,getPtr(fsens_abstol));
    if(flag != CV_SUCCESS) cvodes_error("CVodeSensSStolerances",flag);
    
    // Set optional inputs
    bool errconS = getOption("fsens_err_con");
    flag = CVodeSetSensErrCon(mem_, errconS);
    if(flag != CV_SUCCESS) cvodes_error("CVodeSetSensErrCon",flag);
    
    // Quadrature equations
    if(nq_>0){
      for(vector<N_Vector>::iterator it=qF_.begin(); it!=qF_.end(); ++it) N_VConst(0.0,*it);
      flag = CVodeQuadSensInit(mem_, rhsQS_wrapper, getPtr(qF_));
      if(flag != CV_SUCCESS) cvodes_error("CVodeQuadSensInit",flag);

      // Set tolerances
      flag = CVodeQuadSensSStolerances(mem_,fsens_reltol_,getPtr(fsens_abstol));
      if(flag != CV_SUCCESS) cvodes_error("CVodeQuadSensSStolerances",flag);
    }
  } // enable fsens
    
  // Adjoint sensitivity problem
  whichB_.resize(nadir_);

  // Allocate n-vectors
  xA0_.resize(nadir_,0);
  xA_.resize(nadir_,0);
  for(int i=0; i<nadir_; ++i){
    xA0_[i] = N_VMake_Serial(nx_,&adjSeed(INTEGRATOR_XF,i).front());
    xA_[i] = N_VMake_Serial(nx_,&adjSens(INTEGRATOR_X0,i).front());
  }

  // Allocate n-vectors for quadratures
  qA_.resize(nadir_,0);
  for(int i=0; i<nadir_; ++i){
    casadi_assert(adjSens(INTEGRATOR_P,i).size()==np_);
    qA_[i] = N_VMake_Serial(np_,&adjSens(INTEGRATOR_P,i).front());
  }
  
  if(nadir_>0){
    // Get the number of steos per checkpoint
    int Nd = getOption("steps_per_checkpoint");

    // Get the interpolation type
    int interpType;
    if(getOption("interpolation_type")=="hermite")
      interpType = CV_HERMITE;
    else if(getOption("interpolation_type")=="polynomial")
      interpType = CV_POLYNOMIAL;
    else throw CasadiException("\"interpolation_type\" must be \"hermite\" or \"polynomial\"");
      
    // Initialize adjoint sensitivities
    flag = CVodeAdjInit(mem_, Nd, interpType);
    if(flag != CV_SUCCESS) cvodes_error("CVodeAdjInit",flag);
  }
        
  is_init = true;
  isInitAdj_ = false;
}


void CVodesInternal::initAdj(){
  int flag;

  for(int dir=0; dir<nadir_; ++dir){

      // Create backward problem (use the same lmm and iter)
      flag = CVodeCreateB(mem_, lmm_, iter_, &whichB_[dir]);
      if(flag != CV_SUCCESS) cvodes_error("CVodeCreateB",flag);
      
      // Initialize the backward problem
      double tB0 = tf_;
      flag = CVodeInitB(mem_, whichB_[dir], rhsB_wrapper, tB0, xA0_[dir]);
      if(flag != CV_SUCCESS) cvodes_error("CVodeInitB",flag);

      // Set tolerances
      flag = CVodeSStolerancesB(mem_, whichB_[dir], asens_reltol_, asens_abstol_);
      if(flag!=CV_SUCCESS) cvodes_error("CVodeSStolerancesB",flag);

      // User data
      flag = CVodeSetUserDataB(mem_, whichB_[dir], this);
      if(flag != CV_SUCCESS) cvodes_error("CVodeSetUserDataB",flag);

      // attach linear solver to backward problem
      if(getOption("asens_linear_solver")=="dense"){
      // Dense jacobian
      flag = CVDenseB(mem_, whichB_[dir], nx_);
      if(flag!=CV_SUCCESS) cvodes_error("CVDenseB",flag);
      } else if(getOption("asens_linear_solver")=="banded") {
      // Banded jacobian
      flag = CVBandB(mem_, whichB_[dir], nx_, getOption("asens_upper_bandwidth").toInt(), getOption("asens_lower_bandwidth").toInt());
      if(flag!=CV_SUCCESS) cvodes_error("CVBandB",flag);
    } else if(getOption("asens_linear_solver")=="iterative") {
      // Sparse jacobian
      // Preconditioning type
      int pretype;
      if(getOption("asens_pretype")=="none")               pretype = PREC_NONE;
      else if(getOption("asens_pretype")=="left")          pretype = PREC_LEFT;
      else if(getOption("asens_pretype")=="right")         pretype = PREC_RIGHT;
      else if(getOption("asens_pretype")=="both")          pretype = PREC_BOTH;
      else                                            throw CasadiException("Unknown preconditioning type for backward problem");

      // Attach the sparse solver  
      int maxl = getOption("asens_max_krylov");
      if(getOption("asens_iterative_solver")=="gmres"){
        flag = CVSpgmrB(mem_, whichB_[dir], pretype, maxl);
        if(flag!=CV_SUCCESS) cvodes_error("CVSpgmrB",flag);      
      } else if(getOption("asens_iterative_solver")=="bcgstab") {
        flag = CVSpbcgB(mem_, whichB_[dir], pretype, maxl);
        if(flag!=CV_SUCCESS) cvodes_error("CVSpbcgB",flag);
      } else if(getOption("asens_iterative_solver")=="tfqmr") {
        flag = CVSptfqmrB(mem_, whichB_[dir], pretype, maxl);
        if(flag!=CV_SUCCESS) cvodes_error("CVSptfqmrB",flag);
      } else throw CasadiException("Unknown sparse solver for backward problem");
        
      // Attach functions for jacobian information
    } else throw CasadiException("Unknown linear solver for backward problem");

    // Quadratures for the adjoint problem

    N_VConst(0.0, qA_[dir]);
    flag = CVodeQuadInitB(mem_,whichB_[dir],rhsQB_wrapper,qA_[dir]);
    if(flag!=CV_SUCCESS) cvodes_error("CVodeQuadInitB",flag);
    
    if(getOption("quad_err_con").toInt()){
      flag = CVodeSetQuadErrConB(mem_, whichB_[dir],true);
      if(flag != CV_SUCCESS) cvodes_error("CVodeSetQuadErrConB",flag);
        
      flag = CVodeQuadSStolerancesB(mem_, whichB_[dir], asens_reltol_, asens_abstol_);
      if(flag != CV_SUCCESS) cvodes_error("CVodeQuadSStolerancesB",flag);
    }
  }
  
  // Mark initialized
  isInitAdj_ = true;
}

void CVodesInternal::rhs(double t, const double* x, double* xdot){
  if(monitor_rhs_){
      cout << "CVodesInternal::rhs: begin" << endl;
  }

  // Get time
  time1 = clock();

  // Pass input
  f_.setInput(&t,DAE_T);
  f_.setInput(x,DAE_X);
  f_.setInput(input(INTEGRATOR_P),DAE_P);

  if(monitor_rhs_) {
    cout << "t       = " << t << endl;
    cout << "x       = " << f_.input(DAE_X) << endl;
    cout << "p       = " << f_.input(DAE_P) << endl;
  }
    // Evaluate
  f_.evaluate();

  if(monitor_rhs_) {
    cout << "xdot       = " << f_.output(DAE_ODE)<< endl;
  }
    
  // Get results
  f_.getOutput(xdot);

  // Log time
  time2 = clock();
  t_res += double(time2-time1)/CLOCKS_PER_SEC;

  if(monitor_rhs_){
   cout << "CVodesInternal::rhs: end" << endl;
  }

}

int CVodesInternal::rhs_wrapper(double t, N_Vector x, N_Vector xdot, void *user_data){
try{
    casadi_assert(user_data);
    CVodesInternal *this_ = (CVodesInternal*)user_data;
    this_->rhs(t,NV_DATA_S(x),NV_DATA_S(xdot));
    return 0;
  } catch(exception& e){
    cerr << "rhs failed: " << e.what() << endl;
    return 1;
  }
}
  
void CVodesInternal::reset(int nfdir, int nadir){
  if(monitored("reset")){
    cout << "initial state: " << endl;
    cout << "p = " << input(INTEGRATOR_P) << endl;
    cout << "x0 = " << input(INTEGRATOR_X0) << endl;
  }

  // Reset timers
  t_res = t_fres = t_jac = t_lsolve = t_lsetup_jac = t_lsetup_fac = 0;
  
  fsens_order_ = nfdir>0;
  asens_order_ = nadir>0;
  
  // Get the time horizon
  t_ = t0_;

  // Re-initialize
  int flag = CVodeReInit(mem_, t0_, x0_);
  if(flag!=CV_SUCCESS) cvodes_error("CVodeReInit",flag);
  
  // Re-initialize quadratures
  if(nq_>0){
    N_VConst(0.0,q_);
    flag = CVodeQuadReInit(mem_, q_);
    if(flag != CV_SUCCESS) cvodes_error("CVodeQuadReInit",flag);
  }
  
  // Re-initialize sensitivities
  if(fsens_order_>0){
    flag = CVodeSensReInit(mem_,ism_,getPtr(xF0_));
    if(flag != CV_SUCCESS) cvodes_error("CVodeSensReInit",flag);
    
    if(nq_>0){
      for(vector<N_Vector>::iterator it=qF_.begin(); it!=qF_.end(); ++it) N_VConst(0.0,*it);
      flag = CVodeQuadSensReInit(mem_, getPtr(qF_));
      if(flag != CV_SUCCESS) cvodes_error("CVodeQuadSensReInit",flag);
    }
  } else {
    // Turn of sensitivities
    flag = CVodeSensToggleOff(mem_);
    if(flag != CV_SUCCESS) cvodes_error("CVodeSensToggleOff",flag);
  }
  
  // Set the stop time of the integration -- don't integrate past this point
  if(stop_at_end_) setStopTime(tf_);
}

void CVodesInternal::integrate(double t_out){
   log("CVODES::integrate begin");
  int flag;
    
  // tolerance
  double ttol = 1e-9;
  if(fabs(t_-t_out)<ttol){
    copy(input(INTEGRATOR_X0).begin(),input(INTEGRATOR_X0).end(),output(INTEGRATOR_XF).begin());
    if(fsens_order_>0){
      for(int i=0; i<nfdir_; ++i){
        copy(fwdSeed(INTEGRATOR_X0,i).begin(),fwdSeed(INTEGRATOR_X0,i).end(),fwdSens(INTEGRATOR_XF,i).begin());
      }
    }
    return;
  }
  if(asens_order_>0){
    int ncheck; // number of checkpoints stored so far
    flag = CVodeF(mem_, t_out, x_, &t_, CV_NORMAL,&ncheck);
    if(flag!=CV_SUCCESS && flag!=CV_TSTOP_RETURN) cvodes_error("CVodeF",flag);
    
  } else {
    flag = CVode(mem_, t_out, x_, &t_, CV_NORMAL);
    if(flag!=CV_SUCCESS && flag!=CV_TSTOP_RETURN) cvodes_error("CVode",flag);
  }
  
  if(nq_>0){
    double tret;
    flag = CVodeGetQuad(mem_, &tret, q_);
    if(flag!=CV_SUCCESS) cvodes_error("CVodeGetQuad",flag);
  }
  
  if(fsens_order_>0){
    // Get the sensitivities
    flag = CVodeGetSens(mem_, &t_, getPtr(xF_));
    if(flag != CV_SUCCESS) cvodes_error("CVodeGetSens",flag);
    
    if(nq_>0){
      double tret;
      flag = CVodeGetQuadSens(mem_, &tret, getPtr(qF_));
      if(flag != CV_SUCCESS) cvodes_error("CVodeGetQuadSens",flag);
    }
  }

  
  // Print statistics
  if(getOption("print_stats")) printStats(std::cout);
  
  log("CVODES::integrate end");
}

void CVodesInternal::resetAdj(){
  int flag;
  
  if(isInitAdj_){
    for(int dir=0; dir<nadir_; ++dir){
      flag = CVodeReInitB(mem_, whichB_[dir], tf_, xA0_[dir]);
      if(flag != CV_SUCCESS) cvodes_error("CVodeReInitB",flag);

      N_VConst(0.0,qA_.at(dir));
      flag = CVodeQuadReInitB(mem_,whichB_[dir],qA_[dir]);
      if(flag!=CV_SUCCESS) cvodes_error("CVodeQuadReInitB",flag);
    }
  } else {
    // Initialize the adjoint integration
    initAdj();
  }
}

void CVodesInternal::integrateAdj(double t_out){
  int flag;
  
  // Integrate backwards to t_out
  flag = CVodeB(mem_, t_out, CV_NORMAL);
  if(flag<CV_SUCCESS) cvodes_error("CVodeB",flag);

  // Get the sensitivities
  double tret;
  for(int dir=0; dir<nadir_; ++dir){
    flag = CVodeGetB(mem_, whichB_[dir], &tret, xA_[dir]);
    if(flag!=CV_SUCCESS) cvodes_error("CVodeGetB",flag);

    flag = CVodeGetQuadB(mem_, whichB_[dir], &tret, qA_[dir]);
    if(flag!=CV_SUCCESS) cvodes_error("CVodeGetQuadB",flag);
  }
}

void CVodesInternal::printStats(std::ostream &stream) const{

    long nsteps, nfevals, nlinsetups, netfails;
    int qlast, qcur;
    double hinused, hlast, hcur, tcur;
    
    int flag;

    flag = CVodeGetIntegratorStats(mem_, &nsteps, &nfevals,&nlinsetups, &netfails, &qlast, &qcur,&hinused, &hlast, &hcur, &tcur);
    if(flag!=CV_SUCCESS) cvodes_error("CVodeGetIntegratorStats",flag);


    long nfevalsLS_spils=0, nfevalsLS_dls=0, nfevalsBP=0;
    
    flag = flag = CVSpilsGetNumRhsEvals(mem_, &nfevalsLS_spils);
    if(flag!=CV_SUCCESS) cvodes_error("CVSpilsGetNumRhsEvals",flag);
    
    flag = CVDlsGetNumRhsEvals(mem_, &nfevalsLS_dls);
    if(flag!=CV_SUCCESS) cvodes_error("CVDlsGetNumRhsEvals",flag);
    
    //flag = CVBandPrecGetNumRhsEvals(mem_, &nfevalsBP);
    //if(flag!=CV_SUCCESS) cvodes_error("CVBandPrecGetNumRhsEvals",flag);
    
    stream << "number of steps taken by CVODES: " << nsteps << std::endl;
    stream << "number of calls to the user's f function: " << nfevals + nfevalsLS_spils + nfevalsLS_dls << std::endl;
    stream << "   main solver                      : " << nfevals << std::endl;
    stream << "   finite diff (SPILS linear solver): " << nfevalsLS_spils << std::endl;
    stream << "   finite diff (DLS linear solver)  : " << nfevalsLS_dls << std::endl;
    //stream << "   finite diff (preconditionar)     : " << nfevalsBP << std::endl;
    stream << "number of calls made to the linear solver setup function: " << nlinsetups << std::endl;
    stream << "number of error test failures: " << netfails << std::endl;
    stream << "method order used on the last internal step: " << qlast << std::endl;
    stream << "method order to be used on the next internal step: " << qcur << std::endl;
    stream << "actual value of initial step size: " << hinused << std::endl;
    stream << "step size taken on the last internal step: " << hlast << std::endl;
    stream << "step size to be attempted on the next internal step: " << hcur << std::endl;
    stream << "current internal time reached: " << tcur << std::endl;
    stream << std::endl;

    stream << "Time spent in the ODE residual: " << t_res << " s." << endl;
    stream << "Time spent in the forward sensitivity residual: " << t_fres << " s." << endl;
    stream << "Time spent in the jacobian function or jacobian times vector function: " << t_jac << " s." << endl;
    stream << "Time spent in the linear solver solve function: " << t_lsolve << " s." << endl;
    stream << "Time spent to generate the jacobian in the linear solver setup function: " << t_lsetup_jac << " s." << endl;
    stream << "Time spent to factorize the jacobian in the linear solver setup function: " << t_lsetup_fac << " s." << endl;
    stream << std::endl;

    
    
#if 0
  // Quadrature
  if(ops.quadrature && ocp.hasFunction(LTERM)){
      long nfQevals, nQetfails;
      flag = CVodeGetQuadStats(cvode_mem_[k], &nfQevals, &nQetfails);  
      if(flag != CV_SUCCESS) throw "Error in CVodeGetQuadStats";

      stream << "Quadrature: " << std::endl;
      stream << "number of calls made to the user's quadrature right-hand side function: " << nfQevals << std::endl;
      stream << "number of local error test failures due to quadrature variables: " <<  nQetfails << std::endl;
      stream << std::endl;
}
#endif
}
  
map<int,string> CVodesInternal::calc_flagmap(){
  map<int,string> f;
  f[CV_SUCCESS] = "CV_SUCCESS";
  f[CV_TSTOP_RETURN] = "CV_TSTOP_RETURN";
  f[CV_ROOT_RETURN] = "CV_ROOT_RETURN";
  f[CV_WARNING] = "CV_WARNING";
  f[CV_WARNING] = "CV_WARNING";
  f[CV_TOO_MUCH_WORK] = "CV_TOO_MUCH_WORK";
  f[CV_TOO_MUCH_ACC] = "CV_TOO_MUCH_ACC";
  f[CV_ERR_FAILURE] = "CV_ERR_FAILURE";
  f[CV_CONV_FAILURE] = "CV_CONV_FAILURE";
  f[CV_LINIT_FAIL] = "CV_LINIT_FAIL";
  f[CV_LSETUP_FAIL] = "CV_LSETUP_FAIL";
  f[CV_LSOLVE_FAIL] = "CV_LSOLVE_FAIL";
  f[CV_RHSFUNC_FAIL] = "CV_RHSFUNC_FAIL";
  f[CV_FIRST_RHSFUNC_ERR] = "CV_FIRST_RHSFUNC_ERR";
  f[CV_UNREC_RHSFUNC_ERR] = "CV_UNREC_RHSFUNC_ERR";
  f[CV_RTFUNC_FAIL] = "CV_RTFUNC_FAIL";
  f[CV_MEM_FAIL] = "CV_MEM_FAIL";
  f[CV_ILL_INPUT] = "CV_ILL_INPUT";
  f[CV_NO_MALLOC] = "CV_NO_MALLOC";
  f[CV_BAD_K] = "CV_BAD_K";
  f[CV_BAD_T] = "CV_BAD_T";
  f[CV_BAD_DKY] = "CV_BAD_DKY";
  f[CV_TOO_CLOSE] = "CV_TOO_CLOSE";
  f[CV_QRHSFUNC_FAIL] = "CV_QRHSFUNC_FAIL";
  f[CV_FIRST_QRHSFUNC_ERR] = "CV_FIRST_QRHSFUNC_ERR";
  f[CV_REPTD_QRHSFUNC_ERR] = "CV_REPTD_QRHSFUNC_ERR";
  f[CV_UNREC_QRHSFUNC_ERR] = "CV_UNREC_QRHSFUNC_ERR";
  f[CV_NO_SENS ] = "CV_NO_SENS ";
  f[CV_SRHSFUNC_FAIL] = "CV_SRHSFUNC_FAIL";
  return f;
}
  
map<int,string> CVodesInternal::flagmap = CVodesInternal::calc_flagmap();

void CVodesInternal::cvodes_error(const string& module, int flag){
  // Find the error
  map<int,string>::const_iterator it = flagmap.find(flag);
  
  stringstream ss;
  if(it == flagmap.end()){
    ss << "Unknown error (" << flag << ") from module \"" << module << "\".";
  } else {
    ss << "Module \"" << module << "\" returned flag \"" << it->second << "\".";
  }
  ss << " Consult Cvodes documentation.";
  casadi_error(ss.str());
}
  
void CVodesInternal::ehfun_wrapper(int error_code, const char *module, const char *function, char *msg, void *eh_data){
try{
    casadi_assert(eh_data);
    CVodesInternal *this_ = (CVodesInternal*)eh_data;
    this_->ehfun(error_code,module,function,msg);        
  } catch(exception& e){
    cerr << "ehfun failed: " << e.what() << endl;
  }
}
  
void CVodesInternal::ehfun(int error_code, const char *module, const char *function, char *msg){
  if(!disable_internal_warnings_){
    cerr << msg << endl;
  }
}

void CVodesInternal::rhsS(int Ns, double t, N_Vector x, N_Vector xdot, N_Vector *xF, N_Vector *xdotF, N_Vector tmp1, N_Vector tmp2){
  casadi_assert(Ns==nfdir_);

  // Record the current cpu time
  time1 = clock();
  
    // Pass input
  f_.setInput(&t,DAE_T);
  f_.setInput(NV_DATA_S(x),DAE_X);
  f_.setInput(input(INTEGRATOR_P),DAE_P);

   // Calculate the forward sensitivities, nfdir_f_ directions at a time
   for(int j=0; j<nfdir_; j += nfdir_f_){
     for(int dir=0; dir<nfdir_f_ && j+dir<nfdir_; ++dir){
       // Pass forward seeds 
       f_.fwdSeed(DAE_T,dir).setZero();
       f_.setFwdSeed(NV_DATA_S(xF[j+dir]),DAE_X,dir);
       f_.setFwdSeed(fwdSeed(INTEGRATOR_P,j+dir),DAE_P,dir);
     }

     // Evaluate the AD forward algorithm
     f_.evaluate(nfdir_f_,0);
      
     // Get the output seeds
     for(int dir=0; dir<nfdir_f_ && j+dir<nfdir_; ++dir){
       f_.getFwdSens(NV_DATA_S(xdotF[j+dir]),DAE_ODE,dir);
     }
   }
  
  // Record timings
  time2 = clock();
  t_fres += double(time2-time1)/CLOCKS_PER_SEC;
}

int CVodesInternal::rhsS_wrapper(int Ns, double t, N_Vector x, N_Vector xdot, N_Vector *xF, N_Vector *xdotF, void *user_data, N_Vector tmp1, N_Vector tmp2){
try{
    casadi_assert(user_data);
    CVodesInternal *this_ = (CVodesInternal*)user_data;
    this_->rhsS(Ns,t,x,xdot,xF,xdotF,tmp1,tmp2);
    return 0;
  } catch(exception& e){
    cerr << "fs failed: " << e.what() << endl;
    return 1;
  }
}

void CVodesInternal::rhsS1(int Ns, double t, N_Vector x, N_Vector xdot, int iS, N_Vector xF, N_Vector xdotF, N_Vector tmp1, N_Vector tmp2){
  casadi_assert(Ns==nfdir_);
  
    // Pass input
  f_.setInput(&t,DAE_T);
  f_.setInput(NV_DATA_S(x),DAE_X);
  f_.setInput(input(INTEGRATOR_P),DAE_P);

  // Pass forward seeds
  f_.fwdSeed(DAE_T).setZero();
  f_.setFwdSeed(NV_DATA_S(xF),DAE_X);
  f_.setFwdSeed(fwdSeed(INTEGRATOR_P,iS),DAE_P);
    
  // Evaluate the AD forward algorithm
  f_.evaluate(1,0);
  
  // Get the fwd sensitivities
  f_.getFwdSens(NV_DATA_S(xdotF));
}

int CVodesInternal::rhsS1_wrapper(int Ns, double t, N_Vector x, N_Vector xdot, int iS, N_Vector xF, N_Vector xdotF, void *user_data, N_Vector tmp1, N_Vector tmp2){
try{
    casadi_assert(user_data);
    CVodesInternal *this_ = (CVodesInternal*)user_data;
    this_->rhsS1(Ns,t,x,xdot,iS,xF,xdotF,tmp1,tmp2);
    return 0;
  } catch(exception& e){
    cerr << "fs failed: " << e.what() << endl;
    return 1;
  }
}

int CVodesInternal::rhsQ_wrapper(double t, N_Vector x, N_Vector qdot, void *user_data){
try{
    casadi_assert(user_data);
    CVodesInternal *this_ = (CVodesInternal*)user_data;
    this_->rhsQ(t,NV_DATA_S(x),NV_DATA_S(qdot));
    return 0;
  } catch(exception& e){
    cerr << "rhsQ failed: " << e.what() << endl;;
    return 1;
  }
}

void CVodesInternal::rhsQ(double t, const double* x, double* qdot){
  // Pass input
  f_.setInput(&t,DAE_T);
  f_.setInput(x,DAE_X);
  f_.setInput(input(INTEGRATOR_P),DAE_P);

  // Evaluate
  f_.evaluate();
    
  // Get results
  f_.getOutput(qdot,DAE_QUAD);
}

void CVodesInternal::rhsQS(int Ns, double t, N_Vector x, N_Vector *xF, N_Vector qdot, N_Vector *qdotF, N_Vector tmp1, N_Vector tmp2){
  casadi_assert(Ns==nfdir_);
  
  // Pass input
  f_.setInput(&t,DAE_T);
  f_.setInput(NV_DATA_S(x),DAE_X);
  f_.setInput(input(INTEGRATOR_P),DAE_P);

  for(int i=0; i<nfdir_; ++i){
    // Pass forward seeds
    f_.fwdSeed(DAE_T).setZero();
    f_.setFwdSeed(NV_DATA_S(xF[i]),DAE_X);
    f_.setFwdSeed(fwdSeed(INTEGRATOR_P,i),DAE_P);

    // Evaluate the AD forward algorithm
    f_.evaluate(1,0);
      
    // Get the forward sensitivities
    f_.getFwdSens(NV_DATA_S(qdotF[i]),DAE_QUAD);
  }
}

int CVodesInternal::rhsQS_wrapper(int Ns, double t, N_Vector x, N_Vector *xF, N_Vector qdot, N_Vector *qdotF, void *user_data, N_Vector tmp1, N_Vector tmp2){
try{
//    casadi_assert(user_data);
    CVodesInternal *this_ = (CVodesInternal*)user_data;
    if(!this_){
      // SUNDIALS BUG!!!
      for(int i=0; i<Ns; ++i) N_VConst(0.0,qdotF[i]);
      return 0;
    }
    this_->rhsQS(Ns,t,x,xF,qdot,qdotF,tmp1,tmp2);
    return 0;
  } catch(exception& e){
    cerr << "rhsQS failed: " << e.what() << endl;;
    return 1;
  }
}

void CVodesInternal::rhsB(double t, const double* x, const double *xA, double* xzdotA){
  if(monitor_rhsB_){
    cout << "CVodesInternal::rhsB: begin" << endl;
  }
    
  // Pass input
  f_.setInput(&t,DAE_T);
  f_.setInput(x,DAE_X);
  f_.setInput(input(INTEGRATOR_P),DAE_P);

  // Pass adjoint seeds
  f_.setAdjSeed(xA,DAE_ODE);
  if(nq_>0){
    f_.setAdjSeed(adjSeed(INTEGRATOR_QF),DAE_QUAD);
  }

  if(monitor_rhsB_){
    cout << "t       = " << t << endl;
    cout << "y       = " << f_.input(DAE_X) << endl;
    cout << "p       = " << f_.input(DAE_P) << endl;
    cout << "aseed   = " << f_.adjSeed(DAE_ODE) << endl;
  }
  
  // Evaluate and tape
  f_.evaluate(0,1);

  if(monitor_rhsB_){
    cout << "f_asens = " << f_.adjSens(DAE_X) << endl;
  }
  
  // Save to output
  const vector<double>& fres = f_.adjSens(DAE_X).data();
  for(int i=0; i<nx_; ++i)
    xzdotA[i] = -fres[i];

  if(monitor_rhsB_){
    cout << "CVodesInternal::rhsB: end" << endl;
  }

}

int CVodesInternal::rhsB_wrapper(double t, N_Vector x, N_Vector xA, N_Vector xzdotA, void *user_data){
try{
    casadi_assert(user_data);
    CVodesInternal *this_ = (CVodesInternal*)user_data;
    this_->rhsB(t,NV_DATA_S(x),NV_DATA_S(xA),NV_DATA_S(xzdotA));
    return 0;
  } catch(exception& e){
    cerr << "rhsB failed: " << e.what() << endl;;
    return 1;
  }
}

int CVodesInternal::rhsQB_wrapper(double t, N_Vector x, N_Vector xA, N_Vector qdotA, void *user_data){
  try{
    casadi_assert(user_data);
    CVodesInternal *this_ = (CVodesInternal*)user_data;
    this_->rhsQB(t,NV_DATA_S(x),NV_DATA_S(xA),NV_DATA_S(qdotA));
    return 0;
  } catch(exception& e){
    cerr << "rhsQB failed: " << e.what() << endl;;
    return 1;
  }
}

void CVodesInternal::rhsQB(double t, const double* y, const double* xA, double* qdotA){
  if(monitor_rhsQB_){
    cout << "CVodesInternal::rhsQB: begin" << endl;
  }

  // Pass input
  f_.setInput(&t,DAE_T);
  f_.setInput(y,DAE_X);
  f_.setInput(input(INTEGRATOR_P),DAE_P);
  
  if(monitor_rhs_) {
    cout << "t       = " << t << endl;
    cout << "y       = " << f_.input(DAE_X) << endl;
    cout << "p       = " << f_.input(DAE_P) << endl;
  }

  // Pass adjoint seeds
  f_.setAdjSeed(xA,DAE_ODE);
  if(nq_>0){
    f_.setAdjSeed(adjSeed(INTEGRATOR_QF),DAE_QUAD);
  }
  if(monitor_rhsQB_) {
    cout << "adjSeed       = " << f_.adjSeed(DAE_ODE) << endl;
  }

  // Evaluate
  f_.evaluate(0,1);

  // Save to output
  f_.getAdjSens(qdotA,DAE_P);
  
  if(monitor_rhsQB_) {
    cout << "adjSens       = " << f_.adjSens(DAE_P) << endl;
  }
    
  // Negate as we are integrating backwards
  for(int i=0; i<np_; ++i)
    qdotA[i] *= -1;
    
  if(monitor_rhsQB_){
   cout << "CVodesInternal::rhsQB: end" << endl;
  }
}

int CVodesInternal::jtimes_wrapper(N_Vector v, N_Vector Jv, double t, N_Vector x, N_Vector xdot, void *user_data, N_Vector tmp){
  try{
    casadi_assert(user_data);
    CVodesInternal *this_ = (CVodesInternal*)user_data;
    casadi_assert(this_->f_.fwdSens(DAE_ODE).size() == this_->nx_);
    casadi_assert(NV_LENGTH_S(v) == this_->nx_);
    casadi_assert(NV_LENGTH_S(Jv) == this_->nx_);
    this_->jtimes(NV_DATA_S(v),NV_DATA_S(Jv),t,NV_DATA_S(x),NV_DATA_S(xdot),NV_DATA_S(tmp));
    return 0;
  } catch(exception& e){
    cerr << "jtimes failed: " << e.what() << endl;;
    return 1;
  }
}

void CVodesInternal::jtimes(const double *v, double* Jv, double t, const double* y, const double* xdot, double* tmp){
  // Get time
  time1 = clock();

  // Pass input
  f_.setInput(&t,DAE_T);
  f_.setInput(y,DAE_X);
  f_.setInput(input(INTEGRATOR_P),DAE_P);

  // Pass input seeds
  f_.fwdSeed(DAE_T).setZero();
  f_.setFwdSeed(v,DAE_X);
  fill_n(f_.fwdSeed(DAE_P).begin(),np_,0.0);
  
  // Evaluate
  f_.evaluate(1,0);

  // Get the output seeds
  f_.getFwdSens(Jv,DAE_ODE);
  
  // Log time duration
  time2 = clock();
  t_jac += double(time2-time1)/CLOCKS_PER_SEC;
}

int CVodesInternal::djac_wrapper(SUNDIALS_INT N, double t, N_Vector x, N_Vector xdot, DlsMat Jac, void *user_data,N_Vector tmp1, N_Vector tmp2, N_Vector tmp3){
  try{
    casadi_assert(user_data);
    CVodesInternal *this_ = (CVodesInternal*)user_data;
    this_->djac(N, t, x, xdot, Jac, tmp1, tmp2, tmp3);
    return 0;
  } catch(exception& e){
    cerr << "djac failed: " << e.what() << endl;;
    return 1;
  }
}

void CVodesInternal::djac(SUNDIALS_INT N, double t, N_Vector x, N_Vector xdot, DlsMat Jac, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3){
  // Get time
  time1 = clock();

  // Pass inputs to the jacobian function
  jac_f_.setInput(&t,DAE_T);
  jac_f_.setInput(NV_DATA_S(x),DAE_X);
  jac_f_.setInput(f_.input(DAE_P),DAE_P);

  // Evaluate
  jac_f_.evaluate();
  
  // Get sparsity and non-zero elements
  const vector<int>& rowind = jac_f_.output().rowind();
  const vector<int>& col = jac_f_.output().col();
  const vector<double>& val = jac_f_.output().data();

  // Loop over rows
  for(int i=0; i<rowind.size()-1; ++i){
    // Loop over non-zero entries
    for(int el=rowind[i]; el<rowind[i+1]; ++el){
      // Get column
      int j = col[el];
      
      // Set the element
      DENSE_ELEM(Jac,i,j) = val[el];
    }
  }
  
  // Log time duration
  time2 = clock();
  t_jac += double(time2-time1)/CLOCKS_PER_SEC;
}

int CVodesInternal::bjac_wrapper(SUNDIALS_INT N, SUNDIALS_INT mupper, SUNDIALS_INT mlower, double t, N_Vector x, N_Vector xdot, DlsMat Jac, void *user_data,     
                        N_Vector tmp1, N_Vector tmp2, N_Vector tmp3){
  try{
    casadi_assert(user_data);
    CVodesInternal *this_ = (CVodesInternal*)user_data;
    this_->bjac(N, mupper, mlower, t, x, xdot, Jac, tmp1, tmp2, tmp3);
    return 0;
  } catch(exception& e){
    cerr << "bjac failed: " << e.what() << endl;;
    return 1;
  }
}

void CVodesInternal::bjac(SUNDIALS_INT N, SUNDIALS_INT mupper, SUNDIALS_INT mlower, double t, N_Vector x, N_Vector xdot, DlsMat Jac, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3){
  // Get time
  time1 = clock();

  // Pass inputs to the jacobian function
  jac_f_.setInput(&t,DAE_T);
  jac_f_.setInput(NV_DATA_S(x),DAE_X);
  jac_f_.setInput(f_.input(DAE_P),DAE_P);

  // Evaluate
  jac_f_.evaluate();
  
  // Get sparsity and non-zero elements
  const vector<int>& rowind = jac_f_.output().rowind();
  const vector<int>& col = jac_f_.output().col();
  const vector<double>& val = jac_f_.output().data();

  // Loop over rows
  for(int i=0; i<rowind.size()-1; ++i){
    // Loop over non-zero entries
    for(int el=rowind[i]; el<rowind[i+1]; ++el){
      // Get column
      int j = col[el];
      
      // Set the element
      if(i-j>=-mupper && i-j<=mlower)
        BAND_ELEM(Jac,i,j) = val[el];
    }
  }
  
  // Log time duration
  time2 = clock();
  t_jac += double(time2-time1)/CLOCKS_PER_SEC;
}

void CVodesInternal::setStopTime(double tf){
  // Set the stop time of the integration -- don't integrate past this point
  int flag = CVodeSetStopTime(mem_, tf);
  if(flag != CV_SUCCESS) cvodes_error("CVodeSetStopTime",flag);
}

int CVodesInternal::psolve_wrapper(double t, N_Vector x, N_Vector xdot, N_Vector r, N_Vector z, double gamma, double delta, int lr, void *user_data, N_Vector tmp){
  try{
    CVodesInternal *this_ = (CVodesInternal*)user_data;
    casadi_assert(this_);
    this_->psolve(t, x, xdot, r, z, gamma, delta, lr, tmp);
    return 0;
  } catch(exception& e){
    cerr << "psolve failed: " << e.what() << endl;;
    return 1;
  }
}

int CVodesInternal::psetup_wrapper(double t, N_Vector x, N_Vector xdot, booleantype jok, booleantype *jcurPtr, double gamma, void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3){
  try{
    CVodesInternal *this_ = (CVodesInternal*)user_data;
    casadi_assert(this_);
    this_->psetup(t, x, xdot, jok, jcurPtr, gamma, tmp1, tmp2, tmp3);
    return 0;
  } catch(exception& e){
    cerr << "psetup failed: " << e.what() << endl;;
    return 1;
  }
}

void CVodesInternal::psolve(double t, N_Vector x, N_Vector xdot, N_Vector r, N_Vector z, double gamma, double delta, int lr, N_Vector tmp){
  // Get time
  time1 = clock();

  // Copy input to output, if necessary
  if(r!=z){
    N_VScale(1.0, r, z);
  }

  // Solve the (possibly factorized) system 
  casadi_assert(linsol_.output().size()*nrhs_ == NV_LENGTH_S(z));
  linsol_.solve(NV_DATA_S(z),nrhs_);
  
  // Log time duration
  time2 = clock();
  t_lsolve += double(time2-time1)/CLOCKS_PER_SEC;
}

void CVodesInternal::psetup(double t, N_Vector x, N_Vector xdot, booleantype jok, booleantype *jcurPtr, double gamma, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3){
  // Get time
  time1 = clock();

  // Pass input to the jacobian function
  jac_.setInput(t,M_T);
  jac_.setInput(NV_DATA_S(x),M_Y);
  jac_.setInput(input(INTEGRATOR_P),M_P);
  jac_.setInput(gamma,M_GAMMA);

  // Evaluate jacobian
  jac_.evaluate();
  
  // Log time duration
  time2 = clock();
  t_lsetup_jac += double(time2-time1)/CLOCKS_PER_SEC;

  // Pass non-zero elements, scaled by -gamma, to the linear solver
  linsol_.setInput(jac_.output(),0);

  // Prepare the solution of the linear system (e.g. factorize) -- only if the linear solver inherits from LinearSolver
  linsol_.prepare();

  // Log time duration
  time1 = clock();
  t_lsetup_fac += double(time1-time2)/CLOCKS_PER_SEC;
}

void CVodesInternal::lsetup(CVodeMem cv_mem, int convfail, N_Vector x, N_Vector xdot, booleantype *jcurPtr, N_Vector vtemp1, N_Vector vtemp2, N_Vector vtemp3){
  // Current time
  double t = cv_mem->cv_tn;

  // Scaling factor before J
  double gamma = cv_mem->cv_gamma;

  // Call the preconditioner setup function (which sets up the linear solver)
  psetup(t, x, xdot, FALSE, jcurPtr, gamma, vtemp1, vtemp2, vtemp3);
}

int CVodesInternal::lsetup_wrapper(CVodeMem cv_mem, int convfail, N_Vector x, N_Vector xdot, booleantype *jcurPtr, N_Vector vtemp1, N_Vector vtemp2, N_Vector vtemp3){
  try{
    CVodesInternal *this_ = (CVodesInternal*)(cv_mem->cv_lmem);
    casadi_assert(this_);
    this_->lsetup(cv_mem, convfail, x, xdot, jcurPtr, vtemp1, vtemp2, vtemp3);
    return 0;
  } catch(exception& e){
    cerr << "lsetup failed: " << e.what() << endl;;
    return 1;
  }
}

void CVodesInternal::lsolve(CVodeMem cv_mem, N_Vector b, N_Vector weight, N_Vector x, N_Vector xdot){
  // Current time
  double t = cv_mem->cv_tn;

  // Scaling factor before J
  double gamma = cv_mem->cv_gamma;

  // Accuracy
  double delta = 0.0;
  
  // Left/right preconditioner
  int lr = 1;
  
  // Call the preconditioner solve function (which solves the linear system)
  psolve(t, x, xdot, b, b, gamma, delta, lr, 0);
}

int CVodesInternal::lsolve_wrapper(CVodeMem cv_mem, N_Vector b, N_Vector weight, N_Vector x, N_Vector xdot){
  try{
    CVodesInternal *this_ = (CVodesInternal*)(cv_mem->cv_lmem);
    casadi_assert(this_);
    this_->lsolve(cv_mem, b, weight, x, xdot);
    return 0;
  } catch(exception& e){
    cerr << "lsolve failed: " << e.what() << endl;;
    return 1;
  }
}

void CVodesInternal::initUserDefinedLinearSolver(){
  // Make sure that a Jacobian has been provided
  if(jac_.isNull()) throw CasadiException("CVodesInternal::initUserDefinedLinearSolver(): No Jacobian has been provided.");

  // Make sure that a linear solver has been providided
  if(linsol_.isNull()) throw CasadiException("CVodesInternal::initUserDefinedLinearSolver(): No user defined linear solver has been provided.");

  //  Set fields in the IDA memory
  CVodeMem cv_mem = CVodeMem(mem_);
  cv_mem->cv_lmem   = this;
  cv_mem->cv_lsetup = lsetup_wrapper;
  cv_mem->cv_lsolve = lsolve_wrapper;
  cv_mem->cv_setupNonNull = TRUE;
}

void CVodesInternal::setLinearSolver(const LinearSolver& linsol, const FX& jac){
  linsol_ = linsol;
  jac_ = jac;
}

FX CVodesInternal::getJacobian(){
  return jac_;
}
  
LinearSolver CVodesInternal::getLinearSolver(){
  return linsol_;
}


void CVodesInternal::deepCopyMembers(std::map<SharedObjectNode*,SharedObject>& already_copied){
  SundialsInternal::deepCopyMembers(already_copied);
  jac_f_ = deepcopy(jac_f_,already_copied);
}


} // namespace Sundials
} // namespace CasADi


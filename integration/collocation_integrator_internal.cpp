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

#include "collocation_integrator_internal.hpp"
#include "casadi/stl_vector_tools.hpp"
#include "casadi/matrix/sparsity_tools.hpp"
#include "casadi/matrix/matrix_tools.hpp"
#include "casadi/sx/sx_tools.hpp"
#include "casadi/fx/sx_function.hpp"
#include "casadi/mx/mx_tools.hpp"

using namespace std;
namespace CasADi{

CollocationIntegratorInternal::CollocationIntegratorInternal(const FX& f, const FX& q) : IntegratorInternal(f,q){
  addOption("number_of_finite_elements",     OT_INTEGER,  20, "Number of finite elements");
  addOption("interpolation_order",           OT_INTEGER,  3,  "Order of the interpolating polynomials");
  addOption("collocation_scheme",            OT_STRING,  "radau",  "Collocation scheme (radau or legendre)");
  addOption("implicit_solver",               OT_IMPLICITFUNCTION,  GenericType(), "An implicit function solver");
  addOption("implicit_solver_options",       OT_DICTIONARY, GenericType(), "Options to be passed to the NLP Solver");
  addOption("expand",                        OT_BOOLEAN,  false, "Expand the nonlinear system of equation in an SX graph");
  addOption("hotstart",                      OT_BOOLEAN,  true, "Initialize the trajectory at the previous solution");
}

void CollocationIntegratorInternal::deepCopyMembers(std::map<SharedObjectNode*,SharedObject>& already_copied){
}

CollocationIntegratorInternal::~CollocationIntegratorInternal(){
}

void CollocationIntegratorInternal::init(){
  
  // Init ODE rhs function and quadrature functions, jacobian function
  if(!f_.isInit()) f_.init();
  if(!q_.isNull() && !q_.isInit()) q_.init();
  
  log("CollocationIntegratorInternal::init","functions initialized");
  
  // Check dimensions
  casadi_assert_message(f_.getNumInputs()==DAE_NUM_IN, "CollocationIntegratorInternal: f has wrong number of inputs");
  casadi_assert_message(f_.getNumOutputs()==DAE_NUM_OUT, "CollocationIntegratorInternal: f has wrong number of outputs");
  if(!q_.isNull()){
    casadi_assert_message(q_.getNumInputs()==DAE_NUM_IN, "CollocationIntegratorInternal: q has wrong number of inputs");
    casadi_assert_message(q_.getNumOutputs()==DAE_NUM_OUT, "CollocationIntegratorInternal: q has wrong number of outputs");
  }

  ny_ = f_.input(DAE_Y).numel();
  nq_ = q_.isNull() ? 0 : q_.output().numel();
  
  int np = f_.input(DAE_P).numel();
  setDimensions(ny_+nq_,np);

  // Call the base class init
  IntegratorInternal::init();
  
  // Legendre collocation points
  double legendre_points[][6] = {
    {0},
    {0,0.500000},
    {0,0.211325,0.788675},
    {0,0.112702,0.500000,0.887298},
    {0,0.069432,0.330009,0.669991,0.930568},
    {0,0.046910,0.230765,0.500000,0.769235,0.953090}};
    
  // Radau collocation points
  double radau_points[][6] = {
    {0},
    {0,1.000000},
    {0,0.333333,1.000000},
    {0,0.155051,0.644949,1.000000},
    {0,0.088588,0.409467,0.787659,1.000000},
    {0,0.057104,0.276843,0.583590,0.860240,1.000000}};

  // Read options
  bool use_radau;
  if(getOption("collocation_scheme")=="radau"){
    use_radau = true;
  } else if(getOption("collocation_scheme")=="legendre"){
    use_radau = false;
  } else {
    throw CasadiException("Option \"collocation_scheme\" must be either \"radau\" or \"legendre\"");
  }
  
  // Hotstart?
  hotstart_ = getOption("hotstart");
  
  // Number of finite elements
  int nk = getOption("number_of_finite_elements");
  
  // Interpolation order
  int deg = getOption("interpolation_order");

  // All collocation time points
  double* tau_root = use_radau ? radau_points[deg] : legendre_points[deg];

  // Size of the finite elements
  double h = (tf_-t0_)/nk;
  
  // MX version of the same
  MX h_mx = h;
    
  // Coefficients of the collocation equation
  vector<vector<MX> > C(deg+1,vector<MX>(deg+1));

  // Coefficients of the continuity equation
  vector<MX> D(deg+1);

  // Collocation point
  SXMatrix tau = ssym("tau");

  // For all collocation points
  for(int j=0; j<deg+1; ++j){
    // Construct Lagrange polynomials to get the polynomial basis at the collocation point
    SXMatrix L = 1;
    for(int j2=0; j2<deg+1; ++j2){
      if(j2 != j){
        L *= (tau-tau_root[j2])/(tau_root[j]-tau_root[j2]);
      }
    }
    
    SXFunction lfcn(tau,L);
    lfcn.init();
  
    // Evaluate the polynomial at the final time to get the coefficients of the continuity equation
    lfcn.setInput(1.0);
    lfcn.evaluate();
    D[j] = lfcn.output();

    // Evaluate the time derivative of the polynomial at all collocation points to get the coefficients of the continuity equation
    for(int j2=0; j2<deg+1; ++j2){
      lfcn.setInput(tau_root[j2]);
      lfcn.setFwdSeed(1.0);
      lfcn.evaluate(1,0);
      C[j][j2] = lfcn.fwdSens();
    }
  }
  
  // Free parameters
  MX P("P",np_);

  // Initial state
  MX Y0("Y0",ny_);
  
  // Collocated states
  int nY = nk*(deg+1)*ny_ + ny_;
  
  // Unknowns
  MX V("V",nY-ny_);
  int offset = 0;
  
  // Get collocated states
  vector<vector<MX> > Y(nk+1);
  for(int k=0; k<nk; ++k){
    Y[k].resize(deg+1);

    // Get the expression for the state vector
    for(int j=0; j<deg+1; ++j){
      
      // If it's the first interval
      if(k==0 && j==0){
        Y[k][j] = Y0;
      } else {
        Y[k][j] = V[range(offset,offset+ny_)];
        offset += ny_;
      }
    }
  }
  
  // State at end time
  Y[nk].resize(1);
  Y[nk][0] = V[range(offset,offset+ny_)];
  offset += ny_;
  
  // Check offset for consistency
  casadi_assert(offset==V.size());

  // Constraints
  vector<MX> g;
  g.reserve(nk);
  
  // For all finite elements
  for(int k=0; k<nk; ++k){
  
    // For all collocation points
    for(int j=1; j<deg+1; ++j){
      // Get the time
      MX tk = h*(k + tau_root[j]);
        
      // Get an expression for the state derivative at the collocation point
      MX xp_jk = 0;
      for(int j2=0; j2<deg+1; ++j2){
        xp_jk += C[j2][j]*Y[k][j2];
      }
      
      // Add collocation equations to the NLP
      vector<MX> f_in(DAE_NUM_IN);
      f_in[DAE_T] = tk;
      f_in[DAE_P] = P;
      f_in[DAE_Y] = Y[k][j];
      vector<MX> f_out = f_.call(f_in);
      g.push_back(h_mx*f_out[DAE_RES] - xp_jk);
    }
    
    // Get an expression for the state at the end of the finite element
    MX yf_k = 0;
    for(int j=0; j<deg+1; ++j){
      yf_k += D[j]*Y[k][j];
    }

    // Add continuity equation to NLP
    g.push_back(Y[k+1][0] - yf_k);
  }
  
  // Constraint expression
  MX gv = vertcat(g);
  
  // Make sure that the dimension is consistent with the number of unknowns
  casadi_assert_message(gv.size()==V.size(),"Implicit function unknowns and equations do not match");

  // Nonlinear constraint function input
  vector<MX> gfcn_in(1+2);
  gfcn_in[0] = V;
  gfcn_in[1+0] = Y0;
  gfcn_in[1+1] = P;
  
  // Nonlinear constraint function
  MXFunction gfcn(gfcn_in,gv);
  
  // Expand?
  bool expand = getOption("expand");
  if(expand){
    gfcn.init();
    gfcn_ = SXFunction(gfcn);
  } else {
    gfcn_ = gfcn;
  }

  // Get the NLP creator function
  implicitFunctionCreator implicit_function_creator = getOption("implicit_solver");
  
  // Allocate an NLP solver
  implicit_solver_ = implicit_function_creator(gfcn_);
  
  // Pass options
  if(hasSetOption("implicit_solver_options")){
    const Dictionary& implicit_solver_options = getOption("implicit_solver_options");
    implicit_solver_.setOption(implicit_solver_options);
  }
  
  // Initialize the solver
  implicit_solver_.init();
  
  // Mark the system not yet integrated
  integrated_once_ = false;
}
  
void CollocationIntegratorInternal::initAdj(){
}

void CollocationIntegratorInternal::reset(int nfdir, int nadir){
  // Store the sensitivity directions
  nfdir_ = nfdir;
  nadir_ = nadir;
  
  // Pass the inputs
  const vector<double>& x0 = input(INTEGRATOR_X0).data();
  implicit_solver_.input(0).setArray(&x0.front(),ny_); // only the ny_ first elements
  implicit_solver_.input(1).set(input(INTEGRATOR_P));
  
  for(int dir=0; dir<nfdir_; ++dir){
    // Pass the forward seeds
    const vector<double>& x0 = fwdSeed(INTEGRATOR_X0,dir).data();
    implicit_solver_.fwdSeed(0,dir).setArray(&x0.front(),ny_); // only the ny_ first elements
    implicit_solver_.fwdSeed(1,dir).set(fwdSeed(INTEGRATOR_P,dir));
  }

  for(int dir=0; dir<nadir_; ++dir){
    // Pass the adjoint seeds
    vector<double>& v = implicit_solver_.adjSeed(0,dir).data();
    const vector<double>& xf = adjSeed(INTEGRATOR_XF,dir).data();

    // Set seed to zero
    fill(v.begin(),v.end(),0.0);
    
    // Pass seed for final state
    for(int i=0; i<ny_; ++i){
      v[v.size()-ny_+i] = xf[i];
    }
  }
  
  // Pass solution guess (if this is the first integration or if hotstart is disabled)
  if(!hotstart_ || !integrated_once_){
    vector<double>& v = implicit_solver_.output().data();
    for(int offs=0; offs<v.size(); offs+=ny_){
      for(int i=0; i<ny_; ++i){
        v[i+offs] = x0[i];
      }
    }
  }
    
  // Solve the system of equations
  implicit_solver_.evaluate(nfdir_, nadir_);
  
  // Mark the system integrated at least once
  integrated_once_ = true;
}

void CollocationIntegratorInternal::resetAdj(){
}

void CollocationIntegratorInternal::integrate(double t_out){
  const vector<double>& v = implicit_solver_.output().data();
  vector<double>& xf = output(INTEGRATOR_XF).data();
  
  for(int i=0; i<ny_; ++i){
    xf[i] = v[v.size()-ny_+i];
  }
  for(int i=ny_; i<nx_; ++i){
    xf[i] = 0;
  }

  for(int dir=0; dir<nfdir_; ++dir){
    // Get the forward sensitivities
    const vector<double>& v = implicit_solver_.fwdSens(0,dir).data();
    vector<double>& xf = fwdSens(INTEGRATOR_XF,dir).data();
    
    for(int i=0; i<ny_; ++i){
      xf[i] = v[v.size()-ny_+i];
    }
    for(int i=ny_; i<nx_; ++i){
      xf[i] = 0;
    }
  }

  for(int dir=0; dir<nadir_; ++dir){
    // Get the adjoint sensitivities
    vector<double>& x0 = adjSens(INTEGRATOR_X0,dir).data();
    implicit_solver_.adjSens(0,dir).getArray(&x0.front(),ny_); // only the ny_ first elements
    implicit_solver_.adjSens(1,dir).get(adjSens(INTEGRATOR_P,dir));
    
    for(int i=ny_; i<nx_; ++i){
      x0[i] = 0;
    }
  }
}

void CollocationIntegratorInternal::integrateAdj(double t_out){
}

FX CollocationIntegratorInternal::getJacobian(){
  return FX();
}
  
LinearSolver CollocationIntegratorInternal::getLinearSolver(){
  return LinearSolver();
}

void CollocationIntegratorInternal::setLinearSolver(const LinearSolver& linsol, const FX& jac){
}

void CollocationIntegratorInternal::printStats(std::ostream &stream) const{
}

void CollocationIntegratorInternal::setStopTime(double tf){
}

} // namespace CasADi

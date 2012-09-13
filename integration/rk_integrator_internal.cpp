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

#include "rk_integrator_internal.hpp"
#include "casadi/stl_vector_tools.hpp"
#include "casadi/matrix/sparsity_tools.hpp"
#include "casadi/matrix/matrix_tools.hpp"
#include "casadi/sx/sx_tools.hpp"
#include "casadi/fx/sx_function.hpp"
#include "casadi/mx/mx_tools.hpp"

using namespace std;
namespace CasADi{

RKIntegratorInternal::RKIntegratorInternal(const FX& f, const FX& g) : IntegratorInternal(f,g){
  addOption("number_of_finite_elements",     OT_INTEGER,  20, "Number of finite elements");
  addOption("interpolation_order",           OT_INTEGER,  4,  "Order of the interpolating polynomials");
  addOption("expand_f",                      OT_BOOLEAN,  false, "Expand the ODE/DAE residual function in an SX graph");
  addOption("expand_q",                      OT_BOOLEAN,  false, "Expand the quadrature function in an SX graph");
}

void RKIntegratorInternal::deepCopyMembers(std::map<SharedObjectNode*,SharedObject>& already_copied){
}

RKIntegratorInternal::~RKIntegratorInternal(){
}

void RKIntegratorInternal::init(){
  // Call the base class init
  IntegratorInternal::init();
  casadi_assert_message(nq_==0, "Quadratures not supported.");
  
  // Number of finite elements
  int nk = getOption("number_of_finite_elements");
  
  // Interpolation order
  int deg = getOption("interpolation_order");

  // Assume explicit ODE
  bool explicit_ode = f_.input(DAE_XDOT).size()==0;

  // Expand f?
  bool expand_f = getOption("expand_f");

  // Size of the finite elements
  double h = (tf_-t0_)/nk;
  
  // MX version of the same
  MX h_mx = h;
    
  // Initial state
  MX Y0("Y0",nx_);
  
  // Free parameters
  MX P("P",np_);

  // Current state
  MX Y = Y0;
  
  // Dummy time
  MX T = 0;
  
  // Zero YDOT argument to get an explicit ode
  MX YDOT = MX::zeros(Y.shape());
  
  // Integrate until the end
  for(int k=0; k<nk; ++k){
    
    // Call the ode right hand side function
    vector<MX> f_in(DAE_NUM_IN);
    f_in[DAE_T] = T;
    f_in[DAE_X] = Y;
    f_in[DAE_XDOT] = YDOT;
    f_in[DAE_P] = P;
    vector<MX> f_out = f_.call(f_in);
    MX ode_rhs = f_out[DAE_ODE];
    
    // Explicit Euler step
    Y += h_mx*ode_rhs;
  }
  
  // Create a function which returns the state at the end of the time horizon
  vector<MX> yf_in(2);
  yf_in[0] = Y0;
  yf_in[1] = P;
  MXFunction yf_fun(yf_in,Y);
  
  // Should the function be expanded in elementary operations?
  if(expand_f){
    yf_fun.init();
    yf_fun_ = SXFunction(yf_fun);
  } else {
    yf_fun_ = yf_fun;
  }
  
  // Set number of derivative directions
  yf_fun_.setOption("number_of_fwd_dir",getOption("number_of_fwd_dir"));
  yf_fun_.setOption("number_of_adj_dir",getOption("number_of_adj_dir"));
  
  // Initialize function
  yf_fun_.init();
}
  
void RKIntegratorInternal::initAdj(){
}

void RKIntegratorInternal::evaluate(int nfdir, int nadir){
  // Store the sensitivity directions
  nfdir_ = nfdir;
  nadir_ = nadir;
  
  // Call the base class method
  IntegratorInternal::evaluate(nfdir,nadir);
}

void RKIntegratorInternal::reset(int nsens, int nsensB, int nsensB_store){
  int nfdir = 0; // NOTE: need to update the function below to the new integrator formulation
  
  // Pass the inputs
  yf_fun_.setInput(input(INTEGRATOR_X0),0);
  yf_fun_.setInput(input(INTEGRATOR_P),1);
  
  // Pass the forward seeds
  for(int dir=0; dir<nfdir_; ++dir){
    yf_fun_.setFwdSeed(fwdSeed(INTEGRATOR_X0,dir),0,dir);
    yf_fun_.setFwdSeed(fwdSeed(INTEGRATOR_P,dir),1,dir);
  }
  
  // Pass the adjoint seeds
  for(int dir=0; dir<nadir_; ++dir){
    yf_fun_.setAdjSeed(adjSeed(INTEGRATOR_XF,dir),0,dir);
  }
  
  // Evaluate the function
  yf_fun_.evaluate(nfdir);
  
  // Get the outputs
  yf_fun_.getOutput(output(INTEGRATOR_XF),0);
  
  // Get the forward sensitivities
  for(int dir=0; dir<nfdir_; ++dir){
    yf_fun_.getFwdSens(fwdSens(INTEGRATOR_XF,dir),0,dir);
  }
  
  // Get the adjoint sensitivities
  for(int dir=0; dir<nadir_; ++dir){
    yf_fun_.getAdjSens(adjSens(INTEGRATOR_X0,dir),0,dir);
    yf_fun_.getAdjSens(adjSens(INTEGRATOR_P,dir),1,dir);
  }
}

void RKIntegratorInternal::resetB(){
}

void RKIntegratorInternal::integrate(double t_out){
}

void RKIntegratorInternal::integrateB(double t_out){
}

FX RKIntegratorInternal::getJacobian(){
  return FX();
}
  
LinearSolver RKIntegratorInternal::getLinearSolver(){
  return LinearSolver();
}

void RKIntegratorInternal::setLinearSolver(const LinearSolver& linsol, const FX& jac){
}

void RKIntegratorInternal::printStats(std::ostream &stream) const{
}

void RKIntegratorInternal::setStopTime(double tf){
}

FX RKIntegratorInternal::jacobian(const std::vector<std::pair<int,int> >& jblocks){
  return yf_fun_.jacobian(jblocks);
}

CRSSparsity RKIntegratorInternal::getJacSparsity(int iind, int oind){
  return yf_fun_.jacSparsity(iind, oind, true);
}

} // namespace CasADi

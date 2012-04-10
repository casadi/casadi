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

#include "integrator_internal.hpp"
#include <cassert>
#include "../stl_vector_tools.hpp"
#include "jacobian.hpp"
#include "../matrix/matrix_tools.hpp"
#include "../mx/mx_tools.hpp"
#include "../sx/sx_tools.hpp"
#include "mx_function.hpp"
#include "sx_function.hpp"

INPUTSCHEME(IntegratorInput)
OUTPUTSCHEME(IntegratorOutput)

using namespace std;
namespace CasADi{

IntegratorInternal::IntegratorInternal(const FX& fd, const FX& fq) : fd_(fd), fq_(fq){
  new_design_ = false;
  ctorInit();
}

IntegratorInternal::IntegratorInternal(const FX& f, const FX& g, const FX& h) : f_(f), g_(g), h_(h){
  new_design_ = true;
  ctorInit();
}

void IntegratorInternal::ctorInit(){
  // set default options
  setOption("name","unnamed_integrator"); // name of the function 
  
  // Additional options
  addOption("print_stats",                 OT_BOOLEAN,  false, "Print out statistics after integration");
  addOption("nrhs",                        OT_INTEGER, 1); // number of right hand sides
  addOption("t0",                          OT_REAL, 0.0); // start of the integration
  addOption("tf",                          OT_REAL, 1.0); // end of the integration
  
  // Negative number of parameters for consistancy checking
  np_ = -1;
}

IntegratorInternal::~IntegratorInternal(){ 
}

void IntegratorInternal::setDimensions(int nx, int np){
  nx_ = nx;
  np_ = np;
  
  // Allocate space for inputs
  input_.resize(INTEGRATOR_NUM_IN);
  input(INTEGRATOR_X0)  = DMatrix(nx_,1,0); // initial state value
  input(INTEGRATOR_XP0) = DMatrix(nx_,1,0); // initial state derivative value
  input(INTEGRATOR_P)   = DMatrix(np_,1,0); // parameter
  
  // Allocate space for outputs
  output_.resize(INTEGRATOR_NUM_OUT);
  output(INTEGRATOR_XF) = DMatrix(nx_,1,0);
  output(INTEGRATOR_XPF)= DMatrix(nx_,1,0);
}

void IntegratorInternal::evaluate(int nfdir, int nadir){
  
  // Reset solver
  reset(nfdir, nadir);

  // Integrate forward to the end of the time horizon
  integrate(tf_);

  // If backwards integration is needed
  if(nadir>0){
    
    // Re-initialize backward problem
    resetAdj();

    // Integrate backwards to the beginning
    integrateAdj(t0_);
  }
  
  // Print statistics
  if(getOption("print_stats")) printStats(std::cout);
}

void IntegratorInternal::init(){
  
  // Initialize the functions and get dimensions
  if(new_design_){
    
    // Initialize the functions
    casadi_assert(!f_.isNull());
    
    // Initialize, get and assert dimensions of the forward integration
    if(!f_.isInit()) f_.init();
    nx_ = f_.input(NEW_DAE_X).numel();
    nz_ = f_.input(NEW_DAE_Z).numel();
    np_  = f_.input(NEW_DAE_P).numel();
    nq_ = f_.output(NEW_DAE_QUAD).numel();
    casadi_assert_message(f_.output(NEW_DAE_ODE).numel()==nx_,"Inconsistent dimensions");
    casadi_assert_message(f_.output(NEW_DAE_ALG).numel()==nz_,"Inconsistent dimensions");
    
    // Make sure that both h and g are given, or neither
    casadi_assert_message(h_.isNull()==g_.isNull(),"Either both h and g should be given, or neither of them");
    if(h_.isNull()){
      nrx_ = 0;
      nrq_ = 0;
      nrz_ = 0;
    } else {
      // Initialize, get and assert dimensions of the terminal constraint function
      if(!h_.isInit()) h_.init();
      casadi_assert_message(h_.input(TERM_X).numel()==nx_,"Inconsistent dimensions");
      casadi_assert_message(h_.input(TERM_P).numel()==np_,"Inconsistent dimensions");
      nrx_ = h_.output(TERM_RX).numel();
      
      // Initialize and assert the dimensions of the backward integration
      if(!g_.isInit()) g_.init();
      nrz_ = g_.output(RDAE_ALG).numel();
      nrq_ = g_.output(RDAE_QUAD).numel();
      casadi_assert_message(g_.input(RDAE_X).numel()==nx_,"Inconsistent dimensions");
      casadi_assert_message(g_.input(RDAE_Z).numel()==nz_,"Inconsistent dimensions");
      casadi_assert_message(g_.input(RDAE_RX).numel()==nrx_,"Inconsistent dimensions");
      casadi_assert_message(g_.input(RDAE_RZ).numel()==nrz_,"Inconsistent dimensions");
      casadi_assert_message(g_.input(RDAE_P).numel()==np_,"Inconsistent dimensions");
      casadi_assert_message(g_.output(RDAE_ODE).numel()==nrx_,"Inconsistent dimensions");
    }
    
    // Allocate space for inputs
    input_.resize(NEW_INTEGRATOR_NUM_IN);
    input(NEW_INTEGRATOR_X0) = f_.output(NEW_DAE_ODE);
    input(NEW_INTEGRATOR_P) = f_.input(NEW_DAE_P);
  
    // Allocate space for outputs
    output_.resize(NEW_INTEGRATOR_NUM_OUT);
    output(NEW_INTEGRATOR_XF) = f_.output(NEW_DAE_ODE);
    output(NEW_INTEGRATOR_QF) = f_.output(NEW_DAE_QUAD);
    if(!g_.isNull()){
      output(NEW_INTEGRATOR_RX0) = g_.output(RDAE_ODE);
      output(NEW_INTEGRATOR_RQ0) = g_.output(RDAE_QUAD);
    }
  }
  
  // Make sure that the dimensions have been set
  casadi_assert_message(np_>=0, "\"setDimensions\" has not been called.");
  
  // Call the base class method
  FXInternal::init();

  // read options
  nrhs_ = getOption("nrhs");
  
  // Give an intial value for the time horizon
  t0_ = getOption("t0");
  tf_ = getOption("tf");
}

void IntegratorInternal::deepCopyMembers(std::map<SharedObjectNode*,SharedObject>& already_copied){
  FXInternal::deepCopyMembers(already_copied);
  if(new_design_){
    f_ = deepcopy(f_,already_copied);
    g_ = deepcopy(g_,already_copied);
    h_ = deepcopy(h_,already_copied);
  } else {
    fd_ = deepcopy(fd_,already_copied);
    fq_ = deepcopy(fq_,already_copied);
  }
}

} // namespace CasADi



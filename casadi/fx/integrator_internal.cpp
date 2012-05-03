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

IntegratorInternal::IntegratorInternal(const FX& f, const FX& g) : f_(f), g_(g){
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
  
  // Initialize the functions
  casadi_assert(!f_.isNull());
  
  // Initialize, get and assert dimensions of the forward integration
  if(!f_.isInit()) f_.init();
  casadi_assert_message(f_.getNumInputs()==DAE_NUM_IN,"Wrong number of inputs for the DAE callback function");
  casadi_assert_message(f_.getNumOutputs()==DAE_NUM_OUT,"Wrong number of outputs for the DAE callback function");
  nx_ = f_.input(DAE_X).numel();
  nz_ = f_.input(DAE_Z).numel();
  nq_ = f_.output(DAE_QUAD).numel();
  np_  = f_.input(DAE_P).numel();
  casadi_assert_message(f_.output(DAE_ODE).numel()==nx_,"Inconsistent dimensions");
  casadi_assert_message(f_.output(DAE_ALG).numel()==nz_,"Inconsistent dimensions");
  
  // Initialize, get and assert dimensions of the backwards integration
  if(g_.isNull()){
    // No backwards integration
    nrx_ = nrz_ = nrq_ = 0;
  } else {
    if(!g_.isInit()) g_.init();
    casadi_assert_message(g_.getNumInputs()==RDAE_NUM_IN,"Wrong number of inputs for the backwards DAE callback function");
    casadi_assert_message(g_.getNumOutputs()==RDAE_NUM_OUT,"Wrong number of outputs for the backwards DAE callback function");
    nrx_ = f_.input(RDAE_RX).numel();
    nrz_ = f_.input(RDAE_RZ).numel();
    nrq_ = f_.output(RDAE_QUAD).numel();
    casadi_assert_message(g_.input(RDAE_P).numel()==np_,"Inconsistent dimensions");
    casadi_assert_message(g_.input(RDAE_X).numel()==nx_,"Inconsistent dimensions");
    casadi_assert_message(g_.input(RDAE_Z).numel()==nz_,"Inconsistent dimensions");
    casadi_assert_message(g_.output(RDAE_ODE).numel()==nrx_,"Inconsistent dimensions");
    casadi_assert_message(g_.output(RDAE_ALG).numel()==nrz_,"Inconsistent dimensions");
  }
  
  // Allocate space for inputs
  input_.resize(INTEGRATOR_NUM_IN);
  input(INTEGRATOR_X0)  = f_.input(DAE_X);
  input(INTEGRATOR_P)   = f_.input(DAE_P);
  if(!g_.isNull()){
    input(INTEGRATOR_RX0)  = g_.input(RDAE_RX);
  }
  
  // Allocate space for outputs
  output_.resize(INTEGRATOR_NUM_OUT);
  output(INTEGRATOR_XF) = f_.output(DAE_ODE);
  output(INTEGRATOR_QF) = f_.output(DAE_QUAD);
  if(!g_.isNull()){
    output(INTEGRATOR_RXF)  = g_.output(RDAE_ODE);
    output(INTEGRATOR_RQF)  = g_.output(RDAE_QUAD);
  }
  
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
  f_ = deepcopy(f_,already_copied);
  g_ = deepcopy(g_,already_copied);
}

} // namespace CasADi



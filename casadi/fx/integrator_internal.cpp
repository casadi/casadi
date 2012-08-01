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
  casadi_assert_message(f_.input(DAE_X).dense(),"State vector must be dense in the DAE callback function");
  casadi_assert_message(f_.output(DAE_ODE).dense(),"Right hand side vector must be dense in the DAE callback function");
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

std::pair<FX,FX> IntegratorInternal::getAugmented(int nfwd, int nadj){
  // Currently only handles this special case
  casadi_assert(nfwd==0);
  casadi_assert(nadj==1);
  
  // Augmented DAE
  FX f_aug=f_;
  FX g_aug;
  
  // Handle this special case separately until the general case below is efficient enough
  SXFunction f = shared_cast<SXFunction>(f_);
  SXFunction g = shared_cast<SXFunction>(g_);
  if(f.isNull()==f_.isNull() && g.isNull()==g_.isNull()){
    vector<SXMatrix> arg = f.inputsSX();
    vector<SXMatrix> res = f.outputsSX();
    casadi_assert(arg.size()==DAE_NUM_IN);
    casadi_assert(res.size()==DAE_NUM_OUT);
    SXMatrix x = arg[DAE_X];
    SXMatrix z = arg[DAE_Z];
    SXMatrix p = arg[DAE_P];
    SXMatrix t = arg[DAE_T];
    SXMatrix xdot = arg[DAE_XDOT];
    
    SXMatrix rx = ssym("rx",res[DAE_ODE].sparsity());
    SXMatrix rz = ssym("rz",res[DAE_ALG].sparsity());
    SXMatrix rxdot = ssym("rxdot",res[DAE_ODE].sparsity());
    SXMatrix rp = ssym("rp",res[DAE_QUAD].sparsity());
    
    // Is the ODE part explicit?
    bool ode_is_explict = xdot.size()==0;
    
    // Number of adjoint sweeps needed
    int n_sweep = ode_is_explict ? 1 : 2;
    
    vector<vector<SXMatrix> > dummy;
    vector<vector<SXMatrix> > aseed(n_sweep,res);
    aseed[0][DAE_ODE] = rx;
    aseed[0][DAE_ALG] = rz;
    aseed[0][DAE_QUAD] = rp;
    if(!ode_is_explict){
      aseed[1][DAE_ODE] = -rxdot;
      aseed[1][DAE_ALG].setZero();
      aseed[1][DAE_QUAD].setZero();
    }
    vector<vector<SXMatrix> > asens(n_sweep,arg);
    f.evalSX(arg,res,dummy,dummy,aseed,asens,true);
    
    // Formulate the backwards integration problem
    vector<SXMatrix> rdae_in(RDAE_NUM_IN);
    rdae_in[RDAE_RX] = rx;
    rdae_in[RDAE_RZ] = rz;
    rdae_in[RDAE_RP] = rp;
    rdae_in[RDAE_X] = x;
    rdae_in[RDAE_Z] = z;
    rdae_in[RDAE_P] = p;
    rdae_in[RDAE_T] = t;
    rdae_in[RDAE_RXDOT] = rxdot;
    
    vector<SXMatrix> rdae_out(RDAE_NUM_OUT);
    rdae_out[RDAE_ODE] = asens[0][DAE_X];
    if(ode_is_explict){
      rdae_out[RDAE_ODE] += rxdot;
    } else {
      rdae_out[RDAE_ODE] += asens[1][DAE_XDOT];
    }
    rdae_out[RDAE_ALG] = asens[0][DAE_Z];
    rdae_out[RDAE_QUAD] = asens[0][DAE_P];
    
    g_aug = SXFunction(rdae_in,rdae_out);
  } else { // Not SXFunction

    // Formulate backwards integration problem
    MXFunction f = shared_cast<MXFunction>(f_);
    if(!f.isNull()){
      vector<MX> arg = f.inputsMX();
      vector<MX> res = f.outputsMX();
      casadi_assert(arg.size()==DAE_NUM_IN);
      casadi_assert(res.size()==DAE_NUM_OUT);
    
      MX x = arg[DAE_X];
      MX z = arg[DAE_Z];
      MX p = arg[DAE_P];
      MX t = arg[DAE_T];
      MX xdot = arg[DAE_XDOT];
    
      MX rx = msym("rx",res[DAE_ODE].sparsity());
      MX rz = msym("rz",res[DAE_ALG].sparsity());
      MX rxdot = msym("rxdot",res[DAE_ODE].sparsity());
      MX rp = msym("rp",res[DAE_QUAD].sparsity());
    
      // Is the ODE part explicit?
      bool ode_is_explict = xdot.size()==0;
    
      // Number of adjoint sweeps needed
      int n_sweep = ode_is_explict ? 1 : 2;
    
      vector<vector<MX> > dummy;
      vector<vector<MX> > aseed(n_sweep,res);
      aseed[0][DAE_ODE] = rx;
      aseed[0][DAE_ALG] = rz;
      aseed[0][DAE_QUAD] = rp;
      if(!ode_is_explict){
        aseed[1][DAE_ODE] = -rxdot;
        aseed[1][DAE_ALG] = MX::zeros(rz.size());
        aseed[1][DAE_QUAD] = MX::zeros(rp.size());
      }
      vector<vector<MX> > asens(n_sweep,arg);
      f.evalMX(arg,res,dummy,dummy,aseed,asens,true);
    
      // Formulate the backwards integration problem
      vector<MX> rdae_in(RDAE_NUM_IN);
      rdae_in[RDAE_RX] = rx;
      rdae_in[RDAE_RZ] = rz;
      rdae_in[RDAE_RP] = rp;
      rdae_in[RDAE_X] = x;
      rdae_in[RDAE_Z] = z;
      rdae_in[RDAE_P] = p;
      rdae_in[RDAE_T] = t;
      rdae_in[RDAE_RXDOT] = rxdot;
    
      vector<MX> rdae_out(RDAE_NUM_OUT);
      rdae_out[RDAE_ODE] = asens[0][DAE_X];
      if(ode_is_explict){
        rdae_out[RDAE_ODE] += rxdot;
      } else {
        rdae_out[RDAE_ODE] += asens[1][DAE_XDOT];
      }
      rdae_out[RDAE_ALG] = asens[0][DAE_Z];
      rdae_out[RDAE_QUAD] = asens[0][DAE_P];
      
      g_aug = MXFunction(rdae_in,rdae_out);
    }
  }
  
  return std::pair<FX,FX>(f_aug,g_aug);
}


FX IntegratorInternal::getDerivative(int nfwd, int nadj){
  // Generate augmented DAE
  std::pair<FX,FX> aug_dae = getAugmented(nfwd,nadj);
  
  // Create integrator for augmented DAE
  Integrator ret;
  ret.assignNode(create(aug_dae.first,aug_dae.second));
  return ret;
}


} // namespace CasADi



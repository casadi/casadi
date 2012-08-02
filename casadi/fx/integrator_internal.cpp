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
  // Adjoint derivatives are calculated with source code transformation
  if(nadir>0){
    
    // Get derivative function
    FX dfcn = derivative(0,nadir);
    
    // Pass function values
    dfcn.setInput(input(INTEGRATOR_X0),INTEGRATOR_X0);
    input(INTEGRATOR_P).get(dfcn.input(INTEGRATOR_P).ptr());
    adjSeed(INTEGRATOR_QF).get(dfcn.input(INTEGRATOR_P).ptr() + np_);
    
    // Pass forward seeds
    for(int dir=0; dir<nfdir; ++dir){
      fwdSeed(INTEGRATOR_X0,dir).get(dfcn.fwdSeed(INTEGRATOR_X0,dir));
      dfcn.fwdSeed(INTEGRATOR_P,dir).setZero();
      fwdSeed(INTEGRATOR_P,dir).get(dfcn.fwdSeed(INTEGRATOR_P,dir).ptr());
    }
    
    // Pass adjoint seeds
    casadi_assert(nadir==1);
    dfcn.setInput(adjSeed(INTEGRATOR_XF),INTEGRATOR_RX0);
    
    // Evaluate with forward sensitivities
    dfcn.evaluate(nfdir,0);
    
    // Get results
    dfcn.getOutput(output(INTEGRATOR_XF),INTEGRATOR_XF);
    
    // Get forward sensitivities
    for(int dir=0; dir<nfdir; ++dir){
      dfcn.getFwdSens(fwdSens(INTEGRATOR_XF,dir),INTEGRATOR_XF,dir);
      dfcn.getFwdSens(fwdSens(INTEGRATOR_QF,dir),INTEGRATOR_QF,dir);
    }
    
    // Get adjoint sensitivities
    casadi_assert(nadir==1);
    dfcn.getOutput(adjSens(INTEGRATOR_X0),INTEGRATOR_RXF);
    dfcn.getOutput(adjSens(INTEGRATOR_P),INTEGRATOR_RQF);
    
  } else {
  
    // Reset solver
    reset(nfdir);

    // Integrate forward to the end of the time horizon
    integrate(tf_);

    // If backwards integration is needed
    if(!g_.isNull()){
      
      // Re-initialize backward problem
      resetAdj();

      // Integrate backwards to the beginning
      integrateAdj(t0_);
    }
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
    nrx_ = g_.input(RDAE_RX).numel();
    nrz_ = g_.input(RDAE_RZ).numel();
    nrq_ = g_.output(RDAE_QUAD).numel();
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
  FX f_aug, g_aug;
  
  // Handle this special case separately until the general case below is efficient enough
  SXFunction f = shared_cast<SXFunction>(f_);
  SXFunction g = shared_cast<SXFunction>(g_);
  if(f.isNull()==f_.isNull() && g.isNull()==g_.isNull()){
    vector<SXMatrix> dae_in = f.inputsSX();
    vector<SXMatrix> dae_out = f.outputsSX();
    casadi_assert(dae_in.size()==DAE_NUM_IN);
    casadi_assert(dae_out.size()==DAE_NUM_OUT);
    SXMatrix x = dae_in[DAE_X];
    SXMatrix z = dae_in[DAE_Z];
    SXMatrix p = dae_in[DAE_P];
    SXMatrix t = dae_in[DAE_T];
    SXMatrix xdot = dae_in[DAE_XDOT];
    
    SXMatrix rx = ssym("rx",dae_out[DAE_ODE].sparsity());
    SXMatrix rz = ssym("rz",dae_out[DAE_ALG].sparsity());
    SXMatrix rxdot = ssym("rxdot",dae_out[DAE_ODE].sparsity());
    SXMatrix rp = ssym("rp",dae_out[DAE_QUAD].sparsity());
    
    // Is the ODE part explicit?
    bool ode_is_explict = xdot.size()==0;
    
    // Number of adjoint sweeps needed
    int n_sweep = ode_is_explict ? 1 : 2;
    
    vector<vector<SXMatrix> > dummy;
    vector<vector<SXMatrix> > aseed(n_sweep,dae_out);
    aseed[0][DAE_ODE] = rx;
    aseed[0][DAE_ALG] = rz;
    aseed[0][DAE_QUAD] = rp;
    if(!ode_is_explict){
      aseed[1][DAE_ODE] = -rxdot;
      aseed[1][DAE_ALG].setZero();
      aseed[1][DAE_QUAD].setZero();
    }
    vector<vector<SXMatrix> > asens(n_sweep,dae_in);
    f.evalSX(dae_in,dae_out,dummy,dummy,aseed,asens,true);
    
    // Augment parameter vector
    SXMatrix p_aug = vertcat(p,rp);
    dae_in[DAE_P] = p_aug;
    
    // Formulate the backwards integration problem
    vector<SXMatrix> rdae_in(RDAE_NUM_IN);
    rdae_in[RDAE_RX] = rx;
    rdae_in[RDAE_RZ] = rz;
    rdae_in[RDAE_X] = x;
    rdae_in[RDAE_Z] = z;
    rdae_in[RDAE_P] = p_aug;
    rdae_in[RDAE_T] = t;
    rdae_in[RDAE_XDOT] = xdot;
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
    
    f_aug = SXFunction(dae_in,dae_out);
    g_aug = SXFunction(rdae_in,rdae_out);
    
  } else { // Not SXFunction

    // Formulate backwards integration problem
    MXFunction f = shared_cast<MXFunction>(f_);
    if(!f.isNull()){
      vector<MX> dae_in = f.inputsMX();
      vector<MX> dae_out = f.outputsMX();
      casadi_assert(dae_in.size()==DAE_NUM_IN);
      casadi_assert(dae_out.size()==DAE_NUM_OUT);
    
      MX x = dae_in[DAE_X];
      MX z = dae_in[DAE_Z];
      MX p = dae_in[DAE_P];
      MX t = dae_in[DAE_T];
      MX xdot = dae_in[DAE_XDOT];
    
      MX rx = msym("rx",dae_out[DAE_ODE].sparsity());
      MX rz = msym("rz",dae_out[DAE_ALG].sparsity());
      MX rxdot = msym("rxdot",dae_out[DAE_ODE].sparsity());
      MX rp = msym("rp",dae_out[DAE_QUAD].sparsity());
    
      // Is the ODE part explicit?
      bool ode_is_explict = xdot.size()==0;
    
      // Number of adjoint sweeps needed
      int n_sweep = ode_is_explict ? 1 : 2;
    
      vector<vector<MX> > dummy;
      vector<vector<MX> > aseed(n_sweep,dae_out);
      aseed[0][DAE_ODE] = rx;
      aseed[0][DAE_ALG] = rz;
      aseed[0][DAE_QUAD] = rp;
      if(!ode_is_explict){
        aseed[1][DAE_ODE] = -rxdot;
        aseed[1][DAE_ALG] = MX::zeros(rz.size());
        aseed[1][DAE_QUAD] = MX::zeros(rp.size());
      }
      vector<vector<MX> > asens(n_sweep,dae_in);
      f.evalMX(dae_in,dae_out,dummy,dummy,aseed,asens,true);
    
     // Augment p
      MX p_aug("p_aug",np_+nrq_);
      
      // Variables to replace
      vector<MX> v(2);
      v[0] = p;
      v[1] = rp;
      
      // What to replace with
      vector<MX> vdef(2);
      vdef[0] = p_aug[Slice(0,np_)];
      vdef[1] = p_aug[Slice(np_,np_+nrq_)];
       
      // Replace in expressions
      dae_out = substitute(dae_out,v,vdef);
      
      // Write forward dae in terms of p_aug
      dae_in[DAE_P] = p_aug;
      f_aug = MXFunction(dae_in,dae_out);
      
      // Formulate the backwards integration problem, input...
      vector<MX> rdae_in(RDAE_NUM_IN);
      rdae_in[RDAE_RX] = rx;
      rdae_in[RDAE_RZ] = rz;
      rdae_in[RDAE_X] = x;
      rdae_in[RDAE_Z] = z;
      rdae_in[RDAE_P] = p_aug;
      rdae_in[RDAE_T] = t;
      rdae_in[RDAE_XDOT] = xdot;
      rdae_in[RDAE_RXDOT] = rxdot;
      
      // ... and output
      vector<MX> rdae_out(RDAE_NUM_OUT);
      rdae_out[RDAE_ODE] = asens[0][DAE_X];
      if(ode_is_explict){
        rdae_out[RDAE_ODE] += rxdot;
      } else {
        rdae_out[RDAE_ODE] += asens[1][DAE_XDOT];
      }
      rdae_out[RDAE_ALG] = asens[0][DAE_Z];
      rdae_out[RDAE_QUAD] = asens[0][DAE_P];
      
      // Replace p and rp in expressions
      rdae_out = substitute(rdae_out,v,vdef);
      
      // Create backward integration function
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
  
  // Copy options
  ret.setOption(dictionary());
  
  return ret;
}


} // namespace CasADi



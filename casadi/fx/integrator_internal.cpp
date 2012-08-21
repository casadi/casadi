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
  if(nfdir>0 || nadir==0){
  
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
  
  // Adjoint derivatives are calculated with source code transformation
  if(nadir>0){ 
    // NOTE: The following is a general functionality that should be moved to the base class
    
    // Get derivative function
    FX dfcn = derivative(0,nadir);

    // Pass function values
    for(int i=0; i<INTEGRATOR_NUM_IN; ++i){
      dfcn.setInput(input(i),i);
    }
    
    // Pass adjoint seeds
    for(int dir=0; dir<nadir; ++dir){
      for(int i=0; i<INTEGRATOR_NUM_OUT; ++i){
        dfcn.setInput(adjSeed(i,dir),INTEGRATOR_NUM_IN + dir*INTEGRATOR_NUM_OUT + i);
      }
    }
    
    // Evaluate to get function values and adjoint sensitivities
    dfcn.evaluate();
    
    // Get nondifferentiated results
    if(nfdir==0){
      for(int i=0; i<INTEGRATOR_NUM_OUT; ++i){
        dfcn.getOutput(output(i),i);
      }
    }
    
    // Get adjoint sensitivities 
    for(int dir=0; dir<nadir; ++dir){
      for(int i=0; i<INTEGRATOR_NUM_IN; ++i){
        dfcn.getOutput(adjSens(i,dir),INTEGRATOR_NUM_OUT + dir*INTEGRATOR_NUM_IN + i);
      }
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
  casadi_assert_message(f_.output(DAE_ODE).numel()==nx_,"Inconsistent dimensions. Expecting DAE_ODE output of size " << nx_ << ", but got " << f_.output(DAE_ODE).numel() << " instead.");
  casadi_assert_message(f_.output(DAE_ALG).numel()==nz_,"Inconsistent dimensions. Expecting DAE_ALG output of size " << nz_ << ", but got " << f_.output(DAE_ALG).numel() << " instead.");
  
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
    casadi_assert_message(g_.input(RDAE_P).numel()==np_,"Inconsistent dimensions. Expecting RDAE_P input of size " << np_ << ", but got " << g_.input(RDAE_P).numel() << " instead.");
    casadi_assert_message(g_.input(RDAE_X).numel()==nx_,"Inconsistent dimensions. Expecting RDAE_X input of size " << nx_ << ", but got " << g_.input(RDAE_P).numel() << " instead.");
    casadi_assert_message(g_.input(RDAE_Z).numel()==nz_,"Inconsistent dimensions. Expecting RDAE_Z input of size " << nz_ << ", but got " << g_.input(RDAE_P).numel() << " instead.");
    casadi_assert_message(g_.output(RDAE_ODE).numel()==nrx_,"Inconsistent dimensions. Expecting RDAE_ODE input of size " << nrx_ << ", but got " << g_.input(RDAE_P).numel() << " instead.");
    casadi_assert_message(g_.output(RDAE_ALG).numel()==nrz_,"Inconsistent dimensions. Expecting RDAE_ALG input of size " << nrz_ << ", but got " << g_.input(RDAE_P).numel() << " instead.");
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
    
    vector<SXMatrix> rx = ssym("rx",dae_out[DAE_ODE].sparsity(),nadj);
    vector<SXMatrix> rz = ssym("rz",dae_out[DAE_ALG].sparsity(),nadj);
    vector<SXMatrix> rxdot = ssym("rxdot",dae_out[DAE_ODE].sparsity(),nadj);
    vector<SXMatrix> rp = ssym("rp",dae_out[DAE_QUAD].sparsity(),nadj);
    
    // Is the ODE part explicit?
    bool ode_is_explict = xdot.size()==0;
    
    // Number of adjoint sweeps needed
    int n_sweep = nadj*(ode_is_explict ? 1 : 2);
    
    vector<vector<SXMatrix> > dummy;
    vector<vector<SXMatrix> > aseed(n_sweep,dae_out);
    for(int dir=0; dir<nadj; ++dir){
      aseed[dir][DAE_ODE] = rx[dir];
      aseed[dir][DAE_ALG] = rz[dir];
      aseed[dir][DAE_QUAD] = rp[dir];
      if(!ode_is_explict){
        aseed[nadj+dir][DAE_ODE] = -rxdot[dir];
        aseed[nadj+dir][DAE_ALG].setZero();
        aseed[nadj+dir][DAE_QUAD].setZero();
      }
    }
    vector<vector<SXMatrix> > asens(n_sweep,dae_in);
    f.evalSX(dae_in,dae_out,dummy,dummy,aseed,asens,true);
    
    // Augment parameter vector
    SXMatrix p_aug = vertcat(p,vertcat(rp));
    dae_in[DAE_P] = p_aug;
    
    // Formulate the backwards integration problem
    vector<SXMatrix> rdae_in(RDAE_NUM_IN);
    rdae_in[RDAE_RX] = rx.back();
    rdae_in[RDAE_RZ] = rz.back();
    rdae_in[RDAE_X] = x;
    rdae_in[RDAE_Z] = z;
    rdae_in[RDAE_P] = p_aug;
    rdae_in[RDAE_T] = t;
    rdae_in[RDAE_XDOT] = xdot;
    rdae_in[RDAE_RXDOT] = rxdot.back();
    
    vector<SXMatrix> rdae_out(RDAE_NUM_OUT);
    rdae_out[RDAE_ODE] = asens[0][DAE_X];
    if(ode_is_explict){
      rdae_out[RDAE_ODE] += rxdot.back();
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
    
      vector<MX> rx = msym("rx",dae_out[DAE_ODE].sparsity(),nadj);
      vector<MX> rz = msym("rz",dae_out[DAE_ALG].sparsity(),nadj);
      vector<MX> rxdot = msym("rxdot",dae_out[DAE_ODE].sparsity(),nadj);
      vector<MX> rp = msym("rp",dae_out[DAE_QUAD].sparsity(),nadj);
    
      // Is the ODE part explicit?
      bool ode_is_explict = xdot.size()==0;
    
      // Number of adjoint sweeps needed
      int n_sweep = ode_is_explict ? 1 : 2;
    
      vector<vector<MX> > dummy;
      vector<vector<MX> > aseed(n_sweep,dae_out);
      aseed[0][DAE_ODE] = rx.back();
      aseed[0][DAE_ALG] = rz.back();
      aseed[0][DAE_QUAD] = rp.back();
      if(!ode_is_explict){
        aseed[1][DAE_ODE] = -rxdot.back();
        aseed[1][DAE_ALG] = MX::zeros(rz.back().size());
        aseed[1][DAE_QUAD] = MX::zeros(rp.back().size());
      }
      vector<vector<MX> > asens(n_sweep,dae_in);
      f.evalMX(dae_in,dae_out,dummy,dummy,aseed,asens,true);
    
     // Augment p
      MX p_aug("p_aug",np_+nrq_);
      
      // Variables to replace
      vector<MX> v(2);
      v[0] = p;
      v[1] = rp.back();
      
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
      rdae_in[RDAE_RX] = rx.back();
      rdae_in[RDAE_RZ] = rz.back();
      rdae_in[RDAE_X] = x;
      rdae_in[RDAE_Z] = z;
      rdae_in[RDAE_P] = p_aug;
      rdae_in[RDAE_T] = t;
      rdae_in[RDAE_XDOT] = xdot;
      rdae_in[RDAE_RXDOT] = rxdot.back();
      
      // ... and output
      vector<MX> rdae_out(RDAE_NUM_OUT);
      rdae_out[RDAE_ODE] = asens[0][DAE_X];
      if(ode_is_explict){
        rdae_out[RDAE_ODE] += rxdot.back();
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
  Integrator integrator;
  integrator.assignNode(create(aug_dae.first,aug_dae.second));
  
  // Copy options
  integrator.setOption(dictionary());
  
  // Initialize the integrator since we will call it below
  integrator.init();
  
  // All inputs of the return function
  vector<MX> ret_in;
  ret_in.reserve(INTEGRATOR_NUM_IN*(1+nfwd) + INTEGRATOR_NUM_OUT*nadj);
  
  // Augmented state
  MX x0_aug, p_aug, rx0_aug;
  
  // Temp stringstream
  stringstream ss;
    
  // Inputs or forward/adjoint seeds in one direction
  vector<MX> dd;
  
  // Add nondifferentiated inputs and forward seeds
  dd.resize(INTEGRATOR_NUM_IN);
  for(int dir=-1; dir<nfwd; ++dir){
    
    // Differential state
    ss.clear();
    ss << "x0";
    if(dir>=0) ss << "_" << dir;
    dd[INTEGRATOR_X0] = msym(ss.str(),input(INTEGRATOR_X0).sparsity());
    x0_aug.append(dd[INTEGRATOR_X0]);

    // Parameter
    ss.clear();
    ss << "p";
    if(dir>=0) ss << "_" << dir;
    dd[INTEGRATOR_P] = msym(ss.str(),input(INTEGRATOR_P).sparsity());
    p_aug.append(dd[INTEGRATOR_P]);
    
    // Backward state
    ss.clear();
    ss << "rx0";
    if(dir>=0) ss << "_" << dir;
    dd[INTEGRATOR_RX0] = msym(ss.str(),input(INTEGRATOR_RX0).sparsity());
    rx0_aug.append(dd[INTEGRATOR_RX0]);
    
    // Add to input vector
    ret_in.insert(ret_in.end(),dd.begin(),dd.end());
  }
    
  // Add adjoint seeds
  dd.resize(INTEGRATOR_NUM_OUT);
  for(int dir=0; dir<nadj; ++dir){
    
    // Differential states become backward differential state
    ss.clear();
    ss << "xf" << "_" << dir;
    dd[INTEGRATOR_XF] = msym(ss.str(),output(INTEGRATOR_XF).sparsity());
    rx0_aug.append(dd[INTEGRATOR_XF]);

    // Quadratures become parameters
    ss.clear();
    ss << "qf" << "_" << dir;
    dd[INTEGRATOR_QF] = msym(ss.str(),output(INTEGRATOR_QF).sparsity());
    p_aug.append(dd[INTEGRATOR_QF]);

    // Backward differential states becomes forward differential states
    ss.clear();
    ss << "rxf" << "_" << dir;
    dd[INTEGRATOR_RXF] = msym(ss.str(),output(INTEGRATOR_RXF).sparsity());
    x0_aug.append(dd[INTEGRATOR_RXF]);
    
    // Backward quadratures becomes parameters
    ss.clear();
    ss << "rqf" << "_" << dir;
    dd[INTEGRATOR_RQF] = msym(ss.str(),output(INTEGRATOR_RQF).sparsity());
    p_aug.append(dd[INTEGRATOR_RQF]);
    
    // Add to input vector
    ret_in.insert(ret_in.end(),dd.begin(),dd.end());
  }
  
  // Call the integrator
  vector<MX> integrator_in(INTEGRATOR_NUM_IN);
  integrator_in[INTEGRATOR_X0] = x0_aug;
  integrator_in[INTEGRATOR_P] = p_aug;
  integrator_in[INTEGRATOR_RX0] = rx0_aug;
  vector<MX> integrator_out = integrator.call(integrator_in);
  
  // Augmented results
  MX xf_aug = integrator_out[INTEGRATOR_XF];
  MX rxf_aug = integrator_out[INTEGRATOR_RXF];
  MX qf_aug = integrator_out[INTEGRATOR_QF];
  MX rqf_aug = integrator_out[INTEGRATOR_RQF];
  
  // Offset in each of the above vectors
  int xf_offset = 0, rxf_offset = 0, qf_offset = 0, rqf_offset = 0;
  
  // All outputs of the return function
  vector<MX> ret_out;
  ret_out.reserve(INTEGRATOR_NUM_OUT*(1+nfwd) + INTEGRATOR_NUM_IN*nadj);
  
  // Collect the nondifferentiated results and forward sensitivities
  dd.resize(INTEGRATOR_NUM_OUT);
  fill(dd.begin(),dd.end(),MX());
  for(int dir=-1; dir<nfwd; ++dir){
    dd[INTEGRATOR_XF] = xf_aug[Slice(xf_offset,xf_offset+nx_)]; xf_offset += nx_;
    if(nq_>0){ dd[INTEGRATOR_QF] = qf_aug[Slice(qf_offset,qf_offset+nq_)]; qf_offset += nq_; }
    if(nrx_>0){ dd[INTEGRATOR_RXF] = rxf_aug[Slice(rxf_offset,rxf_offset+nrx_)]; rxf_offset += nrx_; }
    if(nrq_>0){ dd[INTEGRATOR_RQF] = rqf_aug[Slice(rqf_offset,rqf_offset+nrq_)]; rqf_offset += nrq_; }
    ret_out.insert(ret_out.end(),dd.begin(),dd.end());
  }
  
  // Collect the adjoint sensitivities
  dd.resize(INTEGRATOR_NUM_IN);
  fill(dd.begin(),dd.end(),MX());
  for(int dir=0; dir<nadj; ++dir){
    dd[INTEGRATOR_X0] = rxf_aug[Slice(rxf_offset,rxf_offset+nx_)]; rxf_offset += nx_;
    if(np_>0){ dd[INTEGRATOR_P] = rqf_aug[Slice(rqf_offset,rqf_offset+np_)]; rqf_offset += np_; }
    if(nrx_>0){ dd[INTEGRATOR_RX0] = xf_aug[Slice(xf_offset,xf_offset+nrx_)]; xf_offset += nrx_; }
    ret_out.insert(ret_out.end(),dd.begin(),dd.end());
  }
  
  // Create derivative function and return
  return MXFunction(ret_in,ret_out);
}


} // namespace CasADi



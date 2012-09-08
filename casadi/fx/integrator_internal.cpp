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
  if(nfdir==0 && nadir==0){
  
    // Reset solver
    reset();

    // Integrate forward to the end of the time horizon
    integrate(tf_);

    // If backwards integration is needed
    if(!g_.isNull()){
      
      // Re-initialize backward problem
      resetAdj();

      // Integrate backwards to the beginning
      integrateAdj(t0_);
    }
  } else {
    // NOTE: The following is a general functionality that should be moved to the base class
    
    // Get derivative function
    FX dfcn = derivative(nfdir,nadir);

    // Pass function values
    int input_index = 0;
    for(int i=0; i<INTEGRATOR_NUM_IN; ++i){
      dfcn.setInput(input(i),input_index++);
    }
    
    // Pass forward seeds
    for(int dir=0; dir<nfdir; ++dir){
      for(int i=0; i<INTEGRATOR_NUM_IN; ++i){
        dfcn.setInput(fwdSeed(i,dir),input_index++);
      }
    }
    
    // Pass adjoint seeds
    for(int dir=0; dir<nadir; ++dir){
      for(int i=0; i<INTEGRATOR_NUM_OUT; ++i){
        dfcn.setInput(adjSeed(i,dir),input_index++);
      }
    }
    
    // Evaluate to get function values and adjoint sensitivities
    dfcn.evaluate();
    
    // Get nondifferentiated results
    int output_index = 0;
    for(int i=0; i<INTEGRATOR_NUM_OUT; ++i){
      dfcn.getOutput(output(i),output_index++);
    }
    
    // Get forward sensitivities 
    for(int dir=0; dir<nfdir; ++dir){
      for(int i=0; i<INTEGRATOR_NUM_OUT; ++i){
        dfcn.getOutput(fwdSens(i,dir),output_index++);
      }
    }
    
    // Get adjoint sensitivities 
    for(int dir=0; dir<nadir; ++dir){
      for(int i=0; i<INTEGRATOR_NUM_IN; ++i){
        dfcn.getOutput(adjSens(i,dir),output_index++);
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
  if(is_a<SXFunction>(f_)){
    casadi_assert_message(g_.isNull() || is_a<SXFunction>(g_), "Currently, g_ must be of the same type as f_");
    return getAugmentedGen<SXMatrix,SXFunction>(nfwd,nadj);
  } else if(is_a<MXFunction>(f_)){
    casadi_assert_message(g_.isNull() || is_a<MXFunction>(g_), "Currently, g_ must be of the same type as f_");
    return getAugmentedGen<MX,MXFunction>(nfwd,nadj);
  } else {
    throw CasadiException("Currently, f_ must be either SXFunction or MXFunction");
  }
}
  
template<class Mat,class XFunc>
std::pair<FX,FX> IntegratorInternal::getAugmentedGen(int nfwd, int nadj){
    
  // Get derivatived type
  XFunc f = shared_cast<XFunc>(f_);
  XFunc g = shared_cast<XFunc>(g_);
  
  // Take apart forward problem
  vector<Mat> dae_in = f.inputExpr();
  vector<Mat> dae_out = f.outputExpr();
  casadi_assert(dae_in.size()==DAE_NUM_IN);
  casadi_assert(dae_out.size()==DAE_NUM_OUT);
  Mat x = dae_in[DAE_X];
  Mat z = dae_in[DAE_Z];
  Mat p = dae_in[DAE_P];
  Mat t = dae_in[DAE_T];
  Mat xdot = dae_in[DAE_XDOT];
  Mat ode = dae_out[DAE_ODE];
  Mat alg = dae_out[DAE_ALG];
  Mat quad = dae_out[DAE_QUAD];
  
  // Take apart the backwards problem
  vector<Mat> rdae_in(RDAE_NUM_IN), rdae_out(RDAE_NUM_OUT);
  if(!g.isNull()){
    rdae_in = g.inputExpr();
    rdae_out = g.outputExpr();
    // TODO: Assert that rdae_in[RDAE_X]==x, rdae_in[RDAE_Z]==z, rdae_in[RDAE_P]==p, rdae_in[RDAE_XDOT]==xdot
  }
  Mat rx = rdae_in[RDAE_RX];
  Mat rz = rdae_in[RDAE_RZ];
  Mat rxdot = rdae_in[RDAE_RXDOT];
  Mat rode = rdae_out[RDAE_ODE];
  Mat ralg = rdae_out[RDAE_ALG];
  Mat rquad = rdae_out[RDAE_QUAD];
  
  // Allocate adjoint sensitivities
  vector<Mat> adj_ode = Mat::sym("adj_ode",ode.sparsity(),nadj);
  vector<Mat> adj_alg = Mat::sym("adj_alg",alg.sparsity(),nadj);
  vector<Mat> adj_odedot = Mat::sym("adj_odedot",ode.sparsity(),nadj);
  vector<Mat> adj_quad = Mat::sym("adj_quad",quad.sparsity(),nadj);
  
  // Allocate forward sensitivities
  vector<Mat> fwd_x = Mat::sym("fwd_x",x.sparsity(),nfwd);
  vector<Mat> fwd_z = Mat::sym("fwd_z",z.sparsity(),nfwd);
  vector<Mat> fwd_xdot = Mat::sym("fwd_xdot",xdot.sparsity(),nfwd);
  vector<Mat> fwd_p = Mat::sym("fwd_p",p.sparsity(),nfwd);

  // Is the ODE part explicit?
  bool ode_is_explict = xdot.size()==0;
  
  // Forward seeds
  vector<vector<Mat> > fseed(nfwd,vector<Mat>(DAE_NUM_IN));
  for(int dir=0; dir<nfwd; ++dir){
    fseed[dir][DAE_X] = fwd_x[dir];
    fseed[dir][DAE_Z] = fwd_z[dir];
    fseed[dir][DAE_P] = fwd_p[dir];
    if(!t.isNull()) fseed[dir][DAE_T] = Mat(t.sparsity());
    fseed[dir][DAE_XDOT] = fwd_xdot[dir];
  }

  // Adjoint seeds
  vector<vector<Mat> > aseed(nadj*(ode_is_explict ? 1 : 2),vector<Mat>(DAE_NUM_OUT));
  for(int dir=0; dir<nadj; ++dir){
    aseed[dir][DAE_ODE] = adj_ode[dir];
    aseed[dir][DAE_ALG] = adj_alg[dir];
    aseed[dir][DAE_QUAD] = adj_quad[dir];
    if(!ode_is_explict){
      aseed[nadj+dir][DAE_ODE] = -adj_odedot[dir];
      aseed[nadj+dir][DAE_ALG] = Mat(dae_out[DAE_ALG].sparsity());
      aseed[nadj+dir][DAE_QUAD] = Mat(dae_out[DAE_QUAD].sparsity());
    }
  }
  
  // Calculate forward and adjoint sensitivities
  vector<vector<Mat> > fsens(fseed.size(),dae_out);
  vector<vector<Mat> > asens(aseed.size(),dae_in);
  f.eval(dae_in,dae_out,fseed,fsens,aseed,asens,true);
  
  // Augment forward state vectors
  x.append(vertcat(fwd_x));
  z.append(vertcat(fwd_z));
  xdot.append(vertcat(fwd_xdot));
  
  // Augment parameter vector
  p.append(vertcat(fwd_p));
  p.append(vertcat(adj_quad));
  
  // Augment backward state vectors
  rx.append(vertcat(adj_ode));
  rz.append(vertcat(adj_alg));
  rxdot.append(vertcat(adj_odedot));
  
  // Augment forward problem
  for(int dir=0; dir<nfwd; ++dir){
    ode.append(fsens[dir][DAE_ODE]);
    alg.append(fsens[dir][DAE_ALG]);
    quad.append(fsens[dir][DAE_QUAD]);
  }
  
  // Augment backwards problem
  for(int dir=0; dir<nadj; ++dir){
    Mat rode_sens = asens[dir][DAE_X];
    if(ode_is_explict){
      rode_sens += adj_odedot[dir];
    } else {
      rode_sens += asens[nadj+dir][DAE_XDOT];
    }
    rode.append(rode_sens);
    ralg.append(asens[dir][DAE_Z]);
    rquad.append(asens[dir][DAE_P]);
  }
  
  // Make sure that the augmented problem is dense
  makeDense(ode);
  makeDense(alg);
  makeDense(quad);
  makeDense(rode);
  makeDense(ralg);
  makeDense(rquad);
  
  // Update the forward problem inputs ...
  dae_in[DAE_X] = x;
  dae_in[DAE_Z] = z;
  dae_in[DAE_P] = p;
  dae_in[DAE_T] = t;
  dae_in[DAE_XDOT] = xdot;

  // ... and outputs
  dae_out[DAE_ODE] = ode;
  dae_out[DAE_ALG] = alg;
  dae_out[DAE_QUAD] = quad;
  
  // Update the backward problem inputs ...
  rdae_in[RDAE_RX] = rx;
  rdae_in[RDAE_RZ] = rz;
  rdae_in[RDAE_X] = x;
  rdae_in[RDAE_Z] = z;
  rdae_in[RDAE_P] = p;
  rdae_in[RDAE_T] = t;
  rdae_in[RDAE_XDOT] = xdot;
  rdae_in[RDAE_RXDOT] = rxdot;
  
  // ... and outputs
  rdae_out[RDAE_ODE] = rode;
  rdae_out[RDAE_ALG] = ralg;
  rdae_out[RDAE_QUAD] = rquad;
  
  // Create functions for the augmented problems
  XFunc f_aug(dae_in,dae_out);
  XFunc g_aug(rdae_in,rdae_out);

  // Workaround, delete g_aug if its empty
  if(g.isNull() && nadj==0) g_aug = XFunc();
  
  return pair<FX,FX>(f_aug,g_aug);
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



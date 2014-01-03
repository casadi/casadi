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
#include "../matrix/matrix_tools.hpp"
#include "../mx/mx_tools.hpp"
#include "../sx/sx_tools.hpp"
#include "mx_function.hpp"
#include "sx_function.hpp"

INPUTSCHEME(IntegratorInput)
OUTPUTSCHEME(IntegratorOutput)

using namespace std;
namespace CasADi{

  IntegratorInternal::IntegratorInternal(const FX& f, const FX& g) : GenericIntegratorInternal(f,g), f_(f), g_(g){
    // set default options
    setOption("name","unnamed_integrator"); // name of the function 
  
    // Additional options
    addOption("print_stats",              OT_BOOLEAN,     false, "Print out statistics after integration");
    addOption("t0",                       OT_REAL,        0.0, "Beginning of the time horizon"); 
    addOption("tf",                       OT_REAL,        1.0, "End of the time horizon");
    addOption("augmented_options",        OT_DICTIONARY,  GenericType(), "Options to be passed down to the augmented integrator, if one is constructed.");
  
    // Negative number of parameters for consistancy checking
    np_ = -1;
  
    input_.scheme = SCHEME_IntegratorInput;
    output_.scheme = SCHEME_IntegratorOutput;
  }

  IntegratorInternal::~IntegratorInternal(){ 
  }

  void IntegratorInternal::evaluate(){
    // Reset solver
    reset();

    // Integrate forward to the end of the time horizon
    integrate(tf_);

    // If backwards integration is needed
    if(nrx_>0){
      
      // Re-initialize backward problem
      resetB();
      
      // Integrate backwards to the beginning
      integrateB(t0_);
    }
      
    // Print statistics
    if(getOption("print_stats")) printStats(std::cout);
  }

  void IntegratorInternal::init(){
  
    // Initialize the functions
    casadi_assert(!f_.isNull());
  
    // Initialize and get dimensions for the forward integration
    if(!f_.isInit()) f_.init();
    casadi_assert_message(f_.getNumInputs()==DAE_NUM_IN,"Wrong number of inputs for the DAE callback function");
    casadi_assert_message(f_.getNumOutputs()==DAE_NUM_OUT,"Wrong number of outputs for the DAE callback function");
    nx_ = f_.input(DAE_X).size();
    nz_ = f_.input(DAE_Z).size();
    nq_ = f_.output(DAE_QUAD).size();
    np_  = f_.input(DAE_P).size();

    // Initialize and get dimensions for the backward integration
    if(g_.isNull()){
      // No backwards integration
      nrx_ = nrz_ = nrq_ = nrp_ = 0;
    } else {
      if(!g_.isInit()) g_.init();
      casadi_assert_message(g_.getNumInputs()==RDAE_NUM_IN,"Wrong number of inputs for the backwards DAE callback function");
      casadi_assert_message(g_.getNumOutputs()==RDAE_NUM_OUT,"Wrong number of outputs for the backwards DAE callback function");
      nrx_ = g_.input(RDAE_RX).size();
      nrz_ = g_.input(RDAE_RZ).size();
      nrp_ = g_.input(RDAE_RP).size();
      nrq_ = g_.output(RDAE_QUAD).size();
    }

    // Allocate space for inputs
    setNumInputs(INTEGRATOR_NUM_IN);
    input(INTEGRATOR_X0)  = DMatrix::zeros(f_.input(DAE_X).sparsity());
    input(INTEGRATOR_P)   = DMatrix::zeros(f_.input(DAE_P).sparsity());
    if(!g_.isNull()){
      input(INTEGRATOR_RX0)  = DMatrix::zeros(g_.input(RDAE_RX).sparsity());
      input(INTEGRATOR_RP)  = DMatrix::zeros(g_.input(RDAE_RP).sparsity());
    }
  
    // Allocate space for outputs
    setNumOutputs(INTEGRATOR_NUM_OUT);
    output(INTEGRATOR_XF) = input(INTEGRATOR_X0);
    output(INTEGRATOR_QF) = DMatrix::zeros(f_.output(DAE_QUAD).sparsity());
    if(!g_.isNull()){
      output(INTEGRATOR_RXF)  = input(INTEGRATOR_RX0);
      output(INTEGRATOR_RQF)  = DMatrix::zeros(g_.output(RDAE_QUAD).sparsity());
    }

    // Allocate space for algebraic variable
    z_ = DMatrix::zeros(f_.input(DAE_Z).sparsity());
    if(!g_.isNull()){
      rz_ = DMatrix::zeros(g_.input(RDAE_RZ).sparsity());
    }

    // Warn if sparse inputs (was previously an error)
    casadi_assert_warning(f_.input(DAE_X).dense(),"Sparse states in integrators are experimental");

    // Consistency checks
    casadi_assert_message(f_.output(DAE_ODE).shape()==input(INTEGRATOR_X0).shape(),"Inconsistent dimensions. Expecting DAE_ODE output of shape " << input(INTEGRATOR_X0).shape() << ", but got " << f_.output(DAE_ODE).shape() << " instead.");
    casadi_assert(f_.output(DAE_ODE).sparsity()==input(INTEGRATOR_X0).sparsity());
    casadi_assert_message(f_.output(DAE_ALG).shape()==z_.shape(),"Inconsistent dimensions. Expecting DAE_ALG output of shape " << z_.shape() << ", but got " << f_.output(DAE_ALG).shape() << " instead.");
    casadi_assert(f_.output(DAE_ALG).sparsity()==z_.sparsity());
    if(!g_.isNull()){
      casadi_assert(g_.input(RDAE_P).sparsity()==input(INTEGRATOR_P).sparsity());
      casadi_assert(g_.input(RDAE_X).sparsity()==input(INTEGRATOR_X0).sparsity());
      casadi_assert(g_.input(RDAE_Z).sparsity()==z_.sparsity());
      casadi_assert(g_.output(RDAE_ODE).sparsity()==input(INTEGRATOR_RX0).sparsity());
      casadi_assert(g_.output(RDAE_ALG).sparsity()==rz_.sparsity());
    }
  
    // Call the base class method
    FXInternal::init();

    {
      std::stringstream ss;
      ss << "Integrator dimensions: nx=" << nx_ << ", nz="<< nz_ << ", nq=" << nq_ << ", np=" << np_;
      log("IntegratorInternal::init",ss.str());
    }
  
    // read options
    t0_ = getOption("t0");
    tf_ = getOption("tf");
  }

  void IntegratorInternal::deepCopyMembers(std::map<SharedObjectNode*,SharedObject>& already_copied){
    FXInternal::deepCopyMembers(already_copied);
    f_ = deepcopy(f_,already_copied);
    g_ = deepcopy(g_,already_copied);
  }

  std::pair<FX,FX> IntegratorInternal::getAugmented(int nfwd, int nadj, vector<int>& xf_offset, vector<int>& qf_offset, vector<int>& rxf_offset, vector<int>& rqf_offset){
    log("IntegratorInternal::getAugmented","call");
    if(is_a<SXFunction>(f_)){
      casadi_assert_message(g_.isNull() || is_a<SXFunction>(g_), "Currently, g_ must be of the same type as f_");
      return getAugmentedGen<SXMatrix,SXFunction>(nfwd,nadj,xf_offset,qf_offset,rxf_offset,rqf_offset);
    } else if(is_a<MXFunction>(f_)){
      casadi_assert_message(g_.isNull() || is_a<MXFunction>(g_), "Currently, g_ must be of the same type as f_");
      return getAugmentedGen<MX,MXFunction>(nfwd,nadj,xf_offset,qf_offset,rxf_offset,rqf_offset);
    } else {
      throw CasadiException("Currently, f_ must be either SXFunction or MXFunction");
    }
  }
  
  template<class Mat,class XFunc>
  std::pair<FX,FX> IntegratorInternal::getAugmentedGen(int nfwd, int nadj, vector<int>& xf_offset, vector<int>& qf_offset, vector<int>& rxf_offset, vector<int>& rqf_offset){
  
    log("IntegratorInternal::getAugmentedGen","begin");

    // Reset the offset
    xf_offset.clear();    xf_offset.push_back(0);
    qf_offset.clear();    qf_offset.push_back(0);
    rxf_offset.clear();   rxf_offset.push_back(0);
    rqf_offset.clear();   rqf_offset.push_back(0);
  
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
    Mat ode = dae_out[DAE_ODE];
    Mat alg = dae_out[DAE_ALG];
    Mat quad = dae_out[DAE_QUAD];

    // Take apart the backwards problem
    vector<Mat> rdae_in(RDAE_NUM_IN), rdae_out(RDAE_NUM_OUT);
    if(!g.isNull()){
      rdae_in = g.inputExpr();
      rdae_out = g.outputExpr();
      // TODO: Assert that rdae_in[RDAE_X]==x, rdae_in[RDAE_Z]==z, rdae_in[RDAE_P]==p
    } else {
      rdae_in[RDAE_X]=x;
      rdae_in[RDAE_Z]=z;
      rdae_in[RDAE_P]=p;
      rdae_in[RDAE_T]=t;
    }
    Mat rx = rdae_in[RDAE_RX];
    Mat rz = rdae_in[RDAE_RZ];
    Mat rp = rdae_in[RDAE_RP];
    Mat rode = rdae_out[RDAE_ODE];
    Mat ralg = rdae_out[RDAE_ALG];
    Mat rquad = rdae_out[RDAE_QUAD];

    // Get offset for nondifferentiated problem
    if( nx_>0) xf_offset.push_back(x.size1());
    if( nq_>0) qf_offset.push_back(quad.size1());
    if(nrx_>0) rxf_offset.push_back(rx.size1());
    if(nrq_>0) rqf_offset.push_back(rquad.size1());

    // Function evaluating f and g
    vector<Mat> fg_out(DAE_NUM_OUT+RDAE_NUM_OUT);
    copy(dae_out.begin(),dae_out.end(),fg_out.begin());
    copy(rdae_out.begin(),rdae_out.end(),fg_out.begin()+DAE_NUM_OUT);
    XFunc fg(rdae_in,fg_out);
    fg.init();
  
    // Allocate forward sensitivities
    vector<Mat> fwd_x = Mat::sym("fwd_x",x.sparsity(),nfwd);
    vector<Mat> fwd_z = Mat::sym("fwd_z",z.sparsity(),nfwd);
    vector<Mat> fwd_p = Mat::sym("fwd_p",p.sparsity(),nfwd);
    vector<Mat> fwd_rx = Mat::sym("fwd_rx",rx.sparsity(),nfwd);
    vector<Mat> fwd_rz = Mat::sym("fwd_rz",rz.sparsity(),nfwd);
    vector<Mat> fwd_rp = Mat::sym("fwd_rp",rp.sparsity(),nfwd);

    // Allocate adjoint sensitivities
    vector<Mat> adj_ode = Mat::sym("adj_ode",ode.sparsity(),nadj);
    vector<Mat> adj_alg = Mat::sym("adj_alg",alg.sparsity(),nadj);
    vector<Mat> adj_quad = Mat::sym("adj_quad",quad.sparsity(),nadj);
    vector<Mat> adj_rode = Mat::sym("adj_rode",rode.sparsity(),nadj);
    vector<Mat> adj_ralg = Mat::sym("adj_ralg",ralg.sparsity(),nadj);
    vector<Mat> adj_rquad = Mat::sym("adj_rquad",rquad.sparsity(),nadj);
    
    // Forward seeds
    vector<vector<Mat> > fseed(nfwd,vector<Mat>(RDAE_NUM_IN));
    for(int dir=0; dir<nfwd; ++dir){
      fseed[dir][RDAE_X] = fwd_x[dir];
      fseed[dir][RDAE_Z] = fwd_z[dir];
      fseed[dir][RDAE_P] = fwd_p[dir];
      if(!t.isNull()) fseed[dir][RDAE_T] = Mat(t.sparsity());
      fseed[dir][RDAE_RX] = fwd_rx[dir];
      fseed[dir][RDAE_RZ] = fwd_rz[dir];
      fseed[dir][RDAE_RP] = fwd_rp[dir];

      // Save number of rows
      if( nx_>0) xf_offset.push_back(x.size1());
      if( nq_>0) qf_offset.push_back(quad.size1());
      if(nrx_>0) rxf_offset.push_back(rx.size1());
      if(nrq_>0) rqf_offset.push_back(rquad.size1());
    }

    // Adjoint seeds
    vector<vector<Mat> > aseed(nadj,vector<Mat>(DAE_NUM_OUT+RDAE_NUM_OUT));
    for(int dir=0; dir<nadj; ++dir){
      aseed[dir][DAE_ODE] = adj_ode[dir];
      aseed[dir][DAE_ALG] = adj_alg[dir];
      aseed[dir][DAE_QUAD] = adj_quad[dir];
    
      aseed[dir][DAE_NUM_OUT+RDAE_ODE] = adj_rode[dir];
      aseed[dir][DAE_NUM_OUT+RDAE_ALG] = adj_ralg[dir];
      aseed[dir][DAE_NUM_OUT+RDAE_QUAD] = adj_rquad[dir];

      // Save number of rows
      if( nx_>0) rxf_offset.push_back(x.size1());
      if( np_>0) rqf_offset.push_back(p.size1());
      if(nrx_>0) xf_offset.push_back(rx.size1());
      if(nrp_>0) qf_offset.push_back(rp.size1());
    }

    // Get cummulative offsets
    for(int i=1; i<xf_offset.size(); ++i) xf_offset[i] += xf_offset[i-1];
    for(int i=1; i<qf_offset.size(); ++i) qf_offset[i] += qf_offset[i-1];
    for(int i=1; i<rxf_offset.size(); ++i) rxf_offset[i] += rxf_offset[i-1];
    for(int i=1; i<rqf_offset.size(); ++i) rqf_offset[i] += rqf_offset[i-1];
  
    // Calculate forward and adjoint sensitivities
    vector<vector<Mat> > fsens(fseed.size(),fg_out);
    vector<vector<Mat> > asens(aseed.size(),rdae_in);
    fg.eval(rdae_in,fg_out,fseed,fsens,aseed,asens);
  
    // Augment differential state
    x.append(vertcat(fwd_x));
    x.append(vertcat(adj_rode));
  
    // Augment algebraic state
    z.append(vertcat(fwd_z));
    z.append(vertcat(adj_ralg));
  
    // Augment parameter vector
    p.append(vertcat(fwd_p));
    p.append(vertcat(adj_rquad));
  
    // Augment backward differential state
    rx.append(vertcat(fwd_rx));
    rx.append(vertcat(adj_ode));
  
    // Augment backward algebraic state
    rz.append(vertcat(fwd_rz));
    rz.append(vertcat(adj_alg));
  
    // Augment backwards parameter vector
    rp.append(vertcat(fwd_rp));
    rp.append(vertcat(adj_quad));
  
    // Augment forward sensitivity equations to the DAE
    for(int dir=0; dir<nfwd; ++dir){
      ode.append(fsens[dir][DAE_ODE]);
      alg.append(fsens[dir][DAE_ALG]);
      quad.append(fsens[dir][DAE_QUAD]);

      rode.append(fsens[dir][DAE_NUM_OUT+RDAE_ODE]);
      ralg.append(fsens[dir][DAE_NUM_OUT+RDAE_ALG]);
      rquad.append(fsens[dir][DAE_NUM_OUT+RDAE_QUAD]);
    }
  
    // Augment backward sensitivity equations to the DAE
    for(int dir=0; dir<nadj; ++dir){
      rode.append(asens[dir][RDAE_X]);
      ralg.append(asens[dir][RDAE_Z]);
      rquad.append(asens[dir][RDAE_P]);
      ode.append(asens[dir][RDAE_RX]);
      alg.append(asens[dir][RDAE_RZ]);
      quad.append(asens[dir][RDAE_RP]);
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

    // ... and outputs
    dae_out[DAE_ODE] = ode;
    dae_out[DAE_ALG] = alg;
    dae_out[DAE_QUAD] = quad;
  
    // Update the backward problem inputs ...
    rdae_in[RDAE_RX] = rx;
    rdae_in[RDAE_RZ] = rz;
    rdae_in[RDAE_RP] = rp;
    rdae_in[RDAE_X] = x;
    rdae_in[RDAE_Z] = z;
    rdae_in[RDAE_P] = p;
    rdae_in[RDAE_T] = t;
  
    // ... and outputs
    rdae_out[RDAE_ODE] = rode;
    rdae_out[RDAE_ALG] = ralg;
    rdae_out[RDAE_QUAD] = rquad;
  
    // Create functions for the augmented problems
    XFunc f_aug(dae_in,dae_out);
    XFunc g_aug(rdae_in,rdae_out);

    f_aug.init();
  
    casadi_assert_message(f_aug.getFree().size()==0,"IntegratorInternal::getDerivative: Found free variables " << f_aug.getFree() << " while constructing augmented dae. Make sure that gx, gz and gq have a linear dependency on rx, rz and rp. This is a restriction of the implementation.");
  
    // Workaround, delete g_aug if its empty
    if(g.isNull() && nadj==0) g_aug = XFunc();
  
    log("IntegratorInternal::getAugmentedGen","end");
    
    return pair<FX,FX>(f_aug,g_aug);
  }

  void IntegratorInternal::spEvaluate(bool fwd){
    log("IntegratorInternal::spEvaluate","begin");
    /**  This is a bit better than the FXInternal implementation: XF and QF never depend on RX0 and RP, 
     *   i.e. the worst-case structure of the Jacobian is:
     *        x0  p rx0 rp
     *        --------------
     *   xf  | x  x        |
     *   qf  | x  x        |
     *  rxf  | x  x  x  x  |
     *  rqf  | x  x  x  x  |
     *        --------------
     * 
     *  An even better structure of the Jacobian can be obtained by propagating sparsity through the callback functions.
     */
  
    // Variable which depends on all states and parameters
    bvec_t all_depend(0);
  
    if(fwd){
    
      // Have dependency on anything in x0 or p
      for(int k=0; k<2; ++k){
        int iind = k==0 ? INTEGRATOR_X0 : INTEGRATOR_P;
        const DMatrix& m = inputNoCheck(iind);
        const bvec_t* v = reinterpret_cast<const bvec_t*>(m.ptr());
        for(int i=0; i<m.size(); ++i){
          all_depend |= v[i];
        }
      }
    
      // Propagate to xf and qf (that only depend on x0 and p)
      for(int k=0; k<2; ++k){
        int oind = k==0 ? INTEGRATOR_XF : INTEGRATOR_QF;
        DMatrix& m = outputNoCheck(oind);
        bvec_t* v = reinterpret_cast<bvec_t*>(m.ptr());
        for(int i=0; i<m.size(); ++i){
          v[i] = all_depend;
        }
      }
    
      // Add dependency on rx0 or rp
      for(int k=0; k<2; ++k){
        int iind = k==0 ? INTEGRATOR_RX0 : INTEGRATOR_RP;
        const DMatrix& m = inputNoCheck(iind);
        const bvec_t* v = reinterpret_cast<const bvec_t*>(m.ptr());
        for(int i=0; i<m.size(); ++i){
          all_depend |= v[i];
        }
      }
    
      // Propagate to rxf and rqf
      for(int k=0; k<2; ++k){
        int oind = k==0 ? INTEGRATOR_RXF : INTEGRATOR_RQF;
        DMatrix& m = outputNoCheck(oind);
        bvec_t* v = reinterpret_cast<bvec_t*>(m.ptr());
        for(int i=0; i<m.size(); ++i){
          v[i] = all_depend;
        }
      }
    
    } else {
    
      // First find out what influences only rxf and rqf
      for(int k=0; k<2; ++k){
        int oind = k==0 ? INTEGRATOR_RXF : INTEGRATOR_RQF;
        const DMatrix& m = outputNoCheck(oind);
        const bvec_t* v = get_bvec_t(m.data());
        for(int i=0; i<m.size(); ++i){
          all_depend |= v[i];
        }
      }
    
      // Propagate to rx0 and rp
      for(int k=0; k<2; ++k){
        int iind = k==0 ? INTEGRATOR_RX0 : INTEGRATOR_RP;
        DMatrix& m = inputNoCheck(iind);
        bvec_t* v = get_bvec_t(m.data());
        for(int i=0; i<m.size(); ++i){
          v[i] = all_depend;
        }
      }
    
      // Add dependencies to xf and qf
      for(int k=0; k<2; ++k){
        int oind = k==0 ? INTEGRATOR_XF : INTEGRATOR_QF;
        const DMatrix& m = outputNoCheck(oind);
        const bvec_t* v = get_bvec_t(m.data());
        for(int i=0; i<m.size(); ++i){
          all_depend |= v[i];
        }
      }
    
      // Propagate to x0 and p
      for(int k=0; k<2; ++k){
        int iind = k==0 ? INTEGRATOR_X0 : INTEGRATOR_P;
        DMatrix& m = inputNoCheck(iind);
        bvec_t* v = get_bvec_t(m.data());
        for(int i=0; i<m.size(); ++i){
          v[i] = all_depend;
        }
      }
    }
    log("IntegratorInternal::spEvaluate","end");
  }

  FX IntegratorInternal::getDerivative(int nfwd, int nadj){
    log("IntegratorInternal::getDerivative","begin");

    // Generate augmented DAE
    vector<int> xf_offset, qf_offset, rxf_offset, rqf_offset;
    std::pair<FX,FX> aug_dae = getAugmented(nfwd,nadj,xf_offset,qf_offset,rxf_offset,rqf_offset);
  
    // Create integrator for augmented DAE
    Integrator integrator;
    integrator.assignNode(create(aug_dae.first,aug_dae.second));
  
    // Copy options
    integrator.setOption(dictionary());
  
    // Pass down specific options if provided
    if (hasSetOption("augmented_options"))
      integrator.setOption(getOption("augmented_options"));
  
    // Initialize the integrator since we will call it below
    integrator.init();
  
    // All inputs of the return function
    vector<MX> ret_in;
    ret_in.reserve(INTEGRATOR_NUM_IN*(1+nfwd) + INTEGRATOR_NUM_OUT*nadj);
  
    // Augmented state
    MX x0_aug, p_aug, rx0_aug, rp_aug;
  
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

      // Backward parameter
      ss.clear();
      ss << "rp";
      if(dir>=0) ss << "_" << dir;
      dd[INTEGRATOR_RP] = msym(ss.str(),input(INTEGRATOR_RP).sparsity());
      rp_aug.append(dd[INTEGRATOR_RP]);
    
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

      // Quadratures become backward parameters
      ss.clear();
      ss << "qf" << "_" << dir;
      dd[INTEGRATOR_QF] = msym(ss.str(),output(INTEGRATOR_QF).sparsity());
      rp_aug.append(dd[INTEGRATOR_QF]);

      // Backward differential states becomes forward differential states
      ss.clear();
      ss << "rxf" << "_" << dir;
      dd[INTEGRATOR_RXF] = msym(ss.str(),output(INTEGRATOR_RXF).sparsity());
      x0_aug.append(dd[INTEGRATOR_RXF]);
    
      // Backward quadratures becomes (forward) parameters
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
    integrator_in[INTEGRATOR_RP] = rp_aug;
    vector<MX> integrator_out = integrator.call(integrator_in);
  
    // Augmented results
    vector<MX> xf_aug = vertsplit(integrator_out[INTEGRATOR_XF],xf_offset);
    vector<MX> rxf_aug = vertsplit(integrator_out[INTEGRATOR_RXF],rxf_offset);
    vector<MX> qf_aug = vertsplit(integrator_out[INTEGRATOR_QF],qf_offset);
    vector<MX> rqf_aug = vertsplit(integrator_out[INTEGRATOR_RQF],rqf_offset);
    vector<MX>::const_iterator xf_aug_it = xf_aug.begin();
    vector<MX>::const_iterator rxf_aug_it = rxf_aug.begin();
    vector<MX>::const_iterator qf_aug_it = qf_aug.begin();
    vector<MX>::const_iterator rqf_aug_it = rqf_aug.begin();

    // All outputs of the return function
    vector<MX> ret_out;
    ret_out.reserve(INTEGRATOR_NUM_OUT*(1+nfwd) + INTEGRATOR_NUM_IN*nadj);

    // Collect the nondifferentiated results and forward sensitivities
    dd.resize(INTEGRATOR_NUM_OUT);
    fill(dd.begin(),dd.end(),MX());
    for(int dir=-1; dir<nfwd; ++dir){
      if( nx_>0) dd[INTEGRATOR_XF]  = *xf_aug_it++;
      if( nq_>0) dd[INTEGRATOR_QF]  = *qf_aug_it++;
      if(nrx_>0) dd[INTEGRATOR_RXF] = *rxf_aug_it++;
      if(nrq_>0) dd[INTEGRATOR_RQF] = *rqf_aug_it++;
      ret_out.insert(ret_out.end(),dd.begin(),dd.end());
    }
  
    // Collect the adjoint sensitivities
    dd.resize(INTEGRATOR_NUM_IN);
    fill(dd.begin(),dd.end(),MX());
    for(int dir=0; dir<nadj; ++dir){
      if( nx_>0) dd[INTEGRATOR_X0]  = *rxf_aug_it++;
      if( np_>0) dd[INTEGRATOR_P]   = *rqf_aug_it++;
      if(nrx_>0) dd[INTEGRATOR_RX0] = *xf_aug_it++;
      if(nrp_>0) dd[INTEGRATOR_RP]  = *qf_aug_it++;
      ret_out.insert(ret_out.end(),dd.begin(),dd.end());
    }
    log("IntegratorInternal::getDerivative","end");
  
    // Create derivative function and return
    return MXFunction(ret_in,ret_out);
  }

  FX IntegratorInternal::getJacobian(int iind, int oind, bool compact, bool symmetric){
    vector<MX> arg = symbolicInput();
    vector<MX> res = shared_from_this<FX>().call(arg);
    MXFunction f(arg,res);
    f.setOption("ad_mode","forward");
    f.init();
    return f.jacobian(iind,oind,compact,symmetric);
  }

  void IntegratorInternal::reset(){
    log("IntegratorInternal::reset","begin");
    
    // Initialize output (relevant for integration with a zero advance time )
    copy(input(INTEGRATOR_X0).begin(),input(INTEGRATOR_X0).end(),output(INTEGRATOR_XF).begin());
    
    log("IntegratorInternal::reset","end");
  }


} // namespace CasADi



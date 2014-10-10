/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A symbolic framework for dynamic optimization.
 *    Copyright (C) 2010-2014 Joel Andersson, Joris Gillis, Moritz Diehl,
 *                            K.U. Leuven. All rights reserved.
 *    Copyright (C) 2011-2014 Greg Horn
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
#include "../std_vector_tools.hpp"
#include "../matrix/matrix_tools.hpp"
#include "../mx/mx_tools.hpp"
#include "../sx/sx_tools.hpp"
#include "mx_function.hpp"
#include "sx_function.hpp"

INPUTSCHEME(IntegratorInput)
OUTPUTSCHEME(IntegratorOutput)

using namespace std;
namespace casadi {

  IntegratorInternal::IntegratorInternal(const Function& f, const Function& g) : f_(f), g_(g) {
    // set default options
    setOption("name", "unnamed_integrator"); // name of the function

    // Additional options
    addOption("print_stats",              OT_BOOLEAN,     false,
              "Print out statistics after integration");
    addOption("t0",                       OT_REAL,        0.0,
              "Beginning of the time horizon");
    addOption("tf",                       OT_REAL,        1.0,
              "End of the time horizon");
    addOption("augmented_options",        OT_DICTIONARY,  GenericType(),
              "Options to be passed down to the augmented integrator, if one is constructed.");
    addOption("expand_augmented",         OT_BOOLEAN,     true,
              "If DAE callback functions are SXFunction, have augmented"
              " DAE callback function also be SXFunction.");

    // Negative number of parameters for consistancy checking
    np_ = -1;

    input_.scheme = SCHEME_IntegratorInput;
    output_.scheme = SCHEME_IntegratorOutput;
  }

  IntegratorInternal::~IntegratorInternal() {
  }

  void IntegratorInternal::evaluate() {
    // Reset solver
    reset();

    // Integrate forward to the end of the time horizon
    integrate(tf_);

    // If backwards integration is needed
    if (nrx_>0) {

      // Re-initialize backward problem
      resetB();

      // Integrate backwards to the beginning
      integrateB(t0_);
    }

    // Print statistics
    if (print_stats_) printStats(std::cout);
  }

  void IntegratorInternal::init() {

    // Initialize the functions
    casadi_assert(!f_.isNull());

    // Initialize and get dimensions for the forward integration
    if (!f_.isInit()) f_.init();
    casadi_assert_message(f_.getNumInputs()==DAE_NUM_IN,
                          "Wrong number of inputs for the DAE callback function");
    casadi_assert_message(f_.getNumOutputs()==DAE_NUM_OUT,
                          "Wrong number of outputs for the DAE callback function");
    nx_ = f_.input(DAE_X).size();
    nz_ = f_.input(DAE_Z).size();
    nq_ = f_.output(DAE_QUAD).size();
    np_  = f_.input(DAE_P).size();

    // Initialize and get dimensions for the backward integration
    if (g_.isNull()) {
      // No backwards integration
      nrx_ = nrz_ = nrq_ = nrp_ = 0;
    } else {
      if (!g_.isInit()) g_.init();
      casadi_assert_message(g_.getNumInputs()==RDAE_NUM_IN,
                            "Wrong number of inputs for the backwards DAE callback function");
      casadi_assert_message(g_.getNumOutputs()==RDAE_NUM_OUT,
                            "Wrong number of outputs for the backwards DAE callback function");
      nrx_ = g_.input(RDAE_RX).size();
      nrz_ = g_.input(RDAE_RZ).size();
      nrp_ = g_.input(RDAE_RP).size();
      nrq_ = g_.output(RDAE_QUAD).size();
    }

    // Allocate space for inputs
    setNumInputs(INTEGRATOR_NUM_IN);
    x0()  = DMatrix::zeros(f_.input(DAE_X).sparsity());
    p()   = DMatrix::zeros(f_.input(DAE_P).sparsity());
    z0()   = DMatrix::zeros(f_.input(DAE_Z).sparsity());
    if (!g_.isNull()) {
      rx0()  = DMatrix::zeros(g_.input(RDAE_RX).sparsity());
      rp()  = DMatrix::zeros(g_.input(RDAE_RP).sparsity());
      rz0()  = DMatrix::zeros(g_.input(RDAE_RZ).sparsity());
    }

    // Allocate space for outputs
    setNumOutputs(INTEGRATOR_NUM_OUT);
    xf() = x0();
    qf() = DMatrix::zeros(f_.output(DAE_QUAD).sparsity());
    zf() = z0();
    if (!g_.isNull()) {
      rxf()  = rx0();
      rqf()  = DMatrix::zeros(g_.output(RDAE_QUAD).sparsity());
      rzf()  = rz0();
    }

    // Warn if sparse inputs (was previously an error)
    casadi_assert_warning(f_.input(DAE_X).isDense(),
                          "Sparse states in integrators are experimental");

    // Consistency checks
    casadi_assert_message(f_.output(DAE_ODE).shape()==x0().shape(),
                          "Inconsistent dimensions. Expecting DAE_ODE output of shape "
                          << x0().shape() << ", but got "
                          << f_.output(DAE_ODE).shape() << " instead.");
    casadi_assert(f_.output(DAE_ODE).sparsity()==x0().sparsity());
    casadi_assert_message(f_.output(DAE_ALG).shape()==z0().shape(),
                          "Inconsistent dimensions. Expecting DAE_ALG output of shape "
                          << z0().shape() << ", but got "
                          << f_.output(DAE_ALG).shape() << " instead.");
    casadi_assert(f_.output(DAE_ALG).sparsity()==z0().sparsity());
    if (!g_.isNull()) {
      casadi_assert(g_.input(RDAE_P).sparsity()==p().sparsity());
      casadi_assert(g_.input(RDAE_X).sparsity()==x0().sparsity());
      casadi_assert(g_.input(RDAE_Z).sparsity()==z0().sparsity());
      casadi_assert(g_.output(RDAE_ODE).sparsity()==rx0().sparsity());
      casadi_assert(g_.output(RDAE_ALG).sparsity()==rz0().sparsity());
    }

    // Call the base class method
    FunctionInternal::init();

    {
      std::stringstream ss;
      ss << "Integrator dimensions: nx=" << nx_ << ", nz="<< nz_
         << ", nq=" << nq_ << ", np=" << np_;
      log("IntegratorInternal::init", ss.str());
    }

    // read options
    t0_ = getOption("t0");
    tf_ = getOption("tf");
    print_stats_ = getOption("print_stats");

    // Form a linear solver for the sparsity propagation
    linsol_f_ = LinearSolver(spJacF());
    linsol_f_.init();
    if (!g_.isNull()) {
      linsol_g_ = LinearSolver(spJacG());
      linsol_g_.init();
    }
  }

  void IntegratorInternal::deepCopyMembers(
      std::map<SharedObjectNode*, SharedObject>& already_copied) {
    FunctionInternal::deepCopyMembers(already_copied);
    f_ = deepcopy(f_, already_copied);
    g_ = deepcopy(g_, already_copied);
    linsol_f_ = deepcopy(linsol_f_, already_copied);
    linsol_g_ = deepcopy(linsol_g_, already_copied);
  }

  std::pair<Function, Function> IntegratorInternal::getAugmented(int nfwd, int nadj,
                                                                AugOffset& offset) {
    log("IntegratorInternal::getAugmented", "call");

    //    cout << "here" << endl;

    // Return object
    std::pair<Function, Function> ret;

    // Calculate offsets
    offset = getAugOffset(nfwd, nadj);

    // Create augmented problem
    MX aug_t = MX::sym("aug_t", f_.input(DAE_T).sparsity());
    MX aug_x = MX::sym("aug_x", x0().size1(), offset.x.back());
    MX aug_z = MX::sym("aug_z", std::max(z0().size1(), rz0().size1()), offset.z.back());
    MX aug_p = MX::sym("aug_p", std::max(p().size1(), rp().size1()), offset.p.back());
    MX aug_rx = MX::sym("aug_rx", x0().size1(), offset.rx.back());
    MX aug_rz = MX::sym("aug_rz", std::max(z0().size1(), rz0().size1()), offset.rz.back());
    MX aug_rp = MX::sym("aug_rp", std::max(qf().size1(), rp().size1()), offset.rp.back());

    // Split up the augmented vectors
    vector<MX> aug_x_split = horzsplit(aug_x, offset.x);
    vector<MX>::const_iterator aug_x_split_it = aug_x_split.begin();
    vector<MX> aug_z_split = horzsplit(aug_z, offset.z);
    vector<MX>::const_iterator aug_z_split_it = aug_z_split.begin();
    vector<MX> aug_p_split = horzsplit(aug_p, offset.p);
    vector<MX>::const_iterator aug_p_split_it = aug_p_split.begin();
    vector<MX> aug_rx_split = horzsplit(aug_rx, offset.rx);
    vector<MX>::const_iterator aug_rx_split_it = aug_rx_split.begin();
    vector<MX> aug_rz_split = horzsplit(aug_rz, offset.rz);
    vector<MX>::const_iterator aug_rz_split_it = aug_rz_split.begin();
    vector<MX> aug_rp_split = horzsplit(aug_rp, offset.rp);
    vector<MX>::const_iterator aug_rp_split_it = aug_rp_split.begin();

    // Temporary vector
    vector<MX> tmp;

    // Zero with the dimension of t
    MX zero_t = DMatrix::zeros(aug_t.sparsity());

    // The DAE being constructed
    vector<MX> f_ode, f_alg, f_quad, g_ode, g_alg, g_quad;

    // Forward derivatives of f
    Function d = f_.derivative(nfwd, 0);
    vector<MX> f_arg;
    f_arg.reserve(d.getNumInputs());
    tmp.resize(DAE_NUM_IN);
    fill(tmp.begin(), tmp.end(), MX());

    // Collect arguments for calling d
    for (int dir=-1; dir<nfwd; ++dir) {
      tmp[DAE_T] = dir<0 ? aug_t : zero_t;
      if ( nx_>0) tmp[DAE_X] = *aug_x_split_it++;
      if ( nz_>0) tmp[DAE_Z] = *aug_z_split_it++;
      if ( np_>0) tmp[DAE_P] = *aug_p_split_it++;
      f_arg.insert(f_arg.end(), tmp.begin(), tmp.end());
    }

    // Call d
    vector<MX> res = d.call(f_arg);
    vector<MX>::const_iterator res_it = res.begin();

    // Collect right-hand-sides
    tmp.resize(DAE_NUM_OUT);
    fill(tmp.begin(), tmp.end(), MX());
    for (int dir=-1; dir<nfwd; ++dir) {
      copy(res_it, res_it+tmp.size(), tmp.begin());
      res_it += tmp.size();
      if ( nx_>0) f_ode.push_back(tmp[DAE_ODE]);
      if ( nz_>0) f_alg.push_back(tmp[DAE_ALG]);
      if ( nq_>0) f_quad.push_back(tmp[DAE_QUAD]);
    }

    // Consistency check
    casadi_assert(res_it==res.end());

    vector<MX> g_arg;
    if (!g_.isNull()) {

      // Forward derivatives of g
      d = g_.derivative(nfwd, 0);
      g_arg.reserve(d.getNumInputs());
      tmp.resize(RDAE_NUM_IN);
      fill(tmp.begin(), tmp.end(), MX());

      // Reset iterators
      aug_x_split_it = aug_x_split.begin();
      aug_z_split_it = aug_z_split.begin();
      aug_p_split_it = aug_p_split.begin();

      // Collect arguments for calling d
      for (int dir=-1; dir<nfwd; ++dir) {
        tmp[RDAE_T] = dir<0 ? aug_t : zero_t;
        if ( nx_>0) tmp[RDAE_X] = *aug_x_split_it++;
        if ( nz_>0) tmp[RDAE_Z] = *aug_z_split_it++;
        if ( np_>0) tmp[RDAE_P] = *aug_p_split_it++;
        if (nrx_>0) tmp[RDAE_RX] = *aug_rx_split_it++;
        if (nrz_>0) tmp[RDAE_RZ] = *aug_rz_split_it++;
        if (nrp_>0) tmp[RDAE_RP] = *aug_rp_split_it++;
        g_arg.insert(g_arg.end(), tmp.begin(), tmp.end());
      }

      // Call d
      res = d.call(g_arg);
      res_it = res.begin();

      // Collect right-hand-sides
      tmp.resize(RDAE_NUM_OUT);
      fill(tmp.begin(), tmp.end(), MX());
      for (int dir=-1; dir<nfwd; ++dir) {
        copy(res_it, res_it+tmp.size(), tmp.begin());
        res_it += tmp.size();
        if (nrx_>0) g_ode.push_back(tmp[RDAE_ODE]);
        if (nrz_>0) g_alg.push_back(tmp[RDAE_ALG]);
        if (nrq_>0) g_quad.push_back(tmp[RDAE_QUAD]);
      }

      // Consistency check
      casadi_assert(res_it==res.end());
    }

    if (nadj>0) {

      // Adjoint derivatives of f
      d = f_.derivative(0, nadj);
      f_arg.resize(DAE_NUM_IN);
      f_arg.reserve(d.getNumInputs());

      // Collect arguments for calling d
      tmp.resize(DAE_NUM_OUT);
      fill(tmp.begin(), tmp.end(), MX());
      for (int dir=0; dir<nadj; ++dir) {
        if ( nx_>0) tmp[DAE_ODE] = *aug_rx_split_it++;
        if ( nz_>0) tmp[DAE_ALG] = *aug_rz_split_it++;
        if ( nq_>0) tmp[DAE_QUAD] = *aug_rp_split_it++;
        f_arg.insert(f_arg.end(), tmp.begin(), tmp.end());
      }

      // Call der
      res = d.call(f_arg);
      res_it = res.begin() + DAE_NUM_OUT;

      // Record locations in augg for later
      int g_ode_ind = g_ode.size();
      int g_alg_ind = g_alg.size();
      int g_quad_ind = g_quad.size();

      // Collect right-hand-sides
      tmp.resize(DAE_NUM_IN);
      for (int dir=0; dir<nadj; ++dir) {
        copy(res_it, res_it+tmp.size(), tmp.begin());
        res_it += tmp.size();
        if ( nx_>0) g_ode.push_back(tmp[DAE_X]);
        if ( nz_>0) g_alg.push_back(tmp[DAE_Z]);
        if ( np_>0) g_quad.push_back(tmp[DAE_P]);
      }

      // Consistency check
      casadi_assert(res_it==res.end());

      if (!g_.isNull()) {

        // Adjoint derivatives of g
        d = g_.derivative(0, nadj);
        g_arg.resize(RDAE_NUM_IN);
        g_arg.reserve(d.getNumInputs());

        // Collect arguments for calling der
        tmp.resize(RDAE_NUM_OUT);
        fill(tmp.begin(), tmp.end(), MX());
        for (int dir=0; dir<nadj; ++dir) {
          if (nrx_>0) tmp[RDAE_ODE] = *aug_x_split_it++;
          if (nrz_>0) tmp[RDAE_ALG] = *aug_z_split_it++;
          if (nrq_>0) tmp[RDAE_QUAD] = *aug_p_split_it++;
          g_arg.insert(g_arg.end(), tmp.begin(), tmp.end());
        }

        // Call der
        res = d.call(g_arg);
        res_it = res.begin() + RDAE_NUM_OUT;

        // Collect right-hand-sides
        tmp.resize(RDAE_NUM_IN);
        for (int dir=0; dir<nadj; ++dir) {
          copy(res_it, res_it+tmp.size(), tmp.begin());
          res_it += tmp.size();
          if ( nx_>0) g_ode[g_ode_ind++] += tmp[RDAE_X];
          if ( nz_>0) g_alg[g_alg_ind++] += tmp[RDAE_Z];
          if ( np_>0) g_quad[g_quad_ind++] += tmp[RDAE_P];
        }

        // Consistency check
        casadi_assert(g_ode_ind == g_ode.size());
        casadi_assert(g_alg_ind == g_alg.size());
        casadi_assert(g_quad_ind == g_quad.size());

        // Remove the dependency of rx, rz, rp in the forward integration (see Joel's thesis)
        if (nrx_>0) g_arg[RDAE_RX] = MX::zeros(g_arg[RDAE_RX].sparsity());
        if (nrz_>0) g_arg[RDAE_RZ] = MX::zeros(g_arg[RDAE_RZ].sparsity());
        if (nrp_>0) g_arg[RDAE_RP] = MX::zeros(g_arg[RDAE_RP].sparsity());

        // Call der again
        res = d.call(g_arg);
        res_it = res.begin() + RDAE_NUM_OUT;

        // Collect right-hand-sides and add contribution to the forward integration
        tmp.resize(RDAE_NUM_IN);
        for (int dir=0; dir<nadj; ++dir) {
          copy(res_it, res_it+tmp.size(), tmp.begin());
          res_it += tmp.size();
          if (nrx_>0) f_ode.push_back(tmp[RDAE_RX]);
          if (nrz_>0) f_alg.push_back(tmp[RDAE_RZ]);
          if (nrp_>0) f_quad.push_back(tmp[RDAE_RP]);
        }

        // Consistency check
        casadi_assert(res_it==res.end());
      }
    }

    // Do we want to expand MXFunction->SXFunction?
    bool expand = getOption("expand_augmented");

    // Can we expand?
    expand = expand && is_a<SXFunction>(f_) && (g_.isNull() || is_a<SXFunction>(g_));

    // Form the augmented forward integration function
    if (g_.isNull() && nfwd==0) {
      ret.first = f_; // reuse the existing one
    } else {
      vector<MX> f_in(DAE_NUM_IN), f_out(DAE_NUM_OUT);
      f_in[DAE_T] = aug_t;
      f_in[DAE_X] = aug_x;
      f_in[DAE_Z] = aug_z;
      f_in[DAE_P] = aug_p;
      if (!f_ode.empty()) f_out[DAE_ODE] = dense(horzcat(f_ode));
      if (!f_alg.empty()) f_out[DAE_ALG] = dense(horzcat(f_alg));
      if (!f_quad.empty()) f_out[DAE_QUAD] = dense(horzcat(f_quad));
      MXFunction f_mx(f_in, f_out);

      // Expand to SXFuncion?
      if (expand) {
        f_mx.init();
        ret.first = SXFunction(f_mx);
      } else {
        ret.first = f_mx;
      }
    }

    // Form the augmented backward integration function
    if (!g_ode.empty()) {
      vector<MX> g_in(RDAE_NUM_IN), g_out(RDAE_NUM_OUT);
      g_in[RDAE_T] = aug_t;
      g_in[RDAE_X] = aug_x;
      g_in[RDAE_Z] = aug_z;
      g_in[RDAE_P] = aug_p;
      g_in[RDAE_RX] = aug_rx;
      g_in[RDAE_RZ] = aug_rz;
      g_in[RDAE_RP] = aug_rp;
      if (!g_ode.empty()) g_out[RDAE_ODE] = dense(horzcat(g_ode));
      if (!g_alg.empty()) g_out[RDAE_ALG] = dense(horzcat(g_alg));
      if (!g_quad.empty()) g_out[RDAE_QUAD] = dense(horzcat(g_quad));
      MXFunction g_mx(g_in, g_out);

      // Expand to SXFuncion?
      if (expand) {
        g_mx.init();
        ret.second = SXFunction(g_mx);
      } else {
        ret.second = g_mx;
      }
    }

    // Consistency check
    casadi_assert(aug_x_split_it == aug_x_split.end());
    casadi_assert(aug_z_split_it == aug_z_split.end());
    casadi_assert(aug_p_split_it == aug_p_split.end());
    casadi_assert(aug_rx_split_it == aug_rx_split.end());
    casadi_assert(aug_rz_split_it == aug_rz_split.end());
    casadi_assert(aug_rp_split_it == aug_rp_split.end());

    // Return functions
    return ret;
  }

  void IntegratorInternal::spEvaluate(bool fwd) {
    log("IntegratorInternal::spEvaluate", "begin");

    // Temporary vectors
    bvec_t *tmp_f1, *tmp_f2, *tmp_g1= NULL, *tmp_g2=NULL;
    casadi_assert(!linsol_f_.isNull());
    tmp_f1 = reinterpret_cast<bvec_t*>(linsol_f_.input(LINSOL_B).ptr());
    tmp_f2 = reinterpret_cast<bvec_t*>(linsol_f_.output(LINSOL_X).ptr());
    if (!g_.isNull()) {
      casadi_assert(!linsol_g_.isNull());
      tmp_g1 = reinterpret_cast<bvec_t*>(linsol_g_.input(LINSOL_B).ptr());
      tmp_g2 = reinterpret_cast<bvec_t*>(linsol_g_.output(LINSOL_X).ptr());
    }

    // Initialize the callback functions for sparsity propagation
    f_.spInit(fwd);
    if (!g_.isNull()) {
      g_.spInit(fwd);
    }

    if (fwd) {

      // Propagate through the DAE
      f_.input(DAE_T).setZeroBV();
      f_.input(DAE_X).setBV(x0());
      f_.input(DAE_P).setBV(p());
      f_.input(DAE_Z).setZeroBV();
      f_.spEvaluate(true);
      f_.output(DAE_ODE).getArrayBV(tmp_f1, nx_);
      f_.output(DAE_ALG).getArrayBV(tmp_f1+nx_, nz_);

      // Propagate interdependencies
      x0().getArrayBV(tmp_f2, nx_);
      std::fill(tmp_f2+nx_, tmp_f2+nx_+nz_, 0);
      linsol_f_.spSolve(tmp_f2, tmp_f1, true);
      xf().setArrayBV(tmp_f2, nx_);
      zf().setArrayBV(tmp_f2+nx_, nz_);

      // Get influence on the quadratures
      if (nq_>0) {
        f_.input(DAE_X).setArrayBV(tmp_f2, nx_);
        f_.input(DAE_Z).setBV(zf());
        f_.spEvaluate(true);
        f_.output(DAE_QUAD).getBV(qf());
      }

      // Propagate through g
      if (!g_.isNull()) {

        // Propagate through the backward DAE
        g_.input(RDAE_T).setZeroBV();
        g_.input(RDAE_X).setBV(xf());
        g_.input(RDAE_P).setBV(p());
        g_.input(RDAE_Z).setBV(zf());
        g_.input(RDAE_RX).setBV(rx0());
        g_.input(RDAE_RP).setBV(rp());
        g_.input(RDAE_RZ).setZeroBV();
        g_.spEvaluate(true);
        g_.output(RDAE_ODE).getArrayBV(tmp_g1, nrx_);
        g_.output(RDAE_ALG).getArrayBV(tmp_g1+nrx_, nrz_);

        // Propagate interdependencies
        rx0().getArrayBV(tmp_g2, nrx_);
        std::fill(tmp_g2+nrx_, tmp_g2+nrx_+nrz_, 0);
        linsol_g_.spSolve(tmp_g2, tmp_g1, true);
        rxf().setArrayBV(tmp_g2, nrx_);
        rzf().setArrayBV(tmp_g2+nrx_, nrz_);

        // Get influence on the backward quadratures
        if (nrq_>0) {
          g_.input(RDAE_RX).setBV(rxf());
          g_.input(RDAE_RZ).setBV(rzf());
          g_.spEvaluate(true);
          g_.output(RDAE_QUAD).getBV(rqf());
        }
      }
    } else { // reverse mode

      // Clear adjoint sensitivities
      x0().setBV(xf());
      p().setZeroBV();
      z0().setBV(zf());
      rx0().setBV(rxf());
      rp().setZeroBV();
      rz0().setBV(rzf());

      // Propagate through g
      if (!g_.isNull()) {

        // Propagate influence from the backward quadratures
        if (nrq_>0) {
          g_.output(RDAE_ODE).setZeroBV();
          g_.output(RDAE_ALG).setZeroBV();
          g_.output(RDAE_QUAD).setBV(rqf());
          g_.spEvaluate(false);
          rx0().borBV(g_.input(RDAE_RX));
          rz0().borBV(g_.input(RDAE_RZ));
        }

        // Propagate interdependencies
        rx0().getArrayBV(tmp_g2, nrx_);
        rz0().getArrayBV(tmp_g2+nrx_, nrz_);
        std::fill(tmp_g1, tmp_g1+nrx_+nrz_, 0);
        linsol_g_.spSolve(tmp_g1, tmp_g2, false);

        // Propagate through the backward DAE
        g_.output(RDAE_ODE).setArrayBV(tmp_g1, nrx_);
        g_.output(RDAE_ODE).borBV(rx0());
        g_.output(RDAE_ALG).setArrayBV(tmp_g1+nrx_, nrz_);
        g_.output(RDAE_ALG).borBV(rz0());
        g_.spEvaluate(false);
        x0().borBV(g_.input(RDAE_X));
        p().borBV(g_.input(RDAE_P));
        z0().borBV(g_.input(RDAE_Z));
        rx0().borBV(g_.input(RDAE_RX));
        rp().borBV(g_.input(RDAE_RP));

        // No dependency on initial guess
        rz0().setZeroBV();
      }

      // Propagate influence from the quadratures
      if (nq_>0) {
        f_.output(DAE_ODE).setZeroBV();
        f_.output(DAE_ALG).setZeroBV();
        f_.output(DAE_QUAD).setBV(qf());
        f_.spEvaluate(false);
        x0().borBV(f_.input(DAE_X));
        z0().borBV(f_.input(DAE_Z));
      }

      // Propagate interdependencies
      x0().getArrayBV(tmp_f2, nx_);
      z0().getArrayBV(tmp_f2+nx_, nz_);
      std::fill(tmp_f1, tmp_f1+nx_+nz_, 0);
      linsol_f_.spSolve(tmp_f1, tmp_f2, false);

      // Propagate through the DAE
      f_.output(DAE_ODE).setArrayBV(tmp_f1, nx_);
      f_.output(DAE_ODE).borBV(x0());
      f_.output(DAE_ALG).setArrayBV(tmp_f1+nx_, nz_);
      f_.output(DAE_ALG).borBV(z0());
      f_.spEvaluate(false);
      x0().borBV(f_.input(DAE_X));
      p().borBV(f_.input(DAE_P));

      // No dependency on initial guess
      z0().setZeroBV();
    }

    log("IntegratorInternal::spEvaluate", "end");
  }

  IntegratorInternal::AugOffset IntegratorInternal::getAugOffset(int nfwd, int nadj) {
    // Form return object
    AugOffset ret;
    ret.x.resize(1, 0);
    ret.z.resize(1, 0);
    ret.q.resize(1, 0);
    ret.p.resize(1, 0);
    ret.rx.resize(1, 0);
    ret.rz.resize(1, 0);
    ret.rq.resize(1, 0);
    ret.rp.resize(1, 0);

    // Count nondifferentiated and forward sensitivities
    for (int dir=-1; dir<nfwd; ++dir) {
      if ( nx_>0) ret.x.push_back(x0().size2());
      if ( nz_>0) ret.z.push_back(z0().size2());
      if ( nq_>0) ret.q.push_back(qf().size2());
      if ( np_>0) ret.p.push_back(p().size2());
      if (nrx_>0) ret.rx.push_back(rx0().size2());
      if (nrz_>0) ret.rz.push_back(rz0().size2());
      if (nrq_>0) ret.rq.push_back(rqf().size2());
      if (nrp_>0) ret.rp.push_back(rp().size2());
    }

    // Count adjoint sensitivities
    for (int dir=0; dir<nadj; ++dir) {
      if ( nx_>0) ret.rx.push_back(x0().size2());
      if ( nz_>0) ret.rz.push_back(z0().size2());
      if ( np_>0) ret.rq.push_back(p().size2());
      if ( nq_>0) ret.rp.push_back(qf().size2());
      if (nrx_>0) ret.x.push_back(rx0().size2());
      if (nrz_>0) ret.z.push_back(rz0().size2());
      if (nrp_>0) ret.q.push_back(rp().size2());
      if (nrq_>0) ret.p.push_back(rqf().size2());
    }

    // Get cummulative offsets
    for (int i=1; i<ret.x.size(); ++i) ret.x[i] += ret.x[i-1];
    for (int i=1; i<ret.z.size(); ++i) ret.z[i] += ret.z[i-1];
    for (int i=1; i<ret.q.size(); ++i) ret.q[i] += ret.q[i-1];
    for (int i=1; i<ret.p.size(); ++i) ret.p[i] += ret.p[i-1];
    for (int i=1; i<ret.rx.size(); ++i) ret.rx[i] += ret.rx[i-1];
    for (int i=1; i<ret.rz.size(); ++i) ret.rz[i] += ret.rz[i-1];
    for (int i=1; i<ret.rq.size(); ++i) ret.rq[i] += ret.rq[i-1];
    for (int i=1; i<ret.rp.size(); ++i) ret.rp[i] += ret.rp[i-1];

    // Return the offsets
    return ret;
  }

  Function IntegratorInternal::getDerivative(int nfwd, int nadj) {
    log("IntegratorInternal::getDerivative", "begin");

    // Form the augmented DAE
    AugOffset offset;
    std::pair<Function, Function> aug_dae = getAugmented(nfwd, nadj, offset);

    // Create integrator for augmented DAE
    Integrator integrator;
    integrator.assignNode(create(aug_dae.first, aug_dae.second));

    // Set solver specific options
    setDerivativeOptions(integrator, offset);

    // Pass down specific options if provided
    if (hasSetOption("augmented_options"))
      integrator.setOption(getOption("augmented_options"));

    // Initialize the integrator since we will call it below
    integrator.init();

    // All inputs of the return function
    vector<MX> ret_in;
    ret_in.reserve(INTEGRATOR_NUM_IN*(1+nfwd) + INTEGRATOR_NUM_OUT*nadj);

    // Augmented state
    MX x0_aug, p_aug, z0_aug, rx0_aug, rp_aug, rz0_aug;

    // Temp stringstream
    stringstream ss;

    // Inputs or forward/adjoint seeds in one direction
    vector<MX> dd;

    // Add nondifferentiated inputs and forward seeds
    dd.resize(INTEGRATOR_NUM_IN);
    for (int dir=-1; dir<nfwd; ++dir) {

      // Differential state
      ss.clear();
      ss << "x0";
      if (dir>=0) ss << "_" << dir;
      dd[INTEGRATOR_X0] = MX::sym(ss.str(), x0().sparsity());
      x0_aug.appendColumns(dd[INTEGRATOR_X0]);

      // Parameter
      ss.clear();
      ss << "p";
      if (dir>=0) ss << "_" << dir;
      dd[INTEGRATOR_P] = MX::sym(ss.str(), p().sparsity());
      p_aug.appendColumns(dd[INTEGRATOR_P]);

      // Initial guess for algebraic variable
      ss.clear();
      ss << "r0";
      if (dir>=0) ss << "_" << dir;
      dd[INTEGRATOR_Z0] = MX::sym(ss.str(), z0().sparsity());
      z0_aug.appendColumns(dd[INTEGRATOR_Z0]);

      // Backward state
      ss.clear();
      ss << "rx0";
      if (dir>=0) ss << "_" << dir;
      dd[INTEGRATOR_RX0] = MX::sym(ss.str(), rx0().sparsity());
      rx0_aug.appendColumns(dd[INTEGRATOR_RX0]);

      // Backward parameter
      ss.clear();
      ss << "rp";
      if (dir>=0) ss << "_" << dir;
      dd[INTEGRATOR_RP] = MX::sym(ss.str(), rp().sparsity());
      rp_aug.appendColumns(dd[INTEGRATOR_RP]);

      // Initial guess for backward algebraic variable
      ss.clear();
      ss << "rz0";
      if (dir>=0) ss << "_" << dir;
      dd[INTEGRATOR_RZ0] = MX::sym(ss.str(), rz0().sparsity());
      rz0_aug.appendColumns(dd[INTEGRATOR_RZ0]);

      // Add to input vector
      ret_in.insert(ret_in.end(), dd.begin(), dd.end());
    }

    // Add adjoint seeds
    dd.resize(INTEGRATOR_NUM_OUT);
    for (int dir=0; dir<nadj; ++dir) {

      // Differential states become backward differential state
      ss.clear();
      ss << "xf" << "_" << dir;
      dd[INTEGRATOR_XF] = MX::sym(ss.str(), xf().sparsity());
      rx0_aug.appendColumns(dd[INTEGRATOR_XF]);

      // Quadratures become backward parameters
      ss.clear();
      ss << "qf" << "_" << dir;
      dd[INTEGRATOR_QF] = MX::sym(ss.str(), qf().sparsity());
      rp_aug.appendColumns(dd[INTEGRATOR_QF]);

      // Algebraic variables become backward algebraic variables
      ss.clear();
      ss << "zf" << "_" << dir;
      dd[INTEGRATOR_ZF] = MX::sym(ss.str(), zf().sparsity());
      rz0_aug.appendColumns(dd[INTEGRATOR_ZF]);

      // Backward differential states becomes forward differential states
      ss.clear();
      ss << "rxf" << "_" << dir;
      dd[INTEGRATOR_RXF] = MX::sym(ss.str(), rxf().sparsity());
      x0_aug.appendColumns(dd[INTEGRATOR_RXF]);

      // Backward quadratures becomes (forward) parameters
      ss.clear();
      ss << "rqf" << "_" << dir;
      dd[INTEGRATOR_RQF] = MX::sym(ss.str(), rqf().sparsity());
      p_aug.appendColumns(dd[INTEGRATOR_RQF]);

      // Backward differential states becomes forward differential states
      ss.clear();
      ss << "rzf" << "_" << dir;
      dd[INTEGRATOR_RZF] = MX::sym(ss.str(), rzf().sparsity());
      z0_aug.appendColumns(dd[INTEGRATOR_RZF]);

      // Add to input vector
      ret_in.insert(ret_in.end(), dd.begin(), dd.end());
    }

    // Call the integrator
    vector<MX> integrator_in(INTEGRATOR_NUM_IN);
    integrator_in[INTEGRATOR_X0] = x0_aug;
    integrator_in[INTEGRATOR_P] = p_aug;
    integrator_in[INTEGRATOR_Z0] = z0_aug;
    integrator_in[INTEGRATOR_RX0] = rx0_aug;
    integrator_in[INTEGRATOR_RP] = rp_aug;
    integrator_in[INTEGRATOR_RZ0] = rz0_aug;
    vector<MX> integrator_out = integrator.call(integrator_in);

    // Augmented results
    vector<MX> xf_aug = horzsplit(integrator_out[INTEGRATOR_XF], offset.x);
    vector<MX> qf_aug = horzsplit(integrator_out[INTEGRATOR_QF], offset.q);
    vector<MX> zf_aug = horzsplit(integrator_out[INTEGRATOR_ZF], offset.z);
    vector<MX> rxf_aug = horzsplit(integrator_out[INTEGRATOR_RXF], offset.rx);
    vector<MX> rqf_aug = horzsplit(integrator_out[INTEGRATOR_RQF], offset.rq);
    vector<MX> rzf_aug = horzsplit(integrator_out[INTEGRATOR_RZF], offset.rz);
    vector<MX>::const_iterator xf_aug_it = xf_aug.begin();
    vector<MX>::const_iterator qf_aug_it = qf_aug.begin();
    vector<MX>::const_iterator zf_aug_it = zf_aug.begin();
    vector<MX>::const_iterator rxf_aug_it = rxf_aug.begin();
    vector<MX>::const_iterator rqf_aug_it = rqf_aug.begin();
    vector<MX>::const_iterator rzf_aug_it = rzf_aug.begin();

    // All outputs of the return function
    vector<MX> ret_out;
    ret_out.reserve(INTEGRATOR_NUM_OUT*(1+nfwd) + INTEGRATOR_NUM_IN*nadj);

    // Collect the nondifferentiated results and forward sensitivities
    dd.resize(INTEGRATOR_NUM_OUT);
    fill(dd.begin(), dd.end(), MX());
    for (int dir=-1; dir<nfwd; ++dir) {
      if ( nx_>0) dd[INTEGRATOR_XF]  = *xf_aug_it++;
      if ( nq_>0) dd[INTEGRATOR_QF]  = *qf_aug_it++;
      if ( nz_>0) dd[INTEGRATOR_ZF]  = *zf_aug_it++;
      if (nrx_>0) dd[INTEGRATOR_RXF] = *rxf_aug_it++;
      if (nrq_>0) dd[INTEGRATOR_RQF] = *rqf_aug_it++;
      if (nrz_>0) dd[INTEGRATOR_RZF] = *rzf_aug_it++;
      ret_out.insert(ret_out.end(), dd.begin(), dd.end());
    }

    // Collect the adjoint sensitivities
    dd.resize(INTEGRATOR_NUM_IN);
    fill(dd.begin(), dd.end(), MX());
    for (int dir=0; dir<nadj; ++dir) {
      if ( nx_>0) dd[INTEGRATOR_X0]  = *rxf_aug_it++;
      if ( np_>0) dd[INTEGRATOR_P]   = *rqf_aug_it++;
      if ( nz_>0) dd[INTEGRATOR_Z0]  = *rzf_aug_it++;
      if (nrx_>0) dd[INTEGRATOR_RX0] = *xf_aug_it++;
      if (nrp_>0) dd[INTEGRATOR_RP]  = *qf_aug_it++;
      if (nrz_>0) dd[INTEGRATOR_RZ0] = *zf_aug_it++;
      ret_out.insert(ret_out.end(), dd.begin(), dd.end());
    }
    log("IntegratorInternal::getDerivative", "end");

    // Create derivative function and return
    return MXFunction(ret_in, ret_out);
  }

  Function IntegratorInternal::getJacobian(int iind, int oind, bool compact, bool symmetric) {
    vector<MX> arg = symbolicInput();
    vector<MX> res = shared_from_this<Function>().call(arg);
    MXFunction f(arg, res);
    f.setOption("ad_mode", "forward");
    f.init();
    return f.jacobian(iind, oind, compact, symmetric);
  }

  void IntegratorInternal::reset() {
    log("IntegratorInternal::reset", "begin");

    // Go to the start time
    t_ = t0_;

    // Initialize output
    xf().set(x0());
    zf().set(z0());

    // Reset summation states
    qf().set(0.0);

    log("IntegratorInternal::reset", "end");
  }

  void IntegratorInternal::resetB() {
    log("IntegratorInternal::resetB", "begin");

    // Go to the end time
    t_ = tf_;

    // Initialize output
    rxf().set(rx0());
    rzf().set(rz0());

    // Reset summation states
    rqf().set(0.0);

    log("IntegratorInternal::resetB", "end");
  }

  void IntegratorInternal::setDerivativeOptions(Integrator& integrator, const AugOffset& offset) {
    // Copy all options
    integrator.setOption(dictionary());
  }

  Sparsity IntegratorInternal::spJacF() {
    // Start with the sparsity pattern of the ODE part
    Sparsity ret = f_.jacSparsity(DAE_X, DAE_ODE).T();

    // Add diagonal to get interdependencies
    ret = ret.patternUnion(Sparsity::diag(nx_));

    // Quick return if no algebraic variables
    if (nz_==0) return ret;

    // Add contribution from algebraic variables and equations
    Sparsity jac_ode_z = f_.jacSparsity(DAE_Z, DAE_ODE).T();
    Sparsity jac_alg_x = f_.jacSparsity(DAE_X, DAE_ALG).T();
    Sparsity jac_alg_z = f_.jacSparsity(DAE_Z, DAE_ALG).T();
    ret = vertcat(ret, jac_ode_z);
    ret.appendColumns(vertcat(jac_alg_x, jac_alg_z));
    return ret;
  }

  Sparsity IntegratorInternal::spJacG() {
    // Start with the sparsity pattern of the ODE part
    Sparsity ret = g_.jacSparsity(RDAE_RX, RDAE_ODE).T();

    // Add diagonal to get interdependencies
    ret = ret.patternUnion(Sparsity::diag(nrx_));

    // Quick return if no algebraic variables
    if (nrz_==0) return ret;

    // Add contribution from algebraic variables and equations
    Sparsity jac_ode_z = g_.jacSparsity(RDAE_RZ, RDAE_ODE).T();
    Sparsity jac_alg_x = g_.jacSparsity(RDAE_RX, RDAE_ALG).T();
    Sparsity jac_alg_z = g_.jacSparsity(RDAE_RZ, RDAE_ALG).T();
    ret = vertcat(ret, jac_ode_z);
    ret.appendColumns(vertcat(jac_alg_x, jac_alg_z));
    return ret;
  }

  std::map<std::string, IntegratorInternal::Plugin> IntegratorInternal::solvers_;

  const std::string IntegratorInternal::infix_ = "integrator";

  void IntegratorInternal::setStopTime(double tf) {
    casadi_error("IntegratorInternal::setStopTime not defined for class "
                 << typeid(*this).name());
  }

} // namespace casadi



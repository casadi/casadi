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


#include "integrator_impl.hpp"
#include "../std_vector_tools.hpp"

using namespace std;
namespace casadi {

  bool has_integrator(const string& name) {
    return Integrator::hasPlugin(name);
  }

  void load_integrator(const string& name) {
    Integrator::loadPlugin(name);
  }

  string doc_integrator(const string& name) {
    return Integrator::getPlugin(name).doc;
  }

  Function integrator(const string& name, const string& solver,
                            const SXDict& dae, const Dict& opts) {
    return integrator(name, solver, Integrator::map2problem(dae), opts);
  }

  Function integrator(const string& name, const string& solver,
                            const MXDict& dae, const Dict& opts) {
    return integrator(name, solver, Integrator::map2problem(dae), opts);
  }

  Function integrator(const string& name, const string& solver,
                            const Function& dae, const Dict& opts) {
    if (dae.is_a("sxfunction")) {
      SXProblem p = Integrator::fun2problem<SX>(dae);
      return integrator(name, solver, p, opts);
    } else {
      MXProblem p = Integrator::fun2problem<MX>(dae);
      return integrator(name, solver, p, opts);
    }
  }

  Function integrator(const string& name, const string& solver,
                                const pair<Function, Function>& dae,
                                const Dict& opts) {
    if (dae.first.is_a("sxfunction")) {
      SXProblem p = Integrator::fun2problem<SX>(dae.first, dae.second);
      return integrator(name, solver, p, opts);
    } else {
      MXProblem p = Integrator::fun2problem<MX>(dae.first, dae.second);
      return integrator(name, solver, p, opts);
    }
  }

  Function integrator(const string& name, const string& solver,
                                const XProblem& dae, const Dict& opts) {
    Function ret;
    ret.assignNode(Integrator::getPlugin(solver).creator(name, dae));
    ret->construct(opts);
    return ret;
  }

  Function Function::integrator_dae() {
    casadi_assert(!is_null());
    Integrator* n = dynamic_cast<Integrator*>(get());
    casadi_assert_message(n!=0, "Not an integrator");
    return n->f_;
  }

  vector<string> integrator_in() {
    vector<string> ret(integrator_n_in());
    for (size_t i=0; i<ret.size(); ++i) ret[i]=integrator_in(i);
    return ret;
  }

  vector<string> integrator_out() {
    vector<string> ret(integrator_n_out());
    for (size_t i=0; i<ret.size(); ++i) ret[i]=integrator_out(i);
    return ret;
  }

  string integrator_in(int ind) {
    switch (static_cast<IntegratorInput>(ind)) {
    case INTEGRATOR_X0:  return "x0";
    case INTEGRATOR_P:   return "p";
    case INTEGRATOR_Z0:  return "z0";
    case INTEGRATOR_RX0: return "rx0";
    case INTEGRATOR_RP:  return "rp";
    case INTEGRATOR_RZ0: return "rz0";
    case INTEGRATOR_NUM_IN: break;
    }
    return string();
  }

  string integrator_out(int ind) {
    switch (static_cast<IntegratorOutput>(ind)) {
    case INTEGRATOR_XF:  return "xf";
    case INTEGRATOR_QF:  return "qf";
    case INTEGRATOR_ZF:  return "zf";
    case INTEGRATOR_RXF: return "rxf";
    case INTEGRATOR_RQF: return "rqf";
    case INTEGRATOR_RZF: return "rzf";
    case INTEGRATOR_NUM_OUT: break;
    }
    return string();
  }

  int integrator_n_in() {
    return INTEGRATOR_NUM_IN;
  }

  int integrator_n_out() {
    return INTEGRATOR_NUM_OUT;
  }

  Integrator::Integrator(const std::string& name, const XProblem& dae)
    : FunctionInternal(name), dae_(dae) {

    if (dae.is_sx) {
      f_ = get_f<SX>();
      g_ = get_g<SX>();
    } else {
      f_ = get_f<MX>();
      g_ = get_g<MX>();
    }

    // Negative number of parameters for consistancy checking
    np_ = -1;

    // Get the sparsities
    t_ = f_.sparsity_in(DAE_T);
    x_ = f_.sparsity_in(DAE_X);
    z_ = f_.sparsity_in(DAE_Z);
    p_ = f_.sparsity_in(DAE_P);
    q_ = f_.sparsity_out(DAE_QUAD);
    if (!g_.is_null()) {
      rx_ = g_.sparsity_in(RDAE_RX);
      rz_ = g_.sparsity_in(RDAE_RZ);
      rp_ = g_.sparsity_in(RDAE_RP);
      rq_ = g_.sparsity_out(RDAE_QUAD);
    }

    // Default options
    print_stats_ = false;
    output_t0_ = false;
  }

  Integrator::~Integrator() {
  }

  template<typename MatType>
  Function Integrator::get_f() const {
    vector<MatType> dae_in(DAE_NUM_IN), dae_out(DAE_NUM_OUT);
    const Problem<MatType>& dae = this->dae_;
    dae_in[DAE_T]=dae.in[DE_T];
    dae_in[DAE_X]=dae.in[DE_X];
    dae_in[DAE_Z]=dae.in[DE_Z];
    dae_in[DAE_P]=dae.in[DE_P];
    dae_out[DAE_ODE]=dae.out[DE_ODE];
    dae_out[DAE_ALG]=dae.out[DE_ALG];
    dae_out[DAE_QUAD]=dae.out[DE_QUAD];
    return Function("dae", dae_in, dae_out);
  }

  template<typename MatType>
  Function Integrator::get_g() const {
    vector<MatType> rdae_in(RDAE_NUM_IN), rdae_out(RDAE_NUM_OUT);
    const Problem<MatType>& dae = this->dae_;
    rdae_in[RDAE_T]=dae.in[DE_T];
    rdae_in[RDAE_X]=dae.in[DE_X];
    rdae_in[RDAE_Z]=dae.in[DE_Z];
    rdae_in[RDAE_P]=dae.in[DE_P];
    rdae_in[RDAE_RX]=dae.in[DE_RX];
    rdae_in[RDAE_RZ]=dae.in[DE_RZ];
    rdae_in[RDAE_RP]=dae.in[DE_RP];
    rdae_out[RDAE_ODE]=dae.out[DE_RODE];
    rdae_out[RDAE_ALG]=dae.out[DE_RALG];
    rdae_out[RDAE_QUAD]=dae.out[DE_RQUAD];
    return Function("rdae", rdae_in, rdae_out);
  }

  Sparsity Integrator::get_sparsity_in(int ind) const {
    switch (static_cast<IntegratorInput>(ind)) {
    case INTEGRATOR_X0: return x_;
    case INTEGRATOR_P: return p_;
    case INTEGRATOR_Z0: return z_;
    case INTEGRATOR_RX0: return repmat(rx_, 1, ntout_);
    case INTEGRATOR_RP: return repmat(rp_, 1, ntout_);
    case INTEGRATOR_RZ0: return repmat(rz_, 1, ntout_);
    case INTEGRATOR_NUM_IN: break;
    }
    return Sparsity();
  }

  Sparsity Integrator::get_sparsity_out(int ind) const {
    switch (static_cast<IntegratorOutput>(ind)) {
    case INTEGRATOR_XF: return repmat(x_, 1, ntout_);
    case INTEGRATOR_QF: return repmat(q_, 1, ntout_);
    case INTEGRATOR_ZF: return repmat(z_, 1, ntout_);
    case INTEGRATOR_RXF: return rx_;
    case INTEGRATOR_RQF: return rq_;
    case INTEGRATOR_RZF: return rz_;
    case INTEGRATOR_NUM_OUT: break;
    }
    return Sparsity();
  }

  void Integrator::
  eval(Memory& mem, const double** arg, double** res, int* iw, double* w) const {
    IntegratorMemory& m = dynamic_cast<IntegratorMemory&>(mem);
    Integrator* this_ = const_cast<Integrator*>(this);

    // Read inputs
    const double* x0 = arg[INTEGRATOR_X0];
    const double* z0 = arg[INTEGRATOR_Z0];
    const double* p = arg[INTEGRATOR_P];
    const double* rx0 = arg[INTEGRATOR_RX0];
    const double* rz0 = arg[INTEGRATOR_RZ0];
    const double* rp = arg[INTEGRATOR_RP];
    arg += INTEGRATOR_NUM_IN;

    // Read outputs
    double* x = res[INTEGRATOR_XF];
    double* z = res[INTEGRATOR_ZF];
    double* q = res[INTEGRATOR_QF];
    double* rx = res[INTEGRATOR_RXF];
    double* rz = res[INTEGRATOR_RZF];
    double* rq = res[INTEGRATOR_RQF];
    res += INTEGRATOR_NUM_OUT;

    // Setup memory object
    setup(m, arg, res, iw, w);

    // Reset solver, take time to t0
    reset(m, grid_.front(), x0, z0, p);

    // Integrate forward
    for (int k=0; k<grid_.size(); ++k) {
      // Skip t0?
      if (k==0 && !output_t0_) continue;

      // Integrate forward
      this_->advance(m, grid_[k], x, z, q);
      if (x) x += x_.nnz();
      if (z) z += z_.nnz();
      if (q) q += q_.nnz();
    }

    // If backwards integration is needed
    if (nrx_>0) {
      // Integrate backward
      this_->resetB(m, grid_.back(), rx0, rz0, rp);

      // Proceed to t0
      this_->retreat(m, grid_.front(), rx, rz, rq);
    }

    // Print statistics
    if (print_stats_) printStats(m, userOut());
  }

  Options Integrator::options_
  = {{&FunctionInternal::options_},
     {{"print_stats",
       {OT_BOOL,
        "Print out statistics after integration"}},
      {"t0",
       {OT_DOUBLE,
        "Beginning of the time horizon"}},
      {"tf",
       {OT_DOUBLE,
        "End of the time horizon"}},
      {"grid",
       {OT_DOUBLEVECTOR,
        "Time grid"}},
      {"augmented_options",
       {OT_DICT,
        "Options to be passed down to the augmented integrator, if one is constructed."}},
      {"output_t0",
       {OT_BOOL,
        "Output the state at the initial time"}}
     }
  };

  void Integrator::init(const Dict& opts) {
    // Default (temporary) options
    double t0=0, tf=1;

    // Read options
    for (auto&& op : opts) {
      if (op.first=="output_t0") {
        output_t0_ = op.second;
      } else if (op.first=="print_stats") {
        print_stats_ = op.second;
      } else if (op.first=="grid") {
        grid_ = op.second;
      } else if (op.first=="augmented_options") {
        augmented_options_ = op.second;
      } else if (op.first=="t0") {
        t0 = op.second;
      } else if (op.first=="tf") {
        tf = op.second;
      }
    }

    // Store a copy of the options, for creating augmented integrators
    opts_ = opts;

    // If grid unset, default to [t0, tf]
    if (grid_.empty()) {
      grid_ = {t0, tf};
    }

    ngrid_ = grid_.size();
    ntout_ = output_t0_ ? ngrid_ : ngrid_-1;

    // Call the base class method
    FunctionInternal::init(opts);

    // Get dimensions
    nx_ = x().nnz();
    nz_ = z().nnz();
    nq_ = q().nnz();
    np_  = p().nnz();
    nrx_ = rx().nnz();
    nrz_ = rz().nnz();
    nrp_ = rp().nnz();
    nrq_ = rq().nnz();

    // Consistency checks
    casadi_assert(!f_.is_null());
    casadi_assert_message(f_.n_in()==DAE_NUM_IN,
                          "Wrong number of inputs for the DAE callback function");
    casadi_assert_message(f_.n_out()==DAE_NUM_OUT,
                          "Wrong number of outputs for the DAE callback function");
    if (!g_.is_null()) {
      casadi_assert_message(g_.n_in()==RDAE_NUM_IN,
                            "Wrong number of inputs for the backwards DAE callback function");
      casadi_assert_message(g_.n_out()==RDAE_NUM_OUT,
                            "Wrong number of outputs for the backwards DAE callback function");
    }
    casadi_assert_message(f_.size_out(DAE_ODE)==x().size(),
                          "Inconsistent dimensions. Expecting DAE_ODE output of shape "
                          << size_in(INTEGRATOR_X0) << ", but got "
                          << f_.size_out(DAE_ODE) << " instead.");
    casadi_assert(f_.sparsity_out(DAE_ODE)==x());
    casadi_assert_message(f_.size_out(DAE_ALG)==z().size(),
                          "Inconsistent dimensions. Expecting DAE_ALG output of shape "
                          << z().size() << ", but got "
                          << f_.size_out(DAE_ALG) << " instead.");
    casadi_assert(f_.sparsity_out(DAE_ALG)==z());
    if (!g_.is_null()) {
      casadi_assert(g_.sparsity_in(RDAE_P)==p());
      casadi_assert(g_.sparsity_in(RDAE_X)==x());
      casadi_assert(g_.sparsity_in(RDAE_Z)==z());
      casadi_assert(g_.sparsity_out(RDAE_ODE)==rx());
      casadi_assert(g_.sparsity_out(RDAE_ALG)==rz());
    }

    // Warn if sparse inputs (was previously an error)
    casadi_assert_warning(f_.sparsity_in(DAE_X).is_dense(),
                          "Sparse states in integrators are experimental");

    // Form a linear solver for the sparsity propagation
    linsol_f_ = linsol("linsol_f", "none", spJacF(), 1);
    if (!g_.is_null()) {
      linsol_g_ = linsol("linsol_g", "none", spJacG(), 1);
    }

    // Allocate sufficiently large work vectors
    size_t sz_w = f_.sz_w();
    alloc(f_);
    if (!g_.is_null()) {
      alloc(g_);
      sz_w = max(sz_w, g_.sz_w());
    }
    sz_w = max(sz_w, static_cast<size_t>(nx_+nz_));
    sz_w = max(sz_w, static_cast<size_t>(nrx_+nrz_));
    alloc_w(sz_w + nx_ + nz_ + nrx_ + nrz_);
  }

  void Integrator::init_memory(Memory& mem) const {
  }

  template<typename MatType>
  map<string, MatType> Integrator::aug_fwd(int nfwd, AugOffset& offset) {
    log("Integrator::aug_fwd", "call");

    // Return object
    map<string, MatType> ret;

    // Calculate offsets
    offset = getAugOffset(nfwd, 0);

    // Create augmented problem
    MatType aug_t = MatType::sym("aug_t", t());
    MatType aug_x = MatType::sym("aug_x", x().size1(), offset.x.back());
    MatType aug_z = MatType::sym("aug_z", std::max(z().size1(), rz().size1()), offset.z.back());
    MatType aug_p = MatType::sym("aug_p", std::max(p().size1(), rp().size1()), offset.p.back());
    MatType aug_rx = MatType::sym("aug_rx", x().size1(), offset.rx.back());
    MatType aug_rz = MatType::sym("aug_rz", std::max(z().size1(), rz().size1()), offset.rz.back());
    MatType aug_rp = MatType::sym("aug_rp", std::max(q().size1(), rp().size1()), offset.rp.back());

    // Split up the augmented vectors
    vector<MatType> aug_x_split = horzsplit(aug_x, offset.x);
    auto aug_x_split_it = aug_x_split.begin();
    vector<MatType> aug_z_split = horzsplit(aug_z, offset.z);
    auto aug_z_split_it = aug_z_split.begin();
    vector<MatType> aug_p_split = horzsplit(aug_p, offset.p);
    auto aug_p_split_it = aug_p_split.begin();
    vector<MatType> aug_rx_split = horzsplit(aug_rx, offset.rx);
    auto aug_rx_split_it = aug_rx_split.begin();
    vector<MatType> aug_rz_split = horzsplit(aug_rz, offset.rz);
    auto aug_rz_split_it = aug_rz_split.begin();
    vector<MatType> aug_rp_split = horzsplit(aug_rp, offset.rp);
    auto aug_rp_split_it = aug_rp_split.begin();

    // Temporary vector
    vector<MatType> tmp;

    // Zero with the dimension of t
    MatType zero_t = DM::zeros(aug_t.sparsity());

    // The DAE being constructed
    vector<MatType> f_ode, f_alg, f_quad, g_ode, g_alg, g_quad;

    // Forward derivatives of f
    Function d = f_.derivative(nfwd, 0);
    vector<MatType> f_arg;
    f_arg.reserve(d.n_in());
    tmp.resize(DAE_NUM_IN);
    fill(tmp.begin(), tmp.end(), MatType());

    // Collect arguments for calling d
    for (int dir=-1; dir<nfwd; ++dir) {
      tmp[DAE_T] = dir<0 ? aug_t : zero_t;
      if ( nx_>0) tmp[DAE_X] = *aug_x_split_it++;
      if ( nz_>0) tmp[DAE_Z] = *aug_z_split_it++;
      if ( np_>0) tmp[DAE_P] = *aug_p_split_it++;
      f_arg.insert(f_arg.end(), tmp.begin(), tmp.end());
    }

    // Call d
    vector<MatType> res = d(f_arg);
    auto res_it = res.begin();

    // Collect right-hand-sides
    tmp.resize(DAE_NUM_OUT);
    fill(tmp.begin(), tmp.end(), MatType());
    for (int dir=-1; dir<nfwd; ++dir) {
      copy(res_it, res_it+tmp.size(), tmp.begin());
      res_it += tmp.size();
      if ( nx_>0) f_ode.push_back(tmp[DAE_ODE]);
      if ( nz_>0) f_alg.push_back(tmp[DAE_ALG]);
      if ( nq_>0) f_quad.push_back(tmp[DAE_QUAD]);
    }

    // Consistency check
    casadi_assert(res_it==res.end());

    vector<MatType> g_arg;
    if (!g_.is_null()) {

      // Forward derivatives of g
      d = g_.derivative(nfwd, 0);
      g_arg.reserve(d.n_in());
      tmp.resize(RDAE_NUM_IN);
      fill(tmp.begin(), tmp.end(), MatType());

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
      res = d(g_arg);
      res_it = res.begin();

      // Collect right-hand-sides
      tmp.resize(RDAE_NUM_OUT);
      fill(tmp.begin(), tmp.end(), MatType());
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

    // Form the augmented forward integration
    ret["t"] = aug_t;
    ret["x"] = aug_x;
    ret["z"] = aug_z;
    ret["p"] = aug_p;
    if (!f_ode.empty()) ret["ode"] = densify(horzcat(f_ode));
    if (!f_alg.empty()) ret["alg"] = densify(horzcat(f_alg));
    if (!f_quad.empty()) ret["quad"] = densify(horzcat(f_quad));

    // Form the augmented backward integration function
    if (!g_ode.empty()) {
      ret["rx"] = aug_rx;
      ret["rz"] = aug_rz;
      ret["rp"] = aug_rp;
      if (!g_ode.empty()) ret["rode"] = densify(horzcat(g_ode));
      if (!g_alg.empty()) ret["ralg"] = densify(horzcat(g_alg));
      if (!g_quad.empty()) ret["rquad"] = densify(horzcat(g_quad));
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

  template<typename MatType>
  map<string, MatType> Integrator::aug_adj(int nadj, AugOffset& offset) {
    log("Integrator::aug_adj", "call");

    // Return object
    map<string, MatType> ret;

    // Calculate offsets
    offset = getAugOffset(0, nadj);

    // Create augmented problem
    MatType aug_t = MatType::sym("aug_t", t());
    MatType aug_x = MatType::sym("aug_x", x().size1(), offset.x.back());
    MatType aug_z = MatType::sym("aug_z", std::max(z().size1(), rz().size1()), offset.z.back());
    MatType aug_p = MatType::sym("aug_p", std::max(p().size1(), rp().size1()), offset.p.back());
    MatType aug_rx = MatType::sym("aug_rx", size1_in(INTEGRATOR_X0), offset.rx.back());
    MatType aug_rz = MatType::sym("aug_rz", std::max(z().size1(), rz().size1()), offset.rz.back());
    MatType aug_rp = MatType::sym("aug_rp", std::max(q().size1(), rp().size1()), offset.rp.back());

    // Split up the augmented vectors
    vector<MatType> aug_x_split = horzsplit(aug_x, offset.x);
    auto aug_x_split_it = aug_x_split.begin();
    vector<MatType> aug_z_split = horzsplit(aug_z, offset.z);
    auto aug_z_split_it = aug_z_split.begin();
    vector<MatType> aug_p_split = horzsplit(aug_p, offset.p);
    auto aug_p_split_it = aug_p_split.begin();
    vector<MatType> aug_rx_split = horzsplit(aug_rx, offset.rx);
    auto aug_rx_split_it = aug_rx_split.begin();
    vector<MatType> aug_rz_split = horzsplit(aug_rz, offset.rz);
    auto aug_rz_split_it = aug_rz_split.begin();
    vector<MatType> aug_rp_split = horzsplit(aug_rp, offset.rp);
    auto aug_rp_split_it = aug_rp_split.begin();

    // Temporary vector
    vector<MatType> tmp;

    // Zero with the dimension of t
    MatType zero_t = DM::zeros(aug_t.sparsity());

    // The DAE being constructed
    vector<MatType> f_ode, f_alg, f_quad, g_ode, g_alg, g_quad;

    // Forward derivatives of f
    Function d = f_.derivative(0, 0);
    vector<MatType> f_arg;
    f_arg.reserve(d.n_in());
    tmp.resize(DAE_NUM_IN);
    fill(tmp.begin(), tmp.end(), MatType());

    // Collect arguments for calling d
    tmp[DAE_T] = aug_t;
    if ( nx_>0) tmp[DAE_X] = *aug_x_split_it++;
    if ( nz_>0) tmp[DAE_Z] = *aug_z_split_it++;
    if ( np_>0) tmp[DAE_P] = *aug_p_split_it++;
    f_arg.insert(f_arg.end(), tmp.begin(), tmp.end());

    // Call d
    vector<MatType> res = d(f_arg);
    auto res_it = res.begin();

    // Collect right-hand-sides
    tmp.resize(DAE_NUM_OUT);
    fill(tmp.begin(), tmp.end(), MatType());
    copy(res_it, res_it+tmp.size(), tmp.begin());
    res_it += tmp.size();
    if ( nx_>0) f_ode.push_back(tmp[DAE_ODE]);
    if ( nz_>0) f_alg.push_back(tmp[DAE_ALG]);
    if ( nq_>0) f_quad.push_back(tmp[DAE_QUAD]);

    // Consistency check
    casadi_assert(res_it==res.end());

    vector<MatType> g_arg;
    if (!g_.is_null()) {

      // Forward derivatives of g
      d = g_.derivative(0, 0);
      g_arg.reserve(d.n_in());
      tmp.resize(RDAE_NUM_IN);
      fill(tmp.begin(), tmp.end(), MatType());

      // Reset iterators
      aug_x_split_it = aug_x_split.begin();
      aug_z_split_it = aug_z_split.begin();
      aug_p_split_it = aug_p_split.begin();

      // Collect arguments for calling d
      tmp[RDAE_T] = aug_t;
      if ( nx_>0) tmp[RDAE_X] = *aug_x_split_it++;
      if ( nz_>0) tmp[RDAE_Z] = *aug_z_split_it++;
      if ( np_>0) tmp[RDAE_P] = *aug_p_split_it++;
      if (nrx_>0) tmp[RDAE_RX] = *aug_rx_split_it++;
      if (nrz_>0) tmp[RDAE_RZ] = *aug_rz_split_it++;
      if (nrp_>0) tmp[RDAE_RP] = *aug_rp_split_it++;
      g_arg.insert(g_arg.end(), tmp.begin(), tmp.end());

      // Call d
      res = d(g_arg);
      res_it = res.begin();

      // Collect right-hand-sides
      tmp.resize(RDAE_NUM_OUT);
      fill(tmp.begin(), tmp.end(), MatType());
      copy(res_it, res_it+tmp.size(), tmp.begin());
      res_it += tmp.size();
      if (nrx_>0) g_ode.push_back(tmp[RDAE_ODE]);
      if (nrz_>0) g_alg.push_back(tmp[RDAE_ALG]);
      if (nrq_>0) g_quad.push_back(tmp[RDAE_QUAD]);

      // Consistency check
      casadi_assert(res_it==res.end());
    }

    // Adjoint derivatives of f
    d = f_.derivative(0, nadj);
    f_arg.resize(DAE_NUM_IN);
    f_arg.reserve(d.n_in());

    // Collect arguments for calling d
    tmp.resize(DAE_NUM_OUT);
    fill(tmp.begin(), tmp.end(), MatType());
    for (int dir=0; dir<nadj; ++dir) {
      if ( nx_>0) tmp[DAE_ODE] = *aug_rx_split_it++;
      if ( nz_>0) tmp[DAE_ALG] = *aug_rz_split_it++;
      if ( nq_>0) tmp[DAE_QUAD] = *aug_rp_split_it++;
      f_arg.insert(f_arg.end(), tmp.begin(), tmp.end());
    }

    // Call der
    res = d(f_arg);
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

    if (!g_.is_null()) {

      // Adjoint derivatives of g
      d = g_.derivative(0, nadj);
      g_arg.resize(RDAE_NUM_IN);
      g_arg.reserve(d.n_in());

      // Collect arguments for calling der
      tmp.resize(RDAE_NUM_OUT);
      fill(tmp.begin(), tmp.end(), MatType());
      for (int dir=0; dir<nadj; ++dir) {
        if (nrx_>0) tmp[RDAE_ODE] = *aug_x_split_it++;
        if (nrz_>0) tmp[RDAE_ALG] = *aug_z_split_it++;
        if (nrq_>0) tmp[RDAE_QUAD] = *aug_p_split_it++;
        g_arg.insert(g_arg.end(), tmp.begin(), tmp.end());
      }

      // Call der
      res = d(g_arg);
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
      if (nrx_>0) g_arg[RDAE_RX] = MatType::zeros(g_arg[RDAE_RX].sparsity());
      if (nrz_>0) g_arg[RDAE_RZ] = MatType::zeros(g_arg[RDAE_RZ].sparsity());
      if (nrp_>0) g_arg[RDAE_RP] = MatType::zeros(g_arg[RDAE_RP].sparsity());

      // Call der again
      res = d(g_arg);
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

    // Form the augmented forward integration
    ret["t"] = aug_t;
    ret["x"] = aug_x;
    ret["z"] = aug_z;
    ret["p"] = aug_p;
    if (!f_ode.empty()) ret["ode"] = densify(horzcat(f_ode));
    if (!f_alg.empty()) ret["alg"] = densify(horzcat(f_alg));
    if (!f_quad.empty()) ret["quad"] = densify(horzcat(f_quad));

    // Form the augmented backward integration function
    if (!g_ode.empty()) {
      ret["rx"] = aug_rx;
      ret["rz"] = aug_rz;
      ret["rp"] = aug_rp;
      if (!g_ode.empty()) ret["rode"] = densify(horzcat(g_ode));
      if (!g_alg.empty()) ret["ralg"] = densify(horzcat(g_alg));
      if (!g_quad.empty()) ret["rquad"] = densify(horzcat(g_quad));
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

  void Integrator::spFwd(const bvec_t** arg, bvec_t** res, int* iw, bvec_t* w, int mem) {
    log("Integrator::spFwd", "begin");

    // Work vectors
    bvec_t *tmp_x = w; w += nx_;
    bvec_t *tmp_z = w; w += nz_;
    bvec_t *tmp_rx = w; w += nrx_;
    bvec_t *tmp_rz = w; w += nrz_;

    // Propagate through f
    const bvec_t** arg1 = arg+n_in();
    fill(arg1, arg1+DAE_NUM_IN, static_cast<bvec_t*>(0));
    arg1[DAE_X] = arg[INTEGRATOR_X0];
    arg1[DAE_P] = arg[INTEGRATOR_P];
    bvec_t** res1 = res+n_out();
    fill(res1, res1+DAE_NUM_OUT, static_cast<bvec_t*>(0));
    res1[DAE_ODE] = tmp_x;
    res1[DAE_ALG] = tmp_z;
    f_(arg1, res1, iw, w, 0);
    if (arg[INTEGRATOR_X0]) {
      const bvec_t *tmp = arg[INTEGRATOR_X0];
      for (int i=0; i<nx_; ++i) tmp_x[i] |= *tmp++;
    }

    // "Solve" in order to resolve interdependencies (cf. Rootfinder)
    copy(tmp_x, tmp_x+nx_+nz_, w);
    fill_n(tmp_x, nx_+nz_, 0);
    casadi_assert(!linsol_f_.is_null());
    linsol_f_.linsol_spsolve(tmp_x, w, false);

    // Get xf and zf
    if (res[INTEGRATOR_XF])
      copy(tmp_x, tmp_x+nx_, res[INTEGRATOR_XF]);
    if (res[INTEGRATOR_ZF])
      copy(tmp_z, tmp_z+nz_, res[INTEGRATOR_ZF]);

    // Propagate to quadratures
    if (nq_>0 && res[INTEGRATOR_QF]) {
      arg1[DAE_X] = tmp_x;
      arg1[DAE_Z] = tmp_z;
      res1[DAE_ODE] = res1[DAE_ALG] = 0;
      res1[DAE_QUAD] = res[INTEGRATOR_QF];
      f_(arg1, res1, iw, w, 0);
    }

    if (!g_.is_null()) {
      // Propagate through g
      fill(arg1, arg1+RDAE_NUM_IN, static_cast<bvec_t*>(0));
      arg1[RDAE_X] = tmp_x;
      arg1[RDAE_P] = arg[INTEGRATOR_P];
      arg1[RDAE_Z] = tmp_z;
      arg1[RDAE_RX] = arg[INTEGRATOR_X0];
      arg1[RDAE_RX] = arg[INTEGRATOR_RX0];
      arg1[RDAE_RP] = arg[INTEGRATOR_RP];
      fill(res1, res1+RDAE_NUM_OUT, static_cast<bvec_t*>(0));
      res1[RDAE_ODE] = tmp_rx;
      res1[RDAE_ALG] = tmp_rz;
      g_(arg1, res1, iw, w, 0);
      if (arg[INTEGRATOR_RX0]) {
        const bvec_t *tmp = arg[INTEGRATOR_RX0];
        for (int i=0; i<nrx_; ++i) tmp_rx[i] |= *tmp++;
      }

      // "Solve" in order to resolve interdependencies (cf. Rootfinder)
      copy(tmp_rx, tmp_rx+nrx_+nrz_, w);
      fill_n(tmp_rx, nrx_+nrz_, 0);
      casadi_assert(!linsol_g_.is_null());
      linsol_g_.linsol_spsolve(tmp_rx, w, false);

      // Get rxf and rzf
      if (res[INTEGRATOR_RXF])
        copy(tmp_rx, tmp_rx+nrx_, res[INTEGRATOR_RXF]);
      if (res[INTEGRATOR_RZF])
        copy(tmp_rz, tmp_rz+nrz_, res[INTEGRATOR_RZF]);

      // Propagate to quadratures
      if (nrq_>0 && res[INTEGRATOR_RQF]) {
        arg1[RDAE_RX] = tmp_rx;
        arg1[RDAE_RZ] = tmp_rz;
        res1[RDAE_ODE] = res1[RDAE_ALG] = 0;
        res1[RDAE_QUAD] = res[INTEGRATOR_RQF];
        g_(arg1, res1, iw, w, 0);
      }
    }
    log("Integrator::spFwd", "end");
  }

  void Integrator::spAdj(bvec_t** arg, bvec_t** res, int* iw, bvec_t* w, int mem) {
    log("Integrator::spAdj", "begin");

    // Work vectors
    bvec_t** arg1 = arg+n_in();
    bvec_t** res1 = res+n_out();
    bvec_t *tmp_x = w; w += nx_;
    bvec_t *tmp_z = w; w += nz_;

    // Shorthands
    bvec_t* x0 = arg[INTEGRATOR_X0];
    bvec_t* p = arg[INTEGRATOR_P];
    bvec_t* xf = res[INTEGRATOR_XF];
    bvec_t* zf = res[INTEGRATOR_ZF];
    bvec_t* qf = res[INTEGRATOR_QF];

    // Propagate from outputs to state vectors
    if (xf) {
      copy_n(xf, nx_, tmp_x);
      fill_n(xf, nx_, 0);
    } else {
      fill_n(tmp_x, nx_, 0);
    }
    if (zf) {
      copy_n(zf, nz_, tmp_z);
      fill_n(zf, nz_, 0);
    } else {
      fill_n(tmp_z, nz_, 0);
    }

    if (!g_.is_null()) {
      // Work vectors
      bvec_t *tmp_rx = w; w += nrx_;
      bvec_t *tmp_rz = w; w += nrz_;

      // Shorthands
      bvec_t* rx0 = arg[INTEGRATOR_RX0];
      bvec_t* rp = arg[INTEGRATOR_RP];
      bvec_t* rxf = res[INTEGRATOR_RXF];
      bvec_t* rzf = res[INTEGRATOR_RZF];
      bvec_t* rqf = res[INTEGRATOR_RQF];

      // Propagate from outputs to state vectors
      if (rxf) {
        copy_n(rxf, nrx_, tmp_rx);
        fill_n(rxf, nrx_, 0);
      } else {
        fill_n(tmp_rx, nrx_, 0);
      }
      if (rzf) {
        copy_n(rzf, nrz_, tmp_rz);
        fill_n(rzf, nrz_, 0);
      } else {
        fill_n(tmp_rz, nrz_, 0);
      }

      // Get dependencies from backward quadratures
      fill(res1, res1+RDAE_NUM_OUT, static_cast<bvec_t*>(0));
      fill(arg1, arg1+RDAE_NUM_IN, static_cast<bvec_t*>(0));
      res1[RDAE_QUAD] = rqf;
      arg1[RDAE_X] = tmp_x;
      arg1[RDAE_Z] = tmp_z;
      arg1[RDAE_P] = p;
      arg1[RDAE_RX] = tmp_rx;
      arg1[RDAE_RZ] = tmp_rz;
      arg1[RDAE_RP] = rp;
      g_.rev(arg1, res1, iw, w, 0);

      // Propagate interdependencies
      casadi_assert(!linsol_g_.is_null());
      fill_n(w, nrx_+nrz_, 0);
      linsol_g_.linsol_spsolve(w, tmp_rx, true);
      copy(w, w+nrx_+nrz_, tmp_rx);

      // Direct dependency rx0 -> rxf
      if (rx0) for (int i=0; i<nrx_; ++i) rx0[i] |= tmp_rx[i];

      // Indirect dependency via g
      res1[RDAE_ODE] = tmp_rx;
      res1[RDAE_ALG] = tmp_rz;
      res1[RDAE_QUAD] = 0;
      arg1[RDAE_RX] = rx0;
      arg1[RDAE_RZ] = 0; // arg[INTEGRATOR_RZ0] is a guess, no dependency
      g_.rev(arg1, res1, iw, w, 0);
    }

    // Get dependencies from forward quadratures
    fill(res1, res1+DAE_NUM_OUT, static_cast<bvec_t*>(0));
    fill(arg1, arg1+DAE_NUM_IN, static_cast<bvec_t*>(0));
    res1[DAE_QUAD] = qf;
    arg1[DAE_X] = tmp_x;
    arg1[DAE_Z] = tmp_z;
    arg1[DAE_P] = p;
    if (qf && nq_>0) f_.rev(arg1, res1, iw, w, 0);

    // Propagate interdependencies
    casadi_assert(!linsol_f_.is_null());
    fill_n(w, nx_+nz_, 0);
    linsol_f_.linsol_spsolve(w, tmp_x, true);
    copy(w, w+nx_+nz_, tmp_x);

    // Direct dependency x0 -> xf
    if (x0) for (int i=0; i<nx_; ++i) x0[i] |= tmp_x[i];

    // Indirect dependency through f
    res1[DAE_ODE] = tmp_x;
    res1[DAE_ALG] = tmp_z;
    res1[DAE_QUAD] = 0;
    arg1[DAE_X] = x0;
    arg1[DAE_Z] = 0; // arg[INTEGRATOR_Z0] is a guess, no dependency
    f_.rev(arg1, res1, iw, w, 0);

    log("Integrator::spAdj", "end");
  }

  Integrator::AugOffset Integrator::getAugOffset(int nfwd, int nadj) {
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
      if ( nx_>0) ret.x.push_back(x().size2());
      if ( nz_>0) ret.z.push_back(z().size2());
      if ( nq_>0) ret.q.push_back(q().size2());
      if ( np_>0) ret.p.push_back(p().size2());
      if (nrx_>0) ret.rx.push_back(rx().size2());
      if (nrz_>0) ret.rz.push_back(rz().size2());
      if (nrq_>0) ret.rq.push_back(rq().size2());
      if (nrp_>0) ret.rp.push_back(rp().size2());
    }

    // Count adjoint sensitivities
    for (int dir=0; dir<nadj; ++dir) {
      if ( nx_>0) ret.rx.push_back(x().size2());
      if ( nz_>0) ret.rz.push_back(z().size2());
      if ( np_>0) ret.rq.push_back(p().size2());
      if ( nq_>0) ret.rp.push_back(q().size2());
      if (nrx_>0) ret.x.push_back(x().size2());
      if (nrz_>0) ret.z.push_back(rz().size2());
      if (nrp_>0) ret.q.push_back(rp().size2());
      if (nrq_>0) ret.p.push_back(rq().size2());
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

  Function Integrator::get_forward(const std::string& name, int nfwd, Dict& opts) {
    log("Integrator::get_forward", "begin");

    // Integrator options
    Dict aug_opts = getDerivativeOptions(true);
    for (auto&& i : augmented_options_) {
      aug_opts[i.first] = i.second;
    }

    // Temp stringstream
    stringstream ss;
    ss << "aug_f" << nfwd << name_;

    // Create integrator for augmented DAE
    Function aug_int;
    AugOffset offset;
    if (f_.is_a("sxfunction")) {
      SXDict aug_dae = aug_fwd<SX>(nfwd, offset);
      aug_int = integrator(ss.str(), plugin_name(), aug_dae, aug_opts);
    } else {
      casadi_assert(f_.is_a("mxfunction"));
      MXDict aug_dae = aug_fwd<MX>(nfwd, offset);
      aug_int = integrator(ss.str(), plugin_name(), aug_dae, aug_opts);
    }

    // All inputs of the return function
    vector<MX> ret_in;
    ret_in.reserve(INTEGRATOR_NUM_IN*(1+nfwd) + INTEGRATOR_NUM_OUT);

    // Augmented state
    vector<MX> x0_augv, p_augv, z0_augv, rx0_augv, rp_augv, rz0_augv;

    // Inputs or forward/adjoint seeds in one direction
    vector<MX> dd;

    // Add nondifferentiated inputs and forward seeds
    dd.resize(INTEGRATOR_NUM_IN);
    for (int dir=-1; dir<nfwd; ++dir) {

      // Differential state
      ss.clear();
      ss << "x0";
      if (dir>=0) ss << "_" << dir;
      dd[INTEGRATOR_X0] = MX::sym(ss.str(), x());
      x0_augv.push_back(dd[INTEGRATOR_X0]);

      // Parameter
      ss.clear();
      ss << "p";
      if (dir>=0) ss << "_" << dir;
      dd[INTEGRATOR_P] = MX::sym(ss.str(), p());
      p_augv.push_back(dd[INTEGRATOR_P]);

      // Initial guess for algebraic variable
      ss.clear();
      ss << "r0";
      if (dir>=0) ss << "_" << dir;
      dd[INTEGRATOR_Z0] = MX::sym(ss.str(), z());
      z0_augv.push_back(dd[INTEGRATOR_Z0]);

      // Backward state
      ss.clear();
      ss << "rx0";
      if (dir>=0) ss << "_" << dir;
      dd[INTEGRATOR_RX0] = MX::sym(ss.str(), rx());
      rx0_augv.push_back(dd[INTEGRATOR_RX0]);

      // Backward parameter
      ss.clear();
      ss << "rp";
      if (dir>=0) ss << "_" << dir;
      dd[INTEGRATOR_RP] = MX::sym(ss.str(), rp());
      rp_augv.push_back(dd[INTEGRATOR_RP]);

      // Initial guess for backward algebraic variable
      ss.clear();
      ss << "rz0";
      if (dir>=0) ss << "_" << dir;
      dd[INTEGRATOR_RZ0] = MX::sym(ss.str(), rz());
      rz0_augv.push_back(dd[INTEGRATOR_RZ0]);

      // Add to input vector
      ret_in.insert(ret_in.end(), dd.begin(), dd.end());

      // Make space for dummy outputs
      if (dir==-1) ret_in.resize(ret_in.size() + INTEGRATOR_NUM_OUT);
    }

    // Call the integrator
    vector<MX> integrator_in(INTEGRATOR_NUM_IN);
    integrator_in[INTEGRATOR_X0] = horzcat(x0_augv);
    integrator_in[INTEGRATOR_P] = horzcat(p_augv);
    integrator_in[INTEGRATOR_Z0] = horzcat(z0_augv);
    integrator_in[INTEGRATOR_RX0] = horzcat(rx0_augv);
    integrator_in[INTEGRATOR_RP] = horzcat(rp_augv);
    integrator_in[INTEGRATOR_RZ0] = horzcat(rz0_augv);
    vector<MX> integrator_out = aug_int(integrator_in);

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

    // Add dummy inputs (outputs of the nondifferentiated funciton)
    dd.resize(INTEGRATOR_NUM_OUT);
    dd[INTEGRATOR_XF]  = MX::sym("xf_dummy", Sparsity(size_out(INTEGRATOR_XF)));
    dd[INTEGRATOR_QF]  = MX::sym("qf_dummy", Sparsity(q().size()));
    dd[INTEGRATOR_ZF]  = MX::sym("zf_dummy", Sparsity(z().size()));
    dd[INTEGRATOR_RXF]  = MX::sym("rxf_dummy", Sparsity(rx().size()));
    dd[INTEGRATOR_RQF]  = MX::sym("rqf_dummy", Sparsity(rq().size()));
    dd[INTEGRATOR_RZF]  = MX::sym("rzf_dummy", Sparsity(rz().size()));
    std::copy(dd.begin(), dd.end(), ret_in.begin()+INTEGRATOR_NUM_IN);

    // All outputs of the return function
    vector<MX> ret_out;
    ret_out.reserve(INTEGRATOR_NUM_OUT*nfwd);

    // Collect the forward sensitivities
    fill(dd.begin(), dd.end(), MX());
    for (int dir=-1; dir<nfwd; ++dir) {
      if ( nx_>0) dd[INTEGRATOR_XF]  = *xf_aug_it++;
      if ( nq_>0) dd[INTEGRATOR_QF]  = *qf_aug_it++;
      if ( nz_>0) dd[INTEGRATOR_ZF]  = *zf_aug_it++;
      if (nrx_>0) dd[INTEGRATOR_RXF] = *rxf_aug_it++;
      if (nrq_>0) dd[INTEGRATOR_RQF] = *rqf_aug_it++;
      if (nrz_>0) dd[INTEGRATOR_RZF] = *rzf_aug_it++;
      if (dir>=0) // Nondifferentiated output ignored
        ret_out.insert(ret_out.end(), dd.begin(), dd.end());
    }
    log("Integrator::get_forward", "end");

    // Create derivative function and return
    return Function(name, ret_in, ret_out, opts);
  }

  Function Integrator::get_reverse(const std::string& name, int nadj, Dict& opts) {
    log("Integrator::get_reverse", "begin");

    // Integrator options
    Dict aug_opts = getDerivativeOptions(false);
    for (auto&& i : augmented_options_) {
      aug_opts[i.first] = i.second;
    }

    // Temp stringstream
    stringstream ss;
    ss << "aug_r" << nadj << name_;

    // Create integrator for augmented DAE
    Function aug_int;
    AugOffset offset;
    if (f_.is_a("sxfunction")) {
      SXDict aug_dae = aug_adj<SX>(nadj, offset);
      aug_int = integrator(ss.str(), plugin_name(), aug_dae, aug_opts);
    } else {
      casadi_assert(f_.is_a("mxfunction"));
      MXDict aug_dae = aug_adj<MX>(nadj, offset);
      aug_int = integrator(ss.str(), plugin_name(), aug_dae, aug_opts);
    }

    // All inputs of the return function
    vector<MX> ret_in;
    ret_in.reserve(INTEGRATOR_NUM_IN + INTEGRATOR_NUM_OUT*(1+nadj));

    // Augmented state
    vector<MX> x0_augv, p_augv, z0_augv, rx0_augv, rp_augv, rz0_augv;

    // Inputs or forward/adjoint seeds in one direction
    vector<MX> dd;

    // Add nondifferentiated inputs and forward seeds
    dd.resize(INTEGRATOR_NUM_IN);
    fill(dd.begin(), dd.end(), MX());

    // Differential state
    dd[INTEGRATOR_X0] = MX::sym("x0", x());
    x0_augv.push_back(dd[INTEGRATOR_X0]);

    // Parameter
    dd[INTEGRATOR_P] = MX::sym("p", p());
    p_augv.push_back(dd[INTEGRATOR_P]);

    // Initial guess for algebraic variable
    dd[INTEGRATOR_Z0] = MX::sym("r0", z());
    z0_augv.push_back(dd[INTEGRATOR_Z0]);

    // Backward state
    dd[INTEGRATOR_RX0] = MX::sym("rx0", rx());
    rx0_augv.push_back(dd[INTEGRATOR_RX0]);

    // Backward parameter
    dd[INTEGRATOR_RP] = MX::sym("rp", rp());
    rp_augv.push_back(dd[INTEGRATOR_RP]);

    // Initial guess for backward algebraic variable
    dd[INTEGRATOR_RZ0] = MX::sym("rz0", rz());
    rz0_augv.push_back(dd[INTEGRATOR_RZ0]);

    // Add to input vector
    ret_in.insert(ret_in.end(), dd.begin(), dd.end());

    // Add dummy inputs (outputs of the nondifferentiated funciton)
    dd.resize(INTEGRATOR_NUM_OUT);
    dd[INTEGRATOR_XF]  = MX::sym("xf_dummy", Sparsity(x().size()));
    dd[INTEGRATOR_QF]  = MX::sym("qf_dummy", Sparsity(q().size()));
    dd[INTEGRATOR_ZF]  = MX::sym("zf_dummy", Sparsity(z().size()));
    dd[INTEGRATOR_RXF]  = MX::sym("rxf_dummy", Sparsity(rx().size()));
    dd[INTEGRATOR_RQF]  = MX::sym("rqf_dummy", Sparsity(rq().size()));
    dd[INTEGRATOR_RZF]  = MX::sym("rzf_dummy", Sparsity(rz().size()));
    ret_in.insert(ret_in.end(), dd.begin(), dd.end());

    // Add adjoint seeds
    dd.resize(INTEGRATOR_NUM_OUT);
    fill(dd.begin(), dd.end(), MX());
    for (int dir=0; dir<nadj; ++dir) {

      // Differential states become backward differential state
      ss.clear();
      ss << "xf" << "_" << dir;
      dd[INTEGRATOR_XF] = MX::sym(ss.str(), x());
      rx0_augv.push_back(dd[INTEGRATOR_XF]);

      // Quadratures become backward parameters
      ss.clear();
      ss << "qf" << "_" << dir;
      dd[INTEGRATOR_QF] = MX::sym(ss.str(), q());
      rp_augv.push_back(dd[INTEGRATOR_QF]);

      // Algebraic variables become backward algebraic variables
      ss.clear();
      ss << "zf" << "_" << dir;
      dd[INTEGRATOR_ZF] = MX::sym(ss.str(), z());
      rz0_augv.push_back(dd[INTEGRATOR_ZF]);

      // Backward differential states becomes forward differential states
      ss.clear();
      ss << "rxf" << "_" << dir;
      dd[INTEGRATOR_RXF] = MX::sym(ss.str(), rx());
      x0_augv.push_back(dd[INTEGRATOR_RXF]);

      // Backward quadratures becomes (forward) parameters
      ss.clear();
      ss << "rqf" << "_" << dir;
      dd[INTEGRATOR_RQF] = MX::sym(ss.str(), rq());
      p_augv.push_back(dd[INTEGRATOR_RQF]);

      // Backward differential states becomes forward differential states
      ss.clear();
      ss << "rzf" << "_" << dir;
      dd[INTEGRATOR_RZF] = MX::sym(ss.str(), rz());
      z0_augv.push_back(dd[INTEGRATOR_RZF]);

      // Add to input vector
      ret_in.insert(ret_in.end(), dd.begin(), dd.end());
    }

    // Call the integrator
    vector<MX> integrator_in(INTEGRATOR_NUM_IN);
    integrator_in[INTEGRATOR_X0] = horzcat(x0_augv);
    integrator_in[INTEGRATOR_P] = horzcat(p_augv);
    integrator_in[INTEGRATOR_Z0] = horzcat(z0_augv);
    integrator_in[INTEGRATOR_RX0] = horzcat(rx0_augv);
    integrator_in[INTEGRATOR_RP] = horzcat(rp_augv);
    integrator_in[INTEGRATOR_RZ0] = horzcat(rz0_augv);
    vector<MX> integrator_out = aug_int(integrator_in);

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
    ret_out.reserve(INTEGRATOR_NUM_IN*nadj);

    // Collect the nondifferentiated results and forward sensitivities
    dd.resize(INTEGRATOR_NUM_OUT);
    fill(dd.begin(), dd.end(), MX());
    for (int dir=-1; dir<0; ++dir) {
      if ( nx_>0) dd[INTEGRATOR_XF]  = *xf_aug_it++;
      if ( nq_>0) dd[INTEGRATOR_QF]  = *qf_aug_it++;
      if ( nz_>0) dd[INTEGRATOR_ZF]  = *zf_aug_it++;
      if (nrx_>0) dd[INTEGRATOR_RXF] = *rxf_aug_it++;
      if (nrq_>0) dd[INTEGRATOR_RQF] = *rqf_aug_it++;
      if (nrz_>0) dd[INTEGRATOR_RZF] = *rzf_aug_it++;
      //ret_out.insert(ret_out.end(), dd.begin(), dd.end());
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
    log("Integrator::getDerivative", "end");

    // Create derivative function and return
    return Function(name, ret_in, ret_out, opts);
  }

  void Integrator::set_temp(Memory& mem, const double** arg, double** res,
                            int* iw, double* w) const {
    IntegratorMemory& m = dynamic_cast<IntegratorMemory&>(mem);
    m.arg = arg;
    m.res = res;
    m.iw = iw;
    m.w = w;
  }

  Dict Integrator::getDerivativeOptions(bool fwd) {
    // Copy all options
    return opts_;
  }

  Sparsity Integrator::spJacF() {
    // Start with the sparsity pattern of the ODE part
    Sparsity jac_ode_x = f_.sparsity_jac(DAE_X, DAE_ODE);

    // Add diagonal to get interdependencies
    jac_ode_x = jac_ode_x + Sparsity::diag(nx_);

    // Quick return if no algebraic variables
    if (nz_==0) return jac_ode_x;

    // Add contribution from algebraic variables and equations
    Sparsity jac_ode_z = f_.sparsity_jac(DAE_Z, DAE_ODE);
    Sparsity jac_alg_x = f_.sparsity_jac(DAE_X, DAE_ALG);
    Sparsity jac_alg_z = f_.sparsity_jac(DAE_Z, DAE_ALG);
    return blockcat(jac_ode_x, jac_ode_z,
                    jac_alg_x, jac_alg_z);
  }

  Sparsity Integrator::spJacG() {
    // Start with the sparsity pattern of the ODE part
    Sparsity jac_ode_x = g_.sparsity_jac(RDAE_RX, RDAE_ODE);

    // Add diagonal to get interdependencies
    jac_ode_x = jac_ode_x + Sparsity::diag(nrx_);

    // Quick return if no algebraic variables
    if (nrz_==0) return jac_ode_x;

    // Add contribution from algebraic variables and equations
    Sparsity jac_ode_z = g_.sparsity_jac(RDAE_RZ, RDAE_ODE);
    Sparsity jac_alg_x = g_.sparsity_jac(RDAE_RX, RDAE_ALG);
    Sparsity jac_alg_z = g_.sparsity_jac(RDAE_RZ, RDAE_ALG);
    return blockcat(jac_ode_x, jac_ode_z,
                    jac_alg_x, jac_alg_z);
  }

  std::map<std::string, Integrator::Plugin> Integrator::solvers_;

  const std::string Integrator::infix_ = "integrator";

  void Integrator::setStopTime(IntegratorMemory& mem, double tf) const {
    casadi_error("Integrator::setStopTime not defined for class "
                 << typeid(*this).name());
  }

  FixedStepIntegrator::FixedStepIntegrator(const std::string& name, const XProblem& dae)
    : Integrator(name, dae) {

    // Default options
    nk_ = 20;
  }

  FixedStepIntegrator::~FixedStepIntegrator() {
  }

  Options FixedStepIntegrator::options_
  = {{&Integrator::options_},
     {{"number_of_finite_elements",
       {OT_INT,
        "Number of finite elements"}}
     }
  };

  void FixedStepIntegrator::init(const Dict& opts) {
    // Call the base class init
    Integrator::init(opts);

    // Read options
    for (auto&& op : opts) {
      if (op.first=="number_of_finite_elements") {
        nk_ = op.second;
      }
    }

    // Number of finite elements and time steps
    casadi_assert(nk_>0);
    h_ = (grid_.back() - grid_.front())/nk_;

    // Setup discrete time dynamics
    setupFG();

    // Get discrete time dimensions
    nZ_ = F_.nnz_in(DAE_Z);
    nRZ_ =  G_.is_null() ? 0 : G_.nnz_in(RDAE_RZ);
  }

  void FixedStepIntegrator::init_memory(Memory& mem) const {
    Integrator::init_memory(mem);
    FixedStepMemory& m = dynamic_cast<FixedStepMemory&>(mem);

    // Discrete time algebraic variable
    m.Z = DM::zeros(F_.sparsity_in(DAE_Z));
    m.RZ = G_.is_null() ? DM() : DM::zeros(G_.sparsity_in(RDAE_RZ));

    // Allocate tape if backward states are present
    if (nrx_>0) {
      m.x_tape.resize(nk_+1, vector<double>(nx_));
      m.Z_tape.resize(nk_, vector<double>(nZ_));
    }

    // Allocate state
    m.x.resize(nx_);
    m.z.resize(nz_);
    m.p.resize(np_);
    m.q.resize(nq_);
    m.rx.resize(nrx_);
    m.rz.resize(nrz_);
    m.rp.resize(nrp_);
    m.rq.resize(nrq_);
    m.x_prev.resize(nx_);
    m.Z_prev.resize(nZ_);
    m.q_prev.resize(nq_);
    m.rx_prev.resize(nrx_);
    m.RZ_prev.resize(nRZ_);
    m.rq_prev.resize(nrq_);
  }

  void FixedStepIntegrator::advance(IntegratorMemory& mem, double t,
                                    double* x, double* z, double* q) const {
    FixedStepMemory& m = dynamic_cast<FixedStepMemory&>(mem);

    // Get discrete time sought
    int k_out = std::ceil((t - grid_.front())/h_);
    k_out = std::min(k_out, nk_); //  make sure that rounding errors does not result in k_out>nk_
    casadi_assert(k_out>=0);

    // Explicit discrete time dynamics
    const Function& F = getExplicit();

    // Discrete dynamics function inputs ...
    fill_n(m.arg, F.n_in(), nullptr);
    m.arg[DAE_T] = &m.t;
    m.arg[DAE_X] = get_ptr(m.x_prev);
    m.arg[DAE_Z] = get_ptr(m.Z_prev);
    m.arg[DAE_P] = get_ptr(m.p);

    // ... and outputs
    fill_n(m.res, F.n_out(), nullptr);
    m.res[DAE_ODE] = get_ptr(m.x);
    m.res[DAE_ALG] = get_ptr(m.Z);
    m.res[DAE_QUAD] = get_ptr(m.q);

    // Take time steps until end time has been reached
    while (m.k<k_out) {
      // Update the previous step
      casadi_copy(get_ptr(m.x), nx_, get_ptr(m.x_prev));
      casadi_copy(get_ptr(m.Z), nZ_, get_ptr(m.Z_prev));
      casadi_copy(get_ptr(m.q), nq_, get_ptr(m.q_prev));

      // Take step
      F(m.arg, m.res, m.iw, m.w, 0);
      casadi_axpy(nq_, 1., get_ptr(m.q_prev), get_ptr(m.q));

      // Tape
      if (nrx_>0) {
        casadi_copy(get_ptr(m.x), nx_, get_ptr(m.x_tape.at(m.k+1)));
        casadi_copy(get_ptr(m.Z), m.Z.nnz(), get_ptr(m.Z_tape.at(m.k)));
      }

      // Advance time
      m.k++;
      m.t = grid_.front() + m.k*h_;
    }

    // Return to user TODO(@jaeandersson): interpolate
    casadi_copy(get_ptr(m.x), nx_, x);
    casadi_copy(get_ptr(m.Z)+m.Z.nnz()-nz_, nz_, z);
    casadi_copy(get_ptr(m.q), nq_, q);
  }

  void FixedStepIntegrator::retreat(IntegratorMemory& mem, double t,
                                    double* rx, double* rz, double* rq) const {
    FixedStepMemory& m = dynamic_cast<FixedStepMemory&>(mem);

    // Get discrete time sought
    int k_out = std::floor((t - grid_.front())/h_);
    k_out = std::max(k_out, 0); //  make sure that rounding errors does not result in k_out>nk_
    casadi_assert(k_out<=nk_);

    // Explicit discrete time dynamics
    const Function& G = getExplicitB();

    // Discrete dynamics function inputs ...
    fill_n(m.arg, G.n_in(), nullptr);
    m.arg[RDAE_T] = &m.t;
    m.arg[RDAE_P] = get_ptr(m.p);
    m.arg[RDAE_RX] = get_ptr(m.rx_prev);
    m.arg[RDAE_RZ] = get_ptr(m.RZ_prev);
    m.arg[RDAE_RP] = get_ptr(m.rp);

    // ... and outputs
    fill_n(m.res, G.n_out(), nullptr);
    m.res[RDAE_ODE] = get_ptr(m.rx);
    m.res[RDAE_ALG] = get_ptr(m.RZ);
    m.res[RDAE_QUAD] = get_ptr(m.rq);

    // Take time steps until end time has been reached
    while (m.k>k_out) {
      // Advance time
      m.k--;
      m.t = grid_.front() + m.k*h_;

      // Update the previous step
      casadi_copy(get_ptr(m.rx), nrx_, get_ptr(m.rx_prev));
      casadi_copy(get_ptr(m.RZ), nRZ_, get_ptr(m.RZ_prev));
      casadi_copy(get_ptr(m.rq), nrq_, get_ptr(m.rq_prev));

      // Take step
      m.arg[RDAE_X] = get_ptr(m.x_tape.at(m.k));
      m.arg[RDAE_Z] = get_ptr(m.Z_tape.at(m.k));
      G(m.arg, m.res, m.iw, m.w, 0);
      casadi_axpy(nrq_, 1., get_ptr(m.rq_prev), get_ptr(m.rq));
    }

    // Return to user TODO(@jaeandersson): interpolate
    casadi_copy(get_ptr(m.rx), nrx_, rx);
    casadi_copy(get_ptr(m.RZ)+m.RZ.nnz()-nrz_, nrz_, rz);
    casadi_copy(get_ptr(m.rq), nrq_, rq);
  }

  void FixedStepIntegrator::
  reset(IntegratorMemory& mem, double t,
        const double* x, const double* z, const double* p) const {
    FixedStepMemory& m = dynamic_cast<FixedStepMemory&>(mem);

    // Update time
    m.t = t;

    // Set parameters
    casadi_copy(p, np_, get_ptr(m.p));

    // Update the state
    casadi_copy(x, nx_, get_ptr(m.x));
    casadi_copy(z, nz_, get_ptr(m.z));

    // Reset summation states
    casadi_fill(get_ptr(m.q), nq_, 0.);

    // Bring discrete time to the beginning
    m.k = 0;

    // Get consistent initial conditions
    casadi_fill(m.Z.ptr(), m.Z.nnz(), numeric_limits<double>::quiet_NaN());

    // Add the first element in the tape
    if (nrx_>0) {
      casadi_copy(x, nx_, get_ptr(m.x_tape.at(0)));
    }
  }

  void FixedStepIntegrator::resetB(IntegratorMemory& mem, double t, const double* rx,
                                   const double* rz, const double* rp) const {
    FixedStepMemory& m = dynamic_cast<FixedStepMemory&>(mem);

    // Update time
    m.t = t;

    // Set parameters
    casadi_copy(rp, nrp_, get_ptr(m.rp));

    // Update the state
    casadi_copy(rx, nrx_, get_ptr(m.rx));
    casadi_copy(rz, nrz_, get_ptr(m.rz));

    // Reset summation states
    casadi_fill(get_ptr(m.rq), nrq_, 0.);

    // Bring discrete time to the end
    m.k = nk_;

    // Get consistent initial conditions
    casadi_fill(m.RZ.ptr(), m.RZ.nnz(), numeric_limits<double>::quiet_NaN());
  }

  ImplicitFixedStepIntegrator::
  ImplicitFixedStepIntegrator(const std::string& name, const XProblem& dae)
    : FixedStepIntegrator(name, dae) {
  }

  ImplicitFixedStepIntegrator::~ImplicitFixedStepIntegrator() {
  }

  Options ImplicitFixedStepIntegrator::options_
  = {{&FixedStepIntegrator::options_},
     {{"implicit_solver",
       {OT_STRING,
        "An implicit function solver"}},
      {"implicit_solver_options",
       {OT_DICT,
        "Options to be passed to the NLP Solver"}}
     }
  };

  void ImplicitFixedStepIntegrator::init(const Dict& opts) {
    // Call the base class init
    FixedStepIntegrator::init(opts);

    // Default (temporary) options
    std::string implicit_function_name;
    Dict implicit_solver_options;

    // Read options
    for (auto&& op : opts) {
      if (op.first=="implicit_solver") {
        implicit_function_name = op.second.to_string();
      } else if (op.first=="implicit_solver_options") {
        implicit_solver_options = op.second;
      }
    }

    // Complete implicit_solver dictionary
    implicit_solver_options["implicit_input"] = DAE_Z;
    implicit_solver_options["implicit_output"] = DAE_ALG;

    // Allocate a solver
    implicit_solver_ = rootfinder(name_ + "_implicit_solver", implicit_function_name,
                                  F_, implicit_solver_options);
    alloc(implicit_solver_);

    // Allocate a root-finding solver for the backward problem
    if (nRZ_>0) {
      // Options
      Dict backward_implicit_solver_options = implicit_solver_options;
      backward_implicit_solver_options["implicit_input"] = RDAE_RZ;
      backward_implicit_solver_options["implicit_output"] = RDAE_ALG;
      string backward_implicit_function_name = implicit_function_name;

      // Allocate a Newton solver
      backward_implicit_solver_ =
        rootfinder(name_+ "_backward_implicit_solver",
                   backward_implicit_function_name,
                   G_, backward_implicit_solver_options);
      alloc(backward_implicit_solver_);
    }
  }

} // namespace casadi

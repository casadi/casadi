/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A symbolic framework for dynamic optimization.
 *    Copyright (C) 2010-2026 Joel Andersson, Joris Gillis, Moritz Diehl,
 *                            KU Leuven. All rights reserved.
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


#include "rosenbrock.hpp"
#include "rosenbrock_tableau.hpp"

namespace casadi {

  extern "C"
  int CASADI_INTEGRATOR_ROSENBROCK_EXPORT
      casadi_register_integrator_rosenbrock(Integrator::Plugin* plugin) {
    plugin->creator = Rosenbrock::creator;
    plugin->name = "rosenbrock";
    plugin->doc = Rosenbrock::meta_doc.c_str();
    plugin->version = CASADI_VERSION;
    plugin->options = &Rosenbrock::options_;
    plugin->deserialize = &Rosenbrock::deserialize;
    return 0;
  }

  extern "C"
  void CASADI_INTEGRATOR_ROSENBROCK_EXPORT casadi_load_integrator_rosenbrock() {
    Integrator::registerPlugin(casadi_register_integrator_rosenbrock);
  }

  Rosenbrock::Rosenbrock(const std::string& name, const Function& dae, double t0,
      const std::vector<double>& tout) : Integrator(name, dae, t0, tout) {
  }

  Rosenbrock::~Rosenbrock() {
    clear_mem();
  }

  const Options Rosenbrock::options_
  = {{&Integrator::options_},
     {{"abstol",
       {OT_DOUBLE,
        "Absolute tolerance [1e-8]"}},
      {"reltol",
       {OT_DOUBLE,
        "Relative tolerance [1e-6]"}},
      {"max_num_steps",
       {OT_INT,
        "Maximum number of steps per output interval [10000]"}},
      {"step0",
       {OT_DOUBLE,
        "Initial step size [default: estimated]"}},
      {"max_step_size",
       {OT_DOUBLE,
        "Maximum step size [default: unbounded]"}},
      {"min_step_size",
       {OT_DOUBLE,
        "Minimum step size [1e-14]"}},
      {"fixed_steps",
       {OT_INT,
        "Number of equidistant steps per output interval, "
        "no error control [default: 0, adaptive]"}},
      {"quad_err_con",
       {OT_BOOL,
        "Include the quadratures in the error control [false]"}},
      {"init_tol",
       {OT_DOUBLE,
        "Tolerance on the scaled Newton step when making the algebraic states "
        "consistent [1e-3]"}},
      {"init_max_iter",
       {OT_INT,
        "Maximum number of Newton iterations for consistent algebraic states [20]"}},
      {"max_tape_steps",
       {OT_INT,
        "Capacity, in steps, of a tape in the work vector, for adjoint sensitivities "
        "with adaptive steps; needed for code generation [0: tape of unbounded size]. "
        "With fixed_steps, the tape is sized exactly."}}
     }
  };

  void Rosenbrock::init(const Dict& opts) {
    // Call the base class init
    Integrator::init(opts);

    // Default options
    abstol_ = 1e-8;
    reltol_ = 1e-6;
    max_num_steps_ = 10000;
    step0_ = 0;
    max_step_size_ = 0;
    min_step_size_ = 1e-14;
    fixed_steps_ = 0;
    quad_err_con_ = false;
    init_tol_ = 1e-3;
    init_max_iter_ = 20;
    max_tape_steps_ = 0;

    // Read options
    for (auto&& op : opts) {
      if (op.first == "abstol") {
        abstol_ = op.second;
      } else if (op.first == "reltol") {
        reltol_ = op.second;
      } else if (op.first == "max_num_steps") {
        max_num_steps_ = op.second;
      } else if (op.first == "step0") {
        step0_ = op.second;
      } else if (op.first == "max_step_size") {
        max_step_size_ = op.second;
      } else if (op.first == "min_step_size") {
        min_step_size_ = op.second;
      } else if (op.first == "fixed_steps") {
        fixed_steps_ = op.second;
      } else if (op.first == "quad_err_con") {
        quad_err_con_ = op.second;
      } else if (op.first == "init_tol") {
        init_tol_ = op.second;
      } else if (op.first == "init_max_iter") {
        init_max_iter_ = op.second;
      } else if (op.first == "max_tape_steps") {
        max_tape_steps_ = op.second;
      }
    }

    // Events not supported
    casadi_assert(ne_ == 0, "Zero-crossing events are not supported by 'rosenbrock'");

    // Error control on the nominal states only, so derivative integrators share the mesh
    if (derivative_of_.is_null()) {
      nx_err_ = nx1_;
      nz_err_ = nz1_;
      nq_err_ = nq1_;
      nu_err_ = nu1_;
    } else {
      const Rosenbrock* d = derivative_of_.get<Rosenbrock>();
      casadi_assert_dev(d != nullptr);
      nx_err_ = d->nx_err_;
      nz_err_ = d->nz_err_;
      nq_err_ = d->nq_err_;
      nu_err_ = d->nu_err_;
    }

    // Jacobian blocks w.r.t. x, z and time, in the order of casadi_rosenbrock_jac_t
    Function jac = create_function("jacF", {"t", "x", "z", "p", "u"},
      {"jac:ode:x", "jac:alg:x", "jac:quad:x", "jac:ode:z", "jac:alg:z", "jac:quad:z",
       "jac:ode:t", "jac:alg:t", "jac:quad:t"});
    sp_jac_.resize(ROSENBROCK_NUM_JAC);
    for (casadi_int k = 0; k < ROSENBROCK_NUM_JAC; ++k) sp_jac_[k] = jac.sparsity_out(k);

    // Right-hand side [ode; alg; quad] and its vector-Jacobian product, one function each
    create_function("rhsF", dyn_in(), {"ode", "alg", "quad"});
    if (nadj_ > 0) {
      create_function(rdae_, "vjpB", bdyn_in(), {"adj_x", "adj_z", "adj_p", "adj_u"});
      set_function(jac.reverse(nadj_), "jacB");
    }

    // Forward sensitivities of all of these
    if (nfwd_ > 0) {
      create_forward("rhsF", nfwd_);
      create_forward("jacF", nfwd_);
      if (nadj_ > 0) {
        create_forward("vjpB", nfwd_);
        create_forward("jacB", nfwd_);
      }
    }

    // Iteration matrix M/(h gamma) - J of the nominal [x; z]
    sp_iter_ = blockcat(sp_jac_[ROSENBROCK_ODE_X] + Sparsity::diag(nx1_),
      sp_jac_[ROSENBROCK_ODE_Z], sp_jac_[ROSENBROCK_ALG_X], sp_jac_[ROSENBROCK_ALG_Z]);
    // Positions of the nonzeros of the blocks ODE_X, ALG_X, ODE_Z, ALG_Z, then the x-diagonal
    iter_map_.clear();
    const casadi_int blk[4] = {ROSENBROCK_ODE_X, ROSENBROCK_ALG_X, ROSENBROCK_ODE_Z,
      ROSENBROCK_ALG_Z};
    for (casadi_int k = 0; k < 4; ++k) {
      std::vector<casadi_int> r = sp_jac_[blk[k]].get_row(), c = sp_jac_[blk[k]].get_col();
      casadi_int r_off = blk[k] == ROSENBROCK_ALG_X || blk[k] == ROSENBROCK_ALG_Z ? nx1_ : 0;
      casadi_int c_off = blk[k] == ROSENBROCK_ODE_Z || blk[k] == ROSENBROCK_ALG_Z ? nx1_ : 0;
      for (size_t i = 0; i < r.size(); ++i) {
        iter_map_.push_back(sp_iter_.get_nz(r[i] + r_off, c[i] + c_off));
      }
    }
    for (casadi_int i = 0; i < nx1_; ++i) iter_map_.push_back(sp_iter_.get_nz(i, i));

    // Symbolic sparse QR factorizations
    sp_iter_.qr_sparse(sp_iter_v_, sp_iter_r_, iter_prinv_, iter_pc_);
    if (nz1_ > 0) sp_jac_[ROSENBROCK_ALG_Z].qr_sparse(sp_alg_v_, sp_alg_r_, alg_prinv_, alg_pc_);

    // Reverse communication problem
    set_rosenbrock_prob();

    // Forward directions folded into the DAE change the arithmetic of the nominal states:
    // replay the steps of a nominal integration with the functions of the original integrator
    const Rosenbrock* r = root();
    if (r != this && r->nx1_ != nx1_ && fixed_steps_ == 0) {
      set_function(r->get_function("rhsF"), "rhsF_nom");
      set_function(r->get_function("jacF"), "jacF_nom");
    }
    set_nominal_prob();

    // Work vectors of the kernel, and of the nominal integration
    casadi_int sz_arg = 0, sz_res = 0, sz_iw = 0, sz_w = 0;
    casadi_rosenbrock_work(&p_, &sz_arg, &sz_res, &sz_iw, &sz_w);
    if (replay_) casadi_rosenbrock_work(&p_nom_, &sz_arg, &sz_res, &sz_iw, &sz_w);
    alloc_w(sz_w, true);
    // Parameters and controls in generated code
    alloc_w(np_ + nu_, true);
  }

  const Rosenbrock* Rosenbrock::root() const {
    // Integrator at the start of the chain of derivatives
    const Rosenbrock* r = this;
    while (!r->derivative_of_.is_null()) r = r->derivative_of_.get<Rosenbrock>();
    return r;
  }

  void Rosenbrock::set_nominal_prob() {
    // Problem of the original integrator, if its steps are replayed
    replay_ = has_function("rhsF_nom");
    calls_ = task_calls("");
    if (!replay_) return;
    calls_nom_ = root()->task_calls("_nom");
    p_nom_ = root()->p_;
    p_nom_.tape_size = 0;
    casadi_rosenbrock_setup(&p_nom_);
  }

  void Rosenbrock::run(RosenbrockMemory* m) const {
    // Carry out the tasks requested by the kernel, for it or for the nominal integration
    while (casadi_rosenbrock(&m->d)) {
      casadi_rosenbrock_data<double>* a = m->d.active;
      if (a == &m->d ? task(m) : eval_task(m, a, calls_nom_)) a->fail = 1;
    }
  }

  casadi_int Rosenbrock::tape_capacity() const {
    casadi_int ny = nx_ + nz_ + nq_, n_step = fixed_steps_ > 0 ? fixed_steps_ * nt()
      : max_tape_steps_, n_init;
    if (nadj_ == 0 || n_step == 0) return 0;
    // Every output interval may start with a consistent initialization
    n_init = nz_ > 0 ? nt() : 0;
    return n_step * (ny * (rosenbrock_tableau::n_stages + 1) + nu_ + 5) + n_init * (ny + nu_ + 5);
  }

  void Rosenbrock::set_rosenbrock_prob() {
    // Coefficients in original form (rosenbrock_tableau.hpp), Gamma lower triangular
    using namespace rosenbrock_tableau;
    const casadi_int s = n_stages;
    // Inverse of Gamma = (gamma_ij), lower triangular, in long double
    std::vector<long double> Gi(s * s, 0);
    for (casadi_int j = 0; j < s; ++j) {
      Gi[j * s + j] = 1 / gamma;
      for (casadi_int i = j + 1; i < s; ++i) {
        long double acc = 0;
        for (casadi_int k = j; k < i; ++k) acc += gamma_ij[i][k] * Gi[k * s + j];
        Gi[i * s + j] = -acc / gamma;
      }
    }
    // a = A Gamma^-1, c = diag(1/gamma) - Gamma^-1, m = Gamma^-T b, m_err = Gamma^-T (b - bhat)
    a_.assign(s * s, 0);
    c_.assign(s * s, 0);
    alpha_.assign(s, 0);
    gamma_sum_.assign(s, 0);
    m_.assign(s, 0);
    m_err_.assign(s, 0);
    for (casadi_int i = 0; i < s; ++i) {
      long double alpha_i = 0, gamma_i = gamma;
      for (casadi_int j = 0; j < i; ++j) {
        alpha_i += alpha[i][j];
        gamma_i += gamma_ij[i][j];
        long double aij = 0;
        for (casadi_int k = j; k < i; ++k) aij += alpha[i][k] * Gi[k * s + j];
        a_[i * s + j] = static_cast<double>(aij);
        c_[i * s + j] = static_cast<double>(-Gi[i * s + j]);
      }
      alpha_[i] = static_cast<double>(alpha_i);
      gamma_sum_[i] = static_cast<double>(gamma_i);
      long double mi = 0, mhi = 0;
      for (casadi_int k = i; k < s; ++k) {
        mi += b[k] * Gi[k * s + i];
        mhi += bhat[k] * Gi[k * s + i];
      }
      m_[i] = static_cast<double>(mi);
      m_err_[i] = static_cast<double>(mi - mhi);
    }
    // Dense output weights in the same form: m(theta) = Gamma^-T b(theta)
    m_dense_.assign(n_dense * s, 0);
    for (casadi_int q = 0; q < n_dense; ++q) {
      for (casadi_int i = 0; i < s; ++i) {
        long double mi = 0;
        for (casadi_int k = i; k < s; ++k) mi += b_dense[q][k] * Gi[k * s + i];
        m_dense_[q * s + i] = static_cast<double>(mi);
      }
    }

    // Problem structure for the kernel
    p_.nx1 = nx1_;
    p_.nz1 = nz1_;
    p_.nq1 = nq1_;
    p_.np1 = np1_;
    p_.nu1 = nu1_;
    p_.nfwd = nfwd_;
    p_.nadj = nadj_;
    p_.nx_err = nx_err_;
    p_.nz_err = nz_err_;
    p_.nq_err = quad_err_con_ ? nq_err_ : 0;
    p_.nu_err = nu_err_;
    p_.n_stages = s;
    p_.order_err = order_embedded;
    p_.gamma = static_cast<double>(gamma);
    p_.a = get_ptr(a_);
    p_.c = get_ptr(c_);
    p_.alpha = get_ptr(alpha_);
    p_.gamma_sum = get_ptr(gamma_sum_);
    p_.m = get_ptr(m_);
    p_.m_err = get_ptr(m_err_);
    p_.n_dense = n_dense;
    p_.m_dense = get_ptr(m_dense_);
    for (casadi_int k = 0; k < ROSENBROCK_NUM_JAC; ++k) p_.sp_jac[k] = sp_jac_[k];
    p_.sp_iter = sp_iter_;
    p_.iter_map = get_ptr(iter_map_);
    p_.sp_iter_v = sp_iter_v_;
    p_.sp_iter_r = sp_iter_r_;
    p_.iter_prinv = get_ptr(iter_prinv_);
    p_.iter_pc = get_ptr(iter_pc_);
    p_.sp_alg_v = nz1_ > 0 ? static_cast<const casadi_int*>(sp_alg_v_) : nullptr;
    p_.sp_alg_r = nz1_ > 0 ? static_cast<const casadi_int*>(sp_alg_r_) : nullptr;
    p_.alg_prinv = get_ptr(alg_prinv_);
    p_.alg_pc = get_ptr(alg_pc_);
    p_.abstol = abstol_;
    p_.reltol = reltol_;
    p_.safety = 0.9;
    p_.fac_min = 0.2;
    p_.fac_max = 6;
    p_.h_min = min_step_size_;
    p_.h_max = max_step_size_;
    p_.h_init = step0_;
    p_.max_num_steps = max_num_steps_;
    p_.fixed_steps = fixed_steps_;
    p_.init_tol = init_tol_;
    p_.init_max_iter = init_max_iter_;
    p_.tape_size = tape_capacity();
    casadi_rosenbrock_setup(&p_);
  }

  void Rosenbrock::set_work(void* mem, const double**& arg, double**& res,
      casadi_int*& iw, double*& w) const {
    auto m = static_cast<RosenbrockMemory*>(mem);
    Integrator::set_work(mem, arg, res, iw, w);
    m->d.prob = &p_;
    casadi_rosenbrock_set_work(&m->d, &arg, &res, &iw, &w);
    if (replay_) {
      m->d_nom.prob = &p_nom_;
      casadi_rosenbrock_set_work(&m->d_nom, &arg, &res, &iw, &w);
      m->d.nom = &m->d_nom;
    }
  }

  int Rosenbrock::init_mem(void* mem) const {
    if (Integrator::init_mem(mem)) return 1;
    auto m = static_cast<RosenbrockMemory*>(mem);
    casadi_rosenbrock_clear_stats(&m->d);
    return 0;
  }

  void Rosenbrock::reset(IntegratorMemory* mem, bool first_call) const {
    auto m = static_cast<RosenbrockMemory*>(mem);
    casadi_rosenbrock_data<double>* d = &m->d;
    // State y = [x; z; q]
    casadi_copy(m->x, nx_, d->y);
    casadi_copy(m->z, nz_, d->y + nx_);
    casadi_copy(m->q, nq_, d->y + nx_ + nz_);
    casadi_rosenbrock_reset(d, m->t, m->u, first_call);
  }

  int Rosenbrock::advance_noevent(IntegratorMemory* mem) const {
    auto m = static_cast<RosenbrockMemory*>(mem);
    casadi_rosenbrock_data<double>* d = &m->d;
    casadi_rosenbrock_advance(d, m->t_next, m->t_stop, nrx_ > 0);
    run(m);
    if (d->status != ROSENBROCK_SUCCESS) {
      if (verbose_) {
        casadi_warning(name_ + ": integration failed at t = " + str(d->t) + ", status "
          + str(static_cast<casadi_int>(d->status)));
      }
      return 1;
    }
    casadi_copy(d->y_out, nx_, m->x);
    casadi_copy(d->y_out + nx_, nz_, m->z);
    casadi_copy(d->y_out + nx_ + nz_, nq_, m->q);
    return 0;
  }

  void Rosenbrock::resetB(IntegratorMemory* mem) const {
    auto m = static_cast<RosenbrockMemory*>(mem);
    casadi_rosenbrock_reset_backward(&m->d);
  }

  void Rosenbrock::impulseB(IntegratorMemory* mem,
      const double* adj_x, const double* adj_z, const double* adj_q) const {
    auto m = static_cast<RosenbrockMemory*>(mem);
    double* seed = m->d.seed;
    casadi_copy(adj_x, nrx_, seed);
    casadi_copy(adj_z, nrz_, seed + nrx_);
    casadi_copy(adj_q, nrp_, seed + nrx_ + nrz_);
    casadi_rosenbrock_impulse(&m->d, m->t);
    run(m);
    casadi_assert(m->d.status == ROSENBROCK_SUCCESS, name_ + ": adjoint seed at t = "
      + str(m->t) + " failed, status " + str(static_cast<casadi_int>(m->d.status)));
  }

  casadi_int Rosenbrock::next_stop(casadi_int k, const double* u) const {
    // Steps cross output times unless the nominal controls change: derivatives share the mesh
    return casadi_rosenbrock_next_stop(u, nu_, nu_err_, k, nt());
  }

  void Rosenbrock::retreat(IntegratorMemory* mem, const double* u,
      double* adj_x, double* adj_p, double* adj_u) const {
    auto m = static_cast<RosenbrockMemory*>(mem);
    casadi_rosenbrock_data<double>* d = &m->d;
    // The controls of each step are taped
    casadi_rosenbrock_retreat(d, m->t, m->t_next);
    run(m);
    casadi_assert(d->status == ROSENBROCK_SUCCESS,
      name_ + ": adjoint sweep failed, status " + str(static_cast<casadi_int>(d->status)));
    casadi_copy(d->ry, nrx_, adj_x);
    casadi_copy(d->adj_p, nrq_, adj_p);
    casadi_copy(d->adj_u, nuq_, adj_u);
  }

  int Rosenbrock::task(RosenbrockMemory* m) const {
    casadi_rosenbrock_data<double>* d = &m->d;
    if (d->task == ROSENBROCK_TAPE_RESERVE) {
      // The tape persists between evaluations: grow it only when too small
      if (m->tape.size() < static_cast<size_t>(d->tape_need)) {
        m->tape.resize(std::max(static_cast<size_t>(d->tape_need), 2 * m->tape.size()));
      }
      d->tape = get_ptr(m->tape);
      d->tape_cap = m->tape.size();
      return 0;
    }
    return eval_task(m, d, calls_);
  }

  int Rosenbrock::eval_task(RosenbrockMemory* m, casadi_rosenbrock_data<double>* d,
      const std::vector<std::vector<TaskCall>>& calls) const {
    // Function evaluations requested by the kernel, as described by the task table
    for (const TaskCall& c : calls.at(d->task)) {
      std::fill_n(m->arg, c.n_in, nullptr);
      std::fill_n(m->res, c.n_out, nullptr);
      for (const Slot& e : c.arg) m->arg[e.i] = slot_ptr(m, d, e);
      for (const Slot& e : c.res) m->res[e.i] = slot_ptr(m, d, e);
      if (calc_function(m, c.fname)) return 1;
    }
    return 0;
  }

  double* Rosenbrock::slot_ptr(RosenbrockMemory* m, casadi_rosenbrock_data<double>* d,
      const Slot& e) const {
    // Buffer of a function argument or result, as slot_expr in generated code
    double* b = nullptr;
    switch (e.src) {
      case SRC_T: b = &d->t_arg; break;
      case SRC_Y: b = const_cast<double*>(d->y_arg); break;
      case SRC_P: b = m->p; break;
      case SRC_U: b = d->u_step; break;
      case SRC_RES: b = d->res; break;
      case SRC_W: b = const_cast<double*>(d->w); break;
      case SRC_JAC: b = d->jac; break;
      case SRC_JAC_FWD: b = d->jac_fwd; break;
      case SRC_JAC_ADJ: b = d->jac_adj; break;
      case SRC_JAC_ADJ_FWD: b = d->jac_adj_fwd; break;
      case SRC_VJP_X: b = d->vjp_x; break;
      case SRC_VJP_Z: b = d->vjp_z; break;
      case SRC_VJP_P: b = d->vjp_p; break;
      case SRC_VJP_U: b = d->vjp_u; break;
    }
    return b + e.off;
  }

  std::string Rosenbrock::slot_expr(const std::string& d, const Slot& e) {
    // Buffer of a function argument or result in generated code, as slot_ptr
    std::string b;
    switch (e.src) {
      case SRC_T: b = "&" + d + ".t_arg"; break;
      case SRC_Y: b = d + ".y_arg"; break;
      case SRC_P: b = "pv"; break;
      case SRC_U: b = d + ".u_step"; break;
      case SRC_RES: b = d + ".res"; break;
      case SRC_W: b = d + ".w"; break;
      case SRC_JAC: b = d + ".jac"; break;
      case SRC_JAC_FWD: b = d + ".jac_fwd"; break;
      case SRC_JAC_ADJ: b = d + ".jac_adj"; break;
      case SRC_JAC_ADJ_FWD: b = d + ".jac_adj_fwd"; break;
      case SRC_VJP_X: b = d + ".vjp_x"; break;
      case SRC_VJP_Z: b = d + ".vjp_z"; break;
      case SRC_VJP_P: b = d + ".vjp_p"; break;
      case SRC_VJP_U: b = d + ".vjp_u"; break;
    }
    return e.off == 0 ? b : b + "+" + str(e.off);
  }

  std::vector<std::vector<Rosenbrock::TaskCall>> Rosenbrock::task_calls(
      const std::string& suffix) const {
    // For each task, the functions to call and the buffers of their arguments and results
    std::vector<std::vector<TaskCall>> calls(ROSENBROCK_TAPE_RESERVE);
    static_assert(static_cast<int>(BDYN_U) == static_cast<int>(DYN_U)
      && static_cast<int>(BDYN_T) == static_cast<int>(DYN_T),
      "Backward DAE inputs start with t, x, z, p, u");
    const casadi_int nj = ROSENBROCK_NUM_JAC;
    auto call = [&](const std::string& fname) -> TaskCall {
      const Function& f = get_function(fname);
      return {fname + suffix, f.n_in(), f.n_out(), {}, {}};
    };
    // Inputs t, x, z, p, u at offset i; forward seeds are the sensitivity parts
    auto dyn = [&](TaskCall& c, casadi_int i, bool fwd) -> void {
      if (!fwd) c.arg.push_back({i + DYN_T, SRC_T, 0});
      c.arg.push_back({i + DYN_X, SRC_Y, fwd ? nx1_ : 0});
      c.arg.push_back({i + DYN_Z, SRC_Y, nx_ + (fwd ? nz1_ : 0)});
      c.arg.push_back({i + DYN_P, SRC_P, fwd ? np1_ : 0});
      c.arg.push_back({i + DYN_U, SRC_U, fwd ? nu1_ : 0});
    };
    // Jacobian blocks in buffer src with ndir directions, block k at slot i + k
    auto blocks = [&](std::vector<Slot>& v, casadi_int i, Src src, casadi_int ndir) -> void {
      for (casadi_int k = 0; k < nj; ++k) v.push_back({i + k, src, ndir * p_.jac_off[k]});
    };
    // Right-hand side [ode; alg; quad] at (t, y)
    TaskCall rhs = call("rhsF");
    dyn(rhs, 0, false);
    rhs.res = {{0, SRC_RES, 0}, {1, SRC_RES, nx_}, {2, SRC_RES, nx_ + nz_}};
    calls[ROSENBROCK_EVAL_RHS].push_back(rhs);
    if (nfwd_ > 0) {
      TaskCall c = call(forward_name("rhsF", nfwd_));
      dyn(c, 0, false);
      for (casadi_int k = 0; k < 3; ++k) {
        c.arg.push_back({DYN_NUM_IN + k, SRC_RES, rhs.res[k].off});
      }
      dyn(c, DYN_NUM_IN + 3, true);
      c.res = {{0, SRC_RES, nx1_}, {1, SRC_RES, nx_ + nz1_}, {2, SRC_RES, nx_ + nz_ + nq1_}};
      calls[ROSENBROCK_EVAL_RHS].push_back(c);
    }
    // Jacobian blocks; forward derivatives give the augmented off-diagonal blocks
    TaskCall jac = call("jacF");
    dyn(jac, 0, false);
    blocks(jac.res, 0, SRC_JAC, 1);
    calls[ROSENBROCK_EVAL_JAC].push_back(jac);
    if (nfwd_ > 0) {
      TaskCall c = call(forward_name("jacF", nfwd_));
      dyn(c, 0, false);
      blocks(c.arg, DYN_NUM_IN, SRC_JAC, 1);
      dyn(c, DYN_NUM_IN + nj, true);
      blocks(c.res, 0, SRC_JAC_FWD, nfwd_);
      calls[ROSENBROCK_EVAL_JAC].push_back(c);
    }
    if (nadj_ > 0) {
      // [dF/dx; dF/dz; dF/dp; dF/du]^T w, with forward-over-adjoint
      const Src vjp[4] = {SRC_VJP_X, SRC_VJP_Z, SRC_VJP_P, SRC_VJP_U};
      const casadi_int n_vjp[4] = {nrx1_ * nadj_, nrz1_ * nadj_, nrq1_ * nadj_, nuq1_ * nadj_};
      TaskCall vjpb = call("vjpB");
      dyn(vjpb, 0, false);
      vjpb.arg.push_back({BDYN_ADJ_ODE, SRC_W, 0});
      vjpb.arg.push_back({BDYN_ADJ_ALG, SRC_W, nrx_});
      vjpb.arg.push_back({BDYN_ADJ_QUAD, SRC_W, nrx_ + nrz_});
      for (casadi_int k = 0; k < 4; ++k) vjpb.res.push_back({k, vjp[k], 0});
      calls[ROSENBROCK_EVAL_VJP].push_back(vjpb);
      if (nfwd_ > 0) {
        TaskCall c = call(forward_name("vjpB", nfwd_));
        casadi_int i = BDYN_NUM_IN + 4;
        c.arg = vjpb.arg;
        for (casadi_int k = 0; k < 4; ++k) c.arg.push_back({BDYN_NUM_IN + k, vjp[k], 0});
        dyn(c, i, true);
        c.arg.push_back({i + BDYN_ADJ_ODE, SRC_W, nrx1_ * nadj_});
        c.arg.push_back({i + BDYN_ADJ_ALG, SRC_W, nrx_ + nrz1_ * nadj_});
        c.arg.push_back({i + BDYN_ADJ_QUAD, SRC_W, nrx_ + nrz_ + nrp1_ * nadj_});
        for (casadi_int k = 0; k < 4; ++k) c.res.push_back({k, vjp[k], n_vjp[k]});
        calls[ROSENBROCK_EVAL_VJP].push_back(c);
      }
      // d/d(x, z, p, u) <jac_adj, Jacobian blocks>, with forward-over-adjoint
      TaskCall jacb = call("jacB");
      dyn(jacb, 0, false);
      blocks(jacb.arg, DYN_NUM_IN, SRC_JAC, 1);
      blocks(jacb.arg, DYN_NUM_IN + nj, SRC_JAC_ADJ, nadj_);
      for (casadi_int k = 0; k < 4; ++k) jacb.res.push_back({DYN_X + k, vjp[k], 0});
      calls[ROSENBROCK_EVAL_JAC_ADJ].push_back(jacb);
      if (nfwd_ > 0) {
        TaskCall c = call(forward_name("jacB", nfwd_));
        casadi_int i = DYN_NUM_IN + 2 * nj, i_fwd = i + DYN_NUM_IN;
        c.arg = jacb.arg;
        for (casadi_int k = 0; k < 4; ++k) c.arg.push_back({i + DYN_X + k, vjp[k], 0});
        dyn(c, i_fwd, true);
        blocks(c.arg, i_fwd + DYN_NUM_IN + nj, SRC_JAC_ADJ_FWD, nadj_ * nfwd_);
        for (casadi_int k = 0; k < 4; ++k) c.res.push_back({DYN_X + k, vjp[k], n_vjp[k]});
        calls[ROSENBROCK_EVAL_JAC_ADJ].push_back(c);
      }
    }
    return calls;
  }

  bool Rosenbrock::has_codegen() const {
    // Adjoint sensitivities need a tape of bounded size
    return nadj_ == 0 || tape_capacity() > 0;
  }

  void Rosenbrock::codegen_declarations(CodeGenerator& g) const {
    std::vector<std::string> fcn = {"rhsF", "jacF"};
    if (nadj_ > 0) {
      fcn.push_back("vjpB");
      fcn.push_back("jacB");
    }
    for (const std::string& f : fcn) {
      g.add_dependency(get_function(f));
      if (nfwd_ > 0) g.add_dependency(get_function(forward_name(f, nfwd_)));
    }
    if (replay_) {
      g.add_dependency(get_function("rhsF_nom"));
      g.add_dependency(get_function("jacF_nom"));
    }
  }

  void Rosenbrock::codegen_body(CodeGenerator& g) const {
    g.add_auxiliary(CodeGenerator::AUX_ROSENBROCK);
    const std::string nl = ";\n";
    g.local("p", "struct casadi_rosenbrock_prob");
    g.local("d", "struct casadi_rosenbrock_data");
    if (replay_) {
      g.local("pn", "struct casadi_rosenbrock_prob");
      g.local("dn", "struct casadi_rosenbrock_data");
    }
    g.local("sarg", "const casadi_real", "**");
    g.local("sres", "casadi_real", "**");
    g.local("pv", "casadi_real", "*");
    g.local("uv", "casadi_real", "*");
    g.local("t", "casadi_real");
    g.local("k", "casadi_int");
    g.local("k_stop", "casadi_int");

    g.comment("Problem structure");
    codegen_prob(g, "p");
    if (replay_) {
      g.comment("Nominal integration, whose steps are replayed");
      root()->codegen_prob(g, "pn");
      g << "pn.tape_size = 0;\n";
      g << "casadi_rosenbrock_setup(&pn);\n";
    }

    g.comment("Work vectors");
    g << "sarg = arg + " << INTEGRATOR_NUM_IN << "; sres = res + " << INTEGRATOR_NUM_OUT << nl;
    g << "pv = w; w += " << np_ << nl;
    g << "uv = w; w += " << nu_ << nl;
    g << "d.prob = &p" << nl;
    g << "casadi_rosenbrock_set_work(&d, &arg, &res, &iw, &w)" << nl;
    if (replay_) {
      g << "dn.prob = &pn" << nl;
      g << "casadi_rosenbrock_set_work(&dn, &arg, &res, &iw, &w)" << nl;
      g << "d.nom = &dn" << nl;
    }

    g.comment("Initial state, parameters");
    g << g.copy("arg[" + str(INTEGRATOR_X0) + "]", nx_, "d.y") << "\n";
    g << g.copy("arg[" + str(INTEGRATOR_Z0) + "]", nz_, "d.y+" + str(nx_)) << "\n";
    g << g.clear("d.y+" + str(nx_ + nz_), nq_) << "\n";
    g << g.copy("arg[" + str(INTEGRATOR_P) + "]", np_, "pv") << "\n";

    g.comment("Forward integration, output interval by output interval");
    std::string tout = g.constant(tout_);
    std::string u_k = "arg[" + str(INTEGRATOR_U) + "] ? arg[" + str(INTEGRATOR_U) + "] + k*"
      + str(nu_) + " : 0";
    g << "t = " << g.constant(t0_) << nl;
    g << "k_stop = -1" << nl;
    g << "for (k = 0; k < " << nt() << "; ++k) {\n";
    g << g.copy(u_k, nu_, "uv") << "\n";
    g.comment("Reset where the nominal controls change, as Integrator::eval");
    g << "if (k > k_stop) {\n";
    g << "k_stop = casadi_rosenbrock_next_stop(" << u_k << ", " << nu_ << ", " << nu_err_
      << ", k, " << nt() << ")" << nl;
    g << "casadi_rosenbrock_reset(&d, t, uv, k == 0)" << nl;
    g << "}\n";
    g << "casadi_rosenbrock_advance(&d, " << tout << "[k], " << tout << "[k_stop], "
      << (nrx_ > 0 ? 1 : 0) << ")" << nl;
    codegen_loop(g);
    g << "if (d.status != ROSENBROCK_SUCCESS) return 1" << nl;
    g << "if (res[" << INTEGRATOR_XF << "]) " << g.copy("d.y_out", nx_, "res["
      + str(INTEGRATOR_XF) + "] + k*" + str(nx_)) << "\n";
    g << "if (res[" << INTEGRATOR_ZF << "]) " << g.copy("d.y_out+" + str(nx_), nz_, "res["
      + str(INTEGRATOR_ZF) + "] + k*" + str(nz_)) << "\n";
    g << "if (res[" << INTEGRATOR_QF << "]) " << g.copy("d.y_out+" + str(nx_ + nz_), nq_,
      "res[" + str(INTEGRATOR_QF) + "] + k*" + str(nq_)) << "\n";
    g << "t = " << tout << "[k]" << nl;
    g << "}\n";

    if (nrx_ > 0) {
      g.comment("Discrete adjoint, output interval by output interval");
      g << "casadi_rosenbrock_reset_backward(&d)" << nl;
      g << "for (k = " << nt() - 1 << "; k >= 0; --k) {\n";
      g.comment("Impulse from the adjoint seeds at output time k");
      const casadi_int adj_in[3] = {INTEGRATOR_ADJ_XF, INTEGRATOR_ADJ_ZF, INTEGRATOR_ADJ_QF};
      const casadi_int adj_n[3] = {nrx_, nrz_, nrp_};
      const casadi_int adj_off[3] = {0, nrx_, nrx_ + nrz_};
      for (casadi_int i = 0; i < 3; ++i) {
        if (adj_n[i] == 0) continue;
        g << g.copy("arg[" + str(adj_in[i]) + "] ? arg[" + str(adj_in[i]) + "] + k*"
          + str(adj_n[i]) + " : 0", adj_n[i], "d.seed+" + str(adj_off[i])) << "\n";
      }
      g << "casadi_rosenbrock_impulse(&d, " << tout << "[k])" << nl;
      codegen_loop(g);
      g << "if (d.status != ROSENBROCK_SUCCESS) return 1" << nl;
      g << "casadi_rosenbrock_retreat(&d, " << tout << "[k], k > 0 ? " << tout << "[k-1] : "
        << g.constant(t0_) << ")" << nl;
      codegen_loop(g);
      g << "if (d.status != ROSENBROCK_SUCCESS) return 1" << nl;
      g << "if (res[" << INTEGRATOR_ADJ_U << "]) " << g.copy("d.adj_u", nuq_, "res["
        + str(INTEGRATOR_ADJ_U) + "] + k*" + str(nuq_)) << "\n";
      g << "}\n";
      g << "if (res[" << INTEGRATOR_ADJ_X0 << "]) " << g.copy("d.ry", nrx_, "res["
        + str(INTEGRATOR_ADJ_X0) + "]") << "\n";
      g << "if (res[" << INTEGRATOR_ADJ_P << "]) " << g.copy("d.adj_p", nrq_, "res["
        + str(INTEGRATOR_ADJ_P) + "]") << "\n";
      g.comment("Control adjoints per output interval rather than cumulative");
      g << "if (res[" << INTEGRATOR_ADJ_U << "]) for (k = 0; k < " << nt() - 1 << "; ++k) "
        << g.axpy(nuq_, "-1.", "res[" + str(INTEGRATOR_ADJ_U) + "] + (k+1)*" + str(nuq_),
        "res[" + str(INTEGRATOR_ADJ_U) + "] + k*" + str(nuq_)) << "\n";
    }
  }

  void Rosenbrock::codegen_loop(CodeGenerator& g) const {
    // Carry out the tasks requested by the kernel, for it or for the nominal integration
    g << "while (casadi_rosenbrock(&d)) {\n";
    if (replay_) {
      g << "if (d.active == &d) {\n";
      codegen_task(g, "d", calls_);
      g << "} else {\n";
      codegen_task(g, "dn", calls_nom_);
      g << "}\n";
    } else {
      codegen_task(g, "d", calls_);
    }
    g << "}\n";
  }

  void Rosenbrock::codegen_prob(CodeGenerator& g, const std::string& p) const {
    // Problem structure of the kernel in variable p, as set_rosenbrock_prob
    auto set = [&](const std::string& field, const std::string& value) -> void {
      g << p << "." << field << " = " << value << ";\n";
    };
    auto sp = [&](const Sparsity& s) -> std::string { return g.sparsity(s); };
    set("nx1", str(nx1_));
    set("nz1", str(nz1_));
    set("nq1", str(nq1_));
    set("np1", str(np1_));
    set("nu1", str(nu1_));
    set("nfwd", str(nfwd_));
    set("nadj", str(nadj_));
    set("nx_err", str(p_.nx_err));
    set("nz_err", str(p_.nz_err));
    set("nq_err", str(p_.nq_err));
    set("nu_err", str(p_.nu_err));
    set("n_stages", str(p_.n_stages));
    set("order_err", str(p_.order_err));
    set("gamma", g.constant(p_.gamma));
    set("a", g.constant(a_));
    set("c", g.constant(c_));
    set("alpha", g.constant(alpha_));
    set("gamma_sum", g.constant(gamma_sum_));
    set("m", g.constant(m_));
    set("m_err", g.constant(m_err_));
    set("n_dense", str(p_.n_dense));
    set("m_dense", g.constant(m_dense_));
    for (casadi_int k = 0; k < ROSENBROCK_NUM_JAC; ++k) {
      set("sp_jac[" + str(k) + "]", sp(sp_jac_[k]));
    }
    set("sp_iter", sp(sp_iter_));
    set("iter_map", g.constant(iter_map_));
    set("sp_iter_v", sp(sp_iter_v_));
    set("sp_iter_r", sp(sp_iter_r_));
    set("iter_prinv", g.constant(iter_prinv_));
    set("iter_pc", g.constant(iter_pc_));
    set("sp_alg_v", nz1_ > 0 ? sp(sp_alg_v_) : "0");
    set("sp_alg_r", nz1_ > 0 ? sp(sp_alg_r_) : "0");
    set("alg_prinv", nz1_ > 0 ? g.constant(alg_prinv_) : "0");
    set("alg_pc", nz1_ > 0 ? g.constant(alg_pc_) : "0");
    set("abstol", g.constant(p_.abstol));
    set("reltol", g.constant(p_.reltol));
    set("safety", g.constant(p_.safety));
    set("fac_min", g.constant(p_.fac_min));
    set("fac_max", g.constant(p_.fac_max));
    set("h_min", g.constant(p_.h_min));
    set("h_max", g.constant(p_.h_max));
    set("h_init", g.constant(p_.h_init));
    set("max_num_steps", str(p_.max_num_steps));
    set("fixed_steps", str(p_.fixed_steps));
    set("init_tol", g.constant(p_.init_tol));
    set("init_max_iter", str(p_.init_max_iter));
    set("tape_size", str(p_.tape_size));
    g << "casadi_rosenbrock_setup(&" << p << ");\n";
  }

  void Rosenbrock::codegen_task(CodeGenerator& g, const std::string& d,
      const std::vector<std::vector<TaskCall>>& calls) const {
    // Straight-line code for the task table, as eval_task: direct calls, no function pointers
    g << "switch (" << d << ".task) {\n";
    for (casadi_int task = 0; task < static_cast<casadi_int>(calls.size()); ++task) {
      if (calls[task].empty()) continue;
      g << "case " << task_name(task) << ":\n";
      // Slots set earlier in the same case need not be set again
      std::map<casadi_int, std::string> arg_set, res_set;
      for (const TaskCall& c : calls[task]) {
        std::vector<std::string> arg(c.n_in, "0"), res(c.n_out, "0");
        for (const Slot& e : c.arg) arg[e.i] = slot_expr(d, e);
        for (const Slot& e : c.res) res[e.i] = slot_expr(d, e);
        for (casadi_int i = 0; i < c.n_in; ++i) {
          if (arg_set.count(i) && arg_set[i] == arg[i]) continue;
          g << "sarg[" << i << "] = " << arg[i] << ";\n";
          arg_set[i] = arg[i];
        }
        for (casadi_int i = 0; i < c.n_out; ++i) {
          if (res_set.count(i) && res_set[i] == res[i]) continue;
          g << "sres[" << i << "] = " << res[i] << ";\n";
          res_set[i] = res[i];
        }
        // On failure, skip the remaining evaluations of the task
        std::string flag = g(get_function(c.fname), "sarg", "sres", "iw", "w");
        g << "if (" << flag << ") {\n" << d << ".fail = 1;\nbreak;\n}\n";
      }
      g << "break;\n";
    }
    g << "default: " << d << ".fail = 1;\n";
    g << "}\n";
  }

  std::string Rosenbrock::task_name(casadi_int task) {
    switch (task) {
      case ROSENBROCK_EVAL_RHS: return "ROSENBROCK_EVAL_RHS";
      case ROSENBROCK_EVAL_JAC: return "ROSENBROCK_EVAL_JAC";
      case ROSENBROCK_EVAL_VJP: return "ROSENBROCK_EVAL_VJP";
      case ROSENBROCK_EVAL_JAC_ADJ: return "ROSENBROCK_EVAL_JAC_ADJ";
      default: return "ROSENBROCK_TAPE_RESERVE";
    }
  }

  Dict Rosenbrock::get_stats(void* mem) const {
    Dict stats = Integrator::get_stats(mem);
    auto m = static_cast<RosenbrockMemory*>(mem);
    stats["n_steps"] = m->d.n_steps;
    stats["n_reject"] = m->d.n_reject;
    stats["n_rhs"] = m->d.n_rhs;
    stats["n_jac"] = m->d.n_jac;
    stats["n_factor"] = m->d.n_factor;
    stats["n_solve"] = m->d.n_solve;
    stats["n_vjp"] = m->d.n_vjp;
    return stats;
  }

  void Rosenbrock::print_stats(IntegratorMemory* mem) const {
    auto m = static_cast<RosenbrockMemory*>(mem);
    print("Rosenbrock: %lld steps, %lld rejected, %lld rhs, %lld jac, %lld factorizations\n",
      static_cast<long long>(m->d.n_steps), static_cast<long long>(m->d.n_reject),
      static_cast<long long>(m->d.n_rhs), static_cast<long long>(m->d.n_jac),
      static_cast<long long>(m->d.n_factor));
  }

  Rosenbrock::Rosenbrock(DeserializingStream& s) : Integrator(s) {
    s.version("Rosenbrock", 1);
    s.unpack("Rosenbrock::abstol", abstol_);
    s.unpack("Rosenbrock::reltol", reltol_);
    s.unpack("Rosenbrock::step0", step0_);
    s.unpack("Rosenbrock::max_step_size", max_step_size_);
    s.unpack("Rosenbrock::min_step_size", min_step_size_);
    s.unpack("Rosenbrock::init_tol", init_tol_);
    s.unpack("Rosenbrock::max_num_steps", max_num_steps_);
    s.unpack("Rosenbrock::fixed_steps", fixed_steps_);
    s.unpack("Rosenbrock::init_max_iter", init_max_iter_);
    s.unpack("Rosenbrock::max_tape_steps", max_tape_steps_);
    s.unpack("Rosenbrock::quad_err_con", quad_err_con_);
    s.unpack("Rosenbrock::nx_err", nx_err_);
    s.unpack("Rosenbrock::nz_err", nz_err_);
    s.unpack("Rosenbrock::nq_err", nq_err_);
    s.unpack("Rosenbrock::nu_err", nu_err_);
    s.unpack("Rosenbrock::sp_jac", sp_jac_);
    s.unpack("Rosenbrock::sp_iter", sp_iter_);
    s.unpack("Rosenbrock::iter_map", iter_map_);
    s.unpack("Rosenbrock::sp_iter_v", sp_iter_v_);
    s.unpack("Rosenbrock::sp_iter_r", sp_iter_r_);
    s.unpack("Rosenbrock::iter_prinv", iter_prinv_);
    s.unpack("Rosenbrock::iter_pc", iter_pc_);
    s.unpack("Rosenbrock::sp_alg_v", sp_alg_v_);
    s.unpack("Rosenbrock::sp_alg_r", sp_alg_r_);
    s.unpack("Rosenbrock::alg_prinv", alg_prinv_);
    s.unpack("Rosenbrock::alg_pc", alg_pc_);
    set_rosenbrock_prob();
    set_nominal_prob();
  }

  void Rosenbrock::serialize_body(SerializingStream &s) const {
    Integrator::serialize_body(s);
    s.version("Rosenbrock", 1);
    s.pack("Rosenbrock::abstol", abstol_);
    s.pack("Rosenbrock::reltol", reltol_);
    s.pack("Rosenbrock::step0", step0_);
    s.pack("Rosenbrock::max_step_size", max_step_size_);
    s.pack("Rosenbrock::min_step_size", min_step_size_);
    s.pack("Rosenbrock::init_tol", init_tol_);
    s.pack("Rosenbrock::max_num_steps", max_num_steps_);
    s.pack("Rosenbrock::fixed_steps", fixed_steps_);
    s.pack("Rosenbrock::init_max_iter", init_max_iter_);
    s.pack("Rosenbrock::max_tape_steps", max_tape_steps_);
    s.pack("Rosenbrock::quad_err_con", quad_err_con_);
    s.pack("Rosenbrock::nx_err", nx_err_);
    s.pack("Rosenbrock::nz_err", nz_err_);
    s.pack("Rosenbrock::nq_err", nq_err_);
    s.pack("Rosenbrock::nu_err", nu_err_);
    s.pack("Rosenbrock::sp_jac", sp_jac_);
    s.pack("Rosenbrock::sp_iter", sp_iter_);
    s.pack("Rosenbrock::iter_map", iter_map_);
    s.pack("Rosenbrock::sp_iter_v", sp_iter_v_);
    s.pack("Rosenbrock::sp_iter_r", sp_iter_r_);
    s.pack("Rosenbrock::iter_prinv", iter_prinv_);
    s.pack("Rosenbrock::iter_pc", iter_pc_);
    s.pack("Rosenbrock::sp_alg_v", sp_alg_v_);
    s.pack("Rosenbrock::sp_alg_r", sp_alg_r_);
    s.pack("Rosenbrock::alg_prinv", alg_prinv_);
    s.pack("Rosenbrock::alg_pc", alg_pc_);
  }

} // namespace casadi

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


#include "runtime_integrator.hpp"
#include "casadi_misc.hpp"

namespace casadi {

RuntimeIntegrator::RuntimeIntegrator(const std::string& name, const Function& dae, double t0,
    const std::vector<double>& tout) : Integrator(name, dae, t0, tout) {
}

RuntimeIntegrator::~RuntimeIntegrator() {
}

const Options RuntimeIntegrator::options_
= {{&Integrator::options_},
    {{"control_group_size",
      {OT_INTVECTOR,
      "Number of output intervals in each control group, over which the controls, one "
      "column each, are constant: a run of steps [default: 1 for each output interval]"}},
    {"step_size",
      {OT_DOUBLE,
      "Fixed steps: each run of steps, over a control group or up to the end without "
      "controls, is divided into equal steps of at most step_size [0: adaptive steps]"}},
    {"abstol",
      {OT_DOUBLE,
      "Absolute tolerance [1e-8]"}},
    {"reltol",
      {OT_DOUBLE,
      "Relative tolerance [1e-6]"}},
    {"max_num_steps",
      {OT_INT,
      "Adaptive steps: maximum number of steps per output interval [10000]"}},
    {"step0",
      {OT_DOUBLE,
      "Adaptive steps: initial step size [default: estimated]"}},
    {"max_step_size",
      {OT_DOUBLE,
      "Adaptive steps: maximum step size [default: unbounded]"}},
    {"min_step_size",
      {OT_DOUBLE,
      "Adaptive steps: minimum step size [1e-14]"}},
    {"quad_err_con",
      {OT_BOOL,
      "Adaptive steps: include the quadratures in the error control [true]"}},
    {"max_tape_steps",
      {OT_INT,
      "Adaptive steps: capacity, in steps, of a tape in the work vector, for adjoint "
      "sensitivities; needed for code generation [0: tape of unbounded size]. "
      "With fixed steps, the tape is sized exactly."}}
    }
};

void RuntimeIntegrator::init(const Dict& opts) {
  // Control groups, before the base class sizes the controls; a single one without controls
  group_start_ = range(nt() + 1);
  if (oracle_.numel_in(DYN_U) == 0) group_start_ = {0, nt()};
  auto cgs = opts.find("control_group_size");
  if (cgs != opts.end()) {
    group_start_ = {0};
    for (casadi_int n : cgs->second.to_int_vector()) {
      casadi_assert(n > 0, "'control_group_size' entries must be positive");
      group_start_.push_back(group_start_.back() + n);
    }
    casadi_assert(group_start_.back() == nt(), "'control_group_size' must sum to the number "
      "of output times, " + str(nt()) + ", got " + str(group_start_.back()));
  }

  // Call the base class init
  Integrator::init(opts);

  // Controls and control adjoints by output interval, for sparsity propagation
  if (ngroup() < nt()) alloc_w((nu_ + nuq_) * nt(), true);

  // Default options
  step_size_ = 0;
  abstol_ = 1e-8;
  reltol_ = 1e-6;
  max_num_steps_ = 10000;
  step0_ = 0;
  max_step_size_ = 0;
  min_step_size_ = 1e-14;
  quad_err_con_ = true;
  max_tape_steps_ = 0;

  // Read options
  for (auto&& op : opts) {
    if (op.first == "step_size") {
      step_size_ = op.second;
    } else if (op.first == "abstol") {
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
    } else if (op.first == "quad_err_con") {
      quad_err_con_ = op.second;
    } else if (op.first == "max_tape_steps") {
      max_tape_steps_ = op.second;
    }
  }
  casadi_assert(step_size_ >= 0, "'step_size' must be nonnegative");

  // Error control on the nominal states only: derivative integrators share the mesh
  if (derivative_of_.is_null()) {
    nx_err_ = nx1_;
    nz_err_ = nz1_;
    nq_err_ = nq1_;
  } else {
    const RuntimeIntegrator* d = derivative_of_.get<RuntimeIntegrator>();
    casadi_assert_dev(d != nullptr);
    nx_err_ = d->nx_err_;
    nz_err_ = d->nz_err_;
    nq_err_ = d->nq_err_;
  }

  // Right-hand side [ode; alg; quad] and its vector-Jacobian product, one function each
  create_function("rhsF", dyn_in(), {"ode", "alg", "quad"});
  if (nadj_ > 0) {
    create_function(rdae_, "vjpB", bdyn_in(), {"adj_x", "adj_z", "adj_p", "adj_u"});
  }

  // Zero-crossing functions and their time derivative; the jump at an event, its derivatives
  if (ne_ > 0) {
    create_function("zeroF", dyn_in(), {"zero"});
    create_forward("zeroF", 1);
    if (has_function("transition")) {
      Function jump = jump_function();
      set_function(jump, "jumpF");
      if (nfwd_ > 0) create_forward("jumpF", nfwd_);
      if (nadj_ > 0) {
        // Adjoint seeds of x, z, q to adjoints of x, z, p, u
        std::vector<MX> in = jump.mx_in(), out = jump(in), seed;
        for (casadi_int i = 0; i < jump.n_out(); ++i) {
          Sparsity sp = repmat(jump.sparsity_out(i), 1, nadj_);
          seed.push_back(MX::sym("adj_" + jump.name_out(i), sp));
        }
        std::vector<MX> arg = in;
        arg.insert(arg.end(), out.begin(), out.end());
        arg.insert(arg.end(), seed.begin(), seed.end());
        std::vector<MX> adj = jump.reverse(nadj_)(arg);
        in.insert(in.end(), seed.begin(), seed.end());
        set_function(Function("jumpB", in, {adj[1], adj[2], adj[4], adj[5]}), "jumpB");
        if (nfwd_ > 0) create_forward("jumpB", nfwd_);
      }
    }
  }

  // Jacobian blocks w.r.t. x, z and time, in the order of casadi_ivp_jac_t
  sp_jac_.clear();
  if (needs_jacobian()) {
    Function jac = create_function("jacF", {"t", "x", "z", "p", "u"},
      {"jac:ode:x", "jac:alg:x", "jac:quad:x", "jac:ode:z", "jac:alg:z", "jac:quad:z",
       "jac:ode:t", "jac:alg:t", "jac:quad:t"});
    for (casadi_int k = 0; k < IVP_NUM_JAC; ++k) sp_jac_.push_back(jac.sparsity_out(k));
    if (nadj_ > 0) set_function(jac.reverse(nadj_), "jacB");
  }

  // Forward sensitivities of all of these
  if (nfwd_ > 0) {
    create_forward("rhsF", nfwd_);
    if (nadj_ > 0) create_forward("vjpB", nfwd_);
    if (needs_jacobian()) {
      create_forward("jacF", nfwd_);
      if (nadj_ > 0) create_forward("jacB", nfwd_);
    }
  }
}

Function RuntimeIntegrator::jump_function() const {
  // State after the event selected by w as a function of the state before it: moved along
  // the trajectory to the event surface, through the transition and back to t, to second
  // order in the time s to the surface. At an event s vanishes, so the value is that of the
  // transition and the first and second derivatives are those of the event, including the
  // dependence of its time. Algebraic states are consistent along the way: the given ones
  // are only starting points
  MX t = MX::sym("t"), s = MX::sym("s"), w = MX::sym("w", ne_), q = MX::sym("q", nq1_);
  MX x = MX::sym("x", oracle_.sparsity_in(DYN_X)), z = MX::sym("z", oracle_.sparsity_in(DYN_Z));
  MX p = MX::sym("p", oracle_.sparsity_in(DYN_P)), u = MX::sym("u", oracle_.sparsity_in(DYN_U));
  const Function& tr = get_function("transition");
  // Argument i of f, empty if f does not depend on it
  auto arg = [](const Function& f, casadi_int i, const MX& v) -> MX {
    return f.nnz_in(i) == 0 ? MX(f.sparsity_in(i)) : v;
  };
  auto dae = [&](const MX& tt, const MX& xx, const MX& zz) -> std::vector<MX> {
    return oracle_(std::vector<MX>{arg(oracle_, DYN_T, tt), xx, zz, p, u});
  };
  // Algebraic states consistent with (tt, xx), from zz
  Function z_solve;
  if (nz1_ > 0) {
    MX zs = MX::sym("z", z.sparsity()), ts = MX::sym("t"), xs = MX::sym("x", x.sparsity());
    Function alg("alg", {zs, ts, xs, p, u}, {dae(ts, xs, zs)[DYN_ALG]});
    z_solve = rootfinder("z_consistent", "fast_newton", alg,
      Dict{{"abstol", 1e-14}, {"abstolStep", 1e-14}, {"max_iter", 50}});
  }
  auto consistent = [&](const MX& tt, const MX& xx, const MX& zz) -> MX {
    return nz1_ == 0 ? zz : z_solve(std::vector<MX>{zz, tt, xx, p, u}).at(0);
  };
  // Rates of x, q at (tt, xx) and their derivatives along the trajectory
  auto rates = [&](const MX& tt, const MX& xx, const MX& zz, std::vector<MX>& r) -> void {
    std::vector<MX> r0 = dae(tt, xx, consistent(tt, xx, zz));
    MX x1 = xx + s * r0[DYN_ODE];
    std::vector<MX> r1 = dae(tt + s, x1, consistent(tt + s, x1, zz));
    r = MX::substitute({r0[DYN_ODE], r0[DYN_QUAD], MX::jacobian(r1[DYN_ODE], s),
      MX::jacobian(r1[DYN_QUAD], s)}, {s}, {MX(0)});
  };
  // Time to the surface of the selected zero-crossing function
  std::vector<MX> r;
  rates(t, x, z, r);
  MX x_s = x + s * r[0] + 0.5 * s * s * r[2];
  MX g = MX::dot(w, dae(t + s, x_s, consistent(t + s, x_s, z))[DYN_ZERO]);
  MX g_s = MX::jacobian(g, s);
  std::vector<MX> gv = MX::substitute({g, g_s, MX::jacobian(g_s, s)}, {s}, {MX(0)});
  MX s1 = -gv[0] / gv[1];
  MX ss = s1 - 0.5 * gv[2] * s1 * s1 / gv[1];
  // Transition on the surface
  std::vector<double> index(ne_);
  for (casadi_int i = 0; i < ne_; ++i) index[i] = static_cast<double>(i);
  MX x_ev = x + ss * r[0] + 0.5 * ss * ss * r[2];
  MX q_ev = q + ss * r[1] + 0.5 * ss * ss * r[3];
  std::vector<MX> post = tr(std::vector<MX>{arg(tr, EVENT_INDEX, MX::dot(w, MX(DM(index)))),
    arg(tr, EVENT_T, t + ss), x_ev, consistent(t + ss, x_ev, z), arg(tr, EVENT_P, p),
    arg(tr, EVENT_U, u)});
  // Back to t along the trajectory after the event
  MX z_post = consistent(t + ss, post[EVENT_POST_X], post[EVENT_POST_Z]);
  rates(t + ss, post[EVENT_POST_X], z_post, r);
  MX x_post = post[EVENT_POST_X] - ss * r[0] + 0.5 * ss * ss * r[2];
  MX q_post = q_ev - ss * r[1] + 0.5 * ss * ss * r[3];
  return Function("jumpF", {t, x, z, q, p, u, w}, {x_post, z_post, q_post},
    {"t", "x", "z", "q", "p", "u", "w"}, {"post_x", "post_z", "post_q"});
}

const RuntimeIntegrator* RuntimeIntegrator::root() const {
  // Integrator at the start of the chain of derivatives
  const RuntimeIntegrator* r = this;
  while (!r->derivative_of_.is_null()) r = r->derivative_of_.get<RuntimeIntegrator>();
  return r;
}

void RuntimeIntegrator::set_ivp_prob() {
  // Problem structure for the shared part of the kernel
  p_.nx1 = nx1_;
  p_.nz1 = nz1_;
  p_.nq1 = nq1_;
  p_.np1 = np1_;
  p_.nu1 = nu1_;
  p_.nfwd = nfwd_;
  p_.nadj = nadj_;
  p_.n_stages = alpha_.size();
  p_.a = get_ptr(a_);
  p_.alpha = get_ptr(alpha_);
  p_.m = get_ptr(m_);
  p_.n_dense = n_dense_;
  p_.m_dense = get_ptr(m_dense_);
  for (casadi_int k = 0; k < IVP_NUM_JAC; ++k) {
    p_.sp_jac[k] = sp_jac_.empty() ? nullptr : static_cast<const casadi_int*>(sp_jac_[k]);
  }
  p_.tape_size = tape_capacity();
  p_.step_size = step_size_;
  p_.nx_err = nx_err_;
  p_.nz_err = nz_err_;
  p_.nq_err = quad_err_con_ ? nq_err_ : 0;
  p_.abstol = abstol_;
  p_.reltol = reltol_;
  p_.order_err = order_err_;
  p_.m_err = m_err_.empty() ? nullptr : get_ptr(m_err_);
  p_.safety = 0.9;
  p_.fac_min = 0.2;
  p_.fac_max = 6;
  p_.h_min = min_step_size_;
  p_.h_max = max_step_size_;
  p_.h_init = step0_;
  p_.max_num_steps = max_num_steps_;
  p_.ne = ne_;
  p_.has_transition = has_function("transition");
  p_.event_tol = event_tol_;
  p_.event_acceptable_tol = std::isinf(event_acceptable_tol_) ? -1 : event_acceptable_tol_;
  p_.max_event_iter = max_event_iter_;
  p_.max_events = max_events_;
  casadi_ivp_setup(&p_);
  calls_ = task_calls("");
}

void RuntimeIntegrator::set_nominal_prob() {
  // Forward directions folded into the DAE change the arithmetic of the nominal states: with
  // adaptive steps, replay the steps of a nominal integration with the original functions
  const RuntimeIntegrator* r = root();
  if (!has_function("rhsF_nom") && r != this && r->nx1_ != nx1_ && step_size_ == 0) {
    for (auto&& calls : r->task_calls("")) {
      for (auto&& c : calls) {
        std::string fname = c.fname + "_nom";
        if (!has_function(fname)) set_function(r->get_function(c.fname), fname);
      }
    }
  }
  replay_ = has_function("rhsF_nom");
  if (!replay_) return;
  calls_nom_ = root()->task_calls("_nom");
  p_nom_ = root()->p_;
  p_nom_.tape_size = 0;
  casadi_ivp_setup(&p_nom_);
}

void RuntimeIntegrator::alloc_ivp_work() {
  // Shared kernel data, and of the nominal integration
  casadi_int sz_arg = 0, sz_res = 0, sz_iw = 0, sz_w = 0;
  casadi_ivp_work(&p_, &sz_arg, &sz_res, &sz_iw, &sz_w);
  if (replay_) casadi_ivp_work(&p_nom_, &sz_arg, &sz_res, &sz_iw, &sz_w);
  alloc_w(sz_w, true);
  // Parameters and controls in generated code
  alloc_w(np_ + nu_, true);
}

casadi_int RuntimeIntegrator::tape_capacity() const {
  // Steps: a user bound with adaptive steps, a run of steps per control group with fixed steps
  casadi_int ny = nx_ + nz_ + nq_, n_step = max_tape_steps_, n_init, n_event = 0;
  if (step_size_ > 0) {
    n_step = 0;
    for (casadi_int j = 0; j < ngroup(); ++j) {
      double t = j == 0 ? t0_ : tout_[group_start_[j] - 1];
      n_step += casadi_ivp_steps(tout_[group_start_[j + 1] - 1] - t, step_size_);
    }
  }
  if (nadj_ == 0 || n_step == 0) return 0;
  // Every output interval may start with a consistent initialization
  n_init = nz_ > 0 ? nt() : 0;
  if (ne_ > 0) {
    // Event records, room for one more with each step; after an event, a consistent
    // initialization and, with fixed steps, a new run of steps
    n_event = max_events_ + 1;
    if (nz_ > 0) n_init += max_events_;
    if (step_size_ > 0) n_step += max_events_;
  }
  return n_step * (ny * (static_cast<casadi_int>(alpha_.size()) + 1) + nu_ + 5)
    + (n_init + n_event) * (ny + nu_ + 5);
}

void RuntimeIntegrator::set_work(void* mem, const double**& arg, double**& res,
    casadi_int*& iw, double*& w) const {
  auto m = static_cast<RuntimeIntegratorMemory*>(mem);
  Integrator::set_work(mem, arg, res, iw, w);
  m->d.prob = &p_;
  casadi_ivp_set_work(&m->d, &arg, &res, &iw, &w);
  if (replay_) {
    m->d_nom.prob = &p_nom_;
    casadi_ivp_set_work(&m->d_nom, &arg, &res, &iw, &w);
    m->d.nom = &m->d_nom;
  }
}

int RuntimeIntegrator::init_mem(void* mem) const {
  if (Integrator::init_mem(mem)) return 1;
  auto m = static_cast<RuntimeIntegratorMemory*>(mem);
  casadi_ivp_clear_stats(&m->d);
  return 0;
}

void RuntimeIntegrator::run(RuntimeIntegratorMemory* m) const {
  // Carry out the tasks requested by the kernel, for it or for the nominal integration
  while (kernel(m)) {
    casadi_ivp_data<double>* a = m->d.active;
    if (a == &m->d ? task(m) : eval_task(m, a, calls_nom_)) a->fail = 1;
  }
}

void RuntimeIntegrator::reset(IntegratorMemory* mem, bool first_call) const {
  auto m = static_cast<RuntimeIntegratorMemory*>(mem);
  casadi_ivp_data<double>* d = &m->d;
  // State y = [x; z; q]
  casadi_copy(m->x, nx_, d->y);
  casadi_copy(m->z, nz_, d->y + nx_);
  casadi_copy(m->q, nq_, d->y + nx_ + nz_);
  casadi_ivp_reset(d, m->t, m->u, first_call);
}

int RuntimeIntegrator::advance_noevent(IntegratorMemory* mem) const {
  auto m = static_cast<RuntimeIntegratorMemory*>(mem);
  casadi_ivp_data<double>* d = &m->d;
  casadi_ivp_advance(d, m->t_next, m->t_stop, nrx_ > 0);
  run(m);
  casadi_assert(d->status != IVP_TAPE_FAILED, name_ + ": the tape of " + str(max_tape_steps_)
    + " steps is full at t = " + str(d->t) + ", increase max_tape_steps");
  casadi_assert(d->status != IVP_TOO_MANY_EVENTS, name_ + ": more than max_events = "
    + str(max_events_) + " events, the last at t = " + str(d->t_ev));
  casadi_assert(d->status != IVP_EVENT_FAILED, name_ + ": event near t = " + str(d->t_ev)
    + " not located within event_acceptable_tol");
  if (d->status != IVP_SUCCESS) {
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

void RuntimeIntegrator::resetB(IntegratorMemory* mem) const {
  auto m = static_cast<RuntimeIntegratorMemory*>(mem);
  casadi_ivp_reset_backward(&m->d);
}

void RuntimeIntegrator::impulseB(IntegratorMemory* mem,
    const double* adj_x, const double* adj_z, const double* adj_q) const {
  auto m = static_cast<RuntimeIntegratorMemory*>(mem);
  double* seed = m->d.seed;
  casadi_copy(adj_x, nrx_, seed);
  casadi_copy(adj_z, nrz_, seed + nrx_);
  casadi_copy(adj_q, nrp_, seed + nrx_ + nrz_);
  casadi_ivp_impulse(&m->d, m->t);
  run(m);
  casadi_assert(m->d.status == IVP_SUCCESS, name_ + ": adjoint seed at t = "
    + str(m->t) + " failed, status " + str(static_cast<casadi_int>(m->d.status)));
}

Sparsity RuntimeIntegrator::get_sparsity_in(casadi_int i) {
  if (i == INTEGRATOR_U) return Sparsity::dense(nu1_, ngroup() * (1 + nfwd_));
  return Integrator::get_sparsity_in(i);
}

Sparsity RuntimeIntegrator::get_sparsity_out(casadi_int i) {
  if (i == INTEGRATOR_ADJ_U) return Sparsity::dense(nuq1_, nadj_ * (1 + nfwd_) * ngroup());
  return Integrator::get_sparsity_out(i);
}

int RuntimeIntegrator::eval(const double** arg, double** res,
    casadi_int* iw, double* w, void* mem) const {
  auto m = static_cast<RuntimeIntegratorMemory*>(mem);
  const double* u = arg[INTEGRATOR_U];
  const double* adj_xf = arg[INTEGRATOR_ADJ_XF];
  const double* adj_zf = arg[INTEGRATOR_ADJ_ZF];
  const double* adj_qf = arg[INTEGRATOR_ADJ_QF];
  double* xf = res[INTEGRATOR_XF];
  double* zf = res[INTEGRATOR_ZF];
  double* qf = res[INTEGRATOR_QF];
  double* adj_u = res[INTEGRATOR_ADJ_U];
  setup(m, arg + INTEGRATOR_NUM_IN, res + INTEGRATOR_NUM_OUT, iw, w);
  set_q(m, nullptr);
  set_x(m, arg[INTEGRATOR_X0]);
  set_z(m, arg[INTEGRATOR_Z0]);
  set_p(m, arg[INTEGRATOR_P]);
  m->num_events = 0;
  m->t = t0_;
  // Forward integration: a run of steps per control group
  for (casadi_int j = 0; j < ngroup(); ++j) {
    set_u(m, u ? u + j * nu_ : nullptr);
    m->t_stop = tout_[group_start_[j + 1] - 1];
    reset(m, j == 0);
    for (m->k = group_start_[j]; m->k < group_start_[j + 1]; ++m->k) {
      m->t_next = tout_[m->k];
      if (advance_noevent(m)) return 1;
      m->t = m->t_next;
      get_x(m, xf ? xf + m->k * nx_ : nullptr);
      get_z(m, zf ? zf + m->k * nz_ : nullptr);
      get_q(m, qf ? qf + m->k * nq_ : nullptr);
    }
  }
  // Discrete adjoint, output interval by output interval
  if (nrx_ > 0) {
    resetB(m);
    for (casadi_int j = ngroup(); j-- > 0; ) {
      // Control adjoints of this group only
      casadi_clear(m->d.adj_u, nuq_);
      for (m->k = group_start_[j + 1]; m->k-- > group_start_[j]; ) {
        m->t = tout_[m->k];
        impulseB(m, adj_xf ? adj_xf + m->k * nrx_ : nullptr,
          adj_zf ? adj_zf + m->k * nrz_ : nullptr, adj_qf ? adj_qf + m->k * nrp_ : nullptr);
        m->t_next = m->k > 0 ? tout_[m->k - 1] : t0_;
        retreat(m, nullptr, m->k == 0 ? res[INTEGRATOR_ADJ_X0] : nullptr,
          m->k == 0 ? res[INTEGRATOR_ADJ_P] : nullptr,
          adj_u && m->k == group_start_[j] ? adj_u + j * nuq_ : nullptr);
      }
    }
  }
  join_results(m);
  if (print_stats_) print_stats(m);
  return 0;
}

namespace {
  // Dependencies with a column per control group to one with a column per output interval
  void groups_to_intervals(const std::vector<casadi_int>& gs, const bvec_t* g, bvec_t* k,
      casadi_int n) {
    for (size_t j = 0; j + 1 < gs.size(); ++j) {
      for (casadi_int i = gs[j]; i < gs[j + 1]; ++i) std::copy_n(g + j * n, n, k + i * n);
    }
  }
  // Dependencies with a column per output interval added to one with a column per control group
  void intervals_to_groups(const std::vector<casadi_int>& gs, const bvec_t* k, bvec_t* g,
      casadi_int n) {
    for (size_t j = 0; j + 1 < gs.size(); ++j) {
      for (casadi_int i = gs[j]; i < gs[j + 1]; ++i) {
        for (casadi_int l = 0; l < n; ++l) g[j * n + l] |= k[i * n + l];
      }
    }
  }
}  // namespace

int RuntimeIntegrator::sp_forward(const bvec_t** arg, bvec_t** res,
    casadi_int* iw, bvec_t* w, void* mem) const {
  if (ne_ > 0 || ngroup() == nt()) return Integrator::sp_forward(arg, res, iw, w, mem);
  // Controls and control adjoints by output interval
  const bvec_t* u = arg[INTEGRATOR_U];
  bvec_t* adj_u = res[INTEGRATOR_ADJ_U];
  bvec_t* u_k = w; w += nu_ * nt();
  bvec_t* adj_u_k = w; w += nuq_ * nt();
  if (u) groups_to_intervals(group_start_, u, u_k, nu_);
  arg[INTEGRATOR_U] = u ? u_k : nullptr;
  res[INTEGRATOR_ADJ_U] = adj_u ? adj_u_k : nullptr;
  int flag = Integrator::sp_forward(arg, res, iw, w, mem);
  arg[INTEGRATOR_U] = u;
  res[INTEGRATOR_ADJ_U] = adj_u;
  if (adj_u) {
    std::fill_n(adj_u, nuq_ * ngroup(), 0);
    intervals_to_groups(group_start_, adj_u_k, adj_u, nuq_);
  }
  return flag;
}

int RuntimeIntegrator::sp_reverse(bvec_t** arg, bvec_t** res,
    casadi_int* iw, bvec_t* w, void* mem) const {
  if (ne_ > 0 || ngroup() == nt()) return Integrator::sp_reverse(arg, res, iw, w, mem);
  // Controls and control adjoints by output interval
  bvec_t* u = arg[INTEGRATOR_U];
  bvec_t* adj_u = res[INTEGRATOR_ADJ_U];
  bvec_t* u_k = w; w += nu_ * nt();
  bvec_t* adj_u_k = w; w += nuq_ * nt();
  std::fill_n(u_k, nu_ * nt(), 0);
  if (adj_u) {
    groups_to_intervals(group_start_, adj_u, adj_u_k, nuq_);
    std::fill_n(adj_u, nuq_ * ngroup(), 0);
  }
  arg[INTEGRATOR_U] = u ? u_k : nullptr;
  res[INTEGRATOR_ADJ_U] = adj_u ? adj_u_k : nullptr;
  int flag = Integrator::sp_reverse(arg, res, iw, w, mem);
  arg[INTEGRATOR_U] = u;
  res[INTEGRATOR_ADJ_U] = adj_u;
  if (u) intervals_to_groups(group_start_, u_k, u, nu_);
  return flag;
}

void RuntimeIntegrator::retreat(IntegratorMemory* mem, const double* u,
    double* adj_x, double* adj_p, double* adj_u) const {
  auto m = static_cast<RuntimeIntegratorMemory*>(mem);
  casadi_ivp_data<double>* d = &m->d;
  // The controls of each step are taped
  casadi_ivp_retreat(d, m->t, m->t_next);
  run(m);
  casadi_assert(d->status == IVP_SUCCESS,
    name_ + ": adjoint sweep failed, status " + str(static_cast<casadi_int>(d->status)));
  casadi_copy(d->ry, nrx_, adj_x);
  casadi_copy(d->adj_p, nrq_, adj_p);
  casadi_copy(d->adj_u, nuq_, adj_u);
}

int RuntimeIntegrator::task(RuntimeIntegratorMemory* m) const {
  casadi_ivp_data<double>* d = &m->d;
  if (d->task == IVP_TAPE_RESERVE) {
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

int RuntimeIntegrator::eval_task(RuntimeIntegratorMemory* m, casadi_ivp_data<double>* d,
    const std::vector<std::vector<TaskCall>>& calls) const {
  // Adaptive step size control, or function evaluations as described by the task table
  if (d->task == IVP_STEP_CONTROL) {
    casadi_ivp_control(d);
    return 0;
  }
  if (d->task == IVP_EVENT) {
    while (casadi_ivp_event(d)) {
      if (eval_task(m, d, calls)) d->fail = 1;
    }
    return 0;
  }
  if (d->task == IVP_EVENT_ADJ) {
    while (casadi_ivp_event_adj(d)) {
      if (eval_task(m, d, calls)) d->fail = 1;
    }
    return 0;
  }
  for (const TaskCall& c : calls.at(d->task)) {
    std::fill_n(m->arg, c.n_in, nullptr);
    std::fill_n(m->res, c.n_out, nullptr);
    for (const Slot& e : c.arg) m->arg[e.i] = slot_ptr(m, d, e);
    for (const Slot& e : c.res) m->res[e.i] = slot_ptr(m, d, e);
    if (calc_function(m, c.fname)) return 1;
  }
  return 0;
}

double* RuntimeIntegrator::slot_ptr(RuntimeIntegratorMemory* m, casadi_ivp_data<double>* d,
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
    case SRC_EV_ONE: b = &d->ev_one; break;
    case SRC_EV_W: b = d->ev_w; break;
    case SRC_YDOT: b = d->ydot; break;
    case SRC_EDOT: b = d->edot; break;
  }
  return b + e.off;
}

std::string RuntimeIntegrator::slot_expr(const std::string& d, const Slot& e) {
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
    case SRC_EV_ONE: b = "&" + d + ".ev_one"; break;
    case SRC_EV_W: b = d + ".ev_w"; break;
    case SRC_YDOT: b = d + ".ydot"; break;
    case SRC_EDOT: b = d + ".edot"; break;
  }
  return e.off == 0 ? b : b + "+" + str(e.off);
}

std::vector<std::vector<RuntimeIntegrator::TaskCall>> RuntimeIntegrator::task_calls(
    const std::string& suffix) const {
  // For each task, the functions to call and the buffers of their arguments and results
  std::vector<std::vector<TaskCall>> calls(IVP_TAPE_RESERVE);
  static_assert(static_cast<int>(BDYN_U) == static_cast<int>(DYN_U)
    && static_cast<int>(BDYN_T) == static_cast<int>(DYN_T),
    "Backward DAE inputs start with t, x, z, p, u");
  const casadi_int nj = IVP_NUM_JAC;
  // Results of the adjoint functions: adjoints of x, z, p, u and their forward derivatives
  const Src vjp[4] = {SRC_VJP_X, SRC_VJP_Z, SRC_VJP_P, SRC_VJP_U};
  const casadi_int n_vjp[4] = {nrx1_ * nadj_, nrz1_ * nadj_, nrq1_ * nadj_, nuq1_ * nadj_};
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
  calls[IVP_EVAL_RHS].push_back(rhs);
  if (nfwd_ > 0) {
    TaskCall c = call(forward_name("rhsF", nfwd_));
    dyn(c, 0, false);
    for (casadi_int k = 0; k < 3; ++k) {
      c.arg.push_back({DYN_NUM_IN + k, SRC_RES, rhs.res[k].off});
    }
    dyn(c, DYN_NUM_IN + 3, true);
    c.res = {{0, SRC_RES, nx1_}, {1, SRC_RES, nx_ + nz1_}, {2, SRC_RES, nx_ + nz_ + nq1_}};
    calls[IVP_EVAL_RHS].push_back(c);
  }
  if (needs_jacobian()) {
    // Jacobian blocks; forward derivatives give the augmented off-diagonal blocks
    TaskCall jac = call("jacF");
    dyn(jac, 0, false);
    blocks(jac.res, 0, SRC_JAC, 1);
    calls[IVP_EVAL_JAC].push_back(jac);
    if (nfwd_ > 0) {
      TaskCall c = call(forward_name("jacF", nfwd_));
      dyn(c, 0, false);
      blocks(c.arg, DYN_NUM_IN, SRC_JAC, 1);
      dyn(c, DYN_NUM_IN + nj, true);
      blocks(c.res, 0, SRC_JAC_FWD, nfwd_);
      calls[IVP_EVAL_JAC].push_back(c);
    }
  }
  if (ne_ > 0) {
    // Zero-crossing functions; derivative in the direction (1, ydot); forward sensitivities
    TaskCall zero = call("zeroF");
    dyn(zero, 0, false);
    zero.res = {{0, SRC_RES, 0}};
    calls[IVP_EVAL_ZERO].push_back(zero);
    TaskCall c = call(forward_name("zeroF", 1));
    dyn(c, 0, false);
    c.arg.push_back({DYN_NUM_IN, SRC_RES, 0});
    c.arg.push_back({DYN_NUM_IN + 1 + DYN_T, SRC_EV_ONE, 0});
    c.arg.push_back({DYN_NUM_IN + 1 + DYN_X, SRC_YDOT, 0});
    c.arg.push_back({DYN_NUM_IN + 1 + DYN_Z, SRC_YDOT, nx_});
    c.res = {{0, SRC_EDOT, 0}};
    calls[IVP_EVAL_ZERO].push_back(c);
    if (has_function("jumpF")) {
      // Inputs t, x, z, q, p, u, w of the jump at offset i; forward seeds are the sensitivities
      auto jump = [&](TaskCall& c, casadi_int i, bool fwd) -> void {
        if (!fwd) c.arg.push_back({i, SRC_T, 0});
        c.arg.push_back({i + 1, SRC_Y, fwd ? nx1_ : 0});
        c.arg.push_back({i + 2, SRC_Y, nx_ + (fwd ? nz1_ : 0)});
        c.arg.push_back({i + 3, SRC_Y, nx_ + nz_ + (fwd ? nq1_ : 0)});
        c.arg.push_back({i + 4, SRC_P, fwd ? np1_ : 0});
        c.arg.push_back({i + 5, SRC_U, fwd ? nu1_ : 0});
        if (!fwd) c.arg.push_back({i + 6, SRC_EV_W, 0});
      };
      TaskCall jf = call("jumpF");
      jump(jf, 0, false);
      jf.res = {{0, SRC_RES, 0}, {1, SRC_RES, nx_}, {2, SRC_RES, nx_ + nz_}};
      calls[IVP_EVAL_JUMP].push_back(jf);
      if (nfwd_ > 0) {
        c = call(forward_name("jumpF", nfwd_));
        jump(c, 0, false);
        for (casadi_int k = 0; k < 3; ++k) c.arg.push_back({7 + k, SRC_RES, jf.res[k].off});
        jump(c, 10, true);
        c.res = {{0, SRC_RES, nx1_}, {1, SRC_RES, nx_ + nz1_}, {2, SRC_RES, nx_ + nz_ + nq1_}};
        calls[IVP_EVAL_JUMP].push_back(c);
      }
      if (nadj_ > 0) {
        // Adjoint seeds of x, z, q to adjoints of x, z, p, u, with forward-over-adjoint
        TaskCall jb = call("jumpB");
        jump(jb, 0, false);
        jb.arg.push_back({7, SRC_W, 0});
        jb.arg.push_back({8, SRC_W, nrx_});
        jb.arg.push_back({9, SRC_W, nrx_ + nrz_});
        for (casadi_int k = 0; k < 4; ++k) jb.res.push_back({k, vjp[k], 0});
        calls[IVP_EVAL_JUMP_ADJ].push_back(jb);
        if (nfwd_ > 0) {
          c = call(forward_name("jumpB", nfwd_));
          c.arg = jb.arg;
          for (casadi_int k = 0; k < 4; ++k) c.arg.push_back({10 + k, vjp[k], 0});
          jump(c, 14, true);
          c.arg.push_back({21, SRC_W, nrx1_ * nadj_});
          c.arg.push_back({22, SRC_W, nrx_ + nrz1_ * nadj_});
          c.arg.push_back({23, SRC_W, nrx_ + nrz_ + nrp1_ * nadj_});
          for (casadi_int k = 0; k < 4; ++k) c.res.push_back({k, vjp[k], n_vjp[k]});
          calls[IVP_EVAL_JUMP_ADJ].push_back(c);
        }
      }
    }
  }
  if (nadj_ > 0) {
    // [dF/dx; dF/dz; dF/dp; dF/du]^T w, with forward-over-adjoint
    TaskCall vjpb = call("vjpB");
    dyn(vjpb, 0, false);
    vjpb.arg.push_back({BDYN_ADJ_ODE, SRC_W, 0});
    vjpb.arg.push_back({BDYN_ADJ_ALG, SRC_W, nrx_});
    vjpb.arg.push_back({BDYN_ADJ_QUAD, SRC_W, nrx_ + nrz_});
    for (casadi_int k = 0; k < 4; ++k) vjpb.res.push_back({k, vjp[k], 0});
    calls[IVP_EVAL_VJP].push_back(vjpb);
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
      calls[IVP_EVAL_VJP].push_back(c);
    }
    if (needs_jacobian()) {
      // d/d(x, z, p, u) <jac_adj, Jacobian blocks>, with forward-over-adjoint
      TaskCall jacb = call("jacB");
      dyn(jacb, 0, false);
      blocks(jacb.arg, DYN_NUM_IN, SRC_JAC, 1);
      blocks(jacb.arg, DYN_NUM_IN + nj, SRC_JAC_ADJ, nadj_);
      for (casadi_int k = 0; k < 4; ++k) jacb.res.push_back({DYN_X + k, vjp[k], 0});
      calls[IVP_EVAL_JAC_ADJ].push_back(jacb);
      if (nfwd_ > 0) {
        TaskCall c = call(forward_name("jacB", nfwd_));
        casadi_int i = DYN_NUM_IN + 2 * nj, i_fwd = i + DYN_NUM_IN;
        c.arg = jacb.arg;
        for (casadi_int k = 0; k < 4; ++k) c.arg.push_back({i + DYN_X + k, vjp[k], 0});
        dyn(c, i_fwd, true);
        blocks(c.arg, i_fwd + DYN_NUM_IN + nj, SRC_JAC_ADJ_FWD, nadj_ * nfwd_);
        for (casadi_int k = 0; k < 4; ++k) c.res.push_back({DYN_X + k, vjp[k], n_vjp[k]});
        calls[IVP_EVAL_JAC_ADJ].push_back(c);
      }
    }
  }
  return calls;
}

bool RuntimeIntegrator::has_codegen() const {
  // Adjoint sensitivities need a tape of bounded size
  return nadj_ == 0 || tape_capacity() > 0;
}

void RuntimeIntegrator::codegen_declarations(CodeGenerator& g) const {
  // The functions of the task tables
  for (auto* calls : {&calls_, &calls_nom_}) {
    for (auto&& task : *calls) {
      for (auto&& c : task) g.add_dependency(get_function(c.fname));
    }
  }
}

void RuntimeIntegrator::codegen_body(CodeGenerator& g) const {
  casadi_assert(has_codegen(), name_ + ": adjoint sensitivities with adaptive steps need "
    "max_tape_steps for code generation");
  g.add_auxiliary(CodeGenerator::AUX_IVP);
  if (step_size_ == 0) g.add_auxiliary(CodeGenerator::AUX_IVP_ADAPT);
  if (ne_ > 0) g.add_auxiliary(CodeGenerator::AUX_IVP_EVENT);
  const std::string nl = ";\n", method = plugin_name();
  g.local("p", "struct casadi_ivp_prob");
  g.local("d", "struct casadi_ivp_data");
  g.local("pe", "struct casadi_" + method + "_prob");
  g.local("e", "struct casadi_" + method + "_data");
  if (replay_) {
    g.local("pn", "struct casadi_ivp_prob");
    g.local("dn", "struct casadi_ivp_data");
    g.local("pen", "struct casadi_" + method + "_prob");
    g.local("en", "struct casadi_" + method + "_data");
  }
  g.local("sarg", "const casadi_real", "**");
  g.local("sres", "casadi_real", "**");
  g.local("pv", "casadi_real", "*");
  if (nu_ > 0) g.local("uv", "casadi_real", "*");
  g.local("j", "casadi_int");
  g.local("k", "casadi_int");

  g.comment("Problem structure");
  codegen_ivp_prob(g, "p", p_.tape_size);
  if (replay_) {
    g.comment("Nominal integration, whose steps are replayed");
    root()->codegen_ivp_prob(g, "pn", 0);
  }

  g.comment("Work vectors");
  g << "sarg = arg + " << INTEGRATOR_NUM_IN << "; sres = res + " << INTEGRATOR_NUM_OUT << nl;
  g << "pv = w; w += " << np_ << nl;
  if (nu_ > 0) g << "uv = w; w += " << nu_ << nl;
  g << "d.prob = &p" << nl;
  g << "casadi_ivp_set_work(&d, &arg, &res, &iw, &w)" << nl;
  if (replay_) {
    g << "dn.prob = &pn" << nl;
    g << "casadi_ivp_set_work(&dn, &arg, &res, &iw, &w)" << nl;
    g << "d.nom = &dn" << nl;
  }
  g.comment("Method");
  codegen_method(g);

  g.comment("Initial state, parameters");
  g << g.copy("arg[" + str(INTEGRATOR_X0) + "]", nx_, "d.y") << "\n";
  g << g.copy("arg[" + str(INTEGRATOR_Z0) + "]", nz_, "d.y+" + str(nx_)) << "\n";
  g << g.clear("d.y+" + str(nx_ + nz_), nq_) << "\n";
  g << g.copy("arg[" + str(INTEGRATOR_P) + "]", np_, "pv") << "\n";

  g.comment("Forward integration: a run of steps per control group, as RuntimeIntegrator::eval");
  std::string tout = g.constant(tout_), gs = g.constant(group_start_), t0 = g.constant(t0_);
  g << "for (j = 0; j < " << ngroup() << "; ++j) {\n";
  if (nu_ > 0) {
    g << g.copy("arg[" + str(INTEGRATOR_U) + "] ? arg[" + str(INTEGRATOR_U) + "] + j*"
      + str(nu_) + " : 0", nu_, "uv") << "\n";
  }
  g << "casadi_ivp_reset(&d, j ? " << tout << "[" << gs << "[j]-1] : " << t0 << ", "
    << (nu_ > 0 ? "uv" : "0") << ", j == 0)" << nl;
  g << "for (k = " << gs << "[j]; k < " << gs << "[j+1]; ++k) {\n";
  g << "casadi_ivp_advance(&d, " << tout << "[k], " << tout << "[" << gs << "[j+1]-1], "
    << (nrx_ > 0 ? 1 : 0) << ")" << nl;
  codegen_loop(g);
  g << "if (d.status != IVP_SUCCESS) return 1" << nl;
  g << "if (res[" << INTEGRATOR_XF << "]) " << g.copy("d.y_out", nx_, "res["
    + str(INTEGRATOR_XF) + "] + k*" + str(nx_)) << "\n";
  g << "if (res[" << INTEGRATOR_ZF << "]) " << g.copy("d.y_out+" + str(nx_), nz_, "res["
    + str(INTEGRATOR_ZF) + "] + k*" + str(nz_)) << "\n";
  g << "if (res[" << INTEGRATOR_QF << "]) " << g.copy("d.y_out+" + str(nx_ + nz_), nq_,
    "res[" + str(INTEGRATOR_QF) + "] + k*" + str(nq_)) << "\n";
  g << "}\n}\n";

  if (nrx_ > 0) {
    g.comment("Discrete adjoint, output interval by output interval");
    g << "casadi_ivp_reset_backward(&d)" << nl;
    g << "for (j = " << ngroup() - 1 << "; j >= 0; --j) {\n";
    g.comment("Control adjoints of this group only");
    g << g.clear("d.adj_u", nuq_) << "\n";
    g << "for (k = " << gs << "[j+1]-1; k >= " << gs << "[j]; --k) {\n";
    g.comment("Impulse from the adjoint seeds at output time k");
    const casadi_int adj_in[3] = {INTEGRATOR_ADJ_XF, INTEGRATOR_ADJ_ZF, INTEGRATOR_ADJ_QF};
    const casadi_int adj_n[3] = {nrx_, nrz_, nrp_};
    const casadi_int adj_off[3] = {0, nrx_, nrx_ + nrz_};
    for (casadi_int i = 0; i < 3; ++i) {
      if (adj_n[i] == 0) continue;
      g << g.copy("arg[" + str(adj_in[i]) + "] ? arg[" + str(adj_in[i]) + "] + k*"
        + str(adj_n[i]) + " : 0", adj_n[i], "d.seed+" + str(adj_off[i])) << "\n";
    }
    g << "casadi_ivp_impulse(&d, " << tout << "[k])" << nl;
    codegen_loop(g);
    g << "if (d.status != IVP_SUCCESS) return 1" << nl;
    g << "casadi_ivp_retreat(&d, " << tout << "[k], k > 0 ? " << tout << "[k-1] : " << t0
      << ")" << nl;
    codegen_loop(g);
    g << "if (d.status != IVP_SUCCESS) return 1" << nl;
    g << "}\n";
    g << "if (res[" << INTEGRATOR_ADJ_U << "]) " << g.copy("d.adj_u", nuq_, "res["
      + str(INTEGRATOR_ADJ_U) + "] + j*" + str(nuq_)) << "\n";
    g << "}\n";
    g << "if (res[" << INTEGRATOR_ADJ_X0 << "]) " << g.copy("d.ry", nrx_, "res["
      + str(INTEGRATOR_ADJ_X0) + "]") << "\n";
    g << "if (res[" << INTEGRATOR_ADJ_P << "]) " << g.copy("d.adj_p", nrq_, "res["
      + str(INTEGRATOR_ADJ_P) + "]") << "\n";
  }
}

void RuntimeIntegrator::codegen_loop(CodeGenerator& g) const {
  // Carry out the tasks requested by the kernel, for it or for the nominal integration
  g << "while (casadi_" << plugin_name() << "(&d, &e)) {\n";
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

void RuntimeIntegrator::codegen_ivp_prob(CodeGenerator& g, const std::string& p,
    casadi_int tape_size) const {
  // Shared problem structure in variable p, as set_ivp_prob
  auto set = [&](const std::string& field, const std::string& value) -> void {
    g << p << "." << field << " = " << value << ";\n";
  };
  set("nx1", str(nx1_));
  set("nz1", str(nz1_));
  set("nq1", str(nq1_));
  set("np1", str(np1_));
  set("nu1", str(nu1_));
  set("nfwd", str(nfwd_));
  set("nadj", str(nadj_));
  set("n_stages", str(alpha_.size()));
  set("a", g.constant(a_));
  set("alpha", g.constant(alpha_));
  set("m", g.constant(m_));
  set("n_dense", str(n_dense_));
  set("m_dense", g.constant(m_dense_));
  for (casadi_int k = 0; k < IVP_NUM_JAC; ++k) {
    set("sp_jac[" + str(k) + "]", sp_jac_.empty() ? "0" : g.sparsity(sp_jac_[k]));
  }
  set("tape_size", str(tape_size));
  set("step_size", g.constant(p_.step_size));
  set("nx_err", str(p_.nx_err));
  set("nz_err", str(p_.nz_err));
  set("nq_err", str(p_.nq_err));
  set("abstol", g.constant(p_.abstol));
  set("reltol", g.constant(p_.reltol));
  set("order_err", str(p_.order_err));
  set("m_err", m_err_.empty() ? "0" : g.constant(m_err_));
  set("safety", g.constant(p_.safety));
  set("fac_min", g.constant(p_.fac_min));
  set("fac_max", g.constant(p_.fac_max));
  set("h_min", g.constant(p_.h_min));
  set("h_max", g.constant(p_.h_max));
  set("h_init", g.constant(p_.h_init));
  set("max_num_steps", str(p_.max_num_steps));
  set("ne", str(p_.ne));
  if (p_.ne > 0) {
    set("has_transition", str(p_.has_transition));
    set("event_tol", g.constant(p_.event_tol));
    set("event_acceptable_tol", g.constant(p_.event_acceptable_tol));
    set("max_event_iter", str(p_.max_event_iter));
    set("max_events", str(p_.max_events));
  }
  g << "casadi_ivp_setup(&" << p << ");\n";
}

void RuntimeIntegrator::codegen_task(CodeGenerator& g, const std::string& d,
    const std::vector<std::vector<TaskCall>>& calls, bool in_event) const {
  // Straight-line code for the task table, as eval_task: direct calls, no function pointers
  static const char* task_name[] = {"IVP_EVAL_RHS", "IVP_EVAL_JAC", "IVP_EVAL_VJP",
    "IVP_EVAL_JAC_ADJ", "IVP_EVAL_ZERO", "IVP_EVAL_JUMP", "IVP_EVAL_JUMP_ADJ"};
  g << "switch (" << d << ".task) {\n";
  for (casadi_int task = 0; task < static_cast<casadi_int>(calls.size()); ++task) {
    if (calls[task].empty()) continue;
    // Tasks of the events are requested by casadi_ivp_event and casadi_ivp_event_adj only
    bool ev = task >= IVP_EVAL_ZERO && task <= IVP_EVAL_JUMP_ADJ;
    if (in_event != ev) continue;
    g << "case " << task_name[task] << ":\n";
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
  if (step_size_ == 0 && !in_event) {
    g << "case IVP_STEP_CONTROL:\n";
    g << "casadi_ivp_control(&" << d << ");\n";
    g << "break;\n";
  }
  if (ne_ > 0 && !in_event) {
    for (casadi_int adj = 0; adj <= (nadj_ > 0 ? 1 : 0); ++adj) {
      g << "case IVP_EVENT" << (adj ? "_ADJ" : "") << ":\n";
      g << "while (casadi_ivp_event" << (adj ? "_adj" : "") << "(&" << d << ")) {\n";
      codegen_task(g, d, calls, true);
      g << "}\n";
      g << "break;\n";
    }
  }
  g << "default: " << d << ".fail = 1;\n";
  g << "}\n";
}

Dict RuntimeIntegrator::get_stats(void* mem) const {
  Dict stats = Integrator::get_stats(mem);
  auto m = static_cast<RuntimeIntegratorMemory*>(mem);
  stats["n_steps"] = m->d.n_steps;
  stats["n_reject"] = m->d.n_reject;
  stats["n_rhs"] = m->d.n_rhs;
  stats["n_jac"] = m->d.n_jac;
  stats["n_vjp"] = m->d.n_vjp;
  stats["n_events"] = m->d.n_events;
  return stats;
}

void RuntimeIntegrator::print_stats(IntegratorMemory* mem) const {
  auto m = static_cast<RuntimeIntegratorMemory*>(mem);
  print("%s: %lld steps, %lld rejected, %lld rhs, %lld jac, %lld vjp, %lld events\n",
    plugin_name(), static_cast<long long>(m->d.n_steps), static_cast<long long>(m->d.n_reject),
    static_cast<long long>(m->d.n_rhs), static_cast<long long>(m->d.n_jac),
    static_cast<long long>(m->d.n_vjp), static_cast<long long>(m->d.n_events));
}

RuntimeIntegrator::RuntimeIntegrator(DeserializingStream& s) : Integrator(s) {
  s.version("RuntimeIntegrator", 1);
  s.unpack("RuntimeIntegrator::group_start", group_start_);
  s.unpack("RuntimeIntegrator::nx_err", nx_err_);
  s.unpack("RuntimeIntegrator::nz_err", nz_err_);
  s.unpack("RuntimeIntegrator::nq_err", nq_err_);
  s.unpack("RuntimeIntegrator::sp_jac", sp_jac_);
  s.unpack("RuntimeIntegrator::step_size", step_size_);
  s.unpack("RuntimeIntegrator::abstol", abstol_);
  s.unpack("RuntimeIntegrator::reltol", reltol_);
  s.unpack("RuntimeIntegrator::step0", step0_);
  s.unpack("RuntimeIntegrator::max_step_size", max_step_size_);
  s.unpack("RuntimeIntegrator::min_step_size", min_step_size_);
  s.unpack("RuntimeIntegrator::max_num_steps", max_num_steps_);
  s.unpack("RuntimeIntegrator::max_tape_steps", max_tape_steps_);
  s.unpack("RuntimeIntegrator::quad_err_con", quad_err_con_);
}

void RuntimeIntegrator::serialize_body(SerializingStream &s) const {
  Integrator::serialize_body(s);
  s.version("RuntimeIntegrator", 1);
  s.pack("RuntimeIntegrator::group_start", group_start_);
  s.pack("RuntimeIntegrator::nx_err", nx_err_);
  s.pack("RuntimeIntegrator::nz_err", nz_err_);
  s.pack("RuntimeIntegrator::nq_err", nq_err_);
  s.pack("RuntimeIntegrator::sp_jac", sp_jac_);
  s.pack("RuntimeIntegrator::step_size", step_size_);
  s.pack("RuntimeIntegrator::abstol", abstol_);
  s.pack("RuntimeIntegrator::reltol", reltol_);
  s.pack("RuntimeIntegrator::step0", step0_);
  s.pack("RuntimeIntegrator::max_step_size", max_step_size_);
  s.pack("RuntimeIntegrator::min_step_size", min_step_size_);
  s.pack("RuntimeIntegrator::max_num_steps", max_num_steps_);
  s.pack("RuntimeIntegrator::max_tape_steps", max_tape_steps_);
  s.pack("RuntimeIntegrator::quad_err_con", quad_err_con_);
}

} // namespace casadi

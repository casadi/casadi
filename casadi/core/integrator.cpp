/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A symbolic framework for dynamic optimization.
 *    Copyright (C) 2010-2023 Joel Andersson, Joris Gillis, Moritz Diehl,
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


#include "integrator_impl.hpp"
#include "casadi_misc.hpp"

namespace casadi {

std::string to_string(DynIn v) {
  switch (v) {
  case DYN_T: return "t";
  case DYN_X: return "x";
  case DYN_Z: return "z";
  case DYN_P: return "p";
  case DYN_U: return "u";
  default: break;
  }
  return "";
}

std::string to_string(DynOut v) {
  switch (v) {
  case DYN_ODE: return "ode";
  case DYN_ALG: return "alg";
  case DYN_QUAD: return "quad";
  case DYN_ZERO: return "zero";
  default: break;
  }
  return "";
}

std::string to_string(EventIn v) {
  switch (v) {
  case EVENT_INDEX: return "index";
  case EVENT_T: return "t";
  case EVENT_X: return "x";
  case EVENT_Z: return "z";
  case EVENT_P: return "p";
  case EVENT_U: return "u";
  default: break;
  }
  return "";
}

std::string to_string(EventOut v) {
  switch (v) {
  case EVENT_POST_X: return "post_x";
  case EVENT_POST_Z: return "post_z";
  default: break;
  }
  return "";
}

std::string Integrator::bdyn_in(casadi_int i) {
  switch (i) {
    case BDYN_T: return "t";
    case BDYN_X: return "x";
    case BDYN_Z: return "z";
    case BDYN_P: return "p";
    case BDYN_U: return "u";
    case BDYN_OUT_ODE: return "out_ode";
    case BDYN_OUT_ALG: return "out_alg";
    case BDYN_OUT_QUAD: return "out_quad";
    case BDYN_OUT_ZERO: return "out_zero";
    case BDYN_ADJ_ODE: return "adj_ode";
    case BDYN_ADJ_ALG: return "adj_alg";
    case BDYN_ADJ_QUAD: return "adj_quad";
    case BDYN_ADJ_ZERO: return "adj_zero";
  default: break;
  }
  return "";
}

std::vector<std::string> Integrator::bdyn_in() {
  std::vector<std::string> ret(BDYN_NUM_IN);
  for (casadi_int i = 0; i < BDYN_NUM_IN; ++i)
    ret[i] = bdyn_in(i);
  return ret;
}

std::string Integrator::bdyn_out(casadi_int i) {
  switch (i) {
    case BDYN_ADJ_T: return "adj_t";
    case BDYN_ADJ_X: return "adj_x";
    case BDYN_ADJ_Z: return "adj_z";
    case BDYN_ADJ_P: return "adj_p";
    case BDYN_ADJ_U: return "adj_u";
  default: break;
  }
  return "";
}

std::vector<std::string> Integrator::bdyn_out() {
  std::vector<std::string> ret(BDYN_NUM_OUT);
  for (casadi_int i = 0; i < BDYN_NUM_OUT; ++i)
    ret[i] = bdyn_out(i);
  return ret;
}

bool has_integrator(const std::string& name) {
  return Integrator::has_plugin(name);
}

void load_integrator(const std::string& name) {
  Integrator::load_plugin(name);
}

std::string doc_integrator(const std::string& name) {
  return Integrator::getPlugin(name).doc;
}

Function integrator(const std::string& name, const std::string& solver,
    const SXDict& dae, const Dict& opts) {
  return integrator(name, solver, dae, 0.0, std::vector<double>{1.0}, opts);
}

Function integrator(const std::string& name, const std::string& solver,
    const MXDict& dae, const Dict& opts) {
  return integrator(name, solver, dae, 0.0, std::vector<double>{1.0}, opts);
}

Function integrator(const std::string& name, const std::string& solver,
    const Function& dae, const Dict& opts) {
  return integrator(name, solver, dae, 0.0, std::vector<double>{1.0}, opts);
}

Function integrator(const std::string& name, const std::string& solver,
    const SXDict& dae, double t0, const std::vector<double>& tout, const Dict& opts) {
  return integrator(name, solver, Integrator::map2oracle("dae", dae), t0, tout, opts);
}

Function integrator(const std::string& name, const std::string& solver,
    const MXDict& dae, double t0, const std::vector<double>& tout, const Dict& opts) {
  return integrator(name, solver, Integrator::map2oracle("dae", dae), t0, tout, opts);
}

Function integrator(const std::string& name, const std::string& solver,
    const Function& dae, double t0, const std::vector<double>& tout, const Dict& opts) {
  // Make sure that dae is sound
  if (dae.has_free()) {
    casadi_error("Cannot create '" + name + "' since " + str(dae.get_free()) + " are free.");
  }
  Integrator* intg = Integrator::getPlugin(solver).creator(name, dae, t0, tout);
  return intg->create_advanced(opts);
}

Function integrator(const std::string& name, const std::string& solver,
    const SXDict& dae, double t0, double tf, const Dict& opts) {
  return integrator(name, solver, dae, t0, std::vector<double>{tf}, opts);
}

Function integrator(const std::string& name, const std::string& solver,
    const MXDict& dae, double t0, double tf, const Dict& opts) {
  return integrator(name, solver, dae, t0, std::vector<double>{tf}, opts);
}

Function integrator(const std::string& name, const std::string& solver,
    const Function& dae, double t0, double tf, const Dict& opts) {
  return integrator(name, solver, dae, t0, std::vector<double>{tf}, opts);
}

std::vector<std::string> integrator_in() {
  std::vector<std::string> ret(integrator_n_in());
  for (size_t i=0; i<ret.size(); ++i) ret[i]=integrator_in(i);
  return ret;
}

std::vector<std::string> integrator_out() {
  std::vector<std::string> ret(integrator_n_out());
  for (size_t i=0; i<ret.size(); ++i) ret[i]=integrator_out(i);
  return ret;
}

std::string integrator_in(casadi_int ind) {
  switch (static_cast<IntegratorInput>(ind)) {
  case INTEGRATOR_X0:  return "x0";
  case INTEGRATOR_Z0:  return "z0";
  case INTEGRATOR_P:   return "p";
  case INTEGRATOR_U:   return "u";
  case INTEGRATOR_ADJ_XF: return "adj_xf";
  case INTEGRATOR_ADJ_ZF: return "adj_zf";
  case INTEGRATOR_ADJ_QF:  return "adj_qf";
  case INTEGRATOR_NUM_IN: break;
  }
  return std::string();
}

std::string integrator_out(casadi_int ind) {
  switch (static_cast<IntegratorOutput>(ind)) {
  case INTEGRATOR_XF:  return "xf";
  case INTEGRATOR_ZF:  return "zf";
  case INTEGRATOR_QF:  return "qf";
  case INTEGRATOR_ADJ_X0: return "adj_x0";
  case INTEGRATOR_ADJ_Z0: return "adj_z0";
  case INTEGRATOR_ADJ_P: return "adj_p";
  case INTEGRATOR_ADJ_U: return "adj_u";
  case INTEGRATOR_NUM_OUT: break;
  }
  return std::string();
}

casadi_int integrator_n_in() {
  return INTEGRATOR_NUM_IN;
}

casadi_int integrator_n_out() {
  return INTEGRATOR_NUM_OUT;
}

std::vector<std::string> dyn_in() {
  return enum_names<DynIn>();
}

std::vector<std::string> dyn_out() {
  return enum_names<DynOut>();
}

std::string dyn_in(casadi_int ind) {
  return to_string(static_cast<DynIn>(ind));
}

std::string dyn_out(casadi_int ind) {
  return to_string(static_cast<DynOut>(ind));
}

casadi_int dyn_n_in() {
  return DYN_NUM_IN;
}

casadi_int dyn_n_out() {
  return DYN_NUM_OUT;
}

Integrator::Integrator(const std::string& name, const Function& oracle,
    double t0, const std::vector<double>& tout)
    : OracleFunction(name, oracle), t0_(t0), tout_(tout) {

  // Negative number of parameters for consistancy checking
  np_ = -1;

  // Default options
  nfwd_ = 0;
  nadj_ = 0;
  print_stats_ = false;
}

Integrator::~Integrator() {
}

Sparsity Integrator::get_sparsity_in(casadi_int i) {
  switch (static_cast<IntegratorInput>(i)) {
  case INTEGRATOR_X0: return Sparsity::dense(nx1_, 1 + nfwd_);
  case INTEGRATOR_Z0: return Sparsity::dense(nz1_, 1 + nfwd_);
  case INTEGRATOR_P: return Sparsity::dense(np1_, 1 + nfwd_);
  case INTEGRATOR_U: return Sparsity::dense(nu1_, nt() * (1 + nfwd_));
  case INTEGRATOR_ADJ_XF: return Sparsity::dense(nrx1_, nadj_ * (1 + nfwd_) * nt());
  case INTEGRATOR_ADJ_ZF: return Sparsity::dense(nrz1_, nadj_ * (1 + nfwd_) * nt());
  case INTEGRATOR_ADJ_QF: return Sparsity::dense(nrp1_, nadj_ * (1 + nfwd_) * nt());
  case INTEGRATOR_NUM_IN: break;
  }
  return Sparsity();
}

Sparsity Integrator::get_sparsity_out(casadi_int i) {
  switch (static_cast<IntegratorOutput>(i)) {
  case INTEGRATOR_XF: return Sparsity::dense(nx1_, nt() * (1 + nfwd_));
  case INTEGRATOR_ZF: return Sparsity::dense(nz1_, nt() * (1 + nfwd_));
  case INTEGRATOR_QF: return Sparsity::dense(nq1_, nt() * (1 + nfwd_));
  case INTEGRATOR_ADJ_X0: return Sparsity::dense(nrx1_, nadj_ * (1 + nfwd_));
  case INTEGRATOR_ADJ_Z0: return Sparsity(nrz1_, nadj_ * (1 + nfwd_));  // always zero
  case INTEGRATOR_ADJ_P: return Sparsity::dense(nrq1_, nadj_ * (1 + nfwd_));
  case INTEGRATOR_ADJ_U: return Sparsity::dense(nuq1_, nadj_ * (1 + nfwd_) * nt());
  case INTEGRATOR_NUM_OUT: break;
  }
  return Sparsity();
}

bool Integrator::grid_in(casadi_int i) {
  switch (static_cast<IntegratorInput>(i)) {
    case INTEGRATOR_U:
    case INTEGRATOR_ADJ_XF:
    case INTEGRATOR_ADJ_ZF:
    case INTEGRATOR_ADJ_QF:
      return true;
    default: break;
  }
  return false;
}

bool Integrator::grid_out(casadi_int i) {
  switch (static_cast<IntegratorOutput>(i)) {
    case INTEGRATOR_XF:
    case INTEGRATOR_ZF:
    case INTEGRATOR_QF:
    case INTEGRATOR_ADJ_U:
      return true;
    default: break;
  }
  return false;
}

casadi_int Integrator::adjmap_out(casadi_int i) {
  switch (static_cast<IntegratorInput>(i)) {
    case INTEGRATOR_X0: return INTEGRATOR_ADJ_X0;
    case INTEGRATOR_Z0: return INTEGRATOR_ADJ_Z0;
    case INTEGRATOR_P: return INTEGRATOR_ADJ_P;
    case INTEGRATOR_U: return INTEGRATOR_ADJ_U;
    case INTEGRATOR_ADJ_XF: return INTEGRATOR_XF;
    case INTEGRATOR_ADJ_ZF: return INTEGRATOR_ZF;
    case INTEGRATOR_ADJ_QF: return INTEGRATOR_QF;
    default: break;
  }
  return -1;
}

Function Integrator::create_advanced(const Dict& opts) {
  return Function::create(this, opts);
}

int Integrator::eval(const double** arg, double** res,
    casadi_int* iw, double* w, void* mem) const {
  auto m = static_cast<IntegratorMemory*>(mem);

  // Read inputs
  const double* x0 = arg[INTEGRATOR_X0];
  const double* z0 = arg[INTEGRATOR_Z0];
  const double* p = arg[INTEGRATOR_P];
  const double* u = arg[INTEGRATOR_U];
  const double* adj_xf = arg[INTEGRATOR_ADJ_XF];
  const double* rz0 = arg[INTEGRATOR_ADJ_ZF];
  const double* rp = arg[INTEGRATOR_ADJ_QF];
  arg += INTEGRATOR_NUM_IN;

  // Read outputs
  double* x = res[INTEGRATOR_XF];
  double* z = res[INTEGRATOR_ZF];
  double* q = res[INTEGRATOR_QF];
  double* rx = res[INTEGRATOR_ADJ_X0];
  double* rq = res[INTEGRATOR_ADJ_P];
  double* adj_u = res[INTEGRATOR_ADJ_U];
  res += INTEGRATOR_NUM_OUT;

  // Setup memory object
  setup(m, arg, res, iw, w);

  // Reset solver, take time to t0
  m->t = t0_;
  reset(m, u, x0, z0, p);

  // Next stop time due to step change in input
  casadi_int k_stop = next_stop(0, u);

  // Integrate forward
  for (m->k = 0; m->k < nt(); ++m->k) {
    // Update stopping time, if needed
    if (m->k > k_stop) k_stop = next_stop(m->k, u);
    // Advance solution
    m->t_next = tout_[m->k];
    m->t_stop = tout_[k_stop];
    if (verbose_) casadi_message("Integrating forward to output time " + str(m->k) + ": t_next = "
      + str(m->t_next) + ", t_stop = " + str(m->t_stop));
    advance(m, u, x, z, q);
    if (x) x += nx_;
    if (z) z += nz_;
    if (q) q += nq_;
    if (u) u += nu_;
    m->t = m->t_next;
  }

  // Backwards integration, if needed
  if (nrx_ > 0) {
    // Take adj_xf, rz0, rp past the last grid point
    if (adj_xf) adj_xf += nrx_ * nt();
    if (rz0) rz0 += nrz_ * nt();
    if (rp) rp += nrp_ * nt();
    if (adj_u) adj_u += nuq_ * nt();
    // Next stop time due to step change in input
    k_stop = nt();
    // Reset the solver
    resetB(m);
    // Any adjoint seed so far?
    bool any_impulse = false;
    // Integrate backward
    for (m->k = nt(); m->k-- > 0; ) {
      m->t = tout_[m->k];
      // Add impulse to backwards integration
      if (adj_xf) adj_xf -= nrx_;
      if (rz0) rz0 -= nrz_;
      if (rp) rp -= nrp_;
      if (adj_u) adj_u -= nuq_;
      if (u) u -= nu_;
      if (!all_zero(adj_xf, nrx_) || !all_zero(rz0, nrz_) || !all_zero(rp, nrp_)) {
        if (verbose_) casadi_message("Impulse from adjoint seeds at output time " + str(m->k));
        impulseB(m, adj_xf, rz0, rp);
        any_impulse = true;
      }
      // Next output time, or beginning
      casadi_int k_next = m->k - 1;
      m->t_next = k_next < 0 ? t0_ : tout_[k_next];
      // Update integrator stopping time
      if (k_next < k_stop) k_stop = next_stopB(m->k, u);
      m->t_stop = k_stop < 0 ? t0_ : tout_[k_stop];
      // Proceed to the previous time point or t0
      if (any_impulse) {
        if (verbose_) casadi_message("Integrating backward from output time " + str(m->k)
          + ": t_next = " + str(m->t_next) + ", t_stop = " + str(m->t_stop));
        if (m->k > 0) {
          retreat(m, u, 0, 0, adj_u);
        } else {
          retreat(m, u, rx, rq, adj_u);
        }
      } else {
        if (verbose_) casadi_message("No adjoint seeds from output time " + str(m->k)
          + ": t_next = " + str(m->t_next) + ", t_stop = " + str(m->t_stop));
        casadi_clear(adj_u, nuq_);
        if (m->k == 0) {
          casadi_clear(rx, nrx_);
          casadi_clear(rq, nrq_);
        }
      }
    }
    // adj_u should contain the contribution from the grid point, not cumulative
    if (adj_u) {
      for (m->k = 0; m->k < nt() - 1; ++m->k) {
        casadi_axpy(nuq_, -1., adj_u + nuq_, adj_u);
        adj_u += nuq_;
      }
    }
  }

  // Collect oracle statistics
  join_results(m);

  // Print integrator statistics
  if (print_stats_) print_stats(m);

  return 0;
}

const Options Integrator::options_
= {{&OracleFunction::options_},
    {{"expand",
      {OT_BOOL,
      "Replace MX with SX expressions in problem formulation [false]"}},
    {"print_stats",
      {OT_BOOL,
      "Print out statistics after integration"}},
    {"nfwd",
     {OT_INT,
      "Number of forward sensitivities to be calculated [0]"}},
    {"nadj",
     {OT_INT,
      "Number of adjoint sensitivities to be calculated [0]"}},
    {"t0",
      {OT_DOUBLE,
      "[DEPRECATED] Beginning of the time horizon"}},
    {"tf",
      {OT_DOUBLE,
      "[DEPRECATED] End of the time horizon"}},
    {"grid",
      {OT_DOUBLEVECTOR,
      "[DEPRECATED] Time grid"}},
    {"augmented_options",
      {OT_DICT,
      "Options to be passed down to the augmented integrator, if one is constructed."}},
    {"output_t0",
      {OT_BOOL,
      "[DEPRECATED] Output the state at the initial time"}}
    }
};

void Integrator::init(const Dict& opts) {
  // Default (temporary) options
  double t0 = 0, tf = 1;
  bool expand = false;
  bool output_t0 = false;
  std::vector<double> grid;
  bool uses_legacy_options = false;

  // Read options
  for (auto&& op : opts) {
    if (op.first=="expand") {
      expand = op.second;
    } else if (op.first=="output_t0") {
      output_t0 = op.second;
      uses_legacy_options = true;
    } else if (op.first=="print_stats") {
      print_stats_ = op.second;
    } else if (op.first=="nfwd") {
      nfwd_ = op.second;
    } else if (op.first=="nadj") {
      nadj_ = op.second;
    } else if (op.first=="grid") {
      grid = op.second;
      uses_legacy_options = true;
    } else if (op.first=="augmented_options") {
      augmented_options_ = op.second;
    } else if (op.first=="t0") {
      t0 = op.second;
      uses_legacy_options = true;
    } else if (op.first=="tf") {
      tf = op.second;
      uses_legacy_options = true;
    }
  }

  // Replace MX oracle with SX oracle?
  if (expand) this->expand();

  // Store a copy of the options, for creating augmented integrators
  opts_ = opts;

  // Construct t0_ and tout_ gbased on legacy options
  if (uses_legacy_options) {
    static bool first_encounter = true;
    if (first_encounter) {
    // Deprecation warning
    casadi_warning("The options 't0', 'tf', 'grid' and 'output_t0' have been deprecated.\n"
      "The same functionality is provided by providing additional input arguments to "
      "the 'integrator' function, in particular:\n"
      " * Call integrator(..., t0, tf, options) for a single output time, or\n"
      " * Call integrator(..., t0, grid, options) for multiple grid points.\n"
      "The legacy 'output_t0' option can be emulated by including or excluding 't0' in 'grid'.\n"
      "Backwards compatibility is provided in this release only.");
      first_encounter = false;
    }

    // If grid unset, default to [t0, tf]
    if (grid.empty()) grid = {t0, tf};

    // Construct t0 and tout from grid and output_t0
    t0_ = grid.front();
    tout_ = grid;
    if (!output_t0) tout_.erase(tout_.begin());
  }

  // Consistency checks: Sensitivities
  casadi_assert(nfwd_ >= 0, "Number of forward sensitivities must be non-negative");
  casadi_assert(nadj_ >= 0, "Number of adjoint sensitivities must be non-negative");

  // Consistency check: Valid oracle
  casadi_assert(oracle_.n_in() == DYN_NUM_IN, "DAE has wrong number of inputs");
  casadi_assert(oracle_.n_out() == DYN_NUM_OUT, "DAE has wrong number of outputs");

  // Consistency checks, input sparsities
  for (casadi_int i = 0; i < DYN_NUM_IN; ++i) {
    const Sparsity& sp = oracle_.sparsity_in(i);
    if (i == DYN_T) {
      casadi_assert(sp.is_empty() || sp.is_scalar(), "DAE time variable must be empty or scalar. "
        "Got dimension " + str(sp.size()));
    } else {
      casadi_assert(sp.is_vector(), "DAE inputs must be vectors. "
        + dyn_in(i) + " has dimension " + str(sp.size()) + ".");
    }
    casadi_assert(sp.is_dense(), "DAE inputs must be dense . "
      + dyn_in(i) + " is sparse.");
  }

  // Consistency checks, output sparsities
  for (casadi_int i = 0; i < DYN_NUM_OUT; ++i) {
    const Sparsity& sp = oracle_.sparsity_out(i);
    casadi_assert(sp.is_vector(), "DAE outputs must be vectors. "
      + dyn_out(i) + " has dimension " + str(sp.size()));
    casadi_assert(sp.is_dense(), "DAE outputs must be dense . "
      + dyn_out(i) + " is sparse.");
  }

  // Get dimensions (excluding sensitivity equations), forward problem
  nx1_ = oracle_.numel_in(DYN_X);
  nz1_ = oracle_.numel_in(DYN_Z);
  nq1_ = oracle_.numel_out(DYN_QUAD);
  np1_ = oracle_.numel_in(DYN_P);
  nu1_ = oracle_.numel_in(DYN_U);
  ne_ = oracle_.numel_out(DYN_ZERO);

  // Event support not implemented
  if (ne_ > 0) {
    casadi_warning("Event support has not yet been implemented");
  }

  // Consistency checks
  casadi_assert(nx1_ > 0, "Ill-posed ODE - no state");
  casadi_assert(nx1_ == oracle_.numel_out(DYN_ODE), "Dimension mismatch for 'ode'");
  casadi_assert(nz1_ == oracle_.numel_out(DYN_ALG), "Dimension mismatch for 'alg'");

  // Backward problem, if any
  if (nadj_ > 0) {
    // Generate backward DAE
    rdae_ = oracle_.reverse(nadj_);
    // Consistency checks
    casadi_assert(rdae_.n_in() == BDYN_NUM_IN, "Backward DAE has wrong number of inputs");
    casadi_assert(rdae_.n_out() == BDYN_NUM_OUT, "Backward DAE has wrong number of outputs");
    casadi_assert(rdae_.numel_in(BDYN_X) == nx1_, "Dimension mismatch");
    casadi_assert(rdae_.numel_in(BDYN_Z) == nz1_, "Dimension mismatch");
    casadi_assert(rdae_.numel_in(BDYN_P) == np1_, "Dimension mismatch");
    casadi_assert(rdae_.numel_in(BDYN_U) == nu1_, "Dimension mismatch");
    casadi_assert(rdae_.numel_in(BDYN_ADJ_ODE) == nx1_ * nadj_, "Inconsistent dimensions");
    casadi_assert(rdae_.numel_in(BDYN_ADJ_ALG) == nz1_ * nadj_, "Inconsistent dimensions");
    casadi_assert(rdae_.numel_in(BDYN_ADJ_QUAD) == nq1_ * nadj_, "Inconsistent dimensions");
    casadi_assert(rdae_.numel_out(BDYN_ADJ_P) == np1_ * nadj_, "Inconsistent dimensions");
    casadi_assert(rdae_.numel_out(BDYN_ADJ_U) == nu1_ * nadj_, "Inconsistent dimensions");

    // Dimensions (excluding sensitivity equations), backward problem
    nrx1_ = nx1_;
    nrz1_ = nz1_;
    nrp1_ = nq1_;
    nrq1_ = np1_;
    nuq1_ = nu1_;
  } else {
    // No backward problem
    nrx1_ = nrz1_ = nrp1_ = nrq1_ = nuq1_ = 0;
  }

  // Get dimensions (including sensitivity equations)
  nx_ = nx1_ * (1 + nfwd_);
  nz_ = nz1_ * (1 + nfwd_);
  nq_ = nq1_ * (1 + nfwd_);
  np_ = np1_ * (1 + nfwd_);
  nu_ = nu1_ * (1 + nfwd_);
  nrx_ = nrx1_ * nadj_ * (1 + nfwd_);
  nrz_ = nrz1_ * nadj_ * (1 + nfwd_);
  nrp_ = nrp1_ * nadj_ * (1 + nfwd_);
  nrq_ = nrq1_ * nadj_ * (1 + nfwd_);
  nuq_ = nuq1_ * nadj_ * (1 + nfwd_);

  // Call the base class method
  OracleFunction::init(opts);

  // Create problem functions, forward problem
  create_function("daeF", dyn_in(), dae_out());
  if (nq_ > 0) create_function("quadF", dyn_in(), quad_out());
  if (nfwd_ > 0) {
    // one direction to conserve memory, symbolic processing time
    create_forward("daeF", 1);
    if (nq_ > 0) create_forward("quadF", 1);
  }
  // Zero-crossing function
  if (ne_ > 0) create_function("zero", dyn_in(), zero_out());

  // Create problem functions, backward problem
  if (nadj_ > 0) {
    create_function(rdae_, "daeB", bdyn_in(), bdae_out());
    if (nrq_ > 0 || nuq_ > 0) create_function(rdae_, "quadB", bdyn_in(), bquad_out());
    if (nfwd_ > 0) {
      // one direction to conserve memory, symbolic processing time
      create_forward("daeB", 1);
      if (nrq_ > 0 || nuq_ > 0) create_forward("quadB", 1);
    }
  }

  // Nominal values for states
  nom_x_ = oracle_.nominal_in(DYN_X);
  nom_z_ = oracle_.nominal_in(DYN_Z);

  // Get the sparsities of the forward and reverse DAE
  sp_jac_dae_ = sp_jac_dae();
  casadi_assert(!sp_jac_dae_.is_singular(),
    "Jacobian of the forward problem is structurally rank-deficient. "
    "sprank(J)=" + str(sprank(sp_jac_dae_)) + "<" + str(nx_+nz_));
  if (nadj_ > 0) {
    sp_jac_rdae_ = sp_jac_rdae();
    casadi_assert(!sp_jac_rdae_.is_singular(),
      "Jacobian of the backward problem is structurally rank-deficient. "
      "sprank(J)=" + str(sprank(sp_jac_rdae_)) + "<" + str(nrx_+nrz_));
  }

  // Work vectors for sparsity pattern propagation: Can be reused in derived classes
  alloc_w(nx_ + nz_, true); // x, z
  alloc_w(nx_, true); // x_prev
  alloc_w(nrx_ + nrz_, true); // rx, rz
  alloc_w(nrx_, true); // adj_ode
  alloc_w(nrq_, true); // rq
  alloc_w(nx_+nz_);  // Sparsity::sp_solve
  alloc_w(nrx_+nrz_);  // Sparsity::sp_solve
  alloc_w(np_, true); // p
  alloc_w(nu_, true); // u
  alloc_w(nrp_, true); // rp
}

void Integrator::set_work(void* mem, const double**& arg, double**& res,
    casadi_int*& iw, double*& w) const {
  auto m = static_cast<IntegratorMemory*>(mem);

  // Set work in base classes
  OracleFunction::set_work(mem, arg, res, iw, w);

  // Work vectors
  m->x = w; w += nx_;  // doubles as xz
  m->z = w; w += nz_;
  m->x_prev = w; w += nx_;
  m->rx = w; w += nrx_;  // doubles as xz
  m->rz = w; w += nrz_;
  m->adj_ode = w; w += nrx_;
  m->rq = w; w += nrq_;
  m->p = w; w += np_;
  m->u = w; w += nu_;
  m->rp = w; w += nrp_;
}

int Integrator::init_mem(void* mem) const {
  if (OracleFunction::init_mem(mem)) return 1;

  //auto m = static_cast<IntegratorMemory*>(mem);
  return 0;
}

Function Integrator::augmented_dae() const {
  // If no sensitivities, augmented oracle is the oracle itself
  if (nfwd_ == 0) return oracle_;
  // Name of augmented DAE
  std::string aug_name = "fsens" + str(nfwd_) + "_" + oracle_.name();
  // Use function in cache, if available
  Function ret;
  // if (incache(aug_name, ret)) return ret;  // caching disabled while implementing #3047
  // Create new augmented oracle
  try {
    if (oracle_.is_a("SXFunction")) {
      ret = get_forward_dae<SX>(aug_name);
    } else {
      ret = get_forward_dae<MX>(aug_name);
    }
  } catch (std::exception& e) {
    casadi_error("Failed to generate augmented DAE for " + name_ + ":\n" + e.what());
  }
  // Save to Function cache and return
  // tocache(ret);  // caching disabled while implementing #3047
  return ret;
}

template<typename MatType>
Function Integrator::get_forward_dae(const std::string& name) const {
  if (verbose_) casadi_message(name_ + "::get_forward_dae");

  // Get input and output expressions
  std::vector<MatType> arg = MatType::get_input(oracle_);
  std::vector<MatType> res = oracle_(arg);

  // Symbolic expression for augmented DAE
  std::vector<std::vector<MatType>> aug_in(DYN_NUM_IN);
  for (casadi_int i = 0; i < DYN_NUM_IN; ++i) aug_in[i].push_back(arg.at(i));
  std::vector<std::vector<MatType>> aug_out(DYN_NUM_OUT);
  for (casadi_int i = 0; i < DYN_NUM_OUT; ++i) aug_out[i].push_back(res.at(i));

  // Zero of time dimension
  MatType zero_t = MatType::zeros(oracle_.sparsity_in(DYN_T));

  // Augment aug_in with forward sensitivity seeds
  std::vector<std::vector<MatType>> seed(nfwd_, std::vector<MatType>(DYN_NUM_IN));
  for (casadi_int d = 0; d < nfwd_; ++d) {
    // Create expressions for augmented states
    std::string pref = "aug" + str(d) + "_";
    for (casadi_int i = 0; i < DYN_NUM_IN; ++i) {
      if (i == DYN_T) {
        seed[d][i] = zero_t;
      } else {
        seed[d][i] = MatType::sym(pref + dyn_in(i), oracle_.sparsity_in(i));
      }
    }
    // Save to augmented function inputs
    for (casadi_int i = 0; i < DYN_NUM_IN; ++i) {
      if (i != DYN_T) aug_in[i].push_back(seed[d][i]);
    }
  }

  // Calculate directional derivatives
  std::vector<std::vector<MatType>> sens;
  bool always_inline = oracle_.is_a("SXFunction") || oracle_.is_a("MXFunction");
  oracle_->call_forward(arg, res, seed, sens, always_inline, false);

  // Augment aug_out with forward sensitivity equations
  casadi_assert_dev(sens.size() == nfwd_);
  for (casadi_int d = 0; d < nfwd_; ++d) {
    casadi_assert_dev(sens[d].size() == DYN_NUM_OUT);
    for (casadi_int i = 0; i < DYN_NUM_OUT; ++i) {
      aug_out[i].push_back(project(sens[d][i], oracle_.sparsity_out(i)));
    }
  }

  // Concatenate arrays
  for (casadi_int i = 0; i < DYN_NUM_IN; ++i) arg.at(i) = vertcat(aug_in[i]);
  for (casadi_int i = 0; i < DYN_NUM_OUT; ++i) res.at(i) = vertcat(aug_out[i]);

  // Convert to oracle function and return
  return Function(name, arg, res, dyn_in(), dyn_out());
}

int Integrator::fdae_sp_forward(SpForwardMem* m, const bvec_t* x,
    const bvec_t* p, const bvec_t* u, bvec_t* ode, bvec_t* alg) const {
  // Evaluate nondifferentiated
  m->arg[DYN_T] = nullptr;  // t
  m->arg[DYN_X] = x;  // x
  m->arg[DYN_Z] = nullptr;  // z
  m->arg[DYN_P] = p;  // p
  m->arg[DYN_U] = u;  // u
  m->res[DAE_ODE] = ode;  // ode
  m->res[DAE_ALG] = alg;  // alg
  if (calc_sp_forward("daeF", m->arg, m->res, m->iw, m->w)) return 1;
  // Evaluate sensitivities
  for (casadi_int i = 0; i < nfwd_; ++i) {
    m->arg[DYN_NUM_IN + DAE_ODE] = ode;  // out:ode
    m->arg[DYN_NUM_IN + DAE_ALG] = alg;  // out:alg
    m->arg[DYN_NUM_IN + DAE_NUM_OUT + DYN_T] = nullptr;  // fwd:t
    m->arg[DYN_NUM_IN + DAE_NUM_OUT + DYN_X] = x + (i + 1) * nx1_;  // fwd:x
    m->arg[DYN_NUM_IN + DAE_NUM_OUT + DYN_Z] = nullptr;  // fwd:z
    m->arg[DYN_NUM_IN + DAE_NUM_OUT + DYN_P] = p + (i + 1) * np1_;  // fwd:p
    m->arg[DYN_NUM_IN + DAE_NUM_OUT + DYN_U] = u + (i + 1) * nu1_;  // fwd:u
    m->res[DAE_ODE] = ode + (i + 1) * nx1_;  // fwd:ode
    m->res[DAE_ALG] = alg + (i + 1) * nz1_;  // fwd:alg
    if (calc_sp_forward(forward_name("daeF", 1), m->arg, m->res, m->iw, m->w)) return 1;
  }
  return 0;
}

int Integrator::fquad_sp_forward(SpForwardMem* m, const bvec_t* x, const bvec_t* z,
    const bvec_t* p, const bvec_t* u, bvec_t* quad) const {
  // Evaluate nondifferentiated
  m->arg[DYN_T] = nullptr;  // t
  m->arg[DYN_X] = x;  // x
  m->arg[DYN_Z] = z;  // z
  m->arg[DYN_P] = p;  // p
  m->arg[DYN_U] = u;  // u
  m->res[QUAD_QUAD] = quad;  // quad
  if (calc_sp_forward("quadF", m->arg, m->res, m->iw, m->w)) return 1;
  // Evaluate sensitivities
  for (casadi_int i = 0; i < nfwd_; ++i) {
    m->arg[DYN_NUM_IN + QUAD_QUAD] = quad;  // out:quad
    m->arg[DYN_NUM_IN + QUAD_NUM_OUT + DYN_T] = nullptr;  // fwd:t
    m->arg[DYN_NUM_IN + QUAD_NUM_OUT + DYN_X] = x + (i + 1) * nx1_;  // fwd:x
    m->arg[DYN_NUM_IN + QUAD_NUM_OUT + DYN_Z] = z + (i + 1) * nz1_;  // fwd:z
    m->arg[DYN_NUM_IN + QUAD_NUM_OUT + DYN_P] = p + (i + 1) * np1_;  // fwd:p
    m->arg[DYN_NUM_IN + QUAD_NUM_OUT + DYN_U] = u + (i + 1) * nu1_;  // fwd:u
    m->res[QUAD_QUAD] = quad + (i + 1) * nq1_;  // fwd:quad
    if (calc_sp_forward(forward_name("quadF", 1), m->arg, m->res, m->iw, m->w)) return 1;
  }
  return 0;
}

int Integrator::bdae_sp_forward(SpForwardMem* m, const bvec_t* x, const bvec_t* z,
    const bvec_t* p, const bvec_t* u, const bvec_t* adj_ode, const bvec_t* adj_quad,
    bvec_t* adj_x, bvec_t* adj_z) const {
  // Evaluate nondifferentiated
  m->arg[BDYN_T] = nullptr;  // t
  m->arg[BDYN_X] = x;  // x
  m->arg[BDYN_Z] = z;  // z
  m->arg[BDYN_P] = p;  // p
  m->arg[BDYN_U] = u;  // u
  m->arg[BDYN_OUT_ODE] = nullptr;  // out_ode
  m->arg[BDYN_OUT_ALG] = nullptr;  // out_alg
  m->arg[BDYN_OUT_QUAD] = nullptr;  // out_quad
  m->arg[BDYN_OUT_ZERO] = nullptr;  // out_zero
  m->arg[BDYN_ADJ_ODE] = adj_ode;  // adj_ode
  m->arg[BDYN_ADJ_ALG] = nullptr;  // adj_alg
  m->arg[BDYN_ADJ_QUAD] = adj_quad;  // adj_quad
  m->arg[BDYN_ADJ_ZERO] = nullptr;  // adj_zero
  m->res[BDAE_ADJ_X] = adj_x;  // adj_x
  m->res[BDAE_ADJ_Z] = adj_z;  // adj_z
  if (calc_sp_forward("daeB", m->arg, m->res, m->iw, m->w)) return 1;
  // Evaluate sensitivities
  for (casadi_int i = 0; i < nfwd_; ++i) {
    m->arg[BDYN_NUM_IN + BDAE_ADJ_X] = adj_x;  // out:adj_x
    m->arg[BDYN_NUM_IN + BDAE_ADJ_Z] = adj_z;  // out:adj_z
    m->arg[BDYN_NUM_IN + BDAE_NUM_OUT + BDYN_T] = nullptr;  // fwd:t
    m->arg[BDYN_NUM_IN + BDAE_NUM_OUT + BDYN_X] = x + (i + 1) * nx1_;  // fwd:x
    m->arg[BDYN_NUM_IN + BDAE_NUM_OUT + BDYN_Z] = z + (i + 1) * nz1_;  // fwd:z
    m->arg[BDYN_NUM_IN + BDAE_NUM_OUT + BDYN_P] = p + (i + 1) * np1_;  // fwd:p
    m->arg[BDYN_NUM_IN + BDAE_NUM_OUT + BDYN_U] = u + (i + 1) * nu1_;  // fwd:u
    m->arg[BDYN_NUM_IN + BDAE_NUM_OUT + BDYN_OUT_ODE] = nullptr;  // fwd:out_ode
    m->arg[BDYN_NUM_IN + BDAE_NUM_OUT + BDYN_OUT_ALG] = nullptr;  // fwd:out_alg
    m->arg[BDYN_NUM_IN + BDAE_NUM_OUT + BDYN_OUT_QUAD] = nullptr;  // fwd:out_quad
    m->arg[BDYN_NUM_IN + BDAE_NUM_OUT + BDYN_OUT_ZERO] = nullptr;  // fwd:out_zero
    m->arg[BDYN_NUM_IN + BDAE_NUM_OUT + BDYN_ADJ_ODE]
      = adj_ode + (i + 1) * nrx1_ * nadj_;  // fwd:adj_ode
    m->arg[BDYN_NUM_IN + BDAE_NUM_OUT + BDYN_ADJ_ALG] = nullptr;  // fwd:adj_alg
    m->arg[BDYN_NUM_IN + BDAE_NUM_OUT + BDYN_ADJ_QUAD]
      = adj_quad + (i + 1) * nrz1_ * nadj_;  // fwd:adj_quad
    m->arg[BDYN_NUM_IN + BDAE_NUM_OUT + BDYN_ADJ_ZERO] = nullptr;  // fwd:adj_zero
    m->res[BDAE_ADJ_X] = adj_x + (i + 1) * nrx1_ * nadj_;  // fwd:adj_x
    m->res[BDAE_ADJ_Z] = adj_z + (i + 1) * nrz1_ * nadj_;  // fwd:adj_z
    if (calc_sp_forward(forward_name("daeB", 1), m->arg, m->res, m->iw, m->w)) return 1;
  }
  return 0;
}

int Integrator::bquad_sp_forward(SpForwardMem* m, const bvec_t* x, const bvec_t* z,
    const bvec_t* p, const bvec_t* u, const bvec_t* adj_ode, const bvec_t* rz, const bvec_t* adj_quad,
    bvec_t* adj_p, bvec_t* adj_u) const {
  // Evaluate nondifferentiated
  m->arg[BDYN_T] = nullptr;  // t
  m->arg[BDYN_X] = x;  // x
  m->arg[BDYN_Z] = z;  // z
  m->arg[BDYN_P] = p;  // p
  m->arg[BDYN_U] = u;  // u
  m->arg[BDYN_OUT_ODE] = nullptr;  // out_ode
  m->arg[BDYN_OUT_ALG] = nullptr;  // out_alg
  m->arg[BDYN_OUT_QUAD] = nullptr;  // out_quad
  m->arg[BDYN_OUT_ZERO] = nullptr;  // out_zero
  m->arg[BDYN_ADJ_ODE] = adj_ode;  // adj_ode
  m->arg[BDYN_ADJ_ALG] = rz;  // adj_alg
  m->arg[BDYN_ADJ_QUAD] = adj_quad;  // adj_quad
  m->arg[BDYN_ADJ_ZERO] = nullptr;  // adj_zero
  m->res[BQUAD_ADJ_P] = adj_p;  // adj_p
  m->res[BQUAD_ADJ_U] = adj_u;  // adj_u
  if (calc_sp_forward("quadB", m->arg, m->res, m->iw, m->w)) return 1;
  // Evaluate sensitivities
  for (casadi_int i = 0; i < nfwd_; ++i) {
    m->arg[BDYN_NUM_IN + BQUAD_ADJ_P] = adj_p;  // out:adj_p
    m->arg[BDYN_NUM_IN + BQUAD_ADJ_U] = adj_u;  // out:adj_u
    m->arg[BDYN_NUM_IN + BQUAD_NUM_OUT + BDYN_T] = nullptr;  // fwd:t
    m->arg[BDYN_NUM_IN + BQUAD_NUM_OUT + BDYN_X] = x + (i + 1) * nx1_;  // fwd:x
    m->arg[BDYN_NUM_IN + BQUAD_NUM_OUT + BDYN_Z] = z + (i + 1) * nz1_;  // fwd:z
    m->arg[BDYN_NUM_IN + BQUAD_NUM_OUT + BDYN_P] = p + (i + 1) * np1_;  // fwd:p
    m->arg[BDYN_NUM_IN + BQUAD_NUM_OUT + BDYN_U] = u + (i + 1) * nu1_;  // fwd:u
    m->arg[BDYN_NUM_IN + BQUAD_NUM_OUT + BDYN_OUT_ODE] = nullptr;  // fwd:out_ode
    m->arg[BDYN_NUM_IN + BQUAD_NUM_OUT + BDYN_OUT_ALG] = nullptr;  // fwd:out_alg
    m->arg[BDYN_NUM_IN + BQUAD_NUM_OUT + BDYN_OUT_QUAD] = nullptr;  // fwd:out_quad
    m->arg[BDYN_NUM_IN + BQUAD_NUM_OUT + BDYN_OUT_ZERO] = nullptr;  // fwd:out_zero
    m->arg[BDYN_NUM_IN + BQUAD_NUM_OUT + BDYN_ADJ_ODE] =
      adj_ode + (i + 1) * nrx1_ * nadj_;  // fwd:adj_ode
    m->arg[BDYN_NUM_IN + BQUAD_NUM_OUT + BDYN_ADJ_ALG] =
      rz + (i + 1) * nrz1_ * nadj_;  // fwd:adj_alg
    m->arg[BDYN_NUM_IN + BQUAD_NUM_OUT + BDYN_ADJ_QUAD] =
      adj_quad + (i + 1) * nrp1_ * nadj_;  // fwd:adj_quad
    m->arg[BDYN_NUM_IN + BQUAD_NUM_OUT + BDYN_ADJ_ZERO] = nullptr;  // fwd:adj_zero
    m->res[BQUAD_ADJ_P] = adj_p ? adj_p + (i + 1) * nrq1_ * nadj_ : 0;  // fwd:adj_p
    m->res[BQUAD_ADJ_U] = adj_u ? adj_u + (i + 1) * nuq1_ * nadj_: 0;  // fwd:adj_u
    if (calc_sp_forward(forward_name("quadB", 1), m->arg, m->res, m->iw, m->w)) return 1;
  }
  return 0;
}

int Integrator::sp_forward(const bvec_t** arg, bvec_t** res,
    casadi_int* iw, bvec_t* w, void* mem) const {
  if (verbose_) casadi_message(name_ + "::sp_forward");

  // Inputs
  const bvec_t* x0 = arg[INTEGRATOR_X0];
  const bvec_t* p = arg[INTEGRATOR_P];
  const bvec_t* u = arg[INTEGRATOR_U];
  const bvec_t* adj_xf = arg[INTEGRATOR_ADJ_XF];
  const bvec_t* adj_qf = arg[INTEGRATOR_ADJ_QF];
  arg += n_in_;

  // Outputs
  bvec_t* xf = res[INTEGRATOR_XF];
  bvec_t* zf = res[INTEGRATOR_ZF];
  bvec_t* qf = res[INTEGRATOR_QF];
  bvec_t* adj_x0 = res[INTEGRATOR_ADJ_X0];
  bvec_t* adj_p0 = res[INTEGRATOR_ADJ_P];
  bvec_t* adj_u = res[INTEGRATOR_ADJ_U];
  res += n_out_;

  // Work vectors
  bvec_t *ode = w; w += nx_;
  bvec_t *alg = w; w += nz_;
  bvec_t *x_prev = w; w += nx_;
  bvec_t *rx = w; w += nrx_;
  bvec_t *rz = w; w += nrz_;
  bvec_t *adj_ode = w; w += nrx_;
  bvec_t *rq = w; w += nrq_;

  // Memory struct for function calls below
  SpForwardMem m = {arg, res, iw, w};

  // Copy initial guess to x_prev
  std::copy_n(x0, nx_, x_prev);

  // Propagate forward
  for (casadi_int k = 0; k < nt(); ++k) {
    // Propagate through DAE function
    if (fdae_sp_forward(&m, x_prev, p, u, ode, alg)) return 1;
    for (casadi_int i = 0; i < nx_; ++i) ode[i] |= x_prev[i];

    // "Solve" in order to resolve interdependencies (cf. Rootfinder)
    std::copy_n(ode, nx_, w);
    std::copy_n(alg, nz_, w + nx_);
    std::fill_n(ode, nx_ + nz_, 0);
    sp_jac_dae_.spsolve(ode, w, false);

    // Get xf and zf
    if (xf) std::copy_n(ode, nx_, xf);
    if (zf) std::copy_n(alg, nz_, zf);

    // Propagate to quadratures
    if (nq_ > 0 && qf) {
      if (fquad_sp_forward(&m, ode, alg, p, u, qf)) return 1;
    }

    // Shift time
    std::copy_n(ode, nx_, x_prev);
    if (xf) xf += nx_;
    if (zf) zf += nz_;
    if (qf) qf += nq_;
    if (u) u += nu_;
  }

  if (nrx_ > 0) {
    // Clear adj_ode, adj_p0
    std::fill_n(adj_ode, nrx_, 0);
    if (adj_p0) std::fill_n(adj_p0, nrq_, 0);

    // Take adj_xf, rp, adj_u past the last grid point
    if (adj_xf) adj_xf += nrx_ * nt();
    if (adj_qf) adj_qf += nrp_ * nt();
    if (adj_u) adj_u += nuq_ * nt();

    // Integrate backward
    for (casadi_int k = nt(); k-- > 0; ) {
      // Shift time
      if (adj_xf) adj_xf -= nrx_;
      if (adj_qf) adj_qf -= nrp_;
      if (adj_u) adj_u -= nuq_;
      if (u) u -= nu_;

      // Add impulse from adj_xf
      if (adj_xf) {
        for (casadi_int i = 0; i < nrx_; ++i) adj_ode[i] |= adj_xf[i];
      }

      // Propagate through DAE function
      if (bdae_sp_forward(&m, ode, alg, p, u, adj_ode, adj_qf, rx, rz)) return 1;
      for (casadi_int i = 0; i < nrx_; ++i) rx[i] |= adj_ode[i];

      // "Solve" in order to resolve interdependencies (cf. Rootfinder)
      std::copy_n(rx, nrx_ + nrz_, w);
      std::fill_n(rx, nrx_ + nrz_, 0);
      sp_jac_rdae_.spsolve(rx, w, false);

      // Propagate to quadratures
      if ((nrq_ > 0 && adj_p0) || (nuq_ > 0 && adj_u)) {
        if (bquad_sp_forward(&m, ode, alg, p, u, rx, rz, adj_qf, rq, adj_u)) return 1;
        // Sum contributions to adj_p0
        if (adj_p0) {
          for (casadi_int i = 0; i < nrq_; ++i) adj_p0[i] |= rq[i];
        }
      }

      // Update adj_ode
      std::copy_n(rx, nx_, adj_ode);
    }

    // Get adj_x0 at initial time
    if (adj_x0) std::copy_n(rx, nrx_, adj_x0);
  }
  return 0;
}

int Integrator::fdae_sp_reverse(SpReverseMem* m, bvec_t* x,
    bvec_t* p, bvec_t* u, bvec_t* ode, bvec_t* alg) const {
  // Nondifferentiated inputs
  m->arg[DYN_T] = nullptr;  // t
  m->arg[DYN_X] = x;  // x
  m->arg[DYN_Z] = nullptr;  // z
  m->arg[DYN_P] = p;  // p
  m->arg[DYN_U] = u;  // u
  // Propagate through sensitivities
  for (casadi_int i = 0; i < nfwd_; ++i) {
    m->res[DAE_ODE] = ode + (i + 1) * nx1_;  // fwd:ode
    m->res[DAE_ALG] = alg + (i + 1) * nz1_;  // fwd:alg
    m->arg[DYN_NUM_IN + DAE_ODE] = ode;  // out:ode
    m->arg[DYN_NUM_IN + DAE_ALG] = alg;  // out:alg
    m->arg[DYN_NUM_IN + DAE_NUM_OUT + DYN_T] = nullptr;  // fwd:t
    m->arg[DYN_NUM_IN + DAE_NUM_OUT + DYN_X] = x + (i + 1) * nx1_;  // fwd:x
    m->arg[DYN_NUM_IN + DAE_NUM_OUT + DYN_Z] = nullptr;  // fwd:z
    m->arg[DYN_NUM_IN + DAE_NUM_OUT + DYN_P] = p + (i + 1) * np1_;  // fwd:p
    m->arg[DYN_NUM_IN + DAE_NUM_OUT + DYN_U] = u + (i + 1) * nu1_;  // fwd:u
    if (calc_sp_reverse(forward_name("daeF", 1), m->arg, m->res, m->iw, m->w)) return 1;
  }
  // Propagate through nondifferentiated
  m->res[DAE_ODE] = ode;  // ode
  m->res[DAE_ALG] = alg;  // alg
  if (calc_sp_reverse("daeF", m->arg, m->res, m->iw, m->w)) return 1;
  return 0;
}

int Integrator::fquad_sp_reverse(SpReverseMem* m, bvec_t* x, bvec_t* z,
    bvec_t* p, bvec_t* u, bvec_t* quad) const {
  // Nondifferentiated inputs
  m->arg[DYN_T] = nullptr;  // t
  m->arg[DYN_X] = x;  // x
  m->arg[DYN_Z] = z;  // z
  m->arg[DYN_P] = p;  // p
  m->arg[DYN_U] = u;  // u
  // Propagate through sensitivities
  for (casadi_int i = 0; i < nfwd_; ++i) {
    m->res[QUAD_QUAD] = quad + (i + 1) * nq1_;  // fwd:quad
    m->arg[DYN_NUM_IN + QUAD_QUAD] = quad;  // out:quad
    m->arg[DYN_NUM_IN + QUAD_NUM_OUT + DYN_T] = nullptr;  // fwd:t
    m->arg[DYN_NUM_IN + QUAD_NUM_OUT + DYN_X] = x + (i + 1) * nx1_;  // fwd:x
    m->arg[DYN_NUM_IN + QUAD_NUM_OUT + DYN_Z] = z + (i + 1) * nz1_;  // fwd:z
    m->arg[DYN_NUM_IN + QUAD_NUM_OUT + DYN_P] = p + (i + 1) * np1_;  // fwd:p
    m->arg[DYN_NUM_IN + QUAD_NUM_OUT + DYN_U] = u + (i + 1) * nu1_;  // fwd:u
    if (calc_sp_reverse(forward_name("quadF", 1), m->arg, m->res, m->iw, m->w)) return 1;
  }
  // Propagate through nondifferentiated
  m->res[QUAD_QUAD] = quad;  // quad
  if (calc_sp_reverse("quadF", m->arg, m->res, m->iw, m->w)) return 1;
  return 0;
}

int Integrator::bdae_sp_reverse(SpReverseMem* m, bvec_t* x, bvec_t* z,
    bvec_t* p, bvec_t* u, bvec_t* adj_ode, bvec_t* rp,
    bvec_t* adj_x, bvec_t* adj_z) const {
  // Nondifferentiated inputs
  m->arg[BDYN_T] = nullptr;  // t
  m->arg[BDYN_X] = x;  // x
  m->arg[BDYN_Z] = z;  // z
  m->arg[BDYN_P] = p;  // p
  m->arg[BDYN_U] = u;  // u
  m->arg[BDYN_OUT_ODE] = nullptr;  // out_ode
  m->arg[BDYN_OUT_ALG] = nullptr;  // out_alg
  m->arg[BDYN_OUT_QUAD] = nullptr;  // out_quad
  m->arg[BDYN_OUT_ZERO] = nullptr;  // out_zero
  m->arg[BDYN_ADJ_ODE] = adj_ode;  // adj_ode
  m->arg[BDYN_ADJ_ALG] = nullptr;  // adj_alg
  m->arg[BDYN_ADJ_QUAD] = rp;  // adj_quad
  m->arg[BDYN_ADJ_ZERO] = nullptr;  // adj_zero
  // Propagate through sensitivities
  for (casadi_int i = 0; i < nfwd_; ++i) {
    m->res[BDAE_ADJ_X] = adj_x + (i + 1) * nrx1_ * nadj_;  // fwd:adj_x
    m->res[BDAE_ADJ_Z] = adj_z + (i + 1) * nrz1_ * nadj_;  // fwd:adj_z
    m->arg[BDYN_NUM_IN + BDAE_ADJ_X] = adj_x;  // out:adj_x
    m->arg[BDYN_NUM_IN + BDAE_ADJ_Z] = adj_z;  // out:adj_z
    m->arg[BDYN_NUM_IN + BDAE_NUM_OUT + BDYN_T] = nullptr;  // fwd:t
    m->arg[BDYN_NUM_IN + BDAE_NUM_OUT + BDYN_X] = x + (i + 1) * nx1_;  // fwd:x
    m->arg[BDYN_NUM_IN + BDAE_NUM_OUT + BDYN_Z] = z + (i + 1) * nz1_;  // fwd:z
    m->arg[BDYN_NUM_IN + BDAE_NUM_OUT + BDYN_P] = p + (i + 1) * np1_;  // fwd:p
    m->arg[BDYN_NUM_IN + BDAE_NUM_OUT + BDYN_U] = u + (i + 1) * nu1_;  // fwd:u
    m->arg[BDYN_NUM_IN + BDAE_NUM_OUT + BDYN_OUT_ODE] = nullptr;  // fwd:out_ode
    m->arg[BDYN_NUM_IN + BDAE_NUM_OUT + BDYN_OUT_ALG] = nullptr;  // fwd:out_alg
    m->arg[BDYN_NUM_IN + BDAE_NUM_OUT + BDYN_OUT_QUAD] = nullptr;  // fwd:out_quad
    m->arg[BDYN_NUM_IN + BDAE_NUM_OUT + BDYN_OUT_ZERO] = nullptr;  // fwd:out_zero
    m->arg[BDYN_NUM_IN + BDAE_NUM_OUT + BDYN_ADJ_ODE] =
      adj_ode + (i + 1) * nrx1_ * nadj_;  // fwd:adj_ode
    m->arg[BDYN_NUM_IN + BDAE_NUM_OUT + BDYN_ADJ_ALG] = nullptr;  // fwd:adj_alg
    m->arg[BDYN_NUM_IN + BDAE_NUM_OUT + BDYN_ADJ_QUAD] =
      rp + (i + 1) * nrz1_ * nadj_;  // fwd:adj_quad
    m->arg[BDYN_NUM_IN + BDAE_NUM_OUT + BDYN_ADJ_ZERO] = nullptr;  // fwd:adj_zero
    if (calc_sp_reverse(forward_name("daeB", 1), m->arg, m->res, m->iw, m->w)) return 1;
  }
  // Propagate through nondifferentiated
  m->res[BDAE_ADJ_X] = adj_x;  // adj_x
  m->res[BDAE_ADJ_Z] = adj_z;  // adj_z
  if (calc_sp_reverse("daeB", m->arg, m->res, m->iw, m->w)) return 1;
  return 0;
}

int Integrator::bquad_sp_reverse(SpReverseMem* m, bvec_t* x, bvec_t* z,
    bvec_t* p, bvec_t* u, bvec_t* adj_ode, bvec_t* rz, bvec_t* adj_quad,
    bvec_t* adj_p, bvec_t* adj_u) const {
  // Nondifferentiated inputs
  m->arg[BDYN_T] = nullptr;  // t
  m->arg[BDYN_X] = x;  // x
  m->arg[BDYN_Z] = z;  // z
  m->arg[BDYN_P] = p;  // p
  m->arg[BDYN_U] = u;  // u
  m->arg[BDYN_OUT_ODE] = adj_ode;  // out_ode
  m->arg[BDYN_OUT_ALG] = rz;  // out_alg
  m->arg[BDYN_OUT_QUAD] = adj_quad;  // out_quad
  m->arg[BDYN_OUT_ZERO] = nullptr;  // out_zero
  m->arg[BDYN_ADJ_ODE] = adj_ode;  // adj_ode
  m->arg[BDYN_ADJ_ALG] = rz;  // adj_alg
  m->arg[BDYN_ADJ_QUAD] = adj_quad;  // adj_quad
  m->arg[BDYN_ADJ_ZERO] = nullptr;  // adj_zero
  // Propagate through sensitivities
  for (casadi_int i = 0; i < nfwd_; ++i) {
    m->res[BQUAD_ADJ_P] = adj_p ? adj_p + (i + 1) * nrq1_ * nadj_ : 0;  // fwd:adj_p
    m->res[BQUAD_ADJ_U] = adj_u ? adj_u + (i + 1) * nuq1_ * nadj_ : 0;  // fwd:adj_u
    m->arg[BDYN_NUM_IN + BQUAD_ADJ_P] = adj_p;  // out:adj_p
    m->arg[BDYN_NUM_IN + BQUAD_ADJ_U] = adj_u;  // out:adj_u
    m->arg[BDYN_NUM_IN + BQUAD_NUM_OUT + BDYN_T] = nullptr;  // fwd:t
    m->arg[BDYN_NUM_IN + BQUAD_NUM_OUT + BDYN_X] = x + (i + 1) * nx1_;  // fwd:x
    m->arg[BDYN_NUM_IN + BQUAD_NUM_OUT + BDYN_Z] = z + (i + 1) * nz1_;  // fwd:z
    m->arg[BDYN_NUM_IN + BQUAD_NUM_OUT + BDYN_P] = p + (i + 1) * np1_;  // fwd:p
    m->arg[BDYN_NUM_IN + BQUAD_NUM_OUT + BDYN_U] = u + (i + 1) * nu1_;  // fwd:u
    m->arg[BDYN_NUM_IN + BQUAD_NUM_OUT + BDYN_OUT_ODE] = nullptr;  // fwd:out_ode
    m->arg[BDYN_NUM_IN + BQUAD_NUM_OUT + BDYN_OUT_ALG] = nullptr;  // fwd:out_alg
    m->arg[BDYN_NUM_IN + BQUAD_NUM_OUT + BDYN_OUT_QUAD] = nullptr;  // fwd:out_quad
    m->arg[BDYN_NUM_IN + BQUAD_NUM_OUT + BDYN_OUT_ZERO] = nullptr;  // fwd:out_zero
    m->arg[BDYN_NUM_IN + BQUAD_NUM_OUT + BDYN_ADJ_ODE] =
      adj_ode + (i + 1) * nrx1_ * nadj_;  // fwd:adj_ode
    m->arg[BDYN_NUM_IN + BQUAD_NUM_OUT + BDYN_ADJ_ALG] =
      rz + (i + 1) * nrz1_ * nadj_;  // fwd:adj_alg
    m->arg[BDYN_NUM_IN + BQUAD_NUM_OUT + BDYN_ADJ_QUAD] =
      adj_quad + (i + 1) * nrp1_ * nadj_;  // fwd:adj_quad
    m->arg[BDYN_NUM_IN + BQUAD_NUM_OUT + BDYN_ADJ_ZERO] = nullptr;  // fwd:adj_zero
    if (calc_sp_reverse(forward_name("quadB", 1), m->arg, m->res, m->iw, m->w)) return 1;
  }
  // Propagate through nondifferentiated
  m->res[BQUAD_ADJ_P] = adj_p;  // adj_p
  m->res[BQUAD_ADJ_U] = adj_u;  // adj_u
  if (calc_sp_reverse("quadB", m->arg, m->res, m->iw, m->w)) return 1;
  return 0;
}

int Integrator::sp_reverse(bvec_t** arg, bvec_t** res,
    casadi_int* iw, bvec_t* w, void* mem) const {
  if (verbose_) casadi_message(name_ + "::sp_reverse");

  // Inputs
  bvec_t* x0 = arg[INTEGRATOR_X0];
  bvec_t* p = arg[INTEGRATOR_P];
  bvec_t* u = arg[INTEGRATOR_U];
  bvec_t* adj_xf = arg[INTEGRATOR_ADJ_XF];
  bvec_t* rp = arg[INTEGRATOR_ADJ_QF];
  arg += n_in_;

  // Outputs
  bvec_t* xf = res[INTEGRATOR_XF];
  bvec_t* zf = res[INTEGRATOR_ZF];
  bvec_t* qf = res[INTEGRATOR_QF];
  bvec_t* adj_x0 = res[INTEGRATOR_ADJ_X0];
  bvec_t* adj_p0 = res[INTEGRATOR_ADJ_P];
  bvec_t* adj_u = res[INTEGRATOR_ADJ_U];
  res += n_out_;

  // Work vectors
  bvec_t *ode = w; w += nx_;
  bvec_t *alg = w; w += nz_;
  bvec_t *x_prev = w; w += nx_;
  bvec_t *rx = w; w += nrx_;
  bvec_t *rz = w; w += nrz_;
  bvec_t *adj_ode = w; w += nrx_;
  bvec_t *rq = w; w += nrq_;

  // Memory struct for function calls below
  SpReverseMem m = {arg, res, iw, w};

  // Clear state vector
  std::fill_n(ode, nx_, 0);
  std::fill_n(alg, nz_, 0);

  if (nrx_ > 0) {
    // Propagate from adj_x0 initial time
    if (adj_x0) {
      std::copy_n(adj_x0, nrx_, rx);
      std::fill_n(adj_x0, nrx_, 0);
    } else {
      std::fill_n(rx, nrx_, 0);
    }
    // Reset rz
    std::fill_n(rz, nrz_, 0);

    // Save adj_p0: See note below
    if (adj_p0) std::copy_n(adj_p0, nrq_, rq);

    // Step backwards through backward problem
    for (casadi_int k = 0; k < nt(); ++k) {
      // Restore adj_p0: See note below
      if (adj_p0) std::copy_n(rq, nrq_, adj_p0);

      // Add impulse from adj_xf
      if (adj_xf) {
        for (casadi_int i = 0; i < nrx_; ++i) rx[i] |= adj_xf[i];
        std::fill_n(adj_xf, nrx_, 0);
      }

      // Get dependencies from backward quadratures
      if ((nrq_ > 0 && adj_p0) || (nuq_ > 0 && adj_u)) {
        if (bquad_sp_reverse(&m, ode, alg, p, u, rx, rz, rp, adj_p0, adj_u)) return 1;
      }

      // Propagate interdependencies
      std::fill_n(w, nrx_+nrz_, 0);
      sp_jac_rdae_.spsolve(w, rx, true);
      std::copy_n(w, nrx_+nrz_, rx);

      // Direct dependency adj_ode -> rx
      std::copy_n(rx, nrx_, adj_ode);

      // Indirect dependency via g
      if (bdae_sp_reverse(&m, ode, alg, p, u, adj_ode, rp, rx, rz)) return 1;

      // Update rx, rz
      std::copy_n(adj_ode, nrx_, rx);
      std::fill_n(rz, nrz_, 0);

      // Shift time
      if (adj_xf) adj_xf += nrx_;
      if (rp) rp += nrp_;
      if (adj_u) adj_u += nuq_;
      if (u) u += nu_;
    }
  } else {
    // Take u past the last grid point
    if (u) u += nu_ * nt();
  }

  // Take xf, zf, qf past the last grid point
  if (xf) xf += nx_ * nt();
  if (zf) zf += nz_ * nt();
  if (qf) qf += nq_ * nt();

  // Step backwards through forward problem
  for (casadi_int k = nt(); k-- > 0; ) {
    // Shift time
    if (xf) xf -= nx_;
    if (zf) zf -= nz_;
    if (qf) qf -= nq_;
    if (u) u -= nu_;

    // Add impulse from outputs
    if (xf) {
      for (casadi_int i = 0; i < nx_; ++i) ode[i] |= xf[i];
      std::fill_n(xf, nx_, 0);
    }
    if (zf) {
      for (casadi_int i = 0; i < nz_; ++i) alg[i] |= zf[i];
      std::fill_n(zf, nz_, 0);
    }

    // Get dependencies from forward quadratures, if any
    if (nq_ > 0 && qf) {
      if (fquad_sp_reverse(&m, ode, alg, p, u, qf)) return 1;
    }

    // Propagate interdependencies
    std::fill_n(w, nx_ + nz_, 0);
    sp_jac_dae_.spsolve(w, ode, true);
    std::copy_n(w, nx_ + nz_, ode);

    // Direct dependency x_prev -> x
    std::copy_n(ode, nx_, x_prev);

    // Indirect dependency through f
    if (fdae_sp_reverse(&m, x_prev, p, u, ode, alg)) return 1;

    // Update x, z
    std::copy_n(x_prev, nx_, ode);
    std::fill_n(alg, nz_, 0);
  }

  // Direct dependency x0 -> x
  if (x0) {
    for (casadi_int i = 0; i < nx_; ++i) x0[i] |= x_prev[i];
  }

  return 0;
}

Function Integrator::get_forward(casadi_int nfwd, const std::string& name,
    const std::vector<std::string>& inames,
    const std::vector<std::string>& onames,
    const Dict& opts) const {
  if (verbose_) casadi_message(name_ + "::get_forward");

  // Integrator options
  Dict aug_opts = getDerivativeOptions(true);
  for (auto&& i : augmented_options_) {
    aug_opts[i.first] = i.second;
  }

  // Get current DAE, with any existing sensitivity equations augmented
  Function this_dae = augmented_dae();

  // Create integrator for augmented DAE
  std::string aug_prefix = "fsens" + str(nfwd) + "_";
  aug_opts["derivative_of"] = self();
  aug_opts["nfwd"] = nfwd;
  aug_opts["nadj"] = nadj_;
  Function aug_int = integrator(aug_prefix + name_, plugin_name(),
    this_dae, t0_, tout_, aug_opts);

  // All inputs of the return function
  std::vector<MX> ret_in;
  ret_in.reserve(INTEGRATOR_NUM_IN + INTEGRATOR_NUM_OUT + INTEGRATOR_NUM_IN);

  // Add nondifferentiated inputs to ret_in
  for (casadi_int i = 0; i < INTEGRATOR_NUM_IN; ++i) {
    ret_in.push_back(MX::sym(integrator_in(i), sparsity_in(i)));
  }

  // Add nondifferentiated outputs (unused) to ret_in
  for (casadi_int i = 0; i < INTEGRATOR_NUM_OUT; ++i) {
    ret_in.push_back(MX::sym("out_" + integrator_out(i), Sparsity(size_out(i))));
  }

  // Create symbolic expressions for augmented problem, add forward seeds to ret_in
  std::vector<std::vector<MX>> aug_in(INTEGRATOR_NUM_IN);
  std::vector<MX> v(nfwd);
  for (casadi_int i = 0; i < INTEGRATOR_NUM_IN; ++i) {
    for (casadi_int d = 0; d < nfwd; ++d) {
      v[d] = MX::sym("fwd" + str(d) + "_" + integrator_in(i), sparsity_in(i));
      aug_in[i].push_back(v[d]);
    }
    ret_in.push_back(horzcat(v));
  }

  // Call the augmented integrator
  std::vector<MX> integrator_in(INTEGRATOR_NUM_IN);
  for (casadi_int i = 0; i < INTEGRATOR_NUM_IN; ++i) {
    if (size1_in(i) > 0 && grid_in(i) && nt() > 1) {
      // Split nondifferentiated input by grid point
      std::vector<MX> ret_in_split = horzsplit_n(ret_in[i], nt());
      // Split augmented input by grid point
      std::vector<std::vector<MX>> aug_in_split(nfwd);
      for (casadi_int d = 0; d < nfwd; ++d) {
        aug_in_split[d] = horzsplit_n(aug_in[i][d], nt());
      }
      // Reorder columns
      v.clear();
      for (casadi_int k = 0; k < nt(); ++k) {
        v.push_back(ret_in_split.at(k));
        for (casadi_int d = 0; d < nfwd; ++d) {
          v.push_back(aug_in_split[d].at(k));
        }
      }
    } else {
      // No reordering necessary
      v = aug_in[i];
      v.insert(v.begin(), ret_in[i]);
    }
    // Flatten all elements
    for (MX& e : v) e = vec(e);
    integrator_in[i] = horzcat(v);
  }
  std::vector<MX> integrator_out = aug_int(integrator_in);

  // Collect forward sensitivites
  std::vector<MX> ret_out;
  ret_out.reserve(INTEGRATOR_NUM_OUT);
  for (casadi_int i = 0; i < INTEGRATOR_NUM_OUT; ++i) {
    // Split return by grid points and sensitivities
    casadi_int n_grid = grid_out(i) ? nt() : 1;
    std::vector<casadi_int> offset = {0};
    for (casadi_int k = 0; k < n_grid; ++k) {
      for (casadi_int d = 0; d <= nfwd; ++d) {
        offset.push_back(offset.back() + size2_out(i) / n_grid);
      }
    }
    std::vector<MX> integrator_out_split = horzsplit(
      reshape(integrator_out[i], size1_out(i), offset.back()), offset);
    // Collect sensitivity blocks in the right order
    std::vector<MX> ret_out_split;
    ret_out_split.reserve(n_grid * nfwd);
    for (casadi_int d = 0; d < nfwd; ++d) {
      for (casadi_int k = 0; k < n_grid; ++k) {
        ret_out_split.push_back(integrator_out_split.at((nfwd + 1) * k + d + 1));
      }
    }
    ret_out.push_back(horzcat(ret_out_split));
  }

  Dict options = opts;
  options["allow_duplicate_io_names"] = true;

  // Create derivative function and return
  return Function(name, ret_in, ret_out, inames, onames, options);
}

Function Integrator::get_reverse(casadi_int nadj, const std::string& name,
    const std::vector<std::string>& inames,
    const std::vector<std::string>& onames,
    const Dict& opts) const {
  if (verbose_) casadi_message(name_ + "::get_reverse");

  // Integrator options
  Dict aug_opts = getDerivativeOptions(false);
  for (auto&& i : augmented_options_) {
    aug_opts[i.first] = i.second;
  }

  // Get the current oracle, augmented with any existing forward sensitivity equations
  Function this_dae = augmented_dae();

  // Create integrator for augmented DAE
  std::string aug_prefix = "asens" + str(nadj) + "_";
  aug_opts["derivative_of"] = self();
  if (nrx_ == 0) {
    // Add backward problem
    aug_opts["nadj"] = nadj;
    aug_opts["nfwd"] = 0;
  } else {
    // Reformulate as forward-over-reverse
    aug_opts["nfwd"] = nadj;
    aug_opts["nadj"] = nadj_;
  }
  Function aug_int = integrator(aug_prefix + name_, plugin_name(),
    this_dae, t0_, tout_, aug_opts);

  // All inputs of the return function
  std::vector<MX> ret_in;
  ret_in.reserve(INTEGRATOR_NUM_IN + INTEGRATOR_NUM_OUT + INTEGRATOR_NUM_IN);

  // Add nondifferentiated inputs to ret_in
  for (casadi_int i = 0; i < INTEGRATOR_NUM_IN; ++i) {
    ret_in.push_back(MX::sym(integrator_in(i), sparsity_in(i)));
  }

  // Add nondifferentiated outputs (unused) to ret_in
  for (casadi_int i = 0; i < INTEGRATOR_NUM_OUT; ++i) {
    ret_in.push_back(MX::sym("out_" + integrator_out(i), Sparsity(size_out(i))));
  }

  // Create symbolic expressions for augmented problem, add adjoint seeds to ret_in
  std::vector<std::vector<MX>> aug_in(INTEGRATOR_NUM_OUT);
  std::vector<MX> v(nadj);
  for (casadi_int i = 0; i < INTEGRATOR_NUM_OUT; ++i) {
    for (casadi_int d=0; d<nadj; ++d) {
      v[d] = MX::sym("adj" + str(d) + "_" + integrator_out(i), sparsity_out(i));
      aug_in[i].push_back(v[d]);
    }
    ret_in.push_back(horzcat(v));
  }

  // Call the augmented integrator
  std::vector<MX> integrator_in(INTEGRATOR_NUM_IN);
  for (casadi_int i = 0; i < INTEGRATOR_NUM_IN; ++i) {
    // Output index contributing to adjoint seeds
    casadi_int j = adjmap_out(i);
    // Number of grid points for this integrator input
    casadi_int n_grid = grid_in(i) ? nt() : 1;
    // Split input and seeds by grid points, if necessary
    std::vector<MX> ret_in_split;
    std::vector<std::vector<MX>> aug_in_split(nadj);
    if (size1_in(i) > 0 && grid_in(i) && n_grid > 1) {
      // Split nondifferentiated input by grid point
      ret_in_split = horzsplit_n(ret_in[i], nt());
      // Split augmented input by grid point
      for (casadi_int d = 0; d < nadj; ++d) {
        aug_in_split[d] = horzsplit_n(aug_in[j][d], nt());
      }
    } else {
      // No reordering necessary
      ret_in_split = {ret_in[i]};
      for (casadi_int d = 0; d < nadj; ++d) aug_in_split[d] = {aug_in[j][d]};
    }
    // Vectorize all inputs to allow concatenation (unlike forward sensitivities,
    // number of rows for sensitivities may be different from original inputs)
    for (auto&& e : ret_in_split) e = vec(e);
    for (auto&& e1 : aug_in_split) {
      for (auto&& e2 : e1) e2 = vec(e2);
    }
    // Assemble input argument
    v.clear();
    for (casadi_int k = 0; k < ret_in_split.size(); ++k) {
      v.push_back(ret_in_split.at(k));
      for (casadi_int d = 0; d < nadj; ++d) {
        v.push_back(aug_in_split[d].at(k));
      }
    }
    integrator_in[i] = reshape(vertcat(v), aug_int.size_in(i));
  }
  std::vector<MX> integrator_out = aug_int(integrator_in);

  // Collect adjoint sensitivites
  std::vector<MX> ret_out;
  ret_out.reserve(INTEGRATOR_NUM_IN);
  for (casadi_int i = 0; i < INTEGRATOR_NUM_IN; ++i) {
    casadi_int j = adjmap_out(i);
    // Split return by grid points and sensitivities
    casadi_int n_grid = grid_out(j) ? nt() : 1;
    std::vector<casadi_int> offset = {0};
    for (casadi_int k = 0; k < n_grid; ++k) {
      offset.push_back(offset.back() + numel_out(j) / n_grid);
      for (casadi_int d = 0; d < nadj; ++d) {
        offset.push_back(offset.back() + numel_in(i) / n_grid);
      }
    }
    std::vector<MX> integrator_out_split = vertsplit(vec(integrator_out[j]), offset);
    // Collect sensitivity blocks in the right order
    std::vector<MX> ret_out_split;
    ret_out_split.reserve(n_grid * nadj);
    for (casadi_int d = 0; d < nadj; ++d) {
      for (casadi_int k = 0; k < n_grid; ++k) {
        ret_out_split.push_back(reshape(integrator_out_split.at((nadj + 1) * k + d + 1),
          size1_in(i), size2_in(i) / n_grid));
      }
    }
    ret_out.push_back(horzcat(ret_out_split));
  }

  Dict options = opts;
  options["allow_duplicate_io_names"] = true;

  // Create derivative function and return
  return Function(name, ret_in, ret_out, inames, onames, options);
}

Dict Integrator::getDerivativeOptions(bool fwd) const {
  // Copy all options
  return opts_;
}

Sparsity Integrator::sp_jac_aug(const Sparsity& J, const Sparsity& J1) const {
  // Row 1, column 2 in the augmented Jacobian
  Sparsity J12(J.size1(), nfwd_ * J.size2());
  // Row 2, column 1 in the augmented Jacobian
  Sparsity J21 = vertcat(std::vector<Sparsity>(nfwd_, J1));
  // Row 2, column 2 in the augmented Jacobian
  Sparsity J22 = diagcat(std::vector<Sparsity>(nfwd_, J));
  // Form block matrix
  return blockcat(J, J12, J21, J22);
}


Sparsity Integrator::sp_jac_dae() {
  // Get the functions
  const Function& F = get_function("daeF");
  // Sparsity pattern for nonaugmented system
  Sparsity J_xx = F.jac_sparsity(DAE_ODE, DYN_X) + Sparsity::diag(nx1_);
  Sparsity J_xz = F.jac_sparsity(DAE_ODE, DYN_Z);
  Sparsity J_zx = F.jac_sparsity(DAE_ALG, DYN_X);
  Sparsity J_zz = F.jac_sparsity(DAE_ALG, DYN_Z);
  // Augment with sensitivity equations
  if (nfwd_ > 0) {
    const Function& fwd_F = get_function(forward_name("daeF", 1));
    J_xx = sp_jac_aug(J_xx, fwd_F.jac_sparsity(DAE_ODE, DYN_X));
    J_xz = sp_jac_aug(J_xz, fwd_F.jac_sparsity(DAE_ODE, DYN_Z));
    J_zx = sp_jac_aug(J_zx, fwd_F.jac_sparsity(DAE_ALG, DYN_X));
    J_zz = sp_jac_aug(J_zz, fwd_F.jac_sparsity(DAE_ALG, DYN_Z));
  }
  // Assemble the block matrix
  return blockcat(J_xx, J_xz, J_zx, J_zz);
}

Sparsity Integrator::sp_jac_rdae() {
  // Get the functions
  const Function& G = get_function("daeB");
  // Sparsity pattern for nonaugmented system
  Sparsity J_xx = G.jac_sparsity(BDAE_ADJ_X, BDYN_ADJ_ODE) + Sparsity::diag(nrx1_ * nadj_);
  Sparsity J_xz = G.jac_sparsity(BDAE_ADJ_X, BDYN_ADJ_ALG);
  Sparsity J_zx = G.jac_sparsity(BDAE_ADJ_Z, BDYN_ADJ_ODE);
  Sparsity J_zz = G.jac_sparsity(BDAE_ADJ_Z, BDYN_ADJ_ALG);
  // Augment with sensitivity equations
  if (nfwd_ > 0) {
    const Function& fwd_G = get_function(forward_name("daeB", 1));
    J_xx = sp_jac_aug(J_xx, fwd_G.jac_sparsity(BDAE_ADJ_X, BDYN_ADJ_ODE));
    J_xz = sp_jac_aug(J_xz, fwd_G.jac_sparsity(BDAE_ADJ_X, BDYN_ADJ_ALG));
    J_zx = sp_jac_aug(J_zx, fwd_G.jac_sparsity(BDAE_ADJ_Z, BDYN_ADJ_ODE));
    J_zz = sp_jac_aug(J_zz, fwd_G.jac_sparsity(BDAE_ADJ_Z, BDYN_ADJ_ALG));
  }
  // Assemble the block matrix
  return blockcat(J_xx, J_xz, J_zx, J_zz);
}

std::map<std::string, Integrator::Plugin> Integrator::solvers_;

#ifdef CASADI_WITH_THREADSAFE_SYMBOLICS
  std::mutex Integrator::mutex_solvers_;
#endif // CASADI_WITH_THREADSAFE_SYMBOLICS

const std::string Integrator::infix_ = "integrator";

FixedStepIntegrator::FixedStepIntegrator(const std::string& name, const Function& dae,
    double t0, const std::vector<double>& tout) : Integrator(name, dae, t0, tout) {

  // Default options
  nk_target_ = 20;
}

FixedStepIntegrator::~FixedStepIntegrator() {
  clear_mem();
}

const Options FixedStepIntegrator::options_
= {{&Integrator::options_},
    {{"number_of_finite_elements",
      {OT_INT,
      "Target number of finite elements. "
      "The actual number may be higher to accommodate all output times"}},
    {"simplify",
      {OT_BOOL,
      "Implement as MX Function (codegeneratable/serializable) default: false"}},
    {"simplify_options",
      {OT_DICT,
      "Any options to pass to simplified form Function constructor"}}
    }
};


Function FixedStepIntegrator::create_advanced(const Dict& opts) {
  Function temp = Function::create(this, opts);

  // Check if we need to simplify
  bool simplify = false;
  auto it = opts.find("simplify");
  if (it != opts.end()) simplify = it->second;

  if (simplify && nrx_==0 && nt()==1) {
    // Retrieve explicit simulation step (one finite element)
    Function F = get_function("step");

    MX z0 = MX::sym("z0", sparsity_in(INTEGRATOR_Z0));

    // Create symbols
    std::vector<MX> F_in = F.mx_in();

    // Prepare return Function inputs
    std::vector<MX> intg_in(INTEGRATOR_NUM_IN);
    intg_in[INTEGRATOR_X0] = F_in[STEP_X0];
    intg_in[INTEGRATOR_Z0] = z0;
    intg_in[INTEGRATOR_P] = F_in[STEP_P];
    intg_in[INTEGRATOR_U] = F_in[STEP_U];
    F_in[STEP_V0] = algebraic_state_init(intg_in[INTEGRATOR_X0], z0);

    // Number of finite elements and time steps
    double h = (tout_.back() - t0_)/static_cast<double>(disc_.back());

    // Prepare return Function outputs
    std::vector<MX> intg_out(INTEGRATOR_NUM_OUT);
    F_in[STEP_T] = t0_;
    F_in[STEP_H] = h;

    std::vector<MX> F_out;
    // Loop over finite elements
    for (casadi_int k=0; k<disc_.back(); ++k) {
      F_out = F(F_in);

      F_in[STEP_X0] = F_out[STEP_XF];
      F_in[STEP_V0] = F_out[STEP_VF];
      intg_out[INTEGRATOR_QF] = k==0? F_out[STEP_QF] : intg_out[INTEGRATOR_QF]+F_out[STEP_QF];
      F_in[STEP_T] += h;
    }

    intg_out[INTEGRATOR_XF] = F_out[STEP_XF];

    // If-clause needed because rk abuses STEP_VF output for intermediate state output
    if (nz_) {
      intg_out[INTEGRATOR_ZF] = algebraic_state_output(F_out[STEP_VF]);
    }

    // Extract options for Function constructor
    Dict sopts;
    sopts["print_time"] = print_time_;
    auto it = opts.find("simplify_options");
    if (it!=opts.end()) update_dict(sopts, it->second);

    return Function(temp.name(), intg_in, intg_out, integrator_in(), integrator_out(), sopts);
  } else {
    return temp;
  }
}

void FixedStepIntegrator::init(const Dict& opts) {
  // Call the base class init
  Integrator::init(opts);

  // Instantiate functions, forward and backward problem
  set_function(oracle_, "dae");
  if (nadj_ > 0) set_function(rdae_, "rdae");

  // Read options
  for (auto&& op : opts) {
    if (op.first=="number_of_finite_elements") {
      nk_target_ = op.second;
    }
  }

  // Consistency check
  casadi_assert(nk_target_ > 0, "Number of finite elements must be strictly positive");

  // Target interval length
  double h_target = (tout_.back() - t0_) / nk_target_;

  // Number of finite elements for each control interval and in total
  disc_.reserve(1 + nt());
  disc_.push_back(0);
  double t_cur = t0_;
  for (double t_next : tout_) {
    disc_.push_back(disc_.back() + std::ceil((t_next - t_cur) / h_target));
    t_cur = t_next;
  }

  // Setup discrete time dynamics
  setup_step();

  // Get discrete time dimensions
  const Function& F = get_function(has_function("step") ? "step" : "implicit_step");
  nv1_ = F.nnz_out(STEP_VF);
  nrv1_ = nv1_ * nadj_;
  nv_ = nv1_ * (1 + nfwd_);
  nrv_ = nrv1_ * (1 + nfwd_);

  // Work vectors, forward problem
  alloc_w(nv_, true); // v
  alloc_w(nq_, true); // q
  alloc_w(nv_, true); // v_prev
  alloc_w(nq_, true); // q_prev

  // Work vectors, backward problem
  alloc_w(nrv_, true); // rv
  alloc_w(nuq_, true); // adj_u
  alloc_w(nrq_, true); // rq_prev
  alloc_w(nuq_, true); // adj_u_prev

  // Allocate tape if backward states are present
  if (nrx_ > 0) {
    alloc_w((disc_.back() + 1) * nx_, true); // x_tape
    alloc_w(disc_.back() * nv_, true); // v_tape
  }
}

void FixedStepIntegrator::set_work(void* mem, const double**& arg, double**& res,
    casadi_int*& iw, double*& w) const {
  auto m = static_cast<FixedStepMemory*>(mem);

  // Set work in base classes
  Integrator::set_work(mem, arg, res, iw, w);

  // Work vectors, forward problem
  m->v = w; w += nv_;
  m->q = w; w += nq_;
  m->v_prev = w; w += nv_;
  m->q_prev = w; w += nq_;

  // Work vectors, backward problem
  m->rv = w; w += nrv_;
  m->adj_u = w; w += nuq_;
  m->rq_prev = w; w += nrq_;
  m->adj_u_prev = w; w += nuq_;

  // Allocate tape if backward states are present
  if (nrx_ > 0) {
    m->x_tape = w; w += (disc_.back() + 1) * nx_;
    m->v_tape = w; w += disc_.back() * nv_;
  }
}

int FixedStepIntegrator::init_mem(void* mem) const {
  if (Integrator::init_mem(mem)) return 1;
  // auto m = static_cast<FixedStepMemory*>(mem);

  return 0;
}

void FixedStepIntegrator::advance(IntegratorMemory* mem,
    const double* u, double* x, double* z, double* q) const {
  auto m = static_cast<FixedStepMemory*>(mem);

  // Set controls
  casadi_copy(u, nu_, m->u);

  // Number of finite elements and time steps
  casadi_int nj = disc_[m->k + 1] - disc_[m->k];
  double h = (m->t_next - m->t) / nj;

  // Take steps
  for (casadi_int j = 0; j < nj; ++j) {
    // Current time
    double t = m->t + j * h;

    // Update the previous step
    casadi_copy(m->x, nx_, m->x_prev);
    casadi_copy(m->v, nv_, m->v_prev);
    casadi_copy(m->q, nq_, m->q_prev);

    // Take step
    stepF(m, t, h, m->x_prev, m->v_prev, m->x, m->v, m->q);
    casadi_axpy(nq_, 1., m->q_prev, m->q);

    // Save state, if needed
    if (nrx_ > 0) {
      casadi_int tapeind = disc_[m->k] + j;
      casadi_copy(m->x, nx_, m->x_tape + nx_ * (tapeind + 1));
      casadi_copy(m->v, nv_, m->v_tape + nv_ * tapeind);
    }
  }

  // Return to user
  casadi_copy(m->x, nx_, x);
  casadi_copy(m->v + nv_ - nz_, nz_, z);
  casadi_copy(m->q, nq_, q);
}

void FixedStepIntegrator::retreat(IntegratorMemory* mem, const double* u,
    double* rx, double* rq, double* adj_u) const {
  auto m = static_cast<FixedStepMemory*>(mem);

  // Set controls
  casadi_copy(u, nu_, m->u);

  // Number of finite elements and time steps
  casadi_int nj = disc_[m->k + 1] - disc_[m->k];
  double h = (m->t - m->t_next) / nj;

  // Take steps
  for (casadi_int j = nj; j-- > 0; ) {
    // Current time
    double t = m->t_next + j * h;

    // Update the previous step
    casadi_copy(m->rx, nrx_, m->adj_ode);
    casadi_copy(m->rq, nrq_, m->rq_prev);
    casadi_copy(m->adj_u, nuq_, m->adj_u_prev);

    // Take step
    casadi_int tapeind = disc_[m->k] + j;
    stepB(m, t, h,
      m->x_tape + nx_ * tapeind, m->x_tape + nx_ * (tapeind + 1),
      m->v_tape + nv_ * tapeind,
      m->adj_ode, m->rv, m->rx, m->rq, m->adj_u);
    casadi_clear(m->rv, nrv_);
    casadi_axpy(nrq_, 1., m->rq_prev, m->rq);
    casadi_axpy(nuq_, 1., m->adj_u_prev, m->adj_u);
  }

  // Return to user
  casadi_copy(m->rx, nrx_, rx);
  casadi_copy(m->rq, nrq_, rq);
  casadi_copy(m->adj_u, nuq_, adj_u);
}

void FixedStepIntegrator::stepF(FixedStepMemory* m, double t, double h,
    const double* x0, const double* v0, double* xf, double* vf, double* qf) const {
  // Evaluate nondifferentiated
  std::fill(m->arg, m->arg + STEP_NUM_IN, nullptr);
  m->arg[STEP_T] = &t;  // t
  m->arg[STEP_H] = &h;  // h
  m->arg[STEP_X0] = x0;  // x0
  m->arg[STEP_V0] = v0;  // v0
  m->arg[STEP_P] = m->p;  // p
  m->arg[STEP_U] = m->u;  // u
  std::fill(m->res, m->res + STEP_NUM_OUT, nullptr);
  m->res[STEP_XF] = xf;  // xf
  m->res[STEP_VF] = vf;  // vf
  m->res[STEP_QF] = qf;  // qf
  calc_function(m, "step");
  // Evaluate sensitivities
  if (nfwd_ > 0) {
    m->arg[STEP_NUM_IN + STEP_XF] = xf;  // out:xf
    m->arg[STEP_NUM_IN + STEP_VF] = vf;  // out:vf
    m->arg[STEP_NUM_IN + STEP_QF] = qf;  // out:qf
    m->arg[STEP_NUM_IN + STEP_NUM_OUT + STEP_T] = nullptr;  // fwd:t
    m->arg[STEP_NUM_IN + STEP_NUM_OUT + STEP_H] = nullptr;  // fwd:h
    m->arg[STEP_NUM_IN + STEP_NUM_OUT + STEP_X0] = x0 + nx1_;  // fwd:x0
    m->arg[STEP_NUM_IN + STEP_NUM_OUT + STEP_V0] = v0 + nv1_;  // fwd:v0
    m->arg[STEP_NUM_IN + STEP_NUM_OUT + STEP_P] = m->p + np1_;  // fwd:p
    m->arg[STEP_NUM_IN + STEP_NUM_OUT + STEP_U] = m->u + nu1_;  // fwd:u
    m->res[STEP_XF] = xf + nx1_;  // fwd:xf
    m->res[STEP_VF] = vf + nv1_;  // fwd:vf
    m->res[STEP_QF] = qf + nq1_;  // fwd:qf
    calc_function(m, forward_name("step", nfwd_));
  }
}

void FixedStepIntegrator::stepB(FixedStepMemory* m, double t, double h,
    const double* x0, const double* xf, const double* vf,
    const double* adj_xf, const double* rv0,
    double* adj_x0, double* adj_p, double* adj_u) const {
  // Evaluate nondifferentiated
  std::fill(m->arg, m->arg + BSTEP_NUM_IN, nullptr);
  m->arg[BSTEP_T] = &t;  // t
  m->arg[BSTEP_H] = &h;  // h
  m->arg[BSTEP_X0] = x0;  // x0
  m->arg[BSTEP_V0] = nullptr;  // v0
  m->arg[BSTEP_P] = m->p;  // p
  m->arg[BSTEP_U] = m->u;  // u
  m->arg[BSTEP_OUT_XF] = xf;  // out:xf
  m->arg[BSTEP_OUT_VF] = vf;  // out:vf
  m->arg[BSTEP_OUT_QF] = nullptr;  // out:qf
  m->arg[BSTEP_ADJ_XF] = adj_xf;  // adj:xf
  m->arg[BSTEP_ADJ_VF] = rv0;  // adj:vf
  m->arg[BSTEP_ADJ_QF] = m->rp;  // adj:qf
  std::fill(m->res, m->res + BSTEP_NUM_OUT, nullptr);
  m->res[BSTEP_ADJ_T] = nullptr;  // adj:t
  m->res[BSTEP_ADJ_H] = nullptr;  // adj:h
  m->res[BSTEP_ADJ_X0] = adj_x0;  // adj:x0
  m->res[BSTEP_ADJ_V0] = nullptr;  // adj:v0
  m->res[BSTEP_ADJ_P] = adj_p;  // adj:p
  m->res[BSTEP_ADJ_U] = adj_u;  // adj:u
  calc_function(m, reverse_name("step", nadj_));
  // Evaluate sensitivities
  if (nfwd_ > 0) {
    m->arg[BSTEP_NUM_IN + BSTEP_ADJ_T] = nullptr;  // out:adj:t
    m->arg[BSTEP_NUM_IN + BSTEP_ADJ_H] = nullptr;  // out:adj:h
    m->arg[BSTEP_NUM_IN + BSTEP_ADJ_X0] = adj_x0;  // out:adj:x0
    m->arg[BSTEP_NUM_IN + BSTEP_ADJ_V0] = nullptr;  // out:adj:v0
    m->arg[BSTEP_NUM_IN + BSTEP_ADJ_P] = adj_p;  // out:adj:p
    m->arg[BSTEP_NUM_IN + BSTEP_ADJ_U] = adj_u;  // out:adj:u
    m->arg[BSTEP_NUM_IN + BSTEP_NUM_OUT + BSTEP_T] = nullptr;  // fwd:t
    m->arg[BSTEP_NUM_IN + BSTEP_NUM_OUT + BSTEP_H] = nullptr;  // fwd:h
    m->arg[BSTEP_NUM_IN + BSTEP_NUM_OUT + BSTEP_X0] = x0 + nx1_;  // fwd:x0
    m->arg[BSTEP_NUM_IN + BSTEP_NUM_OUT + BSTEP_V0] = nullptr;  // fwd:v0
    m->arg[BSTEP_NUM_IN + BSTEP_NUM_OUT + BSTEP_P] = m->p + np1_;  // fwd:p
    m->arg[BSTEP_NUM_IN + BSTEP_NUM_OUT + BSTEP_U] = m->u + nu1_;  // fwd:u
    m->arg[BSTEP_NUM_IN + BSTEP_NUM_OUT + BSTEP_OUT_XF] = xf + nx1_;  // fwd:out:xf
    m->arg[BSTEP_NUM_IN + BSTEP_NUM_OUT + BSTEP_OUT_VF] = vf + nv1_;  // fwd:out:vf
    m->arg[BSTEP_NUM_IN + BSTEP_NUM_OUT + BSTEP_OUT_QF] = nullptr;  // fwd:out:qf
    m->arg[BSTEP_NUM_IN + BSTEP_NUM_OUT + BSTEP_ADJ_XF] = adj_xf + nrx1_ * nadj_;  // fwd:adj:xf
    m->arg[BSTEP_NUM_IN + BSTEP_NUM_OUT + BSTEP_ADJ_VF] = rv0 + nrv1_;  // fwd:adj:vf
    m->arg[BSTEP_NUM_IN + BSTEP_NUM_OUT + BSTEP_ADJ_QF] = m->rp + nrp1_ * nadj_;  // fwd:adj:qf
    m->res[BSTEP_ADJ_T] = nullptr;  // fwd:adj:t
    m->res[BSTEP_ADJ_H] = nullptr;  // fwd:adj:h
    m->res[BSTEP_ADJ_X0] = adj_x0 + nrx1_ * nadj_;  // fwd:adj_x0
    m->res[BSTEP_ADJ_V0] = nullptr;  // fwd:adj:v0
    m->res[BSTEP_ADJ_P] = adj_p + nrq1_ * nadj_;  // fwd:adj_p
    m->res[BSTEP_ADJ_U] = adj_u + nuq1_ * nadj_;  // fwd:adj_u
    calc_function(m, forward_name(reverse_name("step", nadj_), nfwd_));
  }
}

void FixedStepIntegrator::reset(IntegratorMemory* mem, const double* u, const double* x,
    const double* z, const double* p) const {
  auto m = static_cast<FixedStepMemory*>(mem);

  // Reset the base classes
  Integrator::reset(mem, u, x, z, p);

  // Set parameters
  casadi_copy(p, np_, m->p);

  // Set controls
  casadi_copy(u, nu_, m->u);

  // Update the state
  casadi_copy(x, nx_, m->x);
  casadi_copy(z, nz_, m->z);

  // Reset summation states
  casadi_clear(m->q, nq_);

  // Get consistent initial conditions
  casadi_fill(m->v, nv_, std::numeric_limits<double>::quiet_NaN());

  // Add the first element in the tape
  if (nrx_ > 0) {
    casadi_copy(x, nx_, m->x_tape);
  }
}

void FixedStepIntegrator::resetB(IntegratorMemory* mem) const {
  auto m = static_cast<FixedStepMemory*>(mem);

  // Clear adjoint seeds
  casadi_clear(m->rp, nrp_);
  casadi_clear(m->rx, nrx_);

  // Reset summation states
  casadi_clear(m->rq, nrq_);
  casadi_clear(m->adj_u, nuq_);

  // Update backwards dependent variables
  casadi_clear(m->rv, nrv_);
}

void FixedStepIntegrator::impulseB(IntegratorMemory* mem,
    const double* rx, const double* rz, const double* rp) const {
  auto m = static_cast<FixedStepMemory*>(mem);
  // Add impulse to backward parameters
  casadi_axpy(nrp_, 1., rp, m->rp);

  // Add impulse to state
  casadi_axpy(nrx_, 1., rx, m->rx);

  // Add impulse to backwards dependent variables
  casadi_axpy(nrz_, 1., rz, m->rv + nrv_ - nrz_);
}

ImplicitFixedStepIntegrator::ImplicitFixedStepIntegrator(
    const std::string& name, const Function& dae, double t0, const std::vector<double>& tout)
    : FixedStepIntegrator(name, dae, t0, tout) {
}

ImplicitFixedStepIntegrator::~ImplicitFixedStepIntegrator() {
}

const Options ImplicitFixedStepIntegrator::options_
= {{&FixedStepIntegrator::options_},
    {{"rootfinder",
      {OT_STRING,
      "An implicit function solver"}},
    {"rootfinder_options",
      {OT_DICT,
      "Options to be passed to the NLP Solver"}}
    }
};

void ImplicitFixedStepIntegrator::init(const Dict& opts) {
  // Call the base class init
  FixedStepIntegrator::init(opts);

  // Default (temporary) options
  std::string implicit_function_name = "newton";
  Dict rootfinder_options;

  // Read options
  for (auto&& op : opts) {
    if (op.first=="rootfinder") {
      implicit_function_name = op.second.to_string();
    } else if (op.first=="rootfinder_options") {
      rootfinder_options = op.second;
    }
  }

  // Complete rootfinder dictionary
  rootfinder_options["implicit_input"] = STEP_V0;
  rootfinder_options["implicit_output"] = STEP_VF;

  // Allocate a solver
  Function rf = rootfinder("step", implicit_function_name,
    get_function("implicit_step"), rootfinder_options);
  set_function(rf);
  if (nfwd_ > 0) set_function(rf.forward(nfwd_));

  // Backward integration
  if (nadj_ > 0) {
    Function adj_F = rf.reverse(nadj_);
    set_function(adj_F, adj_F.name(), true);
    if (nfwd_ > 0) {
      create_forward(adj_F.name(), nfwd_);
    }
  }
}

template<typename XType>
Function Integrator::map2oracle(const std::string& name,
    const std::map<std::string, XType>& d) {
  std::vector<XType> de_in(DYN_NUM_IN), de_out(DYN_NUM_OUT);
  for (auto&& i : d) {
    if (i.first=="t") {
      de_in[DYN_T]=i.second;
    } else if (i.first=="x") {
      de_in[DYN_X]=i.second;
    } else if (i.first=="z") {
      de_in[DYN_Z]=i.second;
    } else if (i.first=="p") {
      de_in[DYN_P]=i.second;
    } else if (i.first=="u") {
      de_in[DYN_U]=i.second;
    } else if (i.first=="ode") {
      de_out[DYN_ODE]=i.second;
    } else if (i.first=="alg") {
      de_out[DYN_ALG]=i.second;
    } else if (i.first=="quad") {
      de_out[DYN_QUAD]=i.second;
    } else if (i.first=="zero") {
      de_out[DYN_ZERO]=i.second;
    } else {
      casadi_error("No such field: " + i.first);
    }
  }

  // Consistency checks, input sparsities
  for (casadi_int i = 0; i < DYN_NUM_IN; ++i) {
    const Sparsity& sp = de_in[i].sparsity();
    if (i == DYN_T) {
      casadi_assert(sp.is_empty() || sp.is_scalar(), "DAE time variable must be empty or scalar. "
        "Got dimension " + str(sp.size()));
    } else {
      casadi_assert(sp.is_empty() || sp.is_vector(), "DAE inputs must be empty or vectors. "
        + dyn_in(i) + " has dimension " + str(sp.size()) + ".");
    }
    casadi_assert(sp.is_dense(), "DAE inputs must be dense . "
      + dyn_in(i) + " is sparse.");
    // Convert row vectors to column vectors
    de_in[i] = vec(de_in[i]);
  }

  // Consistency checks, output sparsities
  for (casadi_int i = 0; i < DYN_NUM_OUT; ++i) {
    const Sparsity& sp = de_out[i].sparsity();
    casadi_assert(sp.is_empty() || sp.is_vector(), "DAE outputs must be empty or vectors. "
      + dyn_out(i) + " has dimension " + str(sp.size()));
    // Make sure dense and vector
    de_out[i] = vec(densify(de_out[i]));
  }

  // Construct
  return Function(name, de_in, de_out, dyn_in(), dyn_out());
}

void Integrator::serialize_body(SerializingStream &s) const {
  OracleFunction::serialize_body(s);

  s.version("Integrator", 3);

  s.pack("Integrator::sp_jac_dae", sp_jac_dae_);
  s.pack("Integrator::sp_jac_rdae", sp_jac_rdae_);
  s.pack("Integrator::t0", t0_);
  s.pack("Integrator::tout", tout_);
  s.pack("Integrator::nfwd", nfwd_);
  s.pack("Integrator::nadj", nadj_);
  s.pack("Integrator::rdae", rdae_);

  s.pack("Integrator::nx", nx_);
  s.pack("Integrator::nz", nz_);
  s.pack("Integrator::nq", nq_);
  s.pack("Integrator::nx1", nx1_);
  s.pack("Integrator::nz1", nz1_);
  s.pack("Integrator::nq1", nq1_);

  s.pack("Integrator::nrx", nrx_);
  s.pack("Integrator::nrz", nrz_);
  s.pack("Integrator::nrq", nrq_);
  s.pack("Integrator::nuq", nuq_);
  s.pack("Integrator::nrx1", nrx1_);
  s.pack("Integrator::nrz1", nrz1_);
  s.pack("Integrator::nrq1", nrq1_);
  s.pack("Integrator::nuq1", nuq1_);

  s.pack("Integrator::np", np_);
  s.pack("Integrator::nrp", nrp_);
  s.pack("Integrator::np1", np1_);
  s.pack("Integrator::nrp1", nrp1_);

  s.pack("Integrator::nu", nu_);
  s.pack("Integrator::nu1", nu1_);

  s.pack("Integrator::nom_x", nom_x_);
  s.pack("Integrator::nom_z", nom_z_);

  s.pack("Integrator::augmented_options", augmented_options_);
  s.pack("Integrator::opts", opts_);
  s.pack("Integrator::print_stats", print_stats_);
}

void Integrator::serialize_type(SerializingStream &s) const {
  OracleFunction::serialize_type(s);
  PluginInterface<Integrator>::serialize_type(s);
}

ProtoFunction* Integrator::deserialize(DeserializingStream& s) {
  return PluginInterface<Integrator>::deserialize(s);
}

Integrator::Integrator(DeserializingStream & s) : OracleFunction(s) {
  s.version("Integrator", 3);

  s.unpack("Integrator::sp_jac_dae", sp_jac_dae_);
  s.unpack("Integrator::sp_jac_rdae", sp_jac_rdae_);
  s.unpack("Integrator::t0", t0_);
  s.unpack("Integrator::tout", tout_);
  s.unpack("Integrator::nfwd", nfwd_);
  s.unpack("Integrator::nadj", nadj_);
  s.unpack("Integrator::rdae", rdae_);

  s.unpack("Integrator::nx", nx_);
  s.unpack("Integrator::nz", nz_);
  s.unpack("Integrator::nq", nq_);
  s.unpack("Integrator::nx1", nx1_);
  s.unpack("Integrator::nz1", nz1_);
  s.unpack("Integrator::nq1", nq1_);

  s.unpack("Integrator::nrx", nrx_);
  s.unpack("Integrator::nrz", nrz_);
  s.unpack("Integrator::nrq", nrq_);
  s.unpack("Integrator::nuq", nuq_);
  s.unpack("Integrator::nrx1", nrx1_);
  s.unpack("Integrator::nrz1", nrz1_);
  s.unpack("Integrator::nrq1", nrq1_);
  s.unpack("Integrator::nuq1", nuq1_);

  s.unpack("Integrator::np", np_);
  s.unpack("Integrator::nrp", nrp_);
  s.unpack("Integrator::np1", np1_);
  s.unpack("Integrator::nrp1", nrp1_);

  s.unpack("Integrator::nu", nu_);
  s.unpack("Integrator::nu1", nu1_);

  s.unpack("Integrator::nom_x", nom_x_);
  s.unpack("Integrator::nom_z", nom_z_);

  s.unpack("Integrator::augmented_options", augmented_options_);
  s.unpack("Integrator::opts", opts_);
  s.unpack("Integrator::print_stats", print_stats_);
}

void FixedStepIntegrator::serialize_body(SerializingStream &s) const {
  Integrator::serialize_body(s);

  s.version("FixedStepIntegrator", 3);
  s.pack("FixedStepIntegrator::nk_target", nk_target_);
  s.pack("FixedStepIntegrator::disc", disc_);
  s.pack("FixedStepIntegrator::nv", nv_);
  s.pack("FixedStepIntegrator::nv1", nv1_);
  s.pack("FixedStepIntegrator::nrv", nrv_);
  s.pack("FixedStepIntegrator::nrv1", nrv1_);
}

FixedStepIntegrator::FixedStepIntegrator(DeserializingStream & s) : Integrator(s) {
  s.version("FixedStepIntegrator", 3);
  s.unpack("FixedStepIntegrator::nk_target", nk_target_);
  s.unpack("FixedStepIntegrator::disc", disc_);
  s.unpack("FixedStepIntegrator::nv", nv_);
  s.unpack("FixedStepIntegrator::nv1", nv1_);
  s.unpack("FixedStepIntegrator::nrv", nrv_);
  s.unpack("FixedStepIntegrator::nrv1", nrv1_);
}

void ImplicitFixedStepIntegrator::serialize_body(SerializingStream &s) const {
  FixedStepIntegrator::serialize_body(s);

  s.version("ImplicitFixedStepIntegrator", 2);
}

ImplicitFixedStepIntegrator::ImplicitFixedStepIntegrator(DeserializingStream & s) :
    FixedStepIntegrator(s) {
  s.version("ImplicitFixedStepIntegrator", 2);
}

void Integrator::reset(IntegratorMemory* m, const double* x, const double* z,
    const double* p) const {
  // Update the state
  casadi_copy(x, nx_, m->x);
  casadi_copy(z, nz_, m->z);

  // Set parameters
  casadi_copy(p, np_, m->p);
}

casadi_int Integrator::next_stop(casadi_int k, const double* u) const {
  // Integrate till the end if no input signals
  if (nu_ == 0 || u == 0) return nt() - 1;
  // Find the next discontinuity, if any
  for (; k + 1 < nt(); ++k) {
    // Next control value
    const double *u_next = u + nu_;
    // Check if there is any change in input from k to k + 1
    for (casadi_int i = 0; i < nu_; ++i) {
      // Step change detected: stop integration at k
      if (u[i] != u_next[i]) return k;
    }
    // Shift u
    u = u_next;
  }
  // No step changes detected
  return k;
}

casadi_int Integrator::next_stopB(casadi_int k, const double* u) const {
  // Integrate till the beginning if no input signals
  if (nu_ == 0 || u == 0) return -1;
  // Find the next discontinuity, if any
  for (; k-- > 0; ) {
    // Next control value
    const double *u_next = u - nu_;
    // Check if there is any change in input from k to k + 1
    for (casadi_int i = 0; i < nu_; ++i) {
      // Step change detected: stop integration at k
      if (u[i] != u_next[i]) return k;
    }
    // Shift u
    u = u_next;
  }
  // No step changes detected
  return k;
}

bool Integrator::all_zero(const double* v, casadi_int n) {
  // Quick return if trivially zero
  if (v == 0 || n == 0) return true;
  // Loop over entries
  for (casadi_int i = 0; i < n; ++i) {
    if (v[i] != 0.) return false;
  }
  // All zero if reached here
  return true;
}


} // namespace casadi

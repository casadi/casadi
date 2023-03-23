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
#include "casadi_misc.hpp"

namespace casadi {

std::string to_string(DynIn v) {
  switch (v) {
  case DYN_T: return "t";
  case DYN_X: return "x";
  case DYN_Z: return "z";
  case DYN_P: return "p";
  case DYN_U: return "u";
  case DYN_RX: return "rx";
  case DYN_RZ: return "rz";
  case DYN_RP: return "rp";
  default: break;
  }
  return "";
}

std::string to_string(DynOut v) {
  switch (v) {
  case DYN_ODE: return "ode";
  case DYN_ALG: return "alg";
  case DYN_QUAD: return "quad";
  case DYN_RODE: return "rode";
  case DYN_RALG: return "ralg";
  case DYN_RQUAD: return "rquad";
  case DYN_UQUAD: return "uquad";
  default: break;
  }
  return "";
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
  // Create function oracle and backwards DAE, if any
  Function oracle, rdae;
  oracle = Integrator::map2oracle("dae", dae, &rdae);
  // Create integrator instance
  if (rdae.is_null()) {
    return integrator(name, solver, oracle, t0, tout, opts);
  } else {
    Dict opts2 = opts;
    opts2["rdae"] = rdae;
    return integrator(name, solver, oracle, t0, tout, opts2);
  }
}

Function integrator(const std::string& name, const std::string& solver,
    const MXDict& dae, double t0, const std::vector<double>& tout, const Dict& opts) {
  // Create function oracle and backwards DAE, if any
  Function oracle, rdae;
  oracle = Integrator::map2oracle("dae", dae, &rdae);
  // Create integrator instance
  if (rdae.is_null()) {
    return integrator(name, solver, oracle, t0, tout, opts);
  } else {
    Dict opts2 = opts;
    opts2["rdae"] = rdae;
    return integrator(name, solver, oracle, t0, tout, opts2);
  }
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
  case INTEGRATOR_P:   return "p";
  case INTEGRATOR_U:   return "u";
  case INTEGRATOR_Z0:  return "z0";
  case INTEGRATOR_RX0: return "rx0";
  case INTEGRATOR_RP:  return "rp";
  case INTEGRATOR_RZ0: return "rz0";
  case INTEGRATOR_NUM_IN: break;
  }
  return std::string();
}

std::string integrator_out(casadi_int ind) {
  switch (static_cast<IntegratorOutput>(ind)) {
  case INTEGRATOR_XF:  return "xf";
  case INTEGRATOR_QF:  return "qf";
  case INTEGRATOR_ZF:  return "zf";
  case INTEGRATOR_RXF: return "rxf";
  case INTEGRATOR_RQF: return "rqf";
  case INTEGRATOR_RZF: return "rzf";
  case INTEGRATOR_UQF: return "uqf";
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
  print_stats_ = false;
}

Integrator::~Integrator() {
}

Sparsity Integrator::get_sparsity_in(casadi_int i) {
  switch (static_cast<IntegratorInput>(i)) {
  case INTEGRATOR_X0: return repmat(oracle_.sparsity_in(DYN_X), 1, 1 + nfwd_);
  case INTEGRATOR_P: return repmat(oracle_.sparsity_in(DYN_P), 1, 1 + nfwd_);
  case INTEGRATOR_U: return repmat(oracle_.sparsity_in(DYN_U), 1, nt() * (1 + nfwd_));
  case INTEGRATOR_Z0: return repmat(oracle_.sparsity_in(DYN_Z), 1, 1 + nfwd_);
  case INTEGRATOR_RX0: return repmat(oracle_.sparsity_in(DYN_RX), 1, nt() * (1 + nfwd_));
  case INTEGRATOR_RP: return repmat(oracle_.sparsity_in(DYN_RP), 1, nt() * (1 + nfwd_));
  case INTEGRATOR_RZ0: return repmat(oracle_.sparsity_in(DYN_RZ), 1, nt() * (1 + nfwd_));
  case INTEGRATOR_NUM_IN: break;
  }
  return Sparsity();
}

Sparsity Integrator::get_sparsity_out(casadi_int i) {
  switch (static_cast<IntegratorOutput>(i)) {
  case INTEGRATOR_XF: return repmat(oracle_.sparsity_in(DYN_X), 1, nt() * (1 + nfwd_));
  case INTEGRATOR_QF: return repmat(oracle_.sparsity_out(DYN_QUAD), 1, nt() * (1 + nfwd_));
  case INTEGRATOR_ZF: return repmat(oracle_.sparsity_in(DYN_Z), 1, nt() * (1 + nfwd_));
  case INTEGRATOR_RXF: return repmat(oracle_.sparsity_in(DYN_RX), 1, 1 + nfwd_);
  case INTEGRATOR_RQF: return repmat(oracle_.sparsity_out(DYN_RQUAD), 1, 1 + nfwd_);
  case INTEGRATOR_RZF: return repmat(oracle_.sparsity_in(DYN_RZ), 1, 1 + nfwd_);
  case INTEGRATOR_UQF: return repmat(oracle_.sparsity_out(DYN_UQUAD), 1, nt() * (1 + nfwd_));
  case INTEGRATOR_NUM_OUT: break;
  }
  return Sparsity();
}

bool Integrator::grid_in(casadi_int i) {
  switch (static_cast<IntegratorInput>(i)) {
    case INTEGRATOR_U:
    case INTEGRATOR_RX0:
    case INTEGRATOR_RP:
    case INTEGRATOR_RZ0:
      return true;
    default: break;
  }
  return false;
}

bool Integrator::grid_out(casadi_int i) {
  switch (static_cast<IntegratorOutput>(i)) {
    case INTEGRATOR_XF:
    case INTEGRATOR_QF:
    case INTEGRATOR_ZF:
    case INTEGRATOR_UQF:
      return true;
    default: break;
  }
  return false;
}

casadi_int Integrator::adjmap_in(casadi_int i) {
  switch (static_cast<IntegratorOutput>(i)) {
    case INTEGRATOR_XF: return INTEGRATOR_RX0;
    case INTEGRATOR_QF: return INTEGRATOR_RP;
    case INTEGRATOR_ZF: return INTEGRATOR_RZ0;
    case INTEGRATOR_RXF: return INTEGRATOR_X0;
    case INTEGRATOR_RQF: return INTEGRATOR_P;
    case INTEGRATOR_RZF: return INTEGRATOR_Z0;
    case INTEGRATOR_UQF: return INTEGRATOR_U;
    default: break;
  }
  return -1;
}

casadi_int Integrator::adjmap_out(casadi_int i) {
  switch (static_cast<IntegratorInput>(i)) {
    case INTEGRATOR_X0: return INTEGRATOR_RXF;
    case INTEGRATOR_P: return INTEGRATOR_RQF;
    case INTEGRATOR_U: return INTEGRATOR_UQF;
    case INTEGRATOR_Z0: return INTEGRATOR_RZF;
    case INTEGRATOR_RX0: return INTEGRATOR_XF;
    case INTEGRATOR_RP: return INTEGRATOR_QF;
    case INTEGRATOR_RZ0: return INTEGRATOR_ZF;
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
  double* uq = res[INTEGRATOR_UQF];
  res += INTEGRATOR_NUM_OUT;

  // Setup memory object
  setup(m, arg, res, iw, w);

  // Reset solver, take time to t0
  m->t = t0_;
  reset(m, x0, z0, p);

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
    // Take rx0, rz0, rp past the last grid point
    if (rx0) rx0 += nrx_ * nt();
    if (rz0) rz0 += nrz_ * nt();
    if (rp) rp += nrp_ * nt();
    if (uq) uq += nuq_ * nt();
    // Next stop time due to step change in input
    k_stop = nt();
    // Integrate backward
    for (m->k = nt(); m->k-- > 0; ) {
      m->t = tout_[m->k];
      // Reset the solver, add impulse to backwards integration
      if (rx0) rx0 -= nrx_;
      if (rz0) rz0 -= nrz_;
      if (rp) rp -= nrp_;
      if (uq) uq -= nuq_;
      if (u) u -= nu_;
      if (m->k == nt() - 1) {
        resetB(m, rx0, rz0, rp);
      } else {
        impulseB(m, rx0, rz0, rp);
      }
      // Next output time, or beginning
      casadi_int k_next = m->k - 1;
      m->t_next = k_next < 0 ? t0_ : tout_[k_next];
      // Update integrator stopping time
      if (k_next < k_stop) k_stop = next_stopB(m->k, u);
      m->t_stop = k_stop < 0 ? t0_ : tout_[k_stop];
      // Proceed to the previous time point or t0
      if (verbose_) casadi_message("Integrating backward from output time " + str(m->k)
        + ": t_next = " + str(m->t_next) + ", t_stop = " + str(m->t_stop));
      if (m->k > 0) {
        retreat(m, u, 0, 0, 0, uq);
      } else {
        retreat(m, u, rx, rz, rq, uq);
      }
    }
    // uq should contain the contribution from the grid point, not cumulative
    if (uq) {
      for (m->k = 0; m->k < nt() - 1; ++m->k) {
        casadi_axpy(nuq_, -1., uq + nuq_, uq);
        uq += nuq_;
      }
    }
  }

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
    {"rdae",
      {OT_FUNCTION,
      "Function for evaluating the backwards DAE (transitional option - to be removed)"}},
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
  double t0=0, tf=1;
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
    } else if (op.first=="rdae") {
      rdae_ = op.second;
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
    // Deprecation warning
    casadi_warning("The options 't0', 'tf', 'grid' and 'output_t0' have been deprecated. "
      "Set the time grid by proving additional argument to the 'integrator' call instead.");

    // If grid unset, default to [t0, tf]
    if (grid.empty()) grid = {t0, tf};

    // Construct t0 and tout from grid and output_t0
    t0_ = grid.front();
    tout_ = grid;
    if (!output_t0) tout_.erase(tout_.begin());
  }

  // Call the base class method
  OracleFunction::init(opts);

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

  // Get dimensions (excluding sensitivity equations)
  nx1_ = oracle_.numel_in(DYN_X);
  nz1_ = oracle_.numel_in(DYN_Z);
  nq1_ = oracle_.numel_out(DYN_QUAD);
  np1_ = oracle_.numel_in(DYN_P);
  nu1_ = oracle_.numel_in(DYN_U);
  nrx1_ = oracle_.numel_in(DYN_RX);
  nrz1_ = oracle_.numel_in(DYN_RZ);
  nrp1_ = oracle_.numel_in(DYN_RP);
  nrq1_ = oracle_.numel_out(DYN_RQUAD);
  nuq1_ = oracle_.numel_out(DYN_UQUAD);

  // Get/check dimensions in backwards DAE
  if (rdae_.is_null()) {
    casadi_assert_dev(nrx1_ == 0);
  } else {
    casadi_assert_dev(nrx1_ > 0);
    casadi_assert(nx1_ == rdae_.numel_in(BDYN_X), "Dimension mismatch");
    casadi_assert(nz1_ == rdae_.numel_in(BDYN_Z), "Dimension mismatch");
    casadi_assert(np1_ == rdae_.numel_in(BDYN_P), "Dimension mismatch");
    casadi_assert(nu1_ == rdae_.numel_in(BDYN_U), "Dimension mismatch");
    casadi_assert(nrx1_ == rdae_.numel_in(BDYN_RX), "Dimension mismatch");
    casadi_assert(nrz1_ == rdae_.numel_in(BDYN_RZ), "Dimension mismatch");
    casadi_assert(nrp1_ == rdae_.numel_in(BDYN_RP), "Dimension mismatch");
    casadi_assert(oracle_.numel_out(DYN_RODE) == rdae_.numel_out(BDYN_RODE), "Dimension mismatch");
    casadi_assert(oracle_.numel_out(DYN_RALG) == rdae_.numel_out(BDYN_RALG), "Dimension mismatch");
    casadi_assert(nrq1_ == rdae_.numel_out(BDYN_RQUAD), "Dimension mismatch");
    casadi_assert(nuq1_ == rdae_.numel_out(BDYN_UQUAD), "Dimension mismatch");
  }

  // Consistency checks
  casadi_assert(nx1_ > 0, "Ill-posed ODE - no state");
  casadi_assert(nx1_ == oracle_.numel_out(DYN_ODE), "Dimension mismatch for 'ode'");
  casadi_assert(nz1_ == oracle_.numel_out(DYN_ALG), "Dimension mismatch for 'alg'");
  casadi_assert(nrx1_ == oracle_.numel_out(DYN_RODE), "Dimension mismatch for 'rode'");
  casadi_assert(nrz1_ == oracle_.numel_out(DYN_RALG), "Dimension mismatch for 'ralg'");

  // Get dimensions (including sensitivity equations)
  nx_ = nx1_ * (1 + nfwd_);
  nz_ = nz1_ * (1 + nfwd_);
  nq_ = nq1_ * (1 + nfwd_);
  np_ = np1_ * (1 + nfwd_);
  nu_ = nu1_ * (1 + nfwd_);
  nrx_ = nrx1_ * (1 + nfwd_);
  nrz_ = nrz1_ * (1 + nfwd_);
  nrp_ = nrp1_ * (1 + nfwd_);
  nrq_ = nrq1_ * (1 + nfwd_);
  nuq_ = nuq1_ * (1 + nfwd_);

  // Create problem functions, forward problem
  create_function("daeF", fdyn_in(), fdae_out());
  create_function("quadF", fdyn_in(), fquad_out());
  if (nfwd_ > 0) {
    // one direction to conserve memory, symbolic processing time
    create_forward("daeF", 1);
    create_forward("quadF", 1);
  }

  // Create problem functions, backward problem
  if (nrx1_ > 0) {
    create_function(rdae_, "daeB", bdyn_in(), bdae_out());
    create_function(rdae_, "quadB", bdyn_in(), bquad_out());
    if (nfwd_ > 0) {
      // one direction to conserve memory, symbolic processing time
      create_forward("daeB", 1);
      create_forward("quadB", 1);
    }
  }

  // Get the sparsities of the forward and reverse DAE
  sp_jac_dae_ = sp_jac_dae();
  casadi_assert(!sp_jac_dae_.is_singular(),
    "Jacobian of the forward problem is structurally rank-deficient. "
    "sprank(J)=" + str(sprank(sp_jac_dae_)) + "<" + str(nx_+nz_));
  if (nrx_ > 0) {
    sp_jac_rdae_ = sp_jac_rdae();
    casadi_assert(!sp_jac_rdae_.is_singular(),
      "Jacobian of the backward problem is structurally rank-deficient. "
      "sprank(J)=" + str(sprank(sp_jac_rdae_)) + "<" + str(nrx_+nrz_));
  }

  // Work vectors for sparsity pattern propagation: Can be reused in derived classes
  alloc_w(nx_, true); // x
  alloc_w(nz_, true); // z
  alloc_w(nx_, true); // x_prev
  alloc_w(nrx_, true); // rx
  alloc_w(nrz_, true); // rz
  alloc_w(nrx_, true); // rx_prev
  alloc_w(nrq_, true); // rq
  alloc_w(nx_+nz_);  // Sparsity::sp_solve
  alloc_w(nrx_+nrz_);  // Sparsity::sp_solve
}

int Integrator::init_mem(void* mem) const {
  if (OracleFunction::init_mem(mem)) return 1;

  //auto m = static_cast<IntegratorMemory*>(mem);
  return 0;
}

Function Integrator::augmented_dae(Function* rdae) const {
  // If no sensitivities, augmented oracle is the oracle itself
  if (nfwd_ == 0) {
    if (rdae) *rdae = rdae_;
    return oracle_;
  }
  // Name of augmented DAE
  std::string aug_name = "fsens" + str(nfwd_) + "_" + oracle_.name();
  // Use function in cache, if available
  Function ret;
  // if (incache(aug_name, ret)) return ret;  // caching disabled while implementing #3047
  // Create new augmented oracle
  try {
    if (oracle_.is_a("SXFunction")) {
      ret = get_augmented_dae<SX>(aug_name, rdae);
    } else {
      ret = get_augmented_dae<MX>(aug_name, rdae);
    }
  } catch (std::exception& e) {
    casadi_error("Failed to generate augmented DAE for " + name_ + ":\n" + e.what());
  }
  // Save to Function cache and return
  // tocache(ret);  // caching disabled while implementing #3047
  return ret;
}

template<typename MatType>
Function Integrator::get_augmented_dae(const std::string& name, Function* rdae) const {
  if (verbose_) casadi_message(name_ + "::get_augmented_dae");

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

  // Create backwards DAE
  if (rdae != 0 && !arg[DYN_RX].is_empty()) {
    // Sort expressions by name
    std::map<std::string, MatType> r;
    for (casadi_int i = 0; i < DYN_NUM_IN; ++i) r[dyn_in(i)] = arg[i];
    for (casadi_int i = 0; i < DYN_NUM_OUT; ++i) r[dyn_out(i)] = res[i];
    // Construct backwards DAE
    std::vector<MatType> rdae_in, rdae_out;
    for (auto& n : bdyn_in()) rdae_in.push_back(r.at(n));
    for (auto& n : bdyn_out()) rdae_out.push_back(r.at(n));
    *rdae = Function("rdae", rdae_in, rdae_out, bdyn_in(), bdyn_out());
  }

  // Convert to oracle function and return
  return Function(name, arg, res, dyn_in(), dyn_out());
}

template<typename MatType>
Function Integrator::aug_adj(const Function& this_dae, const Function& this_rdae,
    casadi_int nadj, Function* aug_rdae) const {
  if (verbose_) casadi_message(name_ + "::aug_adj");

  // Get input and output expressions
  std::vector<MatType> arg = MatType::get_input(this_dae);
  std::vector<MatType> res = this_dae(arg);

  // Symbolic expression for augmented DAE
  std::vector<std::vector<MatType>> aug_in(DYN_NUM_IN);
  for (casadi_int i = 0; i < DYN_NUM_IN; ++i) aug_in[i].push_back(arg.at(i));
  std::vector<std::vector<MatType>> aug_out(DYN_NUM_OUT);
  for (casadi_int i = 0; i < DYN_NUM_OUT; ++i) aug_out[i].push_back(res.at(i));

  // Zero of time dimension
  MatType zero_t = MatType::zeros(this_dae.sparsity_in(DYN_T));

  // Reverse mode directional derivatives
  std::vector<std::vector<MatType>> seed(nadj, std::vector<MatType>(DYN_NUM_OUT));
  for (casadi_int d = 0; d < nadj; ++d) {
    std::string pref = "aug" + str(d) + "_";
    for (casadi_int i = 0; i < DYN_NUM_OUT; ++i) {
      seed[d][i] = MatType::sym(pref + dyn_out(i), this_dae.sparsity_out(i));
    }
    aug_in[DYN_RX].push_back(seed[d][DYN_ODE]);
    aug_in[DYN_RZ].push_back(seed[d][DYN_ALG]);
    aug_in[DYN_RP].push_back(seed[d][DYN_QUAD]);
    aug_in[DYN_X].push_back(seed[d][DYN_RODE]);
    aug_in[DYN_Z].push_back(seed[d][DYN_RALG]);
    aug_in[DYN_P].push_back(seed[d][DYN_RQUAD]);
    aug_in[DYN_U].push_back(seed[d][DYN_UQUAD]);
  }

  // Calculate directional derivatives
  std::vector<std::vector<MatType>> sens;
  bool always_inline = this_dae.is_a("SXFunction") || this_dae.is_a("MXFunction");
  this_dae->call_reverse(arg, res, seed, sens, always_inline, false);

  // Collect sensitivity equations
  casadi_assert_dev(sens.size()==nadj);
  for (casadi_int d = 0; d < nadj; ++d) {
    casadi_assert_dev(sens[d].size() == DYN_NUM_IN);
    aug_out[DYN_RODE].push_back(project(sens[d][DYN_X], this_dae.sparsity_in(DYN_X)));
    aug_out[DYN_RALG].push_back(project(sens[d][DYN_Z], this_dae.sparsity_in(DYN_Z)));
    aug_out[DYN_RQUAD].push_back(project(sens[d][DYN_P], this_dae.sparsity_in(DYN_P)));
    aug_out[DYN_UQUAD].push_back(project(sens[d][DYN_U], this_dae.sparsity_in(DYN_U)));
    aug_out[DYN_ODE].push_back(project(sens[d][DYN_RX], this_dae.sparsity_in(DYN_RX)));
    aug_out[DYN_ALG].push_back(project(sens[d][DYN_RZ], this_dae.sparsity_in(DYN_RZ)));
    aug_out[DYN_QUAD].push_back(project(sens[d][DYN_RP], this_dae.sparsity_in(DYN_RP)));
  }

  // Concatenate expressions
  for (casadi_int i = 0; i < DYN_NUM_IN; ++i) arg.at(i) = vertcat(aug_in[i]);
  for (casadi_int i = 0; i < DYN_NUM_OUT; ++i) res.at(i) = vertcat(aug_out[i]);

  // Hack: Make sure that forward problem does not depend on backward states
  // Should become unnecessary once #3047 is complete
  std::map<std::string, MatType> r;
  for (casadi_int i = 0; i < DYN_NUM_IN; ++i) r[dyn_in(i)] = arg[i];
  for (casadi_int i = 0; i < DYN_NUM_OUT; ++i) r[dyn_out(i)] = res[i];
  Function f("f", {r["t"], r["x"], r["z"], r["p"], r["u"]},
                  {r["ode"], r["alg"], r["quad"]}, {{"allow_free", true}});
  if (f.has_free()) {
    // Replace dependencies of rx, rz and rp with zeros
    f = Function("f", {r["t"], r["x"], r["z"], r["p"], r["u"],
                        r["rx"], r["rz"], r["rp"]},
                      {r["ode"], r["alg"], r["quad"]});
    std::vector<MatType> v = {r["t"], r["x"], r["z"], r["p"], r["u"], 0, 0, 0};
    v = f(v);
    res[DYN_ODE] = v.at(0);
    res[DYN_ALG] = v.at(1);
    res[DYN_QUAD] = v.at(2);
  }

  // Create backwards DAE
  if (aug_rdae != 0 && !arg[DYN_RX].is_empty()) {
    std::vector<MatType> rdae_in, rdae_out;
    for (auto& n : bdyn_in()) rdae_in.push_back(r.at(n));
    for (auto& n : bdyn_out()) rdae_out.push_back(r.at(n));
    *aug_rdae = Function("rdae", rdae_in, rdae_out, bdyn_in(), bdyn_out());
  }

  // Convert to oracle function and return
  std::string aug_prefix = "asens" + str(nadj) + "_";
  std::string dae_name = aug_prefix + this_dae.name();
  return Function(dae_name, arg, res, dyn_in(), dyn_out());
}

int Integrator::fdae_sp_forward(SpForwardMem* m, const bvec_t* x,
    const bvec_t* p, const bvec_t* u, bvec_t* ode, bvec_t* alg) const {
  // Evaluate nondifferentiated
  m->arg[FDYN_T] = nullptr;  // t
  m->arg[FDYN_X] = x;  // x
  m->arg[FDYN_Z] = nullptr;  // z
  m->arg[FDYN_P] = p;  // p
  m->arg[FDYN_U] = u;  // u
  m->res[FDAE_ODE] = ode;  // ode
  m->res[FDAE_ALG] = alg;  // alg
  if (calc_sp_forward("daeF", m->arg, m->res, m->iw, m->w)) return 1;
  // Evaluate sensitivities
  for (casadi_int i = 0; i < nfwd_; ++i) {
    m->arg[FDYN_NUM_IN + FDAE_ODE] = ode;  // out:ode
    m->arg[FDYN_NUM_IN + FDAE_ALG] = alg;  // out:alg
    m->arg[FDYN_NUM_IN + FDAE_NUM_OUT + FDYN_T] = nullptr;  // fwd:t
    m->arg[FDYN_NUM_IN + FDAE_NUM_OUT + FDYN_X] = x + (i + 1) * nx1_;  // fwd:x
    m->arg[FDYN_NUM_IN + FDAE_NUM_OUT + FDYN_Z] = nullptr;  // fwd:z
    m->arg[FDYN_NUM_IN + FDAE_NUM_OUT + FDYN_P] = p + (i + 1) * np1_;  // fwd:p
    m->arg[FDYN_NUM_IN + FDAE_NUM_OUT + FDYN_U] = u + (i + 1) * nu1_;  // fwd:u
    m->res[FDAE_ODE] = ode + (i + 1) * nx1_;  // fwd:ode
    m->res[FDAE_ALG] = alg + (i + 1) * nz1_;  // fwd:alg
    if (calc_sp_forward(forward_name("daeF", 1), m->arg, m->res, m->iw, m->w)) return 1;
  }
  return 0;
}

int Integrator::fquad_sp_forward(SpForwardMem* m, const bvec_t* x, const bvec_t* z,
    const bvec_t* p, const bvec_t* u, bvec_t* quad) const {
  // Evaluate nondifferentiated
  m->arg[FDYN_T] = nullptr;  // t
  m->arg[FDYN_X] = x;  // x
  m->arg[FDYN_Z] = z;  // z
  m->arg[FDYN_P] = p;  // p
  m->arg[FDYN_U] = u;  // u
  m->res[FQUAD_QUAD] = quad;  // quad
  if (calc_sp_forward("quadF", m->arg, m->res, m->iw, m->w)) return 1;
  // Evaluate sensitivities
  for (casadi_int i = 0; i < nfwd_; ++i) {
    m->arg[FDYN_NUM_IN + FQUAD_QUAD] = quad;  // out:quad
    m->arg[FDYN_NUM_IN + FQUAD_NUM_OUT + FDYN_T] = nullptr;  // fwd:t
    m->arg[FDYN_NUM_IN + FQUAD_NUM_OUT + FDYN_X] = x + (i + 1) * nx1_;  // fwd:x
    m->arg[FDYN_NUM_IN + FQUAD_NUM_OUT + FDYN_Z] = z + (i + 1) * nz1_;  // fwd:z
    m->arg[FDYN_NUM_IN + FQUAD_NUM_OUT + FDYN_P] = p + (i + 1) * np1_;  // fwd:p
    m->arg[FDYN_NUM_IN + FQUAD_NUM_OUT + FDYN_U] = u + (i + 1) * nu1_;  // fwd:u
    m->res[FQUAD_QUAD] = quad + (i + 1) * nq1_;  // fwd:quad
    if (calc_sp_forward(forward_name("quadF", 1), m->arg, m->res, m->iw, m->w)) return 1;
  }
  return 0;
}

int Integrator::bdae_sp_forward(SpForwardMem* m, const bvec_t* x, const bvec_t* z,
    const bvec_t* p, const bvec_t* u, const bvec_t* rx, const bvec_t* rp,
    bvec_t* rode, bvec_t* ralg) const {
  // Evaluate nondifferentiated
  m->arg[BDYN_T] = nullptr;  // t
  m->arg[BDYN_X] = x;  // x
  m->arg[BDYN_Z] = z;  // z
  m->arg[BDYN_P] = p;  // p
  m->arg[BDYN_U] = u;  // u
  m->arg[BDYN_RX] = rx;  // rx
  m->arg[BDYN_RZ] = nullptr;  // rz
  m->arg[BDYN_RP] = rp;  // rp
  m->res[BDAE_RODE] = rode;  // rode
  m->res[BDAE_RALG] = ralg;  // ralg
  if (calc_sp_forward("daeB", m->arg, m->res, m->iw, m->w)) return 1;
  // Evaluate sensitivities
  for (casadi_int i = 0; i < nfwd_; ++i) {
    m->arg[BDYN_NUM_IN + BDAE_RODE] = rode;  // out:rode
    m->arg[BDYN_NUM_IN + BDAE_RALG] = ralg;  // out:ralg
    m->arg[BDYN_NUM_IN + BDAE_NUM_OUT + BDYN_T] = nullptr;  // fwd:t
    m->arg[BDYN_NUM_IN + BDAE_NUM_OUT + BDYN_X] = x + (i + 1) * nx1_;  // fwd:x
    m->arg[BDYN_NUM_IN + BDAE_NUM_OUT + BDYN_Z] = z + (i + 1) * nz1_;  // fwd:z
    m->arg[BDYN_NUM_IN + BDAE_NUM_OUT + BDYN_P] = p + (i + 1) * np1_;  // fwd:p
    m->arg[BDYN_NUM_IN + BDAE_NUM_OUT + BDYN_U] = u + (i + 1) * nu1_;  // fwd:u
    m->arg[BDYN_NUM_IN + BDAE_NUM_OUT + BDYN_RX] = rx + (i + 1) * nrx1_;  // fwd:rx
    m->arg[BDYN_NUM_IN + BDAE_NUM_OUT + BDYN_RZ] = nullptr;  // fwd:rz
    m->arg[BDYN_NUM_IN + BDAE_NUM_OUT + BDYN_RP] = rp + (i + 1) * nrz1_;  // fwd:rp
    m->res[BDAE_RODE] = rode + (i + 1) * nrx1_;  // fwd:rode
    m->res[BDAE_RALG] = ralg + (i + 1) * nrz1_;  // fwd:ralg
    if (calc_sp_forward(forward_name("daeB", 1), m->arg, m->res, m->iw, m->w)) return 1;
  }
  return 0;
}

int Integrator::bquad_sp_forward(SpForwardMem* m, const bvec_t* x, const bvec_t* z,
    const bvec_t* p, const bvec_t* u, const bvec_t* rx, const bvec_t* rz, const bvec_t* rp,
    bvec_t* rquad, bvec_t* uquad) const {
  // Evaluate nondifferentiated
  m->arg[BDYN_T] = nullptr;  // t
  m->arg[BDYN_X] = x;  // x
  m->arg[BDYN_Z] = z;  // z
  m->arg[BDYN_P] = p;  // p
  m->arg[BDYN_U] = u;  // u
  m->arg[BDYN_RX] = rx;  // rx
  m->arg[BDYN_RZ] = rz;  // rz
  m->arg[BDYN_RP] = rp;  // rp
  m->res[BQUAD_RQUAD] = rquad;  // rquad
  m->res[BQUAD_UQUAD] = uquad;  // uquad
  if (calc_sp_forward("quadB", m->arg, m->res, m->iw, m->w)) return 1;
  // Evaluate sensitivities
  for (casadi_int i = 0; i < nfwd_; ++i) {
    m->arg[BDYN_NUM_IN + BQUAD_RQUAD] = rquad;  // out:rquad
    m->arg[BDYN_NUM_IN + BQUAD_UQUAD] = uquad;  // out:uquad
    m->arg[BDYN_NUM_IN + BQUAD_NUM_OUT + BDYN_T] = nullptr;  // fwd:t
    m->arg[BDYN_NUM_IN + BQUAD_NUM_OUT + BDYN_X] = x + (i + 1) * nx1_;  // fwd:x
    m->arg[BDYN_NUM_IN + BQUAD_NUM_OUT + BDYN_Z] = z + (i + 1) * nz1_;  // fwd:z
    m->arg[BDYN_NUM_IN + BQUAD_NUM_OUT + BDYN_P] = p + (i + 1) * np1_;  // fwd:p
    m->arg[BDYN_NUM_IN + BQUAD_NUM_OUT + BDYN_U] = u + (i + 1) * nu1_;  // fwd:u
    m->arg[BDYN_NUM_IN + BQUAD_NUM_OUT + BDYN_RX] = rx + (i + 1) * nrx1_;  // fwd:rx
    m->arg[BDYN_NUM_IN + BQUAD_NUM_OUT + BDYN_RZ] = rz + (i + 1) * nrz1_;  // fwd:rz
    m->arg[BDYN_NUM_IN + BQUAD_NUM_OUT + BDYN_RP] = rp + (i + 1) * nrp1_;  // fwd:rp
    m->res[BQUAD_RQUAD] = rquad ? rquad + (i + 1) * nrq1_ : 0;  // fwd:rquad
    m->res[BQUAD_UQUAD] = uquad ? uquad + (i + 1) * nuq1_ : 0;  // fwd:uquad
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
  const bvec_t* rx0 = arg[INTEGRATOR_RX0];
  const bvec_t* rp = arg[INTEGRATOR_RP];
  arg += n_in_;

  // Outputs
  bvec_t* xf = res[INTEGRATOR_XF];
  bvec_t* zf = res[INTEGRATOR_ZF];
  bvec_t* qf = res[INTEGRATOR_QF];
  bvec_t* rxf = res[INTEGRATOR_RXF];
  bvec_t* rzf = res[INTEGRATOR_RZF];
  bvec_t* rqf = res[INTEGRATOR_RQF];
  bvec_t* uqf = res[INTEGRATOR_UQF];
  res += n_out_;

  // Work vectors
  bvec_t *x = w; w += nx_;
  bvec_t *z = w; w += nz_;
  bvec_t *x_prev = w; w += nx_;
  bvec_t *rx = w; w += nrx_;
  bvec_t *rz = w; w += nrz_;
  bvec_t *rx_prev = w; w += nrx_;
  bvec_t *rq = w; w += nrq_;

  // Memory struct for function calls below
  SpForwardMem m = {arg, res, iw, w};

  // Copy initial guess to x_prev
  std::copy_n(x0, nx_, x_prev);

  // Propagate forward
  for (casadi_int k = 0; k < nt(); ++k) {
    // Propagate through DAE function
    if (fdae_sp_forward(&m, x_prev, p, u, x, z)) return 1;
    for (casadi_int i = 0; i < nx_; ++i) x[i] |= x_prev[i];

    // "Solve" in order to resolve interdependencies (cf. Rootfinder)
    std::copy_n(x, nx_, w);
    std::copy_n(z, nz_, w + nx_);
    std::fill_n(x, nx_ + nz_, 0);
    sp_jac_dae_.spsolve(x, w, false);

    // Get xf and zf
    if (xf) std::copy_n(x, nx_, xf);
    if (zf) std::copy_n(z, nz_, zf);

    // Propagate to quadratures
    if (nq_ > 0 && qf) {
      if (fquad_sp_forward(&m, x, z, p, u, qf)) return 1;
    }

    // Shift time
    std::copy_n(x, nx_, x_prev);
    if (xf) xf += nx_;
    if (zf) zf += nz_;
    if (qf) qf += nq_;
    if (u) u += nu_;
  }

  if (nrx_ > 0) {
    // Clear rx_prev, rqf
    std::fill_n(rx_prev, nrx_, 0);
    if (rqf) std::fill_n(rqf, nrq_, 0);

    // Take rx0, rp, uqf past the last grid point
    if (rx0) rx0 += nrx_ * nt();
    if (rp) rp += nrp_ * nt();
    if (uqf) uqf += nuq_ * nt();

    // Integrate backward
    for (casadi_int k = nt(); k-- > 0; ) {
      // Shift time
      if (rx0) rx0 -= nrx_;
      if (rp) rp -= nrp_;
      if (uqf) uqf -= nuq_;
      if (u) u -= nu_;

      // Add impulse from rx0
      if (rx0) {
        for (casadi_int i = 0; i < nrx_; ++i) rx_prev[i] |= rx0[i];
      }

      // Propagate through DAE function
      if (bdae_sp_forward(&m, x, z, p, u, rx_prev, rp, rx, rz)) return 1;
      for (casadi_int i = 0; i < nrx_; ++i) rx[i] |= rx_prev[i];

      // "Solve" in order to resolve interdependencies (cf. Rootfinder)
      std::copy_n(rx, nrx_ + nrz_, w);
      std::fill_n(rx, nrx_ + nrz_, 0);
      sp_jac_rdae_.spsolve(rx, w, false);

      // Propagate to quadratures
      if ((nrq_ > 0 && rqf) || (nuq_ > 0 && uqf)) {
        if (bquad_sp_forward(&m, x, z, p, u, rx, rz, rp, rq, uqf)) return 1;
        // Sum contributions to rqf
        if (rqf) {
          for (casadi_int i = 0; i < nrq_; ++i) rqf[i] |= rq[i];
        }
      }

      // Update rx_prev
      std::copy_n(rx, nx_, rx_prev);
    }

    // Get rxf and rzf at initial time
    if (rxf) std::copy_n(rx, nrx_, rxf);
    if (rzf) std::copy_n(rz, nrz_, rzf);
  }
  return 0;
}

int Integrator::fdae_sp_reverse(SpReverseMem* m, bvec_t* x,
    bvec_t* p, bvec_t* u, bvec_t* ode, bvec_t* alg) const {
  // Nondifferentiated inputs
  m->arg[FDYN_T] = nullptr;  // t
  m->arg[FDYN_X] = x;  // x
  m->arg[FDYN_Z] = nullptr;  // z
  m->arg[FDYN_P] = p;  // p
  m->arg[FDYN_U] = u;  // u
  // Propagate through sensitivities
  for (casadi_int i = 0; i < nfwd_; ++i) {
    m->res[FDAE_ODE] = ode + (i + 1) * nx1_;  // fwd:ode
    m->res[FDAE_ALG] = alg + (i + 1) * nz1_;  // fwd:alg
    m->arg[FDYN_NUM_IN + FDAE_ODE] = ode;  // out:ode
    m->arg[FDYN_NUM_IN + FDAE_ALG] = alg;  // out:alg
    m->arg[FDYN_NUM_IN + FDAE_NUM_OUT + FDYN_T] = nullptr;  // fwd:t
    m->arg[FDYN_NUM_IN + FDAE_NUM_OUT + FDYN_X] = x + (i + 1) * nx1_;  // fwd:x
    m->arg[FDYN_NUM_IN + FDAE_NUM_OUT + FDYN_Z] = nullptr;  // fwd:z
    m->arg[FDYN_NUM_IN + FDAE_NUM_OUT + FDYN_P] = p + (i + 1) * np1_;  // fwd:p
    m->arg[FDYN_NUM_IN + FDAE_NUM_OUT + FDYN_U] = u + (i + 1) * nu1_;  // fwd:u
    if (calc_sp_reverse(forward_name("daeF", 1), m->arg, m->res, m->iw, m->w)) return 1;
  }
  // Propagate through nondifferentiated
  m->res[FDAE_ODE] = ode;  // ode
  m->res[FDAE_ALG] = alg;  // alg
  if (calc_sp_reverse("daeF", m->arg, m->res, m->iw, m->w)) return 1;
  return 0;
}

int Integrator::fquad_sp_reverse(SpReverseMem* m, bvec_t* x, bvec_t* z,
    bvec_t* p, bvec_t* u, bvec_t* quad) const {
  // Nondifferentiated inputs
  m->arg[FDYN_T] = nullptr;  // t
  m->arg[FDYN_X] = x;  // x
  m->arg[FDYN_Z] = z;  // z
  m->arg[FDYN_P] = p;  // p
  m->arg[FDYN_U] = u;  // u
  // Propagate through sensitivities
  for (casadi_int i = 0; i < nfwd_; ++i) {
    m->res[FQUAD_QUAD] = quad + (i + 1) * nq1_;  // fwd:quad
    m->arg[FDYN_NUM_IN + FQUAD_QUAD] = quad;  // out:quad
    m->arg[FDYN_NUM_IN + FQUAD_NUM_OUT + FDYN_T] = nullptr;  // fwd:t
    m->arg[FDYN_NUM_IN + FQUAD_NUM_OUT + FDYN_X] = x + (i + 1) * nx1_;  // fwd:x
    m->arg[FDYN_NUM_IN + FQUAD_NUM_OUT + FDYN_Z] = z + (i + 1) * nz1_;  // fwd:z
    m->arg[FDYN_NUM_IN + FQUAD_NUM_OUT + FDYN_P] = p + (i + 1) * np1_;  // fwd:p
    m->arg[FDYN_NUM_IN + FQUAD_NUM_OUT + FDYN_U] = u + (i + 1) * nu1_;  // fwd:u
    if (calc_sp_reverse(forward_name("quadF", 1), m->arg, m->res, m->iw, m->w)) return 1;
  }
  // Propagate through nondifferentiated
  m->res[FQUAD_QUAD] = quad;  // quad
  if (calc_sp_reverse("quadF", m->arg, m->res, m->iw, m->w)) return 1;
  return 0;
}

int Integrator::bdae_sp_reverse(SpReverseMem* m, bvec_t* x, bvec_t* z,
    bvec_t* p, bvec_t* u, bvec_t* rx, bvec_t* rp,
    bvec_t* rode, bvec_t* ralg) const {
  // Nondifferentiated inputs
  m->arg[BDYN_T] = nullptr;  // t
  m->arg[BDYN_X] = x;  // x
  m->arg[BDYN_Z] = z;  // z
  m->arg[BDYN_P] = p;  // p
  m->arg[BDYN_U] = u;  // u
  m->arg[BDYN_RX] = rx;  // rx
  m->arg[BDYN_RZ] = nullptr;  // rz
  m->arg[BDYN_RP] = rp;  // rp
  // Propagate through sensitivities
  for (casadi_int i = 0; i < nfwd_; ++i) {
    m->res[BDAE_RODE] = rode + (i + 1) * nrx1_;  // fwd:rode
    m->res[BDAE_RALG] = ralg + (i + 1) * nrz1_;  // fwd:ralg
    m->arg[BDYN_NUM_IN + BDAE_RODE] = rode;  // out:rode
    m->arg[BDYN_NUM_IN + BDAE_RALG] = ralg;  // out:ralg
    m->arg[BDYN_NUM_IN + BDAE_NUM_OUT + BDYN_T] = nullptr;  // fwd:t
    m->arg[BDYN_NUM_IN + BDAE_NUM_OUT + BDYN_X] = x + (i + 1) * nx1_;  // fwd:x
    m->arg[BDYN_NUM_IN + BDAE_NUM_OUT + BDYN_Z] = z + (i + 1) * nz1_;  // fwd:z
    m->arg[BDYN_NUM_IN + BDAE_NUM_OUT + BDYN_P] = p + (i + 1) * np1_;  // fwd:p
    m->arg[BDYN_NUM_IN + BDAE_NUM_OUT + BDYN_U] = u + (i + 1) * nu1_;  // fwd:u
    m->arg[BDYN_NUM_IN + BDAE_NUM_OUT + BDYN_RX] = rx + (i + 1) * nrx1_;  // fwd:rx
    m->arg[BDYN_NUM_IN + BDAE_NUM_OUT + BDYN_RZ] = nullptr;  // fwd:rz
    m->arg[BDYN_NUM_IN + BDAE_NUM_OUT + BDYN_RP] = rp + (i + 1) * nrz1_;  // fwd:rp
    if (calc_sp_reverse(forward_name("daeB", 1), m->arg, m->res, m->iw, m->w)) return 1;
  }
  // Propagate through nondifferentiated
  m->res[BDAE_RODE] = rode;  // rode
  m->res[BDAE_RALG] = ralg;  // ralg
  if (calc_sp_reverse("daeB", m->arg, m->res, m->iw, m->w)) return 1;
  return 0;
}

int Integrator::bquad_sp_reverse(SpReverseMem* m, bvec_t* x, bvec_t* z,
    bvec_t* p, bvec_t* u, bvec_t* rx, bvec_t* rz, bvec_t* rp,
    bvec_t* rquad, bvec_t* uquad) const {
  // Nondifferentiated inputs
  m->arg[BDYN_T] = nullptr;  // t
  m->arg[BDYN_X] = x;  // x
  m->arg[BDYN_Z] = z;  // z
  m->arg[BDYN_P] = p;  // p
  m->arg[BDYN_U] = u;  // u
  m->arg[BDYN_RX] = rx;  // rx
  m->arg[BDYN_RZ] = rz;  // rz
  m->arg[BDYN_RP] = rp;  // rp
  // Propagate through sensitivities
  for (casadi_int i = 0; i < nfwd_; ++i) {
    m->res[BQUAD_RQUAD] = rquad ? rquad + (i + 1) * nrq1_ : 0;  // fwd:rquad
    m->res[BQUAD_UQUAD] = uquad ? uquad + (i + 1) * nuq1_ : 0;  // fwd:uquad
    m->arg[BDYN_NUM_IN + BQUAD_RQUAD] = rquad;  // out:rquad
    m->arg[BDYN_NUM_IN + BQUAD_UQUAD] = uquad;  // out:uquad
    m->arg[BDYN_NUM_IN + BQUAD_NUM_OUT + BDYN_T] = nullptr;  // fwd:t
    m->arg[BDYN_NUM_IN + BQUAD_NUM_OUT + BDYN_X] = x + (i + 1) * nx1_;  // fwd:x
    m->arg[BDYN_NUM_IN + BQUAD_NUM_OUT + BDYN_Z] = z + (i + 1) * nz1_;  // fwd:z
    m->arg[BDYN_NUM_IN + BQUAD_NUM_OUT + BDYN_P] = p + (i + 1) * np1_;  // fwd:p
    m->arg[BDYN_NUM_IN + BQUAD_NUM_OUT + BDYN_U] = u + (i + 1) * nu1_;  // fwd:u
    m->arg[BDYN_NUM_IN + BQUAD_NUM_OUT + BDYN_RX] = rx + (i + 1) * nrx1_;  // fwd:rx
    m->arg[BDYN_NUM_IN + BQUAD_NUM_OUT + BDYN_RZ] = rz + (i + 1) * nrz1_;  // fwd:rz
    m->arg[BDYN_NUM_IN + BQUAD_NUM_OUT + BDYN_RP] = rp + (i + 1) * nrp1_;  // fwd:rp
    if (calc_sp_reverse(forward_name("quadB", 1), m->arg, m->res, m->iw, m->w)) return 1;
  }
  // Propagate through nondifferentiated
  m->res[BQUAD_RQUAD] = rquad;  // rquad
  m->res[BQUAD_UQUAD] = uquad;  // uquad
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
  bvec_t* rx0 = arg[INTEGRATOR_RX0];
  bvec_t* rp = arg[INTEGRATOR_RP];
  arg += n_in_;

  // Outputs
  bvec_t* xf = res[INTEGRATOR_XF];
  bvec_t* zf = res[INTEGRATOR_ZF];
  bvec_t* qf = res[INTEGRATOR_QF];
  bvec_t* rxf = res[INTEGRATOR_RXF];
  bvec_t* rzf = res[INTEGRATOR_RZF];
  bvec_t* rqf = res[INTEGRATOR_RQF];
  bvec_t* uqf = res[INTEGRATOR_UQF];
  res += n_out_;

  // Work vectors
  bvec_t *x = w; w += nx_;
  bvec_t *z = w; w += nz_;
  bvec_t *x_prev = w; w += nx_;
  bvec_t *rx = w; w += nrx_;
  bvec_t *rz = w; w += nrz_;
  bvec_t *rx_prev = w; w += nrx_;
  bvec_t *rq = w; w += nrq_;

  // Memory struct for function calls below
  SpReverseMem m = {arg, res, iw, w};

  // Clear state vector
  std::fill_n(x, nx_, 0);
  std::fill_n(z, nz_, 0);

  if (nrx_ > 0) {
    // Propagate from rxf and rzf at initial time
    if (rxf) {
      std::copy_n(rxf, nrx_, rx);
      std::fill_n(rxf, nrx_, 0);
    } else {
      std::fill_n(rx, nrx_, 0);
    }
    if (rzf) {
      std::copy_n(rzf, nrz_, rz);
      std::fill_n(rzf, nrz_, 0);
    } else {
      std::fill_n(rz, nrz_, 0);
    }

    // Save rqf: See note below
    if (rqf) std::copy_n(rqf, nrq_, rq);

    // Step backwards through backward problem
    for (casadi_int k = 0; k < nt(); ++k) {
      // Restore rqf: See note below
      if (rqf) std::copy_n(rq, nrq_, rqf);

      // Add impulse from rx0
      if (rx0) {
        for (casadi_int i = 0; i < nrx_; ++i) rx[i] |= rx0[i];
        std::fill_n(rx0, nrx_, 0);
      }

      // Get dependencies from backward quadratures
      if ((nrq_ > 0 && rqf) || (nuq_ > 0 && uqf)) {
        if (bquad_sp_reverse(&m, x, z, p, u, rx, rz, rp, rqf, uqf)) return 1;
      }

      // Propagate interdependencies
      std::fill_n(w, nrx_+nrz_, 0);
      sp_jac_rdae_.spsolve(w, rx, true);
      std::copy_n(w, nrx_+nrz_, rx);

      // Direct dependency rx_prev -> rx
      std::copy_n(rx, nrx_, rx_prev);

      // Indirect dependency via g
      if (bdae_sp_reverse(&m, x, z, p, u, rx_prev, rp, rx, rz)) return 1;

      // Update rx, rz
      std::copy_n(rx_prev, nrx_, rx);
      std::fill_n(rz, nrz_, 0);

      // Shift time
      if (rx0) rx0 += nrx_;
      if (rp) rp += nrp_;
      if (uqf) uqf += nuq_;
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
      for (casadi_int i = 0; i < nx_; ++i) x[i] |= xf[i];
      std::fill_n(xf, nx_, 0);
    }
    if (zf) {
      for (casadi_int i = 0; i < nz_; ++i) z[i] |= zf[i];
      std::fill_n(zf, nz_, 0);
    }

    // Get dependencies from forward quadratures, if any
    if (nq_ > 0 && qf) {
      if (fquad_sp_reverse(&m, x, z, p, u, qf)) return 1;
    }

    // Propagate interdependencies
    std::fill_n(w, nx_ + nz_, 0);
    sp_jac_dae_.spsolve(w, x, true);
    std::copy_n(w, nx_ + nz_, x);

    // Direct dependency x_prev -> x
    std::copy_n(x, nx_, x_prev);

    // Indirect dependency through f
    if (fdae_sp_reverse(&m, x_prev, p, u, x, z)) return 1;

    // Update x, z
    std::copy_n(x_prev, nx_, x);
    std::fill_n(z, nz_, 0);
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
  Function this_dae, this_rdae;
  this_dae = augmented_dae(&this_rdae);

  // Create integrator for augmented DAE
  std::string aug_prefix = "fsens" + str(nfwd) + "_";
  aug_opts["derivative_of"] = self();
  aug_opts["nfwd"] = nfwd;
  if (!this_rdae.is_null()) aug_opts["rdae"] = this_rdae;
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
      std::vector<MX> ret_in_split = horzsplit(ret_in[i], ret_in[i].size2() / nt());
      // Split augmented input by grid point
      std::vector<std::vector<MX>> aug_in_split(nfwd);
      for (casadi_int d = 0; d < nfwd; ++d) {
        aug_in_split[d] = horzsplit(aug_in[i][d], aug_in[i][d].size2() / nt());
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

  // Create derivative function and return
  return Function(name, ret_in, ret_out, inames, onames, opts);
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

  // Get the current oracle, augmented with forward sensitivity equations if any
  Function this_dae, this_rdae;
  this_dae = augmented_dae(&this_rdae);

  // Create integrator for augmented DAE
  Function aug_dae, aug_rdae;
  std::string aug_prefix = "asens" + str(nadj) + "_";
  if (this_dae.is_a("SXFunction")) {
    aug_dae = aug_adj<SX>(this_dae, this_rdae, nadj, &aug_rdae);
  } else {
    aug_dae = aug_adj<MX>(this_dae, this_rdae, nadj, &aug_rdae);
  }
  aug_opts["derivative_of"] = self();
  aug_opts["nfwd"] = 0;
  if (!aug_rdae.is_null()) aug_opts["rdae"] = aug_rdae;
  Function aug_int = integrator(aug_prefix + name_, plugin_name(),
    aug_dae, t0_, tout_, aug_opts);

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
      ret_in_split = horzsplit(ret_in[i], ret_in[i].size2() / nt());
      // Split augmented input by grid point
      for (casadi_int d = 0; d < nadj; ++d) {
        aug_in_split[d] = horzsplit(aug_in[j][d], aug_in[j][d].size2() / nt());
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
    integrator_in[i] = reshape(vertcat(v), -1, n_grid);
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

  // Create derivative function and return
  return Function(name, ret_in, ret_out, inames, onames, opts);
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
  Sparsity J_xx = F.jac_sparsity(FDAE_ODE, FDYN_X) + Sparsity::diag(nx1_);
  Sparsity J_xz = F.jac_sparsity(FDAE_ODE, FDYN_Z);
  Sparsity J_zx = F.jac_sparsity(FDAE_ALG, FDYN_X);
  Sparsity J_zz = F.jac_sparsity(FDAE_ALG, FDYN_Z);
  // Augment with sensitivity equations
  if (nfwd_ > 0) {
    const Function& fwd_F = get_function(forward_name("daeF", 1));
    J_xx = sp_jac_aug(J_xx, fwd_F.jac_sparsity(FDAE_ODE, FDYN_X));
    J_xz = sp_jac_aug(J_xz, fwd_F.jac_sparsity(FDAE_ODE, FDYN_Z));
    J_zx = sp_jac_aug(J_zx, fwd_F.jac_sparsity(FDAE_ALG, FDYN_X));
    J_zz = sp_jac_aug(J_zz, fwd_F.jac_sparsity(FDAE_ALG, FDYN_Z));
  }
  // Assemble the block matrix
  return blockcat(J_xx, J_xz, J_zx, J_zz);
}

Sparsity Integrator::sp_jac_rdae() {
  // Get the functions
  const Function& G = get_function("daeB");
  // Sparsity pattern for nonaugmented system
  Sparsity J_xx = G.jac_sparsity(BDAE_RODE, BDYN_RX) + Sparsity::diag(nrx1_);
  Sparsity J_xz = G.jac_sparsity(BDAE_RODE, BDYN_RZ);
  Sparsity J_zx = G.jac_sparsity(BDAE_RALG, BDYN_RX);
  Sparsity J_zz = G.jac_sparsity(BDAE_RALG, BDYN_RZ);
  // Augment with sensitivity equations
  if (nfwd_ > 0) {
    const Function& fwd_G = get_function(forward_name("daeB", 1));
    J_xx = sp_jac_aug(J_xx, fwd_G.jac_sparsity(BDAE_RODE, BDYN_RX));
    J_xz = sp_jac_aug(J_xz, fwd_G.jac_sparsity(BDAE_RODE, BDYN_RZ));
    J_zx = sp_jac_aug(J_zx, fwd_G.jac_sparsity(BDAE_RALG, BDYN_RX));
    J_zz = sp_jac_aug(J_zz, fwd_G.jac_sparsity(BDAE_RALG, BDYN_RZ));
  }
  // Assemble the block matrix
  return blockcat(J_xx, J_xz, J_zx, J_zz);
}

std::map<std::string, Integrator::Plugin> Integrator::solvers_;

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
    Function F = get_function("stepF");

    MX z0 = MX::sym("z0", sparsity_in(INTEGRATOR_Z0));

    // Create symbols
    std::vector<MX> F_in = F.mx_in();

    // Prepare return Function inputs
    std::vector<MX> intg_in(INTEGRATOR_NUM_IN);
    intg_in[INTEGRATOR_X0] = F_in[FSTEP_X0];
    intg_in[INTEGRATOR_P] = F_in[FSTEP_P];
    intg_in[INTEGRATOR_U] = F_in[FSTEP_U];
    intg_in[INTEGRATOR_Z0] = z0;
    F_in[FSTEP_V0] = algebraic_state_init(intg_in[INTEGRATOR_X0], z0);

    // Number of finite elements and time steps
    double h = (tout_.back() - t0_)/static_cast<double>(disc_.back());

    // Prepare return Function outputs
    std::vector<MX> intg_out(INTEGRATOR_NUM_OUT);
    F_in[FSTEP_T] = t0_;
    F_in[FSTEP_H] = h;

    std::vector<MX> F_out;
    // Loop over finite elements
    for (casadi_int k=0; k<disc_.back(); ++k) {
      F_out = F(F_in);

      F_in[FSTEP_X0] = F_out[FSTEP_XF];
      F_in[FSTEP_V0] = F_out[FSTEP_VF];
      intg_out[INTEGRATOR_QF] = k==0? F_out[FSTEP_QF] : intg_out[INTEGRATOR_QF]+F_out[FSTEP_QF];
      F_in[FSTEP_T] += h;
    }

    intg_out[INTEGRATOR_XF] = F_out[FSTEP_XF];

    // If-clause needed because rk abuses FSTEP_VF output for intermediate state output
    if (nz_) {
      intg_out[INTEGRATOR_ZF] = algebraic_state_output(F_out[FSTEP_VF]);
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

  // Create dynamic functions, forward and backward problem
  create_function("dynF", fdyn_in(), fdyn_out());
  if (nrx1_ > 0) set_function(rdae_);

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
  const Function& F = get_function(has_function("stepF") ? "stepF" : "implicit_stepF");
  nv1_ = F.nnz_in(FSTEP_V0);
  if (nrx1_ > 0) {
    const Function& G = get_function(has_function("stepB") ? "stepB" : "implicit_stepB");
    nrv1_ = G.nnz_in(BSTEP_RV0);
  } else {
    nrv1_ = 0;
  }
  nv_ = nv1_ * (1 + nfwd_);
  nrv_ = nrv1_ * (1 + nfwd_);

  // Work vectors, forward problem
  alloc_w(nv_, true); // v
  alloc_w(np_, true); // p
  alloc_w(nu_, true); // u
  alloc_w(nq_, true); // q
  alloc_w(nv_, true); // v_prev
  alloc_w(nq_, true); // q_prev

  // Work vectors, backward problem
  alloc_w(nrv_, true); // rv
  alloc_w(nrp_, true); // rp
  alloc_w(nuq_, true); // uq
  alloc_w(nrv_, true); // rv_prev
  alloc_w(nrq_, true); // rq_prev
  alloc_w(nuq_, true); // uq_prev

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

  // Work vectors, allocated in base class
  m->x = w; w += nx_;
  m->z = w; w += nz_;
  m->x_prev = w; w += nx_;
  m->rx = w; w += nrx_;
  m->rz = w; w += nrz_;
  m->rx_prev = w; w += nrx_;
  m->rq = w; w += nrq_;

  // Work vectors, forward problem
  m->v = w; w += nv_;
  m->p = w; w += np_;
  m->u = w; w += nu_;
  m->q = w; w += nq_;
  m->v_prev = w; w += nv_;
  m->q_prev = w; w += nq_;

  // Work vectors, backward problem
  m->rv = w; w += nrv_;
  m->rp = w; w += nrp_;
  m->uq = w; w += nuq_;
  m->rv_prev = w; w += nrv_;
  m->rq_prev = w; w += nrq_;
  m->uq_prev = w; w += nuq_;

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
    double* rx, double* rz, double* rq, double* uq) const {
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
    casadi_copy(m->rx, nrx_, m->rx_prev);
    casadi_copy(m->rv, nrv_, m->rv_prev);
    casadi_copy(m->rq, nrq_, m->rq_prev);
    casadi_copy(m->uq, nuq_, m->uq_prev);

    // Take step
    casadi_int tapeind = disc_[m->k] + j;
    stepB(m, t, h, m->x_tape + nx_ * tapeind, m->v_tape + nv_ * tapeind,
      m->rx_prev, m->rv_prev, m->rx, m->rv, m->rq, m->uq);
    casadi_axpy(nrq_, 1., m->rq_prev, m->rq);
    casadi_axpy(nuq_, 1., m->uq_prev, m->uq);
  }

  // Return to user
  casadi_copy(m->rx, nrx_, rx);
  casadi_copy(m->rv + nrv_ - nrz_, nrz_, rz);
  casadi_copy(m->rq, nrq_, rq);
  casadi_copy(m->uq, nuq_, uq);
}

void FixedStepIntegrator::stepF(FixedStepMemory* m, double t, double h,
    const double* x0, const double* v0, double* xf, double* vf, double* qf) const {
  // Evaluate nondifferentiated
  std::fill(m->arg, m->arg + FSTEP_NUM_IN, nullptr);
  m->arg[FSTEP_T] = &t;  // t
  m->arg[FSTEP_H] = &h;  // h
  m->arg[FSTEP_X0] = x0;  // x0
  m->arg[FSTEP_V0] = v0;  // v0
  m->arg[FSTEP_P] = m->p;  // p
  m->arg[FSTEP_U] = m->u;  // u
  std::fill(m->res, m->res + FSTEP_NUM_OUT, nullptr);
  m->res[FSTEP_XF] = xf;  // xf
  m->res[FSTEP_VF] = vf;  // vf
  m->res[FSTEP_QF] = qf;  // qf
  calc_function(m, "stepF");
  // Evaluate sensitivities
  if (nfwd_ > 0) {
    m->arg[FSTEP_NUM_IN + FSTEP_XF] = xf;  // out:xf
    m->arg[FSTEP_NUM_IN + FSTEP_VF] = vf;  // out:vf
    m->arg[FSTEP_NUM_IN + FSTEP_QF] = qf;  // out:qf
    m->arg[FSTEP_NUM_IN + FSTEP_NUM_OUT + FSTEP_T] = nullptr;  // fwd:t
    m->arg[FSTEP_NUM_IN + FSTEP_NUM_OUT + FSTEP_H] = nullptr;  // fwd:h
    m->arg[FSTEP_NUM_IN + FSTEP_NUM_OUT + FSTEP_X0] = x0 + nx1_;  // fwd:x0
    m->arg[FSTEP_NUM_IN + FSTEP_NUM_OUT + FSTEP_V0] = v0 + nv1_;  // fwd:v0
    m->arg[FSTEP_NUM_IN + FSTEP_NUM_OUT + FSTEP_P] = m->p + np1_;  // fwd:p
    m->arg[FSTEP_NUM_IN + FSTEP_NUM_OUT + FSTEP_U] = m->u + nu1_;  // fwd:u
    m->res[FSTEP_XF] = xf + nx1_;  // fwd:xf
    m->res[FSTEP_VF] = vf + nv1_;  // fwd:vf
    m->res[FSTEP_QF] = qf + nq1_;  // fwd:qf
    calc_function(m, forward_name("stepF", nfwd_));
  }
}

void FixedStepIntegrator::stepB(FixedStepMemory* m, double t, double h,
    const double* x, const double* v, const double* rx0, const double* rv0,
    double* rxf, double* rvf, double* rqf, double* uqf) const {
  // Evaluate nondifferentiated
  std::fill(m->arg, m->arg + BSTEP_NUM_IN, nullptr);
  m->arg[BSTEP_T] = &t;  // t
  m->arg[BSTEP_H] = &h;  // h
  m->arg[BSTEP_RX0] = rx0;  // rx0
  m->arg[BSTEP_RV0] = rv0;  // rv0
  m->arg[BSTEP_RP] = m->rp;  // rp
  m->arg[BSTEP_X] = x;  // x
  m->arg[BSTEP_V] = v;  // v
  m->arg[BSTEP_P] = m->p;  // p
  m->arg[BSTEP_U] = m->u;  // u
  std::fill(m->res, m->res + BSTEP_NUM_OUT, nullptr);
  m->res[BSTEP_RXF] = rxf;  // rxf
  m->res[BSTEP_RVF] = rvf;  // rvf
  m->res[BSTEP_RQF] = rqf;  // rqf
  m->res[BSTEP_UQF] = uqf;  // uqf
  calc_function(m, "stepB");
  // Evaluate sensitivities
  if (nfwd_ > 0) {
    m->arg[BSTEP_NUM_IN + BSTEP_RXF] = rxf;  // out:rxf
    m->arg[BSTEP_NUM_IN + BSTEP_RVF] = rvf;  // out:rvf
    m->arg[BSTEP_NUM_IN + BSTEP_RQF] = rqf;  // out:rqf
    m->arg[BSTEP_NUM_IN + BSTEP_UQF] = uqf;  // out:uqf
    m->arg[BSTEP_NUM_IN + BSTEP_NUM_OUT + BSTEP_T] = nullptr;  // fwd:t
    m->arg[BSTEP_NUM_IN + BSTEP_NUM_OUT + BSTEP_H] = nullptr;  // fwd:h
    m->arg[BSTEP_NUM_IN + BSTEP_NUM_OUT + BSTEP_RX0] = rx0 + nrx1_;  // fwd:rx0
    m->arg[BSTEP_NUM_IN + BSTEP_NUM_OUT + BSTEP_RV0] = rv0 + nrv1_;  // fwd:rv0
    m->arg[BSTEP_NUM_IN + BSTEP_NUM_OUT + BSTEP_RP] = m->rp + nrp1_;  // fwd:rp
    m->arg[BSTEP_NUM_IN + BSTEP_NUM_OUT + BSTEP_X] = x + nx1_;  // fwd:x
    m->arg[BSTEP_NUM_IN + BSTEP_NUM_OUT + BSTEP_V] = v + nv1_;  // fwd:v
    m->arg[BSTEP_NUM_IN + BSTEP_NUM_OUT + BSTEP_P] = m->p + np1_;  // fwd:p
    m->arg[BSTEP_NUM_IN + BSTEP_NUM_OUT + BSTEP_U] = m->u + nu1_;  // fwd:u
    m->res[BSTEP_RXF] = rxf + nrx1_;  // fwd:rxf
    m->res[BSTEP_RVF] = rvf + nrv1_;  // fwd:rvf
    m->res[BSTEP_RQF] = rqf + nrq1_;  // fwd:rqf
    m->res[BSTEP_UQF] = uqf + nuq1_;  // fwd:uqf
    calc_function(m, forward_name("stepB", nfwd_));
  }
}

void FixedStepIntegrator::reset(IntegratorMemory* mem, const double* x, const double* z,
    const double* p) const {
  auto m = static_cast<FixedStepMemory*>(mem);

  // Set parameters
  casadi_copy(p, np_, m->p);

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

void FixedStepIntegrator::resetB(IntegratorMemory* mem,
    const double* rx, const double* rz, const double* rp) const {
  auto m = static_cast<FixedStepMemory*>(mem);

  // Set parameters
  casadi_copy(rp, nrp_, m->rp);

  // Update the state
  casadi_copy(rx, nrx_, m->rx);
  casadi_copy(rz, nrz_, m->rz);

  // Reset summation states
  casadi_clear(m->rq, nrq_);
  casadi_clear(m->uq, nuq_);

  // Get consistent initial conditions
  casadi_fill(m->rv, nrv_, std::numeric_limits<double>::quiet_NaN());
}

void FixedStepIntegrator::impulseB(IntegratorMemory* mem,
    const double* rx, const double* rz, const double* rp) const {
  auto m = static_cast<FixedStepMemory*>(mem);
  // Add impulse to backward parameters
  casadi_axpy(nrp_, 1., rp, m->rp);

  // Add impulse to state
  casadi_axpy(nrx_, 1., rx, m->rx);
  casadi_axpy(nrz_, 1., rz, m->rz);
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
  rootfinder_options["implicit_input"] = FSTEP_V0;
  rootfinder_options["implicit_output"] = FSTEP_VF;

  // Allocate a solver
  Function rf = rootfinder("stepF", implicit_function_name,
    get_function("implicit_stepF"), rootfinder_options);
  set_function(rf);
  if (nfwd_ > 0) set_function(rf.forward(nfwd_));

  // Allocate a root-finding solver for the backward problem
  if (nrv1_ > 0) {
    // Options
    Dict backward_rootfinder_options = rootfinder_options;
    backward_rootfinder_options["implicit_input"] = BSTEP_RV0;
    backward_rootfinder_options["implicit_output"] = BSTEP_RVF;
    std::string backward_implicit_function_name = implicit_function_name;

    // Allocate a Newton solver
    Function brf = rootfinder("stepB", backward_implicit_function_name,
      get_function("implicit_stepB"), backward_rootfinder_options);
    set_function(brf);
    if (nfwd_ > 0) set_function(brf.forward(nfwd_));
  }
}

template<typename XType>
Function Integrator::map2oracle(const std::string& name,
    const std::map<std::string, XType>& d, Function *rdae) {
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
    } else if (i.first=="rx") {
      de_in[DYN_RX]=i.second;
    } else if (i.first=="rz") {
      de_in[DYN_RZ]=i.second;
    } else if (i.first=="rp") {
      de_in[DYN_RP]=i.second;
    } else if (i.first=="ode") {
      de_out[DYN_ODE]=i.second;
    } else if (i.first=="alg") {
      de_out[DYN_ALG]=i.second;
    } else if (i.first=="quad") {
      de_out[DYN_QUAD]=i.second;
    } else if (i.first=="rode") {
      de_out[DYN_RODE]=i.second;
    } else if (i.first=="ralg") {
      de_out[DYN_RALG]=i.second;
    } else if (i.first=="rquad") {
      de_out[DYN_RQUAD]=i.second;
    } else if (i.first=="uquad") {
      de_out[DYN_UQUAD]=i.second;
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

  // Backwards DAE, if any
  if (rdae && de_in[DYN_RX].numel() > 0) {
    // Sort expressions by name
    std::map<std::string, XType> r;
    for (casadi_int i = 0; i < DYN_NUM_IN; ++i) r[dyn_in(i)] = de_in[i];
    for (casadi_int i = 0; i < DYN_NUM_OUT; ++i) r[dyn_out(i)] = de_out[i];
    // Create new expressions corresponding to the rdae IO scheme
    std::vector<XType> rdae_in, rdae_out;
    for (auto& n : bdyn_in()) rdae_in.push_back(r.at(n));
    for (auto& n : bdyn_out()) rdae_out.push_back(r.at(n));
    // Form rdae instance
    *rdae = Function("rdae", rdae_in, rdae_out, bdyn_in(), bdyn_out());
  }

  // Construct
  return Function(name, de_in, de_out, dyn_in(), dyn_out());
}

void Integrator::serialize_body(SerializingStream &s) const {
  OracleFunction::serialize_body(s);

  s.version("Integrator", 2);

  s.pack("Integrator::sp_jac_dae", sp_jac_dae_);
  s.pack("Integrator::sp_jac_rdae", sp_jac_rdae_);
  s.pack("Integrator::t0", t0_);
  s.pack("Integrator::tout", tout_);
  s.pack("Integrator::nfwd", nfwd_);
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
  s.version("Integrator", 2);

  s.unpack("Integrator::sp_jac_dae", sp_jac_dae_);
  s.unpack("Integrator::sp_jac_rdae", sp_jac_rdae_);
  s.unpack("Integrator::t0", t0_);
  s.unpack("Integrator::tout", tout_);
  s.unpack("Integrator::nfwd", nfwd_);
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

  s.unpack("Integrator::augmented_options", augmented_options_);
  s.unpack("Integrator::opts", opts_);
  s.unpack("Integrator::print_stats", print_stats_);
}

void FixedStepIntegrator::serialize_body(SerializingStream &s) const {
  Integrator::serialize_body(s);

  s.version("FixedStepIntegrator", 2);
  s.pack("FixedStepIntegrator::nk_target", nk_target_);
  s.pack("FixedStepIntegrator::disc", disc_);
  s.pack("FixedStepIntegrator::nv", nv_);
  s.pack("FixedStepIntegrator::nv1", nv1_);
  s.pack("FixedStepIntegrator::nrv", nrv_);
  s.pack("FixedStepIntegrator::nrv1", nrv1_);
}

FixedStepIntegrator::FixedStepIntegrator(DeserializingStream & s) : Integrator(s) {
  s.version("FixedStepIntegrator", 2);
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


} // namespace casadi

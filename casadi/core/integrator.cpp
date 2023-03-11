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
  return integrator(name, solver, dae, 0, {1}, opts);
}

Function integrator(const std::string& name, const std::string& solver,
    const MXDict& dae, const Dict& opts) {
  return integrator(name, solver, dae, 0, {1}, opts);
}

Function integrator(const std::string& name, const std::string& solver,
    const Function& dae, const Dict& opts) {
  return integrator(name, solver, dae, 0, {1}, opts);
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
  print_stats_ = false;
}

Integrator::~Integrator() {
}

Sparsity Integrator::get_sparsity_in(casadi_int i) {
  switch (static_cast<IntegratorInput>(i)) {
  case INTEGRATOR_X0: return x();
  case INTEGRATOR_P: return p();
  case INTEGRATOR_U: return repmat(u(), 1, nt());
  case INTEGRATOR_Z0: return z();
  case INTEGRATOR_RX0: return repmat(rx(), 1, nt());
  case INTEGRATOR_RP: return repmat(rp(), 1, nt());
  case INTEGRATOR_RZ0: return repmat(rz(), 1, nt());
  case INTEGRATOR_NUM_IN: break;
  }
  return Sparsity();
}

Sparsity Integrator::get_sparsity_out(casadi_int i) {
  switch (static_cast<IntegratorOutput>(i)) {
  case INTEGRATOR_XF: return repmat(x(), 1, nt());
  case INTEGRATOR_QF: return repmat(q(), 1, nt());
  case INTEGRATOR_ZF: return repmat(z(), 1, nt());
  case INTEGRATOR_RXF: return rx();
  case INTEGRATOR_RQF: return rq();
  case INTEGRATOR_RZF: return rz();
  case INTEGRATOR_UQF: return repmat(uq(), 1, nt());
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

int Integrator::
eval(const double** arg, double** res, casadi_int* iw, double* w, void* mem) const {
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
    if (grid.empty()) {
      grid = {t0, tf};
    }

    // Construct t0 and tout from grid and output_t0
    t0_ = grid.front();
    tout_ = grid;
    if (!output_t0) tout_.erase(tout_.begin());
  }

  // Call the base class method
  OracleFunction::init(opts);

  // For sparsity pattern propagation
  alloc(oracle_);

  // Error if sparse input
  casadi_assert(x().is_dense(), "Sparse DAE not supported");
  casadi_assert(z().is_dense(), "Sparse DAE not supported");
  casadi_assert(p().is_dense(), "Sparse DAE not supported");
  casadi_assert(rx().is_dense(), "Sparse DAE not supported");
  casadi_assert(rz().is_dense(), "Sparse DAE not supported");
  casadi_assert(rp().is_dense(), "Sparse DAE not supported");

  // Get dimensions (excluding sensitivity equations)
  nx1_ = x().size1();
  nz1_ = z().size1();
  nq1_ = q().size1();
  np1_  = p().size1();
  nu1_  = u().size1();
  nrx1_ = rx().size1();
  nrz1_ = rz().size1();
  nrp1_ = rp().size1();
  nrq1_ = rq().size1();
  nuq1_ = uq().size1();

  // Get dimensions (including sensitivity equations)
  nx_ = x().nnz();
  nz_ = z().nnz();
  nq_ = q().nnz();
  np_  = p().nnz();
  nu_  = u().nnz();
  nrx_ = rx().nnz();
  nrz_ = rz().nnz();
  nrp_ = rp().nnz();
  nrq_ = rq().nnz();
  nuq_ = uq().nnz();

  // Number of sensitivities
  ns_ = x().size2()-1;

  // Get the sparsities of the forward and reverse DAE
  sp_jac_dae_ = sp_jac_dae();
  casadi_assert(!sp_jac_dae_.is_singular(),
                        "Jacobian of the forward problem is structurally rank-deficient. "
                        "sprank(J)=" + str(sprank(sp_jac_dae_)) + "<"
                        + str(nx_+nz_));
  if (nrx_>0) {
    sp_jac_rdae_ = sp_jac_rdae();
    casadi_assert(!sp_jac_rdae_.is_singular(),
                          "Jacobian of the backward problem is structurally rank-deficient. "
                          "sprank(J)=" + str(sprank(sp_jac_rdae_)) + "<"
                          + str(nrx_+nrz_));
  }

  // Consistency check

  // Allocate sufficiently large work vectors
  alloc_w(nx_+nz_);
  alloc_w(nrx_+nrz_);
  alloc_w(nx_ + nz_ + nrx_ + nrz_, true);
}

int Integrator::init_mem(void* mem) const {
  if (OracleFunction::init_mem(mem)) return 1;

  //auto m = static_cast<IntegratorMemory*>(mem);
  return 0;
}

template<typename MatType>
std::map<std::string, MatType> Integrator::aug_fwd(casadi_int nfwd) const {
  if (verbose_) casadi_message(name_ + "::aug_fwd");

  // Get input expressions
  std::vector<MatType> arg = MatType::get_input(oracle_);
  std::vector<MatType> aug_x, aug_z, aug_p, aug_u, aug_rx, aug_rz, aug_rp;
  MatType aug_t = arg.at(DYN_T);
  aug_x.push_back(vec(arg.at(DYN_X)));
  aug_z.push_back(vec(arg.at(DYN_Z)));
  aug_p.push_back(vec(arg.at(DYN_P)));
  aug_u.push_back(vec(arg.at(DYN_U)));
  aug_rx.push_back(vec(arg.at(DYN_RX)));
  aug_rz.push_back(vec(arg.at(DYN_RZ)));
  aug_rp.push_back(vec(arg.at(DYN_RP)));

  // Get output expressions
  std::vector<MatType> res = oracle_(arg);
  std::vector<MatType> aug_ode, aug_alg, aug_quad, aug_rode, aug_ralg, aug_rquad, aug_uquad;
  aug_ode.push_back(vec(res.at(DYN_ODE)));
  aug_alg.push_back(vec(res.at(DYN_ALG)));
  aug_quad.push_back(vec(res.at(DYN_QUAD)));
  aug_rode.push_back(vec(res.at(DYN_RODE)));
  aug_ralg.push_back(vec(res.at(DYN_RALG)));
  aug_rquad.push_back(vec(res.at(DYN_RQUAD)));
  aug_uquad.push_back(vec(res.at(DYN_UQUAD)));

  // Zero of time dimension
  MatType zero_t = MatType::zeros(t());

  // Forward directional derivatives
  std::vector<std::vector<MatType>> seed(nfwd, std::vector<MatType>(DYN_NUM_IN));
  for (casadi_int d=0; d<nfwd; ++d) {
    seed[d][DYN_T] = zero_t;
    std::string pref = "aug" + str(d) + "_";
    aug_x.push_back(vec(seed[d][DYN_X] = MatType::sym(pref + "x", x())));
    aug_z.push_back(vec(seed[d][DYN_Z] = MatType::sym(pref + "z", z())));
    aug_p.push_back(vec(seed[d][DYN_P] = MatType::sym(pref + "p", p())));
    aug_u.push_back(vec(seed[d][DYN_U] = MatType::sym(pref + "u", u())));
    aug_rx.push_back(vec(seed[d][DYN_RX] = MatType::sym(pref + "rx", rx())));
    aug_rz.push_back(vec(seed[d][DYN_RZ] = MatType::sym(pref + "rz", rz())));
    aug_rp.push_back(vec(seed[d][DYN_RP] = MatType::sym(pref + "rp", rp())));
  }

  // Calculate directional derivatives
  std::vector<std::vector<MatType>> sens;
  oracle_->call_forward(arg, res, seed, sens, true, false);

  // Collect sensitivity equations
  casadi_assert_dev(sens.size()==nfwd);
  for (casadi_int d=0; d<nfwd; ++d) {
    casadi_assert_dev(sens[d].size()==DYN_NUM_OUT);
    aug_ode.push_back(vec(project(sens[d][DYN_ODE], x())));
    aug_alg.push_back(vec(project(sens[d][DYN_ALG], z())));
    aug_quad.push_back(vec(project(sens[d][DYN_QUAD], q())));
    aug_rode.push_back(vec(project(sens[d][DYN_RODE], rx())));
    aug_ralg.push_back(vec(project(sens[d][DYN_RALG], rz())));
    aug_rquad.push_back(vec(project(sens[d][DYN_RQUAD], rq())));
    aug_uquad.push_back(vec(project(sens[d][DYN_UQUAD], uq())));
  }

  // Construct return object
  std::map<std::string, MatType> ret;
  ret["t"] = aug_t;
  ret["x"] = horzcat(aug_x);
  ret["z"] = horzcat(aug_z);
  ret["p"] = horzcat(aug_p);
  ret["u"] = horzcat(aug_u);
  ret["ode"] = horzcat(aug_ode);
  ret["alg"] = horzcat(aug_alg);
  ret["quad"] = horzcat(aug_quad);
  ret["rx"] = horzcat(aug_rx);
  ret["rz"] = horzcat(aug_rz);
  ret["rp"] = horzcat(aug_rp);
  ret["rode"] = horzcat(aug_rode);
  ret["ralg"] = horzcat(aug_ralg);
  ret["rquad"] = horzcat(aug_rquad);
  ret["uquad"] = horzcat(aug_uquad);

  return ret;
}

template<typename MatType>
std::map<std::string, MatType> Integrator::aug_adj(casadi_int nadj) const {
  if (verbose_) casadi_message(name_ + "::aug_adj");

  // Get input expressions
  std::vector<MatType> arg = MatType::get_input(oracle_);
  std::vector<MatType> aug_x, aug_z, aug_p, aug_u, aug_rx, aug_rz, aug_rp;
  MatType aug_t = arg.at(DYN_T);
  aug_x.push_back(vec(arg.at(DYN_X)));
  aug_z.push_back(vec(arg.at(DYN_Z)));
  aug_p.push_back(vec(arg.at(DYN_P)));
  aug_u.push_back(vec(arg.at(DYN_U)));
  aug_rx.push_back(vec(arg.at(DYN_RX)));
  aug_rz.push_back(vec(arg.at(DYN_RZ)));
  aug_rp.push_back(vec(arg.at(DYN_RP)));

  // Get output expressions
  std::vector<MatType> res = oracle_(arg);
  std::vector<MatType> aug_ode, aug_alg, aug_quad, aug_rode, aug_ralg, aug_rquad, aug_uquad;
  aug_ode.push_back(vec(res.at(DYN_ODE)));
  aug_alg.push_back(vec(res.at(DYN_ALG)));
  aug_quad.push_back(vec(res.at(DYN_QUAD)));
  aug_rode.push_back(vec(res.at(DYN_RODE)));
  aug_ralg.push_back(vec(res.at(DYN_RALG)));
  aug_rquad.push_back(vec(res.at(DYN_RQUAD)));
  aug_uquad.push_back(vec(res.at(DYN_UQUAD)));

  // Zero of time dimension
  MatType zero_t = MatType::zeros(t());

  // Reverse mode directional derivatives
  std::vector<std::vector<MatType>> seed(nadj, std::vector<MatType>(DYN_NUM_OUT));
  for (casadi_int d=0; d<nadj; ++d) {
    std::string pref = "aug" + str(d) + "_";
    aug_rx.push_back(vec(seed[d][DYN_ODE] = MatType::sym(pref + "ode", x())));
    aug_rz.push_back(vec(seed[d][DYN_ALG] = MatType::sym(pref + "alg", z())));
    aug_rp.push_back(vec(seed[d][DYN_QUAD] = MatType::sym(pref + "quad", q())));
    aug_x.push_back(vec(seed[d][DYN_RODE] = MatType::sym(pref + "rode", rx())));
    aug_z.push_back(vec(seed[d][DYN_RALG] = MatType::sym(pref + "ralg", rz())));
    aug_p.push_back(vec(seed[d][DYN_RQUAD] = MatType::sym(pref + "rquad", rq())));
    aug_u.push_back(vec(seed[d][DYN_UQUAD] = MatType::sym(pref + "uquad", uq())));
  }

  // Calculate directional derivatives
  std::vector<std::vector<MatType>> sens;
  oracle_->call_reverse(arg, res, seed, sens, true, false);

  // Collect sensitivity equations
  casadi_assert_dev(sens.size()==nadj);
  for (casadi_int d=0; d<nadj; ++d) {
    casadi_assert_dev(sens[d].size()==DYN_NUM_IN);
    aug_rode.push_back(vec(project(sens[d][DYN_X], x())));
    aug_ralg.push_back(vec(project(sens[d][DYN_Z], z())));
    aug_rquad.push_back(vec(project(sens[d][DYN_P], p())));
    aug_uquad.push_back(vec(project(sens[d][DYN_U], u())));
    aug_ode.push_back(vec(project(sens[d][DYN_RX], rx())));
    aug_alg.push_back(vec(project(sens[d][DYN_RZ], rz())));
    aug_quad.push_back(vec(project(sens[d][DYN_RP], rp())));
  }

  // Construct return object
  std::map<std::string, MatType> ret;
  ret["t"] = aug_t;
  ret["x"] = vertcat(aug_x);
  ret["z"] = vertcat(aug_z);
  ret["p"] = vertcat(aug_p);
  ret["u"] = vertcat(aug_u);
  ret["ode"] = vertcat(aug_ode);
  ret["alg"] = vertcat(aug_alg);
  ret["quad"] = vertcat(aug_quad);
  ret["rx"] = vertcat(aug_rx);
  ret["rz"] = vertcat(aug_rz);
  ret["rp"] = vertcat(aug_rp);
  ret["rode"] = vertcat(aug_rode);
  ret["ralg"] = vertcat(aug_ralg);
  ret["rquad"] = vertcat(aug_rquad);
  ret["uquad"] =  vertcat(aug_uquad);

  // Make sure that forward problem does not depend on backward states
  Function f("f", {ret["t"], ret["x"], ret["z"], ret["p"], ret["u"]},
                  {ret["ode"], ret["alg"], ret["quad"]}, {{"allow_free", true}});
  if (f.has_free()) {
    // Replace dependencies of rx, rz and rp with zeros
    f = Function("f", {ret["t"], ret["x"], ret["z"], ret["p"], ret["u"],
                        ret["rx"], ret["rz"], ret["rp"]},
                      {ret["ode"], ret["alg"], ret["quad"]});
    std::vector<MatType> v = {ret["t"], ret["x"], ret["z"], ret["p"], ret["u"], 0, 0, 0};
    v = f(v);
    ret["ode"] = v.at(0);
    ret["alg"] = v.at(1);
    ret["quad"] = v.at(2);
  }

  return ret;
}

int Integrator::
sp_forward(const bvec_t** arg, bvec_t** res, casadi_int* iw, bvec_t* w, void* mem) const {
  if (verbose_) casadi_message(name_ + "::sp_forward");

  // Work vectors
  bvec_t *tmp_x = w; w += nx_;
  bvec_t *tmp_z = w; w += nz_;
  bvec_t *tmp_rx = w; w += nrx_;
  bvec_t *tmp_rz = w; w += nrz_;

  // Propagate forward
  const bvec_t** arg1 = arg+n_in_;
  std::fill_n(arg1, static_cast<size_t>(DYN_NUM_IN), nullptr);
  arg1[DYN_X] = arg[INTEGRATOR_X0];
  arg1[DYN_P] = arg[INTEGRATOR_P];
  arg1[DYN_U] = arg[INTEGRATOR_U];
  bvec_t** res1 = res+n_out_;
  std::fill_n(res1, static_cast<size_t>(DYN_NUM_OUT), nullptr);
  res1[DYN_ODE] = tmp_x;
  res1[DYN_ALG] = tmp_z;
  oracle_(arg1, res1, iw, w, 0);
  if (arg[INTEGRATOR_X0]) {
    const bvec_t *tmp = arg[INTEGRATOR_X0];
    for (casadi_int i=0; i<nx_; ++i) tmp_x[i] |= *tmp++;
  }

  // "Solve" in order to resolve interdependencies (cf. Rootfinder)
  std::copy_n(tmp_x, nx_+nz_, w);
  std::fill_n(tmp_x, nx_+nz_, 0);
  sp_jac_dae_.spsolve(tmp_x, w, false);

  // Get xf and zf
  if (res[INTEGRATOR_XF]) std::copy_n(tmp_x, nx_, res[INTEGRATOR_XF]);
  if (res[INTEGRATOR_ZF]) std::copy_n(tmp_z, nz_, res[INTEGRATOR_ZF]);

  // Propagate to quadratures
  if (nq_>0 && res[INTEGRATOR_QF]) {
    arg1[DYN_X] = tmp_x;
    arg1[DYN_Z] = tmp_z;
    res1[DYN_ODE] = res1[DYN_ALG] = nullptr;
    res1[DYN_QUAD] = res[INTEGRATOR_QF];
    if (oracle_(arg1, res1, iw, w, 0)) return 1;
  }

  if (nrx_>0) {
    // Propagate through g
    std::fill_n(arg1, static_cast<size_t>(DYN_NUM_IN), nullptr);
    arg1[DYN_X] = tmp_x;
    arg1[DYN_P] = arg[INTEGRATOR_P];
    arg1[DYN_U] = arg[INTEGRATOR_U];
    arg1[DYN_Z] = tmp_z;
    arg1[DYN_RX] = arg[INTEGRATOR_X0];
    arg1[DYN_RX] = arg[INTEGRATOR_RX0];
    arg1[DYN_RP] = arg[INTEGRATOR_RP];
    std::fill_n(res1, static_cast<size_t>(DYN_NUM_OUT), nullptr);
    res1[DYN_RODE] = tmp_rx;
    res1[DYN_RALG] = tmp_rz;
    oracle_(arg1, res1, iw, w, 0);
    if (arg[INTEGRATOR_RX0]) {
      const bvec_t *tmp = arg[INTEGRATOR_RX0];
      for (casadi_int i=0; i<nrx_; ++i) tmp_rx[i] |= *tmp++;
    }

    // "Solve" in order to resolve interdependencies (cf. Rootfinder)
    std::copy_n(tmp_rx, nrx_+nrz_, w);
    std::fill_n(tmp_rx, nrx_+nrz_, 0);
    sp_jac_rdae_.spsolve(tmp_rx, w, false);

    // Get rxf and rzf
    if (res[INTEGRATOR_RXF]) std::copy_n(tmp_rx, nrx_, res[INTEGRATOR_RXF]);
    if (res[INTEGRATOR_RZF]) std::copy_n(tmp_rz, nrz_, res[INTEGRATOR_RZF]);

    // Propagate to quadratures
    if (nrq_>0 && res[INTEGRATOR_RQF]) {
      arg1[DYN_RX] = tmp_rx;
      arg1[DYN_RZ] = tmp_rz;
      res1[DYN_RODE] = res1[DYN_RALG] = nullptr;
      res1[DYN_RQUAD] = res[INTEGRATOR_RQF];
      if (oracle_(arg1, res1, iw, w, 0)) return 1;
    }
  }
  return 0;
}

int Integrator::sp_reverse(bvec_t** arg, bvec_t** res,
    casadi_int* iw, bvec_t* w, void* mem) const {
  if (verbose_) casadi_message(name_ + "::sp_reverse");

  // Work vectors
  bvec_t** arg1 = arg+n_in_;
  bvec_t** res1 = res+n_out_;
  bvec_t *tmp_x = w; w += nx_;
  bvec_t *tmp_z = w; w += nz_;

  // Shorthands
  bvec_t* x0 = arg[INTEGRATOR_X0];
  bvec_t* p = arg[INTEGRATOR_P];
  bvec_t* u = arg[INTEGRATOR_U];
  bvec_t* xf = res[INTEGRATOR_XF];
  bvec_t* zf = res[INTEGRATOR_ZF];
  bvec_t* qf = res[INTEGRATOR_QF];

  // Propagate from outputs to state vectors
  if (xf) {
    std::copy_n(xf, nx_, tmp_x);
    std::fill_n(xf, nx_, 0);
  } else {
    std::fill_n(tmp_x, nx_, 0);
  }
  if (zf) {
    std::copy_n(zf, nz_, tmp_z);
    std::fill_n(zf, nz_, 0);
  } else {
    std::fill_n(tmp_z, nz_, 0);
  }

  if (nrx_>0) {
    // Work vectors
    bvec_t *tmp_rx = w; w += nrx_;
    bvec_t *tmp_rz = w; w += nrz_;

    // Shorthands
    bvec_t* rx0 = arg[INTEGRATOR_RX0];
    bvec_t* rp = arg[INTEGRATOR_RP];
    bvec_t* rxf = res[INTEGRATOR_RXF];
    bvec_t* rzf = res[INTEGRATOR_RZF];
    bvec_t* rqf = res[INTEGRATOR_RQF];
    bvec_t* uqf = res[INTEGRATOR_UQF];

    // Propagate from outputs to state vectors
    if (rxf) {
      std::copy_n(rxf, nrx_, tmp_rx);
      std::fill_n(rxf, nrx_, 0);
    } else {
      std::fill_n(tmp_rx, nrx_, 0);
    }
    if (rzf) {
      std::copy_n(rzf, nrz_, tmp_rz);
      std::fill_n(rzf, nrz_, 0);
    } else {
      std::fill_n(tmp_rz, nrz_, 0);
    }

    // Get dependencies from backward quadratures
    std::fill_n(res1, static_cast<size_t>(DYN_NUM_OUT), nullptr);
    std::fill_n(arg1, static_cast<size_t>(DYN_NUM_IN), nullptr);
    res1[DYN_RQUAD] = rqf;
    res1[DYN_UQUAD] = uqf;
    arg1[DYN_X] = tmp_x;
    arg1[DYN_Z] = tmp_z;
    arg1[DYN_P] = p;
    arg1[DYN_U] = u;
    arg1[DYN_RX] = tmp_rx;
    arg1[DYN_RZ] = tmp_rz;
    arg1[DYN_RP] = rp;
    if (oracle_.rev(arg1, res1, iw, w, 0)) return 1;

    // Propagate interdependencies
    std::fill_n(w, nrx_+nrz_, 0);
    sp_jac_rdae_.spsolve(w, tmp_rx, true);
    std::copy_n(w, nrx_+nrz_, tmp_rx);

    // Direct dependency rx0 -> rxf
    if (rx0) for (casadi_int i=0; i<nrx_; ++i) rx0[i] |= tmp_rx[i];

    // Indirect dependency via g
    res1[DYN_RODE] = tmp_rx;
    res1[DYN_RALG] = tmp_rz;
    res1[DYN_RQUAD] = nullptr;
    res1[DYN_UQUAD] = nullptr;
    arg1[DYN_RX] = rx0;
    arg1[DYN_RZ] = nullptr; // arg[INTEGRATOR_RZ0] is a guess, no dependency
    if (oracle_.rev(arg1, res1, iw, w, 0)) return 1;
  }

  // Get dependencies from forward quadratures
  std::fill_n(res1, static_cast<size_t>(DYN_NUM_OUT), nullptr);
  std::fill_n(arg1, static_cast<size_t>(DYN_NUM_IN), nullptr);
  res1[DYN_QUAD] = qf;
  arg1[DYN_X] = tmp_x;
  arg1[DYN_Z] = tmp_z;
  arg1[DYN_P] = p;
  arg1[DYN_U] = u;
  if (qf && nq_>0) {
    if (oracle_.rev(arg1, res1, iw, w, 0)) return 1;
  }

  // Propagate interdependencies
  std::fill_n(w, nx_+nz_, 0);
  sp_jac_dae_.spsolve(w, tmp_x, true);
  std::copy_n(w, nx_+nz_, tmp_x);

  // Direct dependency x0 -> xf
  if (x0) for (casadi_int i=0; i<nx_; ++i) x0[i] |= tmp_x[i];

  // Indirect dependency through f
  res1[DYN_ODE] = tmp_x;
  res1[DYN_ALG] = tmp_z;
  res1[DYN_QUAD] = nullptr;
  arg1[DYN_X] = x0;
  arg1[DYN_Z] = nullptr; // arg[INTEGRATOR_Z0] is a guess, no dependency
  if (oracle_.rev(arg1, res1, iw, w, 0)) return 1;
  return 0;
}

Function Integrator::
get_forward(casadi_int nfwd, const std::string& name,
            const std::vector<std::string>& inames,
            const std::vector<std::string>& onames,
            const Dict& opts) const {
  if (verbose_) casadi_message(name_ + "::get_forward");

  // Integrator options
  Dict aug_opts = getDerivativeOptions(true);
  for (auto&& i : augmented_options_) {
    aug_opts[i.first] = i.second;
  }

  // Create integrator for augmented DAE
  Function aug_dae;
  std::string aug_prefix = "fsens" + str(nfwd) + "_";
  std::string dae_name = aug_prefix + oracle_.name();
  Dict dae_opts = {{"derivative_of", oracle_}};
  if (oracle_.is_a("SXFunction")) {
    aug_dae = map2oracle(dae_name, aug_fwd<SX>(nfwd));
  } else {
    aug_dae = map2oracle(dae_name, aug_fwd<MX>(nfwd));
  }
  aug_opts["derivative_of"] = self();
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
      integrator_in[i] = horzcat(v);
    } else {
      // No reordering necessary
      v = aug_in[i];
      v.insert(v.begin(), ret_in[i]);
      integrator_in[i] = horzcat(v);
    }
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
    std::vector<MX> integrator_out_split = horzsplit(integrator_out[i], offset);
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

Function Integrator::
get_reverse(casadi_int nadj, const std::string& name,
            const std::vector<std::string>& inames,
            const std::vector<std::string>& onames,
            const Dict& opts) const {
  if (verbose_) casadi_message(name_ + "::get_reverse");

  // Integrator options
  Dict aug_opts = getDerivativeOptions(false);
  for (auto&& i : augmented_options_) {
    aug_opts[i.first] = i.second;
  }

  // Create integrator for augmented DAE
  Function aug_dae;
  std::string aug_prefix = "asens" + str(nadj) + "_";
  std::string dae_name = aug_prefix + oracle_.name();
  if (oracle_.is_a("SXFunction")) {
    aug_dae = map2oracle(dae_name, aug_adj<SX>(nadj));
  } else {
    aug_dae = map2oracle(dae_name, aug_adj<MX>(nadj));
  }
  aug_opts["derivative_of"] = self();
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
      offset.push_back(offset.back() + size1_out(j));
      for (casadi_int d = 0; d < nadj; ++d) {
        offset.push_back(offset.back() + size1_in(i));
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

Sparsity Integrator::sp_jac_dae() {
  // Start with the sparsity pattern of the ODE part
  Sparsity jac_ode_x = oracle_.jac_sparsity(DYN_ODE, DYN_X);

  // Add diagonal to get interdependencies
  jac_ode_x = jac_ode_x + Sparsity::diag(nx_);

  // Quick return if no algebraic variables
  if (nz_==0) return jac_ode_x;

  // Add contribution from algebraic variables and equations
  Sparsity jac_ode_z = oracle_.jac_sparsity(DYN_ODE, DYN_Z);
  Sparsity jac_alg_x = oracle_.jac_sparsity(DYN_ALG, DYN_X);
  Sparsity jac_alg_z = oracle_.jac_sparsity(DYN_ALG, DYN_Z);
  return blockcat(jac_ode_x, jac_ode_z,
                  jac_alg_x, jac_alg_z);
}

Sparsity Integrator::sp_jac_rdae() {
  // Start with the sparsity pattern of the ODE part
  Sparsity jac_ode_x = oracle_.jac_sparsity(DYN_RODE, DYN_RX);

  // Add diagonal to get interdependencies
  jac_ode_x = jac_ode_x + Sparsity::diag(nrx_);

  // Quick return if no algebraic variables
  if (nrz_==0) return jac_ode_x;

  // Add contribution from algebraic variables and equations
  Sparsity jac_ode_z = oracle_.jac_sparsity(DYN_RODE, DYN_RZ);
  Sparsity jac_alg_x = oracle_.jac_sparsity(DYN_RALG, DYN_RX);
  Sparsity jac_alg_z = oracle_.jac_sparsity(DYN_RALG, DYN_RZ);
  return blockcat(jac_ode_x, jac_ode_z,
                  jac_alg_x, jac_alg_z);
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
    Function F = getExplicit();

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
    F_in[FSTEP_T0] = t0_;
    F_in[FSTEP_H] = h;

    std::vector<MX> F_out;
    // Loop over finite elements
    for (casadi_int k=0; k<disc_.back(); ++k) {
      F_out = F(F_in);

      F_in[FSTEP_X0] = F_out[FSTEP_XF];
      F_in[FSTEP_V0] = F_out[FSTEP_VF];
      intg_out[INTEGRATOR_QF] = k==0? F_out[FSTEP_QF] : intg_out[INTEGRATOR_QF]+F_out[FSTEP_QF];
      F_in[FSTEP_T0] += h;
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
  setupFG();

  // Get discrete time dimensions
  nv_ = F_.nnz_in(FSTEP_V0);
  nrv_ =  G_.is_null() ? 0 : G_.nnz_in(BSTEP_RV0);

  // Work vectors, forward problem
  alloc_w(nx_, true); // x
  alloc_w(nz_, true); // z
  alloc_w(np_, true); // p
  alloc_w(nu_, true); // u
  alloc_w(nq_, true); // q

  // Work vectors, backward problem
  alloc_w(nrx_, true); // rx
  alloc_w(nrz_, true); // rz
  alloc_w(nrp_, true); // rp
  alloc_w(nrq_, true); // rq
  alloc_w(nuq_, true); // uq

  /// Work vectors, dependent variables
  alloc_w(nv_, true); // v
  alloc_w(nrv_, true); // rv
}

void FixedStepIntegrator::set_work(void* mem, const double**& arg, double**& res,
    casadi_int*& iw, double*& w) const {
  auto m = static_cast<FixedStepMemory*>(mem);

  // Set work in base classes
  Integrator::set_work(mem, arg, res, iw, w);

  // Work vectors, forward problem
  m->x = w; w += nx_;
  m->z = w; w += nz_;
  m->p = w; w += np_;
  m->u = w; w += nu_;
  m->q = w; w += nq_;

  // Work vectors, backward problem
  m->rx = w; w += nrx_;
  m->rz = w; w += nrz_;
  m->rp = w; w += nrp_;
  m->rq = w; w += nrq_;
  m->uq = w; w += nuq_;

  /// Work vectors, dependent variables
  m->v = w; w += nv_;
  m->rv = w; w += nrv_;
}

int FixedStepIntegrator::init_mem(void* mem) const {
  if (Integrator::init_mem(mem)) return 1;
  auto m = static_cast<FixedStepMemory*>(mem);

  // Allocate tape if backward states are present
  if (nrx_>0) {
    m->x_tape.resize((disc_.back() + 1) * nx_);
    m->v_tape.resize(disc_.back() * nv_);
  }

  // Allocate state
  m->x_prev.resize(nx_);
  m->v_prev.resize(nv_);
  m->q_prev.resize(nq_);
  m->rx_prev.resize(nrx_);
  m->rv_prev.resize(nrv_);
  m->rq_prev.resize(nrq_);
  m->uq_prev.resize(nuq_);

  return 0;
}


void FixedStepIntegrator::advance(IntegratorMemory* mem,
    const double* u, double* x, double* z, double* q) const {
  auto m = static_cast<FixedStepMemory*>(mem);

  // Set controls
  casadi_copy(u, nu_, m->u);

  // Explicit discrete time dynamics
  const Function& F = getExplicit();

  // Number of finite elements and time steps
  casadi_int nj = disc_[m->k + 1] - disc_[m->k];
  double h = (m->t_next - m->t) / nj;

  // Current time
  double t;

  // Discrete dynamics function inputs ...
  std::fill_n(m->arg, F.n_in(), nullptr);
  m->arg[FSTEP_T0] = &t;
  m->arg[FSTEP_H] = &h;
  m->arg[FSTEP_X0] = get_ptr(m->x_prev);
  m->arg[FSTEP_V0] = get_ptr(m->v_prev);
  m->arg[FSTEP_P] = m->p;
  m->arg[FSTEP_U] = m->u;

  // ... and outputs
  std::fill_n(m->res, F.n_out(), nullptr);
  m->res[FSTEP_XF] = m->x;
  m->res[FSTEP_VF] = m->v;
  m->res[FSTEP_QF] = m->q;

  // Take steps
  for (casadi_int j = 0; j < nj; ++j) {
    // Update time
    t = m->t + j * h;

    // Update the previous step
    casadi_copy(m->x, nx_, get_ptr(m->x_prev));
    casadi_copy(m->v, nv_, get_ptr(m->v_prev));
    casadi_copy(m->q, nq_, get_ptr(m->q_prev));

    // Take step
    F(m->arg, m->res, m->iw, m->w);
    casadi_axpy(nq_, 1., get_ptr(m->q_prev), m->q);

    // Save state, if needed
    if (nrx_ > 0) {
      casadi_int tapeind = disc_[m->k] + j;
      casadi_copy(m->x, nx_, get_ptr(m->x_tape) + nx_ * (tapeind + 1));
      casadi_copy(m->v, nv_, get_ptr(m->v_tape) + nv_ * tapeind);
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

  // Explicit discrete time dynamics
  const Function& G = getExplicitB();

  // Number of finite elements and time steps
  casadi_int nj = disc_[m->k + 1] - disc_[m->k];
  double h = (m->t - m->t_next) / nj;

  // Current time
  double t;

  // Discrete dynamics function inputs ...
  std::fill_n(m->arg, G.n_in(), nullptr);
  m->arg[BSTEP_T0] = &t;
  m->arg[BSTEP_H] = &h;
  m->arg[BSTEP_P] = m->p;
  m->arg[BSTEP_U] = m->u;
  m->arg[BSTEP_RX0] = get_ptr(m->rx_prev);
  m->arg[BSTEP_RV0] = get_ptr(m->rv_prev);
  m->arg[BSTEP_RP] = m->rp;

  // ... and outputs
  std::fill_n(m->res, G.n_out(), nullptr);
  m->res[BSTEP_RXF] = m->rx;
  m->res[BSTEP_RVF] = m->rv;
  m->res[BSTEP_RQF] = m->rq;
  m->res[BSTEP_UQF] = m->uq;

  // Take steps
  for (casadi_int j = nj; j-- > 0; ) {
    // Update time
    t = m->t_next + j * h;

    // Update the previous step
    casadi_copy(m->rx, nrx_, get_ptr(m->rx_prev));
    casadi_copy(m->rv, nrv_, get_ptr(m->rv_prev));
    casadi_copy(m->rq, nrq_, get_ptr(m->rq_prev));
    casadi_copy(m->uq, nuq_, get_ptr(m->uq_prev));

    // Take step
    casadi_int tapeind = disc_[m->k] + j;
    m->arg[BSTEP_X] = get_ptr(m->x_tape) + nx_ * tapeind;
    m->arg[BSTEP_V] = get_ptr(m->v_tape) + nv_ * tapeind;
    G(m->arg, m->res, m->iw, m->w);
    casadi_axpy(nrq_, 1., get_ptr(m->rq_prev), m->rq);
    casadi_axpy(nuq_, 1., get_ptr(m->uq_prev), m->uq);
  }

  // Return to user
  casadi_copy(m->rx, nrx_, rx);
  casadi_copy(m->rv + nrv_ - nrz_, nrz_, rz);
  casadi_copy(m->rq, nrq_, rq);
  casadi_copy(m->uq, nuq_, uq);
}

void FixedStepIntegrator::
reset(IntegratorMemory* mem,
      const double* x, const double* z, const double* p) const {
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
    casadi_copy(x, nx_, get_ptr(m->x_tape));
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
  rootfinder_ = rootfinder(name_ + "_rootfinder", implicit_function_name,
    F_, rootfinder_options);
  alloc(rootfinder_);

  // Allocate a root-finding solver for the backward problem
  if (nrv_>0) {
    // Options
    Dict backward_rootfinder_options = rootfinder_options;
    backward_rootfinder_options["implicit_input"] = BSTEP_RV0;
    backward_rootfinder_options["implicit_output"] = BSTEP_RVF;
    std::string backward_implicit_function_name = implicit_function_name;

    // Allocate a Newton solver
    backward_rootfinder_ =
      rootfinder(name_+ "_backward_rootfinder",
                  backward_implicit_function_name,
                  G_, backward_rootfinder_options);
    alloc(backward_rootfinder_);
  }
}

template<typename XType>
Function Integrator::map2oracle(const std::string& name,
  const std::map<std::string, XType>& d, const Dict& opts) {
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

  // Make sure x and ode exist
  casadi_assert(!de_in[DYN_X].is_empty(), "Ill-posed ODE - no state");

  // Number of right-hand-sides
  casadi_int nrhs = de_in[DYN_X].size2();

  // Make sure consistent number of right-hand-sides
  for (bool b : {true, false}) {
    for (auto&& e : b ? de_in : de_out) {
      // Skip time
      if (&e == &de_in[DYN_T]) continue;
      // Number of rows
      casadi_int nr = e.size1();
      // Make sure no change in number of elements
      casadi_assert(e.numel()==nr*nrhs, "Inconsistent number of rhs");
      e = reshape(e, nr, nrhs);
    }
  }

  // Consistent sparsity for x
  casadi_assert(de_in[DYN_X].size()==de_out[DYN_ODE].size(),
    "Dimension mismatch for 'ode'");
  de_out[DYN_ODE] = project(de_out[DYN_ODE], de_in[DYN_X].sparsity());

  // Consistent sparsity for z
  casadi_assert(de_in[DYN_Z].size()==de_out[DYN_ALG].size(),
    "Dimension mismatch for 'alg'");
  de_out[DYN_ALG] = project(de_out[DYN_ALG], de_in[DYN_Z].sparsity());

  // Consistent sparsity for rx
  casadi_assert(de_in[DYN_RX].size()==de_out[DYN_RODE].size(),
    "Dimension mismatch for 'rode'");
  de_out[DYN_RODE] = project(de_out[DYN_RODE], de_in[DYN_RX].sparsity());

  // Consistent sparsity for rz
  casadi_assert(de_in[DYN_RZ].size()==de_out[DYN_RALG].size(),
    "Dimension mismatch for 'ralg'");
  de_out[DYN_RALG] = project(de_out[DYN_RALG], de_in[DYN_RZ].sparsity());

  // Construct
  return Function(name, de_in, de_out, dyn_in(), dyn_out(), opts);
}

void Integrator::serialize_body(SerializingStream &s) const {
  OracleFunction::serialize_body(s);

  s.version("Integrator", 2);
  s.pack("Integrator::sp_jac_dae", sp_jac_dae_);
  s.pack("Integrator::sp_jac_rdae", sp_jac_rdae_);
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
  s.pack("Integrator::nu", nu_);
  s.pack("Integrator::nrp", nrp_);
  s.pack("Integrator::np1", np1_);
  s.pack("Integrator::nu1", nu1_);
  s.pack("Integrator::nrp1", nrp1_);
  s.pack("Integrator::ns", ns_);
  s.pack("Integrator::augmented_options", augmented_options_);
  s.pack("Integrator::opts", opts_);
  s.pack("Integrator::print_stats", print_stats_);
  s.pack("Integrator::t0", t0_);
  s.pack("Integrator::tout", tout_);
}

void Integrator::serialize_type(SerializingStream &s) const {
  OracleFunction::serialize_type(s);
  PluginInterface<Integrator>::serialize_type(s);
}

ProtoFunction* Integrator::deserialize(DeserializingStream& s) {
  return PluginInterface<Integrator>::deserialize(s);
}

Integrator::Integrator(DeserializingStream & s) : OracleFunction(s) {
  int version = s.version("Integrator", 1, 2);
  s.unpack("Integrator::sp_jac_dae", sp_jac_dae_);
  s.unpack("Integrator::sp_jac_rdae", sp_jac_rdae_);
  s.unpack("Integrator::nx", nx_);
  s.unpack("Integrator::nz", nz_);
  s.unpack("Integrator::nq", nq_);
  s.unpack("Integrator::nx1", nx1_);
  s.unpack("Integrator::nz1", nz1_);
  s.unpack("Integrator::nq1", nq1_);
  s.unpack("Integrator::nrx", nrx_);
  s.unpack("Integrator::nrz", nrz_);
  s.unpack("Integrator::nrq", nrq_);
  s.unpack("Integrator::nrx1", nrx1_);
  s.unpack("Integrator::nrz1", nrz1_);
  s.unpack("Integrator::nrq1", nrq1_);
  s.unpack("Integrator::np", np_);
  s.unpack("Integrator::nrp", nrp_);
  s.unpack("Integrator::np1", np1_);
  s.unpack("Integrator::nrp1", nrp1_);
  s.unpack("Integrator::ns", ns_);
  s.unpack("Integrator::augmented_options", augmented_options_);
  s.unpack("Integrator::opts", opts_);
  s.unpack("Integrator::print_stats", print_stats_);
  if (version >= 2) {
    s.unpack("Integrator::t0", t0_);
    s.unpack("Integrator::tout", tout_);
    s.unpack("Integrator::nu", nu_);
    s.unpack("Integrator::nu1", nu1_);
    s.unpack("Integrator::nuq", nuq_);
    s.unpack("Integrator::nuq1", nuq1_);
  } else {
    // Time grid
    std::vector<double> grid;
    s.unpack("Integrator::grid", grid);
    // Is the first time point in output?
    bool output_t0;
    s.unpack("Integrator::output_t0", output_t0);
    // Construct t0 and tout from grid and output_t0
    t0_ = grid.front();
    tout_ = grid;
    if (!output_t0) tout_.erase(tout_.begin());
    // Controls did not exist in version 1
    nu_ = nu1_ = nuq_ = nuq1_ = 0;
  }
}

void FixedStepIntegrator::serialize_body(SerializingStream &s) const {
  Integrator::serialize_body(s);

  s.version("FixedStepIntegrator", 1);
  s.pack("FixedStepIntegrator::F", F_);
  s.pack("FixedStepIntegrator::G", G_);
  //s.pack("FixedStepIntegrator::nk", nk_);
  //s.pack("FixedStepIntegrator::h", h_);
  s.pack("FixedStepIntegrator::nZ", nv_);
  s.pack("FixedStepIntegrator::nRZ", nrv_);
}

FixedStepIntegrator::FixedStepIntegrator(DeserializingStream & s) : Integrator(s) {
  s.version("FixedStepIntegrator", 1);
  s.unpack("FixedStepIntegrator::F", F_);
  s.unpack("FixedStepIntegrator::G", G_);
  //s.unpack("FixedStepIntegrator::nk", nk_);
  //s.unpack("FixedStepIntegrator::h", h_);
  s.unpack("FixedStepIntegrator::nZ", nv_);
  s.unpack("FixedStepIntegrator::nRZ", nrv_);
}

void ImplicitFixedStepIntegrator::serialize_body(SerializingStream &s) const {
  FixedStepIntegrator::serialize_body(s);

  s.version("ImplicitFixedStepIntegrator", 1);
  s.pack("ImplicitFixedStepIntegrator::rootfinder", rootfinder_);
  s.pack("ImplicitFixedStepIntegrator::backward_rootfinder", backward_rootfinder_);
}

ImplicitFixedStepIntegrator::ImplicitFixedStepIntegrator(DeserializingStream & s) :
    FixedStepIntegrator(s) {
  s.version("ImplicitFixedStepIntegrator", 1);
  s.unpack("ImplicitFixedStepIntegrator::rootfinder", rootfinder_);
  s.unpack("ImplicitFixedStepIntegrator::backward_rootfinder", backward_rootfinder_);
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

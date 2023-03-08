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


#include "simulator_impl.hpp"
#include "casadi_misc.hpp"

namespace casadi {

bool has_simulator(const std::string& name) {
  return Simulator::has_plugin(name);
}

void load_simulator(const std::string& name) {
  Simulator::load_plugin(name);
}

std::string doc_simulator(const std::string& name) {
  return Simulator::getPlugin(name).doc;
}

Function simulator(const std::string& name, const std::string& solver, const Function& dae,
    const std::vector<double>& grid, const Dict& opts) {
  // Make sure that dae is sound
  if (dae.has_free()) {
    casadi_error("Cannot create '" + name + "' since " + str(dae.get_free()) + " are free.");
  }
  Simulator* intg = Simulator::getPlugin(solver).creator(name, dae, grid);
  return intg->create_advanced(opts);
}

Function simulator(const std::string& name, const std::string& solver, const SXDict& dae,
    const std::vector<double>& grid, const Dict& opts) {
  return simulator(name, solver, Simulator::map2oracle("dae", dae), grid, opts);
}

Function simulator(const std::string& name, const std::string& solver, const MXDict& dae,
    const std::vector<double>& grid, const Dict& opts) {
  return simulator(name, solver, Simulator::map2oracle("dae", dae), grid, opts);
}

std::vector<std::string> simulator_in() {
  std::vector<std::string> ret(simulator_n_in());
  for (size_t i=0; i<ret.size(); ++i) ret[i]=simulator_in(i);
  return ret;
}

std::vector<std::string> simulator_out() {
  std::vector<std::string> ret(simulator_n_out());
  for (size_t i=0; i<ret.size(); ++i) ret[i]=simulator_out(i);
  return ret;
}

std::string simulator_in(casadi_int ind) {
  switch (static_cast<SimulatorInput>(ind)) {
  case SIMULATOR_X0:  return "x0";
  case SIMULATOR_U:   return "u";
  case SIMULATOR_Z0:  return "z0";
  case SIMULATOR_P:   return "p";
  case SIMULATOR_NUM_IN: break;
  }
  return std::string();
}

std::string simulator_out(casadi_int ind) {
  switch (static_cast<SimulatorOutput>(ind)) {
  case SIMULATOR_X:  return "x";
  case SIMULATOR_Z:  return "z";
  case SIMULATOR_NUM_OUT: break;
  }
  return std::string();
}

Simulator::Simulator(const std::string& name, const Function& oracle,
    const std::vector<double>& grid) : OracleFunction(name, oracle), grid_(grid) {
  // Negative number of parameters for consistancy checking
  np_ = -1;

  // Default options
  print_stats_ = false;
  nondiff_ = true;
  nfwd_ = 0;
}

Simulator::~Simulator() {
}

size_t Simulator::get_n_in() {
  size_t ret = 0;
  // Regular inputs
  ret += SIMULATOR_NUM_IN;
  // Nondifferentiated outputs
  if (!nondiff_ && nfwd_ > 0) ret += SIMULATOR_NUM_OUT;
  // Forward seeds
  if (nfwd_ > 0) ret += SIMULATOR_NUM_IN;
  return ret;
}

size_t Simulator::get_n_out() {
  size_t ret = 0;
  // Regular outputs
  if (nondiff_) ret += SIMULATOR_NUM_OUT;
  // Forward sensitivities
  if (nfwd_ > 0) ret += SIMULATOR_NUM_OUT;
  return ret;
}

Sparsity Simulator::get_sparsity_in(casadi_int i) {
  // Regular inputs
  switch (i) {
    case SIMULATOR_X0: return x();
    case SIMULATOR_U: return repmat(u(), 1, ng_ - 1);
    case SIMULATOR_Z0: return z();
    case SIMULATOR_P: return p();
  }
  i -= SIMULATOR_NUM_IN;
  // Nondifferentiated outputs
  if (!nondiff_ && nfwd_ > 0) {
    switch (i) {
      case SIMULATOR_X: return repmat(Sparsity(x().size()), 1, ng_);
      case SIMULATOR_Z: return repmat(Sparsity(z().size()), 1, ng_);
    }
    i -= SIMULATOR_NUM_OUT;
  }
  // Forward seeds
  if (nfwd_ > 0) {
    switch (i) {
      case SIMULATOR_X0: return repmat(x(), 1, nfwd_);
      case SIMULATOR_U: return repmat(u(), 1, (ng_ - 1) * nfwd_);
      case SIMULATOR_Z0: return repmat(z(), 1, nfwd_);
      case SIMULATOR_P: return repmat(p(), 1, nfwd_);
    }
    i -= SIMULATOR_NUM_IN;
  }
  // Default return
  return Sparsity();
}

Sparsity Simulator::get_sparsity_out(casadi_int i) {
  // Regular outputs
  if (nondiff_) {
    switch (i) {
      case SIMULATOR_X: return repmat(x(), 1, ng_);
      case SIMULATOR_Z: return repmat(z(), 1, ng_);
    }
    i -= SIMULATOR_NUM_OUT;
  }
  // Forward sensitivities
  if (nfwd_ > 0) {
    switch (i) {
      case SIMULATOR_X: return repmat(x(), 1, ng_ * nfwd_);
      case SIMULATOR_Z: return repmat(z(), 1, ng_ * nfwd_);
    }
    i -= SIMULATOR_NUM_OUT;
  }
  // Default return
  return Sparsity();
}

std::string Simulator::get_name_in(casadi_int i) {
  // Regular inputs
  if (i < SIMULATOR_NUM_IN) return simulator_in(i);
  i -= SIMULATOR_NUM_IN;
  // Nondifferentiated outputs
  if (!nondiff_ && nfwd_ > 0) {
    if (i < SIMULATOR_NUM_OUT) return "out_" + simulator_out(i);
    i -= SIMULATOR_NUM_OUT;
  }
  // Forward seeds
  if (nfwd_ > 0) {
    if (i < SIMULATOR_NUM_IN) return "fwd_" + simulator_in(i);
    i -= SIMULATOR_NUM_IN;
  }
  // Default return
  return std::string();
}

std::string Simulator::get_name_out(casadi_int i) {
  // Regular outputs
  if (nondiff_) {
    if (i < SIMULATOR_NUM_OUT) return simulator_out(i);
    i -= SIMULATOR_NUM_OUT;
  }
  // Forward sensitivities
  if (nfwd_ > 0) {
    if (i < SIMULATOR_NUM_OUT) return "fwd_" + simulator_out(i);
    i -= SIMULATOR_NUM_OUT;
  }
  // Default return
  return std::string();
}

Function Simulator::create_advanced(const Dict& opts) {
  return Function::create(this, opts);
}

void Simulator::set_work(void* mem, const double**& arg, double**& res,
    casadi_int*& iw, double*& w) const {
  auto m = static_cast<SimulatorMemory*>(mem);
  // Call the base class method
  OracleFunction::set_work(m, arg, res, iw, w);
  // Inputs
  m->x0 = arg[SIMULATOR_X0];
  m->u = arg[SIMULATOR_U];
  m->z0 = arg[SIMULATOR_Z0];
  m->p = arg[SIMULATOR_P];
  arg += SIMULATOR_NUM_IN;
  // Nondifferentiated outputs, if given
  if (!nondiff_ && nfwd_ > 0) {
    m->out_x = arg[SIMULATOR_X];
    m->out_z = arg[SIMULATOR_Z];
    arg += SIMULATOR_NUM_OUT;
  } else {
    m->out_x = m->out_z = 0;
  }
  // Forward seeds
  if (nfwd_ > 0) {
    m->fwd_x0 = arg[SIMULATOR_X0];
    m->fwd_u = arg[SIMULATOR_U];
    m->fwd_z0 = arg[SIMULATOR_Z0];
    m->fwd_p = arg[SIMULATOR_P];
    arg += SIMULATOR_NUM_IN;
  } else {
    m->fwd_x0 = m->fwd_u = m->fwd_z0 = m->fwd_p = 0;
  }
  // Outputs, if requested
  if (nondiff_) {
    m->x = res[SIMULATOR_X];
    m->z = res[SIMULATOR_Z];
    res += SIMULATOR_NUM_OUT;
  } else {
    m->x = m->z = 0;
  }
  // Forward sensitivities
  if (nfwd_ > 0) {
    m->fwd_x = res[SIMULATOR_X];
    m->fwd_z = res[SIMULATOR_Z];
    res += SIMULATOR_NUM_OUT;
  } else {
    m->fwd_x = m->fwd_z = 0;
  }
  // Current state
  m->xk = w;  w += nx_;
  m->zk = w;  w += nz_;
  m->fwd_xk = w;  w += nfwd_ * nx_;
  m->fwd_zk = w;  w += nfwd_ * nz_;
}

int Simulator::eval(const double** arg, double** res, casadi_int* iw, double* w, void* mem) const {
  auto m = static_cast<SimulatorMemory*>(mem);
  // Setup memory object
  setup(m, arg, res, iw, w);
  // Get initial state, algebraic variable guess
  casadi_copy(m->x0, nx_, m->xk);
  casadi_copy(m->z0, nz_, m->zk);
  // Forward seeds w.r.t. initial state
  casadi_copy(m->fwd_x0, nfwd_ * nx_, m->fwd_xk);
  casadi_copy(m->fwd_z0, nfwd_ * nz_, m->fwd_zk);
  // Reset solver, take time to t0, calculate outputs at t0
  m->t = grid_.front();
  reset(m);
  // Get state
  casadi_copy(m->xk, nx_, m->x);
  casadi_copy(m->zk, nz_, m->z);
  // Get state sensitivities
  for (casadi_int i = 0; i < nfwd_; ++i) {
    casadi_copy(m->fwd_xk + nx_ * i, nx_, m->fwd_x + nx_ * ng_ * i);
    casadi_copy(m->fwd_zk + nz_ * i, nz_, m->fwd_z + nz_ * ng_ * i);
  }
  // Advance output time
  if (m->x) m->x += nx_;
  if (m->z) m->z += nz_;
  // Next stop time due to step change in input
  casadi_int k_stop = next_stop(1, m->u);
  // Integrate forward
  for (casadi_int k = 1; k < ng_; ++k) {
    // Update stopping time, if needed
    if (k > k_stop) k_stop = next_stop(k, m->u);
    // Integrate forward
    if (verbose_) casadi_message("Integrating to " + str(grid_[k])
      + ", stopping time " + str(grid_[k_stop]));
    advance(m, grid_[k], grid_[k_stop]);
    casadi_copy(m->xk, nx_, m->x);
    casadi_copy(m->zk, nz_, m->z);
    // Collect state sensitivities
    for (casadi_int i = 0; i < nfwd_; ++i) {
      casadi_copy(m->fwd_xk + nx_ * i, nx_, m->fwd_x + nx_ * (ng_ * i + k));
      casadi_copy(m->fwd_zk + nz_ * i, nz_, m->fwd_z + nz_ * (ng_ * i + k));
    }
    // Advance output time
    if (m->x) m->x += nx_;
    if (m->u) m->u += nu_;
    if (m->fwd_u) m->fwd_u += nu_;
    if (m->z) m->z += nz_;
  }
  // Print stats
  if (print_stats_) print_stats(m);
  return 0;
}

const Options Simulator::options_
= {{&OracleFunction::options_},
   {{"print_stats",
     {OT_BOOL,
      "Print out statistics after integration"}},
    {"nondiff",
     {OT_BOOL,
      "Output nondifferentiated"}},
    {"nfwd",
     {OT_INT,
      "Number of forward sensitivities"}}
   }
};

void Simulator::init(const Dict& opts) {
  // Read options
  for (auto&& op : opts) {
    if (op.first=="print_stats") {
      print_stats_ = op.second;
    } else if (op.first=="nondiff") {
      nondiff_ = op.second;
    } else if (op.first=="nfwd") {
      nfwd_ = op.second;
    }
  }

  // Number of grid points
  ng_ = grid_.size();

  // Consistency checks
  casadi_assert(nondiff_ || nfwd_ > 0, "Inconsistent options");
  casadi_assert(ng_ >= 2, "Need at least two grid points for simulation");

  // Call the base class method
  OracleFunction::init(opts);

  // Oracle can be evaluated directly
  set_function(oracle_, "dae");

  // Error if sparse input
  casadi_assert(x().is_dense(), "Sparse DAE not supported");
  casadi_assert(u().is_dense(), "Sparse DAE not supported");
  casadi_assert(z().is_dense(), "Sparse DAE not supported");
  casadi_assert(p().is_dense(), "Sparse DAE not supported");

  // Get dimensions (including sensitivity equations)
  nx_ = x().nnz();
  nu_ = u().nnz();
  nz_ = z().nnz();
  np_ = p().nnz();

  // Sensitivities not implemented
  if (nfwd_ > 0) casadi_warning("Forward sensitivities experimental");

  // Nominal values for states
  nom_x_ = oracle_.nominal_in(DYN_X);
  nom_z_ = oracle_.nominal_in(DYN_Z);

  // Get the sparsities of the forward and reverse DAE
  sp_jac_dae_ = sp_jac_dae();
  casadi_assert(!sp_jac_dae_.is_singular(),
    "Jacobian of the forward problem is structurally rank-deficient. "
    "sprank(J)=" + str(sprank(sp_jac_dae_)) + "<" + str(nx_+nz_));

  // Work vectors
  alloc_w(nx_, true);  // xk
  alloc_w(nz_, true);  // zk
  alloc_w(nx_ * nfwd_, true);  // fwd_xk
  alloc_w(nz_ * nfwd_, true);  // fwd_zk
}

int Simulator::init_mem(void* mem) const {
  if (OracleFunction::init_mem(mem)) return 1;

  //auto m = static_cast<SimulatorMemory*>(mem);
  return 0;
}

Sparsity Simulator::sp_jac_dae() {
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

std::map<std::string, Simulator::Plugin> Simulator::solvers_;

const std::string Simulator::infix_ = "simulator";

template<typename XType>
Function Simulator::map2oracle(const std::string& name,
    const std::map<std::string, XType>& d, const Dict& opts) {
  // Gather symbolic inputs and outputs
  std::vector<XType> de_in(enum_traits<DynIn>::n_enum), de_out(enum_traits<DynOut>::n_enum);
  for (auto&& i : d) {
    if (has_enum<DynIn>(i.first)) {
      de_in[to_enum<DynIn>(i.first)] = i.second;
    } else if (has_enum<DynOut>(i.first)) {
      de_out[to_enum<DynOut>(i.first)] = i.second;
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
  // Construct
  return Function(name, de_in, de_out, dyn_in(), dyn_out(), opts);
}

casadi_int Simulator::next_stop(casadi_int k, const double* u) const {
  // Integrate till the end if no input signals
  if (nu_ == 0 || u == 0) return ng_ - 1;
  // Find the next discontinuity, if any
  for (; k + 1 < ng_; ++k) {
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

} // namespace casadi

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

std::string to_string(DynIn v) {
  switch (v) {
  case DYN_T: return "t";
  case DYN_X: return "x";
  case DYN_U: return "u";
  case DYN_Z: return "z";
  case DYN_P: return "p";
  default: break;
  }
  return "";
}

std::string to_string(DynOut v) {
  switch (v) {
  case DYN_ODE: return "ode";
  case DYN_ALG: return "alg";
  case DYN_Y: return "y";
  default: break;
  }
  return "";
}

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
  case SIMULATOR_Y:  return "y";
  case SIMULATOR_NUM_OUT: break;
  }
  return std::string();
}

casadi_int simulator_n_in() {
  return SIMULATOR_NUM_IN;
}

casadi_int simulator_n_out() {
  return SIMULATOR_NUM_OUT;
}

Simulator::Simulator(const std::string& name, const Function& oracle,
    const std::vector<double>& grid) : OracleFunction(name, oracle), grid_(grid) {
  // Negative number of parameters for consistancy checking
  np_ = -1;

  // Default options
  print_stats_ = false;
}

Simulator::~Simulator() {
}

Sparsity Simulator::get_sparsity_in(casadi_int i) {
  switch (static_cast<SimulatorInput>(i)) {
  case SIMULATOR_X0: return x();
  case SIMULATOR_U: return repmat(u(), 1, grid_.size() - 1);
  case SIMULATOR_Z0: return z();
  case SIMULATOR_P: return p();
  case SIMULATOR_NUM_IN: break;
  }
  return Sparsity();
}

Sparsity Simulator::get_sparsity_out(casadi_int i) {
  switch (static_cast<SimulatorOutput>(i)) {
  case SIMULATOR_X: return repmat(x(), 1, grid_.size());
  case SIMULATOR_Z: return repmat(z(), 1, grid_.size());
  case SIMULATOR_Y: return repmat(y(), 1, grid_.size());
  case SIMULATOR_NUM_OUT: break;
  }
  return Sparsity();
}

Function Simulator::create_advanced(const Dict& opts) {
  return Function::create(this, opts);
}

int Simulator::
eval(const double** arg, double** res, casadi_int* iw, double* w, void* mem) const {
  auto m = static_cast<SimulatorMemory*>(mem);
  // Read inputs
  const double* x0 = arg[SIMULATOR_X0];
  const double* u = arg[SIMULATOR_U];
  const double* z0 = arg[SIMULATOR_Z0];
  const double* p = arg[SIMULATOR_P];
  arg += SIMULATOR_NUM_IN;
  // Read outputs
  double* x = res[SIMULATOR_X];
  double* z = res[SIMULATOR_Z];
  double* y = res[SIMULATOR_Y];
  res += SIMULATOR_NUM_OUT;
  // Setup memory object
  setup(m, arg, res, iw, w);
  // Reset solver, take time to t0, calculate outputs at t0
  reset(m, grid_.front(), x0, u, z0, p, y);
  if (x) x += nx_;
  if (z) z += nz_;
  if (y) y += ny_;
  // Integrate forward
  for (casadi_int k = 1; k < grid_.size(); ++k) {
    // Integrate forward
    advance(m, grid_[k], x, u, z, p, y);
    if (x) x += nx_;
    if (u) u += nu_;
    if (z) z += nz_;
    if (y) y += ny_;
  }
  if (print_stats_) print_stats(m);

  return 0;
}

const Options Simulator::options_
= {{&OracleFunction::options_},
   {{"print_stats",
     {OT_BOOL,
      "Print out statistics after integration"}}
   }
};

void Simulator::init(const Dict& opts) {
  // Read options
  for (auto&& op : opts) {
    if (op.first=="print_stats") {
      print_stats_ = op.second;
    }
  }

  // Store a copy of the options, for creating augmented simulators
  opts_ = opts;

  // Consistency check
  casadi_assert(grid_.size() >= 2, "Need at least two grid points for simulation");

  // Call the base class method
  OracleFunction::init(opts);

  // For sparsity pattern propagation
  set_function(oracle_);

  // Error if sparse input
  casadi_assert(x().is_dense(), "Sparse DAE not supported");
  casadi_assert(u().is_dense(), "Sparse DAE not supported");
  casadi_assert(z().is_dense(), "Sparse DAE not supported");
  casadi_assert(p().is_dense(), "Sparse DAE not supported");

  // Get dimensions (excluding sensitivity equations)
  nx1_ = x().size1();
  nu1_ = u().size1();
  nz1_ = z().size1();
  np1_ = p().size1();
  ny1_ = y().size1();

  // Get dimensions (including sensitivity equations)
  nx_ = x().nnz();
  nu_ = u().nnz();
  nz_ = z().nnz();
  np_ = p().nnz();
  ny_ = y().nnz();

  // Number of sensitivities
  ns_ = x().size2()-1;

  // Get the sparsities of the forward and reverse DAE
  sp_jac_dae_ = sp_jac_dae();
  casadi_assert(!sp_jac_dae_.is_singular(),
                        "Jacobian of the forward problem is structurally rank-deficient. "
                        "sprank(J)=" + str(sprank(sp_jac_dae_)) + "<"
                        + str(nx_+nz_));

  // Allocate sufficiently large work vectors
  alloc_w(nx_ + nz_);
  alloc_w(nx_ + nz_, true);
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

void Simulator::setStopTime(SimulatorMemory* mem, double tf) const {
  casadi_error("setStopTime not defined for class " + class_name());
}

template<typename XType>
Function Simulator::map2oracle(const std::string& name,
    const std::map<std::string, XType>& d, const Dict& opts) {
  // Gather symbolic inputs and outputs
  std::vector<XType> de_in(DYN_NUM_IN), de_out(DYN_NUM_OUT);
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
  return Function(name, de_in, de_out, enum_names<DynIn>(), enum_names<DynOut>(), opts);
}

void Simulator::eval_y(SimulatorMemory* mem, double t, const double* x, const double* u,
    const double* z, const double* p, double* y) const {
  // Calculate outputs
  std::fill_n(mem->arg, DYN_NUM_IN, nullptr);
  mem->arg[DYN_T] = &t;
  mem->arg[DYN_X] = x;
  mem->arg[DYN_U] = u;
  mem->arg[DYN_Z] = z;
  mem->arg[DYN_P] = p;
  std::fill_n(mem->res, DYN_NUM_OUT, nullptr);
  mem->res[DYN_Y] = y;
  if (calc_function(mem, "dae")) casadi_error("'dae' calculation failed");
}

} // namespace casadi

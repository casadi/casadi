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


#include "fmu_impl.hpp"
#include "fmu_function.hpp"
#include "dae_builder_internal.hpp"

#ifdef WITH_FMI2
#include "fmu2.hpp"
#endif  // WITH_FMI2

#ifdef WITH_FMI3
#include "fmu3.hpp"
#endif  // WITH_FMI3

namespace casadi {

// Throw informative error message
#define THROW_ERROR(FNAME, WHAT) \
throw CasadiException("Error in Fmu::" FNAME " for '" + this->name() + "' "\
  "[" + this->class_name() + "] at " + CASADI_WHERE + ":\n"\
  + std::string(WHAT));


Fmu::Fmu() {
}

Fmu::Fmu(const std::string& name, FmuApi api, const DaeBuilderInternal* dae,
    const std::vector<std::string>& scheme_in,
    const std::vector<std::string>& scheme_out,
    const std::map<std::string, std::vector<size_t>>& scheme,
    const std::vector<std::string>& aux) {
  if (api == FmuApi::FMI2) {
#ifdef WITH_FMI2
  // Create
  own(new Fmu2(name, scheme_in, scheme_out, scheme, aux));
#else  // WITH_FMI2
  // No compilation support
  casadi_error("CasADi was not compiled with WITH_FMI2=ON.");
#endif  // WITH_FMI2
  } else if (api == FmuApi::FMI3) {
#ifdef WITH_FMI3
  // Create
  own(new Fmu3(name, scheme_in, scheme_out, scheme, aux));
#else  // WITH_FMI3
  // No compilation support
  casadi_error("CasADi was not compiled with WITH_FMI3=ON.");
#endif  // WITH_FMI3
  } else {
    // Not supported
    casadi_error("Unsupported FMU API: " + to_string(api));
  }
  // Initialize
  try {
    (*this)->init(dae);
    (*this)->finalize();
  } catch(std::exception& e) {
    THROW_ERROR("init", e.what());
  }
}

FmuInternal* Fmu::operator->() {
  return static_cast<FmuInternal*>(SharedObject::operator->());
}

const FmuInternal* Fmu::operator->() const {
  return static_cast<const FmuInternal*>(SharedObject::operator->());
}

FmuInternal* Fmu::get() const {
  return static_cast<FmuInternal*>(SharedObject::get());
}

const std::string& Fmu::name() const {
  if (is_null()) {
    static std::string null = "null";
    return null;
  } else {
    return (*this)->name_;
  }
}

const std::string& Fmu::instance_name() const {
  if (is_null()) {
    static std::string null = "null";
    return null;
  } else {
    return (*this)->instance_name_;
  }
}

size_t Fmu::n_in() const {
  try {
    return (*this)->n_in();
  } catch(std::exception& e) {
    THROW_ERROR("n_in", e.what());
  }
}

size_t Fmu::n_out() const {
  try {
    return (*this)->n_out();
  } catch(std::exception& e) {
    THROW_ERROR("n_out", e.what());
  }
}

size_t Fmu::index_in(const std::string& n) const {
  try {
    return (*this)->index_in(n);
  } catch(std::exception& e) {
    THROW_ERROR("index_in", e.what());
  }
}

size_t Fmu::index_out(const std::string& n) const {
  try {
    return (*this)->index_out(n);
  } catch(std::exception& e) {
    THROW_ERROR("index_out", e.what());
  }
}

const std::vector<size_t>& Fmu::ired(size_t ind) const {
  try {
    return (*this)->ired_.at(ind);
  } catch(std::exception& e) {
    THROW_ERROR("ired", e.what());
  }
}

const std::vector<size_t>& Fmu::ored(size_t ind) const {
  try {
    return (*this)->ored_.at(ind);
  } catch(std::exception& e) {
    THROW_ERROR("ored", e.what());
  }
}

double Fmu::nominal_in(size_t ind) const {
  try {
    return (*this)->nominal_in_.at(ind);
  } catch(std::exception& e) {
    THROW_ERROR("nominal_in", e.what());
  }
}

double Fmu::nominal_out(size_t ind) const {
  try {
    return (*this)->nominal_out_.at(ind);
  } catch(std::exception& e) {
    THROW_ERROR("nominal_out", e.what());
  }
}

double Fmu::min_in(size_t ind) const {
  try {
    return (*this)->min_in_.at(ind);
  } catch(std::exception& e) {
    THROW_ERROR("min_in", e.what());
  }
}

double Fmu::max_in(size_t ind) const {
  try {
    return (*this)->max_in_.at(ind);
  } catch(std::exception& e) {
    THROW_ERROR("max_in", e.what());
  }
}

std::vector<double> Fmu::all_nominal_in(size_t ind) const {
  try {
    return (*this)->all_nominal_in(ind);
  } catch(std::exception& e) {
    THROW_ERROR("all_nominal_in", e.what());
  }
}

std::vector<double> Fmu::all_nominal_out(size_t ind) const {
  try {
    return (*this)->all_nominal_out(ind);
  } catch(std::exception& e) {
    THROW_ERROR("all_nominal_out", e.what());
  }
}

std::string Fmu::desc_in(FmuMemory* m, size_t id, bool more) const {
  try {
    return (*this)->desc_in(m, id, more);
  } catch(std::exception& e) {
    THROW_ERROR("desc_in", e.what());
  }
}

bool Fmu::provides_directional_derivatives() const {
  try {
    return (*this)->provides_directional_derivatives_;
  } catch(std::exception& e) {
    THROW_ERROR("provides_directional_derivatives", e.what());
  }
}

bool Fmu::provides_adjoint_derivatives() const {
  try {
    return (*this)->provides_adjoint_derivatives_;
  } catch(std::exception& e) {
    THROW_ERROR("provides_adjoint_derivatives", e.what());
  }
}

Sparsity Fmu::jac_sparsity(const std::vector<size_t>& osub,
    const std::vector<size_t>& isub) const {
  try {
    return (*this)->jac_sparsity(osub, isub);
  } catch(std::exception& e) {
    THROW_ERROR("jac_sparsity", e.what());
  }
}

Sparsity Fmu::hess_sparsity(const std::vector<size_t>& r, const std::vector<size_t>& c) const {
  try {
    return (*this)->hess_sparsity(r, c);
  } catch(std::exception& e) {
    THROW_ERROR("hess_sparsity", e.what());
  }
}

int Fmu::init_mem(FmuMemory* m) const {
  try {
    return (*this)->init_mem(m);
  } catch(std::exception& e) {
    THROW_ERROR("init_mem", e.what());
    return 1;
  }
}

void Fmu::free_instance(void* instance) const {
  try {
    return (*this)->free_instance(instance);
  } catch(std::exception& e) {
    THROW_ERROR("free_instance", e.what());
  }
}

void Fmu::set(FmuMemory* m, size_t ind, const double* value) const {
  try {
    return (*this)->set(m, ind, value);
  } catch(std::exception& e) {
    THROW_ERROR("set", e.what());
  }
}

void Fmu::request(FmuMemory* m, size_t ind) const {
  try {
    return (*this)->request(m, ind);
  } catch(std::exception& e) {
    THROW_ERROR("request", e.what());
  }
}

int Fmu::eval(FmuMemory* m) const {
  try {
    return (*this)->eval(m);
  } catch(std::exception& e) {
    THROW_ERROR("eval", e.what());
  }
}

void Fmu::get(FmuMemory* m, size_t id, double* value) const {
  try {
    return (*this)->get(m, id, value);
  } catch(std::exception& e) {
    THROW_ERROR("get", e.what());
  }
}

void Fmu::set_fwd(FmuMemory* m, casadi_int nseed, const casadi_int* id, const double* v) const {
  try {
    return (*this)->set_fwd(m, nseed, id, v);
  } catch(std::exception& e) {
    THROW_ERROR("set_fwd", e.what());
  }
}

void Fmu::set_fwd(FmuMemory* m, size_t ind, const double* v) const {
  try {
    return (*this)->set_fwd(m, ind, v);
  } catch(std::exception& e) {
    THROW_ERROR("set_fwd", e.what());
  }
}

void Fmu::request_fwd(FmuMemory* m, casadi_int nsens, const casadi_int* id,
    const casadi_int* wrt_id) const {
  try {
    return (*this)->request_fwd(m, nsens, id, wrt_id);
  } catch(std::exception& e) {
    THROW_ERROR("request_fwd", e.what());
  }
}

void Fmu::request_fwd(FmuMemory* m, casadi_int ind) const {
  try {
    return (*this)->request_fwd(m, ind);
  } catch(std::exception& e) {
    THROW_ERROR("request_fwd", e.what());
  }
}

int Fmu::eval_fwd(FmuMemory* m, bool independent_seeds) const {
  try {
    return (*this)->eval_fwd(m, independent_seeds);
  } catch(std::exception& e) {
    THROW_ERROR("eval_fwd", e.what());
  }
}

void Fmu::get_fwd(FmuMemory* m, casadi_int nsens, const casadi_int* id, double* v) const {
  try {
    return (*this)->get_fwd(m, nsens, id, v);
  } catch(std::exception& e) {
    THROW_ERROR("get_fwd", e.what());
  }
}

void Fmu::get_fwd(FmuMemory* m, size_t ind, double* v) const {
  try {
    return (*this)->get_fwd(m, ind, v);
  } catch(std::exception& e) {
    THROW_ERROR("get_fwd", e.what());
  }
}

void Fmu::set_adj(FmuMemory* m, casadi_int nseed, const casadi_int* id, const double* v) const {
  try {
    return (*this)->set_adj(m, nseed, id, v);
  } catch(std::exception& e) {
    THROW_ERROR("set_adj", e.what());
  }
}

void Fmu::set_adj(FmuMemory* m, size_t ind, const double* v) const {
  try {
    return (*this)->set_adj(m, ind, v);
  } catch(std::exception& e) {
    THROW_ERROR("set_adj", e.what());
  }
}

void Fmu::request_adj(FmuMemory* m, casadi_int nsens, const casadi_int* id,
    const casadi_int* wrt_id) const {
  try {
    return (*this)->request_adj(m, nsens, id, wrt_id);
  } catch(std::exception& e) {
    THROW_ERROR("request_adj", e.what());
  }
}

void Fmu::request_adj(FmuMemory* m, casadi_int ind) const {
  try {
    return (*this)->request_adj(m, ind);
  } catch(std::exception& e) {
    THROW_ERROR("request_adj", e.what());
  }
}

int Fmu::eval_adj(FmuMemory* m) const {
  try {
    return (*this)->eval_adj(m);
  } catch(std::exception& e) {
    THROW_ERROR("eval_adj", e.what());
  }
}

void Fmu::get_adj(FmuMemory* m, casadi_int nsens, const casadi_int* id, double* v) const {
  try {
    return (*this)->get_adj(m, nsens, id, v);
  } catch(std::exception& e) {
    THROW_ERROR("get_adj", e.what());
  }
}

void Fmu::get_adj(FmuMemory* m, size_t ind, double* v) const {
  try {
    return (*this)->get_adj(m, ind, v);
  } catch(std::exception& e) {
    THROW_ERROR("get_adj", e.what());
  }
}

void Fmu::get_stats(FmuMemory* m, Dict* stats,
    const std::vector<std::string>& name_in, const InputStruct* in) const {
  try {
    return (*this)->get_stats(m, stats, name_in, in);
  } catch(std::exception& e) {
    THROW_ERROR("get_stats", e.what());
  }
}

FmuInternal::FmuInternal(const std::string& name,
    const std::vector<std::string>& scheme_in,
    const std::vector<std::string>& scheme_out,
    const std::map<std::string, std::vector<size_t>>& scheme,
    const std::vector<std::string>& aux)
    : name_(name), scheme_in_(scheme_in), scheme_out_(scheme_out), scheme_(scheme), aux_(aux) {
}

FmuInternal::~FmuInternal() {
}

void FmuInternal::init(const DaeBuilderInternal* dae) {
  // Copy info from DaeBuilder
  fmutol_ = dae->fmutol_;
  instance_name_ = dae->model_identifier_;
  instantiation_token_ = dae->instantiation_token_;
  logging_on_ = dae->debug_;
  number_of_event_indicators_ = dae->number_of_event_indicators_;
  provides_directional_derivatives_ = dae->provides_directional_derivatives_;
  provides_adjoint_derivatives_ = dae->provides_adjoint_derivatives_;

  // Path to resource directory
  resource_loc_ = "file://" + dae->path_ + "/resources";

  // Load DLL
  std::string instance_name_no_dot = dae->model_identifier_;
  std::replace(instance_name_no_dot.begin(), instance_name_no_dot.end(), '.', '_');
  std::string dll_path = dae->path_ + "/binaries/" + system_infix()
    + "/" + instance_name_no_dot + dll_suffix();
  li_ = Importer(dll_path, "dll");

  // Mark input indices
  size_t numel = 0;
  std::vector<bool> lookup(dae->n_variables(), false);
  for (auto&& n : scheme_in_) {
    for (size_t i : scheme_.at(n)) {
      casadi_assert(!lookup.at(i), "Duplicate variable: " + dae->variable(i).name);
      lookup.at(i) = true;
      numel++;
    }
  }
  // Input mappings
  iind_.reserve(numel);
  iind_map_.reserve(lookup.size());
  for (size_t k = 0; k < lookup.size(); ++k) {
    if (lookup[k]) {
      iind_map_.push_back(iind_.size());
      iind_.push_back(k);
    } else {
      iind_map_.push_back(-1);
    }
  }
  // Mark output indices
  numel = 0;
  std::fill(lookup.begin(), lookup.end(), false);
  for (auto&& n : scheme_out_) {
    for (size_t i : scheme_.at(n)) {
      casadi_assert(!lookup.at(i), "Duplicate variable: " + dae->variable(i).name);
      lookup.at(i) = true;
      numel++;
    }
  }
  // Construct mappings
  oind_.reserve(numel);
  oind_map_.reserve(lookup.size());
  for (size_t k = 0; k < lookup.size(); ++k) {
    if (lookup[k]) {
      oind_map_.push_back(oind_.size());
      oind_.push_back(k);
    } else {
      oind_map_.push_back(-1);
    }
  }
  // Inputs
  ired_.resize(scheme_in_.size());
  for (size_t i = 0; i < ired_.size(); ++i) {
    auto&& s = scheme_.at(scheme_in_[i]);
    ired_[i].resize(s.size());
    for (size_t k = 0; k < s.size(); ++k) {
      ired_[i][k] = iind_map_.at(s[k]);
    }
  }
  // Outputs
  ored_.resize(scheme_out_.size());
  for (size_t i = 0; i < ored_.size(); ++i) {
    auto&& s = scheme_.at(scheme_out_[i]);
    ored_[i].resize(s.size());
    for (size_t k = 0; k < s.size(); ++k) {
      ored_[i][k] = oind_map_.at(s[k]);
    }
  }

  // Collect meta information for inputs
  nominal_in_.reserve(iind_.size());
  min_in_.reserve(iind_.size());
  max_in_.reserve(iind_.size());
  vn_in_.reserve(iind_.size());
  vr_in_.reserve(iind_.size());
  for (size_t i : iind_) {
    const Variable& v = dae->variable(i);
    nominal_in_.push_back(v.nominal);
    min_in_.push_back(v.min);
    max_in_.push_back(v.max);
    vn_in_.push_back(v.name);
    vr_in_.push_back(v.value_reference);
  }
  // Collect meta information for outputs
  nominal_out_.reserve(oind_.size());
  min_out_.reserve(oind_.size());
  max_out_.reserve(oind_.size());
  vn_out_.reserve(oind_.size());
  vr_out_.reserve(oind_.size());
  for (size_t i : oind_) {
    const Variable& v = dae->variable(i);
    nominal_out_.push_back(v.nominal);
    min_out_.push_back(v.min);
    max_out_.push_back(v.max);
    vn_out_.push_back(v.name);
    vr_out_.push_back(v.value_reference);
  }

  // Numerical values for inputs
  value_in_.resize(iind_.size());

  // Get Jacobian sparsity information
  jac_sp_ = dae->jac_sparsity(oind_, iind_);

  // Get Hessian sparsity information
  hess_sp_ = dae->hess_sparsity(oind_, iind_);
}

int FmuInternal::get_adjoint_derivative(void* instance, const unsigned int* vr_out, size_t n_out,
    const unsigned int* vr_in, size_t n_in, const double* seed, size_t n_seed,
    double* sensitivity, size_t n_sensitivity) const {
  casadi_error("Adjoint derivatives not supported for " + class_name());
  return 1;
}

void FmuInternal::finalize() {
  // Get FMI C functions
  load_functions();

  // Create a temporary instance
  void* c = instantiate();
  // Set all values
  if (set_values(c)) {
    casadi_error("FmuInternal::set_values failed");
  }
  // Initialization mode begins
  if (enter_initialization_mode(c)) {
    casadi_error("FmuInternal::enter_initialization_mode failed");
  }
  // Get input values
  if (get_real(c, get_ptr(vr_in_), vr_in_.size(), get_ptr(value_in_), value_in_.size())) {
    casadi_error("FmuInternal::get_in failed");
  }
  // Get auxilliary variables
  if (get_aux(c)) {
    casadi_error("FmuInternal::get_aux failed");
  }
  // Free memory
  free_instance(c);
}

void FmuInternal::disp(std::ostream& stream, bool more) const {
  (void)more;  // unused
  stream << name_ << " " << class_name();
}

std::string to_string(FmuApi v) {
  switch (v) {
  case FmuApi::FMI2: return "fmi2";
  case FmuApi::FMI3: return "fmi3";
  default: break;
  }
  return "";
}

size_t FmuInternal::index_in(const std::string& n) const {
  // Linear search for the input
  for (size_t i = 0; i < scheme_in_.size(); ++i) {
    if (scheme_in_[i] == n) return i;
  }
  // Not found
  casadi_error("No such input: " + n);
  return -1;
}

size_t FmuInternal::index_out(const std::string& n) const {
  // Linear search for the input
  for (size_t i = 0; i < scheme_out_.size(); ++i) {
    if (scheme_out_[i] == n) return i;
  }
  // Not found
  casadi_error("No such output: " + n);
  return -1;
}

Sparsity FmuInternal::jac_sparsity(const std::vector<size_t>& osub,
    const std::vector<size_t>& isub) const {
  // Convert to casadi_int type
  std::vector<casadi_int> osub1(osub.begin(), osub.end());
  std::vector<casadi_int> isub1(isub.begin(), isub.end());
  // Index mapping (not used)
  std::vector<casadi_int> mapping;
  // Get selection
  return jac_sp_.sub(osub1, isub1, mapping);
}

Sparsity FmuInternal::hess_sparsity(const std::vector<size_t>& r,
    const std::vector<size_t>& c) const {
  // Convert to casadi_int type
  std::vector<casadi_int> r1(r.begin(), r.end());
  std::vector<casadi_int> c1(c.begin(), c.end());
  // Index mapping (not used)
  std::vector<casadi_int> mapping;
  // Get selection
  return hess_sp_.sub(r1, c1, mapping);
}

std::vector<double> FmuInternal::all_nominal_in(size_t i) const {
  auto&& ind = ired_.at(i);
  std::vector<double> n;
  n.reserve(ind.size());
  for (size_t k : ind) n.push_back(nominal_in_.at(k));
  return n;
}

std::vector<double> FmuInternal::all_nominal_out(size_t i) const {
  auto&& ind = ored_.at(i);
  std::vector<double> n;
  n.reserve(ind.size());
  for (size_t k : ind) n.push_back(nominal_out_.at(k));
  return n;
}

std::string FmuInternal::dll_suffix() {
#if defined(_WIN32)
  // Windows system
  return ".dll";
#elif defined(__APPLE__)
  // OSX
  return ".dylib";
#else
  // Linux
  return ".so";
#endif
}

std::string FmuInternal::desc_in(FmuMemory* m, size_t id, bool more) const {
  // Create description
  if (more) {
    // Detailed description
    std::stringstream ss;
    ss << vn_in_[id] << " = " << m->ibuf_[id] << " (nominal " << nominal_in_[id]
    << ", min " << min_in_[id] << ", max " << max_in_[id] << ")";
    return ss.str();
  } else {
    return vn_in_[id];
  }
}

int FmuInternal::eval_fwd(FmuMemory* m, bool independent_seeds) const {
  // Gather input and output indices
  gather_fwd(m);
  // Calculate derivatives using FMU directional derivative support
  if (m->self.uses_directional_derivatives_) {
    // Evaluate using AD
    if (eval_ad(m)) return 1;
  }
  // Calculate derivatives using finite differences
  if (!m->self.uses_directional_derivatives_ || m->self.validate_forward_) {
    // Evaluate using FD
    if (eval_fd(m, independent_seeds)) return 1;
  }
  return 0;
}

int FmuInternal::eval_adj(FmuMemory* m) const {
  // Gather input and output indices
  gather_adj(m);
  // Quick return if nothing to be calculated
  if (m->id_in_.size() == 0) return 0;
  // Evaluate adjoint derivatives
  if (get_adjoint_derivative(m->instance,
      get_ptr(m->vr_out_), m->id_out_.size(),
      get_ptr(m->vr_in_), m->id_in_.size(),
      get_ptr(m->d_out_), m->id_out_.size(),
      get_ptr(m->d_in_), m->id_in_.size())) {
    casadi_warning("FMU adjoint derivative failed");
    return 1;
  }
  // Collect requested variables
  auto it = m->d_in_.begin();
  for (size_t id : m->id_in_) {
    m->isens_[id] = *it++;
  }
  // Successful return
  return 0;
}

int FmuInternal::eval_ad(FmuMemory* m) const {
  // Number of inputs and outputs
  size_t n_known = m->id_in_.size();
  size_t n_unknown = m->id_out_.size();
  // Quick return if nothing to be calculated
  if (n_unknown == 0) return 0;
  // Evalute (should not be necessary)
  if (get_real(m->instance, get_ptr(m->vr_out_), n_unknown, get_ptr(m->v_out_), n_unknown)) {
    casadi_warning("FMU evaluation failed");
    return 1;
  }
  // Evaluate directional derivatives
  if (get_directional_derivative(m->instance,
      get_ptr(m->vr_out_), n_unknown,
      get_ptr(m->vr_in_), n_known,
      get_ptr(m->d_in_), n_known,
      get_ptr(m->d_out_), n_unknown)) {
    casadi_warning("FMU directional derivative failed");
    return 1;
  }
  // Collect requested variables
  auto it = m->d_out_.begin();
  for (size_t id : m->id_out_) {
    m->osens_[id] = *it++;
  }
  // Successful return
  return 0;
}

int FmuInternal::eval_fd(FmuMemory* m, bool independent_seeds) const {
  // Number of inputs and outputs
  size_t n_known = m->id_in_.size();
  size_t n_unknown = m->id_out_.size();
  // Quick return if nothing to be calculated
  if (n_unknown == 0) return 0;
  // Evalute (should not be necessary)
  if (get_real(m->instance, get_ptr(m->vr_out_), n_unknown, get_ptr(m->v_out_), n_unknown)) {
    casadi_warning("Evaluating FMU failed");
    return 1;
  }
  // Make outputs dimensionless
  for (size_t k = 0; k < n_unknown; ++k) m->v_out_[k] /= nominal_out_[m->id_out_[k]];
  // Number of points in FD stencil
  casadi_int n_points = n_fd_points(m->self.fd_);
  // Offset for points
  casadi_int offset = fd_offset(m->self.fd_);
  // Memory for perturbed outputs
  m->fd_out_.resize(n_points * n_unknown);
  // Which inputs are in bounds
  m->in_bounds_.resize(n_known);
  // Memory for perturbed inputs
  m->v_pert_.resize(n_known);
  // Do any any inputs need flipping?
  m->flip_.resize(n_known);
  size_t first_flip = -1;
  for (size_t i = 0; i < n_known; ++i) {
    // Try to take step
    double test = m->v_in_[i] + m->self.step_ * m->d_in_[i];
    // Check if in bounds
    size_t id = m->id_in_[i];
    if (test >= min_in_[id] && test <= max_in_[id]) {
      // Positive perturbation is fine
      m->flip_[i] = false;
    } else {
      // Try negative direction instead
      test = m->v_in_[i] - m->self.step_ * m->d_in_[i];
      casadi_assert(test >= min_in_[id] && test <= max_in_[id],
        "Cannot perturb " + vn_in_[id] + " at " + str(m->v_in_[i]) + ", min " + str(min_in_[id])
        + ", max " + str(max_in_[id]) + ", nominal " + str(nominal_in_[id]));
      m->flip_[i] = true;
      if (first_flip == size_t(-1)) first_flip = i;
    }
  }
  // If seeds are not independent, we have to flip the sign for all of the seeds or none
  if (first_flip != size_t(-1) && !independent_seeds) {
    // Flip the rest of the seeds
    for (size_t i = 0; i < n_known; ++i) {
      if (!m->flip_[i]) {
        // Test negative direction
        double test = m->v_in_[i] - m->self.step_ * m->d_in_[i];
        size_t id = m->id_in_[i];
        casadi_assert(test >= min_in_[id] && test <= max_in_[id],
          "Cannot perturb both " + vn_in_[id] + " and " + vn_in_[first_flip]);
        // Flip it too
        m->flip_[i] = true;
      }
    }
  }
  // All perturbed outputs
  const double* yk_all[5] = {0};

  // Calculate all perturbed outputs
  for (casadi_int k = 0; k < n_points; ++k) {
    // Where to save the perturbed outputs
    double* yk = &m->fd_out_[n_unknown * k];
    casadi_assert_dev(k < 5);
    yk_all[k] = yk;
    // If unperturbed output, quick return
    if (k == offset) {
      casadi_copy(get_ptr(m->v_out_), n_unknown, yk);
      continue;
    }
    // Perturbation size
    double pert = (k - offset) * m->self.step_;
    // Perturb inputs, if allowed
    for (size_t i = 0; i < n_known; ++i) {
      // Try to take step
      double sign = m->flip_[i] ? -1 : 1;
      double test = m->v_in_[i] + pert * sign * m->d_in_[i];
      // Check if in bounds
      size_t id = m->id_in_[i];
      m->in_bounds_[i] = test >= min_in_[id] && test <= max_in_[id];
      // Take step, if allowed
      m->v_pert_[i] = m->in_bounds_[i] ? test : m->v_in_[i];
    }
    // Pass perturbed inputs to FMU
    if (set_real(m->instance, get_ptr(m->vr_in_), n_known, get_ptr(m->v_pert_), n_known)) {
      casadi_warning("Setting FMU variables failed");
      return 1;
    }
    // Evaluate perturbed FMU
    if (get_real(m->instance, get_ptr(m->vr_out_), n_unknown, yk, n_unknown)) {
      casadi_warning("Evaluation failed");
      return 1;
    }
    // Post-process yk
    for (size_t i = 0; i < n_unknown; ++i) {
      // Variable id
      size_t id = m->id_out_[i];
      // Differentiation with respect to what variable
      size_t wrt_id = m->wrt_.at(id);
      // Find the corresponding input variable
      size_t wrt_i;
      for (wrt_i = 0; wrt_i < n_known; ++wrt_i) {
        if (m->id_in_[wrt_i] == wrt_id) break;
      }
      // Check if in bounds
      if (m->in_bounds_.at(wrt_i)) {
        // Input was in bounds: Keep output, make dimensionless
        yk[i] /= nominal_out_[m->id_out_[i]];
      } else {
        // Input was out of bounds: Discard output
        yk[i] = nan;
      }
    }
  }
  // Restore FMU inputs
  if (set_real(m->instance, get_ptr(m->vr_in_), n_known, get_ptr(m->v_in_), n_known)) {
    casadi_warning("Setting FMU variables failed");
    return 1;
  }
  // Step size
  double h = m->self.step_;

  // Calculate FD approximation
  finite_diff(m->self.fd_, yk_all, get_ptr(m->d_out_), h, n_unknown, eps);

  // Collect requested variables
  for (size_t ind = 0; ind < m->id_out_.size(); ++ind) {
    // Variable id
    size_t id = m->id_out_[ind];
    // With respect to what variable
    size_t wrt = m->wrt_[id];
    // Find the corresponding input variable
    size_t wrt_i;
    for (wrt_i = 0; wrt_i < n_known; ++wrt_i) {
      if (m->id_in_[wrt_i] == wrt) break;
    }
    // Nominal value
    double n = nominal_out_[id];
    // Get the value
    double d_fd = m->d_out_[ind] * n;
    // Correct sign, if necessary
    if (m->flip_[wrt_i]) d_fd = -d_fd;
    // Use FD instead of AD or to compare with AD
    if (m->self.validate_forward_) {
      // Value to compare with
      double d_ad = m->osens_[id];
      // Nominal value used as seed
      d_ad /= nominal_in_[wrt];
      d_fd /= nominal_in_[wrt];
      // Is it a not a number?
      bool d_is_nan = d_ad != d_ad;
      // Magnitude of derivatives
      double d_max = std::fmax(std::fabs(d_fd), std::fabs(d_ad));
      // Check if NaN or error exceeds thresholds
      if (d_is_nan || (d_max > n * m->self.abstol_
          && std::fabs(d_ad - d_fd) > d_max * m->self.reltol_)) {
        // Offset for printing the stencil
        double off = m->fd_out_.at(ind + offset * n_unknown);
        // Warning or add to file
        std::stringstream ss;
        if (m->self.validate_ad_file_.empty()) {
          // Issue warning
          ss << (d_is_nan ? "NaN" : "Inconsistent") << " derivatives of " << vn_out_[id]
            << " w.r.t. " << desc_in(m, wrt) << ", got " << d_ad
            << " for AD vs. " << d_fd << " for FD[" << to_string(m->self.fd_) << "].";
          // Print the stencil:
          ss << "\nValues for step size " << h << ": " << (n * off) << " + [";
          for (casadi_int k = 0; k < n_points; ++k) {
            if (k > 0) ss << ", ";
            ss << (n * (m->fd_out_.at(ind + k * n_unknown) - off));
          }
          ss << "]";
          // Issue warning
          casadi_warning(ss.str());
        } else {
          // Output
          ss << vn_out_[id] << " ";
          // Input
          ss << vn_in_[wrt] << " ";
          // Value
          ss << m->ibuf_[wrt] << " ";
          // Noninal
          ss << nominal_in_[wrt] << " ";
          // Min
          ss << min_in_[wrt] << " ";
          // Max
          ss << max_in_[wrt] << " ";
          // AD
          ss << d_ad << " ";
          // FD
          ss << d_fd << " ";
          // Step
          ss << h << " ";
          // Offset
          ss << off << " ";
          // Stencil
          ss << "[";
          for (casadi_int k = 0; k < n_points; ++k) {
            if (k > 0) ss << ",";
            ss << (n * (m->fd_out_.at(ind + k * n_unknown) - off));
          }
          ss << "]" << std::endl;
          // Append to file
          std::ofstream valfile;
          valfile.open(m->self.validate_ad_file_, std::ios_base::app);
          valfile << ss.str();
        }
      }
    } else {
      // Use instead of AD
      m->osens_[id] = d_fd;
    }
  }
  // Successful return
  return 0;
}

void FmuInternal::get_fwd(FmuMemory* m, casadi_int nsens, const casadi_int* id, double* v) const {
  for (casadi_int i = 0; i < nsens; ++i) {
    *v++ = m->osens_.at(*id++);
  }
}

void FmuInternal::get_fwd(FmuMemory* m, size_t ind, double* v) const {
  // Quick return if not needed
  if (!v) return;
  // Retrieve all sensitivities FIXME(@jaeandersson): should use compatible types
  for (size_t id : ored_[ind]) {
    casadi_int id2 = id;
    get_fwd(m, 1, &id2, v++);
  }
}

void FmuInternal::set_adj(FmuMemory* m, casadi_int nseed,
    const casadi_int* id, const double* v) const {
  for (casadi_int i = 0; i < nseed; ++i) {
    m->osens_.at(*id) = *v++;
    m->omarked_.at(*id) = true;
    id++;
  }
}

void FmuInternal::set_adj(FmuMemory* m, size_t ind, const double* v) const {
  // If seeds are zero, no need to add to seed buffers
  if (!v) return;
  // Pass all seeds FIXME(@jaeandersson): should use compatible types
  for (size_t id : ored_[ind]) {
    casadi_int id2 = id;
    set_adj(m, 1, &id2, v++);
  }
}

void FmuInternal::request_adj(FmuMemory* m, casadi_int nsens, const casadi_int* id,
    const casadi_int* wrt_id) const {
  for (casadi_int i = 0; i < nsens; ++i) {
    m->imarked_.at(*id) = true;
    m->wrt_.at(*id) = *wrt_id++;
    id++;
  }
}

void FmuInternal::request_adj(FmuMemory* m, casadi_int ind) const {
  // Request all sensitivities FIXME(@jaeandersson): should use compatible types
  casadi_int wrt_id = -1;
  for (size_t id : ired_[ind]) {
    casadi_int id2 = id;
    request_adj(m, 1, &id2, &wrt_id);
  }
}

void FmuInternal::get_adj(FmuMemory* m, casadi_int nsens, const casadi_int* id, double* v) const {
  for (casadi_int i = 0; i < nsens; ++i) {
    *v++ = m->isens_.at(*id++);
  }
}

void FmuInternal::get_adj(FmuMemory* m, size_t ind, double* v) const {
  // Quick return if not needed
  if (!v) return;
  // Retrieve all sensitivities FIXME(@jaeandersson): should use compatible types
  for (size_t id : ired_[ind]) {
    casadi_int id2 = id;
    get_adj(m, 1, &id2, v++);
  }
}

int FmuInternal::discrete_states_iter(void* instance) const {
  // Quick return if no event indicators
  if (number_of_event_indicators_ == 0) return 0;
  // Helper function: update_discrete_states
  EventMemory eventmem;
  const size_t max_update_iter = 10;
  for (size_t update_iter = 0; update_iter < max_update_iter; ++update_iter) {
    if (update_discrete_states(instance, &eventmem)) {
      casadi_warning("update_discrete_states");
      return 1;
    }
    // Not implemented
    if (eventmem.discrete_states_need_update)
      casadi_warning("Discrete state update not implemented");
    if (eventmem.terminate_simulation)
      casadi_warning("Terminate solution not implemented");
    if (eventmem.nominals_of_continuous_states_changed)
      casadi_warning("Update of nominals of states not implemented");
    if (eventmem.values_of_continuous_states_changed)
      casadi_warning("Update of values of continuous states not implemented");
    if (eventmem.next_event_time_defined)
      casadi_warning("Ignoring next time defined: " + std::to_string(eventmem.next_event_time));
    // Successful return
    if (!eventmem.discrete_states_need_update) {
      casadi_warning("Discrete state update successful");
      return 0;
    }
  }
  // Too many iterations
  casadi_warning("Discrete state update failed");
  return 1;
}

int FmuInternal::init_mem(FmuMemory* m) const {
  // Ensure not already instantiated
  casadi_assert(m->instance == 0, "Already instantiated");
  // Create instance
  m->instance = instantiate();
  // Set all values
  if (set_values(m->instance)) {
    casadi_warning("FmuInternal::set_values failed");
    return 1;
  }
  // Initialization mode begins
  if (enter_initialization_mode(m->instance)) return 1;
  // Initialization mode ends
  if (exit_initialization_mode(m->instance)) return 1;
  // Iteration to update discrete states
  if (discrete_states_iter(m->instance)) return 1;
  // Allocate/reset input buffer
  m->ibuf_.resize(iind_.size());
  std::fill(m->ibuf_.begin(), m->ibuf_.end(), casadi::nan);
  // Allocate/reset output buffer
  m->obuf_.resize(oind_.size());
  std::fill(m->obuf_.begin(), m->obuf_.end(), casadi::nan);
  // Maximum input or output
  size_t max_io = std::max(iind_.size(), oind_.size());
  // Allocate/reset seeds
  m->isens_.resize(max_io);
  std::fill(m->isens_.begin(), m->isens_.end(), 0);
  // Allocate/reset sensitivities
  m->osens_.resize(max_io);
  std::fill(m->osens_.begin(), m->osens_.end(), 0);
  // Allocate/reset changed
  m->imarked_.resize(max_io);
  std::fill(m->imarked_.begin(), m->imarked_.end(), false);
  // Allocate/reset requested
  m->omarked_.resize(max_io);
  std::fill(m->omarked_.begin(), m->omarked_.end(), false);
  // Also allocate memory for corresponding Jacobian entry (for debugging)
  m->wrt_.resize(max_io);
  // Successful return
  return 0;
}

void FmuInternal::set(FmuMemory* m, size_t ind, const double* value) const {
  if (value) {
    // Argument is given
    for (size_t id : ired_[ind]) {
      if (*value != m->ibuf_.at(id)) {
        m->ibuf_.at(id) = *value;
        m->imarked_.at(id) = true;
      }
      value++;
    }
  } else {
    // Argument is null - all zeros
    for (size_t id : ired_[ind]) {
      if (0 != m->ibuf_.at(id)) {
        m->ibuf_.at(id) = 0;
        m->imarked_.at(id) = true;
      }
    }
  }
}

void FmuInternal::request(FmuMemory* m, size_t ind) const {
  for (size_t id : ored_[ind]) {
    // Mark as requested
    m->omarked_.at(id) = true;
    // Also log corresponding input index
    m->wrt_.at(id) = -1;
  }
}

int FmuInternal::eval(FmuMemory* m) const {
  // Gather inputs and outputs
  gather_io(m);
  // Number of inputs and outputs
  size_t n_set = m->id_in_.size();
  size_t n_out = m->id_out_.size();
  // Set all variables
  if (set_real(m->instance, get_ptr(m->vr_in_), n_set, get_ptr(m->v_in_), n_set)) {
    casadi_warning("Setting FMU variables failed");
    return 1;
  }
  // Quick return if nothing requested
  if (n_out == 0) return 0;
  // Calculate all variables
  m->v_out_.resize(n_out);
  if (get_real(m->instance, get_ptr(m->vr_out_), n_out, get_ptr(m->v_out_), n_out)) {
    casadi_warning("Evaluation failed");
    return 1;
  }
  // Collect requested variables
  auto it = m->v_out_.begin();
  for (size_t id : m->id_out_) {
    m->obuf_[id] = *it++;
  }
  // Successful return
  return 0;
}

void FmuInternal::get(FmuMemory* m, size_t ind, double* value) const {
  // Save to return
  for (size_t id : ored_[ind]) {
    *value++ = m->obuf_.at(id);
  }
}

void FmuInternal::set_fwd(FmuMemory* m, casadi_int nseed,
    const casadi_int* id, const double* v) const {
  for (casadi_int i = 0; i < nseed; ++i) {
    m->isens_.at(*id) = *v++;
    m->imarked_.at(*id) = true;
    id++;
  }
}

void FmuInternal::set_fwd(FmuMemory* m, size_t ind, const double* v) const {
  // If seeds are zero, no need to add to seed buffers
  if (!v) return;
  // Pass all seeds FIXME(@jaeandersson): should use compatible types
  for (size_t id : ired_[ind]) {
    casadi_int id2 = id;
    set_fwd(m, 1, &id2, v++);
  }
}

void FmuInternal::request_fwd(FmuMemory* m, casadi_int nsens, const casadi_int* id,
    const casadi_int* wrt_id) const {
  for (casadi_int i = 0; i < nsens; ++i) {
    m->omarked_.at(*id) = true;
    m->wrt_.at(*id) = *wrt_id++;
    id++;
  }
}

void FmuInternal::request_fwd(FmuMemory* m, casadi_int ind) const {
  // Request all sensitivities FIXME(@jaeandersson): should use compatible types
  casadi_int wrt_id = -1;
  for (size_t id : ored_[ind]) {
    casadi_int id2 = id;
    request_fwd(m, 1, &id2, &wrt_id);
  }
}

void FmuInternal::gather_io(FmuMemory* m) const {
  // Collect input indices and corresponding value references and values
  m->id_in_.clear();
  m->vr_in_.clear();
  m->v_in_.clear();
  for (size_t id = 0; id < m->imarked_.size(); ++id) {
    if (m->imarked_[id]) {
      m->id_in_.push_back(id);
      m->vr_in_.push_back(vr_in_[id]);
      m->v_in_.push_back(m->ibuf_[id]);
      m->imarked_[id] = false;
    }
  }
  // Collect output indices, corresponding value references
  m->id_out_.clear();
  m->vr_out_.clear();
  for (size_t id = 0; id < m->omarked_.size(); ++id) {
    if (m->omarked_[id]) {
      m->id_out_.push_back(id);
      m->vr_out_.push_back(vr_out_[id]);
      m->omarked_[id] = false;
    }
  }
}

void FmuInternal::gather_fwd(FmuMemory* m) const {
  // Gather input and output indices
  gather_io(m);
  // Number of inputs and outputs
  size_t n_known = m->id_in_.size();
  size_t n_unknown = m->id_out_.size();
  // Get/clear seeds
  m->d_in_.clear();
  for (size_t id : m->id_in_) {
    m->d_in_.push_back(m->isens_[id]);
    m->isens_[id] = 0;
  }
  // Ensure at least one seed
  casadi_assert(n_known != 0, "No seeds");
  // Allocate result vectors
  m->v_out_.resize(n_unknown);
  m->d_out_.resize(n_unknown);
}

void FmuInternal::gather_adj(FmuMemory* m) const {
  // Gather input and output indices
  gather_io(m);
  // Number of inputs and outputs
  size_t n_known = m->id_in_.size();
  size_t n_unknown = m->id_out_.size();
  // Get/clear seeds
  m->d_out_.clear();
  for (size_t id : m->id_out_) {
    m->d_out_.push_back(m->osens_[id]);
    m->osens_[id] = 0;
  }
  // Ensure at least one seed
  casadi_assert(n_unknown != 0, "No seeds");
  // Allocate result vectors
  m->v_in_.resize(n_known);
  m->d_in_.resize(n_known);
}

void Fmu::serialize(SerializingStream &s) const {
  return (*this)->serialize(s);
}

Fmu Fmu::deserialize(DeserializingStream& s) {
  return Fmu::create(FmuInternal::deserialize(s));
}

Fmu Fmu::create(FmuInternal* node) {
  Fmu ret;
  ret.own(node);
  return ret;
}

void FmuInternal::serialize(SerializingStream& s) const {
  serialize_type(s);
  serialize_body(s);
}

void FmuInternal::serialize_type(SerializingStream& s) const {
  s.pack("FmuInternal::type", class_name());
}

void FmuInternal::serialize_body(SerializingStream& s) const {
  s.version("FmuInternal", 2);
  s.pack("FmuInternal::name", name_);
  s.pack("FmuInternal::scheme_in", scheme_in_);
  s.pack("FmuInternal::scheme_out", scheme_out_);
  s.pack("FmuInternal::scheme", scheme_);
  s.pack("FmuInternal::aux", aux_);
  s.pack("FmuInternal::li", li_);
  s.pack("FmuInternal::iind", iind_);
  s.pack("FmuInternal::iind_map", iind_map_);
  s.pack("FmuInternal::oind", oind_);
  s.pack("FmuInternal::oind_map", oind_map_);
  s.pack("FmuInternal::nominal_in", nominal_in_);
  s.pack("FmuInternal::nominal_out", nominal_out_);
  s.pack("FmuInternal::min_in", min_in_);
  s.pack("FmuInternal::min_out", min_out_);
  s.pack("FmuInternal::max_in", max_in_);
  s.pack("FmuInternal::max_out", max_out_);
  s.pack("FmuInternal::vn_in", vn_in_);
  s.pack("FmuInternal::vn_out", vn_out_);
  s.pack("FmuInternal::vr_in", vr_in_);
  s.pack("FmuInternal::vr_out", vr_out_);

  s.pack("FmuInternal::value_in", value_in_);
  s.pack("FmuInternal::ired", ired_);
  s.pack("FmuInternal::ored", ored_);
  s.pack("FmuInternal::jac_sp", jac_sp_);
  s.pack("FmuInternal::hess_sp", hess_sp_);

  s.pack("FmuInternal::resource_loc", resource_loc_);
  s.pack("FmuInternal::fmutol", fmutol_);
  s.pack("FmuInternal::instance_name", instance_name_);
  s.pack("FmuInternal::instantiation_token", instantiation_token_);
  s.pack("FmuInternal::logging_on", logging_on_);
  s.pack("FmuInternal::number_of_event_indicators", number_of_event_indicators_);
  s.pack("FmuInternal::provides_directional_derivatives", provides_directional_derivatives_);
  s.pack("FmuInternal::provides_adjoint_derivatives", provides_adjoint_derivatives_);
}

FmuInternal::FmuInternal(DeserializingStream& s) {
  s.version("FmuInternal", 2);
  s.unpack("FmuInternal::name", name_);
  s.unpack("FmuInternal::scheme_in", scheme_in_);
  s.unpack("FmuInternal::scheme_out", scheme_out_);
  s.unpack("FmuInternal::scheme", scheme_);
  s.unpack("FmuInternal::aux", aux_);
  s.unpack("FmuInternal::li", li_);
  s.unpack("FmuInternal::iind", iind_);
  s.unpack("FmuInternal::iind_map", iind_map_);
  s.unpack("FmuInternal::oind", oind_);
  s.unpack("FmuInternal::oind_map", oind_map_);
  s.unpack("FmuInternal::nominal_in", nominal_in_);
  s.unpack("FmuInternal::nominal_out", nominal_out_);
  s.unpack("FmuInternal::min_in", min_in_);
  s.unpack("FmuInternal::min_out", min_out_);
  s.unpack("FmuInternal::max_in", max_in_);
  s.unpack("FmuInternal::max_out", max_out_);
  s.unpack("FmuInternal::vn_in", vn_in_);
  s.unpack("FmuInternal::vn_out", vn_out_);
  s.unpack("FmuInternal::vr_in", vr_in_);
  s.unpack("FmuInternal::vr_out", vr_out_);

  s.unpack("FmuInternal::value_in", value_in_);
  s.unpack("FmuInternal::ired", ired_);
  s.unpack("FmuInternal::ored", ored_);
  s.unpack("FmuInternal::jac_sp", jac_sp_);
  s.unpack("FmuInternal::hess_sp", hess_sp_);

  s.unpack("FmuInternal::resource_loc", resource_loc_);
  s.unpack("FmuInternal::fmutol", fmutol_);
  s.unpack("FmuInternal::instance_name", instance_name_);
  s.unpack("FmuInternal::instantiation_token", instantiation_token_);
  s.unpack("FmuInternal::logging_on", logging_on_);
  s.unpack("FmuInternal::number_of_event_indicators", number_of_event_indicators_);
  s.unpack("FmuInternal::provides_directional_derivatives", provides_directional_derivatives_);
  s.unpack("FmuInternal::provides_adjoint_derivatives", provides_adjoint_derivatives_);
}

FmuInternal* FmuInternal::deserialize(DeserializingStream& s) {
  std::string class_name;
  s.unpack("FmuInternal::type", class_name);
  if (class_name=="Fmu2") {
#ifdef WITH_FMI2
    return Fmu2::deserialize(s);
#else
    casadi_error("CasADi was not compiled with WITH_FMI2=ON.");
#endif // WITH_FMI2
  } else if (class_name=="Fmu3") {
#ifdef WITH_FMI3
    casadi_error("Not implemented");
#else
    casadi_error("CasADi was not compiled with WITH_FMI2=ON.");
#endif // WITH_FMI3
  } else {
    casadi_error("Cannot deserialize type '" + class_name + "'");
  }
}

} // namespace casadi

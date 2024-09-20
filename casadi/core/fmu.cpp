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

bool Fmu::has_fwd() const {
  try {
    return (*this)->has_fwd();
  } catch(std::exception& e) {
    THROW_ERROR("has_fwd", e.what());
  }
}

bool Fmu::has_adj() const {
  try {
    return (*this)->has_adj();
  } catch(std::exception& e) {
    THROW_ERROR("has_adj", e.what());
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

void Fmu::set_seed(FmuMemory* m, casadi_int nseed, const casadi_int* id, const double* v) const {
  try {
    return (*this)->set_seed(m, nseed, id, v);
  } catch(std::exception& e) {
    THROW_ERROR("set_seed", e.what());
  }
}

void Fmu::request_sens(FmuMemory* m, casadi_int nsens, const casadi_int* id,
    const casadi_int* wrt_id) const {
  try {
    return (*this)->request_sens(m, nsens, id, wrt_id);
  } catch(std::exception& e) {
    THROW_ERROR("request_sens", e.what());
  }
}

int Fmu::eval_derivative(FmuMemory* m, bool independent_seeds) const {
  try {
    return (*this)->eval_derivative(m, independent_seeds);
  } catch(std::exception& e) {
    THROW_ERROR("eval_derivative", e.what());
  }
}

void Fmu::get_sens(FmuMemory* m, casadi_int nsens, const casadi_int* id, double* v) const {
  try {
    return (*this)->get_sens(m, nsens, id, v);
  } catch(std::exception& e) {
    THROW_ERROR("get_sens", e.what());
  }
}

void Fmu::set_fwd(FmuMemory* m, size_t ind, const double* v) const {
  try {
    return (*this)->set_fwd(m, ind, v);
  } catch(std::exception& e) {
    THROW_ERROR("set_fwd", e.what());
  }
}

void Fmu::request_fwd(FmuMemory* m, casadi_int ind) const {
  try {
    return (*this)->request_fwd(m, ind);
  } catch(std::exception& e) {
    THROW_ERROR("request_fwd", e.what());
  }
}

void Fmu::get_fwd(FmuMemory* m, size_t ind, double* v) const {
  try {
    return (*this)->get_fwd(m, ind, v);
  } catch(std::exception& e) {
    THROW_ERROR("get_fwd", e.what());
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

int FmuInternal::eval_derivative(FmuMemory* m, bool independent_seeds) const {
  // Gather input and output indices
  gather_sens(m);
  // Calculate derivatives using FMU directional derivative support
  if (m->self.enable_ad_) {
    // Evaluate using AD
    if (eval_ad(m)) return 1;
  }
  // Calculate derivatives using finite differences
  if (!m->self.enable_ad_ || m->self.validate_ad_) {
    // Evaluate using FD
    if (eval_fd(m, independent_seeds)) return 1;
  }
  return 0;
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
  // Allocate/reset input buffer
  m->ibuf_.resize(iind_.size());
  std::fill(m->ibuf_.begin(), m->ibuf_.end(), casadi::nan);
  // Allocate/reset output buffer
  m->obuf_.resize(oind_.size());
  std::fill(m->obuf_.begin(), m->obuf_.end(), casadi::nan);
  // Allocate/reset seeds
  m->seed_.resize(iind_.size());
  std::fill(m->seed_.begin(), m->seed_.end(), 0);
  // Allocate/reset sensitivities
  m->sens_.resize(oind_.size());
  std::fill(m->sens_.begin(), m->sens_.end(), 0);
  // Allocate/reset changed
  m->changed_.resize(iind_.size());
  std::fill(m->changed_.begin(), m->changed_.end(), false);
  // Allocate/reset requested
  m->requested_.resize(oind_.size());
  std::fill(m->requested_.begin(), m->requested_.end(), false);
  // Also allocate memory for corresponding Jacobian entry (for debugging)
  m->wrt_.resize(oind_.size());
  // Successful return
  return 0;
}

void FmuInternal::set(FmuMemory* m, size_t ind, const double* value) const {
  if (value) {
    // Argument is given
    for (size_t id : ired_[ind]) {
      if (*value != m->ibuf_.at(id)) {
        m->ibuf_.at(id) = *value;
        m->changed_.at(id) = true;
      }
      value++;
    }
  } else {
    // Argument is null - all zeros
    for (size_t id : ired_[ind]) {
      if (0 != m->ibuf_.at(id)) {
        m->ibuf_.at(id) = 0;
        m->changed_.at(id) = true;
      }
    }
  }
}

void FmuInternal::request(FmuMemory* m, size_t ind) const {
  for (size_t id : ored_[ind]) {
    // Mark as requested
    m->requested_.at(id) = true;
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

void FmuInternal::set_seed(FmuMemory* m, casadi_int nseed,
    const casadi_int* id, const double* v) const {
  for (casadi_int i = 0; i < nseed; ++i) {
    m->seed_.at(*id) = *v++;
    m->changed_.at(*id) = true;
    id++;
  }
}

void FmuInternal::request_sens(FmuMemory* m, casadi_int nsens, const casadi_int* id,
    const casadi_int* wrt_id) const {
  for (casadi_int i = 0; i < nsens; ++i) {
    m->requested_.at(*id) = true;
    m->wrt_.at(*id) = *wrt_id++;
    id++;
  }
}

void FmuInternal::get_sens(FmuMemory* m, casadi_int nsens, const casadi_int* id, double* v) const {
  for (casadi_int i = 0; i < nsens; ++i) {
    *v++ = m->sens_.at(*id++);
  }
}

void FmuInternal::set_fwd(FmuMemory* m, size_t ind, const double* v) const {
  // If seeds are zero, no need to add to seed buffers
  if (!v) return;
  // Pass all seeds FIXME(@jaeandersson): should use compatible types
  for (size_t id : ired_[ind]) {
    casadi_int id2 = id;
    set_seed(m, 1, &id2, v++);
  }
}

void FmuInternal::request_fwd(FmuMemory* m, casadi_int ind) const {
  // Request all sensitivities FIXME(@jaeandersson): should use compatible types
  casadi_int wrt_id = -1;
  for (size_t id : ored_[ind]) {
    casadi_int id2 = id;
    request_sens(m, 1, &id2, &wrt_id);
  }
}

void FmuInternal::get_fwd(FmuMemory* m, size_t ind, double* v) const {
  // Quick return if not needed
  if (!v) return;
  // Retrieve all sensitivities FIXME(@jaeandersson): should use compatible types
  for (size_t id : ored_[ind]) {
    casadi_int id2 = id;
    get_sens(m, 1, &id2, v++);
  }
}

void FmuInternal::gather_io(FmuMemory* m) const {
  // Collect input indices and corresponding value references and values
  m->id_in_.clear();
  m->vr_in_.clear();
  m->v_in_.clear();
  for (size_t id = 0; id < m->changed_.size(); ++id) {
    if (m->changed_[id]) {
      m->id_in_.push_back(id);
      m->vr_in_.push_back(vr_in_[id]);
      m->v_in_.push_back(m->ibuf_[id]);
      m->changed_[id] = false;
    }
  }
  // Collect output indices, corresponding value references
  m->id_out_.clear();
  m->vr_out_.clear();
  for (size_t id = 0; id < m->requested_.size(); ++id) {
    if (m->requested_[id]) {
      m->id_out_.push_back(id);
      m->vr_out_.push_back(vr_out_[id]);
      m->requested_[id] = false;
    }
  }
}

void FmuInternal::gather_sens(FmuMemory* m) const {
  // Gather input and output indices
  gather_io(m);
  // Number of inputs and outputs
  size_t n_known = m->id_in_.size();
  size_t n_unknown = m->id_out_.size();
  // Get/clear seeds
  m->d_in_.clear();
  for (size_t id : m->id_in_) {
    m->d_in_.push_back(m->seed_[id]);
    m->seed_[id] = 0;
  }
  // Ensure at least one seed
  casadi_assert(n_known != 0, "No seeds");
  // Allocate result vectors
  m->v_out_.resize(n_unknown);
  m->d_out_.resize(n_unknown);
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
  s.version("FmuInternal", 1);
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


}

FmuInternal::FmuInternal(DeserializingStream& s) {
  s.version("FmuInternal", 1);
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

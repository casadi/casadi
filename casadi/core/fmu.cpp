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

namespace casadi {

// Throw informative error message
#define THROW_ERROR(FNAME, WHAT) \
throw CasadiException("Error in Fmu::" FNAME " for '" + this->name() + "' "\
  "[" + this->class_name() + "] at " + CASADI_WHERE + ":\n"\
  + std::string(WHAT));

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
  } else {
    // Not supported
    casadi_error("Unsupported FMU API: " + to_string(api));
  }
  // Initialize
  try {
    (*this)->init(dae);
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

bool Fmu::has_ad() const {
  try {
    return (*this)->has_ad();
  } catch(std::exception& e) {
    THROW_ERROR("has_ad", e.what());
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

void Fmu::free_instance(void* c) const {
  try {
    return (*this)->free_instance(c);
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

void FmuInternal::disp(std::ostream& stream, bool more) const {
  (void)more;  // unused
  stream << name_ << " " << class_name();
}

std::string to_string(FmuApi v) {
  switch (v) {
  case FmuApi::FMI2: return "fmi2";
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

} // namespace casadi

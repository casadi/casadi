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


#include "dae_builder_internal.hpp"

#include <cctype>
#include <ctime>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <algorithm>

#include "casadi_misc.hpp"
#include "exception.hpp"
#include "code_generator.hpp"
#include "calculus.hpp"
#include "xml_file.hpp"
#include "external.hpp"
#include "integrator.hpp"

namespace casadi {

// Throw informative error message
#define THROW_ERROR(FNAME, WHAT) \
throw CasadiException("Error in DaeBuilder::" FNAME " for '" + this->name() \
  + "'  at " + CASADI_WHERE + ":\n" + std::string(WHAT));

DaeBuilder::DaeBuilder() {
}

DaeBuilder::DaeBuilder(const std::string& name, const std::string& path, const Dict& opts) {
  own(new DaeBuilderInternal(name, path, opts));
  if (!path.empty()) load_fmi_description(path + "/modelDescription.xml");
}

const std::string& DaeBuilder::name() const {
  return (*this)->name_;
}

const MX& DaeBuilder::time() const {
  try {
    return (*this)->time();
  } catch (std::exception& e) {
    THROW_ERROR("time", e.what());
    static const MX t;
    return t;  // never reached
  }
}

std::vector<std::string> DaeBuilder::y() const {
  try {
    return (*this)->name((*this)->outputs_);
  } catch (std::exception& e) {
    THROW_ERROR("y", e.what());
    return {};  // never reached
  }
}

std::vector<MX> DaeBuilder::ode() const {
  try {
    return (*this)->output(OutputCategory::ODE);
  } catch (std::exception& e) {
    THROW_ERROR("ode", e.what());
    return {};  // never reached
  }
}

std::vector<MX> DaeBuilder::alg() const {
  try {
    return (*this)->output(OutputCategory::ALG);
  } catch (std::exception& e) {
    THROW_ERROR("alg", e.what());
    return {};  // never reached
  }
}

std::vector<MX> DaeBuilder::quad() const {
  try {
    return (*this)->output(OutputCategory::QUAD);
  } catch (std::exception& e) {
    THROW_ERROR("quad", e.what());
    return {};  // never reached
  }
}

std::vector<MX> DaeBuilder::zero() const {
  try {
    return (*this)->output(OutputCategory::ZERO);
  } catch (std::exception& e) {
    THROW_ERROR("zero", e.what());
  }
}

std::vector<MX> DaeBuilder::ydef() const {
  try {
    return (*this)->output(OutputCategory::Y);
  } catch (std::exception& e) {
    THROW_ERROR("ydef", e.what());
    return {};  // never reached
  }
}

std::vector<MX> DaeBuilder::cdef() const {
  try {
    return (*this)->cdef();
  } catch (std::exception& e) {
    THROW_ERROR("cdef", e.what());
    return {};  // never reached
  }
}

std::vector<MX> DaeBuilder::ddef() const {
  try {
    return (*this)->output(OutputCategory::D);
  } catch (std::exception& e) {
    THROW_ERROR("ddef", e.what());
    return {};  // never reached
  }
}

std::vector<MX> DaeBuilder::wdef() const {
  try {
    return (*this)->output(OutputCategory::W);
  } catch (std::exception& e) {
    THROW_ERROR("wdef", e.what());
    return {};  // never reached
  }
}

std::vector<MX> DaeBuilder::init_lhs() const {
  return (*this)->init_lhs();
}

std::vector<MX> DaeBuilder::init_rhs() const {
  return (*this)->init_rhs();
}

std::vector<std::string> DaeBuilder::outputs() const {
  try {
    return (*this)->name((*this)->outputs_);
  } catch (std::exception& e) {
    THROW_ERROR("outputs", e.what());
    return {};  // never reached
  }
}

std::vector<std::string> DaeBuilder::derivatives() const {
  try {
    return (*this)->name((*this)->derivatives_);
  } catch (std::exception& e) {
    THROW_ERROR("derivatives", e.what());
    return {};  // never reached
  }
}

std::vector<std::string> DaeBuilder::initial_unknowns() const {
  try {
    return (*this)->name((*this)->initial_unknowns_);
  } catch (std::exception& e) {
    THROW_ERROR("initial_unknowns", e.what());
    return {};  // never reached
  }
}

bool DaeBuilder::has_t() const {
  try {
    return (*this)->has_t();
  } catch (std::exception& e) {
    THROW_ERROR("has_t", e.what());
    return false;  // never reached
  }
}

casadi_int DaeBuilder::nx() const {
  return (*this)->size(Category::X);
}

casadi_int DaeBuilder::nz() const {
  return (*this)->size(Category::Z);
}

casadi_int DaeBuilder::nq() const {
  return (*this)->size(Category::Q);
}

casadi_int DaeBuilder::nzero() const {
  return (*this)->event_indicators_.size();
}

casadi_int DaeBuilder::ny() const {
  return (*this)->outputs_.size();
}

casadi_int DaeBuilder::nu() const {
  return (*this)->size(Category::U);
}

casadi_int DaeBuilder::np() const {
  return (*this)->size(Category::P);
}

casadi_int DaeBuilder::nc() const {
  return (*this)->size(Category::C);
}

casadi_int DaeBuilder::nd() const {
  return (*this)->size(Category::D);
}

casadi_int DaeBuilder::nw() const {
  return (*this)->size(Category::W);
}

void DaeBuilder::load_fmi_description(const std::string& filename) {
  try {
    (*this)->load_fmi_description(filename);
  } catch (std::exception& e) {
    THROW_ERROR("load_fmi_description", e.what());
  }
}

bool DaeBuilder::provides_directional_derivatives() const {
  try {
    casadi_assert(!(*this)->symbolic_, "Functionality only applies to imported standard FMUs");
    return (*this)->provides_directional_derivatives_;
  } catch (std::exception& e) {
    THROW_ERROR("provides_directional_derivatives", e.what());
    return false;
  }
}

std::vector<std::string> DaeBuilder::export_fmu(const Dict& opts) {
  try {
    return (*this)->export_fmu(opts);
  } catch (std::exception& e) {
    THROW_ERROR("export_fmu", e.what());
  }
}

void DaeBuilder::eliminate_quad() {
  try {
    (*this)->eliminate_quad();
  } catch (std::exception& e) {
    THROW_ERROR("eliminate_quad", e.what());
  }
}

void DaeBuilder::sort_d() {
  try {
    (*this)->sort_d();
  } catch (std::exception& e) {
    THROW_ERROR("sort_d", e.what());
  }
}

void DaeBuilder::sort_w() {
  try {
    (*this)->sort_w();
  } catch (std::exception& e) {
    THROW_ERROR("sort_w", e.what());
  }
}

void DaeBuilder::sort_z(const std::vector<std::string>& z_order) {
  try {
    (*this)->sort_z(z_order);
  } catch (std::exception& e) {
    THROW_ERROR("sort_z", e.what());
  }
}

void DaeBuilder::prune(bool prune_p, bool prune_u) {
  try {
    (*this)->prune(prune_p, prune_u);
  } catch (std::exception& e) {
    THROW_ERROR("prune", e.what());
  }
}

void DaeBuilder::tear() {
  try {
    (*this)->tear();
  } catch (std::exception& e) {
    THROW_ERROR("tear", e.what());
  }
}

bool DaeBuilder::has(const std::string& name) const {
  try {
    return (*this)->has(name);
  } catch (std::exception& e) {
    THROW_ERROR("has", e.what());
    return false;  // never reached
  }
}

std::vector<std::string> DaeBuilder::all() const {
  try {
    return (*this)->all();
  } catch (std::exception& e) {
    THROW_ERROR("all", e.what());
    return {};  // never reached
  }
}

std::vector<std::string> DaeBuilder::all(const std::string& cat) const {
  try {
    return (*this)->all(to_enum<Category>(cat));
  } catch (std::exception& e) {
    THROW_ERROR("all", e.what());
    return {};  // never reached
  }
}

#ifdef WITH_DEPRECATED_FEATURES
Variable& DaeBuilder::new_variable(const std::string& name, casadi_int numel) {
  try {
    return (*this)->new_variable(name, {numel});
  } catch (std::exception& e) {
    THROW_ERROR("new_variable", e.what());
  }
}

Variable& DaeBuilder::variable(const std::string& name) {
  try {
    return (*this)->variable(name);
  } catch (std::exception& e) {
    THROW_ERROR("variable", e.what());
  }
}

const Variable& DaeBuilder::variable(const std::string& name) const {
  try {
    return (*this)->variable(name);
  } catch (std::exception& e) {
    THROW_ERROR("variable", e.what());
  }
}

Variable& DaeBuilder::variable(size_t ind) {
  try {
    return (*this)->variable(ind);
  } catch (std::exception& e) {
    THROW_ERROR("variable", e.what());
  }
}

const Variable& DaeBuilder::variable(size_t ind) const {
  try {
    return (*this)->variable(ind);
  } catch (std::exception& e) {
    THROW_ERROR("variable", e.what());
  }
}

size_t DaeBuilder::find(const std::string& name) const {
  try {
    return (*this)->find(name);
  } catch (std::exception& e) {
    THROW_ERROR("find", e.what());
    return -1; // never reached
  }
}

std::vector<size_t> DaeBuilder::find(const std::vector<std::string>& name) const {
  try {
    return (*this)->find(name);
  } catch (std::exception& e) {
    THROW_ERROR("find", e.what());
    return {}; // never reached
  }
}

const std::string& DaeBuilder::name(size_t ind) const {
  try {
    return (*this)->name(ind);
  } catch (std::exception& e) {
    THROW_ERROR("name", e.what());
    static std::string dummy;
    return dummy; // never reached
  }
}

std::vector<std::string> DaeBuilder::name(const std::vector<size_t>& ind) const {
  try {
    return (*this)->name(ind);
  } catch (std::exception& e) {
    THROW_ERROR("name", e.what());
    return {}; // never reached
  }
}
const MX& DaeBuilder::var(size_t ind) const {
  try {
    return (*this)->var(ind);
  } catch (std::exception& e) {
    THROW_ERROR("var", e.what());
    static MX dummy;
    return dummy; // never reached
  }
}

std::vector<MX> DaeBuilder::var(const std::vector<size_t>& ind) const {
  try {
    return (*this)->var(ind);
  } catch (std::exception& e) {
    THROW_ERROR("var", e.what());
    return {}; // never reached
  }
}

MX DaeBuilder::add_variable(const std::string& name, casadi_int n) {
  return add_variable(name, Sparsity::dense(n));
}

MX DaeBuilder::add_variable(const std::string& name, const Sparsity& sp) {
  Variable& v = new_variable(name);
  v.v = MX::sym(name, sp);
  return v.v;
}

void DaeBuilder::add_variable(const MX& new_v) {
  Variable& v = new_variable(new_v.name());
  v.v = new_v;
}

size_t DaeBuilder::add_variable_new(const std::string& name, casadi_int n) {
  return add_variable_new(name, Sparsity::dense(n));
}

size_t DaeBuilder::add_variable_new(const std::string& name, const Sparsity& sp) {
  Variable& v = new_variable(name);
  v.v = MX::sym(name, sp);
  return v.index;
}

size_t DaeBuilder::add_variable_new(const MX& new_v) {
  Variable& v = new_variable(new_v.name());
  v.v = new_v;
  return v.index;
}

void DaeBuilder::register_t(const std::string& name) {
  // Save to class
  casadi_assert(!has_t(), "'t' already defined");
  (*this)->indices(Category::T).push_back(find(name));
}

void DaeBuilder::register_p(const std::string& name) {
  (*this)->indices(Category::P).push_back(find(name));
}

void DaeBuilder::register_u(const std::string& name) {
  (*this)->indices(Category::U).push_back(find(name));
}

void DaeBuilder::register_x(const std::string& name) {
  (*this)->indices(Category::X).push_back(find(name));
}

void DaeBuilder::register_z(const std::string& name) {
  (*this)->indices(Category::Z).push_back(find(name));
}

void DaeBuilder::register_q(const std::string& name) {
  (*this)->indices(Category::Q).push_back(find(name));
}

void DaeBuilder::register_c(const std::string& name) {
  (*this)->indices(Category::C).push_back(find(name));
}

void DaeBuilder::register_d(const std::string& name) {
  (*this)->indices(Category::D).push_back(find(name));
}

void DaeBuilder::register_w(const std::string& name) {
  (*this)->indices(Category::W).push_back(find(name));
}

void DaeBuilder::register_y(const std::string& name) {
  (*this)->outputs_.push_back(find(name));
}

void DaeBuilder::register_e(const std::string& name) {
  (*this)->event_indicators_.push_back(find(name));
}

void DaeBuilder::clear_all(const std::string& v) {
  try {
    (*this)->clear_cache_ = true;  // Clear cache after this
    (*this)->indices(to_enum<Category>(v)).clear();
  } catch (std::exception& e) {
    THROW_ERROR("clear_all", e.what());
  }
}

void DaeBuilder::set_all(const std::string& v, const std::vector<std::string>& name) {
  try {
    (*this)->clear_cache_ = true;  // Clear cache after this
    (*this)->indices(to_enum<Category>(v)) = (*this)->find(name);
  } catch (std::exception& e) {
    THROW_ERROR("set_all", e.what());
  }
}

#endif // WITH_DEPRECATED_FEATURES

void DaeBuilder::reorder(const std::string& cat, const std::vector<std::string>& v) {
  try {
    auto vind = (*this)->find(v);
    if (cat == "y") {
      // Reorder outputs
      (*this)->reorder("y", (*this)->outputs_, vind);
    } else {
      // Reorder inputs
      (*this)->reorder(to_enum<Category>(cat), vind);
    }
  } catch (std::exception& e) {
    THROW_ERROR("reorder", e.what());
  }
}

MX DaeBuilder::add(const std::string& name, const std::string& causality,
    const std::string& variability, const Dict& opts) {
  try {
    return (*this)->add(name, to_enum<Causality>(causality),
      to_enum<Variability>(variability), opts).v;
  } catch (std::exception& e) {
    THROW_ERROR("add", e.what());
    return MX();
  }
}

MX DaeBuilder::add(const std::string& name, const std::string& causality, const Dict& opts) {
  try {
    return (*this)->add(name, to_enum<Causality>(causality), opts).v;
  } catch (std::exception& e) {
    THROW_ERROR("add", e.what());
    return MX();
  }
}

MX DaeBuilder::add(const std::string& name, const Dict& opts) {
  try {
    return (*this)->add(name, opts).v;
  } catch (std::exception& e) {
    THROW_ERROR("add", e.what());
    return MX();
  }
}

void DaeBuilder::add(const std::string& name, const std::string& causality,
  const std::string& variability, const MX& expr, const Dict& opts) {
  try {
    // Ensure pure symbolic expression
    casadi_assert(expr.is_symbolic(), "Expression must be symbolic");
    // Make sure name matches expression
    casadi_assert(name == expr.name(), "Name must match expression");
    // Add variable
    (*this)->add(name, to_enum<Causality>(causality), to_enum<Variability>(variability),
      expr, opts);
  } catch (std::exception& e) {
    THROW_ERROR("add", e.what());
  }
}

#ifdef WITH_DEPRECATED_FEATURES
MX DaeBuilder::add_t(const std::string& name) {
  return add(name, "independent");
}

MX DaeBuilder::add_p(const std::string& name) {
  casadi_assert(!name.empty(), "Variable name is required");
  return add(name, "parameter", "tunable");
}

MX DaeBuilder::add_u(const std::string& name) {
  casadi_assert(!name.empty(), "Variable name is required");
  return add(name, "input");
}

MX DaeBuilder::add_x(const std::string& name) {
  casadi_assert(!name.empty(), "Variable name is required");
  return add(name);
}

MX DaeBuilder::add_z(const std::string& name) {
  casadi_assert(!name.empty(), "Variable name is required");
  return add(name);
}

MX DaeBuilder::add_q(const std::string& name) {
  casadi_assert(!name.empty(), "Variable name is required");
  return add(name);
}

MX DaeBuilder::add_c(const std::string& name, const MX& new_cdef) {
  MX v = add(name, "local", "constant");
  set_beq(name, new_cdef);
  return v;
}

MX DaeBuilder::add_d(const std::string& name, const MX& new_ddef) {
  MX v = add(name, "calculatedParameter", "fixed");
  set_beq(name, new_ddef);
  return v;
}

MX DaeBuilder::add_w(const std::string& name, const MX& new_wdef) {
  MX v = add(name);
  eq(v, new_wdef);
  return v;
}

MX DaeBuilder::add_y(const std::string& name, const MX& new_ydef) {
  MX v = add(name, "output");
  eq(v, new_ydef);
  return v;
}

void DaeBuilder::set_beq(const std::string& name, const MX& val) {
  try {
    eq(var(name), val);
  } catch (std::exception& e) {
    THROW_ERROR("set_beq", e.what());
  }
}

#endif  // WITH_DEPRECATED_FEATURES

void DaeBuilder::eq(const MX& lhs, const MX& rhs, const Dict& opts) {
  try {
    (*this)->eq(lhs, rhs, opts);
  } catch (std::exception& e) {
    THROW_ERROR("eq", e.what());
  }
}

void DaeBuilder::when(const MX& cond, const std::vector<std::string>& eqs, const Dict& opts) {
  try {
    (*this)->when(cond, eqs, opts);
  } catch (std::exception& e) {
    THROW_ERROR("when", e.what());
  }
}

std::string DaeBuilder::assign(const std::string& name, const MX& val) {
  try {
    return (*this)->assign(name, val).name;
  } catch (std::exception& e) {
    THROW_ERROR("assign", e.what());
    return std::string();  // never reached
  }
}

std::string DaeBuilder::reinit(const std::string& name, const MX& val) {
  try {
    return (*this)->reinit(name, val).name;
  } catch (std::exception& e) {
    THROW_ERROR("reinit", e.what());
    return std::string();  // never reached
  }
}

void DaeBuilder::set_init(const std::string& name, const MX& init_rhs) {
  try {
    (*this)->set_init(name, init_rhs);
  } catch (std::exception& e) {
    THROW_ERROR("set_init", e.what());
  }
}

void DaeBuilder::sanity_check() const {
  try {
    (*this)->sanity_check();
  } catch (std::exception& e) {
    THROW_ERROR("sanity_check", e.what());
  }
}

MX DaeBuilder::var(const std::string& name) const {
  try {
    return (*this)->variable(name).v;
  } catch (std::exception& e) {
    THROW_ERROR("var", e.what());
    return MX();  // never reached
  }
}

std::string DaeBuilder::der(const std::string& name) const {
  try {
    // Get variable index
    size_t ind = (*this)->find(name);
    // Differentiate
    ind = (*this)->variable(ind).der;
    casadi_assert(ind != size_t(-1), "No derivative expression for " + name);
    // Return name
    return (*this)->variable(ind).name;
  } catch (std::exception& e) {
    THROW_ERROR("der", e.what());
    return std::string();  // never reached
  }
}

std::string DaeBuilder::pre(const std::string& name) const {
  try {
    // Not implemented
    static bool warned = false;
    if (!warned) {
      casadi_warning("DaeBuilder::pre has not been implemented: Returning identity mapping");
      warned = true;
    }
    return name;
  } catch (std::exception& e) {
    THROW_ERROR("pre", e.what());
    return std::string();  // never reached
  }
}

MX DaeBuilder::pre(const MX& v) const {
  try {
    // Not implemented
    static bool warned = false;
    if (!warned) {
      casadi_warning("DaeBuilder::pre has not been implemented: Returning identity mapping");
      warned = true;
    }
    return v;
  } catch (std::exception& e) {
    THROW_ERROR("pre", e.what());
    return MX();  // never reached
  }
}

std::vector<std::string> DaeBuilder::der(const std::vector<std::string>& name) const {
  try {
    std::vector<std::string> r(name.size());
    for (size_t i = 0; i < r.size(); ++i) r[i] = der(name[i]);
    return r;
  } catch (std::exception& e) {
    THROW_ERROR("der", e.what());
    return {};  // never reached
  }
}

std::vector<std::string> DaeBuilder::pre(const std::vector<std::string>& name) const {
  try {
    std::vector<std::string> r(name.size());
    for (size_t i = 0; i < r.size(); ++i) r[i] = pre(name[i]);
    return r;
  } catch (std::exception& e) {
    THROW_ERROR("pre", e.what());
    return {};  // never reached
  }
}

bool DaeBuilder::has_beq(const std::string& name) const {
  try {
    return !(*this)->variable(name).has_beq();
  } catch (std::exception& e) {
    THROW_ERROR("has_beq", e.what());
    return false;  // never reached
  }
}

MX DaeBuilder::beq(const std::string& name) const {
  try {
    const Variable& v = (*this)->variable(name);
    return (*this)->variable(v.bind).v;
  } catch (std::exception& e) {
    THROW_ERROR("beq", e.what());
    return MX();  // never reached
  }
}

void DaeBuilder::eliminate_d() {
  try {
    (*this)->eliminate_d();
  } catch (std::exception& e) {
    THROW_ERROR("eliminate_d", e.what());
  }
}
void DaeBuilder::eliminate_w() {
  try {
    (*this)->eliminate_w();
  } catch (std::exception& e) {
    THROW_ERROR("eliminate_w", e.what());
  }
}

void DaeBuilder::lift(bool lift_shared, bool lift_calls) {
  try {
    (*this)->lift(lift_shared, lift_calls);
  } catch (std::exception& e) {
    THROW_ERROR("lift", e.what());
  }
}

casadi_int DaeBuilder::value_reference(const std::string& name) const {
  return (*this)->variable(name).value_reference;
}

void DaeBuilder::set_value_reference(const std::string& name, casadi_int val) {
  (*this)->variable(name).value_reference = val;
}

std::string DaeBuilder::description(const std::string& name) const {
  return (*this)->variable(name).description;
}

void DaeBuilder::set_description(const std::string& name, const std::string& val) {
  (*this)->variable(name).description = val;
}

std::string DaeBuilder::type(const std::string& name, casadi_int fmi_version) const {
  // Check version
  casadi_assert(fmi_version == 2 || fmi_version == 3, "Only FMI version 2 or 3 supported");
  // Handle FMI 2
  if (fmi_version == 2) {
    return to_string(to_fmi2((*this)->variable(name).type));
  }
  // Assume FMI 3
  return to_string((*this)->variable(name).type);
}

void DaeBuilder::set_type(const std::string& name, const std::string& val) {
  // Fallback to FMI 2, if necessary
  if (has_enum<TypeFmi2>(val) && !has_enum<Type>(val)) {
    (*this)->variable(name).type = from_fmi2(to_enum<TypeFmi2>(val));
  }
  // Assume FMI 3
  (*this)->variable(name).type = to_enum<Type>(val);
}

std::string DaeBuilder::causality(const std::string& name) const {
  try {
    return to_string((*this)->causality((*this)->find(name)));
  } catch (std::exception& e) {
    THROW_ERROR("causality", e.what());
    return std::string();  // never reached
  }
}

void DaeBuilder::set_causality(const std::string& name, const std::string& val) {
  try {
    (*this)->set_causality((*this)->find(name), to_enum<Causality>(val));
  } catch (std::exception& e) {
    THROW_ERROR("set_causality", e.what());
  }
}

std::string DaeBuilder::variability(const std::string& name) const {
  try {
    return to_string((*this)->variability((*this)->find(name)));
  } catch (std::exception& e) {
    THROW_ERROR("variability", e.what());
    return std::string();  // never reached
  }
}

void DaeBuilder::set_variability(const std::string& name, const std::string& val) {
  try {
    (*this)->set_variability((*this)->find(name), to_enum<Variability>(val));
  } catch (std::exception& e) {
    THROW_ERROR("set_variability", e.what());
  }
}

std::string DaeBuilder::category(const std::string& name) const {
  try {
    return to_string((*this)->category((*this)->find(name)));
  } catch (std::exception& e) {
    THROW_ERROR("category", e.what());
    return std::string();  // never reached
  }
}

void DaeBuilder::set_category(const std::string& name, const std::string& val) {
  try {
    (*this)->set_category((*this)->find(name), to_enum<Category>(val));
  } catch (std::exception& e) {
    THROW_ERROR("set_category", e.what());
  }
}

std::string DaeBuilder::initial(const std::string& name) const {
  return to_string((*this)->variable(name).initial);
}

void DaeBuilder::set_initial(const std::string& name, const std::string& val) {
  (*this)->variable(name).initial = to_enum<Initial>(val);
}

std::string DaeBuilder::unit(const std::string& name) const {
  return (*this)->variable(name).unit;
}

void DaeBuilder::set_unit(const std::string& name, const std::string& val) {
  (*this)->variable(name).unit = val;
}

std::string DaeBuilder::display_unit(const std::string& name) const {
  return (*this)->variable(name).display_unit;
}

void DaeBuilder::set_display_unit(const std::string& name, const std::string& val) {
  (*this)->variable(name).display_unit = val;
}

casadi_int DaeBuilder::numel(const std::string& name) const {
  return (*this)->variable(name).numel;
}

std::vector<casadi_int> DaeBuilder::dimension(const std::string& name) const {
  return (*this)->variable(name).dimension;
}

double DaeBuilder::start_time() const {
  try {
    return (*this)->start_time_;
  } catch (std::exception& e) {
    THROW_ERROR("start_time", e.what());
    return nan;
  }
}

void DaeBuilder::set_start_time(double val) {
  try {
    (*this)->start_time_ = val;
  } catch (std::exception& e) {
    THROW_ERROR("set_start_time", e.what());
  }
}

double DaeBuilder::stop_time() const {
  try {
    return (*this)->stop_time_;
  } catch (std::exception& e) {
    THROW_ERROR("stop_time", e.what());
    return nan;
  }
}

void DaeBuilder::set_stop_time(double val) {
  try {
    (*this)->stop_time_ = val;
  } catch (std::exception& e) {
    THROW_ERROR("set_stop_time", e.what());
  }
}

double DaeBuilder::tolerance() const {
  try {
    return (*this)->tolerance_;
  } catch (std::exception& e) {
    THROW_ERROR("tolerance", e.what());
    return nan;
  }
}

void DaeBuilder::set_tolerance(double val) {
  try {
    (*this)->tolerance_ = val;
  } catch (std::exception& e) {
    THROW_ERROR("set_tolerance", e.what());
  }
}

double DaeBuilder::step_size() const {
  try {
    return (*this)->step_size_;
  } catch (std::exception& e) {
    THROW_ERROR("step_size", e.what());
    return nan;
  }
}

void DaeBuilder::set_step_size(double val) {
  try {
    (*this)->step_size_ = val;
  } catch (std::exception& e) {
    THROW_ERROR("set_step_size", e.what());
  }
}

void DaeBuilder::add_lc(const std::string& name,
    const std::vector<std::string>& f_out) {
  try {
    (*this)->add_lc(name, f_out);
  } catch (std::exception& e) {
    THROW_ERROR("add_lc", e.what());
  }
}

Function DaeBuilder::create(const std::string& fname,
    const std::vector<std::string>& name_in,
    const std::vector<std::string>& name_out, bool sx, bool lifted_calls) const {
  try {
    return (*this)->create(fname, name_in, name_out, Dict(), sx, lifted_calls);
  } catch (std::exception& e) {
    THROW_ERROR("create", e.what());
    return Function();  // never reached
  }
}

Function DaeBuilder::create(const std::string& fname,
    const std::vector<std::string>& name_in,
    const std::vector<std::string>& name_out, const Dict& opts) const {
  try {
    return (*this)->create(fname, name_in, name_out, opts, false, false);
  } catch (std::exception& e) {
    THROW_ERROR("create", e.what());
    return Function();  // never reached
  }
}

Function DaeBuilder::create(const std::string& name, const Dict& opts) const {
  try {
    return (*this)->create(name, dyn_in(), dyn_out(), opts, false, false);
  } catch (std::exception& e) {
    THROW_ERROR("create", e.what());
    return Function(); // never reached
  }
}

Function DaeBuilder::add_fun(const Function& f) {
  try {
    return (*this)->add_fun(f);
  } catch (std::exception& e) {
    THROW_ERROR("add_fun", e.what());
    return Function();  // never reached
  }
}

Function DaeBuilder::add_fun(const std::string& name, const std::vector<std::string>& arg,
    const std::vector<std::string>& res, const Dict& opts) {
  try {
    return (*this)->add_fun(name, arg, res, opts);
  } catch (std::exception& e) {
    THROW_ERROR("add_fun", e.what());
    return Function();  // never reached
  }
}

Function DaeBuilder::add_fun(const std::string& name, const Importer& compiler, const Dict& opts) {
  casadi_assert(!has_fun(name), "Function '" + name + "' already exists");
  return add_fun(external(name, compiler, opts));
}

bool DaeBuilder::has_fun(const std::string& name) const {
  try {
    return (*this)->has_fun(name);
  } catch (std::exception& e) {
    THROW_ERROR("has_fun", e.what());
    return false;  // never reached
  }
}

Function DaeBuilder::fun(const std::string& name) const {
  try {
    return (*this)->fun(name);
  } catch (std::exception& e) {
    THROW_ERROR("fun", e.what());
    return Function();  // never reached
  }
}

void DaeBuilder::gather_fun(casadi_int max_depth) {
  try {
    // Get a function corresponding to all equations (no inputs)
    Function all_eq = (*this)->gather_eq();
    // Gather all functions
    std::vector<Function> allfun = all_eq.find_functions(max_depth);
    // Add to list of functions
    for (const Function& f : allfun) {
      if (has_fun(f.name())) {
        // Skip functions with duplicate names
        casadi_warning("Duplicate function: '" + f.name() + "', ignored");
      } else {
        // Add to list of functions
        add_fun(f);
      }
    }
  } catch (std::exception& e) {
    THROW_ERROR("gather_fun", e.what());
  }
}

std::vector<Function> DaeBuilder::fun() const {
  return (*this)->fun_;
}

Function DaeBuilder::oracle(bool sx, bool elim_w, bool lifted_calls) const {
  try {
    return (*this)->oracle(sx, elim_w, lifted_calls);
  } catch (std::exception& e) {
    THROW_ERROR("oracle", e.what());
    return Function(); // never reached
  }
}

Function DaeBuilder::dependent_fun(const std::string& fname,
    const std::vector<std::string>& s_in,
    const std::vector<std::string>& s_out) const {
  try {
    return (*this)->dependent_fun(fname, s_in, s_out);
  } catch (std::exception& e) {
    THROW_ERROR("dependent_fun", e.what());
    return Function(); // never reached
  }
}

Function DaeBuilder::transition(const std::string& fname, casadi_int index) const {
  try {
    return (*this)->transition(fname, index);
  } catch (std::exception& e) {
    THROW_ERROR("transition", e.what());
    return Function(); // never reached
  }
}

Function DaeBuilder::transition(const std::string& fname) const {
  try {
    return (*this)->transition(fname);
  } catch (std::exception& e) {
    THROW_ERROR("transition", e.what());
    return Function(); // never reached
  }
}

bool DaeBuilder::test_cast(const SharedObjectInternal* ptr) {
  return dynamic_cast<const DaeBuilderInternal*>(ptr) != nullptr;
}

DaeBuilderInternal* DaeBuilder::operator->() {
  return static_cast<DaeBuilderInternal*>(SharedObject::operator->());
}

const DaeBuilderInternal* DaeBuilder::operator->() const {
  return static_cast<const DaeBuilderInternal*>(SharedObject::operator->());
}

MX DaeBuilder::der(const MX& v) const {
  try {
    return (*this)->der(v);
  } catch (std::exception& e) {
    THROW_ERROR("der", e.what());
    return MX();  // never reached
  }
}

double DaeBuilder::attribute(const std::string& a, const std::string& name) const {
  try {
    return (*this)->attribute(to_enum<Attribute>(a), name);
  } catch (std::exception& e) {
    THROW_ERROR("attribute", e.what());
    return 0; // never reached
  }
}

std::vector<double> DaeBuilder::attribute(const std::string& a,
    const std::vector<std::string>& name) const {
  try {
    return (*this)->attribute(to_enum<Attribute>(a), name);
  } catch (std::exception& e) {
    THROW_ERROR("attribute", e.what());
    return {}; // never reached
  }
}

void DaeBuilder::set_attribute(const std::string& a, const std::string& name, double val) {
  try {
    (*this)->set_attribute(to_enum<Attribute>(a), name, val);
  } catch (std::exception& e) {
    THROW_ERROR("set_attribute", e.what());
  }
}

void DaeBuilder::set_attribute(const std::string& a, const std::vector<std::string>& name,
    const std::vector<double>& val) {
  try {
    (*this)->set_attribute(to_enum<Attribute>(a), name, val);
  } catch (std::exception& e) {
    THROW_ERROR("set_attribute", e.what());
  }
}

double DaeBuilder::min(const std::string& name) const {
  try {
    return (*this)->variable(name).min;
  } catch (std::exception& e) {
    THROW_ERROR("min", e.what());
    return 0; // never reached
  }
}

std::vector<double> DaeBuilder::min(const std::vector<std::string>& name) const {
  try {
    return (*this)->attribute(Attribute::MIN, name);
  } catch (std::exception& e) {
    THROW_ERROR("min", e.what());
    return {}; // never reached
  }
}

void DaeBuilder::set_min(const std::string& name, double val) {
  try {
    (*this)->variable(name).min = val;
  } catch (std::exception& e) {
    THROW_ERROR("set_min", e.what());
  }
}

void DaeBuilder::set_min(const std::vector<std::string>& name, const std::vector<double>& val) {
  try {
    (*this)->set_attribute(Attribute::MIN, name, val);
  } catch (std::exception& e) {
    THROW_ERROR("set_min", e.what());
  }
}

double DaeBuilder::max(const std::string& name) const {
  try {
    return (*this)->variable(name).max;
  } catch (std::exception& e) {
    THROW_ERROR("max", e.what());
    return 0; // never reached
  }
}

std::vector<double> DaeBuilder::max(const std::vector<std::string>& name) const {
  try {
    return (*this)->attribute(Attribute::MAX, name);
  } catch (std::exception& e) {
    THROW_ERROR("max", e.what());
    return {}; // never reached
  }
}

void DaeBuilder::set_max(const std::string& name, double val) {
  try {
    (*this)->variable(name).max = val;
  } catch (std::exception& e) {
    THROW_ERROR("set_max", e.what());
  }
}

void DaeBuilder::set_max(const std::vector<std::string>& name, const std::vector<double>& val) {
  try {
    (*this)->set_attribute(Attribute::MAX, name, val);
  } catch (std::exception& e) {
    THROW_ERROR("set_max", e.what());
  }
}

double DaeBuilder::nominal(const std::string& name) const {
  try {
    return (*this)->variable(name).nominal;
  } catch (std::exception& e) {
    THROW_ERROR("nominal", e.what());
    return 0; // never reached
  }
}

std::vector<double> DaeBuilder::nominal(const std::vector<std::string>& name) const {
  try {
    return (*this)->attribute(Attribute::NOMINAL, name);
  } catch (std::exception& e) {
    THROW_ERROR("nominal", e.what());
    return {}; // never reached
  }
}

void DaeBuilder::set_nominal(const std::string& name, double val) {
  try {
    (*this)->variable(name).nominal = val;
  } catch (std::exception& e) {
    THROW_ERROR("set_nominal", e.what());
  }
}

void DaeBuilder::set_nominal(const std::vector<std::string>& name, const std::vector<double>& val) {
  try {
    (*this)->set_attribute(Attribute::NOMINAL, name, val);
  } catch (std::exception& e) {
    THROW_ERROR("set_mininal", e.what());
  }
}

std::vector<double> DaeBuilder::start(const std::string& name) const {
  try {
    return (*this)->attribute(Attribute::START, std::vector<std::string>{name});
  } catch (std::exception& e) {
    THROW_ERROR("start", e.what());
    return {}; // never reached
  }
}

std::vector<double> DaeBuilder::start(const std::vector<std::string>& name) const {
  try {
    return (*this)->attribute(Attribute::START, name);
  } catch (std::exception& e) {
    THROW_ERROR("start", e.what());
    return {}; // never reached
  }
}

void DaeBuilder::set_start(const std::string& name, double val) {
  try {
    (*this)->set_attribute(Attribute::START, name, val);
  } catch (std::exception& e) {
    THROW_ERROR("set_start", e.what());
  }
}

void DaeBuilder::set_start(const std::string& name, const std::vector<double>& val) {
  try {
    (*this)->set_attribute(Attribute::START, std::vector<std::string>{name}, val);
  } catch (std::exception& e) {
    THROW_ERROR("set_start", e.what());
  }
}

void DaeBuilder::set_start(const std::vector<std::string>& name, const std::vector<double>& val) {
  try {
    (*this)->set_attribute(Attribute::START, name, val);
  } catch (std::exception& e) {
    THROW_ERROR("set_start", e.what());
  }
}

void DaeBuilder::reset() {
  try {
    (*this)->reset();
  } catch (std::exception& e) {
    THROW_ERROR("reset", e.what());
  }
}

void DaeBuilder::set(const std::string& name, double val) {
  try {
    (*this)->set_attribute(Attribute::VALUE, name, val);
  } catch (std::exception& e) {
    THROW_ERROR("set", e.what());
  }
}

void DaeBuilder::set(const std::string& name, const std::string& val) {
  try {
    (*this)->set_string_attribute(Attribute::STRINGVALUE, name, val);
  } catch (std::exception& e) {
    THROW_ERROR("set", e.what());
  }
}

void DaeBuilder::set(const std::vector<std::string>& name, const std::vector<double>& val) {
  try {
    (*this)->set_attribute(Attribute::VALUE, name, val);
  } catch (std::exception& e) {
    THROW_ERROR("set", e.what());
  }
}

void DaeBuilder::set(const std::vector<std::string>& name,
    const std::vector<std::string>& val) {
  try {
    (*this)->set_string_attribute(Attribute::STRINGVALUE, name, val);
  } catch (std::exception& e) {
    THROW_ERROR("set", e.what());
  }
}

GenericType DaeBuilder::get(const std::string& name) const {
  return get(std::vector<std::string>{name}).front();
}

std::vector<GenericType> DaeBuilder::get(const std::vector<std::string>& name) const {
  try {
    // Create a temporary FmuFunction instance
    Function f = create(this->name() + "_get", {}, {}, Dict{{"aux", name}});
    // Get the stats
    Dict stats = f.stats().at("aux");
    // Return in the same order as inputs
    std::vector<GenericType> ret;
    ret.reserve(name.size());
    for (const std::string& n : name) ret.push_back(stats.at(n));
    return ret;
  } catch (std::exception& e) {
    THROW_ERROR("get", e.what());
    return {};
  }
}

Sparsity DaeBuilder::jac_sparsity(const std::vector<std::string>& onames,
    const std::vector<std::string>& inames) const {
  try {
    return (*this)->jac_sparsity((*this)->find(onames), (*this)->find(inames));
  } catch (std::exception& e) {
    THROW_ERROR("jac_sparsity", e.what());
    return Sparsity();  // never reached
  }
}

} // namespace casadi

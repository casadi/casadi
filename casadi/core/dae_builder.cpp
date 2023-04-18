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

const MX& DaeBuilder::t() const {
  return var((*this)->t_.at(0));
}

std::vector<std::string> DaeBuilder::x() const {
  return name((*this)->x_);
}

std::vector<MX> DaeBuilder::ode() const {
  return (*this)->ode();
}

std::vector<std::string> DaeBuilder::z() const {
  return name((*this)->z_);
}

std::vector<MX> DaeBuilder::alg() const {
  return (*this)->alg();
}

std::vector<std::string> DaeBuilder::q() const {
  return name((*this)->q_);
}

std::vector<MX> DaeBuilder::quad() const {
  return (*this)->quad();
}

std::vector<std::string> DaeBuilder::y() const {
  return name((*this)->y_);
}

std::vector<MX> DaeBuilder::ydef() const {
  return (*this)->ydef();
}

std::vector<std::string> DaeBuilder::u() const {
  return name((*this)->u_);
}

std::vector<std::string> DaeBuilder::p() const {
  return name((*this)->p_);
}

std::vector<std::string> DaeBuilder::c() const {
  return name((*this)->c_);
}

std::vector<MX> DaeBuilder::cdef() const {
  return (*this)->cdef();
}

std::vector<std::string> DaeBuilder::d() const {
  return name((*this)->d_);
}

std::vector<MX> DaeBuilder::ddef() const {
  return (*this)->ddef();
}

std::vector<std::string> DaeBuilder::w() const {
  return name((*this)->w_);
}

std::vector<MX> DaeBuilder::wdef() const {
  return (*this)->wdef();
}

const std::vector<MX>& DaeBuilder::aux() const {
  return (*this)->aux_;
}

const std::vector<MX>& DaeBuilder::init_lhs() const {
  return (*this)->init_lhs_;
}

const std::vector<MX>& DaeBuilder::init_rhs() const {
  return (*this)->init_rhs_;
}

const std::vector<MX>& DaeBuilder::when_cond() const {
  return (*this)->when_cond_;
}

const std::vector<MX>& DaeBuilder::when_lhs() const {
  return (*this)->when_lhs_;
}

const std::vector<MX>& DaeBuilder::when_rhs() const {
  return (*this)->when_rhs_;
}

std::vector<std::string> DaeBuilder::outputs() const {
  try {
    return name((*this)->outputs_);
  } catch (std::exception& e) {
    THROW_ERROR("outputs", e.what());
    return {};  // never reached
  }
}

std::vector<std::string> DaeBuilder::derivatives() const {
  try {
    return name((*this)->derivatives_);
  } catch (std::exception& e) {
    THROW_ERROR("derivatives", e.what());
    return {};  // never reached
  }
}

std::vector<std::string> DaeBuilder::initial_unknowns() const {
  try {
    return name((*this)->initial_unknowns_);
  } catch (std::exception& e) {
    THROW_ERROR("initial_unknowns", e.what());
    return {};  // never reached
  }
}

bool DaeBuilder::has_t() const {
  return !(*this)->t_.empty();
}

casadi_int DaeBuilder::nx() const {
  return (*this)->x_.size();
}

casadi_int DaeBuilder::nz() const {
  return (*this)->z_.size();
}

casadi_int DaeBuilder::nq() const {
  return (*this)->q_.size();
}

casadi_int DaeBuilder::ny() const {
  return (*this)->y_.size();
}

casadi_int DaeBuilder::nu() const {
  return (*this)->u_.size();
}

casadi_int DaeBuilder::np() const {
  return (*this)->p_.size();
}

casadi_int DaeBuilder::nc() const {
  return (*this)->c_.size();
}

casadi_int DaeBuilder::nd() const {
  return (*this)->d_.size();
}

casadi_int DaeBuilder::nw() const {
  return (*this)->w_.size();
}

void DaeBuilder::load_fmi_description(const std::string& filename) {
  try {
    (*this)->load_fmi_description(filename);
  } catch (std::exception& e) {
    THROW_ERROR("load_fmi_description", e.what());
  }
}

bool DaeBuilder::provides_directional_derivative() const {
  try {
    casadi_assert(!(*this)->symbolic_, "Functionality only applies to imported standard FMUs");
    return (*this)->provides_directional_derivative_;
  } catch (std::exception& e) {
    THROW_ERROR("provides_directional_derivative", e.what());
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

bool DaeBuilder::has_variable(const std::string& name) const {
  try {
    return (*this)->has_variable(name);
  } catch (std::exception& e) {
    THROW_ERROR("has_variable", e.what());
    return false;  // never reached
  }
}

std::vector<std::string> DaeBuilder::all_variables() const {
  try {
    return (*this)->all_variables();
  } catch (std::exception& e) {
    THROW_ERROR("all_variables", e.what());
    return {};  // never reached
  }
}

Variable& DaeBuilder::new_variable(const std::string& name, casadi_int numel) {
  try {
    return (*this)->new_variable(name, numel);
  } catch (std::exception& e) {
    THROW_ERROR("new_variable", e.what());
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
  (*this)->t_.push_back(find(name));
}

void DaeBuilder::register_p(const std::string& name) {
  (*this)->p_.push_back(find(name));
}

void DaeBuilder::register_u(const std::string& name) {
  (*this)->u_.push_back(find(name));
}

void DaeBuilder::register_x(const std::string& name) {
  (*this)->x_.push_back(find(name));
}

void DaeBuilder::register_z(const std::string& name) {
  (*this)->z_.push_back(find(name));
}

void DaeBuilder::register_q(const std::string& name) {
  (*this)->q_.push_back(find(name));
}

void DaeBuilder::register_c(const std::string& name) {
  (*this)->c_.push_back(find(name));
}

void DaeBuilder::register_d(const std::string& name) {
  (*this)->d_.push_back(find(name));
}

void DaeBuilder::register_w(const std::string& name) {
  (*this)->w_.push_back(find(name));
}

void DaeBuilder::register_y(const std::string& name) {
  (*this)->y_.push_back(find(name));
}

void DaeBuilder::set_z(const std::vector<std::string>& name, const std::vector<std::string>& alg) {
  try {
    // Consistency check
    if (!alg.empty()) {
      casadi_assert(alg.size() == name.size(), "Inconsistent number of algebraic variables");
    }

    // Update z
    set_all("z", name);

    // Update algebraic equations
    if (!alg.empty()) {
      for (size_t i = 0; i < name.size(); ++i) variable(name[i]).alg = find(alg[i]);
    }
  } catch (std::exception& e) {
    THROW_ERROR("set_z", e.what());
  }
}

void DaeBuilder::clear_all(const std::string& v) {
  try {
    (*this)->clear_all(v);
  } catch (std::exception& e) {
    THROW_ERROR("clear_all", e.what());
  }
}

void DaeBuilder::set_all(const std::string& v, const std::vector<std::string>& name) {
  try {
    (*this)->set_all(v, name);
  } catch (std::exception& e) {
    THROW_ERROR("set_all", e.what());
  }
}

MX DaeBuilder::add_t(const std::string& name) {
  casadi_assert((*this)->t_.empty(), "'t' already defined");
  size_t new_t = add_variable_new(name);
  (*this)->t_.push_back(new_t);
  return var(new_t);
}

MX DaeBuilder::add_p(const std::string& name) {
  try {
    return (*this)->add_p(name);
  } catch (std::exception& e) {
    THROW_ERROR("add_p", e.what());
    return MX();
  }
}

MX DaeBuilder::add_u(const std::string& name) {
  try {
    return (*this)->add_u(name);
  } catch (std::exception& e) {
    THROW_ERROR("add_u", e.what());
    return MX();
  }
}

MX DaeBuilder::add_x(const std::string& name) {
  try {
    return (*this)->add_x(name);
  } catch (std::exception& e) {
    THROW_ERROR("add_x", e.what());
    return MX();
  }
}

MX DaeBuilder::add_z(const std::string& name) {
  try {
    return (*this)->add_z(name);
  } catch (std::exception& e) {
    THROW_ERROR("add_z", e.what());
    return MX();
  }
}

MX DaeBuilder::add_q(const std::string& name) {
  try {
    return (*this)->add_q(name);
  } catch (std::exception& e) {
    THROW_ERROR("add_q", e.what());
    return MX();
  }
}

MX DaeBuilder::add_c(const std::string& name, const MX& new_cdef) {
  try {
    return (*this)->add_c(name, new_cdef);
  } catch (std::exception& e) {
    THROW_ERROR("add_c", e.what());
    return MX();
  }
}

MX DaeBuilder::add_d(const std::string& name, const MX& new_ddef) {
  try {
    return (*this)->add_d(name, new_ddef);
  } catch (std::exception& e) {
    THROW_ERROR("add_d", e.what());
    return MX();
  }
}

MX DaeBuilder::add_w(const std::string& name, const MX& new_wdef) {
  try {
    return (*this)->add_w(name, new_wdef);
  } catch (std::exception& e) {
    THROW_ERROR("add_w", e.what());
    return MX();
  }
}

MX DaeBuilder::add_y(const std::string& name, const MX& new_ydef) {
  try {
    return (*this)->add_y(name, new_ydef);
  } catch (std::exception& e) {
    THROW_ERROR("add_y", e.what());
    return MX();
  }
}

MX DaeBuilder::add_aux(const std::string& name, casadi_int n) {
  MX new_aux = add_variable(name, n);
  (*this)->aux_.push_back(new_aux);
  return new_aux;
}

void DaeBuilder::add_init(const MX& lhs, const MX& rhs) {
  (*this)->init_lhs_.push_back(lhs);
  (*this)->init_rhs_.push_back(rhs);
}

void DaeBuilder::add_when(const MX& cond, const MX& lhs, const MX& rhs) {
  (*this)->when_cond_.push_back(cond);
  (*this)->when_lhs_.push_back(lhs);
  (*this)->when_rhs_.push_back(rhs);
}

void DaeBuilder::set_ode(const std::string& name, const MX& ode_rhs) {
  try {
    (*this)->set_ode(name, ode_rhs);
  } catch (std::exception& e) {
    THROW_ERROR("set_ode", e.what());
  }
}

void DaeBuilder::set_alg(const std::string& name, const MX& alg_rhs) {
  try {
    (*this)->set_alg(name, alg_rhs);
  } catch (std::exception& e) {
    THROW_ERROR("set_alg", e.what());
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
    return variable(name).v;
  } catch (std::exception& e) {
    THROW_ERROR("var", e.what());
    return MX();  // never reached
  }
}

std::string DaeBuilder::der(const std::string& name) const {
  try {
    // Get variable index
    size_t ind = find(name);
    // Differentiate
    ind = variable(ind).der;
    casadi_assert(ind != size_t(-1), "No derivative expression for " + name);
    // Return name
    return variable(ind).name;
  } catch (std::exception& e) {
    THROW_ERROR("der", e.what());
    return std::string();  // never reached
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

MX DaeBuilder::beq(const std::string& name) const {
  try {
    return variable(name).beq;
  } catch (std::exception& e) {
    THROW_ERROR("beq", e.what());
    return MX();  // never reached
  }
}

void DaeBuilder::set_beq(const std::string& name, const MX& val) {
  try {
    variable(name).beq = val;
  } catch (std::exception& e) {
    THROW_ERROR("set_beq", e.what());
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
  return variable(name).value_reference;
}

void DaeBuilder::set_value_reference(const std::string& name, casadi_int val) {
  variable(name).value_reference = val;
}

std::string DaeBuilder::description(const std::string& name) const {
  return variable(name).description;
}

void DaeBuilder::set_description(const std::string& name, const std::string& val) {
  variable(name).description = val;
}

std::string DaeBuilder::type(const std::string& name, casadi_int fmi_version) const {
  // Check version
  casadi_assert(fmi_version == 2 || fmi_version == 3, "Only FMI version 2 or 3 supported");
  // Handle FMI 2
  if (fmi_version == 2) {
    return to_string(to_fmi2(variable(name).type));
  }
  // Assume FMI 3
  return to_string(variable(name).type);
}

void DaeBuilder::set_type(const std::string& name, const std::string& val) {
  // Fallback to FMI 2, if necessary
  if (has_enum<TypeFmi2>(val) && !has_enum<Type>(val)) {
    variable(name).type = from_fmi2(to_enum<TypeFmi2>(val));
  }
  // Assume FMI 3
  variable(name).type = to_enum<Type>(val);
}

std::string DaeBuilder::causality(const std::string& name) const {
  return to_string(variable(name).causality);
}

void DaeBuilder::set_causality(const std::string& name, const std::string& val) {
  variable(name).causality = to_enum<Causality>(val);
}

std::string DaeBuilder::variability(const std::string& name) const {
  return to_string(variable(name).variability);
}

void DaeBuilder::set_variability(const std::string& name, const std::string& val) {
  variable(name).variability = to_enum<Variability>(val);
}

std::string DaeBuilder::initial(const std::string& name) const {
  return to_string(variable(name).initial);
}

void DaeBuilder::set_initial(const std::string& name, const std::string& val) {
  variable(name).initial = to_enum<Initial>(val);
}

std::string DaeBuilder::unit(const std::string& name) const {
  return variable(name).unit;
}

void DaeBuilder::set_unit(const std::string& name, const std::string& val) {
  variable(name).unit = val;
}

std::string DaeBuilder::display_unit(const std::string& name) const {
  return variable(name).display_unit;
}

void DaeBuilder::set_display_unit(const std::string& name, const std::string& val) {
  variable(name).display_unit = val;
}

casadi_int DaeBuilder::numel(const std::string& name) const {
  return variable(name).numel;
}

std::vector<casadi_int> DaeBuilder::dimension(const std::string& name) const {
  return variable(name).dimension;
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

bool DaeBuilder::test_cast(const SharedObjectInternal* ptr) {
  return dynamic_cast<const DaeBuilderInternal*>(ptr) != nullptr;
}

DaeBuilderInternal* DaeBuilder::operator->() {
  return static_cast<DaeBuilderInternal*>(SharedObject::operator->());
}

const DaeBuilderInternal* DaeBuilder::operator->() const {
  return static_cast<const DaeBuilderInternal*>(SharedObject::operator->());
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
    return variable(ind).name;
  } catch (std::exception& e) {
    THROW_ERROR("name", e.what());
    static std::string dummy;
    return dummy; // never reached
  }
}

std::vector<std::string> DaeBuilder::name(const std::vector<size_t>& ind) const {
  try {
    std::vector<std::string> r(ind.size());
    for (size_t i = 0; i < r.size(); ++i) r[i] = name(ind[i]);
    return r;
  } catch (std::exception& e) {
    THROW_ERROR("name", e.what());
    return {}; // never reached
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
    return variable(name).min;
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
    variable(name).min = val;
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
    return variable(name).max;
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
    variable(name).max = val;
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
    return variable(name).nominal;
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
    variable(name).nominal = val;
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

double DaeBuilder::start(const std::string& name) const {
  try {
    casadi_assert(numel(name) == 1, "Variable " + name + " is not scalar");
    return variable(name).start.front();
  } catch (std::exception& e) {
    THROW_ERROR("start", e.what());
    return 0; // never reached
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
    return (*this)->jac_sparsity(find(onames), find(inames));
  } catch (std::exception& e) {
    THROW_ERROR("jac_sparsity", e.what());
    return Sparsity();  // never reached
  }
}

} // namespace casadi

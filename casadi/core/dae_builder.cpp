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

namespace casadi {

// Throw informative error message
#define THROW_ERROR(FNAME, WHAT) \
throw CasadiException("Error in DaeBuilder::" FNAME " for '" + this->name() \
  + "'  at " + CASADI_WHERE + ":\n" + std::string(WHAT));

DaeBuilder::DaeBuilder() {
}

DaeBuilder::DaeBuilder(const std::string& name) {
  own(new DaeBuilderInternal(name));
}

const std::string& DaeBuilder::name() const {
  return (*this)->name_;
}

const MX& DaeBuilder::t() const {
  casadi_assert(has_t(), "No time variable");
  return (*this)->t();
}

const std::vector<MX>& DaeBuilder::x() const {
  return (*this)->x_;
}

const std::vector<MX>& DaeBuilder::ode() const {
  return (*this)->ode_;
}

const std::vector<MX>& DaeBuilder::z() const {
  return (*this)->z_;
}

const std::vector<MX>& DaeBuilder::alg() const {
  return (*this)->alg_;
}

const std::vector<MX>& DaeBuilder::q() const {
  return (*this)->q_;
}

const std::vector<MX>& DaeBuilder::quad() const {
  return (*this)->quad_;
}

const std::vector<MX>& DaeBuilder::y() const {
  return (*this)->y_;
}

std::vector<MX> DaeBuilder::ydef() const {
  return (*this)->ydef_;
}

const std::vector<MX>& DaeBuilder::u() const {
  return (*this)->u_;
}

const std::vector<MX>& DaeBuilder::p() const {
  return (*this)->p_;
}

const std::vector<MX>& DaeBuilder::c() const {
  return (*this)->c_;
}

std::vector<MX> DaeBuilder::cdef() const {
  return (*this)->cdef_;
}

const std::vector<MX>& DaeBuilder::d() const {
  return (*this)->d_;
}

std::vector<MX> DaeBuilder::ddef() const {
  return (*this)->ddef_;
}

const std::vector<MX>& DaeBuilder::w() const {
  return (*this)->w_;
}

std::vector<MX> DaeBuilder::wdef() const {
  return (*this)->wdef_;
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

bool DaeBuilder::has_t() const {
  return (*this)->has_t();
}

void DaeBuilder::parse_fmi(const std::string& filename) {
  try {
    (*this)->parse_fmi(filename);
  } catch (std::exception& e) {
    THROW_ERROR("parse_fmi", e.what());
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

bool DaeBuilder::has_variable(const std::string& name) const {
  try {
    return (*this)->has_variable(name);
  } catch (std::exception& e) {
    THROW_ERROR("has_variable", e.what());
    return false;  // never reached
  }
}

void DaeBuilder::add_variable(const std::string& name, const Variable& var) {
  try {
    return (*this)->add_variable(name, var);
  } catch (std::exception& e) {
    THROW_ERROR("add_variable", e.what());
  }
}

MX DaeBuilder::add_variable(const std::string& name, casadi_int n) {
  return add_variable(name, Sparsity::dense(n));
}

MX DaeBuilder::add_variable(const std::string& name, const Sparsity& sp) {
  Variable v(name);
  v.v = MX::sym(name, sp);
  add_variable(name, v);
  return v.v;
}

void DaeBuilder::add_variable(const MX& new_v) {
  Variable v(new_v.name());
  v.v = new_v;
  add_variable(new_v.name(), v);
}

void DaeBuilder::register_p(const MX& new_p) {
  // Consistency checks
  casadi_assert(has_variable(new_p.name()), "No such variable: " + new_p.name());
  // Add to list
  (*this)->p_.push_back(new_p);
}

void DaeBuilder::register_u(const MX& new_u) {
  // Consistency checks
  casadi_assert(has_variable(new_u.name()), "No such variable: " + new_u.name());
  // Add to list
  (*this)->u_.push_back(new_u);
}

void DaeBuilder::register_x(const MX& new_x) {
  // Consistency checks
  casadi_assert(has_variable(new_x.name()), "No such variable: " + new_x.name());
  // Add to list
  (*this)->x_.push_back(new_x);
}

void DaeBuilder::register_z(const MX& new_z) {
  // Consistency checks
  casadi_assert(has_variable(new_z.name()), "No such variable: " + new_z.name());
  // Add to list
  (*this)->z_.push_back(new_z);
}

void DaeBuilder::register_t(const MX& new_t) {
  try {
    (*this)->register_t(new_t);
  } catch (std::exception& e) {
    THROW_ERROR("register_t", e.what());
  }
}

void DaeBuilder::register_c(const MX& new_c, const MX& new_cdef) {
  // Consistency check
  casadi_assert(new_c.sparsity() == new_cdef.sparsity(), "Mismatching sparsity");
  casadi_assert(has_variable(new_c.name()), "No such variable: " + new_c.name());
  // Add to lists
  (*this)->c_.push_back(new_c);
  (*this)->cdef_.push_back(new_cdef);
}

void DaeBuilder::register_d(const MX& new_d, const MX& new_ddef) {
  // Consistency checks
  casadi_assert(new_d.sparsity() == new_ddef.sparsity(), "Mismatching sparsity");
  casadi_assert(has_variable(new_d.name()), "No such variable: " + new_d.name());
  // Add to lists
  (*this)->d_.push_back(new_d);
  (*this)->ddef_.push_back(new_ddef);
}

void DaeBuilder::register_w(const MX& new_w, const MX& new_wdef) {
  // Consistency checks
  casadi_assert(new_w.sparsity() == new_wdef.sparsity(), "Mismatching sparsity");
  casadi_assert(has_variable(new_w.name()), "No such variable: " + new_w.name());
  // Add to lists
  (*this)->w_.push_back(new_w);
  (*this)->wdef_.push_back(new_wdef);
}

void DaeBuilder::register_y(const MX& new_y, const MX& new_ydef) {
  // Consistency checks
  casadi_assert(new_y.sparsity() == new_ydef.sparsity(), "Mismatching sparsity");
  casadi_assert(has_variable(new_y.name()), "No such variable: " + new_y.name());
  // Add to lists
  (*this)->y_.push_back(new_y);
  (*this)->ydef_.push_back(new_ydef);
}

void DaeBuilder::clear_in(const std::string& v) {
  try {
    (*this)->clear_in(v);
  } catch (std::exception& e) {
    THROW_ERROR("clear_in", e.what());
  }
}

void DaeBuilder::clear_out(const std::string& v) {
  try {
    (*this)->clear_out(v);
  } catch (std::exception& e) {
    THROW_ERROR("clear_out", e.what());
  }
}

MX DaeBuilder::add_t(const std::string& name) {
  try {
    return (*this)->add_t(name);
  } catch (std::exception& e) {
    THROW_ERROR("add_t", e.what());
    return MX();
  }
}

MX DaeBuilder::add_p(const std::string& name, casadi_int n) {
  try {
    if (name.empty()) return add_p("p" + str(p().size()), n);
    return (*this)->add_p(name, n);
  } catch (std::exception& e) {
    THROW_ERROR("add_p", e.what());
    return MX();
  }
}

MX DaeBuilder::add_u(const std::string& name, casadi_int n) {
  try {
    if (name.empty()) return add_u("u" + str(u().size()), n);
    return (*this)->add_p(name, n);
  } catch (std::exception& e) {
    THROW_ERROR("add_u", e.what());
    return MX();
  }
}

MX DaeBuilder::add_x(const std::string& name, casadi_int n) {
  try {
    if (name.empty()) return add_x("x" + str(x().size()), n);
    return (*this)->add_x(name, n);
  } catch (std::exception& e) {
    THROW_ERROR("add_x", e.what());
    return MX();
  }
}

MX DaeBuilder::add_z(const std::string& name, casadi_int n) {
  try {
    if (name.empty()) return add_z("z" + str(z().size()), n);
    return (*this)->add_z(name, n);
  } catch (std::exception& e) {
    THROW_ERROR("add_z", e.what());
    return MX();
  }
}

MX DaeBuilder::add_q(const std::string& name, casadi_int n) {
  try {
    if (name.empty()) return add_q("q" + str(q().size()), n);
    return (*this)->add_x(name, n);
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
  if (name.empty()) return add_aux("aux" + str(aux().size()), n);
  MX new_aux = add_variable(name, n);
  (*this)->aux_.push_back(new_aux);
  return new_aux;
}

void DaeBuilder::add_init(const MX& lhs, const MX& rhs) {
  (*this)->init_lhs_.push_back(lhs);
  (*this)->init_rhs_.push_back(rhs);
}

void DaeBuilder::add_ode(const std::string& name, const MX& new_ode) {
  (*this)->ode_.push_back(new_ode);
  (*this)->clear_cache_ = true;
}

void DaeBuilder::add_alg(const std::string& name, const MX& new_alg) {
  (*this)->alg_.push_back(new_alg);
  (*this)->clear_cache_ = true;
}

void DaeBuilder::add_quad(const std::string& name, const MX& new_quad) {
  (*this)->quad_.push_back(new_quad);
  (*this)->clear_cache_ = true;
}

void DaeBuilder::sanity_check() const {
  try {
    (*this)->sanity_check();
  } catch (std::exception& e) {
    THROW_ERROR("sanity_check", e.what());
  }
}

MX DaeBuilder::var(const std::string& name) const {
  return variable(name).v;
}

MX DaeBuilder::der(const std::string& name) const {
  return (*this)->variables_.at(variable(name).derivative).v;
}

MX DaeBuilder::der(const MX& var) const {
  casadi_assert_dev(var.is_column() && var.is_symbolic());
  return der(var.name());
}

void DaeBuilder::eliminate_w() {
  try {
    (*this)->eliminate_w();
  } catch (std::exception& e) {
    THROW_ERROR("eliminate_w", e.what());
  }
}

void DaeBuilder::lift(bool lift_shared, bool lift_calls, bool inline_calls) {
  try {
    (*this)->lift(lift_shared, lift_calls, inline_calls);
  } catch (std::exception& e) {
    THROW_ERROR("lift", e.what());
  }
}

std::string DaeBuilder::description(const std::string& name) const {
  return variable(name).description;
}

void DaeBuilder::set_description(const std::string& name, const std::string& val) {
  variable(name).description = val;
}

std::string DaeBuilder::type(const std::string& name) const {
  return to_string(variable(name).type);
}

void DaeBuilder::set_type(const std::string& name, const std::string& val) {
  variable(name).type = to_enum<Variable::Type>(val);
}

std::string DaeBuilder::causality(const std::string& name) const {
  return to_string(variable(name).causality);
}

void DaeBuilder::set_causality(const std::string& name, const std::string& val) {
  variable(name).causality = to_enum<Variable::Causality>(val);
}

std::string DaeBuilder::variability(const std::string& name) const {
  return to_string(variable(name).variability);
}

void DaeBuilder::set_variability(const std::string& name, const std::string& val) {
  variable(name).variability = to_enum<Variable::Variability>(val);
}

std::string DaeBuilder::initial(const std::string& name) const {
  return to_string(variable(name).initial);
}

void DaeBuilder::set_initial(const std::string& name, const std::string& val) {
  variable(name).initial = to_enum<Variable::Initial>(val);
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

MX DaeBuilder::nominal(const std::string& name) const {
  return variable(name).nominal;
}

void DaeBuilder::set_nominal(const std::string& name, const MX& val) {
  variable(name).nominal = val;
}

MX DaeBuilder::min(const std::string& name) const {
  return variable(name).min;
}

void DaeBuilder::set_min(const std::string& name, const MX& val) {
  variable(name).min = val;
}

MX DaeBuilder::max(const std::string& name) const {
  return variable(name).max;
}

void DaeBuilder::set_max(const std::string& name, const MX& val) {
  variable(name).max = val;
}

MX DaeBuilder::start(const std::string& name) const {
  return variable(name).start;
}

void DaeBuilder::set_start(const std::string& name, const MX& val) {
  variable(name).start = val;
}

const casadi::MX& DaeBuilder::binding_equation(const std::string& name) const {
  return variable(name).beq;
}

void DaeBuilder::set_binding_equation(const std::string& name, const MX& val) {
  variable(name).beq = val;
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
    const std::vector<std::string>& s_in,
    const std::vector<std::string>& s_out, bool sx, bool lifted_calls) const {
  try {
    return (*this)->create(fname, s_in, s_out, sx, lifted_calls);
  } catch (std::exception& e) {
    THROW_ERROR("create", e.what());
    return Function();  // never reached
  }
}

Function DaeBuilder::add_fun(const Function& f) {
  casadi_assert(!has_fun(f.name()), "Function '" + f.name() + "' already exists");
  (*this)->fun_.push_back(f);
  return f;
}

Function DaeBuilder::add_fun(const std::string& name,
                             const std::vector<std::string>& arg,
                             const std::vector<std::string>& res,
                             const Dict& opts) {
  casadi_assert(!has_fun(name), "Function '" + name + "' already exists");

  // Get inputs
  std::vector<MX> arg_ex, res_ex;
  for (auto&& s : arg) arg_ex.push_back(var(s));
  for (auto&& s : res) {
    // Find the binding expression FIXME(@jaeandersson)
    casadi_int v_ind;
    for (v_ind=0; v_ind < (*this)->w_.size(); ++v_ind) {
      if (s == (*this)->w_.at(v_ind).name()) {
        res_ex.push_back((*this)->wdef_.at(v_ind));
        break;
      }
    }
    casadi_assert(v_ind < (*this)->w_.size(), "Cannot find dependent '" + s + "'");
  }
  Function ret(name, arg_ex, res_ex, arg, res, opts);
  return add_fun(ret);
}

Function DaeBuilder::add_fun(const std::string& name, const Importer& compiler,
                             const Dict& opts) {
  casadi_assert(!has_fun(name), "Function '" + name + "' already exists");
  return add_fun(external(name, compiler, opts));
}

bool DaeBuilder::has_fun(const std::string& name) const {
  for (const Function& f : (*this)->fun_) {
    if (f.name()==name) return true;
  }
  return false;
}

Function DaeBuilder::fun(const std::string& name) const {
  casadi_assert(has_fun(name), "No such function: '" + name + "'");
  for (const Function& f : (*this)->fun_) {
    if (f.name()==name) return f;
  }
  return Function();
}

void DaeBuilder::gather_fun(casadi_int max_depth) {
  try {
    // Get a function corresponding to all equations (no inputs)
    Function all_eq = (*this)->gather_eq();
    // Gather all functions
    std::vector<Function> allfun = all_eq.find(max_depth);
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

Function DaeBuilder::attribute_fun(const std::string& fname,
    const std::vector<std::string>& s_in,
    const std::vector<std::string>& s_out) const {
  try {
    return (*this)->attribute_fun(fname, s_in, s_out);
  } catch (std::exception& e) {
    THROW_ERROR("attribute_fun", e.what());
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

void DaeBuilder::prune_d() {
  try {
    (*this)->prune_d();
  } catch (std::exception& e) {
    THROW_ERROR("prune_d", e.what());
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

DaeBuilderInternal* DaeBuilder::operator->() {
  return static_cast<DaeBuilderInternal*>(SharedObject::operator->());
}

const DaeBuilderInternal* DaeBuilder::operator->() const {
  return static_cast<const DaeBuilderInternal*>(SharedObject::operator->());
}

} // namespace casadi

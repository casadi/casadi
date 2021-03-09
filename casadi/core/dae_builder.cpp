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


#include "dae_builder.hpp"

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

std::string to_string(Variable::Type v) {
  switch (v) {
  case Variable::REAL: return "real";
  case Variable::INTEGER: return "integer";
  case Variable::BOOLEAN: return "boolean";
  case Variable::STRING: return "string";
  case Variable::ENUM: return "enum";
  default: break;
  }
  return "";
}

std::string to_string(Variable::Causality v) {
  switch (v) {
  case Variable::PARAMETER: return "parameter";
  case Variable::CALCULATED_PARAMETER: return "calculatedParameter";
  case Variable::INPUT: return "input";
  case Variable::OUTPUT: return "output";
  case Variable::LOCAL: return "local";
  case Variable::INDEPENDENT: return "independent";
  default: break;
  }
  return "";
}

std::string to_string(Variable::Variability v) {
  switch (v) {
  case Variable::CONSTANT: return "constant";
  case Variable::FIXED: return "fixed";
  case Variable::TUNABLE: return "tunable";
  case Variable::DISCRETE: return "discrete";
  case Variable::CONTINUOUS: return "continuous";
  default: break;
  }
  return "";
}

std::string to_string(Variable::Initial v) {
  switch (v) {
  case Variable::EXACT: return "exact";
  case Variable::APPROX: return "approx";
  case Variable::CALCULATED: return "calculated";
  case Variable::INITIAL_NA: return "initial_na";
  default: break;
  }
  return "";
}

CASADI_EXPORT std::string to_string(Variable::Attribute v) {
  switch (v) {
  case Variable::MIN: return "min";
  case Variable::MAX: return "max";
  case Variable::NOMINAL: return "nominal";
  case Variable::START: return "start";
  default: break;
  }
  return "";
}

Variable::Initial Variable::default_initial(Variable::Causality causality,
    Variable::Variability variability) {
  // According to table in FMI 2.0.2 specification, section 2.2.7
  switch (variability) {
  case CONSTANT:
    if (causality == OUTPUT || causality == LOCAL)
      return EXACT;
    break;
  case FIXED:
    // Fall-through
  case TUNABLE:
    if (causality == PARAMETER)
      return EXACT;
    else if (causality == CALCULATED_PARAMETER || causality == LOCAL)
      return CALCULATED;
    break;
  case DISCRETE:
  // Fall-through
  case CONTINUOUS:
    if (causality == OUTPUT || causality == LOCAL)
      return CALCULATED;
    break;
  default: break;
  }
  // Initial value not available
  return INITIAL_NA;
}

Variable::Variable(const std::string& name) : name(name),
    value_reference(-1), description(""),
    type(REAL), causality(LOCAL), variability(CONTINUOUS),
    unit(""), display_unit(""),
    min(-std::numeric_limits<double>::infinity()), max(std::numeric_limits<double>::infinity()),
    nominal(1.0), start(0.0), derivative(-1), antiderivative(-1), dependency(false) {
}

void Variable::disp(std::ostream &stream, bool more) const {
  stream << name;
}

DaeBuilder::DaeBuilder() {
  this->t = MX::sym("t");
  clear_cache_ = false;
}

void DaeBuilder::parse_fmi(const std::string& filename) {

  // Load
  XmlFile xml_file("tinyxml");
  XmlNode document = xml_file.parse(filename);

  // Number of variables before adding new ones
  size_t n_vars_before = variables_.size();

  // **** Add model variables ****
  {
    // Get a reference to the ModelVariables node
    const XmlNode& modvars = document[0]["ModelVariables"];

    // Add variables
    for (casadi_int i = 0; i < modvars.size(); ++i) {
      // Get a reference to the variable
      const XmlNode& vnode = modvars[i];

      // Name of variable, ensure unique
      std::string name = vnode.attribute<std::string>("name");
      casadi_assert(varind_.find(name) == varind_.end(), "Duplicate variable: " + name);

      // Create new variable
      Variable var(name);
      var.v = MX::sym(name);

      // Read common attributes, cf. FMI 2.0.2 specification, 2.2.7
      var.value_reference = vnode.attribute<casadi_int>("valueReference");
      var.description = vnode.attribute<std::string>("description", "");
      var.causality = to_enum<Variable::Causality>(
        vnode.attribute<std::string>("causality", "local"));
      var.variability = to_enum<Variable::Variability>(
        vnode.attribute<std::string>("variability", "continuous"));
      std::string initial_str = vnode.attribute<std::string>("initial", "");
      if (initial_str.empty()) {
        // Default value
        var.initial = Variable::default_initial(var.causality, var.variability);
      } else {
        // Consistency check
        casadi_assert(var.causality != Variable::INPUT && var.causality != Variable::INDEPENDENT,
          "The combination causality = '" + to_string(var.causality) + "', "
          "initial = '" + initial_str + "' is not allowed per FMI 2.0 specification.");
        // Value specified
        var.initial = to_enum<Variable::Initial>(initial_str);
      }
      // Other properties
      if (vnode.has_child("Real")) {
        const XmlNode& props = vnode["Real"];
        var.unit = props.attribute<std::string>("unit", var.unit);
        var.display_unit = props.attribute<std::string>("displayUnit", var.display_unit);
        var.min = props.attribute<double>("min", -inf);
        var.max = props.attribute<double>("max", inf);
        var.nominal = props.attribute<double>("nominal", 1.);
        var.start = props.attribute<double>("start", 0.);
        var.derivative = props.attribute<casadi_int>("derivative", var.derivative);
      }
      // Add to list of variables
      add_variable(name, var);
    }
    // Handle derivatives/antiderivatives
    for (auto it = variables_.begin() + n_vars_before; it != variables_.end(); ++it) {
      if (it->derivative >= 0) {
        // Add variable offset, make index 1
        it->derivative += n_vars_before - 1;
        // Set antiderivative
        variables_.at(it->derivative).antiderivative = it - variables_.begin();
      }
    }
  }

  // **** Process model structure ****
  if (document[0].has_child("ModelStructure")) {
    // Get a reference to the ModelVariables node
    const XmlNode& modst = document[0]["ModelStructure"];
    // Test both Outputs and Derivatives
    for (const char* dtype : {"Outputs", "Derivatives"}) {
      // Derivative variables
      if (modst.has_child(dtype)) {
        const XmlNode& outputs = modst[dtype];
        for (casadi_int i = 0; i < outputs.size(); ++i) {
          // Get a reference to the output
          const XmlNode& onode = outputs[i];
          // Read attribute
          casadi_int index = onode.attribute<casadi_int>("index", -1);
          casadi_assert(index >= 1, "Non-positive output index");
          // Convert to index in variables list
          index += n_vars_before - 1;
          // Get dependencies
          std::vector<casadi_int> dependencies = onode.attribute<std::vector<casadi_int>>(
            "dependencies", {});
          // Convert to indices in variables list
          for (casadi_int& d : dependencies) {
            // Consistency check, add offset
            casadi_assert(d >= 1, "Non-positive dependency index");
            d += n_vars_before - 1;
            // Mark corresponding variable as dependency
            variables_.at(d).dependency = true;
          }
        }
      }
    }
  }

  // **** Postprocess / sort variables ****
  for (auto it = variables_.begin() + n_vars_before; it != variables_.end(); ++it) {
    // Sort by types
    if (it->causality == Variable::INDEPENDENT) {
      // Independent (time) variable
      this->t = it->v;
    } else if (it->causality == Variable::INPUT) {
      this->u.push_back(it->v);
    } else if (it->variability == Variable::CONSTANT) {
      // Named constant
      this->c.push_back(it->v);
      this->cdef.push_back(it->start);
    } else if (it->variability == Variable::FIXED || it->variability == Variable::TUNABLE) {
      this->p.push_back(it->v);
    } else if (it->variability == Variable::CONTINUOUS) {
      if (it->antiderivative >= 0) {
        // Is the variable needed to calculate other states, algebraic variables?
        if (it->dependency) {
          // Add to list of differential equations
          this->x.push_back(it->v);
          add_ode("ode_" + it->name, variables_.at(it->antiderivative).v);
        } else {
          // Add to list of quadrature equations
          this->q.push_back(it->v);
          add_quad("quad_" + it->name, variables_.at(it->antiderivative).v);
        }
      } else if (it->dependency || it->derivative >= 0) {
        // Add to list of algebraic equations
        this->z.push_back(it->v);
        add_alg("alg_" + it->name, it->v - nan);
      }
      // Is it (also) an output variable?
      if (it->causality == Variable::OUTPUT) {
        this->y.push_back(it->v);
        this->ydef.push_back(it->v);
      }
    } else if (it->dependency) {
      casadi_warning("Cannot sort " + it->name);
    }
  }
}

Variable& DaeBuilder::read_variable(const XmlNode& node) {
  // Qualified name
  std::string qn = qualified_name(node);

  // Find and return the variable
  return variable(qn);
}

MX DaeBuilder::read_expr(const XmlNode& node) {
  const std::string& fullname = node.name();
  if (fullname.find("exp:")== std::string::npos) {
    casadi_error("DaeBuilder::read_expr: unknown - expression is supposed to "
                 "start with 'exp:' , got " + fullname);
  }

  // Chop the 'exp:'
  std::string name = fullname.substr(4);

  // The switch below is alphabetical, and can be thus made more efficient,
  // for example by using a switch statement of the first three letters,
  // if it would ever become a bottleneck
  if (name=="Add") {
    return read_expr(node[0]) + read_expr(node[1]);
  } else if (name=="Acos") {
    return acos(read_expr(node[0]));
  } else if (name=="Asin") {
    return asin(read_expr(node[0]));
  } else if (name=="Atan") {
    return atan(read_expr(node[0]));
  } else if (name=="Cos") {
    return cos(read_expr(node[0]));
  } else if (name=="Der") {
    return variables_.at(read_variable(node[0]).derivative).v;
  } else if (name=="Div") {
    return read_expr(node[0]) / read_expr(node[1]);
  } else if (name=="Exp") {
    return exp(read_expr(node[0]));
  } else if (name=="Identifier") {
    return read_variable(node).v;
  } else if (name=="IntegerLiteral") {
    casadi_int val;
    node.getText(val);
    return val;
  } else if (name=="Instant") {
    double val;
    node.getText(val);
    return val;
  } else if (name=="Log") {
    return log(read_expr(node[0]));
  } else if (name=="LogLeq") { // Logical less than equal
    return read_expr(node[0]) <= read_expr(node[1]);
  } else if (name=="LogGeq") { // Logical greater than equal
    return read_expr(node[0]) >= read_expr(node[1]);
  } else if (name=="LogLt") { // Logical less than
    return read_expr(node[0]) < read_expr(node[1]);
  } else if (name=="LogGt") { // Logical greater than
    return read_expr(node[0]) > read_expr(node[1]);
  } else if (name=="Max") {
    return fmax(read_expr(node[0]), read_expr(node[1]));
  } else if (name=="Min") {
    return fmin(read_expr(node[0]), read_expr(node[1]));
  } else if (name=="Mul") { // Multiplication
    return read_expr(node[0]) * read_expr(node[1]);
  } else if (name=="Neg") {
    return -read_expr(node[0]);
  } else if (name=="NoEvent") {
    // NOTE: This is a workaround, we assume that whenever NoEvent occurs,
    // what is meant is a switch
    casadi_int n = node.size();

    // Default-expression
    MX ex = read_expr(node[n-1]);

    // Evaluate ifs
    for (casadi_int i=n-3; i>=0; i -= 2) {
      ex = if_else(read_expr(node[i]), read_expr(node[i+1]), ex);
    }

    return ex;
  } else if (name=="Pow") {
    return pow(read_expr(node[0]), read_expr(node[1]));
  } else if (name=="RealLiteral") {
    double val;
    node.getText(val);
    return val;
  } else if (name=="Sin") {
    return sin(read_expr(node[0]));
  } else if (name=="Sqrt") {
    return sqrt(read_expr(node[0]));
  } else if (name=="StringLiteral") {
    throw CasadiException(node.getText());
  } else if (name=="Sub") {
    return read_expr(node[0]) - read_expr(node[1]);
  } else if (name=="Tan") {
    return tan(read_expr(node[0]));
  } else if (name=="Time") {
    return t;
  } else if (name=="TimedVariable") {
    return read_variable(node[0]).v;
  }

  // throw error if reached this point
  throw CasadiException(std::string("DaeBuilder::read_expr: Unknown node: ") + name);

}

void DaeBuilder::disp(std::ostream& stream, bool more) const {
  // Assert correctness
  if (more) sanity_check();

  // Print dimensions
  stream << "nx = " << this->x.size() << ", "
         << "nz = " << this->z.size() << ", "
         << "nq = " << this->q.size() << ", "
         << "ny = " << this->y.size() << ", "
         << "np = " << this->p.size() << ", "
         << "nc = " << this->c.size() << ", "
         << "nd = " << this->d.size() << ", "
         << "nw = " << this->w.size() << ", "
         << "nu = " << this->u.size();

  // Quick return?
  if (!more) return;
  stream << std::endl;

  // Print the functions
  if (!fun_.empty()) {
    stream << "Functions" << std::endl;
    for (const Function& f : fun_) {
      stream << "  " << f << std::endl;
    }
  }

  // Print the variables
  stream << "Variables" << std::endl;
  stream << "  t = " << str(this->t) << std::endl;
  if (!this->c.empty()) stream << "  c = " << str(this->c) << std::endl;
  if (!this->p.empty()) stream << "  p = " << str(this->p) << std::endl;
  if (!this->d.empty()) stream << "  d = " << str(this->d) << std::endl;
  if (!this->x.empty()) stream << "  x = " << str(this->x) << std::endl;
  if (!this->z.empty()) stream << "  z = " << str(this->z) << std::endl;
  if (!this->q.empty()) stream << "  q = " << str(this->q) << std::endl;
  if (!this->y.empty()) stream << "  y = " << str(this->y) << std::endl;
  if (!this->w.empty()) stream << "  w = " << str(this->w) << std::endl;
  if (!this->u.empty()) stream << "  u = " << str(this->u) << std::endl;

  if (!this->c.empty()) {
    stream << "Constants" << std::endl;
    for (casadi_int i=0; i<this->c.size(); ++i)
      stream << "  " << str(this->c[i]) << " == " << str(this->cdef[i]) << std::endl;
  }

  if (!this->d.empty()) {
    stream << "Dependent parameters" << std::endl;
    for (casadi_int i=0; i<this->d.size(); ++i)
      stream << "  " << str(this->d[i]) << " == " << str(this->ddef[i]) << std::endl;
  }

  if (!this->w.empty()) {
    stream << "Dependent variables" << std::endl;
    for (casadi_int i=0; i<this->w.size(); ++i)
      stream << "  " << str(this->w[i]) << " == " << str(this->wdef[i]) << std::endl;
  }

  if (!this->x.empty()) {
    stream << "Differential equations" << std::endl;
    for (casadi_int k=0; k<this->x.size(); ++k) {
      stream << "  " << str(der(this->x[k])) << " == " << str(this->ode[k]) << std::endl;
    }
  }

  if (!this->alg.empty()) {
    stream << "Algebraic equations" << std::endl;
    for (casadi_int k=0; k<this->z.size(); ++k) {
      stream << "  0 == " << str(this->alg[k]) << std::endl;
    }
  }

  if (!this->q.empty()) {
    stream << "Quadrature equations" << std::endl;
    for (casadi_int k=0; k<this->q.size(); ++k) {
      stream << "  " << str(der(this->q[k])) << " == " << str(this->quad[k]) << std::endl;
    }
  }

  if (!this->init_lhs.empty()) {
    stream << "Initial equations" << std::endl;
    for (casadi_int k=0; k<this->init_lhs.size(); ++k) {
      stream << "  " << str(this->init_lhs.at(k)) << " == " << str(this->init_rhs.at(k))
        << std::endl;
    }
  }

  if (!this->y.empty()) {
    stream << "Output variables" << std::endl;
    for (casadi_int i=0; i<this->y.size(); ++i) {
      stream << "  " << str(this->y[i]) << " == " << str(this->ydef[i]) << std::endl;
    }
  }
}

void DaeBuilder::eliminate_quad() {
  // Move all the quadratures to the list of differential states
  this->x.insert(this->x.end(), this->q.begin(), this->q.end());
  this->q.clear();
}

void DaeBuilder::sort_d() {
  sort_dependent(this->d, this->ddef);
}

void DaeBuilder::sort_w() {
  sort_dependent(this->w, this->wdef);
}

void DaeBuilder::sort_z(const std::vector<std::string>& z_order) {
  // Make sure lengths agree
  casadi_assert(z_order.size() == this->z.size(), "Dimension mismatch");
  // Mark existing components in z
  std::vector<bool> old_z(variables_.size(), false);
  for (size_t i = 0; i < this->z.size(); ++i) {
    std::string s = this->z.at(i).name();
    auto it = varind_.find(s);
    casadi_assert(it != varind_.end(), "No such variable: \"" + s + "\".");
    old_z.at(it->second) = true;
  }
  // New vector of z
  std::vector<MX> new_z;
  new_z.reserve(z_order.size());
  for (const std::string& s : z_order) {
    auto it = varind_.find(s);
    casadi_assert(it != varind_.end(), "No such variable: \"" + s + "\".");
    casadi_assert(old_z.at(it->second), "Variable \"" + s + "\" is not an algebraic variable.");
    new_z.push_back(variables_.at(it->second).v);
  }
  // Success: Update z
  std::copy(new_z.begin(), new_z.end(), this->z.begin());
}

void DaeBuilder::prune(bool prune_p, bool prune_u) {
  // Function inputs and outputs
  std::vector<MX> f_in, f_out, v;
  std::vector<std::string> f_in_name, f_out_name;
  // Collect all DAE input variables with at least one entry, skip u
  for (casadi_int i = 0; i != DAE_BUILDER_NUM_IN; ++i) {
    if (prune_p && i == DAE_BUILDER_P) continue;
    if (prune_u && i == DAE_BUILDER_U) continue;
    v = input(static_cast<DaeBuilderIn>(i));
    if (!v.empty()) {
      f_in.push_back(vertcat(v));
      f_in_name.push_back(to_string(static_cast<DaeBuilderIn>(i)));
    }
  }
  // Collect all DAE output variables with at least one entry
  for (casadi_int i = 0; i != DAE_BUILDER_NUM_OUT; ++i) {
    v = output(static_cast<DaeBuilderOut>(i));
    if (!v.empty()) {
      f_out.push_back(vertcat(v));
      f_out_name.push_back(to_string(static_cast<DaeBuilderOut>(i)));
    }
  }
  // Create a function
  Function f("prune_fcn", f_in, f_out, f_in_name, f_out_name);
  // Mark which variables are free
  std::vector<bool> free_variables(variables_.size(), false);
  for (const std::string& s : f.get_free()) {
    auto it = varind_.find(s);
    casadi_assert(it != varind_.end(), "No such variable: \"" + s + "\".");
    free_variables.at(it->second) = true;
  }
  // Prune p
  if (prune_p) {
    size_t np = 0;
    for (size_t i = 0; i < this->p.size(); ++i) {
      std::string s = this->p.at(i).name();
      auto it = varind_.find(s);
      casadi_assert(it != varind_.end(), "No such variable: \"" + s + "\".");
      if (!free_variables.at(it->second)) this->p.at(np++) = this->p.at(i);
    }
    this->p.resize(np);
  }
  // Prune u
  if (prune_u) {
    size_t nu = 0;
    for (size_t i = 0; i < this->u.size(); ++i) {
      std::string s = this->u.at(i).name();
      auto it = varind_.find(s);
      casadi_assert(it != varind_.end(), "No such variable: \"" + s + "\".");
      if (!free_variables.at(it->second)) this->u.at(nu++) = this->u.at(i);
    }
    this->u.resize(nu);
  }
}

const Variable& DaeBuilder::variable(const std::string& name) const {
  return const_cast<DaeBuilder*>(this)->variable(name);
}

Variable& DaeBuilder::variable(const std::string& name) {
  // Find the variable
  auto it = varind_.find(name);
  if (it == varind_.end()) casadi_error("No such variable: \"" + name + "\".");

  // Return the variable
  return variables_.at(it->second);
}

bool DaeBuilder::has_variable(const std::string& name) const {
  return varind_.find(name) != varind_.end();
}

void DaeBuilder::add_variable(const std::string& name, const Variable& var) {
  // Try to find the component
  casadi_assert(!has_variable(name), "Variable \"" + name + "\" has already been added.");
  // Add to the map of all variables
  varind_[name] = variables_.size();
  variables_.push_back(var);
  // Clear cache
  clear_cache_ = true;
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

MX DaeBuilder::add_x(const std::string& name, casadi_int n) {
  if (name.empty()) return add_x("x" + str(this->x.size()), n);
  MX new_x = add_variable(name, n);
  this->x.push_back(new_x);
  return new_x;
}

void DaeBuilder::register_p(const MX& new_p) {
  // Consistency checks
  casadi_assert(has_variable(new_p.name()), "No such variable: " + new_p.name());
  // Add to list
  this->p.push_back(new_p);
}

void DaeBuilder::register_u(const MX& new_u) {
  // Consistency checks
  casadi_assert(has_variable(new_u.name()), "No such variable: " + new_u.name());
  // Add to list
  this->u.push_back(new_u);
}

void DaeBuilder::register_x(const MX& new_x) {
  // Consistency checks
  casadi_assert(has_variable(new_x.name()), "No such variable: " + new_x.name());
  // Add to list
  this->x.push_back(new_x);
}

void DaeBuilder::register_z(const MX& new_z) {
  // Consistency checks
  casadi_assert(has_variable(new_z.name()), "No such variable: " + new_z.name());
  // Add to list
  this->z.push_back(new_z);
}

void DaeBuilder::register_t(const MX& new_t) {
  // Save to class
  this->t = new_t;
}

void DaeBuilder::register_c(const MX& new_c, const MX& new_cdef) {
  // Consistency check
  casadi_assert(new_c.sparsity() == new_cdef.sparsity(), "Mismatching sparsity");
  casadi_assert(has_variable(new_c.name()), "No such variable: " + new_c.name());
  // Add to lists
  this->c.push_back(new_c);
  this->cdef.push_back(new_cdef);
}

void DaeBuilder::register_d(const MX& new_d, const MX& new_ddef) {
  // Consistency checks
  casadi_assert(new_d.sparsity() == new_ddef.sparsity(), "Mismatching sparsity");
  casadi_assert(has_variable(new_d.name()), "No such variable: " + new_d.name());
  // Add to lists
  this->d.push_back(new_d);
  this->ddef.push_back(new_ddef);
  this->lam_ddef.push_back(MX::sym("lam_" + new_d.name(), new_d.sparsity()));
}

void DaeBuilder::register_w(const MX& new_w, const MX& new_wdef) {
  // Consistency checks
  casadi_assert(new_w.sparsity() == new_wdef.sparsity(), "Mismatching sparsity");
  casadi_assert(has_variable(new_w.name()), "No such variable: " + new_w.name());
  // Add to lists
  this->w.push_back(new_w);
  this->wdef.push_back(new_wdef);
  this->lam_wdef.push_back(MX::sym("lam_" + new_w.name(), new_w.sparsity()));
}

void DaeBuilder::register_y(const MX& new_y, const MX& new_ydef) {
  // Consistency checks
  casadi_assert(new_y.sparsity() == new_ydef.sparsity(), "Mismatching sparsity");
  casadi_assert(has_variable(new_y.name()), "No such variable: " + new_y.name());
  // Add to lists
  this->y.push_back(new_y);
  this->ydef.push_back(new_ydef);
  this->lam_ydef.push_back(MX::sym("lam_" + new_y.name(), new_y.sparsity()));
}

MX DaeBuilder::add_q(const std::string& name, casadi_int n) {
  if (name.empty()) return add_q("q" + str(this->q.size()), n);
  MX new_q = add_variable(name, n);
  this->q.push_back(new_q);
  return new_q;
}

MX DaeBuilder::add_z(const std::string& name, casadi_int n) {
  if (name.empty()) return add_z("z" + str(this->z.size()), n);
  MX new_z = add_variable(name, n);
  this->z.push_back(new_z);
  return new_z;
}

MX DaeBuilder::add_p(const std::string& name, casadi_int n) {
  if (name.empty()) return add_p("p" + str(this->p.size()), n);
  MX new_p = add_variable(name, n);
  this->p.push_back(new_p);
  return new_p;
}

MX DaeBuilder::add_u(const std::string& name, casadi_int n) {
  if (name.empty()) return add_u("u" + str(this->u.size()), n);
  MX new_u = add_variable(name, n);
  this->u.push_back(new_u);
  return new_u;
}

MX DaeBuilder::add_aux(const std::string& name, casadi_int n) {
  if (name.empty()) return add_aux("aux" + str(this->aux.size()), n);
  MX new_aux = add_variable(name, n);
  this->aux.push_back(new_aux);
  return new_aux;
}

void DaeBuilder::add_init(const MX& lhs, const MX& rhs) {
  this->init_lhs.push_back(lhs);
  this->init_rhs.push_back(rhs);
}

MX DaeBuilder::add_d(const std::string& name, const MX& new_ddef) {
  MX new_d = add_variable(name, new_ddef.sparsity());
  this->d.push_back(new_d);
  this->ddef.push_back(new_ddef);
  this->lam_ddef.push_back(MX::sym("lam_" + name, new_ddef.sparsity()));
  return new_d;
}

MX DaeBuilder::add_w(const std::string& name, const MX& new_wdef) {
  MX new_w = add_variable(name, new_wdef.sparsity());
  this->w.push_back(new_w);
  this->wdef.push_back(new_wdef);
  this->lam_wdef.push_back(MX::sym("lam_" + name, new_wdef.sparsity()));
  return new_w;
}

MX DaeBuilder::add_y(const std::string& name, const MX& new_ydef) {
  MX new_y = add_variable(name, new_ydef.sparsity());
  this->y.push_back(new_y);
  this->ydef.push_back(new_ydef);
  this->lam_ydef.push_back(MX::sym("lam_" + name, new_ydef.sparsity()));
  return new_y;
}

void DaeBuilder::add_ode(const std::string& name, const MX& new_ode) {
  this->ode.push_back(new_ode);
  this->lam_ode.push_back(MX::sym("lam_" + name, new_ode.sparsity()));
  clear_cache_ = true;
}

void DaeBuilder::add_alg(const std::string& name, const MX& new_alg) {
  this->alg.push_back(new_alg);
  this->lam_alg.push_back(MX::sym("lam_" + name, new_alg.sparsity()));
  clear_cache_ = true;
}

void DaeBuilder::add_quad(const std::string& name, const MX& new_quad) {
  this->quad.push_back(new_quad);
  this->lam_quad.push_back(MX::sym("lam_" + name, new_quad.sparsity()));
  clear_cache_ = true;
}

void DaeBuilder::sanity_check() const {
  // Time
  casadi_assert(this->t.is_symbolic(), "Non-symbolic time t");
  casadi_assert(this->t.is_scalar(), "Non-scalar time t");

  // Differential states
  casadi_assert(this->x.size()==this->ode.size(),
                        "x and ode have different lengths");
  for (casadi_int i=0; i<this->x.size(); ++i) {
    casadi_assert(this->x[i].size()==this->ode[i].size(),
                          "ode has wrong dimensions");
    casadi_assert(this->x[i].is_symbolic(), "Non-symbolic state x");
  }

  // Algebraic variables/equations
  casadi_assert(this->z.size()==this->alg.size(),
                        "z and alg have different lengths");
  for (casadi_int i=0; i<this->z.size(); ++i) {
    casadi_assert(this->z[i].is_symbolic(), "Non-symbolic algebraic variable z");
    casadi_assert(this->z[i].size()==this->alg[i].size(),
                          "alg has wrong dimensions");
  }

  // Quadrature states/equations
  casadi_assert(this->q.size()==this->quad.size(), "q and quad have different lengths");
  for (casadi_int i=0; i<this->q.size(); ++i) {
    casadi_assert(this->q[i].is_symbolic(), "Non-symbolic quadrature state q");
    casadi_assert(this->q[i].size()==this->quad[i].size(),
                          "quad has wrong dimensions");
  }

  // Dependent parameters
  casadi_assert(this->d.size()==this->ddef.size(), "d and ddef have different lengths");
  for (casadi_int i=0; i<this->d.size(); ++i) {
    casadi_assert(this->d[i].is_symbolic(), "Non-symbolic dependent parameter d");
    casadi_assert(this->d[i].size()==this->ddef[i].size(), "ddef has wrong dimensions");
  }

  // Dependent variables
  casadi_assert(this->w.size()==this->wdef.size(), "w and wdef have different lengths");
  for (casadi_int i=0; i<this->w.size(); ++i) {
    casadi_assert(this->w[i].is_symbolic(), "Non-symbolic dependent parameter v");
    casadi_assert(this->w[i].size()==this->wdef[i].size(), "wdef has wrong dimensions");
  }

  // Output equations
  casadi_assert(this->y.size()==this->ydef.size(), "y and ydef have different lengths");
  for (casadi_int i=0; i<this->y.size(); ++i) {
    casadi_assert(this->y[i].is_symbolic(), "Non-symbolic output y");
    casadi_assert(this->y[i].size()==this->ydef[i].size(), "ydef has wrong dimensions");
  }

  // Control
  for (casadi_int i=0; i<this->u.size(); ++i) {
    casadi_assert(this->u[i].is_symbolic(), "Non-symbolic control u");
  }

  // Parameter
  for (casadi_int i=0; i<this->p.size(); ++i) {
    casadi_assert(this->p[i].is_symbolic(), "Non-symbolic parameter p");
  }

  // Initial equations
  casadi_assert(this->init_lhs.size() == this->init_rhs.size(),
    "init_lhs and init_rhs have different lengths");
}

std::string DaeBuilder::qualified_name(const XmlNode& nn) {
  // Stringstream to assemble name
  std::stringstream qn;

  for (casadi_int i=0; i<nn.size(); ++i) {
    // Add a dot
    if (i!=0) qn << ".";

    // Get the name part
    qn << nn[i].attribute<std::string>("name");

    // Get the index, if any
    if (nn[i].size()>0) {
      casadi_int ind;
      nn[i]["exp:ArraySubscripts"]["exp:IndexExpression"]["exp:IntegerLiteral"].getText(ind);
      qn << "[" << ind << "]";
    }
  }

  // Return the name
  return qn.str();
}

MX DaeBuilder::var(const std::string& name) const {
  return variable(name).v;
}

MX DaeBuilder::der(const std::string& name) const {
  return variables_.at(variable(name).derivative).v;
}

MX DaeBuilder::der(const MX& var) const {
  casadi_assert_dev(var.is_column() && var.is_symbolic());
  return der(var.name());
}

void DaeBuilder::eliminate_w() {
  // Quick return if no w
  if (this->w.empty()) return;
  // Ensure variables are sorted
  sort_w();
  // Expressions where the variables are also being used
  std::vector<MX> ex;
  ex.insert(ex.end(), this->alg.begin(), this->alg.end());
  ex.insert(ex.end(), this->ode.begin(), this->ode.end());
  ex.insert(ex.end(), this->quad.begin(), this->quad.end());
  ex.insert(ex.end(), this->ydef.begin(), this->ydef.end());
  // Perform elimination
  substitute_inplace(this->w, this->wdef, ex);
  // Clear lists
  this->w.clear();
  this->wdef.clear();
  // Get algebraic equations
  auto it = ex.begin();
  std::copy(it, it + this->alg.size(), this->alg.begin());
  it += this->alg.size();
  // Get differential equations
  std::copy(it, it + this->ode.size(), this->ode.begin());
  it += this->ode.size();
  // Get quadrature equations
  std::copy(it, it + this->quad.size(), this->quad.begin());
  it += this->quad.size();
  // Get output equations
  std::copy(it, it + this->ydef.size(), this->ydef.begin());
  it += this->ydef.size();
  // Consistency check
  casadi_assert_dev(it == ex.end());
}

void DaeBuilder::lift(bool lift_shared, bool lift_calls) {
  // Partially implemented
  if (x.size() > 0) casadi_warning("Only lifting algebraic variables");
  // Lift algebraic expressions
  std::vector<MX> new_w, new_wdef;
  Dict opts{{"lift_shared", lift_shared}, {"lift_calls", lift_calls},
    {"prefix", "w_"}, {"suffix", ""}, {"offset", static_cast<casadi_int>(this->w.size())}};
  extract(this->alg, new_w, new_wdef, opts);
  // Register as dependent variables
  for (size_t i = 0; i < new_w.size(); ++i) {
    add_variable(new_w.at(i));
    register_w(new_w.at(i), new_wdef.at(i));
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

std::string DaeBuilder::unit(const MX& var) const {
  casadi_assert(!var.is_column() && var.is_valid_input(),
                        "DaeBuilder::unit: Argument must be a symbolic vector");
  if (var.is_empty()) {
    return "n/a";
  } else {
    std::vector<MX> prim = var.primitives();
    std::string ret = unit(prim.at(0).name());
    for (casadi_int i=1; i<prim.size(); ++i) {
      casadi_assert(ret == unit(prim.at(i).name()),
                            "DaeBuilder::unit: Argument has mixed units");
    }
    return ret;
  }
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

std::string to_string(DaeBuilder::DaeBuilderIn v) {
  switch (v) {
  case DaeBuilder::DAE_BUILDER_T: return "t";
  case DaeBuilder::DAE_BUILDER_C: return "c";
  case DaeBuilder::DAE_BUILDER_P: return "p";
  case DaeBuilder::DAE_BUILDER_D: return "d";
  case DaeBuilder::DAE_BUILDER_W: return "w";
  case DaeBuilder::DAE_BUILDER_U: return "u";
  case DaeBuilder::DAE_BUILDER_X: return "x";
  case DaeBuilder::DAE_BUILDER_Z: return "z";
  case DaeBuilder::DAE_BUILDER_Q: return "q";
  case DaeBuilder::DAE_BUILDER_Y: return "y";
  default: break;
  }
  return "";
}

std::string to_string(DaeBuilder::DaeBuilderOut v) {
  switch (v) {
  case DaeBuilder::DAE_BUILDER_DDEF: return "ddef";
  case DaeBuilder::DAE_BUILDER_WDEF: return "wdef";
  case DaeBuilder::DAE_BUILDER_ODE: return "ode";
  case DaeBuilder::DAE_BUILDER_ALG: return "alg";
  case DaeBuilder::DAE_BUILDER_QUAD: return "quad";
  case DaeBuilder::DAE_BUILDER_YDEF: return "ydef";
  default: break;
  }
  return "";
}

std::vector<MX> DaeBuilder::input(DaeBuilderIn ind) const {
  switch (ind) {
  case DAE_BUILDER_T: return std::vector<MX>(1, this->t);
  case DAE_BUILDER_C: return this->c;
  case DAE_BUILDER_P: return this->p;
  case DAE_BUILDER_D: return this->d;
  case DAE_BUILDER_W: return this->w;
  case DAE_BUILDER_U: return this->u;
  case DAE_BUILDER_X: return this->x;
  case DAE_BUILDER_Z: return this->z;
  case DAE_BUILDER_Q: return this->q;
  case DAE_BUILDER_Y: return this->y;
  default: return std::vector<MX>();
  }
}

std::vector<MX> DaeBuilder::input(const std::vector<DaeBuilderIn>& ind) const {
  std::vector<MX> ret(ind.size());
  for (casadi_int i=0; i<ind.size(); ++i) {
    ret[i] = vertcat(input(ind[i]));
  }
  return ret;
}

std::vector<MX> DaeBuilder::output(DaeBuilderOut ind) const {
  switch (ind) {
  case DAE_BUILDER_DDEF: return this->ddef;
  case DAE_BUILDER_WDEF: return this->wdef;
  case DAE_BUILDER_ODE: return this->ode;
  case DAE_BUILDER_ALG: return this->alg;
  case DAE_BUILDER_QUAD: return this->quad;
  case DAE_BUILDER_YDEF: return this->ydef;
  default: return std::vector<MX>();
  }
}

std::vector<MX> DaeBuilder::output(const std::vector<DaeBuilderOut>& ind) const {
  std::vector<MX> ret(ind.size());
  for (casadi_int i=0; i<ind.size(); ++i) {
    ret[i] = vertcat(output(ind[i]));
  }
  return ret;
}

void DaeBuilder::add_lc(const std::string& name,
                      const std::vector<std::string>& f_out) {
  // Make sure object valid
  sanity_check();

  // Make sure name is valid
  casadi_assert(!name.empty(), "DaeBuilder::add_lc: \"name\" is empty");
  for (std::string::const_iterator i=name.begin(); i!=name.end(); ++i) {
    casadi_assert(isalnum(*i),
                          "DaeBuilder::add_lc: \"name\" must be alphanumeric");
  }

  // Consistency checks
  casadi_assert(!f_out.empty(), "DaeBuilder::add_lc: Linear combination is empty");
  std::vector<bool> in_use(DAE_BUILDER_NUM_OUT, false);
  for (casadi_int i=0; i<f_out.size(); ++i) {
    DaeBuilderOut oind = to_enum<DaeBuilderOut>(f_out[i]);
    casadi_assert(!in_use[oind], "DaeBuilder::add_lc: Duplicate expression " + f_out[i]);
    in_use[oind] = true;
  }

  std::vector<std::string>& ret1 = lc_[name];
  if (!ret1.empty()) casadi_warning("DaeBuilder::add_lc: Overwriting " << name);
  ret1 = f_out;
}

Function DaeBuilder::create(const std::string& fname,
    const std::vector<std::string>& s_in,
    const std::vector<std::string>& s_out, bool sx, bool lifted_calls) const {
  // Are there any '_' in the names?
  bool with_underscore = false;
  for (auto s_io : {&s_in, &s_out}) {
    for (const std::string& s : *s_io) {
      with_underscore = with_underscore || std::count(s.begin(), s.end(), '_');
    }
  }
  // Replace '_' with ':', if needed
  if (with_underscore) {
    std::vector<std::string> s_in_mod(s_in), s_out_mod(s_out);
    for (auto s_io : {&s_in_mod, &s_out_mod}) {
      for (std::string& s : *s_io) replace(s.begin(), s.end(), '_', ':');
    }
    // Recursive call
    return create(fname, s_in_mod, s_out_mod, sx, lifted_calls);
  }
  // Check if dependent variables are given and needed
  bool elim_w = false;
  if (!this->w.empty()) {
    // Dependent variables exists, eliminate unless v is given
    elim_w = true;
    for (const std::string& s : s_in) {
      if (s == "w") {
        // Dependent variables are given
        elim_w = false;
        break;
      }
    }
  }
  // Are lifted calls really needed?
  if (lifted_calls) {
    // Consistency check
    casadi_assert(!elim_w, "Lifted calls cannot be used if dependent variables are eliminated");
    // Only lift calls if really needed
    lifted_calls = false;
    for (const MX& vdef_comp : this->wdef) {
      if (vdef_comp.is_output()) {
        // There are indeed function calls present
        lifted_calls = true;
        break;
      }
    }
  }
  // Call factory without lifted calls
  std::string fname_nocalls = lifted_calls ? fname + "_nocalls" : fname;
  Function ret = oracle(sx, elim_w, lifted_calls).factory(fname_nocalls, s_in, s_out, lc_);
  // If no lifted calls, done
  if (!lifted_calls) return ret;
  // MX expressions for ret without lifted calls
  std::vector<MX> ret_in = ret.mx_in();
  std::vector<MX> ret_out = ret(ret_in);
  // Offsets in v
  std::vector<casadi_int> h_offsets = offset(this->w);
  // Split "w", "lam_wdef" into components
  std::vector<MX> v_in, lam_vdef_in;
  for (size_t i = 0; i < s_in.size(); ++i) {
    if (ret.name_in(i) == "w") {
      v_in = vertsplit(ret_in[i], h_offsets);
    } else if (ret.name_in(i) == "lam_wdef") {
      lam_vdef_in = vertsplit(ret_in[i], h_offsets);
    }
  }
  // Map dependent variables into index in vector
  std::map<MXNode*, size_t> v_map;
  for (size_t i = 0; i < this->w.size(); ++i) {
    v_map[this->w.at(i).get()] = i;
  }
  // Collect all the call nodes
  std::map<MXNode*, CallIO> call_nodes;
  for (size_t vdefind = 0; vdefind < this->wdef.size(); ++vdefind) {
    // Current element handled
    const MX& vdefref = this->wdef.at(vdefind);
    // Handle function call nodes
    if (vdefref.is_output()) {
      // Get function call node
      MX c = vdefref.dep(0);
      // Find the corresponding call node in the map
      auto call_it = call_nodes.find(c.get());
      // If first time this call node is encountered
      if (call_it == call_nodes.end()) {
        // Create new CallIO struct
        CallIO cio;
        // Save function instance
        cio.f = c.which_function();
        // Expressions for function call inputs
        cio.v.resize(c.n_dep(), -1);
        cio.arg.resize(cio.v.size());
        for (casadi_int i = 0; i < cio.v.size(); ++i) {
          if (c.dep(i).is_constant()) {
            cio.arg.at(i) = c.dep(i);
          } else {
            size_t v_ind = v_map.at(c.dep(i).get());
            cio.v.at(i) = v_ind;
            cio.arg.at(i) = v_in.at(v_ind);
          }
        }
        // Allocate memory for function call outputs
        cio.vdef.resize(c.n_out(), -1);
        cio.res.resize(cio.vdef.size());
        // Allocate memory for adjoint seeds, if any
        if (!lam_vdef_in.empty()) cio.adj1_arg.resize(c.n_out());
        // Save to map and update iterator
        call_it = call_nodes.insert(std::make_pair(c.get(), cio)).first;
      }
      // Which output of the function are we calculating?
      casadi_int oind = vdefref.which_output();
      // Save output expression to structure
      call_it->second.vdef.at(oind) = vdefind;
      call_it->second.res.at(oind) = v_in.at(vdefind);
      // Save adjoint seed to structure, if any
      if (!lam_vdef_in.empty()) call_it->second.adj1_arg.at(oind) = lam_vdef_in.at(vdefind);
    }
  }
  // Additional term in jac_vdef_v
  for (size_t i = 0; i < ret_out.size(); ++i) {
    if (ret.name_out(i) == "jac_wdef_w") {
      ret_out.at(i) += jac_vdef_v_from_calls(call_nodes, h_offsets);
    }
  }
  // Additional term in hess_?_v_v where ? is any linear combination containing vdef
  MX extra_hess_v_v;  // same for all linear combinations, if multiple
  for (auto&& e : lc_) {
    // Find out of vdef is part of the linear combination
    bool has_vdef = false;
    for (const std::string& r : e.second) {
      if (r == "wdef") {
        has_vdef = true;
        break;
      }
    }
    // Skip if linear combination does not depend on vdef
    if (!has_vdef) continue;
    // Search for matching function outputs
    for (size_t i = 0; i < ret_out.size(); ++i) {
      if (ret.name_out(i) == "hess_" + e.first + "_w_w") {
        // Calculate contribution to hess_?_v_v
        if (extra_hess_v_v.is_empty())
          extra_hess_v_v = hess_v_v_from_calls(call_nodes, h_offsets);
        // Add contribution to output
        ret_out.at(i) += extra_hess_v_v;
      }
    }
  }
  // Assemble modified return function and return
  ret = Function(fname, ret_in, ret_out, ret.name_in(), ret.name_out());
  return ret;
}

MX DaeBuilder::jac_vdef_v_from_calls(std::map<MXNode*, CallIO>& call_nodes,
    const std::vector<casadi_int>& h_offsets) const {
  // Calculate all Jacobian expressions
  for (auto call_it = call_nodes.begin(); call_it != call_nodes.end(); ++call_it) {
    call_it->second.calc_jac();
  }
  // Row offsets in jac_vdef_v
  casadi_int voffset_begin = 0, voffset_end = 0, voffset_last = 0;
  // Vertical and horizontal slices of jac_vdef_v
  std::vector<MX> vblocks, hblocks;
  // All blocks for this block row
  std::map<size_t, MX> jac_brow;
  // Collect all Jacobian blocks
  for (size_t vdefind = 0; vdefind < this->wdef.size(); ++vdefind) {
    // Current element handled
    const MX& vdefref = this->wdef.at(vdefind);
    // Update vertical offset
    voffset_begin = voffset_end;
    voffset_end += vdefref.numel();
    // Handle function call nodes
    if (vdefref.is_output()) {
      // Which output of the function are we calculating?
      casadi_int oind = vdefref.which_output();
      // Get function call node
      MX c = vdefref.dep(0);
      // Find data about inputs and outputs
      auto call_it = call_nodes.find(c.get());
      casadi_assert_dev(call_it != call_nodes.end());
      // Collect all blocks for this block row
      jac_brow.clear();
      for (casadi_int iind = 0; iind < call_it->second.arg.size(); ++iind) {
        size_t vind = call_it->second.v.at(iind);
        if (vind != size_t(-1)) {
          jac_brow[vind] = call_it->second.jac(oind, iind);
        }
      }
      // Add empty rows to vblocks, if any
      if (voffset_last != voffset_begin) {
        vblocks.push_back(MX(voffset_begin - voffset_last, h_offsets.back()));
      }
      // Collect horizontal blocks
      hblocks.clear();
      casadi_int hoffset = 0;
      for (auto e : jac_brow) {
        // Add empty block before Jacobian block, if needed
        if (hoffset < h_offsets.at(e.first))
          hblocks.push_back(MX(vdefref.numel(), h_offsets.at(e.first) - hoffset));
        // Add Jacobian block
        hblocks.push_back(e.second);
        // Update offsets
        hoffset = h_offsets.at(e.first + 1);
      }
      // Add trailing empty block, if needed
      if (hoffset < h_offsets.back())
        hblocks.push_back(MX(vdefref.numel(), h_offsets.back() - hoffset));
      // Add new block row to vblocks
      vblocks.push_back(horzcat(hblocks));
      // Keep track of the offset handled in jac_brow
      voffset_last = voffset_end;
    }
  }
  // Add empty trailing row to vblocks, if any
  if (voffset_last != voffset_end) {
    vblocks.push_back(MX(voffset_end - voffset_last, h_offsets.back()));
  }
  // Return additional term in jac_vdef_v
  return vertcat(vblocks);
}

MX DaeBuilder::hess_v_v_from_calls(std::map<MXNode*, CallIO>& call_nodes,
    const std::vector<casadi_int>& h_offsets) const {
  // Calculate all Hessian expressions
  for (auto&& call_ref : call_nodes) call_ref.second.calc_hess();
  // Row offsets in hess_v_v
  casadi_int voffset_begin = 0, voffset_end = 0, voffset_last = 0;
  // Vertical and horizontal slices of hess_v_v
  std::vector<MX> vblocks, hblocks;
  // All blocks for a block row
  std::map<size_t, MX> hess_brow;
  // Loop over block rows
  for (size_t vind1 = 0; vind1 < this->w.size(); ++vind1) {
    // Current element handled
    const MX& vref = this->w.at(vind1);
    // Update vertical offset
    voffset_begin = voffset_end;
    voffset_end += vref.numel();
    // Collect all blocks for this block row
    hess_brow.clear();
    for (auto&& call_ref : call_nodes) {
      // Locate the specific index
      for (size_t iind1 = 0; iind1 < call_ref.second.v.size(); ++iind1) {
        if (call_ref.second.v.at(iind1) == vind1) {
          // Add contribution to block row
          for (size_t iind2 = 0; iind2 < call_ref.second.v.size(); ++iind2) {
            // Corresponding index in v
            size_t vind2 = call_ref.second.v[iind2];
            if (vind2 == size_t(-1)) continue;
            // Hessian contribution
            MX H_contr = call_ref.second.hess(iind1, iind2);
            // Insert new block or add to existing one
            auto it = hess_brow.find(vind2);
            if (it != hess_brow.end()) {
              it->second += H_contr;
            } else {
              hess_brow[vind2] = H_contr;
            }
          }
          // An index can only appear once
          break;
        }
      }
    }
    // If no blocks, skip row
    if (hess_brow.empty()) continue;
    // Add empty rows to vblocks, if any
    if (voffset_last != voffset_begin) {
      vblocks.push_back(MX(voffset_begin - voffset_last, h_offsets.back()));
    }
    // Collect horizontal blocks
    hblocks.clear();
    casadi_int hoffset = 0;
    for (auto e : hess_brow) {
      // Add empty block before Jacobian block, if needed
      if (hoffset < h_offsets.at(e.first))
        hblocks.push_back(MX(vref.numel(), h_offsets.at(e.first) - hoffset));
      // Add Jacobian block
      hblocks.push_back(e.second);
      // Update offsets
      hoffset = h_offsets.at(e.first + 1);
    }
    // Add trailing empty block, if needed
    if (hoffset < h_offsets.back())
      hblocks.push_back(MX(vref.numel(), h_offsets.back() - hoffset));
    // Add new block row to vblocks
    vblocks.push_back(horzcat(hblocks));
    // Keep track of the offset handled in jac_brow
    voffset_last = voffset_end;
  }
  // Add empty trailing row to vblocks, if any
  if (voffset_last != voffset_end) {
    vblocks.push_back(MX(voffset_end - voffset_last, h_offsets.back()));
  }
  // Return additional term in jac_vdef_v
  return vertcat(vblocks);
}

Function DaeBuilder::add_fun(const Function& f) {
  casadi_assert(!has_fun(f.name()), "Function '" + f.name() + "' already exists");
  fun_.push_back(f);
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
    for (v_ind=0; v_ind<this->w.size(); ++v_ind) {
      if (s==this->w.at(v_ind).name()) {
        res_ex.push_back(this->wdef.at(v_ind));
        break;
      }
    }
    casadi_assert(v_ind<this->w.size(), "Cannot find dependent '" + s + "'");
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
  for (const Function& f : fun_) {
    if (f.name()==name) return true;
  }
  return false;
}

Function DaeBuilder::fun(const std::string& name) const {
  casadi_assert(has_fun(name), "No such function: '" + name + "'");
  for (const Function& f : fun_) {
    if (f.name()==name) return f;
  }
  return Function();
}

void DaeBuilder::clear_cache() const {
  for (bool sx : {false, true}) {
    for (bool elim_w : {false, true}) {
      for (bool lifted_calls : {false, true}) {
        Function& fref = oracle_[sx][elim_w][lifted_calls];
        if (!fref.is_null()) fref = Function();
      }
    }
  }
  clear_cache_ = false;
}

const Function& DaeBuilder::oracle(bool sx, bool elim_w, bool lifted_calls) const {
  // Clear cache now, if necessary
  if (clear_cache_) clear_cache();
  // Create an MX oracle, if needed
  if (oracle_[false][elim_w][lifted_calls].is_null()) {
    // Oracle function inputs and outputs
    std::vector<MX> f_in, f_out, v;
    std::vector<std::string> f_in_name, f_out_name;
    // Index for vdef
    casadi_int vdef_ind = -1;
    // Options consistency check
    casadi_assert(!(elim_w && lifted_calls), "Incompatible options");
    // Do we need to substitute out v
    bool subst_v = false;
    // Collect all DAE input variables with at least one entry
    for (casadi_int i = 0; i != DAE_BUILDER_NUM_IN; ++i) {
      v = input(static_cast<DaeBuilderIn>(i));
      if (!v.empty()) {
        if (elim_w && i == DAE_BUILDER_W) {
          subst_v = true;
        } else {
          f_in.push_back(vertcat(v));
          f_in_name.push_back(to_string(static_cast<DaeBuilderIn>(i)));
        }
      }
    }
    // Collect all DAE output variables with at least one entry
    for (casadi_int i = 0; i != DAE_BUILDER_NUM_OUT; ++i) {
      v = output(static_cast<DaeBuilderOut>(i));
      if (!v.empty()) {
        if (i == DAE_BUILDER_WDEF) vdef_ind = f_out.size();
        f_out.push_back(vertcat(v));
        f_out_name.push_back(to_string(static_cast<DaeBuilderOut>(i)));
      }
    }
    // Eliminate v from inputs
    if (subst_v) {
      // Make a copy of dependent variable definitions to avoid modifying member variable
      std::vector<MX> vdef(this->wdef);
      // Perform in-place substitution
      substitute_inplace(this->w, vdef, f_out, false);
    } else if (lifted_calls && vdef_ind >= 0) {
      // Make a copy of dependent variable definitions to avoid modifying member variable
      std::vector<MX> vdef(this->wdef);
      // Remove references to call nodes
      for (MX& vdefref : vdef) {
        if (vdefref.is_output()) vdefref = MX::zeros(vdefref.sparsity());
      }
      // Save to oracle outputs
      f_out.at(vdef_ind) = vertcat(vdef);
    }
    // Create oracle
    oracle_[false][elim_w][lifted_calls]
      = Function("mx_oracle", f_in, f_out, f_in_name, f_out_name);
  }
  // Return MX oracle, if requested
  if (!sx) return oracle_[false][elim_w][lifted_calls];
  // Create SX oracle, if needed
  Function& sx_oracle = oracle_[true][elim_w][lifted_calls];
  if (sx_oracle.is_null()) sx_oracle = oracle_[false][elim_w][lifted_calls].expand("sx_oracle");
  // Return SX oracle reference
  return sx_oracle;
}

void DaeBuilder::CallIO::calc_jac() {
  // Consistency checks
  for (casadi_int i = 0; i < this->f.n_in(); ++i) {
    casadi_assert(this->f.size_in(i) == this->arg.at(i).size(), "Call input not provided");
  }
  for (casadi_int i = 0; i < this->f.n_out(); ++i) {
    casadi_assert(this->f.size_out(i) == this->res.at(i).size(), "Call output not provided");
  }
  // Get/generate the (cached) Jacobian function
  // casadi_message("Retrieving the Jacobian of " + str(this->f));
  this->J = this->f.jacobian();
  // casadi_message("Retrieving the Jacobian of " + str(this->f) + " done");
  // Input expressions for the call to J
  std::vector<MX> call_in = this->arg;
  call_in.insert(call_in.end(), this->res.begin(), this->res.end());
  // Create expressions for Jacobian blocks and save to struct
  this->jac_res = this->J(call_in);
}

void DaeBuilder::CallIO::calc_grad() {
  // Consistency checks
  for (casadi_int i = 0; i < this->f.n_in(); ++i) {
    casadi_assert(this->f.size_in(i) == this->arg.at(i).size(), "Call input not provided");
  }
  casadi_assert(this->adj1_arg.size() == this->res.size(), "Input 'lam_vdef' not provided");
  for (casadi_int i = 0; i < this->f.n_out(); ++i) {
    casadi_assert(this->f.size_out(i) == this->res.at(i).size(), "Call output not provided");
    casadi_assert(this->adj1_arg.at(i).size() == this->res.at(i).size(),
      "Call adjoint seed not provided");
  }
  // We should make use of the Jacobian blocks here, if available
  if (!this->jac_res.empty())
    casadi_warning("Jacobian blocks currently not reused for gradient calculation");
  // Get/generate the (cached) adjoint function
  // casadi_message("Retrieving the gradient of " + str(this->f));
  this->adj1_f = this->f.reverse(1);
  // casadi_message("Retrieving the gradient of " + str(this->f) + " done");
  // Input expressions for the call to adj1_f
  std::vector<MX> call_in = this->arg;
  call_in.insert(call_in.end(), this->res.begin(), this->res.end());
  call_in.insert(call_in.end(), this->adj1_arg.begin(), this->adj1_arg.end());
  // Create expressions for adjoint sweep and save to struct
  this->adj1_res = this->adj1_f(call_in);
}

void DaeBuilder::CallIO::calc_hess() {
  // Calculate gradient, if needed
  if (this->adj1_f.is_null()) calc_grad();
  // Get/generate the (cached) Hessian function
  // casadi_message("Retrieving the Hessian of " + str(this->f));
  this->H = this->adj1_f.jacobian();
  // casadi_message("Retrieving the Hessian of " + str(this->f) + " done");
  // Input expressions for the call to H
  std::vector<MX> call_in = this->arg;
  call_in.insert(call_in.end(), this->res.begin(), this->res.end());
  call_in.insert(call_in.end(), this->adj1_arg.begin(), this->adj1_arg.end());
  call_in.insert(call_in.end(), this->adj1_res.begin(), this->adj1_res.end());
  // Create expressions for Hessian blocks and save to struct
  this->hess_res = this->H(call_in);
}

const MX& DaeBuilder::CallIO::jac(casadi_int oind, casadi_int iind) const {
  // Flat index
  casadi_int ind = iind + oind * this->arg.size();
  // Return reference
  return this->jac_res.at(ind);
}

const MX& DaeBuilder::CallIO::hess(casadi_int iind1, casadi_int iind2) const {
  // Flat index
  casadi_int ind = iind1 + iind1 * this->adj1_arg.size();
  // Return reference
  return this->hess_res.at(ind);
}

void DaeBuilder::sort_dependent(std::vector<MX>& v, std::vector<MX>& vdef) {
  // Calculate sparsity pattern of dvdef/dv
  Function vfcn("vfcn", {vertcat(v)}, {vertcat(vdef)}, {"v"}, {"vdef"});
  Sparsity Jv = vfcn.jac_sparsity(0, 0);
  // Add diagonal (equation is v-vdef = 0)
  Jv = Jv + Sparsity::diag(Jv.size1());
  // If lower triangular, nothing to do
  if (Jv.is_triu()) return;
  // Perform a Dulmage-Mendelsohn decomposition
  std::vector<casadi_int> rowperm, colperm, rowblock, colblock, coarse_rowblock, coarse_colblock;
  casadi_int nz = Jv.btf(rowperm, colperm, rowblock, colblock, coarse_rowblock, coarse_colblock);
  // Reorder the variables
  std::vector<MX> tmp(v.size());
  for (size_t k = 0; k < v.size(); ++k) tmp[k] = v.at(colperm.at(k));
  std::copy(tmp.begin(), tmp.end(), v.begin());
  // Reorder the equations
  for (size_t k = 0; k < v.size(); ++k) tmp[k] = vdef.at(rowperm.at(k));
  std::copy(tmp.begin(), tmp.end(), vdef.begin());
}

MX Variable::attribute(Attribute att) const {
  switch (att) {
  case MIN:
    return this->min;
  case MAX:
    return this->max;
  case NOMINAL:
    return this->nominal;
  case START:
    return this->start;
  default:
    casadi_error("Cannot process attribute '" + to_string(att) + "'");
    return MX();
  }
}

Function DaeBuilder::attribute_fun(const std::string& fname,
    const std::vector<std::string>& s_in,
    const std::vector<std::string>& s_out) const {
  // Convert inputs to enums
  std::vector<DaeBuilderIn> v_in;
  v_in.reserve(v_in.size());
  for (const std::string& s : s_in) v_in.push_back(to_enum<DaeBuilderIn>(s));
  // Convert outputs to enums
  std::vector<Variable::Attribute> a_out;
  std::vector<DaeBuilderIn> v_out;
  a_out.reserve(s_out.size());
  v_out.reserve(s_out.size());
  for (const std::string& s : s_out) {
    // Locate the underscore divider
    size_t pos = s.find('_');
    casadi_assert(pos < s.size(), "Cannot process \"" + s + "\"");
    // Get attribute
    a_out.push_back(to_enum<Variable::Attribute>(s.substr(0, pos)));
    // Get variable
    v_out.push_back(to_enum<DaeBuilderIn>(s.substr(pos + 1, std::string::npos)));
  }
  // Collect input expressions
  std::vector<MX> f_in;
  f_in.reserve(s_in.size());
  for (DaeBuilderIn v : v_in) f_in.push_back(vertcat(input(v)));
  // Collect output expressions
  std::vector<MX> f_out;
  f_out.reserve(s_out.size());
  for (size_t i = 0; i < s_out.size(); ++i) {
    // Get expressions for which attributes are requested
    std::vector<MX> vars = input(v_out.at(i));
    // Expressions for the attributes
    std::vector<MX> attr;
    attr.reserve(vars.size());
    // Collect attributes
    for (const MX& vi : vars)
      attr.push_back(variable(vi.name()).attribute(a_out.at(i)));
    // Add to output expressions
    f_out.push_back(vertcat(attr));
  }
  // Assemble return function
  return Function(fname, f_in, f_out, s_in, s_out);
}

Function DaeBuilder::dependent_fun(const std::string& fname,
    const std::vector<std::string>& s_in,
    const std::vector<std::string>& s_out) const {
  // Are we calculating d and/or w
  bool calc_d = false, calc_w = false;
  // Convert outputs to enums
  std::vector<DaeBuilderIn> v_out;
  v_out.reserve(v_out.size());
  for (const std::string& s : s_out) {
    DaeBuilderIn e = to_enum<DaeBuilderIn>(s);
    if (e == DAE_BUILDER_D) {
      calc_d = true;
    } else if (e == DAE_BUILDER_W) {
      calc_w = true;
    } else {
      casadi_error("Can only calculate d and/or w");
    }
    v_out.push_back(e);
  }
  // Consistency check
  casadi_assert(calc_d || calc_w, "Nothing to calculate");
  // Convert inputs to enums
  std::vector<DaeBuilderIn> v_in;
  v_in.reserve(v_in.size());
  for (const std::string& s : s_in) {
    DaeBuilderIn e = to_enum<DaeBuilderIn>(s);
    if (calc_d && e == DAE_BUILDER_D) casadi_error("'d' cannot be both input and output");
    if (calc_w && e == DAE_BUILDER_W) casadi_error("'w' cannot be both input and output");
    v_in.push_back(e);
  }
  // Collect input expressions
  std::vector<MX> f_in;
  f_in.reserve(s_in.size());
  for (DaeBuilderIn v : v_in) f_in.push_back(vertcat(input(v)));
  // Collect output expressions
  std::vector<MX> f_out;
  f_out.reserve(s_out.size());
  for (DaeBuilderIn v : v_out) f_out.push_back(vertcat(input(v)));
  // Variables to be substituted
  std::vector<MX> dw, dwdef;
  if (calc_d) {
    dw.insert(dw.end(), this->d.begin(), this->d.end());
    dwdef.insert(dwdef.end(), this->ddef.begin(), this->ddef.end());
  }
  if (calc_w) {
    dw.insert(dw.end(), this->w.begin(), this->w.end());
    dwdef.insert(dwdef.end(), this->wdef.begin(), this->wdef.end());
  }
  // Perform elimination
  substitute_inplace(dw, dwdef, f_out);
  // Assemble return function
  return Function(fname, f_in, f_out, s_in, s_out);
}

void DaeBuilder::prune_d() {
  // If no d, quick return
  if (this->d.empty()) return;
  // Create a dependent function with all inputs except for d itself
  Function dfun = dependent_fun("dfun", {"c"}, {"d"});
  // If no free variables, all good
  if (!dfun.has_free()) return;
  // Print progress
  casadi_message("Eliminating " + str(dfun.get_free()));
  // Variables to be eliminated
  MX elim = vertcat(dfun.free_mx());
  // Create a function for identifying which d cannot be calculated
  dfun = Function("dfun", {vertcat(this->d), elim}, {vertcat(this->ddef)}, {"d", "elim"}, {"ddef"});
  // Seed all elim
  std::vector<bvec_t> elim_sp(dfun.nnz_in("elim"), static_cast<bvec_t>(1));
  // Do not seed d
  std::vector<bvec_t> d_sp(dfun.nnz_in("d"), static_cast<bvec_t>(0));
  // Propagate dependencies to ddef
  std::vector<bvec_t> ddef_sp(dfun.nnz_out("ddef"), static_cast<bvec_t>(0));
  dfun({&d_sp.front(), &elim_sp.front()}, {&ddef_sp.front()});
  // Get vertical offsets in d vector
  std::vector<casadi_int> d_off = offset(this->d);
  // Eliminate dependent variables that depend on any free variable
  std::vector<MX> d_new, ddef_new;
  d_new.reserve(this->d.size());
  ddef_new.reserve(this->ddef.size());
  for (size_t k = 0; k < this->d.size(); ++k) {
    // Does the dependence enter in the definition of the variable?
    bvec_t ddef_sp_any(0);
    for (casadi_int i = d_off.at(k); i < d_off.at(k + 1); ++i) {
      ddef_sp_any = ddef_sp_any | ddef_sp.at(i);
    }
    // If there is a dependence?
    if (ddef_sp_any) {
      // Yes: Eliminate
      casadi_message("Eliminating 'd' that depends on free variables: " + this->d.at(k).name());
    } else {
      // No: Keep
      d_new.push_back(this->d.at(k));
      ddef_new.push_back(this->ddef.at(k));
    }
  }
  // Update d, ddef
  this->d.resize(d_new.size());
  std::copy(d_new.begin(), d_new.end(), this->d.begin());
  this->ddef.resize(ddef_new.size());
  std::copy(ddef_new.begin(), ddef_new.end(), this->ddef.begin());
  // Tail recursive call to handle dependencies of these eliminated dependent parameters
  prune_d();
}

} // namespace casadi

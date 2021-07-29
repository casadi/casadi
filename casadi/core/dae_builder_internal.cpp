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
#include "fmu_function.hpp"

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

DaeBuilderInternal::~DaeBuilderInternal() {
#ifdef WITH_FMU
  if (fmu_) delete fmu_;
#endif // WITH_FMU
}

DaeBuilderInternal::DaeBuilderInternal(const std::string& name, const std::string& path)
    : name_(name), path_(path) {
  clear_cache_ = false;
  number_of_event_indicators_ = 0;
  provides_directional_derivative_ = 0;
  fmu_ = 0;
}

void DaeBuilderInternal::load_fmi_description(const std::string& filename) {
  // Ensure no variables already
  casadi_assert(variables_.empty(), "Instance already has variables");

  // Parse XML file
  XmlFile xml_file("tinyxml");
  XmlNode fmi_desc = xml_file.parse(filename)[0];  // One child; fmiModelDescription

  // Read attributes
  fmi_version_ = fmi_desc.attribute<std::string>("fmiVersion", "");
  model_name_ = fmi_desc.attribute<std::string>("modelName", "");
  guid_ = fmi_desc.attribute<std::string>("guid", "");
  description_ = fmi_desc.attribute<std::string>("description", "");
  author_ = fmi_desc.attribute<std::string>("author", "");
  copyright_ = fmi_desc.attribute<std::string>("copyright", "");
  license_ = fmi_desc.attribute<std::string>("license", "");
  generation_tool_ = fmi_desc.attribute<std::string>("generationTool", "");
  generation_date_and_time_ = fmi_desc.attribute<std::string>("generationDateAndTime", "");
  variable_naming_convention_ = fmi_desc.attribute<std::string>("variableNamingConvention", "");
  number_of_event_indicators_ = fmi_desc.attribute<casadi_int>("numberOfEventIndicators", 0);

  // Process ModelExchange
  if (fmi_desc.has_child("ModelExchange"))
    import_model_exchange(fmi_desc["ModelExchange"]);

  // Process ModelVariables
  casadi_assert(fmi_desc.has_child("ModelVariables"), "Missing 'ModelVariables'");
  import_model_variables(fmi_desc["ModelVariables"]);

  // Process model structure
  if (fmi_desc.has_child("ModelStructure"))
    import_model_structure(fmi_desc["ModelStructure"]);

  // Postprocess / sort variables
  for (size_t k = 0; k < variables_.size(); ++k) {
    // Get reference to the variable
    Variable& v = variables_[k];
    // Sort by types
    if (v.causality == Variable::INDEPENDENT) {
      // Independent (time) variable
      t_.push_back(k);
    } else if (v.causality == Variable::INPUT) {
      u_.push_back(k);
    } else if (v.variability == Variable::CONSTANT) {
      // Named constant
      c_.push_back(k);
      v.beq = v.start;
    } else if (v.variability == Variable::FIXED || v.variability == Variable::TUNABLE) {
      p_.push_back(k);
    } else if (v.variability == Variable::CONTINUOUS) {
      if (v.antiderivative >= 0) {
        // Is the variable needed to calculate other states, algebraic variables?
        if (v.dependency) {
          // Add to list of differential equations
          x_.push_back(k);
          ode_.push_back(v.antiderivative);
        } else {
          // Add to list of quadrature equations
          q_.push_back(k);
          quad_.push_back(v.antiderivative);
        }
      } else if (v.dependency && v.derivative < 0) {
        // Add to list of algebraic equations
        z_.push_back(k);
      }
      // Is it (also) an output variable?
      if (v.causality == Variable::OUTPUT) {
        y_.push_back(k);
        v.beq = v.v;
      }
    } else if (v.dependency) {
      casadi_warning("Cannot sort " + v.name);
    }
  }

  // **** Add binding equations ****
  if (fmi_desc.has_child("equ:BindingEquations")) {
    // Get a reference to the BindingEquations node
    const XmlNode& bindeqs = fmi_desc["equ:BindingEquations"];
    // Loop over binding equations
    for (casadi_int i = 0; i < bindeqs.size(); ++i) {
      // Reference to the binding equation
      const XmlNode& beq_node = bindeqs[i];
      // Get the variable and binding expression
      Variable& var = read_variable(beq_node[0]);
      if (beq_node[1].size() == 1) {
        // Regular expression
        var.beq = read_expr(beq_node[1][0]);
      } else {
        // OpenModelica 1.17 occationally generates integer values without type specifier (bug?)
        casadi_assert(beq_node[1].size() == 0, "Not implemented");
        casadi_int val;
        beq_node[1].get(&val);
        casadi_warning(var.name + " has binding equation without type specifier: " + str(val));
        var.beq = val;
      }
    }
  }

  // **** Add dynamic equations, initial equations ****
  for (bool init_eq : {true, false}) {
    const char* equ = init_eq ? "equ:InitialEquations" : "equ:DynamicEquations";
    if (fmi_desc.has_child(equ)) {
      // Get a reference to the DynamicEquations node
      const XmlNode& dyneqs = fmi_desc[equ];
      // Add equations
      for (casadi_int i = 0; i < dyneqs.size(); ++i) {
        // Get a reference to the variable
        const XmlNode& n = dyneqs[i];
        try {
          // Consistency checks
          casadi_assert_dev(n.name == "equ:Equation");
          casadi_assert_dev(n.size() == 1 && n[0].name == "exp:Sub");
          // Ensure not empty
          if (n[0].size() == 0) {
            casadi_warning(str(equ) + "#" + str(i) + " is empty, ignored.");
            continue;
          }
          // Get the left-hand-sides and right-hand-sides
          const XmlNode& lhs = n[0][0];
          const XmlNode& rhs = n[0][1];
          // Left-hand-side needs to be a variable
          Variable& v = read_variable(lhs);
          // Right-hand-side is the binding equation
          MX beq = read_expr(rhs);
          // Set the equation
          w_.push_back(find(v.name));
          v.beq = beq;
          // Also add to list of initial equations
          if (init_eq) {
            init_lhs_.push_back(v.v);
            init_rhs_.push_back(beq);
          }
        } catch (std::exception& e) {
          uerr() << "Failed to read " << equ << "#" << i << ": " << e.what() << std::endl;
        }
      }
    }
  }
}

void DaeBuilderInternal::load_fmi_functions(const std::string& path) {
  // All expressions
  std::vector<MX> var_xd;
  // All value references
  std::vector<size_t> id_xd, id_yd;
  // Get the IDs and expressions of the initial unknowns
  id_xd.reserve(u_.size() + x_.size());
  var_xd.reserve(u_.size() + x_.size());
  for (auto& v : {u_, x_}) {
    for (size_t k : v) {
      id_xd.push_back(k);
      var_xd.push_back(variable(k).v);
    }
  }
  // Get the IDs of the output variables
  id_yd.reserve(outputs_.size() + derivatives_.size());
  for (auto& v : {outputs_, derivatives_})
    for (size_t k : v) id_yd.push_back(variable(k).value_reference);
  // Create an FMU function
  Dict opts = {
    {"enable_fd", !provides_directional_derivative_},
    {"fd_method", "smoothing"}};
  Function fmu = fmu_fun(name_, {id_xd}, {id_yd}, {"xd"}, {"yd"}, opts);
  add_fun(fmu);
  // Auxiliary variables for xd
  Variable xd(name_ + "_xd");
  xd.v = MX::sym(xd.name, id_xd.size());
  xd.variability = Variable::CONTINUOUS;
  xd.beq = vertcat(var_xd);
  w_.push_back(add_variable(xd.name, xd));
  // Create function call to fmu
  std::vector<MX> fmu_lhs = fmu(std::vector<MX>{xd.v});
  // Auxiliary variables for yd
  Variable yd(name_ + "_yd");
  yd.v = MX::sym(yd.name, id_yd.size());
  yd.variability = Variable::CONTINUOUS;
  yd.beq = fmu_lhs.at(0);
  w_.push_back(add_variable(yd.name, yd));
  // Split up yd.v
  std::vector<MX> yd_split = vertsplit(yd.v);
  // Update binding equations for yd, add to list of dependent variables
  casadi_assert_dev(yd_split.size() == outputs_.size() + derivatives_.size());
  for (size_t k = 0; k < outputs_.size(); ++k) {
    variable(outputs_[k]).beq = yd_split.at(k);
    w_.push_back(outputs_[k]);
  }
  for (size_t k = 0; k < derivatives_.size(); ++k) {
    variable(derivatives_[k]).beq = yd_split.at(outputs_.size() + k);
    w_.push_back(derivatives_[k]);
  }
}

Variable& DaeBuilderInternal::read_variable(const XmlNode& node) {
  // Qualified name
  std::string qn = qualified_name(node);

  // Find and return the variable
  return variable(qn);
}

MX DaeBuilderInternal::read_expr(const XmlNode& node) {
  const std::string& fullname = node.name;
  if (fullname.find("exp:")== std::string::npos) {
    casadi_error("DaeBuilderInternal::read_expr: unknown - expression is supposed to "
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
  } else if (name=="IntegerLiteral" || name=="BooleanLiteral") {
    casadi_int val;
    node.get(&val);
    return val;
  } else if (name=="Instant") {
    double val;
    node.get(&val);
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
    node.get(&val);
    return val;
  } else if (name=="Sin") {
    return sin(read_expr(node[0]));
  } else if (name=="Sqrt") {
    return sqrt(read_expr(node[0]));
  } else if (name=="StringLiteral") {
    casadi_error(node.text);
  } else if (name=="Sub") {
    return read_expr(node[0]) - read_expr(node[1]);
  } else if (name=="Tan") {
    return tan(read_expr(node[0]));
  } else if (name=="Time") {
    return var(t_.at(0));
  } else if (name=="TimedVariable") {
    return read_variable(node[0]).v;
  } else if (name=="FunctionCall") {
    // Get the name of the function
    std::string fname = qualified_name(node["exp:Name"]);
    casadi_warning("Function call to '" + fname + "' incomplete");
    // Collect the arguments
    const XmlNode& args = node["exp:Arguments"];
    std::vector<MX> farg(args.size());
    for (casadi_int i = 0; i < args.size(); ++i) {
      // Lift input arguments
      Variable v("w_" + str(w_.size()));
      v.v = MX::sym(v.name);
      // Add to list of variables
      w_.push_back(add_variable(v.name, v));
      // Set binding expression
      v.beq = read_expr(args[i]);
      // Add to list of function arguments
      farg[i] = v.v;
    }
    // Return argument (scalar for now)
    Variable r("w_" + str(w_.size()));
    r.v = MX::sym(r.name);
    // Add to list of variables
    w_.push_back(add_variable(r.name, r));
    // Return output variable
    return r.v;
  } else if (name=="Array") {
    // Array of arguments
    std::vector<MX> v(node.size());
    for (casadi_int i = 0; i < v.size(); ++i) v[i] = read_expr(node[i]);
    return vertcat(v);
  }

  // throw error if reached this point
  casadi_error("Unknown node: " + name);
}

void DaeBuilderInternal::disp(std::ostream& stream, bool more) const {
  // Assert correctness
  if (more) sanity_check();

  // Print dimensions
  stream << "nx = " << x_.size() << ", "
         << "nz = " << z_.size() << ", "
         << "nq = " << q_.size() << ", "
         << "ny = " << y_.size() << ", "
         << "np = " << p_.size() << ", "
         << "nc = " << c_.size() << ", "
         << "nd = " << d_.size() << ", "
         << "nw = " << w_.size() << ", "
         << "nu = " << u_.size();

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
  if (!t_.empty()) stream << "  t = " << var(t_.at(0)) << std::endl;
  if (!c_.empty()) stream << "  c = " << var(c_) << std::endl;
  if (!p_.empty()) stream << "  p = " << var(p_) << std::endl;
  if (!d_.empty()) stream << "  d = " << var(d_) << std::endl;
  if (!x_.empty()) stream << "  x = " << var(x_) << std::endl;
  if (!z_.empty()) stream << "  z = " << var(z_) << std::endl;
  if (!q_.empty()) stream << "  q = " << var(q_) << std::endl;
  if (!y_.empty()) stream << "  y = " << var(y_) << std::endl;
  if (!w_.empty()) stream << "  w = " << var(w_) << std::endl;
  if (!u_.empty()) stream << "  u = " << var(u_) << std::endl;

  if (!c_.empty()) {
    stream << "Constants" << std::endl;
    for (size_t c : c_) {
      stream << "  " << var(c) << " == " << variable(c).beq << std::endl;
    }
  }

  if (!d_.empty()) {
    stream << "Dependent parameters" << std::endl;
    for (size_t d : d_) {
      stream << "  " << var(d) << " == " << variable(d).beq << std::endl;
    }
  }

  if (!w_.empty()) {
    stream << "Dependent equations" << std::endl;
    for (size_t w : w_) {
      stream << "  " << var(w) << " == " << variable(w).beq << std::endl;
    }
  }

  if (!x_.empty()) {
    stream << "State derivatives" << std::endl;
    for (casadi_int k = 0; k < x_.size(); ++k) {
      stream << "  " << var(x_[k]) << ": " << var(ode_[k]) << std::endl;
    }
  }

  if (!alg_.empty()) {
    stream << "Algebraic equations" << std::endl;
    for (casadi_int k = 0; k < alg_.size(); ++k) {
      stream << "  0 == " << var(alg_[k]) << std::endl;
    }
  }

  if (!q_.empty()) {
    stream << "Quadrature derivatives" << std::endl;
    for (casadi_int k = 0; k < q_.size(); ++k) {
      stream << "  " << var(q_[k]) << ": " << var(quad_[k]) << std::endl;
    }
  }

  if (!init_lhs_.empty()) {
    stream << "Initial equations" << std::endl;
    for (casadi_int k=0; k < init_lhs_.size(); ++k) {
      stream << "  " << str(init_lhs_.at(k)) << " == " << str(init_rhs_.at(k))
        << std::endl;
    }
  }

  if (!when_cond_.empty()) {
    stream << "When statements" << std::endl;
    for (casadi_int k = 0; k < when_cond_.size(); ++k) {
      stream << "  when " << str(when_cond_.at(k)) << " < 0: " << str(when_lhs_.at(k))
        << " := " << str(when_rhs_.at(k)) << std::endl;
    }
  }

  if (!y_.empty()) {
    stream << "Output variables" << std::endl;
    for (size_t y : y_) {
      stream << "  " << var(y) << std::endl;
    }
  }
}

void DaeBuilderInternal::eliminate_quad() {
  // Move all the quadratures to the list of differential states
  x_.insert(x_.end(), q_.begin(), q_.end());
  q_.clear();
}

void DaeBuilderInternal::sort_d() {
  std::vector<MX> d = var(d_), ddef = this->ddef();
  sort_dependent(d, ddef);
  d_.clear();
  for (const MX& e : d) d_.push_back(find(e.name()));
}

void DaeBuilderInternal::sort_w() {
  std::vector<MX> w = var(w_), wdef = this->wdef();
  sort_dependent(w, wdef);
  w_.clear();
  for (const MX& e : w) w_.push_back(find(e.name()));
}

void DaeBuilderInternal::sort_z(const std::vector<std::string>& z_order) {
  // Make sure lengths agree
  casadi_assert(z_order.size() == z_.size(), "Dimension mismatch");
  // Mark existing components in z
  std::vector<bool> old_z(variables_.size(), false);
  for (size_t i : z_) old_z.at(i) = true;
  // New vector of z
  std::vector<size_t> new_z;
  new_z.reserve(z_order.size());
  for (const std::string& s : z_order) {
    size_t i = find(s);
    casadi_assert(old_z.at(i), "Variable \"" + s + "\" is not an algebraic variable.");
    new_z.push_back(i);
  }
  // Success: Update z
  std::copy(new_z.begin(), new_z.end(), z_.begin());
}

void DaeBuilderInternal::clear_in(const std::string& v) {
  switch (to_enum<DaeBuilderInternalIn>(v)) {
  case DAE_BUILDER_T: return t_.clear();
  case DAE_BUILDER_P: return p_.clear();
  case DAE_BUILDER_U: return u_.clear();
  case DAE_BUILDER_X: return x_.clear();
  case DAE_BUILDER_Z: return z_.clear();
  case DAE_BUILDER_Q: return q_.clear();
  case DAE_BUILDER_C: return c_.clear();
  case DAE_BUILDER_D: return d_.clear();
  case DAE_BUILDER_W: return w_.clear();
  case DAE_BUILDER_Y: return y_.clear();
  default: break;
  }
  casadi_error("Cannot clear input: " + v);
}

void DaeBuilderInternal::clear_out(const std::string& v) {
  switch (to_enum<DaeBuilderInternalOut>(v)) {
  case DAE_BUILDER_ODE: return ode_.clear();
  case DAE_BUILDER_ALG: return alg_.clear();
  case DAE_BUILDER_QUAD: return quad_.clear();
  default: break;
  }
  casadi_error("Cannot clear output: " + v);
}

void DaeBuilderInternal::prune(bool prune_p, bool prune_u) {
  // Function inputs and outputs
  std::vector<MX> f_in, f_out, v;
  std::vector<std::string> f_in_name, f_out_name;
  // Collect all DAE input variables with at least one entry, skip u
  for (casadi_int i = 0; i != DAE_BUILDER_NUM_IN; ++i) {
    if (prune_p && i == DAE_BUILDER_P) continue;
    if (prune_u && i == DAE_BUILDER_U) continue;
    v = input(static_cast<DaeBuilderInternalIn>(i));
    if (!v.empty()) {
      f_in.push_back(vertcat(v));
      f_in_name.push_back(to_string(static_cast<DaeBuilderInternalIn>(i)));
    }
  }
  // Collect all DAE output variables with at least one entry
  for (casadi_int i = 0; i != DAE_BUILDER_NUM_OUT; ++i) {
    v = output(static_cast<DaeBuilderInternalOut>(i));
    if (!v.empty()) {
      f_out.push_back(vertcat(v));
      f_out_name.push_back(to_string(static_cast<DaeBuilderInternalOut>(i)));
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
    for (size_t i = 0; i < p_.size(); ++i) {
      if (!free_variables.at(p_.at(i))) p_.at(np++) = p_.at(i);
    }
    p_.resize(np);
  }
  // Prune u
  if (prune_u) {
    size_t nu = 0;
    for (size_t i = 0; i < u_.size(); ++i) {
      if (!free_variables.at(u_.at(i))) u_.at(nu++) = u_.at(i);
    }
    u_.resize(nu);
  }
}

void DaeBuilderInternal::tear() {
  // Prefix
  const std::string res_prefix = "res__";
  // Get residual variables, iteration variables
  std::vector<std::string> res, iv, iv_on_hold;
  tearing_variables(&res, &iv, &iv_on_hold);
  // All iteration variables
  std::set<std::string> iv_set;
  for (auto& e : iv) iv_set.insert(e);
  for (auto& e : iv_on_hold) iv_set.insert(e);
  // Remove any (held or not held) iteration variables, equations from z and alg
  size_t sz = 0;
  casadi_assert(z_.size() == alg_.size(), "z and alg have different lengths");
  for (size_t k = 0; k < z_.size(); ++k) {
    if (!iv_set.count(variable(z_[k]).name)) {
      // Non-iteration variable: Keep
      z_.at(k) = z_.at(sz);
      alg_.at(k) = alg_.at(sz);
      sz++;
    }
  }
  z_.resize(sz);
  alg_.resize(sz);
  // Remove any (held or not held) iteration variables, equations from u
  sz = 0;
  for (size_t k = 0; k < u_.size(); ++k) {
    if (!iv_set.count(variable(u_[k]).name)) {
      // Non-iteration variable: Keep
      u_.at(k) = u_.at(sz++);
    }
  }
  u_.resize(sz);
  // Add algebraic variables
  for (auto& e : iv) z_.push_back(find(e));
  // Add residual variables
  for (auto& e : res) alg_.push_back(find(e));
  // Add output variables
  for (auto& e : iv_on_hold) u_.push_back(find(e));
}

void DaeBuilderInternal::tearing_variables(std::vector<std::string>* res,
    std::vector<std::string>* iv, std::vector<std::string>* iv_on_hold) const {
  // Clear output
  if (res) res->clear();
  if (iv) iv->clear();
  if (iv_on_hold) iv_on_hold->clear();
  // Prefix
  const std::string res_prefix = "res__";
  // Collect hold indices
  std::vector<MX> r_hold, iv_hold;
  // Any hold variable?
  bool any_hold = false;
  // Collect residual variables, iteration variables, expression for hold indices, if any
  for (const Variable& v : variables_) {
    // Residual variables are specified with a "res__" prefix
    if (v.name.rfind(res_prefix, 0) == 0) {
      // Process iteration variable name, names of hold markers
      std::string iv_name, res_hold_name, iv_hold_name;
      // Iteration variable, hold markers are contained in the remainder of the name
      try {
        size_t pos = res_prefix.size();
        // Find the next "__", if any
        size_t end = v.name.find("__", pos);
        if (end == std::string::npos) end = v.name.size();
        // Look up iteration variable
        iv_name = v.name.substr(pos, end - pos);
        // Ensure that the variable exists
        casadi_assert(has_variable(iv_name), "No such variable: " + iv_name);
        // Get hold indices, if any
        if (end != v.name.size()) {
          // Find next "__", read hold index for residual variable
          pos = end + 2;
          end = v.name.find("__", pos);
          if (end == std::string::npos) end = v.name.size();
          res_hold_name = v.name.substr(pos, end - pos);
          // Ensure that the variable exists
          casadi_assert(has_variable(res_hold_name), "No such variable: " + res_hold_name);
          // The remainder of the name contains iv_hold_name
          if (end != v.name.size()) {
            iv_hold_name = v.name.substr(end + 2);
            casadi_assert(has_variable(iv_hold_name), "No such variable: " + iv_hold_name);
          }
        }
      } catch (std::exception& e) {
        // Generate warning
        casadi_warning("Cannot process residual variable: " + v.name + ":" + std::string(e.what()));
        continue;
      }
      // Add residual variable, corresponding hold variable
      if (res_hold_name.empty()) {
        r_hold.push_back(false);
      } else {
        any_hold = true;
        r_hold.push_back(variable(res_hold_name).v);
        casadi_assert(r_hold.back().is_scalar(), "Non-scalar hold variable for " + res_hold_name);
      }
      if (res) res->push_back(v.name);
      // Add iteration variable, corresponding hold variable
      if (iv_hold_name.empty()) {
        iv_hold.push_back(false);
      } else {
        any_hold = true;
        iv_hold.push_back(variable(iv_hold_name).v);
        casadi_assert(iv_hold.back().is_scalar(), "Non-scalar hold variable for " + iv_hold_name);
      }
      if (iv) iv->push_back(iv_name);
    }
  }
  // Evaluate hold variables, if needed
  if (any_hold) {
    try {
      // Get start attributes for p
      Function startfun_p = attribute_fun("startfun_p", {}, {"start_p"});
      if (startfun_p.has_free()) {
        casadi_error("startfun has free variables: " + str(startfun_p.get_free()));
      }
      DM p0 = startfun_p(std::vector<DM>{}).at(0);
      // Create function to evaluate the hold attributes
      Function holdfun("holdfun", {vertcat(var(p_))},
        {vertcat(r_hold), vertcat(iv_hold)}, {"p"}, {"r_hold", "iv_hold"});
      if (holdfun.has_free()) {
        casadi_error("holdfun has free variables: " + str(holdfun.get_free()));
      }
      // Evaluate holdfun to get hold attributes
      std::vector<DM> hold0 = holdfun(std::vector<DM>{p0});
      std::vector<double> r_hold0 = hold0.at(0).nonzeros();
      std::vector<double> iv_hold0 = hold0.at(1).nonzeros();
      casadi_assert_dev(r_hold0.size() == res->size());
      casadi_assert_dev(iv_hold0.size() == iv->size());
      // Remove hold variables from residual variables
      size_t sz = 0;
      if (res) {
        for (size_t k = 0; k < res->size(); ++k) {
          if (!static_cast<bool>(r_hold0.at(k))) {
            res->at(sz++) = res->at(k);
          }
        }
        res->resize(sz);
      }
      // Remove hold variables from iteration variables
      sz = 0;
      for (size_t k = 0; k < iv->size(); ++k) {
        if (!static_cast<bool>(iv_hold0.at(k))) {
          if (iv) iv->at(sz++) = iv->at(k);
        } else {
          if (iv_on_hold) iv_on_hold->push_back(iv->at(k));
        }
      }
      if (iv) iv->resize(sz);
    } catch (std::exception& e) {
      // Warning instead of error
      casadi_warning("Failed to evaluate hold variables: " + std::string(e.what()));
    }
  }
}

bool DaeBuilderInternal::has_variable(const std::string& name) const {
  return varind_.find(name) != varind_.end();
}

size_t DaeBuilderInternal::add_variable(const std::string& name, const Variable& var) {
  // Try to find the component
  casadi_assert(!has_variable(name), "Variable \"" + name + "\" has already been added.");
  // Index of the variable
  size_t ind = variables_.size();
  // Add to the map of all variables
  varind_[name] = ind;
  variables_.push_back(var);
  // Clear cache
  clear_cache_ = true;
  // Return index of new variable
  return ind;
}

void DaeBuilderInternal::sanity_check() const {
  // Time
  if (!t_.empty()) {
    casadi_assert(t_.size() == 1, "At most one time variable allowed");
    casadi_assert(var(t_[0]).is_scalar(), "Non-scalar time t");
  }

  // Differential states
  casadi_assert(x_.size() == ode_.size(), "x and ode have different lengths");
  for (casadi_int i = 0; i < x_.size(); ++i) {
    casadi_assert(var(x_[i]).size() == var(ode_[i]).size(), "ode has wrong dimensions");
  }

  // Algebraic variables/equations
  casadi_assert(z_.size() == alg_.size(), "z and alg have different lengths");
  for (casadi_int i = 0; i < z_.size(); ++i) {
    casadi_assert(var(z_[i]).size() == var(alg_[i]).size(), "alg has wrong dimensions");
  }

  // Quadrature states/equations
  casadi_assert(q_.size() == quad_.size(), "q and quad have different lengths");
  for (casadi_int i = 0; i < q_.size(); ++i) {
    casadi_assert(var(q_[i]).size() == var(quad_[i]).size(), "quad has wrong dimensions");
  }

  // Initial equations
  casadi_assert(init_lhs_.size() == init_rhs_.size(),
    "init_lhs and init_rhs have different lengths");

  // When statements
  casadi_assert(when_cond_.size() == when_lhs_.size() && when_lhs_.size() == when_rhs_.size(),
    "when_cond, when_lhs and when_rhs must all have the the same length");
}

std::string DaeBuilderInternal::qualified_name(const XmlNode& nn) {
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
      nn[i]["exp:ArraySubscripts"]["exp:IndexExpression"]["exp:IntegerLiteral"].get(&ind);
      qn << "[" << ind << "]";
    }
  }

  // Return the name
  return qn.str();
}

const MX& DaeBuilderInternal::var(const std::string& name) const {
  return variable(name).v;
}

MX DaeBuilderInternal::der(const std::string& name) const {
  return variables_.at(variable(name).derivative).v;
}

MX DaeBuilderInternal::der(const MX& var) const {
  casadi_assert_dev(var.is_column() && var.is_symbolic());
  return der(var.name());
}

void DaeBuilderInternal::eliminate_w() {
  // Quick return if no w
  if (w_.empty()) return;
  // Clear cache after this
  clear_cache_ = true;
  // Ensure variables are sorted
  sort_w();
  // Expressions where the variables are also being used
  std::vector<MX> ex;
  for (const Variable& v : variables_) {
    if (!v.min.is_constant()) ex.push_back(v.min);
    if (!v.max.is_constant()) ex.push_back(v.max);
    if (!v.nominal.is_constant()) ex.push_back(v.nominal);
    if (!v.start.is_constant()) ex.push_back(v.start);
    if (!v.beq.is_constant()) ex.push_back(v.beq);
  }
  // Perform elimination
  std::vector<MX> w = var(w_);
  std::vector<MX> wdef = this->wdef();
  substitute_inplace(w, wdef, ex);
  // Clear list of dependent variables
  w_.clear();
  // Get variable attributes
  auto it = ex.begin();
  for (Variable& v : variables_) {
    if (!v.min.is_constant()) v.min = *it++;
    if (!v.max.is_constant()) v.max = *it++;
    if (!v.nominal.is_constant()) v.nominal = *it++;
    if (!v.start.is_constant()) v.start = *it++;
    if (!v.beq.is_constant()) v.beq = *it++;
  }
  // Consistency check
  casadi_assert_dev(it == ex.end());
}

void DaeBuilderInternal::lift(bool lift_shared, bool lift_calls) {
  // Not tested if w is non-empty before
  if (!w_.empty()) casadi_warning("'w' already has entries");
  // Expressions where the variables are also being used
  std::vector<MX> ex;
  for (auto& vv : {ode_, alg_, quad_, y_})
    for (size_t v : vv) ex.push_back(variable(v).beq);
  // Lift expressions
  std::vector<MX> new_w, new_wdef;
  Dict opts{{"lift_shared", lift_shared}, {"lift_calls", lift_calls},
    {"prefix", "w_"}, {"suffix", ""}, {"offset", static_cast<casadi_int>(w_.size())}};
  extract(ex, new_w, new_wdef, opts);
  // Register as dependent variables
  for (size_t i = 0; i < new_w.size(); ++i) {
    Variable v(new_w.at(i).name());
    v.v = new_w.at(i);
    v.beq = new_wdef.at(i);
    w_.push_back(add_variable(v.name, v));
  }
  // Get expressions
  auto it = ex.begin();
  for (auto& vv : {ode_, alg_, quad_, y_})
    for (size_t v : vv) variable(v).beq = *it++;
  // Consistency check
  casadi_assert_dev(it == ex.end());
}

std::string to_string(DaeBuilderInternal::DaeBuilderInternalIn v) {
  switch (v) {
  case DaeBuilderInternal::DAE_BUILDER_T: return "t";
  case DaeBuilderInternal::DAE_BUILDER_P: return "p";
  case DaeBuilderInternal::DAE_BUILDER_U: return "u";
  case DaeBuilderInternal::DAE_BUILDER_X: return "x";
  case DaeBuilderInternal::DAE_BUILDER_Z: return "z";
  case DaeBuilderInternal::DAE_BUILDER_Q: return "q";
  case DaeBuilderInternal::DAE_BUILDER_C: return "c";
  case DaeBuilderInternal::DAE_BUILDER_D: return "d";
  case DaeBuilderInternal::DAE_BUILDER_W: return "w";
  case DaeBuilderInternal::DAE_BUILDER_Y: return "y";
  default: break;
  }
  return "";
}

std::string to_string(DaeBuilderInternal::DaeBuilderInternalOut v) {
  switch (v) {
  case DaeBuilderInternal::DAE_BUILDER_ODE: return "ode";
  case DaeBuilderInternal::DAE_BUILDER_ALG: return "alg";
  case DaeBuilderInternal::DAE_BUILDER_QUAD: return "quad";
  case DaeBuilderInternal::DAE_BUILDER_DDEF: return "ddef";
  case DaeBuilderInternal::DAE_BUILDER_WDEF: return "wdef";
  case DaeBuilderInternal::DAE_BUILDER_YDEF: return "ydef";
  default: break;
  }
  return "";
}

std::vector<MX> DaeBuilderInternal::input(DaeBuilderInternalIn ind) const {
  switch (ind) {
  case DAE_BUILDER_T: return var(t_);
  case DAE_BUILDER_C: return var(c_);
  case DAE_BUILDER_P: return var(p_);
  case DAE_BUILDER_D: return var(d_);
  case DAE_BUILDER_W: return var(w_);
  case DAE_BUILDER_U: return var(u_);
  case DAE_BUILDER_X: return var(x_);
  case DAE_BUILDER_Z: return var(z_);
  case DAE_BUILDER_Q: return var(q_);
  case DAE_BUILDER_Y: return var(y_);
  default: return std::vector<MX>{};
  }
}

std::vector<MX> DaeBuilderInternal::input(const std::vector<DaeBuilderInternalIn>& ind) const {
  std::vector<MX> ret(ind.size());
  for (casadi_int i=0; i<ind.size(); ++i) {
    ret[i] = vertcat(input(ind[i]));
  }
  return ret;
}

std::vector<MX> DaeBuilderInternal::output(DaeBuilderInternalOut ind) const {
  switch (ind) {
  case DAE_BUILDER_ODE: return ode();
  case DAE_BUILDER_ALG: return alg();
  case DAE_BUILDER_QUAD: return quad();
  case DAE_BUILDER_DDEF: return ddef();
  case DAE_BUILDER_WDEF: return wdef();
  case DAE_BUILDER_YDEF: return ydef();
  default: return std::vector<MX>();
  }
}

std::vector<MX> DaeBuilderInternal::output(const std::vector<DaeBuilderInternalOut>& ind) const {
  std::vector<MX> ret(ind.size());
  for (casadi_int i=0; i<ind.size(); ++i) {
    ret[i] = vertcat(output(ind[i]));
  }
  return ret;
}

void DaeBuilderInternal::add_lc(const std::string& name, const std::vector<std::string>& f_out) {
  // Make sure object valid
  sanity_check();

  // Make sure name is valid
  casadi_assert(!name.empty(), "DaeBuilderInternal::add_lc: \"name\" is empty");
  for (std::string::const_iterator i=name.begin(); i!=name.end(); ++i) {
    casadi_assert(isalnum(*i),
                          "DaeBuilderInternal::add_lc: \"name\" must be alphanumeric");
  }

  // Consistency checks
  casadi_assert(!f_out.empty(), "DaeBuilderInternal::add_lc: Linear combination is empty");
  std::vector<bool> in_use(DAE_BUILDER_NUM_OUT, false);
  for (casadi_int i=0; i<f_out.size(); ++i) {
    DaeBuilderInternalOut oind = to_enum<DaeBuilderInternalOut>(f_out[i]);
    casadi_assert(!in_use[oind], "DaeBuilderInternal::add_lc: Duplicate expression " + f_out[i]);
    in_use[oind] = true;
  }

  std::vector<std::string>& ret1 = lc_[name];
  if (!ret1.empty()) casadi_warning("DaeBuilderInternal::add_lc: Overwriting " << name);
  ret1 = f_out;
}

Function DaeBuilderInternal::create(const std::string& fname,
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
      for (std::string& s : *s_io) std::replace(s.begin(), s.end(), '_', ':');
    }
    // Recursive call
    return create(fname, s_in_mod, s_out_mod, sx, lifted_calls);
  }
  // Check if dependent variables are given and needed
  bool elim_w = false;
  if (!w_.empty()) {
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
    for (const MX& vdef_comp : wdef()) {
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
  std::vector<casadi_int> h_offsets = offset(var(w_));
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
  for (size_t i = 0; i < w_.size(); ++i) {
    v_map[var(w_.at(i)).get()] = i;
  }
  // Definitions of w
  std::vector<MX> wdef = this->wdef();
  // Collect all the call nodes
  std::map<MXNode*, CallIO> call_nodes;
  for (size_t vdefind = 0; vdefind < wdef.size(); ++vdefind) {
    // Current element handled
    const MX& vdefref = wdef.at(vdefind);
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

MX DaeBuilderInternal::jac_vdef_v_from_calls(std::map<MXNode*, CallIO>& call_nodes,
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
  // Definitions of w
  std::vector<MX> wdef = this->wdef();
  // Collect all Jacobian blocks
  for (size_t vdefind = 0; vdefind < wdef.size(); ++vdefind) {
    // Current element handled
    const MX& vdefref = wdef.at(vdefind);
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

MX DaeBuilderInternal::hess_v_v_from_calls(std::map<MXNode*, CallIO>& call_nodes,
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
  for (size_t vind1 = 0; vind1 < w_.size(); ++vind1) {
    // Current element handled
    const MX& vref = var(w_.at(vind1));
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

void DaeBuilderInternal::clear_cache() const {
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

const Function& DaeBuilderInternal::oracle(bool sx, bool elim_w, bool lifted_calls) const {
  // Clear cache now, if necessary
  if (clear_cache_) clear_cache();
  // Create an MX oracle, if needed
  if (oracle_[false][elim_w][lifted_calls].is_null()) {
    // Oracle function inputs and outputs
    std::vector<MX> f_in, f_out, v;
    std::vector<std::string> f_in_name, f_out_name;
    // Index for wdef
    casadi_int wdef_ind = -1;
    // Options consistency check
    casadi_assert(!(elim_w && lifted_calls), "Incompatible options");
    // Do we need to substitute out v
    bool subst_v = false;
    // Collect all DAE input variables with at least one entry
    for (casadi_int i = 0; i != DAE_BUILDER_NUM_IN; ++i) {
      if (i == DAE_BUILDER_Y) continue;  // fixme
      v = input(static_cast<DaeBuilderInternalIn>(i));
      if (!v.empty()) {
        if (elim_w && i == DAE_BUILDER_W) {
          subst_v = true;
        } else {
          f_in.push_back(vertcat(v));
          f_in_name.push_back(to_string(static_cast<DaeBuilderInternalIn>(i)));
        }
      }
    }
    // Collect all DAE output variables with at least one entry
    for (casadi_int i = 0; i != DAE_BUILDER_NUM_OUT; ++i) {
      v = output(static_cast<DaeBuilderInternalOut>(i));
      if (!v.empty()) {
        if (i == DAE_BUILDER_WDEF) wdef_ind = f_out.size();
        f_out.push_back(vertcat(v));
        f_out_name.push_back(to_string(static_cast<DaeBuilderInternalOut>(i)));
      }
    }
    // Eliminate v from inputs
    if (subst_v) {
      // Dependent variable definitions
      std::vector<MX> wdef = this->wdef();
      // Perform in-place substitution
      substitute_inplace(var(w_), wdef, f_out, false);
    } else if (lifted_calls && wdef_ind >= 0) {
      // Dependent variable definitions
      std::vector<MX> wdef = this->wdef();
      // Remove references to call nodes
      for (MX& wdefref : wdef) {
        if (wdefref.is_output()) wdefref = MX::zeros(wdefref.sparsity());
      }
      // Save to oracle outputs
      f_out.at(wdef_ind) = vertcat(wdef);
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

void DaeBuilderInternal::CallIO::calc_jac() {
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

void DaeBuilderInternal::CallIO::calc_grad() {
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

void DaeBuilderInternal::CallIO::calc_hess() {
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

const MX& DaeBuilderInternal::CallIO::jac(casadi_int oind, casadi_int iind) const {
  // Flat index
  casadi_int ind = iind + oind * this->arg.size();
  // Return reference
  return this->jac_res.at(ind);
}

const MX& DaeBuilderInternal::CallIO::hess(casadi_int iind1, casadi_int iind2) const {
  // Flat index
  casadi_int ind = iind1 + iind1 * this->adj1_arg.size();
  // Return reference
  return this->hess_res.at(ind);
}

void DaeBuilderInternal::sort_dependent(std::vector<MX>& v, std::vector<MX>& vdef) {
  // Form function to evaluate dependent variables
  Function vfcn("vfcn", {vertcat(v)}, {vertcat(vdef)}, {"v"}, {"vdef"});
  // Is any variable vector-valued?
  bool any_vector_valued = false;
  for (const MX& v_i : v) {
    casadi_assert(!v_i.is_empty(), "Cannot have zero-dimension dependent variables");
    if (!v_i.is_scalar()) {
      any_vector_valued = true;
      break;
    }
  }
  // If vector-valued variables exists, collapse them
  if (any_vector_valued) {
    // New v corresponding to one scalar input per v argument
    std::vector<MX> vfcn_in(v), vfcn_arg(v);
    for (size_t i = 0; i < v.size(); ++i) {
      if (!v.at(i).is_scalar()) {
        vfcn_in.at(i) = MX::sym(v.at(i).name());
        vfcn_arg.at(i) = repmat(vfcn_in.at(i), v.at(i).size1());
      }
    }
    // Wrap vfcn
    std::vector<MX> vfcn_out = vfcn(vertcat(vfcn_arg));
    vfcn_out = vertsplit(vfcn_out.at(0), offset(v));
    // Collapse vector-valued outputs
    for (size_t i = 0; i < v.size(); ++i) {
      if (!v.at(i).is_scalar()) {
        vfcn_out.at(i) = dot(vfcn_out.at(i), vfcn_out.at(i));
      }
    }
    // Recreate vfcn with smaller dimensions
    vfcn = Function(vfcn.name(), {vertcat(vfcn_in)}, {vertcat(vfcn_out)},
      vfcn.name_in(), vfcn.name_out());
  }
  // Calculate sparsity pattern of dvdef/dv
  Sparsity Jv = vfcn.jac_sparsity(0, 0);
  // Add diagonal (equation is v-vdef = 0)
  Jv = Jv + Sparsity::diag(Jv.size1());
  // If lower triangular, nothing to do
  if (Jv.is_triu()) return;
  // Perform a Dulmage-Mendelsohn decomposition
  std::vector<casadi_int> rowperm, colperm, rowblock, colblock, coarse_rowblock, coarse_colblock;
  (void)Jv.btf(rowperm, colperm, rowblock, colblock, coarse_rowblock, coarse_colblock);
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

Function DaeBuilderInternal::attribute_fun(const std::string& fname,
    const std::vector<std::string>& s_in,
    const std::vector<std::string>& s_out) const {
  // Convert inputs to enums
  std::vector<DaeBuilderInternalIn> v_in;
  v_in.reserve(v_in.size());
  for (const std::string& s : s_in) v_in.push_back(to_enum<DaeBuilderInternalIn>(s));
  // Convert outputs to enums
  std::vector<Variable::Attribute> a_out;
  std::vector<DaeBuilderInternalIn> v_out;
  a_out.reserve(s_out.size());
  v_out.reserve(s_out.size());
  for (const std::string& s : s_out) {
    // Locate the underscore divider
    size_t pos = s.find('_');
    casadi_assert(pos < s.size(), "Cannot process \"" + s + "\"");
    // Get attribute
    a_out.push_back(to_enum<Variable::Attribute>(s.substr(0, pos)));
    // Get variable
    v_out.push_back(to_enum<DaeBuilderInternalIn>(s.substr(pos + 1, std::string::npos)));
  }
  // Collect input expressions
  std::vector<MX> f_in;
  f_in.reserve(s_in.size());
  for (DaeBuilderInternalIn v : v_in) f_in.push_back(vertcat(input(v)));
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
    // Handle dependency on w
    if (depends_on(f_out.back(), vertcat(var(w_)))) {
      // Make copies of w and wdef and sort
      std::vector<MX> w_sorted = var(w_), wdef = this->wdef();
      sort_dependent(w_sorted, wdef);
      // Eliminate dependency on w
      substitute_inplace(w_sorted, wdef, attr, false);
      // Update f_out
      f_out.back() = vertcat(attr);
    }
    // Handle interdependencies
    if (depends_on(f_out.back(), vertcat(vars))) {
      // Make copies of vars and attr and sort
      std::vector<MX> vars_sorted = vars, attr_sorted = attr;
      sort_dependent(vars_sorted, attr_sorted);
      // Eliminate interdependencies
      substitute_inplace(vars_sorted, attr_sorted, attr, false);
      // Update f_out
      f_out.back() = vertcat(attr);
    }
  }
  // Assemble return function
  return Function(fname, f_in, f_out, s_in, s_out);
}

Function DaeBuilderInternal::dependent_fun(const std::string& fname,
    const std::vector<std::string>& s_in,
    const std::vector<std::string>& s_out) const {
  // Are we calculating d and/or w
  bool calc_d = false, calc_w = false;
  // Convert outputs to enums
  std::vector<DaeBuilderInternalIn> v_out;
  v_out.reserve(v_out.size());
  for (const std::string& s : s_out) {
    DaeBuilderInternalIn e = to_enum<DaeBuilderInternalIn>(s);
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
  std::vector<DaeBuilderInternalIn> v_in;
  v_in.reserve(v_in.size());
  for (const std::string& s : s_in) {
    DaeBuilderInternalIn e = to_enum<DaeBuilderInternalIn>(s);
    if (calc_d && e == DAE_BUILDER_D) casadi_error("'d' cannot be both input and output");
    if (calc_w && e == DAE_BUILDER_W) casadi_error("'w' cannot be both input and output");
    v_in.push_back(e);
  }
  // Collect input expressions
  std::vector<MX> f_in;
  f_in.reserve(s_in.size());
  for (DaeBuilderInternalIn v : v_in) f_in.push_back(vertcat(input(v)));
  // Collect output expressions
  std::vector<MX> f_out;
  f_out.reserve(s_out.size());
  for (DaeBuilderInternalIn v : v_out) f_out.push_back(vertcat(input(v)));
  // Variables to be substituted
  std::vector<MX> dw, dwdef;
  if (calc_d) {
    std::vector<MX> d = var(d_);
    dw.insert(dw.end(), d.begin(), d.end());
    std::vector<MX> ddef = this->ddef();
    dwdef.insert(dwdef.end(), ddef.begin(), ddef.end());
  }
  if (calc_w) {
    std::vector<MX> w = var(w_);
    dw.insert(dw.end(), w.begin(), w.end());
    std::vector<MX> wdef = this->wdef();
    dwdef.insert(dwdef.end(), wdef.begin(), wdef.end());
  }
  // Perform elimination
  substitute_inplace(dw, dwdef, f_out);
  // Assemble return function
  return Function(fname, f_in, f_out, s_in, s_out);
}

Function DaeBuilderInternal::fmu_fun(const std::string& name,
    const std::vector<std::vector<size_t>>& id_in,
    const std::vector<std::vector<size_t>>& id_out,
    const std::vector<std::string>& name_in,
    const std::vector<std::string>& name_out,
    const Dict& opts) const {
#ifdef WITH_FMU
  return Function::create(new FmuFunction(name, shared_from_this<DaeBuilder>(),
    id_in, id_out, name_in, name_out), opts);
#else  // WITH_FMU
  casadi_error("FMU support not enabled. Recompile CasADi with 'WITH_FMU=ON'");
  return Function();
#endif  // WITH_FMU
}

Function DaeBuilderInternal::gather_eq() const {
  // Output expressions
  std::vector<MX> f_out;
  // Names of outputs
  std::vector<std::string> f_out_name;
  // Get all expressions
  for (casadi_int i = 0; i != DAE_BUILDER_NUM_OUT; ++i) {
    std::vector<MX> v = output(static_cast<DaeBuilderInternalOut>(i));
    if (!v.empty()) {
      f_out.push_back(vertcat(v));
      f_out_name.push_back(to_string(static_cast<DaeBuilderInternalOut>(i)));
    }
  }
  // Construct function
  return Function("all_eq", {}, f_out, {}, f_out_name);
}

std::vector<MX> DaeBuilderInternal::cdef() const {
  std::vector<MX> ret;
  ret.reserve(c_.size());
  for (size_t c : c_) ret.push_back(variable(c).beq);
  return ret;
}

std::vector<MX> DaeBuilderInternal::ddef() const {
  std::vector<MX> ret;
  ret.reserve(d_.size());
  for (size_t d : d_) ret.push_back(variable(d).beq);
  return ret;
}

std::vector<MX> DaeBuilderInternal::wdef() const {
  std::vector<MX> ret;
  ret.reserve(w_.size());
  for (size_t w : w_) ret.push_back(variable(w).beq);
  return ret;
}

std::vector<MX> DaeBuilderInternal::ydef() const {
  std::vector<MX> ret;
  ret.reserve(y_.size());
  for (size_t v : y_) ret.push_back(variable(v).beq);
  return ret;
}

std::vector<MX> DaeBuilderInternal::ode() const {
  std::vector<MX> ret;
  ret.reserve(ode_.size());
  for (size_t v : ode_) ret.push_back(variable(v).beq);
  return ret;
}

std::vector<MX> DaeBuilderInternal::alg() const {
  std::vector<MX> ret;
  ret.reserve(alg_.size());
  for (size_t v : alg_) ret.push_back(variable(v).beq);
  return ret;
}

std::vector<MX> DaeBuilderInternal::quad() const {
  std::vector<MX> ret;
  ret.reserve(quad_.size());
  for (size_t v : quad_) ret.push_back(variable(v).beq);
  return ret;
}

MX DaeBuilderInternal::add_t(const std::string& name) {
  casadi_assert(t_.empty(), "'t' already defined");
  Variable v(name);
  v.v = MX::sym(name);
  v.causality = Variable::INDEPENDENT;
  t_.push_back(add_variable(name, v));
  return v.v;
}

MX DaeBuilderInternal::add_p(const std::string& name, casadi_int n) {
  Variable v(name);
  v.v = MX::sym(name, n);
  v.variability = Variable::FIXED;
  v.causality = Variable::INPUT;
  p_.push_back(add_variable(name, v));
  return v.v;
}

MX DaeBuilderInternal::add_u(const std::string& name, casadi_int n) {
  Variable v(name);
  v.v = MX::sym(name, n);
  v.variability = Variable::CONTINUOUS;
  v.causality = Variable::INPUT;
  u_.push_back(add_variable(name, v));
  return v.v;
}

MX DaeBuilderInternal::add_x(const std::string& name, casadi_int n) {
  Variable v(name);
  v.v = MX::sym(name, n);
  v.variability = Variable::CONTINUOUS;
  v.causality = Variable::LOCAL;
  x_.push_back(add_variable(name, v));
  return v.v;
}

MX DaeBuilderInternal::add_z(const std::string& name, casadi_int n) {
  Variable v(name);
  v.v = MX::sym(name, n);
  v.variability = Variable::CONTINUOUS;
  v.causality = Variable::LOCAL;
  z_.push_back(add_variable(name, v));
  return v.v;
}

MX DaeBuilderInternal::add_q(const std::string& name, casadi_int n) {
  Variable v(name);
  v.v = MX::sym(name, n);
  v.variability = Variable::CONTINUOUS;
  v.causality = Variable::LOCAL;
  q_.push_back(add_variable(name, v));
  return v.v;
}

MX DaeBuilderInternal::add_c(const std::string& name, const MX& new_cdef) {
  Variable v(name);
  v.v = MX::sym(name);
  v.variability = Variable::CONSTANT;
  v.beq = new_cdef;
  c_.push_back(add_variable(name, v));
  return v.v;
}

MX DaeBuilderInternal::add_d(const std::string& name, const MX& new_ddef) {
  Variable v(name);
  v.v = MX::sym(name);
  v.variability = Variable::FIXED;
  v.causality = Variable::CALCULATED_PARAMETER;
  v.beq = new_ddef;
  d_.push_back(add_variable(name, v));
  return v.v;
}

MX DaeBuilderInternal::add_w(const std::string& name, const MX& new_wdef) {
  Variable v(name);
  v.v = MX::sym(name);
  v.variability = Variable::CONTINUOUS;
  v.beq = new_wdef;
  w_.push_back(add_variable(name, v));
  return v.v;
}

MX DaeBuilderInternal::add_y(const std::string& name, const MX& new_ydef) {
  Variable v(name);
  v.v = MX::sym(name);
  v.causality = Variable::OUTPUT;
  v.beq = new_ydef;
  y_.push_back(add_variable(name, v));
  return v.v;
}

MX DaeBuilderInternal::add_ode(const std::string& name, const MX& new_ode) {
  Variable v(name);
  v.v = MX::sym(name);
  v.causality = Variable::OUTPUT;
  v.beq = new_ode;
  ode_.push_back(add_variable(name, v));
  return v.v;
}

MX DaeBuilderInternal::add_alg(const std::string& name, const MX& new_alg) {
  Variable v(name);
  v.v = MX::sym(name);
  v.causality = Variable::OUTPUT;
  v.beq = new_alg;
  alg_.push_back(add_variable(name, v));
  return v.v;
}

MX DaeBuilderInternal::add_quad(const std::string& name, const MX& new_quad) {
  Variable v(name);
  v.v = MX::sym(name);
  v.causality = Variable::OUTPUT;
  v.beq = new_quad;
  quad_.push_back(add_variable(name, v));
  return v.v;
}

template<typename T>
std::vector<T> read_list(const XmlNode& n) {
  // Number of elements
  size_t sz = n.size();
  // Read the elements
  std::vector<T> r;
  r.reserve(sz);
  for (size_t i = 0; i < sz; ++i) {
    r.push_back(T(n[i]));
  }
  return r;
}

void DaeBuilderInternal::import_model_exchange(const XmlNode& n) {
  // Read attributes
  provides_directional_derivative_
    = n.attribute<bool>("providesDirectionalDerivative", false);
  model_identifier_ = n.attribute<std::string>("modelIdentifier");
  // Get list of source files
  if (n.has_child("SourceFiles")) {
    for (const XmlNode& sf : n["SourceFiles"].children) {
      source_files_.push_back(sf.attribute<std::string>("name"));
    }
  }
}

void DaeBuilderInternal::import_model_variables(const XmlNode& modvars) {
  // Add variables
  for (casadi_int i = 0; i < modvars.size(); ++i) {
    // Get a reference to the variable
    const XmlNode& vnode = modvars[i];

    // Name of variable
    std::string name = vnode.attribute<std::string>("name");

    // Ignore duplicate variables
    if (varind_.find(name) != varind_.end()) {
      casadi_warning("Duplicate variable '" + name + "' ignored");
      continue;
    }

    // Create new variable
    Variable var(name);
    var.v = MX::sym(name);

    // Read common attributes, cf. FMI 2.0.2 specification, 2.2.7
    var.value_reference = vnode.attribute<casadi_int>("valueReference");
    var.description = vnode.attribute<std::string>("description", "");
    std::string causality_str = vnode.attribute<std::string>("causality", "local");
    if (causality_str == "internal") causality_str = "local";  // FMI 1.0 -> FMI 2.0
    var.causality = to_enum<Variable::Causality>(causality_str);
    std::string variability_str = vnode.attribute<std::string>("variability", "continuous");
    if (variability_str == "parameter") variability_str = "fixed";  // FMI 1.0 -> FMI 2.0
    var.variability = to_enum<Variable::Variability>(variability_str);
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
    (void)add_variable(name, var);
  }
  // Handle derivatives/antiderivatives
  for (auto it = variables_.begin(); it != variables_.end(); ++it) {
    if (it->derivative >= 0) {
      // Add variable offset, make index 1
      it->derivative -= 1;
      // Set antiderivative
      variables_.at(it->derivative).antiderivative = it - variables_.begin();
    }
  }
}

void DaeBuilderInternal::import_model_structure(const XmlNode& n) {
  // Outputs
  if (n.has_child("Outputs")) {
    for (auto& e : n["Outputs"].children) {
      // Get index
      outputs_.push_back(e.attribute<casadi_int>("index", 0) - 1);
      // Get dependencies
      for (casadi_int d : e.attribute<std::vector<casadi_int>>("dependencies", {})) {
        variables_.at(d - 1).dependency = true;
      }
    }
  }
  // Derivatives
  if (n.has_child("Derivatives")) {
    for (auto& e : n["Derivatives"].children) {
      // Get index
      derivatives_.push_back(e.attribute<casadi_int>("index", 0) - 1);
      // Get dependencies
      for (casadi_int d : e.attribute<std::vector<casadi_int>>("dependencies", {})) {
        variables_.at(d - 1).dependency = true;
      }
    }
  }
  // Initial unknowns
  if (n.has_child("InitialUnknowns")) {
    for (auto& e : n["InitialUnknowns"].children) {
      // Get index
      initial_unknowns_.push_back(e.attribute<casadi_int>("index", 0) - 1);
      // Get dependencies
      for (casadi_int d : e.attribute<std::vector<casadi_int>>("dependencies", {})) {
        variables_.at(d - 1).dependency = true;
      }
    }
  }
}

const MX& DaeBuilderInternal::var(size_t ind) const {
  return variables_.at(ind).v;
}

std::vector<MX> DaeBuilderInternal::var(const std::vector<size_t>& ind) const {
  std::vector<MX> ret;
  ret.reserve(ind.size());
  for (size_t i : ind) ret.push_back(var(i));
  return ret;
}

size_t DaeBuilderInternal::find(const std::string& name) const {
  auto it = varind_.find(name);
  casadi_assert(it != varind_.end(), "No such variable: \"" + name + "\".");
  return it->second;
}

Function DaeBuilderInternal::add_fun(const Function& f) {
  casadi_assert(!has_fun(f.name()), "Function '" + f.name() + "' already exists");
  fun_.push_back(f);
  return f;
}

Function DaeBuilderInternal::add_fun(const std::string& name,
    const std::vector<std::string>& arg,
    const std::vector<std::string>& res, const Dict& opts) {
  casadi_assert(!has_fun(name), "Function '" + name + "' already exists");

  // Dependent variable definitions
  std::vector<MX> wdef = this->wdef();
  // Get inputs
  std::vector<MX> arg_ex, res_ex;
  for (auto&& s : arg) arg_ex.push_back(var(s));
  for (auto&& s : res) {
    // Find the binding expression FIXME(@jaeandersson)
    casadi_int v_ind;
    for (v_ind = 0; v_ind < w_.size(); ++v_ind) {
      if (s == variable(w_.at(v_ind)).name) {
        res_ex.push_back(wdef.at(v_ind));
        break;
      }
    }
    casadi_assert(v_ind < w_.size(), "Cannot find dependent '" + s + "'");
  }
  Function ret(name, arg_ex, res_ex, arg, res, opts);
  return add_fun(ret);
}

bool DaeBuilderInternal::has_fun(const std::string& name) const {
  for (const Function& f : fun_) {
    if (f.name()==name) return true;
  }
  return false;
}

Function DaeBuilderInternal::fun(const std::string& name) const {
  casadi_assert(has_fun(name), "No such function: '" + name + "'");
  for (const Function& f : fun_) {
    if (f.name()==name) return f;
  }
  return Function();
}

std::string DaeBuilderInternal::system_infix() {
#if defined(_WIN32)
  // Windows system
#ifdef _WIN64
  return "win64";
#else
  return "win32";
#endif
#elif defined(__APPLE__)
  // OSX
  return sizeof(void*) == 4 ? "darwin32" : "darwin64";
#else
  // Linux
  return sizeof(void*) == 4 ? "linux32" : "linux64";
#endif
}

std::string DaeBuilderInternal::dll_suffix() {
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

void DaeBuilderInternal::init_fmu() const {
#ifdef WITH_FMU
  // Free existing instance, if any
  if (fmu_) {
    delete fmu_;
    fmu_ = 0;
  }
  // Allocate new instance
  fmu_ = new Fmu(*this);
  // Directory where the DLL is stored, per the FMI specification
  std::string instance_name_no_dot = model_identifier_;
  std::replace(instance_name_no_dot.begin(), instance_name_no_dot.end(), '.', '_');
  std::string dll_path = path_ + "/binaries/" + system_infix()
    + "/" + instance_name_no_dot + dll_suffix();
  // Load DLL
  fmu_->li_ = Importer(dll_path, "dll");
  // Load functions
  fmu_->init();
#else  // WITH_FMU
  casadi_error("FMU support not enabled. Recompile CasADi with 'WITH_FMU=ON'");
#endif  // WITH_FMU
}

#ifdef WITH_FMU

Fmu::Fmu(const DaeBuilderInternal& self) : self_(self) {
  // Initialize to null pointers
  instantiate_ = 0;
  free_instance_ = 0;
  reset_ = 0;
  setup_experiment_ = 0;
  enter_initialization_mode_ = 0;
  exit_initialization_mode_ = 0;
  enter_continuous_time_mode_ = 0;
  set_real_ = 0;
  set_boolean_ = 0;
  get_real_ = 0;
  get_directional_derivative_ = 0;
}

void Fmu::init() {
  instantiate_ = reinterpret_cast<fmi2InstantiateTYPE*>(get_function("fmi2Instantiate"));
  free_instance_ = reinterpret_cast<fmi2FreeInstanceTYPE*>(get_function("fmi2FreeInstance"));
  reset_ = reinterpret_cast<fmi2ResetTYPE*>(get_function("fmi2Reset"));
  setup_experiment_ = reinterpret_cast<fmi2SetupExperimentTYPE*>(
    get_function("fmi2SetupExperiment"));
  enter_initialization_mode_ = reinterpret_cast<fmi2EnterInitializationModeTYPE*>(
    get_function("fmi2EnterInitializationMode"));
  exit_initialization_mode_ = reinterpret_cast<fmi2ExitInitializationModeTYPE*>(
    get_function("fmi2ExitInitializationMode"));
  enter_continuous_time_mode_ = reinterpret_cast<fmi2EnterContinuousTimeModeTYPE*>(
    get_function("fmi2EnterContinuousTimeMode"));
  set_real_ = reinterpret_cast<fmi2SetRealTYPE*>(get_function("fmi2SetReal"));
  set_boolean_ = reinterpret_cast<fmi2SetBooleanTYPE*>(get_function("fmi2SetBoolean"));
  get_real_ = reinterpret_cast<fmi2GetRealTYPE*>(get_function("fmi2GetReal"));
  if (self_.provides_directional_derivative_) {
    get_directional_derivative_ = reinterpret_cast<fmi2GetDirectionalDerivativeTYPE*>(
      get_function("fmi2GetDirectionalDerivative"));
  }

  // Callback functions
  functions_.logger = logger;
  functions_.allocateMemory = calloc;
  functions_.freeMemory = free;
  functions_.stepFinished = 0;
  functions_.componentEnvironment = 0;

  // Path to resource directory
  resource_loc_ = "file://" + self_.path_ + "/resources";
}

signal_t Fmu::get_function(const std::string& symname) {
  // Load the function
  signal_t f = li_.get_function(symname);
  // Ensure that it was found
  casadi_assert(f != 0, "Cannot retrieve '" + symname + "'");
  // Return function to be type converted
  return f;
}

void Fmu::logger(fmi2ComponentEnvironment componentEnvironment,
    fmi2String instanceName,
    fmi2Status status,
    fmi2String category,
    fmi2String message, ...) {
  // Variable number of arguments
  va_list args;
  va_start(args, message);
  // Static & dynamic buffers
  char buf[256];
  size_t buf_sz = sizeof(buf);
  char* buf_dyn = nullptr;
  // Try to print with a small buffer
  int n = vsnprintf(buf, buf_sz, message, args);
  // Need a larger buffer?
  if (n > buf_sz) {
    buf_sz = n + 1;
    buf_dyn = new char[buf_sz];
    n = vsnprintf(buf_dyn, buf_sz, message, args);
  }
  // Print buffer content
  if (n >= 0) {
    uout() << "[" << instanceName << ":" << category << "] "
      << (buf_dyn ? buf_dyn : buf) << std::endl;
  }
  // Cleanup
  delete[] buf_dyn;
  va_end(args);
  // Throw error if failure
  casadi_assert(n>=0, "Print failure while processing '" + std::string(message) + "'");
}

int Fmu::instantiate() {
  // Create instance
  fmi2String instanceName = self_.model_identifier_.c_str();
  fmi2Type fmuType = fmi2ModelExchange;
  fmi2String fmuGUID = self_.guid_.c_str();
  fmi2String fmuResourceLocation = resource_loc_.c_str();
  fmi2Boolean visible = fmi2False;
  fmi2Boolean loggingOn = fmi2False;
  fmi2Component c = instantiate_(instanceName, fmuType, fmuGUID, fmuResourceLocation,
    &functions_, visible, loggingOn);
  if (c == 0) casadi_error("fmi2Instantiate failed");
  // Reuse an element of the memory pool
  for (int ind = 0; ind < mem_.size(); ++ind) {
    if (mem_[ind] == 0) {
      // Reuse
      mem_[ind] = c;
      return ind;
    }
  }
  // New element in the memory pool
  mem_.push_back(c);
  return mem_.size() - 1;
}

fmi2Component Fmu::mem(int ind) {
  return mem_.at(ind);
}

void Fmu::free_instance(int ind) {
  if (ind >= 0) {
    casadi_assert_dev(free_instance_ != 0);
    // Remove component from memory pool
    fmi2Component c = mem_.at(ind);
    mem_.at(ind) = 0;
    // Free memory
    free_instance_(c);
  }
}

int Fmu::setup_experiment(int ind) {
  fmi2Status status = setup_experiment_(mem(ind), fmi2False, 0.0, 0., fmi2True, 1.);
  if (status != fmi2OK) {
    casadi_warning("fmi2SetupExperiment failed");
    return 1;
  }
  return 0;
}

int Fmu::reset(int ind) {
  fmi2Status status = reset_(mem(ind));
  if (status != fmi2OK) {
    casadi_warning("fmi2Reset failed");
    return 1;
  }
  return 0;
}

#endif  // WITH_FMU

} // namespace casadi

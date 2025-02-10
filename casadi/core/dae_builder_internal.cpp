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
#include "fmu_function.hpp"
#include "integrator.hpp"

// Throw informative error message
#define THROW_ERROR_NODE(FNAME, NODE, WHAT) \
throw CasadiException("Error in DaeBuilderInternal::" FNAME " for '" + this->name_ \
  + "', node '" + NODE.name + "' (line " + str(NODE.line) + ") at " \
  + CASADI_WHERE + ":\n" + std::string(WHAT));

namespace casadi {

Type from_fmi2(TypeFmi2 v) {
  switch (v) {
  case TypeFmi2::REAL: return Type::FLOAT64;
  case TypeFmi2::INTEGER: return Type::INT32;
  case TypeFmi2::BOOLEAN: return Type::BOOLEAN;
  case TypeFmi2::STRING: return Type::STRING;
  case TypeFmi2::ENUM: return Type::ENUMERATION;
  default: break;
  }
  return Type::NUMEL;
}

TypeFmi2 to_fmi2(Type v) {
  switch (v) {
  case Type::FLOAT64: return TypeFmi2::REAL;
  case Type::INT32: return TypeFmi2::INTEGER;
  case Type::BOOLEAN: return TypeFmi2::BOOLEAN;
  case Type::STRING: return TypeFmi2::STRING;
  case Type::ENUMERATION: return TypeFmi2::ENUM;
  case Type::FLOAT32:  // fall-through
  case Type::INT8:  // fall-through
  case Type::UINT8:  // fall-through
  case Type::INT16:  // fall-through
  case Type::UINT16:  // fall-through
  case Type::UINT32:  // fall-through
  case Type::INT64:  // fall-through
  case Type::UINT64:  // fall-through
  case Type::BINARY:  // fall-through
  case Type::CLOCK:  // fall-through
    casadi_error(to_string(v) + " cannot be converted to FMI 2");
  default: break;
  }
  return TypeFmi2::NUMEL;
}

std::string to_string(TypeFmi2 v) {
  switch (v) {
  case TypeFmi2::REAL: return "real";
  case TypeFmi2::INTEGER: return "integer";
  case TypeFmi2::BOOLEAN: return "boolean";
  case TypeFmi2::STRING: return "string";
  case TypeFmi2::ENUM: return "enum";
  default: break;
  }
  return "";
}

std::string to_string(Type v) {
  switch (v) {
  case Type::FLOAT32: return "Float32";
  case Type::FLOAT64: return "Float64";
  case Type::INT8: return "Int8";
  case Type::UINT8: return "UInt8";
  case Type::INT16: return "Int16";
  case Type::UINT16: return "UInt16";
  case Type::INT32: return "Int32";
  case Type::UINT32: return "UInt32";
  case Type::INT64: return "Int64";
  case Type::UINT64: return "UInt64";
  case Type::BOOLEAN: return "Boolean";
  case Type::STRING: return "String";
  case Type::BINARY: return "Binary";
  case Type::ENUMERATION: return "Enumeration";
  case Type::CLOCK: return "Clock";
  default: break;
  }
  return "";
}

std::string to_string(Causality v) {
  switch (v) {
  case Causality::PARAMETER: return "parameter";
  case Causality::CALCULATED_PARAMETER: return "calculatedParameter";
  case Causality::INPUT: return "input";
  case Causality::OUTPUT: return "output";
  case Causality::LOCAL: return "local";
  case Causality::INDEPENDENT: return "independent";
  default: break;
  }
  return "";
}

std::string to_string(Variability v) {
  switch (v) {
  case Variability::CONSTANT: return "constant";
  case Variability::FIXED: return "fixed";
  case Variability::TUNABLE: return "tunable";
  case Variability::DISCRETE: return "discrete";
  case Variability::CONTINUOUS: return "continuous";
  default: break;
  }
  return "";
}

std::string to_string(Initial v) {
  switch (v) {
  case Initial::EXACT: return "exact";
  case Initial::APPROX: return "approx";
  case Initial::CALCULATED: return "calculated";
  case Initial::NA: return "na";
  default: break;
  }
  return "";
}

CASADI_EXPORT std::string to_string(Attribute v) {
  switch (v) {
  case Attribute::MIN: return "min";
  case Attribute::MAX: return "max";
  case Attribute::NOMINAL: return "nominal";
  case Attribute::START: return "start";
  case Attribute::VALUE: return "value";
  case Attribute::STRINGVALUE: return "stringvalue";
  default: break;
  }
  return "";
}

CASADI_EXPORT std::string to_string(DependenciesKind v) {
  switch (v) {
  case DependenciesKind::DEPENDENT: return "dependent";
  case DependenciesKind::CONSTANT: return "constant";
  case DependenciesKind::FIXED: return "fixed";
  case DependenciesKind::TUNABLE: return "tunable";
  case DependenciesKind::DISCRETE: return "discrete";
  default: break;
  }
  return "";
}

std::string to_string(Category v) {
  switch (v) {
  case Category::T: return "t";
  case Category::P: return "p";
  case Category::U: return "u";
  case Category::X: return "x";
  case Category::Z: return "z";
  case Category::Q: return "q";
  case Category::C: return "c";
  case Category::D: return "d";
  case Category::W: return "w";
  case Category::Y: return "y";
  case Category::E: return "e";
  default: break;
  }
  return "";
}


double Variable::attribute(Attribute a) const {
  switch (a) {
    case Attribute::MIN:
      return min;
    case Attribute::MAX:
      return max;
    case Attribute::NOMINAL:
      return nominal;
    case Attribute::START:
      casadi_assert(numel == 1, "Not a scalar variable");
      return start.front();
    case Attribute::VALUE:
      casadi_assert(numel == 1, "Not a scalar variable");
      return value.front();
    default:
      break;
  }
  casadi_error("Cannot handle: " + to_string(a));
  return 0;
}

void Variable::set_attribute(Attribute a, double val) {
  switch (a) {
    case Attribute::MIN:
      min = val;
      return;
    case Attribute::MAX:
      max = val;
      return;
    case Attribute::NOMINAL:
      nominal = val;
      return;
    case Attribute::START:
      std::fill(start.begin(), start.end(), val);
      return;
    case Attribute::VALUE:
      std::fill(value.begin(), value.end(), val);
      return;
    default:
      break;
  }
  casadi_error("Cannot handle: " + to_string(a));
}

std::string Variable::string_attribute(Attribute a) const {
  switch (a) {
    case Attribute::STRINGVALUE:
      return stringvalue;
    default:
      break;
  }
  casadi_error("Cannot handle: " + to_string(a));
  return std::string();
}

void Variable::set_string_attribute(Attribute a, const std::string& val) {
  switch (a) {
    case Attribute::STRINGVALUE:
      stringvalue = val;
      return;
    default:
      break;
  }
}

Initial Variable::default_initial(Causality causality, Variability variability) {
  // According to table in FMI 2.0.2 specification, section 2.2.7
  switch (variability) {
  case Variability::CONSTANT:
    if (causality == Causality::OUTPUT || causality == Causality::LOCAL)
      return Initial::EXACT;
    break;
  case Variability::FIXED:
    // Fall-through
  case Variability::TUNABLE:
    if (causality == Causality::PARAMETER)
      return Initial::EXACT;
    else if (causality == Causality::CALCULATED_PARAMETER || causality == Causality::LOCAL)
      return Initial::CALCULATED;
    break;
  case Variability::DISCRETE:
  // Fall-through
  case Variability::CONTINUOUS:
    if (causality == Causality::OUTPUT || causality == Causality::LOCAL)
      return Initial::CALCULATED;
    break;
  default: break;
  }
  // Initial value not available
  return Initial::NA;
}

Variable::Variable(casadi_int index, const std::string& name,
    const std::vector<casadi_int>& dimension)
    : index(index), name(name), dimension(dimension) {
  // Consistency checks
  casadi_assert(dimension.size() > 0, "Variable must have at least one dimension");
  for (casadi_int d : dimension) casadi_assert(d > 0, "Dimensions must be positive");
  // Number of elements
  numel = 1;
  for (casadi_int d : dimension) numel *= d;
  // Symbolic expression (always flattened)
  this->v = MX::sym(name, numel);
  // Default arguments
  this->value_reference = index;
  this->type = Type::FLOAT64;
  this->causality = Causality::LOCAL;
  this->variability = Variability::CONTINUOUS;
  this->category = Category::NUMEL;
  this->min = -inf;
  this->max = inf,
  this->nominal = 1.0;
  this->start.resize(numel, 0.0);
  this->der_of = -1;
  this->der = -1;
  this->alg = -1;
  this->value.resize(numel, nan);
  this->dependency = false;
}

XmlNode Variable::export_xml(const DaeBuilderInternal& self) const {
  // Create new XmlNode
  XmlNode r;
  r.name = to_string(type);
  // Name of variable
  r.set_attribute("name", name);
  // Value reference
  r.set_attribute("valueReference", static_cast<casadi_int>(value_reference));
  // Description, if any
  if (!description.empty()) r.set_attribute("description", description);
  // Causality
  if (causality != Causality::LOCAL)
    r.set_attribute("causality", to_string(causality));
  // Variability (real variables are continuous by default)
  if (!(is_real() && variability == Variability::CONTINUOUS))
    r.set_attribute("variability", to_string(variability));
  // Minimum attribute
  if (min != -inf) {
    if (is_real()) {
      r.set_attribute("min", min);
    } else {
      r.set_attribute("min", static_cast<casadi_int>(min));
    }
  }
  // Maximum attribute
  if (max != inf) {
    if (is_real()) {
      r.set_attribute("max", max);
    } else {
      r.set_attribute("max", static_cast<casadi_int>(max));
    }
  }
  // Unit
  if (!unit.empty()) r.set_attribute("unit", unit);
  // Display unit
  if (!display_unit.empty()) r.set_attribute("displayUnit", display_unit);
  // Nominal value, only for floats
  if (is_real() && nominal != 1.) r.set_attribute("nominal", nominal);
  // Start attribute, if any
  if (has_start()) {
    if (type == Type::BINARY || type == Type::STRING) {
      casadi_warning("Start attribute for String, Binary not implemented.");
    } else {
      // Convert to string
      std::stringstream ss;
      for (size_t i = 0; i < start.size(); ++i) {
        if (i > 0) ss << " ";
        if (is_real()) {
          ss << start.at(i);
        } else {
          ss << static_cast<casadi_int>(start.at(i));
        }
      }
      r.set_attribute("start", ss.str());
    }
  }
  // Derivative attribute, if any
  if (der_of >= 0) {
      r.set_attribute("derivative",
        static_cast<casadi_int>(self.variable(der_of).value_reference));
  }
  // Return XML representation
  return r;
}


bool Variable::has_start() const {
  // Rules, according to the FMI 3.0 specification, Section 2.4.7.5.
  if (initial == Initial::EXACT || initial == Initial::APPROX) return true;
  if (initial == Initial::CALCULATED || causality == Causality::INDEPENDENT) return false;
  if (causality == Causality::PARAMETER) return true;
  if (causality == Causality::INPUT) return true;
  if (variability == Variability::CONSTANT) return true;
  return false;
}

bool Variable::needs_der() const {
  // Only continuous variables can have derivatives
  if (variability != Variability::CONTINUOUS) return false;
  // Independent variables have trivial derivatives (1)
  if (causality == Causality::INDEPENDENT) return false;
  // Inputs are assumed piecewise constant
  if (causality == Causality::INPUT) return false;
  // Other variables may have derivatives
  return true;
}

MX Variable::dot(const DaeBuilderInternal& self) const {
  if (causality == Causality::INDEPENDENT) {
    // Time derivative of independent variable is 1
    return 1;
  } else if (needs_der()) {
    // Should have a derviative expression
    return self.variable(der).v;
  } else {
    // Constant or piecewise constant variable
    return MX::zeros(v.sparsity());
  }
}

DaeBuilderInternal::~DaeBuilderInternal() {
  for (Variable* v : variables_) {
    if (v) delete v;
  }
}

DaeBuilderInternal::DaeBuilderInternal(const std::string& name, const std::string& path,
    const Dict& opts) : name_(name), path_(path) {
  clear_cache_ = false;
  number_of_event_indicators_ = 0;
  provides_directional_derivatives_ = false;
  provides_adjoint_derivatives_ = false;
  can_be_instantiated_only_once_per_process_ = false;
  symbolic_ = true;
  start_time_ = nan;
  stop_time_ = nan;
  tolerance_ = nan;
  step_size_ = nan;
  // Default options
  debug_ = false;
  fmutol_ = 0;
  ignore_time_ = false;
  // Read options
  for (auto&& op : opts) {
    if (op.first=="debug") {
      debug_ = op.second;
    } else if (op.first=="fmutol") {
      fmutol_ = op.second;
    } else if (op.first=="ignore_time") {
      ignore_time_ = op.second;
    } else {
      casadi_error("No such option: " + op.first);
    }
  }
}

void DaeBuilderInternal::load_fmi_description(const std::string& filename) {
  // Ensure no variables already
  casadi_assert(n_variables() == 0, "Instance already has variables");

  // Parse XML file
  XmlFile xml_file("tinyxml");
  XmlNode fmi_desc = xml_file.parse(filename)[0];  // One child; fmiModelDescription

  // Read FMU version
  fmi_version_ = fmi_desc.attribute<std::string>("fmiVersion", "");
  if (fmi_version_.rfind("1.", 0) == 0) {
    // FMI 1.0: Only for symbolic FMUs
    fmi_major_ = 1;
  } else if (fmi_version_.rfind("2.", 0) == 0) {
    // FMI 2.0
    fmi_major_ = 2;
  } else if (fmi_version_.rfind("3.", 0) == 0) {
    // FMI 3
    fmi_major_ = 3;
  } else {
    casadi_error("Unknown FMI version: " + fmi_version_);
  }

  // Read attributes
  model_name_ = fmi_desc.attribute<std::string>("modelName", "");
  instantiation_token_ = fmi_desc.attribute<std::string>(
    fmi_major_ >= 3 ? "instantiationToken" : "guid");
  description_ = fmi_desc.attribute<std::string>("description", "");
  author_ = fmi_desc.attribute<std::string>("author", "");
  copyright_ = fmi_desc.attribute<std::string>("copyright", "");
  license_ = fmi_desc.attribute<std::string>("license", "");
  generation_tool_ = fmi_desc.attribute<std::string>("generationTool", "");
  generation_date_and_time_ = fmi_desc.attribute<std::string>("generationDateAndTime", "");
  variable_naming_convention_ = fmi_desc.attribute<std::string>("variableNamingConvention", "");
  if (fmi_major_ >= 3) {
    number_of_event_indicators_ = 0;  // In FMI 3: Obtain from binary
  } else {
    number_of_event_indicators_ = fmi_desc.attribute<casadi_int>("numberOfEventIndicators", 0);
  }

  // Process DefaultExperiment
  if (fmi_desc.has_child("DefaultExperiment")) {
    import_default_experiment(fmi_desc["DefaultExperiment"]);
  }

  // Process ModelExchange
  bool has_model_exchange = false;
  if (fmi_desc.has_child("ModelExchange")) {
    has_model_exchange = true;
    import_model_exchange(fmi_desc["ModelExchange"]);
  }

  // Process ModelVariables
  casadi_assert(fmi_desc.has_child("ModelVariables"), "Missing 'ModelVariables'");
  import_model_variables(fmi_desc["ModelVariables"]);

  // Process model structure
  if (fmi_desc.has_child("ModelStructure")) {
    import_model_structure(fmi_desc["ModelStructure"]);
  }

  // Is a symbolic representation available?
  symbolic_ = false;  // use DLL by default

  // Add symbolic binding equations
  if (fmi_desc.has_child("equ:BindingEquations")) {
    symbolic_ = true;
    import_binding_equations(fmi_desc["equ:BindingEquations"]);
  }

  // Add symbolic initial equations
  if (fmi_desc.has_child("equ:InitialEquations")) {
    symbolic_ = true;
    import_initial_equations(fmi_desc["equ:InitialEquations"]);
  }

  // Add symbolic dynamic equations
  if (fmi_desc.has_child("equ:DynamicEquations")) {
    symbolic_ = true;
    import_dynamic_equations(fmi_desc["equ:DynamicEquations"]);
  }

  // Ensure model equations are available, binary or symbolic
  casadi_assert(symbolic_ || has_model_exchange,
    "FMU must be of ModelExchange type or be symbolic (FMUX)");
}

std::string DaeBuilderInternal::generate_build_description(
    const std::vector<std::string>& cfiles) const {
  // Default arguments
  int fmi_major = 3;
  int fmi_minor = 0;
  std::string model_name = name_;
  // Construct XML file
  XmlNode r;
  // Preamble
  r.name = "fmiBuildDescription";
  r.set_attribute("fmiVersion", std::to_string(fmi_major) + "." + std::to_string(fmi_minor));
  // Set of source files
  XmlNode source_file_set;
  source_file_set.name = "SourceFileSet";
  for (auto&& f : cfiles) {
    XmlNode source_file;
    source_file.name = "SourceFile";
    source_file.set_attribute("name", f);
    source_file_set.children.push_back(source_file);
  }
  // Build configurations
  XmlNode bc;
  bc.name = "BuildConfiguration";
  bc.set_attribute("modelIdentifier", model_name);
  bc.children.push_back(source_file_set);
  r.children.push_back(bc);
  // XML file name
  std::string xml_filename = "buildDescription.xml";
  // Construct ModelDescription
  XmlNode build_description;
  build_description.children.push_back(r);
  // Export to file
  XmlFile xml_file("tinyxml");
  xml_file.dump(xml_filename, build_description);
  return xml_filename;
}

std::string DaeBuilderInternal::generate_model_description(const std::string& guid) const {
  // Default arguments
  int fmi_major = 3;
  int fmi_minor = 0;
  std::string model_name = name_;
  std::string description;  // none
  std::string author;  // none
  std::string version;  // none
  std::string copyright;  // none
  std::string license;  // none
  // Construct XML file
  XmlNode r;
  // Preamble
  r.name = "fmiModelDescription";
  r.set_attribute("fmiVersion", std::to_string(fmi_major) + "." + std::to_string(fmi_minor));
  r.set_attribute("modelName", model_name);
  r.set_attribute(fmi_major >= 3 ? "instantiationToken" : "guid", guid);
  if (!description.empty()) r.set_attribute("description", description);
  if (!author.empty()) r.set_attribute("author", author);
  if (!version.empty()) r.set_attribute("version", version);
  if (!copyright.empty()) r.set_attribute("copyright", copyright);
  if (!license.empty()) r.set_attribute("license", license);
  r.set_attribute("generationTool", "CasADi");
  r.set_attribute("generationDateAndTime", iso_8601_time());
  r.set_attribute("variableNamingConvention", "structured");  // flat better?
  if (fmi_major < 3) r.set_attribute("numberOfEventIndicators", "0");
  // Model exchange marker
  XmlNode me;
  me.name = "ModelExchange";
  me.set_attribute("modelIdentifier", model_name);  // sanitize name?
  r.children.push_back(me);
  // Model variables
  r.children.push_back(generate_model_variables());
  // Model structure
  r.children.push_back(generate_model_structure());
  // XML file name
  std::string xml_filename = "modelDescription.xml";
  // Construct ModelDescription
  XmlNode model_description;
  model_description.children.push_back(r);
  // Export to file
  XmlFile xml_file("tinyxml");
  xml_file.dump(xml_filename, model_description);
  return xml_filename;
}


XmlNode DaeBuilderInternal::generate_model_variables() const {
  XmlNode r;
  r.name = "ModelVariables";
  for (auto&& v : variables_) {
    r.children.push_back(v->export_xml(*this));
  }
  return r;
}

XmlNode DaeBuilderInternal::generate_model_structure() const {
  XmlNode r;
  r.name = "ModelStructure";
  // Add outputs
  for (size_t i : y_) {
    const Variable& y = variable(i);
    XmlNode c;
    c.name = "Output";
    c.set_attribute("valueReference", static_cast<casadi_int>(y.value_reference));
    c.set_attribute("dependencies", y.dependencies);
 r.children.push_back(c);
  }
  // Add state derivatives
  for (size_t i : x_) {
    const Variable& xdot = variable(variable(i).der);
    XmlNode c;
    c.name = "ContinuousStateDerivative";
    c.set_attribute("valueReference", static_cast<casadi_int>(xdot.value_reference));
    c.set_attribute("dependencies", xdot.dependencies);
    r.children.push_back(c);
  }
  // Add initial unknowns: Outputs
  for (size_t i : y_) {
    const Variable& y = variable(i);
    XmlNode c;
    c.name = "InitialUnknown";
    c.set_attribute("valueReference", static_cast<casadi_int>(y.value_reference));
    c.set_attribute("dependencies", y.dependencies);
    r.children.push_back(c);
  }
  // Add initial unknowns: State derivative
  for (size_t i : x_) {
    const Variable& xdot = variable(variable(i).der);
    XmlNode c;
    c.name = "InitialUnknown";
    c.set_attribute("valueReference", static_cast<casadi_int>(xdot.value_reference));
    c.set_attribute("dependencies", xdot.dependencies);
    r.children.push_back(c);
  }
  return r;
}

void DaeBuilderInternal::update_dependencies() const {
  // Get oracle function
  const Function& oracle = this->oracle();
  // Dependendencies of the ODE right-hand-side
  Sparsity dode_dxT = oracle.jac_sparsity(oracle.index_out("ode"), oracle.index_in("x")).T();
  Sparsity dode_duT = oracle.jac_sparsity(oracle.index_out("ode"), oracle.index_in("u")).T();
  for (casadi_int i = 0; i < x_.size(); ++i) {
    // Get output variable
    const Variable& xdot = variable(variable(x_.at(i)).der);
    // Clear dependencies
    xdot.dependencies.clear();
    // Dependencies on states
    for (casadi_int k = dode_dxT.colind(i); k < dode_dxT.colind(i + 1); ++k) {
      casadi_int j = dode_dxT.row(k);
      xdot.dependencies.push_back(variable(x_.at(j)).value_reference);
    }
    // Dependencies on controls
    for (casadi_int k = dode_duT.colind(i); k < dode_duT.colind(i + 1); ++k) {
      casadi_int j = dode_duT.row(k);
      xdot.dependencies.push_back(variable(u_.at(j)).value_reference);
    }
  }
  // Dependendencies of the output function
  Sparsity dydef_dxT = oracle.jac_sparsity(oracle.index_out("ydef"), oracle.index_in("x")).T();
  Sparsity dydef_duT = oracle.jac_sparsity(oracle.index_out("ydef"), oracle.index_in("u")).T();
  for (casadi_int i = 0; i < y_.size(); ++i) {
    // Get output variable
    const Variable& y = variable(y_.at(i));
    // Clear dependencies
    y.dependencies.clear();
    // Dependencies on states
    for (casadi_int k = dydef_dxT.colind(i); k < dydef_dxT.colind(i + 1); ++k) {
      casadi_int j = dydef_dxT.row(k);
      y.dependencies.push_back(variable(x_.at(j)).value_reference);
    }
    // Dependencies on controls
    for (casadi_int k = dydef_duT.colind(i); k < dydef_duT.colind(i + 1); ++k) {
      casadi_int j = dydef_duT.row(k);
      y.dependencies.push_back(variable(u_.at(j)).value_reference);
    }
  }
}

std::vector<std::string> DaeBuilderInternal::export_fmu(const Dict& opts) const {
  // Default options
  bool no_warning = false;
  for (auto&& op : opts) {
    if (op.first == "no_warning") {
      no_warning = op.second;
    }
  }
  // Feature incomplete
  if (!no_warning) casadi_warning("FMU generation is experimental and incomplete")
  // Return object
  std::vector<std::string> ret;
  // GUID
  std::string guid = generate_guid();
  // Generate model function
  std::string dae_filename = name_;
  Function dae = shared_from_this<DaeBuilder>().create(dae_filename,
    {"t", "x", "p", "u"}, {"ode", "ydef"});
  // Generate C code for model equations
  Dict codegen_opts;
  codegen_opts["with_header"] = true;
  CodeGenerator gen(dae_filename, codegen_opts);
  gen.add(dae);
  gen.add(dae.forward(1));
  gen.add(dae.reverse(1));
  ret.push_back(gen.generate());
  ret.push_back(dae_filename + ".h");
  // Make sure dependencies are up-to-date
  update_dependencies();
  // Generate FMU wrapper file
  std::string wrapper_filename = generate_wrapper(guid, gen);
  ret.push_back(wrapper_filename);
  // Generate build description
  ret.push_back(generate_build_description(ret));
  // Generate modelDescription file
  ret.push_back(generate_model_description(guid));
  // Return list of files
  return ret;
}

std::string DaeBuilderInternal::generate(const std::vector<size_t>& v) {
  std::stringstream ss;
  ss << "{";
  bool first = true;
  for (double e : v) {
    // Separator
    if (!first) ss << ", ";
    first = false;
    // Print element
    ss << e;
  }
  ss << "}";
  return ss.str();
}

std::string DaeBuilderInternal::generate(const std::vector<double>& v) {
  std::stringstream ss;
  ss << "{";
  bool first = true;
  for (double e : v) {
    // Separator
    if (!first) ss << ", ";
    first = false;
    // Print element
    ss << std::scientific << std::setprecision(std::numeric_limits<double>::digits10 + 1) << e;
  }
  ss << "}";
  return ss.str();
}

std::vector<double> DaeBuilderInternal::start_all() const {
  std::vector<double> r;
  for (const Variable* v : variables_) {
    for (double s : v->start) r.push_back(s);
  }
  return r;
}

std::string DaeBuilderInternal::generate_wrapper(const std::string& guid,
    const CodeGenerator& gen) const {
  // Create file
  std::string wrapper_filename = name_ + "_wrap.c";
  std::ofstream f;
  CodeGenerator::file_open(f, wrapper_filename, false);

  // Add includes
  f << "#include <fmi3Functions.h>\n"
    << "#include \"" << name_ << ".h\"\n"
    << "\n";

  // Total number of variables
  f << "#define N_VAR " << n_variables() << "\n";

  // Memory size
  f << "#define SZ_MEM " << n_mem() << "\n";

  // Work vectors sizes
  size_t sz_arg, sz_res, sz_iw, sz_w;
  gen.sz_work(sz_arg, sz_res, sz_iw, sz_w);
  f << "#define SZ_ARG " << sz_arg << "\n"
    << "#define SZ_RES " << sz_res << "\n"
    << "#define SZ_IW " << sz_iw << "\n"
    << "#define SZ_W " << sz_w << "\n";

  // Memory offsets
  f << "const size_t var_offset[N_VAR + 1] = {0";
  size_t mem_ind = 0;
  for (const Variable* v : variables_) {
    mem_ind += v->numel;
    f << ", " << mem_ind;
  }
  f << "};\n\n";

  // Start attributes
  f << "casadi_real start[SZ_MEM] = " << generate(start_all()) << ";\n\n";

  // States
  f << "#define N_X " << x_.size() << "\n"
    << "fmi3ValueReference x_vr[N_X] = " << generate(x_) << ";\n"
    << "\n";

  // Controls
  f << "#define N_U " << u_.size() << "\n"
    << "fmi3ValueReference u_vr[N_U] = " << generate(u_) << ";\n"
    << "\n";

  // Parameters
  f << "#define N_P " << p_.size() << "\n"
    << "fmi3ValueReference p_vr[N_P] = " << generate(p_) << ";\n"
    << "\n";

  // State derivatives
  std::vector<size_t> xdot;
  for (size_t v : x_) xdot.push_back(variable(v).der);
  f << "fmi3ValueReference xdot_vr[N_X] = " << generate(xdot) << ";\n"
    << "\n";

  // Outputs
  f << "#define N_Y " << y_.size() << "\n"
    << "fmi3ValueReference y_vr[N_Y] = " << generate(y_) << ";\n"
    << "\n";

  // Memory structure
  f << CodeGenerator::fmu_helpers(name_);

  // Finalize file
  CodeGenerator::file_close(f, false);
  return wrapper_filename;
}

Variable& DaeBuilderInternal::read_variable(const XmlNode& node, Attribute* att) {
  try {
    // Qualified name
    std::string qn = qualified_name(node, att);

    return variable(qn);
  } catch (std::exception& e) {
    THROW_ERROR_NODE("read_variable", node, e.what());
    //return {};
  }
}

MX DaeBuilderInternal::read_identifier(const XmlNode& node) {
  Attribute att;  // attribute
  Variable& v = read_variable(node, &att);
  if (att == Attribute::VALUE) {
    return v.v;
  } else if (att == Attribute::START) {
    return v.start;
  } else {
    casadi_error("Cannot read attribute " + to_string(att));
    return MX();
  }
}

MX DaeBuilderInternal::read_expr(const XmlNode& node) {
  try {
    const std::string& fullname = node.name;

    if (fullname == "fun:If") {
      // Expression for condition, if and else
      const XmlNode& cond = node["fun:Condition"];
      const XmlNode& then_stmt = node["fun:Statements"];
      const XmlNode& else_stmt = node["fun:Else"];
      // Only scalar expression implemented
      casadi_assert(cond.size() == 1, "Only one condition in if expression supported");
      casadi_assert(then_stmt.size() == 1, "Only one then statement in if expression supported");
      casadi_assert(else_stmt.size() == 1, "Only one else statement in if expression supported");
      return if_else(read_expr(cond[0]), read_expr(then_stmt[0]), read_expr(else_stmt[0]));
    }

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
      return variable(read_variable(node[0]).der_of).v;
    } else if (name=="Div") {
      return read_expr(node[0]) / read_expr(node[1]);
    } else if (name=="Exp") {
      return exp(read_expr(node[0]));
    } else if (name=="Identifier") {
      return read_identifier(node);
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
    } else if (name=="Not") { // Logical not
      return !read_expr(node[0]);
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
    } else if (name=="Pre") {
      casadi_warning("Ignoring pre attribute");
      return read_expr(node[0]);
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
        Variable& v = new_variable("w_" + str(w_.size()));
        // Add to list of variables
        w_.push_back(v.index);
        // Set binding expression
        v.beq = read_expr(args[i]);
        // Add to list of function arguments
        farg[i] = v.v;
      }
      // Return argument (scalar for now)
      Variable& r = new_variable("w_" + str(w_.size()));
      // Add to list of variables
      w_.push_back(r.index);
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
  } catch (std::exception& e) {
    THROW_ERROR_NODE("read_expr", node, e.what());
    return {};
  }
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
    stream << "Differential equations" << std::endl;
    for (size_t k : x_) {
      const Variable& x = variable(k);
      stream << "  \\dot{" << x.name << "} == ";
      if (x.der >= 0) {
        // Derivative variable
        const Variable& xdot = variable(x.der);
        stream << xdot.name;
        if (!xdot.beq.is_empty()) stream << " := " << xdot.beq;
      } else if (x.variability == Variability::DISCRETE) {
        // Discrete variable - derivative is zero
        stream << 0;
      } else {
        // Missing ODE?
        casadi_error("Missing derivative for " + str(x.name));
      }
      stream << std::endl;
    }
  }

  if (!z_.empty()) {
    stream << "Algebraic equations" << std::endl;
    for (size_t k : z_) {
      const Variable& z = variable(k);
      casadi_assert(z.alg >= 0, "No residual variable for " + z.name);
      const Variable& alg = variable(z.alg);
      stream << "  0 == " << alg.beq;
      if (!alg.beq.is_empty()) stream << " := " << alg.beq;
      stream << std::endl;
    }
  }

  if (!q_.empty()) {
    stream << "Quadrature equations" << std::endl;
    for (size_t k : q_) {
      const Variable& q = variable(k);
      casadi_assert(q.der >= 0, "No derivative variable for " + q.name);
      const Variable& qdot = variable(q.der);
      stream << "  \\dot{" << q.name << "} == " << qdot.name;
      if (!qdot.beq.is_empty()) stream << " := " << qdot.beq;
      stream << std::endl;
    }
  }

  if (!init_.empty()) {
    stream << "Initial equations" << std::endl;
    for (size_t k : init_) {
      const Variable& v = variable(k);
      stream << "  " << v.name;
      if (!v.ieq.is_empty()) stream << " := " << v.ieq;
      stream << std::endl;
    }
  }

  if (!when_cond_.empty()) {
    stream << "When statements" << std::endl;
    for (casadi_int k = 0; k < when_cond_.size(); ++k) {
      stream << "  when " << str(when_cond_.at(k)) << " > 0 : " << str(when_lhs_.at(k))
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
  std::vector<bool> old_z(n_variables(), false);
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

std::vector<size_t>& DaeBuilderInternal::ind_in(const std::string& v) {
  switch (to_enum<Category>(v)) {
  case Category::T: return t_;
  case Category::P: return p_;
  case Category::U: return u_;
  case Category::X: return x_;
  case Category::Z: return z_;
  case Category::Q: return q_;
  case Category::C: return c_;
  case Category::D: return d_;
  case Category::W: return w_;
  case Category::Y: return y_;
  default: break;
  }
  // Unsuccessful
  casadi_error("Cannot access input indices for " + v);
  // Just to resolve warnings
  static std::vector<size_t> dummy;
  return dummy;
}

const std::vector<size_t>& DaeBuilderInternal::ind_in(const std::string& v) const {
  return const_cast<DaeBuilderInternal*>(this)->ind_in(v);
}

void DaeBuilderInternal::clear_all(const std::string& v) {
  clear_cache_ = true;  // Clear cache after this
  ind_in(v).clear();
}

void DaeBuilderInternal::set_all(const std::string& v, const std::vector<std::string>& name) {
  clear_cache_ = true;  // Clear cache after this
  ind_in(v) = find(name);
}

void DaeBuilderInternal::prune(bool prune_p, bool prune_u) {
  // Function inputs and outputs
  std::vector<MX> f_in, f_out, v;
  std::vector<std::string> f_in_name, f_out_name;
  // Collect all DAE input variables with at least one entry, skip u
  for (casadi_int i = 0; i != enum_traits<Category>::n_enum; ++i) {
    auto cat = static_cast<Category>(i);
    if (prune_p && cat == Category::P) continue;
    if (prune_u && cat == Category::U) continue;
    v = input(cat);
    if (!v.empty()) {
      f_in.push_back(vertcat(v));
      f_in_name.push_back(to_string(cat));
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
  std::vector<bool> free_variables(n_variables(), false);
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
  for (size_t k = 0; k < z_.size(); ++k) {
    if (!iv_set.count(variable(z_[k]).name)) {
      // Non-iteration variable: Keep
      z_.at(k) = z_.at(sz);
      sz++;
    }
  }
  z_.resize(sz);
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
  for (const Variable* v : variables_) {
    // Residual variables are specified with a "res__" prefix
    if (v->name.rfind(res_prefix, 0) == 0) {
      // Process iteration variable name, names of hold markers
      std::string iv_name, res_hold_name, iv_hold_name;
      // Iteration variable, hold markers are contained in the remainder of the name
      try {
        size_t pos = res_prefix.size();
        // Find the next "__", if any
        size_t end = v->name.find("__", pos);
        if (end == std::string::npos) end = v->name.size();
        // Look up iteration variable
        iv_name = v->name.substr(pos, end - pos);
        // Ensure that the variable exists
        casadi_assert(has(iv_name), "No such variable: " + iv_name);
        // Get hold indices, if any
        if (end != v->name.size()) {
          // Find next "__", read hold index for residual variable
          pos = end + 2;
          end = v->name.find("__", pos);
          if (end == std::string::npos) end = v->name.size();
          res_hold_name = v->name.substr(pos, end - pos);
          // Ensure that the variable exists
          casadi_assert(has(res_hold_name), "No such variable: " + res_hold_name);
          // The remainder of the name contains iv_hold_name
          if (end != v->name.size()) {
            iv_hold_name = v->name.substr(end + 2);
            casadi_assert(has(iv_hold_name), "No such variable: " + iv_hold_name);
          }
        }
      } catch (std::exception& e) {
        // Generate warning
        casadi_warning("Cannot process residual variable: " + v->name + ":" +
                       std::string(e.what()));
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
      if (res) res->push_back(v->name);
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
      // Code below needs to be refactored
      casadi_error("not implemented");
#if 0
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
#endif
    } catch (std::exception& e) {
      // Warning instead of error
      casadi_warning("Failed to evaluate hold variables: " + std::string(e.what()));
    }
  }
}

bool DaeBuilderInternal::has(const std::string& name) const {
  return varind_.find(name) != varind_.end();
}

std::vector<std::string> DaeBuilderInternal::all() const {
  std::vector<std::string> r;
  r.reserve(n_variables());
  for (const Variable* v : variables_) r.push_back(v->name);
  return r;
}

std::vector<std::string> DaeBuilderInternal::all(const std::string& cat) const {
  return name(ind_in(cat));
}

size_t DaeBuilderInternal::n_mem() const {
  size_t n = 0;
  for (const Variable* v : variables_) n += v->numel;
  return n;
}

Variable& DaeBuilderInternal::new_variable(const std::string& name,
    const std::vector<casadi_int>& dimension) {
  // Name check
  casadi_assert(!name.empty(), "Name is empty string");
  // Try to find the component
  casadi_assert(!has(name), "Variable \"" + name + "\" already exists.");
  // Index of the variable
  size_t ind = n_variables();
  // Add to the map of all variables
  varind_[name] = ind;
  variables_.push_back(new Variable(ind, name, dimension));
  // Clear cache
  clear_cache_ = true;
  // Return reference to new variable
  return *variables_.back();
}

void DaeBuilderInternal::sanity_check() const {
  // Time
  if (!t_.empty()) {
    casadi_assert(t_.size() == 1, "At most one time variable allowed");
    casadi_assert(var(t_[0]).is_scalar(), "Non-scalar time t");
  }

  // When statements
  casadi_assert(when_cond_.size() == when_lhs_.size() && when_lhs_.size() == when_rhs_.size(),
    "when_cond, when_lhs and when_rhs must all have the the same length");
}

std::string DaeBuilderInternal::qualified_name(const XmlNode& nn, Attribute* att) {
  // std::stringstream to assemble name
  std::stringstream qn;
  bool first_part = true;
  if (att) *att = Attribute::VALUE;  // value attribute by default

  // Loop over name parts
  for (casadi_int i=0; i<nn.size(); ++i) {
    // Get the name part
    std::string np = nn[i].attribute<std::string>("name");

    // Check if an attribute
    if (np == "$START") {
      if (att) {
        *att = Attribute::START;
      } else {
        casadi_error("Ignoring attribute " + np);
      }
      continue;
    } else if (np == "$PRE") {
      casadi_warning("$PRE attribute has not been implemented, ignoring");
      continue;
    }

    // Add the name part to the variable name
    if (!first_part) qn << ".";
    qn << np;

    // Get the index, if any
    if (nn[i].size()>0) {
      casadi_int ind;
      nn[i]["exp:ArraySubscripts"]["exp:IndexExpression"]["exp:IntegerLiteral"].get(&ind);
      qn << "[" << ind << "]";
    }

    // Dot prefix if a additional parts
    first_part = false;
  }

  // Return the name
  return qn.str();
}

const MX& DaeBuilderInternal::var(const std::string& name) const {
  return variable(name).v;
}

MX DaeBuilderInternal::der(const MX& var) const {
  // Must be a vector
  casadi_assert(var.is_column(), "Input expression must be a vector");
  // Quick return if symbolic variable
  if (var.is_symbolic()) return der(find(var));
  // If a vertical concatenation
  if (var.is_valid_input()) {
    // Differentiate each primitive
    auto var_split = var.primitives();
    for (MX& s : var_split) s = der(s);
    // Return the concatenation
    return var.join_primitives(var_split);
  }
  // Handle general case: Get dependent symbolic primitives
  std::vector<MX> dep = symvar(var);
  // Get derivatives of dependent symbolic primitives
  std::vector<MX> dep_der;
  for (size_t ind : find(dep)) dep_der.push_back(der(ind));
  // Forward directional derivative to get time derivative:
  // dot(var) = d_var/d_dep * dot(dep)
  std::vector<std::vector<MX>> r = {dep_der};
  r = forward(std::vector<MX>{var}, dep, r);
  casadi_assert_dev(r.size() == 1);
  return vertcat(r.at(0));
}

void DaeBuilderInternal::eliminate_d() {
  // Quick return if no d
  if (d_.empty()) return;
  // Clear cache after this
  clear_cache_ = true;
  // Ensure variables are sorted
  sort_d();
  // Expressions where the variables are also being used
  std::vector<MX> ex;
  for (const Variable* v : variables_) {
    if (!v->beq.is_constant()) ex.push_back(v->beq);
  }
  // Perform elimination
  std::vector<MX> d = var(d_);
  std::vector<MX> ddef = this->ddef();
  substitute_inplace(d, ddef, ex);
  // Clear list of dependent parameters
  d_.clear();
  // Get binding equations
  auto it = ex.begin();
  for (Variable* v : variables_) {
    if (!v->beq.is_constant()) v->beq = *it++;
  }
  // Consistency check
  casadi_assert_dev(it == ex.end());
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
  for (const Variable* v : variables_) {
    if (!v->beq.is_constant()) ex.push_back(v->beq);
  }
  // Perform elimination
  std::vector<MX> w = var(w_);
  std::vector<MX> wdef = this->wdef();
  substitute_inplace(w, wdef, ex);
  // Clear list of dependent variables
  w_.clear();
  // Get binding equations
  auto it = ex.begin();
  for (Variable* v : variables_) {
    if (!v->beq.is_constant()) v->beq = *it++;
  }
  // Consistency check
  casadi_assert_dev(it == ex.end());
}

void DaeBuilderInternal::lift(bool lift_shared, bool lift_calls) {
  // Not tested if w is non-empty before
  if (!w_.empty()) casadi_warning("'w' already has entries");
  // Expressions where the variables are also being used
  std::vector<MX> ex;
  for (size_t v : x_) ex.push_back(variable(variable(v).der).beq);
  for (size_t v : q_) ex.push_back(variable(variable(v).der).beq);
  for (size_t v : z_) ex.push_back(variable(variable(v).alg).beq);
  for (size_t v : y_) ex.push_back(variable(v).beq);
  // Lift expressions
  std::vector<MX> new_w, new_wdef;
  Dict opts{{"lift_shared", lift_shared}, {"lift_calls", lift_calls},
    {"prefix", "w_"}, {"suffix", ""}, {"offset", static_cast<casadi_int>(w_.size())}};
  extract(ex, new_w, new_wdef, opts);
  // Register as dependent variables
  for (size_t i = 0; i < new_w.size(); ++i) {
    Variable& v = new_variable(new_w.at(i).name());
    v.v = new_w.at(i);
    v.beq = new_wdef.at(i);
    w_.push_back(v.index);
  }
  // Get expressions
  auto it = ex.begin();
  for (size_t v : x_) variable(variable(v).der).beq = *it++;
  for (size_t v : q_) variable(variable(v).der).beq = *it++;
  for (size_t v : z_) variable(variable(v).alg).beq = *it++;
  for (size_t v : y_) variable(v).beq = *it++;
  // Consistency check
  casadi_assert_dev(it == ex.end());
}

std::string to_string(DaeBuilderInternal::DaeBuilderInternalOut v) {
  switch (v) {
  case DaeBuilderInternal::DAE_BUILDER_ODE: return "ode";
  case DaeBuilderInternal::DAE_BUILDER_ALG: return "alg";
  case DaeBuilderInternal::DAE_BUILDER_QUAD: return "quad";
  case DaeBuilderInternal::DAE_BUILDER_ZERO: return "zero";
  case DaeBuilderInternal::DAE_BUILDER_DDEF: return "ddef";
  case DaeBuilderInternal::DAE_BUILDER_WDEF: return "wdef";
  case DaeBuilderInternal::DAE_BUILDER_YDEF: return "ydef";
  default: break;
  }
  return "";
}

std::vector<MX> DaeBuilderInternal::input(Category ind) const {
  switch (ind) {
  case Category::T: return var(t_);
  case Category::C: return var(c_);
  case Category::P: return var(p_);
  case Category::D: return var(d_);
  case Category::W: return var(w_);
  case Category::U: return var(u_);
  case Category::X: return var(x_);
  case Category::Z: return var(z_);
  case Category::Q: return var(q_);
  case Category::Y: return var(y_);
  case Category::E: return var(e_);
  default: return std::vector<MX>{};
  }
}

std::vector<MX> DaeBuilderInternal::input(const std::vector<Category>& ind) const {
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
  case DAE_BUILDER_ZERO: return zero();
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
    const std::vector<std::string>& s_out, const Dict& opts, bool sx, bool lifted_calls) const {
  // Are there any '_' in the names?
  bool with_underscore = false;
  for (auto s_io : {&s_in, &s_out}) {
    for (const std::string& s : *s_io) {
      with_underscore = with_underscore || std::count(s.begin(), s.end(), '_');
    }
  }
  // Model equations in DLL
  if (!symbolic_) {
    // Cannot lift calls in an FMU
    casadi_assert(!lifted_calls, "Lifting requires a symbolic representation");
    // Cannot convert to SX
    casadi_assert(!sx, "SX expansion requires a symbolic representation");
    // Redirect to FmuFunction creation
    return fmu_fun(fname, s_in, s_out, opts);
  }
  // Replace '_' with ':', if needed
  if (with_underscore) {
    std::vector<std::string> s_in_mod(s_in), s_out_mod(s_out);
    for (auto s_io : {&s_in_mod, &s_out_mod}) {
      for (std::string& s : *s_io) std::replace(s.begin(), s.end(), '_', ':');
    }
    // Recursive call
    return create(fname, s_in_mod, s_out_mod, opts, sx, lifted_calls);
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
    // Collect all DAE input variables
    for (size_t i = 0; i != enum_traits<Category>::n_enum; ++i) {
      auto cat = static_cast<Category>(i);
      if (cat == Category::Y) continue;  // fixme2
      f_in_name.push_back(to_string(cat));
      v = input(cat);
      if (v.empty()) {
        f_in.push_back(MX(0, 1));
      } else {
        if (elim_w && cat == Category::W) {
          subst_v = true;
        } else {
          f_in.push_back(vertcat(v));
        }
      }
    }
    // Collect all DAE output variables
    for (size_t i = 0; i != DAE_BUILDER_NUM_OUT; ++i) {
      f_out_name.push_back(to_string(static_cast<DaeBuilderInternalOut>(i)));
      v = output(static_cast<DaeBuilderInternalOut>(i));
      if (v.empty()) {
        f_out.push_back(MX(0, 1));
      } else {
        if (i == DAE_BUILDER_WDEF) wdef_ind = f_out.size();
        f_out.push_back(vertcat(v));
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
  Function vfcn("vfcn", {vertcat(v)}, {vertcat(vdef)}, {"v"}, {"vdef"},
    Dict{{"allow_free", true}});
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
      vfcn.name_in(), vfcn.name_out(), {{"allow_free", true}});
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

Function DaeBuilderInternal::dependent_fun(const std::string& fname,
    const std::vector<std::string>& s_in,
    const std::vector<std::string>& s_out) const {
  // Are we calculating d and/or w
  bool calc_d = false, calc_w = false;
  // Convert outputs to enums
  std::vector<Category> v_out;
  v_out.reserve(v_out.size());
  for (const std::string& s : s_out) {
    Category e = to_enum<Category>(s);
    if (e == Category::D) {
      calc_d = true;
    } else if (e == Category::W) {
      calc_w = true;
    } else {
      casadi_error("Can only calculate d and/or w");
    }
    v_out.push_back(e);
  }
  // Consistency check
  casadi_assert(calc_d || calc_w, "Nothing to calculate");
  // Convert inputs to enums
  std::vector<Category> v_in;
  v_in.reserve(v_in.size());
  for (const std::string& s : s_in) {
    Category e = to_enum<Category>(s);
    if (calc_d && e == Category::D) casadi_error("'d' cannot be both input and output");
    if (calc_w && e == Category::W) casadi_error("'w' cannot be both input and output");
    v_in.push_back(e);
  }
  // Collect input expressions
  std::vector<MX> f_in;
  f_in.reserve(s_in.size());
  for (Category v : v_in) f_in.push_back(vertcat(input(v)));
  // Collect output expressions
  std::vector<MX> f_out;
  f_out.reserve(s_out.size());
  for (Category v : v_out) f_out.push_back(vertcat(input(v)));
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

Function DaeBuilderInternal::event_transition(const std::string& fname, casadi_int index,
    bool dummy_index_input) const {

  // Make sure that the index is valid
  casadi_assert(index >= 0 && index < when_cond_.size(), "Illegal event index");

  // Get input expressions for the oracle
  std::vector<MX> oracle_in = oracle().mx_in();

  // Input expressions for the event functions, without the index
  std::vector<MX> ret_in(DYN_NUM_IN);
  ret_in[DYN_T] = oracle_in.at(static_cast<size_t>(Category::T));
  ret_in[DYN_X] = oracle_in.at(static_cast<size_t>(Category::X));
  ret_in[DYN_Z] = oracle_in.at(static_cast<size_t>(Category::Z));
  ret_in[DYN_P] = oracle_in.at(static_cast<size_t>(Category::P));
  ret_in[DYN_U] = oracle_in.at(static_cast<size_t>(Category::U));

  // Expressions for x and z after event
  std::vector<MX> ret_out = {ret_in[DYN_X], ret_in[DYN_Z]};
  ret_out = MX::substitute(ret_out,
    std::vector<MX>{when_lhs_.at(index)},
    std::vector<MX>{when_rhs_.at(index)});

  // Check if a dummy index input needes to be included
  if (dummy_index_input) {
    // Create a function with the event_transition input signature
    ret_in.insert(ret_in.begin(), MX());
    return Function(fname, ret_in, ret_out, event_in(), event_out());
  } else {
    // Create a function with the DAE function input signature
    return Function(fname, ret_in, ret_out, dyn_in(), event_out());
  }
}

Function DaeBuilderInternal::event_transition(const std::string& fname) const {
  // If no events, return null
  if (when_cond_.empty()) return Function();

  // If just a single event, create an event function with a dummy index input
  if (when_cond_.size() == 1) return event_transition(fname, 0, true);

  // Create separate transition functions for each event
  std::vector<Function> f_all;
  for (casadi_int i = 0; i < when_cond_.size(); ++i) {
    f_all.push_back(event_transition(fname + "_" + str(i), i));
  }

  // Make the last function the default value in the switch
  Function f_def = f_all.back();
  f_all.pop_back();

  // Create a switch function
  return Function::conditional(fname, f_all, f_def);
}


Function DaeBuilderInternal::fmu_fun(const std::string& name,
    const std::vector<std::string>& name_in,
    const std::vector<std::string>& name_out,
    const Dict& opts) const {
  // Iterator for options lookup
  Dict::const_iterator it;
  // Scheme inputs
  std::vector<std::string> scheme_in;
  bool has_in = false;
  if ((it = opts.find("scheme_in")) != opts.end()) {
    try {
      scheme_in = it->second;
    } catch (std::exception& e) {
      casadi_error(std::string("Cannot read 'scheme_in': ") + e.what());
    }
    has_in = true;
  }
  // Scheme outputs
  std::vector<std::string> scheme_out;
  bool has_out = false;
  if ((it = opts.find("scheme_out")) != opts.end()) {
    try {
      scheme_out = it->second;
      has_out = true;
    } catch (std::exception& e) {
      casadi_error(std::string("Cannot read 'scheme_out': ") + e.what());
    }
  }
  // If scheme_in and/or scheme_out not provided, identify from name_in, name_out
  if (!has_in || !has_out) {
    FmuFunction::identify_io(has_in ? 0 : &scheme_in, has_out ? 0 : &scheme_out, name_in, name_out);
  }
  // IO scheme
  std::map<std::string, std::vector<size_t>> scheme;
  if ((it = opts.find("scheme")) != opts.end()) {
    try {
      // Argument is a Dict
      Dict scheme_dict = it->second;
      // Convert indices
      for (auto&& e : scheme_dict) {
        std::vector<std::string> v = e.second;
        scheme[e.first] = find(v);
      }
    } catch (std::exception& e) {
      casadi_error(std::string("Cannot read 'scheme': ") + e.what());
    }
  } else {
    // Initialize all scheme entries
    for (auto&& s : dyn_in()) scheme[s] = std::vector<size_t>();
    for (auto&& s : dyn_out()) scheme[s] = std::vector<size_t>();
    // Default IO scheme
    scheme["t"] = ind_in("t");
    scheme["x"] = ind_in("x");
    scheme["u"] = ind_in("u");
    scheme["z"] = ind_in("z");
    scheme["p"] = ind_in("p");
    scheme["ode"] = x_;
    for (size_t& i : scheme["ode"]) i = variable(i).der;
    scheme["alg"] = z_;
    casadi_assert(z_.empty(), "Not implemented)");
    scheme["ydef"] = y_;
  }
  // Auxilliary variables, if any
  std::vector<std::string> aux;
  if ((it = opts.find("aux")) != opts.end()) {
    try {
      aux = it->second;
    } catch (std::exception& e) {
      casadi_error(std::string("Cannot read 'aux': ") + e.what());
    }
  }
  // New FMU instance (to be shared between derivative functions)
  Fmu fmu(name, fmi_major_ >= 3 ? FmuApi::FMI3 : FmuApi::FMI2, this,
    scheme_in, scheme_out, scheme, aux);

  // Crete new function
  return Function::create(new FmuFunction(name, fmu, name_in, name_out), opts);
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
  return Function("all_eq", {}, f_out, {}, f_out_name, {{"allow_free", true}});
}

const MX& DaeBuilderInternal::time() const {
  casadi_assert(has_t(), "No explicit time variable");
  return var(t_.at(0));
}

bool DaeBuilderInternal::has_t() const {
  return !t_.empty();
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
  ret.reserve(x_.size());
  for (size_t v : x_) {
    const Variable& x = variable(v);
    if (x.der >= 0) {
      // Derivative variable
      ret.push_back(variable(x.der).beq);
    } else if (x.variability == Variability::DISCRETE) {
      // Discrete variable - derivative is zero
      ret.push_back(MX::zeros(x.v.sparsity()));
    } else {
      // Missing ODE?
      casadi_error("Missing derivative for " + str(x.name));
    }
  }
  return ret;
}

std::vector<MX> DaeBuilderInternal::alg() const {
  std::vector<MX> ret;
  ret.reserve(z_.size());
  for (size_t v : z_) {
    const Variable& z = variable(v);
    casadi_assert(z.alg >= 0, "No residual variable for " + z.name);
    const Variable& alg = variable(z.alg);
    ret.push_back(alg.beq);
  }
  return ret;
}

std::vector<MX> DaeBuilderInternal::quad() const {
  std::vector<MX> ret;
  ret.reserve(q_.size());
  for (size_t v : q_) {
    const Variable& q = variable(v);
    casadi_assert(q.der >= 0, "No derivative variable for " + q.name);
    const Variable& qdot = variable(q.der);
    ret.push_back(qdot.beq);
  }
  return ret;
}


std::vector<MX> DaeBuilderInternal::zero() const {
  std::vector<MX> ret;
  ret.reserve(e_.size());
  for (size_t v : e_) ret.push_back(variable(v).beq);
  return ret;
}

std::vector<MX> DaeBuilderInternal::init_lhs() const {
  std::vector<MX> ret;
  ret.reserve(init_.size());
  for (size_t ind : init_) {
    ret.push_back(variable(ind).v);
  }
  return ret;
}

std::vector<MX> DaeBuilderInternal::init_rhs() const {
  std::vector<MX> ret;
  ret.reserve(init_.size());
  for (size_t ind : init_) {
    ret.push_back(variable(ind).ieq);
  }
  return ret;
}

MX DaeBuilderInternal::add(const std::string& name, Causality causality,
    Variability variability, const Dict& opts) {
  // Default options
  std::string cat;
  std::vector<casadi_int> dimension = {1};
  std::string description;
  casadi::MX def;
  // Read options
  for (auto&& op : opts) {
    if (op.first=="cat") {
      cat = op.second.to_string();
    } else if (op.first=="dimension") {
      dimension = op.second.to_int_vector();
    } else if (op.first=="description") {
      description = op.second.to_string();
    } else {
      casadi_error("No such option: " + op.first);
    }
  }
  // Create a new variable
  Variable& v = new_variable(name, dimension);
  v.description = description;
  // Handle different categories
  if (cat.empty()) {
    // Use provided causality and variability
    v.causality = causality;
    v.variability = variability;
    // Auto-detect category (default)
    switch (causality) {
      case Causality::INDEPENDENT:
        // Independent variable
        casadi_assert(t_.empty(), "'t' already defined");
        v.category = Category::T;
        t_.push_back(v.index);
        break;
      case Causality::INPUT:
        // Control
        v.category = Category::U;
        u_.push_back(v.index);
        break;
      case Causality::OUTPUT:
        // Output
        v.category = Category::Y;
        y_.push_back(v.index);
        break;
      case Causality::PARAMETER:
        // Parameter
        if (variability == Variability::TUNABLE) {
          v.category = Category::P;
          p_.push_back(v.index);
        }
        break;
      case Causality::CALCULATED_PARAMETER:
        // Calculated parameter: No categorization for now
        break;
      case Causality::LOCAL:
        // Type determined by providing equation
        break;
      default:
        casadi_error("Unknown causality");
    }
    // Also create a derivative variable, if needed
    if (v.needs_der()) {
      Variable& der_v = new_variable("der_" + name, dimension);
      der_v.der_of = v.index;
      v.der = der_v.index;
    }
  } else {
    // Explicit category (legacy, to be removed)
    switch (to_enum<Category>(cat)) {
      case Category::T:
        // Independent variable
        casadi_assert(t_.empty(), "'t' already defined");
        casadi_assert(dimension.size() == 1 && dimension[0] == 1,
          "Independent variable must be scalar");
        v.causality = Causality::INDEPENDENT;
        t_.push_back(v.index);
        break;
      case Category::P:
        // Parameter
        v.variability = Variability::FIXED;
        v.causality = Causality::INPUT;
        p_.push_back(v.index);
        break;
      case Category::U:
        // Control
        v.variability = Variability::CONTINUOUS;
        v.causality = Causality::INPUT;
        u_.push_back(v.index);
        break;
      case Category::X:
        // State
        v.variability = Variability::CONTINUOUS;
        v.causality = Causality::LOCAL;
        x_.push_back(v.index);
        break;
      case Category::Z:
        // Algebraic variable
        v.variability = Variability::CONTINUOUS;
        v.causality = Causality::LOCAL;
        z_.push_back(v.index);
        break;
      case Category::Q:
        // Quadrature variable
        v.variability = Variability::CONTINUOUS;
        v.causality = Causality::LOCAL;
        q_.push_back(v.index);
        break;
      case Category::C:
        // Constant
        v.variability = Variability::CONSTANT;
        c_.push_back(v.index);
        break;
      case Category::D:
        // Calculated parameter
        v.variability = Variability::FIXED;
        v.causality = Causality::CALCULATED_PARAMETER;
        d_.push_back(v.index);
        break;
      case Category::W:
        // Dependent variable
        v.variability = Variability::CONTINUOUS;
        v.causality = Causality::LOCAL;
        w_.push_back(v.index);
        break;
      case Category::Y:
        // Output
        v.causality = Causality::OUTPUT;
        y_.push_back(v.index);
        break;
      case Category::E:
        // Zero-crossing
        v.causality = Causality::OUTPUT;
        e_.push_back(v.index);
        break;
      default: break;
    }
  }
  // Return variable expression
  return v.v;
}

MX DaeBuilderInternal::add(const std::string& name, Causality causality, const Dict& opts) {
  // Default variability per FMI 3.0.2, section 2.4.7.4
  switch (causality) {
    // "The default for variables of causality parameter, structural parameter or
    // calculated parameter is fixed."
    case Causality::PARAMETER:  // fall-through
    case Causality::CALCULATED_PARAMETER:
      return add(name, causality, Variability::FIXED, opts);
    // "The default for variables of type Float32 and Float64 and causality other
    // than parameter, structuralParameter or calculatedParameter is continuous"
    case Causality::INPUT:  // fall-through
    case Causality::OUTPUT:  // fall-through
    case Causality::LOCAL:  // fall-through
    case Causality::INDEPENDENT:
      return add(name, causality, Variability::CONTINUOUS, opts);
    default:
      casadi_error("Unknown causality");
  }
}

void DaeBuilderInternal::eq(const MX& lhs, const MX& rhs) {
  // Make sure vectors
  casadi_assert(lhs.is_column(), "Left-hand-side must be a column vector");
  casadi_assert(rhs.is_column(), "Right-hand-side must be a column vector");
  // Make sure dense
  if (!lhs.is_dense()) return eq(densify(lhs), rhs);
  if (!rhs.is_dense()) return eq(lhs, densify(rhs));
  // Make sure dimensions agree
  if (lhs.size1() != rhs.size1()) {
    // Handle mismatching dimnensions by recursion
    if (lhs.size1() == 1 && rhs.size1() > 1) {
      return eq(repmat(lhs, rhs.size1()), rhs);
    } else if (lhs.size1() > 1 && rhs.size1() == 1) {
      return eq(lhs, repmat(rhs, lhs.size1()));
    } else {
      casadi_error("Mismatched dimensions: " + str(lhs.size1()) + " vs " + str(rhs.size1()));
    }
  }
  // Make sure right-hand-side only depends on known model variables
  std::vector<size_t> rhs_vars = find(symvar(rhs));
  // Try to honor a == b as an explicit equation
  if (lhs.is_valid_input()) {
    // Explicit equation
    if (lhs.is_symbolic()) {
      // Regular symbolic: Find the variable
      Variable& v = variable(lhs);
      // Set the binding equation
      if (v.has_beq()) {
        // TODO(@jaeandersson): Treat as implicit equation
        casadi_error("Already a binding equation for " + v.name);
      } else {
        // Set the binding equation
        v.beq = rhs;
      }
      // (Re)classify variables
      if (v.der_of >= 0) {
        // Derivative variable is being set - find the corresponding state variable
        Variable& x = variable(v.der_of);
        // Update category
        if (x.category == Category::NUMEL) {
          x.category = Category::X;
          // Insert into list of state variables: Keep ordered
          size_t loc = x_.size();
          for (size_t i = 0; i < x_.size(); ++i) {
            if (variable(x_[i]).index >= x.index) {
              loc = i;
              break;
            }
          }
          x_.insert(x_.begin() + loc, x.index);
        }
      } else {
        casadi_error("Only ODEs implemented");
      }
    } else {
      // Concatenation: Split into primitives
      auto lhs_split = lhs.primitives();
      std::vector<MX> rhs_split = lhs.split_primitives(rhs);
      // Call recursively
      for (size_t k = 0; k < lhs_split.size(); ++k) {
        eq(lhs_split.at(k), rhs_split.at(k));
      }
      return;
    }
  } else {
    // Implicit equation
    casadi_error("Implicit equations not yet supported");
  }
}

void DaeBuilderInternal::set_ode(const std::string& name, const MX& ode_rhs) {
  // Find the state variable
  const Variable& x = variable(name);
  // Check if derivative exists
  if (x.der < 0) {
    // New derivative variable
    Variable& xdot = new_variable("der_" + name);
    xdot.causality = Causality::LOCAL;
    xdot.der_of = find(name);
    xdot.beq = ode_rhs;
    variable(name).der = xdot.index;
  } else {
    // Variable exists: Update binding equation
    variable(x.der).beq = ode_rhs;
  }
}

void DaeBuilderInternal::set_alg(const std::string& name, const MX& alg_rhs) {
  // Find the algebraic variable
  const Variable& z = variable(name);
  // Check if residual exists
  if (z.alg < 0) {
    // New derivative variable
    Variable& alg = new_variable("alg_" + name);
    alg.causality = Causality::OUTPUT;
    alg.beq = alg_rhs;
    variable(name).alg = alg.index;
  } else {
    // Variable exists: Update binding equation
    variable(z.alg).beq = alg_rhs;
  }
}

void DaeBuilderInternal::set_init(const std::string& name, const MX& init_rhs) {
  // Find the algebraic variable
  Variable& v = variable(name);
  // If variable already has an initial binding equation, remove it
  if (!v.ieq.is_empty()) {
    // Remove from list of initial equations
    auto old_loc = std::find(init_.begin(), init_.end(), v.index);
    if (old_loc == init_.end()) casadi_error("Corrupted list of initial equations");
    init_.erase(old_loc);
    v.ieq = MX();
  }
  // If right-hand-side is empty, just erase
  if (init_rhs.is_empty()) return;

  // Make sure not already in list of initial equations
  if (std::find(init_.begin(), init_.end(), v.index) != init_.end()) {
    casadi_error("Initial equation for " + name + " has already been set");
  }
  // Add to list of initial equations
  init_.push_back(v.index);
  v.ieq = init_rhs;
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

void DaeBuilderInternal::import_default_experiment(const XmlNode& n) {
  start_time_ = n.attribute<double>("startTime", nan);
  stop_time_ = n.attribute<double>("stopTime", nan);
  tolerance_ = n.attribute<double>("tolerance", nan);
  step_size_ = n.attribute<double>("stepSize", nan);
}

void DaeBuilderInternal::import_model_exchange(const XmlNode& n) {
  // Read attributes
  provides_directional_derivatives_ = n.attribute<bool>(
    fmi_major_ >= 3 ? "providesDirectionalDerivatives" : "providesDirectionalDerivative", false);
  provides_adjoint_derivatives_
    = n.attribute<bool>("providesAdjointDerivatives", false);
  model_identifier_ = n.attribute<std::string>("modelIdentifier");
  can_be_instantiated_only_once_per_process_ =
    n.attribute<bool>("canBeInstantiatedOnlyOncePerProcess", false);
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
    Variable& var = new_variable(name);

    // Get type (FMI 3)
    if (fmi_major_ >= 3) var.type = to_enum<Type>(vnode.name);

    // Read common attributes, cf. FMI 3.0 specification, 2.4.7.4
    var.value_reference = static_cast<unsigned int>(vnode.attribute<casadi_int>("valueReference"));
    vrmap_[var.value_reference] = var.index;
    var.description = vnode.attribute<std::string>("description", "");
    std::string causality_str = vnode.attribute<std::string>("causality", "local");
    std::string variability_str = vnode.attribute<std::string>("variability", "continuous");
    if (fmi_major_ == 1) {
        // FMI 1.0 -> FMI 2.0+
      if (causality_str == "internal") causality_str = "local";
      if (variability_str == "parameter") variability_str = "fixed";
    }
    var.causality = to_enum<Causality>(causality_str);
    var.variability = to_enum<Variability>(variability_str);
    std::string initial_str = vnode.attribute<std::string>("initial", "");
    if (initial_str.empty()) {
      // Default value
      var.initial = Variable::default_initial(var.causality, var.variability);
    } else {
      // Consistency check
      casadi_assert(var.causality != Causality::INPUT && var.causality != Causality::INDEPENDENT,
        "The combination causality = '" + to_string(var.causality) + "', "
        "initial = '" + initial_str + "' is not allowed per the FMI specification.");
      // Value specified
      var.initial = to_enum<Initial>(initial_str);
    }
    // Type specific properties
    if (fmi_major_ >= 3) {
      // FMI 3.0: Type information in the same node
      switch (var.type) {
      case Type::FLOAT32:  // fall-through
      case Type::FLOAT64:
        // Floating point valued variables
        var.unit = vnode.attribute<std::string>("unit", var.unit);
        var.display_unit = vnode.attribute<std::string>("displayUnit", var.display_unit);
        var.min = vnode.attribute<double>("min", -inf);
        var.max = vnode.attribute<double>("max", inf);
        var.nominal = vnode.attribute<double>("nominal", 1.);
        var.set_attribute(Attribute::START, vnode.attribute<double>("start", 0.));
        var.der_of = vnode.attribute<casadi_int>("derivative", var.der_of);
        break;
      case Type::INT8:  // fall-through
      case Type::UINT8:  // fall-through
      case Type::INT16:  // fall-through
      case Type::UINT16:  // fall-through
      case Type::INT32:  // fall-through
      case Type::UINT32:  // fall-through
      case Type::INT64:  // fall-through
      case Type::UINT64:  // fall-through
        // Integer valued variables
        var.min = vnode.attribute<double>("min", -inf);
        var.max = vnode.attribute<double>("max", inf);
        break;
      default:
        break;
      }
    } else {
      // FMI 1.0 / 2.0: Type information in a separate node
      if (vnode.has_child("Real")) {
        const XmlNode& props = vnode["Real"];
        var.unit = props.attribute<std::string>("unit", var.unit);
        var.display_unit = props.attribute<std::string>("displayUnit", var.display_unit);
        var.min = props.attribute<double>("min", -inf);
        var.max = props.attribute<double>("max", inf);
        var.nominal = props.attribute<double>("nominal", 1.);
        var.set_attribute(Attribute::START, props.attribute<double>("start", 0.));
        var.der_of = props.attribute<casadi_int>("derivative", var.der_of);
      } else if (vnode.has_child("Integer")) {
        const XmlNode& props = vnode["Integer"];
        var.type = Type::INT32;
        var.min = props.attribute<double>("min", -inf);
        var.max = props.attribute<double>("max", inf);
      } else if (vnode.has_child("Boolean")) {
        var.type = Type::BOOLEAN;
      } else if (vnode.has_child("String")) {
        var.type = Type::STRING;
      } else if (vnode.has_child("Enumeration")) {
        var.type = Type::ENUMERATION;
      } else {
        casadi_warning("Unknown type for " + name);
      }
    }

    // Initial classification of variables (states/outputs to be added later)
    if (var.causality == Causality::INDEPENDENT) {
      // Independent (time) variable
      if (!ignore_time_) t_.push_back(var.index);
    } else if (var.causality == Causality::INPUT) {
      // Check if description starts with PARAMETER:
      if (var.description.rfind("PARAMETER:", 0) == 0) {
        // Make tunable parameter
        var.variability = Variability::TUNABLE;
        p_.push_back(var.index);
      } else {
        u_.push_back(var.index);
      }
    } else if (var.variability == Variability::TUNABLE) {
      p_.push_back(var.index);
    }
  }

  // Handle derivatives
  for (size_t i = 0; i < n_variables(); ++i) {
    if (variable(i).der_of >= 0) {
      if (fmi_major_ >= 3) {
        // Value reference is given: Find corresponding variable
        variable(i).der_of = vrmap_.at(static_cast<unsigned int>(variable(i).der_of));
      } else {
        // Variable given with index-1, make index 0
        variable(i).der_of -= 1;
      }
      // Set der
      variable(variable(i).der_of).der = i;
    }
  }
}

std::vector<casadi_int> DaeBuilderInternal::read_dependencies(const XmlNode& n) {
  // Default behaviour should be no known structure
  casadi_assert(n.has_attribute("dependencies"),
    "Default 'dependencies' not implemented");
  // Read list of dependencies
  std::vector<casadi_int> r = n.attribute<std::vector<casadi_int>>("dependencies", {});
  // Get corresponding variable index
  for (casadi_int& e : r) {
    if (fmi_major_ >= 3) {
      // Value reference is given
      e = vrmap_.at(static_cast<unsigned int>(e));
    } else {
      // Index-1 is given
      e--;
    }
  }
  // Return list of dependencies
  return r;
}

std::vector<DependenciesKind> DaeBuilderInternal::read_dependencies_kind(
    const XmlNode& n, size_t ndep) {
  // Quick return if node not provided
  if (!n.has_attribute("dependenciesKind")) {
    // No structure known, assume general dependency
    return std::vector<DependenciesKind>(ndep, DependenciesKind::DEPENDENT);
  } else {
    // Read list of strings
    auto dk_str = n.attribute<std::vector<std::string>>("dependenciesKind");
    // Make sure expected length
    casadi_assert(dk_str.size() == ndep, "Mismatching 'dependenciesKind'");
    // Convert to enums
    std::vector<DependenciesKind> r(ndep);
    for (size_t i = 0; i < ndep; ++i) {
      r[i] = to_enum<DependenciesKind>(dk_str[i]);
    }
    return r;
  }
}

void DaeBuilderInternal::import_model_structure(const XmlNode& n) {
  if (fmi_major_ >= 3) {
    // Loop over ModelStructure elements
    for (casadi_int i = 0; i < n.size(); ++i) {
      const XmlNode& e = n[i];
      // Get a reference to the variable
      if (e.name == "Output") {
        // Get index
        outputs_.push_back(vrmap_.at(e.attribute<size_t>("valueReference")));
        // Corresponding variable
        Variable& v = variable(outputs_.back());
        // Add to y, unless state
        if (v.der < 0) {
          y_.push_back(v.index);
          v.beq = v.v;
        }
        // Get dependencies
        v.dependencies = read_dependencies(e);
        v.dependenciesKind = read_dependencies_kind(e, v.dependencies.size());
        // Mark interdependencies
        for (casadi_int d : v.dependencies) variable(d).dependency = true;
      } else if (e.name == "ContinuousStateDerivative") {
        // Get index
        derivatives_.push_back(vrmap_.at(e.attribute<size_t>("valueReference")));
        // Corresponding variable
        Variable& v = variable(derivatives_.back());
        // Add to list of states
        casadi_assert(v.der_of >= 0, "Error processing derivative info for " + v.name);
        x_.push_back(v.der_of);
        // Make sure der field is consistent
        variable(x_.back()).der = derivatives_.back();
        // Get dependencies
        v.dependencies = read_dependencies(e);
        v.dependenciesKind = read_dependencies_kind(e, v.dependencies.size());
        // Mark interdependencies
        for (casadi_int d : v.dependencies) variable(d).dependency = true;
      } else if (e.name == "ClockedState") {
        // Clocked state
        casadi_message("ClockedState not implemented, ignoring");
      } else if (e.name == "InitialUnknown") {
        // Get index
        initial_unknowns_.push_back(vrmap_.at(e.attribute<size_t>("valueReference")));
        // Get dependencies
        for (casadi_int d : read_dependencies(e)) variable(d).dependency = true;
      } else if (e.name == "EventIndicator") {
        // Event indicator
        number_of_event_indicators_++;
      } else {
        // Unknown
        casadi_error("Unknown ModelStructure element: " + e.name);
      }
    }
  } else {
    // Derivatives
    if (n.has_child("Derivatives")) {
      for (auto& e : n["Derivatives"].children) {
        // Get index
        derivatives_.push_back(e.attribute<casadi_int>("index", 0) - 1);
        // Corresponding variable
        Variable& v = variable(derivatives_.back());
        // Add to list of states
        casadi_assert(v.der_of >= 0, "Error processing derivative info for " + v.name);
        x_.push_back(v.der_of);
        // Make sure der field is consistent
        variable(x_.back()).der = derivatives_.back();
      }
    }

    // What if dependencies attributed is missing from Outputs,Derivatives?
    // Depends on x_ having been populated
    std::vector<casadi_int> default_dependencies;
    default_dependencies.insert(default_dependencies.begin(), t_.begin(), t_.end());
    default_dependencies.insert(default_dependencies.begin(), x_.begin(), x_.end());
    default_dependencies.insert(default_dependencies.begin(), u_.begin(), u_.end());
    std::sort(default_dependencies.begin(), default_dependencies.end());

    // FMI2 standard does not allow parameters to be included in dependencies attribute
    //  (Table in section 2.2.8)
    // FMI3 standard does (Table in section 2.4.8)
    //
    // A uniform treatment is unconditionally add dense dependencies to all parameters
    // when loading FMI2

    std::vector<casadi_int> additional_dependencies;
    additional_dependencies.insert(additional_dependencies.begin(), p_.begin(), p_.end());
    std::vector<DependenciesKind> additional_dependencies_kind(p_.size(),
                                                               DependenciesKind::DEPENDENT);

    // Derivatives
    if (n.has_child("Derivatives")) {
      // Separate pass for dependencies
      for (auto& e : n["Derivatives"].children) {
        casadi_int index = e.attribute<casadi_int>("index", 0)-1;

        // Corresponding variable
        Variable& v = variable(index);

        // Get dependencies
        if (e.has_attribute("dependencies")) {
          v.dependencies = read_dependencies(e);
        } else {
          v.dependencies = default_dependencies;
        }
        v.dependenciesKind = read_dependencies_kind(e, v.dependencies.size());

        v.dependencies.insert(v.dependencies.end(),
          additional_dependencies.begin(), additional_dependencies.end());
        v.dependenciesKind.insert(v.dependenciesKind.end(),
          additional_dependencies_kind.begin(), additional_dependencies_kind.end());

        // Mark interdependencies
        for (casadi_int d : v.dependencies) variable(d).dependency = true;
      }
    }
    // Outputs
    if (n.has_child("Outputs")) {
      for (auto& e : n["Outputs"].children) {
        // Get index
        outputs_.push_back(e.attribute<casadi_int>("index", 0) - 1);
        // Corresponding variable
        Variable& v = variable(outputs_.back());
        // Add to y, unless state
        if (v.der < 0) {
          y_.push_back(v.index);
          v.beq = v.v;
        }
        // Get dependencies
        if (e.has_attribute("dependencies")) {
          v.dependencies = read_dependencies(e);
        } else {
          v.dependencies = default_dependencies;
        }
        v.dependenciesKind = read_dependencies_kind(e, v.dependencies.size());

        v.dependencies.insert(v.dependencies.end(),
          additional_dependencies.begin(), additional_dependencies.end());
        v.dependenciesKind.insert(v.dependenciesKind.end(),
          additional_dependencies_kind.begin(), additional_dependencies_kind.end());

        // Mark interdependencies
        for (casadi_int d : v.dependencies) variable(d).dependency = true;
      }
    }
    // What if dependencies is missing from InitialUnknowns?
    // Depends on x_ having been populated
    default_dependencies.clear();
    default_dependencies.insert(default_dependencies.begin(), t_.begin(), t_.end());
    for (const Variable* v : variables_) {
      if (v->initial == Initial::EXACT) default_dependencies.push_back(v->index);
    }
    default_dependencies.insert(default_dependencies.begin(), u_.begin(), u_.end());
    std::sort(default_dependencies.begin(), default_dependencies.end());

    // Initial unknowns
    if (n.has_child("InitialUnknowns")) {
      for (auto& e : n["InitialUnknowns"].children) {
        // Get index
        initial_unknowns_.push_back(e.attribute<casadi_int>("index", 0) - 1);

        std::vector<casadi_int> dependencies;
        // Get dependencies
        if (e.has_attribute("dependencies")) {
          dependencies = read_dependencies(e);
        } else {
          dependencies = default_dependencies;
        }

        dependencies.insert(dependencies.end(),
          additional_dependencies.begin(), additional_dependencies.end());

        // Get dependencies
        for (casadi_int d : dependencies) variable(d).dependency = true;
      }
    }
  }
}

void DaeBuilderInternal::import_binding_equations(const XmlNode& eqs) {
  // Loop over binding equations
  for (casadi_int i = 0; i < eqs.size(); ++i) {
    const XmlNode& eq = eqs[i];
    std::string eq_name = "beq_" + str(i);
    // Try to read the node
    try {
      // Get the variable and binding expression
      Variable& var = read_variable(eq[0]);
      if (eq[1].size() == 1) {
        // Regular expression
        var.beq = read_expr(eq[1][0]);
      } else {
        // OpenModelica 1.17 occationally generates integer values without type specifier (bug?)
        casadi_assert(eq[1].size() == 0, "Not implemented");
        casadi_int val;
        eq[1].get(&val);
        casadi_warning(var.name + " has binding equation without type specifier: " + str(val));
        var.beq = val;
      }
      // Add to list of dependent parameters
      d_.push_back(var.index);
    } catch (std::exception& e) {
      casadi_error("Failed to read " + eq_name + ":" + str(e.what()));
    }
  }
}

void DaeBuilderInternal::import_dynamic_equations(const XmlNode& eqs) {
  // Add discrete states to x_
  // Note: Make generic, should also apply to regular FMUs
  for (Variable* v : variables_) {
    if (v->variability == Variability::DISCRETE) x_.push_back(v->index);
  }
  // Add equations
  for (casadi_int i = 0; i < eqs.size(); ++i) {
    const XmlNode& eq = eqs[i];
    std::string eq_name = "dyneq_" + str(i);
    // Try to read the node
    try {
      if (eq.name == "equ:When") {  // When equation
        // Nodes for condition, equations
        const XmlNode& n_cond = eq["equ:Condition"];
        const XmlNode& n_equ = eq["equ:Equation"];
        // Consistency checks - only working for very simple expressions
        casadi_assert(n_cond.size() == 1, "Only one condition in when equation supported");
        casadi_assert(n_equ.size() == 1, "Only one equation in when equation supported");
        // Get expression for condition
        Variable& cond = read_variable(n_cond[0]);
        // Left-hand-side and right-hand-side
        MX lhs, rhs;
        // Handle different types of equations
        if (n_equ[0].name == "exp:Sub") {
          // Assume equation is an assignment
          lhs = read_expr(n_equ[0][0]);
          rhs = read_expr(n_equ[0][1]);
        } else if (n_equ[0].name == "exp:Reinit") {
          // Reinitialization
          lhs = read_identifier(n_equ[0][0]);
          rhs = read_expr(n_equ[0][1]);
        } else {
          // Not implemented
          casadi_error(n_equ[0].name + " in when equation not supported");
        }
        // Turn non-snooth zero-crossing expression into a smooth zero-crossing expression
        MX zc;
        switch (cond.beq.op()) {
          case OP_LT:
            // x1 < x2 <=> x2 - x1 > 0
            zc = cond.beq.dep(1) - cond.beq.dep(0);
            break;
          default:
            casadi_error("Cannot turn " + str(cond.beq) + " into a smooth expression");
        }
        set_init(cond.name, MX());  // remove initial conditions, if any
        auto w_it = std::find(w_.begin(), w_.end(), cond.index);
        if (w_it != w_.end()) w_.erase(w_it);  // remove from dependent equations
        auto x_it = std::find(x_.begin(), x_.end(), cond.index);
        if (x_it != x_.end()) x_.erase(x_it);  // remove from states

        // Create event indicator
        (void)add(cond.name + "_smooth", {{"cat", "e"}});
        variable(cond.name + "_smooth").beq = zc;
        // Add to list of when equations
        when_cond_.push_back(zc);
        when_lhs_.push_back(lhs);
        when_rhs_.push_back(rhs);
      } else if (eq.name == "equ:Equation") {  // Residual equation
        // Consistency checks
        casadi_assert_dev(eq.size() == 1 && eq[0].name == "exp:Sub");
        // Skip if empty
        if (eq[0].size() == 0) {
          casadi_warning(eq_name + " is empty, ignored.");
          continue;
        }
        // Get the left-hand-sides and right-hand-sides
        const XmlNode& lhs = eq[0][0];
        const XmlNode& rhs = eq[0][1];
        // Right-hand-side is the binding equation
        MX beq = read_expr(rhs);
        // Left-hand-side is a variable or derivative
        if (lhs.name == "exp:Der") {
          // Differentiated variable
          Variable& v = read_variable(lhs[0]);
          // Corresponding time derivative
          Variable& dot_v = variable("der(" + v.name + ")");
          // Map to each other
          v.der = dot_v.index;
          dot_v.der_of = v.index;
          // Mark as state
          x_.push_back(v.index);
          // Set binding equation to derivative variable
          dot_v.beq = beq;
        } else {
          // Left-hand-side is a variable
          Variable& v = read_variable(lhs);
          // Set the equation
          w_.push_back(find(v.name));
          v.beq = beq;
        }
      } else {
        casadi_error("Unknown dynamic equation type, got:" + eq.name);
      }
    } catch (std::exception& e) {
      casadi_error("Failed to read " + eq_name + ":" + str(e.what()));
    }
  }
}

void DaeBuilderInternal::import_initial_equations(const XmlNode& eqs) {
  // Hack: OpenModelica XML may contain duplicate expressions, ignore these
  std::set<std::string> already_added;
  // Add equations
  for (casadi_int i = 0; i < eqs.size(); ++i) {
    const XmlNode& eq = eqs[i];
    // Try to read the node
    std::string eq_name = "initeq_" + str(i);
    try {
      if (eq.name == "equ:Equation") {  // Residual equation
        // Consistency checks
        casadi_assert_dev(eq.size() == 1 && eq[0].name == "exp:Sub");
        // Ensure not empty
        if (eq[0].size() == 0) {
          casadi_warning(eq_name + " is empty, ignored.");
          continue;
        }
        // Get the left-hand-sides and right-hand-sides
        const XmlNode& lhs = eq[0][0];
        const XmlNode& rhs = eq[0][1];
        // Hack: Ignore expressions that just set a $PRE variable
        if (lhs.size() > 0 && lhs[0].attribute<std::string>("name") == "$PRE") {
          casadi_warning(eq_name + " defines a pre-variable, ignored");
          continue;
        }
        // Right-hand-side is the binding equation
        MX beq = read_expr(rhs);
        // Left-hand-side is a variable
        Variable& v = read_variable(lhs);
        // Hack: ignore initial equations for tunable parameters
        if (v.variability == Variability::TUNABLE) {
          casadi_warning(eq_name + " defines a tunable parameter, ignored")
          continue;
        }
        // Hack: avoid duplicate equations
        std::string eq_str = str(v.v) + " == " + str(beq);
        auto it = already_added.insert(eq_str);
        if (!it.second) {
          casadi_warning(eq_name + " duplicate of previous equation " + eq_str + ", ignored")
          continue;
        }
        // Set initial condition
        set_init(v.name, beq);
      } else {
        casadi_error("Unknown initial equation type, got:" + eq.name);
      }
    } catch (std::exception& e) {
      casadi_error("Failed to read " + eq_name + ":" + str(e.what()));
    }
  }
}

const MX& DaeBuilderInternal::var(size_t ind) const {
  return variable(ind).v;
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

size_t DaeBuilderInternal::find(const MX& v) const {
  // Make sure it is a valid input
  casadi_assert(v.is_symbolic(), "Variable must be symbolic");
  // Find the prospective variable
  size_t ind = find(v.name());
  // Make sure that the expression (not just the name) is correct
  casadi_assert(is_equal(v, variables_.at(ind)->v),
    "Variable \"" + v.name() + "\" has mismatching symbolic expression");
  // Return index
  return ind;
}

std::vector<size_t> DaeBuilderInternal::find(const std::vector<std::string>& name) const {
  std::vector<size_t> r(name.size());
  for (size_t i = 0; i < r.size(); ++i) r[i] = find(name[i]);
  return r;
}

std::vector<size_t> DaeBuilderInternal::find(const std::vector<MX>& v) const {
  std::vector<size_t> r(v.size());
  for (size_t i = 0; i < r.size(); ++i) r[i] = find(v[i]);
  return r;
}

const std::string& DaeBuilderInternal::name(size_t ind) const {
  return variable(ind).name;
}

std::vector<std::string> DaeBuilderInternal::name(const std::vector<size_t>& ind) const {
  std::vector<std::string> r(ind.size());
  for (size_t i = 0; i < r.size(); ++i) r[i] = name(ind[i]);
  return r;
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

void DaeBuilderInternal::reset() {
  for (Variable* v : variables_) {
    std::fill(v->value.begin(), v->value.end(), nan);
    v->stringvalue = std::string();
  }
}

double DaeBuilderInternal::attribute(Attribute a, const std::string& name) const {
  return variable(name).attribute(a);
}

std::vector<double> DaeBuilderInternal::attribute(Attribute a,
    const std::vector<std::string>& name) const {
  std::vector<double> r;
  r.reserve(name.size());
  for (auto& n : name) r.push_back(variable(n).attribute(a));
  return r;
}

void DaeBuilderInternal::set_attribute(Attribute a, const std::string& name, double val) {
  variable(name).set_attribute(a, val);
}

void DaeBuilderInternal::set_attribute(Attribute a, const std::vector<std::string>& name,
    const std::vector<double>& val) {
  casadi_assert(name.size() == val.size(), "Dimension mismatch");
  for (size_t k = 0; k < name.size(); ++k) variable(name[k]).set_attribute(a, val[k]);
}

std::string DaeBuilderInternal::string_attribute(Attribute a,
    const std::string& name) const {
  return variable(name).string_attribute(a);
}

std::vector<std::string> DaeBuilderInternal::string_attribute(Attribute a,
    const std::vector<std::string>& name) const {
  std::vector<std::string> r;
  r.reserve(name.size());
  for (auto& n : name) r.push_back(variable(n).string_attribute(a));
  return r;
}

void DaeBuilderInternal::set_string_attribute(Attribute a, const std::string& name,
    const std::string& val) {
  variable(name).set_string_attribute(a, val);
}

void DaeBuilderInternal::set_string_attribute(Attribute a,
    const std::vector<std::string>& name, const std::vector<std::string>& val) {
  casadi_assert(name.size() == val.size(), "Dimension mismatch");
  for (size_t k = 0; k < name.size(); ++k) variable(name[k]).set_string_attribute(a, val[k]);
}

Sparsity DaeBuilderInternal::jac_sparsity(const std::vector<size_t>& oind,
    const std::vector<size_t>& iind) const {
  // Mark inputs
  std::vector<casadi_int> lookup(n_variables(), -1);
  for (size_t i = 0; i < iind.size(); ++i)
    lookup.at(iind[i]) = i;
  // Sparsity pattern for the Jacobian block
  std::vector<casadi_int> row, col;
  // Loop over output nonzeros
  for (casadi_int j = 0; j < oind.size(); ++j) {
    for (casadi_int d : variable(oind[j]).dependencies) {
      casadi_int i = lookup.at(d);
      if (i >= 0) {
        row.push_back(j);  // Note: May not be sorted in ascending order
        col.push_back(i);
      }
    }
  }
  // Assemble sparsity in triplet format
  return Sparsity::triplet(oind.size(), iind.size(), row, col);
}

Sparsity DaeBuilderInternal::hess_sparsity(const std::vector<size_t>& oind,
    const std::vector<size_t>& iind) const {
  // Mark inputs
  std::vector<casadi_int> lookup(n_variables(), -1);
  for (size_t i = 0; i < iind.size(); ++i) lookup.at(iind[i]) = i;
  // Which variables enter as a nonlinear dependency in any variable in oind
  std::vector<bool> nonlin(iind.size(), false);
  // List of nonlinearly entering variables for the specific output
  std::vector<casadi_int> nonlin_list;
  // Rows and columns of the Hessian
  std::vector<casadi_int> row, col;
  // Loop over output variables
  for (casadi_int j = 0; j < oind.size(); ++j) {
    const Variable& v = variable(oind[j]);
    // Loop over dependencies
    for (size_t k = 0; k < v.dependencies.size(); ++k) {
      if (v.dependenciesKind.empty() || v.dependenciesKind.at(k) == DependenciesKind::DEPENDENT) {
        casadi_int i = lookup.at(v.dependencies[k]);
        if (i >= 0 && !nonlin.at(i)) {
          // Add to list
          nonlin_list.push_back(i);
          nonlin.at(i) = true;
        }
      }
    }
    // Add all combinations to sparsity pattern
    for (casadi_int k1 : nonlin_list) {
      for (casadi_int k2 : nonlin_list) {
        row.push_back(k1);
        col.push_back(k2);
      }
    }
    // If row/col vectors grow too large, remove duplicates
    if (col.size() > 2 * iind.size() * iind.size()) {
      Sparsity r = Sparsity::triplet(iind.size(), iind.size(), row, col);
      row = r.get_row();
      col = r.get_col();
    }
    // Reset nonlin, nonlin_list for next iteration
    for (casadi_int k : nonlin_list) nonlin[k] = false;
    nonlin_list.clear();
  }
  // Create sparsity pattern
  return Sparsity::triplet(iind.size(), iind.size(), row, col);
}

std::string DaeBuilderInternal::iso_8601_time() {
  // Get current time
  auto now = std::chrono::system_clock::now();
  std::time_t tt = std::chrono::system_clock::to_time_t(now);
  auto local_tm = *std::localtime(&tt);  // NOLINT(runtime/threadsafe_fn)
  // Convert to ISO 8601 (YYYY-MM-DDThh:mm:ssZ) format and return
  std::stringstream ss;
  ss << local_tm.tm_year + 1900 << '-';  // YYYY-
  ss << std::setfill('0') << std::setw(2) << local_tm.tm_mon + 1 << '-';  // MM-
  ss << std::setfill('0') << std::setw(2) << local_tm.tm_mday << 'T';  // DDT
  ss << std::setfill('0') << std::setw(2) << local_tm.tm_hour << ':';  // hh:
  ss << std::setfill('0') << std::setw(2) << local_tm.tm_min << ':';  // mm:
  ss << std::setfill('0') << std::setw(2) << local_tm.tm_sec << 'Z'; // ssZ
  return ss.str();
}

std::string DaeBuilderInternal::generate_guid() {
  // Initialize random seed
  static bool initialized = false;
  if (!initialized) {
      srand(::time(nullptr));  // NOLINT(runtime/threadsafe_fn)
    initialized = true;
  }
  // Possible characters
  const char h[] = "0123456789abcdef";
  // Length of GUID
  const size_t len = 32;
  // Generate random hex string
  std::vector<char> buf(len);
  for (size_t i = 0; i < len; ++i)
    buf[i] = h[rand() % 16];  // NOLINT(runtime/threadsafe_fn)
  return std::string(&buf.front(), len);
}

} // namespace casadi

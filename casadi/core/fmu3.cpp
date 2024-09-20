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


#include "fmu3.hpp"
#include "fmu_function.hpp"
#include "dae_builder_internal.hpp"

namespace casadi {

Fmu3::~Fmu3() {
}

std::string Fmu3::system_infix() {  
  // Architecture
  std::string arch;
#ifdef __arm__
  if (sizeof(void*) == 4) {
    // ARM 32-bit Architecture
    arch = "aarch32";
  } else {
    // ARM 64-bit Architecture
    arch = "aarch64";
  }
#else
  if (sizeof(void*) == 4) {
    // Intel/AMD x86 32-bit
    arch = "x86";
  } else {
    // Intel/AMD x86 64-bit
    arch = "x86_64";
  }
#endif  
  // Operating system
  std::string sys;
#if defined(_WIN32)
  // Microsoft Windows
  sys = "windows";
#elif defined(__APPLE__)
  // Darwin (macOS, iOS, watchOS, tvOS, audioOS)
  sys = "darwin";
#else
  // Linux
  sys = "linux";
#endif
  // Return platform tuple, according to Section 2.5.1.4.1. of the FMI 3 specification
  return arch + "-" + sys;
}

std::string Fmu3::dll_suffix() {
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

void Fmu3::init(const DaeBuilderInternal* dae) {
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

  // Collect input and parameter values
  vr_real_.clear();
  vr_integer_.clear();
  vr_boolean_.clear();
  vr_string_.clear();
  init_real_.clear();
  init_integer_.clear();
  init_boolean_.clear();
  init_string_.clear();
  for (size_t i = 0; i < dae->n_variables(); ++i) {
    const Variable& v = dae->variable(i);
    casadi_assert(v.numel == 1, "Vector variable support not implemented");
    // Skip if the wrong type
    if (v.causality != Causality::PARAMETER && v.causality != Causality::INPUT) continue;
    // Variable has not been set - keep default value
    if (!v.is_set()) continue;
    // Value reference
    fmi3ValueReference vr = v.value_reference;
    // Get value
    switch (to_fmi2(v.type)) {
      case TypeFmi2::REAL:
        init_real_.push_back(static_cast<fmi3Float64>(v.value.front()));
        vr_real_.push_back(vr);
        break;
      case TypeFmi2::INTEGER:
      case TypeFmi2::ENUM:
        init_integer_.push_back(static_cast<fmi3Int32>(v.value.front()));
        vr_integer_.push_back(vr);
        break;
      case TypeFmi2::BOOLEAN:
        init_boolean_.push_back(static_cast<fmi3Boolean>(v.value.front()));
        vr_boolean_.push_back(vr);
        break;
      case TypeFmi2::STRING:
        init_string_.push_back(v.stringvalue);
        vr_string_.push_back(vr);
        break;
      default:
        casadi_warning("Ignoring " + v.name + ", type: " + to_string(v.type));
    }
  }

  // Collect auxilliary variables
  vn_aux_real_.clear();
  vn_aux_integer_.clear();
  vn_aux_boolean_.clear();
  vn_aux_string_.clear();
  vr_aux_real_.clear();
  vr_aux_integer_.clear();
  vr_aux_boolean_.clear();
  vr_aux_string_.clear();
  for (auto&& s : aux_) {
    const Variable& v = dae->variable(s);
    // Convert to expected type
    fmi3ValueReference vr = v.value_reference;
    // Sort by type
    switch (to_fmi2(v.type)) {
      case TypeFmi2::REAL:
        // Real
        vn_aux_real_.push_back(v.name);
        vr_aux_real_.push_back(vr);
        break;
      case TypeFmi2::INTEGER:
      case TypeFmi2::ENUM:
        // Integer or enum
        vn_aux_integer_.push_back(v.name);
        vr_aux_integer_.push_back(vr);
        break;
      case TypeFmi2::BOOLEAN:
        // Boolean
        vn_aux_boolean_.push_back(v.name);
        vr_aux_boolean_.push_back(vr);
        break;
      case TypeFmi2::STRING:
        // String
        vn_aux_string_.push_back(v.name);
        vr_aux_string_.push_back(vr);
        break;
      default:
        casadi_warning("Ignoring " + v.name + ", type: " + to_string(v.type));
    }
  }

  /// Allocate numerical values for initial auxilliary variables
  aux_value_.v_real.resize(vn_aux_real_.size());
  aux_value_.v_integer.resize(vn_aux_integer_.size());
  aux_value_.v_boolean.resize(vn_aux_boolean_.size());
  aux_value_.v_string.resize(vn_aux_string_.size());

  // Get Jacobian sparsity information
  jac_sp_ = dae->jac_sparsity(oind_, iind_);

  // Get Hessian sparsity information
  hess_sp_ = dae->hess_sparsity(oind_, iind_);

  // Load DLL
  std::string instance_name_no_dot = dae->model_identifier_;
  std::replace(instance_name_no_dot.begin(), instance_name_no_dot.end(), '.', '_');
  std::string dll_path = dae->path_ + "/binaries/" + system_infix()
    + "/" + instance_name_no_dot + dll_suffix();
  li_ = Importer(dll_path, "dll");

  // Types of analytic AD, if any
  provides_directional_derivatives_ = dae->provides_directional_derivatives_;
  provides_adjoint_derivatives_ = dae->provides_adjoint_derivatives_;

  // Path to resource directory
  resource_loc_ = "file://" + dae->path_ + "/resources";

  // Copy info from DaeBuilder
  fmutol_ = dae->fmutol_;
  instance_name_ = dae->model_identifier_;
  guid_ = dae->instantiation_token_;
  logging_on_ = dae->debug_;
}

void Fmu3::finalize() {
  // Get FMI C functions
  instantiate_model_exchange_ = load_function<fmi3InstantiateModelExchangeTYPE>(
    "fmi3InstantiateModelExchange");
  free_instance_ = load_function<fmi3FreeInstanceTYPE>("fmi3FreeInstance");
  reset_ = load_function<fmi3ResetTYPE>("fmi3Reset");
  enter_initialization_mode_ = load_function<fmi3EnterInitializationModeTYPE>(
    "fmi3EnterInitializationMode");
  exit_initialization_mode_ = load_function<fmi3ExitInitializationModeTYPE>(
    "fmi3ExitInitializationMode");
  enter_continuous_time_mode_ = load_function<fmi3EnterContinuousTimeModeTYPE>(
    "fmi3EnterContinuousTimeMode");
  get_float64_ = load_function<fmi3GetFloat64TYPE>("fmi3GetFloat64");
  set_float64_ = load_function<fmi3SetFloat64TYPE>("fmi3SetFloat64");
  get_int32_ = load_function<fmi3GetInt32TYPE>("fmi3GetInt32");
  set_int32_ = load_function<fmi3SetInt32TYPE>("fmi3SetInt32");
  get_boolean_ = load_function<fmi3GetBooleanTYPE>("fmi3GetBoolean");
  set_boolean_ = load_function<fmi3SetBooleanTYPE>("fmi3SetBoolean");
  get_string_ = load_function<fmi3GetStringTYPE>("fmi3GetString");
  set_string_ = load_function<fmi3SetStringTYPE>("fmi3SetString");
  if (provides_directional_derivatives_) {
    get_directional_derivative_ =
      load_function<fmi3GetDirectionalDerivativeTYPE>("fmi3GetDirectionalDerivative");
  }
  if (provides_adjoint_derivatives_) {
    get_adjoint_derivative_ =
      load_function<fmi3GetAdjointDerivativeTYPE>("fmi3GetAdjointDerivative");
  }

  // Create a temporary instance
  fmi3Instance c = instantiate();
  // Set all values
  if (set_values(c)) {
    casadi_error("Fmu3::set_values failed");
  }
  // Initialization mode begins
  if (enter_initialization_mode(c)) {
    casadi_error("Fmu3::enter_initialization_mode failed");
  }
  // Get input values
  if (get_in(c, &value_in_)) {
    casadi_error("Fmu3::get_in failed");
  }
  // Get auxilliary variables
  if (get_aux(c, &aux_value_)) {
    casadi_error("Fmu3::get_aux failed");
  }
  // Free memory
  free_instance(c);
}

void Fmu3::log_message_callback(fmi3InstanceEnvironment instanceEnvironment,
    fmi3Status status, fmi3String category, fmi3String message) {
  // Print message content
  uout() << "[" << category << "] " << message << std::endl;
}

fmi3Instance Fmu3::instantiate() const {
  // Instantiate FMU
  fmi3String instanceName = instance_name_.c_str();
  fmi3String instantiationToken = guid_.c_str();
  fmi3String resourcePath = resource_loc_.c_str();
  fmi3Boolean visible = fmi3False;
  fmi3InstanceEnvironment instanceEnvironment = 0;
  fmi3Instance c = instantiate_model_exchange_(instanceName, instantiationToken, resourcePath,
    visible, logging_on_, instanceEnvironment, log_message_callback);
  if (c == 0) casadi_error("fmi3Instantiate failed");
  return c;
}

void Fmu3::free_instance(void* c) const {
  if (free_instance_) {
    free_instance_(static_cast<fmi3Instance>(c));
  } else {
    casadi_warning("No free_instance function pointer available");
  }
}

int Fmu3::init_mem(FmuMemory* m) const {
  // Ensure not already instantiated
  casadi_assert(m->instance == 0, "Already instantiated");
  // Create instance
  m->instance = instantiate();
  // Set all values
  if (set_values(m->instance)) {
    casadi_warning("Fmu3::set_values failed");
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

int Fmu3::reset(fmi3Instance c) {
  fmi3Status status = reset_(c);
  if (status != fmi3OK) {
    casadi_warning("fmi3Reset failed");
    return 1;
  }
  return 0;
}

int Fmu3::enter_initialization_mode(fmi3Instance c) const {
  fmi3Status status = enter_initialization_mode_(c, fmutol_ > 0, fmutol_, 0., fmi3True, 1.);
  if (status != fmi3OK) {
    casadi_warning("fmi3EnterInitializationMode failed: " + str(status));
    return 1;
  }
  return 0;
}

int Fmu3::exit_initialization_mode(fmi3Instance c) const {
  fmi3Status status = exit_initialization_mode_(c);
  if (status != fmi3OK) {
    casadi_warning("fmi3ExitInitializationMode failed");
    return 1;
  }
  return 0;
}

int Fmu3::set_values(fmi3Instance c) const {
  // Pass real values before initialization
  if (!vr_real_.empty()) {
    fmi3Status status = set_float64_(c, get_ptr(vr_real_), vr_real_.size(),
      get_ptr(init_real_), init_real_.size());
    if (status != fmi3OK) {
      casadi_warning("fmi3SetFloat64 failed");
      return 1;
    }
  }
  // Pass integer values before initialization (also enums)
  if (!vr_integer_.empty()) {
    fmi3Status status = set_int32_(c, get_ptr(vr_integer_), vr_integer_.size(),
      get_ptr(init_integer_), init_integer_.size());
    if (status != fmi3OK) {
      casadi_warning("fmi3SetInt32 failed");
      return 1;
    }
  }
  // Pass boolean values before initialization
  if (!vr_boolean_.empty()) {
    casadi_error("Broken");
    // nullptr -> get_ptr(init_boolean_)
    fmi3Status status = set_boolean_(c, get_ptr(vr_boolean_), vr_boolean_.size(),
      nullptr, init_boolean_.size());
    if (status != fmi3OK) {
      casadi_warning("fmi3SetBoolean failed");
      return 1;
    }
  }
  // Pass string values before initialization
  for (size_t k = 0; k < vr_string_.size(); ++k) {
    fmi3ValueReference vr = vr_string_[k];
    fmi3String value = init_string_[k].c_str();
    fmi3Status status = set_string_(c, &vr, 1, &value, 1);
    if (status != fmi3OK) {
      casadi_error("fmi3SetString failed for value reference " + str(vr));
    }
  }
  // Successful return
  return 0;
}

int Fmu3::get_in(fmi3Instance c, std::vector<fmi3Float64>* v) const {
  if (!vr_in_.empty()) {
    fmi3Status status = get_float64_(c, get_ptr(vr_in_), vr_in_.size(),
      get_ptr(*v), vr_in_.size());
    if (status != fmi3OK) {
      casadi_warning("fmi3GetFloat64 failed");
      return 1;
    }
  }
  // Successful return
  return 0;
}

int Fmu3::get_aux(fmi3Instance c, Value* v) const {
  // Get real auxilliary variables
  if (!vr_aux_real_.empty()) {
    fmi3Status status = get_float64_(c, get_ptr(vr_aux_real_), vr_aux_real_.size(),
      get_ptr(v->v_real), vr_aux_real_.size());
    if (status != fmi3OK) {
      casadi_warning("fmi3GetFloat64 failed");
      return 1;
    }
  }
  // Get integer/enum auxilliary variables
  if (!vr_aux_integer_.empty()) {
    fmi3Status status = get_int32_(c, get_ptr(vr_aux_integer_), vr_aux_integer_.size(),
      get_ptr(v->v_integer), vr_aux_integer_.size());
    if (status != fmi3OK) {
      casadi_warning("fmi3GetInt32 failed");
      return 1;
    }
  }
  // Get boolean auxilliary variables
  if (!vr_aux_boolean_.empty()) {
    casadi_error("Broken");
    // nullptr -> get_ptr(v->v_boolean)
    fmi3Status status = get_boolean_(c, get_ptr(vr_aux_boolean_), vr_aux_boolean_.size(),
      nullptr, vr_aux_boolean_.size());
    if (status != fmi3OK) {
      casadi_warning("fmi3GetBoolean failed");
      return 1;
    }
  }
  // Get string auxilliary variables
  for (size_t k = 0; k < vr_aux_string_.size(); ++k) {
    fmi3ValueReference vr = vr_aux_string_[k];
    fmi3String value = v->v_string.at(k).c_str();
    fmi3Status status = set_string_(c, &vr, 1, &value, 1);
    if (status != fmi3OK) {
      casadi_error("fmi3GetString failed for value reference " + str(vr));
    }
  }
  // Successful return
  return 0;
}

void Fmu3::get_stats(FmuMemory* m, Dict* stats,
    const std::vector<std::string>& name_in, const InputStruct* in) const {
  // To do: Use auxillary variables from last evaluation
  (void)m;  // unused
  // Auxilliary values to be copied
  const Value& v = aux_value_;
  // Collect auxilliary variables
  Dict aux;
  // Real
  for (size_t k = 0; k < vn_aux_real_.size(); ++k) {
    aux[vn_aux_real_[k]] = static_cast<double>(v.v_real[k]);
  }
  // Integer
  for (size_t k = 0; k < vn_aux_integer_.size(); ++k) {
    aux[vn_aux_integer_[k]] = static_cast<casadi_int>(v.v_integer[k]);
  }
  // Boolean
  for (size_t k = 0; k < vn_aux_boolean_.size(); ++k) {
    aux[vn_aux_boolean_[k]] = static_cast<bool>(v.v_boolean[k]);
  }
  // String
  for (size_t k = 0; k < vn_aux_string_.size(); ++k) {
    aux[vn_aux_string_[k]] = v.v_string[k];
  }
  // Copy to stats
  (*stats)["aux"] = aux;
  // Loop over input variables
  for (size_t k = 0; k < name_in.size(); ++k) {
    // Only consider regular inputs
    if (in[k].type == InputType::REG) {
      // Get the indices
      const std::vector<size_t>& iind = ired_.at(in[k].ind);
      // Collect values
      std::vector<double> v(iind.size());
      for (size_t i = 0; i < v.size(); ++i) v[i] = value_in_.at(iind[i]);
      // Save to stats
      (*stats)[name_in[k]] = v;
    }
  }
}

int Fmu3::eval(FmuMemory* m) const {
  // Gather inputs and outputs
  gather_io(m);
  // Number of inputs and outputs
  size_t n_set = m->id_in_.size();
  size_t n_out = m->id_out_.size();
  // Fmi return flag
  fmi3Status status;
  // Set all variables
  status = set_float64_(m->instance, get_ptr(m->vr_in_), n_set, get_ptr(m->v_in_), n_set);
  if (status != fmi3OK) {
    casadi_warning("fmi3SetFloat64 failed");
    return 1;
  }
  // Quick return if nothing requested
  if (n_out == 0) return 0;
  // Calculate all variables
  m->v_out_.resize(n_out);
  status = get_float64_(m->instance, get_ptr(m->vr_out_), n_out, get_ptr(m->v_out_), n_out);
  if (status != fmi3OK) {
    casadi_warning("fmi3GetFloat64 failed");
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

int Fmu3::eval_ad(FmuMemory* m) const {
  // Number of inputs and outputs
  size_t n_known = m->id_in_.size();
  size_t n_unknown = m->id_out_.size();
  // Quick return if nothing to be calculated
  if (n_unknown == 0) return 0;
  // Evalute (should not be necessary)
  fmi3Status status = get_float64_(m->instance, get_ptr(m->vr_out_), n_unknown,
    get_ptr(m->v_out_), n_unknown);
  if (status != fmi3OK) {
    casadi_warning("fmi3GetFloat64 failed");
    return 1;
  }
  // Evaluate directional derivatives
  status = get_directional_derivative_(m->instance, get_ptr(m->vr_out_), n_unknown,
    get_ptr(m->vr_in_), n_known, get_ptr(m->d_in_), n_known, get_ptr(m->d_out_), n_unknown);
  if (status != fmi3OK) {
    casadi_warning("fmi3GetDirectionalDerivative failed");
    return 1;
  }
  // Collect requested variables
  auto it = m->d_out_.begin();
  for (size_t id : m->id_out_) {
    m->sens_[id] = *it++;
  }
  // Successful return
  return 0;
}

int Fmu3::eval_fd(FmuMemory* m, bool independent_seeds) const {
  // Number of inputs and outputs
  size_t n_known = m->id_in_.size();
  size_t n_unknown = m->id_out_.size();
  // Quick return if nothing to be calculated
  if (n_unknown == 0) return 0;
  // Evalute (should not be necessary)
  fmi3Status status = get_float64_(m->instance, get_ptr(m->vr_out_), n_unknown,
    get_ptr(m->v_out_), n_unknown);
  if (status != fmi3OK) {
    casadi_warning("fmi3GetFloat64 failed");
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
  const fmi3Float64* yk_all[5] = {0};

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
    status = set_float64_(m->instance, get_ptr(m->vr_in_), n_known, get_ptr(m->v_pert_), n_known);
    if (status != fmi3OK) {
      casadi_warning("fmi3SetFloat64 failed");
      return 1;
    }
    // Evaluate perturbed FMU
    status = get_float64_(m->instance, get_ptr(m->vr_out_), n_unknown, yk, n_unknown);
    if (status != fmi3OK) {
      casadi_warning("fmi3GetFloat64 failed");
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
  status = set_float64_(m->instance, get_ptr(m->vr_in_), n_known, get_ptr(m->v_in_), n_known);
  if (status != fmi3OK) {
    casadi_warning("fmi3SetFloat64 failed");
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
    if (m->self.validate_ad_) {
      // Value to compare with
      double d_ad = m->sens_[id];
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
      m->sens_[id] = d_fd;
    }
  }
  // Successful return
  return 0;
}

Fmu3::Fmu3(const std::string& name,
    const std::vector<std::string>& scheme_in,
    const std::vector<std::string>& scheme_out,
    const std::map<std::string, std::vector<size_t>>& scheme,
    const std::vector<std::string>& aux)
    : FmuInternal(name, scheme_in, scheme_out, scheme, aux) {
  instantiate_model_exchange_ = 0;
  free_instance_ = 0;
  reset_ = 0;
  enter_initialization_mode_ = 0;
  exit_initialization_mode_ = 0;
  enter_continuous_time_mode_ = 0;
  set_float64_ = 0;
  set_boolean_ = 0;
  get_float64_ = 0;
  get_directional_derivative_ = 0;
  get_adjoint_derivative_ = 0;
}

} // namespace casadi

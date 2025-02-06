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

std::string Fmu3::system_infix() const {
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

void Fmu3::init(const DaeBuilderInternal* dae) {
  // Initialize base classes
  FmuInternal::init(dae);

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
}

void Fmu3::load_functions() {
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
  set_time_ = load_function<fmi3SetTimeTYPE>("fmi3SetTime");
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
  if (true) {
    update_discrete_states_ =
      load_function<fmi3UpdateDiscreteStatesTYPE>("fmi3UpdateDiscreteStates");
  }
}

void Fmu3::log_message_callback(fmi3InstanceEnvironment instanceEnvironment,
    fmi3Status status, fmi3String category, fmi3String message) {
  // Print message content
  uout() << "[" << category << "] " << message << std::endl;
}

void* Fmu3::instantiate() const {
  // Instantiate FMU
  fmi3String instanceName = instance_name_.c_str();
  fmi3String instantiationToken = instantiation_token_.c_str();
  fmi3String resourcePath = resource_loc_.c_str();
  fmi3Boolean visible = fmi3False;
  fmi3InstanceEnvironment instanceEnvironment = 0;
  fmi3Instance c = instantiate_model_exchange_(instanceName, instantiationToken, resourcePath,
    visible, logging_on_, instanceEnvironment, log_message_callback);
  if (c == 0) casadi_error("fmi3Instantiate failed");
  return c;
}

void Fmu3::free_instance(void* instance) const {
  if (free_instance_) {
    auto c = static_cast<fmi3Instance>(instance);
    free_instance_(c);
  } else {
    casadi_warning("No free_instance function pointer available");
  }
}

int Fmu3::reset(void* instance) {
  auto c = static_cast<fmi3Instance>(instance);
  fmi3Status status = reset_(c);
  if (status != fmi3OK) {
    casadi_warning("fmi3Reset failed");
    return 1;
  }
  return 0;
}

int Fmu3::enter_initialization_mode(void* instance) const {
  auto c = static_cast<fmi3Instance>(instance);
  fmi3Status status = enter_initialization_mode_(c, fmutol_ > 0, fmutol_, 0., fmi3True, 1.);
  if (status != fmi3OK) {
    casadi_warning("fmi3EnterInitializationMode failed: " + str(status));
    return 1;
  }
  return 0;
}

int Fmu3::exit_initialization_mode(void* instance) const {
  auto c = static_cast<fmi3Instance>(instance);
  fmi3Status status = exit_initialization_mode_(c);
  if (status != fmi3OK) {
    casadi_warning("fmi3ExitInitializationMode failed");
    return 1;
  }
  return 0;
}

int Fmu3::update_discrete_states(void* instance, EventMemory* eventmem) const {
  auto c = static_cast<fmi3Instance>(instance);
  // Return arguments in FMI types
  fmi3Boolean discreteStatesNeedUpdate, terminateSimulation, nominalsOfContinuousStatesChanged,
    valuesOfContinuousStatesChanged, nextEventTimeDefined;
  fmi3Float64 nextEventTime;
  // Call FMU
  fmi3Status status = update_discrete_states_(c, &discreteStatesNeedUpdate, &terminateSimulation,
    &nominalsOfContinuousStatesChanged, &valuesOfContinuousStatesChanged, &nextEventTimeDefined,
    &nextEventTime);
  // Pass to event iteration memory
  eventmem->discrete_states_need_update = discreteStatesNeedUpdate;
  eventmem->terminate_simulation = terminateSimulation;
  eventmem->nominals_of_continuous_states_changed = nominalsOfContinuousStatesChanged;
  eventmem->values_of_continuous_states_changed = nominalsOfContinuousStatesChanged;
  eventmem->next_event_time_defined = nextEventTimeDefined;
  eventmem->next_event_time = nextEventTime;
  return status != fmi3OK;
}

int Fmu3::set_real(void* instance, const unsigned int* vr, size_t n_vr,
    const double* values, size_t n_values) const {
  // Set time variable, if any
  if (has_independent_ && n_vr > 0 && *vr == vr_in_[0]) {
    // Update FMU time
    fmi3Status status = set_time_(instance, *values);
    if (status != fmi3OK) return 1;
    // Skip when setting remaining variables
    vr++;
    n_vr--;
    values++;
    n_values--;
  }
  // Set remaining variables
  fmi3Status status = set_float64_(instance, vr, n_vr, values, n_values);
  return status != fmi3OK;
}

int Fmu3::get_real(void* instance, const unsigned int* vr, size_t n_vr,
    double* values, size_t n_values) const {
  fmi3Status status = get_float64_(instance, vr, n_vr, values, n_values);
  return status != fmi3OK;
}

int Fmu3::get_directional_derivative(void* instance, const unsigned int* vr_out, size_t n_out,
    const unsigned int* vr_in, size_t n_in, const double* seed, size_t n_seed,
    double* sensitivity, size_t n_sensitivity) const {
  fmi3Status status = get_directional_derivative_(instance, vr_out, n_out, vr_in, n_in,
    seed, n_seed, sensitivity, n_sensitivity);
  return status != fmi3OK;
}

int Fmu3::get_adjoint_derivative(void* instance, const unsigned int* vr_out, size_t n_out,
    const unsigned int* vr_in, size_t n_in, const double* seed, size_t n_seed,
    double* sensitivity, size_t n_sensitivity) const {
  fmi3Status status = get_adjoint_derivative_(instance, vr_out, n_out, vr_in, n_in,
    seed, n_seed, sensitivity, n_sensitivity);
  return status != fmi3OK;
}

int Fmu3::set_values(void* instance) const {
  auto c = static_cast<fmi3Instance>(instance);
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

int Fmu3::get_aux(void* instance) {
  auto c = static_cast<fmi3Instance>(instance);
  // Get real auxilliary variables
  if (!vr_aux_real_.empty()) {
    fmi3Status status = get_float64_(c, get_ptr(vr_aux_real_), vr_aux_real_.size(),
      get_ptr(aux_value_.v_real), vr_aux_real_.size());
    if (status != fmi3OK) {
      casadi_warning("fmi3GetFloat64 failed");
      return 1;
    }
  }
  // Get integer/enum auxilliary variables
  if (!vr_aux_integer_.empty()) {
    fmi3Status status = get_int32_(c, get_ptr(vr_aux_integer_), vr_aux_integer_.size(),
      get_ptr(aux_value_.v_integer), vr_aux_integer_.size());
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
    fmi3String value = aux_value_.v_string.at(k).c_str();
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
  set_time_ = 0;
  set_float64_ = 0;
  set_boolean_ = 0;
  get_float64_ = 0;
  get_directional_derivative_ = 0;
  get_adjoint_derivative_ = 0;
  update_discrete_states_ = 0;
}

} // namespace casadi

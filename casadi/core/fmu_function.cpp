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


#include "fmu_function.hpp"
#include "casadi_misc.hpp"
#include "serializing_stream.hpp"
#include "dae_builder_internal.hpp"

#include <fstream>
#include <iostream>
#include <sstream>

namespace casadi {

#ifdef WITH_FMU

std::string Fmu::system_infix() {
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

std::string Fmu::dll_suffix() {
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

void Fmu::init() {
  // DaeBuilder instance
  auto dae = this->dae();

  // Load DLL
  std::string instance_name_no_dot = dae->model_identifier_;
  std::replace(instance_name_no_dot.begin(), instance_name_no_dot.end(), '.', '_');
  std::string dll_path = dae->path_ + "/binaries/" + system_infix()
    + "/" + instance_name_no_dot + dll_suffix();
  li_ = Importer(dll_path, "dll");

  // Get FMI C functions
  instantiate_ = reinterpret_cast<fmi2InstantiateTYPE*>(load_function("fmi2Instantiate"));
  free_instance_ = reinterpret_cast<fmi2FreeInstanceTYPE*>(load_function("fmi2FreeInstance"));
  reset_ = reinterpret_cast<fmi2ResetTYPE*>(load_function("fmi2Reset"));
  setup_experiment_ = reinterpret_cast<fmi2SetupExperimentTYPE*>(
    load_function("fmi2SetupExperiment"));
  enter_initialization_mode_ = reinterpret_cast<fmi2EnterInitializationModeTYPE*>(
    load_function("fmi2EnterInitializationMode"));
  exit_initialization_mode_ = reinterpret_cast<fmi2ExitInitializationModeTYPE*>(
    load_function("fmi2ExitInitializationMode"));
  enter_continuous_time_mode_ = reinterpret_cast<fmi2EnterContinuousTimeModeTYPE*>(
    load_function("fmi2EnterContinuousTimeMode"));
  get_real_ = reinterpret_cast<fmi2GetRealTYPE*>(load_function("fmi2GetReal"));
  set_real_ = reinterpret_cast<fmi2SetRealTYPE*>(load_function("fmi2SetReal"));
  get_integer_ = reinterpret_cast<fmi2GetIntegerTYPE*>(load_function("fmi2GetInteger"));
  set_integer_ = reinterpret_cast<fmi2SetIntegerTYPE*>(load_function("fmi2SetInteger"));
  get_boolean_ = reinterpret_cast<fmi2GetBooleanTYPE*>(load_function("fmi2GetBoolean"));
  set_boolean_ = reinterpret_cast<fmi2SetBooleanTYPE*>(load_function("fmi2SetBoolean"));
  get_string_ = reinterpret_cast<fmi2GetStringTYPE*>(load_function("fmi2GetString"));
  set_string_ = reinterpret_cast<fmi2SetStringTYPE*>(load_function("fmi2SetString"));
  if (dae->provides_directional_derivative_) {
    get_directional_derivative_ = reinterpret_cast<fmi2GetDirectionalDerivativeTYPE*>(
      load_function("fmi2GetDirectionalDerivative"));
  }

  // Callback functions
  functions_.logger = logger;
  functions_.allocateMemory = calloc;
  functions_.freeMemory = free;
  functions_.stepFinished = 0;
  functions_.componentEnvironment = 0;

  // Path to resource directory
  resource_loc_ = "file://" + dae->path_ + "/resources";


  // Collect variables
  vr_real_.clear();
  vr_integer_.clear();
  vr_boolean_.clear();
  vr_string_.clear();
  init_real_.clear();
  init_integer_.clear();
  init_boolean_.clear();
  init_string_.clear();

  // Collect input and parameter values
  for (const Variable& v : dae->variables_) {
    // Skip if the wrong type
    if (v.causality != Causality::PARAMETER && v.causality != Causality::INPUT) continue;
    // If nan - variable has not been set - keep default value
    if (std::isnan(v.value)) continue;
    // Value reference
    fmi2ValueReference vr = v.value_reference;
    // Get value
    switch (v.type) {
      case Type::REAL:
        init_real_.push_back(static_cast<fmi2Real>(v.value));
        vr_real_.push_back(vr);
        break;
      case Type::INTEGER:
      case Type::ENUM:
        init_integer_.push_back(static_cast<fmi2Integer>(v.value));
        vr_integer_.push_back(vr);
        break;
      case Type::BOOLEAN:
        init_boolean_.push_back(static_cast<fmi2Boolean>(v.value));
        vr_boolean_.push_back(vr);
        break;
      case Type::STRING:
        init_string_.push_back(v.stringvalue);
        vr_string_.push_back(vr);
        break;
      default:
        casadi_warning("Ignoring " + v.name + ", type: " + to_string(v.type));
    }
  }
}

signal_t Fmu::load_function(const std::string& symname) {
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

fmi2Component FmuFunction::instantiate() const {
  // DaeBuilder instance
  auto dae = fmu_->dae();
  // Instantiate FMU
  fmi2String instanceName = dae->model_identifier_.c_str();
  fmi2Type fmuType = fmi2ModelExchange;
  fmi2String fmuGUID = dae->guid_.c_str();
  fmi2String fmuResourceLocation = fmu_->resource_loc_.c_str();
  fmi2Boolean visible = fmi2False;
  fmi2Boolean loggingOn = dae->debug_;
  fmi2Component c = fmu_->instantiate_(instanceName, fmuType, fmuGUID, fmuResourceLocation,
    &fmu_->functions_, visible, loggingOn);
  if (c == 0) casadi_error("fmi2Instantiate failed");
  return c;
}

int FmuFunction::init_mem(void* mem) const {
  FmuMemory* m = static_cast<FmuMemory*>(mem);
  casadi_assert(m != 0, "Memory is null");
  // Instantiate base classes
  if (FunctionInternal::init_mem(mem)) return 1;
  // DaeBuilder instance
  auto dae = fmu_->dae();
  // Ensure not already instantiated
  casadi_assert(m->c == 0, "Already instantiated");
  // Create instance
  m->c = instantiate();
  // Reset solver
  setup_experiment(m);
  // Set all values
  if (fmu_->set_values(m)) {
    casadi_warning("Fmu::set_values failed");
    return 1;
  }
  // Allocate/reset value buffer
  m->buffer_.resize(dae->variables_.size());
  std::fill(m->buffer_.begin(), m->buffer_.end(), casadi::nan);
  // Allocate/reset sensitivities
  m->sens_.resize(dae->variables_.size());
  std::fill(m->sens_.begin(), m->sens_.end(), 0);
  // Allocate/reset changed
  m->changed_.resize(dae->variables_.size());
  std::fill(m->changed_.begin(), m->changed_.end(), false);
  // Allocate/reset requested
  m->requested_.resize(dae->variables_.size());
  std::fill(m->requested_.begin(), m->requested_.end(), false);
  // Also allocate memory for corresponding Jacobian entry (for debugging)
  m->wrt_.resize(dae->variables_.size());
  // Initialization mode begins
  if (fmu_->enter_initialization_mode(m)) return 1;
  // Get all values
  if (fmu_->get_values(m)) {
    casadi_warning("Fmu::get_values failed");
    return 1;
  }
  // Initialization mode ends
  if (fmu_->exit_initialization_mode(m)) return 1;
  // Successful return
  return 0;
}

void FmuFunction::free_mem(void *mem) const {
  if (mem == 0) {
    casadi_warning("Memory is null");
    return;
  }
  // Free FMI memory
  FmuMemory* m = static_cast<FmuMemory*>(mem);
  if (m->c) {
    if (fmu_->free_instance_) {
      fmu_->free_instance_(m->c);
      m->c = nullptr;
    } else {
      casadi_warning("No free_instance function pointer available");
    }
  }
  delete m;
}

void FmuFunction::setup_experiment(FmuMemory* m) const {
  // DaeBuilder instance
  auto dae = fmu_->dae();
  // Call fmi2SetupExperiment
  fmi2Status status = fmu_->setup_experiment_(m->c, dae->fmutol_ > 0, dae->fmutol_, 0.,
    fmi2True, 1.);
  casadi_assert(status == fmi2OK, "fmi2SetupExperiment failed");
}

int Fmu::reset(FmuMemory* m) {
  fmi2Status status = reset_(m->c);
  if (status != fmi2OK) {
    casadi_warning("fmi2Reset failed");
    return 1;
  }
  return 0;
}

int Fmu::enter_initialization_mode(FmuMemory* m) const {
  fmi2Status status = enter_initialization_mode_(m->c);
  if (status != fmi2OK) {
    casadi_warning("fmi2EnterInitializationMode failed: " + str(status));
    return 1;
  }
  return 0;
}

int Fmu::exit_initialization_mode(FmuMemory* m) const {
  fmi2Status status = exit_initialization_mode_(m->c);
  if (status != fmi2OK) {
    casadi_warning("fmi2ExitInitializationMode failed");
    return 1;
  }
  return 0;
}

int Fmu::set_values(FmuMemory* m) const {
  // Pass real values before initialization
  if (!vr_real_.empty()) {
    fmi2Status status = set_real_(m->c, get_ptr(vr_real_), vr_real_.size(), get_ptr(init_real_));
    if (status != fmi2OK) {
      casadi_warning("fmi2SetReal failed");
      return 1;
    }
  }
  // Pass integer values before initialization (also enums)
  if (!vr_integer_.empty()) {
    fmi2Status status = set_integer_(m->c, get_ptr(vr_integer_), vr_integer_.size(),
      get_ptr(init_integer_));
    if (status != fmi2OK) {
      casadi_warning("fmi2SetInteger failed");
      return 1;
    }
  }
  // Pass boolean values before initialization
  if (!vr_boolean_.empty()) {
    fmi2Status status = set_boolean_(m->c, get_ptr(vr_boolean_), vr_boolean_.size(),
      get_ptr(init_boolean_));
    if (status != fmi2OK) {
      casadi_warning("fmi2SetBoolean failed");
      return 1;
    }
  }
  // Pass string valeus before initialization
  for (size_t k = 0; k < vr_string_.size(); ++k) {
    fmi2ValueReference vr = vr_string_[k];
    fmi2String value = init_string_[k].c_str();
    fmi2Status status = set_string_(m->c, &vr, 1, &value);
    if (status != fmi2OK) {
      casadi_error("fmi2SetString failed for value reference " + str(vr));
    }
  }
  // Successful return
  return 0;
}

int Fmu::get_values(FmuMemory* m) const {
  // DaeBuilder instance
  auto dae = this->dae();
  // Retrieve values
  for (Variable& v : dae->variables_) {
    // Convert to expected type
    fmi2ValueReference vr = v.value_reference;
    // Get value
    switch (v.type) {
      // Skip if the wrong type
      if (v.causality == Causality::PARAMETER || v.causality == Causality::INPUT) continue;
      // Get by type
      case Type::REAL:
        {
          // Real
          fmi2Real value;
          fmi2Status status = get_real_(m->c, &vr, 1, &value);
          if (status != fmi2OK) {
            casadi_warning("fmi2GetReal failed for " + v.name);
            return 1;
          }
          v.value = value;
          break;
        }
      case Type::INTEGER:
      case Type::ENUM:
        {
          // Integer: Convert to double
          fmi2Integer value;
          fmi2Status status = get_integer_(m->c, &vr, 1, &value);
          if (status != fmi2OK) {
            casadi_warning("fmi2GetInteger failed for " + v.name);
            return 1;
          }
          v.value = value;
          break;
        }
      case Type::BOOLEAN:
        {
          // Boolean: Convert to double
          fmi2Boolean value;
          fmi2Status status = get_boolean_(m->c, &vr, 1, &value);
          if (status != fmi2OK) {
            casadi_warning("fmi2GetBoolean failed for " + v.name);
            return 1;
          }
          v.value = value;
          break;
        }
      case Type::STRING:
        {
          // String
          fmi2String value;
          fmi2Status status = get_string_(m->c, &vr, 1, &value);
          if (status != fmi2OK) {
            casadi_warning("fmi2GetString failed for " + v.name);
            return 1;
          }
          v.stringvalue = value;
          v.value = 0;
          break;
        }
      default:
        casadi_warning("Ignoring " + v.name + ", type: " + to_string(v.type));
    }
  }
  // Successful return
  return 0;
}

void FmuFunction::set(FmuMemory* m, size_t id, double value) const {
  // Update buffer
  if (value != m->buffer_.at(id)) {
    m->buffer_.at(id) = value;
    m->changed_.at(id) = true;
  }
}

void FmuFunction::set_seed(FmuMemory* m, size_t id, double value) const {
  // Update buffer
  if (value != 0) {
    m->sens_.at(id) = value;
    m->changed_.at(id) = true;
  }
}

void FmuFunction::request(FmuMemory* m, size_t id, size_t wrt_id) const {
  // Mark as requested
  m->requested_.at(id) = true;
  // Also log corresponding input index
  m->wrt_.at(id) = wrt_id;
}

int FmuFunction::eval(FmuMemory* m) const {
  // Gather inputs and outputs
  gather_io(m);
  // Number of inputs and outputs
  size_t n_set = m->id_in_.size();
  size_t n_out = m->id_out_.size();
  // Fmi return flag
  fmi2Status status;
  // Set all variables
  status = fmu_->set_real_(m->c, get_ptr(m->vr_in_), n_set, get_ptr(m->v_in_));
  if (status != fmi2OK) {
    casadi_warning("fmi2SetReal failed");
    return 1;
  }
  // Quick return if nothing requested
  if (n_out == 0) return 0;
  // Calculate all variables
  m->v_out_.resize(n_out);
  status = fmu_->get_real_(m->c, get_ptr(m->vr_out_), n_out, get_ptr(m->v_out_));
  if (status != fmi2OK) {
    casadi_warning("fmi2GetReal failed");
    return 1;
  }
  // Collect requested variables
  auto it = m->v_out_.begin();
  for (size_t id : m->id_out_) {
    m->buffer_[id] = *it++;
  }
  // Successful return
  return 0;
}

void FmuFunction::get(FmuMemory* m, size_t id, double* value) const {
  // Save to return
  *value = m->buffer_.at(id);
}

void FmuFunction::gather_io(FmuMemory* m) const {
  // DaeBuilder instance
  auto dae = fmu_->dae();
  // Collect input indices and corresponding value references and values
  m->id_in_.clear();
  m->vr_in_.clear();
  m->v_in_.clear();
  for (size_t id = 0; id < m->changed_.size(); ++id) {
    if (m->changed_[id]) {
      const Variable& v = dae->variable(id);
      m->id_in_.push_back(id);
      m->vr_in_.push_back(v.value_reference);
      m->v_in_.push_back(m->buffer_[id]);
      m->changed_[id] = false;
    }
  }
  // Collect output indices, corresponding value references
  m->id_out_.clear();
  m->vr_out_.clear();
  for (size_t id = 0; id < m->requested_.size(); ++id) {
    if (m->requested_[id]) {
      const Variable& v = dae->variable(id);
      m->id_out_.push_back(id);
      m->vr_out_.push_back(v.value_reference);
      m->requested_[id] = false;
    }
  }
}

void FmuFunction::gather_sens(FmuMemory* m) const {
  // Gather input and output indices
  gather_io(m);
  // Number of inputs and outputs
  size_t n_known = m->id_in_.size();
  size_t n_unknown = m->id_out_.size();
  // Get/clear seeds
  m->d_in_.clear();
  for (size_t id : m->id_in_) {
    m->d_in_.push_back(m->sens_[id]);
    m->sens_[id] = 0;
  }
  // Ensure at least one seed
  casadi_assert(n_known != 0, "No seeds");
  // Allocate result vectors
  m->v_out_.resize(n_unknown);
  m->d_out_.resize(n_unknown);
}

int FmuFunction::eval_ad(FmuMemory* m) const {
  // Number of inputs and outputs
  size_t n_known = m->id_in_.size();
  size_t n_unknown = m->id_out_.size();
  // Quick return if nothing to be calculated
  if (n_unknown == 0) return 0;
  // Evalute (should not be necessary)
  fmi2Status status = fmu_->get_real_(m->c, get_ptr(m->vr_out_), n_unknown, get_ptr(m->v_out_));
  if (status != fmi2OK) {
    casadi_warning("fmi2GetReal failed");
    return 1;
  }
  // Evaluate directional derivatives
  status = fmu_->get_directional_derivative_(m->c, get_ptr(m->vr_out_), n_unknown,
    get_ptr(m->vr_in_), n_known, get_ptr(m->d_in_), get_ptr(m->d_out_));
  if (status != fmi2OK) {
    casadi_warning("fmi2GetDirectionalDerivative failed");
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

int FmuFunction::eval_fd(FmuMemory* m) const {
  // DaeBuilder instance
  auto dae = fmu_->dae();
  // Number of inputs and outputs
  size_t n_known = m->id_in_.size();
  size_t n_unknown = m->id_out_.size();
  // Quick return if nothing to be calculated
  if (n_unknown == 0) return 0;
  // Evalute (should not be necessary)
  fmi2Status status = fmu_->get_real_(m->c, get_ptr(m->vr_out_), n_unknown, get_ptr(m->v_out_));
  if (status != fmi2OK) {
    casadi_warning("fmi2GetReal failed");
    return 1;
  }
  // Get nominal values for outputs
  m->nominal_out_.clear();
  for (size_t id : m->id_out_) m->nominal_out_.push_back(dae->variable(id).nominal);
  // Make outputs dimensionless
  for (size_t k = 0; k < n_unknown; ++k) m->v_out_[k] /= m->nominal_out_[k];
  // Perturbed outputs
  double* yk[4];
  // Error estimate used to update step size
  double u = nan;
  // Initial step size
  double h = step_;
  // For backward, negate h
  if (fd_ == BACKWARD) h = -h;
  // Number of perturbations
  casadi_int n_pert = this->n_pert();
  // Allocate memory for perturbed outputs
  m->fd_out_.resize(n_pert * n_unknown);
  // Perform finite difference algorithm with different step sizes
  for (casadi_int iter = 0; iter < 1 + h_iter_; ++iter) {
    // Calculate all perturbed outputs
    for (casadi_int k = 0; k < n_pert; ++k) {
      // Where to save the perturbed outputs
      yk[k] = &m->fd_out_[n_unknown * k];
      // Get perturbation expression
      double pert;
      switch (fd_) {
      case FORWARD:
      case BACKWARD:
        pert = h;
        break;
      case CENTRAL:
        pert = (2 * static_cast<double>(k) - 1) * h;
        break;
      case SMOOTHING:
        pert = static_cast<double>((2*(k/2)-1) * (k%2+1)) * h;
        break;
      default: casadi_error("Not implemented");
      }
      // Perturb inputs, if allowed
      m->v_pert_.resize(n_known);
      m->in_bounds_.clear();
      for (size_t i = 0; i < n_known; ++i) {
        // Try to take step
        double test = m->v_in_[i] + pert * m->d_in_[i];
        // Check if in bounds
        const Variable& v = dae->variable(m->id_in_[i]);
        bool in_bounds = test >= v.min && test <= v.max;
        // Take step, if allowed
        m->v_pert_[i] = in_bounds ? test : m->v_in_[i];
        // Keep track for later
        m->in_bounds_.push_back(in_bounds);
      }
      // Pass perturbed inputs to FMU
      status = fmu_->set_real_(m->c, get_ptr(m->vr_in_), n_known, get_ptr(m->v_pert_));
      if (status != fmi2OK) {
        casadi_warning("fmi2SetReal failed");
        return 1;
      }
      // Evaluate perturbed FMU
      status = fmu_->get_real_(m->c, get_ptr(m->vr_out_), n_unknown, yk[k]);
      if (status != fmi2OK) {
        casadi_warning("fmi2GetReal failed");
        return 1;
      }
      // Post-process yk[k]
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
          yk[k][i] /= m->nominal_out_[i];
        } else {
          // Input was out of bounds: Discard output
          yk[k][i] = nan;
        }
      }
    }
    // Restore FMU inputs
    status = fmu_->set_real_(m->c, get_ptr(m->vr_in_), n_known, get_ptr(m->v_in_));
    if (status != fmi2OK) {
      casadi_warning("fmi2SetReal failed");
      return 1;
    }
    // FD memory
    casadi_finite_diff_mem<double> fd_mem;
    fd_mem.reltol = reltol_;
    fd_mem.abstol = abstol_;
    fd_mem.smoothing = eps;
    switch (fd_) {
      case FORWARD:
      case BACKWARD:
        u = casadi_forward_diff(yk, get_ptr(m->v_out_),
          get_ptr(m->d_out_), h, n_unknown, &fd_mem);
        break;
      case CENTRAL:
        u = casadi_central_diff(yk, get_ptr(m->v_out_), get_ptr(m->d_out_),
          h, n_unknown, &fd_mem);
        break;
      case SMOOTHING:
        u = casadi_smoothing_diff(yk, get_ptr(m->v_out_), get_ptr(m->d_out_),
          h, n_unknown, &fd_mem);
        break;
      default: casadi_error("Not implemented");
    }
    // Stop, if no more stepsize iterations
    if (iter == h_iter_) break;
    // Update step size
    if (u < 0) {
      // Perturbation failed, try a smaller step size
      h /= u_aim_;
    } else {
      // Update h to get u near the target ratio
      h *= sqrt(u_aim_ / fmax(1., u));
    }
    // Make sure h stays in the range [h_min_,h_max_]
    h = fmin(fmax(h, h_min_), h_max_);
  }
  // Collect requested variables
  for (size_t ind = 0; ind < m->id_out_.size(); ++ind) {
    // Variable id
    size_t id = m->id_out_[ind];
    // Nominal value
    double n = m->nominal_out_[ind];
    // Get the value
    double d_fd = m->d_out_[ind] * n;
    // Use FD instead of AD or to compare with AD
    if (validate_ad_) {
      // Access output variable
      const Variable& v = dae->variable(id);
      // With respect to what variable
      size_t wrt_id = m->wrt_[id];
      const Variable& wrt = dae->variable(wrt_id);
      // Nominal value for input
      double wrt_nom = wrt.nominal;
      // Value to compare with
      double d_ad = m->sens_[id];
      // Magnitude of derivatives
      double d_max = std::fmax(std::fabs(d_fd), std::fabs(d_ad));
      // Check if error exceeds thresholds
      if (d_max > wrt_nom * n * abstol_ && std::fabs(d_ad - d_fd) > d_max * reltol_) {
        // Which element in the vector
        size_t wrt_ind = 0;
        for (; wrt_ind < m->id_in_.size(); ++wrt_ind)
          if (m->id_in_[wrt_ind] == wrt_id) break;
        casadi_assert(wrt_ind < m->id_in_.size(), "Inconsistent variable index for validation");
        // Issue warning
        std::stringstream ss;
        ss << "Inconsistent derivatives of " << v.name << " w.r.t. " << wrt.name << "\n"
          << "At " << m->v_in_[wrt_ind] << ", nominal " << wrt.nominal << ", min " << wrt.min
          << ", max " << wrt.max << ", got " << d_ad
          << " for AD vs. " << d_fd << " for FD[" << to_string(fd_) << "].\n";
        // Also print the stencil:
        std::vector<double> stencil;
        switch (fd_) {
          case FORWARD:
            stencil = {m->v_out_[ind], yk[0][ind]};
            break;
          case BACKWARD:
            stencil = {yk[0][ind], m->v_out_[ind]};
            break;
          case CENTRAL:
            stencil = {yk[0][ind], m->v_out_[ind], yk[1][ind]};
            break;
          case SMOOTHING:
            stencil = {yk[1][ind], yk[0][ind], m->v_out_[ind], yk[2][ind], yk[3][ind]};
            break;
          default: casadi_error("Not implemented");
        }
        // Scale by nominal value
        for (double& s : stencil) s *= n;
        ss << "Values for step size " << h << ", error ratio " << u << ": " << stencil;
        casadi_warning(ss.str());
      }
    } else {
      // Use instead of AD
      m->sens_[id] = d_fd;
    }
  }
  // Successful return
  return 0;
}

int FmuFunction::eval_derivative(FmuMemory* m) const {
  // Gather input and output indices
  gather_sens(m);
  // Calculate derivatives using FMU directional derivative support
  if (enable_ad_) {
    // Evaluate using AD
    if (eval_ad(m)) return 1;
  }
  // Calculate derivatives using finite differences
  if (!enable_ad_ || validate_ad_) {
    // Evaluate using FD
    if (eval_fd(m)) return 1;
  }
  return 0;
}

double FmuFunction::get_sens(FmuMemory* m, size_t id) const {
  return m->sens_.at(id);
}

Fmu::Fmu(const DaeBuilderInternal* dae,
    const std::vector<std::string>& name_in,
    const std::vector<std::string>& name_out,
    const std::map<std::string, std::vector<size_t>>& scheme,
    const std::map<std::string, std::vector<size_t>>& lc)
    : scheme_(scheme), lc_(lc) {
  dae_ = dae->shared_from_this<DaeBuilder>();
  counter_ = 0;
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
  // Mark input indices
  name_in_.clear();
  size_t numel = 0;
  std::vector<bool> lookup(dae->variables_.size(), false);
  for (auto&& n : name_in) {
    auto it = scheme_.find(n);
    if (it == scheme_.end()) continue;
    name_in_.push_back(n);
    for (size_t i : it->second) {
      casadi_assert(!lookup.at(i), "Duplicate variable: " + dae->variables_.at(i).name);
      lookup.at(i) = true;
    }
    numel += it->second.size();
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
  name_out_.clear();
  numel = 0;
  std::fill(lookup.begin(), lookup.end(), false);
  for (auto&& n : name_out) {
    auto it = scheme_.find(n);
    if (it == scheme_.end()) continue;
    name_out_.push_back(n);
    for (size_t i : it->second) {
      casadi_assert(!lookup.at(i), "Duplicate variable: " + dae->variables_.at(i).name);
      lookup.at(i) = true;
    }
    numel += it->second.size();
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
  in_.resize(name_in_.size());
  for (size_t i = 0; i < in_.size(); ++i) {
    auto&& s = scheme_.at(name_in_[i]);
    in_[i].resize(s.size());
    for (size_t k = 0; k < s.size(); ++k) {
      in_[i][k] = iind_map_.at(s[k]);
    }
  }
  // Outputs
  out_.resize(name_out_.size());
  for (size_t i = 0; i < out_.size(); ++i) {
    auto&& s = scheme_.at(name_out_[i]);
    out_[i].resize(s.size());
    for (size_t k = 0; k < s.size(); ++k) {
      out_[i][k] = oind_map_.at(s[k]);
    }
  }
}

size_t Fmu::index_in(const std::string& n) const {
  // Linear search for the input
  for (size_t i = 0; i < name_in_.size(); ++i) {
    if (name_in_[i] == n) return i;
  }
  // Not found
  casadi_error("No such input: " + n);
  return -1;
}

size_t Fmu::index_out(const std::string& n) const {
  // Linear search for the input
  for (size_t i = 0; i < name_out_.size(); ++i) {
    if (name_out_[i] == n) return i;
  }
  // Not found
  casadi_error("No such output: " + n);
  return -1;
}

DaeBuilderInternal* Fmu::dae() const {
  // Get a pointer to the DaeBuilder class
  casadi_assert(dae_.alive(), "DaeBuilder instance has been deleted");
  return static_cast<DaeBuilderInternal*>(dae_->raw_);
}

FmuFunction::FmuFunction(const std::string& name, Fmu* fmu,
    const std::vector<std::string>& name_in,
    const std::vector<std::string>& name_out)
    : FunctionInternal(name), fmu_(fmu) {
  // Parse input IDs
  in_.resize(name_in.size());
  for (size_t k = 0; k < name_in.size(); ++k)
    parse_input(&in_[k], name_in[k]);
  // Parse output IDs
  out_.resize(name_out.size());
  for (size_t k = 0; k < name_out.size(); ++k)
    parse_output(&out_[k], name_out[k]);
  // Set input/output names
  name_in_ = name_in;
  name_out_ = name_out;
  // Default options
  enable_ad_ = fmu_->get_directional_derivative_ != 0;
  validate_ad_ = false;
  step_ = 1e-6;
  abstol_ = 1e-3;
  reltol_ = 1e-3;
  u_aim_ = 100.;
  h_iter_ = 0;
  h_min_ = 0;
  h_max_ = inf;
}

FmuFunction::~FmuFunction() {
  // Free memory
  clear_mem();
  // Decrease reference pointer to Fmu instance
  if (fmu_ && fmu_->counter_-- == 0) delete fmu_;
}

const Options FmuFunction::options_
= {{&FunctionInternal::options_},
   {{"enable_ad",
     {OT_BOOL,
      "Calculate first order derivatives using FMU directional derivative support"}},
    {"validate_ad",
     {OT_BOOL,
      "Compare analytic derivatives with finite differences for validation"}},
    {"step",
     {OT_DOUBLE,
      "Step size, scaled by nominal value"}},
    {"abstol",
     {OT_DOUBLE,
      "Absolute error tolerance, scaled by nominal value"}},
    {"reltol",
     {OT_DOUBLE,
      "Relative error tolerance"}},
    {"h_iter",
     {OT_INT,
      "Number of step size iterations"}},
    {"u_aim",
     {OT_DOUBLE,
      "Target ratio of truncation error to roundoff error"}},
    {"h_min",
     {OT_DOUBLE,
      "Minimum step size"}},
    {"h_max",
     {OT_DOUBLE,
      "Maximum step size"}}
   }
};

void FmuFunction::init(const Dict& opts) {
  // Read options
  for (auto&& op : opts) {
    if (op.first=="enable_ad") {
      enable_ad_ = op.second;
    } else if (op.first=="validate_ad") {
      validate_ad_ = op.second;
    } else if (op.first=="step") {
      step_ = op.second;
    } else if (op.first=="abstol") {
      abstol_ = op.second;
    } else if (op.first=="reltol") {
      reltol_ = op.second;
    } else if (op.first=="h_iter") {
      h_iter_ = op.second;
    } else if (op.first=="u_aim") {
      u_aim_ = op.second;
    } else if (op.first=="h_min") {
      h_min_ = op.second;
    } else if (op.first=="h_max") {
      h_max_ = op.second;
    }
  }

  // Call the initialization method of the base class
  FunctionInternal::init(opts);

  // Read FD mode
  fd_ = to_enum<FmuFunction::FdMode>(fd_method_, "forward");

  // Get a pointer to the DaeBuilder class
  auto dae = fmu_->dae();

  // Consistency checks
  if (enable_ad_) casadi_assert(dae->provides_directional_derivative_,
    "FMU does not provide support for analytic derivatives");
  if (validate_ad_ && !enable_ad_) casadi_error("Inconsistent options");

  // Number of variables
  nv_ = dae->variables_.size();

  // Collect all inputs in any Jacobian or adjoint block
  std::vector<bool> in_jac(nv_, false);
  for (auto&& i : out_) {
    if (i.type == JAC_OUTPUT) {
      for (size_t j : fmu_->in_[i.wrt]) in_jac[fmu_->iind_[j]] = true;
    } else if (i.type == ADJ_SENS) {
      for (size_t j : fmu_->in_[i.ind]) in_jac[fmu_->iind_[j]] = true;
    }
  }
  jac_in_.clear();
  for (size_t k = 0; k < in_jac.size(); ++k) {
    if (in_jac[k]) jac_in_.push_back(k);
  }

  // Do the outputs include Jacobian blocks and/or adjoint directional derivatives
  has_adj_ = has_jac_ = false;

  // Collect all outputs in any Jacobian or adjoint block
  std::fill(in_jac.begin(), in_jac.end(), false);
  for (auto&& i : out_) {
    if (i.type == JAC_OUTPUT) {
      for (size_t j : fmu_->out_[i.ind]) in_jac[fmu_->oind_[j]] = true;
      has_jac_ = true;
    }
  }
  for (auto&& i : in_) {
    if (i.type == ADJ_SEED) {
      for (size_t j : fmu_->out_[i.ind]) in_jac[fmu_->oind_[j]] = true;
      has_adj_ = true;
    }
  }
  jac_out_.clear();
  for (size_t k = 0; k < in_jac.size(); ++k) {
    if (in_jac[k]) jac_out_.push_back(k);
  }

  // Get sparsity pattern for extended Jacobian
  sp_ext_ = dae->jac_sparsity(jac_out_, jac_in_);

  // Calculate graph coloring
  coloring_ = sp_ext_.uni_coloring();
  if (verbose_) casadi_message("Graph coloring: " + str(sp_ext_.size2())
    + " -> " + str(coloring_.size2()) + " directions");

  // Work vectors
  if (has_adj_) alloc_w(nv_, false);  // adjw
}

void FmuFunction::parse_input(InputStruct* s, const std::string& n) const {
  // Look for prefix
  if (has_prefix(n)) {
    // Get the prefix
    std::string pref, rem;
    pref = pop_prefix(n, &rem);
    if (pref == "out") {
      // Nondifferentiated function output (unused)
      s->type = DUMMY_OUTPUT;
      s->ind = fmu_->index_out(rem);
    } else if (pref == "adj") {
      // Adjoint seed
      s->type = ADJ_SEED;
      s->ind = fmu_->index_out(rem);
    } else {
      // No such prefix
      casadi_error("No such prefix: " + pref);
    }
  } else {
    // No prefix - regular input
    s->type = REG_INPUT;
    s->ind = fmu_->index_in(n);
  }
}

void FmuFunction::parse_output(OutputStruct* s, const std::string& n) const {
  // Look for prefix
  if (has_prefix(n)) {
    // Get the prefix
    std::string part1, rem;
    part1 = pop_prefix(n, &rem);
    if (part1 == "jac") {
      // Jacobian block
      casadi_assert(has_prefix(rem), "Two arguments expected for Jacobian block");
      std::string part2 = pop_prefix(rem, &rem);
      s->type = JAC_OUTPUT;
      s->ind = fmu_->index_out(part2);
      s->wrt = fmu_->index_in(rem);
    } else if (part1 == "adj") {
      // Adjoint sensitivity
      s->type = ADJ_SENS;
      s->ind = fmu_->index_in(rem);
    } else {
      // No such prefix
      casadi_error("No such prefix: " + part1);
    }
  } else {
    // No prefix - regular output
    s->type = REG_OUTPUT;
    s->ind = fmu_->index_out(n);
  }
}

Sparsity FmuFunction::get_sparsity_in(casadi_int i) {
  switch (in_.at(i).type) {
    case REG_INPUT:
      return Sparsity::dense(fmu_->in_[in_.at(i).ind].size(), 1);
    case ADJ_SEED:
      return Sparsity::dense(fmu_->out_[in_.at(i).ind].size(), 1);
    case DUMMY_OUTPUT:
      return Sparsity(fmu_->out_[in_.at(i).ind].size(), 1);
  }
  return Sparsity();
}

Sparsity FmuFunction::get_sparsity_out(casadi_int i) {
  switch (out_.at(i).type) {
    case REG_OUTPUT:
      return Sparsity::dense(fmu_->out_[out_.at(i).ind].size(), 1);
    case ADJ_SENS:
      return Sparsity::dense(fmu_->in_[out_.at(i).ind].size(), 1);
    case JAC_OUTPUT:
      return fmu_->dae()->jac_sparsity(fmu_->get_out(out_.at(i).ind),
        fmu_->get_in(out_.at(i).wrt));
  }
  return Sparsity();
}

int FmuFunction::eval(const double** arg, double** res, casadi_int* iw, double* w,
    void* mem) const {
  // Get memory struct
  FmuMemory* m = static_cast<FmuMemory*>(mem);
  casadi_assert(m != 0, "Memory is null");

  // Work vector holding the adjoint seeds and sensitivities
  double* adjw = 0;
  if (has_adj_) {
    adjw = w;
    w += nv_;
  }

  // DaeBuilder instance
  auto dae = fmu_->dae();
  // Pass all regular inputs
  for (size_t k = 0; k < in_.size(); ++k) {
    if (in_[k].type == REG_INPUT) {
      const std::vector<size_t>& iind = fmu_->in_[in_[k].ind];
      for (size_t i = 0; i < iind.size(); ++i) {
        set(m, fmu_->iind_[iind[i]], arg[k] ? arg[k][i] : 0);
      }
    }
  }
  // Request all regular outputs to be evaluated
  for (size_t k = 0; k < out_.size(); ++k) {
    if (res[k] && out_[k].type == REG_OUTPUT) {
      const std::vector<size_t>& oind = fmu_->out_[out_[k].ind];
      for (size_t i = 0; i < oind.size(); ++i) {
        request(m, fmu_->oind_[oind[i]]);
      }
    }
  }
  // Evaluate
  if (eval(m)) return 1;
  // Get regular outputs
  for (size_t k = 0; k < out_.size(); ++k) {
    if (res[k] && out_[k].type == REG_OUTPUT) {
      const std::vector<size_t>& oind = fmu_->out_[out_[k].ind];
      for (size_t i = 0; i < oind.size(); ++i) {
        get(m, fmu_->oind_[oind[i]], &res[k][i]);
      }
    }
  }

  // What other blocks are there?
  bool need_jac = false;
  for (size_t k = 0; k < out_.size(); ++k) {
    if (res[k] && (out_[k].type == JAC_OUTPUT || out_[k].type == ADJ_SENS)) {
     need_jac = true;
    }
  }

  // Evalute Jacobian blocks
  if (need_jac) {
    // Copy adjoint seeds to adjw
    if (has_adj_) {
      std::fill(adjw, adjw + nv_, 0);
      for (size_t i = 0; i < in_.size(); ++i) {
        if (arg[i] && in_[i].type == ADJ_SEED) {
          const std::vector<size_t>& oind = fmu_->out_[in_[i].ind];
          for (size_t k = 0; k < oind.size(); ++k) adjw[fmu_->oind_[oind[k]]] = arg[i][k];
        }
      }
    }
    // Loop over colors
    for (casadi_int c = 0; c < coloring_.size2(); ++c) {
      // Loop over input indices for color
      for (casadi_int kc = coloring_.colind(c); kc < coloring_.colind(c + 1); ++kc) {
        casadi_int vin = coloring_.row(kc);
        // Differentiation with respect to what variable
        size_t Jc = jac_in_.at(vin);
        // Nominal value
        double nom = dae->variable(Jc).nominal;
        // Set seed for column
        set_seed(m, Jc, nom);
        // Request corresponding outputs
        for (casadi_int Jk = sp_ext_.colind(vin); Jk < sp_ext_.colind(vin + 1); ++Jk) {
          casadi_int vout = sp_ext_.row(Jk);
          request(m, jac_out_.at(vout), Jc);
        }
      }
      // Calculate derivatives
      if (eval_derivative(m)) return 1;
      // Loop over input indices for color
      for (casadi_int kc = coloring_.colind(c); kc < coloring_.colind(c + 1); ++kc) {
        casadi_int vin = coloring_.row(kc);
        // Differentiation with respect to what variable
        size_t Jc = jac_in_.at(vin);
        // Inverse of nominal value
        double inv_nom = 1. / dae->variable(Jc).nominal;
        // Fetch Jacobian blocks
        for (size_t k = 0; k < out_.size(); ++k) {
          if (res[k] && out_[k].type == JAC_OUTPUT) {
            // Find input index
            const std::vector<size_t>& iind = fmu_->in_[out_[k].wrt];
            for (size_t Bc = 0; Bc < iind.size(); ++Bc) {
              if (fmu_->iind_[iind[Bc]] == Jc) {
                // Column exists in Jacobian block
                const Sparsity& sp = sparsity_out(k);
                const std::vector<size_t>& oind = fmu_->out_[out_[k].ind];
                for (casadi_int Bk = sp.colind(Bc); Bk < sp.colind(Bc + 1); ++Bk) {
                  // Save Jacobian nonzero, scaled by nominal value factor
                  res[k][Bk] = inv_nom * get_sens(m, fmu_->oind_[oind.at(sp.row(Bk))]);
                }
              }
            }
          }
        }
        // Propagate adjoint derivatives
        if (has_adj_) {
          for (casadi_int Jk = sp_ext_.colind(vin); Jk < sp_ext_.colind(vin + 1); ++Jk) {
            casadi_int vout = sp_ext_.row(Jk);
            size_t Jr = jac_out_.at(vout);
            adjw[Jc] += inv_nom * get_sens(m, Jr) * adjw[Jr];
          }
        }
      }
    }
    // Collect adjoint sensitivities
    if (has_adj_) {
      for (size_t i = 0; i < out_.size(); ++i) {
        if (res[i] && out_[i].type == ADJ_SENS) {
          const std::vector<size_t>& iind = fmu_->in_[out_[i].ind];
          for (size_t k = 0; k < iind.size(); ++k) res[i][k] = adjw[fmu_->iind_[iind[k]]];
        }
      }
    }
  }
  // Successful return
  return 0;
}

casadi_int FmuFunction::n_pert() const {
  switch (fd_) {
  case FmuFunction::FORWARD:
  case FmuFunction::BACKWARD:
    return 1;
  case FmuFunction::CENTRAL:
    return 2;
  case FmuFunction::SMOOTHING:
    return 4;
  default: break;
  }
  casadi_error("Not implemented");
}

std::string to_string(FmuFunction::FdMode v) {
  switch (v) {
  case FmuFunction::FORWARD: return "forward";
  case FmuFunction::BACKWARD: return "backward";
  case FmuFunction::CENTRAL: return "central";
  case FmuFunction::SMOOTHING: return "smoothing";
  default: break;
  }
  return "";
}

bool FmuFunction::has_prefix(const std::string& s) {
  return s.find('_') < s.size();
}

std::string FmuFunction::pop_prefix(const std::string& s, std::string* rem) {
  // Get prefix
  casadi_assert_dev(!s.empty());
  size_t pos = s.find('_');
  casadi_assert(pos < s.size(), "Cannot process \"" + s + "\"");
  // Get prefix
  std::string r = s.substr(0, pos);
  // Remainder, if requested (note that rem == &s is possible)
  if (rem) *rem = s.substr(pos+1, std::string::npos);
  // Return prefix
  return r;
}

Function FmuFunction::get_jacobian(const std::string& name, const std::vector<std::string>& inames,
    const std::vector<std::string>& onames, const Dict& opts) const {
  // Return value
  Function ret;
  ret.own(new FmuFunction(name, fmu_, inames, onames));
  // Hack: Manually enable finite differenting (pending implementation in class)
  Dict opts2 = opts;
  opts2["enable_fd"] = true;
  ret->construct(opts2);
  return ret;
}

Function FmuFunction::get_reverse(casadi_int nadj, const std::string& name,
    const std::vector<std::string>& inames,
    const std::vector<std::string>& onames,
    const Dict& opts) const {
  // Only single directional derivative implemented
  casadi_assert(nadj == 1, "Not implemented");
  // Return value
  Function ret;
  ret.own(new FmuFunction(name, fmu_, inames, onames));
  // Hack: Manually enable finite differenting (pending implementation in class)
  Dict opts2 = opts;
  opts2["enable_fd"] = true;
  ret->construct(opts2);
  return ret;
}

bool FmuFunction::has_jac_sparsity(casadi_int oind, casadi_int iind) const {
  return in_.at(iind).type == REG_INPUT && out_.at(oind).type == REG_OUTPUT;
}

Sparsity FmuFunction::get_jac_sparsity(casadi_int oind, casadi_int iind,
    bool symmetric) const {
  // DaeBuilder instance
  auto dae = fmu_->dae();
  // Get the indices
  std::vector<size_t> oi = fmu_->out_[out_.at(oind).ind];
  std::vector<size_t> ii = fmu_->in_[in_.at(iind).ind];
  // Convert to indices in FMU
  for (size_t& i : oi) i = fmu_->oind_.at(i);
  for (size_t& i : ii) i = fmu_->iind_.at(i);
  // Get the Jacobian block
  return dae->jac_sparsity(oi, ii);
}

#endif  // WITH_FMU

} // namespace casadi

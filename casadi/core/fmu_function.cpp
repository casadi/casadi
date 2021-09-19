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

Fmu::~Fmu() {
  // Free all memory objects
  for (int id = 0; id < mem_.size(); ++id) {
    if (mem_[id].c && free_instance_) free_instance_(mem_[id].c);
  }
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
  get_real_ = reinterpret_cast<fmi2GetRealTYPE*>(get_function("fmi2GetReal"));
  set_real_ = reinterpret_cast<fmi2SetRealTYPE*>(get_function("fmi2SetReal"));
  get_integer_ = reinterpret_cast<fmi2GetIntegerTYPE*>(get_function("fmi2GetInteger"));
  set_integer_ = reinterpret_cast<fmi2SetIntegerTYPE*>(get_function("fmi2SetInteger"));
  get_boolean_ = reinterpret_cast<fmi2GetBooleanTYPE*>(get_function("fmi2GetBoolean"));
  set_boolean_ = reinterpret_cast<fmi2SetBooleanTYPE*>(get_function("fmi2SetBoolean"));
  get_string_ = reinterpret_cast<fmi2GetStringTYPE*>(get_function("fmi2GetString"));
  set_string_ = reinterpret_cast<fmi2SetStringTYPE*>(get_function("fmi2SetString"));
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

fmi2Component Fmu::instantiate() {
  fmi2String instanceName = self_.model_identifier_.c_str();
  fmi2Type fmuType = fmi2ModelExchange;
  fmi2String fmuGUID = self_.guid_.c_str();
  fmi2String fmuResourceLocation = resource_loc_.c_str();
  fmi2Boolean visible = fmi2False;
  fmi2Boolean loggingOn = self_.debug_;
  fmi2Component c = instantiate_(instanceName, fmuType, fmuGUID, fmuResourceLocation,
    &functions_, visible, loggingOn);
  if (c == 0) casadi_error("fmi2Instantiate failed");
  return c;
}

int Fmu::checkout() {
  // Reuse an element of the memory pool if possible
  int mem;
  for (mem = 0; mem < mem_.size(); ++mem) {
    if (!mem_[mem].in_use) break;
  }
  // Add new element necessary
  if (mem == mem_.size()) mem_.push_back(Memory());
  // Memory object
  Memory& m = mem_.at(mem);
  // Mark as in use
  casadi_assert(!m.in_use, "Memory object is already in use");
  m.in_use = true;
  // Make sure not already instantiated
  casadi_assert(m.c == 0, "Already instantiated");
  // Create instance
  m.c = instantiate();
  // Reset solver
  setup_experiment(mem);
  // Always initialize before first call
  m.need_init = true;
  // Allocate/reset value buffer
  m.buffer_.resize(self_.variables_.size());
  std::fill(m.buffer_.begin(), m.buffer_.end(), casadi::nan);
  // Allocate/reset sensitivities
  m.sens_.resize(self_.variables_.size());
  std::fill(m.sens_.begin(), m.sens_.end(), 0);
  // Allocate/reset changed
  m.changed_.resize(self_.variables_.size());
  std::fill(m.changed_.begin(), m.changed_.end(), false);
  // Allocate/reset requested
  m.requested_.resize(self_.variables_.size());
  std::fill(m.requested_.begin(), m.requested_.end(), false);
  // Also allocate memory for corresponding Jacobian entry (for debugging)
  m.wrt_.resize(self_.variables_.size());
  // Return memory object
  return mem;
}

fmi2Component Fmu::memory(int mem) {
  return mem_.at(mem).c;
}

fmi2Component Fmu::pop_memory(int mem) {
  fmi2Component c = mem_.at(mem).c;
  mem_.at(mem).c = 0;
  return c;
}

void Fmu::release(int mem) {
  // Mark as in not in use
  if (mem >= 0) {
    if (!mem_.at(mem).in_use) casadi_warning("Memory object not in use");
    if (mem_.at(mem).c && free_instance_) {
      free_instance_(mem_.at(mem).c);
      mem_.at(mem).c = nullptr;
    }
    mem_.at(mem).in_use = false;
  }
}

void Fmu::setup_experiment(int mem) {
  fmi2Status status = setup_experiment_(memory(mem), self_.fmutol_ > 0, self_.fmutol_, 0.,
    fmi2True, 1.);
  casadi_assert(status == fmi2OK, "fmi2SetupExperiment failed");
}

int Fmu::reset(int mem) {
  fmi2Status status = reset_(memory(mem));
  if (status != fmi2OK) {
    casadi_warning("fmi2Reset failed");
    return 1;
  }
  return 0;
}

int Fmu::enter_initialization_mode(int mem) {
  fmi2Status status = enter_initialization_mode_(memory(mem));
  if (status != fmi2OK) {
    casadi_warning("fmi2EnterInitializationMode failed: " + str(status));
    return 1;
  }
  return 0;
}

int Fmu::exit_initialization_mode(int mem) {
  fmi2Status status = exit_initialization_mode_(memory(mem));
  if (status != fmi2OK) {
    casadi_warning("fmi2ExitInitializationMode failed");
    return 1;
  }
  return 0;
}

void Fmu::set(int mem, size_t id, double value) {
  // Get memory
  Memory& m = mem_.at(mem);
  // Update buffer
  if (value != m.buffer_.at(id)) {
    m.buffer_.at(id) = value;
    m.changed_.at(id) = true;
  }
}

void Fmu::set_seed(int mem, size_t id, double value) {
  // Get memory
  Memory& m = mem_.at(mem);
  // Update buffer
  if (value != 0) {
    m.sens_.at(id) = value;
    m.changed_.at(id) = true;
  }
}

void Fmu::request(int mem, size_t id, size_t wrt_id) {
  // Get memory
  Memory& m = mem_.at(mem);
  // Mark as requested
  m.requested_.at(id) = true;
  // Also log corresponding input index
  m.wrt_.at(id) = wrt_id;
}

int Fmu::eval(int mem, const FmuFunction& f) {
  // Gather inputs and outputs
  gather_io(mem);
  // Get memory
  Memory& m = mem_.at(mem);
  // Number of inputs and outputs
  size_t n_set = m.id_in_.size();
  size_t n_out = m.id_out_.size();
  // Fmi return flag
  fmi2Status status;
  // Initialize, if needed
  if (m.need_init) {
    // DaeBuilder instance
    casadi_assert(f.dae_.alive(), "DaeBuilder instance has been deleted");
    auto dae = static_cast<DaeBuilderInternal*>(f.dae_->raw_);
    // Set all variables before initialization
    for (const Variable& v : dae->variables_) {
      continue;
      // If nan - variable has not been set - keep default value
      if (std::isnan(v.value)) continue;
      // Convert to expected type
      fmi2ValueReference vr = v.value_reference;
      // Get value
      switch (v.type) {
        case Variable::REAL:
          {
            // Real
            fmi2Real value = v.value;
            status = set_real_(m.c, &vr, 1, &value);
            if (status != fmi2OK) {
              casadi_warning("fmi2SetReal failed for " + v.name + ", value: " + str(v.value));
              return 1;
            }
            break;
          }
        case Variable::INTEGER:
        case Variable::ENUM:
          {
            // Integer: Convert to double
            fmi2Integer value = static_cast<casadi_int>(v.value);
            status = set_integer_(m.c, &vr, 1, &value);
            if (status != fmi2OK) {
              casadi_warning("fmi2SetInteger failed for " + v.name + ", value: " + str(v.value));
              return 1;
            }
            break;
          }
        case Variable::BOOLEAN:
          {
            // Boolean: Convert to double
            fmi2Boolean value = static_cast<casadi_int>(v.value);
            status = set_boolean_(m.c, &vr, 1, &value);
            if (status != fmi2OK) {
              casadi_warning("fmi2SetBoolean failed for " + v.name + ", value: " + str(v.value));
              return 1;
            }
            break;
          }
        case Variable::STRING:
          {
            // String
            fmi2String value = v.stringvalue.c_str();
            status = set_string_(m.c, &vr, 1, &value);
            if (status != fmi2OK) {
              casadi_warning("fmi2SetString failed for " + v.name + ", value: " + v.stringvalue);
              return 1;
            }
            break;
          }
        default:
          casadi_warning("Ignoring " + v.name + ", type: " + to_string(v.type));
      }
    }
    // Set all variables before initialization
    status = set_real_(m.c, get_ptr(m.vr_in_), n_set, get_ptr(m.v_in_));
    if (status != fmi2OK) {
      casadi_warning("fmi2SetReal failed");
      return 1;
    }
    // Initialization mode begins
    if (enter_initialization_mode(mem)) return 1;
    // Get all values
    for (Variable& v : dae->variables_) {
      // Convert to expected type
      fmi2ValueReference vr = v.value_reference;
      // Get value
      switch (v.type) {
        case Variable::REAL:
          {
            // Real
            fmi2Real value;
            status = get_real_(m.c, &vr, 1, &value);
            if (status != fmi2OK) {
              casadi_warning("fmi2GetReal failed for " + v.name);
              return 1;
            }
            v.value = value;
            break;
          }
        case Variable::INTEGER:
        case Variable::ENUM:
          {
            // Integer: Convert to double
            fmi2Integer value;
            status = get_integer_(m.c, &vr, 1, &value);
            if (status != fmi2OK) {
              casadi_warning("fmi2GetInteger failed for " + v.name);
              return 1;
            }
            v.value = value;
            break;
          }
        case Variable::BOOLEAN:
          {
            // Boolean: Convert to double
            fmi2Boolean value;
            status = get_boolean_(m.c, &vr, 1, &value);
            if (status != fmi2OK) {
              casadi_warning("fmi2GetBoolean failed for " + v.name);
              return 1;
            }
            v.value = value;
            break;
          }
        case Variable::STRING:
          {
            // String
            fmi2String value;
            status = get_string_(m.c, &vr, 1, &value);
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
    // Initialization mode ends
    if (exit_initialization_mode(mem)) return 1;
    // Initialized
    m.need_init = false;
  }
  // Set all variables
  status = set_real_(m.c, get_ptr(m.vr_in_), n_set, get_ptr(m.v_in_));
  if (status != fmi2OK) {
    casadi_warning("fmi2SetReal failed");
    return 1;
  }
  // Quick return if nothing requested
  if (n_out == 0) return 0;
  // Calculate all variables
  m.v_out_.resize(n_out);
  status = get_real_(m.c, get_ptr(m.vr_out_), n_out, get_ptr(m.v_out_));
  if (status != fmi2OK) {
    casadi_warning("fmi2GetReal failed");
    return 1;
  }
  // Collect requested variables
  auto it = m.v_out_.begin();
  for (size_t id : m.id_out_) {
    m.buffer_[id] = *it++;
  }
  // Successful return
  return 0;
}

void Fmu::get(int mem, size_t id, double* value) {
  // Get memory
  Memory& m = mem_.at(mem);
  // Save to return
  *value = m.buffer_.at(id);
}

void Fmu::gather_io(int mem) {
  // Get memory
  Memory& m = mem_.at(mem);
  // Collect input indices and corresponding value references and values
  m.id_in_.clear();
  m.vr_in_.clear();
  m.v_in_.clear();
  for (size_t id = 0; id < m.changed_.size(); ++id) {
    if (m.changed_[id]) {
      const Variable& v = self_.variable(id);
      m.id_in_.push_back(id);
      m.vr_in_.push_back(v.value_reference);
      m.v_in_.push_back(m.buffer_[id]);
      m.changed_[id] = false;
    }
  }
  // Collect output indices, corresponding value references
  m.id_out_.clear();
  m.vr_out_.clear();
  for (size_t id = 0; id < m.requested_.size(); ++id) {
    if (m.requested_[id]) {
      const Variable& v = self_.variable(id);
      m.id_out_.push_back(id);
      m.vr_out_.push_back(v.value_reference);
      m.requested_[id] = false;
    }
  }
}

void Fmu::gather_sens(int mem) {
  // Gather input and output indices
  gather_io(mem);
  // Get memory
  Memory& m = mem_.at(mem);
  // Number of inputs and outputs
  size_t n_known = m.id_in_.size();
  size_t n_unknown = m.id_out_.size();
  // Get/clear seeds
  m.d_in_.clear();
  for (size_t id : m.id_in_) {
    m.d_in_.push_back(m.sens_[id]);
    m.sens_[id] = 0;
  }
  // Ensure at least one seed
  casadi_assert(n_known != 0, "No seeds");
  // Allocate result vectors
  m.v_out_.resize(n_unknown);
  m.d_out_.resize(n_unknown);
}

int Fmu::eval_ad(int mem, const FmuFunction& f) {
  // Get memory
  Memory& m = mem_.at(mem);
  // Number of inputs and outputs
  size_t n_known = m.id_in_.size();
  size_t n_unknown = m.id_out_.size();
  // Quick return if nothing to be calculated
  if (n_unknown == 0) return 0;
  // Evalute (should not be necessary)
  fmi2Status status = get_real_(m.c, get_ptr(m.vr_out_), n_unknown, get_ptr(m.v_out_));
  if (status != fmi2OK) {
    casadi_warning("fmi2GetReal failed");
    return 1;
  }
  // Evaluate directional derivatives
  status = get_directional_derivative_(m.c, get_ptr(m.vr_out_), n_unknown,
    get_ptr(m.vr_in_), n_known, get_ptr(m.d_in_), get_ptr(m.d_out_));
  if (status != fmi2OK) {
    casadi_warning("fmi2GetDirectionalDerivative failed");
    return 1;
  }
  // Collect requested variables
  auto it = m.d_out_.begin();
  for (size_t id : m.id_out_) {
    m.sens_[id] = *it++;
  }
  // Successful return
  return 0;
}

int Fmu::eval_fd(int mem, const FmuFunction& f) {
  // Get memory
  Memory& m = mem_.at(mem);
  // Number of inputs and outputs
  size_t n_known = m.id_in_.size();
  size_t n_unknown = m.id_out_.size();
  // Quick return if nothing to be calculated
  if (n_unknown == 0) return 0;
  // Evalute (should not be necessary)
  fmi2Status status = get_real_(m.c, get_ptr(m.vr_out_), n_unknown, get_ptr(m.v_out_));
  if (status != fmi2OK) {
    casadi_warning("fmi2GetReal failed");
    return 1;
  }
  // Get nominal values for outputs
  m.nominal_out_.clear();
  for (size_t id : m.id_out_) m.nominal_out_.push_back(self_.variable(id).nominal);
  // Make outputs dimensionless
  for (size_t k = 0; k < n_unknown; ++k) m.v_out_[k] /= m.nominal_out_[k];
  // Perturbed outputs
  double* yk[4];
  // Error estimate used to update step size
  double u = nan;
  // Initial step size
  double h = f.step_;
  // For backward, negate h
  if (f.fd_ == BACKWARD) h = -h;
  // Number of perturbations
  casadi_int n_pert = f.n_pert();
  // Allocate memory for perturbed outputs
  m.fd_out_.resize(n_pert * n_unknown);
  // Perform finite difference algorithm with different step sizes
  for (casadi_int iter = 0; iter < 1 + f.h_iter_; ++iter) {
    // Calculate all perturbed outputs
    for (casadi_int k = 0; k < n_pert; ++k) {
      // Where to save the perturbed outputs
      yk[k] = &m.fd_out_[n_unknown * k];
      // Get perturbation expression
      double pert;
      switch (f.fd_) {
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
      m.v_pert_.resize(n_known);
      m.in_bounds_.clear();
      for (size_t i = 0; i < n_known; ++i) {
        // Try to take step
        double test = m.v_in_[i] + pert * m.d_in_[i];
        // Check if in bounds
        const Variable& v = self_.variable(m.id_in_[i]);
        bool in_bounds = test >= v.min && test <= v.max;
        // Take step, if allowed
        m.v_pert_[i] = in_bounds ? test : m.v_in_[i];
        // Keep track for later
        m.in_bounds_.push_back(in_bounds);
      }
      // Pass perturbed inputs to FMU
      status = set_real_(m.c, get_ptr(m.vr_in_), n_known, get_ptr(m.v_pert_));
      if (status != fmi2OK) {
        casadi_warning("fmi2SetReal failed");
        return 1;
      }
      // Evaluate perturbed FMU
      status = get_real_(m.c, get_ptr(m.vr_out_), n_unknown, yk[k]);
      if (status != fmi2OK) {
        casadi_warning("fmi2GetReal failed");
        return 1;
      }
      // Post-process yk[k]
      for (size_t i = 0; i < n_unknown; ++i) {
        // Variable id
        size_t id = m.id_out_[i];
        // Differentiation with respect to what variable
        size_t wrt_id = m.wrt_.at(id);
        // Find the corresponding input variable
        size_t wrt_i;
        for (wrt_i = 0; wrt_i < n_known; ++wrt_i) {
          if (m.id_in_[wrt_i] == wrt_id) break;
        }
        // Check if in bounds
        if (m.in_bounds_.at(wrt_i)) {
          // Input was in bounds: Keep output, make dimensionless
          yk[k][i] /= m.nominal_out_[i];
        } else {
          // Input was out of bounds: Discard output
          yk[k][i] = nan;
        }
      }
    }
    // Restore FMU inputs
    status = set_real_(m.c, get_ptr(m.vr_in_), n_known, get_ptr(m.v_in_));
    if (status != fmi2OK) {
      casadi_warning("fmi2SetReal failed");
      return 1;
    }
    // FD memory
    casadi_finite_diff_mem<double> fd_mem;
    fd_mem.reltol = f.reltol_;
    fd_mem.abstol = f.abstol_;
    fd_mem.smoothing = eps;
    switch (f.fd_) {
      case FORWARD:
      case BACKWARD:
        u = casadi_forward_diff(yk, get_ptr(m.v_out_),
          get_ptr(m.d_out_), h, n_unknown, &fd_mem);
        break;
      case CENTRAL:
        u = casadi_central_diff(yk, get_ptr(m.v_out_), get_ptr(m.d_out_),
          h, n_unknown, &fd_mem);
        break;
      case SMOOTHING:
        u = casadi_smoothing_diff(yk, get_ptr(m.v_out_), get_ptr(m.d_out_),
          h, n_unknown, &fd_mem);
        break;
      default: casadi_error("Not implemented");
    }
    // Stop, if no more stepsize iterations
    if (iter == f.h_iter_) break;
    // Update step size
    if (u < 0) {
      // Perturbation failed, try a smaller step size
      h /= f.u_aim_;
    } else {
      // Update h to get u near the target ratio
      h *= sqrt(f.u_aim_ / fmax(1., u));
    }
    // Make sure h stays in the range [h_min_,h_max_]
    h = fmin(fmax(h, f.h_min_), f.h_max_);
  }
  // Collect requested variables
  for (size_t ind = 0; ind < m.id_out_.size(); ++ind) {
    // Variable id
    size_t id = m.id_out_[ind];
    // Nominal value
    double n = m.nominal_out_[ind];
    // Get the value
    double d_fd = m.d_out_[ind] * n;
    // Use FD instead of AD or to compare with AD
    if (f.validate_ad_) {
      // Access output variable
      const Variable& v = self_.variable(id);
      // With respect to what variable
      size_t wrt_id = m.wrt_[id];
      const Variable& wrt = self_.variable(wrt_id);
      // Nominal value for input
      double wrt_nom = wrt.nominal;
      // Value to compare with
      double d_ad = m.sens_[id];
      // Magnitude of derivatives
      double d_max = std::fmax(std::fabs(d_fd), std::fabs(d_ad));
      // Check if error exceeds thresholds
      if (d_max > wrt_nom * n * f.abstol_ && std::fabs(d_ad - d_fd) > d_max * f.reltol_) {
        // Which element in the vector
        size_t wrt_ind = 0;
        for (; wrt_ind < m.id_in_.size(); ++wrt_ind)
          if (m.id_in_[wrt_ind] == wrt_id) break;
        casadi_assert(wrt_ind < m.id_in_.size(), "Inconsistent variable index for validation");
        // Issue warning
        std::stringstream ss;
        ss << "Inconsistent derivatives of " << v.name << " w.r.t. " << wrt.name << "\n"
          << "At " << m.v_in_[wrt_ind] << ", nominal " << wrt.nominal << ", min " << wrt.min
          << ", max " << wrt.max << ", got " << d_ad
          << " for AD vs. " << d_fd << " for FD[" << to_string(f.fd_) << "].\n";
        // Also print the stencil:
        std::vector<double> stencil;
        switch (f.fd_) {
          case FORWARD:
            stencil = {m.v_out_[ind], yk[0][ind]};
            break;
          case BACKWARD:
            stencil = {yk[0][ind], m.v_out_[ind]};
            break;
          case CENTRAL:
            stencil = {yk[0][ind], m.v_out_[ind], yk[1][ind]};
            break;
          case SMOOTHING:
            stencil = {yk[1][ind], yk[0][ind], m.v_out_[ind], yk[2][ind], yk[3][ind]};
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
      m.sens_[id] = d_fd;
    }
  }
  // Successful return
  return 0;
}

int Fmu::eval_derivative(int mem, const FmuFunction& f) {
  // Gather input and output indices
  gather_sens(mem);
  // Calculate derivatives using FMU directional derivative support
  if (f.enable_ad_) {
    // Evaluate using AD
    if (eval_ad(mem, f)) return 1;
  }
  // Calculate derivatives using finite differences
  if (!f.enable_ad_ || f.validate_ad_) {
    // Evaluate using FD
    if (eval_fd(mem, f)) return 1;
  }
  return 0;
}

void Fmu::get_sens(int mem, size_t id, double* value) {
  // Get memory
  Memory& m = mem_.at(mem);
  // Save to return
  *value = m.sens_.at(id);
}

int Fmu::eval(int mem, const double** arg, double** res, const FmuFunction& f) {
  // Pass all regular inputs
  for (size_t k = 0; k < f.in_.size(); ++k) {
    if (f.in_[k]->is_reg()) {
      for (size_t i = 0; i < f.in_[k]->size(); ++i) {
        set(mem, f.in_[k]->ind(i), arg[k] ? arg[k][i] : 0);
      }
    }
  }
  // Request all regular outputs to be evaluated
  for (size_t k = 0; k < f.out_.size(); ++k) {
    if (res[k] && f.out_[k]->is_reg()) {
      for (size_t i = 0; i < f.out_[k]->size(); ++i) {
        request(mem, f.out_[k]->ind(i));
      }
    }
  }
  // Evaluate
  if (eval(mem, f)) return 1;
  // Get regular outputs
  for (size_t k = 0; k < f.out_.size(); ++k) {
    if (res[k] && f.out_[k]->is_reg()) {
      for (size_t i = 0; i < f.out_[k]->size(); ++i) {
        get(mem, f.out_[k]->ind(i), &res[k][i]);
      }
    }
  }
  // What other blocks are there?
  bool any_jac = false;
  for (size_t k = 0; k < f.out_.size(); ++k) {
    if (res[k] && f.out_[k]->is_jac()) {
      any_jac = true;
    }
  }
  // Evalute Jacobian blocks
  if (any_jac) {
    // Loop over colors
    for (casadi_int c = 0; c < f.coloring_.size2(); ++c) {
      // Loop over input indices for color
      for (casadi_int kc = f.coloring_.colind(c); kc < f.coloring_.colind(c + 1); ++kc) {
        casadi_int vind = f.coloring_.row(kc);
        // Differentiation with respect to what variable
        size_t Jc = f.jac_in_.at(vind);
        // Nominal value
        double nom = self_.variable(Jc).nominal;
        // Set seed for column
        set_seed(mem, Jc, nom);
        // Request corresponding outputs
        for (casadi_int Jk = f.sp_ext_.colind(vind); Jk < f.sp_ext_.colind(vind + 1); ++Jk) {
          casadi_int Jr = f.sp_ext_.row(Jk);
          request(mem, f.jac_out_.at(Jr), Jc);
        }
      }
      // Calculate derivatives
      if (eval_derivative(mem, f)) return 1;
      // Loop over input indices for color
      for (casadi_int kc = f.coloring_.colind(c); kc < f.coloring_.colind(c + 1); ++kc) {
        casadi_int vind = f.coloring_.row(kc);
        // Differentiation with respect to what variable
        size_t Jc = f.jac_in_.at(vind);
        // Inverse of nominal value
        double inv_nom = 1. / self_.variable(Jc).nominal;
        // Fetch Jacobian blocks
        for (size_t k = 0; k < f.out_.size(); ++k) {
          if (res[k] && f.out_[k]->is_jac()) {
            // Find input index
            const std::vector<size_t>& ind2 = f.out_[k]->ind2();
            for (size_t Bc = 0; Bc < ind2.size(); ++Bc) {
              if (ind2[Bc] == Jc) {
                // Column exists in Jacobian block
                const Sparsity& sp = f.sparsity_out(k);
                const std::vector<size_t>& ind1 = f.out_[k]->ind1();
                for (casadi_int Bk = sp.colind(Bc); Bk < sp.colind(Bc + 1); ++Bk) {
                  // Get the Jacobian nonzero
                  double J_nz;
                  get_sens(mem, ind1.at(sp.row(Bk)), &J_nz);
                  // Remove nominal value factor
                  J_nz *= inv_nom;
                  // Save to output
                  res[k][Bk] = J_nz;
                }
              }
            }
          }
        }
      }
    }
  }
  // Successful return
  return 0;
}

FmuInput::~FmuInput() {
}

size_t FmuInput::ind(size_t k) const {
  casadi_error("ind(k) not implemented for " + class_name());
  return -1;
}

const std::vector<size_t>& FmuInput::ind() const {
  casadi_error("ind not implemented for " + class_name());
  static const std::vector<size_t> dummy;
  return dummy;
}

size_t FmuInput::size() const {
  casadi_error("size() not implemented for " + class_name());
  return -1;
}

RegInput::~RegInput() {
}

DummyInput::~DummyInput() {
}

FmuOutput::~FmuOutput() {
}

RegOutput::~RegOutput() {
}

size_t FmuOutput::ind(size_t k) const {
  casadi_error("ind(k) not implemented for " + class_name());
  return -1;
}

const std::vector<size_t>& FmuOutput::ind() const {
  casadi_error("ind() not implemented for " + class_name());
  static const std::vector<size_t> dummy;
  return dummy;
}

const std::vector<size_t>& FmuOutput::ind1() const {
  casadi_error("ind1() not implemented for " + class_name());
  static const std::vector<size_t> dummy;
  return dummy;
}

const std::vector<size_t>& FmuOutput::ind2() const {
  casadi_error("ind2() not implemented for " + class_name());
  static const std::vector<size_t> dummy;
  return dummy;
}

size_t FmuOutput::size() const {
  casadi_error("size() not implemented for " + class_name());
  return -1;
}

JacOutput::~JacOutput() {
}

Sparsity JacOutput::sparsity(const FmuFunction& f) const {
  // DaeBuilder instance
  casadi_assert(f.dae_.alive(), "DaeBuilder instance has been deleted");
  auto dae = static_cast<const DaeBuilderInternal*>(f.dae_->raw_);
  // Get the Jacobian block
  return dae->jac_sparsity(ind1_, ind2_);
}

FmuFunction::FmuFunction(const std::string& name, const DaeBuilder& dae,
    const std::vector<std::string>& name_in,
    const std::vector<std::string>& name_out,
    const std::map<std::string, std::vector<size_t>>& scheme,
    const std::map<std::string, std::vector<size_t>>& lc)
    : FunctionInternal(name), dae_(dae), scheme_(scheme), lc_(lc) {
  // Get input IDs
  in_.resize(name_in.size(), nullptr);
  for (size_t k = 0; k < name_in.size(); ++k) {
    // Look for prefix
    if (has_prefix(name_in[k])) {
      // Get the prefix
      std::string pref, rem;
      pref = pop_prefix(name_in[k], &rem);
      if (pref == "out") {
        // Nondifferentiated function output (unused)
        in_[k] = new DummyInput(scheme.at(rem).size());
      } else {
        // No such prefix
        casadi_error("No such prefix: " + pref);
      }
    } else {
      // No prefix - regular input
      in_[k] = new RegInput(scheme.at(name_in[k]));
    }
  }
  // Get input IDs
  out_.resize(name_out.size(), nullptr);
  for (size_t k = 0; k < name_out.size(); ++k) {
    // Look for prefix
    if (has_prefix(name_out[k])) {
      // Get the prefix
      std::string part1, rem;
      part1 = pop_prefix(name_out[k], &rem);
      if (part1 == "jac") {
        // Jacobian block
        casadi_assert(has_prefix(rem), "Two arguments expected for Jacobian block");
        std::string part2 = pop_prefix(rem, &rem);
        out_[k] = new JacOutput(scheme.at(part2), scheme.at(rem));
      } else {
        // No such prefix
        casadi_error("No such prefix: " + part1);
      }
    } else {
      // No prefix - regular output
      out_[k] = new RegOutput(scheme.at(name_out[k]));
    }
  }
  // Set input/output names
  name_in_ = name_in;
  name_out_ = name_out;
  // Default options
  enable_ad_ = dae->provides_directional_derivative_;
  validate_ad_ = false;
  step_ = 1e-6;
  abstol_ = 1e-3;
  reltol_ = 1e-3;
  u_aim_ = 100.;
  h_iter_ = 0;
  h_min_ = 0;
  h_max_ = inf;
  // No memory object allocated
  mem_ = -1;
}

FmuFunction::~FmuFunction() {
  // Free memory
  clear_mem();
  // Free input and output memory structures
  for (FmuInput* e : in_) if (e) delete(e);
  for (FmuOutput* e : out_) if (e) delete(e);
  // Release memory if DaeBuilder instance still exists
  if (mem_ >= 0 && dae_.alive()) {
    // Get a pointer to the DaeBuilder class
    auto dae = static_cast<const DaeBuilderInternal*>(dae_->raw_);
    // Release memory object
    dae->fmu_->release(mem_);
  }
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
  fd_ = to_enum<Fmu::FdMode>(fd_method_, "forward");

  // Get a pointer to the DaeBuilder class
  casadi_assert(dae_.alive(), "DaeBuilder instance has been deleted");
  auto dae = static_cast<const DaeBuilderInternal*>(dae_->raw_);

  // Consistency checks
  if (enable_ad_) casadi_assert(dae->provides_directional_derivative_,
    "FMU does not provide support for analytic derivatives");
  if (validate_ad_ && !enable_ad_) casadi_error("Inconsistent options");

  // Collect all inputs in any Jacobian block
  std::vector<bool> in_jac(dae->variables_.size(), false);
  for (FmuOutput* i : out_) {
    if (i->is_jac()) {
      for (size_t j : i->ind2()) in_jac[j] = true;
    }
  }
  jac_in_.clear();
  for (size_t k = 0; k < in_jac.size(); ++k) {
    if (in_jac[k]) jac_in_.push_back(k);
  }
  // Collect all outputs in any Jacobian block
  std::fill(in_jac.begin(), in_jac.end(), false);
  for (FmuOutput* i : out_) {
    if (i->is_jac()) {
      for (size_t j : i->ind1()) in_jac[j] = true;
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

  // Load on first encounter
  if (dae->fmu_ == 0) dae->init_fmu();

  // Create instance
  mem_ = dae->fmu_->checkout();
}

int FmuFunction::eval(const double** arg, double** res, casadi_int* iw, double* w,
    void* mem) const {
  // DaeBuilder instance
  casadi_assert(dae_.alive(), "DaeBuilder instance has been deleted");
  auto dae = static_cast<const DaeBuilderInternal*>(dae_->raw_);
  // Evaluate fmu
  int flag = dae->fmu_->eval(mem_, arg, res, *this);
  // Return error flag
  return flag;
}

casadi_int FmuFunction::n_pert() const {
  switch (fd_) {
  case Fmu::FORWARD:
  case Fmu::BACKWARD:
    return 1;
  case Fmu::CENTRAL:
    return 2;
  case Fmu::SMOOTHING:
    return 4;
  default: break;
  }
  casadi_error("Not implemented");
}

std::string to_string(Fmu::FdMode v) {
  switch (v) {
  case Fmu::FORWARD: return "forward";
  case Fmu::BACKWARD: return "backward";
  case Fmu::CENTRAL: return "central";
  case Fmu::SMOOTHING: return "smoothing";
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
  // DaeBuilder instance
  casadi_assert(dae_.alive(), "DaeBuilder instance has been deleted");
  DaeBuilder dae = shared_cast<DaeBuilder>(const_cast<FmuFunction*>(this)->dae_.shared());
  // Return value
  Function ret;
  ret.own(new FmuFunction(name, dae, inames, onames, scheme_, lc_));
  // Hack: Manually enable finite differenting (pending implementation in class)
  Dict opts2 = opts;
  opts2["enable_fd"] = true;
  ret->construct(opts2);
  return ret;
}

bool FmuFunction::has_jac_sparsity(casadi_int oind, casadi_int iind) const {
  return in_.at(iind)->is_reg() && out_.at(oind)->is_reg();
}

Sparsity FmuFunction::get_jac_sparsity(casadi_int oind, casadi_int iind,
    bool symmetric) const {
  // DaeBuilder instance
  casadi_assert(dae_.alive(), "DaeBuilder instance has been deleted");
  auto dae = static_cast<const DaeBuilderInternal*>(dae_->raw_);
  // Get the Jacobian block
  return dae->jac_sparsity(out_.at(oind)->ind(), in_.at(iind)->ind());
}

#endif  // WITH_FMU

} // namespace casadi

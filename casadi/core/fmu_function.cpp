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
  // Create instance if one is not already allocated
  if (m.c == 0) m.c = instantiate();
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
    mem_.at(mem).in_use = false;
  }
}

int Fmu::setup_experiment(int mem, const FmuFunction& f) {
  fmi2Status status = setup_experiment_(memory(mem), f.fmutol_ > 0, f.fmutol_, 0., fmi2True, 1.);
  if (status != fmi2OK) {
    casadi_warning("fmi2SetupExperiment failed");
    return 1;
  }
  return 0;
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

void Fmu::request(int mem, size_t id) {
  // Get memory
  Memory& m = mem_.at(mem);
  // Mark as requested
  m.requested_.at(id) = true;
}

int Fmu::eval(int mem) {
  // Gather inputs and outputs
  gather_io(mem);
  // Get memory
  Memory& m = mem_.at(mem);
  // Number of inputs and outputs
  size_t n_set = m.id_in_.size();
  size_t n_out = m.id_out_.size();
  // Set all variables
  fmi2Status status = set_real_(m.c, get_ptr(m.vr_in_), n_set, get_ptr(m.v_in_));
  if (status != fmi2OK) {
    casadi_warning("fmi2SetReal failed");
    return 1;
  }
  // Initialization mode begins
  if (enter_initialization_mode(mem)) return 1;
  // Initialization mode ends
  if (exit_initialization_mode(mem)) return 1;
  // Set all variables again (state variables get reset during initialization?)
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
  // Collect input indices
  m.id_in_.clear();
  for (size_t id = 0; id < m.changed_.size(); ++id) {
    if (m.changed_[id]) {
      m.id_in_.push_back(id);
      m.changed_[id] = false;
    }
  }
  // Get corresponding value references
  m.vr_in_.clear();
  for (size_t id : m.id_in_) {
    m.vr_in_.push_back(self_.variable(id).value_reference);
  }
  // Collect output indices
  m.id_out_.clear();
  for (size_t id = 0; id < m.requested_.size(); ++id) {
    if (m.requested_[id]) {
      m.id_out_.push_back(id);
      m.requested_[id] = false;
    }
  }
  // Get corresponding value references
  m.vr_out_.clear();
  for (size_t id : m.id_out_) {
    m.vr_out_.push_back(self_.variable(id).value_reference);
  }
  // Get function inputs
  m.v_in_.clear();
  for (size_t id : m.id_in_) m.v_in_.push_back(m.buffer_[id]);
}

int Fmu::eval_derivative(int mem, const FmuFunction& f) {
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
  // Quick return if nothing to be calculated
  if (n_unknown == 0) return 0;
  // Allocate result
  m.d_out_.resize(n_unknown);
  // Calculate derivatives using FMU directional derivative support
  if (f.enable_ad_) {
    fmi2Status status = get_directional_derivative_(m.c, get_ptr(m.vr_out_), n_unknown,
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
  }
  // Calculate derivatives using finite differences
  if (!f.enable_ad_ || f.validate_ad_) {
    // Average nominal value
    double nom = 0;
    for (size_t id : m.id_in_) {
      nom += double(self_.variable(id).nominal);
    }
    nom /= n_known;
    // Step size (fixed for now)
    double h = nom * f.step_;
    // For backward, negate h
    if (f.fd_ == BACKWARD) h = -h;
    // Number of perturbations
    casadi_int n_pert = f.fd_ == FORWARD || f.fd_ == BACKWARD ? 1 :
      f.fd_ == CENTRAL ? 2 : 4;
    // Allocate memory for perturbed outputs
    m.fd_out_.resize(n_pert * n_unknown);
    // Get unperturbed outputs
    m.v_out_.resize(n_unknown);
    fmi2Status status = get_real_(m.c, get_ptr(m.vr_out_), n_unknown, get_ptr(m.v_out_));
    if (status != fmi2OK) {
      casadi_warning("fmi2GetReal failed");
      return 1;
    }
    // Calculate all perturbed outputs
    double* yk[4];
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
      }
      // Pass perturbed inputs to FMU
      casadi_axpy(n_known, pert, get_ptr(m.d_in_), get_ptr(m.v_in_));
      status = set_real_(m.c, get_ptr(m.vr_in_), n_known, get_ptr(m.v_in_));
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
      // Remove purburbation
      casadi_axpy(n_known, -pert, get_ptr(m.d_in_), get_ptr(m.v_in_));
    }
    // Restore FMU inputs
    status = set_real_(m.c, get_ptr(m.vr_in_), n_known, get_ptr(m.v_in_));
    if (status != fmi2OK) {
      casadi_warning("fmi2SetReal failed");
      return 1;
    }
    // FD memory
    casadi_finite_diff_mem<double> fd_mem;
    fd_mem.reltol = nan;
    fd_mem.abstol = nan;
    fd_mem.smoothing = eps;
    switch (f.fd_) {
      case FORWARD:
      case BACKWARD:
        (void)casadi_forward_diff(yk, get_ptr(m.v_out_),
          get_ptr(m.d_out_), h, n_unknown, &fd_mem);
        break;
      case CENTRAL:
        (void)casadi_central_diff(yk, get_ptr(m.v_out_), get_ptr(m.d_out_),
          h, n_unknown, &fd_mem);
        break;
      case SMOOTHING:
        (void)casadi_smoothing_diff(yk, get_ptr(m.v_out_), get_ptr(m.d_out_),
          h, n_unknown, &fd_mem);
        break;
    }
    // Collect requested variables
    auto it = m.d_out_.begin();
    for (size_t id : m.id_out_) {
      // Get the value
      double d = *it++;
      // Use FD instead of AD or to compare with AD
      if (f.validate_ad_) {
        // Access variable
        const Variable& v = self_.variable(id);
        // Get the nominal value
        double nom_out = double(v.nominal);
        // Maximum error
        double etol = std::fabs(d) * f.reltol_ + nom_out * f.abstol_;
        // Check relative error
        if (std::fabs(m.sens_[id] - d) > etol) {
          // Issue warning
          std::stringstream ss;
          ss << "Inconsistent derivatives of " << v.name << " w.r.t. ";
          for (size_t j = 0; j < n_known; ++j) {
            if (j != 0) ss << ", ";
            ss << self_.variable(m.id_in_[j]).name;
          }
          ss << ". Got " << m.sens_[id] << " for AD vs. " << d << " for FD" << std::endl;
          casadi_warning(ss.str());
        }
      } else {
        // Use instead of AD
        m.sens_[id] = d;
      }
    }
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
  // Set inputs
  for (size_t k = 0; k < f.id_in_.size(); ++k) {
    for (size_t i = 0; i < f.id_in_[k].size(); ++i) {
      set(mem, f.id_in_[k][i], arg[k] ? arg[k][i] : 0);
    }
  }
  // Request outputs to be evaluated
  for (size_t k = 0; k < f.id_out_.size(); ++k) {
    if (res[k]) {
      for (size_t i = 0; i < f.id_out_[k].size(); ++i) {
        request(mem, f.id_out_[k][i]);
      }
    }
  }
  // Reset solver
  if (setup_experiment(mem, f)) return 1;
  // Evaluate
  if (eval(mem)) return 1;
  // Reset solver
  if (reset(mem)) return 1;
  // Get outputs
  for (size_t k = 0; k < f.id_out_.size(); ++k) {
    if (res[k]) {
      for (size_t i = 0; i < f.id_out_[k].size(); ++i) {
        get(mem, f.id_out_[k][i], &res[k][i]);
      }
    }
  }
  // Successful return
  return 0;
}

int Fmu::eval_jac(int mem, const double** arg, double** res, const FmuFunction& f,
    const std::vector<Sparsity>& sp_jac) {
  // Set inputs
  for (size_t k = 0; k < f.id_in_.size(); ++k) {
    for (size_t i = 0; i < f.id_in_[k].size(); ++i) {
      set(mem, f.id_in_[k][i], arg[k] ? arg[k][i] : 0);
    }
  }
  // Reset solver
  if (setup_experiment(mem, f)) return 1;
  // Evaluate
  if (eval(mem)) return 1;
  // Loop over function inputs
  for (size_t i1 = 0; i1 < f.id_in_.size(); ++i1) {
    // Calculate Jacobian, one column at a time
    for (casadi_int i2 = 0; i2 < f.id_in_[i1].size(); ++i2) {
      // Set seed for column
      set_seed(mem, f.id_in_[i1][i2], 1.);
      // Loop over function output
      for (size_t j1 = 0; j1 < f.id_out_.size(); ++j1) {
        // Index of the Jacobian block
        size_t res_ind = j1 * f.id_in_.size() + i1;
        // Only calculate outputs that are requested
        if (res[res_ind]) {
          // Get sparsity
          const casadi_int* colind = sp_jac[res_ind].colind();
          const casadi_int* row = sp_jac[res_ind].row();
          // Request all nonzero elements of the column
          for (size_t k = colind[i2]; k < colind[i2 + 1]; ++k) {
            request(mem, f.id_out_[j1][row[k]]);
          }
        }
      }
      // Calculate derivatives
      if (eval_derivative(mem, f)) return 1;
      // Loop over function outputs
      for (size_t j1 = 0; j1 < f.id_out_.size(); ++j1) {
        // Index of the Jacobian block
        size_t res_ind = j1 * f.id_in_.size() + i1;
        // Only calculate outputs that are requested
        if (res[res_ind]) {
          // Get sparsity
          const casadi_int* colind = sp_jac[res_ind].colind();
          const casadi_int* row = sp_jac[res_ind].row();
          // Collect all nonzero elements of the column
          for (size_t k = colind[i2]; k < colind[i2 + 1]; ++k) {
            get_sens(mem, f.id_out_[j1][row[k]], &res[res_ind][k]);
          }
        }
      }
    }
  }
  // Reset solver
  if (reset(mem)) return 1;
  // Successful return
  return 0;
}

int Fmu::eval_adj(int mem, const double** arg, double** res, const FmuFunction& f) {
  // Set inputs
  for (size_t k = 0; k < f.id_in_.size(); ++k) {
    for (size_t i = 0; i < f.id_in_[k].size(); ++i) {
      set(mem, f.id_in_[k][i], arg[k] ? arg[k][i] : 0);
    }
  }
  // Setup experiment
  if (setup_experiment(mem, f)) return 1;
  // Evaluate
  if (eval(mem)) return 1;
  // Loop over function inputs
  for (size_t i1 = 0; i1 < f.id_in_.size(); ++i1) {
    // Sensitivities to be calculated
    double* sens = res[i1];
    // Skip if not requested
    if (sens == 0) continue;
    // Calculate Jacobian, one column at a time
    for (casadi_int i2 = 0; i2 < f.id_in_[i1].size(); ++i2) {
      // Initialize return to zero
      sens[i2] = 0;
      // Set seed for column i
      set_seed(mem, f.id_in_[i1][i2], 1.);
      // Request all elements of the column, unless corresponding seed is zero
      for (size_t j1 = 0; j1 < f.id_out_.size(); ++j1) {
        if (arg[f.id_in_.size() + f.id_out_.size() + j1]) {
          for (size_t j2 = 0; j2 < f.id_out_[j1].size(); ++j2) {
            request(mem, f.id_out_[j1][j2]);
          }
        }
      }
      // Calculate derivatives
      if (eval_derivative(mem, f)) return 1;
      // Get sensitivities
      for (size_t j1 = 0; j1 < f.id_out_.size(); ++j1) {
        const double* seed = arg[f.id_in_.size() + f.id_out_.size() + j1];
        if (seed) {
          for (size_t j2 = 0; j2 < f.id_out_[j1].size(); ++j2) {
            double J_ij;
            get_sens(mem, f.id_out_[j1][j2], &J_ij);
            sens[i2] += seed[j2] * J_ij;
          }
        }
      }
    }
  }
  // Reset solver
  if (reset(mem)) return 1;
  // Successful return
  return 0;
}

FmuFunction::FmuFunction(const std::string& name, const DaeBuilder& dae,
    const std::vector<std::vector<size_t>>& id_in,
    const std::vector<std::vector<size_t>>& id_out,
    const std::vector<std::string>& name_in,
    const std::vector<std::string>& name_out)
  : FunctionInternal(name), dae_(dae), id_in_(id_in), id_out_(id_out) {
  // Names of inputs
  if (!name_in.empty()) {
    casadi_assert(id_in.size() == name_in.size(), "Mismatching number of input names");
    name_in_ = name_in;
  }
  // Names of outputs
  if (!name_out.empty()) {
    casadi_assert(id_out.size() == name_out.size(), "Mismatching number of output names");
    name_out_ = name_out;
  }

  // Default options
  enable_ad_ = dae->provides_directional_derivative_;
  validate_ad_ = false;
  step_ = 1e-6;
  abstol_ = 1e-3;
  reltol_ = 1e-3;
  fmutol_ = 0;
}

FmuFunction::~FmuFunction() {
  // Free memory
  clear_mem();
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
    {"fmutol",
     {OT_DOUBLE,
      "Tolerance to be passed to the fmu (0 if not defined)"}}
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
    } else if (op.first=="fmutol") {
      fmutol_ = op.second;
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

  // Load on first encounter
  if (dae->fmu_ == 0) dae->init_fmu();
}

Sparsity FmuFunction::get_sparsity_in(casadi_int i) {
  return Sparsity::dense(id_in_.at(i).size(), 1);
}

Sparsity FmuFunction::get_sparsity_out(casadi_int i) {
  return Sparsity::dense(id_out_.at(i).size(), 1);
}

Sparsity FmuFunction::get_jac_sparsity(casadi_int oind, casadi_int iind,
    bool symmetric) const {
  // DaeBuilder instance
  casadi_assert(dae_.alive(), "DaeBuilder instance has been deleted");
  auto dae = static_cast<const DaeBuilderInternal*>(dae_->raw_);
  // Lookup for inputs
  std::vector<casadi_int> lookup(dae->variables_.size(), -1);
  for (casadi_int i = 0; i < id_in_.at(iind).size(); ++i)
    lookup.at(id_in_.at(iind).at(i)) = i;
  // Nonzeros of the Jacobian
  std::vector<casadi_int> row, col;
  // Loop over output nonzeros
  for (casadi_int j = 0; j < id_out_.at(oind).size(); ++j) {
    // Loop over dependencies
    for (casadi_int d : dae->variables_.at(id_out_.at(oind).at(j)).dependencies) {
      casadi_int i = lookup.at(d);
      if (i >= 0) {
        row.push_back(j);
        col.push_back(i);
      }
    }
  }
  // Assemble sparsity pattern
  return Sparsity::triplet(id_out_.at(oind).size(), id_in_.at(iind).size(), row, col);
}

int FmuFunction::eval(const double** arg, double** res, casadi_int* iw, double* w,
    void* mem) const {
  // DaeBuilder instance
  casadi_assert(dae_.alive(), "DaeBuilder instance has been deleted");
  auto dae = static_cast<const DaeBuilderInternal*>(dae_->raw_);
  // Create instance
  int m = dae->fmu_->checkout();
  // Evaluate fmu
  int flag = dae->fmu_->eval(m, arg, res, *this);
  // Release memory object
  dae->fmu_->release(m);
  // Return error flag
  return flag;
}

Function FmuFunction::get_jacobian(const std::string& name, const std::vector<std::string>& inames,
    const std::vector<std::string>& onames, const Dict& opts) const {
  Function ret;
  ret.own(new FmuFunctionJac(name));
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
  casadi_assert(nadj == 1, "Not supported");
  Function ret;
  ret.own(new FmuFunctionAdj(name));
  // Hack: Manually enable finite differenting (pending implementation in class)
  Dict opts2 = opts;
  opts2["enable_fd"] = true;
  ret->construct(opts2);
  return ret;
}

FmuFunctionJac::~FmuFunctionJac() {
  // Free memory
  clear_mem();
}

int FmuFunctionJac::eval(const double** arg, double** res, casadi_int* iw, double* w,
    void* mem) const {
  // Non-differentiated class
  auto self = derivative_of_.get<FmuFunction>();
  // DaeBuilder instance
  casadi_assert(self->dae_.alive(), "DaeBuilder instance has been deleted");
  auto dae = static_cast<const DaeBuilderInternal*>(self->dae_->raw_);
  // Create instance
  int m = dae->fmu_->checkout();
  // Evaluate fmu
  int flag = dae->fmu_->eval_jac(m, arg, res, *self, sparsity_out_);
  // Release memory object
  dae->fmu_->release(m);
  // Return error flag
  return flag;
}

FmuFunctionAdj::~FmuFunctionAdj() {
  // Free memory
  clear_mem();
}

int FmuFunctionAdj::eval(const double** arg, double** res, casadi_int* iw, double* w,
    void* mem) const {
  // Non-differentiated class
  auto self = derivative_of_.get<FmuFunction>();
  // DaeBuilder instance
  casadi_assert(self->dae_.alive(), "DaeBuilder instance has been deleted");
  auto dae = static_cast<const DaeBuilderInternal*>(self->dae_->raw_);
  // Create instance
  int m = dae->fmu_->checkout();
  // Evaluate fmu
  int flag = dae->fmu_->eval_adj(m, arg, res, *self);
  // Release memory object
  dae->fmu_->release(m);
  // Return error flag
  return flag;
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

#endif  // WITH_FMU

} // namespace casadi

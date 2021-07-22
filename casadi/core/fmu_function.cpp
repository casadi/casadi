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

#include <fstream>
#include <iostream>
#include <sstream>

namespace casadi {

#ifdef WITH_FMU

FmuFunction::FmuFunction(const DaeBuilderInternal& dae,
    const std::string& name,
    const std::vector<std::vector<casadi_int>>& id_in,
    const std::vector<std::vector<casadi_int>>& id_out,
    const std::vector<std::string>& name_in,
    const std::vector<std::string>& name_out)
  : FunctionInternal(name) {
  // Cast id_in to right type
  id_in_.resize(id_in.size());
  for (size_t k = 0; k < id_in.size(); ++k) {
    id_in_[k].resize(id_in[k].size());
    std::copy(id_in[k].begin(), id_in[k].end(), id_in_[k].begin());
  }
  // Cast id_out to right type
  id_out_.resize(id_out.size());
  for (size_t k = 0; k < id_out.size(); ++k) {
    id_out_[k].resize(id_out[k].size());
    std::copy(id_out[k].begin(), id_out[k].end(), id_out_[k].begin());
  }
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
  // Information from DaeBuilder instance
  path_ = dae.path_;
  guid_ = dae.guid_;
  instance_name_ = dae.model_identifier_;
  provides_directional_derivative_ = dae.provides_directional_derivative_;
  // Initialize to null pointers
  instantiate_ = 0;
  free_instance_ = 0;
  setup_experiment_ = 0;
  enter_initialization_mode_ = 0;
  exit_initialization_mode_ = 0;
  enter_continuous_time_mode_ = 0;
  set_real_ = 0;
  set_boolean_ = 0;
  get_real_ = 0;
  get_directional_derivative_ = 0;
  c_ = 0;
}

FmuFunction::~FmuFunction() {
  // Free memory
  if (c_ && free_instance_) free_instance_(c_);
  clear_mem();
}

void FmuFunction::init(const Dict& opts) {
  // Call the initialization method of the base class
  FunctionInternal::init(opts);
  // Directory where the DLL is stored, per the FMI specification
  std::string instance_name_no_dot = instance_name_;
  std::replace(instance_name_no_dot.begin(), instance_name_no_dot.end(), '.', '_');
  std::string dll_path = path_ + "/binaries/" + system_infix()
    + "/" + instance_name_no_dot + dll_suffix();

  // Path to resource directory
  resource_loc_ = "file://" + path_ + "/resources";

  // Load the DLL
  li_ = Importer(dll_path, "dll");

  // Load functions
  instantiate_ = reinterpret_cast<fmi2InstantiateTYPE*>(get_function("fmi2Instantiate"));
  free_instance_ = reinterpret_cast<fmi2FreeInstanceTYPE*>(get_function("fmi2FreeInstance"));
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
  if (provides_directional_derivative_) {
    get_directional_derivative_ = reinterpret_cast<fmi2GetDirectionalDerivativeTYPE*>(
      get_function("fmi2GetDirectionalDerivative"));
  }

  // Callback functions
  functions_.logger = logger;
  functions_.allocateMemory = calloc;
  functions_.freeMemory = free;
  functions_.stepFinished = 0;
  functions_.componentEnvironment = 0;

  // Create instance
  fmi2String instanceName = instance_name_.c_str();
  fmi2Type fmuType = fmi2ModelExchange;
  fmi2String fmuGUID = guid_.c_str();
  fmi2String fmuResourceLocation = resource_loc_.c_str();
  fmi2Boolean visible = fmi2False;
  fmi2Boolean loggingOn = fmi2False;
  c_ = instantiate_(instanceName, fmuType, fmuGUID, fmuResourceLocation,
    &functions_, visible, loggingOn);
  if (c_ == 0) casadi_error("fmi2Instantiate failed");

  // Reset solver
  fmi2Status status = setup_experiment_(c_, fmi2False, 0.0, 0., fmi2True, 1.);
  if (status != fmi2OK) casadi_error("fmi2SetupExperiment failed");

  // Initialization mode begins
  status = enter_initialization_mode_(c_);
  if (status != fmi2OK) casadi_error("fmi2EnterInitializationMode failed: " + str(status));

  // Initialization mode ends
  status = exit_initialization_mode_(c_);
  if (status != fmi2OK) casadi_error("fmi2ExitInitializationMode failed");

  // Continuous time mode only for event mode
  // status = enter_continuous_time_mode_(c_);
  // if (status != fmi2OK) casadi_error("fmi2EnterContinuousTimeMode failed: " + str(status));
  // uout() << "enter_continuous_time_mode_ done" << std::endl;
}

std::string FmuFunction::system_infix() {
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

std::string FmuFunction::dll_suffix() {
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

signal_t FmuFunction::get_function(const std::string& symname) {
  // Load the function
  signal_t f = li_.get_function(symname);
  // Ensure that it was found
  casadi_assert(f != 0, "Cannot retrieve '" + symname + "'");
  // Return function to be type converted
  return f;
}

Sparsity FmuFunction::get_sparsity_in(casadi_int i) {
  return Sparsity::dense(id_in_.at(i).size(), 1);
}

Sparsity FmuFunction::get_sparsity_out(casadi_int i) {
  return Sparsity::dense(id_out_.at(i).size(), 1);
}

int FmuFunction::set_inputs(const double** x) const {
  for (size_t k = 0; k < id_in_.size(); ++k) {
    // Set input k
    if (x[k]) {
      // Input is available
      fmi2Status status = set_real_(c_, get_ptr(id_in_[k]), id_in_[k].size(), x[k]);
      if (status != fmi2OK) {
        casadi_warning("fmi2SetReal failed for input '" + name_in_[k] + "'");
        return 1;
      }
    } else {
      // Input is null (zero)
      static const double zero = 0;
      for (size_t i = 0; i < id_in_[k].size(); ++i) {
        fmi2Status status = set_real_(c_, &id_in_[k][i], 1, &zero);
        if (status != fmi2OK) {
          casadi_warning("fmi2SetReal failed for input '" + name_in_[k] + "', component " + str(i)
            + " (value reference " + str(id_in_[k][i]) + ")" );
          return 1;
        }
      }
    }
  }
  // Success
  return 0;
}

int FmuFunction::eval(const double** arg, double** res, casadi_int* iw, double* w,
    void* mem) const {
  // Return flag
  fmi2Status status;

  // Initialization mode begins
  // status = enter_initialization_mode_(c_);
  // if (status != fmi2OK) casadi_error("fmi2EnterInitializationMode failed: " + str(status));

  // Pass inputs
  if (set_inputs(arg)) return 1;

  // Evaluate
  if (res[0]) {
    status = get_real_(c_, get_ptr(id_out_[0]), id_out_[0].size(), res[0]);
    if (status != fmi2OK) casadi_error("fmi2GetReal failed for output '" + name_out_[0] + "'");
  }

  // Initialization mode ends
  // status = exit_initialization_mode_(c_);
  // if (status != fmi2OK) casadi_error("fmi2ExitInitializationMode failed");

  return 0;
}

int FmuFunction::eval_jac(const double** arg, double** res, casadi_int* iw, double* w,
    void* mem) const {
  // Dimensions
  casadi_int n_xd = nnz_in(0);
  casadi_int n_yd = nnz_out(0);
  // Outputs
  double* jac = res[0];
  // Forward seed, sensitivity
  double* fwd_xd = w; w += n_xd;
  double* fwd_yd = w; w += n_yd;
  // FMI return flag
  fmi2Status status;
  // Pass inputs
  if (set_inputs(arg)) return 1;
  // Clear seeds
  casadi_clear(fwd_xd, n_xd);
  // Calculate Jacobian, one column at a time
  for (casadi_int i = 0; i < n_xd; ++i) {
    // Set seed for column i
    fwd_xd[i] = 1.;
    // Calculate directional derivative
    status = get_directional_derivative_(c_, get_ptr(id_out_[0]), id_out_[0].size(),
      get_ptr(id_in_[0]), id_in_[0].size(), fwd_xd, fwd_yd);
    if (status != fmi2OK) {
      casadi_warning("fmi2GetDirectionalDerivative failed");
      return 1;
    }
    // Copy column to Jacobian
    casadi_copy(fwd_yd, n_yd, jac);
    jac += n_yd;
    // Remove seed
    fwd_xd[i] = 0;
  }
  // Successful return
  return 0;
}

int FmuFunction::eval_adj(const double** arg, double** res, casadi_int* iw, double* w,
    void* mem) const {
  // Dimensions
  casadi_int n_xd = nnz_in(0);
  casadi_int n_yd = nnz_out(0);
  // Inputs
  // const double* xd = arg[0];
  // const double* xn = arg[1];  // not implemented
  // const double* out_yd = arg[2];  // not used
  // const double* out_yn = arg[3];  // not used
  const double* adj_yd = arg[4];
  // const double* adj_yn = arg[5];  // non-differentiable
  // Outputs
  double* adj_xd = res[0];
  // double* adj_xn = res[1];  // non-differentiable, not implemented
  // Forward seed, sensitivity for calculating columns of the Jacobian
  double* fwd_xd = w; w += n_xd;
  double* fwd_yd = w; w += n_yd;
  // FMI return flag
  fmi2Status status;
  // Pass inputs
  if (set_inputs(arg)) return 1;
  // Reset results
  casadi_clear(adj_xd, n_xd);
  // Clear seeds
  casadi_clear(fwd_xd, n_xd);
  // Calculate Jacobian, one column at a time
  for (casadi_int i = 0; i < n_xd; ++i) {
    // Set seed for column i
    fwd_xd[i] = 1.;
    // Calculate directional derivative
    status = get_directional_derivative_(c_, get_ptr(id_out_[0]), id_out_[0].size(),
      get_ptr(id_in_[0]), id_in_[0].size(), fwd_xd, fwd_yd);
    if (status != fmi2OK) {
      casadi_warning("fmi2GetDirectionalDerivative failed");
      return 1;
    }
    // Add contribution from first seed
    adj_xd[i] += casadi_dot(n_yd, fwd_yd, adj_yd);
    // Remove seed
    fwd_xd[i] = 0;
  }
  // Successful return
  return 0;
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

void FmuFunction::logger(fmi2ComponentEnvironment componentEnvironment,
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

FmuFunctionJac::~FmuFunctionJac() {
  clear_mem();
}

void FmuFunctionJac::init(const Dict& opts) {
  // Call the base class initializer
  FunctionInternal::init(opts);
  // Work vectors
  alloc_w(derivative_of_.nnz_in(0), true);
  alloc_w(derivative_of_.nnz_out(0), true);
}

int FmuFunctionJac::eval(const double** arg, double** res, casadi_int* iw, double* w,
    void* mem) const {
  // Redirect to non-differentiated class
  auto m = derivative_of_.get<FmuFunction>();
  return m->eval_jac(arg, res, iw, w, mem);
}

FmuFunctionAdj::~FmuFunctionAdj() {
  clear_mem();
}

void FmuFunctionAdj::init(const Dict& opts) {
  // Call the base class initializer
  FunctionInternal::init(opts);
  // Work vectors
  alloc_w(derivative_of_.nnz_in(0), true);
  alloc_w(derivative_of_.nnz_out(0), true);
}

int FmuFunctionAdj::eval(const double** arg, double** res, casadi_int* iw, double* w,
    void* mem) const {
  // Redirect to non-differentiated class
  auto m = derivative_of_.get<FmuFunction>();
  return m->eval_adj(arg, res, iw, w, mem);
}

#endif  // WITH_FMU

} // namespace casadi

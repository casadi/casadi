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


#include "fmu_function_impl.hpp"
#include "casadi_misc.hpp"
#include "serializing_stream.hpp"

#include <fstream>
#include <iostream>
#include <sstream>

namespace casadi {

Function fmu_function(const std::string& name, const std::string& path,
    const std::vector<std::vector<casadi_int>>& id_in,
    const std::vector<std::vector<casadi_int>>& id_out,
    const std::string& guid, const Dict& opts) {
  return Function::create(new FmuFunction(name, path, id_in, id_out, guid), opts);
}

FmuFunction::FmuFunction(const std::string& name, const std::string& path,
    const std::vector<std::vector<casadi_int>>& id_in,
    const std::vector<std::vector<casadi_int>>& id_out,
    const std::string& guid)
  : FunctionInternal(name),
    path_(path), id_in_(id_in), id_out_(id_out), guid_(guid) {
#ifdef WITH_FMU
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
#endif  // WITH_FMU
}

FmuFunction::~FmuFunction() {
#ifdef WITH_FMU
  // Free memory
  if (c_ && free_instance_) free_instance_(c_);
#endif  // WITH_FMU
  clear_mem();
}

void FmuFunction::init(const Dict& opts) {
  // Consistency checks
  casadi_assert(id_in_.size() == 2,
    "Expected two input lists: differentiated and nondifferentiated variables");
  casadi_assert(id_out_.size() == 2,
    "Expected two output lists: differentiated and nondifferentiated variables");

  // Call the initialization method of the base class
  FunctionInternal::init(opts);

#ifdef WITH_FMU
  // Input indices
  xd_.resize(id_in_[0].size());
  std::copy(id_in_[0].begin(), id_in_[0].end(), xd_.begin());
  xn_.resize(id_in_[1].size());
  std::copy(id_in_[1].begin(), id_in_[1].end(), xn_.begin());

  // Output indices
  yd_.resize(id_out_[0].size());
  std::copy(id_out_[0].begin(), id_out_[0].end(), yd_.begin());
  yn_.resize(id_out_[1].size());
  std::copy(id_out_[1].begin(), id_out_[1].end(), yn_.begin());

  // Directory where the DLL is stored, per the FMI specification
  std::string dll_path = path_ + "/binaries/" + system_infix() + "/" + name_ + dll_suffix();

  // Path to resource directory
  resource_loc_ = "file:" + path_ + "/resources/";

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
  get_directional_derivative_ = reinterpret_cast<fmi2GetDirectionalDerivativeTYPE*>(
    get_function("fmi2GetDirectionalDerivative"));

  // Callback functions
  fmi2CallbackFunctions functions;
  functions.componentEnvironment = 0;
  functions.allocateMemory = calloc;
  functions.freeMemory = free;
  functions.stepFinished = 0;
  functions.componentEnvironment = 0;

  // Create instance
  fmi2String instanceName = name_.c_str();
  fmi2Type fmuType = fmi2ModelExchange;
  fmi2String fmuGUID = guid_.c_str();
  fmi2String fmuResourceLocation = resource_loc_.c_str();
  fmi2Boolean visible = fmi2False;
  fmi2Boolean loggingOn = fmi2False;
  c_ = instantiate_(instanceName, fmuType, fmuGUID, fmuResourceLocation,
    &functions, visible, loggingOn);
  if (c_ == 0) casadi_error("fmi2Instantiate failed");

  // Reset solver
  fmi2Status status = setup_experiment_(c_, fmi2False, 0.0, 0., fmi2True, 1.);
  if (status != fmi2OK) casadi_error("fmi2SetupExperiment failed");

#else  // WITH_FMU
  casadi_error("FMU support not enabled. Recompile CasADi with 'WITH_FMU=ON'");
#endif  // WITH_FMU
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

int FmuFunction::eval(const double** arg, double** res, casadi_int* iw, double* w,
    void* mem) const {
#ifdef WITH_FMU
  // Return flag
  fmi2Status status;

  // Initialization mode begins
  status = enter_initialization_mode_(c_);
  if (status != fmi2OK) casadi_error("fmi2EnterInitializationMode failed");

  // Pass differentiable inputs
  status = set_real_(c_, get_ptr(xd_), xd_.size(), arg[0]);
  if (status != fmi2OK) casadi_error("fmi2SetReal failed");

  // Evaluate
  status = get_real_(c_, get_ptr(yd_), yd_.size(), res[0]);
  if (status != fmi2OK) casadi_error("fmi2GetReal failed");

  // Initialization mode ends
  status = exit_initialization_mode_(c_);
  if (status != fmi2OK) casadi_error("fmi2ExitInitializationMode failed");
#endif  // WITH_FMU

  return 0;
}

} // namespace casadi

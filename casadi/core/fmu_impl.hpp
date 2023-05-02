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

#ifndef CASADI_FMU_IMPL_HPP
#define CASADI_FMU_IMPL_HPP

#ifdef WITH_FMU

#include "fmu.hpp"
#include "importer.hpp"
#include "shared_object_internal.hpp"

#include <fmi2Functions.h>

/// \cond INTERNAL

namespace casadi {

// Forward declarations
class DaeBuilderInternal;
struct FmuMemory;
struct InputStruct;

/** \brief Interface to binary FMU

    Internal API.

    \author Joel Andersson
    \date 2023
*/
class CASADI_EXPORT Fmu2 : public SharedObjectInternal {
 public:
  // Constructor
  Fmu2(const std::string& name,
    const std::vector<std::string>& scheme_in, const std::vector<std::string>& scheme_out,
    const std::map<std::string, std::vector<size_t>>& scheme, const std::vector<std::string>& aux);

  /// Destructor
  ~Fmu2() override;

  /** \brief Get type name */
  std::string class_name() const override { return "Fmu2";}

  /** \brief Print */
  void disp(std::ostream& stream, bool more) const override;

  // Initialize
  void init(const DaeBuilderInternal* dae);

  /// Instance name
  std::string name_;

  // DLL
  Importer li_;

  // IO scheme
  std::vector<std::string> scheme_in_, scheme_out_;
  std::map<std::string, std::vector<size_t>> scheme_;

  // Auxilliary outputs
  std::vector<std::string> aux_;

  // Mapping from scheme variable to and from FMU variable indices
  std::vector<size_t> iind_, iind_map_, oind_, oind_map_;

  // Meta information about the input/output variable subsets
  std::vector<double> nominal_in_, nominal_out_;
  std::vector<double> min_in_, min_out_;
  std::vector<double> max_in_, max_out_;
  std::vector<std::string> vn_in_, vn_out_;
  std::vector<fmi2ValueReference> vr_in_, vr_out_;

  // Numerical values for inputs
  std::vector<fmi2Real> value_in_;

  // Reduced space indices for all inputs and outputs
  std::vector<std::vector<size_t>> ired_, ored_;

  // FMU C API function prototypes. Cf. FMI specification 2.0.2
  fmi2InstantiateTYPE* instantiate_;
  fmi2FreeInstanceTYPE* free_instance_;
  fmi2ResetTYPE* reset_;
  fmi2SetupExperimentTYPE* setup_experiment_;
  fmi2EnterInitializationModeTYPE* enter_initialization_mode_;
  fmi2ExitInitializationModeTYPE* exit_initialization_mode_;
  fmi2EnterContinuousTimeModeTYPE* enter_continuous_time_mode_;
  fmi2GetRealTYPE* get_real_;
  fmi2SetRealTYPE* set_real_;
  fmi2GetBooleanTYPE* get_boolean_;
  fmi2SetBooleanTYPE* set_boolean_;
  fmi2GetIntegerTYPE* get_integer_;
  fmi2SetIntegerTYPE* set_integer_;
  fmi2GetStringTYPE* get_string_;
  fmi2SetStringTYPE* set_string_;
  fmi2GetDirectionalDerivativeTYPE* get_directional_derivative_;

  // Callback functions
  fmi2CallbackFunctions functions_;

  // Path to the FMU resource directory
  std::string resource_loc_;

  // Tolerance
  double fmutol_;

  // Instance name
  std::string instance_name_;

  // GUID
  std::string guid_;

  // Logging?
  bool logging_on_;

  // Collection of variable values, all types
  struct Value {
    std::vector<fmi2Real> v_real;
    std::vector<fmi2Integer> v_integer;
    std::vector<fmi2Boolean> v_boolean;
    std::vector<std::string> v_string;
  };

  // Variables used for initialization, by type
  std::vector<fmi2ValueReference> vr_real_, vr_integer_, vr_boolean_, vr_string_;
  std::vector<fmi2Real> init_real_;
  std::vector<fmi2Integer> init_integer_;
  std::vector<fmi2Boolean> init_boolean_;
  std::vector<std::string> init_string_;

  // Auxilliary variables, by type
  std::vector<std::string> vn_aux_real_, vn_aux_integer_, vn_aux_boolean_, vn_aux_string_;
  std::vector<fmi2ValueReference> vr_aux_real_, vr_aux_integer_, vr_aux_boolean_, vr_aux_string_;
  Value aux_value_;

  // Sparsity pattern for extended Jacobian, Hessian
  Sparsity jac_sp_, hess_sp_;

  /** \brief Get the number of scheme inputs */
  size_t n_in() const { return iind_.size();}

  /** \brief Get the number of scheme outputs */
  size_t n_out() const { return oind_.size();}

  // Index lookup for input
  size_t index_in(const std::string& n) const;

  // Index lookup for output
  size_t index_out(const std::string& n) const;

  // Does the FMU support analytic derivatives?
  bool has_ad() const { return get_directional_derivative_ != nullptr; }

  // Load an FMI function
  template<typename T>
  T* load_function(const std::string& symname);

  // New memory object
  fmi2Component instantiate() const;

  // Free FMU instance
  void free_instance(void* c) const;

  // Reset solver
  int reset(fmi2Component c);

  // Setup experiment
  void setup_experiment(fmi2Component c) const;

  // Enter initialization mode
  int enter_initialization_mode(fmi2Component c) const;

  // Exit initialization mode
  int exit_initialization_mode(fmi2Component c) const;

  // Copy values set in DaeBuilder to FMU
  int set_values(fmi2Component c) const;

  // Retrieve input variable values from FMU
  int get_in(fmi2Component c, std::vector<fmi2Real>* v) const;

  // Retrieve auxilliary variables from FMU
  int get_aux(fmi2Component c, Value* v) const;

  /** \brief Get stats */
  void get_stats(FmuMemory* m, Dict* stats,
    const std::vector<std::string>& name_in, const InputStruct* in) const;

  /** \brief Initalize memory block */
  int init_mem(FmuMemory* m) const;

  // Set value
  void set(FmuMemory* m, size_t ind, const double* value) const;

  // Request the calculation of a variable
  void request(FmuMemory* m, size_t ind) const;

  // Gather user inputs and outputs
  void gather_io(FmuMemory* m) const;

  // Calculate all requested variables
  int eval(FmuMemory* m) const;

  // Get a calculated variable
  void get(FmuMemory* m, size_t id, double* value) const;

  // Set seed
  void set_seed(FmuMemory* m, casadi_int nseed, const casadi_int* id, const double* v) const;

  // Request the calculation of a sensitivity
  void request_sens(FmuMemory* m, casadi_int nsens, const casadi_int* id,
    const casadi_int* wrt_id) const;

  // Get calculated derivatives
  void get_sens(FmuMemory* m, casadi_int nsens, const casadi_int* id, double* v) const;

  // Gather user sensitivities
  void gather_sens(FmuMemory* m) const;

  // Calculate directional derivatives using AD
  int eval_ad(FmuMemory* m) const;

  // Calculate directional derivatives using FD
  int eval_fd(FmuMemory* m, bool independent_seeds) const;

  // Calculate directional derivatives
  int eval_derivative(FmuMemory* m, bool independent_seeds) const;

  // Get Jacobian sparsity for a subset of inputs and outputs
  Sparsity jac_sparsity(const std::vector<size_t>& osub, const std::vector<size_t>& isub) const;

  // Get Hessian sparsity for a subset of inputs
  Sparsity hess_sparsity(const std::vector<size_t>& r, const std::vector<size_t>& c) const;

  /// @{
  /** \brief Retreive nominal values */
  std::vector<double> all_nominal_in(size_t i) const;
  std::vector<double> all_nominal_out(size_t i) const;
  /// @}

  // Print description of an input
  std::string desc_in(FmuMemory* m, size_t id, bool more = true) const;

  // Name of system, per the FMI specification
  static std::string system_infix();

  // DLL suffix, per the FMI specification
  static std::string dll_suffix();

  // Process message
  static void logger(fmi2ComponentEnvironment componentEnvironment,
    fmi2String instanceName,
    fmi2Status status,
    fmi2String category,
    fmi2String message, ...);
};

} // namespace casadi

/// \endcond

#endif  // WITH_FMU

#endif // CASADI_FMU_IMPL_HPP

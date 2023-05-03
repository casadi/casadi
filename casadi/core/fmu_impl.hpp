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
class CASADI_EXPORT FmuInternal : public SharedObjectInternal {
  friend class Fmu;
 public:
  // Constructor
  FmuInternal(const std::string& name,
    const std::vector<std::string>& scheme_in, const std::vector<std::string>& scheme_out,
    const std::map<std::string, std::vector<size_t>>& scheme, const std::vector<std::string>& aux);

  /// Destructor
  ~FmuInternal() override;

  // Initialize
  virtual void init(const DaeBuilderInternal* dae) = 0;

  /** \brief Print */
  void disp(std::ostream& stream, bool more) const override;

  /** \brief Get the number of scheme inputs */
  size_t n_in() const { return iind_.size();}

  /** \brief Get the number of scheme outputs */
  size_t n_out() const { return oind_.size();}

  // Index lookup for input
  size_t index_in(const std::string& n) const;

  // Index lookup for output
  size_t index_out(const std::string& n) const;

  // Does the FMU support analytic derivatives?
  virtual bool has_ad() const = 0;

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

  // Load an FMI function
  template<typename T>
  T* load_function(const std::string& symname);

  /** \brief Initalize memory block */
  virtual int init_mem(FmuMemory* m) const = 0;

  // Free FMU instance
  virtual void free_instance(void* c) const = 0;

  // Set value
  virtual void set(FmuMemory* m, size_t ind, const double* value) const = 0;

  // Request the calculation of a variable
  virtual void request(FmuMemory* m, size_t ind) const = 0;

  // Calculate all requested variables
  virtual int eval(FmuMemory* m) const = 0;

  // Get a calculated variable
  virtual void get(FmuMemory* m, size_t id, double* value) const = 0;

  // Set seed
  virtual void set_seed(FmuMemory* m, casadi_int nseed,
    const casadi_int* id, const double* v) const = 0;

  // Request the calculation of a sensitivity
  virtual void request_sens(FmuMemory* m, casadi_int nsens, const casadi_int* id,
    const casadi_int* wrt_id) const = 0;

  // Calculate directional derivatives
  virtual int eval_derivative(FmuMemory* m, bool independent_seeds) const = 0;

  // Get calculated derivatives
  virtual void get_sens(FmuMemory* m, casadi_int nsens,
    const casadi_int* id, double* v) const = 0;

  // Gather user sensitivities
  virtual void gather_sens(FmuMemory* m) const = 0;

  /** \brief Get stats */
  virtual void get_stats(FmuMemory* m, Dict* stats,
    const std::vector<std::string>& name_in, const InputStruct* in) const = 0;

 protected:

  /// Instance name
  std::string name_;

  // IO scheme
  std::vector<std::string> scheme_in_, scheme_out_;
  std::map<std::string, std::vector<size_t>> scheme_;

  // Auxilliary outputs
  std::vector<std::string> aux_;

  /// DLL
  Importer li_;

  // Mapping from scheme variable to and from FMU variable indices
  std::vector<size_t> iind_, iind_map_, oind_, oind_map_;

  // Meta information about the input/output variable subsets
  std::vector<double> nominal_in_, nominal_out_;
  std::vector<double> min_in_, min_out_;
  std::vector<double> max_in_, max_out_;
  std::vector<std::string> vn_in_, vn_out_;
  std::vector<unsigned int> vr_in_, vr_out_;

  // Numerical values for inputs
  std::vector<double> value_in_;

  // Reduced space indices for all inputs and outputs
  std::vector<std::vector<size_t>> ired_, ored_;

  // Sparsity pattern for extended Jacobian, Hessian
  Sparsity jac_sp_, hess_sp_;
};

template<typename T>
T* FmuInternal::load_function(const std::string& symname) {
  // Load the function
  signal_t f = li_.get_function(symname);
  // Ensure that it was found
  casadi_assert(f != 0, "Cannot retrieve '" + symname + "'");
  // Return function with the right type
  return reinterpret_cast<T*>(f);
}

/** \brief Interface to binary FMU

    Internal API.

    \author Joel Andersson
    \date 2023
*/
class CASADI_EXPORT Fmu2 : public FmuInternal {
 public:
  // Constructor
  Fmu2(const std::string& name,
    const std::vector<std::string>& scheme_in, const std::vector<std::string>& scheme_out,
    const std::map<std::string, std::vector<size_t>>& scheme, const std::vector<std::string>& aux);

  /// Destructor
  ~Fmu2() override;

  /** \brief Get type name */
  std::string class_name() const override { return "Fmu2";}

  // Initialize
  void init(const DaeBuilderInternal* dae) override;

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

  // Does the FMU support analytic derivatives?
  bool has_ad() const override { return get_directional_derivative_ != nullptr; }

  // New memory object
  fmi2Component instantiate() const;

  // Free FMU instance
  void free_instance(void* c) const override;

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
    const std::vector<std::string>& name_in, const InputStruct* in) const override;

  /** \brief Initalize memory block */
  int init_mem(FmuMemory* m) const override;

  // Set value
  void set(FmuMemory* m, size_t ind, const double* value) const override;

  // Request the calculation of a variable
  void request(FmuMemory* m, size_t ind) const override;

  // Gather user inputs and outputs
  void gather_io(FmuMemory* m) const;

  // Calculate all requested variables
  int eval(FmuMemory* m) const override;

  // Get a calculated variable
  void get(FmuMemory* m, size_t id, double* value) const override;

  // Set seed
  void set_seed(FmuMemory* m, casadi_int nseed,
    const casadi_int* id, const double* v) const override;

  // Request the calculation of a sensitivity
  void request_sens(FmuMemory* m, casadi_int nsens, const casadi_int* id,
    const casadi_int* wrt_id) const override;

  // Get calculated derivatives
  void get_sens(FmuMemory* m, casadi_int nsens,
    const casadi_int* id, double* v) const override;

  // Calculate directional derivatives using AD
  int eval_ad(FmuMemory* m) const;

  // Calculate directional derivatives using FD
  int eval_fd(FmuMemory* m, bool independent_seeds) const;

  // Calculate directional derivatives
  int eval_derivative(FmuMemory* m, bool independent_seeds) const override;

  // Gather user sensitivities
  void gather_sens(FmuMemory* m) const override;

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

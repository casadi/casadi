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

#include "fmu.hpp"
#include "importer.hpp"
#include "shared_object.hpp"

/// \cond INTERNAL

namespace casadi {

// Memory for event iteration
struct EventMemory {
  bool discrete_states_need_update;
  bool terminate_simulation;
  bool nominals_of_continuous_states_changed;
  bool values_of_continuous_states_changed;
  bool next_event_time_defined;
  double next_event_time;
};

// Forward declarations
class DaeBuilderInternal;
struct FmuMemory;
struct InputStruct;

/** \brief Interface to binary FMU

    Internal API.

    \author Joel Andersson
    \date 2023

    \identifier{26z} */
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
  virtual void init(const DaeBuilderInternal* dae);

  // Set C API functions
  virtual void load_functions() = 0;

  // Enter initialization mode
  virtual int enter_initialization_mode(void* instance) const = 0;

  // Exit initialization mode
  virtual int exit_initialization_mode(void* instance) const = 0;

  // Enter continuous-time mode
  virtual int enter_continuous_time_mode(void* instance) const = 0;

  // Update discrete states
  virtual int update_discrete_states(void* instance, EventMemory* eventmem) const = 0;

  // Set real values
  virtual int set_real(void* instance, const unsigned int* vr, size_t n_vr,
    const double* values, size_t n_values) const = 0;

  // Get/evaluate real values
  virtual int get_real(void* instance, const unsigned int* vr, size_t n_vr,
    double* values, size_t n_values) const = 0;

  // Forward mode AD
  virtual int get_directional_derivative(void* instance, const unsigned int* vr_out, size_t n_out,
    const unsigned int* vr_in, size_t n_in, const double* seed, size_t n_seed,
    double* sensitivity, size_t n_sensitivity) const = 0;

  // Reverse mode AD
  virtual int get_adjoint_derivative(void* instance, const unsigned int* vr_out, size_t n_out,
    const unsigned int* vr_in, size_t n_in, const double* seed, size_t n_seed,
    double* sensitivity, size_t n_sensitivity) const;

  // Copy values set in DaeBuilder to FMU
  virtual int set_values(void* instance) const = 0;

  // Retrieve auxilliary variables from FMU
  virtual int get_aux(void* instance) = 0;

  // Finalize
  void finalize();

  /** \brief Print

      \identifier{26m} */
  void disp(std::ostream& stream, bool more) const override;

  /** \brief Get the number of scheme inputs

      \identifier{270} */
  size_t n_in() const { return iind_.size();}

  /** \brief Get the number of scheme outputs

      \identifier{26o} */
  size_t n_out() const { return oind_.size();}

  // Index lookup for input
  size_t index_in(const std::string& n) const;

  // Index lookup for output
  size_t index_out(const std::string& n) const;

  // Get Jacobian sparsity for a subset of inputs and outputs
  Sparsity jac_sparsity(const std::vector<size_t>& osub, const std::vector<size_t>& isub) const;

  // Get Hessian sparsity for a subset of inputs
  Sparsity hess_sparsity(const std::vector<size_t>& r, const std::vector<size_t>& c) const;

  /// @{
  /** \brief Retreive nominal values

      \identifier{26r} */
  std::vector<double> all_nominal_in(size_t i) const;
  std::vector<double> all_nominal_out(size_t i) const;
  /// @}

  // Print description of an input
  std::string desc_in(FmuMemory* m, size_t id, bool more = true) const;

  // Name of system, per the FMI specification
  virtual std::string system_infix() const = 0;

  // DLL suffix, per the FMI specification
  static std::string dll_suffix();

  // Load an FMI function
  template<typename T>
  T* load_function(const std::string& symname);

  // Iteration to update discrete states
  int discrete_states_iter(void* instance) const;

  /** \brief Initalize memory block

      \identifier{271} */
  int init_mem(FmuMemory* m) const;

  // New memory object
  virtual void* instantiate() const = 0;

  // Free FMU instance
  virtual void free_instance(void* c) const = 0;

  // Set value
  void set(FmuMemory* m, size_t ind, const double* value) const;

  // Request the calculation of a variable
  void request(FmuMemory* m, size_t ind) const;

  // Calculate all requested variables
  int eval(FmuMemory* m) const;

  // Get a calculated variable
  void get(FmuMemory* m, size_t id, double* value) const;

  // Set forward seeds
  void set_fwd(FmuMemory* m, casadi_int nseed,
    const casadi_int* id, const double* v) const;

  // Set all forward seeds for a single input
  void set_fwd(FmuMemory* m, size_t ind, const double* v) const;

  // Request the calculation of forward sensitivities
  void request_fwd(FmuMemory* m, casadi_int nsens, const casadi_int* id,
    const casadi_int* wrt_id) const;

  // Request the calculation of all forward sensitivities for an output
  void request_fwd(FmuMemory* m, casadi_int ind) const;

  // Calculate forward directional derivatives
  int eval_fwd(FmuMemory* m, bool independent_seeds) const;

  // Calculate forward directional derivatives using AD
  int eval_ad(FmuMemory* m) const;

  // Calculate forward directional derivatives using FD
  int eval_fd(FmuMemory* m, bool independent_seeds) const;

  // Get forward  sensitivities
  void get_fwd(FmuMemory* m, casadi_int nsens,
    const casadi_int* id, double* v) const;

  // Get the forward sensitivities for a single output
  void get_fwd(FmuMemory* m, size_t ind, double* v) const;

  // Set adjoint seeds
  void set_adj(FmuMemory* m, casadi_int nseed,
    const casadi_int* id, const double* v) const;

  // Set all adjoint seeds for a single output
  void set_adj(FmuMemory* m, size_t ind, const double* v) const;

  // Request the calculation of adjoint sensitivities
  void request_adj(FmuMemory* m, casadi_int nsens, const casadi_int* id,
    const casadi_int* wrt_id) const;

  // Request the calculation of all adjoint sensitivities for an input
  void request_adj(FmuMemory* m, casadi_int ind) const;

  // Calculate adjoint sensitivities
  int eval_adj(FmuMemory* m) const;

  // Get adjoint sensitivities
  void get_adj(FmuMemory* m, casadi_int nsens,
    const casadi_int* id, double* v) const;

  // Get the adjoint sensitivities for a single input
  void get_adj(FmuMemory* m, size_t ind, double* v) const;

  // Gather forward sensitivities
  void gather_fwd(FmuMemory* m) const;

  // Gather adjoint sensitivities
  void gather_adj(FmuMemory* m) const;

  // Gather user inputs and outputs
  void gather_io(FmuMemory* m) const;

  /** \brief Get stats

      \identifier{272} */
  virtual void get_stats(FmuMemory* m, Dict* stats,
    const std::vector<std::string>& name_in, const InputStruct* in) const = 0;

  void serialize(SerializingStream& s) const;

  virtual void serialize_type(SerializingStream& s) const;
  virtual void serialize_body(SerializingStream& s) const;

  static FmuInternal* deserialize(DeserializingStream& s);

 protected:
  explicit FmuInternal(DeserializingStream& s);

  /// Instance name
  std::string name_;

  // IO scheme
  std::vector<std::string> scheme_in_, scheme_out_;
  std::map<std::string, std::vector<size_t>> scheme_;

  // Auxilliary outputs
  std::vector<std::string> aux_;

  // Path to the FMU resource directory
  std::string resource_loc_;

  // Tolerance
  double fmutol_;

  // Instance name
  std::string instance_name_;

  // GUID / instantiation_token
  std::string instantiation_token_;

  // Logging?
  bool logging_on_;

  // Number of event indicators
  casadi_int number_of_event_indicators_;

  // Does the FMU declare analytic derivatives support?
  bool provides_directional_derivatives_, provides_adjoint_derivatives_;

  // Does the FMU declare restrictions on instantiation?
  bool can_be_instantiated_only_once_per_process_;

  /// DLL
  Importer li_;

  // Mapping from scheme variable to and from FMU variable indices
  std::vector<size_t> iind_, iind_map_, oind_, oind_map_;

  // Is there an independent variable?
  bool has_independent_;

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

} // namespace casadi

/// \endcond

#endif // CASADI_FMU_IMPL_HPP

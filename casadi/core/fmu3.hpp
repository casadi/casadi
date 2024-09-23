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

#ifndef CASADI_FMU3_HPP
#define CASADI_FMU3_HPP

#include "fmu_impl.hpp"

#include <fmi3Functions.h>

/// \cond INTERNAL

namespace casadi {

/** \brief Interface to a binary FMU, adhering to FMI version 2.0.

    \author Joel Andersson
    \date 2024
  */
class CASADI_EXPORT Fmu3 : public FmuInternal {
 public:
  // Constructor
  Fmu3(const std::string& name,
    const std::vector<std::string>& scheme_in, const std::vector<std::string>& scheme_out,
    const std::map<std::string, std::vector<size_t>>& scheme, const std::vector<std::string>& aux);

  /// Destructor
  ~Fmu3() override;

  /** \brief Get type name */
  std::string class_name() const override { return "Fmu3";}

  // Initialize
  void init(const DaeBuilderInternal* dae) override;

  // Set C API functions
  void load_functions() override;

  // Variables used for initialization, by type
  std::vector<fmi3ValueReference> vr_real_, vr_integer_, vr_boolean_, vr_string_;
  std::vector<fmi3Float64> init_real_;
  std::vector<fmi3Int32> init_integer_;
  std::vector<fmi3Boolean> init_boolean_;
  std::vector<std::string> init_string_;

  // Auxilliary variables, by type
  std::vector<std::string> vn_aux_real_, vn_aux_integer_, vn_aux_boolean_, vn_aux_string_;
  std::vector<fmi3ValueReference> vr_aux_real_, vr_aux_integer_, vr_aux_boolean_, vr_aux_string_;

  // Following members set in finalize

  // FMU C API function prototypes. Cf. FMI specification 2.0.2
  fmi3InstantiateModelExchangeTYPE* instantiate_model_exchange_;
  fmi3FreeInstanceTYPE* free_instance_;
  fmi3ResetTYPE* reset_;
  fmi3EnterInitializationModeTYPE* enter_initialization_mode_;
  fmi3ExitInitializationModeTYPE* exit_initialization_mode_;
  fmi3EnterContinuousTimeModeTYPE* enter_continuous_time_mode_;
  fmi3GetFloat64TYPE* get_float64_;
  fmi3SetFloat64TYPE* set_float64_;
  fmi3GetBooleanTYPE* get_boolean_;
  fmi3SetBooleanTYPE* set_boolean_;
  fmi3GetInt32TYPE* get_int32_;
  fmi3SetInt32TYPE* set_int32_;
  fmi3GetStringTYPE* get_string_;
  fmi3SetStringTYPE* set_string_;
  fmi3GetDirectionalDerivativeTYPE* get_directional_derivative_;
  fmi3GetAdjointDerivativeTYPE* get_adjoint_derivative_;

  // Collection of variable values, all types
  struct Value {
    std::vector<fmi3Float64> v_real;
    std::vector<fmi3Int32> v_integer;
    std::vector<fmi3Boolean> v_boolean;
    std::vector<std::string> v_string;
  };

  Value aux_value_;

  // Name of system, per the FMI specification
  std::string system_infix() const override;

  // New memory object
  void* instantiate() const override;

  // Free FMU instance
  void free_instance(void* instance) const override;

  // Reset solver
  int reset(void* instance);

  // Enter initialization mode
  int enter_initialization_mode(void* instance) const override;

  // Exit initialization mode
  int exit_initialization_mode(void* instance) const override;

  // Set real values
  int set_real(void* instance, const unsigned int* vr, size_t n_vr,
    const double* values, size_t n_values) const override;

  // Get/evaluate real values
  int get_real(void* instance, const unsigned int* vr, size_t n_vr,
    double* values, size_t n_values) const override;

  // Forward mode AD
  int get_directional_derivative(void* instance, const unsigned int* vr_out, size_t n_out,
    const unsigned int* vr_in, size_t n_in, const double* seed, size_t n_seed,
    double* sensitivity, size_t n_sensitivity) const override;

  // Reverse mode AD
  int get_adjoint_derivative(void* instance, const unsigned int* vr_out, size_t n_out,
    const unsigned int* vr_in, size_t n_in, const double* seed, size_t n_seed,
    double* sensitivity, size_t n_sensitivity) const override;

  // Copy values set in DaeBuilder to FMU
  int set_values(void* instance) const override;

  // Retrieve auxilliary variables from FMU
  int get_aux(void* instance) override;

  /** \brief Get stats */
  void get_stats(FmuMemory* m, Dict* stats,
    const std::vector<std::string>& name_in, const InputStruct* in) const override;

  // Process message
  static void log_message_callback(fmi3InstanceEnvironment instanceEnvironment,
    fmi3Status status, fmi3String category, fmi3String message);
};

} // namespace casadi

/// \endcond

#endif // CASADI_FMU3_HPP

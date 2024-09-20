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

#ifndef CASADI_FMU2_HPP
#define CASADI_FMU2_HPP

#include "fmu_impl.hpp"

#include <fmi2Functions.h>

/// \cond INTERNAL

namespace casadi {

/** \brief Interface to a binary FMU, adhering to FMI version 2.0.

    \author Joel Andersson
    \date 2023

    \identifier{273} */
class CASADI_EXPORT Fmu2 : public FmuInternal {
 public:
  // Constructor
  Fmu2(const std::string& name,
    const std::vector<std::string>& scheme_in, const std::vector<std::string>& scheme_out,
    const std::map<std::string, std::vector<size_t>>& scheme, const std::vector<std::string>& aux);

  /// Destructor
  ~Fmu2() override;

  /** \brief Get type name

      \identifier{274} */
  std::string class_name() const override { return "Fmu2";}

  // Initialize
  void init(const DaeBuilderInternal* dae) override;

  // Finalize
  void finalize() override;

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

  // Variables used for initialization, by type
  std::vector<fmi2ValueReference> vr_real_, vr_integer_, vr_boolean_, vr_string_;
  std::vector<fmi2Real> init_real_;
  std::vector<fmi2Integer> init_integer_;
  std::vector<fmi2Boolean> init_boolean_;
  std::vector<std::string> init_string_;

  // Auxilliary variables, by type
  std::vector<std::string> vn_aux_real_, vn_aux_integer_, vn_aux_boolean_, vn_aux_string_;
  std::vector<fmi2ValueReference> vr_aux_real_, vr_aux_integer_, vr_aux_boolean_, vr_aux_string_;

  // Does the FMU declare analytic derivatives support?
  bool declared_ad_;

  // Following members set in finalize

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

  // Collection of variable values, all types
  struct Value {
    std::vector<fmi2Real> v_real;
    std::vector<fmi2Integer> v_integer;
    std::vector<fmi2Boolean> v_boolean;
    std::vector<std::string> v_string;
  };

  Value aux_value_;

  // Does the FMU support analytic forward derivatives?
  bool has_fwd() const override { return get_directional_derivative_ != nullptr; }

  // New memory object
  void* instantiate() const;

  // Free FMU instance
  void free_instance(void* instance) const override;

  // Reset solver
  int reset(void* instance);

  // Enter initialization mode
  int enter_initialization_mode(void* instance) const;

  // Exit initialization mode
  int exit_initialization_mode(void* instance) const;

  // Copy values set in DaeBuilder to FMU
  int set_values(void* instance) const;

  // Retrieve input variable values from FMU
  int get_in(void* instance, std::vector<fmi2Real>* v) const;

  // Retrieve auxilliary variables from FMU
  int get_aux(void* instance, Value* v) const;

  /** \brief Get stats

      \identifier{275} */
  void get_stats(FmuMemory* m, Dict* stats,
    const std::vector<std::string>& name_in, const InputStruct* in) const override;

  /** \brief Initalize memory block

      \identifier{276} */
  int init_mem(FmuMemory* m) const override;

  // Calculate all requested variables
  int eval(FmuMemory* m) const override;

  // Calculate directional derivatives using AD
  int eval_ad(FmuMemory* m) const override;

  // Calculate directional derivatives using FD
  int eval_fd(FmuMemory* m, bool independent_seeds) const override;

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

  void serialize_body(SerializingStream& s) const override;

  static Fmu2* deserialize(DeserializingStream& s);

  protected:
    explicit Fmu2(DeserializingStream& s);
};

} // namespace casadi

/// \endcond

#endif // CASADI_FMU2_HPP

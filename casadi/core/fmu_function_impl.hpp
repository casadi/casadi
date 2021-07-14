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


#ifndef CASADI_FMU_FUNCTION_IMPL_HPP
#define CASADI_FMU_FUNCTION_IMPL_HPP

#include "fmu_function.hpp"
#include "function_internal.hpp"
#include "importer.hpp"

#ifdef WITH_FMU
//#include <fmi2FunctionTypes.h>
#include <fmi2Functions.h>
//#include <fmi2TypesPlatform.h>
#endif  // WITH_FMU

#ifdef WITH_DL
#ifdef _WIN32 // also for 64-bit
#ifndef NOMINMAX
#define NOMINMAX
#endif
#ifndef _WIN32_WINNT
#define _WIN32_WINNT 0x0502
#endif
#include <windows.h>
#else // _WIN32
#include <dlfcn.h>
#endif // _WIN32
#endif // WITH_DL

/// \cond INTERNAL

namespace casadi {

#ifdef WITH_FMU

class CASADI_EXPORT FmuFunction : public FunctionInternal {
 protected:
  // Path to the unpacked FMU
  std::string path_;

  // Value reference to the inputs and outputs
  std::vector<std::vector<casadi_int>> id_in_, id_out_;

  // Global unique identifier
  const std::string guid_;

  /** \brief Information about the library */
  Importer li_;

  // Callback functions
  fmi2CallbackFunctions functions_;

  // Path to the FMU resource directory
  std::string resource_loc_;

  // User-set options
  bool provides_directional_derivative_;
  std::string instance_name_;

  // Component references
  std::vector<fmi2ValueReference> xd_, xn_, yd_, yn_;

  // FMU C API function prototypes. Cf. FMI specification 2.0.2
  fmi2InstantiateTYPE* instantiate_;
  fmi2FreeInstanceTYPE* free_instance_;
  fmi2SetupExperimentTYPE* setup_experiment_;
  fmi2EnterInitializationModeTYPE* enter_initialization_mode_;
  fmi2ExitInitializationModeTYPE* exit_initialization_mode_;
  fmi2EnterContinuousTimeModeTYPE* enter_continuous_time_mode_;
  fmi2SetRealTYPE* set_real_;
  fmi2SetBooleanTYPE* set_boolean_;
  fmi2GetRealTYPE* get_real_;
  fmi2GetDirectionalDerivativeTYPE* get_directional_derivative_;

  // Pointer to instance
  fmi2Component c_;

 public:

  /** \brief Constructor */
  FmuFunction(const std::string& name, const std::string& path,
      const std::vector<std::vector<casadi_int>>& id_in,
      const std::vector<std::vector<casadi_int>>& id_out,
      const std::string& guid);

  /** \brief Destructor */
  ~FmuFunction() override;

  /** \brief Get type name */
  std::string class_name() const override { return "FmuFunction";}

  ///@{
  /** \brief Options */
  static const Options options_;
  const Options& get_options() const override { return options_;}
  ///@}

  /// Initialize
  void init(const Dict& opts) override;

  ///@{
  /** \brief Number of function inputs and outputs */
  size_t get_n_in() override { return 2;}
  size_t get_n_out() override {return 2;}
  ///@}

  ///@{
  /** \brief Names of function input and outputs */
  std::string get_name_in(casadi_int i) override { return i == 0 ? "xd" : "xn";}
  std::string get_name_out(casadi_int i) override { return i == 0 ? "yd" : "yn";}
  /// @}

  /// @{
  /** \brief Retreive differentiability */
  bool get_diff_in(casadi_int i) override {return i == 0;}
  bool get_diff_out(casadi_int i) override {return i == 0;}
  /// @}

  /// @{
  /** \brief Retreive sparsities */
  Sparsity get_sparsity_in(casadi_int i) override;
  Sparsity get_sparsity_out(casadi_int i) override;
  /// @}

  // Evaluate numerically
  int eval(const double** arg, double** res, casadi_int* iw, double* w, void* mem) const override;

  // Evaluate Jacobian numerically
  int eval_jac(const double** arg, double** res, casadi_int* iw, double* w, void* mem) const;

  // Evaluate adjoint numerically
  int eval_adj(const double** arg, double** res, casadi_int* iw, double* w, void* mem) const;

  ///@{
  /** \brief Full Jacobian */
  bool has_jacobian() const override { return provides_directional_derivative_;}
  Function get_jacobian(const std::string& name,
    const std::vector<std::string>& inames,
    const std::vector<std::string>& onames,
    const Dict& opts) const override;
  ///@}

  ///@{
  /** \brief Adjoint */
  bool has_reverse(casadi_int nadj) const override {
    return provides_directional_derivative_ && nadj == 1;
  }
  ///@}
  Function get_reverse(casadi_int nadj, const std::string& name,
    const std::vector<std::string>& inames,
    const std::vector<std::string>& onames,
    const Dict& opts) const override;

  // Name of system, per the FMI specification
  static std::string system_infix();

  // DLL suffix, per the FMI specification
  static std::string dll_suffix();

  // Load an FMI function
  signal_t get_function(const std::string& symname);

  // Process message
  static void logger(fmi2ComponentEnvironment componentEnvironment,
    fmi2String instanceName,
    fmi2Status status,
    fmi2String category,
    fmi2String message, ...);
};

/** Jacobian */
class CASADI_EXPORT FmuFunctionJac : public FunctionInternal {
 public:
  /// Constructor
  FmuFunctionJac(const std::string& name) : FunctionInternal(name) {}

  /// Destructor
  ~FmuFunctionJac() override;

  /// Initialize
  void init(const Dict& opts) override;

  /** \brief Get type name */
  std::string class_name() const override { return "FmuFunctionJac";}

  /// Evaluate numerically
  int eval(const double** arg, double** res, casadi_int* iw, double* w, void* mem) const override;
};

/** Adjoint */
class CASADI_EXPORT FmuFunctionAdj : public FunctionInternal {
 public:
  /// Constructor
  FmuFunctionAdj(const std::string& name) : FunctionInternal(name) {}

  /// Destructor
  ~FmuFunctionAdj() override;

  /// Initialize
  void init(const Dict& opts) override;

  /** \brief Get type name */
  std::string class_name() const override { return "FmuFunctionAdj";}

  /// Evaluate numerically
  int eval(const double** arg, double** res, casadi_int* iw, double* w, void* mem) const override;
};

#endif  // WITH_FMU

} // namespace casadi
/// \endcond

#endif // CASADI_FMU_FUNCTION_IMPL_HPP

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

#include "function_internal.hpp"
#include "dae_builder_internal.hpp"


/// \cond INTERNAL

namespace casadi {

#ifdef WITH_FMU

struct CASADI_EXPORT FmuFunctionMemory : public FunctionMemory {
  // Pointer to instance
  fmi2Component c;

  // First run
  bool first_run;
};


class CASADI_EXPORT FmuFunction : public FunctionInternal {
 protected:
  // DaeBuilder instance (non-owning reference to avoid circular dependency)
  WeakRef dae_;

  // Variable references for inputs and outputs
  std::vector<std::vector<size_t>> id_in_, id_out_;

  // Value reference to the inputs and outputs
  std::vector<std::vector<fmi2ValueReference>> vref_in_, vref_out_;

  // Callback functions
  fmi2CallbackFunctions functions_;

  // Path to the FMU resource directory
  std::string resource_loc_;

 public:

  /** \brief Constructor */
  FmuFunction(const std::string& name, const DaeBuilder& dae,
      const std::vector<std::vector<size_t>>& id_in,
      const std::vector<std::vector<size_t>>& id_out,
      const std::vector<std::string>& name_in,
      const std::vector<std::string>& name_out);

  /** \brief Destructor */
  ~FmuFunction() override;

  /** \brief Get type name */
  std::string class_name() const override { return "FmuFunction";}

  /// Initialize
  void init(const Dict& opts) override;

  /** \brief Create memory block */
  void* alloc_mem() const override { return new FmuFunctionMemory();}

  /** \brief Initalize memory block */
  int init_mem(void* mem) const override;

  /** \brief Free memory block */
  void free_mem(void *mem) const override;

  ///@{
  /** \brief Number of function inputs and outputs */
  size_t get_n_in() override { return id_in_.size();}
  size_t get_n_out() override {return id_out_.size();}
  ///@}

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
  bool has_jacobian() const override;
  Function get_jacobian(const std::string& name,
    const std::vector<std::string>& inames,
    const std::vector<std::string>& onames,
    const Dict& opts) const override;
  ///@}

  ///@{
  /** \brief Adjoint */
  bool has_reverse(casadi_int nadj) const override;
  ///@}
  Function get_reverse(casadi_int nadj, const std::string& name,
    const std::vector<std::string>& inames,
    const std::vector<std::string>& onames,
    const Dict& opts) const override;

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

  // Pass inputs to FMU
  int set_inputs(FmuFunctionMemory* m, const double** x) const;

  // Get/calculate outputs
  int get_outputs(FmuFunctionMemory* m, double** r) const;
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

  /** \brief Create memory block */
  void* alloc_mem() const override { return new FmuFunctionMemory();}

  /** \brief Initalize memory block */
  int init_mem(void* mem) const override;

  /** \brief Free memory block */
  void free_mem(void *mem) const override;

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

  /** \brief Create memory block */
  void* alloc_mem() const override { return new FmuFunctionMemory();}

  /** \brief Initalize memory block */
  int init_mem(void* mem) const override;

  /** \brief Free memory block */
  void free_mem(void *mem) const override;

  /** \brief Get type name */
  std::string class_name() const override { return "FmuFunctionAdj";}

  /// Evaluate numerically
  int eval(const double** arg, double** res, casadi_int* iw, double* w, void* mem) const override;
};

#endif  // WITH_FMU

} // namespace casadi
/// \endcond

#endif // CASADI_FMU_FUNCTION_IMPL_HPP

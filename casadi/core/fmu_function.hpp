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

class CASADI_EXPORT FmuFunction : public FunctionInternal {
 public:
  // DaeBuilder instance (non-owning reference to avoid circular dependency)
  WeakRef dae_;

  // Variable references for inputs and outputs
  std::vector<std::vector<size_t>> id_in_, id_out_;

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

  // Evaluate, non-differentated
  int eval_fmu(const DaeBuilderInternal* dae, int mem, const double** arg, double** res,
    const std::vector<std::vector<size_t>>& id_in,
    const std::vector<std::vector<size_t>>& id_out) const;

  // Evaluate Jacobian numerically
  int eval_jac(const DaeBuilderInternal* dae, int mem, const double** arg, double** res,
    const std::vector<std::vector<size_t>>& id_in,
    const std::vector<std::vector<size_t>>& id_out) const;

  // Evaluate adjoint numerically
  int eval_adj(const DaeBuilderInternal* dae, int mem, const double** arg, double** res,
    const std::vector<std::vector<size_t>>& id_in,
    const std::vector<std::vector<size_t>>& id_out) const;

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
};

/** Jacobian */
class CASADI_EXPORT FmuFunctionJac : public FunctionInternal {
 public:
  /// Constructor
  FmuFunctionJac(const std::string& name) : FunctionInternal(name) {}

  /// Destructor
  ~FmuFunctionJac() override;

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

  /** \brief Get type name */
  std::string class_name() const override { return "FmuFunctionAdj";}

  /// Evaluate numerically
  int eval(const double** arg, double** res, casadi_int* iw, double* w, void* mem) const override;
};

#endif  // WITH_FMU

} // namespace casadi
/// \endcond

#endif // CASADI_FMU_FUNCTION_IMPL_HPP

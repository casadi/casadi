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


#ifndef CASADI_FMU_HPP
#define CASADI_FMU_HPP

#include "shared_object.hpp"
#include "printable.hpp"
#include "sparsity.hpp"

/// \cond INTERNAL

namespace casadi {

// Forward declarations
class DaeBuilderInternal;
class FmuFunction;
struct FmuMemory;
struct InputStruct;

// Forward declaration of internal class
class FmuInternal;

/// Which C API
enum class FmuApi {FMI2, FMI3, NUMEL};

/// Convert to string
CASADI_EXPORT std::string to_string(FmuApi v);

/** \brief Interface to binary FMU

    Can be shared between multiple CasADi functions.

    \author Joel Andersson
    \date 2023

    \identifier{26s} */
class CASADI_EXPORT Fmu
  : public SharedObject,
    public SWIG_IF_ELSE(PrintableCommon, Printable<Fmu>) {
 public:
  /** \brief Get type name

      \identifier{26t} */
  static std::string type_name() {return "Fmu";}

  /// Default constructor
  Fmu();

  /// Importer factory
  Fmu(const std::string& name, FmuApi api, const DaeBuilderInternal* dae,
    const std::vector<std::string>& scheme_in,
    const std::vector<std::string>& scheme_out,
    const std::map<std::string, std::vector<size_t>>& scheme,
    const std::vector<std::string>& aux);

  /** \brief Name of the instance

      \identifier{26u} */
  const std::string& name() const;

  /** \brief Name of the FMU

      \identifier{2ab} */
  const std::string& instance_name() const;

  ///@{
  /// Access functions of the node
  FmuInternal* operator->();
  const FmuInternal* operator->() const;
  FmuInternal* get() const;
  ///@}

  /** \brief Get the number of scheme inputs

      \identifier{26v} */
  size_t n_in() const;

  /** \brief Get the number of scheme outputs

      \identifier{26w} */
  size_t n_out() const;

  // Index lookup for input
  size_t index_in(const std::string& n) const;

  // Index lookup for output
  size_t index_out(const std::string& n) const;

  // Get all reduced space indices for an input
  const std::vector<size_t>& ired(size_t ind) const;

  // Get all reduced space indices for an output
  const std::vector<size_t>& ored(size_t ind) const;

  // Get nominal value for reduced-space input
  double nominal_in(size_t ind) const;

  // Get nominal value for reduced-space output
  double nominal_out(size_t ind) const;

  // Get minimum value for reduced-space input
  double min_in(size_t ind) const;

  // Get maximum value for reduced-space input
  double max_in(size_t ind) const;

  // Get all nominal values for an input
  std::vector<double> all_nominal_in(size_t ind) const;

  // Get all nominal values for an output
  std::vector<double> all_nominal_out(size_t ind) const;

  // Description of an input
  std::string desc_in(FmuMemory* m, size_t id, bool more = true) const;

  /// Does the FMU provide support for forward directional derivatives
  bool provides_directional_derivatives() const;

  /// Does the FMU provide support for adjoint directional derivatives
  bool provides_adjoint_derivatives() const;

  /// Does the FMU declare restrictions on instantiation?
  bool can_be_instantiated_only_once_per_process() const;

  // Get Jacobian sparsity for a subset of inputs and outputs
  Sparsity jac_sparsity(const std::vector<size_t>& osub, const std::vector<size_t>& isub) const;

  // Get Jacobian sparsity for an output/input pair
  Sparsity jac_sparsity(size_t oind, size_t iind) const {
    return jac_sparsity(ored(oind), ired(iind));
  }

  // Get Hessian sparsity for a subset of inputs
  Sparsity hess_sparsity(const std::vector<size_t>& r, const std::vector<size_t>& c) const;

  // Get Jacobian sparsity for an input/input pair
  Sparsity hess_sparsity(size_t r, size_t c) const {
    return hess_sparsity(ired(r), ired(c));
  }

  /** \brief Create memory block

      \identifier{2dn} */
  FmuMemory* alloc_mem(const FmuFunction& f) const;

  /** \brief Initalize memory block

      \identifier{26x} */
  int init_mem(FmuMemory* m) const;

  /** \brief Free memory block

      \identifier{2dt} */
  void free_mem(void *mem) const;

  // Free FMU instance
  void free_instance(void* instance) const;

  // Set value
  void set(FmuMemory* m, size_t ind, const double* value) const;

  // Request the calculation of a variable
  void request(FmuMemory* m, size_t ind) const;

  // Calculate all requested variables
  int eval(FmuMemory* m) const;

  // Get a calculated variable
  void get(FmuMemory* m, size_t ind, double* value) const;

  // Set forward seed by variable index
  void set_fwd(FmuMemory* m, casadi_int nseed, const casadi_int* id, const double* v) const;

  // Set all forward seeds for a single input
  void set_fwd(FmuMemory* m, size_t ind, const double* v) const;

  // Request the calculation of forward sensitivities
  void request_fwd(FmuMemory* m, casadi_int nsens, const casadi_int* id,
    const casadi_int* wrt_id) const;

  // Request the calculation of all forward sensitivities for an output
  void request_fwd(FmuMemory* m, casadi_int ind) const;

  // Calculate forward directional derivatives
  int eval_fwd(FmuMemory* m, bool independent_seeds) const;

  // Get forward sensitivities
  void get_fwd(FmuMemory* m, casadi_int nsens, const casadi_int* id, double* v) const;

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

  // Calculate adjoint directional derivatives
  int eval_adj(FmuMemory* m) const;

  // Get adjoint sensitivities
  void get_adj(FmuMemory* m, casadi_int nsens,
    const casadi_int* id, double* v) const;

  // Get the adjoint sensitivities for a single input
  void get_adj(FmuMemory* m, size_t ind, double* v) const;

  /** \brief Get stats

      \identifier{26y} */
  void get_stats(FmuMemory* m, Dict* stats,
    const std::vector<std::string>& name_in, const InputStruct* in) const;

  /// \cond INTERNAL
#ifndef SWIG
    /** \brief  Create from node

        \identifier{27c} */
    static Fmu create(FmuInternal* node);
#endif // SWIG
  /// \endcond

  /** \brief Serialize an object

      \identifier{27d} */
  void serialize(SerializingStream &s) const;

  /** \brief Deserialize with type disambiguation

      \identifier{27e} */
  static Fmu deserialize(DeserializingStream& s);
};

} // namespace casadi

/// \endcond

#endif // CASADI_FMU_HPP

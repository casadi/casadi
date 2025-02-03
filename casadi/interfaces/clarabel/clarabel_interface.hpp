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

#ifndef CASADI_CLARABEL_INTERFACE_HPP
#define CASADI_CLARABEL_INTERFACE_HPP

#include "casadi/core/conic_impl.hpp"
#include <casadi/interfaces/clarabel/casadi_conic_clarabel_export.h>
#include <Clarabel.h>
#include <string>

namespace casadi {
  #include "clarabel_runtime.hpp"
}

/** \defgroup plugin_Conic_clarabel Title
    \par
      Interface to the Clarabel solver for sparse Quadratic Programs.
      See the Clarabel documentation (e.g. clarabel.dev) for more details.
*/

/** \pluginsection{Conic,clarabel} */

/// \cond INTERNAL

namespace casadi {

  struct CASADI_CONIC_CLARABEL_EXPORT ClarabelMemory : public ConicMemory {
    // Problem data structure (Clarabel-specific)
    casadi_clarabel_data<double> d;
  };

  /** \brief \pluginbrief{Conic,clarabel}
      @copydoc Clarabel_doc
      @copydoc plugin_Conic_clarabel

      \author Joris Gillis
      \date 2025
  */
  class CASADI_CONIC_CLARABEL_EXPORT ClarabelInterface : public Conic {
  public:
    /// Creator function for the plugin
    static Conic* creator(const std::string& name,
                          const std::map<std::string, Sparsity>& st) {
      return new ClarabelInterface(name, st);
    }

    /// Constructor using sparsity patterns
    explicit ClarabelInterface(const std::string& name,
                               const std::map<std::string, Sparsity>& st);

    /// Destructor
    ~ClarabelInterface() override;

    /// Get the name of the plugin
    const char* plugin_name() const override { return "clarabel"; }

    /// Get the class name
    std::string class_name() const override { return "ClarabelInterface"; }

    /// Options for Clarabel
    static const Options options_;
    const Options& get_options() const override { return options_; }

    void set_clarabel_prob();
    void set_clarabel_prob(CodeGenerator& g) const;

    /// Generate code for the function body
    void codegen_body(CodeGenerator& g) const override;

    /// Codegen for initializing memory
    void codegen_init_mem(CodeGenerator& g) const override;

    /// Codegen for freeing memory
    void codegen_free_mem(CodeGenerator& g) const override;

    /// Thread-local memory object type
    std::string codegen_mem_type() const override { return "struct casadi_clarabel_data"; }

    /// Initialize the solver
    void init(const Dict& opts) override;

    /// Initialize dependent (read-only) members
    void init_dependent();

    /// Create memory block
    void* alloc_mem() const override { return new ClarabelMemory(); }

    /// Initialize the memory block
    int init_mem(void* mem) const override;

    /// Free the memory block
    void free_mem(void* mem) const override;

    /// Set the (persistent) work vectors
    void set_work(void* mem, const double**& arg, double**& res,
                  casadi_int*& iw, double*& w) const override;

    /// Get all statistics
    Dict get_stats(void* mem) const override;

    /// Solve the QP
    int solve(const double** arg, double** res,
              casadi_int* iw, double* w, void* mem) const override;

    /// A documentation string
    static const std::string meta_doc;

    /// Options passed to Clarabel (distinct from options in the solver)
    Dict opts_;

    void serialize_body(SerializingStream &s) const override;

    /** \brief Deserialize with type disambiguation */
    static ProtoFunction* deserialize(DeserializingStream& s) { return new ClarabelInterface(s); }

    /// Clarabel does not (in this example) support discrete variables.
    bool integer_support() const override { return false; }

  protected:
    /** \brief Deserializing constructor */
    explicit ClarabelInterface(DeserializingStream& s);

  private:
    // Problem data (constructed from the CasADi QP data)
    casadi_clarabel_prob<double> p_;

    // Sparsity information for the quadratic cost (P) and constraints (A)
    std::vector<int> colindp_, rowp_;
    std::vector<int> colinda_, rowa_;

    // The underlying QP data from CasADi
    casadi_qp_prob<double> p_qp_;
    casadi_qp_data<double> d_qp_;
  };
} // end namespace casadi
/// \endcond
#endif // CASADI_CLARABEL_INTERFACE_HPP

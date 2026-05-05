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

#ifndef CASADI_XPRESS_INTERFACE_HPP
#define CASADI_XPRESS_INTERFACE_HPP

#include "casadi/core/conic_impl.hpp"
#include <casadi/interfaces/xpress/casadi_conic_xpress_export.h>

#include <xprs.h>
#include <string>

namespace casadi {
  #include "xpress_runtime.hpp"
}

/** \defgroup plugin_Conic_xpress Title
    \par

      Interface to FICO Xpress Optimizer for sparse Linear, Quadratic and
      Mixed-Integer Quadratic Programs.  See https://www.fico.com/en/products/fico-xpress-optimization
      for more information and the Xpress reference manual for the list of
      controls (options) accepted under the ``xpress`` option dict.

    \identifier{xpress} */

/** \pluginsection{Conic,xpress} */

/// \cond INTERNAL

namespace casadi {

  struct CASADI_CONIC_XPRESS_EXPORT XpressMemory : public ConicMemory {
    // Problem data structure
    casadi_xpress_data<double> d;
  };

  /** \brief \pluginbrief{Conic,xpress}

      @copydoc Xpress_doc
      @copydoc plugin_Conic_xpress

      \author Joris Gillis
      \date 2026
  */
  class CASADI_CONIC_XPRESS_EXPORT XpressInterface : public Conic {
  public:
    /** \brief  Create a new QP Solver */
    static Conic* creator(const std::string& name,
                          const std::map<std::string, Sparsity>& st) {
      return new XpressInterface(name, st);
    }

    /// Constructor using sparsity patterns
    explicit XpressInterface(const std::string& name,
                             const std::map<std::string, Sparsity>& st);

    /// Destructor
    ~XpressInterface() override;

    // Get name of the plugin
    const char* plugin_name() const override { return "xpress";}

    // Get name of the class
    std::string class_name() const override { return "XpressInterface";}

    ///@{
    /** \brief Options */
    static const Options options_;
    const Options& get_options() const override { return options_;}
    ///@}

    void set_xpress_prob();
    void set_xpress_prob(CodeGenerator& g) const;

    /** \brief Generate code for the function body */
    void codegen_body(CodeGenerator& g) const override;

    /** \brief Codegen init for init_mem */
    void codegen_init_mem(CodeGenerator& g) const override;

    /** \brief Codegen for free_mem */
    void codegen_free_mem(CodeGenerator& g) const override;

    /** \brief Thread-local memory object type */
    std::string codegen_mem_type() const override { return "struct casadi_xpress_data"; }

    /** \brief Is thread-local memory object needed? */
    bool codegen_needs_mem() const override { return true; }

    // Initialize the solver
    void init(const Dict& opts) override;

    // Initialize dependant read-only members of class
    void init_dependent();

    /** \brief Create memory block */
    void* alloc_mem() const override { return new XpressMemory();}

    /** \brief Initalize memory block */
    int init_mem(void* mem) const override;

    /** \brief Free memory block */
    void free_mem(void *mem) const override;

    /** \brief Set the (persistent) work vectors */
    void set_work(void* mem, const double**& arg, double**& res,
                          casadi_int*& iw, double*& w) const override;

    /// Get all statistics
    Dict get_stats(void* mem) const override;

    // Solve the QP
    int solve(const double** arg, double** res,
      casadi_int* iw, double* w, void* mem) const override;

    /// A documentation string
    static const std::string meta_doc;

    /// All Xpress options
    Dict opts_;

    void serialize_body(SerializingStream &s) const override;

    /** \brief Deserialize with type disambiguation */
    static ProtoFunction* deserialize(DeserializingStream& s) { return new XpressInterface(s); }

    /// Can discrete variables be treated
    bool integer_support() const override { return true; }

    /// SOCP / PSD-formulated cone constraint support
    bool psd_support() const override { return true; }

  protected:
    /** \brief Deserializing constructor */
    explicit XpressInterface(DeserializingStream& s);

  private:

    // Memory structure
    casadi_xpress_prob<double> p_;

    // Linear constraint matrix in CSC form (as int, since Xpress takes int)
    std::vector<int> colinda_, rowa_;

    // Upper-triangular triplet form of the Hessian H_ in column-major order.
    // Entry k corresponds to H_->row()[k_full], H_->colind()[..] elements
    // filtered to row<=col, in the order they appear in H_'s nonzero array.
    // The mapping from H_-nz-index to triplet-index is identity for the
    // upper-triangular nonzeros (we drop the strict lower triangle since
    // CasADi's H_ is typically symmetric and we expect both triangles, or
    // we could expect upper-triangular only; here we accept both and pick
    // upper).  qobj_nz_idx_[k] gives the H_ nz index for the k-th triplet.
    std::vector<int> qobj_col1_, qobj_col2_;
    std::vector<int> qobj_nz_idx_;

    // Per-column type ('I' for integer/discrete, 'C' for continuous).
    // Empty if no discrete variables.
    std::vector<char> coltype_;

    // SOS data, flattened from sos_groups/sos_weights/sos_types options.
    // settype_[k] is '1' or '2'; setstart_[k]..setstart_[k+1] indexes into
    // setind_/refval_ for the k-th group.
    std::vector<char> sos_settype_;
    std::vector<int> sos_setstart_, sos_setind_;
    std::vector<double> sos_refval_;

    // SOCP support: precomputed mapping between CasADi's Q/P input pair and
    // per-cone constraint blocks.  Populated when Q_.nnz() > 0.
    bool has_socp_;
    casadi_socp_prob<double> socp_;
    // Backing storage for socp_ (must outlive socp_)
    std::vector<casadi_int> socp_r_;
    std::vector<casadi_int> socp_mq_colind_, socp_mq_row_, socp_mq_data_;
    std::vector<casadi_int> socp_map_P_;

    // Build socp_ from the result of Conic::sdp_to_socp_init
    void build_socp_config();
  };
} // end namespace casadi
/// \endcond
#endif // CASADI_XPRESS_INTERFACE_HPP

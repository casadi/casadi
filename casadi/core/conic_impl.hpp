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


#ifndef CASADI_CONIC_IMPL_HPP
#define CASADI_CONIC_IMPL_HPP

// blas.hpp must come BEFORE any include that pulls in casadi_runtime.hpp /
// casadi_condensing.hpp.  Phase-1 unqualified lookup of the
// casadi_blas_mtimes calls in condensing.hpp's templated body happens
// when those headers are first parsed; if blas.hpp isn't visible by then,
// the call has no phase-1 candidate, and phase-2 ADL on primitive types
// finds nothing -- compile error at every call site.  (And without the
// declaration, even if the body did parse, conic.cpp would emit a weak
// primary instantiation of casadi_blas_mtimes<double> that races the
// strong spec in blas.cpp -- silently dropping BLAS dispatch.)
#include "blas.hpp"
#include "conic.hpp"
#include "function_internal.hpp"
#include "plugin_interface.hpp"
#include "im.hpp"
// casadi_ocp_block / casadi_condensing_prob / casadi_condensing_data come in
// via casadi_runtime.hpp which is transitively included by im.hpp -> matrix.

/// \cond INTERNAL
namespace casadi {

  struct CASADI_EXPORT ConicMemory : public FunctionMemory {
    // Problem data structure
    casadi_qp_data<double> d_qp;

    // Condensing data struct (populated by Conic::set_work when
    // condense_=True; left null otherwise).  Per-call primal/dual
    // outputs land in d_cond.x / lam_x / lam_a (no separate stash).
    casadi_condensing_data<double> d_cond;
  };

  /// Internal class
  class CASADI_EXPORT Conic : public FunctionInternal, public PluginInterface<Conic> {
  public:
    // Memory structure
    casadi_qp_prob<double> p_qp_;

    // Constructor
    Conic(const std::string& name, const std::map<std::string, Sparsity> &st);

    // Destructor
    ~Conic() override = 0;

    ///@{
    /** \brief Number of function inputs and outputs

        \identifier{1j6} */
    size_t get_n_in() override { return CONIC_NUM_IN;}
    size_t get_n_out() override { return CONIC_NUM_OUT;}
    ///@}

    /// @{
    /** \brief Sparsities of function inputs and outputs

        \identifier{1j7} */
    Sparsity get_sparsity_in(casadi_int i) override;
    Sparsity get_sparsity_out(casadi_int i) override;
    /// @}

    ///@{
    /** \brief Names of function input and outputs

        \identifier{1j8} */
    std::string get_name_in(casadi_int i) override { return conic_in(i);}
    std::string get_name_out(casadi_int i) override { return conic_out(i);}
    /// @}

    ///@{
    /** \brief Options

        \identifier{1j9} */
    static const Options options_;
    const Options& get_options() const override { return options_;}
    ///@}

    /// Solve the QP
    int eval(const double** arg, double** res, casadi_int* iw, double* w, void* mem) const final;

    /// Solve the QP
    virtual int solve(const double** arg, double** res,
      casadi_int* iw, double* w, void* mem) const = 0;

    // Initialize
    void init(const Dict& opts) override;

    // Finalize the object creation
    void finalize() override;

    /** \brief Initalize memory block

        \identifier{1ja} */
    int init_mem(void* mem) const override;

    /** \brief Set the (persistent) work vectors

        \identifier{1jb} */
    void set_work(void* mem, const double**& arg, double**& res,
                          casadi_int*& iw, double*& w) const override;

    /// \brief Check if the numerical values of the supplied bounds make sense
    virtual void check_inputs(const double* lbx, const double* ubx,
                             const double* lba, const double* uba) const;

    /** Generate native code in the interfaced language for debugging */
    virtual void generateNativeCode(std::ostream& file) const;

    // Creator function for internal class
    typedef Conic* (*Creator)(const std::string& name,
                              const std::map<std::string, Sparsity>& st);

    // No static functions exposed
    struct Exposed{ };

    /// Collection of solvers
    static std::map<std::string, Plugin> solvers_;

#ifdef CASADI_WITH_THREADSAFE_SYMBOLICS
    static std::mutex mutex_solvers_;
#endif // CASADI_WITH_THREADSAFE_SYMBOLICS

    /// Infix
    static const std::string infix_;

    /// Short name
    static std::string shortname() { return "conic";}

    /** \brief Check if the function is of a particular type

        \identifier{1jc} */
    bool is_a(const std::string& type, bool recursive) const override;

    /** \brief Get default input value

        \identifier{1jd} */
    double get_default_in(casadi_int ind) const override;

    /// Can discrete variables be treated
    virtual bool integer_support() const { return false;}

    /// Can psd constraints be treated
    virtual bool psd_support() const { return false;}

    /// Get all statistics
    Dict get_stats(void* mem) const override;

    /** \brief Generate code for the function body

        \identifier{24y} */
    void qp_codegen_body(CodeGenerator& g) const;

    /** \brief Emit cleanup/lift after the plugin's solve in codegen */
    void qp_codegen_post(CodeGenerator& g) const;

  protected:
    /// Options
    std::vector<bool> discrete_;
    std::vector<bool> equality_;
    bool print_problem_;
    bool solver_version_check_;
    bool debug_;
    casadi_int condensed_block_count_;
    std::string condense_partition_strategy_;

    /// Partial-condensing wrap (around derived solve)
    bool condense_;
    enum CondenseStructureDetection {
      COND_STRUCT_NONE, COND_STRUCT_AUTO, COND_STRUCT_MANUAL
    };
    CondenseStructureDetection condense_structure_detection_;
    casadi_int N_;            // OCP horizon (original)
    casadi_int N_hat_;        // condensed horizon
    std::vector<casadi_int> nxs_, nus_, ngs_;       // length N+1 (nu[N]=0)
    std::vector<casadi_int> M_;                     // partition, length N_hat+1
    std::vector<casadi_int> M_user_;                // user-supplied M (empty if none)
    std::vector<casadi_int> nx_hat_, nu_hat_, ng_hat_;  // length N_hat+1
    // Block descriptors packed as 4-tuples (offset_r, offset_c, rows, cols)
    std::vector<casadi_int> AB_blocks_packed_;
    std::vector<casadi_int> CD_blocks_packed_;
    std::vector<casadi_int> RSQ_blocks_packed_;
    std::vector<casadi_int> AB_hat_blocks_packed_;
    std::vector<casadi_int> CD_hat_blocks_packed_;
    std::vector<casadi_int> RSQ_hat_blocks_packed_;
    std::vector<casadi_int> AB_offsets_, CD_offsets_, RSQ_offsets_;
    std::vector<casadi_int> AB_hat_offsets_, CD_hat_offsets_, RSQ_hat_offsets_;
    Sparsity RSQsp_, ABsp_, CDsp_;
    Sparsity H_hat_sp_, A_hat_sp_;
    casadi_int nx_total_hat_, na_total_hat_;
    casadi_int nx_max_, nu_max_block_, nxu_max_block_;

    // Native casadi_ocp_block storage (unpacked from *_blocks_packed_)
    // and the condensing prob struct that points into our vectors.
    // These are filled in build_condense_blocks() and
    // finalize_condense_prob() (called from init / deserialize).
    mutable std::vector<casadi_ocp_block> AB_blocks_;
    mutable std::vector<casadi_ocp_block> CD_blocks_;
    mutable std::vector<casadi_ocp_block> RSQ_blocks_;
    mutable std::vector<casadi_ocp_block> AB_hat_blocks_;
    mutable std::vector<casadi_ocp_block> CD_hat_blocks_;
    mutable std::vector<casadi_ocp_block> RSQ_hat_blocks_;
    mutable casadi_condensing_prob<double> p_cond_;
    void finalize_condense_prob();
    /// Unpack a flat [off_r, off_c, rows, cols]*N vector into ocp_block records
    static void unpack_ocp_blocks(const std::vector<casadi_int>& packed,
                                  std::vector<casadi_ocp_block>& dst);
    /// Prepend a length scalar to a vector (used to length-prefix codegen consts)
    static std::vector<casadi_int> len_prefixed(casadi_int n,
                                                const std::vector<casadi_int>& v);

    /// Run structure detection (auto) or use provided (manual)
    void detect_condense_structure();
    /// Build packed block descriptors and condensed sparsities
    void build_condense_blocks();

    /// Pick M_ from condense_partition_strategy_ + condensed_block_count_.
    /// Called from build_condense_blocks() when M_ is empty.
    void derive_condense_partition();
    /// Frison closed-form: pick n_hat for uniform problems if 0,
    /// then emit a uniform-length partition into M_.
    void frison_uniform_partition();
    /// O(K_max * N^2) DP picking M to minimize sum_K (nx[M_K] +
    /// sum nu[j])^3.  If condensed_block_count_ == 0, scans K=1..K_max
    /// and picks the argmin; otherwise fixes K = condensed_block_count_.
    void dp_optimal_partition();

    /// Problem structure
    Sparsity H_, A_, Q_, P_;

    /// Number of decision variables
    casadi_int nx_;

    /// The number of constraints (counting both equality and inequality) == A.size1()
    casadi_int na_;

    /// The shape of psd constraint matrix
    casadi_int np_;

    /// SDP to SOCP conversion memory
    struct SDPToSOCPMem {
      // Block partition vector for SOCP (block i runs from r[i] to r[i+1])
      std::vector<casadi_int> r;

      // Tranpose of A, and corresponding mapping
      Sparsity AT;
      std::vector<casadi_int> A_mapping;

      // Aggregate SOCP helper constraints (lhs)
      IM map_Q;

      // Aggregate SOCP helper constraints (rhs)
      std::vector<casadi_int> map_P;

      // Maximum size of ind/val vectors
      casadi_int indval_size;
    };

    /// SDP to SOCP conversion initialization
    void sdp_to_socp_init(SDPToSOCPMem& mem) const;

    void serialize(SerializingStream &s, const SDPToSOCPMem& m) const;
    void deserialize(DeserializingStream &s, SDPToSOCPMem& m);

  public:
      /** \brief Serialize an object without type information

          \identifier{1je} */
    void serialize_body(SerializingStream &s) const override;
    /** \brief Serialize type information

        \identifier{1jf} */
    void serialize_type(SerializingStream &s) const override;

    /** \brief String used to identify the immediate FunctionInternal subclass

        \identifier{1jg} */
    std::string serialize_base_function() const override { return "Conic"; }
    /** \brief Deserialize with type disambiguation

        \identifier{1jh} */
    static ProtoFunction* deserialize(DeserializingStream& s);

  protected:

    /** \brief Deserializing constructor

        \identifier{1ji} */
    explicit Conic(DeserializingStream& s);
  private:
    void set_qp_prob();
  };


} // namespace casadi
/// \endcond
#endif // CASADI_CONIC_IMPL_HPP

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


#ifndef CASADI_MULTIPLICATION_HPP
#define CASADI_MULTIPLICATION_HPP

#include "mx_node.hpp"

/// \cond INTERNAL

namespace casadi {
  /** \brief An MX atomic for matrix-matrix product,

             note that the first factor must be provided transposed
      \author Joel Andersson
      \date 2010

      The base implementation handles arbitrary sparsity via casadi_mtimes.
      Subclasses specialize the kernel for structured operands (all-dense,
      dense * sparse, compactible / pseudo-dense). Each subclass overrides
      eval_kernel + generate + serialize_type; everything else is shared.

      \identifier{11h} */
  class CASADI_EXPORT Multiplication : public MXNode {
  public:

    /** \brief  Factory: dispatch to the most specific subclass for the given operands

        Detection ladder is most-specific-first; each branch is at most O(nnz) on
        the operand sparsities. Callers should use this rather than picking a
        subclass directly.

        `blas` selects the dense matrix-multiply backend (e.g. "reference",
        "classic", "blasfeo"). It is honored only by the dense and pseudo-
        dense branches; sparse-only nodes ignore it but preserve the choice
        across serialization. */
    static MX create(const MX& z, const MX& x, const MX& y,
                     const std::string& blas = "reference");

    /** \brief  Constructor

        \identifier{11i} */
    Multiplication(const MX& z, const MX& x, const MX& y,
                   const std::string& blas = "reference");

    /** \brief  Destructor

        \identifier{11j} */
    ~Multiplication() override {}

    /** \brief  Print expression

        \identifier{11k} */
    std::string disp(const std::vector<std::string>& arg) const override;

    /** \brief Generate code for the operation

        \identifier{11l} */
    void generate(CodeGenerator& g,
                  const std::vector<casadi_int>& arg,
                  const std::vector<casadi_int>& res,
                  const std::vector<bool>& arg_is_ref,
                  std::vector<bool>& res_is_ref) const override;

    /** \brief Subclass hook: mathematical kernel z += x*y on the input buffers

        Default implementation calls casadi_mtimes (general sparse).
        Two non-template overloads avoid virtual templates while still
        sharing the eval_gen wrapper across types. */
    virtual void eval_kernel(const double** arg, double** res, double* w) const;
    virtual void eval_kernel(const SXElem** arg, SXElem** res, SXElem* w) const;

    /// Evaluate the function (template) — copy z, then run subclass kernel
    template<typename T>
    int eval_gen(const T** arg, T** res, casadi_int* iw, T* w) const {
      if (arg[0]!=res[0]) std::copy(arg[0], arg[0]+dep(0).nnz(), res[0]);
      eval_kernel(arg, res, w);
      return 0;
    }

    /// Evaluate the function numerically
    int eval(const double** arg, double** res, casadi_int* iw, double* w) const override {
      return eval_gen<double>(arg, res, iw, w);
    }

    /// Evaluate the function symbolically (SX)
    int eval_sx(const SXElem** arg, SXElem** res, casadi_int* iw, SXElem* w) const override {
      return eval_gen<SXElem>(arg, res, iw, w);
    }

    /** \brief  Evaluate symbolically (MX)

        \identifier{11m} */
    void eval_mx(const std::vector<MX>& arg, std::vector<MX>& res,
        const std::vector<bool>& unique={}) const override;

    /** \brief Evaluate the MX node on a const/linear/nonlinear partition

        \identifier{28a} */
    void eval_linear(const std::vector<std::array<MX, 3> >& arg,
                        std::vector<std::array<MX, 3> >& res) const override;

    /** \brief Calculate forward mode directional derivatives

        \identifier{11n} */
    void ad_forward(const std::vector<std::vector<MX> >& fseed,
                         std::vector<std::vector<MX> >& fsens) const override;

    /** \brief Calculate reverse mode directional derivatives

        \identifier{11o} */
    void ad_reverse(const std::vector<std::vector<MX> >& aseed,
                         std::vector<std::vector<MX> >& asens) const override;

    /** \brief  Propagate sparsity forward

        \identifier{11p} */
    int sp_forward(const bvec_t** arg, bvec_t** res, casadi_int* iw, bvec_t* w) const override;

    /** \brief  Propagate sparsity backwards

        \identifier{11q} */
    int sp_reverse(bvec_t** arg, bvec_t** res, casadi_int* iw, bvec_t* w) const override;

    /** \brief Get the operation

        \identifier{11r} */
    casadi_int op() const override { return OP_MTIMES;}

    /// Can the operation be performed inplace (i.e. overwrite the result)
    casadi_int n_inplace() const override { return 1;}

    /** \brief Check if two nodes are equivalent up to a given depth

        \identifier{11s} */
    bool is_equal(const MXNode* node, casadi_int depth) const override {
      return sameOpAndDeps(node, depth) && dynamic_cast<const Multiplication*>(node)!=nullptr;
    }

    /** \brief Get required length of w field

        \identifier{11t} */
    size_t sz_w() const override { return sparsity().size1();}

    /** \brief Serialize specific part of node

        \identifier{11u} */
    void serialize_type(SerializingStream& s) const override;

    /** \brief Serialize body (BLAS plugin name)

        Stored on the base so all subclasses round-trip the choice,
        even those that don't currently consume it. */
    void serialize_body(SerializingStream& s) const override;

    /** \brief Deserialize with type disambiguation

        \identifier{11v} */
    static MXNode* deserialize(DeserializingStream& s);

  protected:
    /** \brief Deserializing constructor.

        Reads MXNode state, then unpacks the BLAS plugin name into
        blas_shorthand_. legacy=true is used by the pre-3.8 wire-format
        path (no blas field on the wire) and forces shorthand 0 ("reference").

        \identifier{11w} */
    explicit Multiplication(DeserializingStream& s, bool legacy = false);

    /// Cached BLAS plugin shorthand. 0 == "reference". Derived state, never
    /// serialized directly; the underlying name string is what crosses
    /// serialization boundaries.
    casadi_int blas_shorthand_;
  };


  /** \brief Dense * Dense -> Dense matrix product

      \identifier{11x} */
  class CASADI_EXPORT DenseMultiplication : public Multiplication{
  public:
    /// Returns a fresh node iff x, y, z are all dense; otherwise nullptr.
    static MXNode* try_create(const MX& z, const MX& x, const MX& y,
                              const std::string& blas = "reference");

    DenseMultiplication(const MX& z, const MX& x, const MX& y,
                        const std::string& blas = "reference")
        : Multiplication(z, x, y, blas) {}
    ~DenseMultiplication() override {}

    void eval_kernel(const double** arg, double** res, double* w) const override;
    void eval_kernel(const SXElem** arg, SXElem** res, SXElem* w) const override;

    void generate(CodeGenerator& g,
                  const std::vector<casadi_int>& arg,
                  const std::vector<casadi_int>& res,
                  const std::vector<bool>& arg_is_ref,
                  std::vector<bool>& res_is_ref) const override;

    void serialize_type(SerializingStream& s) const override;
    explicit DenseMultiplication(DeserializingStream& s, bool legacy = false)
        : Multiplication(s, legacy) {}
  };


  /** \brief Dense * Sparse -> Dense matrix product

      Iterate columns of y; per column, walk y's nonzeros and AXPY a column of x into z. */
  class CASADI_EXPORT DenseSparseMultiplication : public Multiplication {
  public:
    /// Returns a fresh node iff x, z are dense; otherwise nullptr.
    static MXNode* try_create(const MX& z, const MX& x, const MX& y,
                              const std::string& blas = "reference");

    DenseSparseMultiplication(const MX& z, const MX& x, const MX& y,
                              const std::string& blas = "reference")
        : Multiplication(z, x, y, blas) {}
    ~DenseSparseMultiplication() override {}

    void eval_kernel(const double** arg, double** res, double* w) const override;
    void eval_kernel(const SXElem** arg, SXElem** res, SXElem* w) const override;

    void generate(CodeGenerator& g,
                  const std::vector<casadi_int>& arg,
                  const std::vector<casadi_int>& res,
                  const std::vector<bool>& arg_is_ref,
                  std::vector<bool>& res_is_ref) const override;

    void serialize_type(SerializingStream& s) const override;
    explicit DenseSparseMultiplication(DeserializingStream& s) : Multiplication(s) {}
  };


  /** \brief Compactible * Compactible -> Compactible product

      When the nonzero patterns of x, y, z are each Cartesian products
      whose connecting index sets agree (col(x)==row(y), row(x)==row(z),
      col(y)==col(z)), the CCS nonzero buffers are already laid out as
      column-major dense matrices of compact dimensions. We can call
      casadi_mtimes_dense directly on the nz buffers — no gather/scatter.

      The compact dimensions are stored at construction time. */
  class CASADI_EXPORT PseudoDenseMultiplication : public Multiplication {
  public:
    /// Returns a fresh node iff x, y, z are all compactible with matching
    /// connecting index sets; otherwise nullptr.
    static MXNode* try_create(const MX& z, const MX& x, const MX& y,
                              const std::string& blas = "reference");

    PseudoDenseMultiplication(const MX& z, const MX& x, const MX& y,
                              casadi_int nrow_x_compact, casadi_int ncol_x_compact,
                              casadi_int ncol_y_compact,
                              const std::string& blas = "reference");
    ~PseudoDenseMultiplication() override {}

    void eval_kernel(const double** arg, double** res, double* w) const override;
    void eval_kernel(const SXElem** arg, SXElem** res, SXElem* w) const override;

    void generate(CodeGenerator& g,
                  const std::vector<casadi_int>& arg,
                  const std::vector<casadi_int>& res,
                  const std::vector<bool>& arg_is_ref,
                  std::vector<bool>& res_is_ref) const override;

    void serialize_type(SerializingStream& s) const override;
    void serialize_body(SerializingStream& s) const override;
    explicit PseudoDenseMultiplication(DeserializingStream& s);

  private:
    // Compact (dense) dimensions seen by casadi_mtimes_dense
    casadi_int a_, b_, c_;  // x: a x b, y: b x c, z: a x c
  };


} // namespace casadi
/// \endcond

#endif // CASADI_MULTIPLICATION_HPP

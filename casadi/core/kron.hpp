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


#ifndef CASADI_KRON_HPP
#define CASADI_KRON_HPP

#include "mx_node.hpp"

/// \cond INTERNAL

namespace casadi {

  /** \brief Kronecker product

      F = A \otimes B, shape (mA*mB, nA*nB), with
      F[i*mB+r, j*nB+s] = A[i,j] * B[r,s].

      The base implementation handles arbitrary sparsity via casadi_kron.
      Subclasses specialize the kernel for dense / dense-sparse / sparse-dense
      operands. Each subclass overrides eval_kernel + generate + serialize_type;
      AD, sparsity propagation, and everything else is shared.

      \author Joris Gillis
      \date 2026

      \identifier{2gn} */
  class CASADI_EXPORT Kron : public MXNode {
  public:

    /// Factory: dispatch to the most specific subclass for the given operands
    static MX create(const MX& a, const MX& b);

    /// Constructor
    Kron(const MX& a, const MX& b);

    /// Destructor
    ~Kron() override {}

    /// Subclass hook: r = kron(a, b) on the input buffers
    virtual void eval_kernel(const double** arg, double** res) const;
    virtual void eval_kernel(const SXElem** arg, SXElem** res) const;

    /// Evaluate the function (template) — dispatches to eval_kernel
    template<typename T>
    int eval_gen(const T** arg, T** res, casadi_int* iw, T* w) const {
      eval_kernel(arg, res);
      return 0;
    }

    /** \brief  Evaluate numerically

        \identifier{2go} */
    int eval(const double** arg, double** res, casadi_int* iw, double* w) const override {
      return eval_gen<double>(arg, res, iw, w);
    }

    /** \brief  Evaluate symbolically (SX)

        \identifier{2gp} */
    int eval_sx(const SXElem** arg, SXElem** res, casadi_int* iw, SXElem* w) const override {
      return eval_gen<SXElem>(arg, res, iw, w);
    }

    /** \brief  Evaluate symbolically (MX)

        \identifier{2gq} */
    void eval_mx(const std::vector<MX>& arg, std::vector<MX>& res,
        const std::vector<bool>& unique={}) const override;

    /** \brief  Propagate sparsity forward

        \identifier{2gr} */
    int sp_forward(const bvec_t** arg, bvec_t** res, casadi_int* iw, bvec_t* w) const override;

    /** \brief  Propagate sparsity backwards

        \identifier{2gs} */
    int sp_reverse(bvec_t** arg, bvec_t** res, casadi_int* iw, bvec_t* w) const override;

    /** \brief Calculate forward mode directional derivatives

        \identifier{2gt} */
    void ad_forward(const std::vector<std::vector<MX> >& fseed,
                    std::vector<std::vector<MX> >& fsens) const override;

    /** \brief Calculate reverse mode directional derivatives

        \identifier{2gu} */
    void ad_reverse(const std::vector<std::vector<MX> >& aseed,
                    std::vector<std::vector<MX> >& asens) const override;

    /** \brief Generate code for the operation

        \identifier{2gv} */
    void generate(CodeGenerator& g,
                  const std::vector<casadi_int>& arg,
                  const std::vector<casadi_int>& res,
                  const std::vector<bool>& arg_is_ref,
                  std::vector<bool>& res_is_ref) const override;

    /** \brief  Print expression

        \identifier{2gw} */
    std::string disp(const std::vector<std::string>& arg) const override;

    /** \brief Get the operation

        \identifier{2gx} */
    casadi_int op() const override { return OP_KRON;}

    /** \brief Serialize specific part of node

        \identifier{2gy} */
    void serialize_type(SerializingStream& s) const override;

    /** \brief Deserialize with type disambiguation

        \identifier{2gz} */
    static MXNode* deserialize(DeserializingStream& s);

  protected:

    /** \brief Deserializing constructor

        \identifier{2h0} */
    explicit Kron(DeserializingStream& s) : MXNode(s) {}
  };


  /** \brief Kron specialization: both operands dense

      \identifier{2h1} */
  class CASADI_EXPORT DenseKron : public Kron {
  public:
    static MXNode* try_create(const MX& a, const MX& b);

    DenseKron(const MX& a, const MX& b) : Kron(a, b) {}
    ~DenseKron() override {}

    void eval_kernel(const double** arg, double** res) const override;
    void eval_kernel(const SXElem** arg, SXElem** res) const override;

    void generate(CodeGenerator& g,
                  const std::vector<casadi_int>& arg,
                  const std::vector<casadi_int>& res,
                  const std::vector<bool>& arg_is_ref,
                  std::vector<bool>& res_is_ref) const override;

    void serialize_type(SerializingStream& s) const override;
    explicit DenseKron(DeserializingStream& s) : Kron(s) {}
  };


  /** \brief Kron specialization: dense a + sparse b

      \identifier{2h2} */
  class CASADI_EXPORT DenseSparseKron : public Kron {
  public:
    static MXNode* try_create(const MX& a, const MX& b);

    DenseSparseKron(const MX& a, const MX& b) : Kron(a, b) {}
    ~DenseSparseKron() override {}

    void eval_kernel(const double** arg, double** res) const override;
    void eval_kernel(const SXElem** arg, SXElem** res) const override;

    void generate(CodeGenerator& g,
                  const std::vector<casadi_int>& arg,
                  const std::vector<casadi_int>& res,
                  const std::vector<bool>& arg_is_ref,
                  std::vector<bool>& res_is_ref) const override;

    void serialize_type(SerializingStream& s) const override;
    explicit DenseSparseKron(DeserializingStream& s) : Kron(s) {}
  };


  /** \brief Kron specialization: sparse a + dense b

      \identifier{2h3} */
  class CASADI_EXPORT SparseDenseKron : public Kron {
  public:
    static MXNode* try_create(const MX& a, const MX& b);

    SparseDenseKron(const MX& a, const MX& b) : Kron(a, b) {}
    ~SparseDenseKron() override {}

    void eval_kernel(const double** arg, double** res) const override;
    void eval_kernel(const SXElem** arg, SXElem** res) const override;

    void generate(CodeGenerator& g,
                  const std::vector<casadi_int>& arg,
                  const std::vector<casadi_int>& res,
                  const std::vector<bool>& arg_is_ref,
                  std::vector<bool>& res_is_ref) const override;

    void serialize_type(SerializingStream& s) const override;
    explicit SparseDenseKron(DeserializingStream& s) : Kron(s) {}
  };


  /** \brief Kronecker contraction

      Viewing M (shape mA*mB x nA*nB) as a 4-tensor M[i, r, j, s], this is a
      double contraction with a 2-tensor X. Two modes:

      - `inner = true`:  Y[i, j] = sum over (r, s) of M[i*mB+r, j*nB+s] * X[r, s]
                         (X shape mB x nB; Y shape mA x nA)
      - `inner = false`: Y[r, s] = sum over (i, j) of X[i, j] * M[i*mB+r, j*nB+s]
                         (X shape mA x nA; Y shape mB x nB)

      Together with Kron, this closes the AD algebra: each AD rule emits a
      constant number of nodes regardless of operand dimensions.

      \author Joris Gillis
      \date 2026

      \identifier{2h4} */
  class CASADI_EXPORT KronContract : public MXNode {
  public:

    /// Factory: dispatch to the most specific subclass for the given operands
    static MX create(const MX& m, const MX& x, bool inner);

    /// Constructor (auto-derives output sparsity)
    KronContract(const MX& m, const MX& x, bool inner);

    /// Destructor
    ~KronContract() override {}

    /// Subclass hook: y = kron_contract(m, x, inner_) on the input buffers
    /// w is scratch (sz_w bytes worth of T) for use by the kernel.
    virtual void eval_kernel(const double** arg, double** res, double* w) const;
    virtual void eval_kernel(const SXElem** arg, SXElem** res, SXElem* w) const;

    /// Evaluate the function (template) — dispatches to eval_kernel
    template<typename T>
    int eval_gen(const T** arg, T** res, casadi_int* iw, T* w) const {
      eval_kernel(arg, res, w);
      return 0;
    }

    /** \brief  Evaluate numerically

        \identifier{2h5} */
    int eval(const double** arg, double** res, casadi_int* iw, double* w) const override {
      return eval_gen<double>(arg, res, iw, w);
    }

    /** \brief  Evaluate symbolically (SX)

        \identifier{2h6} */
    int eval_sx(const SXElem** arg, SXElem** res, casadi_int* iw, SXElem* w) const override {
      return eval_gen<SXElem>(arg, res, iw, w);
    }

    /** \brief  Evaluate symbolically (MX)

        \identifier{2h7} */
    void eval_mx(const std::vector<MX>& arg, std::vector<MX>& res,
        const std::vector<bool>& unique={}) const override;

    /** \brief  Propagate sparsity forward

        \identifier{2h8} */
    int sp_forward(const bvec_t** arg, bvec_t** res, casadi_int* iw, bvec_t* w) const override;

    /** \brief  Propagate sparsity backwards

        \identifier{2h9} */
    int sp_reverse(bvec_t** arg, bvec_t** res, casadi_int* iw, bvec_t* w) const override;

    /** \brief Calculate forward mode directional derivatives

        \identifier{2ha} */
    void ad_forward(const std::vector<std::vector<MX> >& fseed,
                    std::vector<std::vector<MX> >& fsens) const override;

    /** \brief Calculate reverse mode directional derivatives

        \identifier{2hb} */
    void ad_reverse(const std::vector<std::vector<MX> >& aseed,
                    std::vector<std::vector<MX> >& asens) const override;

    /** \brief Get required length of w field

        \identifier{2hc} */
    size_t sz_w() const override;

    /** \brief Generate code for the operation

        \identifier{2hd} */
    void generate(CodeGenerator& g,
                  const std::vector<casadi_int>& arg,
                  const std::vector<casadi_int>& res,
                  const std::vector<bool>& arg_is_ref,
                  std::vector<bool>& res_is_ref) const override;

    /** \brief  Print expression

        \identifier{2he} */
    std::string disp(const std::vector<std::string>& arg) const override;

    /** \brief Get the operation

        \identifier{2hf} */
    casadi_int op() const override { return OP_KRON_CONTRACT;}

    /** \brief Serialize an object without type information

        \identifier{2hg} */
    void serialize_body(SerializingStream& s) const override;

    /** \brief Serialize specific part of node

        \identifier{2hh} */
    void serialize_type(SerializingStream& s) const override;

    /** \brief Deserialize with type disambiguation

        \identifier{2hi} */
    static MXNode* deserialize(DeserializingStream& s);

    /// Which axes are contracted (true = inner (mB, nB); false = outer (mA, nA))
    bool inner_;

  protected:

    /** \brief Deserializing constructor

        \identifier{2hj} */
    explicit KronContract(DeserializingStream& s);
  };


  /** \brief KronContract specialization: M dense, X dense (=> Y dense)

      \identifier{2hk} */
  class CASADI_EXPORT DenseKronContract : public KronContract {
  public:
    static MXNode* try_create(const MX& m, const MX& x, bool inner);

    DenseKronContract(const MX& m, const MX& x, bool inner) : KronContract(m, x, inner) {}
    ~DenseKronContract() override {}

    void eval_kernel(const double** arg, double** res, double* w) const override;
    void eval_kernel(const SXElem** arg, SXElem** res, SXElem* w) const override;

    // sz_w stays at base size (dep(1).numel()) because sp_forward uses it
    // even though our specialized eval_kernel doesn't.

    void generate(CodeGenerator& g,
                  const std::vector<casadi_int>& arg,
                  const std::vector<casadi_int>& res,
                  const std::vector<bool>& arg_is_ref,
                  std::vector<bool>& res_is_ref) const override;

    void serialize_type(SerializingStream& s) const override;
    explicit DenseKronContract(DeserializingStream& s) : KronContract(s) {}
  };


  /** \brief KronContract specialization: M dense, X sparse (=> Y dense)

      \identifier{2hl} */
  class CASADI_EXPORT DenseSparseKronContract : public KronContract {
  public:
    static MXNode* try_create(const MX& m, const MX& x, bool inner);

    DenseSparseKronContract(const MX& m, const MX& x, bool inner) : KronContract(m, x, inner) {}
    ~DenseSparseKronContract() override {}

    void eval_kernel(const double** arg, double** res, double* w) const override;
    void eval_kernel(const SXElem** arg, SXElem** res, SXElem* w) const override;

    // sz_w stays at base size for sp_forward's densify step.

    void generate(CodeGenerator& g,
                  const std::vector<casadi_int>& arg,
                  const std::vector<casadi_int>& res,
                  const std::vector<bool>& arg_is_ref,
                  std::vector<bool>& res_is_ref) const override;

    void serialize_type(SerializingStream& s) const override;
    explicit DenseSparseKronContract(DeserializingStream& s) : KronContract(s) {}
  };


  /** \brief KronContract specialization: M sparse, X dense

      \identifier{2hm} */
  class CASADI_EXPORT SparseDenseKronContract : public KronContract {
  public:
    static MXNode* try_create(const MX& m, const MX& x, bool inner);

    SparseDenseKronContract(const MX& m, const MX& x, bool inner) : KronContract(m, x, inner) {}
    ~SparseDenseKronContract() override {}

    void eval_kernel(const double** arg, double** res, double* w) const override;
    void eval_kernel(const SXElem** arg, SXElem** res, SXElem* w) const override;

    // sz_w stays at base size for sp_forward's densify step.

    void generate(CodeGenerator& g,
                  const std::vector<casadi_int>& arg,
                  const std::vector<casadi_int>& res,
                  const std::vector<bool>& arg_is_ref,
                  std::vector<bool>& res_is_ref) const override;

    void serialize_type(SerializingStream& s) const override;
    explicit SparseDenseKronContract(DeserializingStream& s) : KronContract(s) {}
  };

} // namespace casadi
/// \endcond

#endif // CASADI_KRON_HPP

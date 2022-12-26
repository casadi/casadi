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


#ifndef CASADI_MULTIPLICATION_HPP
#define CASADI_MULTIPLICATION_HPP

#include "mx_node.hpp"

/// \cond INTERNAL

namespace casadi {
  /** \brief An MX atomic for matrix-matrix product,

             note that the first factor must be provided transposed
      \author Joel Andersson
      \date 2010

      \identifier{11h} */
  class CASADI_EXPORT Multiplication : public MXNode {
  public:

    /** \brief  Constructor

        \identifier{11i} */
    Multiplication(const MX& z, const MX& x, const MX& y, const Dict& opts=Dict());

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
                  const std::vector<casadi_int>& res, bool prefer_inline=false) const override;

    /// Evaluate the function (template)
    template<typename T>
    int eval_gen(const T** arg, T** res, casadi_int* iw, T* w) const;

    /// Evaluate the function numerically
    int eval(const double** arg, double** res, casadi_int* iw, double* w) const override;

    /// Evaluate the function symbolically (SX)
    int eval_sx(const SXElem** arg, SXElem** res, casadi_int* iw, SXElem* w) const override;

    /** \brief  Evaluate symbolically (MX)

        \identifier{11m} */
    void eval_mx(const std::vector<MX>& arg, std::vector<MX>& res) const override;

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
    casadi_int n_inplace() const override { return n_dep()==3;}

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

    /** \brief Deserialize with type disambiguation

        \identifier{11v} */
    static MXNode* deserialize(DeserializingStream& s);

  protected:
    /** \brief Deserializing constructor

        \identifier{11w} */
    explicit Multiplication(DeserializingStream& s) : MXNode(s) {}
  };


  /** \brief An MX atomic for matrix-matrix product,

             note that the factor must be provided transposed
      \author Joel Andersson
      \date 2010

      \identifier{11x} */
  class CASADI_EXPORT DenseMultiplication : public Multiplication{
  public:

    /** \brief  Constructor

        \identifier{11y} */
    DenseMultiplication(const MX& z, const MX& x, const MX& y, const Dict& opts=Dict())
        : Multiplication(z, x, y, opts) {}

    /** \brief  Destructor

        \identifier{11z} */
    ~DenseMultiplication() override {}

    /** \brief Generate code for the operation

        \identifier{120} */
    void generate(CodeGenerator& g,
                  const std::vector<casadi_int>& arg,
                  const std::vector<casadi_int>& res, bool prefer_inline=false) const override;

    /** \brief Serialize specific part of node

        \identifier{121} */
    void serialize_type(SerializingStream& s) const override;

    /** \brief Deserializing constructor

        \identifier{122} */
    explicit DenseMultiplication(DeserializingStream& s) : Multiplication(s) {}
  };


} // namespace casadi
/// \endcond

#endif // CASADI_MULTIPLICATION_HPP

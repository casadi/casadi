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


#ifndef CASADI_UNARY_MX_HPP
#define CASADI_UNARY_MX_HPP

#include "mx_node.hpp"
/// \cond INTERNAL

namespace casadi {
  /** \brief Represents a general unary operation on an MX

      \author Joel Andersson
      \date 2010

      \identifier{176} */
  class CASADI_EXPORT UnaryMX : public MXNode {
  public:

    /** \brief  Constructor is private, use "create" below

        \identifier{177} */
    UnaryMX(Operation op, MX x);

    /** \brief  Destructor

        \identifier{178} */
    ~UnaryMX() override {}

    /** \brief  Print expression

        \identifier{179} */
    std::string disp(const std::vector<std::string>& arg) const override;

    /// Evaluate the function numerically
    int eval(const double** arg, double** res, casadi_int* iw, double* w) const override;

    /// Evaluate the function symbolically (SX)
    int eval_sx(const SXElem** arg, SXElem** res, casadi_int* iw, SXElem* w) const override;

    /** \brief  Evaluate symbolically (MX)

        \identifier{17a} */
    void eval_mx(const std::vector<MX>& arg, std::vector<MX>& res) const override;

    /** \brief Evaluate the MX node on a const/linear/nonlinear partition

        \identifier{28d} */
    void eval_linear(const std::vector<std::array<MX, 3> >& arg,
                        std::vector<std::array<MX, 3> >& res) const override;

    /** \brief Calculate forward mode directional derivatives

        \identifier{17b} */
    void ad_forward(const std::vector<std::vector<MX> >& fseed,
                         std::vector<std::vector<MX> >& fsens) const override;

    /** \brief Calculate reverse mode directional derivatives

        \identifier{17c} */
    void ad_reverse(const std::vector<std::vector<MX> >& aseed,
                         std::vector<std::vector<MX> >& asens) const override;

    /** \brief  Propagate sparsity forward

        \identifier{17d} */
    int sp_forward(const bvec_t** arg, bvec_t** res, casadi_int* iw, bvec_t* w) const override;

    /** \brief  Propagate sparsity backwards

        \identifier{17e} */
    int sp_reverse(bvec_t** arg, bvec_t** res, casadi_int* iw, bvec_t* w) const override;

    /** \brief Check if unary operation

        \identifier{17f} */
    bool is_unary() const override { return true;}

    /** \brief Get the operation

        \identifier{17g} */
    casadi_int op() const override { return op_;}

    /** \brief Generate code for the operation

        \identifier{17h} */
    void generate(CodeGenerator& g,
                  const std::vector<casadi_int>& arg,
                  const std::vector<casadi_int>& res) const override;

    /// Can the operation be performed inplace (i.e. overwrite the result)
    casadi_int n_inplace() const override { return 1;}

    /// Get a unary operation
    MX get_unary(casadi_int op) const override;

    /// Get a binary operation operation
    MX _get_binary(casadi_int op, const MX& y, bool scX, bool scY) const override;

    /** \brief Check if two nodes are equivalent up to a given depth

        \identifier{17i} */
    bool is_equal(const MXNode* node, casadi_int depth) const override {
      return sameOpAndDeps(node, depth);
    }

    /** \brief Serialize an object without type information

        \identifier{17j} */
    void serialize_body(SerializingStream& s) const override;

    /** \brief Deserialize without type information

        \identifier{17k} */
    static MXNode* deserialize(DeserializingStream& s) { return new UnaryMX(s); }

    //! \brief operation
    Operation op_;

  protected:
    /** \brief Deserializing constructor

        \identifier{17l} */
    explicit UnaryMX(DeserializingStream& s);

  };

} // namespace casadi

/// \endcond

#endif // CASADI_UNARY_MX_HPP

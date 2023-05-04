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


#ifndef CASADI_RESHAPE_HPP
#define CASADI_RESHAPE_HPP

#include "mx_node.hpp"
#include <map>
#include <stack>

/// \cond INTERNAL

namespace casadi {
  /** \brief Reshape an expression

      \author Joel Andersson
      \date 2013

      \identifier{1dj} */
  class CASADI_EXPORT Reshape : public MXNode {
  public:

    /// Constructor
    Reshape(const MX& x, Sparsity sp);

    /// Destructor
    ~Reshape() override {}

    /// Evaluate the function (template)
    template<typename T>
    int eval_gen(const T** arg, T** res, casadi_int* iw, T* w) const;

    /// Evaluate the function numerically
    int eval(const double** arg, double** res, casadi_int* iw, double* w) const override;

    /// Evaluate the function symbolically (SX)
    int eval_sx(const SXElem** arg, SXElem** res, casadi_int* iw, SXElem* w) const override;

    /** \brief  Evaluate symbolically (MX)

        \identifier{1dk} */
    void eval_mx(const std::vector<MX>& arg, std::vector<MX>& res) const override;

    /** \brief Calculate forward mode directional derivatives

        \identifier{1dl} */
    void ad_forward(const std::vector<std::vector<MX> >& fseed,
                         std::vector<std::vector<MX> >& fsens) const override;

    /** \brief Calculate reverse mode directional derivatives

        \identifier{1dm} */
    void ad_reverse(const std::vector<std::vector<MX> >& aseed,
                         std::vector<std::vector<MX> >& asens) const override;

    /** \brief  Propagate sparsity forward

        \identifier{1dn} */
    int sp_forward(const bvec_t** arg, bvec_t** res, casadi_int* iw, bvec_t* w) const override;

    /** \brief  Propagate sparsity backwards

        \identifier{1do} */
    int sp_reverse(bvec_t** arg, bvec_t** res, casadi_int* iw, bvec_t* w) const override;

    /** \brief  Print expression

        \identifier{1dp} */
    std::string disp(const std::vector<std::string>& arg) const override;

    /** \brief Generate code for the operation

        \identifier{1dq} */
    void generate(CodeGenerator& g,
                  const std::vector<casadi_int>& arg,
                  const std::vector<casadi_int>& res) const override;

    /** \brief Get the operation

        \identifier{1dr} */
    casadi_int op() const override { return OP_RESHAPE;}

    /// Can the operation be performed inplace (i.e. overwrite the result)
    casadi_int n_inplace() const override { return 1;}

    /// Reshape
    MX get_reshape(const Sparsity& sp) const override;

    /** \brief Check if two nodes are equivalent up to a given depth

        \identifier{1ds} */
    bool is_equal(const MXNode* node, casadi_int depth) const override
    { return sameOpAndDeps(node, depth) && sparsity()==node->sparsity();}

    /// Transpose (if a dimension is one)
    MX get_transpose() const override;

    /** \brief  Check if valid function input

        \identifier{1dt} */
    bool is_valid_input() const override;

    /** \brief Get the number of symbolic primitives

        \identifier{1du} */
    casadi_int n_primitives() const override;

    /** \brief Get symbolic primitives

        \identifier{1dv} */
    void primitives(std::vector<MX>::iterator& it) const override;

    /** \brief Split up an expression along symbolic primitives

        \identifier{1dw} */
    void split_primitives(const MX& x, std::vector<MX>::iterator& it) const override;

    /** \brief Join an expression along symbolic primitives

        \identifier{1dx} */
    MX join_primitives(std::vector<MX>::const_iterator& it) const override;

    /** \brief Detect duplicate symbolic expressions

        \identifier{1dy} */
    bool has_duplicates() const override;

    /** \brief Reset the marker for an input expression

        \identifier{1dz} */
    void reset_input() const override;

    /** \brief Deserialize without type information

        \identifier{1e0} */
    static MXNode* deserialize(DeserializingStream& s) { return new Reshape(s); }
  protected:
    /** \brief Deserializing constructor

        \identifier{1e1} */
    explicit Reshape(DeserializingStream& s) : MXNode(s) {}
  };

} // namespace casadi
/// \endcond

#endif // CASADI_RESHAPE_HPP

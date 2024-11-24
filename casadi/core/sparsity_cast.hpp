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


#ifndef CASADI_SPARSITY_CAST_HPP
#define CASADI_SPARSITY_CAST_HPP

#include "mx_node.hpp"
#include <map>
#include <stack>

/// \cond INTERNAL

namespace casadi {

  class CASADI_EXPORT SparsityCast : public MXNode {
  public:

    /// Constructor
    SparsityCast(const MX& x, Sparsity sp);

    /// Destructor
    ~SparsityCast() override {}

    /// Evaluate the function (template)
    template<typename T>
    int eval_gen(const T** arg, T** res, casadi_int* iw, T* w) const;

    /// Evaluate the function numerically
    int eval(const double** arg, double** res, casadi_int* iw, double* w) const override;

    /// Evaluate the function symbolically (SX)
    int eval_sx(const SXElem** arg, SXElem** res, casadi_int* iw, SXElem* w) const override;

    /** \brief  Evaluate symbolically (MX)

        \identifier{24i} */
    void eval_mx(const std::vector<MX>& arg, std::vector<MX>& res) const override;

    /** \brief Calculate forward mode directional derivatives

        \identifier{24j} */
    void ad_forward(const std::vector<std::vector<MX> >& fseed,
                         std::vector<std::vector<MX> >& fsens) const override;

    /** \brief Calculate reverse mode directional derivatives

        \identifier{24k} */
    void ad_reverse(const std::vector<std::vector<MX> >& aseed,
                         std::vector<std::vector<MX> >& asens) const override;

    /** \brief  Propagate sparsity forward

        \identifier{24l} */
    int sp_forward(const bvec_t** arg, bvec_t** res, casadi_int* iw, bvec_t* w) const override;

    /** \brief  Propagate sparsity backwards

        \identifier{24m} */
    int sp_reverse(bvec_t** arg, bvec_t** res, casadi_int* iw, bvec_t* w) const override;

    /** \brief  Print expression

        \identifier{24n} */
    std::string disp(const std::vector<std::string>& arg) const override;

    /** \brief Generate code for the operation

        \identifier{24o} */
    void generate(CodeGenerator& g,
                  const std::vector<casadi_int>& arg,
                  const std::vector<casadi_int>& res,
                  const std::vector<bool>& arg_is_ref,
                  std::vector<bool>& res_is_ref) const override;

    /** \brief Get the operation

        \identifier{24p} */
    casadi_int op() const override { return OP_SPARSITY_CAST;}

    /// Can the operation be performed inplace (i.e. overwrite the result)
    casadi_int n_inplace() const override { return 1;}

    /// SparsityCast
    MX get_reshape(const Sparsity& sp) const override;

    /** \brief Get the nonzeros of matrix

        \identifier{24q} */
    MX get_nzref(const Sparsity& sp, const std::vector<casadi_int>& nz) const override;

    /// SparsityCast
    MX get_sparsity_cast(const Sparsity& sp) const override;

    /** \brief Check if two nodes are equivalent up to a given depth

        \identifier{24r} */
    bool is_equal(const MXNode* node, casadi_int depth) const override
    { return sameOpAndDeps(node, depth) && sparsity()==node->sparsity();}

    /// Transpose (if a dimension is one)
    MX get_transpose() const override;

    /** \brief  Check if valid function input

        \identifier{24s} */
    bool is_valid_input() const override;

    /** \brief Get the number of symbolic primitives

        \identifier{27u} */
    casadi_int n_primitives() const override;

    /** \brief Get symbolic primitives

        \identifier{27v} */
    void primitives(std::vector<MX>::iterator& it) const override;

    /// Split up an expression along primitives (template)
    template<typename T>
    void split_primitives_gen(const T& x, typename std::vector<T>::iterator& it) const;

    /// @{
    /** \brief Split up an expression along symbolic primitives

        \identifier{27w} */
    void split_primitives(const MX& x, std::vector<MX>::iterator& it) const override;
    void split_primitives(const SX& x, std::vector<SX>::iterator& it) const override;
    void split_primitives(const DM& x, std::vector<DM>::iterator& it) const override;
    /// @}

    /// Join an expression along symbolic primitives (template)
    template<typename T>
    T join_primitives_gen(typename std::vector<T>::const_iterator& it) const;

    /// @{
    /** \brief Join an expression along symbolic primitives

        \identifier{27x} */
    MX join_primitives(std::vector<MX>::const_iterator& it) const override;
    SX join_primitives(std::vector<SX>::const_iterator& it) const override;
    DM join_primitives(std::vector<DM>::const_iterator& it) const override;
    /// @}

    /** \brief Detect duplicate symbolic expressions

        \identifier{24t} */
    bool has_duplicates() const override;

    /** \brief Reset the marker for an input expression

        \identifier{24u} */
    void reset_input() const override;

    /** \brief Deserialize without type information

        \identifier{24v} */
    static MXNode* deserialize(DeserializingStream& s) { return new SparsityCast(s); }
  protected:
    /** \brief Deserializing constructor

        \identifier{24w} */
    explicit SparsityCast(DeserializingStream& s) : MXNode(s) {}
  };

} // namespace casadi
/// \endcond

#endif // CASADI_SPARSITY_CAST_HPP

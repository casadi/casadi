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


#ifndef CASADI_BINARY_MX_HPP
#define CASADI_BINARY_MX_HPP

#include "mx_node.hpp"

/// \cond INTERNAL

namespace casadi {
  /** \brief Represents any binary operation that involves two matrices
      \author Joel Andersson
      \date 2010
  */
  template<bool ScX, bool ScY>
  class CASADI_EXPORT BinaryMX : public MXNode {
  public:
    /** \brief  Constructor */
    BinaryMX(Operation op, const MX& x, const MX& y);

    /** \brief  Destructor */
    ~BinaryMX() override;

    /** \brief  Print expression */
    std::string print(const std::vector<std::string>& arg) const override;

    /** \brief Get the operation */
    int op() const override { return op_;}

    /** \brief Check if binary operation */
    bool is_binary() const override { return true;}

    /** \brief  Evaluate symbolically (MX) */
    void eval_mx(const std::vector<MX>& arg, std::vector<MX>& res) const override;

    /** \brief Calculate forward mode directional derivatives */
    void ad_forward(const std::vector<std::vector<MX> >& fseed,
                         std::vector<std::vector<MX> >& fsens) const override;

    /** \brief Calculate reverse mode directional derivatives */
    void ad_reverse(const std::vector<std::vector<MX> >& aseed,
                         std::vector<std::vector<MX> >& asens) const override;

    /** \brief  Propagate sparsity forward */
    int sp_forward(const bvec_t** arg, bvec_t** res, int* iw, bvec_t* w) const override;

    /** \brief  Propagate sparsity backwards */
    int sp_reverse(bvec_t** arg, bvec_t** res, int* iw, bvec_t* w) const override;

    /// Evaluate the function (template)
    template<typename T>
    int eval_gen(const T* const* arg, T* const* res, int* iw, T* w) const;

    /// Evaluate the function numerically
    int eval(const double** arg, double** res, int* iw, double* w) const override;

    /// Evaluate the function symbolically (SX)
    int eval_sx(const SXElem** arg, SXElem** res, int* iw, SXElem* w) const override;

    /// Can the operation be performed inplace (i.e. overwrite the result)
    int n_inplace() const override { return 2;}

    /** \brief Generate code for the operation */
    void generate(CodeGenerator& g,
                  const std::vector<int>& arg, const std::vector<int>& res) const override;

    /// Get a unary operation
    MX get_unary(int op) const override;

    /// Get a binary operation operation
    MX _get_binary(int op, const MX& y, bool scX, bool scY) const override;

    /** \brief Check if two nodes are equivalent up to a given depth */
    bool is_equal(const MXNode* node, int depth) const override {
      if (op_==node->op()) {
        if (MX::is_equal(dep(0), node->dep(0), depth-1)
            && MX::is_equal(dep(1), node->dep(1), depth-1)) {
          // If arguments are equal
          return true;
        } else {
          // If arguments are flipped
          return operation_checker<CommChecker>(op_)
            && MX::is_equal(dep(1), node->dep(0), depth-1)
            && MX::is_equal(dep(0), node->dep(1), depth-1);
        }
      } else {
        return false;
      }
    }

    //! \brief Operation
    Operation op_;

  };

} // namespace casadi
/// \endcond

#endif // CASADI_BINARY_MX_HPP

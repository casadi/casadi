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


#ifndef CASADI_BINARY_MX_HPP
#define CASADI_BINARY_MX_HPP

#include "mx_node.hpp"

/// \cond INTERNAL

namespace casadi {
  /** \brief Represents any binary operation that involves two matrices

      \author Joel Andersson
      \date 2010

      \identifier{1fn} */
  template<bool ScX, bool ScY>
  class CASADI_EXPORT BinaryMX : public MXNode {
  public:
    /** \brief  Constructor

        \identifier{1fo} */
    BinaryMX(Operation op, const MX& x, const MX& y);

    /** \brief  Destructor

        \identifier{1fp} */
    ~BinaryMX() override;

    /** \brief  Print expression

        \identifier{1fq} */
    std::string disp(const std::vector<std::string>& arg) const override;

    /** \brief Get the operation

        \identifier{1fr} */
    casadi_int op() const override { return op_;}

    /** \brief Check if binary operation

        \identifier{1fs} */
    bool is_binary() const override { return true;}

    /** \brief  Evaluate symbolically (MX)

        \identifier{1ft} */
    void eval_mx(const std::vector<MX>& arg, std::vector<MX>& res) const override;

    /** \brief Calculate forward mode directional derivatives

        \identifier{1fu} */
    void ad_forward(const std::vector<std::vector<MX> >& fseed,
                         std::vector<std::vector<MX> >& fsens) const override;

    /** \brief Calculate reverse mode directional derivatives

        \identifier{1fv} */
    void ad_reverse(const std::vector<std::vector<MX> >& aseed,
                         std::vector<std::vector<MX> >& asens) const override;

    /** \brief  Propagate sparsity forward

        \identifier{1fw} */
    int sp_forward(const bvec_t** arg, bvec_t** res, casadi_int* iw, bvec_t* w) const override;

    /** \brief  Propagate sparsity backwards

        \identifier{1fx} */
    int sp_reverse(bvec_t** arg, bvec_t** res, casadi_int* iw, bvec_t* w) const override;

    /// Evaluate the function (template)
    template<typename T>
    int eval_gen(const T* const* arg, T* const* res, casadi_int* iw, T* w) const;

    /// Evaluate the function numerically
    int eval(const double** arg, double** res, casadi_int* iw, double* w) const override;

    /// Evaluate the function symbolically (SX)
    int eval_sx(const SXElem** arg, SXElem** res, casadi_int* iw, SXElem* w) const override;

    /// Can the operation be performed inplace (i.e. overwrite the result)
    casadi_int n_inplace() const override { return 2;}

    /** \brief Generate code for the operation

        \identifier{1fy} */
    void generate(CodeGenerator& g,
                  const std::vector<casadi_int>& arg,
                  const std::vector<casadi_int>& res) const override;

    /** \brief Serialize an object without type information

        \identifier{1fz} */
    void serialize_body(SerializingStream& s) const override;

    /** \brief Serialize type information

        \identifier{1g0} */
    void serialize_type(SerializingStream& s) const override;

    /** \brief Deserialize with type disambiguation

        \identifier{1g1} */
    static MXNode* deserialize(DeserializingStream& s);

    /// Get a unary operation
    MX get_unary(casadi_int op) const override;

    /// Get a binary operation operation
    MX _get_binary(casadi_int op, const MX& y, bool scX, bool scY) const override;

    /** \brief Check if two nodes are equivalent up to a given depth

        \identifier{1g2} */
    bool is_equal(const MXNode* node, casadi_int depth) const override {
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

    /// Solve a system of linear equations, upper triangular A
    MX get_solve_triu(const MX& r, bool tr) const override;

    /// Solve a system of linear equations, lower triangular A
    MX get_solve_tril(const MX& r, bool tr) const override;

    //! \brief Operation
    Operation op_;

    /** \brief Deserializing constructor

        \identifier{1g3} */
    explicit BinaryMX(DeserializingStream& s);

  };

} // namespace casadi
/// \endcond

#endif // CASADI_BINARY_MX_HPP

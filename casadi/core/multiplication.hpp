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
  */
  class CASADI_EXPORT Multiplication : public MXNode {
  public:

    /** \brief  Constructor */
    Multiplication(const MX& z, const MX& x, const MX& y);

    /** \brief  Destructor */
    ~Multiplication() override {}

    /** \brief  Print expression */
    std::string print(const std::vector<std::string>& arg) const override;

    /** \brief Generate code for the operation */
    void generate(CodeGenerator& g, const std::string& mem,
                          const std::vector<int>& arg, const std::vector<int>& res) const override;

    /// Evaluate the function (template)
    template<typename T>
    void evalGen(const T** arg, T** res, int* iw, T* w, int mem) const;

    /// Evaluate the function numerically
    void eval(const double** arg, double** res, int* iw, double* w, int mem) const override;

    /// Evaluate the function symbolically (SX)
    void eval_sx(const SXElem** arg, SXElem** res, int* iw, SXElem* w, int mem) const override;

    /** \brief  Evaluate symbolically (MX) */
    void eval_mx(const std::vector<MX>& arg, std::vector<MX>& res) const override;

    /** \brief Calculate forward mode directional derivatives */
    void eval_forward(const std::vector<std::vector<MX> >& fseed,
                         std::vector<std::vector<MX> >& fsens) const override;

    /** \brief Calculate reverse mode directional derivatives */
    void eval_reverse(const std::vector<std::vector<MX> >& aseed,
                         std::vector<std::vector<MX> >& asens) const override;

    /** \brief  Propagate sparsity forward */
    void sp_fwd(const bvec_t** arg, bvec_t** res, int* iw, bvec_t* w, int mem) const override;

    /** \brief  Propagate sparsity backwards */
    void sp_rev(bvec_t** arg, bvec_t** res, int* iw, bvec_t* w, int mem) const override;

    /** \brief Get the operation */
    int op() const override { return OP_MTIMES;}

    /// Can the operation be performed inplace (i.e. overwrite the result)
    int numInplace() const override { return 1;}

    /** \brief Check if two nodes are equivalent up to a given depth */
    bool is_equal(const MXNode* node, int depth) const override {
      return sameOpAndDeps(node, depth) && dynamic_cast<const Multiplication*>(node)!=0;
    }

    /** \brief Get required length of w field */
    size_t sz_w() const override { return sparsity().size1();}
  };


  /** \brief An MX atomic for matrix-matrix product,
             note that the factor must be provided transposed
      \author Joel Andersson
      \date 2010
  */
  class CASADI_EXPORT DenseMultiplication : public Multiplication{
  public:

    /** \brief  Constructor */
    DenseMultiplication(const MX& z, const MX& x, const MX& y)
        : Multiplication(z, x, y) {}

    /** \brief  Destructor */
    ~DenseMultiplication() override {}

    /** \brief Generate code for the operation */
    void generate(CodeGenerator& g, const std::string& mem,
                          const std::vector<int>& arg, const std::vector<int>& res) const override;
  };


} // namespace casadi
/// \endcond

#endif // CASADI_MULTIPLICATION_HPP

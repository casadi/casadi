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


#ifndef CASADI_EINSTEIN_HPP
#define CASADI_EINSTEIN_HPP

#include "mx_node.hpp"

/// \cond INTERNAL

namespace casadi {
  /** \brief An MX atomic for an Einstein product,
      \author Joris Gillis
      \date 2016
  */
  class CASADI_EXPORT Einstein : public MXNode {
  public:

    /** \brief  Constructor */
    Einstein(const MX& C, const MX& A, const MX& B,
      const std::vector<int>& dim_c, const std::vector<int>& dim_a, const std::vector<int>& dim_b,
      const std::vector<int>& c, const std::vector<int>& a, const std::vector<int>& b);

    /** \brief  Destructor */
    virtual ~Einstein() {}

    /** \brief  Print expression */
    virtual std::string print(const std::vector<std::string>& arg) const;

    /** \brief Generate code for the operation */
    virtual void generate(CodeGenerator& g, const std::string& mem,
                          const std::vector<int>& arg, const std::vector<int>& res) const;

    /// Evaluate the function (template)
    template<typename T>
    void evalGen(const T** arg, T** res, int* iw, T* w, int mem) const;

    /// Evaluate the function numerically
    virtual void eval(const double** arg, double** res, int* iw, double* w, int mem) const;

    /// Evaluate the function symbolically (SX)
    virtual void eval_sx(const SXElem** arg, SXElem** res, int* iw, SXElem* w, int mem) const;

    /** \brief  Evaluate symbolically (MX) */
    virtual void eval_mx(const std::vector<MX>& arg, std::vector<MX>& res) const;

    /** \brief Calculate forward mode directional derivatives */
    virtual void eval_forward(const std::vector<std::vector<MX> >& fseed,
                         std::vector<std::vector<MX> >& fsens) const;

    /** \brief Calculate reverse mode directional derivatives */
    virtual void eval_reverse(const std::vector<std::vector<MX> >& aseed,
                         std::vector<std::vector<MX> >& asens) const;

    /** \brief  Propagate sparsity forward */
    virtual void sp_fwd(const bvec_t** arg, bvec_t** res, int* iw, bvec_t* w, int mem) const;

    /** \brief  Propagate sparsity backwards */
    virtual void sp_rev(bvec_t** arg, bvec_t** res, int* iw, bvec_t* w, int mem) const;

    /** \brief Get the operation */
    virtual int op() const { return OP_EINSTEIN;}

    /// Can the operation be performed inplace (i.e. overwrite the result)
    virtual int numInplace() const { return 1;}

    /** \brief Check if two nodes are equivalent up to a given depth */
    virtual bool is_equal(const MXNode* node, int depth) const {
      return sameOpAndDeps(node, depth) && dynamic_cast<const Einstein*>(node)!=0;
    }

    /** \brief Get required length of w field */
    virtual size_t sz_w() const { return sparsity().size1();}

    /// Dimensions of tensors A B C
    std::vector<int> dim_c_, dim_a_, dim_b_;
    /// Einstein indices
    std::vector<int> c_, a_, b_;

    std::vector<int> iter_dims_;

    std::vector<int> strides_a_;
    std::vector<int> strides_b_;
    std::vector<int> strides_c_;

    int n_iter_;

  };


} // namespace casadi
/// \endcond

#endif // CASADI_EINSTEIN_HPP

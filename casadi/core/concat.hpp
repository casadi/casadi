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


#ifndef CASADI_CONCAT_HPP
#define CASADI_CONCAT_HPP

#include "mx_node.hpp"
#include <map>
#include <stack>

/// \cond INTERNAL

namespace casadi {
  /** \brief Concatenation: Join multiple expressions stacking the nonzeros
      \author Joel Andersson
      \date 2014
  */
  class CASADI_EXPORT Concat : public MXNode {
  public:

    /// Constructor
    Concat(const std::vector<MX>& x);

    /// Destructor
    ~Concat() override = 0;

    /// Evaluate the function (template)
    template<typename T>
    void evalGen(const T* const* arg, T* const* res, int* iw, T* w) const;

    /// Evaluate the function numerically
    void eval(const double** arg, double** res, int* iw, double* w, int mem) const override;

    /// Evaluate the function symbolically (SX)
    void eval_sx(const SXElem** arg, SXElem** res, int* iw, SXElem* w, int mem) const override;

    /** \brief  Propagate sparsity forward */
    void sp_fwd(const bvec_t** arg, bvec_t** res, int* iw, bvec_t* w, int mem) const override;

    /** \brief  Propagate sparsity backwards */
    void sp_rev(bvec_t** arg, bvec_t** res, int* iw, bvec_t* w, int mem) const override;

    /** \brief Generate code for the operation */
    void generate(CodeGenerator& g, const std::string& mem,
                          const std::vector<int>& arg, const std::vector<int>& res) const override;

    /// Get the nonzeros of matrix
    MX getGetNonzeros(const Sparsity& sp, const std::vector<int>& nz) const override;

    /** \brief Check if two nodes are equivalent up to a given depth */
    bool is_equal(const MXNode* node, int depth) const override {
      return sameOpAndDeps(node, depth);
    }

    /** \brief  Check if valid function input */
    bool is_valid_input() const override;

    /** \brief Get the number of symbolic primitives */
    int n_primitives() const override;

    /** \brief Get symbolic primitives */
    void primitives(std::vector<MX>::iterator& it) const override;

    /** \brief Detect duplicate symbolic expressions */
    bool has_duplicates() const override;

    /** \brief Reset the marker for an input expression */
    void reset_input() const override;
  };


  /** \brief Horizontal concatenation
      \author Joel Andersson
      \date 2013
  */
  class CASADI_EXPORT Horzcat : public Concat {
  public:

    /// Constructor
    Horzcat(const std::vector<MX>& x);

    /// Destructor
    ~Horzcat() override {}

    /** \brief  Print expression */
    std::string print(const std::vector<std::string>& arg) const override;

    /** \brief  Evaluate symbolically (MX) */
    void eval_mx(const std::vector<MX>& arg, std::vector<MX>& res) const override;

    /** \brief Calculate forward mode directional derivatives */
    void eval_forward(const std::vector<std::vector<MX> >& fseed,
                         std::vector<std::vector<MX> >& fsens) const override;

    /** \brief Calculate reverse mode directional derivatives */
    void eval_reverse(const std::vector<std::vector<MX> >& aseed,
                         std::vector<std::vector<MX> >& asens) const override;

    /** \brief Get the operation */
    int op() const override { return OP_HORZCAT;}

    /** \brief Split up an expression along symbolic primitives */
    void split_primitives(const MX& x, std::vector<MX>::iterator& it) const override;

    /** \brief Join an expression along symbolic primitives */
    MX join_primitives(std::vector<MX>::const_iterator& it) const override;

    /** \brief Get offsets for split */
    std::vector<int> off() const;
  };

  /** \brief Vertical concatenation of vectors
      \author Joel Andersson
      \date 2014
  */
  class CASADI_EXPORT Vertcat : public Concat {
  public:

    /// Constructor
    Vertcat(const std::vector<MX>& x);

    /// Destructor
    ~Vertcat() override {}

    /** \brief  Print expression */
    std::string print(const std::vector<std::string>& arg) const override;

    /** \brief  Evaluate symbolically (MX) */
    void eval_mx(const std::vector<MX>& arg, std::vector<MX>& res) const override;

    /** \brief Calculate forward mode directional derivatives */
    void eval_forward(const std::vector<std::vector<MX> >& fseed,
                         std::vector<std::vector<MX> >& fsens) const override;

    /** \brief Calculate reverse mode directional derivatives */
    void eval_reverse(const std::vector<std::vector<MX> >& aseed,
                         std::vector<std::vector<MX> >& asens) const override;

    /** \brief Get the operation */
    int op() const override { return OP_VERTCAT;}

    /** \brief Split up an expression along symbolic primitives */
    void split_primitives(const MX& x, std::vector<MX>::iterator& it) const override;

    /** \brief Join an expression along symbolic primitives */
    MX join_primitives(std::vector<MX>::const_iterator& it) const override;

    /** \brief Get offsets for split */
    std::vector<int> off() const;
  };

  /** \brief Diagonal concatenation of matrices
      \author Joris Gillis
      \date 2014
  */
  class CASADI_EXPORT Diagcat : public Concat {
  public:

    /// Constructor
    Diagcat(const std::vector<MX>& x);

    /// Destructor
    ~Diagcat() override {}

    /** \brief  Print expression */
    std::string print(const std::vector<std::string>& arg) const override;

    /** \brief  Evaluate symbolically (MX) */
    void eval_mx(const std::vector<MX>& arg, std::vector<MX>& res) const override;

    /** \brief Calculate forward mode directional derivatives */
    void eval_forward(const std::vector<std::vector<MX> >& fseed,
                         std::vector<std::vector<MX> >& fsens) const override;

    /** \brief Calculate reverse mode directional derivatives */
    void eval_reverse(const std::vector<std::vector<MX> >& aseed,
                         std::vector<std::vector<MX> >& asens) const override;

    /** \brief Get the operation */
    int op() const override { return OP_DIAGCAT;}

    /** \brief Split up an expression along symbolic primitives */
    void split_primitives(const MX& x, std::vector<MX>::iterator& it) const override;

    /** \brief Join an expression along symbolic primitives */
    MX join_primitives(std::vector<MX>::const_iterator& it) const override;

    /** \brief Get offsets for split */
    std::pair<std::vector<int>, std::vector<int> > off() const;
  };

} // namespace casadi
/// \endcond

#endif // CASADI_CONCAT_HPP

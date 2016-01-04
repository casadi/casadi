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
    virtual ~Concat() = 0;

    /// Evaluate the function (template)
    template<typename T>
    void evalGen(const T* const* arg, T* const* res, int* iw, T* w) const;

    /// Evaluate the function numerically
    virtual void eval(const double** arg, double** res, int* iw, double* w, int mem) const;

    /// Evaluate the function symbolically (SX)
    virtual void eval_sx(const SXElem** arg, SXElem** res, int* iw, SXElem* w, int mem);

    /** \brief  Propagate sparsity forward */
    virtual void spFwd(const bvec_t** arg, bvec_t** res, int* iw, bvec_t* w, int mem);

    /** \brief  Propagate sparsity backwards */
    virtual void spAdj(bvec_t** arg, bvec_t** res, int* iw, bvec_t* w, int mem);

    /** \brief Generate code for the operation */
    virtual void generate(CodeGenerator& g, const std::string& mem,
                          const std::vector<int>& arg, const std::vector<int>& res) const;

    /// Get the nonzeros of matrix
    virtual MX getGetNonzeros(const Sparsity& sp, const std::vector<int>& nz) const;

    /** \brief Check if two nodes are equivalent up to a given depth */
    virtual bool is_equal(const MXNode* node, int depth) const {
      return sameOpAndDeps(node, depth);
    }

    /** \brief  Check if valid function input */
    virtual bool is_valid_input() const;

    /** \brief Get the number of symbolic primitives */
    virtual int n_primitives() const;

    /** \brief Get symbolic primitives */
    virtual void primitives(std::vector<MX>::iterator& it) const;

    /** \brief Detect duplicate symbolic expressions */
    virtual bool has_duplicates();

    /** \brief Reset the marker for an input expression */
    virtual void resetInput();
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
    virtual ~Horzcat() {}

    /** \brief  Print expression */
    virtual std::string print(const std::vector<std::string>& arg) const;

    /** \brief  Evaluate symbolically (MX) */
    virtual void eval_mx(const std::vector<MX>& arg, std::vector<MX>& res);

    /** \brief Calculate forward mode directional derivatives */
    virtual void evalFwd(const std::vector<std::vector<MX> >& fseed,
                         std::vector<std::vector<MX> >& fsens);

    /** \brief Calculate reverse mode directional derivatives */
    virtual void evalAdj(const std::vector<std::vector<MX> >& aseed,
                         std::vector<std::vector<MX> >& asens);

    /** \brief Get the operation */
    virtual int op() const { return OP_HORZCAT;}

    /** \brief Split up an expression along symbolic primitives */
    virtual void split_primitives(const MX& x, std::vector<MX>::iterator& it) const;

    /** \brief Join an expression along symbolic primitives */
    virtual MX join_primitives(std::vector<MX>::const_iterator& it) const;

    /** \brief Get offsets for split */
    std::vector<int> offset() const;
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
    virtual ~Vertcat() {}

    /** \brief  Print expression */
    virtual std::string print(const std::vector<std::string>& arg) const;

    /** \brief  Evaluate symbolically (MX) */
    virtual void eval_mx(const std::vector<MX>& arg, std::vector<MX>& res);

    /** \brief Calculate forward mode directional derivatives */
    virtual void evalFwd(const std::vector<std::vector<MX> >& fseed,
                         std::vector<std::vector<MX> >& fsens);

    /** \brief Calculate reverse mode directional derivatives */
    virtual void evalAdj(const std::vector<std::vector<MX> >& aseed,
                         std::vector<std::vector<MX> >& asens);

    /** \brief Get the operation */
    virtual int op() const { return OP_VERTCAT;}

    /** \brief Split up an expression along symbolic primitives */
    virtual void split_primitives(const MX& x, std::vector<MX>::iterator& it) const;

    /** \brief Join an expression along symbolic primitives */
    virtual MX join_primitives(std::vector<MX>::const_iterator& it) const;

    /** \brief Get offsets for split */
    std::vector<int> offset() const;
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
    virtual ~Diagcat() {}

    /** \brief  Print expression */
    virtual std::string print(const std::vector<std::string>& arg) const;

    /** \brief  Evaluate symbolically (MX) */
    virtual void eval_mx(const std::vector<MX>& arg, std::vector<MX>& res);

    /** \brief Calculate forward mode directional derivatives */
    virtual void evalFwd(const std::vector<std::vector<MX> >& fseed,
                         std::vector<std::vector<MX> >& fsens);

    /** \brief Calculate reverse mode directional derivatives */
    virtual void evalAdj(const std::vector<std::vector<MX> >& aseed,
                         std::vector<std::vector<MX> >& asens);

    /** \brief Get the operation */
    virtual int op() const { return OP_DIAGCAT;}

    /** \brief Split up an expression along symbolic primitives */
    virtual void split_primitives(const MX& x, std::vector<MX>::iterator& it) const;

    /** \brief Join an expression along symbolic primitives */
    virtual MX join_primitives(std::vector<MX>::const_iterator& it) const;

    /** \brief Get offsets for split */
    std::pair<std::vector<int>, std::vector<int> > offset() const;
  };

} // namespace casadi
/// \endcond

#endif // CASADI_CONCAT_HPP

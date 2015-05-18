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
  */
  class CASADI_EXPORT Reshape : public MXNode {
  public:

    /// Constructor
    Reshape(const MX& x, Sparsity sp);

    /// Clone function
    virtual Reshape* clone() const;

    /// Destructor
    virtual ~Reshape() {}

    /// Evaluate the function (template)
    template<typename T>
    void evalGen(const T** arg, T** res, int* iw, T* w);

    /// Evaluate the function numerically
    virtual void evalD(const double** arg, double** res, int* iw, double* w);

    /// Evaluate the function symbolically (SX)
    virtual void evalSX(const SXElement** arg, SXElement** res, int* iw, SXElement* w);

    /** \brief  Evaluate symbolically (MX) */
    virtual void evalMX(const std::vector<MX>& arg, std::vector<MX>& res);

    /** \brief Calculate forward mode directional derivatives */
    virtual void evalFwd(const std::vector<std::vector<MX> >& fseed,
                         std::vector<std::vector<MX> >& fsens);

    /** \brief Calculate reverse mode directional derivatives */
    virtual void evalAdj(const std::vector<std::vector<MX> >& aseed,
                         std::vector<std::vector<MX> >& asens);

    /** \brief  Propagate sparsity forward */
    virtual void spFwd(const bvec_t** arg, bvec_t** res, int* iw, bvec_t* w);

    /** \brief  Propagate sparsity backwards */
    virtual void spAdj(bvec_t** arg, bvec_t** res, int* iw, bvec_t* w);

    /** \brief  Print expression */
    virtual std::string print(const std::vector<std::string>& arg) const;

    /** \brief Generate code for the operation */
    virtual void generate(const std::vector<int>& arg, const std::vector<int>& res,
                          CodeGenerator& g) const;

    /** \brief Get the operation */
    virtual int getOp() const { return OP_RESHAPE;}

    /// Can the operation be performed inplace (i.e. overwrite the result)
    virtual int numInplace() const { return 1;}

    /// Reshape
    virtual MX getReshape(const Sparsity& sp) const;

    /** \brief Check if two nodes are equivalent up to a given depth */
    virtual bool zz_isEqual(const MXNode* node, int depth) const
    { return sameOpAndDeps(node, depth) && sparsity()==node->sparsity();}

    /// Transpose (if a dimension is one)
    virtual MX getTranspose() const;

    /** \brief  Check if valid function input */
    virtual bool isValidInput() const;

    /** \brief Get the number of symbolic primitives */
    virtual int numPrimitives() const;

    /** \brief Get symbolic primitives */
    virtual void getPrimitives(std::vector<MX>::iterator& it) const;

    /** \brief Split up an expression along symbolic primitives */
    virtual void splitPrimitives(const MX& x, std::vector<MX>::iterator& it) const;

    /** \brief Join an expression along symbolic primitives */
    virtual MX joinPrimitives(std::vector<MX>::const_iterator& it) const;

    /** \brief Detect duplicate symbolic expressions */
    virtual bool hasDuplicates();

    /** \brief Reset the marker for an input expression */
    virtual void resetInput();
  };

} // namespace casadi
/// \endcond

#endif // CASADI_RESHAPE_HPP

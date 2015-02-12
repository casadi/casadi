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
    void evaluateGen(const T* const* input, T** output, int* itmp, T* rtmp);

    /// Evaluate the function numerically
    virtual void evaluateD(const double* const* input, double** output, int* itmp, double* rtmp);

    /// Evaluate the function symbolically (SX)
    virtual void evaluateSX(const SXElement* const* input, SXElement** output,
                            int* itmp, SXElement* rtmp);

    /// Evaluate the function symbolically (MX)
    virtual void eval(const MXPtrV& input, MXPtrV& output);

    /** \brief Calculate forward mode directional derivatives */
    virtual void evalFwd(const MXPtrVV& fwdSeed, MXPtrVV& fwdSens);

    /** \brief Calculate reverse mode directional derivatives */
    virtual void evalAdj(MXPtrVV& adjSeed, MXPtrVV& adjSens);

    /// Propagate sparsity
    virtual void propagateSparsity(double** input, double** output, bool fwd);

    /// Print a part of the expression */
    virtual void printPart(std::ostream &stream, int part) const;

    /** \brief Generate code for the operation */
    virtual void generateOperation(std::ostream &stream, const std::vector<int>& arg,
                                   const std::vector<int>& res, CodeGenerator& gen) const;

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
  };

} // namespace casadi
/// \endcond

#endif // CASADI_RESHAPE_HPP

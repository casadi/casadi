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


#ifndef CASADI_INNER_PROD_HPP
#define CASADI_INNER_PROD_HPP

#include "mx_node.hpp"
#include <map>
#include <stack>

/// \cond INTERNAL

namespace casadi {
  /** \brief Matrix inner_prod
      \author Joel Andersson
      \date 2013
  */
  class CASADI_CORE_EXPORT InnerProd : public MXNode {
  public:

    /// Constructor
    InnerProd(const MX& x, const MX& y);

    /// Clone function
    virtual InnerProd* clone() const { return new InnerProd(*this);}

    /// Destructor
    virtual ~InnerProd() {}

    /** \brief  Evaluate the function numerically */
    virtual void evaluateD(const DMatrixPtrV& input, DMatrixPtrV& output, std::vector<int>& itmp,
                           std::vector<double>& rtmp);

    /** \brief  Evaluate the function symbolically (SX) */
    virtual void evaluateSX(const SXPtrV& input, SXPtrV& output, std::vector<int>& itmp,
                            std::vector<SXElement>& rtmp);

    /// Evaluate the function (template)
    template<typename T, typename MatV, typename MatVV>
    void evaluateGen(const MatV& input, MatV& output, std::vector<int>& itmp, std::vector<T>& rtmp);

    /** \brief  Propagate sparsity */
    virtual void propagateSparsity(DMatrixPtrV& input, DMatrixPtrV& output, bool fwd);

    /// Evaluate the function symbolically (MX)
    virtual void evaluateMX(const MXPtrV& input, MXPtrV& output, const MXPtrVV& fwdSeed,
                            MXPtrVV& fwdSens, const MXPtrVV& adjSeed, MXPtrVV& adjSens,
                            bool output_given);

    /** \brief Generate code for the operation */
    void generateOperation(std::ostream &stream, const std::vector<std::string>& arg,
                           const std::vector<std::string>& res, CodeGenerator& gen) const;

    /// Print a part of the expression */
    virtual void printPart(std::ostream &stream, int part) const;

    /** \brief Get the operation */
    virtual int getOp() const { return OP_INNER_PROD;}
  };


} // namespace casadi
/// \endcond

#endif // CASADI_INNER_PROD_HPP

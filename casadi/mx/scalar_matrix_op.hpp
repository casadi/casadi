/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A symbolic framework for dynamic optimization.
 *    Copyright (C) 2010 by Joel Andersson, Moritz Diehl, K.U.Leuven. All rights reserved.
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

#ifndef SCALAR_MATRIX_OP_HPP
#define SCALAR_MATRIX_OP_HPP

#include "mx_node.hpp"

namespace CasADi{
/** Represents any binary operation that involves a matrix and a scalar
  \author Joel Andersson 
  \date 2010
*/
class ScalarMatrixOp : public MXNode{
  public:

    /** \brief  Constructor */
    ScalarMatrixOp(Operation op, const MX& x, MX y);

    /** \brief  Destructor */
    virtual ~ScalarMatrixOp(){}

    /** \brief  Clone function */
    virtual ScalarMatrixOp * clone() const;

    /** \brief  Print */
    virtual void print(std::ostream &stream, const std::vector<std::string>& args) const;

    /** \brief  Evaluate the function and store the result in the node */
    virtual void evaluate(const VDptr& input, Dptr& output, const VVDptr& fwdSeed, VDptr& fwdSens, const VDptr& adjSeed, VVDptr& adjSens, int nfwd, int nadj);

    /// Symbolic forward sensitivities
    virtual MX adFwd(const std::vector<MX>& jx);

    /** \brief  Evaluate symbolically (SX) */
    virtual void evaluateSX(const std::vector<SXMatrix*> &input, SXMatrix& output);

    //! \brief Operation
    Operation op;
};


} // namespace CasADi


#endif // SCALAR_MATRIX_OP_HPP

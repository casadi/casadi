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

#ifndef MATRIX_MATRIX_OP_HPP
#define MATRIX_MATRIX_OP_HPP

#include "mx_node.hpp"

namespace CasADi{
/** \brief Represents any binary operation that involves two matrices 
  \author Joel Andersson 
  \date 2010	
*/
class MatrixMatrixOp : public MXNode{
  public:

    /** \brief  Constructor */
    MatrixMatrixOp (OPERATION op, const MX& x, const MX& y);

    /** \brief  Clone function */
    virtual MatrixMatrixOp * clone() const;

    /** \brief  Print */
    virtual void print(std::ostream &stream=std::cout) const;

    /** \brief  Evaluate the function and store the result in the node */
    virtual void evaluate(const VDptr& input, Dptr& output, const VVDptr& fwdSeed, VDptr& fwdSens, const VDptr& adjSeed, VVDptr& adjSens, int nfwd, int nadj);

  protected:
    //! \brief Operation
    OPERATION op;
    
    /// Does the two matrices have the same sparsity
    bool same_sparsity_;
    
};

} // namespace CasADi


#endif // MATRIX_MATRIX_OP_HPP

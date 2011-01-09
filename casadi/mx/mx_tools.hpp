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

#ifndef MX_TOOLS_HPP
#define MX_TOOLS_HPP

#include "mx.hpp"

namespace CasADi{

//@{
/** \brief  concatenate */
MX vertcat(const std::vector<MX>& comp);
MX horzcat(const std::vector<MX>& comp);
MX vertcat(const MX& a, const MX& b);
MX horzcat(const MX& a, const MX& b);
//@}

/** \brief  Take the 2-norm of a MX
Internally represented by Norm2
*/
MX norm_2(const MX &x);
/** \brief  Take the 1-norm of a MX
Internally represented by Norm1
*/
MX norm_1(const MX &x);
/** \brief  Take the infinity-norm of a MX
Internally represented by NormInf
*/
MX norm_inf(const MX &x);
/** \brief  Take the transpose of a MX 
Internally represented by Transpose
*/
MX trans(const MX &x); // transpose
/** \brief  Take the matrix product of 2 MX objects */
MX prod(const MX &x, const MX &y); // matrix product
/** \brief  Take the inner product of two vectors 
        Equals
        \code
        trans(x)*y
        \endcode
        with x and y vectors
*/
MX inner_prod(const MX &x, const MX &y); // 
/** \brief  Take the outer product of two vectors 
        Equals
        \code
        x*trans(y)
        \endcode
         with x and y vectors
*/
MX outer_prod(const MX &x, const MX &y); // x*trans(y) with x and y vectors

/** \brief Branching on MX nodes
Ternary operator, "cond ? if_true : if_false"
Internally represented by IfElseNode.
*/
MX if_else(const MX &cond, const MX &if_true, const MX &if_false); 

#ifndef SWIG
//! \brief Returns a reshaped version of the MX
MX reshape(const MX &x, int n, int m);
#endif // SWIG

//! \brief Returns a reshaped version of the MX, dimensions as a vector
MX reshape(const MX &x, const std::vector<int> sz);

/** \brief Returns a flattened version of the MX
    Flattening is a cheap (non-copying) operation
    Same as reshape(x, x.numel(),1)
*/
MX flatten(const MX &x);
  
} // namespace CasADi

#endif // MX_TOOLS_HPP

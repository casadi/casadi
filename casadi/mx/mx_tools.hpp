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

/** \brief  concatenate vertically */
MX vertcat(const std::vector<MX>& comp);

/** \brief  concatenate horizontally */
MX horzcat(const std::vector<MX>& comp);

/** \brief  concatenate vertically while vectorizing all arguments with getNZ */
MX veccat(const std::vector<MX>& comp);

#ifndef SWIG
/** \brief  concatenate vertically, two matrices */
MX vertcat(const MX& a, const MX& b);

/** \brief  concatenate horizontally, two matrices */
MX horzcat(const MX& a, const MX& b);
#endif // SWIG

#ifndef SWIG
/**
Apply a function f to each element in a vector
*/
std::vector<MX> applymap(MX (*f)(const MX& ),const std::vector<MX>&);

/**
Apply a function f to each element in a vector
*/
void applymap(void (*f)(MX&), std::vector<MX>&);
#endif // SWIG

/** \brief  Take the 2-norm of a MX
Internally represented by Norm2
*/
MX norm_2(const MX &x);

/** \brief  Take the 2-norm squared of a MX
Internally represented by Norm2
*/
MX norm_22(const MX &x);

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
*/
MX if_else(const MX &cond, const MX &if_true, const MX &if_false); 

/** \brief Conditional evaluation
cond ? if_true : 0
*/
MX if_else_zero(const MX &cond, const MX &if_true); 

#ifndef SWIG
//! \brief Returns a reshaped version of the MX
MX reshape(const MX &x, int n, int m);
#endif // SWIG

//! \brief Returns a reshaped version of the MX, dimensions as a vector
MX reshape(const MX &x, const std::vector<int> sz);

//! \brief Reshape the MX
MX reshape(const MX &x, const CRSSparsity& sp);

/** \brief Returns a flattened version of the MX
    Flattening is a cheap (non-copying) operation
    Same as reshape(x, x.numel(),1)
*/
MX vec(const MX &x);

/** \brief Returns a flattened version of the MX, preserving only nonzeros
*/
MX vecNZ(const MX &x);

/** \brief  Unite two matrices no overlapping sparsity */
MX unite(const MX& A, const MX& B);

/** \brief  check if symbolic */
bool isSymbolic(const MX& ex);

/** \brief  check if identity */
bool isIdentity(const MX& ex);

/** \brief  check if zero */
bool isZero(const MX& ex);

MX trace(const MX& A);

/** \brief Repeat matrix A n times vertically and m times horizontally */
MX repmat(const MX &A, int n, int m); 

/** \brief create a clipped view into a matrix
Create a sparse matrix from a dense matrix A, with sparsity pattern sp
**/
//MX clip(const MX& A, const CRSSparsity& sp);

/** \brief Lift a node */
MX lift(const MX& x);

/** \brief  Make the matrix dense */
void makeDense(MX& x);

/** \brief  Create a parent MX on which all given MX's will depend.

 In some sense, this function is the inverse of 
 
 \param deps  Must all be symbolic matrices.
 */
MX createParent(std::vector<MX> &deps);

/** \brief  Create a parent MX on which a bunch of MX's (sizes given as argument) will depend
 */
std::pair<MX, std::vector<MX> > createParent(const std::vector<MX> &deps);


/** \brief  Create a parent MX on which a bunch of MX's (sizes given as argument) will depend
 */
std::pair<MX, std::vector<MX> > createParent(const std::vector<CRSSparsity> &deps);

// Equality
MX operator==(const MX& a, const MX& b);
MX operator>=(const MX& a, const MX& b);
MX operator<=(const MX& a, const MX& b);
#ifndef SWIG
MX operator!(const MX& a);
#endif // SWIG

/** \brief  Get the diagonal */
MX diag(const MX& x);

} // namespace CasADi

#ifdef SWIG
// Template instantiations
%template(Pair_MX_MXVector) std::pair<CasADi::MX, std::vector<CasADi::MX> >;
#endif // SWIG


#endif // MX_TOOLS_HPP

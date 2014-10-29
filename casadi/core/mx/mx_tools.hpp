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


#ifndef CASADI_MX_TOOLS_HPP
#define CASADI_MX_TOOLS_HPP

#include "mx.hpp"

#include "../matrix/generic_matrix_tools.hpp"
#include "../matrix/generic_expression_tools.hpp"
#include "../function/linear_solver.hpp"

namespace casadi {

/**
\ingroup expression_tools
@{
*/

  /** \brief  concatenate horizontally
  *
  *  horzcat(horzsplit(x, ...)) = x
  */
  CASADI_CORE_EXPORT MX horzcat(const std::vector<MX>& x);

  /** \brief  split horizontally, retaining groups of cols
  * \param output_offset List of all start cols for each group
  *      the last col group will run to the end.
  *
  *   horzcat(horzsplit(x, ...)) = x
  */
  CASADI_CORE_EXPORT std::vector<MX> horzsplit(const MX& x,
                                                   const std::vector<int>& output_offset);

  /** \brief  split horizontally, retaining fixed-sized groups of cols
  * \param incr Size of each group of cols
  *
  *  horzcat(horzsplit(x, ...)) = x
  */
  CASADI_CORE_EXPORT std::vector<MX> horzsplit(const MX& x, int incr=1);

  /** \brief  concatenate vertically
  *
  *   vertcat(vertsplit(x, ...)) = x
  */
  CASADI_CORE_EXPORT MX vertcat(const std::vector<MX>& comp);

  /** \brief  split vertically, retaining groups of rows
  * \param output_offset List of all start rows for each group
  *      the last row group will run to the end.
  *
  *   vertcat(vertsplit(x, ...)) = x
  */
  CASADI_CORE_EXPORT std::vector<MX> vertsplit(const MX& x,
                                                   const std::vector<int>& output_offset);

  /** \brief  split vertically, retaining fixed-sized groups of rows
  * \param incr Size of each group of rows
  *
  *   vertcat(vertsplit(x, ...)) = x
  */
  CASADI_CORE_EXPORT std::vector<MX> vertsplit(const MX& x, int incr=1);

  /** \brief Construct a matrix from a list of list of blocks.
  *
  *   blockcat(blocksplit(x,..., ...)) = x
  */
  CASADI_CORE_EXPORT MX blockcat(const std::vector< std::vector<MX > > &v);

  /** \brief  chop up into blocks
  * \brief vert_offset Defines the boundaries of the block cols
  * \brief horz_offset Defines the boundaries of the block rows
  *
  *   blockcat(blocksplit(x,..., ...)) = x
  */
  CASADI_CORE_EXPORT std::vector< std::vector<MX > >
    blocksplit(const MX& x, const std::vector<int>& vert_offset,
               const std::vector<int>& horz_offset);

  /** \brief  chop up into blocks
  * \brief vert_incr Defines the increment for block boundaries in col dimension
  * \brief horz_incr Defines the increment for block boundaries in row dimension
  *
  *   blockcat(blocksplit(x,..., ...)) = x
  */
  CASADI_CORE_EXPORT std::vector< std::vector<MX > >
    blocksplit(const MX& x, int vert_incr = 1, int horz_incr = 1);

#ifndef SWIG
  /** \brief Construct a matrix from a list of list of blocks.*/
  CASADI_CORE_EXPORT MX blockcat(const MX &A, const MX &B, const MX &C, const MX &D);
#endif // SWIG

  /** \brief  concatenate diagonally
  *
  *  diagcat(diagsplit(x, ...)) = x
  */
  CASADI_CORE_EXPORT MX diagcat(const std::vector<MX>& x);

/// \cond INTERNAL
#ifndef SWIG
  /** \brief  split diagonally, retaining matrices
  * \param output_offset1 List of all start locations (row) for each matrix
  * \param output_offset1 List of all start locations (column) for each matrix
  *      the last matrix will run to the end.
  *   diagcat(diagsplit(x, ...)) = x
  */
  CASADI_CORE_EXPORT std::vector<MX> diagsplitNative(const MX& x,
                                                   const std::vector<int>& output_offset1,
                                                   const std::vector<int>& output_offset2);
#endif // SWIG
/// \endcond

  /** \brief Concatenate vertically while vectorizing all arguments */
  CASADI_CORE_EXPORT MX veccat(const std::vector<MX>& comp);

  /** \brief  concatenate vertically while vecing all arguments with vecNZ */
  CASADI_CORE_EXPORT MX vecNZcat(const std::vector<MX>& comp);

#ifndef SWIG
  /** \brief  concatenate vertically, two matrices */
  CASADI_CORE_EXPORT MX horzcat(const MX& a, const MX& b);

  /** \brief  concatenate horizontally, two matrices */
  CASADI_CORE_EXPORT MX vertcat(const MX& a, const MX& b);
#endif // SWIG

  /** \brief  Frobenius norm  */
  CASADI_CORE_EXPORT MX norm_F(const MX &x);

  /** \brief  2-norm  */
  CASADI_CORE_EXPORT MX norm_2(const MX &x);

  /** \brief 1-norm  */
  CASADI_CORE_EXPORT MX norm_1(const MX &x);

  /** \brief Infinity-norm */
  CASADI_CORE_EXPORT MX norm_inf(const MX &x);

  /** \brief Transpose an expression */
  CASADI_CORE_EXPORT MX transpose(const MX &x);

  /** \brief  Take the matrix product of 2 MX objects
  *
  * With optional sp_z you can specify the sparsity of the result
  * A typical use case might be where the product is only constructed to
  * inspect the trace of it. sp_z diagonal will be more efficient then.
  *
  */
  CASADI_CORE_EXPORT MX mul(const MX &x, const MX &y, const Sparsity& sp_z=Sparsity());

  /** \brief  Take the matrix product of n MX objects */
  CASADI_CORE_EXPORT MX mul(const std::vector< MX > &x);

  /** \brief  Take the inner product of two vectors
      Equals
      \code
      trans(x)*y
      \endcode
      with x and y vectors
  */
  CASADI_CORE_EXPORT MX inner_prod(const MX &x, const MX &y);

  /** \brief  Take the outer product of two vectors
      Equals
      \code
      x*trans(y)
      \endcode
      with x and y vectors
  */
  CASADI_CORE_EXPORT MX outer_prod(const MX &x, const MX &y);

  /** \brief Branching on MX nodes
      Ternary operator, "cond ? if_true : if_false"
  */
  CASADI_CORE_EXPORT MX if_else(const MX &cond, const MX &if_true, const MX &if_false);

#ifndef SWIG
  //! \brief Returns a reshaped version of the MX
  CASADI_CORE_EXPORT MX reshape(const MX &x, int nrow, int ncol);
#endif // SWIG

  //! \brief Returns a reshaped version of the MX, dimensions as a vector
  CASADI_CORE_EXPORT MX reshape(const MX &x, std::pair<int, int> rc);

  //! \brief Reshape the MX
  CASADI_CORE_EXPORT MX reshape(const MX &x, const Sparsity& sp);

  /** \brief Returns a vectorized version of the MX
      Same as reshape(x, x.numel(), 1)

      a c
      b d

      turns into

      a
      b
      c
      d

  */
  CASADI_CORE_EXPORT MX vec(const MX &x);

  /** \brief Returns a vectorized version of the MX, preserving only nonzeros
  */
  CASADI_CORE_EXPORT MX vecNZ(const MX &x);

  /** \brief  Unite two matrices no overlapping sparsity */
  CASADI_CORE_EXPORT MX unite(const MX& A, const MX& B);

  /** \brief  Simplify an expression */
  CASADI_CORE_EXPORT void simplify(MX& ex);

  /** \brief Matrix trace */
  CASADI_CORE_EXPORT MX trace(const MX& A);

  /** \brief Repeat matrix A n times vertically and m times horizontally */
  CASADI_CORE_EXPORT MX repmat(const MX &A, int n, int m);

  /** \brief create a clipped view into a matrix
      Create a sparse matrix from a dense matrix A, with sparsity pattern sp
  **/
  //MX clip(const MX& A, const Sparsity& sp);

#ifndef SWIGOCTAVE
  /** \brief  Make the matrix dense if not already */
  CASADI_CORE_EXPORT MX dense(const MX& x);
#endif // SWIGOCTAVE

  /** \brief  Create a parent MX on which all given MX's will depend.

      In some sense, this function is the inverse of

      \param deps  Must all be symbolic matrices.
  */
  CASADI_CORE_EXPORT MX createParent(std::vector<MX> &deps);

  /** \brief  Create a parent MX on which a bunch of MX's (sizes given as argument) will depend
   */
  CASADI_CORE_EXPORT MX createParent(const std::vector<MX> &deps,
                                         std::vector<MX>& SWIG_OUTPUT(children));


  /** \brief  Create a parent MX on which a bunch of MX's (sizes given as argument) will depend
   */
  CASADI_CORE_EXPORT MX createParent(const std::vector<Sparsity> &deps,
                                         std::vector<MX>& SWIG_OUTPUT(children));

  /** Count number of nodes */
  CASADI_CORE_EXPORT int countNodes(const MX& A);

  /** \brief  Get the diagonal of a matrix or construct a diagonal

      When the input is square, the diagonal elements are returned.
      If the input is vector-like, a diagonal matrix is constructed with it.
  */
  CASADI_CORE_EXPORT MX diag(const MX& x);

  /** \brief   Construct a matrix with given blocks on the diagonal */
  CASADI_CORE_EXPORT MX blkdiag(const std::vector<MX> &A);

#ifndef SWIG
  /** \brief   Construct a matrix with given blocks on the diagonal */
  CASADI_CORE_EXPORT MX blkdiag(const MX &A, const MX& B);
#endif // SWIG

  /** \brief Return a col-wise summation of elements */
  CASADI_CORE_EXPORT MX sumCols(const MX &x);

  /** \brief Return a row-wise summation of elements */
  CASADI_CORE_EXPORT MX sumRows(const MX &x);

  /// Return summation of all elements
  CASADI_CORE_EXPORT MX sumAll(const MX &x);

  /** \brief  Evaluate a polynomial with coefficients p in x */
  CASADI_CORE_EXPORT MX polyval(const MX& p, const MX& x);

  /** \brief Get a string representation for a binary MX, using custom arguments */
  CASADI_CORE_EXPORT std::string
    getOperatorRepresentation(const MX& xb, const std::vector<std::string>& args);

  /** \brief  Substitute variable v with expression vdef in an expression ex */
  CASADI_CORE_EXPORT MX substitute(const MX &ex, const MX& v, const MX& vdef);

  /** \brief  Substitute variable var with expression expr in multiple expressions */
  CASADI_CORE_EXPORT std::vector<MX>
    substitute(const std::vector<MX> &ex, const std::vector<MX> &v, const std::vector<MX> &vdef);

  /** \brief  Substitute variable v with expression vdef in an expression ex, preserving nodes */
  CASADI_CORE_EXPORT MX
     graph_substitute(const MX &ex, const std::vector<MX> &v, const std::vector<MX> &vdef);

  /** \brief  Substitute variable var with expression expr in
   * multiple expressions, preserving nodes  */
  CASADI_CORE_EXPORT std::vector<MX>
    graph_substitute(const std::vector<MX> &ex, const std::vector<MX> &v,
                     const std::vector<MX> &vdef);

  /** \brief Inplace substitution
   * Substitute variables v out of the expressions vdef sequentially
   */
#ifndef SWIG
  CASADI_CORE_EXPORT void
    substituteInPlace(const std::vector<MX>& v, std::vector<MX>& vdef, bool reverse=false);
#else // SWIG
  CASADI_CORE_EXPORT void
    substituteInPlace(const std::vector<MX>& v, std::vector<MX>& INOUT, bool reverse=false);
#endif // SWIG

  /** \brief Inplace substitution with piggyback expressions
   * Substitute variables v out of the expressions vdef sequentially,
   * as well as out of a number of other expressions piggyback */
#ifndef SWIG
  CASADI_CORE_EXPORT void
    substituteInPlace(const std::vector<MX>& v, std::vector<MX>& vdef,
                      std::vector<MX>& ex, bool reverse=false);
#else // SWIG
  CASADI_CORE_EXPORT void
    substituteInPlace(const std::vector<MX>& v, std::vector<MX>& INOUT,
                      std::vector<MX>& INOUT, bool reverse=false);
#endif // SWIG

  /** \brief Extract shared subexpressions from an set of expressions */
  CASADI_CORE_EXPORT void extractShared(std::vector<MX>& ex,
                     std::vector<MX>& v, std::vector<MX>& vdef,
                     const std::string& v_prefix="v_", const std::string& v_suffix="");

  /** \brief Print compact, introducing new variables for shared subexpressions */
  CASADI_CORE_EXPORT void printCompact(const MX& ex, std::ostream &stream=std::cout);

  ///@{
  /** \brief Calculate jacobian via source code transformation

      Uses casadi::MXFunction::jac
  */
  CASADI_CORE_EXPORT MX jacobian(const MX &ex, const MX &arg);
  CASADI_CORE_EXPORT MX gradient(const MX &ex, const MX &arg);
  CASADI_CORE_EXPORT MX tangent(const MX &ex, const MX &arg);
  ///@}

  /** \brief Computes the nullspace of a matrix A
  *
  * Finds Z m-by-(m-n) such that AZ = 0
  * with A n-by-m with m > n
  *
  * Assumes A is full rank
  *
  * Inspired by Numerical Methods in Scientific Computing by Ake Bjorck
  */
  CASADI_CORE_EXPORT MX nullspace(const MX& A);

  /** \brief Matrix determinant (experimental) */
  CASADI_CORE_EXPORT MX det(const MX& A);

  /** \brief Matrix inverse (experimental) */
  CASADI_CORE_EXPORT MX inv(const MX& A);

  /** \brief Get all symbols contained in the supplied expression
  * Get all symbols on which the supplied expression depends
  * \see MXFunction::getFree()
  */
  CASADI_CORE_EXPORT std::vector<MX> getSymbols(const MX& e);

  /** \brief Get all symbols contained in the supplied expression
  * Get all symbols on which the supplied expression depends
  * \see MXFunction::getFree()
  */
  CASADI_CORE_EXPORT std::vector<MX> getSymbols(const std::vector<MX>& e);

  /** \brief Check if expression depends on any of the arguments
    The arguments must be symbolic
  */
  CASADI_CORE_EXPORT bool dependsOn(const MX& ex, const std::vector<MX> &arg);

  /** \brief Expand MX graph to SXFunction call
  *
  *  Expand the given expression e, optionally
  *  supplying expressions contained in it at which expansion should stop.
  *
  */
  CASADI_CORE_EXPORT MX
    matrix_expand(const MX& e, const std::vector<MX> &boundary = std::vector<MX>());

  /** \brief Expand MX graph to SXFunction call
  *
  *  Expand the given expression e, optionally
  *  supplying expressions contained in it at which expansion should stop.
  *
  */
  CASADI_CORE_EXPORT std::vector<MX>
    matrix_expand(const std::vector<MX>& e, const std::vector<MX> &boundary = std::vector<MX>());

  /** \brief Kronecker tensor product
  *
  * Creates a block matrix in which each element (i, j) is a_ij*b
  */
  CASADI_CORE_EXPORT MX kron(const MX& a, const MX& b);

  /** \brief Solve a system of equations: A*x = b
  */
  CASADI_CORE_EXPORT MX
  solve(const MX& A, const MX& b, const std::string& lsolver = "symbolicqr",
        const Dictionary& dict = Dictionary());

  /** \brief Computes the Moore-Penrose pseudo-inverse
  *
  * If the matrix A is fat (size1>size2), mul(A, pinv(A)) is unity.
  * If the matrix A is slender (size2<size1), mul(pinv(A), A) is unity.
  *
  */
  CASADI_CORE_EXPORT MX pinv(const MX& A, const std::string& lsolver,
                                 const Dictionary& dict = Dictionary());

/** @}
*/

} // namespace casadi

#ifdef SWIG
// Template instantiations
%template(Pair_MX_MXVector) std::pair<casadi::MX, std::vector<casadi::MX> >;
#endif // SWIG



#endif // CASADI_MX_TOOLS_HPP

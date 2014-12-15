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

  /** \brief  split horizontally, retaining groups of cols
  * \param output_offset List of all start cols for each group
  *      the last col group will run to the end.
  *
  *   horzcat(horzsplit(x, ...)) = x
  */
  inline std::vector<MX> horzsplit(const MX& x, const std::vector<int>& output_offset) {
    return x.zz_horzsplit(output_offset);
  }

  /** \brief  split horizontally, retaining fixed-sized groups of cols
  * \param incr Size of each group of cols
  *
  *  horzcat(horzsplit(x, ...)) = x
  */
  inline std::vector<MX> horzsplit(const MX& x, int incr=1) { return x.zz_horzsplit(incr);}

  /** \brief  split vertically, retaining groups of rows
  * \param output_offset List of all start rows for each group
  *      the last row group will run to the end.
  *
  *   vertcat(vertsplit(x, ...)) = x
  */
  inline std::vector<MX> vertsplit(const MX& x, const std::vector<int>& output_offset) {
    return x.zz_vertsplit(output_offset);
  }

  /** \brief  split vertically, retaining fixed-sized groups of rows
  * \param incr Size of each group of rows
  *
  *   vertcat(vertsplit(x, ...)) = x
  */
  inline std::vector<MX> vertsplit(const MX& x, int incr=1) { return x.zz_vertsplit(incr);}

  /** \brief Construct a matrix from a list of list of blocks.
  *
  *   blockcat(blocksplit(x,..., ...)) = x
  */
  inline MX blockcat(const std::vector< std::vector<MX > > &v) { return MX::zz_blockcat(v);}

  /** \brief  chop up into blocks
  * \brief vert_offset Defines the boundaries of the block cols
  * \brief horz_offset Defines the boundaries of the block rows
  *
  *   blockcat(blocksplit(x,..., ...)) = x
  */
  inline std::vector< std::vector<MX > >
  blocksplit(const MX& x, const std::vector<int>& vert_offset,
             const std::vector<int>& horz_offset) {
    return x.zz_blocksplit(vert_offset, horz_offset);
  }

  /** \brief  chop up into blocks
  * \brief vert_incr Defines the increment for block boundaries in col dimension
  * \brief horz_incr Defines the increment for block boundaries in row dimension
  *
  *   blockcat(blocksplit(x,..., ...)) = x
  */
  inline std::vector< std::vector<MX > >
  blocksplit(const MX& x, int vert_incr = 1, int horz_incr = 1) {
    return x.zz_blocksplit(vert_incr, horz_incr);
  }

#ifndef SWIG
  /** \brief Construct a matrix from a list of list of blocks.*/
  inline MX blockcat(const MX &A, const MX &B, const MX &C, const MX &D) {
    return MX::zz_blockcat(A, B, C, D);
  }
#endif // SWIG

  /** \brief  concatenate diagonally
  *
  *  diagcat(diagsplit(x, ...)) = x
  */
  inline MX diagcat(const std::vector<MX>& x) { return MX::zz_diagcat(x);}

/// \cond INTERNAL
#ifndef SWIG
  /** \brief  split diagonally, retaining matrices
  * \param output_offset1 List of all start locations (row) for each matrix
  * \param output_offset1 List of all start locations (column) for each matrix
  *      the last matrix will run to the end.
  *   diagcat(diagsplit(x, ...)) = x
  */
  inline std::vector<MX> diagsplitNative(const MX& x,
                                         const std::vector<int>& output_offset1,
                                         const std::vector<int>& output_offset2) {
    return x.zz_diagsplitNative(output_offset1, output_offset2);
  }
#endif // SWIG
/// \endcond

  /** \brief Concatenate vertically while vectorizing all arguments */
  inline MX veccat(const std::vector<MX>& comp) { return MX::zz_veccat(comp);}

  /** \brief  concatenate vertically while vecing all arguments with vecNZ */
  inline MX vecNZcat(const std::vector<MX>& comp) { return MX::zz_vecNZcat(comp);}

  /** \brief  Frobenius norm  */
  inline MX norm_F(const MX &x) { return x.zz_norm_F();}

  /** \brief  2-norm  */
  inline MX norm_2(const MX &x) { return x.zz_norm_2();}

  /** \brief 1-norm  */
  inline MX norm_1(const MX &x) { return x.zz_norm_1();}

  /** \brief Infinity-norm */
  inline MX norm_inf(const MX &x) { return x.zz_norm_inf();}

  /** \brief  Take the inner product of two vectors
      Equals
      \code
      trans(x)*y
      \endcode
      with x and y vectors
  */
  inline MX inner_prod(const MX &x, const MX &y) { return x.inner_prod(y);}

  /** \brief  Take the outer product of two vectors
      Equals
      \code
      x*trans(y)
      \endcode
      with x and y vectors
  */
  inline MX outer_prod(const MX &x, const MX &y) { return x.outer_prod(y);}

  /** \brief Branching on MX nodes
      Ternary operator, "cond ? if_true : if_false"
  */
  inline MX if_else(const MX &cond, const MX &if_true, const MX &if_false) {
    return cond.zz_if_else(if_true, if_false);
  }

  /** \brief  Unite two matrices no overlapping sparsity */
  inline MX unite(const MX& A, const MX& B) { return A.zz_unite(B);}

  /** \brief  Simplify an expression */
  inline void simplify(MX& ex) { ex.zz_simplify();}

  /** \brief Repeat matrix A n times vertically and m times horizontally */
  inline MX repmat(const MX &A, int n, int m) { return A.zz_repmat(n, m);}

  /** \brief  Make the matrix dense if not already */
  inline MX dense(const MX& x) { return x.zz_dense();}

  /** \brief  Create a parent MX on which all given MX's will depend.

      In some sense, this function is the inverse of

      \param deps  Must all be symbolic matrices.
  */
  inline MX createParent(std::vector<MX> &deps) { return MX::zz_createParent(deps);}

  /** \brief  Create a parent MX on which a bunch of MX's (sizes given as argument) will depend
   */
  inline MX createParent(const std::vector<MX> &deps,
                         std::vector<MX>& SWIG_OUTPUT(children)) {
    return MX::zz_createParent(deps, children);
  }

  /** \brief  Create a parent MX on which a bunch of MX's (sizes given as argument) will depend
   */
  inline MX createParent(const std::vector<Sparsity> &deps,
                         std::vector<MX>& SWIG_OUTPUT(children)) {
    return MX::zz_createParent(deps, children);
  }

  /** Count number of nodes */
  inline int countNodes(const MX& A) { return A.zz_countNodes();}

  /** \brief  Get the diagonal of a matrix or construct a diagonal

      When the input is square, the diagonal elements are returned.
      If the input is vector-like, a diagonal matrix is constructed with it.
  */
  inline MX diag(const MX& x) { return x.zz_diag();}

  /** \brief Return a col-wise summation of elements */
  inline MX sumCols(const MX &x) { return x.zz_sumCols();}

  /** \brief Return a row-wise summation of elements */
  inline MX sumRows(const MX &x) { return x.zz_sumRows();}

  /// Return summation of all elements
  inline MX sumAll(const MX &x) { return x.zz_sumAll();}

  /** \brief  Evaluate a polynomial with coefficients p in x */
  inline MX polyval(const MX& p, const MX& x) { return p.zz_polyval(x);}

  /** \brief Get a string representation for a binary MX, using custom arguments */
  inline std::string
  getOperatorRepresentation(const MX& xb, const std::vector<std::string>& args) {
    return xb.zz_getOperatorRepresentation(args);
  }

  /** \brief  Substitute variable v with expression vdef in an expression ex */
  inline MX substitute(const MX &ex, const MX& v, const MX& vdef) {
    return ex.zz_substitute(v, vdef);
  }

  /** \brief  Substitute variable var with expression expr in multiple expressions */
  inline std::vector<MX>
  substitute(const std::vector<MX> &ex, const std::vector<MX> &v, const std::vector<MX> &vdef) {
    return MX::zz_substitute(ex, v, vdef);
  }

  /** \brief  Substitute variable v with expression vdef in an expression ex, preserving nodes */
  inline MX
  graph_substitute(const MX &ex, const std::vector<MX> &v, const std::vector<MX> &vdef) {
    return ex.zz_graph_substitute(v, vdef);
  }

  /** \brief  Substitute variable var with expression expr in
   * multiple expressions, preserving nodes  */
  inline std::vector<MX>
  graph_substitute(const std::vector<MX> &ex, const std::vector<MX> &v,
                   const std::vector<MX> &vdef) {
    return MX::zz_graph_substitute(ex, v, vdef);
  }

  /** \brief Inplace substitution
   * Substitute variables v out of the expressions vdef sequentially
   */
  inline void
  substituteInPlace(const std::vector<MX>& v,
                    std::vector<MX>& SWIG_INOUT(vdef), bool reverse=false) {
    return MX::zz_substituteInPlace(v, vdef, reverse);
  }

  /** \brief Inplace substitution with piggyback expressions
   * Substitute variables v out of the expressions vdef sequentially,
   * as well as out of a number of other expressions piggyback */
  inline void
  substituteInPlace(const std::vector<MX>& v,
                    std::vector<MX>& SWIG_INOUT(vdef),
                    std::vector<MX>& SWIG_INOUT(ex), bool reverse=false) {
    return MX::zz_substituteInPlace(v, vdef, ex, reverse);
  }

  /** \brief Extract shared subexpressions from an set of expressions */
  inline void extractShared(std::vector<MX>& ex,
                            std::vector<MX>& v, std::vector<MX>& vdef,
                            const std::string& v_prefix="v_",
                            const std::string& v_suffix="") {
    return MX::zz_extractShared(ex, v, vdef, v_prefix, v_suffix);
  }

  /** \brief Print compact, introducing new variables for shared subexpressions */
  inline void printCompact(const MX& ex, std::ostream &stream=CASADI_COUT) {
    return ex.zz_printCompact(stream);
  }

  ///@{
  /** \brief Calculate jacobian via source code transformation

      Uses casadi::MXFunction::jac
  */
  inline MX jacobian(const MX &ex, const MX &arg) { return ex.zz_jacobian(arg);}
  inline MX gradient(const MX &ex, const MX &arg) { return ex.zz_gradient(arg);}
  inline MX tangent(const MX &ex, const MX &arg) { return ex.zz_tangent(arg);}
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
  inline MX nullspace(const MX& A) { return A.zz_nullspace();}

  /** \brief Get all symbols contained in the supplied expression
  * Get all symbols on which the supplied expression depends
  * \see MXFunction::getFree()
  */
  inline std::vector<MX> getSymbols(const MX& e) { return e.zz_getSymbols();}

  /** \brief Get all symbols contained in the supplied expression
  * Get all symbols on which the supplied expression depends
  * \see MXFunction::getFree()
  */
  inline std::vector<MX> getSymbols(const std::vector<MX>& e) {
    return MX::zz_getSymbols(e);
  }

  /** \brief Check if expression depends on any of the arguments
    The arguments must be symbolic
  */
  inline bool dependsOn(const MX& ex, const std::vector<MX> &arg) { return ex.zz_dependsOn(arg);}

  /** \brief Expand MX graph to SXFunction call
  *
  *  Expand the given expression e, optionally
  *  supplying expressions contained in it at which expansion should stop.
  *
  */
  inline MX
  matrix_expand(const MX& e, const std::vector<MX> &boundary = std::vector<MX>()) {
    return MX::zz_matrix_expand(e, boundary);
  }

  /** \brief Expand MX graph to SXFunction call
  *
  *  Expand the given expression e, optionally
  *  supplying expressions contained in it at which expansion should stop.
  *
  */
  inline std::vector<MX>
  matrix_expand(const std::vector<MX>& e, const std::vector<MX> &boundary = std::vector<MX>()) {
    return MX::zz_matrix_expand(e, boundary);
  }

  /** \brief Kronecker tensor product
  *
  * Creates a block matrix in which each element (i, j) is a_ij*b
  */
  inline MX kron(const MX& a, const MX& b) { return a.zz_kron(b);}

  /** \brief Solve a system of equations: A*x = b
  */
  inline MX
  solve(const MX& A, const MX& b, const std::string& lsolver = "symbolicqr",
        const Dictionary& dict = Dictionary()) { return A.zz_solve(b, lsolver, dict);}

  /** \brief Computes the Moore-Penrose pseudo-inverse
  *
  * If the matrix A is fat (size1>size2), mul(A, pinv(A)) is unity.
  * If the matrix A is slender (size2<size1), mul(pinv(A), A) is unity.
  *
  */
  inline MX pinv(const MX& A, const std::string& lsolver,
                 const Dictionary& dict = Dictionary()) { return A.zz_pinv(lsolver, dict);}
/** @}
*/

} // namespace casadi

#endif // CASADI_MX_TOOLS_HPP

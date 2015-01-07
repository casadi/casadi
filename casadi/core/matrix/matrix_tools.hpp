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


#ifndef CASADI_MATRIX_TOOLS_HPP
#define CASADI_MATRIX_TOOLS_HPP

#include "matrix.hpp"
#include <algorithm>

#include "../options_functionality.hpp"

#include "sparsity_tools.hpp"
#include "generic_matrix_tools.hpp"

namespace casadi {

/**
\ingroup expression_tools
@{
*/
  /** \brief Construct a matrix from a list of list of blocks.
   */
  template<typename DataType>
  Matrix<DataType> blockcat(const std::vector< std::vector<Matrix<DataType> > > &v) {
    return Matrix<DataType>::zz_blockcat(v);
  }

#ifndef SWIG
  /** \brief Construct a matrix from 4 blocks
   */
  template<typename DataType>
  Matrix<DataType> blockcat(const Matrix<DataType> &A, const Matrix<DataType> &B,
                            const Matrix<DataType> &C, const Matrix<DataType> &D) {
    return Matrix<DataType>::zz_blockcat(A, B, C, D);
  }
#endif // SWIG

  /** \brief  split horizontally, retaining groups of columns
   * \param offset List of all start columns for each group
   *      the last column group will run to the end.
   *
   *   horzcat(horzsplit(x, ...)) = x
   */
  template<typename DataType>
  std::vector<Matrix<DataType> > horzsplit(const Matrix<DataType> &v,
                                           const std::vector<int>& offset) {
    return v.zz_horzsplit(offset);
  }

  /** \brief  split horizontally, retaining fixed-sized groups of columns
   * \param incr Size of each group of columns
   *
   *   horzcat(horzsplit(x, ...)) = x
   */
  template<typename DataType>
  std::vector<Matrix<DataType> > horzsplit(const Matrix<DataType> &v, int incr=1) {
    return v.zz_horzsplit(incr);
  }

  /** \brief  split vertically, retaining groups of rows
   * \param output_offset List of all start rows for each group
   *      the last row group will run to the end.
   *
   *   vertcat(vertsplit(x, ...)) = x
   */
  template<typename DataType>
  std::vector<Matrix<DataType> > vertsplit(const Matrix<DataType> &v,
                                           const std::vector<int>& offset) {
    return v.zz_vertsplit(offset);
  }

  /** \brief  split vertically, retaining fixed-sized groups of rows
   * \param incr Size of each group of rows
   *
   *   vertcat(vertsplit(x, ...)) = x
   */
  template<typename DataType>
  std::vector<Matrix<DataType> > vertsplit(const Matrix<DataType> &v, int incr=1) {
    return v.zz_vertsplit(incr);
  }

  /** \brief  chop up into blocks
   * \param vert_offset Defines the boundaries of the block rows
   * \param horz_offset Defines the boundaries of the block columns
   *
   *   blockcat(blocksplit(x,..., ...)) = x
   */
  template<typename DataType>
  std::vector< std::vector< Matrix<DataType> > > blocksplit(const Matrix<DataType>& x,
                                                            const std::vector<int>& vert_offset,
                                                            const std::vector<int>& horz_offset) {
    return x.zz_blocksplit(vert_offset, horz_offset);
  }

  /** \brief  chop up into blocks
   * \param vert_incr Defines the increment for block boundaries in row dimension
   * \param horz_incr Defines the increment for block boundaries in column dimension
   *
   *   blockcat(blocksplit(x,..., ...)) = x
   */
  template<typename DataType>
  std::vector< std::vector< Matrix<DataType> > > blocksplit(const Matrix<DataType>& x,
                                                            int vert_incr=1,
                                                            int horz_incr=1) {
    return x.zz_blocksplit(vert_incr, horz_incr);
  }

/// \cond INTERNAL
#ifndef SWIG
  /** \brief  split diagonally, retaining matrices
  * \param output_offset1 List of all start locations (row) for each matrix
  *      the last matrix will run to the end.
  *   diagcat(diagsplit(x, ...)) = x
  */
  template<typename DataType>
  std::vector< Matrix<DataType> > diagsplitNative(const Matrix<DataType>& x,
                                                  const std::vector<int>& output_offset1,
                                                  const std::vector<int>& output_offset2) {
    return x.zz_diagsplitNative(output_offset1, output_offset2);
  }
#endif // SWIG
/// \endcond

  template<typename DataType>
  /** \brief  concatenate vertically while vectorizing all arguments with vec */
  Matrix<DataType> veccat(const std::vector< Matrix<DataType> >& x) {
    return Matrix<DataType>::zz_veccat(x);
  }

  template<typename DataType>
  /** \brief  concatenate vertically while vectorizing all arguments with vecNZ */
  Matrix<DataType> vecNZcat(const std::vector< Matrix<DataType> >& x) {
    return Matrix<DataType>::zz_vecNZcat(x);
  }

  /** \brief Inner product of two matrices
      Equals
      \code
      sumAll(x*y)
      \endcode
      with x and y matrices of the same dimension
  */
  // inner product
  template<typename DataType>
  Matrix<DataType> inner_prod(const Matrix<DataType> &x, const Matrix<DataType> &y) {
    return x.zz_inner_prod(y);
  }

  /** \brief Outer product of two vectors
      Equals
      \code
      x*transpose(y)
      \endcode
      with x and y vectors
  */
  template<typename DataType>
  Matrix<DataType> outer_prod(const Matrix<DataType> &x, const Matrix<DataType> &y) {
    return x.zz_outer_prod(y);
  }

  /** \brief  QR factorization using the modified Gram-Schmidt algorithm
   * More stable than the classical Gram-Schmidt, but may break down if the rows of A
   * are nearly linearly dependent
   * See J. Demmel: Applied Numerical Linear Algebra (algorithm 3.1.).
   * Note that in SWIG, Q and R are returned by value. */
  template<typename DataType>
  void qr(const Matrix<DataType>& A,
          Matrix<DataType>& SWIG_OUTPUT(Q),
          Matrix<DataType>& SWIG_OUTPUT(R)) {
    return A.zz_qr(Q, R);
  }

  /** \brief Computes the nullspace of a matrix A
  *
  * Finds Z m-by-(m-n) such that AZ = 0
  * with A n-by-m with m > n
  *
  * Assumes A is full rank
  *
  * Inspired by Numerical Methods in Scientific Computing by Ake Bjorck
  */
  template<typename DataType>
  Matrix<DataType> nullspace(const Matrix<DataType>& A) { return A.zz_nullspace();}

  /** \brief  Solve a system of equations: A*x = b
      The solve routine works similar to Matlab's backslash when A is square and nonsingular.
      The algorithm used is the following:
      1. A simple forward or backward substitution if A is upper or lower triangular
      2. If the linear system is at most 3-by-3, form the inverse via minor expansion and multiply
      3. Permute the variables and equations as to get a (structurally) nonzero diagonal,
      then perform a QR factorization without pivoting and solve the factorized system.

      Note 1: If there are entries of the linear system known to be zero, these will be removed.
      Elements that are very small, or will evaluate to be zero, can still cause numerical errors,
      due to the lack of pivoting (which is not possible since cannot compare the size of entries)

      Note 2: When permuting the linear system, a BLT (block lower triangular) transformation is
      formed. Only the permutation part of this is however used. An improvement would be to solve
      block-by-block if there are multiple BLT blocks.

  */
  template<typename DataType>
  Matrix<DataType> solve(const Matrix<DataType>& A, const Matrix<DataType>& b) {
    return A.zz_solve(b);
  }

  /** \brief Computes the Moore-Penrose pseudo-inverse
  *
  * If the matrix A is fat (size2>size1), mul(A, pinv(A)) is unity.
  * If the matrix A is slender (size1<size2), mul(pinv(A), A) is unity.
  *
  */
  template<typename DataType>
  Matrix<DataType> pinv(const Matrix<DataType>& A) { return A.zz_pinv();}

  /** \brief Solve a system of equations: A*x = b
  */
  CASADI_EXPORT Matrix<double> solve(const Matrix<double>& A, const Matrix<double>& b,
                                          const std::string& lsolver,
                                          const Dictionary& dict = Dictionary());


  /** \brief Computes the Moore-Penrose pseudo-inverse
  *
  * If the matrix A is fat (size1>size2), mul(A, pinv(A)) is unity.
  * If the matrix A is slender (size2<size1), mul(pinv(A), A) is unity.
  *
  */
  CASADI_EXPORT Matrix<double> pinv(const Matrix<double>& A, const std::string& lsolver,
                                    const Dictionary& dict = Dictionary());


  /** Inf-norm of a Matrix-matrix product, no memory allocation
  *   mul(x, y)
  *
  * \param Dwork  A double work vector that you must allocate
  *               Minimum size: y.size1()
  * \param Iwork  A integer work vector that you must allocate
  *               Minimum size: y.size1()+x.size2()+1
  */
  CASADI_EXPORT double norm_inf_mul_nn(const Matrix<double> &x,
                              const Matrix<double> &y,
                              std::vector<double>& Dwork,
                              std::vector<int>& Iwork);

  /** 0-norm (nonzero count) of a Matrix-matrix product, no memory allocation
  *   mul(x, y)
  *
  * \param Bwork  A boolean work vector that you must allocate
  *               Minimum size: y.size1()
  * \param Iwork  A integer work vector that you must allocate
  *               Minimum size: y.size1()+x.size2()+1
  */
  template<typename DataType>
  int norm_0_mul_nn(const Matrix<DataType> &x,
                    const Matrix<DataType> &y,
                    std::vector<bool>& Bwork,
                    std::vector<int>& Iwork) {
    return x.zz_norm_0_mul_nn(y, Bwork, Iwork);
  }

  /** \brief Kronecker tensor product
  *
  * Creates a block matrix in which each element (i, j) is a_ij*b
  */
  template<typename DataType>
  Matrix<DataType> kron(const Matrix<DataType>& a, const Matrix<DataType>& b) {
    return a.zz_kron(b);
  }

  /** \brief  Frobenius norm  */
  template<typename DataType>
  Matrix<DataType> norm_F(const Matrix<DataType> &x) { return x.zz_norm_F();}

  /** \brief  2-norm  */
  template<typename DataType>
  Matrix<DataType> norm_2(const Matrix<DataType> &x) { return x.zz_norm_2();}

  /** \brief 1-norm  */
  template<typename DataType>
  Matrix<DataType> norm_1(const Matrix<DataType> &x) { return x.zz_norm_1();}

  /** \brief Infinity-norm */
  template<typename DataType>
  Matrix<DataType> norm_inf(const Matrix<DataType> &x) { return x.zz_norm_inf();}

  /// Return summation of all elements
  template<typename DataType>
  Matrix<DataType> sumAll(const Matrix<DataType> &x) { return x.zz_sumAll();}

  /** \brief Return a col-wise summation of elements */
  template<typename DataType>
  Matrix<DataType> sumCols(const Matrix<DataType> &x) { return x.zz_sumCols();}

  /** \brief Return a row-wise summation of elements */
  template<typename DataType>
  Matrix<DataType> sumRows(const Matrix<DataType> &x) { return x.zz_sumRows();}

  /// Returns true only if every element in the matrix is true
  template<typename DataType>
  Matrix<DataType> all(const Matrix<DataType> &x) { return x.zz_all();}

  /// Returns true if any element in the matrix is true
  template<typename DataType>
  Matrix<DataType> any(const Matrix<DataType> &x) { return x.zz_any();}

  /** \brief Repeat matrix A n times vertically and m times horizontally */
  template<typename DataType>
  Matrix<DataType> repmat(const Matrix<DataType> &A, int n, int m) { return A.zz_repmat(n, m);}

  /** \brief  Evaluate a polynomial with coefficients p in x */
  template<typename DataType>
  Matrix<DataType> polyval(const Matrix<DataType>& p, const Matrix<DataType>& x) {
    return p.zz_polyval(x);
  }

  /** \brief   Get the diagonal of a matrix or construct a diagonal
      When the input is square, the diagonal elements are returned.
      If the input is vector-like, a diagonal matrix is constructed with it. */
  template<typename DataType>
  Matrix<DataType> diag(const Matrix<DataType> &A) { return A.zz_diag();}

  /** \brief  Unite two matrices no overlapping sparsity */
  template<typename DataType>
  Matrix<DataType> unite(const Matrix<DataType>& A, const Matrix<DataType>& B) {
    return A.zz_unite(B);
  }

  /** \brief  Make a matrix dense */
  template<typename DataType>
  Matrix<DataType> dense(const Matrix<DataType>& A);

  /** \brief  Make a matrix sparse by removing numerical zeros*/
  template<typename DataType>
  Matrix<DataType> sparse(const Matrix<DataType>& A, double tol=0);

  /// same as: res += mul(A, v)
  template<typename DataType>
  void addMultiple(const Matrix<DataType>& A, const std::vector<DataType>& v,
                   std::vector<DataType>& res, bool trans_A=false) {
    return A.zz_addMultiple(v, res, trans_A);
  }

  /// \cond INTERNAL
  /// Get a pointer to the data contained in the vector
  template<typename DataType>
  DataType* getPtr(Matrix<DataType> &v);

  /// Get a pointer to the data contained in the vector
  template<typename DataType>
  const DataType* getPtr(const Matrix<DataType> &v);
  /// \endcond

  /** \brief Create a new matrix with a given sparsity pattern but with the
   * nonzeros taken from an existing matrix */
  template<typename DataType>
  Matrix<DataType> project(const Matrix<DataType>& A, const Sparsity& sp) {
    return A.zz_project(sp);
  }

  /// Obtain the structural rank of a sparsity-pattern
  template<typename DataType>
  int sprank(const Matrix<DataType>& A) { return A.zz_sprank();}
/*
@}
*/

} // namespace casadi

// Global namespace

#ifndef SWIG
#include <iterator>

namespace casadi {
  // Implementations
  template<typename DataType>
  DataType* getPtr(Matrix<DataType> &v) {
    if (v.isEmpty())
      return 0;
    else
      return &v.front();
  }

  template<typename DataType>
  const DataType* getPtr(const Matrix<DataType> &v) {
    if (v.isEmpty())
      return 0;
    else
      return &v.front();
  }

  template<typename DataType>
  Matrix<DataType> dense(const Matrix<DataType>& A) {
    Matrix<DataType> ret = A;
    ret.densify();
    return ret;
  }

  template<typename DataType>
  Matrix<DataType> sparse(const Matrix<DataType>& A, double tol) {
    Matrix<DataType> ret = A;
    ret.sparsify(tol);
    return ret;
  }

} // namespace casadi

#endif //SWIG

#endif // CASADI_MATRIX_TOOLS_HPP

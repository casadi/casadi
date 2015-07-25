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


/** \brief Calculate quadratic form X^T A X*/
inline SWIG_FRIEND MatType quad_form(const MatType &X, const MatType &A) {
  return X.zz_quad_form(A);
}

/** \brief Calculate quadratic form X^T X*/
inline SWIG_FRIEND MatType quad_form(const MatType &X) {
  return X.zz_quad_form();
}

/** \brief Calculate some of squares: sum_ij X_ij^2  */
inline SWIG_FRIEND MatType sum_square(const MatType &X) {
  return X.zz_sum_square();
}

/** \brief Matlab's \c linspace command
 */
inline SWIG_FRIEND MatType linspace(const MatType &a, const MatType &b, int nsteps) {
  return a.zz_linspace(b, nsteps);
}

/** \brief Matlab's \c cross command
 */
inline SWIG_FRIEND MatType cross(const MatType &a, const MatType &b, int dim = -1) {
  return a.zz_cross(b, dim);
}

/** \brief Matrix determinant (experimental) */
inline SWIG_FRIEND MatType det(const MatType& A) { return A.zz_det();}

/** \brief Matrix inverse (experimental) */
inline SWIG_FRIEND MatType inv(const MatType& A) { return A.zz_inv();}

/** \brief Matrix trace */
inline SWIG_FRIEND MatType trace(const MatType& a) { return a.zz_trace();}

/** \brief Convert a lower triangular matrix to a symmetric one
 */
inline SWIG_FRIEND MatType tril2symm(const MatType &a) { return a.zz_tril2symm();}

/** \brief Convert a upper triangular matrix to a symmetric one
 */
inline SWIG_FRIEND MatType triu2symm(const MatType &a) { return a.zz_triu2symm();}

/** \brief  Frobenius norm  */
inline SWIG_FRIEND MatType norm_F(const MatType &x) { return x.zz_norm_F();}

/** \brief  2-norm  */
inline SWIG_FRIEND MatType norm_2(const MatType &x) { return x.zz_norm_2();}

/** \brief 1-norm  */
inline SWIG_FRIEND MatType norm_1(const MatType &x) { return x.zz_norm_1();}

/** \brief Infinity-norm */
inline SWIG_FRIEND MatType norm_inf(const MatType &x) { return x.zz_norm_inf();}

/** \brief Return a col-wise summation of elements */
inline SWIG_FRIEND MatType sumCols(const MatType &x) { return x.zz_sumCols();}

/** \brief Return a row-wise summation of elements */
inline SWIG_FRIEND MatType sumRows(const MatType &x) { return x.zz_sumRows();}

/** \brief Inner product of two matrices
    with x and y matrices of the same dimension
*/
inline SWIG_FRIEND MatType inner_prod(const MatType &x, const MatType &y) {
  return x.zz_inner_prod(y);
}

/** \brief  Take the outer product of two vectors
    Equals
    \code
    x*y.T()
    \endcode
    with x and y vectors
*/
inline SWIG_FRIEND MatType outer_prod(const MatType &x, const MatType &y) {
  return x.zz_outer_prod(y);
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
inline SWIG_FRIEND MatType nullspace(const MatType& A) { return A.zz_nullspace();}

/** \brief  Evaluate a polynomial with coefficients p in x */
inline SWIG_FRIEND MatType polyval(const MatType& p, const MatType& x) {
  return p.zz_polyval(x);
}

/** \brief   Get the diagonal of a matrix or construct a diagonal
    When the input is square, the diagonal elements are returned.
    If the input is vector-like, a diagonal matrix is constructed with it. */
inline SWIG_FRIEND MatType diag(const MatType &A) { return A.zz_diag();}

/** \brief  Unite two matrices no overlapping sparsity */
inline SWIG_FRIEND MatType unite(const MatType& A, const MatType& B) {
  return A.zz_unite(B);
}

/** \brief  Make the matrix dense if not already */
inline SWIG_FRIEND MatType densify(const MatType& x) { return x.zz_densify();}

/** \brief Create a new matrix with a given sparsity pattern but with the
 * nonzeros taken from an existing matrix */
inline SWIG_FRIEND MatType project(const MatType& A, const Sparsity& sp,
                              bool intersect=false) {
  return A.zz_project(sp, intersect);
}

/** \brief Check if expression depends on the argument
    The argument must be symbolic
*/
//inline SWIG_FRIEND bool dependsOn(const MatType& f, const MatType &arg) {
//return f.zz_dependsOn(arg);
//}

/** \brief Branching on MX nodes
    Ternary operator, "cond ? if_true : if_false"
*/
inline SWIG_FRIEND MatType if_else(const MatType &cond, const MatType &if_true,
                              const MatType &if_false, bool short_circuit=true) {
  return cond.zz_if_else(if_true, if_false, short_circuit);
}

/** \brief Create a switch
 *
 * If the condition \param ind evaluates to the integer k, where 0<=k<f.size(),
 * then x[k] will be returned, otherwise \param x_default will be returned.
 */
inline SWIG_FRIEND MatType conditional(const MatType& ind, const std::vector<MatType> &x,
                                  const MatType &x_default, bool short_circuit=true) {
  return ind.zz_conditional(x, x_default, short_circuit);
}

/** \brief Check if expression depends on the argument
    The argument must be symbolic
*/
inline SWIG_FRIEND bool dependsOn(const MatType& f, const MatType &arg) {
  return f.zz_dependsOn(arg);
}

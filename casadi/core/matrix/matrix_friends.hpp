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


#ifndef CASADI_MATRIX_FRIENDS_HPP
#define CASADI_MATRIX_FRIENDS_HPP

/** \brief Matrix adjoint */
inline SWIG_FRIEND Matrix<DataType> adj(const Matrix<DataType>& A) { return A.zz_adj();}

/** \brief Get the (i,j) minor matrix */
inline SWIG_FRIEND Matrix<DataType> getMinor(const Matrix<DataType> &x, int i, int j) {
  return x.zz_getMinor(i, j);
}

/** \brief Get the (i,j) cofactor matrix */
inline SWIG_FRIEND Matrix<DataType> cofactor(const Matrix<DataType> &x, int i, int j) {
  return x.zz_cofactor(i, j);
}

/** \brief  QR factorization using the modified Gram-Schmidt algorithm
 * More stable than the classical Gram-Schmidt, but may break down if the rows of A
 * are nearly linearly dependent
 * See J. Demmel: Applied Numerical Linear Algebra (algorithm 3.1.).
 * Note that in SWIG, Q and R are returned by value. */
inline SWIG_FRIEND void qr(const Matrix<DataType>& A, Matrix<DataType>& Q, Matrix<DataType>& R) {
  return A.zz_qr(Q, R);
}

/** \brief Obtain a Cholesky factorisation of a matrix
 *
 *  Returns an upper triangular R such that R'R = A.
 *  Matrix A must be positive definite.
 *
 * At the moment, the algorithm is dense (Cholesky-Banachiewicz).
 * There is an open ticket #1212 to make it sparse.
 */
inline SWIG_FRIEND Matrix<DataType> chol(const Matrix<DataType>& A) {
  return A.zz_chol();
}

/// Returns true only if every element in the matrix is true
inline SWIG_FRIEND Matrix<DataType> all(const Matrix<DataType> &x) { return x.zz_all();}

/// Returns true if any element in the matrix is true
inline SWIG_FRIEND Matrix<DataType> any(const Matrix<DataType> &x) { return x.zz_any();}

/** Inf-norm of a Matrix-Matrix product */
inline SWIG_FRIEND Matrix<DataType>
norm_inf_mul(const Matrix<DataType> &x, const Matrix<DataType> &y) {
  return x.zz_norm_inf_mul(y);
}

/** \brief  Make a matrix sparse by removing numerical zeros*/
inline SWIG_FRIEND Matrix<DataType> sparsify(const Matrix<DataType>& A, double tol=0) {
  return A.zz_sparsify(tol);
}

#endif // CASADI_MATRIX_FRIENDS_HPP

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


#ifndef CASADI_GENERIC_MATRIX_TOOLS_HPP
#define CASADI_GENERIC_MATRIX_TOOLS_HPP

#include "slice.hpp"
#include "submatrix.hpp"
#include "nonzeros.hpp"
#include "sparsity.hpp"
#include "../casadi_math.hpp"
#include "../casadi_exception.hpp"

namespace casadi {
  // Helper function
  template<typename MatType>
  const MatType& mat(const GenericMatrix<MatType>& x) {
    return static_cast<const MatType&>(x);
  }

/**
\ingroup expression_tools
@{
*/
  /** \brief Calculate quadratic form X^T A X*/
  template<typename MatType>
  MatType quad_form(const GenericMatrix<MatType> &X, const GenericMatrix<MatType> &A);

  /** \brief Calculate quadratic form X^T X*/
  template<typename MatType>
  MatType quad_form(const GenericMatrix<MatType> &X);

  /** \brief Calculate some of squares: sum_ij X_ij^2  */
  template<typename MatType>
  MatType sum_square(const GenericMatrix<MatType> &X);

  /** \brief Matlab's \c linspace command
   */
  template<typename MatType>
  MatType linspace(const GenericMatrix<MatType> &a, const GenericMatrix<MatType> &b, int nsteps);

  /** \brief Matlab's \c cross command
   */
  template<typename MatType>
  MatType cross(const GenericMatrix<MatType> &a, const GenericMatrix<MatType> &b, int dim = -1);

  /** \brief Convert a lower triangular matrix to a symmetric one
   */
  template<typename MatType>
  MatType tril2symm(const GenericMatrix<MatType> &a);

  /** \brief Convert a upper triangular matrix to a symmetric one
   */
  template<typename MatType>
  MatType triu2symm(const GenericMatrix<MatType> &a);

  /** \brief Get the upper triangular part of a matrix
   */
  template<typename MatType>
  MatType triu(const GenericMatrix<MatType> &a);

  /** \brief Get the lower triangular part of a matrix
   */
  template<typename MatType>
  MatType tril(const GenericMatrix<MatType> &a);

  /** \brief Check if two expressions are equal, assuming that they are comparable */
  template<typename MatType>
  bool isEqual(const GenericMatrix<MatType>& x, const GenericMatrix<MatType>& y) {
    return mat(x).isEqual(mat(y));
  }

  /** \brief  split diagonally, retaining groups of square matrices
  * \param incr Size of each matrix
  *
  *  diagsplit(diagsplit(x, ...)) = x
  */
  template<typename MatType>
  std::vector< MatType > diagsplit(const GenericMatrix<MatType>& x, int incr=1);

  /** \brief  split diagonally, retaining fixed-sized matrices
  * \param incr1 Row dimension of each matrix
  * \param incr2 Column dimension of each matrix
  *
  *  diagsplit(diagsplit(x, ...)) = x
  */
  template<typename MatType>
  std::vector< MatType > diagsplit(const GenericMatrix<MatType>& x, int incr1, int incr2);

  /** \brief  split diagonally, retaining square matrices
  * \param output_offset List of all start locations for each group
  *      the last matrix will run to the end.
  *
  *   diagcat(diagsplit(x, ...)) = x
  */
  template<typename MatType>
  std::vector< MatType > diagsplit(const GenericMatrix<MatType>& x,
                                   const std::vector<int>& output_offset);

  /** \brief  split diagonally, retaining square matrices
  * \param output_offset1 List of all start locations (row) for each group
  *      the last matrix will run to the end.
  * \param output_offset2 List of all start locations (row) for each group
  *      the last matrix will run to the end.
  *
  *   diagcat(diagsplit(x, ...)) = x
  */
  template<typename MatType>
  std::vector< MatType > diagsplit(const GenericMatrix<MatType>& x,
                                   const std::vector<int>& output_offset1,
                                   const std::vector<int>& output_offset2);

  /** \brief Matrix determinant (experimental) */
  template<typename MatType>
  MatType det(const GenericMatrix<MatType>& A) { return mat(A).zz_det();}

  /** \brief Matrix inverse (experimental) */
  template<typename MatType>
  MatType inv(const GenericMatrix<MatType>& A) { return mat(A).zz_inv();}

  /** \brief Matrix adjoint */
  template<typename MatType>
  MatType adj(const GenericMatrix<MatType>& A) { return mat(A).zz_adj();}

  /** \brief Get the (i,j) minor matrix */
  template<typename MatType>
  MatType getMinor(const GenericMatrix<MatType> &x, int i, int j) {
    return mat(x).zz_getMinor(i, j);
  }

  /** \brief Get the (i,j) cofactor matrix */
  template<typename MatType>
  MatType cofactor(const GenericMatrix<MatType> &x, int i, int j) {
    return mat(x).zz_cofactor(i, j);
  }

  //! \brief Returns a reshaped version of the matrix
  template<typename MatType>
  MatType reshape(const GenericMatrix<MatType>& a, int nrow, int ncol) {
    return mat(a).zz_reshape(nrow, ncol);
  }

  //! \brief Returns a reshaped version of the matrix, dimensions as a vector
  template<typename MatType>
  MatType reshape(const GenericMatrix<MatType>& a, std::pair<int, int> rc) {
    return mat(a).zz_reshape(rc);
  }

  //! \brief Reshape the matrix
  template<typename MatType>
  MatType reshape(const GenericMatrix<MatType>& a,
                  const Sparsity& sp) { return mat(a).zz_reshape(sp);}

  /** \brief  make a vector
      Reshapes/vectorizes the matrix such that the shape becomes (expr.numel(), 1).
      Columns are stacked on top of each other.
      Same as reshape(expr, expr.numel(), 1)

      a c \n
      b d \n

      turns into

      a \n
      b \n
      c \n
      d \n

  */
  template<typename MatType>
  MatType vec(const GenericMatrix<MatType>& a) { return mat(a).zz_vec();}

  /** \brief Returns a flattened version of the matrix, preserving only nonzeros
   */
  template<typename MatType>
  MatType vecNZ(const GenericMatrix<MatType>& a) { return mat(a).zz_vecNZ();}

  /** \brief Matrix trace */
  template<typename MatType>
  MatType trace(const GenericMatrix<MatType>& a) { return mat(a).zz_trace();}

#ifndef SWIG
  template<typename MatType>
  MatType linspace(const GenericMatrix<MatType> &a_, const GenericMatrix<MatType> &b_, int nsteps) {
    const MatType& a = mat(a_);
    const MatType& b = mat(b_);
    std::vector<MatType> ret(nsteps);
    ret[0] = a;
    MatType step = (b-a)/(nsteps-1);

    for (int i=1; i<nsteps-1; ++i)
      ret[i] = ret[i-1] + step;

    ret[nsteps-1] = b;
    return vertcat(ret);
  }
#endif // SWIG

#ifndef SWIG
  // Implementations

  template<typename MatType>
  MatType cross(const GenericMatrix<MatType> &a, const GenericMatrix<MatType> &b, int dim) {
    casadi_assert_message(a.size1()==b.size1() && a.size2()==b.size2(),
                          "cross(a, b): Inconsistent dimensions. Dimension of a ("
                          << a.dimString() << " ) must equal that of b ("
                          << b.dimString() << ").");

    casadi_assert_message(a.size1()==3 || a.size2()==3,
                          "cross(a, b): One of the dimensions of a should have length 3, but got "
                          << a.dimString() << ".");
    casadi_assert_message(dim==-1 || dim==1 || dim==2,
                          "cross(a, b, dim): Dim must be 1, 2 or -1 (automatic).");

    std::vector<MatType> ret(3);

    bool t = a.size1()==3;

    if (dim==1) t = true;
    if (dim==2) t = false;

    MatType a1 = t ? a(0, ALL) : a(ALL, 0);
    MatType a2 = t ? a(1, ALL) : a(ALL, 1);
    MatType a3 = t ? a(2, ALL) : a(ALL, 2);

    MatType b1 = t ? b(0, ALL) : b(ALL, 0);
    MatType b2 = t ? b(1, ALL) : b(ALL, 1);
    MatType b3 = t ? b(2, ALL) : b(ALL, 2);

    ret[0] = a2*b3-a3*b2;
    ret[1] = a3*b1-a1*b3;
    ret[2] = a1*b2-a2*b1;

    return t ? vertcat(ret) : horzcat(ret);

  }

  template<typename MatType>
  MatType tril2symm(const GenericMatrix<MatType> &a_) {
    const MatType& a = mat(a_);
    casadi_assert_message(a.isSquare(),
                          "Shape error in tril2symm. Expecting square shape but got "
                          << a.dimString());
    casadi_assert_message(a.sizeU()-a.sizeD()==0,
                          "Sparsity error in tril2symm. Found above-diagonal entries in argument: "
                          << a.dimString());
    return a +  a.T() - diag(diag(a));
  }


  template<typename MatType>
  MatType triu2symm(const GenericMatrix<MatType> &a_) {
    const MatType& a = mat(a_);
    casadi_assert_message(a.isSquare(),
                          "Shape error in triu2symm. Expecting square shape but got "
                          << a.dimString());
    casadi_assert_message(a.sizeL()-a.sizeD()==0,
                          "Sparsity error in triu2symm. Found below-diagonal entries in argument: "
                          << a.dimString());
    return a + a.T() - diag(diag(a));
  }

  template<typename MatType>
  MatType triu(const GenericMatrix<MatType> &a_) {
    const MatType& a = mat(a_);
    return a.setSparse(a.sparsity().getTriu());
  }

  template<typename MatType>
  MatType tril(const GenericMatrix<MatType> &a_) {
    const MatType& a = mat(a_);
    return a.setSparse(a.sparsity().getTril());
  }

  template<typename MatType>
  std::vector< MatType > diagsplit(const GenericMatrix<MatType>& x_, int incr) {
    const MatType& x = mat(x_);
    casadi_assert(incr>=1);
    casadi_assert_message(x.isSquare(),
      "diagsplit(x,incr)::input must be square but got " << x.dimString()  << ".");
    std::vector<int> offset2 = range(0, x.size2(), incr);
    offset2.push_back(x.size2());
    return diagsplit(x, offset2);
  }

  template<typename MatType>
  std::vector< MatType > diagsplit(const GenericMatrix<MatType>& x_, int incr1, int incr2) {
    const MatType& x = mat(x_);
    casadi_assert(incr1>=1);
    casadi_assert(incr2>=1);
    std::vector<int> offset1 = range(0, x.size1(), incr1);
    offset1.push_back(x.size1());
    std::vector<int> offset2 = range(0, x.size2(), incr2);
    offset2.push_back(x.size2());
    return diagsplit(x, offset1, offset2);
  }

  template<typename MatType>
  std::vector< MatType > diagsplit(const GenericMatrix<MatType>& x_,
    const std::vector<int>& output_offset1,
    const std::vector<int>& output_offset2) {
    const MatType& x = mat(x_);
    return diagsplitNative(x, output_offset1, output_offset2);
  }

  template<typename MatType>
  std::vector< MatType > diagsplit(const GenericMatrix<MatType>& x_,
      const std::vector<int>& output_offset) {
    const MatType& x = mat(x_);
    casadi_assert_message(x.isSquare(),
      "diagsplit(x,incr)::input must be square but got " << x.dimString()  << ".");
    return diagsplit(x, output_offset, output_offset);
  }

  template<typename MatType>
  MatType quad_form(const GenericMatrix<MatType> &X_, const GenericMatrix<MatType> &A_) {
    const MatType& X = mat(X_);
    const MatType& A = mat(A_);
    return mul(X.T(), mul(A, X));
  }

  template<typename MatType>
  MatType quad_form(const GenericMatrix<MatType> &X_) {
    const MatType& X = mat(X_);
    return mul(X.T(), X);
  }

  template<typename MatType>
  MatType sum_square(const GenericMatrix<MatType> &X_) {
    const MatType& X = mat(X_);
    return sumAll(X*X);
  }
#endif // SWIG

/**
* @}
*/

} // namespace casadi

#endif // CASADI_GENERIC_MATRIX_TOOLS_HPP

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
  /** \brief Convert a lower triangular matrix to a symmetric one
   */
  template<typename MatType>
  MatType tril2symm(const GenericMatrix<MatType> &a);

  /** \brief Convert a upper triangular matrix to a symmetric one
   */
  template<typename MatType>
  MatType triu2symm(const GenericMatrix<MatType> &a);

  /** \brief Check if two expressions are equal, assuming that they are comparable */
  template<typename MatType>
  bool isEqual(const GenericMatrix<MatType>& x, const GenericMatrix<MatType>& y) {
    return mat(x).isEqual(mat(y));
  }

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

  /** \brief Matrix trace */
  template<typename MatType>
  MatType trace(const GenericMatrix<MatType>& a) { return mat(a).zz_trace();}

#ifndef SWIG
  // Implementations
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
#endif // SWIG

/**
* @}
*/

} // namespace casadi

#endif // CASADI_GENERIC_MATRIX_TOOLS_HPP

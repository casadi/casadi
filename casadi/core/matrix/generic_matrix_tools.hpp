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

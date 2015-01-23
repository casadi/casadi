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

namespace casadi {

/**
\ingroup expression_tools
@{
*/
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
/*
@}
*/

} // namespace casadi

#endif // CASADI_MATRIX_TOOLS_HPP

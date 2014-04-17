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

#ifndef CSPARSE_CHOLESKY_HPP
#define CSPARSE_CHOLESKY_HPP

#include "casadi/symbolic/function/linear_solver.hpp"
#include <casadi/interfaces/csparse/casadi_csparse_interface_export.h>

namespace casadi{


/** \brief  Forward declaration of internal class */
class CSparseCholeskyInternal;

/** \brief  LinearSolver with CSparseCholesky Interface
*
 @copydoc LinearSolver_doc
*
* CSparseCholesky is an casadi::Function mapping from 2 inputs
* [ A (matrix),b (vector)] to one output [x (vector)].
*
*  A = LL'
*    Ax = b
*    LL'x = b
*    L'x = L^-1 b
*
* The usual procedure to use CSparseCholesky is: \n
*  -# init()
*  -# set the first input (A)
*  -# prepare()
*  -# set the second input (b)
*  -# solve()
*  -# Repeat steps 4 and 5 to work with other b vectors.
*
* The method evaluate() combines the prepare() and solve()
* step and is therefore more expensive if A is invariant.
*
*/
class CASADI_CSPARSE_INTERFACE_EXPORT CSparseCholesky : public LinearSolver{
public:

  /// Default (empty) constructor
  CSparseCholesky();

  /// Create a linear solver given a sparsity pattern
  CSparseCholesky(const Sparsity& sp, int nrhs=1);

  /** \brief  Access internal functions and data members */
  CSparseCholeskyInternal* operator->();

  /** \brief  Access internal functions and data members */
  const CSparseCholeskyInternal* operator->() const;

  /// Check if the node is pointing to the right type of object
  virtual bool checkNode() const;

  /// Solve the system of equations Lx = b
  void solveL(double* x, int nrhs, bool transpose);

  /// Obtain a symbolic Cholesky factorization
  Sparsity getFactorizationSparsity(bool transpose=false) const;

  /// Obtain a numeric Cholesky factorization
  DMatrix getFactorization(bool transpose=false) const;

  /// Static creator function
  #ifdef SWIG
  %callback("%s_cb");
  #endif
  static LinearSolver creator(const Sparsity& sp, int rhs){ return CSparseCholesky(sp, rhs);}
  #ifdef SWIG
  %nocallback;
  #endif

};

} // namespace casadi

#endif //CSPARSE_CHOLESKY_HPP

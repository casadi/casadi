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

#ifndef LINEAR_SOLVER_HPP
#define LINEAR_SOLVER_HPP

#include "fx.hpp"

/** \defgroup LinearSolver_doc 
 * 
 * Solves the linear system X*A = B or X*A^T = B for X
 * with A square and non-singular
 *
 *  If A is structurally singular, an error will be thrown during init.
 *  If A is numerically singular, the prepare step will fail.
 *
 *
 * Note: the transposed form is equivalent to A X^T = B^T
 *       which is the same as A x = b  with x = X^T, b = B^T
 */

namespace CasADi{
  
/// Input arguments of a linear solver [linsolIn]
enum LinsolInput{
  /// The square matrix A: sparse, (n x n). [A]
  LINSOL_A,
  /// The right-hand-side matrix b: dense,  (m x n) [B]
  LINSOL_B,
  LINSOL_NUM_IN};

/// Output arguments of a linear solver [linsolOut]
enum LinsolOutput{
  /// Solution to the linear system of equations [X]
  LINSOL_X,
  LINSOL_NUM_OUT};

  // Forward declaration of internal class
  class LinearSolverInternal;

  /** Abstract base class for the linear solver classes
   *  @copydoc LinearSolver_doc
   \author Joel Andersson
   \date 2010-2013
  */
  class LinearSolver : public FX{
  public:
  
    /// Access functions of the node
    LinearSolverInternal* operator->();

    /// Const access functions of the node
    const LinearSolverInternal* operator->() const;

    /// Factorize the matrix
    void prepare();

    /// Solve the system of equations, internal vector
    void solve(bool transpose=false);

#ifndef SWIG
    /// Solve the factorized system of equations
    void solve(double* x, int nrhs=1, bool transpose=false);

    /// Propagate sparsity through a linear solve
    void spSolve(bvec_t* X, bvec_t* B, bool transpose=false) const;

#endif // SWIG

    /// Create a solve node
    MX solve(const MX& A, const MX& B, bool transpose=false);

    /// Check if prepared
    bool prepared() const;
  
    /// Check if the node is pointing to the right type of object
    virtual bool checkNode() const;
  };

} // namespace CasADi

#endif //LINEAR_SOLVER_HPP


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

#include "function.hpp"

/** \defgroup LinearSolver_doc 
 * 
 * Solves the linear system A*X = B or A^T*X = B for X
 * with A square and non-singular
 *
 *  If A is structurally singular, an error will be thrown during init.
 *  If A is numerically singular, the prepare step will fail.
 *
 */

namespace CasADi{
  
/// Input arguments of a linear solver [linsolIn]
enum LinsolInput{
  /// The square matrix A: sparse, (n x n). [A]
  LINSOL_A,
  /// The right-hand-side matrix b: dense,  (n x m) [B]
  LINSOL_B,
  LINSOL_NUM_IN};

/// Output arguments of a linear solver [linsolOut]
enum LinsolOutput{
  /// Solution to the linear system of equations [X]
  LINSOL_X,
  LINSOL_NUM_OUT};

  // Forward declaration of internal class
  class LinearSolverInternal;

  /** Base class for the linear solver classes
   *  @copydoc LinearSolver_doc
   \author Joel Andersson
   \date 2010-2013
  */
  class LinearSolver : public Function{
  public:
    
    /// \cond INTERNAL
    /// Default (empty) constructor
    LinearSolver();
    /// \endcond
  
    /// Create a linear solver given a sparsity pattern (creates a dummy solver only)
    explicit LinearSolver(const Sparsity& sp, int nrhs=1);

    /// \cond INTERNAL
    /// Access functions of the node
    LinearSolverInternal* operator->();

    /// Const access functions of the node
    const LinearSolverInternal* operator->() const;
    /// \endcond

    /// Factorize the matrix
    void prepare();

    /// Solve the system of equations, internal vector
    void solve(bool transpose=false);

/// \cond INTERNAL
#ifndef SWIG
    /// Solve the factorized system of equations
    void solve(double* x, int nrhs=1, bool transpose=false);

    //@{
    /// Propagate sparsity through a linear solve
    void spSolve(bvec_t* X, const bvec_t* B, bool transpose=false) const;
    void spSolve(DMatrix& X, const DMatrix& B, bool transpose=false) const;
    //@}

#endif // SWIG
/// \endcond

    /// Create a solve node
    MX solve(const MX& A, const MX& B, bool transpose=false);

    /// Check if prepared
    bool prepared() const;
  
    /// Check if the node is pointing to the right type of object
    virtual bool checkNode() const;
  };

} // namespace CasADi

#endif //LINEAR_SOLVER_HPP


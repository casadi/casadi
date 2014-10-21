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


#ifndef CASADI_LINEAR_SOLVER_HPP
#define CASADI_LINEAR_SOLVER_HPP

#include "function.hpp"

/** \defgroup LinearSolver_doc
 *
 * Solves the linear system A*X = B or A^T*X = B for X
 * with A square and non-singular
 *
 *  If A is structurally singular, an error will be thrown during init.
 *  If A is numerically singular, the prepare step will fail.
 *
 *  The usual procedure to use LinearSolver is: \n
 *  -# init()
 *  -# set the first input (A)
 *  -# prepare()
 *  -# set the second input (b)
 *  -# solve()
 *  -# Repeat steps 4 and 5 to work with other b vectors.
 *
 * The method evaluate() combines the prepare() and solve() step and is
 * therefore more expensive if A is invariant.
 *
*/

namespace casadi {

  /// Input arguments of a linear solver [linsolIn]
  enum LinsolInput {
    /// The square matrix A: sparse, (n x n). [A]
    LINSOL_A,
    /// The right-hand-side matrix b: dense,  (n x m) [B]
    LINSOL_B,
    LINSOL_NUM_IN};

  /// Output arguments of a linear solver [linsolOut]
  enum LinsolOutput {
    /// Solution to the linear system of equations [X]
    LINSOL_X,
    LINSOL_NUM_OUT};

  // Forward declaration of internal class
  class LinearSolverInternal;

  /** \brief Base class for the linear solver classes

   *  @copydoc LinearSolver_doc

      \generalsection{LinearSolver}
      \pluginssection{LinearSolver}

   \author Joel Andersson
   \date 2010-2013
  */
  class CASADI_CORE_EXPORT LinearSolver : public Function {
  public:

    /// \cond INTERNAL
    /// Default (empty) constructor
    LinearSolver();
    /// \endcond

    /// Create a linear solver given a sparsity pattern (creates a dummy solver only)
    explicit LinearSolver(const Sparsity& sp, int nrhs=1);

    /** \brief Create a linear solver given a sparsity pattern
    * \param name \pluginargument{LinearSolver}
    */
    LinearSolver(const std::string& name, const Sparsity& sp, int nrhs=1);

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

    ///@{
    /// Propagate sparsity through a linear solve
    void spSolve(bvec_t* X, const bvec_t* B, bool transpose=false) const;
    void spSolve(DMatrix& X, const DMatrix& B, bool transpose=false) const;
    ///@}

#endif // SWIG
    /// \endcond

    /// Create a solve node
    MX solve(const MX& A, const MX& B, bool transpose=false);

    /// Check if prepared
    bool prepared() const;

    /** \brief Solve the system of equations <tt>Lx = b</tt>
        Only when a Cholesky factorization is available
    */
    void solveL(double* x, int nrhs, bool transpose);

    /** \brief Obtain a symbolic Cholesky factorization
        Only for Cholesky solvers
    */
    Sparsity getFactorizationSparsity(bool transpose=false) const;

    /** \brief Obtain a numeric Cholesky factorization
        Only for Cholesky solvers
     */
    DMatrix getFactorization(bool transpose=false) const;

    /// Check if a particular cast is allowed
    static bool testCast(const SharedObjectNode* ptr);

    /// Check if a plugin is available
    static bool hasPlugin(const std::string& name);

    /// Explicitly load a plugin dynamically
    static void loadPlugin(const std::string& name);

    /// Get solver specific documentation
    static std::string doc(const std::string& name);
  };

} // namespace casadi

#endif // CASADI_LINEAR_SOLVER_HPP


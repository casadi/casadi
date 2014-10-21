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


#ifndef CASADI_QCQP_SOLVER_HPP
#define CASADI_QCQP_SOLVER_HPP

#include "function.hpp"

/** \defgroup QcqpSolver_doc

    Solves the following strictly convex problem:

    \verbatim
    min          1/2 x' H x + g' x
    x

    subject to
    1/2 x' Pi x  +  qi' x + ri  <= 0   for i=0..nq-1
    LBA <= A x <= UBA
    LBX <= x   <= UBX

    with :
    H, Pi sparse (n x n) positive definite
    g, qi dense  (n x 1)
    ri scalar

    n: number of decision variables (x)
    nc: number of linear constraints (A)
    nq: number of quadratic constraints

    \endverbatim

    If H, Pi is not positive-definite, the solver should throw an error.

*/

namespace casadi {

  /// Input arguments of a QP problem [qcqpIn]
  enum QcqpSolverInput {
    /// The square matrix H: sparse, (n x n). Only the lower triangular part is actually used.
    /// The matrix is assumed to be symmetrical. [h]
    QCQP_SOLVER_H,
    /// The vector g: dense,  (n x 1) [g]
    QCQP_SOLVER_G,
    /// The horizontal stack of all Pi. Each Pi is sparse (n x n). Only the lower
    /// triangular part is actually used. The matrix is assumed to be symmetrical. [p]
    QCQP_SOLVER_P,
    /// The vertical stack of all qi: dense,  (nq n x 1) [q]
    QCQP_SOLVER_Q,
    /// The vertical stack of all scalars ri (nq x 1)  [r]
    QCQP_SOLVER_R,
    /// The matrix A: sparse, (nc x n) - product with x must be dense. [a]
    QCQP_SOLVER_A,
    /// dense, (nc x 1) [lba]
    QCQP_SOLVER_LBA,
    /// dense, (nc x 1) [uba]
    QCQP_SOLVER_UBA,
    /// dense, (n x 1) [lbx]
    QCQP_SOLVER_LBX,
    /// dense, (n x 1) [ubx]
    QCQP_SOLVER_UBX,
    /// dense, (n x 1) [x0]
    QCQP_SOLVER_X0,
    /// dense [lam_x0]
    QCQP_SOLVER_LAM_X0,
    QCQP_SOLVER_NUM_IN};

  /// Output arguments of an QP Solver [qcqpOut]
  enum QcqpSolverOutput {
    /// The primal solution [x]
    QCQP_SOLVER_X,
    /// The optimal cost [cost]
    QCQP_SOLVER_COST,
    /// The dual solution corresponding to linear bounds [lam_a]
    QCQP_SOLVER_LAM_A,
    /// The dual solution corresponding to simple bounds [lam_x]
    QCQP_SOLVER_LAM_X,
    QCQP_SOLVER_NUM_OUT};


  /// Structure specification of a QP [qcqpStruct]
  enum QCQPStruct {
    /// The square matrix H: sparse, (n x n). Only the lower triangular part is actually used.
    /// The matrix is assumed to be symmetrical. [h]
    QCQP_STRUCT_H,
    /// The horizontal stack of all Pi. Each Pi is sparse (n x n). Only the lower
    /// triangular part is actually used. The matrix is assumed to be symmetrical. [p]
    QCQP_STRUCT_P,
    /// The matrix A: sparse, (nc x n) - product with x must be dense. [a]
    QCQP_STRUCT_A,
    QCQP_STRUCT_NUM};

  // Forward declaration of internal class
  class QcqpSolverInternal;

  /** \brief QcqpSolver


      @copydoc QcqpSolver_doc

      \generalsection{QcqpSolver}
      \pluginssection{QcqpSolver}

      \author Joris Gillis
      \date 2013
  */
  class CASADI_CORE_EXPORT QcqpSolver : public Function {
  public:

    /// Default constructor
    QcqpSolver();

    /** \brief Constructor
     *  \param name \pluginargument{QcqpSolver}
     *  \param st \structargument{QCQP}
     */
    QcqpSolver(const std::string& name, const QCQPStructure& st);

    /// Access functions of the node
    QcqpSolverInternal* operator->();
    const QcqpSolverInternal* operator->() const;

    /// Check if a particular cast is allowed
    static bool testCast(const SharedObjectNode* ptr);

    /// Check if a plugin is available
    static bool hasPlugin(const std::string& name);

    /// Explicitly load a plugin dynamically
    static void loadPlugin(const std::string& name);

    /// Get solver specific documentation
    static std::string doc(const std::string& name);

    /// Set options that make the QP solver more suitable for solving LPs
    void setQPOptions();
  };

} // namespace casadi

#endif // CASADI_QCQP_SOLVER_HPP

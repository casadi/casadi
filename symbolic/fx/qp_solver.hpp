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

#ifndef QP_SOLVER_HPP
#define QP_SOLVER_HPP

#include "fx.hpp"


/** \defgroup QPSolver_doc

  Solves the following strictly convex problem:
  
  \verbatim
  min          1/2 x'.H.x + G'.x 
   x
  
  subject to
              LBA <= A.x <= UBA
              LBX <= x   <= UBX
              
      with H positive definite
              
      nx: number of decision variables (x)
      nc: number of constraints (A)
      
  \endverbatim
  
  If H is not positive-definite, the solver should throw an error.
  
*/
      
namespace CasADi{
  
/// Input arguments of a QP problem [qpIn]
enum QPInput{
  /// The square matrix H: sparse, (nx x nx). Only the lower triangular part is actually used. The matrix is assumed to be symmetrical. [h]
  QP_H,
  /// The vector G: dense,  (nx x 1) [g]
  QP_G,
  /// The matrix A: sparse, (nc x nx) - product with x must be dense. [a]
  QP_A,
  /// dense, (nc x 1) [lba]
  QP_LBA,
  /// dense, (nc x 1) [uba]
  QP_UBA,
  /// dense, (nx x 1) [lbx]
  QP_LBX,
  /// dense, (nx x 1) [ubx]
  QP_UBX,
  /// dense, (nx x 1) [x_init]
  QP_X_INIT,
  /// dense [lambda_init]
  QP_LAMBDA_INIT,
  QP_NUM_IN};

/// Output arguments of an QP Solver [qpOut]
enum QPOutput{
  /// The primal solution [primal]
  QP_PRIMAL,
  /// The optimal cost [cost]
  QP_COST,
  /// The dual solution corresponding to linear bounds [lambda_a]
  QP_LAMBDA_A,
  /// The dual solution corresponding to simple bounds [lambda_x]
  QP_LAMBDA_X,
  QP_NUM_OUT};

// Forward declaration of internal class
class QPSolverInternal;

/** \brief QPSolver


@copydoc QPSolver_doc

  \author Joel Andersson 
  \date 2010
*/
class QPSolver : public FX{
  public:

  /// Default constructor
  QPSolver();
  
  /// Access functions of the node
  QPSolverInternal* operator->();
  const QPSolverInternal* operator->() const;

  /// Check if the node is pointing to the right type of object
  virtual bool checkNode() const;
};

} // namespace CasADi

#endif // QP_SOLVER_HPP


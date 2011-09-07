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

namespace CasADi{
  
  /** Quadratic programming solver in the following form
  min          x'Hx + G'x 
  
  subject to
              LBA <= Ax <= UBA
              LBX <= x  <= UBX
              
      nx: number of decision variables (x)
      nc: number of constraints (A)

*/
enum QPInput{
  /// The square matrix H: sparse, (nx x nx). Only the lower triangular part is actually used. The matrix is assumed to be symmetrical.
  QP_H,
  /// The column vector G: dense,  (nx x 1)
  QP_G,
  /// The matrix A: sparse, (nc x nx) - product with x must be dense.
  QP_A,
  /// dense, (nc x 1)
  QP_LBA,
  /// dense, (nc x 1)
  QP_UBA,
  /// dense, (nx x 1)
  QP_LBX,
  /// dense, (nx x 1)
  QP_UBX,
  /// dense, (nx x 1)
  QP_X_INIT,
  QP_LAMBDA_INIT,QP_NUM_IN};

/// Outputs arguments of an QP Solver
enum QPOutput{
  /// The primal solution
  QP_PRIMAL,
  /// The optimal cost
  QP_COST,
  /// The dual solution corresponding to linear bounds
  QP_DUAL_A,
  /// The dual solution corresponding to simple bounds
  QP_DUAL_X,
  
  QP_NUM_OUT};

// Forward declaration of internal class
class QPSolverInternal;

/** \brief QPSolver

Input arguments of an QP Solver CasADi::QPInput: QP_X_INIT,QP_LBX,QP_UBX,QP_LBG,QP_UBG,QP_LAMBDA_INIT\n
Output arguments of an QP Solver CasADi::QPOutput: QP_X_OPT,QP_COST,QP_LAMBDA_OPT,QP_LAMBDA_LBX,QP_LAMBDA_UBX\n

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


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

#ifndef SDQP_SOLVER_HPP
#define SDQP_SOLVER_HPP

#include "fx.hpp"


/** \defgroup SDQPSolver_doc

  Same as an SDPSolver, but with a quadratic objective 1/2 x' H x
  
*/
      
namespace CasADi{
  
/// Input arguments of a SDQP problem [sdqpIn]
enum SDQPInput{
  /// The matrix H: sparse ( n x n) [h]
  SDQP_SOLVER_H,
  /// The vector c: ( n x 1) [c]
  SDQP_SOLVER_C,
  /// The vertical stack of all matrices F_i: ( nm x m) [f]
  SDQP_SOLVER_F,
  /// The matrix G: ( m x m) [g]
  SDQP_SOLVER_G,
  /// The matrix A: ( nc x n) [a]
  SDQP_SOLVER_A,
  /// Lower bounds on Ax ( nc x 1) [lba]
  SDQP_SOLVER_LBA,
  /// Upper bounds on Ax  ( nc x 1) [uba]
  SDQP_SOLVER_UBA,
  /// Lower bounds on x ( n x 1 ) [lbx]
  SDQP_SOLVER_LBX,
  /// Upper bounds on x ( n x 1 ) [ubx]
  SDQP_SOLVER_UBX,
  SDQP_SOLVER_NUM_IN};

/// Output arguments of an SDQP Solver [sdqpOut]
enum SDQPOutput{
  /// The primal solution (n x 1) - may be used as initial guess [x]
  SDQP_SOLVER_X,
  /// The solution P (m x m) - may be used as initial guess [p]
  SDQP_SOLVER_P,
  /// The dual solution (m x m) - may be used as initial guess [dual]
  SDQP_SOLVER_DUAL,
  /// The primal optimal cost (1 x 1) [cost]
  SDQP_SOLVER_COST,
  /// The dual optimal cost (1 x 1) [dual_cost]
  SDQP_SOLVER_DUAL_COST,
  /// The dual solution corresponding to the linear constraints  (nc x 1) [lam_a]
  SDQP_SOLVER_LAM_A,
  /// The dual solution corresponding to simple bounds  (n x 1) [lam_x]
  SDQP_SOLVER_LAM_X,
  SDQP_SOLVER_NUM_OUT};
  
/// Structure specification of an SDQP [sdqpStruct]
enum SDQPStruct{
  /// The matrix H: sparse ( n x n) [h]
  SDQP_STRUCT_H,
  /// The vertical stack of all matrices F_i: ( nm x m) [f]
  SDQP_STRUCT_F,
  /// The matrix G: ( m x m) [g]
  SDQP_STRUCT_G,
  /// The matrix A: ( nc x n) [a]
  SDQP_STRUCT_A,
  SDQP_STRUCT_NUM};
  
// Forward declaration of internal class
class SDQPSolverInternal;

/** \brief SDQPSolver


@copydoc SDQPSolver_doc

  \author Joel Andersson 
  \date 2010
*/
class SDQPSolver : public FX{
  public:

  /// Default constructor
  SDQPSolver();
  
  /// Access functions of the node
  SDQPSolverInternal* operator->();
  const SDQPSolverInternal* operator->() const;

  /// Check if the node is pointing to the right type of object
  virtual bool checkNode() const;
  
  /// Set options that make the SDQP solver more suitable for solving SOCPs
  void setSOCQPOptions();
  
};

} // namespace CasADi

#endif // SDQP_SOLVER_HPP


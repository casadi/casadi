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

#ifndef SDP_SOLVER_HPP
#define SDP_SOLVER_HPP

#include "fx.hpp"


/** \defgroup SDPSolver_doc

  Solves an SDP problem in standard form.
  See http://sdpa.indsys.chuo-u.ac.jp/sdpa/files/sdpa-c.6.2.0.manual.pdf
  
  Primal:

  \verbatim
  min          b' x 
   x 
  subject to
                P = Sum_m A_i x_i - C 
                P positive semidefinite   
              
      with x ( m x 1)
           b ( m x 1 )
           C, A_i  sparse symmetric (n x n)
           X dense symmetric ( n x n )
      
  \endverbatim
  
  This formulation is chosen as primal, because it does not call for a large decision variable space.
  
  Dual:
  
  \verbatim
  max          trace(C Y)
   Y 
  
  subject to
              trace(A_i Y) = b_i
              Y positive semidefinite
              
      with Y dense symmetric ( n x n)

  \endverbatim
  
  On generality: you might have formulation with block partitioning:
  
  Primal:
  
  \verbatim
  min          b' x 
   x 
  subject to
                Pj = Sum_m A_ij x_i - Cj   for all j
                Pj positive semidefinite   for all j
              
      with x ( m x 1)
           b ( m x 1 )
           C, A_i  sparse symmetric (n x n)
           X dense symmetric ( n x n )
      
  \endverbatim
  
  Dual:
  \verbatim
  max          Sum_j trace(Cj Yj)
   Yj 
  
  subject to
              Sum_j trace(A_ij Yj) = b_i   for all j
              Yj positive semidefinite     for all j
              
      with Y dense symmetric ( n x n)

  \endverbatim
  
  You can cast this into the standard form with:
    C  = blkdiag(Cj for all j)
    Ai = blkdiag(A_ij for all j)
    
  Implementations of SDPSolver are encouraged to exploit this block structure.
  
*/
      
namespace CasADi{
  
/// Input arguments of a SDP problem [sdpIn]
enum SDPInput{
  /// The vertical stack of all matrices A_i: ( nm x n) [a]
  SDP_A,
  /// The vector b: ( m x 1) [b]
  SDP_B,
  /// The matrix C: ( n x n) [c]
  SDP_C,
  SDP_NUM_IN};

/// Output arguments of an SDP Solver [sdpOut]
enum SDPOutput{
  /// The primal solution (m x 1) - may be used as initial guess [primal]
  SDP_PRIMAL,
  /// The solution P (n x n) - may be used as initial guess [p]
  SDP_PRIMAL_P,
  /// The dual solution (n x n) - may be used as initial guess [dual]
  SDP_DUAL,
  /// The primal optimal cost (1 x 1) [primal_cost]
  SDP_PRIMAL_COST,
  /// The dual optimal cost (1 x 1) [dual_cost]
  SDP_DUAL_COST,
  SDP_NUM_OUT};

// Forward declaration of internal class
class SDPSolverInternal;

/** \brief SDPSolver


@copydoc SDPSolver_doc

  \author Joel Andersson 
  \date 2010
*/
class SDPSolver : public FX{
  public:

  /// Default constructor
  SDPSolver();
  
  /// Access functions of the node
  SDPSolverInternal* operator->();
  const SDPSolverInternal* operator->() const;

  /// Check if the node is pointing to the right type of object
  virtual bool checkNode() const;
};

} // namespace CasADi

#endif // SDP_SOLVER_HPP


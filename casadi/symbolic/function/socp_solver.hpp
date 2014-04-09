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

#ifndef SOCP_SOLVER_HPP
#define SOCP_SOLVER_HPP

#include "function.hpp"


/** \defgroup SOCPSolver_doc

  Solves an Second Order Cone Programming (SOCP) problem in standard form.

  Primal:

  \verbatim
  min          c' x
   x
  subject to
                || Gi' x + hi ||_2 <= ei' x + fi  i = 1..m

              LBA <= A x <= UBA
              LBX <= x   <= UBX

      with x ( n x 1)
           c ( n x 1 )
           Gi  sparse (n x ni)
           hi  dense (ni x 1)
           ei  dense (n x 1)
           fi  dense (1 x 1)
           N = Sum_i^m ni
           A sparse (nc x n)
           LBA, UBA dense vector (nc x 1)
           LBX, UBX dense vector (n x 1)

  \endverbatim

*/

namespace casadi{

/// Input arguments of a SOCP problem [socpIn]
enum SOCPInput{
  /// The horizontal stack of all matrices Gi: ( n x N) [g]
  SOCP_SOLVER_G,
  /// The vertical stack of all vectors hi: ( N x 1) [h]
  SOCP_SOLVER_H,
  /// The vertical stack of all vectors ei: ( nm x 1) [e]
  SOCP_SOLVER_E,
  /// The vertical stack of all scalars fi: ( m x 1) [f]
  SOCP_SOLVER_F,
  /// The vector c: ( n x 1) [c]
  SOCP_SOLVER_C,
  /// The matrix A: ( nc x n) [a]
  SOCP_SOLVER_A,
  /// Lower bounds on Ax ( nc x 1) [lba]
  SOCP_SOLVER_LBA,
  /// Upper bounds on Ax  ( nc x 1) [uba]
  SOCP_SOLVER_UBA,
  /// Lower bounds on x ( n x 1 ) [lbx]
  SOCP_SOLVER_LBX,
  /// Upper bounds on x ( n x 1 ) [ubx]
  SOCP_SOLVER_UBX,
  SOCP_SOLVER_NUM_IN};

/// Output arguments of an SOCP Solver [socpOut]
enum SOCPOutput{
  /// The primal solution (n x 1) [x]
  SOCP_SOLVER_X,
  /// The primal optimal cost (1 x 1) [cost]
  SOCP_SOLVER_COST,
  /// The dual solution corresponding to the linear constraints  (nc x 1) [lam_a]
  SOCP_SOLVER_LAM_A,
  /// The dual solution corresponding to simple bounds  (n x 1) [lam_x]
  SOCP_SOLVER_LAM_X,
  SOCP_SOLVER_NUM_OUT};

/// Structure specification of an SOCP [socpStruct]
enum SOCPStruct{
  /// The horizontal stack of all matrices Gi: ( n x N) [g]
  SOCP_STRUCT_G,
  /// The matrix A: ( nc x n) [a]
  SOCP_STRUCT_A,
  SOCP_STRUCT_NUM};

// Forward declaration of internal class
class SOCPSolverInternal;

/** \brief SOCPSolver


@copydoc SOCPSolver_doc

  \author Joris Gillis
  \date 2013
*/
class CASADI_SYMBOLIC_EXPORT SOCPSolver : public Function{
  public:

  /// Default constructor
  SOCPSolver();

  /// Access functions of the node
  SOCPSolverInternal* operator->();
  const SOCPSolverInternal* operator->() const;

  /// Check if the node is pointing to the right type of object
  virtual bool checkNode() const;
};

} // namespace casadi

#endif // SOCP_SOLVER_HPP


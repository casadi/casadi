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

#ifndef LP_SOLVER_HPP
#define LP_SOLVER_HPP

#include "fx.hpp"


/** \defgroup LPSolver_doc

  Solves the following linear problem:
  
  \verbatim
  min          c' x 
   x
  
  subject to
              LBA <= A x <= UBA
              LBX <= x   <= UBX
              
      with x ( n x 1)
           c ( n x 1 )
           A sparse matrix ( nc x n)
           LBA, UBA dense vector (nc x 1)
           LBX, UBX dense vector (n x 1)
              
      n: number of decision variables (x)
      nc: number of constraints (A)
      
  \endverbatim
  

*/
      
namespace CasADi{
  
/// Input arguments of a LP problem [lpIn]
enum LPSolverInput{
  /// The vector c: dense (n x 1) [c]
  LP_SOLVER_C,
  /// The matrix A: sparse, (nc x n) - product with x must be dense. [a]
  LP_SOLVER_A,
  /// dense, (nc x 1) [lba]
  LP_SOLVER_LBA,
  /// dense, (nc x 1) [uba]
  LP_SOLVER_UBA,
  /// dense, (n x 1) [lbx]
  LP_SOLVER_LBX,
  /// dense, (n x 1) [ubx]
  LP_SOLVER_UBX,
  LP_SOLVER_NUM_IN};

/// Output arguments of an LP Solver [lpOut]
enum LPSolverOutput{
  /// The primal solution [x]
  LP_SOLVER_X,
  /// The optimal cost [cost]
  LP_SOLVER_COST,
  /// The dual solution corresponding to linear bounds [lam_a]
  LP_SOLVER_LAM_A,
  /// The dual solution corresponding to simple bounds [lam_x]
  LP_SOLVER_LAM_X,
  LP_SOLVER_NUM_OUT};
  
/// Structure specification of an LP [lpStruct]
enum LPStruct{
  /// The matrix A: sparse [a]
  LP_STRUCT_A,
  LP_STRUCT_NUM};

// Forward declaration of internal class
class LPSolverInternal;

/** \brief LPSolver


@copydoc LPSolver_doc

  \author Joris Gillis 
  \date 2013
*/
class LPSolver : public FX{
  public:

  /// Default constructor
  LPSolver();
  
  /// Access functions of the node
  LPSolverInternal* operator->();
  const LPSolverInternal* operator->() const;

  /// Check if the node is pointing to the right type of object
  virtual bool checkNode() const;
};

} // namespace CasADi

#endif // LP_SOLVER_HPP


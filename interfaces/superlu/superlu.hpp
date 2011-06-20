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

#ifndef SUPERLU_HPP
#define SUPERLU_HPP

#include "casadi/fx/linear_solver.hpp"

namespace CasADi{
  
/** \brief  Forward declaration of internal class */
class SuperLUInternal;

/** \brief  LU LinearSolver with SuperLU Interface
*
* This class solves the linear system A.x=b by making an LU factorization of A: \n
* A = L.U, with L lower and U upper triangular
* 
* SuperLU is an CasADi::FX mapping from 2 inputs [ A (matrix),b (vector)] to one output [x (vector)].
*
* The usual procedure to use SuperLU is: \n
*  -# init()
*  -# set the first input (A)
*  -# prepare()
*  -# set the second input (b)
*  -# solve()
*  -# Repeat steps 4 and 5 to work with other b vectors.
*
* The method evaluate() combines the prepare() and solve() step and is therefore more expensive if A is invariant.
*
*/
class SuperLU : public LinearSolver{
public:

  /// Default (empty) constructor
  SuperLU();
  
  /// Create a linear solver given a sparsity pattern
  SuperLU(const CRSSparsity& sparsity, int nrhs=1);
    
  /** \brief  Access functions of the node */
  SuperLUInternal* operator->();
  const SuperLUInternal* operator->() const;
};

} // namespace CasADi

#endif //SUPERLU_HPP


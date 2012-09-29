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

#ifndef OOQP_SOLVER_HPP
#define OOQP_SOLVER_HPP

#include "symbolic/fx/qp_solver.hpp"

namespace CasADi {
  
  
// Forward declaration of internal class 
class OOQPInternal;

/** OOQP Solver for quadratic programming:

      @copydoc QPSolver_doc

      The current implementation assumes that OOQP is configured with the MA27 sparse linear solver.
      
      NOTE: when doing multiple calls to evaluate(), check if you need to reInit();
*/
class OOQPSolver : public QPSolver {
public:

  /** \brief  Default constructor */
  OOQPSolver();
  
  OOQPSolver(const CRSSparsity& H, const CRSSparsity& A);
  
  /** \brief  Access functions of the node */
  OOQPInternal* operator->();
  const OOQPInternal* operator->() const;

  /// Check if the node is pointing to the right type of object
  virtual bool checkNode() const;
  
  /**
   * \brief Reinitialize the problem 
   * This method needs to be called before evaluate() whenever the nature of any constraint has changed. This occurs when: \n
   *  - Any of LBA, UBA, LBX, UBX changes to/from (+-)infinity  \n
   *  - An entry of LBA becomes equal/unequal to UBA: this indicates that an inequality becomes an equality or visa versa. \n
   * 
   * You do not need to call this method before doing the very first evaluate() run
   */
  void reInit();
  
  /// Static creator function
  #ifdef SWIG
  %callback("%s_cb");
  #endif
  static QPSolver creator(const CRSSparsity& H, const CRSSparsity& A){ return OOQPSolver(H,A);}
  #ifdef SWIG
  %nocallback;
  #endif

  
};


} // namespace CasADi

#endif //OOQP_SOLVER_HPP


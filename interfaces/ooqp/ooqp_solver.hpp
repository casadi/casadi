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

#include "symbolic/function/qp_solver.hpp"

namespace CasADi {
  
  
// Forward declaration of internal class 
class OOQPInternal;

/** \brief Interface to the OOQP Solver for quadratic programming:

      @copydoc QPSolver_doc

      The current implementation assumes that OOQP is configured with the MA27 sparse linear solver.
      
      NOTE: when doing multiple calls to evaluate(), check if you need to reInit();
*/
class OOQPSolver : public QPSolver {
public:

  /** \brief  Default constructor */
  OOQPSolver();
  
  /** \brief Constructor
  *  \param st Problem structure
  *  \copydoc scheme_QPStruct
  */
  OOQPSolver(const QPStructure &st);
  
  /** \brief  Access functions of the node */
  OOQPInternal* operator->();
  const OOQPInternal* operator->() const;

  /// Check if the node is pointing to the right type of object
  virtual bool checkNode() const;
  
  /// Static creator function
  #ifdef SWIG
  %callback("%s_cb");
  #endif
  static QPSolver creator(const QPStructure &st){ return OOQPSolver(st);}
  #ifdef SWIG
  %nocallback;
  #endif

  
};


} // namespace CasADi

#endif //OOQP_SOLVER_HPP


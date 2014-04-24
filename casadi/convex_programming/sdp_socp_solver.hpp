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

#ifndef SDP_SOCP_SOLVER_HPP
#define SDP_SOCP_SOLVER_HPP

#include "casadi/core/function/socp_solver.hpp"

#include <casadi/convex_programming/casadi_convex_programming_export.h>

namespace casadi {


// Forward declaration of internal class
class SDPSOCPInternal;

  /** \brief SOCP Solver for quadratic programming

   @copydoc SOCPSolver_doc

   \author Joris Gillis
   \date 2013
  */
class CASADI_CONVEX_PROGRAMMING_EXPORT SDPSOCPSolver : public SOCPSolver {
public:

  /** \brief  Default constructor */
  SDPSOCPSolver();


  /** \brief Constructor
  *  \param st Problem structure
  *  \copydoc scheme_SOCPStruct
  */
  explicit SDPSOCPSolver(const SOCPStructure & st);

  /** \brief  Access functions of the node */
  SDPSOCPInternal* operator->();
  const SDPSOCPInternal* operator->() const;

  /// Check if the node is pointing to the right type of object
  virtual bool checkNode() const;

  /// Static creator function
  #ifdef SWIG
  %callback("%s_cb");
  #endif
  static SOCPSolver creator(const SOCPStructure & st) { return SDPSOCPSolver(st);}
  #ifdef SWIG
  %nocallback;
  #endif

  /// Access underlying SDP solver
  SDPSolver & getSolver();

};


} // namespace casadi

#endif //SDP_SOCP_SOLVER_HPP


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

#ifndef SDP_SDQP_SOLVER_HPP
#define SDP_SDQP_SOLVER_HPP

#include "symbolic/fx/qcqp_solver.hpp"

namespace CasADi {
  
  
// Forward declaration of internal class 
class SDPSDQPInternal;

  /** \brief SDP SDQP Solver for quadratic programming
   *
   *  Note: this implementation relies on Cholesky decomposition:  Chol(H) = L ->  H = LL' with L lower triangular
   *   This requires Pi, H to be positive definite. Positive semi-definite is not sufficient.
   *    Notably, H==0  will not work.
   *
   *  A better implementation would rely on matrix square root, but we need singular value decomposition to implement that.
   *
   *
   @copydoc SDQPSolver_doc
      
   \author Joris Gillis
   \date 2013
  */
class SDPSDQPSolver : public SDQPSolver {
public:

  /** \brief  Default constructor */
  SDPSDQPSolver();
  
  /** \brief Constructor
  *  \param st Problem structure
  *  \copydoc scheme_SDQPStruct
  */
  explicit SDPSDQPSolver(const SDQPStructure & st);
  
  /** \brief  Access functions of the node */
  SDPSDQPInternal* operator->();
  const SDPSDQPInternal* operator->() const;

  /// Check if the node is pointing to the right type of object
  virtual bool checkNode() const;
  
  /// Static creator function
  #ifdef SWIG
  %callback("%s_cb");
  #endif
  static SDQPSolver creator(const SDQPStructure & st){ return SDPSDQPSolver(st);}
  #ifdef SWIG
  %nocallback;
  #endif
  
  /// Access underlying SDP solver
  SDPSolver & getSolver();

};


} // namespace CasADi

#endif //SDP_SDQP_SOLVER_HPP


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

#ifndef QP_LP_SOLVER_HPP
#define QP_LP_SOLVER_HPP

#include "symbolic/fx/lp_solver.hpp"

namespace CasADi {
  
  
// Forward declaration of internal class 
class QPLPInternal;

  /** \brief IPOPT QP Solver for quadratic programming

   @copydoc LPSolver_doc
      
   \author Joris Gillis
   \date 2013
  */
class QPLPSolver : public LPSolver {
public:

  /** \brief  Default constructor */
  QPLPSolver();
  
  
  /** \brief Constructor
  *  \param st Problem structure
  *  \copydoc scheme_LPStruct
  */
  explicit QPLPSolver(const LPStructure & st);
  
  /** \brief  Access functions of the node */
  QPLPInternal* operator->();
  const QPLPInternal* operator->() const;

  /// Check if the node is pointing to the right type of object
  virtual bool checkNode() const;
  
  /// Static creator function
  #ifdef SWIG
  %callback("%s_cb");
  #endif
  static LPSolver creator(const LPStructure & st){ return QPLPSolver(st);}
  #ifdef SWIG
  %nocallback;
  #endif

};


} // namespace CasADi

#endif //QP_LP_SOLVER_HPP


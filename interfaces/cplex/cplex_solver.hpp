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
#ifndef CPLEX_SOLVER_HPP
#define CPLEX_SOLVER_HPP

#include "symbolic/function/qp_solver.hpp"

namespace CasADi{


// Forward declaration of internal class
class CplexInternal;

  /** \brief Interface to Cplex solver for sparse Quadratic Programs
   @copydoc QPSolver_doc
   \author Attila Kozma, Joel Andersson
   \date 2012
   */
class CplexSolver : public QPSolver{
public:

  /** \brief Default constructor  */
  CplexSolver();

  /** \brief Constructor
  *  \param st Problem structure
  *  \copydoc scheme_QPStruct
  */
  CplexSolver(const QPStructure &st);
  
  CplexInternal* operator->();
  const CplexInternal* operator->() const;

 
  virtual bool checkNode() const;

  /// Static creator function 
  #ifdef SWIG
  %callback("%s_cb");
  #endif
  static QPSolver creator(const QPStructure &st){ return CplexSolver(st);}
  #ifdef SWIG
  %nocallback;
  #endif
};


} // end namespace CasADi

#endif //CPLEX_SOLVER_HPP

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

#ifndef NEWTON_IMPLICIT_SOLVER_HPP
#define NEWTON_IMPLICIT_SOLVER_HPP

#include "casadi/symbolic/function/implicit_function.hpp"

#include <casadi/nonlinear_programming/casadi_nonlinear_programming_export.h>

namespace casadi {


// Forward declaration of internal class
class NewtonImplicitInternal;

  /** \brief Implements simple newton iterations to solve an implicit function.

   @copydoc ImplicitFunction_doc

   \author Joris Gillis
   \date 2012
  */
class CASADI_NONLINEAR_PROGRAMMING_EXPORT NewtonImplicitSolver : public ImplicitFunction {
public:

  /** \brief Default constructor */
  NewtonImplicitSolver();

  /** \brief Create a solver instance */
  explicit NewtonImplicitSolver(const Function& f, const Function& jac=Function(), const LinearSolver& linsol=LinearSolver());

  /** \brief  Access functions of the node */
  NewtonImplicitInternal* operator->();
  const NewtonImplicitInternal* operator->() const;

  /// Check if the node is pointing to the right type of object
  virtual bool checkNode() const;

  /// Static creator function
  #ifdef SWIG
  %callback("%s_cb");
  #endif
  static ImplicitFunction creator(const Function& f, const Function& jac, const LinearSolver& linsol){ return NewtonImplicitSolver(f,jac,linsol);}
  #ifdef SWIG
  %nocallback;
  #endif

};


} // namespace casadi

#endif //NEWTON_IMPLICIT_SOLVER_HPP


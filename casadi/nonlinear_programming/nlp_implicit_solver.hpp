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

#ifndef NLP_IMPLICIT_SOLVER_HPP
#define NLP_IMPLICIT_SOLVER_HPP

#include "casadi/symbolic/function/implicit_function.hpp"

#include <casadi/nonlinear_programming/casadi_nonlinear_programming_export.h>

namespace casadi {


// Forward declaration of internal class
class NLPImplicitInternal;

  /** \brief Use an NLPSolver as ImplicitFunction solver

   @copydoc ImplicitFunction_doc

   \author Joris Gillis
   \date 2012
  */
class CASADI_NONLINEAR_PROGRAMMING_EXPORT NLPImplicitSolver : public ImplicitFunction {
public:

  /** \brief  Default constructor */
  NLPImplicitSolver();

  /** \brief Create a new solver instance */
  explicit NLPImplicitSolver(const Function& f, const Function& jac=Function(),
                             const LinearSolver& linsol=LinearSolver());

  /** \brief  Access functions of the node */
  NLPImplicitInternal* operator->();
  const NLPImplicitInternal* operator->() const;

  /// Check if the node is pointing to the right type of object
  virtual bool checkNode() const;

  /// Access NLP solver
  NLPSolver& getNLPSolver();

  /// Static creator function
  #ifdef SWIG
  %callback("%s_cb");
  #endif
  static ImplicitFunction creator(const Function& f, const Function& jac,
                                  const LinearSolver& linsol)
  { return NLPImplicitSolver(f,jac,linsol);}
  #ifdef SWIG
  %nocallback;
  #endif

};


} // namespace casadi

#endif //NLP_IMPLICIT_SOLVER_HPP


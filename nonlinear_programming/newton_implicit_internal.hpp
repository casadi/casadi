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

#ifndef NEWTON_IMPLICIT_INTERNAL_HPP
#define NEWTON_IMPLICIT_INTERNAL_HPP

#include "newton_implicit_solver.hpp"
#include "symbolic/function/implicit_function_internal.hpp"
#include "symbolic/function/nlp_solver.hpp"
#include "symbolic/function/linear_solver.hpp"

/// \cond INTERNAL
namespace CasADi{

  /** \brief Internal class for NewtonImplicitInternal
   * 
   @copydoc ImplicitFunction_doc
   * */
  class NewtonImplicitInternal : public ImplicitFunctionInternal {
    friend class NewtonImplicitSolver;
  public:
    /** \brief  Constructor */
    explicit NewtonImplicitInternal(const Function& f, const Function& jac, const LinearSolver& linsol);

    /** \brief  Destructor */
    virtual ~NewtonImplicitInternal();

    /** \brief  Clone */
    virtual NewtonImplicitInternal* clone() const{ return new NewtonImplicitInternal(*this);}
  
    /** \brief  Create a new ImplicitFunctionInternal */
    virtual ImplicitFunctionInternal* create(const Function& f, const Function& jac, const LinearSolver& linsol) const { return new NewtonImplicitInternal(f,jac,linsol);}

    /** \brief  Initialize */
    virtual void init();

    /** \brief  Solve the nonlinear system of equations */ 
    virtual void solveNonLinear();
  
  protected:
    /// Maximum number of Newton iterations
    int max_iter_;
    
    /// Absolute tolerance that should be met on residual
    double abstol_;
    
    /// Absolute tolerance that should be met on step
    double abstolStep_;
  };

} // namespace CasADi
/// \endcond
#endif //NEWTON_IMPLICIT_INTERNAL_HPP


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

#ifndef NLP_IMPLICIT_INTERNAL_HPP
#define NLP_IMPLICIT_INTERNAL_HPP

#include "nlp_implicit_solver.hpp"
#include "casadi/symbolic/function/implicit_function_internal.hpp"
#include "casadi/symbolic/function/nlp_solver.hpp"
#include "casadi/symbolic/function/linear_solver.hpp"

/// \cond INTERNAL
namespace casadi{

  /** \brief Internal class for NLPImplicitInternal
   *
   @copydoc ImplicitFunction_doc
   * */
  class CASADI_NONLINEAR_PROGRAMMING_EXPORT NLPImplicitInternal : public ImplicitFunctionInternal {
    friend class NLPImplicitSolver;
  public:
    /** \brief  Constructor */
    explicit NLPImplicitInternal(const Function& f, const Function& jac,
                                 const LinearSolver& linsol);

    /** \brief  Destructor */
    virtual ~NLPImplicitInternal();

    /** \brief  Clone */
    virtual NLPImplicitInternal* clone() const{ return new NLPImplicitInternal(*this);}

    /** \brief  Deep copy data members */
    virtual void deepCopyMembers(std::map<SharedObjectNode*,SharedObject>& already_copied);

    /** \brief  Create a new ImplicitFunctionInternal */
    virtual NLPImplicitInternal* create(const Function& f, const Function& jac,
                                        const LinearSolver& linsol) const
    { return new NLPImplicitInternal(f,jac,linsol);}

    /** \brief  Initialize */
    virtual void init();

    /** \brief  Solve the nonlinear system of equations */
    virtual void solveNonLinear();

    // NLP solver instance
    NLPSolver nlp_solver_;

  };

} // namespace casadi
/// \endcond
#endif //NLP_IMPLICIT_INTERNAL_HPP

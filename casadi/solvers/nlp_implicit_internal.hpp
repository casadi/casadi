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

#ifndef CASADI_NLP_IMPLICIT_INTERNAL_HPP
#define CASADI_NLP_IMPLICIT_INTERNAL_HPP

#include "casadi/core/function/implicit_function_internal.hpp"
#include "casadi/core/function/nlp_solver.hpp"
#include "casadi/core/function/linear_solver.hpp"

#include <casadi/solvers/casadi_implicitfunction_nlp_export.h>

/** \defgroup plugin_ImplicitFunction_nlp
  Use an NlpSolver as ImplicitFunction solver
*/
/** \pluginsection{ImplicitFunction,nlp} */

/// \cond INTERNAL
namespace casadi {

  /** \brief  \pluginbrief{ImplicitFunction,nlp}

   @copydoc ImplicitFunction_doc
   @copydoc plugin_ImplicitFunction_nlp
  
   \author Joris Gillis
   \date 2012
  */
  class CASADI_IMPLICITFUNCTION_NLP_EXPORT NLPImplicitInternal : public ImplicitFunctionInternal {
    friend class NLPImplicitSolver;
  public:
    /** \brief  Constructor */
    explicit NLPImplicitInternal(const Function& f, const Function& jac,
                                 const LinearSolver& linsol);

    /** \brief  Destructor */
    virtual ~NLPImplicitInternal();

    /** \brief  Clone */
    virtual NLPImplicitInternal* clone() const { return new NLPImplicitInternal(*this);}

    /** \brief  Deep copy data members */
    virtual void deepCopyMembers(std::map<SharedObjectNode*, SharedObject>& already_copied);

    /** \brief  Create a new ImplicitFunctionInternal */
    virtual NLPImplicitInternal* create(const Function& f, const Function& jac,
                                        const LinearSolver& linsol) const
    { return new NLPImplicitInternal(f, jac, linsol);}

    /** \brief  Create a new ImplicitFunction */
    static ImplicitFunctionInternal* creator(const Function& f, const Function& jac,
                                             const LinearSolver& linsol)
    { return new NLPImplicitInternal(f, jac, linsol);}

    /** \brief  Initialize */
    virtual void init();

    /** \brief  Solve the nonlinear system of equations */
    virtual void solveNonLinear();

    // NLP solver instance
    NlpSolver nlp_solver_;

    /// A documentation string
    static const std::string meta_doc;

  };

} // namespace casadi
/// \endcond
#endif // CASADI_NLP_IMPLICIT_INTERNAL_HPP

/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A symbolic framework for dynamic optimization.
 *    Copyright (C) 2010-2014 Joel Andersson, Joris Gillis, Moritz Diehl,
 *                            K.U. Leuven. All rights reserved.
 *    Copyright (C) 2011-2014 Greg Horn
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


#ifndef CASADI_NEWTON_HPP
#define CASADI_NEWTON_HPP

#include "casadi/core/function/implicit_function_internal.hpp"
#include "casadi/core/function/nlp_solver.hpp"
#include "casadi/core/function/linear_solver.hpp"

#include <casadi/solvers/casadi_implicitfunction_newton_export.h>

/** \defgroup plugin_ImplicitFunction_newton
     Implements simple newton iterations to solve an implicit function.
*/

/** \pluginsection{ImplicitFunction,newton} */

/// \cond INTERNAL
namespace casadi {

  /** \brief \pluginbrief{ImplicitFunction,newton}

      @copydoc ImplicitFunction_doc
      @copydoc plugin_ImplicitFunction_newton

      \author Joris Gillis
      \date 2012
  */
  class CASADI_IMPLICITFUNCTION_NEWTON_EXPORT Newton
      : public ImplicitFunctionInternal {
  public:
    /** \brief  Constructor */
    explicit Newton(const Function& f);

    /** \brief  Destructor */
    virtual ~Newton();

    /** \brief  Clone */
    virtual Newton* clone() const { return new Newton(*this);}

    /** \brief  Create a new ImplicitFunctionInternal */
    virtual ImplicitFunctionInternal* create(const Function& f) const
    { return new Newton(f);}

    /** \brief  Create a new ImplicitFunction */
    static ImplicitFunctionInternal* creator(const Function& f)
    { return new Newton(f);}

    /** \brief  Initialize */
    virtual void init();

    /** \brief  Solve the nonlinear system of equations */
    virtual void solveNonLinear();

    /// A documentation string
    static const std::string meta_doc;

  protected:
    /// Maximum number of Newton iterations
    int max_iter_;

    /// Absolute tolerance that should be met on residual
    double abstol_;

    /// Absolute tolerance that should be met on step
    double abstolStep_;

    /// If true, each iteration will be printed
    bool print_iteration_;

    /// Print iteration header
    void printIteration(std::ostream &stream);

    /// Print iteration
    void printIteration(std::ostream &stream, int iter, double abstol, double abstolStep);

  };

} // namespace casadi
/// \endcond
#endif // CASADI_NEWTON_HPP

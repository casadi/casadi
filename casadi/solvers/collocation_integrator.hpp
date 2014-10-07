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


#ifndef CASADI_COLLOCATION_INTEGRATOR_HPP
#define CASADI_COLLOCATION_INTEGRATOR_HPP

#include "implicit_fixed_step_integrator.hpp"
#include "casadi/core/function/mx_function.hpp"
#include "casadi/core/function/implicit_function.hpp"
#include "casadi/core/misc/integration_tools.hpp"
#include <casadi/solvers/casadi_integrator_collocation_export.h>

/** \defgroup plugin_Integrator_collocation

     Fixed-step implicit Runge-Kutta integrator
     ODE/DAE integrator based on collocation schemes

     The method is still under development

*/

/** \pluginsection{Integrator,collocation} */

/// \cond INTERNAL
namespace casadi {

  /**
     \brief \pluginbrief{Integrator,collocation}

     @copydoc DAE_doc
     @copydoc plugin_Integrator_collocation

     \author Joel Andersson
     \date 2014
  */
  class CASADI_INTEGRATOR_COLLOCATION_EXPORT CollocationIntegrator :
        public ImplicitFixedStepIntegrator {
  public:

    /// Constructor
    explicit CollocationIntegrator(const Function& f, const Function& g);

    /// Deep copy data members
    virtual void deepCopyMembers(std::map<SharedObjectNode*, SharedObject>& already_copied);

    /// Clone
    virtual CollocationIntegrator* clone() const
    { return new CollocationIntegrator(*this);}

    /// Create a new integrator
    virtual CollocationIntegrator* create(const Function& f, const Function& g) const
    { return new CollocationIntegrator(f, g);}

    /** \brief  Create a new integrator */
    static IntegratorInternal* creator(const Function& f, const Function& g)
    { return new CollocationIntegrator(f, g);}

    /// Destructor
    virtual ~CollocationIntegrator();

    /// Initialize stage
    virtual void init();

    /// Setup F and G
    virtual void setupFG();

    // Return zero if smaller than machine epsilon
    static double zeroIfSmall(double x);

    /// Get initial guess for the algebraic variable
    virtual void calculateInitialConditions();

    /// Get initial guess for the algebraic variable (backward problem)
    virtual void calculateInitialConditionsB();

    // Interpolation order
    int deg_;

    /// A documentation string
    static const std::string meta_doc;

  };

} // namespace casadi
/// \endcond
#endif // CASADI_COLLOCATION_INTEGRATOR_HPP

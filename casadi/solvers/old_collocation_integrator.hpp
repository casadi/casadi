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


#ifndef CASADI_OLD_COLLOCATION_INTEGRATOR_HPP
#define CASADI_OLD_COLLOCATION_INTEGRATOR_HPP

#include "casadi/core/function/integrator_internal.hpp"
#include "casadi/core/function/mx_function.hpp"
#include "casadi/core/function/implicit_function.hpp"
#include "casadi/core/misc/integration_tools.hpp"
#include <casadi/solvers/casadi_integrator_oldcollocation_export.h>

/** \defgroup plugin_Integrator_oldcollocation

     Collocation integrator
     ODE/DAE integrator based on collocation

     The method is still under development

*/
/** \pluginsection{Integrator,oldcollocation} */

/// \cond INTERNAL
namespace casadi {

  /**
     \brief \pluginbrief{Integrator,oldcollocation}

     @copydoc DAE_doc
     @copydoc plugin_Integrator_oldcollocation

     \author Joel Andersson
     \date 2011
  */
  class CASADI_INTEGRATOR_OLDCOLLOCATION_EXPORT
  OldCollocationIntegrator : public IntegratorInternal {

  public:

    /// Constructor
    explicit OldCollocationIntegrator(const Function& f, const Function& g);

    /// Deep copy data members
    virtual void deepCopyMembers(std::map<SharedObjectNode*, SharedObject>& already_copied);

    /// Clone
    virtual OldCollocationIntegrator* clone() const
    { return new OldCollocationIntegrator(*this);}

    /// Create a new integrator
    virtual OldCollocationIntegrator* create(const Function& f, const Function& g) const
    { return new OldCollocationIntegrator(f, g);}

    /** \brief  Create a new integrator */
    static IntegratorInternal* creator(const Function& f, const Function& g)
    { return new OldCollocationIntegrator(f, g);}

    /// Destructor
    virtual ~OldCollocationIntegrator();

    /// Initialize stage
    virtual void init();

    /// Reset the forward problem and bring the time back to t0
    virtual void reset();

    /// Reset the backward problem and take time to tf
    virtual void resetB() {}

    ///  Integrate until a specified time point
    virtual void integrate(double t_out);

    /// Integrate backwards in time until a specified time point
    virtual void integrateB(double t_out);

    // Startup integrator (generates an initial trajectory guess)
    Integrator startup_integrator_;

    // Implicit function solver
    ImplicitFunction implicit_solver_;

    // With hotstart
    bool hotstart_;

    // Has the system been integrated once
    bool integrated_once_;

    // Collocated times
    std::vector<std::vector<double> > coll_time_;

    /// A documentation string
    static const std::string meta_doc;

  };

} // namespace casadi
/// \endcond
#endif // CASADI_OLD_COLLOCATION_INTEGRATOR_HPP

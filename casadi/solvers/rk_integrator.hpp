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


#ifndef CASADI_RK_INTEGRATOR_HPP
#define CASADI_RK_INTEGRATOR_HPP

#include "fixed_step_integrator.hpp"
#include <casadi/solvers/casadi_integrator_rk_export.h>

/** \defgroup plugin_Integrator_rk
      Fixed-step explicit Runge-Kutta integrator for ODEs
      Currently implements RK4.

      The method is still under development
*/
/** \pluginsection{Integrator,rk} */

/// \cond INTERNAL
namespace casadi {

  /** \brief \pluginbrief{Integrator,rk}


      @copydoc DAE_doc
      @copydoc plugin_Integrator_rk

      \author Joel Andersson
      \date 2011-2014
  */
  class CASADI_INTEGRATOR_RK_EXPORT RkIntegrator : public FixedStepIntegrator {
  public:

    /// Constructor
    explicit RkIntegrator(const Function& f, const Function& g);

    /// Deep copy data members
    virtual void deepCopyMembers(std::map<SharedObjectNode*, SharedObject>& already_copied);

    /// Clone
    virtual RkIntegrator* clone() const { return new RkIntegrator(*this);}

    /// Create a new integrator
    virtual RkIntegrator* create(const Function& f, const Function& g) const
    { return new RkIntegrator(f, g);}

    /** \brief  Create a new integrator */
    static IntegratorInternal* creator(const Function& f, const Function& g)
    { return new RkIntegrator(f, g);}

    /// Destructor
    virtual ~RkIntegrator();

    /// Initialize stage
    virtual void init();

    /// Setup F and G
    virtual void setupFG();

    /// A documentation string
    static const std::string meta_doc;

  };

} // namespace casadi

/// \endcond
#endif // CASADI_RK_INTEGRATOR_HPP

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

#ifndef RK_INTEGRATOR_HPP
#define RK_INTEGRATOR_HPP

#include "fixed_step_integrator.hpp"

namespace casadi{

  class RKIntegratorInternal;

  /** \brief Fixed-step explicit Runge-Kutta integrator for ODEs
      Currently implements RK4.

      The method is still under development

      \author Joel Andersson
      \date 2011-2014
  */
  class CASADI_INTEGRATION_EXPORT RKIntegrator : public FixedStepIntegrator {
  public:
    /** \brief  Default constructor */
    RKIntegrator();

    /** \brief  Create an integrator for explicit ODEs
     *   \param f dynamical system
     * \copydoc scheme_DAEInput
     * \copydoc scheme_DAEOutput
     *   \param g backwards system
     * \copydoc scheme_RDAEInput
     * \copydoc scheme_RDAEOutput
     */
    explicit RKIntegrator(const Function& f, const Function& g=Function());

    //@{
    /// Access functions of the node
    RKIntegratorInternal* operator->();
    const RKIntegratorInternal* operator->() const;
    //@}

    /// Check if the node is pointing to the right type of object
    virtual bool checkNode() const;

    /// Static creator function
#ifdef SWIG
    %callback("%s_cb");
#endif
    static Integrator creator(const Function& f, const Function& g){ return RKIntegrator(f,g);}
#ifdef SWIG
    %nocallback;
#endif

  };

} // namespace casadi

#endif //RK_INTEGRATOR_HPP

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

#ifndef OLD_COLLOCATION_INTEGRATOR_HPP
#define OLD_COLLOCATION_INTEGRATOR_HPP

#include "casadi/symbolic/function/integrator.hpp"
#include <casadi/integration/casadi_integration_export.h>

namespace casadi{

  class OldCollocationIntegratorInternal;

  /**
     \brief Collocation integrator
     ODE/DAE integrator based on collocation

     The method is still under development

     @copydoc DAE_doc

     \author Joel Andersson
     \date 2011
  */
  class CASADI_INTEGRATION_EXPORT OldCollocationIntegrator : public Integrator {
  public:
    /** \brief  Default constructor */
    OldCollocationIntegrator();

    /** \brief  Create an integrator for explicit ODEs
     *   \param f dynamical system
     * \copydoc scheme_DAEInput
     * \copydoc scheme_DAEOutput
     *   \param g backwards system
     * \copydoc scheme_RDAEInput
     * \copydoc scheme_RDAEOutput
     */
    explicit OldCollocationIntegrator(const Function& f, const Function& g=Function());

    //@{
    /// Access functions of the node
    OldCollocationIntegratorInternal* operator->();
    const OldCollocationIntegratorInternal* operator->() const;
    //@}

    /// Check if the node is pointing to the right type of object
    virtual bool checkNode() const;

    /// Static creator function
#ifdef SWIG
    %callback("%s_cb");
#endif
    static Integrator creator(const Function& f, const Function& g){ return OldCollocationIntegrator(f,g);}
#ifdef SWIG
    %nocallback;
#endif
  };

} // namespace casadi

#endif //OLD_COLLOCATION_INTEGRATOR_HPP

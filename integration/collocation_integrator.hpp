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

#ifndef COLLOCATION_INTEGRATOR_HPP
#define COLLOCATION_INTEGRATOR_HPP

#include "symbolic/fx/integrator.hpp"

namespace CasADi{
  
  class CollocationIntegratorInternal;
  
  /**
     \brief Collocation integrator
     ODE/DAE integrator based on collocation
  
     The method is still under development
  
     @copydoc DAE_doc
  
     \author Joel Andersson
     \date 2011
  */
  class CollocationIntegrator : public Integrator {
  public:
    /** \brief  Default constructor */
    CollocationIntegrator();
    
    /** \brief  Create an integrator for explicit ODEs
     *   \param f dynamical system
     * \copydoc scheme_DAEInput
     * \copydoc scheme_DAEOutput
     *   \param g backwards system
     * \copydoc scheme_RDAEInput
     * \copydoc scheme_RDAEOutput
     */
    explicit CollocationIntegrator(const FX& f, const FX& g=FX());

    //@{
    /// Access functions of the node
    CollocationIntegratorInternal* operator->();
    const CollocationIntegratorInternal* operator->() const;
    //@}

    /// Check if the node is pointing to the right type of object
    virtual bool checkNode() const;

    /// Static creator function
#ifdef SWIG
    %callback("%s_cb");
#endif
    static Integrator creator(const FX& f, const FX& g){ return CollocationIntegrator(f,g);}
#ifdef SWIG
    %nocallback;
#endif
  };

} // namespace CasADi

#endif //COLLOCATION_INTEGRATOR_HPP

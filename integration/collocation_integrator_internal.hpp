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

#ifndef COLLOCATION_INTEGRATOR_INTERNAL_HPP
#define COLLOCATION_INTEGRATOR_INTERNAL_HPP

#include "collocation_integrator.hpp"
#include "symbolic/fx/integrator_internal.hpp"
#include "symbolic/fx/mx_function.hpp"
#include "symbolic/fx/implicit_function.hpp"
#include "integration_tools.hpp"

namespace CasADi{
    
  class CollocationIntegratorInternal : public IntegratorInternal{
  
  public:
  
    /// Constructor
    explicit CollocationIntegratorInternal(const FX& f, const FX& g);

    /// Deep copy data members
    virtual void deepCopyMembers(std::map<SharedObjectNode*,SharedObject>& already_copied);

    /// Clone
    virtual CollocationIntegratorInternal* clone() const{ return new CollocationIntegratorInternal(*this);}

    /// Create a new integrator
    virtual CollocationIntegratorInternal* create(const FX& f, const FX& g) const{ return new CollocationIntegratorInternal(f,g);}
  
    /// Destructor
    virtual ~CollocationIntegratorInternal();

    /// Initialize stage
    virtual void init();
  
    /// Reset the forward problem and bring the time back to t0
    virtual void reset();

    /// Reset the backward problem and take time to tf
    virtual void resetB(){}

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
  
  };

} // namespace CasADi

#endif //COLLOCATION_INTEGRATOR_INTERNAL_HPP

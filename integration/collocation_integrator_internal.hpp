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
#include "implicit_fixed_step_integrator_internal.hpp"
#include "symbolic/function/mx_function.hpp"
#include "symbolic/function/implicit_function.hpp"
#include "integration_tools.hpp"

/// \cond INTERNAL
namespace CasADi{
    
  class CollocationIntegratorInternal : public ImplicitFixedStepIntegratorInternal{
  public:
  
    /// Constructor
    explicit CollocationIntegratorInternal(const Function& f, const Function& g);

    /// Deep copy data members
    virtual void deepCopyMembers(std::map<SharedObjectNode*,SharedObject>& already_copied);

    /// Clone
    virtual CollocationIntegratorInternal* clone() const{ return new CollocationIntegratorInternal(*this);}

    /// Create a new integrator
    virtual CollocationIntegratorInternal* create(const Function& f, const Function& g) const{ return new CollocationIntegratorInternal(f,g);}
  
    /// Destructor
    virtual ~CollocationIntegratorInternal();

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
  };

} // namespace CasADi
/// \endcond
#endif //COLLOCATION_INTEGRATOR_INTERNAL_HPP

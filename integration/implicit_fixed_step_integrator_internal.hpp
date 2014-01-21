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

#ifndef IMPLICIT_FIXED_STEP_INTEGRATOR_INTERNAL_HPP
#define IMPLICIT_FIXED_STEP_INTEGRATOR_INTERNAL_HPP

#include "implicit_fixed_step_integrator.hpp"
#include "fixed_step_integrator_internal.hpp"
#include "symbolic/fx/implicit_function.hpp"

namespace CasADi{
    
  class ImplicitFixedStepIntegratorInternal : public FixedStepIntegratorInternal{
  public:
  
    /// Constructor
    explicit ImplicitFixedStepIntegratorInternal(const FX& f, const FX& g);

    /// Deep copy data members
    virtual void deepCopyMembers(std::map<SharedObjectNode*,SharedObject>& already_copied);

    /// Clone
    virtual ImplicitFixedStepIntegratorInternal* clone() const = 0;

    /// Create a new integrator
    virtual ImplicitFixedStepIntegratorInternal* create(const FX& f, const FX& g) const = 0;
  
    /// Destructor
    virtual ~ImplicitFixedStepIntegratorInternal();

    /// Initialize stage
    virtual void init();

    /// Get explicit dynamics
    virtual FX& getExplicit(){ return implicit_solver_;}

    /// Get explicit dynamics (backward problem)
    virtual FX& getExplicitB(){ return backward_implicit_solver_;}

    // Implicit function solver
    ImplicitFunction implicit_solver_, backward_implicit_solver_;
  };

} // namespace CasADi

#endif //IMPLICIT_FIXED_STEP_INTEGRATOR_INTERNAL_HPP

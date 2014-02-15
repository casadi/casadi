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

#ifndef FIXED_STEP_INTEGRATOR_HPP
#define FIXED_STEP_INTEGRATOR_HPP

#include "symbolic/fx/integrator.hpp"

namespace CasADi{
  
  class FixedStepIntegratorInternal;
  
  /**
     \brief Base class for fixed step integrators
  
     \author Joel Andersson
     \date 2014
  */
  class FixedStepIntegrator : public Integrator {
  public:
    /** \brief  Default constructor */
    FixedStepIntegrator();
    
    /// Access functions of the node
    FixedStepIntegratorInternal* operator->();

    /// Access functions of the node
    const FixedStepIntegratorInternal* operator->() const;

    /// Check if the node is pointing to the right type of object
    virtual bool checkNode() const;
    
  };

} // namespace CasADi

#endif //FIXED_STEP_INTEGRATOR_HPP

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

#ifndef IMPLICIT_FIXED_STEP_INTEGRATOR_HPP
#define IMPLICIT_FIXED_STEP_INTEGRATOR_HPP

#include "fixed_step_integrator.hpp"

namespace casadi {

  class ImplicitFixedStepIntegratorInternal;

  /**
     \brief Base class for implicit fixed step integrators

     \author Joel Andersson
     \date 2014
  */
  class CASADI_INTEGRATION_EXPORT ImplicitFixedStepIntegrator : public FixedStepIntegrator {
  public:
    /** \brief  Default constructor */
    ImplicitFixedStepIntegrator();

    /// Access functions of the node
    ImplicitFixedStepIntegratorInternal* operator->();

    /// Access functions of the node
    const ImplicitFixedStepIntegratorInternal* operator->() const;

    /// Check if the node is pointing to the right type of object
    virtual bool checkNode() const;

  };

} // namespace casadi

#endif //IMPLICIT_FIXED_STEP_INTEGRATOR_HPP

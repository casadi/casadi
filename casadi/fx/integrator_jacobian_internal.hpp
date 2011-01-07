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

#ifndef INTEGRATOR_JACOBIAN_INTERNAL_HPP
#define INTEGRATOR_JACOBIAN_INTERNAL_HPP

#include "integrator_jacobian.hpp"
#include "integrator.hpp"
#include "fx_internal.hpp"

namespace CasADi{

class IntegratorJacobianInternal : public FXInternal{
public:
  /** \brief  Constructor */
  explicit IntegratorJacobianInternal(const Integrator& integrator);

  /** \brief  Destructor */
  virtual ~IntegratorJacobianInternal();

  /** \brief  Clone */
  virtual IntegratorJacobianInternal* clone() const;
    
  /** \brief  evaluate */
  virtual void evaluate(int fsens_order, int asens_order);

  /** \brief  Initialize */
  virtual void init();

  protected:
    // Integrator integrating the ODE/DAE augmented with forward sensitivity equations
    Integrator integrator_;

    // Number of sensitivities
    int ns_;
  
    // Number of states
    int nx_;
    
    // Mapping between augmented dae states and jacobian
    std::vector<int> jacmap_;
    
    // Intitial value for the augmented dae states
    std::vector<double> jacinit_;
    
};
  
} // namespace CasADi

#endif // INTEGRATOR_JACOBIAN_INTERNAL_HPP

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

#ifndef INTEGRATOR_JACOBIAN_HPP
#define INTEGRATOR_JACOBIAN_HPP

#include "fx.hpp"

namespace CasADi{

/// Output arguments of an integratorJacobian
enum IntegratorJacobianOutput{ 
  /** Jacobian matrix */
  INTEGRATORJACOBIAN_J,   
  /**  Differential state at tf */
  INTEGRATORJACOBIAN_XF, 
  /**  Differential state derivative at tf */ 
  INTEGRATORJACOBIAN_XPF, 
  /** Algebraic state at tf*/
  INTEGRATORJACOBIAN_ZF,  
  INTEGRATORJACOBIAN_NUM_OUT
};

// Forward declaration of internal class
class IntegratorJacobianInternal;

// Forward declaration of Integrator
class Integrator;

/** Evaluate the Jacobian of an integrator by means of integrating the forward sensitivity equations.

  This class maps an integrator which is able to calculate sensitivities as a Jacobian
  This will enable a uniform treatment of all kinds of function, whether SXFunction, MXFunction, etc.
  
  IntegratorJacobian is an CasADi::FX mapping from CasADi::IntegratorInput to CasADi::IntegratorJacobianOutput.  
  
  \author Joel Andersson
  \date 2010
*/

class IntegratorJacobian : public FX{
public:
  /// Default constructor
  IntegratorJacobian();

  /// Constructor
  explicit IntegratorJacobian(const Integrator& integrator);
  
  /// Access functions of the node
  IntegratorJacobianInternal* operator->();

  /// Const access functions of the node
  const IntegratorJacobianInternal* operator->() const;
  
  /// Check if the node is pointing to the right type of object
  virtual bool checkNode() const;
};

} // namespace CasADi

#endif //INTEGRATOR_JACOBIAN_HPP

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

#ifndef RK_INTEGRATOR_INTERNAL_HPP
#define RK_INTEGRATOR_INTERNAL_HPP

#include "rk_integrator.hpp"
#include "symbolic/fx/integrator_internal.hpp"
#include "symbolic/fx/mx_function.hpp"
#include "symbolic/fx/implicit_function.hpp"

namespace CasADi{
    
class RKIntegratorInternal : public IntegratorInternal{

public:
  
  /// Constructor
  explicit RKIntegratorInternal(const FX& f, const FX& g);

  /// Deep copy data members
  virtual void deepCopyMembers(std::map<SharedObjectNode*,SharedObject>& already_copied);

  /// Clone
  virtual RKIntegratorInternal* clone() const{ return new RKIntegratorInternal(*this);}

  /// Create a new integrator
  virtual RKIntegratorInternal* create(const FX& f, const FX& g) const{ return new RKIntegratorInternal(f,g);}
  
  /// Destructor
  virtual ~RKIntegratorInternal();

  /// Initialize stage
  virtual void init();
  
  /// Initialize the adjoint problem (can only be called after the first integration)
  virtual void initAdj();

  /// Reset the forward problem and bring the time back to t0
  virtual void reset();

  /// Reset the backward problem and take time to tf
  virtual void resetB();

  ///  Integrate until a specified time point
  virtual void integrate(double t_out);

  /// Integrate backward in time until a specified time point
  virtual void integrateB(double t_out);

  /// Generate a function that calculates a Jacobian function
  virtual FX getJacobian(int iind, int oind, bool compact, bool symmetric);

  /// Generate the sparsity of a Jacobian block
  virtual CRSSparsity getJacSparsity(int iind, int oind);

  // Function which returns the state at the final time
  FX yf_fun_;  
};

} // namespace CasADi

#endif //RK_INTEGRATOR_INTERNAL_HPP

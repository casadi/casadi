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
#include "casadi/fx/integrator_internal.hpp"
#include "casadi/fx/mx_function.hpp"
#include "casadi/fx/implicit_function.hpp"

namespace CasADi{
    
class CollocationIntegratorInternal : public IntegratorInternal{

public:
  
  /// Constructor
  explicit CollocationIntegratorInternal(const FX& f, const FX& q);

  /// Deep copy data members
  virtual void deepCopyMembers(std::map<SharedObjectNode*,SharedObject>& already_copied);

  /// Clone
  virtual CollocationIntegratorInternal* clone() const{ return new CollocationIntegratorInternal(*this);}

  /// Create a new integrator
  virtual CollocationIntegratorInternal* create(const FX& f, const FX& q) const{ return new CollocationIntegratorInternal(f,q);}
  
  /// Destructor
  virtual ~CollocationIntegratorInternal();

  /// Initialize stage
  virtual void init();
  
  /// Initialize the adjoint problem (can only be called after the first integration)
  virtual void initAdj();

  /// Reset the solver and bring the time back to t0
  virtual void reset(int nfdir, int nadir);

  /// Reset the solver of the adjoint problem and take time to tf
  virtual void resetAdj();

  ///  Integrate until a specified time point
  virtual void integrate(double t_out);

  /// Integrate backwards in time until a specified time point
  virtual void integrateAdj(double t_out);

  /// Get the jacobian in the nonlinear iteration
  virtual FX getJacobian();
  
  /// Get the Linear solver
  virtual LinearSolver getLinearSolver();

  /// Set linear solver
  virtual void setLinearSolver(const LinearSolver& linsol, const FX& jac);

  /// Print statistics
  virtual void printStats(std::ostream &stream) const;

  /// Set the stop time of the forward integration
  virtual void setStopTime(double tf);

  // Nonlinear system of equations
  FX gfcn_;

  // Implicit function solver
  ImplicitFunction implicit_solver_;
  
  // With hotstart
  bool hotstart_;
  
  // Has the system been integrated once
  bool integrated_once_;
  
  // Number of sensitivity directions
  int nfdir_, nadir_;
};

} // namespace CasADi

#endif //COLLOCATION_INTEGRATOR_INTERNAL_HPP

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

#ifndef INTEGRATOR_INTERNAL_HPP
#define INTEGRATOR_INTERNAL_HPP

#include "integrator.hpp"
#include "fx_internal.hpp"

namespace CasADi{

/** \brief Internal storage for integrator related data
  \author Joel Andersson 
  \date 2010
*/
class IntegratorInternal : public FXInternal{
public:
  /** \brief  Constructor */
  IntegratorInternal();

  /** \brief  Destructor */
  virtual ~IntegratorInternal()=0;

  /** \brief  Clone */
  virtual IntegratorInternal* clone() const=0;
  
  /** \brief  Set linear solver */
  virtual void setLinearSolver(const LinearSolver& linsol, const FX& jac)=0;
  
  /** \brief  Print solver statistics */
  virtual void printStats(std::ostream &stream) const = 0;

    /** \brief  Reset the solver and bring the time back to t0 */
  virtual void reset(int fsens_order, int asens_order) = 0;

    /** \brief  Reset the solver of the adjoint problem and take time to tf */
  virtual void resetAdj() = 0;

  /** \brief  Integrate until a specified time point */
  virtual void integrate(double t_out) = 0;

  /** \brief  Integrate backwards in time until a specified time point */
  virtual void integrateAdj(double t_out) = 0;

  /** \brief  Set stop time for the integration */
  virtual void setStopTime(double tf) = 0;
  
  /** \brief  evaluate */
  virtual void evaluate(int fsens_order, int asens_order);

  /** \brief  Initialize */
  virtual void init();

  /** \brief Create an integrator which integrates the ODE/DAE augmented with the forward sensitivity equations */
  virtual Integrator jac(int iind=0, int oind=0) = 0;

  /** \brief Jacobian of output oind with respect to input iind */
  virtual FX jacobian(int iind=0, int oind=0);
 
  /// Get the Jacobian
  virtual FX getJacobian() = 0;
  
  /// Get the Linear solver
  virtual LinearSolver getLinearSolver() = 0;

  /// Number of differential states (including quadrature states)
  int nx_;
  
  /// Number of parameters
  int np_;
  
  /// Number of algebraic states
  int nz_;
  
  /// Number of right hand sides
  int nrhs_;
  
  /// Current time
  double t_;

  // Do not integrate past the end point
  bool stop_at_end_;
  
  //@{
  /// options
  bool exact_jacobian_;
  double abstol_, reltol_;
  double fsens_abstol_, fsens_reltol_;
  double asens_abstol_, asens_reltol_;
  int max_num_steps_;
  bool finite_difference_fsens_;  
  //@}
  
  protected:
    // Set dimensions
    void setDimensions(int nx, int np, int nz);
};
  
} // namespace CasADi

#endif // INTEGRATOR_INTERNAL_HPP

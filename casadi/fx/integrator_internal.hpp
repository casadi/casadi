/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A minimalistic computer algebra system with automatic differentiation 
 *    and framework for dynamic optimization.
 *    Copyright (C) 2010 by Joel Andersson, Moritz Diehl et al., K.U.Leuven. All rights reserved.
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

namespace CasADi{

/** \brief Internal storage for integrator related data
  \author Joel Andersson 
  \date 2010
*/
class IntegratorInternal : public FXNode{
public:
  /** \brief  Constructor */
  IntegratorInternal(int nx, int np);

  /** \brief  Destructor */
  virtual ~IntegratorInternal()=0;
  
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

  /** \brief  Set linear solver */
  void setLinearSolver(LinearSolver::Creator creator);

  /** Lenght of the state vector 
   (also includes "states" evaluated from quadrature formulas)
  */
  int nx_;
  
  /// Number of parameters
  int np_;
  
  /// Current time
  double t_;

  //@{
  /// options
  bool exact_jacobian_;
  double abstol_, reltol_;
  double fsens_abstol_, fsens_reltol_;
  double asens_abstol_, asens_reltol_;
  int max_num_steps_;
  bool finite_difference_fsens_;  
  //@}
  
  /// Linear solver creator
  LinearSolver::Creator linsolve_creator_;
  
  /// Create linear solver
  LinearSolver createLinearSolver(int nrow, int ncol, const std::vector<int>& rowind, const std::vector<int>& col, int nrhs=1) const;
  
};
  
} // namespace CasADi

#endif // INTEGRATOR_INTERNAL_HPP

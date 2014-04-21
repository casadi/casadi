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

#ifndef SUNDIALS_INTERNAL_HPP
#define SUNDIALS_INTERNAL_HPP

#include "sundials_integrator.hpp"
#include "casadi/symbolic/function/integrator_internal.hpp"

#include <nvector/nvector_serial.h>
#include <sundials/sundials_dense.h>
#include <sundials/sundials_iterative.h>
#include <sundials/sundials_types.h>

/// \cond INTERNAL
namespace casadi{

class CASADI_SUNDIALS_INTERFACE_EXPORT SundialsInternal : public IntegratorInternal{
public:
  /** \brief  Constructor */
  SundialsInternal(const Function& f, const Function& g);

  /** \brief  Destructor */
  virtual ~SundialsInternal()=0;

  /** \brief  Initialize */
  virtual void init();

  /** \brief  Reset the forward problem and bring the time back to t0 */
  virtual void reset() = 0;

  /** \brief  Deep copy data members */
  virtual void deepCopyMembers(std::map<SharedObjectNode*,SharedObject>& already_copied);

  /** \brief  Set stop time for the integration */
  virtual void setStopTime(double tf) = 0;

  /// Linear solver forward, backward
  LinearSolver linsol_, linsolB_;

  ///@{
  /// options
  bool exact_jacobian_, exact_jacobianB_;
  double abstol_, reltol_;
  double fsens_abstol_, fsens_reltol_;
  double abstolB_, reltolB_;
  int max_num_steps_;
  bool finite_difference_fsens_;
  bool stop_at_end_;
  ///@}

  /// number of checkpoints stored so far
  int ncheck_;

  /// Supported linear solvers in Sundials
  enum LinearSolverType{SD_USER_DEFINED, SD_DENSE, SD_BANDED, SD_ITERATIVE};

  /// Supported iterative solvers in Sundials
  enum IterativeSolverType{SD_GMRES,SD_BCGSTAB,SD_TFQMR};

  /// Linear solver data (dense)
  struct LinSolDataDense{};

  /// Linear solver
  LinearSolverType linsol_f_, linsol_g_;

  /// Iterative solver
  IterativeSolverType itsol_f_, itsol_g_;

  /// Preconditioning
  int pretype_f_, pretype_g_;

  /// Max Krylov size
  int max_krylov_, max_krylovB_;

  /// Use preconditioning
  bool use_preconditioner_, use_preconditionerB_;

  // Jacobian of the DAE with respect to the state and state derivatives
  Function jac_, jacB_;

  // Jacobian times vector functions
  Function f_fwd_, g_fwd_;

  /** \brief  Get the integrator Jacobian for the forward problem */
  virtual Function getJac()=0;

  /** \brief  Get the integrator Jacobian for the backward problem */
  virtual Function getJacB()=0;

};

} // namespace casadi

/// \endcond
#endif // SUNDIALS_INTERNAL_HPP

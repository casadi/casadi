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
#include "casadi/fx/integrator_internal.hpp"

#include <nvector/nvector_serial.h>
#include <sundials/sundials_dense.h>
#include <sundials/sundials_iterative.h>
#include <sundials/sundials_types.h>

#ifdef WITH_SUNDIALS_2_5
typedef long SUNDIALS_INT;
#else  // WITH_SUNDIALS_2_5
typedef int SUNDIALS_INT;
#endif // WITH_SUNDIALS_2_5

namespace CasADi{
namespace Sundials{

class SundialsInternal : public IntegratorInternal{
public:
  /** \brief  Constructor */
  SundialsInternal(const FX& f, const FX& g);

  /** \brief  Destructor */
  virtual ~SundialsInternal()=0;
    
  /** \brief  Initialize */
  virtual void init();

  /** \brief  Deep copy data members */
  virtual void deepCopyMembers(std::map<SharedObjectNode*,SharedObject>& already_copied);
  
  /// Get the Jacobian
  virtual FX getJacobian() = 0;
  
  /// Get the Linear solver
  virtual LinearSolver getLinearSolver() = 0;

  /** \brief Create an integrator which integrates the ODE/DAE augmented with the forward sensitivity equations */
  virtual SundialsIntegrator jac(bool with_x, bool with_p);

  /** \brief Calculate the jacobian of a number of function outputs with respect to a number of function inputs, optionally include the function outputs */
  virtual FX jacobian(const std::vector<std::pair<int,int> >& jblocks);

  /// Generate the sparsity of a Jacobian block
  virtual CRSSparsity getJacSparsity(int iind, int oind);

  /** \brief  Set stop time for the integration */
  virtual void setStopTime(double tf) = 0;
  
  /// Set initial time (to be removed)
  void setInitialTime(double t0);

  /// Set final time (to be removed)
  void setFinalTime(double tf);

  /// Jacobian of the ODE/DAE with respect to the state and state derivatives (to be removed)
  FX jac_;

  /// Linear solver
  LinearSolver linsol_;
  
  /// Backwards integration function (to be replaced by g_)
  FX g_new_;
  
  //@{
  /// options
  bool exact_jacobian_;
  double abstol_, reltol_;
  double fsens_abstol_, fsens_reltol_;
  double asens_abstol_, asens_reltol_;
  int max_num_steps_;
  bool finite_difference_fsens_;  
  bool stop_at_end_;
  //@}
  
  /// Current time (to be removed)
  double t_;
  
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
  
  /// Use preconditioning
  bool use_preconditioner_;
  
};
  
} // namespace Sundials
} // namespace CasADi

#endif // SUNDIALS_INTERNAL_HPP

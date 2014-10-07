/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A symbolic framework for dynamic optimization.
 *    Copyright (C) 2010-2014 Joel Andersson, Joris Gillis, Moritz Diehl,
 *                            K.U. Leuven. All rights reserved.
 *    Copyright (C) 2011-2014 Greg Horn
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


#ifndef GSL_INTERNAL_HPP
#define GSL_INTERNAL_HPP

#include "gsl_integrator.hpp"
#include "core/function/integrator_internal.hpp"
#include "core/function/linear_solver.hpp"
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>

#include <ctime>

namespace casadi{
namespace GSL {
  
/**
@copydoc DAE_doc
*/
class GslInternal : public IntegratorInternal{
  friend class GslIntegrator;
public:
  /** \brief  Constructor */
  explicit GslInternal(const Function& f, const Function& q);

  /** \brief  Clone */
  virtual GslInternal* clone() const;
  
  /** \brief  Create a new integrator */
  virtual GslInternal* create(const Function& f, const Function& q) const{ return new GslInternal(f,q);}

  /** \brief  Destructor */
  virtual ~GslInternal();

  /** \brief  Initialize stage */
  virtual void init();

  /** \brief  Reset the solver and bring the time back to t0 */
  virtual void reset(int fsens_order, int asens_order);

  /** \brief  Integrate until a specified time point */
  virtual void integrate(double t_out);

  /** \brief  Set the stop time of the forward integration */
  virtual void setStopTime(double tf);
  
  /** \brief  Reset the solver of the adjoint problem and take time to tf */
  virtual void resetAdj();

  /** \brief  Integrate backwards in time until a specified time point */
  virtual void integrateAdj(double t_out);

  /// Nothing to see here
  virtual Function getJacobian();
  
  /// Get the Linear solver
  virtual LinearSolver getLinearSolver();

  // Set linear solver
  virtual void setLinearSolver(const LinearSolver& linsol, const Function& jac);
  
  protected:

  // GSL callback functions
  void rhs (double t, const double y[], double f[]);
  void jac (double t, const double y[], double *dfdy, double dfdt[]);
  
  // Static wrappers to be passed to GSL
  static int rhs_wrapper (double t, const double y[], double f[], void *userdata);
  static int jac_wrapper (double t, const double y[], double *dfdy, double dfdt[], void *userdata);
  
  // Current state of the integrator
  DMatrix state;
  
  virtual void printStats(std::ostream &stream) const;
  
  // Throw error
  static void gsl_error(const std::string& module, int flag);
  
  bool is_init;
  
  // Calculate the error message map
  static std::map<int,std::string> calc_flagmap();
  
  // Error message map
  static std::map<int,std::string> flagmap;
  
  // The jacobian of the ODE rhs fcn
  Function jac_f_;
  
  // The time derivative of the rhs fcn
  Function dt_f_;
  
  // For timings
  clock_t time1, time2;
  
  // Accummulated time since last reset:
  double t_res; // time spent in the DAE residual
  double t_fres; // time spent in the forward sensitivity residual
  double t_jac; // time spent in the jacobian, or jacobian times vector function
  double t_lsolve; // preconditioner/linear solver solve function
  double t_lsetup_jac; // preconditioner/linear solver setup function, generate jacobian
  double t_lsetup_fac; // preconditioner setup function, factorize jacobian

  // Type of integrator
  const gsl_odeiv_step_type *type_ptr;
  
  gsl_odeiv_step *step_ptr;
  gsl_odeiv_control *control_ptr;
  gsl_odeiv_evolve *evolve_ptr;

  gsl_odeiv_system my_system;	// structure with the rhs function, etc. 
  
  
};


} // namespace GSL
} // namespace casadi

#endif //GSL_INTERNAL_HPP


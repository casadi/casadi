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

#ifndef CVODES_INTERNAL_HPP
#define CVODES_INTERNAL_HPP

#include "cvodes_integrator.hpp"
#include "sundials_internal.hpp"
#include "casadi/fx/linear_solver.hpp"
#include <cvodes/cvodes.h>            /* prototypes for CVode fcts. and consts. */
#include <cvodes/cvodes_dense.h>
#include <cvodes/cvodes_band.h> 
#include <cvodes/cvodes_spgmr.h>
#include <cvodes/cvodes_spbcgs.h>
#include <cvodes/cvodes_sptfqmr.h>
#include <cvodes/cvodes_impl.h> /* Needed for the provided linear solver */
#include <ctime>

namespace CasADi{
namespace Sundials{
  
/**
  @copydoc ODE_doc
*/
class CVodesInternal : public SundialsInternal{
  friend class CVodesIntegrator;
public:
  /** \brief  Constructor */
  explicit CVodesInternal(const FX& f, const FX& g);

  /** \brief  Deep copy data members */
  virtual void deepCopyMembers(std::map<SharedObjectNode*,SharedObject>& already_copied);

  /** \brief  Clone */
  virtual CVodesInternal* clone() const;
  
  /** \brief  Create a new integrator */
  virtual CVodesInternal* create(const FX& f, const FX& g) const{ return new CVodesInternal(f,g);}

  /** \brief  Destructor */
  virtual ~CVodesInternal();

  /** \brief  Free all CVodes memory */
  virtual void freeCVodes();

  /** \brief  Initialize stage */
  virtual void init();

  /** \brief  Update the number of sensitivity directions during or after initialization */
  virtual void updateNumSens(bool recursive);

  /** \brief Initialize the adjoint problem (can only be called after the first integration) */
  virtual void initAdj();

  /** \brief  Reset the forward problem and bring the time back to t0 */
  virtual void reset(int nsens, int nsensB, int nsensB_store);

  /** \brief  Reset the backward problem and take time to tf */
  virtual void resetB();

  /** \brief  Integrate forward until a specified time point */
  virtual void integrate(double t_out);

  /** \brief  Integrate backward until a specified time point */
  virtual void integrateB(double t_out);

  /** \brief  Set the stop time of the forward integration */
  virtual void setStopTime(double tf);

  /** \brief Get the jacobian in the nonlinear iteration
  * The result is an CasADi::FX mapping from CasADi::Sundials::MInput to CasADi::Sundials::MOutput
  */
  virtual FX getJacobian();
  
  /// Get the Linear solver
  virtual LinearSolver getLinearSolver();

  protected:

  // Sundials callback functions
  void rhs(double t, const double* x, double* xdot);
  void ehfun(int error_code, const char *module, const char *function, char *msg);
  void rhsS(int Ns, double t, N_Vector x, N_Vector xdot, N_Vector *xF, N_Vector *xdotF, N_Vector tmp1, N_Vector tmp2);
  void rhsS1(int Ns, double t, N_Vector x, N_Vector xdot, int iS, N_Vector xF, N_Vector xdotF, N_Vector tmp1, N_Vector tmp2);
  void rhsQ(double t, const double* x, double* qdot);
  void rhsQS(int Ns, double t, N_Vector x, N_Vector *xF, N_Vector qdot, N_Vector *qFdot, N_Vector tmp1, N_Vector tmp2);
  void rhsB(double t, const double* x, const double *rx, double* rxdot);
  void rhsBS(double t, N_Vector x, N_Vector *xF, N_Vector xA, N_Vector xdotA);
  void rhsQB(double t, const double* x, const double* rx, double* rqdot);
  void jtimes(const double *v, double* Jv, double t, const double* x, const double* fy, double* tmp);
  void djac(SUNDIALS_INT N, double t, N_Vector x, N_Vector fy, DlsMat Jac, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
  void bjac(SUNDIALS_INT N, SUNDIALS_INT mupper, SUNDIALS_INT mlower, double t, N_Vector x, N_Vector fy, DlsMat Jac, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
  void psolve(double t, N_Vector x, N_Vector fy, N_Vector r, N_Vector z, double gamma, double delta, int lr, N_Vector tmp);
  void psetup(double t, N_Vector x, N_Vector fy, booleantype jok, booleantype *jcurPtr, double gamma, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
  void lsetup(CVodeMem cv_mem, int convfail, N_Vector ypred, N_Vector fpred, booleantype *jcurPtr, N_Vector vtemp1, N_Vector vtemp2, N_Vector vtemp3);
  void lsolve(CVodeMem cv_mem, N_Vector b, N_Vector weight, N_Vector ycur, N_Vector fcur);
  
  // Static wrappers to be passed to Sundials
  static int rhs_wrapper(double t, N_Vector x, N_Vector xdot, void *user_data);
  static void ehfun_wrapper(int error_code, const char *module, const char *function, char *msg, void *eh_data);
  static int rhsS_wrapper(int Ns, double t, N_Vector x, N_Vector xdot, N_Vector *xF, N_Vector *xdotF, void *user_data, N_Vector tmp1, N_Vector tmp2);
  static int rhsS1_wrapper(int Ns, double t, N_Vector x, N_Vector xdot, int iS, N_Vector xF, N_Vector xdotF, void *user_data, N_Vector tmp1, N_Vector tmp2);
  static int rhsQ_wrapper(double t, N_Vector x, N_Vector qdot, void *user_data);
  static int rhsQS_wrapper(int Ns, double t, N_Vector x, N_Vector *xF, N_Vector qdot, N_Vector *qdotF, void *user_data, N_Vector tmp1, N_Vector tmp2);
  static int rhsB_wrapper(double t, N_Vector x, N_Vector xA, N_Vector xdotA, void *user_data);
  static int rhsBS_wrapper(double t, N_Vector x, N_Vector *xF, N_Vector xA, N_Vector xdotA, void *user_data);
  static int rhsQB_wrapper(double t, N_Vector x, N_Vector xA, N_Vector qdotA, void *user_data);
  static int jtimes_wrapper(N_Vector v, N_Vector Jv, double t, N_Vector x, N_Vector xdot, void *user_data, N_Vector tmp);
  static int djac_wrapper(SUNDIALS_INT N, double t, N_Vector x, N_Vector xdot, DlsMat Jac, void *user_data,N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
  static int bjac_wrapper(SUNDIALS_INT N, SUNDIALS_INT mupper, SUNDIALS_INT mlower, double t, N_Vector x, N_Vector xdot, DlsMat Jac, void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
  static int psolve_wrapper(double t, N_Vector x, N_Vector xdot, N_Vector r, N_Vector z, double gamma, double delta, int lr, void *user_data, N_Vector tmp);
  static int psetup_wrapper(double t, N_Vector x, N_Vector xdot, booleantype jok, booleantype *jcurPtr, double gamma, void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
  static int lsetup_wrapper(CVodeMem cv_mem, int convfail, N_Vector x, N_Vector xdot, booleantype *jcurPtr, N_Vector vtemp1, N_Vector vtemp2, N_Vector vtemp3);
  static int lsolve_wrapper(CVodeMem cv_mem, N_Vector b, N_Vector weight, N_Vector x, N_Vector xdot);
  
  virtual void printStats(std::ostream &stream) const;
  
  // CVodes memory block
  void* mem_;
  
  // The jacobian of the ODE rhs fcn
  FX jac_f_;
  
  // For timings
  clock_t time1, time2;
  
  // Accummulated time since last reset:
  double t_res; // time spent in the DAE residual
  double t_fres; // time spent in the forward sensitivity residual
  double t_jac; // time spent in the jacobian, or jacobian times vector function
  double t_lsolve; // preconditioner/linear solver solve function
  double t_lsetup_jac; // preconditioner/linear solver setup function, generate jacobian
  double t_lsetup_fac; // preconditioner setup function, factorize jacobian
  
  // N-vectors for the forward integration
  N_Vector x0_, x_, q_;
  
  // N-vectors for the backward integration
  N_Vector rx0_, rx_, rq_;

  // N-vectors for the forward sensitivities
  std::vector<N_Vector> xF0_, xF_, qF_;

  bool is_init;
  bool isInitAdj_;

  int ism_;
  
  // Calculate the error message map
  static std::map<int,std::string> calc_flagmap();
  
    // Error message map
  static std::map<int,std::string> flagmap;
 
  // Throw error
  static void cvodes_error(const std::string& module, int flag);

  // Ids of backward problem
  int whichB_;

  // Number of forward directions for the functions f and g
  int nfdir_f_, nfdir_g_;

  // Initialize the dense linear solver
  void initDenseLinearSolver();
  
  // Initialize the banded linear solver
  void initBandedLinearSolver();
  
  // Initialize the iterative linear solver
  void initIterativeLinearSolver();
  
  // Initialize the user defined linear solver
  void initUserDefinedLinearSolver();
  
  // Initialize the dense linear solver (backward integration)
  void initDenseLinearSolverB();
  
  // Initialize the banded linear solver (backward integration)
  void initBandedLinearSolverB();
  
  // Initialize the iterative linear solver (backward integration)
  void initIterativeLinearSolverB();
  
  // Initialize the user defined linear solver (backward integration)
  void initUserDefinedLinearSolverB();

  // Set linear solver
  virtual void setLinearSolver(const LinearSolver& linsol, const FX& jac);

  int lmm_; // linear multistep method
  int iter_; // nonlinear solver iteration

  bool monitor_rhsB_;
  bool monitor_rhs_;
  bool monitor_rhsQB_;
  
  bool disable_internal_warnings_;
  
};


} // namespace Sundials
} // namespace CasADi

#endif //CVODES_INTERNAL_HPP


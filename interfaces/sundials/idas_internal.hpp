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

#ifndef IDAS_INTERNAL_HPP
#define IDAS_INTERNAL_HPP

#include "idas_integrator.hpp"
#include "sundials_internal.hpp"
#include "symbolic/fx/linear_solver.hpp"
#include <idas/idas.h>            /* prototypes for CVODE fcts. and consts. */
#include <idas/idas_dense.h>
#include <idas/idas_band.h> 
#include <idas/idas_spgmr.h>
#include <idas/idas_spbcgs.h>
#include <idas/idas_sptfqmr.h>
#include <idas/idas_impl.h> /* Needed for the provided linear solver */
#include <ctime>

namespace CasADi{
  
/**
@copydoc IdasIntegrator_doc
*/
class IdasInternal : public SundialsInternal{
  friend class IdasIntegrator;

  public:
  
  /** \brief  Constructor */
  explicit IdasInternal(const FX& f, const FX& g);

  /** \brief  Copy constructor */
//  IdasInternal(const IdasInternal& integrator);

  /** \brief  Clone */
  virtual IdasInternal* clone() const;
  
  /** \brief  Create a new integrator */
  virtual IdasInternal* create(const FX& f, const FX& g) const{ return new IdasInternal(f,g);}

  /** \brief  Deep copy data members */
  virtual void deepCopyMembers(std::map<SharedObjectNode*,SharedObject>& already_copied);

  /** \brief  Destructor */
  virtual ~IdasInternal();

  /** \brief  Free all IDAS memory */
  virtual void freeIDAS();

  /** \brief  Initialize */
  virtual void init();
  
  /** \brief  Update the number of sensitivity directions during or after initialization */
  virtual void updateNumSens(bool recursive);

  /** \brief Initialize the taping */
  virtual void initTaping();
  
  /** \brief Initialize the backward problem (can only be called after the first integration) */
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
  
  /** \brief  Print solver statistics */  
  virtual void printStats(std::ostream &stream) const;
  
  /** \brief  Get the integrator Jacobian for the forward problem */
  FX getJacobian();
  
  /** \brief  Get the integrator Jacobian for the forward problem (generic) */
  template<typename FunctionType>
  FunctionType getJacobianGen();
  
  /** \brief  Get the integrator Jacobian for the backward problem */
  FX getJacobianB();
  
  /** \brief  Get the integrator Jacobian for the backward problem (generic) */
  template<typename FunctionType>
  FunctionType getJacobianGenB();
  
  /// Get the Linear solver
  virtual LinearSolver getLinearSolver();

  /// Correct the initial conditions, i.e. calculate
  void correctInitialConditions();
  
  protected:

  // Sundials callback functions
  void res(double t, const double* xz, const double* xzdot, double* rr);
  void ehfun(int error_code, const char *module, const char *function, char *msg);
  void jtimes(double t, const double *xz, const double *xzdot, const double *rr, const double *v, double *Jv, double cj, double *tmp1, double *tmp2);
  void resS(int Ns, double t, const double* xz, const double* xzdot, const double *resval, N_Vector *xzF, N_Vector* xzdotF, N_Vector *rrF, double *tmp1, double *tmp2, double *tmp3);
  void rhsQ(double t, const double* xz, const double* xzdot, double* qdot);
  void rhsQS(int Ns, double t, N_Vector xz, N_Vector xzdot, N_Vector *xzF, N_Vector *xzdotF, N_Vector rrQ, N_Vector *qdotF, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
  void resB(double t, const double* y, const double* xzdot, const double* xA, const double* xzdotA, double* rrA);
  void rhsQB(double t, const double* y, const double* xzdot, const double* xA, const double* xzdotA, double *qdotA);
  void psolve(double t, N_Vector xz, N_Vector xzdot, N_Vector rr, N_Vector rvec, N_Vector zvec, double cj, double delta, N_Vector tmp);
  void psetup(double t, N_Vector xz, N_Vector xzdot, N_Vector rr, double cj, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
  void djac(long Neq, double t, double cj, N_Vector xz, N_Vector xzdot, N_Vector rr, DlsMat Jac, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
  void bjac(long Neq, long mupper, long mlower, double tt, double cj, N_Vector xz, N_Vector xzdot, N_Vector rr, DlsMat Jac, N_Vector tmp1, N_Vector tmp2,N_Vector tmp3);
  void lsetup(IDAMem IDA_mem, N_Vector xzp, N_Vector xzdotp, N_Vector resp, N_Vector vtemp1, N_Vector vtemp2, N_Vector vtemp3);
  void lsolve(IDAMem IDA_mem, N_Vector b, N_Vector weight, N_Vector xzcur, N_Vector xzdotcur, N_Vector rescur);

  // Static wrappers to be passed to Sundials
  static int res_wrapper(double t, N_Vector xz, N_Vector xzdot, N_Vector rr, void *user_data);
  static void ehfun_wrapper(int error_code, const char *module, const char *function, char *msg, void *eh_data);
  static int jtimes_wrapper(double t, N_Vector xz, N_Vector xzdot, N_Vector rr, N_Vector v, N_Vector Jv, double cj, void *user_data, N_Vector tmp1, N_Vector tmp2);
  static int resS_wrapper(int Ns, double t, N_Vector xz, N_Vector xzdot, N_Vector resval, N_Vector *xzF, N_Vector *xzdotF, N_Vector *resF, void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
  static int rhsQ_wrapper(double t, N_Vector xz, N_Vector xzdot, N_Vector qdot, void *user_data);
  static int rhsQS_wrapper(int Ns, double t, N_Vector xz, N_Vector xzdot, N_Vector *xzF, N_Vector *xzdotF, N_Vector rrQ, N_Vector *qdotF, void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
  static int resB_wrapper(double t, N_Vector y, N_Vector xzdot, N_Vector xA, N_Vector xzdotA, N_Vector resA, void *user_dataB);
  static int rhsQB_wrapper(double t, N_Vector y, N_Vector xzdot, N_Vector xA, N_Vector xzdotA, N_Vector qdotA, void *user_dataB);
  static int psolve_wrapper(double t, N_Vector xz, N_Vector xzdot, N_Vector rr, N_Vector rvec, N_Vector zvec, double cj, double delta, void *user_data, N_Vector tmp);
  static int psetup_wrapper(double t, N_Vector xz, N_Vector xzdot, N_Vector rr, double cj, void* user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
  static int djac_wrapper(long Neq, double t, double cj, N_Vector xz, N_Vector xzdot, N_Vector rr, DlsMat Jac, void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
  static int bjac_wrapper(long Neq, long mupper, long mlower, double tt, double cj, N_Vector xz, N_Vector xzdot, N_Vector rr, DlsMat Jac, void *user_data, N_Vector tmp1, N_Vector tmp2,N_Vector tmp3);
  static int lsetup_wrapper(IDAMem IDA_mem, N_Vector xz, N_Vector xzdot, N_Vector resp, N_Vector vtemp1, N_Vector vtemp2, N_Vector vtemp3);
  static int lsolve_wrapper(IDAMem IDA_mem, N_Vector b, N_Vector weight, N_Vector ycur, N_Vector xzdotcur, N_Vector rescur);
  
 public:

  // Idas memory block
  void* mem_;

  // N-vectors for the forward integration
  N_Vector xz_, xzdot_, q_;
  
  // N-vectors for the backward integration
  N_Vector rxz_, rxzdot_, rq_;

  // N-vectors for the forward sensitivities
  std::vector<N_Vector> xzF_, xzdotF_, qF_;

  // sensitivity method
  int ism_;
  
  // Calculate the error message map
  static std::map<int,std::string> calc_flagmap();
  
  // Error message map
  static std::map<int,std::string> flagmap;
 
  // Throw error
  static void idas_error(const std::string& module, int flag);
  
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
  
  // Ids of backward problem
  int whichB_;

  // For timings
  clock_t time1, time2;
  
  // Accummulated time since last reset:
  double t_res; // time spent in the DAE residual
  double t_fres; // time spent in the forward sensitivity residual
  double t_jac; // time spent in the jacobian, or jacobian times vector function
  double t_lsolve; // preconditioner/linear solver solve function
  double t_lsetup_jac; // preconditioner/linear solver setup function, generate jacobian
  double t_lsetup_fac; // preconditioner setup function, factorize jacobian
    
  // Has the adjoint problem been initialized
  bool isInitAdj_;
  bool isInitTaping_;
  
  // Number of forward seeds for the function f
  int nfdir_f_;
  
  // Scaling of cj
  bool cj_scaling_;

  // Set linear solver
  virtual void setLinearSolver(const LinearSolver& linsol, const FX& jac);
  
  // Disable IDAS internal warning messages
  bool disable_internal_warnings_;
  
  //  Initial values for xdot and z
  std::vector<double> init_z_, init_xdot_;
  
  // Jacobian of the DAE with respect to the state and state derivatives
  FX jac_, jacB_;
};

} // namespace CasADi

#endif //IDAS_INTERNAL_HPP


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
#include "casadi/fx/integrator_internal.hpp"
#include "casadi/fx/linear_solver.hpp"
#include <nvector/nvector_serial.h>   /* serial N_Vector types, fcts., and macros */
#include <sundials/sundials_dense.h>  /* definitions DlsMat DENSE_ELEM */
#include <sundials/sundials_types.h>  /* definition of type double */
#include <idas/idas.h>            /* prototypes for CVODE fcts. and consts. */
#include <idas/idas_dense.h>
#include <idas/idas_band.h> 
#include <idas/idas_spgmr.h>
#include <idas/idas_spbcgs.h>
#include <idas/idas_sptfqmr.h>
#include <idas/idas_impl.h> /* Needed for the provided linear solver */
#include <ctime>

namespace CasADi{
namespace Sundials{
  
class IdasInternal : public IntegratorInternal{
  friend class IdasIntegrator;

  public:
  
  /** \brief  Constructor */
  explicit IdasInternal(const FX& f, const FX& q);

  /** \brief  Clone */
  virtual IdasInternal* clone() const;

  /** \brief  Destructor */
  virtual ~IdasInternal();

  /** \brief  Initialize */
  virtual void init();

  /** \brief Initialize the adjoint problem (can only be called after the first integration) */
  virtual void initAdj();
  
  /** \brief  Reset the solver and bring the time back to t0 */
  virtual void reset(int fsens_order, int asens_order);

  /** \brief  Reset the solver of the adjoint problem and take time to tf */
  virtual void resetAdj();

  /** \brief  Integrate until a specified time point */
  virtual void integrate(double t_out);

  /** \brief  Integrate backwards in time until a specified time point */
  virtual void integrateAdj(double t_out);

  /** \brief  Set the stop time of the forward integration */
  virtual void setStopTime(double tf);
  
  /** \brief  Print solver statistics */  
  virtual void printStats(std::ostream &stream) const;
  
  protected:

  // Sundials callback functions
  void res(double t, const double* yy, const double* yp, double* rr);
  void ehfun(int error_code, const char *module, const char *function, char *msg);
  void jtimes(double t, const double *yy, const double *yp, const double *rr, const double *v, double *Jv, double cj, double *tmp1, double *tmp2);
  void resS(int Ns, double t, const double* yy, const double* yp, const double *resval, N_Vector *yS, N_Vector* ypS, N_Vector *resvalS, double *tmp1, double *tmp2, double *tmp3);
  void rhsQ(double t, const double* yy, const double* yp, double* rhsQ);
  void rhsQS(double t, const double* yy, const double* yp, const double *yyS, const double *ypS, const double* rrQ, double* rhsvalQS, double* tmp1, double* tmp2, double* tmp3);
  void rhsQS(int Ns, double t, N_Vector yy, N_Vector yp, N_Vector *yyS, N_Vector *ypS, N_Vector rrQ, N_Vector *rhsvalQS, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
  void resB(double t, const double* y, const double* yp, const double* yB, const double* ypB, double* resvalB);
  void rhsQB(double t, const double* y, const double* yp, const double* yB, const double* ypB, double *rhsvalBQ);
  void psolve(double t, N_Vector yy, N_Vector yp, N_Vector rr, N_Vector rvec, N_Vector zvec, double cj, double delta, N_Vector tmp);
  void psetup(double t, N_Vector yy, N_Vector yp, N_Vector rr, double cj, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
  void djac(int Neq, double t, double cj, N_Vector yy, N_Vector yp, N_Vector rr, DlsMat Jac, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
  void bjac(int Neq, int mupper, int mlower, double tt, double cj, N_Vector yy, N_Vector yp, N_Vector rr, DlsMat Jac, N_Vector tmp1, N_Vector tmp2,N_Vector tmp3);
  void lsetup(IDAMem IDA_mem, N_Vector yyp, N_Vector ypp, N_Vector resp, N_Vector vtemp1, N_Vector vtemp2, N_Vector vtemp3);
  void lsolve(IDAMem IDA_mem, N_Vector b, N_Vector weight, N_Vector ycur, N_Vector ypcur, N_Vector rescur);

  // Static wrappers to be passed to Sundials
  static int res_wrapper(double t, N_Vector yy, N_Vector yp, N_Vector rr, void *user_data);
  static void ehfun_wrapper(int error_code, const char *module, const char *function, char *msg, void *eh_data);
  static int jtimes_wrapper(double t, N_Vector yy, N_Vector yp, N_Vector rr, N_Vector v, N_Vector Jv, double cj, void *user_data, N_Vector tmp1, N_Vector tmp2);
  static int resS_wrapper(int Ns, double t, N_Vector yy, N_Vector yp, N_Vector resval, N_Vector *yS, N_Vector *ypS, N_Vector *resvalS, void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
  static int rhsQ_wrapper(double t, N_Vector yy, N_Vector yp, N_Vector rhsQ, void *user_data);
  static int rhsQS_wrapper(int Ns, double t, N_Vector yy, N_Vector yp, N_Vector *yyS, N_Vector *ypS, N_Vector rrQ, N_Vector *rhsvalQS, void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
  static int resB_wrapper(double t, N_Vector y, N_Vector yp, N_Vector yB, N_Vector ypB, N_Vector resvalB, void *user_dataB);
  static int rhsQB_wrapper(double t, N_Vector y, N_Vector yp, N_Vector yB, N_Vector ypB, N_Vector rhsvalBQ, void *user_dataB);
  static int psolve_wrapper(double t, N_Vector yy, N_Vector yp, N_Vector rr, N_Vector rvec, N_Vector zvec, double cj, double delta, void *user_data, N_Vector tmp);
  static int psetup_wrapper(double t, N_Vector yy, N_Vector yp, N_Vector rr, double cj, void* user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
  static int djac_wrapper(int Neq, double t, double cj, N_Vector yy, N_Vector yp, N_Vector rr, DlsMat Jac, void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
  static int bjac_wrapper(int Neq, int mupper, int mlower, double tt, double cj, N_Vector yy, N_Vector yp, N_Vector rr, DlsMat Jac, void *user_data, N_Vector tmp1, N_Vector tmp2,N_Vector tmp3);
  static int lsetup_wrapper(IDAMem IDA_mem, N_Vector yyp, N_Vector ypp, N_Vector resp, N_Vector vtemp1, N_Vector vtemp2, N_Vector vtemp3);
  static int lsolve_wrapper(IDAMem IDA_mem, N_Vector b, N_Vector weight, N_Vector ycur, N_Vector ypcur, N_Vector rescur);
  
 public:

  // Idas memory block
  void* mem_;

  // ODE rhs
  FX f_;

  // Quadrature function
  FX q_;

  // N-vectors for the DAE integration
  N_Vector  y0_,  y_;
  N_Vector yp0_, yp_; 
  N_Vector yQ0_, yQ_;

  // N-vectors for the forward and adjoint sensitivities
  std::vector<N_Vector> yS0_, yS_, ypS0_, ypS_, yQS0_, yQS_;

  // N-vectors for the adjoint sensitivities
  std::vector<N_Vector> yB0_, yB_, ypB0_, ypB_, /*yQB0_, */ yQB_;

  // dimensions
  int ny_; // number of ode states
  int nq_; // number of quadratures
  
  bool is_init;
  
  // sensitivity method
  int ism_;
  
  // Calculate the error message map
  static map<int,std::string> calc_flagmap();
  
  // Error message map
  static std::map<int,std::string> flagmap;
 
  // Throw error
  static void idas_error(const std::string& module, int flag);
  
  // Auxiliary
  static int getNX(const FX& f, const FX& q); // count the total number of states
  static int getNP(const FX& f); // count the number of parameters

  // Set the user defined linear solver
  void initUserDefinedLinearSolver();
  
  // Ids of backward problem
  std::vector<int> whichB_;

  int fsens_order_, asens_order_;
  
  // Calculate consistent initial conditions
  bool calc_ic_;
  
  // Jacobian of the ODE with respect to the state and state derivatives (separately)
  FX jacx_, jacxdot_;
  
  // Jacobian of the ODE with respect to the state and state derivatives
  FX jac_;

  // Jacobian of the ODE with respect to the parameters
  FX jacp_;
  
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
  
  // Number of forward and adjoint seeds for the functions f and q
  int nfdir_f_, nadir_f_, nfdir_q_, nadir_q_;
  
  // Linear solver
  LinearSolver linsol_;  

  // Can calcic be used?
  bool calc_ic_ok_;
  
  // Scaling of cj
  bool cj_scaling_;
};


} // namespace Sundials
} // namespace CasADi

#endif //IDAS_INTERNAL_HPP


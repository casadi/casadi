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

/*
TODO:
  To make the class useful in a dynamic optimization setting, it is necessary to make a clear division of the
  state vector into a differential and an algebraic part. The differential part contains all the states that
  appear differentiated in the DAE.
  
  The following changes will apper:
  
  * The algebraic state will be included into DAEInput as DAS_Z and JACInput as JAC_Z, user functions must be changed accordingly
  * The option "is_differential" is no longer needed and should be removed
  * IdasInternal will allocate its own state vectors using N_VNew_Serial instead of N_VMake_Serial
  * A user interface must be provided to allow access to the algebraic state as well as the time derivative of the differential state,
    as these will no longer accessable via getInput/setInput
  
  Open questions to be answered before programming starts:
  * Should the differential state and algebraic state be forced to follow in that exact order? If so, how to get e.g. banded Jacobians to work well? 
    If not, how to provide an ordering?
  * How many state vectors should be stored, one or two?
  * Where to store the consistent initial values for the forward and backward integration?
    - In the integrator class: Good if the integrator will be called several times with similar arguments, e.g. one integrator per shooting node
    - Outside the integrator: Good if the integrator will be called with very different arguments, e.g. one integrator for all shooting nodes
    
  If outside the integrator, where to put it?
    - In a new class? If so, all the changes above might not be needed, in particular internal allocation (and copying!) of the state vector
    Is this a problem which can be relevant to explicit integrators as well?
    - In the simulator class? If so, need to give the simulator class extra functionality, in particular dealing with controls
    
  Conclusion: The most user-friendly alternative is probably to move this functionality to the simulator class, keeping integrator as low level as possible
  This will also probably save a lot of memory, as only one integrator needs to be allocated. On the other hand, on parallell machines, multiple
  integrators are necessary for load balancing - and in this case, it is better not to use the simulator at all.  

*/



namespace CasADi{
namespace Sundials{
  
class IdasInternal : public IntegratorInternal{
  friend class IdasIntegrator;

  public:
  
  /** \brief  Constructor */
  explicit IdasInternal(const FX& f, const FX& q);

  /** \brief  Copy constructor */
//  IdasInternal(const IdasInternal& integrator);

  /** \brief  Clone */
  virtual IdasInternal* clone() const;

  /** \brief  Destructor */
  virtual ~IdasInternal();

  /** \brief  Initialize */
  virtual void init();

  /** \brief Initialize the adjoint problem (can only be called after the first integration) */
  virtual void initAdj();
  
  /** \brief  Reset the solver and bring the time back to t0 */
  virtual void reset(int nfsens=0, int nasens=0);

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
  
  /** \brief Create an integrator which integrates the ODE/DAE augmented with the forward sensitivity equations */
  virtual Integrator jac(int iind=0, int oind=0);

  /// Get the Jacobian
  virtual FX getJacobian();
  
  /// Get the Linear solver
  virtual LinearSolver getLinearSolver();

  protected:

  // Sundials callback functions
  void res(double t, const double* yz, const double* yp, double* rr);
  void ehfun(int error_code, const char *module, const char *function, char *msg);
  void jtimes(double t, const double *yz, const double *yp, const double *rr, const double *v, double *Jv, double cj, double *tmp1, double *tmp2);
  void resS(int Ns, double t, const double* yz, const double* yp, const double *resval, N_Vector *yS, N_Vector* ypS, N_Vector *resvalS, double *tmp1, double *tmp2, double *tmp3);
  void rhsQ(double t, const double* yz, const double* yp, double* rhsQ);
  void rhsQS(double t, const double* yz, const double* yp, const double *yzS, const double *ypS, const double* rrQ, double* rhsvalQS, double* tmp1, double* tmp2, double* tmp3);
  void rhsQS(int Ns, double t, N_Vector yz, N_Vector yp, N_Vector *yzS, N_Vector *ypS, N_Vector rrQ, N_Vector *rhsvalQS, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
  void resB(double t, const double* y, const double* yp, const double* yB, const double* ypB, double* resvalB);
  void rhsQB(double t, const double* y, const double* yp, const double* yB, const double* ypB, double *rhsvalBQ);
  void psolve(double t, N_Vector yz, N_Vector yp, N_Vector rr, N_Vector rvec, N_Vector zvec, double cj, double delta, N_Vector tmp);
  void psetup(double t, N_Vector yz, N_Vector yp, N_Vector rr, double cj, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
  void djac(int Neq, double t, double cj, N_Vector yz, N_Vector yp, N_Vector rr, DlsMat Jac, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
  void bjac(int Neq, int mupper, int mlower, double tt, double cj, N_Vector yz, N_Vector yp, N_Vector rr, DlsMat Jac, N_Vector tmp1, N_Vector tmp2,N_Vector tmp3);
  void lsetup(IDAMem IDA_mem, N_Vector yzp, N_Vector ypp, N_Vector resp, N_Vector vtemp1, N_Vector vtemp2, N_Vector vtemp3);
  void lsolve(IDAMem IDA_mem, N_Vector b, N_Vector weight, N_Vector ycur, N_Vector ypcur, N_Vector rescur);

  // Static wrappers to be passed to Sundials
  static int res_wrapper(double t, N_Vector yz, N_Vector yp, N_Vector rr, void *user_data);
  static void ehfun_wrapper(int error_code, const char *module, const char *function, char *msg, void *eh_data);
  static int jtimes_wrapper(double t, N_Vector yz, N_Vector yp, N_Vector rr, N_Vector v, N_Vector Jv, double cj, void *user_data, N_Vector tmp1, N_Vector tmp2);
  static int resS_wrapper(int Ns, double t, N_Vector yz, N_Vector yp, N_Vector resval, N_Vector *yS, N_Vector *ypS, N_Vector *resvalS, void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
  static int rhsQ_wrapper(double t, N_Vector yz, N_Vector yp, N_Vector rhsQ, void *user_data);
  static int rhsQS_wrapper(int Ns, double t, N_Vector yz, N_Vector yp, N_Vector *yzS, N_Vector *ypS, N_Vector rrQ, N_Vector *rhsvalQS, void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
  static int resB_wrapper(double t, N_Vector y, N_Vector yp, N_Vector yB, N_Vector ypB, N_Vector resvalB, void *user_dataB);
  static int rhsQB_wrapper(double t, N_Vector y, N_Vector yp, N_Vector yB, N_Vector ypB, N_Vector rhsvalBQ, void *user_dataB);
  static int psolve_wrapper(double t, N_Vector yz, N_Vector yp, N_Vector rr, N_Vector rvec, N_Vector zvec, double cj, double delta, void *user_data, N_Vector tmp);
  static int psetup_wrapper(double t, N_Vector yz, N_Vector yp, N_Vector rr, double cj, void* user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
  static int djac_wrapper(int Neq, double t, double cj, N_Vector yz, N_Vector yp, N_Vector rr, DlsMat Jac, void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
  static int bjac_wrapper(int Neq, int mupper, int mlower, double tt, double cj, N_Vector yz, N_Vector yp, N_Vector rr, DlsMat Jac, void *user_data, N_Vector tmp1, N_Vector tmp2,N_Vector tmp3);
  static int lsetup_wrapper(IDAMem IDA_mem, N_Vector yzp, N_Vector ypp, N_Vector resp, N_Vector vtemp1, N_Vector vtemp2, N_Vector vtemp3);
  static int lsolve_wrapper(IDAMem IDA_mem, N_Vector b, N_Vector weight, N_Vector ycur, N_Vector ypcur, N_Vector rescur);
  
 public:

  // Idas memory block
  void* mem_;

  // ODE rhs
  FX f_;

  // Quadrature function
  FX q_;

  // N-vectors for the DAE integration
  N_Vector  yz_, yP_, yQ_;

  // N-vectors for the forward and adjoint sensitivities
  std::vector<N_Vector> yzS_, yPS_, yQS_;

  // N-vectors for the adjoint sensitivities
  std::vector<N_Vector> yzB_, yPB_, yBB_;

  // Which components are differential
  N_Vector id_;
  
  // dimensions
  int ny_; // number of differential states
  int nq_; // number of quadratures
  int nyz_; // number of states seen by IDA (differential states and algebraic states)
  
  bool is_init;
  
  // sensitivity method
  int ism_;
  
  // Calculate the error message map
  static map<int,std::string> calc_flagmap();
  
  // Error message map
  static std::map<int,std::string> flagmap;
 
  // Throw error
  static void idas_error(const std::string& module, int flag);
  
  // Set the user defined linear solver
  void initUserDefinedLinearSolver();
  
  // Ids of backward problem
  std::vector<int> whichB_;

  int fsens_order_, asens_order_;
  
  // Jacobian of the ODE with respect to the state and state derivatives (separately)
  FX jacx_, jacxdot_, jacz_;
  
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
  
  // number of checkpoints stored so far
  int ncheck_; 
  
  // Has the adjoint problem been initialized
  bool isInitAdj_;
  
  // Number of forward and adjoint seeds for the functions f and q
  int nfdir_f_, nadir_f_, nfdir_q_, nadir_q_;
  
  // Linear solver
  LinearSolver linsol_;  

  // Scaling of cj
  bool cj_scaling_;

  // Set linear solver
  virtual void setLinearSolver(const LinearSolver& linsol, const FX& jac);
  
  // Get the initial state
  void getInitialState();
  
  // Set the final state
  void setFinalState();
  
  // Get the forward seeds
  void getForwardSeeds();
  
  // Set the forward sensitivities
  void setForwardSensitivities();
  
  // Get the adjoint seeds
  void getAdjointSeeds();
  
  // Set the adjoint sensitivities
  void setAdjointSensitivities();
  
};


} // namespace Sundials
} // namespace CasADi

#endif //IDAS_INTERNAL_HPP


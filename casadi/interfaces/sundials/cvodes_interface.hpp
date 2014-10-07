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


#ifndef CASADI_CVODES_INTERFACE_HPP
#define CASADI_CVODES_INTERFACE_HPP

#include <casadi/interfaces/sundials/casadi_integrator_cvodes_export.h>
#include "sundials_interface.hpp"
#include "casadi/core/function/linear_solver.hpp"
#include <cvodes/cvodes.h>            /* prototypes for CVode fcts. and consts. */
#include <cvodes/cvodes_dense.h>
#include <cvodes/cvodes_band.h>
#include <cvodes/cvodes_spgmr.h>
#include <cvodes/cvodes_spbcgs.h>
#include <cvodes/cvodes_sptfqmr.h>
#include <cvodes/cvodes_impl.h> /* Needed for the provided linear solver */
#include <ctime>

/** \defgroup plugin_Integrator_cvodes

      Interface to CVodes from the Sundials suite.

      A call to evaluate will integrate to the end.

      You can retrieve the entire state trajectory as follows, after the evaluate call:
      Call reset. Then call integrate(t_i) and getOuput for a series of times t_i.

      Note: depending on the dimension and structure of your problem,
      you may experience a dramatic speed-up by using a sparse linear solver:

      \verbatim
       intg.setOption("linear_solver","csparse")
       intg.setOption("linear_solver_type","user_defined")
      \endverbatim
*/

/** \pluginsection{Integrator,cvodes} */

/// \cond INTERNAL
namespace casadi {

  /** \brief \pluginbrief{Integrator,cvodes}

      @copydoc DAE_doc
      @copydoc plugin_Integrator_cvodes

  */
  class CASADI_INTEGRATOR_CVODES_EXPORT CvodesInterface : public SundialsInterface {
  public:
    /** \brief  Constructor */
    explicit CvodesInterface(const Function& f, const Function& g);

    /** \brief  Deep copy data members */
    virtual void deepCopyMembers(std::map<SharedObjectNode*, SharedObject>& already_copied);

    /** \brief  Clone */
    virtual CvodesInterface* clone() const;

    /** \brief  Create a new integrator */
    virtual CvodesInterface* create(const Function& f, const Function& g) const
    { return new CvodesInterface(f, g);}

    /** \brief  Create a new integrator */
    static IntegratorInternal* creator(const Function& f, const Function& g)
    { return new CvodesInterface(f, g);}

    /** \brief  Destructor */
    virtual ~CvodesInterface();

    /** \brief  Free all CVodes memory */
    virtual void freeCVodes();

    /** \brief  Initialize stage */
    virtual void init();

    /** \brief Initialize the adjoint problem (can only be called after the first integration) */
    virtual void initAdj();

    /** \brief  Reset the forward problem and bring the time back to t0 */
    virtual void reset();

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

    /** \brief  Get the integrator Jacobian for the forward problem (generic) */
    template<typename FunctionType>
      FunctionType getJacGen();

    /** \brief  Get the integrator Jacobian for the backward problem (generic) */
    template<typename FunctionType>
      FunctionType getJacGenB();

    /** \brief  Get the integrator Jacobian for the forward problem */
    virtual Function getJac();

    /** \brief  Get the integrator Jacobian for the backward problem */
    virtual Function getJacB();

    /// A documentation string
    static const std::string meta_doc;

  protected:

    // Sundials callback functions
    void rhs(double t, const double* x, double* xdot);
    void ehfun(int error_code, const char *module, const char *function, char *msg);
    void rhsS(int Ns, double t, N_Vector x, N_Vector xdot, N_Vector *xF, N_Vector *xdotF,
              N_Vector tmp1, N_Vector tmp2);
    void rhsS1(int Ns, double t, N_Vector x, N_Vector xdot, int iS, N_Vector xF, N_Vector xdotF,
               N_Vector tmp1, N_Vector tmp2);
    void rhsQ(double t, const double* x, double* qdot);
    void rhsQS(int Ns, double t, N_Vector x, N_Vector *xF, N_Vector qdot, N_Vector *qFdot,
               N_Vector tmp1, N_Vector tmp2);
    void rhsB(double t, const double* x, const double *rx, double* rxdot);
    void rhsBS(double t, N_Vector x, N_Vector *xF, N_Vector xB, N_Vector xdotB);
    void rhsQB(double t, const double* x, const double* rx, double* rqdot);
    void jtimes(N_Vector v, N_Vector Jv, double t, N_Vector x, N_Vector xdot, N_Vector tmp);
    void jtimesB(N_Vector vB, N_Vector JvB, double t, N_Vector x, N_Vector xB,
                 N_Vector xdotB, N_Vector tmpB);
    void djac(long N, double t, N_Vector x, N_Vector xdot, DlsMat Jac,
              N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
    void djacB(long NeqB, double t, N_Vector x, N_Vector xB, N_Vector xdotB, DlsMat JacB,
               N_Vector tmp1B, N_Vector tmp2B, N_Vector tmp3B);
    void bjac(long N, long mupper, long mlower, double t, N_Vector x, N_Vector xdot, DlsMat Jac,
              N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
    void bjacB(long NeqB, long mupperB, long mlowerB, double t, N_Vector x, N_Vector xB,
               N_Vector xdotB, DlsMat JacB, N_Vector tmp1B, N_Vector tmp2B, N_Vector tmp3B);
    /// <tt>z = M^(-1).r</tt>
    void psolve(double t, N_Vector x, N_Vector xdot, N_Vector r, N_Vector z,
                double gamma, double delta, int lr, N_Vector tmp);
    void psolveB(double t, N_Vector x, N_Vector xB, N_Vector xdotB, N_Vector rvecB, N_Vector zvecB,
                 double gammaB, double deltaB, int lr, N_Vector tmpB);
    /// <tt>M = I-gamma*df/dx</tt>, factorize
    void psetup(double t, N_Vector x, N_Vector xdot, booleantype jok, booleantype *jcurPtr,
                double gamma,
                N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
    void psetupB(double t, N_Vector x, N_Vector xB, N_Vector xdotB,
                 booleantype jokB, booleantype *jcurPtrB, double gammaB,
                 N_Vector tmp1B, N_Vector tmp2B, N_Vector tmp3B);
    /// <tt>M = I-gamma*df/dx</tt>, factorize
    void lsetup(CVodeMem cv_mem, int convfail, N_Vector ypred, N_Vector fpred, booleantype *jcurPtr,
                N_Vector vtemp1, N_Vector vtemp2, N_Vector vtemp3);
    void lsetupB(double t, double gamma, int convfail, N_Vector x, N_Vector xB, N_Vector xdotB,
                 booleantype *jcurPtr,
                 N_Vector vtemp1, N_Vector vtemp2, N_Vector vtemp3);
    /// <tt>b = M^(-1).b</tt>
    void lsolve(CVodeMem cv_mem, N_Vector b, N_Vector weight, N_Vector ycur, N_Vector fcur);
    void lsolveB(double t, double gamma, N_Vector b, N_Vector weight, N_Vector x,
                 N_Vector xB, N_Vector xdotB);

    // Static wrappers to be passed to Sundials
    static int rhs_wrapper(double t, N_Vector x, N_Vector xdot, void *user_data);
    static void ehfun_wrapper(int error_code, const char *module, const char *function, char *msg,
                              void *user_data);
    static int rhsS_wrapper(int Ns, double t, N_Vector x, N_Vector xdot, N_Vector *xF,
                            N_Vector *xdotF, void *user_data, N_Vector tmp1, N_Vector tmp2);
    static int rhsS1_wrapper(int Ns, double t, N_Vector x, N_Vector xdot, int iS,
                             N_Vector xF, N_Vector xdotF, void *user_data,
                             N_Vector tmp1, N_Vector tmp2);
    static int rhsQ_wrapper(double t, N_Vector x, N_Vector qdot, void *user_data);
    static int rhsQS_wrapper(int Ns, double t, N_Vector x, N_Vector *xF,
                             N_Vector qdot, N_Vector *qdotF,
                             void *user_data, N_Vector tmp1, N_Vector tmp2);
    static int rhsB_wrapper(double t, N_Vector x, N_Vector xB, N_Vector xdotB, void *user_data);
    static int rhsBS_wrapper(double t, N_Vector x, N_Vector *xF, N_Vector xB, N_Vector xdotB,
                             void *user_data);
    static int rhsQB_wrapper(double t, N_Vector x, N_Vector xB, N_Vector qdotB, void *user_data);
    static int jtimes_wrapper(N_Vector v, N_Vector Jv, double t, N_Vector x, N_Vector xdot,
                              void *user_data, N_Vector tmp);
    static int jtimesB_wrapper(N_Vector vB, N_Vector JvB, double t, N_Vector x, N_Vector xB,
                               N_Vector xdotB, void *user_data , N_Vector tmpB);
    static int djac_wrapper(long N, double t, N_Vector x, N_Vector xdot, DlsMat Jac,
                            void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
    static int djacB_wrapper(long NeqB, double t, N_Vector x, N_Vector xB, N_Vector xdotB,
                             DlsMat JacB, void *user_data,
                             N_Vector tmp1B, N_Vector tmp2B, N_Vector tmp3B);
    static int bjac_wrapper(long N, long mupper, long mlower, double t, N_Vector x, N_Vector xdot,
                            DlsMat Jac, void *user_data,
                            N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
    static int bjacB_wrapper(long NeqB, long mupperB, long mlowerB, double t,
                             N_Vector x, N_Vector xB,
                             N_Vector xdotB, DlsMat JacB, void *user_data,
                             N_Vector tmp1B, N_Vector tmp2B, N_Vector tmp3B);
    static int psolve_wrapper(double t, N_Vector x, N_Vector xdot, N_Vector r, N_Vector z,
                              double gamma, double delta, int lr, void *user_data, N_Vector tmp);
    static int psolveB_wrapper(double t, N_Vector x, N_Vector xB, N_Vector xdotB, N_Vector rvecB,
                               N_Vector zvecB, double gammaB, double deltaB,
                               int lr, void *user_data, N_Vector tmpB);
    static int psetup_wrapper(double t, N_Vector x, N_Vector xdot, booleantype jok,
                              booleantype *jcurPtr, double gamma, void *user_data,
                              N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
    static int psetupB_wrapper(double t, N_Vector x, N_Vector xB, N_Vector xdotB,
                               booleantype jokB, booleantype *jcurPtrB, double gammaB,
                               void *user_data, N_Vector tmp1B, N_Vector tmp2B, N_Vector tmp3B);
    static int lsetup_wrapper(CVodeMem cv_mem, int convfail, N_Vector x, N_Vector xdot,
                              booleantype *jcurPtr,
                              N_Vector vtemp1, N_Vector vtemp2, N_Vector vtemp3);
    static int lsolve_wrapper(CVodeMem cv_mem, N_Vector b, N_Vector weight, N_Vector x,
                              N_Vector xdot);
    static int lsetupB_wrapper(CVodeMem cv_mem, int convfail, N_Vector x, N_Vector xdot,
                               booleantype *jcurPtr,
                               N_Vector vtemp1, N_Vector vtemp2, N_Vector vtemp3);
    static int lsolveB_wrapper(CVodeMem cv_mem, N_Vector b, N_Vector weight,
                               N_Vector x, N_Vector xdot);

    // CVodes memory block
    void* mem_;

    // For timings
    clock_t time1, time2;

    // Accumulated time since last reset:
    double t_res; // time spent in the DAE residual
    double t_fres; // time spent in the forward sensitivity residual
    double t_jac; // time spent in the Jacobian, or Jacobian times vector function
    double t_lsolve; // preconditioner/linear solver solve function
    double t_lsetup_jac; // preconditioner/linear solver setup function, generate Jacobian
    double t_lsetup_fac; // preconditioner setup function, factorize Jacobian

    // N-vectors for the forward integration
    N_Vector x0_, x_, q_;

    // N-vectors for the backward integration
    N_Vector rx0_, rx_, rq_;

    // N-vectors for the forward sensitivities
    std::vector<N_Vector> xF0_, xF_, qF_;

    bool isInitAdj_;

    int ism_;

    // Calculate the error message map
    static std::map<int, std::string> calc_flagmap();

    // Error message map
    static std::map<int, std::string> flagmap;

    // Throw error
    static void cvodes_error(const std::string& module, int flag);

    // Ids of backward problem
    int whichB_;

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

    int lmm_; // linear multistep method
    int iter_; // nonlinear solver iteration

    bool monitor_rhsB_;
    bool monitor_rhs_;
    bool monitor_rhsQB_;

    bool disable_internal_warnings_;

  };

} // namespace casadi

/// \endcond
#endif // CASADI_CVODES_INTERFACE_HPP


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


#ifndef CASADI_IDAS_INTERFACE_HPP
#define CASADI_IDAS_INTERFACE_HPP

#include <casadi/interfaces/sundials/casadi_integrator_idas_export.h>
#include "sundials_interface.hpp"
#include "casadi/core/function/linear_solver.hpp"
#include <idas/idas.h>            /* prototypes for CVODE fcts. and consts. */
#include <idas/idas_dense.h>
#include <idas/idas_band.h>
#include <idas/idas_spgmr.h>
#include <idas/idas_spbcgs.h>
#include <idas/idas_sptfqmr.h>
#include <idas/idas_impl.h> /* Needed for the provided linear solver */
#include <ctime>

/** \defgroup plugin_Integrator_idas
      Interface to IDAS from the Sundials suite.

      Note: depending on the dimension and structure of your problem,
      you may experience a dramatic speed-up by using a sparse linear solver:

      \verbatim
       intg.setOption("linear_solver","csparse")
       intg.setOption("linear_solver_type","user_defined")
      \endverbatim
*/

/** \pluginsection{Integrator,idas} */

/// \cond INTERNAL

namespace casadi {

  /** \brief \pluginbrief{Integrator,idas}

      @copydoc IdasIntegrator_doc
      @copydoc plugin_Integrator_idas

      \author Joel Andersson
      \date 2010
  */
  class CASADI_INTEGRATOR_IDAS_EXPORT IdasInterface : public SundialsInterface {

  public:

    /** \brief  Constructor */
    explicit IdasInterface(const Function& f, const Function& g);

    /** \brief  Copy constructor */
    //  IdasInterface(const IdasInterface& integrator);

    /** \brief  Clone */
    virtual IdasInterface* clone() const;

    /** \brief  Create a new integrator */
    virtual IdasInterface* create(const Function& f, const Function& g) const
    { return new IdasInterface(f, g);}

    /** \brief  Create a new integrator */
    static IntegratorInternal* creator(const Function& f, const Function& g)
    { return new IdasInterface(f, g);}

    /** \brief  Deep copy data members */
    virtual void deepCopyMembers(std::map<SharedObjectNode*, SharedObject>& already_copied);

    /** \brief  Destructor */
    virtual ~IdasInterface();

    /** \brief  Free all IDAS memory */
    virtual void freeIDAS();

    /** \brief  Initialize */
    virtual void init();

    /** \brief Initialize the taping */
    virtual void initTaping();

    /** \brief Initialize the backward problem (can only be called after the first integration) */
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

    /** \brief  Get the integrator Jacobian for the backward problem (generic)
     *   Structure:
     *
     * \verbatim
     *   | diff(gx, rx) + cj*diff(gx, dot(rx))  |   diff(gx, rz) |
     *   | diff(gz, rx)                        |   diff(gz, rz) |
     * \endverbatim
     */
    template<typename FunctionType>
      FunctionType getJacGenB();

    /** \brief  Get the integrator Jacobian for the forward problem */
    virtual Function getJac();

    /** \brief  Get the integrator Jacobian for the backward problem */
    virtual Function getJacB();

    /// Correct the initial conditions, i.e. calculate
    void correctInitialConditions();

    /// A documentation string
    static const std::string meta_doc;

  protected:

    // Sundials callback functions
    void res(double t, const double* xz, const double* xzdot, double* rr);
    void ehfun(int error_code, const char *module, const char *function, char *msg);
    void jtimes(double t, const double *xz, const double *xzdot, const double *rr,
                const double *v, double *Jv, double cj, double *tmp1, double *tmp2);
    void jtimesB(double t, const double *xz, const double *xzdot, const double *xzB,
                 const double *xzdotB, const double *resvalB, const double *vB, double *JvB,
                 double cjB, double * tmp1B, double * tmp2B);
    void resS(int Ns, double t, const double* xz, const double* xzdot, const double *resval,
              N_Vector *xzF, N_Vector* xzdotF, N_Vector *rrF,
              double *tmp1, double *tmp2, double *tmp3);
    void rhsQ(double t, const double* xz, const double* xzdot, double* qdot);
    void rhsQS(int Ns, double t, N_Vector xz, N_Vector xzdot, N_Vector *xzF, N_Vector *xzdotF,
               N_Vector rrQ, N_Vector *qdotF,
               N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
    void resB(double t, const double* y, const double* xzdot, const double* xA,
              const double* xzdotB, double* rrB);
    void rhsQB(double t, const double* y, const double* xzdot, const double* xA,
               const double* xzdotB, double *qdotA);
    void psolve(double t, N_Vector xz, N_Vector xzdot, N_Vector rr, N_Vector rvec, N_Vector zvec,
                double cj, double delta, N_Vector tmp);
    void psolveB(double t, N_Vector xz, N_Vector xzdot, N_Vector xzB, N_Vector xzdotB,
                 N_Vector resvalB, N_Vector rvecB, N_Vector zvecB, double cjB,
                 double deltaB, N_Vector tmpB);
    void psetup(double t, N_Vector xz, N_Vector xzdot, N_Vector rr, double cj,
                N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
    void psetupB(double t, N_Vector xz, N_Vector xzdot, N_Vector xzB, N_Vector xzdotB,
                 N_Vector resvalB, double cjB,
                 N_Vector tmp1B, N_Vector tmp2B, N_Vector tmp3B);
    void djac(long Neq, double t, double cj, N_Vector xz, N_Vector xzdot, N_Vector rr, DlsMat Jac,
              N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
    void djacB(long NeqB, double t, double cjB, N_Vector xz, N_Vector xzdot, N_Vector xzB,
               N_Vector xzdotB, N_Vector rrB, DlsMat JacB,
               N_Vector tmp1B, N_Vector tmp2B, N_Vector tmp3B);
    void bjac(long Neq, long mupper, long mlower, double tt, double cj, N_Vector xz,
              N_Vector xzdot, N_Vector rr, DlsMat Jac,
              N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
    void bjacB(long NeqB, long mupperB, long mlowerB, double tt, double cjB, N_Vector xz,
               N_Vector xzdot, N_Vector xzB, N_Vector xzdotB, N_Vector resvalB, DlsMat JacB,
               N_Vector tmp1B, N_Vector tmp2B, N_Vector tmp3B);
    void lsetup(IDAMem IDA_mem, N_Vector xz, N_Vector xzdot, N_Vector resp,
                N_Vector vtemp1, N_Vector vtemp2, N_Vector vtemp3);
    void lsetupB(double t, double cj, N_Vector xz, N_Vector xzdot, N_Vector xzB, N_Vector xzdotB,
                 N_Vector resp, N_Vector vtemp1, N_Vector vtemp2, N_Vector vtemp3);
    void lsolve(IDAMem IDA_mem, N_Vector b, N_Vector weight, N_Vector xz, N_Vector xzdot,
                N_Vector rr);
    void lsolveB(double t, double cj, double cjratio, N_Vector b, N_Vector weight,
                 N_Vector xz, N_Vector xzdot, N_Vector xzB, N_Vector xzdotB, N_Vector rr);

    // Static wrappers to be passed to Sundials
    static int res_wrapper(double t, N_Vector xz, N_Vector xzdot, N_Vector rr, void *user_data);
    static int resB_wrapper(double t, N_Vector xz, N_Vector xzdot, N_Vector xzB, N_Vector xzdotB,
                            N_Vector rrB, void *user_data);
    static void ehfun_wrapper(int error_code, const char *module, const char *function, char *msg,
                              void *eh_data);
    static int jtimes_wrapper(double t, N_Vector xz, N_Vector xzdot, N_Vector rr, N_Vector v,
                              N_Vector Jv, double cj, void *user_data,
                              N_Vector tmp1, N_Vector tmp2);
    static int jtimesB_wrapper(double t, N_Vector xz, N_Vector xzdot, N_Vector xzB, N_Vector xzdotB,
                               N_Vector resvalB, N_Vector vB, N_Vector JvB, double cjB,
                               void *user_data, N_Vector tmp1B, N_Vector tmp2B);
    static int resS_wrapper(int Ns, double t, N_Vector xz, N_Vector xzdot, N_Vector resval,
                            N_Vector *xzF, N_Vector *xzdotF, N_Vector *resF, void *user_data,
                            N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
    static int rhsQ_wrapper(double t, N_Vector xz, N_Vector xzdot, N_Vector qdot, void *user_data);
    static int rhsQB_wrapper(double t, N_Vector xz, N_Vector xzdot, N_Vector xzB, N_Vector xzdotB,
                             N_Vector qdotA, void *user_data);
    static int rhsQS_wrapper(int Ns, double t, N_Vector xz, N_Vector xzdot,
                             N_Vector *xzF, N_Vector *xzdotF, N_Vector rrQ, N_Vector *qdotF,
                             void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
    static int psolve_wrapper(double t, N_Vector xz, N_Vector xzdot, N_Vector rr, N_Vector rvec,
                              N_Vector zvec, double cj, double delta, void *user_data,
                              N_Vector tmp);
    static int psetup_wrapper(double t, N_Vector xz, N_Vector xzdot, N_Vector rr, double cj,
                              void* user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
    static int psolveB_wrapper(double t,
                               N_Vector xz, N_Vector xzdot, N_Vector xzB, N_Vector xzdotB,
                               N_Vector resvalB, N_Vector rvecB, N_Vector zvecB, double cjB,
                               double deltaB, void *user_dataB, N_Vector tmpB);
    static int psetupB_wrapper(double t,
                               N_Vector xz, N_Vector xzdot, N_Vector xzB, N_Vector xzdotB,
                               N_Vector resvalB, double cjB, void *user_dataB,
                               N_Vector tmp1B, N_Vector tmp2B, N_Vector tmp3B);
    static int djac_wrapper(long Neq, double t, double cj, N_Vector xz, N_Vector xzdot,
                            N_Vector rr, DlsMat Jac, void *user_data,
                            N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
    static int djacB_wrapper(long NeqB, double t, double cjB, N_Vector xz, N_Vector xzdot,
                             N_Vector xzB, N_Vector xzdotB, N_Vector rrB, DlsMat JacB,
                             void *user_data, N_Vector tmp1B, N_Vector tmp2B, N_Vector tmp3B);
    static int bjac_wrapper(long Neq, long mupper, long mlower, double t, double cj, N_Vector xz,
                            N_Vector xzdot, N_Vector rr, DlsMat Jac, void *user_data,
                            N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
    static int bjacB_wrapper(long NeqB, long mupperB, long mlowerB, double t, double cjB,
                             N_Vector xz, N_Vector xzdot, N_Vector xzB, N_Vector xzdotB,
                             N_Vector resvalB, DlsMat JacB, void *user_data,
                             N_Vector tmp1B, N_Vector tmp2B, N_Vector tmp3B);
    static int lsetup_wrapper(IDAMem IDA_mem, N_Vector xz, N_Vector xzdot, N_Vector resp,
                              N_Vector vtemp1, N_Vector vtemp2, N_Vector vtemp3);
    static int lsolve_wrapper(IDAMem IDA_mem, N_Vector b, N_Vector weight, N_Vector ycur,
                              N_Vector xzdotcur, N_Vector rescur);
    static int lsetupB_wrapper(IDAMem IDA_mem, N_Vector xz, N_Vector xzdot, N_Vector resp,
                               N_Vector vtemp1, N_Vector vtemp2, N_Vector vtemp3);
    static int lsolveB_wrapper(IDAMem IDA_mem, N_Vector b, N_Vector weight, N_Vector ycur,
                               N_Vector xzdotcur, N_Vector rescur);
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

    // Accumulated time since last reset:
    double t_res; // time spent in the DAE residual
    double t_fres; // time spent in the forward sensitivity residual
    double t_jac, t_jacB; // time spent in the Jacobian, or Jacobian times vector function
    double t_lsolve; // preconditioner/linear solver solve function
    double t_lsetup_jac; // preconditioner/linear solver setup function, generate Jacobian
    double t_lsetup_fac; // preconditioner setup function, factorize Jacobian

    // Has the adjoint problem been initialized
    bool isInitAdj_;
    bool isInitTaping_;

    // Options
    bool cj_scaling_;
    bool calc_ic_;
    bool calc_icB_;

    // Disable IDAS internal warning messages
    bool disable_internal_warnings_;

    //  Initial values for \p xdot and \p z
    std::vector<double> init_xdot_;

  };

} // namespace casadi

/// \endcond
#endif // CASADI_IDAS_INTERFACE_HPP


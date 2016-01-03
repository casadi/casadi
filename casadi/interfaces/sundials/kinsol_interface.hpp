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


#ifndef CASADI_KINSOL_INTERFACE_HPP
#define CASADI_KINSOL_INTERFACE_HPP

#include <casadi/interfaces/sundials/casadi_rootfinder_kinsol_export.h>
#include "casadi/core/function/rootfinder.hpp"
#include <nvector/nvector_serial.h>   /* serial N_Vector types, fcts., and macros */
#include <sundials/sundials_dense.h>  /* definitions DlsMat DENSE_ELEM */
#include <sundials/sundials_types.h>  /* definition of type double */
#include <kinsol/kinsol.h>            /* prototypes for CVode fcts. and consts. */
#include <kinsol/kinsol_dense.h>
#include <kinsol/kinsol_band.h>
#include <kinsol/kinsol_spgmr.h>
#include <kinsol/kinsol_spbcgs.h>
#include <kinsol/kinsol_sptfqmr.h>
#include <kinsol/kinsol_impl.h> /* Needed for the provided linear solver */
#include <ctime>

/** \defgroup plugin_Rootfinder_kinsol
 KINSOL interface from the Sundials suite
*/
/** \pluginsection{Rootfinder,kinsol} */

/// \cond INTERNAL
namespace casadi {

  /** \brief \pluginbrief{Rootfinder,kinsol}
  *
  * @copydoc Rootfinder_doc
  * @copydoc plugin_Rootfinder_kinsol
  */
  class CASADI_ROOTFINDER_KINSOL_EXPORT KinsolInterface : public Rootfinder {
  public:
    /** \brief  Constructor */
    explicit KinsolInterface(const std::string& name, const Function& f);

    /** \brief  Destructor */
    virtual ~KinsolInterface();

    /** \brief  Create a new Rootfinder */
    static Rootfinder* creator(const std::string& name, const Function& f) {
      return new KinsolInterface(name, f);
    }

    /** \brief  Initialize stage */
    virtual void init();

    /// Solve the system of equations and calculate derivatives
    virtual void eval(const double** arg, double** res, int* iw, double* w, void* mem);

    // Get name of the plugin
    virtual const char* plugin_name() const { return "kinsol";}

    // Scaling
    N_Vector u_scale_, f_scale_;

    /// Globalization strategy
    int strategy_;

    // Should KINSOL internal warning messages be ignored
    bool disable_internal_warnings_;

    // Maximum number of iterations
    int max_iter_;

    // Use exact Jacobian?
    bool exact_jacobian_;

    // Type of linear solver
    enum LinsolType { DENSE, BANDED, ITERATIVE, USER_DEFINED};
    LinsolType linear_solver_type_;

    // Bandwidth (for banded solvers)
    int upper_bandwidth_, lower_bandwidth_;

    // Krylov subspace size (for iterative solvers)
    int maxl_;

    // Iterative solver
    enum IterativeSolver { GMRES, BCGSTAB, TFQMR};
    IterativeSolver iterative_solver_;

    // Should a preconditioner be used
    bool use_preconditioner_;

    // Absolute tolerance
    double abstol_;

    // Jacobian times vector function
    Function f_fwd_;

    // Raise an error specific to KinSol
    void kinsol_error(const std::string& module, int flag, bool fatal=true);

    /// A documentation string
    static const std::string meta_doc;

    /** \brief Allocate memory block */
    virtual Memory2* alloc_mem();

    /** \brief Free allocated memory block */
    virtual void free_mem(void* mem);
  };

  // Memory
  struct CASADI_ROOTFINDER_KINSOL_EXPORT KinsolMemory : public Memory2 {
    /// Shared memory
    KinsolInterface& self;

    /// Constructor
    KinsolMemory(KinsolInterface& s);

    /// Destructor
    virtual ~KinsolMemory();

    /// KINSOL memory block
    void* mem_;

    // Linear solver memory
    Memory linsol_mem_;

    /// Function arguments
    const double** arg_;
    double** res_;
    int* iw_;
    double* w_;

    /// Variable
    N_Vector u_;

    /// For timings
    clock_t time1_, time2_;

    /// Accumulated time since last reset:
    double t_func_; // time spent in the residual function
    double t_jac_; // time spent in the Jacobian function
    double t_lsolve_; // preconditioner/linear solver solve function
    double t_lsetup_jac_; // preconditioner/linear solver setup function, generate Jacobian
    double t_lsetup_fac_; // preconditioner setup function, factorize Jacobian

    /** \brief Callback functions */
    void func(N_Vector u, N_Vector fval);
    void djac(long N, N_Vector u, N_Vector fu, DlsMat J, N_Vector tmp1, N_Vector tmp2);
    void bjac(long N, long mupper, long mlower, N_Vector u, N_Vector fu, DlsMat J,
              N_Vector tmp1, N_Vector tmp2);
    void jtimes(N_Vector v, N_Vector Jv, N_Vector u, int* new_u);
    void psetup(N_Vector u, N_Vector uscale, N_Vector fval, N_Vector fscale,
                N_Vector tmp1, N_Vector tmp2);
    void psolve(N_Vector u, N_Vector uscale, N_Vector fval, N_Vector fscale, N_Vector v,
                N_Vector tmp);
    void lsetup(KINMem kin_mem);
    void lsolve(KINMem kin_mem, N_Vector x, N_Vector b, double *res_norm);
    void ehfun(int error_code, const char *module, const char *function, char *msg);

    /** \brief Wrappers to callback functions*/
    static int func_wrapper(N_Vector u, N_Vector fval, void *user_data);
    static int djac_wrapper(long N, N_Vector u, N_Vector fu, DlsMat J, void *user_data,
                            N_Vector tmp1, N_Vector tmp2);
    static int bjac_wrapper(long N, long mupper, long mlower, N_Vector u, N_Vector fu, DlsMat J,
                            void *user_data, N_Vector tmp1, N_Vector tmp2);
    static int jtimes_wrapper(N_Vector v, N_Vector Jv, N_Vector u, int* new_u, void *user_data);
    static int psetup_wrapper(N_Vector u, N_Vector uscale, N_Vector fval, N_Vector fscale,
                              void* user_data, N_Vector tmp1, N_Vector tmp2);
    static int psolve_wrapper(N_Vector u, N_Vector uscale, N_Vector fval, N_Vector fscale,
                              N_Vector v, void* user_data, N_Vector tmp);
    static int lsetup_wrapper(KINMem kin_mem);
    static int lsolve_wrapper(KINMem kin_mem, N_Vector x, N_Vector b, double *res_norm);
    static void ehfun_wrapper(int error_code, const char *module, const char *function,
                              char *msg, void *eh_data);
  };

} // namespace casadi

/// \endcond
#endif // CASADI_KINSOL_INTERFACE_HPP


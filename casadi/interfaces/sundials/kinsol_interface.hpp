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

#include <casadi/interfaces/sundials/casadi_implicitfunction_kinsol_export.h>
#include "casadi/core/function/implicit_function_internal.hpp"
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

/** \defgroup plugin_ImplicitFunction_kinsol
 KINSOL interface from the Sundials suite
*/
/** \pluginsection{ImplicitFunction,kinsol} */

/// \cond INTERNAL
namespace casadi {

  typedef std::pair< std::string, std::string> Message;

  /** \brief \pluginbrief{ImplicitFunction,kinsol}
  *
  * @copydoc ImplicitFunction_doc
  * @copydoc plugin_ImplicitFunction_kinsol
  */
  class CASADI_IMPLICITFUNCTION_KINSOL_EXPORT KinsolInterface : public ImplicitFunctionInternal {
  public:
    /** \brief  Constructor */
    explicit KinsolInterface(const Function& f);

    /** \brief  Destructor */
    virtual ~KinsolInterface();

    /** \brief  Clone */
    virtual KinsolInterface* clone() const;

    /** \brief  Create a new ImplicitFunctionInternal */
    virtual ImplicitFunctionInternal* create(const Function& f) const
    { return new KinsolInterface(f);}

    /** \brief  Create a new ImplicitFunction */
    static ImplicitFunctionInternal* creator(const Function& f)
    { return new KinsolInterface(f);}

    /** \brief  Initialize stage */
    virtual void init();

    /** \brief  Solve the nonlinear system of equations */
    virtual void solveNonLinear();

    /** \brief Residual */
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

    /** \brief Wrappers */
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

    /// KINSOL memory block
    void* mem_;

    /// Variable
    N_Vector u_;

    // Scaling
    N_Vector u_scale_, f_scale_;

    /// For timings
    clock_t time1_, time2_;

    /// Accumulated time since last reset:
    double t_func_; // time spent in the residual function
    double t_jac_; // time spent in the Jacobian function
    double t_lsolve_; // preconditioner/linear solver solve function
    double t_lsetup_jac_; // preconditioner/linear solver setup function, generate Jacobian
    double t_lsetup_fac_; // preconditioner setup function, factorize Jacobian

    /// Globalization strategy
    int strategy_;

    // Should KINSOL internal warning messages be ignored
    bool disable_internal_warnings_;

    // Jacobian times vector function
    Function f_fwd_;

    // Calculate the error message map
    static std::map<int, Message > calc_flagmap();

    // Error message map
    static std::map<int, Message> flagmap;


    // Raise an error specific to KinSol
    void kinsol_error(const std::string& module, int flag, bool fatal=true);

    /// A documentation string
    static const std::string meta_doc;

  };

} // namespace casadi

/// \endcond
#endif // CASADI_KINSOL_INTERFACE_HPP


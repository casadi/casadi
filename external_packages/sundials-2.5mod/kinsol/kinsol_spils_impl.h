/*
 * -----------------------------------------------------------------
 * $Revision: 1.5 $
 * $Date: 2010/12/01 22:43:33 $
 * -----------------------------------------------------------------
 * Programmer(s): Radu Serban and Aaron Collier @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2002, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see the LICENSE file.
 * -----------------------------------------------------------------
 * Common implementation header file for the scaled, preconditioned
 * linear solver modules.
 * -----------------------------------------------------------------
 */

#ifndef _KINSPILS_IMPL_H
#define _KINSPILS_IMPL_H

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

#include <kinsol/kinsol_spils.h>
#include "kinsol_impl.h"

/* Types of iterative linear solvers */

#define SPILS_SPGMR   1
#define SPILS_SPBCG   2
#define SPILS_SPTFQMR 3

/*
 * -----------------------------------------------------------------
 * keys for KINPrintInfo (do not use 1 -> conflict with PRNT_RETVAL)
 * -----------------------------------------------------------------
 */

#define PRNT_NLI   101
#define PRNT_EPS   102

/*
 * -----------------------------------------------------------------
 * Types : struct KINSpilsMemRec and struct *KINSpilsMem
 * -----------------------------------------------------------------
 * A variable declaration of type struct *KINSpilsMem denotes a
 * pointer to a data structure of type struct KINSpilsMemRec. The
 * KINSpilsMemRec structure contains fields that must be accessible
 * by KINSPILS/SPGMR solver module routines.
 * -----------------------------------------------------------------
 */

typedef struct KINSpilsMemRec {

  int s_type;           /* type of scaled preconditioned iterative LS          */

  /* problem specification data */

  int  s_maxl;          /* maximum allowable dimension of Krylov subspace      */     
  int  s_pretype;       /* preconditioning type: PREC_NONE, PREC_RIGHT,
			   PREC_LEFT or PREC_BOTH (used by SPGMR module and
			   defined in sundials_iterative.h)                    */
  int  s_gstype;        /* Gram-Schmidt orthogonalization procedure:
			   CLASSICAL_GS or MODIFIED_GS (used by SPGMR module
			   and defined in sundials_iterative.h)                */
  booleantype s_new_uu; /* flag indicating if the iterate has been updated -
			   Jacobian must be updated/reevaluated (meant to be
			   used by user-supplied jtimes function)              */
  int s_maxlrst;        /* maximum number of times the SPGMR linear solver
			   can be restarted                                    */

  /* counters */

  long int s_nli;     /* number of linear iterations performed                 */
  long int s_npe;     /* number of preconditioner evaluations                  */
  long int s_nps;     /* number of calls to preconditioner solve fun.          */
  long int s_ncfl;    /* number of linear convergence failures                 */
  long int s_nfes;    /* number of evaluations of the system function F(u) or
			 number of calls made to func routine                  */    
  long int s_njtimes; /* number of times the matrix-vector product J(u)*v
			 was computed or number of calls made to jtimes
			 routine                                               */

  /* miscellaneous */

  void *s_spils_mem;    /* pointer to generic linear solver memory block       */

  long int s_last_flag; /* last flag returned                                  */


   /* Preconditioner computation
   * (a) user-provided:
   *     - P_data == user_data
   *     - pfree == NULL (the user dealocates memory for user_data)
   * (b) internal preconditioner module
   *     - P_data == kinsol_mem
   *     - pfree == set by the prec. module and called in KINSpilsFree
   */
 
  KINSpilsPrecSetupFn s_pset;     
  KINSpilsPrecSolveFn s_psolve; 
  void (*s_pfree)(KINMem kin_mem);
  void *s_P_data;

  /* Jacobian times vector compuation
   * (a) jtimes function provided by the user:
   *     - J_data == user_data
   *     - jtimesDQ == FALSE
   * (b) internal jtimes
   *     - J_data == kinsol_mem
   *     - jtimesDQ == TRUE
   */

  booleantype s_jtimesDQ;
  KINSpilsJacTimesVecFn s_jtimes;
  void *s_J_data;

} *KINSpilsMem;


/*
 * -----------------------------------------------------------------
 * Prototypes of internal functions
 * -----------------------------------------------------------------
 */

/* KINSpgmr Atimes and PSolve routines called by generic SPGMR solver */

int KINSpilsAtimes(void *kinsol_mem, N_Vector v, N_Vector z);
int KINSpilsPSolve(void *kinsol_mem, N_Vector r, N_Vector z, int lr);

/* difference quotient approximation for jacobian times vector */

int KINSpilsDQJtimes(N_Vector v, N_Vector Jv,
                     N_Vector u, booleantype *new_u, 
                     void *data);

/*
 * -----------------------------------------------------------------
 * KINSPILS error messages
 * -----------------------------------------------------------------
 */

#define MSGS_KINMEM_NULL "KINSOL memory is NULL."
#define MSGS_MEM_FAIL    "A memory request failed."
#define MSGS_BAD_NVECTOR "A required vector operation is not implemented."
#define MSGS_LMEM_NULL   "Linear solver memory is NULL."
#define MSGS_NEG_MAXRS   "maxrs < 0 illegal."


/*
 * -----------------------------------------------------------------
 * KINSPILS info messages
 * -----------------------------------------------------------------
 */

#define INFO_NLI  "nli_inc = %d"

#if defined(SUNDIALS_EXTENDED_PRECISION)

#define INFO_EPS  "residual norm = %12.3Lg  eps = %12.3Lg"

#elif defined(SUNDIALS_DOUBLE_PRECISION)

#define INFO_EPS  "residual norm = %12.3lg  eps = %12.3lg"

#else

#define INFO_EPS  "residual norm = %12.3g  eps = %12.3g"

#endif



#ifdef __cplusplus
}
#endif

#endif

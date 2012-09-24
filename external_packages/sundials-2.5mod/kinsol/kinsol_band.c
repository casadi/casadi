/*
 * -----------------------------------------------------------------
 * $Revision: 1.11 $
 * $Date: 2011/03/23 23:37:59 $
 * ----------------------------------------------------------------- 
 * Programmer(s): Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2002, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see the LICENSE file.
 * -----------------------------------------------------------------
 * This is the implementation file for the KINBAND linear solver.
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>

#include <kinsol/kinsol_band.h>
#include "kinsol_direct_impl.h"
#include "kinsol_impl.h"

#include <sundials/sundials_math.h>

/* Constants */

#define ZERO         RCONST(0.0)
#define ONE          RCONST(1.0)
#define TWO          RCONST(2.0)

/* 
 * =================================================================
 * PROTOTYPES FOR PRIVATE FUNCTIONS
 * =================================================================
 */

/* KINBAND linit, lsetup, lsolve, and lfree routines */
static int kinBandInit(KINMem kin_mem);
static int kinBandSetup(KINMem kin_mem);
static int kinBandsolve(KINMem kin_mem, N_Vector x, N_Vector b,
                        realtype *res_norm);
static void kinBandFree(KINMem kin_mem);

/*
 * =================================================================
 * READIBILITY REPLACEMENTS
 * =================================================================
 */

#define lrw1           (kin_mem->kin_lrw1)
#define liw1           (kin_mem->kin_liw1)
#define func           (kin_mem->kin_func)
#define printfl        (kin_mem->kin_printfl)
#define linit          (kin_mem->kin_linit)
#define lsetup         (kin_mem->kin_lsetup)
#define lsolve         (kin_mem->kin_lsolve)
#define lfree          (kin_mem->kin_lfree)
#define lmem           (kin_mem->kin_lmem)
#define inexact_ls     (kin_mem->kin_inexact_ls)
#define uu             (kin_mem->kin_uu)
#define fval           (kin_mem->kin_fval)
#define uscale         (kin_mem->kin_uscale)
#define fscale         (kin_mem->kin_fscale)
#define sqrt_relfunc   (kin_mem->kin_sqrt_relfunc)
#define sJpnorm        (kin_mem->kin_sJpnorm)
#define sfdotJp        (kin_mem->kin_sfdotJp)
#define errfp          (kin_mem->kin_errfp)
#define infofp         (kin_mem->kin_infofp)
#define setupNonNull   (kin_mem->kin_setupNonNull)
#define vtemp1         (kin_mem->kin_vtemp1)
#define vec_tmpl       (kin_mem->kin_vtemp1)
#define vtemp2         (kin_mem->kin_vtemp2)

#define mtype          (kindls_mem->d_type)
#define n              (kindls_mem->d_n)
#define ml             (kindls_mem->d_ml)
#define mu             (kindls_mem->d_mu)
#define smu            (kindls_mem->d_smu)
#define jacDQ          (kindls_mem->d_jacDQ)
#define bjac           (kindls_mem->d_bjac)
#define J              (kindls_mem->d_J)
#define lpivots        (kindls_mem->d_lpivots)
#define nje            (kindls_mem->d_nje)
#define nfeDQ          (kindls_mem->d_nfeDQ)
#define J_data         (kindls_mem->d_J_data)
#define last_flag      (kindls_mem->d_last_flag)

/* 
 * =================================================================
 * EXPORTED FUNCTIONS
 * =================================================================
 */

/*
 * -----------------------------------------------------------------
 * KINBand
 * -----------------------------------------------------------------
 * This routine initializes the memory record and sets various function
 * fields specific to the band linear solver module.  KINBand first calls
 * the existing lfree routine if this is not NULL.  It then sets the
 * cv_linit, cv_lsetup, cv_lsolve, and cv_lfree fields in (*cvode_mem)
 * to be kinBandInit, kinBandSetup, kinBandsolve, and kinBandFree,
 * respectively.  It allocates memory for a structure of type
 * KINDlsMemRec and sets the cv_lmem field in (*cvode_mem) to the
 * address of this structure.  It sets setupNonNull in (*cvode_mem) to be
 * TRUE, b_mu to be mupper, b_ml to be mlower, and the bjac field to be 
 * kinDlsBandDQJac.
 * Finally, it allocates memory for M, savedJ, and pivot.
 *
 * NOTE: The band linear solver assumes a serial implementation
 *       of the NVECTOR package. Therefore, KINBand will first 
 *       test for compatible a compatible N_Vector internal
 *       representation by checking that the function 
 *       N_VGetArrayPointer exists.
 * -----------------------------------------------------------------
 */
                  
int KINBand(void *kinmem, long int N, long int mupper, long int mlower)
{
  KINMem kin_mem;
  KINDlsMem kindls_mem;

  /* Return immediately if kinmem is NULL */
  if (kinmem == NULL) {
    KINProcessError(NULL, KINDLS_MEM_NULL, "KINBAND", "KINBand", MSGD_KINMEM_NULL);
    return(KINDLS_MEM_NULL);
  }
  kin_mem = (KINMem) kinmem;

  /* Test if the NVECTOR package is compatible with the BAND solver */
  if (vec_tmpl->ops->nvgetarraypointer == NULL) {
    KINProcessError(kin_mem, KINDLS_ILL_INPUT, "KINBAND", "KINBand", MSGD_BAD_NVECTOR);
    return(KINDLS_ILL_INPUT);
  }

  if (lfree != NULL) lfree(kin_mem);

  /* Set four main function fields in kin_mem */  
  linit  = kinBandInit;
  lsetup = kinBandSetup;
  lsolve = kinBandsolve;
  lfree  = kinBandFree;
  
  /* Get memory for KINDlsMemRec */
  kindls_mem = NULL;
  kindls_mem = (KINDlsMem) malloc(sizeof(struct KINDlsMemRec));
  if (kindls_mem == NULL) {
    KINProcessError(kin_mem, KINDLS_MEM_FAIL, "KINBAND", "KINBand", MSGD_MEM_FAIL);
    return(KINDLS_MEM_FAIL);
  }

  /* Set matrix type */
  mtype = SUNDIALS_BAND;  

  /* Set default Jacobian routine and Jacobian data */
  jacDQ  = TRUE;
  bjac   = NULL;
  J_data = NULL;
  last_flag = KINDLS_SUCCESS;

  setupNonNull = TRUE;
  
  /* Load problem dimension */
  n = N;

  /* Load half-bandwiths in kindls_mem */
  ml = mlower;
  mu = mupper;

  /* Test ml and mu for legality */
  if ((ml < 0) || (mu < 0) || (ml >= N) || (mu >= N)) {
    KINProcessError(kin_mem, KINDLS_MEM_FAIL, "KINBAND", "KINBand", MSGD_MEM_FAIL);
    free(kindls_mem); kindls_mem = NULL;
    return(KINDLS_ILL_INPUT);
  }

  /* Set extended upper half-bandwith for M (required for pivoting) */
  smu = MIN(N-1, mu + ml);

  /* Allocate memory for J and pivot array */
  J = NULL;
  J = NewBandMat(N, mu, ml, smu);
  if (J == NULL) {
    KINProcessError(kin_mem, KINDLS_MEM_FAIL, "KINBAND", "KINBand", MSGD_MEM_FAIL);
    free(kindls_mem); kindls_mem = NULL;
    return(KINDLS_MEM_FAIL);
  }

  lpivots = NULL;
  lpivots = NewLintArray(N);
  if (lpivots == NULL) {
    KINProcessError(kin_mem, KINDLS_MEM_FAIL, "KINBAND", "KINBand", MSGD_MEM_FAIL);
    DestroyMat(J);
    free(kindls_mem); kindls_mem = NULL;
    return(KINDLS_MEM_FAIL);
  }

  /* This is a direct linear solver */
  inexact_ls = FALSE;

  /* Attach linear solver memory to integrator memory */
  lmem = kindls_mem;

  return(KINDLS_SUCCESS);
}

/* 
 * =================================================================
 *  PRIVATE FUNCTIONS
 * =================================================================
 */

/*
 * -----------------------------------------------------------------
 * kinBandInit
 * -----------------------------------------------------------------
 * This routine does remaining initializations specific to the band
 * linear solver.
 * -----------------------------------------------------------------
 */

static int kinBandInit(KINMem kin_mem)
{
  KINDlsMem kindls_mem;

  kindls_mem = (KINDlsMem) lmem;

  nje   = 0;
  nfeDQ = 0;

  if (jacDQ) {
    bjac = kinDlsBandDQJac;
    J_data = kin_mem;
  } else {
    J_data = kin_mem->kin_user_data;
  }

  last_flag = KINDLS_SUCCESS;
  return(0);
}

/*
 * -----------------------------------------------------------------
 * kinBandSetup
 * -----------------------------------------------------------------
 * This routine does the setup operations for the band linear solver.
 * It makes a decision whether or not to call the Jacobian evaluation
 * routine based on various state variables, and if not it uses the 
 * saved copy.  In any case, it constructs the Newton matrix 
 * M = I - gamma*J, updates counters, and calls the band LU 
 * factorization routine.
 * -----------------------------------------------------------------
 */

static int kinBandSetup(KINMem kin_mem)
{
  KINDlsMem kindls_mem;
  int retval;
  long int ier;

  kindls_mem = (KINDlsMem) lmem;

  nje++;
  SetToZero(J); 
  retval = bjac(n, mu, ml, uu, fval, J, J_data, vtemp1, vtemp2);
  if (retval != 0) {
    last_flag = -1;
    return(-1);
  }
  
  /* Do LU factorization of J */
  ier = BandGBTRF(J, lpivots);

  /* Return 0 if the LU was complete; otherwise return -1 */
  last_flag = ier;
  if (ier > 0) return(-1);

  return(0);
}

/*
 * -----------------------------------------------------------------
 * kinBandsolve
 * -----------------------------------------------------------------
 * This routine handles the solve operation for the band linear solver
 * by calling the band backsolve routine.  The return value is 0.
 * -----------------------------------------------------------------
 */

static int kinBandsolve(KINMem kin_mem, N_Vector x, N_Vector b, realtype *res_norm)
{
  KINDlsMem kindls_mem;
  realtype *xd;

  kindls_mem = (KINDlsMem) lmem;

  /* Copy the right-hand side into x */

  N_VScale(ONE, b, x);
  
  xd = N_VGetArrayPointer(x);

  /* Back-solve and get solution in x */

  BandGBTRS(J, lpivots, xd);

  /* Compute the terms Jpnorm and sfdotJp for use in the global strategy
     routines and in KINForcingTerm. Both of these terms are subsequently
     corrected if the step is reduced by constraints or the line search.
     
     sJpnorm is the norm of the scaled product (scaled by fscale) of
     the current Jacobian matrix J and the step vector p.

     sfdotJp is the dot product of the scaled f vector and the scaled
     vector J*p, where the scaling uses fscale. */

  sJpnorm = N_VWL2Norm(b,fscale);
  N_VProd(b, fscale, b);
  N_VProd(b, fscale, b);
  sfdotJp = N_VDotProd(fval, b);

  last_flag = KINDLS_SUCCESS;

  return(0);
}

/*
 * -----------------------------------------------------------------
 * kinBandFree
 * -----------------------------------------------------------------
 * This routine frees memory specific to the band linear solver.
 * -----------------------------------------------------------------
 */

static void kinBandFree(KINMem kin_mem)
{
  KINDlsMem kindls_mem;

  kindls_mem = (KINDlsMem) lmem;

  DestroyMat(J);
  DestroyArray(lpivots);
  free(kindls_mem); kindls_mem = NULL;
}


/*
 * -----------------------------------------------------------------
 * $Revision: 1.11 $
 * $Date: 2011/02/16 22:43:28 $
 * ----------------------------------------------------------------- 
 * Programmer: Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2006, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see the LICENSE file.
 * -----------------------------------------------------------------
 * This is the implementation file for a KINSOL dense linear solver
 * using BLAS and LAPACK functions.
 * -----------------------------------------------------------------
 */

/* 
 * =================================================================
 * IMPORTED HEADER FILES
 * =================================================================
 */

#include <stdio.h>
#include <stdlib.h>

#include <kinsol/kinsol_lapack.h>
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

/* KINLAPACK DENSE linit, lsetup, lsolve, and lfree routines */ 
static int kinLapackDenseInit(KINMem kin_mem);
static int kinLapackDenseSetup(KINMem kin_mem);
static int kinLapackDenseSolve(KINMem kin_mem, N_Vector x, N_Vector b, realtype *res_norm);
static void kinLapackDenseFree(KINMem kin_mem);

/* KINLAPACK BAND linit, lsetup, lsolve, and lfree routines */ 
static int kinLapackBandInit(KINMem kin_mem);
static int kinLapackBandSetup(KINMem kin_mem);
static int kinLapackBandSolve(KINMem kin_mem, N_Vector x, N_Vector b, realtype *res_norm);
static void kinLapackBandFree(KINMem kin_mem);

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
#define djac           (kindls_mem->d_djac)
#define bjac           (kindls_mem->d_bjac)
#define J              (kindls_mem->d_J)
#define pivots         (kindls_mem->d_pivots)
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
 * KINLapackDense
 * -----------------------------------------------------------------
 * This routine initializes the memory record and sets various function
 * fields specific to the linear solver module.  KINLapackDense first
 * calls the existing lfree routine if this is not NULL.  Then it sets
 * the kin_linit, kin_lsetup, kin_lsolve, kin_lfree fields in (*kinmem)
 * to be kinLapackDenseInit, kinLapackDenseSetup, kinLapackDenseSolve, 
 * and kinLapackDenseFree, respectively.  It allocates memory for a 
 * structure of type KINDlsMemRec and sets the kin_lmem field in 
 * (*kinmem) to the address of this structure.  It sets lsetup_exists 
 * in (*kinmem) to TRUE, and the djac field to the default 
 * kinLapackDenseDQJac. Finally, it allocates memory for M, pivots, and 
 * (if needed) savedJ.
 *
 * NOTE: The dense linear solver assumes a serial implementation
 *       of the NVECTOR package. Therefore, KINLapackDense will first 
 *       test for a compatible N_Vector internal representation 
 *       by checking that N_VGetArrayPointer and N_VSetArrayPointer 
 *       exist.
 * -----------------------------------------------------------------
 */
int KINLapackDense(void *kinmem, int N)
{
  KINMem kin_mem;
  KINDlsMem kindls_mem;

  /* Return immediately if kinmem is NULL */
  if (kinmem == NULL) {
    KINProcessError(NULL, KINDLS_MEM_NULL, "KINLAPACK", "KINLapackDense", MSGD_KINMEM_NULL);
    return(KINDLS_MEM_NULL);
  }
  kin_mem = (KINMem) kinmem;

  /* Test if the NVECTOR package is compatible with the DENSE solver */
  if (vec_tmpl->ops->nvgetarraypointer == NULL ||
      vec_tmpl->ops->nvsetarraypointer == NULL) {
    KINProcessError(kin_mem, KINDLS_ILL_INPUT, "KINLAPACK", "KINLapackDense", MSGD_BAD_NVECTOR);
    return(KINDLS_ILL_INPUT);
  }

  if (lfree !=NULL) lfree(kin_mem);

  /* Set four main function fields in kin_mem */
  linit  = kinLapackDenseInit;
  lsetup = kinLapackDenseSetup;
  lsolve = kinLapackDenseSolve;
  lfree  = kinLapackDenseFree;

  /* Get memory for KINDlsMemRec */
  kindls_mem = NULL;
  kindls_mem = (KINDlsMem) malloc(sizeof(struct KINDlsMemRec));
  if (kindls_mem == NULL) {
    KINProcessError(kin_mem, KINDLS_MEM_FAIL, "KINLAPACK", "KINLapackDense", MSGD_MEM_FAIL);
    return(KINDLS_MEM_FAIL);
  }

  /* Set matrix type */
  mtype = SUNDIALS_DENSE;  

  /* Set default Jacobian routine and Jacobian data */
  jacDQ  = TRUE;
  djac   = NULL;
  J_data = NULL;

  last_flag = KINDLS_SUCCESS;

  setupNonNull = TRUE;

  /* Set problem dimension */
  n = (long int) N;

  /* Allocate memory for J and pivot array */
  
  J = NULL;
  J = NewDenseMat(n, n);
  if (J == NULL) {
    KINProcessError(kin_mem, KINDLS_MEM_FAIL, "KINLAPACK", "KINLapackDense", MSGD_MEM_FAIL);
    free(kindls_mem); kindls_mem = NULL;
    return(KINDLS_MEM_FAIL);
  }

  pivots = NULL;
  pivots = NewIntArray(N);
  if (pivots == NULL) {
    KINProcessError(kin_mem, KINDLS_MEM_FAIL, "KINLAPACK", "KINLapackDense", MSGD_MEM_FAIL);
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
 * -----------------------------------------------------------------
 * KINLapackBand
 * -----------------------------------------------------------------
 * This routine initializes the memory record and sets various function
 * fields specific to the band linear solver module. It first calls
 * the existing lfree routine if this is not NULL.  It then sets the
 * kin_linit, kin_lsetup, kin_lsolve, and kin_lfree fields in (*kinmem)
 * to be kinLapackBandInit, kinLapackBandSetup, kinLapackBandSolve, 
 * and kinLapackBandFree, respectively.  It allocates memory for a 
 * structure of type KINLapackBandMemRec and sets the kin_lmem field in 
 * (*kinmem) to the address of this structure.  It sets lsetup_exists 
 * in (*kinmem) to be TRUE, mu to be mupper, ml to be mlower, and 
 * the bjac field to kinDlsBandDQJac
 * Finally, it allocates memory for M, pivots, and (if needed) savedJ.  
 * The KINLapackBand return value is KINDLS_SUCCESS = 0, 
 * KINDLS_MEM_FAIL = -1, or KINDLS_ILL_INPUT = -2.
 *
 * NOTE: The KINLAPACK linear solver assumes a serial implementation
 *       of the NVECTOR package. Therefore, KINLapackBand will first 
 *       test for compatible a compatible N_Vector internal
 *       representation by checking that the function 
 *       N_VGetArrayPointer exists.
 * -----------------------------------------------------------------
 */                  
int KINLapackBand(void *kinmem, int N, int mupper, int mlower)
{
  KINMem kin_mem;
  KINDlsMem kindls_mem;

  /* Return immediately if kinmem is NULL */
  if (kinmem == NULL) {
    KINProcessError(NULL, KINDLS_MEM_NULL, "KINLAPACK", "KINLapackBand", MSGD_KINMEM_NULL);
    return(KINDLS_MEM_NULL);
  }
  kin_mem = (KINMem) kinmem;

  /* Test if the NVECTOR package is compatible with the BAND solver */
  if (vec_tmpl->ops->nvgetarraypointer == NULL) {
    KINProcessError(kin_mem, KINDLS_ILL_INPUT, "KINLAPACK", "KINLapackBand", MSGD_BAD_NVECTOR);
    return(KINDLS_ILL_INPUT);
  }

  if (lfree != NULL) lfree(kin_mem);

  /* Set four main function fields in kin_mem */  
  linit  = kinLapackBandInit;
  lsetup = kinLapackBandSetup;
  lsolve = kinLapackBandSolve;
  lfree  = kinLapackBandFree;
  
  /* Get memory for KINDlsMemRec */
  kindls_mem = NULL;
  kindls_mem = (KINDlsMem) malloc(sizeof(struct KINDlsMemRec));
  if (kindls_mem == NULL) {
    KINProcessError(kin_mem, KINDLS_MEM_FAIL, "KINLAPACK", "KINLapackBand", MSGD_MEM_FAIL);
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
  n = (long int) N;

  /* Load half-bandwidths in kindls_mem */
  ml = (long int) mlower;
  mu = (long int) mupper;

  /* Test ml and mu for legality */
  if ((ml < 0) || (mu < 0) || (ml >= n) || (mu >= n)) {
    KINProcessError(kin_mem, KINDLS_MEM_FAIL, "KINLAPACK", "KINLapackBand", MSGD_MEM_FAIL);
    free(kindls_mem); kindls_mem = NULL;
    return(KINDLS_ILL_INPUT);
  }

  /* Set extended upper half-bandwith for M (required for pivoting) */
  smu = MIN(n-1, mu + ml);

  /* Allocate memory for J and pivot array */
  J = NULL;
  J = NewBandMat(n, mu, ml, smu);
  if (J == NULL) {
    KINProcessError(kin_mem, KINDLS_MEM_FAIL, "KINLAPACK", "KINLapackBand", MSGD_MEM_FAIL);
    free(kindls_mem); kindls_mem = NULL;
    return(KINDLS_MEM_FAIL);
  }

  pivots = NULL;
  pivots = NewIntArray(N);
  if (pivots == NULL) {
    KINProcessError(kin_mem, KINDLS_MEM_FAIL, "KINLAPACK", "KINLapackBand", MSGD_MEM_FAIL);
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
 *  PRIVATE FUNCTIONS FOR SOLUTION WITH DENSE JACOBIANS
 * =================================================================
 */

/*
 * -----------------------------------------------------------------
 * kinLapackDenseInit
 * -----------------------------------------------------------------
 * This routine does remaining initializations specific to the dense
 * linear solver.
 * -----------------------------------------------------------------
 */

static int kinLapackDenseInit(KINMem kin_mem)
{
  KINDlsMem kindls_mem;

  kindls_mem = (KINDlsMem) lmem;
  
  nje   = 0;
  nfeDQ = 0;
  
  if (jacDQ) {
    djac = kinDlsDenseDQJac;
    J_data = kin_mem;
  } else {
    J_data = kin_mem->kin_user_data;
  }

  last_flag = KINDLS_SUCCESS;
  return(0);
}

/*
 * -----------------------------------------------------------------
 * kinLapackDenseSetup
 * -----------------------------------------------------------------
 * This routine does the setup operations for the dense linear solver.
 * It calls the dense LU factorization routine.
 * -----------------------------------------------------------------
 */

static int kinLapackDenseSetup(KINMem kin_mem)
{
  KINDlsMem kindls_mem;
  int ier, retval;
  int intn;

  kindls_mem = (KINDlsMem) lmem;

  intn = (int) n;

  nje++;
  SetToZero(J); 
  retval = djac(n, uu, fval, J, J_data, vtemp1, vtemp2);
  if (retval != 0) {
    last_flag = -1;
    return(-1);
  }

  /* Do LU factorization of J */
  dgetrf_f77(&intn, &intn, J->data, &intn, pivots, &ier);

  /* Return 0 if the LU was complete; otherwise return -1 */
  last_flag = (long int) ier;
  if (ier > 0) return(-1);

  return(0);
}

/*
 * -----------------------------------------------------------------
 * kinLapackDenseSolve
 * -----------------------------------------------------------------
 * This routine handles the solve operation for the dense linear solver
 * by calling the dense backsolve routine.  The returned value is 0.
 * -----------------------------------------------------------------
 */

static int kinLapackDenseSolve(KINMem kin_mem, N_Vector x, N_Vector b, realtype *res_norm)
{
  KINDlsMem kindls_mem;
  realtype *xd;
  int ier, one = 1;
  int intn;

  kindls_mem = (KINDlsMem) lmem;

  intn = (int) n;

  /* Copy the right-hand side into x */
  N_VScale(ONE, b, x);
  xd = N_VGetArrayPointer(x);

  /* Back-solve and get solution in x */
  dgetrs_f77("N", &intn, &one, J->data, &intn, pivots, xd, &intn, &ier, 1); 
  if (ier > 0) return(-1);

  /* Compute the terms Jpnorm and sfdotJp for use in the global strategy
   * routines and in KINForcingTerm. Both of these terms are subsequently
   * corrected if the step is reduced by constraints or the line search.
   *
   * sJpnorm is the norm of the scaled product (scaled by fscale) of
   * the current Jacobian matrix J and the step vector p.
   *
   * sfdotJp is the dot product of the scaled f vector and the scaled
   * vector J*p, where the scaling uses fscale. 
   */
  sJpnorm = N_VWL2Norm(b,fscale);
  N_VProd(b, fscale, b);
  N_VProd(b, fscale, b);
  sfdotJp = N_VDotProd(fval, b);

  last_flag = KINDLS_SUCCESS;

  return(0);
}

/*
 * -----------------------------------------------------------------
 * kinLapackDenseFree
 * -----------------------------------------------------------------
 * This routine frees memory specific to the dense linear solver.
 * -----------------------------------------------------------------
 */

static void kinLapackDenseFree(KINMem kin_mem)
{
  KINDlsMem  kindls_mem;

  kindls_mem = (KINDlsMem) lmem;
  
  DestroyMat(J);
  DestroyArray(pivots);
  free(kindls_mem); kindls_mem = NULL;
}


/* 
 * =================================================================
 *  PRIVATE FUNCTIONS FOR SOLUTION WITH BANDED JACOBIANS
 * =================================================================
 */

/*
 * -----------------------------------------------------------------
 * kinLapackBandInit
 * -----------------------------------------------------------------
 * This routine does remaining initializations specific to the band
 * linear solver.
 * -----------------------------------------------------------------
 */

static int kinLapackBandInit(KINMem kin_mem)
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
 * kinLapackBandSetup
 * -----------------------------------------------------------------
 * This routine does the setup operations for the band linear solver.
 * It makes a decision whether or not to call the Jacobian evaluation
 * routine based on various state variables, and if not it uses the 
 * saved copy.  In any case, it constructs the Newton matrix 
 * M = I - gamma*J, updates counters, and calls the band LU 
 * factorization routine.
 * -----------------------------------------------------------------
 */

static int kinLapackBandSetup(KINMem kin_mem)
{
  KINDlsMem kindls_mem;
  int ier, retval;
  int intn, iml, imu, ldmat;

  kindls_mem = (KINDlsMem) lmem;

  intn = (int) n;
  iml = (int) ml;
  imu = (int) mu;
  ldmat = J->ldim;

  nje++;
  SetToZero(J); 
  retval = bjac(n, mu, ml, uu, fval, J, J_data, vtemp1, vtemp2);
  if (retval != 0) {
    last_flag = -1;
    return(-1);
  }
  
  /* Do LU factorization of J */
  dgbtrf_f77(&intn, &intn, &iml, &imu, J->data, &ldmat, pivots, &ier);

  /* Return 0 if the LU was complete; otherwise return -1 */
  last_flag = (long int) ier;
  if (ier > 0) return(-1);

  return(0);
}

/*
 * -----------------------------------------------------------------
 * kinLapackBandSolve
 * -----------------------------------------------------------------
 * This routine handles the solve operation for the band linear solver
 * by calling the band backsolve routine.  The return value is 0.
 * -----------------------------------------------------------------
 */

static int kinLapackBandSolve(KINMem kin_mem, N_Vector x, N_Vector b, realtype *res_norm)
{
  KINDlsMem kindls_mem;
  realtype *xd;
  int ier, one = 1;
  int intn, iml, imu, ldmat;

  kindls_mem = (KINDlsMem) lmem;

  intn = (int) n;
  iml = (int) ml;
  imu = (int) mu;
  ldmat = J->ldim;

  /* Copy the right-hand side into x */
  N_VScale(ONE, b, x);
  xd = N_VGetArrayPointer(x);

  /* Back-solve and get solution in x */
  dgbtrs_f77("N", &intn, &iml, &imu, &one, J->data, &ldmat, pivots, xd, &intn, &ier, 1);
  if (ier > 0) return(-1);

  /* Compute the terms Jpnorm and sfdotJp for use in the global strategy
   * routines and in KINForcingTerm. Both of these terms are subsequently
   * corrected if the step is reduced by constraints or the line search.
   * 
   * sJpnorm is the norm of the scaled product (scaled by fscale) of
   * the current Jacobian matrix J and the step vector p.
   *
   * sfdotJp is the dot product of the scaled f vector and the scaled
   * vector J*p, where the scaling uses fscale. 
   */
  sJpnorm = N_VWL2Norm(b,fscale);
  N_VProd(b, fscale, b);
  N_VProd(b, fscale, b);
  sfdotJp = N_VDotProd(fval, b);

  last_flag = KINDLS_SUCCESS;

  return(0);
}

/*
 * -----------------------------------------------------------------
 * kinLapackBandFree
 * -----------------------------------------------------------------
 * This routine frees memory specific to the band linear solver.
 * -----------------------------------------------------------------
 */

static void kinLapackBandFree(KINMem kin_mem)
{
  KINDlsMem kindls_mem;

  kindls_mem = (KINDlsMem) lmem;

  DestroyMat(J);
  DestroyArray(pivots);
  free(kindls_mem); kindls_mem = NULL;
}

/*
 * -----------------------------------------------------------------
 * $Revision: 1.15 $
 * $Date: 2011/03/23 22:58:46 $
 * ----------------------------------------------------------------- 
 * Programmer: Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2006, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see the LICENSE file.
 * -----------------------------------------------------------------
 * This is the implementation file for a dense or banded CVODES 
 * linear solver using BLAS and LAPACK functions.
 * -----------------------------------------------------------------
 */

/* 
 * =================================================================
 * IMPORTED HEADER FILES
 * =================================================================
 */

#include <stdio.h>
#include <stdlib.h>

#include <cvodes/cvodes_lapack.h>
#include "cvodes_direct_impl.h"
#include "cvodes_impl.h"

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

/* CVSLAPACK DENSE linit, lsetup, lsolve, and lfree routines */ 
static int cvLapackDenseInit(CVodeMem cv_mem);
static int cvLapackDenseSetup(CVodeMem cv_mem, int convfail, 
                              N_Vector yP, N_Vector fctP, 
                              booleantype *jcurPtr,
                              N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
static int cvLapackDenseSolve(CVodeMem cv_mem, N_Vector b, N_Vector weight,
                              N_Vector yC, N_Vector fctC);
static void cvLapackDenseFree(CVodeMem cv_mem);

/* CVSLAPACK BAND linit, lsetup, lsolve, and lfree routines */ 
static int cvLapackBandInit(CVodeMem cv_mem);
static int cvLapackBandSetup(CVodeMem cv_mem, int convfail, 
                             N_Vector yP, N_Vector fctP, 
                             booleantype *jcurPtr,
                             N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
static int cvLapackBandSolve(CVodeMem cv_mem, N_Vector b, N_Vector weight,
                             N_Vector yC, N_Vector fctC);
static void cvLapackBandFree(CVodeMem cv_mem);

/* CVSLAPACK lfreeB functions */
static void cvLapackDenseFreeB(CVodeBMem cvB_mem);
static void cvLapackBandFreeB(CVodeBMem cvB_mem);

/* 
 * ================================================================
 *
 *                   PART I - forward problems
 *
 * ================================================================
 */

/* Readability Replacements */

#define lmm            (cv_mem->cv_lmm)
#define f              (cv_mem->cv_f)
#define uround         (cv_mem->cv_uround)
#define nst            (cv_mem->cv_nst)
#define tn             (cv_mem->cv_tn)
#define h              (cv_mem->cv_h)
#define gamma          (cv_mem->cv_gamma)
#define gammap         (cv_mem->cv_gammap)
#define gamrat         (cv_mem->cv_gamrat)
#define ewt            (cv_mem->cv_ewt)

#define linit          (cv_mem->cv_linit)
#define lsetup         (cv_mem->cv_lsetup)
#define lsolve         (cv_mem->cv_lsolve)
#define lfree          (cv_mem->cv_lfree)
#define lmem           (cv_mem->cv_lmem)
#define tempv          (cv_mem->cv_tempv)
#define setupNonNull   (cv_mem->cv_setupNonNull)

#define mtype          (cvdls_mem->d_type)
#define n              (cvdls_mem->d_n)
#define jacDQ          (cvdls_mem->d_jacDQ)
#define djac           (cvdls_mem->d_djac)
#define bjac           (cvdls_mem->d_bjac)
#define M              (cvdls_mem->d_M)
#define mu             (cvdls_mem->d_mu)
#define ml             (cvdls_mem->d_ml)
#define smu            (cvdls_mem->d_smu)
#define pivots         (cvdls_mem->d_pivots)
#define savedJ         (cvdls_mem->d_savedJ)
#define nstlj          (cvdls_mem->d_nstlj)
#define nje            (cvdls_mem->d_nje)
#define nfeDQ          (cvdls_mem->d_nfeDQ)
#define J_data         (cvdls_mem->d_J_data)
#define last_flag      (cvdls_mem->d_last_flag)

/*
 * -----------------------------------------------------------------
 * CVLapackDense
 * -----------------------------------------------------------------
 * This routine initializes the memory record and sets various function
 * fields specific to the linear solver module.  CVLapackDense first
 * calls the existing lfree routine if this is not NULL.  Then it sets
 * the cv_linit, cv_lsetup, cv_lsolve, cv_lfree fields in (*cvode_mem)
 * to be cvLapackDenseInit, cvLapackDenseSetup, cvLapackDenseSolve, 
 * and cvLapackDenseFree, respectively.  It allocates memory for a 
 * structure of type CVDlsMemRec and sets the cv_lmem field in 
 * (*cvode_mem) to the address of this structure.  It sets setupNonNull 
 * in (*cvode_mem) to TRUE, and the d_jac field to the default 
 * cvDlsDenseDQJac. Finally, it allocates memory for M, pivots, and 
 * (if needed) savedJ.
 * The return value is SUCCESS = 0, or LMEM_FAIL = -1.
 *
 * NOTE: The dense linear solver assumes a serial implementation
 *       of the NVECTOR package. Therefore, CVLapackDense will first 
 *       test for a compatible N_Vector internal representation 
 *       by checking that N_VGetArrayPointer and N_VSetArrayPointer 
 *       exist.
 * -----------------------------------------------------------------
 */
int CVLapackDense(void *cvode_mem, int N)
{
  CVodeMem cv_mem;
  CVDlsMem cvdls_mem;

  /* Return immediately if cvode_mem is NULL */
  if (cvode_mem == NULL) {
    cvProcessError(NULL, CVDLS_MEM_NULL, "CVSLAPACK", "CVLapackDense", MSGD_CVMEM_NULL);
    return(CVDLS_MEM_NULL);
  }
  cv_mem = (CVodeMem) cvode_mem;

  /* Test if the NVECTOR package is compatible with the CVSLAPACK solver */
  if (tempv->ops->nvgetarraypointer == NULL ||
      tempv->ops->nvsetarraypointer == NULL) {
    cvProcessError(cv_mem, CVDLS_ILL_INPUT, "CVSLAPACK", "CVLapackDense", MSGD_BAD_NVECTOR);
    return(CVDLS_ILL_INPUT);
  }

  if (lfree !=NULL) lfree(cv_mem);

  /* Set four main function fields in cv_mem */
  linit  = cvLapackDenseInit;
  lsetup = cvLapackDenseSetup;
  lsolve = cvLapackDenseSolve;
  lfree  = cvLapackDenseFree;

  /* Get memory for CVDlsMemRec */
  cvdls_mem = NULL;
  cvdls_mem = (CVDlsMem) malloc(sizeof(struct CVDlsMemRec));
  if (cvdls_mem == NULL) {
    cvProcessError(cv_mem, CVDLS_MEM_FAIL, "CVSLAPACK", "CVLapackDense", MSGD_MEM_FAIL);
    return(CVDLS_MEM_FAIL);
  }

  /* Set matrix type */
  mtype = SUNDIALS_DENSE;

  /* Initialize Jacobian-related data */
  jacDQ  = TRUE;
  djac   = NULL;
  J_data = NULL;

  last_flag = CVDLS_SUCCESS;
  setupNonNull = TRUE;

  /* Set problem dimension */
  n = (long int) N;

  /* Allocate memory for M, pivot array, and (if needed) savedJ */
  M = NULL;
  pivots = NULL;
  savedJ = NULL;

  M = NewDenseMat(n, n);
  if (M == NULL) {
    cvProcessError(cv_mem, CVDLS_MEM_FAIL, "CVSLAPACK", "CVLapackDense", MSGD_MEM_FAIL);
    free(cvdls_mem); cvdls_mem = NULL;
    return(CVDLS_MEM_FAIL);
  }
  pivots = NewIntArray(N);
  if (pivots == NULL) {
    cvProcessError(cv_mem, CVDLS_MEM_FAIL, "CVSLAPACK", "CVLapackDense", MSGD_MEM_FAIL);
    DestroyMat(M);
    free(cvdls_mem); cvdls_mem = NULL;
    return(CVDLS_MEM_FAIL);
  }
  savedJ = NewDenseMat(n, n);
  if (savedJ == NULL) {
    cvProcessError(cv_mem, CVDLS_MEM_FAIL, "CVSLAPACK", "CVLapackDense", MSGD_MEM_FAIL);
    DestroyMat(M);
    DestroyArray(pivots);
    free(cvdls_mem); cvdls_mem = NULL;
    return(CVDLS_MEM_FAIL);
  }

  /* Attach linear solver memory to integrator memory */
  lmem = cvdls_mem;

  return(CVDLS_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * CVLapackBand
 * -----------------------------------------------------------------
 * This routine initializes the memory record and sets various function
 * fields specific to the band linear solver module. It first calls
 * the existing lfree routine if this is not NULL.  It then sets the
 * cv_linit, cv_lsetup, cv_lsolve, and cv_lfree fields in (*cvode_mem)
 * to be cvLapackBandInit, cvLapackBandSetup, cvLapackBandSolve, 
 * and cvLapackBandFree, respectively.  It allocates memory for a 
 * structure of type CVLapackBandMemRec and sets the cv_lmem field in 
 * (*cvode_mem) to the address of this structure.  It sets setupNonNull 
 * in (*cvode_mem) to be TRUE, mu to be mupper, ml to be mlower, and 
 * the jacE and jacI field to NULL.
 * Finally, it allocates memory for M, pivots, and (if needed) savedJ.  
 * The CVLapackBand return value is CVDLS_SUCCESS = 0, 
 * CVDLS_MEM_FAIL = -1, or CVDLS_ILL_INPUT = -2.
 *
 * NOTE: The CVSLAPACK linear solver assumes a serial implementation
 *       of the NVECTOR package. Therefore, CVLapackBand will first 
 *       test for compatible a compatible N_Vector internal
 *       representation by checking that the function 
 *       N_VGetArrayPointer exists.
 * -----------------------------------------------------------------
 */                  
int CVLapackBand(void *cvode_mem, int N, int mupper, int mlower)
{
  CVodeMem cv_mem;
  CVDlsMem cvdls_mem;

  /* Return immediately if cvode_mem is NULL */
  if (cvode_mem == NULL) {
    cvProcessError(NULL, CVDLS_MEM_NULL, "CVSLAPACK", "CVLapackBand", MSGD_CVMEM_NULL);
    return(CVDLS_MEM_NULL);
  }
  cv_mem = (CVodeMem) cvode_mem;

  /* Test if the NVECTOR package is compatible with the BAND solver */
  if (tempv->ops->nvgetarraypointer == NULL) {
    cvProcessError(cv_mem, CVDLS_ILL_INPUT, "CVSLAPACK", "CVLapackBand", MSGD_BAD_NVECTOR);
    return(CVDLS_ILL_INPUT);
  }

  if (lfree != NULL) lfree(cv_mem);

  /* Set four main function fields in cv_mem */  
  linit  = cvLapackBandInit;
  lsetup = cvLapackBandSetup;
  lsolve = cvLapackBandSolve;
  lfree  = cvLapackBandFree;
  
  /* Get memory for CVDlsMemRec */
  cvdls_mem = NULL;
  cvdls_mem = (CVDlsMem) malloc(sizeof(struct CVDlsMemRec));
  if (cvdls_mem == NULL) {
    cvProcessError(cv_mem, CVDLS_MEM_FAIL, "CVSLAPACK", "CVLapackBand", MSGD_MEM_FAIL);
    return(CVDLS_MEM_FAIL);
  }

  /* Set matrix type */
  mtype = SUNDIALS_BAND;

  /* Initialize Jacobian-related data */
  jacDQ  = TRUE;
  bjac   = NULL;
  J_data = NULL;

  last_flag = CVDLS_SUCCESS;
  setupNonNull = TRUE;
  
  /* Load problem dimension */
  n = (long int) N;

  /* Load half-bandwiths in cvdls_mem */
  ml = (long int) mlower;
  mu = (long int) mupper;

  /* Test ml and mu for legality */
  if ((ml < 0) || (mu < 0) || (ml >= n) || (mu >= n)) {
    cvProcessError(cv_mem, CVDLS_ILL_INPUT, "CVSLAPACK", "CVLapackBand", MSGD_BAD_SIZES);
    free(cvdls_mem); cvdls_mem = NULL;
    return(CVDLS_ILL_INPUT);
  }

  /* Set extended upper half-bandwith for M (required for pivoting) */
  smu = MIN(n-1, mu + ml);

  /* Allocate memory for M, savedJ, and pivot arrays */
  M = NULL;
  pivots = NULL;
  savedJ = NULL;

  M = NewBandMat(n, mu, ml, smu);
  if (M == NULL) {
    cvProcessError(cv_mem, CVDLS_MEM_FAIL, "CVSLAPACK", "CVLapackBand", MSGD_MEM_FAIL);
    free(cvdls_mem); cvdls_mem = NULL;
    return(CVDLS_MEM_FAIL);
  }  
  pivots = NewIntArray(N);
  if (pivots == NULL) {
    cvProcessError(cv_mem, CVDLS_MEM_FAIL, "CVSLAPACK", "CVLapackBand", MSGD_MEM_FAIL);
    DestroyMat(M);
    free(cvdls_mem); cvdls_mem = NULL;
    return(CVDLS_MEM_FAIL);
  }
  savedJ = NewBandMat(n, mu, ml, smu);
  if (savedJ == NULL) {
    cvProcessError(cv_mem, CVDLS_MEM_FAIL, "CVSLAPACK", "CVLapackBand", MSGD_MEM_FAIL);
    DestroyMat(M);
    DestroyArray(pivots);
    free(cvdls_mem); cvdls_mem = NULL;
    return(CVDLS_MEM_FAIL);
  }

  /* Attach linear solver memory to integrator memory */
  lmem = cvdls_mem;

  return(CVDLS_SUCCESS);
}

/* 
 * =================================================================
 *  PRIVATE FUNCTIONS FOR IMPLICIT INTEGRATION WITH DENSE JACOBIANS
 * =================================================================
 */

/*
 * cvLapackDenseInit does remaining initializations specific to the dense
 * linear solver.
 */
static int cvLapackDenseInit(CVodeMem cv_mem)
{
  CVDlsMem cvdls_mem;

  cvdls_mem = (CVDlsMem) lmem;
  
  nje   = 0;
  nfeDQ = 0;
  nstlj = 0;

  /* Set Jacobian function and data, depending on jacDQ */
  if (jacDQ) {
    djac = cvDlsDenseDQJac;
    J_data = cv_mem;
  } else {
    J_data = cv_mem->cv_user_data;
  }

  last_flag = CVDLS_SUCCESS;
  return(0);
}

/*
 * cvLapackDenseSetup does the setup operations for the dense linear solver.
 * It makes a decision whether or not to call the Jacobian evaluation
 * routine based on various state variables, and if not it uses the 
 * saved copy. In any case, it constructs the Newton matrix M = I - gamma*J
 * updates counters, and calls the dense LU factorization routine.
 */
static int cvLapackDenseSetup(CVodeMem cv_mem, int convfail,
                              N_Vector yP, N_Vector fctP,
                              booleantype *jcurPtr,
                              N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  CVDlsMem cvdls_mem;
  realtype dgamma, fact;
  booleantype jbad, jok;
  int ier, retval, one = 1;
  int intn, lenmat;

  cvdls_mem = (CVDlsMem) lmem;
  intn = (int) n;
  lenmat = M->ldata;

  /* Use nst, gamma/gammap, and convfail to set J eval. flag jok */
  dgamma = ABS((gamma/gammap) - ONE);
  jbad = (nst == 0) || (nst > nstlj + CVD_MSBJ) ||
    ((convfail == CV_FAIL_BAD_J) && (dgamma < CVD_DGMAX)) ||
    (convfail == CV_FAIL_OTHER);
  jok = !jbad;
  
  if (jok) {
    
    /* If jok = TRUE, use saved copy of J */
    *jcurPtr = FALSE;
    dcopy_f77(&lenmat, savedJ->data, &one, M->data, &one);
    
  } else {
    
    /* If jok = FALSE, call jac routine for new J value */
    nje++;
    nstlj = nst;
    *jcurPtr = TRUE;
    SetToZero(M);

    retval = djac(n, tn, yP, fctP, M, J_data, tmp1, tmp2, tmp3);
    if (retval == 0) {
      dcopy_f77(&lenmat, M->data, &one, savedJ->data, &one);
    } else if (retval < 0) {
      cvProcessError(cv_mem, CVDLS_JACFUNC_UNRECVR, "CVSLAPACK", "cvLapackDenseSetup", MSGD_JACFUNC_FAILED);
      last_flag = CVDLS_JACFUNC_UNRECVR;
      return(-1);
    } else if (retval > 0) {
      last_flag = CVDLS_JACFUNC_RECVR;
      return(1);
    }
    
  }

  /* Scale J by - gamma */
  fact = -gamma;
  dscal_f77(&lenmat, &fact, M->data, &one);
  
  /* Add identity to get M = I - gamma*J*/
  AddIdentity(M);

  /* Do LU factorization of M */
  dgetrf_f77(&intn, &intn, M->data, &intn, pivots, &ier);

  /* Return 0 if the LU was complete; otherwise return 1 */
  last_flag = (long int) ier;
  if (ier > 0) return(1);
  return(0);
}

/*
 * cvLapackDenseSolve handles the solve operation for the dense linear solver
 * by calling the dense backsolve routine.
 */
static int cvLapackDenseSolve(CVodeMem cv_mem, N_Vector b, N_Vector weight,
                              N_Vector yC, N_Vector fctC)
{
  CVDlsMem cvdls_mem;
  realtype *bd, fact;
  int ier, one = 1;
  int intn;

  cvdls_mem = (CVDlsMem) lmem;

  intn = (int) n;

  bd = N_VGetArrayPointer(b);

  dgetrs_f77("N", &intn, &one, M->data, &intn, pivots, bd, &intn, &ier, 1); 
  if (ier > 0) return(1);

  /* For BDF, scale the correction to account for change in gamma */
  if ((lmm == CV_BDF) && (gamrat != ONE)) {
    fact = TWO/(ONE + gamrat);
    dscal_f77(&intn, &fact, bd, &one); 
  }
  
  last_flag = CVDLS_SUCCESS;
  return(0);
}

/*
 * cvLapackDenseFree frees memory specific to the dense linear solver.
 */
static void cvLapackDenseFree(CVodeMem cv_mem)
{
  CVDlsMem  cvdls_mem;

  cvdls_mem = (CVDlsMem) lmem;
  
  DestroyMat(M);
  DestroyArray(pivots);
  DestroyMat(savedJ);
  free(cvdls_mem); 
  cvdls_mem = NULL;
}

/* 
 * =================================================================
 *  PRIVATE FUNCTIONS FOR IMPLICIT INTEGRATION WITH BAND JACOBIANS
 * =================================================================
 */

/*
 * cvLapackBandInit does remaining initializations specific to the band
 * linear solver.
 */
static int cvLapackBandInit(CVodeMem cv_mem)
{
  CVDlsMem cvdls_mem;

  cvdls_mem = (CVDlsMem) lmem;

  nje   = 0;
  nfeDQ = 0;
  nstlj = 0;

  /* Set Jacobian function and data, depending on jacDQ */
  if (jacDQ) {
    bjac = cvDlsBandDQJac;
    J_data = cv_mem;
  } else {
    J_data = cv_mem->cv_user_data;
  }

  last_flag = CVDLS_SUCCESS;
  return(0);
}

/*
 * cvLapackBandSetup does the setup operations for the band linear solver.
 * It makes a decision whether or not to call the Jacobian evaluation
 * routine based on various state variables, and if not it uses the 
 * saved copy. In any case, it constructs the Newton matrix M = I - gamma*J, 
 * updates counters, and calls the band LU factorization routine.
 */
static int cvLapackBandSetup(CVodeMem cv_mem, int convfail, 
                             N_Vector yP, N_Vector fctP, 
                             booleantype *jcurPtr,
                             N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  CVDlsMem cvdls_mem;
  realtype dgamma, fact;
  booleantype jbad, jok;
  int ier, retval, one = 1;
  int intn, iml, imu, lenmat, ldmat;

  cvdls_mem = (CVDlsMem) lmem;

  intn = (int) n;
  iml = (int) ml;
  imu = (int) mu;
  lenmat = M->ldata;
  ldmat = M->ldim;

  /* Use nst, gamma/gammap, and convfail to set J eval. flag jok */
  dgamma = ABS((gamma/gammap) - ONE);
  jbad = (nst == 0) || (nst > nstlj + CVD_MSBJ) ||
    ((convfail == CV_FAIL_BAD_J) && (dgamma < CVD_DGMAX)) ||
    (convfail == CV_FAIL_OTHER);
  jok = !jbad;
  
  if (jok) {
    
    /* If jok = TRUE, use saved copy of J */
    *jcurPtr = FALSE;
    dcopy_f77(&lenmat, savedJ->data, &one, M->data, &one);
    
  } else {
    
    /* If jok = FALSE, call jac routine for new J value */
    nje++;
    nstlj = nst;
    *jcurPtr = TRUE;
    SetToZero(M); 

    retval = bjac(n, mu, ml, tn, yP, fctP, M, J_data, tmp1, tmp2, tmp3);
    if (retval == 0) {
      dcopy_f77(&lenmat, M->data, &one, savedJ->data, &one);
    } else if (retval < 0) {
      cvProcessError(cv_mem, CVDLS_JACFUNC_UNRECVR, "CVSLAPACK", "cvLapackBandSetup", MSGD_JACFUNC_FAILED);
      last_flag = CVDLS_JACFUNC_UNRECVR;
      return(-1);
    } else if (retval > 0) {
      last_flag = CVDLS_JACFUNC_RECVR;
      return(1);
    }
    
  }
  
  /* Scale J by - gamma */
  fact = -gamma;
  dscal_f77(&lenmat, &fact, M->data, &one);
  
  /* Add identity to get M = I - gamma*J*/
  AddIdentity(M);
  
  /* Do LU factorization of M */
  dgbtrf_f77(&intn, &intn, &iml, &imu, M->data, &ldmat, pivots, &ier);

  /* Return 0 if the LU was complete; otherwise return 1 */
  last_flag = (long int) ier;
  if (ier > 0) return(1);
  return(0);

}

/*
 * cvLapackBandSolve handles the solve operation for the band linear solver
 * by calling the band backsolve routine.
 */
static int cvLapackBandSolve(CVodeMem cv_mem, N_Vector b, N_Vector weight,
                             N_Vector yC, N_Vector fctC)
{
  CVDlsMem cvdls_mem;
  realtype *bd, fact;
  int ier, one = 1;
  int intn, iml, imu, ldmat;

  cvdls_mem = (CVDlsMem) lmem;

  intn = (int) n;
  iml = (int) ml;
  imu = (int) mu;
  ldmat = M->ldim;

  bd = N_VGetArrayPointer(b);

  dgbtrs_f77("N", &intn, &iml, &imu, &one, M->data, &ldmat, pivots, bd, &intn, &ier, 1);
  if (ier > 0) return(1);

  /* For BDF, scale the correction to account for change in gamma */
  if ((lmm == CV_BDF) && (gamrat != ONE)) {
    fact = TWO/(ONE + gamrat);
    dscal_f77(&intn, &fact, bd, &one); 
  }

  last_flag = CVDLS_SUCCESS;
  return(0);
}

/*
 * cvLapackBandFree frees memory specific to the band linear solver.
 */
static void cvLapackBandFree(CVodeMem cv_mem)
{
  CVDlsMem  cvdls_mem;

  cvdls_mem = (CVDlsMem) lmem;
  
  DestroyMat(M);
  DestroyArray(pivots);
  DestroyMat(savedJ);
  free(cvdls_mem); 
  cvdls_mem = NULL;
}

/* 
 * ================================================================
 *
 *                   PART II - backward problems
 *
 * ================================================================
 */

/*
 * CVLapackDenseB is a wraper around CVLapackDense. It attaches the
 * dense CVSLAPACK linear solver to the backward problem memory block.
 */

int CVLapackDenseB(void *cvode_mem, int which, int nB)
{
  CVodeMem cv_mem;
  CVadjMem ca_mem;
  CVodeBMem cvB_mem;
  void *cvodeB_mem;
  CVDlsMemB cvdlsB_mem;
  int flag;

  /* Check if cvode_mem exists */
  if (cvode_mem == NULL) {
    cvProcessError(NULL, CVDLS_MEM_NULL, "CVSLAPACK", "CVLapackDenseB", MSGD_CVMEM_NULL);
    return(CVDLS_MEM_NULL);
  }
  cv_mem = (CVodeMem) cvode_mem;

  /* Was ASA initialized? */
  if (cv_mem->cv_adjMallocDone == FALSE) {
    cvProcessError(cv_mem, CVDLS_NO_ADJ, "CVSLAPACK", "CVLapackDenseB", MSGD_NO_ADJ);
    return(CVDLS_NO_ADJ);
  } 
  ca_mem = cv_mem->cv_adj_mem;

  /* Check which */
  if ( which >= ca_mem->ca_nbckpbs ) {
    cvProcessError(cv_mem, CVDLS_ILL_INPUT, "CVSLAPACK", "CVLapackDenseB", MSGCV_BAD_WHICH);
    return(CVDLS_ILL_INPUT);
  }

  /* Find the CVodeBMem entry in the linked list corresponding to which */
  cvB_mem = ca_mem->cvB_mem;
  while (cvB_mem != NULL) {
    if ( which == cvB_mem->cv_index ) break;
    cvB_mem = cvB_mem->cv_next;
  }

  cvodeB_mem = (void *) (cvB_mem->cv_mem);

  /* Get memory for CVDlsMemRecB */
  cvdlsB_mem = (CVDlsMemB) malloc(sizeof(struct CVDlsMemRecB));
  if (cvdlsB_mem == NULL) {
    cvProcessError(cv_mem, CVDLS_MEM_FAIL, "CVSLAPACK", "CVLapackDenseB", MSGD_MEM_FAIL);
    return(CVDLS_MEM_FAIL);
  }

  /* set matrix type */
  cvdlsB_mem->d_typeB = SUNDIALS_DENSE;

  /* initialize Jacobian function */
  cvdlsB_mem->d_djacB = NULL;

  /* attach lmemB and lfreeB */
  cvB_mem->cv_lmem = cvdlsB_mem;
  cvB_mem->cv_lfree = cvLapackDenseFreeB;

  flag = CVLapackDense(cvodeB_mem, nB);

  if (flag != CVDLS_SUCCESS) {
    free(cvdlsB_mem);
    cvdlsB_mem = NULL;
  }

  return(flag);
}

/*
 * cvLapackDenseFreeB frees the memory associated with the dense CVSLAPACK
 * linear solver for backward integration.
 */

static void cvLapackDenseFreeB(CVodeBMem cvB_mem)
{
  CVDlsMemB cvdlsB_mem;

  cvdlsB_mem = (CVDlsMemB) (cvB_mem->cv_lmem);

  free(cvdlsB_mem);
}

/*
 * CVLapackBandB is a wraper around CVLapackBand. It attaches the band
 * CVSLAPACK linear solver to the backward problem memory block.
 */

int CVLapackBandB(void *cvode_mem, int which,
                  int nB, int mupperB, int mlowerB)
{
  CVodeMem cv_mem;
  CVadjMem ca_mem;
  CVodeBMem cvB_mem;
  void *cvodeB_mem;
  CVDlsMemB cvdlsB_mem;
  int flag;

  /* Check if cvode_mem exists */
  if (cvode_mem == NULL) {
    cvProcessError(NULL, CVDLS_MEM_NULL, "CVSLAPACK", "CVLapackBandB", MSGD_CVMEM_NULL);
    return(CVDLS_MEM_NULL);
  }
  cv_mem = (CVodeMem) cvode_mem;

  /* Was ASA initialized? */
  if (cv_mem->cv_adjMallocDone == FALSE) {
    cvProcessError(cv_mem, CVDLS_NO_ADJ, "CVSLAPACK", "CVLapackBandB", MSGD_NO_ADJ);
    return(CVDLS_NO_ADJ);
  } 
  ca_mem = cv_mem->cv_adj_mem;

  /* Check which */
  if ( which >= ca_mem->ca_nbckpbs ) {
    cvProcessError(cv_mem, CVDLS_ILL_INPUT, "CVSLAPACK", "CVLapackBandB", MSGCV_BAD_WHICH);
    return(CVDLS_ILL_INPUT);
  }

  /* Find the CVodeBMem entry in the linked list corresponding to which */
  cvB_mem = ca_mem->cvB_mem;
  while (cvB_mem != NULL) {
    if ( which == cvB_mem->cv_index ) break;
    cvB_mem = cvB_mem->cv_next;
  }

  cvodeB_mem = (void *) (cvB_mem->cv_mem);

  /* Get memory for CVDlsMemRecB */
  cvdlsB_mem = (CVDlsMemB) malloc(sizeof(struct CVDlsMemRecB));
  if (cvdlsB_mem == NULL) {
    cvProcessError(cv_mem, CVDLS_MEM_FAIL, "CVSLAPACK", "CVLapackBandB", MSGD_MEM_FAIL);
    return(CVDLS_MEM_FAIL);
  }

  /* set matrix type */
  cvdlsB_mem->d_typeB = SUNDIALS_BAND;

  /* initialize Jacobian function */
  cvdlsB_mem->d_bjacB = NULL;

  /* attach lmemB and lfreeB */
  cvB_mem->cv_lmem = cvdlsB_mem;
  cvB_mem->cv_lfree = cvLapackBandFreeB;

  flag = CVLapackBand(cvodeB_mem, nB, mupperB, mlowerB);

  if (flag != CVDLS_SUCCESS) {
    free(cvdlsB_mem);
    cvdlsB_mem = NULL;
  }

  return(flag);
}

/*
 * cvLapackBandFreeB frees the memory associated with the band CVSLAPACK
 * linear solver for backward integration.
 */

static void cvLapackBandFreeB(CVodeBMem cvB_mem)
{
  CVDlsMemB cvdlsB_mem;

  cvdlsB_mem = (CVDlsMemB) (cvB_mem->cv_lmem);

  free(cvdlsB_mem);
}

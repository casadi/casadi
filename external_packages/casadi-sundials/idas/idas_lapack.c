/*
 * -----------------------------------------------------------------
 * $Revision: 1.12 $
 * $Date: 2011/03/23 23:25:35 $
 * ----------------------------------------------------------------- 
 * Programmer: Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2006, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see the LICENSE file.
 * -----------------------------------------------------------------
 * This is the implementation file for a IDAS dense linear solver
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

#include <idas/idas_lapack.h>
#include "idas_direct_impl.h"
#include "idas_impl.h"

#include <sundials/sundials_math.h>

/* 
 * =================================================================
 * FUNCTION SPECIFIC CONSTANTS
 * =================================================================
 */

#define ZERO         RCONST(0.0)
#define ONE          RCONST(1.0)
#define TWO          RCONST(2.0)

/* 
 * =================================================================
 * PROTOTYPES FOR PRIVATE FUNCTIONS
 * =================================================================
 */

/* IDALAPACK DENSE linit, lsetup, lsolve, and lfree routines */ 
static int idaLapackDenseInit(IDAMem IDA_mem);
static int idaLapackDenseSetup(IDAMem IDA_mem,
                               N_Vector yP, N_Vector ypP, N_Vector fctP, 
                               N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
static int idaLapackDenseSolve(IDAMem IDA_mem, N_Vector b, N_Vector weight,
                               N_Vector yC, N_Vector ypC, N_Vector fctC);
static int idaLapackDenseFree(IDAMem IDA_mem);

/* IDALAPACK BAND linit, lsetup, lsolve, and lfree routines */ 
static int idaLapackBandInit(IDAMem IDA_mem);
static int idaLapackBandSetup(IDAMem IDA_mem,
                              N_Vector yP, N_Vector ypP, N_Vector fctP, 
                              N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
static int idaLapackBandSolve(IDAMem IDA_mem, N_Vector b, N_Vector weight,
                              N_Vector yC, N_Vector ypC, N_Vector fctC);
static int idaLapackBandFree(IDAMem IDA_mem);

/*
 * =================================================================
 * READIBILITY REPLACEMENTS
 * =================================================================
 */

#define res            (IDA_mem->ida_res)
#define nst            (IDA_mem->ida_nst)
#define tn             (IDA_mem->ida_tn)
#define hh             (IDA_mem->ida_hh)
#define cj             (IDA_mem->ida_cj)
#define cjratio        (IDA_mem->ida_cjratio)
#define ewt            (IDA_mem->ida_ewt)
#define constraints    (IDA_mem->ida_constraints)

#define linit          (IDA_mem->ida_linit)
#define lsetup         (IDA_mem->ida_lsetup)
#define lsolve         (IDA_mem->ida_lsolve)
#define lfree          (IDA_mem->ida_lfree)
#define lperf          (IDA_mem->ida_lperf)
#define lmem           (IDA_mem->ida_lmem)
#define tempv          (IDA_mem->ida_tempv1)
#define setupNonNull   (IDA_mem->ida_setupNonNull)

#define mtype          (idadls_mem->d_type)
#define n              (idadls_mem->d_n)
#define ml             (idadls_mem->d_ml)
#define mu             (idadls_mem->d_mu)
#define smu            (idadls_mem->d_smu)
#define jacDQ          (idadls_mem->d_jacDQ)
#define djac           (idadls_mem->d_djac)
#define bjac           (idadls_mem->d_bjac)
#define JJ             (idadls_mem->d_J)
#define pivots         (idadls_mem->d_pivots)
#define nje            (idadls_mem->d_nje)
#define nreDQ          (idadls_mem->d_nreDQ)
#define J_data         (idadls_mem->d_J_data)
#define last_flag      (idadls_mem->d_last_flag)

/* 
 * =================================================================
 * EXPORTED FUNCTIONS FOR IMPLICIT INTEGRATION
 * =================================================================
 */
              
/*
 * -----------------------------------------------------------------
 * IDALapackDense
 * -----------------------------------------------------------------
 * This routine initializes the memory record and sets various function
 * fields specific to the linear solver module.  IDALapackDense first
 * calls the existing lfree routine if this is not NULL.  Then it sets
 * the ida_linit, ida_lsetup, ida_lsolve, ida_lfree fields in (*ida_mem)
 * to be idaLapackDenseInit, idaLapackDenseSetup, idaLapackDenseSolve, 
 * and idaLapackDenseFree, respectively.  It allocates memory for a 
 * structure of type IDADlsMemRec and sets the ida_lmem field in 
 * (*ida_mem) to the address of this structure.  It sets setupNonNull 
 * in (*ida_mem) to TRUE, and the d_jac field to the default 
 * idaLapackDenseDQJac. Finally, it allocates memory for M, pivots.
 *
 * The return value is SUCCESS = 0, or LMEM_FAIL = -1.
 *
 * NOTE: The dense linear solver assumes a serial implementation
 *       of the NVECTOR package. Therefore, IDALapackDense will first 
 *       test for a compatible N_Vector internal representation 
 *       by checking that N_VGetArrayPointer and N_VSetArrayPointer 
 *       exist.
 * -----------------------------------------------------------------
 */
int IDALapackDense(void *ida_mem, int N)
{
  IDAMem IDA_mem;
  IDADlsMem idadls_mem;

  /* Return immediately if ida_mem is NULL */
  if (ida_mem == NULL) {
    IDAProcessError(NULL, IDADLS_MEM_NULL, "IDASLAPACK", "IDALapackDense", MSGD_IDAMEM_NULL);
    return(IDADLS_MEM_NULL);
  }
  IDA_mem = (IDAMem) ida_mem;

  /* Test if the NVECTOR package is compatible with the LAPACK solver */
  if (tempv->ops->nvgetarraypointer == NULL ||
      tempv->ops->nvsetarraypointer == NULL) {
    IDAProcessError(IDA_mem, IDADLS_ILL_INPUT, "IDASLAPACK", "IDALapackDense", MSGD_BAD_NVECTOR);
    return(IDADLS_ILL_INPUT);
  }

  if (lfree !=NULL) lfree(IDA_mem);

  /* Set four main function fields in IDA_mem */
  linit  = idaLapackDenseInit;
  lsetup = idaLapackDenseSetup;
  lsolve = idaLapackDenseSolve;
  lperf  = NULL;
  lfree  = idaLapackDenseFree;

  /* Get memory for IDADlsMemRec */
  idadls_mem = NULL;
  idadls_mem = (IDADlsMem) malloc(sizeof(struct IDADlsMemRec));
  if (idadls_mem == NULL) {
    IDAProcessError(IDA_mem, IDADLS_MEM_FAIL, "IDASLAPACK", "IDALapackDense", MSGD_MEM_FAIL);
    return(IDADLS_MEM_FAIL);
  }

  /* Set matrix type */
  mtype = SUNDIALS_DENSE;

  /* Set default Jacobian routine and Jacobian data */
  jacDQ  = TRUE;
  djac   = NULL;
  J_data = NULL;

  last_flag = IDADLS_SUCCESS;
  setupNonNull = TRUE;

  /* Set problem dimension */
  n = (long int) N;

  /* Allocate memory for JJ and pivot array */
  JJ = NULL;
  pivots = NULL;

  JJ = NewDenseMat(n, n);
  if (JJ == NULL) {
    IDAProcessError(IDA_mem, IDADLS_MEM_FAIL, "IDASLAPACK", "IDALapackDense", MSGD_MEM_FAIL);
    free(idadls_mem); idadls_mem = NULL;
    return(IDADLS_MEM_FAIL);
  }
  pivots = NewIntArray(N);
  if (pivots == NULL) {
    IDAProcessError(IDA_mem, IDADLS_MEM_FAIL, "IDASLAPACK", "IDALapackDense", MSGD_MEM_FAIL);
    DestroyMat(JJ);
    free(idadls_mem); idadls_mem = NULL;
    return(IDADLS_MEM_FAIL);
  }

  /* Attach linear solver memory to integrator memory */
  lmem = idadls_mem;

  return(IDADLS_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * IDALapackBand
 * -----------------------------------------------------------------
 * This routine initializes the memory record and sets various function
 * fields specific to the band linear solver module. It first calls
 * the existing lfree routine if this is not NULL.  It then sets the
 * ida_linit, ida_lsetup, ida_lsolve, and ida_lfree fields in (*ida_mem)
 * to be idaLapackBandInit, idaLapackBandSetup, idaLapackBandSolve, 
 * and idaLapackBandFree, respectively.  It allocates memory for a 
 * structure of type IDALapackBandMemRec and sets the ida_lmem field in 
 * (*ida_mem) to the address of this structure.  It sets setupNonNull 
 * in (*ida_mem) to be TRUE, mu to be mupper, ml to be mlower, and 
 * the jacE and jacI field to NULL.
 * Finally, it allocates memory for M and pivots.
 * The IDALapackBand return value is IDADLS_SUCCESS = 0, 
 * IDADLS_MEM_FAIL = -1, or IDADLS_ILL_INPUT = -2.
 *
 * NOTE: The IDALAPACK linear solver assumes a serial implementation
 *       of the NVECTOR package. Therefore, IDALapackBand will first 
 *       test for compatible a compatible N_Vector internal
 *       representation by checking that the function 
 *       N_VGetArrayPointer exists.
 * -----------------------------------------------------------------
 */                  
int IDALapackBand(void *ida_mem, int N, int mupper, int mlower)
{
  IDAMem IDA_mem;
  IDADlsMem idadls_mem;

  /* Return immediately if ida_mem is NULL */
  if (ida_mem == NULL) {
    IDAProcessError(NULL, IDADLS_MEM_NULL, "IDASLAPACK", "IDALapackBand", MSGD_IDAMEM_NULL);
    return(IDADLS_MEM_NULL);
  }
  IDA_mem = (IDAMem) ida_mem;

  /* Test if the NVECTOR package is compatible with the BAND solver */
  if (tempv->ops->nvgetarraypointer == NULL) {
    IDAProcessError(IDA_mem, IDADLS_ILL_INPUT, "IDASLAPACK", "IDALapackBand", MSGD_BAD_NVECTOR);
    return(IDADLS_ILL_INPUT);
  }

  if (lfree != NULL) lfree(IDA_mem);

  /* Set four main function fields in IDA_mem */  
  linit  = idaLapackBandInit;
  lsetup = idaLapackBandSetup;
  lsolve = idaLapackBandSolve;
  lperf  = NULL;
  lfree  = idaLapackBandFree;
  
  /* Get memory for IDADlsMemRec */
  idadls_mem = NULL;
  idadls_mem = (IDADlsMem) malloc(sizeof(struct IDADlsMemRec));
  if (idadls_mem == NULL) {
    IDAProcessError(IDA_mem, IDADLS_MEM_FAIL, "IDASLAPACK", "IDALapackBand", MSGD_MEM_FAIL);
    return(IDADLS_MEM_FAIL);
  }

  /* Set matrix type */
  mtype = SUNDIALS_BAND;

  /* Set default Jacobian routine and Jacobian data */
  jacDQ  = TRUE;
  bjac   = NULL;
  J_data = NULL;

  last_flag = IDADLS_SUCCESS;
  setupNonNull = TRUE;
  
  /* Load problem dimension */
  n = (long int) N;

  /* Load half-bandwiths in idadls_mem */
  ml = (long int) mlower;
  mu = (long int) mupper;

  /* Test ml and mu for legality */
  if ((ml < 0) || (mu < 0) || (ml >= n) || (mu >= n)) {
    IDAProcessError(IDA_mem, IDADLS_ILL_INPUT, "IDASLAPACK", "IDALapackBand", MSGD_BAD_SIZES);
    free(idadls_mem); idadls_mem = NULL;
    return(IDADLS_ILL_INPUT);
  }

  /* Set extended upper half-bandwith for M (required for pivoting) */
  smu = MIN(n-1, mu + ml);

  /* Allocate memory for JJ and pivot arrays */
  JJ = NULL;
  pivots = NULL;

  JJ = NewBandMat(n, mu, ml, smu);
  if (JJ == NULL) {
    IDAProcessError(IDA_mem, IDADLS_MEM_FAIL, "IDASLAPACK", "IDALapackBand", MSGD_MEM_FAIL);
    free(idadls_mem); idadls_mem = NULL;
    return(IDADLS_MEM_FAIL);
  }  
  pivots = NewIntArray(N);
  if (pivots == NULL) {
    IDAProcessError(IDA_mem, IDADLS_MEM_FAIL, "IDASLAPACK", "IDALapackBand", MSGD_MEM_FAIL);
    DestroyMat(JJ);
    free(idadls_mem); idadls_mem = NULL;
    return(IDADLS_MEM_FAIL);
  }

  /* Attach linear solver memory to integrator memory */
  lmem = idadls_mem;

  return(IDADLS_SUCCESS);
}

/* 
 * =================================================================
 *  PRIVATE FUNCTIONS FOR IMPLICIT INTEGRATION WITH DENSE JACOBIANS
 * =================================================================
 */

/*
 * idaLapackDenseInit does remaining initializations specific to the dense
 * linear solver.
 */
static int idaLapackDenseInit(IDAMem IDA_mem)
{
  IDADlsMem idadls_mem;

  idadls_mem = (IDADlsMem) lmem;
  
  nje   = 0;
  nreDQ = 0;
  
  if (jacDQ) {
    djac = idaDlsDenseDQJac;
    J_data = IDA_mem;
  } else {
    J_data = IDA_mem->ida_user_data;
  }

  last_flag = IDADLS_SUCCESS;
  return(0);
}

/*
 * idaLapackDenseSetup does the setup operations for the dense linear solver. 
 * It calls the Jacobian function to obtain the Newton matrix M = F_y + c_j*F_y', 
 * updates counters, and calls the dense LU factorization routine.
 */
static int idaLapackDenseSetup(IDAMem IDA_mem,
                               N_Vector yP, N_Vector ypP, N_Vector fctP,
                               N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  IDADlsMem idadls_mem;
  int ier, retval;
  int intn;

  idadls_mem = (IDADlsMem) lmem;

  intn = (int) n;

  /* Call Jacobian function */
  nje++;
  SetToZero(JJ);
  retval = djac(n, tn, cj, yP, ypP, fctP, JJ, J_data, tmp1, tmp2, tmp3);
  if (retval < 0) {
    IDAProcessError(IDA_mem, IDADLS_JACFUNC_UNRECVR, "IDASLAPACK", "idaLapackDenseSetup", MSGD_JACFUNC_FAILED);
    last_flag = IDADLS_JACFUNC_UNRECVR;
    return(-1);
  } else if (retval > 0) {
    last_flag = IDADLS_JACFUNC_RECVR;
    return(1);
  }
  
  /* Do LU factorization of M */
  dgetrf_f77(&intn, &intn, JJ->data, &intn, pivots, &ier);

  /* Return 0 if the LU was complete; otherwise return 1 */
  last_flag = (long int) ier;
  if (ier > 0) return(1);
  return(0);
}

/*
 * idaLapackDenseSolve handles the solve operation for the dense linear solver
 * by calling the dense backsolve routine.
 */
static int idaLapackDenseSolve(IDAMem IDA_mem, N_Vector b, N_Vector weight,
                               N_Vector yC, N_Vector ypC, N_Vector fctC)
{
  IDADlsMem idadls_mem;
  realtype *bd, fact;
  int ier, one = 1;
  int intn;

  idadls_mem = (IDADlsMem) lmem;

  intn = (int) n;

  bd = N_VGetArrayPointer(b);

  dgetrs_f77("N", &intn, &one, JJ->data, &intn, pivots, bd, &intn, &ier, 1); 
  if (ier > 0) return(1);

  /* Scale the correction to account for change in cj. */
  if (cjratio != ONE) {
    fact = TWO/(ONE + cjratio);
    dscal_f77(&intn, &fact, bd, &one); 
  }

  last_flag = IDADLS_SUCCESS;
  return(0);
}

/*
 * idaLapackDenseFree frees memory specific to the dense linear solver.
 */
static int idaLapackDenseFree(IDAMem IDA_mem)
{
  IDADlsMem  idadls_mem;

  idadls_mem = (IDADlsMem) lmem;
  
  DestroyMat(JJ);
  DestroyArray(pivots);
  free(idadls_mem); 
  idadls_mem = NULL;

  return(0);
}

/* 
 * =================================================================
 *  PRIVATE FUNCTIONS FOR IMPLICIT INTEGRATION WITH BAND JACOBIANS
 * =================================================================
 */

/*
 * idaLapackBandInit does remaining initializations specific to the band
 * linear solver.
 */
static int idaLapackBandInit(IDAMem IDA_mem)
{
  IDADlsMem idadls_mem;

  idadls_mem = (IDADlsMem) lmem;

  nje   = 0;
  nreDQ = 0;

  if (jacDQ) {
    bjac = idaDlsBandDQJac;
    J_data = IDA_mem;
  } else {
    J_data = IDA_mem->ida_user_data;
  }

  last_flag = IDADLS_SUCCESS;
  return(0);
}

/*
 * idaLapackBandSetup does the setup operations for the band linear solver.
 * It calls the Jacobian function to obtain the Newton matrix M = F_y + c_j*F_y', 
 * updates counters, and calls the band LU factorization routine.
 */
static int idaLapackBandSetup(IDAMem IDA_mem,
                              N_Vector yP, N_Vector ypP, N_Vector fctP, 
                              N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  IDADlsMem idadls_mem;
  int ier, retval;
  int intn, iml, imu, ldmat;

  idadls_mem = (IDADlsMem) lmem;

  intn = (int) n;
  iml = (int) ml;
  imu = (int) mu;
  ldmat = JJ->ldim;

  /* Call Jacobian function */
  nje++;
  SetToZero(JJ);
  retval = bjac(n, mu, ml, tn, cj, yP, ypP, fctP, JJ, J_data, tmp1, tmp2, tmp3);
  if (retval < 0) {
    IDAProcessError(IDA_mem, IDADLS_JACFUNC_UNRECVR, "IDASLAPACK", "idaLapackBandSetup", MSGD_JACFUNC_FAILED);
    last_flag = IDADLS_JACFUNC_UNRECVR;
    return(-1);
  } else if (retval > 0) {
    last_flag = IDADLS_JACFUNC_RECVR;
    return(+1);
  }
  
  /* Do LU factorization of M */
  dgbtrf_f77(&intn, &intn, &iml, &imu, JJ->data, &ldmat, pivots, &ier);

  /* Return 0 if the LU was complete; otherwise return 1 */
  last_flag = (long int) ier;
  if (ier > 0) return(1);
  return(0);

}

/*
 * idaLapackBandSolve handles the solve operation for the band linear solver
 * by calling the band backsolve routine.
 */
static int idaLapackBandSolve(IDAMem IDA_mem, N_Vector b, N_Vector weight,
                              N_Vector yC, N_Vector ypC, N_Vector fctC)
{
  IDADlsMem idadls_mem;
  realtype *bd, fact;
  int ier, one = 1;
  int intn, iml, imu, ldmat;

  idadls_mem = (IDADlsMem) lmem;

  intn = (int) n;
  iml = (int) ml;
  imu = (int) mu;
  ldmat = JJ->ldim;

  bd = N_VGetArrayPointer(b);

  dgbtrs_f77("N", &intn, &iml, &imu, &one, JJ->data, &ldmat, pivots, bd, &intn, &ier, 1);
  if (ier > 0) return(1);

  /* For BDF, scale the correction to account for change in cj */
  if (cjratio != ONE) {
    fact = TWO/(ONE + cjratio);
    dscal_f77(&intn, &fact, bd, &one); 
  }

  last_flag = IDADLS_SUCCESS;
  return(0);
}

/*
 * idaLapackBandFree frees memory specific to the band linear solver.
 */
static int idaLapackBandFree(IDAMem IDA_mem)
{
  IDADlsMem  idadls_mem;

  idadls_mem = (IDADlsMem) lmem;
  
  DestroyMat(JJ);
  DestroyArray(pivots);
  free(idadls_mem); 
  idadls_mem = NULL;

  return(0);
}


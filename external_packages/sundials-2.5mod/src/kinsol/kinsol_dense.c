/*
 * -----------------------------------------------------------------
 * $Revision: 1.11 $
 * $Date: 2010/12/01 22:43:33 $
 * ----------------------------------------------------------------- 
 * Programmer(s): Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2002, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see the LICENSE file.
 * -----------------------------------------------------------------
 * This is the implementation file for the KINDENSE linear solver.
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>

#include <kinsol/kinsol_dense.h>
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

/* KINDENSE linit, lsetup, lsolve, and lfree routines */ 
static int kinDenseInit(KINMem kin_mem);
static int kinDenseSetup(KINMem kin_mem);
static int kinDenseSolve(KINMem kin_mem, N_Vector x, N_Vector b,
                         realtype *res_norm);
static void kinDenseFree(KINMem kin_mem);

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
 * KINDense
 * -----------------------------------------------------------------
 * This routine initializes the memory record and sets various function
 * fields specific to the dense linear solver module. 
 * KINDense sets the kin_linit, kin_lsetup, kin_lsolve, kin_lfree fields 
 * in *kinmem to be kinDenseInit, kinDenseSetup, kinDenseSolve, and 
 * kinDenseFree, respectively.  
 * It allocates memory for a structure of type KINDlsMemRec and sets 
 * the kin_lmem field in *kinmem to the address of this structure.  
 * It sets setupNonNull in *kinmem to TRUE, and the djac field to the 
 * default kinDlsDenseDQJac.
 * Finally, it allocates memory for J and lpivots.
 *
 * NOTE: The dense linear solver assumes a serial implementation
 *       of the NVECTOR package. Therefore, KINDense will first 
 *       test for compatible a compatible N_Vector internal
 *       representation by checking that N_VGetArrayPointer and
 *       N_VSetArrayPointer exist.
 * -----------------------------------------------------------------
 */

int KINDense(void *kinmem, long int N)
{
  KINMem kin_mem;
  KINDlsMem kindls_mem;

  /* Return immediately if kinmem is NULL */
  if (kinmem == NULL) {
    KINProcessError(NULL, KINDLS_MEM_NULL, "KINDENSE", "KINDense", MSGD_KINMEM_NULL);
    return(KINDLS_MEM_NULL);
  }
  kin_mem = (KINMem) kinmem;

  /* Test if the NVECTOR package is compatible with the DENSE solver */
  if (vec_tmpl->ops->nvgetarraypointer == NULL ||
      vec_tmpl->ops->nvsetarraypointer == NULL) {
    KINProcessError(kin_mem, KINDLS_ILL_INPUT, "KINDENSE", "KINDense", MSGD_BAD_NVECTOR);
    return(KINDLS_ILL_INPUT);
  }

  if (lfree !=NULL) lfree(kin_mem);

  /* Set four main function fields in kin_mem */
  linit  = kinDenseInit;
  lsetup = kinDenseSetup;
  lsolve = kinDenseSolve;
  lfree  = kinDenseFree;

  /* Get memory for KINDlsMemRec */
  kindls_mem = NULL;
  kindls_mem = (KINDlsMem) malloc(sizeof(struct KINDlsMemRec));
  if (kindls_mem == NULL) {
    KINProcessError(kin_mem, KINDLS_MEM_FAIL, "KINDENSE", "KINDense", MSGD_MEM_FAIL);
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
  n = N;

  /* Allocate memory for J and pivot array */
  
  J = NULL;
  J = NewDenseMat(N, N);
  if (J == NULL) {
    KINProcessError(kin_mem, KINDLS_MEM_FAIL, "KINDENSE", "KINDense", MSGD_MEM_FAIL);
    free(kindls_mem); kindls_mem = NULL;
    return(KINDLS_MEM_FAIL);
  }

  lpivots = NULL;
  lpivots = NewLintArray(N);
  if (lpivots == NULL) {
    KINProcessError(kin_mem, KINDLS_MEM_FAIL, "KINDENSE", "KINDense", MSGD_MEM_FAIL);
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
 * kinDenseInit
 * -----------------------------------------------------------------
 * This routine does remaining initializations specific to the dense
 * linear solver.
 * -----------------------------------------------------------------
 */

static int kinDenseInit(KINMem kin_mem)
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
 * kinDenseSetup
 * -----------------------------------------------------------------
 * This routine does the setup operations for the dense linear solver.
 * It calls the dense LU factorization routine.
 * -----------------------------------------------------------------
 */

static int kinDenseSetup(KINMem kin_mem)
{
  KINDlsMem kindls_mem;
  int retval;
  long int ier;

  kindls_mem = (KINDlsMem) lmem;
 
  nje++;
  SetToZero(J); 
  retval = djac(n, uu, fval, J, J_data, vtemp1, vtemp2);
  if (retval != 0) {
    last_flag = -1;
    return(-1);
  }

  /* Do LU factorization of J */
  ier = DenseGETRF(J, lpivots); 

  /* Return 0 if the LU was complete; otherwise return -1 */
  last_flag = ier;
  if (ier > 0) return(-1);

  return(0);
}

/*
 * -----------------------------------------------------------------
 * kinDenseSolve
 * -----------------------------------------------------------------
 * This routine handles the solve operation for the dense linear solver
 * by calling the dense backsolve routine.  The returned value is 0.
 * -----------------------------------------------------------------
 */

static int kinDenseSolve(KINMem kin_mem, N_Vector x, N_Vector b, realtype *res_norm)
{
  KINDlsMem kindls_mem;
  realtype *xd;

  kindls_mem = (KINDlsMem) lmem;

  /* Copy the right-hand side into x */

  N_VScale(ONE, b, x);
  
  xd = N_VGetArrayPointer(x);

  /* Back-solve and get solution in x */
  
  DenseGETRS(J, lpivots, xd);

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
 * kinDenseFree
 * -----------------------------------------------------------------
 * This routine frees memory specific to the dense linear solver.
 * -----------------------------------------------------------------
 */

static void kinDenseFree(KINMem kin_mem)
{
  KINDlsMem  kindls_mem;

  kindls_mem = (KINDlsMem) lmem;
  
  DestroyMat(J);
  DestroyArray(lpivots);
  free(kindls_mem); kindls_mem = NULL;
}


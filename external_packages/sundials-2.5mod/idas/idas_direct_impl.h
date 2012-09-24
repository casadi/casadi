/*
 * -----------------------------------------------------------------
 * $Revision: 1.7 $
 * $Date: 2010/12/01 22:39:18 $
 * ----------------------------------------------------------------- 
 * Programmer: Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2006, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see the LICENSE file.
 * -----------------------------------------------------------------
 * Implementation header file for the IDADLS linear solvers.
 * -----------------------------------------------------------------
 */

#ifndef _IDASDLS_IMPL_H
#define _IDASDLS_IMPL_H

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

#include <idas/idas_direct.h>

/*
 * =================================================================
 * I D A S D I R E C T    I N T E R N A L    C O N S T A N T S
 * =================================================================
 */

/*
 * =================================================================
 * PART I:  F O R W A R D    P R O B L E M S
 * =================================================================
 */

/*
 * -----------------------------------------------------------------
 * Types : IDADlsMemRec, IDADlsMem                             
 * -----------------------------------------------------------------
 * IDADlsMem is pointer to a IDADlsMemRec structure.
 * -----------------------------------------------------------------
 */

typedef struct IDADlsMemRec {

  int d_type;               /* Type of Jacobians (DENSE or BAND)             */

  long int d_n;             /* problem dimension                             */

  long int d_ml;            /* b_ml = lower bandwidth of savedJ              */
  long int d_mu;            /* b_mu = upper bandwidth of savedJ              */ 
  long int d_smu;           /* upper bandwith of M = MIN(N-1,b_mu+b_ml)      */

  booleantype d_jacDQ;      /* TRUE if using internal DQ Jacobian approx.    */
  IDADlsDenseJacFn d_djac;  /* dense Jacobian routine to be called           */
  IDADlsBandJacFn d_bjac;   /* band Jacobian routine to be called            */
  void *d_J_data;           /* J_data is passed to djac or bjac              */

  DlsMat d_J;               /* J = dF/dy + cj*dF/dy'                         */

  int *d_pivots;            /* pivots = int pivot array for PM = LU          */
  long int *d_lpivots;      /* lpivots = long int pivot array for PM = LU    */
  
  long int d_nje;           /* nje = no. of calls to jac                     */

  long int d_nreDQ;         /* no. of calls to res due to DQ Jacobian approx.*/

  long int d_last_flag;     /* last error return flag                        */
  
} *IDADlsMem;

/*
 * -----------------------------------------------------------------
 * Prototypes of internal functions
 * -----------------------------------------------------------------
 */
  
int idaDlsDenseDQJac(long int N, realtype tt, realtype c_j,
		     N_Vector yy, N_Vector yp, N_Vector rr, 
		     DlsMat Jac, void *data,
		     N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
  
int idaDlsBandDQJac(long int N, long int mupper, long int mlower,
		    realtype tt, realtype c_j, 
		    N_Vector yy, N_Vector yp, N_Vector rr,
		    DlsMat Jac, void *data,
		    N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

/*
 * =================================================================
 * PART II:  B A C K W A R D    P R O B L E M S
 * =================================================================
 */

/*
 * -----------------------------------------------------------------
 * Types : IDADlsMemRecB, IDADlsMemB       
 * -----------------------------------------------------------------
 * An IDADLS linear solver's specification function attaches such
 * a structure to the lmemB filed of IDABMem
 * -----------------------------------------------------------------
 */

typedef struct IDADlsMemRecB {

  int d_typeB;

  IDADlsDenseJacFnB d_djacB;
  IDADlsBandJacFnB  d_bjacB;

} *IDADlsMemB;


/*
 * =================================================================
 * E R R O R   M E S S A G E S
 * =================================================================
 */

#define MSGD_IDAMEM_NULL "Integrator memory is NULL."
#define MSGD_BAD_NVECTOR "A required vector operation is not implemented."
#define MSGD_BAD_SIZES "Illegal bandwidth parameter(s). Must have 0 <=  ml, mu <= N-1."
#define MSGD_MEM_FAIL "A memory request failed."
#define MSGD_LMEM_NULL "Linear solver memory is NULL."
#define MSGD_JACFUNC_FAILED "The Jacobian routine failed in an unrecoverable manner."

#define MSGD_CAMEM_NULL "idaadj_mem = NULL illegal."
#define MSGD_LMEMB_NULL "Linear solver memory is NULL for the backward integration."
#define MSGD_BAD_T "Bad t for interpolation."
#define MSGD_BAD_WHICH "Illegal value for which."
#define MSGD_NO_ADJ "Illegal attempt to call before calling IDAAdjInit."

#ifdef __cplusplus
}
#endif

#endif

/*
 * -----------------------------------------------------------------
 * $Revision: 1.8 $
 * $Date: 2010/12/01 22:30:42 $
 * ----------------------------------------------------------------- 
 * Programmer(s): Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2005, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see the LICENSE file.
 * -----------------------------------------------------------------
 * Implementation header file for the CVBBDPRE module.
 * -----------------------------------------------------------------
 */

#ifndef _CVSBBDPRE_IMPL_H
#define _CVSBBDPRE_IMPL_H

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

#include <cvodes/cvodes_bbdpre.h>
#include <sundials/sundials_band.h>

/*
 * -----------------------------------------------------------------
 * Type: CVBBDPrecData
 * -----------------------------------------------------------------
 */

typedef struct CVBBDPrecDataRec {

  /* passed by user to CVBBDPrecAlloc and used by PrecSetup/PrecSolve */

  long int mudq, mldq, mukeep, mlkeep;
  realtype dqrely;
  CVLocalFn gloc;
  CVCommFn cfn;

  /* set by CVBBDPrecSetup and used by CVBBDPrecSolve */

  DlsMat savedJ;
  DlsMat savedP;
  long int *lpivots;

  /* set by CVBBDPrecAlloc and used by CVBBDPrecSetup */

  long int n_local;

  /* available for optional output */

  long int rpwsize;
  long int ipwsize;
  long int nge;

  /* pointer to cvode_mem */

  void *cvode_mem;

} *CVBBDPrecData;


/*
 * -----------------------------------------------------------------
 * Type: CVBBDPrecDataB
 * -----------------------------------------------------------------
 */

typedef struct CVBBDPrecDataRecB {

  /* BBD user functions (glocB and cfnB) for backward run */
  CVLocalFnB glocB;
  CVCommFnB  cfnB;

} *CVBBDPrecDataB;

/*
 * -----------------------------------------------------------------
 * CVBBDPRE error messages
 * -----------------------------------------------------------------
 */

#define MSGBBD_MEM_NULL    "Integrator memory is NULL."
#define MSGBBD_LMEM_NULL   "Linear solver memory is NULL. One of the SPILS linear solvers must be attached."
#define MSGBBD_MEM_FAIL    "A memory request failed."
#define MSGBBD_BAD_NVECTOR "A required vector operation is not implemented."
#define MSGBBD_PMEM_NULL   "BBD peconditioner memory is NULL. CVBBDPrecInit must be called."
#define MSGBBD_FUNC_FAILED "The gloc or cfn routine failed in an unrecoverable manner."

#define MSGBBD_NO_ADJ      "Illegal attempt to call before calling CVodeAdjInit."
#define MSGBBD_BAD_WHICH   "Illegal value for the which parameter."
#define MSGBBD_PDATAB_NULL "BBD preconditioner memory is NULL for the backward integration."
#define MSGBBD_BAD_TINTERP "Bad t for interpolation."


#ifdef __cplusplus
}
#endif

#endif

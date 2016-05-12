/*
 * -----------------------------------------------------------------
 * $Revision: 1.10 $
 * $Date: 2010/12/01 22:30:42 $
 * ----------------------------------------------------------------- 
 * Programmer(s): Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2005, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see the LICENSE file.
 * -----------------------------------------------------------------
 * Implementation header file for the CVBANDPRE module.
 * -----------------------------------------------------------------
 */

#ifndef _CVSBANDPRE_IMPL_H
#define _CVSBANDPRE_IMPL_H

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

#include <cvodes/cvodes_bandpre.h>
#include <sundials/sundials_band.h>
#include <sundials/sundials_direct.h>

/*
 * -----------------------------------------------------------------
 * Type: CVBandPrecData
 * -----------------------------------------------------------------
 */

typedef struct CVBandPrecDataRec {

  /* Data set by user in CVBandPrecInit */

  long int N;
  long int ml, mu;

  /* Data set by CVBandPrecSetup */

  DlsMat savedJ;
  DlsMat savedP;
  long int *lpivots;

  /* Rhs calls */

  long int nfeBP;

  /* Pointer to cvode_mem */

  void *cvode_mem;

} *CVBandPrecData;

/*
 * -----------------------------------------------------------------
 * CVBANDPRE error messages
 * -----------------------------------------------------------------
 */

#define MSGBP_MEM_NULL       "Integrator memory is NULL."
#define MSGBP_LMEM_NULL      "Linear solver memory is NULL. One of the SPILS linear solvers must be attached."
#define MSGBP_MEM_FAIL       "A memory request failed."
#define MSGBP_BAD_NVECTOR    "A required vector operation is not implemented."
#define MSGBP_PMEM_NULL      "Band preconditioner memory is NULL. CVBandPrecInit must be called."
#define MSGBP_RHSFUNC_FAILED "The right-hand side routine failed in an unrecoverable manner."

#define MSGBP_NO_ADJ         "Illegal attempt to call before calling CVodeAdjInit."
#define MSGBP_BAD_WHICH      "Illegal value for parameter which."

#ifdef __cplusplus
}
#endif

#endif

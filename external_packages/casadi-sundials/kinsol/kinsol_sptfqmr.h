/*
 * -----------------------------------------------------------------
 * $Revision: 1.2 $
 * $Date: 2006/11/29 00:05:07 $
 * -----------------------------------------------------------------
 * Programmer(s): Aaron Collier @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2005, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see the LICENSE file.
 * -----------------------------------------------------------------
 * This is the public header file for the KINSOL scaled preconditioned
 * TFQMR linear solver module, KINSPTFQMR.
 * -----------------------------------------------------------------
 */

#ifndef _KINSPTFQMR_H
#define _KINSPTFQMR_H

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

#include <kinsol/kinsol_spils.h>
#include <sundials/sundials_sptfqmr.h>

/*
 * -----------------------------------------------------------------
 * Function : KINSptfqmr
 * -----------------------------------------------------------------
 * KINSptfqmr links the main KINSOL solver module with the SPTFQMR
 * linear solver module. The routine establishes the inter-module
 * interface by setting the generic KINSOL pointers linit, lsetup,
 * lsolve, and lfree to KINSptfqmrInit, KINSptfqmrSetup, KINSptfqmrSolve,
 * and KINSptfqmrFree, respectively.
 *
 *  kinmem  pointer to an internal memory block allocated during a
 *          prior call to KINCreate
 *
 *  maxl  maximum allowable dimension of Krylov subspace (passing
 *        a value of 0 (zero) will cause the default value
 *        KINSPTFQMR_MAXL (predefined constant) to be used)
 *
 * If successful, KINSptfqmr returns KINSPTFQMR_SUCCESS. If an error
 * occurs, then KINSptfqmr returns an error code (negative integer
 * value).
 *
 * -----------------------------------------------------------------
 * KINSptfqmr Return Values
 * -----------------------------------------------------------------
 * The possible return values for the KINSptfqmr subroutine are the
 * following:
 *
 * KINSPTFQMR_SUCCESS : means the KINSPTFQMR linear solver module
 *                      (implementation of the TFQMR method) was
 *                      successfully initialized - allocated system
 *                      memory and set shared variables to default
 *                      values [0]
 *
 * KINSPTFQMR_MEM_NULL : means a NULL KINSOL memory block pointer
 *                       was given (must call the KINCreate and
 *                       KINMalloc memory allocation subroutines
 *                       prior to calling KINSptfqmr) [-1]
 *
 * KINSPTFQMR_MEM_FAIL : means either insufficient system resources
 *                       were available to allocate memory for the
 *                       main KINSPTFQMR data structure (type
 *                       KINSptfqmrMemRec), or the SptfqmrMalloc
 *                       subroutine failed (unable to allocate enough
 *                       system memory for vector storate and/or the
 *                       main SPTFQMR data structure
 *                       (type SptfqmrMemRec)) [-4]
 *
 * KINSPTFQMR_ILL_INPUT : means either a supplied parameter was invalid,
 *                        or the NVECTOR implementation is NOT
 *                        compatible [-3]
 *
 * The above constants are defined in kinsol_spils.h
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT int KINSptfqmr(void *kinmem, int maxl);


#ifdef __cplusplus
}
#endif

#endif

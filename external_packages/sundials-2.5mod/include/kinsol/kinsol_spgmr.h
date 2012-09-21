/*
 * -----------------------------------------------------------------
 * $Revision: 1.2 $
 * $Date: 2006/11/29 00:05:07 $
 * -----------------------------------------------------------------
 * Programmer(s): Allan Taylor, Alan Hindmarsh, Radu Serban, and
 *                Aaron Collier @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2002, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see the LICENSE file.
 * -----------------------------------------------------------------
 * This is the header file for the KINSOL Scaled Preconditioned GMRES     
 * linear solver module, KINSPGMR.    
 * -----------------------------------------------------------------
 */

#ifndef _KINSPGMR_H
#define _KINSPGMR_H

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

#include <kinsol/kinsol_spils.h>
#include <sundials/sundials_spgmr.h>

/*
 * -----------------------------------------------------------------
 * Function : KINSpgmr
 * -----------------------------------------------------------------
 * KINSpgmr links the main KINSOL solver module with the SPGMR
 * linear solver module. The routine establishes the inter-module
 * interface by setting the generic KINSOL pointers linit, lsetup,
 * lsolve, and lfree to KINSpgmrInit, KINSpgmrSetup, KINSpgmrSolve,
 * and KINSpgmrFree, respectively.
 *
 *  kinmem  pointer to an internal memory block allocated during a
 *          prior call to KINCreate
 *
 *  maxl  maximum allowable dimension of Krylov subspace (passing
 *        a value of 0 (zero) will cause the default value
 *        KINSPILS_MAXL (predefined constant) to be used)
 *
 * -----------------------------------------------------------------
 * KINSpgmr Return Values
 * -----------------------------------------------------------------
 *
 * The possible return values for the KINSpgmr subroutine are the
 * following:
 *
 * KINSPILS_SUCCESS : means the KINSPGMR linear solver module
 *                    (implementation of the GMRES method) was
 *                    successfully initialized - allocated system
 *                    memory and set shared variables to default
 *                    values [0]
 *
 * KINSPILS_MEM_NULL : means a NULL KINSOL memory block pointer was
 *                     given (must call the KINCreate and KINMalloc
 *                     memory allocation subroutines prior to
 *                     calling KINSpgmr) [-1]
 *
 * KINSPILS_MEM_FAIL : means either insufficient system resources
 *                     were available to allocate memory for the main
 *                     KINSPGMR data structure (type KINSpgmrMemRec),
 *                     or the SpgmrMalloc subroutine failed (unable
 *                     to allocate enough system memory for vector
 *                     storage and/or the main SPGMR data structure
 *                     (type SpgmrMemRec)) [-4]
 *
 * KINSPILS_ILL_INPUT : means a supplied parameter was invalid
 *                      (check error message) [-3]
 *
 * The above constants are defined in kinsol_spils.h
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT int KINSpgmr(void *kinmem, int maxl);


#ifdef __cplusplus
}
#endif

#endif

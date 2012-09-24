/*
 * -----------------------------------------------------------------
 * $Revision: 1.5 $
 * $Date: 2010/12/01 22:16:17 $
 * -----------------------------------------------------------------
 * Programmer(s): Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2002, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see the LICENSE file.
 * -----------------------------------------------------------------
 * This is the header file for the KINSOL dense linear solver module, 
 * KINDENSE.
 * -----------------------------------------------------------------
 */

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

#ifndef _KINDENSE_H
#define _KINDENSE_H

#include <kinsol/kinsol_direct.h>
#include <sundials/sundials_dense.h>

/*
 * -----------------------------------------------------------------
 * Function : KINDense
 * -----------------------------------------------------------------
 * A call to the KINDense function links the main solver with the
 * KINDENSE linear solver. Its arguments are as follows:
 *
 * kinmem - pointer to an internal memory block allocated during a
 *          prior call to KINCreate
 *
 * N      - problem size
 *
 * The return value of KINDense is one of:
 *    KINDLS_SUCCESS   if successful
 *    KINDLS_MEM_NULL  if the kinsol memory was NULL
 *    KINDLS_MEM_FAIL  if there was a memory allocation failure
 *    KINDLS_ILL_INPUT if a required vector operation is missing
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT int KINDense(void *kinmem, long int N);

#endif

#ifdef __cplusplus
}
#endif

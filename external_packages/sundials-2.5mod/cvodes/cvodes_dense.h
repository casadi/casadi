/*
 * -----------------------------------------------------------------
 * $Revision: 1.8 $
 * $Date: 2010/12/01 22:13:10 $
 * ----------------------------------------------------------------- 
 * Programmer(s): Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2005, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see the LICENSE file.
 * -----------------------------------------------------------------
 * This is the header file for the dense linear solver CVSDENSE.
 * -----------------------------------------------------------------
 */

#ifndef _CVSDENSE_H
#define _CVSDENSE_H

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

#include <cvodes/cvodes_direct.h>
#include <sundials/sundials_dense.h>

/*
 * -----------------------------------------------------------------
 * Function: CVDense
 * -----------------------------------------------------------------
 * A call to the CVDense function links the main integrator with
 * the CVSDENSE linear solver.
 *
 * cvode_mem is the pointer to the integrator memory returned by
 *           CVodeCreate.
 *
 * N is the size of the ODE system.
 *
 * The return value of CVDense is one of:
 *    CVDLS_SUCCESS   if successful
 *    CVDLS_MEM_NULL  if the cvode memory was NULL
 *    CVDLS_MEM_FAIL  if there was a memory allocation failure
 *    CVDLS_ILL_INPUT if a required vector operation is missing
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT int CVDense(void *cvode_mem, long int N);

/*
 * -----------------------------------------------------------------
 * Function: CVDenseB
 * -----------------------------------------------------------------
 * CVDenseB links the main CVODE integrator with the CVSDENSE
 * linear solver for the backward integration.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT int CVDenseB(void *cvode_mem, int which, long int nB);

#ifdef __cplusplus
}
#endif

#endif

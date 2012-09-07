/*
 * -----------------------------------------------------------------
 * $Revision: 1.6 $
 * $Date: 2008/04/18 19:42:37 $
 * ----------------------------------------------------------------- 
 * Programmer: Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2006, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see the LICENSE file.
 * -----------------------------------------------------------------
 * Header file for the CVODES dense linear solver CVSLAPACK.
 * -----------------------------------------------------------------
 */

#ifndef _CVSLAPACK_H
#define _CVSLAPACK_H

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

#include <cvodes/cvodes_direct.h>
#include <sundials/sundials_lapack.h>

/*
 * -----------------------------------------------------------------
 * Function: CVLapackDense
 * -----------------------------------------------------------------
 * A call to the CVLapackDense function links the main integrator
 * with the CVSLAPACK linear solver using dense Jacobians.
 *
 * cvode_mem is the pointer to the integrator memory returned by
 *           CVodeCreate.
 *
 * N is the size of the ODE system.
 *
 * The return value of CVLapackDense is one of:
 *    CVDLS_SUCCESS   if successful
 *    CVDLS_MEM_NULL  if the CVODES memory was NULL
 *    CVDLS_MEM_FAIL  if there was a memory allocation failure
 *    CVDLS_ILL_INPUT if a required vector operation is missing
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT int CVLapackDense(void *cvode_mem, int N);

/*
 * -----------------------------------------------------------------
 * Function: CVLapackBand
 * -----------------------------------------------------------------
 * A call to the CVLapackBand function links the main integrator
 * with the CVSLAPACK linear solver using banded Jacobians. 
 *
 * cvode_mem is the pointer to the integrator memory returned by
 *           CVodeCreate.
 *
 * N is the size of the ODE system.
 *
 * mupper is the upper bandwidth of the band Jacobian approximation.
 *
 * mlower is the lower bandwidth of the band Jacobian approximation.
 *
 * The return value of CVLapackBand is one of:
 *    CVDLS_SUCCESS   if successful
 *    CVDLS_MEM_NULL  if the CVODES memory was NULL
 *    CVDLS_MEM_FAIL  if there was a memory allocation failure
 *    CVDLS_ILL_INPUT if a required vector operation is missing or
 *                       if a bandwidth has an illegal value.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT int CVLapackBand(void *cvode_mem, int N, int mupper, int mlower);

/*
 * -----------------------------------------------------------------
 * Function: CVLapackDenseB
 * -----------------------------------------------------------------
 * CVLapackDenseB links the main CVODE integrator with the dense
 * CVSLAPACK linear solver for the backward integration.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT int CVLapackDenseB(void *cvode_mem, int which, int nB);

/*
 * -----------------------------------------------------------------
 * Function: CVLapackBandB
 * -----------------------------------------------------------------
 * CVLapackBandB links the main CVODE integrator with the band
 * CVSLAPACK linear solver for the backward integration.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT int CVLapackBandB(void *cvode_mem, int which,
                                  int nB, int mupperB, int mlowerB);

#ifdef __cplusplus
}
#endif

#endif

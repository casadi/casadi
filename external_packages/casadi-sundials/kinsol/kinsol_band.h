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
 * This is the header file for the KINSOL band linear solver, KINBAND.
 * -----------------------------------------------------------------
 */

#ifndef _KINBAND_H
#define _KINBAND_H

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

#include <kinsol/kinsol_direct.h>
#include <sundials/sundials_band.h>

/*
 * -----------------------------------------------------------------
 * Function : KINBand
 * -----------------------------------------------------------------
 * A call to the KINBand function links the main solver with the 
 * KINBAND linear solver. Its arguments are as follows:
 *
 * kinmem - pointer to the integrator memory returned by KINCreate.
 *
 * N      - problem size
 *
 * mupper - upper bandwidth of the band Jacobian
 *
 * mlower - lower bandwidth of the band Jacobian
 *
 * The return value of KINBand is one of:
 *    KINDLS_SUCCESS   if successful
 *    KINDLS_MEM_NULL  if the kinsol memory was NULL
 *    KINDLS_MEM_FAIL  if there was a memory allocation failure
 *    KINDLS_ILL_INPUT if a required vector operation is missing
 *                        or if a bandwidth has an illegal value.
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT int KINBand(void *kinmem, long int N, long int mupper, long int mlower);

#ifdef __cplusplus
}
#endif

#endif

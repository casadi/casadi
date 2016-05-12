/*
 * -----------------------------------------------------------------
 * $Revision: 1.4 $
 * $Date: 2007/03/22 18:05:50 $
 * ----------------------------------------------------------------- 
 * Programmer(s): Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2005, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see the LICENSE file.
 * -----------------------------------------------------------------
 * This is the header file for the CVODES scaled preconditioned 
 * GMRES linear solver, CVSPGMR.
 *
 * Part I contains function prototypes for using CVSPGMR on forward 
 * problems (IVP integration and/or FSA)
 *
 * Part II contains function prototypes for using CVSPGMR on adjoint 
 * (backward) problems
 * -----------------------------------------------------------------
 */

#ifndef _CVSSPGMR_H
#define _CVSSPGMR_H

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

#include <cvodes/cvodes_spils.h>
#include <sundials/sundials_spgmr.h>

/* 
 * -----------------------------------------------------------------
 * PART I - forward problems
 * -----------------------------------------------------------------
 */

/*
 * -----------------------------------------------------------------
 * Function : CVSpgmr
 * -----------------------------------------------------------------
 * A call to the CVSpgmr function links the main CVODE integrator
 * with the CVSPGMR linear solver.
 *
 * cvode_mem is the pointer to the integrator memory returned by
 *           CVodeCreate.
 *
 * pretype   is the type of user preconditioning to be done.
 *           This must be one of the four enumeration constants
 *           NONE, LEFT, RIGHT, or BOTH defined in iterative.h.
 *           These correspond to no preconditioning,
 *           left preconditioning only, right preconditioning
 *           only, and both left and right preconditioning,
 *           respectively.
 *
 * maxl      is the maximum Krylov dimension. This is an
 *           optional input to the CVSPGMR solver. Pass 0 to
 *           use the default value CVSPILS_MAXL=5.
 *
 * The return value of CVSpgmr is one of:
 *    CVSPILS_SUCCESS   if successful
 *    CVSPILS_MEM_NULL  if the cvode memory was NULL
 *    CVSPILS_MEM_FAIL  if there was a memory allocation failure
 *    CVSPILS_ILL_INPUT if a required vector operation is missing
 * The above constants are defined in cvodes_spils.h
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT int CVSpgmr(void *cvode_mem, int pretype, int maxl);

/* 
 * -----------------------------------------------------------------
 * PART II - backward problems
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT int CVSpgmrB(void *cvode_mem, int which,
                             int pretypeB, int maxlB);
  


#ifdef __cplusplus
}
#endif

#endif

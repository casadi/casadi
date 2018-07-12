#include <math.h>
#include <float.h>
/* 
 * This file contains properties configurable by the user
 */
#ifndef GLOBALS_H
#define GLOBALS_H
    /*
     * ---------------------------------
     * constants used with double data type
     * ---------------------------------
     */
    #define real_t double
    /* data types have different absolute value functions */ 
    #define ABS(x) fabs(x)
    #define MAX(x,y) fmax(x,y)
    #define MIN(x,y) fmin(x,y)

    /* Machine accuracy of IEEE double */ 
    #define MACHINE_ACCURACY DBL_EPSILON
    /* large number use with things like indicator functions */ 
    #define LARGE 10000000000

    /* return values for failure and success of function, the unix way */
    #define FAILURE 1
    #define SUCCESS 0

    /* define the 2 boolean states */
    #define TRUE 1
    #define FALSE 0

    /* minimum amount of steps Panoc always should execute */
    #ifndef PANOC_MIN_STEPS
        #define PANOC_MIN_STEPS 0
    #endif

    /* 
    * ---------------------------------
    * Proximal gradient descent definitions
    * ---------------------------------
    */

    /* constant delta used to estimate lipschitz constant  */
    #ifndef PROXIMAL_GRAD_DESC_SAFETY_VALUE
        #define PROXIMAL_GRAD_DESC_SAFETY_VALUE 0.05
    #endif

    /* ---------------------------------
    * lipschitz etimator definitions
    * ---------------------------------
    */
    #ifndef DELTA_LIPSCHITZ
        #define DELTA_LIPSCHITZ (1e-12)
    #endif
    #ifndef DELTA_LIPSCHITZ_SAFETY_VALUE
        #define DELTA_LIPSCHITZ_SAFETY_VALUE (1e-6)
    #endif

    #ifndef FBE_LINESEARCH_MAX_ITERATIONS
        #define FBE_LINESEARCH_MAX_ITERATIONS 10
    #endif

    #ifndef LBGFS_SAFETY_SMALL_VALUE
        #define LBGFS_SAFETY_SMALL_VALUE (10e-12)
    #endif

#endif
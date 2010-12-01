#ifndef INTEGRATOR_C_H 
#define INTEGRATOR_C_H 

#ifdef __cplusplus
  extern "C" {
#endif 

#include "fx_c.h"

/* Print to cerr: redirect to print to a file */
int casadi_integrator_integrate(fx_ref ref, double t_out);

/* Reset the solver and bring the time back to t0 */
int casadi_integrator_reset(fx_ref ref, int with_sens);

#ifdef __cplusplus
  }
#endif 

#endif /* INTEGRATOR_C_H  */

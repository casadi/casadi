#ifndef SIMULATOR_C_H 
#define SIMULATOR_C_H 

#ifdef __cplusplus
  extern "C" {
#endif 

#include "integrator_c.h"

/* Create a function for numerical evaluation */
int casadi_simulator(fx_ref fcn, fx_ref integrator, fx_ref output_fcn, int ngrid);

#ifdef __cplusplus
  }
#endif 

#endif /* SIMULATOR_C_H */

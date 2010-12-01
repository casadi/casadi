#ifndef SUNDIALS_INTERFACE_C_H 
#define SUNDIALS_INTERFACE_C_H

#ifdef __cplusplus
  extern "C" {
#endif 

#include "../casadi/c_interface/fx_c.h"

/* Create a CVodes integrator */
int casadi_cvodes_integrator(fx_ref fcn, fx_ref ffcn);

/* Create an IDAS integrator */
int casadi_idas_integrator(fx_ref fcn, fx_ref ffcn);

#ifdef __cplusplus
  }
#endif 

#endif /* SUNDIALS_INTERFACE_C_H */

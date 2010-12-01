#ifndef IPOPT_SOLVER_C_H 
#define IPOPT_SOLVER_C_H 

#ifdef __cplusplus
  extern "C" {
#endif 

#include "../casadi/c_interface/fx_c.h"

/* Create a function for numerical evaluation */
int casadi_ipopt_solver(fx_ref fcn, fx_ref f, fx_ref g, fx_ref h, fx_ref j);

#ifdef __cplusplus
  }
#endif 

#endif /* IPOPT_SOLVER_C_H */

#ifndef MX_FUNCTION_C_H 
#define MX_FUNCTION_C_H 

#ifdef __cplusplus
  extern "C" {
#endif 

#include "mx_c.h"
#include "fx_c.h"

/* Create a function for numerical evaluation */
int casadi_mx_function(fx_ref fcn, mx_vec iv, mx_vec ov);

#ifdef __cplusplus
  }
#endif 

#endif /* MX_FUNCTION_C_H */

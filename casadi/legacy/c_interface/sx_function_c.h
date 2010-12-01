#ifndef SX_FUNCTION_C_H 
#define SX_FUNCTION_C_H 

#ifdef __cplusplus
  extern "C" {
#endif 

#include "sx_matrix_c.h"
#include "fx_c.h"

/* Create a function for symbolical and numerical evaluation */
int casadi_sx_function(fx_ref fcn, sx_matrix_vec iv, sx_matrix_vec ov);

/* Evaluate a function symbolically (change -> should return a vector) */
int casadi_sx_function_evaluate_symbolically(fx_ref fcn, sx_matrix_vec iv, sx_matrix_ref res);

#ifdef __cplusplus
  }
#endif 

#endif /* SX_FUNCTION_C_H */

#ifndef CASADI_C_INTERFACE_H
#define CASADI_C_INTERFACE_H

#ifdef __cplusplus
extern "C" {
#endif

#include "mex.h"

double test(double aa);

void* convert_to_swig(mxArray *array_ptr);

mxArray* convert_from_swig(void* proxy_ptr);

void convert_input_mx(void* proxy_ptr);

void convert_input_swig(void* proxy_ptr);

#ifdef __cplusplus
}
#endif

#endif /* CASADI_C_INTERFACE_H */


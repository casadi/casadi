#ifndef CASADI_C_INTERFACE_H
#define CASADI_C_INTERFACE_H

#ifdef __cplusplus
extern "C" {
#endif

#include "mex.h"

double test(double aa);

int test_mex(const mxArray *array_ptr);

#ifdef __cplusplus
}
#endif

#endif /* CASADI_C_INTERFACE_H */


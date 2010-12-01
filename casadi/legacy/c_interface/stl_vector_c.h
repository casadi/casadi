#ifndef STL_VECTOR_C_H
#define STL_VECTOR_C_H

#ifdef __cplusplus
  extern "C" {
#endif 

#include <stdio.h>

/* STL vector class */
typedef void* vector_ptr;

/* Allocate STL string */
vector_ptr casadi_vector_new(void);

/* Delete vector */
int casadi_vector_delete(vector_ptr str);

/* Resize */
int casadi_vector_resize(vector_ptr v, int len);

/* Get length of array */
int casadi_vector_size(vector_ptr v);

/* Set the value of the array */
int casadi_vector_set(vector_ptr v, const double* val);

/* Get the value of the array */
int casadi_vector_get(vector_ptr v, double* val);

/* Get double array */
const double* casadi_vector_get_ptr(vector_ptr v);

#ifdef __cplusplus
  }
#endif 

#endif /* STL_VECTOR_C_H */

#ifndef SX_MATRIX_C_H
#define SX_MATRIX_C_H

#ifdef __cplusplus
  extern "C" {
#endif 

#include <stdio.h>
#include "stl_string_c.h"
#include "sx_c.h"

/* A pointer to an sx_matrix */
typedef void* sx_matrix_ref;

/* Vector of sx_matrix:es */
typedef void* sx_matrix_vec;

/* Create a matrix expression */
sx_matrix_ref casadi_sx_matrix_new();

/* Assign symbolic variable */
int casadi_sx_matrix_symbol(sx_matrix_ref ref, const char* name, int nrow, int ncol);

/* Assign constant */
int casadi_sx_matrix_constant(sx_matrix_ref ref, const double* data, int nrow, int ncol, char order);

/* Constructor that takes an symbolic scalar */
int casadi_sx_matrix_sx(sx_matrix_ref ref, sx_ptr scalar);

/* Constructor that takes a vector of symbolic scalars */
int casadi_sx_matrix_sx_vec(sx_matrix_ref ref, sx_vec_ptr v);

/* Delete a variable */
int casadi_sx_matrix_delete(sx_matrix_ref ref);

/* Print */
int casadi_sx_matrix_print(sx_matrix_ref ref, string_ptr str);
int casadi_sx_matrix_print_cout(sx_matrix_ref ref); /* remove */
int casadi_sx_matrix_print_cerr(sx_matrix_ref ref); /* remove */

/* Binary operation */
int casadi_sx_matrix_binary(sx_matrix_ref r, int op, sx_matrix_ref x, sx_matrix_ref y);

/* Unary operation */
int casadi_sx_matrix_unary(sx_matrix_ref r, int op, sx_matrix_ref x);

/* Matrix multiplication */
int casadi_sx_matrix_prod(sx_matrix_ref r, sx_matrix_ref x, sx_matrix_ref y);

/* Vertical concatenation */
int casadi_sx_matrix_vertcat(sx_matrix_ref r, sx_matrix_vec v);

/* Horizontal concatenation */
int casadi_sx_matrix_horzcat(sx_matrix_ref r, sx_matrix_vec v);

/* Create an sx_matrix vector */
sx_matrix_vec casadi_sx_matrix_vec_new(void);

/* delete an sx_matrix vector */
int casadi_sx_matrix_vec_delete(sx_matrix_vec v);

/* Add an element */
int casadi_sx_matrix_vec_push_back(sx_matrix_vec v, sx_matrix_ref ref);

/* Get the length of the the array */
int casadi_sx_matrix_vec_size(sx_matrix_vec v);

/* Print */
int casadi_sx_matrix_vec_print(sx_matrix_ref ref, string_ptr str);
int casadi_sx_matrix_vec_print_cout(sx_matrix_ref v);

/* Get the size of the matrix */
int casadi_sx_matrix_size(sx_matrix_ref ref, int *sz);
int casadi_sx_matrix_size1(sx_matrix_ref ref, int *sz);
int casadi_sx_matrix_size2(sx_matrix_ref ref, int *sz);

/* Transpose of a matrix */
int casadi_sx_matrix_transpose(sx_matrix_ref res, sx_matrix_ref ref);


#ifdef __cplusplus
  }
#endif 

#endif /* SX_MATRIX_C_H */

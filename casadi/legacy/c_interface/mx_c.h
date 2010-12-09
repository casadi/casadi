/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A symbolic framework for dynamic optimization.
 *    Copyright (C) 2010 by Joel Andersson, Moritz Diehl, K.U.Leuven. All rights reserved.
 *
 *    CasADi is free software; you can redistribute it and/or
 *    modify it under the terms of the GNU Lesser General Public
 *    License as published by the Free Software Foundation; either
 *    version 3 of the License, or (at your option) any later version.
 *
 *    CasADi is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *    Lesser General Public License for more details.
 *
 *    You should have received a copy of the GNU Lesser General Public
 *    License along with CasADi; if not, write to the Free Software
 *    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 *
 */

#ifndef MX_C_H
#define MX_C_H

#ifdef __cplusplus
  extern "C" {
#endif 

#include <stdio.h>
#include "stl_string_c.h"

/* A pointer to an mx */
typedef void* mx_ref;

/* Vector of mx:es */
typedef void* mx_vec;

/* Create a matrix expression */
mx_ref casadi_mx_new(void);

/* Assign symbolic variable */
int casadi_mx_symbol(mx_ref ref, const char* name, int nrow, int ncol);

/* Assign constant */
int casadi_mx_constant(mx_ref ref, const double* data, int nrow, int ncol, char order);

/* Delete a variable */
int casadi_mx_delete(mx_ref ref);

/* Print to a string */
int casadi_mx_print_string(mx_ref ref, string_ptr str);

/* Print to cout: redirect to print to a file */
int casadi_mx_print_cout(mx_ref ref);

/* Print to cerr: redirect to print to a file */
int casadi_mx_print_cerr(mx_ref ref);

/* Binary operation */
int casadi_mx_binary(mx_ref r, int op, mx_ref x, mx_ref y);

/* Unary operation */
int casadi_mx_unary(mx_ref r, int op, mx_ref x);

/* Matrix multiplication */
int casadi_mx_prod(mx_ref r, mx_ref x, mx_ref y);

/* Vertical concatenation */
int casadi_mx_vertcat(mx_ref r, mx_vec v);

/* Horizontal concatenation */
int casadi_mx_horzcat(mx_ref r, mx_vec v);

/* Create an mx vector */
mx_vec casadi_mx_vec_new(void);

/* delete an mx vector */
int casadi_mx_vec_delete(mx_vec v);

/* Add an element */
int casadi_mx_vec_push_back(mx_vec v, mx_ref ref);

/* Get the length of the the array */
int casadi_mx_vec_size(mx_vec v);

/* Print to array cout: redirect to print to a file */
int casadi_mx_vec_print_cout(mx_vec v);

/* Get the size of the matrix */
int casadi_mx_size(mx_ref ref, int *sz);
int casadi_mx_size1(mx_ref ref, int *sz);
int casadi_mx_size2(mx_ref ref, int *sz);

/* Transpose of a matrix */
int casadi_mx_transpose(mx_ref res, mx_ref ref);

#ifdef __cplusplus
  }
#endif 

#endif /* MX_C_H */

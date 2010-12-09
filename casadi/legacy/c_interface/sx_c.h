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

#ifndef SX_C_H
#define SX_C_H

#ifdef __cplusplus
  extern "C" {
#endif 

#include <stdio.h>
#include "stl_string_c.h"

/* A pointer to an sx */
typedef void* sx_ptr;

/* Vector of sx:es */
typedef void* sx_vec_ptr;

/* Create a matrix expression */
sx_ptr casadi_sx_new(void);

/* Assign symbolic variable */
int casadi_sx_symbol(sx_ptr ptr, const char* name);

/* Assign constant */
int casadi_sx_constant(sx_ptr ptr, double value);

/* Delete a variable */
int casadi_sx_delete(sx_ptr ptr);

/* Print to a string */
int casadi_sx_print(sx_ptr ptr, string_ptr str);

/* Binary (or unary) operation */
int casadi_sx_binary(sx_ptr r, int op, sx_ptr x, sx_ptr y);

/* Unary operation */
int casadi_sx_unary(sx_ptr r, int op, sx_ptr x);

/* Create an sx vector */
sx_vec_ptr casadi_sx_vec_new(void);

/* delete an sx vector */
int casadi_sx_vec_delete(sx_vec_ptr v);

/* Add an element */
int casadi_sx_vec_push_back(sx_vec_ptr v, sx_ptr ptr);

/* Get the length of the the array */
int casadi_sx_vec_size(sx_vec_ptr v);

/* Print to array cout: redirect to print to a file */
int casadi_sx_vec_print(sx_vec_ptr ptr, string_ptr str);

#ifdef __cplusplus
  }
#endif 

#endif /* SX_C_H */

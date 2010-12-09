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

/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A minimalistic computer algebra system with automatic differentiation 
 *    and framework for dynamic optimization.
 *    Copyright (C) 2010 by Joel Andersson, Moritz Diehl et al., K.U.Leuven. All rights reserved.
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

#ifndef STL_STRING_C_H
#define STL_STRING_C_H

#ifdef __cplusplus
  extern "C" {
#endif 

#include <stdio.h>

/* STL string class */
typedef void* string_ptr;

/* Allocate STL string */
string_ptr casadi_string_new(void);

/* Delete STL string */
int casadi_string_delete(string_ptr str);

/* Assign STL string */
int casadi_string_assign(string_ptr str, const char* s);

/* Get char array */
const char* casadi_string_get(string_ptr str);

#ifdef __cplusplus
  }
#endif 

#endif /* STL_STRING_C_H */

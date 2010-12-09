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

#ifndef FX_C_H
#define FX_C_H

#ifdef __cplusplus
  extern "C" {
#endif 

#include <stdio.h>
#include "stl_string_c.h"
#include "stl_vector_c.h"

/* A pointer to an mx */
typedef void* fx_ref;

/* Matlab work around */
int casadi_ptr_to_int(void* ref, int *ref_int); // 32 bit systems
int casadi_ptr_to_long(void* ref, long *ref_int); // 64  bit systems

/* Create a matrix expression */
fx_ref casadi_fx_new(void);

/* Delete a variable */
int casadi_fx_delete(fx_ref ref);

/* Print to cout: redirect to print to a file */
int casadi_fx_print_cout(fx_ref ref);

/* Print to cerr: redirect to print to a file */
int casadi_fx_print_cerr(fx_ref ref);

/* Print to an internal string */
int casadi_fx_print_string(fx_ref ref, string_ptr str);

/* Set an string valued option */
int casadi_fx_setoption_string(fx_ref ref, const char* name, const char* opval);

/* Set a double valued option */
int casadi_fx_setoption_double(fx_ref ref, const char* name, const double* opval, int n);

/* Get an string valued option */
int casadi_fx_getoption_string(fx_ref ref, const char* name, string_ptr str);

/* Get a double valued option */
int casadi_fx_getoption_double(fx_ref ref, const char* name, vector_ptr str);

/* Get the option type */
int casadi_fx_option_is_string(fx_ref ref, const char* name, int *is_string);

/* Print an option */
int casadi_fx_print_option(fx_ref ref, const char* name, string_ptr str);

/* Print all options */
int casadi_fx_print_options(fx_ref ref);

/* Initialize the expression */
int casadi_fx_init(fx_ref ref);

/* Get input size */
int casadi_fx_input_size(fx_ref ref, int ind, int *sz);

/* Get input value */
int casadi_fx_getinput(fx_ref ref, int ind, int ord, double* val);

/* Set input value */
int casadi_fx_setinput(fx_ref ref, int ind, int ord, const double* val);

/* Get output size */
int casadi_fx_output_size(fx_ref ref, int ind, int *sz);

/* Get outputvalue */
int casadi_fx_getoutput(fx_ref ref, int ind, int ord, double* val);

/* Set output value */
int casadi_fx_setoutput(fx_ref ref, int ind, int ord, const double* val);

/* Evaluate */
int casadi_fx_evaluate(fx_ref ref);

/* Evaluate forward AD*/
// int casadi_fx_evaluate_fwd(fx_ref ref);

/* Evaluate adjoint AD */
// int casadi_fx_evaluate_adj(fx_ref ref);

/* Jacobian of an expression */
int casadi_fx_jacobian(fx_ref ref, fx_ref res, int iind, int oind);

/* Hessian of an expression */
int casadi_fx_hessian(fx_ref ref, fx_ref res, int iind, int oind);



#ifdef __cplusplus
  }
#endif 

#endif /* FX_C_H */

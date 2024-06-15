/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A symbolic framework for dynamic optimization.
 *    Copyright (C) 2010-2023 Joel Andersson, Joris Gillis, Moritz Diehl,
 *                            KU Leuven. All rights reserved.
 *    Copyright (C) 2011-2014 Greg Horn
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


#ifndef CASADI_C_H
#define CASADI_C_H

#ifdef __cplusplus
extern "C" {
#endif

#ifndef CASADI_CASADI_HPP

/* Used to manage symbol visibility */
#ifndef CASADI_EXPORT
#define CASADI_EXPORT
#endif /* CASADI_EXPORT */

/* Integer type */
#ifndef casadi_int
#define casadi_int long long int
#endif /* casadi_int */

#endif /* CASADI_CASADI_HPP */

/* ==================================
*   Managing CasADi Functions
*  ================================== */

/** \brief Read from a serialized CasADi file
 * 
 * Pushes to an internal buffer any Functions found
 * 
 * Not thread-safe
 * Return 0 when successful
*/

CASADI_EXPORT int casadi_c_push_file(const char *filename);

/** \brief Unloads the last batch of added Functions
* Not thread-safe
*/
CASADI_EXPORT void casadi_c_pop(void);

/** \brief Unloads all Functions
 * 
 * Not thread-safe
 */
CASADI_EXPORT void casadi_c_clear(void);

/** \brief Show the amount of Functions loaded */
CASADI_EXPORT int casadi_c_n_loaded(void);


/** \brief Get a Function id by name
 *
 * id is really just a linear index into the loaded Functions vector 
*/
CASADI_EXPORT int casadi_c_id(const char* funname);
CASADI_EXPORT int casadi_c_activate(int id);


/** \brief Get name of Function */
CASADI_EXPORT const char* casadi_c_name();
CASADI_EXPORT const char* casadi_c_name_id(int id);

/** \brief Get width of integer type casadi_int */
CASADI_EXPORT int casadi_c_int_width();
/** \brief Get width of real type double */
CASADI_EXPORT int casadi_c_real_width();

/* ===================================================
*   Codegen-like API for evaluating CasADi Functions
*  =================================================== */


CASADI_EXPORT void casadi_c_incref(void);
CASADI_EXPORT void casadi_c_decref(void);
CASADI_EXPORT int casadi_c_checkout(void);
CASADI_EXPORT void casadi_c_release(int mem);
CASADI_EXPORT double casadi_c_default_in(casadi_int i);
CASADI_EXPORT casadi_int casadi_c_n_in(void);
CASADI_EXPORT casadi_int casadi_c_n_out(void);
CASADI_EXPORT const char* casadi_c_name_in(casadi_int i);
CASADI_EXPORT const char* casadi_c_name_out(casadi_int i);
CASADI_EXPORT const casadi_int* casadi_c_sparsity_in(casadi_int i);
CASADI_EXPORT const casadi_int* casadi_c_sparsity_out(casadi_int i);
CASADI_EXPORT int casadi_c_work(casadi_int *sz_arg, casadi_int* sz_res,
  casadi_int *sz_iw, casadi_int *sz_w);
CASADI_EXPORT int casadi_c_eval(const double** arg, double** res,
  casadi_int* iw, double* w, int mem);


CASADI_EXPORT void casadi_c_incref_id(int id);
CASADI_EXPORT void casadi_c_decref_id(int id);
CASADI_EXPORT int casadi_c_checkout_id(int id);
CASADI_EXPORT void casadi_c_release_id(int id, int mem);
CASADI_EXPORT double casadi_c_default_in_id(int id, casadi_int i);
CASADI_EXPORT casadi_int casadi_c_n_in_id(int id);
CASADI_EXPORT casadi_int casadi_c_n_out_id(int id);
CASADI_EXPORT const char* casadi_c_name_in_id(int id, casadi_int i);
CASADI_EXPORT const char* casadi_c_name_out_id(int id, casadi_int i);
CASADI_EXPORT const casadi_int* casadi_c_sparsity_in_id(int id, casadi_int i);
CASADI_EXPORT const casadi_int* casadi_c_sparsity_out_id(int id, casadi_int i);
CASADI_EXPORT int casadi_c_work_id(int id, casadi_int *sz_arg, casadi_int* sz_res,
  casadi_int *sz_iw, casadi_int *sz_w);
CASADI_EXPORT int casadi_c_eval_id(int id, const double** arg, double** res,
  casadi_int* iw, double* w, int mem);

CASADI_EXPORT void casadi_c_logger_write(const char* msg, int num);
CASADI_EXPORT void casadi_c_logger_flush(void);

#ifdef __cplusplus
}
#endif

#endif /* CASADI_C_H */

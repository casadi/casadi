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

/* ===================================================
*   Stats, as generated code with option "stats"
*  =================================================== */

/* Standalone copy of runtime/casadi_stats_sink_callback.hpp (C89) */
#ifndef CASADI_STATS_SINK
#define CASADI_STATS_SINK
struct casadi_stats_sink {
  unsigned char* p;
  casadi_int cap;
  casadi_int needed;
  unsigned char* (*reserve)(struct casadi_stats_sink* s, casadi_int len, casadi_int* pos);
  void (*commit)(struct casadi_stats_sink* s, casadi_int pos);
  void* data;
};
#endif /* CASADI_STATS_SINK */

/** \brief Sink on p[0, cap), thread-safe */
CASADI_EXPORT void casadi_c_stats_init_sink(struct casadi_stats_sink* s, unsigned char* p,
  casadi_int cap);
/** \brief Drop all records */
CASADI_EXPORT void casadi_c_stats_clear(struct casadi_stats_sink* s);
/** \brief Bytes the stream needed, including dropped records */
CASADI_EXPORT casadi_int casadi_c_stats_nbytes(const struct casadi_stats_sink* s);
/** \brief Were records dropped? */
CASADI_EXPORT int casadi_c_stats_truncated(const struct casadi_stats_sink* s);
/** \brief Recorded bytes held in memory: *n of them */
CASADI_EXPORT const unsigned char* casadi_c_stats_data(const struct casadi_stats_sink* s,
  casadi_int* n);
/** \brief Store where each node ends, after the last write: queries skip subtrees */
CASADI_EXPORT void casadi_c_stats_index(struct casadi_stats_sink* s);
/** \brief Next call after cursor (-1: first) at any depth below scope (-1: root), or -1

    A pattern containing ':' matches the full id, otherwise the function name; null: any */
CASADI_EXPORT casadi_int casadi_c_stats_find_function(const struct casadi_stats_sink* s,
  casadi_int scope, casadi_int cursor, const char* pattern);
/** \brief As casadi_c_stats_find_function, among the direct children of parent */
CASADI_EXPORT casadi_int casadi_c_stats_select_function(const struct casadi_stats_sink* s,
  casadi_int parent, casadi_int cursor, const char* pattern);
/** \brief Next iteration of parent after cursor (-1: first) with index (-1: any), or -1 */
CASADI_EXPORT casadi_int casadi_c_stats_select_iteration(const struct casadi_stats_sink* s,
  casadi_int parent, casadi_int cursor, casadi_int index);
/** \brief Iteration of parent recorded last, or -1 */
CASADI_EXPORT casadi_int casadi_c_stats_select_last_iteration(
  const struct casadi_stats_sink* s, casadi_int parent);
/** \brief Next section of parent after cursor (-1: first) named name (null: any), or -1 */
CASADI_EXPORT casadi_int casadi_c_stats_select_section(const struct casadi_stats_sink* s,
  casadi_int parent, casadi_int cursor, const char* name);
/** \brief Value of key of node; returns 0 if found

    Reserved keys: id, name, mem and flag of a call, name of a section, index of an
    iteration */
CASADI_EXPORT int casadi_c_stats_get_int(const struct casadi_stats_sink* s, casadi_int node,
  const char* key, casadi_int* v);
CASADI_EXPORT int casadi_c_stats_get_real(const struct casadi_stats_sink* s, casadi_int node,
  const char* key, double* v);
/** \brief At most cap-1 chars in v, null-terminated if cap>0

    *len: full length without the null, truncated if *len >= cap (len may be null) */
CASADI_EXPORT int casadi_c_stats_get_text(const struct casadi_stats_sink* s, casadi_int node,
  const char* key, char* v, casadi_int cap, casadi_int* len);
/** \brief First cap elements in v; *len: full length (len may be null) */
CASADI_EXPORT int casadi_c_stats_get_reals(const struct casadi_stats_sink* s, casadi_int node,
  const char* key, double* v, casadi_int cap, casadi_int* len);

/** \brief Evaluate, recording to sink under parent (-1 for a root call) */
CASADI_EXPORT int casadi_c_eval_with_stats(const double** arg, double** res,
  casadi_int* iw, double* w, int mem, struct casadi_stats_sink* sink, casadi_int parent);
CASADI_EXPORT int casadi_c_eval_with_stats_id(int id, const double** arg, double** res,
  casadi_int* iw, double* w, int mem, struct casadi_stats_sink* sink, casadi_int parent);

CASADI_EXPORT void casadi_c_logger_write(const char* msg, int num);
CASADI_EXPORT void casadi_c_logger_flush(void);

#ifdef __cplusplus
}
#endif

#endif /* CASADI_C_H */

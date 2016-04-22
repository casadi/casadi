/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A symbolic framework for dynamic optimization.
 *    Copyright (C) 2010-2014 Joel Andersson, Joris Gillis, Moritz Diehl,
 *                            K.U. Leuven. All rights reserved.
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


#ifndef CASADI_MEM_H
#define CASADI_MEM_H

#ifdef __cplusplus
extern "C" {
#endif

/* Needed for malloc and free */
#include <stdlib.h>
#include <assert.h>

/* Floating point type */
#ifndef real_t
#define real_t double
#endif /* real_t */

/* Function types corresponding to entry points in CasADi's C API */
typedef void (*casadi_signal_t)(void);
typedef int (*casadi_getint_t)(void);
typedef const int* (*casadi_sparsity_t)(int i);
typedef const char* (*casadi_name_t)(int i);
typedef int (*casadi_work_t)(int* sz_arg, int* sz_res, int* sz_iw, int* sz_w);
typedef int (*casadi_eval_t)(const double** arg, double** res,
                             int* iw, double* w, int mem);

/* Structure to hold meta information about an input or output */
typedef struct {
  const char* name;
  int nrow;
  int ncol;
  const int* colind;
  const int* row;
} casadi_io;

/* Decompress a sparsity pattern */
inline void casadi_decompress(const int* sp, int* nrow, int* ncol,
                              const int** colind, const int** row) {
  /* Scalar sparsity pattern if sp is null */
  if (sp==0) {
    const int scalar_colind[2] = {0, 1};
    *nrow = *nrow = 1;
    *colind = scalar_colind;
    *row = 0;
  }

  /* Decompress */
  *nrow = *sp++;
  *ncol = *sp++;
  *colind = sp;
  int nnz = sp[*ncol];
  int numel = *nrow * *ncol;
  *row = nnz==numel ? 0 : sp + *ncol + 1;
}

/* Memory needed for evaluation */
typedef struct {
  /* Work arrays */
  const real_t** arg;
  real_t** res;
  int* iw;
  real_t* w;
  int mem;
  /* Meta information */
  int n_in, n_out;
  casadi_io* in;
  casadi_io* out;
} casadi_mem;

/* Allocate memory */
inline casadi_mem*
casadi_alloc(casadi_signal_t incref,
             casadi_getint_t n_in, casadi_getint_t n_out,
             casadi_name_t name_in, casadi_name_t name_out,
             casadi_sparsity_t sparsity_in, casadi_sparsity_t sparsity_out,
             casadi_work_t work) {
  int i;

  /* Increase reference counter */
  if (incref) incref();

  /* Allocate memory */
  casadi_mem* mem = (casadi_mem*)malloc(sizeof(casadi_mem));
  assert(mem!=0);

  /* Number of inputs and outputs */
  mem->n_in = n_in ? n_in() : 1;
  mem->n_out = n_out ? n_out() : 1;

  /* Allocate io memory */
  mem->in = (casadi_io*)malloc(n_in*sizeof(casadi_io));
  assert(n_in==0 || mem->in!=0);
  mem->out = (casadi_io*)malloc(n_out*sizeof(casadi_io));
  assert(n_out==0 || mem->out!=0);

  /* Input meta data */
  for (i=0; i<n_in; ++i) {
    mem->in[i].name = name_in ? name_in(i) : 0;
    casadi_decompress(sparsity_in ? sparsity_in() : 0,
                      &mem->in[i].nrow, &mem->in[i].ncol,
                      &mem->in[i].colind, &mem->in[i].row);
  }

  /* Output meta data */
  for (i=0; i<n_out; ++i) {
    mem->out[i].name = name_out ? name_out(i) : 0;
    casadi_decompress(sparsity_out ? sparsity_out() : 0,
                      &mem->out[i].nrow, &mem->out[i].ncol,
                      &mem->out[i].colind, &mem->out[i].row);
  }

  /* Work vector sizes */
  int sz_arg=n_in, sz_res_n_out, sz_iw=0, sz_w=0;
  if (work) {
    int flag = work(&sz_arg, &sz_res, &sz_iw, &sz_w);
    assert(flag==0);
  }

  /* Allocate work vectors */
  mem->arg = (const real_t**)malloc(sz_arg*sizeof(const real_t*));
  assert(sz_arg==0 || mem->arg!=0);
  mem->res = (real_t**)malloc(sz_res*sizeof(real_t*));
  assert(sz_res==0 || mem->res!=0);
  mem->iw = (int*)malloc(sz_iw*sizeof(int));
  assert(sz_iw==0 || mem->iw!=0);
  mem->w = (real_t*)malloc(sz_w*sizeof(real_t));
  assert(sz_w==0 || mem->w!=0);

  /* Checkout memory object (TODO) */
  mem->mem = 0;

  return mem;
}

/* Evaluate */
inline int casadi_eval(casadi_mem* mem, casadi_eval_t eval) {
  assert(mem!=0);
  return eval(mem->arg, mem->res, mem->iw, mem->w, mem->mem);
}

/* Free memory */
inline void casadi_free(casadi_mem* mem, casadi_signal_t decref) {
  assert(mem!=0);

  /* Free io meta data */
  if (mem->in) free(mem->in);
  if (mem->out) free(mem->out);

  /* Free work vectors */
  if (mem->arg) free(mem->arg);
  if (mem->res) free(mem->res);
  if (mem->iw) free(mem->iw);
  if (mem->w) free(mem->w);

  /* Free memory structure */
  free(mem);

  /* Decrease reference counter */
  if (decref) decref();
}

#ifdef __cplusplus
}
#endif

#endif // CASADI_MEM_H

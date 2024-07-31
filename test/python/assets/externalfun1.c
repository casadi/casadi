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

#ifdef __cplusplus
extern "C" {
#endif

/* How to prefix internal symbols */
#ifdef CASADI_CODEGEN_PREFIX
  #define CASADI_NAMESPACE_CONCAT(NS, ID) _CASADI_NAMESPACE_CONCAT(NS, ID)
  #define _CASADI_NAMESPACE_CONCAT(NS, ID) NS ## ID
  #define CASADI_PREFIX(ID) CASADI_NAMESPACE_CONCAT(CODEGEN_PREFIX, ID)
#else
  #define CASADI_PREFIX(ID) f_ ## ID
#endif

#include <math.h>
#include <stdio.h>
#include <stdarg.h>
#include <string.h>

#ifndef casadi_real
#define casadi_real double
#endif

#ifndef casadi_int
#define casadi_int long long int
#endif

#define CASADI_MAX_NUM_THREADS 2

/* Symbol visibility in DLLs */
#ifndef CASADI_SYMBOL_EXPORT
  #if defined(_WIN32) || defined(__WIN32__) || defined(__CYGWIN__)
    #if defined(STATIC_LINKED)
      #define CASADI_SYMBOL_EXPORT
    #else
      #define CASADI_SYMBOL_EXPORT __declspec(dllexport)
    #endif
  #elif defined(__GNUC__) && defined(GCC_HASCLASSVISIBILITY)
    #define CASADI_SYMBOL_EXPORT __attribute__ ((visibility ("default")))
  #else
    #define CASADI_SYMBOL_EXPORT
  #endif
#endif

#define casadi_s1 CASADI_PREFIX(s1)

static const casadi_int casadi_s1[5] = {1, 1, 0, 1, 0};

struct demo_external_local {
  casadi_int foo;
};

struct demo_external_global {
  casadi_int n_flags;
};

static int casadi_F_mem_counter = 0;
static int casadi_F_unused_stack_counter = -1;
static int casadi_F_unused_stack[CASADI_MAX_NUM_THREADS];
static struct demo_external_local casadi_F_mem[CASADI_MAX_NUM_THREADS];


static int demo_external_global_counter = 0;
static struct demo_external_global demo_external_global_data;

static int demo_external_argc = 0;
static const char** demo_external_argv = 0;


CASADI_SYMBOL_EXPORT int F(const casadi_real** arg, casadi_real** res, casadi_int* iw, casadi_real* w, int mem){
  struct demo_external_global* g = &demo_external_global_data;
  struct demo_external_local* m = &casadi_F_mem[mem];
  if (res[0]) res[0][0] = g->n_flags;
  if (res[1]) res[1][0] = m->foo;
  if (res[2]) res[2][0] = mem;
  return 0;
}

CASADI_SYMBOL_EXPORT int F_alloc_mem(void) {
  return casadi_F_mem_counter++;
}

CASADI_SYMBOL_EXPORT int F_init_mem(int mem) {
  struct demo_external_local* m = &casadi_F_mem[mem];
  m->foo = mem;
  return 0;
}

CASADI_SYMBOL_EXPORT void F_free_mem(int mem) {

}


CASADI_SYMBOL_EXPORT int F_checkout(void) {
  int mid;
  printf("checkout()\n");
  if (casadi_F_unused_stack_counter>=0) {
    return casadi_F_unused_stack[casadi_F_unused_stack_counter--];
  } else {
    if (casadi_F_mem_counter==CASADI_MAX_NUM_THREADS) return -1;
    mid = F_alloc_mem();
    if (mid<0) return -1;
    if(F_init_mem(mid)) return -1;
    return mid;
  }
}

CASADI_SYMBOL_EXPORT void F_release(int mem) {
  printf("release(%d)\n", mem);
  casadi_F_unused_stack[++casadi_F_unused_stack_counter] = mem;
}

CASADI_SYMBOL_EXPORT int F_config(int argc, const char** argv) {
  demo_external_argc = argc;
  demo_external_argv = argv;
  if (demo_external_argc>=3) return 0; /* success */
  return 1;
}
 
CASADI_SYMBOL_EXPORT void F_incref(void) {
  struct demo_external_global* g = &demo_external_global_data;
  if (demo_external_global_counter==0) {
    g->n_flags = demo_external_argc;
  }
  demo_external_global_counter++;
}

CASADI_SYMBOL_EXPORT void F_decref(void) {
  demo_external_global_counter--;
}

CASADI_SYMBOL_EXPORT casadi_int F_n_in(void) {
  return 1;
}

CASADI_SYMBOL_EXPORT casadi_int F_n_out(void) {
  return 3;
}

CASADI_SYMBOL_EXPORT const char* F_name_in(casadi_int i) {
  return "x";
}

CASADI_SYMBOL_EXPORT const char* F_name_out(casadi_int i) {
  if (i==0) {
      return "flags";
  } else if (i==1) {
      return "mem_obj";
  } else {
      return "mem";
  }
}

CASADI_SYMBOL_EXPORT const casadi_int* F_sparsity_in(casadi_int i) {
  return casadi_s1;
}

CASADI_SYMBOL_EXPORT const casadi_int* F_sparsity_out(casadi_int i) {
  return casadi_s1;
}


#ifdef __cplusplus
} /* extern "C" */
#endif

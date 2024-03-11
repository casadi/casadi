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
#include <engine.h>
#include <stdio.h>
#include <stdarg.h>
#include <string.h>

#ifndef casadi_real
#define casadi_real double
#endif

#ifndef casadi_int
#define casadi_int long long int
#endif

#ifndef CASADI_MAX_NUM_THREADS
#define CASADI_MAX_NUM_THREADS 1
#endif

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

/* Structure to hold meta information about an input or output */
typedef struct {
  char* name;
  casadi_int nrow;
  casadi_int ncol;
  casadi_int nnz;
  casadi_int* compressed;
} casadi_io;

struct matlab_external_local {
  mxArray** arg;
  mxArray** res;

  mxArray* args; // cell
};

struct matlab_external_global {
  Engine *ep;
  casadi_int n_in;
  casadi_int n_out;
  casadi_io* in;
  casadi_io* out;
};

static int casadi_F_mem_counter = 0;
static int casadi_F_unused_stack_counter = -1;
static int casadi_F_unused_stack[CASADI_MAX_NUM_THREADS];
static struct matlab_external_local casadi_F_mem[CASADI_MAX_NUM_THREADS];


static int matlab_external_global_counter = 0;
static struct matlab_external_global matlab_external_global_data;

static int matlab_external_argc = 0;
static const char** matlab_external_argv = 0;

casadi_int to_val(mxArray* m, casadi_int* val) {
      if (mxGetNumberOfElements(m)>1) return 1;
      switch (mxGetClassID(m)) {
        case mxINT8_CLASS:
          if (val) *val = (casadi_int)(*(int8_T*) mxGetData(m));
          break;
        case mxUINT8_CLASS:
          if(val) *val = (casadi_int)(*(uint8_T*) mxGetData(m));
          break;
        case mxINT16_CLASS:
          if(val) *val = (casadi_int)(*(int16_T*) mxGetData(m));
          break;
        case mxUINT16_CLASS:
          if(val) *val = (casadi_int)(*(uint16_T*) mxGetData(m));
          break;
        case mxINT32_CLASS:
          if(val) *val = (casadi_int)(*(int32_T*) mxGetData(m));
          break;
        case mxUINT32_CLASS:
          if(val) *val = (casadi_int)(*(uint32_T*) mxGetData(m));
          break;
        case mxINT64_CLASS:
          if(val) *val = (casadi_int)(*(int64_T*) mxGetData(m));
          break;
        case mxUINT64_CLASS: 
          if(val) *val = (casadi_int)(*(uint64_T*) mxGetData(m));
          break;
        case mxDOUBLE_CLASS:
        {
          double v = mxGetScalar(m);
          if (v!=floor(v)) return 1;
          if (val) *val = (casadi_int) v;
          break;
        }
        default:
          return 1;
    }
    return 0;
}

void casadi_copy(const double* x, casadi_int n, double* y) {
  casadi_int i;
  if (y) {
    if (x) {
      for (i=0; i<n; ++i) *y++ = *x++;
    } else {
      for (i=0; i<n; ++i) *y++ = 0.;
    }
  }
}

CASADI_SYMBOL_EXPORT int F(const casadi_real** arg, casadi_real** res, casadi_int* iw, casadi_real* w, int mem){
  struct matlab_external_global* g = &matlab_external_global_data;
  struct matlab_external_local* m = &casadi_F_mem[mem];

  for (int i=0;i<g->n_in;++i) {
    casadi_copy(arg[i], g->in[i].nnz, (double*) mxGetData(m->arg[i]));
  }

  engPutVariable(g->ep, "temp", m->args);
  engEvalString(g->ep, "ret = cb.eval(temp);");

  mxArray* ret = engGetVariable(g->ep, "ret");
  if (!ret) return 1;

  for (int i=0;i<g->n_out;++i) {
    mxArray* mm = mxGetCell(ret, i);
    if (mxGetClassID(mm)!=mxDOUBLE_CLASS) return 1;
    casadi_copy((double*) mxGetData(mm), g->out[i].nnz, res[i]);
  }

  mxDestroyArray(ret);

  return 0;
}

CASADI_SYMBOL_EXPORT int F_alloc_mem(void) {
  return casadi_F_mem_counter++;
}

CASADI_SYMBOL_EXPORT int F_init_mem(int mem) {
  struct matlab_external_global* g = &matlab_external_global_data;
  struct matlab_external_local* m = &casadi_F_mem[mem];
  m->arg = (mxArray**) malloc(g->n_in*sizeof(mxArray*));
  m->args = mxCreateCellMatrix(1, g->n_in);

  for (int i=0;i<g->n_in;++i) {
    m->arg[i] = mxCreateDoubleMatrix(g->in[i].nrow, g->in[i].ncol, mxREAL);
    // Cell takes ownership
    mxSetCell(m->args, i, m->arg[i]);
  }

  return 0;
}

CASADI_SYMBOL_EXPORT void F_free_mem(int mem) {
  struct matlab_external_local* m = &casadi_F_mem[mem];
  mxDestroyArray(m->args);
  free(m->arg);
}


CASADI_SYMBOL_EXPORT int F_checkout(void) {
  int mid;
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
  casadi_F_unused_stack[++casadi_F_unused_stack_counter] = mem;

}

mxArray* engCall(Engine *e, const char* cmd) {
  char ct[1000] = {'t','e','m','p','='};

  snprintf(ct+5, sizeof(ct)-5, "%s", cmd);

  if (engEvalString(e, ct)) return 0;
  return engGetVariable(e, "temp");
}

mxArray* engCall_d(Engine *e, const char* cmd, casadi_int i) {
  char ct[1000] = {'t','e','m','p','='};

  snprintf(ct+5, sizeof(ct)-5, cmd, i);

  if (engEvalString(e, ct)) return 0;
  return engGetVariable(e, "temp");
}

mxArray* engCall_var(Engine *e, const char* cmd, ...) {
  char ct[1000] = {'t','e','m','p','='};

  va_list args;

  // Start variadic argument processing
  va_start(args, cmd);
  snprintf(ct+5, sizeof(ct)-5, cmd, args);
  // End variadic argument processing
  va_end(args);

  if (engEvalString(e, ct)) return 0;
  return engGetVariable(e, "temp");
}

int engCallTo_casadi_int(Engine *e, casadi_int* ret, const char* cmd) {
  mxArray* temp = engCall(e, cmd);
  if (!temp) return 1;
  to_val(temp, ret);
  mxDestroyArray(temp);
  return 0;
}

int engCallTo_char(Engine *e, char** ret, const char* cmd, int i) {
  mxArray* temp = engCall_d(e, cmd, i);
  if (!temp) return 1;
  if (!mxIsChar(temp)) {
    mxDestroyArray(temp);
    return 1;
  }
  casadi_int n = mxGetN(temp);
  *ret = malloc(n+1);
  mxGetString(temp, *ret, (n+1)*sizeof(char));
  mxDestroyArray(temp);
  return 0;
}


int engCallTo_casadi_io(Engine *e, casadi_io* ret, const char* cmd, int i) {
  mxArray* temp = engCall_d(e, cmd, i);
  if (!temp) return 1;
  if (mxGetClassID(temp)!=mxINT64_CLASS) {
    mxDestroyArray(temp);
    return 1;
  }

  const int64_T* sp = (const int64_T*) mxGetData(temp);

  int m = mxGetM(temp);
  int n = mxGetN(temp);

  if (mxGetNumberOfDimensions(temp)!=2 || (m!=1 && n!=1) ) {
    mxDestroyArray(temp);
    return 1;
  }

  if (m>n) n = m;

  ret->nrow = sp[0];
  ret->ncol = sp[1];
  if (n==2) {
    ret->nnz = sp[0]*sp[1];
    ret->compressed = malloc(sizeof(casadi_int)*3);
    ret->compressed[0] = sp[0];
    ret->compressed[1] = sp[1];
    ret->compressed[2] = 1; // Dense indicator
  } else {
    ret->nnz = sp[2+ret->ncol];
    ret->compressed = malloc(sizeof(casadi_int)*(2+ret->ncol+1+ret->nnz));
    for (int i=0;i<2+ret->ncol+1+ret->nnz;++i) ret->compressed[i] = sp[i];
  }

  mxDestroyArray(temp);
  return 0;
}

CASADI_SYMBOL_EXPORT int F_config(int argc_, const char** argv_) {
  matlab_external_argc = argc_;
  matlab_external_argv = argv_;
  printf("config %s\n", argv_[0]);
  if (matlab_external_argc>=3) return 0; // success
  return 1;
}
 
CASADI_SYMBOL_EXPORT void F_incref(void) {
  struct matlab_external_global* g = &matlab_external_global_data;
  if (matlab_external_global_counter==0) {

    int len = 10;
    for (int i = 2; i < matlab_external_argc; i++) {
      len += strlen(matlab_external_argv[i])+3;
    }
    char ct[strlen(matlab_external_argv[1])+12];
    char ct2[len];

    printf(":%s:\n",matlab_external_argv[1]);

    snprintf(ct, sizeof(ct), "addpath('%s')", matlab_external_argv[1]);
    snprintf(ct2, sizeof(ct2), "cb = %s(", matlab_external_argv[2]);

    // Append additional arguments to the function call
    for (int i = 3; i < matlab_external_argc; i++) {
        strncat(ct2, "'", sizeof(ct2) - strlen(ct2) - 1);
        strncat(ct2, matlab_external_argv[i], sizeof(ct2) - strlen(ct2) - 1);
        strncat(ct2, "'", sizeof(ct2) - strlen(ct2) - 1);
        if (i < matlab_external_argc - 1) {
            strncat(ct2, ",", sizeof(ct2) - strlen(ct2) - 1);
        }
    }
    strncat(ct2, ")", sizeof(ct2) - strlen(ct2) - 1);

    printf("%s\n",ct);
    printf("%s\n",ct2);

    g->ep = engOpen("\0");
    if (!g->ep) {
      printf("Can't start MATLAB engine\n");
    }


    engEvalString(g->ep, ct);
    engEvalString(g->ep, ct2);
    engCallTo_casadi_int(g->ep, &g->n_in, "cb.get_n_in()");
    engCallTo_casadi_int(g->ep, &g->n_out, "cb.get_n_out()");
    g->in = (casadi_io*) malloc(g->n_in*sizeof(casadi_io));
    g->out = (casadi_io*) malloc(g->n_out*sizeof(casadi_io));

    for (int i=0;i<g->n_in;++i) {
      g->in[i].name = 0;
      engCallTo_casadi_io(g->ep, g->in+i, "cb.get_sparsity_in(%d)", i);
      engCallTo_char(g->ep, &g->in[i].name, "cb.get_name_in(%d)", i);
    }
    for (int i=0;i<g->n_out;++i) {
      g->out[i].name = 0;
      engCallTo_casadi_io(g->ep, g->out+i, "cb.get_sparsity_out(%d)", i);
      engCallTo_char(g->ep, &g->out[i].name, "cb.get_name_out(%d)", i);
    }
  }
  matlab_external_global_counter++;
}

CASADI_SYMBOL_EXPORT void F_decref(void) {
  int i;
  struct matlab_external_global* g = &matlab_external_global_data;
  matlab_external_global_counter--;
  if (matlab_external_global_counter==0) {
    engClose(g->ep);
    for (int i=0;i<g->n_in;++i) {
      free(g->in[i].compressed);
      free(g->in[i].name);
    }
    for (int i=0;i<g->n_out;++i) {
      free(g->out[i].compressed);
      free(g->out[i].name);
    }
    free(g->in);
    free(g->out);
    for (i=0;i<casadi_F_mem_counter;++i) {
      F_free_mem(i);
    }
  }
}

CASADI_SYMBOL_EXPORT casadi_int F_n_in(void) {
  struct matlab_external_global* g = &matlab_external_global_data;
  return g->n_in;
}

CASADI_SYMBOL_EXPORT casadi_int F_n_out(void) {
  struct matlab_external_global* g = &matlab_external_global_data;
  return g->n_out;
}

CASADI_SYMBOL_EXPORT const char* F_name_in(casadi_int i) {
  struct matlab_external_global* g = &matlab_external_global_data;
  return g->in[i].name;
}

CASADI_SYMBOL_EXPORT const char* F_name_out(casadi_int i) {
  struct matlab_external_global* g = &matlab_external_global_data;
  return g->out[i].name;
}

CASADI_SYMBOL_EXPORT const casadi_int* F_sparsity_in(casadi_int i) {
  struct matlab_external_global* g = &matlab_external_global_data;
  return g->in[i].compressed;
}

CASADI_SYMBOL_EXPORT const casadi_int* F_sparsity_out(casadi_int i) {
  struct matlab_external_global* g = &matlab_external_global_data;
  return g->out[i].compressed;
}


#ifdef __cplusplus
} /* extern "C" */
#endif
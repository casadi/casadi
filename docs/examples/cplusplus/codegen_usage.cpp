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


/**
 *  Example program demonstrating the usage of C-code generated from CasADi
 *  Generated code is encapsulated in self-contained code with the entry points
 *  defined below.
 *  Note that other usage, e.g. accessing the internal data structures in the
 *  generated files is not recommended and subject to change.
 *
 *  We show how the generated code can be used from C (or C++), without requiring
 *  any CasADi classes as well as how to use it from C++ using CasADi's external.
 *
 *  Joel Andersson, 2013-2015
 */

#include <stdio.h>
#include <dlfcn.h>

// Usage from C
int usage_c(){
  printf("---\n");
  printf("Standalone usage from C/C++:\n");
  printf("\n");

  /* Handle to the dll */
  void* handle;

  /* Load the dll */
  handle = dlopen("./f.so", RTLD_LAZY);
  if(handle==0){
    printf("Cannot open f.so, error %s\n", dlerror());
    return 1;
  }

  /* Reset error */
  dlerror();

  typedef long long int casadi_int;

  /* Typedefs */
  typedef void (*signal_t)(void);
  typedef casadi_int (*getint_t)(void);
  typedef int (*work_t)(casadi_int* sz_arg, casadi_int* sz_res, casadi_int* sz_iw, casadi_int* sz_w);
  typedef const casadi_int* (*sparsity_t)(casadi_int ind);
  typedef int (*eval_t)(const double** arg, double** res, casadi_int* iw, double* w, int mem);
  typedef int (*casadi_checkout_t)(void);
  typedef void (*casadi_release_t)(int);

  /* Memory management -- increase reference counter */
  signal_t incref = (signal_t)dlsym(handle, "f_incref");
  if(dlerror()) dlerror(); // No such function, reset error flags

  /* Memory management -- decrease reference counter */
  signal_t decref = (signal_t)dlsym(handle, "f_decref");
  if(dlerror()) dlerror(); // No such function, reset error flags

  /* Thread-local memory management -- checkout memory */
  casadi_checkout_t checkout = (casadi_checkout_t)dlsym(handle, "f_checkout");
  if(dlerror()) dlerror(); // No such function, reset error flags

  /* Thread-local memory management -- release memory */
  casadi_release_t release = (casadi_release_t)dlsym(handle, "f_release");
  if(dlerror()) dlerror(); // No such function, reset error flags

  /* Number of inputs */
  getint_t n_in_fcn = (getint_t)dlsym(handle, "f_n_in");
  if (dlerror()) return 1;
  casadi_int n_in = n_in_fcn();

  /* Number of outputs */
  getint_t n_out_fcn = (getint_t)dlsym(handle, "f_n_out");
  if (dlerror()) return 1;
  casadi_int n_out = n_out_fcn();

  /* Get sizes of the required work vectors */
  casadi_int sz_arg=n_in, sz_res=n_out, sz_iw=0, sz_w=0;
  work_t work = (work_t)dlsym(handle, "f_work");
  if(dlerror()) dlerror(); // No such function, reset error flags
  if (work && work(&sz_arg, &sz_res, &sz_iw, &sz_w)) return 1;
  printf("Work vector sizes:\n");
  printf("sz_arg = %lld, sz_res = %lld, sz_iw = %lld, sz_w = %lld\n\n",
         sz_arg, sz_res, sz_iw, sz_w);

  /* Input sparsities */
  sparsity_t sparsity_in = (sparsity_t)dlsym(handle, "f_sparsity_in");
  if (dlerror()) return 1;

  /* Output sparsities */
  sparsity_t sparsity_out = (sparsity_t)dlsym(handle, "f_sparsity_out");
  if (dlerror()) return 1;

  /* Print the sparsities of the inputs and outputs */
  casadi_int i;
  for(i=0; i<n_in + n_out; ++i){
    // Retrieve the sparsity pattern - CasADi uses column compressed storage (CCS)
    const casadi_int *sp_i;
    if (i<n_in) {
      printf("Input %lld\n", i);
      sp_i = sparsity_in(i);
    } else {
      printf("Output %lld\n", i-n_in);
      sp_i = sparsity_out(i-n_in);
    }
    if (sp_i==0) return 1;
    casadi_int nrow = *sp_i++; /* Number of rows */
    casadi_int ncol = *sp_i++; /* Number of columns */
    const casadi_int *colind = sp_i; /* Column offsets */
    const casadi_int *row = sp_i + ncol+1; /* Row nonzero */
    casadi_int nnz = sp_i[ncol]; /* Number of nonzeros */

    /* Print the pattern */
    printf("  Dimension: %lld-by-%lld (%lld nonzeros)\n", nrow, ncol, nnz);
    printf("  Nonzeros: {");
    casadi_int rr,cc,el;
    for(cc=0; cc<ncol; ++cc){                    /* loop over columns */
      for(el=colind[cc]; el<colind[cc+1]; ++el){ /* loop over the nonzeros entries of the column */
        if(el!=0) printf(", ");                  /* Separate the entries */
        rr = row[el];                            /* Get the row */
        printf("{%lld,%lld}",rr,cc);                 /* Print the nonzero */
      }
    }
    printf("}\n\n");
  }

  /* Function for numerical evaluation */
  eval_t eval = (eval_t)dlsym(handle, "f");
  if(dlerror()){
    printf("Failed to retrieve \"f\" function.\n");
    return 1;
  }

  /* Allocate input/output buffers and work vectors*/
  const double *arg[sz_arg];
  double *res[sz_res];
  casadi_int iw[sz_iw];
  double w[sz_w];

  /* Function input and output */
  const double x_val[] = {1,2,3,4};
  const double y_val = 5;
  double res0;
  double res1[4];

  // Allocate memory (thread-safe)
  incref();

  /* Evaluate the function */
  arg[0] = x_val;
  arg[1] = &y_val;
  res[0] = &res0;
  res[1] = res1;

  // Checkout thread-local memory (not thread-safe)
  // Note MAX_NUM_THREADS
  int mem = checkout();

  // Evaluation is thread-safe
  if (eval(arg, res, iw, w, mem)) return 1;

  // Release thread-local (not thread-safe)
  release(mem);

  /* Print result of evaluation */
  printf("result (0): %g\n",res0);
  printf("result (1): [%g,%g;%g,%g]\n",res1[0],res1[1],res1[2],res1[3]);

  /* Free memory (thread-safe) */
  decref();

  /* Free the handle */
  dlclose(handle);

  return 0;
}

// C++ (and CasADi) from here on
#include <casadi/casadi.hpp>
using namespace casadi;
using namespace std;

void usage_cplusplus(){
  cout << "---" << endl;
  cout << "Usage from CasADi C++:" << endl;
  cout << endl;

  // Use CasADi's "external" to load the compiled function
  Function f = external("f");

  // Use like any other CasADi function
  vector<double> x = {1, 2, 3, 4};
  vector<DM> arg = {reshape(DM(x), 2, 2), 5};
  vector<DM> res = f(arg);

  cout << "result (0): " << res.at(0) << endl;
  cout << "result (1): " << res.at(1) << endl;
}


// Usage from C using mem.h
#include <casadi/mem.h>
int usage_c_with_mem(){
  printf("---\n");
  printf("Usage from C/C++ with casadi/mem.h:\n");
  printf("\n");

  /* Handle to the dll */
  void* handle;

  /* Load the dll */
  handle = dlopen("./f_with_mem.so", RTLD_LAZY);
  if(handle==0){
    printf("Cannot open f_with_mem.so, error %s\n", dlerror());
    return 1;
  }

  /* Reset error */
  dlerror();

  /* Typedefs */
  typedef casadi_functions* (*functions_t)(void);

  /* mem.h interface */
  functions_t functions = (functions_t)dlsym(handle, "f_functions");
  if(dlerror()) dlerror(); // No such function, reset error flags

  casadi_functions* f = functions();

  casadi_mem* mem = casadi_alloc(f);

  /* Function input and output */
  const double x_val[] = {1,2,3,4};
  const double y_val = 5;
  double res0;
  double res1[4];

  /* Evaluate the function */
  mem->arg[0] = x_val;
  mem->arg[1] = &y_val;
  mem->res[0] = &res0;
  mem->res[1] = res1;

  casadi_eval(mem);

  /* Print result of evaluation */
  printf("result (0): %g\n",res0);
  printf("result (1): [%g,%g;%g,%g]\n",res1[0],res1[1],res1[2],res1[3]);

  casadi_free(mem);

  /* Success */
  return 0;
}


int main(){

  // Variables
  SX x = SX::sym("x", 2, 2);
  SX y = SX::sym("y");

  // Simple function
  Function f("f", {x, y}, {sqrt(y)-1, sin(x)-y});

  // Generate C-code
  f.generate("f");

  // Compile the C-code to a shared library
  string compile_command = "gcc -fPIC -shared -O3 f.c -o f.so";
  int flag = system(compile_command.c_str());
  casadi_assert(flag==0, "Compilation failed");

  // Usage from C
  flag = usage_c();
  casadi_assert(flag==0, "Example failed");

  // Usage from C++
  usage_cplusplus();

  // Generate C-code
  f.generate("f_with_mem", {{"with_mem", true}});

  // Compile the C-code to a shared library
  compile_command = "gcc -fPIC -I"+ std::string(INCLUDE_DIR) + " -shared -g f_with_mem.c -o f_with_mem.so";
  flag = system(compile_command.c_str());
  casadi_assert(flag==0, "Compilation failed");

  // Usage from c with mem.h
  usage_c_with_mem();

  return 0;
}

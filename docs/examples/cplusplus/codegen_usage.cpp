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
 *  any CasADi classes as well as how to use it from C++ using CasADi's Function::external.
 *
 *  Joel Andersson, 2013-2015
 */

#include <stdio.h>
#include <dlfcn.h>

// Usage from C
int usage_c(){
  printf("---\n");
  printf("Usage from C.\n");
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

  /* Initialize and get the number of inputs and outputs */
  typedef int (*init_t)(void **mem, int *n_in, int *n_out, void *data);
  init_t init = (init_t)dlsym(handle, "f_init");
  if(dlerror()){
    printf("Failed to retrieve \"init\" function.\n");
    return 1;
  }
  int n_in, n_out;
  void *mem;
  if (init(&mem, &n_in, &n_out, 0)) return 1;
  printf("n_in = %d, n_out = %d\n", n_in, n_out);

  /* Function for retrieving sparsities */
  typedef int (*sparsity_t)(void* mem, int ind, int *n_col, int *n_row,
                             const int** colind, const int** row);
  sparsity_t sparsity = (sparsity_t)dlsym(handle, "f_sparsity");
  if(dlerror()){
    printf("Failed to retrieve \"sparsity\" function.\n");
    return 1;
  }

  /* Get sparsities of the inputs and outptus */
  int ind;
  for(ind=0; ind<n_in + n_out; ++ind){
    if(ind<n_in){
      printf("Input %d\n",ind);
    } else {
      printf("Output %d\n",ind-n_in);
    }

    int nrow,ncol;
    const int *colind, *row;
    if (sparsity(mem, ind, &nrow, &ncol, &colind, &row)) return 1;

    printf("  Dimension: %d-by-%d\n", nrow, ncol);
    printf("  Nonzeros: {");
    int rr,cc,el;
    for(cc=0; cc<ncol; ++cc){ /* loop over rows */
      for(el=colind[cc]; el<colind[cc+1]; ++el){ /* loop over nonzeros */
        /* Separate the entries */
        if(el!=0) printf(", ");

        /* Get the column */
        rr = row[el]; 
        
        /* Print the nonzero */
        printf("{%d,%d}",rr,cc);
      }
    }
    printf("}\n\n");
  }

  /* Get sizes of the required work vectors */
  int sz_arg, sz_res, sz_iw, sz_w;
  typedef int (*work_t)(void* mem, int* sz_arg, int* sz_res, int* sz_iw, int* sz_w);
  work_t work = (work_t)dlsym(handle, "f_work");
  if(dlerror()){
    // No work vectors by default
    sz_arg = n_in;
    sz_res = n_out;
    sz_iw = sz_w = 0;
    // Reset error flags
    dlerror();
  } else {
    // Read from functions
    if (work(mem, &sz_arg, &sz_res, &sz_iw, &sz_w)) return 1;
  }
  printf("sz_arg = %d, sz_res = %d, sz_iw = %d, sz_w = %d\n", sz_arg, sz_res, sz_iw, sz_w);

  /* Function for numerical evaluation */
  typedef int (*eval_t)(void* mem, const double** arg, double** res, int* iw, double* w);
  eval_t eval = (eval_t)dlsym(handle, "f");
  if(dlerror()){
    printf("Failed to retrieve \"f\" function.\n");
    return 1;
  }

  /* Function for freeing memory */
  typedef int (*freemem_t)(void *mem);
  freemem_t freemem = (freemem_t)dlsym(handle, "f_free");
  if(dlerror()){
    // No free function
    freemem = 0;
    dlerror(); // Reset error flags
  }

  /* Allocate input/output buffers and work vectors*/
  const double *arg[sz_arg];
  double *res[sz_res];
  int iw[sz_iw];
  double w[sz_w];

  /* Function input and output */
  const double x_val[] = {1,2,3,4};
  const double y_val = 5;
  double res0;
  double res1[4];

  /* Evaluate the function */
  arg[0] = x_val;
  arg[1] = &y_val;
  res[0] = &res0;
  res[1] = res1;
  if (eval(mem, arg, res, iw, w)) return 1;

  /* Print result of evaluation */
  printf("result (0): %g\n",res0);
  printf("result (1): [%g,%g;%g,%g]\n",res1[0],res1[1],res1[2],res1[3]);

  /* Free memory */
  freemem(mem);

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
  cout << "Usage from C++" << endl;
  cout << endl;

  // Use CasADi's "Function::external" to load the compiled function
  Function f = Function::external("f");

  // Use like any other CasADi function
  vector<double> x = {1, 2, 3, 4};
  vector<DMatrix> arg = {reshape(DMatrix(x), 2, 2), 5};
  vector<DMatrix> res = f(arg);

  cout << "result (0): " << res.at(0) << endl;
  cout << "result (1): " << res.at(1) << endl;
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
  casadi_assert_message(flag==0, "Compilation failed");

  // Usage from C
  flag = usage_c();
  casadi_assert_message(flag==0, "Example failed");

  // Usage from C++
  usage_cplusplus();

  return 0;
}

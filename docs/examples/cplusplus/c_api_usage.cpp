/*
 *    MIT No Attribution
 *
 *    Copyright 2023 Joel Andersson, Joris Gillis, Moritz Diehl, KU Leuven.
 *
 *    Permission is hereby granted, free of charge, to any person obtaining a copy of this
 *    software and associated documentation files (the "Software"), to deal in the Software
 *    without restriction, including without limitation the rights to use, copy, modify,
 *    merge, publish, distribute, sublicense, and/or sell copies of the Software, and to
 *    permit persons to whom the Software is furnished to do so.
 *
 *    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
 *    INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
 *    PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 *    HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
 *    OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
 *    SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 *
 */

#include <stdio.h>

#include <casadi/casadi_c.h>

// Usage from C
int usage_c(){
  printf("---\n");
  printf("Standalone usage from C/C++:\n");
  printf("\n");

  // Sanity-check on integer type
  if (casadi_c_int_width()!=sizeof(casadi_int)) {
    printf("Mismatch in integer size\n");
    return -1;
  }
  if (casadi_c_real_width()!=sizeof(double)) {
    printf("Mismatch in double size\n");
    return -1;
  }

  // Push Function(s) to a stack
  int ret = casadi_c_push_file("f.casadi");
  if (ret) {
    printf("Failed to load file 'f.casadi'.\n");
    return -1;
  }
  ret = casadi_c_push_file("gh.casadi");
  if (ret) {
    printf("Failed to load file 'gh.casadi'.\n");
    return -1;
  }

  printf("Loaded number of functions: %d\n", casadi_c_n_loaded());

  // Identify a Function by name
  int id = casadi_c_id("g");


  casadi_int n_in = casadi_c_n_in_id(id);
  casadi_int n_out = casadi_c_n_out_id(id);

  casadi_int sz_arg=n_in, sz_res=n_out, sz_iw=0, sz_w=0;

  casadi_c_work_id(id, &sz_arg, &sz_res, &sz_iw, &sz_w);
  printf("Work vector sizes:\n");
  printf("sz_arg = %lld, sz_res = %lld, sz_iw = %lld, sz_w = %lld\n\n",
         sz_arg, sz_res, sz_iw, sz_w);

  /* Print the sparsities of the inputs and outputs */
  casadi_int i;
  for(i=0; i<n_in + n_out; ++i){
    // Retrieve the sparsity pattern - CasADi uses column compressed storage (CCS)
    const casadi_int *sp_i;
    if (i<n_in) {
      printf("Input %lld\n", i);
      sp_i = casadi_c_sparsity_in_id(id, i);
    } else {
      printf("Output %lld\n", i-n_in);
      sp_i = casadi_c_sparsity_out_id(id, i-n_in);
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
  casadi_c_incref_id(id);

  /* Evaluate the function */
  arg[0] = x_val;
  arg[1] = &y_val;
  res[0] = &res0;
  res[1] = res1;

  // Checkout thread-local memory (not thread-safe)
  int mem = casadi_c_checkout_id(id);

  // Evaluation is thread-safe
  if (casadi_c_eval_id(id, arg, res, iw, w, mem)) return 1;

  // Release thread-local (not thread-safe)
  casadi_c_release_id(id, mem);

  /* Print result of evaluation */
  printf("result (0): %g\n",res0);
  printf("result (1): [%g,%g;%g,%g]\n",res1[0],res1[1],res1[2],res1[3]);

  /* Free memory (thread-safe) */
  casadi_c_decref_id(id);

  // Clear the last loaded Function(s) from the stack
  casadi_c_pop();

  return 0;
}

// C++ (and CasADi) from here on
#include <casadi/casadi.hpp>
using namespace casadi;

int main(){

  // Variables
  SX x = SX::sym("x", 2, 2);
  SX y = SX::sym("y");

  // Simple function
  Function f("f", {x, y}, {x*y});

  // Mode 1: Function::save
  f.save("f.casadi");

  // More simple functions
  Function g("g", {x, y}, {sqrt(y)-1, sin(x)-y});
  Function h("h", {x}, {y*y});

  // Mode 2: FileSerializer (allows packing a list of Functions)
  {
    FileSerializer gh("gh.casadi",{{"debug", true}});
    gh.pack(std::vector<Function>{g, h});
  }

  // Usage from C
  usage_c();

  return 0;
}

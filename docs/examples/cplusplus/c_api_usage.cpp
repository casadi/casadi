/*
 *    MIT No Attribution
 *
 *    Copyright (C) 2010-2023 Joel Andersson, Joris Gillis, Moritz Diehl, KU Leuven.
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
#include <math.h>

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

// Solver stats from C
int usage_c_stats(){
  printf("---\n");
  printf("Stats from C:\n");
  printf("\n");

  if (casadi_c_push_file("solver.casadi")) {
    printf("No 'solver.casadi' (Ipopt not available): skipped.\n");
    return 0;
  }
  int id = casadi_c_id("solver");

  casadi_int sz_arg, sz_res, sz_iw, sz_w;
  casadi_c_work_id(id, &sz_arg, &sz_res, &sz_iw, &sz_w);
  const double *arg[sz_arg];
  double *res[sz_res];
  casadi_int iw[sz_iw];
  double w[sz_w];
  casadi_int i;
  for (i=0; i<sz_arg; ++i) arg[i] = 0;
  for (i=0; i<sz_res; ++i) res[i] = 0;

  /* Rosenbrock; null inputs are zero, not their default: pass the bounds */
  const double inf = HUGE_VAL;
  const double x0[2] = {-1, 1}, p = 100, lbx[2] = {-inf, -inf}, ubx[2] = {inf, inf};
  const double lbg = -inf, ubg = 1;
  double x[2];
  arg[0] = x0;
  arg[1] = &p;
  arg[2] = lbx;
  arg[3] = ubx;
  arg[4] = &lbg;
  arg[5] = &ubg;
  res[0] = x;

  /* Sink on a caller-owned buffer */
  static unsigned char buf[1 << 16];
  struct casadi_stats_sink s = casadi_c_stats_make_sink(buf, sizeof(buf));

  casadi_c_incref_id(id);
  int mem = casadi_c_checkout_id(id);
  int k;
  for (k=0; k<2; ++k) {
    /* -1: root call; both calls are recorded */
    if (casadi_c_eval_with_stats_id(id, arg, res, iw, w, mem, &s, -1)) return 1;
  }
  casadi_c_release_id(id, mem);

  /* Optional, after the last write: queries skip whole subtrees */
  casadi_c_stats_index(&s);

  /* Nodes are stream offsets, -1 for none or the root; a cursor continues after a node.
     First call of solver at any depth; "#5:solver" would match the full id */
  casadi_int c = casadi_c_stats_find_function(&s, -1, -1, "solver");
  char sid[64];
  casadi_int iter_count;
  if (casadi_c_stats_get_text(&s, c, "id", sid, sizeof(sid), 0)) return 1;
  if (casadi_c_stats_get_int(&s, c, "iter_count", &iter_count)) return 1;
  printf("%s: %lld iterations, x = [%g, %g] (%lld bytes of stats)\n",
         sid, iter_count, x[0], x[1], casadi_c_stats_nbytes(&s));

  /* Iteration 5: objective, and the flag of the nlp_jac_g call made in it */
  casadi_int it = casadi_c_stats_select_iteration(&s, c, -1, 5);
  casadi_int g = casadi_c_stats_select_function(&s, it, -1, "nlp_jac_g");
  double obj;
  casadi_int flag;
  if (casadi_c_stats_get_real(&s, it, "obj", &obj)) return 1;
  if (casadi_c_stats_get_int(&s, g, "flag", &flag)) return 1;
  printf("iteration 5: obj = %g, nlp_jac_g flag %lld\n", obj, flag);

  /* Next match: the second call; inf_pr of its last iteration */
  c = casadi_c_stats_find_function(&s, -1, c, "solver");
  double inf_pr;
  if (casadi_c_stats_get_int(&s, c, "iter_count", &iter_count)) return 1;
  it = casadi_c_stats_select_iteration(&s, c, -1, iter_count);
  if (casadi_c_stats_get_real(&s, it, "inf_pr", &inf_pr)) return 1;
  printf("second call: inf_pr = %g\n", inf_pr);

  /* Calls after the iterations are in section post */
  casadi_int post = casadi_c_stats_select_section(&s, c, -1, "post");
  printf("post: %s\n", casadi_c_stats_select_function(&s, post, -1, 0) >= 0 ? "calls" : "none");

  casadi_c_decref_id(id);
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
  Function h("h", {y}, {y*y});

  // Mode 2: FileSerializer (allows packing a list of Functions)
  {
    FileSerializer gh("gh.casadi",{{"debug", true}});
    gh.pack(std::vector<Function>{g, h});
  }

  if (has_nlpsol("ipopt")) {
    MX xs = MX::sym("x", 2);
    MX ps = MX::sym("p");
    MX x1 = xs(0), x2 = xs(1);
    MXDict nlp = {{"x", xs}, {"p", ps}, {"g", x1 + x2}, {"f", sq(1 - x1) + ps * sq(x2 - sq(x1))}};
    Function solver = nlpsol("solver", "ipopt", nlp,
      {{"ipopt.print_level", 0}, {"ipopt.sb", "yes"}, {"print_time", false}});
    solver.save("solver.casadi");
  }

  // Usage from C
  usage_c();

  return usage_c_stats();
}

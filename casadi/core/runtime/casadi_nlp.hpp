//
//    MIT No Attribution
//
//    Copyright (C) 2010-2023 Joel Andersson, Joris Gillis, Moritz Diehl, KU Leuven.
//
//    Permission is hereby granted, free of charge, to any person obtaining a copy of this
//    software and associated documentation files (the "Software"), to deal in the Software
//    without restriction, including without limitation the rights to use, copy, modify,
//    merge, publish, distribute, sublicense, and/or sell copies of the Software, and to
//    permit persons to whom the Software is furnished to do so.
//
//    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
//    INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
//    PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
//    HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
//    OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
//    SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
//


// SYMBOL "nlpsol_detect_bounds_prob"
template<typename T1>
struct casadi_nlpsol_detect_bounds_prob {
  casadi_int sz_arg;
  casadi_int sz_res;
  casadi_int sz_iw;
  casadi_int sz_w;
  // Original number og constraints
  casadi_int ng;
  // Number of bounds
  casadi_int nb;
  casadi_int *target_x;
  casadi_int *target_g;
  char *is_simple;
};
// C-REPLACE "casadi_nlpsol_detect_bounds_prob<T1>" "struct casadi_nlpsol_detect_bounds_prob"

// SYMBOL "nlpsol_prob"
template<typename T1>
struct casadi_nlpsol_prob {
  casadi_int nx, ng, np;

  casadi_nlpsol_detect_bounds_prob<T1> detect_bounds;
};
// C-REPLACE "casadi_nlpsol_prob<T1>" "struct casadi_nlpsol_prob"

// SYMBOL "nlpsol_detect_bounds_data"
template<typename T1>
struct casadi_nlpsol_detect_bounds_data {
  // Work vectors
  const double** arg;
  double** res;
  casadi_int* iw;
  double* w;

  // Simple bounds g(x)=A(b)x+b(p);
  // a[i]*x[target_x[i]]+b[i]
  T1* a;
  T1* b;
  casadi_int* target_l;
  casadi_int* target_u;
  T1* lam_xl;
  T1* lam_xu;
};
// C-REPLACE "casadi_nlpsol_detect_bounds_data<T1>" "struct casadi_nlpsol_detect_bounds_data"

// SYMBOL "nlpsol_data"
template<typename T1>
struct casadi_nlpsol_data {
  // Problem structure
  const casadi_nlpsol_prob<T1>* prob;
  // Variable bounds
  T1 *lbz, *ubz;
  // Current primal solution
  T1 *z;
  // Current dual solution
  T1 *lam;
  // Storage for f when null pointer passed
  T1 objective;

  // NLP data, pointers to arg (no allocations needed)
  const T1 *p, *lbx, *ubx, *lbg, *ubg, *x0, *lam_x0, *lam_g0;
  // NLP results, pointers to res (no allocations needed)
  T1 *f, *x, *g, *lam_x, *lam_g, *lam_p;

  casadi_nlpsol_detect_bounds_data<T1> detect_bounds;
};
// C-REPLACE "casadi_nlpsol_data<T1>" "struct casadi_nlpsol_data"

// SYMBOL "nlpsol_work"
template<typename T1>
void casadi_nlpsol_work(const casadi_nlpsol_prob<T1>* p, casadi_int* sz_arg, casadi_int* sz_res,
    casadi_int* sz_iw, casadi_int* sz_w) {
  // Reset sz_arg, sz_res
  *sz_arg = *sz_res = 0;
  // Reset sz_w, sz_iw
  *sz_w = *sz_iw = 0;
  *sz_w += p->nx + p->ng; // z
  *sz_w += p->nx + p->ng; // lbz
  *sz_w += p->nx + p->ng; // ubz
  *sz_w += p->nx + p->ng; // lam

  if (p->detect_bounds.ng) {
    *sz_arg += p->detect_bounds.sz_arg;
    *sz_res += p->detect_bounds.sz_res;
    *sz_iw += p->detect_bounds.sz_iw;
    *sz_w += p->detect_bounds.sz_w;

    *sz_w += p->detect_bounds.nb; // a;
    *sz_w += p->detect_bounds.nb; // b;
    *sz_iw += p->nx; // target_l
    *sz_iw += p->nx; // target_u
    *sz_w += p->nx; // lam_xl;
    *sz_w += p->nx; // lam_xu;
  }

}


// SYMBOL "nlpsol_init"
template<typename T1>
void casadi_nlpsol_init(casadi_nlpsol_data<T1>* d, const T1*** arg, T1*** res,
    casadi_int** iw, T1** w) {
  // Local variables
  casadi_int nx, ng;
  const casadi_nlpsol_prob<T1>* p = d->prob;
  nx = p->nx;
  ng = p->ng;
  // Get matrix number of nonzeros
  d->z = *w; *w += nx + ng;
  d->lbz = *w; *w += nx + ng;
  d->ubz = *w; *w += nx + ng;
  d->lam = *w; *w += nx + ng;

  if (p->detect_bounds.ng) {
    d->detect_bounds.arg = *arg; *arg += p->detect_bounds.sz_arg;
    d->detect_bounds.res = *res; *res += p->detect_bounds.sz_res;
    d->detect_bounds.iw = *iw; *iw += p->detect_bounds.sz_iw;
    d->detect_bounds.w = *w; *w += p->detect_bounds.sz_w;

    d->detect_bounds.a = *w; *w += p->detect_bounds.nb;
    d->detect_bounds.b = *w; *w += p->detect_bounds.nb;
    d->detect_bounds.target_l = *iw; *iw += p->nx;
    d->detect_bounds.target_u = *iw; *iw += p->nx;
    d->detect_bounds.lam_xl = *w; *w += nx;
    d->detect_bounds.lam_xu = *w; *w += nx;
  }

}

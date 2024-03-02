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


// C-REPLACE "casadi_nlpsol_prob<T1>" "struct casadi_nlpsol_prob"

// SYMBOL "sqpmethod_prob"
template<typename T1>
struct casadi_sqpmethod_prob {
  const casadi_nlpsol_prob<T1>* nlp;
  // Sparsity patterns
  const casadi_int *sp_h, *sp_a, *sp_hr;
  casadi_int merit_memsize;
  casadi_int max_iter_ls;
};
// C-REPLACE "casadi_sqpmethod_prob<T1>" "struct casadi_sqpmethod_prob"


// SYMBOL "sqpmethod_data"
template<typename T1>
struct casadi_sqpmethod_data {
  // Problem structure
  const casadi_sqpmethod_prob<T1>* prob;

  T1* z_cand;
  // Lagrange gradient in the next iterate
  T1 *gLag, *gLag_old;
  // Gradient of the objective
  T1 *gf;
  // Bounds of the QP
  T1 *lbdz, *ubdz;
  // QP solution
  T1 *dx, *dlam;
  // Hessian approximation
  T1 *Bk;
  // Jacobian
  T1* Jk;
  // merit_mem
  T1* merit_mem;
  // temp_mem
  T1* temp_mem;
  // temp_sol
  T1* temp_sol;

  const T1** arg;
  T1** res;
  casadi_int* iw;
  T1* w;
};
// C-REPLACE "casadi_sqpmethod_data<T1>" "struct casadi_sqpmethod_data"


// SYMBOL "sqpmethod_work"
template<typename T1>
void casadi_sqpmethod_work(const casadi_sqpmethod_prob<T1>* p,
    casadi_int* sz_iw, casadi_int* sz_w, int elastic_mode, int so_corr) {
  // Local variables
  casadi_int nnz_h, nnz_a, nx, ng;
  nnz_h = p->sp_h[2+p->sp_h[1]];
  nnz_a = p->sp_a[2+p->sp_a[1]];
  nx = p->nlp->nx;
  ng = p->nlp->ng;

  // Reset sz_w, sz_iw
  *sz_w = *sz_iw = 0;
  if (p->max_iter_ls>0 || so_corr) *sz_w += nx + ng; // z_cand
  // Lagrange gradient in the next iterate
  *sz_w += nx; // gLag
  *sz_w += nx; // gLag_old
  // Gradient of the objective
  *sz_w += nx; // gf
  // Bounds of the QP
  *sz_w += nx + ng; // lbdz
  *sz_w += nx + ng; // ubdz
  // QP solution
  *sz_w += nx; // dx
  *sz_w += nx + ng; // dlam
  // Hessian approximation
  *sz_w += nnz_h; // Bk
  // Jacobian
  *sz_w += nnz_a; // Jk
  // merit_mem
  if (p->max_iter_ls>0 || so_corr) *sz_w += p->merit_memsize;

  if (elastic_mode) {
    // Additional work for larger objective gradient
    *sz_w += 2*ng; // gf
    // Additional work for the larger bounds
    *sz_w += 2*ng; // lbdz
    *sz_w += 2*ng; // ubdz
    // Additional work for larger solution
    *sz_w += 2*ng; // dx
    *sz_w += 2*ng; // dlam
    // Additional work for larger jacobian
    *sz_w += 2*ng; // Jk
    // Additional work for temp memory
    *sz_w += ng;
  }

  if (so_corr) *sz_w += nx+nx+ng; // Temp memory for failing soc
}

// SYMBOL "sqpmethod_init"
template<typename T1>
void casadi_sqpmethod_init(casadi_sqpmethod_data<T1>* d,
    const T1*** arg, T1*** res, casadi_int** iw, T1** w,
    int elastic_mode, int so_corr) {
  // Local variables
  casadi_int nnz_h, nnz_a, nx, ng;
  const casadi_sqpmethod_prob<T1>* p = d->prob;
  // Get matrix number of nonzeros
  nnz_h = p->sp_h[2+p->sp_h[1]];
  nnz_a = p->sp_a[2+p->sp_a[1]];
  nx = p->nlp->nx;
  ng = p->nlp->ng;
  if (p->max_iter_ls>0 || so_corr) {
    d->z_cand = *w; *w += nx + ng;
  }
  // Lagrange gradient in the next iterate
  d->gLag = *w; *w += nx;
  d->gLag_old = *w; *w += nx;
  // Hessian approximation
  d->Bk = *w; *w += nnz_h;
  // merit_mem
  if (p->max_iter_ls>0 || so_corr) {
    d->merit_mem = *w; *w += p->merit_memsize;
  }

  if (so_corr) {
    d->temp_sol = *w; *w += nx+nx+ng;
  }

  if (elastic_mode) {
    // Gradient of the objective
    d->gf = *w; *w += nx + 2*ng;
    // Bounds of the QP
    d->lbdz = *w; *w += nx + 3*ng;
    d->ubdz = *w; *w += nx + 3*ng;
    // QP solution
    d->dx = *w; *w += nx + 2*ng;
    d->dlam = *w; *w += nx + 3*ng;
    // Jacobian
    d->Jk = *w; *w += nnz_a + 2*ng;
    // temp mem
    d->temp_mem = *w; *w += ng;
  } else {
    // Gradient of the objective
    d->gf = *w; *w += nx;
    // Bounds of the QP
    d->lbdz = *w; *w += nx + ng;
    d->ubdz = *w; *w += nx + ng;
    // QP solution
    d->dx = *w; *w += nx;
    d->dlam = *w; *w += nx + ng;
    // Jacobian
    d->Jk = *w; *w += nnz_a;
  }
  d->arg = *arg;
  d->res = *res;
  d->iw = *iw;
  d->w = *w;
}

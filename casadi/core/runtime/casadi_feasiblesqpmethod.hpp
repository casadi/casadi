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

// SYMBOL "feasiblesqpmethod_prob"
template<typename T1>
struct casadi_feasiblesqpmethod_prob {
  const casadi_nlpsol_prob<T1>* nlp;
  // Sparsity patterns
  const casadi_int *sp_h, *sp_a, *sp_hr;
  // casadi_int merit_memsize;
  // casadi_int max_iter_ls;
};
// C-REPLACE "casadi_feasiblesqpmethod_prob<T1>" "struct casadi_feasiblesqpmethod_prob"


// SYMBOL "feasiblesqpmethod_data"
template<typename T1>
struct casadi_feasiblesqpmethod_data {
  // Problem structure
  const casadi_feasiblesqpmethod_prob<T1>* prob;

  // Lagrange gradient in the next iterate
  T1 *gLag, *gLag_old;
  // Gradient of the objective
  T1 *gf;
  // Bounds of the QP
  T1 *lbdz, *ubdz;
  // QP solution
  T1 *dx, *dlam;
  // Feasibility QP solution
  T1 *dx_feas;
  T1 *dlam_feas;
  // Feasibility iterate
  T1 *z_feas;
  T1 *gf_feas;
  T1 *lbdz_feas, *ubdz_feas;
  T1* z_tmp;
  T1 *tr_scale_vector;
  casadi_int *tr_mask;
  // Hessian approximation
  T1 *Bk;
  // Jacobian
  T1* Jk;
  // Anderson vectors
  T1* anderson_memory_step;
  T1* anderson_memory_iterate;
  T1* gamma;

  // Function value of feasibility iterate
  T1 f_feas;
  // merit_mem
  // T1* merit_mem;
  // temp_mem
  // T1* temp_mem;
  // temp_sol
  // T1* temp_sol;
};
// C-REPLACE "casadi_feasiblesqpmethod_data<T1>" "struct casadi_feasiblesqpmethod_data"


// SYMBOL "feasiblesqpmethod_work"
template<typename T1>
void casadi_feasiblesqpmethod_work(const casadi_feasiblesqpmethod_prob<T1>* p,
    casadi_int* sz_iw, casadi_int* sz_w, int sz_anderson_memory) {
  // Local variables
  casadi_int nnz_h, nnz_a, nx, ng;
  nnz_h = p->sp_h[2+p->sp_h[1]];
  nnz_a = p->sp_a[2+p->sp_a[1]];
  nx = p->nlp->nx;
  ng = p->nlp->ng;

  // Reset sz_w, sz_iw
  *sz_w = *sz_iw = 0;
  // if (p->max_iter_ls>0) *sz_w += nx + ng; // z_cand
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
  // Feasibility QP solution
  *sz_w += nx; // dx_feas
  *sz_w += nx + ng; // dlam_feas
  // Feasibility iterate
  *sz_w += nx + ng; // x_feas + g_feas
  *sz_w += nx; // gf_feas
  *sz_w += nx + ng; // lower bounds feasibile QP
  *sz_w += nx + ng; // upper bounds feasible QP
  *sz_w += nx+ng; // x tmp feasible QP
  *sz_w += nx; // tr_scale_vector
  *sz_iw += nx; // tr_mask
  // Hessian approximation
  *sz_w += nnz_h; // Bk
  // Jacobian
  *sz_w += nnz_a; // Jk
  // merit_mem
  // if (p->max_iter_ls>0) *sz_w += p->merit_memsize;

  if (sz_anderson_memory > 0) {
    // for step (dx)
    *sz_w += sz_anderson_memory*nx;
    // for x
    *sz_w += sz_anderson_memory*nx;
    // for gamma
    *sz_w += sz_anderson_memory;
  }

  // if (elastic_mode) {
  //   // Additional work for larger objective gradient
  //   *sz_w += 2*ng; // gf
  //   // Additional work for the larger bounds
  //   *sz_w += 2*ng; // lbdz
  //   *sz_w += 2*ng; // ubdz
  //   // Additional work for larger solution
  //   *sz_w += 2*ng; // dx
  //   *sz_w += 2*ng; // dlam
  //   // Additional work for larger jacobian
  //   *sz_w += 2*ng; // Jk
  //   // Additional work for temp memory
  //   *sz_w += ng;
  // }

  // if (so_corr) *sz_w += nx+nx+ng; // Temp memory for failing soc
}

// SYMBOL "feasiblesqpmethod_init"
template<typename T1>
void casadi_feasiblesqpmethod_init(casadi_feasiblesqpmethod_data<T1>* d,
    casadi_int** iw, T1** w, int sz_anderson_memory) {
  // Local variables
  casadi_int nnz_h, nnz_a, nx, ng;
  const casadi_feasiblesqpmethod_prob<T1>* p = d->prob;
  // Get matrix number of nonzeros
  nnz_h = p->sp_h[2+p->sp_h[1]];
  nnz_a = p->sp_a[2+p->sp_a[1]];
  nx = p->nlp->nx;
  ng = p->nlp->ng;
  // if (p->max_iter_ls>0) {
  //   d->z_cand = *w; *w += nx + ng;
  // }
  // Lagrange gradient in the next iterate
  d->gLag = *w; *w += nx;
  d->gLag_old = *w; *w += nx;
  // Hessian approximation
  d->Bk = *w; *w += nnz_h;
  // merit_mem
  // if (p->max_iter_ls>0) {
  //   d->merit_mem = *w; *w += p->merit_memsize;
  // }

  // if (so_corr) {
  //   d->temp_sol = *w; *w += nx+nx+ng;
  // }

  // Gradient of the objective
  d->gf = *w; *w += nx;
  // Bounds of the QP
  d->lbdz = *w; *w += nx + ng;
  d->ubdz = *w; *w += nx + ng;
  // QP solution
  d->dx = *w; *w += nx;
  d->dlam = *w; *w += nx + ng;
  // Feasible QP solution
  d->dx_feas = *w; *w += nx;
  d->dlam_feas = *w; *w += nx + ng;
  // feasibility iterate
  d->z_feas = *w; *w += nx + ng;
  d->gf_feas = *w; *w += nx;
  // Bounds of the feasibility QPs
  d->lbdz_feas = *w; *w += nx + ng;
  d->ubdz_feas = *w; *w += nx + ng;
  // x tmp for QPs
  d->z_tmp = *w; *w += sz_anderson_memory*nx+ng;
  // trust-region scale vector
  d->tr_scale_vector = *w; *w += nx;
  d->tr_mask = *iw; *iw += nx;
  // Jacobian
  d->Jk = *w; *w += nnz_a;
  // Anderson vector
  d->anderson_memory_step = *w; *w += sz_anderson_memory*nx;
  d->anderson_memory_iterate = *w; *w += sz_anderson_memory*nx;
  d->gamma = *w; *w += sz_anderson_memory;

}

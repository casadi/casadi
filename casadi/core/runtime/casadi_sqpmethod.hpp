// NOLINT(legal/copyright)

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
};
// C-REPLACE "casadi_sqpmethod_data<T1>" "struct casadi_sqpmethod_data"


// SYMBOL "sqpmethod_work"
template<typename T1>
void casadi_sqpmethod_work(const casadi_sqpmethod_prob<T1>* p,
    casadi_int* sz_iw, casadi_int* sz_w) {
  // Local variables
  casadi_int nnz_h, nnz_a, nx, ng;
  nnz_h = p->sp_h[2+p->sp_h[1]];
  nnz_a = p->sp_a[2+p->sp_a[1]];
  nx = p->nlp->nx;
  ng = p->nlp->ng;

  // Reset sz_w, sz_iw
  *sz_w = *sz_iw = 0;
  if (p->max_iter_ls>0) *sz_w += nx + ng; // z_cand
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
  if (p->max_iter_ls>0) *sz_w += p->merit_memsize;
}

// SYMBOL "sqpmethod_init"
template<typename T1>
void casadi_sqpmethod_init(casadi_sqpmethod_data<T1>* d, casadi_int** iw, T1** w) {
  // Local variables
  casadi_int nnz_h, nnz_a, nx, ng;
  const casadi_sqpmethod_prob<T1>* p = d->prob;
  // Get matrix number of nonzeros
  nnz_h = p->sp_h[2+p->sp_h[1]];
  nnz_a = p->sp_a[2+p->sp_a[1]];
  nx = p->nlp->nx;
  ng = p->nlp->ng;
  if (p->max_iter_ls>0) {
    d->z_cand = *w;
    *w += nx + ng;
  }
  // Lagrange gradient in the next iterate
  d->gLag = *w; *w += nx;
  d->gLag_old = *w; *w += nx;
  // Gradient of the objective
  d->gf = *w; *w += nx;
  // Bounds of the QP
  d->lbdz = *w; *w += nx + ng;
  d->ubdz = *w; *w += nx + ng;
  // QP solution
  d->dx = *w; *w += nx;
  d->dlam = *w; *w += nx + ng;
  // Hessian approximation
  d->Bk = *w; *w += nnz_h;
  // Jacobian
  d->Jk = *w; *w += nnz_a;
  // merit_mem
  if (p->max_iter_ls>0) {
    d->merit_mem = *w;
    *w += p->merit_memsize;
  }
}

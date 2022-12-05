// NOLINT(legal/copyright)

// SYMBOL "nlpsol_prob"
template<typename T1>
struct casadi_nlpsol_prob {
  casadi_int nx, ng, np;
};
// C-REPLACE "casadi_nlpsol_prob<T1>" "struct casadi_nlpsol_prob"

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
};
// C-REPLACE "casadi_nlpsol_data<T1>" "struct casadi_nlpsol_data"

// SYMBOL "nlpsol_work"
template<typename T1>
void casadi_nlpsol_work(const casadi_nlpsol_prob<T1>* p, casadi_int* sz_iw, casadi_int* sz_w) {
  // Reset sz_w, sz_iw
  *sz_w = *sz_iw = 0;
  *sz_w += p->nx + p->ng; // z
  *sz_w += p->nx + p->ng; // lbz
  *sz_w += p->nx + p->ng; // ubz
  *sz_w += p->nx + p->ng; // lam
}


// SYMBOL "nlpsol_init"
template<typename T1>
void casadi_nlpsol_init(casadi_nlpsol_data<T1>* d, casadi_int** iw, T1** w) {
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
}

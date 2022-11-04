// NOLINT(legal/copyright)

// C-REPLACE "fmin" "casadi_fmin"
// C-REPLACE "fmax" "casadi_fmax"
// C-REPLACE "std::numeric_limits<T1>::min()" "casadi_real_min"
// C-REPLACE "std::numeric_limits<T1>::infinity()" "casadi_inf"
// C-REPLACE "static_cast<int>" "(int) "

// SYMBOL "qp_prob"
template<typename T1>
struct casadi_qp_prob {
  // Sparsity patterns
  const casadi_int *sp_a, *sp_h;
  // Dimensions
  casadi_int nx, na, nz;
};
// C-REPLACE "casadi_qp_prob<T1>" "struct casadi_qp_prob"

// SYMBOL "qp_setup"
template<typename T1>
void casadi_qp_setup(casadi_qp_prob<T1>* p) {
  p->na = p->sp_a[0];
  p->nx = p->sp_a[1];
  p->nz = p->na+p->nx;
}

// SYMBOL "qp_data"
template<typename T1>
struct casadi_qp_data {
  // Problem structure
  const casadi_qp_prob<T1>* prob;

  // QP data, pointers to arg (no allocations needed)
  const T1 *a, *h, *g, *lbx, *ubx, *lba, *uba, *x0, *lam_x0, *lam_a0;
  // QP results, pointers to res (no allocations needed)
  T1 *f, *x, *lam_x, *lam_a;

};
// C-REPLACE "casadi_qp_data<T1>" "struct casadi_qp_data"

// SYMBOL "qp_init"
template<typename T1>
void casadi_qp_init(casadi_qp_data<T1>* d, casadi_int** iw, T1** w) {
  // Local variables
  //const casadi_qp_prob<T1>* p = d->prob;
}

// SYMBOL "qp_work"
template<typename T1>
void casadi_qp_work(const casadi_qp_prob<T1>* p, casadi_int* sz_iw, casadi_int* sz_w) {
  // Reset sz_w, sz_iw
  *sz_w = *sz_iw = 0;
}

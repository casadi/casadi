// NOLINT(legal/copyright)

// C-REPLACE "fmin" "casadi_fmin"
// C-REPLACE "fmax" "casadi_fmax"
// C-REPLACE "std::numeric_limits<T1>::min()" "casadi_real_min"
// C-REPLACE "std::numeric_limits<T1>::infinity()" "casadi_inf"
// C-REPLACE "static_cast<int>" "(int) "
// SYMBOL "ipqp_prob"
template<typename T1>
struct casadi_ipqp_prob {
  // Sparsity patterns
  const casadi_int *sp_a, *sp_h, *sp_at, *sp_kkt;
  // Symbolic QR factorization
  const casadi_int *prinv, *pc, *sp_v, *sp_r;
  // Dimensions
  casadi_int nx, na, nz;
  // Smallest nonzero number
  T1 dmin;
  // Infinity
  T1 inf;
  // Smallest multiplier treated as inactive for the initial active set
  T1 min_lam;
  // Maximum number of iterations
  casadi_int max_iter;
  // Primal and dual error tolerance
  T1 constr_viol_tol, dual_inf_tol;
};
// C-REPLACE "casadi_ipqp_prob<T1>" "struct casadi_ipqp_prob"

// SYMBOL "ipqp_setup"
template<typename T1>
void casadi_ipqp_setup(casadi_ipqp_prob<T1>* p) {
  p->na = p->sp_a[0];
  p->nx = p->sp_a[1];
  p->nz = p->nx + p->na;
  p->dmin = std::numeric_limits<T1>::min();
  p->inf = std::numeric_limits<T1>::infinity();
  p->min_lam = 0;
  p->max_iter = 1000;
  p->constr_viol_tol = 1e-8;
  p->dual_inf_tol = 1e-8;
}

// SYMBOL "ipqp_work"
template<typename T1>
void casadi_ipqp_work(const casadi_ipqp_prob<T1>* p, casadi_int* sz_iw, casadi_int* sz_w) {
  // Local variables
  casadi_int nnz_a, nnz_kkt, nnz_v, nnz_r;
  // Get matrix number of nonzeros
  nnz_a = p->sp_a[2+p->sp_a[1]];
  nnz_kkt = p->sp_kkt[2+p->sp_kkt[1]];
  nnz_v = p->sp_v[2+p->sp_v[1]];
  nnz_r = p->sp_r[2+p->sp_r[1]];
  // Reset sz_w, sz_iw
  *sz_w = *sz_iw = 0;
  // Temporary work vectors
  *sz_w = casadi_max(*sz_w, p->nz); // casadi_project, tau memory
  *sz_iw = casadi_max(*sz_iw, p->nz); // casadi_trans, tau type, allzero
  *sz_w = casadi_max(*sz_w, 2*p->nz); // casadi_qr
  // Persistent work vectors
  *sz_w += nnz_kkt; // kkt
  *sz_w += p->nz; // z=[xk,gk]
  *sz_w += p->nz; // lbz
  *sz_w += p->nz; // ubz
  *sz_w += p->nz; // lam
  *sz_w += nnz_a; // trans(a)
  *sz_w += p->nz; // dz
  *sz_w += p->nz; // dlam
  *sz_w += p->nx; // infeas
  *sz_w += p->nx; // tinfeas
  *sz_w += p->nz; // sens
  *sz_iw += p->nz; // neverzero
  *sz_iw += p->nz; // neverupper
  *sz_iw += p->nz; // neverlower
  *sz_iw += p->nz; // lincomb
  *sz_w += casadi_max(nnz_v+nnz_r, nnz_kkt); // [v,r] or trans(kkt)
  *sz_w += p->nz; // beta
  *sz_w += p->nz; // D
  *sz_w += p->nz; // S
  *sz_w += p->nz; // lam_lbz
  *sz_w += p->nz; // lam_ubz
  *sz_w += p->nz; // dlam_lbz
  *sz_w += p->nz; // dlam_ubz
  *sz_w += p->nz; // rz
  *sz_w += p->nz; // rlam
  *sz_w += p->nz; // rlam_lbz
  *sz_w += p->nz; // rlam_ubz
  *sz_w += p->nz; // dinv_lbz
  *sz_w += p->nz; // dinv_ubz
}

// SYMBOL "ipqp_flag_t"
typedef enum {
  IPQP_SUCCESS,
  IPQP_MAX_ITER,
  IPQP_NO_SEARCH_DIR,
  IPQP_PRINTING_ERROR
} casadi_ipqp_flag_t;

// SYMBOL "ipqp_task_t"
typedef enum {
  IPQP_MV,
  IPQP_PROGRESS,
  IPQP_FACTOR,
  IPQP_SOLVE} casadi_ipqp_task_t;

// SYMBOL "ipqp_task_t"
typedef enum {
  IPQP_INIT,
  IPQP_RESIDUAL,
  IPQP_NEWITER,
  IPQP_PREPARE,
  IPQP_PREDICTOR,
  IPQP_CORRECTOR} casadi_ipqp_next_t;

// SYMBOL "ipqp_blocker_t"
typedef enum {
  IPQP_NONE = 0x0,
  IPQP_UPPER = 0x1,
  IPQP_LOWER = 0x2,
  IPQP_PRIMAL = 0x4,
  IPQP_DUAL = 0x8
} casadi_blocker_t;

// SYMBOL "ipqp_data"
template<typename T1>
struct casadi_ipqp_data {
  // Problem structure
  const casadi_ipqp_prob<T1>* prob;
  // Solver status
  casadi_ipqp_flag_t status;
  // Cost
  T1 f;
  // QP data
  const T1 *nz_a, *nz_h, *g;
  // Vectors
  T1 *z, *lbz, *ubz, *infeas, *tinfeas, *sens, *lam, *w, *dz, *dlam;
  casadi_int *iw, *neverzero, *neverlower, *neverupper, *lincomb;
  // Numeric QR factorization
  T1 *nz_at, *nz_kkt, *beta, *nz_v, *nz_r;
  // Message buffer
  const char *msg;
  // Message index
  casadi_int msg_ind;
  // Stepsize
  T1 tau;
  // Singularity
  casadi_int sing;
  // Do we already have a search direction?
  int has_search_dir;
  // Smallest diagonal value for the QR factorization
  T1 mina;
  casadi_int imina;
  // Primal and dual error, corresponding index
  T1 pr, du, epr, edu;
  casadi_int ipr, idu;
  // Pending active-set change
  casadi_int index, sign;
  // Feasibility restoration active-set change
  casadi_int r_index, r_sign;
  // Iteration
  casadi_int iter;
  // Diagonal entries
  T1* D;
  // Scaling factor
  T1* S;
  // lam_lbx, lam_ubz
  T1* lam_lbz;
  T1* lam_ubz;
  // dlam_lbx, dlam_ubz
  T1* dlam_lbz;
  T1* dlam_ubz;
  // Residual
  T1* rz;
  T1* rlam;
  T1* rlam_lbz;
  T1* rlam_ubz;
  // Inverse of margin to bounds (0 if no bound)
  T1* dinv_lbz;
  T1* dinv_ubz;
  // Complementarity measure
  T1 mu;
  // Complementarity constraint error and corresponding index
  T1 co;
  casadi_int ico;
  // Number of finite constraints
  casadi_int n_con;
  // User task
  casadi_ipqp_task_t task;
  // Next step
  casadi_ipqp_next_t next;
  // Linear system
  T1* linsys;
};
// C-REPLACE "casadi_ipqp_data<T1>" "struct casadi_ipqp_data"

// SYMBOL "ipqp_init"
template<typename T1>
void casadi_ipqp_init(casadi_ipqp_data<T1>* d, casadi_int** iw, T1** w) {
  // Local variables
  casadi_int nnz_a, nnz_kkt, nnz_v, nnz_r;
  const casadi_ipqp_prob<T1>* p = d->prob;
  // Get matrix number of nonzeros
  nnz_a = p->sp_a[2+p->sp_a[1]];
  nnz_kkt = p->sp_kkt[2+p->sp_kkt[1]];
  nnz_v = p->sp_v[2+p->sp_v[1]];
  nnz_r = p->sp_r[2+p->sp_r[1]];
  d->nz_kkt = *w; *w += nnz_kkt;
  d->z = *w; *w += p->nz;
  d->lbz = *w; *w += p->nz;
  d->ubz = *w; *w += p->nz;
  d->lam = *w; *w += p->nz;
  d->dz = *w; *w += p->nz;
  d->dlam = *w; *w += p->nz;
  d->nz_v = *w; *w += casadi_max(nnz_v+nnz_r, nnz_kkt);
  d->nz_r = d->nz_v + nnz_v;
  d->beta = *w; *w += p->nz;
  d->nz_at = *w; *w += nnz_a;
  d->infeas = *w; *w += p->nx;
  d->tinfeas = *w; *w += p->nx;
  d->sens = *w; *w += p->nz;
  d->neverzero = *iw; *iw += p->nz;
  d->neverupper = *iw; *iw += p->nz;
  d->neverlower = *iw; *iw += p->nz;
  d->lincomb = *iw; *iw += p->nz;
  d->D = *w; *w += p->nz;
  d->S = *w; *w += p->nz;
  d->lam_lbz = *w; *w += p->nz;
  d->lam_ubz = *w; *w += p->nz;
  d->dlam_lbz = *w; *w += p->nz;
  d->dlam_ubz = *w; *w += p->nz;
  d->rz = *w; *w += p->nz;
  d->rlam = *w; *w += p->nz;
  d->rlam_lbz = *w; *w += p->nz;
  d->rlam_ubz = *w; *w += p->nz;
  d->dinv_lbz = *w; *w += p->nz;
  d->dinv_ubz = *w; *w += p->nz;
  d->w = *w;
  d->iw = *iw;
}

// The following routines require stdio
#ifndef CASADI_PRINTF

#endif  // CASADI_PRINTF

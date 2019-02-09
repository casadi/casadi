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
// C-REPLACE "casadi_qp_prob<T1>" "struct casadi_qp_prob"

// SYMBOL "qp_setup"
template<typename T1>
void casadi_qp_setup(casadi_qp_prob<T1>* p) {
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

// SYMBOL "qp_work"
template<typename T1>
void casadi_qp_work(const casadi_qp_prob<T1>* p, casadi_int* sz_iw, casadi_int* sz_w) {
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
}

// SYMBOL "qp_flag_t"
typedef enum {
  QP_SUCCESS,
  QP_MAX_ITER,
  QP_NO_SEARCH_DIR,
  QP_PRINTING_ERROR
} casadi_qp_flag_t;

// SYMBOL "qp_data"
template<typename T1>
struct casadi_qp_data {
  // Problem structure
  const casadi_qp_prob<T1>* prob;
  // Solver status
  casadi_qp_flag_t status;
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
};
// C-REPLACE "casadi_qp_data<T1>" "struct casadi_qp_data"

// SYMBOL "qp_init"
template<typename T1>
void casadi_qp_init(casadi_qp_data<T1>* d, casadi_int** iw, T1** w) {
  // Local variables
  casadi_int nnz_a, nnz_kkt, nnz_v, nnz_r;
  const casadi_qp_prob<T1>* p = d->prob;
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
  d->w = *w;
  d->iw = *iw;
}

// SYMBOL "qp_reset"
template<typename T1>
int casadi_qp_reset(casadi_qp_data<T1>* d) {
  // Local variables
  casadi_int i;
  const casadi_qp_prob<T1>* p = d->prob;
  // Reset variables corresponding to previous iteration
  d->msg = 0;
  d->tau = 0.;
  d->sing = 0;
  // Correct lam if needed, determine permitted signs
  for (i=0; i<p->nz; ++i) {
    // Permitted signs for lam
    d->neverzero[i] = d->lbz[i] == d->ubz[i];
    d->neverupper[i] = d->ubz[i] == p->inf;
    d->neverlower[i] = d->lbz[i] == -p->inf;
    if (d->neverzero[i] && d->neverupper[i] && d->neverlower[i]) return 1;
    // Small enough lambdas are treated as inactive
    if (!d->neverzero[i] && fabs(d->lam[i]) < p->min_lam) d->lam[i] = 0.;
    // Prevent illegal active sets
    if (d->neverzero[i] && d->lam[i] == 0.) {
      d->lam[i] = d->neverupper[i]
                || d->z[i]-d->lbz[i] <= d->ubz[i]-d->z[i] ? -p->dmin : p->dmin;
    } else if (d->neverupper[i] && d->lam[i]>0.) {
      d->lam[i] = d->neverzero[i] ? -p->dmin : 0.;
    } else if (d->neverlower[i] && d->lam[i]<0.) {
      d->lam[i] = d->neverzero[i] ? p->dmin : 0.;
    }
  }
  // Transpose A
  casadi_trans(d->nz_a, p->sp_a, d->nz_at, p->sp_at, d->iw);
  // No pending active-set change
  d->index = -2;
  d->sign = 0;
  // No restoration index
  d->r_index = -2;
  d->r_sign = 0;
  // Reset iteration counter
  d->iter = 0;
  return 0;
}

// SYMBOL "qp_pr"
template<typename T1>
void casadi_qp_pr(casadi_qp_data<T1>* d) {
  // Calculate largest constraint violation
  casadi_int i;
  const casadi_qp_prob<T1>* p = d->prob;
  d->pr = 0;
  d->ipr = -1;
  for (i=0; i<p->nz; ++i) {
    if (d->z[i] > d->ubz[i]+d->pr) {
      d->pr = d->z[i]-d->ubz[i];
      d->ipr = i;
    } else if (d->z[i] < d->lbz[i]-d->pr) {
      d->pr = d->lbz[i]-d->z[i];
      d->ipr = i;
    }
  }
}

// SYMBOL "qp_du"
template<typename T1>
void casadi_qp_du(casadi_qp_data<T1>* d) {
  // Calculate largest constraint violation
  casadi_int i;
  const casadi_qp_prob<T1>* p = d->prob;
  d->du = 0;
  d->idu = -1;
  for (i=0; i<p->nx; ++i) {
    if (d->infeas[i] > d->du) {
      d->du = d->infeas[i];
      d->idu = i;
    } else if (d->infeas[i] < -d->du) {
      d->du = -d->infeas[i];
      d->idu = i;
    }
  }
}

// SYMBOL "qp_du_check"
template<typename T1>
int casadi_qp_du_check(casadi_qp_data<T1>* d, casadi_int i) {
  // Local variables
  casadi_int k;
  T1 new_du;
  const casadi_int *at_colind, *at_row;
  const casadi_qp_prob<T1>* p = d->prob;
  // AT sparsity
  at_colind = p->sp_at + 2;
  at_row = at_colind + p->na + 1;
  // Maximum infeasibility from setting from setting lam[i]=0
  if (i<p->nx) {
    new_du = fabs(d->infeas[i]-d->lam[i]);
  } else {
    new_du = 0.;
    for (k=at_colind[i-p->nx]; k<at_colind[i-p->nx+1]; ++k) {
      new_du = fmax(new_du, fabs(d->infeas[at_row[k]]-d->nz_at[k]*d->lam[i]));
    }
  }
  return new_du <= d->du;
}

// SYMBOL "qp_du_index"
template<typename T1>
void casadi_qp_du_index(casadi_qp_data<T1>* d) {
  // Try to improve dual feasibility by removing a constraint
  // Local variables
  casadi_int i, s;
  T1 best_sens;
  const casadi_qp_prob<T1>* p = d->prob;
  // Find the best lam[i] to make zero
  d->index = -1;
  best_sens = -1;
  for (i = 0; i < p->nz; ++i) {
    // Skip if no dual infeasibility sensitivity
    if (d->sens[i] == 0.) continue;
    // Is the constraint enforced?
    if (d->lam[i] == 0) {
      // We're enforcing constraints
      s = d->sens[i] > 0 ? 1 : -1;
      // Make sure that enforcing the constraint is possible
      if (s > 0 ? d->neverupper[i] : d->neverlower[i]) continue;
    } else {
      // We're removing constraints
      s = 0;
      // Make sure that it's a constraint that can be removed
      if (d->neverzero[i]) continue;
      // If variable influences du, make sure sign is right
      if (d->lam[i] > 0. ? d->sens[i] > 0. : d->sens[i] < 0.) continue;
      // Skip if maximum infeasibility increases
      if (!casadi_qp_du_check(d, i)) continue;
    }
    // Check if best so far
    if (fabs(d->sens[i]) > best_sens) {
      best_sens = fabs(d->sens[i]);
      d->index = i;
      d->sign = s;
    }
  }
  // Accept, if any
  if (d->index >= 0) {
    if (d->sign > 0) {
      d->msg = "Enforced ubz to reduce |du|";
    } else if (d->sign < 0) {
      d->msg = "Enforced lbz to reduce |du|";
    } else if (d->lam[d->index] > 0) {
      d->msg = "Dropped ubz to reduce |du|";
    } else {
      d->msg = "Dropped lbz to reduce |du|";
    }
    d->msg_ind = d->index;
  }
}

// SYMBOL "qp_pr_index"
template<typename T1>
void casadi_qp_pr_index(casadi_qp_data<T1>* d) {
  // Try to improve primal feasibility by adding a constraint
  if (d->lam[d->ipr] == 0.) {
    // Add the most violating constraint
    if (d->z[d->ipr] < d->lbz[d->ipr]) {
      d->sign = -1;
      d->msg = "Added lbz to reduce |pr|";
    } else {
      d->sign = 1;
      d->msg = "Added ubz to reduce |pr|";
    }
    d->msg_ind = d->ipr;
    d->index = d->ipr;
  } else {
    // No improvement possible
    d->index = -1;
  }
}

// SYMBOL "qp_kkt"
template<typename T1>
void casadi_qp_kkt(casadi_qp_data<T1>* d) {
  // Local variables
  casadi_int i, k;
  const casadi_int *h_colind, *h_row, *a_colind, *a_row, *at_colind, *at_row,
                   *kkt_colind, *kkt_row;
  const casadi_qp_prob<T1>* p = d->prob;
  // Extract sparsities
  a_row = (a_colind = p->sp_a+2) + p->nx + 1;
  at_row = (at_colind = p->sp_at+2) + p->na + 1;
  h_row = (h_colind = p->sp_h+2) + p->nx + 1;
  kkt_row = (kkt_colind = p->sp_kkt+2) + p->nz + 1;
  // Reset w to zero
  casadi_clear(d->w, p->nz);
  // Loop over rows of the (transposed) KKT
  for (i=0; i<p->nz; ++i) {
    // Copy row of KKT to w
    if (i<p->nx) {
      if (d->lam[i]==0) {
        for (k=h_colind[i]; k<h_colind[i+1]; ++k) d->w[h_row[k]] = d->nz_h[k];
        for (k=a_colind[i]; k<a_colind[i+1]; ++k) d->w[p->nx+a_row[k]] = d->nz_a[k];
      } else {
        d->w[i] = 1.;
      }
    } else {
      if (d->lam[i]==0) {
        d->w[i] = -1.;
      } else {
        for (k=at_colind[i-p->nx]; k<at_colind[i-p->nx+1]; ++k) {
          d->w[at_row[k]] = d->nz_at[k];
        }
      }
    }
    // Copy row to KKT, zero out w
    for (k=kkt_colind[i]; k<kkt_colind[i+1]; ++k) {
      d->nz_kkt[k] = d->w[kkt_row[k]];
      d->w[kkt_row[k]] = 0;
    }
  }
}

// SYMBOL "qp_kkt_vector"
template<typename T1>
void casadi_qp_kkt_vector(casadi_qp_data<T1>* d, T1* kkt_i, casadi_int i) {
  // Local variables
  casadi_int k;
  const casadi_int *h_colind, *h_row, *a_colind, *a_row, *at_colind, *at_row;
  const casadi_qp_prob<T1>* p = d->prob;
  // Extract sparsities
  a_row = (a_colind = p->sp_a+2) + p->nx + 1;
  at_row = (at_colind = p->sp_at+2) + p->na + 1;
  h_row = (h_colind = p->sp_h+2) + p->nx + 1;
  // Reset kkt_i to zero
  casadi_clear(kkt_i, p->nz);
  // Copy sparse entries
  if (i<p->nx) {
    for (k=h_colind[i]; k<h_colind[i+1]; ++k) kkt_i[h_row[k]] = d->nz_h[k];
    for (k=a_colind[i]; k<a_colind[i+1]; ++k) kkt_i[p->nx+a_row[k]] = d->nz_a[k];
  } else {
    for (k=at_colind[i-p->nx]; k<at_colind[i-p->nx+1]; ++k) {
      kkt_i[at_row[k]] = -d->nz_at[k];
    }
  }
  // Add diagonal entry
  kkt_i[i] -= 1.;
}

// SYMBOL "qp_kkt_dot"
template<typename T1>
T1 casadi_qp_kkt_dot(casadi_qp_data<T1>* d, const T1* v, casadi_int i) {
  // Local variables
  casadi_int k;
  const casadi_int *h_colind, *h_row, *a_colind, *a_row, *at_colind, *at_row;
  T1 r;
  const casadi_qp_prob<T1>* p = d->prob;
  // Extract sparsities
  a_row = (a_colind = p->sp_a + 2) + p->nx + 1;
  at_row = (at_colind = p->sp_at + 2) + p->na + 1;
  h_row = (h_colind = p->sp_h + 2) + p->nx + 1;
  // Scalar product with the diagonal
  r = v[i];
  // Scalar product with the sparse entries
  if (i < p->nx) {
    for (k=h_colind[i]; k<h_colind[i+1]; ++k) r -= v[h_row[k]] * d->nz_h[k];
    for (k=a_colind[i]; k<a_colind[i+1]; ++k) r -= v[p->nx+a_row[k]] * d->nz_a[k];
  } else {
    for (k=at_colind[i-p->nx]; k<at_colind[i-p->nx+1]; ++k) {
      r += v[at_row[k]] * d->nz_at[k];
    }
  }
  return r;
}

// SYMBOL "qp_kkt_residual"
template<typename T1>
void casadi_qp_kkt_residual(casadi_qp_data<T1>* d, T1* r) {
  casadi_int i;
  const casadi_qp_prob<T1>* p = d->prob;
  for (i=0; i<p->nz; ++i) {
    if (d->lam[i]>0.) {
      r[i] = d->ubz[i]-d->z[i];
    } else if (d->lam[i]<0.) {
      r[i] = d->lbz[i]-d->z[i];
    } else if (i<p->nx) {
      r[i] = d->lam[i]-d->infeas[i];
    } else {
      r[i] = d->lam[i];
    }
  }
}

// SYMBOL "qp_zero_blocking"
template<typename T1>
int casadi_qp_zero_blocking(casadi_qp_data<T1>* d) {
  // Local variables
  casadi_int i;
  T1 dz_max = 0;
  const casadi_qp_prob<T1>* p = d->prob;
  // Look for violated constraints that are not improving
  for (i = 0; i < p->nz; ++i) {
    if (d->dz[i] < -dz_max && d->lbz[i] - d->z[i] >= d->epr) {
      dz_max = -d->dz[i];
      d->index = i;
      d->sign = -1;
      d->msg = "lbz violated with zero step";
      d->msg_ind = d->index;
    } else if (d->dz[i] > dz_max && d->z[i] - d->ubz[i] >= d->epr) {
      dz_max = d->dz[i];
      d->index = i;
      d->sign = 1;
      d->msg = "ubz violated with zero step";
      d->msg_ind = d->index;
    }
  }
  return dz_max > 0;
}

// SYMBOL "qp_primal_blocking"
template<typename T1>
void casadi_qp_primal_blocking(casadi_qp_data<T1>* d) {
  // Local variables
  casadi_int i;
  T1 trial_z;
  const casadi_qp_prob<T1>* p = d->prob;
  // Check if violation with tau=0 and not improving
  if (casadi_qp_zero_blocking(d)) {
    d->tau = 0.;
    return;
  }
  // Loop over all primal variables
  for (i = 0; i < p->nz; ++i) {
    if (d->dz[i] == 0.) continue; // Skip zero steps
    // Trial primal step
    trial_z = d->z[i] + d->tau * d->dz[i];
    if (d->dz[i] < 0 && trial_z < d->lbz[i] - d->epr) {
      // Trial would increase maximum infeasibility
      d->tau = (d->lbz[i] - d->epr - d->z[i]) / d->dz[i];
      d->index = d->lam[i] < 0. ? -1 : i;
      d->sign = -1;
      d->msg = "Enforcing lbz";
      d->msg_ind = i;
    } else if (d->dz[i] > 0 && trial_z > d->ubz[i] + d->epr) {
      // Trial would increase maximum infeasibility
      d->tau = (d->ubz[i] + d->epr - d->z[i]) / d->dz[i];
      d->index = d->lam[i] > 0. ? -1 : i;
      d->sign = 1;
      d->msg = "Enforcing ubz";
      d->msg_ind = i;
    }
    if (d->tau <= 0) return;
  }
}

// SYMBOL "qp_dual_breakpoints"
template<typename T1>
casadi_int casadi_qp_dual_breakpoints(casadi_qp_data<T1>* d, T1* tau_list,
                                      casadi_int* ind_list, T1 tau) {
  // Local variables
  casadi_int i, n_tau, loc, next_ind, tmp_ind, j;
  T1 trial_lam, new_tau, next_tau, tmp_tau;
  const casadi_qp_prob<T1>* p = d->prob;
  // Dual feasibility is piecewise linear. Start with one interval [0,tau]:
  tau_list[0] = tau;
  ind_list[0] = -1; // no associated index
  n_tau = 1;
  // Find the taus corresponding to lam crossing zero and insert into list
  for (i=0; i<p->nz; ++i) {
    if (d->dlam[i]==0.) continue; // Skip zero steps
    if (d->lam[i]==0.) continue; // Skip inactive constraints
    // Trial dual step
    trial_lam = d->lam[i] + tau*d->dlam[i];
    // Skip if no sign change
    if (d->lam[i]>0 ? trial_lam>=0 : trial_lam<=0) continue;
    // Location of the sign change
    new_tau = -d->lam[i]/d->dlam[i];
    // Where to insert the w[i]
    for (loc=0; loc<n_tau-1; ++loc) {
      if (new_tau<tau_list[loc]) break;
    }
    // Insert element
    n_tau++;
    next_tau=new_tau;
    next_ind=i;
    for (j=loc; j<n_tau; ++j) {
      tmp_tau = tau_list[j];
      tau_list[j] = next_tau;
      next_tau = tmp_tau;
      tmp_ind = ind_list[j];
      ind_list[j] = next_ind;
      next_ind = tmp_ind;
    }
  }
  return n_tau;
}

// SYMBOL "qp_dual_blocking"
template<typename T1>
casadi_int casadi_qp_dual_blocking(casadi_qp_data<T1>* d) {
  // Local variables
  casadi_int i, n_tau, j, k, du_index;
  T1 tau_k, dtau, new_infeas, tau1, infeas, tinfeas;
  const casadi_int *at_colind, *at_row;
  const casadi_qp_prob<T1>* p = d->prob;
  // Extract sparsities
  at_row = (at_colind = p->sp_at+2) + p->na + 1;
  // Dual feasibility is piecewise linear in tau. Get the intervals:
  n_tau = casadi_qp_dual_breakpoints(d, d->w, d->iw, d->tau);
  // No dual blocking yet
  du_index = -1;
  // How long step can we take without exceeding e?
  tau_k = 0.;
  for (j=0; j<n_tau; ++j) {
    // Distance to the next tau (may be zero)
    dtau = d->w[j] - tau_k;
    // Check if maximum dual infeasibilty gets exceeded
    for (k=0; k<p->nx; ++k) {
      // Get infeasibility and infeasibility tangent
      infeas  = d->infeas[k];
      tinfeas  = d->tinfeas[k];
      // Make sure tinfeas>0
      if (fabs(tinfeas)<1e-14) {
        // Skip
        continue;
      } else if (tinfeas<0) {
        // Switch signs
        infeas *= -1;
        tinfeas *= -1;
      }
      // Tentative new infeasibility
      new_infeas = infeas + dtau*tinfeas;
      // Does infeasibility get exceeded
      if (new_infeas > d->edu) {
        // Sign change and exceeded
        tau1 = fmax(tau_k, tau_k + (d->edu - infeas)/tinfeas);
        if (tau1 < d->tau) {
          // Enforce dual blocking constraint
          d->tau = tau1;
          du_index = k;
        }
      }
    }
    // Update infeasibility
    casadi_axpy(p->nx, fmin(d->tau - tau_k, dtau), d->tinfeas, d->infeas);
    // Stop here if dual blocking constraint
    if (du_index>=0) return du_index;
    // Continue to the next tau
    tau_k = d->w[j];
    // Get component, break if last
    i = d->iw[j];
    if (i<0) break;
    // Update sign or tinfeas
    if (!d->neverzero[i]) {
      // lam becomes zero, update the infeasibility tangent
      if (i<p->nx) {
        // Set a lam_x to zero
        d->tinfeas[i] -= d->dlam[i];
      } else {
        // Set a lam_a to zero
        for (k=at_colind[i-p->nx]; k<at_colind[i-p->nx+1]; ++k) {
          d->tinfeas[at_row[k]] -= d->nz_at[k]*d->dlam[i];
        }
      }
    }
  }
  return du_index;
}

// SYMBOL "qp_take_step"
template<typename T1>
void casadi_qp_take_step(casadi_qp_data<T1>* d) {
  // Local variables
  casadi_int i;
  const casadi_qp_prob<T1>* p = d->prob;
  // Get current sign
  for (i=0; i<p->nz; ++i) d->iw[i] = d->lam[i]>0. ? 1 : d->lam[i]<0 ? -1 : 0;
  // Take primal-dual step
  casadi_axpy(p->nz, d->tau, d->dz, d->z);
  casadi_axpy(p->nz, d->tau, d->dlam, d->lam);
  // Update sign
  for (i=0; i<p->nz; ++i) {
    // Allow sign changes for certain components
    if (d->neverzero[i] && (d->iw[i]<0 ? d->lam[i]>0 : d->lam[i]<0)) {
      d->iw[i]=-d->iw[i];
    }
    // Ensure correct sign
    switch (d->iw[i]) {
      case -1: d->lam[i] = fmin(d->lam[i], -p->dmin); break;
      case  1: d->lam[i] = fmax(d->lam[i],  p->dmin); break;
      case  0: d->lam[i] = 0.; break;
    }
  }
}

// SYMBOL "qp_flip_check"
template<typename T1>
int casadi_qp_flip_check(casadi_qp_data<T1>* d) {
  const casadi_qp_prob<T1>* p = d->prob;
  // Calculate the difference between unenforced and enforced column index
  casadi_qp_kkt_vector(d, d->dlam, d->index);
  // Calculate the difference between old and new column index
  if (d->sign == 0) casadi_scal(p->nz, -1., d->dlam);
  // Try to find a linear combination of the new columns
  casadi_qr_solve(d->dlam, 1, 0, p->sp_v, d->nz_v, p->sp_r, d->nz_r, d->beta,
    p->prinv, p->pc, d->w);
  // If dlam[index]!=1, new columns must be linearly independent
  if (fabs(d->dlam[d->index]-1.) >= 1e-12) return 0;
  // Next, find a linear combination of the new rows
  casadi_clear(d->dz, p->nz);
  d->dz[d->index] = 1;
  casadi_qr_solve(d->dz, 1, 1, p->sp_v, d->nz_v, p->sp_r, d->nz_r, d->beta,
    p->prinv, p->pc, d->w);
  // Normalize dlam, dz
  casadi_scal(p->nz, 1./sqrt(casadi_dot(p->nz, d->dlam, d->dlam)), d->dlam);
  casadi_scal(p->nz, 1./sqrt(casadi_dot(p->nz, d->dz, d->dz)), d->dz);
  // KKT system will be singular
  return 1;
}

// SYMBOL "qp_factorize"
template<typename T1>
void casadi_qp_factorize(casadi_qp_data<T1>* d) {
  const casadi_qp_prob<T1>* p = d->prob;
  // Do we already have a search direction due to lost singularity?
  if (d->has_search_dir) {
    d->sing = 1;
    return;
  }
  // Construct the KKT matrix
  casadi_qp_kkt(d);
  // QR factorization
  casadi_qr(p->sp_kkt, d->nz_kkt, d->w, p->sp_v, d->nz_v, p->sp_r,
            d->nz_r, d->beta, p->prinv, p->pc);
  // Check singularity
  d->sing = casadi_qr_singular(&d->mina, &d->imina, d->nz_r, p->sp_r, p->pc, 1e-12);
}

// SYMBOL "qp_expand_step"
template<typename T1>
void casadi_qp_expand_step(casadi_qp_data<T1>* d) {
  // Local variables
  casadi_int i;
  const casadi_qp_prob<T1>* p = d->prob;
  // Calculate change in Lagrangian gradient
  casadi_clear(d->dlam, p->nx);
  casadi_mv(d->nz_h, p->sp_h, d->dz, d->dlam, 0); // gradient of the objective
  casadi_mv(d->nz_a, p->sp_a, d->dz + p->nx, d->dlam, 1); // gradient of the Lagrangian
  // Step in lam[:nx]
  casadi_scal(p->nx, -1., d->dlam);
  // For inactive constraints, lam(x) step is zero
  for (i = 0; i < p->nx; ++i) if (d->lam[i] == 0.) d->dlam[i] = 0.;
  // Step in lam[nx:]
  casadi_copy(d->dz+p->nx, p->na, d->dlam + p->nx);
  // Step in z[nx:]
  casadi_clear(d->dz + p->nx, p->na);
  casadi_mv(d->nz_a, p->sp_a, d->dz, d->dz + p->nx, 0);
  // Avoid steps that are nonzero due to numerics
  for (i = 0; i < p->nz; ++i) if (fabs(d->dz[i]) < 1e-14) d->dz[i] = 0.;
  // Tangent of the dual infeasibility at tau=0
  casadi_clear(d->tinfeas, p->nx);
  casadi_mv(d->nz_h, p->sp_h, d->dz, d->tinfeas, 0);
  casadi_mv(d->nz_a, p->sp_a, d->dlam + p->nx, d->tinfeas, 1);
  casadi_axpy(p->nx, 1., d->dlam, d->tinfeas);
}

// SYMBOL "casadi_qp_pr_direction"
template<typename T1>
int casadi_qp_pr_direction(casadi_qp_data<T1>* d) {
  casadi_int i;
  const casadi_qp_prob<T1>* p = d->prob;
  for (i=0; i<p->nz; ++i) {
    if (d->lbz[i] - d->z[i] >= d->epr) {
      // Prevent further violation of lower bound
      if (d->dz[i] < 0 || d->dlam[i] > 0) return 1;
    } else if (d->z[i] - d->ubz[i] >= d->epr) {
      // Prevent further violation of upper bound
      if (d->dz[i] > 0 || d->dlam[i] < 0) return 1;
    }
  }
  return 0;
}

// SYMBOL "casadi_qp_du_direction"
template<typename T1>
int casadi_qp_du_direction(casadi_qp_data<T1>* d) {
  casadi_int i;
  const casadi_qp_prob<T1>* p = d->prob;
  for (i=0; i<p->nx; ++i) {
    // Prevent further increase in dual infeasibility
    if (d->infeas[i] <= -d->edu && d->tinfeas[i] < -1e-12) {
      return 1;
    } else if (d->infeas[i] >= d->edu && d->tinfeas[i] > 1e-12) {
      return 1;
    }
  }
  return 0;
}

// SYMBOL "qp_enforceable"
template<typename T1>
int casadi_qp_enforceable(casadi_qp_data<T1>* d, casadi_int i, casadi_int s) {
  // Local variables
  casadi_int k;
  const casadi_int *at_colind, *at_row;
  const casadi_qp_prob<T1>* p = d->prob;
  // Can always enforce if not at bound
  if (fabs(d->infeas[i]) < d->edu) return 1;
  // AT sparsity
  at_colind = p->sp_at + 2;
  at_row = at_colind + p->na + 1;
  // Can we set lam[i] := s*DMIN without exceeding edu?
  if (i<p->nx) {
    return (s < 0) == (d->infeas[i] > 0);
  } else {
    for (k=at_colind[i-p->nx]; k<at_colind[i-p->nx+1]; ++k) {
      if (d->nz_at[k] > 0) {
        if ((s > 0) == (d->infeas[at_row[k]] > 0)) return 0;
      } else if (d->nz_at[k] < 0) {
        if ((s < 0) == (d->infeas[at_row[k]] > 0)) return 0;
      }
    }
    return 1;
  }
}

// SYMBOL "qp_singular_step"
// C-REPLACE "static_cast<T1*>(0)" "0"
template<typename T1>
int casadi_qp_singular_step(casadi_qp_data<T1>* d) {
  // Local variables
  T1 tau_test, tau;
  casadi_int nnz_kkt, nk, k, i, best_k, best_neg, neg;
  const casadi_qp_prob<T1>* p = d->prob;
  // Find the columns that take part in any linear combination
  for (i = 0; i < p->nz; ++i) d->lincomb[i] = 0;
  for (k = 0; k < d->sing; ++k) {
    if (!d->has_search_dir) {
      casadi_qr_colcomb(d->dlam, d->nz_r, p->sp_r, p->pc, 1e-12, k);
    }
    for (i = 0; i < p->nz; ++i) if (fabs(d->dlam[i]) >= 1e-12) d->lincomb[i]++;
  }

  if (d->has_search_dir) {
    // One, given search direction
    nk = 1;
  } else {
    // QR factorization of the transpose
    casadi_trans(d->nz_kkt, p->sp_kkt, d->nz_v, p->sp_kkt, d->iw);
    nnz_kkt = p->sp_kkt[2+p->nz]; // kkt_colind[nz]
    casadi_copy(d->nz_v, nnz_kkt, d->nz_kkt);
    casadi_qr(p->sp_kkt, d->nz_kkt, d->w, p->sp_v, d->nz_v, p->sp_r, d->nz_r,
              d->beta, p->prinv, p->pc);
    // For all nullspace vectors
    nk = casadi_qr_singular(static_cast<T1*>(0), 0, d->nz_r, p->sp_r, p->pc, 1e-12);
  }
  // Best flip
  best_k = best_neg = -1;
  tau = p->inf;
  for (k=0; k<nk; ++k) {
    if (!d->has_search_dir) {
      // Get a linear combination of the rows in kkt
      casadi_qr_colcomb(d->dz, d->nz_r, p->sp_r, p->pc, 1e-12, k);
    }
    // Which constraints can be flipped in order to increase rank?
    for (i=0; i<p->nz; ++i) {
      d->iw[i] = d->lincomb[i] && fabs(casadi_qp_kkt_dot(d, d->dz, i)) > 1e-12;
    }
    // Calculate step, dz and dlam
    casadi_qp_expand_step(d);
    // Try both positive and negative direction
    for (neg = 0; neg < 2; ++neg) {
      // Negate direction
      if (neg) {
        casadi_scal(p->nz, -1., d->dz);
        casadi_scal(p->nz, -1., d->dlam);
        casadi_scal(p->nx, -1., d->tinfeas);
      }
      // Make sure primal infeasibility doesn't exceed limits
      if (casadi_qp_pr_direction(d)) continue;
      // Make sure dual infeasibility doesn't exceed limits
      if (casadi_qp_du_direction(d)) continue;
      // Loop over potential active set changes
      for (i=0; i<p->nz; ++i) {
        // Skip if no rank increase
        if (!d->iw[i]) continue;
        // Enforced or not?
        if (d->lam[i]==0.) {
          if (d->z[i] <= d->ubz[i] && (d->z[i] >= d->lbz[i] ?
              d->dz[i] < -1e-12 : d->dz[i] > 1e-12)) {
            // Enforce lower bound?
            if (!d->neverlower[i]
                && (tau_test = (d->lbz[i] - d->z[i]) / d->dz[i]) < tau
                && casadi_qp_enforceable(d, i, -1)) {
              tau = tau_test;
              d->r_index = i;
              d->r_sign = -1;
              best_k = k;
              best_neg = neg;
            }
          } else if (d->z[i] >= d->lbz[i] && (d->z[i] <= d->ubz[i] ?
              d->dz[i] > 1e-12 : d->dz[i] < -1e-12)) {
            // Enforce upper bound?
            if (!d->neverupper[i]
                && (tau_test = (d->ubz[i] - d->z[i]) / d->dz[i]) < tau
                && casadi_qp_enforceable(d, i, 1)) {
              tau = tau_test;
              d->r_index = i;
              d->r_sign = 1;
              best_k = k;
              best_neg = neg;
            }
          }
        } else if (!d->neverzero[i]) {
          // Drop a constraint?
          if (d->lam[i] > 0 ? d->dlam[i] < -1e-12 : d->dlam[i] > 1e-12) {
            if ((tau_test = -d->lam[i] / d->dlam[i]) < tau) {
              tau = tau_test;
              d->r_index = i;
              d->r_sign = 0;
              best_k = k;
              best_neg = neg;
            }
          }
        }
      }
    }
  }
  // Can we restore feasibility?
  if (d->r_index < 0) return 1;
  // Recalculate direction, if needed
  if (--k != best_k) {
    // Need to recalculate direction
    casadi_qr_colcomb(d->dz, d->nz_r, p->sp_r, p->pc, 1e-12, best_k);
    casadi_qp_expand_step(d);
    if (best_neg) tau *= -1;
  } else if (--neg != best_neg) {
    // No need to recalculate, but opposite direction
    tau *= -1;
  }
  // Scale step so that that tau=1 corresponds to a full step
  casadi_scal(p->nz, tau, d->dz);
  casadi_scal(p->nz, tau, d->dlam);
  casadi_scal(p->nx, tau, d->tinfeas);
  return 0;
}

// SYMBOL "qp_calc_step"
template<typename T1>
int casadi_qp_calc_step(casadi_qp_data<T1>* d) {
  // Local variables
  const casadi_qp_prob<T1>* p = d->prob;
  // Reset returns
  d->r_index = -1;
  d->r_sign = 0;
  // Handle singularity
  if (d->sing) return casadi_qp_singular_step(d);
  // Negative KKT residual
  casadi_qp_kkt_residual(d, d->dz);
  // Solve to get step in z[:nx] and lam[nx:]
  casadi_qr_solve(d->dz, 1, 1, p->sp_v, d->nz_v, p->sp_r, d->nz_r, d->beta,
                  p->prinv, p->pc, d->w);
  // Have step in dz[:nx] and dlam[nx:]. Calculate complete dz and dlam
  casadi_qp_expand_step(d);
  // Successful return
  return 0;
}

// SYMBOL "qp_calc_sens"
template<typename T1>
void casadi_qp_calc_sens(casadi_qp_data<T1>* d, casadi_int i) {
  // Local variables
  const casadi_qp_prob<T1>* p = d->prob;
  // Calculate sensitivities in decreasing dual infeasibility index i
  casadi_clear(d->sens, p->nz);
  if (i >= 0) {
    d->sens[i] = d->infeas[i] > 0 ? -1. : 1.;
    casadi_mv(d->nz_a, p->sp_a, d->sens, d->sens + p->nx, 0);
  }
}

// SYMBOL "qp_calc_dependent"
template<typename T1>
void casadi_qp_calc_dependent(casadi_qp_data<T1>* d) {
  // Local variables
  casadi_int i;
  T1 r;
  const casadi_qp_prob<T1>* p = d->prob;
  // Calculate f
  d->f = casadi_bilin(d->nz_h, p->sp_h, d->z, d->z)/2.
       + casadi_dot(p->nx, d->z, d->g);
  // Calculate z[nx:]
  casadi_clear(d->z+p->nx, p->na);
  casadi_mv(d->nz_a, p->sp_a, d->z, d->z+p->nx, 0);
  // Calculate gradient of the Lagrangian
  casadi_copy(d->g, p->nx, d->infeas);
  casadi_mv(d->nz_h, p->sp_h, d->z, d->infeas, 0);
  casadi_mv(d->nz_a, p->sp_a, d->lam+p->nx, d->infeas, 1);
  // Calculate lam[:nx] without changing the sign accidentally, dual infeasibility
  for (i=0; i<p->nx; ++i) {
    // No change if zero
    if (d->lam[i]==0) continue;
    // lam[i] with no sign restrictions
    r = -d->infeas[i];
    if (d->lam[i]>0) {
      if (d->neverzero[i] && !d->neverlower[i]) {
        d->lam[i] = r==0 ? p->dmin : r; // keep sign if r==0
      } else {
        d->lam[i] = fmax(r, p->dmin); // no sign change
      }
    } else {
      if (d->neverzero[i] && !d->neverupper[i]) {
        d->lam[i] = r==0 ? -p->dmin : r; // keep sign if r==0
      } else {
        d->lam[i] = fmin(r, -p->dmin); // no sign change
      }
    }
    // Update dual infeasibility
    d->infeas[i] += d->lam[i];
  }
  // Calculate primal and dual error
  casadi_qp_pr(d);
  casadi_qp_du(d);
  // Acceptable primal and dual error
  d->epr = fmax(d->pr, (0.5 * p->constr_viol_tol / p->dual_inf_tol) * d->du);
  d->edu = fmax(d->du, (0.5 * p->dual_inf_tol / p->constr_viol_tol) * d->pr);
  // Sensitivity in decreasing |du|
  casadi_qp_calc_sens(d, d->idu);
}

// SYMBOL "qp_linesearch"
template<typename T1>
void casadi_qp_linesearch(casadi_qp_data<T1>* d) {
  // Local variables
  casadi_int du_index;
  // Start with a full step and no active set change
  d->sign = 0;
  d->index = -1;
  d->tau = 1.;
  // Find largest possible step without exceeding acceptable |pr|
  casadi_qp_primal_blocking(d);
  // Find largest possible step without exceeding acceptable |du|
  du_index = casadi_qp_dual_blocking(d);
  // Take primal-dual step, avoiding accidental sign changes for lam
  casadi_qp_take_step(d);
  // Handle dual blocking constraints
  if (du_index >= 0) {
    // Sensititivity in decreasing du_index
    casadi_qp_calc_sens(d, du_index);
    // Find corresponding index
    casadi_qp_du_index(d);
  }
}

// SYMBOL "qp_flip"
template<typename T1>
void casadi_qp_flip(casadi_qp_data<T1>* d) {
  // Local variables
  const casadi_qp_prob<T1>* p = d->prob;
  // Try to restore regularity if possible
  if (d->index == -1 && d->r_index >= 0) {
    if (d->r_sign != 0 || casadi_qp_du_check(d, d->r_index)) {
      d->index = d->r_index;
      d->sign = d->r_sign;
      if (d->sign > 0) {
        d->msg = "Enforced ubz for regularity";
      } else if (d->sign < 0) {
        d->msg = "Enforced lbz for regularity";
      } else if (d->lam[d->index] > 0) {
        d->msg = "Dropped ubz for regularity";
      } else {
        d->msg = "Dropped lbz for regularity";
      }
      d->msg_ind = d->index;
    }
  }
  // If nonsingular and nonzero error, try to flip a constraint
  if (!d->sing && d->index == -1) {
    if (d->pr * p->dual_inf_tol >= p->constr_viol_tol * d->du) {
      // Improve primal feasibility if dominating
      if (d->pr >= p->constr_viol_tol) casadi_qp_pr_index(d);
    } else {
      // Improve dual feasibility if dominating
      if (d->du >= p->dual_inf_tol) casadi_qp_du_index(d);
    }
  }
  // No search direction given by default
  d->has_search_dir = 0;
  // If a constraint was added
  if (d->index >= 0) {
    // Detect singularity before it happens and get nullspace vectors
    if (!d->sing) d->has_search_dir = casadi_qp_flip_check(d);
    // Perform the active-set change
    d->lam[d->index] = d->sign==0 ? 0 : d->sign > 0 ? p->dmin : -p->dmin;
    // Recalculate primal and dual infeasibility
    casadi_qp_calc_dependent(d);
  }
}

// SYMBOL "qp_prepare"
template<typename T1>
int casadi_qp_prepare(casadi_qp_data<T1>* d) {
  // Local variables
  const casadi_qp_prob<T1>* p = d->prob;
  // Calculate dependent quantities
  casadi_qp_calc_dependent(d);
  // Make an active set change
  casadi_qp_flip(d);
  // Form and factorize the KKT system
  casadi_qp_factorize(d);
  // Termination message
  if (!d->sing && d->index == -1) {
    d->status = QP_SUCCESS;
    d->msg = "Converged";
    d->msg_ind = -2;
    return 1;
  } else if (d->iter >= p->max_iter) {
    d->status = QP_MAX_ITER;
    d->msg = "Max iter";
    d->msg_ind = -2;
    return 1;
  } else if (!d->sing && d->ipr < 0 && d->idu < 0) {
    d->status = QP_SUCCESS;
    d->msg = "No primal or dual error";
    d->msg_ind = -2;
    return 1;
  } else {
    // Keep iterating
    return 0;
  }
}

// SYMBOL "qp_iterate"
template<typename T1>
int casadi_qp_iterate(casadi_qp_data<T1>* d) {
  // Reset message flag
  d->msg = 0;
  // Start a new iteration
  d->iter++;
  // Calculate search direction
  if (casadi_qp_calc_step(d)) {
    d->status = QP_NO_SEARCH_DIR;
    return 1;
  }
  // Line search in the calculated direction
  casadi_qp_linesearch(d);
  // Keep iterating
  return 0;
}

// The following routines require stdio
#ifndef CASADI_PRINTF

// SYMBOL "qp_print_header"
template<typename T1>
int casadi_qp_print_header(casadi_qp_data<T1>* d, char* buf, size_t buf_sz) {
  int flag;
  // Print to string
  flag = snprintf(buf, buf_sz, "%5s %5s %9s %9s %5s %9s %5s %9s %5s %9s  %4s",
          "Iter", "Sing", "fk", "|pr|", "con", "|du|", "var",
          "min_R", "con", "last_tau", "Note");
  // Check if error
  if (flag < 0) {
    d->status = QP_PRINTING_ERROR;
    return 1;
  }
  // Successful return
  return 0;
}

// SYMBOL "qp_print_colcomb"
template<typename T1>
int casadi_qp_print_colcomb(casadi_qp_data<T1>* d, char* buf, size_t buf_sz, casadi_int j) {
  casadi_int num_size, n_print, i, k, buf_offset, val;
  size_t b;
  const casadi_qp_prob<T1>* p = d->prob;
  casadi_qr_colcomb(d->dlam, d->nz_r, p->sp_r, p->pc, 1e-12, j);

  // Determine max printing size
  num_size = 1;
  val = p->nz-1;
  while (val) {
    val/=10;
    num_size++;
  }

  if (buf_sz<=4) return 1;

  // How many numbers can be printed?
  // Need some extra space for '...'
  // and null
  n_print = (buf_sz-4)/num_size;

  // Clear buffer
  for (b=0;b<buf_sz;++b) buf[b]=' ';

  buf_offset = 0;
  for (i=0;i<p->nz;++i) {
    if (fabs(d->dlam[i]) >= 1e-12) {
      if (n_print==0) {
        buf[buf_sz-4] = '.';
        buf[buf_sz-3] = '.';
        buf[buf_sz-2] = '.';
        buf[buf_sz-1] = '\0';
        return 1;
      }
      n_print--;
      snprintf(buf+buf_offset, num_size, "%d", static_cast<int>(i));
      // Clear null chars
      for (k=0;k<num_size;++k) {
        if (buf[buf_offset+k]=='\0') buf[buf_offset+k] = ' ';
      }
      buf_offset += num_size;
    }
  }
  buf[buf_sz-1] = '\0';

  // Successful return
  return 0;
}

// SYMBOL "qp_print_iteration"
template<typename T1>
int casadi_qp_print_iteration(casadi_qp_data<T1>* d, char* buf, int buf_sz) {
  int flag;
  // Print iteration data without note to string
  flag = snprintf(buf, buf_sz,
    "%5d %5d %9.2g %9.2g %5d %9.2g %5d %9.2g %5d %9.2g  ",
    static_cast<int>(d->iter), static_cast<int>(d->sing), d->f, d->pr,
    static_cast<int>(d->ipr), d->du, static_cast<int>(d->idu),
    d->mina, static_cast<int>(d->imina), d->tau);
  // Check if error
  if (flag < 0) {
    d->status = QP_PRINTING_ERROR;
    return 1;
  }
  // Rest of buffer reserved for iteration note
  buf += flag;
  buf_sz -= flag;
  // Print iteration note, if any
  if (d->msg) {
    if (d->msg_ind > -2) {
      flag = snprintf(buf, buf_sz, "%s, i=%d", d->msg, static_cast<int>(d->msg_ind));
    } else {
      flag = snprintf(buf, buf_sz, "%s", d->msg);
    }
    // Check if error
    if (flag < 0) {
      d->status = QP_PRINTING_ERROR;
      return 1;
    }
  }
  // Successful return
  return 0;
}

#endif  // CASADI_PRINTF

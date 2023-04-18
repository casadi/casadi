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


// C-REPLACE "fmin" "casadi_fmin"
// C-REPLACE "fmax" "casadi_fmax"
// C-REPLACE "std::numeric_limits<T1>::min()" "casadi_real_min"
// C-REPLACE "std::numeric_limits<T1>::infinity()" "casadi_inf"
// C-REPLACE "static_cast<int>" "(int) "
// SYMBOL "ipqp_prob"
template<typename T1>
struct casadi_ipqp_prob {
  // Dimensions
  casadi_int nx, na, nz;
  // Smallest nonzero number
  T1 dmin;
  // Infinity
  T1 inf;
  // Maximum number of iterations
  casadi_int max_iter;
  // Error tolerance
  T1 pr_tol, du_tol, co_tol, mu_tol;
};
// C-REPLACE "casadi_ipqp_prob<T1>" "struct casadi_ipqp_prob"

// SYMBOL "ipqp_setup"
template<typename T1>
void casadi_ipqp_setup(casadi_ipqp_prob<T1>* p, casadi_int nx, casadi_int na) {
  p->nx = nx;
  p->na = na;
  p->nz = nx + na;
  p->dmin = std::numeric_limits<T1>::min();
  p->inf = std::numeric_limits<T1>::infinity();
  p->max_iter = 100;
  p->pr_tol = 1e-8;
  p->du_tol = 1e-8;
  p->co_tol = 1e-8;
  p->mu_tol = 1e-8;
}

// SYMBOL "ipqp_flag_t"
typedef enum {
  IPQP_SUCCESS,
  IPQP_MAX_ITER,
  IPQP_NO_SEARCH_DIR,
  IPQP_MV_ERROR,
  IPQP_FACTOR_ERROR,
  IPQP_SOLVE_ERROR,
  IPQP_PROGRESS_ERROR
} casadi_ipqp_flag_t;

// SYMBOL "ipqp_task_t"
typedef enum {
  IPQP_MV,
  IPQP_PROGRESS,
  IPQP_FACTOR,
  IPQP_SOLVE} casadi_ipqp_task_t;

// SYMBOL "ipqp_task_t"
typedef enum {
  IPQP_RESET,
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
  // QP data
  const T1 *g;
  // Number of finite constraints
  casadi_int n_con;
  // Solver status
  casadi_ipqp_flag_t status;
  // User task
  casadi_ipqp_task_t task;
  // Next step
  casadi_ipqp_next_t next;
  // Linear system
  T1* linsys;
  // Message buffer
  const char *msg;
  // Complementarity measure
  T1 mu;
  // Stepsize
  T1 tau;
  // Primal and dual error, complementarity error, corresponding index
  T1 pr, du, co;
  casadi_int ipr, idu, ico;
  // Iteration
  casadi_int iter;
  // Bounds
  T1 *lbz, *ubz;
  // Current solution
  T1 *z, *lam, *lam_lbz, *lam_ubz;
  // Step
  T1 *dz, *dlam, *dlam_lbz, *dlam_ubz;
  // Residual
  T1 *rz, *rlam, *rlam_lbz, *rlam_ubz;
  // Diagonal entries
  T1* D;
  // Scaling factor
  T1* S;
  // Inverse of margin to bounds (0 if no bound)
  T1 *dinv_lbz, *dinv_ubz;
};
// C-REPLACE "casadi_ipqp_data<T1>" "struct casadi_ipqp_data"

// SYMBOL "ipqp_sz_w"
template<typename T1>
casadi_int casadi_ipqp_sz_w(const casadi_ipqp_prob<T1>* p) {
  // Return value
  casadi_int sz_w = 0;
  // Persistent work vectors
  sz_w += p->nz; // lbz
  sz_w += p->nz; // ubz
  sz_w += p->nz; // z
  sz_w += p->nz; // lam
  sz_w += p->nz; // lam_lbz
  sz_w += p->nz; // lam_ubz
  sz_w += p->nz; // dz
  sz_w += p->nz; // dlam
  sz_w += p->nz; // dlam_lbz
  sz_w += p->nz; // dlam_ubz
  sz_w += p->nz; // rz
  sz_w += p->nz; // rlam
  sz_w += p->nz; // rlam_lbz
  sz_w += p->nz; // rlam_ubz
  sz_w += p->nz; // D
  sz_w += p->nz; // S
  sz_w += p->nz; // dinv_lbz
  sz_w += p->nz; // dinv_ubz
  return sz_w;
}

// SYMBOL "ipqp_init"
template<typename T1>
void casadi_ipqp_init(casadi_ipqp_data<T1>* d, casadi_int** iw, T1** w) {
  const casadi_ipqp_prob<T1>* p = d->prob;
  // Assign memory
  d->lbz = *w; *w += p->nz;
  d->ubz = *w; *w += p->nz;
  d->z = *w; *w += p->nz;
  d->lam = *w; *w += p->nz;
  d->lam_lbz = *w; *w += p->nz;
  d->lam_ubz = *w; *w += p->nz;
  d->dz = *w; *w += p->nz;
  d->dlam = *w; *w += p->nz;
  d->dlam_lbz = *w; *w += p->nz;
  d->dlam_ubz = *w; *w += p->nz;
  d->rz = *w; *w += p->nz;
  d->rlam = *w; *w += p->nz;
  d->rlam_lbz = *w; *w += p->nz;
  d->rlam_ubz = *w; *w += p->nz;
  d->D = *w; *w += p->nz;
  d->S = *w; *w += p->nz;
  d->dinv_lbz = *w; *w += p->nz;
  d->dinv_ubz = *w; *w += p->nz;
  // New QP
  d->next = IPQP_RESET;
}

// SYMBOL "ipqp_bounds"
template<typename T1>
void casadi_ipqp_bounds(casadi_ipqp_data<T1>* d, const T1* g,
    const T1* lbx, const T1* ubx, const T1* lba, const T1* uba) {
  // Local variables
  const casadi_ipqp_prob<T1>* p = d->prob;
  // Pass pointer to linear term in objective
  d->g = g;
  // Pass bounds on z
  casadi_copy(lbx, p->nx, d->lbz);
  casadi_copy(lba, p->na, d->lbz + p->nx);
  casadi_copy(ubx, p->nx, d->ubz);
  casadi_copy(uba, p->na, d->ubz + p->nx);
}

// SYMBOL "ipqp_guess"
template<typename T1>
void casadi_ipqp_guess(casadi_ipqp_data<T1>* d,
    const T1* x0, const T1* lam_x0, const T1* lam_a0) {
  // Local variables
  const casadi_ipqp_prob<T1>* p = d->prob;
  // Pass initial guess
  casadi_copy(x0, p->nx, d->z);
  casadi_fill(d->z + p->nx, p->na, 0.);
  casadi_copy(lam_x0, p->nx, d->lam);
  casadi_copy(lam_a0, p->na, d->lam + p->nx);
  casadi_fill(d->lam_lbz, p->nz, 0.);
  casadi_fill(d->lam_ubz, p->nz, 0.);
}

// SYMBOL "ipqp_reset"
template<typename T1>
void casadi_ipqp_reset(casadi_ipqp_data<T1>* d) {
  // Local variables
  casadi_int k;
  T1 margin, mid;
  const casadi_ipqp_prob<T1>* p = d->prob;
  // Required margin to constraints
  margin = .1;
  // Reset constraint count
  d->n_con = 0;
  // Initialize constraints to zero
  for (k = p->nx; k < p->nz; ++k) d->z[k] = 0;
  // Find interior point
  for (k = 0; k < p->nz; ++k) {
    if (d->lbz[k] > -p->inf) {
      if (d->ubz[k] < p->inf) {
        // Both upper and lower bounds
        mid = .5 * (d->lbz[k] + d->ubz[k]);
        // Ensure margin to boundary, without crossing midpoint
        if (d->z[k] < mid) {
          d->z[k] = fmin(fmax(d->z[k], d->lbz[k] + margin), mid);
        } else if (d->z[k] > mid) {
          d->z[k] = fmax(fmin(d->z[k], d->ubz[k] - margin), mid);
        }
        if (d->ubz[k] > d->lbz[k] + p->dmin) {
          d->lam_lbz[k] = 1;
          d->lam_ubz[k] = 1;
          d->n_con += 2;
        }
      } else {
        // Only lower bound
        d->z[k] = fmax(d->z[k], d->lbz[k] + margin);
        d->lam_lbz[k] = 1;
        d->n_con++;
      }
    } else {
      if (d->ubz[k] < p->inf) {
        // Only upper bound
        d->z[k] = fmin(d->z[k], d->ubz[k] - margin);
        d->lam_ubz[k] = 1;
        d->n_con++;
      }
    }
  }
  // Clear residual
  casadi_clear(d->rz, p->nz);
  // Reset iteration counter
  d->iter = 0;
  // Reset iteration variables
  d->msg = 0;
  d->tau = -1;
}

// SYMBOL "ipqp_diag"
template<typename T1>
void casadi_ipqp_diag(casadi_ipqp_data<T1>* d) {
  // Local variables
  casadi_int k;
  const casadi_ipqp_prob<T1>* p = d->prob;
  // Diagonal entries corresponding to variables
  for (k = 0; k < p->nx; ++k) {
    if (d->ubz[k] <= d->lbz[k] + p->dmin) {
      // Fixed variable (eliminate)
      d->D[k] = -1;
    } else {
      d->D[k] = d->lam_lbz[k] * d->dinv_lbz[k]
        + d->lam_ubz[k] * d->dinv_ubz[k];
    }
  }
  // Diagonal entries corresponding to constraints
  for (; k < p->nz; ++k) {
    if (d->lbz[k] <= -p->inf && d->ubz[k] >= p->inf) {
      // Unconstrained (eliminate)
      d->D[k] = -1;
    } else if (d->ubz[k] <= d->lbz[k] + p->dmin) {
      // Equality constrained
      d->D[k] = 0;
    } else {
      d->D[k] = 1. / (d->lam_lbz[k] * d->dinv_lbz[k]
        + d->lam_ubz[k] * d->dinv_ubz[k]);
    }
  }
  // Scale diagonal entries
  for (k = 0; k < p->nz; ++k) {
    if (d->D[k] < 0) {
      // Eliminate
      d->S[k] = 0;
      d->D[k] = 1;
    } else {
      // Scale
      d->S[k] = fmin(1., std::sqrt(1. / d->D[k]));
      d->D[k] = fmin(1., d->D[k]);
    }
  }
}

// SYMBOL "ipqp_newiter"
template<typename T1>
int casadi_ipqp_newiter(casadi_ipqp_data<T1>* d) {
  // Local variables
  const casadi_ipqp_prob<T1>* p = d->prob;
  // Converged?
  if (d->pr < p->pr_tol && d->du < p->du_tol && d->co < p->co_tol
      && d->mu < p->mu_tol) {
    d->status = IPQP_SUCCESS;
    return 1;
  }
  // Max number of iterations reached
  if (d->iter >= p->max_iter) {
    d->status = IPQP_MAX_ITER;
    return 1;
  }
  // Start new iteration
  d->iter++;
  // Calculate diagonal entries and scaling factors
  casadi_ipqp_diag(d);
  // Success
  return 0;
}

// SYMBOL "ipqp_residual"
template<typename T1>
void casadi_ipqp_residual(casadi_ipqp_data<T1>* d) {
  // Local variables
  casadi_int k;
  T1 bdiff, viol;
  const casadi_ipqp_prob<T1>* p = d->prob;
  // Gradient of the Lagrangian
  casadi_axpy(p->nx, 1., d->g, d->rz);
  for (k = 0; k < p->nx; ++k) {
    if (d->ubz[k] <= d->lbz[k] + p->dmin) {
      // Fixed variable: Solve to get multiplier explicitly
      d->lam[k] = -d->rz[k];
      d->rz[k] = 0;
    } else {
      // Residual
      d->rz[k] += d->lam[k];
    }
  }
  // Constraint violation (only possible for linear constraints)
  d->ipr = -1;
  d->pr = 0;
  for (k = p->nx; k < p->nz; ++k) {
    if (d->lbz[k] <= -p->inf && d->ubz[k] >= p->inf) {
      // Unconstrained: Solve to get g explicitly
      d->z[k] = d->rz[k];
      d->lam[k] = 0;
    } else {
      // Check constraint violation
      if (d->rz[k] + d->pr < d->lbz[k]) {
        d->pr = d->lbz[k] - d->rz[k];
        d->ipr = k;
      } else if (d->rz[k] - d->pr > d->ubz[k]) {
        d->pr = d->rz[k] - d->ubz[k];
        d->ipr = k;
      }
    }
  }
  // Dual infeasibility
  d->idu = -1;
  d->du = 0;
  for (k = 0; k < p->nx; ++k) {
    if (fabs(d->rz[k]) > d->du) {
      d->du = fabs(d->rz[k]);
      d->idu = k;
    }
  }
  // Linear constraint
  casadi_axpy(p->na, -1., d->z + p->nx, d->rz + p->nx);
  // Multiplier consistency
  for (k = 0; k < p->nz; ++k) {
    if (d->ubz[k] <= d->lbz[k] + p->dmin) {
      // Fixed variable: Solve to get lam_lbz, lam_ubz
      d->lam_ubz[k] = fmax(d->lam[k], 0.);
      d->lam_lbz[k] = fmax(-d->lam[k], 0.);
      d->rlam[k] = 0;
    } else {
      // Residual
      d->rlam[k] = d->lam_ubz[k] - d->lam_lbz[k] - d->lam[k];
    }
  }
  // Complementarity conditions, mu
  d->mu = 0;
  d->ico = -1;
  d->co = 0;
  for (k = 0; k < p->nz; ++k) {
    // Lower bound
    if (d->lbz[k] > -p->inf && d->ubz[k] > d->lbz[k] + p->dmin) {
      // Inequality constraint
      bdiff = d->z[k] - d->lbz[k];
      d->mu += d->rlam_lbz[k] = d->lam_lbz[k] * bdiff;
      d->dinv_lbz[k] = 1. / bdiff;
      // Constraint violation
      viol = bdiff * fmax(-d->lam[k], 0.);
      if (viol > d->co) {
        d->co = viol;
        d->ico = k;
      }
    } else {
      // No bound or equality constraint
      d->rlam_lbz[k] = 0;
      d->dinv_lbz[k] = 0;
    }
    // Upper bound
    if (d->ubz[k] < p->inf && d->ubz[k] > d->lbz[k] + p->dmin) {
      // Inequality constraint
      bdiff = d->ubz[k] - d->z[k];
      d->mu += d->rlam_ubz[k] = d->lam_ubz[k] * bdiff;
      d->dinv_ubz[k] = 1. / bdiff;
      // Constraint violation
      viol = bdiff * fmax(d->lam[k], 0.);
      if (viol > d->co) {
        d->co = viol;
        d->ico = k;
      }
    } else {
      // No bound or equality constraint
      d->rlam_ubz[k] = 0;
      d->dinv_ubz[k] = 0;
    }
  }
  // Divide mu by total number of finite constraints
  if (d->n_con > 0) d->mu /= d->n_con;
}

// SYMBOL "ipqp_predictor_prepare"
template<typename T1>
void casadi_ipqp_predictor_prepare(casadi_ipqp_data<T1>* d) {
  // Local variables
  casadi_int k;
  const casadi_ipqp_prob<T1>* p = d->prob;
  // Store r_lam - dinv_lbz * rlam_lbz + dinv_ubz * rlam_ubz in dz
  casadi_copy(d->rlam, p->nz, d->dz);
  for (k=0; k<p->nz; ++k) d->dz[k] += d->dinv_lbz[k] * d->rlam_lbz[k];
  for (k=0; k<p->nz; ++k) d->dz[k] -= d->dinv_ubz[k] * d->rlam_ubz[k];
  // Finish calculating x-component of right-hand-side and store in dz[:nx]
  for (k=0; k<p->nx; ++k) d->dz[k] += d->rz[k];
  // Copy tilde{r}_lam to dlam[nx:] (needed to calculate step in g later)
  for (k=p->nx; k<p->nz; ++k) d->dlam[k] = d->dz[k];
  // Finish calculating g-component of right-hand-side and store in dz[nx:]
  for (k=p->nx; k<p->nz; ++k) {
    if (d->S[k] == 0.) {
      // Eliminate
      d->dz[k] = 0;
    } else {
      d->dz[k] *= d->D[k] / (d->S[k] * d->S[k]);
      d->dz[k] += d->rz[k];
    }
  }
  // Scale and negate right-hand-side
  for (k=0; k<p->nz; ++k) d->dz[k] *= -d->S[k];
  // dlam_lbz := -rlam_lbz, dlam_ubz := -rlam_ubz
  for (k=0; k<p->nz; ++k) d->dlam_lbz[k] = -d->rlam_lbz[k];
  for (k=0; k<p->nz; ++k) d->dlam_ubz[k] = -d->rlam_ubz[k];
  // dlam_x := rlam_x
  for (k=0; k<p->nx; ++k) d->dlam[k] = d->rlam[k];
  // Solve to get step
  d->linsys = d->dz;
}

// SYMBOL "ipqp_maxstep"
template<typename T1>
int casadi_ipqp_maxstep(casadi_ipqp_data<T1>* d, T1* alpha, casadi_int* ind) {
  // Local variables
  T1 test;
  casadi_int k, blocking_k;
  int flag;
  const casadi_ipqp_prob<T1>* p = d->prob;
  // Reset variables
  blocking_k = -1;
  flag = IPQP_NONE;
  // Maximum step size is 1
  *alpha = 1.;
  // Primal step
  for (k=0; k<p->nz; ++k) {
    if (d->dz[k] < 0 && d->lbz[k] > -p->inf) {
      if ((test = (d->lbz[k] - d->z[k]) / d->dz[k]) < *alpha) {
        *alpha = test;
        blocking_k = k;
        flag = IPQP_PRIMAL | IPQP_LOWER;
      }
    }
    if (d->dz[k] > 0 && d->ubz[k] < p->inf) {
      if ((test = (d->ubz[k] - d->z[k]) / d->dz[k]) < *alpha) {
        *alpha = test;
        blocking_k = k;
        flag = IPQP_PRIMAL | IPQP_UPPER;
      }
    }
  }
  // Dual step
  for (k=0; k<p->nz; ++k) {
    if (d->dlam_lbz[k] < 0.) {
      if ((test = -d->lam_lbz[k] / d->dlam_lbz[k]) < *alpha) {
        *alpha = test;
        blocking_k = k;
        flag = IPQP_DUAL | IPQP_LOWER;
      }
    }
    if (d->dlam_ubz[k] < 0.) {
      if ((test = -d->lam_ubz[k] / d->dlam_ubz[k]) < *alpha) {
        *alpha = test;
        blocking_k = k;
        flag = IPQP_DUAL | IPQP_UPPER;
      }
    }
  }
  // Return information about blocking constraints
  if (ind) *ind = blocking_k;
  return flag;
}

// SYMBOL "ipqp_predictor"
template<typename T1>
void casadi_ipqp_predictor(casadi_ipqp_data<T1>* d) {
  // Local variables
  casadi_int k;
  T1 t, alpha, sigma;
  const casadi_ipqp_prob<T1>* p = d->prob;
  // Scale results
  for (k=0; k<p->nz; ++k) d->dz[k] *= d->S[k];
  // Calculate step in z(g), lam(g)
  for (k=p->nx; k<p->nz; ++k) {
    if (d->S[k] == 0.) {
      // Eliminate
      d->dlam[k] = d->dz[k] = 0;
    } else {
      t = d->D[k] / (d->S[k] * d->S[k]) * (d->dz[k] - d->dlam[k]);
      d->dlam[k] = d->dz[k];
      d->dz[k] = t;
    }
  }
  // Finish calculation in dlam_lbz, dlam_ubz
  for (k=0; k<p->nz; ++k) {
    d->dlam_lbz[k] -= d->lam_lbz[k] * d->dz[k];
    d->dlam_lbz[k] *= d->dinv_lbz[k];
  }
  for (k=0; k<p->nz; ++k) {
    d->dlam_ubz[k] += d->lam_ubz[k] * d->dz[k];
    d->dlam_ubz[k] *= d->dinv_ubz[k];
  }
  // Finish calculation of dlam(x)
  for (k=0; k<p->nx; ++k) d->dlam[k] += d->dlam_ubz[k] - d->dlam_lbz[k];
  // Maximum primal and dual step
  (void)casadi_ipqp_maxstep(d, &alpha, 0);
  // Calculate sigma
  sigma = casadi_ipqp_sigma(d, alpha);
  // Prepare corrector step
  casadi_ipqp_corrector_prepare(d, -sigma * d->mu);
  // Solve to get step
  d->linsys = d->rz;
}

// SYMBOL "ipqp_step"
template<typename T1>
void casadi_ipqp_step(casadi_ipqp_data<T1>* d, T1 alpha_pr, T1 alpha_du) {
  // Local variables
  casadi_int k;
  const casadi_ipqp_prob<T1>* p = d->prob;
  // Primal step
  for (k=0; k<p->nz; ++k) d->z[k] += alpha_pr * d->dz[k];
  // Dual step
  for (k=0; k<p->nz; ++k) d->lam[k] += alpha_du * d->dlam[k];
  for (k=0; k<p->nz; ++k) d->lam_lbz[k] += alpha_du * d->dlam_lbz[k];
  for (k=0; k<p->nz; ++k) d->lam_ubz[k] += alpha_du * d->dlam_ubz[k];
}

// SYMBOL "ipqp_mu"
template<typename T1>
T1 casadi_ipqp_mu(casadi_ipqp_data<T1>* d, T1 alpha) {
  // Local variables
  T1 mu;
  casadi_int k;
  const casadi_ipqp_prob<T1>* p = d->prob;
  // Quick return if no inequalities
  if (d->n_con == 0) return 0;
  // Calculate projected mu (and save to sigma variable)
  mu = 0;
  for (k = 0; k < p->nz; ++k) {
    // Lower bound
    if (d->lbz[k] > -p->inf && d->ubz[k] > d->lbz[k] + p->dmin) {
      mu += (d->lam_lbz[k] + alpha * d->dlam_lbz[k])
        * (d->z[k] - d->lbz[k] + alpha * d->dz[k]);
    }
    // Upper bound
    if (d->ubz[k] < p->inf && d->ubz[k] > d->lbz[k] + p->dmin) {
      mu += (d->lam_ubz[k] + alpha * d->dlam_ubz[k])
        * (d->ubz[k] - d->z[k] - alpha * d->dz[k]);
    }
  }
  // Divide mu by total number of finite constraints
  mu /= d->n_con;
  return mu;
}

// SYMBOL "ipqp_sigma"
template<typename T1>
T1 casadi_ipqp_sigma(casadi_ipqp_data<T1>* d, T1 alpha) {
  // Local variables
  T1 sigma;
  // Quick return if no inequalities
  if (d->n_con == 0) return 0;
  // Calculate projected mu (and save to sigma variable)
  sigma = casadi_ipqp_mu(d, alpha);
  // Finish calculation of sigma := (mu_aff / mu)^3
  sigma /= d->mu;
  sigma *= sigma * sigma;
  return sigma;
}

// SYMBOL "ipqp_corrector_prepare"
template<typename T1>
void casadi_ipqp_corrector_prepare(casadi_ipqp_data<T1>* d, T1 shift) {
  // Local variables
  casadi_int k;
  const casadi_ipqp_prob<T1>* p = d->prob;
  // Modified residual in lam_lbz, lam_ubz
  for (k=0; k<p->nz; ++k) d->rlam_lbz[k] = d->dlam_lbz[k] * d->dz[k] + shift;
  for (k=0; k<p->nz; ++k) d->rlam_ubz[k] = -d->dlam_ubz[k] * d->dz[k] + shift;
  // Difference in tilde(r)_x, tilde(r)_lamg
  for (k=0; k<p->nz; ++k)
    d->rz[k] = d->dinv_lbz[k] * d->rlam_lbz[k]
      - d->dinv_ubz[k] * d->rlam_ubz[k];
  // Difference in tilde(r)_g
  for (k=p->nx; k<p->nz; ++k) {
    if (d->S[k] == 0.) {
      // Eliminate
      d->rlam[k] = d->rz[k] = 0;
    } else {
      d->rlam[k] = d->rz[k];
      d->rz[k] *= d->D[k] / (d->S[k] * d->S[k]);
    }
  }
  // Scale and negate right-hand-side
  for (k=0; k<p->nz; ++k) d->rz[k] *= -d->S[k];
}

// SYMBOL "ipqp_corrector"
template<typename T1>
void casadi_ipqp_corrector(casadi_ipqp_data<T1>* d) {
  // Local variables
  T1 t, mu_test, primal_slack, primal_step, dual_slack, dual_step, max_tau;
  casadi_int k;
  int flag;
  const casadi_ipqp_prob<T1>* p = d->prob;
  // Scale results
  for (k=0; k<p->nz; ++k) d->rz[k] *= d->S[k];
  // Calculate step in z(g), lam(g)
  for (k=p->nx; k<p->nz; ++k) {
    if (d->S[k] == 0.) {
      // Eliminate
      d->rlam[k] = d->rz[k] = 0;
    } else {
      t = d->D[k] / (d->S[k] * d->S[k]) * (d->rz[k] - d->rlam[k]);
      d->rlam[k] = d->rz[k];
      d->rz[k] = t;
    }
  }
  // Update step in dz, dlam
  for (k=0; k<p->nz; ++k) d->dz[k] += d->rz[k];
  for (k=p->nx; k<p->nz; ++k) d->dlam[k] += d->rlam[k];
  // Update step in lam_lbz
  for (k=0; k<p->nz; ++k) {
    t = d->dinv_lbz[k] * (-d->rlam_lbz[k] - d->lam_lbz[k] * d->rz[k]);
    d->dlam_lbz[k] += t;
    if (k<p->nx) d->dlam[k] -= t;
  }
  // Update step in lam_ubz
  for (k=0; k<p->nz; ++k) {
    t = d->dinv_ubz[k] * (-d->rlam_ubz[k] + d->lam_ubz[k] * d->rz[k]);
    d->dlam_ubz[k] += t;
    if (k<p->nx) d->dlam[k] += t;
  }
  // Find the largest step size, keeping track of blocking constraints
  flag = casadi_ipqp_maxstep(d, &max_tau, &k);
  // Handle blocking constraints using Mehrotra's heuristic
  if (flag == IPQP_NONE) {
    // No blocking constraints
    d->tau = 1;
  } else {
    // Calculate mu for maximum step
    mu_test = casadi_ipqp_mu(d, max_tau);
    // Get distance to constraints for blocking variable
    if (flag & IPQP_UPPER) {
      primal_slack = d->ubz[k] - d->z[k];
      primal_step = -d->dz[k];
      dual_slack = d->lam_ubz[k];
      dual_step = d->dlam_ubz[k];
    } else {
      primal_slack = d->z[k] - d->lbz[k];
      primal_step = d->dz[k];
      dual_slack = d->lam_lbz[k];
      dual_step = d->dlam_lbz[k];
    }
    // Mehrotra's heuristic as in in OOQP per communication with S. Wright
    if (flag & IPQP_PRIMAL) {
      d->tau = (0.01 * mu_test / (dual_slack + max_tau * dual_step)
        - primal_slack) / primal_step;
    } else {
      d->tau = (0.01 * mu_test / (primal_slack + max_tau * primal_step)
        - dual_slack) / dual_step;
    }
    d->tau = fmax(d->tau, 0.99 * max_tau);
  }
  // Take step
  casadi_ipqp_step(d, d->tau, d->tau);
  // Clear residual
  casadi_clear(d->rz, p->nz);
}

// SYMBOL "ipqp"
template<typename T1>
int casadi_ipqp(casadi_ipqp_data<T1>* d) {
  switch (d->next) {
    case IPQP_RESET:
      casadi_ipqp_reset(d);
      d->task = IPQP_MV;
      d->next = IPQP_RESIDUAL;
      return 1;
    case IPQP_RESIDUAL:
      // Calculate residual
      if (d->status == IPQP_MV_ERROR) break;
      casadi_ipqp_residual(d);
      d->task = IPQP_PROGRESS;
      d->next = IPQP_NEWITER;
      return 1;
    case IPQP_NEWITER:
      // New iteration
      if (d->status == IPQP_PROGRESS_ERROR) break;
      if (casadi_ipqp_newiter(d)) break;
      d->task = IPQP_FACTOR;
      d->next = IPQP_PREPARE;
      return 1;
    case IPQP_PREPARE:
      // Prepare predictor step
      if (d->status == IPQP_FACTOR_ERROR) break;
      casadi_ipqp_predictor_prepare(d);
      d->task = IPQP_SOLVE;
      d->next = IPQP_PREDICTOR;
      return 1;
    case IPQP_PREDICTOR:
      // Complete predictor step
      if (d->status == IPQP_SOLVE_ERROR) break;
      casadi_ipqp_predictor(d);
      d->task = IPQP_SOLVE;
      d->next = IPQP_CORRECTOR;
      return 1;
    case IPQP_CORRECTOR:
      // Complete predictor step
      if (d->status == IPQP_SOLVE_ERROR) break;
      casadi_ipqp_corrector(d);
      d->task = IPQP_MV;
      d->next = IPQP_RESIDUAL;
      return 1;
    default:
      break;
  }
  // Done iterating
  d->next = IPQP_RESET;
  return 0;
}

// SYMBOL "ipqp_return_status"
inline
const char* casadi_ipqp_return_status(casadi_ipqp_flag_t status) {
  switch (status) {
    case IPQP_SUCCESS: return "success";
    case IPQP_MAX_ITER: return "Maximum number of iterations reached";
    case IPQP_NO_SEARCH_DIR: return "Failed to calculate search direction";
    case IPQP_MV_ERROR: return "Matrix-vector evaluation error";
    case IPQP_FACTOR_ERROR: return "Linear solver factorization error";
    case IPQP_SOLVE_ERROR: return "Linear solver solution error";
    case IPQP_PROGRESS_ERROR: return "Printing error";
  }
  return 0;
}

// SYMBOL "ipqp_solution"
template<typename T1>
void casadi_ipqp_solution(casadi_ipqp_data<T1>* d, T1* x, T1* lam_x, T1* lam_a) {
  // Local variables
  const casadi_ipqp_prob<T1>* p = d->prob;
  // Copy solution
  casadi_copy(d->z, p->nx, x);
  casadi_copy(d->lam, p->nx, lam_x);
  casadi_copy(d->lam + p->nx, p->na, lam_a);
}

// SYMBOL "ipqp_print_header"
template<typename T1>
int casadi_ipqp_print_header(casadi_ipqp_data<T1>* d, char* buf, size_t buf_sz) {
#ifdef CASADI_SNPRINTF
  int flag;
  // Print to string
  flag = CASADI_SNPRINTF(buf, buf_sz, "%5s %9s %9s %5s %9s %5s "
          "%9s %5s %9s %4s",
          "Iter", "mu", "|pr|", "con", "|du|", "var", "|co|", "con",
          "last_tau", "Note");
  // Check if error
  if (flag < 0) {
    d->status = IPQP_PROGRESS_ERROR;
    return 1;
  }
#else
  if (buf_sz) buf[0] = '\0';
#endif
  // Successful return
  return 0;
}

// SYMBOL "ipqp_print_iteration"
template<typename T1>
int casadi_ipqp_print_iteration(casadi_ipqp_data<T1>* d, char* buf, int buf_sz) {
#ifdef CASADI_SNPRINTF
  int flag;
  // Print iteration data without note to string
  flag = CASADI_SNPRINTF(buf, buf_sz,
    "%5d %9.2g %9.2g %5d %9.2g %5d %9.2g %5d %9.2g  ",
    static_cast<int>(d->iter), d->mu,
    d->pr, static_cast<int>(d->ipr),
    d->du, static_cast<int>(d->idu),
    d->co, static_cast<int>(d->ico),
    d->tau);
  // Check if error
  if (flag < 0) {
    d->status = IPQP_PROGRESS_ERROR;
    return 1;
  }
  // Rest of buffer reserved for iteration note
  buf += flag;
  buf_sz -= flag;
  // Print iteration note, if any
  if (d->msg) {
    flag = CASADI_SNPRINTF(buf, buf_sz, "%s", d->msg);
    // Check if error
    if (flag < 0) {
      d->status = IPQP_PROGRESS_ERROR;
      return 1;
    }
  }
#else
  if (buf_sz) buf[0] = '\0';
#endif
  // Successful return
  return 0;
}

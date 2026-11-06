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

// C-REPLACE "casadi_qp_prob<T1>" "struct casadi_qp_prob"
// C-REPLACE "casadi_qp_data<T1>" "struct casadi_qp_data"

// C-REPLACE "reinterpret_cast<int**>" "(int**) "
// C-REPLACE "reinterpret_cast<int*>" "(int*) "
// C-REPLACE "const_cast<DAQPSettings*>" "(DAQPSettings*) "
// C-REPLACE "static_cast<const T1*>" "(const casadi_real*) "
// C-REPLACE "const_cast<T1*>" "(casadi_real*) "

template<typename T1>
struct casadi_daqp_prob {
  const casadi_qp_prob<T1>* qp;

  DAQPSettings settings;
  int warm_start;
  int warm_start_previous;
  const int *integrality;
};
// C-REPLACE "casadi_daqp_prob<T1>" "struct casadi_daqp_prob"

// SYMBOL "daqp_setup"
template<typename T1>
void casadi_daqp_setup(casadi_daqp_prob<T1>* p) {

}



// SYMBOL "daqp_data"
template<typename T1>
struct casadi_daqp_data {
  // Problem structure
  const casadi_daqp_prob<T1>* prob;
  // Problem structure
  casadi_qp_data<T1>* qp;

  DAQPWorkspace work;
  DAQPProblem daqp;
  DAQPResult res;
  // Settings and cached problem data must outlive a generated evaluation.
  DAQPSettings settings;
  T1* cache;
  int* sense;
  int initialized;
  int workspace_reused;
  int update_mask;

  int return_status;
  int nodecount;
  int bnb_itercount;
};
// C-REPLACE "casadi_daqp_data<T1>" "struct casadi_daqp_data"

// C-REPLACE "DAQPWorkspace()" "(DAQPWorkspace) {0}"
// C-REPLACE "reinterpret_cast<T1*>" "(casadi_real*) "

// SYMBOL "daqp_init_mem"
template<typename T1>
int casadi_daqp_init_mem(casadi_daqp_data<T1>* d) {
  d->work = DAQPWorkspace();
  d->cache = 0;
  d->sense = 0;
  d->initialized = 0;
  d->workspace_reused = d->update_mask = 0;
  d->return_status = d->nodecount = d->bnb_itercount = 0;
  return 0;
}

// SYMBOL "daqp_free_mem"
template<typename T1>
void casadi_daqp_free_mem(casadi_daqp_data<T1>* d) {
  d->work.settings = 0;
  free_daqp_workspace(&d->work);
  free_daqp_ldp(&d->work);
  free(d->cache);
  free(d->sense);
  d->cache = 0;
  d->sense = 0;
  d->initialized = 0;
  d->work = DAQPWorkspace();
}

// SYMBOL "daqp_work"
template<typename T1>
void casadi_daqp_work(const casadi_daqp_prob<T1>* p, casadi_int* sz_arg, casadi_int* sz_res, casadi_int* sz_iw, casadi_int* sz_w) {
  casadi_qp_work(p->qp, sz_arg, sz_res, sz_iw, sz_w);
  // Only multipliers need a per-evaluation scratch buffer.
  *sz_w = p->qp->nz;
  *sz_iw = 0;
}

// SYMBOL "daqp_set_work"
template<typename T1>
void casadi_daqp_set_work(casadi_daqp_data<T1>* d, const T1*** arg, T1*** res,
    casadi_int** iw, T1** w) {
  d->res.lam = *w; *w += d->prob->qp->nz;
}

// C-REPLACE "fabs" "casadi_fabs"

// C-REPLACE "SOLVER_RET_UNKNOWN" "1"
// C-REPLACE "SOLVER_RET_INFEASIBLE" "4"
// C-REPLACE "SOLVER_RET_SUCCESS" "0"
// C-REPLACE "SOLVER_RET_LIMITED" "2"

// SYMBOL "daqp_update_data"
template<typename T1>
int casadi_daqp_update_data(const T1* src, casadi_int n, T1* dst) {
  casadi_int i;
  int changed = 0;
  for (i=0; i<n; ++i) {
    T1 value = src ? src[i] : 0;
    if (dst[i] != value) {
      changed = 1;
      dst[i] = value;
    }
  }
  return changed;
}

// SYMBOL "daqp_update_matrix"
template<typename T1>
int casadi_daqp_update_matrix(const T1* src, const casadi_int* sp, T1* dst,
    int transpose) {
  casadi_int col, k;
  const casadi_int* colind = sp+2;
  const casadi_int* row = sp+3+sp[1];
  // The sparsity pattern is fixed. Compare only stored entries against the
  // persistent dense matrix; no temporary dense matrix or sparse cache is needed.
  for (col=0; col<sp[1]; ++col) {
    for (k=colind[col]; k<colind[col+1]; ++k) {
      casadi_int index = transpose ? col+row[k]*sp[1] : row[k]+col*sp[0];
      if (dst[index] != (src ? src[k] : 0)) {
        casadi_densify(src, sp, dst, transpose);
        return 1;
      }
    }
  }
  return 0;
}

// SYMBOL "daqp_alloc_data"
template<typename T1>
int casadi_daqp_alloc_data(casadi_daqp_data<T1>* d) {
  const casadi_qp_prob<T1>* p_qp = d->prob->qp;
  T1* buffer;
  if (!d->cache) {
    d->cache = reinterpret_cast<T1*>(calloc(p_qp->nx*p_qp->nx +
      p_qp->na*p_qp->nx + p_qp->nx + 2*p_qp->nz, sizeof(T1)));
    d->sense = reinterpret_cast<int*>(calloc(p_qp->nz, sizeof(int)));
    if (!d->cache || !d->sense) {
      casadi_daqp_free_mem(d);
      return 1;
    }
  }
  buffer = d->cache;
  d->daqp.H = buffer; buffer += p_qp->nx*p_qp->nx;
  d->daqp.A = buffer; buffer += p_qp->na*p_qp->nx;
  d->daqp.f = buffer; buffer += p_qp->nx;
  d->daqp.blower = buffer; buffer += p_qp->nz;
  d->daqp.bupper = buffer;
  d->daqp.sense = d->sense;
  return 0;
}

// SYMBOL "daqp_update_inputs"
template<typename T1>
int casadi_daqp_update_inputs(casadi_daqp_data<T1>* d, int cold_start) {
  int mask = DAQP_UPDATE_unconstrained;
  const casadi_daqp_prob<T1>* p = d->prob;
  const casadi_qp_prob<T1>* p_qp = p->qp;
  casadi_qp_data<T1>* d_qp = d->qp;
  // Recheck bounds after resetting sense, even when their values are unchanged.
  if (cold_start) mask |= DAQP_UPDATE_sense | DAQP_UPDATE_d | DAQP_UPDATE_eliminate;
  if (casadi_daqp_update_matrix(d_qp->h, p_qp->sp_h, d->daqp.H, 0)) {
    mask |= DAQP_UPDATE_Rinv;
  }
  if (casadi_daqp_update_matrix(d_qp->a, p_qp->sp_a, d->daqp.A, 1)) {
    mask |= DAQP_UPDATE_M;
  }
  if (casadi_daqp_update_data(d_qp->g, p_qp->nx, d->daqp.f)) {
    mask |= DAQP_UPDATE_v;
  }
  if (casadi_daqp_update_data(d_qp->lbx, p_qp->nx, d->daqp.blower) |
      casadi_daqp_update_data(d_qp->lba, p_qp->na, d->daqp.blower+p_qp->nx)) {
    mask |= DAQP_UPDATE_d;
  }
  if (casadi_daqp_update_data(d_qp->ubx, p_qp->nx, d->daqp.bupper) |
      casadi_daqp_update_data(d_qp->uba, p_qp->na, d->daqp.bupper+p_qp->nx)) {
    mask |= DAQP_UPDATE_d;
  }
  return mask;
}

// SYMBOL "daqp_update_sense"
template<typename T1>
int casadi_daqp_update_sense(casadi_daqp_data<T1>* d, int mask, int cold_start) {
  casadi_int i;
  const casadi_qp_prob<T1>* p_qp = d->prob->qp;
  // DAQP skips immutable constraints in its bound check. Restore input sense
  // on bound changes so equalities can become inequalities again.
  if (mask & DAQP_UPDATE_d) mask |= DAQP_UPDATE_sense;
  // Retain only active-set flags; DAQP derives immutable flags from fresh bounds.
  if (d->initialized && !cold_start) {
    daqp_eq_restore(&d->work);
    for (i=0; i<p_qp->nz; ++i) {
      d->sense[i] |= d->work.sense[i] & (DAQP_ACTIVE | DAQP_LOWER);
    }
  }
  // Rebuild the active-set factorization after H or A changes.
  if (mask & (DAQP_UPDATE_Rinv | DAQP_UPDATE_M)) mask |= DAQP_UPDATE_sense;
  if (mask & DAQP_UPDATE_sense) mask |= DAQP_UPDATE_d;
  // The new equality set is inferred inside DAQP. Rebuild full constraints
  // before that inference when a previous reduction may have omitted rows.
  if ((mask & DAQP_UPDATE_sense) && d->work.eq) mask |= DAQP_UPDATE_M;
  return mask;
}

// SYMBOL "daqp_check_inputs"
template<typename T1>
int casadi_daqp_check_inputs(casadi_daqp_data<T1>* d) {
  casadi_int i;
  const casadi_daqp_prob<T1>* p = d->prob;
  const casadi_qp_prob<T1>* p_qp = p->qp;
  casadi_qp_data<T1>* d_qp = d->qp;
  // DAQP classifies equalities from the bounds during its update.
  for (i=0; i<p_qp->nz; ++i) d->sense[i] = 0;

  if (p->integrality) {
    for (casadi_int j = 0; j < p_qp->nx; ++j) {
      if (!p->integrality[j]) continue;

      double lb = d_qp->lbx[j];
      double ub = d_qp->ubx[j];

      int binary = fabs(lb - 0.0) < 1e-9 && fabs(ub - 1.0) < 1e-9;

      // Inputs have already been copied into the cache. A rejected call must
      // not leave those values ahead of the workspace's matrix factors.
      if (!binary) casadi_daqp_free_mem(d);
      casadi_assert(binary, "DAQP only supports binary variables with bounds [0,1], but variable " + str(j) + " has bounds [" + str(lb) + ", " + str(ub) + "]."); // NOLINT(whitespace/line_length)
      if (!binary) return 1;

      d->sense[j] |= DAQP_BINARY;  // mark as binary
    }
  }
  // DAQP 0.10 can omit an immutable zero row without checking consistency.
  for (i=0; i<p_qp->na; ++i) {
    if (d->daqp.blower[p_qp->nx+i] <= p->settings.primal_tol &&
        d->daqp.bupper[p_qp->nx+i] >= -p->settings.primal_tol) continue;
    casadi_int j;
    for (j=0; j<p_qp->nx; ++j) {
      if (d->daqp.A[i*p_qp->nx+j] != 0) break;
    }
    if (j==p_qp->nx) {
      d->return_status = DAQP_EXIT_INFEASIBLE;
      d_qp->unified_return_status = SOLVER_RET_INFEASIBLE;
      casadi_daqp_free_mem(d);
      return 1;
    }
  }
  return 0;
}

// SYMBOL "daqp_init_start"
template<typename T1>
void casadi_daqp_init_start(casadi_daqp_data<T1>* d, int cold_start) {
  const casadi_daqp_prob<T1>* p = d->prob;
  const casadi_qp_prob<T1>* p_qp = p->qp;
  casadi_qp_data<T1>* d_qp = d->qp;
  if (d->initialized && cold_start) {
    daqp_eq_restore(&d->work);
    reset_daqp_workspace(&d->work);
    // Proximal iterations use the primal vector as their initial center.
    casadi_copy(static_cast<const T1*>(0), p_qp->nx, d->work.u);
    casadi_copy(static_cast<const T1*>(0), p_qp->nx, d->work.x);
    casadi_copy(static_cast<const T1*>(0), p_qp->nx, d->work.xold);
    d->work.state &= ~DAQP_STATE_INCUMBENT;
    // DAQP otherwise implicitly loads the previous root active set.
    if (d->work.bnb) d->work.bnb->n_root_WS = 0;
  }
  if (p->warm_start) {
    if (d_qp->lam_x0 || d_qp->lam_a0) {
      casadi_copy(d_qp->lam_x0, p_qp->nx, d->res.lam);
      casadi_copy(d_qp->lam_a0, p_qp->na, d->res.lam+p_qp->nx);
      daqp_dual_init_active(&d->daqp, d->res.lam);
    } else if (d_qp->x0) {
      daqp_primal_init_active(&d->daqp, const_cast<T1*>(d_qp->x0));
    }
  }
}

// SYMBOL "daqp_store_result"
// DAQP exit flags (see daqp/constants.h):
//   2: soft optimal, 1: optimal
//  -1: infeasible, -2: cycle, -3: unbounded, -4: iteration limit
//  -5: nonconvex, -6: overdetermined initial active set
//  -7: time limit, -8: unsupported problem
template<typename T1>
void casadi_daqp_store_result(casadi_daqp_data<T1>* d) {
  const casadi_qp_prob<T1>* p_qp = d->prob->qp;
  casadi_qp_data<T1>* d_qp = d->qp;
  casadi_copy(d->res.lam, p_qp->nx, d_qp->lam_x);
  casadi_copy(d->res.lam+p_qp->nx, p_qp->na, d_qp->lam_a);
  if (d->work.bnb) {
    d->nodecount = d->work.bnb->nodecount;
    d->bnb_itercount = d->work.bnb->itercount;
  } else {
    d->nodecount = 0;
    d->bnb_itercount = 0;
  }
  if (d_qp->f) *d_qp->f = d->res.fval;
  if (d->res.exitflag < 0) casadi_daqp_free_mem(d);

  d->return_status = d->res.exitflag;
  d_qp->iter_count = d->res.iter;
  d_qp->success = d->res.exitflag == DAQP_EXIT_OPTIMAL;
  if (d_qp->success) {
    d_qp->unified_return_status = SOLVER_RET_SUCCESS;
  } else if (d->res.exitflag == DAQP_EXIT_ITERLIMIT ||
             d->res.exitflag == DAQP_EXIT_TIMELIMIT) {
    d_qp->unified_return_status = SOLVER_RET_LIMITED;
  } else if (d->res.exitflag == DAQP_EXIT_INFEASIBLE) {
    d_qp->unified_return_status = SOLVER_RET_INFEASIBLE;
  }
}

// SYMBOL "daqp_solve"
template<typename T1>
int casadi_daqp_solve(casadi_daqp_data<T1>* d, const double** arg, double** res,
    casadi_int* iw, double* w) {
  int flag, mask, cold_start;
  const casadi_daqp_prob<T1>* p = d->prob;
  casadi_qp_data<T1>* d_qp = d->qp;

  d->return_status = DAQP_EXIT_UNSUPPORTED;
  d->nodecount = d->bnb_itercount = 0;
  d->workspace_reused = d->update_mask = 0;
  d_qp->success = 0;
  d_qp->iter_count = -1;
  d_qp->unified_return_status = SOLVER_RET_UNKNOWN;
  d->res.exitflag = d->res.iter = 0;
  d->daqp.problem_type = d->daqp.nh = 0;
  d->daqp.break_points = 0;
  d->daqp.n = p->qp->nx;
  d->daqp.m = p->qp->nz;
  d->daqp.ms = p->qp->nx;
  d->settings = p->settings;
  d->work.settings = &d->settings;
  d->res.x = d_qp->x;

  if (casadi_daqp_alloc_data(d)) return 1;
  cold_start = !p->warm_start_previous ||
    (p->warm_start && (d_qp->x0 || d_qp->lam_x0 || d_qp->lam_a0));
  mask = casadi_daqp_update_inputs(d, cold_start);
  if (casadi_daqp_check_inputs(d)) return 1;
  mask = casadi_daqp_update_sense(d, mask, cold_start);
  casadi_daqp_init_start(d, cold_start);

  d->workspace_reused = d->initialized;
  d->update_mask = mask;
  if (d->initialized) {
    // Equality elimination saves Hessian pointers invalidated by H updates.
    if (mask & DAQP_UPDATE_Rinv) free_daqp_eq(&d->work);
    flag = daqp_update_ldp(mask, &d->work, &d->daqp);
  } else {
    flag = setup_daqp_main(&d->daqp, &d->work, &d->res.setup_time, mask);
  }
  if (flag<0) {
    d->return_status = flag;
    casadi_daqp_free_mem(d);
    return 1;
  }
  d->initialized = 1;
  if (p->warm_start && d_qp->x0) {
    daqp_set_primal_start(&d->work, const_cast<T1*>(d_qp->x0));
  }
  daqp_solve(&d->res, &d->work);
  casadi_daqp_store_result(d);
  return 0;
}

//
//    MIT No Attribution
//
//    Copyright (C) 2010-2026 Joel Andersson, Joris Gillis, Moritz Diehl, KU Leuven.
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
// C-REPLACE "casadi_socp_prob<T1>" "struct casadi_socp_prob"
// C-REPLACE "casadi_socp_data<T1>" "struct casadi_socp_data"

// C-REPLACE "reinterpret_cast<int**>" "(int**) "
// C-REPLACE "reinterpret_cast<int*>" "(int*) "
// C-REPLACE "reinterpret_cast<char*>" "(char*) "
// C-REPLACE "const_cast<int*>" "(int*) "
// C-REPLACE "static_cast<int>" "(int) "
// C-REPLACE "static_cast<MSKboundkeye>" "(MSKboundkeye) "
// C-REPLACE "static_cast<MSKvariabletypee>" "(MSKvariabletypee) "

template<typename T1>
struct casadi_mosek_prob {
  const casadi_qp_prob<T1>* qp;

  // Linear constraint matrix in CSC form (cast to int, since Mosek takes int)
  const int *colinda, *rowa;
  // Quadratic objective in lower-triangular triplet form.  Mosek's
  // MSK_putqobj wants entries with row >= col (lower triangle), and
  // interprets the objective as 0.5 x'Qx + c'x with Q symmetric --
  // matching CasADi's convention.  Pass the H[i,j] values as-is.
  const int *qobj_row, *qobj_col;
  // qobj_nz_idx[k] = index into d_qp->h[] of the k-th lower-tri entry
  const int *qobj_nz_idx;
  int nquad;
  // Discrete variable flags ('I' for integer, 'C' for continuous), or NULL
  const char *coltype;

  // SOCP support.  Null if no SOC cones.
  const casadi_socp_prob<T1> *socp;
};
// C-REPLACE "casadi_mosek_prob<T1>" "struct casadi_mosek_prob"

// SYMBOL "mosek_setup"
template<typename T1>
void casadi_mosek_setup(casadi_mosek_prob<T1>* p) {

}

// SYMBOL "mosek_data"
template<typename T1>
struct casadi_mosek_data {
  // Problem structure
  const casadi_mosek_prob<T1>* prob;
  // QP runtime data
  casadi_qp_data<T1>* qp;

  // Mosek environment + task handles.  The env is a per-instance handle
  // (created in init_mem) -- Mosek allows multiple envs per process.
  MSKenv_t  env;
  MSKtask_t task;

  // Workspace owned by the C runtime (sized via casadi_mosek_work)
  int      *col_idx;     // 0..nx-1, scratch for MSK_putvartypelist
  int      *vtype;       // var type code int[nx], scratch for MSK_putvartypelist
  T1       *qobj_val;    // scaled Q triplet values, length nquad
  T1       *dual_scratch; // length nx, scratch for sux retrieval
  int      *cone_idx;    // scratch for MSK_appendcone variable list (length max_block)

  // SOCP runtime data (only meaningful when prob->socp != NULL)
  casadi_socp_data<T1> socp;

  // 0 until the first solve has appended vars/cons/cones to the task;
  // gates the one-shot model-building path inside solve.
  int model_built;

  // Status / statistics
  int return_status;
  int prob_status;
  int sol_status;
  int simplex_iter;
  int barrier_iter;
  int mip_nodes;
  T1  obj_val;
};
// C-REPLACE "casadi_mosek_data<T1>" "struct casadi_mosek_data"

// SYMBOL "mosek_init_mem"
template<typename T1>
int casadi_mosek_init_mem(casadi_mosek_data<T1>* d) {
  d->env = 0;
  d->task = 0;
  d->model_built = 0;
  if (MSK_makeenv(&d->env, 0) != MSK_RES_OK) return 1;
  if (MSK_maketask(d->env, 0, 0, &d->task) != MSK_RES_OK) return 1;
  return 0;
}

// SYMBOL "mosek_free_mem"
template<typename T1>
void casadi_mosek_free_mem(casadi_mosek_data<T1>* d) {
  if (d->task) {
    MSK_deletetask(&d->task);
    d->task = 0;
  }
  if (d->env) {
    MSK_deleteenv(&d->env);
    d->env = 0;
  }
}

// SYMBOL "mosek_work"
template<typename T1>
void casadi_mosek_work(const casadi_mosek_prob<T1>* p,
    casadi_int* sz_arg, casadi_int* sz_res, casadi_int* sz_iw, casadi_int* sz_w) {
  casadi_qp_work(p->qp, sz_arg, sz_res, sz_iw, sz_w);
  // int scratch (overlaid on iw slots)
  *sz_iw += p->qp->nx;               // col_idx
  *sz_iw += p->qp->nx;               // vtype
  // double scratch
  *sz_w  += p->nquad;                // qobj_val
  *sz_w  += p->qp->nx;               // dual_scratch (sux retrieval)

  if (p->socp) {
    casadi_socp_work(p->socp, sz_iw, sz_w);
    // Cone variable index scratch, sized to the largest cone block.
    *sz_iw += p->socp->max_block;
  }
}

// SYMBOL "mosek_init"
template<typename T1>
void casadi_mosek_init(casadi_mosek_data<T1>* d,
    const T1*** arg, T1*** res, casadi_int** iw, T1** w) {
  const casadi_mosek_prob<T1>* p = d->prob;
  d->col_idx      = reinterpret_cast<int*>(*iw);  *iw += p->qp->nx;
  d->vtype        = reinterpret_cast<int*>(*iw);  *iw += p->qp->nx;
  d->qobj_val     = *w; *w += p->nquad;
  d->dual_scratch = *w; *w += p->qp->nx;
  d->cone_idx     = 0;

  if (p->socp) {
    d->socp.prob = p->socp;
    casadi_socp_init(&d->socp, iw, w);
    d->cone_idx = reinterpret_cast<int*>(*iw);  *iw += p->socp->max_block;
  }
}


// C-REPLACE "SOLVER_RET_SUCCESS" "0"
// C-REPLACE "SOLVER_RET_LIMITED" "2"
// C-REPLACE "std::numeric_limits<T1>::infinity()" "casadi_inf"


// SYMBOL "mosek_solve"
template<typename T1>
int casadi_mosek_solve(casadi_mosek_data<T1>* d,
    const double** arg, double** res, casadi_int* iw, double* w) {

  const casadi_mosek_prob<T1>* p = d->prob;
  const casadi_qp_prob<T1>* p_qp = p->qp;
  casadi_qp_data<T1>* d_qp = d->qp;

  int i, k, j;
  int has_mip = p->coltype ? 1 : 0;
  MSKtask_t task = d->task;
  MSKrescodee r;

  // Total variable count: nx + (lifted, if SOCP)
  int n_lifted = (p->socp) ? static_cast<int>(p->socp->n_lifted) : 0;
  int nvar = p_qp->nx + n_lifted;
  int ncon = p_qp->na + ((p->socp) ? static_cast<int>(p->socp->n_eq) : 0);

  if (!d->model_built) {
    if (MSK_appendvars(task, nvar) != MSK_RES_OK) return 1;
    if (MSK_appendcons(task, ncon) != MSK_RES_OK) return 1;
    MSK_putobjsense(task, MSK_OBJECTIVE_SENSE_MINIMIZE);
  }

  // Linear objective coefficient: c = g (for original variables)
  for (j = 0; j < p_qp->nx; ++j) {
    MSK_putcj(task, j, d_qp->g[j]);
  }
  // Lifted variables get c = 0 (handled by default; appendvars zero-inits).

  // Variable bounds: original variables use lbx/ubx
  for (j = 0; j < p_qp->nx; ++j) {
    T1 lo = d_qp->lbx[j], up = d_qp->ubx[j];
    int lo_inf = (lo <= -std::numeric_limits<T1>::infinity());
    int up_inf = (up >=  std::numeric_limits<T1>::infinity());
    MSKboundkeye bk;
    if (lo_inf && up_inf) { bk = MSK_BK_FR; }
    else if (lo_inf)      { bk = MSK_BK_UP; lo = -std::numeric_limits<T1>::infinity(); }
    else if (up_inf)      { bk = MSK_BK_LO; up =  std::numeric_limits<T1>::infinity(); }
    else if (lo == up)    { bk = MSK_BK_FX; }
    else                  { bk = MSK_BK_RA; }
    MSK_putvarbound(task, j, bk, lo, up);
  }

  // Linear A: column-slice from CSC
  for (j = 0; j < p_qp->nx; ++j) {
    int kbeg = p->colinda[j];
    int kend = p->colinda[j + 1];
    if (kend > kbeg) {
      // Allow MSK_putaijlist or per-column put.  We use MSK_putacol.
      MSK_putacol(task, j, kend - kbeg, p->rowa + kbeg, d_qp->a + kbeg);
    }
  }

  // Constraint bounds: [con bk][lba/uba]
  for (i = 0; i < p_qp->na; ++i) {
    T1 lo = d_qp->lba[i], up = d_qp->uba[i];
    int lo_inf = (lo <= -std::numeric_limits<T1>::infinity());
    int up_inf = (up >=  std::numeric_limits<T1>::infinity());
    MSKboundkeye bk;
    if (lo_inf && up_inf) { bk = MSK_BK_FR; }
    else if (lo_inf)      { bk = MSK_BK_UP; lo = -std::numeric_limits<T1>::infinity(); }
    else if (up_inf)      { bk = MSK_BK_LO; up =  std::numeric_limits<T1>::infinity(); }
    else if (lo == up)    { bk = MSK_BK_FX; }
    else                  { bk = MSK_BK_RA; }
    MSK_putconbound(task, i, bk, lo, up);
  }

  // Quadratic objective: lower triangle, values as-is.  Mosek interprets
  // the qobj as 0.5 x'Qx with Q symmetric -- matches CasADi's H.
  for (k = 0; k < p->nquad; ++k) {
    d->qobj_val[k] = d_qp->h[p->qobj_nz_idx[k]];
  }
  if (p->nquad > 0) {
    if (MSK_putqobj(task, p->nquad, p->qobj_row, p->qobj_col, d->qobj_val)
        != MSK_RES_OK) return 1;
  }

  // SOCP: add lifted variables, equality rows, then native cones.
  if (p->socp) {
    casadi_socp_data<T1>* sd = &d->socp;
    const casadi_socp_prob<T1>* sp = sd->prob;
    casadi_int b;

    casadi_socp_build(sd);

    // Lifted variable bounds (X free, Z >= 0).
    for (j = 0; j < n_lifted; ++j) {
      int jvar = p_qp->nx + j;
      T1 lo = sd->lb_lift[j], up = sd->ub_lift[j];
      int lo_inf = (lo <= -std::numeric_limits<T1>::infinity());
      int up_inf = (up >=  std::numeric_limits<T1>::infinity());
      MSKboundkeye bk;
      if (lo_inf && up_inf) { bk = MSK_BK_FR; }
      else if (lo_inf)      { bk = MSK_BK_UP; lo = -std::numeric_limits<T1>::infinity(); }
      else if (up_inf)      { bk = MSK_BK_LO; up =  std::numeric_limits<T1>::infinity(); }
      else if (lo == up)    { bk = MSK_BK_FX; }
      else                  { bk = MSK_BK_RA; }
      MSK_putvarbound(task, jvar, bk, lo, up);
    }

    // Equality rows Q*[x; lifted] = -P, placed at constraint indices na..na+n_eq-1.
    // eq_colind[] are full-variable-space indices already (0..nx+n_lifted-1).
    int eq_offset = p_qp->na;
    for (i = 0; i < (int) sp->n_eq; ++i) {
      int kbeg = sd->eq_start[i];
      int kend = sd->eq_start[i + 1];
      MSK_putarow(task, eq_offset + i, kend - kbeg,
                  sd->eq_colind + kbeg, sd->eq_coef + kbeg);
      MSK_putconbound(task, eq_offset + i, MSK_BK_FX,
                      sd->eq_rhs[i], sd->eq_rhs[i]);
    }

    // Native quadratic cones: per cone block b, the lifted layout is
    // [X_0, X_1, ..., X_{bs-2}, Z].  Mosek's MSK_CT_QUAD expects
    // [t, x_0, x_1, ...] with t the cone tip.  Permute Z to position 0.
    // Cones are model structure (not data) -- append once.
    if (!d->model_built) {
      for (b = 0; b < sp->n_blocks; ++b) {
        casadi_int bs = sp->r[b + 1] - sp->r[b];
        int base = static_cast<int>(p_qp->nx + sp->r[b]);
        d->cone_idx[0] = base + static_cast<int>(bs - 1);
        for (j = 0; j < (int) bs - 1; ++j) d->cone_idx[1 + j] = base + j;
        MSK_appendcone(task, MSK_CT_QUAD, 0.0, static_cast<int>(bs), d->cone_idx);
      }
    }
  }

  // MIP integer flags (model structure, set once)
  if (has_mip && !d->model_built) {
    for (j = 0; j < p_qp->nx; ++j) {
      d->vtype[j] = (p->coltype[j] == 'I')
                    ? MSK_VAR_TYPE_INT : MSK_VAR_TYPE_CONT;
    }
    for (j = 0; j < p_qp->nx; ++j) {
      MSK_putvartype(task, j,
          static_cast<MSKvariabletypee>(d->vtype[j]));
    }
  }

  d->model_built = 1;

  // Solve.
  MSKrescodee trm = MSK_RES_OK;
  r = MSK_optimizetrm(task, &trm);
  d->return_status = static_cast<int>(r);
  if (r != MSK_RES_OK) {
    d_qp->success = 0;
    return 1;
  }

  // Pick the best available solution slice.  After Mosek's IPM, basis
  // identification fills SOL_BAS (sharper at vertices) for LP/QP without
  // SOC; for QPs and SOCPs only SOL_ITR is meaningful; for MIP, SOL_ITG.
  MSKsoltypee soltype;
  if (has_mip) {
    soltype = MSK_SOL_ITG;
  } else {
    MSKbooleant bas_def = 0;
    MSK_solutiondef(task, MSK_SOL_BAS, &bas_def);
    soltype = bas_def ? MSK_SOL_BAS : MSK_SOL_ITR;
  }

  // Retrieve solution status
  MSKsolstae solsta = MSK_SOL_STA_UNKNOWN;
  MSKprostae prosta = MSK_PRO_STA_UNKNOWN;
  MSK_getsolsta(task, soltype, &solsta);
  MSK_getprosta(task, soltype, &prosta);
  d->sol_status  = static_cast<int>(solsta);
  d->prob_status = static_cast<int>(prosta);

  d_qp->success = (solsta == MSK_SOL_STA_OPTIMAL ||
                   solsta == MSK_SOL_STA_INTEGER_OPTIMAL ||
                   solsta == MSK_SOL_STA_PRIM_AND_DUAL_FEAS);
  if (d_qp->success) d_qp->unified_return_status = SOLVER_RET_SUCCESS;

  // Solution: x for original variables only (drop lifted SOCP slack vars)
  if (d_qp->x) {
    MSK_getxxslice(task, soltype, 0, p_qp->nx, d_qp->x);
  }

  // Objective
  MSK_getprimalobj(task, soltype, &d->obj_val);
  if (d_qp->f) *d_qp->f = d->obj_val;

  // Duals
  if (has_mip) {
    if (d_qp->lam_x) for (i = 0; i < p_qp->nx; ++i) d_qp->lam_x[i] = 0;
    if (d_qp->lam_a) for (i = 0; i < p_qp->na; ++i) d_qp->lam_a[i] = 0;
  } else {
    // Constraint duals: y = slc - suc.  CasADi sign convention: lam = -y.
    if (d_qp->lam_a && p_qp->na > 0) {
      MSK_getyslice(task, soltype, 0, p_qp->na, d_qp->lam_a);
      for (i = 0; i < p_qp->na; ++i) d_qp->lam_a[i] = -d_qp->lam_a[i];
    }
    // Variable duals: slx - sux.  Use combined snx for QP; for LP fall back to slx-sux.
    if (d_qp->lam_x && p_qp->nx > 0) {
      // Pull slx and sux separately; lam_x = -(slx - sux) per CasADi convention.
      MSK_getslxslice(task, soltype, 0, p_qp->nx, d_qp->lam_x);
      MSK_getsuxslice(task, soltype, 0, p_qp->nx, d->dual_scratch);
      for (i = 0; i < p_qp->nx; ++i)
        d_qp->lam_x[i] = -(d_qp->lam_x[i] - d->dual_scratch[i]);
    }
  }

  return 0;
}

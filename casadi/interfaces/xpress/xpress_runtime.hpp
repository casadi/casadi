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
// C-REPLACE "casadi_socp_prob<T1>" "struct casadi_socp_prob"
// C-REPLACE "casadi_socp_data<T1>" "struct casadi_socp_data"

// C-REPLACE "reinterpret_cast<int**>" "(int**) "
// C-REPLACE "reinterpret_cast<int*>" "(int*) "
// C-REPLACE "reinterpret_cast<char*>" "(char*) "
// C-REPLACE "const_cast<int*>" "(int*) "
// C-REPLACE "static_cast<int>" "(int) "

template<typename T1>
struct casadi_xpress_prob {
  const casadi_qp_prob<T1>* qp;

  // Linear constraint matrix in CSC form
  const int *colinda, *rowa;
  // Quadratic objective in upper-triangular triplet form
  // (Xpress expects only one triangle of the symmetric Q matrix;
  //  diagonal entries pass full value, off-diagonal entries pass
  //  half their actual value -- this is handled at solve time)
  const int *qobj_col1, *qobj_col2;
  // qobj_nz_idx[k] = index into the d_qp->h[] nonzero array of H_
  // (CasADi stores the full symmetric H; we picked upper-tri entries here)
  const int *qobj_nz_idx;
  int nquad;
  // Discrete variable flags ('I' for integer, 'C' for continuous), or NULL
  const char *coltype;

  // SOS constraints (empty if no SOS).  Stored in CSR-like form:
  // sos_setstart[k]..sos_setstart[k+1]-1 indexes into sos_setind/sos_refval
  // for the k-th SOS group.  sos_settype[k] is '1' or '2'.
  int n_sos_sets;
  int n_sos_elems;
  const char *sos_settype;
  const int *sos_setstart;
  const int *sos_setind;
  const double *sos_refval;

  // SOCP support.  Null if no SOC cones.
  const casadi_socp_prob<T1> *socp;
};
// C-REPLACE "casadi_xpress_prob<T1>" "struct casadi_xpress_prob"

// SYMBOL "xpress_setup"
template<typename T1>
void casadi_xpress_setup(casadi_xpress_prob<T1>* p) {

}


// SYMBOL "xpress_data"
template<typename T1>
struct casadi_xpress_data {
  // Problem structure
  const casadi_xpress_prob<T1>* prob;
  // Problem structure
  casadi_qp_data<T1>* qp;

  // Xpress problem handle
  XPRSprob xprob;

  // Workspace owned by the C runtime (sized via casadi_xpress_work)
  char *qrtype;
  int *col_idx;       // 0..nx-1, scratch for XPRSchgcoltype
  T1 *rhs;
  T1 *rng;
  T1 *qobj_val;

  // SOCP runtime data (only meaningful when prob->socp != NULL)
  casadi_socp_data<T1> socp;

  // Status / statistics
  int return_status;
  int lp_status;
  int mip_status;
  int simplex_iter;
  int barrier_iter;
  int mip_nodes;
  T1 obj_val;
};
// C-REPLACE "casadi_xpress_data<T1>" "struct casadi_xpress_data"

// SYMBOL "xpress_init_mem"
template<typename T1>
int casadi_xpress_init_mem(casadi_xpress_data<T1>* d) {
  d->xprob = 0;
  if (XPRScreateprob(&d->xprob)) return 1;
  return 0;
}

// SYMBOL "xpress_free_mem"
template<typename T1>
void casadi_xpress_free_mem(casadi_xpress_data<T1>* d) {
  if (d->xprob) {
    XPRSdestroyprob(d->xprob);
    d->xprob = 0;
  }
}

// SYMBOL "xpress_work"
template<typename T1>
void casadi_xpress_work(const casadi_xpress_prob<T1>* p,
    casadi_int* sz_arg, casadi_int* sz_res, casadi_int* sz_iw, casadi_int* sz_w) {
  casadi_qp_work(p->qp, sz_arg, sz_res, sz_iw, sz_w);
  // Workspace for row sense, rhs, range, and quadratic-objective values
  *sz_iw += p->qp->na;          // qrtype (chars packed in iw)
  *sz_iw += p->qp->nx;          // scratch int[] for XPRSchgcoltype indices
  *sz_w  += p->qp->na;          // rhs
  *sz_w  += p->qp->na;          // rng
  *sz_w  += p->nquad;           // qobj_val (scaled Q triplet values)

  if (p->socp) casadi_socp_work(p->socp, sz_iw, sz_w);
}

// SYMBOL "xpress_init"
template<typename T1>
void casadi_xpress_init(casadi_xpress_data<T1>* d,
    const T1*** arg, T1*** res, casadi_int** iw, T1** w) {
  const casadi_xpress_prob<T1>* p = d->prob;
  d->qrtype   = reinterpret_cast<char*>(*iw); *iw += p->qp->na;
  d->col_idx  = reinterpret_cast<int*>(*iw);  *iw += p->qp->nx;
  d->rhs      = *w; *w += p->qp->na;
  d->rng      = *w; *w += p->qp->na;
  d->qobj_val = *w; *w += p->nquad;

  if (p->socp) {
    d->socp.prob = p->socp;
    casadi_socp_init(&d->socp, iw, w);
  }
}


// C-REPLACE "SOLVER_RET_SUCCESS" "0"
// C-REPLACE "SOLVER_RET_LIMITED" "2"
// C-REPLACE "std::numeric_limits<T1>::infinity()" "casadi_inf"

// SYMBOL "xpress_solve"
template<typename T1>
int casadi_xpress_solve(casadi_xpress_data<T1>* d,
    const double** arg, double** res, casadi_int* iw, double* w) {

  const casadi_xpress_prob<T1>* p = d->prob;
  const casadi_qp_prob<T1>* p_qp = p->qp;
  casadi_qp_data<T1>* d_qp = d->qp;

  int i, k;
  int has_mip = p->coltype ? 1 : 0;

  // Build per-row sense/rhs/range from (lba, uba)
  for (i = 0; i < p_qp->na; ++i) {
    T1 lo = d_qp->lba[i];
    T1 up = d_qp->uba[i];
    int lo_inf = (lo <= -std::numeric_limits<T1>::infinity());
    int up_inf = (up >=  std::numeric_limits<T1>::infinity());
    if (!lo_inf && !up_inf) {
      if (lo == up) {
        d->qrtype[i] = 'E'; d->rhs[i] = up; d->rng[i] = 0;
      } else {
        d->qrtype[i] = 'R'; d->rhs[i] = up; d->rng[i] = up - lo;
      }
    } else if (!up_inf) {
      d->qrtype[i] = 'L'; d->rhs[i] = up; d->rng[i] = 0;
    } else if (!lo_inf) {
      d->qrtype[i] = 'G'; d->rhs[i] = lo; d->rng[i] = 0;
    } else {
      // Free row: model as <= +inf
      d->qrtype[i] = 'L'; d->rhs[i] = std::numeric_limits<T1>::infinity(); d->rng[i] = 0;
    }
  }

  // CasADi: f = g'x + 0.5 x'Hx with H symmetric.  Xpress: 0.5 x'Qx; pass
  // ONE triangle of Q only (Xpress auto-symmetrizes).  We pass H[i,j]
  // unmodified for both diagonal and off-diagonal entries.
  for (k = 0; k < p->nquad; ++k) {
    d->qobj_val[k] = d_qp->h[p->qobj_nz_idx[k]];
  }

  // Load the (Q)P into Xpress
  int rc = XPRSloadqp(d->xprob, "casadi_qp",
      p_qp->nx, p_qp->na,
      d->qrtype, d->rhs, d->rng,
      d_qp->g,
      p->colinda, 0, p->rowa, d_qp->a,
      d_qp->lbx, d_qp->ubx,
      p->nquad, p->qobj_col1, p->qobj_col2, d->qobj_val);
  if (rc) return 1;

  // SOCP: lift the model with helper variables, link them to the
  // original variables via Q*[x; lifted] = -P, and add a "<= 0"
  // row per cone that we make quadratic via XPRSaddqmatrix.
  if (p->socp) {
    casadi_socp_data<T1>* sd = &d->socp;
    const casadi_socp_prob<T1>* sp = sd->prob;
    casadi_int b;
    int n_rows_after_eq;

    casadi_socp_build(sd);

    if (XPRSaddcols(d->xprob, static_cast<int>(sp->n_lifted), 0,
        sd->obj_lift, sd->lift_start, 0, 0,
        sd->lb_lift, sd->ub_lift)) return 1;

    if (XPRSaddrows(d->xprob, static_cast<int>(sp->n_eq),
        static_cast<int>(sp->eq_nnz),
        sd->eq_type, sd->eq_rhs, sd->eq_rng,
        sd->eq_start, sd->eq_colind, sd->eq_coef)) return 1;

    XPRSgetintattrib(d->xprob, XPRS_ROWS, &n_rows_after_eq);

    if (XPRSaddrows(d->xprob, static_cast<int>(sp->n_blocks), 0,
        sd->cone_type, sd->cone_rhs, sd->cone_rng,
        sd->cone_start, 0, 0)) return 1;

    for (b = 0; b < sp->n_blocks; ++b) {
      casadi_int bs = casadi_socp_cone_build(sd, b);
      if (XPRSaddqmatrix(d->xprob, n_rows_after_eq + static_cast<int>(b),
          static_cast<int>(bs), sd->qcol1, sd->qcol2, sd->qcoef)) return 1;
    }
  }

  if (has_mip) {
    // Mark column types -- pass the index list and types in lockstep so
    // coltype[k] applies to colind[k] (a subtle XPRSchgcoltype API contract).
    for (i = 0; i < p_qp->nx; ++i) d->col_idx[i] = i;
    if (XPRSchgcoltype(d->xprob, p_qp->nx, d->col_idx, p->coltype)) return 1;
    if (p->n_sos_sets > 0) {
      if (XPRSaddsets(d->xprob, p->n_sos_sets, p->n_sos_elems,
                      p->sos_settype, p->sos_setstart, p->sos_setind,
                      p->sos_refval)) return 1;
    }
    rc = XPRSmipoptimize(d->xprob, "");
  } else {
    rc = XPRSlpoptimize(d->xprob, "");
  }
  if (rc) return 1;

  // Retrieve status
  if (has_mip) {
    XPRSgetintattrib(d->xprob, XPRS_MIPSTATUS, &d->mip_status);
    XPRSgetdblattrib(d->xprob, XPRS_MIPOBJVAL, &d->obj_val);
    d->return_status = d->mip_status;
    d_qp->success = (d->mip_status == XPRS_MIP_OPTIMAL ||
                     d->mip_status == XPRS_MIP_SOLUTION);
    if (d->mip_status == XPRS_MIP_OPTIMAL)
      d_qp->unified_return_status = SOLVER_RET_SUCCESS;
  } else {
    XPRSgetintattrib(d->xprob, XPRS_LPSTATUS, &d->lp_status);
    XPRSgetdblattrib(d->xprob, XPRS_LPOBJVAL, &d->obj_val);
    d->return_status = d->lp_status;
    d_qp->success = (d->lp_status == XPRS_LP_OPTIMAL);
    if (d->lp_status == XPRS_LP_OPTIMAL)
      d_qp->unified_return_status = SOLVER_RET_SUCCESS;
    if (d->lp_status == XPRS_LP_UNFINISHED)
      d_qp->unified_return_status = SOLVER_RET_LIMITED;
  }

  // Retrieve solution.  XPRSgetsolution works for both LP and MIP.
  // Duals/reduced costs are only meaningful for the LP case.
  XPRSgetsolution(d->xprob, 0, d_qp->x, 0, p_qp->nx - 1);

  if (has_mip) {
    if (d_qp->lam_x) for (i = 0; i < p_qp->nx; ++i) d_qp->lam_x[i] = 0;
    if (d_qp->lam_a) for (i = 0; i < p_qp->na; ++i) d_qp->lam_a[i] = 0;
  } else {
    // CasADi sign convention: lam = -dual returned by Xpress
    if (d_qp->lam_a && p_qp->na > 0) {
      XPRSgetduals(d->xprob, 0, d_qp->lam_a, 0, p_qp->na - 1);
      for (i = 0; i < p_qp->na; ++i) d_qp->lam_a[i] = -d_qp->lam_a[i];
    }
    if (d_qp->lam_x && p_qp->nx > 0) {
      XPRSgetredcosts(d->xprob, 0, d_qp->lam_x, 0, p_qp->nx - 1);
      for (i = 0; i < p_qp->nx; ++i) d_qp->lam_x[i] = -d_qp->lam_x[i];
    }
  }

  if (d_qp->f) *d_qp->f = d->obj_val;

  XPRSgetintattrib(d->xprob, XPRS_SIMPLEXITER, &d->simplex_iter);
  XPRSgetintattrib(d->xprob, XPRS_BARITER,     &d->barrier_iter);
  XPRSgetintattrib(d->xprob, XPRS_NODES,       &d->mip_nodes);

  return 0;
}

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

// C-REPLACE "SOLVER_RET_SUCCESS" "0"
// C-REPLACE "SOLVER_RET_INFEASIBLE" "4"

// C-REPLACE "casadi_nlpsol_prob<T1>" "struct casadi_nlpsol_prob"
// C-REPLACE "std::numeric_limits<T1>::infinity()" "casadi_inf"
// C-REPLACE "casadi_nlpsol_data<T1>" "struct casadi_nlpsol_data"
// C-REPLACE "static_cast<casadi_int>" "(casadi_int) "

// Lift: constant +-1 maps between the caller's problem and the lifted one, sizes _l lifted.
// A slack column whose rows span several stages becomes one helper state m, constant over
// the horizon, entering every lower side it relaxes with +1 and every upper side with -1;
// m_0, its stage-0 copy, is bounded 0 <= m_0 <= ubs and priced like the slack.
// SYMBOL "ipmc_rewrite_prob"
template<typename T1>
struct casadi_ipmc_rewrite_prob {
  casadi_int ne;              // helper components, 0 without lift
  casadi_int nx, ng;          // the caller's sizes
  casadi_int nnz_au, nnz_hu;  // nonzeros of the caller's Jacobian, Hessian
  const casadi_int *m0;       // [ne] stage-0 copy of each component, bounded and priced
  const casadi_int *col;      // [ne] the column of S each helper is
  const casadi_int *Px_sp, *Pg_sp, *C_sp, *B_sp, *Blo_sp, *Bup_sp, *Gz_sp, *H_sp;
  const T1 *Px;               // [nx_l x nx] caller variables, x_u = Px' x
  const T1 *Pg;               // [ng_l x ng] caller rows, g = Pg g_u + C x
  const T1 *C;                // [ng_l x nx_l] helper terms of the rows
  const T1 *B;                // [nz_l x nz] caller z entry carried, lam_u = B' lam
  const T1 *Blo, *Bup;        // B where the lower, upper bound is kept
  const T1 *Gz;               // [ng x nz_l] constraint values back
  const T1 *H;                // [nx_l x ne] helper component each variable copies
  const T1 *lbz0, *ubz0;      // [nz_l] +-inf where nothing is carried, else 0
};
// C-REPLACE "casadi_ipmc_rewrite_prob<T1>" "struct casadi_ipmc_rewrite_prob"

// SYMBOL "ipmc_rewrite_data"
template<typename T1>
struct casadi_ipmc_rewrite_data {
  casadi_nlpsol_data<T1>* user;  // the caller's problem
  T1 *m_start;                   // [ne] helper start values
};
// C-REPLACE "casadi_ipmc_rewrite_data<T1>" "struct casadi_ipmc_rewrite_data"

// SYMBOL "ipmc_rewrite_work"
template<typename T1>
void casadi_ipmc_rewrite_work(const casadi_ipmc_rewrite_prob<T1>* p, casadi_int* sz_w) {
  if (!p->ne) return;
  *sz_w += p->ne;  // m_start
}

// SYMBOL "ipmc_rewrite_set_work"
template<typename T1>
void casadi_ipmc_rewrite_set_work(casadi_ipmc_rewrite_data<T1>* d,
    const casadi_ipmc_rewrite_prob<T1>* p, T1** w) {
  if (!p->ne) return;
  d->m_start = *w; *w += p->ne;
}

// SYMBOL "ipmc_rewrite_helpers"
// Helper bounds 0 <= m_0 <= ubs, start s0 clipped to them; ipmc pins m_0 at ubs = 0,
// which hardens the rows the column relaxes
template<typename T1>
void casadi_ipmc_rewrite_helpers(const casadi_ipmc_rewrite_prob<T1>* p,
    casadi_ipmc_rewrite_data<T1>* d, casadi_nlpsol_data<T1>* v,
    const T1* ubs, const T1* s0) {
  casadi_int e, col, a;
  T1 m0;
  for (e=0;e<p->ne;++e) {
    col = p->col[e];
    a = p->m0[e];
    v->lbz[a] = 0;
    v->ubz[a] = ubs[col];
    m0 = 0;
    if (s0) {
      m0 = s0[col];
      if (!(m0>0)) m0 = 0;
      if (m0>ubs[col]) m0 = ubs[col];
    }
    d->m_start[e] = m0;
  }
  casadi_mv(p->H, p->H_sp, d->m_start, v->z, 0);
}

// SYMBOL "ipmc_rewrite_expand"
// Caller's bounds and x0 onto the lifted problem; ipmc takes no multiplier guess
template<typename T1>
void casadi_ipmc_rewrite_expand(const casadi_ipmc_rewrite_prob<T1>* p,
    casadi_ipmc_rewrite_data<T1>* d, casadi_nlpsol_data<T1>* v,
    const T1* ubs, const T1* s0) {
  casadi_int nzl;
  casadi_nlpsol_data<T1>* u = d->user;
  if (!p->ne) return;
  nzl = v->prob->nx + v->prob->ng;
  v->p = u->p;
  casadi_copy(p->lbz0, nzl, v->lbz);
  casadi_mv(p->Blo, p->Blo_sp, u->lbz, v->lbz, 0);
  casadi_copy(p->ubz0, nzl, v->ubz);
  casadi_mv(p->Bup, p->Bup_sp, u->ubz, v->ubz, 0);
  casadi_clear(v->z, v->prob->nx);
  casadi_mv(p->Px, p->Px_sp, u->z, v->z, 0);
  casadi_clear(v->lam, nzl);
  casadi_ipmc_rewrite_helpers(p, d, v, ubs, s0);
}

// SYMBOL "ipmc_rewrite_collect"
// Lifted solution back to the caller: x, lam = B' lam, and s, lam_s from m_0
template<typename T1>
void casadi_ipmc_rewrite_collect(const casadi_ipmc_rewrite_prob<T1>* p,
    casadi_ipmc_rewrite_data<T1>* d, casadi_nlpsol_data<T1>* v,
    T1* s, T1* lam_s) {
  casadi_int e;
  casadi_nlpsol_data<T1>* u = d->user;
  if (!p->ne) return;
  casadi_clear(u->z, p->nx);
  casadi_mv(p->Px, p->Px_sp, v->z, u->z, 1);
  casadi_clear(u->lam, p->nx+p->ng);
  casadi_mv(p->B, p->B_sp, v->lam, u->lam, 1);
  for (e=0;e<p->ne;++e) {
    s[p->col[e]] = v->z[p->m0[e]];
    lam_s[p->col[e]] = v->lam[p->m0[e]];
  }
}

// SYMBOL "ipmc_rewrite_collect_g"
// Lifted constraint values of an intermediate iterate back to the caller's g
template<typename T1>
void casadi_ipmc_rewrite_collect_g(const casadi_ipmc_rewrite_prob<T1>* p,
    casadi_ipmc_rewrite_data<T1>* d, casadi_nlpsol_data<T1>* v) {
  casadi_nlpsol_data<T1>* u = d->user;
  if (!p->ne) return;
  casadi_clear(u->z+p->nx, p->ng);
  casadi_mv(p->Gz, p->Gz_sp, v->z, u->z+p->nx, 0);
}

// C-REPLACE "casadi_ocp_block" "struct casadi_ocp_block"

// The handed-over problem as ipmc sees it, fixed when the solver is built
// SYMBOL "ipmc_prob"
template<typename T1>
struct casadi_ipmc_prob {
  const casadi_nlpsol_prob<T1>* nlp;
  // Stage partition of the handed-over problem, k = 0..N
  const casadi_int *nx, *nu;
  casadi_int N;
  // Jacobian blocks over [x_k; u_k]: AB[k] dynamics rows, CD[k] path rows
  casadi_ocp_block *AB, *CD;

  // Slacks: ipmc's n_soft stage-local slack variables, stage after stage
  casadi_int n_soft;
  const casadi_int *slack_perm;  // [n_soft] the column of S each ipmc slack is
  const T1 *fs_z;                // [ns] penalty gradient at s=0, per column
  const T1 *fs_Z;                // [ns] penalty Hessian diagonal, per column

  // Lift (rewrite.ne==0: none); everything above describes the lifted problem
  casadi_ipmc_rewrite_prob<T1> rewrite;

  // ipmc's inequality rows as z-space indices, stage after stage
  casadi_int n_ineq;
  const casadi_int *ineq_z;
  // [n_ineq] ipmc slack relaxing the lower / upper side of the row, -1 if hard
  const casadi_int *ineq_lo, *ineq_up;

  // Pack tables: code c>0 reads src[c-1], c<0 -src[-c-1]; column c to column col[c] of blk[c]
  const casadi_int *bat_sp, *bat_code, *bat_blk, *bat_col;  // BAt
  const casadi_int *rsq_sp, *rsq_code, *rsq_blk, *rsq_col;  // RSQ, full diagonal
  const casadi_int *rsqs_sp, *rsqs_code;                    // column k: RSQ_slack[k]
  const casadi_int *gi_sp, *gi_code, *gi_blk, *gi_col;      // Gt_ineq
  // Jacobian checks: au[ichk[i]] == 1; cchk (nonzero, value, stage, state) per constant state
  casadi_int n_ichk, n_cchk;
  const casadi_int *ichk, *cchk;
};
// C-REPLACE "casadi_ipmc_prob<T1>" "struct casadi_ipmc_prob"

// ipmc's stage vector is [u_k; x_k], casadi's [x_k; u_k] at CD[k].offset_c

// SYMBOL "ipmc_read_primal_data"
template<typename T1>
void casadi_ipmc_read_primal_data(const casadi_ipmc_prob<T1>* p,
    const double* primal_data, T1* x, const struct IpmcLayout *s) {
  casadi_int k;
  for (k=0;k<s->K;++k) {
    casadi_copy(primal_data+s->ux_offs[k], s->nu[k], x+p->CD[k].offset_c+p->nx[k]);
    casadi_copy(primal_data+s->ux_offs[k]+s->nu[k], p->nx[k], x+p->CD[k].offset_c);
  }
}

// SYMBOL "ipmc_write_primal_data"
template<typename T1>
void casadi_ipmc_write_primal_data(const casadi_ipmc_prob<T1>* p,
    const double* x, T1* primal_data, const struct IpmcLayout *s) {
  casadi_int k;
  for (k=0;k<s->K;++k) {
    casadi_copy(x+p->CD[k].offset_c+p->nx[k], s->nu[k], primal_data+s->ux_offs[k]);
    casadi_copy(x+p->CD[k].offset_c, p->nx[k], primal_data+s->ux_offs[k]+s->nu[k]);
  }
}

// SYMBOL "ipmc_error_t"
typedef enum {
  CASADI_IPMC_OK,
  CASADI_IPMC_GAP_BOUNDS,    // dynamics row (index: z) with bounds not equal and finite
  CASADI_IPMC_UBS,           // ubs <= 0 on a stage-local slack relaxing a finite bound
  CASADI_IPMC_PENALTY,       // ipmc refused the slack penalty (index: the IpmcError)
  CASADI_IPMC_IDENTITY,      // gap-closing entry au[ichk[index]] is not 1
  CASADI_IPMC_NXC            // constant state cchk[4*index..] has another dynamics row
} casadi_ipmc_error_t;

// SYMBOL "ipmc_data"
template<typename T1>
struct casadi_ipmc_data {
  const casadi_ipmc_prob<T1>* prob;
  casadi_nlpsol_data<T1>* nlp;  // the handed-over problem (lifted if any)
  casadi_ipmc_rewrite_data<T1> rewrite;

  const T1** arg;
  T1** res;
  casadi_int* iw;
  T1* w;

  int unified_return_status;
  int success;
  int return_status;
  casadi_ipmc_error_t error;
  casadi_int error_index;

  // Handed-over problem: point, constraint values or gradient, multipliers
  T1 *x, *g, *lam;
  // Caller oracle's buffers, aliases without lift; au, hu end in the constants the tables read
  T1 *xu, *gu, *au, *hu, *lamu;

  // Native slacks: NlpsolMemory::slack_* ([ns], never null if ns>0)
  T1 *slack_s, *slack_lam_s, *slack_ubs;
  // ipmc's slacks and their bound multipliers, [n_soft] each
  T1 *sv, *zs_lo, *zs_up;
  int nxc_checked;    // constant-state rows checked this solve

  // Zero the blocks at the first evaluation of a solve; ipmc modifies them in place
  int jac_fresh, hess_fresh;

  // Scratch for ipmc_set_bounds/_initial/_soft_penalty/_initial_slack
  T1 *set_lower, *set_upper, *set_ux0, *set_pen, *set_s0;

  // g for the iteration callback, and ipmc's point it was evaluated at
  T1 *cb_g, *cb_x;

  struct IpmcStats stats;

  // Built once per memory object in memory casadi owns
  struct IpmcSolver *solver;
  const struct IpmcLayout *layout;
};
// C-REPLACE "casadi_ipmc_data<T1>" "struct casadi_ipmc_data"

// SYMBOL "ipmc_init_mem"
// The memory object's own state; init_mem builds the solver next to it
template<typename T1>
int casadi_ipmc_init_mem(casadi_ipmc_data<T1>* d) {
  d->solver = 0;
  d->layout = 0;
  d->error = CASADI_IPMC_OK;
  d->nxc_checked = 0;
  return 0;
}

// Lifting the caller oracle's outputs; no-ops without lift

// SYMBOL "ipmc_penalty"
// *acc += sum_e (z_e m_e + 1/2 Z_e m_e^2), m_e the stage-0 helpers in x
template<typename T1>
void casadi_ipmc_penalty(const casadi_ipmc_prob<T1>* p, const T1* x, T1* acc) {
  casadi_int e, col;
  T1 m0;
  const casadi_ipmc_rewrite_prob<T1>* q = &p->rewrite;
  for (e=0;e<q->ne;++e) {
    col = q->col[e];
    m0 = x[q->m0[e]];
    *acc += p->fs_z[col]*m0 + 0.5*p->fs_Z[col]*m0*m0;
  }
}

// SYMBOL "ipmc_penalty_grad"
// g += scale * gradient of the penalty
template<typename T1>
void casadi_ipmc_penalty_grad(const casadi_ipmc_prob<T1>* p, const T1* x, T1 scale, T1* g) {
  casadi_int e, col;
  const casadi_ipmc_rewrite_prob<T1>* q = &p->rewrite;
  for (e=0;e<q->ne;++e) {
    col = q->col[e];
    g[q->m0[e]] += scale*(p->fs_z[col] + p->fs_Z[col]*x[q->m0[e]]);
  }
}

// SYMBOL "ipmc_lift_x"
// xu = Px' x
template<typename T1>
void casadi_ipmc_lift_x(const casadi_ipmc_prob<T1>* p, casadi_ipmc_data<T1>* d) {
  const casadi_ipmc_rewrite_prob<T1>* q = &p->rewrite;
  if (!q->ne) return;
  casadi_clear(d->xu, q->nx);
  casadi_mv(q->Px, q->Px_sp, d->x, d->xu, 1);
}

// SYMBOL "ipmc_set_x"
// ipmc's point x into d->x and the caller oracle's first two inputs
template<typename T1>
void casadi_ipmc_set_x(casadi_ipmc_data<T1>* d, const double* x) {
  casadi_ipmc_read_primal_data(d->prob, x, d->x, d->layout);
  casadi_ipmc_lift_x(d->prob, d);
  d->arg[0] = d->xu;
  d->arg[1] = d->nlp->p;
}

// SYMBOL "ipmc_lift_obj"
// f += penalty
template<typename T1>
void casadi_ipmc_lift_obj(const casadi_ipmc_prob<T1>* p, casadi_ipmc_data<T1>* d, T1* obj) {
  casadi_ipmc_penalty(p, d->x, obj);
}

// SYMBOL "ipmc_lift_grad_f"
// g = Px gu + grad penalty
template<typename T1>
void casadi_ipmc_lift_grad_f(const casadi_ipmc_prob<T1>* p, casadi_ipmc_data<T1>* d) {
  const casadi_ipmc_rewrite_prob<T1>* q = &p->rewrite;
  if (!q->ne) return;
  casadi_clear(d->g, p->nlp->nx);
  casadi_mv(q->Px, q->Px_sp, d->gu, d->g, 0);
  casadi_ipmc_penalty_grad(p, d->x, 1., d->g);
}

// SYMBOL "ipmc_lift_g"
// g = Pg gu + C x
template<typename T1>
void casadi_ipmc_lift_g(const casadi_ipmc_prob<T1>* p, casadi_ipmc_data<T1>* d) {
  const casadi_ipmc_rewrite_prob<T1>* q = &p->rewrite;
  if (!q->ne) return;
  casadi_clear(d->g, p->nlp->ng);
  casadi_mv(q->Pg, q->Pg_sp, d->gu, d->g, 0);
  casadi_mv(q->C, q->C_sp, d->x, d->g, 0);
}

// SYMBOL "ipmc_lift_lam"
// lamu = Pg' lam
template<typename T1>
void casadi_ipmc_lift_lam(const casadi_ipmc_prob<T1>* p, casadi_ipmc_data<T1>* d) {
  const casadi_ipmc_rewrite_prob<T1>* q = &p->rewrite;
  if (!q->ne) return;
  casadi_clear(d->lamu, q->ng);
  casadi_mv(q->Pg, q->Pg_sp, d->lam, d->lamu, 1);
}

// SYMBOL "ipmc_lift_hess_l"
// Lagrangian gradient g = Px gu + C' lam + obj_scale grad penalty; the Hessian needs none
template<typename T1>
void casadi_ipmc_lift_hess_l(const casadi_ipmc_prob<T1>* p, casadi_ipmc_data<T1>* d,
    T1 obj_scale) {
  const casadi_ipmc_rewrite_prob<T1>* q = &p->rewrite;
  if (!q->ne) return;
  casadi_clear(d->g, p->nlp->nx);
  casadi_mv(q->Px, q->Px_sp, d->gu, d->g, 0);
  casadi_mv(q->C, q->C_sp, d->lam, d->g, 1);
  casadi_ipmc_penalty_grad(p, d->x, obj_scale, d->g);
}

// SYMBOL "ipmc_zval"
// [d->x; d->g][z]
template<typename T1>
T1 casadi_ipmc_zval(const casadi_ipmc_data<T1>* d, casadi_int z) {
  casadi_int nx = d->prob->nlp->nx;
  return z<nx ? d->x[z] : d->g[z-nx];
}

// SYMBOL "ipmc_snapshot_g"
// Copy d->g (shared scratch) and ipmc's point for the iteration callback
template<typename T1>
void casadi_ipmc_snapshot_g(casadi_ipmc_data<T1>* d, const struct IpmcLayout* s,
            const double* primal_data) {
  const casadi_ipmc_prob<T1>* p = d->prob;
  casadi_int n_ux = s->ux_offs[s->K-1] + s->nu[s->K-1] + s->nx[s->K-1];

  casadi_copy(d->g, p->nlp->ng, d->cb_g);
  casadi_copy(primal_data, n_ux, d->cb_x);
}

// SYMBOL "ipmc_stage_values"
// Stage k: gin = z on inequalities, gap = lbg - g on dynamics
template<typename T1>
void casadi_ipmc_stage_values(const casadi_ipmc_data<T1>* d, casadi_int k,
    T1* gin, T1* gap) {
  casadi_int i, z;
  const casadi_ipmc_prob<T1>* p = d->prob;
  const struct IpmcLayout* s = d->layout;
  const T1* lbz = d->nlp->lbz;
  for (i=0;i<s->ng_ineq[k];++i) {
    gin[i] = casadi_ipmc_zval(d, p->ineq_z[s->ineq_offs[k]+i]);
  }
  if (k==p->N) return;
  for (i=0;i<p->nx[k+1];++i) {
    z = p->nlp->nx+p->AB[k].offset_r+i;
    gap[i] = lbz[z]-casadi_ipmc_zval(d, z);
  }
}

// SYMBOL "ipmc_pack_constr_viol"
// Constraint values (already in d->g) into ipmc's layout
template<typename T1>
void casadi_ipmc_pack_constr_viol(casadi_ipmc_data<T1>* d, const struct IpmcLayout* s,
            double* res) {
  casadi_int k;
  for (k=0;k<s->K;++k) {
    casadi_ipmc_stage_values(d, k, res+s->g_ineq_offs[k], k<s->K-1 ? res+s->dyn_eq_offs[k] : 0);
  }
}


// SYMBOL "ipmc_read_lam"
// ipmc's row multipliers into d->lam; dynamics flip sign, ipmc's rows are x_{k+1} - F
template<typename T1>
void casadi_ipmc_read_lam(casadi_ipmc_data<T1>* d, const struct IpmcLayout* s,
            const double* lam_data) {
  casadi_int k, i, z, nx;
  const casadi_ipmc_prob<T1>* p = d->prob;

  // Rows only, every one an ipmc row; simple bounds enter in casadi_ipmc_pack_lag_hess
  nx = p->nlp->nx;
  for (k=0;k<s->K;++k) {
    for (i=0;i<s->ng_ineq[k];++i) {
      z = p->ineq_z[s->ineq_offs[k]+i];
      if (z>=nx) d->lam[z-nx] = lam_data[s->g_ineq_offs[k]+i];
    }
  }
  for (k=0;k<s->K-1;++k) {
    casadi_scaled_copy(-1.0, lam_data+s->dyn_eq_offs[k], p->nx[k+1], d->lam+p->AB[k].offset_r);
  }
}

// SYMBOL "ipmc_scatter"
// Execute a pack table (see casadi_ipmc_prob)
template<typename T1>
void casadi_ipmc_scatter(const casadi_int* sp, const casadi_int* code,
    const casadi_int* blk, const casadi_int* col, const T1* src, struct blasfeo_dmat* M) {
  casadi_int ncol, c, el, s;
  const casadi_int *colind, *row;
  struct blasfeo_dmat* B;
  ncol = sp[1];
  colind = sp+2;
  row = colind+ncol+1;
  for (c=0;c<ncol;++c) {
    B = M+blk[c];
    for (el=colind[c];el<colind[c+1];++el) {
      s = code[el];
      BLASFEO_DMATEL(B, row[el], col[c]) = s>0 ? src[s-1] : -src[-s-1];
    }
  }
}

// SYMBOL "ipmc_pack_lag_hess"
// Hessian and Lagrangian gradient into RSQ, RSQ_slack, rq; rq adds x_{k+1} and bound terms
template<typename T1>
void casadi_ipmc_pack_lag_hess(casadi_ipmc_data<T1>* d, const struct IpmcLayout* s,
            const double* lam_data, T1 obj_scale, struct blasfeo_dmat* RSQ_p,
            struct blasfeo_dvec* RSQ_slack_p, struct blasfeo_dvec* rq_p) {
  casadi_int k, i, e, c, z;
  const casadi_int *colind, *row;
  const casadi_ipmc_prob<T1>* p = d->prob;
  T1* src = d->hu;

  // x_{k+1} term of the dynamics rows
  for (k=0;k<s->K-1;++k) {
    casadi_axpy(p->nx[k+1], 1.0, lam_data+s->dyn_eq_offs[k], d->g+p->CD[k+1].offset_c);
  }
  // Simple-bound multipliers
  for (k=0;k<s->K;++k) {
    for (i=0;i<s->ng_ineq[k];++i) {
      z = p->ineq_z[s->ineq_offs[k]+i];
      if (z<p->nlp->nx) d->g[z] += lam_data[s->g_ineq_offs[k]+i];
    }
  }
  // rq is [u_k; x_k]
  for (k=0;k<s->K;++k) {
    blasfeo_pack_dvec(p->nx[k], d->g+p->CD[k].offset_c, 1, rq_p+k, p->nu[k]);
    blasfeo_pack_dvec(p->nu[k], d->g+p->CD[k].offset_c+p->nx[k], 1, rq_p+k, 0);
  }

  // Constants after the oracle's nonzeros: 0, then the penalty curvature of each helper
  src[p->rewrite.nnz_hu] = 0;
  for (e=0;e<p->rewrite.ne;++e) {
    src[p->rewrite.nnz_hu+1+e] = obj_scale*p->fs_Z[p->rewrite.col[e]];
  }
  if (d->hess_fresh) {
    for (k=0;k<s->K;++k) blasfeo_dgese(RSQ_p[k].m, RSQ_p[k].n, 0.0, RSQ_p+k, 0, 0);
    d->hess_fresh = 0;
  }
  casadi_ipmc_scatter(p->rsq_sp, p->rsq_code, p->rsq_blk, p->rsq_col, src, RSQ_p);
  colind = p->rsqs_sp+2;
  row = colind+p->rsqs_sp[1]+1;
  for (k=0;k<p->rsqs_sp[1];++k) {
    for (i=colind[k];i<colind[k+1];++i) {
      c = p->rsqs_code[i];
      BLASFEO_DVECEL(RSQ_slack_p+k, row[i]) = c>0 ? src[c-1] : -src[-c-1];
    }
  }
}

// SYMBOL "ipmc_pack_constr_jac"
// Jacobian (d->au) into BAt, Gt_ineq and b, g_ineq; 1 with d->error on a broken promise
template<typename T1>
int casadi_ipmc_pack_constr_jac(casadi_ipmc_data<T1>* d, const struct IpmcLayout* s,
            struct blasfeo_dmat* BAt_p, struct blasfeo_dvec* b_p,
            struct blasfeo_dmat* Gt_ineq_p, struct blasfeo_dvec* g_ineq_p) {
  casadi_int i, k;
  const casadi_ipmc_prob<T1>* p = d->prob;

  d->au[p->rewrite.nnz_au] = 1.0;
  if (d->jac_fresh) {
    for (k=0;k<s->K;++k) {
      if (k<s->K-1) blasfeo_dgese(BAt_p[k].m, BAt_p[k].n, 0.0, BAt_p+k, 0, 0);
      blasfeo_dgese(Gt_ineq_p[k].m, Gt_ineq_p[k].n, 0.0, Gt_ineq_p+k, 0, 0);
    }
    d->jac_fresh = 0;
  }
  casadi_ipmc_scatter(p->bat_sp, p->bat_code, p->bat_blk, p->bat_col, d->au, BAt_p);
  casadi_ipmc_scatter(p->gi_sp, p->gi_code, p->gi_blk, p->gi_col, d->au, Gt_ineq_p);

  for (i=0;i<p->n_ichk;++i) {
    if (d->au[p->ichk[i]]!=1.0) {
      d->error = CASADI_IPMC_IDENTITY;
      d->error_index = i;
      return 1;
    }
  }
  if (!d->nxc_checked) {
    for (i=0;i<p->n_cchk;++i) {
      if (d->au[p->cchk[4*i]]!=p->cchk[4*i+1]) {
        d->error = CASADI_IPMC_NXC;
        d->error_index = i;
        return 1;
      }
    }
    d->nxc_checked = 1;
  }

  for (k=0;k<s->K;++k) {
    casadi_ipmc_stage_values(d, k, g_ineq_p[k].pa, k<s->K-1 ? b_p[k].pa : 0);
  }
  return 0;
}

// SYMBOL "ipmc_work"
template<typename T1>
void casadi_ipmc_work(const casadi_ipmc_prob<T1>* p, casadi_int* sz_arg, casadi_int* sz_res,
    casadi_int* sz_iw, casadi_int* sz_w) {
  casadi_nlpsol_work(p->nlp, sz_arg, sz_res, sz_iw, sz_w);

  *sz_w += p->nlp->nx;                          // x
  *sz_w += p->nlp->nx+p->nlp->ng;               // lam
  *sz_w += casadi_max(p->nlp->nx, p->nlp->ng);  // g
  *sz_w += p->rewrite.nnz_au + 1;                   // au, then 1
  *sz_w += p->rewrite.nnz_hu + 1 + p->rewrite.ne;   // hu, then 0 and the helper curvatures

  *sz_w += 2*p->n_ineq;          // set_lower, set_upper
  *sz_w += p->nlp->nx;           // set_ux0
  *sz_w += 4*p->n_soft;          // set_pen (Z, z, ubs), set_s0
  *sz_w += p->nlp->ng;           // cb_g
  *sz_w += p->nlp->nx;           // cb_x, ipmc's primal length
  *sz_w += 3*p->n_soft;          // sv, zs_lo, zs_up

  // Caller oracle buffers, lift only
  if (p->rewrite.ne) {
    *sz_w += p->rewrite.nx;                             // xu
    *sz_w += casadi_max(p->rewrite.nx, p->rewrite.ng);  // gu
    *sz_w += p->rewrite.ng;                             // lamu
  }
  casadi_ipmc_rewrite_work(&p->rewrite, sz_w);
}

// SYMBOL "ipmc_set_work"
template<typename T1>
void casadi_ipmc_set_work(casadi_ipmc_data<T1>* d, const T1*** arg, T1*** res,
    casadi_int** iw, T1** w) {
  const casadi_ipmc_prob<T1>* p = d->prob;

  d->x = *w; *w += p->nlp->nx;
  d->lam = *w; *w += p->nlp->nx+p->nlp->ng;
  d->g = *w; *w += casadi_max(p->nlp->nx, p->nlp->ng);
  d->au = *w; *w += p->rewrite.nnz_au + 1;
  d->hu = *w; *w += p->rewrite.nnz_hu + 1 + p->rewrite.ne;

  d->set_lower = *w; *w += p->n_ineq;
  d->set_upper = *w; *w += p->n_ineq;
  d->set_ux0 = *w;   *w += p->nlp->nx;
  d->set_pen = *w;   *w += 3*p->n_soft;
  d->set_s0 = *w;    *w += p->n_soft;
  d->cb_g = *w;      *w += p->nlp->ng;
  d->cb_x = *w;      *w += p->nlp->nx;
  d->sv = *w;        *w += p->n_soft;
  d->zs_lo = *w;     *w += p->n_soft;
  d->zs_up = *w;     *w += p->n_soft;

  // Caller oracle buffers; aliases without lift
  if (p->rewrite.ne) {
    d->xu = *w;   *w += p->rewrite.nx;
    d->gu = *w;   *w += casadi_max(p->rewrite.nx, p->rewrite.ng);
    d->lamu = *w; *w += p->rewrite.ng;
  } else {
    d->xu = d->x;
    d->gu = d->g;
    d->lamu = d->lam;
  }
  casadi_ipmc_rewrite_set_work(&d->rewrite, &p->rewrite, w);

  d->arg = *arg;
  d->res = *res;
  d->iw = *iw;
  d->w = *w;
}

// SYMBOL "ipmc_hand_over"
// Bounds, x0, slack penalty, ubs and s0 onto the solver; 1 with d->error if they do not fit
template<typename T1>
int casadi_ipmc_hand_over(casadi_ipmc_data<T1>* d) {
  casadi_int k, i, z, col;
  T1 lo, up;
  IpmcError err;
  const casadi_ipmc_prob<T1>* p = d->prob;
  casadi_nlpsol_data<T1>* d_nlp = d->nlp;
  const T1 inf = std::numeric_limits<T1>::infinity();
  T1 *Z, *zl, *ubs;
  d->error = CASADI_IPMC_OK;
  d->nxc_checked = 0;
  d->jac_fresh = 1;
  d->hess_fresh = 1;
  // d->sv: the stage-local slack relaxes a finite bound; ubs <= 0 elsewhere is inert
  casadi_clear(d->sv, p->n_soft);
  for (i=0;i<p->n_ineq;++i) {
    if (p->ineq_lo[i]>=0 && d_nlp->lbz[p->ineq_z[i]]>-inf) d->sv[p->ineq_lo[i]] = 1;
    if (p->ineq_up[i]>=0 && d_nlp->ubz[p->ineq_z[i]]<inf) d->sv[p->ineq_up[i]] = 1;
  }
  for (i=0;i<p->n_soft;++i) {
    col = p->slack_perm[i];
    if (!(d->slack_ubs[col]>0) && d->sv[i]) {
      d->error = CASADI_IPMC_UBS;
      d->error_index = col;
      return 1;
    }
  }
  for (k=0;k<p->N;++k) {
    for (i=0;i<p->nx[k+1];++i) {
      z = p->nlp->nx+p->AB[k].offset_r+i;
      lo = d_nlp->lbz[z];
      up = d_nlp->ubz[z];
      if (!(lo==up && lo>-inf && lo<inf)) {
        d->error = CASADI_IPMC_GAP_BOUNDS;
        d->error_index = z;
        return 1;
      }
    }
  }
  for (i=0;i<p->n_ineq;++i) {
    d->set_lower[i] = d_nlp->lbz[p->ineq_z[i]];
    d->set_upper[i] = d_nlp->ubz[p->ineq_z[i]];
  }
  ipmc_set_bounds(d->solver, d->set_lower, d->set_upper);
  casadi_ipmc_write_primal_data(p, d_nlp->z, d->set_ux0, d->layout);
  ipmc_set_initial(d->solver, d->set_ux0);
  if (p->n_soft==0) return 0;
  Z = d->set_pen;
  zl = Z + p->n_soft;
  ubs = zl + p->n_soft;
  for (i=0;i<p->n_soft;++i) {
    col = p->slack_perm[i];
    Z[i] = p->fs_Z[col];
    zl[i] = p->fs_z[col];
    ubs[i] = d->slack_ubs[col]>0 ? d->slack_ubs[col] : 1;
    d->set_s0[i] = d->slack_s[col];
  }
  err = ipmc_set_soft_penalty(d->solver, Z, zl, ubs);
  if (err!=IPMC_OK) {
    d->error = CASADI_IPMC_PENALTY;
    d->error_index = static_cast<casadi_int>(err);
    return 1;
  }
  ipmc_set_initial_slack(d->solver, d->set_s0);
  return 0;
}

// SYMBOL "ipmc_read_dual"
// ipmc's multipliers into d->nlp->lam; every z entry is an ipmc row
template<typename T1>
void casadi_ipmc_read_dual(casadi_ipmc_data<T1>* d, const struct IpmcLayout* str,
                              const double* dual_data) {
  casadi_int k, i;
  const casadi_ipmc_prob<T1>* p = d->prob;
  const casadi_nlpsol_prob<T1>* p_nlp = p->nlp;
  casadi_nlpsol_data<T1>* d_nlp = d->nlp;
  for (k=0;k<str->K;++k) {
    for (i=0;i<str->ng_ineq[k];++i) {
      d_nlp->lam[p->ineq_z[str->ineq_offs[k]+i]] = dual_data[str->g_ineq_offs[k]+i];
    }
  }
  // Dynamics
  for (k=0;k<str->K-1;++k) {
    casadi_scaled_copy(-1.0, dual_data+str->dyn_eq_offs[k], p->nx[k+1],
                       d_nlp->lam+p_nlp->nx+p->AB[k].offset_r);
  }
}

// SYMBOL "ipmc_read_slacks"
// ipmc's slacks into the columns of S; lam_s = zs_up - zs_lo, negative at s = 0
template<typename T1>
void casadi_ipmc_read_slacks(casadi_ipmc_data<T1>* d) {
  casadi_int i, col;
  const casadi_ipmc_prob<T1>* p = d->prob;
  if (p->n_soft==0) return;
  ipmc_get_slack(d->solver, d->sv);
  ipmc_get_slack_dual(d->solver, d->zs_lo, d->zs_up);
  for (i=0;i<p->n_soft;++i) {
    col = p->slack_perm[i];
    d->slack_s[col] = d->sv[i];
    d->slack_lam_s[col] = d->zs_up[i] - d->zs_lo[i];
    // An inert slack pinned by ubs <= 0 is zero
    if (!(d->slack_ubs[col]>0)) d->slack_s[col] = d->slack_lam_s[col] = 0;
  }
}

// SYMBOL "ipmc_report_iterate"
// An iterate (x, g, lam, s, f) into the caller's z-space for Nlpsol::callback
template<typename T1>
void casadi_ipmc_report_iterate(casadi_ipmc_data<T1>* d, const struct IpmcLayout* str,
    const double* primal_data, const double* dual_data, T1 f) {
  T1 pen;
  const casadi_ipmc_prob<T1>* p = d->prob;
  const casadi_nlpsol_prob<T1>* p_nlp = p->nlp;
  casadi_nlpsol_data<T1>* d_nlp = d->nlp;
  casadi_ipmc_read_primal_data(p, primal_data, d_nlp->z, str);
  casadi_copy(d->cb_g, p_nlp->ng, d_nlp->z+p_nlp->nx);
  if (dual_data) casadi_ipmc_read_dual(d, str, dual_data);
  casadi_ipmc_read_slacks(d);
  pen = 0;
  casadi_ipmc_penalty(p, d_nlp->z, &pen);
  d->rewrite.user->objective = f - pen;
  casadi_ipmc_rewrite_collect(&p->rewrite, &d->rewrite, d_nlp,
                                d->slack_s, d->slack_lam_s);
  casadi_ipmc_rewrite_collect_g(&p->rewrite, &d->rewrite, d_nlp);
}

// SYMBOL "ipmc_finish"
// After the solve loop: status, stats, solution and slacks
template<typename T1>
void casadi_ipmc_finish(casadi_ipmc_data<T1>* d) {
  const casadi_ipmc_prob<T1>* p = d->prob;
  casadi_nlpsol_data<T1>* d_nlp = d->nlp;
  const struct IpmcLayout* str = d->layout;
  const double *primal_data, *dual_data;
  ipmc_int ret;

  ret = ipmc_get_status(d->solver);

  d->return_status = ret;
  if (ret==IPMC_SOLVED) {
    d->unified_return_status = SOLVER_RET_SUCCESS;
    d->success = 1;
  }

  if (ret==IPMC_INFEASIBLE) {
    // Restoration converged to an infeasible point
    d->unified_return_status = SOLVER_RET_INFEASIBLE;
  }

  primal_data = ipmc_get_primal(d->solver);
  dual_data = ipmc_get_dual(d->solver);

  d->stats = *ipmc_get_stats(d->solver);

  casadi_ipmc_read_primal_data(p, primal_data, d_nlp->z, str);
  casadi_ipmc_read_dual(d, str, dual_data);
  casadi_ipmc_read_slacks(d);
}

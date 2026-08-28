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
// C-REPLACE "SOLVER_RET_UNKNOWN" "1"
// C-REPLACE "SOLVER_RET_LIMITED" "2"
// C-REPLACE "SOLVER_RET_INFEASIBLE" "4"
// C-REPLACE "SOLVER_RET_EXCEPTION" "5"

// C-REPLACE "casadi_nlpsol_prob<T1>" "struct casadi_nlpsol_prob"
// C-REPLACE "std::numeric_limits<T1>::infinity()" "casadi_inf"
// C-REPLACE "casadi_nlpsol_data<T1>" "struct casadi_nlpsol_data"

// C-REPLACE "reinterpret_cast<int**>" "(int**) "
// C-REPLACE "reinterpret_cast<int*>" "(int*) "
// C-REPLACE "const_cast<int*>" "(int*) "

/* ==== problem rewrite: the caller's z-space onto a lifted problem ========

   Nlpsol owns the caller's decision-vector space, so a plugin that
   reformulates the problem into a larger equivalent one (more variables,
   more rows) carries a second casadi_nlpsol_prob/_data and moves the
   caller's runtime inputs (bounds, x0, multiplier guesses) onto it on every
   solve.  This block is that move: a pure vector scatter driven by the index
   maps below, which are construction-time constants; only the loops run per
   solve.  One implementation serves both the C++ and the generated-C path,
   like casadi_nlpsol_detect_bounds_before/_after in casadi_nlp.hpp.

   The rewrite is "slack lifting": a slack column that is not local to one
   stage becomes ne HELPER VARIABLES, replicated over K stages and tied
   together by equality rows, with every softened row turned into a hard
   one-sided row against its helper.  The plugin chooses the layout; the
   semantics below are fixed.

     ne        : helper components; 0 means nothing was rewritten and none
                 of this is read
     nx, ng    : the CALLER's nx and ng
     K         : stages the helper chain is replicated over
     x[i]      : caller variable i -> rewritten variable
     xfree[i]  : its simple bound moved to a path row (softened)
     row[r]    : rewritten row r -> the caller z-space index it came from
                 (variable i of a softened bound, or nx + constraint row);
                 -1 for a continuity row.  Pre-offset at construction so one
                 loop body serves a split row and a split bound.
     kind[r]   : 0 plain row, 1 lower / 2 upper half of a split row,
                 3 helper continuity row
     ent[r]    : the helper component r relaxes, -1 if none
     mflat[r]  : flat helper index (stage*ne + component) of the helper the
                 row touches, -1 if none.  m[mflat[r]] is the helper VARIABLE
                 at r's own stage; a continuity row ties it to the next
                 stage's copy m[mflat[r]+ne].
     m[k*ne+e] : helper component e of stage k -> rewritten variable
     slot[e]   : helper component e -> slack slot (side*ns + column)

   The ORACLE stays the caller's: the rewritten problem's f, g, gradients,
   jacobian and Hessian are produced numerically at solve time by the lift_*
   routines further down, out of the caller oracle's outputs.  Two more
   construction-time maps drive the sparse ones (nnz_au / nnz_hu are the
   caller's nonzero counts):
     a_src[i]  : rewritten jacobian nonzero i -> caller nonzero, or the
                 constant -1 = +1.0 / -2 = -1.0 (a helper, split-bound or
                 continuity coefficient)
     h_src[i]  : rewritten hessian nonzero i -> caller nonzero, or -2-e for
                 the stage-0 penalty diagonal of helper component e */

// SYMBOL "ipmc_rewrite_prob"
template<typename T1>
struct casadi_ipmc_rewrite_prob {
  casadi_int ne, nx, ng, K;
  casadi_int nnz_au, nnz_hu;
  const casadi_int *x, *row, *ent, *mflat, *m, *slot;
  const casadi_int *a_src, *h_src;
  const char *xfree, *kind;
};
// C-REPLACE "casadi_ipmc_rewrite_prob<T1>" "struct casadi_ipmc_rewrite_prob"

// SYMBOL "ipmc_rewrite_data"
template<typename T1>
struct casadi_ipmc_rewrite_data {
  /* the CALLER's problem; the plugin's own d->nlp is the rewritten one */
  casadi_nlpsol_data<T1>* user;
  /* [p->ne] 1 once a non-vacuous row was found for that helper component */
  casadi_int* used;
};
// C-REPLACE "casadi_ipmc_rewrite_data<T1>" "struct casadi_ipmc_rewrite_data"

// SYMBOL "ipmc_rewrite_work"
template<typename T1>
void casadi_ipmc_rewrite_work(const casadi_ipmc_rewrite_prob<T1>* p, casadi_int* sz_iw) {
  *sz_iw += p->ne;  // used
}

// SYMBOL "ipmc_rewrite_set_work"
template<typename T1>
void casadi_ipmc_rewrite_set_work(casadi_ipmc_rewrite_data<T1>* d,
    const casadi_ipmc_rewrite_prob<T1>* p, casadi_int** iw) {
  d->used = *iw; *iw += p->ne;
}

// SYMBOL "ipmc_rewrite_expand"
template<typename T1>
void casadi_ipmc_rewrite_expand(const casadi_ipmc_rewrite_prob<T1>* p,
    casadi_ipmc_rewrite_data<T1>* d, casadi_nlpsol_data<T1>* v,
    const T1* ubs, const T1* s0, const T1* lam_s0) {
  casadi_int i, k, e, r, a, src, kind, slot, ne, K, nxu, nxl, ngl;
  T1 lo, up, m0, lam0;
  casadi_nlpsol_data<T1>* u = d->user;
  if (!p->ne) return;
  nxu = p->nx;
  nxl = v->prob->nx;
  ngl = v->prob->ng;
  ne = p->ne;
  K = p->K;
  v->p = u->p;
  /* caller variables, in their new places.  A softened simple bound is
     reimposed as a path row, so the variable's own bound is opened. */
  for (i=0;i<nxu;++i) {
    a = p->x[i];
    v->z[a] = u->z[i];
    if (p->xfree[i]) {
      v->lbz[a] = -std::numeric_limits<T1>::infinity();
      v->ubz[a] = std::numeric_limits<T1>::infinity();
      v->lam[a] = 0;
    } else {
      v->lbz[a] = u->lbz[i];
      v->ubz[a] = u->ubz[i];
      v->lam[a] = u->lam[i];
    }
  }
  /* helper components: free everywhere but on stage 0, where the bound
     below pins them */
  for (i=0;i<K*ne;++i) {
    a = p->m[i];
    v->z[a] = 0;
    v->lbz[a] = -std::numeric_limits<T1>::infinity();
    v->ubz[a] = std::numeric_limits<T1>::infinity();
    v->lam[a] = 0;
  }
  for (e=0;e<ne;++e) d->used[e] = 0;
  /* rows.  A split row is one-sided; the multiplier guess is divided over
     the two halves by sign. */
  for (r=0;r<ngl;++r) {
    src = p->row[r];
    kind = p->kind[r];
    lo = 0;
    up = 0;
    lam0 = 0;
    if (kind!=3) {
      lo = kind==2 ? -std::numeric_limits<T1>::infinity() : u->lbz[src];
      up = kind==1 ? std::numeric_limits<T1>::infinity() : u->ubz[src];
      lam0 = u->lam[src];
      if (kind==1 && lam0>0) lam0 = 0;
      if (kind==2 && lam0<0) lam0 = 0;
    }
    v->lbz[nxl+r] = lo;
    v->ubz[nxl+r] = up;
    v->lam[nxl+r] = lam0;
    e = p->ent[r];
    /* a row vacuous once split (no finite bound on this side) leaves its
       helper component with nothing to do */
    if (e>=0 && !(lo==-std::numeric_limits<T1>::infinity()
               && up==std::numeric_limits<T1>::infinity())) {
      d->used[e] = 1;
    }
  }
  /* 0 <= m <= ubs, on the stage-0 copy.  An unused component, or one with
     ubs==0, becomes a FIXED variable instead: 0 <= m <= 0 has no barrier
     interior, and pinning rather than deleting keeps the row coupled to m
     and its multiplier a shadow price. */
  for (e=0;e<ne;++e) {
    slot = p->slot[e];
    a = p->m[e];
    m0 = 0;
    v->lbz[a] = 0;
    v->ubz[a] = 0;
    if (d->used[e] && ubs[slot]>0) {
      v->ubz[a] = ubs[slot];
      if (s0) {
        m0 = s0[slot];
        if (!(m0>0)) m0 = 0;
        if (m0>ubs[slot]) m0 = ubs[slot];
      }
      if (lam_s0) v->lam[a] = lam_s0[slot];
    }
    /* start the whole chain m_0 = ... = m_N at one value: zero initial gap */
    for (k=0;k<K;++k) v->z[p->m[k*ne+e]] = m0;
  }
}

// SYMBOL "ipmc_rewrite_collect"
template<typename T1>
void casadi_ipmc_rewrite_collect(const casadi_ipmc_rewrite_prob<T1>* p,
    casadi_ipmc_rewrite_data<T1>* d, casadi_nlpsol_data<T1>* v,
    T1* s, T1* lam_s) {
  casadi_int i, e, r, a, slot, ne, nxu, ngu, nxl, ngl;
  casadi_nlpsol_data<T1>* u = d->user;
  if (!p->ne) return;
  nxu = p->nx;
  ngu = p->ng;
  nxl = v->prob->nx;
  ngl = v->prob->ng;
  ne = p->ne;
  for (i=0;i<nxu;++i) {
    a = p->x[i];
    u->z[i] = v->z[a];
    u->lam[i] = p->xfree[i] ? 0 : v->lam[a];
  }
  for (i=0;i<ngu;++i) u->lam[nxu+i] = 0;
  /* at most one half of a split row can be active, so summing the halves is
     the same statement the reference expansion makes */
  for (r=0;r<ngl;++r) {
    if (p->kind[r]!=3) u->lam[p->row[r]] += v->lam[nxl+r];
  }
  /* m_0 IS the slack value, and its simple-bound multiplier already carries
     casadi's sign convention */
  for (e=0;e<ne;++e) {
    slot = p->slot[e];
    a = p->m[e];
    s[slot] = v->z[a];
    lam_s[slot] = v->lam[a];
  }
}

// SYMBOL "ipmc_rewrite_collect_g"
/* The rewritten CONSTRAINT VALUES back into the caller's g (the tail of
   u->z).  An untouched row IS the caller's row.  A split half carries the
   helper it was relaxed against (+m lower, -m upper), so taking that off
   recovers the caller's value; either half works.  A softened simple bound
   and a continuity row have no caller g row.  The helper is read at the
   row's own stage (p->mflat): the copies agree only where the continuity
   rows hold, and a reported iterate has to be the point it actually is.
   Called only for intermediate iterates; at the end of a solve Nlpsol
   evaluates the caller's own oracle instead. */
template<typename T1>
void casadi_ipmc_rewrite_collect_g(const casadi_ipmc_rewrite_prob<T1>* p,
    casadi_ipmc_rewrite_data<T1>* d, casadi_nlpsol_data<T1>* v) {
  casadi_int r, src, kind, nxl, ngl;
  casadi_nlpsol_data<T1>* u = d->user;
  if (!p->ne) return;
  nxl = v->prob->nx;
  ngl = v->prob->ng;
  for (r=0;r<ngl;++r) {
    kind = p->kind[r];
    src = p->row[r];
    if (kind==3 || src<p->nx) continue;
    if (kind==1) {
      u->z[src] = v->z[nxl+r] - v->z[p->m[p->mflat[r]]];
    } else if (kind==2) {
      u->z[src] = v->z[nxl+r] + v->z[p->m[p->mflat[r]]];
    } else {
      u->z[src] = v->z[nxl+r];
    }
  }
}

// SYMBOL "ipmc_mproject"
template<typename T1>
void casadi_ipmc_mproject(T1 factor, const T1* x, const casadi_int* sp_x,
    T1* y, const casadi_int* sp_y, T1* w) {
  casadi_int ncol_y = sp_y[1];
  const casadi_int* colind_y = sp_y+2;
  casadi_project(x, sp_x, y, sp_y, w);
  casadi_scal(colind_y[ncol_y], factor, y);
}

// C-REPLACE "casadi_ocp_block" "struct casadi_ocp_block"

template<typename T1>
struct casadi_ipmc_prob {
  const casadi_nlpsol_prob<T1>* nlp;
  const casadi_int *nx, *nu, *ng;
  /* trailing states of every stage whose dynamics are the identity; the
     same number on every stage, 0 when the user declared none */
  casadi_int nxc;
  casadi_int nx_max, nu_max, nxu_max;
  // Sparsity patterns
  const casadi_int *sp_h, *sp_a;

  casadi_int nnz_h, nnz_a;

  // Sparsities
  const casadi_int *ABsp;
  const casadi_int *AB_offsets;
  const casadi_int *CDsp;
  const casadi_int *CD_offsets;
  const casadi_int *RSQsp;
  const casadi_int *RSQ_offsets;
  const casadi_int *Isp;
  const casadi_int *I_offsets;

  casadi_int N;
  casadi_ocp_block *AB, *CD, *RSQ, *I;

  /* ---- native soft constraints (slacks); all zero when ns==0 ------------
       slack_g[i]  : casadi slack column of user constraint row g[i], -1 = hard
       slack_x[i]  : casadi slack column of user variable x[i], -1 = hard
       slack_perm  : ipmc slack index -> casadi slack column
       slack_offs  : [N+2] prefix sums, slack_offs[k+1]-slack_offs[k] columns
                     live on stage k
     casadi's slack vector is [s_l; s_u], each 'ns' long, so the upper-side
     quantity of column j sits at index ns+j. */
  casadi_int ns;
  const casadi_int *slack_g;
  const casadi_int *slack_x;
  const casadi_int *slack_perm;
  const casadi_int *slack_offs;
  /* f_s = c + z's + 1/2 s'Zs, validated at construction to be a separable
     quadratic independent of p.  fs_z is its gradient at s=0 and fs_Z its
     Hessian diagonal, both dense of length 2*ns, evaluated once by
     IpmcInterface::init.  Nothing here calls f_s. */
  const T1 *fs_z;
  const T1 *fs_Z;

  /* The interface may hand over a problem that is not the caller's but a
     larger equivalent (IpmcInterface lifts a cross-stage slack column that
     way).  Everything above then describes the REWRITTEN problem, and this
     carries the caller's z-space onto it and back.  rewrite.ne==0 means
     nothing was rewritten; see the rewrite block above. */
  casadi_ipmc_rewrite_prob<T1> rewrite;

  IpmcWrite write;
  IpmcFlush flush;
};
// C-REPLACE "casadi_ipmc_prob<T1>" "struct casadi_ipmc_prob"

/* The ipmc ux vector of stage k is [u_k ; x_k].  These two helpers translate
   between that layout and casadi's decision vector, whose stage-k offset is
   CD[k].offset_c. */

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

// SYMBOL "ipmc_setup"
template<typename T1>
void casadi_ipmc_setup(casadi_ipmc_prob<T1>* p) {
  casadi_int k;
  if (p->sp_h) {
    p->nnz_h = p->sp_h[2+p->sp_h[1]];
  } else {
    p->nnz_h = 0;
  }
  p->nnz_a = p->sp_a[2+p->sp_a[1]];

  p->nx_max = 0;
  for (k=0;k<p->N+1;++k) {
    if (p->nx[k]>p->nx_max) p->nx_max = p->nx[k];
  }
  p->nu_max = 0;
  p->nxu_max = 0;
  for (k=0;k<p->N;++k) {
    if (p->nu[k]>p->nu_max) p->nu_max = p->nu[k];
    if (p->nu[k]+p->nx[k]>p->nxu_max) p->nxu_max = p->nu[k]+p->nx[k];
  }
}

// SYMBOL "ipmc_data"
template<typename T1>
struct casadi_ipmc_data {
  // Problem structure
  const casadi_ipmc_prob<T1>* prob;
  /* The REWRITTEN problem when p->rewrite.ne (rewrite.user is then the
     caller's); the same pointer otherwise. */
  casadi_nlpsol_data<T1>* nlp;
  casadi_ipmc_rewrite_data<T1> rewrite;

  T1 *AB, *CD, *RSQ, *I;

  casadi_int *a_eq, *a_ineq, *a_eq_idx, *a_ineq_idx;
  casadi_int *x_eq, *x_ineq, *x_eq_idx, *x_ineq_idx;

  const T1** arg;
  T1** res;
  casadi_int* iw;
  T1* w;

  int unified_return_status;
  int success;
  int return_status;

  T1 *pv, *x, *a, *g, *h, *lam;

  /* The caller oracle's buffers on the lifted path: point, values/gradient,
     jacobian and hessian nonzeros, folded multipliers.  Aliases of
     x/g/a/h/lam when nothing was rewritten, so the solve loop uses them
     unconditionally and the lift_* routines are then no-ops. */
  T1 *xu, *gu, *au, *hu, *lamu;

  /* ---- native soft constraints ----------------------------------------
     slack_s/slack_lam_s/slack_ubs point straight at NlpsolMemory::slack_*
     (the codegen mirror uses the nlp_slack_* locals).  Never null when
     p->ns>0; 2*ns doubles each, laid out [lower; upper].  m->s / m->lam_s
     are NOT usable here: they may be null. */
  T1 *slack_s, *slack_lam_s, *slack_ubs;
  T1 *sl, *su, *zsl, *zsu;
  /* Rebuilt on every presolve: they depend on ubs and on the eq/ineq
     classification, both solve-time inputs. */
  casadi_int *soft_ns, *soft_offs, *soft_map, *soft_local;
  casadi_int *soft_idxs;
  /* The problem STRUCTURE, in the shape ipmc_create takes it: nu, nx, ng,
     ng_ineq and ns per stage, then the row -> slack map.  `structure` is
     what this presolve wants, `structure_built` what the existing solver was
     built from; they differ exactly when the solver has to be rebuilt.
     Nothing else forces a rebuild: bounds, starting point, penalty and
     initial slacks are set on the live solver and re-read at every
     ipmc_start.  ipmc_int because ipmc_create reads them; malloc'ed rather
     than carved out of iw because iw is handed out fresh and zero-filled on
     every call (alloc_iw's persistent=true does not survive between calls).
     Both live on the memory object next to `solver`, freed in free_mem. */
  ipmc_int *structure, *structure_built;
  casadi_int sz_structure;
  casadi_int n_soft;
  int soft_error;
  /* 1 once the identity dynamics PROMISED for the constant states (nxc is a
     declaration, not a detection) have been checked against the Jacobian
     this solve.  Checked once per solve, at the first Jacobian; a value that
     changes later is a nonlinearity the first iterate already shows.
     Cleared by presolve. */
  int nxc_checked;
  /* set by presolve: 1 iff it created a new solver object.  IpmcInterface
     ::solve() pushes the plugin options exactly when this is set -- reading
     `d->solver == 0` BEFORE presolve is wrong, because presolve may destroy
     and re-create the solver, which would then never see the options. */
  int solver_created;

  struct blasfeo_dvec v, r;
  struct blasfeo_dmat R;

  /* Scratch for the numbers handed to ipmc_set_bounds / _initial /
     _soft_penalty / _initial_slack: plain arrays in ipmc's own layout,
     built stage after stage in presolve. */
  T1 *set_lower, *set_upper, *set_ux0, *set_pen, *set_sl0, *set_su0;

  /* Constraint values snapshotted for the progress hook, with the point
     they were evaluated at (ipmc may step back to an earlier iterate).
     ipmc's own callback data cannot yield casadi's g: it carries the
     VIOLATION, folded into its stagewise layout.  Both are quantities of
     the problem that was HANDED OVER (cb_g is p->nlp->ng long) and reach
     the caller's z-space through casadi_ipmc_report_iterate.  Only written
     while a hook is installed: the C++ interface points it at
     Nlpsol::callback, generated C has no Function to call back into and
     leaves it null. */
  T1 *cb_g, *cb_x;

  struct IpmcStats stats;

  struct IpmcSolver *solver;
  /* ipmc's vector layout, valid for as long as `solver` is */
  const struct IpmcLayout *layout;
};
// C-REPLACE "casadi_ipmc_data<T1>" "struct casadi_ipmc_data"

// SYMBOL "ipmc_init_mem"
template<typename T1>
int casadi_ipmc_init_mem(casadi_ipmc_data<T1>* d) {
  d->solver = 0;
  d->soft_error = 0;
  d->n_soft = 0;
  d->structure = 0;
  d->structure_built = 0;
  d->sz_structure = 0;
  d->layout = 0;
  d->solver_created = 0;
  d->nxc_checked = 0;
  return 0;
}

// SYMBOL "ipmc_free_mem"
template<typename T1>
void casadi_ipmc_free_mem(casadi_ipmc_data<T1>* d) {
  if (d->solver) {
    ipmc_destroy(d->solver);
    d->solver = 0;
    d->layout = 0;
  }
  if (d->structure) {
    free(d->structure);
    d->structure = 0;
    d->structure_built = 0;
    d->sz_structure = 0;
  }
}
// C-REPLACE "static_cast< casadi_ipmc_data<T1>* >" "(struct casadi_ipmc_data*)"
// C-REPLACE "casadi_error" "//casadi_error"

/* ==== lifting the caller oracle's outputs ================================
   On the lifted path the oracle is the CALLER's, so every request of the
   solve loop gathers the caller's point out of the rewritten one, calls the
   caller's function, and lifts the output into the rewritten problem's
   quantity.  Each routine below is one such lift; all are no-ops when
   nothing was rewritten (and the u-buffers then alias the lifted ones). */

// SYMBOL "ipmc_lift_x"
/* the caller's point, gathered out of the rewritten one (d->x -> d->xu) */
template<typename T1>
void casadi_ipmc_lift_x(const casadi_ipmc_prob<T1>* p, casadi_ipmc_data<T1>* d) {
  casadi_int i;
  const casadi_ipmc_rewrite_prob<T1>* q = &p->rewrite;
  if (!q->ne) return;
  for (i=0;i<q->nx;++i) d->xu[i] = d->x[q->x[i]];
}

// SYMBOL "ipmc_lift_obj"
/* the penalty on the stage-0 helpers at the point in d->x, onto the
   caller's objective value; the scaling happens outside */
template<typename T1>
void casadi_ipmc_lift_obj(const casadi_ipmc_prob<T1>* p, casadi_ipmc_data<T1>* d, T1* obj) {
  casadi_int e, slot;
  T1 m0;
  const casadi_ipmc_rewrite_prob<T1>* q = &p->rewrite;
  for (e=0;e<q->ne;++e) {
    slot = q->slot[e];
    m0 = d->x[q->m[e]];
    *obj += p->fs_z[slot]*m0 + 0.5*p->fs_Z[slot]*m0*m0;
  }
}

// SYMBOL "ipmc_lift_grad_f"
/* grad f of the rewritten problem out of the caller's (d->gu -> d->g): the
   caller's entries in their new places, the penalty gradient on the stage-0
   helpers, zero on every other helper copy */
template<typename T1>
void casadi_ipmc_lift_grad_f(const casadi_ipmc_prob<T1>* p, casadi_ipmc_data<T1>* d) {
  casadi_int i, e, slot;
  const casadi_ipmc_rewrite_prob<T1>* q = &p->rewrite;
  if (!q->ne) return;
  for (i=0;i<q->nx;++i) d->g[q->x[i]] = d->gu[i];
  for (i=0;i<q->K*q->ne;++i) d->g[q->m[i]] = 0;
  for (e=0;e<q->ne;++e) {
    slot = q->slot[e];
    d->g[q->m[e]] = p->fs_z[slot] + p->fs_Z[slot]*d->x[q->m[e]];
  }
}

// SYMBOL "ipmc_lift_g"
/* g of the rewritten problem out of the caller's (d->gu -> d->g), at the
   point in d->x: a split half carries +-m, a softened simple bound becomes
   x +- m, a continuity row is m_{k+1} - m_k */
template<typename T1>
void casadi_ipmc_lift_g(const casadi_ipmc_prob<T1>* p, casadi_ipmc_data<T1>* d) {
  casadi_int r, src, kind, fl, ngl;
  const casadi_ipmc_rewrite_prob<T1>* q = &p->rewrite;
  if (!q->ne) return;
  ngl = p->nlp->ng;
  for (r=0;r<ngl;++r) {
    kind = q->kind[r];
    fl = q->mflat[r];
    if (kind==3) {
      d->g[r] = d->x[q->m[fl+q->ne]] - d->x[q->m[fl]];
    } else {
      src = q->row[r];
      d->g[r] = src>=q->nx ? d->gu[src-q->nx] : d->x[q->x[src]];
      if (kind==1) d->g[r] += d->x[q->m[fl]];
      else if (kind==2) d->g[r] -= d->x[q->m[fl]];
    }
  }
}

// SYMBOL "ipmc_lift_jac_g"
/* the rewritten jacobian nonzeros out of the caller's (d->au -> d->a),
   through the construction-time map; the negative entries are the +-1
   coefficients of the helper, split-bound and continuity columns */
template<typename T1>
void casadi_ipmc_lift_jac_g(const casadi_ipmc_prob<T1>* p, casadi_ipmc_data<T1>* d) {
  casadi_int i, s;
  const casadi_ipmc_rewrite_prob<T1>* q = &p->rewrite;
  if (!q->ne) return;
  for (i=0;i<p->nnz_a;++i) {
    s = q->a_src[i];
    d->a[i] = s>=0 ? d->au[s] : (s==-1 ? 1.0 : -1.0);
  }
}

// SYMBOL "ipmc_lift_lam"
/* the caller's constraint multipliers folded out of the rewritten ones in
   d->lam (what nlp_hess_l wants for lam:g): both halves of a split row sum
   onto the caller's row; bound and continuity rows have no caller row */
template<typename T1>
void casadi_ipmc_lift_lam(const casadi_ipmc_prob<T1>* p, casadi_ipmc_data<T1>* d) {
  casadi_int r, src, ngl;
  const casadi_ipmc_rewrite_prob<T1>* q = &p->rewrite;
  if (!q->ne) return;
  for (r=0;r<q->ng;++r) d->lamu[r] = 0;
  ngl = p->nlp->ng;
  for (r=0;r<ngl;++r) {
    src = q->row[r];
    if (q->kind[r]!=3 && src>=q->nx) d->lamu[src-q->nx] += d->lam[r];
  }
}

// SYMBOL "ipmc_lift_hess_l"
/* the rewritten Lagrangian gradient and Hessian out of the caller's
   (d->gu / d->hu -> d->g / d->h).  d->gu is grad:gamma:x of the CALLER's
   problem at the folded multipliers, so the caller part is done; added here
   is what the caller's oracle cannot know: the multipliers of the rows with
   no caller g row on the variables they touch, the +-m columns of the split
   rows, and the penalty.  The Hessian gains only the penalty diagonal --
   every helper term is linear. */
template<typename T1>
void casadi_ipmc_lift_hess_l(const casadi_ipmc_prob<T1>* p, casadi_ipmc_data<T1>* d,
    T1 obj_scale) {
  casadi_int i, r, e, src, kind, fl, slot, s, ngl;
  T1 lam;
  const casadi_ipmc_rewrite_prob<T1>* q = &p->rewrite;
  if (!q->ne) return;
  for (i=0;i<q->nx;++i) d->g[q->x[i]] = d->gu[i];
  for (i=0;i<q->K*q->ne;++i) d->g[q->m[i]] = 0;
  ngl = p->nlp->ng;
  for (r=0;r<ngl;++r) {
    kind = q->kind[r];
    if (kind==0) continue;
    lam = d->lam[r];
    fl = q->mflat[r];
    if (kind==3) {
      d->g[q->m[fl+q->ne]] += lam;
      d->g[q->m[fl]] -= lam;
    } else {
      src = q->row[r];
      if (src<q->nx) d->g[q->x[src]] += lam;
      if (kind==1) d->g[q->m[fl]] += lam;
      else d->g[q->m[fl]] -= lam;
    }
  }
  for (e=0;e<q->ne;++e) {
    slot = q->slot[e];
    d->g[q->m[e]] += obj_scale*(p->fs_z[slot] + p->fs_Z[slot]*d->x[q->m[e]]);
  }
  for (i=0;i<p->nnz_h;++i) {
    s = q->h_src[i];
    d->h[i] = s>=0 ? d->hu[s] : obj_scale*p->fs_Z[q->slot[-2-s]];
  }
}

// SYMBOL "ipmc_snapshot_g"
/* The constraint values put aside for the iteration callback, tagged with
   the point they were evaluated at.  d->g is shared scratch (nlp_grad_f and
   nlp_jac_g write it too), so g has to be copied while it is still there;
   cb_x lets the post-iteration path tell whether the snapshot belongs to
   the iterate being reported.  Both callers evaluate g right before this. */
template<typename T1>
void casadi_ipmc_snapshot_g(casadi_ipmc_data<T1>* d, const struct IpmcLayout* s,
            const double* primal_data) {
  const casadi_ipmc_prob<T1>* p = d->prob;
  casadi_int n_ux = s->ux_offs[s->K-1] + s->nu[s->K-1] + s->nx[s->K-1];

  casadi_copy(d->g, p->nlp->ng, d->cb_g);
  casadi_copy(primal_data, n_ux, d->cb_x);
}

// SYMBOL "ipmc_pack_constr_viol"
/* The violation vector ipmc reads, out of constraint values an nlp_g has
   ALREADY written.  Pack, not eval: the nlp_g call is the caller's.  The
   rows ipmc holds are not the caller's rows: a two-sided casadi row shows
   up once, an equality row is measured against its own bound and a simple
   bound on x is a row of its own. */
template<typename T1>
void casadi_ipmc_pack_constr_viol(casadi_ipmc_data<T1>* d, const struct IpmcLayout* s,
            double* res) {
  casadi_int i,k,column;
  const casadi_ipmc_prob<T1>* p = d->prob;
  casadi_nlpsol_data<T1>* d_nlp = d->nlp;

  for (k=0;k<s->K;++k) {
    column = 0;
    for (i=d->a_ineq_idx[k];i<d->a_ineq_idx[k+1];++i) {
      res[s->g_ineq_offs[k]+column] = d->g[d->a_ineq[i]];
      column++;
    }
    for (i=d->x_ineq_idx[k];i<d->x_ineq_idx[k+1];++i) {
      res[s->g_ineq_offs[k]+column] = d->x[d->x_ineq[i]];
      column++;
    }
    column = 0;
    for (i=d->a_eq_idx[k];i<d->a_eq_idx[k+1];++i) {
      res[s->g_offs[k]+column] = d->g[d->a_eq[i]]-d_nlp->lbz[p->nlp->nx+d->a_eq[i]];
      column++;
    }
    for (i=d->x_eq_idx[k];i<d->x_eq_idx[k+1];++i) {
      res[s->g_offs[k]+column] = d->x[d->x_eq[i]]-d_nlp->lbz[d->x_eq[i]];
      column++;
    }
  }
  for (k=0;k<s->K-1;++k) {
    const T1* lbg_k = d_nlp->lbz+p->nlp->nx+p->AB[k].offset_r;
    const T1* g_k = d->g+p->AB[k].offset_r;
    casadi_copy(lbg_k, p->nx[k+1], res+s->dyn_eq_offs[k]);
    casadi_axpy(p->nx[k+1], -1.0, g_k, res+s->dyn_eq_offs[k]);
  }
}

// C-REPLACE "const_cast<T1*>" "(T1*)"

// SYMBOL "ipmc_pack_BAbt"
template<typename T1>
void casadi_ipmc_pack_BAbt(casadi_ipmc_data<T1>* d, struct blasfeo_dmat *res,
    const ipmc_int k) {
  const casadi_ipmc_prob<T1>* p = d->prob;
  casadi_nlpsol_data<T1>* d_nlp = d->nlp;
  const T1* lbg_k = d_nlp->lbz+p->nlp->nx+p->AB[k].offset_r;
  const T1* g_k = d->g+p->AB[k].offset_r;
  casadi_int i, j;
  /* The identity block ipmc exploits, the same on every transition and
     trailing-aligned at both ends.  nxd1 is what is left to factorize. */
  casadi_int nc1 = p->nxc;
  casadi_int nxd1 = p->nx[k+1]-nc1;
  blasfeo_pack_tran_dmat(nxd1, p->nx[k], d->AB+p->AB_offsets[k], p->nx[k+1], res, p->nu[k], 0);
  blasfeo_pack_tran_dmat(nxd1, p->nu[k], d->AB+p->AB_offsets[k]+p->nx[k]*p->nx[k+1], p->nx[k+1], res, 0, 0);

  /* The gap row stays full width: the constant rows have their own dynamics
     residual x_{k+1}-x_k, which is zero only at a feasible point. */
  for (i=0; i<p->nx[k+1]; ++i) {
    BLASFEO_DMATEL(res, p->nu[k]+p->nx[k], i) = lbg_k[i]-g_k[i];
  }

  if (nc1>0 && !d->nxc_checked) {
    /* Row nxd1+i of this block must be e_{nxd0+i}: the declared constant
       states have to be the ones the problem actually holds constant. */
    casadi_int nxd0 = p->nx[k]-nc1;
    const T1* AB_k = d->AB+p->AB_offsets[k];
    for (i=0; i<nc1; ++i) {
      for (j=0; j<p->nx[k]+p->nu[k]; ++j) {
        if (AB_k[j*p->nx[k+1]+nxd1+i] != (j==nxd0+i ? 1.0 : 0.0)) {
          /* The raise stays on ONE line: codegen rewrites it into a line
             comment, and a continuation line would survive as bare code. */
          casadi_error("nxc: stage " + str(k+1) + " state " + str(nxd1+i) + " was declared constant, but its dynamics row is not x_{k+1}=x_k.");
        }
      }
    }
  }
}

// SYMBOL "ipmc_pack_RSQrqt"
template<typename T1>
void casadi_ipmc_pack_RSQrqt(
    casadi_ipmc_data<T1>* d,
    struct blasfeo_dmat *res,
    const ipmc_int k) {
  const casadi_ipmc_prob<T1>* p = d->prob;

  int n = p->nx[k]+p->nu[k];
  blasfeo_pack_dmat(p->nx[k], p->nx[k],
    d->RSQ+p->RSQ_offsets[k], n, res, p->nu[k], p->nu[k]);
  blasfeo_pack_dmat(p->nu[k], p->nu[k],
    d->RSQ+p->RSQ_offsets[k]+p->nx[k]*n+p->nx[k], n, res, 0, 0);
  blasfeo_pack_dmat(p->nu[k], p->nx[k],
    d->RSQ+p->RSQ_offsets[k]+p->nx[k], n, res, 0, p->nu[k]);
  blasfeo_pack_dmat(p->nx[k], p->nu[k],
    d->RSQ+p->RSQ_offsets[k]+p->nx[k]*n, n, res, p->nu[k], 0);

  blasfeo_pack_dmat(1, p->nx[k], d->g+p->CD[k].offset_c, 1, res, p->nu[k]+p->nx[k], p->nu[k]);
  blasfeo_pack_dmat(1, p->nu[k], d->g+p->CD[k].offset_c+p->nx[k], 1, res, p->nu[k]+p->nx[k], 0);
}

// SYMBOL "ipmc_read_lam"
/* ipmc's stage-by-stage multipliers, scattered into the caller's row order
   in d->lam -- what nlp_hess_l wants for its lam argument, so this runs
   BEFORE the oracle call.  The dynamics multipliers change sign on the way:
   ipmc's rows are the gaps x_{k+1}-F(x_k,u_k), casadi's are
   F(x_k,u_k)-x_{k+1}. */
template<typename T1>
void casadi_ipmc_read_lam(casadi_ipmc_data<T1>* d, const struct IpmcLayout* s,
            const double* lam_data) {
  casadi_int k, column, i;
  const casadi_ipmc_prob<T1>* p = d->prob;

  for (k=0;k<s->K;++k) {
    column = 0;
    for (i=d->a_ineq_idx[k];i<d->a_ineq_idx[k+1];++i) {
      d->lam[d->a_ineq[i]] = lam_data[s->g_ineq_offs[k]+column];
      column++;
    }
    column = 0;
    for (i=d->a_eq_idx[k];i<d->a_eq_idx[k+1];++i) {
      d->lam[d->a_eq[i]] = lam_data[s->g_offs[k]+column];
      column++;
    }
  }
  for (k=0;k<s->K-1;++k) {
    casadi_scaled_copy(-1.0, lam_data+s->dyn_eq_offs[k], p->nx[k+1], d->lam+p->AB[k].offset_r);
  }
}

// SYMBOL "ipmc_pack_lag_hess"
/* The RSQrqt blocks ipmc reads, out of a Hessian and a Lagrangian gradient
   an nlp_hess_l has ALREADY written into d->h and d->g.  Pack, not eval.
   The gradient needs finishing first: the Lagrangian casadi differentiates
   carries lam*F(x_k,u_k) without the x_{k+1} term, and the multipliers of
   the rows ipmc holds as simple bounds are not in it at all, so both are
   added back here before d->g is the last row of RSQrqt. */
template<typename T1>
void casadi_ipmc_pack_lag_hess(casadi_ipmc_data<T1>* d, const struct IpmcLayout* s,
            const double* lam_data, struct blasfeo_dmat* res) {
  casadi_int k, column, i;
  const casadi_ipmc_prob<T1>* p = d->prob;

  casadi_project(d->h, p->sp_h, d->RSQ, p->RSQsp, d->pv);

  /* the missing x_{k+1} term of the dynamics rows */
  for (k=0;k<s->K-1;++k) {
    casadi_axpy(p->nx[k+1], 1.0, lam_data+s->dyn_eq_offs[k], d->g+p->CD[k+1].offset_c);
  }
  /* the simple-bound multipliers */
  for (k=0;k<s->K;++k) {
    column = d->a_ineq_idx[k+1]-d->a_ineq_idx[k];
    for (i=d->x_ineq_idx[k];i<d->x_ineq_idx[k+1];++i) {
      d->g[d->x_ineq[i]] += lam_data[s->g_ineq_offs[k]+column];
      column++;
    }
    column = d->a_eq_idx[k+1]-d->a_eq_idx[k];
    for (i=d->x_eq_idx[k];i<d->x_eq_idx[k+1];++i) {
      d->g[d->x_eq[i]] += lam_data[s->g_offs[k]+column];
      column++;
    }
  }

  for (k = 0; k < s->K; k++) {
    casadi_ipmc_pack_RSQrqt(d, res + k, k);
  }
}

// SYMBOL "ipmc_pack_Ggt"
template<typename T1>
void casadi_ipmc_pack_Ggt(
      casadi_ipmc_data<T1>* d,
      struct blasfeo_dmat *res,
      const ipmc_int k) {
  const casadi_ipmc_prob<T1>* p = d->prob;
  casadi_int i, column;

  int n_a_eq = d->a_eq_idx[k+1]-d->a_eq_idx[k];
  int n_x_eq = d->x_eq_idx[k+1]-d->x_eq_idx[k];
  int ng_eq = n_a_eq+n_x_eq;

  blasfeo_dgese(p->nu[k]+p->nx[k]+1, ng_eq, 0.0, res, 0, 0);

  column = 0;
  for (i=d->a_eq_idx[k];i<d->a_eq_idx[k+1];++i) {
    blasfeo_pack_tran_dmat(1, p->nx[k],
      d->CD+p->CD_offsets[k]+(d->a_eq[i]-p->CD[k].offset_r),
      p->CD[k].rows, res, p->nu[k], column);
    blasfeo_pack_tran_dmat(1, p->nu[k],
      d->CD+p->CD_offsets[k]+(d->a_eq[i]-p->CD[k].offset_r)+p->nx[k]*p->CD[k].rows,
      p->CD[k].rows, res, 0,        column);
    BLASFEO_DMATEL(res, p->nu[k]+p->nx[k], column) = d->g[d->a_eq[i]]-d->nlp->lbz[p->nlp->nx+d->a_eq[i]];
    column++;
  }
  for (i=d->x_eq_idx[k];i<d->x_eq_idx[k+1];++i) {
    int j = d->x_eq[i]-p->CD[k].offset_c;
    if (j>=p->nx[k]) {
      j -= p->nx[k];
    } else {
      j += p->nu[k];
    }
    BLASFEO_DMATEL(res, j, column) = 1;
    BLASFEO_DMATEL(res, p->nu[k]+p->nx[k], column) = d->x[d->x_eq[i]]-d->nlp->lbz[d->x_eq[i]];
    column++;
  }
}

// SYMBOL "ipmc_pack_Ggt_ineq"
template<typename T1>
void casadi_ipmc_pack_Ggt_ineq(
    casadi_ipmc_data<T1>* d,
    struct blasfeo_dmat *res,
    const ipmc_int k) {
  const casadi_ipmc_prob<T1>* p = d->prob;
  casadi_int i, column;

  int n_a_ineq = d->a_ineq_idx[k+1]-d->a_ineq_idx[k];
  int n_x_ineq = d->x_ineq_idx[k+1]-d->x_ineq_idx[k];
  int ng_ineq = n_a_ineq+n_x_ineq;

  blasfeo_dgese(p->nu[k]+p->nx[k]+1, ng_ineq, 0.0, res, 0, 0);

  column = 0;
  for (i=d->a_ineq_idx[k];i<d->a_ineq_idx[k+1];++i) {
    blasfeo_pack_tran_dmat(1, p->nx[k],
      d->CD+p->CD_offsets[k]+(d->a_ineq[i]-p->CD[k].offset_r),
      p->CD[k].rows, res, p->nu[k], column);
    blasfeo_pack_tran_dmat(1, p->nu[k],
      d->CD+p->CD_offsets[k]+(d->a_ineq[i]-p->CD[k].offset_r)+p->nx[k]*p->CD[k].rows,
      p->CD[k].rows, res, 0,        column);
    BLASFEO_DMATEL(res, p->nu[k]+p->nx[k], column) = d->g[d->a_ineq[i]];
    column++;
  }
  for (i=d->x_ineq_idx[k];i<d->x_ineq_idx[k+1];++i) {
    int j = d->x_ineq[i]-p->CD[k].offset_c;
    if (j>=p->nx[k]) {
      j -= p->nx[k];
    } else {
      j += p->nu[k];
    }
    BLASFEO_DMATEL(res, j, column) = 1;
    BLASFEO_DMATEL(res, p->nu[k]+p->nx[k], column) = d->x[d->x_ineq[i]];
    column++;
  }
}

// SYMBOL "ipmc_pack_constr_jac"
/* The constraint jacobian into the block layout ipmc reads, out of nonzeros
   an nlp_jac_g has ALREADY written -- d->x, d->g and d->a are one and the
   same point on entry.  Pack, not eval: the nonzeros are scattered into the
   dense AB/CD/I blocks, then packed stage by stage into BAbt, Ggt and
   Ggt_ineq. */
template<typename T1>
void casadi_ipmc_pack_constr_jac(casadi_ipmc_data<T1>* d, const struct IpmcLayout* s,
            struct blasfeo_dmat* BAbt_p, struct blasfeo_dmat* Ggt_p,
            struct blasfeo_dmat* Ggt_ineq_p) {
  casadi_int i, k;
  const casadi_ipmc_prob<T1>* p = d->prob;

  casadi_ipmc_mproject(-1.0, d->a, p->sp_a, d->AB, p->ABsp, d->pv);
  casadi_project(d->a, p->sp_a, d->CD, p->CDsp, d->pv);
  casadi_project(d->a, p->sp_a, d->I, p->Isp, d->pv);

  for (i=0;i<p->Isp[2+p->Isp[1]];++i) {
    if (d->I[i]!=1.0) {
      casadi_error("Structure mismatch: gap-closing constraints must be like this: x_{k+1}-F(xk,uk).");
    }
  }

  for (k = 0; k < s->K-1; k++) {
    casadi_ipmc_pack_BAbt(d, BAbt_p + k, k);
  }
  d->nxc_checked = 1;

  for (k = 0; k < s->K; k++) {
    if (s->ng[k]>0) casadi_ipmc_pack_Ggt(d, Ggt_p + k, k);
  }

  for (k = 0; k < s->K; k++) {
    if (s->ng_ineq[k]>0) casadi_ipmc_pack_Ggt_ineq(d, Ggt_ineq_p + k, k);
  }
}

// SYMBOL "ipmc_work"
template<typename T1>
void casadi_ipmc_work(const casadi_ipmc_prob<T1>* p, casadi_int* sz_arg, casadi_int* sz_res, casadi_int* sz_iw, casadi_int* sz_w) {
  casadi_nlpsol_work(p->nlp, sz_arg, sz_res, sz_iw, sz_w);

  // Temporary work vectors
  *sz_w = casadi_max(*sz_w, 2*(p->nlp->nx+p->nlp->ng)); // pv

  // Persistent work vectors
  *sz_w += casadi_sp_nnz(p->ABsp); // AB
  *sz_w += casadi_sp_nnz(p->CDsp); // CD
  *sz_w += casadi_sp_nnz(p->RSQsp); // RSQ
  *sz_w += casadi_sp_nnz(p->Isp); // I
  *sz_w += p->nlp->nx;  // x
  *sz_w += p->nlp->nx+p->nlp->ng; // lam
  *sz_w += casadi_sp_nnz(p->sp_a);  // a
  *sz_w += casadi_sp_nnz(p->sp_h);  // h
  *sz_w += casadi_max(p->nlp->nx,p->nlp->ng);  // g
  *sz_w += blasfeo_memsize_dvec(p->nxu_max+1)+64; // v
  *sz_w += blasfeo_memsize_dvec(p->nx_max+p->nlp->ng)+64; // r
  *sz_w += blasfeo_memsize_dmat(p->nxu_max, p->nxu_max)+64; // R

  *sz_iw += p->N+2; // a_eq_idx
  *sz_iw += p->N+2; // a_ineq_idx
  *sz_iw += p->N+2; // x_eq_idx
  *sz_iw += p->N+2; // x_ineq_idx

  *sz_iw += p->nlp->ng; // a_eq
  *sz_iw += p->nlp->ng; // a_ineq
  *sz_iw += p->nlp->nx; // x_eq
  *sz_iw += p->nlp->nx; // x_ineq

  casadi_int n_ineq_slots = p->nlp->ng + p->nlp->nx;

  /* What is handed to ipmc_set_bounds / _initial / _soft_penalty /
     _initial_slack: plain arrays in ipmc's own layout */
  *sz_w += 2*n_ineq_slots;        // set_lower, set_upper
  *sz_w += p->nlp->nx;            // set_ux0
  *sz_w += 8*p->ns;               // set_pen, set_sl0, set_su0
  *sz_w += p->nlp->ng;            // cb_g
  *sz_w += p->nlp->nx;            // cb_x, ipmc's primal length

  // Native soft constraints
  *sz_w += 4*p->ns;               // sl, su, zsl, zsu
  *sz_iw += p->ns>0 ? p->N+1 : 0; // soft_ns
  *sz_iw += p->ns>0 ? p->N+2 : 0; // soft_offs
  *sz_iw += p->ns;                // soft_map
  *sz_iw += p->ns;                // soft_local
  *sz_iw += p->ns>0 ? n_ineq_slots : 0; // soft_idxs

  // The caller oracle's own buffers, only on the lifted path
  if (p->rewrite.ne) {
    *sz_w += p->rewrite.nx;                             // xu
    *sz_w += casadi_max(p->rewrite.nx, p->rewrite.ng);  // gu
    *sz_w += p->rewrite.nnz_au;                         // au
    *sz_w += p->rewrite.nnz_hu;                         // hu
    *sz_w += p->rewrite.ng;                             // lamu
  }
  casadi_ipmc_rewrite_work(&p->rewrite, sz_iw);
  /* d->structure is NOT here: half of it has to survive between calls, see
     the comment on casadi_ipmc_data */
}

// SYMBOL "ipmc_set_work"
template<typename T1>
void casadi_ipmc_set_work(casadi_ipmc_data<T1>* d, const T1*** arg, T1*** res, casadi_int** iw, T1** w) {
  // Problem structure
  const casadi_ipmc_prob<T1>* p = d->prob;

  d->AB = *w; *w += casadi_sp_nnz(p->ABsp);
  d->CD = *w; *w += casadi_sp_nnz(p->CDsp);
  d->RSQ = *w; *w += casadi_sp_nnz(p->RSQsp);
  d->I = *w; *w += casadi_sp_nnz(p->Isp);
  d->x = *w; *w += p->nlp->nx;
  d->lam = *w; *w += p->nlp->nx+p->nlp->ng;
  d->a = *w; *w += casadi_sp_nnz(p->sp_a);
  d->h = *w; *w += casadi_sp_nnz(p->sp_h);
  d->g = *w; *w += casadi_max(p->nlp->nx,p->nlp->ng);
  blasfeo_create_dvec(p->nxu_max+1, &d->v, (void*) (((unsigned long long) (*w)+63)/64*64));
  *w += blasfeo_memsize_dvec(p->nxu_max+1)+64;
  blasfeo_create_dvec(p->nx_max+p->nlp->ng, &d->r, (void*) (((unsigned long long) (*w)+63)/64*64));
  *w += blasfeo_memsize_dvec(p->nx_max+p->nlp->ng)+64;
  blasfeo_create_dmat(p->nxu_max, p->nxu_max, &d->R, (void*) (((unsigned long long) (*w)+63)/64*64));
  *w += blasfeo_memsize_dmat(p->nxu_max, p->nxu_max)+64;

  d->a_eq_idx = *iw;   *iw += p->N+2;
  d->a_ineq_idx = *iw; *iw += p->N+2;
  d->x_eq_idx = *iw;   *iw += p->N+2;
  d->x_ineq_idx = *iw; *iw += p->N+2;

  casadi_int n_ineq_slots = p->nlp->ng + p->nlp->nx;

  d->a_eq = *iw;   *iw += p->nlp->ng;
  d->a_ineq = *iw; *iw += p->nlp->ng;
  d->x_eq = *iw;  *iw += p->nlp->nx;
  d->x_ineq = *iw; *iw += p->nlp->nx;

  d->set_lower = *w; *w += n_ineq_slots;
  d->set_upper = *w; *w += n_ineq_slots;
  d->set_ux0 = *w;   *w += p->nlp->nx;
  d->set_pen = *w;   *w += 6*p->ns;
  d->set_sl0 = *w;   *w += p->ns;
  d->set_su0 = *w;   *w += p->ns;
  d->cb_g = *w;      *w += p->nlp->ng;
  d->cb_x = *w;      *w += p->nlp->nx;

  /* Native soft constraints -- must consume w/iw in exactly the order
     casadi_ipmc_work accounts for them */
  d->sl = *w;      *w += p->ns;
  d->su = *w;      *w += p->ns;
  d->zsl = *w;     *w += p->ns;
  d->zsu = *w;     *w += p->ns;
  d->soft_ns = *iw;   *iw += p->ns>0 ? p->N+1 : 0;
  d->soft_offs = *iw; *iw += p->ns>0 ? p->N+2 : 0;
  d->soft_map = *iw;  *iw += p->ns;
  d->soft_local = *iw; *iw += p->ns;
  d->soft_idxs = *iw; *iw += p->ns>0 ? n_ineq_slots : 0;

  // The caller oracle's own buffers; aliases when nothing was rewritten, so
  // the solve loop uses them unconditionally
  if (p->rewrite.ne) {
    d->xu = *w;   *w += p->rewrite.nx;
    d->gu = *w;   *w += casadi_max(p->rewrite.nx, p->rewrite.ng);
    d->au = *w;   *w += p->rewrite.nnz_au;
    d->hu = *w;   *w += p->rewrite.nnz_hu;
    d->lamu = *w; *w += p->rewrite.ng;
  } else {
    d->xu = d->x;
    d->gu = d->g;
    d->au = d->a;
    d->hu = d->h;
    d->lamu = d->lam;
  }
  casadi_ipmc_rewrite_set_work(&d->rewrite, &p->rewrite, iw);

  d->pv = *w;

  d->arg = *arg;
  d->res = *res;
  d->iw = *iw;
  d->w = *w;
}

// SYMBOL "ipmc_presolve"
template<typename T1>
void casadi_ipmc_presolve(casadi_ipmc_data<T1>* d) {
  casadi_int k, i, start, stop, nx, j, col, n_soft, ns_k, n_ineq;
  casadi_int K, n_ineq_slots, column;
  ipmc_int *st_nu, *st_nx, *st_ng, *st_ngi, *st_ns, *st_idx;
  int soft;
  T1 lo, up;
  // Problem structure
  const casadi_ipmc_prob<T1>* p = d->prob;
  const casadi_nlpsol_prob<T1>* p_nlp = p->nlp;
  casadi_nlpsol_data<T1>* d_nlp = d->nlp;

  d->solver_created = 0;
  d->nxc_checked = 0;

  nx = p_nlp->nx;

  /* ---- native soft constraints, part 1: reset the per-solve maps -------
     (z and Z need nothing here: p->fs_z / p->fs_Z are constants evaluated
     once by IpmcInterface::init) */
  d->soft_error = 0;
  d->n_soft = 0;
  if (p->ns>0) {
    for (j=0;j<p->ns;++j) d->soft_local[j] = -1;
  }

  d->a_eq_idx[0] = 0;
  d->a_ineq_idx[0] = 0;
  d->x_eq_idx[0] = 0;
  d->x_ineq_idx[0] = 0;

  /* Loop over CD blocks.  A SOFTENED row is forced into the inequality set
     even when lbz==ubz: once relaxed it is genuinely two-sided
     (lbz - s_l <= z <= ubz + s_u), and ipmc can only soften inequality
     rows.  A softened row with both bounds infinite is still dropped: it is
     vacuous, relaxed or not. */
  for (k=0;k<p->N+1;++k) {
    d->a_eq_idx[k+1] = d->a_eq_idx[k];
    d->a_ineq_idx[k+1] = d->a_ineq_idx[k];
    start = p->CD[k].offset_r;
    stop  = p->CD[k].offset_r+p->CD[k].rows;
    for (i=start;i<stop;++i) {
      col = p->ns>0 ? p->slack_g[i] : -1;
      soft = col>=0;
      lo = d_nlp->lbz[nx+i];
      up = d_nlp->ubz[nx+i];
      if (lo==up && !soft) {
        d->a_eq[d->a_eq_idx[k+1]++] = i;
      } else {
        if (lo==-std::numeric_limits<T1>::infinity() && up==std::numeric_limits<T1>::infinity()) continue;
        d->a_ineq[d->a_ineq_idx[k+1]++] = i;
      }
    }
    d->x_eq_idx[k+1] = d->x_eq_idx[k];
    d->x_ineq_idx[k+1] = d->x_ineq_idx[k];
    start = p->CD[k].offset_c;
    stop  = p->CD[k].offset_c+p->CD[k].cols;

    for (i=start;i<stop;++i) {
      col = p->ns>0 ? p->slack_x[i] : -1;
      soft = col>=0;
      lo = d_nlp->lbz[i];
      up = d_nlp->ubz[i];
      if (lo==up && !soft) {
        d->x_eq[d->x_eq_idx[k+1]++] = i;
      } else {
        if (lo==-std::numeric_limits<T1>::infinity() && up==std::numeric_limits<T1>::infinity()) continue;
        d->x_ineq[d->x_ineq_idx[k+1]++] = i;
      }
    }
  }

  /* ---- native soft constraints, part 2 --------------------------------- */
  if (p->ns>0) {
    /* which columns are actually handed to ipmc?  ubs==0 pins the slack to
       zero, i.e. the row is hard; ipmc refuses a non-positive slack upper
       bound, so drop the column instead. */
    n_soft = 0;
    for (k=0;k<p->N+1;++k) {
      d->soft_offs[k] = n_soft;
      ns_k = 0;
      for (j=p->slack_offs[k];j<p->slack_offs[k+1];++j) {
        col = p->slack_perm[j];
        if (d->slack_ubs[col]<=0 && d->slack_ubs[p->ns+col]<=0) {
          d->soft_local[col] = -1;
          continue;
        }
        if (d->slack_ubs[col]<=0 || d->slack_ubs[p->ns+col]<=0) d->soft_error = 5;
        d->soft_local[col] = ns_k;
        d->soft_map[n_soft+ns_k] = col;
        ns_k++;
      }
      d->soft_ns[k] = ns_k;
      n_soft += ns_k;
    }
    d->soft_offs[p->N+1] = n_soft;
    d->n_soft = n_soft;

    /* idxs_rev, stage after stage, in the order the inequality rows are
       handed over */
    n_ineq = 0;
    for (k=0;k<p->N+1;++k) {
      for (i=d->a_ineq_idx[k];i<d->a_ineq_idx[k+1];++i) {
        col = p->slack_g[d->a_ineq[i]];
        d->soft_idxs[n_ineq++] = col<0 ? -1 : d->soft_local[col];
      }
      for (i=d->x_ineq_idx[k];i<d->x_ineq_idx[k+1];++i) {
        col = p->slack_x[d->x_ineq[i]];
        d->soft_idxs[n_ineq++] = col<0 ? -1 : d->soft_local[col];
      }
    }
  }

  if (d->soft_error) {
    /* nothing usable to hand over; IpmcInterface::solve reports it */
    if (d->solver) {
      ipmc_destroy(d->solver);
      d->solver = 0;
      d->layout = 0;
    }
    return;
  }

  /* ---- the problem STRUCTURE, in the shape ipmc_create takes it --------
     The only thing that can force the cached solver to be rebuilt: bounds,
     starting point, penalty and initial slacks are set on the live solver
     below and re-read at every ipmc_start. */
  K = p->N+1;
  n_ineq_slots = p_nlp->ng + p_nlp->nx;
  if (!d->structure) {
    /* one malloc for both halves; structure_built has to outlive the call,
       see casadi_ipmc_data.  5 per-stage rows, the row -> slack map, and
       one trailing scalar, the horizon-wide nxc. */
    d->sz_structure = 5*K + n_ineq_slots + 1;
    d->structure = (ipmc_int*) malloc(2*d->sz_structure*sizeof(ipmc_int));
    if (d->structure) {
      d->structure_built = d->structure + d->sz_structure;
      /* a value no real structure can take: the first compare always rebuilds */
      for (i=0;i<d->sz_structure;++i) d->structure_built[i] = -2;
    } else {
      d->sz_structure = 0;
    }
  }
  if (!d->structure) {
    d->soft_error = 7;
    if (d->solver) {
      ipmc_destroy(d->solver);
      d->solver = 0;
      d->layout = 0;
    }
    return;
  }
  st_nu = d->structure;
  st_nx = st_nu + K;
  st_ng = st_nx + K;
  st_ngi = st_ng + K;
  st_ns = st_ngi + K;
  st_idx = st_ns + K;
  n_ineq = 0;
  for (k=0;k<K;++k) {
    st_nu[k] = p->nu[k];
    st_nx[k] = p->nx[k];
    st_ng[k] = d->a_eq_idx[k+1] - d->a_eq_idx[k]
             + d->x_eq_idx[k+1] - d->x_eq_idx[k];
    st_ngi[k] = d->a_ineq_idx[k+1] - d->a_ineq_idx[k]
              + d->x_ineq_idx[k+1] - d->x_ineq_idx[k];
    st_ns[k] = p->ns>0 ? d->soft_ns[k] : 0;
    n_ineq += st_ngi[k];
  }
  for (i=0;i<n_ineq;++i) st_idx[i] = p->ns>0 ? d->soft_idxs[i] : -1;
  /* the horizon-wide nxc trails the map; n_ineq is pinned by the ng_ineq
     rows above, so a changed n_ineq is caught before this entry is ever
     read at a shifted position */
  d->structure[5*K+n_ineq] = p->nxc;

  if (d->solver) {
    int changed = 0;
    for (i=0;i<5*K+n_ineq+1;++i) {
      if (d->structure[i]!=d->structure_built[i]) { changed = 1; break; }
    }
    if (changed) {
      ipmc_destroy(d->solver);
      d->solver = 0;
      d->layout = 0;
    }
  }
  if (!d->solver) {
    struct IpmcProblem desc;
    IpmcError err = IPMC_OK;
    desc.K = K;
    desc.nu = st_nu;
    desc.nx = st_nx;
    desc.ng = st_ng;
    desc.ng_ineq = st_ngi;
    desc.nxc = p->nxc;
    desc.ns = p->ns>0 ? st_ns : 0;
    desc.soft_idx = p->ns>0 ? st_idx : 0;
    d->solver = ipmc_create(&desc, &err);
    if (!d->solver) {
      d->soft_error = err ? (int)err : -1;
      d->layout = 0;
      return;
    }
    d->solver_created = 1;
    d->layout = ipmc_get_layout(d->solver);
    ipmc_set_output(d->solver, p->write, p->flush);
    for (i=0;i<5*K+n_ineq+1;++i) d->structure_built[i] = d->structure[i];
  }

  /* ---- the NUMBERS.  Handed over on every solve, so a cached solver can
     never keep solving the previous call's bounds or starting point. ----- */
  column = 0;
  for (k=0;k<K;++k) {
    for (i=d->a_ineq_idx[k];i<d->a_ineq_idx[k+1];++i) {
      d->set_lower[column] = d_nlp->lbz[nx+d->a_ineq[i]];
      d->set_upper[column] = d_nlp->ubz[nx+d->a_ineq[i]];
      column++;
    }
    for (i=d->x_ineq_idx[k];i<d->x_ineq_idx[k+1];++i) {
      d->set_lower[column] = d_nlp->lbz[d->x_ineq[i]];
      d->set_upper[column] = d_nlp->ubz[d->x_ineq[i]];
      column++;
    }
  }
  ipmc_set_bounds(d->solver, d->set_lower, d->set_upper);

  /* the starting point, in ipmc's (u_0, x_0, u_1, x_1, ...) layout */
  column = 0;
  for (k=0;k<K;++k) {
    casadi_copy(d_nlp->z+p->CD[k].offset_c+p->nx[k], p->nu[k], d->set_ux0+column);
    column += p->nu[k];
    casadi_copy(d_nlp->z+p->CD[k].offset_c, p->nx[k], d->set_ux0+column);
    column += p->nx[k];
  }
  ipmc_set_initial(d->solver, d->set_ux0);

  /* the penalty: ipmc adds zl*sl + 1/2 Zl sl^2 + zu*su + 1/2 Zu su^2 to the
     objective -- the Taylor expansion of the (validated) separable
     quadratic f_s around s=0 */
  if (d->n_soft>0) {
    IpmcError err;
    T1* Zl = d->set_pen;
    T1* Zu = Zl + d->n_soft;
    T1* zl = Zu + d->n_soft;
    T1* zu = zl + d->n_soft;
    T1* ubsl = zu + d->n_soft;
    T1* ubsu = ubsl + d->n_soft;
    for (k=0;k<K;++k) {
      for (j=0;j<d->soft_ns[k];++j) {
        i = d->soft_offs[k]+j;
        col = d->soft_map[i];
        Zl[i] = p->fs_Z[col];
        Zu[i] = p->fs_Z[p->ns+col];
        zl[i] = p->fs_z[col];
        zu[i] = p->fs_z[p->ns+col];
        ubsl[i] = d->slack_ubs[col];
        ubsu[i] = d->slack_ubs[p->ns+col];
      }
    }
    err = ipmc_set_soft_penalty(d->solver, Zl, Zu, zl, zu, ubsl, ubsu);
    if (err!=IPMC_OK) {
      d->soft_error = (int)err;
      ipmc_destroy(d->solver);
      d->solver = 0;
      d->layout = 0;
      return;
    }
    /* the user's s0 */
    if (d->slack_s) {
      for (k=0;k<K;++k) {
        for (j=0;j<d->soft_ns[k];++j) {
          i = d->soft_offs[k]+j;
          col = d->soft_map[i];
          d->set_sl0[i] = d->slack_s[col];
          d->set_su0[i] = d->slack_s[p->ns+col];
        }
      }
      ipmc_set_initial_slack(d->solver, d->set_sl0, d->set_su0);
    }
  }
}

// SYMBOL "ipmc_read_dual"
/* ipmc's multipliers, laid out stage after stage, back into the lam vector
   of the problem that was HANDED OVER -- the rewritten one when a slack
   column was lifted, and then casadi_ipmc_rewrite_collect takes it the rest
   of the way.  Used when the solve is over and, when a progress hook is
   installed, at every iteration boundary. */
template<typename T1>
void casadi_ipmc_read_dual(casadi_ipmc_data<T1>* d, const struct IpmcLayout* str,
                              const double* dual_data) {
  casadi_int k, i, column;
  const casadi_ipmc_prob<T1>* p = d->prob;
  const casadi_nlpsol_prob<T1>* p_nlp = p->nlp;
  casadi_nlpsol_data<T1>* d_nlp = d->nlp;
  // Inequalities
  for (k=0;k<str->K;++k) {
    column = 0;
    for (i=d->a_ineq_idx[k];i<d->a_ineq_idx[k+1];++i) {
      d_nlp->lam[p_nlp->nx+d->a_ineq[i]] = dual_data[str->g_ineq_offs[k]+column];
      column++;
    }
    for (i=d->x_ineq_idx[k];i<d->x_ineq_idx[k+1];++i) {
      d_nlp->lam[d->x_ineq[i]] = dual_data[str->g_ineq_offs[k]+column];
      column++;
    }
  }
  // Equalities
  for (k=0;k<str->K;++k) {
    column = 0;
    for (i=d->a_eq_idx[k];i<d->a_eq_idx[k+1];++i) {
      d_nlp->lam[p_nlp->nx+d->a_eq[i]] = dual_data[str->g_offs[k]+column];
      column++;
    }
    for (i=d->x_eq_idx[k];i<d->x_eq_idx[k+1];++i) {
      d_nlp->lam[d->x_eq[i]] = dual_data[str->g_offs[k]+column];
      column++;
    }
  }
  // Dynamics
  for (k=0;k<str->K-1;++k) {
    casadi_scaled_copy(-1.0, dual_data+str->dyn_eq_offs[k], p->nx[k+1], d_nlp->lam+p_nlp->nx+p->AB[k].offset_r);
  }
}

// SYMBOL "ipmc_lift_penalty"
/* The penalty the rewritten objective puts on the lifted helper columns:
   sum_e (z m_0 + 1/2 Z m_0^2), the same formula casadi_ipmc_lift_obj adds
   at evaluation time, here read at the reported iterate in d->nlp->z.  The
   caller's objective is the one WITHOUT that penalty, so this comes back
   off whatever is reported. */
template<typename T1>
T1 casadi_ipmc_lift_penalty(casadi_ipmc_data<T1>* d) {
  casadi_int e, slot;
  T1 m0, pen;
  const casadi_ipmc_prob<T1>* p = d->prob;
  pen = 0;
  for (e=0;e<p->rewrite.ne;++e) {
    slot = p->rewrite.slot[e];
    m0 = d->nlp->z[p->rewrite.m[e]];
    pen += p->fs_z[slot]*m0 + 0.5*p->fs_Z[slot]*m0*m0;
  }
  return pen;
}

// SYMBOL "ipmc_report_iterate"
/* One iterate, in the CALLER's z-space, ready for Nlpsol::callback.
   Everything ipmc holds is in the space of the problem that was handed
   over, the REWRITTEN one whenever a slack column was lifted, so all four
   of x, g, lam and f make the same trip the end-of-solve read-back makes.
   Order matters: the rewritten iterate has to be complete before _collect /
   _collect_g gather out of it, and the penalty is read off the helper
   variables it puts there.  When nothing was lifted every step is a no-op
   or writes the caller's own vector in place. */
template<typename T1>
void casadi_ipmc_report_iterate(casadi_ipmc_data<T1>* d, const struct IpmcLayout* str,
    const double* primal_data, const double* dual_data, T1 f) {
  const casadi_ipmc_prob<T1>* p = d->prob;
  const casadi_nlpsol_prob<T1>* p_nlp = p->nlp;
  casadi_nlpsol_data<T1>* d_nlp = d->nlp;
  casadi_ipmc_read_primal_data(p, primal_data, d_nlp->z, str);
  casadi_copy(d->cb_g, p_nlp->ng, d_nlp->z+p_nlp->nx);
  if (dual_data) casadi_ipmc_read_dual(d, str, dual_data);
  d->rewrite.user->objective = f - casadi_ipmc_lift_penalty(d);
  casadi_ipmc_rewrite_collect(&p->rewrite, &d->rewrite, d_nlp,
                                d->slack_s, d->slack_lam_s);
  casadi_ipmc_rewrite_collect_g(&p->rewrite, &d->rewrite, d_nlp);
}

// SYMBOL "ipmc_finish"
/* Everything a solve does once its loop is over: status, statistics, the
   primal and dual read-back, and the slacks.  No oracle call is reachable
   from here, which is why it is ONE routine while the loop driving ipmc is
   written out twice (C++ in IpmcInterface::solve, emitted C in its
   codegen_body); each calls it where the loop it just ran ends. */
template<typename T1>
void casadi_ipmc_finish(casadi_ipmc_data<T1>* d) {
  casadi_int i, col;
  const casadi_ipmc_prob<T1>* p = d->prob;
  casadi_nlpsol_data<T1>* d_nlp = d->nlp;
  const struct IpmcLayout* str = d->layout;
  ipmc_int ret;

  ret = ipmc_get_status(d->solver);

  d->return_status = ret;
  if (ret==IPMC_SOLVED) {
    d->unified_return_status = SOLVER_RET_SUCCESS;
    d->success = 1;
  }

  if (ret==IPMC_INFEASIBLE) {
    /* the restoration phase converged to a point still infeasible for the
       original problem -- ipopt's Restoration_Failed.  NOT a solution. */
    d->unified_return_status = SOLVER_RET_INFEASIBLE;
  }

  const double* primal_data = ipmc_get_primal(d->solver);
  const double* dual_data = ipmc_get_dual(d->solver);
  const struct IpmcStats* stats = ipmc_get_stats(d->solver);

  d->stats.compute_sd_time = stats->compute_sd_time;
  d->stats.duinf_time = stats->duinf_time;
  d->stats.eval_hess_time = stats->eval_hess_time;
  d->stats.eval_jac_time = stats->eval_jac_time;
  d->stats.eval_cv_time = stats->eval_cv_time;
  d->stats.eval_grad_time = stats->eval_grad_time;
  d->stats.eval_obj_time = stats->eval_obj_time;
  d->stats.initialization_time = stats->initialization_time;
  d->stats.time_total = stats->time_total;
  d->stats.eval_hess_count = stats->eval_hess_count;
  d->stats.eval_jac_count = stats->eval_jac_count;
  d->stats.eval_cv_count = stats->eval_cv_count;
  d->stats.eval_grad_count = stats->eval_grad_count;
  d->stats.eval_obj_count = stats->eval_obj_count;
  d->stats.iterations_count = stats->iterations_count;
  d->stats.restoration_iterations_count = stats->restoration_iterations_count;
  d->stats.return_flag = stats->return_flag;

  casadi_ipmc_read_primal_data(p, primal_data, d_nlp->z, str);
  casadi_ipmc_read_dual(d, str, dual_data);

  /* Slacks.  ipmc lays them out stage after stage; soft_map takes us back
     to casadi's own column order.  A column that was dropped (ubs==0 on
     both sides) keeps s=0, the only value it could have taken anyway.
     Sign: ipmc returns zsl>=0, the multiplier of sl>=0, and casadi's
     convention is lam<0 at an active LOWER bound, so lam_s = -zs. */
  if (p->ns>0) {
    for (i=0;i<2*p->ns;++i) {
      d->slack_s[i] = 0;
      d->slack_lam_s[i] = 0;
    }
    if (d->n_soft>0) {
      ipmc_get_slack(d->solver, d->sl, d->su);
      ipmc_get_slack_dual(d->solver, d->zsl, d->zsu);
      for (i=0;i<d->n_soft;++i) {
        col = d->soft_map[i];
        d->slack_s[col] = d->sl[i];
        d->slack_s[p->ns+col] = d->su[i];
        d->slack_lam_s[col] = -d->zsl[i];
        d->slack_lam_s[p->ns+col] = -d->zsu[i];
      }
    }
  }
}

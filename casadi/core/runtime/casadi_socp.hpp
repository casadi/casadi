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

// Reusable runtime building block that converts CasADi's Q/P input pair
// (SDP-style cone description, see Conic::sdp_to_socp_init) into the
// per-cone SOCP material that solver plugins need:
//
//   * lifted helper variables [X_1; Z_1; X_2; Z_2; ...] with bounds
//     X in (-inf, +inf), Z in [0, +inf)
//
//   * a block of linear equality constraints  Q*[x; lifted] = -P,
//     materialized as CSR-style arrays (start/colind/coef + rowtype/rhs/rng)
//
//   * a block of "<= 0" rows, one per cone, that the consumer turns
//     quadratic via its native API (e.g. XPRSaddqmatrix, GRBaddqconstr).
//     The cone's diagonal quadratic triplet (1, 1, ..., 1, -1) is
//     filled into a per-cone scratch by casadi_socp_cone_build.
//
// Typical consumer plugin flow (in solve):
//   d->q = arg[CONIC_Q]; d->p = arg[CONIC_P];
//   casadi_socp_build(d);
//   solver_add_cols(n_lifted, d->obj_lift, d->lift_start,
//                   d->lb_lift, d->ub_lift);
//   solver_add_rows(p->n_eq, p->eq_nnz, d->eq_type, d->eq_rhs, d->eq_rng,
//                   d->eq_start, d->eq_colind, d->eq_coef);
//   solver_add_rows(p->n_blocks, 0, d->cone_type, d->cone_rhs, d->cone_rng,
//                   d->cone_start, NULL, NULL);
//   for (b = 0; b < p->n_blocks; ++b) {
//     casadi_int bs = casadi_socp_cone_build(d, b);
//     solver_attach_qmatrix(cone_row_offset + b, bs,
//                           d->qcol1, d->qcol2, d->qcoef);
//   }

// SYMBOL "socp_prob"
template<typename T1>
struct casadi_socp_prob {
  // Number of cone blocks; 0 if no SOC constraints (entire component inert).
  casadi_int n_blocks;
  // Block boundaries in lifted space, length n_blocks+1.  Block b owns
  // lifted indices r[b]..r[b+1]-1.  The last entry of each block is the
  // cone's "Z" scalar; preceding entries are the cone's "X" vector.
  const casadi_int *r;
  // Total lifted variable count = r[n_blocks]
  casadi_int n_lifted;
  // Number of original (un-lifted) decision variables.  Used to
  // distinguish original-variable rows from lifted-variable rows in
  // map_Q.
  casadi_int nx;
  // map_Q in CSC form.  Rows = (nx + n_lifted), cols = n_lifted.
  // Each column is one equality constraint.
  const casadi_int *mq_colind;   // length n_lifted+1
  const casadi_int *mq_row;      // length eq_nnz
  // For nz k, if mq_row[k] < nx the entry's coefficient comes from
  // q[mq_data[k]] (an index into the user-supplied Q value array);
  // otherwise the coefficient is -1 (the -I block).
  const casadi_int *mq_data;     // length eq_nnz
  // map_P, length n_lifted.  rhs of equality row i is
  //   -p[map_P[i]]  if map_P[i] != -1
  //   0             otherwise
  const casadi_int *map_P;

  // Derived (filled by casadi_socp_setup)
  casadi_int n_eq;
  casadi_int eq_nnz;
  casadi_int max_block;
};
// C-REPLACE "casadi_socp_prob<T1>" "struct casadi_socp_prob"
// C-REPLACE "reinterpret_cast<int*>" "(int*) "
// C-REPLACE "reinterpret_cast<char*>" "(char*) "
// C-REPLACE "static_cast<int>" "(int) "
// C-REPLACE "std::numeric_limits<T1>::infinity()" "casadi_inf"

// SYMBOL "socp_setup"
template<typename T1>
void casadi_socp_setup(casadi_socp_prob<T1>* p) {
  casadi_int b;
  if (p->n_blocks == 0) {
    p->n_eq = 0;
    p->eq_nnz = 0;
    p->max_block = 0;
    return;
  }
  p->n_eq = p->n_lifted;
  p->eq_nnz = p->mq_colind[p->n_eq];
  p->max_block = 0;
  for (b = 0; b < p->n_blocks; ++b) {
    casadi_int bs = p->r[b + 1] - p->r[b];
    if (bs > p->max_block) p->max_block = bs;
  }
}

// SYMBOL "socp_data"
template<typename T1>
struct casadi_socp_data {
  const casadi_socp_prob<T1>* prob;

  // Inputs (pointed at user-supplied Q and P value arrays).  May be NULL,
  // in which case the corresponding entries are filled with 0.
  const T1 *q;
  const T1 *p;

  // Workspace populated in casadi_socp_init.  All sized to fit one
  // "build" pass; callers must not retain pointers across builds.

  // Lifted-variable column data (n_lifted entries each)
  T1 *obj_lift;
  T1 *lb_lift;
  T1 *ub_lift;
  int *lift_start;       // n_lifted+1, all zero (no entries in existing rows)

  // Linear equality rows (n_eq rows, eq_nnz nonzeros)
  int *eq_start;         // n_eq+1
  int *eq_colind;        // eq_nnz
  T1 *eq_coef;           // eq_nnz
  T1 *eq_rhs;            // n_eq
  T1 *eq_rng;            // n_eq, all zero
  char *eq_type;         // n_eq, all 'E'

  // Per-cone "<= 0" placeholder rows (n_blocks rows, no linear nonzeros)
  int *cone_start;       // n_blocks+1, all zero
  T1 *cone_rhs;          // n_blocks, all zero
  T1 *cone_rng;          // n_blocks, all zero
  char *cone_type;       // n_blocks, all 'L'

  // Per-cone qmatrix scratch (max_block entries each).  Re-filled per
  // cone by casadi_socp_cone_build.
  int *qcol1;
  int *qcol2;
  T1 *qcoef;
};
// C-REPLACE "casadi_socp_data<T1>" "struct casadi_socp_data"

// SYMBOL "socp_work"
template<typename T1>
void casadi_socp_work(const casadi_socp_prob<T1>* p,
    casadi_int* sz_iw, casadi_int* sz_w) {
  if (p->n_blocks == 0) return;
  // doubles
  *sz_w  += 3 * p->n_lifted;       // obj_lift, lb_lift, ub_lift
  *sz_w  += p->eq_nnz;             // eq_coef
  *sz_w  += 2 * p->n_eq;           // eq_rhs, eq_rng
  *sz_w  += 2 * p->n_blocks;       // cone_rhs, cone_rng
  *sz_w  += p->max_block;          // qcoef
  // ints / chars (overlaid on iw slots)
  *sz_iw += p->n_lifted + 1;       // lift_start
  *sz_iw += p->n_eq + 1;           // eq_start
  *sz_iw += p->eq_nnz;             // eq_colind
  *sz_iw += p->n_eq;               // eq_type
  *sz_iw += p->n_blocks + 1;       // cone_start
  *sz_iw += p->n_blocks;           // cone_type
  *sz_iw += 2 * p->max_block;      // qcol1, qcol2
}

// SYMBOL "socp_init"
template<typename T1>
void casadi_socp_init(casadi_socp_data<T1>* d, casadi_int** iw, T1** w) {
  const casadi_socp_prob<T1>* p = d->prob;
  if (p->n_blocks == 0) return;
  d->obj_lift   = *w;  *w += p->n_lifted;
  d->lb_lift    = *w;  *w += p->n_lifted;
  d->ub_lift    = *w;  *w += p->n_lifted;
  d->eq_coef    = *w;  *w += p->eq_nnz;
  d->eq_rhs     = *w;  *w += p->n_eq;
  d->eq_rng     = *w;  *w += p->n_eq;
  d->cone_rhs   = *w;  *w += p->n_blocks;
  d->cone_rng   = *w;  *w += p->n_blocks;
  d->qcoef      = *w;  *w += p->max_block;
  d->lift_start = reinterpret_cast<int*>(*iw);  *iw += p->n_lifted + 1;
  d->eq_start   = reinterpret_cast<int*>(*iw);  *iw += p->n_eq + 1;
  d->eq_colind  = reinterpret_cast<int*>(*iw);  *iw += p->eq_nnz;
  d->eq_type    = reinterpret_cast<char*>(*iw); *iw += p->n_eq;
  d->cone_start = reinterpret_cast<int*>(*iw);  *iw += p->n_blocks + 1;
  d->cone_type  = reinterpret_cast<char*>(*iw); *iw += p->n_blocks;
  d->qcol1      = reinterpret_cast<int*>(*iw);  *iw += p->max_block;
  d->qcol2      = reinterpret_cast<int*>(*iw);  *iw += p->max_block;
}

// SYMBOL "socp_build"
// Populate all lifted-variable + equality-row + cone-row arrays from
// d->q / d->p.  After this returns, the consumer can hand the arrays
// straight to its solver's add-cols / add-rows APIs.
//
// The per-cone qmatrix scratch (d->qcol1, d->qcol2, d->qcoef) is NOT
// filled by this function -- call casadi_socp_cone_build per block.
// C-REPLACE "std::numeric_limits<T1>::infinity()" "casadi_inf"
template<typename T1>
void casadi_socp_build(casadi_socp_data<T1>* d) {
  const casadi_socp_prob<T1>* p = d->prob;
  casadi_int b, j, k, kk, kbeg, kend;
  if (p->n_blocks == 0) return;

  // Lifted column metadata
  for (j = 0; j < p->n_lifted; ++j) {
    d->obj_lift[j] = 0;
    d->ub_lift[j] = std::numeric_limits<T1>::infinity();
  }
  for (j = 0; j < p->n_lifted + 1; ++j) d->lift_start[j] = 0;
  for (b = 0; b < p->n_blocks; ++b) {
    casadi_int bs = p->r[b + 1] - p->r[b];
    for (j = 0; j < bs - 1; ++j) {
      d->lb_lift[p->r[b] + j] = -std::numeric_limits<T1>::infinity();
    }
    d->lb_lift[p->r[b] + bs - 1] = 0;
  }

  // Linear equality rows
  kk = 0;
  for (j = 0; j < p->n_eq; ++j) {
    d->eq_start[j] = static_cast<int>(kk);
    d->eq_type[j] = 'E';
    casadi_int idx = p->map_P[j];
    d->eq_rhs[j] = (d->p && idx >= 0) ? -d->p[idx] : 0;
    d->eq_rng[j] = 0;
    kbeg = p->mq_colind[j];
    kend = p->mq_colind[j + 1];
    for (k = kbeg; k < kend; ++k) {
      casadi_int row = p->mq_row[k];
      d->eq_colind[kk] = static_cast<int>(row);
      d->eq_coef[kk] = (d->q && row < p->nx) ? d->q[p->mq_data[k]] : -1;
      kk++;
    }
  }
  d->eq_start[p->n_eq] = static_cast<int>(kk);

  // Per-cone "<= 0" placeholder rows (no linear coefficients).
  for (b = 0; b < p->n_blocks; ++b) {
    d->cone_type[b] = 'L';
    d->cone_rhs[b] = 0;
    d->cone_rng[b] = 0;
    d->cone_start[b] = 0;
  }
  d->cone_start[p->n_blocks] = 0;
}

// SYMBOL "socp_cone_build"
// Fill the per-cone qmatrix scratch (d->qcol1, d->qcol2, d->qcoef) for
// cone block b.  Returns block_size.
//
// The constraint is  X_1^2 + ... + X_{bs-1}^2 - Z^2 <= 0,
// expressed as the diagonal quadratic [1, 1, ..., 1, -1] over the
// lifted variables.
template<typename T1>
casadi_int casadi_socp_cone_build(casadi_socp_data<T1>* d, casadi_int b) {
  const casadi_socp_prob<T1>* p = d->prob;
  casadi_int bs = p->r[b + 1] - p->r[b];
  casadi_int j;
  for (j = 0; j < bs; ++j) {
    int col = static_cast<int>(p->nx + p->r[b] + j);
    d->qcol1[j] = col;
    d->qcol2[j] = col;
    d->qcoef[j] = (j == bs - 1) ? -1 : 1;
  }
  return bs;
}

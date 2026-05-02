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

// casadi_blas_mtimes is the BLAS-dispatch wrapper around casadi_mtimes_dense.
// vm: C++ template in blas.hpp; <double> specialization in blas.cpp routes
//     through Blas::mtimes(default_, ...).
// codegen: emitted by AUX_BLAS_MTIMES as a thin static wrapper around
//     casadi_mtimes_dense (which the plugin AUX hook may override with
//     cblas_dgemm).  Same call sites in source and generated C.

// Plain (M-step) condensing of OCP-structured QPs.
//
// Given an OCP with horizon N, stage dimensions nx[k], nu[k] (with nu[N]=0),
// path-constraint dimensions ng[k], and a partition M[0]=0 < M[1] < ... <
// M[N_hat] = N, condense each block [M[K], M[K+1]) into a single condensed
// stage K by eliminating the inner states.  The output is again an OCP, with
// the same data layout convention, with horizon N_hat and:
//
//   nx_hat[K]  = nx[M[K]]                              (only "boundary" state)
//   nu_hat[K]  = sum_{j in [M[K], M[K+1])} nu[j]       (concatenated controls)
//   ng_hat[K]  = sum_{j in [M[K], M[K+1])} ng[j]
//              + sum_{j in [M[K]+1, M[K+1])} nx[j]     (lifted state bounds)
//   nx_hat[N_hat] = nx[N], nu_hat[N_hat] = 0,
//   ng_hat[N_hat] = ng[N]
//
// Layout convention (must match the caller's producer):
//   AB block:  rows = nx[k+1], cols = nx[k] + nu[k];
//              column-major dense; first nx[k] columns are A_k, the next
//              nu[k] columns are B_k.
//   RSQ block: rows = cols = nx[k] + nu[k]; symmetric Hessian laid out
//              column-major as [Qxx Qxu; Qxu' Quu] -- i.e. variable order
//              is [x_k; u_k] within a stage.
//   CD block:  rows = ng[k], cols = nx[k] + nu[k];
//              column-major; first nx[k] columns are C_k, the next nu[k]
//              columns are D_k.
//   Stage gradient qr[k]: length nx[k]+nu[k], packed [q_k; r_k].
//   Affine offset b[k]:   length nx[k+1] (dynamics RHS).
//   Bounds lbx[k], ubx[k]: length nx[k]; lbu[k], ubu[k]: length nu[k];
//   lbg[k], ubg[k]:      length ng[k].
//
// The caller stores per-block data flat: AB_val, CD_val, RSQ_val are flat
// arrays indexed via AB_offsets[K] etc. (matching the fatrop interface
// convention).  Per-stage vectors (b, qr, lbx, ubx, lbu, ubu, lbg, ubg) are
// concatenated stage-by-stage; the helper computes offsets internally from
// nx, nu, ng.
//
// Algorithm (per condensed block K spanning original stages [k_a, k_b)):
//   Initialise Phi = I_{nx[k_a]}, Gamma = (zero columns), phi = 0
//   Initialise Hxx = 0, Hxu = 0, Huu = 0, hx = 0, hu = 0
//   For j = 0, 1, ..., M-1 with k = k_a + j:
//     read R_tilde = RSQ[k] = [Qxx Qxu; Qxu' Quu], qr = [q_k; r_k]
//     read A_k, B_k (from AB[k]), b_k (affine), lbu_k, ubu_k, lbg_k, ubg_k
//
//     Add stage-cost contribution to (Hxx, Hxu, Huu, hx, hu) using the
//     substitution x_k = Phi*xi + Gamma*nu_<j + phi.
//
//     If j > 0 lift the state bound on x_{k_a+j} to a path inequality on
//     (xi, nu_<j) and append to CD_hat[K].  Append the original path
//     inequality C_k x_k + D_k u_k <= d_k similarly (column block for u_k
//     gets D_k, others zero).
//
//     Append lbu_k, ubu_k to lbu_hat[K], ubu_hat[K].
//
//     Advance:
//       phi   <- A_k * phi + b_k
//       Gamma <- [A_k * Gamma  B_k]    (one new column block of width nu[k])
//       Phi   <- A_k * Phi
//   End for.
//
//   Output for block K:
//     AB_hat[K]   = [Phi_M | Gamma_M]
//     b_hat[K]    = phi_M
//     RSQ_hat[K]  = [Hxx Hxu; Hxu' Huu]
//     qr_hat[K]   = [hx; hu]
//     CD_hat[K]   = stacked path-ineq + lifted-state-bound rows
//     lbx_hat[K]  = lbx[k_a]; ubx_hat[K] = ubx[k_a]   (boundary state)
//
// Complexity per condensed block of length M, with nx_hat = O(n_x) and
// per-stage nu = O(n_u): O(M * n_x^3 + M^2 * n_x^2 * n_u + M^2 * n_x * n_u^2).
// Summed over all blocks (sum M = N): O(N * M_avg * n_x^2 * (n_x + n_u)) when
// M_avg = max M.  Linear in N once M is fixed; one M factor better than the
// naive "rebuild Phi, Gamma from scratch per j" formulation.

// C-REPLACE "casadi_ocp_block" "struct casadi_ocp_block"
// C-REPLACE "casadi_condensing_prob<T1>" "struct casadi_condensing_prob"
// C-REPLACE "casadi_condensing_data<T1>" "struct casadi_condensing_data"
// C-REPLACE "static_cast<int>" "(int) "
// C-REPLACE "reinterpret_cast<casadi_ocp_block*>" "(casadi_ocp_block*) "
// C-REPLACE "std::numeric_limits<T1>::infinity()" "casadi_inf"

// SYMBOL "condensing_prob"
template<typename T1>
struct casadi_condensing_prob {
  // Original OCP descriptors -- caller-owned, just pointers
  const casadi_int *nx, *nu, *ng;     // length N+1
  casadi_int N;
  const casadi_ocp_block *AB, *CD, *RSQ;
  // Flat-storage offsets (length N+1; AB has length N -- last entry sentinel)
  const casadi_int *AB_offsets, *CD_offsets, *RSQ_offsets;

  // Partition: condensed stage K covers original stages [M[K], M[K+1])
  const casadi_int *M;
  casadi_int N_hat;

  // Note: per-stage derived arrays (nx_hat[], nu_hat[], ng_hat[],
  // AB_hat[], CD_hat[], RSQ_hat[], *_hat_offsets[]) live on
  // casadi_condensing_data and are populated by casadi_condensing_set_work.
  // Only *scalar* derived metadata stays on prob (below).

  // Maxima for workspace sizing (filled by condensing_setup)
  casadi_int nx_max;          // max over k of nx[k]
  casadi_int nu_max_block;    // max over K of nu_hat[K]
  casadi_int nxu_max_block;   // max over K of (nx_hat[K] + nu_hat[K])

  // Cumulative sizes (filled by condensing_setup) -- used by
  // condensing_work / condensing_init to claim all per-stage flat buffers
  // and condensed-CSC scratch from the iw / w workspace.
  casadi_int nnz_RSQ;          // sum_k (nx[k]+nu[k])^2          [== RSQ_offsets[N+1]]
  casadi_int nnz_AB;           // sum_{k<N} nx[k+1]*(nx[k]+nu[k]) [== AB_offsets[N]]
  casadi_int nnz_CD;           // sum_k ng[k]*(nx[k]+nu[k])      [== CD_offsets[N+1]]
  casadi_int total_x;          // sum_k nx[k]
  casadi_int total_u;          // sum_k nu[k]
  casadi_int total_g;          // sum_k ng[k]
  casadi_int total_b;          // sum_{k<N} nx[k+1]
  casadi_int total_qr;         // sum_k (nx[k]+nu[k])  (== nx_total)
  casadi_int nnz_RSQ_hat;
  casadi_int nnz_AB_hat;
  casadi_int nnz_CD_hat;
  casadi_int total_x_hat, total_u_hat, total_g_hat, total_b_hat;
  casadi_int nx_total_hat;     // sum_K (nx_hat[K]+nu_hat[K])
  casadi_int na_total_hat;     // sum_{K<N_hat} nx_hat[K+1] + sum_K ng_hat[K]
  casadi_int nnz_a_hat_csc;    // CSC nnz of the condensed A matrix
};

// SYMBOL "condensing_setup"
//   Computes the *scalar* derived metadata on prob (cumulative sizes,
//   maxima) used to size workspace.  All per-stage derived arrays
//   (nx_hat[], nu_hat[], ng_hat[], hat block descriptors and
//   *_hat_offsets[]) live on data and are populated by
//   casadi_condensing_set_work from the same recurrence.  The two stay
//   consistent because both walk prob->nx/nu/ng/M.
template<typename T1>
void casadi_condensing_setup(casadi_condensing_prob<T1>* p) {
  casadi_int K, j, k_a, k_b, M_, sum_nu, sum_ng_path, sum_ng_lift;
  casadi_int nx_K, nxu_K;
  casadi_int off_AB = 0, off_CD = 0, off_RSQ = 0;

  p->nx_max = 0;
  p->nu_max_block = 0;
  p->nxu_max_block = 0;

  // Pass 1: maxima, hat dimensions (locally), hat offsets (locally).
  // Only writes scalars onto prob.
  for (K = 0; K <= p->N_hat; ++K) {
    k_a = p->M[K];
    if (p->nx[k_a] > p->nx_max) p->nx_max = p->nx[k_a];
  }
  casadi_int nnz_RSQ_hat = 0, nnz_AB_hat = 0, nnz_CD_hat = 0;
  casadi_int tot_x_hat = 0, tot_u_hat = 0, tot_g_hat = 0, tot_b_hat = 0;
  casadi_int nx_total_hat = 0, na_total_hat = 0;
  casadi_int nnz_a = 0;
  casadi_int nu_hat_K, ng_hat_K, nx_hat_K, nx_hat_Kp1;
  for (K = 0; K < p->N_hat; ++K) {
    k_a = p->M[K];
    k_b = p->M[K + 1];
    M_ = k_b - k_a;
    sum_nu = 0;
    sum_ng_path = 0;
    sum_ng_lift = 0;
    for (j = 0; j < M_; ++j) {
      sum_nu += p->nu[k_a + j];
      sum_ng_path += p->ng[k_a + j];
      if (j > 0) sum_ng_lift += p->nx[k_a + j];
    }
    nu_hat_K = sum_nu;
    ng_hat_K = sum_ng_path + sum_ng_lift;
    nx_hat_K = p->nx[k_a];
    nx_hat_Kp1 = p->nx[p->M[K + 1]];
    if (nu_hat_K > p->nu_max_block) p->nu_max_block = nu_hat_K;
    nxu_K = nx_hat_K + nu_hat_K;
    if (nxu_K > p->nxu_max_block) p->nxu_max_block = nxu_K;

    off_AB  += nx_hat_Kp1 * nxu_K;
    off_CD  += ng_hat_K * nxu_K;
    off_RSQ += nxu_K * nxu_K;

    tot_x_hat += nx_hat_K;
    tot_u_hat += nu_hat_K;
    tot_g_hat += ng_hat_K;
    tot_b_hat += nx_hat_Kp1;
    nx_total_hat += nxu_K;
    na_total_hat += ng_hat_K + nx_hat_Kp1;
    if (K >= 1) nnz_a += nx_hat_K;
    nnz_a += nxu_K * nx_hat_Kp1;
    nnz_a += nxu_K * ng_hat_K;
  }
  // Terminal stage
  nx_K = p->nx[p->M[p->N_hat]];
  ng_hat_K = p->ng[p->N];
  if (nx_K > p->nxu_max_block) p->nxu_max_block = nx_K;
  off_RSQ += nx_K * nx_K;
  nnz_CD_hat = off_CD + ng_hat_K * nx_K;
  tot_x_hat += nx_K;
  tot_g_hat += ng_hat_K;
  nx_total_hat += nx_K;
  na_total_hat += ng_hat_K;
  if (p->N_hat >= 1) nnz_a += nx_K;
  nnz_a += nx_K * ng_hat_K;
  nnz_RSQ_hat = off_RSQ;
  nnz_AB_hat = off_AB;

  // Cumulative sizes over the *original* horizon (don't depend on M)
  {
    casadi_int kk;
    casadi_int nnz_RSQ = 0, nnz_AB = 0, nnz_CD = 0;
    casadi_int tot_x = 0, tot_u = 0, tot_g = 0, tot_b = 0;
    for (kk = 0; kk <= p->N; ++kk) {
      casadi_int nxu_k = p->nx[kk] + p->nu[kk];
      nnz_RSQ += nxu_k * nxu_k;
      nnz_CD  += p->ng[kk] * nxu_k;
      tot_x   += p->nx[kk];
      tot_u   += p->nu[kk];
      tot_g   += p->ng[kk];
      if (kk < p->N) {
        nnz_AB += p->nx[kk + 1] * nxu_k;
        tot_b  += p->nx[kk + 1];
      }
    }
    p->nnz_RSQ = nnz_RSQ;
    p->nnz_AB = nnz_AB;
    p->nnz_CD = nnz_CD;
    p->total_x = tot_x;
    p->total_u = tot_u;
    p->total_g = tot_g;
    p->total_b = tot_b;
    p->total_qr = tot_x + tot_u;
  }

  p->nnz_RSQ_hat = nnz_RSQ_hat;
  p->nnz_AB_hat = nnz_AB_hat;
  p->nnz_CD_hat = nnz_CD_hat;
  p->total_x_hat = tot_x_hat;
  p->total_u_hat = tot_u_hat;
  p->total_g_hat = tot_g_hat;
  p->total_b_hat = tot_b_hat;
  p->nx_total_hat = nx_total_hat;
  p->na_total_hat = na_total_hat;
  p->nnz_a_hat_csc = nnz_a;
}

// SYMBOL "condensing_data"
//
// All of the per-stage flat input/output buffers and the condensed-CSC +
// bound + primal/dual scratch are claimed from `w` by condensing_init.
// The caller writes into the input pointers (AB_val, RSQ_val, ...) via
// the extract helpers (or casadi_project for sparse projection); the
// solver writes into the condensed primal/dual fields; condensing_lift
// reads them back.
template<typename T1>
struct casadi_condensing_data {
  const casadi_condensing_prob<T1> *prob;

  // Per-condensed-stage derived dimensions and block descriptors.
  // Length N_hat+1.  Claimed from iw and populated in
  // casadi_condensing_set_work (reads prob->nx, ->nu, ->ng, ->M).
  casadi_int *nx_hat, *nu_hat, *ng_hat;
  casadi_ocp_block *AB_hat, *CD_hat, *RSQ_hat;
  casadi_int *AB_hat_offsets, *CD_hat_offsets, *RSQ_hat_offsets;

  // Per-stage flat input buffers (filled by caller before eval)
  T1 *AB_val;        // length p->nnz_AB
  T1 *CD_val;        // length p->nnz_CD
  T1 *RSQ_val;       // length p->nnz_RSQ
  T1 *b_val;         // length p->total_b
  T1 *qr_val;        // length p->total_qr
  T1 *lbx_val, *ubx_val;  // length p->total_x
  T1 *lbu_val, *ubu_val;  // length p->total_u
  T1 *lbg_val, *ubg_val;  // length p->total_g

  // Per-stage hat output buffers (written by eval)
  T1 *AB_hat_val, *CD_hat_val, *RSQ_hat_val;
  T1 *b_hat_val, *qr_hat_val;
  T1 *lbx_hat_val, *ubx_hat_val;
  T1 *lbu_hat_val, *ubu_hat_val;
  T1 *lbg_hat_val, *ubg_hat_val;

  // Condensed CSC + bounds + primal/dual scratch -- written by the
  // assemble/pack helpers and consumed by the inner solver.
  T1 *h_hat_csc;     // length p->nnz_RSQ_hat
  T1 *a_hat_csc;     // length p->nnz_a_hat_csc
  T1 *lbx, *ubx;  // length p->nx_total_hat
  T1 *lba, *uba;  // length p->na_total_hat
  T1 *x;
  T1 *lam_x;
  T1 *lam_a;

  // Lifted (full-horizon) primal/dual -- written by casadi_condensing_lift.
  // Caller copies these into res[CONIC_X/LAM_X/LAM_A] after lift returns.
  T1 *x_lifted;       // length p->total_qr
  T1 *lam_x_lifted;   // length p->total_qr
  T1 *lam_a_lifted;   // length p->total_b + p->total_g

  // User-side input snapshots -- caller assigns these before eval.
  // Eval reads from them and demultiplexes into the per-stage flat
  // input layout above.
  const T1 *g_orig;              // gradient (length sum_k nx[k]+nu[k])
  const T1 *lbx_orig, *ubx_orig; // bounds on z, interleaved [x_k; u_k]
  const T1 *lba_orig, *uba_orig; // bounds on Az, per stage [gap; path]

  // Internal eval scratch -- buffer pairs avoid aliasing in propagation.
  T1 *Phi;        // running transition matrix (nx_max x nx_max), col-major
  T1 *Phi_new;
  T1 *Gamma;      // running controllability (nx_max x nu_max_block), col-major
  T1 *Gamma_new;
  T1 *phi;        // running affine offset (nx_max)
  T1 *phi_new;
  T1 *Qxx_c;      // contiguous Qxx extract from RSQ[k] (nx_max x nx_max)
  T1 *Qxu_c;      // contiguous Qxu extract                (nx_max x nu_max_block)
  T1 *Quu_c;      // contiguous Quu extract                (nu_max_block x nu_max_block)
  T1 *gemm_xx;    // gemm scratch (nx_max x nx_max) for Q*Phi etc.
  T1 *gemm_xu;    // gemm scratch (nx_max x nu_max_block) for Q*Gamma, C*Phi, ...
  T1 *tmp_v;      // vector scratch (nx_max + nu_max_block)
};

// Inline offset helpers -- forward-declared because eval uses them and
// the C codegen output (after sanitize_source) requires definitions to
// appear before use.  Bodies are at the bottom of this file.  These read
// from data (not prob), since the per-stage *_hat arrays live on data.
template<typename T1>
static casadi_int casadi_condensing_off_lbx(const casadi_condensing_data<T1>* d,
                                            casadi_int K);
template<typename T1>
static casadi_int casadi_condensing_off_lbu(const casadi_condensing_data<T1>* d,
                                            casadi_int K);
template<typename T1>
static casadi_int casadi_condensing_off_lbg(const casadi_condensing_data<T1>* d,
                                            casadi_int K);

// SYMBOL "condensing_work"
//   Bumps sz_w (and sz_iw) by ALL the workspace this transform needs --
//   per-stage flat input/output buffers, per-stage hat outputs, internal
//   eval scratch, condensed-CSC scratch, and condensed bounds + primal/
//   dual scratch.  The caller does not need to know any of these sizes.
template<typename T1>
void casadi_condensing_work(const casadi_condensing_prob<T1>* p,
                            casadi_int* sz_iw, casadi_int* sz_w) {
  casadi_int nx2 = p->nx_max * p->nx_max;
  casadi_int nxnu = p->nx_max * p->nu_max_block;
  casadi_int nu2 = p->nu_max_block * p->nu_max_block;
  // Per-condensed-stage derived dim/offset/block-descriptor arrays
  // (populated by casadi_condensing_set_work).
  *sz_iw += 3 * (p->N_hat + 1);          // nx_hat, nu_hat, ng_hat
  *sz_iw += 3 * (p->N_hat + 1);          // *_hat_offsets
  *sz_iw += 4 * p->N_hat;                // AB_hat: 4 ints per ocp_block
  *sz_iw += 4 * (p->N_hat + 1);          // CD_hat
  *sz_iw += 4 * (p->N_hat + 1);          // RSQ_hat
  // Per-stage flat inputs
  *sz_w += p->nnz_RSQ;
  *sz_w += p->nnz_AB;
  *sz_w += p->nnz_CD;
  *sz_w += p->total_qr;
  *sz_w += p->total_b;
  *sz_w += 2 * p->total_x;
  *sz_w += 2 * p->total_u;
  *sz_w += 2 * p->total_g;
  // Per-stage hat outputs
  *sz_w += p->nnz_RSQ_hat;
  *sz_w += p->nnz_AB_hat;
  *sz_w += p->nnz_CD_hat;
  *sz_w += p->nx_total_hat;            // qr_hat
  *sz_w += p->total_b_hat;
  *sz_w += 2 * p->total_x_hat;
  *sz_w += 2 * p->total_u_hat;
  *sz_w += 2 * p->total_g_hat;
  // Internal eval scratch
  *sz_w += 2 * nx2;                    // Phi, Phi_new
  *sz_w += 2 * nxnu;                   // Gamma, Gamma_new
  *sz_w += 2 * p->nx_max;              // phi, phi_new
  *sz_w += nx2;                        // Qxx_c
  *sz_w += nxnu;                       // Qxu_c
  *sz_w += nu2;                        // Quu_c
  *sz_w += nx2;                        // gemm_xx
  *sz_w += nxnu;                       // gemm_xu
  *sz_w += p->nx_max + p->nu_max_block;
  // Condensed CSC + bounds + primal/dual scratch
  *sz_w += p->nnz_RSQ_hat;             // h_hat_csc
  *sz_w += p->nnz_a_hat_csc;
  *sz_w += 2 * p->nx_total_hat;        // lbx, ubx
  *sz_w += 2 * p->na_total_hat;        // lba, uba
  *sz_w += p->nx_total_hat;            // x
  *sz_w += p->nx_total_hat;            // lam_x
  *sz_w += p->na_total_hat;            // lam_a
  // x_lifted / lam_x_lifted / lam_a_lifted alias onto the dead per-stage
  // bounds region (lbx_val..ubg_val, contiguous, ends right before the
  // *_hat_val region in w).  No extra claim needed -- see set_work.
}

// SYMBOL "condensing_set_work"
//   Claims ALL workspace from iw and w (per-condensed-stage derived
//   dim/offset/block-descriptor arrays, per-stage flat in/out, internal
//   scratch, condensed CSC, condensed bounds, primal/dual) AND populates
//   the derived per-stage arrays from the prob recurrence.  After this
//   call d->RSQ_val, d->AB_val, ..., d->h_hat_csc, d->x, d->nx_hat, ...
//   are all ready to use.
template<typename T1>
void casadi_condensing_set_work(casadi_condensing_data<T1>* d,
    const T1*** arg, T1*** res, casadi_int** iw, T1** w) {
  (void)arg; (void)res;
  const casadi_condensing_prob<T1>* p = d->prob;
  casadi_int nx2 = p->nx_max * p->nx_max;
  casadi_int nxnu = p->nx_max * p->nu_max_block;
  casadi_int nu2 = p->nu_max_block * p->nu_max_block;

  // Per-condensed-stage derived arrays.  Claim from iw, then fill via
  // the same per-K recurrence that casadi_condensing_setup uses for the
  // scalars.  Block-descriptor structs (4 casadi_int fields each) ride
  // in iw via reinterpret_cast.
  d->nx_hat = *iw; *iw += p->N_hat + 1;
  d->nu_hat = *iw; *iw += p->N_hat + 1;
  d->ng_hat = *iw; *iw += p->N_hat + 1;
  d->AB_hat_offsets  = *iw; *iw += p->N_hat + 1;
  d->CD_hat_offsets  = *iw; *iw += p->N_hat + 1;
  d->RSQ_hat_offsets = *iw; *iw += p->N_hat + 1;
  d->AB_hat  = reinterpret_cast<casadi_ocp_block*>(*iw); *iw += 4 * p->N_hat;
  d->CD_hat  = reinterpret_cast<casadi_ocp_block*>(*iw); *iw += 4 * (p->N_hat + 1);
  d->RSQ_hat = reinterpret_cast<casadi_ocp_block*>(*iw); *iw += 4 * (p->N_hat + 1);
  {
    casadi_int K, j, k_a, k_b, M_, sum_nu, sum_ng_path, sum_ng_lift;
    casadi_int nx_K, nxu_K, off_AB = 0, off_CD = 0, off_RSQ = 0;
    for (K = 0; K <= p->N_hat; ++K) d->nx_hat[K] = p->nx[p->M[K]];
    for (K = 0; K < p->N_hat; ++K) {
      k_a = p->M[K];
      k_b = p->M[K + 1];
      M_ = k_b - k_a;
      sum_nu = 0;
      sum_ng_path = 0;
      sum_ng_lift = 0;
      for (j = 0; j < M_; ++j) {
        sum_nu += p->nu[k_a + j];
        sum_ng_path += p->ng[k_a + j];
        if (j > 0) sum_ng_lift += p->nx[k_a + j];
      }
      d->nu_hat[K] = sum_nu;
      d->ng_hat[K] = sum_ng_path + sum_ng_lift;
      nx_K = d->nx_hat[K];
      nxu_K = nx_K + sum_nu;
      d->AB_hat[K].offset_r = 0;
      d->AB_hat[K].offset_c = 0;
      d->AB_hat[K].rows = d->nx_hat[K + 1];
      d->AB_hat[K].cols = nxu_K;
      d->AB_hat_offsets[K] = off_AB;
      off_AB += d->AB_hat[K].rows * d->AB_hat[K].cols;
      d->CD_hat[K].offset_r = 0;
      d->CD_hat[K].offset_c = 0;
      d->CD_hat[K].rows = d->ng_hat[K];
      d->CD_hat[K].cols = nxu_K;
      d->CD_hat_offsets[K] = off_CD;
      off_CD += d->CD_hat[K].rows * d->CD_hat[K].cols;
      d->RSQ_hat[K].offset_r = 0;
      d->RSQ_hat[K].offset_c = 0;
      d->RSQ_hat[K].rows = nxu_K;
      d->RSQ_hat[K].cols = nxu_K;
      d->RSQ_hat_offsets[K] = off_RSQ;
      off_RSQ += nxu_K * nxu_K;
    }
    // Terminal stage K = N_hat: only path constraints.
    d->nu_hat[p->N_hat] = 0;
    d->ng_hat[p->N_hat] = p->ng[p->N];
    nx_K = d->nx_hat[p->N_hat];
    d->AB_hat_offsets[p->N_hat] = off_AB;
    d->CD_hat[p->N_hat].offset_r = 0;
    d->CD_hat[p->N_hat].offset_c = 0;
    d->CD_hat[p->N_hat].rows = d->ng_hat[p->N_hat];
    d->CD_hat[p->N_hat].cols = nx_K;
    d->CD_hat_offsets[p->N_hat] = off_CD;
    d->RSQ_hat[p->N_hat].offset_r = 0;
    d->RSQ_hat[p->N_hat].offset_c = 0;
    d->RSQ_hat[p->N_hat].rows = nx_K;
    d->RSQ_hat[p->N_hat].cols = nx_K;
    d->RSQ_hat_offsets[p->N_hat] = off_RSQ;
  }

  // Per-stage flat inputs
  d->RSQ_val = *w; *w += p->nnz_RSQ;
  d->AB_val  = *w; *w += p->nnz_AB;
  d->CD_val  = *w; *w += p->nnz_CD;
  d->qr_val  = *w; *w += p->total_qr;
  d->b_val   = *w; *w += p->total_b;
  d->lbx_val = *w; *w += p->total_x;
  d->ubx_val = *w; *w += p->total_x;
  d->lbu_val = *w; *w += p->total_u;
  d->ubu_val = *w; *w += p->total_u;
  d->lbg_val = *w; *w += p->total_g;
  d->ubg_val = *w; *w += p->total_g;

  // Per-stage hat outputs
  d->RSQ_hat_val = *w; *w += p->nnz_RSQ_hat;
  d->AB_hat_val  = *w; *w += p->nnz_AB_hat;
  d->CD_hat_val  = *w; *w += p->nnz_CD_hat;
  d->qr_hat_val  = *w; *w += p->nx_total_hat;
  d->b_hat_val   = *w; *w += p->total_b_hat;
  d->lbx_hat_val = *w; *w += p->total_x_hat;
  d->ubx_hat_val = *w; *w += p->total_x_hat;
  d->lbu_hat_val = *w; *w += p->total_u_hat;
  d->ubu_hat_val = *w; *w += p->total_u_hat;
  d->lbg_hat_val = *w; *w += p->total_g_hat;
  d->ubg_hat_val = *w; *w += p->total_g_hat;

  // Internal eval scratch
  d->Phi       = *w; *w += nx2;
  d->Phi_new   = *w; *w += nx2;
  d->Gamma     = *w; *w += nxnu;
  d->Gamma_new = *w; *w += nxnu;
  d->phi       = *w; *w += p->nx_max;
  d->phi_new   = *w; *w += p->nx_max;
  d->Qxx_c     = *w; *w += nx2;
  d->Qxu_c     = *w; *w += nxnu;
  d->Quu_c     = *w; *w += nu2;
  d->gemm_xx   = *w; *w += nx2;
  d->gemm_xu   = *w; *w += nxnu;
  d->tmp_v     = *w; *w += p->nx_max + p->nu_max_block;

  // Condensed CSC + bounds + scratch primal/dual
  d->h_hat_csc  = *w; *w += p->nnz_RSQ_hat;
  d->a_hat_csc  = *w; *w += p->nnz_a_hat_csc;
  d->lbx   = *w; *w += p->nx_total_hat;
  d->ubx   = *w; *w += p->nx_total_hat;
  d->lba   = *w; *w += p->na_total_hat;
  d->uba   = *w; *w += p->na_total_hat;
  d->x     = *w; *w += p->nx_total_hat;
  d->lam_x = *w; *w += p->nx_total_hat;
  d->lam_a = *w; *w += p->na_total_hat;

  // Lift output ptrs alias onto the per-stage bounds region (dead at
  // lift time; lift reads only x/lam_x/lam_a + {AB,b,RSQ,qr}_val +
  // Phi_new/Gamma_new scratch).  The contiguous w-region from lbx_val
  // onward spans 2*(total_x+total_u+total_g) bytes for the bounds plus
  // the entire *_hat_val region after that, so 2*total_qr + total_b +
  // total_g fits unconditionally.
  d->x_lifted     = d->lbx_val;
  d->lam_x_lifted = d->lbx_val + p->total_qr;
  d->lam_a_lifted = d->lbx_val + 2 * p->total_qr;
}

// Helper: zero a buffer of length n.
template<typename T1>
static void casadi_condensing_zero(T1* x, casadi_int n) {
  casadi_int i;
  for (i = 0; i < n; ++i) x[i] = 0;
}

// Helper: copy n elements from src to dst.
template<typename T1>
static void casadi_condensing_copy(const T1* src, casadi_int n, T1* dst) {
  casadi_int i;
  for (i = 0; i < n; ++i) dst[i] = src[i];
}

// Helper: y(nrow, ncol) := alpha * y; column-major, full dense.
template<typename T1>
static void casadi_condensing_scale(T1 alpha, T1* y, casadi_int n) {
  casadi_int i;
  for (i = 0; i < n; ++i) y[i] *= alpha;
}

// SYMBOL "condensing_eval"
//   Run one condensing pass: read user inputs from d->g_orig,
//   d->lbx_orig, d->ubx_orig, d->lba_orig, d->uba_orig (caller assigns before
//   calling), copy/demultiplex into d's per-stage flat layout,
//   condense, and pack into the flat-dense Conic-input form
//   (h_hat_csc, a_hat_csc, lbx/ubx, lba/uba) ready for the
//   inner solver.  RSQ_val / AB_val / CD_val must also be pre-filled
//   by the caller (typically via casadi_project from the user's H, A
//   inputs).  Returns 0 on success.
template<typename T1>
int casadi_condensing_eval(casadi_condensing_data<T1>* d) {
  const casadi_condensing_prob<T1>* p = d->prob;
  casadi_int K, j, k, k_a, k_b, M_;
  casadi_int nx_K, nu_K, nxu_K, nx_k, nu_k, ng_k, nx_kp1;
  casadi_int cum_nu, off_lbu_block, off_lbg_block, off_lbx, off_lbu, off_lbg;
  casadi_int off_b = 0, off_qr = 0;            // running input vector offsets

  // Copy/demultiplex Conic-format user inputs -> d's per-stage flat layout.
  {
    casadi_int kk, ii, src = 0, dx = 0, du = 0, ob = 0, og = 0, sa = 0;
    casadi_int total_qr = p->total_qr;
    for (ii = 0; ii < total_qr; ++ii) d->qr_val[ii] = d->g_orig[ii];
    // lbx/ubx_orig are interleaved [x_k; u_k]; demultiplex.
    for (kk = 0; kk <= p->N; ++kk) {
      for (ii = 0; ii < p->nx[kk]; ++ii) {
        d->lbx_val[dx + ii] = d->lbx_orig[src + ii];
        d->ubx_val[dx + ii] = d->ubx_orig[src + ii];
      }
      src += p->nx[kk];  dx  += p->nx[kk];
      for (ii = 0; ii < p->nu[kk]; ++ii) {
        d->lbu_val[du + ii] = d->lbx_orig[src + ii];
        d->ubu_val[du + ii] = d->ubx_orig[src + ii];
      }
      src += p->nu[kk];  du  += p->nu[kk];
    }
    // lba/uba_orig: per-stage (nx[k+1] gap + ng[k] path), terminal ng[N].
    for (kk = 0; kk < p->N; ++kk) {
      for (ii = 0; ii < p->nx[kk + 1]; ++ii) {
        d->b_val[ob + ii] = -d->lba_orig[sa + ii];
      }
      sa += p->nx[kk + 1];  ob += p->nx[kk + 1];
      for (ii = 0; ii < p->ng[kk]; ++ii) {
        d->lbg_val[og + ii] = d->lba_orig[sa + ii];
        d->ubg_val[og + ii] = d->uba_orig[sa + ii];
      }
      sa += p->ng[kk];  og += p->ng[kk];
    }
    for (ii = 0; ii < p->ng[p->N]; ++ii) {
      d->lbg_val[og + ii] = d->lba_orig[sa + ii];
      d->ubg_val[og + ii] = d->uba_orig[sa + ii];
    }
  }
  casadi_int row_off;                          // running row offset within CD_hat
  T1 *Phi, *Gamma, *phi, *phi_new, *tmp_swap;
  const T1 *A, *B, *Qxx, *Qxu, *Quu, *q_k, *r_k, *C, *D, *b_k;
  T1 *RSQ_hat, *Hxx, *Hxu, *Huu, *hx, *hu, *AB_hat, *CD_hat, *qr_hat, *b_hat;
  T1 *Phi_out, *Gamma_out;
  casadi_int i, jj, ii, kk;

  // Pre-compute cumulative offsets for input vectors keyed by stage.
  // We advance them in-place as we walk stages.
  off_lbx = 0;  // accumulates sum nx[k] up to current stage
  off_lbu = 0;  // accumulates sum nu[k] up to current stage
  off_lbg = 0;  // accumulates sum ng[k] up to current stage

  for (K = 0; K < p->N_hat; ++K) {
    k_a = p->M[K];
    k_b = p->M[K + 1];
    M_ = k_b - k_a;
    nx_K = d->nx_hat[K];
    nu_K = d->nu_hat[K];
    nxu_K = nx_K + nu_K;

    // Initialise Phi = I_{nx_K}, Gamma = 0 (no columns yet), phi = 0
    Phi   = d->Phi;
    Gamma = d->Gamma;
    phi   = d->phi;
    casadi_condensing_zero(Phi, nx_K * nx_K);
    for (i = 0; i < nx_K; ++i) Phi[i + i * nx_K] = 1;
    casadi_condensing_zero(phi, nx_K);
    // Gamma: no columns at j=0, nothing to initialise.

    // Output buffers for this condensed block
    AB_hat  = d->AB_hat_val  + d->AB_hat_offsets[K];
    CD_hat  = d->CD_hat_val  + d->CD_hat_offsets[K];
    RSQ_hat = d->RSQ_hat_val + d->RSQ_hat_offsets[K];
    qr_hat  = d->qr_hat_val;        // walked stage-by-stage; offset below
    b_hat   = d->b_hat_val;         // walked likewise

    // Find offsets for output qr_hat[K] and b_hat[K].
    // These mirror the per-stage concatenation convention.
    {
      casadi_int t_off_qr = 0, t_off_b = 0;
      for (i = 0; i < K; ++i) {
        t_off_qr += d->nx_hat[i] + d->nu_hat[i];
        t_off_b  += d->nx_hat[i + 1];
      }
      qr_hat = d->qr_hat_val + t_off_qr;
      b_hat  = d->b_hat_val  + t_off_b;
    }

    // RSQ_hat block: row-major view as Hxx (nx_K x nx_K), Hxu (nx_K x nu_K),
    // Huu (nu_K x nu_K) within a single (nxu_K x nxu_K) col-major matrix.
    casadi_condensing_zero(RSQ_hat, nxu_K * nxu_K);
    casadi_condensing_zero(qr_hat, nxu_K);
    Hxx = RSQ_hat;                                  // top-left  nx_K x nx_K
    Hxu = RSQ_hat + nx_K * nxu_K;                   // top-right (in xu cols)
    Huu = RSQ_hat + nx_K * nxu_K + nx_K;            // bottom-right
    hx  = qr_hat;
    hu  = qr_hat + nx_K;

    // CD_hat: zero whole block; we write contiguous row groups
    casadi_condensing_zero(CD_hat, d->ng_hat[K] * nxu_K);

    // Bounds: copy xi-bounds from stage k_a; collect lbu/ubu/lbg (and lifted
    // state bounds) as we walk j.
    casadi_condensing_copy(d->lbx_val + off_lbx, nx_K, d->lbx_hat_val + casadi_condensing_off_lbx(d, K));
    casadi_condensing_copy(d->ubx_val + off_lbx, nx_K, d->ubx_hat_val + casadi_condensing_off_lbx(d, K));
    // (We use a small inline helper below to compute lbx_hat offsets.)

    cum_nu = 0;          // running width of Gamma at start of step j
    row_off = 0;         // running row offset into CD_hat for path-ineq + lifts
    off_lbu_block = 0;   // running offset within lbu_hat[K] / ubu_hat[K]
    off_lbg_block = 0;   // running offset within lbg_hat[K] / ubg_hat[K]

    // Walk inner stages j = 0, 1, ..., M-1
    for (j = 0; j < M_; ++j) {
      k = k_a + j;
      nx_k   = p->nx[k];
      nu_k   = p->nu[k];
      ng_k   = p->ng[k];
      nx_kp1 = p->nx[k + 1];

      A   = d->AB_val  + p->AB_offsets[k];                     // nx_kp1 x nx_k
      B   = A + nx_kp1 * nx_k;                                 // nx_kp1 x nu_k
      Qxx = d->RSQ_val + p->RSQ_offsets[k];                    // nx_k x nx_k
      Qxu = Qxx + nx_k * (nx_k + nu_k);                        // nx_k x nu_k (in (nx_k+nu_k) col-stride)
      Quu = Qxx + nx_k * (nx_k + nu_k) + nx_k;                 // nu_k x nu_k
      // Wait -- correction: the RSQ block is col-major (nx_k+nu_k)x(nx_k+nu_k).
      // Top-left Qxx occupies rows [0,nx_k), cols [0,nx_k); each col has stride
      // (nx_k+nu_k).  Using a flat pointer with stride = nx_k+nu_k requires
      // explicit indexing below; we cannot use casadi_mtimes_dense on Qxx as
      // a contiguous (nx_k x nx_k).  See "rebuild Q-blocks" below.
      (void)Qxx; (void)Qxu; (void)Quu;  // we re-extract these below
      C   = d->CD_val  + p->CD_offsets[k];                     // ng_k x nx_k
      D   = C + ng_k * nx_k;                                   // ng_k x nu_k
      q_k = d->qr_val + off_qr;                                // length nx_k
      r_k = q_k + nx_k;                                        // length nu_k
      b_k = d->b_val + off_b;                                  // length nx_kp1

      // Extract Qxx, Qxu, Quu into contiguous buffers because the RSQ
      // block has column-stride (nx_k+nu_k) and our gemm primitive
      // expects contiguous column-major operands.
      {
        const casadi_int sQ = nx_k + nu_k;
        const T1* RSQ_full = d->RSQ_val + p->RSQ_offsets[k];
        T1 *Qxx_c = d->Qxx_c;
        T1 *Qxu_c = d->Qxu_c;
        T1 *Quu_c = d->Quu_c;
        for (jj = 0; jj < nx_k; ++jj) {
          for (ii = 0; ii < nx_k; ++ii)
            Qxx_c[ii + jj * nx_k] = RSQ_full[ii + jj * sQ];
        }
        for (jj = 0; jj < nu_k; ++jj) {
          for (ii = 0; ii < nx_k; ++ii)
            Qxu_c[ii + jj * nx_k] = RSQ_full[ii + (nx_k + jj) * sQ];
        }
        for (jj = 0; jj < nu_k; ++jj) {
          for (ii = 0; ii < nu_k; ++ii)
            Quu_c[ii + jj * nu_k] = RSQ_full[(nx_k + ii) + (nx_k + jj) * sQ];
        }

        // ---- Cost contribution from stage k = k_a + j ----
        //
        // Substitute x_k = Phi*xi + Gamma*nu_<j + phi.
        // Let Z_x := Phi (nx_k x nx_K), Z_v := Gamma (nx_k x cum_nu),
        // Z_p := phi (nx_k).  Then in (xi, nu) coords the stage cost adds:
        //   Hxx += Phi'  Qxx Phi
        //   Hxu_left += Phi' Qxx Gamma          -> goes to Hxu cols [0, cum_nu)
        //   Hxu_right += Phi' Qxu               -> goes to Hxu cols [cum_nu, cum_nu+nu_k)
        //   Huu_TL += Gamma' Qxx Gamma          -> Huu rows/cols [0, cum_nu)
        //   Huu_TR += Gamma' Qxu                -> Huu rows [0,cum_nu), cols [cum_nu, +nu_k)
        //   Huu_BR += Quu                       -> Huu rows/cols [cum_nu, +nu_k)
        //   hx  += Phi'   * (Qxx*phi + q_k)
        //   hu_left += Gamma' * (Qxx*phi + q_k)
        //   hu_right += Qxu' * phi + r_k
        //
        // Compute s := Qxx*phi + q_k     (length nx_k)
        // and    t := Qxu' * phi + r_k    (length nu_k)
        {
          T1 *s = d->tmp_v;                    // length nx_k
          T1 *t = d->tmp_v + nx_k;             // length nu_k
          for (ii = 0; ii < nx_k; ++ii) s[ii] = q_k[ii];
          for (ii = 0; ii < nu_k; ++ii) t[ii] = r_k[ii];
          // s += Qxx * phi
          casadi_blas_mtimes(Qxx_c, nx_k, nx_k, phi, 1, s, 0);
          // t += Qxu' * phi
          casadi_blas_mtimes(Qxu_c, nx_k, nu_k, phi, 1, t, 1);

          // hx += Phi' * s   (Phi is nx_k x nx_K, s is nx_k, hx is nx_K)
          casadi_blas_mtimes(Phi, nx_k, nx_K, s, 1, hx, 1);
          // hu_left += Gamma' * s   (Gamma is nx_k x cum_nu)
          if (cum_nu > 0) {
            casadi_blas_mtimes(Gamma, nx_k, cum_nu, s, 1, hu, 1);
          }
          // hu_right += t   (write into hu[cum_nu .. cum_nu+nu_k-1])
          for (ii = 0; ii < nu_k; ++ii) hu[cum_nu + ii] += t[ii];
        }

        // Hessian contributions.  Use d->gemm_xx and d->gemm_xu as scratch.
        // (1) Hxx += Phi' Qxx Phi.  Hxx is the top-left sub-block of
        // RSQ_hat with column stride nxu_K, so we write one column at a
        // time (column kk of Hxx is Phi' * (Q*Phi)[:,kk]).
        {
          T1 *tmp = d->gemm_xx;       // (nx_k x nx_K)
          casadi_condensing_zero(tmp, nx_k * nx_K);
          casadi_blas_mtimes(Qxx_c, nx_k, nx_k, Phi, nx_K, tmp, 0);
          for (kk = 0; kk < nx_K; ++kk) {
            T1 *hxx_col = Hxx + kk * nxu_K;
            casadi_blas_mtimes(Phi, nx_k, nx_K, tmp + kk * nx_k, 1, hxx_col, 1);
          }
        }

        // (2) Hxu_left  (cols [0, cum_nu))      += Phi' Qxx Gamma
        // (3) Hxu_right (cols [cum_nu, +nu_k))  += Phi' Qxu
        {
          T1 *tmp = d->gemm_xu;       // (nx_k x cum_nu) when used
          if (cum_nu > 0) {
            casadi_condensing_zero(tmp, nx_k * cum_nu);
            casadi_blas_mtimes(Qxx_c, nx_k, nx_k, Gamma, cum_nu, tmp, 0);
            // Hxu_left columns are not contiguous in RSQ_hat (col-stride nxu_K,
            // length nx_K).  Write one column at a time.
            for (kk = 0; kk < cum_nu; ++kk) {
              T1 *hxu_col = Hxu + kk * nxu_K;
              casadi_blas_mtimes(Phi, nx_k, nx_K, tmp + kk * nx_k, 1, hxu_col, 1);
            }
          }
          // Hxu_right += Phi' Qxu  -- write one column at a time
          for (kk = 0; kk < nu_k; ++kk) {
            T1 *hxu_col = Hxu + (cum_nu + kk) * nxu_K;
            casadi_blas_mtimes(Phi, nx_k, nx_K, Qxu_c + kk * nx_k, 1, hxu_col, 1);
          }
        }

        // (4) Huu_TL += Gamma' Qxx Gamma
        // (5) Huu_TR += Gamma' Qxu
        // (6) Huu_BR += Quu
        // Huu is the bottom-right (nu_K x nu_K) sub-block of RSQ_hat;
        // Huu[i,j] = RSQ_hat[(nx_K + i) + (nx_K + j) * nxu_K].
        if (cum_nu > 0) {
          T1 *tmp = d->gemm_xu;
          casadi_condensing_zero(tmp, nx_k * cum_nu);
          casadi_blas_mtimes(Qxx_c, nx_k, nx_k, Gamma, cum_nu, tmp, 0);
          for (kk = 0; kk < cum_nu; ++kk) {
            T1 *huu_col = Huu + kk * nxu_K;
            casadi_blas_mtimes(Gamma, nx_k, cum_nu, tmp + kk * nx_k, 1, huu_col, 1);
          }
          for (kk = 0; kk < nu_k; ++kk) {
            T1 *huu_col = Huu + (cum_nu + kk) * nxu_K;
            casadi_blas_mtimes(Gamma, nx_k, cum_nu, Qxu_c + kk * nx_k, 1, huu_col, 1);
          }
        }
        // (6): Huu_BR += Quu   (nu_k x nu_k) into rows/cols [cum_nu, +nu_k)
        for (kk = 0; kk < nu_k; ++kk) {
          T1 *huu_col = Huu + (cum_nu + kk) * nxu_K;  // length nxu_K
          for (ii = 0; ii < nu_k; ++ii) {
            huu_col[cum_nu + ii] += Quu_c[ii + kk * nu_k];
          }
        }

        // Symmetrize Huu_TR into Huu_BL (mirror to ensure full symmetry of
        // RSQ_hat block; some downstream solvers assume symmetric storage).
        // Skip if cum_nu == 0.
        if (cum_nu > 0) {
          for (kk = 0; kk < nu_k; ++kk) {        // dst col (cum_nu+kk)
            for (ii = 0; ii < cum_nu; ++ii) {
              T1 v = Huu[ii + (cum_nu + kk) * nxu_K];
              Huu[(cum_nu + kk) + ii * nxu_K] = v;
            }
          }
        }
      }

      // ---- Path inequality contribution at stage k ----
      // Original constraint: lbg <= C_k x_k + D_k u_k <= ubg
      // After substitution:  C_k Phi xi + C_k Gamma nu_<j + D_k u_k
      //                       in [lbg - C_k phi, ubg - C_k phi].
      // Append ng_k rows starting at row_off in CD_hat[K].
      if (ng_k > 0) {
        // Compute C_k * phi (length ng_k)
        T1 *Cphi = d->tmp_v;
        casadi_condensing_zero(Cphi, ng_k);
        casadi_blas_mtimes(C, ng_k, nx_k, phi, 1, Cphi, 0);

        // (a) Columns [0, nx_K) of CD_hat block:  C_k * Phi
        // tmp (ng_k x nx_K) = C_k * Phi
        {
          T1 *tmp = d->gemm_xu;            // borrow (nx_max x nu_max_block); plenty
          casadi_condensing_zero(tmp, ng_k * nx_K);
          casadi_blas_mtimes(C, ng_k, nx_k, Phi, nx_K, tmp, 0);
          // Write into CD_hat at rows [row_off, row_off+ng_k), cols [0, nx_K)
          for (kk = 0; kk < nx_K; ++kk) {
            T1 *col = CD_hat + kk * d->ng_hat[K];   // ng_hat-strided col-major
            for (ii = 0; ii < ng_k; ++ii)
              col[row_off + ii] = tmp[ii + kk * ng_k];
          }
        }
        // (b) Columns [nx_K, nx_K+cum_nu) of CD_hat:  C_k * Gamma
        if (cum_nu > 0) {
          T1 *tmp = d->gemm_xu;
          casadi_condensing_zero(tmp, ng_k * cum_nu);
          casadi_blas_mtimes(C, ng_k, nx_k, Gamma, cum_nu, tmp, 0);
          for (kk = 0; kk < cum_nu; ++kk) {
            T1 *col = CD_hat + (nx_K + kk) * d->ng_hat[K];
            for (ii = 0; ii < ng_k; ++ii)
              col[row_off + ii] = tmp[ii + kk * ng_k];
          }
        }
        // (c) Columns [nx_K+cum_nu, nx_K+cum_nu+nu_k) of CD_hat: D_k
        for (kk = 0; kk < nu_k; ++kk) {
          T1 *col = CD_hat + (nx_K + cum_nu + kk) * d->ng_hat[K];
          for (ii = 0; ii < ng_k; ++ii)
            col[row_off + ii] = D[ii + kk * ng_k];
        }
        // RHS adjustments: lbg_hat[row_off+i] = lbg[i] - Cphi[i]; ubg likewise.
        for (ii = 0; ii < ng_k; ++ii) {
          d->lbg_hat_val[casadi_condensing_off_lbg(d, K) + off_lbg_block + ii]
              = d->lbg_val[off_lbg + ii] - Cphi[ii];
          d->ubg_hat_val[casadi_condensing_off_lbg(d, K) + off_lbg_block + ii]
              = d->ubg_val[off_lbg + ii] - Cphi[ii];
        }
        off_lbg_block += ng_k;
        row_off += ng_k;
      }

      // ---- Lifted state bound contribution for j > 0 ----
      // Original: lbx[k_a+j] <= x_{k_a+j} <= ubx[k_a+j]
      // After substitution: lbx - phi <= Phi xi + Gamma nu_<j <= ubx - phi
      // Add nx_k rows.  Rows write [Phi | Gamma | 0 ...] columns and shifted
      // bounds.
      if (j > 0) {
        // Columns [0, nx_K): copy Phi rows
        for (kk = 0; kk < nx_K; ++kk) {
          T1 *col = CD_hat + kk * d->ng_hat[K];
          for (ii = 0; ii < nx_k; ++ii)
            col[row_off + ii] = Phi[ii + kk * nx_k];
        }
        // Columns [nx_K, nx_K + cum_nu): copy Gamma rows
        for (kk = 0; kk < cum_nu; ++kk) {
          T1 *col = CD_hat + (nx_K + kk) * d->ng_hat[K];
          for (ii = 0; ii < nx_k; ++ii)
            col[row_off + ii] = Gamma[ii + kk * nx_k];
        }
        // Remaining u-columns are zero (already from initial zeroing of CD_hat).
        // RHS shifts.
        for (ii = 0; ii < nx_k; ++ii) {
          d->lbg_hat_val[casadi_condensing_off_lbg(d, K) + off_lbg_block + ii]
              = d->lbx_val[off_lbx + ii] - phi[ii];
          d->ubg_hat_val[casadi_condensing_off_lbg(d, K) + off_lbg_block + ii]
              = d->ubx_val[off_lbx + ii] - phi[ii];
        }
        off_lbg_block += nx_k;
        row_off += nx_k;
      }

      // ---- Control bounds: append lbu_k, ubu_k ----
      for (ii = 0; ii < nu_k; ++ii) {
        d->lbu_hat_val[casadi_condensing_off_lbu(d, K) + off_lbu_block + ii]
            = d->lbu_val[off_lbu + ii];
        d->ubu_hat_val[casadi_condensing_off_lbu(d, K) + off_lbu_block + ii]
            = d->ubu_val[off_lbu + ii];
      }
      off_lbu_block += nu_k;

      // ---- Advance Phi, Gamma, phi to step j+1 ----
      // Pick the *other* buffer than the one currently pointed to, so the
      // matmul never aliases input and output.  After writing, the local
      // pointer follows.
      Phi_out   = (Phi   == d->Phi)   ? d->Phi_new   : d->Phi;
      Gamma_out = (Gamma == d->Gamma) ? d->Gamma_new : d->Gamma;
      phi_new   = (phi   == d->phi)   ? d->phi_new   : d->phi;

      // phi_new = A * phi + b_k
      casadi_condensing_copy(b_k, nx_kp1, phi_new);
      casadi_blas_mtimes(A, nx_kp1, nx_k, phi, 1, phi_new, 0);

      // Phi_new = A * Phi  (nx_kp1 x nx_K)
      casadi_condensing_zero(Phi_out, nx_kp1 * nx_K);
      casadi_blas_mtimes(A, nx_kp1, nx_k, Phi, nx_K, Phi_out, 0);

      // Gamma_new = [A * Gamma   B_k]   (nx_kp1 x (cum_nu + nu_k))
      casadi_condensing_zero(Gamma_out, nx_kp1 * (cum_nu + nu_k));
      if (cum_nu > 0) {
        casadi_blas_mtimes(A, nx_kp1, nx_k, Gamma, cum_nu, Gamma_out, 0);
      }
      // Append B_k as the last nu_k columns
      for (kk = 0; kk < nu_k; ++kk) {
        const T1 *src = B + kk * nx_kp1;
        T1 *dst = Gamma_out + (cum_nu + kk) * nx_kp1;
        for (ii = 0; ii < nx_kp1; ++ii) dst[ii] = src[ii];
      }

      // Local pointers now follow the freshly-written buffers.
      Phi   = Phi_out;
      Gamma = Gamma_out;
      phi   = phi_new;
      (void)tmp_swap;

      // Advance input vector offsets (b: nx_kp1, qr: nx_k+nu_k, lbx: nx_k,
      // lbu: nu_k, lbg: ng_k)
      off_b  += nx_kp1;
      off_qr += nx_k + nu_k;
      off_lbx += nx_k;
      off_lbu += nu_k;
      off_lbg += ng_k;
      cum_nu += nu_k;
    }   // end inner loop over j

    // After the M-step recursion, Phi (nx_{k_b} x nx_K), Gamma (nx_{k_b} x nu_K),
    // phi (nx_{k_b}) hold the M-step transition.  Write them out.
    // AB_hat[K] columns [0, nx_K): Phi
    // AB_hat[K] columns [nx_K, nx_K+nu_K): Gamma
    {
      casadi_int rows_out = d->nx_hat[K + 1];
      for (kk = 0; kk < nx_K; ++kk) {
        for (ii = 0; ii < rows_out; ++ii)
          AB_hat[ii + kk * rows_out] = Phi[ii + kk * rows_out];
      }
      for (kk = 0; kk < nu_K; ++kk) {
        for (ii = 0; ii < rows_out; ++ii)
          AB_hat[ii + (nx_K + kk) * rows_out] = Gamma[ii + kk * rows_out];
      }
      // b_hat[K] = phi
      for (ii = 0; ii < rows_out; ++ii) b_hat[ii] = phi[ii];
    }

    // Symmetrize RSQ_hat[K]: mirror the upper triangle into the lower
    // triangle so the output Hessian is fully symmetric.  We have
    // populated:
    //   Hxx (top-left, nx_K x nx_K, all entries) -- written via Phi'QPhi
    //   Hxu (top-right, nx_K x nu_K, all entries) -- written via Phi'(QGamma|Qxu)
    //   Huu (bottom-right, nu_K x nu_K), incl. internal symmetrization.
    // Mirror Hxu into Hux (bottom-left).  Hxx and Huu are already symmetric
    // because the underlying R_tilde block is symmetric and we accumulate
    // R_tilde^{1/2}-type forms.
    for (jj = 0; jj < nu_K; ++jj) {
      for (ii = 0; ii < nx_K; ++ii) {
        T1 v = RSQ_hat[ii + (nx_K + jj) * nxu_K];   // Hxu[ii, jj]
        RSQ_hat[(nx_K + jj) + ii * nxu_K] = v;       // Hux[jj, ii]
      }
    }
  }   // end loop over K

  // -------- Terminal stage K = N_hat --------
  // Only path constraints (CD) and terminal-cost RSQ (just Qxx + qx).
  {
    casadi_int K2 = p->N_hat;
    casadi_int n = d->nx_hat[K2];
    const casadi_int sQ = p->nx[p->N] + p->nu[p->N];
    const T1* RSQ_full = d->RSQ_val + p->RSQ_offsets[p->N];
    T1 *RSQ_hatN = d->RSQ_hat_val + d->RSQ_hat_offsets[K2];
    T1 *qr_hatN, *CD_hatN;
    casadi_int t_off_qr = 0;
    for (i = 0; i < K2; ++i) t_off_qr += d->nx_hat[i] + d->nu_hat[i];
    qr_hatN = d->qr_hat_val + t_off_qr;
    CD_hatN = d->CD_hat_val + d->CD_hat_offsets[K2];

    // RSQ_hat[N_hat] = Qxx of stage N (top-left n x n block of full RSQ[N])
    for (jj = 0; jj < n; ++jj) {
      for (ii = 0; ii < n; ++ii)
        RSQ_hatN[ii + jj * n] = RSQ_full[ii + jj * sQ];
    }
    // qr_hat[N_hat] = q of stage N (length n)
    for (ii = 0; ii < n; ++ii) qr_hatN[ii] = d->qr_val[off_qr + ii];

    // CD_hat[N_hat] = C of stage N (ng[N] x n)
    {
      casadi_int ng = p->ng[p->N];
      const T1* CN = d->CD_val + p->CD_offsets[p->N];
      for (jj = 0; jj < n; ++jj) {
        for (ii = 0; ii < ng; ++ii)
          CD_hatN[ii + jj * ng] = CN[ii + jj * ng];
      }
      // lbg, ubg: copy through
      for (ii = 0; ii < ng; ++ii) {
        d->lbg_hat_val[casadi_condensing_off_lbg(d, K2) + ii] = d->lbg_val[off_lbg + ii];
        d->ubg_hat_val[casadi_condensing_off_lbg(d, K2) + ii] = d->ubg_val[off_lbg + ii];
      }
    }
    // xi-bounds: copy lbx[N], ubx[N]
    for (ii = 0; ii < n; ++ii) {
      d->lbx_hat_val[casadi_condensing_off_lbx(d, K2) + ii] = d->lbx_val[off_lbx + ii];
      d->ubx_hat_val[casadi_condensing_off_lbx(d, K2) + ii] = d->ubx_val[off_lbx + ii];
    }
  }

  // ---- Pack the condensed-OCP form into the flat-dense Conic-input form
  // (h_hat_csc, a_hat_csc, lbx/ubx, lba/uba).  These outputs
  // are fields of d that are claimed by casadi_condensing_set_work.
  // H_hat is block-diagonal: CSC values are RSQ_hat_val concatenated.
  for (i = 0; i < p->nnz_RSQ_hat; ++i) d->h_hat_csc[i] = d->RSQ_hat_val[i];
  // A_hat CSC: walk z_hat columns, emit per-col contributions in
  // ascending row order (-I from K-1's gap, AB_hat[K], CD_hat[K]).
  {
    casadi_int K, j_local, ai_idx = 0;
    for (K = 0; K <= p->N_hat; ++K) {
      casadi_int nx_K = d->nx_hat[K];
      casadi_int nxu_K = nx_K + d->nu_hat[K];
      casadi_int ng_K = d->ng_hat[K];
      casadi_int nxp1 = (K < p->N_hat) ? d->nx_hat[K+1] : 0;
      for (j_local = 0; j_local < nxu_K; ++j_local) {
        if (K >= 1 && j_local < nx_K) {
          d->a_hat_csc[ai_idx++] = -1;
        }
        if (K < p->N_hat) {
          const T1* col = d->AB_hat_val + d->AB_hat_offsets[K] + j_local * nxp1;
          for (i = 0; i < nxp1; ++i) d->a_hat_csc[ai_idx++] = col[i];
        }
        if (ng_K > 0) {
          const T1* col = d->CD_hat_val + d->CD_hat_offsets[K] + j_local * ng_K;
          for (i = 0; i < ng_K; ++i) d->a_hat_csc[ai_idx++] = col[i];
        }
      }
    }
  }
  // lbx/ubx: interleave xi-bounds and u-bounds per K
  {
    casadi_int K, dst = 0, off_x = 0, off_u = 0;
    for (K = 0; K <= p->N_hat; ++K) {
      casadi_int nxK = d->nx_hat[K];
      casadi_int nuK = d->nu_hat[K];
      for (i = 0; i < nxK; ++i) {
        d->lbx[dst] = d->lbx_hat_val[off_x + i];
        d->ubx[dst] = d->ubx_hat_val[off_x + i];
        ++dst;
      }
      off_x += nxK;
      for (i = 0; i < nuK; ++i) {
        d->lbx[dst] = d->lbu_hat_val[off_u + i];
        d->ubx[dst] = d->ubu_hat_val[off_u + i];
        ++dst;
      }
      off_u += nuK;
    }
  }
  // lba/uba: per K, gap rows (lba=uba=-b_hat) then path rows (lbg/ubg)
  {
    casadi_int K, dst = 0, off_b = 0, off_g = 0;
    for (K = 0; K < p->N_hat; ++K) {
      casadi_int nxp1 = d->nx_hat[K+1];
      for (i = 0; i < nxp1; ++i) {
        d->lba[dst] = -d->b_hat_val[off_b + i];
        d->uba[dst] = -d->b_hat_val[off_b + i];
        ++dst;
      }
      off_b += nxp1;
      casadi_int ngK = d->ng_hat[K];
      for (i = 0; i < ngK; ++i) {
        d->lba[dst] = d->lbg_hat_val[off_g + i];
        d->uba[dst] = d->ubg_hat_val[off_g + i];
        ++dst;
      }
      off_g += ngK;
    }
    casadi_int ngN = d->ng_hat[p->N_hat];
    for (i = 0; i < ngN; ++i) {
      d->lba[dst] = d->lbg_hat_val[off_g + i];
      d->uba[dst] = d->ubg_hat_val[off_g + i];
      ++dst;
    }
  }

  return 0;
}

// SYMBOL "condensing_lift"
//
// Reconstruct the full original-problem primal and dual from a condensed-
// problem solution.  Output layout (canonical interleaved OCP order):
//
//   z_full       : [x_0; u_0; x_1; u_1; ...; x_{N-1}; u_{N-1}; x_N]
//                  total length = sum_{k=0..N} nx[k] + sum_{k=0..N-1} nu[k]
//   lam_x_full   : same layout as z_full -- per-variable bound duals
//   lam_a_full   : [gap_0; gap_1; ...; gap_{N-1}; path_0; path_1; ...]
//                  gap rows length nx[k+1] each, path rows ng[k] each.
//                  Caller layout decision: we emit gap rows for k=0..N-1
//                  followed by path rows for k=0..N (current API supports
//                  ng=0; path-ineq lifting is documented as TODO).
//
// Inputs (read from d):
//   d->x       : condensed primal       (length p->nx_total_hat)
//   d->lam_x   : condensed bound duals  (same length)
//   d->lam_a   : condensed linear-constraint duals (length p->na_total_hat)
//
// Outputs (written into d):
//   d->x_lifted     : full primal       (length p->total_qr)
//   d->lam_x_lifted : full bound duals  (length p->total_qr)
//   d->lam_a_lifted : full linear-constraint duals
//                     (length p->total_b + p->total_g)
//
// Caller copies d->*_lifted into res[CONIC_X/LAM_X/LAM_A] explicitly.
//
// Workspace: requires d->Phi/Gamma/phi/Phi_new/Gamma_new/phi_new
// already allocated by casadi_condensing_set_work.  Lift uses Phi_new (nx_max
// scratch) for the backward Riccati accumulator.
//
// Limitations (v1):
//   - Assumes ng[k] = 0 (no path inequalities).  Path-ineq dual mapping
//     is straightforward but not implemented yet.
template<typename T1>
void casadi_condensing_lift(casadi_condensing_data<T1>* d) {
  const casadi_condensing_prob<T1>* p = d->prob;
  const T1* zhat = d->x;
  const T1* lam_x_hat = d->lam_x;
  const T1* lam_a_hat = d->lam_a;
  T1* z_full = d->x_lifted;
  T1* lam_x_full = d->lam_x_lifted;
  T1* lam_a_full = d->lam_a_lifted;
  casadi_int K, j, k, k_a, k_b, M_blk, i, jj, c, kk;
  casadi_int xi_off, nu_off;     // walking offsets in zhat / lam_x_hat
  casadi_int row_cur;            // walking offset in lam_a_hat
  casadi_int z_off, gap_off;     // walking offsets in z_full / lam_a_full

  // ---- Pass 1: primal lift (reconstruct x_k, u_k by walking dynamics) ----
  // Initialise running offsets through the interleaved layouts.
  z_off = 0;
  xi_off = 0;
  nu_off = 0;
  for (K = 0; K <= p->N_hat; ++K) {
    casadi_int nxK = d->nx_hat[K];
    nu_off += nxK;   // skip xi_K block to position at nu_K start
    // (nu_off updated below per K)
    if (K < p->N_hat) nu_off += d->nu_hat[K];
  }
  // The above approach is too tangled; recompute cleanly per K below.

  // For each block K, x_{M[K]} = xi_K, then walk dynamics with the
  // controls slice of nu_K to get x_{M[K]+1}, ..., x_{M[K+1]}.
  z_off = 0;
  xi_off = 0;   // running offset in zhat for current xi_K
  nu_off = 0;   // running offset in zhat for current nu_K
  for (K = 0; K < p->N_hat; ++K) {
    k_a = p->M[K];
    k_b = p->M[K + 1];
    M_blk = k_b - k_a;
    casadi_int nxK = d->nx_hat[K];

    // xi_K is at zhat[xi_off ... xi_off + nxK)
    // Walk to find positions: xi_off updated outside K loop below.
    // Compute it explicitly here:
    {
      casadi_int o = 0, KK;
      for (KK = 0; KK < K; ++KK) o += d->nx_hat[KK] + d->nu_hat[KK];
      xi_off = o;
      nu_off = o + nxK;
    }
    // Write x_{k_a} = xi_K into z_full at z_off
    for (i = 0; i < nxK; ++i) z_full[z_off + i] = zhat[xi_off + i];

    // Walk inner stages: append u_k followed by x_{k+1}
    casadi_int u_within = 0;
    for (j = 0; j < M_blk; ++j) {
      k = k_a + j;
      casadi_int nx_k = p->nx[k];
      casadi_int nu_k = p->nu[k];
      casadi_int nx_kp1 = p->nx[k + 1];
      casadi_int z_x_off = z_off;     // x_k is at z_off
      casadi_int z_u_off = z_off + nx_k;   // u_k follows
      casadi_int z_x1_off = z_off + nx_k + nu_k;  // x_{k+1}

      // u_k slice from nu_K
      for (i = 0; i < nu_k; ++i)
        z_full[z_u_off + i] = zhat[nu_off + u_within + i];
      u_within += nu_k;

      // x_{k+1} = A_k x_k + B_k u_k + b_k
      const T1* AB = d->AB_val + p->AB_offsets[k];   // nx_kp1 x (nx_k + nu_k), col-major
      const T1* b_k = d->b_val;
      // Compute b_k offset (sum nx[k'+1] for k'<k)
      casadi_int b_pos = 0;
      for (kk = 0; kk < k; ++kk) b_pos += p->nx[kk + 1];
      b_k += b_pos;

      for (i = 0; i < nx_kp1; ++i) {
        T1 s = b_k[i];
        for (jj = 0; jj < nx_k; ++jj)
          s += AB[i + jj * nx_kp1] * z_full[z_x_off + jj];
        for (jj = 0; jj < nu_k; ++jj)
          s += AB[i + (nx_k + jj) * nx_kp1] * z_full[z_u_off + jj];
        z_full[z_x1_off + i] = s;
      }

      // Advance z_off to past u_k (next iter's z_x_off = z_x1_off)
      z_off = z_x1_off;
    }
    // After block K, z_off points at x_{k_b} which equals xi_{K+1}.
  }
  // Terminal block already populated by walking the last block above.

  // ---- Pass 2: lift bound duals (lam_x) ----
  // (a) Boundary state duals + control duals -- direct mappings.
  z_off = 0;
  for (K = 0; K <= p->N_hat; ++K) {
    k_a = p->M[K];
    casadi_int nxK = d->nx_hat[K];
    // xi_off / nu_off in zhat / lam_x_hat at start of block K
    casadi_int o = 0, KK;
    for (KK = 0; KK < K; ++KK) o += d->nx_hat[KK] + d->nu_hat[KK];
    xi_off = o;
    nu_off = o + nxK;

    // Boundary state bound dual: lam_x_full[x_{k_a}] = lam_x_hat[xi_K]
    // Find z_off for x_{k_a}: sum nx[k'] + nu[k'] for k' < k_a, plus nx_a partial
    casadi_int z_xa = 0;
    for (kk = 0; kk < k_a; ++kk) z_xa += p->nx[kk] + p->nu[kk];
    for (i = 0; i < nxK; ++i)
      lam_x_full[z_xa + i] = lam_x_hat[xi_off + i];

    // Control bound duals for stages [k_a, k_b)
    if (K < p->N_hat) {
      casadi_int u_within = 0;
      k_b = p->M[K + 1];
      for (j = 0; j < k_b - k_a; ++j) {
        k = k_a + j;
        casadi_int nu_k = p->nu[k];
        // z position of u_k = (sum nx[k'] + nu[k'] for k'<k) + nx[k]
        casadi_int z_uk = p->nx[k];
        for (kk = 0; kk < k; ++kk) z_uk += p->nx[kk] + p->nu[kk];
        for (i = 0; i < nu_k; ++i)
          lam_x_full[z_uk + i] = lam_x_hat[nu_off + u_within + i];
        u_within += nu_k;
      }
    }
  }

  // ---- Pass 3: linear-constraint duals (lam_a) and inner-state bound duals ----
  // Walk lam_a_hat: per K, gap rows (nx_hat[K+1]) then lifted-bound rows
  // (ng_hat[K]).
  row_cur = 0;
  gap_off = 0;
  for (K = 0; K < p->N_hat; ++K) {
    k_a = p->M[K];
    k_b = p->M[K + 1];
    M_blk = k_b - k_a;
    casadi_int nx_kp1 = p->nx[k_b];

    // Last gap-closing dual in this block: lam_a_full at row offset
    // sum_{k'<k_b-1} nx[k'+1] (== sum_{k'=0..k_b-2} nx[k'+1])
    casadi_int gap_last_off = 0;
    for (kk = 0; kk < k_b - 1; ++kk) gap_last_off += p->nx[kk + 1];
    for (i = 0; i < nx_kp1; ++i)
      lam_a_full[gap_last_off + i] = lam_a_hat[row_cur + i];
    row_cur += nx_kp1;

    // Lifted-bound duals for inner stages (j=1..M-1): map to lam_x_full
    // at x_{k_a+j} positions.  Each inner stage contributes nx[k_a+j] rows.
    casadi_int lift_within = 0;
    for (j = 1; j < M_blk; ++j) {
      k = k_a + j;
      casadi_int nx_k = p->nx[k];
      casadi_int z_xk = 0;
      for (kk = 0; kk < k; ++kk) z_xk += p->nx[kk] + p->nu[kk];
      for (i = 0; i < nx_k; ++i)
        lam_x_full[z_xk + i] = lam_a_hat[row_cur + lift_within + i];
      lift_within += nx_k;
    }
    row_cur += d->ng_hat[K];

    // Backward Riccati: recover lam^gap_{k-1} for k = k_b-1, k_b-2, ..., k_a+1
    //   lam^gap_{k-1} = A_k' lam^gap_k - R_xx x_k - R_xu u_k - q_k - mu^x_k
    // where mu^x_k is the bound dual for x_k (already lifted above).
    // Use d->Phi_new as scratch for lam_k (size nx_max).
    if (M_blk > 1) {
      T1* lam_k = d->Phi_new;          // borrow (nx_max scratch)
      // Initialise lam_k = lam^gap_{k_b - 1}
      for (i = 0; i < nx_kp1; ++i) lam_k[i] = lam_a_full[gap_last_off + i];

      for (j = M_blk - 1; j >= 1; --j) {
        k = k_a + j;
        casadi_int nx_k = p->nx[k];
        casadi_int nu_k = p->nu[k];
        casadi_int sQ = nx_k + nu_k;

        // Read original-stage data
        const T1* AB = d->AB_val + p->AB_offsets[k];        // nx[k+1] x (nx_k+nu_k)
        const T1* RSQ = d->RSQ_val + p->RSQ_offsets[k];     // sQ x sQ
        const T1* qr = d->qr_val;
        casadi_int qr_pos = 0;
        for (kk = 0; kk < k; ++kk) qr_pos += p->nx[kk] + p->nu[kk];
        const T1* q_k = qr + qr_pos;                        // length nx_k

        // x_k, u_k positions in z_full
        casadi_int z_xk = 0;
        for (kk = 0; kk < k; ++kk) z_xk += p->nx[kk] + p->nu[kk];
        const T1* x_k = z_full + z_xk;
        const T1* u_k = z_full + z_xk + nx_k;
        const T1* mu_x_k = lam_x_full + z_xk;   // bound dual for x_k

        // Compute lam_km1 = A_k' lam_k - R_xx*x_k - R_xu*u_k - q_k - mu_x_k
        // lam_km1 has length nx[k] (= nx_k assuming nx[k_a..k_b-1] == nx_k for inner stages)
        // Reuse Gamma_new's first nx_max entries as scratch for lam_km1
        T1* lam_km1 = d->Gamma_new;

        // A_k is nx[k+1] x nx_k.  A_k' has shape nx_k x nx[k+1].
        // (A_k' lam_k)[i] = sum_p A_k[p,i] * lam_k[p] = sum_p AB[p + i*nx[k+1]] * lam_k[p]
        casadi_int nx_kp_loc = p->nx[k + 1];
        for (i = 0; i < nx_k; ++i) {
          T1 s = 0;
          for (c = 0; c < nx_kp_loc; ++c)
            s += AB[c + i * nx_kp_loc] * lam_k[c];
          lam_km1[i] = s;
        }
        // CasADi dual convention: H z + g + lam_x + A^T lam_a = 0 at the
        // optimum.  For an interior stage, the components for x_k yield
        //   lam_gap_(k-1) = R_xx x_k + R_xu u_k + q_k + lam_x_k + A_k^T lam_gap_k
        // so all three additional terms are added (NOT subtracted).
        for (i = 0; i < nx_k; ++i) {
          T1 s = 0;
          for (c = 0; c < nx_k; ++c) s += RSQ[i + c * sQ] * x_k[c];
          lam_km1[i] += s;
        }
        for (i = 0; i < nx_k; ++i) {
          T1 s = 0;
          for (c = 0; c < nu_k; ++c) s += RSQ[i + (nx_k + c) * sQ] * u_k[c];
          lam_km1[i] += s;
        }
        for (i = 0; i < nx_k; ++i) lam_km1[i] += q_k[i] + mu_x_k[i];

        // Write into lam_a_full at gap row k-1 offset
        casadi_int gap_km1_off = 0;
        for (kk = 0; kk < k - 1; ++kk) gap_km1_off += p->nx[kk + 1];
        for (i = 0; i < nx_k; ++i) lam_a_full[gap_km1_off + i] = lam_km1[i];

        // Advance: lam_k <- lam_km1 for next iteration (j-1)
        for (i = 0; i < nx_k; ++i) lam_k[i] = lam_km1[i];
      }
    }
  }

  // (gap_off unused; suppress warning)
  (void)gap_off;
}

// Inline offset helpers: cumulative sums of condensed-stage dimensions.
// Read from data (the *_hat arrays live there now).
template<typename T1>
static casadi_int casadi_condensing_off_lbx(const casadi_condensing_data<T1>* d,
                                            casadi_int K) {
  casadi_int s = 0, i;
  for (i = 0; i < K; ++i) s += d->nx_hat[i];
  return s;
}
template<typename T1>
static casadi_int casadi_condensing_off_lbu(const casadi_condensing_data<T1>* d,
                                            casadi_int K) {
  casadi_int s = 0, i;
  for (i = 0; i < K; ++i) s += d->nu_hat[i];
  return s;
}
template<typename T1>
static casadi_int casadi_condensing_off_lbg(const casadi_condensing_data<T1>* d,
                                            casadi_int K) {
  casadi_int s = 0, i;
  for (i = 0; i < K; ++i) s += d->ng_hat[i];
  return s;
}

// (extract_g / extract_bounds / extract_a_rhs have been folded into
//  casadi_condensing_eval -- it now takes the user g/lbx/ubx/lba/uba
//  pointers and demultiplexes them internally.  assemble_a / pack_lbx
//  / pack_lba are also folded in -- the flat-dense Conic-input form
//  lands in d->h_hat_csc / a_hat_csc / lbx/ubx / lba/
//  uba as part of eval.)

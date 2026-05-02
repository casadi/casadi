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

// Partial condensing eliminates states inside each partition interval [M[K], M[K+1]).
// Stage data is column-major: AB = [A B], CD = [C D], RSQ = [Q S; S' R], qr = [q; r].
// Condensed stages retain boundary states and concatenate controls and lifted state bounds.

// C-REPLACE "casadi_blas_mtimes" "casadi_mtimes_dense"
// C-REPLACE "casadi_ocp_block" "struct casadi_ocp_block"
// C-REPLACE "casadi_condensing_prob<T1>" "struct casadi_condensing_prob"
// C-REPLACE "casadi_condensing_data<T1>" "struct casadi_condensing_data"
// C-REPLACE "static_cast<int>" "(int) "
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


  // Maxima for workspace sizing (filled by condensing_setup)
  casadi_int nx_max;          // max over k of nx[k]
  casadi_int nu_max_block;    // max over K of nu_hat[K]
  casadi_int ng_max, sz_gemm, sz_v;
  const casadi_int* gap_nz;
  casadi_int nxu_max_block;   // max over K of (nx_hat[K] + nu_hat[K])

  // Cumulative sizes computed by setup
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
template<typename T1>
void casadi_condensing_setup(casadi_condensing_prob<T1>* p) {
  casadi_int K, j, k_a, k_b, M_, sum_nu, sum_ng_path, sum_ng_lift;
  casadi_int nx_K, nxu_K;
  casadi_int off_AB = 0, off_CD = 0, off_RSQ = 0;

  p->nx_max = 0;
  p->nu_max_block = 0;
  p->nxu_max_block = 0;

  p->ng_max = 0;
  for (K = 0; K <= p->N; ++K) {
    if (p->nx[K] > p->nx_max) p->nx_max = p->nx[K];
    if (p->ng[K] > p->ng_max) p->ng_max = p->ng[K];
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
  p->sz_gemm = p->nx_max * p->nu_max_block;
  if (p->ng_max * p->nx_max > p->sz_gemm) p->sz_gemm = p->ng_max * p->nx_max;
  if (p->ng_max * p->nu_max_block > p->sz_gemm)
    p->sz_gemm = p->ng_max * p->nu_max_block;
  p->sz_v = p->nx_max + p->nu_max_block;
  if (p->ng_max > p->sz_v) p->sz_v = p->ng_max;
}

// SYMBOL "condensing_data"
template<typename T1>
struct casadi_condensing_data {
  const casadi_condensing_prob<T1> *prob;
  T1 cost;

  // Condensed dimensions and block descriptors, populated by set_work
  casadi_int *nx_hat, *nu_hat, *ng_hat;
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

  // Condensed solver inputs and outputs
  T1 *h_hat_csc;     // length p->nnz_RSQ_hat
  T1 *a_hat_csc;     // length p->nnz_a_hat_csc
  T1 *lbx, *ubx;  // length p->nx_total_hat
  T1 *lba, *uba;  // length p->na_total_hat
  T1 *x;
  T1 *lam_x;
  T1 *lam_a;

  // Lifted primal and dual solution
  T1 *x_lifted;       // length p->total_qr
  T1 *lam_x_lifted;   // length p->total_qr
  T1 *lam_a_lifted;   // length p->total_b + p->total_g

  // Original inputs, assigned by the caller before eval
  const T1 *x0_orig, *lam_x0_orig, *lam_a0_orig;
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

// Cumulative offsets into condensed-stage vectors
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
template<typename T1>
void casadi_condensing_work(const casadi_condensing_prob<T1>* p,
                            casadi_int* sz_iw, casadi_int* sz_w) {
  casadi_int nx2 = p->nx_max * p->nx_max;
  casadi_int nxnu = p->nx_max * p->nu_max_block;
  casadi_int nu2 = p->nu_max_block * p->nu_max_block;
  // Condensed dimensions, offsets and block descriptors
  *sz_iw += 3 * (p->N_hat + 1);          // nx_hat, nu_hat, ng_hat
  *sz_iw += 3 * (p->N_hat + 1);          // *_hat_offsets
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
  *sz_w += p->sz_gemm;                 // gemm_xu
  *sz_w += p->sz_v;
  // Condensed CSC + bounds + primal/dual scratch
  *sz_w += p->nnz_RSQ_hat;             // h_hat_csc
  *sz_w += p->nnz_a_hat_csc;
  *sz_w += 2 * p->nx_total_hat;        // lbx, ubx
  *sz_w += 2 * p->na_total_hat;        // lba, uba
  *sz_w += p->nx_total_hat;            // x
  *sz_w += p->nx_total_hat;            // lam_x
  *sz_w += p->na_total_hat;            // lam_a
  *sz_w += 2 * p->total_qr + p->total_b + p->total_g;
}

// SYMBOL "condensing_set_work"
template<typename T1>
void casadi_condensing_set_work(casadi_condensing_data<T1>* d,
    const T1*** arg, T1*** res, casadi_int** iw, T1** w) {
  (void)arg; (void)res;
  const casadi_condensing_prob<T1>* p = d->prob;
  casadi_int nx2 = p->nx_max * p->nx_max;
  casadi_int nxnu = p->nx_max * p->nu_max_block;
  casadi_int nu2 = p->nu_max_block * p->nu_max_block;

  // Condensed dimensions and offsets
  d->nx_hat = *iw; *iw += p->N_hat + 1;
  d->nu_hat = *iw; *iw += p->N_hat + 1;
  d->ng_hat = *iw; *iw += p->N_hat + 1;
  d->AB_hat_offsets  = *iw; *iw += p->N_hat + 1;
  d->CD_hat_offsets  = *iw; *iw += p->N_hat + 1;
  d->RSQ_hat_offsets = *iw; *iw += p->N_hat + 1;
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
      d->AB_hat_offsets[K] = off_AB;
      off_AB += d->nx_hat[K+1] * nxu_K;
      d->CD_hat_offsets[K] = off_CD;
      off_CD += d->ng_hat[K] * nxu_K;
      d->RSQ_hat_offsets[K] = off_RSQ;
      off_RSQ += nxu_K * nxu_K;
    }
    // Terminal stage K = N_hat: only path constraints.
    d->nu_hat[p->N_hat] = 0;
    d->ng_hat[p->N_hat] = p->ng[p->N];
    nx_K = d->nx_hat[p->N_hat];
    d->AB_hat_offsets[p->N_hat] = off_AB;
    d->CD_hat_offsets[p->N_hat] = off_CD;
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
  d->gemm_xu   = *w; *w += p->sz_gemm;
  d->tmp_v     = *w; *w += p->sz_v;

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

  d->x_lifted = *w; *w += p->total_qr;
  d->lam_x_lifted = *w; *w += p->total_qr;
  d->lam_a_lifted = *w; *w += p->total_b + p->total_g;
}


// SYMBOL "condensing_check"
template<typename T1>
int casadi_condensing_check(const casadi_condensing_prob<T1>* p,
    const T1* a, const T1* lba, const T1* uba) {
  casadi_int k, i, nz = 0;
  for (k = 0; k < p->N; ++k) {
    for (i = 0; i < p->nx[k+1]; ++i, ++nz) {
      casadi_int row = p->AB[k].offset_r + i;
      T1 lb = lba ? lba[row] : -std::numeric_limits<T1>::infinity();
      T1 ub = uba ? uba[row] : std::numeric_limits<T1>::infinity();
      if (!(lb == ub && lb > -std::numeric_limits<T1>::infinity()
          && ub < std::numeric_limits<T1>::infinity())) return 1;
      if (!a || a[p->gap_nz[nz]] != -1) return 2;
    }
  }
  return 0;
}

// SYMBOL "condensing_eval"
// Project H and A into RSQ_val, AB_val and CD_val before calling.
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
    for (ii = 0; ii < total_qr; ++ii) d->qr_val[ii] = d->g_orig ? d->g_orig[ii] : 0;
    // lbx/ubx_orig are interleaved [x_k; u_k]; demultiplex.
    for (kk = 0; kk <= p->N; ++kk) {
      for (ii = 0; ii < p->nx[kk]; ++ii) {
        d->lbx_val[dx + ii] = (d->lbx_orig ? d->lbx_orig[src + ii] : -std::numeric_limits<T1>::infinity());
        d->ubx_val[dx + ii] = (d->ubx_orig ? d->ubx_orig[src + ii] : std::numeric_limits<T1>::infinity());
      }
      src += p->nx[kk];  dx  += p->nx[kk];
      for (ii = 0; ii < p->nu[kk]; ++ii) {
        d->lbu_val[du + ii] = (d->lbx_orig ? d->lbx_orig[src + ii] : -std::numeric_limits<T1>::infinity());
        d->ubu_val[du + ii] = (d->ubx_orig ? d->ubx_orig[src + ii] : std::numeric_limits<T1>::infinity());
      }
      src += p->nu[kk];  du  += p->nu[kk];
    }
    // lba/uba_orig: per-stage (nx[k+1] gap + ng[k] path), terminal ng[N].
    for (kk = 0; kk < p->N; ++kk) {
      for (ii = 0; ii < p->nx[kk + 1]; ++ii) {
        d->b_val[ob + ii] = -(d->lba_orig ? d->lba_orig[sa + ii] : -std::numeric_limits<T1>::infinity());
      }
      sa += p->nx[kk + 1];  ob += p->nx[kk + 1];
      for (ii = 0; ii < p->ng[kk]; ++ii) {
        d->lbg_val[og + ii] = (d->lba_orig ? d->lba_orig[sa + ii] : -std::numeric_limits<T1>::infinity());
        d->ubg_val[og + ii] = (d->uba_orig ? d->uba_orig[sa + ii] : std::numeric_limits<T1>::infinity());
      }
      sa += p->ng[kk];  og += p->ng[kk];
    }
    for (ii = 0; ii < p->ng[p->N]; ++ii) {
      d->lbg_val[og + ii] = (d->lba_orig ? d->lba_orig[sa + ii] : -std::numeric_limits<T1>::infinity());
      d->ubg_val[og + ii] = (d->uba_orig ? d->uba_orig[sa + ii] : std::numeric_limits<T1>::infinity());
    }
  }
  casadi_int row_off;                          // running row offset within CD_hat
  T1 *Phi, *Gamma, *phi, *phi_new;
  const T1 *A, *B, *q_k, *r_k, *C, *D, *b_k;
  T1 *RSQ_hat, *Hxx, *Hxu, *Huu, *hx, *hu, *AB_hat, *CD_hat, *qr_hat, *b_hat;
  T1 *Phi_out, *Gamma_out;
  casadi_int i, jj, ii, kk;

  // Running offsets into the original stage vectors
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
    casadi_clear(Phi, nx_K * nx_K);
    for (i = 0; i < nx_K; ++i) Phi[i + i * nx_K] = 1;
    casadi_clear(phi, nx_K);
    // Gamma: no columns at j=0, nothing to initialise.

    // Output buffers for this condensed block
    AB_hat  = d->AB_hat_val  + d->AB_hat_offsets[K];
    CD_hat  = d->CD_hat_val  + d->CD_hat_offsets[K];
    RSQ_hat = d->RSQ_hat_val + d->RSQ_hat_offsets[K];
    qr_hat  = d->qr_hat_val;        // walked stage-by-stage; offset below
    b_hat   = d->b_hat_val;         // walked likewise

    // Offsets into condensed gradient and dynamics vectors
    {
      casadi_int t_off_qr = 0, t_off_b = 0;
      for (i = 0; i < K; ++i) {
        t_off_qr += d->nx_hat[i] + d->nu_hat[i];
        t_off_b  += d->nx_hat[i + 1];
      }
      qr_hat = d->qr_hat_val + t_off_qr;
      b_hat  = d->b_hat_val  + t_off_b;
    }

    // Hessian subblocks use column stride nxu_K.
    casadi_clear(RSQ_hat, nxu_K * nxu_K);
    casadi_clear(qr_hat, nxu_K);
    Hxx = RSQ_hat;                                  // top-left  nx_K x nx_K
    Hxu = RSQ_hat + nx_K * nxu_K;                   // top-right (in xu cols)
    Huu = RSQ_hat + nx_K * nxu_K + nx_K;            // bottom-right
    hx  = qr_hat;
    hu  = qr_hat + nx_K;

    // CD_hat: zero whole block; we write contiguous row groups
    casadi_clear(CD_hat, d->ng_hat[K] * nxu_K);

    // Copy boundary-state bounds; collect the remaining bounds below.
    casadi_copy(d->lbx_val + off_lbx, nx_K, d->lbx_hat_val + casadi_condensing_off_lbx(d, K));
    casadi_copy(d->ubx_val + off_lbx, nx_K, d->ubx_hat_val + casadi_condensing_off_lbx(d, K));
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
      C   = d->CD_val  + p->CD_offsets[k];                     // ng_k x nx_k
      D   = C + ng_k * nx_k;                                   // ng_k x nu_k
      q_k = d->qr_val + off_qr;                                // length nx_k
      r_k = q_k + nx_k;                                        // length nu_k
      b_k = d->b_val + off_b;                                  // length nx_kp1

      // Extract contiguous Hessian subblocks for BLAS.
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

        // Substitute x_k = Phi*xi + Gamma*u + phi into the stage cost.
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

        // Hxx += Phi' Qxx Phi; write strided output one column at a time.
        {
          T1 *tmp = d->gemm_xx;       // (nx_k x nx_K)
          casadi_clear(tmp, nx_k * nx_K);
          casadi_blas_mtimes(Qxx_c, nx_k, nx_k, Phi, nx_K, tmp, 0);
          for (kk = 0; kk < nx_K; ++kk) {
            T1 *hxx_col = Hxx + kk * nxu_K;
            casadi_blas_mtimes(Phi, nx_k, nx_K, tmp + kk * nx_k, 1, hxx_col, 1);
          }
        }

        // Hxu += Phi' [Qxx Gamma, Qxu]
        {
          T1 *tmp = d->gemm_xu;       // (nx_k x cum_nu) when used
          if (cum_nu > 0) {
            casadi_clear(tmp, nx_k * cum_nu);
            casadi_blas_mtimes(Qxx_c, nx_k, nx_k, Gamma, cum_nu, tmp, 0);
            // Write Hxu one column at a time (stride nxu_K).
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

        // Huu += [Gamma' Qxx Gamma, Gamma' Qxu; Qxu' Gamma, Quu]
        if (cum_nu > 0) {
          T1 *tmp = d->gemm_xu;
          casadi_clear(tmp, nx_k * cum_nu);
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

        // Mirror the off-diagonal control Hessian block.
        if (cum_nu > 0) {
          for (kk = 0; kk < nu_k; ++kk) {        // dst col (cum_nu+kk)
            for (ii = 0; ii < cum_nu; ++ii) {
              T1 v = Huu[ii + (cum_nu + kk) * nxu_K];
              Huu[(cum_nu + kk) + ii * nxu_K] = v;
            }
          }
        }
      }

      // Substitute the state propagation into the path constraints.
      if (ng_k > 0) {
        // Compute C_k * phi (length ng_k)
        T1 *Cphi = d->tmp_v;
        casadi_clear(Cphi, ng_k);
        casadi_blas_mtimes(C, ng_k, nx_k, phi, 1, Cphi, 0);

        // State columns: C_k * Phi
        {
          T1 *tmp = d->gemm_xu;
          casadi_clear(tmp, ng_k * nx_K);
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
          casadi_clear(tmp, ng_k * cum_nu);
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

      // Lift interior-state bounds: lbx - phi <= Phi*xi + Gamma*u <= ubx - phi.
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
        // Shift the bounds by the affine state offset.
        for (ii = 0; ii < nx_k; ++ii) {
          d->lbg_hat_val[casadi_condensing_off_lbg(d, K) + off_lbg_block + ii]
              = d->lbx_val[off_lbx + ii] - phi[ii];
          d->ubg_hat_val[casadi_condensing_off_lbg(d, K) + off_lbg_block + ii]
              = d->ubx_val[off_lbx + ii] - phi[ii];
        }
        off_lbg_block += nx_k;
        row_off += nx_k;
      }

      // Control bounds: append lbu_k, ubu_k
      for (ii = 0; ii < nu_k; ++ii) {
        d->lbu_hat_val[casadi_condensing_off_lbu(d, K) + off_lbu_block + ii]
            = d->lbu_val[off_lbu + ii];
        d->ubu_hat_val[casadi_condensing_off_lbu(d, K) + off_lbu_block + ii]
            = d->ubu_val[off_lbu + ii];
      }
      off_lbu_block += nu_k;

      // Alternate propagation buffers to avoid aliasing BLAS inputs and outputs.
      Phi_out   = (Phi   == d->Phi)   ? d->Phi_new   : d->Phi;
      Gamma_out = (Gamma == d->Gamma) ? d->Gamma_new : d->Gamma;
      phi_new   = (phi   == d->phi)   ? d->phi_new   : d->phi;

      // phi_new = A * phi + b_k
      casadi_copy(b_k, nx_kp1, phi_new);
      casadi_blas_mtimes(A, nx_kp1, nx_k, phi, 1, phi_new, 0);

      // Phi_new = A * Phi  (nx_kp1 x nx_K)
      casadi_clear(Phi_out, nx_kp1 * nx_K);
      casadi_blas_mtimes(A, nx_kp1, nx_k, Phi, nx_K, Phi_out, 0);

      // Gamma_new = [A * Gamma   B_k]   (nx_kp1 x (cum_nu + nu_k))
      casadi_clear(Gamma_out, nx_kp1 * (cum_nu + nu_k));
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

      // Advance original-stage vector offsets
      off_b  += nx_kp1;
      off_qr += nx_k + nu_k;
      off_lbx += nx_k;
      off_lbu += nu_k;
      off_lbg += ng_k;
      cum_nu += nu_k;
    }   // end inner loop over j

    // Store the condensed dynamics AB_hat = [Phi Gamma], b_hat = phi.
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

    // Mirror Hxu into Hux.
    for (jj = 0; jj < nu_K; ++jj) {
      for (ii = 0; ii < nx_K; ++ii) {
        T1 v = RSQ_hat[ii + (nx_K + jj) * nxu_K];   // Hxu[ii, jj]
        RSQ_hat[(nx_K + jj) + ii * nxu_K] = v;       // Hux[jj, ii]
      }
    }
  }   // end loop over K

  // Copy terminal cost, constraints and bounds.
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

  // The block-diagonal Hessian is already in CSC nonzero order.
  for (i = 0; i < p->nnz_RSQ_hat; ++i) d->h_hat_csc[i] = d->RSQ_hat_val[i];
  // Pack A columns in row order: previous gap, current dynamics, path constraints.
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

  // Map warm starts into the condensed variables and constraints.
  {
    casadi_int xhat = 0, ahat = 0;
    for (K = 0; K <= p->N_hat; ++K) {
      k_a = p->M[K];
      k_b = K < p->N_hat ? p->M[K+1] : p->N;
      casadi_int x = p->RSQ[k_a].offset_r;
      for (i = 0; i < p->nx[k_a]; ++i, ++xhat) {
        d->x[xhat] = d->x0_orig ? d->x0_orig[x+i] : 0;
        d->lam_x[xhat] = d->lam_x0_orig ? d->lam_x0_orig[x+i] : 0;
      }
      if (K < p->N_hat) {
        casadi_int row = p->AB[k_b-1].offset_r;
        for (i = 0; i < p->nx[k_b]; ++i, ++ahat)
          d->lam_a[ahat] = d->lam_a0_orig ? d->lam_a0_orig[row+i] : 0;
      }
      for (k = k_a; k < k_b || (K == p->N_hat && k == k_b); ++k) {
        x = p->RSQ[k].offset_r;
        for (i = 0; i < p->nu[k]; ++i, ++xhat) {
          d->x[xhat] = d->x0_orig ? d->x0_orig[x+p->nx[k]+i] : 0;
          d->lam_x[xhat] = d->lam_x0_orig ? d->lam_x0_orig[x+p->nx[k]+i] : 0;
        }
        casadi_int row = p->CD[k].offset_r;
        for (i = 0; i < p->ng[k]; ++i, ++ahat)
          d->lam_a[ahat] = d->lam_a0_orig ? d->lam_a0_orig[row+i] : 0;
        if (k > k_a) {
          for (i = 0; i < p->nx[k]; ++i, ++ahat)
            d->lam_a[ahat] = d->lam_x0_orig ? d->lam_x0_orig[x+i] : 0;
        }
      }
    }
  }

  return 0;
}

// SYMBOL "condensing_lift"
template<typename T1>
void casadi_condensing_lift(casadi_condensing_data<T1>* d) {
  const casadi_condensing_prob<T1>* p = d->prob;
  casadi_int K, k, i, j, ka, kb, xhat = 0, ahat = 0, b = 0;
  casadi_clear(d->lam_x_lifted, p->total_qr);
  casadi_clear(d->lam_a_lifted, p->total_b + p->total_g);
  for (K = 0; K <= p->N_hat; ++K) {
    ka = p->M[K];
    kb = K < p->N_hat ? p->M[K+1] : p->N;
    casadi_int x = p->RSQ[ka].offset_r;
    casadi_int uhat = xhat + p->nx[ka];
    casadi_copy(d->x + xhat, p->nx[ka], d->x_lifted + x);
    casadi_copy(d->lam_x + xhat, p->nx[ka], d->lam_x_lifted + x);
    xhat += d->nx_hat[K] + d->nu_hat[K];
    if (K < p->N_hat) {
      casadi_copy(d->lam_a + ahat, p->nx[kb],
        d->lam_a_lifted + p->AB[kb-1].offset_r);
      ahat += p->nx[kb];
    }
    // Condensed path rows contain each stage's paths followed by its lifted state bounds.
    for (k = ka; k < kb || (K == p->N_hat && k == kb); ++k) {
      x = p->RSQ[k].offset_r;
      casadi_copy(d->lam_a + ahat, p->ng[k], d->lam_a_lifted + p->CD[k].offset_r);
      ahat += p->ng[k];
      if (k > ka) {
        casadi_copy(d->lam_a + ahat, p->nx[k], d->lam_x_lifted + x);
        ahat += p->nx[k];
      }
      if (k == p->N) break;
      casadi_copy(d->x + uhat, p->nu[k], d->x_lifted + x + p->nx[k]);
      casadi_copy(d->lam_x + uhat, p->nu[k], d->lam_x_lifted + x + p->nx[k]);
      uhat += p->nu[k];
      const T1* AB = d->AB_val + p->AB_offsets[k];
      T1* xn = d->x_lifted + p->RSQ[k+1].offset_r;
      for (i = 0; i < p->nx[k+1]; ++i) {
        T1 v = d->b_val[b++];
        for (j = 0; j < p->nx[k] + p->nu[k]; ++j)
          v += AB[i + j * p->nx[k+1]] * d->x_lifted[x+j];
        xn[i] = v;
      }
    }
    // Recover interior dynamics multipliers from original-problem stationarity.
    for (k = kb-1; k > ka; --k) {
      x = p->RSQ[k].offset_r;
      const T1* AB = d->AB_val + p->AB_offsets[k];
      const T1* CD = d->CD_val + p->CD_offsets[k];
      const T1* H = d->RSQ_val + p->RSQ_offsets[k];
      casadi_int n = p->nx[k] + p->nu[k];
      for (i = 0; i < p->nx[k]; ++i) {
        T1 v = d->qr_val[x+i] + d->lam_x_lifted[x+i];
        for (j = 0; j < n; ++j) v += H[i+j*n] * d->x_lifted[x+j];
        for (j = 0; j < p->nx[k+1]; ++j)
          v += AB[j+i*p->nx[k+1]] * d->lam_a_lifted[p->AB[k].offset_r+j];
        for (j = 0; j < p->ng[k]; ++j)
          v += CD[j+i*p->ng[k]] * d->lam_a_lifted[p->CD[k].offset_r+j];
        d->lam_a_lifted[p->AB[k-1].offset_r+i] = v;
      }
    }
  }
  // Evaluate the original objective, including the affine-elimination constant.
  d->cost = 0;
  for (k = 0; k <= p->N; ++k) {
    casadi_int n = p->nx[k] + p->nu[k];
    const T1* x = d->x_lifted + p->RSQ[k].offset_r;
    const T1* H = d->RSQ_val + p->RSQ_offsets[k];
    const T1* g = d->qr_val + p->RSQ[k].offset_r;
    for (i = 0; i < n; ++i) {
      T1 v = 0;
      for (j = 0; j < n; ++j) v += H[i+j*n] * x[j];
      d->cost += x[i] * (g[i] + 0.5*v);
    }
  }
}

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

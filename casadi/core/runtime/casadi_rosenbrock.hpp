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

// Rosenbrock-Wanner integrator for M y' = F(t, y), see casadi_ivp.hpp
// Step: (M/(h gamma) - J) U_i = F(t0 + alpha_i h, y0 + sum_j a_ij U_j) + sum_j (c_ij/h) M U_j
//         + gamma_i h dF/dt, y1 = y0 + sum_i m_i U_i, error sum_i m_err_i U_i
// Steps depend on the nominal entries only: sensitivities are exact derivatives of the scheme

// C-REPLACE "fabs" "casadi_fabs"
// C-REPLACE "static_cast<T1>" "(casadi_real) "
// C-REPLACE "static_cast<casadi_int*>(0)" "0"
// C-REPLACE "casadi_ivp_prob<T1>" "struct casadi_ivp_prob"
// C-REPLACE "casadi_ivp_data<T1>" "struct casadi_ivp_data"

// SYMBOL "rosenbrock_prob"
template<typename T1>
struct casadi_rosenbrock_prob {
  // Shared problem structure
  const casadi_ivp_prob<T1>* ivp;
  // Diagonal coefficient
  T1 gamma;
  // Method coefficients, transformed form (a, alpha, m, m_err: see ivp): c strictly lower
  // triangular, row-major, sums gamma_i of the rows of Gamma
  const T1 *c, *gamma_sum;
  // Iteration matrix M/(h gamma) - J: sparsity, positions of ODE_X, ALG_X, ODE_Z, ALG_Z, diag
  const casadi_int *sp_iter, *iter_map;
  // Sparse QR of the iteration matrix and of d(alg)/dz
  const casadi_int *sp_iter_v, *sp_iter_r, *iter_prinv, *iter_pc;
  const casadi_int *sp_alg_v, *sp_alg_r, *alg_prinv, *alg_pc;
  // Consistent initialization: tolerance on scaled Newton step, iteration limit
  T1 init_tol;
  casadi_int init_max_iter;
  // Derived: nonzeros of the iteration matrix and of the QR factors
  casadi_int nnz_iter, nnz_iter_v, nnz_iter_r, nnz_alg_v, nnz_alg_r;
  // Derived: length of the work vectors
  casadi_int nv, n_qr_w;
};
// C-REPLACE "casadi_rosenbrock_prob<T1>" "struct casadi_rosenbrock_prob"

// SYMBOL "rosenbrock_setup"
template<typename T1>
void casadi_rosenbrock_setup(casadi_rosenbrock_prob<T1>* rp) {
  const casadi_ivp_prob<T1>* p = rp->ivp;
  casadi_int n1 = p->nx1 + p->nz1;
  rp->nnz_iter = casadi_ivp_nnz(rp->sp_iter);
  rp->nnz_iter_v = casadi_ivp_nnz(rp->sp_iter_v);
  rp->nnz_iter_r = casadi_ivp_nnz(rp->sp_iter_r);
  rp->nnz_alg_v = p->nz1 > 0 ? casadi_ivp_nnz(rp->sp_alg_v) : 0;
  rp->nnz_alg_r = p->nz1 > 0 ? casadi_ivp_nnz(rp->sp_alg_r) : 0;
  // Work vectors: forward and adjoint vectors with all directions
  rp->nv = p->ny;
  if (p->nr + p->nrp + p->nru > rp->nv) rp->nv = p->nr + p->nrp + p->nru;
  if (n1 * p->nadj * (1 + p->nfwd) > rp->nv) rp->nv = n1 * p->nadj * (1 + p->nfwd);
  // Work vector of the sparse QR: extended number of rows
  rp->n_qr_w = rp->sp_iter_v[0];
  if (p->nz1 > 0 && rp->sp_alg_v[0] > rp->n_qr_w) rp->n_qr_w = rp->sp_alg_v[0];
}

// SYMBOL "rosenbrock_next_t"
typedef enum {
  // Phases: each one consumes the result of the previously requested task
  ROSENBROCK_INIT_BEGIN = IVP_NUM_ENTRY,
  ROSENBROCK_INIT_LINEARIZE,
  ROSENBROCK_INIT_NEWTON_STEP,
  ROSENBROCK_INIT_SENS_RHS,
  ROSENBROCK_INIT_SENS,
  ROSENBROCK_INIT_END,
  ROSENBROCK_STEP_BEGIN,
  ROSENBROCK_STEP_SIZE,
  ROSENBROCK_STEP_INITIAL_SIZE,
  ROSENBROCK_STEP_LINEARIZE,
  ROSENBROCK_STEP_FACTORIZE,
  ROSENBROCK_STEP_REPLAY,
  ROSENBROCK_STAGE_SOLVE,
  ROSENBROCK_STEP_CHECK,
  ROSENBROCK_STEP_EVENT,
  ROSENBROCK_OUTPUT_CORRECT,
  ROSENBROCK_BACK_BEGIN,
  ROSENBROCK_BACK_INIT_SOLVE,
  ROSENBROCK_BACK_INIT_END,
  ROSENBROCK_BACK_ADJ_BEGIN,
  ROSENBROCK_BACK_ADJ_STAGE,
  ROSENBROCK_BACK_ADJ_SPREAD,
  ROSENBROCK_BACK_STEP_END,
  ROSENBROCK_IMPULSE_LINEARIZE,
  ROSENBROCK_IMPULSE_CORRECT,
  ROSENBROCK_IMPULSE_END
} casadi_rosenbrock_next_t;

// C-REPLACE "casadi_rosenbrock_data<T1>" "struct casadi_rosenbrock_data"
// SYMBOL "rosenbrock_data"
template<typename T1>
struct casadi_rosenbrock_data {
  // Problem structure
  const casadi_rosenbrock_prob<T1>* prob;
  // Method data of the nominal integration whose steps are replayed (or null)
  casadi_rosenbrock_data<T1>* nom;
  // h * gamma, time derivative of F
  T1 h_gamma;
  T1* dfdt;
  // Iteration matrix and its QR factors, QR factors of d(alg)/dz
  T1 *iter, *iter_v, *iter_r, *iter_beta, *alg_v, *alg_r, *alg_beta, *qr_w;
  // Work vectors of length nv
  T1 *v1, *v2;
  // Consistent initialization pending, Newton iteration
  int init_pending;
  casadi_int init_iter;
};

// SYMBOL "rosenbrock_work"
template<typename T1>
void casadi_rosenbrock_work(const casadi_rosenbrock_prob<T1>* rp, casadi_int* sz_arg,
    casadi_int* sz_res, casadi_int* sz_iw, casadi_int* sz_w) {
  const casadi_ivp_prob<T1>* p = rp->ivp;
  casadi_int n1 = p->nx1 + p->nz1;
  // dF/dt
  *sz_w += p->ny;
  // Iteration matrix, QR factors of it and of d(alg)/dz, QR work vector
  *sz_w += rp->nnz_iter + rp->nnz_iter_v + rp->nnz_iter_r + n1;
  *sz_w += rp->nnz_alg_v + rp->nnz_alg_r + p->nz1 + rp->n_qr_w;
  // Work vectors
  *sz_w += 2 * rp->nv;
}

// SYMBOL "rosenbrock_set_work"
template<typename T1>
void casadi_rosenbrock_set_work(casadi_rosenbrock_data<T1>* r, const T1*** arg, T1*** res,
    casadi_int** iw, T1** w) {
  const casadi_rosenbrock_prob<T1>* rp = r->prob;
  const casadi_ivp_prob<T1>* p = rp->ivp;
  casadi_int n1 = p->nx1 + p->nz1;
  r->dfdt = *w; *w += p->ny;
  r->iter = *w; *w += rp->nnz_iter;
  r->iter_v = *w; *w += rp->nnz_iter_v;
  r->iter_r = *w; *w += rp->nnz_iter_r;
  r->iter_beta = *w; *w += n1;
  r->alg_v = *w; *w += rp->nnz_alg_v;
  r->alg_r = *w; *w += rp->nnz_alg_r;
  r->alg_beta = *w; *w += p->nz1;
  r->qr_w = *w; *w += rp->n_qr_w;
  r->v1 = *w; *w += rp->nv;
  r->v2 = *w; *w += rp->nv;
  r->nom = 0;
}

// ---------------------------------------------------------------------------------------
// Linear algebra of the sensitivity-augmented iteration matrix
// ---------------------------------------------------------------------------------------

// SYMBOL "rosenbrock_gather"
template<typename T1>
void casadi_rosenbrock_gather(const T1* x, const T1* z, casadi_int nx1, casadi_int nz1,
    T1* v) {
  // v = [x; z]
  casadi_copy(x, nx1, v);
  casadi_copy(z, nz1, v + nx1);
}

// SYMBOL "rosenbrock_scatter"
template<typename T1>
void casadi_rosenbrock_scatter(const T1* v, casadi_int nx1, casadi_int nz1, T1* x, T1* z) {
  // [x; z] = v
  casadi_copy(v, nx1, x);
  casadi_copy(v + nx1, nz1, z);
}

// SYMBOL "rosenbrock_dfdt"
template<typename T1>
void casadi_rosenbrock_dfdt(casadi_ivp_data<T1>* d, casadi_rosenbrock_data<T1>* r) {
  // Time derivative of [ode; alg; quad], with sensitivities, from the time columns
  const casadi_ivp_prob<T1>* p = d->prob;
  casadi_int k, dir, off, n1;
  T1 one = 1;
  casadi_clear(r->dfdt, p->ny);
  for (k = 0; k < 3; ++k) {
    // Empty if the right-hand side does not depend on time
    if (p->sp_jac[IVP_ODE_T + k][1] == 0) continue;
    off = k == 0 ? 0 : k == 1 ? p->nx : p->nx + p->nz;
    n1 = k == 0 ? p->nx1 : k == 1 ? p->nz1 : p->nq1;
    for (dir = -1; dir < p->nfwd; ++dir) {
      casadi_ivp_jac_mv(d, IVP_ODE_T + k, dir, &one, r->dfdt + off + (dir + 1) * n1, 0);
    }
  }
}

// SYMBOL "rosenbrock_factorize"
template<typename T1>
int casadi_rosenbrock_factorize(casadi_ivp_data<T1>* d, casadi_rosenbrock_data<T1>* r) {
  // QR of M/(h gamma) - J for the nominal [x; z]; returns 1 if singular
  const casadi_ivp_prob<T1>* p = d->prob;
  const casadi_rosenbrock_prob<T1>* rp = r->prob;
  casadi_int k, i, n, pos = 0;
  const casadi_int blk[4] = {IVP_ODE_X, IVP_ALG_X, IVP_ODE_Z, IVP_ALG_Z};
  const T1* J;
  T1 rmin = 0;
  casadi_clear(r->iter, rp->nnz_iter);
  for (k = 0; k < 4; ++k) {
    J = casadi_ivp_jac_nz(d, d->jac, blk[k], 1, 0);
    n = casadi_ivp_nnz(p->sp_jac[blk[k]]);
    for (i = 0; i < n; ++i) r->iter[rp->iter_map[pos++]] -= J[i];
  }
  for (i = 0; i < p->nx1; ++i) r->iter[rp->iter_map[pos++]] += 1 / r->h_gamma;
  casadi_qr(rp->sp_iter, r->iter, r->qr_w, rp->sp_iter_v, r->iter_v, rp->sp_iter_r, r->iter_r,
    r->iter_beta, rp->iter_prinv, rp->iter_pc);
  casadi_qr_singular(&rmin, static_cast<casadi_int*>(0), r->iter_r, rp->sp_iter_r,
    rp->iter_pc, 0.);
  return !(rmin > 0) || !casadi_ivp_isfinite(rmin);
}

// SYMBOL "rosenbrock_iter_solve"
template<typename T1>
void casadi_rosenbrock_iter_solve(casadi_rosenbrock_data<T1>* r, T1* v, casadi_int nrhs,
    int tr) {
  // Solve with the factorized iteration matrix (or its transpose), nrhs right-hand sides
  const casadi_rosenbrock_prob<T1>* rp = r->prob;
  casadi_qr_solve(v, nrhs, tr, rp->sp_iter_v, r->iter_v, rp->sp_iter_r, r->iter_r, r->iter_beta,
    rp->iter_prinv, rp->iter_pc, r->qr_w);
}

// SYMBOL "rosenbrock_solve"
template<typename T1>
void casadi_rosenbrock_solve(casadi_ivp_data<T1>* d, casadi_rosenbrock_data<T1>* r, T1* v) {
  // A U_0 = r_0, A U_d = r_d + dJ_d U_0 (A = M/hg - J), quadratures U_q = hg (r_q + J_qy U)
  const casadi_ivp_prob<T1>* p = d->prob;
  casadi_int dir, n1 = p->nx1 + p->nz1;
  T1 *u0 = r->v1, *ud, *vx = v, *vz = v + p->nx, *vq = v + p->nx + p->nz, *vqd;
  casadi_rosenbrock_gather(vx, vz, p->nx1, p->nz1, u0);
  casadi_rosenbrock_iter_solve(r, u0, 1, 0);
  for (dir = 0; dir < p->nfwd; ++dir) {
    ud = r->v2 + dir * n1;
    casadi_rosenbrock_gather(vx + (dir + 1) * p->nx1, vz + (dir + 1) * p->nz1,
      p->nx1, p->nz1, ud);
    casadi_ivp_jac_mv(d, IVP_ODE_X, dir, u0, ud, 0);
    casadi_ivp_jac_mv(d, IVP_ODE_Z, dir, u0 + p->nx1, ud, 0);
    casadi_ivp_jac_mv(d, IVP_ALG_X, dir, u0, ud + p->nx1, 0);
    casadi_ivp_jac_mv(d, IVP_ALG_Z, dir, u0 + p->nx1, ud + p->nx1, 0);
  }
  if (p->nfwd > 0) casadi_rosenbrock_iter_solve(r, r->v2, p->nfwd, 0);
  for (dir = 0; dir < p->nfwd; ++dir) {
    ud = r->v2 + dir * n1;
    vqd = vq + (dir + 1) * p->nq1;
    casadi_ivp_jac_mv(d, IVP_QUAD_X, -1, ud, vqd, 0);
    casadi_ivp_jac_mv(d, IVP_QUAD_Z, -1, ud + p->nx1, vqd, 0);
    casadi_ivp_jac_mv(d, IVP_QUAD_X, dir, u0, vqd, 0);
    casadi_ivp_jac_mv(d, IVP_QUAD_Z, dir, u0 + p->nx1, vqd, 0);
    casadi_scal(p->nq1, r->h_gamma, vqd);
    casadi_rosenbrock_scatter(ud, p->nx1, p->nz1, vx + (dir + 1) * p->nx1,
      vz + (dir + 1) * p->nz1);
  }
  casadi_ivp_jac_mv(d, IVP_QUAD_X, -1, u0, vq, 0);
  casadi_ivp_jac_mv(d, IVP_QUAD_Z, -1, u0 + p->nx1, vq, 0);
  casadi_scal(p->nq1, r->h_gamma, vq);
  casadi_rosenbrock_scatter(u0, p->nx1, p->nz1, vx, vz);
}

// SYMBOL "rosenbrock_solve_trans"
template<typename T1>
void casadi_rosenbrock_solve_trans(casadi_ivp_data<T1>* d, casadi_rosenbrock_data<T1>* r,
    T1* v) {
  // A^T mu = l + J_qy^T hg l_q, then A^T dmu_d = dl_d + J_qy^T hg dl_q + dJ_d^T [mu; hg l_q]
  const casadi_ivp_prob<T1>* p = d->prob;
  casadi_int a, dir, n1 = p->nx1 + p->nz1;
  casadi_int nrx1 = p->nx1 * p->nadj, nrz1 = p->nz1 * p->nadj, nrq1 = p->nq1 * p->nadj;
  T1 *vx = v, *vz = v + p->nrx, *vq = v + p->nrx + p->nrz;
  T1 *mu, *mq, *dmu, *dmq, *dvx, *dvz;
  casadi_scal(p->nrq, r->h_gamma, vq);
  for (a = 0; a < p->nadj; ++a) {
    mu = r->v1 + a * n1;
    mq = vq + a * p->nq1;
    casadi_rosenbrock_gather(vx + a * p->nx1, vz + a * p->nz1, p->nx1, p->nz1, mu);
    casadi_ivp_jac_mv(d, IVP_QUAD_X, -1, mq, mu, 1);
    casadi_ivp_jac_mv(d, IVP_QUAD_Z, -1, mq, mu + p->nx1, 1);
  }
  casadi_rosenbrock_iter_solve(r, r->v1, p->nadj, 1);
  for (dir = 0; dir < p->nfwd; ++dir) {
    for (a = 0; a < p->nadj; ++a) {
      mu = r->v1 + a * n1;
      mq = vq + a * p->nq1;
      dmu = r->v2 + (dir * p->nadj + a) * n1;
      dmq = vq + (dir + 1) * nrq1 + a * p->nq1;
      dvx = vx + (dir + 1) * nrx1 + a * p->nx1;
      dvz = vz + (dir + 1) * nrz1 + a * p->nz1;
      casadi_rosenbrock_gather(dvx, dvz, p->nx1, p->nz1, dmu);
      casadi_ivp_jac_mv(d, IVP_QUAD_X, -1, dmq, dmu, 1);
      casadi_ivp_jac_mv(d, IVP_QUAD_Z, -1, dmq, dmu + p->nx1, 1);
      casadi_ivp_jac_mv(d, IVP_QUAD_X, dir, mq, dmu, 1);
      casadi_ivp_jac_mv(d, IVP_QUAD_Z, dir, mq, dmu + p->nx1, 1);
      casadi_ivp_jac_mv(d, IVP_ODE_X, dir, mu, dmu, 1);
      casadi_ivp_jac_mv(d, IVP_ALG_X, dir, mu + p->nx1, dmu, 1);
      casadi_ivp_jac_mv(d, IVP_ODE_Z, dir, mu, dmu + p->nx1, 1);
      casadi_ivp_jac_mv(d, IVP_ALG_Z, dir, mu + p->nx1, dmu + p->nx1, 1);
    }
  }
  if (p->nfwd > 0) casadi_rosenbrock_iter_solve(r, r->v2, p->nadj * p->nfwd, 1);
  for (a = 0; a < p->nadj; ++a) {
    casadi_rosenbrock_scatter(r->v1 + a * n1, p->nx1, p->nz1, vx + a * p->nx1,
      vz + a * p->nz1);
    for (dir = 0; dir < p->nfwd; ++dir) {
      casadi_rosenbrock_scatter(r->v2 + (dir * p->nadj + a) * n1, p->nx1, p->nz1,
        vx + (dir + 1) * nrx1 + a * p->nx1, vz + (dir + 1) * nrz1 + a * p->nz1);
    }
  }
}

// SYMBOL "rosenbrock_factorize_alg"
template<typename T1>
int casadi_rosenbrock_factorize_alg(casadi_ivp_data<T1>* d, casadi_rosenbrock_data<T1>* r) {
  // QR of d(alg)/dz; returns 1 if singular
  const casadi_rosenbrock_prob<T1>* rp = r->prob;
  T1 rmin = 0;
  casadi_qr(d->prob->sp_jac[IVP_ALG_Z], casadi_ivp_jac_nz(d, d->jac, IVP_ALG_Z, 1, 0),
    r->qr_w, rp->sp_alg_v, r->alg_v, rp->sp_alg_r, r->alg_r, r->alg_beta, rp->alg_prinv,
    rp->alg_pc);
  casadi_qr_singular(&rmin, static_cast<casadi_int*>(0), r->alg_r, rp->sp_alg_r, rp->alg_pc,
    0.);
  return !(rmin > 0) || !casadi_ivp_isfinite(rmin);
}

// SYMBOL "rosenbrock_alg_solve"
template<typename T1>
void casadi_rosenbrock_alg_solve(casadi_rosenbrock_data<T1>* r, T1* v, casadi_int nrhs,
    int tr) {
  // Solve with the factorized d(alg)/dz (or its transpose), nrhs right-hand sides
  const casadi_rosenbrock_prob<T1>* rp = r->prob;
  casadi_qr_solve(v, nrhs, tr, rp->sp_alg_v, r->alg_v, rp->sp_alg_r, r->alg_r, r->alg_beta,
    rp->alg_prinv, rp->alg_pc, r->qr_w);
}

// SYMBOL "rosenbrock_newton_step"
template<typename T1>
void casadi_rosenbrock_newton_step(casadi_ivp_data<T1>* d, casadi_rosenbrock_data<T1>* r,
    T1* v) {
  // Augmented Newton step: g_z w_0 = v_0, g_z w_d = v_d - dg_z,d w_0
  const casadi_ivp_prob<T1>* p = d->prob;
  casadi_int dir;
  T1* vd;
  casadi_rosenbrock_alg_solve(r, v, 1, 0);
  for (dir = 0; dir < p->nfwd; ++dir) {
    vd = v + (dir + 1) * p->nz1;
    casadi_scal(p->nz1, -1., vd);
    casadi_ivp_jac_mv(d, IVP_ALG_Z, dir, v, vd, 0);
    casadi_scal(p->nz1, -1., vd);
  }
  if (p->nfwd > 0) casadi_rosenbrock_alg_solve(r, v + p->nz1, p->nfwd, 0);
}

// SYMBOL "rosenbrock_newton_step_trans"
template<typename T1>
void casadi_rosenbrock_newton_step_trans(casadi_ivp_data<T1>* d,
    casadi_rosenbrock_data<T1>* r, T1* v) {
  // Transposed: g_z^T mu = v, g_z^T dmu_d = dv_d - dg_z,d^T mu, for each adjoint direction
  const casadi_ivp_prob<T1>* p = d->prob;
  casadi_int a, dir, nrz1 = p->nz1 * p->nadj;
  T1 *mu, *dmu;
  casadi_rosenbrock_alg_solve(r, v, p->nadj, 1);
  for (dir = 0; dir < p->nfwd; ++dir) {
    for (a = 0; a < p->nadj; ++a) {
      mu = v + a * p->nz1;
      dmu = v + (dir + 1) * nrz1 + a * p->nz1;
      casadi_scal(p->nz1, -1., dmu);
      casadi_ivp_jac_mv(d, IVP_ALG_Z, dir, mu, dmu, 1);
      casadi_scal(p->nz1, -1., dmu);
    }
  }
  if (p->nfwd > 0) casadi_rosenbrock_alg_solve(r, v + nrz1, p->nadj * p->nfwd, 1);
}

// SYMBOL "rosenbrock_accum_jac_adj"
template<typename T1>
void casadi_rosenbrock_accum_jac_adj(casadi_ivp_data<T1>* d, const T1* w, const T1* u,
    T1 coef) {
  // jac_adj_a += w_a u^T (x, z columns), coef w_a (time column), and forward derivatives
  const casadi_ivp_prob<T1>* p = d->prob;
  casadi_int a, row, c, k, dir, n_row, n_col;
  const T1 *wa, *dwa, *uc;
  T1 one = 1;
  for (a = 0; a < p->nadj; ++a) {
    for (row = 0; row < 3; ++row) {
      // Rows: ode, alg, quad
      n_row = row == 0 ? p->nx1 : row == 1 ? p->nz1 : p->nq1;
      wa = (row == 0 ? w : row == 1 ? w + p->nrx : w + p->nrx + p->nrz) + a * n_row;
      for (c = 0; c < 3; ++c) {
        // Columns: x, z, t
        k = 3 * c + row;
        if (casadi_ivp_nnz(p->sp_jac[k]) == 0) continue;
        n_col = c == 0 ? p->nx1 : p->nz1;
        uc = c == 0 ? u : u + p->nx;
        casadi_rank1(casadi_ivp_jac_nz(d, d->jac_adj, k, p->nadj, a), p->sp_jac[k],
          c < 2 ? one : coef, wa, c < 2 ? uc : &one);
        for (dir = 0; dir < p->nfwd; ++dir) {
          dwa = wa + (dir + 1) * n_row * p->nadj;
          casadi_rank1(casadi_ivp_jac_nz(d, d->jac_adj_fwd, k, p->nadj * p->nfwd,
            dir * p->nadj + a), p->sp_jac[k], c < 2 ? one : coef, dwa, c < 2 ? uc : &one);
          if (c < 2) {
            casadi_rank1(casadi_ivp_jac_nz(d, d->jac_adj_fwd, k, p->nadj * p->nfwd,
              dir * p->nadj + a), p->sp_jac[k], one, wa, uc + (dir + 1) * n_col);
          }
        }
      }
    }
  }
}

// ---------------------------------------------------------------------------------------
// Consistent initialization
// ---------------------------------------------------------------------------------------

// SYMBOL "rosenbrock_newton_norm"
template<typename T1>
T1 casadi_rosenbrock_newton_norm(const casadi_ivp_data<T1>* d, const T1* dz) {
  // Weighted RMS norm of a Newton step for the error controlled algebraic states
  const casadi_ivp_prob<T1>* p = d->prob;
  casadi_int i;
  T1 ri, sum = 0;
  if (p->nz_err == 0) return 0;
  for (i = 0; i < p->nz_err; ++i) {
    ri = dz[i] / (p->abstol + p->reltol * fabs(d->y[p->nx + i]));
    sum += ri * ri;
  }
  return sqrt(sum / static_cast<T1>(p->nz_err));
}

// SYMBOL "rosenbrock_stage_rhs"
template<typename T1>
void casadi_rosenbrock_stage_rhs(casadi_ivp_data<T1>* d, casadi_rosenbrock_data<T1>* r,
    casadi_int i) {
  // U_i += sum_{j<i} (c_ij/h) M U_j + gamma_i h dF/dt, with U_i holding F on entry
  const casadi_ivp_prob<T1>* p = d->prob;
  const casadi_rosenbrock_prob<T1>* rp = r->prob;
  casadi_int j;
  T1 *Ui = d->U + i * p->ny, cij;
  for (j = 0; j < i; ++j) {
    cij = rp->c[i * p->n_stages + j] / d->h;
    // M is the identity on x and q, zero on z
    casadi_axpy(p->nx, cij, d->U + j * p->ny, Ui);
    casadi_axpy(p->nq, cij, d->U + j * p->ny + p->nx + p->nz, Ui + p->nx + p->nz);
  }
  casadi_axpy(p->ny, rp->gamma_sum[i] * d->h, r->dfdt, Ui);
}

// SYMBOL "rosenbrock_stage_request"
template<typename T1>
int casadi_rosenbrock_stage_request(casadi_ivp_data<T1>* d) {
  // Right-hand side at the argument of the current stage
  const casadi_ivp_prob<T1>* p = d->prob;
  casadi_ivp_stage_argument(d, d->stage);
  d->t_arg = d->t0 + p->alpha[d->stage] * d->h;
  d->y_arg = d->y_stage;
  d->res = d->U + d->stage * p->ny;
  return casadi_ivp_request(d, IVP_EVAL_RHS, ROSENBROCK_STAGE_SOLVE);
}

// ---------------------------------------------------------------------------------------
// Forward integration and discrete adjoint
// ---------------------------------------------------------------------------------------

// SYMBOL "rosenbrock"
template<typename T1>
int casadi_rosenbrock(casadi_ivp_data<T1>* d, casadi_rosenbrock_data<T1>* r);

// SYMBOL "rosenbrock_forward"
template<typename T1>
int casadi_rosenbrock_forward(casadi_ivp_data<T1>* d, casadi_rosenbrock_data<T1>* r) {
  // Forward integration; returns 1 if a task is requested, 0 when done or failed
  const casadi_ivp_prob<T1>* p = d->prob;
  const casadi_rosenbrock_prob<T1>* rp = r->prob;
  T1 err, *v;
  int ret;
  switch (d->next) {
  case IVP_ADVANCE:
    if (d->reset) {
      // New state or controls: consistent algebraic states, a new run of fixed steps
      casadi_ivp_begin_run(d);
      r->init_pending = p->nz > 0;
    }
    d->next = r->init_pending ? ROSENBROCK_INIT_BEGIN : ROSENBROCK_STEP_BEGIN;
    return casadi_rosenbrock_forward(d, r);
  // Consistent initialization: Newton's method for alg(t, x, z) = 0, x fixed
  case ROSENBROCK_INIT_BEGIN:
    if (d->fail) return casadi_ivp_stop(d, IVP_TAPE_FAILED);
    if ((ret = casadi_ivp_reserve(d, IVP_RECORD_INIT, ROSENBROCK_INIT_BEGIN)) >= 0) return ret;
    casadi_copy(d->u, p->nu, d->u_step);
    r->init_iter = 0;
    d->t_arg = d->t;
    d->y_arg = d->y;
    d->res = d->U;
    return casadi_ivp_request(d, IVP_EVAL_RHS, ROSENBROCK_INIT_LINEARIZE);
  case ROSENBROCK_INIT_LINEARIZE:
    if (d->fail) return casadi_ivp_stop(d, IVP_INIT_FAILED);
    d->n_rhs++;
    return casadi_ivp_request(d, IVP_EVAL_JAC, ROSENBROCK_INIT_NEWTON_STEP);
  case ROSENBROCK_INIT_NEWTON_STEP:
    if (d->fail || casadi_rosenbrock_factorize_alg(d, r)) {
      return casadi_ivp_stop(d, IVP_INIT_FAILED);
    }
    d->n_jac++;
    // z <- z - g_z^{-1} g, with sensitivities
    v = d->U + p->nx;
    casadi_rosenbrock_newton_step(d, r, v);
    casadi_axpy(p->nz, -1., v, d->y + p->nx);
    err = casadi_rosenbrock_newton_norm(d, v);
    if (++r->init_iter > rp->init_max_iter || !casadi_ivp_isfinite(err)) {
      return casadi_ivp_stop(d, IVP_INIT_FAILED);
    }
    if (err > rp->init_tol) {
      // Next Newton iteration
      d->res = d->U;
      return casadi_ivp_request(d, IVP_EVAL_RHS, ROSENBROCK_INIT_LINEARIZE);
    }
    if (p->nfwd == 0) {
      d->next = ROSENBROCK_INIT_END;
      return casadi_rosenbrock_forward(d, r);
    }
    // Consistent algebraic sensitivities at the solution (implicit function theorem)
    return casadi_ivp_request(d, IVP_EVAL_JAC, ROSENBROCK_INIT_SENS_RHS);
  case ROSENBROCK_INIT_SENS_RHS:
    if (d->fail || casadi_rosenbrock_factorize_alg(d, r)) {
      return casadi_ivp_stop(d, IVP_INIT_FAILED);
    }
    d->n_jac++;
    d->res = d->U;
    return casadi_ivp_request(d, IVP_EVAL_RHS, ROSENBROCK_INIT_SENS);
  case ROSENBROCK_INIT_SENS:
    if (d->fail) return casadi_ivp_stop(d, IVP_INIT_FAILED);
    d->n_rhs++;
    // dz_d -= g_z^{-1} (g_x dx_d + g_z dz_d + g_p dp_d + g_u du_d)
    v = d->U + p->nx + p->nz1;
    casadi_rosenbrock_alg_solve(r, v, p->nfwd, 0);
    casadi_axpy(p->nz1 * p->nfwd, -1., v, d->y + p->nx + p->nz1);
    d->next = ROSENBROCK_INIT_END;
    return casadi_rosenbrock_forward(d, r);
  case ROSENBROCK_INIT_END:
    r->init_pending = 0;
    if (d->taping) {
      d->record = IVP_RECORD_INIT;
      casadi_ivp_tape_push(d, d->y, static_cast<T1>(0), d->t, d->t);
    }
    d->next = ROSENBROCK_STEP_BEGIN;
    return casadi_rosenbrock_forward(d, r);
  // Integrator step
  case ROSENBROCK_STEP_BEGIN:
    if (d->fail) return casadi_ivp_stop(d, IVP_TAPE_FAILED);
    if (d->t >= d->t_out) {
      if (d->pause) return casadi_ivp_stop(d, IVP_SUCCESS);
      if (!casadi_ivp_output(d) || p->nz == 0) return casadi_ivp_stop(d, IVP_SUCCESS);
      // Output inside the last step, algebraic states: Newton step on alg(t_out, x, z) = 0,
      // d(alg)/dz at the step start
      if (casadi_rosenbrock_factorize_alg(d, r)) {
        return casadi_ivp_stop(d, IVP_EVAL_FAILED);
      }
      d->t_arg = d->t_out;
      d->y_arg = d->y_out;
      d->res = d->y_stage;
      return casadi_ivp_request(d, IVP_EVAL_RHS, ROSENBROCK_OUTPUT_CORRECT);
    }
    if (d->jump) {
      // State changed by an event: restart
      d->reset = 1;
      d->next = IVP_ADVANCE;
      return casadi_rosenbrock_forward(d, r);
    }
    if (casadi_ivp_step_limit(d)) return casadi_ivp_stop(d, IVP_MAX_NUM_STEPS);
    if ((ret = casadi_ivp_reserve(d, IVP_RECORD_STEP, ROSENBROCK_STEP_BEGIN)) >= 0) return ret;
    // Controls of the new step
    casadi_copy(d->u, p->nu, d->u_step);
    if (d->nom) {
      // Step of the nominal integration
      d->next = ROSENBROCK_STEP_REPLAY;
      return casadi_rosenbrock_forward(d, r);
    }
    if (p->step_size > 0) {
      d->h = d->run_h;
      d->next = ROSENBROCK_STEP_LINEARIZE;
      return casadi_rosenbrock_forward(d, r);
    }
    return casadi_ivp_request_control(d, IVP_CONTROL_SIZE, ROSENBROCK_STEP_SIZE);
  case ROSENBROCK_STEP_SIZE:
    if (d->h > 0) {
      d->next = ROSENBROCK_STEP_LINEARIZE;
      return casadi_rosenbrock_forward(d, r);
    }
    // Estimate from the right-hand side at the initial point
    d->t_arg = d->t;
    d->y_arg = d->y;
    d->res = d->U;
    return casadi_ivp_request(d, IVP_EVAL_RHS, ROSENBROCK_STEP_INITIAL_SIZE);
  case ROSENBROCK_STEP_INITIAL_SIZE:
    if (d->fail) return casadi_ivp_stop(d, IVP_EVAL_FAILED);
    d->n_rhs++;
    return casadi_ivp_request_control(d, IVP_CONTROL_INIT, ROSENBROCK_STEP_SIZE);
  case ROSENBROCK_STEP_LINEARIZE:
    // Linearize at the beginning of the step
    d->t0 = d->t;
    casadi_copy(d->y, p->ny, d->y0);
    d->t_arg = d->t0;
    d->y_arg = d->y0;
    return casadi_ivp_request(d, IVP_EVAL_JAC, ROSENBROCK_STEP_FACTORIZE);
  case ROSENBROCK_STEP_REPLAY:
    // Advance the nominal integration by one step; its tasks are handed to the caller
    if (casadi_rosenbrock(d->nom, r->nom)) {
      d->active = d->nom;
      return 1;
    }
    d->active = d;
    if (d->nom->status != IVP_SUCCESS) return casadi_ivp_stop(d, d->nom->status);
    if (d->nom->next == IVP_DONE) return casadi_ivp_stop(d, IVP_EVAL_FAILED);
    // Same step size, same end time
    d->h = d->nom->h;
    d->next = ROSENBROCK_STEP_LINEARIZE;
    return casadi_rosenbrock_forward(d, r);
  case ROSENBROCK_STEP_FACTORIZE:
    if (d->fail) return casadi_ivp_stop(d, IVP_EVAL_FAILED);
    d->n_jac++;
    casadi_rosenbrock_dfdt(d, r);
    r->h_gamma = d->h * rp->gamma;
    if (casadi_rosenbrock_factorize(d, r)) {
      // Singular iteration matrix: retry with a smaller step, if adaptive
      if (casadi_ivp_prescribed(d)) return casadi_ivp_stop(d, IVP_EVAL_FAILED);
      return casadi_ivp_request_control(d, IVP_CONTROL_FAIL, ROSENBROCK_STEP_CHECK);
    }
    d->stage = 0;
    return casadi_rosenbrock_stage_request(d);
  case ROSENBROCK_STAGE_SOLVE:
    if (d->fail) {
      // Evaluation failed: retry with a smaller step, if adaptive
      if (casadi_ivp_prescribed(d)) return casadi_ivp_stop(d, IVP_EVAL_FAILED);
      return casadi_ivp_request_control(d, IVP_CONTROL_FAIL, ROSENBROCK_STEP_CHECK);
    }
    d->n_rhs++;
    casadi_rosenbrock_stage_rhs(d, r, d->stage);
    casadi_rosenbrock_solve(d, r, d->U + d->stage * p->ny);
    if (++d->stage < p->n_stages) return casadi_rosenbrock_stage_request(d);
    // Candidate solution y = y0 + sum_i m_i U_i
    casadi_copy(d->y0, p->ny, d->y);
    casadi_mv_dense(d->U, p->ny, p->n_stages, p->m, d->y, 0);
    if (!casadi_ivp_prescribed(d)) {
      return casadi_ivp_request_control(d, IVP_CONTROL_TEST, ROSENBROCK_STEP_CHECK);
    }
    d->accepted = 1;
    d->next = ROSENBROCK_STEP_CHECK;
    return casadi_rosenbrock_forward(d, r);
  case ROSENBROCK_STEP_CHECK:
    if (d->status != IVP_SUCCESS) return casadi_ivp_stop(d, d->status);
    if (!d->accepted) {
      // Retry from t0
      d->next = ROSENBROCK_STEP_BEGIN;
      return casadi_rosenbrock_forward(d, r);
    }
    casadi_ivp_accept(d);
    if (p->ne > 0) return casadi_ivp_request(d, IVP_EVENT, ROSENBROCK_STEP_EVENT);
    d->next = ROSENBROCK_STEP_EVENT;
    return casadi_rosenbrock_forward(d, r);
  case ROSENBROCK_STEP_EVENT:
    if (d->status != IVP_SUCCESS) return casadi_ivp_stop(d, d->status);
    d->next = ROSENBROCK_STEP_BEGIN;
    if (d->pause) {
      // Hand the step to the integration that replays it
      d->status = IVP_SUCCESS;
      return 0;
    }
    return casadi_rosenbrock_forward(d, r);
  case ROSENBROCK_OUTPUT_CORRECT:
    if (d->fail) return casadi_ivp_stop(d, IVP_EVAL_FAILED);
    d->n_rhs++;
    v = d->y_stage + p->nx;
    casadi_rosenbrock_newton_step(d, r, v);
    casadi_axpy(p->nz, -1., v, d->y_out + p->nx);
    return casadi_ivp_stop(d, IVP_SUCCESS);
  default:
    return 0;
  }
}

// SYMBOL "rosenbrock_backward"
template<typename T1>
int casadi_rosenbrock_backward(casadi_ivp_data<T1>* d, casadi_rosenbrock_data<T1>* r) {
  // Discrete adjoint of the taped steps; returns 1 if a task is requested, 0 when done
  const casadi_ivp_prob<T1>* p = d->prob;
  const casadi_rosenbrock_prob<T1>* rp = r->prob;
  casadi_int i, j, s = p->n_stages;
  T1 cij, *rUi, *rUj;
  switch (d->next) {
  case IVP_RETREAT:
  case ROSENBROCK_BACK_BEGIN:
    if (!casadi_ivp_back_begin(d)) return 0;
    if (d->record == IVP_RECORD_EVENT) {
      return casadi_ivp_request(d, IVP_EVENT_ADJ, ROSENBROCK_BACK_BEGIN);
    }
    // Linearize at the beginning of the record
    d->t_arg = d->t0;
    d->y_arg = d->y0;
    return casadi_ivp_request(d, IVP_EVAL_JAC,
      d->record == IVP_RECORD_INIT ? ROSENBROCK_BACK_INIT_SOLVE : ROSENBROCK_BACK_ADJ_BEGIN);
  case ROSENBROCK_BACK_INIT_SOLVE:
    if (d->fail || casadi_rosenbrock_factorize_alg(d, r)) break;
    // Implicit function theorem: mu = (d(alg)/dz)^{-T} zbar, xbar -= (d(alg)/dx)^T mu
    casadi_clear(r->v1, p->nr);
    casadi_copy(d->ry + p->nrx, p->nrz, r->v1 + p->nrx);
    casadi_rosenbrock_newton_step_trans(d, r, r->v1 + p->nrx);
    casadi_scal(p->nrz, -1., r->v1 + p->nrx);
    casadi_ivp_clear_vjp(d);
    d->w = r->v1;
    return casadi_ivp_request(d, IVP_EVAL_VJP, ROSENBROCK_BACK_INIT_END);
  case ROSENBROCK_BACK_INIT_END:
    if (d->fail) break;
    d->n_vjp++;
    casadi_ivp_add_vjp(d, d->ry);
    // The algebraic adjoint has been propagated into the differential states
    casadi_clear(d->ry + p->nrx, p->nrz);
    d->t = d->t0;
    d->next = ROSENBROCK_BACK_BEGIN;
    return casadi_rosenbrock_backward(d, r);
  case ROSENBROCK_BACK_ADJ_BEGIN:
    if (d->fail) break;
    r->h_gamma = d->h * rp->gamma;
    if (casadi_rosenbrock_factorize(d, r)) break;
    casadi_ivp_back_step_begin(d);
    d->next = ROSENBROCK_BACK_ADJ_STAGE;
    return casadi_rosenbrock_backward(d, r);
  case ROSENBROCK_BACK_ADJ_SPREAD:
    // Adjoint of y0 + sum_j a_ij U_j and sum_j (c_ij/h) M U_j to y0 and earlier stages
    if (d->fail) break;
    d->n_vjp++;
    i = d->stage;
    rUi = d->rU + i * p->nr;
    casadi_ivp_add_vjp(d, d->ry);
    for (j = 0; j < i; ++j) {
      rUj = d->rU + j * p->nr;
      casadi_axpy(p->nrx, p->a[i * s + j], d->vjp_x, rUj);
      casadi_axpy(p->nrz, p->a[i * s + j], d->vjp_z, rUj + p->nrx);
      cij = rp->c[i * s + j] / d->h;
      casadi_axpy(p->nrx, cij, rUi, rUj);
      casadi_axpy(p->nrq, cij, rUi + p->nrx + p->nrz, rUj + p->nrx + p->nrz);
    }
    d->next = ROSENBROCK_BACK_ADJ_STAGE;
    return casadi_rosenbrock_backward(d, r);
  case ROSENBROCK_BACK_ADJ_STAGE:
    if (d->stage == 0) {
      // Contribution through J(t0, y0), dF/dt(t0, y0)
      casadi_ivp_clear_vjp(d);
      d->t_arg = d->t0;
      d->y_arg = d->y0;
      return casadi_ivp_request(d, IVP_EVAL_JAC_ADJ, ROSENBROCK_BACK_STEP_END);
    }
    // Previous stage: rbar_i = A^{-T} Ubar_i, then through J, dF/dt at y0 and F at the stage
    i = --d->stage;
    rUi = d->rU + i * p->nr;
    casadi_rosenbrock_solve_trans(d, r, rUi);
    casadi_rosenbrock_accum_jac_adj(d, rUi, d->U + i * p->ny, rp->gamma_sum[i] * d->h);
    casadi_ivp_stage_argument(d, i);
    casadi_ivp_clear_vjp(d);
    d->t_arg = d->t0 + p->alpha[i] * d->h;
    d->y_arg = d->y_stage;
    d->w = rUi;
    return casadi_ivp_request(d, IVP_EVAL_VJP, ROSENBROCK_BACK_ADJ_SPREAD);
  case ROSENBROCK_BACK_STEP_END:
    if (d->fail) break;
    casadi_ivp_add_vjp(d, d->ry);
    casadi_clear(d->jac_adj, p->nnz_jac * p->nadj * (1 + p->nfwd));
    d->t = d->t0;
    d->next = ROSENBROCK_BACK_BEGIN;
    return casadi_rosenbrock_backward(d, r);
  // Adjoint seed of the output at t_seed
  case IVP_IMPULSE:
    if (!casadi_ivp_impulse_begin(d)) return 0;
    // Output inside the pending step
    if (p->nz == 0) {
      casadi_ivp_seed_stages(d);
      return casadi_ivp_stop(d, IVP_SUCCESS);
    }
    // Adjoint of the Newton step on the algebraic states, as in the forward integration
    d->t_arg = d->t0;
    d->y_arg = d->y0;
    return casadi_ivp_request(d, IVP_EVAL_JAC, ROSENBROCK_IMPULSE_LINEARIZE);
  case ROSENBROCK_IMPULSE_LINEARIZE:
    if (d->fail || casadi_rosenbrock_factorize_alg(d, r)) break;
    d->t_arg = d->t_seed;
    d->y_arg = d->y_out;
    d->res = d->y_stage;
    return casadi_ivp_request(d, IVP_EVAL_RHS, ROSENBROCK_IMPULSE_CORRECT);
  case ROSENBROCK_IMPULSE_CORRECT:
    // z_out = z - G^{-1} alg with G = d(alg)/dz at the step start: dz = G^{-1} alg,
    // mu = G^{-T} zbar, Gbar = mu dz^T, (x, z, p, u)bar -= d(alg)/d(x, z, p, u)^T mu
    if (d->fail) break;
    casadi_rosenbrock_newton_step(d, r, d->y_stage + p->nx);
    casadi_clear(r->v1, p->nr);
    casadi_copy(d->seed + p->nrx, p->nrz, r->v1 + p->nrx);
    casadi_rosenbrock_newton_step_trans(d, r, r->v1 + p->nrx);
    casadi_clear(r->v2, p->ny);
    casadi_copy(d->y_stage + p->nx, p->nz, r->v2 + p->nx);
    casadi_rosenbrock_accum_jac_adj(d, r->v1, r->v2, 0.);
    casadi_scal(p->nrz, -1., r->v1 + p->nrx);
    casadi_ivp_clear_vjp(d);
    d->w = r->v1;
    return casadi_ivp_request(d, IVP_EVAL_VJP, ROSENBROCK_IMPULSE_END);
  case ROSENBROCK_IMPULSE_END:
    if (d->fail) break;
    d->n_vjp++;
    casadi_axpy(p->nrx, 1., d->vjp_x, d->seed);
    casadi_axpy(p->nrz, 1., d->vjp_z, d->seed + p->nrx);
    casadi_axpy(p->nrp, 1., d->vjp_p, d->adj_p);
    // Control adjoints belong to the output interval of the step
    casadi_axpy(p->nru, 1., d->vjp_u, d->ru_step);
    casadi_ivp_seed_stages(d);
    return casadi_ivp_stop(d, IVP_SUCCESS);
  default:
    return 0;
  }
  return casadi_ivp_stop(d, IVP_EVAL_FAILED);
}

template<typename T1>
int casadi_rosenbrock(casadi_ivp_data<T1>* d, casadi_rosenbrock_data<T1>* r) {
  // Reverse communication: returns 1 while a task in d->task is to be carried out
  if (d->next == IVP_IMPULSE || d->next == IVP_RETREAT || d->next >= ROSENBROCK_BACK_BEGIN) {
    return casadi_rosenbrock_backward(d, r);
  }
  return casadi_rosenbrock_forward(d, r);
}

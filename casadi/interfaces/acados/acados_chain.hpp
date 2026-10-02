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

// Chaining of acados sim solves over nt output intervals, with first and second order
// sensitivities w.r.t. w = [x0; u_0; ...; u_{nt-1}; p_d], where the differentiable
// parameters p_d enter acados as extra controls: acados u = [u_k; p_d] on interval k.
// Shared by the VM and generated code.

// C-REPLACE "casadi_acados_chain_prob<T1>" "struct casadi_acados_chain_prob"
// C-REPLACE "casadi_acados_chain_data<T1>" "struct casadi_acados_chain_data"
// C-REPLACE "casadi_acados_chain_solve<T1>" "casadi_acados_chain_solve"

template<typename T1>
struct casadi_acados_chain_prob {
  // Dimensions: states, controls (per interval), parameters, differentiable parameters
  // (0 or np), intervals
  casadi_int nx, nu, np, npd, nt;
  // Interval lengths
  const T1* T;
  // Work for first (jacobian), second order (hessian) sensitivities
  int jac, hess;
  // Setup-derived: acados w (nx + nu + npd), chain w (nx + nu*nt + npd), zero vector length
  casadi_int nw, nW, nzero;
};

template<typename T1>
struct casadi_acados_chain_data {
  const casadi_acados_chain_prob<T1>* prob;
  // One acados solve over an interval (any of xf, S_forw, S_adj, S_hess may be null)
  int (*solve)(void* ctx, T1 T, const T1* x0, const T1* p, const T1* u, T1* xf,
    T1* S_forw, const T1* seed_adj, T1* S_adj, T1* S_hess);
  void* ctx;
  // Work: states x_0..x_nt, acados controls, S_forw, S_adj, S_hess of one interval,
  // J = d[x_1; ...; x_nt]/dw, H = Hessian w.r.t. w, Hessian times W, adjoint, zeros
  T1 *xs, *ua, *S_forw, *S_adj, *S_hess, *J, *H, *tmp, *lam, *zero;
};

// SYMBOL "acados_chain_setup"
template<typename T1>
void casadi_acados_chain_setup(casadi_acados_chain_prob<T1>* p) {
  p->nw = p->nx + p->nu + p->npd;
  p->nW = p->nx + p->nu * p->nt + p->npd;
  p->nzero = p->nw;
  if (p->nu * p->nt > p->nzero) p->nzero = p->nu * p->nt;
  if (p->np > p->nzero) p->nzero = p->np;
}

// SYMBOL "acados_chain_work"
template<typename T1>
void casadi_acados_chain_work(const casadi_acados_chain_prob<T1>* p, casadi_int* sz_w) {
  casadi_int nx = p->nx, nw = p->nw, nt = p->nt, nW = p->nW;
  *sz_w += nx * (nt + 1) + (nw - nx) + nw + nx + p->nzero;
  if (p->jac) *sz_w += nx * nw + nx * nt * nW;
  if (p->hess) *sz_w += nw * nw + nW * nW + nw * nW;
}

// SYMBOL "acados_chain_set_work"
template<typename T1>
void casadi_acados_chain_set_work(casadi_acados_chain_data<T1>* d, const T1*** arg, T1*** res,
    casadi_int** iw, T1** w) {
  const casadi_acados_chain_prob<T1>* p = d->prob;
  casadi_int nx = p->nx, nw = p->nw, nt = p->nt, nW = p->nW;
  d->xs = *w; *w += nx * (nt + 1);
  d->ua = *w; *w += nw - nx;
  d->S_adj = *w; *w += nw;
  d->lam = *w; *w += nx;
  d->S_forw = d->J = d->S_hess = d->H = d->tmp = 0;
  if (p->jac) {
    d->S_forw = *w; *w += nx * nw;
    d->J = *w; *w += nx * nt * nW;
  }
  if (p->hess) {
    d->S_hess = *w; *w += nw * nw;
    d->H = *w; *w += nW * nW;
    d->tmp = *w; *w += nw * nW;
  }
  d->zero = *w; *w += p->nzero;
  casadi_clear(d->zero, p->nzero);
}

// C (nr-by-nc, leading dimension ldc) = A (nr-by-nk) * B (nk-by-nc, leading dimension ldb)
// SYMBOL "acados_chain_mm"
template<typename T1>
void casadi_acados_chain_mm(const T1* A, casadi_int nr, casadi_int nk, const T1* B,
    casadi_int ldb, casadi_int nc, T1* C, casadi_int ldc) {
  casadi_int r, c, i;
  for (c = 0; c < nc; ++c) {
    for (r = 0; r < nr; ++r) C[r + c * ldc] = 0;
    for (i = 0; i < nk; ++i) {
      for (r = 0; r < nr; ++r) C[r + c * ldc] += A[r + i * nr] * B[i + c * ldb];
    }
  }
}

// One acados solve over interval k, with acados controls [u_k; p_d]
// SYMBOL "acados_chain_solve"
template<typename T1>
int casadi_acados_chain_solve(casadi_acados_chain_data<T1>* d, casadi_int k, const T1* p,
    const T1* u, T1* xf, T1* S_forw, const T1* seed_adj, T1* S_adj, T1* S_hess) {
  const casadi_acados_chain_prob<T1>* q = d->prob;
  casadi_copy(u + k * q->nu, q->nu, d->ua);
  casadi_copy(p, q->npd, d->ua + q->nu);
  return d->solve(d->ctx, q->T[k], d->xs + k * q->nx, p, d->ua, xf, S_forw, seed_adj, S_adj,
    S_hess);
}

// Row block k of J from S_forw = [A B C] of interval k (blocks 0..k-1 already filled)
// SYMBOL "acados_chain_jac_step"
template<typename T1>
void casadi_acados_chain_jac_step(casadi_acados_chain_data<T1>* d, casadi_int k) {
  const casadi_acados_chain_prob<T1>* p = d->prob;
  casadi_int nx = p->nx, nu = p->nu, npd = p->npd, nr = p->nx * p->nt, nW = p->nW, c, r;
  casadi_int pc = nW - npd;
  const T1 *A = d->S_forw, *B = d->S_forw + nx * nx, *C = d->S_forw + nx * (nx + nu);
  T1* Jk = d->J + k * nx;
  if (k == 0) {
    // d x_1 / d (x0, p) = (A, C)
    for (c = 0; c < nx; ++c) for (r = 0; r < nx; ++r) Jk[r + c * nr] = A[r + c * nx];
    for (c = 0; c < npd; ++c) for (r = 0; r < nx; ++r) Jk[r + (pc + c) * nr] = C[r + c * nx];
  } else {
    // d x_{k+1} / d (x0, u_0..u_{k-1}) = A * d x_k / d (x0, u_0..u_{k-1})
    casadi_acados_chain_mm(A, nx, nx, Jk - nx, nr, nx + k * nu, Jk, nr);
    // d x_{k+1} / d p = A * d x_k / d p + C
    casadi_acados_chain_mm(A, nx, nx, Jk - nx + pc * nr, nr, npd, Jk + pc * nr, nr);
    for (c = 0; c < npd; ++c) for (r = 0; r < nx; ++r) Jk[r + (pc + c) * nr] += C[r + c * nx];
  }
  // d x_{k+1} / d u_k = B, independent of u_{k+1}..
  for (c = 0; c < nu; ++c) {
    for (r = 0; r < nx; ++r) Jk[r + (nx + k * nu + c) * nr] = B[r + c * nx];
  }
  for (c = nx + (k + 1) * nu; c < pc; ++c) for (r = 0; r < nx; ++r) Jk[r + c * nr] = 0;
}

// Nominal: xf (nx-by-nt, may be null) = [x_1, ..., x_nt]
// SYMBOL "acados_chain_nom"
template<typename T1>
int casadi_acados_chain_nom(casadi_acados_chain_data<T1>* d, const T1* x0, const T1* p,
    const T1* u, T1* xf) {
  casadi_int nx = d->prob->nx, k;
  casadi_copy(x0, nx, d->xs);
  for (k = 0; k < d->prob->nt; ++k) {
    if (casadi_acados_chain_solve<T1>(d, k, p, u, d->xs + (k + 1) * nx, 0, 0, 0, 0)) return 1;
  }
  casadi_copy(d->xs + nx, nx * d->prob->nt, xf);
  return 0;
}

// States x_0..x_nk, and if jac: row blocks 0..nk-1 of J
// SYMBOL "acados_chain_fwd"
template<typename T1>
int casadi_acados_chain_fwd(casadi_acados_chain_data<T1>* d, const T1* x0, const T1* p,
    const T1* u, casadi_int nk, int jac) {
  casadi_int nx = d->prob->nx, k;
  casadi_copy(x0, nx, d->xs);
  for (k = 0; k < nk; ++k) {
    if (casadi_acados_chain_solve<T1>(d, k, p, u, d->xs + (k + 1) * nx, jac ? d->S_forw : 0,
      0, 0, 0)) return 1;
    if (jac) casadi_acados_chain_jac_step(d, k);
  }
  return 0;
}

// Jacobian: J = [d xf / d x0, d xf / d u, d xf / d p_d], (nx*nt)-by-nW in d->J
// SYMBOL "acados_chain_jac"
template<typename T1>
int casadi_acados_chain_jac(casadi_acados_chain_data<T1>* d, const T1* x0, const T1* p,
    const T1* u) {
  return casadi_acados_chain_fwd(d, x0, p, u, d->prob->nt, 1);
}

// One adjoint direction, given the states x_0..x_{nt-1}: seed (nx-by-nt, may be null)
// -> adj_x0 (nx), adj_u (nu-by-nt), adj_p (npd); any output may be null
// SYMBOL "acados_chain_adj"
template<typename T1>
int casadi_acados_chain_adj(casadi_acados_chain_data<T1>* d, const T1* p, const T1* u,
    const T1* seed, T1* adj_x0, T1* adj_u, T1* adj_p) {
  casadi_int nx = d->prob->nx, nu = d->prob->nu, npd = d->prob->npd, k;
  casadi_clear(d->lam, nx);
  casadi_clear(adj_p, npd);
  for (k = d->prob->nt; k-- > 0; ) {
    if (seed) casadi_axpy(nx, 1., seed + k * nx, d->lam);
    if (casadi_acados_chain_solve<T1>(d, k, p, u, 0, 0, d->lam, d->S_adj, 0)) return 1;
    casadi_copy(d->S_adj, nx, d->lam);
    casadi_copy(d->S_adj + nx, nu, adj_u ? adj_u + k * nu : 0);
    if (adj_p) casadi_axpy(npd, 1., d->S_adj + nx + nu, adj_p);
  }
  casadi_copy(d->lam, nx, adj_x0);
  return 0;
}

// Hessian of sum_k seed_k^T x_{k+1} w.r.t. w in d->H, nW-by-nW,
// and the Jacobian J (for the derivatives w.r.t. the seed)
// SYMBOL "acados_chain_hess"
template<typename T1>
int casadi_acados_chain_hess(casadi_acados_chain_data<T1>* d, const T1* x0, const T1* p,
    const T1* u, const T1* seed) {
  casadi_int nx = d->prob->nx, nu = d->prob->nu, nt = d->prob->nt, nw = d->prob->nw;
  casadi_int nr = nx * nt, nW = d->prob->nW, pc = nW - d->prob->npd, k, i, r, c, uc;
  T1 s;
  const T1* Hk = d->S_hess;
  // States and Jacobian rows, except for the last interval
  if (casadi_acados_chain_fwd(d, x0, p, u, nt - 1, 1)) return 1;
  casadi_clear(d->H, nW * nW);
  casadi_clear(d->lam, nx);
  for (k = nt; k-- > 0; ) {
    // Weight on x_{k+1}: seed plus adjoint from the later intervals
    if (seed) casadi_axpy(nx, 1., seed + k * nx, d->lam);
    if (casadi_acados_chain_solve<T1>(d, k, p, u, 0, d->S_forw, d->lam, d->S_adj, d->S_hess))
      return 1;
    if (k == nt - 1) casadi_acados_chain_jac_step(d, k);
    casadi_copy(d->S_adj, nx, d->lam);
    if (nt == 1) {
      // W = I
      casadi_copy(Hk, nw * nw, d->H);
      continue;
    }
    // W = d(x_k, u_k, p_d)/dw: x part from J (identity for k == 0), u part selects u_k,
    // p part selects p_d
    uc = nx + k * nu;
    // tmp = Hk * W, nw-by-nW
    for (c = 0; c < nW; ++c) {
      for (r = 0; r < nw; ++r) {
        s = 0;
        if (k == 0) {
          if (c < nx) s = Hk[r + c * nw];
        } else {
          for (i = 0; i < nx; ++i) s += Hk[r + i * nw] * d->J[(k - 1) * nx + i + c * nr];
        }
        if (c >= uc && c < uc + nu) s += Hk[r + (nx + c - uc) * nw];
        if (c >= pc) s += Hk[r + (nx + nu + c - pc) * nw];
        d->tmp[r + c * nw] = s;
      }
    }
    // H += W^T * tmp
    for (c = 0; c < nW; ++c) {
      for (r = 0; r < nW; ++r) {
        s = 0;
        if (k == 0) {
          if (r < nx) s = d->tmp[r + c * nw];
        } else {
          for (i = 0; i < nx; ++i) s += d->J[(k - 1) * nx + i + r * nr] * d->tmp[i + c * nw];
        }
        if (r >= uc && r < uc + nu) s += d->tmp[nx + r - uc + c * nw];
        if (r >= pc) s += d->tmp[nx + nu + r - pc + c * nw];
        d->H[r + c * nW] += s;
      }
    }
  }
  return 0;
}

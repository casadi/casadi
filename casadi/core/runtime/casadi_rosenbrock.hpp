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

// Adaptive Rosenbrock-Wanner integrator for M y' = F(t, y), y = [x; z; q], M = diag(I, 0, I)
// Step: (M/(h gamma) - J) U_i = F(t0 + alpha_i h, y0 + sum_j a_ij U_j) + sum_j (c_ij/h) M U_j
//         + gamma_i h dF/dt, y1 = y0 + sum_i m_i U_i, error sum_i m_err_i U_i
// Steps depend on the nominal entries only: sensitivities are exact derivatives of the scheme
// Layouts: x = [x_0; ...; x_nfwd] (nominal, forward), rx = [rx_0 (nx1 x nadj); ...; rx_nfwd]

// C-REPLACE "fmin" "casadi_fmin"
// C-REPLACE "fmax" "casadi_fmax"
// C-REPLACE "fabs" "casadi_fabs"
// C-REPLACE "static_cast<T1>" "(casadi_real) "
// C-REPLACE "static_cast<casadi_int*>(0)" "0"
// C-REPLACE "static_cast<int>" "(int) "

// SYMBOL "rosenbrock_jac_t"
typedef enum {
  // Jacobian blocks w.r.t. x, z and t, in the order of evaluation by the caller
  ROSENBROCK_ODE_X, ROSENBROCK_ALG_X, ROSENBROCK_QUAD_X,
  ROSENBROCK_ODE_Z, ROSENBROCK_ALG_Z, ROSENBROCK_QUAD_Z,
  ROSENBROCK_ODE_T, ROSENBROCK_ALG_T, ROSENBROCK_QUAD_T,
  ROSENBROCK_NUM_JAC
} casadi_rosenbrock_jac_t;

// SYMBOL "rosenbrock_prob"
template<typename T1>
struct casadi_rosenbrock_prob {
  // Nominal block sizes: differential, algebraic, quadrature states, parameters, controls
  casadi_int nx1, nz1, nq1, np1, nu1;
  // Number of forward and adjoint directions
  casadi_int nfwd, nadj;
  // Leading entries of x, z, q subject to error control
  casadi_int nx_err, nz_err, nq_err;
  // Leading controls whose change triggers a consistent initialization
  casadi_int nu_err;
  // Number of stages, order of the error estimate
  casadi_int n_stages, order_err;
  // Diagonal coefficient
  T1 gamma;
  // Method coefficients, transformed form: a, c are strictly lower triangular, row-major
  const T1 *a, *c, *alpha, *gamma_sum, *m, *m_err;
  // Dense output weights m_i(theta) = sum_q m_dense[q * n_stages + i] theta^(q+1)
  casadi_int n_dense;
  const T1* m_dense;
  // Sparsity patterns of the Jacobian blocks, see casadi_rosenbrock_jac_t
  const casadi_int* sp_jac[ROSENBROCK_NUM_JAC];
  // Iteration matrix M/(h gamma) - J: sparsity, positions of ODE_X, ALG_X, ODE_Z, ALG_Z, diag
  const casadi_int *sp_iter, *iter_map;
  // Sparse QR of the iteration matrix and of d(alg)/dz
  const casadi_int *sp_iter_v, *sp_iter_r, *iter_prinv, *iter_pc;
  const casadi_int *sp_alg_v, *sp_alg_r, *alg_prinv, *alg_pc;
  // Error tolerances
  T1 abstol, reltol;
  // Step size controller
  T1 safety, fac_min, fac_max, h_min, h_max, h_init;
  // Maximum number of steps per output interval
  casadi_int max_num_steps;
  // Number of equidistant steps per output interval (0: adaptive)
  casadi_int fixed_steps;
  // Consistent initialization: tolerance on scaled Newton step, iteration limit
  T1 init_tol;
  casadi_int init_max_iter;
  // Capacity of the tape in the work vector (0: the caller provides the tape)
  casadi_int tape_size;
  // Derived: forward and adjoint state dimensions, controls, vector lengths
  casadi_int nx, nz, nq, nu, nrx, nrz, nrq, nrp, nru, ny, nr;
  // Derived: offsets of the Jacobian blocks in a buffer with one direction
  casadi_int jac_off[ROSENBROCK_NUM_JAC + 1];
  // Derived: nonzeros of all Jacobian blocks, of the iteration matrix and its factors
  casadi_int nnz_jac, nnz_iter, nnz_iter_v, nnz_iter_r, nnz_alg_v, nnz_alg_r;
  // Derived: length of the work vectors
  casadi_int nv, n_qr_w;
};
// C-REPLACE "casadi_rosenbrock_prob<T1>" "struct casadi_rosenbrock_prob"

// SYMBOL "rosenbrock_nnz"
inline
casadi_int casadi_rosenbrock_nnz(const casadi_int* sp) {
  // Number of nonzeros of a sparsity pattern
  return sp[2 + sp[1]];
}

// SYMBOL "rosenbrock_isfinite"
template<typename T1>
int casadi_rosenbrock_isfinite(T1 x) {
  // False for +-inf and NaN
  return x - x == 0;
}

// SYMBOL "rosenbrock_setup"
template<typename T1>
void casadi_rosenbrock_setup(casadi_rosenbrock_prob<T1>* p) {
  casadi_int k, n1;
  // Forward states including sensitivities
  p->nx = p->nx1 * (1 + p->nfwd);
  p->nz = p->nz1 * (1 + p->nfwd);
  p->nq = p->nq1 * (1 + p->nfwd);
  p->nu = p->nu1 * (1 + p->nfwd);
  // Adjoint states and adjoint parameters, controls, including forward sensitivities
  p->nrx = p->nx1 * p->nadj * (1 + p->nfwd);
  p->nrz = p->nz1 * p->nadj * (1 + p->nfwd);
  p->nrq = p->nq1 * p->nadj * (1 + p->nfwd);
  p->nrp = p->np1 * p->nadj * (1 + p->nfwd);
  p->nru = p->nu1 * p->nadj * (1 + p->nfwd);
  p->ny = p->nx + p->nz + p->nq;
  p->nr = p->nrx + p->nrz + p->nrq;
  // Jacobian blocks
  p->jac_off[0] = 0;
  for (k = 0; k < ROSENBROCK_NUM_JAC; ++k) {
    p->jac_off[k + 1] = p->jac_off[k] + casadi_rosenbrock_nnz(p->sp_jac[k]);
  }
  p->nnz_jac = p->jac_off[ROSENBROCK_NUM_JAC];
  p->nnz_iter = casadi_rosenbrock_nnz(p->sp_iter);
  p->nnz_iter_v = casadi_rosenbrock_nnz(p->sp_iter_v);
  p->nnz_iter_r = casadi_rosenbrock_nnz(p->sp_iter_r);
  p->nnz_alg_v = p->nz1 > 0 ? casadi_rosenbrock_nnz(p->sp_alg_v) : 0;
  p->nnz_alg_r = p->nz1 > 0 ? casadi_rosenbrock_nnz(p->sp_alg_r) : 0;
  // Work vectors: forward and adjoint vectors with all directions
  n1 = p->nx1 + p->nz1;
  p->nv = p->ny;
  if (p->nr + p->nrp + p->nru > p->nv) p->nv = p->nr + p->nrp + p->nru;
  if (n1 * p->nadj * (1 + p->nfwd) > p->nv) p->nv = n1 * p->nadj * (1 + p->nfwd);
  // Work vector of the sparse QR: extended number of rows
  p->n_qr_w = p->sp_iter_v[0];
  if (p->nz1 > 0 && p->sp_alg_v[0] > p->n_qr_w) p->n_qr_w = p->sp_alg_v[0];
}

// SYMBOL "rosenbrock_task_t"
typedef enum {
  // res = F(t_arg, y_arg), with forward sensitivities
  ROSENBROCK_EVAL_RHS,
  // jac, jac_fwd = Jacobian blocks at (t_arg, y_arg), with forward sensitivities
  ROSENBROCK_EVAL_JAC,
  // vjp = dF/d(x,z,p,u)^T w at (t_arg, y_arg), with forward sensitivities
  ROSENBROCK_EVAL_VJP,
  // vjp = d/d(x,z,p,u) <jac_adj, Jacobian blocks> at (t_arg, y_arg), w. forward sensitivities
  ROSENBROCK_EVAL_JAC_ADJ,
  // Enlarge the tape to at least tape_need entries, keeping its contents; set tape, tape_cap
  ROSENBROCK_TAPE_RESERVE
} casadi_rosenbrock_task_t;

// SYMBOL "rosenbrock_next_t"
typedef enum {
  // Phases: each one consumes the result of the previously requested task
  ROSENBROCK_DONE,
  ROSENBROCK_INIT_BEGIN,
  ROSENBROCK_INIT_LINEARIZE,
  ROSENBROCK_INIT_NEWTON_STEP,
  ROSENBROCK_INIT_SENS_RHS,
  ROSENBROCK_INIT_SENS,
  ROSENBROCK_INIT_END,
  ROSENBROCK_STEP_BEGIN,
  ROSENBROCK_STEP_INITIAL_SIZE,
  ROSENBROCK_STEP_FACTORIZE,
  ROSENBROCK_STEP_REPLAY,
  ROSENBROCK_STAGE_SOLVE,
  ROSENBROCK_OUTPUT_CORRECT,
  ROSENBROCK_BACK_BEGIN,
  ROSENBROCK_BACK_INIT_SOLVE,
  ROSENBROCK_BACK_INIT_END,
  ROSENBROCK_BACK_ADJ_BEGIN,
  ROSENBROCK_BACK_ADJ_STAGE,
  ROSENBROCK_BACK_ADJ_SPREAD,
  ROSENBROCK_BACK_STEP_END,
  ROSENBROCK_IMPULSE_BEGIN,
  ROSENBROCK_IMPULSE_LINEARIZE,
  ROSENBROCK_IMPULSE_CORRECT,
  ROSENBROCK_IMPULSE_END
} casadi_rosenbrock_next_t;

// SYMBOL "rosenbrock_flag_t"
typedef enum {
  ROSENBROCK_SUCCESS,
  ROSENBROCK_MAX_NUM_STEPS,
  ROSENBROCK_STEP_TOO_SMALL,
  ROSENBROCK_INIT_FAILED,
  ROSENBROCK_EVAL_FAILED,
  ROSENBROCK_TAPE_FAILED
} casadi_rosenbrock_flag_t;

// SYMBOL "rosenbrock_record_t"
typedef enum {
  // Integrator step from t0 to t0 + h
  ROSENBROCK_RECORD_STEP,
  // Consistent initialization of the algebraic states at t0
  ROSENBROCK_RECORD_INIT
} casadi_rosenbrock_record_t;

// C-REPLACE "casadi_rosenbrock_data<T1>" "struct casadi_rosenbrock_data"
// SYMBOL "rosenbrock_data"
template<typename T1>
struct casadi_rosenbrock_data {
  // Problem structure
  const casadi_rosenbrock_prob<T1>* prob;
  // Nominal integration whose steps are replayed (or null), data whose task is requested
  casadi_rosenbrock_data<T1> *nom, *active;
  // Stop after every accepted step and skip the outputs (nominal integration)
  int pause;
  // Task for the caller, where to resume, solver status
  casadi_rosenbrock_task_t task;
  casadi_rosenbrock_next_t next;
  casadi_rosenbrock_flag_t status;
  // Set by the caller if a task could not be carried out
  int fail;
  // Task arguments and results
  T1 t_arg;
  const T1 *y_arg, *w;
  T1 *res;
  // Current time and state, output time, latest time to step to, start of a retreat
  T1 t, t_out, t_stop, t_start;
  T1* y;
  // State at the output time
  T1* y_out;
  // Controls set by the caller, controls of the current or popped step
  const T1* u;
  T1* u_step;
  // Current step: start time, size, h * gamma, start state, stage counter
  T1 t0, h, h_gamma;
  // End time and end of the output interval of a popped record
  T1 t1, t_record;
  T1* y0;
  casadi_int stage;
  // Stage increments U_i (n_stages x ny), stage argument, time derivative of F
  T1 *U, *y_stage, *dfdt;
  // Adjoint state, stage adjoints (n_stages x nr), adjoint of the start of the popped step
  T1 *ry, *rU, *ry0;
  // Adjoint seed of an output at time t_seed, dense output weights
  T1 *seed, *w_dense;
  T1 t_seed;
  // Control adjoints of the popped step from outputs inside it
  T1* ru_step;
  // Jacobian blocks and their forward derivatives, nonzeros in casadi_rosenbrock_jac_t order
  T1 *jac, *jac_fwd;
  // Adjoint seeds of the Jacobian blocks, with forward derivatives
  T1 *jac_adj, *jac_adj_fwd;
  // Iteration matrix and its QR factors, QR factors of d(alg)/dz
  T1 *iter, *iter_v, *iter_r, *iter_beta, *alg_v, *alg_r, *alg_beta, *qr_w;
  // Results of EVAL_VJP and EVAL_JAC_ADJ: adjoints of x, z, p, u
  T1 *vjp_x, *vjp_z, *vjp_p, *vjp_u;
  // Accumulated adjoints of the parameters and controls
  T1 *adj_p, *adj_u;
  // Work vectors of length nv
  T1 *v1, *v2;
  // Tape: records [U (steps only), y0, u, h, t0, t1, t_out, type]; entries used, capacity, needed
  T1* tape;
  casadi_int tape_pos, tape_cap, tape_need;
  // Step size controller: proposed next step, previous error
  T1 h_next, err_prev;
  int last_rejected, last_clipped, last_step;
  // Steps taken in the current output interval
  casadi_int n_steps_interval;
  // Consistent initialization pending, Newton iteration
  int init_pending;
  casadi_int init_iter;
  // Nominal controls at the last consistent initialization
  T1* u_init;
  // Taping requested, record type, popped record that belongs to an earlier interval
  int taping;
  casadi_rosenbrock_record_t record;
  int record_pending;
  // Statistics
  casadi_int n_steps, n_reject, n_rhs, n_jac, n_factor, n_solve, n_vjp;
};

// SYMBOL "rosenbrock_work"
template<typename T1>
void casadi_rosenbrock_work(const casadi_rosenbrock_prob<T1>* p, casadi_int* sz_arg,
    casadi_int* sz_res, casadi_int* sz_iw, casadi_int* sz_w) {
  casadi_int n1 = p->nx1 + p->nz1;
  // State: y, y_out, y0, stages, stage argument, dF/dt, dense output weights
  *sz_w += p->ny * (p->n_stages + 5) + p->n_stages;
  // Adjoint: adjoint state, stage adjoints, adjoint of the step start, seed
  *sz_w += p->nr * (p->n_stages + 3);
  // Controls of a step, their adjoints from outputs inside it
  *sz_w += p->nu + p->nru;
  // Jacobian blocks, forward derivatives, adjoint seeds and their forward derivatives
  *sz_w += p->nnz_jac * (1 + p->nfwd) * (1 + p->nadj);
  // Iteration matrix, QR factors of it and of d(alg)/dz, QR work vector
  *sz_w += p->nnz_iter + p->nnz_iter_v + p->nnz_iter_r + n1;
  *sz_w += p->nnz_alg_v + p->nnz_alg_r + p->nz1 + p->n_qr_w;
  // Function adjoints, accumulated adjoints
  *sz_w += p->nrx + p->nrz + 2 * (p->nrp + p->nru);
  // Work vectors, nominal controls at the last consistent initialization
  *sz_w += 2 * p->nv + p->nu_err;
  // Tape
  *sz_w += p->tape_size;
}

// SYMBOL "rosenbrock_set_work"
template<typename T1>
void casadi_rosenbrock_set_work(casadi_rosenbrock_data<T1>* d, const T1*** arg, T1*** res,
    casadi_int** iw, T1** w) {
  const casadi_rosenbrock_prob<T1>* p = d->prob;
  casadi_int n1 = p->nx1 + p->nz1;
  d->y = *w; *w += p->ny;
  d->y_out = *w; *w += p->ny;
  d->y0 = *w; *w += p->ny;
  d->U = *w; *w += p->ny * p->n_stages;
  d->y_stage = *w; *w += p->ny;
  d->dfdt = *w; *w += p->ny;
  d->w_dense = *w; *w += p->n_stages;
  d->ry = *w; *w += p->nr;
  d->rU = *w; *w += p->nr * p->n_stages;
  d->ry0 = *w; *w += p->nr;
  d->seed = *w; *w += p->nr;
  d->u_step = *w; *w += p->nu;
  d->ru_step = *w; *w += p->nru;
  d->jac = *w; *w += p->nnz_jac;
  d->jac_fwd = *w; *w += p->nnz_jac * p->nfwd;
  d->jac_adj = *w; *w += p->nnz_jac * p->nadj;
  d->jac_adj_fwd = *w; *w += p->nnz_jac * p->nadj * p->nfwd;
  d->iter = *w; *w += p->nnz_iter;
  d->iter_v = *w; *w += p->nnz_iter_v;
  d->iter_r = *w; *w += p->nnz_iter_r;
  d->iter_beta = *w; *w += n1;
  d->alg_v = *w; *w += p->nnz_alg_v;
  d->alg_r = *w; *w += p->nnz_alg_r;
  d->alg_beta = *w; *w += p->nz1;
  d->qr_w = *w; *w += p->n_qr_w;
  d->vjp_x = *w; *w += p->nrx;
  d->vjp_z = *w; *w += p->nrz;
  d->vjp_p = *w; *w += p->nrp;
  d->vjp_u = *w; *w += p->nru;
  d->adj_p = *w; *w += p->nrp;
  d->adj_u = *w; *w += p->nru;
  d->v1 = *w; *w += p->nv;
  d->v2 = *w; *w += p->nv;
  d->u_init = *w; *w += p->nu_err;
  d->tape = *w; *w += p->tape_size;
  d->tape_cap = p->tape_size;
  d->tape_pos = 0;
  d->next = ROSENBROCK_DONE;
  d->status = ROSENBROCK_SUCCESS;
  d->fail = 0;
  d->taping = 0;
  d->record_pending = 0;
  d->nom = 0;
  d->active = d;
  d->pause = 0;
}

// SYMBOL "rosenbrock_clear_stats"
template<typename T1>
void casadi_rosenbrock_clear_stats(casadi_rosenbrock_data<T1>* d) {
  d->n_steps = d->n_reject = d->n_rhs = d->n_jac = d->n_factor = d->n_solve = d->n_vjp = 0;
}

// SYMBOL "rosenbrock_reset"
template<typename T1>
void casadi_rosenbrock_reset(casadi_rosenbrock_data<T1>* d, T1 t, const T1* u,
    int first_call) {
  // Start at time t, controls u; the state d->y is set by the caller
  const casadi_rosenbrock_prob<T1>* p = d->prob;
  casadi_int i;
  int init_alg = first_call;
  d->t = t;
  d->u = u;
  if (first_call) {
    // Every evaluation starts from the same controller state
    d->h_next = 0;
    d->err_prev = 1;
    d->last_rejected = 0;
    d->tape_pos = 0;
    casadi_rosenbrock_clear_stats(d);
  }
  // Consistent algebraic states at the start and when the nominal controls change
  for (i = 0; i < p->nu_err; ++i) init_alg = init_alg || u[i] != d->u_init[i];
  if (init_alg) casadi_copy(u, p->nu_err, d->u_init);
  d->init_pending = init_alg && p->nz > 0;
  if (d->nom) {
    // Initial nominal states: the leading entries of x, z and q; afterwards, the nominal
    // integration continues its own trajectory
    const casadi_rosenbrock_prob<T1>* pn = d->nom->prob;
    if (first_call) {
      casadi_copy(d->y, pn->nx, d->nom->y);
      casadi_copy(d->y + p->nx, pn->nz, d->nom->y + pn->nx);
      casadi_copy(d->y + p->nx + p->nz, pn->nq, d->nom->y + pn->nx + pn->nz);
    }
    casadi_rosenbrock_reset(d->nom, t, u, first_call);
  }
}

// SYMBOL "rosenbrock_advance"
template<typename T1>
void casadi_rosenbrock_advance(casadi_rosenbrock_data<T1>* d, T1 t_out, T1 t_stop,
    int taping) {
  // Output at t_out, stepping no further than t_stop; fixed steps end at every output time
  d->t_out = t_out;
  d->t_stop = d->prob->fixed_steps > 0 ? t_out : t_stop;
  d->taping = taping;
  d->n_steps_interval = 0;
  d->status = ROSENBROCK_SUCCESS;
  d->fail = 0;
  d->next = d->init_pending ? ROSENBROCK_INIT_BEGIN : ROSENBROCK_STEP_BEGIN;
  if (d->nom) {
    casadi_rosenbrock_advance(d->nom, t_out, t_stop, 0);
    d->nom->pause = 1;
  }
}

// SYMBOL "rosenbrock_next_stop"
template<typename T1>
casadi_int casadi_rosenbrock_next_stop(const T1* u, casadi_int nu, casadi_int nu_nom,
    casadi_int k, casadi_int nt) {
  // Last output interval from k on before the nominal controls change; u is column k
  casadi_int i;
  if (u == 0) return nt - 1;
  for (; k + 1 < nt; ++k) {
    for (i = 0; i < nu_nom; ++i) if (u[i] != u[nu + i]) return k;
    u += nu;
  }
  return k;
}

// SYMBOL "rosenbrock_reset_backward"
template<typename T1>
void casadi_rosenbrock_reset_backward(casadi_rosenbrock_data<T1>* d) {
  const casadi_rosenbrock_prob<T1>* p = d->prob;
  casadi_clear(d->ry, p->nr);
  casadi_clear(d->adj_p, p->nrp);
  casadi_clear(d->adj_u, p->nru);
  casadi_clear(d->jac_adj, p->nnz_jac * p->nadj * (1 + p->nfwd));
  d->record_pending = 0;
}

// SYMBOL "rosenbrock_retreat"
template<typename T1>
void casadi_rosenbrock_retreat(casadi_rosenbrock_data<T1>* d, T1 t_start, T1 t_out) {
  // Adjoint from t_start back to t_out; skips records of output intervals ending after t_start
  d->t = d->t_start = t_start;
  d->t_out = t_out;
  d->status = ROSENBROCK_SUCCESS;
  d->next = ROSENBROCK_BACK_BEGIN;
}

// SYMBOL "rosenbrock_request"
template<typename T1>
int casadi_rosenbrock_request(casadi_rosenbrock_data<T1>* d, casadi_rosenbrock_task_t task,
    casadi_rosenbrock_next_t next) {
  // Hand a task to the caller, resume at phase next
  d->active = d;
  d->task = task;
  d->next = next;
  d->fail = 0;
  return 1;
}

// SYMBOL "rosenbrock_stop"
template<typename T1>
int casadi_rosenbrock_stop(casadi_rosenbrock_data<T1>* d, casadi_rosenbrock_flag_t status) {
  d->status = status;
  d->next = ROSENBROCK_DONE;
  return 0;
}

// ---------------------------------------------------------------------------------------
// Linear algebra of the sensitivity-augmented iteration matrix
// ---------------------------------------------------------------------------------------

// SYMBOL "rosenbrock_jac_nz"
template<typename T1>
T1* casadi_rosenbrock_jac_nz(const casadi_rosenbrock_data<T1>* d, T1* jac, casadi_int k,
    casadi_int ndir, casadi_int dir) {
  // Nonzeros of Jacobian block k, direction dir, in a buffer with ndir directions per block
  const casadi_int* off = d->prob->jac_off;
  return jac + ndir * off[k] + dir * (off[k + 1] - off[k]);
}

// SYMBOL "rosenbrock_jac_mv"
template<typename T1>
void casadi_rosenbrock_jac_mv(const casadi_rosenbrock_data<T1>* d, casadi_int k,
    casadi_int dir, const T1* v, T1* r, int tr) {
  // r += J_k v (or J_k^T v): Jacobian block k (dir < 0) or its forward derivative
  const casadi_rosenbrock_prob<T1>* p = d->prob;
  const T1* nz = dir < 0 ? casadi_rosenbrock_jac_nz(d, d->jac, k, 1, 0)
    : casadi_rosenbrock_jac_nz(d, d->jac_fwd, k, p->nfwd, dir);
  casadi_mv(nz, p->sp_jac[k], v, r, tr);
}

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
void casadi_rosenbrock_dfdt(casadi_rosenbrock_data<T1>* d) {
  // Time derivative of [ode; alg; quad], with sensitivities, from the time columns
  const casadi_rosenbrock_prob<T1>* p = d->prob;
  casadi_int k, dir, off, n1;
  T1 one = 1;
  casadi_clear(d->dfdt, p->ny);
  for (k = 0; k < 3; ++k) {
    // Empty if the right-hand side does not depend on time
    if (p->sp_jac[ROSENBROCK_ODE_T + k][1] == 0) continue;
    off = k == 0 ? 0 : k == 1 ? p->nx : p->nx + p->nz;
    n1 = k == 0 ? p->nx1 : k == 1 ? p->nz1 : p->nq1;
    for (dir = -1; dir < p->nfwd; ++dir) {
      casadi_rosenbrock_jac_mv(d, ROSENBROCK_ODE_T + k, dir, &one,
        d->dfdt + off + (dir + 1) * n1, 0);
    }
  }
}

// SYMBOL "rosenbrock_factorize"
template<typename T1>
int casadi_rosenbrock_factorize(casadi_rosenbrock_data<T1>* d) {
  // QR of M/(h gamma) - J for the nominal [x; z]; returns 1 if singular
  const casadi_rosenbrock_prob<T1>* p = d->prob;
  casadi_int k, i, n, pos = 0;
  const casadi_int blk[4] = {ROSENBROCK_ODE_X, ROSENBROCK_ALG_X, ROSENBROCK_ODE_Z,
    ROSENBROCK_ALG_Z};
  const T1* J;
  T1 rmin = 0;
  casadi_clear(d->iter, p->nnz_iter);
  for (k = 0; k < 4; ++k) {
    J = casadi_rosenbrock_jac_nz(d, d->jac, blk[k], 1, 0);
    n = casadi_rosenbrock_nnz(p->sp_jac[blk[k]]);
    for (i = 0; i < n; ++i) d->iter[p->iter_map[pos++]] -= J[i];
  }
  for (i = 0; i < p->nx1; ++i) d->iter[p->iter_map[pos++]] += 1 / d->h_gamma;
  casadi_qr(p->sp_iter, d->iter, d->qr_w, p->sp_iter_v, d->iter_v, p->sp_iter_r, d->iter_r,
    d->iter_beta, p->iter_prinv, p->iter_pc);
  casadi_qr_singular(&rmin, static_cast<casadi_int*>(0), d->iter_r, p->sp_iter_r,
    p->iter_pc, 0.);
  d->n_factor++;
  return !(rmin > 0) || !casadi_rosenbrock_isfinite(rmin);
}

// SYMBOL "rosenbrock_iter_solve"
template<typename T1>
void casadi_rosenbrock_iter_solve(casadi_rosenbrock_data<T1>* d, T1* v, casadi_int nrhs,
    int tr) {
  // Solve with the factorized iteration matrix (or its transpose), nrhs right-hand sides
  const casadi_rosenbrock_prob<T1>* p = d->prob;
  casadi_qr_solve(v, nrhs, tr, p->sp_iter_v, d->iter_v, p->sp_iter_r, d->iter_r, d->iter_beta,
    p->iter_prinv, p->iter_pc, d->qr_w);
}

// SYMBOL "rosenbrock_solve"
template<typename T1>
void casadi_rosenbrock_solve(casadi_rosenbrock_data<T1>* d, T1* v) {
  // A U_0 = r_0, A U_d = r_d + dJ_d U_0 (A = M/hg - J), quadratures U_q = hg (r_q + J_qy U)
  const casadi_rosenbrock_prob<T1>* p = d->prob;
  casadi_int dir, n1 = p->nx1 + p->nz1;
  T1 *u0 = d->v1, *ud, *vx = v, *vz = v + p->nx, *vq = v + p->nx + p->nz, *vqd;
  casadi_rosenbrock_gather(vx, vz, p->nx1, p->nz1, u0);
  casadi_rosenbrock_iter_solve(d, u0, 1, 0);
  for (dir = 0; dir < p->nfwd; ++dir) {
    ud = d->v2 + dir * n1;
    casadi_rosenbrock_gather(vx + (dir + 1) * p->nx1, vz + (dir + 1) * p->nz1,
      p->nx1, p->nz1, ud);
    casadi_rosenbrock_jac_mv(d, ROSENBROCK_ODE_X, dir, u0, ud, 0);
    casadi_rosenbrock_jac_mv(d, ROSENBROCK_ODE_Z, dir, u0 + p->nx1, ud, 0);
    casadi_rosenbrock_jac_mv(d, ROSENBROCK_ALG_X, dir, u0, ud + p->nx1, 0);
    casadi_rosenbrock_jac_mv(d, ROSENBROCK_ALG_Z, dir, u0 + p->nx1, ud + p->nx1, 0);
  }
  if (p->nfwd > 0) casadi_rosenbrock_iter_solve(d, d->v2, p->nfwd, 0);
  for (dir = 0; dir < p->nfwd; ++dir) {
    ud = d->v2 + dir * n1;
    vqd = vq + (dir + 1) * p->nq1;
    casadi_rosenbrock_jac_mv(d, ROSENBROCK_QUAD_X, -1, ud, vqd, 0);
    casadi_rosenbrock_jac_mv(d, ROSENBROCK_QUAD_Z, -1, ud + p->nx1, vqd, 0);
    casadi_rosenbrock_jac_mv(d, ROSENBROCK_QUAD_X, dir, u0, vqd, 0);
    casadi_rosenbrock_jac_mv(d, ROSENBROCK_QUAD_Z, dir, u0 + p->nx1, vqd, 0);
    casadi_scal(p->nq1, d->h_gamma, vqd);
    casadi_rosenbrock_scatter(ud, p->nx1, p->nz1, vx + (dir + 1) * p->nx1,
      vz + (dir + 1) * p->nz1);
  }
  casadi_rosenbrock_jac_mv(d, ROSENBROCK_QUAD_X, -1, u0, vq, 0);
  casadi_rosenbrock_jac_mv(d, ROSENBROCK_QUAD_Z, -1, u0 + p->nx1, vq, 0);
  casadi_scal(p->nq1, d->h_gamma, vq);
  casadi_rosenbrock_scatter(u0, p->nx1, p->nz1, vx, vz);
  d->n_solve++;
}

// SYMBOL "rosenbrock_solve_trans"
template<typename T1>
void casadi_rosenbrock_solve_trans(casadi_rosenbrock_data<T1>* d, T1* v) {
  // A^T mu = l + J_qy^T hg l_q, then A^T dmu_d = dl_d + J_qy^T hg dl_q + dJ_d^T [mu; hg l_q]
  const casadi_rosenbrock_prob<T1>* p = d->prob;
  casadi_int a, dir, n1 = p->nx1 + p->nz1;
  casadi_int nrx1 = p->nx1 * p->nadj, nrz1 = p->nz1 * p->nadj, nrq1 = p->nq1 * p->nadj;
  T1 *vx = v, *vz = v + p->nrx, *vq = v + p->nrx + p->nrz;
  T1 *mu, *mq, *dmu, *dmq, *dvx, *dvz;
  casadi_scal(p->nrq, d->h_gamma, vq);
  for (a = 0; a < p->nadj; ++a) {
    mu = d->v1 + a * n1;
    mq = vq + a * p->nq1;
    casadi_rosenbrock_gather(vx + a * p->nx1, vz + a * p->nz1, p->nx1, p->nz1, mu);
    casadi_rosenbrock_jac_mv(d, ROSENBROCK_QUAD_X, -1, mq, mu, 1);
    casadi_rosenbrock_jac_mv(d, ROSENBROCK_QUAD_Z, -1, mq, mu + p->nx1, 1);
  }
  casadi_rosenbrock_iter_solve(d, d->v1, p->nadj, 1);
  for (dir = 0; dir < p->nfwd; ++dir) {
    for (a = 0; a < p->nadj; ++a) {
      mu = d->v1 + a * n1;
      mq = vq + a * p->nq1;
      dmu = d->v2 + (dir * p->nadj + a) * n1;
      dmq = vq + (dir + 1) * nrq1 + a * p->nq1;
      dvx = vx + (dir + 1) * nrx1 + a * p->nx1;
      dvz = vz + (dir + 1) * nrz1 + a * p->nz1;
      casadi_rosenbrock_gather(dvx, dvz, p->nx1, p->nz1, dmu);
      casadi_rosenbrock_jac_mv(d, ROSENBROCK_QUAD_X, -1, dmq, dmu, 1);
      casadi_rosenbrock_jac_mv(d, ROSENBROCK_QUAD_Z, -1, dmq, dmu + p->nx1, 1);
      casadi_rosenbrock_jac_mv(d, ROSENBROCK_QUAD_X, dir, mq, dmu, 1);
      casadi_rosenbrock_jac_mv(d, ROSENBROCK_QUAD_Z, dir, mq, dmu + p->nx1, 1);
      casadi_rosenbrock_jac_mv(d, ROSENBROCK_ODE_X, dir, mu, dmu, 1);
      casadi_rosenbrock_jac_mv(d, ROSENBROCK_ALG_X, dir, mu + p->nx1, dmu, 1);
      casadi_rosenbrock_jac_mv(d, ROSENBROCK_ODE_Z, dir, mu, dmu + p->nx1, 1);
      casadi_rosenbrock_jac_mv(d, ROSENBROCK_ALG_Z, dir, mu + p->nx1, dmu + p->nx1, 1);
    }
  }
  if (p->nfwd > 0) casadi_rosenbrock_iter_solve(d, d->v2, p->nadj * p->nfwd, 1);
  for (a = 0; a < p->nadj; ++a) {
    casadi_rosenbrock_scatter(d->v1 + a * n1, p->nx1, p->nz1, vx + a * p->nx1,
      vz + a * p->nz1);
    for (dir = 0; dir < p->nfwd; ++dir) {
      casadi_rosenbrock_scatter(d->v2 + (dir * p->nadj + a) * n1, p->nx1, p->nz1,
        vx + (dir + 1) * nrx1 + a * p->nx1, vz + (dir + 1) * nrz1 + a * p->nz1);
    }
  }
  d->n_solve++;
}

// SYMBOL "rosenbrock_factorize_alg"
template<typename T1>
int casadi_rosenbrock_factorize_alg(casadi_rosenbrock_data<T1>* d) {
  // QR of d(alg)/dz; returns 1 if singular
  const casadi_rosenbrock_prob<T1>* p = d->prob;
  T1 rmin = 0;
  casadi_qr(p->sp_jac[ROSENBROCK_ALG_Z], casadi_rosenbrock_jac_nz(d, d->jac, ROSENBROCK_ALG_Z,
    1, 0), d->qr_w, p->sp_alg_v, d->alg_v, p->sp_alg_r, d->alg_r, d->alg_beta, p->alg_prinv,
    p->alg_pc);
  casadi_qr_singular(&rmin, static_cast<casadi_int*>(0), d->alg_r, p->sp_alg_r, p->alg_pc, 0.);
  return !(rmin > 0) || !casadi_rosenbrock_isfinite(rmin);
}

// SYMBOL "rosenbrock_alg_solve"
template<typename T1>
void casadi_rosenbrock_alg_solve(casadi_rosenbrock_data<T1>* d, T1* v, casadi_int nrhs,
    int tr) {
  // Solve with the factorized d(alg)/dz (or its transpose), nrhs right-hand sides
  const casadi_rosenbrock_prob<T1>* p = d->prob;
  casadi_qr_solve(v, nrhs, tr, p->sp_alg_v, d->alg_v, p->sp_alg_r, d->alg_r, d->alg_beta,
    p->alg_prinv, p->alg_pc, d->qr_w);
}

// SYMBOL "rosenbrock_newton_step"
template<typename T1>
void casadi_rosenbrock_newton_step(casadi_rosenbrock_data<T1>* d, T1* v) {
  // Augmented Newton step: g_z w_0 = v_0, g_z w_d = v_d - dg_z,d w_0
  const casadi_rosenbrock_prob<T1>* p = d->prob;
  casadi_int dir;
  T1* vd;
  casadi_rosenbrock_alg_solve(d, v, 1, 0);
  for (dir = 0; dir < p->nfwd; ++dir) {
    vd = v + (dir + 1) * p->nz1;
    casadi_scal(p->nz1, -1., vd);
    casadi_rosenbrock_jac_mv(d, ROSENBROCK_ALG_Z, dir, v, vd, 0);
    casadi_scal(p->nz1, -1., vd);
  }
  if (p->nfwd > 0) casadi_rosenbrock_alg_solve(d, v + p->nz1, p->nfwd, 0);
}

// SYMBOL "rosenbrock_newton_step_trans"
template<typename T1>
void casadi_rosenbrock_newton_step_trans(casadi_rosenbrock_data<T1>* d, T1* v) {
  // Transposed: g_z^T mu = v, g_z^T dmu_d = dv_d - dg_z,d^T mu, for each adjoint direction
  const casadi_rosenbrock_prob<T1>* p = d->prob;
  casadi_int a, dir, nrz1 = p->nz1 * p->nadj;
  T1 *mu, *dmu;
  casadi_rosenbrock_alg_solve(d, v, p->nadj, 1);
  for (dir = 0; dir < p->nfwd; ++dir) {
    for (a = 0; a < p->nadj; ++a) {
      mu = v + a * p->nz1;
      dmu = v + (dir + 1) * nrz1 + a * p->nz1;
      casadi_scal(p->nz1, -1., dmu);
      casadi_rosenbrock_jac_mv(d, ROSENBROCK_ALG_Z, dir, mu, dmu, 1);
      casadi_scal(p->nz1, -1., dmu);
    }
  }
  if (p->nfwd > 0) casadi_rosenbrock_alg_solve(d, v + nrz1, p->nadj * p->nfwd, 1);
}

// SYMBOL "rosenbrock_accum_jac_adj"
template<typename T1>
void casadi_rosenbrock_accum_jac_adj(casadi_rosenbrock_data<T1>* d, const T1* w,
    const T1* u, T1 coef) {
  // jac_adj_a += w_a u^T (x, z columns), coef w_a (time column), and forward derivatives
  const casadi_rosenbrock_prob<T1>* p = d->prob;
  casadi_int a, r, c, k, dir, n_row, n_col;
  const T1 *wa, *dwa, *uc;
  T1 one = 1;
  for (a = 0; a < p->nadj; ++a) {
    for (r = 0; r < 3; ++r) {
      // Rows: ode, alg, quad
      n_row = r == 0 ? p->nx1 : r == 1 ? p->nz1 : p->nq1;
      wa = (r == 0 ? w : r == 1 ? w + p->nrx : w + p->nrx + p->nrz) + a * n_row;
      for (c = 0; c < 3; ++c) {
        // Columns: x, z, t
        k = 3 * c + r;
        if (casadi_rosenbrock_nnz(p->sp_jac[k]) == 0) continue;
        n_col = c == 0 ? p->nx1 : p->nz1;
        uc = c == 0 ? u : u + p->nx;
        casadi_rank1(casadi_rosenbrock_jac_nz(d, d->jac_adj, k, p->nadj, a), p->sp_jac[k],
          c < 2 ? one : coef, wa, c < 2 ? uc : &one);
        for (dir = 0; dir < p->nfwd; ++dir) {
          dwa = wa + (dir + 1) * n_row * p->nadj;
          casadi_rank1(casadi_rosenbrock_jac_nz(d, d->jac_adj_fwd, k, p->nadj * p->nfwd,
            dir * p->nadj + a), p->sp_jac[k], c < 2 ? one : coef, dwa, c < 2 ? uc : &one);
          if (c < 2) {
            casadi_rank1(casadi_rosenbrock_jac_nz(d, d->jac_adj_fwd, k, p->nadj * p->nfwd,
              dir * p->nadj + a), p->sp_jac[k], one, wa, uc + (dir + 1) * n_col);
          }
        }
      }
    }
  }
}

// SYMBOL "rosenbrock_add_vjp"
template<typename T1>
void casadi_rosenbrock_add_vjp(casadi_rosenbrock_data<T1>* d, T1* r) {
  // Add x, z adjoints of EVAL_VJP/EVAL_JAC_ADJ to r, p, u adjoints to the running sums
  const casadi_rosenbrock_prob<T1>* p = d->prob;
  if (r) {
    casadi_axpy(p->nrx, 1., d->vjp_x, r);
    casadi_axpy(p->nrz, 1., d->vjp_z, r + p->nrx);
  }
  casadi_axpy(p->nrp, 1., d->vjp_p, d->adj_p);
  casadi_axpy(p->nru, 1., d->vjp_u, d->adj_u);
}

// SYMBOL "rosenbrock_clear_vjp"
template<typename T1>
void casadi_rosenbrock_clear_vjp(casadi_rosenbrock_data<T1>* d) {
  const casadi_rosenbrock_prob<T1>* p = d->prob;
  casadi_clear(d->vjp_x, p->nrx);
  casadi_clear(d->vjp_z, p->nrz);
  casadi_clear(d->vjp_p, p->nrp);
  casadi_clear(d->vjp_u, p->nru);
}

// ---------------------------------------------------------------------------------------
// Tape: records [U (step records only), y0, u, h, t0, t1, end of the output interval, type]
// ---------------------------------------------------------------------------------------

// SYMBOL "rosenbrock_tape_full"
template<typename T1>
int casadi_rosenbrock_tape_full(casadi_rosenbrock_data<T1>* d, casadi_rosenbrock_record_t record) {
  // No room for a record of this type? Sets the capacity needed
  const casadi_rosenbrock_prob<T1>* p = d->prob;
  casadi_int n = (record == ROSENBROCK_RECORD_STEP ? p->ny * p->n_stages : 0) + p->ny + p->nu + 5;
  d->tape_need = d->tape_pos + n;
  return n > d->tape_cap - d->tape_pos;
}

// SYMBOL "rosenbrock_tape_push"
template<typename T1>
void casadi_rosenbrock_tape_push(casadi_rosenbrock_data<T1>* d) {
  // Append the current record, for which room has been reserved
  const casadi_rosenbrock_prob<T1>* p = d->prob;
  casadi_int n_u = d->record == ROSENBROCK_RECORD_STEP ? p->ny * p->n_stages : 0;
  T1* r = d->tape + d->tape_pos;
  casadi_copy(d->U, n_u, r);
  r += n_u;
  casadi_copy(d->y0, p->ny, r);
  casadi_copy(d->u_step, p->nu, r + p->ny);
  r += p->ny + p->nu;
  r[0] = d->h;
  r[1] = d->t0;
  r[2] = d->t;
  r[3] = d->t_out;
  r[4] = d->record;
  d->tape_pos += n_u + p->ny + p->nu + 5;
}

// SYMBOL "rosenbrock_tape_pop"
template<typename T1>
int casadi_rosenbrock_tape_pop(casadi_rosenbrock_data<T1>* d) {
  // Restore the last record, clear the seeds of its stages; returns 1 if the tape is empty
  const casadi_rosenbrock_prob<T1>* p = d->prob;
  casadi_int n_u;
  T1* r;
  if (d->tape_pos < p->ny + p->nu + 5) return 1;
  r = d->tape + d->tape_pos;
  d->record = static_cast<int>(r[-1]) == ROSENBROCK_RECORD_STEP ? ROSENBROCK_RECORD_STEP
    : ROSENBROCK_RECORD_INIT;
  d->t_record = r[-2];
  d->t1 = r[-3];
  d->t0 = r[-4];
  d->h = r[-5];
  n_u = d->record == ROSENBROCK_RECORD_STEP ? p->ny * p->n_stages : 0;
  r -= p->ny + p->nu + 5;
  casadi_copy(r, p->ny, d->y0);
  casadi_copy(r + p->ny, p->nu, d->u_step);
  casadi_copy(r - n_u, n_u, d->U);
  d->tape_pos -= n_u + p->ny + p->nu + 5;
  casadi_clear(d->rU, p->nr * p->n_stages);
  casadi_clear(d->ry0, p->nr);
  casadi_clear(d->ru_step, p->nru);
  return 0;
}

// ---------------------------------------------------------------------------------------
// Step size control
// ---------------------------------------------------------------------------------------

// SYMBOL "rosenbrock_error_norm"
template<typename T1>
T1 casadi_rosenbrock_error_norm(casadi_rosenbrock_data<T1>* d) {
  // Weighted RMS norm of the error estimate over the error controlled entries
  const casadi_rosenbrock_prob<T1>* p = d->prob;
  casadi_int k, i, n, offset[3], n_err[3];
  T1 r, sum;
  T1* e = d->y_stage;
  // e = sum_j m_err_j U_j
  casadi_clear(e, p->ny);
  casadi_mv_dense(d->U, p->ny, p->n_stages, p->m_err, e, 0);
  offset[0] = 0; offset[1] = p->nx; offset[2] = p->nx + p->nz;
  n_err[0] = p->nx_err; n_err[1] = p->nz_err; n_err[2] = p->nq_err;
  sum = 0;
  n = 0;
  for (k = 0; k < 3; ++k) {
    for (i = offset[k]; i < offset[k] + n_err[k]; ++i) {
      r = e[i] / (p->abstol + p->reltol * fmax(fabs(d->y0[i]), fabs(d->y[i])));
      sum += r * r;
      n++;
    }
  }
  return n > 0 ? sqrt(sum / static_cast<T1>(n)) : 0;
}

// SYMBOL "rosenbrock_newton_norm"
template<typename T1>
T1 casadi_rosenbrock_newton_norm(const casadi_rosenbrock_data<T1>* d, const T1* dz) {
  // Weighted RMS norm of a Newton step for the error controlled algebraic states
  const casadi_rosenbrock_prob<T1>* p = d->prob;
  casadi_int i;
  T1 r, sum = 0;
  if (p->nz_err == 0) return 0;
  for (i = 0; i < p->nz_err; ++i) {
    r = dz[i] / (p->abstol + p->reltol * fabs(d->y[p->nx + i]));
    sum += r * r;
  }
  return sqrt(sum / static_cast<T1>(p->nz_err));
}

// SYMBOL "rosenbrock_initial_step"
template<typename T1>
T1 casadi_rosenbrock_initial_step(const casadi_rosenbrock_data<T1>* d, const T1* f0) {
  // Initial step size from |y0| and |f(y0)| over the error controlled states
  const casadi_rosenbrock_prob<T1>* p = d->prob;
  casadi_int i;
  T1 sc, d0 = 0, d1 = 0, h;
  for (i = 0; i < p->nx_err; ++i) {
    sc = p->abstol + p->reltol * fabs(d->y[i]);
    d0 += (d->y[i] / sc) * (d->y[i] / sc);
    d1 += (f0[i] / sc) * (f0[i] / sc);
  }
  d0 = sqrt(d0 / static_cast<T1>(p->nx_err));
  d1 = sqrt(d1 / static_cast<T1>(p->nx_err));
  h = (d0 < 1e-5 || d1 < 1e-5) ? 1e-6 : 0.01 * d0 / d1;
  h = fmin(h, d->t_stop - d->t);
  return fmax(h, p->h_min);
}

// SYMBOL "rosenbrock_controller"
template<typename T1>
T1 casadi_rosenbrock_controller(casadi_rosenbrock_data<T1>* d, T1 err) {
  // PI controller: factor for the next step size after an accepted step
  const casadi_rosenbrock_prob<T1>* p = d->prob;
  T1 k = static_cast<T1>(p->order_err + 1), fac;
  err = fmax(err, 1e-10);
  fac = p->safety * pow(err, -0.7 / k) * pow(d->err_prev, 0.4 / k);
  fac = fmin(p->fac_max, fmax(p->fac_min, fac));
  if (d->last_rejected) fac = fmin(fac, 1.);
  d->err_prev = fmax(err, 1e-4);
  return fac;
}

// SYMBOL "rosenbrock_reject"
template<typename T1>
int casadi_rosenbrock_reject(casadi_rosenbrock_data<T1>* d, T1 fac) {
  // Restore state, shrink step size; returns 1 if the step size became too small
  const casadi_rosenbrock_prob<T1>* p = d->prob;
  casadi_copy(d->y0, p->ny, d->y);
  d->n_reject++;
  d->last_rejected = 1;
  d->h_next = d->h * fac;
  d->next = ROSENBROCK_STEP_BEGIN;
  d->fail = 0;
  if (p->fixed_steps > 0 || d->nom) {
    // The step sizes are prescribed
    casadi_rosenbrock_stop(d, ROSENBROCK_EVAL_FAILED);
    return 1;
  }
  if (d->h_next < p->h_min) {
    casadi_rosenbrock_stop(d, ROSENBROCK_STEP_TOO_SMALL);
    return 1;
  }
  return 0;
}

// SYMBOL "rosenbrock_stage_argument"
template<typename T1>
void casadi_rosenbrock_stage_argument(casadi_rosenbrock_data<T1>* d, casadi_int i) {
  // y_stage = y0 + sum_{j<i} a_ij U_j
  const casadi_rosenbrock_prob<T1>* p = d->prob;
  casadi_copy(d->y0, p->ny, d->y_stage);
  casadi_mv_dense(d->U, p->ny, i, p->a + i * p->n_stages, d->y_stage, 0);
}

// SYMBOL "rosenbrock_stage_rhs"
template<typename T1>
void casadi_rosenbrock_stage_rhs(casadi_rosenbrock_data<T1>* d, casadi_int i) {
  // U_i += sum_{j<i} (c_ij/h) M U_j + gamma_i h dF/dt, with U_i holding F on entry
  const casadi_rosenbrock_prob<T1>* p = d->prob;
  casadi_int j;
  T1 *Ui = d->U + i * p->ny, cij;
  for (j = 0; j < i; ++j) {
    cij = p->c[i * p->n_stages + j] / d->h;
    // M is the identity on x and q, zero on z
    casadi_axpy(p->nx, cij, d->U + j * p->ny, Ui);
    casadi_axpy(p->nq, cij, d->U + j * p->ny + p->nx + p->nz, Ui + p->nx + p->nz);
  }
  casadi_axpy(p->ny, p->gamma_sum[i] * d->h, d->dfdt, Ui);
}

// SYMBOL "rosenbrock_stage_request"
template<typename T1>
int casadi_rosenbrock_stage_request(casadi_rosenbrock_data<T1>* d) {
  // Right-hand side at the argument of the current stage
  const casadi_rosenbrock_prob<T1>* p = d->prob;
  casadi_rosenbrock_stage_argument(d, d->stage);
  d->t_arg = d->t0 + p->alpha[d->stage] * d->h;
  d->y_arg = d->y_stage;
  d->res = d->U + d->stage * p->ny;
  return casadi_rosenbrock_request(d, ROSENBROCK_EVAL_RHS, ROSENBROCK_STAGE_SOLVE);
}

// SYMBOL "rosenbrock_dense_weights"
template<typename T1>
void casadi_rosenbrock_dense_weights(casadi_rosenbrock_data<T1>* d, T1 theta) {
  // w_i = m_i(theta), the weights of y(t0 + theta h) = y0 + sum_i m_i(theta) U_i
  const casadi_rosenbrock_prob<T1>* p = d->prob;
  casadi_int q;
  T1 th = theta;
  casadi_clear(d->w_dense, p->n_stages);
  for (q = 0; q < p->n_dense; ++q) {
    casadi_axpy(p->n_stages, th, p->m_dense + q * p->n_stages, d->w_dense);
    th *= theta;
  }
}

// SYMBOL "rosenbrock_interpolate"
template<typename T1>
void casadi_rosenbrock_interpolate(casadi_rosenbrock_data<T1>* d, T1 t, T1* y) {
  // Dense output of the current step: y = y0 + sum_i m_i(theta) U_i, weights in d->w_dense
  const casadi_rosenbrock_prob<T1>* p = d->prob;
  casadi_rosenbrock_dense_weights(d, (t - d->t0) / d->h);
  casadi_copy(d->y0, p->ny, y);
  casadi_mv_dense(d->U, p->ny, p->n_stages, d->w_dense, y, 0);
}

// SYMBOL "rosenbrock_seed_stages"
template<typename T1>
void casadi_rosenbrock_seed_stages(casadi_rosenbrock_data<T1>* d) {
  // Adjoint of the dense output with seed d->seed, to the stages and the start of the step
  const casadi_rosenbrock_prob<T1>* p = d->prob;
  casadi_int i;
  for (i = 0; i < p->n_stages; ++i) {
    casadi_axpy(p->nr, d->w_dense[i], d->seed, d->rU + i * p->nr);
  }
  casadi_axpy(p->nr, 1., d->seed, d->ry0);
}

// SYMBOL "rosenbrock_impulse"
template<typename T1>
void casadi_rosenbrock_impulse(casadi_rosenbrock_data<T1>* d, T1 t) {
  // Adjoint seed d->seed of the output at time t; tasks follow, see casadi_rosenbrock
  d->t_seed = t;
  d->status = ROSENBROCK_SUCCESS;
  d->next = ROSENBROCK_IMPULSE_BEGIN;
}

// ---------------------------------------------------------------------------------------
// Forward integration and discrete adjoint
// ---------------------------------------------------------------------------------------

// SYMBOL "rosenbrock"
template<typename T1>
int casadi_rosenbrock(casadi_rosenbrock_data<T1>* d);

// SYMBOL "rosenbrock_forward"
template<typename T1>
int casadi_rosenbrock_forward(casadi_rosenbrock_data<T1>* d) {
  // Forward integration; returns 1 if a task is requested, 0 when done or failed
  const casadi_rosenbrock_prob<T1>* p = d->prob;
  T1 err, fac, *v;
  switch (d->next) {
  // Consistent initialization: Newton's method for alg(t, x, z) = 0, x fixed
  case ROSENBROCK_INIT_BEGIN:
    if (d->fail) return casadi_rosenbrock_stop(d, ROSENBROCK_TAPE_FAILED);
    if (d->taping && casadi_rosenbrock_tape_full(d, ROSENBROCK_RECORD_INIT)) {
      if (p->tape_size > 0) return casadi_rosenbrock_stop(d, ROSENBROCK_TAPE_FAILED);
      return casadi_rosenbrock_request(d, ROSENBROCK_TAPE_RESERVE, ROSENBROCK_INIT_BEGIN);
    }
    casadi_copy(d->u, p->nu, d->u_step);
    d->init_iter = 0;
    d->t_arg = d->t;
    d->y_arg = d->y;
    d->res = d->U;
    return casadi_rosenbrock_request(d, ROSENBROCK_EVAL_RHS, ROSENBROCK_INIT_LINEARIZE);
  case ROSENBROCK_INIT_LINEARIZE:
    if (d->fail) return casadi_rosenbrock_stop(d, ROSENBROCK_INIT_FAILED);
    d->n_rhs++;
    return casadi_rosenbrock_request(d, ROSENBROCK_EVAL_JAC, ROSENBROCK_INIT_NEWTON_STEP);
  case ROSENBROCK_INIT_NEWTON_STEP:
    if (d->fail || casadi_rosenbrock_factorize_alg(d)) {
      return casadi_rosenbrock_stop(d, ROSENBROCK_INIT_FAILED);
    }
    d->n_jac++;
    // z <- z - g_z^{-1} g, with sensitivities
    v = d->U + p->nx;
    casadi_rosenbrock_newton_step(d, v);
    casadi_axpy(p->nz, -1., v, d->y + p->nx);
    err = casadi_rosenbrock_newton_norm(d, v);
    if (++d->init_iter > p->init_max_iter || !casadi_rosenbrock_isfinite(err)) {
      return casadi_rosenbrock_stop(d, ROSENBROCK_INIT_FAILED);
    }
    if (err > p->init_tol) {
      // Next Newton iteration
      d->res = d->U;
      return casadi_rosenbrock_request(d, ROSENBROCK_EVAL_RHS, ROSENBROCK_INIT_LINEARIZE);
    }
    if (p->nfwd == 0) {
      d->next = ROSENBROCK_INIT_END;
      return casadi_rosenbrock_forward(d);
    }
    // Consistent algebraic sensitivities at the solution (implicit function theorem)
    return casadi_rosenbrock_request(d, ROSENBROCK_EVAL_JAC, ROSENBROCK_INIT_SENS_RHS);
  case ROSENBROCK_INIT_SENS_RHS:
    if (d->fail || casadi_rosenbrock_factorize_alg(d)) {
      return casadi_rosenbrock_stop(d, ROSENBROCK_INIT_FAILED);
    }
    d->n_jac++;
    d->res = d->U;
    return casadi_rosenbrock_request(d, ROSENBROCK_EVAL_RHS, ROSENBROCK_INIT_SENS);
  case ROSENBROCK_INIT_SENS:
    if (d->fail) return casadi_rosenbrock_stop(d, ROSENBROCK_INIT_FAILED);
    d->n_rhs++;
    // dz_d -= g_z^{-1} (g_x dx_d + g_z dz_d + g_p dp_d + g_u du_d)
    v = d->U + p->nx + p->nz1;
    casadi_rosenbrock_alg_solve(d, v, p->nfwd, 0);
    casadi_axpy(p->nz1 * p->nfwd, -1., v, d->y + p->nx + p->nz1);
    d->next = ROSENBROCK_INIT_END;
    return casadi_rosenbrock_forward(d);
  case ROSENBROCK_INIT_END:
    d->init_pending = 0;
    if (d->taping) {
      d->record = ROSENBROCK_RECORD_INIT;
      d->t0 = d->t;
      d->h = 0;
      casadi_copy(d->y, p->ny, d->y0);
      casadi_rosenbrock_tape_push(d);
    }
    d->next = ROSENBROCK_STEP_BEGIN;
    return casadi_rosenbrock_forward(d);
  // Integrator step
  case ROSENBROCK_STEP_BEGIN:
    if (d->fail) return casadi_rosenbrock_stop(d, ROSENBROCK_TAPE_FAILED);
    if (d->t >= d->t_out) {
      if (d->pause) return casadi_rosenbrock_stop(d, ROSENBROCK_SUCCESS);
      if (d->t == d->t_out) {
        casadi_copy(d->y, p->ny, d->y_out);
        return casadi_rosenbrock_stop(d, ROSENBROCK_SUCCESS);
      }
      // Output inside the last step
      casadi_rosenbrock_interpolate(d, d->t_out, d->y_out);
      if (p->nz == 0) return casadi_rosenbrock_stop(d, ROSENBROCK_SUCCESS);
      // Algebraic states: Newton step on alg(t_out, x, z) = 0, d(alg)/dz at the step start
      if (casadi_rosenbrock_factorize_alg(d)) {
        return casadi_rosenbrock_stop(d, ROSENBROCK_EVAL_FAILED);
      }
      d->t_arg = d->t_out;
      d->y_arg = d->y_out;
      d->res = d->y_stage;
      return casadi_rosenbrock_request(d, ROSENBROCK_EVAL_RHS, ROSENBROCK_OUTPUT_CORRECT);
    }
    if (d->n_steps_interval >= p->max_num_steps) {
      return casadi_rosenbrock_stop(d, ROSENBROCK_MAX_NUM_STEPS);
    }
    if (d->taping && casadi_rosenbrock_tape_full(d, ROSENBROCK_RECORD_STEP)) {
      if (p->tape_size > 0) return casadi_rosenbrock_stop(d, ROSENBROCK_TAPE_FAILED);
      return casadi_rosenbrock_request(d, ROSENBROCK_TAPE_RESERVE, ROSENBROCK_STEP_BEGIN);
    }
    // Controls of the new step
    casadi_copy(d->u, p->nu, d->u_step);
    d->last_clipped = 0;
    d->last_step = 0;
    if (d->nom && p->fixed_steps == 0) {
      // Step of the nominal integration
      d->next = ROSENBROCK_STEP_REPLAY;
      return casadi_rosenbrock_forward(d);
    } else if (p->fixed_steps > 0) {
      // Equidistant steps
      d->h = (d->t_out - d->t) / static_cast<T1>(p->fixed_steps - d->n_steps_interval);
      d->last_step = d->n_steps_interval + 1 == p->fixed_steps;
    } else {
      if (d->h_next <= 0) {
        if (p->h_init > 0) {
          d->h_next = p->h_init;
        } else {
          // Estimate from the right-hand side at the initial point
          d->t_arg = d->t;
          d->y_arg = d->y;
          d->res = d->U;
          return casadi_rosenbrock_request(d, ROSENBROCK_EVAL_RHS, ROSENBROCK_STEP_INITIAL_SIZE);
        }
      }
      d->h = p->h_max > 0 ? fmin(d->h_next, p->h_max) : d->h_next;
      if (d->t + 1.05 * d->h >= d->t_stop) {
        // Hit the stop time exactly
        d->last_clipped = d->t + d->h > d->t_stop;
        d->h = d->t_stop - d->t;
        d->last_step = 1;
      } else if (d->t + 2 * d->h > d->t_stop) {
        // Avoid a tiny last step
        d->h = 0.5 * (d->t_stop - d->t);
      }
    }
    // Linearize at the beginning of the step
    d->t0 = d->t;
    casadi_copy(d->y, p->ny, d->y0);
    d->t_arg = d->t0;
    d->y_arg = d->y0;
    return casadi_rosenbrock_request(d, ROSENBROCK_EVAL_JAC, ROSENBROCK_STEP_FACTORIZE);
  case ROSENBROCK_STEP_REPLAY:
    // Advance the nominal integration by one step; its tasks are handed to the caller
    if (casadi_rosenbrock(d->nom)) {
      d->active = d->nom;
      return 1;
    }
    d->active = d;
    if (d->nom->status != ROSENBROCK_SUCCESS) return casadi_rosenbrock_stop(d, d->nom->status);
    if (d->nom->next == ROSENBROCK_DONE) return casadi_rosenbrock_stop(d, ROSENBROCK_EVAL_FAILED);
    // Same step size, same end time
    d->h = d->nom->h;
    d->t0 = d->t;
    casadi_copy(d->y, p->ny, d->y0);
    d->t_arg = d->t0;
    d->y_arg = d->y0;
    return casadi_rosenbrock_request(d, ROSENBROCK_EVAL_JAC, ROSENBROCK_STEP_FACTORIZE);
  case ROSENBROCK_STEP_INITIAL_SIZE:
    if (d->fail) return casadi_rosenbrock_stop(d, ROSENBROCK_EVAL_FAILED);
    d->n_rhs++;
    d->h_next = casadi_rosenbrock_initial_step(d, d->U);
    d->next = ROSENBROCK_STEP_BEGIN;
    return casadi_rosenbrock_forward(d);
  case ROSENBROCK_STEP_FACTORIZE:
    if (d->fail) return casadi_rosenbrock_stop(d, ROSENBROCK_EVAL_FAILED);
    d->n_jac++;
    casadi_rosenbrock_dfdt(d);
    d->h_gamma = d->h * p->gamma;
    if (casadi_rosenbrock_factorize(d)) {
      // Singular iteration matrix: retry with a smaller step
      if (casadi_rosenbrock_reject(d, 0.25)) return 0;
      return casadi_rosenbrock_forward(d);
    }
    d->stage = 0;
    return casadi_rosenbrock_stage_request(d);
  case ROSENBROCK_STAGE_SOLVE:
    if (d->fail) {
      // Evaluation failed: retry with a smaller step
      if (casadi_rosenbrock_reject(d, 0.25)) return 0;
      return casadi_rosenbrock_forward(d);
    }
    d->n_rhs++;
    casadi_rosenbrock_stage_rhs(d, d->stage);
    casadi_rosenbrock_solve(d, d->U + d->stage * p->ny);
    if (++d->stage < p->n_stages) return casadi_rosenbrock_stage_request(d);
    // Candidate solution y = y0 + sum_i m_i U_i
    casadi_copy(d->y0, p->ny, d->y);
    casadi_mv_dense(d->U, p->ny, p->n_stages, p->m, d->y, 0);
    if (p->fixed_steps == 0 && !d->nom) {
      err = casadi_rosenbrock_error_norm(d);
      if (!(err <= 1)) {
        // Reject (also if not a number)
        fac = 0.25;
        if (casadi_rosenbrock_isfinite(err)) {
          fac = fmax(p->fac_min, p->safety * pow(err, -1. / (p->order_err + 1)));
        }
        if (casadi_rosenbrock_reject(d, fac)) return 0;
        return casadi_rosenbrock_forward(d);
      }
      fac = casadi_rosenbrock_controller(d, err);
      // A step clipped at the stop time does not shrink the proposal
      d->h_next = d->last_clipped ? fmax(d->h_next, d->h * fac) : d->h * fac;
    }
    d->last_rejected = 0;
    if (d->nom && p->fixed_steps == 0) {
      d->t = d->nom->t;
    } else {
      d->t = d->last_step ? d->t_stop : d->t0 + d->h;
    }
    if (d->taping) {
      d->record = ROSENBROCK_RECORD_STEP;
      casadi_rosenbrock_tape_push(d);
    }
    d->n_steps++;
    d->n_steps_interval++;
    d->next = ROSENBROCK_STEP_BEGIN;
    if (d->pause) {
      // Hand the step to the integration that replays it
      d->status = ROSENBROCK_SUCCESS;
      return 0;
    }
    return casadi_rosenbrock_forward(d);
  case ROSENBROCK_OUTPUT_CORRECT:
    if (d->fail) return casadi_rosenbrock_stop(d, ROSENBROCK_EVAL_FAILED);
    d->n_rhs++;
    v = d->y_stage + p->nx;
    casadi_rosenbrock_newton_step(d, v);
    casadi_axpy(p->nz, -1., v, d->y_out + p->nx);
    return casadi_rosenbrock_stop(d, ROSENBROCK_SUCCESS);
  default:
    return 0;
  }
}

// SYMBOL "rosenbrock_backward"
template<typename T1>
int casadi_rosenbrock_backward(casadi_rosenbrock_data<T1>* d) {
  // Discrete adjoint of the taped steps; returns 1 if a task is requested, 0 when done
  const casadi_rosenbrock_prob<T1>* p = d->prob;
  casadi_int i, j, s = p->n_stages;
  T1 cij, *rUi, *rUj;
  switch (d->next) {
  case ROSENBROCK_BACK_BEGIN:
    // Next record, unless one from an earlier output interval is pending
    if (!d->record_pending && casadi_rosenbrock_tape_pop(d)) {
      // Tape exhausted
      return casadi_rosenbrock_stop(d, ROSENBROCK_SUCCESS);
    }
    // Record belongs to a later output interval, skipped with zero adjoint?
    if (d->t_record > d->t_start) {
      d->next = ROSENBROCK_BACK_BEGIN;
      return casadi_rosenbrock_backward(d);
    }
    // Record belongs to an earlier output interval?
    if (d->t_record < d->t_start) {
      d->record_pending = 1;
      d->t = d->t_out;
      return casadi_rosenbrock_stop(d, ROSENBROCK_SUCCESS);
    }
    d->record_pending = 0;
    // Linearize at the beginning of the record
    d->t_arg = d->t0;
    d->y_arg = d->y0;
    return casadi_rosenbrock_request(d, ROSENBROCK_EVAL_JAC,
      d->record == ROSENBROCK_RECORD_INIT ? ROSENBROCK_BACK_INIT_SOLVE
      : ROSENBROCK_BACK_ADJ_BEGIN);
  case ROSENBROCK_BACK_INIT_SOLVE:
    if (d->fail || casadi_rosenbrock_factorize_alg(d)) break;
    // Implicit function theorem: mu = (d(alg)/dz)^{-T} zbar, xbar -= (d(alg)/dx)^T mu
    casadi_clear(d->v1, p->nr);
    casadi_copy(d->ry + p->nrx, p->nrz, d->v1 + p->nrx);
    casadi_rosenbrock_newton_step_trans(d, d->v1 + p->nrx);
    casadi_scal(p->nrz, -1., d->v1 + p->nrx);
    casadi_rosenbrock_clear_vjp(d);
    d->w = d->v1;
    return casadi_rosenbrock_request(d, ROSENBROCK_EVAL_VJP, ROSENBROCK_BACK_INIT_END);
  case ROSENBROCK_BACK_INIT_END:
    if (d->fail) break;
    d->n_vjp++;
    casadi_rosenbrock_add_vjp(d, d->ry);
    // The algebraic adjoint has been propagated into the differential states
    casadi_clear(d->ry + p->nrx, p->nrz);
    d->t = d->t0;
    d->next = ROSENBROCK_BACK_BEGIN;
    return casadi_rosenbrock_backward(d);
  case ROSENBROCK_BACK_ADJ_BEGIN:
    if (d->fail) break;
    d->h_gamma = d->h * p->gamma;
    if (casadi_rosenbrock_factorize(d)) break;
    // Stage adjoints from y1 = y0 + sum_i m_i U_i, added to the seeds of dense outputs
    for (i = 0; i < s; ++i) casadi_axpy(p->nr, p->m[i], d->ry, d->rU + i * p->nr);
    casadi_axpy(p->nr, 1., d->ry0, d->ry);
    casadi_axpy(p->nru, 1., d->ru_step, d->adj_u);
    d->stage = s;
    d->next = ROSENBROCK_BACK_ADJ_STAGE;
    return casadi_rosenbrock_backward(d);
  case ROSENBROCK_BACK_ADJ_SPREAD:
    // Adjoint of y0 + sum_j a_ij U_j and sum_j (c_ij/h) M U_j to y0 and earlier stages
    if (d->fail) break;
    d->n_vjp++;
    i = d->stage;
    rUi = d->rU + i * p->nr;
    casadi_rosenbrock_add_vjp(d, d->ry);
    for (j = 0; j < i; ++j) {
      rUj = d->rU + j * p->nr;
      casadi_axpy(p->nrx, p->a[i * s + j], d->vjp_x, rUj);
      casadi_axpy(p->nrz, p->a[i * s + j], d->vjp_z, rUj + p->nrx);
      cij = p->c[i * s + j] / d->h;
      casadi_axpy(p->nrx, cij, rUi, rUj);
      casadi_axpy(p->nrq, cij, rUi + p->nrx + p->nrz, rUj + p->nrx + p->nrz);
    }
    d->next = ROSENBROCK_BACK_ADJ_STAGE;
    return casadi_rosenbrock_backward(d);
  case ROSENBROCK_BACK_ADJ_STAGE:
    if (d->stage == 0) {
      // Contribution through J(t0, y0), dF/dt(t0, y0)
      casadi_rosenbrock_clear_vjp(d);
      d->t_arg = d->t0;
      d->y_arg = d->y0;
      return casadi_rosenbrock_request(d, ROSENBROCK_EVAL_JAC_ADJ, ROSENBROCK_BACK_STEP_END);
    }
    // Previous stage: rbar_i = A^{-T} Ubar_i, then through J, dF/dt at y0 and F at the stage
    i = --d->stage;
    rUi = d->rU + i * p->nr;
    casadi_rosenbrock_solve_trans(d, rUi);
    casadi_rosenbrock_accum_jac_adj(d, rUi, d->U + i * p->ny, p->gamma_sum[i] * d->h);
    casadi_rosenbrock_stage_argument(d, i);
    casadi_rosenbrock_clear_vjp(d);
    d->t_arg = d->t0 + p->alpha[i] * d->h;
    d->y_arg = d->y_stage;
    d->w = rUi;
    return casadi_rosenbrock_request(d, ROSENBROCK_EVAL_VJP, ROSENBROCK_BACK_ADJ_SPREAD);
  case ROSENBROCK_BACK_STEP_END:
    if (d->fail) break;
    casadi_rosenbrock_add_vjp(d, d->ry);
    casadi_clear(d->jac_adj, p->nnz_jac * p->nadj * (1 + p->nfwd));
    d->t = d->t0;
    d->next = ROSENBROCK_BACK_BEGIN;
    return casadi_rosenbrock_backward(d);
  // Adjoint seed of the output at t_seed
  case ROSENBROCK_IMPULSE_BEGIN:
    // The record containing t_seed; records of later output intervals had zero seeds
    while (!d->record_pending && !casadi_rosenbrock_tape_pop(d)) {
      d->record_pending = d->t_record <= d->t_seed;
    }
    if (!d->record_pending || d->record != ROSENBROCK_RECORD_STEP || d->t_seed == d->t1) {
      // Seed of the state at the end of the pending record
      casadi_axpy(p->nr, 1., d->seed, d->ry);
      return casadi_rosenbrock_stop(d, ROSENBROCK_SUCCESS);
    }
    // Output inside the pending step
    casadi_rosenbrock_interpolate(d, d->t_seed, d->y_out);
    if (p->nz == 0) {
      casadi_rosenbrock_seed_stages(d);
      return casadi_rosenbrock_stop(d, ROSENBROCK_SUCCESS);
    }
    // Adjoint of the Newton step on the algebraic states, as in the forward integration
    d->t_arg = d->t0;
    d->y_arg = d->y0;
    return casadi_rosenbrock_request(d, ROSENBROCK_EVAL_JAC, ROSENBROCK_IMPULSE_LINEARIZE);
  case ROSENBROCK_IMPULSE_LINEARIZE:
    if (d->fail || casadi_rosenbrock_factorize_alg(d)) break;
    d->t_arg = d->t_seed;
    d->y_arg = d->y_out;
    d->res = d->y_stage;
    return casadi_rosenbrock_request(d, ROSENBROCK_EVAL_RHS, ROSENBROCK_IMPULSE_CORRECT);
  case ROSENBROCK_IMPULSE_CORRECT:
    // z_out = z - G^{-1} alg with G = d(alg)/dz at the step start: dz = G^{-1} alg,
    // mu = G^{-T} zbar, Gbar = mu dz^T, (x, z, p, u)bar -= d(alg)/d(x, z, p, u)^T mu
    if (d->fail) break;
    casadi_rosenbrock_newton_step(d, d->y_stage + p->nx);
    casadi_clear(d->v1, p->nr);
    casadi_copy(d->seed + p->nrx, p->nrz, d->v1 + p->nrx);
    casadi_rosenbrock_newton_step_trans(d, d->v1 + p->nrx);
    casadi_clear(d->v2, p->ny);
    casadi_copy(d->y_stage + p->nx, p->nz, d->v2 + p->nx);
    casadi_rosenbrock_accum_jac_adj(d, d->v1, d->v2, 0.);
    casadi_scal(p->nrz, -1., d->v1 + p->nrx);
    casadi_rosenbrock_clear_vjp(d);
    d->w = d->v1;
    return casadi_rosenbrock_request(d, ROSENBROCK_EVAL_VJP, ROSENBROCK_IMPULSE_END);
  case ROSENBROCK_IMPULSE_END:
    if (d->fail) break;
    d->n_vjp++;
    casadi_axpy(p->nrx, 1., d->vjp_x, d->seed);
    casadi_axpy(p->nrz, 1., d->vjp_z, d->seed + p->nrx);
    casadi_axpy(p->nrp, 1., d->vjp_p, d->adj_p);
    // Control adjoints belong to the output interval of the step
    casadi_axpy(p->nru, 1., d->vjp_u, d->ru_step);
    casadi_rosenbrock_seed_stages(d);
    return casadi_rosenbrock_stop(d, ROSENBROCK_SUCCESS);
  default:
    return 0;
  }
  return casadi_rosenbrock_stop(d, ROSENBROCK_EVAL_FAILED);
}

template<typename T1>
int casadi_rosenbrock(casadi_rosenbrock_data<T1>* d) {
  // Reverse communication: returns 1 while a task in d->task is to be carried out
  if (d->next >= ROSENBROCK_BACK_BEGIN) return casadi_rosenbrock_backward(d);
  return casadi_rosenbrock_forward(d);
}

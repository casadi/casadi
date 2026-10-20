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

// Shared part of the one-step integrators for M y' = F(t, y), y = [x; z; q], M = diag(I, 0, I):
// stages, dense output, tape, output times, stops and the record-wise discrete adjoint
// Layouts: x = [x_0; ...; x_nfwd] (nominal, forward), rx = [rx_0 (nx1 x nadj); ...; rx_nfwd]

// C-REPLACE "static_cast<T1>" "(casadi_real) "
// C-REPLACE "static_cast<casadi_int>" "(casadi_int) "

// SYMBOL "ivp_task_t"
typedef enum {
  // res = F(t_arg, y_arg), with forward sensitivities
  IVP_EVAL_RHS,
  // vjp = dF/d(x,z,p,u)^T w at (t_arg, y_arg), with forward sensitivities
  IVP_EVAL_VJP,
  // Enlarge the tape to at least tape_need entries, keeping its contents; set tape, tape_cap
  IVP_TAPE_RESERVE
} casadi_ivp_task_t;

// SYMBOL "ivp_entry_t"
typedef enum {
  // Entry points of a method; its own phases start at IVP_NUM_ENTRY
  IVP_DONE,
  IVP_ADVANCE,
  IVP_IMPULSE,
  IVP_RETREAT,
  IVP_NUM_ENTRY
} casadi_ivp_entry_t;

// SYMBOL "ivp_flag_t"
typedef enum {
  IVP_SUCCESS,
  IVP_MAX_NUM_STEPS,
  IVP_EVAL_FAILED,
  IVP_TAPE_FAILED
} casadi_ivp_flag_t;

// SYMBOL "ivp_prob"
template<typename T1>
struct casadi_ivp_prob {
  // Nominal block sizes: differential, algebraic, quadrature states, parameters, controls
  casadi_int nx1, nz1, nq1, np1, nu1;
  // Number of forward and adjoint directions
  casadi_int nfwd, nadj;
  // Stages: increments U_i, arguments y0 + sum_j a_ij U_j at t0 + alpha_i h (a strictly lower
  // triangular, row-major), solution y1 = y0 + sum_i m_i U_i
  casadi_int n_stages;
  const T1 *a, *alpha, *m;
  // Dense output weights m_i(theta) = sum_q m_dense[q * n_stages + i] theta^(q+1)
  casadi_int n_dense;
  const T1* m_dense;
  // Capacity of the tape in the work vector (0: the caller provides the tape)
  casadi_int tape_size;
  // Fixed steps: each run of steps up to a stop has equal steps of at most step_size
  T1 step_size;
  // Derived: forward and adjoint state dimensions, controls, vector lengths
  casadi_int nx, nz, nq, nu, nrx, nrz, nrq, nrp, nru, ny, nr;
};
// C-REPLACE "casadi_ivp_prob<T1>" "struct casadi_ivp_prob"

// SYMBOL "ivp_setup"
template<typename T1>
void casadi_ivp_setup(casadi_ivp_prob<T1>* p) {
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
}

// C-REPLACE "casadi_ivp_data<T1>" "struct casadi_ivp_data"
// SYMBOL "ivp_data"
template<typename T1>
struct casadi_ivp_data {
  // Problem structure
  const casadi_ivp_prob<T1>* prob;
  // Task for the caller, where to resume (casadi_ivp_entry_t or a phase of the method), status
  casadi_ivp_task_t task;
  casadi_int next;
  casadi_ivp_flag_t status;
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
  // Current step: start time, size, start state, stage counter
  T1 t0, h;
  // End time and end of the output interval of a popped record
  T1 t1, t_record;
  T1* y0;
  casadi_int stage;
  // Stage increments U_i (n_stages x ny), stage argument, dense output weights
  T1 *U, *y_stage, *w_dense;
  // Adjoint state, stage adjoints (n_stages x nr), adjoint of the start of the popped step
  T1 *ry, *rU, *ry0;
  // Adjoint seed of an output at time t_seed
  T1* seed;
  T1 t_seed;
  // Control adjoints of the popped step from outputs inside it
  T1* ru_step;
  // Results of EVAL_VJP: adjoints of x, z, p, u
  T1 *vjp_x, *vjp_z, *vjp_p, *vjp_u;
  // Accumulated adjoints of the parameters and controls
  T1 *adj_p, *adj_u;
  // Tape: records [U, y0, u, h, t0, t1, t_out]; entries used, capacity, needed
  T1* tape;
  casadi_int tape_pos, tape_cap, tape_need;
  // Reset pending, taping requested
  int reset, taping;
  // Popped record that belongs to an earlier interval
  int record_pending;
  // Fixed steps: current run start, stop, step size, number of steps, steps taken
  T1 run_t0, run_t1, run_h;
  casadi_int run_n, run_k;
  // Statistics
  casadi_int n_steps, n_rhs, n_vjp;
};

// SYMBOL "ivp_work"
template<typename T1>
void casadi_ivp_work(const casadi_ivp_prob<T1>* p, casadi_int* sz_arg,
    casadi_int* sz_res, casadi_int* sz_iw, casadi_int* sz_w) {
  // State: y, y_out, y0, stages, stage argument, dense output weights
  *sz_w += p->ny * (p->n_stages + 4) + p->n_stages;
  // Adjoint: adjoint state, stage adjoints, adjoint of the step start, seed
  *sz_w += p->nr * (p->n_stages + 3);
  // Controls of a step, their adjoints from outputs inside it
  *sz_w += p->nu + p->nru;
  // Function adjoints, accumulated adjoints
  *sz_w += p->nrx + p->nrz + 2 * (p->nrp + p->nru);
  // Tape
  *sz_w += p->tape_size;
}

// SYMBOL "ivp_set_work"
template<typename T1>
void casadi_ivp_set_work(casadi_ivp_data<T1>* d, const T1*** arg, T1*** res,
    casadi_int** iw, T1** w) {
  const casadi_ivp_prob<T1>* p = d->prob;
  d->y = *w; *w += p->ny;
  d->y_out = *w; *w += p->ny;
  d->y0 = *w; *w += p->ny;
  d->U = *w; *w += p->ny * p->n_stages;
  d->y_stage = *w; *w += p->ny;
  d->w_dense = *w; *w += p->n_stages;
  d->ry = *w; *w += p->nr;
  d->rU = *w; *w += p->nr * p->n_stages;
  d->ry0 = *w; *w += p->nr;
  d->seed = *w; *w += p->nr;
  d->u_step = *w; *w += p->nu;
  d->ru_step = *w; *w += p->nru;
  d->vjp_x = *w; *w += p->nrx;
  d->vjp_z = *w; *w += p->nrz;
  d->vjp_p = *w; *w += p->nrp;
  d->vjp_u = *w; *w += p->nru;
  d->adj_p = *w; *w += p->nrp;
  d->adj_u = *w; *w += p->nru;
  d->tape = *w; *w += p->tape_size;
  d->tape_cap = p->tape_size;
  d->tape_pos = 0;
  d->next = IVP_DONE;
  d->status = IVP_SUCCESS;
  d->fail = 0;
  d->reset = 0;
  d->taping = 0;
  d->record_pending = 0;
}

// SYMBOL "ivp_clear_stats"
template<typename T1>
void casadi_ivp_clear_stats(casadi_ivp_data<T1>* d) {
  d->n_steps = d->n_rhs = d->n_vjp = 0;
}

// SYMBOL "ivp_reset"
template<typename T1>
void casadi_ivp_reset(casadi_ivp_data<T1>* d, T1 t, const T1* u, int first_call) {
  // Start at time t, controls u; the state d->y is set by the caller
  d->t = t;
  d->u = u;
  d->reset = 1;
  if (first_call) {
    d->tape_pos = 0;
    casadi_ivp_clear_stats(d);
  }
}

// SYMBOL "ivp_steps"
template<typename T1>
casadi_int casadi_ivp_steps(T1 len, T1 step_size) {
  // Number of equal steps covering len, of at most step_size up to rounding
  casadi_int n = static_cast<casadi_int>(ceil(len / step_size - 1e-9));
  return n < 1 ? 1 : n;
}

// SYMBOL "ivp_begin_run"
template<typename T1>
void casadi_ivp_begin_run(casadi_ivp_data<T1>* d) {
  // Handle a pending reset: a new run of steps up to the stop time
  const casadi_ivp_prob<T1>* p = d->prob;
  d->run_t0 = d->t;
  d->run_t1 = d->t_stop;
  d->run_n = casadi_ivp_steps(d->t_stop - d->t, p->step_size);
  d->run_h = (d->run_t1 - d->run_t0) / static_cast<T1>(d->run_n);
  d->run_k = 0;
  d->reset = 0;
}

// SYMBOL "ivp_advance"
template<typename T1>
void casadi_ivp_advance(casadi_ivp_data<T1>* d, T1 t_out, T1 t_stop, int taping) {
  // Output at t_out, stepping no further than t_stop; tasks follow, see the method
  d->t_out = t_out;
  d->t_stop = t_stop;
  d->taping = taping;
  d->status = IVP_SUCCESS;
  d->fail = 0;
  d->next = IVP_ADVANCE;
}

// SYMBOL "ivp_reset_backward"
template<typename T1>
void casadi_ivp_reset_backward(casadi_ivp_data<T1>* d) {
  const casadi_ivp_prob<T1>* p = d->prob;
  casadi_clear(d->ry, p->nr);
  casadi_clear(d->adj_p, p->nrp);
  casadi_clear(d->adj_u, p->nru);
  d->record_pending = 0;
}

// SYMBOL "ivp_retreat"
template<typename T1>
void casadi_ivp_retreat(casadi_ivp_data<T1>* d, T1 t_start, T1 t_out) {
  // Adjoint from t_start back to t_out; skips records of output intervals ending after t_start
  d->t = d->t_start = t_start;
  d->t_out = t_out;
  d->status = IVP_SUCCESS;
  d->next = IVP_RETREAT;
}

// SYMBOL "ivp_impulse"
template<typename T1>
void casadi_ivp_impulse(casadi_ivp_data<T1>* d, T1 t) {
  // Adjoint seed d->seed of the output at time t; tasks follow, see the method
  d->t_seed = t;
  d->status = IVP_SUCCESS;
  d->next = IVP_IMPULSE;
}

// SYMBOL "ivp_request"
template<typename T1>
int casadi_ivp_request(casadi_ivp_data<T1>* d, casadi_ivp_task_t task, casadi_int next) {
  // Hand a task to the caller, resume at phase next
  d->task = task;
  d->next = next;
  d->fail = 0;
  return 1;
}

// SYMBOL "ivp_stop"
template<typename T1>
int casadi_ivp_stop(casadi_ivp_data<T1>* d, casadi_ivp_flag_t status) {
  d->status = status;
  d->next = IVP_DONE;
  return 0;
}

// SYMBOL "ivp_add_vjp"
template<typename T1>
void casadi_ivp_add_vjp(casadi_ivp_data<T1>* d, T1* r) {
  // Add x, z adjoints of EVAL_VJP to r, p, u adjoints to the running sums
  const casadi_ivp_prob<T1>* p = d->prob;
  if (r) {
    casadi_axpy(p->nrx, 1., d->vjp_x, r);
    casadi_axpy(p->nrz, 1., d->vjp_z, r + p->nrx);
  }
  casadi_axpy(p->nrp, 1., d->vjp_p, d->adj_p);
  casadi_axpy(p->nru, 1., d->vjp_u, d->adj_u);
}

// SYMBOL "ivp_clear_vjp"
template<typename T1>
void casadi_ivp_clear_vjp(casadi_ivp_data<T1>* d) {
  const casadi_ivp_prob<T1>* p = d->prob;
  casadi_clear(d->vjp_x, p->nrx);
  casadi_clear(d->vjp_z, p->nrz);
  casadi_clear(d->vjp_p, p->nrp);
  casadi_clear(d->vjp_u, p->nru);
}

// ---------------------------------------------------------------------------------------
// Tape: records [U, y0, u, h, t0, t1, end of the output interval]
// ---------------------------------------------------------------------------------------

// SYMBOL "ivp_tape_full"
template<typename T1>
int casadi_ivp_tape_full(casadi_ivp_data<T1>* d) {
  // No room for a record? Sets the capacity needed
  const casadi_ivp_prob<T1>* p = d->prob;
  casadi_int n = p->ny * (p->n_stages + 1) + p->nu + 4;
  d->tape_need = d->tape_pos + n;
  return n > d->tape_cap - d->tape_pos;
}

// SYMBOL "ivp_reserve"
template<typename T1>
int casadi_ivp_reserve(casadi_ivp_data<T1>* d, casadi_int next) {
  // Room on the tape for a record, if taping? -1 if so, else the value for the method to return
  if (!d->taping || !casadi_ivp_tape_full(d)) return -1;
  if (d->prob->tape_size > 0) return casadi_ivp_stop(d, IVP_TAPE_FAILED);
  return casadi_ivp_request(d, IVP_TAPE_RESERVE, next);
}

// SYMBOL "ivp_tape_push"
template<typename T1>
void casadi_ivp_tape_push(casadi_ivp_data<T1>* d) {
  // Append the current record, for which room has been reserved
  const casadi_ivp_prob<T1>* p = d->prob;
  casadi_int n_u = p->ny * p->n_stages;
  T1* rec = d->tape + d->tape_pos;
  casadi_copy(d->U, n_u, rec);
  rec += n_u;
  casadi_copy(d->y0, p->ny, rec);
  casadi_copy(d->u_step, p->nu, rec + p->ny);
  rec += p->ny + p->nu;
  rec[0] = d->h;
  rec[1] = d->t0;
  rec[2] = d->t;
  rec[3] = d->t_out;
  d->tape_pos += n_u + p->ny + p->nu + 4;
}

// SYMBOL "ivp_tape_pop"
template<typename T1>
int casadi_ivp_tape_pop(casadi_ivp_data<T1>* d) {
  // Restore the last record, clear the seeds of its stages; returns 1 if the tape is empty
  const casadi_ivp_prob<T1>* p = d->prob;
  casadi_int n_u = p->ny * p->n_stages;
  T1* rec;
  if (d->tape_pos < n_u + p->ny + p->nu + 4) return 1;
  rec = d->tape + d->tape_pos;
  d->t_record = rec[-1];
  d->t1 = rec[-2];
  d->t0 = rec[-3];
  d->h = rec[-4];
  rec -= p->ny + p->nu + 4;
  casadi_copy(rec, p->ny, d->y0);
  casadi_copy(rec + p->ny, p->nu, d->u_step);
  casadi_copy(rec - n_u, n_u, d->U);
  d->tape_pos -= n_u + p->ny + p->nu + 4;
  casadi_clear(d->rU, p->nr * p->n_stages);
  casadi_clear(d->ry0, p->nr);
  casadi_clear(d->ru_step, p->nru);
  return 0;
}

// ---------------------------------------------------------------------------------------
// Stages and dense output
// ---------------------------------------------------------------------------------------

// SYMBOL "ivp_stage_argument"
template<typename T1>
void casadi_ivp_stage_argument(casadi_ivp_data<T1>* d, casadi_int i) {
  // y_stage = y0 + sum_{j<i} a_ij U_j, skipping zero coefficients
  const casadi_ivp_prob<T1>* p = d->prob;
  casadi_int j;
  const T1* a = p->a + i * p->n_stages;
  casadi_copy(d->y0, p->ny, d->y_stage);
  for (j = 0; j < i; ++j) if (a[j] != 0) casadi_axpy(p->ny, a[j], d->U + j * p->ny, d->y_stage);
}

// SYMBOL "ivp_dense_weights"
template<typename T1>
void casadi_ivp_dense_weights(casadi_ivp_data<T1>* d, T1 theta) {
  // w_i = m_i(theta), the weights of y(t0 + theta h) = y0 + sum_i m_i(theta) U_i
  const casadi_ivp_prob<T1>* p = d->prob;
  casadi_int q;
  T1 th = theta;
  casadi_clear(d->w_dense, p->n_stages);
  for (q = 0; q < p->n_dense; ++q) {
    casadi_axpy(p->n_stages, th, p->m_dense + q * p->n_stages, d->w_dense);
    th *= theta;
  }
}

// SYMBOL "ivp_interpolate"
template<typename T1>
void casadi_ivp_interpolate(casadi_ivp_data<T1>* d, T1 t, T1* y) {
  // Dense output of the current step: y = y0 + sum_i m_i(theta) U_i, weights in d->w_dense
  const casadi_ivp_prob<T1>* p = d->prob;
  casadi_ivp_dense_weights(d, (t - d->t0) / d->h);
  casadi_copy(d->y0, p->ny, y);
  casadi_mv_dense(d->U, p->ny, p->n_stages, d->w_dense, y, 0);
}

// SYMBOL "ivp_seed_stages"
template<typename T1>
void casadi_ivp_seed_stages(casadi_ivp_data<T1>* d) {
  // Adjoint of the dense output with seed d->seed, to the stages and the start of the step
  const casadi_ivp_prob<T1>* p = d->prob;
  casadi_int i;
  for (i = 0; i < p->n_stages; ++i) {
    casadi_axpy(p->nr, d->w_dense[i], d->seed, d->rU + i * p->nr);
  }
  casadi_axpy(p->nr, 1., d->seed, d->ry0);
}

// ---------------------------------------------------------------------------------------
// Step sizes
// ---------------------------------------------------------------------------------------

// SYMBOL "ivp_step_limit"
template<typename T1>
int casadi_ivp_step_limit(const casadi_ivp_data<T1>* d) {
  // No further step allowed?
  return d->run_k >= d->run_n;
}

// ---------------------------------------------------------------------------------------
// Building blocks of the forward integration and of the discrete adjoint
// ---------------------------------------------------------------------------------------

// SYMBOL "ivp_output"
template<typename T1>
void casadi_ivp_output(casadi_ivp_data<T1>* d) {
  // State at t_out <= t, in or at the end of the last step
  const casadi_ivp_prob<T1>* p = d->prob;
  if (d->t == d->t_out) {
    casadi_copy(d->y, p->ny, d->y_out);
  } else {
    casadi_ivp_interpolate(d, d->t_out, d->y_out);
  }
}

// SYMBOL "ivp_accept"
template<typename T1>
void casadi_ivp_accept(casadi_ivp_data<T1>* d) {
  // The step from t0, with y the solution at its end, is accepted
  d->run_k++;
  d->t = d->run_k == d->run_n ? d->run_t1 : d->run_t0 + static_cast<T1>(d->run_k) * d->run_h;
  if (d->taping) casadi_ivp_tape_push(d);
  d->n_steps++;
}

// SYMBOL "ivp_back_begin"
template<typename T1>
int casadi_ivp_back_begin(casadi_ivp_data<T1>* d) {
  // Record to process next; returns 0 when the retreat is done (tape exhausted, or the
  // record belongs to an earlier output interval and becomes pending)
  for (;;) {
    if (!d->record_pending && casadi_ivp_tape_pop(d)) return casadi_ivp_stop(d, IVP_SUCCESS);
    // Records of later output intervals are skipped with zero adjoint
    if (d->t_record > d->t_start) {
      d->record_pending = 0;
      continue;
    }
    if (d->t_record < d->t_start) {
      d->record_pending = 1;
      d->t = d->t_out;
      return casadi_ivp_stop(d, IVP_SUCCESS);
    }
    d->record_pending = 0;
    return 1;
  }
}

// SYMBOL "ivp_back_step_begin"
template<typename T1>
void casadi_ivp_back_step_begin(casadi_ivp_data<T1>* d) {
  // Stage adjoints from y1 = y0 + sum_i m_i U_i, added to the seeds of dense outputs
  const casadi_ivp_prob<T1>* p = d->prob;
  casadi_int i;
  for (i = 0; i < p->n_stages; ++i) casadi_axpy(p->nr, p->m[i], d->ry, d->rU + i * p->nr);
  casadi_axpy(p->nr, 1., d->ry0, d->ry);
  casadi_axpy(p->nru, 1., d->ru_step, d->adj_u);
  d->stage = p->n_stages;
}

// SYMBOL "ivp_impulse_begin"
template<typename T1>
int casadi_ivp_impulse_begin(casadi_ivp_data<T1>* d) {
  // Pops the record containing t_seed (records of later intervals had zero seeds); returns
  // 1 if t_seed is inside its step, with the dense output at t_seed in y_out
  const casadi_ivp_prob<T1>* p = d->prob;
  while (!d->record_pending && !casadi_ivp_tape_pop(d)) {
    d->record_pending = d->t_record <= d->t_seed;
  }
  if (!d->record_pending || d->t_seed == d->t1) {
    // Seed of the state at the end of the pending record
    casadi_axpy(p->nr, 1., d->seed, d->ry);
    return casadi_ivp_stop(d, IVP_SUCCESS);
  }
  casadi_ivp_interpolate(d, d->t_seed, d->y_out);
  return 1;
}

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

// C-REPLACE "static_cast<int>" "(int) "
// C-REPLACE "static_cast<T1>" "(casadi_real) "
// C-REPLACE "static_cast<casadi_int>" "(casadi_int) "
// C-REPLACE "static_cast<casadi_ivp_record_t>" "(casadi_ivp_record_t) "

// SYMBOL "ivp_jac_t"
typedef enum {
  // Jacobian blocks w.r.t. x, z and t, in the order of evaluation by the caller
  IVP_ODE_X, IVP_ALG_X, IVP_QUAD_X,
  IVP_ODE_Z, IVP_ALG_Z, IVP_QUAD_Z,
  IVP_ODE_T, IVP_ALG_T, IVP_QUAD_T,
  IVP_NUM_JAC
} casadi_ivp_jac_t;

// SYMBOL "ivp_task_t"
typedef enum {
  // res = F(t_arg, y_arg), with forward sensitivities
  IVP_EVAL_RHS,
  // jac, jac_fwd = Jacobian blocks at (t_arg, y_arg), with forward sensitivities
  IVP_EVAL_JAC,
  // vjp = dF/d(x,z,p,u)^T w at (t_arg, y_arg), with forward sensitivities
  IVP_EVAL_VJP,
  // vjp = d/d(x,z,p,u) <jac_adj, Jacobian blocks> at (t_arg, y_arg), w. forward sensitivities
  IVP_EVAL_JAC_ADJ,
  // res = zero-crossing functions at (t_arg, y_arg), edot = derivative in direction (1, ydot)
  IVP_EVAL_ZERO,
  // res = y after the event selected by ev_w at (t_arg, y_arg), with forward sensitivities
  IVP_EVAL_JUMP,
  // vjp = adjoint of IVP_EVAL_JUMP with seed w, with forward sensitivities
  IVP_EVAL_JUMP_ADJ,
  // Enlarge the tape to at least tape_need entries, keeping its contents; set tape, tape_cap
  IVP_TAPE_RESERVE,
  // Adaptive step size control, see casadi_ivp_control
  IVP_STEP_CONTROL,
  // Events in the accepted step, see casadi_ivp_event
  IVP_EVENT,
  // Adjoint of the popped event record, see casadi_ivp_event_adj
  IVP_EVENT_ADJ
} casadi_ivp_task_t;

// SYMBOL "ivp_control_t"
typedef enum {
  // Size h of a step from t; h = 0 if an initial estimate from F(t, y) is needed first
  IVP_CONTROL_SIZE,
  // Initial estimate from F(t, y) in U, then as IVP_CONTROL_SIZE
  IVP_CONTROL_INIT,
  // Error test of the candidate y: accepted, or rejected with y restored
  IVP_CONTROL_TEST,
  // Evaluation failed: rejected with y restored
  IVP_CONTROL_FAIL
} casadi_ivp_control_t;

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
  IVP_STEP_TOO_SMALL,
  IVP_INIT_FAILED,
  IVP_EVAL_FAILED,
  IVP_TAPE_FAILED,
  IVP_TOO_MANY_EVENTS,
  IVP_EVENT_FAILED
} casadi_ivp_flag_t;

// SYMBOL "ivp_record_t"
typedef enum {
  // Integrator step from t0 to t0 + h
  IVP_RECORD_STEP,
  // Consistent initialization of the algebraic states at t0
  IVP_RECORD_INIT,
  // Event at t0: state before it, index in place of h
  IVP_RECORD_EVENT
} casadi_ivp_record_t;

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
  // Sparsity patterns of the Jacobian blocks (null: not used), see casadi_ivp_jac_t
  const casadi_int* sp_jac[IVP_NUM_JAC];
  // Capacity of the tape in the work vector (0: the caller provides the tape)
  casadi_int tape_size;
  // Fixed steps: each run of steps up to a stop has equal steps of at most step_size
  // (0: adaptive steps)
  T1 step_size;
  // Error weights: leading entries of x, z, q subject to error control, tolerances
  casadi_int nx_err, nz_err, nq_err;
  T1 abstol, reltol;
  // Adaptive steps: order and weights of the error estimate sum_i m_err_i U_i, step size
  // controller, maximum number of steps per output interval
  casadi_int order_err;
  const T1* m_err;
  T1 safety, fac_min, fac_max, h_min, h_max, h_init;
  casadi_int max_num_steps;
  // Zero-crossing events: number, transition given, tolerances on the event time (acceptable:
  // negative if none), iterations to locate an event, events per evaluation
  casadi_int ne;
  int has_transition;
  T1 event_tol, event_acceptable_tol;
  casadi_int max_event_iter, max_events;
  // Derived: forward and adjoint state dimensions, controls, vector lengths
  casadi_int nx, nz, nq, nu, nrx, nrz, nrq, nrp, nru, ny, nr;
  // Derived: offsets of the Jacobian blocks in a buffer with one direction, total nonzeros
  casadi_int jac_off[IVP_NUM_JAC + 1];
  casadi_int nnz_jac;
};
// C-REPLACE "casadi_ivp_prob<T1>" "struct casadi_ivp_prob"

// SYMBOL "ivp_nnz"
inline
casadi_int casadi_ivp_nnz(const casadi_int* sp) {
  // Number of nonzeros of a sparsity pattern, zero if null
  return sp ? sp[2 + sp[1]] : 0;
}

// SYMBOL "ivp_setup"
template<typename T1>
void casadi_ivp_setup(casadi_ivp_prob<T1>* p) {
  casadi_int k;
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
  for (k = 0; k < IVP_NUM_JAC; ++k) {
    p->jac_off[k + 1] = p->jac_off[k] + casadi_ivp_nnz(p->sp_jac[k]);
  }
  p->nnz_jac = p->jac_off[IVP_NUM_JAC];
}

// C-REPLACE "casadi_ivp_data<T1>" "struct casadi_ivp_data"
// SYMBOL "ivp_data"
template<typename T1>
struct casadi_ivp_data {
  // Problem structure
  const casadi_ivp_prob<T1>* prob;
  // Data whose task is requested: this one, or one the method delegates to
  casadi_ivp_data<T1>* active;
  // Nominal integration whose steps are replayed (or null), stop after every accepted step
  casadi_ivp_data<T1>* nom;
  int pause;
  // Task for the caller, where to resume (casadi_ivp_entry_t or a phase of the method), status
  casadi_ivp_task_t task;
  // Step size control requested, step accepted by it
  casadi_ivp_control_t control;
  int accepted;
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
  // Jacobian blocks and their forward derivatives, nonzeros in casadi_ivp_jac_t order
  T1 *jac, *jac_fwd;
  // Adjoint seeds of the Jacobian blocks, with forward derivatives
  T1 *jac_adj, *jac_adj_fwd;
  // Results of EVAL_VJP and EVAL_JAC_ADJ: adjoints of x, z, p, u
  T1 *vjp_x, *vjp_z, *vjp_p, *vjp_u;
  // Accumulated adjoints of the parameters and controls
  T1 *adj_p, *adj_u;
  // Tape: records [U (steps only), y0, u, h, t0, t1, t_out, type]; entries used, capacity, needed
  T1* tape;
  casadi_int tape_pos, tape_cap, tape_need;
  // Reset pending, taping requested
  int reset, taping;
  // Record type, popped record that belongs to an earlier interval
  casadi_ivp_record_t record;
  int record_pending;
  // Fixed steps: current run start, stop, step size, number of steps, steps taken
  T1 run_t0, run_t1, run_h;
  casadi_int run_n, run_k;
  // Adaptive steps: proposed next step size, previous error, last step rejected, clipped at
  // the stop or ending at the stop
  T1 h_next, err_prev;
  int last_rejected, last_clipped, last_step;
  // Steps taken in the current output interval
  casadi_int n_steps_interval;
  // Events: zero-crossing functions at the bracket ends (e0 also at the step start) and at a
  // trial time, candidates, time derivatives at these, selection of the event, state derivative
  T1 *e0, *e1, *em, *e_mask, *edot0, *edot1, *edot, *ev_w, *ydot;
  // Event phase, e0 known at t, state changed by an event (algebraic states inconsistent),
  // restart after an event
  int ev_phase, e_valid, jump, ev_restart;
  // Bracket, trial time, event time and its accuracy, the constant 1
  T1 ev_ta, ev_tb, ev_tm, t_ev, ev_dt, ev_one;
  // Iterations, event in the last step (-1: none), events in this evaluation
  casadi_int ev_iter, ev_fired, n_events;
  // Statistics
  casadi_int n_steps, n_reject, n_rhs, n_jac, n_vjp;
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
  // Jacobian blocks, forward derivatives, adjoint seeds and their forward derivatives
  *sz_w += p->nnz_jac * (1 + p->nfwd) * (1 + p->nadj);
  // Function adjoints, accumulated adjoints
  *sz_w += p->nrx + p->nrz + 2 * (p->nrp + p->nru);
  // Tape
  *sz_w += p->tape_size;
  // Events
  if (p->ne > 0) *sz_w += 8 * p->ne + p->ny;
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
  d->jac = *w; *w += p->nnz_jac;
  d->jac_fwd = *w; *w += p->nnz_jac * p->nfwd;
  d->jac_adj = *w; *w += p->nnz_jac * p->nadj;
  d->jac_adj_fwd = *w; *w += p->nnz_jac * p->nadj * p->nfwd;
  d->vjp_x = *w; *w += p->nrx;
  d->vjp_z = *w; *w += p->nrz;
  d->vjp_p = *w; *w += p->nrp;
  d->vjp_u = *w; *w += p->nru;
  d->adj_p = *w; *w += p->nrp;
  d->adj_u = *w; *w += p->nru;
  d->tape = *w; *w += p->tape_size;
  if (p->ne > 0) {
    d->e0 = *w; *w += p->ne;
    d->e1 = *w; *w += p->ne;
    d->em = *w; *w += p->ne;
    d->e_mask = *w; *w += p->ne;
    d->edot0 = *w; *w += p->ne;
    d->edot1 = *w; *w += p->ne;
    d->edot = *w; *w += p->ne;
    d->ev_w = *w; *w += p->ne;
    d->ydot = *w; *w += p->ny;
  }
  d->tape_cap = p->tape_size;
  d->tape_pos = 0;
  d->next = IVP_DONE;
  d->status = IVP_SUCCESS;
  d->fail = 0;
  d->reset = 0;
  d->taping = 0;
  d->record_pending = 0;
  d->active = d;
  d->nom = 0;
  d->pause = 0;
  d->ev_phase = 0;
  d->e_valid = 0;
  d->jump = 0;
  d->ev_restart = 0;
  d->ev_dt = 0;
  d->ev_one = 1;
  d->ev_fired = -1;
  d->n_events = 0;
}

// SYMBOL "ivp_clear_stats"
template<typename T1>
void casadi_ivp_clear_stats(casadi_ivp_data<T1>* d) {
  d->n_steps = d->n_reject = d->n_rhs = d->n_jac = d->n_vjp = 0;
}

// SYMBOL "ivp_reset"
template<typename T1>
void casadi_ivp_reset(casadi_ivp_data<T1>* d, T1 t, const T1* u, int first_call) {
  // Start at time t, controls u; the state d->y is set by the caller
  const casadi_ivp_prob<T1>* p = d->prob;
  const casadi_ivp_prob<T1>* pn;
  d->t = t;
  d->u = u;
  d->reset = 1;
  d->e_valid = 0;
  if (first_call) {
    d->n_events = 0;
    d->ev_restart = 0;
    d->tape_pos = 0;
    casadi_ivp_clear_stats(d);
    // Every evaluation starts from the same controller state
    d->h_next = 0;
    d->err_prev = 1;
    d->last_rejected = 0;
  }
  if (d->nom) {
    // Initial nominal states: the leading entries of x, z and q; afterwards, the nominal
    // integration continues its own trajectory
    pn = d->nom->prob;
    if (first_call) {
      casadi_copy(d->y, pn->nx, d->nom->y);
      casadi_copy(d->y + p->nx, pn->nz, d->nom->y + pn->nx);
      casadi_copy(d->y + p->nx + p->nz, pn->nq, d->nom->y + pn->nx + pn->nz);
    }
    casadi_ivp_reset(d->nom, t, u, first_call);
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
  // Handle a pending reset: with fixed steps, a new run of steps up to the stop time
  const casadi_ivp_prob<T1>* p = d->prob;
  if (p->step_size > 0) {
    d->run_t0 = d->t;
    d->run_t1 = d->t_stop;
    d->run_n = casadi_ivp_steps(d->t_stop - d->t, p->step_size);
    d->run_h = (d->run_t1 - d->run_t0) / static_cast<T1>(d->run_n);
    d->run_k = 0;
  }
  d->reset = 0;
  d->jump = 0;
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
  d->n_steps_interval = 0;
  if (d->nom) {
    // The nominal integration stops after every step, to be replayed
    casadi_ivp_advance(d->nom, t_out, t_stop, 0);
    d->nom->pause = 1;
  }
}

// SYMBOL "ivp_reset_backward"
template<typename T1>
void casadi_ivp_reset_backward(casadi_ivp_data<T1>* d) {
  const casadi_ivp_prob<T1>* p = d->prob;
  casadi_clear(d->ry, p->nr);
  casadi_clear(d->adj_p, p->nrp);
  casadi_clear(d->adj_u, p->nru);
  casadi_clear(d->jac_adj, p->nnz_jac * p->nadj * (1 + p->nfwd));
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
  d->active = d;
  d->task = task;
  d->next = next;
  d->fail = 0;
  return 1;
}

// SYMBOL "ivp_request_control"
template<typename T1>
int casadi_ivp_request_control(casadi_ivp_data<T1>* d, casadi_ivp_control_t control,
    casadi_int next) {
  // Hand step size control to the caller, resume at phase next
  d->control = control;
  return casadi_ivp_request(d, IVP_STEP_CONTROL, next);
}

// SYMBOL "ivp_stop"
template<typename T1>
int casadi_ivp_stop(casadi_ivp_data<T1>* d, casadi_ivp_flag_t status) {
  d->status = status;
  d->next = IVP_DONE;
  return 0;
}

// SYMBOL "ivp_jac_nz"
template<typename T1>
T1* casadi_ivp_jac_nz(const casadi_ivp_data<T1>* d, T1* jac, casadi_int k,
    casadi_int ndir, casadi_int dir) {
  // Nonzeros of Jacobian block k, direction dir, in a buffer with ndir directions per block
  const casadi_int* off = d->prob->jac_off;
  return jac + ndir * off[k] + dir * (off[k + 1] - off[k]);
}

// SYMBOL "ivp_jac_mv"
template<typename T1>
void casadi_ivp_jac_mv(const casadi_ivp_data<T1>* d, casadi_int k,
    casadi_int dir, const T1* v, T1* r, int tr) {
  // r += J_k v (or J_k^T v): Jacobian block k (dir < 0) or its forward derivative
  const casadi_ivp_prob<T1>* p = d->prob;
  const T1* nz = dir < 0 ? casadi_ivp_jac_nz(d, d->jac, k, 1, 0)
    : casadi_ivp_jac_nz(d, d->jac_fwd, k, p->nfwd, dir);
  casadi_mv(nz, p->sp_jac[k], v, r, tr);
}

// SYMBOL "ivp_add_vjp"
template<typename T1>
void casadi_ivp_add_vjp(casadi_ivp_data<T1>* d, T1* r) {
  // Add x, z adjoints of EVAL_VJP/EVAL_JAC_ADJ to r, p, u adjoints to the running sums
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
// Tape: records [U (step records only), y0, u, h, t0, t1, end of the output interval, type]
// ---------------------------------------------------------------------------------------

// SYMBOL "ivp_tape_full"
template<typename T1>
int casadi_ivp_tape_full(casadi_ivp_data<T1>* d, casadi_ivp_record_t record) {
  // No room for a record of this type (a step with events: and an event record)? Sets the
  // capacity needed
  const casadi_ivp_prob<T1>* p = d->prob;
  casadi_int n = (record == IVP_RECORD_STEP ? p->ny * p->n_stages : 0) + p->ny + p->nu + 5;
  if (record == IVP_RECORD_STEP && p->ne > 0) n += p->ny + p->nu + 5;
  d->tape_need = d->tape_pos + n;
  return n > d->tape_cap - d->tape_pos;
}

// SYMBOL "ivp_reserve"
template<typename T1>
int casadi_ivp_reserve(casadi_ivp_data<T1>* d, casadi_ivp_record_t record, casadi_int next) {
  // Room on the tape for a record, if taping? -1 if so, else the value for the method to return
  if (!d->taping || !casadi_ivp_tape_full(d, record)) return -1;
  if (d->prob->tape_size > 0) return casadi_ivp_stop(d, IVP_TAPE_FAILED);
  return casadi_ivp_request(d, IVP_TAPE_RESERVE, next);
}

// SYMBOL "ivp_tape_push"
template<typename T1>
void casadi_ivp_tape_push(casadi_ivp_data<T1>* d, const T1* y0, T1 h, T1 t0, T1 t1) {
  // Append a record of type d->record from t0 to t1, for which room has been reserved
  const casadi_ivp_prob<T1>* p = d->prob;
  casadi_int n_u = d->record == IVP_RECORD_STEP ? p->ny * p->n_stages : 0;
  T1* rec = d->tape + d->tape_pos;
  casadi_copy(d->U, n_u, rec);
  rec += n_u;
  casadi_copy(y0, p->ny, rec);
  casadi_copy(d->u_step, p->nu, rec + p->ny);
  rec += p->ny + p->nu;
  rec[0] = h;
  rec[1] = t0;
  rec[2] = t1;
  rec[3] = d->t_out;
  rec[4] = d->record;
  d->tape_pos += n_u + p->ny + p->nu + 5;
}

// SYMBOL "ivp_tape_pop"
template<typename T1>
int casadi_ivp_tape_pop(casadi_ivp_data<T1>* d) {
  // Restore the last record, clear the seeds of its stages; returns 1 if the tape is empty
  const casadi_ivp_prob<T1>* p = d->prob;
  casadi_int n_u;
  T1* rec;
  if (d->tape_pos < p->ny + p->nu + 5) return 1;
  rec = d->tape + d->tape_pos;
  d->record = static_cast<casadi_ivp_record_t>(static_cast<int>(rec[-1]));
  d->t_record = rec[-2];
  d->t1 = rec[-3];
  d->t0 = rec[-4];
  d->h = rec[-5];
  n_u = d->record == IVP_RECORD_STEP ? p->ny * p->n_stages : 0;
  rec -= p->ny + p->nu + 5;
  casadi_copy(rec, p->ny, d->y0);
  casadi_copy(rec + p->ny, p->nu, d->u_step);
  casadi_copy(rec - n_u, n_u, d->U);
  d->tape_pos -= n_u + p->ny + p->nu + 5;
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
// Step sizes; adaptive control is in casadi_ivp_adapt.hpp
// ---------------------------------------------------------------------------------------

// SYMBOL "ivp_isfinite"
template<typename T1>
int casadi_ivp_isfinite(T1 x) {
  // False for +-inf and NaN
  return x - x == 0;
}

// SYMBOL "ivp_step_limit"
template<typename T1>
int casadi_ivp_step_limit(const casadi_ivp_data<T1>* d) {
  // No further step allowed?
  const casadi_ivp_prob<T1>* p = d->prob;
  if (p->step_size > 0 && !d->nom) return d->run_k >= d->run_n;
  return d->n_steps_interval >= p->max_num_steps;
}

// SYMBOL "ivp_prescribed"
template<typename T1>
int casadi_ivp_prescribed(const casadi_ivp_data<T1>* d) {
  // Step sizes prescribed: fixed steps, or replayed from the nominal integration
  return d->prob->step_size > 0 || d->nom != 0;
}

// ---------------------------------------------------------------------------------------
// Building blocks of the forward integration and of the discrete adjoint
// ---------------------------------------------------------------------------------------

// SYMBOL "ivp_output"
template<typename T1>
int casadi_ivp_output(casadi_ivp_data<T1>* d) {
  // State at t_out <= t, in or at the end of the last step; returns 1 if interpolated
  const casadi_ivp_prob<T1>* p = d->prob;
  if (d->t == d->t_out) {
    casadi_copy(d->y, p->ny, d->y_out);
    return 0;
  }
  casadi_ivp_interpolate(d, d->t_out, d->y_out);
  return 1;
}

// SYMBOL "ivp_accept"
template<typename T1>
void casadi_ivp_accept(casadi_ivp_data<T1>* d) {
  // The step from t0, with y the solution at its end, is accepted
  const casadi_ivp_prob<T1>* p = d->prob;
  if (d->nom) {
    d->t = d->nom->t;
  } else if (p->step_size > 0) {
    d->run_k++;
    d->t = d->run_k == d->run_n ? d->run_t1 : d->run_t0 + static_cast<T1>(d->run_k) * d->run_h;
  } else {
    d->t = d->last_step ? d->t_stop : d->t0 + d->h;
  }
  d->last_rejected = 0;
  d->n_steps_interval++;
  if (d->taping) {
    d->record = IVP_RECORD_STEP;
    casadi_ivp_tape_push(d, d->y0, d->h, d->t0, d->t);
  }
  d->n_steps++;
}

// SYMBOL "ivp_back_begin"
template<typename T1>
int casadi_ivp_back_begin(casadi_ivp_data<T1>* d) {
  // Record to process next; returns 0 when the retreat is done (tape exhausted, or the
  // record belongs to an earlier output interval and becomes pending). An event belongs to the
  // output interval that contains it
  T1 r;
  if (d->status != IVP_SUCCESS) return casadi_ivp_stop(d, d->status);
  for (;;) {
    if (!d->record_pending && casadi_ivp_tape_pop(d)) return casadi_ivp_stop(d, IVP_SUCCESS);
    r = d->record == IVP_RECORD_EVENT ? d->t0 : d->t_record;
    // Records of later output intervals are skipped with zero adjoint
    if (r > d->t_start) {
      d->record_pending = 0;
      continue;
    }
    if (d->record == IVP_RECORD_EVENT ? r <= d->t_out : r < d->t_start) {
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
    d->record_pending = (d->record == IVP_RECORD_EVENT ? d->t0 : d->t_record) <= d->t_seed;
  }
  if (!d->record_pending || d->record != IVP_RECORD_STEP || d->t_seed == d->t1) {
    // Seed of the state at the end of the pending record
    casadi_axpy(p->nr, 1., d->seed, d->ry);
    return casadi_ivp_stop(d, IVP_SUCCESS);
  }
  casadi_ivp_interpolate(d, d->t_seed, d->y_out);
  return 1;
}

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

// Explicit Runge-Kutta integrator for y' = F(t, y), y = [x; q], see casadi_ivp.hpp
// Stages U_i = h F(t0 + alpha_i h, y0 + sum_{j<i} a_ij U_j), y1 = y0 + sum_i m_i U_i
// Fixed or adaptive steps; with the first-same-as-last property, the last stage of a step is the
// first of the next

// C-REPLACE "casadi_ivp_prob<T1>" "struct casadi_ivp_prob"
// C-REPLACE "casadi_ivp_data<T1>" "struct casadi_ivp_data"

// SYMBOL "erk_prob"
template<typename T1>
struct casadi_erk_prob {
  // Last stage at t0 + h, y1 (first same as last)
  int fsal;
};
// C-REPLACE "casadi_erk_prob<T1>" "struct casadi_erk_prob"

// SYMBOL "erk_next_t"
typedef enum {
  // Phases: each one consumes the result of the previously requested task
  ERK_STEP_BEGIN = IVP_NUM_ENTRY,
  ERK_STEP_SIZE,
  ERK_STEP_INITIAL_SIZE,
  ERK_STEP_REPLAY,
  ERK_STAGE,
  ERK_STEP_CHECK,
  ERK_STEP_EVENT,
  ERK_BACK_STAGE,
  ERK_BACK_SPREAD
} casadi_erk_next_t;

// C-REPLACE "casadi_erk_data<T1>" "struct casadi_erk_data"
// SYMBOL "erk_data"
template<typename T1>
struct casadi_erk_data {
  // Problem structure
  const casadi_erk_prob<T1>* prob;
  // Method data of the nominal integration whose steps are replayed (or null)
  casadi_erk_data<T1>* nom;
  // Stage holding h_f F(t, y) for the next step (-1: none), h_f
  casadi_int f_stage;
  T1 f_scale;
};

// SYMBOL "erk"
template<typename T1>
int casadi_erk(casadi_ivp_data<T1>* d, casadi_erk_data<T1>* e);

// SYMBOL "erk_stage_request"
template<typename T1>
int casadi_erk_stage_request(casadi_ivp_data<T1>* d) {
  // Right-hand side at the argument of the current stage
  const casadi_ivp_prob<T1>* p = d->prob;
  casadi_ivp_stage_argument(d, d->stage);
  d->t_arg = d->t0 + p->alpha[d->stage] * d->h;
  d->y_arg = d->y_stage;
  d->res = d->U + d->stage * p->ny;
  return casadi_ivp_request(d, IVP_EVAL_RHS, ERK_STAGE);
}

// SYMBOL "erk_step_start"
template<typename T1>
int casadi_erk_step_start(casadi_ivp_data<T1>* d, casadi_erk_data<T1>* e) {
  // Step of size h from t, reusing F(t, y) if available
  const casadi_ivp_prob<T1>* p = d->prob;
  d->t0 = d->t;
  casadi_copy(d->y, p->ny, d->y0);
  d->stage = 0;
  if (e->f_stage >= 0) {
    if (e->f_stage > 0) casadi_copy(d->U + e->f_stage * p->ny, p->ny, d->U);
    casadi_scal(p->ny, d->h / e->f_scale, d->U);
    e->f_stage = -1;
    d->stage = 1;
  }
  return casadi_erk_stage_request(d);
}

// SYMBOL "erk_forward"
template<typename T1>
int casadi_erk_forward(casadi_ivp_data<T1>* d, casadi_erk_data<T1>* e) {
  // Forward integration; returns 1 if a task is requested, 0 when done or failed
  const casadi_ivp_prob<T1>* p = d->prob;
  int ret;
  switch (d->next) {
  case IVP_ADVANCE:
    if (d->reset) {
      // New state or controls: F(t, y) is not known
      casadi_ivp_begin_run(d);
      e->f_stage = -1;
    }
    d->next = ERK_STEP_BEGIN;
    return casadi_erk_forward(d, e);
  case ERK_STEP_BEGIN:
    if (d->fail) return casadi_ivp_stop(d, IVP_TAPE_FAILED);
    if (d->t >= d->t_out) {
      if (!d->pause) casadi_ivp_output(d);
      return casadi_ivp_stop(d, IVP_SUCCESS);
    }
    if (d->jump) {
      // State changed by an event: restart
      d->reset = 1;
      d->next = IVP_ADVANCE;
      return casadi_erk_forward(d, e);
    }
    if (casadi_ivp_step_limit(d)) return casadi_ivp_stop(d, IVP_MAX_NUM_STEPS);
    if ((ret = casadi_ivp_reserve(d, IVP_RECORD_STEP, ERK_STEP_BEGIN)) >= 0) return ret;
    casadi_copy(d->u, p->nu, d->u_step);
    if (d->nom) {
      // Step of the nominal integration
      d->next = ERK_STEP_REPLAY;
      return casadi_erk_forward(d, e);
    }
    if (d->prob->step_size > 0) {
      d->h = d->run_h;
      return casadi_erk_step_start(d, e);
    }
    return casadi_ivp_request_control(d, IVP_CONTROL_SIZE, ERK_STEP_SIZE);
  case ERK_STEP_SIZE:
    if (d->h > 0) return casadi_erk_step_start(d, e);
    // Estimate from the right-hand side at the initial point
    d->t_arg = d->t;
    d->y_arg = d->y;
    d->res = d->U;
    return casadi_ivp_request(d, IVP_EVAL_RHS, ERK_STEP_INITIAL_SIZE);
  case ERK_STEP_INITIAL_SIZE:
    if (d->fail) return casadi_ivp_stop(d, IVP_EVAL_FAILED);
    d->n_rhs++;
    e->f_stage = 0;
    e->f_scale = 1;
    return casadi_ivp_request_control(d, IVP_CONTROL_INIT, ERK_STEP_SIZE);
  case ERK_STEP_REPLAY:
    // Advance the nominal integration by one step; its tasks are handed to the caller
    if (casadi_erk(d->nom, e->nom)) {
      d->active = d->nom;
      return 1;
    }
    d->active = d;
    if (d->nom->status != IVP_SUCCESS) return casadi_ivp_stop(d, d->nom->status);
    if (d->nom->next == IVP_DONE) return casadi_ivp_stop(d, IVP_EVAL_FAILED);
    d->h = d->nom->h;
    return casadi_erk_step_start(d, e);
  case ERK_STAGE:
    // F(t0, y0), if evaluated, holds for a retry from the same point
    e->f_stage = d->stage > 0 ? 0 : -1;
    e->f_scale = d->h;
    if (d->fail) {
      // Evaluation failed: retry with a smaller step, if adaptive
      if (casadi_ivp_prescribed(d)) return casadi_ivp_stop(d, IVP_EVAL_FAILED);
      return casadi_ivp_request_control(d, IVP_CONTROL_FAIL, ERK_STEP_CHECK);
    }
    d->n_rhs++;
    casadi_scal(p->ny, d->h, d->U + d->stage * p->ny);
    if (++d->stage < p->n_stages) return casadi_erk_stage_request(d);
    if (e->prob->fsal) {
      // y1 is the argument of the last stage
      casadi_copy(d->y_stage, p->ny, d->y);
    } else {
      // y = y0 + sum_i m_i U_i
      casadi_copy(d->y0, p->ny, d->y);
      casadi_mv_dense(d->U, p->ny, p->n_stages, p->m, d->y, 0);
    }
    if (!casadi_ivp_prescribed(d)) {
      return casadi_ivp_request_control(d, IVP_CONTROL_TEST, ERK_STEP_CHECK);
    }
    d->accepted = 1;
    d->next = ERK_STEP_CHECK;
    return casadi_erk_forward(d, e);
  case ERK_STEP_CHECK:
    if (d->status != IVP_SUCCESS) return casadi_ivp_stop(d, d->status);
    if (!d->accepted) {
      // Retry from t0
      d->next = ERK_STEP_BEGIN;
      return casadi_erk_forward(d, e);
    }
    e->f_stage = e->prob->fsal ? p->n_stages - 1 : -1;
    casadi_ivp_accept(d);
    if (p->ne > 0) return casadi_ivp_request(d, IVP_EVENT, ERK_STEP_EVENT);
    d->next = ERK_STEP_EVENT;
    return casadi_erk_forward(d, e);
  case ERK_STEP_EVENT:
    if (d->status != IVP_SUCCESS) return casadi_ivp_stop(d, d->status);
    d->next = ERK_STEP_BEGIN;
    if (d->pause) {
      // Hand the step to the integration that replays it
      d->status = IVP_SUCCESS;
      return 0;
    }
    return casadi_erk_forward(d, e);
  default:
    return 0;
  }
}

// SYMBOL "erk_backward"
template<typename T1>
int casadi_erk_backward(casadi_ivp_data<T1>* d, casadi_erk_data<T1>* e) {
  // Discrete adjoint of the taped steps; returns 1 if a task is requested, 0 when done
  const casadi_ivp_prob<T1>* p = d->prob;
  casadi_int i, j, s = p->n_stages;
  T1* rUi;
  switch (d->next) {
  case IVP_RETREAT:
    if (!casadi_ivp_back_begin(d)) return 0;
    if (d->record == IVP_RECORD_EVENT) return casadi_ivp_request(d, IVP_EVENT_ADJ, IVP_RETREAT);
    casadi_ivp_back_step_begin(d);
    d->next = ERK_BACK_STAGE;
    return casadi_erk_backward(d, e);
  case ERK_BACK_STAGE:
    if (d->stage == 0) {
      // Next record
      d->t = d->t0;
      d->next = IVP_RETREAT;
      return casadi_erk_backward(d, e);
    }
    // Previous stage: U_i = h F at the stage argument
    i = --d->stage;
    rUi = d->rU + i * p->nr;
    casadi_scal(p->nr, d->h, rUi);
    casadi_ivp_stage_argument(d, i);
    casadi_ivp_clear_vjp(d);
    d->t_arg = d->t0 + p->alpha[i] * d->h;
    d->y_arg = d->y_stage;
    d->w = rUi;
    return casadi_ivp_request(d, IVP_EVAL_VJP, ERK_BACK_SPREAD);
  case ERK_BACK_SPREAD:
    // Adjoint of y0 + sum_j a_ij U_j to y0 and earlier stages
    if (d->fail) return casadi_ivp_stop(d, IVP_EVAL_FAILED);
    d->n_vjp++;
    i = d->stage;
    casadi_ivp_add_vjp(d, d->ry);
    for (j = 0; j < i; ++j) {
      if (p->a[i * s + j] != 0) casadi_axpy(p->nrx, p->a[i * s + j], d->vjp_x, d->rU + j * p->nr);
    }
    d->next = ERK_BACK_STAGE;
    return casadi_erk_backward(d, e);
  case IVP_IMPULSE:
    // Adjoint seed of the output at t_seed
    if (casadi_ivp_impulse_begin(d)) {
      casadi_ivp_seed_stages(d);
      return casadi_ivp_stop(d, IVP_SUCCESS);
    }
    return 0;
  default:
    return 0;
  }
}

template<typename T1>
int casadi_erk(casadi_ivp_data<T1>* d, casadi_erk_data<T1>* e) {
  // Reverse communication: returns 1 while a task in d->task is to be carried out
  if (d->next == IVP_IMPULSE || d->next == IVP_RETREAT || d->next >= ERK_BACK_STAGE) {
    return casadi_erk_backward(d, e);
  }
  return casadi_erk_forward(d, e);
}

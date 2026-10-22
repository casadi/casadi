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

// Zero-crossing events of the one-step integrators, see casadi_ivp.hpp: the kernel requests
// task IVP_EVENT after each accepted step, so only integrators with events need this file
// Event i occurs where e_i crosses from positive to negative; it is located on the dense output
// of the step, and the integration restarts from the state after the transition
// The jump (task IVP_EVAL_JUMP) is a function of the state before the event whose derivatives
// are those of the event, including the dependence of its time: forward, adjoint (task
// IVP_EVENT_ADJ for event records of the tape) and forward-over-adjoint sensitivities

// C-REPLACE "static_cast<T1>" "(casadi_real) "
// C-REPLACE "static_cast<casadi_int>" "(casadi_int) "
// C-REPLACE "fabs" "casadi_fabs"
// C-REPLACE "casadi_ivp_prob<T1>" "struct casadi_ivp_prob"
// C-REPLACE "casadi_ivp_data<T1>" "struct casadi_ivp_data"

// SYMBOL "ivp_event_phase_t"
typedef enum {
  // Phases of casadi_ivp_event: each one consumes the result of the previously requested task
  IVP_EV_BEGIN,
  IVP_EV_START,
  IVP_EV_END,
  IVP_EV_TRIAL,
  IVP_EV_POST,
  IVP_EV_BACK_SEED
} casadi_ivp_event_phase_t;

// SYMBOL "ivp_event_request"
template<typename T1>
int casadi_ivp_event_request(casadi_ivp_data<T1>* d, casadi_ivp_task_t task, int phase) {
  // Hand a task to the caller, resume at phase
  d->task = task;
  d->ev_phase = phase;
  d->fail = 0;
  return 1;
}

// SYMBOL "ivp_event_done"
template<typename T1>
int casadi_ivp_event_done(casadi_ivp_data<T1>* d, casadi_ivp_flag_t status) {
  d->ev_phase = IVP_EV_BEGIN;
  if (status != IVP_SUCCESS) d->status = status;
  return 0;
}

// SYMBOL "ivp_dense_derivative"
template<typename T1>
void casadi_ivp_dense_derivative(casadi_ivp_data<T1>* d, T1 t, T1* ydot) {
  // Time derivative of the dense output of the current step at t
  const casadi_ivp_prob<T1>* p = d->prob;
  casadi_int q;
  T1 theta = (t - d->t0) / d->h, th = 1. / d->h;
  casadi_clear(d->w_dense, p->n_stages);
  for (q = 0; q < p->n_dense; ++q) {
    casadi_axpy(p->n_stages, static_cast<T1>(q + 1) * th, p->m_dense + q * p->n_stages,
      d->w_dense);
    th *= theta;
  }
  casadi_clear(ydot, p->ny);
  casadi_mv_dense(d->U, p->ny, p->n_stages, d->w_dense, ydot, 0);
}

// SYMBOL "ivp_event_hermite"
template<typename T1>
casadi_int casadi_ivp_event_hermite(casadi_ivp_data<T1>* d, T1* t) {
  // First crossing in the bracket of the cubic Hermite interpolants of the candidates, from
  // values and derivatives at its ends; its time in t; -1 if none
  const casadi_ivp_prob<T1>* p = d->prob;
  casadi_int i, k, ind = -1;
  T1 H = d->ev_tb - d->ev_ta, c[4], sa, sb, s, v, t_i;
  for (i = 0; i < p->ne; ++i) {
    if (d->e_mask[i] == 0 || !(d->e1[i] < 0)) continue;
    // p(s) = c0 + c1 s + c2 s^2 + c3 s^3, p(0) >= 0 > p(1)
    c[0] = d->e0[i];
    c[1] = H * d->edot0[i];
    c[2] = 3 * (d->e1[i] - d->e0[i]) - H * (2 * d->edot0[i] + d->edot1[i]);
    c[3] = 2 * (d->e0[i] - d->e1[i]) + H * (d->edot0[i] + d->edot1[i]);
    // First sign change on a uniform grid, then bisection
    sa = 0;
    sb = 1;
    for (k = 1; k <= 8; ++k) {
      s = static_cast<T1>(k) / 8;
      if (c[0] + s * (c[1] + s * (c[2] + s * c[3])) < 0) {
        sb = s;
        break;
      }
      sa = s;
    }
    for (k = 0; k < 60 && sb - sa > 1e-15; ++k) {
      s = 0.5 * (sa + sb);
      v = c[0] + s * (c[1] + s * (c[2] + s * c[3]));
      if (v < 0) {
        sb = s;
      } else {
        sa = s;
      }
    }
    t_i = d->ev_ta + sb * H;
    if (ind < 0 || t_i < *t) {
      *t = t_i;
      ind = i;
    }
  }
  return ind;
}

// SYMBOL "ivp_event_eval"
template<typename T1>
int casadi_ivp_event_eval(casadi_ivp_data<T1>* d, T1 t, int phase) {
  // Zero-crossing functions and their time derivatives at t on the dense output of the step
  const casadi_ivp_prob<T1>* p = d->prob;
  d->ev_tm = t;
  if (t == d->t) {
    casadi_copy(d->y, p->ny, d->y_out);
  } else {
    casadi_ivp_interpolate(d, t, d->y_out);
  }
  casadi_ivp_dense_derivative(d, t, d->ydot);
  d->t_arg = t;
  d->y_arg = d->y_out;
  d->res = d->em;
  return casadi_ivp_event_request(d, IVP_EVAL_ZERO, phase);
}

// SYMBOL "ivp_event_restart"
template<typename T1>
int casadi_ivp_event_restart(casadi_ivp_data<T1>* d) {
  // Continue from the state after the event, adaptive steps from a new initial step size; the
  // method restarts once the outputs before the event, on the step that contains it, are done
  d->t = d->t_ev;
  d->e_valid = 0;
  d->ev_restart = 1;
  d->jump = 1;
  d->h_next = 0;
  d->err_prev = 1;
  d->last_rejected = 0;
  return casadi_ivp_event_done(d, IVP_SUCCESS);
}

// SYMBOL "ivp_event_select"
template<typename T1>
void casadi_ivp_event_select(casadi_ivp_data<T1>* d, casadi_int ind) {
  // Argument of the jump selecting event ind
  casadi_clear(d->ev_w, d->prob->ne);
  d->ev_w[ind] = 1;
}

// SYMBOL "ivp_event_post"
template<typename T1>
int casadi_ivp_event_post(casadi_ivp_data<T1>* d) {
  // State after the event in y: event record (room reserved with the step), restart
  if (d->taping) {
    // State before the event, index in place of the step size
    d->record = IVP_RECORD_EVENT;
    casadi_ivp_tape_push(d, d->y_out, static_cast<T1>(d->ev_fired), d->t_ev, d->t_ev);
  }
  return casadi_ivp_event_restart(d);
}

// SYMBOL "ivp_event_fire"
template<typename T1>
int casadi_ivp_event_fire(casadi_ivp_data<T1>* d) {
  // Event ev_fired at t_ev: state before it in y_out, then the jump
  const casadi_ivp_prob<T1>* p = d->prob;
  if (!d->nom && ++d->n_events > p->max_events) {
    return casadi_ivp_event_done(d, IVP_TOO_MANY_EVENTS);
  }
  if (!d->nom && d->t_ev == d->t) {
    casadi_copy(d->y, p->ny, d->y_out);
  } else {
    casadi_ivp_interpolate(d, d->t_ev, d->y_out);
  }
  casadi_copy(d->y_out, p->ny, d->y);
  if (!p->has_transition) return casadi_ivp_event_post(d);
  casadi_ivp_event_select(d, d->ev_fired);
  d->t_arg = d->t_ev;
  d->y_arg = d->y_out;
  d->res = d->y;
  return casadi_ivp_event_request(d, IVP_EVAL_JUMP, IVP_EV_POST);
}

// SYMBOL "ivp_event"
template<typename T1>
int casadi_ivp_event(casadi_ivp_data<T1>* d) {
  // Events in the accepted step from t0 to t (task IVP_EVENT): returns 1 while a task in
  // d->task is to be carried out; afterwards ev_fired >= 0 and jump set if one occurred
  const casadi_ivp_prob<T1>* p = d->prob;
  casadi_int i, ind;
  int any;
  T1 tn = 0, dt;
  switch (d->ev_phase) {
  case IVP_EV_BEGIN:
    d->ev_fired = -1;
    if (d->nom) {
      // Replayed step: the event of the nominal integration, if any
      if (d->nom->ev_fired < 0) return casadi_ivp_event_done(d, IVP_SUCCESS);
      d->ev_fired = d->nom->ev_fired;
      d->t_ev = d->nom->t_ev;
      return casadi_ivp_event_fire(d);
    }
    casadi_clear(d->e_mask, p->ne);
    if (!d->e_valid) return casadi_ivp_event_eval(d, d->t0, IVP_EV_START);
    return casadi_ivp_event_eval(d, d->t, IVP_EV_END);
  case IVP_EV_START:
    if (d->fail) return casadi_ivp_event_done(d, IVP_EVAL_FAILED);
    casadi_copy(d->em, p->ne, d->e0);
    casadi_copy(d->edot, p->ne, d->edot0);
    if (d->ev_restart) {
      // On the zero surface after an event, up to the accuracy of the event time: positive if
      // rising, so that an event cannot repeat until its function has been positive
      d->ev_restart = 0;
      for (i = 0; i < p->ne; ++i) {
        if (fabs(d->e0[i]) > 10 * fabs(d->edot0[i]) * d->ev_dt) continue;
        d->e_mask[i] = d->edot0[i] > 0;
        d->e0[i] = 0;
      }
    }
    return casadi_ivp_event_eval(d, d->t, IVP_EV_END);
  case IVP_EV_END:
    if (d->fail) return casadi_ivp_event_done(d, IVP_EVAL_FAILED);
    casadi_copy(d->em, p->ne, d->e1);
    casadi_copy(d->edot, p->ne, d->edot1);
    // Candidates: positive at the start of the step, negative at its end
    any = 0;
    for (i = 0; i < p->ne; ++i) {
      d->e_mask[i] = (d->e0[i] > 0 || d->e_mask[i] != 0) && d->e1[i] < 0;
      if (d->e_mask[i] != 0) any = 1;
    }
    if (!any) {
      casadi_copy(d->e1, p->ne, d->e0);
      casadi_copy(d->edot1, p->ne, d->edot0);
      d->e_valid = 1;
      return casadi_ivp_event_done(d, IVP_SUCCESS);
    }
    d->ev_ta = d->t0;
    d->ev_tb = d->t;
    d->ev_iter = 0;
    casadi_ivp_event_hermite(d, &tn);
    return casadi_ivp_event_eval(d, tn, IVP_EV_TRIAL);
  case IVP_EV_TRIAL:
    if (d->fail) return casadi_ivp_event_done(d, IVP_EVAL_FAILED);
    // Shrink the bracket: a candidate crossed before the trial time or not
    any = 0;
    for (i = 0; i < p->ne; ++i) if (d->e_mask[i] != 0 && d->em[i] < 0) any = 1;
    if (any) {
      d->ev_tb = d->ev_tm;
      casadi_copy(d->em, p->ne, d->e1);
      casadi_copy(d->edot, p->ne, d->edot1);
    } else {
      d->ev_ta = d->ev_tm;
      casadi_copy(d->em, p->ne, d->e0);
      casadi_copy(d->edot, p->ne, d->edot0);
    }
    // Next estimate; its distance to the trial time measures the error
    ind = casadi_ivp_event_hermite(d, &tn);
    dt = fabs(tn - d->ev_tm);
    d->ev_iter++;
    if (dt <= p->event_tol || d->ev_tb - d->ev_ta <= p->event_tol
        || d->ev_iter >= p->max_event_iter) {
      if (p->event_acceptable_tol >= 0 && dt >= p->event_acceptable_tol) {
        return casadi_ivp_event_done(d, IVP_EVENT_FAILED);
      }
      d->ev_fired = ind;
      d->t_ev = tn;
      d->ev_dt = dt > p->event_tol ? dt : p->event_tol;
      return casadi_ivp_event_fire(d);
    }
    return casadi_ivp_event_eval(d, tn, IVP_EV_TRIAL);
  case IVP_EV_POST:
    if (d->fail) return casadi_ivp_event_done(d, IVP_EVAL_FAILED);
    return casadi_ivp_event_post(d);
  default:
    return casadi_ivp_event_done(d, IVP_EVENT_FAILED);
  }
}

// SYMBOL "ivp_event_adj"
template<typename T1>
int casadi_ivp_event_adj(casadi_ivp_data<T1>* d) {
  // Adjoint of the popped event record (task IVP_EVENT_ADJ), then the stages of the step it
  // cut short are seeded through the dense output; returns 1 while a task in d->task is to be
  // carried out
  const casadi_ivp_prob<T1>* p = d->prob;
  switch (d->ev_phase) {
  case IVP_EV_BEGIN:
    d->t_ev = d->t0;
    if (p->has_transition) {
      casadi_ivp_event_select(d, static_cast<casadi_int>(d->h));
      casadi_ivp_clear_vjp(d);
      d->t_arg = d->t0;
      d->y_arg = d->y0;
      d->w = d->ry;
      return casadi_ivp_event_request(d, IVP_EVAL_JUMP_ADJ, IVP_EV_BACK_SEED);
    }
    d->ev_phase = IVP_EV_BACK_SEED;
    return casadi_ivp_event_adj(d);
  case IVP_EV_BACK_SEED:
    if (d->fail) return casadi_ivp_event_done(d, IVP_EVAL_FAILED);
    if (p->has_transition) {
      // Adjoint of the state before the event; quadratures pass through
      casadi_clear(d->ry, p->nrx + p->nrz);
      casadi_ivp_add_vjp(d, d->ry);
    }
    // The step that contains the event, ending there
    if (casadi_ivp_tape_pop(d) || d->record != IVP_RECORD_STEP) {
      return casadi_ivp_event_done(d, IVP_EVENT_FAILED);
    }
    d->record_pending = 1;
    if (d->t_ev != d->t1) {
      casadi_copy(d->ry, p->nr, d->seed);
      casadi_clear(d->ry, p->nr);
      casadi_ivp_interpolate(d, d->t_ev, d->y_out);
      casadi_ivp_seed_stages(d);
    }
    d->t = d->t_ev;
    return casadi_ivp_event_done(d, IVP_SUCCESS);
  default:
    return casadi_ivp_event_done(d, IVP_EVENT_FAILED);
  }
}

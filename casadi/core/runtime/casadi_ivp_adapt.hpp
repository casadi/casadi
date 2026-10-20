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

// Adaptive step size control of the one-step integrators, see casadi_ivp.hpp: the kernel
// requests it as task IVP_STEP_CONTROL, so only integrators with adaptive steps need this file

// C-REPLACE "static_cast<T1>" "(casadi_real) "
// C-REPLACE "fmin" "casadi_fmin"
// C-REPLACE "fmax" "casadi_fmax"
// C-REPLACE "fabs" "casadi_fabs"
// C-REPLACE "casadi_ivp_prob<T1>" "struct casadi_ivp_prob"
// C-REPLACE "casadi_ivp_data<T1>" "struct casadi_ivp_data"

// SYMBOL "ivp_step_size"
template<typename T1>
void casadi_ivp_step_size(casadi_ivp_data<T1>* d) {
  // Size h of a step from t; h = 0 if an initial estimate is needed first
  const casadi_ivp_prob<T1>* p = d->prob;
  d->last_clipped = 0;
  d->last_step = 0;
  if (d->h_next <= 0) {
    if (p->h_init <= 0) {
      d->h = 0;
      return;
    }
    d->h_next = p->h_init;
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

// SYMBOL "ivp_initial_step"
template<typename T1>
void casadi_ivp_initial_step(casadi_ivp_data<T1>* d, const T1* f0) {
  // Initial step size from |y0| and |f(y0)| over the error controlled states
  const casadi_ivp_prob<T1>* p = d->prob;
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
  d->h_next = fmax(h, p->h_min);
}

// SYMBOL "ivp_error_norm"
template<typename T1>
T1 casadi_ivp_error_norm(casadi_ivp_data<T1>* d) {
  // Weighted RMS norm of the error estimate over the error controlled entries
  const casadi_ivp_prob<T1>* p = d->prob;
  casadi_int k, i, n, offset[3], n_err[3];
  T1 ri, sum;
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
      ri = e[i] / (p->abstol + p->reltol * fmax(fabs(d->y0[i]), fabs(d->y[i])));
      sum += ri * ri;
      n++;
    }
  }
  return n > 0 ? sqrt(sum / static_cast<T1>(n)) : 0;
}

// SYMBOL "ivp_reject"
template<typename T1>
void casadi_ivp_reject(casadi_ivp_data<T1>* d, T1 fac) {
  // Restore the state, shrink the step by fac; status IVP_STEP_TOO_SMALL below h_min
  const casadi_ivp_prob<T1>* p = d->prob;
  casadi_copy(d->y0, p->ny, d->y);
  d->n_reject++;
  d->last_rejected = 1;
  d->h_next = d->h * fac;
  d->accepted = 0;
  if (d->h_next < p->h_min) d->status = IVP_STEP_TOO_SMALL;
}

// SYMBOL "ivp_error_test"
template<typename T1>
void casadi_ivp_error_test(casadi_ivp_data<T1>* d) {
  // Error test of the candidate y, PI controller for the next step size
  const casadi_ivp_prob<T1>* p = d->prob;
  T1 err, k, fac;
  err = casadi_ivp_error_norm(d);
  k = static_cast<T1>(p->order_err + 1);
  if (!(err <= 1)) {
    // Reject (also if not a number)
    fac = 0.25;
    if (casadi_ivp_isfinite(err)) fac = fmax(p->fac_min, p->safety * pow(err, -1. / k));
    casadi_ivp_reject(d, fac);
    return;
  }
  // PI controller
  err = fmax(err, 1e-10);
  fac = p->safety * pow(err, -0.7 / k) * pow(d->err_prev, 0.4 / k);
  fac = fmin(p->fac_max, fmax(p->fac_min, fac));
  if (d->last_rejected) fac = fmin(fac, 1.);
  d->err_prev = fmax(err, 1e-4);
  // A step clipped at the stop time does not shrink the proposal
  d->h_next = d->last_clipped ? fmax(d->h_next, d->h * fac) : d->h * fac;
  d->accepted = 1;
}

// SYMBOL "ivp_control"
template<typename T1>
void casadi_ivp_control(casadi_ivp_data<T1>* d) {
  // Carry out the step size control requested by the kernel (task IVP_STEP_CONTROL)
  switch (d->control) {
  case IVP_CONTROL_INIT:
    casadi_ivp_initial_step(d, d->U);
    casadi_ivp_step_size(d);
    break;
  case IVP_CONTROL_SIZE:
    casadi_ivp_step_size(d);
    break;
  case IVP_CONTROL_TEST:
    casadi_ivp_error_test(d);
    break;
  case IVP_CONTROL_FAIL:
    casadi_ivp_reject(d, 0.25);
    break;
  }
}

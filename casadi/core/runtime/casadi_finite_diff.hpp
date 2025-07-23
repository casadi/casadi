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

// C-REPLACE "std::numeric_limits<T1>::quiet_NaN()" "NAN"
// C-REPLACE "fmax" "casadi_fmax"

// SYMBOL "forward_diff"
template<typename T1>
void casadi_forward_diff(const T1** yk, T1* J, T1 h, casadi_int n_y) {
  // Local variables
  casadi_int i;
  const double *yf, *yc;
  T1 hinv;
  // Inverse of step size
  hinv = 1. / h;
  // Get stencil
  yc = yk[0];
  yf = yk[1];
  // Calculate FD approximation
  for (i = 0; i < n_y; ++i) J[i] = hinv * (yf[i] - yc[i]);
}

// SYMBOL "central_diff"
template<typename T1>
void casadi_central_diff(const T1** yk, T1* J, T1 h, casadi_int n_y) {
  // Local variables
  casadi_int i;
  const T1 *yf, *yc, *yb;
  T1 hinv;
  // Inverse of step size
  hinv = 1. / h;
  // Get stencil
  yb = yk[0];
  yc = yk[1];
  yf = yk[2];
  // Set u and stencils to zero (also supresses warnings)
  for (i = 0; i < n_y; ++i) {
    if (isfinite(yb[i])) {
      if (isfinite(yf[i])) {
        // Both forward and backward allowed
        J[i] = 0.5 * hinv * (yf[i] - yb[i]);
      } else {
        // Backward but not forward allowed
        J[i] = hinv * (yc[i] - yb[i]);
      }
    } else {
      if (isfinite(yf[i])) {
        // Forward but not backward allowed
        J[i] = hinv * (yf[i] - yc[i]);
      } else {
        // Neither forward nor backward possible
        J[i] = std::numeric_limits<T1>::quiet_NaN();
      }
    }
  }
}

// SYMBOL "central_diff_err"
template<typename T1>
T1 casadi_central_diff_err(const T1** yk, T1 h, casadi_int n_y, casadi_int i,
    T1 abstol, T1 reltol) {
  // Local variables
  const T1 *yf, *yc, *yb;
  T1 err_trunc, err_round;
  // Get stencil
  yb = yk[0];
  yc = yk[1];
  yf = yk[2];
  // Only consider points where both forward and backward allowed
  if (isfinite(yb[i]) && isfinite(yf[i])) {
    // Truncation error
    err_trunc = yf[i] - 2*yc[i] + yb[i];
    // Roundoff error
    err_round = reltol / h * fmax(fabs(yf[i] - yc[i]), fabs(yc[i] - yb[i])) + abstol;
    // Error quotient estimate
    return err_trunc / err_round;
  } else {
    // Cannot be calculated
    return std::numeric_limits<T1>::quiet_NaN();;
  }
}

// SYMBOL "smoothing_diff_weights"
template<typename T1>
T1 casadi_smoothing_diff_weights(casadi_int k, T1 yb, T1 yc, T1 yf, T1 *J) {
  // Calculate shifted finite difference approximation, weights
  if (k == 0) {
    // Backward shifted
    // 7.10 in Conte & Carl de Boor: Elementary Numerical Analysis (1972)
    // and 25.3.4 in Abramowitz and Stegun, Handbook of Mathematical Functions (1964)
    if (J) *J = 3 * yf - 4 * yc + yb;
    // Relative weight is 1
    return 1;
  } else if (k == 1) {
    // Central
    // We give this the relative weight 4 since if all weights are equal,
    // this would amount to a five point formula for the derivative
    // (yb2 - 8*yb + 8*yf - y_f2)/(12*h)
    // cf. 25.3.6 in Abramowitz and Stegun, Handbook of Mathematical Functions (1964)
    if (J) *J = yf - yb;
    // Relative weight is 4
    return 4;
  } else {
    // Forward shifted, cf. backward shifted above
    if (J) *J = -3 * yb + 4 * yc - yf;
    // Relative weight is 1
    return 1;
  }
}

// SYMBOL "smoothing_diff"
template<typename T1>
void casadi_smoothing_diff(const T1** yk, T1* J, T1 h, casadi_int n_y, T1 smoothing) {
  // Stencil
  T1 yb, yc, yf;
  // Local variables
  T1 Jk, wk, sw, ui, sm;
  casadi_int i, k;
  // Set stencils to zero (also supresses warnings)
  yf = yc = yb = 0;
  for (i = 0; i < n_y; ++i) {
    // Reset derivative estimate, sum of weights, error estimate
    J[i] = sw = ui = 0;
    // For backward shifted, central and forward shifted
    for (k = 0; k < 3; ++k) {
      // Get stencil
      yb = yk[k][i];
      yc = yk[k + 1][i];
      yf = yk[k + 2][i];
      // No contribuation if any value is infinite
      if (!isfinite(yb) || !isfinite(yc) || !isfinite(yf)) continue;
      // Calculate weights
      wk = casadi_smoothing_diff_weights(k, yb, yc, yf, &Jk);
      // Smoothness measure (second order derivative)
      sm = yf - 2*yc + yb;
      sm /= h*h;
      // Modify the weight according to smoothness
      wk /= sm*sm + smoothing;
      sw += wk;
      // Added weighted contribution to weight and error
      J[i] += wk * Jk;
    }
    // If sw is 0, no stencil worked
    if (sw == 0) {
      // Set component to 0, return -1
      J[i] = std::numeric_limits<T1>::quiet_NaN();
    } else {
      // Finalize estimate using the sum of weights and the step length
      J[i] /= 2*h*sw;
    }
  }
}

// SYMBOL "smoothing_diff_err"
template<typename T1>
T1 casadi_smoothing_diff_err(const T1** yk, T1 h, casadi_int n_y, casadi_int i,
    T1 abstol, T1 reltol, T1 smoothing) {
  // Stencil
  T1 yb, yc, yf;
  // Local variables
  T1 wk, sw, ui, err_trunc, err_round, sm;
  casadi_int k;
  // Set stencils to zero (also supresses warnings)
  yf = yc = yb = 0;
  // Reset derivative estimate, sum of weights, error estimate
  sw = ui = 0;
  // For backward shifted, central and forward shifted
  for (k = 0; k < 3; ++k) {
    // Get stencil
    yb = yk[k][i];
    yc = yk[k + 1][i];
    yf = yk[k + 2][i];
    // No contribuation if any value is infinite
    if (!isfinite(yb) || !isfinite(yc) || !isfinite(yf)) continue;
    // Calculate weights
    wk = casadi_smoothing_diff_weights(k, yb, yc, yf, static_cast<T1*>(0));
    // Truncation error
    err_trunc = yf - 2*yc + yb;
    // Roundoff error
    err_round = reltol/h*fmax(fabs(yf - yc), fabs(yc - yb)) + abstol;
    // We use the second order derivative as a smoothness measure
    sm = err_trunc/(h*h);
    // Modify the weight according to smoothness
    wk /= sm*sm + smoothing;
    sw += wk;
    // Added weighted contribution to weight and error
    ui += wk * fabs(err_trunc / err_round);
  }
  // If sw is 0, no stencil worked
  if (sw == 0) {
    // Cannot be calculated
    return std::numeric_limits<T1>::quiet_NaN();
  } else {
    // Finalize estimate using the sum of weights and the step length
    return ui / sw;
  }
}

// SYMBOL "finite_diff_mem"
template<typename T1>
struct casadi_finite_diff_mem {
  // Input precision
  T1 reltol;
  // Output precision
  T1 abstol;
  // Smoothness parameter
  // Smaller epsilon: More discontinuity rejecting
  // Larger epsilon: More accurate (higher order) if smooth
  T1 smoothing;
};

// C-REPLACE "casadi_finite_diff_mem<T1>" "struct casadi_finite_diff_mem"

// SYMBOL "forward_diff_old"
template<typename T1>
T1 casadi_forward_diff_old(T1** yk, T1* y0, T1* J,
    T1 h, casadi_int n_y, const casadi_finite_diff_mem<T1>* m) {
  casadi_int i;
  for (i=0; i<n_y; ++i) {
    J[i] = (yk[0][i]-y0[i])/h;
  }
  return -1;
}

// SYMBOL "central_diff_old"
template<typename T1>
T1 casadi_central_diff_old(T1** yk, T1* y0, T1* J,
    T1 h, casadi_int n_y, const casadi_finite_diff_mem<T1>* m) {
  // Return value
  T1 u;
  // Stencil
  T1 yf, yc, yb;
  // Local variables
  T1 err_trunc, err_round;
  casadi_int i;
  // Set u and stencils to zero (also supresses warnings)
  yf = yc = yb = u = 0;
  for (i=0; i<n_y; ++i) {
    // Copy to local variables, return -1 if invalid entry
    if (!isfinite((yf=yk[1][i])) || !isfinite((yc=y0[i])) || !isfinite((yb=yk[0][i]))) {
      J[i] = std::numeric_limits<T1>::quiet_NaN();
      u = -1;
      continue;
    }
    // Central difference approximation
    J[i] = (yf - yb)/(2*h);
    // Truncation error
    err_trunc = yf - 2*yc + yb;
    // Roundoff error
    err_round = m->reltol/h*fmax(fabs(yf-yc), fabs(yc-yb)) + m->abstol;
    // Update error estimate
    if (u>=0) u = fmax(u, fabs(err_trunc/err_round));
  }
  return u;
}

// SYMBOL "smoothing_diff_old"
template<typename T1>
T1 casadi_smoothing_diff_old(T1** yk, T1* y0, T1* J,
    T1 h, casadi_int n_y, const casadi_finite_diff_mem<T1>* m) {
  // Return value
  T1 u;
  // Stencil
  T1 yb, yc, yf;
  // Local variables
  T1 Jk, wk, sw, ui, err_trunc, err_round, sm;
  casadi_int i, k;
  // Set u and stencils to zero (also supresses warnings)
  yf = yc = yb = u = 0;
  for (i=0; i<n_y; ++i) {
    // Reset derivative estimate, sum of weights, error estimate
    J[i] = sw = ui = 0;
    // For backward shifted, central and forward shifted
    for (k=0; k<3; ++k) {
      // Calculate shifted finite difference approximation
      if (k==0) {
        // Backward shifted
        // 7.10 in Conte & Carl de Boor: Elementary Numerical Analysis (1972)
        // and 25.3.4 in Abramowitz and Stegun, Handbook of Mathematical Functions (1964)
        if (!isfinite((yc=yk[0][i]))) continue;
        if (!isfinite((yb=yk[1][i]))) continue;
        yf = y0[i];
        Jk = 3*yf - 4*yc + yb;
        wk = 1;
      } else if (k==1) {
        // Central
        // We give this the "nominal weight" 4 since if all weights are equal,
        // this would amount to a five point formula for the derivative
        // (yb2 - 8*yb + 8*yf - y_f2)/(12*h)
        // cf. 25.3.6 in Abramowitz and Stegun, Handbook of Mathematical Functions (1964)
        if (!isfinite((yf=yk[2][i]))) continue;
        if (!isfinite((yb=yc))) continue;
        yc = y0[i];
        Jk = yf-yb;
        wk = 4;
      } else {
        // Forward shifted
        if (!isfinite((yc=yf))) continue;
        if (!isfinite((yf=yk[3][i]))) continue;
        yb = y0[i];
        Jk = -3*yb + 4*yc - yf;
        wk = 1;
      }
      // Truncation error
      err_trunc = yf - 2*yc + yb;
      // Roundoff error
      err_round = m->reltol/h*fmax(fabs(yf-yc), fabs(yc-yb)) + m->abstol;
      // We use the second order derivative as a smoothness measure
      sm = err_trunc/(h*h);
      // Modify the weight according to smoothness
      wk /= sm*sm + m->smoothing;
      sw += wk;
      // Added weighted contribution to weight and error
      J[i] += wk * Jk;
      ui += wk * fabs(err_trunc/err_round);
    }
    // If sw is 0, no stencil worked
    if (sw==0) {
      // Set component to 0, return -1
      J[i] = std::numeric_limits<T1>::quiet_NaN();
      u = -1;
    } else {
      // Finalize estimate using the sum of weights and the step length
      J[i] /= 2*h*sw;
      if (u>=0) u = fmax(u, ui/sw);
    }
  }
  return u;
}

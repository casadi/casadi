// NOLINT(legal/copyright)
template<typename T1>
struct casadi_central_diff_mem {
  // Dimensions
  int n_x, n_r;
  // Perturbation sizes
  T1 *h;
  // Perturbation size range
  T1 h_max, eps, eps1;
  // Ratio of roundoff error to truncation error
  T1 u_min, u_aim, u_max;
  // (Current) differentiable inputs
  T1 *x, *x0;
  // (Current) differentiable outputs
  T1 *r, *r0;
  // Jacobian
  T1* J;
  // Control
  int status, i;
};

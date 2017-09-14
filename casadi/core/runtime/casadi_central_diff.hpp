// NOLINT(legal/copyright)
// SYMBOL "central_diff_mem"
template<typename T1>
struct casadi_central_diff_mem {
  // Dimension
  int n_r;
  // Perturbation size
  T1 h;
  // Perturbation size range
  T1 h_max, eps, eps1;
  // Target ratio of roundoff error to truncation error
  T1 u_aim;
  // (Current) differentiable inputs
  T1 x, x0;
  // (Current) differentiable outputs
  T1 *r, *r0;
  // Jacobian
  T1* J;
  // Control
  int status;
};

// C-REPLACE "casadi_central_diff_mem<T1>" "struct casadi_central_diff_mem"
// SYMBOL "central_diff"
template<typename T1>
int casadi_central_diff(casadi_central_diff_mem<T1>* m) {
  // Control loop
  const int wait_for_pos_pert=1, wait_for_neg_pert=2;
  switch (m->status) {
  case 0:
    // Backup r
    m->x0 = m->x;
    casadi_copy(m->r, m->n_r, m->r0);
    // Perturb x in positive direction and call function
    m->x += m->h/2;
    m->status = wait_for_pos_pert;
    return 1;
  case wait_for_pos_pert:
    casadi_copy(m->r, m->n_r, m->J);
    // Perturb x in negative direction and call function
    m->x = m->x0 - m->h/2;
    m->status = wait_for_neg_pert;
    return 1;
  case wait_for_neg_pert:
    m->x = m->x0;
    // Calculate finite difference approximation
    casadi_axpy(m->n_r, -1., m->r, m->J);
    casadi_scal(m->n_r, 1/m->h, m->J);
    // Successful return, restore r
    casadi_copy(m->r0, m->n_r, m->r);
    m->status = 0;
  }
  return 0;
}

// NOLINT(legal/copyright)
#define central_diff_mem CASADI_PREFIX(central_diff_mem)<T1>
template<typename T1>
int CASADI_PREFIX(central_diff)(central_diff_mem* m) {
  // Control loop
  const int wait_for_pos_pert=1, wait_for_neg_pert=2;
  switch (m->status) {
    case 0: break;
    case wait_for_pos_pert: goto back_from_pos_pert;
    case wait_for_neg_pert: goto back_from_neg_pert;
    default: return 0;
  }
  // Backup x and f
  CASADI_PREFIX(copy)(m->x, m->n_x, m->x0);
  CASADI_PREFIX(copy)(m->f, m->n_f, m->f0);
  // Loop over directions
  for (m->i=0; m->i<m->n_x; ++m->i) {
    // Perturb x in positive direction and call function
    m->x[m->i] += m->h/2;
    m->status = wait_for_pos_pert;
    return 1;
    back_from_pos_pert:
    CASADI_PREFIX(copy)(m->f, m->n_f, m->J + m->i*m->n_f);
    // Perturb x in negative direction and call function
    m->x[m->i] = m->x0[m->i] - m->h/2;
    m->status = wait_for_neg_pert;
    return 1;
    back_from_neg_pert:
    m->x[m->i] = m->x0[m->i];
    // Calculate finite difference approximation
    CASADI_PREFIX(axpy)(m->n_f, -1., m->f, m->J + m->i*m->n_f);
    CASADI_PREFIX(scal)(m->n_f, 1/m->h, m->J + m->i*m->n_f);
  }
  // Successful return, restore f
  CASADI_PREFIX(copy)(m->f0, m->n_f, m->f);
  m->status = 0;
  return 0;
}
#undef central_diff_mem

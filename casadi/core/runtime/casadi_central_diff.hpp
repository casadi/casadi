// NOLINT(legal/copyright)
#define central_diff_mem CASADI_PREFIX(central_diff_mem)<T1>
template<typename T1>
int CASADI_PREFIX(central_diff)(central_diff_mem* m) {
  // Define cases
  const int INIT=0, POS_PERT=1, NEG_PERT=2, CALC_FD=3;
  // Call user function?
  const int DONE=0, CALL=1;
  // Control loop
  switch (m->next) {
    case INIT:
    // Backup x and f
    CASADI_PREFIX(copy)(m->x, m->n_x, m->x0);
    CASADI_PREFIX(copy)(m->f, m->n_f, m->f0);
    // Reset direction counter
    m->i = 0;
    m->next = POS_PERT;
    // No break or return
    case POS_PERT:
    // Perturb x, positive direction
    m->x[m->i] += m->h/2;
    m->next = NEG_PERT;
    return CALL;
    case NEG_PERT:
    // Save result, perturb in negative direction
    CASADI_PREFIX(copy)(m->f, m->n_f, m->J + m->i*m->n_f);
    m->x[m->i] = m->x0[m->i] - m->h/2;
    m->next = CALC_FD;
    return CALL;
    case CALC_FD:
    // Reset x
    m->x[m->i] = m->x0[m->i];
    // Calculate finite difference approximation
    CASADI_PREFIX(axpy)(m->n_f, -1., m->f, m->J + m->i*m->n_f);
    CASADI_PREFIX(scal)(m->n_f, 1/m->h, m->J + m->i*m->n_f);
    // Procede to next direction
    m->i++;
    if (m->i < m->n_x) {
      // Continue to the next direction
      m->next = POS_PERT;
      return CALL;
    }
    // Restore f
    CASADI_PREFIX(copy)(m->x0, m->n_x, m->x);
  }
  return DONE;
}
#undef central_diff_mem

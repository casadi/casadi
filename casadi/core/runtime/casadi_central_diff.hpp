// NOLINT(legal/copyright)
template<typename T1>
bool CASADI_PREFIX(central_diff)(CASADI_PREFIX(central_diff_mem)<T1>* m) {
  bool ret = true;
  switch (m->n_calls) {
    case 0:
    // Backup x and f
    CASADI_PREFIX(copy)(m->x, m->n_x, m->x0);
    CASADI_PREFIX(copy)(m->f, m->n_f, m->f0);
    // Perturb x, positive direction
    CASADI_PREFIX(axpy)(m->n_x, m->h/2, m->v, m->x);
    break;
    case 1:
    // Save result, perturb in negative direction
    CASADI_PREFIX(copy)(m->f, m->n_f, m->Jv);
    CASADI_PREFIX(copy)(m->x0, m->n_x, m->x);
    CASADI_PREFIX(axpy)(m->n_x, -m->h/2, m->v, m->x);
    break;
    case 2:
    // Calculate finite difference approximation
    CASADI_PREFIX(axpy)(m->n_f, -1., m->f, m->Jv);
    CASADI_PREFIX(scal)(m->n_f, 1/m->h, m->Jv);
    // Restore x and f
    CASADI_PREFIX(copy)(m->x0, m->n_x, m->x);
    CASADI_PREFIX(copy)(m->f0, m->n_f, m->f);
    ret = false;
  }
  // Increase function call counter
  if (ret) m->n_calls++;
  return ret;
}

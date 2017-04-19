template<typename real_t>
real_t CASADI_PREFIX(norm_1)(int n, const real_t* x) {
  real_t ret = 0;
  int i;
  if (x) {
    for (i=0; i<n; ++i) ret += fabs(*x++);
  }
  return ret;
}

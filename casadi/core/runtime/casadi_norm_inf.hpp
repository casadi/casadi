template<typename real_t>
real_t CASADI_PREFIX(norm_inf)(int n, const real_t* x) {
  real_t ret = 0;
  int i;
  for (i=0; i<n; ++i) ret = fmax(ret, fabs(*x++));
  return ret;
}

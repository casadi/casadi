template<typename real_t>
real_t CASADI_PREFIX(dot)(int n, const real_t* x, const real_t* y) {
  real_t r = 0;
  int i;
  for (i=0; i<n; ++i) r += *x++ * *y++;
  return r;
}

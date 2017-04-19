template<typename real_t>
real_t CASADI_PREFIX(norm_2)(int n, const real_t* x) {
  return sqrt(CASADI_PREFIX(dot)(n, x, x));
}

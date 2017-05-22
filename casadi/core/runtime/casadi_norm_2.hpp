// NOLINT(legal/copyright)
template<typename T1>
T1 CASADI_PREFIX(norm_2)(int n, const T1* x) {
  return sqrt(CASADI_PREFIX(dot)(n, x, x));
}

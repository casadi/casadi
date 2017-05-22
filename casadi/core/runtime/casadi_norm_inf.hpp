// NOLINT(legal/copyright)
template<typename T1>
T1 CASADI_PREFIX(norm_inf)(int n, const T1* x) {
  T1 ret = 0;
  int i;
  for (i=0; i<n; ++i) ret = fmax(ret, fabs(*x++));
  return ret;
}

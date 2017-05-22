// NOLINT(legal/copyright)
template<typename T1>
T1 CASADI_PREFIX(norm_1)(int n, const T1* x) {
  T1 ret = 0;
  int i;
  if (x) {
    for (i=0; i<n; ++i) ret += fabs(*x++);
  }
  return ret;
}

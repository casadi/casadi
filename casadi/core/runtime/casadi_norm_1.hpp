// NOLINT(legal/copyright)
// SYMBOL "norm_1"
template<typename T1>
T1 casadi_norm_1(int n, const T1* x) {
  T1 ret = 0;
  int i;
  if (x) {
    for (i=0; i<n; ++i) ret += fabs(*x++);
  }
  return ret;
}

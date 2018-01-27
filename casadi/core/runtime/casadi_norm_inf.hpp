// NOLINT(legal/copyright)
// SYMBOL "norm_inf"
template<typename T1>
T1 casadi_norm_inf(casadi_int n, const T1* x) {
  T1 ret = 0;
  casadi_int i;
  for (i=0; i<n; ++i) ret = fmax(ret, fabs(*x++));
  return ret;
}

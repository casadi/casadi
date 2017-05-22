// NOLINT(legal/copyright)
template<typename T1>
T1 CASADI_PREFIX(max_viol)(int n, const T1* x, const T1* lb, const T1* ub) {
  T1 r = 0;
  const T1 zero = 0;
  int i;
  for (i=0; i<n; ++i) {
    T1 x_i = x ? *x++ : zero;
    T1 lb_i = lb ? *lb++ : zero;
    T1 ub_i = ub ? *ub++ : zero;
    r = fmax(r, fmax(x_i-ub_i, zero));
    r = fmax(r, fmax(lb_i-x_i, zero));
  }
  return r;
}

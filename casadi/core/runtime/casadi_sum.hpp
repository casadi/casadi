// NOLINT(legal/copyright)
// SYMBOL "sum"
template<typename T1>
T1 casadi_sum(const T1* x, casadi_int n) {
  casadi_int i;
  T1 r = 0;

  r = 0;
  for (i=0; i<n; ++i) {
    r += x[i];
  }
  return r;
}

// NOLINT(legal/copyright)
// SYMBOL "axpy"
template<typename T1>
void casadi_axpy(casadi_int n, T1 alpha, const T1* x, T1* y) {
  if (!x || !y) return;
  casadi_int i;
  for (i=0; i<n; ++i) *y++ += alpha**x++;
}

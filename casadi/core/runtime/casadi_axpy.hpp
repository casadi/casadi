// NOLINT(legal/copyright)
template<typename T1>
void CASADI_PREFIX(axpy)(int n, T1 alpha, const T1* x, T1* y) {
  int i;
  for (i=0; i<n; ++i) *y++ += alpha**x++;
}

// NOLINT(legal/copyright)
template<typename T1>
void CASADI_PREFIX(scal)(int n, T1 alpha, T1* x) {
  int i;
  for (i=0; i<n; ++i) *x++ *= alpha;
}

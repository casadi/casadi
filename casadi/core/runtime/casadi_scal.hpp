// NOLINT(legal/copyright)
// SYMBOL "scal"
template<typename T1>
void casadi_scal(int n, T1 alpha, T1* x) {
  if (!x) return;
  int i;
  for (i=0; i<n; ++i) *x++ *= alpha;
}

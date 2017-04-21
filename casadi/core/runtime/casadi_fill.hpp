// NOLINT(legal/copyright)
template<typename T1>
void CASADI_PREFIX(fill)(T1* x, int n, T1 alpha) {
  int i;
  if (x) {
    for (i=0; i<n; ++i) *x++ = alpha;
  }
}

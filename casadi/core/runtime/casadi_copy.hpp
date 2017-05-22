// NOLINT(legal/copyright)
template<typename T1>
void CASADI_PREFIX(copy)(const T1* x, int n, T1* y) {
  int i;
  if (y) {
    if (x) {
      for (i=0; i<n; ++i) *y++ = *x++;
    } else {
      for (i=0; i<n; ++i) *y++ = 0.;
    }
  }
}

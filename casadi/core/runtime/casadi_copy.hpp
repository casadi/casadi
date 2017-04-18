template<typename real_t>
void CASADI_PREFIX(copy)(const real_t* x, int n, real_t* y) {
  int i;
  if (y) {
    if (x) {
      for (i=0; i<n; ++i) *y++ = *x++;
    } else {
      for (i=0; i<n; ++i) *y++ = 0.;
    }
  }
}

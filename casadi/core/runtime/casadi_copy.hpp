// NOLINT(legal/copyright)
// SYMBOL "copy"
template<typename T1>
void casadi_copy(const T1* x, casadi_int n, T1* y) {
  casadi_int i;
  if (y) {
    if (x) {
      //for (i=0; i<n; ++i) *y++ = *x++;
      memcpy(y, x, n*sizeof(T1));
    } else {
      memset(y, 0, n*sizeof(T1));
      //for (i=0; i<n; ++i) *y++ = 0.;
    }
  }
}

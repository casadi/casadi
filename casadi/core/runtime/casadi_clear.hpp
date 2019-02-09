// NOLINT(legal/copyright)
// SYMBOL "clear"
template<typename T1>
void casadi_clear(T1* x, casadi_int n) {
  casadi_int i;
  if (x) {
    for (i=0; i<n; ++i) *x++ = 0;
  }
}

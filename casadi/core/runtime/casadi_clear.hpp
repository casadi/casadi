// NOLINT(legal/copyright)
// SYMBOL "clear"
template<typename T1>
void casadi_clear(T1* x, casadi_int n) {
  casadi_int i;
  if (x) {
    memset(x, 0, n*sizeof(T1));
    //for (i=0; i<n; ++i) *x++ = 0;
  }
}

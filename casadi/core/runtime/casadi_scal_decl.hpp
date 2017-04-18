template<typename real_t>
void CASADI_PREFIX(scal)(int n, real_t alpha, real_t* x) {
  int i;
  for (i=0; i<n; ++i) *x++ *= alpha;
}

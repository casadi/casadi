template<typename real_t>
void CASADI_PREFIX(axpy)(int n, real_t alpha, const real_t* x, real_t* y) {
  int i;
  for (i=0; i<n; ++i) *y++ += alpha**x++;
}

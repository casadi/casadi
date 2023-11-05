// NOLINT(legal/copyright)
// SYMBOL "axpy"
template<typename T1>
void casadi_axpy(casadi_int n, T1 alpha, const T1* x, T1* y) {
  casadi_int i;
  if (!x || !y) return;
#if defined(CASADI_WITH_CBLAS_DOUBLE)
  cblas_daxpy(n, alpha, x, 1, y, 1);
#elif defined(CASADI_WITH_CBLAS_SINGLE)
  cblas_saxpy(n, alpha, x, 1, y, 1);
#else
  for (i=0; i<n; ++i) *y++ += alpha**x++;
#endif
}

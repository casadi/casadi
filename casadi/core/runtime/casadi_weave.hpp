// NOLINT(legal/copyright)
// SYMBOL "weave"
template<typename T1>
void casadi_weave(const T1* x, casadi_int nnz, casadi_int m, casadi_int n, T1* y) {
  casadi_int i, j, stride;
  const T1 * src;

  if (y) {
    if (x) {
      stride = nnz*m;
      for (j=0;j<m;++j) {
        src = x;
        for (i=0;i<n;++i) {
          casadi_copy(src, nnz, y);
          y += nnz;
          src += stride;
        }
        x += nnz;
      }
    } else {
      casadi_clear(y, nnz*m*n);
    }
  }
}

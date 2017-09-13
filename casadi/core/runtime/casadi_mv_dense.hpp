// NOLINT(legal/copyright)
// SYMBOL "mv_dense"
template<typename T1>
void casadi_mv_dense(const T1* x, int nrow_x, int ncol_x, const T1* y, T1* z, int tr) {
  if (!x || !y || !z) return;
  int i, j;
  if (tr) {
    for (i=0; i<ncol_x; ++i) {
      for (j=0; j<nrow_x; ++j) {
        z[i] += *x++ * y[j];
      }
    }
  } else {
    for (i=0; i<ncol_x; ++i) {
      for (j=0; j<nrow_x; ++j) {
        z[j] += *x++ * y[i];
      }
    }
  }
}

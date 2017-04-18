template<typename real1_t, typename real2_t>
void CASADI_PREFIX(sparsify)(const real1_t* x, real2_t* y, const int* sp_y, int tr) {
  int nrow_y = sp_y[0], ncol_y = sp_y[1];
  const int *colind_y = sp_y+2, *row_y = sp_y+ncol_y+3;
  int i, el;
  if (tr) {
    for (i=0; i<ncol_y; ++i) {
      for (el=colind_y[i]; el!=colind_y[i+1]; ++el) {
        *y++ = CASADI_CAST(real2_t, x[i + row_y[el]*ncol_y]);
      }
    }
  } else {
    for (i=0; i<ncol_y; ++i) {
      for (el=colind_y[i]; el!=colind_y[i+1]; ++el) {
        *y++ = CASADI_CAST(real2_t, x[row_y[el] + i*nrow_y]);
      }
    }
  }
}

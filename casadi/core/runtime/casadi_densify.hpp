template<typename real1_t, typename real2_t>
void CASADI_PREFIX(densify)(const real1_t* x, const int* sp_x, real2_t* y, int tr) {
  /* Quick return - output ignored */
  if (!y) return;
  int nrow_x = sp_x[0], ncol_x = sp_x[1];
  const int *colind_x = sp_x+2, *row_x = sp_x+ncol_x+3;
  /* Zero out return value */
  CASADI_PREFIX(fill)(y, nrow_x*ncol_x, CASADI_CAST(real2_t, 0));
  /* Quick return - input is zero */
  if (!x) return;
  /* Copy nonzeros */
  int i, el;
  if (tr) {
    for (i=0; i<ncol_x; ++i) {
      for (el=colind_x[i]; el!=colind_x[i+1]; ++el) {
        y[i + row_x[el]*ncol_x] = CASADI_CAST(real2_t, *x++);
      }
    }
  } else {
    for (i=0; i<ncol_x; ++i) {
      for (el=colind_x[i]; el!=colind_x[i+1]; ++el) {
        y[row_x[el] + i*nrow_x] = CASADI_CAST(real2_t, *x++);
      }
    }
  }
}

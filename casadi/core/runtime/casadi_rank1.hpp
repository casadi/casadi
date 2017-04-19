template<typename real_t>
void CASADI_PREFIX(rank1)(real_t* A, const int* sp_A, real_t alpha, const real_t* x, const real_t* y) {
  /* Get sparsities */
  int ncol_A = sp_A[1];
  const int *colind_A = sp_A+2, *row_A = sp_A + 2 + ncol_A+1;

  /* Loop over the columns of A */
  int cc, rr, el;
  for (cc=0; cc<ncol_A; ++cc) {
    /* Loop over the nonzeros of A */
    for (el=colind_A[cc]; el<colind_A[cc+1]; ++el) {
      /* Get row */
      rr = row_A[el];

      /* Add the multiple */
      A[el] += alpha*x[rr]*y[cc];
    }
  }
}

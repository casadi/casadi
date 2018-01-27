// NOLINT(legal/copyright)
// SYMBOL "ldl"
// Calculate the nonzeros of the transposed L factor (strictly lower entries only)
// as well as D for an LDL^T factorization
// len[w] >= n
template<typename T1>
void casadi_ldl_new(const int* sp_a, const int* sp_lt,
                    const T1* a, T1* lt, T1* d, T1* w) {
  // Extract sparsities
  int n=sp_lt[1];
  const int *lt_colind=sp_lt+2, *lt_row=sp_lt+2+n+1;
  const int *a_colind=sp_a+2, *a_row=sp_a+2+n+1;
  // Local variables
  int r, c, k, k2;
  // Clear w
  for (r=0; r<n; ++r) w[r] = 0;
  // Sparse copy of A to L and D
  for (c=0; c<n; ++c) {
    for (k=a_colind[c]; k<a_colind[c+1]; ++k) w[a_row[k]] = a[k];
    for (k=lt_colind[c]; k<lt_colind[c+1]; ++k) lt[k] = w[lt_row[k]];
    d[c] = w[c];
    for (k=a_colind[c]; k<a_colind[c+1]; ++k) w[a_row[k]] = 0;
  }
  // Loop over columns of L
  for (c=0; c<n; ++c) {
    for (k=lt_colind[c]; k<lt_colind[c+1]; ++k) {
      r = lt_row[k];
      // Calculate l(r,c) with r<c
      for (k2=lt_colind[r]; k2<lt_colind[r+1]; ++k2) {
        lt[k] -= lt[k2] * w[lt_row[k2]];
      }
      w[r] = lt[k];
      lt[k] /= d[r];
      // Update d(c)
      d[c] -= w[r]*lt[k];
    }
    // Clear w
    for (k=lt_colind[c]; k<lt_colind[c+1]; ++k) w[lt_row[k]] = 0;
  }
}

// SYMBOL "ldl_trs"
// Solve for (I+R) with R an optionally transposed strictly upper triangular matrix.
template<typename T1>
void casadi_ldl_trs_new(const int* sp_r, const T1* nz_r, T1* x, int tr) {
  // Extract sparsity
  int ncol=sp_r[1];
  const int *colind=sp_r+2, *row=sp_r+2+ncol+1;
  // Local variables
  int c, k;
  if (tr) {
    // Forward substitution
    for (c=0; c<ncol; ++c) {
      for (k=colind[c]; k<colind[c+1]; ++k) {
        x[c] -= nz_r[k]*x[row[k]];
      }
    }
  } else {
    // Backward substitution
    for (c=ncol-1; c>=0; --c) {
      for (k=colind[c+1]-1; k>=colind[c]; --k) {
        x[row[k]] -= nz_r[k]*x[c];
      }
    }
  }
}

// SYMBOL "ldl_solve"
// Linear solve using an LDL factorized linear system
template<typename T1>
void casadi_ldl_solve_new(T1* x, int nrhs, const int* sp_lt, const T1* lt,
                      const T1* d) {
  int n = sp_lt[1];
  int i, k;
  for (k=0; k<nrhs; ++k) {
    //      LDL'x = b <=> x = L'\D\L\b
    //  Solve for L
    casadi_ldl_trs_new(sp_lt, lt, x, 1);
    // Divide by D
    for (i=0; i<n; ++i) x[i] /= d[i];
    // Solve for L'
    casadi_ldl_trs_new(sp_lt, lt, x, 0);
    // Next rhs
    x += n;
  }
}

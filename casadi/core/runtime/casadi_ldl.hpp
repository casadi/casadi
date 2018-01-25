// NOLINT(legal/copyright)
// SYMBOL "ldl"
// Calculate the nonzeros of the L factor (strictly lower entries only)
// as well as D for an LDL^T factorization
// Ref: User Guide for LDL by Tim Davis
// len[iw] >= 2*n
// len[w] >= n
/* Only declaration, no implementation due to license restrictions.
  Cf. CasADi issue #2158
 */
template<typename T1>
void casadi_ldl_new(const int* sp_a, const int* parent, const int* sp_lt,
                 const T1* a, T1* lt, T1* d, int *iw, T1* w) {
  // Extract sparsities
  int n=sp_lt[1];
  const int *lt_colind=sp_lt+2, *lt_row=sp_lt+2+n+1;
  const int *a_colind=sp_a+2, *a_row=sp_a+2+n+1;
  // Local variables
  int r, c, k;
  // Clear w
  for (r=0; r<n; ++r) w[r] = 0;
  // Loop over columns of A
  for (c=0; c<n; ++c) {
    // Copy column of a to w
    for (k=a_colind[c]; k<a_colind[c+1]; ++k) w[a_row[k]] = a[k];
    // Calculate d using d_cc = a_cc - sum_{k<c} d_k * l_ck^2
    d[c] = w[c];
    for (k=lt_colind[c]; k<lt_colind[c+1]; ++k) {
      r = lt_row[k];
      d[c] -= d[r]*lt[k]*lt[k]
    }
    // Calculate l using l_rc = (a_rc - sum_{k<c} d_k * l_rk * l_kc) / d_c
    for (k=lt_colind[c]; k<lt_colind[c+1]; ++k) {
      r = lt_row[k];
      // FIXME: Looping over r will destroy algorithm complexity, need to
      // keep track of visited nodes
    }
    // Clear w
    for (k=a_colind[c]; k<a_colind[c+1]; ++k) w[a_row[k]] = 0;
  }
}

// SYMBOL "ldl_trs"
// Solve for (I+L) with L an optionally transposed strictly lower triangular matrix.
template<typename T1>
void casadi_ldl_trs(const int* sp_l, const T1* nz_l, T1* x, int tr) {
  // Extract sparsity
  int ncol=sp_l[1];
  const int *colind=sp_l+2, *row=sp_l+2+ncol+1;
  // Local variables
  int c, k;
  if (tr) {
    // Backward substitution
    for (c=ncol-1; c>=0; --c) {
      for (k=colind[c]; k<colind[c+1]; ++k) {
        x[c] -= nz_l[k]*x[row[k]];
      }
    }
  } else {
    // Forward substitution
    for (c=0; c<ncol; ++c) {
      for (k=colind[c]; k<colind[c+1]; ++k) {
        x[row[k]] -= nz_l[k]*x[c];
      }
    }
  }
}

// SYMBOL "ldl_solve"
// Linear solve using an LDL factorized linear system
template<typename T1>
void casadi_ldl_solve(T1* x, int nrhs, const int* sp_l, const T1* l,
                      const T1* d) {
  int n = sp_l[1];
  int i, k;
  for (k=0; k<nrhs; ++k) {
    //      LDL'x = b <=> x = L'\D\L\b
    //  Solve for L
    casadi_ldl_trs(sp_l, l, x, 0);
    // Divide by D
    for (i=0; i<n; ++i) x[i] /= d[i];
    // Solve for L'
    casadi_ldl_trs(sp_l, l, x, 1);
    // Next rhs
    x += n;
  }
}

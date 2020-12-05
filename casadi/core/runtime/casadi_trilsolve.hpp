// NOLINT(legal/copyright)

// SYMBOL "trilsolve"
template<typename T1>
void casadi_trilsolve(const casadi_int* sp_a, const T1* nz_a, T1* x, int tr, int unity,
    casadi_int nrhs) {
  // Local variables
  casadi_int nrow, ncol, r, c, k, rhs;
  const casadi_int *colind, *row;
  // Extract sparsity
  nrow = sp_a[0];
  ncol = sp_a[1];
  colind = sp_a + 2;
  row = colind + ncol + 1;
  // For all right hand sides
  for (rhs = 0; rhs < nrhs; ++rhs) {
    if (unity) {
      if (tr) {
        // Backward substitution
        for (c = ncol; c-- > 0; ) {
          for (k = colind[c + 1]; k-- > colind[c]; ) {
            x[c] += nz_a[k] * x[row[k]];
          }
        }
      } else {
        // Forward substitution
        for (c = 0; c < ncol; ++c) {
          for (k = colind[c]; k < colind[c+1]; ++k) {
            x[row[k]] += nz_a[k] * x[c];
          }
        }
      }
    } else {
      if (tr) {
        // Backward substitution
        for (c = ncol; c-- > 0; ) {
          for (k = colind[c + 1]; k-- > colind[c]; ) {
            r = row[k];
            if (r == c) {
              x[c] /= nz_a[k];
            } else {
              x[c] -= nz_a[k] * x[r];
            }
          }
        }
      } else {
        // Forward substitution
        for (c = 0; c < ncol; ++c) {
          for (k = colind[c]; k < colind[c+1]; ++k) {
            r = row[k];
            if (r == c) {
              x[r] /= nz_a[k];
            } else {
              x[r] -= nz_a[k] * x[c];
            }
          }
        }
      }
    }
    // Next right-hand-side
    x += nrow;
  }
}

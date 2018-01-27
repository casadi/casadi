// NOLINT(legal/copyright)
// SYMBOL "house"
// Householder reflection
// Ref: Chapter 5, Direct Methods for Sparse Linear Systems by Tim Davis
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wconversion"
template<typename T1>
T1 casadi_house(T1* x, T1* beta, casadi_int n) {
  // Local variable
  casadi_int i;
  // Calculate norm
  T1 x0 = x[0]; // Save x0 (overwritten below)
  T1 sigma=0;
  for (i=1; i<n; ++i) sigma += x[i]*x[i];
  T1 s = sqrt(x0*x0 + sigma); // s = norm(x)
  // Calculate consistently with symbolic datatypes (SXElem)
  T1 sigma_is_zero = sigma==0;
  T1 x0_nonpos = x0<=0;
  x[0] = if_else(sigma_is_zero, 1,
                 if_else(x0_nonpos, x0-s, -sigma/(x0+s)));
  *beta = if_else(sigma_is_zero, 2*x0_nonpos, -1/(s*x[0]));
  return s;
}

// SYMBOL "qr"
// Numeric QR factorization
// Ref: Chapter 5, Direct Methods for Sparse Linear Systems by Tim Davis
// Note: nrow <= nrow_ext <= nrow+ncol
// len[iw] = nrow_ext + ncol
// len[x] = nrow_ext
// sp_v = [nrow_ext, ncol, 0, 0, ...] len[3 + ncol + nnz_v]
// len[v] nnz_v
// sp_r = [nrow_ext, ncol, 0, 0, ...] len[3 + ncol + nnz_r]
// len[r] nnz_r
// len[beta] ncol
/* Only declaration, no implementation due to license restrictions.
  Cf. CasADi issue #2158
 */
template<typename T1>
void casadi_qr(const casadi_int* sp_a, const T1* nz_a, casadi_int* iw, T1* x,
               const casadi_int* sp_v, T1* nz_v, const casadi_int* sp_r, T1* nz_r, T1* beta,
               const casadi_int* leftmost, const casadi_int* parent, const casadi_int* pinv);

// SYMBOL "qr_mv"
// Multiply QR Q matrix from the right with a vector, with Q represented
// by the Householder vectors V and beta
// x = Q*x or x = Q'*x
// with Q = (I-beta(1)*v(:,1)*v(:,1)')*...*(I-beta(n)*v(:,n)*v(:,n)')
// len[x] >= nrow_ext
template<typename T1>
void casadi_qr_mv(const casadi_int* sp_v, const T1* v, const T1* beta, T1* x,
                  casadi_int tr) {
  // Extract sparsity
  casadi_int ncol=sp_v[1];
  const casadi_int *colind=sp_v+2, *row=sp_v+2+ncol+1;
  // Local variables
  casadi_int c, c1, k;
  T1 alpha;
  // Loop over vectors
  for (c1=0; c1<ncol; ++c1) {
    // Forward order for transpose, otherwise backwards
    c = tr ? c1 : ncol-1-c1;
    // Calculate scalar factor alpha = beta(c)*dot(v(:,c), x)
    alpha=0;
    for (k=colind[c]; k<colind[c+1]; ++k) alpha += v[k]*x[row[k]];
    alpha *= beta[c];
    // x -= alpha*v(:,c)
    for (k=colind[c]; k<colind[c+1]; ++k) x[row[k]] -= alpha*v[k];
  }
}

// SYMBOL "qr_trs"
// Solve for an (optionally transposed) upper triangular matrix R
template<typename T1>
void casadi_qr_trs(const casadi_int* sp_r, const T1* nz_r, T1* x, casadi_int tr) {
  // Extract sparsity
  casadi_int ncol=sp_r[1];
  const casadi_int *colind=sp_r+2, *row=sp_r+2+ncol+1;
  // Local variables
  casadi_int r, c, k;
  if (tr) {
    // Forward substitution
    for (c=0; c<ncol; ++c) {
      for (k=colind[c]; k<colind[c+1]; ++k) {
        r = row[k];
        if (r==c) {
          x[c] /= nz_r[k];
        } else {
          x[c] -= nz_r[k]*x[r];
        }
      }
    }
  } else {
    // Backward substitution
    for (c=ncol-1; c>=0; --c) {
      for (k=colind[c+1]-1; k>=colind[c]; --k) {
        r=row[k];
        if (r==c) {
          x[r] /= nz_r[k];
        } else {
          x[r] -= nz_r[k]*x[c];
        }
      }
    }
  }
}

// SYMBOL "qr_solve"
// Solve a factorized linear system
// len[w] >= max(ncol, nrow_ext)
template<typename T1>
void casadi_qr_solve(T1* x, casadi_int nrhs, casadi_int tr,
                     const casadi_int* sp_v, const T1* v, const casadi_int* sp_r, const T1* r,
                     const T1* beta, const casadi_int* pinv, T1* w) {
  casadi_int k, c;
  casadi_int nrow_ext = sp_v[0], ncol = sp_v[1];
  for (k=0; k<nrhs; ++k) {
    if (tr) {
      // ('P'Q R)' x = R'Q'P x = b <-> x = P' Q R' \ b
      // Copy to w
      casadi_copy(x, ncol, w);
      //  Solve for R'
      casadi_qr_trs(sp_r, r, w, 1);
      // Multiply by Q
      casadi_qr_mv(sp_v, v, beta, w, 0);
      // Multiply by P'
      for (c=0; c<ncol; ++c) x[c] = w[pinv[c]];
    } else {
      //P'Q R x = b <-> x = R \ Q' P b
      // Multiply with P
      // C-REPLACE "T1(0)" "0"
      casadi_fill(w, nrow_ext, T1(0));
      for (c=0; c<ncol; ++c) w[pinv[c]] = x[c];
      // Multiply with Q'
      casadi_qr_mv(sp_v, v, beta, w, 1);
      //  Solve for R
      casadi_qr_trs(sp_r, r, w, 0);
      // Copy to x
      casadi_copy(w, ncol, x);
    }
    x += ncol;
  }
}
#pragma GCC diagnostic pop

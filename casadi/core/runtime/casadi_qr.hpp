// NOLINT(legal/copyright)
// SYMBOL "house"
// Householder reflection
// Ref: Chapter 5, Direct Methods for Sparse Linear Systems by Tim Davis
template<typename T1>
T1 casadi_house(T1* x, T1* beta, casadi_int n) {
  // Local variable
  casadi_int i;
  T1 x0, sigma, s, sigma_is_zero, x0_nonpos;
  // Calculate norm
  x0 = x[0]; // Save x0 (overwritten below)
  sigma=0;
  for (i=1; i<n; ++i) sigma += x[i]*x[i];
  s = sqrt(x0*x0 + sigma); // s = norm(x)
  // Calculate consistently with symbolic datatypes (SXElem)
  sigma_is_zero = sigma==0;
  x0_nonpos = x0<=0;
  // C-REPLACE "if_else" "casadi_if_else"
  x[0] = if_else(sigma_is_zero, 1,
                 if_else(x0_nonpos, x0-s, -sigma/(x0+s)));
  *beta = if_else(sigma_is_zero, 2*x0_nonpos, -1/(s*x[0]));
  return s;
}

// SYMBOL "qr"
// Numeric QR factorization
// Ref: Chapter 5, Direct Methods for Sparse Linear Systems by Tim Davis
// len[x] = nrow
// sp_v = [nrow, ncol, 0, 0, ...] len[3 + ncol + nnz_v]
// len[v] nnz_v
// sp_r = [nrow, ncol, 0, 0, ...] len[3 + ncol + nnz_r]
// len[r] nnz_r
// len[beta] ncol
template<typename T1>
void casadi_qr(const casadi_int* sp_a, const T1* nz_a, T1* x,
               const casadi_int* sp_v, T1* nz_v, const casadi_int* sp_r, T1* nz_r, T1* beta,
               const casadi_int* prinv, const casadi_int* pc) {
   // Local variables
   casadi_int ncol, nrow, r, c, k, k1;
   T1 alpha;
   const casadi_int *a_colind, *a_row, *v_colind, *v_row, *r_colind, *r_row;
   // Extract sparsities
   ncol = sp_a[1];
   a_colind=sp_a+2; a_row=sp_a+2+ncol+1;
   nrow = sp_v[0];
   v_colind=sp_v+2; v_row=sp_v+2+ncol+1;
   r_colind=sp_r+2; r_row=sp_r+2+ncol+1;
   // Clear work vector
   for (r=0; r<nrow; ++r) x[r] = 0;
   // Loop over columns of R, A and V
   for (c=0; c<ncol; ++c) {
     // Copy (permuted) column of A to x
     for (k=a_colind[pc[c]]; k<a_colind[pc[c]+1]; ++k) x[prinv[a_row[k]]] = nz_a[k];
     // Use the equality R = (I-betan*vn*vn')*...*(I-beta1*v1*v1')*A to get
     // strictly upper triangular entries of R
     for (k=r_colind[c]; k<r_colind[c+1] && (r=r_row[k])<c; ++k) {
       // Calculate scalar factor alpha = beta(r)*dot(v(:,r), x)
       alpha = 0;
       for (k1=v_colind[r]; k1<v_colind[r+1]; ++k1) alpha += nz_v[k1]*x[v_row[k1]];
       alpha *= beta[r];
       // x -= alpha*v(:,r)
       for (k1=v_colind[r]; k1<v_colind[r+1]; ++k1) x[v_row[k1]] -= alpha*nz_v[k1];
       // Get r entry
       *nz_r++ = x[r];
       // Strictly upper triangular entries in x no longer needed
       x[r] = 0;
     }
     // Get V column
     for (k=v_colind[c]; k<v_colind[c+1]; ++k) {
       nz_v[k] = x[v_row[k]];
       // Lower triangular entries of x no longer needed
       x[v_row[k]] = 0;
     }
     // Get diagonal entry of R, normalize V column
     *nz_r++ = casadi_house(nz_v + v_colind[c], beta + c, v_colind[c+1] - v_colind[c]);
   }
 }

// SYMBOL "qr_mv"
// Multiply QR Q matrix from the right with a vector, with Q represented
// by the Householder vectors V and beta
// x = Q*x or x = Q'*x
// with Q = (I-beta(1)*v(:,1)*v(:,1)')*...*(I-beta(n)*v(:,n)*v(:,n)')
// len[x] >= nrow_ext
template<typename T1>
void casadi_qr_mv(const casadi_int* sp_v, const T1* v, const T1* beta, T1* x,
                  casadi_int tr) {
  // Local variables
  casadi_int ncol, c, c1, k;
  T1 alpha;
  const casadi_int *colind, *row;
  // Extract sparsity
  ncol=sp_v[1];
  colind=sp_v+2; row=sp_v+2+ncol+1;
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
  // Local variables
  casadi_int ncol, r, c, k;
  const casadi_int *colind, *row;
  // Extract sparsity
  ncol=sp_r[1];
  colind=sp_r+2; row=sp_r+2+ncol+1;
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
                     const T1* beta, const casadi_int* prinv, const casadi_int* pc, T1* w) {
  casadi_int k, c, nrow_ext, ncol;
  nrow_ext = sp_v[0]; ncol = sp_v[1];
  for (k=0; k<nrhs; ++k) {
    if (tr) {
      // (PR' Q R PC)' x = PC' R' Q' PR x = b <-> x = PR' Q R' \ PC b
      // Multiply by PC
      for (c=0; c<ncol; ++c) w[c] = x[pc[c]];
      //  Solve for R'
      casadi_qr_trs(sp_r, r, w, 1);
      // Multiply by Q
      casadi_qr_mv(sp_v, v, beta, w, 0);
      // Multiply by PR'
      for (c=0; c<ncol; ++c) x[c] = w[prinv[c]];
    } else {
      //PR' Q R PC x = b <-> x = PC' R \ Q' PR b
      // Multiply with PR
      for (c=0; c<nrow_ext; ++c) w[c] = 0;
      for (c=0; c<ncol; ++c) w[prinv[c]] = x[c];
      // Multiply with Q'
      casadi_qr_mv(sp_v, v, beta, w, 1);
      //  Solve for R
      casadi_qr_trs(sp_r, r, w, 0);
      // Multiply with PC'
      for (c=0; c<ncol; ++c) x[pc[c]] = w[c];
    }
    x += ncol;
  }
}

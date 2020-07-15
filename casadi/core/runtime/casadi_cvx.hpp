// NOLINT(legal/copyright)
// C-REPLACE "fmax" "casadi_fmax"
// C-REPLACE "nullptr" "0"

// SYMBOL "cvx_house"
/// Computes Householder vector
/// beta: scalar
/// v: vector of length nv
/// Returns 2-norm of v
///
/// Ref: Golub & Van Loan Alg 5.1.1
template<typename T1>
T1 casadi_cvx_house(T1* v, T1* beta, casadi_int nv) {
  // Local variable
  casadi_int i;
  T1 v0, sigma, s, v02;
  // Calculate norm
  v0 = v[0]; // Save v0 (overwritten below)
  sigma=0;
  for (i=1; i<nv; ++i) sigma += v[i]*v[i];
  s = sqrt(v0*v0 + sigma); // s = norm(v)
  if (sigma==0) {
    *beta = 0;
  } else {
    if (v0<=0) {
      v0 -= s;
    } else {
      v0 = -sigma/(v0+s);
    }
    v02 = v0*v0;
    *beta = 2*v02/(sigma+v02);
    v[0] = 1;
    for (i=1;i<nv;++i) v[i] /= v0;
  }
  return s;
}


// SYMBOL "cvx_house_apply_symm"
// Apply householder transform to a symmetric submatrix
// on dense A m-by-n matrix
//
// A is modified in-place
//
// s : stride
//     normally equal to m
//     when A is a submatrix of a bigger matrix, set equal to latter's number of rows
// v : compact housholder factorisation (length m)
//     First element (always one) is used to store beta
// p : length n
//
// Reference: Golub & Van Loan, Alg. 8.3.1
template<typename T1>
void casadi_cvx_house_apply_symm(casadi_int n, casadi_int k, T1* A, T1* p, T1* v) {
  casadi_int i, j, stride, N;
  T1 *a;
  T1 beta = v[0];
  v[0] = 1;
  stride = k+1;
  A+= k+1+n*k;
  N = n-k-1;

  // p <- beta * A(k+1:n,k+1:n) v
  casadi_clear(p, N);
  a = A+n;

  // Loop over columns
  for (i=0;i<N;++i) {
    p[i] += beta*(*a++)*v[i];
    // Loop over rows
    for (j=i+1;j<N;++j) {
      p[j] += beta*(*a)*v[i];
      p[i] += beta*(*a++)*v[j];
    }
    a += stride+i+1;
  }

  // p <- p - (beta p'v/2) v
  casadi_axpy(N, -beta*casadi_dot(N, p, v)/2, v, p);

  // Rank-2 update
  a = A+n;
  // Loop over columns
  for (i=0;i<N;++i) {
    *a++ -= 2*v[i]*p[i];
    // Loop over rows
    for (j=i+1;j<N;++j) {
      *a++ -= v[i]*p[j]+v[j]*p[i];
    }
    a += stride+i+1;
  }
  v[0] = beta;
}


// SYMBOL "cvx_tri"
// Tri-diagonalize a symmetric matrix in-place
// Results are in lower-triangular part
//
// Upper triangular part contains compact housholder factorisations
//
// A: n-by-n dense
// p: work vector; length n
template<typename T1>
void casadi_cvx_tri(T1* A, casadi_int n, T1* p) {
  T1 pp[1000];
  casadi_int k, N;
  T1 *A_base, *v;
  T1 beta;
  for (k=0;k<n-2;++k) {
    A_base = A+k+1+n*k;
    N = n-k-1;

    v = A+N*n;

    // Compute Householder transformation
    casadi_copy(A_base, N, v);

    // Assign 2-norm
    *A_base = casadi_cvx_house(v, &beta, N);

    v[0] = beta;
    casadi_cvx_house_apply_symm(n, k, A, pp, v);

  }
}


// SYMBOL "cvx_givens"
// Ref: Golub & Van Loan Alg. 5.1.3
template<typename T1>
void casadi_cvx_givens(T1 a, T1 b, T1* c, T1* s) {
  T1 r;
  if (b==0) {
    *c = 1;
    *s = 0;
  } else {
    if (fabs(b)>fabs(a)) {
      r = -a/b;
      *s = 1/sqrt(1+r*r);
      *c = (*s)*r;
    } else {
      r = -b/a;
      *c = 1/sqrt(1+r*r);
      *s = (*c)*r;
    }
  }
}


// SYMBOL "cvx_implicit_qr"
// Implicit Symmetric QR step with Wilkinson shift
//
// Tri-diagonal n-by-n matrix
// Diagonal: t_diag (length n)
// Off-diagonal: t_off (length n-1)
// cs: [c0 s0 c1 s1 ...] length 2*(n-1)
// Golub & Van Loan Alg. 8.3.2
template<typename T1>
void casadi_cvx_implicit_qr(casadi_int n, T1* t_diag, T1* t_off, T1* cs) {
  T1 d, mu, to2, x, z, c, s, t1, t2, d0, d1, o0, o1, sd;
  casadi_int i;
  d = 0.5*(t_diag[n-2]-t_diag[n-1]);
  to2 = t_off[n-2]*t_off[n-2];
  sd = 1;
  if (d<0) sd = -1;
  mu = t_diag[n-1]-to2/(d+sd*sqrt(d*d+to2));
  x = t_diag[0]-mu;
  z = t_off[0];
  for (i=0;i<n-1;++i) {
    // Compute Givens transformation
    casadi_cvx_givens(x, z, &c, &s);
    // T = G'TG (worked out with scalars)
    d0 = t_diag[i];
    d1 = t_diag[i+1];
    o0 = t_off[i];
    o1 = t_off[i+1];
    t1 = d0*c-o0*s;
    t2 = o0*c-d1*s;
    t_diag[i]   = c*t1-s*t2;
    t_off[i]    = s*t1+c*t2;
    t_diag[i+1] = d0*s*s+2*s*o0*c+d1*c*c;
    t_off[i+1] *= c;
    if (i>0) {
      t_off[i-1] = t_off[i-1]*c-z*s;
    }
    x = t_off[i];
    z = -s*o1;
    if (cs) {
      *cs++ = c;
      *cs++ = s;
    }
  }
}


// SYMBOL "cvx_symm_schur"
// Eigen-decomposition Q'TQ = D
// T tri-diagonal, with:
//   - t_diag the diagonal vector (length n)
//   - t_off the off-diagonal vector (length n-1)
//
// Eigenvalues can be read from returned t_diag
//
// tolerance greater than machine precision
//
// trace_meta: length 1+3*n_iter
// trace: length 2*(n-1)*n_iter
//
/// Golub & Van Loan Alg. 8.3.3
template<typename T1>
int casadi_cvx_symm_schur(casadi_int n, T1* t_diag, T1* t_off, T1 tol, casadi_int max_iter,
    casadi_int* trace_meta, T1* trace) {
  casadi_int i, p, q, sp, sq, trace_offset, nn;
  casadi_int* n_iter;
  n_iter = trace_meta++;

  trace_offset = 0;
  q = 0;
  *n_iter = 0;

  while (q<n) {
    if (*n_iter==max_iter) return 1;
    // Clip converged entries
    for (i=0;i<n-1;++i) {
      if (fabs(t_off[i])<=tol*(fabs(t_diag[i])+fabs(t_diag[i+1]))) {
        t_off[i] = 0;
      }
    }

    // Determine p, q
    p = 0;
    q = 0;
    sp = 0;
    sq = 0;
    for (i=0;i<n-1;++i) {
      if (t_off[n-i-2]==0 && sq==0) {
        q++;
      } else {
        sq = 1;
      }
      if (t_off[i]==0 && sp==0) {
        p++;
      } else {
        sp = 1;
      }
      if (q==n-1) {
        q = n;
        p = 0;
      }
    }

    nn = n-q-p;
    if (q<n) {
      casadi_cvx_implicit_qr(nn, t_diag+p, t_off+p, trace ? trace+trace_offset : nullptr);
      trace_offset += 2*(nn-1);

      if (trace_meta) {
        *trace_meta++ = nn;
        *trace_meta++ = p;
        *trace_meta++ = trace_offset;
      }
      (*n_iter)++;
    }
  }
  return 0;
}



// SYMBOL "cvx_givens_apply"
template<typename T1>
void casadi_cvx_givens_apply(casadi_int n, T1* q, T1 c, T1 s, casadi_int p) {
  T1 t1, t2, t3, t4, a, b;
  casadi_int i;
  // Update rows
  T1 *m = q;
  m += p;
  for (i=0;i<p;++i) {
    a = m[0];
    b = m[1];
    m[0] = c*a+s*b;
    m[1] = c*b-s*a;
    m+=n;
  }
  // Update central patch
  t1 = c*m[0]+s*m[1];
  t2 = c*m[1]+s*m[n+1];
  t3 = c*m[1]-s*m[0];
  t4 = c*m[n+1]-s*m[1];
  m[0] = c*t1+s*t2;
  m[1] = c*t2-s*t1;
  m[n+1] = c*t4-s*t3;
  // Update columns
  m = q+n*p+p+2;
  for (i=0;i<n-p-2;++i) {
    a = m[0];
    b = m[n];
    m[0] = c*a+s*b;
    m[n] = c*b-s*a;
    m++;
  }
}


// SYMBOL "cvx_house_apply"
/// Apply householder transform
/// on dense A m-by-n matrix
///
/// A is modified in-place
///
/// s : stride
///     normally equal to m
///     when A is a submatrix of a bigger matrix, set equal to latter's number of rows
/// v : compact housholder factorisation (length m)
///     First element (always one) is used to store beta
/// p : length n
///
template<typename T1>
void casadi_cvx_house_apply(casadi_int n, casadi_int m, casadi_int s, T1* A,
    T1* p, const T1* v) {
  casadi_int i, j;
  T1 *a;
  T1 beta;
  beta = v[0];

  // pi <- beta Aji vj
  casadi_clear(p, n);
  a = A;

  // Loop over columns
  for (i=0;i<n;++i) {
    p[i] += beta*a[0];
    // Loop over rows
    for (j=1;j<m;++j) {
      p[i] += beta*a[j]*v[j];
    }
    a += s;
  }

  a = A;
  // Loop over columns
  for (i=0;i<n;++i) {
    a[0] -= p[i];
    // Loop over rows
    for (j=1;j<m;++j) {
      a[j] -= v[j]*p[i];
    }
    a += s;
  }
}

template<typename T1>
T1 casadi_cvx_scalar(T1 epsilon, casadi_int reflect, T1 eig) {
  return fmax(epsilon, reflect ? fabs(eig) : eig);
}


// SYMBOL "cvx"
// Convexify a dense symmetric Hessian
//
// w real work vector: length max(n,2*(n-1)*n_iter)
// iw integer work vector: 1+3*n_iter
//
// tol:     tolerance for symmetric schur
// epsilon: minimum magnitude of eigenvalues
// reflect: when nonzero, reflect negative eigenvalues
template<typename T1>
int casadi_cvx(casadi_int n, T1 *A, T1 epsilon, T1 tol, casadi_int reflect, casadi_int max_iter,
    T1* w, casadi_int* iw) {
  casadi_int i, j, k, n_iter, nn, p, trace_offset;
  casadi_int *t_meta;
  T1 c, s, t_off0;
  T1 *cs, *t_diag, *t_off;

  // Short-circuit for empty matrices
  if (n==0) return 0;

  // Short-circuit for scalar matrices
  if (n==1) {
    A[0] = casadi_cvx_scalar(epsilon, reflect, A[0]);
    return 0;
  }

  casadi_cvx_tri(A, n, w);

  for (i=0;i<n;++i) {
    for (j=0;j<n;++j) {
      if (i-j>=2) {
        A[i+j*n] = 0;
      }
    }
  }

  // Represent tri-diagonal as vector pair (t_diag, t_off)
  t_off0 = A[1];
  t_diag = A;
  t_off = A+n;
  for (i=1;i<n;++i) {
    t_diag[i] = A[i+n*i];
  }
  t_off[0] = t_off0;
  for (i=1;i<n-1;++i) {
    t_off[i] = A[i+1+n*i];
  }

  // Diagonalize matrix by Symmetric QR
  if (casadi_cvx_symm_schur(n, t_diag, t_off, tol, max_iter, iw, w)) return 1;

  // Retain diagonals (eigenvalues)
  for (i=0;i<n;++i) {
    A[i+n*i] = casadi_cvx_scalar(epsilon, reflect, t_diag[i]);
  }

  // Reset other elements
  for (i=0;i<n;++i) {
    for (j=i+1;j<n;++j) A[j+i*n] = 0;
  }

  // Undo Symmetric QR
  n_iter = iw[0];
  t_meta = iw+3*(n_iter-1)+1;

  for (i=0;i<n_iter;++i) {
    nn = *t_meta++;
    p = *t_meta++;
    trace_offset = *t_meta++;
    cs = w+trace_offset;
    t_meta-= 6;
    for (j=0;j<nn-1;j++) {
      s = *--cs;
      c = *--cs;
      casadi_cvx_givens_apply(n, A, c, s, p+nn-j-2);
    }
  }

  // Undo triangularization
  for (k = n-3; k>=0; --k) {
    casadi_int N = n-k-1;
    T1 *v = A+N*n;
    casadi_cvx_house_apply_symm(n, k, A, w, v);
    casadi_cvx_house_apply(k+1, N, n, A+k+1, w, v);
  }

  return 0;
}

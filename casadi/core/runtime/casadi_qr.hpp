// NOLINT(legal/copyright)
// SYMBOL "etree"
// Calculate the elimination tree for a matrix
// len[w] >= ata ? ncol + nrow : ncol
// len[parent] == ncol
// Ref: Chapter 4, Direct Methods for Sparse Linear Systems by Tim Davis
inline
void casadi_etree(const int* sp, int* parent, int *w, int ata) {
  int r, c, k, rnext;
  // Extract sparsity
  int nrow = *sp++, ncol = *sp++;
  const int *colind = sp, *row = sp+ncol+1;
  // Highest known ascestor of a node
  int *ancestor=w;
  // Path for A'A
  int *prev;
  if (ata) {
    prev=w+ncol;
    for (r=0; r<nrow; ++r) prev[r] = -1;
  }
  // Loop over columns
  for (c=0; c<ncol; ++c) {
    parent[c] = -1; // No parent yet
    ancestor[c] = -1; // No ancestor
    // Loop over nonzeros
    for (k=colind[c]; k<colind[c+1]; ++k) {
      r = row[k];
      if (ata) r = prev[r];
      // Traverse from r to c
      while (r!=-1 && r<c) {
        rnext = ancestor[r];
        ancestor[r] = c;
        if (rnext==-1) parent[r] = c;
        r = rnext;
      }
      if (ata) prev[row[k]] = c;
    }
  }
}

// SYMBOL "postorder_dfs"
// Traverse an elimination tree using depth first search
// Ref: Chapter 4, Direct Methods for Sparse Linear Systems by Tim Davis
inline
int casadi_postorder_dfs(int j, int k, int* head, int* next,
                         int* post, int* stack) {
  int i, p, top=0;
  stack[0] = j;
  while (top>=0) {
    p = stack[top];
    i = head[p];
    if (i==-1) {
      // No children
      top--;
      post[k++] = p;
    } else {
      // Add to stack
      head[p] = next[i];
      stack[++top] = i;
    }
  }
  return k;
}

// SYMBOL "postorder"
// Calculate the postorder permuation
// Ref: Chapter 4, Direct Methods for Sparse Linear Systems by Tim Davis
// len[w] >= 3*n
// len[post] == n
inline
void casadi_postorder(const int* parent, int n, int* post, int* w) {
  int j, k=0;
  // Work vectors
  int *head, *next, *stack;
  head=w; w+=n;
  next=w; w+=n;
  stack=w; w+=n;
  // Empty linked lists
  for (j=0; j<n; ++j) head[j] = -1;
  // Traverse nodes in reverse order
  for (j=n-1; j>=0; --j) {
    if (parent[j]!=-1) {
      next[j] = head[parent[j]];
      head[parent[j]] = j;
    }
  }
  for (j=0; j<n; j++) {
    if (parent[j]==-1) {
      k = casadi_postorder_dfs(j, k, head, next, post, stack);
    }
  }
}

// SYMBOL "leaf"
// Needed by casadi_qr_colind
// Ref: Chapter 4, Direct Methods for Sparse Linear Systems by Tim Davis
inline
int casadi_leaf(int i, int j, const int* first, int* maxfirst,
                 int* prevleaf, int* ancestor, int* jleaf) {
  int q, s, sparent, jprev;
  *jleaf = 0;
  // Quick return if j is not a leaf
  if (i<=j || first[j]<=maxfirst[i]) return -1;
  // Update max first[j] seen so far
  maxfirst[i] = first[j];
  // Previous leaf of ith subtree
  jprev = prevleaf[i];
  prevleaf[i] = j;
  // j is first or subsequent leaf
  *jleaf = (jprev == -1) ? 1 : 2;
  // if first leaf, q is root of ith subtree
  if (*jleaf==1) return i;
  // Path compression
  for (q=jprev; q!=ancestor[q]; q=ancestor[q]) {}
  for (s=jprev; s!=q; s=sparent) {
    sparent = ancestor[s];
    ancestor[s] = q;
  }
  // Return least common ancestor
  return q;
}

// SYMBOL "imin"
// Smallest of two integers
inline
int casadi_imin(int a, int b) {
  return a<b ? a : b;
}

// SYMBOL "imax"
// Largest of two integers
inline
int casadi_imax(int a, int b) {
  return a>b ? a : b;
}

// SYMBOL "qr_colind"
// Calculate the column offsets for the QR R matrix
// Ref: Chapter 4, Direct Methods for Sparse Linear Systems by Tim Davis
// len[counts] = ncol
// len[w] >= 5*ncol + nrow + 1
inline
int casadi_qr_counts(const int* tr_sp, const int* parent,
                     const int* post, int* counts, int* w) {
  int ncol = *tr_sp++, nrow = *tr_sp++;
  const int *rowind=tr_sp, *col=tr_sp+nrow+1;
  int i, j, k, J, p, q, jleaf, *maxfirst, *prevleaf,
    *ancestor, *head=0, *next=0, *first;
  // Work vectors
  ancestor=w; w+=ncol;
  maxfirst=w; w+=ncol;
  prevleaf=w; w+=ncol;
  first=w; w+=ncol;
  head=w; w+=ncol+1;
  next=w; w+=nrow;
  // Find first [j]
  for (k=0; k<ncol; ++k) first[k]=-1;
  for (k=0; k<ncol; ++k) {
    j=post[k];
    // counts[j]=1 if j is a leaf
    counts[j] = (first[j]==-1) ? 1 : 0;
    for (; j!=-1 && first[j]==-1; j=parent[j]) first[j]=k;
  }
  // Invert post (use ancestor as work vector)
  for (k=0; k<ncol; ++k) ancestor[post[k]] = k;
  for (k=0; k<ncol+1; ++k) head[k]=-1;
  for (i=0; i<nrow; ++i) {
    for (k=ncol, p=rowind[i]; p<rowind[i+1]; ++p) {
      k = casadi_imin(k, ancestor[col[p]]);
    }
    // Place row i in linked list k
    next[i] = head[k];
    head[k] = i;
  }

  // Clear workspace
  for (k=0; k<ncol; ++k) maxfirst[k]=-1;
  for (k=0; k<ncol; ++k) prevleaf[k]=-1;
  // Each node in its own set
  for (i=0; i<ncol; ++i) ancestor[i]=i;
  for (k=0; k<ncol; ++k) {
    // j is the kth node in the postordered etree
    j=post[k];
    if (parent[j]!=-1) counts[parent[j]]--; // j is not a root
    J=head[k];
    while (J!=-1) { // J=j for LL' = A case
      for (p=rowind[J]; p<rowind[J+1]; ++p) {
        i=col[p];
        q = casadi_leaf(i, j, first, maxfirst, prevleaf, ancestor, &jleaf);
        if (jleaf>=1) counts[j]++; // A(i,j) is in skeleton
        if (jleaf==2) counts[q]--; // account for overlap in q
      }
      J = next[J];
    }
    if (parent[j]!=-1) ancestor[j]=parent[j];
  }
  // Sum up counts of each child
  for (j=0; j<ncol; ++j) {
    if (parent[j]!=-1) counts[parent[j]] += counts[j];
  }

  // Sum of counts
  int sum_counts = 0;
  for (j=0; j<ncol; ++j) sum_counts += counts[j];
  return sum_counts;
}

// SYMBOL "qr_nnz"
// Calculate the number of nonzeros in the QR V matrix
// Ref: Chapter 5, Direct Methods for Sparse Linear Systems by Tim Davis
// len[w] >= nrow + 3*ncol
// len[pinv] == nrow + ncol
// len[leftmost] == nrow
inline
int casadi_qr_nnz(const int* sp, int* pinv, int* leftmost,
                  const int* parent, int* nrow_ext, int* w) {
  // Extract sparsity
  int nrow = sp[0], ncol = sp[1];
  const int *colind=sp+2, *row=sp+2+ncol+1;
  // Work vectors
  int *next=w; w+=nrow;
  int *head=w; w+=ncol;
  int *tail=w; w+=ncol;
  int *nque=w; w+=ncol;
  // Local variables
  int r, c, k, pa;
  // Clear queue
  for (c=0; c<ncol; ++c) head[c] = -1;
  for (c=0; c<ncol; ++c) tail[c] = -1;
  for (c=0; c<ncol; ++c) nque[c] = 0;
  for (r=0; r<nrow; ++r) leftmost[r] = -1;
  // leftmost[r] = min(find(A(r,:)))
  for (c=ncol-1; c>=0; --c) {
    for (k=colind[c]; k<colind[c+1]; ++k) {
      leftmost[row[k]] = c;
    }
  }
  // Scan rows in reverse order
  for (r=nrow-1; r>=0; --r) {
    pinv[r] = -1; // row r not yet ordered
    c=leftmost[r];

    if (c==-1) continue; // row r is empty
    if (nque[c]++ == 0) tail[c]=r; // first row in queue c
    next[r] = head[c]; // put r at head of queue c
    head[c] = r;
  }
  // Find row permutation and nnz(V)
  int v_nnz = 0;
  int nrow_new = nrow;
  for (c=0; c<ncol; ++c) {
    r = head[c]; // remove r from queue c
    v_nnz++; // count V(c,c) as nonzero
    if (r<0) r=nrow_new++; // add a fictitious row
    pinv[r] = c; // associate row r with V(:,c)
    if (--nque[c]<=0) continue; // skip if V(c+1,nrow,c) is empty
    v_nnz += nque[c]; // nque[c] is nnz(V(c+1:nrow, c))
    if ((pa=parent[c]) != -1) {
      // Move all rows to parent of c
      if (nque[pa]==0) tail[pa] = tail[c];
      next[tail[c]] = head[pa];
      head[pa] = next[r];
      nque[pa] += nque[c];
    }
  }
  for (r=0; r<nrow; ++r) if (pinv[r]<0) pinv[r] = c++;
  if (nrow_ext) *nrow_ext = nrow_new;
  return v_nnz;
}

// SYMBOL "qr_init"
// Setup QP solver
// Ref: Chapter 5, Direct Methods for Sparse Linear Systems by Tim Davis
// len[w] >= nrow + 7*ncol + 1
// len[pinv] == nrow + ncol
// len[leftmost] == nrow
inline
void casadi_qr_init(const int* sp, const int* sp_tr,
                    int* leftmost, int* parent, int* pinv,
                    int* nrow_ext, int* v_nnz, int* r_nnz, int* w) {
  // Extract sparsity
  int ncol = sp[1];
  // Calculate elimination tree for A'A
  casadi_etree(sp, parent, w, 1); // len[w] >= nrow+ncol
  // Calculate postorder
  int* post = w; w += ncol;
  casadi_postorder(parent, ncol, post, w); // len[w] >= 3*ncol
  // Calculate nnz in R
  *r_nnz = casadi_qr_counts(sp_tr, parent, post, w, w+ncol);
  // Calculate nnz in V
  *v_nnz = casadi_qr_nnz(sp, pinv, leftmost, parent, nrow_ext, w);
}

// SYMBOL "qr_row"
// Get the row indices for V and R in QR factorization
// Ref: Chapter 5, Direct Methods for Sparse Linear Systems by Tim Davis
// Note: nrow <= nrow_ext <= nrow+ncol
// len[iw] = nrow_ext + ncol
// len[x] = nrow_ext
// sp_v = [nrow_ext, ncol, 0, 0, ...] len[3 + ncol + nnz_v]
// len[v] nnz_v
// sp_r = [nrow_ext, ncol, 0, 0, ...] len[3 + ncol + nnz_r]
// len[r] nnz_r
// len[beta] ncol
inline
void casadi_qr_sparsities(const int* sp_a, int nrow_ext, int* sp_v, int* sp_r,
                          const int* leftmost, const int* parent, const int* pinv,
                          int* iw) {
  // Extract sparsities
  int ncol = sp_a[1];
  const int *colind=sp_a+2, *row=sp_a+2+ncol+1;
  int *v_colind=sp_v+2, *v_row=sp_v+2+ncol+1;
  int *r_colind=sp_r+2, *r_row=sp_r+2+ncol+1;
  // Specify dimensions of V and R
  sp_v[0] = sp_r[0] = nrow_ext;
  sp_v[1] = sp_r[1] = ncol;
  // Work vectors
  int* s = iw; iw += ncol;
  // Local variables
  int r, c, k, k1, top, len, k2, r2;
  // Clear w to mark nodes
  for (r=0; r<nrow_ext; ++r) iw[r] = -1;
  // Number of nonzeros in v and r
  int nnz_r=0, nnz_v=0;
  // Compute V and R
  for (c=0; c<ncol; ++c) {
    // R(:,c) starts here
    r_colind[c] = nnz_r;
    // V(:, c) starts here
    v_colind[c] = k1 = nnz_v;
    // Add V(c,c) to pattern of V
    iw[c] = c;
    v_row[nnz_v++] = c;
    top = ncol;
    for (k=colind[c]; k<colind[c+1]; ++k) {
      r = leftmost[row[k]]; // r = min(find(A(r,:))
      // Traverse up c
      for (len=0; iw[r]!=c; r=parent[r]) {
        s[len++] = r;
        iw[r] = c;
      }
      while (len>0) s[--top] = s[--len]; // push path on stack
      r = pinv[row[k]]; // r = permuted row of A(:,c)
      if (r>c && iw[r]<c) {
        v_row[nnz_v++] = r; // add r to pattern of V(:,c)
        iw[r] = c;
      }
    }
    // For each r in pattern of R(:,c)
    for (k = top; k<ncol; ++k) {
      // R(r,c) is nonzero
      r = s[k];
      // Apply (V(r), beta(r)) to x: x -= v*beta*v'*x
      r_row[nnz_r++] = r;
      if (parent[r]==c) {
        for (k2=v_colind[r]; k2<v_colind[r+1]; ++k2) {
          r2 = v_row[k2];
          if (iw[r2]<c) {
            iw[r2] = c;
            v_row[nnz_v++] = r2;
          }
        }
      }
    }
    // R(c,c) = norm(x)
    r_row[nnz_r++] = c;
  }
  // Finalize R, V
  r_colind[ncol] = nnz_r;
  v_colind[ncol] = nnz_v;
}

// SYMBOL "house"
// Householder reflection
// Ref: Chapter 5, Direct Methods for Sparse Linear Systems by Tim Davis
template<typename T1>
T1 casadi_house(T1* x, T1* beta, int n) {
  // Local variable
  int i;
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
template<typename T1>
void casadi_qr(const int* sp_a, const T1* nz_a, int* iw, T1* x,
               const int* sp_v, T1* nz_v, const int* sp_r, T1* nz_r, T1* beta,
               const int* leftmost, const int* parent, const int* pinv) {
  // Extract sparsities
  int ncol = sp_a[1];
  const int *colind=sp_a+2, *row=sp_a+2+ncol+1;
  int nrow_ext = sp_v[0];
  const int *v_colind=sp_v+2, *v_row=sp_v+2+ncol+1;
  // Work vectors
  int* s = iw; iw += ncol;
  // Local variables
  int r, c, k, k1, top, len, k2, r2;
  T1 tau;
  // Clear workspace x
  for (r=0; r<nrow_ext; ++r) x[r] = 0;
  // Clear w to mark nodes
  for (r=0; r<nrow_ext; ++r) iw[r] = -1;
  // Number of nonzeros in v and r
  int nnz_r=0, nnz_v=0;
  // Compute V and R
  for (c=0; c<ncol; ++c) {
    // V(:, c) starts here
    k1 = nnz_v;
    // Add V(c,c) to pattern of V
    iw[c] = c;
    nnz_v++;
    top = ncol;
    for (k=colind[c]; k<colind[c+1]; ++k) {
      r = leftmost[row[k]]; // r = min(find(A(r,:))
      // Traverse up c
      for (len=0; iw[r]!=c; r=parent[r]) {
        s[len++] = r;
        iw[r] = c;
      }
      while (len>0) s[--top] = s[--len]; // push path on stack
      r = pinv[row[k]]; // r = permuted row of A(:,c)
      x[r] = nz_a[k]; // x(r) = A(:,c)
      if (r>c && iw[r]<c) {
        nnz_v++; // add r to pattern of V(:,c)
        iw[r] = c;
      }
    }
    // For each r in pattern of R(:,c)
    for (k = top; k<ncol; ++k) {
      // R(r,c) is nonzero
      r = s[k];
      // Apply (V(r), beta(r)) to x: x -= v*beta*v'*x
      tau=0;
      for (k2=v_colind[r]; k2<v_colind[r+1]; ++k2) tau += nz_v[k2] * x[v_row[k2]];
      tau *= beta[r];
      for (k2=v_colind[r]; k2<v_colind[r+1]; ++k2) x[v_row[k2]] -= nz_v[k2]*tau;
      nz_r[nnz_r++] = x[r];
      x[r] = 0;
      if (parent[r]==c) {
        for (k2=v_colind[r]; k2<v_colind[r+1]; ++k2) {
          r2 = v_row[k2];
          if (iw[r2]<c) {
            iw[r2] = c;
            nnz_v++;
          }
        }
      }
    }
    // Gather V(:,c) = x
    for (k=k1; k<nnz_v; ++k) {
      nz_v[k] = x[v_row[k]];
      x[v_row[k]] = 0;
    }
    // R(c,c) = norm(x)
    nz_r[nnz_r++] = casadi_house(nz_v + k1, beta + c, nnz_v-k1);
  }
}

// SYMBOL "qr_mv"
// Multiply QR Q matrix from the right with a vector, with Q represented
// by the Householder vectors V and beta
// x = Q*x or x = Q'*x
// with Q = (I-beta(1)*v(:,1)*v(:,1)')*...*(I-beta(n)*v(:,n)*v(:,n)')
// len[x] >= nrow_ext
template<typename T1>
void casadi_qr_mv(const int* sp_v, const T1* v, const T1* beta, T1* x,
                  int tr) {
  // Extract sparsity
  int ncol=sp_v[1];
  const int *colind=sp_v+2, *row=sp_v+2+ncol+1;
  // Local variables
  int c, c1, k;
  T1 alpha;
  // Loop over vectors
  for (c1=0; c1<ncol; ++c1) {
    // Forward order for transpose, otherwards backwards
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
void casadi_qr_trs(const int* sp_r, const T1* nz_r, T1* x, int tr) {
  // Extract sparsity
  int ncol=sp_r[1];
  const int *colind=sp_r+2, *row=sp_r+2+ncol+1;
  // Local variables
  int r, c, k;
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
void casadi_qr_solve(T1* x, int nrhs, int tr,
                     const int* sp_v, const T1* v, const int* sp_r, const T1* r,
                     const T1* beta, const int* pinv, T1* w) {
  int k, c;
  int nrow_ext = sp_v[0], ncol = sp_v[1];
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

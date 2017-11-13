// NOLINT(legal/copyright)
// SYMBOL "ldl_colind"
// Calculate the column offsets for the L factor (strictly lower entries only)
// for an LDL^T factorization
// Ref: User Guide for LDL by Tim Davis
// len[colind] = ncol+1
// len[parent] = ncol
// len[w] >= ncol
inline
void casadi_ldl_colind(const int* sp, int* parent,
                        int* l_colind, int* w) {
  int n = sp[0];
  const int *colind=sp+2, *row=sp+2+n+1;
  // Local variables
  int r, c, k;
  // Work vectors
  int* visited=w; w+=n;
  // Loop over columns
  for (c=0; c<n; ++c) {
    // L(c,:) pattern: all nodes reachable in etree from nz in A(0:c-1,c)
    parent[c] = -1; // parent of c is not yet known
    visited[c] = c; // mark node c as visited
    l_colind[1+c] = 0; // count of nonzeros in column c of L
    // Loop over strictly upper triangular entries A
    for (k=colind[c]; k<colind[c+1] && (r=row[k])<c; ++k) {
      // Follow path from r to root of etree, stop at visited node
      while (visited[r]!=c) {
        // Find parent of r if not yet determined
        if (parent[r]==-1) parent[r]=c;
        l_colind[1+r]++; // L(c,r) is nonzero
        visited[r] = c; // mark r as visited
        r=parent[r]; // proceed to parent row
      }
    }
  }
  // Cumsum
  l_colind[0] = 0;
  for (c=0; c<n; ++c) l_colind[c+1] += l_colind[c];
}

// SYMBOL "ldl_row"
// Calculate the row indices for the L factor (strictly lower entries only)
// for an LDL^T factorization
// Ref: User Guide for LDL by Tim Davis
// len[w] >= n
inline
void casadi_ldl_row(const int* sp, const int* parent, int* l_colind, int* l_row, int *w) {
  // Extract sparsity
  int n = sp[0];
  const int *colind = sp+2, *row = sp+n+3;
  // Work vectors
  int *visited=w; w+=n;
  // Local variables
  int r, c, k;
  // Compute nonzero pattern of kth row of L
  for (c=0; c<n; ++c) {
    // Not yet visited
    visited[c] = c;
    // Loop over nonzeros in upper triangular half
    for (k=colind[c]; k<colind[c+1] && (r=row[k])<c; ++k) {
      // Loop over dependent rows
      while (visited[r]!=c) {
        l_row[l_colind[r]++] = c; // L(c,r) is nonzero
        visited[r] = c; // mark r as visited
        r=parent[r]; // proceed to parent row
      }
    }
  }
  // Restore l_colind by shifting it forward
  k=0;
  for (c=0; c<n; ++c) {
    r=l_colind[c];
    l_colind[c]=k;
    k=r;
  }
}

// SYMBOL "ldl"
// Calculate the nonzeros of the L factor (strictly lower entries only)
// as well as D for an LDL^T factorization
// Ref: User Guide for LDL by Tim Davis
// len[iw] >= 2*n
// len[w] >= n
template<typename T1>
void casadi_ldl(const int* sp_a, const int* parent, const int* sp_l,
                 const T1* a, T1* l, T1* d, int *iw, T1* w) {
  // Extract sparsities
  int n = sp_a[0];
  const int *colind = sp_a+2, *row = sp_a+n+3;
  const int *l_colind = sp_l+2, *l_row = sp_l+n+3;
  // Work vectors
  int *visited=iw; iw+=n;
  int *currcol=iw; iw+=n;
  T1* y = w; w+=n;
  // Local variables
  int r, c, k, k2;
  T1 yr;
  // Keep track of current nonzero for each column of L
  for (c=0; c<n; ++c) currcol[c] = l_colind[c];
  // Compute nonzero pattern of kth row of L
  for (c=0; c<n; ++c) {
    // Not yet visited
    visited[c] = c;
    // Get nonzeros of column c in a dense vector
    y[c]=0; // Make sure y is all-zero until index c
    for (k=colind[c]; k<colind[c+1] && (r=row[k])<=c; ++k) y[r] = a[k];
    // Get D(c,c) and clear Y(c)
    d[c] = y[c];
    y[c] = 0;
    // Loop over matching entries in L(:,c), i.e. L(c,:)
    for (k=colind[c]; k<colind[c+1] && (r=row[k])<c; ++k) {
      while (visited[r]!=c) {
        // Get and clear y(r)
        yr = y[r];
        y[r] = 0;
        for (k2=l_colind[r]; k2<l_colind[r+1]; ++k2) y[l_row[k2]] -= l[k2]*yr;
        // The nonzero entry L(c,r)
        d[c] -= (l[currcol[r]++]=yr/d[r]) * yr;
        visited[r] = c; // mark r as visited
        r=parent[r]; // proceed to parent row
      }
    }
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
    //      LDL'x = b <=> x = L\D\L'\b
    //  Solve for L'
    casadi_ldl_trs(sp_l, l, x, 0);
    // Divide by D
    for (i=0; i<n; ++i) x[i] /= d[i];
    // Solve for L
    casadi_ldl_trs(sp_l, l, x, 1);
    // Next rhs
    x += n;
  }
}

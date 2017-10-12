// NOLINT(legal/copyright)
// SYMBOL "etree"
// Calculate the elimination tree for a matrix
// len[w] >= ata ? ncol + nrow : ncol
// len[parent] == ncol
// Ref: Section 4.1, Direct Methods for Sparse Linear Systems by Tim Davis
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
// Ref: Section 4.3, Direct Methods for Sparse Linear Systems by Tim Davis
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
// Ref: Section 4.3, Direct Methods for Sparse Linear Systems by Tim Davis
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

// SYMBOL "chol_colind"
// Calculate the column offsets for the L factor (strictly lower entries only)
// for an LDL^T factorization
// Ref: User Guide for LDL by Tim Davis
// len[colind] = ncol+1
// len[parent] = ncol
// len[w] >= ncol
inline
void casadi_chol_colind(const int* sp, int* parent,
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

// SYMBOL "chol_row"
// Calculate the row indices for the L factor (strictly lower entries only)
// for an LDL^T factorization
// Ref: User Guide for LDL by Tim Davis
// len[w] >= n
inline
void casadi_chol_row(const int* sp, const int* parent, int* l_colind, int* l_row, int *w) {
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

// SYMBOL "leaf"
// Needed by casadi_lu_colind
// Ref: Section 4.4, Direct Methods for Sparse Linear Systems by Tim Davis
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

// SYMBOL "lu_colind"
// Calculate the column offsets for the LU L factor
// Ref: Section 4.5, Direct Methods for Sparse Linear Systems by Tim Davis
// len[colind] = ncol+1
// len[w] >= 5*ncol + nrow + 1
// C-REPLACE "std::min" "casadi_min"
inline
void casadi_lu_colind(const int* tr_sp, const int* parent,
                      const int* post, int* l_colind, int* w) {
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
    // l_colind[j]=1 if j is a leaf
    l_colind[1+j] = (first[j]==-1) ? 1 : 0;
    for (; j!=-1 && first[j]==-1; j=parent[j]) first[j]=k;
  }
  // Invert post (use ancestor as work vector)
  for (k=0; k<ncol; ++k) ancestor[post[k]] = k;
  for (k=0; k<ncol+1; ++k) head[k]=-1;
  for (i=0; i<nrow; ++i) {
    for (k=ncol, p=rowind[i]; p<rowind[i+1]; ++p) {
      k = std::min(k, ancestor[col[p]]);
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
    if (parent[j]!=-1) l_colind[1+parent[j]]--; // j is not a root
    J=head[k];
    while (J!=-1) { // J=j for LL' = A case
      for (p=rowind[J]; p<rowind[J+1]; ++p) {
        i=col[p];
        q = casadi_leaf(i, j, first, maxfirst, prevleaf, ancestor, &jleaf);
        if (jleaf>=1) l_colind[1+j]++; // A(i,j) is in skeleton
        if (jleaf==2) l_colind[1+q]--; // account for overlap in q
      }
      J = next[J];
    }
    if (parent[j]!=-1) ancestor[j]=parent[j];
  }
  // Sum up counts of each child
  for (j=0; j<ncol; ++j) {
    if (parent[j]!=-1) l_colind[1+parent[j]] += l_colind[1+j];
  }

  // Cumsum
  l_colind[0] = 0;
  for (j=0; j<ncol; ++j) {
    l_colind[j+1] += l_colind[j];
  }
}

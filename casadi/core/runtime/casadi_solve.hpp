// NOLINT(legal/copyright)
// SYMBOL "etree"
// Calculate the elimination tree for a matrix
// len[w] >= col ? ncol + nrow : ncol
// len[parent] == ncol
// Ref: Section 4.1, Direct Methods for Sparse Linear Systems by Tim Davis
inline
void casadi_etree(const int* sp, int* parent, int *w, int col) {
  int r, c, k, rnext;
  // Extract sparsity
  int nrow = *sp++, ncol = *sp++;
  const int *colind = sp, *row = sp+ncol+1;
  // Highest known ascestor of a node
  int *ancestor=w;
  // Path for A'A
  int *prev;
  if (col) {
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
      if (col) r = prev[r];
      // Traverse from r to c
      while (r!=-1 && r<c) {
        rnext = ancestor[r];
        ancestor[r] = c;
        if (rnext==-1) parent[r] = c;
        r = rnext;
      }
      if (col) prev[row[k]] = c;
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
  int *head=w, *next=w+n, *stack=w+2*n;
  for (j=0; j<n; ++j) head[j] = -1;
  for (j=n-1; j>=0; --j) {
    if (parent[j]!=-1) {
      next[j] = head[parent[j]];
      head[parent[j]] = j;
    }
  }
  for (j=0; j<n; j++) {
    if (parent[j]!=-1) {
      k = casadi_postorder_dfs(j, k, head, next, post, stack);
    }
  }
}

// SYMBOL "leaf"
// Needed by casadi_symbfact
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
  for (q=jprev; q!=ancestor[q]; q=ancestor[q]);
  for (s=jprev; s!=q; s=sparent) {
    sparent = ancestor[s];
    ancestor[s] = q;
  }
  // Return least common ancestor
  return q;
}

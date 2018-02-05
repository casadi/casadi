// NOLINT(legal/copyright)
// SYMBOL "lb_eig"
// Use Gershgorin to finds upper and lower bounds on the eigenvalues
template<typename T1>
double casadi_lb_eig(const casadi_int* sp_h, const T1* h) {
  // Get sparsity
  casadi_int ncol = sp_h[1];
  const casadi_int *colind = sp_h+2, *row = sp_h+ncol+3;
  // Return value
  T1 lb_eig;
  // Local variables
  int c, k;
  T1 center, radius;
  for (c=0; c<ncol; ++c) {
    // Calculate Gershgorin discs
    center = 0;
    radius = 0;
    for (k=colind[c]; k<colind[c+1]; ++k) {
      if (row[k]==c) {
        center = h[k];
      } else {
        radius += std::fabs(h[k]);
      }
    }
    // Update the eigenvalue estimates
    if (c==0) {
      lb_eig = center - radius;
    } else {
      lb_eig = std::fmin(lb_eig, center - radius);
    }
  }
  return lb_eig;
}

// SYMBOL "regularize"
// Add a multiple of the identity matrix to the diagonal
template<typename T1>
void casadi_regularize(const casadi_int* sp_h, T1* h, T1 reg) {
  // Get sparsity
  casadi_int ncol = sp_h[1];
  const casadi_int *colind = sp_h+2, *row = sp_h+ncol+3;
  // Local variables
  int c, k;
  // Shift diagonal entries
  for (c=0; c<ncol; ++c) {
    for (k=colind[c]; k<colind[c+1]; ++k) {
      if (row[k]==c) h[k] += reg;
    }
  }
}

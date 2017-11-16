// NOLINT(legal/copyright)
// SYMBOL "interpn_grad"
template<typename T1>
void casadi_interpn_grad(T1* grad, int ndim, const T1* grid, const int* offset, const T1* values, const T1* x, const int* lookup_mode, int m, int* iw, T1* w) { // NOLINT(whitespace/line_length)
  // Quick return
  if (!grad) return;
  // Work vectors
  T1* alpha = w; w += ndim;
  T1* coeff = w; w += ndim;
  T1* v = w; w+= m;
  int* index = iw; iw += ndim;
  int* corner = iw; iw += ndim;

  // Left index and fraction of interval
  casadi_interpn_weights(ndim, grid, offset, x, alpha, index, lookup_mode);
  // Loop over all corners, add contribution to output
  casadi_fill_int(corner, ndim, 0);
  casadi_fill(grad, ndim*m, 0.);
  do {
    // Get coefficients
    casadi_fill(v, m, 0.);
    casadi_interpn_interpolate(v, ndim, offset, values,
      alpha, index, corner, coeff, m);
    // Propagate to alpha
    int i, j;
    for (i=ndim-1; i>=0; --i) {
      if (corner[i]) {
        for (j=0; j<m; ++j) {
          grad[i*m+j] += v[j]*coeff[i];
          v[j] *= alpha[i];
        }
      } else {
        for (j=0; j<m; ++j) {
          grad[i*m+j] -= v[j]*coeff[i];
          v[j] *= 1-alpha[i];
        }
      }
    }
  } while (casadi_flip(corner, ndim));
  // Propagate to x
  int i, k;
  for (i=0; i<ndim; ++i) {
    const T1* g = grid + offset[i];
    int j = index[i];
    T1 delta =  g[j+1]-g[j];
    for (k=0;k<m;++k) grad[k] /= delta;
    grad += m;
  }
}

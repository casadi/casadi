// NOLINT(legal/copyright)
// SYMBOL "interpn_weights"
template<typename T1>
void casadi_interpn_weights(int ndim, const T1* grid, const int* offset, const T1* x, T1* alpha, int* index, const int* lookup_mode) { // NOLINT(whitespace/line_length)
  // Left index and fraction of interval
  int i;
  for (i=0; i<ndim; ++i) {
    // Grid point
    T1 xi = x ? x[i] : 0;
    // Grid
    const T1* g = grid + offset[i];
    int ng = offset[i+1]-offset[i];
    // Find left index
    int j = index[i] = casadi_low(xi, g, ng, lookup_mode[i]);
    // Get interpolation/extrapolation alpha
    alpha[i] = (xi-g[j])/(g[j+1]-g[j]);
  }
}

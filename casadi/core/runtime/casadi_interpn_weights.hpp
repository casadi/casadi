template<typename real_t>
void CASADI_PREFIX(interpn_weights)(int ndim, const real_t* grid, const int* offset, const real_t* x, real_t* alpha, int* index, const int* lookup_mode) {
  /* Left index and fraction of interval */
  int i;
  for (i=0; i<ndim; ++i) {
    /* Grid point */
    real_t xi = x ? x[i] : 0;
    /* Grid */
    const real_t* g = grid + offset[i];
    int ng = offset[i+1]-offset[i];
    /* Find left index */
    int j = index[i] = CASADI_PREFIX(low)(xi, g, ng, lookup_mode[i]);
    /* Get interpolation/extrapolation alpha */
    alpha[i] = (xi-g[j])/(g[j+1]-g[j]);
  }
}

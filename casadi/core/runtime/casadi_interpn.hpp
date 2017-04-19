template<typename real_t>
real_t CASADI_PREFIX(interpn)(int ndim, const real_t* grid, const int* offset, const real_t* values, const real_t* x, const int* lookup_mode, int* iw, real_t* w) {
  /* Work vectors */
  real_t* alpha = w; w += ndim;
  int* index = iw; iw += ndim;
  int* corner = iw; iw += ndim;
  /* Left index and fraction of interval */
  CASADI_PREFIX(interpn_weights)(ndim, grid, offset, x, alpha, index, lookup_mode);
  /* Loop over all corners, add contribution to output */
  CASADI_PREFIX(fill_int)(corner, ndim, 0);
  real_t ret = 0;
  do {
    real_t* coeff = 0;
    ret += CASADI_PREFIX(interpn_interpolate)(ndim, offset, values,
      alpha, index, corner, coeff);
  } while (CASADI_PREFIX(flip)(corner, ndim));
  return ret;
}

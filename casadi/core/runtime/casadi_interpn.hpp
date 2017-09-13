// NOLINT(legal/copyright)
// SYMBOL "interpn"
template<typename T1>
T1 casadi_interpn(int ndim, const T1* grid, const int* offset, const T1* values, const T1* x, const int* lookup_mode, int* iw, T1* w) { // NOLINT(whitespace/line_length)
  // Work vectors
  T1* alpha = w; w += ndim;
  int* index = iw; iw += ndim;
  int* corner = iw; iw += ndim;
  // Left index and fraction of interval
  casadi_interpn_weights(ndim, grid, offset, x, alpha, index, lookup_mode);
  // Loop over all corners, add contribution to output
  casadi_fill_int(corner, ndim, 0);
  T1 ret = 0;
  do {
    T1* coeff = 0;
    ret += casadi_interpn_interpolate(ndim, offset, values,
      alpha, index, corner, coeff);
  } while (casadi_flip(corner, ndim));
  return ret;
}

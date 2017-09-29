// NOLINT(legal/copyright)
// SYMBOL "interpn_interpolate"
template<typename T1>
T1 casadi_interpn_interpolate(int ndim, const int* offset, const T1* values, const T1* alpha, const int* index, const int* corner, T1* coeff) { // NOLINT(whitespace/line_length)
  // Get weight and value for corner
  T1 c=1;
  int ld=1; // leading dimension
  int i;
  for (i=0; i<ndim; ++i) {
    if (coeff) *coeff++ = c;
    if (corner[i]) {
      c *= alpha[i];
    } else {
      c *= 1-alpha[i];
    }
    values += (index[i]+corner[i])*ld;
    ld *= offset[i+1]-offset[i];
  }
  if (coeff) {
    return *values;
  } else {
    return c**values;
  }
}

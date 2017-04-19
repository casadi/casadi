template<typename real_t>
real_t CASADI_PREFIX(interpn_interpolate)(int ndim, const int* offset, const real_t* values, const real_t* alpha, const int* index, const int* corner, real_t* coeff) {
  /* Get weight and value for corner */
  real_t c=1;
  int ld=1; /* leading dimension */
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

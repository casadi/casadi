// NOLINT(legal/copyright)

// SYMBOL "vfmax"
// elementwise fmin of two vectors
template<typename T1>
void casadi_vector_fmin(casadi_int n, const T1* x, const T1* y, T1* z) {
  casadi_int i;
  // C-REPLACE "fmax" "casadi_fmax"
  for (i=0; i<n; ++i) z[i] = fmin(x[i], y[i]);
}

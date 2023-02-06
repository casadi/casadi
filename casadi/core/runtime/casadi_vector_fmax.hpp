// NOLINT(legal/copyright)

// SYMBOL "vector_fmax"
// elementwise fmax of two vectors
template<typename T1>
void casadi_vector_fmax(casadi_int n, const T1* x, const T1* y, T1* z) {
  casadi_int i;
  // C-REPLACE "fmax" "casadi_fmax"
  for (i=0; i<n; ++i) z[i] = fmax(x[i], y[i]);
}

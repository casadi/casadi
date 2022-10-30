// NOLINT(legal/copyright)

// SYMBOL "vfmax"
// elementwise fmin of two vectors
template<typename T1>
void casadi_elem_vfmin(casadi_int n, const T1* x, T1*y) {
  casadi_int i;
  // C-REPLACE "fmax" "casadi_fmax"
  for (i=0; i<n; ++i) y[i] = fmin(x[i], y[i]);
}

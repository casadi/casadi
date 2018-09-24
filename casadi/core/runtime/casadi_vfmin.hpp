// NOLINT(legal/copyright)

// SYMBOL "vfmin"
template<typename T1>
T1 casadi_vfmin(const T1* x, casadi_int n, T1 r) {
  casadi_int i;
  // C-REPLACE "fmin" "casadi_fmin"
  for (i=0; i<n; ++i) r = fmin(r, x[i]);
  return r;
}

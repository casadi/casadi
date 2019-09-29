// NOLINT(legal/copyright)

// SYMBOL "vfmax"
template<typename T1>
T1 casadi_vfmax(const T1* x, casadi_int n, T1 r) {
  casadi_int i;
  // C-REPLACE "fmax" "casadi_fmax"
  for (i=0; i<n; ++i) r = fmax(r, x[i]);
  return r;
}

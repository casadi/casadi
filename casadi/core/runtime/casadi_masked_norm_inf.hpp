// NOLINT(legal/copyright)
// SYMBOL "masked_norm_inf"
template<typename T1>
T1 casadi_masked_norm_inf(casadi_int n, const T1* x, const casadi_int* mask) {
  casadi_int i;
  T1 ret = 0;
// C-REPLACE "fmax" "casadi_fmax"
  for (i=0; i<n; ++i) {
    if (mask[i]) {
      ret = fmax(ret, fabs(x[i]));
    }
  }

  return ret;
}

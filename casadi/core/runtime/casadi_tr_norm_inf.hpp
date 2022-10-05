// NOLINT(legal/copyright)
// SYMBOL "norm_inf"
template<typename T1>
T1 casadi_tr_norm_inf(casadi_int n, const T1* x, const T1* binary_vector) {
  casadi_int i;
  T1 ret = 0;
// C-REPLACE "fmax" "casadi_fmax"
  for (i=0; i<n; ++i){
    if (binary_vector[i] != 0){
      ret = fmax(ret, fabs(x[i]));
    }
  }
    
  return ret;
}

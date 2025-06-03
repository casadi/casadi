// NOLINT(legal/copyright)
// SYMBOL "relayout"
template<typename T1>
void casadi_relayout(const T1* arg, T1* res, const casadi_int* source, const casadi_int* perm, const casadi_int* target, casadi_int* iw) {
  casadi_int nnz, i, j, k, kk;
  casadi_int *counter, *target_strides_perm;
  const casadi_int *dims, *target_dims, *strides, *target_strides;
  casadi_int n_dims = source[1];
  nnz = source[2];
  dims = source+3;
  target_dims = target+3;
  strides = dims+2*n_dims;
  target_strides = target_dims+2*n_dims;
  
  counter = iw;
  target_strides_perm = iw+n_dims;
  for (j=0;j<n_dims;++j) {
    target_strides_perm[perm[j]] = target_strides[j];
    counter[j] = 0;
  }
  for (i=0;i<nnz;++i) {
    k = 0;
    for (j=0;j<n_dims;++j) {
      k += strides[j]*counter[j];
    }
    kk = 0;
    for (j=0;j<n_dims;++j) {
      kk += target_strides_perm[j]*counter[j];
    }
    res[kk] = arg[k];
    for (j=0;j<n_dims;++j) {
      counter[j]++;
      if (counter[j]<dims[j]) break;
      counter[j]=0;
    }
  }
}
// NOLINT(legal/copyright)
template<typename T1>
mxArray* CASADI_PREFIX(to_mex)(const int* sp, const T1* x) {
  int nrow = *sp++, ncol = *sp++, nnz = sp[ncol];
  mxArray* p = mxCreateSparse(nrow, ncol, nnz, mxREAL);
  int i;
  mwIndex* j;
  for (i=0, j=mxGetJc(p); i<=ncol; ++i) *j++ = *sp++;
  for (i=0, j=mxGetIr(p); i<nnz; ++i) *j++ = *sp++;
  if (x) {
    double* d = (double*)mxGetData(p);
    for (i=0; i<nnz; ++i) *d++ = to_double(*x++);
  }
  return p;
}

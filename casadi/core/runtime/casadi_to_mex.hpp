// NOLINT(legal/copyright)
template<typename T1>
mxArray* CASADI_PREFIX(to_mex)(const int* sp, const T1* x) {
  int nrow = sp[0], ncol = sp[1], nnz = sp[ncol+2];

  int dense = 0;

#ifdef CASASI_MEX_ALLOW_DENSE
  dense = nrow*ncol==nnz;
#endif
#ifdef CASASI_MEX_ALWAYS_DENSE
  dense = 1;
#endif

  mxArray* p;
  if (dense) {
    p = mxCreateDoubleMatrix(nrow, ncol, mxREAL);
    densify(x, sp, (double*)mxGetData(p), 0);
  } else {
    p = mxCreateSparse(nrow, ncol, nnz, mxREAL);
    int i;
    mwIndex* j;
    sp += 2;
    for (i=0, j=mxGetJc(p); i<=ncol; ++i) *j++ = *sp++;
    for (i=0, j=mxGetIr(p); i<nnz; ++i) *j++ = *sp++;
    if (x) {
      double* d = (double*)mxGetData(p);
      for (i=0; i<nnz; ++i) *d++ = to_double(*x++);
    }
  }
  return p;
}

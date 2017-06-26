// NOLINT(legal/copyright)
template<typename T1>
mxArray* CASADI_PREFIX(to_mex)(const int* sp, const T1* x) {
  int nrow = *sp++, ncol = *sp++, nnz = sp[ncol];
  const int *colind = sp, *row = sp+ncol+1;
#ifndef CASADI_MEX_NO_SPARSE
  if (nnz!=nrow*ncol) {
    mxArray*p = mxCreateSparse(nrow, ncol, nnz, mxREAL);
    int i;
    mwIndex* j;
    for (i=0, j=mxGetJc(p); i<=ncol; ++i) *j++ = *colind++;
    for (i=0, j=mxGetIr(p); i<nnz; ++i) *j++ = *row++;
    if (x) {
      double* d = (double*)mxGetData(p);
      for (i=0; i<nnz; ++i) *d++ = to_double(*x++);
    }
    return p;
  }
#endif /* CASADI_MEX_NO_SPARSE */
  mxArray* p = mxCreateDoubleMatrix(nrow, ncol, mxREAL);
  if (x) {
    double* d = (double*)mxGetData(p);
    int c, k;
    for (c=0; c<ncol; ++c) {
      for (k=colind[c]; k<colind[c+1]; ++k) {
        d[row[k]+c*nrow] = to_double(*x++);
      }
    }
  }
  return p;
}

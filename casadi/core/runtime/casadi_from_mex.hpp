// NOLINT(legal/copyright)
template<typename T1>
T1* CASADI_PREFIX(from_mex)(const mxArray* p, T1* y, const int* sp, T1* w) {
  if (!mxIsDouble(p) || mxGetNumberOfDimensions(p)!=2)
    mexErrMsgIdAndTxt("Casadi:RuntimeError","\"from_mex\" failed: Not a two-dimensional matrix of double precision.");
  size_t nrow = *sp++;
  size_t ncol = *sp++;
  size_t nnz = sp[ncol];
  const int *colind=sp, *row=sp+ncol+1;
  size_t p_nrow = mxGetM(p), p_ncol = mxGetN(p);
  const double* p_data = (const double*)mxGetData(p);
  int is_sparse = mxIsSparse(p);
  mwIndex *Jc = is_sparse ? mxGetJc(p) : 0;
  mwIndex *Ir = is_sparse ? mxGetIr(p) : 0;
  if (p_nrow==1 && p_ncol==1) {
    double v = is_sparse && Jc[1]==0 ? 0 : *p_data;
    fill(y, nnz, v);
  } else {
    int tr = false;
    if (nrow!=p_nrow || ncol!=p_ncol) {
      tr = nrow==p_ncol && ncol==p_nrow && (nrow==1 || ncol==1);
      if (!tr) mexErrMsgIdAndTxt("Casadi:RuntimeError","\"from_mex\" failed: Dimension mismatch.");
    }
    size_t c,k;
    if (is_sparse) {
      if (tr) {
        for (c=0; c<ncol; ++c)
          for (k=(size_t) colind[c]; k<(size_t) colind[c+1]; ++k) w[row[k]+c*nrow]=0;
        for (c=0; c<p_ncol; ++c)
          for (k=Jc[c]; k<Jc[c+1]; ++k) w[c+Ir[k]*p_ncol] = p_data[k];
        for (c=0; c<ncol; ++c)
          for (k=(size_t) colind[c]; k<(size_t) colind[c+1]; ++k) y[k] = w[row[k]+c*nrow];
      } else {
        for (c=0; c<ncol; ++c) {
          for (k=(size_t) colind[c]; k<(size_t) colind[c+1]; ++k) w[row[k]]=0;
          for (k=Jc[c]; k<Jc[c+1]; ++k) w[Ir[k]]=p_data[k];
          for (k=(size_t) colind[c]; k<(size_t) colind[c+1]; ++k) y[k]=w[row[k]];
        }
      }
    } else {
      for (c=0; c<ncol; ++c) {
        for (k=(size_t) colind[c]; k<(size_t) colind[c+1]; ++k) {
          y[k] = p_data[row[k]+c*nrow];
        }
      }
    }
  }
  return y;
}

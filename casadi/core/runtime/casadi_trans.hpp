// NOLINT(legal/copyright)
// SYMBOL "trans"
template<typename T1>
void casadi_trans(const T1* x, const int* sp_x, T1* y, const int* sp_y, int* tmp) {
  int ncol_x = sp_x[1];
  int nnz_x = sp_x[2 + ncol_x];
  const int* row_x = sp_x + 2 + ncol_x+1;
  int ncol_y = sp_y[1];
  const int* colind_y = sp_y+2;
  int k;
  for (k=0; k<ncol_y; ++k) tmp[k] = colind_y[k];
  for (k=0; k<nnz_x; ++k) {
    y[tmp[row_x[k]]++] = x[k];
  }
}

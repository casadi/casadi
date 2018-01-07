// NOLINT(legal/copyright)
// SYMBOL "trans"
template<typename T1>
void casadi_trans(const T1* x, const casadi_int* sp_x, T1* y,
    const casadi_int* sp_y, casadi_int* tmp) {
  casadi_int ncol_x = sp_x[1];
  casadi_int nnz_x = sp_x[2 + ncol_x];
  const casadi_int* row_x = sp_x + 2 + ncol_x+1;
  casadi_int ncol_y = sp_y[1];
  const casadi_int* colind_y = sp_y+2;
  casadi_int k;
  for (k=0; k<ncol_y; ++k) tmp[k] = colind_y[k];
  for (k=0; k<nnz_x; ++k) {
    y[tmp[row_x[k]]++] = x[k];
  }
}

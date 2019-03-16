// NOLINT(legal/copyright)
// SYMBOL "tri_project"
template<typename T1>
void casadi_tri_project(const T1* x, const casadi_int* sp_x, T1* y, casadi_int lower) {
  casadi_int ncol_x, j, k;
  const casadi_int *colind_x, *row_x;
  ncol_x = sp_x[1];
  colind_x = sp_x+2; row_x = sp_x + 2 + ncol_x+1;

  for (j=0; j<ncol_x; ++j) {
    for (k=colind_x[j]; k<colind_x[j+1]; ++k) {
      if (lower) {
        if (row_x[k]>=j) *y++ = x ? x[k] : 0;
      } else {
        if (row_x[k]<=j) *y++ = x ? x[k] : 0;
      }
    }
  }
}

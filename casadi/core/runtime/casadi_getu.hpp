// NOLINT(legal/copyright)
// SYMBOL "getu"
template<typename T1>
void casadi_getu(const T1* x, const casadi_int* sp_x, T1* v) {
  // Get sparsities
  casadi_int ncol_x = sp_x[1];
  const casadi_int *colind_x = sp_x+2, *row_x = sp_x + 2 + ncol_x+1;
  // Loop over the columns of x
  casadi_int cc, el;
  for (cc=0; cc<ncol_x; ++cc) {
    // Loop over the nonzeros of x
    for (el=colind_x[cc]; el<colind_x[cc+1] && row_x[el]<=cc; ++el) {
      *v++ = x[el];
    }
  }
}

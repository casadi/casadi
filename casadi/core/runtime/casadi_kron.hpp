// NOLINT(legal/copyright)
// SYMBOL "kron"
template<typename T1>
void casadi_kron(const T1* a, const casadi_int* sp_a, const T1* b, const casadi_int* sp_b, T1* r) {
    casadi_int a_ncol, b_ncol, k;
    const casadi_int *a_colind, *b_colind;
    T1 a_v, b_v;

    k = 0;

    a_ncol = sp_a[1];
    a_colind = sp_a+2;
    b_ncol = sp_b[1];
    b_colind = sp_b+2;

    // Loop over the columns
    for (casadi_int a_cc=0; a_cc<a_ncol; ++a_cc) {
      // Loop over the columns
      for (casadi_int b_cc=0; b_cc<b_ncol; ++b_cc) {
        // Loop over existing nonzeros
        for (casadi_int a_el=a_colind[a_cc]; a_el<a_colind[a_cc+1]; ++a_el) {
          a_v = a[a_el];
          // Loop over existing nonzeros
          for (casadi_int b_el=b_colind[b_cc]; b_el<b_colind[b_cc+1]; ++b_el) {
            b_v = b[b_el];
            r[k++] = a_v*b_v;
          }
        }
      }
    }

}

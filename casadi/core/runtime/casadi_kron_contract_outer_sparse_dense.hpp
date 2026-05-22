//
//    MIT No Attribution
//
//    Copyright (C) 2010-2023 Joel Andersson, Joris Gillis, Moritz Diehl, KU Leuven.
//
//    Permission is hereby granted, free of charge, to any person obtaining a copy of this
//    software and associated documentation files (the "Software"), to deal in the Software
//    without restriction, including without limitation the rights to use, copy, modify,
//    merge, publish, distribute, sublicense, and/or sell copies of the Software, and to
//    permit persons to whom the Software is furnished to do so.
//
//    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
//    INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
//    PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
//    HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
//    OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
//    SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
//

// SYMBOL "kron_contract_outer_sparse_dense"
// y[r,s] = sum_{i,j} a[i,j] * m[i*mB+r, j*nB+s]
// where m has CSC sparsity sp_m (mA*mB, nA*nB), a is dense (mA, nA),
// y has CSC sparsity sp_y (mB, nB). Walks m's nonzeros only.
template<typename T1>
void casadi_kron_contract_outer_sparse_dense(const T1* m, const casadi_int* sp_m,
                                             const T1* a, casadi_int mA, casadi_int nA,
                                             T1* y, const casadi_int* sp_y) {
    casadi_int nB = sp_y[1];
    const casadi_int* y_colind = sp_y+2;
    const casadi_int* y_row    = sp_y+2+nB+1;
    casadi_int m_ncol = sp_m[1];
    casadi_int mB = sp_y[0];
    const casadi_int* m_colind = sp_m+2;
    const casadi_int* m_row    = sp_m+2+m_ncol+1;
    casadi_int cc, j, s, el, rr, i, r, y_el, y_col_start, y_col_end, k;
    T1 a_val;
    for (k=0; k<y_colind[nB]; ++k) y[k] = 0;
    // Nested (j, s) loop avoids per-column cc / nB and cc % nB.
    for (j=0; j<nA; ++j) {
      for (s=0; s<nB; ++s) {
        cc = j*nB + s;
        y_col_start = y_colind[s];
        y_col_end = y_colind[s+1];
        if (y_col_start == y_col_end) continue;
        for (el=m_colind[cc]; el<m_colind[cc+1]; ++el) {
          rr = m_row[el];
          i = rr / mB;
          r = rr % mB;
          a_val = a[j*mA + i];
          for (y_el=y_col_start; y_el<y_col_end; ++y_el) {
            if (y_row[y_el] == r) { y[y_el] += a_val * m[el]; break; }
            if (y_row[y_el] > r) break;
          }
        }
      }
    }
}

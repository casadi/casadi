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

// SYMBOL "kron_contract_outer_dense_sparse"
// y[r,s] = sum over (i,j) in a's nonzeros of a[i,j] * m[i*mB+r, j*nB+s]
// where m is dense (mA*mB, nA*nB), a is sparse (sp_a, mA x nA), y is dense (mB, nB).
template<typename T1>
void casadi_kron_contract_outer_dense_sparse(const T1* m, casadi_int mB, casadi_int nB,
                                             const T1* a, const casadi_int* sp_a,
                                             T1* y) {
    casadi_int mA = sp_a[0], nA = sp_a[1];
    const casadi_int* a_colind = sp_a+2;
    const casadi_int* a_row    = sp_a+2+nA+1;
    casadi_int j, a_el, i, s, rr, k;
    T1 a_val;
    const T1* m_col, *m_block;
    T1* y_col;
    for (k=0; k<mB*nB; ++k) y[k] = 0;
    for (j=0; j<nA; ++j) {
      for (a_el=a_colind[j]; a_el<a_colind[j+1]; ++a_el) {
        i = a_row[a_el];
        a_val = a[a_el];
        for (s=0; s<nB; ++s) {
          m_col = m + (j*nB+s)*(mA*mB);
          m_block = m_col + i*mB;
          y_col = y + s*mB;
          for (rr=0; rr<mB; ++rr) y_col[rr] += a_val * m_block[rr];
        }
      }
    }
}

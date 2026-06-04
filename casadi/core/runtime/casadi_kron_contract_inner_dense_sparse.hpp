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

// SYMBOL "kron_contract_inner_dense_sparse"
// y[i,j] = sum over (r,s) in b's nonzeros of m[i*mB+r, j*nB+s] * b[r,s]
// where m is dense (mA*mB, nA*nB), b is sparse (sp_b, mB x nB), y is dense (mA, nA).
template<typename T1>
void casadi_kron_contract_inner_dense_sparse(const T1* m, casadi_int mA, casadi_int nA,
                                             const T1* b, const casadi_int* sp_b,
                                             T1* y) {
    casadi_int mB = sp_b[0], nB = sp_b[1];
    const casadi_int* b_colind = sp_b+2;
    const casadi_int* b_row    = sp_b+2+nB+1;
    casadi_int s, b_el, rr, j, i, k;
    T1 b_val;
    const T1* m_col;
    for (k=0; k<mA*nA; ++k) y[k] = 0;
    for (s=0; s<nB; ++s) {
      for (b_el=b_colind[s]; b_el<b_colind[s+1]; ++b_el) {
        rr = b_row[b_el];
        b_val = b[b_el];
        for (j=0; j<nA; ++j) {
          m_col = m + (j*nB+s)*(mA*mB);
          for (i=0; i<mA; ++i) y[j*mA + i] += b_val * m_col[i*mB + rr];
        }
      }
    }
}

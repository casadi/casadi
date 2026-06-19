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

// SYMBOL "kron_contract_outer_dense"
// y[r,s] = sum_{i,j} a[i,j] * m[i*mB+r, j*nB+s] where m, a, y are all dense
// column-major.
template<typename T1>
void casadi_kron_contract_outer_dense(const T1* m, casadi_int mB, casadi_int nB,
                                      const T1* a, casadi_int mA, casadi_int nA,
                                      T1* y) {
    casadi_int j, s, i, rr, k;
    const T1 *m_col, *m_block;
    T1 *y_col;
    T1 a_val;
    for (k=0; k<mB*nB; ++k) y[k] = 0;
    for (j=0; j<nA; ++j) {
      for (s=0; s<nB; ++s) {
        m_col = m + (j*nB+s)*(mA*mB);
        y_col = y + s*mB;
        for (i=0; i<mA; ++i) {
          a_val = a[j*mA + i];
          m_block = m_col + i*mB;
          for (rr=0; rr<mB; ++rr) y_col[rr] += a_val * m_block[rr];
        }
      }
    }
}

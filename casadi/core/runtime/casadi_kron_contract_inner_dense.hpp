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

// SYMBOL "kron_contract_inner_dense"
// y[i,j] = sum_{r,s} m[i*mB+r, j*nB+s] * b[r,s] where m, b, y are all dense
// column-major.
template<typename T1>
void casadi_kron_contract_inner_dense(const T1* m, casadi_int mA, casadi_int nA,
                                      const T1* b, casadi_int mB, casadi_int nB,
                                      T1* y) {
    casadi_int j, s, i, rr, k;
    const T1 *m_col, *m_block, *b_col;
    T1 acc;
    for (k=0; k<mA*nA; ++k) y[k] = 0;
    for (j=0; j<nA; ++j) {
      for (s=0; s<nB; ++s) {
        m_col = m + (j*nB+s)*(mA*mB);
        b_col = b + s*mB;
        for (i=0; i<mA; ++i) {
          m_block = m_col + i*mB;
          acc = 0;
          for (rr=0; rr<mB; ++rr) acc += m_block[rr] * b_col[rr];
          y[j*mA + i] += acc;
        }
      }
    }
}

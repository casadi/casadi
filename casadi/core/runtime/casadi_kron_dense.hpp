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

// SYMBOL "kron_dense"
// r = kron(a, b) where a is (mA, nA) dense column-major, b is (mB, nB) dense
// column-major. r is (mA*mB, nA*nB) dense column-major.
template<typename T1>
void casadi_kron_dense(const T1* a, casadi_int mA, casadi_int nA,
                       const T1* b, casadi_int mB, casadi_int nB,
                       T1* r) {
    casadi_int j, s, i, rr, k;
    T1 a_val;
    k = 0;
    for (j=0; j<nA; ++j) {
      for (s=0; s<nB; ++s) {
        for (i=0; i<mA; ++i) {
          a_val = a[j*mA + i];
          for (rr=0; rr<mB; ++rr) {
            r[k++] = a_val * b[s*mB + rr];
          }
        }
      }
    }
}

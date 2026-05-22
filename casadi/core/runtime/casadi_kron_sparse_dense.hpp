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

// SYMBOL "kron_sparse_dense"
// r = kron(a, b) where a has CSC sparsity sp_a (mA, nA),
// b is (mB, nB) dense column-major.
// r is laid out per the kron sparsity Sparsity::kron(sp_a, dense(mB,nB)).
template<typename T1>
void casadi_kron_sparse_dense(const T1* a, const casadi_int* sp_a,
                              const T1* b, casadi_int mB, casadi_int nB,
                              T1* r) {
    casadi_int nA = sp_a[1];
    const casadi_int* a_colind = sp_a+2;
    casadi_int a_cc, b_cc, a_el, rr, k;
    T1 a_val;
    k = 0;
    for (a_cc=0; a_cc<nA; ++a_cc) {
      for (b_cc=0; b_cc<nB; ++b_cc) {
        for (a_el=a_colind[a_cc]; a_el<a_colind[a_cc+1]; ++a_el) {
          a_val = a[a_el];
          for (rr=0; rr<mB; ++rr) {
            r[k++] = a_val * b[b_cc*mB + rr];
          }
        }
      }
    }
}

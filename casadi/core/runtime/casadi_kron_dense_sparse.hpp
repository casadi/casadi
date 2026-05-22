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

// SYMBOL "kron_dense_sparse"
// r = kron(a, b) where a is (mA, nA) dense column-major,
// b has CSC sparsity sp_b (mB, nB).
// r is laid out per the kron sparsity Sparsity::kron(dense(mA,nA), sp_b).
template<typename T1>
void casadi_kron_dense_sparse(const T1* a, casadi_int mA, casadi_int nA,
                              const T1* b, const casadi_int* sp_b,
                              T1* r) {
    casadi_int nB = sp_b[1];
    const casadi_int* b_colind = sp_b+2;
    casadi_int j, s, i, b_el, k;
    T1 a_val;
    k = 0;
    for (j=0; j<nA; ++j) {
      for (s=0; s<nB; ++s) {
        for (i=0; i<mA; ++i) {
          a_val = a[j*mA + i];
          for (b_el=b_colind[s]; b_el<b_colind[s+1]; ++b_el) {
            r[k++] = a_val * b[b_el];
          }
        }
      }
    }
}

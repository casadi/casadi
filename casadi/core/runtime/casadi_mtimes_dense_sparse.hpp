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

// SYMBOL "mtimes_dense_sparse"
// z (m x p, dense col-major) += x (m x n, dense col-major) * y (n x p, CCS sparse)
template<typename T1>
void casadi_mtimes_dense_sparse(const T1* x, casadi_int nrow_x,
    const T1* y, const casadi_int* sp_y, T1* z) {
  casadi_int j, kk, i, ncol_y;
  const casadi_int *colind_y, *row_y;
  ncol_y = sp_y[1];
  colind_y = sp_y + 2;
  row_y = sp_y + 2 + ncol_y + 1;
  for (j = 0; j < ncol_y; ++j) {
    T1* zj = z + j * nrow_x;
    for (kk = colind_y[j]; kk < colind_y[j + 1]; ++kk) {
      casadi_int rr = row_y[kk];
      T1 yk = y[kk];
      const T1* x_col = x + rr * nrow_x;
      for (i = 0; i < nrow_x; ++i) zj[i] += x_col[i] * yk;
    }
  }
}

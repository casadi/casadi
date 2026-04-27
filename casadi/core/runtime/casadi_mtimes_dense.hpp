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

// SYMBOL "mtimes_dense"
template<typename T1>
void casadi_mtimes_dense(const T1* x, casadi_int nrow_x, casadi_int ncol_x,
    const T1* y, casadi_int ncol_y, T1* z, casadi_int tr) {
  casadi_int i, j, k;
  if (tr) {
    // z(ncol_x, ncol_y) += x(nrow_x, ncol_x)^T * y(nrow_x, ncol_y)
    for (i=0; i<ncol_y; ++i) {
      for (j=0; j<ncol_x; ++j) {
        for (k=0; k<nrow_x; ++k) {
          z[j + i*ncol_x] += x[k + j*nrow_x] * y[k + i*nrow_x];
        }
      }
    }
  } else {
    // z(nrow_x, ncol_y) += x(nrow_x, ncol_x) * y(ncol_x, ncol_y)
    for (i=0; i<ncol_y; ++i) {
      for (j=0; j<nrow_x; ++j) {
        for (k=0; k<ncol_x; ++k) {
          z[j + i*nrow_x] += x[j + k*nrow_x] * y[k + i*ncol_x];
        }
      }
    }
  }
}

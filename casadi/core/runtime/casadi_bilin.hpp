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

// SYMBOL "bilin"
template<typename T1>
T1 casadi_bilin(const T1* A, const casadi_int* sp_A, const T1* x, const T1* y) {
  casadi_int ncol_A, cc, rr, el;
  const casadi_int *colind_A, *row_A;
  T1 ret;
  // Get sparsities
  ncol_A = sp_A[1];
  colind_A = sp_A+2; row_A = sp_A + 2 + ncol_A+1;
  // Return value
  ret=0;
  // Loop over the columns of A
  for (cc=0; cc<ncol_A; ++cc) {
    // Loop over the nonzeros of A
    for (el=colind_A[cc]; el<colind_A[cc+1]; ++el) {
      // Get the row
      rr=row_A[el];
      // Add contribution
      ret += x[rr]*A[el]*y[cc];
    }
  }
  return ret;
}

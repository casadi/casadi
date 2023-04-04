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

// SYMBOL "mv"
template<typename T1>
void casadi_mv(const T1* x, const casadi_int* sp_x, const T1* y, T1* z, casadi_int tr) {
  casadi_int ncol_x, i, el;
  const casadi_int *colind_x, *row_x;
  if (!x || !y || !z) return;
  // Get sparsities
  ncol_x = sp_x[1];
  colind_x = sp_x+2; row_x = sp_x + 2 + ncol_x+1;
  if (tr) {
    // loop over the columns of x
    for (i=0; i<ncol_x; ++i) {
      // loop over the non-zeros of x
      for (el=colind_x[i]; el<colind_x[i+1]; ++el) {
        z[i] += x[el] * y[row_x[el]];
      }
    }
  } else {
    // loop over the columns of x
    for (i=0; i<ncol_x; ++i) {
      // loop over the non-zeros of x
      for (el=colind_x[i]; el<colind_x[i+1]; ++el) {
        z[row_x[el]] += x[el] * y[i];
      }
    }
  }
}

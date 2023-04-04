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

// SYMBOL "sparsify"
template<typename T1, typename T2>
void casadi_sparsify(const T1* x, T2* y, const casadi_int* sp_y, casadi_int tr) {
  casadi_int nrow_y, ncol_y, i, el;
  const casadi_int *colind_y, *row_y;
  nrow_y = sp_y[0];
  ncol_y = sp_y[1];
  colind_y = sp_y+2; row_y = sp_y+ncol_y+3;
  if (tr) {
    for (i=0; i<ncol_y; ++i) {
      for (el=colind_y[i]; el!=colind_y[i+1]; ++el) {
        *y++ = CASADI_CAST(T2, x[i + row_y[el]*ncol_y]);
      }
    }
  } else {
    for (i=0; i<ncol_y; ++i) {
      for (el=colind_y[i]; el!=colind_y[i+1]; ++el) {
        *y++ = CASADI_CAST(T2, x[row_y[el]]);
      }
      x += nrow_y;
    }
  }
}

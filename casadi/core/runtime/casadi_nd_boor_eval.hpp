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

// SYMBOL "nd_boor_eval"
template<typename T1>
void casadi_nd_boor_eval(T1* ret, casadi_int n_dims, const T1* all_knots, const casadi_int* offset, const casadi_int* all_degree, const casadi_int* strides, const T1* c, casadi_int m, const T1* all_x, const casadi_int* lookup_mode, casadi_int* iw, T1* w) { // NOLINT(whitespace/line_length)
  casadi_int k;
  casadi_int *boor_offset, *starts;
  T1 *all_boor;

  boor_offset = iw; iw+=n_dims+1;
  starts = iw; iw+=n_dims;

  all_boor = w;

  boor_offset[0] = 0;

  for (k=0;k<n_dims;++k) {
    T1 *boor;
    const T1* knots;
    T1 x;
    casadi_int degree, n_knots, n_b, L, start;
    boor = all_boor+boor_offset[k];

    degree = all_degree[k];
    knots = all_knots + offset[k];
    n_knots = offset[k+1]-offset[k];
    n_b = n_knots-degree-1;

    x = all_x[k];
    L = casadi_low(x, knots+degree, n_knots-2*degree, lookup_mode[k]);

    start = L;
    if (start>n_b-degree-1) start = n_b-degree-1;

    starts[k] = start;

    casadi_clear(boor, 2*degree+1);
    if (x>=knots[0] && x<=knots[n_knots-1]) {
      if (x==knots[1]) {
        casadi_fill(boor, degree+1, 1.0);
      } else if (x==knots[n_knots-1]) {
        boor[degree] = 1;
      } else if (knots[L+degree]==x) {
        boor[degree-1] = 1;
      } else {
        boor[degree] = 1;
      }
    }
    casadi_de_boor(x, knots+start, 2*degree+2, degree, boor);
    boor+= degree+1;
    boor_offset[k+1] = boor_offset[k] + degree+1;
  }

  casadi_tensor_ttv(ret, n_dims-1, n_dims, all_boor, boor_offset,
      starts, strides, c, m, 1.0, 0);
}

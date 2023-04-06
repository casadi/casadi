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

// SYMBOL "interpn"
template<typename T1>
void casadi_interpn(T1* res, casadi_int ndim, const T1* grid, const casadi_int* offset, const T1* values, const T1* x, const casadi_int* lookup_mode, casadi_int m, casadi_int* iw, T1* w) { // NOLINT(whitespace/line_length)
  // Work vectors
  T1* alpha;
  casadi_int *index, *corner;
  alpha = w; w += ndim;
  index = iw; iw += ndim;
  corner = iw; iw += ndim;
  // Left index and fraction of interval
  casadi_interpn_weights(ndim, grid, offset, x, alpha, index, lookup_mode);
  // Loop over all corners, add contribution to output
  casadi_clear_casadi_int(corner, ndim);
  casadi_clear(res, m);
  do {
    T1* coeff = 0;
    casadi_interpn_interpolate(res, ndim, offset, values,
      alpha, index, corner, coeff, m);
  } while (casadi_flip(corner, ndim));
}

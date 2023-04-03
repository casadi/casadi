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

// SYMBOL "interpn_weights"
template<typename T1>
void casadi_interpn_weights(casadi_int ndim, const T1* grid, const casadi_int* offset, const T1* x, T1* alpha, casadi_int* index, const casadi_int* lookup_mode) { // NOLINT(whitespace/line_length)
  // Left index and fraction of interval
  casadi_int i;
  for (i=0; i<ndim; ++i) {
    casadi_int ng, j;
    T1 xi;
    const T1* g;
    // Grid point
    xi = x ? x[i] : 0;
    // Grid
    g = grid + offset[i];
    ng = offset[i+1]-offset[i];
    // Find left index
    j = index[i] = casadi_low(xi, g, ng, lookup_mode[i]);
    // Get interpolation/extrapolation alpha
    alpha[i] = (xi-g[j])/(g[j+1]-g[j]);
  }
}

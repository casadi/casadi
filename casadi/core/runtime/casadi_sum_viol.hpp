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

// SYMBOL "sum_viol"
template<typename T1>
T1 casadi_sum_viol(casadi_int n, const T1* x, const T1* lb, const T1* ub) {
  T1 r;
  casadi_int i;
  const T1 zero = 0;
  r = 0;
  for (i=0; i<n; ++i) {
    T1 x_i, lb_i, ub_i;
    x_i = x ? *x++ : zero;
    lb_i = lb ? *lb++ : zero;
    ub_i = ub ? *ub++ : zero;
    // C-REPLACE "fmax" "casadi_fmax"
    r += fmax(x_i-ub_i, zero);
    r += fmax(lb_i-x_i, zero);
  }
  return r;
}

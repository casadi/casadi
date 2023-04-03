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

// SYMBOL "logsumexp"
template<typename T1>
T1 casadi_logsumexp(const T1* x, casadi_int n) {
  casadi_int i, max_ind;
  T1 max, r;
  if (n==1) return x[0];
  max_ind = 0;
  max = x[0];
  // Determine max, argmax
  for (i=1; i<n; ++i) {
    if (x[i]>x[0]) {
      max = x[i];
      max_ind = i;
    }
  }
  r = 0;
  for (i=0; i<n; ++i) {
    if (i!=max_ind) r += exp(x[i]-max);
  }
  // C-REPLACE "log1p" "casadi_log1p"
  return log1p(r)+max;
}

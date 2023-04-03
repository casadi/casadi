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


// C-REPLACE "fmin" "casadi_fmin"

// SYMBOL "lb_eig"
// Use Gershgorin to finds upper and lower bounds on the eigenvalues
template<typename T1>
T1 casadi_lb_eig(const casadi_int* sp_h, const T1* h) {
  // Local variables
  casadi_int ncol, c, k;
  T1 center, radius;
  const casadi_int *colind, *row;
  // Return value
  T1 lb_eig = 0;
  // Get sparsity
  ncol = sp_h[1];
  colind = sp_h+2; row = sp_h+ncol+3;
  for (c=0; c<ncol; ++c) {
    // Calculate Gershgorin discs
    center = 0;
    radius = 0;
    for (k=colind[c]; k<colind[c+1]; ++k) {
      if (row[k]==c) {
        center = h[k];
      } else {
        radius += fabs(h[k]);
      }
    }
    // Update the eigenvalue estimates
    if (c==0) {
      lb_eig = center - radius;
    } else {
      lb_eig = fmin(lb_eig, center - radius);
    }
  }
  return lb_eig;
}

// SYMBOL "regularize"
// Add a multiple of the identity matrix to the diagonal
template<typename T1>
void casadi_regularize(const casadi_int* sp_h, T1* h, T1 reg) {
  // Local variables
  casadi_int ncol, c, k;
  const casadi_int *colind, *row;
  // Get sparsity
  ncol = sp_h[1];
  colind = sp_h+2; row = sp_h+ncol+3;
  // Shift diagonal entries
  for (c=0; c<ncol; ++c) {
    for (k=colind[c]; k<colind[c+1]; ++k) {
      if (row[k]==c) h[k] += reg;
    }
  }
}

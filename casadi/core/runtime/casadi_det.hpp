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

// SYMBOL "det"
// Determinant of a matrix from its sparse QR factors (as produced by casadi_qr).
// det = det(Q) * det(R) = prod_c (1 - beta(c)*v(:,c)'v(:,c)) * prod_c R(c,c).
// The caller multiplies by the (data-independent) sign of the row/column pivoting
// to recover the determinant of the original, unpermuted matrix.
template<typename T1>
T1 casadi_det(const casadi_int* sp_v, const T1* nz_v,
              const casadi_int* sp_r, const T1* nz_r, const T1* beta) {
  // Local variables
  casadi_int ncol, c, k;
  const casadi_int *v_colind, *r_colind;
  T1 det, vtv;
  // Extract sparsities
  ncol = sp_v[1];
  v_colind = sp_v+2; r_colind = sp_r+2;
  // Product of the diagonal entries of R
  det = 1;
  for (c=0; c<ncol; ++c) det *= nz_r[r_colind[c+1]-1];
  // Determinant of Q (product of Householder reflectors)
  for (c=0; c<ncol; ++c) {
    vtv = 0;
    for (k=v_colind[c]; k<v_colind[c+1]; ++k) vtv += nz_v[k]*nz_v[k];
    det *= 1 - beta[c]*vtv;
  }
  return det;
}

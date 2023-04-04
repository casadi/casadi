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


// SYMBOL "kkt"
template<typename T1>
void casadi_kkt(const casadi_int* sp_kkt, T1* nz_kkt,
    const casadi_int* sp_h, const T1* nz_h,
    const casadi_int* sp_a, const T1* nz_a,
    const T1* S, const T1* D, T1* w, casadi_int* iw) {
  // Local variables
  casadi_int i, k, j, nx, nz;
  const casadi_int *h_colind, *h_row, *a_colind, *a_row,
                   *kkt_colind, *kkt_row;
  // Extract sparsities
  a_row = (a_colind = sp_a + 2) + (nx = sp_a[1]) + 1;
  h_row = (h_colind = sp_h + 2) + nx + 1;
  kkt_row = (kkt_colind = sp_kkt + 2) + (nz = sp_kkt[1]) + 1;
  // Running indices for each row of A
  for (i = nx; i < nz; ++i) iw[i - nx] = kkt_colind[i];
  // Reset w to zero
  casadi_clear(w, nz);
  // Loop over columns of [H + D_x; A]
  for (i=0; i<nx; ++i) {
    // Copy scaled column of H to w
    for (k=h_colind[i]; k<h_colind[i+1]; ++k) {
      j = h_row[k];
      w[j] = nz_h[k] * S[i] * S[j];
    }
    // Copy scaled column of A to w
    for (k=a_colind[i]; k<a_colind[i+1]; ++k) {
      j = a_row[k] + nx;
      w[j] = nz_a[k] * S[i] * S[j];
    }
    // Add D_x to diagonal
    w[i] += D[i];
    // Copy column to KKT
    for (k=kkt_colind[i]; k<kkt_colind[i+1]; ++k) {
      j = kkt_row[k];
      nz_kkt[k] = w[j];
      if (j >= nx) nz_kkt[iw[j - nx]++] = w[j];
    }
    // Zero out w
    for (k=h_colind[i]; k<h_colind[i+1]; ++k) w[h_row[k]] = 0;
    for (k=a_colind[i]; k<a_colind[i+1]; ++k) w[a_row[k] + nx] = 0;
  }
  // Copy -D_g to diagonal
  for (i=nx; i<nz; ++i) {
    nz_kkt[iw[i - nx]++] = -D[i];
  }
}

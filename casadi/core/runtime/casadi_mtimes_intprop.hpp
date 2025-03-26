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

// SYMBOL "mtimes_intprop"
template<typename T1>
void casadi_mtimes_intprop(const T1* xL, const T1* xR, const casadi_int* sp_x,
    const T1* yL, const T1* yR, const casadi_int* sp_y,
    T1* zL, T1* zR, const casadi_int* sp_z, T1* w, casadi_int tr) { // NOLINT(whitespace/line_length)
  casadi_int ncol_x, ncol_y, ncol_z, cc;
  const casadi_int *colind_x, *row_x, *colind_y, *row_y, *colind_z, *row_z;
  T1 *wL, *wR;

  // Get sparsities
  ncol_x = sp_x[1];
  colind_x = sp_x+2; row_x = sp_x + 2 + ncol_x+1;
  ncol_y = sp_y[1];
  colind_y = sp_y+2; row_y = sp_y + 2 + ncol_y+1;
  ncol_z = sp_z[1];
  colind_z = sp_z+2; row_z = sp_z + 2 + ncol_z+1;

  wL = w;
  wR = w + sp_z[0];

  if (tr) {
    // Loop over the columns of y and z
    for (cc=0; cc<ncol_z; ++cc) {
      casadi_int kk;
      // Get the dense column of y
      for (kk=colind_y[cc]; kk<colind_y[cc+1]; ++kk) {
        wL[row_y[kk]] = yL[kk];
        wR[row_y[kk]] = yR[kk];
      }
      // Loop over the nonzeros of z
      for (kk=colind_z[cc]; kk<colind_z[cc+1]; ++kk) {
        casadi_int kk1;
        casadi_int rr = row_z[kk];
        // Loop over corresponding columns of x
        for (kk1=colind_x[rr]; kk1<colind_x[rr+1]; ++kk1) {
          T1 aL = xL[kk1];
          T1 bL = wL[row_x[kk1]];
          T1 aR = xR[kk1];
          T1 bR = wR[row_x[kk1]];
          zL[kk] += fmin(fmin(fmin(aL*bL, aL*bR), aR*bL), aR*bR);
          zR[kk] += fmax(fmax(fmax(aL*bL, aL*bR), aR*bL), aR*bR);
        }
      }
    }
  } else {
    // Loop over the columns of y and z
    for (cc=0; cc<ncol_y; ++cc) {
      casadi_int kk;
      // Get the dense column of z
      for (kk=colind_z[cc]; kk<colind_z[cc+1]; ++kk) {
        wL[row_z[kk]] = zL[kk];
        wR[row_z[kk]] = zR[kk];
      }
      // Loop over the nonzeros of y
      for (kk=colind_y[cc]; kk<colind_y[cc+1]; ++kk) {
        casadi_int kk1;
        casadi_int rr = row_y[kk];
        // Loop over corresponding columns of x
        for (kk1=colind_x[rr]; kk1<colind_x[rr+1]; ++kk1) {
          T1 aL = xL[kk1];
          T1 bL = yL[kk];
          T1 aR = xR[kk1];
          T1 bR = yR[kk];
          wL[row_x[kk1]] += fmin(fmin(fmin(aL*bL, aL*bR), aR*bL), aR*bR);
          wR[row_x[kk1]] += fmax(fmax(fmax(aL*bL, aL*bR), aR*bL), aR*bR);
        }
      }
      // Get the sparse column of z
      for (kk=colind_z[cc]; kk<colind_z[cc+1]; ++kk) {
        zL[kk] = wL[row_z[kk]];
        zR[kk] = wR[row_z[kk]];
      }
    }
  }
}

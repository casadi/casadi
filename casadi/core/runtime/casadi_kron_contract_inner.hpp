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

// SYMBOL "kron_contract_inner"
// Y[i,j] = sum over (r,s) of M[i*mB+r, j*nB+s] * B[r,s]
// where M has CSC sparsity sp_m of size (mA*mB)x(nA*nB),
//       B has CSC sparsity sp_b of size mB x nB,
//       Y has CSC sparsity sp_y of size mA x nA.
// w is scratch of length mB*nB used to densify b.
template<typename T1>
void casadi_kron_contract_inner(const T1* m, const casadi_int* sp_m,
                     const T1* b, const casadi_int* sp_b,
                     T1* y, const casadi_int* sp_y,
                     T1* w) {
    casadi_int nA = sp_y[1];
    const casadi_int* y_colind = sp_y+2;
    const casadi_int* y_row    = sp_y+2+nA+1;
    casadi_int mB = sp_b[0];
    casadi_int nB = sp_b[1];
    const casadi_int* b_colind = sp_b+2;
    const casadi_int* b_row    = sp_b+2+nB+1;
    casadi_int m_ncol = sp_m[1];
    const casadi_int* m_colind = sp_m+2;
    const casadi_int* m_row    = sp_m+2+m_ncol+1;

    casadi_int k, cc, el, j, s, rr, i, r, y_el, y_col_start, y_col_end;
    T1 b_val;

    // Densify b into w (column-major, mB*nB entries)
    for (k=0; k<mB*nB; ++k) w[k] = 0;
    for (cc=0; cc<nB; ++cc) {
      for (el=b_colind[cc]; el<b_colind[cc+1]; ++el) {
        w[cc*mB + b_row[el]] = b[el];
      }
    }

    // Zero y
    for (k=0; k<y_colind[nA]; ++k) y[k] = 0;

    // Walk M block-column by block-column; (j, s) ARE the column decomposition
    // (no per-column / or %). Row decomposition still uses / and % per nonzero.
    for (j=0; j<nA; ++j) {
      y_col_start = y_colind[j];
      y_col_end   = y_colind[j+1];
      if (y_col_start == y_col_end) continue;
      for (s=0; s<nB; ++s) {
        cc = j*nB + s;
        for (el=m_colind[cc]; el<m_colind[cc+1]; ++el) {
          rr = m_row[el];
          i = rr / mB;
          r = rr % mB;
          b_val = w[s*mB + r];
          // Linear scan y's column j for row i. y_row is sorted ascending.
          for (y_el=y_col_start; y_el<y_col_end; ++y_el) {
            if (y_row[y_el] == i) { y[y_el] += m[el] * b_val; break; }
            if (y_row[y_el] > i) break;
          }
        }
      }
    }
}

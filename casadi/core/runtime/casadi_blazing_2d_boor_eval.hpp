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

// C-REPLACE "casadi_blazing_boor_init<T1>" "casadi_blazing_boor_init"
// C-REPLACE "casadi_blazing_dbasis<T1>" "casadi_blazing_dbasis"
// C-REPLACE "casadi_blazing_d2basis<T1>" "casadi_blazing_d2basis"
// C-REPLACE "casadi_blazing_tensor_ttv2<T1>" "casadi_blazing_tensor_ttv2"

// SYMBOL "blazing_2d_boor_eval"
template<typename T1>
void casadi_blazing_2d_boor_eval(T1* f, T1* J, T1* H, const T1* all_knots, const T1* all_inv, const casadi_int* offset, const T1* c, const T1* dc, const T1* ddc, const T1* all_x, const casadi_int* lookup_mode, casadi_int* iw, T1* w) { // NOLINT(whitespace/line_length)
  casadi_int *starts;
  iw+=2+1;
  starts = iw;

  casadi_int stride1 = offset[1]-offset[0]-4;

  simde__m256d zero = simde_mm256_set1_pd(0.0);

  // Per-dimension de Boor evaluation
  simde__m256d d0[2], d1[2], d2[2];
  const T1* inv2[2] = {0, 0}; const T1* inv3[2] = {0, 0};

  starts[0] = casadi_blazing_boor_init<T1>(all_x[0], all_knots, all_inv,
      offset[0], offset[1], lookup_mode[0], &d0[0], &d1[0], &d2[0], &inv2[0], &inv3[0]);
  starts[1] = casadi_blazing_boor_init<T1>(all_x[1], all_knots, all_inv,
      offset[1], offset[2], lookup_mode[1], &d0[1], &d1[1], &d2[1], &inv2[1], &inv3[1]);

  // Load coefficient sub-tensor (same C used for value and all NPC derivatives)
  simde__m256d C[4];
  for (int j=0;j<4;++j) {
    C[j] = simde_mm256_loadu_pd(c+(starts[1]+j)*stride1+starts[0]);
  }

  if (f) {
    f[0] = casadi_blazing_tensor_ttv2<T1>(C, d0[0], d0[1]);
  }

  if (J || H) {
    const T1* t0 = all_knots + offset[0] + starts[0];
    const T1* t1 = all_knots + offset[1] + starts[1];

    // Effective first-derivative basis: raw for precomputed, summation-by-parts for NPC
    simde__m256d J0 = dc ? d1[0] : casadi_blazing_dbasis<T1>(d1[0], t0, inv3[0]);
    simde__m256d J1 = dc ? d1[1] : casadi_blazing_dbasis<T1>(d1[1], t1, inv3[1]);

    // First derivatives
    if (J) {
      if (dc) {
        stride1 = offset[1]-offset[0]-4-1;
        for (int j=0;j<4;++j) {
          C[j] = simde_mm256_loadu_pd(dc+(starts[1]+j)*stride1+starts[0]-1);
        }
        dc += stride1*(offset[2]-offset[1]-4);
      }
      J[0] = casadi_blazing_tensor_ttv2<T1>(C, J0, d0[1]);

      if (dc) {
        stride1 = offset[1]-offset[0]-4;
        for (int j=0;j<4;++j) {
          if (j==0) {
            C[j] = zero;
          } else {
            C[j] = simde_mm256_loadu_pd(dc+(starts[1]+j-1)*stride1+starts[0]);
          }
        }
      }
      J[1] = casadi_blazing_tensor_ttv2<T1>(C, d0[0], J1);
    }

    // Second derivatives
    if (H) {
      // Effective second-derivative basis
      simde__m256d H0 = ddc ? d2[0] : casadi_blazing_d2basis<T1>(d2[0], t0, inv2[0], inv3[0]);
      simde__m256d H1 = ddc ? d2[1] : casadi_blazing_d2basis<T1>(d2[1], t1, inv2[1], inv3[1]);

      // Diagonal
      if (ddc) {
        stride1 = offset[1]-offset[0]-4-2;
        for (int j=0;j<4;++j)
          C[j] = simde_mm256_loadu_pd(ddc+(starts[1]+j)*stride1+starts[0]-2);
        ddc += stride1*(offset[2]-offset[1]-4);
      }
      H[0] = casadi_blazing_tensor_ttv2<T1>(C, H0, d0[1]);

      if (ddc) {
        stride1 = offset[1]-offset[0]-4;
        for (int j=0;j<4;++j) {
          if (j<=1) {
            C[j] = zero;
          } else {
            C[j] = simde_mm256_loadu_pd(ddc+(starts[1]+j-2)*stride1+starts[0]);
          }
        }
        ddc += stride1*(offset[2]-offset[1]-4-2);
      }
      H[3] = casadi_blazing_tensor_ttv2<T1>(C, d0[0], H1);

      // Off-diagonal
      if (ddc) {
        stride1 = offset[1]-offset[0]-5;
        for (int j=0;j<4;++j) {
          if (j==0) {
            C[j] = zero;
          } else {
            C[j] = simde_mm256_loadu_pd(ddc+(starts[1]+j-1)*stride1+starts[0]-1);
          }
        }
      }
      H[1] = H[2] = casadi_blazing_tensor_ttv2<T1>(C, J0, J1);
    }
  }
}

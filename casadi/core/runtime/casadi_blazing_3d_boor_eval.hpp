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
// C-REPLACE "casadi_blazing_tensor_ttv3<T1>" "casadi_blazing_tensor_ttv3"

// SYMBOL "blazing_3d_boor_eval"
template<typename T1>
void casadi_blazing_3d_boor_eval(T1* f, T1* J, T1* H, const T1* all_knots, const T1* all_inv, const casadi_int* offset, const T1* c, const T1* dc, const T1* ddc, const T1* all_x, const casadi_int* lookup_mode, casadi_int* iw, T1* w) { // NOLINT(whitespace/line_length)
  casadi_int *starts;
  iw+=3+1;
  starts = iw;

  casadi_int stride1 = offset[1]-offset[0]-4;
  casadi_int stride2 = (offset[2]-offset[1]-4)*stride1;

  simde__m256d zero = simde_mm256_set1_pd(0.0);

  // Per-dimension de Boor evaluation
  simde__m256d d0[3], d1[3], d2[3];
  const T1* inv2[3] = {0, 0, 0}; const T1* inv3[3] = {0, 0, 0};
  for (int i = 0; i < 3; ++i) {
    starts[i] = casadi_blazing_boor_init<T1>(all_x[i], all_knots, all_inv,
        offset[i], offset[i+1], lookup_mode[i], &d0[i], &d1[i], &d2[i], &inv2[i], &inv3[i]);
  }

  // Load coefficient sub-tensor (same C used for value and all NPC derivatives)
  simde__m256d C[16];
  for (int j=0;j<4;++j) {
      for (int k=0;k<4;++k) {
          C[j+4*k] = simde_mm256_loadu_pd(c+(starts[1]+j)*stride1+(starts[2]+k)*stride2+starts[0]);
      }
  }

  if (f) {
    f[0] = casadi_blazing_tensor_ttv3<T1>(C, d0[0], d0[1], d0[2]);
  }

  // Jacobian and Hessian share knot pointers and J bases
  if (J || H) {
    const T1* t[3];
    for (int i = 0; i < 3; ++i) {
      t[i] = all_knots + offset[i] + starts[i];
    }

    // Effective first-derivative basis: raw for precomputed, summation-by-parts for NPC
    simde__m256d Ji[3];
    for (int i = 0; i < 3; ++i) {
      Ji[i] = dc ? d1[i] : casadi_blazing_dbasis<T1>(d1[i], t[i], inv3[i]);
    }

    // First derivatives
    if (J) {
      // J[0]: d/dx0
      if (dc) {
        stride1 = offset[1]-offset[0]-4-1;
        stride2 = (offset[2]-offset[1]-4)*stride1;
        for (int j=0;j<4;++j) {
          for (int k=0;k<4;++k) {
            C[j+4*k] = simde_mm256_loadu_pd(
                        dc+(starts[1]+j)*stride1+(starts[2]+k)*stride2+starts[0]-1);
          }
        }
        dc += stride2*(offset[3]-offset[2]-4);
      }
      J[0] = casadi_blazing_tensor_ttv3<T1>(C, Ji[0], d0[1], d0[2]);

      // J[1]: d/dx1
      if (dc) {
        stride1 = offset[1]-offset[0]-4;
        stride2 = (offset[2]-offset[1]-4-1)*stride1;
        for (int j=0;j<4;++j) {
          for (int k=0;k<4;++k) {
            if (j==0) {
              C[j+4*k] = zero;
            } else {
              C[j+4*k] = simde_mm256_loadu_pd(
                          dc+(starts[1]+j-1)*stride1+(starts[2]+k)*stride2+starts[0]);
            }
          }
        }
        dc += stride2*(offset[3]-offset[2]-4);
      }
      J[1] = casadi_blazing_tensor_ttv3<T1>(C, d0[0], Ji[1], d0[2]);

      // J[2]: d/dx2
      if (dc) {
        stride1 = offset[1]-offset[0]-4;
        stride2 = (offset[2]-offset[1]-4)*stride1;
        for (int j=0;j<4;++j) {
          for (int k=0;k<4;++k) {
            if (k==0) {
              C[j+4*k] = zero;
            } else {
              C[j+4*k] = simde_mm256_loadu_pd(
                          dc+(starts[1]+j)*stride1+(starts[2]+k-1)*stride2+starts[0]);
            }
          }
        }
      }
      J[2] = casadi_blazing_tensor_ttv3<T1>(C, d0[0], d0[1], Ji[2]);
    }

    // Second derivatives
    if (H) {
      // Effective second-derivative basis
      simde__m256d Hi[3];
      for (int i = 0; i < 3; ++i) {
        Hi[i] = ddc ? d2[i] : casadi_blazing_d2basis<T1>(d2[i], t[i], inv2[i], inv3[i]);
      }

      // Diagonal: H[i*(n+1)]
      // H[0] = d2/dx0^2
      if (ddc) {
        stride1 = offset[1]-offset[0]-4-2;
        stride2 = (offset[2]-offset[1]-4)*stride1;
        for (int j=0;j<4;++j) {
          for (int k=0;k<4;++k) {
            C[j+4*k] = simde_mm256_loadu_pd(
                        ddc+(starts[1]+j)*stride1+(starts[2]+k)*stride2+starts[0]-2);
          }
        }
        ddc += stride2*(offset[3]-offset[2]-4);
      }
      H[0] = casadi_blazing_tensor_ttv3<T1>(C, Hi[0], d0[1], d0[2]);

      // H[4] = d2/dx1^2
      if (ddc) {
        stride1 = offset[1]-offset[0]-4;
        stride2 = (offset[2]-offset[1]-4-2)*stride1;
        for (int j=0;j<4;++j) {
          for (int k=0;k<4;++k) {
            if (j<=1) {
              C[j+4*k] = zero;
            } else {
              C[j+4*k] = simde_mm256_loadu_pd(
                          ddc+(starts[1]+j-2)*stride1+(starts[2]+k)*stride2+starts[0]);
            }
          }
        }
        ddc += stride2*(offset[3]-offset[2]-4);
      }
      H[4] = casadi_blazing_tensor_ttv3<T1>(C, d0[0], Hi[1], d0[2]);

      // H[8] = d2/dx2^2
      if (ddc) {
        stride1 = offset[1]-offset[0]-4;
        stride2 = (offset[2]-offset[1]-4)*stride1;
        for (int j=0;j<4;++j) {
          for (int k=0;k<4;++k) {
            if (k<=1) {
              C[j+4*k] = zero;
            } else {
              C[j+4*k] = simde_mm256_loadu_pd(
                          ddc+(starts[1]+j)*stride1+(starts[2]+k-2)*stride2+starts[0]);
            }
          }
        }
        ddc += stride2*(offset[3]-offset[2]-4-2);
      }
      H[8] = casadi_blazing_tensor_ttv3<T1>(C, d0[0], d0[1], Hi[2]);

      // Off-diagonal
      // H[1] = H[3] = d2/dx0 dx1
      if (ddc) {
        stride1 = offset[1]-offset[0]-5;
        stride2 = (offset[2]-offset[1]-5)*stride1;
        for (int j=0;j<4;++j) {
          for (int k=0;k<4;++k) {
            if (j==0) {
              C[j+4*k] = zero;
            } else {
              C[j+4*k] = simde_mm256_loadu_pd(
                          ddc+(starts[1]+j-1)*stride1+(starts[2]+k)*stride2+starts[0]-1);
            }
          }
        }
        ddc += stride2*(offset[3]-offset[2]-4);
      }
      H[1] = H[3] = casadi_blazing_tensor_ttv3<T1>(C, Ji[0], Ji[1], d0[2]);

      // H[5] = H[7] = d2/dx1 dx2
      if (ddc) {
        stride1 = offset[1]-offset[0]-4;
        stride2 = (offset[2]-offset[1]-5)*stride1;
        for (int j=0;j<4;++j) {
          for (int k=0;k<4;++k) {
            if (k==0) {
              C[j+4*k] = zero;
            } else {
              C[j+4*k] = simde_mm256_loadu_pd(
                          ddc+(starts[1]+j-1)*stride1+(starts[2]+k-1)*stride2+starts[0]);
            }
          }
        }
        ddc += stride2*(offset[3]-offset[2]-5);
      }
      H[5] = H[7] = casadi_blazing_tensor_ttv3<T1>(C, d0[0], Ji[1], Ji[2]);

      // H[2] = H[6] = d2/dx0 dx2
      if (ddc) {
        stride1 = offset[1]-offset[0]-5;
        stride2 = (offset[2]-offset[1]-4)*stride1;
        for (int j=0;j<4;++j) {
          for (int k=0;k<4;++k) {
            if (k==0) {
              C[j+4*k] = zero;
            } else {
              C[j+4*k] = simde_mm256_loadu_pd(
                          ddc+(starts[1]+j)*stride1+(starts[2]+k-1)*stride2+starts[0]-1);
            }
          }
        }
      }
      H[2] = H[6] = casadi_blazing_tensor_ttv3<T1>(C, Ji[0], d0[1], Ji[2]);
    }
  }
}

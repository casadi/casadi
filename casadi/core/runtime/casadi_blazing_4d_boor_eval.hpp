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
// C-REPLACE "casadi_blazing_tensor_ttv4<T1>" "casadi_blazing_tensor_ttv4"

// SYMBOL "blazing_4d_boor_eval"
template<typename T1>
void casadi_blazing_4d_boor_eval(T1* f, T1* J, T1* H, const T1* all_knots, const T1* all_inv, const casadi_int* offset, const T1* c, const T1* dc, const T1* ddc, const T1* all_x, const casadi_int* lookup_mode, casadi_int* iw, T1* w) { // NOLINT(whitespace/line_length)
  casadi_int *starts;
  iw+=4+1;
  starts = iw;

  casadi_int n_b[4];

  simde__m256d d0[4], d1[4], d2[4];
  const T1* inv2[4] = {0, 0, 0, 0}; const T1* inv3[4] = {0, 0, 0, 0};

  // Per-dimension de Boor evaluation
  for (int i = 0; i < 4; ++i) {
    starts[i] = casadi_blazing_boor_init<T1>(all_x[i], all_knots, all_inv,
        offset[i], offset[i+1], lookup_mode[i], &d0[i], &d1[i], &d2[i], &inv2[i], &inv3[i]);
    n_b[i] = offset[i+1] - offset[i] - 3 - 1;
  }

  // Compute strides and base pointer for 4D coefficient tensor
  casadi_int s1 = n_b[0];
  casadi_int s2 = n_b[0]*n_b[1];
  casadi_int s3 = n_b[0]*n_b[1]*n_b[2];
  const T1* base = c + starts[0] + starts[1]*s1 + starts[2]*s2 + starts[3]*s3;

  if (f) {
    f[0] = casadi_blazing_tensor_ttv4<T1>(base, s1, s2, s3, d0[0], d0[1], d0[2], d0[3]);
  }

  // Jacobian and Hessian share knot pointers and J bases
  if (J || H) {
    const T1* t[4];
    simde__m256d dJ[4];
    for (int i = 0; i < 4; ++i) {
      t[i] = all_knots + offset[i] + starts[i];
      dJ[i] = casadi_blazing_dbasis<T1>(d1[i], t[i], inv3[i]);
    }

    if (J) {
      J[0] = casadi_blazing_tensor_ttv4<T1>(base, s1, s2, s3, dJ[0], d0[1], d0[2], d0[3]);
      J[1] = casadi_blazing_tensor_ttv4<T1>(base, s1, s2, s3, d0[0], dJ[1], d0[2], d0[3]);
      J[2] = casadi_blazing_tensor_ttv4<T1>(base, s1, s2, s3, d0[0], d0[1], dJ[2], d0[3]);
      J[3] = casadi_blazing_tensor_ttv4<T1>(base, s1, s2, s3, d0[0], d0[1], d0[2], dJ[3]);
    }

    if (H) {
      simde__m256d dH[4];
      for (int i = 0; i < 4; ++i)
        dH[i] = casadi_blazing_d2basis<T1>(d2[i], t[i], inv2[i], inv3[i]);

      // Diagonal
      H[0]  = casadi_blazing_tensor_ttv4<T1>(base, s1, s2, s3, dH[0], d0[1], d0[2], d0[3]);
      H[5]  = casadi_blazing_tensor_ttv4<T1>(base, s1, s2, s3, d0[0], dH[1], d0[2], d0[3]);
      H[10] = casadi_blazing_tensor_ttv4<T1>(base, s1, s2, s3, d0[0], d0[1], dH[2], d0[3]);
      H[15] = casadi_blazing_tensor_ttv4<T1>(base, s1, s2, s3, d0[0], d0[1], d0[2], dH[3]);

      // Off-diagonal
      H[1] = H[4]  = casadi_blazing_tensor_ttv4<T1>(base, s1, s2, s3, dJ[0], dJ[1], d0[2], d0[3]);
      H[2] = H[8]  = casadi_blazing_tensor_ttv4<T1>(base, s1, s2, s3, dJ[0], d0[1], dJ[2], d0[3]);
      H[3] = H[12] = casadi_blazing_tensor_ttv4<T1>(base, s1, s2, s3, dJ[0], d0[1], d0[2], dJ[3]);
      H[6] = H[9]  = casadi_blazing_tensor_ttv4<T1>(base, s1, s2, s3, d0[0], dJ[1], dJ[2], d0[3]);
      H[7] = H[13] = casadi_blazing_tensor_ttv4<T1>(base, s1, s2, s3, d0[0], dJ[1], d0[2], dJ[3]);
      H[11] = H[14] = casadi_blazing_tensor_ttv4<T1>(base, s1, s2, s3, d0[0], d0[1], dJ[2], dJ[3]);
    }
  }
}

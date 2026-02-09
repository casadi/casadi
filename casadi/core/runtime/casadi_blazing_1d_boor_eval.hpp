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

// SYMBOL "blazing_1d_boor_eval"
template<typename T1>
void casadi_blazing_1d_boor_eval(T1* f, T1* J, T1* H, const T1* all_knots, const T1* all_inv, const casadi_int* offset, const T1* c, const T1* dc, const T1* ddc, const T1* all_x, const casadi_int* lookup_mode, casadi_int* iw, T1* w) { // NOLINT(whitespace/line_length)
  casadi_int *starts;
  iw+=1+1;
  starts = iw;

  simde__m256d boor0_d2, boor0_d1, boor0_d0;
  const T1* inv2_0 = 0; const T1* inv3_0 = 0;
  starts[0] = casadi_blazing_boor_init<T1>(all_x[0], all_knots, all_inv,
      offset[0], offset[1], lookup_mode[0], &boor0_d0, &boor0_d1, &boor0_d2, &inv2_0, &inv3_0);

  const T1* C = c+starts[0];
  if (f) {
    double boor0_d0v[4];
    simde_mm256_storeu_pd(boor0_d0v, boor0_d0);

    f[0] = 0;
    for (casadi_int i=0;i<4;++i) {
      f[0] += boor0_d0v[i]*C[i];
    }
  }

  // First derivative
  if (J) {
    if (dc) {
      // Precomputed derivative coefficients
      const T1* Cdc = dc+starts[0];

      double boor0_d1v[4];
      simde_mm256_storeu_pd(boor0_d1v, boor0_d1);

      J[0] = 0;
      for (casadi_int i=0;i<3;++i) {
        J[0] += boor0_d1v[i+1]*Cdc[i];
      }
    } else {
      // NPC: summation by parts
      const T1* t0 = all_knots + offset[0] + starts[0];
      simde__m256d boor0_J = casadi_blazing_dbasis<T1>(boor0_d1, t0, inv3_0);

      double boor0_Jv[4];
      simde_mm256_storeu_pd(boor0_Jv, boor0_J);

      J[0] = 0;
      for (casadi_int i=0;i<4;++i) {
        J[0] += boor0_Jv[i]*C[i];
      }
    }
  }

  // Second derivative
  if (H) {
    if (ddc) {
      // Precomputed derivative coefficients
      const T1* Cddc = ddc+starts[0];

      double boor0_d2v[4];
      simde_mm256_storeu_pd(boor0_d2v, boor0_d2);

      H[0] = 0;
      for (casadi_int i=0;i<2;++i) {
        H[0] += boor0_d2v[i+2]*Cddc[i];
      }
    } else {
      // NPC: summation by parts (applied twice)
      const T1* t0 = all_knots + offset[0] + starts[0];
      simde__m256d boor0_H = casadi_blazing_d2basis<T1>(boor0_d2, t0,
          inv2_0, inv3_0);

      double boor0_Hv[4];
      simde_mm256_storeu_pd(boor0_Hv, boor0_H);

      H[0] = 0;
      for (casadi_int i=0;i<4;++i) {
        H[0] += boor0_Hv[i]*C[i];
      }
    }
  }
}

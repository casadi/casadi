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

// SYMBOL "blazing_1d_boor_eval"
template<typename T1>
void casadi_blazing_1d_boor_eval(T1* f, T1* J, T1* H, const T1* all_knots, const casadi_int* offset, const T1* c, const T1* dc, const T1* ddc, const T1* all_x, const casadi_int* lookup_mode, casadi_int* iw, T1* w) { // NOLINT(whitespace/line_length)
  casadi_int n_dims = 1;
  casadi_int m = 1;
  casadi_int n_iter, i, pivot;
  casadi_int *boor_offset, *starts, *index, *coeff_offset;
  T1 *cumprod;
  boor_offset = iw; iw+=n_dims+1;
  starts = iw; iw+=n_dims;
  index = iw; iw+=n_dims;
  coeff_offset = iw;
  cumprod = w; w+= n_dims+1;
  boor_offset[0] = 0;
  cumprod[n_dims] = 1;
  coeff_offset[n_dims] = 0;

  casadi_int stride1 = offset[1]-offset[0]-4;

  simde__m256d zero = simde_mm256_set1_pd(0.0);

  simde__m256d boor_start_0000 = zero;
  simde__m256d boor_start_1111 = simde_mm256_set1_pd(1.0);
  simde__m256d boor_start_0001 = simde_mm256_set_pd(1.0, 0.0, 0.0, 0.0);
  simde__m256d boor_start_0010 = simde_mm256_set_pd(0.0, 1.0, 0.0, 0.0);

  simde__m256d boor0_d3;
  simde__m256d boor0_d2;
  simde__m256d boor0_d1;
  simde__m256d boor0_d0;

    const T1* knots;
    T1 x;
    casadi_int degree, n_knots, n_b, L, start;
    degree = 3;
    knots = all_knots + offset[0];
    n_knots = offset[0+1]-offset[0];
    n_b = n_knots-degree-1;
    x = all_x[0];
    L = casadi_low(x, knots+degree, n_knots-2*degree, lookup_mode[0]);
    start = L;
    if (start>n_b-degree-1) start = n_b-degree-1;
    starts[0] = start;
    boor0_d3 = boor_start_0000;
    if (x>=knots[0] && x<=knots[n_knots-1]) {
      if (x==knots[1]) {
        boor0_d3 = boor_start_1111;
      } else if (x==knots[n_knots-1]) {
        boor0_d3 = boor_start_0001;
      } else if (knots[L+degree]==x) {
        boor0_d3 = boor_start_0010;
      } else {
        boor0_d3 = boor_start_0001;
      }
    }
    casadi_blazing_de_boor(x, knots+start, &boor0_d0, &boor0_d1, &boor0_d2, &boor0_d3);

  double boor0_d0v[4];
  simde_mm256_storeu_pd(boor0_d0v, boor0_d0);

  const T1* C = c+starts[0];
  if (f) {
    f[0] = 0;
    for (casadi_int i=0;i<4;++i) {
      f[0] += boor0_d0v[i]*C[i];
    }
  }

  // First derivative
  if (dc && J) {
    C = dc+starts[0];

    double boor0_d1v[4];
    simde_mm256_storeu_pd(boor0_d1v, boor0_d1);

    J[0] = 0;
    for (casadi_int i=0;i<3;++i) {
      J[0] += boor0_d1v[i+1]*C[i];
    }
  }

  if (ddc && H) {
    C = ddc+starts[0];

    double boor0_d2v[4];
    simde_mm256_storeu_pd(boor0_d2v, boor0_d2);

    H[0] = 0;
    for (casadi_int i=0;i<2;++i) {
      H[0] += boor0_d2v[i+2]*C[i];
    }
  }
}

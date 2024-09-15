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

// SYMBOL "blazing_de_boor_dim3"
template<typename T1>
void casadi_blazing_de_boor_dim3(T1 x, const T1* knots, simde__m256d* boor_d0, simde__m256d* boor_d1, simde__m256d* boor_d2, simde__m256d boor_d3) {
  simde__m256d x_ = simde_mm256_set1_pd(x);
  simde__m256d zero = simde_mm256_set1_pd(0.0);
  simde__m256d mask_end = simde_mm256_set_pd(0.0, -1.0, -1.0, -1.0);
  

  simde__m256d boor_d3i = boor_d3;

  simde__m256d boor_d3i_1 = simde_mm256_permute4x64_pd(boor_d3i, SIMDE_MM_SHUFFLE(3, 3, 2, 1)); // shift one up
  boor_d3i_1 = simde_mm256_blendv_pd(zero, boor_d3i_1, mask_end);

  simde__m256d knotsi = simde_mm256_loadu_pd(knots);
  simde__m256d knotsi_1 = simde_mm256_loadu_pd(knots+1);
  simde__m256d knotsi_2 = simde_mm256_loadu_pd(knots+2);
  simde__m256d knotsi_3 = simde_mm256_loadu_pd(knots+3);
  simde__m256d knotsi_4 = simde_mm256_loadu_pd(knots+4);

  simde__m256d bottom = simde_mm256_sub_pd(knotsi_1, knotsi); // bottom = knots[i + 1] - knots[i];
  simde__m256d bottom_mask = simde_mm256_cmp_pd(bottom, zero, SIMDE_CMP_EQ_OQ); // if (bottom)

  simde__m256d r = simde_mm256_div_pd(simde_mm256_sub_pd(x_, knotsi), bottom); // (x - knots[i]) / bottom;
  r = simde_mm256_blendv_pd(r, zero, bottom_mask);
  *boor_d2 = simde_mm256_mul_pd(r, boor_d3i);

  *boor_d2 = simde_mm256_blendv_pd(*boor_d2, zero, bottom_mask);

  bottom = simde_mm256_sub_pd(knotsi_2, knotsi_1); // bottom = knots[i + 2] - knots[i + 1];
  bottom_mask = simde_mm256_cmp_pd(bottom, zero, SIMDE_CMP_EQ_OQ);
  r = simde_mm256_div_pd(simde_mm256_sub_pd(knotsi_2, x_), bottom); // (knots[i + 2] - x) / bottom
  r = simde_mm256_blendv_pd(r, zero, bottom_mask);

  *boor_d2 = simde_mm256_fmadd_pd(r, boor_d3i_1, *boor_d2);



  simde__m256d boor_d2i_1 = simde_mm256_permute4x64_pd(*boor_d2, SIMDE_MM_SHUFFLE(3, 3, 2, 1)); // shift one up
  boor_d2i_1 = simde_mm256_blendv_pd(zero, boor_d2i_1, mask_end);

  bottom = simde_mm256_sub_pd(knotsi_2, knotsi); // bottom = knots[i + 2] - knots[i];
  bottom_mask = simde_mm256_cmp_pd(bottom, zero, SIMDE_CMP_EQ_OQ); // if (bottom)

  r = simde_mm256_div_pd(simde_mm256_sub_pd(x_, knotsi), bottom); // (x - knots[i]) / bottom;
  r = simde_mm256_blendv_pd(r, zero, bottom_mask);
  *boor_d1 = simde_mm256_mul_pd(r, *boor_d2);

  *boor_d1 = simde_mm256_blendv_pd(*boor_d1, zero, bottom_mask);

  bottom = simde_mm256_sub_pd(knotsi_3, knotsi_1); // bottom = knots[i + 3] - knots[i + 1];
  bottom_mask = simde_mm256_cmp_pd(bottom, zero, SIMDE_CMP_EQ_OQ);
  r = simde_mm256_div_pd(simde_mm256_sub_pd(knotsi_3, x_), bottom);  // (knots[i + 3] - x) / bottom
  r = simde_mm256_blendv_pd(r, zero, bottom_mask);

  *boor_d1 = simde_mm256_fmadd_pd(r, boor_d2i_1, *boor_d1);


  simde__m256d boor_d1i_1 = simde_mm256_permute4x64_pd(*boor_d1, SIMDE_MM_SHUFFLE(3, 3, 2, 1)); // shift one up
  boor_d1i_1 = simde_mm256_blendv_pd(zero, boor_d1i_1, mask_end);

  bottom = simde_mm256_sub_pd(knotsi_3, knotsi); // bottom = knots[i + 3] - knots[i];
  bottom_mask = simde_mm256_cmp_pd(bottom, zero, SIMDE_CMP_EQ_OQ); // if (bottom)

  r = simde_mm256_div_pd(simde_mm256_sub_pd(x_, knotsi), bottom); // (x - knots[i]) / bottom;
  r = simde_mm256_blendv_pd(r, zero, bottom_mask);
  *boor_d0 = simde_mm256_mul_pd(r, *boor_d1);

  *boor_d0 = simde_mm256_blendv_pd(*boor_d0, zero, bottom_mask);

  bottom = simde_mm256_sub_pd(knotsi_4, knotsi_1); // bottom = knots[i + 4] - knots[i + 1];
  bottom_mask = simde_mm256_cmp_pd(bottom, zero, SIMDE_CMP_EQ_OQ);
  r = simde_mm256_div_pd(simde_mm256_sub_pd(knotsi_4, x_), bottom);  // (knots[i + 4] - x) / bottom
  r = simde_mm256_blendv_pd(r, zero, bottom_mask);

  *boor_d0 = simde_mm256_fmadd_pd(r, boor_d1i_1, *boor_d0);
}

// SYMBOL "blazing_printvec"
template<typename T1>
void casadi_blazing_printvec(simde__m256d e) {
  double elements[4];
  simde_mm256_store_pd(elements, e);
  printf("mm256d: <%.4f %.4f %.4f %.4f>\n", elements[0], elements[1], elements[2], elements[3]);
}

// SYMBOL "blazing_nd_boor_eval_dim3"
template<typename T1>
void casadi_blazing_nd_boor_eval_dim3(T1* f, T1* J, const T1* all_knots, const casadi_int* offset, const T1* c, const T1* dc, const T1* all_x, const casadi_int* lookup_mode, casadi_int* iw, T1* w) {
  casadi_int n_dims = 3;
  casadi_int m = 1;
  casadi_int n_iter, k, i, pivot;
  casadi_int *boor_offset, *starts, *index, *coeff_offset;
  T1 *cumprod, *all_boor;
  boor_offset = iw; iw+=n_dims+1;
  starts = iw; iw+=n_dims;
  index = iw; iw+=n_dims;
  coeff_offset = iw;
  cumprod = w; w+= n_dims+1;
  all_boor = w;
  boor_offset[0] = 0;
  cumprod[n_dims] = 1;
  coeff_offset[n_dims] = 0;

  casadi_int stride1 = offset[1]-offset[0]-4;
  casadi_int stride2 = (offset[2]-offset[1]-4)*stride1;

  simde__m256d zero = simde_mm256_set1_pd(0.0);

  simde__m256d boor_start_0000 = zero;
  simde__m256d boor_start_1111 = simde_mm256_set1_pd(1.0);
  simde__m256d boor_start_0001 = simde_mm256_set_pd(1.0, 0.0, 0.0, 0.0);
  simde__m256d boor_start_0010 = simde_mm256_set_pd(0.0, 1.0, 0.0, 0.0);

  simde__m256d boor0_d3;
  simde__m256d boor0_d2;
  simde__m256d boor0_d1;
  simde__m256d boor0_d0;

  simde__m256d boor1_d3;
  simde__m256d boor1_d2;
  simde__m256d boor1_d1;
  simde__m256d boor1_d0;

  simde__m256d boor2_d3;
  simde__m256d boor2_d2;
  simde__m256d boor2_d1;
  simde__m256d boor2_d0;

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
    casadi_blazing_de_boor_dim3(x, knots+start, &boor0_d0, &boor0_d1, &boor0_d2, boor0_d3);

    knots = all_knots + offset[1];
    n_knots = offset[1+1]-offset[1];
    n_b = n_knots-degree-1;
    x = all_x[1];
    L = casadi_low(x, knots+degree, n_knots-2*degree, lookup_mode[1]);
    start = L;
    if (start>n_b-degree-1) start = n_b-degree-1;
    starts[1] = start;
    boor1_d3 = boor_start_0000;
    if (x>=knots[0] && x<=knots[n_knots-1]) {
      if (x==knots[1]) {
        boor1_d3 = boor_start_1111;
      } else if (x==knots[n_knots-1]) {
        boor1_d3 = boor_start_0001;
      } else if (knots[L+degree]==x) {
        boor1_d3 = boor_start_0010;
      } else {
        boor1_d3 = boor_start_0001;
      }
    }
    casadi_blazing_de_boor_dim3(x, knots+start, &boor1_d0, &boor1_d1, &boor1_d2, boor1_d3);

    knots = all_knots + offset[2];
    n_knots = offset[2+1]-offset[2];
    n_b = n_knots-degree-1;
    x = all_x[2];
    L = casadi_low(x, knots+degree, n_knots-2*degree, lookup_mode[2]);
    start = L;
    if (start>n_b-degree-1) start = n_b-degree-1;
    starts[2] = start;
    boor2_d3 = boor_start_0000;
    if (x>=knots[0] && x<=knots[n_knots-1]) {
      if (x==knots[1]) {
        boor2_d3 = boor_start_1111;
      } else if (x==knots[n_knots-1]) {
        boor2_d3 = boor_start_0001;
      } else if (knots[L+degree]==x) {
        boor2_d3 = boor_start_0010;
      } else {
        boor2_d3 = boor_start_0001;
      }
    }
    casadi_blazing_de_boor_dim3(x, knots+start, &boor2_d0, &boor2_d1, &boor2_d2, boor2_d3);

  simde__m256d C[16];

  for (int j=0;j<4;++j) {
      for (int k=0;k<4;++k) {
          C[j+4*k] = simde_mm256_loadu_pd(c+(starts[1]+j)*stride1+(starts[2]+k)*stride2+starts[0]);
      }
  }

  simde__m256d a, b0, b1, b2, b3, c0, c1, c2, c3, r;
  simde__m256d ab[4], cab[4];
  simde__m128d r0, r1;

  a = boor0_d0;
  b0 = simde_mm256_permute4x64_pd(boor1_d0, SIMDE_MM_SHUFFLE(0, 0, 0, 0));
  b1 = simde_mm256_permute4x64_pd(boor1_d0, SIMDE_MM_SHUFFLE(1, 1, 1, 1));
  b2 = simde_mm256_permute4x64_pd(boor1_d0, SIMDE_MM_SHUFFLE(2, 2, 2, 2));
  b3 = simde_mm256_permute4x64_pd(boor1_d0, SIMDE_MM_SHUFFLE(3, 3, 3, 3));

  c0 = simde_mm256_permute4x64_pd(boor2_d0, SIMDE_MM_SHUFFLE(0, 0, 0, 0));
  c1 = simde_mm256_permute4x64_pd(boor2_d0, SIMDE_MM_SHUFFLE(1, 1, 1, 1));
  c2 = simde_mm256_permute4x64_pd(boor2_d0, SIMDE_MM_SHUFFLE(2, 2, 2, 2));
  c3 = simde_mm256_permute4x64_pd(boor2_d0, SIMDE_MM_SHUFFLE(3, 3, 3, 3));

  // Need to compute sum_abc C_abc A_a B_b C_c

  // Step 1: Outer product a b: A_a B_b 
  ab[0] = simde_mm256_mul_pd(a, b0);
  ab[1] = simde_mm256_mul_pd(a, b1);
  ab[2] = simde_mm256_mul_pd(a, b2);
  ab[3] = simde_mm256_mul_pd(a, b3);

  // Sum over b axis: sum_b C_abc * (A_a B_b)_b 
  // cab <- cab + ab[i]*C[i]
  for (int i=0;i<4;++i) {
    cab[i] = simde_mm256_set1_pd(0);
    cab[i] = simde_mm256_fmadd_pd(ab[0], C[4*i+0], cab[i]);
    cab[i] = simde_mm256_fmadd_pd(ab[1], C[4*i+1], cab[i]);
    cab[i] = simde_mm256_fmadd_pd(ab[2], C[4*i+2], cab[i]);
    cab[i] = simde_mm256_fmadd_pd(ab[3], C[4*i+3], cab[i]);
  }

  if (f) {
    // Reduce over the c direction
    r = simde_mm256_set1_pd(0);
    r = simde_mm256_fmadd_pd(cab[0], c0, r);
    r = simde_mm256_fmadd_pd(cab[1], c1, r);
    r = simde_mm256_fmadd_pd(cab[2], c2, r);
    r = simde_mm256_fmadd_pd(cab[3], c3, r);

    // Sum all r entries
    r0  = simde_mm256_castpd256_pd128(r);
    r1 = simde_mm256_extractf128_pd(r, 1);
    r0  = simde_mm_add_pd(r0, r1);
    f[0] = simde_mm_cvtsd_f64(simde_mm_add_sd(r0, simde_mm_unpackhi_pd(r0, r0)));
  }

  // First derivative
  if (dc) {
    stride1 = offset[1]-offset[0]-4-1;
    stride2 = (offset[2]-offset[1]-4)*stride1;
    for (int j=0;j<4;++j) {
        for (int k=0;k<4;++k) {
            C[j+4*k] = simde_mm256_loadu_pd(dc+(starts[1]+j)*stride1+(starts[2]+k)*stride2+starts[0]-1);
        }
    }
    dc += stride2*(offset[3]-offset[2]-4);

    a = boor0_d1;
    ab[0] = simde_mm256_mul_pd(a, b0);
    ab[1] = simde_mm256_mul_pd(a, b1);
    ab[2] = simde_mm256_mul_pd(a, b2);
    ab[3] = simde_mm256_mul_pd(a, b3);
    // Sum over b axis: sum_b C_abc * (A_a B_b)_b 
    // cab <- cab + ab[i]*C[i]
    for (int i=0;i<4;++i) {
      cab[i] = simde_mm256_set1_pd(0);
      cab[i] = simde_mm256_fmadd_pd(ab[0], C[4*i+0], cab[i]);
      cab[i] = simde_mm256_fmadd_pd(ab[1], C[4*i+1], cab[i]);
      cab[i] = simde_mm256_fmadd_pd(ab[2], C[4*i+2], cab[i]);
      cab[i] = simde_mm256_fmadd_pd(ab[3], C[4*i+3], cab[i]);
    }

    // Reduce over the c direction
    r = simde_mm256_set1_pd(0);
    r = simde_mm256_fmadd_pd(cab[0], c0, r);
    r = simde_mm256_fmadd_pd(cab[1], c1, r);
    r = simde_mm256_fmadd_pd(cab[2], c2, r);
    r = simde_mm256_fmadd_pd(cab[3], c3, r);

    // Sum all r entries
    r0  = simde_mm256_castpd256_pd128(r);
    r1 = simde_mm256_extractf128_pd(r, 1);
    r0  = simde_mm_add_pd(r0, r1);
    J[0] = simde_mm_cvtsd_f64(simde_mm_add_sd(r0, simde_mm_unpackhi_pd(r0, r0)));


    stride1 = offset[1]-offset[0]-4;
    stride2 = (offset[2]-offset[1]-4-1)*stride1;
    for (int j=0;j<4;++j) {
        for (int k=0;k<4;++k) {
          if (j==0) {
            C[j+4*k] = zero;
          } else {
            C[j+4*k] = simde_mm256_loadu_pd(dc+(starts[1]+j-1)*stride1+(starts[2]+k)*stride2+starts[0]);
          }
        }
    }
    dc += stride2*(offset[3]-offset[2]-4);

    a = boor0_d0;

    b0 = simde_mm256_permute4x64_pd(boor1_d1, SIMDE_MM_SHUFFLE(0, 0, 0, 0));
    b1 = simde_mm256_permute4x64_pd(boor1_d1, SIMDE_MM_SHUFFLE(1, 1, 1, 1));
    b2 = simde_mm256_permute4x64_pd(boor1_d1, SIMDE_MM_SHUFFLE(2, 2, 2, 2));
    b3 = simde_mm256_permute4x64_pd(boor1_d1, SIMDE_MM_SHUFFLE(3, 3, 3, 3));

    ab[0] = simde_mm256_mul_pd(a, b0);
    ab[1] = simde_mm256_mul_pd(a, b1);
    ab[2] = simde_mm256_mul_pd(a, b2);
    ab[3] = simde_mm256_mul_pd(a, b3);

    // Sum over b axis: sum_b C_abc * (A_a B_b)_b 
    // cab <- cab + ab[i]*C[i]
    for (int i=0;i<4;++i) {
      cab[i] = simde_mm256_set1_pd(0);
      cab[i] = simde_mm256_fmadd_pd(ab[0], C[4*i+0], cab[i]);
      cab[i] = simde_mm256_fmadd_pd(ab[1], C[4*i+1], cab[i]);
      cab[i] = simde_mm256_fmadd_pd(ab[2], C[4*i+2], cab[i]);
      cab[i] = simde_mm256_fmadd_pd(ab[3], C[4*i+3], cab[i]);
    }

    // Reduce over the c direction
    r = simde_mm256_set1_pd(0);
    r = simde_mm256_fmadd_pd(cab[0], c0, r);
    r = simde_mm256_fmadd_pd(cab[1], c1, r);
    r = simde_mm256_fmadd_pd(cab[2], c2, r);
    r = simde_mm256_fmadd_pd(cab[3], c3, r);

    // Sum all r entries
    r0  = simde_mm256_castpd256_pd128(r);
    r1 = simde_mm256_extractf128_pd(r, 1);
    r0  = simde_mm_add_pd(r0, r1);
    J[1] = simde_mm_cvtsd_f64(simde_mm_add_sd(r0, simde_mm_unpackhi_pd(r0, r0)));

    stride1 = offset[1]-offset[0]-4;
    stride2 = (offset[2]-offset[1]-4)*stride1;
    for (int j=0;j<4;++j) {
        for (int k=0;k<4;++k) {
          if (k==0) {
            C[j+4*k] = zero;
          } else {
            C[j+4*k] = simde_mm256_loadu_pd(dc+(starts[1]+j)*stride1+(starts[2]+k-1)*stride2+starts[0]);
          }
        }
    }

    b0 = simde_mm256_permute4x64_pd(boor1_d0, SIMDE_MM_SHUFFLE(0, 0, 0, 0));
    b1 = simde_mm256_permute4x64_pd(boor1_d0, SIMDE_MM_SHUFFLE(1, 1, 1, 1));
    b2 = simde_mm256_permute4x64_pd(boor1_d0, SIMDE_MM_SHUFFLE(2, 2, 2, 2));
    b3 = simde_mm256_permute4x64_pd(boor1_d0, SIMDE_MM_SHUFFLE(3, 3, 3, 3));

    c0 = simde_mm256_permute4x64_pd(boor2_d1, SIMDE_MM_SHUFFLE(0, 0, 0, 0));
    c1 = simde_mm256_permute4x64_pd(boor2_d1, SIMDE_MM_SHUFFLE(1, 1, 1, 1));
    c2 = simde_mm256_permute4x64_pd(boor2_d1, SIMDE_MM_SHUFFLE(2, 2, 2, 2));
    c3 = simde_mm256_permute4x64_pd(boor2_d1, SIMDE_MM_SHUFFLE(3, 3, 3, 3));
    
    ab[0] = simde_mm256_mul_pd(a, b0);
    ab[1] = simde_mm256_mul_pd(a, b1);
    ab[2] = simde_mm256_mul_pd(a, b2);
    ab[3] = simde_mm256_mul_pd(a, b3);

    // Sum over b axis: sum_b C_abc * (A_a B_b)_b 
    // cab <- cab + ab[i]*C[i]
    for (int i=0;i<4;++i) {
      cab[i] = simde_mm256_set1_pd(0);
      cab[i] = simde_mm256_fmadd_pd(ab[0], C[4*i+0], cab[i]);
      cab[i] = simde_mm256_fmadd_pd(ab[1], C[4*i+1], cab[i]);
      cab[i] = simde_mm256_fmadd_pd(ab[2], C[4*i+2], cab[i]);
      cab[i] = simde_mm256_fmadd_pd(ab[3], C[4*i+3], cab[i]);
    }

    // Reduce over the c direction
    r = simde_mm256_set1_pd(0);
    r = simde_mm256_fmadd_pd(cab[0], c0, r);
    r = simde_mm256_fmadd_pd(cab[1], c1, r);
    r = simde_mm256_fmadd_pd(cab[2], c2, r);
    r = simde_mm256_fmadd_pd(cab[3], c3, r);

    // Sum all r entries
    r0  = simde_mm256_castpd256_pd128(r);
    r1 = simde_mm256_extractf128_pd(r, 1);
    r0  = simde_mm_add_pd(r0, r1);
    J[2] = simde_mm_cvtsd_f64(simde_mm_add_sd(r0, simde_mm_unpackhi_pd(r0, r0)));
    
  }
}
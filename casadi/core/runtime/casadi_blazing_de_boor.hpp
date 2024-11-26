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

// SYMBOL "blazing_printvec"
template<typename T1>
void casadi_blazing_printvec(const simde__m256d* e) {
  double elements[4];
  simde_mm256_storeu_pd(elements, *e);
  printf("mm256d: <%.4f %.4f %.4f %.4f>\n", elements[0], elements[1], elements[2], elements[3]);
}

// SYMBOL "blazing_de_boor"
template<typename T1>
void casadi_blazing_de_boor(T1 x, const T1* knots, simde__m256d* boor_d0, simde__m256d* boor_d1, simde__m256d* boor_d2, const simde__m256d* boor_d3) { // NOLINT(whitespace/line_length)
  simde__m256d x_ = simde_mm256_set1_pd(x);
  simde__m256d zero = simde_mm256_set1_pd(0.0);
  simde__m256d mask_end = simde_mm256_set_pd(0.0, -1.0, -1.0, -1.0);

  // shift one up
  simde__m256d boor_d3i_1 = simde_mm256_permute4x64_pd(*boor_d3, SIMDE_MM_SHUFFLE(3, 3, 2, 1));
  boor_d3i_1 = simde_mm256_blendv_pd(zero, boor_d3i_1, mask_end);

  simde__m256d knotsi = simde_mm256_loadu_pd(knots);
  simde__m256d knotsi_1 = simde_mm256_loadu_pd(knots+1);
  simde__m256d knotsi_2 = simde_mm256_loadu_pd(knots+2);
  simde__m256d knotsi_3 = simde_mm256_loadu_pd(knots+3);
  simde__m256d knotsi_4 = simde_mm256_loadu_pd(knots+4);

  simde__m256d bottom = simde_mm256_sub_pd(knotsi_1, knotsi); // bottom = knots[i + 1] - knots[i];
  simde__m256d bottom_mask = simde_mm256_cmp_pd(bottom, zero, SIMDE_CMP_EQ_OQ); // if (bottom)

  // (x - knots[i]) / bottom;
  simde__m256d r = simde_mm256_div_pd(simde_mm256_sub_pd(x_, knotsi), bottom);
  r = simde_mm256_blendv_pd(r, zero, bottom_mask);
  *boor_d2 = simde_mm256_mul_pd(r, *boor_d3);

  *boor_d2 = simde_mm256_blendv_pd(*boor_d2, zero, bottom_mask);

  bottom = simde_mm256_sub_pd(knotsi_2, knotsi_1); // bottom = knots[i + 2] - knots[i + 1];
  bottom_mask = simde_mm256_cmp_pd(bottom, zero, SIMDE_CMP_EQ_OQ);
  r = simde_mm256_div_pd(simde_mm256_sub_pd(knotsi_2, x_), bottom); // (knots[i + 2] - x) / bottom
  r = simde_mm256_blendv_pd(r, zero, bottom_mask);

  *boor_d2 = simde_mm256_fmadd_pd(r, boor_d3i_1, *boor_d2);

  // shift one up
  simde__m256d boor_d2i_1 = simde_mm256_permute4x64_pd(*boor_d2, SIMDE_MM_SHUFFLE(3, 3, 2, 1));
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

  // shift one up
  simde__m256d boor_d1i_1 = simde_mm256_permute4x64_pd(*boor_d1, SIMDE_MM_SHUFFLE(3, 3, 2, 1));
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

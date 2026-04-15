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

// C-REPLACE "casadi_blazing_boor_der<T1>" "casadi_blazing_boor_der"
// C-REPLACE "casadi_blazing_shift_left<T1>" "casadi_blazing_shift_left"
// C-REPLACE "casadi_blazing_knot_scale<T1>" "casadi_blazing_knot_scale"
// C-REPLACE "casadi_blazing_boor_init<T1>" "casadi_blazing_boor_init"
// C-REPLACE "casadi_blazing_dbasis<T1>" "casadi_blazing_dbasis"
// C-REPLACE "casadi_blazing_d2basis<T1>" "casadi_blazing_d2basis"
// C-REPLACE "casadi_blazing_hsum<T1>" "casadi_blazing_hsum"
// C-REPLACE "casadi_blazing_tensor_ttv2<T1>" "casadi_blazing_tensor_ttv2"
// C-REPLACE "casadi_blazing_tensor_ttv3<T1>" "casadi_blazing_tensor_ttv3"
// C-REPLACE "casadi_blazing_tensor_ttv4<T1>" "casadi_blazing_tensor_ttv4"
// C-REPLACE "casadi_blazing_tensor_ttv5<T1>" "casadi_blazing_tensor_ttv5"

// ===== De Boor evaluation =====

// SYMBOL "blazing_printvec"
template<typename T1>
void casadi_blazing_printvec(const simde__m256d* e) {
  double elements[4];
  simde_mm256_storeu_pd(elements, *e);
  printf("mm256d: <%.4f %.4f %.4f %.4f>\n", elements[0], elements[1], elements[2], elements[3]);
}

// SYMBOL "blazing_de_boor"
template<typename T1>
void casadi_blazing_de_boor(T1 x, const T1* knots, const T1* inv1, const T1* inv2, const T1* inv3, simde__m256d* boor_d0, simde__m256d* boor_d1, simde__m256d* boor_d2, const simde__m256d* boor_d3) { // NOLINT(whitespace/line_length)
  simde__m256d x_ = simde_mm256_set1_pd(x);
  simde__m256d zero = simde_mm256_set1_pd(0.0);
  simde__m256d mask_end = simde_mm256_set_pd(0.0, -1.0, -1.0, -1.0);
  simde__m256d r;

  // shift one up
  simde__m256d boor_d3i_1 = simde_mm256_permute4x64_pd(*boor_d3, SIMDE_MM_SHUFFLE(3, 3, 2, 1));
  boor_d3i_1 = simde_mm256_blendv_pd(zero, boor_d3i_1, mask_end);

  simde__m256d knotsi   = simde_mm256_loadu_pd(knots);
  simde__m256d knotsi_2 = simde_mm256_loadu_pd(knots+2);
  simde__m256d knotsi_3 = simde_mm256_loadu_pd(knots+3);
  simde__m256d knotsi_4 = simde_mm256_loadu_pd(knots+4);

  if (inv1) {
    // ---- Reciprocal path: multiply by precomputed 1/span ----
    // d3 -> d2 (span-1)
    r = simde_mm256_mul_pd(simde_mm256_sub_pd(x_, knotsi), simde_mm256_loadu_pd(inv1));
    *boor_d2 = simde_mm256_mul_pd(r, *boor_d3);
    r = simde_mm256_mul_pd(simde_mm256_sub_pd(knotsi_2, x_), simde_mm256_loadu_pd(inv1 + 1));
    *boor_d2 = simde_mm256_fmadd_pd(r, boor_d3i_1, *boor_d2);

    // shift d2
    simde__m256d boor_d2i_1 = simde_mm256_permute4x64_pd(*boor_d2, SIMDE_MM_SHUFFLE(3, 3, 2, 1));
    boor_d2i_1 = simde_mm256_blendv_pd(zero, boor_d2i_1, mask_end);

    // d2 -> d1 (span-2)
    r = simde_mm256_mul_pd(simde_mm256_sub_pd(x_, knotsi), simde_mm256_loadu_pd(inv2));
    *boor_d1 = simde_mm256_mul_pd(r, *boor_d2);
    r = simde_mm256_mul_pd(simde_mm256_sub_pd(knotsi_3, x_), simde_mm256_loadu_pd(inv2 + 1));
    *boor_d1 = simde_mm256_fmadd_pd(r, boor_d2i_1, *boor_d1);

    // shift d1
    simde__m256d boor_d1i_1 = simde_mm256_permute4x64_pd(*boor_d1, SIMDE_MM_SHUFFLE(3, 3, 2, 1));
    boor_d1i_1 = simde_mm256_blendv_pd(zero, boor_d1i_1, mask_end);

    // d1 -> d0 (span-3)
    r = simde_mm256_mul_pd(simde_mm256_sub_pd(x_, knotsi), simde_mm256_loadu_pd(inv3));
    *boor_d0 = simde_mm256_mul_pd(r, *boor_d1);
    r = simde_mm256_mul_pd(simde_mm256_sub_pd(knotsi_4, x_), simde_mm256_loadu_pd(inv3 + 1));
    *boor_d0 = simde_mm256_fmadd_pd(r, boor_d1i_1, *boor_d0);
  } else {
    // ---- Original division path ----
    simde__m256d knotsi_1 = simde_mm256_loadu_pd(knots+1);
    simde__m256d bottom, bottom_mask;

    bottom = simde_mm256_sub_pd(knotsi_1, knotsi);
    bottom_mask = simde_mm256_cmp_pd(bottom, zero, SIMDE_CMP_EQ_OQ);
    r = simde_mm256_div_pd(simde_mm256_sub_pd(x_, knotsi), bottom);
    r = simde_mm256_blendv_pd(r, zero, bottom_mask);
    *boor_d2 = simde_mm256_mul_pd(r, *boor_d3);
    *boor_d2 = simde_mm256_blendv_pd(*boor_d2, zero, bottom_mask);

    bottom = simde_mm256_sub_pd(knotsi_2, knotsi_1);
    bottom_mask = simde_mm256_cmp_pd(bottom, zero, SIMDE_CMP_EQ_OQ);
    r = simde_mm256_div_pd(simde_mm256_sub_pd(knotsi_2, x_), bottom);
    r = simde_mm256_blendv_pd(r, zero, bottom_mask);
    *boor_d2 = simde_mm256_fmadd_pd(r, boor_d3i_1, *boor_d2);

    simde__m256d boor_d2i_1 = simde_mm256_permute4x64_pd(*boor_d2, SIMDE_MM_SHUFFLE(3, 3, 2, 1));
    boor_d2i_1 = simde_mm256_blendv_pd(zero, boor_d2i_1, mask_end);

    bottom = simde_mm256_sub_pd(knotsi_2, knotsi);
    bottom_mask = simde_mm256_cmp_pd(bottom, zero, SIMDE_CMP_EQ_OQ);
    r = simde_mm256_div_pd(simde_mm256_sub_pd(x_, knotsi), bottom);
    r = simde_mm256_blendv_pd(r, zero, bottom_mask);
    *boor_d1 = simde_mm256_mul_pd(r, *boor_d2);
    *boor_d1 = simde_mm256_blendv_pd(*boor_d1, zero, bottom_mask);

    bottom = simde_mm256_sub_pd(knotsi_3, knotsi_1);
    bottom_mask = simde_mm256_cmp_pd(bottom, zero, SIMDE_CMP_EQ_OQ);
    r = simde_mm256_div_pd(simde_mm256_sub_pd(knotsi_3, x_), bottom);
    r = simde_mm256_blendv_pd(r, zero, bottom_mask);
    *boor_d1 = simde_mm256_fmadd_pd(r, boor_d2i_1, *boor_d1);

    simde__m256d boor_d1i_1 = simde_mm256_permute4x64_pd(*boor_d1, SIMDE_MM_SHUFFLE(3, 3, 2, 1));
    boor_d1i_1 = simde_mm256_blendv_pd(zero, boor_d1i_1, mask_end);

    bottom = simde_mm256_sub_pd(knotsi_3, knotsi);
    bottom_mask = simde_mm256_cmp_pd(bottom, zero, SIMDE_CMP_EQ_OQ);
    r = simde_mm256_div_pd(simde_mm256_sub_pd(x_, knotsi), bottom);
    r = simde_mm256_blendv_pd(r, zero, bottom_mask);
    *boor_d0 = simde_mm256_mul_pd(r, *boor_d1);
    *boor_d0 = simde_mm256_blendv_pd(*boor_d0, zero, bottom_mask);

    bottom = simde_mm256_sub_pd(knotsi_4, knotsi_1);
    bottom_mask = simde_mm256_cmp_pd(bottom, zero, SIMDE_CMP_EQ_OQ);
    r = simde_mm256_div_pd(simde_mm256_sub_pd(knotsi_4, x_), bottom);
    r = simde_mm256_blendv_pd(r, zero, bottom_mask);
    *boor_d0 = simde_mm256_fmadd_pd(r, boor_d1i_1, *boor_d0);
  }
}

// ===== Basis derivative helpers =====

// SYMBOL "blazing_boor_der"
// Summation by parts: convert derivative from coefficient-space to basis-space.
//
// Given boor = [b0, b1, b2, b3] and scale = [s0, s1, s2, s3]:
//   sb[i] = b[i]*s[i]
//   result[i] = sb[i-1] - sb[i]   (with sb[-1]=0)
// i.e. [-s0*b0, s0*b0-s1*b1, s1*b1-s2*b2, s2*b2-s3*b3]
//
// To compute the derivative of f = sum_i c[i]*N_{i,d}(x),
// use boor = shift_left(d_{d-1} basis) and appropriate knot scales.
//
// For cubic Jacobian: boor_J = blazing_boor_der(shift_left(d1), s1)
//   where shift_left([0,N0,N1,N2]) = [N0,N1,N2,0]
//   and s1[j] = 3/(t[j+start+4]-t[j+start+1])
//
// For cubic Hessian: inner = blazing_boor_der(shift_left(d2), s2)
//                    boor_H = blazing_boor_der(shift_left(inner), s1)
//   where s2[j] = 2/(t[j+start+3]-t[j+start+1])
template<typename T1>
simde__m256d casadi_blazing_boor_der(simde__m256d boor, simde__m256d scale) {
    simde__m256d sb = simde_mm256_mul_pd(boor, scale);
    simde__m256d shifted = simde_mm256_permute4x64_pd(sb, SIMDE_MM_SHUFFLE(2, 1, 0, 0));
    shifted = simde_mm256_blend_pd(simde_mm256_setzero_pd(), shifted, 0xE);
    return simde_mm256_sub_pd(shifted, sb);
}

// SYMBOL "blazing_shift_left"
// Shift AVX vector one position to the left, filling position 3 with zero.
// [a, b, c, d] -> [b, c, d, 0]
template<typename T1>
simde__m256d casadi_blazing_shift_left(simde__m256d v) {
    simde__m256d shifted = simde_mm256_permute4x64_pd(v, SIMDE_MM_SHUFFLE(3, 3, 2, 1));
    return simde_mm256_blend_pd(shifted, simde_mm256_setzero_pd(), 0x8);
}

// SYMBOL "blazing_knot_scale"
// Compute knot scale with safe division: degree/(t_hi - t_lo), returning 0 where span is 0.
// This avoids NaN from 0/0 at knot boundaries where both the basis function and
// the knot span are zero.
template<typename T1>
simde__m256d casadi_blazing_knot_scale(simde__m256d degree, simde__m256d t_hi, simde__m256d t_lo) {
    simde__m256d zero = simde_mm256_setzero_pd();
    simde__m256d denom = simde_mm256_sub_pd(t_hi, t_lo);
    simde__m256d denom_mask = simde_mm256_cmp_pd(denom, zero, SIMDE_CMP_EQ_OQ);
    simde__m256d scale = simde_mm256_div_pd(degree, denom);
    return simde_mm256_blendv_pd(scale, zero, denom_mask);
}

// SYMBOL "blazing_boor_init"
// Per-dimension de Boor setup: knot lookup, boundary-case init, and de Boor evaluation.
// Returns the start index and fills d0, d1, d2 with basis function intermediates.
// all_inv can be 0 (NULL) to use the division path.
// Layout (dim-D-K): per-dimension [inv1[n_k], inv2[n_k], inv3[n_k]].
// If inv2_out/inv3_out are non-NULL, they receive pre-positioned derivative
// base pointers that can be passed directly to dbasis/d2basis.
template<typename T1>
casadi_int casadi_blazing_boor_init(
    T1 x, const T1* all_knots, const T1* all_inv,
    casadi_int knot_offset, casadi_int knot_offset_next,
    casadi_int lookup_mode,
    simde__m256d* d0, simde__m256d* d1, simde__m256d* d2,
    const T1** inv2_out, const T1** inv3_out) {
    casadi_int degree = 3;
    const T1* knots = all_knots + knot_offset;
    casadi_int n_knots = knot_offset_next - knot_offset;
    casadi_int n_b = n_knots - degree - 1;
    casadi_int L = casadi_low(x, knots + degree, n_knots - 2*degree, lookup_mode);
    casadi_int start = L;
    if (start > n_b - degree - 1) start = n_b - degree - 1;

    simde__m256d d3 = simde_mm256_setzero_pd();
    if (x >= knots[0] && x <= knots[n_knots-1]) {
      if (x == knots[1]) {
        d3 = simde_mm256_set1_pd(1.0);
      } else if (x == knots[n_knots-1]) {
        d3 = simde_mm256_set_pd(1.0, 0.0, 0.0, 0.0);
      } else if (knots[L+degree] == x) {
        d3 = simde_mm256_set_pd(0.0, 1.0, 0.0, 0.0);
      } else {
        d3 = simde_mm256_set_pd(1.0, 0.0, 0.0, 0.0);
      }
    }
    const T1 *inv1 = 0, *inv2 = 0, *inv3 = 0;
    if (all_inv) {
      const T1* dim_base = all_inv + 3 * knot_offset;
      inv1 = dim_base + start;
      inv2 = dim_base + n_knots + start;
      inv3 = dim_base + 2*n_knots + start;
      if (inv2_out) *inv2_out = dim_base + n_knots + start + 1;
      if (inv3_out) *inv3_out = dim_base + 2*n_knots + start + 1;
    } else {
      if (inv2_out) *inv2_out = 0;
      if (inv3_out) *inv3_out = 0;
    }
    casadi_blazing_de_boor(x, knots + start, inv1, inv2, inv3, d0, d1, d2, &d3);
    return start;
}

// SYMBOL "blazing_dbasis"
// Compute 1st-derivative basis functions for one dimension.
// t points to knots at starts[i] for this dimension.
// inv3 can be 0 (NULL) to use the division path; otherwise 1/(t[k+3]-t[k]).
template<typename T1>
simde__m256d casadi_blazing_dbasis(simde__m256d boor_d1, const T1* t, const T1* inv3) {
    simde__m256d three = simde_mm256_set1_pd(3.0);
    simde__m256d s1;
    if (inv3) {
      s1 = simde_mm256_mul_pd(three, simde_mm256_loadu_pd(inv3));
    } else {
      s1 = casadi_blazing_knot_scale<T1>(three,
          simde_mm256_loadu_pd(t + 4), simde_mm256_loadu_pd(t + 1));
    }
    return casadi_blazing_boor_der<T1>(
        casadi_blazing_shift_left<T1>(boor_d1), s1);
}

// SYMBOL "blazing_d2basis"
// Compute 2nd-derivative basis functions for one dimension.
// t points to knots at starts[i] for this dimension.
// inv2, inv3 can be 0 (NULL) to use the division path;
// otherwise 1/(t[k+2]-t[k]) and 1/(t[k+3]-t[k]).
template<typename T1>
simde__m256d casadi_blazing_d2basis(simde__m256d boor_d2, const T1* t, const T1* inv2, const T1* inv3) {// NOLINT(whitespace/line_length)
    simde__m256d three = simde_mm256_set1_pd(3.0);
    simde__m256d two = simde_mm256_set1_pd(2.0);
    simde__m256d s1, s2;
    if (inv3) {
      s1 = simde_mm256_mul_pd(three, simde_mm256_loadu_pd(inv3));
      s2 = simde_mm256_mul_pd(two, simde_mm256_loadu_pd(inv2));
    } else {
      s1 = casadi_blazing_knot_scale<T1>(three,
          simde_mm256_loadu_pd(t + 4), simde_mm256_loadu_pd(t + 1));
      s2 = casadi_blazing_knot_scale<T1>(two,
          simde_mm256_loadu_pd(t + 3), simde_mm256_loadu_pd(t + 1));
    }
    simde__m256d inner = casadi_blazing_boor_der<T1>(
        casadi_blazing_shift_left<T1>(boor_d2), s2);
    return casadi_blazing_boor_der<T1>(
        casadi_blazing_shift_left<T1>(inner), s1);
}

// ===== Tensor-times-vector contractions =====

// AVX2 horizontal sum: reduce 4-wide __m256d to scalar double
// SYMBOL "blazing_hsum"
template<typename T1>
T1 casadi_blazing_hsum(simde__m256d r) {
  simde__m128d r0 = simde_mm256_castpd256_pd128(r);
  simde__m128d r1 = simde_mm256_extractf128_pd(r, 1);
  r0 = simde_mm_add_pd(r0, r1);
  return simde_mm_cvtsd_f64(simde_mm_add_sd(r0, simde_mm_unpackhi_pd(r0, r0)));
}

// AVX2 tensor-times-vector for 2D cubic B-splines (degree 3, m=1).
//
//   result = sum_{i,j} a[i] * b[j] * C[j][i]
//
// where a[i] are in AVX lanes, b[j] are in AVX lanes (broadcast internally),
// and C[4] holds coefficient vectors (4 doubles along dim 0 each).
//
// SYMBOL "blazing_tensor_ttv2"
template<typename T1>
T1 casadi_blazing_tensor_ttv2(const simde__m256d C[4],
    simde__m256d a, simde__m256d b) {
  simde__m256d r;
  // Broadcast dim 1 weights and form outer product with dim 0
  simde__m256d ab0 = simde_mm256_mul_pd(a,
      simde_mm256_permute4x64_pd(b, SIMDE_MM_SHUFFLE(0, 0, 0, 0)));
  simde__m256d ab1 = simde_mm256_mul_pd(a,
      simde_mm256_permute4x64_pd(b, SIMDE_MM_SHUFFLE(1, 1, 1, 1)));
  simde__m256d ab2 = simde_mm256_mul_pd(a,
      simde_mm256_permute4x64_pd(b, SIMDE_MM_SHUFFLE(2, 2, 2, 2)));
  simde__m256d ab3 = simde_mm256_mul_pd(a,
      simde_mm256_permute4x64_pd(b, SIMDE_MM_SHUFFLE(3, 3, 3, 3)));
  // Contract dim 1: r = sum_j ab[j] * C[j]
  r = simde_mm256_mul_pd(ab0, C[0]);
  r = simde_mm256_fmadd_pd(ab1, C[1], r);
  r = simde_mm256_fmadd_pd(ab2, C[2], r);
  r = simde_mm256_fmadd_pd(ab3, C[3], r);
  // Horizontal sum contracts dim 0
  return casadi_blazing_hsum<T1>(r);
}

// AVX2 tensor-times-vector for 3D cubic B-splines (degree 3, m=1).
//
//   result = sum_{i,j,k} a[i] * b[j] * c[k] * C[j + 4*k][i]
//
// where a[i] are in AVX lanes, b and c hold 4 weights each (broadcast internally),
// and C[16] holds 4x4 coefficient vectors (4 doubles along dim 0 each).
//
// SYMBOL "blazing_tensor_ttv3"
template<typename T1>
T1 casadi_blazing_tensor_ttv3(const simde__m256d C[16],
    simde__m256d a, simde__m256d b, simde__m256d c) {
  simde__m256d ab[4], cab[4], r;
  int i;
  // Broadcast dim 1 weights and form outer product with dim 0
  simde__m256d b0 = simde_mm256_permute4x64_pd(b, SIMDE_MM_SHUFFLE(0, 0, 0, 0));
  simde__m256d b1 = simde_mm256_permute4x64_pd(b, SIMDE_MM_SHUFFLE(1, 1, 1, 1));
  simde__m256d b2 = simde_mm256_permute4x64_pd(b, SIMDE_MM_SHUFFLE(2, 2, 2, 2));
  simde__m256d b3 = simde_mm256_permute4x64_pd(b, SIMDE_MM_SHUFFLE(3, 3, 3, 3));
  ab[0] = simde_mm256_mul_pd(a, b0);
  ab[1] = simde_mm256_mul_pd(a, b1);
  ab[2] = simde_mm256_mul_pd(a, b2);
  ab[3] = simde_mm256_mul_pd(a, b3);
  // Contract dim 1: cab[k] = sum_j ab[j] * C[j + 4*k]
  for (i = 0; i < 4; ++i) {
    cab[i] = simde_mm256_mul_pd(ab[0], C[4*i+0]);
    cab[i] = simde_mm256_fmadd_pd(ab[1], C[4*i+1], cab[i]);
    cab[i] = simde_mm256_fmadd_pd(ab[2], C[4*i+2], cab[i]);
    cab[i] = simde_mm256_fmadd_pd(ab[3], C[4*i+3], cab[i]);
  }
  // Broadcast dim 2 weights and contract: r = sum_k cab[k] * c_k
  r = simde_mm256_mul_pd(cab[0],
      simde_mm256_permute4x64_pd(c, SIMDE_MM_SHUFFLE(0, 0, 0, 0)));
  r = simde_mm256_fmadd_pd(cab[1],
      simde_mm256_permute4x64_pd(c, SIMDE_MM_SHUFFLE(1, 1, 1, 1)), r);
  r = simde_mm256_fmadd_pd(cab[2],
      simde_mm256_permute4x64_pd(c, SIMDE_MM_SHUFFLE(2, 2, 2, 2)), r);
  r = simde_mm256_fmadd_pd(cab[3],
      simde_mm256_permute4x64_pd(c, SIMDE_MM_SHUFFLE(3, 3, 3, 3)), r);
  // Horizontal sum contracts dim 0
  return casadi_blazing_hsum<T1>(r);
}

// AVX2 tensor-times-vector for 4D cubic B-splines (degree 3, m=1).
//
//   result = sum_{i,j,k,l} a[i] * b[j] * c[k] * d[l]
//            * coeffs[i + s1*j + s2*k + s3*l]
//
// Contracts dim 3 from memory into C[16] registers, then delegates to ttv3.
// a[i] in AVX lanes (dim 0), b/c/d hold 4 weights each.
// coeffs points to the base of the 4x4x4x4 sub-tensor, with strides s1, s2, s3.
//
// SYMBOL "blazing_tensor_ttv4"
template<typename T1>
T1 casadi_blazing_tensor_ttv4(const T1* coeffs,
    casadi_int s1, casadi_int s2, casadi_int s3,
    simde__m256d a, simde__m256d b, simde__m256d c, simde__m256d d) {
  simde__m256d C[16];
  int j, k;
  // Broadcast dim 3 weights
  simde__m256d d0 = simde_mm256_permute4x64_pd(d, SIMDE_MM_SHUFFLE(0, 0, 0, 0));
  simde__m256d d1 = simde_mm256_permute4x64_pd(d, SIMDE_MM_SHUFFLE(1, 1, 1, 1));
  simde__m256d d2 = simde_mm256_permute4x64_pd(d, SIMDE_MM_SHUFFLE(2, 2, 2, 2));
  simde__m256d d3 = simde_mm256_permute4x64_pd(d, SIMDE_MM_SHUFFLE(3, 3, 3, 3));
  // Contract dim 3 from memory: C[j+4*k] = sum_l d[l] * coeffs[... + s3*l]
  for (j = 0; j < 4; ++j) {
    for (k = 0; k < 4; ++k) {
      const T1* p = coeffs + s1*j + s2*k;
      C[j+4*k] = simde_mm256_mul_pd(
          simde_mm256_loadu_pd(p), d0);
      C[j+4*k] = simde_mm256_fmadd_pd(
          simde_mm256_loadu_pd(p + s3), d1, C[j+4*k]);
      C[j+4*k] = simde_mm256_fmadd_pd(
          simde_mm256_loadu_pd(p + 2*s3), d2, C[j+4*k]);
      C[j+4*k] = simde_mm256_fmadd_pd(
          simde_mm256_loadu_pd(p + 3*s3), d3, C[j+4*k]);
    }
  }
  // Dims 0-2 handled by ttv3
  return casadi_blazing_tensor_ttv3<T1>(C, a, b, c);
}

// AVX2 tensor-times-vector for 5D cubic B-splines (degree 3, m=1).
//
//   result = sum_{i,j,k,l,m} a[i] * b[j] * c[k] * d[l] * e[m]
//            * coeffs[i + s1*j + s2*k + s3*l + s4*m]
//
// Contracts dims 3+4 fused from memory into C[16] registers,
// then delegates to ttv3. a[i] in AVX lanes (dim 0),
// b/c/d/e hold 4 weights each.
//
// SYMBOL "blazing_tensor_ttv5"
template<typename T1>
T1 casadi_blazing_tensor_ttv5(const T1* coeffs,
    casadi_int s1, casadi_int s2, casadi_int s3, casadi_int s4,
    simde__m256d a, simde__m256d b, simde__m256d c,
    simde__m256d d, simde__m256d e) {
  simde__m256d C[16];
  simde__m256d de[16];
  int j, k, l, m;
  // Pre-broadcast dim 3 and dim 4 weights (compile-time constant immediates)
  simde__m256d dbr[4] = {
    simde_mm256_permute4x64_pd(d, SIMDE_MM_SHUFFLE(0, 0, 0, 0)),
    simde_mm256_permute4x64_pd(d, SIMDE_MM_SHUFFLE(1, 1, 1, 1)),
    simde_mm256_permute4x64_pd(d, SIMDE_MM_SHUFFLE(2, 2, 2, 2)),
    simde_mm256_permute4x64_pd(d, SIMDE_MM_SHUFFLE(3, 3, 3, 3))
  };
  simde__m256d ebr[4] = {
    simde_mm256_permute4x64_pd(e, SIMDE_MM_SHUFFLE(0, 0, 0, 0)),
    simde_mm256_permute4x64_pd(e, SIMDE_MM_SHUFFLE(1, 1, 1, 1)),
    simde_mm256_permute4x64_pd(e, SIMDE_MM_SHUFFLE(2, 2, 2, 2)),
    simde_mm256_permute4x64_pd(e, SIMDE_MM_SHUFFLE(3, 3, 3, 3))
  };
  // Precompute outer product d x e (16 broadcast weights)
  for (l = 0; l < 4; ++l) {
    for (m = 0; m < 4; ++m) {
      de[l+4*m] = simde_mm256_mul_pd(dbr[l], ebr[m]);
    }
  }
  // Contract dims 3+4 from memory: C[j+4*k] = sum_{l,m} de[l+4*m] * coeffs[...]
  for (j = 0; j < 4; ++j) {
    for (k = 0; k < 4; ++k) {
      const T1* p = coeffs + s1*j + s2*k;
      C[j+4*k] = simde_mm256_setzero_pd();
      for (l = 0; l < 4; ++l) {
        for (m = 0; m < 4; ++m) {
          C[j+4*k] = simde_mm256_fmadd_pd(
              simde_mm256_loadu_pd(p + s3*l + s4*m),
              de[l+4*m], C[j+4*k]);
        }
      }
    }
  }
  // Dims 0-2 handled by ttv3
  return casadi_blazing_tensor_ttv3<T1>(C, a, b, c);
}

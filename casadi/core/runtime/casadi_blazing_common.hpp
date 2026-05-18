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

// C-REPLACE "casadi_blazing_low<T1>" "casadi_blazing_low"
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

// ===== Knot-span lookup =====

// SYMBOL "blazing_low"
// Forked from casadi_low for the blazing kernel. Same semantics
// (returns the largest i in [0, ng-2] such that grid[i] <= x < grid[i+1],
// with clamping at both ends), but three faster code paths:
//   lookup_mode = 0 (linear): AVX2-vectorised scan, 4 grid points per
//                             SIMD iteration; ~4x fewer branches than scalar.
//   lookup_mode = 1 (exact):  if grid_inv is non-null, treats grid_inv as a
//                             2-element {intercept, slope} pair and computes
//                             i = floor(x*slope + intercept) with a single
//                             FMA --- no FP divide, no extra subtract.
//                             Falls back to the original divide when null.
//   lookup_mode = 2 (binary): branchless bisection via cmov (data-dependent
//                             updates only, no mispredictable branches inside
//                             the loop).
template<typename T1>
casadi_int casadi_blazing_low(T1 x, const T1* grid, casadi_int ng,
                              casadi_int lookup_mode, const T1* grid_inv) {
  switch (lookup_mode) {
    case 1:
      {
        // 'exact' --- uniform-grid direct lookup
        casadi_int ret;
        if (grid_inv) {
          // Pre-baked at codegen time:
          //   slope     = (ng-1) / (grid[ng-1] - grid[0])
          //   intercept = -grid[0] * slope
          // so x*slope + intercept lands directly in the lookup index.
          // One FMA, one truncate.
          T1 intercept = grid_inv[0];
          T1 slope     = grid_inv[1];
          ret = (casadi_int) (x * slope + intercept);
        } else {
          // Fallback: original formula with the FP divide.
          T1 g0 = grid[0];
          T1 dg = grid[ng-1] - g0;
          ret = (casadi_int) ((x - g0) * (ng-1) / dg);
        }
        if (ret < 0) ret = 0;
        if (ret > ng-2) ret = ng-2;
        return ret;
      }
    case 2:
      {
        // 'binary' --- Skarupke-style branchless lower_bound
        // (https://mhdm.dev/posts/sb_lower_bound/), then converted to
        // casadi_low semantics (largest i with grid[i] <= x, clamped to
        // [0, ng-2]).
        //
        // The `if (cmp) lo += len + 1` form lets the compiler fuse the
        // compare-cmovae-add into a tight chain; measurably faster than
        // the `lo = cmp ? probe : lo; len -= half;` style. Benchmarked
        // at n=100 against std::lower_bound, Skarupke, sb, sbm, bb, sbp;
        // Skarupke is the consistent winner across n=100..10000 in both
        // fixed-x and random-x regimes.
        casadi_int lo = 0;
        casadi_int len = ng;
        while (len > 0) {
          len /= 2;
          if (grid[lo + len] < x) lo += len + 1;
        }
        // Convert lower_bound -> "largest i s.t. grid[i] <= x".
        if (lo > 0) lo--;
        if (lo > ng - 2) lo = ng - 2;
        return lo;
      }
    default:
      {
        // 'linear' --- AVX2-vectorised forward scan.
        // Compare 4 grid points to x per iteration; tzcnt picks the first
        // lane with x < grid[i+1]. Falls back to scalar for the tail.
        casadi_int n = ng - 2;
        if (n <= 0) return 0;
        simde__m256d xv = simde_mm256_set1_pd((double) x);
        casadi_int i = 0;
        for (; i + 4 <= n; i += 4) {
          simde__m256d gv  = simde_mm256_loadu_pd((const double*)(grid + i + 1));
          simde__m256d cmp = simde_mm256_cmp_pd(xv, gv, SIMDE_CMP_LT_OQ);
          int m = simde_mm256_movemask_pd(cmp);
          if (m) return i + __builtin_ctz((unsigned) m);
        }
        for (; i < n; ++i) {
          if (x < grid[i+1]) return i;
        }
        return n;
      }
  }
}

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
// dim_cache points at THIS dimension's slice of the global cache (or NULL).
// The caller is responsible for advancing the pointer between dims.
//
// Per-dimension cache layout (always emitted unconditionally when the cache
// is present, even when lookup_mode != 1):
//     [intercept, slope, inv1[n_k], inv2[n_k], inv3[n_k]]
// where slope     = (ng-1) / (grid[ng-1] - grid[0])
//       intercept = -grid[0] * slope
// and ng = n_k - 2*degree. The first two entries feed 'exact' lookup as a
// single FMA; the inv1/inv2/inv3 spans feed the de Boor recurrence.
//
// If inv2_out/inv3_out are non-NULL, they receive pre-positioned derivative
// base pointers that can be passed directly to dbasis/d2basis.
template<typename T1>
casadi_int casadi_blazing_boor_init(
    T1 x, const T1* all_knots, const T1* dim_cache,
    casadi_int knot_offset, casadi_int knot_offset_next,
    casadi_int lookup_mode,
    simde__m256d* d0, simde__m256d* d1, simde__m256d* d2,
    const T1** inv2_out, const T1** inv3_out) {
    casadi_int degree = 3;
    const T1* knots = all_knots + knot_offset;
    casadi_int n_knots = knot_offset_next - knot_offset;
    casadi_int n_b = n_knots - degree - 1;
    // First two entries of dim_cache are {intercept, slope} for 'exact' lookup.
    const T1* grid_inv = dim_cache;
    casadi_int L = casadi_blazing_low<T1>(x, knots + degree,
                                          n_knots - 2*degree,
                                          lookup_mode, grid_inv);
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
    if (dim_cache) {
      // Skip the {intercept, slope} prefix to reach the inv1/inv2/inv3 spans.
      const T1* inv_base = dim_cache + 2;
      inv1 = inv_base + start;
      inv2 = inv_base + n_knots + start;
      inv3 = inv_base + 2*n_knots + start;
      if (inv2_out) *inv2_out = inv_base + n_knots + start + 1;
      if (inv3_out) *inv3_out = inv_base + 2*n_knots + start + 1;
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
  int j, k, l;
  // Broadcast dim 3 weights into an array indexable by l.
  simde__m256d dbr[4] = {
    simde_mm256_permute4x64_pd(d, SIMDE_MM_SHUFFLE(0, 0, 0, 0)),
    simde_mm256_permute4x64_pd(d, SIMDE_MM_SHUFFLE(1, 1, 1, 1)),
    simde_mm256_permute4x64_pd(d, SIMDE_MM_SHUFFLE(2, 2, 2, 2)),
    simde_mm256_permute4x64_pd(d, SIMDE_MM_SHUFFLE(3, 3, 3, 3))
  };
  // Contract dim 3 from memory: C[j+4*k] = sum_l d[l] * coeffs[... + s3*l]
  // Loop order: large stride s3 outer, so the 4x4 (j,k) gather stays in
  // L1d; reversed order conflict-thrashes when s3 is a power-of-2 stride.
  for (int idx = 0; idx < 16; ++idx) C[idx] = simde_mm256_setzero_pd();
  for (l = 0; l < 4; ++l) {
    const T1* base = coeffs + s3*l;
    simde__m256d dl = dbr[l];
    for (j = 0; j < 4; ++j) {
      for (k = 0; k < 4; ++k) {
        C[j+4*k] = simde_mm256_fmadd_pd(
            simde_mm256_loadu_pd(base + s1*j + s2*k), dl, C[j+4*k]);
      }
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
  // Loop order: large strides s3, s4 outer, so the 4x4 (j,k) gather stays
  // in L1d; reversed order conflict-thrashes on power-of-2 strides.
  for (int idx = 0; idx < 16; ++idx) C[idx] = simde_mm256_setzero_pd();
  for (l = 0; l < 4; ++l) {
    for (m = 0; m < 4; ++m) {
      const T1* base = coeffs + s3*l + s4*m;
      simde__m256d delm = de[l+4*m];
      for (j = 0; j < 4; ++j) {
        for (k = 0; k < 4; ++k) {
          C[j+4*k] = simde_mm256_fmadd_pd(
              simde_mm256_loadu_pd(base + s1*j + s2*k),
              delm, C[j+4*k]);
        }
      }
    }
  }
  // Dims 0-2 handled by ttv3
  return casadi_blazing_tensor_ttv3<T1>(C, a, b, c);
}

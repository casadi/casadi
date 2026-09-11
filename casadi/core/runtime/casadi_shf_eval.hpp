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

// SYMBOL "de_casteljau"
// Evaluates the Bernstein polynomial with control points b[0..n] at t, in place.
template<typename T1>
T1 casadi_de_casteljau(casadi_int n, T1* b, T1 t) {
  casadi_int i, j;
  T1 u;
  u = 1-t;
  for (j=1; j<=n; ++j)
    for (i=0; i<=n-j; ++i) b[i] = b[i]*u + b[i+1]*t;
  return b[0];
}

// SYMBOL "shf_bernstein"
// p-th derivative of s_k(t), the Bernstein polynomial of degree n=2k+1 of
// R(t)=max(0,2t-1): p-fold forward differencing of the control points, then de
// Casteljau at degree n-p, scaled by n!/(n-p)!. Needs beta[n+1] scratch.
template<typename T1>
T1 casadi_shf_bernstein(casadi_int k, casadi_int p, T1 t, T1* beta) {
  casadi_int n, i, d;
  T1 f;
  n = 2*k+1;
  if (p > n) return 0;
  for (i=0; i<=n; ++i) beta[i] = (i<=k) ? 0.0 : (2.0*i-n)/n;
  for (d=0; d<p; ++d)
    for (i=0; i<n-d; ++i) beta[i] = beta[i+1]-beta[i];
  f = 1;
  for (d=0; d<p; ++d) f *= (T1)(n-d);
  return f*casadi_de_casteljau(n-p, beta, t);
}

// SYMBOL "shf_axis"
// Order-p weights of one axis: writes 'width' values into w and the stencil start
// into *start. Three regimes: bulk (two nonzero weights, padded with a zero on
// whichever side keeps the stencil in range) and epsilon-ball around an interior
// grid point (three weights). Outside the grid the bulk branch extrapolates
// linearly, because t is never clamped -- only the interval index is. In the bulk
// the map is affine, so every derivative of order two and up vanishes identically.
template<typename T1>
void casadi_shf_axis(casadi_int k, const T1* g, casadi_int ng, const T1* inv_h,
    T1 x, T1 epsilon, casadi_int p, casadi_int width, casadi_int lookup_mode,
    casadi_int* start, T1* w, T1* beta) {
  casadi_int j, c, i;
  T1 e, t, ta, sa, sb, A, B, sc, hm, hp;
  j = casadi_low(x, g, ng, lookup_mode);
  // Locate a smoothing centre, if any. epsilon is an absolute width. The j>0 and
  // j+1<ng-1 guards keep the ends unsmoothed: with linear extrapolation the first
  // and last grid points carry no kink, so there is nothing there to round off.
  c = -1;
  e = epsilon;
  if (width > 2 && e > 0 && j > 0 && x - g[j] <= e) c = j;
  if (width > 2 && c < 0 && e > 0 && j+1 < ng-1 && g[j+1] - x <= e) c = j+1;
  for (i=0; i<width; ++i) w[i] = 0;
  if (c < 0) {
    if (width > 2 && j > 0) {
      *start = j-1;
      if (p == 0) { t = (x - g[j])*inv_h[j]; w[1] = 1-t; w[2] = t; }
      else if (p == 1) { w[1] = -inv_h[j]; w[2] = inv_h[j]; }
    } else {
      *start = j;
      if (p == 0) { t = (x - g[j])*inv_h[j]; w[0] = 1-t; w[1] = t; }
      else if (p == 1) { w[0] = -inv_h[j]; w[1] = inv_h[j]; }
    }
  } else {
    *start = c-1;
    hm = g[c]-g[c-1];
    hp = g[c+1]-g[c];
    ta = 0.5 + (g[c]-x)/(2*e);
    // One polynomial evaluation supplies both mirror weights: differentiating
    // s_k(t) - s_k(1-t) = 2t-1 gives s^(p)(1-t) = (-1)^p (s^(p)(t) - l_p) with
    // l_0 = 2t-1, l_1 = 2 and l_p = 0 beyond.
    sa = casadi_shf_bernstein(k, p, ta, beta);
    sb = sa - (p == 0 ? 2*ta-1 : (p == 1 ? 2 : 0));
    if (p % 2) sb = -sb;
    if (p == 0) {
      A = e*sa/hm; B = e*sb/hp;
      w[0] = A; w[1] = 1-A-B; w[2] = B;
    } else {
      // chain rule: d(ta)/dx = -1/(2e), d(tb)/dx = +1/(2e)
      sc = 1;
      for (i=0; i<p; ++i) sc /= (2*e);
      A = e*sa/hm; B = e*sb/hp;
      if (p % 2) A = -A;
      A *= sc; B *= sc;
      w[0] = A; w[1] = -A-B; w[2] = B;
    }
  }
}

// SYMBOL "shf_ttv_multi"
// Tensor-times-vector carrying one accumulator per requested derivative
// multi-index instead of a scalar weight. A single traversal of the coefficient
// stencil therefore yields every mixed partial in the set -- the coefficients are
// never transformed and the cost does not grow with the derivative order.
//   all_w  : per axis, weight vectors of order 0..P, stride 3
//   multi  : nb-by-ndim table of multi-indices
// wofs[dim*nb+b] is the precomputed offset into all_w of the weight vector that
// accumulator b needs on axis dim. Resolving it once per evaluation instead of
// per stencil visit keeps the innermost loop free of index arithmetic.
template<typename T1>
void casadi_shf_ttv_multi(T1* ret, casadi_int dim, casadi_int ndim, const T1* all_w,
    const casadi_int* width, const casadi_int* starts, const casadi_int* strides,
    const casadi_int* wofs, casadi_int nb,
    const T1* c, casadi_int m, casadi_int offset, T1* W) {
  casadi_int i, j, b, n_w;
  const casadi_int* ofs;
  const T1 *Win;
  T1 *Wout;
  n_w = width[dim];
  ofs = wofs + dim*nb;
  Win = W + (dim+1)*nb;
  Wout = W + dim*nb;
  if (dim == 0) {
    // At the leaf the accumulator is dead after use, so keep it in a register
    // instead of round-tripping through W. Fusing the two loops also keeps gcc
    // from vectorising a pair of very short runtime-bounded loops, which costs
    // more in setup than it recovers.
    const T1* coeff = c + offset + starts[0]*m;
    for (i=0; i<n_w; ++i) {
      for (b=0; b<nb; ++b) {
        T1 wb = Win[b]*all_w[ofs[b]+i];
        for (j=0; j<m; ++j) ret[b*m+j] += wb*coeff[i*m+j];
      }
    }
  } else {
    for (i=0; i<n_w; ++i) {
      for (b=0; b<nb; ++b) Wout[b] = Win[b]*all_w[ofs[b]+i];
      casadi_shf_ttv_multi(ret, dim-1, ndim, all_w, width, starts, strides,
        wofs, nb, c, m, offset + (starts[dim]+i)*strides[dim], W);
    }
  }
}

// SYMBOL "shf_eval_multi"
// ret is m-by-nb: column b holds the mixed partial d^multi[b] f. One call covers
// value (nb=1, multi=0), jacobian, hessian and any higher order alike.
// epsilon holds one half-width per axis (eps_stride=1) or a single one shared by
// all axes (eps_stride=0); a null pointer means no smoothing at all.
template<typename T1>
void casadi_shf_eval_multi(T1* ret, casadi_int ndim, const T1* grid,
    const casadi_int* offset, const T1* inv_h, const casadi_int* width,
    const casadi_int* strides, const T1* c, casadi_int m, const T1* x,
    const T1* epsilon, casadi_int eps_stride,
    casadi_int k, const casadi_int* multi, casadi_int nb, casadi_int P,
    const casadi_int* lookup_mode, casadi_int* iw, T1* w) {
  casadi_int r, q, b, ng;
  casadi_int *starts, *wofs;
  T1 *all_w, *beta, *W;
  starts = iw; iw += ndim;
  wofs = iw; iw += ndim*nb;
  all_w = w; w += ndim*(P+1)*3;
  W = w; w += (ndim+1)*nb;
  beta = w;
  for (r=0; r<ndim; ++r) {
    ng = offset[r+1]-offset[r];
    for (q=0; q<=P; ++q) {
      casadi_shf_axis(k, grid+offset[r], ng, inv_h+offset[r]-r, x[r],
        epsilon ? epsilon[r*eps_stride] : 0, q,
        width[r], lookup_mode[r], starts+r, all_w + (r*(P+1)+q)*3, beta);
    }
  }
  for (r=0; r<ndim; ++r)
    for (b=0; b<nb; ++b) wofs[r*nb+b] = (r*(P+1) + multi[b*ndim+r])*3;
  if (nb == 1) {
    // Single output (the value, or any one mixed partial): the accumulator is a
    // lone scalar, so hand the contraction to casadi_tensor_ttv, which carries it
    // in a register rather than a memory cell. Compact the selected weight vectors
    // into W first, since tensor_ttv wants them contiguous per axis.
    casadi_int off = 0;
    for (r=0; r<ndim; ++r) {
      wofs[ndim+r] = off;
      for (q=0; q<width[r]; ++q) W[off+q] = all_w[wofs[r]+q];
      off += width[r];
    }
    wofs[ndim+ndim] = off;
    casadi_tensor_ttv(ret, ndim-1, ndim, W, wofs+ndim, starts, strides, c, m, 1.0, 0);
    return;
  }
  for (b=0; b<nb; ++b) W[ndim*nb+b] = 1;
  casadi_shf_ttv_multi(ret, ndim-1, ndim, all_w, width, starts, strides,
    wofs, nb, c, m, 0, W);
}

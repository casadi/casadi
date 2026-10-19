// Derivation of the coefficients of the Rosenbrock-Wanner method used by the
// 'rosenbrock' integrator plugin. Standalone: plain C++ and libquadmath.
//
//   g++ -O2 -std=c++17 -pthread rosenbrock_coefficients.cpp -lquadmath -o rosenbrock_coefficients
//   ./rosenbrock_coefficients [n_starts | -start] [seed] [conditions] [optimization iterations]
// writes rosenbrock_tableau.hpp. The tableau in casadi/solvers was obtained with
//   ./rosenbrock_coefficients -1980 0 7 300
// i.e. start 1980 of a 2560-start search, chosen among the A-stable solutions (main and
// embedded method) for its efficiency on a set of stiff ODE, DAE and Prothero-Robinson tests.
//
// The method: s = 8 stages, stiffly accurate, order 5 for ODEs and index-1 DAEs, with a
// stiffly accurate embedded method of order 4 (s = 7), in the original form
//   (M - h gamma J) k_i = h F(t0 + alpha_i h, y0 + sum_j alpha_ij k_j)
//                         + h J sum_j gamma_ij k_j + h^2 gamma_i F_t,
//   y1 = y0 + sum_i b_i k_i.
// Notation: B = (beta_ij), beta_ij = alpha_ij + gamma_ij, beta_ii = gamma, W = B^{-1}.
// The structure follows the construction principles of G. Steinebach, "Construction of
// Rosenbrock-Wanner method Rodas5P and numerical benchmarks within the Julia Differential
// Equations package", BIT 63:27 (2023):
//   b_i = beta_8i, b_8 = gamma; bhat_i = beta_7i, bhat_7 = gamma;
//   alpha_8i = beta_7i, alpha_7i = beta_6i, alpha_6 = alpha_7 = alpha_8 = 1;
//   alpha_21 = 3 gamma, beta_21 = 0.
// The order conditions are the ones of Tables 2-4 of that paper: ODE conditions are
// generated from rooted trees, the DAE (Roche), index-2, W-method and Prothero-Robinson
// conditions are spelled out below. The free coefficients are found by a multistart
// Levenberg-Marquardt search in double precision. The solutions form a manifold, on which
// the principal error coefficients are minimized, keeping the main and the embedded method
// A-stable. A final Gauss-Newton refinement is done in quadruple precision.

#include <quadmath.h>
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <functional>
#include <mutex>
#include <random>
#include <thread>
#include <vector>

typedef __float128 quad;

// ---------------------------------------------------------------- dual numbers
// Value and gradient with respect to the NV free coefficients
constexpr int NV = 42;
template<typename T>
struct Dual {
  T v;
  T d[NV];
  Dual() : v(0) { std::fill(d, d + NV, T(0)); }
  Dual(double c) : v(c) { std::fill(d, d + NV, T(0)); }
  static Dual var(T value, int i) { Dual r; r.v = value; r.d[i] = T(1); return r; }
};
template<typename T> Dual<T> operator+(const Dual<T>& a, const Dual<T>& b) {
  Dual<T> r; r.v = a.v + b.v; for (int i = 0; i < NV; ++i) r.d[i] = a.d[i] + b.d[i]; return r;
}
template<typename T> Dual<T> operator-(const Dual<T>& a, const Dual<T>& b) {
  Dual<T> r; r.v = a.v - b.v; for (int i = 0; i < NV; ++i) r.d[i] = a.d[i] - b.d[i]; return r;
}
template<typename T> Dual<T> operator-(const Dual<T>& a) { return Dual<T>(0.) - a; }
template<typename T> Dual<T> operator*(const Dual<T>& a, const Dual<T>& b) {
  Dual<T> r; r.v = a.v * b.v;
  for (int i = 0; i < NV; ++i) r.d[i] = a.d[i] * b.v + a.v * b.d[i];
  return r;
}
template<typename T> Dual<T> operator/(const Dual<T>& a, const Dual<T>& b) {
  Dual<T> r; r.v = a.v / b.v;
  for (int i = 0; i < NV; ++i) r.d[i] = (a.d[i] - r.v * b.d[i]) / b.v;
  return r;
}
template<typename T> Dual<T> operator*(double c, const Dual<T>& a) { return Dual<T>(c) * a; }
template<typename T> Dual<T>& operator+=(Dual<T>& a, const Dual<T>& b) { return a = a + b; }

template<typename R> using Vec = std::vector<R>;

// ------------------------------------------------------------- method structure
constexpr int S = 8;

// Free coefficients: gamma, beta_ij (2 < i <= 8, j < i, except beta_21), alpha_ij (i = 3..6)
//   gamma: 1, beta: 27, alpha: 2 + 3 + 4 + 5 = 14  -> NV = 42
template<typename R>
struct Method {
  R gamma;
  R A[S][S], B[S][S];
  Vec<R> b, bh;
};

template<typename R>
Method<R> assemble(const R* x) {
  Method<R> m;
  int k = 0;
  m.gamma = x[k++];
  for (int i = 0; i < S; ++i) {
    for (int j = 0; j < S; ++j) { m.A[i][j] = R(0.); m.B[i][j] = R(0.); }
    m.B[i][i] = m.gamma;
  }
  m.A[1][0] = 3. * m.gamma;  // alpha_21 = 3 gamma, beta_21 = 0
  for (int i = 2; i < S; ++i) for (int j = 0; j < i; ++j) m.B[i][j] = x[k++];
  for (int i = 2; i < 6; ++i) for (int j = 0; j < i; ++j) m.A[i][j] = x[k++];
  for (int j = 0; j < 6; ++j) m.A[6][j] = m.B[5][j];   // alpha_7i = beta_6i
  for (int j = 0; j < 7; ++j) m.A[7][j] = m.B[6][j];   // alpha_8i = beta_7i
  m.b.resize(S);
  m.bh.resize(S);
  for (int j = 0; j < S; ++j) { m.b[j] = m.B[7][j]; m.bh[j] = j < 7 ? m.B[6][j] : R(0.); }
  return m;
}

// ----------------------------------------------------------------- vector algebra
template<typename R> Vec<R> ones(int n) { return Vec<R>(n, R(1.)); }
template<typename R> Vec<R> mul(const R M[S][S], const Vec<R>& v, int n) {
  Vec<R> r(n, R(0.));
  for (int i = 0; i < n; ++i) for (int j = 0; j <= i; ++j) r[i] += M[i][j] * v[j];
  return r;
}
template<typename R> Vec<R> times(const Vec<R>& a, const Vec<R>& b) {
  Vec<R> r(a.size()); for (size_t i = 0; i < a.size(); ++i) r[i] = a[i] * b[i]; return r;
}
template<typename R> R dot(const Vec<R>& a, const Vec<R>& b) {
  R r(0.); for (size_t i = 0; i < a.size(); ++i) r += a[i] * b[i]; return r;
}
template<typename R> Vec<R> rowsum(const R M[S][S], int n) { return mul(M, ones<R>(n), n); }
template<typename R> Vec<R> pw(const Vec<R>& a, int e) {
  Vec<R> r = ones<R>(a.size()); for (int k = 0; k < e; ++k) r = times(r, a); return r;
}
// W v with W = B^{-1}: forward substitution
template<typename R> Vec<R> wsolve(const R B[S][S], const Vec<R>& v, int n) {
  Vec<R> r(n);
  for (int i = 0; i < n; ++i) {
    R acc = v[i];
    for (int j = 0; j < i; ++j) acc = acc - B[i][j] * r[j];
    r[i] = acc / B[i][i];
  }
  return r;
}

// ------------------------------------------------------------------ rooted trees
// A tree is the list of its subtrees
struct Tree { std::vector<Tree> c; };

bool tree_less(const Tree& a, const Tree& b);
int tree_nodes(const Tree& t) { int n = 1; for (auto& c : t.c) n += tree_nodes(c); return n; }
bool tree_less(const Tree& a, const Tree& b) {
  if (tree_nodes(a) != tree_nodes(b)) return tree_nodes(a) < tree_nodes(b);
  if (a.c.size() != b.c.size()) return a.c.size() < b.c.size();
  for (size_t i = 0; i < a.c.size(); ++i) {
    if (tree_less(a.c[i], b.c[i])) return true;
    if (tree_less(b.c[i], a.c[i])) return false;
  }
  return false;
}

// All rooted trees with n nodes (children sorted, so each tree appears once)
std::vector<Tree> trees(int n) {
  std::vector<Tree> out;
  if (n == 1) { out.push_back(Tree()); return out; }
  // Distribute n-1 nodes over children with non-increasing (size, tree) order
  std::function<void(int, Tree&, const Tree*)> rec = [&](int left, Tree& cur, const Tree* maxc) {
    if (left == 0) { out.push_back(cur); return; }
    for (int sz = left; sz >= 1; --sz) {
      for (auto& c : trees(sz)) {
        if (maxc && tree_less(*maxc, c)) continue;
        cur.c.push_back(c);
        rec(left - sz, cur, &cur.c.back());
        cur.c.pop_back();
      }
    }
  };
  Tree t;
  rec(n - 1, t, nullptr);
  return out;
}

// Symmetry factor: product over groups of m identical subtrees of m! sigma(subtree)^m
int tree_symmetry(const Tree& t) {
  int r = 1;
  for (size_t i = 0; i < t.c.size();) {
    size_t j = i;
    while (j < t.c.size() && !tree_less(t.c[i], t.c[j]) && !tree_less(t.c[j], t.c[i])) ++j;
    int m = j - i, f = 1;
    for (int k = 2; k <= m; ++k) f *= k;
    r *= f;
    for (int k = 0; k < m; ++k) r *= tree_symmetry(t.c[i]);
    i = j;
  }
  return r;
}

int tree_density(const Tree& t) {
  int r = tree_nodes(t);
  for (auto& c : t.c) r *= tree_density(c);
  return r;
}

// Elementary weights: a node with one child is linked by B, otherwise by A
template<typename R> Vec<R> phi(const Tree& t, const Method<R>& m, int n) {
  if (t.c.empty()) return ones<R>(n);
  if (t.c.size() == 1) return mul(m.B, phi(t.c[0], m, n), n);
  Vec<R> r = ones<R>(n);
  for (auto& c : t.c) r = times(r, mul(m.A, phi(c, m, n), n));
  return r;
}

// ----------------------------------------------------------------- the conditions
// Exact rational right-hand sides
struct Ratio { int num, den; };
template<typename R> R ratio(const Ratio& r) { return R(double(r.num)) / R(double(r.den)); }

template<typename R>
void ode_conditions(const Method<R>& m, const Vec<R>& b, int n, int pmax, Vec<R>& res) {
  for (int p = 1; p <= pmax; ++p) {
    for (auto& t : trees(p)) res.push_back(dot(b, phi(t, m, n)) - R(1.) / R(double(tree_density(t))));
  }
}

// Index-1 DAE (No 18-40), index-2 (41-44) and W-method (45-47) condition vectors
template<typename R>
void dae_vectors(const Method<R>& m, int n, std::vector<std::pair<Vec<R>, Ratio>>& V) {
  Vec<R> al = rowsum(m.A, n), be = rowsum(m.B, n);
  auto A = [&](const Vec<R>& v) { return mul(m.A, v, n); };
  auto Bm = [&](const Vec<R>& v) { return mul(m.B, v, n); };
  auto W = [&](const Vec<R>& v) { return wsolve(m.B, v, n); };
  Vec<R> Ab = A(be), AWa2 = A(W(pw(al, 2))), AWWa2 = A(W(W(pw(al, 2))));
  V.resize(48);
  V[18] = {W(pw(al, 2)), {1, 1}};
  V[19] = {times(al, A(W(pw(al, 2)))), {1, 4}};
  V[20] = {W(pw(al, 3)), {1, 1}};
  V[21] = {W(times(al, Ab)), {1, 2}};
  V[22] = {W(times(al, AWa2)), {1, 1}};
  V[23] = {times(al, A(W(pw(al, 3)))), {1, 5}};
  V[24] = {times(al, A(W(times(al, Ab)))), {1, 10}};
  V[25] = {times(al, A(W(times(al, AWa2)))), {1, 5}};
  V[26] = {times(Ab, AWa2), {1, 10}};
  V[27] = {times(pw(al, 2), AWa2), {1, 5}};
  V[28] = {Bm(times(al, AWa2)), {1, 20}};
  V[29] = {times(AWa2, AWa2), {1, 5}};
  V[30] = {W(pw(al, 4)), {1, 1}};
  V[31] = {W(times(pw(al, 2), Ab)), {1, 2}};
  V[32] = {W(times(pw(al, 2), AWa2)), {1, 1}};
  V[33] = {W(times(al, A(Bm(be)))), {1, 6}};
  V[34] = {W(times(al, A(pw(al, 2)))), {1, 3}};
  V[35] = {W(times(al, A(W(pw(al, 3))))), {1, 1}};
  V[36] = {W(times(al, A(W(times(al, Ab))))), {1, 2}};
  V[37] = {W(times(al, A(W(times(al, AWa2))))), {1, 1}};
  V[38] = {W(times(Ab, Ab)), {1, 4}};
  V[39] = {W(times(Ab, AWa2)), {1, 2}};
  V[40] = {W(times(AWa2, AWa2)), {1, 1}};
  V[41] = {W(W(pw(al, 2))), {2, 1}};
  V[42] = {times(al, AWWa2), {2, 3}};
  V[43] = {W(times(al, AWWa2)), {2, 1}};
  V[44] = {times(AWWa2, AWWa2), {4, 3}};
  Vec<R> bmal(n); for (int i = 0; i < n; ++i) bmal[i] = be[i] - al[i];
  V[45] = {bmal, {0, 1}};
  V[46] = {A(W(al)), {1, 2}};
  V[47] = {W(al), {1, 1}};
}

// Prothero-Robinson conditions No 48 (A_0..A_s) and No 49 (B_0..B_{s-1})
template<typename R>
void prothero_robinson(const Method<R>& m, const Vec<R>& b, int n, Vec<R>& Ac, Vec<R>& Bc) {
  // Strictly lower part of B and its row sums
  R Bs[S][S];
  for (int i = 0; i < S; ++i) for (int j = 0; j < S; ++j) Bs[i][j] = j < i ? m.B[i][j] : R(0.);
  Vec<R> bs = rowsum(Bs, n), al = rowsum(m.A, n);
  auto chain = [&](int nu, Vec<R> v) { for (int k = 0; k < nu; ++k) v = mul(Bs, v, n); return v; };
  auto M = [&](int nu) { return nu < 0 ? dot(b, ones<R>(n)) : dot(b, chain(nu, bs)); };
  auto N = [&](int sig, int nu) {
    R f(sig == 2 ? 2. : 6.);
    return (nu < 0 ? dot(b, ones<R>(n)) : dot(b, chain(nu, pw(al, sig)))) / f;
  };
  R g = m.gamma;
  Ac.clear(); Bc.clear();
  Ac.push_back(-N(2, -1) + g * M(-1) + M(0));
  for (int i = 1; i < n; ++i) Ac.push_back(-N(2, i - 1) + 2. * g * M(i - 1) + g * g * M(i - 2) + M(i));
  Ac.push_back(g * g * M(n - 2));
  Bc.push_back(-N(3, -1) + N(2, 0));
  for (int i = 1; i < n - 1; ++i) Bc.push_back(-N(3, i - 1) + g * N(2, i - 1) + N(2, i));
  Bc.push_back(-N(3, n - 2) + g * N(2, n - 2));
}

// Condition sets
enum { WITH_INDEX2 = 1, WITH_W = 2, WITH_PR = 4 };

template<typename R>
Vec<R> residuals(const R* x, int flags) {
  Method<R> m = assemble(x);
  Vec<R> res;
  std::vector<std::pair<Vec<R>, Ratio>> V, Vh;
  dae_vectors(m, S, V);
  dae_vectors(m, S, Vh);
  // Main method: ODE order 5 and index-1 DAE order 5
  ode_conditions(m, m.b, S, 5, res);
  for (int no = 18; no <= 40; ++no) res.push_back(dot(m.b, V[no].first) - ratio<R>(V[no].second));
  // Embedded method: order 4
  ode_conditions(m, m.bh, S, 4, res);
  for (int no = 18; no <= 22; ++no) res.push_back(dot(m.bh, V[no].first) - ratio<R>(V[no].second));
  if (flags & WITH_INDEX2) {
    for (int no = 41; no <= 44; ++no) res.push_back(dot(m.b, V[no].first) - ratio<R>(V[no].second));
    res.push_back(dot(m.bh, V[41].first) - ratio<R>(V[41].second));
  }
  if (flags & WITH_W) {
    for (int no = 45; no <= 47; ++no) res.push_back(dot(m.b, V[no].first) - ratio<R>(V[no].second));
  }
  if (flags & WITH_PR) {
    Vec<R> Ac, Bc;
    prothero_robinson(m, m.b, S, Ac, Bc);
    for (auto& e : Ac) res.push_back(e);
    for (auto& e : Bc) res.push_back(e);
    Vec<R> bh7(m.bh.begin(), m.bh.begin() + 7);
    prothero_robinson(m, bh7, 7, Ac, Bc);
    for (auto& e : Ac) res.push_back(e);
  }
  // Structure: alpha_6 = 1, alpha_7 = beta_6 = 1
  R a6(0.), b6(0.);
  for (int j = 0; j < S; ++j) { a6 += m.A[5][j]; b6 += m.B[5][j]; }
  res.push_back(a6 - R(1.));
  res.push_back(b6 - R(1.));
  return res;
}

// Order-6 ODE residuals of the main method (error constants)
template<typename R>
Vec<R> error_terms(const R* x) {
  Method<R> m = assemble(x);
  Vec<R> res;
  for (auto& t : trees(6)) res.push_back(dot(m.b, phi(t, m, S)) - R(1.) / R(double(tree_density(t))));
  return res;
}

// ------------------------------------------------------------ linear least squares
// Minimum-norm least-squares solution of min ||J d + r||^2 + mu ||d||^2 via normal equations
// in Cholesky form (small, dense)
template<typename T>
bool lm_step(const std::vector<std::vector<T>>& J, const Vec<T>& r, T mu, Vec<T>& d) {
  int m = J.size(), n = NV;
  std::vector<std::vector<T>> H(n, Vec<T>(n, T(0)));
  Vec<T> g(n, T(0));
  for (int k = 0; k < m; ++k) {
    for (int i = 0; i < n; ++i) {
      g[i] += J[k][i] * r[k];
      for (int j = 0; j <= i; ++j) H[i][j] += J[k][i] * J[k][j];
    }
  }
  for (int i = 0; i < n; ++i) H[i][i] += mu;
  // Cholesky H = L L^T
  for (int j = 0; j < n; ++j) {
    T s = H[j][j];
    for (int k = 0; k < j; ++k) s -= H[j][k] * H[j][k];
    if (!(s > T(0))) return false;
    T Ljj = sqrtq(s);
    H[j][j] = Ljj;
    for (int i = j + 1; i < n; ++i) {
      T t = H[i][j];
      for (int k = 0; k < j; ++k) t -= H[i][k] * H[j][k];
      H[i][j] = t / Ljj;
    }
  }
  d.assign(n, T(0));
  for (int i = 0; i < n; ++i) {
    T t = -g[i];
    for (int k = 0; k < i; ++k) t -= H[i][k] * d[k];
    d[i] = t / H[i][i];
  }
  for (int i = n - 1; i >= 0; --i) {
    T t = d[i];
    for (int k = i + 1; k < n; ++k) t -= H[k][i] * d[k];
    d[i] = t / H[i][i];
  }
  return true;
}

template<typename T> T sq_norm(const Vec<T>& v) { T s(0); for (auto& e : v) s += e * e; return s; }

// Residuals and their Jacobian
template<typename T>
void eval(const Vec<T>& x, int flags, Vec<T>& r, std::vector<std::vector<T>>& J) {
  std::vector<Dual<T>> xd(NV);
  for (int i = 0; i < NV; ++i) xd[i] = Dual<T>::var(x[i], i);
  Vec<Dual<T>> rd = residuals(xd.data(), flags);
  r.resize(rd.size());
  J.assign(rd.size(), Vec<T>(NV));
  for (size_t k = 0; k < rd.size(); ++k) {
    r[k] = rd[k].v;
    for (int i = 0; i < NV; ++i) J[k][i] = rd[k].d[i];
  }
}

// Levenberg-Marquardt
template<typename T>
T levenberg_marquardt(Vec<T>& x, int flags, int max_iter, T tol) {
  Vec<T> r, rn, d, xn(NV);
  std::vector<std::vector<T>> J, Jn;
  eval(x, flags, r, J);
  T f = sq_norm(r), mu(1e-3);
  for (int it = 0; it < max_iter && f > tol * tol; ++it) {
    if (!lm_step(J, r, mu, d)) { mu *= T(10); continue; }
    for (int i = 0; i < NV; ++i) xn[i] = x[i] + d[i];
    eval(xn, flags, rn, Jn);
    T fn = sq_norm(rn);
    if (fn < f) {
      x = xn; r = rn; J = Jn; f = fn;
      mu = std::max(mu / T(3), T(1e-30));
    } else {
      mu *= T(4);
      if (mu > T(1e20)) break;
    }
  }
  return sqrtq(f);
}

// Gauss-Newton refinement for a consistent, rank-deficient system R(x) = 0: each step is
// the basic solution of J d = -r from Householder QR with column pivoting
// (Businger-Golub), rank decided relative to the largest pivot
template<typename T>
T gauss_newton(Vec<T>& x, int flags, int max_iter, T rank_tol) {
  Vec<T> r;
  std::vector<std::vector<T>> J;
  T nr(0), best_nr(-1);
  Vec<T> best_x = x;
  for (int it = 0; it < max_iter; ++it) {
    eval(x, flags, r, J);
    nr = sqrtq(sq_norm(r));
    if (best_nr < T(0) || nr < best_nr) { best_nr = nr; best_x = x; }
    if (nr == T(0)) break;
    int m = J.size(), n = NV, rank = 0;
    std::vector<int> perm(n);
    for (int j = 0; j < n; ++j) perm[j] = j;
    Vec<T> rhs(m);
    for (int i = 0; i < m; ++i) rhs[i] = -r[i];
    T pivot0(0);
    for (int k = 0; k < std::min(m, n); ++k) {
      // Column with largest remaining norm
      int jmax = k;
      T best(-1);
      for (int j = k; j < n; ++j) {
        T s(0);
        for (int i = k; i < m; ++i) s += J[i][j] * J[i][j];
        if (s > best) { best = s; jmax = j; }
      }
      if (k == 0) pivot0 = sqrtq(best);
      if (sqrtq(best) <= rank_tol * pivot0) break;
      for (int i = 0; i < m; ++i) std::swap(J[i][k], J[i][jmax]);
      std::swap(perm[k], perm[jmax]);
      // Householder reflection zeroing J[k+1:, k]
      T alpha = sqrtq(best);
      if (J[k][k] > 0) alpha = -alpha;
      Vec<T> v(m, T(0));
      for (int i = k; i < m; ++i) v[i] = J[i][k];
      v[k] -= alpha;
      T vv(0);
      for (int i = k; i < m; ++i) vv += v[i] * v[i];
      if (vv > T(0)) {
        for (int j = k; j < n; ++j) {
          T s(0);
          for (int i = k; i < m; ++i) s += v[i] * J[i][j];
          s = 2 * s / vv;
          for (int i = k; i < m; ++i) J[i][j] -= s * v[i];
        }
        T s(0);
        for (int i = k; i < m; ++i) s += v[i] * rhs[i];
        s = 2 * s / vv;
        for (int i = k; i < m; ++i) rhs[i] -= s * v[i];
      }
      rank = k + 1;
    }
    // Back substitution on the leading rank x rank triangle
    Vec<T> y(n, T(0));
    for (int i = rank - 1; i >= 0; --i) {
      T s = rhs[i];
      for (int j = i + 1; j < rank; ++j) s -= J[i][j] * y[j];
      y[i] = s / J[i][i];
    }
    for (int j = 0; j < n; ++j) x[perm[j]] += y[j];
    printf("  refinement %2d: |R| = %.3e, numerical rank %d of %d\n", it, double(nr), rank, n);
  }
  x = best_x;
  return best_nr;
}

// --------------------------------------------------------------- dense output
// Continuous extension y(t0 + theta h) = y0 + sum_i b_i(theta) k_i with polynomial weights
// b_i(theta) = sum_{q=1..ND} beta_qi theta^q: ODE order conditions up to order 4 for every
// theta, b(1) = b; the remaining freedom minimizes the integral over theta in [0, 1] of the
// squared order-5 ODE residuals and index-1 DAE residuals (No 18-22)
constexpr int ND = 4;

// Indices of linearly independent rows of E: Householder QR of E^T with column pivoting
template<typename T>
std::vector<int> independent_rows(const std::vector<Vec<T>>& E, T rank_tol) {
  int m = E.size(), n = E[0].size();
  std::vector<Vec<T>> M(n, Vec<T>(m));
  for (int i = 0; i < m; ++i) for (int j = 0; j < n; ++j) M[j][i] = E[i][j];
  std::vector<int> perm(m);
  for (int i = 0; i < m; ++i) perm[i] = i;
  T pivot0(0);
  int rank = 0;
  for (int k = 0; k < std::min(m, n); ++k) {
    int jmax = k;
    T best(-1);
    for (int j = k; j < m; ++j) {
      T s(0);
      for (int i = k; i < n; ++i) s += M[i][j] * M[i][j];
      if (s > best) { best = s; jmax = j; }
    }
    if (k == 0) pivot0 = sqrtq(best);
    if (sqrtq(best) <= rank_tol * pivot0) break;
    for (int i = 0; i < n; ++i) std::swap(M[i][k], M[i][jmax]);
    std::swap(perm[k], perm[jmax]);
    T alpha = sqrtq(best);
    if (M[k][k] > 0) alpha = -alpha;
    Vec<T> v(n, T(0));
    for (int i = k; i < n; ++i) v[i] = M[i][k];
    v[k] -= alpha;
    T vv(0);
    for (int i = k; i < n; ++i) vv += v[i] * v[i];
    if (vv > T(0)) {
      for (int j = k; j < m; ++j) {
        T s(0);
        for (int i = k; i < n; ++i) s += v[i] * M[i][j];
        s = 2 * s / vv;
        for (int i = k; i < n; ++i) M[i][j] -= s * v[i];
      }
    }
    rank = k + 1;
  }
  return std::vector<int>(perm.begin(), perm.begin() + rank);
}

// Solve M x = r by Gaussian elimination with partial pivoting
template<typename T>
Vec<T> gauss_solve(std::vector<Vec<T>> M, Vec<T> r) {
  int n = r.size();
  for (int k = 0; k < n; ++k) {
    int p = k;
    for (int i = k + 1; i < n; ++i) if (fabsq(M[i][k]) > fabsq(M[p][k])) p = i;
    std::swap(M[k], M[p]);
    std::swap(r[k], r[p]);
    for (int i = k + 1; i < n; ++i) {
      T f = M[i][k] / M[k][k];
      for (int j = k; j < n; ++j) M[i][j] -= f * M[k][j];
      r[i] -= f * r[k];
    }
  }
  Vec<T> x(n);
  for (int i = n - 1; i >= 0; --i) {
    T s = r[i];
    for (int j = i + 1; j < n; ++j) s -= M[i][j] * x[j];
    x[i] = s / M[i][i];
  }
  return x;
}

// Dense output weights beta[q][i], q = 0..ND-1 for theta^(q+1)
template<typename T>
std::vector<Vec<T>> dense_output(const Method<T>& m) {
  int n = ND * S;
  // Conditions: vector, right-hand side, power of theta
  struct Cond { Vec<T> v; T r; int o; };
  std::vector<Cond> eq, obj;
  for (int p = 1; p <= 4; ++p) {
    for (auto& t : trees(p)) eq.push_back({phi(t, m, S), T(1.) / T(double(tree_density(t))), p});
  }
  for (auto& t : trees(5)) obj.push_back({phi(t, m, S), T(1.) / T(double(tree_density(t))), 5});
  std::vector<std::pair<Vec<T>, Ratio>> V;
  dae_vectors(m, S, V);
  const int dae_order[5] = {3, 4, 4, 4, 4};
  for (int no = 18; no <= 22; ++no) obj.push_back({V[no].first, ratio<T>(V[no].second), dae_order[no - 18]});
  // Equality constraints E x = e: each condition for every power of theta, and b(1) = b
  std::vector<Vec<T>> E;
  Vec<T> e;
  for (auto& c : eq) {
    for (int q = 1; q <= ND; ++q) {
      Vec<T> row(n, T(0));
      for (int i = 0; i < S; ++i) row[(q - 1) * S + i] = c.v[i];
      E.push_back(row);
      e.push_back(q == c.o ? c.r : T(0));
    }
  }
  for (int i = 0; i < S; ++i) {
    Vec<T> row(n, T(0));
    for (int q = 0; q < ND; ++q) row[q * S + i] = T(1);
    E.push_back(row);
    e.push_back(m.b[i]);
  }
  std::vector<int> keep = independent_rows(E, T(1e-20));
  // Objective: exact integral over [0, 1] of the squared residuals, plus a small regularization
  std::vector<Vec<T>> H(n, Vec<T>(n, T(0)));
  Vec<T> g(n, T(0));
  for (auto& c : obj) {
    for (int q = 1; q <= ND; ++q) {
      for (int i = 0; i < S; ++i) {
        g[(q - 1) * S + i] += c.r * c.v[i] / T(double(q + c.o + 1));
        for (int qq = 1; qq <= ND; ++qq) {
          for (int j = 0; j < S; ++j) {
            H[(q - 1) * S + i][(qq - 1) * S + j] += c.v[i] * c.v[j] / T(double(q + qq + 1));
          }
        }
      }
    }
  }
  for (int i = 0; i < n; ++i) H[i][i] += T(1e-8);
  // KKT system [H E^T; E 0] [x; mu] = [g; e] with the independent rows of E
  int nk = n + keep.size();
  std::vector<Vec<T>> K(nk, Vec<T>(nk, T(0)));
  Vec<T> rk(nk, T(0));
  for (int i = 0; i < n; ++i) { for (int j = 0; j < n; ++j) K[i][j] = H[i][j]; rk[i] = g[i]; }
  for (size_t k = 0; k < keep.size(); ++k) {
    for (int j = 0; j < n; ++j) { K[n + k][j] = E[keep[k]][j]; K[j][n + k] = E[keep[k]][j]; }
    rk[n + k] = e[keep[k]];
  }
  Vec<T> x = gauss_solve(K, rk);
  // Check all equality constraints
  T worst(0);
  for (size_t k = 0; k < E.size(); ++k) worst = std::max(worst, fabsq(dot(E[k], Vec<T>(x.begin(), x.begin() + n)) - e[k]));
  printf("dense output: %zu of %zu conditions independent, residual %.3e\n", keep.size(), E.size(), double(worst));
  std::vector<Vec<T>> beta(ND, Vec<T>(S));
  for (int q = 0; q < ND; ++q) for (int i = 0; i < S; ++i) beta[q][i] = x[q * S + i];
  return beta;
}

// ----------------------------------------------------------------- stability
// R(z) = 1 + z b^T (I - z B)^{-1} 1
double stability_on_imaginary_axis(const Method<double>& m, const Vec<double>& b, int n) {
  double worst = 0;
  for (int k = -4000; k <= 4000; ++k) {
    double y = std::copysign(std::pow(10., std::abs(k) / 1000.), k) * (k == 0 ? 0 : 1);
    // Solve (I - z B) v = 1 in complex arithmetic, z = i y
    double vr[S], vi[S];
    for (int i = 0; i < n; ++i) {
      double ar = 1, ai = 0;  // rhs
      for (int j = 0; j < i; ++j) {
        // + z B_ij v_j, z = i y: (i y)(vr + i vi) = -y vi + i y vr
        ar += m.B[i][j] * (-y * vi[j]);
        ai += m.B[i][j] * (y * vr[j]);
      }
      // divide by 1 - i y B_ii
      double dr = 1, di = -y * m.B[i][i], den = dr * dr + di * di;
      vr[i] = (ar * dr + ai * di) / den;
      vi[i] = (ai * dr - ar * di) / den;
    }
    double sr = 0, si = 0;
    for (int i = 0; i < n; ++i) { sr += b[i] * vr[i]; si += b[i] * vi[i]; }
    // R = 1 + i y (sr + i si)
    double Rr = 1 - y * si, Ri = y * sr;
    worst = std::max(worst, std::sqrt(Rr * Rr + Ri * Ri));
  }
  return worst;
}

// --------------------------------------------------------- optimization on the manifold
// Orthonormal basis of the null space of J (m x NV), from Householder QR of J^T
// with column pivoting; rank decided relative to the largest pivot
std::vector<Vec<double>> null_space(std::vector<std::vector<double>> J, double rank_tol) {
  int m = J.size(), n = NV;
  // M = J^T (n x m), Q accumulated as n x n
  std::vector<Vec<double>> M(n, Vec<double>(m)), Q(n, Vec<double>(n, 0.));
  for (int i = 0; i < n; ++i) { Q[i][i] = 1; for (int k = 0; k < m; ++k) M[i][k] = J[k][i]; }
  int rank = 0;
  double pivot0 = 0;
  for (int k = 0; k < std::min(m, n); ++k) {
    int jmax = k;
    double best = -1;
    for (int j = k; j < m; ++j) {
      double s = 0;
      for (int i = k; i < n; ++i) s += M[i][j] * M[i][j];
      if (s > best) { best = s; jmax = j; }
    }
    if (k == 0) pivot0 = std::sqrt(best);
    if (std::sqrt(best) <= rank_tol * pivot0) break;
    for (int i = 0; i < n; ++i) std::swap(M[i][k], M[i][jmax]);
    double alpha = M[k][k] > 0 ? -std::sqrt(best) : std::sqrt(best);
    Vec<double> v(n, 0.);
    for (int i = k; i < n; ++i) v[i] = M[i][k];
    v[k] -= alpha;
    double vv = 0;
    for (int i = k; i < n; ++i) vv += v[i] * v[i];
    if (vv > 0) {
      for (int j = k; j < m; ++j) {
        double s = 0;
        for (int i = k; i < n; ++i) s += v[i] * M[i][j];
        s = 2 * s / vv;
        for (int i = k; i < n; ++i) M[i][j] -= s * v[i];
      }
      // Q <- Q H
      for (int i = 0; i < n; ++i) {
        double s = 0;
        for (int l = k; l < n; ++l) s += Q[i][l] * v[l];
        s = 2 * s / vv;
        for (int l = k; l < n; ++l) Q[i][l] -= s * v[l];
      }
    }
    rank = k + 1;
  }
  std::vector<Vec<double>> N;
  for (int j = rank; j < n; ++j) {
    Vec<double> col(n);
    for (int i = 0; i < n; ++i) col[i] = Q[i][j];
    N.push_back(col);
  }
  return N;
}

// Error coefficients and their Jacobian
void eval_error(const Vec<double>& x, Vec<double>& e, std::vector<std::vector<double>>& Je) {
  std::vector<Dual<double>> xd(NV);
  for (int i = 0; i < NV; ++i) xd[i] = Dual<double>::var(x[i], i);
  Vec<Dual<double>> ed = error_terms(xd.data());
  e.resize(ed.size());
  Je.assign(ed.size(), Vec<double>(NV));
  for (size_t k = 0; k < ed.size(); ++k) {
    e[k] = ed[k].v;
    for (int i = 0; i < NV; ++i) Je[k][i] = ed[k].d[i];
  }
}

double stability_on_imaginary_axis(const Method<double>& m, const Vec<double>& b, int n);

// Minimize |error coefficients| over the solution manifold, keeping A-stability
double optimize_on_manifold(Vec<double>& x, int flags, int max_iter) {
  Vec<double> r, e, xn;
  std::vector<std::vector<double>> J, Je;
  auto a_stable = [&](const Vec<double>& y) {
    Method<double> m = assemble(y.data());
    return stability_on_imaginary_axis(m, m.b, S) <= 1 + 1e-12
      && stability_on_imaginary_axis(m, m.bh, 7) <= 1 + 1e-12;
  };
  eval_error(x, e, Je);
  double fe = std::sqrt(sq_norm(e)), radius = 0.05;
  for (int it = 0; it < max_iter && radius > 1e-10; ++it) {
    eval(x, flags, r, J);
    std::vector<Vec<double>> N = null_space(J, 1e-9);
    int nn = N.size();
    // Reduced Gauss-Newton: min |e + Je N z|, damped to the trust radius
    std::vector<Vec<double>> G(e.size(), Vec<double>(nn, 0.));
    for (size_t k = 0; k < e.size(); ++k)
      for (int j = 0; j < nn; ++j)
        for (int i = 0; i < NV; ++i) G[k][j] += Je[k][i] * N[j][i];
    Vec<double> z(nn, 0.), g(nn, 0.);
    std::vector<Vec<double>> H(nn, Vec<double>(nn, 0.));
    for (size_t k = 0; k < e.size(); ++k)
      for (int i = 0; i < nn; ++i) {
        g[i] -= G[k][i] * e[k];
        for (int j = 0; j < nn; ++j) H[i][j] += G[k][i] * G[k][j];
      }
    double mu = 1e-8;
    for (int tries = 0; tries < 60; ++tries, mu *= 10) {
      // Solve (H + mu I) z = g by Gaussian elimination (small)
      std::vector<Vec<double>> Hm = H;
      Vec<double> rhs = g;
      for (int i = 0; i < nn; ++i) Hm[i][i] += mu;
      for (int k = 0; k < nn; ++k)
        for (int i = k + 1; i < nn; ++i) {
          double f = Hm[i][k] / Hm[k][k];
          for (int j = k; j < nn; ++j) Hm[i][j] -= f * Hm[k][j];
          rhs[i] -= f * rhs[k];
        }
      for (int i = nn - 1; i >= 0; --i) {
        z[i] = rhs[i];
        for (int j = i + 1; j < nn; ++j) z[i] -= Hm[i][j] * z[j];
        z[i] /= Hm[i][i];
      }
      if (std::sqrt(sq_norm(z)) <= radius) break;
    }
    xn = x;
    for (int j = 0; j < nn; ++j) for (int i = 0; i < NV; ++i) xn[i] += N[j][i] * z[j];
    // Restore feasibility
    double res = levenberg_marquardt<double>(xn, flags, 200, 1e-14);
    Vec<double> en;
    std::vector<std::vector<double>> Jen;
    eval_error(xn, en, Jen);
    double fn = std::sqrt(sq_norm(en));
    if (res < 1e-12 && fn < fe && a_stable(xn)) {
      x = xn; e = en; Je = Jen; fe = fn;
      radius = std::min(1.0, 2 * radius);
    } else {
      radius /= 4;
    }
  }
  return fe;
}

// ----------------------------------------------------------------------- driver
std::string q2s(quad v) { char buf[64]; quadmath_snprintf(buf, sizeof(buf), "%+.33Qe", v); return buf; }

int main(int argc, char** argv) {
  int n_starts = argc > 1 ? atoi(argv[1]) : 256;
  // Single start: argv[1] < 0 selects seed -argv[1]
  int only = n_starts < 0 ? -n_starts : -1;
  if (only >= 0) n_starts = only + 1;
  int seed = argc > 2 ? atoi(argv[2]) : 0;
  int flags = argc > 3 ? atoi(argv[3]) : (WITH_INDEX2 | WITH_W | WITH_PR);
  int opt_iter = argc > 4 ? atoi(argv[4]) : 300;
  std::mutex mtx;
  struct Candidate { double res, err6, stab, stab_emb; Vec<double> x; int k; };
  std::vector<Candidate> found;
  int n_threads = std::max(1u, std::thread::hardware_concurrency());
  std::vector<std::thread> pool;
  for (int t = 0; t < n_threads; ++t) {
    pool.emplace_back([&, t]() {
      for (int k = t; k < n_starts; k += n_threads) {
        if (only >= 0 && k != only) continue;
        std::mt19937_64 rng(seed + k);
        std::normal_distribution<double> nd(0, 0.5);
        std::uniform_real_distribution<double> ud(0.15, 0.35);
        Vec<double> x(NV);
        for (auto& e : x) e = nd(rng);
        x[0] = ud(rng);
        double res = levenberg_marquardt<double>(x, flags, 2000, 1e-14);
        if (!(res < 1e-11)) continue;
        Method<double> m = assemble(x.data());
        if (opt_iter > 0 && stability_on_imaginary_axis(m, m.b, S) <= 1 + 1e-12
            && stability_on_imaginary_axis(m, m.bh, 7) <= 1 + 1e-12) {
          // Use the remaining freedom to reduce the error coefficients
          optimize_on_manifold(x, flags, opt_iter);
          res = levenberg_marquardt<double>(x, flags, 200, 1e-14);
          m = assemble(x.data());
        }
        double e6 = std::sqrt(sq_norm(error_terms(x.data())));
        double st = stability_on_imaginary_axis(m, m.b, S);
        double ste = stability_on_imaginary_axis(m, m.bh, 7);
        std::lock_guard<std::mutex> lock(mtx);
        found.push_back({res, e6, st, ste, x, k});
        printf("start %4d: residual %.1e gamma %.6f |R(iy)|max %.6f embedded %.6f err6 %.3e\n",
          k, res, x[0], st, ste, e6);
        fflush(stdout);
      }
    });
  }
  for (auto& th : pool) th.join();
  // Best A-stable candidate (main and embedded), smallest error constants
  const Candidate* best = nullptr;
  for (auto& c : found) {
    if (c.stab > 1 + 1e-12 || c.stab_emb > 1 + 1e-12) continue;
    if (!best || c.err6 < best->err6) best = &c;
  }
  printf("%zu solutions, %s A-stable one selected\n", found.size(), best ? "an" : "no");
  if (!best) return 1;
  int best_k = best->k;
  // Refinement in quadruple precision
  Vec<quad> xq(NV);
  for (int i = 0; i < NV; ++i) xq[i] = best->x[i];
  quad res = gauss_newton<quad>(xq, flags, 12, quad(1e-12));
  printf("quad residual %s\n", q2s(res).c_str());
  Method<quad> m = assemble(xq.data());
  // Header with the coefficients in the original form
  FILE* out = fopen("rosenbrock_tableau.hpp", "w");
  auto L = [](quad v) { return q2s(v) + "L"; };
  // License header of casadi/solvers
  fputs(
    "/*\n"
    " *    This file is part of CasADi.\n"
    " *\n"
    " *    CasADi -- A symbolic framework for dynamic optimization.\n"
    " *    Copyright (C) 2010-2026 Joel Andersson, Joris Gillis, Moritz Diehl,\n"
    " *                            KU Leuven. All rights reserved.\n"
    " *    Copyright (C) 2011-2014 Greg Horn\n"
    " *\n"
    " *    CasADi is free software; you can redistribute it and/or\n"
    " *    modify it under the terms of the GNU Lesser General Public\n"
    " *    License as published by the Free Software Foundation; either\n"
    " *    version 3 of the License, or (at your option) any later version.\n"
    " *\n"
    " *    CasADi is distributed in the hope that it will be useful,\n"
    " *    but WITHOUT ANY WARRANTY; without even the implied warranty of\n"
    " *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU\n"
    " *    Lesser General Public License for more details.\n"
    " *\n"
    " *    You should have received a copy of the GNU Lesser General Public\n"
    " *    License along with CasADi; if not, write to the Free Software\n"
    " *    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA\n"
    " *\n"
    " */\n"
    "\n"
    "\n", out);
  fprintf(out, "// Generated by misc/rosenbrock_coefficients.cpp (start %d, seed %d): do not edit.\n", best_k, seed);
  fprintf(out, "// Order conditions residual %.3e, |R(iy)| <= 1 (embedded: <= 1),\n", double(res));
  fprintf(out, "// principal error coefficients (2-norm) %.4e\n", double(std::sqrt(sq_norm(error_terms(best->x.data())))));
  fprintf(out, "#ifndef CASADI_ROSENBROCK_TABLEAU_HPP\n#define CASADI_ROSENBROCK_TABLEAU_HPP\n\n");
  fprintf(out, "namespace casadi {\nnamespace rosenbrock_tableau {\n\n");
  fprintf(out, "// Number of stages, order of the embedded method\n");
  fprintf(out, "const casadi_int n_stages = %d;\nconst casadi_int order_embedded = 4;\n\n", S);
  fprintf(out, "// Diagonal coefficient\nconst long double gamma = %s;\n\n", L(m.gamma).c_str());
  const char* names[2] = {"alpha", "gamma_ij"};
  for (int w = 0; w < 2; ++w) {
    fprintf(out, "// %s, strictly lower triangular\nconst long double %s[%d][%d] = {\n",
      w == 0 ? "Stage coefficients alpha_ij" : "Jacobian coefficients gamma_ij", names[w], S, S);
    for (int i = 0; i < S; ++i) {
      fprintf(out, "  {");
      for (int j = 0; j < S; ++j) {
        quad v = j < i ? (w == 0 ? m.A[i][j] : m.B[i][j] - m.A[i][j]) : quad(0);
        fprintf(out, "%s%s", j < i ? L(v).c_str() : "0", j + 1 < S ? (j + 1 < i ? ",\n   " : ", ") : "");
      }
      fprintf(out, "}%s\n", i + 1 < S ? "," : "");
    }
    fprintf(out, "};\n\n");
  }
  fprintf(out, "// Weights of the method and of the embedded method\nconst long double b[%d] = {\n", S);
  for (int i = 0; i < S; ++i) fprintf(out, "  %s%s\n", L(m.b[i]).c_str(), i + 1 < S ? "," : "");
  fprintf(out, "};\nconst long double bhat[%d] = {\n", S);
  for (int i = 0; i < S; ++i) fprintf(out, "  %s%s\n", L(m.bh[i]).c_str(), i + 1 < S ? "," : "");
  std::vector<Vec<quad>> beta = dense_output(m);
  fprintf(out, "};\n\n// Dense output y(t0 + theta h) = y0 + sum_i b_i(theta) k_i of order 4 for ODEs,\n"
    "// b_i(theta) = sum_q b_dense[q][i] theta^(q+1), b_i(1) = b_i\n");
  fprintf(out, "const casadi_int n_dense = %d;\nconst long double b_dense[%d][%d] = {\n", ND, ND, S);
  for (int q = 0; q < ND; ++q) {
    fprintf(out, "  {");
    for (int i = 0; i < S; ++i) fprintf(out, "%s%s", L(beta[q][i]).c_str(), i + 1 < S ? ",\n   " : "");
    fprintf(out, "}%s\n", q + 1 < ND ? "," : "");
  }
  fprintf(out, "};\n\n} // namespace rosenbrock_tableau\n} // namespace casadi\n\n#endif // CASADI_ROSENBROCK_TABLEAU_HPP\n");
  fclose(out);
  printf("gamma = %s, written to rosenbrock_tableau.hpp\n", q2s(m.gamma).c_str());
  return 0;
}

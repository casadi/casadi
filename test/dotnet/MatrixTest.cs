using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;

namespace CasADi.Tests {

/// <summary>Port of test/java MatrixTest (test/python/matrix.py)</summary>
public class MatrixTest : CasadiTestCase {
  const long NONE_START = long.MinValue, NONE_STOP = long.MaxValue;

  static DM idx(params long[] i) => i;

  static DM range(Sparsity sp) {
    var v = new double[sp.nnz()];
    for (int k = 0; k < v.Length; k++) v[k] = k;
    return ca.DM(sp, v);
  }

  static DM row(params double[] x) => ((DM)x).T();

  long[] permutation(int n) {
    var p = new long[n];
    for (int i = 0; i < n; i++) p[i] = i;
    for (int i = n - 1; i > 0; i--) {
      int j = random.Next(i + 1);
      (p[i], p[j]) = (p[j], p[i]);
    }
    return p;
  }

  static double npDet(DM m) {
    double[,] a = m.ToArray();
    int n = a.GetLength(0);
    double d = 1;
    for (int k = 0; k < n; k++) {
      int p = k;
      for (int i = k + 1; i < n; i++) if (Math.Abs(a[i, k]) > Math.Abs(a[p, k])) p = i;
      if (a[p, k] == 0) return 0;
      if (p != k) {
        for (int j = 0; j < n; j++) (a[p, j], a[k, j]) = (a[k, j], a[p, j]);
        d = -d;
      }
      d *= a[k, k];
      for (int i = k + 1; i < n; i++) {
        double f = a[i, k] / a[k, k];
        for (int j = k; j < n; j++) a[i, j] -= f * a[k, j];
      }
    }
    return d;
  }

  static double[,] npSolve(double[,] A0, double[,] B0) {
    int n = A0.GetLength(0), m = B0.GetLength(1);
    var a = (double[,])A0.Clone();
    var b = (double[,])B0.Clone();
    for (int k = 0; k < n; k++) {
      int p = k;
      for (int i = k + 1; i < n; i++) if (Math.Abs(a[i, k]) > Math.Abs(a[p, k])) p = i;
      for (int j = 0; j < n; j++) (a[p, j], a[k, j]) = (a[k, j], a[p, j]);
      for (int j = 0; j < m; j++) (b[p, j], b[k, j]) = (b[k, j], b[p, j]);
      for (int i = k + 1; i < n; i++) {
        double f = a[i, k] / a[k, k];
        for (int j = k; j < n; j++) a[i, j] -= f * a[k, j];
        for (int j = 0; j < m; j++) b[i, j] -= f * b[k, j];
      }
    }
    var x = new double[n, m];
    for (int i = n - 1; i >= 0; i--)
      for (int j = 0; j < m; j++) {
        double s = b[i, j];
        for (int k = i + 1; k < n; k++) s -= a[i, k] * x[k, j];
        x[i, j] = s / a[i, i];
      }
    return x;
  }

  static double[,] npKron(double[,] a, double[,] b) {
    int n = a.GetLength(0), m = a.GetLength(1), p = b.GetLength(0), q = b.GetLength(1);
    var r = new double[n * p, m * q];
    for (int i = 0; i < n; i++)
      for (int j = 0; j < m; j++)
        for (int k = 0; k < p; k++)
          for (int l = 0; l < q; l++) r[i * p + k, j * q + l] = a[i, j] * b[k, l];
    return r;
  }

  double[,] npRandom(int n, int m) {
    var r = new double[n, m];
    for (int i = 0; i < n; i++)
      for (int j = 0; j < m; j++) r[i, j] = random.NextDouble();
    return r;
  }

  /// <summary>Python helper check_codegen: generate C, compile with the shell importer, compare; optionally serialize</summary>
  void checkCodegen(Function f, List<DM> inputs, bool checkSerialize) {
    if (checkSerialize) {
      Function g = Function.deserialize(f.serialize());
      List<DM> a = f.call(inputs), b = g.call(inputs);
      for (int i = 0; i < a.Count; i++) checkarray(a[i], b[i], 1e-10, "serialize output " + i);
    }
    if (!RunSlow) return;
    string dir = Directory.CreateTempSubdirectory("casadi_cg").FullName;
    string name = "cg_" + f.name() + "_" + random.Next().ToString("x");
    var cg = ca.CodeGenerator(name + ".c");
    cg.add(f);
    string src = cg.generate(dir + "/");
    Function F = ca.external(f.name(), ca.Importer(src, "shell"));
    List<DM> r0 = f.call(inputs), r1 = F.call(inputs);
    for (int i = 0; i < r0.Count; i++) checkarray(r0[i], r1[i], 1e-10, "codegen output " + i);
  }

  public void test_constructorlol() {
    message("List of list constructor");
    DM a = new double[,] { { 1, 2, 3 }, { 4, 5, 6 }, { 7, 8, 9 } };
    var b = ca.DM(3, 3);
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++) b[i, j] = 3 * i + j + 1;
    checkarray(a, b, "List of list constructor");
  }

  public void test_sum() {
    message("sum");
    DM D = ca.DM(new double[,] { { 1, 2, 3 }, { 4, 5, 6 }, { 7, 8, 9 } });
    checkarray(ca.sum1(D), new double[,] { { 12, 15, 18 } });
    checkarray(ca.sum2(D), new double[,] { { 6 }, { 15 }, { 24 } });
  }

  public void test_inv() {
    message("Matrix inverse");
    DM a = ca.DM(new double[,] { { 1, 2 }, { 1, 3 } });
    checkarray(ca.mtimes(ca.inv(a), a), DM.eye(2), "DM inverse");
  }

  public void test_trans() {
    message("trans");
    var a = ca.DM(0, 1);
    DM b = a.T();
    assertEquals(1L, b.size1());
    assertEquals(0L, b.size2());
  }

  public void test_vertcat() {
    message("vertcat");
    DM A = DM.ones(2, 3);
    var B = ca.DM(4, 3);
    DM C = ca.vertcat(A, B);
    assertEquals((6L, 3L), C.size(), "vertcat shape");
    assertEquals(A.nnz(), C.nnz(), "vertcat size");
    Check.Throws<Exception>(() => ca.horzcat(A, B));
  }

  public void test_horzcat() {
    message("horcat");
    DM A = DM.ones(3, 2);
    var B = ca.DM(3, 4);
    DM C = ca.horzcat(A, B);
    assertEquals((3L, 6L), C.size(), "horzcat shape");
    assertEquals(A.nnz(), C.nnz(), "vertcat size");
    Check.Throws<Exception>(() => ca.vertcat(A, B));
  }

  /// <summary>First test_diagcat of matrix.py, shadowed there by the later definition</summary>
  public void test_diagcat_shadowed() {
    MX x = MX.sym("x", 2, 2);
    MX y = MX.sym("y", Sparsity.lower(3));
    MX z = MX.sym("z", 4, 2);
    var L = new List<MX> { x, y, z };
    var fMX = ca.Function("fMX", L, new[] { ca.diagcat(L) });
    var LSX = L.Select(i => SX.sym("", i.sparsity())).ToList();
    var fSX = ca.Function("fSX", LSX, new[] { ca.diagcat(LSX) });
    var inputs = new List<DM>();
    for (long i = 0; i < 3; i++) inputs.Add(range(fMX.sparsity_in(i)));
    checkfunction(fMX, fSX, inputs);
  }

  public void test_veccat() {
    message("vecccat");
    var A = ca.DM(2, 3);
    A[0, 1] = 2;
    A[1, 0] = 1;
    A[1, 2] = 3;
    var B = ca.DM(3, 1);
    B[0, 0] = 4;
    B[1, 0] = 5;
    B[2, 0] = 6;
    DM C = ca.veccat(A, B);
    assertEquals((9L, 1L), C.size(), "veccat shape");
    assertEquals(A.nnz() + B.nnz(), C.nnz(), "veccat size");
    Check.True(new double[] { 1, 2, 3, 4, 5, 6 }.SequenceEqual(C.nonzeros()), "numbers shape");
  }

  /// <summary>Python list slice a1[start:stop:step], with NONE_START/NONE_STOP for omitted bounds</summary>
  static DM pySlice(double[] a, long start, long stop, long step) {
    int n = a.Length;
    long s, e;
    if (step > 0) {
      s = start == NONE_START ? 0 : start < 0 ? Math.Max(0, start + n) : Math.Min(start, n);
      e = stop == NONE_STOP ? n : stop < 0 ? Math.Max(0, stop + n) : Math.Min(stop, n);
    } else {
      s = start == NONE_START ? n - 1 : start < 0 ? Math.Max(-1, start + n) : Math.Min(start, n - 1);
      e = stop == NONE_STOP ? -1 : stop < 0 ? Math.Max(-1, stop + n) : Math.Min(stop, n - 1);
    }
    var r = new List<double>();
    for (long i = s; step > 0 ? i < e : i > e; i += step) r.Add(a[i]);
    return r.ToArray();
  }

  public void test_slicestepnegative() {
    message("Slice step negative");
    double[] a1 = { 1, 2, 3, 4, 5 };
    DM a2 = a1;
    long[,] cases = { { 0, 4, -1 }, { 4, 0, -1 }, { 0, 4, -2 }, { 4, 0, -2 }, { 1, 4, -2 }, { 4, 1, -2 }, { 0, 3, -2 },
                      { 3, 0, -2 }, { NONE_START, NONE_STOP, -1 }, { NONE_START, NONE_STOP, 1 }, { 2, NONE_STOP, -1 },
                      { NONE_START, 2, -1 }, { -1, NONE_STOP, -1 }, { NONE_START, -1, -1 } };
    for (int k = 0; k < cases.GetLength(0); k++) {
      long c0 = cases[k, 0], c1 = cases[k, 1], c2 = cases[k, 2];
      checkarray(a2[ca.Slice(c0, c1, c2), 0], pySlice(a1, c0, c1, c2), "[" + c0 + ":" + c1 + ":" + c2 + "]");
    }
  }

  public void test_indexingOutOfBounds() {
    message("Indexing out of bounds");
    DM y = DM.zeros(4, 5);
    Check.Throws<Exception>(() => _ = y[12, 0]);
    Check.Throws<Exception>(() => _ = y[12, 12]);
    Check.Throws<Exception>(() => _ = y[0, 12]);
    Check.Throws<Exception>(() => _ = y[12, ..]);
    Check.Throws<Exception>(() => _ = y[12..15, 0]);
    Check.Throws<Exception>(() => _ = y[.., 12]);
    Check.Throws<Exception>(() => _ = y[0, 12..15]);
    _ = y[-1, 2];
    Check.Throws<Exception>(() => _ = y[-12, 2]);
    _ = y[^3..^1, 2];
    Check.Throws<Exception>(() => _ = y[^12..^9, 2]);

    Check.Throws<Exception>(() => y[12, 0] = 0);
    Check.Throws<Exception>(() => y[12, 12] = 0);
    Check.Throws<Exception>(() => y[0, 12] = 0);
    Check.Throws<Exception>(() => y[12, ..] = 0);
    Check.Throws<Exception>(() => y[12..15, 0] = 0);
    Check.Throws<Exception>(() => y[.., 12] = 0);
    Check.Throws<Exception>(() => y[0, 12..15] = 0);
    y[-1, 2] = 0;
    Check.Throws<Exception>(() => y[-12, 2] = 0);
    y[^3..^1, 2] = 0;
    Check.Throws<Exception>(() => y[^12..^9, 2] = 0);
  }

  // huge_slice: not a test method in Python (no test_ prefix)

  public void test_nonmonotonous_indexing() {
    message("non-monotonous indexing");
    DM A = ca.DM(new double[,] { { 1, 2, 3 }, { 4, 5, 6 }, { 7, 8, 9 } });
    DM B = A.get(idx(0, 2, 1), ca.Slice(0));
    checkarray(ca.vertcat(1, 7, 4), B, "non-monotonous");
    B = A.get(ca.Slice(0), idx(0, 2, 1));
    checkarray(row(1, 3, 2), B, "non-monotonous");
  }

  public void test_IM_indexing() {
    message("IM");
    var A = ca.DM(2, 2);
    A[0, 0] = 1;
    A[1, 1] = 3;
    A[0, 1] = 2;
    A[1, 0] = 4;
    DM B = ca.vertcat(1, 2, 3, 4, 5);
    DM B_ = B.get(A);
    checkarray(B_, ca.DM(new double[,] { { 2, 3 }, { 5, 4 } }), "matrix indexing");
    B.set(A, ca.DM(new double[,] { { 1, 2 }, { 3, 4 } }));
    checkarray(B, ca.vertcat(1, 1, 2, 4, 3), "matrix indexing assignement");

    B = ca.vertcat(1, 2, 3, 4, 5);
    B_ = B.get(A);
    checkarray(B_.ToArray(), ca.DM(new double[,] { { 2, 3 }, { 5, 4 } }), "Imatrix indexing");
    B.set(A, ca.DM(new double[,] { { 1, 2 }, { 3, 4 } }));
    checkarray(B.ToArray(), ca.vertcat(1, 1, 2, 4, 3), "Imatrix indexing assignement");

    var B2 = ca.DM(5, 1);
    Check.Throws<Exception>(() => B2.get_nz(false, A));
  }

  public void test_sparsity_indexing() {
    message("sparsity");
    DM B = ca.DM(new double[,] { { 1, 2, 3, 4, 5 }, { 6, 7, 8, 9, 10 } });
    DM A = ca.sparsify(ca.DM(new double[,] { { 1, 1, 0, 0, 0 }, { 0, 0, 1, 0, 0 } }));
    Sparsity sp = A.sparsity();
    checkarray(B.get(sp), ca.DM(new double[,] { { 1, 2, 0, 0, 0 }, { 0, 0, 8, 0, 0 } }), "sparsity indexing");
    B.set(sp, -4);
    checkarray(B, ca.DM(new double[,] { { -4, -4, 3, 4, 5 }, { 6, 7, -4, 9, 10 } }), "sparsity indexing assignement");
    DM B2 = ca.DM(new double[,] { { 1, 2, 3, 4, 5 }, { 6, 7, 8, 9, 10 } });
    B2.set(sp, B2 * 2);
    checkarray(B2, ca.DM(new double[,] { { 2, 4, 3, 4, 5 }, { 6, 7, 16, 9, 10 } }), "Imatrix indexing assignement");
    Check.Throws<Exception>(() => B2.get(Sparsity.dense(4, 4)));
  }

  public void test_index_setting() {
    message("index setting");
    DM B = ca.vertcat(1, 2, 3, 4, 5);
    B[0] = 8;
    checkarray(B, ca.vertcat(8, 2, 3, 4, 5), "index setting");
    B[1, 0] = 4;
    checkarray(B, ca.vertcat(8, 4, 3, 4, 5), "index setting");
    B[.., 0] = 7;
    checkarray(B, ca.vertcat(7, 7, 7, 7, 7), "index setting");
  }

  public void test_issue298() {
    message("Issue #298");
    var a = ca.DM(4, 1);
    object b = ca.reshape(a, 2, 2);
    assertEquals(a.GetType(), b.GetType());
  }

  public void test_det() {
    message("Determinant");
    var a = ca.DM(1, 1);
    a[0, 0] = 5;
    checkarray(ca.det(a), npDet(a));

    a = ca.DM(5, 5);
    for (int i = 0; i < 5; i++) a[i, i] = i + 1;
    checkarray(ca.det(a), npDet(a));

    a = ca.DM(5, 5);
    for (int i = 0; i < 4; i++) a[i, i] = i + 1;
    a[0, 4] = 3;
    a[4, 0] = 7;
    checkarray(ca.det(a), npDet(a));

    a = ca.DM(5, 5);
    for (int i = 0; i < 5; i++)
      for (int j = 0; j < 5; j++) a[i, j] = i + j;
    checkarray(ca.det(a), npDet(a));

    foreach (var (s0, s1) in new[] { (4, 5), (5, 4) }) {
      a = ca.DM(5, 5);
      for (int i = 0; i < s0; i++)
        for (int j = 0; j < s1; j++) a[i, j] = i + j;
      checkarray(ca.det(a), npDet(a));
    }
    foreach (var (e0, e1) in new[] { (4L, 1L), (1L, 4L), (4L, 2L), (2L, 4L) }) {
      a = ca.DM(5, 5);
      int ni = e0 == 4 ? 4 : 5, nj = e0 == 4 ? 5 : 4;
      for (int i = 0; i < ni; i++)
        for (int j = 0; j < nj; j++) a[i, j] = i + j;
      a[e0, e1] = 12;
      checkarray(ca.det(a), npDet(a));
    }

    a = ca.DM(50, 50);
    for (int i = 0; i < 50; i++) a[i, i] = i + 1;
    checkarray(ca.det(a) / npDet(a), 1);
  }

  public void test_inv_sparsity() {
    Check.Assume(GlobalOptions.getSimplificationOnTheFly(), "simplification on the fly disabled");
    message("sparsity pattern of inverse");
    int n = 8;
    Sparsity sp = Sparsity.lower(n);
    var syms = new List<SX>();
    for (int i = 0; i < sp.nnz(); i++) syms.Add(SX.sym("a" + i));
    var x = ca.SX(sp, ca.vertcat(syms));
    DM x_ = DM.ones(x.sparsity());
    DM I_ = DM.ones(ca.inv(x).sparsity());
    checkarray(x_, I_, "inv");

    sp = Sparsity.lower(n);
    x = SX.sym("a", sp);
    x[0, n - 1] = 1;
    I_ = DM.ones(ca.inv(x).sparsity());
    checkarray(DM.ones(n, n), I_, "inv");

    x = SX.sym("a", sp);
    x[0, n / 2] = 1;
    DM s_ = DM.ones(sp);
    s_[.., 0..(n / 2 + 1)] = 1;
    I_ = DM.ones(ca.inv_minor(x).sparsity());
    s_ = ca.densify(s_);
    ca.densify(I_);
    checkarray(s_, I_, "inv");
  }

  public void test_mtimes() {
    DM A = DM.ones((4, 3));
    DM B = DM.ones((3, 8));
    DM C = DM.ones((8, 7));
    Check.Throws<Exception>(() => ca.mtimes(new List<DM>()));
    DM D = ca.mtimes(new[] { A });
    assertEquals((4L, 3L), D.size());
    D = ca.mtimes(new[] { A, B });
    assertEquals((4L, 8L), D.size());
    D = ca.mtimes(new[] { A, B, C });
    assertEquals((4L, 7L), D.size());
  }

  public void test_remove() {
    message("remove");
    DM B = ca.DM(new double[,] { { 1, 2, 3, 4 }, { 5, 6, 7, 8 }, { 9, 10, 11, 12 }, { 13, 14, 15, 16 }, { 17, 18, 19, 20 } });
    var A = ca.DM(B);
    A.remove(new long[] { }, new long[] { });
    checkarray(A, B, "remove nothing");
    A = ca.DM(B);
    A.remove(new long[] { }, new long[] { 1 });
    checkarray(A, ca.DM(new double[,] { { 1, 3, 4 }, { 5, 7, 8 }, { 9, 11, 12 }, { 13, 15, 16 }, { 17, 19, 20 } }),
               "remove a column");
    A = ca.DM(B);
    A.remove(new long[] { 0, 3 }, new long[] { 1 });
    checkarray(A, ca.DM(new double[,] { { 5, 7, 8 }, { 9, 11, 12 }, { 17, 19, 20 } }), "remove a column and two rows ");
  }

  /// <summary>The six comparisons of A against c, both operand orders; reference = {&lt;=,&lt;,&gt;,&gt;=,==,!=} of A op c</summary>
  static void checkComparisons(DM A, double c, double[][,] reference) {
    var C = ca.DM(c);
    var lhsD = new List<Func<DM, double, DM>> { (a, b) => ca.le(a, b), (a, b) => ca.lt(a, b), (a, b) => ca.gt(a, b),
                                                (a, b) => ca.ge(a, b), (a, b) => ca.eq(a, b), (a, b) => ca.ne(a, b) };
    var lhsM = new List<Func<DM, DM, DM>> { (a, b) => a <= b, (a, b) => a < b, (a, b) => a > b,
                                            (a, b) => a >= b, (a, b) => ca.eq(a, b), (a, b) => ca.ne(a, b) };
    var rhsD = new List<Func<double, DM, DM>> { (b, a) => ca.ge(b, a), (b, a) => ca.gt(b, a), (b, a) => ca.lt(b, a),
                                                (b, a) => ca.le(b, a), (b, a) => ca.eq(b, a), (b, a) => ca.ne(b, a) };
    var rhsM = new List<Func<DM, DM, DM>> { (b, a) => b >= a, (b, a) => b > a, (b, a) => b < a,
                                            (b, a) => b <= a, (b, a) => ca.eq(b, a), (b, a) => ca.ne(b, a) };
    var lhsI = new List<Func<DM, double, DM>> { (a, b) => a <= b, (a, b) => a < b, (a, b) => a > b,
                                                (a, b) => a >= b, (a, b) => ca.eq(a, b), (a, b) => ca.ne(a, b) };
    string[] names = { "<=", "<", ">", ">=", "==", "!=" };
    for (int k = 0; k < 6; k++) {
      DM r = ca.vertcat(reference[k]);
      checkarray(lhsD[k](A, c), r, 1e-10, "A" + names[k] + c);
      checkarray(lhsM[k](A, C), r, 1e-10, "A" + names[k] + "DM(" + c + ")");
      checkarray(lhsI[k](A, c), r, 1e-10, "A." + names[k] + "(" + c + ")");
      checkarray(rhsD[k](c, A), r, 1e-10, c + " reversed " + names[k]);
      checkarray(rhsM[k](C, A), r, 1e-10, "DM(" + c + ") reversed " + names[k]);
    }
  }

  public void test_comparisons() {
    DM A = ca.DM(new double[,] { { 5, 4 }, { 2, 1 } });
    double[,] T = { { 1, 1 }, { 1, 1 } }, F = { { 0, 0 }, { 0, 0 } };
    // python-only: np.array / np.matrix operands
    foreach (double c in new double[] { 6, 6.0 }) checkComparisons(A, c, new[] { T, T, F, F, F, T });
    foreach (double c in new double[] { 5, 5.0 })
      checkComparisons(A, c, new[] { T, new double[,] { { 0, 1 }, { 1, 1 } }, F, new double[,] { { 1, 0 }, { 0, 0 } },
                                     new double[,] { { 1, 0 }, { 0, 0 } }, new double[,] { { 0, 1 }, { 1, 1 } } });
    foreach (double c in new double[] { 4, 4.0 })
      checkComparisons(A, c, new[] { new double[,] { { 0, 1 }, { 1, 1 } }, new double[,] { { 0, 0 }, { 1, 1 } },
                                     new double[,] { { 1, 0 }, { 0, 0 } }, new double[,] { { 1, 1 }, { 0, 0 } },
                                     new double[,] { { 0, 1 }, { 0, 0 } }, new double[,] { { 1, 0 }, { 1, 1 } } });
    foreach (double c in new double[] { 1, 1.0 })
      checkComparisons(A, c, new[] { new double[,] { { 0, 0 }, { 0, 1 } }, F, new double[,] { { 1, 1 }, { 1, 0 } }, T,
                                     new double[,] { { 0, 0 }, { 0, 1 } }, new double[,] { { 1, 1 }, { 1, 0 } } });
    checkComparisons(A, 0, new[] { F, F, T, T, F, T });
  }

  /// <summary>Python bool(x)</summary>
  static bool truth(DM x) => x.nonzero();

  public void test_truth() {
    assertTrue(truth(ca.vertcat(1)));
    assertFalse(truth(ca.vertcat(0)));
    assertTrue(truth(ca.vertcat(0.2)));
    assertTrue(truth(ca.vertcat(-0.2)));
    Check.Throws<Exception>(() => truth(ca.vertcat(2.0, 3)));
    Check.Throws<Exception>(() => truth(ca.DM()));
  }

  static List<List<long>> permutations(List<long> b) {
    var r = new List<List<long>>();
    if (b.Count == 0) {
      r.Add(new List<long>());
      return r;
    }
    for (int i = 0; i < b.Count; i++) {
      var rest = new List<long>(b);
      long h = rest[i];
      rest.RemoveAt(i);
      foreach (var p in permutations(rest)) {
        p.Insert(0, h);
        r.Add(p);
      }
    }
    return r;
  }

  static void checkListslice(DM d, long[] rowbase, long[] colbase) {
    foreach (var col in permutations(colbase.ToList()))
      foreach (var rw in permutations(rowbase.ToList())) {
        DM r = DM.zeros(rw.Count, col.Count);
        for (int i = 0; i < rw.Count; i++)
          for (int j = 0; j < col.Count; j++) r[i, j] = d[rw[i], col[j]];
        checkarray(d.get(idx(rw.ToArray()), idx(col.ToArray())), r, 1e-10,
                   d + "[" + string.Join(",", rw) + "," + string.Join(",", col) + "]");
      }
  }

  public void test_listslice() {
    checkListslice(range(Sparsity.dense(3, 3)), new long[] { 0, 1, 2 }, new long[] { 0, 1, 2 });
    checkListslice(range(Sparsity.dense(4, 4)), new long[] { 0, 1, 3 }, new long[] { 0, 2, 3 });
    checkListslice(range(Sparsity.dense(3, 3)), new long[] { 0, 0, 1 }, new long[] { 0, 0, 1 });
    checkListslice(range(Sparsity.dense(3, 3)), new long[] { 0, 0, 2 }, new long[] { 0, 0, 2 });
    checkListslice(range(Sparsity.dense(3, 3)), new long[] { 1, 1, 2 }, new long[] { 1, 1, 2 });
    DM d = range(Sparsity.lower(4));
    checkListslice(d, new long[] { 0, 1, 3 }, new long[] { 0, 2, 3 });
    checkListslice(d.T(), new long[] { 0, 1, 3 }, new long[] { 0, 2, 3 });
    d = range(Sparsity.rowcol(new long[] { 0, 1, 2 }, new long[] { 0, 1 }, 4, 4));
    checkListslice(d, new long[] { 0, 3 }, new long[] { 0, 2 });
    checkListslice(range(Sparsity.dense(2, 2)), new long[] { 0, 0, 0 }, new long[] { 0, 0, 0 });
    checkListslice(range(Sparsity.dense(2, 2)), new long[] { 0, 0, 1 }, new long[] { 0, 0, 1 });
    checkListslice(range(Sparsity.dense(2, 2)), new long[] { 1, 1, 0 }, new long[] { 1, 1, 0 });
    checkListslice(range(Sparsity.dense(2, 2)), new long[] { 1, 1, 1 }, new long[] { 1, 1, 1 });
    d = range(Sparsity.lower(3));
    checkListslice(d, new long[] { 0, 1, 2 }, new long[] { 0, 1, 2 });
    checkListslice(d.T(), new long[] { 0, 1, 2 }, new long[] { 0, 1, 2 });
    d = range(Sparsity.rowcol(new long[] { 0, 2 }, new long[] { 0, 1 }, 4, 4));
    checkListslice(d, new long[] { 0, 1, 3 }, new long[] { 0, 2, 3 });
  }

  // test_sparsesym: returns immediately in Python (feature removed in 73f407e)

  public void test_diagcat() {
    message("diagcat");
    DM C = ca.diagcat(ca.DM(new double[,] { { -1.4, -3.2 }, { -3.2, -28 } }),
                   ca.DM(new double[,] { { 15, -12, 2.1 }, { -12, 16, -3.8 }, { 2.1, -3.8, 15 } }), ca.DM(1.8), ca.DM(-4.0));
    DM r = ca.sparsify(ca.DM(new double[,] { { -1.4, -3.2, 0, 0, 0, 0, 0 }, { -3.2, -28, 0, 0, 0, 0, 0 },
                                       { 0, 0, 15, -12, 2.1, 0, 0 }, { 0, 0, -12, 16, -3.8, 0, 0 },
                                       { 0, 0, 2.1, -3.8, 15, 0, 0 }, { 0, 0, 0, 0, 0, 1.8, 0 },
                                       { 0, 0, 0, 0, 0, 0, -4 } }));
    checkarray(C, r);
  }

  public void test_diag_sparse() {
    message("diag sparse");
    foreach (double[] n in new[] { new double[] { 0, 1, 0, 0, 2, 3, 4, 5, 6, 0 }, new double[] { 1, 2, 3, 0 },
                                   new double[] { 0, 1, 2, 3 } }) {
      DM d = ca.sparsify(n);
      DM D = n;
      DM m = ca.diag(d);
      DM M = ca.sparsify(ca.diag(D));
      Check.True(M.sparsity().colind().SequenceEqual(m.sparsity().colind()), "colind");
      Check.True(M.sparsity().row().SequenceEqual(m.sparsity().row()), "row");
    }
  }

  public void test_sprank() {
    message("sprank");
    assertEquals(3L, ca.sprank(ca.sparsify(ca.DM(new double[,] { { 1, 0, 0 }, { 0, 1, 0 }, { 0, 0, 1 } }))));
    assertEquals(2L, ca.sprank(ca.sparsify(ca.DM(new double[,] { { 1, 0, 0 }, { 0, 0, 0 }, { 0, 0, 1 } }))));
    assertEquals(1L, ca.sprank(ca.sparsify(ca.DM(new double[,] { { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 1 } }))));
    assertEquals(0L, ca.sprank(ca.sparsify(ca.DM(new double[,] { { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 } }))));
    assertEquals(1L, ca.sprank(DM.ones(1, 3)));
    assertEquals(1L, ca.sprank(DM.ones(3, 1)));
    assertEquals(2L, ca.sprank(DM.ones(2, 3)));
    assertEquals(2L, ca.sprank(DM.ones(3, 2)));
    assertEquals(3L, ca.sprank(DM.ones(3, 3)));
    assertEquals(3L, ca.sprank(DM.ones(3, 3)));
    var A = ca.DM(6, 4);
    A[0, 0] = 1;
    A[1, 2] = 1;
    A[2, 2] = 1;
    A[5, 3] = 1;
    assertEquals(3L, ca.sprank(A));
  }

  public void test_cross() {
    message("cross products");
    checkarray(ca.cross(ca.vertcat(1, 0, 0), ca.vertcat(0, 1, 0)), ca.vertcat(0, 0, 1));
    checkarray(ca.cross(ca.vertcat(1.1, 1.3, 1.7), ca.vertcat(2, 3, 13)), ca.vertcat(11.8, -10.9, 0.7));
    checkarray(ca.cross(ca.vertcat(1.1, 1.3, 1.7).T(), ca.vertcat(2, 3, 13).T()), ca.vertcat(11.8, -10.9, 0.7).T());
    DM a = ca.DM(new double[,] { { 1.1, 1.3, 1.7 }, { 1, 0, 0 }, { 0, 0, 1 }, { 4, 5, 6 } });
    DM b = ca.DM(new double[,] { { 2, 3, 13 }, { 0, 1, 0 }, { 0, 0, 1 }, { 1, 0, 1 } });
    DM r = ca.DM(new double[,] { { 11.8, -10.9, 0.7 }, { 0, 0, 1 }, { 0, 0, 0 }, { 5, 2, -5 } });
    checkarray(ca.cross(a, b), r);
    checkarray(ca.cross(a.T(), b.T()), r.T());
    checkarray(ca.cross(a, b, 2), r);
    checkarray(ca.cross(a.T(), b.T(), 1), r.T());
  }

  public void test_is_regular() {
    assertTrue(ca.vertcat(1, 2).is_regular());
    assertFalse(ca.vertcat(1, double.PositiveInfinity).is_regular());
    assertFalse(DM.nan(2).is_regular());
  }

  public void test_sizes() {
    assertEquals(10L, Sparsity.diag(10).nnz_diag());
    assertEquals(10L, Sparsity.diag(10).nnz_upper());
    assertEquals(10L, Sparsity.diag(10).nnz_lower());
    assertEquals(10L * 11 / 2, Sparsity.dense(10, 10).nnz_lower());
    assertEquals(10L * 11 / 2, Sparsity.dense(10, 10).nnz_upper());
    assertEquals(10L, Sparsity.dense(10, 10).nnz_diag());
    DM a = ca.sparsify(ca.DM(new double[,] { { 1, 1, 0 }, { 1, 0, 1 }, { 0, 0, 0 } }));
    assertEquals(1L, a.nnz_diag());
    assertEquals(2L, a.nnz_lower());
    assertEquals(3L, a.nnz_upper());
  }

  public void test_tril2symm() {
    DM a = range(Sparsity.upper(3)).T();
    DM s = ca.tril2symm(a);
    checkarray(s, ca.DM(new double[,] { { 0, 1, 3 }, { 1, 2, 4 }, { 3, 4, 5 } }));
    Check.Throws<Exception>(() => ca.tril2symm(DM.ones(5, 3)));
    Console.WriteLine(DM.ones(5, 5).nnz_upper() - DM.ones(5, 5).nnz_diag());
    Check.Throws<Exception>(() => ca.tril2symm(DM.ones(5, 5)));
  }

  public void test_not_null() {
    MX x = MX.sym("x", 3, 1);
    Sparsity sp = Sparsity.upper(2);
    ca.MX(sp, x);
  }

  public void test_segfault() {
    MX x = MX.sym("x", 10, 1);
    Sparsity sp = Sparsity.upper(2);
    MX y = ca.triu2symm(ca.MX(sp, x[1..4]));
    ca.Function("f", new[] { x }, new[] { y });
  }

  public void test_append_empty() {
    DM a = ca.vertcat(ca.DM(0, 0), ca.DM(0, 2));
    assertEquals(0L, a.size1());
    assertEquals(2L, a.size2());
    a = ca.vertcat(ca.DM(0, 0), ca.DM(2, 0), ca.DM(3, 0));
    assertEquals(5L, a.size1());
    assertEquals(0L, a.size2());
  }

  public void test_vertcat_empty() {
    var a = ca.DM(0, 2);
    DM v = ca.vertcat(a, a);
    assertEquals(0L, v.size1());
    assertEquals(2L, v.size2());
    a = ca.DM(2, 0);
    v = ca.vertcat(a, a);
    assertEquals(4L, v.size1());
    assertEquals(0L, v.size2());
  }

  public void test_vertsplit() {
    DM a = range(Sparsity.upper(5)).T();
    List<DM> v = ca.vertsplit(a, new long[] { 0, 2, 4, 5 });
    assertEquals(3, v.Count);
    checkarray(v[0], ca.DM(new double[,] { { 0, 0, 0, 0, 0 }, { 1, 2, 0, 0, 0 } }));
    checkarray(v[1], ca.DM(new double[,] { { 3, 4, 5, 0, 0 }, { 6, 7, 8, 9, 0 } }));
    checkarray(v[2], ca.DM(new double[,] { { 10, 11, 12, 13, 14 } }));

    v = ca.vertsplit(a);
    assertEquals(a.size1(), (long)v.Count);
    checkarray(v[0], ca.DM(new double[,] { { 0, 0, 0, 0, 0 } }));
    checkarray(v[1], ca.DM(new double[,] { { 1, 2, 0, 0, 0 } }));
    checkarray(v[2], ca.DM(new double[,] { { 3, 4, 5, 0, 0 } }));
    checkarray(v[3], ca.DM(new double[,] { { 6, 7, 8, 9, 0 } }));
    checkarray(v[4], ca.DM(new double[,] { { 10, 11, 12, 13, 14 } }));

    v = ca.vertsplit(a, new long[] { 0, 2, 4, 5 });
    assertEquals(3, v.Count);
    checkarray(v[0], ca.DM(new double[,] { { 0, 0, 0, 0, 0 }, { 1, 2, 0, 0, 0 } }));
    checkarray(v[1], ca.DM(new double[,] { { 3, 4, 5, 0, 0 }, { 6, 7, 8, 9, 0 } }));
    checkarray(v[2], ca.DM(new double[,] { { 10, 11, 12, 13, 14 } }));

    v = ca.vertsplit(a, new long[] { 0, 0, 3, 5 });
    assertEquals(3, v.Count);
    assertEquals(0L, v[0].size1());
    assertEquals(5L, v[0].size2());
    checkarray(v[1], ca.DM(new double[,] { { 0, 0, 0, 0, 0 }, { 1, 2, 0, 0, 0 }, { 3, 4, 5, 0, 0 } }));
    checkarray(v[2], ca.DM(new double[,] { { 6, 7, 8, 9, 0 }, { 10, 11, 12, 13, 14 } }));
  }

  public void test_horzsplit() {
    DM a = range(Sparsity.upper(5)).T();
    List<DM> v = ca.horzsplit(a, new long[] { 0, 2, 4, 5 });
    assertEquals(3, v.Count);
    checkarray(v[0], ca.DM(new double[,] { { 0, 0 }, { 1, 2 }, { 3, 4 }, { 6, 7 }, { 10, 11 } }));
    checkarray(v[1], ca.DM(new double[,] { { 0, 0 }, { 0, 0 }, { 5, 0 }, { 8, 9 }, { 12, 13 } }));
    checkarray(v[2], ca.DM(new double[,] { { 0 }, { 0 }, { 0 }, { 0 }, { 14 } }));

    v = ca.horzsplit(a);
    assertEquals(a.size1(), (long)v.Count);
    checkarray(v[0], ca.vertcat(0, 1, 3, 6, 10));
    checkarray(v[1], ca.vertcat(0, 2, 4, 7, 11));
    checkarray(v[2], ca.vertcat(0, 0, 5, 8, 12));
    checkarray(v[3], ca.vertcat(0, 0, 0, 9, 13));
    checkarray(v[4], ca.vertcat(0, 0, 0, 0, 14));

    v = ca.horzsplit(a, new long[] { 0, 2, 4, 5 });
    assertEquals(3, v.Count);
    checkarray(v[0], ca.DM(new double[,] { { 0, 0 }, { 1, 2 }, { 3, 4 }, { 6, 7 }, { 10, 11 } }));
    checkarray(v[1], ca.DM(new double[,] { { 0, 0 }, { 0, 0 }, { 5, 0 }, { 8, 9 }, { 12, 13 } }));
    checkarray(v[2], ca.DM(new double[,] { { 0 }, { 0 }, { 0 }, { 0 }, { 14 } }));

    v = ca.horzsplit(a, new long[] { 0, 0, 3, 5 });
    assertEquals(3, v.Count);
    assertEquals(5L, v[0].size1());
    assertEquals(0L, v[0].size2());
    checkarray(v[1], ca.DM(new double[,] { { 0, 0, 0 }, { 1, 2, 0 }, { 3, 4, 5 }, { 6, 7, 8 }, { 10, 11, 12 } }));
    checkarray(v[2], ca.DM(new double[,] { { 0, 0 }, { 0, 0 }, { 0, 0 }, { 9, 0 }, { 13, 14 } }));
  }

  public void test_blocksplit() {
    DM a = range(Sparsity.upper(5)).T();
    List<List<DM>> v = ca.blocksplit(a, new long[] { 0, 2, 4, 5 }, new long[] { 0, 1, 3, 5 });
    checkarray(v[0][0], ca.vertcat(0, 1));
    checkarray(v[0][1], ca.DM(new double[,] { { 0, 0 }, { 2, 0 } }));
    checkarray(v[1][0], ca.vertcat(3, 6));
    checkarray(ca.blockcat(v), a);
  }

  public void test_solve() {
    var spA = new List<Sparsity> { Sparsity.dense(1, 1) };
    for (long n = 2; n < 5; n++) {
      spA.AddRange(new[] {
        Sparsity.diag(n), Sparsity.dense(n, n), Sparsity.lower(n), Sparsity.lower(n).T(), Sparsity.banded(n, 1),
        ca.diagcat(new[] { Sparsity.diag(n), Sparsity.dense(n, n) }), ca.diagcat(new[] { Sparsity.diag(n), Sparsity.lower(n) }),
        ca.diagcat(new[] { Sparsity.diag(n), Sparsity.lower(n).T() }),
        ca.diagcat(new[] { Sparsity.lower(n), Sparsity.lower(n).T() }),
        Sparsity.diag(n).plus(Sparsity.rowcol(new long[] { 0 }, new long[] { n - 1 }, n, n)),
        Sparsity.diag(n).plus(Sparsity.rowcol(new long[] { 0, n - 1 }, new long[] { n - 1, 0 }, n, n)),
        Sparsity.diag(n).plus(Sparsity.triplet(n, n, new long[] { 0 }, new long[] { n - 1 })),
        Sparsity.diag(n).plus(Sparsity.triplet(n, n, new long[] { 0, n - 1 }, new long[] { n - 1, 0 })) });
    }
    foreach (Sparsity sA in spA) {
      var a = ca.DM(sA, randomValues(sA.nnz()));
      SX A = SX.sym("a", a.sparsity());
      long n1 = a.size1();
      foreach (Sparsity sB in new[] { Sparsity.dense(n1, 1), ca.vertcat(new[] { Sparsity.dense(1, 1), ca.Sparsity(n1 - 1, 1) }),
                                      Sparsity.lower(n1), Sparsity.lower(n1).T() }) {
        var b = ca.DM(sB, randomValues(sB.nnz()));
        SX B = SX.sym("B", b.sparsity());
        SX C = ca.solve(A, B);
        var f = ca.Function("f", new[] { A, B }, new[] { C });
        DM c = f.Call(a, b);
        DM c_ref = ca.sparsify(ca.vertcat(npSolve(a.ToArray(), b.ToArray())));
        Console.WriteLine(sA.dim() + " " + sB.dim());
        Console.WriteLine("foo" + c + c_ref);
        checkarray(c, c_ref);
      }
    }
  }

  double[] randomValues(long n) {
    var v = new double[n];
    for (int k = 0; k < n; k++) v[k] = random.NextDouble();
    return v;
  }

  public void test_kron() {
    DM a = ca.sparsify(ca.DM(new double[,] { { 1, 0, 6 }, { 2, 7, 0 } }));
    DM b = ca.sparsify(ca.DM(new double[,] { { 1, 0, 0 }, { 2, 3, 7 }, { 0, 0, 9 }, { 1, 12, 13 } }));
    DM c_ = ca.kron(a, b);
    assertEquals(a.size1() * b.size1(), c_.size1());
    assertEquals(a.size2() * b.size2(), c_.size2());
    assertEquals(a.nnz() * b.nnz(), c_.nnz());
    checkarray(c_, npKron(a.ToArray(), b.ToArray()));
  }

  static Function jacobian_old(Function f, long i, long j) {
    var outs = new List<string> { "jac:" + f.name_out(j) + ":" + f.name_in(i) };
    outs.AddRange(f.name_out());
    return f.factory(f.name() + "_jac", f.name_in(), outs);
  }

  public void test_vec_kron() {
    SX A = SX.sym("A", 2, 3);
    SX B = SX.sym("B", 4, 5);
    SX P = SX.sym("P", A.size2(), B.size1());
    var f = ca.Function("f", new[] { ca.vec(P.T()), A, B }, new[] { ca.vec(ca.mtimes(new[] { A, P, B }).T()) });
    Function J = jacobian_old(f, 0, 0);
    var J_in = new List<DM> { ca.vertcat(npRandom((int)ca.vec(P.T()).size1(), 1)), ca.vertcat(npRandom(2, 3)), ca.vertcat(npRandom(4, 5)) };
    DM res = J.call(J_in)[0];
    DM reference = ca.kron(J_in[1], J_in[2].T());
    checkarray(res, reference);
  }

  public void test_repmat() {
    DM a = ca.DM(new double[,] { { 1, 2 }, { 3, 4 }, { 5, 6 } });
    checkarray(ca.repmat(a, 2, 3), ca.kron(DM.ones(2, 3), a));
  }

  public void test_triu() {
    DM a = ca.DM(new double[,] { { 1, 2 }, { 3, 4 } });
    checkarray(ca.triu(a), ca.DM(new double[,] { { 1, 2 }, { 0, 4 } }));
  }

  public void test_tril() {
    DM a = ca.DM(new double[,] { { 1, 2 }, { 3, 4 } });
    checkarray(ca.tril(a), ca.DM(new double[,] { { 1, 0 }, { 3, 4 } }));
  }

  public void test_nz() {
    DM a = ca.sparsify(ca.DM(new double[,] { { 1, 2 }, { 0, 0 }, { 3, 4 } }));
    checkarray(a.get_nz(false, ca.Slice()), ca.vertcat(1, 3, 2, 4));
    assertEquals(4L, a.nnz());
    checkarray(a.get_nz(false, ca.Slice(0, -1)), ca.vertcat(1, 3, 2));
    checkarray(a.get_nz(false, ca.Slice(0)), ca.vertcat(1));
  }

  // get_nz/set_nz without the ind1 flag, like get/set
  public void test_nz_without_ind1() {
    DM a = ca.sparsify(ca.DM(new double[,] { { 1, 2 }, { 0, 0 }, { 3, 4 } }));
    checkarray(a.get_nz(ca.Slice(0, -1)), ca.vertcat(1, 3, 2));
    checkarray(a.get_nz(ca.vertcat(3, 0)), ca.vertcat(4, 1));
    checkarray(a.get_nz(ca.Slice(2)), ca.vertcat(2));
    a.set_nz(0..2, ca.vertcat(7, 8));
    a.set_nz(ca.Slice(3), 9);
    a.set_nz(ca.vertcat(2), 6);
    checkarray(a.nonzeros().ToArray(), ca.vertcat(7, 8, 6, 9));
    MX x = MX.sym("x", 4);
    var f = ca.Function("f", new[] { x }, new[] { x.get_nz(1..3) });
    checkarray(f.call(new[] { ca.vertcat(1, 2, 3, 4) })[0], ca.vertcat(2, 3));
    SX y = SX.sym("y", 3);
    y.set_nz(ca.Slice(1), 0);
    assertEquals("[y_0, 0, y_2]", y.ToString());
  }

  void zeroRandomEntries(double[,] A) {
    for (int i = 0; i < 5; i++) A[random.Next(A.GetLength(0)), random.Next(A.GetLength(1))] = 0;
  }

  public void test_norm_inf_mul() {
    foreach (var (s0, s1, s2) in new[] { (10, 2, 8), (8, 2, 10) }) {
      double[,] A = npRandom(s0, s1), B = npRandom(s1, s2);
      checkarray(ca.norm_inf_mul(ca.DM(A), ca.DM(B)), ca.norm_inf(ca.mtimes(ca.DM(A), ca.DM(B))));
      checkarray(ca.norm_0_mul(ca.DM(A), ca.DM(B)), ca.mtimes(ca.DM(A), ca.DM(B)).nnz());
      zeroRandomEntries(A);
      zeroRandomEntries(B);
      DM As = ca.sparsify(ca.DM(A)), Bs = ca.sparsify(ca.DM(B));
      checkarray(ca.norm_inf_mul(As, Bs), ca.norm_inf(ca.mtimes(As, Bs)));
      checkarray(ca.norm_0_mul(As, Bs), ca.mtimes(As, Bs).nnz());
    }
  }

  public void test_mul3_issue_1465() {
    SX w2 = SX.sym("w", 2, 1);
    SX w3 = SX.sym("w", 3, 1);
    DM Q = DM.eye(2);
    ca.mtimes(w2.T(), Q);
    assertInException("Matrix product with incompatible dimensions. Lhs is 1x3 and rhs is 2x2.",
                      () => ca.mtimes(w3.T(), Q));
    // Java's instance forms w.T().mtimes(Q) are the same call as mtimes(w.T(), Q) here
    var Qs = ca.SX(Q);
    ca.mtimes(new[] { w2.T(), Qs, w2 });
    assertInException("Matrix product with incompatible dimensions. Lhs is 1x3 and rhs is 2x2.",
                      () => ca.mtimes(new[] { w3.T(), Qs, w2 }));
    assertInException("Matrix product with incompatible dimensions. Lhs is 1x2 and rhs is 3x1.",
                      () => ca.mtimes(new[] { w2.T(), Qs, w3 }));
    assertInException("Matrix product with incompatible dimensions. Lhs is 1x3 and rhs is 2x2.",
                      () => ca.mtimes(new[] { w3.T(), Qs, w3 }));
    // python-only: @ None and mtimes(x, y, non-string) do not compile
  }

  public void test_chol() {
    for (int i = 0; i < 4; i++) {
      DM A = ca.vertcat(npRandom(3, 3));
      DM H = ca.mtimes(A, A.T());
      DM R = ca.chol(H);
      assertTrue(R.is_triu());
      checkarray(ca.mtimes(R.T(), R), H);
    }
  }

  public void test_skew() {
    DM x = ca.vertcat(1, 7, 13);
    checkarray(ca.inv_skew(ca.skew(x)), x);
    DM y = ca.vertcat(0.2, 0.9, 0.4);
    checkarray(ca.mtimes(ca.skew(x), y), ca.cross(x, y));
  }

  public void test_nz_overflow() {
    DM d = ca.vertcat(2, 3);
    DM r = d.get_nz(false, ca.Slice());
    checkarray(r, d);
  }

  public void test_DMcrash() {
    // Python DM([DM, DM]) builds a DM from a list of matrices; C# has no such constructor: vertcat is the equivalent
    var a = ca.vertcat(ca.DM(1), ca.DM(2));
    checkarray(a, ca.vertcat(1, 2));
  }


  public void test_sparsity_operation() {
    var L = new List<DM> { ca.DM(1), ca.DM(ca.Sparsity(1, 1), 1), ca.DM(ca.Sparsity(2, 1), 1),
                           ca.DM(Sparsity.dense(2, 1), 1) };
    foreach (DM a in L)
      foreach (DM b in L) {
        DM c = a * b;
        if (a.nnz() == 0 || b.nnz() == 0) assertTrue(c.nnz() == 0);
        else assertTrue(c.nnz() > 0);
      }
    assertTrue(ca.sum2(ca.DM(ca.Sparsity(1, 1), 1)).nnz() == 0);
  }

  static double[,] npCumsum(double[,] A, int axis) {
    int n = A.GetLength(0), m = A.GetLength(1);
    var r = new double[n, m];
    for (int i = 0; i < n; i++)
      for (int j = 0; j < m; j++)
        r[i, j] = A[i, j] + (axis == 0 ? (i > 0 ? r[i - 1, j] : 0) : (j > 0 ? r[i, j - 1] : 0));
    return r;
  }

  public void test_matlab_operations() {
    var data = new[] { new double[,] { { 1, 3 }, { 11, 17 } }, new double[,] { { 1, 3 } }, new double[,] { { 1 }, { 3 } },
                       new double[,] { { 3 } } };
    foreach (double[,] A in data) {
      int nr = A.GetLength(0), nc = A.GetLength(1);
      DM B = ca.reshape(ca.DM(A), (nr, nc));
      checkarray(ca.cumsum(B, 0), npCumsum(A, 0));
      checkarray(ca.cumsum(B, 1), npCumsum(A, 1));
      MX Bs = MX.sym("B", B.size());
      DM A0 = ca.cumsum(B, 0);
      DM A1 = ca.cumsum(B, 1);
      checkarray(ca.evalf(ca.cumsum(ca.MX(B), 0)), npCumsum(A, 0));
      checkarray(ca.evalf(ca.cumsum(B, 1)), npCumsum(A, 1));
      var f = ca.Function("f", new[] { Bs }, new[] { ca.MX(A0), ca.MX(A1) });
      checkfunction(f, f.expand(), new List<DM> { B });
      checkCodegen(f, new List<DM> { B }, true);
      var d = new double[nr - 1, nc];
      for (int i = 0; i + 1 < nr; i++)
        for (int j = 0; j < nc; j++) d[i, j] = A[i + 1, j] - A[i, j];
      DM D = ca.diff(B, 1, 0);
      assertEquals((long)(nr - 1), D.size1());
      assertEquals((long)nc, D.size2());
      if (nr - 1 > 0) checkarray(D, d);
    }
  }

  public void test_singular_repmat() {
    foreach (int n_b in new[] { 0, 2 })
      foreach (int m_b in new[] { 0, 2 })
        foreach (int n in new[] { 0, 3 })
          foreach (int m in new[] { 0, 3 }) {
            var reference = ((long)n_b * n, (long)m_b * m);
            assertEquals(reference, ca.repmat(ca.DM(n_b, m_b), n, m).size());
            assertEquals(reference, ca.repmat(ca.SX(n_b, m_b), n, m).size());
            assertEquals(reference, ca.repmat(ca.MX(n_b, m_b), n, m).size());
            assertEquals(reference, ca.repmat(ca.Sparsity(n_b, m_b), n, m).size());
          }
  }

  public void test_serialize() {
    foreach (DM a in new[] { ca.DM(), ca.DM(2), ca.vertcat(1, 2), ca.DM(new double[,] { { 1, 2 }, { 3, 4 }, { 5, 6 } }) }) {
      DM b = DM.deserialize(a.serialize());
      checkarray(a, b);
    }
  }

  // test_iterable: python-only (Python iteration protocol)

  public void test_toarray() {
    // C# analog of toarray(): ToArray() / elements(); numpy layout flags are python-only
    DM A = ca.DM(new double[,] { { 1, 2, 0 }, { 0, 4, 5 } });
    foreach (DM a in new[] { A, ca.sparsify(A), ca.tril(ca.DM(new double[,] { { 1, 2 }, { 3, 4 } })), ca.vertcat(1, 2, 3), ca.vertcat(1, 2, 3).T(),
                             ca.DM(3), ca.DM(3, 0), ca.DM(0, 3) }) {
      double[,] r = a.ToArray();
      assertEquals(a.size1(), (long)r.GetLength(0));
      assertEquals(a.size2(), (long)r.GetLength(1));
      for (int i = 0; i < r.GetLength(0); i++)
        for (int j = 0; j < r.GetLength(1); j++) assertEquals((double)a[i, j], r[i, j]);
    }
    Check.True(new double[] { 1, 2, 3 }.SequenceEqual(ca.vertcat(1, 2, 3).elements()), "elements");
    assertEquals(3.0, ca.DM(3).to_double());
  }

  public void test_from_file() {
    string dir = Directory.CreateTempSubdirectory("casadi_from_file").FullName;
    string mtx = Path.Combine(dir, "test.mtx"), txt = Path.Combine(dir, "test.txt");
    string txt2 = Path.Combine(dir, "test2.txt");
    DM a = ca.sparsify(ca.DM(new double[,] { { 1, 0, -6, 5, 0 }, { 4, 0, -4e-301, 9.3e-18, 0 } }));
    a.print_dense();
    a[0, 4] = 0;
    a[1, 1] = 0;
    a.to_file(mtx);
    assertInException("foobar.mtx", () => DM.from_file("foobar.mtx"));
    DM b = DM.from_file(mtx);
    assertTrue(a.sparsity().is_equal(b.sparsity()));
    checkarray(a, b);
    a.to_file(txt);
    b = DM.from_file(txt);
    checkarray(a, b);
    foreach (string s in new[] { "1 00 -6 5 0\n4 0 -4e-301 9.3e-18 00\n", "1 00 -6 5 0\n4 0 -4e-301 9.3e-18 00",
                                 "1 00 -6 5 0 \n4 0 -4e-301 9.3e-18 00   ", "1\t00\t-6\t5\t0\t\n4 0 -4e-301 9.3e-18 00\n",
                                 "1 00 -6 5 0\n%foobar\n4 0 -4e-301 9.3e-18 00\n" }) {
      File.WriteAllText(txt2, s);
      b = DM.from_file(txt2);
      assertTrue(a.sparsity().is_equal(b.sparsity()));
      assertEquals(0.0, (double)ca.norm_1(a - b));
    }
    File.WriteAllText(txt, "1 aa 6\n4 8 9e-18\n");
    assertInException("Parsing", () => DM.from_file(txt));
    a = ca.vertcat(3, double.PositiveInfinity, double.NegativeInfinity, double.NaN);
    a.to_file(txt);
    b = DM.from_file(txt);
    checkarray(a, b);
  }

  public void test_norm_2() {
    double[] a = { 1, 2, 3 };
    double r = Math.Sqrt(14);
    checkarray(ca.evalf(ca.norm_2((DM)a)), r);
    checkarray(ca.evalf(ca.norm_2(((DM)a).T())), r);
    checkarray(ca.evalf(ca.norm_2((MX)a)), r);
    checkarray(ca.evalf(ca.norm_2(((MX)a).T())), r);
    // ca.SX(a) does not compile: ambiguous between SX(IEnumerable<double>) and SX(DM)
    checkarray(ca.evalf(ca.norm_2((SX)a)), r);
    checkarray(ca.evalf(ca.norm_2(((SX)a).T())), r);
  }

  public void test_ldl() {
    DM.rng(1);
    DM H = ca.diagcat(DM.rand(5, 5), DM.rand(5, 5));
    H = H + H.T() + DM.eye(10) * 2;
    DM p0 = idx(permutation(10));
    H = H.get(p0, p0);
    var (D, Lt, p) = ca.ldl(H);
    DM P = DM.eye(10).get(ca.Slice(), idx(p.ToArray()));
    DM F = ca.mtimes(new[] { ca.sqrt(ca.diag(D)), DM.eye(10) + Lt, P.T() });
    Console.WriteLine(H - ca.mtimes(F.T(), F));
    assertTrue((double)ca.norm_fro(H - ca.mtimes(F.T(), F)) <= 1e-14);
  }

  public void test_im_bugs() {
    object a = ca.vertcat(1, 2);
    assertTrue(a is DM);
    checkarray(ca.linspace(1, 3, 10), ca.linspace(1.0, 3.0, 10));
  }

  static double[] npLinspace(double a, double b, int n) {
    var r = new double[n];
    double step = (b - a) / (n - 1);
    for (int i = 0; i < n; i++) r[i] = step == 0 ? i / (double)(n - 1) * (b - a) + a : i * step + a;
    if (n > 1) r[n - 1] = b;
    return r;
  }

  public void test_linspace() {
    double[,] cases = { { 1.0, 3.0 }, { -2.5, 7.25 }, { 0.0, 1.0 }, { 1e10, 1e-10 }, { -1.0, 1.0 } };
    int[] sizes = { 2, 3, 5, 10, 100, 1001 };
    for (int c = 0; c < cases.GetLength(0); c++) {
      double a = cases[c, 0], b = cases[c, 1];
      foreach (int n in sizes) {
        double[] reference = npLinspace(a, b, n);
        string w = " a=" + a + " b=" + b + " n=" + n;
        Check.True(reference.SequenceEqual(ca.linspace(ca.DM(a), ca.DM(b), n).nonzeros()), "DM linspace mismatch" + w);
        MX x = MX.sym("x"), y = MX.sym("y");
        var f = ca.Function("f", new[] { x, y }, new[] { ca.linspace(x, y, n) });
        Check.True(reference.SequenceEqual(f.Call(a, b).nonzeros()), "MX linspace mismatch" + w);
        SX xs = SX.sym("x"), ys = SX.sym("y");
        f = ca.Function("f", new[] { xs, ys }, new[] { ca.linspace(xs, ys, n) });
        Check.True(reference.SequenceEqual(f.Call(a, b).nonzeros()), "SX linspace mismatch" + w);
      }
    }
  }

  public void test_linspace_mx_node_count() {
    MX x = MX.sym("x"), y = MX.sym("y");
    long n_small = ca.n_nodes(ca.linspace(x, y, 5));
    long n_large = ca.n_nodes(ca.linspace(x, y, 1000));
    assertTrue(n_large < n_small + 5, "MX linspace node count grew with nsteps: " + n_small + " -> " + n_large);
  }

  public void test_permutation() {
    int n = 10;
    long[] p = permutation(n);
    Sparsity S = Sparsity.permutation(p);
    Check.True(p.SequenceEqual(S.permutation_vector()), "permutation_vector");
    assertTrue(S.is_permutation());
    DM v = DM.rand(n);
    // Python S @ v converts the Sparsity to DM.ones(S)
    checkarray(ca.mtimes(DM.ones(S), v), v[p]);
    S = Sparsity.permutation(p, true);
    checkarray(ca.mtimes(DM.ones(S), v)[p], v);
  }

  public void test_sparsity_orthonormality() {
    var alltests = new List<(string, Func<Sparsity, bool>)> {
      ("is_orthonormal_rows(true)", s => s.is_orthonormal_rows(true)),
      ("is_orthonormal_columns(true)", s => s.is_orthonormal_columns(true)),
      ("is_selection(true)", s => s.is_selection(true)),
      ("is_orthonormal_rows(false)", s => s.is_orthonormal_rows(false)),
      ("is_orthonormal_columns(false)", s => s.is_orthonormal_columns(false)),
      ("is_selection(false)", s => s.is_selection(false)),
      ("is_permutation()", s => s.is_permutation()),
      ("is_orthonormal(false)", s => s.is_orthonormal(false)),
      ("is_orthonormal(true)", s => s.is_orthonormal(true)),
    };
    void check(Sparsity S, IList<string> prop) {
      foreach (string p in prop) assertTrue(alltests.Any(t => t.Item1 == p));
      foreach (var (name, test) in alltests) {
        bool quality = prop.Contains(name);
        assertEquals(quality, test(S), "Call " + name + " expected to return " + quality);
      }
    }
    Sparsity sp(double[,] m) => ca.sparsify(ca.DM(m)).sparsity();
    var sel = new List<string> { "is_orthonormal_rows(true)", "is_orthonormal_columns(true)", "is_selection(true)",
                                 "is_orthonormal(true)" };
    foreach (double[,] S in new[] {
               new double[,] { { 0, 0, 0, 0, 0 }, { 0, 0, 0, 0, 0 }, { 0, 0, 0, 0, 0 }, { 0, 0, 0, 0, 0 } },
               new double[,] { { 0, 0, 0, 0, 0 }, { 0, 1, 0, 0, 0 }, { 0, 0, 0, 0, 0 }, { 0, 0, 0, 0, 0 } },
               new double[,] { { 0, 0, 0, 0, 0 }, { 0, 1, 0, 0, 0 }, { 0, 0, 0, 1, 0 }, { 0, 0, 0, 0, 0 } } })
      check(sp(S), sel);
    foreach (double[,] S in new[] {
               new double[,] { { 0, 0, 0, 0, 0 }, { 0, 1, 1, 0, 0 }, { 0, 0, 0, 0, 0 }, { 0, 0, 0, 0, 0 } },
               new double[,] { { 0, 1, 0, 0 }, { 1, 0, 0, 0 }, { 0, 1, 1, 0 }, { 0, 0, 0, 0 }, { 0, 0, 0, 1 } },
               new double[,] { { 0, 0, 0, 0, 0 }, { 0, 1, 0, 0, 0 }, { 0, 1, 0, 0, 0 }, { 0, 0, 0, 0, 0 } } })
      check(sp(S), new List<string>());
    check(sp(new double[,] { { 0, 1, 0, 0, 0 }, { 1, 0, 0, 0, 0 }, { 0, 0, 1, 0, 0 }, { 0, 0, 0, 1, 0 } }),
          new List<string> { "is_orthonormal_rows(true)", "is_orthonormal_columns(true)", "is_selection(true)",
                             "is_orthonormal_rows(false)", "is_selection(false)", "is_orthonormal(true)" });
    check(sp(new double[,] { { 0, 1, 0, 0 }, { 1, 0, 0, 0 }, { 0, 0, 1, 0 }, { 0, 0, 0, 1 } }),
          alltests.Select(t => t.Item1).ToList());
    check(sp(new double[,] { { 0, 1, 0, 0 }, { 1, 0, 0, 0 }, { 0, 0, 1, 0 }, { 0, 0, 0, 0 }, { 0, 0, 0, 1 } }),
          new List<string> { "is_orthonormal_rows(true)", "is_orthonormal_columns(true)", "is_selection(true)",
                             "is_orthonormal_columns(false)", "is_orthonormal(true)" });
  }

  static void generateAndDelete(Function f) {
    f.generate("f.c");
    File.Delete("f.c");
  }

  public void test_permutation_solve() {
    int n = 10;
    DM.rng(1);
    long[] p = permutation(n);
    Sparsity P = Sparsity.permutation(p);
    SX xs = SX.sym("x", P), bs = SX.sym("b", n);
    assertTrue(xs.sparsity().is_orthonormal());
    ca.sparsity_cast(xs, Sparsity.dense(n));
    var f0 = ca.Function("f", new[] { xs, bs }, new[] { ca.solve(xs, bs) });
    generateAndDelete(f0);
    MX xm = MX.sym("x", P), bm = MX.sym("b", n);
    assertTrue(xm.sparsity().is_orthonormal());
    ca.sparsity_cast(xm, Sparsity.dense(n));
    var f1 = ca.Function("f", new[] { xm, bm }, new[] { ca.solve(xm, bm) });
    generateAndDelete(f1);
    DM x0 = DM.rand(xm.sparsity());
    DM b0 = DM.rand(bm.sparsity());
    checkarray(f0.Call(x0, b0), f1.Call(x0, b0));
  }

  public void test_permutation_inv() {
    int n = 10;
    DM.rng(1);
    long[] p = permutation(n);
    Sparsity P = Sparsity.permutation(p);
    SX xs = SX.sym("x", P);
    var f0 = ca.Function("f", new[] { xs }, new[] { ca.solve(xs, SX.eye(n)) });
    generateAndDelete(f0);
    MX xm = MX.sym("x", P);
    var f1 = ca.Function("f", new[] { xm }, new[] { ca.solve(xm, MX.eye(n)) });
    generateAndDelete(f1);
    var x0 = ca.DM(P, DM.rand(n));
    checkarray(f0.Call(x0), f1.Call(x0));

    f0 = ca.Function("f", new[] { xs }, new[] { ca.inv(xs) });
    f0.disp(true);
    f1 = ca.Function("f", new[] { xm }, new[] { ca.inv(xm) });
    f1.disp(true);
    x0 = ca.DM(P, DM.rand(n));
    checkarray(f0.Call(x0), f1.Call(x0));
  }

  public void test_mac_simplify() {
    int n = 10;
    DM.rng(1);
    long[] p = permutation(n);
    Sparsity P = Sparsity.permutation(p);
    SX bs = SX.sym("b", n), xs = SX.sym("x", P);
    var f0 = ca.Function("f", new[] { xs, bs }, new[] { ca.mtimes(xs, bs) });
    generateAndDelete(f0);
    MX bm = MX.sym("b", n), xm = MX.sym("x", P);
    var f1 = ca.Function("f", new[] { xm, bm }, new[] { ca.mtimes(xm, bm) });
    generateAndDelete(f1);
    var x0 = ca.DM(P, DM.rand(n));
    DM b0 = DM.rand(bm.sparsity());
    checkarray(f0.Call(x0, b0), f1.Call(x0, b0));
  }

  public void test_logsumexp_empty() {
    assertInException("Argument must be nonempty", () => ca.logsumexp(ca.DM(0, 1)));
  }

  public void test_linspace_nsteps() {
    foreach (long n in new long[] { 0, -1 })
      assertInException("must be at least 1", () => ca.linspace(ca.DM(0), ca.DM(1), n));
  }

  public void test_densify_value() {
    checkarray(ca.densify(ca.DM(Sparsity.diag(2), 1), 5), ca.DM(new double[,] { { 1, 5 }, { 5, 1 } }));
  }

  public void test_polyval_sparse() {
    checkarray(ca.polyval(ca.vertcat(1, 2), ca.DM(Sparsity.diag(2), 1)), ca.DM(new double[,] { { 3, 2 }, { 2, 3 } }));
    MX x = MX.sym("x");
    var f = ca.Function("f", new[] { x }, new[] { ca.polyval(ca.MX(ca.DM(new double[,] { { 1, 2 } })), x) });
    checkarray(f.Call(3), ca.DM(5));
  }

  public void test_det_empty_column() {
    checkarray(ca.det(ca.DM(Sparsity.triplet(2, 2, new long[] { 0, 1 }, new long[] { 0, 0 }), ca.vertcat(1, 2))), ca.DM(0));
  }

  public void test_negative_size_checks() {
    assertInException("must be a multiple of n", () => ca.repsum(ca.vertcat(1, 2, 3), 2));
    assertInException("Dimensions must be non-negative", () => MX.sym("x", -1, 1));
  }

  public void test_inplace() {
    int n = 10;
    MX x = MX.sym("x", n);
    var y = ca.MX(Sparsity.diag(n), x);
    var f = ca.Function("f", new[] { x }, new[] { y });
    assertTrue(f.sz_w() == n);
  }

  public void test_split_n_zero() {
    DM x = DM.ones(2, 4);
    assertInException("must be a multiple of n", () => ca.horzsplit_n(x, 0));
    assertInException("must be a multiple of n", () => ca.vertsplit_n(x.T(), 0));
  }

  public void test_diagsplit_offset_lengths() {
    assertInException("offset1 and offset2 must have the same length", () => ca.diagsplit(DM.ones(4, 4), 1, 2));
  }

  public void test_vertcat_columns() {
    assertInException("Mismatching number of columns", () => ca.vertcat(DM.ones(1, 2), DM.ones(1, 3)));
    assertInException("Mismatching number of columns",
                      () => ca.blockcat(new[] { new[] { ca.DM(1), ca.DM(2) }, new[] { ca.DM(3) } }));
    assertEquals((2L, 2L), ca.blockcat(new[] { new[] { ca.DM(1), ca.DM(2) }, new[] { DM.ones(1, 2) } }).size());
  }

  public void test_split_user_errors() {
    DM x = DM.ones(2, 4);
    foreach (long[] offset in new[] { new long[] { 1, 4 }, new long[] { 0, 3, 2, 4 } })
      assertInException("must start at 0, be nondecreasing and end at 4", () => ca.horzsplit(x, offset));
    foreach (long incr in new long[] { 0, -1 }) assertInException("must be positive", () => ca.horzsplit(x, incr));
    assertInException("must be positive", () => ca.blocksplit(DM.ones(4, 4), 0, 2));
    assertInException("must start at 0, be nondecreasing and end at 4",
                      () => ca.diagsplit(DM.ones(4, 4), new long[] { 0, 1, 3 }, new long[] { 0, 2, 4 }));
    assertInException("Reshape mismatch", () => ca.reshape(DM.eye(2), Sparsity.dense(4, 2)));
    assertInException("sparsity_cast", () => ca.sparsity_cast(DM.eye(2), Sparsity.dense(3, 1)));
  }

  public void test_repmat_negative() {
    assertInException("must be non-negative", () => ca.repmat(ca.vertcat(1, 2), -1));
  }

  public void test_blockcat_empty_diagsplit_offblock() {
    assertEquals((0L, 0L), ca.blockcat(new List<List<DM>>()).size());
    assertEquals((0L, 0L), ca.blockcat(new[] { new List<DM>() }).size());
    assertEquals((1L, 1L), ca.blockcat(new[] { new List<DM>(), new List<DM> { ca.DM(1) } }).size());
    DM X1 = ca.DM(new double[,] { { 1, 2 }, { 3, 4 } });
    assertInException("outside the diagonal blocks", () => ca.diagsplit(X1, 1));
    assertInException("outside the diagonal blocks", () => ca.diagsplit(SX.sym("x", 2, 2), 1));
    assertInException("outside the diagonal blocks", () => ca.diagsplit(Sparsity.dense(2, 2), 1));
    assertInException("outside the diagonal blocks is unsupported", () => ca.diagsplit(MX.sym("x", 2, 2), 1));
  }

  public void test_triplet_empty() {
    assertEquals((0L, 0L), DM.triplet(new long[] { }, new long[] { }, ca.DM(0, 1)).size());
  }

  public void test_triplet_default_size() {
    DM A = DM.triplet(new long[] { 0, 1 }, new long[] { 1, 0 }, ca.vertcat(5, 6));
    assertEquals((2L, 2L), A.size());
    checkarray(A, ca.DM(new double[,] { { 0, 5 }, { 6, 0 } }));
  }

  public void test_eye_negative() {
    assertInException("Dimensions must be non-negative", () => DM.eye(-1));
    assertInException("Dimensions must be non-negative", () => SX.eye(-1));
    assertInException("Dimensions must be non-negative", () => MX.eye(-1));
  }

  public void test_horzsplit_n() {
    assertTrue(ca.horzsplit_n(DM.rand(1, 6), 2).Count == 2);
    assertTrue(ca.vertsplit_n(DM.rand(6), 2).Count == 2);
    assertTrue(ca.horzsplit_n(ca.DM(0, 0), 2).Count == 2);
    assertTrue(ca.vertsplit_n(ca.DM(0, 0), 2).Count == 2);
    Check.Throws<Exception>(() => ca.horzsplit_n(DM.rand(1, 5), 2));
    Check.Throws<Exception>(() => ca.vertsplit_n(DM.rand(5), 2));
  }
}

}

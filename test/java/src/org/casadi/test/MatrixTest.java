package org.casadi.test;

import static org.junit.jupiter.api.Assertions.assertArrayEquals;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Collections;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.function.BiFunction;
import java.util.function.Function;
import java.util.function.Predicate;

import org.casadi.DM;
import org.casadi.GlobalOptions;
import org.casadi.Importer;
import org.casadi.MX;
import org.casadi.SX;
import org.casadi.Slice;
import org.casadi.Sparsity;
import org.casadi.ca;
import org.junit.jupiter.api.Assumptions;
import org.junit.jupiter.api.Test;

/** Port of test/python/matrix.py */
public class MatrixTest extends CasadiTestCase {

  private static final long NONE_START = Long.MIN_VALUE, NONE_STOP = Long.MAX_VALUE;

  private static DM idx(long... i) {
    double[] r = new double[i.length];
    for (int k = 0; k < i.length; k++) r[k] = i[k];
    return new DM(r);
  }

  private static DM range(Sparsity sp) {
    double[] v = new double[(int) sp.nnz()];
    for (int k = 0; k < v.length; k++) v[k] = k;
    return new DM(sp, new DM(v));
  }

  private static DM row(double... x) {
    return new DM(x).T();
  }

  private long[] permutation(int n) {
    List<Long> p = new ArrayList<>();
    for (long i = 0; i < n; i++) p.add(i);
    Collections.shuffle(p, random);
    return toLongs(p);
  }

  private static double npDet(DM m) {
    double[][] a = m.full();
    int n = a.length;
    double d = 1;
    for (int k = 0; k < n; k++) {
      int p = k;
      for (int i = k + 1; i < n; i++) if (Math.abs(a[i][k]) > Math.abs(a[p][k])) p = i;
      if (a[p][k] == 0) return 0;
      if (p != k) {
        double[] t = a[p];
        a[p] = a[k];
        a[k] = t;
        d = -d;
      }
      d *= a[k][k];
      for (int i = k + 1; i < n; i++) {
        double f = a[i][k] / a[k][k];
        for (int j = k; j < n; j++) a[i][j] -= f * a[k][j];
      }
    }
    return d;
  }

  private static double[][] npSolve(double[][] A0, double[][] B0) {
    int n = A0.length, m = B0[0].length;
    double[][] a = new double[n][], b = new double[n][];
    for (int i = 0; i < n; i++) {
      a[i] = A0[i].clone();
      b[i] = B0[i].clone();
    }
    for (int k = 0; k < n; k++) {
      int p = k;
      for (int i = k + 1; i < n; i++) if (Math.abs(a[i][k]) > Math.abs(a[p][k])) p = i;
      double[] t = a[p];
      a[p] = a[k];
      a[k] = t;
      t = b[p];
      b[p] = b[k];
      b[k] = t;
      for (int i = k + 1; i < n; i++) {
        double f = a[i][k] / a[k][k];
        for (int j = k; j < n; j++) a[i][j] -= f * a[k][j];
        for (int j = 0; j < m; j++) b[i][j] -= f * b[k][j];
      }
    }
    double[][] x = new double[n][m];
    for (int i = n - 1; i >= 0; i--)
      for (int j = 0; j < m; j++) {
        double s = b[i][j];
        for (int k = i + 1; k < n; k++) s -= a[i][k] * x[k][j];
        x[i][j] = s / a[i][i];
      }
    return x;
  }

  private static double[][] npKron(double[][] a, double[][] b) {
    int n = a.length, m = a[0].length, p = b.length, q = b[0].length;
    double[][] r = new double[n * p][m * q];
    for (int i = 0; i < n; i++)
      for (int j = 0; j < m; j++)
        for (int k = 0; k < p; k++)
          for (int l = 0; l < q; l++) r[i * p + k][j * q + l] = a[i][j] * b[k][l];
    return r;
  }

  private double[][] npRandom(int n, int m) {
    double[][] r = new double[n][m];
    for (int i = 0; i < n; i++)
      for (int j = 0; j < m; j++) r[i][j] = random.nextDouble();
    return r;
  }

  /** Python helper check_codegen: generate C, compile with the shell importer, compare; optionally serialize */
  private void checkCodegen(org.casadi.Function f, List<DM> inputs, boolean checkSerialize) {
    if (checkSerialize) {
      org.casadi.Function g = org.casadi.Function.deserialize(f.serialize());
      List<DM> a = f.call(inputs), b = g.call(inputs);
      for (int i = 0; i < a.size(); i++) checkarray(a.get(i), b.get(i), "serialize output " + i);
    }
    if (!RUN_SLOW) return;
    try {
      Path dir = Files.createTempDirectory("casadi_cg");
      String name = "cg_" + f.name() + "_" + Long.toHexString(random.nextLong() & Long.MAX_VALUE);
      org.casadi.CodeGenerator cg = new org.casadi.CodeGenerator(name + ".c");
      cg.add(f);
      String src = cg.generate(dir.toString() + "/");
      org.casadi.Function F = ca.external(f.name(), new Importer(src, "shell"));
      List<DM> a = f.call(inputs), b = F.call(inputs);
      for (int i = 0; i < a.size(); i++) checkarray(a.get(i), b.get(i), "codegen output " + i);
    } catch (IOException e) {
      throw new RuntimeException(e);
    }
  }

  @Test
  public void test_constructorlol() {
    message("List of list constructor");
    DM a = new DM(new double[][] {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}});
    DM b = new DM(3, 3);
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++) b.set(i, j, 3 * i + j + 1);
    checkarray(a, b, "List of list constructor");
  }

  @Test
  public void test_sum() {
    message("sum");
    DM D = dm(new double[][] {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}});
    checkarray(ca.sum1(D), new double[][] {{12, 15, 18}});
    checkarray(ca.sum2(D), new double[][] {{6}, {15}, {24}});
  }

  @Test
  public void test_inv() {
    message("Matrix inverse");
    DM a = dm(new double[][] {{1, 2}, {1, 3}});
    checkarray(ca.mtimes(ca.inv(a), a), DM.eye(2), "DM inverse");
  }

  @Test
  public void test_trans() {
    message("trans");
    DM a = new DM(0, 1);
    DM b = a.T();
    assertEquals(1, b.size1());
    assertEquals(0, b.size2());
  }

  @Test
  public void test_vertcat() {
    message("vertcat");
    DM A = DM.ones(2, 3);
    DM B = new DM(4, 3);
    DM C = ca.vertcat(List.of(A, B));
    assertArrayEquals(new long[] {6, 3}, C.size(), "vertcat shape");
    assertEquals(A.nnz(), C.nnz(), "vertcat size");
    assertThrows(RuntimeException.class, () -> ca.horzcat(List.of(A, B)));
  }

  @Test
  public void test_horzcat() {
    message("horcat");
    DM A = DM.ones(3, 2);
    DM B = new DM(3, 4);
    DM C = ca.horzcat(List.of(A, B));
    assertArrayEquals(new long[] {3, 6}, C.size(), "horzcat shape");
    assertEquals(A.nnz(), C.nnz(), "vertcat size");
    assertThrows(RuntimeException.class, () -> ca.vertcat(List.of(A, B)));
  }

  /** First test_diagcat of matrix.py, shadowed there by the later definition */
  @Test
  public void test_diagcat_shadowed() {
    MX x = MX.sym("x", 2, 2);
    MX y = MX.sym("y", Sparsity.lower(3));
    MX z = MX.sym("z", 4, 2);
    List<MX> L = List.of(x, y, z);
    org.casadi.Function fMX = new org.casadi.Function("fMX", L, List.of(ca.diagcat(L)));
    List<SX> LSX = new ArrayList<>();
    for (MX i : L) LSX.add(SX.sym("", i.sparsity()));
    org.casadi.Function fSX = new org.casadi.Function("fSX", LSX, List.of(ca.diagcat(LSX)));
    List<DM> in = new ArrayList<>();
    for (long i = 0; i < 3; i++) in.add(range(fMX.sparsity_in(i)));
    checkfunction(fMX, fSX, in);
  }

  @Test
  public void test_veccat() {
    message("vecccat");
    DM A = new DM(2, 3);
    A.set(0, 1, 2);
    A.set(1, 0, 1);
    A.set(1, 2, 3);
    DM B = new DM(3, 1);
    B.set(0, 0, 4);
    B.set(1, 0, 5);
    B.set(2, 0, 6);
    DM C = ca.veccat(List.of(A, B));
    assertArrayEquals(new long[] {9, 1}, C.size(), "veccat shape");
    assertEquals(A.nnz() + B.nnz(), C.nnz(), "veccat size");
    assertArrayEquals(new double[] {1, 2, 3, 4, 5, 6}, C.nonzeros(), "numbers shape");
  }

  /** Python list slice a1[start:stop:step], with NONE_START/NONE_STOP for omitted bounds */
  private static DM pySlice(double[] a, long start, long stop, long step) {
    int n = a.length;
    long s, e;
    if (step > 0) {
      s = start == NONE_START ? 0 : start < 0 ? Math.max(0, start + n) : Math.min(start, n);
      e = stop == NONE_STOP ? n : stop < 0 ? Math.max(0, stop + n) : Math.min(stop, n);
    } else {
      s = start == NONE_START ? n - 1 : start < 0 ? Math.max(-1, start + n) : Math.min(start, n - 1);
      e = stop == NONE_STOP ? -1 : stop < 0 ? Math.max(-1, stop + n) : Math.min(stop, n - 1);
    }
    List<Double> r = new ArrayList<>();
    for (long i = s; step > 0 ? i < e : i > e; i += step) r.add(a[(int) i]);
    return new DM(r.stream().mapToDouble(Double::doubleValue).toArray());
  }

  @Test
  public void test_slicestepnegative() {
    message("Slice step negative");
    double[] a1 = {1, 2, 3, 4, 5};
    DM a2 = new DM(a1);
    long[][] cases = {{0, 4, -1}, {4, 0, -1}, {0, 4, -2}, {4, 0, -2}, {1, 4, -2}, {4, 1, -2}, {0, 3, -2}, {3, 0, -2},
        {NONE_START, NONE_STOP, -1}, {NONE_START, NONE_STOP, 1}, {2, NONE_STOP, -1}, {NONE_START, 2, -1},
        {-1, NONE_STOP, -1}, {NONE_START, -1, -1}};
    for (long[] c : cases)
      checkarray(a2.get(new Slice(c[0], c[1], c[2]), 0), pySlice(a1, c[0], c[1], c[2]),
          "[" + c[0] + ":" + c[1] + ":" + c[2] + "]");
  }

  @Test
  public void test_indexingOutOfBounds() {
    message("Indexing out of bounds");
    DM y = DM.zeros(4, 5);
    assertThrows(RuntimeException.class, () -> y.get(12, 0));
    assertThrows(RuntimeException.class, () -> y.get(12, 12));
    assertThrows(RuntimeException.class, () -> y.get(0, 12));
    assertThrows(RuntimeException.class, () -> y.get(12, new Slice()));
    assertThrows(RuntimeException.class, () -> y.get(new Slice(12, 15), 0));
    assertThrows(RuntimeException.class, () -> y.get(new Slice(), 12));
    assertThrows(RuntimeException.class, () -> y.get(0, new Slice(12, 15)));
    y.get(-1, 2);
    assertThrows(RuntimeException.class, () -> y.get(-12, 2));
    y.get(new Slice(-3, -1), 2);
    assertThrows(RuntimeException.class, () -> y.get(new Slice(-12, -9), 2));

    assertThrows(RuntimeException.class, () -> y.set(12, 0, 0));
    assertThrows(RuntimeException.class, () -> y.set(12, 12, 0));
    assertThrows(RuntimeException.class, () -> y.set(0, 12, 0));
    assertThrows(RuntimeException.class, () -> y.set(12, new Slice(), 0));
    assertThrows(RuntimeException.class, () -> y.set(new Slice(12, 15), 0, 0));
    assertThrows(RuntimeException.class, () -> y.set(new Slice(), 12, 0));
    assertThrows(RuntimeException.class, () -> y.set(0, new Slice(12, 15), 0));
    y.set(-1, 2, 0);
    assertThrows(RuntimeException.class, () -> y.set(-12, 2, 0));
    y.set(new Slice(-3, -1), 2, 0);
    assertThrows(RuntimeException.class, () -> y.set(new Slice(-12, -9), 2, 0));
  }

  // huge_slice: not a test method in Python (no test_ prefix)

  @Test
  public void test_nonmonotonous_indexing() {
    message("non-monotonous indexing");
    DM A = dm(new double[][] {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}});
    DM B = A.get(idx(0, 2, 1), 0);
    checkarray(dm(1, 7, 4), B, "non-monotonous");
    B = A.get(0, idx(0, 2, 1));
    checkarray(row(1, 3, 2), B, "non-monotonous");
  }

  @Test
  public void test_IM_indexing() {
    message("IM");
    DM A = new DM(2, 2);
    A.set(0, 0, 1);
    A.set(1, 1, 3);
    A.set(0, 1, 2);
    A.set(1, 0, 4);
    DM B = dm(1, 2, 3, 4, 5);
    DM B_ = B.get(A);
    checkarray(B_, dm(new double[][] {{2, 3}, {5, 4}}), "matrix indexing");
    B.set(A, dm(new double[][] {{1, 2}, {3, 4}}));
    checkarray(B, dm(1, 1, 2, 4, 3), "matrix indexing assignement");

    B = dm(1, 2, 3, 4, 5);
    B_ = B.get(A);
    checkarray(new DM(B_.full()), dm(new double[][] {{2, 3}, {5, 4}}), "Imatrix indexing");
    B.set(A, dm(new double[][] {{1, 2}, {3, 4}}));
    checkarray(new DM(B.full()), dm(1, 1, 2, 4, 3), "Imatrix indexing assignement");

    DM B2 = new DM(5, 1);
    assertThrows(RuntimeException.class, () -> B2.get_nz(false, A));
  }

  @Test
  public void test_sparsity_indexing() {
    message("sparsity");
    DM B = dm(new double[][] {{1, 2, 3, 4, 5}, {6, 7, 8, 9, 10}});
    DM A = ca.sparsify(dm(new double[][] {{1, 1, 0, 0, 0}, {0, 0, 1, 0, 0}}));
    Sparsity sp = A.sparsity();
    checkarray(B.get(sp), dm(new double[][] {{1, 2, 0, 0, 0}, {0, 0, 8, 0, 0}}), "sparsity indexing");
    B.set(sp, -4);
    checkarray(B, dm(new double[][] {{-4, -4, 3, 4, 5}, {6, 7, -4, 9, 10}}), "sparsity indexing assignement");
    DM B2 = dm(new double[][] {{1, 2, 3, 4, 5}, {6, 7, 8, 9, 10}});
    B2.set(sp, B2.times(2));
    checkarray(B2, dm(new double[][] {{2, 4, 3, 4, 5}, {6, 7, 16, 9, 10}}), "Imatrix indexing assignement");
    assertThrows(RuntimeException.class, () -> B2.get(Sparsity.dense(4, 4)));
  }

  @Test
  public void test_index_setting() {
    message("index setting");
    DM B = dm(1, 2, 3, 4, 5);
    B.set(0, 8);
    checkarray(B, dm(8, 2, 3, 4, 5), "index setting");
    B.set(1, 0, 4);
    checkarray(B, dm(8, 4, 3, 4, 5), "index setting");
    B.set(new Slice(), 0, 7);
    checkarray(B, dm(7, 7, 7, 7, 7), "index setting");
  }

  @Test
  public void test_issue298() {
    message("Issue #298");
    DM a = new DM(4, 1);
    Object b = ca.reshape(a, 2, 2);
    assertEquals(a.getClass(), b.getClass());
  }

  @Test
  public void test_det() {
    message("Determinant");
    DM a = new DM(1, 1);
    a.set(0, 0, 5);
    checkarray(ca.det(a), npDet(a));

    a = new DM(5, 5);
    for (int i = 0; i < 5; i++) a.set(i, i, i + 1);
    checkarray(ca.det(a), npDet(a));

    a = new DM(5, 5);
    for (int i = 0; i < 4; i++) a.set(i, i, i + 1);
    a.set(0, 4, 3);
    a.set(4, 0, 7);
    checkarray(ca.det(a), npDet(a));

    a = new DM(5, 5);
    for (int i = 0; i < 5; i++)
      for (int j = 0; j < 5; j++) a.set(i, j, i + j);
    checkarray(ca.det(a), npDet(a));

    int[][] shapes = {{4, 5}, {5, 4}};
    for (int[] s : shapes) {
      a = new DM(5, 5);
      for (int i = 0; i < s[0]; i++)
        for (int j = 0; j < s[1]; j++) a.set(i, j, i + j);
      checkarray(ca.det(a), npDet(a));
    }
    long[][] extra = {{4, 1}, {1, 4}, {4, 2}, {2, 4}};
    for (long[] e : extra) {
      a = new DM(5, 5);
      int ni = e[0] == 4 ? 4 : 5, nj = e[0] == 4 ? 5 : 4;
      for (int i = 0; i < ni; i++)
        for (int j = 0; j < nj; j++) a.set(i, j, i + j);
      a.set(e[0], e[1], 12);
      checkarray(ca.det(a), npDet(a));
    }

    a = new DM(50, 50);
    for (int i = 0; i < 50; i++) a.set(i, i, i + 1);
    checkarray(ca.det(a).div(npDet(a)), 1);
  }

  @Test
  public void test_inv_sparsity() {
    Assumptions.assumeTrue(GlobalOptions.getSimplificationOnTheFly());
    message("sparsity pattern of inverse");
    int n = 8;
    Sparsity sp = Sparsity.lower(n);
    List<SX> syms = new ArrayList<>();
    for (int i = 0; i < sp.nnz(); i++) syms.add(SX.sym("a" + i));
    SX x = new SX(sp, ca.vertcat(syms));
    DM x_ = DM.ones(x.sparsity());
    DM I_ = DM.ones(ca.inv(x).sparsity());
    checkarray(x_, I_, "inv");

    sp = Sparsity.lower(n);
    x = SX.sym("a", sp);
    x.set(0, n - 1, 1);
    I_ = DM.ones(ca.inv(x).sparsity());
    checkarray(DM.ones(n, n), I_, "inv");

    x = SX.sym("a", sp);
    x.set(0, n / 2, 1);
    DM s_ = DM.ones(sp);
    s_.set(new Slice(), new Slice(0, n / 2 + 1), 1);
    I_ = DM.ones(ca.inv_minor(x).sparsity());
    s_ = ca.densify(s_);
    ca.densify(I_);
    checkarray(s_, I_, "inv");
  }

  @Test
  public void test_mtimes() {
    DM A = DM.ones(new long[] {4, 3});
    DM B = DM.ones(new long[] {3, 8});
    DM C = DM.ones(new long[] {8, 7});
    assertThrows(RuntimeException.class, () -> ca.mtimes(new ArrayList<DM>()));
    DM D = ca.mtimes(List.of(A));
    assertEquals(4, D.size()[0]);
    assertEquals(3, D.size()[1]);
    D = ca.mtimes(List.of(A, B));
    assertEquals(4, D.size()[0]);
    assertEquals(8, D.size()[1]);
    D = ca.mtimes(List.of(A, B, C));
    assertEquals(4, D.size()[0]);
    assertEquals(7, D.size()[1]);
  }

  @Test
  public void test_remove() {
    message("remove");
    DM B = dm(new double[][] {{1, 2, 3, 4}, {5, 6, 7, 8}, {9, 10, 11, 12}, {13, 14, 15, 16}, {17, 18, 19, 20}});
    DM A = new DM(B);
    A.remove(new long[] {}, new long[] {});
    checkarray(A, B, "remove nothing");
    A = new DM(B);
    A.remove(new long[] {}, new long[] {1});
    checkarray(A, dm(new double[][] {{1, 3, 4}, {5, 7, 8}, {9, 11, 12}, {13, 15, 16}, {17, 19, 20}}),
        "remove a column");
    A = new DM(B);
    A.remove(new long[] {0, 3}, new long[] {1});
    checkarray(A, dm(new double[][] {{5, 7, 8}, {9, 11, 12}, {17, 19, 20}}), "remove a column and two rows ");
  }

  /** The six comparisons of A against c, both operand orders; ref = {<=,<,>,>=,==,!=} of A op c */
  private static void checkComparisons(DM A, double c, double[][][] ref) {
    DM C = new DM(c);
    List<BiFunction<DM, Double, DM>> lhsD = List.of((a, b) -> ca.le(a, b), (a, b) -> ca.lt(a, b),
        (a, b) -> ca.gt(a, b), (a, b) -> ca.ge(a, b), (a, b) -> ca.eq(a, b), (a, b) -> ca.ne(a, b));
    List<BiFunction<DM, DM, DM>> lhsM = List.of((a, b) -> a.le(b), (a, b) -> a.lt(b), (a, b) -> a.gt(b),
        (a, b) -> a.ge(b), (a, b) -> a.eq(b), (a, b) -> a.ne(b));
    List<BiFunction<Double, DM, DM>> rhsD = List.of((b, a) -> ca.ge(b, a), (b, a) -> ca.gt(b, a),
        (b, a) -> ca.lt(b, a), (b, a) -> ca.le(b, a), (b, a) -> ca.eq(b, a), (b, a) -> ca.ne(b, a));
    List<BiFunction<DM, DM, DM>> rhsM = List.of((b, a) -> ca.ge(b, a), (b, a) -> ca.gt(b, a),
        (b, a) -> ca.lt(b, a), (b, a) -> ca.le(b, a), (b, a) -> ca.eq(b, a), (b, a) -> ca.ne(b, a));
    List<BiFunction<DM, Double, DM>> lhsI = List.of((a, b) -> a.le(b), (a, b) -> a.lt(b), (a, b) -> a.gt(b),
        (a, b) -> a.ge(b), (a, b) -> a.eq(b), (a, b) -> a.ne(b));
    String[] names = {"<=", "<", ">", ">=", "==", "!="};
    for (int k = 0; k < 6; k++) {
      checkarray(lhsD.get(k).apply(A, c), new DM(ref[k]), "A" + names[k] + c);
      checkarray(lhsM.get(k).apply(A, C), new DM(ref[k]), "A" + names[k] + "DM(" + c + ")");
      checkarray(lhsI.get(k).apply(A, c), new DM(ref[k]), "A." + names[k] + "(" + c + ")");
      checkarray(rhsD.get(k).apply(c, A), new DM(ref[k]), c + " reversed " + names[k]);
      checkarray(rhsM.get(k).apply(C, A), new DM(ref[k]), "DM(" + c + ") reversed " + names[k]);
    }
  }

  @Test
  public void test_comparisons() {
    DM A = dm(new double[][] {{5, 4}, {2, 1}});
    double[][] T = {{1, 1}, {1, 1}}, F = {{0, 0}, {0, 0}};
    // python-only: np.array / np.matrix operands
    for (double c : new double[] {6, 6.0}) checkComparisons(A, c, new double[][][] {T, T, F, F, F, T});
    for (double c : new double[] {5, 5.0})
      checkComparisons(A, c, new double[][][] {T, {{0, 1}, {1, 1}}, F, {{1, 0}, {0, 0}}, {{1, 0}, {0, 0}},
          {{0, 1}, {1, 1}}});
    for (double c : new double[] {4, 4.0})
      checkComparisons(A, c, new double[][][] {{{0, 1}, {1, 1}}, {{0, 0}, {1, 1}}, {{1, 0}, {0, 0}},
          {{1, 1}, {0, 0}}, {{0, 1}, {0, 0}}, {{1, 0}, {1, 1}}});
    for (double c : new double[] {1, 1.0})
      checkComparisons(A, c, new double[][][] {{{0, 0}, {0, 1}}, F, {{1, 1}, {1, 0}}, T, {{0, 0}, {0, 1}},
          {{1, 1}, {1, 0}}});
    checkComparisons(A, 0, new double[][][] {F, F, T, T, F, T});
  }

  /** Python bool(x) */
  private static boolean truth(Object x) {
    if (x instanceof org.casadi.DM) return ((org.casadi.DM) x).nonzero();
    if (x instanceof org.casadi.SX) return ((org.casadi.SX) x).nonzero();
    return ((org.casadi.MX) x).nonzero();
  }

  @Test
  public void test_truth() {
    assertTrue(truth(dm(1)));
    assertFalse(truth(dm(0)));
    assertTrue(truth(dm(0.2)));
    assertTrue(truth(dm(-0.2)));
    assertThrows(RuntimeException.class, () -> truth(dm(2.0, 3)));
    assertThrows(RuntimeException.class, () -> truth(new DM()));
  }

  private static List<List<Long>> permutations(List<Long> base) {
    List<List<Long>> r = new ArrayList<>();
    if (base.isEmpty()) {
      r.add(new ArrayList<>());
      return r;
    }
    for (int i = 0; i < base.size(); i++) {
      List<Long> rest = new ArrayList<>(base);
      Long h = rest.remove(i);
      for (List<Long> p : permutations(rest)) {
        p.add(0, h);
        r.add(p);
      }
    }
    return r;
  }

  private static void checkListslice(DM d, long[] rowbase, long[] colbase) {
    List<Long> rb = new ArrayList<>(), cb = new ArrayList<>();
    for (long i : rowbase) rb.add(i);
    for (long i : colbase) cb.add(i);
    for (List<Long> col : permutations(cb))
      for (List<Long> row : permutations(rb)) {
        DM r = DM.zeros(row.size(), col.size());
        for (int i = 0; i < row.size(); i++)
          for (int j = 0; j < col.size(); j++) r.set(i, j, d.get(row.get(i), col.get(j)));
        checkarray(d.get(idx(toLongs(row)), idx(toLongs(col))), r, d + "[" + row + "," + col + "]");
      }
  }

  @Test
  public void test_listslice() {
    checkListslice(range(Sparsity.dense(3, 3)), new long[] {0, 1, 2}, new long[] {0, 1, 2});
    checkListslice(range(Sparsity.dense(4, 4)), new long[] {0, 1, 3}, new long[] {0, 2, 3});
    checkListslice(range(Sparsity.dense(3, 3)), new long[] {0, 0, 1}, new long[] {0, 0, 1});
    checkListslice(range(Sparsity.dense(3, 3)), new long[] {0, 0, 2}, new long[] {0, 0, 2});
    checkListslice(range(Sparsity.dense(3, 3)), new long[] {1, 1, 2}, new long[] {1, 1, 2});
    DM d = range(Sparsity.lower(4));
    checkListslice(d, new long[] {0, 1, 3}, new long[] {0, 2, 3});
    checkListslice(d.T(), new long[] {0, 1, 3}, new long[] {0, 2, 3});
    d = range(Sparsity.rowcol(new long[] {0, 1, 2}, new long[] {0, 1}, 4, 4));
    checkListslice(d, new long[] {0, 3}, new long[] {0, 2});
    checkListslice(range(Sparsity.dense(2, 2)), new long[] {0, 0, 0}, new long[] {0, 0, 0});
    checkListslice(range(Sparsity.dense(2, 2)), new long[] {0, 0, 1}, new long[] {0, 0, 1});
    checkListslice(range(Sparsity.dense(2, 2)), new long[] {1, 1, 0}, new long[] {1, 1, 0});
    checkListslice(range(Sparsity.dense(2, 2)), new long[] {1, 1, 1}, new long[] {1, 1, 1});
    d = range(Sparsity.lower(3));
    checkListslice(d, new long[] {0, 1, 2}, new long[] {0, 1, 2});
    checkListslice(d.T(), new long[] {0, 1, 2}, new long[] {0, 1, 2});
    d = range(Sparsity.rowcol(new long[] {0, 2}, new long[] {0, 1}, 4, 4));
    checkListslice(d, new long[] {0, 1, 3}, new long[] {0, 2, 3});
  }

  // test_sparsesym: returns immediately in Python (feature removed in 73f407e)

  @Test
  public void test_diagcat() {
    message("diagcat");
    DM C = ca.diagcat(List.of(dm(new double[][] {{-1.4, -3.2}, {-3.2, -28}}),
        dm(new double[][] {{15, -12, 2.1}, {-12, 16, -3.8}, {2.1, -3.8, 15}}), new DM(1.8), new DM(-4.0)));
    DM r = ca.sparsify(dm(new double[][] {{-1.4, -3.2, 0, 0, 0, 0, 0}, {-3.2, -28, 0, 0, 0, 0, 0},
        {0, 0, 15, -12, 2.1, 0, 0}, {0, 0, -12, 16, -3.8, 0, 0}, {0, 0, 2.1, -3.8, 15, 0, 0},
        {0, 0, 0, 0, 0, 1.8, 0}, {0, 0, 0, 0, 0, 0, -4}}));
    checkarray(C, r);
  }

  @Test
  public void test_diag_sparse() {
    message("diag sparse");
    for (double[] n : new double[][] {{0, 1, 0, 0, 2, 3, 4, 5, 6, 0}, {1, 2, 3, 0}, {0, 1, 2, 3}}) {
      DM d = ca.sparsify(new DM(n));
      DM D = new DM(n);
      DM m = ca.diag(d);
      DM M = ca.sparsify(ca.diag(D));
      assertArrayEquals(M.sparsity().colind(), m.sparsity().colind());
      assertArrayEquals(M.sparsity().row(), m.sparsity().row());
    }
  }

  @Test
  public void test_sprank() {
    message("sprank");
    assertEquals(3, ca.sprank(ca.sparsify(dm(new double[][] {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}))));
    assertEquals(2, ca.sprank(ca.sparsify(dm(new double[][] {{1, 0, 0}, {0, 0, 0}, {0, 0, 1}}))));
    assertEquals(1, ca.sprank(ca.sparsify(dm(new double[][] {{0, 0, 0}, {0, 0, 0}, {0, 0, 1}}))));
    assertEquals(0, ca.sprank(ca.sparsify(dm(new double[][] {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}}))));
    assertEquals(1, ca.sprank(DM.ones(1, 3)));
    assertEquals(1, ca.sprank(DM.ones(3, 1)));
    assertEquals(2, ca.sprank(DM.ones(2, 3)));
    assertEquals(2, ca.sprank(DM.ones(3, 2)));
    assertEquals(3, ca.sprank(DM.ones(3, 3)));
    assertEquals(3, ca.sprank(DM.ones(3, 3)));
    DM A = new DM(6, 4);
    A.set(0, 0, 1);
    A.set(1, 2, 1);
    A.set(2, 2, 1);
    A.set(5, 3, 1);
    assertEquals(3, ca.sprank(A));
  }

  @Test
  public void test_cross() {
    message("cross products");
    checkarray(ca.cross(dm(1, 0, 0), dm(0, 1, 0)), dm(0, 0, 1));
    checkarray(ca.cross(dm(1.1, 1.3, 1.7), dm(2, 3, 13)), dm(11.8, -10.9, 0.7));
    checkarray(ca.cross(dm(1.1, 1.3, 1.7).T(), dm(2, 3, 13).T()), dm(11.8, -10.9, 0.7).T());
    DM a = dm(new double[][] {{1.1, 1.3, 1.7}, {1, 0, 0}, {0, 0, 1}, {4, 5, 6}});
    DM b = dm(new double[][] {{2, 3, 13}, {0, 1, 0}, {0, 0, 1}, {1, 0, 1}});
    DM r = dm(new double[][] {{11.8, -10.9, 0.7}, {0, 0, 1}, {0, 0, 0}, {5, 2, -5}});
    checkarray(ca.cross(a, b), r);
    checkarray(ca.cross(a.T(), b.T()), r.T());
    checkarray(ca.cross(a, b, 2), r);
    checkarray(ca.cross(a.T(), b.T(), 1), r.T());
  }

  @Test
  public void test_is_regular() {
    assertTrue(dm(1, 2).is_regular());
    assertFalse(dm(1, Double.POSITIVE_INFINITY).is_regular());
    assertFalse(DM.nan(2).is_regular());
  }

  @Test
  public void test_sizes() {
    assertEquals(10, Sparsity.diag(10).nnz_diag());
    assertEquals(10, Sparsity.diag(10).nnz_upper());
    assertEquals(10, Sparsity.diag(10).nnz_lower());
    assertEquals(10 * 11 / 2, Sparsity.dense(10, 10).nnz_lower());
    assertEquals(10 * 11 / 2, Sparsity.dense(10, 10).nnz_upper());
    assertEquals(10, Sparsity.dense(10, 10).nnz_diag());
    DM a = ca.sparsify(dm(new double[][] {{1, 1, 0}, {1, 0, 1}, {0, 0, 0}}));
    assertEquals(1, a.nnz_diag());
    assertEquals(2, a.nnz_lower());
    assertEquals(3, a.nnz_upper());
  }

  @Test
  public void test_tril2symm() {
    DM a = range(Sparsity.upper(3)).T();
    DM s = ca.tril2symm(a);
    checkarray(s, dm(new double[][] {{0, 1, 3}, {1, 2, 4}, {3, 4, 5}}));
    assertThrows(RuntimeException.class, () -> ca.tril2symm(DM.ones(5, 3)));
    System.out.println(DM.ones(5, 5).nnz_upper() - DM.ones(5, 5).nnz_diag());
    assertThrows(RuntimeException.class, () -> ca.tril2symm(DM.ones(5, 5)));
  }

  @Test
  public void test_not_null() {
    MX x = MX.sym("x", 3, 1);
    Sparsity sp = Sparsity.upper(2);
    new MX(sp, x);
  }

  @Test
  public void test_segfault() {
    MX x = MX.sym("x", 10, 1);
    Sparsity sp = Sparsity.upper(2);
    MX y = ca.triu2symm(new MX(sp, x.get(new Slice(1, 4))));
    new org.casadi.Function("f", List.of(x), List.of(y));
  }

  @Test
  public void test_append_empty() {
    DM a = ca.vertcat(List.of(new DM(0, 0), new DM(0, 2)));
    assertEquals(0, a.size1());
    assertEquals(2, a.size2());
    a = ca.vertcat(List.of(new DM(0, 0), new DM(2, 0), new DM(3, 0)));
    assertEquals(5, a.size1());
    assertEquals(0, a.size2());
  }

  @Test
  public void test_vertcat_empty() {
    DM a = new DM(0, 2);
    DM v = ca.vertcat(List.of(a, a));
    assertEquals(0, v.size1());
    assertEquals(2, v.size2());
    a = new DM(2, 0);
    v = ca.vertcat(List.of(a, a));
    assertEquals(4, v.size1());
    assertEquals(0, v.size2());
  }

  @Test
  public void test_vertsplit() {
    DM a = range(Sparsity.upper(5)).T();
    List<DM> v = ca.vertsplit(a, new long[] {0, 2, 4, 5});
    assertEquals(3, v.size());
    checkarray(v.get(0), dm(new double[][] {{0, 0, 0, 0, 0}, {1, 2, 0, 0, 0}}));
    checkarray(v.get(1), dm(new double[][] {{3, 4, 5, 0, 0}, {6, 7, 8, 9, 0}}));
    checkarray(v.get(2), dm(new double[][] {{10, 11, 12, 13, 14}}));

    v = ca.vertsplit(a);
    assertEquals(a.size1(), v.size());
    checkarray(v.get(0), dm(new double[][] {{0, 0, 0, 0, 0}}));
    checkarray(v.get(1), dm(new double[][] {{1, 2, 0, 0, 0}}));
    checkarray(v.get(2), dm(new double[][] {{3, 4, 5, 0, 0}}));
    checkarray(v.get(3), dm(new double[][] {{6, 7, 8, 9, 0}}));
    checkarray(v.get(4), dm(new double[][] {{10, 11, 12, 13, 14}}));

    v = ca.vertsplit(a, new long[] {0, 2, 4, 5});
    assertEquals(3, v.size());
    checkarray(v.get(0), dm(new double[][] {{0, 0, 0, 0, 0}, {1, 2, 0, 0, 0}}));
    checkarray(v.get(1), dm(new double[][] {{3, 4, 5, 0, 0}, {6, 7, 8, 9, 0}}));
    checkarray(v.get(2), dm(new double[][] {{10, 11, 12, 13, 14}}));

    v = ca.vertsplit(a, new long[] {0, 0, 3, 5});
    assertEquals(3, v.size());
    assertEquals(0, v.get(0).size1());
    assertEquals(5, v.get(0).size2());
    checkarray(v.get(1), dm(new double[][] {{0, 0, 0, 0, 0}, {1, 2, 0, 0, 0}, {3, 4, 5, 0, 0}}));
    checkarray(v.get(2), dm(new double[][] {{6, 7, 8, 9, 0}, {10, 11, 12, 13, 14}}));
  }

  @Test
  public void test_horzsplit() {
    DM a = range(Sparsity.upper(5)).T();
    List<DM> v = ca.horzsplit(a, new long[] {0, 2, 4, 5});
    assertEquals(3, v.size());
    checkarray(v.get(0), dm(new double[][] {{0, 0}, {1, 2}, {3, 4}, {6, 7}, {10, 11}}));
    checkarray(v.get(1), dm(new double[][] {{0, 0}, {0, 0}, {5, 0}, {8, 9}, {12, 13}}));
    checkarray(v.get(2), dm(new double[][] {{0}, {0}, {0}, {0}, {14}}));

    v = ca.horzsplit(a);
    assertEquals(a.size1(), v.size());
    checkarray(v.get(0), dm(0, 1, 3, 6, 10));
    checkarray(v.get(1), dm(0, 2, 4, 7, 11));
    checkarray(v.get(2), dm(0, 0, 5, 8, 12));
    checkarray(v.get(3), dm(0, 0, 0, 9, 13));
    checkarray(v.get(4), dm(0, 0, 0, 0, 14));

    v = ca.horzsplit(a, new long[] {0, 2, 4, 5});
    assertEquals(3, v.size());
    checkarray(v.get(0), dm(new double[][] {{0, 0}, {1, 2}, {3, 4}, {6, 7}, {10, 11}}));
    checkarray(v.get(1), dm(new double[][] {{0, 0}, {0, 0}, {5, 0}, {8, 9}, {12, 13}}));
    checkarray(v.get(2), dm(new double[][] {{0}, {0}, {0}, {0}, {14}}));

    v = ca.horzsplit(a, new long[] {0, 0, 3, 5});
    assertEquals(3, v.size());
    assertEquals(5, v.get(0).size1());
    assertEquals(0, v.get(0).size2());
    checkarray(v.get(1), dm(new double[][] {{0, 0, 0}, {1, 2, 0}, {3, 4, 5}, {6, 7, 8}, {10, 11, 12}}));
    checkarray(v.get(2), dm(new double[][] {{0, 0}, {0, 0}, {0, 0}, {9, 0}, {13, 14}}));
  }

  @Test
  public void test_blocksplit() {
    DM a = range(Sparsity.upper(5)).T();
    List<List<DM>> v = ca.blocksplit(a, new long[] {0, 2, 4, 5}, new long[] {0, 1, 3, 5});
    checkarray(v.get(0).get(0), dm(0, 1));
    checkarray(v.get(0).get(1), dm(new double[][] {{0, 0}, {2, 0}}));
    checkarray(v.get(1).get(0), dm(3, 6));
    checkarray(ca.blockcat(v), a);
  }

  private static Sparsity spDiagcat(Sparsity... s) {
    return ca.diagcat(s);
  }

  private static Sparsity spVertcat(Sparsity... s) {
    return ca.vertcat(s);
  }

  @Test
  public void test_solve() {
    List<Sparsity> spA = new ArrayList<>(List.of(Sparsity.dense(1, 1)));
    for (long n = 2; n < 5; n++) {
      spA.addAll(List.of(Sparsity.diag(n), Sparsity.dense(n, n), Sparsity.lower(n), Sparsity.lower(n).T(),
          Sparsity.banded(n, 1), spDiagcat(Sparsity.diag(n), Sparsity.dense(n, n)),
          spDiagcat(Sparsity.diag(n), Sparsity.lower(n)),
          spDiagcat(Sparsity.diag(n), Sparsity.lower(n).T()),
          spDiagcat(Sparsity.lower(n), Sparsity.lower(n).T()),
          Sparsity.diag(n).plus(Sparsity.rowcol(new long[] {0}, new long[] {n - 1}, n, n)),
          Sparsity.diag(n).plus(Sparsity.rowcol(new long[] {0, n - 1}, new long[] {n - 1, 0}, n, n)),
          Sparsity.diag(n).plus(Sparsity.triplet(n, n, new long[] {0}, new long[] {n - 1})),
          Sparsity.diag(n).plus(Sparsity.triplet(n, n, new long[] {0, n - 1}, new long[] {n - 1, 0}))));
    }
    for (Sparsity sA : spA) {
      DM a = new DM(sA, new DM(randomValues(sA.nnz())));
      SX A = SX.sym("a", a.sparsity());
      long n1 = a.size1();
      for (Sparsity sB : List.of(Sparsity.dense(n1, 1), spVertcat(Sparsity.dense(1, 1), new Sparsity(n1 - 1, 1)),
          Sparsity.lower(n1), Sparsity.lower(n1).T())) {
        DM b = new DM(sB, new DM(randomValues(sB.nnz())));
        SX B = SX.sym("B", b.sparsity());
        SX C = ca.solve(A, B);
        org.casadi.Function f = new org.casadi.Function("f", List.of(A, B), List.of(C));
        DM c = (DM) f.invoke(a, b);
        DM c_ref = ca.sparsify(new DM(npSolve(a.full(), b.full())));
        System.out.println(sA.dim() + " " + sB.dim());
        System.out.println("foo" + c + c_ref);
        checkarray(c, c_ref);
      }
    }
  }

  private double[] randomValues(long n) {
    double[] v = new double[(int) n];
    for (int k = 0; k < n; k++) v[k] = random.nextDouble();
    return v;
  }

  @Test
  public void test_kron() {
    DM a = ca.sparsify(dm(new double[][] {{1, 0, 6}, {2, 7, 0}}));
    DM b = ca.sparsify(dm(new double[][] {{1, 0, 0}, {2, 3, 7}, {0, 0, 9}, {1, 12, 13}}));
    DM c_ = ca.kron(a, b);
    assertEquals(a.size1() * b.size1(), c_.size1());
    assertEquals(a.size2() * b.size2(), c_.size2());
    assertEquals(a.nnz() * b.nnz(), c_.nnz());
    checkarray(c_, npKron(a.full(), b.full()));
  }

  private static org.casadi.Function jacobian_old(org.casadi.Function f, long i, long j) {
    List<String> out = new ArrayList<>();
    out.add("jac:" + f.name_out(j) + ":" + f.name_in(i));
    out.addAll(f.name_out());
    return f.factory(f.name() + "_jac", f.name_in(), out);
  }

  @Test
  @SuppressWarnings("unchecked")
  public void test_vec_kron() {
    SX A = SX.sym("A", 2, 3);
    SX B = SX.sym("B", 4, 5);
    SX P = SX.sym("P", A.size2(), B.size1());
    org.casadi.Function f = new org.casadi.Function("f", List.of(ca.vec(P.T()), A, B),
        List.of(ca.vec(ca.mtimes(List.of(A, P, B)).T())));
    org.casadi.Function J = jacobian_old(f, 0, 0);
    List<DM> J_in = List.of(new DM(npRandom((int) ca.vec(P.T()).size1(), 1)), new DM(npRandom(2, 3)),
        new DM(npRandom(4, 5)));
    List<Object> r = (List<Object>) J.invoke(J_in.get(0), J_in.get(1), J_in.get(2));
    DM res = (DM) r.get(0);
    DM ref = ca.kron(J_in.get(1), J_in.get(2).T());
    checkarray(res, ref);
  }

  @Test
  public void test_repmat() {
    DM a = dm(new double[][] {{1, 2}, {3, 4}, {5, 6}});
    checkarray(ca.repmat(a, 2, 3), ca.kron(DM.ones(2, 3), a));
  }

  @Test
  public void test_triu() {
    DM a = dm(new double[][] {{1, 2}, {3, 4}});
    checkarray(ca.triu(a), dm(new double[][] {{1, 2}, {0, 4}}));
  }

  @Test
  public void test_tril() {
    DM a = dm(new double[][] {{1, 2}, {3, 4}});
    checkarray(ca.tril(a), dm(new double[][] {{1, 0}, {3, 4}}));
  }

  @Test
  public void test_nz() {
    DM a = ca.sparsify(dm(new double[][] {{1, 2}, {0, 0}, {3, 4}}));
    checkarray(a.get_nz(false, new Slice()), dm(1, 3, 2, 4));
    assertEquals(4, a.nnz());
    checkarray(a.get_nz(false, new Slice(0, -1)), dm(1, 3, 2));
    checkarray(a.get_nz(false, 0), dm(1));
  }

  // Java: get_nz/set_nz without the ind1 flag, like get/set
  @Test
  public void test_nz_without_ind1() {
    DM a = ca.sparsify(dm(new double[][] {{1, 2}, {0, 0}, {3, 4}}));
    checkarray(a.get_nz(new Slice(0, -1)), dm(1, 3, 2));
    checkarray(a.get_nz(dm(3, 0)), dm(4, 1));
    checkarray(a.get_nz(2), dm(2));
    a.set_nz(new Slice(0, 2), dm(7, 8));
    a.set_nz(3, new DM(9));
    a.set_nz(dm(2), new DM(6));
    checkarray(new DM(a.nonzeros()), dm(7, 8, 6, 9));
    MX x = MX.sym("x", 4);
    org.casadi.Function f = new org.casadi.Function("f", List.of(x), List.of(x.get_nz(new Slice(1, 3))));
    checkarray(f.call(List.of(dm(1, 2, 3, 4))).get(0), dm(2, 3));
    SX y = SX.sym("y", 3);
    y.set_nz(1, new SX(0));
    assertEquals("[y_0, 0, y_2]", y.toString());
  }

  private void zeroRandomEntries(double[][] A) {
    for (int i = 0; i < 5; i++) A[random.nextInt(A.length)][random.nextInt(A[0].length)] = 0;
  }

  @Test
  public void test_norm_inf_mul() {
    for (int[] s : new int[][] {{10, 2, 8}, {8, 2, 10}}) {
      double[][] A = npRandom(s[0], s[1]), B = npRandom(s[1], s[2]);
      checkarray(ca.norm_inf_mul(new DM(A), new DM(B)), ca.norm_inf(ca.mtimes(new DM(A), new DM(B))));
      checkarray(new DM(ca.norm_0_mul(new DM(A), new DM(B))), ca.mtimes(new DM(A), new DM(B)).nnz());
      zeroRandomEntries(A);
      zeroRandomEntries(B);
      DM As = ca.sparsify(new DM(A)), Bs = ca.sparsify(new DM(B));
      checkarray(ca.norm_inf_mul(As, Bs), ca.norm_inf(ca.mtimes(As, Bs)));
      checkarray(new DM(ca.norm_0_mul(As, Bs)), ca.mtimes(As, Bs).nnz());
    }
  }

  @Test
  public void test_mul3_issue_1465() {
    SX w2 = SX.sym("w", 2, 1);
    SX w3 = SX.sym("w", 3, 1);
    DM Q = DM.eye(2);
    ca.mtimes(w2.T(), Q);
    w2.T().mtimes(Q);
    assertInException("Matrix product with incompatible dimensions. Lhs is 1x3 and rhs is 2x2.",
        () -> ca.mtimes(w3.T(), Q));
    assertInException("Matrix product with incompatible dimensions. Lhs is 1x3 and rhs is 2x2.",
        () -> w3.T().mtimes(Q));
    SX Qs = new SX(Q);
    ca.mtimes(List.of(w2.T(), Qs, w2));
    assertInException("Matrix product with incompatible dimensions. Lhs is 1x3 and rhs is 2x2.",
        () -> ca.mtimes(List.of(w3.T(), Qs, w2)));
    assertInException("Matrix product with incompatible dimensions. Lhs is 1x2 and rhs is 3x1.",
        () -> ca.mtimes(List.of(w2.T(), Qs, w3)));
    assertInException("Matrix product with incompatible dimensions. Lhs is 1x3 and rhs is 2x2.",
        () -> ca.mtimes(List.of(w3.T(), Qs, w3)));
    // python-only: @ None and mtimes(x, y, non-string) are compile errors in Java
  }

  @Test
  public void test_chol() {
    for (int i = 0; i < 4; i++) {
      DM A = new DM(npRandom(3, 3));
      DM H = ca.mtimes(A, A.T());
      DM R = ca.chol(H);
      assertTrue(R.is_triu());
      checkarray(ca.mtimes(R.T(), R), H);
    }
  }

  @Test
  public void test_skew() {
    DM x = dm(1, 7, 13);
    checkarray(ca.inv_skew(ca.skew(x)), x);
    DM y = dm(0.2, 0.9, 0.4);
    checkarray(ca.mtimes(ca.skew(x), y), ca.cross(x, y));
  }

  @Test
  public void test_nz_overflow() {
    DM d = dm(2, 3);
    DM r = d.get_nz(false, new Slice());
    checkarray(r, d);
  }

  @Test
  public void test_DMcrash() {
    // Python DM([DM, DM]) list conversion; nearest Java constructor takes List<SX>
    assertThrows(RuntimeException.class, () -> new DM(List.of(new SX(dm(1, 2)), new SX(dm(1, 2)))));
    DM a = new DM(List.of(new SX(dm(1)), new SX(dm(2))));
    checkarray(a, dm(1, 2));
  }

  @Test
  public void test_sparsity_operation() {
    List<DM> L = List.of(new DM(1), new DM(new Sparsity(1, 1), 1), new DM(new Sparsity(2, 1), 1),
        new DM(Sparsity.dense(2, 1), 1));
    for (DM a : L)
      for (DM b : L) {
        DM c = a.times(b);
        if (a.nnz() == 0 || b.nnz() == 0) assertTrue(c.nnz() == 0);
        else assertTrue(c.nnz() > 0);
      }
    assertTrue(ca.sum2(new DM(new Sparsity(1, 1), 1)).nnz() == 0);
  }

  private static double[][] npCumsum(double[][] A, int axis) {
    int n = A.length, m = A[0].length;
    double[][] r = new double[n][m];
    for (int i = 0; i < n; i++)
      for (int j = 0; j < m; j++)
        r[i][j] = A[i][j] + (axis == 0 ? (i > 0 ? r[i - 1][j] : 0) : (j > 0 ? r[i][j - 1] : 0));
    return r;
  }

  @Test
  public void test_matlab_operations() {
    double[][][] data = {{{1, 3}, {11, 17}}, {{1, 3}}, {{1}, {3}}, {{3}}};
    for (double[][] A : data) {
      DM B = ca.reshape(new DM(A), new long[] {A.length, A[0].length});
      checkarray(ca.cumsum(B, 0), npCumsum(A, 0));
      checkarray(ca.cumsum(B, 1), npCumsum(A, 1));
      MX Bs = MX.sym("B", B.size());
      DM A0 = ca.cumsum(B, 0);
      DM A1 = ca.cumsum(B, 1);
      checkarray(ca.evalf(ca.cumsum(new MX(B), 0)), npCumsum(A, 0));
      checkarray(ca.evalf(ca.cumsum(B, 1)), npCumsum(A, 1));
      org.casadi.Function f = new org.casadi.Function("f", List.of(Bs), List.of(new MX(A0), new MX(A1)));
      checkfunction(f, f.expand(), List.of(B));
      checkCodegen(f, List.of(B), true);
      double[][] d = new double[A.length - 1][A[0].length];
      for (int i = 0; i + 1 < A.length; i++)
        for (int j = 0; j < A[0].length; j++) d[i][j] = A[i + 1][j] - A[i][j];
      DM D = ca.diff(B, 1, 0);
      assertEquals(d.length, D.size1());
      assertEquals(A[0].length, D.size2());
      if (d.length > 0) checkarray(D, d);
    }
  }

  @Test
  public void test_singular_repmat() {
    for (int n_b : new int[] {0, 2})
      for (int m_b : new int[] {0, 2})
        for (int n : new int[] {0, 3})
          for (int m : new int[] {0, 3}) {
            long[] ref = {n_b * n, m_b * m};
            assertArrayEquals(ref, ca.repmat(new DM(n_b, m_b), n, m).size());
            assertArrayEquals(ref, ca.repmat(new SX(n_b, m_b), n, m).size());
            assertArrayEquals(ref, ca.repmat(new MX(n_b, m_b), n, m).size());
            assertArrayEquals(ref, ca.repmat(new Sparsity(n_b, m_b), n, m).size());
          }
  }

  @Test
  public void test_serialize() {
    for (DM a : List.of(new DM(), new DM(2), dm(1, 2), dm(new double[][] {{1, 2}, {3, 4}, {5, 6}}))) {
      DM b = DM.deserialize(a.serialize());
      checkarray(a, b);
    }
  }

  // test_iterable: python-only (Python iteration protocol)

  @Test
  public void test_toarray() {
    // Java analog of toarray(): full() / elements(); numpy layout flags are python-only
    DM A = dm(new double[][] {{1, 2, 0}, {0, 4, 5}});
    for (DM a : List.of(A, ca.sparsify(A), ca.tril(dm(new double[][] {{1, 2}, {3, 4}})), dm(1, 2, 3),
        dm(1, 2, 3).T(), new DM(3), new DM(3, 0), new DM(0, 3))) {
      double[][] r = a.full();
      assertEquals(a.size1(), r.length);
      for (double[] ri : r) assertEquals(a.size2(), ri.length);
      for (int i = 0; i < r.length; i++)
        for (int j = 0; j < r[i].length; j++) assertEquals(a.get(i, j).to_double(), r[i][j]);
    }
    assertArrayEquals(new double[] {1, 2, 3}, dm(1, 2, 3).elements());
    assertEquals(3.0, new DM(3).to_double());
  }

  @Test
  public void test_from_file() throws IOException {
    Path dir = Files.createTempDirectory("casadi_from_file");
    String mtx = dir.resolve("test.mtx").toString(), txt = dir.resolve("test.txt").toString();
    Path txt2 = dir.resolve("test2.txt");
    DM a = ca.sparsify(dm(new double[][] {{1, 0, -6, 5, 0}, {4, 0, -4e-301, 9.3e-18, 0}}));
    a.print_dense();
    a.set(0, 4, 0);
    a.set(1, 1, 0);
    a.to_file(mtx);
    assertInException("foobar.mtx", () -> DM.from_file("foobar.mtx"));
    DM b = DM.from_file(mtx);
    assertTrue(a.sparsity().is_equal(b.sparsity()));
    checkarray(a, b);
    a.to_file(txt);
    b = DM.from_file(txt);
    checkarray(a, b);
    for (String s : new String[] {"1 00 -6 5 0\n4 0 -4e-301 9.3e-18 00\n", "1 00 -6 5 0\n4 0 -4e-301 9.3e-18 00",
        "1 00 -6 5 0 \n4 0 -4e-301 9.3e-18 00   ", "1\t00\t-6\t5\t0\t\n4 0 -4e-301 9.3e-18 00\n",
        "1 00 -6 5 0\n%foobar\n4 0 -4e-301 9.3e-18 00\n"}) {
      Files.write(txt2, s.getBytes(StandardCharsets.UTF_8));
      b = DM.from_file(txt2.toString());
      assertTrue(a.sparsity().is_equal(b.sparsity()));
      assertEquals(0.0, ca.norm_1(a.minus(b)).to_double());
    }
    Files.write(Path.of(txt), "1 aa 6\n4 8 9e-18\n".getBytes(StandardCharsets.UTF_8));
    assertInException("Parsing", () -> DM.from_file(txt));
    a = dm(3, Double.POSITIVE_INFINITY, Double.NEGATIVE_INFINITY, Double.NaN);
    a.to_file(txt);
    b = DM.from_file(txt);
    checkarray(a, b);
  }

  @Test
  public void test_norm_2() {
    double[] a = {1, 2, 3};
    double r = Math.sqrt(14);
    checkarray(ca.evalf(ca.norm_2(new DM(a))), r);
    checkarray(ca.evalf(ca.norm_2(new DM(a).T())), r);
    checkarray(ca.evalf(ca.norm_2(new MX(new DM(a)))), r);
    checkarray(ca.evalf(ca.norm_2(new MX(new DM(a)).T())), r);
    checkarray(ca.evalf(ca.norm_2(new SX(a))), r);
    checkarray(ca.evalf(ca.norm_2(new SX(a).T())), r);
  }

  @Test
  public void test_ldl() {
    DM.rng(1);
    DM H = ca.diagcat(List.of(DM.rand(5, 5), DM.rand(5, 5)));
    H = H.plus(H.T()).plus(DM.eye(10).times(2));
    DM p0 = idx(permutation(10));
    H = H.get(p0, p0);
    List<Object> r = ca.ldl(H);
    DM D = (DM) r.get(0), Lt = (DM) r.get(1);
    long[] p = (long[]) r.get(2);
    DM P = DM.eye(10).get(new Slice(), idx(p));
    DM F = ca.mtimes(List.of(ca.sqrt(ca.diag(D)), DM.eye(10).plus(Lt), P.T()));
    System.out.println(H.minus(ca.mtimes(F.T(), F)));
    assertTrue(ca.norm_fro(H.minus(ca.mtimes(F.T(), F))).to_double() <= 1e-14);
  }

  @Test
  public void test_im_bugs() {
    // Python vertcat(1, 2): Java vertcat takes matrices only
    Object a = ca.vertcat(List.of(new DM(1), new DM(2)));
    assertTrue(a instanceof DM);
    checkarray(ca.linspace(1, 3, 10), ca.linspace(1.0, 3.0, 10));
  }

  private static double[] npLinspace(double a, double b, int n) {
    double[] r = new double[n];
    double step = (b - a) / (n - 1);
    for (int i = 0; i < n; i++) r[i] = step == 0 ? i / (double) (n - 1) * (b - a) + a : i * step + a;
    if (n > 1) r[n - 1] = b;
    return r;
  }

  @Test
  public void test_linspace() {
    double[][] cases = {{1.0, 3.0}, {-2.5, 7.25}, {0.0, 1.0}, {1e10, 1e-10}, {-1.0, 1.0}};
    int[] sizes = {2, 3, 5, 10, 100, 1001};
    for (double[] ab : cases) {
      double a = ab[0], b = ab[1];
      for (int n : sizes) {
        double[] ref = npLinspace(a, b, n);
        String w = " a=" + a + " b=" + b + " n=" + n;
        assertArrayEquals(ref, ca.linspace(new DM(a), new DM(b), n).nonzeros(), "DM linspace mismatch" + w);
        MX x = MX.sym("x"), y = MX.sym("y");
        org.casadi.Function f = new org.casadi.Function("f", List.of(x, y), List.of(ca.linspace(x, y, n)));
        assertArrayEquals(ref, ((DM) f.invoke(a, b)).nonzeros(), "MX linspace mismatch" + w);
        SX xs = SX.sym("x"), ys = SX.sym("y");
        f = new org.casadi.Function("f", List.of(xs, ys), List.of(ca.linspace(xs, ys, n)));
        assertArrayEquals(ref, ((DM) f.invoke(a, b)).nonzeros(), "SX linspace mismatch" + w);
      }
    }
  }

  @Test
  public void test_linspace_mx_node_count() {
    MX x = MX.sym("x"), y = MX.sym("y");
    long n_small = ca.n_nodes(ca.linspace(x, y, 5));
    long n_large = ca.n_nodes(ca.linspace(x, y, 1000));
    assertTrue(n_large < n_small + 5, "MX linspace node count grew with nsteps: " + n_small + " -> " + n_large);
  }

  @Test
  public void test_permutation() {
    int n = 10;
    long[] p = permutation(n);
    Sparsity S = Sparsity.permutation(p);
    assertArrayEquals(p, S.permutation_vector());
    assertTrue(S.is_permutation());
    DM v = DM.rand(n);
    // Python S @ v converts the Sparsity to DM.ones(S)
    checkarray(ca.mtimes(DM.ones(S), v), v.get(idx(p)));
    S = Sparsity.permutation(p, true);
    checkarray(ca.mtimes(DM.ones(S), v).get(idx(p)), v);
  }

  @Test
  public void test_sparsity_orthonormality() {
    Map<String, Predicate<Sparsity>> alltests = new LinkedHashMap<>();
    alltests.put("is_orthonormal_rows(true)", s -> s.is_orthonormal_rows(true));
    alltests.put("is_orthonormal_columns(true)", s -> s.is_orthonormal_columns(true));
    alltests.put("is_selection(true)", s -> s.is_selection(true));
    alltests.put("is_orthonormal_rows(false)", s -> s.is_orthonormal_rows(false));
    alltests.put("is_orthonormal_columns(false)", s -> s.is_orthonormal_columns(false));
    alltests.put("is_selection(false)", s -> s.is_selection(false));
    alltests.put("is_permutation()", s -> s.is_permutation());
    alltests.put("is_orthonormal(false)", s -> s.is_orthonormal(false));
    alltests.put("is_orthonormal(true)", s -> s.is_orthonormal(true));
    BiFunction<Sparsity, List<String>, Void> check = (S, prop) -> {
      for (String p : prop) assertTrue(alltests.containsKey(p));
      for (Map.Entry<String, Predicate<Sparsity>> t : alltests.entrySet()) {
        boolean quality = prop.contains(t.getKey());
        assertEquals(quality, t.getValue().test(S), "Call " + t.getKey() + " expected to return " + quality);
      }
      return null;
    };
    Function<double[][], Sparsity> sp = m -> ca.sparsify(new DM(m)).sparsity();
    List<String> sel = List.of("is_orthonormal_rows(true)", "is_orthonormal_columns(true)", "is_selection(true)",
        "is_orthonormal(true)");
    for (double[][] S : new double[][][] {
        {{0, 0, 0, 0, 0}, {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0}},
        {{0, 0, 0, 0, 0}, {0, 1, 0, 0, 0}, {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0}},
        {{0, 0, 0, 0, 0}, {0, 1, 0, 0, 0}, {0, 0, 0, 1, 0}, {0, 0, 0, 0, 0}}})
      check.apply(sp.apply(S), sel);
    for (double[][] S : new double[][][] {
        {{0, 0, 0, 0, 0}, {0, 1, 1, 0, 0}, {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0}},
        {{0, 1, 0, 0}, {1, 0, 0, 0}, {0, 1, 1, 0}, {0, 0, 0, 0}, {0, 0, 0, 1}},
        {{0, 0, 0, 0, 0}, {0, 1, 0, 0, 0}, {0, 1, 0, 0, 0}, {0, 0, 0, 0, 0}}})
      check.apply(sp.apply(S), List.of());
    check.apply(sp.apply(new double[][] {{0, 1, 0, 0, 0}, {1, 0, 0, 0, 0}, {0, 0, 1, 0, 0}, {0, 0, 0, 1, 0}}),
        List.of("is_orthonormal_rows(true)", "is_orthonormal_columns(true)", "is_selection(true)",
            "is_orthonormal_rows(false)", "is_selection(false)", "is_orthonormal(true)"));
    check.apply(sp.apply(new double[][] {{0, 1, 0, 0}, {1, 0, 0, 0}, {0, 0, 1, 0}, {0, 0, 0, 1}}),
        new ArrayList<>(alltests.keySet()));
    check.apply(sp.apply(new double[][] {{0, 1, 0, 0}, {1, 0, 0, 0}, {0, 0, 1, 0}, {0, 0, 0, 0}, {0, 0, 0, 1}}),
        List.of("is_orthonormal_rows(true)", "is_orthonormal_columns(true)", "is_selection(true)",
            "is_orthonormal_columns(false)", "is_orthonormal(true)"));
  }

  private static void generateAndDelete(org.casadi.Function f) {
    f.generate("f.c");
    try {
      Files.deleteIfExists(Path.of("f.c"));
    } catch (IOException e) {
      throw new RuntimeException(e);
    }
  }

  @Test
  public void test_permutation_solve() {
    int n = 10;
    DM.rng(1);
    long[] p = permutation(n);
    Sparsity P = Sparsity.permutation(p);
    SX xs = SX.sym("x", P), bs = SX.sym("b", n);
    assertTrue(xs.sparsity().is_orthonormal());
    ca.sparsity_cast(xs, Sparsity.dense(n));
    org.casadi.Function f0 = new org.casadi.Function("f", List.of(xs, bs), List.of(ca.solve(xs, bs)));
    generateAndDelete(f0);
    MX xm = MX.sym("x", P), bm = MX.sym("b", n);
    assertTrue(xm.sparsity().is_orthonormal());
    ca.sparsity_cast(xm, Sparsity.dense(n));
    org.casadi.Function f1 = new org.casadi.Function("f", List.of(xm, bm), List.of(ca.solve(xm, bm)));
    generateAndDelete(f1);
    DM x0 = DM.rand(xm.sparsity());
    DM b0 = DM.rand(bm.sparsity());
    checkarray((DM) f0.invoke(x0, b0), (DM) f1.invoke(x0, b0));
  }

  @Test
  public void test_permutation_inv() {
    int n = 10;
    DM.rng(1);
    long[] p = permutation(n);
    Sparsity P = Sparsity.permutation(p);
    SX xs = SX.sym("x", P);
    org.casadi.Function f0 = new org.casadi.Function("f", List.of(xs), List.of(ca.solve(xs, SX.eye(n))));
    generateAndDelete(f0);
    MX xm = MX.sym("x", P);
    org.casadi.Function f1 = new org.casadi.Function("f", List.of(xm), List.of(ca.solve(xm, MX.eye(n))));
    generateAndDelete(f1);
    DM x0 = new DM(P, DM.rand(n));
    checkarray((DM) f0.invoke(x0), (DM) f1.invoke(x0));

    f0 = new org.casadi.Function("f", List.of(xs), List.of(ca.inv(xs)));
    f0.disp(true);
    f1 = new org.casadi.Function("f", List.of(xm), List.of(ca.inv(xm)));
    f1.disp(true);
    x0 = new DM(P, DM.rand(n));
    checkarray((DM) f0.invoke(x0), (DM) f1.invoke(x0));
  }

  @Test
  public void test_mac_simplify() {
    int n = 10;
    DM.rng(1);
    long[] p = permutation(n);
    Sparsity P = Sparsity.permutation(p);
    SX bs = SX.sym("b", n), xs = SX.sym("x", P);
    org.casadi.Function f0 = new org.casadi.Function("f", List.of(xs, bs), List.of(ca.mtimes(xs, bs)));
    generateAndDelete(f0);
    MX bm = MX.sym("b", n), xm = MX.sym("x", P);
    org.casadi.Function f1 = new org.casadi.Function("f", List.of(xm, bm), List.of(ca.mtimes(xm, bm)));
    generateAndDelete(f1);
    DM x0 = new DM(P, DM.rand(n));
    DM b0 = DM.rand(bm.sparsity());
    checkarray((DM) f0.invoke(x0, b0), (DM) f1.invoke(x0, b0));
  }

  @Test
  public void test_logsumexp_empty() {
    assertInException("Argument must be nonempty", () -> ca.logsumexp(new DM(0, 1)));
  }

  @Test
  public void test_linspace_nsteps() {
    for (long n : new long[] {0, -1})
      assertInException("must be at least 1", () -> ca.linspace(new DM(0), new DM(1), n));
  }

  @Test
  public void test_densify_value() {
    checkarray(ca.densify(new DM(Sparsity.diag(2), 1), 5), dm(new double[][] {{1, 5}, {5, 1}}));
  }

  @Test
  public void test_polyval_sparse() {
    checkarray(ca.polyval(dm(1, 2), new DM(Sparsity.diag(2), 1)), dm(new double[][] {{3, 2}, {2, 3}}));
    MX x = MX.sym("x");
    org.casadi.Function f = new org.casadi.Function("f", List.of(x),
        List.of(ca.polyval(new MX(dm(new double[][] {{1, 2}})), x)));
    checkarray((DM) f.invoke(3), new DM(5));
  }

  @Test
  public void test_det_empty_column() {
    checkarray(ca.det(new DM(Sparsity.triplet(2, 2, new long[] {0, 1}, new long[] {0, 0}), dm(1, 2))), new DM(0));
  }

  @Test
  public void test_negative_size_checks() {
    assertInException("must be a multiple of n", () -> ca.repsum(dm(1, 2, 3), 2));
    assertInException("Dimensions must be non-negative", () -> MX.sym("x", -1, 1));
  }

  @Test
  public void test_inplace() {
    int n = 10;
    MX x = MX.sym("x", n);
    MX y = new MX(Sparsity.diag(n), x);
    org.casadi.Function f = new org.casadi.Function("f", List.of(x), List.of(y));
    assertTrue(f.sz_w() == n);
  }

  @Test
  public void test_split_n_zero() {
    DM x = DM.ones(2, 4);
    assertInException("must be a multiple of n", () -> ca.horzsplit_n(x, 0));
    assertInException("must be a multiple of n", () -> ca.vertsplit_n(x.T(), 0));
  }

  @Test
  public void test_diagsplit_offset_lengths() {
    assertInException("offset1 and offset2 must have the same length", () -> ca.diagsplit(DM.ones(4, 4), 1, 2));
  }

  @Test
  public void test_vertcat_columns() {
    assertInException("Mismatching number of columns", () -> ca.vertcat(List.of(DM.ones(1, 2), DM.ones(1, 3))));
    assertInException("Mismatching number of columns",
        () -> ca.blockcat(List.of(List.of(new DM(1), new DM(2)), List.of(new DM(3)))));
    assertArrayEquals(new long[] {2, 2},
        ca.blockcat(List.of(List.of(new DM(1), new DM(2)), List.of(DM.ones(1, 2)))).size());
  }

  @Test
  public void test_split_user_errors() {
    DM x = DM.ones(2, 4);
    for (long[] offset : new long[][] {{1, 4}, {0, 3, 2, 4}})
      assertInException("must start at 0, be nondecreasing and end at 4", () -> ca.horzsplit(x, offset));
    for (long incr : new long[] {0, -1}) assertInException("must be positive", () -> ca.horzsplit(x, incr));
    assertInException("must be positive", () -> ca.blocksplit(DM.ones(4, 4), 0, 2));
    assertInException("must start at 0, be nondecreasing and end at 4",
        () -> ca.diagsplit(DM.ones(4, 4), new long[] {0, 1, 3}, new long[] {0, 2, 4}));
    assertInException("Reshape mismatch", () -> ca.reshape(DM.eye(2), Sparsity.dense(4, 2)));
    assertInException("sparsity_cast", () -> ca.sparsity_cast(DM.eye(2), Sparsity.dense(3, 1)));
  }

  @Test
  public void test_repmat_negative() {
    assertInException("must be non-negative", () -> ca.repmat(dm(1, 2), -1));
  }

  @Test
  public void test_blockcat_empty_diagsplit_offblock() {
    assertArrayEquals(new long[] {0, 0}, ca.blockcat(new ArrayList<List<DM>>()).size());
    assertArrayEquals(new long[] {0, 0}, ca.blockcat(List.of(new ArrayList<DM>())).size());
    assertArrayEquals(new long[] {1, 1}, ca.blockcat(List.of(new ArrayList<DM>(), List.of(new DM(1)))).size());
    DM X1 = dm(new double[][] {{1, 2}, {3, 4}});
    assertInException("outside the diagonal blocks", () -> ca.diagsplit(X1, 1));
    assertInException("outside the diagonal blocks", () -> ca.diagsplit(SX.sym("x", 2, 2), 1));
    assertInException("outside the diagonal blocks", () -> ca.diagsplit(Sparsity.dense(2, 2), 1));
    assertInException("outside the diagonal blocks is unsupported", () -> ca.diagsplit(MX.sym("x", 2, 2), 1));
  }

  @Test
  public void test_triplet_empty() {
    assertArrayEquals(new long[] {0, 0}, DM.triplet(new long[] {}, new long[] {}, new DM(0, 1)).size());
  }

  @Test
  public void test_triplet_default_size() {
    DM A = DM.triplet(new long[] {0, 1}, new long[] {1, 0}, dm(5, 6));
    assertArrayEquals(new long[] {2, 2}, A.size());
    checkarray(A, dm(new double[][] {{0, 5}, {6, 0}}));
  }

  @Test
  public void test_eye_negative() {
    assertInException("Dimensions must be non-negative", () -> DM.eye(-1));
    assertInException("Dimensions must be non-negative", () -> SX.eye(-1));
    assertInException("Dimensions must be non-negative", () -> MX.eye(-1));
  }

  @Test
  public void test_horzsplit_n() {
    assertTrue(ca.horzsplit_n(DM.rand(1, 6), 2).size() == 2);
    assertTrue(ca.vertsplit_n(DM.rand(6), 2).size() == 2);
    assertTrue(ca.horzsplit_n(new DM(0, 0), 2).size() == 2);
    assertTrue(ca.vertsplit_n(new DM(0, 0), 2).size() == 2);
    assertThrows(RuntimeException.class, () -> ca.horzsplit_n(DM.rand(1, 5), 2));
    assertThrows(RuntimeException.class, () -> ca.vertsplit_n(DM.rand(5), 2));
  }
}

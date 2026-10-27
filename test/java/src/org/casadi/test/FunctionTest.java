package org.casadi.test;

import static org.junit.jupiter.api.Assertions.assertArrayEquals;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import org.casadi.Callback;
import org.casadi.CodeGenerator;
import org.casadi.DM;
import org.casadi.Function;
import org.casadi.GenericMatrix;
import org.casadi.GlobalOptions;
import org.casadi.MX;
import org.casadi.SX;
import org.casadi.Slice;
import org.casadi.Sparsity;
import org.casadi.ca;
import org.junit.jupiter.api.Assumptions;
import org.junit.jupiter.api.Test;

/** Port of test/python/function.py */
public class FunctionTest extends CasadiTestCase {

  private static final long END = Long.MAX_VALUE;

  private static Function jacobianOld(Function f, long i, long j) {
    List<String> out = new ArrayList<>();
    out.add("jac:" + f.name_out(j) + ":" + f.name_in(i));
    out.addAll(f.name_out());
    return f.factory(f.name() + "_jac", f.name_in(), out);
  }

  private static Function hessianOld(Function f, long i, long j) {
    List<String> out = new ArrayList<>();
    out.add("hess:" + f.name_out(j) + ":" + f.name_in(i) + ":" + f.name_in(i));
    out.add("grad:" + f.name_out(j) + ":" + f.name_in(i));
    out.addAll(f.name_out());
    return f.factory(f.name() + "_hess", f.name_in(), out);
  }

  private static Function toSXFun(Function f) {
    List<SX> ins = f.sx_in();
    return new Function("f", ins, f.call(ins));
  }

  /** Inputs as DM: numbers, double[] (column), double[][] (rows) or DM */
  private static List<DM> ins(Object... a) {
    List<DM> r = new ArrayList<>();
    for (Object o : a) {
      if (o instanceof DM) r.add((DM) o);
      else if (o instanceof Number) r.add(new DM(((Number) o).doubleValue()));
      else if (o instanceof double[]) r.add(new DM((double[]) o));
      else r.add(new DM((double[][]) o));
    }
    return r;
  }

  /** As np.random.random: uniform nonzeros on sp */
  private DM rnd(Sparsity sp) {
    double[] v = new double[(int) sp.nnz()];
    for (int k = 0; k < v.length; k++) v[k] = random.nextDouble();
    return new DM(sp, new DM(v));
  }

  private DM rnd(long n, long m) {
    return rnd(Sparsity.dense(n, m));
  }

  private static double[] linspace(double a, double b, int n) {
    double[] r = new double[n];
    for (int i = 0; i < n; i++) r[i] = n == 1 ? a : a + (b - a) * i / (n - 1);
    return r;
  }

  private static DM d(Object o) {
    return (DM) o;
  }

  private static MX m(Object o) {
    return (MX) o;
  }

  private static SX s(Object o) {
    return (SX) o;
  }

  @SuppressWarnings("unchecked")
  private static <T> List<T> l(Object o) {
    return (List<T>) o;
  }

  private static <T> List<T> cat(List<T> a, List<T> b) {
    List<T> r = new ArrayList<>(a);
    r.addAll(b);
    return r;
  }

  @SafeVarargs
  private static <T> List<T> cat(List<T>... a) {
    List<T> r = new ArrayList<>();
    for (List<T> e : a) r.addAll(e);
    return r;
  }

  private static String captureStderr(Runnable r) {
    PrintStream old = System.err;
    ByteArrayOutputStream b = new ByteArrayOutputStream();
    System.setErr(new PrintStream(b, true, StandardCharsets.UTF_8));
    try {
      r.run();
    } finally {
      System.err.flush();
      System.setErr(old);
    }
    return b.toString(StandardCharsets.UTF_8);
  }

  private void checkfunctionLight(Function f, Function g, List<DM> inputs) {
    checkfun(f, g, inputs, false, false, false, 1e-9);
  }

  /** checkfunction with the jacobian/forward/reverse parts individually selectable */
  private void checkfun(Function f, Function g, List<DM> inputs, boolean jac, boolean fwd, boolean adj, double tol) {
    assertEquals(g.n_in(), f.n_in(), "n_in");
    assertEquals(g.n_out(), f.n_out(), "n_out");
    List<DM> a = f.call(inputs), b = g.call(inputs);
    for (int i = 0; i < a.size(); i++) checkarray(a.get(i), b.get(i), tol, "output " + i);
    if (jac) {
      List<DM> va = f.jacobian().call(cat(inputs, a)), vb = g.jacobian().call(cat(inputs, b));
      for (int i = 0; i < va.size(); i++) checkarray(va.get(i), vb.get(i), tol, "jacobian block " + i);
    }
    if (fwd) {
      List<DM> s = new ArrayList<>(inputs);
      s.addAll(a);
      for (long i = 0; i < f.n_in(); i++) s.add(rnd(f.sparsity_in(i)));
      List<DM> fa = f.forward(1).call(s), fb = g.forward(1).call(s);
      for (int i = 0; i < fa.size(); i++) checkarray(fa.get(i), fb.get(i), tol, "forward " + i);
    }
    if (adj) {
      List<DM> s = new ArrayList<>(inputs);
      s.addAll(a);
      for (long i = 0; i < f.n_out(); i++) s.add(rnd(f.sparsity_out(i)));
      List<DM> ra = f.reverse(1).call(s), rb = g.reverse(1).call(s);
      for (int i = 0; i < ra.size(); i++) checkarray(ra.get(i), rb.get(i), tol, "reverse " + i);
    }
  }

  private void checkSerialize(Function f, List<DM> inputs) {
    Function f2 = Function.deserialize(f.serialize(Map.of("debug", true)));
    List<DM> a = f.call(inputs), b = f2.call(inputs);
    for (int i = 0; i < f.n_out(); i++) checkarray(a.get(i), b.get(i), 1e-15, "serialize output " + i);
  }

  private void checkSerialize(Function f, Map<String, DM> inputs) {
    Function f2 = Function.deserialize(f.serialize(Map.of("debug", true)));
    Map<String, DM> a = f.call(inputs), b = f2.call(inputs);
    assertEquals(f.name_out(), f2.name_out());
    for (String k : f.name_out()) checkarray(a.get(k), b.get(k), 1e-15, "serialize output " + k);
  }

  private static int cgCounter = 0;

  private static void sh(String cmd, File dir, File in, File out) {
    try {
      ProcessBuilder pb = new ProcessBuilder("sh", "-c", cmd).directory(dir);
      if (in != null) pb.redirectInput(in);
      if (out != null) pb.redirectOutput(out);
      else pb.inheritIO();
      int rc = pb.start().waitFor();
      assertEquals(0, rc, "command failed: " + cmd);
    } catch (IOException | InterruptedException e) {
      throw new RuntimeException(e);
    }
  }

  private Function checkCodegen(Function f, List<DM> inputs) {
    return checkCodegen(f, inputs, Map.of(), Map.of(), false, "c89", false);
  }

  /** As helpers.check_codegen: only active with RUN_SLOW (needs gcc); returns the loaded external or null */
  private Function checkCodegen(Function f, List<DM> inputs, Map<String, ?> opts, Map<String, ?> extOpts,
      boolean withJacSparsity, String std, boolean main) {
    if (!RUN_SLOW) return null;
    try {
      File dir = Files.createTempDirectory("casadi_java_cg").toFile();
      String name = "codegen_" + (cgCounter++) + "_" + Math.abs(random.nextLong());
      Map<String, Object> o = new LinkedHashMap<>(opts);
      if (main) o.put("main", true);
      CodeGenerator cg = new CodeGenerator(name, o);
      cg.add(f, withJacSparsity);
      cg.generate(dir.getPath() + "/");
      String lib = GlobalOptions.getCasadiPath(), inc = GlobalOptions.getCasadiIncludePath();
      String flags = "gcc -pedantic -std=" + std + " -fPIC -Wall -Werror -Wextra -I" + inc + " -I" + inc
          + "/highs -Wno-dangling-pointer -Wno-unknown-warning-option -Wno-unknown-pragmas -Wno-long-long"
          + " -Wno-unused-parameter -O3 " + name + ".c";
      String libs = " -L" + lib + " -Wl,-rpath," + lib + " -Wl,-rpath,.";
      sh(flags + " -shared -o " + name + ".so" + libs, dir, null, null);
      List<DM> ref = f.call(inputs);
      if (main) {
        sh(flags + " -o " + name + libs + " -lm", dir, null, null);
        File in = new File(dir, f.name() + "_in.txt"), out = new File(dir, f.name() + "_out.txt");
        f.generate_in(in.getPath(), inputs);
        sh("./" + name + " " + f.name(), dir, in, out);
        List<DM> outputs = f.generate_out(out.getPath());
        for (int i = 0; i < f.n_out(); i++) checkarray(ref.get(i), outputs.get(i), 1e-14, "main output " + i);
      }
      Function f2 = ca.external(f.name(), dir.getPath() + "/" + name + ".so", extOpts);
      List<DM> res = f2.call(inputs);
      for (int i = 0; i < f.n_out(); i++) checkarray(ref.get(i), res.get(i), 1e-14, "codegen output " + i);
      return f2;
    } catch (IOException e) {
      throw new RuntimeException(e);
    }
  }

  private void checkCodegen(Function f, Map<String, DM> inputs, boolean main) {
    checkCodegen(f, f.convert_in(inputs), Map.of(), Map.of(), false, "c89", main);
  }

  @Test
  public void test_call_empty() {
    SX xs = SX.sym("x", 2);
    Function fsx = new Function("fsx", List.of(xs, new SX()), List.of(xs));
    MX x = MX.sym("x", 2);
    Function fmx1 = new Function("fmx1", List.of(x, new MX()), List.of(x));
    Function fmx2 = new Function("fmx2", List.of(x, new MX()), List.of(x));

    for (Function f : List.of(fsx, fmx1, fmx2)) {
      f.invoke(0, 0);
      MX X = MX.sym("X", 2);
      MX F = m(f.invoke(X, new MX()));
      Function g = new Function("g", List.of(X), List.of(F));
      g.invoke(0);
    }

    xs = SX.sym("x", 2);
    fsx = new Function("fsx", List.of(xs), List.of(xs, new SX()));
    x = MX.sym("x", 2);
    fmx1 = new Function("fmx1", List.of(x), List.of(x, new MX()));
    fmx2 = new Function("fmx2", List.of(x), List.of(x, new MX()));

    for (Function f : List.of(fsx, fmx1)) {
      f.invoke(0);
      MX X = MX.sym("X", 2);
      List<MX> F = l(f.invoke(X));
      Function g = new Function("g", List.of(X), F);
      g.invoke(0);
    }
  }

  @Test
  public void test_MX_funSeed() {
    message("MX_funSeed");
    MX x1 = MX.sym("x", 2), y1 = MX.sym("y"), x2 = MX.sym("x", 2), y2 = MX.sym("y");
    Function p = new Function("p", List.of(x1, y1, x2, y2), List.of(ca.sin(x1).plus(y1), ca.sin(x2).plus(y2)));
    DM n1 = dm(4, 5), n2 = dm(5, 7);
    double N1 = 3, N2 = 8;
    List<DM> out = l(p.invoke(n1, N1, n2, N2));
    checkarray(ca.sin(n1).plus(N1), out.get(0), "output");
    checkarray(ca.sin(n2).plus(N2), out.get(1), "output");
  }

  @Test
  public void test_segfault() {
    Function f = new Function();
    assertThrows(RuntimeException.class, () -> f.stats());
  }

  @Test
  public void test_issue304() {
    message("regression test for #304");
    SX x = SX.sym("x");
    Function f = new Function("f", List.of(x), List.of(x.power(2.0), x.power(3.0)));
    MX X = MX.sym("X");
    List<MX> z = l(f.invoke(X));
    new Function("g", List.of(X), z).expand();
  }

  @Test
  public void test_jacobian() {
    SX x = SX.sym("x", 3, 1), y = SX.sym("y", 2, 1);
    Function f = new Function("f", List.of(x, y), List.of(x.power(2.0), y, x.times(y.get(0))));
    Function g = jacobianOld(f, 0, 0);
    assertEquals(f.n_in(), g.n_in());
    assertEquals(f.n_out() + 1, g.n_out());
  }

  @Test
  public void test_xfunction() {
    SX x = SX.sym("x", 3, 1), y = SX.sym("y", 2, 1);
    Function f = new Function("f", List.of(x, y), List.of(x.power(2.0), y, x.times(y.get(0))));
    MX X = MX.sym("x", 3, 1), Y = MX.sym("y", 2, 1);
    Function F = new Function("F", List.of(X, Y), List.of(X.power(2.0), Y, X.times(Y.get(0))));
    checkfun(f, F, ins(new double[] {0.1, 0.7, 1.3}, new double[] {7.1, 2.9}), true, false, false, 1e-9);
  }

  private static void jacobiansTest(Sparsity sp) {
    SX x = SX.sym("x", sp.size2());
    Sparsity sp2 = ca.jacobian(ca.mtimes(DM.ones(sp), x), x).sparsity();
    assertArrayEquals(sp.row(), sp2.row());
    assertArrayEquals(sp.colind(), sp2.colind());
  }

  @Test
  public void test_jacobians() {
    SX x = SX.sym("x");
    assertEquals(0, ca.jacobian(5.0, x).nnz());

    for (int i = 0; i < 5; i++) {
      jacobiansTest(Sparsity.lower(i));
      jacobiansTest(Sparsity.lower(i).T());
      jacobiansTest(Sparsity.dense(i, i));
      jacobiansTest(Sparsity.diag(i));
    }

    for (int i : new int[] {63, 64, 65, 127, 128, 129}) {
      Sparsity d = Sparsity.diag(i);
      jacobiansTest(d);
      jacobiansTest(d.plus(Sparsity.rowcol(new long[] {0}, new long[] {5}, i, i)));
      Sparsity b = Sparsity.band(i, -1).plus(Sparsity.band(i, 1));
      jacobiansTest(b.plus(Sparsity.rowcol(new long[] {0}, new long[] {5}, i, i)));
    }

    DM m = DM.ones(Sparsity.diag(129));
    m.set(new Slice(0, 50), 0, 1.0);
    m.set(new Slice(60, END), 0, 1.0);
    m.set(new Slice(6, 9), 6, 1.0);
    m.set(9, new Slice(9, 12), 1.0);

    Sparsity sp = m.get(new Slice(), new Slice(0, 120)).sparsity();
    jacobiansTest(sp);

    m = DM.ones(Sparsity.diag(64));
    m.set(new Slice(0, 50), 0, 1.0);
    m.set(new Slice(60, END), 0, 1.0);

    for (DM mm : List.of(m.T(), m))
      for (int k : new int[] {40, 20}) {
        sp = mm.get(new Slice(0, k), new Slice()).sparsity();
        jacobiansTest(sp);
        jacobiansTest(sp.T());
      }

    for (int i : new int[] {63, 64, 65, 127, 128, 129}) {
      jacobiansTest(Sparsity.lower(i));
      jacobiansTest(Sparsity.lower(i).T());
    }

    int[] sizes = RUN_SLOW ? new int[] {63, 64, 65, 127, 128, 129} : new int[] {63, 64, 65};
    for (int n : sizes)
      for (int mm : sizes) {
        sp = Sparsity.dense(n, mm);
        jacobiansTest(sp);
        DM I = DM.ones(sp);
        for (int i = 0; i < n; i++)
          for (int j = 0; j < mm; j++)
            if (random.nextDouble() < 0.5) I.set(i, j, 0.0);
        I = ca.sparsify(I);
        Sparsity spHoles = I.sparsity();
        jacobiansTest(spHoles);
        DM z = new DM(spHoles.size1(), spHoles.size2());
        List<DM> v = new ArrayList<>();
        for (int r = 0; r < 5; r++) {
          List<DM> h = new ArrayList<>(java.util.Collections.nCopies(5, z));
          h.set(r, I);
          v.add(ca.horzcat(h));
        }
        jacobiansTest(ca.vertcat(v).sparsity());
      }
  }

  private static void hessiansTest(Sparsity sp) {
    SX x = SX.sym("x", sp.size2());
    assertTrue(sp.is_equal(sp.T()));
    Function f = new Function("f", List.of(x), List.of(ca.mtimes(List.of(x.T(), new SX(DM.ones(sp)), x))));
    Function J = hessianOld(f, 0, 0);
    Sparsity sp2 = J.sparsity_out(0);
    assertArrayEquals(sp.row(), sp2.row());
    assertArrayEquals(sp.colind(), sp2.colind());
  }

  @Test
  public void test_hessians() {
    double[][][] As = {
      {{1, 1, 0, 0, 0, 0}, {1, 1, 1, 0, 1, 1}, {0, 1, 1, 1, 0, 0}, {0, 0, 1, 1, 0, 1}, {0, 1, 0, 0, 1, 0}, {0, 1, 0, 1, 0, 1}},
      {{1, 0, 0, 0, 0, 0}, {0, 1, 1, 0, 1, 1}, {0, 1, 1, 1, 0, 0}, {0, 0, 1, 1, 0, 1}, {0, 1, 0, 0, 1, 0}, {0, 1, 0, 1, 0, 1}},
      {{1, 0, 0, 0, 0, 0}, {0, 1, 0, 0, 1, 1}, {0, 0, 1, 1, 0, 0}, {0, 0, 1, 1, 0, 1}, {0, 1, 0, 0, 1, 0}, {0, 1, 0, 1, 0, 1}},
      {{0, 0, 0, 0, 0, 0}, {0, 1, 0, 0, 1, 1}, {0, 0, 1, 1, 0, 0}, {0, 0, 1, 1, 0, 1}, {0, 1, 0, 0, 1, 0}, {0, 1, 0, 1, 0, 1}},
      {{0, 0, 0, 0, 0, 0}, {0, 1, 0, 0, 1, 0}, {0, 0, 1, 1, 0, 0}, {0, 0, 1, 1, 0, 1}, {0, 1, 0, 0, 1, 0}, {0, 0, 0, 1, 0, 1}}};
    for (double[][] A : As) hessiansTest(ca.sparsify(new DM(A)).sparsity());

    for (int i : new int[] {63, 64, 65, 100, 127, 128, 129}) {
      Sparsity d = Sparsity.diag(i);
      hessiansTest(d);
      Sparsity rc = Sparsity.rowcol(new long[] {0}, new long[] {5}, i, i).plus(Sparsity.rowcol(new long[] {5}, new long[] {0}, i, i));
      hessiansTest(d.plus(rc));
      Sparsity b = Sparsity.band(i, -1).plus(Sparsity.band(i, 1));
      hessiansTest(b);
      hessiansTest(b.plus(rc));
      hessiansTest(Sparsity.dense(i, i));
      long[] zeros = new long[i], range = new long[i];
      for (int k = 0; k < i; k++) range[k] = k;
      hessiansTest(Sparsity.diag(i).plus(Sparsity.triplet(i, i, zeros, range)).plus(Sparsity.triplet(i, i, range, zeros)));

      Sparsity sp = Sparsity.dense(i, i);
      DM I = DM.ones(sp);
      for (int ii = 0; ii < i; ii++)
        for (int jj = 0; jj < i; jj++)
          if (random.nextDouble() < 0.5) {
            I.set(ii, jj, 0.0);
            I.set(jj, ii, 0.0);
          }
      I = ca.sparsify(I);
      Sparsity spHoles = I.sparsity();
      hessiansTest(spHoles);
      DM z = new DM(spHoles.size1(), spHoles.size2());
      List<DM> v = new ArrayList<>();
      for (int r = 0; r < 5; r++) {
        List<DM> h = new ArrayList<>(java.util.Collections.nCopies(5, z));
        h.set(r, I);
        v.add(ca.horzcat(h));
      }
      hessiansTest(ca.vertcat(v).sparsity());
    }
  }

  @Test
  public void test_customIO() {
    SX x = SX.sym("x");
    Function f = new Function("f", List.of(x), List.of(x.times(x), x), List.of("i0"), List.of("foo", "bar"));
    Map<String, GenericMatrix> ret = f.invoke(Map.of("i0", 12));
    checkarray(dm(144), d(ret.get("foo")));
    checkarray(dm(12), d(ret.get("bar")));
    // f_out["baz"]: python NameError on an undefined variable, nothing to port
    ret = f.invoke(Map.of("i0", new SX(12)));
    checkarray(new DM(s(ret.get("foo"))), dm(144));
    checkarray(new DM(s(ret.get("bar"))), dm(12));
    assertFalse(ret.containsKey("baz"));
  }

  @Test
  public void test_derivative_simplifications() {
    int n = 1;
    SX x = SX.sym("x", n);
    double[] r = new double[n];
    for (int i = 0; i < n; i++) r[i] = i;
    Function M = new Function("M", List.of(x), List.of(ca.mtimes(x.minus(new DM(r)), x.T())));
    MX P = MX.sym("P", n, n), X = MX.sym("X", n);
    MX M_X = m(M.invoke(X));
    Function Pf = new Function("P", List.of(X, P), List.of(ca.mtimes(M_X, P)));
    Function P_P = jacobianOld(Pf, 1, 0);
    assertFalse(P_P.toString().contains("derivative"));
  }

  @Test
  public void test_issue1464() {
    int n = 6, N = 9;
    SX x = SX.sym("x", n), u = SX.sym("u");
    Function rk4 = new Function("f", List.of(x, u), List.of(x.plus(u)));
    long[] split = {0, (N + 1) * n, (N + 1) * n + N};
    {
      List<SX> g = new ArrayList<>(), g2 = new ArrayList<>();
      SX V = SX.sym("V", (N + 1) * n + N);
      List<SX> VV = ca.vertsplit(V, split);
      List<SX> VXk = ca.vertsplit(VV.get(0), n), VUk = ca.vertsplit(VV.get(1), 1);
      for (int k = 0; k < N; k++) {
        SX xf = s(rk4.invoke(VXk.get(k), VUk.get(k)));
        List<SX> xfp = ca.vertsplit(xf, n / 2), vp = ca.vertsplit(VXk.get(k + 1), n / 2);
        g.add(xfp.get(0).minus(vp.get(0)));
        g.add(xfp.get(1).minus(vp.get(1)));
        g2.add(xf.minus(VXk.get(k + 1)));
      }
      for (int i = 0; i < 2; i++) {
        Function f = new Function("nlp", List.of(V), List.of(ca.vertcat(g)), Map.of("ad_weight_sp", i));
        assertEquals(162, f.jac_sparsity(0, 0).nnz());
        Function f2 = new Function("nlp", List.of(V), List.of(ca.vertcat(g2)), Map.of("ad_weight_sp", i));
        assertEquals(162, f2.jac_sparsity(0, 0).nnz());
      }
    }
    {
      List<MX> g = new ArrayList<>(), g2 = new ArrayList<>();
      MX V = MX.sym("V", (N + 1) * n + N);
      List<MX> VV = ca.vertsplit(V, split);
      List<MX> VXk = ca.vertsplit(VV.get(0), n), VUk = ca.vertsplit(VV.get(1), 1);
      for (int k = 0; k < N; k++) {
        MX xf = m(rk4.invoke(VXk.get(k), VUk.get(k)));
        List<MX> xfp = ca.vertsplit(xf, n / 2), vp = ca.vertsplit(VXk.get(k + 1), n / 2);
        g.add(xfp.get(0).minus(vp.get(0)));
        g.add(xfp.get(1).minus(vp.get(1)));
        g2.add(xf.minus(VXk.get(k + 1)));
      }
      for (int i = 0; i < 2; i++) {
        Function f = new Function("nlp", List.of(V), List.of(ca.vertcat(g)), Map.of("ad_weight_sp", i));
        assertEquals(162, f.jac_sparsity(0, 0).nnz());
        Function f2 = new Function("nlp", List.of(V), List.of(ca.vertcat(g2)), Map.of("ad_weight_sp", i));
        assertEquals(162, f2.jac_sparsity(0, 0).nnz());
      }
    }
  }

  @Test
  public void test_callback() {
    final class MyCallback extends Callback {
      MyCallback(String name, Map<String, ?> opts) {
        construct(name, opts);
      }

      @Override
      public List<DM> eval(List<DM> argin) {
        return List.of(argin.get(0).power(2.0));
      }
    }
    MyCallback foo = new MyCallback("my_f", Map.of());
    MX x = MX.sym("x");
    MX y = m(foo.invoke(x));
    Function f = new Function("f", List.of(x), List.of(y));
    checkarray(d(f.invoke(5)), 25);
  }

  @Test
  public void test_callback_buffer() {
    final class MyCallback extends Callback {
      MyCallback(String name) {
        construct(name, Map.of());
      }

      @Override
      public boolean has_eval_buffer() {
        return true;
      }

      @Override
      public int eval_buffer(List<java.nio.DoubleBuffer> arg, List<java.nio.DoubleBuffer> res) {
        java.nio.DoubleBuffer a = arg.get(0), r = res.get(0);
        for (int i = 0; i < a.capacity(); i++) r.put(i, a.get(i) * a.get(i));
        return 0;
      }
    }
    MyCallback foo = new MyCallback("my_f");
    MX x = MX.sym("x");
    MX y = m(foo.invoke(x));
    Function f = new Function("f", List.of(x), List.of(y));
    checkarray(d(f.invoke(5)), 25);

    final class MyCallback2 extends Callback {
      MyCallback2(String name) {
        construct(name, Map.of());
      }

      @Override
      public boolean has_eval_buffer() {
        return true;
      }

      @Override
      public long get_n_in() {
        return 3;
      }

      @Override
      public long get_n_out() {
        return 2;
      }

      @Override
      public Sparsity get_sparsity_in(long i) {
        return i == 0 ? Sparsity.dense(1, 1) : i == 1 ? Sparsity.dense(3, 1) : Sparsity.dense(3, 3);
      }

      @Override
      public Sparsity get_sparsity_out(long i) {
        return i == 0 ? Sparsity.dense(3, 1) : Sparsity.dense(3, 3);
      }

      @Override
      public int eval_buffer(List<java.nio.DoubleBuffer> arg, List<java.nio.DoubleBuffer> res) {
        // column-major, as np.frombuffer(..).reshape(order='F')
        java.nio.DoubleBuffer a = arg.get(0), b = arg.get(1), c = arg.get(2), r0 = res.get(0), r1 = res.get(1);
        for (int i = 0; i < 3; i++) {
          double v = 0;
          for (int j = 0; j < 3; j++) v += a.get(0) * c.get(i + 3 * j) * b.get(j);
          r0.put(i, v);
        }
        for (int k = 0; k < 9; k++) r1.put(k, c.get(k) * c.get(k));
        return 0;
      }
    }
    MyCallback2 foo2 = new MyCallback2("my_f");
    DM a = new DM(3), b = dm(1, 2, 3), c = dm(new double[][] {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}});
    List<?> res = (List<?>) foo2.invoke(a, b, c);
    checkarray((DM) res.get(0), ca.mtimes(a.times(c), b));
    checkarray((DM) res.get(1), ca.power(c, 2));
  }

  @Test
  public void test_callback_errors() {
    final class MyCallback extends Callback {
      MyCallback(String name) {
        construct(name, Map.of());
      }

      @Override
      public List<DM> eval(List<DM> argin) {
        throw new RuntimeException("foobar");
      }
    }
    MyCallback foo = new MyCallback("my_f");
    MX x = MX.sym("x");
    MX y = m(foo.invoke(x));
    Function f = new Function("f", List.of(x), List.of(y));
    try {
      f.invoke(3);
    } catch (RuntimeException e) {
      assertTrue(e.getMessage().contains("foobar"), e.getMessage());
    }
  }

  private static final class MySVD extends Callback {
    final long n, m, k;
    final boolean fd;

    MySVD(String name, long n, long m, boolean fd) {
      this.n = n;
      this.m = m;
      this.fd = fd;
      this.k = Math.min(n, m);
      construct(name);
    }

    @Override
    public long get_n_in() {
      return 1;
    }

    @Override
    public long get_n_out() {
      return 3;
    }

    @Override
    public Sparsity get_sparsity_in(long i) {
      return Sparsity.dense(n, m);
    }

    @Override
    public Sparsity get_sparsity_out(long i) {
      if (i == 0) return Sparsity.dense(n, k);
      else if (i == 1) return Sparsity.dense(k);
      else if (i == 2) return Sparsity.dense(k, m);
      throw new IllegalArgumentException("No such index");
    }
    // eval (numpy svd) is never evaluated by the test
  }

  private static final class Convolution extends Callback {
    final DM A;
    final long n, m;
    final boolean flagTransp;
    final List<Convolution> data = new ArrayList<>();

    Convolution(DM A, boolean flagTransp) {
      this.A = A;
      if (flagTransp) {
        m = A.size1();
        n = A.size2();
      } else {
        n = A.size1();
        m = A.size2();
      }
      this.flagTransp = flagTransp;
      construct("Convolution");
    }

    Convolution getTranspose() {
      Convolution out = new Convolution(A, !flagTransp);
      data.add(out);
      return out;
    }

    @Override
    public Sparsity get_sparsity_in(long i) {
      return Sparsity.dense(m, 1);
    }

    @Override
    public Sparsity get_sparsity_out(long i) {
      return Sparsity.dense(n, 1);
    }

    @Override
    public List<DM> eval(List<DM> arg) {
      DM x = arg.get(0);
      return List.of(flagTransp ? ca.mtimes(A.T(), x) : ca.mtimes(A, x));
    }

    @Override
    public long get_n_in() {
      return 1;
    }

    @Override
    public long get_n_out() {
      return 1;
    }

    @Override
    public boolean has_forward(long nfwd) {
      return true;
    }

    @Override
    public boolean has_reverse(long nadj) {
      return true;
    }

    @Override
    public Function get_forward(long nfwd, String name, List<String> inames, List<String> onames, Map<String, ?> opts) {
      Function mf = map(nfwd, "serial");
      MX fwdIn = MX.sym("x", m, nfwd);
      MX fwdOut = m(mf.invoke(fwdIn));
      MX arg = MX.sym("x", m, 1);
      MX dummy = MX.sym("x", new Sparsity(n, 1));
      return new Function(name, List.of(arg, dummy, fwdIn), List.of(fwdOut), inames, onames, opts);
    }

    @Override
    public Function get_reverse(long nadj, String name, List<String> inames, List<String> onames, Map<String, ?> opts) {
      Function mf = getTranspose().map(nadj, "serial");
      MX adjIn = MX.sym("x", n, nadj);
      MX adjOut = m(mf.invoke(adjIn));
      MX arg = MX.sym("x", m, 1);
      MX dummy = new MX(n, 1);
      return new Function(name, List.of(arg, dummy, adjIn), List.of(adjOut), inames, onames, opts);
    }

    @Override
    public boolean has_jacobian() {
      return false;
    }
  }

  @Test
  public void test_callback_convolution() {
    int n = 3, mm = 3;
    MySVD foo = new MySVD("foo", n, mm, true);
    DM.rng(1);
    DM.rand(n, mm);
    MX X = MX.sym("x", n, mm);
    foo.invoke(X);
    DM A = DM.rand(5, 5);
    Convolution conv = new Convolution(A, false);
    conv.forward(1);
    conv.jacobian();
  }

  /** The fun, X, Y, Z, V of the map tests */
  private static final class MapSetup {
    final Function fun;
    final List<MX> X = new ArrayList<>(), Y = new ArrayList<>(), Z = new ArrayList<>(), V = new ArrayList<>();

    MapSetup(int n) {
      SX x = SX.sym("x"), y = SX.sym("y", 2), z = SX.sym("z", 2, 2), v = SX.sym("z", Sparsity.upper(3));
      fun = new Function("f", List.of(x, y, z, v), List.of(ca.mtimes(z, y).plus(x), ca.sin(y.times(x)).T(), v.div(x)));
      for (int i = 0; i < n; i++) {
        X.add(MX.sym("x"));
        Y.add(MX.sym("y", 2));
        Z.add(MX.sym("z", 2, 2));
        V.add(MX.sym("z", Sparsity.upper(3)));
      }
    }

    List<MX> all() {
      return cat(X, Y, Z, V);
    }

    List<MX> hcats() {
      return List.of(ca.horzcat(X), ca.horzcat(Y), ca.horzcat(Z), ca.horzcat(V));
    }
  }

  private List<DM> rndLike(List<MX> xs) {
    List<DM> r = new ArrayList<>();
    for (MX x : xs) r.add(rnd(x.sparsity()));
    return r;
  }

  private static List<MX> sinAll(List<MX> xs) {
    List<MX> r = new ArrayList<>();
    for (MX x : xs) r.add(ca.sin(x));
    return r;
  }

  @Test
  public void test_mapdict() {
    SX x = SX.sym("x"), y = SX.sym("y", 2), z = SX.sym("z", 2, 2), v = SX.sym("z", Sparsity.upper(3));
    Map<String, SX> dict = new LinkedHashMap<>();
    dict.put("x", x);
    dict.put("y", y);
    dict.put("z", z);
    dict.put("v", v);
    dict.put("I", ca.mtimes(z, y).plus(x));
    dict.put("II", ca.sin(y.times(x)).T());
    dict.put("III", v.div(x));
    Function fun = new Function("f", dict, List.of("x", "y", "z", "v"), List.of("I", "II", "III"));

    MapSetup s = new MapSetup(2);
    Map<String, MX> arg = new LinkedHashMap<>();
    arg.put("x", ca.horzcat(s.X));
    arg.put("y", ca.horzcat(s.Y));
    arg.put("z", ca.horzcat(s.Z));
    arg.put("v", ca.horzcat(s.V));
    Map<String, MX> res = fun.map(2).call(arg);
    List<MX> res2 = fun.map(2).call(s.hcats());

    Function F = new Function("F", s.all(), res2);
    Function F2 = new Function("F", s.all(), List.of(res.get("I"), res.get("II"), res.get("III")));
    checkfunctionLight(F, F2, rndLike(s.all()));
  }

  private void mapNode(List<String> parallelizations, boolean light) {
    MapSetup s = new MapSetup(2);
    for (String parallelization : parallelizations) {
      System.out.println(parallelization);
      List<MX> res = s.fun.map(2, parallelization).call(s.hcats());
      Function F = new Function("F", s.all(), sinAll(res));
      List<List<MX>> resref = new ArrayList<>();
      for (int i = 0; i < s.fun.n_out(); i++) resref.add(new ArrayList<>());
      for (int k = 0; k < 2; k++) {
        List<MX> e = sinAll(s.fun.call(List.of(s.X.get(k), s.Y.get(k), s.Z.get(k), s.V.get(k))));
        for (int i = 0; i < e.size(); i++) resref.get(i).add(e.get(i));
      }
      List<MX> hs = new ArrayList<>();
      for (List<MX> r : resref) hs.add(ca.horzcat(r));
      Function Fref = new Function("F", s.all(), hs);
      List<DM> inputs = rndLike(s.all());
      for (Function f : List.of(F, F.expand("expand_" + F.name()))) {
        if (light) checkfunctionLight(f, Fref, inputs);
        else checkfunction(f, Fref, inputs);
      }
    }
  }

  @Test
  public void test_map_node() {
    mapNode(RUN_SLOW ? List.of("serial", "openmp", "unroll", "inline", "thread") : List.of("serial"), false);
  }

  @Test
  public void test_map_node_light() {
    mapNode(List.of("serial", "openmp", "unroll", "inline", "thread"), true);
  }

  @Test
  public void test_map_node_n_threads() {
    MapSetup s = new MapSetup(1);
    Function fun = s.fun;
    List<DM> X_ = new ArrayList<>(), Y_ = new ArrayList<>(), Z_ = new ArrayList<>(), V_ = new ArrayList<>();
    for (int i = 0; i < 10; i++) {
      X_.add(rnd(fun.sparsity_in(0)));
      Y_.add(rnd(fun.sparsity_in(1)));
      Z_.add(rnd(fun.sparsity_in(2)));
      V_.add(rnd(fun.sparsity_in(3)));
    }
    System.out.println(fun.map(3, "thread", 2));
    long[][] cases = {{2, 1}, {3, 1}, {3, 2}, {4, 2}, {4, 5}};
    for (long[] c : cases) {
      int k = (int) c[0];
      List<DM> inputs = List.of(ca.horzcat(X_.subList(0, k)), ca.horzcat(Y_.subList(0, k)), ca.horzcat(Z_.subList(0, k)),
          ca.horzcat(V_.subList(0, k)));
      checkfunctionLight(fun.map(k, "thread", c[1]), fun.map(k), inputs);
    }
  }

  // test_map_fenv: python-only (switches the FPU rounding mode through ctypes/libm)

  @Test
  public void test_mapsum() {
    MapSetup s = new MapSetup(2);
    Function fun = s.fun;
    for (List<MX> Z_alt : List.of(s.Z, List.of(new MX(), new MX(), new MX()))) {
      for (String parallelization : List.of("serial", "openmp", "unroll", "thread")) {
        List<MX> res = fun.mapsum(List.of(ca.horzcat(s.X), ca.horzcat(s.Y), ca.horzcat(Z_alt), ca.horzcat(s.V)), parallelization);
        for (int ad_weight_sp : new int[] {0, 1}) {
          Function F = new Function("F", s.all(), sinAll(res), Map.of("ad_weight", 0, "ad_weight_sp", ad_weight_sp));
          List<MX> resref = new ArrayList<>();
          for (int i = 0; i < fun.n_out(); i++) resref.add(new MX(0.0));
          for (int k = 0; k < 2; k++) {
            List<MX> e = fun.call(List.of(s.X.get(k), s.Y.get(k), Z_alt.get(k), s.V.get(k)));
            for (int i = 0; i < e.size(); i++) resref.set(i, resref.get(i).plus(e.get(i)));
          }
          Function Fref = new Function("F", s.all(), sinAll(resref));
          List<DM> inputs = rndLike(s.all());
          if (!parallelization.equals("thread")) checkCodegen(F, inputs);
          for (Function f : List.of(F, toSXFun(F))) checkfunction(f, Fref, inputs);
        }
      }
    }
  }

  @Test
  public void test_mapsum2() {
    MapSetup s = new MapSetup(2);
    Function fun = s.fun;
    int n = 2;
    MX Z = MX.sym("z", 2, 2), V = MX.sym("z", Sparsity.upper(3));
    for (String parallelization : List.of("true_map_sum", "serial", "openmp", "unroll", "thread"))
      for (int ad_weight_sp : new int[] {0, 1})
        for (int ad_weight : new int[] {0, 1}) {
          Map<String, Object> opts = Map.of("ad_weight_sp", ad_weight_sp, "ad_weight", ad_weight);
          Function F;
          if (parallelization.equals("true_map_sum"))
            F = fun.map(n, new boolean[] {false, false, true, true}, new boolean[] {true, false, false}, opts);
          else F = fun.map("map", parallelization, n, new long[] {2, 3}, new long[] {0}, opts);

          MX acc = new MX(0.0);
          List<MX> bl = new ArrayList<>(), cl = new ArrayList<>();
          for (int k = 0; k < n; k++) {
            List<MX> abc = l(fun.invoke(s.X.get(k), s.Y.get(k), Z, V));
            acc = acc.plus(abc.get(0));
            bl.add(abc.get(1));
            cl.add(abc.get(2));
          }
          Function Fref = new Function("F", List.of(ca.horzcat(s.X), ca.horzcat(s.Y), Z, V),
              List.of(acc, ca.horzcat(bl), ca.horzcat(cl)));
          List<DM> inputs = List.of(ca.horzcat(rndLike(s.X)), ca.horzcat(rndLike(s.Y)), rnd(Z.sparsity()), rnd(V.sparsity()));
          if (!parallelization.equals("thread")) checkCodegen(F, inputs);
          for (Function f : List.of(F, toSXFun(F))) checkfunction(f, Fref, inputs);
          checkSerialize(F, inputs);
        }
  }

  @Test
  public void test_issue4274() {
    MX x = MX.sym("x");
    Function c = new Function("x", List.of(x), List.of(x.printme(4).power(2.0)), List.of("x"), List.of("out"));
    Map<String, Double> localInputs = Map.of("x", 5.0);
    Function F = c.map(2, "thread");
    Map<String, DM> wideInputs = new LinkedHashMap<>();
    for (Map.Entry<String, Double> e : localInputs.entrySet())
      wideInputs.put(e.getKey(), ca.horzcat(List.of(new DM(e.getValue()), new DM(e.getValue() * 1.2))));
    System.out.println(wideInputs);
    F.invoke(wideInputs);
  }

  @Test
  public void test_repmatnode() {
    MX x = MX.sym("x", 2);
    MX y = ca.sin(ca.repmat(x.power(2.0), 1, 3));
    MX z = MX.sym("y", 2, 2);
    Function F = new Function("f", List.of(x, z), List.of(ca.sum2(ca.sum1(y))));
    SX xs = SX.sym("x", 2);
    SX ys = ca.sin(ca.repmat(xs.power(2.0), 1, 3));
    SX zs = SX.sym("y", 2, 2);
    Function Fref = new Function("f", List.of(xs, zs), List.of(ca.sum2(ca.sum1(ys))));
    List<DM> inputs = List.of(dm(1, 7), dm(new double[][] {{3, 0}, {2, 4}}));
    checkCodegen(F, inputs);
    checkfunction(F, Fref, inputs);
  }

  @Test
  public void test_repsumnode() {
    MX x = MX.sym("x", 2), z = MX.sym("y", 2, 2);
    Function F = new Function("f", List.of(x, z),
        List.of(ca.sin(ca.repsum(x.power(2.0).T(), 1, 2)), ca.cos(x.power(2.0)).times(2.0).times(x).T()));
    SX xs = SX.sym("x", 2), zs = SX.sym("y", 2, 2);
    Function Fref = new Function("f", List.of(xs, zs),
        List.of(ca.sin(ca.repsum(xs.power(2.0).T(), 1, 2)), ca.cos(xs.power(2.0)).times(2.0).times(xs).T()));
    List<DM> inputs = List.of(dm(1, 7), dm(new double[][] {{3, 0}, {2, 4}}));
    checkCodegen(F, inputs);
    checkfunction(F, Fref, inputs);
  }

  @Test
  public void test_unknown_options() {
    SX x = SX.sym("x");
    // python spells the constructor SXFunction (undefined, so trivially raises); the intent is ported
    assertThrows(RuntimeException.class, () -> new Function("f", List.of(x), List.of(x), Map.of("fooo", false)));
    assertThrows(RuntimeException.class, () -> new Function("f", List.of(x), List.of(x), Map.of("ad_weight", "foo")));
  }

  @Test
  public void test_unknown_options_stringvector() {
    Assumptions.assumeTrue(hasNlpsol("ipopt"));
    SX x = SX.sym("x");
    ca.nlpsol("mysolver", "ipopt", Map.of("x", x, "f", x.power(2.0)), Map.of("monitor", List.of("nlp_f")));
    String err = captureStderr(() -> ca.nlpsol("mysolver", "ipopt", Map.of("x", x, "f", x.power(2.0)), Map.of("monitor", List.of("abc"))));
    assertTrue(err.contains("Ignoring monitor 'abc'. Available functions: nlp_f"), err);
  }

  @Test
  public void test_mapaccum() {
    SX x = SX.sym("x", 2), y = SX.sym("y"), z = SX.sym("z", 2, 2), v = SX.sym("v", Sparsity.upper(3));
    Function fun = new Function("f", List.of(x, y, z, v), List.of(ca.mtimes(z, x).plus(y), ca.sin(y.times(x)).T(), v.div(y)));
    int n = 2;
    MX X = MX.sym("x", x.sparsity());
    List<MX> Y = new ArrayList<>(), Z = new ArrayList<>(), V = new ArrayList<>();
    for (int i = 0; i < n; i++) {
      Y.add(MX.sym("y", y.sparsity()));
      Z.add(MX.sym("z", z.sparsity()));
      V.add(MX.sym("v", v.sparsity()));
    }
    DM X_ = rnd(x.sparsity());
    List<DM> Y_ = rndLike(Y), Z_ = rndLike(Z), V_ = rndLike(V);

    for (int ad_weight = 0; ad_weight < 2; ad_weight++)
      for (int ad_weight_sp = 0; ad_weight_sp < 2; ad_weight_sp++) {
        Function F = fun.mapaccum("map", n, new long[] {0}, new long[] {0}, Map.of("ad_weight_sp", ad_weight_sp, "ad_weight", ad_weight));
        F.forward(2);
        MX XP = X;
        List<MX> Y0s = new ArrayList<>(), Y1s = new ArrayList<>(), Xps = new ArrayList<>();
        for (int k = 0; k < n; k++) {
          List<MX> r = l(fun.invoke(XP, Y.get(k), Z.get(k), V.get(k)));
          XP = r.get(0);
          Y0s.add(r.get(1));
          Y1s.add(r.get(2));
          Xps.add(XP);
        }
        Function Fref = new Function("f", List.of(X, ca.horzcat(Y), ca.horzcat(Z), ca.horzcat(V)),
            List.of(ca.horzcat(Xps), ca.horzcat(Y0s), ca.horzcat(Y1s)));
        List<DM> inputs = List.of(X_, ca.horzcat(Y_), ca.horzcat(Z_), ca.horzcat(V_));
        for (Function f : List.of(F, toSXFun(F))) {
          checkfunction(f, Fref, inputs);
          checkCodegen(f, inputs);
        }
      }

    fun = new Function("f", List.of(y, x, z, v),
        List.of(ca.mtimes(z, x).plus(y).plus(ca.trace(v).power(2.0)), ca.sin(y.times(x)).T(), v.div(y)));
    for (int ad_weight = 0; ad_weight < 2; ad_weight++)
      for (int ad_weight_sp = 0; ad_weight_sp < 2; ad_weight_sp++) {
        Function F = fun.mapaccum("map", n, new long[] {1, 3}, new long[] {0, 2}, Map.of("ad_weight_sp", ad_weight_sp, "ad_weight", ad_weight));
        MX XP = X, VP = V.get(0);
        List<MX> Y0s = new ArrayList<>(), Xps = new ArrayList<>(), Vps = new ArrayList<>();
        for (int k = 0; k < n; k++) {
          List<MX> r = l(fun.invoke(Y.get(k), XP, Z.get(k), VP));
          XP = r.get(0);
          VP = r.get(2);
          Y0s.add(r.get(1));
          Xps.add(XP);
          Vps.add(VP);
        }
        Function Fref = new Function("f", List.of(ca.horzcat(Y), X, ca.horzcat(Z), V.get(0)),
            List.of(ca.horzcat(Xps), ca.horzcat(Y0s), ca.horzcat(Vps)));
        List<DM> inputs = List.of(ca.horzcat(Y_), X_, ca.horzcat(Z_), V_.get(0));
        for (Function f : List.of(F, toSXFun(F))) {
          checkfunction(f, Fref, inputs);
          checkCodegen(f, inputs);
        }
      }
  }

  @Test
  public void test_mapaccum_schemes() {
    SX x = SX.sym("x", 2), y = SX.sym("y"), z = SX.sym("z", 2, 2), v = SX.sym("v", Sparsity.upper(3));
    Function fun = new Function("f", List.of(y, z, x, v), List.of(ca.mtimes(z, x).plus(y), ca.sin(y.times(x)).T(), v.div(y)),
        List.of("y", "z", "x", "v"), List.of("out0", "out1", "out2"));
    int n = 2;
    Function F = fun.mapaccum("map", n, new long[] {2}, new long[] {0});
    List<String> schemeInFun = fun.name_in(), schemeOutFun = fun.name_out();
    List<String> schemeInF = F.name_in(), schemeOutF = F.name_out();
    for (int rep = 0; rep < 2; rep++) {
      // python: assertTrue(len(a), len(b)) only checks the first is nonzero
      assertTrue(schemeInFun.size() > 0);
      assertTrue(schemeOutFun.size() > 0);
      for (int i = 0; i < Math.min(schemeInFun.size(), schemeInF.size()); i++) assertEquals(schemeInFun.get(i), schemeInF.get(i));
      for (int i = 0; i < Math.min(schemeOutFun.size(), schemeOutF.size()); i++) assertEquals(schemeOutFun.get(i), schemeOutF.get(i));
      fun = new Function("f", List.of(x, y, z, v), List.of(ca.mtimes(z, x).plus(y), ca.sin(y.times(x)).T(), v.div(y)),
          List.of("x", "y", "z", "v"), List.of("out0", "out1", "out2"));
      F = fun.mapaccum("map", n);
    }
  }

  @Test
  public void test_depends_on() {
    SX x = SX.sym("x");
    SX y = x.power(2.0);
    String s = "";
    try {
      ca.depends_on(x, y);
    } catch (RuntimeException e) {
      s = e.getMessage();
    }
    assertTrue(s.contains("not symbolic"), s);
    try {
      new Function("f", List.of(y), List.of(x));
    } catch (RuntimeException e) {
      s = e.getMessage();
    }
    assertTrue(s.contains("not symbolic"), s);
  }

  private static boolean same(Object a, double b) {
    return Math.abs(d(a).to_double() - b) < 1e-8;
  }

  private void checkCodegenSerialize(Function f, List<DM> inputs) {
    Function f2 = checkCodegen(f, inputs);
    if (f2 != null) checkSerialize(f2, inputs);
  }

  @Test
  public void test_1d_interpolant() {
    double[][] grid = {{0, 1, 1.5, 2, 3}};
    double[] values = {0, 1, 2, 5, 3};
    Function F = ca.interpolant("F", "linear", grid, values);
    double[][] pairs = {{3.4, 3 - 0.4 * 2}, {2.4, 5 - 0.4 * 2}, {1.6, 2 + 3 * 0.1 / 0.5}, {1.4, 1 + 0.4 / 0.5}, {0.4, 0.4}, {-.6, -0.6}};
    MX X = MX.sym("x");
    new Function("F", List.of(X), List.of(m(F.invoke(X))));
    for (double[] p : pairs) {
      assertTrue(same(F.invoke(p[0]), p[1]));
      checkCodegenSerialize(F, ins(p[0]));
      checkSerialize(F, ins(p[0]));
    }
    System.out.println(F.stats());

    X = MX.sym("x");
    Function J = new Function("F", List.of(X), List.of(ca.jacobian(m(F.invoke(X)), X)));
    double[][] jpairs = {{3.4, -2}, {2.4, -2}, {1.6, 6}, {1.4, 2}, {0.4, 1}, {-.6, 1}, {1, 2}, {0.99, 1}};
    for (double[] p : jpairs) {
      assertTrue(same(J.invoke(p[0]), p[1]));
      checkCodegen(J, ins(p[0]));
      checkSerialize(J, ins(p[0]));
    }
  }

  @Test
  public void test_2d_interpolant() {
    double[][] grid = {{0, 1, 4, 5}, {0, 2, 3}};
    double[] values = {0, 1, 8, 3, 10, -11, 12, 13, 20, 31, -42, 53};
    Function F = ca.interpolant("F", "linear", grid, values);
    double a0 = -11 + 0.4 * (31 + 11), a1 = 12 + 0.4 * (-42 - 12);
    double[][] pairs = {{1, 2, -11}, {1, 3, 31}, {4, 2, 12}, {4, 3, -42}, {1, 2.4, a0}, {4, 2.4, a1},
      {3, 2, -11 + 2.0 / 3 * (12 + 11)}, {3, 3, 31 + 2.0 / 3 * (-42 - 31)}, {3, 2.4, a0 + 2.0 / 3 * (a1 - a0)}};
    for (double[] p : pairs) {
      DM a = dm(p[0], p[1]);
      checkarray(d(F.invoke(a)), p[2]);
      checkCodegen(F, List.of(a));
    }

    MX X = MX.sym("x", 2);
    Function J = new Function("J", List.of(X), List.of(ca.jacobian(m(F.invoke(X)), X)));
    double jx0 = (12 + 11) / 3.0, jx1 = (-42 - 31) / 3.0, jx2 = (13 - 12), jx3 = (53 + 42);
    double jy0 = 31 + 11, jy1 = -42 - 12;
    double[][] jpairs = {{1, 2, jx0, jy0}, {1, 3, jx1, jy0}, {4, 2, jx2, jy1}, {4, 3, jx3, jy1},
      {1, 2.4, jx0 + (jx1 - jx0) * 0.4, 31 + 11}, {4, 2.4, jx2 + (jx3 - jx2) * 0.4, -42 - 12},
      {3, 2, jx0, jy0 + (jy1 - jy0) * 2.0 / 3}, {3, 3, jx1, jy0 + (jy1 - jy0) * 2.0 / 3},
      {3, 2.4, jx0 + (jx1 - jx0) * 0.4, jy0 + (jy1 - jy0) * 2.0 / 3}};
    for (double[] p : jpairs) {
      DM a = dm(p[0], p[1]);
      checkarray(d(J.invoke(a)).T(), dm(p[2], p[3]));
      checkCodegen(J, List.of(a));
      checkSerialize(J, List.of(a));
    }
  }

  @Test
  public void test_1d_interpolant_uniform() {
    double[][] grid = {{0, 1, 2}};
    double[] values = {0, 1, 2};
    for (String mode : List.of("linear", "exact", "binary")) {
      Function F = ca.interpolant("F", "linear", grid, values, Map.of("lookup_mode", List.of(mode)));
      assertTrue(same(F.invoke(2.4), 2.4));
      assertTrue(same(F.invoke(1.4), 1.4));
      assertTrue(same(F.invoke(0), 0));
      assertTrue(same(F.invoke(1), 1));
      assertTrue(same(F.invoke(2), 2));
      assertTrue(same(F.invoke(6), 6));
      assertTrue(same(F.invoke(0.4), 0.4));
      assertTrue(same(F.invoke(-.6), -.6));
      ca.interpolant("F", "linear", new double[][] {linspace(0, 1, 7)}, new double[] {0, 1, 2, 3, 4, 5, 6}, Map.of("lookup_mode", List.of("exact")));
    }
    grid = new double[][] {{2, 4, 6}};
    values = new double[] {10, 7, 1};
    for (String mode : List.of("linear", "exact", "binary")) {
      Function F = ca.interpolant("F", "linear", grid, values, Map.of("lookup_mode", List.of(mode)));
      assertTrue(same(F.invoke(1), 11.5));
      assertTrue(same(F.invoke(2), 10));
      assertTrue(same(F.invoke(3), 8.5));
      assertTrue(same(F.invoke(4), 7));
      assertTrue(same(F.invoke(5), 4));
      assertTrue(same(F.invoke(6), 1));
      assertTrue(same(F.invoke(7), -2));
      ca.interpolant("F", "linear", new double[][] {linspace(0, 1, 7)}, new double[] {0, 1, 2, 3, 4, 5, 6}, Map.of("lookup_mode", List.of("exact")));
    }
  }

  @Test
  public void test_2d_interpolant_uniform() {
    double[][] grid = {{0, 1, 2}, {0, 1, 2}};
    double[] values = {0, 1, 2, 10, 11, 12, 20, 21, 22};
    for (String mode : List.of("linear", "exact", "binary")) {
      Function F = ca.interpolant("F", "linear", grid, values, Map.of("lookup_mode", List.of(mode, mode)));
      double[][] cases = {{2.4, 0.5, 7.4}, {1.4, 0.5, 6.4}, {0.4, 0.5, 5.4}, {1, 0.5, 6}, {1, 0, 1}, {0, 0, 0}, {0.4, 1, 10.4},
        {-.6, 0.5, 4.4}, {-.6, 1.5, 14.4}, {-.6, 2.5, 24.4}, {-.6, 3.5, 34.4}};
      for (double[] c : cases) assertTrue(same(F.invoke(new double[] {c[0], c[1]}), c[2]));
    }
  }

  /** Multilinear interpolation on a rectilinear grid, data in column-major (first dimension fastest) */
  private static double interpn(double[][] grid, double[] d, double[] x) {
    int n = grid.length;
    int[] lo = new int[n];
    double[] w = new double[n];
    for (int i = 0; i < n; i++) {
      double[] g = grid[i];
      int k = 0;
      while (k < g.length - 2 && x[i] >= g[k + 1]) k++;
      lo[i] = k;
      w[i] = (x[i] - g[k]) / (g[k + 1] - g[k]);
    }
    double r = 0;
    for (int c = 0; c < (1 << n); c++) {
      double wc = 1;
      int idx = 0, stride = 1;
      for (int i = 0; i < n; i++) {
        int b = (c >> i) & 1;
        wc *= b == 1 ? w[i] : 1 - w[i];
        idx += (lo[i] + b) * stride;
        stride *= grid[i].length;
      }
      r += wc * d[idx];
    }
    return r;
  }

  @Test
  public void test_nd_linear() {
    int N = 4;
    for (int n_dim = 1; n_dim < 5; n_dim++) {
      double[][] grid = new double[n_dim][];
      double[] x0 = new double[n_dim];
      for (int i = 0; i < n_dim; i++) {
        double[] g = new double[N];
        for (int k = 0; k < N; k++) g[k] = random.nextDouble();
        Arrays.sort(g);
        grid[i] = g;
        x0[i] = Arrays.stream(g).average().getAsDouble();
      }
      double[] d = new double[(int) Math.pow(N, n_dim)];
      for (int k = 0; k < d.length; k++) d[k] = random.nextDouble();
      double xref = interpn(grid, d, x0);
      Function F = ca.interpolant("F", "linear", grid, d);
      checkarray(d(F.invoke(x0)), xref);
    }
  }

  @Test
  public void test_2d_bspline() {
    double[][] d_knots = {linspace(0, 1, 5), linspace(0, 1, 6)};
    double[] d_flat = rnd(5 * 6, 1).nonzeros();
    Function LUT = ca.interpolant("name", "bspline", d_knots, d_flat);
    jacobianOld(LUT, 0, 0);
    hessianOld(LUT, 0, 0);
    checkCodegen(LUT, List.of(dm(0.2, 0.3)));
    checkSerialize(LUT, List.of(dm(0.2, 0.3)));
    // comparison against scipy.interpolate.RectBivariateSpline: python-only (scipy reference)
  }

  @Test
  public void test_1d_bspline_call_fun() {
    double[][] d_knots = {linspace(0, 1, 5)};
    double[] d_flat = rnd(5, 1).nonzeros();
    Function LUT = ca.interpolant("name", "bspline", d_knots, d_flat);
    SX xs = SX.sym("x");
    Function f0 = new Function("f", List.of(xs), List.of(s(LUT.invoke(ca.sin(xs))).times(2.0)));
    MX xm = MX.sym("x");
    Function f1 = new Function("f", List.of(xm), List.of(m(LUT.invoke(ca.sin(xm))).times(2.0)));
    checkfunction(f0, f1, ins(0.32));
  }

  @Test
  public void test_1d_bspline() {
    double[][] d_knots = {linspace(0, 1, 5)};
    double[] d_flat = rnd(5, 1).nonzeros();
    Function LUT = ca.interpolant("name", "bspline", d_knots, d_flat);
    checkCodegen(LUT, ins(0.2));
    jacobianOld(LUT, 0, 0);
    hessianOld(LUT, 0, 0);
    // comparison against scipy.interpolate.InterpolatedUnivariateSpline: python-only (scipy reference)
  }

  @Test
  public void test_bspline_even_degree() {
    double[] x = linspace(0, 3, 12), y = new double[12];
    for (int i = 0; i < 12; i++) y[i] = Math.sin(x[i]);
    Function LUT = ca.interpolant("LUT", "bspline", new double[][] {x}, y, Map.of("degree", new long[] {2}));
    // comparison against scipy make_interp_spline: python-only (scipy reference)
    checkCodegen(LUT, ins(0.2));
    for (long k : new long[] {4, 6}) {
      LUT = ca.interpolant("LUT", "bspline", new double[][] {x}, y, Map.of("degree", new long[] {k}));
      checkarray(d(LUT.invoke(new DM(x).T())), new DM(y).T(), 1e-10, "");
      checkarray(d(LUT.invoke(1.234)), new DM(Math.sin(1.234)), 1e-4, "");
    }
    double[] xg = linspace(0, 1, 6), yg = linspace(0, 2, 7);
    double[] V = new double[6 * 7];
    for (int j = 0; j < 7; j++)
      for (int i = 0; i < 6; i++) V[i + 6 * j] = Math.cos(xg[i]) * yg[j];
    LUT = ca.interpolant("LUT", "bspline", new double[][] {xg, yg}, V, Map.of("degree", new long[] {2, 3}));
    for (int[] ij : new int[][] {{0, 0}, {2, 3}, {5, 6}, {4, 1}})
      checkarray(d(LUT.invoke(dm(xg[ij[0]], yg[ij[1]]))), new DM(V[ij[0] + 6 * ij[1]]), 1e-10, "");
  }

  /** sin(x+3*y) */
  private static class SinFun extends Callback {
    @Override
    public long get_n_in() {
      return 2;
    }

    @Override
    public long get_n_out() {
      return 1;
    }

    @Override
    public List<DM> eval(List<DM> arg) {
      DM x = arg.get(0), y = arg.get(1);
      DM z0 = y.times(3.0);
      DM z1 = x.plus(z0);
      return List.of(ca.sin(z1));
    }
  }

  @Test
  public void test_Callback_Jacobian() {
    MX x = MX.sym("x"), y = MX.sym("y");
    List<DM> numInputs = ins(0.2, 0.7);
    Function g = new Function("g", List.of(x, y), List.of(ca.sin(x.plus(y.times(3.0)))));
    final class Fun extends SinFun {
      Fun() {
        construct("Fun", Map.of());
      }

      @Override
      public boolean has_forward(long nfwd) {
        return false;
      }

      @Override
      public boolean has_reverse(long nadj) {
        return false;
      }

      @Override
      public boolean has_jacobian() {
        return true;
      }

      @Override
      public Function get_jacobian(String name, List<String> inames, List<String> onames, Map<String, ?> opts) {
        SX x = SX.sym("x"), y = SX.sym("y");
        SX out_g = SX.sym("out_g", new Sparsity(1, 1));
        SX c = ca.cos(x.plus(y.times(3.0)));
        return new Function(name, List.of(x, y, out_g), List.of(c, c.times(3.0)), inames, onames, opts);
      }
    }
    Fun f = new Fun();
    checkfunction(f, g, numInputs);
  }

  // test_Callback_errors: python-only (overrides get_sparsity_in with a wrong signature returning an int)

  @Test
  public void test_Callback_unstringable_exception() {
    final class Unstringable extends RuntimeException {
      private static final long serialVersionUID = 1L;

      @Override
      public String getMessage() {
        throw new RuntimeException("boom in getMessage");
      }

      @Override
      public String toString() {
        throw new RuntimeException("boom in toString");
      }
    }
    final class Fun extends Callback {
      Fun() {
        construct("Fun", Map.of());
      }

      @Override
      public List<DM> eval(List<DM> arg) {
        throw new Unstringable();
      }
    }
    Fun f = new Fun();
    assertThrows(RuntimeException.class, () -> f.invoke(0));
  }

  @Test
  public void test_Callback() {
    MX x = MX.sym("x"), y = MX.sym("y");
    List<DM> numInputs = ins(0.2, 0.7);
    Function g = new Function("g", List.of(x, y), List.of(ca.sin(x.plus(y.times(3.0)))));
    final class Fun extends SinFun {
      Fun() {
        construct("Fun", Map.of());
      }
    }
    for (boolean indirect : new boolean[] {true, false}) {
      Fun cb = new Fun();
      Function f = indirect ? new Function("f", List.of(x, y), List.of(m(cb.invoke(x, y)))) : cb;
      checkfunctionLight(f, g, numInputs);
      // f.gradient(): python AttributeError (Function has no gradient method), nothing to port
      assertThrows(RuntimeException.class, () -> jacobianOld(f, 0, 0));
      assertThrows(RuntimeException.class, () -> f.forward(1));
      assertThrows(RuntimeException.class, () -> f.reverse(1));
    }
  }

  @Test
  public void test_Callback_dimcheck() {
    final class Fun extends Callback {
      final List<DM> ret;

      Fun(List<DM> ret) {
        this.ret = ret;
        construct("Fun");
      }

      @Override
      public long get_n_in() {
        return 2;
      }

      @Override
      public long get_n_out() {
        return 1;
      }

      @Override
      public List<DM> eval(List<DM> arg) {
        return ret;
      }
    }
    Fun f = new Fun(List.of(new DM(2), new DM(1)));
    String s = "";
    try {
      f.invoke(2);
    } catch (RuntimeException e) {
      s = e.getMessage();
    }
    assertTrue(s.contains("Incorrect number of inputs"), s);
    f = new Fun(List.of(new DM(2), new DM(1)));
    s = "";
    try {
      f.invoke(2, 3);
    } catch (RuntimeException e) {
      s = e.getMessage();
    }
    assertTrue(s.contains("Expected 1 output"), s);
    s = "";
    f = new Fun(List.of(DM.zeros(2, 2)));
    try {
      f.invoke(2, 3);
    } catch (RuntimeException e) {
      s = e.getMessage();
    }
    assertTrue(s.contains("Shape mismatch"), s);
  }

  private static final class ForwardFun extends Callback {
    ForwardFun(String name) {
      construct(name, Map.of("verbose", true));
    }

    @Override
    public long get_n_in() {
      return 2 + 1 + 2;
    }

    @Override
    public long get_n_out() {
      return 1;
    }

    @Override
    public List<DM> eval(List<DM> arg) {
      DM x = arg.get(0), y = arg.get(1);
      DM z1 = x.plus(y.times(3.0));
      List<DM> ret = new ArrayList<>();
      for (int i = 3; i < arg.size(); i += 2) {
        DM dx = arg.get(i), dy = arg.get(i + 1);
        DM dz1 = dx.plus(dy.times(3.0));
        ret.add(ca.cos(z1).times(dz1));
      }
      return ret;
    }
  }

  private static final class BackwardFun extends Callback {
    BackwardFun(String name) {
      construct(name, Map.of("verbose", true));
    }

    @Override
    public long get_n_in() {
      return 2 + 1 + 1;
    }

    @Override
    public long get_n_out() {
      return 2;
    }

    @Override
    public List<DM> eval(List<DM> arg) {
      DM x = arg.get(0), y = arg.get(1);
      DM z1 = x.plus(y.times(3.0));
      List<DM> ret = new ArrayList<>();
      for (int i = 3; i < arg.size(); i++) {
        DM bz2 = arg.get(i);
        DM bz1 = bz2.times(ca.cos(z1));
        DM bx = bz1, bz0 = bz1;
        DM by = bz0.times(3.0);
        ret.add(bx);
        ret.add(by);
      }
      return ret;
    }
  }

  private static final class SensFun extends SinFun {
    final boolean hasFwd, hasAdj;
    ForwardFun cbFwd;
    BackwardFun cbRev;

    SensFun(boolean hasFwd, boolean hasAdj, Map<String, ?> opts) {
      this.hasFwd = hasFwd;
      this.hasAdj = hasAdj;
      construct("Fun", opts);
    }

    @Override
    public boolean has_forward(long nfwd) {
      return hasFwd && nfwd == 1;
    }

    @Override
    public Function get_forward(long nfwd, String name, List<String> inames, List<String> onames, Map<String, ?> opts) {
      assertEquals(1, nfwd);
      cbFwd = new ForwardFun(name);
      return cbFwd;
    }

    @Override
    public boolean has_reverse(long nadj) {
      return hasAdj && nadj == 1;
    }

    @Override
    public Function get_reverse(long nadj, String name, List<String> inames, List<String> onames, Map<String, ?> opts) {
      assertEquals(1, nadj);
      cbRev = new BackwardFun(name);
      return cbRev;
    }
  }

  @Test
  public void test_Callback_sens() {
    MX x = MX.sym("x"), y = MX.sym("y");
    List<DM> numInputs = ins(0.2, 0.7);
    Function g = new Function("g", List.of(x, y), List.of(ca.sin(x.plus(y.times(3.0)))));
    List<SensFun> keep = new ArrayList<>();
    for (boolean indirect : new boolean[] {true, false}) {
      boolean[][] variants = {{true, true}, {true, false}, {false, true}, {true, false}, {false, true}};
      for (boolean[] v : variants) {
        SensFun cb = new SensFun(v[0], v[1], Map.of("verbose", true));
        keep.add(cb);
        Function f = indirect ? new Function("f", List.of(x, y), List.of(m(cb.invoke(x, y)))) : cb;
        checkfunction(f, g, numInputs);
      }
    }
  }

  @Test
  public void test_callback_jacobian_sparsity() {
    MX x = MX.sym("x", 2);
    double h = 1e-7;
    for (boolean withJacSparsity : new boolean[] {true, false}) {
      List<DM> calls = new ArrayList<>();
      final class Fun extends Callback {
        Fun() {
          construct("Fun", Map.of("enable_fd", true, "fd_method", "forward", "fd_options", Map.of("h", h, "h_iter", false)));
        }

        @Override
        public long get_n_in() {
          return 1;
        }

        @Override
        public long get_n_out() {
          return 1;
        }

        @Override
        public Sparsity get_sparsity_in(long i) {
          return Sparsity.dense(2, 1);
        }

        @Override
        public Sparsity get_sparsity_out(long i) {
          return Sparsity.dense(2, 1);
        }

        @Override
        public List<DM> eval(List<DM> arg) {
          DM x = arg.get(0);
          calls.add(x.minus(dm(5, 7)).div(h));
          return List.of(x.power(2.0));
        }

        @Override
        public boolean has_forward(long nfwd) {
          return false;
        }

        @Override
        public boolean has_reverse(long nadj) {
          return false;
        }

        @Override
        public boolean has_jac_sparsity(long oind, long iind) {
          return withJacSparsity;
        }

        @Override
        public Sparsity get_jac_sparsity(long oind, long iind, boolean symmetric) {
          assertEquals(0, oind);
          assertEquals(0, iind);
          return Sparsity.diag(2);
        }
      }
      Fun f = new Fun();
      MX y = ca.jacobian(m(f.invoke(x)), x);
      Function F = new Function("F", List.of(x), List.of(y));
      DM[] J = new DM[1];
      captureStdout(() -> J[0] = d(F.invoke(new double[] {5, 7})));
      DM callsM = ca.horzcat(calls);
      DM J_ref = dm(new double[][] {{5 * 2, 0}, {0, 7 * 2}});
      if (withJacSparsity) {
        checkarray(callsM, dm(new double[][] {{0, 0}, {1, 1}}).T(), 1e-5, "");
        J_ref = ca.sparsify(J_ref);
      } else {
        checkarray(callsM, dm(new double[][] {{0, 0}, {1, 0}, {0, 1}}).T(), 1e-5, "");
      }
      checkarray(J[0], J_ref, 1e-5, "");
      assertTrue(J[0].sparsity().is_equal(J_ref.sparsity()));
    }
  }

  @Test
  public void test_common_specific_options() {
    Assumptions.assumeTrue(hasNlpsol("ipopt"));
    SX x = SX.sym("x");
    Map<String, SX> nlp = Map.of("x", x, "f", x.power(4.0));
    Map<String, Double> x0 = Map.of("x0", 0.1);

    String out = captureStdout(() -> ca.nlpsol("solver", "ipopt", nlp));
    assertFalse(out.contains("nlp_f"));
    Function solver = ca.nlpsol("solver", "ipopt", nlp, Map.of("common_options", Map.of("final_options", Map.of("print_in", true))));
    out = captureStdout(() -> solver.invoke(x0));
    System.out.println(out);
    assertTrue(out.contains("Function nlp_f"));
    assertTrue(out.contains("Function nlp_grad_f"));
    Function solver2 = ca.nlpsol("solver", "ipopt", nlp, Map.of("specific_options", Map.of("nlp_f", Map.of("final_options", Map.of("print_in", true)))));
    out = captureStdout(() -> solver2.invoke(x0));
    assertTrue(out.contains("Function nlp_f"));
    assertFalse(out.contains("Function nlp_grad_f"));
    Function solver3 = ca.nlpsol("solver", "ipopt", nlp, Map.of("common_options", Map.of("final_options", Map.of("print_in", true)),
        "specific_options", Map.of("nlp_f", Map.of("final_options", Map.of("print_in", false)))));
    out = captureStdout(() -> solver3.invoke(x0));
    assertFalse(out.contains("Function nlp_f"));
    assertTrue(out.contains("Function nlp_grad_f"));
    Function solver4 = ca.nlpsol("solver", "ipopt", nlp, Map.of("common_options", Map.of("final_options", Map.of("print_in", false)),
        "specific_options", Map.of("nlp_f", Map.of("final_options", Map.of("print_in", true)))));
    out = captureStdout(() -> solver4.invoke(x0));
    assertTrue(out.contains("Function nlp_f"));
    assertFalse(out.contains("Function nlp_grad_f"));

    String err = captureStderr(() -> ca.nlpsol("solver", "ipopt", nlp));
    assertEquals(0, err.length(), err);
    err = captureStderr(() -> ca.nlpsol("solver", "ipopt", nlp, Map.of("specific_options", Map.of("nlp_foo", Map.of("helper_options", Map.of("verbose", true))))));
    assertTrue(("Ignoring" + err).length() > 0);
    assertTrue(err.contains("nlp_g"), err);
    assertThrows(RuntimeException.class, () -> ca.nlpsol("solver", "ipopt", nlp, Map.of("specific_options", Map.of("nlp_foo", 3))));
  }

  @Test
  public void test_oracle_options() {
    Assumptions.assumeTrue(hasNlpsol("ipopt"));
    DM.set_precision(16);
    int N = 5;
    MX x = MX.sym("x", N, 1);
    Map<String, Object> options = new LinkedHashMap<>();
    options.put("ipopt", Map.of("hessian_approximation", "limited-memory"));
    options.put("oracle_options", Map.of("enable_fd", true, "enable_forward", false, "enable_reverse", false, "print_in", true, "fd_method", "central"));
    Function solver = ca.nlpsol("solver", "ipopt", Map.of("x", x, "f", ca.dot(x, x)), options);
    String out = captureStdout(() -> solver.invoke(Map.of("x0", 1)));
    DM.set_precision(6);
    assertTrue(out.contains("[1, 1, 1, 1, 1]"));
  }

  @Test
  public void test_expm() {
    Assumptions.assumeTrue(ca.has_expm("slicot"));
    Assumptions.assumeTrue(hasIntegrator("cvodes"));
    int n = 4;
    DM Anum = rnd(n, n);
    MX A = MX.sym("A", n, n), t = MX.sym("t");
    Function fr = new Function("fr", List.of(A, t), List.of(expmIntegrator(A.times(t), n)));
    Function f = new Function("f", List.of(A, t), List.of(ca.expm(A.times(t))));
    checkfunction(fr, f, ins(Anum, 1.1), 1e-8);
    fr = new Function("fr", List.of(t), List.of(expmIntegrator(new MX(Anum).times(t), n)));
    f = new Function("f", List.of(t), List.of(ca.expm_const(Anum, t)));
    checkfunction(fr, f, ins(1.1), 1e-8);

    MX JA = ca.jacobian(ca.expm(A.times(t)), A), Jt = ca.jacobian(ca.expm(A.times(t)), t);
    assertEquals(n * n * n * n, JA.nnz());
    assertEquals(n * n, Jt.nnz());
    JA = ca.jacobian(ca.expm_const(A, t), A);
    Jt = ca.jacobian(ca.expm_const(A, t), t);
    assertEquals(0, JA.nnz());
    assertEquals(n * n, Jt.nnz());
  }

  private static MX expmIntegrator(MX Ain, int n) {
    MX x = MX.sym("x", n), As = MX.sym("A", n, n);
    Map<String, MX> dae = Map.of("x", x, "p", ca.vec(As), "ode", ca.mtimes(As, x));
    Function intg = ca.integrator("intg", "cvodes", dae, Map.of("reltol", 1e-14, "abstol", 1e-14));
    Function Intg = intg.map("identity", "serial", n, new long[] {1}, new long[] {});
    Map<String, MX> out = Intg.call(Map.of("x0", new MX(DM.eye(n)), "p", ca.vec(As)));
    Function expmF = new Function("expm", List.of(As), List.of(out.get("xf")));
    return m(expmF.invoke(Ain));
  }

  @Test
  public void test_conditional() {
    MX x = MX.sym("x", 2, 2);
    MX sp1 = MX.sym("y", Sparsity.lower(2)), sp2 = MX.sym("z", Sparsity.diag(2));
    Function f1 = new Function("f", List.of(sp2, x), List.of(x.power(2.0), x.times(sp2)));
    Function f2 = new Function("f", List.of(sp1, x), List.of(x.power(2.0).times(2.0), ca.sin(sp1)));
    Function f3 = new Function("f", List.of(sp1, sp2), List.of(sp1.times(sp2), sp1.plus(sp2)));
    Function F = Function.conditional("test", List.of(f1, f2), f3);
    Function Fsx = F.expand();
    double[][] A = rnd(2, 2).full(), B = rnd(2, 2).full();
    for (int i = -1; i < 3; i++) {
      checkfunction(F, Fsx, ins(i, A, B));
      checkCodegen(F, ins(i, A, B));
    }
  }

  // test_max_num_dir (first definition): shadowed in python by the second definition of the same name

  @Test
  public void test_max_num_dir() {
    MX x = MX.sym("x");
    Function f = new Function("ffff", List.of(x), List.of(ca.sin(x)));
    Function fm = f.mapaccum("mapaccum", 100, 1, Map.of("base", 4));
    CodeGenerator c = new CodeGenerator("me");
    c.add(fm);
    String code = c.dump();
    assertTrue(code.contains("ffff_acc4_acc4_acc4"));
  }

  @Test
  public void test_codegen_with_jac_sparsity() {
    Assumptions.assumeTrue(RUN_SLOW);
    MX x = MX.sym("x", 3), y = MX.sym("y", 3, 3);
    Function f = new Function("f", List.of(x, y), List.of(x.power(3.0), ca.mtimes(y, x)));
    CodeGenerator c = new CodeGenerator("me");
    c.add(f, true);
    Function fJ = f.jacobian();
    DM.rng(0);
    DM x0 = DM.rand(x.sparsity()), y0 = DM.rand(y.sparsity());
    Function F2 = checkCodegen(f, List.of(x0, y0), Map.of(), Map.of("enable_fd", true), true, "c89", false);
    Function FJ = F2.jacobian();
    for (int i = 0; i < fJ.n_out(); i++) assertTrue(fJ.sparsity_out(i).is_equal(FJ.sparsity_out(i)));
    F2 = checkCodegen(f, List.of(x0, y0), Map.of(), Map.of("enable_fd", true), false, "c89", false);
    FJ = F2.jacobian();
    for (int i = 0; i < fJ.n_out(); i++) assertTrue(FJ.sparsity_out(i).is_dense());
  }

  /** column-major flat data of two random 5x6 tables, interleaved as in the python tests */
  private double[][] twoTables() {
    double[] d0 = rnd(30, 1).nonzeros(), d1 = rnd(30, 1).nonzeros(), dflat = new double[60];
    for (int i = 0; i < 30; i++) {
      dflat[2 * i] = d0[i];
      dflat[2 * i + 1] = d1[i];
    }
    return new double[][] {d0, d1, dflat};
  }

  private void multiout(String solver, boolean all) {
    double[][] d_knots = {linspace(0, 1, 5), linspace(0, 1, 6)};
    double[][] t = twoTables();
    Function LUT0 = ca.interpolant("name", solver, d_knots, t[0]);
    Function LUT1 = ca.interpolant("name", solver, d_knots, t[1]);
    double[] d_flat = t[2];
    Function LUTf = ca.interpolant("name", solver, d_knots, d_flat);
    MX x = MX.sym("x"), y = MX.sym("y");
    MX xy = ca.vertcat(List.of(x, y));
    Function LUT_sep = new Function("f", List.of(x, y), List.of(ca.vertcat(List.of(m(LUT0.invoke(xy)), m(LUT1.invoke(xy))))));
    Function LUT = new Function("f", List.of(x, y), List.of(m(LUTf.invoke(xy))));
    List<DM> inputs = ins(0.2, 0.333);
    checkfunction(LUT, LUT_sep, inputs);
    checkCodegen(LUT, inputs);
    checkSerialize(LUT, inputs);

    Function LUT_param = ca.interpolant("name", solver, d_knots, 2);
    Function f = new Function("LUTp", List.of(x, y), List.of(m(LUT_param.invoke(xy, new DM(d_flat)))));
    checkfunction(LUT, f, inputs);
    checkCodegen(f, inputs);
    checkSerialize(f, inputs);
    if (!all) return;

    double[] cat = new double[11];
    System.arraycopy(d_knots[0], 0, cat, 0, 5);
    System.arraycopy(d_knots[1], 0, cat, 5, 6);
    DM d_knots_cat = new DM(cat);

    LUT_param = ca.interpolant("name", solver, new long[] {5, 6}, 2);
    f = new Function("LUTp", List.of(x, y), List.of(m(LUT_param.invoke(xy, d_knots_cat, new DM(d_flat)))));
    checkfunction(LUT, f, inputs);
    checkCodegen(f, inputs);
    checkSerialize(f, inputs);

    LUT_param = ca.interpolant("name", solver, new long[] {5, 6}, d_flat);
    f = new Function("LUTp", List.of(x, y), List.of(m(LUT_param.invoke(xy, d_knots_cat))));
    checkfunction(LUT, f, inputs);
    checkCodegen(f, inputs);
    checkSerialize(f, inputs);
  }

  @Test
  public void test_2d_linear_multiout() {
    multiout("linear", true);
  }

  @Test
  public void test_2d_bspline_multiout() {
    multiout("bspline", false);
  }

  /** B-spline sum_i c_i B_{i,p}(x), right-continuous at knots except at the right end (as scipy BSpline) */
  private static double bsplineRef(double[] t, double[] c, int p, double x) {
    int n = c.length;
    int mu = p;
    while (mu < n - 1 && x >= t[mu + 1]) mu++;
    double[] dd = new double[p + 1];
    for (int j = 0; j <= p; j++) dd[j] = c[mu - p + j];
    for (int r = 1; r <= p; r++)
      for (int j = p; j >= r; j--) {
        int i = mu - p + j;
        double a = (x - t[i]) / (t[i + p - r + 1] - t[i]);
        dd[j] = (1 - a) * dd[j - 1] + a * dd[j];
      }
    return dd[p];
  }

  @Test
  public void test_bspline_repeated_knots() {
    Object[][] cases = {{3, new double[] {0, 0, 0, 0, 1, 1, 2, 3, 3, 3, 4, 4, 4, 4}}, {2, new double[] {0, 0, 0, 1, 1, 2, 2, 3, 3, 3}}};
    for (Object[] cs : cases) {
      int degree = (Integer) cs[0];
      double[] knots = (double[]) cs[1];
      int n = knots.length - degree - 1;
      double[] c = rnd(n, 1).nonzeros();
      MX x = MX.sym("x"), C = MX.sym("C", n);
      Function f = new Function("f", List.of(x), List.of(ca.bspline(x, new DM(c), new double[][] {knots}, new long[] {degree}, 1)));
      Function F = new Function("F", List.of(x, C), List.of(ca.bspline(x, C, new double[][] {knots}, new long[] {degree}, 1)));
      for (double xe = 0; xe < knots[knots.length - 1] + 0.25; xe += 0.5) {
        double ref = bsplineRef(knots, c, degree, xe);
        checkarray(d(f.invoke(xe)), new DM(ref), 1e-12, "x=" + xe);
        checkarray(d(F.invoke(xe, c)), new DM(ref), 1e-12, "x=" + xe);
      }
      checkCodegen(f, ins(2));
      checkCodegen(F, ins(3, c));
    }
  }

  @Test
  public void test_parametric_bspline() {
    double[][] knots = {{0, 0, 0, 0, 0.2, 0.5, 0.8, 1, 1, 1, 1}, {0, 0, 0, 0.1, 0.5, 0.9, 1, 1, 1}};
    MX x = MX.sym("x", 2);
    DM data = rnd(7 * 6 * 2, 1);
    MX y = ca.bspline(x, data, knots, new long[] {3, 2}, 2);
    MX C = MX.sym("C", data.size1(), 1);
    MX Y = ca.bspline(x, C, knots, new long[] {3, 2}, 2);
    Function f = new Function("f", List.of(x), List.of(y));
    Function F = new Function("f", List.of(x, C), List.of(Y));
    F = new Function("f", List.of(x), List.of(m(F.invoke(x, data))));
    List<DM> inputs = List.of(dm(0.3, 0.4));
    checkfunction(f, F, inputs);
    checkCodegen(F, inputs);
    checkSerialize(F, inputs);

    List<MX> K = List.of(MX.sym("K0", knots[0].length), MX.sym("K1", knots[1].length));
    MX Y_pk = ca.bspline(x, C, K, new long[] {3, 2}, 2);
    List<MX> knotsDm = List.of(new MX(new DM(knots[0])), new MX(new DM(knots[1])));
    MX Y_sub = ca.substitute(List.of(Y_pk), K, knotsDm).get(0);
    Function F_pk = new Function("f_pk", List.of(x, C), List.of(Y_sub));
    F_pk = new Function("f_pk", List.of(x), List.of(m(F_pk.invoke(x, data))));
    checkfunction(f, F_pk, inputs);
    checkCodegen(F_pk, inputs, Map.of(), Map.of(), false, "c99", false);
  }

  /** np.random.seed(0) data of test_smooth_linear, which pins a reference value */
  private static final double[] SMOOTH_LINEAR_DATA = {0.5488135039273248, 0.26455561210462697, 0.4375872112626925,
    0.6120957227224214, 0.5680445610939323, 0.6976311959272649, 0.7781567509498505, 0.31542835092418386, 0.11827442586893322,
    0.2088767560948347, 0.7151893663724195, 0.7742336894342167, 0.8917730007820798, 0.6169339968747569, 0.925596638292661,
    0.06022547162926983, 0.8700121482468192, 0.3637107709426226, 0.6399210213275238, 0.16130951788499626, 0.6027633760716439,
    0.45615033221654855, 0.9636627605010293, 0.9437480785146242, 0.07103605819788694, 0.6667667154456677, 0.978618342232764,
    0.5701967704178796, 0.1433532874090464, 0.6531083254653984, 0.5448831829968969, 0.5684339488686485, 0.3834415188257777,
    0.6818202991034834, 0.08712929970154071, 0.6706378696181594, 0.7991585642167236, 0.43860151346232035, 0.9446689170495839,
    0.2532916025397821, 0.4236547993389047, 0.018789800436355142, 0.7917250380826646, 0.359507900573786, 0.02021839744032572,
    0.2103825610738409, 0.46147936225293185, 0.9883738380592262, 0.5218483217500717, 0.4663107728563063, 0.6458941130666561,
    0.6176354970758771, 0.5288949197529045, 0.43703195379934145, 0.832619845547938, 0.1289262976548533, 0.7805291762864555,
    0.10204481074802807, 0.4146619399905236, 0.24442559200160274};

  @Test
  public void test_smooth_linear() {
    double[][] d_knots = {linspace(0, 1, 5), linspace(0, 1, 6)};
    double[] d_flat = SMOOTH_LINEAR_DATA;
    Function LUT_linear = ca.interpolant("name", "linear", d_knots, d_flat);
    Function LUT = ca.interpolant("name", "bspline", d_knots, d_flat, Map.of("algorithm", "smooth_linear", "smooth_linear_frac", 0.1));
    Function LUT2 = ca.interpolant("name", "bspline", d_knots, d_flat, Map.of("algorithm", "smooth_linear", "smooth_linear_frac", 0.05));

    double[] p1 = {0.2, 0.333}, p2 = {0.25, 0.4};
    double diff1 = ca.norm_1(d(LUT.invoke(p1)).minus(d(LUT_linear.invoke(p1)))).to_double();
    double diff2 = ca.norm_2(d(LUT.invoke(p2)).minus(d(LUT_linear.invoke(p2)))).to_double();
    double diff22 = ca.norm_2(d(LUT2.invoke(p2)).minus(d(LUT_linear.invoke(p2)))).to_double();
    assertTrue(diff1 <= 1e-10);
    assertTrue(diff2 >= 100 * diff1);
    assertTrue(diff22 <= diff2 * 0.51);
    assertTrue(diff22 >= diff2 * 0.49);
    checkarray(d(LUT2.invoke(new double[] {0.26, 0.39})), dm(0.9261362392504342, 0.9157189108791507));

    MX xy = MX.sym("xy", 2);
    Function LUT_param = ca.interpolant("name", "bspline", d_knots, 2, Map.of("algorithm", "smooth_linear", "smooth_linear_frac", 0.1));
    Function f = new Function("LUTp", List.of(xy), List.of(m(LUT_param.invoke(xy, new DM(d_flat)))));
    List<DM> inputs = List.of(dm(0.2, 0.333));
    checkfunction(LUT, f, inputs);
    checkCodegen(f, inputs, Map.of(), Map.of(), false, "c89", true);
    checkSerialize(f, inputs);
  }

  @Test
  public void test_codegen_avoid_stack() {
    SX x = SX.sym("x", 3, 3);
    Function f = new Function("f", List.of(x), List.of(ca.det(x)));
    checkCodegen(f, List.of(rnd(3, 3)));
    checkCodegen(f, List.of(rnd(3, 3)), Map.of("avoid_stack", true), Map.of(), false, "c89", false);
  }

  @Test
  public void test_serialize() {
    for (Map<String, Object> opts : List.<Map<String, Object>>of(Map.of("debug", true), Map.of())) {
      SX x = SX.sym("x");
      SX z = ca.sin(x.plus(3.0));
      Function f = new Function("f", List.of(x), List.of(z));
      Function fs = Function.deserialize(f.serialize(opts));
      checkfunction(f, fs, ins(2));

      x = SX.sym("x");
      z = ca.sin(x.plus(3.0));
      f = new Function("f", List.of(x), List.of(z, new SX(Double.NaN), new SX(Double.NEGATIVE_INFINITY), new SX(Double.POSITIVE_INFINITY)));
      fs = Function.deserialize(f.serialize(opts));
      checkfunction(f, fs, ins(2));

      x = SX.sym("x");
      SX y = SX.sym("y", Sparsity.lower(3));
      z = x.plus(y);
      SX z1 = ca.sparsify(ca.vertcat(List.of(z.get(0), new SX(0), z.get(1))));
      SX z2 = z.T();
      f = new Function("f", List.of(x, y), List.of(z1, z2, x.power(2.0)), List.of("x", "y"), List.of("a", "b", "c"));
      fs = Function.deserialize(f.serialize(opts));
      assertEquals("x", fs.name_in(0));
      assertEquals("a", fs.name_out(0));
      assertEquals("f", fs.name());
      List<DM> in = ins(3.7, new double[][] {{1, 0, 0}, {2, 3, 0}, {4, 5, 6}});
      checkfunction(f, fs, in);
      fs = Function.deserialize(f.serialize(opts));
      checkfunction(f, fs, in);

      x = SX.sym("x");
      SX p = SX.sym("p");
      new Function("f", List.of(x), List.of(p), Map.of("allow_free", true));
      // pickle round trip of the SXFunction with free parameters: python-only (pickle)

      MX xm = MX.sym("x");
      f = new Function("f", List.of(xm), List.of(xm.power(2.0)));
      fs = Function.deserialize(f.serialize(opts));
      checkfunction(f, fs, ins(3.7));

      xm = MX.sym("x");
      MX ym = MX.sym("y", 2);
      MX w = ca.if_else(xm, ca.atan2(ca.norm_fro(ym).times(3.0).times(ym), xm), xm.minus(ym), true);
      MX zm = ca.sin(xm.times(2.0)).times(w.get(0)).plus(1.0);
      Function g = new Function("g", List.of(xm, ym), List.of(w.minus(xm)));
      Function gmap = g.map(2, "thread", 2);
      Function gmapsx = gmap.expand();
      MX q = m(gmap.invoke(ca.horzcat(List.of(xm.times(2.0), xm.minus(ym.get(1)))), ca.horzcat(List.of(zm.plus(ym), ca.cos(zm.plus(ym))))))
          .plus(ca.rdivide(1.0, m(gmapsx.invoke(ca.horzcat(List.of(xm.times(2.0), xm.minus(ym.get(1)))), ca.repmat(zm.plus(ym), 1, 2)))));
      q = ca.solve(q, ym.times(2.0), "qr");
      q = q.plus(ca.bilin(dm(new double[][] {{1, 3}, {7, 8}}), q, q.times(2.0)));
      f = new Function("f", List.of(xm, ym), List.of(q.plus(1.0), ca.jacobian(q, ca.vertcat(List.of(xm, ym)))));
      fs = Function.deserialize(f.serialize(opts));
      checkfunction(f, fs, List.of(new DM(1.1), dm(2.7, 3)));
    }
  }

  @Test
  public void test_serialize_recursion_limit() {
    SX xs = SX.sym("x");
    SX ys = xs;
    for (int i = 0; i < 10000; i++) ys = ca.sin(ys);
    Function.deserialize(new Function("foo", List.of(xs), List.of(ys)).serialize());
    MX xm = MX.sym("x");
    MX ym = xm;
    for (int i = 0; i < 10000; i++) ym = ca.sin(ym);
    Function.deserialize(new Function("foo", List.of(xm), List.of(ym)).serialize());
  }

  @Test
  public void test_string() {
    MX x = MX.sym("x"), y = MX.sym("y");
    Function f = new Function("f", List.of(x), List.<MX>of(), List.of("x"), List.of());
    assertTrue(f.toString().contains("(x)->()"));
    f = new Function("f", List.<MX>of(), List.of(x), List.of(), List.of("y"), Map.of("allow_free", true));
    assertTrue(f.toString().contains("()->(y)"));
    f = new Function("f", List.of(x), List.of(x.power(2.0)), List.of("x"), List.of("y"));
    assertTrue(f.toString().contains("(x)->(y)"));
    f = new Function("f", List.of(x, y), List.of(x.power(2.0), x.times(y)), List.of("x", "y"), List.of("w", "z"));
    assertTrue(f.toString().contains("(x,y)->(w,z)"));
  }

  @Test
  public void test_fold() {
    SX x = SX.sym("x", 2, 2);
    Function f = new Function("f", List.of(x), List.of(ca.sin(x)));
    Function F = f.fold(10);
    SX x0 = x;
    for (int i = 0; i < 10; i++) x = s(f.invoke(x));
    Function Fref = new Function("f", List.of(x0), List.of(x));
    checkfunction(F, Fref, List.of(dm(new double[][] {{1, 2}, {3, 7}})));
  }

  @Test
  public void test_thread_safety() {
    MX x = MX.sym("x"), y = MX.sym("y");
    Function f = new Function("f", List.of(x, y), List.of(x.power(2.0).minus(y)));
    double[] ones = new double[200];
    Arrays.fill(ones, 1);
    double[] lin = linspace(0, 10, 200);
    for (String rf : List.of("newton", "fast_newton", "bisection")) {
      Function finv = ca.rootfinder("finv", rf, f);
      Function finvPar = finv.map(50, "unroll").map(4, "thread", 4);
      DM res = d(finvPar.invoke(ones, lin));
      checkarray(ca.norm_inf(res.T().minus(ca.sqrt(new DM(lin)))), new DM(0), 1e-5, "");
      finvPar = finv.map(200, "thread", 4);
      res = d(finvPar.invoke(ones, lin));
      checkarray(ca.norm_inf(res.T().minus(ca.sqrt(new DM(lin)))), new DM(0), 1e-5, "");
    }
  }

  @Test
  public void test_threadsafe_symbolics_weakref() throws Exception {
    int nthreads = 24, ntrials = 8000;
    Function[] box = new Function[1];
    org.casadi.WeakRef[] refs = new org.casadi.WeakRef[nthreads];
    boolean[] stop = {false};
    java.util.concurrent.CyclicBarrier start = new java.util.concurrent.CyclicBarrier(nthreads + 1);
    java.util.concurrent.CyclicBarrier done = new java.util.concurrent.CyclicBarrier(nthreads + 1);
    List<Thread> threads = new ArrayList<>();
    for (int t = 0; t < nthreads; t++) {
      final int i = t;
      threads.add(new Thread(() -> {
        try {
          while (true) {
            start.await();
            if (stop[0]) return;
            refs[i] = new org.casadi.WeakRef(box[0]);
            done.await();
          }
        } catch (InterruptedException | java.util.concurrent.BrokenBarrierException e) {
          return;
        }
      }));
    }
    for (Thread t : threads) t.start();
    int orphans = 0;
    try {
      for (int trial = 0; trial < ntrials; trial++) {
        MX x = MX.sym("x", 2);
        box[0] = new Function("f", List.of(x), List.of(ca.sin(x)));
        start.await();
        done.await();
        box[0].delete();
        box[0] = null;
        for (org.casadi.WeakRef r : refs)
          if (r.alive()) orphans++;
        for (org.casadi.WeakRef r : refs) r.delete();
      }
    } finally {
      stop[0] = true;
      for (Thread t : threads) t.interrupt();
      for (Thread t : threads) t.join();
    }
    assertEquals(0, orphans, "weak() double-create race: " + orphans + " orphaned weak references survived their owner's destruction");
  }

  @Test
  public void test_mapped_eval() {
    SX x = SX.sym("x"), y = SX.sym("y", 2);
    Function f = new Function("f", List.of(x, y), List.of(ca.sin(x).times(y)), List.of("x", "y"), List.of("r"));
    DM r1 = d(f.invoke(1, dm(3, 4))), r2 = d(f.invoke(2, dm(3, 4))), r3 = d(f.invoke(3, dm(3, 4)));
    DM r_all = d(f.invoke(dm(new double[][] {{1, 2, 3}}), dm(3, 4)));
    checkarray(r_all, ca.horzcat(List.of(r1, r2, r3)), "Mapped evaluation (DM)");
    MX z = MX.sym("z", 1, 3);
    MX rz = m(f.invoke(z, dm(3, 4)));
    Function F = new Function("F", List.of(z), List.of(rz));
    DM r_mx = d(F.invoke(dm(new double[][] {{1, 2, 3}})));
    checkarray(r_all, r_mx, "Mapped evaluation (MX)");
  }

  @Test
  public void test_default_arg() {
    MX x = MX.sym("x"), y = MX.sym("y"), z = MX.sym("z");
    Function f = new Function("f", List.of(x, y, z), List.of(x, y, z), List.of("x", "y", "z"), List.of("a", "b", "c"),
        Map.of("default_in", new double[] {1, 2, 3}));
    checkCodegen(f, Map.of("x", new DM(5), "z", new DM(3)), true);
  }

  @Test
  public void test_factory_inherit_options() {
    MX x = MX.sym("x", 5);
    for (String op : List.of("grad:f:x", "jac:f:x", "triu:hess:f:x:x")) {
      Function f = new Function("f", List.of(x), List.of(ca.dot(x, x)), List.of("x"), List.of("f"), Map.of("verbose", true));
      String out = captureStdout(() -> f.factory("fgrad", List.of("x"), List.of(op)));
      assertTrue(out.contains("::init"));

      Function f2 = new Function("f", List.of(x), List.of(ca.dot(x, x)), List.of("x"), List.of("f"),
          Map.of("enable_fd", true, "enable_forward", false, "enable_reverse", false, "print_in", true, "fd_method", "central",
              "fd_options", Map.of("h", 1e-7, "h_iter", false)));
      Function fgrad = f2.factory("fgrad", List.of("x"), List.of(op));
      out = captureStdout(() -> fgrad.invoke(0));
      assertTrue(out.contains("-1e-07,") || out.contains("-1e-007,"), out);
      assertTrue(out.contains("1e-07,") || out.contains("1e-007,"), out);
    }
  }

  private static boolean hasImporter(String name) {
    return org.casadi.Importer.has_plugin(name);
  }

  @Test
  public void test_jit_codegen_options() throws IOException {
    Assumptions.assumeTrue(hasImporter("shell"));
    Assumptions.assumeTrue(RUN_SLOW);
    MX x = MX.sym("x");
    String dir = Files.createTempDirectory("casadi_java_jit").toString();
    Map<String, Object> opts = new LinkedHashMap<>();
    opts.put("jit", true);
    opts.put("compiler", "shell");
    opts.put("jit_cleanup", false);
    opts.put("jit_temp_suffix", false);
    opts.put("jit_name", "jit_codegen_options");
    opts.put("jit_options", Map.of("directory", dir));
    opts.put("codegen_options", Map.of("prefix", "cgopt"));
    Function f = new Function("f", List.of(x), List.of(x.power(2.0)), opts);
    for (Function g : List.of(f, Function.deserialize(f.serialize()))) {
      checkarray(d(g.invoke(3)), 9);
      String src = Files.readString(new File(dir, "jit_codegen_options.c").toPath()).replace("\r\n", "\n");
      assertTrue(src.contains("cgopt_"));
    }
  }

  private static java.util.Set<String> compiled(String out) {
    java.util.Set<String> found = new java.util.HashSet<>();
    java.util.regex.Matcher mt = java.util.regex.Pattern.compile("Compiling function '(\\w+?)'").matcher(out);
    while (mt.find()) found.add(mt.group(1));
    return found;
  }

  @Test
  public void test_inherit_jit_options() {
    Assumptions.assumeTrue(hasNlpsol("ipopt"));
    Assumptions.assumeTrue(hasImporter("shell"));
    Assumptions.assumeTrue(RUN_SLOW);
    MX x = MX.sym("x");
    Function f = new Function("Q", List.of(x), List.of(x.power(2.0)), Map.of("jit", true, "compiler", "shell", "verbose", true));
    MX g = f.call(List.of(ca.sin(x)), false, true).get(0);
    MX gc = ca.cos(g);
    Function[] solver = new Function[1];
    String out = captureStdout(() -> solver[0] = ca.nlpsol("solver", "ipopt", Map.of("x", x, "f", x.power(2.0), "g", gc)));
    java.util.Set<String> found = compiled(out);
    assertEquals(1, found.size(), found.toString());
    assertTrue(found.contains("fwd1_Q"));

    out = captureStdout(() -> checkSerialize(solver[0], Map.of("x0", new DM(2))));
    found = compiled(out);
    assertEquals(2, found.size(), found.toString());
    assertTrue(found.contains("Q"));
    assertTrue(found.contains("fwd1_Q"));
  }

  @Test
  public void test_codegen_external_name_collision() {
    Assumptions.assumeTrue(hasImporter("shell"));
    Assumptions.assumeTrue(RUN_SLOW);
    MX x = MX.sym("x");
    Function f = new Function("f", List.of(x), List.of(x.plus(1.0)));
    f.generate("external_name_collision.c");
    org.casadi.Importer compiler = new org.casadi.Importer("external_name_collision.c", "shell");
    Function ext = ca.external("f", compiler);
    checkarray(d(ext.invoke(2)), 3);
    Function wrapper = new Function("wrapper", List.of(x), List.of(m(ext.invoke(x))));
    assertInException("name conflicts with an external dependency", () -> ext.generate("external_name_collision_wrapped.c"));
    for (List<Function> functions : List.of(List.of(f, wrapper), List.of(wrapper, f))) {
      CodeGenerator cg = new CodeGenerator("external_name_collision_module.c");
      cg.add(functions.get(0));
      assertInException("name conflicts with an external dependency", () -> cg.add(functions.get(1)));
    }
    CodeGenerator cg = new CodeGenerator("external_name_collision_valid.c");
    cg.add(wrapper);
    cg.add(new Function("other", List.of(x), List.of(x.times(2.0))));
    assertTrue(cg.dump().contains("wrapper"));
  }

  @Test
  public void test_generate_dependencies_jit() {
    Assumptions.assumeTrue(hasNlpsol("ipopt"));
    Assumptions.assumeTrue(hasImporter("shell"));
    Assumptions.assumeTrue(RUN_SLOW);
    MX x = MX.sym("x", 2);
    Function solver = ca.nlpsol("solver", "ipopt",
        Map.of("x", x, "f", ca.sumsqr(x), "g", x.get(0).plus(x.get(1)).minus(10.0)),
        Map.of("jit", true, "compiler", "shell", "print_time", false, "ipopt.print_level", 0));
    solver.generate_dependencies("generate_dependencies_jit.c");
    org.casadi.Importer compiler = new org.casadi.Importer("generate_dependencies_jit.c", "shell");
    for (String name : List.of("nlp_f", "nlp_g"))
      checkfunctionLight(solver.get_function(name), ca.external(name, compiler), List.of(dm(1, 2), new DM(0)));
  }
}

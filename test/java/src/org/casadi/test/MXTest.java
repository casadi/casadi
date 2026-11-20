package org.casadi.test;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertNotEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.io.File;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.function.BinaryOperator;
import java.util.function.UnaryOperator;

import org.casadi.CasadiException;
import org.casadi.CasadiMeta;
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
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;

/** Port of test/python/mx.py (first part, up to test_project_simplify) */
public class MXTest extends CasadiTestCase {

  private static final boolean KNOWN_BUGS = System.getenv("CASADI_JAVA_KNOWN_BUGS") != null;

  /** FunctionPool entry: casadi operation on MX, reference on DM */
  private static final class Op1 {
    final String name;
    final boolean nozero;
    final UnaryOperator<MX> cas;
    final UnaryOperator<DM> ref;

    Op1(String name, boolean nozero, UnaryOperator<MX> cas, UnaryOperator<DM> ref) {
      this.name = name;
      this.nozero = nozero;
      this.cas = cas;
      this.ref = ref;
    }
  }

  private static final class Op2 {
    final String name;
    final BinaryOperator<MX> cas;
    final BinaryOperator<DM> ref;

    Op2(String name, BinaryOperator<MX> cas, BinaryOperator<DM> ref) {
      this.name = name;
      this.cas = cas;
      this.ref = ref;
    }
  }

  private List<Op1> pool, Jpool;
  private List<Op2> matrixbinarypool;

  @BeforeEach
  public void setUp() {
    pool = new ArrayList<>();
    pool.add(new Op1("sqrt", false, ca::sqrt, ca::sqrt));
    pool.add(new Op1("sin", false, ca::sin, ca::sin));
    pool.add(new Op1("cos", false, ca::cos, ca::cos));
    pool.add(new Op1("tan", false, ca::tan, ca::tan));
    pool.add(new Op1("arctan", false, ca::atan, ca::atan));
    pool.add(new Op1("arcsin", false, ca::asin, ca::asin));
    pool.add(new Op1("arccos", false, ca::acos, ca::acos));
    pool.add(new Op1("exp", false, ca::exp, ca::exp));
    pool.add(new Op1("log", true, ca::log, ca::log));
    pool.add(new Op1("x^0", true, x -> x.power(0.0), x -> x.power(0.0)));
    pool.add(new Op1("^1", false, x -> x.power(1.0), x -> x.power(1.0)));
    pool.add(new Op1("^-2", true, x -> x.power(-2.0), x -> x.power(-2.0)));
    pool.add(new Op1("^0.3", false, x -> x.power(0.3), x -> x.power(0.3)));
    pool.add(new Op1("floor", false, ca::floor, ca::floor));
    pool.add(new Op1("ceil", false, ca::ceil, ca::ceil));
    Jpool = new ArrayList<>();
    Jpool.add(new Op1("sqrt", false, ca::sqrt, x -> ca.diag(ca.rdivide(1.0, ca.sqrt(x).times(2.0)))));
    Jpool.add(new Op1("sin", false, ca::sin, x -> ca.diag(ca.cos(x))));
    Jpool.add(new Op1("cos", false, ca::cos, x -> ca.diag(ca.sin(x).unaryMinus())));
    Jpool.add(new Op1("tan", false, ca::tan, x -> ca.diag(ca.rdivide(1.0, ca.cos(x).power(2.0)))));
    Jpool.add(new Op1("arctan", false, ca::atan, x -> ca.diag(ca.rdivide(1.0, x.power(2.0).plus(1.0)))));
    Jpool.add(new Op1("arcsin", false, ca::asin, x -> ca.diag(ca.rdivide(1.0, ca.sqrt(ca.minus(1.0, x.power(2.0)))))));
    Jpool.add(new Op1("arccos", false, ca::acos, x -> ca.diag(ca.rdivide(-1.0, ca.sqrt(ca.minus(1.0, x.power(2.0)))))));
    Jpool.add(new Op1("exp", false, ca::exp, x -> ca.diag(ca.exp(x))));
    Jpool.add(new Op1("log", false, ca::log, x -> ca.diag(ca.rdivide(1.0, x))));
    Jpool.add(new Op1("x^0", false, x -> x.power(0.0), x -> ca.diag(DM.zeros(x.size1(), x.size2()))));
    Jpool.add(new Op1("^1", false, x -> x.power(1.0), x -> ca.diag(DM.ones(x.size1(), x.size2()))));
    Jpool.add(new Op1("^-2", false, x -> x.power(-2.0), x -> ca.diag(ca.rdivide(-2.0, x.power(3.0)))));
    Jpool.add(new Op1("^0.3", false, x -> x.power(0.3), x -> ca.diag(ca.rdivide(0.3, x.power(0.7)))));
    matrixbinarypool = new ArrayList<>();
    matrixbinarypool.add(new Op2("Matrix+Matrix", ca::plus, ca::plus));
    matrixbinarypool.add(new Op2("Matrix-Matrix", ca::minus, ca::minus));
    matrixbinarypool.add(new Op2("Matrix*Matrix", ca::times, ca::times));
    matrixbinarypool.add(new Op2("fmin", ca::fmax, ca::fmax));
    matrixbinarypool.add(new Op2("fmax", ca::fmin, ca::fmin));
    matrixbinarypool.add(new Op2("mtimes(Matrix,Matrix.T)", (a, b) -> ca.mtimes(a, b.T()), (a, b) -> ca.mtimes(a, b.T())));
    matrixbinarypool.add(new Op2("arctan2", ca::atan2, ca::atan2));
    matrixbinarypool.add(new Op2("mtimes(Matrix,Matrix.T)", (a, b) -> ca.mtimes(a, b.T()), (a, b) -> ca.mtimes(a, b.T())));
  }

  // ---------------------------------------------------------------- helpers

  private static DM dmsp(Sparsity sp, double... v) {
    return new DM(sp, new DM(v));
  }

  private static DM idx(double... i) {
    return new DM(i);
  }

  private static double[] range(int n) {
    double[] r = new double[n];
    for (int i = 0; i < n; i++) r[i] = i;
    return r;
  }

  private static DM row(double... v) {
    return new DM(new double[][] {v});
  }

  private static DM ev(Function f, Object... args) {
    return (DM) f.invoke(args);
  }

  @SuppressWarnings("unchecked")
  private static List<Object> evs(Function f, Object... args) {
    return (List<Object>) f.invoke(args);
  }

  private static Function fun(String name, List<MX> in, List<MX> out) {
    return new Function(name, in, out);
  }

  private static Function jacobianOld(Function f, int i, int j) {
    List<String> out = new ArrayList<>();
    out.add("jac:" + f.name_out(j) + ":" + f.name_in(i));
    out.addAll(f.name_out());
    return f.factory(f.name() + "_jac", f.name_in(), out);
  }

  private static double[][] npT(double[][] a) {
    double[][] r = new double[a[0].length][a.length];
    for (int i = 0; i < a.length; i++)
      for (int j = 0; j < a[0].length; j++) r[j][i] = a[i][j];
    return r;
  }

  /** numpy.reshape, row-major */
  private static double[][] npReshape(double[][] a, int n, int m) {
    double[] flat = new double[a.length * a[0].length];
    int k = 0;
    for (double[] r : a)
      for (double v : r) flat[k++] = v;
    double[][] r = new double[n][m];
    for (int i = 0; i < n; i++)
      for (int j = 0; j < m; j++) r[i][j] = flat[i * m + j];
    return r;
  }

  private static double[][] npScale(double[][] a, double s) {
    double[][] r = new double[a.length][a[0].length];
    for (int i = 0; i < a.length; i++)
      for (int j = 0; j < a[0].length; j++) r[i][j] = a[i][j] * s;
    return r;
  }

  private static double[][] npVstack(double[][] a, double[][] b) {
    double[][] r = new double[a.length + b.length][];
    System.arraycopy(a, 0, r, 0, a.length);
    System.arraycopy(b, 0, r, a.length, b.length);
    return r;
  }

  private static double[][] npHstack(double[][] a, double[][] b) {
    double[][] r = new double[a.length][a[0].length + b[0].length];
    for (int i = 0; i < a.length; i++) {
      System.arraycopy(a[i], 0, r[i], 0, a[i].length);
      System.arraycopy(b[i], 0, r[i], a[i].length, b[i].length);
    }
    return r;
  }

  private static void checkMXoperations(UnaryOperator<MX> ztf, UnaryOperator<double[][]> zrf, String name) {
    MX x = MX.sym("x", 1, 3);
    List<MX> zs = new ArrayList<>();
    for (int i = 0; i < 8; i++) zs.add(x.times(i + 1.0));
    MX z = ca.vertcat(zs);
    Function f = fun("f", List.of(x), List.of(ztf.apply(z)));
    double[] L = {1, 2, 3};
    DM fOut = ev(f, (Object) L);
    double[][] zr = new double[8][];
    for (int i = 0; i < 8; i++) zr[i] = new double[] {L[0] * (i + 1), L[1] * (i + 1), L[2] * (i + 1)};
    checkarray(fOut, new DM(zrf.apply(zr)), name);
  }

  private static void checkMXoperations2(UnaryOperator<MX> ztf, UnaryOperator<double[][]> zrf, String name) {
    MX x = MX.sym("x", 3, 1);
    List<MX> zs = new ArrayList<>();
    for (int i = 0; i < 8; i++) zs.add(x.times((double) i));
    MX z = ca.horzcat(zs);
    Function f = fun("f", List.of(x), List.of(ztf.apply(z)));
    double[] L = {1, 2, 3};
    DM fOut = ev(f, (Object) L);
    double[][] zr = new double[8][];
    for (int i = 0; i < 8; i++) zr[i] = new double[] {L[0] * i, L[1] * i, L[2] * i};
    checkarray(fOut, new DM(zrf.apply(npT(zr))), name);
  }

  private static void checkMXoperations3(UnaryOperator<MX> ztf, UnaryOperator<double[][]> zrf, String name) {
    MX x = MX.sym("x", 3, 1);
    MX p = ca.horzcat(List.of(x.get(0, 0), x.get(1, 0), x.get(2, 0)));
    List<MX> zs = new ArrayList<>();
    for (int i = 0; i < 8; i++) zs.add(p.times((double) i));
    MX z = ca.vertcat(zs);
    Function f = fun("f", List.of(x), List.of(ztf.apply(z)));
    double[] L = {1, 2, 3};
    DM fOut = ev(f, (Object) L);
    double[][] zr = new double[8][];
    for (int i = 0; i < 8; i++) zr[i] = new double[] {L[0] * i, L[1] * i, L[2] * i};
    checkarray(fOut, new DM(zrf.apply(zr)), name);
  }

  /** helpers.py evaluationCheck/numpyEvaluationCheck: Function(x, ft(x)) at setx0 against the reference */
  private void numpyEvaluationCheck(UnaryOperator<MX> ft, DM ref, List<MX> x, List<DM> setx0, String name,
      UnaryOperator<Function> fmod) {
    message(":" + name);
    Function f = fun("f", x, List.of(ft.apply(x.get(0))));
    if (fmod != null) f = fmod.apply(f);
    List<DM> fIn = new ArrayList<>(setx0);
    while (fIn.size() < f.n_in()) fIn.add(new DM());
    DM zt = new DM(f.call(fIn).get(0).full());
    checkarray(zt, ref, name);
  }

  private void numpyEvaluationCheckPool(List<Op1> pool, MX x, DM x0, String name, UnaryOperator<Function> fmod,
      DM setx0, boolean excludeNozero) {
    for (Op1 op : pool) {
      if (excludeNozero && op.nozero) continue;
      numpyEvaluationCheck(op.cas, op.ref.apply(ca.densify(x0)), List.of(x), List.of(setx0 == null ? x0 : setx0),
          name + ":" + op.name, fmod);
    }
  }

  private void numpyEvaluationCheckPool2(List<Op2> pool, MX x, MX y, DM x0, DM y0, String name) {
    for (Op2 op : pool) {
      message(":" + name + ":" + op.name);
      Function f = fun("f", List.of(x, y), List.of(op.cas.apply(x, y)));
      DM zt = new DM(f.call(List.of(x0, y0)).get(0).full());
      checkarray(zt, op.ref.apply(ca.densify(x0), ca.densify(y0)), name + ":" + op.name);
    }
  }

  /** helpers.py checkfunction_light: outputs agree */
  private static void checkfunctionLight(Function f, Function g, List<DM> inputs) {
    List<DM> a = f.call(inputs), b = g.call(inputs);
    assertEquals(a.size(), b.size());
    for (int i = 0; i < a.size(); i++) checkarray(a.get(i), b.get(i), 1e-9, "output " + i);
  }

  private static void checkSerialize(Function f, List<DM> inputs) {
    Function f2 = Function.deserialize(f.serialize(Map.of("debug", true)));
    List<DM> a = f.call(inputs), b = f2.call(inputs);
    for (int i = 0; i < f.n_out(); i++) checkarray(a.get(i), b.get(i), 1e-15, "serialize output " + i);
  }

  private void checkCodegen(Function f, List<DM> inputs) {
    checkCodegen(f, inputs, new HashMap<>(), List.of());
  }

  /** helpers.py check_codegen: generate, compile with gcc, load as external, compare (RUN_SLOW only) */
  private void checkCodegen(Function f, List<DM> inputs, Map<String, Object> opts, List<String> extralibs) {
    if (!RUN_SLOW) return;
    try {
      Path dir = Files.createTempDirectory("casadi_cg");
      String name = "codegen_" + Long.toHexString(random.nextLong() & Long.MAX_VALUE);
      CodeGenerator cg = new CodeGenerator(name, opts);
      cg.add(f);
      cg.generate(dir.toString() + File.separator);
      String libdir = GlobalOptions.getCasadiPath();
      String inc = GlobalOptions.getCasadiIncludePath();
      List<String> cmd = new ArrayList<>(Arrays.asList("gcc", "-pedantic", "-std=c89", "-fPIC", "-shared", "-Wall",
          "-Werror", "-Wextra", "-I" + inc, "-I" + inc + "/highs", "-Wno-unknown-pragmas", "-Wno-long-long",
          "-Wno-unused-parameter", "-O3", name + ".c", "-o", name + ".so", "-L" + libdir, "-Wl,-rpath," + libdir));
      for (String lib : extralibs) {
        if (lib.startsWith("-") || lib.contains(".")) cmd.add(lib);
        else cmd.add("-l" + lib);
      }
      Process p = new ProcessBuilder(cmd).directory(dir.toFile()).inheritIO().start();
      assertEquals(0, p.waitFor(), "compile " + String.join(" ", cmd));
      Function f2 = ca.external(f.name(), dir.resolve(name + ".so").toString());
      List<DM> a = f.call(inputs), b = f2.call(inputs);
      for (int i = 0; i < f.n_out(); i++) checkarray(a.get(i), b.get(i), 1e-15, "codegen output " + i);
    } catch (IOException | InterruptedException e) {
      throw new RuntimeException(e);
    }
  }

  private static Path tmpDir() {
    try {
      return Files.createTempDirectory("casadi_mx");
    } catch (IOException e) {
      throw new RuntimeException(e);
    }
  }

  private static String readFile(Path p) {
    try {
      return new String(Files.readAllBytes(p), StandardCharsets.UTF_8).replace("\r\n", "\n");
    } catch (IOException e) {
      throw new RuntimeException(e);
    }
  }

  // ---------------------------------------------------------------- tests

  @Test
  public void test_MX1() {
    message("MX constructor");
    MX x = MX.sym("x", 2, 3);
    assertEquals(2, x.size1(), "MX fails to indicate its size1");
    assertEquals(3, x.size2(), "MX fails to indicate its size2");
  }

  @Test
  public void test_MXvertcat() {
    message("MX vertcat");
    MX x = MX.sym("x", 1, 3);
    MX y = MX.sym("y", 1, 3);
    MX z = ca.vertcat(List.of(x, y));
    assertEquals(2, z.size1(), "MX fails to indicate its size1");
    assertEquals(3, z.size2(), "MX fails to indicate its size2");
  }

  @Test
  public void test_MX_fun1() {
    message("MXFunction single input, single output");
    MX x = MX.sym("x");
    MX y = ca.times(2.0, x);
    Function f = fun("f", List.of(x), List.of(y));
    assertEquals(1, f.n_in(), "Function fails to indicate correct number of inputs");
    assertEquals(1, f.n_out(), "Function fails to indicate correct number of outputs");
    DM fOut = ev(f, 3);
    double[] yt = fOut.nonzeros();
    assertEquals(1, yt.length, "Output of Function was tuple of floats, as expected, but length is incorrect.");
    checkarray(yt[0], 2 * 3, 1e-10);
  }

  @Test
  public void test_MXfunction2() {
    message("Function multi input, multi output");
    MX x = MX.sym("x");
    MX y = MX.sym("y");
    Function f = fun("f", List.of(x, y), List.of(x.plus(y), y.times(x)));
    assertEquals(2, f.n_in(), "Function fails to indicate correct number of inputs");
    assertEquals(2, f.n_out(), "Function fails to indicate correct number of outputs");
    List<Object> fOut = evs(f, 3, 7);
    double[] zt1 = ((DM) fOut.get(0)).nonzeros();
    double[] zt2 = ((DM) fOut.get(1)).nonzeros();
    assertEquals(1, zt1.length);
    assertEquals(1, zt2.length);
    checkarray(zt2[0], 21, 1e-10);
    checkarray(zt1[0], 10, 1e-10);
  }

  @Test
  public void test_MXfunction3() {
    message("Function single input, multi output (1)");
    MX xy = MX.sym("xy", 2);
    Function f = fun("f", List.of(xy), List.of(xy.get(0).plus(xy.get(1)), xy.get(0).times(xy.get(1))));
    assertEquals(1, f.n_in(), "Function fails to indicate correct number of inputs");
    assertEquals(2, f.n_out(), "Function fails to indicate correct number of outputs");
    List<Object> fOut = evs(f, (Object) new double[] {3, 7});
    double[] zt1 = ((DM) fOut.get(0)).nonzeros();
    double[] zt2 = ((DM) fOut.get(1)).nonzeros();
    assertEquals(1, zt1.length);
    assertEquals(1, zt2.length);
    checkarray(zt2[0], 21, 1e-10);
    checkarray(zt1[0], 10, 1e-10);
  }

  @Test
  public void test_MXfunction3b() {
    message("Function single input, multi output (2)");
    MX xy = MX.sym("xy", 1, 2);
    Function f = fun("f", List.of(xy), List.of(xy.get(0, 0).plus(xy.get(0, 1)), xy.get(0, 0).times(xy.get(0, 1))));
    assertEquals(1, f.n_in(), "Function fails to indicate correct number of inputs");
    assertEquals(2, f.n_out(), "Function fails to indicate correct number of outputs");
    List<Object> fOut = evs(f, (Object) new double[] {3, 7});
    double[][] zt1 = ((DM) fOut.get(0)).full();
    double[][] zt2 = ((DM) fOut.get(1)).full();
    assertEquals(1, zt1.length, "Output of Function is of wrong shape.");
    assertEquals(1, zt1[0].length, "Output of Function is of wrong shape.");
    assertEquals(1, zt2.length, "Output of Function is of wrong shape.");
    assertEquals(1, zt2[0].length, "Output of Function is of wrong shape.");
    checkarray(zt2[0][0], 21, 1e-10);
    checkarray(zt1[0][0], 10, 1e-10);
  }

  @Test
  public void test_MXfunction4() {
    message("Function single input, single output , using vertcat");
    MX xy = MX.sym("xy", 2);
    MX z = ca.vertcat(List.of(xy.get(0).plus(xy.get(1)), xy.get(0).times(xy.get(1))));
    Function f = fun("f", List.of(xy), List.of(z));
    assertEquals(1, f.n_in(), "Function fails to indicate correct number of inputs");
    assertEquals(1, f.n_out(), "Function fails to indicate correct number of outputs");
    double[][] zt = ev(f, (Object) new double[] {3, 7}).full();
    assertEquals(2, zt.length, "Output of Function is of wrong shape.");
    assertEquals(1, zt[0].length, "Output of Function is of wrong shape.");
    checkarray(zt[1][0], 21, 1e-10);
    checkarray(zt[0][0], 10, 1e-10);
  }

  @Test
  public void test_MXfunction5() {
    message("Function single input, single output , using horzcat");
    MX xy = MX.sym("xy", 2);
    MX z = ca.horzcat(List.of(xy.get(0).plus(xy.get(1)), xy.get(0).times(xy.get(1))));
    Function f = fun("f", List.of(xy), List.of(z));
    assertEquals(1, f.n_in(), "Function fails to indicate correct number of inputs");
    assertEquals(1, f.n_out(), "Function fails to indicate correct number of outputs");
    double[][] zt = ev(f, (Object) new double[] {3, 7}).full();
    assertEquals(1, zt.length, "Output of Function is of wrong shape.");
    assertEquals(2, zt[0].length, "Output of Function is of wrong shape.");
    checkarray(zt[0][1], 21, 1e-10);
    checkarray(zt[0][0], 10, 1e-10);
  }

  @Test
  public void test_which_depends_empty() {
    SX xs = SX.sym("x");
    MX xm = MX.sym("x");
    for (boolean tr : new boolean[] {true, false}) {
      for (int i = 0; i < 3; i++) {
        assertEquals(tr ? 1 : 0, ca.which_depends(xs, new SX(0, 1), i, tr).length);
        assertFalse(anyTrue(ca.which_depends(xs, new SX(0, 1), i, tr)));
        assertEquals(tr ? 0 : 1, ca.which_depends(new SX(0, 1), xs, i, tr).length);
        assertFalse(anyTrue(ca.which_depends(new SX(0, 1), xs, i, tr)));
        assertEquals(0, ca.which_depends(new SX(0, 1), new SX(0, 1), i, tr).length);

        assertEquals(tr ? 1 : 0, ca.which_depends(xm, new MX(0, 1), i, tr).length);
        assertFalse(anyTrue(ca.which_depends(xm, new MX(0, 1), i, tr)));
        assertEquals(tr ? 0 : 1, ca.which_depends(new MX(0, 1), xm, i, tr).length);
        assertFalse(anyTrue(ca.which_depends(new MX(0, 1), xm, i, tr)));
        assertEquals(0, ca.which_depends(new MX(0, 1), new MX(0, 1), i, tr).length);
      }
    }
  }

  private static boolean anyTrue(boolean[] b) {
    for (boolean v : b) if (v) return true;
    return false;
  }

  @Test
  public void test_issue83() {
    MX x = MX.sym("x");
    MX y = MX.sym("y");
    MX z = x.plus(y);
    Function f = fun("f", List.of(x, y), List.of(z));
    MX fc = (MX) f.invoke(new MX(3), y);
    Function g = fun("g", List.of(y), List.of(fc));
    DM gOut = ev(g, 7);
    checkarray(gOut.to_double(), 10, 1e-10);

    fc = (MX) f.invoke(x, new MX(7));
    g = fun("g", List.of(x), List.of(fc));
    gOut = ev(g, 3);
    checkarray(gOut.to_double(), 10, 1e-10);
  }

  @Test
  public void test_identitySX() {
    message("identity SXFunction");
    SX x = SX.sym("x");
    Function f = new Function("f", List.of(x), List.of(x));
    DM fOut = ev(f, 3);
    checkarray(fOut.get(0, 0).to_double(), 3, 1e-10);
  }

  @Test
  public void test_identityMX() {
    message("identity Function");
    MX x = MX.sym("x");
    Function f = fun("f", List.of(x), List.of(x));
    DM fOut = ev(f, 3);
    checkarray(fOut.get(0, 0).to_double(), 3, 1e-10);
  }

  @Test
  public void test_MXorder() {
    message("Function order of non-zero elements");
    MX x = MX.sym("x", 2, 3);
    Function f = fun("f", List.of(x), List.of(x.plus(x)));
    assertEquals(1, f.n_in(), "Function fails to indicate correct number of inputs");
    assertEquals(1, f.n_out(), "Function fails to indicate correct number of outputs");
    double[] L = {1, 2, 3, 4, 5, 6};
    DM fIn = new DM(f.sparsity_in(0), new DM(L));
    double[][] zt = ev(f, fIn).full();
    assertEquals(2, zt.length, "Output of Function is of wrong shape.");
    assertEquals(3, zt[0].length, "Output of Function is of wrong shape.");
    for (int i = 0; i < 2; i++)
      for (int j = 0; j < 3; j++) checkarray(L[i + 2 * j] * 2, zt[i][j], 1e-10);
  }

  @Test
  public void test_trans() {
    message("trans");
    MX a = new MX(0, 1);
    MX b = a.T();
    assertEquals(1, b.size1());
    assertEquals(0, b.size2());
  }

  @Test
  public void test_MXtrans() {
    message("trans(MX)");
    MX x = MX.sym("x", 2, 3);
    MX z = x.T();
    assertEquals(3, z.size1(), "Vec returns MX of wrong dimension");
    assertEquals(2, z.size2(), "Vec returns MX of wrong dimension");
    Function f = fun("f", List.of(x), List.of(z));
    assertEquals(1, f.n_in(), "Function fails to indicate correct number of inputs");
    assertEquals(1, f.n_out(), "Function fails to indicate correct number of outputs");
    double[] L = {1, 2, 3, 4, 5, 6};
    DM fIn = new DM(f.sparsity_in(0), new DM(L));
    double[][] ztr = ev(f, fIn).full();
    for (int i = 0; i < 2; i++)
      for (int j = 0; j < 3; j++) checkarray(L[i + 2 * j], ztr[j][i], 1e-10);
  }

  @Test
  public void test_MXvec() {
    double[][] u0 = new double[4][3];
    for (int j = 0; j < 4; j++)
      for (int i = 0; i < 3; i++) u0[j][i] = 10 * j + i;
    DM u = new DM(u0);
    MX U = MX.sym("u", u.size());
    Function f = fun("f", List.of(U), List.of(ca.vec(U)));
    DM fOut = ev(f, u);
    checkarray(ca.vec(u), fOut, "vec");
  }

  @Test
  public void test_MXreshape() {
    message("reshape(MX)");
    MX x = MX.sym("x", 2, 3);
    MX z = ca.reshape(x, new long[] {1, 6});
    assertEquals(1, z.size1(), "Vec returns MX of wrong dimension");
    assertEquals(6, z.size2(), "Vec returns MX of wrong dimension");
    Function f = fun("f", List.of(x), List.of(z));
    assertEquals(1, f.n_in(), "Function fails to indicate correct number of inputs");
    assertEquals(1, f.n_out(), "Function fails to indicate correct number of outputs");
    double[] L = {1, 2, 3, 4, 5, 6};
    DM fIn = new DM(f.sparsity_in(0), new DM(L));
    double[][] zt = ev(f, fIn).full();
    for (int i = 0; i < L.length; i++) checkarray(L[i], zt[0][i], 1e-10);
  }

  @Test
  public void test_sparsity_cast() {
    Sparsity spSource = ca.sparsify(dm(new double[][] {{1, 0, 1}, {0, 0, 1}, {1, 0, 0}})).sparsity();
    MX x = MX.sym("x", spSource);
    SX xsx = SX.sym("x", spSource);
    Sparsity sp = ca.sparsify(dm(new double[][] {{1, 0, 1}, {1, 0, 1}})).sparsity();
    assertInException("mismatch", () -> ca.reshape(x, sp));
    assertInException("Mismatching", () -> ca.sparsity_cast(x, ca.horzcat(List.of(DM.ones(sp), DM.ones(sp))).sparsity()));  // horzcat(List<Sparsity>) does not compile
    MX y = ca.sparsity_cast(x, sp);
    SX ysx = ca.sparsity_cast(xsx, sp);
    Function fsx = new Function("fsx", List.of(xsx), List.of(ysx));
    Function f = fun("f", List.of(x), List.of(y));
    DM inp = dm(new double[][] {{1, 0, 2}, {0, 0, 3}, {4, 0, 0}});
    checkfunction(f, fsx, List.of(inp));
    checkCodegen(f, List.of(inp));
    DM y2 = ev(fsx, inp);
    checkarray(y2, dm(new double[][] {{1, 0, 2}, {4, 0, 3}}));

    Sparsity sp2 = ca.sparsify(dm(new double[][] {{1, 0}, {0, 1}})).sparsity();
    MX y3 = ca.sparsity_cast(MX.sym("y", 2), sp2);
    ca.vec(y3);
  }

  @Test
  public void test_MXcompose() {
    message("compositions of vec, trans, reshape with vertcat");
    checkMXoperations(x -> x, x -> x, "vertcat");
    checkMXoperations(x -> x.T(), MXTest::npT, "trans(vertcat)");
    checkMXoperations(x -> x.T().T(), x -> x, "trans(trans(vertcat))");
    checkMXoperations(x -> ca.vec(x.T()), x -> npReshape(x, x.length * x[0].length, 1), "vec(trans(vertcat))");
    checkMXoperations(x -> ca.vec(x).T(), x -> npT(npReshape(npT(x), x.length * x[0].length, 1)), "vec(trans(vertcat))");
    checkMXoperations(x -> ca.reshape(x.T(), 6, 4).T(), x -> npReshape(x, 4, 6), "reshape(vertcat)");
    checkMXoperations(x -> ca.reshape(x, 6, 4).T(), x -> npReshape(npT(x), 4, 6), "reshape(trans(vertcat))");
    checkMXoperations(x -> ca.reshape(x.T(), 6, 4), x -> npT(npReshape(x, 4, 6)), "trans(reshape(vertcat))");
  }

  @Test
  public void test_MXcompose2() {
    message("compositions of vec, trans, reshape with horzcat");
    checkMXoperations2(x -> x, x -> x, "horzcat");
    checkMXoperations2(x -> x.T(), MXTest::npT, "trans(horzcat)");
    checkMXoperations2(x -> x.T().T(), x -> x, "trans(trans(horzcat))");
    checkMXoperations2(x -> ca.vec(x.T()), x -> npReshape(x, x.length * x[0].length, 1), "vec(trans(horzcat))");
    checkMXoperations2(x -> ca.vec(x).T(), x -> npT(npReshape(npT(x), x.length * x[0].length, 1)), "vec(trans(horzcat))");
    checkMXoperations2(x -> ca.reshape(x.T(), 6, 4).T(), x -> npReshape(x, 4, 6), "reshape(horzcat)");
    checkMXoperations2(x -> ca.reshape(x, 6, 4).T(), x -> npReshape(npT(x), 4, 6), "reshape(trans(horzcat))");
    checkMXoperations2(x -> ca.reshape(x.T(), 6, 4), x -> npT(npReshape(x, 4, 6)), "trans(reshape(horzcat))");
  }

  @Test
  public void test_MXcompose3() {
    message("compositions of vec, trans, reshape with vertcat");
    checkMXoperations3(x -> x, x -> x, "snippet");
    checkMXoperations3(x -> x.T(), MXTest::npT, "trans(snippet)");
    checkMXoperations3(x -> x.T().T(), x -> x, "trans(trans(snippet))");
    checkMXoperations3(x -> ca.vec(x.T()), x -> npReshape(x, x.length * x[0].length, 1), "vec(trans(snippet))");
    checkMXoperations3(x -> ca.vec(x).T(), x -> npT(npReshape(npT(x), x.length * x[0].length, 1)), "vec(trans(snippet))");
    checkMXoperations3(x -> ca.reshape(x.T(), 6, 4).T(), x -> npReshape(x, 4, 6), "reshape(snippet)");
    checkMXoperations3(x -> ca.reshape(x, 6, 4).T(), x -> npReshape(npT(x), 4, 6), "reshape(trans(snippet))");
    checkMXoperations3(x -> ca.reshape(x.T(), 6, 4), x -> npT(npReshape(x, 4, 6)), "trans(reshape(snippet))");
  }

  @Test
  public void test_MXcompose4() {
    message("compositions of horzcat + vertcat");
    checkMXoperations(x -> ca.vertcat(List.of(x)), x -> x, "vertcat(*vertcat)");
    checkMXoperations(x -> ca.vertcat(List.of(x, x.times(2.0))), x -> npVstack(x, npScale(x, 2)), "vertcat(*vertcat,vertcat)");
    checkMXoperations(x -> ca.horzcat(List.of(x)), x -> x, "horzcat(*vertcat)");
    checkMXoperations(x -> ca.horzcat(List.of(x, x.times(2.0))), x -> npHstack(x, npScale(x, 2)), "horzcat(*vertcat,vertcat)");

    checkMXoperations2(x -> ca.vertcat(List.of(x)), x -> x, "vertcat(*horzcat)");
    checkMXoperations2(x -> ca.vertcat(List.of(x, x.times(2.0))), x -> npVstack(x, npScale(x, 2)), "vertcat(*horzcat,horzcat)");
    checkMXoperations2(x -> ca.horzcat(List.of(x)), x -> x, "horzcat(*horzcat)");
    checkMXoperations2(x -> ca.horzcat(List.of(x, x.times(2.0))), x -> npHstack(x, npScale(x, 2)), "horzcat(*horzcat,horzcat)");

    checkMXoperations3(x -> ca.vertcat(List.of(x)), x -> x, "vertcat(*snippet)");
    checkMXoperations3(x -> ca.vertcat(List.of(x, x.times(2.0))), x -> npVstack(x, npScale(x, 2)), "vertcat(*snippet,snippet)");
    checkMXoperations3(x -> ca.horzcat(List.of(x)), x -> x, "horzcat(*snippet)");
    checkMXoperations3(x -> ca.horzcat(List.of(x, x.times(2.0))), x -> npHstack(x, npScale(x, 2)), "horzcat(*snippet,snippet)");
  }

  @Test
  public void test_MXslicingnew() {
    Assumptions.assumeTrue(KNOWN_BUGS, "known bug");
    message("MX slicing new");
    message(":dense");
    MX x = MX.sym("x", 3, 2);
    DM x0 = dm(new double[][] {{0.738, 0.2}, {0.1, 0.39}, {0.99, 0.999999}});
    List<MX> X = List.of(x);
    List<DM> X0 = List.of(x0);
    numpyEvaluationCheck(a -> a.get(1, 0), x0.get(1, 0), X, X0, "x[1,0]", null);
    numpyEvaluationCheck(a -> a.get(0, 1), x0.get(0, 1), X, X0, "x[0,1]", null);
    numpyEvaluationCheck(a -> a.get(0, -1), x0.get(0, -1), X, X0, "x[0,-1]", null);
    numpyEvaluationCheck(a -> a.get(new Slice(), 0), x0.get(new Slice(), 0), X, X0, "x[:,0]", null);
    numpyEvaluationCheck(a -> a.get(new Slice(), 1), x0.get(new Slice(), 1), X, X0, "x[:,1]", null);
    numpyEvaluationCheck(a -> a.get(1, new Slice()), x0.get(1, new Slice()), X, X0, "x[1,:]", null);
    numpyEvaluationCheck(a -> a.get(0, new Slice()), x0.get(0, new Slice()), X, X0, "x[0,:]", null);
    numpyEvaluationCheck(a -> a.get(-1, new Slice()), x0.get(-1, new Slice()), X, X0, "x[-1,:]", null);
    numpyEvaluationCheck(a -> a.get(new Slice(), -2), x0.get(new Slice(), -2), X, X0, "x[:,-2]", null);
    numpyEvaluationCheck(a -> a.get(new Slice(0, -2), new Slice(0, -1)), x0.get(new Slice(0, -2), new Slice(0, -1)), X, X0, "x[0:-2,0:-1]", null);
    numpyEvaluationCheck(a -> a.get(new Slice(0, 2), new Slice(0, 2)), x0.get(new Slice(0, 2), new Slice(0, 2)), X, X0, "x[0:2,0:2]", null);
    numpyEvaluationCheck(a -> a.get(idx(0, 1), new Slice(0, 2)), x0.get(idx(0, 1), new Slice(0, 2)), X, X0, "x[[0,1],0:2]", null);
    numpyEvaluationCheck(a -> a.get_nz(false, idx(0, 2, 3)), dm(0.738, 0.99, 0.2), X, X0, "x[[0,2,3]]", null);
    numpyEvaluationCheck(a -> a.get_nz(false, 1), new DM(0.1), X, X0, "x[1] on dense matrix", null);
    numpyEvaluationCheck(a -> a.get_nz(false, -1), new DM(0.999999), X, X0, "x[-1] on dense matrix", null);
    numpyEvaluationCheck(a -> a.get(idx(0, 1), new Slice(0, 1)), x0.get(idx(0, 1), new Slice(0, 1)), X, X0, "x[:,0:1]", null);
    numpyEvaluationCheck(a -> a.get(new Slice(0, 1), idx(0, 1)), x0.get(new Slice(0, 1), idx(0, 1)), X, X0, "x[0:1,:]", null);

    message(":sparse");
    Sparsity sp = new Sparsity(4, 3, new long[] {0, 2, 2, 3}, new long[] {1, 2, 1});
    x = MX.sym("X", sp);
    X = List.of(x);
    DM sx0 = dm(0.738, 0.39, 0.99);
    x0 = new DM(dmsp(sp, 0.738, 0.39, 0.99).full());
    List<DM> S0 = List.of(sx0);
    numpyEvaluationCheck(a -> a.get(0, 0), x0.get(0, 0), X, S0, "x[0,0]", null);
    numpyEvaluationCheck(a -> a.get(1, 0), x0.get(1, 0), X, S0, "x[1,0]", null);
    numpyEvaluationCheck(a -> a.get(0, 1), x0.get(0, 1), X, S0, "x[1,0]", null);
    numpyEvaluationCheck(a -> a.get(0, -1), x0.get(0, -1), X, S0, "x[0,-1]", null);
    numpyEvaluationCheck(a -> a.get(new Slice(), 0), x0.get(new Slice(), 0), X, S0, "x[:,0]", null);
    numpyEvaluationCheck(a -> a.get(new Slice(), 1), x0.get(new Slice(), 1), X, S0, "x[:,1]", null);
    numpyEvaluationCheck(a -> a.get(1, new Slice()), x0.get(1, new Slice()), X, S0, "x[1,:]", null);
    numpyEvaluationCheck(a -> a.get(0, new Slice()), x0.get(0, new Slice()), X, S0, "x[0,:]", null);
    numpyEvaluationCheck(a -> a.get(-1, new Slice()), x0.get(-1, new Slice()), X, S0, "x[-1,:]", null);
    numpyEvaluationCheck(a -> a.get(new Slice(), -2), x0.get(new Slice(), -2), X, S0, "x[:,-2]", null);
    numpyEvaluationCheck(a -> a.get(new Slice(0, -2), new Slice(0, -1)), x0.get(new Slice(0, -2), new Slice(0, -1)), X, S0, "x[0:-2,0:-1]", null);
    numpyEvaluationCheck(a -> a.get(new Slice(0, 2), new Slice(0, 2)), x0.get(new Slice(0, 2), new Slice(0, 2)), X, S0, "x[0:2,0:2]", null);
    numpyEvaluationCheck(a -> a.get(idx(0, 1), new Slice(0, 2)), x0.get(idx(0, 1), new Slice(0, 2)), X, S0, "x[[0,1],0:2]", null);
    numpyEvaluationCheck(a -> a.get_nz(false, idx(2, 1)), ca.vertcat(List.of(x0.get(1, 2), x0.get(2, 0))), X, List.of(x0), "x[[2,1]]", null);
    numpyEvaluationCheck(a -> a.get_nz(false, new Slice(0, 2)), dm(0.738, 0.39), X, List.of(x0), "x[0:2] on dense matrix", null);
    numpyEvaluationCheck(a -> a.get_nz(false, 1), new DM(0.39), X, S0, "x[1]", null);
    numpyEvaluationCheck(a -> a.get_nz(false, -1), new DM(0.99), X, S0, "x[-1]", null);
  }

  @Test
  public void test_mx_in() {
    message("mx_out/mx_in");
    MX x = MX.sym("x", 2, 3);
    Function f = fun("f", List.of(x), List.of(ca.times(3.0, x)));
    List<MX> xIn = f.mx_in();
    List<MX> xOut = f.call(xIn);
    Function g = fun("g", List.of(xIn.get(0)), List.of(ca.times(6.0, xOut.get(0))));
    double[] n = {1, 2, 3, 4, 5, 6};
    DM fOut = ev(f, new DM(f.sparsity_in(0), new DM(n)));
    DM gOut = ev(g, new DM(g.sparsity_in(0), new DM(n)));
    checkarray(gOut, ca.times(6.0, fOut), "slicing(trans)");
  }

  @Test
  public void test_scalarMX() {
    MX x = MX.sym("x");
    DM x0 = new DM(0.738);
    numpyEvaluationCheckPool(pool, x, x0, "scalarMX", null, null, false);
    // matrixpool is empty in mx.py
  }

  @Test
  public void test_MXJacobian() {
    message("MX(1,1) unary operation, jacobian");
    Jpool = new ArrayList<>();
    message("SX(1,1) unary operation, jacobian");
    MX x = MX.sym("x");
    DM x0 = new DM(new double[][] {{0.738}});
    numpyEvaluationCheckPool(Jpool, x, x0, "MX unary operations, jacobian", f -> jacobianOld(f, 0, 0), null, false);
  }

  @Test
  public void test_MXJacobians() {
    message("MX(3,1) unary operation, jacobian");
    MX x = MX.sym("x", 3, 1);
    DM x0 = dm(0.738, 0.9, 0.3);
    numpyEvaluationCheckPool(Jpool, x, x0, "MX unary operations, jacobian", f -> jacobianOld(f, 0, 0), null, false);
  }

  @Test
  public void test_MXbinary() {
    message("MX binary operations");
    MX x = MX.sym("x", 3, 2);
    MX y = MX.sym("x", 3, 2);
    DM x0 = dm(new double[][] {{0.738, 0.2}, {0.1, 0.39}, {0.99, 0.999999}});
    DM y0 = dm(new double[][] {{1.738, 0.6}, {0.7, 12}, {0, -6}});
    numpyEvaluationCheckPool2(matrixbinarypool, x, y, x0, y0, "MX");
  }

  @Test
  public void test_MXSparse() {
    message("MX unary operations, sparse");
    Sparsity sp = new Sparsity(4, 3, new long[] {0, 2, 2, 3}, new long[] {1, 2, 1});
    MX x = MX.sym("x", sp);
    DM x0 = dmsp(sp, 0.738, 0.1, 0.99);
    numpyEvaluationCheckPool(pool, x, x0, "MX", null, x0, true);
  }

  @Test
  public void test_MXbinarySparse() {
    message("SX binary operations");
    Sparsity spx = new Sparsity(4, 3, new long[] {0, 2, 2, 3}, new long[] {1, 2, 1});
    Sparsity spy = new Sparsity(4, 3, new long[] {0, 2, 2, 3}, new long[] {0, 2, 3});
    MX xx = MX.sym("x", spx);
    MX yy = MX.sym("y", spy);
    DM x0 = dmsp(spx, 0.738, 0.1, 0.99);
    DM y0 = dmsp(spy, 1.738, 0.7, -6);
    numpyEvaluationCheckPool2(matrixbinarypool, xx, yy, x0, y0, "MX");
  }

  @Test
  public void test_symbolcheck() {
    message("Check if non-symbolic inputs are caught");
    assertThrows(RuntimeException.class, () -> fun("f", List.of(new MX(0)), List.of(MX.sym("x"))));
  }

  @Test
  public void test_unite() {
    message("unite operation");
    double[][] xn = new double[3][4];
    for (double[] r : xn)
      for (int j = 0; j < 4; j++) r[j] = random.nextDouble();
    MX x = new MX(3, 4);
    MX y = MX.sym("x", 3, 4);
    MX z = ca.unite(x, y);
    Function f = fun("f", List.of(y), List.of(z));
    checkarray(ev(f, (Object) xn), xn);

    Sparsity spx = new Sparsity(4, 3, new long[] {0, 2, 2, 3}, new long[] {1, 2, 1});
    Sparsity spy = new Sparsity(4, 3, new long[] {0, 1, 2, 3}, new long[] {0, 2, 2});
    DM nx = DM.zeros(spx);
    for (long k = 0; k < nx.nnz(); k++) nx.set_nz(random.nextDouble(), false, k);
    DM ny = DM.zeros(spy);
    for (long k = 0; k < nx.nnz(); k++) ny.set_nz(random.nextDouble(), false, k);
    double[][] nxn = nx.full(), nyn = ny.full();
    double[][] sum = new double[4][3];
    for (int i = 0; i < 4; i++)
      for (int j = 0; j < 3; j++) sum[i][j] = nxn[i][j] + nyn[i][j];
    MX xs = MX.sym("x", spx);
    MX ys = MX.sym("y", spy);
    z = ca.unite(xs, ys);
    f = fun("f", List.of(xs, ys), List.of(z));
    checkarray(ev(f, nx, ny), sum);
  }

  @Test
  public void test_imatrix_index() {
    message("IM indexing");
    MX X = MX.sym("x", 2, 2);
    MX Y = X.get_nz(false, new DM(new double[][] {{0, 2}, {1, 1}, {3, 3}}));
    Function f = fun("f", List.of(X), List.of(Y));
    DM fOut = ev(f, new DM(f.sparsity_in(0), dm(1, 2, 3, 4)));
    checkarray(fOut, new double[][] {{1, 3}, {2, 2}, {4, 4}});

    MX Y2 = X.get(new Slice(), new Slice());
    Y2.set_nz(dm(new double[][] {{9, 8}}), false, new DM(new double[][] {{0, 2}}));
    f = fun("f", List.of(X), List.of(Y2));
    fOut = ev(f, new DM(f.sparsity_in(0), dm(1, 2, 3, 4)));
    checkarray(fOut, new double[][] {{9, 8}, {2, 4}});
  }

  @Test
  public void test_subsass() {
    message("Check subscripted assignment");
    MX X = MX.sym("x", 2, 2);
    X.set(0, 0, new MX(5));
    X.set(0, 0, 5.0);
    X.set(new Slice(), 0, 8.0);

    MX x = MX.sym("X", 3, 4);
    double[][] xn = new double[3][4];
    for (double[] r : xn)
      for (int j = 0; j < 4; j++) r[j] = random.nextDouble();
    double[][] r = new double[7][8];
    long[] cols = {2, 4, 6, 7};
    MX y = MX.zeros(7, 8);
    y.set(new Slice(1, 4), idx(2, 4, 6, 7), x);
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 4; j++) r[1 + i][(int) cols[j]] = xn[i][j];
    Function fy = fun("fy", List.of(x), List.of(y));
    checkarray(ev(fy, (Object) xn), r);

    y = new MX(7, 8);
    y.set(new Slice(1, 4), idx(2, 4, 6, 7), x);
    fy = fun("fy", List.of(x), List.of(y));
    checkarray(ev(fy, (Object) xn), r);

    long[] kl = {2, 4, 5, 8};
    Sparsity s = y.sparsity();
    for (long k : kl) r[(int) s.row()[(int) k]][(int) s.get_col()[(int) k]] = 1.0;
    y.set_nz(new MX(1), false, idx(2, 4, 5, 8));
    fy = fun("fy", List.of(x), List.of(y));
    checkarray(ev(fy, (Object) xn), r);

    y.set_nz(x.get_nz(false, idx(0, 1, 2, 3)), false, idx(2, 4, 5, 8));
    s = y.sparsity();
    Sparsity sx = x.sparsity();
    int cnt = 0;
    for (long k : kl) {
      r[(int) s.row()[(int) k]][(int) s.get_col()[(int) k]] = xn[(int) sx.row()[cnt]][(int) sx.get_col()[cnt]];
      cnt++;
    }
    fy = fun("fy", List.of(x), List.of(y));
    checkarray(ev(fy, (Object) xn), r);
  }

  @Test
  public void test_erase() {
    message("Erase function");
    message(":dense");
    MX y = MX.sym("Y", 7, 8);
    MX z = y.times(2.0);
    z.erase(new long[] {1, 2, 3}, new long[] {2, 4, 6, 7});
    Function f = fun("f", List.of(y), List.of(z));
    double[] ones = new double[56];
    Arrays.fill(ones, 1);
    new DM(f.sparsity_in(0), new DM(ones));
    DM fOut = ev(f, 0);
    DM e = fOut;
    checkarray(fOut, e, "erase");
    message(":sparse");
  }

  private void algebraCheck(DM A_, DM b_, DM C_, DM D_, DM e_, DM x_, MX A, MX b, MX C, MX D, MX e, MX x) {
    MX Axb = ca.mtimes(A, x).plus(b);
    MX Dxe = ca.mtimes(D, x).plus(e);
    MX a = ca.mtimes(ca.mtimes(Axb.T(), C), Dxe);
    Function f = fun("f", List.of(x, A, b, C, D, e), List.of(a));
    DM fOut = ev(f, x_, A_, b_, C_, D_, e_);
    DM Axb_ = ca.mtimes(A_, x_).plus(b_);
    DM Dxe_ = ca.mtimes(D_, x_).plus(e_);
    DM f_ = ca.mtimes(ca.mtimes(Axb_.T(), C_), Dxe_);
    checkarray(fOut, f_, "evaluation");
    DM J_ = ca.mtimes(ca.mtimes(Dxe_.T(), C_.T()), A_).plus(ca.mtimes(ca.mtimes(Axb_.T(), C_), D_));
    for (int w = 0; w < 2; w++) {
      f = new Function("f", List.of(x, A, b, C, D, e), List.of(a), Map.of("ad_weight", w, "ad_weight_sp", w));
      Function J = jacobianOld(f, 0, 0);
      List<DM> JOut = J.call(List.of(x_, A_, b_, C_, D_, e_));
      checkarray(JOut.get(0), J_, "evaluation");
    }
  }

  @Test
  public void test_MXalgebraDense() {
    message("Test some dense algebraic properties of matrices");
    int n = 3, m = 4;
    algebraCheck(randDM(m, n), randDM(m, 1), randDM(m, m), randDM(m, n), randDM(m, 1), randDM(n, 1),
        MX.sym("A", m, n), MX.sym("b", m, 1), MX.sym("C", m, m), MX.sym("D", m, n), MX.sym("e", m, 1),
        MX.sym("x", n, 1));
  }

  private Sparsity randsparsity(int m, int n, boolean skipMid) {
    Sparsity sp = new Sparsity(m, n);
    for (int k = 0; k < (n * m) / 2; k++) {
      int i = random.nextInt(m), j = random.nextInt(n);
      if (!skipMid || (i != m / 2 && (n == 1 || j != n / 2))) sp.add_nz(i, j);
    }
    return sp;
  }

  private DM gentest(int m, int n, boolean skipMid, MX[] sym) {
    Sparsity As = randsparsity(m, n, skipMid);
    DM A_ = DM.zeros(As);
    for (long k = 0; k < As.nnz(); k++) A_.set_nz(random.nextDouble(), false, k);
    sym[0] = MX.sym("A", As);
    return A_;
  }

  private void sparseAlgebra(int n, int m, boolean skipMid) {
    MX[] A = new MX[1], b = new MX[1], C = new MX[1], D = new MX[1], e = new MX[1];
    DM A_ = gentest(m, n, skipMid, A), b_ = gentest(m, 1, skipMid, b), C_ = gentest(m, m, skipMid, C),
        D_ = gentest(m, n, skipMid, D), e_ = gentest(m, 1, skipMid, e);
    double[][] x0 = new double[n][1];
    for (double[] r : x0) r[0] = random.nextDouble();
    algebraCheck(A_, b_, C_, D_, e_, new DM(x0), A[0], b[0], C[0], D[0], e[0], MX.sym("x", n, 1));
  }

  @Test
  public void test_MXalgebraSparse() {
    message("Test some sparse algebraic properties of matrices");
    sparseAlgebra(3, 4, false);
  }

  @Test
  public void test_MXalgebraSparseSparse() {
    message("Test some sparse algebraic properties of matrices");
    sparseAlgebra(8, 10, true);
  }

  @Test
  public void test_chaining() {
    message("Chaining SX and MX together");
    SX x = SX.sym("x");
    SX y = x.power(3.0);
    Function f = new Function("f", List.of(x), List.of(y));
    Function J = jacobianOld(f, 0, 0);
    MX X = MX.sym("X");
    Function F = fun("F", List.of(X), J.call(List.of(X)));
    double x_ = 1.7;
    List<Object> FOut = evs(F, x_);
    checkarray((DM) FOut.get(0), 3 * x_ * x_);
  }

  @Test
  public void test_issue107() {
    message("Regression test for issue 107: +=");
    MX x = MX.sym("x");
    MX y = MX.sym("y");
    MX z = x;
    z = z.plus(y);
    assertTrue(x.is_symbolic());
    assertFalse(z.is_symbolic());
  }

  @Test
  public void test_MXd_trivial() {
    message("symbolic variables and constants jac");
    MX X = MX.sym("X", 10);
    MX V = MX.sym("V");
    MX J = ca.jacobian(X, X);
    assertEquals(10, J.nnz());
    assertEquals(10, J.size1());
    assertEquals(10, J.size2());

    Function g = fun("g", List.of(), List.of(J));
    checkarray(g.call(List.<DM>of()).get(0), DM.eye(10), "unit matrix");
    g = fun("g", List.of(), List.of(ca.jacobian(MX.eye(3), X)));
    checkarray(g.call(List.<DM>of()).get(0), DM.zeros(9, 10), "zero matrix");
    g = fun("g", List.of(), List.of(ca.jacobian(X, V)));
    checkarray(g.call(List.<DM>of()).get(0), DM.zeros(10, 1), "zero matrix");
    g = fun("g", List.of(), List.of(ca.jacobian(MX.eye(3), V)));
    checkarray(g.call(List.<DM>of()).get(0), DM.zeros(9, 1), "zero matrix");
  }

  @Test
  public void test_MXd_substractionl() {
    message("substraction jac");
    MX V = MX.sym("V");
    MX X = MX.sym("X");
    Function g = fun("g", List.of(), List.of(ca.jacobian(X.minus(V), X)));
    checkarray(g.call(List.<DM>of()).get(0), DM.ones(1, 1), "one");
    g = fun("g", List.of(), List.of(ca.jacobian(X.minus(V), V)));
    checkarray(g.call(List.<DM>of()).get(0), DM.ones(1, 1).unaryMinus(), "one");
    g = fun("g", List.of(), List.of(ca.jacobian(V.minus(X), X)));
    checkarray(g.call(List.<DM>of()).get(0), DM.ones(1, 1).unaryMinus(), "one");
    g = fun("g", List.of(), List.of(ca.jacobian(V.minus(X), V)));
    checkarray(g.call(List.<DM>of()).get(0), DM.ones(1, 1), "one");
  }

  @Test
  public void test_MXd_mapping() {
    message("mapping jac");
    MX X = MX.sym("X", 3);
    MX Y = MX.sym("Y", 2);
    MX J = ca.jacobian(ca.vertcat(List.of(X, Y)), X);
    DM JJ = DM.ones(J.sparsity());
    checkarray(JJ, ca.vertcat(List.of(DM.eye(3), DM.zeros(2, 3))), "diag");
    J = ca.jacobian(ca.vertcat(List.of(X, Y)), Y);
    JJ = DM.ones(J.sparsity());
    checkarray(JJ, ca.vertcat(List.of(DM.zeros(3, 2), DM.eye(2))), "diag");
  }

  @Test
  public void test_null() {
    message("Function null");
    MX x = MX.sym("x");
    Function f = fun("f", List.of(x), List.of(x.power(2.0), new MX()));
    assertEquals(0, f.size1_out(1));
    assertEquals(0, f.size2_out(1));
    evs(f, 0);

    f = fun("f", List.of(x, new MX()), List.of(x.power(2.0), new MX()));
    assertEquals(0, f.size1_out(1));
    assertEquals(0, f.size2_out(1));
    evs(f, 0, 0);

    List<Object> r = evs(f, x, new MX());
    assertTrue(((MX) r.get(1)).is_empty(true));

    r = evs(f, new MX(), new MX());
    assertTrue(((MX) r.get(1)).is_empty(true));
  }

  @Test
  public void test_issue184() {
    message("Regression test issue #184");
    MX x = MX.sym("x", 3);
    MX y = x.get(new Slice(0, 0));
    assertEquals(0, y.nnz());
  }

  @Test
  public void test_indexingOutOfBounds() {
    message("Indexing out of bounds");
    DM y = DM.zeros(4, 5);
    assertThrows(CasadiException.class, () -> y.get(12, 0));
    assertThrows(CasadiException.class, () -> y.get(12, 12));
    assertThrows(CasadiException.class, () -> y.get(0, 12));
    assertThrows(CasadiException.class, () -> y.get(12, new Slice()));
    assertThrows(CasadiException.class, () -> y.get(new Slice(12, 15), 0));
    assertThrows(CasadiException.class, () -> y.get(new Slice(), 12));
    assertThrows(CasadiException.class, () -> y.get(0, new Slice(12, 15)));
    y.get(-1, 2);
    assertThrows(CasadiException.class, () -> y.get(-12, 2));
    y.get(new Slice(-3, -1), 2);
    assertThrows(CasadiException.class, () -> y.get(new Slice(-12, -9), 2));

    assertThrows(CasadiException.class, () -> y.set(12, 0, 0.0));
    assertThrows(CasadiException.class, () -> y.set(12, 12, 0.0));
    assertThrows(CasadiException.class, () -> y.set(0, 12, 0.0));
    assertThrows(CasadiException.class, () -> y.set(12, new Slice(), 0.0));
    assertThrows(CasadiException.class, () -> y.set(new Slice(12, 15), 0, 0.0));
    assertThrows(CasadiException.class, () -> y.set(new Slice(), 12, 0.0));
    assertThrows(CasadiException.class, () -> y.set(0, new Slice(12, 15), 0.0));
    y.set(-1, 2, 0.0);
    assertThrows(CasadiException.class, () -> y.set(-12, 2, 0.0));
    y.set(new Slice(-3, -1), 2, 0.0);
    assertThrows(CasadiException.class, () -> y.set(new Slice(-12, -9), 2, 0.0));
  }

  @Test
  public void test_indexinglimits() {
    message("Limits of indexing");
    MX y = MX.sym("y", 3);
    assertThrows(CasadiException.class, () -> y.get(idx(0, 5)));
    assertThrows(CasadiException.class, () -> y.set(idx(0, 5), MX.sym("a")));
    y.get(idx(0, 2));
    y.set(idx(0, 2), MX.sym("a"));
  }

  @Test
  public void test_mtimes() {
    MX A = new MX(DM.ones(new long[] {4, 3}));
    MX B = new MX(DM.ones(new long[] {3, 8}));
    MX C = new MX(DM.ones(new long[] {8, 7}));
    assertThrows(CasadiException.class, () -> ca.mtimes(List.<MX>of()));
    MX D = ca.mtimes(List.of(A));
    assertEquals(4, D.size1());
    assertEquals(3, D.size2());
    D = ca.mtimes(List.of(A, B));
    assertEquals(4, D.size1());
    assertEquals(8, D.size2());
    D = ca.mtimes(List.of(A, B, C));
    assertEquals(4, D.size1());
    assertEquals(7, D.size2());
  }

  private static Sparsity parseSp(String art) {
    List<double[]> rows = new ArrayList<>();
    for (String r : art.strip().split("\n")) {
      r = r.strip();
      if (r.isEmpty()) continue;
      double[] v = new double[r.length()];
      for (int i = 0; i < r.length(); i++) v[i] = r.charAt(i) == 'X' ? 1 : 0;
      rows.add(v);
    }
    return ca.sparsify(new DM(rows.toArray(new double[0][]))).sparsity();
  }

  private static List<String> blasPlugins() {
    List<String> plugins = new ArrayList<>(List.of("reference"));
    for (String opt : new String[] {"classic", "blasfeo"}) {
      try {
        ca.load_blas(opt);
        plugins.add(opt);
      } catch (RuntimeException e) {
        // not available
      }
    }
    return plugins;
  }

  private static List<String> classicLibs() {
    List<String> r = new ArrayList<>();
    for (String t : CasadiMeta.lapack_libraries().split(";")) if (!t.isEmpty()) r.add(t);
    return r;
  }

  private static boolean compactibleMatch(List<Object> a, List<Object> b, List<Object> z, Sparsity zsp) {
    return (Boolean) a.get(0) && (Boolean) b.get(0) && (Boolean) z.get(0) && zsp.nnz() > 0
        && Arrays.equals((long[]) a.get(2), (long[]) b.get(1)) && Arrays.equals((long[]) a.get(1), (long[]) z.get(1))
        && Arrays.equals((long[]) b.get(2), (long[]) z.get(2));
  }

  private static String predictCategory(Sparsity A, Sparsity B, Sparsity Z) {
    if (A.is_dense() && B.is_dense() && Z.is_dense()) return "dense";
    if (A.is_dense() && Z.is_dense()) return "dense_sparse";
    if (compactibleMatch(A.is_compactible(), B.is_compactible(), Z.is_compactible(), Z)) return "dense";
    return "generic";
  }

  private static String generateCode(Function f, String name, Map<String, Object> opts, Path dir) {
    CodeGenerator cg = new CodeGenerator(name, opts);
    cg.add(f);
    cg.generate(dir.toString() + File.separator);
    return readFile(dir.resolve(name + ".c"));
  }

  @Test
  public void test_mtimes_blas_plugin() {
    Path dir = tmpDir();
    DM.rng(1);
    int m = 3, k = 4, n = 5;
    String[] aLabels = {"dense", "compactible", "sparse"};
    Sparsity[] aPats = {parseSp("XXXX\nXXXX\nXXXX"), parseSp("X.X.\n....\nX.X."), parseSp("X.XX\n.X.X\nXX..")};
    String[] bLabels = {"dense", "compactible", "sparse"};
    Sparsity[] bPats = {parseSp("XXXXX\nXXXXX\nXXXXX\nXXXXX"), parseSp(".X.X.\n.....\n.X.X.\n....."),
        parseSp("X.X.X\n.X..X\nXX.X.\n.X.XX")};
    DM AdenseVal = DM.rand(m, k);
    DM BdenseVal = DM.rand(k, n);
    List<String> plugins = blasPlugins();
    Map<String, String> denseKernel = Map.of("reference", "casadi_mtimes_dense(", "classic", "CASADI_BLAS_DGEMM(",
        "blasfeo", "blasfeo_blas_dgemm(");
    Map<String, List<String>> pluginExtralibs = Map.of("reference", List.of(), "classic", classicLibs(), "blasfeo",
        List.of("blasfeo"));
    MX AdenseSym = MX.sym("A", m, k);
    MX BdenseSym = MX.sym("B", k, n);
    MX ZdenseSym = MX.sym("Z", m, n);
    DM ZdenseVal = randDM(m, n);

    for (String name : plugins) {
      for (int ia = 0; ia < 3; ia++) {
        for (int ib = 0; ib < 3; ib++) {
          for (String opLabel : new String[] {"mtimes", "mac"}) {
            Sparsity aSp = aPats[ia], bSp = bPats[ib];
            String tag = opLabel + "/" + name + "/" + aLabels[ia] + "/" + bLabels[ib];
            MX A = ca.project(AdenseSym, aSp);
            MX B = ca.project(BdenseSym, bSp);
            DM Amask = ca.densify(new DM(aSp, 1.0));
            DM Bmask = ca.densify(new DM(bSp, 1.0));
            MX outF, outRef;
            List<MX> symIn;
            List<DM> valIn;
            Sparsity zSp;
            if (opLabel.equals("mtimes")) {
              outF = ca.densify(ca.mtimes(A, B, name));
              outRef = ca.mtimes(AdenseSym.times(new MX(Amask)), BdenseSym.times(new MX(Bmask)), "reference");
              symIn = List.of(AdenseSym, BdenseSym);
              valIn = List.of(AdenseVal, BdenseVal);
              zSp = ca.mtimes(new DM(aSp, 1.0), new DM(bSp, 1.0)).sparsity();
            } else {
              outF = ca.mac(A, B, ZdenseSym, name);
              outRef = ca.mac(AdenseSym.times(new MX(Amask)), BdenseSym.times(new MX(Bmask)), ZdenseSym, "reference");
              symIn = List.of(AdenseSym, BdenseSym, ZdenseSym);
              valIn = List.of(AdenseVal, BdenseVal, ZdenseVal);
              zSp = Sparsity.dense(m, n);
            }
            Function f = fun("f", symIn, List.of(outF));
            Function fRef = fun("f_ref", symIn, List.of(outRef));
            checkfunction(f, fRef, valIn);

            DM Adm = ca.project(AdenseVal, aSp);
            DM Bdm = ca.project(BdenseVal, bSp);
            DM dmOut = opLabel.equals("mtimes") ? ca.mtimes(Adm, Bdm, name) : ca.mac(Adm, Bdm, ZdenseVal, name);
            checkarray(ca.densify(dmOut), fRef.call(valIn).get(0), tag + " (DM-level)");

            String cat = predictCategory(aSp, bSp, zSp);
            String expected = cat.equals("dense") ? denseKernel.get(name)
                : cat.equals("dense_sparse") ? "casadi_mtimes_dense_sparse(" : "casadi_mtimes(";

            String code = generateCode(f, "f", new HashMap<>(), dir);
            assertTrue(code.contains(expected), tag + ": expected kernel " + expected + " in generated code");

            String fpath = dir.resolve("f.casadi").toString();
            f.save(fpath);
            Function fdeser = Function.load(fpath);
            code = generateCode(fdeser, "fdeser", new HashMap<>(), dir);
            assertTrue(code.contains(expected), tag + ": expected kernel " + expected + " in deserialized code");

            checkfunctionLight(f, fdeser, valIn);

            List<String> extralibs = pluginExtralibs.getOrDefault(name, List.of());
            if (name.equals("reference") || !extralibs.isEmpty()) {
              if (expected.equals(denseKernel.get(name)) || name.equals("reference"))
                checkCodegen(f, valIn, new HashMap<>(), extralibs);
            }
          }
        }
      }
    }
  }

  @Test
  public void test_l1_blas_plugin() {
    Path dir = tmpDir();
    DM.rng(2);
    List<String> plugins = blasPlugins();
    Map<String, List<String>> pluginExtralibs = Map.of("reference", List.of(), "classic", classicLibs(), "blasfeo",
        List.of("blasfeo"));
    Map<String, String> expectedKernel = new HashMap<>();
    expectedKernel.put("reference/dot", "casadi_real casadi_dot(");
    expectedKernel.put("reference/norm_1", "casadi_real casadi_norm_1(");
    expectedKernel.put("reference/norm_2", "casadi_real casadi_norm_2(");
    expectedKernel.put("reference/axpy", "void casadi_axpy(");
    expectedKernel.put("reference/copy", "memcpy(y, x");
    expectedKernel.put("classic/dot", "CASADI_BLAS_DDOT");
    expectedKernel.put("classic/norm_1", "CASADI_BLAS_DASUM");
    expectedKernel.put("classic/norm_2", "CASADI_BLAS_DNRM2");
    expectedKernel.put("classic/axpy", "CASADI_BLAS_DAXPY");
    expectedKernel.put("classic/copy", "memcpy(y, x");
    expectedKernel.put("blasfeo/dot", "blasfeo_blas_ddot");
    expectedKernel.put("blasfeo/norm_1", "casadi_real casadi_norm_1(");
    expectedKernel.put("blasfeo/norm_2", "casadi_real casadi_norm_2(");
    expectedKernel.put("blasfeo/axpy", "blasfeo_blas_daxpy");
    expectedKernel.put("blasfeo/copy", "memcpy(y, x");

    int n = 7;
    MX xSym = MX.sym("x", n);
    MX ySym = MX.sym("y", n);
    DM xVal = DM.rand(n);
    DM yVal = DM.rand(n);
    GlobalOptions.setDefaultBlas("reference");
    List<String> opNames = List.of("dot", "norm_1", "norm_2");
    Function fRef = new Function("f_ref", List.of(xSym, ySym),
        List.of(ca.dot(xSym, ySym), ca.norm_1(xSym), ca.norm_2(xSym)), List.of("x", "y"), opNames);
    List<DM> refOut = fRef.call(List.of(xVal, yVal));

    MX innerX = MX.sym("x", n);
    Function inner = fun("inner", List.of(innerX), List.of(ca.times(2.0, innerX)));
    Function fAxpyRef = inner.map(4, new boolean[] {false}, new boolean[] {true});
    List<DM> fAxpyInputs = List.of(ca.horzcat(List.of(xVal, xVal, xVal, xVal)));
    List<DM> fAxpyRefOut = fAxpyRef.call(fAxpyInputs);

    try {
      for (String name : plugins) {
        GlobalOptions.setDefaultBlas(name);
        String tagBase = "blas=" + name;
        Function f = new Function("f", List.of(xSym, ySym),
            List.of(ca.dot(xSym, ySym), ca.norm_1(xSym), ca.norm_2(xSym)), List.of("x", "y"), opNames);
        List<DM> got = f.call(List.of(xVal, yVal));
        for (int i = 0; i < opNames.size(); i++)
          checkarray(got.get(i), refOut.get(i), tagBase + " op=" + opNames.get(i) + " numerical agreement");

        Function fAxpy = inner.map(4, new boolean[] {false}, new boolean[] {true});
        List<DM> gotAxpy = fAxpy.call(fAxpyInputs);
        checkarray(gotAxpy.get(0), fAxpyRefOut.get(0), tagBase + " op=axpy numerical agreement");

        Map<String, Function> perOpCarriers = new LinkedHashMap<>();
        perOpCarriers.put("dot", fun("g_dot", List.of(xSym, ySym), List.of(ca.dot(xSym, ySym))));
        perOpCarriers.put("norm_1", fun("g_norm_1", List.of(xSym), List.of(ca.norm_1(xSym))));
        perOpCarriers.put("norm_2", fun("g_norm_2", List.of(xSym), List.of(ca.norm_2(xSym))));
        perOpCarriers.put("axpy", inner.map(4, new boolean[] {false}, new boolean[] {true}));
        perOpCarriers.put("copy", fun("g_copy", List.of(xSym), List.of(xSym)));
        for (Map.Entry<String, Function> e : perOpCarriers.entrySet()) {
          Map<String, Object> o = new HashMap<>();
          o.put("l1_blas", true);
          String code = generateCode(e.getValue(), "l1", o, dir);
          String exp = expectedKernel.get(name + "/" + e.getKey());
          assertTrue(code.contains(exp), tagBase + " op=" + e.getKey() + ": expected " + exp + " in generated code");
        }

        String fpath = dir.resolve("l1.casadi").toString();
        f.save(fpath);
        Function fdeser = Function.load(fpath);
        checkfunctionLight(f, fdeser, List.of(xVal, yVal));

        List<String> extralibs = pluginExtralibs.getOrDefault(name, List.of());
        if (name.equals("reference") || !extralibs.isEmpty()) {
          Map<String, Object> o = new HashMap<>();
          o.put("l1_blas", true);
          checkCodegen(f, List.of(xVal, yVal), o, extralibs);
          checkCodegen(fAxpy, fAxpyInputs, o, extralibs);
        }
      }
    } finally {
      GlobalOptions.setDefaultBlas("reference");
    }
  }

  @Test
  public void test_truth() {
    message("Truth values");
    assertThrows(RuntimeException.class, () -> MX.sym("x").nonzero());
    assertTrue(new MX(1).nonzero());
    assertFalse(new MX(0).nonzero());
    assertTrue(new MX(0.2).nonzero());
    assertTrue(new MX(-0.2).nonzero());
    assertThrows(RuntimeException.class, () -> new MX(dm(2.0, 3)).nonzero());
    assertThrows(RuntimeException.class, () -> new MX().nonzero());
  }

  @Test
  public void test_MXbool() {
    message("bool");
    MX xy = MX.sym("x", 2);
    MX x = xy.get(0);
    MX y = xy.get(1);
    Function f = fun("f", List.of(xy), List.of(ca.vertcat(List.of(ca.logic_and(x, y), ca.logic_or(x, y), ca.logic_not(x)))));
    for (int t1 = 0; t1 < 2; t1++) {
      for (int t2 = 0; t2 < 2; t2++) {
        boolean T1 = t1 != 0, T2 = t2 != 0;
        DM fOut = ev(f, (Object) new double[] {t1, t2});
        checkarray(fOut, dm(T1 && T2 ? 1 : 0, T1 || T2 ? 1 : 0, !T1 ? 1 : 0), "bool(" + t1 + "," + t2 + "): " + fOut);
      }
    }
  }

  @Test
  public void test_short_circuiting_codegen() {
    MX x = MX.sym("x", 5);
    MX y = ca.if_else(ca.gt(x, 0.0), 0.0, x);
    Function f = new Function("f", List.of(x), List.of(y), List.of("x"), List.of("y"));
    checkCodegen(f, List.of(dm(-200, 200, -100, 100, 0)));

    MX x8 = MX.sym("x", 8);
    MX y8 = MX.sym("y", 8);
    MX z = ca.logic_and(x8, y8);
    MX z2 = ca.logic_or(x8, y8);
    f = fun("f", List.of(x8, y8), List.of(z, z2));
    checkCodegen(f, List.of(dm(0, 1, 0, 1, 0, 1, 0, 1), dm(0, 0, 1, 1, 0, 0, 1, 1)));
  }

  @Test
  public void test_MXineq() {
    message("SX ineq");
    MX xy = MX.sym("x", 2);
    MX x = xy.get(0);
    MX y = xy.get(1);
    Function f = fun("f", List.of(xy),
        List.of(ca.vertcat(List.of(ca.lt(x, y), ca.le(x, y), ca.ge(x, y), ca.eq(x, y), ca.ne(x, y)))));
    double[] vals = {-10, 0.1, 0, 1, 10};
    for (double t1 : vals) {
      for (double t2 : vals) {
        DM fOut = ev(f, (Object) new double[] {t1, t2});
        checkarray(fOut, dm(t1 < t2 ? 1 : 0, t1 <= t2 ? 1 : 0, t1 >= t2 ? 1 : 0, t1 == t2 ? 1 : 0, t1 != t2 ? 1 : 0),
            "ineq(" + t1 + "," + t2 + ")");
      }
    }
  }

  @Test
  public void test_if_else_zero() {
    MX x = MX.sym("x");
    MX y = ca.if_else(x, 5.0, 0.0);
    Function f = fun("f", List.of(x), List.of(y));
    DM fOut = ev(f, 1);
    assertTrue(fOut.to_double() == 5, "if_else_zero " + fOut);
    fOut = ev(f, 0);
    assertTrue(fOut.to_double() == 0, "if_else_zero");
  }

  @Test
  public void test_if_else() {
    MX x = MX.sym("x");
    MX y = ca.if_else(x, 1.0, 2.0);
    Function f = fun("f", List.of(x), List.of(y));
    assertTrue(ev(f, 1).to_double() == 1, "if_else");
    assertTrue(ev(f, 0).to_double() == 2, "if_else");
  }

  private static boolean allFinite(DM r) {
    for (double[] row : r.full())
      for (double v : row) if (!Double.isFinite(v)) return false;
    return true;
  }

  @Test
  public void test_if_else_call_reverse_4387() {
    double[][] masks = {{0, 0}, {1, 0}, {0, 1}, {1, 1}};
    for (double[] maskValue : masks) {
      MX x = MX.sym("x", 3);
      MX mask = MX.sym("mask", 2);
      Function f = new Function("f", List.of(x), List.of(ca.vertcat(List.of(ca.sumsqr(x), ca.sum1(ca.sin(x))))),
          Map.of("never_inline", true));
      MX z = ca.if_else(ca.gt(mask, 0.0), (MX) f.invoke(x), 0.0);
      Function g = fun("g", List.of(x, mask), List.of(ca.gradient(ca.sumsqr(z), x)));
      double[] value = {1, 2, 3};
      double vv = 0, ss = 0;
      for (double v : value) {
        vv += v * v;
        ss += Math.sin(v);
      }
      double[] expected = new double[3];
      for (int i = 0; i < 3; i++)
        expected[i] = maskValue[0] * 4 * vv * value[i] + maskValue[1] * 2 * ss * Math.cos(value[i]);
      checkarray(ev(g, value, maskValue), dm(expected), "rectangular nonlinear " + Arrays.toString(maskValue));
    }

    double[][] nonbooleanMasks = {{0, 0}, {1, 0}, {0, 1}, {1, 1}, {-1, 0}, {0, -2}, {2, -3}};
    String[] names = {"square", "scalar input", "sparse mask"};
    DM[] as = {dm(new double[][] {{1, 1}, {1, -1}}), dm(1, 2), dm(new double[][] {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}})};
    Sparsity[] sps = {Sparsity.dense(2, 1), Sparsity.dense(2, 1), Sparsity.triplet(3, 1, new long[] {0, 2}, new long[] {0, 0})};
    double[][][] maskValuess = {nonbooleanMasks, nonbooleanMasks, {{0, 0, 0}, {1, 0, 0}, {0, 0, 1}, {1, 0, 1}}};
    for (int c = 0; c < 3; c++) {
      DM a = as[c];
      Sparsity sp = sps[c];
      for (double[] maskValue : maskValuess[c]) {
        for (int ndir = 1; ndir <= 2; ndir++) {
          String sub = names[c] + " mask=" + Arrays.toString(maskValue) + " ndir=" + ndir;
          MX x = MX.sym("x", a.size2());
          MX mask = MX.sym("mask", sp);
          Function f = new Function("f", List.of(x), List.of(ca.mtimes(a, x)), Map.of("never_inline", true));
          MX z = ca.if_else(mask, (MX) f.invoke(x), 0.0);
          double[] r2 = new double[(int) a.size1()];
          for (int i = 0; i < r2.length; i++) r2[i] = i + 2;
          List<DM> seeds = List.of(DM.ones(a.size1()), dm(r2)).subList(0, ndir);
          List<List<MX>> v = new ArrayList<>();
          for (DM seed : seeds) v.add(List.of(new MX(seed)));
          List<List<MX>> adj = ca.reverse(List.of(z), List.of(x), v);
          List<MX> cols = new ArrayList<>();
          for (List<MX> d : adj) cols.add(d.get(0));
          Function g = fun("g", List.of(x, mask), List.of(ca.horzcat(cols)));
          List<DM> exp = new ArrayList<>();
          for (DM seed : seeds) exp.add(ca.mtimes(a.T(), ca.ne(dm(maskValue), 0.0).times(seed)));
          checkarray(ev(g, DM.ones(a.size2()), ca.project(dm(maskValue), sp)), ca.horzcat(exp), sub);
        }
      }
    }

    Object[][] cases = {{true, new double[] {0}, new double[] {0, 0, 0}}, {true, new double[] {1}, new double[] {0.5, 1, 1}},
        {false, new double[] {0, 0}, new double[] {0, 0, 0}}, {false, new double[] {1, 0}, new double[] {0.5, 0, 0}}};
    for (Object[] cs : cases) {
      boolean scalarMask = (Boolean) cs[0];
      double[] maskValue = (double[]) cs[1], expected = (double[]) cs[2];
      MX x = MX.sym("x", 3);
      MX mask = MX.sym("mask", scalarMask ? 1 : 2);
      Function f = new Function("f", List.of(x), List.of(ca.vertcat(List.of(ca.sqrt(x.get(0)), x.get(1).plus(x.get(2))))),
          Map.of("never_inline", true));
      MX z = ca.if_else(mask, (MX) f.invoke(x), 0.0);
      Function g = fun("g", List.of(x, mask), List.of(ca.gradient(ca.sum1(z), x)));
      boolean active = maskValue[0] != 0;
      DM result = ev(g, new double[] {active ? 1 : 0, 2, 3}, maskValue);
      assertTrue(allFinite(result));
      checkarray(result, dm(expected), "singular derivative " + scalarMask + " " + Arrays.toString(maskValue));
    }
  }

  @Test
  public void test_if_else_reverse_partial_nan() {
    Assumptions.assumeTrue(KNOWN_BUGS, "known bug");
    for (boolean inline : new boolean[] {false, true}) {
      MX x = MX.sym("x", 3);
      MX mask = MX.sym("mask", 2);
      Function f = new Function("f", List.of(x), List.of(ca.vertcat(List.of(ca.sqrt(x.get(0)), x.get(1).plus(x.get(2))))),
          Map.of("always_inline", inline, "never_inline", !inline));
      MX z = ca.if_else(mask, (MX) f.invoke(x), 0.0);
      Function g = fun("g", List.of(x, mask), List.of(ca.gradient(ca.sum1(z), x)));
      DM result = ev(g, new double[] {0, 2, 3}, new double[] {0, 1});
      assertTrue(allFinite(result));
      checkarray(result, dm(0, 1, 1));
    }
  }

  @Test
  public void test_regression491() {
    message("regression #491");
    SX u = SX.sym("u");
    SX x = SX.sym("x");
    Function F = new Function("F", List.of(u, x), List.of(u.plus(ca.rdivide(1.0, x))));
    MX U = MX.sym("U");
    MX X = (MX) F.invoke(U, U);
    MX G = (MX) F.invoke(U, X);
    for (int kk = 0; kk < 2; kk++) {
      Function gfcn;
      if (kk == 0) gfcn = fun("gfcn", List.of(U), List.of(G)).expand("e_gfcn", Map.of("ad_weight", 1));
      else gfcn = new Function("gfcn", List.of(U), List.of(G), Map.of("ad_weight", 1));
      Function J = jacobianOld(gfcn, 0, 0);
      List<DM> JOut = J.call(List.of(new DM(1)));
      checkarray(JOut.get(0).to_double(), 1, 1e-9);
    }
  }

  // test_ticket: python-only (Python list + MX operator)

  @Test
  public void test_jacobian_tools() {
    message("jacobian");
    MX X = MX.sym("X");
    MX Y = ca.jacobian(X.power(2.0), X);
    Function f = fun("f", List.of(X), List.of(Y));
    checkarray(ev(f, 2.3).to_double(), 4.6, 1e-7);
  }

  @Test
  public void test_reshape() {
    message("reshape");
    MX X = MX.sym("X", 10);
    DM i = new DM(Sparsity.lower(3), dm(range(6)));
    i.print_dense();
    System.out.println(i.T().get_nz(false, new Slice()));
    MX T = X.get_nz(false, i);
    MX q = T.T().get_nz(false, new Slice()).power(2.0);
    Function f = fun("f", List.of(X), List.of(q));
    DM fOut = ev(f, (Object) range(10));
    checkarray(dm(0, 1, 9, 4, 16, 25), fOut);

    MX Y = MX.sym("Y", 10);
    Function ff = fun("ff", List.of(Y), f.call(List.of(Y), true));
    checkarray(dm(0, 1, 9, 4, 16, 25), ev(ff, (Object) range(10)));

    Function J = fun("J", List.of(X), List.of(ca.jacobian(q, X)));
    DM JOut = ev(J, (Object) range(10));
    DM i2 = ca.horzcat(List.of(ca.diag(dm(0, 2, 4, 6, 8, 10)), DM.zeros(6, 4)));
    i2.set(idx(2, 3), new Slice(), i2.get(idx(3, 2), new Slice()));
    checkarray(i2, JOut);

    q = T.T().get_nz(false, new Slice()).power(2.0);
    J = fun("J", List.of(X), List.of(ca.jacobian(q, X)));
    JOut = ev(J, (Object) range(10));
    i2 = ca.horzcat(List.of(ca.diag(dm(0, 2, 4, 6, 8, 10)), DM.zeros(6, 4)));
    i2.set(idx(2, 3), new Slice(), i2.get(idx(3, 2), new Slice()));
    checkarray(i2, JOut);
  }

  @Test
  public void test_vertcat() {
    message("vertcat");
    MX X = MX.sym("X", 10);
    MX T = ca.vertcat(List.of(X.get(4), X.get(2)));
    MX q = T.power(2.0);
    Function f = fun("f", List.of(X), List.of(q));
    List<DM> fOut = f.call(List.of(dm(range(10))));
    checkarray(dm(16, 4), fOut.get(0));

    MX Y = MX.sym("Y", 10);
    Function ff = fun("ff", List.of(Y), f.call(List.of(Y), true));
    checkarray(dm(16, 4), ev(ff, (Object) range(10)));
    Function J = fun("J", List.of(X), List.of(ca.jacobian(q, X)));
    DM JOut = ev(J, (Object) range(10));
    DM i = DM.zeros(2, 10);
    i.set(0, 4, 8.0);
    i.set(1, 2, 4.0);
    checkarray(i, JOut);
    q = T.power(2.0);
    J = fun("J", List.of(X), List.of(ca.jacobian(q, X)));
    JOut = ev(J, (Object) range(10));
    checkarray(i, JOut);
  }

  @Test
  public void test_blockcat() {
    MX x = MX.sym("x");
    MX y = ca.blockcat(List.of(List.of(x, ca.times(2.0, x)), List.of(ca.times(3.0, x), ca.times(4.0, x))));
    Function f = fun("f", List.of(x), List.of(y));
    checkarray(ev(f, 3), dm(new double[][] {{3, 6}, {9, 12}}));
  }

  @Test
  public void test_veccats() {
    MX x = MX.sym("x", 2);
    assertEquals(x.hash(), ca.vec(x).hash());
  }

  @Test
  public void test_constmxmtimes() {
    ca.times(0.1, MX.ones(2));
  }

  @Test
  public void test_is_regular() {
    assertTrue(new MX(dm(0, 1)).is_regular());
    assertFalse(new MX(dm(0, Double.POSITIVE_INFINITY)).is_regular());
    assertThrows(RuntimeException.class, () -> MX.sym("x", 2).is_regular());
  }

  @Test
  public void test_diagcat() {
    MX C = ca.diagcat(List.of(new MX(dm(new double[][] {{-1.4, -3.2}, {-3.2, -28}})),
        new MX(dm(new double[][] {{15, -12, 2.1}, {-12, 16, -3.8}, {2.1, -3.8, 15}})), new MX(1.8), new MX(-4.0)));
    DM r = dm(new double[][] {{-1.4, -3.2, 0, 0, 0, 0, 0}, {-3.2, -28, 0, 0, 0, 0, 0}, {0, 0, 15, -12, 2.1, 0, 0},
        {0, 0, -12, 16, -3.8, 0, 0}, {0, 0, 2.1, -3.8, 15, 0, 0}, {0, 0, 0, 0, 0, 1.8, 0}, {0, 0, 0, 0, 0, 0, -4}});
    r = ca.sparsify(r);
    Function f = fun("f", List.of(), List.of(C));
    checkarray(f.call(List.<DM>of()).get(0), r);

    int cx = 3, cy = 4;
    for (int nx = 2; nx >= 0; nx--) {
      for (int ny = 2; ny >= 0; ny--) {
        System.out.println(nx + " " + ny);
        DM rd = ca.diagcat(List.of(new DM(nx, ny), DM.ones(cx, cy), new DM(nx, ny)));
        assertEquals(cx + 2 * nx, rd.size1());
        assertEquals(cy + 2 * ny, rd.size2());
        SX rs = ca.diagcat(List.of(new SX(nx, ny), SX.ones(cx, cy), new SX(nx, ny)));
        assertEquals(cx + 2 * nx, rs.size1());
        assertEquals(cy + 2 * ny, rs.size2());
        MX rm = ca.diagcat(List.of(new MX(nx, ny), MX.ones(cx, cy), new MX(nx, ny)));
        assertEquals(cx + 2 * nx, rm.size1());
        assertEquals(cy + 2 * ny, rm.size2());
      }
    }
  }

  @Test
  public void test_tril2symm() {
    MX x = MX.sym("x", Sparsity.lower(3));
    Function f = fun("f", List.of(x), List.of(ca.tril2symm(x)));
    DM fOut = ev(f, new DM(f.sparsity_in(0), dm(range(6))));
    checkarray(fOut, dm(new double[][] {{0, 1, 2}, {1, 3, 4}, {2, 4, 5}}));
  }

  @Test
  public void test_sparsity_indexing() {
    message("sparsity");
    DM B_ = dm(new double[][] {{1, 2, 3, 4, 5}, {6, 7, 8, 9, 10}});
    MX B = MX.sym("B", 2, 5);
    DM A = ca.sparsify(dm(new double[][] {{1, 1, 0, 0, 0}, {0, 0, 1, 0, 0}}));
    Sparsity sp = A.sparsity();
    java.util.function.Function<MX, DM> meval = m -> ev(fun("f", List.of(B), List.of(m)), B_);

    checkarray(meval.apply(B.get(sp)), dm(new double[][] {{1, 2, 0, 0, 0}, {0, 0, 8, 0, 0}}), "sparsity indexing");

    MX Bmod = new MX(B);
    Bmod.set(sp, -4.0);
    checkarray(meval.apply(Bmod), dm(new double[][] {{-4, -4, 3, 4, 5}, {6, 7, -4, 9, 10}}), "sparsity indexing assignement");

    Bmod = new MX(B);
    Bmod.set(sp, ca.times(2.0, B));
    checkarray(meval.apply(Bmod), dm(new double[][] {{2, 4, 3, 4, 5}, {6, 7, 16, 9, 10}}), "Imatrix indexing assignement");

    assertThrows(RuntimeException.class, () -> B.get(Sparsity.dense(4, 4)));
  }

  @Test
  public void test_symvar() {
    MX a = MX.sym("a");
    MX b = MX.sym("b");
    MX c = MX.sym("c");
    MX e = ca.cos(a.times(b)).plus(c);
    List<MX> w = ca.symvar(e);
    assertEquals(3, w.size());
    if (GlobalOptions.getSimplificationOnTheFly()) {
      assertTrue(ca.is_equal(w.get(0), a));
      assertTrue(ca.is_equal(w.get(1), b));
      assertTrue(ca.is_equal(w.get(2), c));
    }
  }

  /** Elementwise algebra on SX, MX or DM by name, so one expression list serves all types */
  private static final class Alg {
    final Class<?> cls;

    Alg(Class<?> cls) {
      this.cls = cls;
    }

    Object c(Object v) {
      try {
        if (v instanceof Number) return cls.getConstructor(double.class).newInstance(((Number) v).doubleValue());
        if (v instanceof DM && cls != DM.class) return cls.getConstructor(DM.class).newInstance(v);
        return v;
      } catch (ReflectiveOperationException e) {
        throw new RuntimeException(e);
      }
    }

    Object f(String name, Object... a) {
      Object[] b = new Object[a.length];
      Class<?>[] t = new Class<?>[a.length];
      for (int i = 0; i < a.length; i++) {
        b[i] = c(a[i]);
        t[i] = cls;
      }
      try {
        return ca.class.getMethod(name, t).invoke(null, b);
      } catch (java.lang.reflect.InvocationTargetException e) {
        throw (RuntimeException) e.getCause();
      } catch (ReflectiveOperationException e) {
        throw new RuntimeException(e);
      }
    }

    Object add(Object a, Object b) { return f("plus", a, b); }
    Object sub(Object a, Object b) { return f("minus", a, b); }
    Object mul(Object a, Object b) { return f("times", a, b); }
    Object div(Object a, Object b) { return f("rdivide", a, b); }
    Object pow(Object a, Object b) { return f("power", a, b); }
    Object neg(Object a) {
      try {
        return cls.getMethod("unaryMinus").invoke(c(a));
      } catch (ReflectiveOperationException e) {
        throw new RuntimeException(e);
      }
    }
  }

  private interface Gen {
    Object apply(Alg o, Object x, Object y, Object z, Object X, Object Y);
  }

  private static final class Simp {
    final boolean refcount, onTheFly;
    final Gen genA, genB;

    Simp(boolean refcount, boolean onTheFly, Gen genA, Gen genB) {
      this.refcount = refcount;
      this.onTheFly = onTheFly;
      this.genA = genA;
      this.genB = genB;
    }
  }

  private static Simp s(boolean refcount, boolean onTheFly, Gen genA, Gen genB) {
    return new Simp(refcount, onTheFly, genA, genB);
  }

  /** Python str() of the expected value; a String stands for a Python number literal */
  private static String pyStr(Object o) {
    return o.toString();
  }

  private static List<Simp> simplificationCases(Object dx) {
    final double inf = Double.POSITIVE_INFINITY;
    // "X is ca.MX" is always False in mx.py: X is rebound to sin(A) before the list is built
    final boolean isMX = false;
    return List.of(
        s(false, true, (o, x, y, z, X, Y) -> o.add(o.mul(y, x), o.mul(y, o.sub(1, x))), (o, x, y, z, X, Y) -> y),
        s(false, isMX, (o, x, y, z, X, Y) -> o.sub(o.mul(2, x), x), (o, x, y, z, X, Y) -> x),
        s(false, true, (o, x, y, z, X, Y) -> o.sub(o.mul(4, x), o.mul(3, x)), (o, x, y, z, X, Y) -> x),
        s(false, true, (o, x, y, z, X, Y) -> o.add(o.mul(0.2, x), o.mul(0.8, x)), (o, x, y, z, X, Y) -> x),
        s(false, true, (o, x, y, z, X, Y) -> o.add(x, 0), (o, x, y, z, X, Y) -> x),
        s(false, isMX, (o, x, y, z, X, Y) -> o.add(0, x), (o, x, y, z, X, Y) -> x),
        s(false, true, (o, x, y, z, X, Y) -> o.add(x, o.neg(y)), (o, x, y, z, X, Y) -> o.sub(x, y)),
        s(false, isMX, (o, x, y, z, X, Y) -> o.add(o.neg(x), y), (o, x, y, z, X, Y) -> o.sub(y, x)),
        s(false, true, (o, x, y, z, X, Y) -> o.add(o.mul(0.5, x), o.mul(0.5, x)), (o, x, y, z, X, Y) -> x),
        s(false, true, (o, x, y, z, X, Y) -> o.add(o.div(x, 2), o.div(x, 2)), (o, x, y, z, X, Y) -> x),
        s(false, true, (o, x, y, z, X, Y) -> o.add(o.sub(x, y), y), (o, x, y, z, X, Y) -> x),
        s(false, true, (o, x, y, z, X, Y) -> o.add(y, o.sub(x, y)), (o, x, y, z, X, Y) -> x),
        s(false, true, (o, x, y, z, X, Y) -> o.add(o.pow(o.f("sin", x), 2), o.pow(o.f("cos", x), 2)), (o, x, y, z, X, Y) -> "1"),
        s(false, true, (o, x, y, z, X, Y) -> o.sub(x, 0), (o, x, y, z, X, Y) -> x),
        s(false, isMX, (o, x, y, z, X, Y) -> o.sub(0, x), (o, x, y, z, X, Y) -> o.neg(x)),
        s(false, true, (o, x, y, z, X, Y) -> o.sub(x, x), (o, x, y, z, X, Y) -> "0"),
        s(false, true, (o, x, y, z, X, Y) -> o.sub(x, o.neg(y)), (o, x, y, z, X, Y) -> o.add(x, y)),
        s(false, true, (o, x, y, z, X, Y) -> o.sub(o.add(x, y), y), (o, x, y, z, X, Y) -> x),
        s(false, true, (o, x, y, z, X, Y) -> o.sub(o.add(x, y), x), (o, x, y, z, X, Y) -> y),
        s(false, true, (o, x, y, z, X, Y) -> o.sub(x, o.add(y, x)), (o, x, y, z, X, Y) -> o.neg(y)),
        s(false, true, (o, x, y, z, X, Y) -> o.sub(x, o.add(x, y)), (o, x, y, z, X, Y) -> o.neg(y)),
        s(false, true, (o, x, y, z, X, Y) -> o.sub(o.neg(x), y), (o, x, y, z, X, Y) -> o.neg(o.add(x, y))),
        s(false, true, (o, x, y, z, X, Y) -> o.mul(x, x), (o, x, y, z, X, Y) -> o.pow(x, 2)),
        s(false, true, (o, x, y, z, X, Y) -> o.mul(x, 3), (o, x, y, z, X, Y) -> o.mul(3, x)),
        s(false, true, (o, x, y, z, X, Y) -> o.mul(x, 0), (o, x, y, z, X, Y) -> "0"),
        s(false, true, (o, x, y, z, X, Y) -> o.mul(0, x), (o, x, y, z, X, Y) -> "0"),
        s(false, true, (o, x, y, z, X, Y) -> o.mul(x, 1), (o, x, y, z, X, Y) -> x),
        s(false, isMX, (o, x, y, z, X, Y) -> o.mul(1, x), (o, x, y, z, X, Y) -> x),
        s(false, true, (o, x, y, z, X, Y) -> o.mul(x, -1), (o, x, y, z, X, Y) -> o.neg(x)),
        s(false, isMX, (o, x, y, z, X, Y) -> o.mul(-1, x), (o, x, y, z, X, Y) -> o.neg(x)),
        s(false, true, (o, x, y, z, X, Y) -> o.mul(x, o.div(1, y)), (o, x, y, z, X, Y) -> o.div(x, y)),
        s(false, isMX, (o, x, y, z, X, Y) -> o.mul(o.div(1, x), y), (o, x, y, z, X, Y) -> o.div(y, x)),
        s(false, true, (o, x, y, z, X, Y) -> o.mul(5, o.mul(0.2, x)), (o, x, y, z, X, Y) -> x),
        s(false, true, (o, x, y, z, X, Y) -> o.mul(0.5, o.mul(2, x)), (o, x, y, z, X, Y) -> x),
        s(false, true, (o, x, y, z, X, Y) -> o.mul(5, o.div(x, 5)), (o, x, y, z, X, Y) -> x),
        s(false, true, (o, x, y, z, X, Y) -> o.mul(o.div(2, x), x), (o, x, y, z, X, Y) -> "2"),
        s(false, true, (o, x, y, z, X, Y) -> o.div(o.mul(2, x), x), (o, x, y, z, X, Y) -> "2"),
        s(false, true, (o, x, y, z, X, Y) -> o.div(o.mul(5, x), x), (o, x, y, z, X, Y) -> "5"),
        s(false, true, (o, x, y, z, X, Y) -> o.mul(x, o.div(2, x)), (o, x, y, z, X, Y) -> "2"),
        s(false, isMX, (o, x, y, z, X, Y) -> o.mul(o.neg(x), y), (o, x, y, z, X, Y) -> o.neg(o.mul(x, y))),
        s(false, true, (o, x, y, z, X, Y) -> o.mul(x, o.neg(y)), (o, x, y, z, X, Y) -> o.neg(o.mul(x, y))),
        s(false, isMX, (o, x, y, z, X, Y) -> o.mul(2, o.mul(0.5, x)), (o, x, y, z, X, Y) -> x),
        s(false, isMX, (o, x, y, z, X, Y) -> o.mul(2, o.mul(x, 0.5)), (o, x, y, z, X, Y) -> x),
        s(false, isMX, (o, x, y, z, X, Y) -> o.mul(2, o.mul(0.5, X)), (o, x, y, z, X, Y) -> X),
        s(false, isMX, (o, x, y, z, X, Y) -> o.mul(2, o.mul(X, 0.5)), (o, x, y, z, X, Y) -> X),
        s(false, true, (o, x, y, z, X, Y) -> o.div(0, x), (o, x, y, z, X, Y) -> "0"),
        s(false, true, (o, x, y, z, X, Y) -> o.div(x, 1), (o, x, y, z, X, Y) -> x),
        s(false, true, (o, x, y, z, X, Y) -> o.div(x, -1), (o, x, y, z, X, Y) -> o.neg(x)),
        s(false, true, (o, x, y, z, X, Y) -> o.div(x, x), (o, x, y, z, X, Y) -> "1"),
        s(false, true, (o, x, y, z, X, Y) -> o.div(o.mul(2, x), 2), (o, x, y, z, X, Y) -> x),
        s(false, true, (o, x, y, z, X, Y) -> o.div(o.mul(x, y), x), (o, x, y, z, X, Y) -> y),
        s(false, isMX, (o, x, y, z, X, Y) -> o.div(o.mul(x, y), y), (o, x, y, z, X, Y) -> x),
        s(false, isMX, (o, x, y, z, X, Y) -> o.div(1, x), (o, x, y, z, X, Y) -> o.pow(x, -1)),
        s(false, true, (o, x, y, z, X, Y) -> o.div(x, o.div(1, y)), (o, x, y, z, X, Y) -> o.mul(x, y)),
        s(false, true, (o, x, y, z, X, Y) -> o.div(o.mul(2, x), o.mul(2, y)), (o, x, y, z, X, Y) -> o.div(x, y)),
        s(false, true, (o, x, y, z, X, Y) -> o.div(o.div(x, 5), 0.2), (o, x, y, z, X, Y) -> x),
        s(false, true, (o, x, y, z, X, Y) -> o.div(x, o.mul(2, x)), (o, x, y, z, X, Y) -> "0.5"),
        s(false, true, (o, x, y, z, X, Y) -> o.div(o.neg(x), x), (o, x, y, z, X, Y) -> "-1"),
        s(false, true, (o, x, y, z, X, Y) -> o.div(x, o.neg(x)), (o, x, y, z, X, Y) -> "-1"),
        s(false, true, (o, x, y, z, X, Y) -> o.div(o.neg(x), o.neg(x)), (o, x, y, z, X, Y) -> "1"),
        s(false, true, (o, x, y, z, X, Y) -> o.div(o.div(x, y), x), (o, x, y, z, X, Y) -> o.div(1, y)),
        s(false, isMX, (o, x, y, z, X, Y) -> o.div(o.neg(x), y), (o, x, y, z, X, Y) -> o.neg(o.div(x, y))),
        s(false, true, (o, x, y, z, X, Y) -> o.div(x, o.neg(y)), (o, x, y, z, X, Y) -> o.neg(o.div(x, y))),
        s(false, true, (o, x, y, z, X, Y) -> o.pow(x, 0), (o, x, y, z, X, Y) -> "1"),
        s(false, true, (o, x, y, z, X, Y) -> o.pow(x, 2), (o, x, y, z, X, Y) -> o.mul(x, x)),
        s(false, true, (o, x, y, z, X, Y) -> o.pow(x, 3), (o, x, y, z, X, Y) -> o.mul(x, o.mul(x, x))),
        s(false, true, (o, x, y, z, X, Y) -> o.pow(x, -3), (o, x, y, z, X, Y) -> o.div(1, o.pow(x, 3))),
        s(false, true, (o, x, y, z, X, Y) -> o.pow(x, 4), (o, x, y, z, X, Y) -> o.mul(o.pow(x, 2), o.pow(x, 2))),
        s(false, true, (o, x, y, z, X, Y) -> o.pow(x, 0.5), (o, x, y, z, X, Y) -> o.pow(x, 0.5)),
        s(false, false, (o, x, y, z, X, Y) -> o.pow(x, y), (o, x, y, z, X, Y) -> o.pow(x, y)),
        s(false, true, (o, x, y, z, X, Y) -> o.f("ge", o.pow(x, 2), 0), (o, x, y, z, X, Y) -> "1"),
        s(false, true, (o, x, y, z, X, Y) -> o.f("ge", o.mul(2, o.pow(x, 2)), o.pow(x, 2)), (o, x, y, z, X, Y) -> "1"),
        s(false, true, (o, x, y, z, X, Y) -> o.f("fmin", x, inf), (o, x, y, z, X, Y) -> x),
        s(false, true, (o, x, y, z, X, Y) -> o.f("fmin", inf, x), (o, x, y, z, X, Y) -> x),
        s(false, true, (o, x, y, z, X, Y) -> o.f("fmin", -inf, x), (o, x, y, z, X, Y) -> "-inf"),
        s(false, true, (o, x, y, z, X, Y) -> o.f("fmin", x, x), (o, x, y, z, X, Y) -> x),
        s(false, true, (o, x, y, z, X, Y) -> o.f("fmax", -inf, x), (o, x, y, z, X, Y) -> x),
        s(false, true, (o, x, y, z, X, Y) -> o.f("fmax", x, -inf), (o, x, y, z, X, Y) -> x),
        s(false, true, (o, x, y, z, X, Y) -> o.f("fmax", x, inf), (o, x, y, z, X, Y) -> "inf"),
        s(false, true, (o, x, y, z, X, Y) -> o.f("fmax", x, x), (o, x, y, z, X, Y) -> x),
        s(false, true, (o, x, y, z, X, Y) -> o.f("lt", o.pow(x, 2), 0), (o, x, y, z, X, Y) -> "0"),
        s(false, true, (o, x, y, z, X, Y) -> o.f("eq", x, x), (o, x, y, z, X, Y) -> "1"),
        s(false, true, (o, x, y, z, X, Y) -> o.f("ne", x, x), (o, x, y, z, X, Y) -> "0"),
        s(false, true, (o, x, y, z, X, Y) -> o.pow(o.f("sqrt", x), 2), (o, x, y, z, X, Y) -> x),
        s(false, true, (o, x, y, z, X, Y) -> o.pow(o.f("sqrt", X), 2), (o, x, y, z, X, Y) -> X),
        s(false, true, (o, x, y, z, X, Y) -> o.pow(o.neg(x), 2), (o, x, y, z, X, Y) -> o.pow(x, 2)),
        s(false, isMX, (o, x, y, z, X, Y) -> o.f("abs", o.f("abs", x)), (o, x, y, z, X, Y) -> o.f("abs", x)),
        s(false, isMX, (o, x, y, z, X, Y) -> o.f("abs", o.f("sqrt", x)), (o, x, y, z, X, Y) -> o.f("sqrt", x)),
        s(false, isMX, (o, x, y, z, X, Y) -> o.f("abs", o.f("exp", x)), (o, x, y, z, X, Y) -> o.f("exp", x)),
        s(false, isMX, (o, x, y, z, X, Y) -> o.f("log", o.f("exp", x)), (o, x, y, z, X, Y) -> x),
        s(false, isMX, (o, x, y, z, X, Y) -> o.div(1, o.div(1, x)), (o, x, y, z, X, Y) -> x),
        s(false, true, (o, x, y, z, X, Y) -> o.f("sqrt", o.pow(x, 2)), (o, x, y, z, X, Y) -> o.f("abs", x)),
        s(false, isMX, (o, x, y, z, X, Y) -> o.f("abs", o.neg(x)), (o, x, y, z, X, Y) -> o.f("abs", x)),
        s(false, true, (o, x, y, z, X, Y) -> o.pow(o.f("abs", x), 2), (o, x, y, z, X, Y) -> o.pow(x, 2)),
        s(false, isMX, (o, x, y, z, X, Y) -> o.f("cos", o.neg(x)), (o, x, y, z, X, Y) -> o.f("cos", x)),
        s(false, isMX, (o, x, y, z, X, Y) -> o.f("cos", o.f("abs", x)), (o, x, y, z, X, Y) -> o.f("cos", x)),
        s(false, false, (o, x, y, z, X, Y) -> o.neg(o.neg(x)), (o, x, y, z, X, Y) -> x),
        s(false, true, (o, x, y, z, X, Y) -> o.f("cosh", o.mul(x, 0)), (o, x, y, z, X, Y) -> o.add(o.mul(x, 0), 1)),
        s(false, true, (o, x, y, z, X, Y) -> o.div(x, 0.5), (o, x, y, z, X, Y) -> o.mul(2, x)),
        s(false, true, (o, x, y, z, X, Y) -> o.add(x, x), (o, x, y, z, X, Y) -> o.mul(2, x)),
        s(false, true, (o, x, y, z, X, Y) -> o.f("lt", x, x), (o, x, y, z, X, Y) -> "0"),
        s(false, true, (o, x, y, z, X, Y) -> o.add(o.mul(2, x), x), (o, x, y, z, X, Y) -> o.mul(3, x)),
        s(false, true, (o, x, y, z, X, Y) -> o.add(x, o.mul(2, x)), (o, x, y, z, X, Y) -> o.mul(3, x)),
        s(false, isMX, (o, x, y, z, X, Y) -> o.sub(o.mul(2, x), x), (o, x, y, z, X, Y) -> x),
        s(false, true, (o, x, y, z, X, Y) -> o.sub(x, o.mul(2, x)), (o, x, y, z, X, Y) -> o.neg(x)),
        s(true, false, (o, x, y, z, X, Y) -> o.neg(o.sub(x, y)), (o, x, y, z, X, Y) -> o.sub(y, x)),
        s(true, false, (o, x, y, z, X, Y) -> o.add(o.neg(o.sub(x, y)), o.pow(o.sub(x, y), 2)),
            (o, x, y, z, X, Y) -> o.add(o.neg(dx), o.pow(dx, 2))),
        s(true, false, (o, x, y, z, X, Y) -> o.mul(o.add(x, o.div(y, z)), z), (o, x, y, z, X, Y) -> o.add(o.mul(x, z), y)),
        s(true, false, (o, x, y, z, X, Y) -> o.mul(o.add(o.div(x, z), y), z), (o, x, y, z, X, Y) -> o.add(x, o.mul(y, z))),
        s(true, false, (o, x, y, z, X, Y) -> o.mul(z, o.add(o.div(x, z), y)), (o, x, y, z, X, Y) -> o.add(x, o.mul(z, y))),
        s(true, false, (o, x, y, z, X, Y) -> o.mul(z, o.add(x, o.div(y, z))), (o, x, y, z, X, Y) -> o.add(o.mul(z, x), y)),
        s(true, false, (o, x, y, z, X, Y) -> o.mul(o.sub(x, o.div(y, z)), z), (o, x, y, z, X, Y) -> o.sub(o.mul(x, z), y)),
        s(true, false, (o, x, y, z, X, Y) -> o.mul(o.sub(o.div(x, z), y), z), (o, x, y, z, X, Y) -> o.sub(x, o.mul(y, z))),
        s(true, false, (o, x, y, z, X, Y) -> o.mul(z, o.sub(o.div(x, z), y)), (o, x, y, z, X, Y) -> o.sub(x, o.mul(z, y))),
        s(true, false, (o, x, y, z, X, Y) -> o.mul(z, o.sub(x, o.div(y, z))), (o, x, y, z, X, Y) -> o.sub(o.mul(z, x), y)));
  }

  @Test
  public void test_simplifications() {
    for (Class<?> cls : List.of(SX.class, MX.class)) {
      Alg o = new Alg(cls);
      Alg od = new Alg(DM.class);
      Object a, b, c, A, B;
      if (cls == SX.class) {
        a = SX.sym("a");
        b = SX.sym("b");
        c = SX.sym("c");
        A = SX.sym("A", 2, 2);
        B = SX.sym("B", 2, 2);
      } else {
        a = MX.sym("a");
        b = MX.sym("b");
        c = MX.sym("c");
        A = MX.sym("A", 2, 2);
        B = MX.sym("B", 2, 2);
      }
      List<GenericMatrix> args = List.of((GenericMatrix) a, (GenericMatrix) b, (GenericMatrix) c, (GenericMatrix) A,
          (GenericMatrix) B);
      Object x = o.f("sin", a), y = o.f("sin", b), z = o.f("sin", c);
      Object X = o.f("sin", A), Y = o.f("sin", B);
      DM.rng(1);
      DM a0 = DM.rand(sparsityOf(a));
      DM b0 = DM.rand(sparsityOf(b));
      DM c0 = DM.rand(sparsityOf(c));
      DM A0 = DM.rand(sparsityOf(A));
      DM B0 = DM.rand(sparsityOf(B));
      Object x0 = ca.sin(a0), y0 = ca.sin(b0), z0 = ca.sin(c0), X0 = ca.sin(A0), Y0 = ca.sin(B0);
      Object dx = o.sub(x, y);
      int count = 1477;
      Object Alast = null;
      for (Simp cs : simplificationCases(dx)) {
        System.out.println(count);
        count++;
        Object Ae = cs.genA.apply(o, x, y, z, X, Y);
        Object Be = cs.genB.apply(o, x, y, z, X, Y);
        Alast = Ae;
        Function f = new Function("f", args, List.of((GenericMatrix) Ae));
        if (cs.refcount) {
          assertNotEquals(pyStr(Be), pyStr(Ae));
          System.out.println(Ae);
          f = f.transform();
          Function fref = new Function("f", args, List.of((GenericMatrix) Be));
          assertEquals(fref.call(args, true, false).toString(), f.call(args, true, false).toString());
        } else {
          assertEquals(pyStr(Be), pyStr(Ae));
        }
        if (cs.onTheFly) {
          GlobalOptions.setSimplificationOnTheFly(false);
          Object Ae2;
          try {
            Ae2 = cs.genA.apply(o, x, y, z, X, Y);
          } finally {
            GlobalOptions.setSimplificationOnTheFly(true);
          }
          assertNotEquals(pyStr(Be), pyStr(Ae2));
        }
        DM Aval = f.call(List.of(a0, b0, c0, A0, B0)).get(0);
        DM Aref = (DM) cs.genA.apply(od, x0, y0, z0, X0, Y0);
        System.out.println(Aval);
        System.out.println(Aref);
        checkarray(Aval, Aref);
      }
      System.out.println(new DM(Sparsity.upper(3), 0).is_nonnegative());
      System.out.println(new DM(Sparsity.upper(3), 0.5).is_half());
      System.out.println(new MX(new DM(Sparsity.upper(3), 0.5)).is_half());
      System.out.println(new MX(new DM(Sparsity.dense(3, 3), 0.5)).is_half());
      o.mul(2, o.mul(dm(new double[][] {{2, 3}, {6, 8}}), Alast));
    }
  }

  private static Sparsity sparsityOf(Object x) {
    if (x instanceof SX) return ((SX) x).sparsity();
    if (x instanceof MX) return ((MX) x).sparsity();
    return ((DM) x).sparsity();
  }

  @Test
  public void test_vertcat_empty() {
    Assumptions.assumeTrue(KNOWN_BUGS, "known bug");
    MX a = new MX(new DM(0, 2));
    MX v = ca.vertcat(List.of(a, a));
    assertEquals(0, v.size1());
    assertEquals(2, v.size2());
    a = new MX(new DM(2, 0));
    v = ca.vertcat(List.of(a, a));
    assertEquals(4, v.size1());
    assertEquals(0, v.size2());
  }

  @Test
  public void test_jacobian_empty() {
    MX x = MX.sym("x", 3);
    MX s = ca.jacobian(new DM(0, 0), x);
    assertEquals(0, s.size1());
    assertEquals(3, s.size2());
    s = ca.jacobian(x, MX.sym("x", 0, 4));
    assertEquals(3, s.size1());
    assertEquals(0, s.size2());
  }

  @Test
  public void test_mul_sparsity() {
    int N = 10;
    MX x = MX.sym("x", N, N);
    MX y = MX.sym("y", N, N);
    DM x_ = randDM(N, N);
    DM y_ = randDM(N, N);
    Sparsity filt = Sparsity.diag(N).plus(Sparsity.triplet(N, N, new long[] {1}, new long[] {3}));
    Function f = fun("f", List.of(x, y), List.of(ca.mtimes(x, y)));
    Function g = fun("g", List.of(x, y), List.of(ca.mac(x, y, MX.zeros(filt))));
    DM fOut = ev(f, x_, y_);
    DM gOut = ev(g, x_, y_);
    checkarray(DM.ones(filt), DM.ones(g.sparsity_out(0)));
    checkarray(fOut.get(filt), gOut);
  }

  @Test
  public void test_mul_zero_wrong() {
    assertThrows(CasadiException.class, () -> ca.mtimes(MX.sym("X", 4, 5), MX.zeros(3, 2)));
  }

  private static DM lowerRange(Function f) {
    return new DM(f.sparsity_in(0), dm(range(15)));
  }

  @Test
  public void test_vertsplit() {
    MX a = MX.sym("X", Sparsity.lower(5));
    List<MX> v = ca.vertsplit(a, new long[] {0, 2, 4, 5});
    Function f = fun("f", List.of(a), v);
    List<DM> V = f.call(List.of(lowerRange(f)));
    assertEquals(3, V.size());
    checkarray(V.get(0), dm(new double[][] {{0, 0, 0, 0, 0}, {1, 5, 0, 0, 0}}));
    checkarray(V.get(1), dm(new double[][] {{2, 6, 9, 0, 0}, {3, 7, 10, 12, 0}}));
    checkarray(V.get(2), dm(new double[][] {{4, 8, 11, 13, 14}}));

    v = ca.vertsplit(a);
    f = fun("f", List.of(a), v);
    V = f.call(List.of(lowerRange(f)));
    assertEquals(a.size1(), V.size());
    checkarray(V.get(0), dm(new double[][] {{0, 0, 0, 0, 0}}));
    checkarray(V.get(1), dm(new double[][] {{1, 5, 0, 0, 0}}));
    checkarray(V.get(2), dm(new double[][] {{2, 6, 9, 0, 0}}));
    checkarray(V.get(3), dm(new double[][] {{3, 7, 10, 12, 0}}));
    checkarray(V.get(4), dm(new double[][] {{4, 8, 11, 13, 14}}));

    v = ca.vertsplit(a, 2);
    f = fun("f", List.of(a), v);
    V = f.call(List.of(lowerRange(f)));
    assertEquals(3, V.size());
    checkarray(V.get(0), dm(new double[][] {{0, 0, 0, 0, 0}, {1, 5, 0, 0, 0}}));
    checkarray(V.get(1), dm(new double[][] {{2, 6, 9, 0, 0}, {3, 7, 10, 12, 0}}));
    checkarray(V.get(2), dm(new double[][] {{4, 8, 11, 13, 14}}));

    v = ca.vertsplit(a, new long[] {0, 0, 3, a.size1()});
    f = fun("f", List.of(a), v);
    List<Object> W = evs(f, lowerRange(f));
    assertEquals(3, v.size());
    assertEquals(0, v.get(0).size1());
    assertEquals(5, v.get(0).size2());
    checkarray((DM) W.get(1), dm(new double[][] {{0, 0, 0, 0, 0}, {1, 5, 0, 0, 0}, {2, 6, 9, 0, 0}}));
    checkarray((DM) W.get(2), dm(new double[][] {{3, 7, 10, 12, 0}, {4, 8, 11, 13, 14}}));
  }

  @Test
  public void test_horzsplit() {
    MX a = MX.sym("X", Sparsity.lower(5));
    List<MX> v = ca.horzsplit(a, new long[] {0, 2, 4, 5});
    Function f = fun("f", List.of(a), v);
    List<Object> V = evs(f, lowerRange(f));
    assertEquals(3, V.size());
    checkarray((DM) V.get(0), dm(new double[][] {{0, 0}, {1, 5}, {2, 6}, {3, 7}, {4, 8}}));
    checkarray((DM) V.get(1), dm(new double[][] {{0, 0}, {0, 0}, {9, 0}, {10, 12}, {11, 13}}));
    checkarray((DM) V.get(2), dm(new double[][] {{0}, {0}, {0}, {0}, {14}}));

    v = ca.horzsplit(a);
    f = fun("f", List.of(a), v);
    V = evs(f, lowerRange(f));
    assertEquals(a.size1(), V.size());
    checkarray((DM) V.get(0), dm(0, 1, 2, 3, 4));
    checkarray((DM) V.get(1), dm(0, 5, 6, 7, 8));
    checkarray((DM) V.get(2), dm(0, 0, 9, 10, 11));
    checkarray((DM) V.get(3), dm(0, 0, 0, 12, 13));
    checkarray((DM) V.get(4), dm(0, 0, 0, 0, 14));

    v = ca.horzsplit(a, 2);
    f = fun("f", List.of(a), v);
    V = evs(f, lowerRange(f));
    assertEquals(3, V.size());
    checkarray((DM) V.get(0), dm(new double[][] {{0, 0}, {1, 5}, {2, 6}, {3, 7}, {4, 8}}));
    checkarray((DM) V.get(1), dm(new double[][] {{0, 0}, {0, 0}, {9, 0}, {10, 12}, {11, 13}}));
    checkarray((DM) V.get(2), dm(new double[][] {{0}, {0}, {0}, {0}, {14}}));

    v = ca.horzsplit(a, new long[] {0, 0, 3, a.size2()});
    f = fun("f", List.of(a), v);
    V = evs(f, lowerRange(f));
    assertEquals(3, v.size());
    assertEquals(5, v.get(0).size1());
    assertEquals(0, v.get(0).size2());
    checkarray((DM) V.get(1), dm(new double[][] {{0, 0, 0}, {1, 5, 0}, {2, 6, 9}, {3, 7, 10}, {4, 8, 11}}));
    checkarray((DM) V.get(2), dm(new double[][] {{0, 0}, {0, 0}, {0, 0}, {12, 0}, {13, 14}}));
  }

  @Test
  public void test_blocksplit() {
    MX a = MX.sym("X", Sparsity.lower(5));
    List<List<MX>> v = ca.blocksplit(a, new long[] {0, 2, 4, 5}, new long[] {0, 1, 3, 5});
    List<Function> fs = new ArrayList<>();
    for (List<MX> vr : v) fs.add(fun("fs", List.of(a), vr));
    List<List<DM>> w = new ArrayList<>();
    for (int i = 0; i < 3; i++) w.add(fs.get(i).call(List.of(lowerRange(fs.get(i)))));
    checkarray(w.get(0).get(0), dm(0, 1));
    checkarray(w.get(0).get(1), dm(new double[][] {{0, 0}, {5, 0}}));
    checkarray(w.get(1).get(0), dm(2, 3));
    checkarray(ca.blockcat(w), lowerRange(fs.get(0)));
  }

  private static boolean is(Object x, String method, Object... a) {
    try {
      Class<?>[] t = new Class<?>[a.length];
      for (int i = 0; i < a.length; i++) t[i] = double.class;
      return (Boolean) x.getClass().getMethod(method, t).invoke(x, a);
    } catch (java.lang.reflect.InvocationTargetException e) {
      throw (RuntimeException) e.getCause();
    } catch (ReflectiveOperationException e) {
      throw new RuntimeException(e);
    }
  }

  @Test
  public void test_is_one_etc() {
    double inf = Double.POSITIVE_INFINITY;
    List<java.util.function.Function<DM, Object>> ops = List.of(d -> d, MX::new, SX::new);
    for (java.util.function.Function<DM, Object> op : ops) {
      assertTrue(is(op.apply(new DM(Sparsity.upper(3), 0)), "is_zero"));
      assertTrue(is(op.apply(new DM(Sparsity.dense(3, 4), 0)), "is_zero"));
      assertTrue(is(op.apply(new DM(3, 4)), "is_zero"));
      assertFalse(is(op.apply(dm(new double[][] {{0, 0}, {0, 7}})), "is_zero"));

      assertFalse(is(op.apply(new DM(Sparsity.upper(3), 1)), "is_one"));
      assertTrue(is(op.apply(new DM(Sparsity.dense(3, 4), 1)), "is_one"));
      assertFalse(is(op.apply(dm(new double[][] {{1, 1}, {1, 7}})), "is_one"));

      assertFalse(is(op.apply(new DM(Sparsity.upper(3), 0.5)), "is_half"));
      assertTrue(is(op.apply(new DM(Sparsity.dense(3, 4), 0.5)), "is_half"));
      assertFalse(is(op.apply(dm(new double[][] {{0.5, 0.5}, {0.5, 7}})), "is_half"));

      assertFalse(is(op.apply(new DM(Sparsity.upper(3), -1)), "is_minus_one"));
      assertTrue(is(op.apply(new DM(Sparsity.dense(3, 4), -1)), "is_minus_one"));
      assertFalse(is(op.apply(dm(new double[][] {{-1, -1}, {-1, 7}})), "is_minus_one"));

      assertFalse(is(op.apply(new DM(Sparsity.upper(3), inf)), "is_inf"));
      assertTrue(is(op.apply(new DM(Sparsity.dense(3, 4), inf)), "is_inf"));
      assertFalse(is(op.apply(dm(new double[][] {{inf, inf}, {inf, 7}})), "is_inf"));

      assertFalse(is(op.apply(new DM(Sparsity.upper(3), -inf)), "is_minus_inf"));
      assertTrue(is(op.apply(new DM(Sparsity.dense(3, 4), -inf)), "is_minus_inf"));
      assertFalse(is(op.apply(dm(new double[][] {{-inf, -inf}, {-inf, 7}})), "is_minus_inf"));

      assertTrue(is(op.apply(new DM(Sparsity.upper(3), 3)), "is_integer"));
      assertTrue(is(op.apply(new DM(Sparsity.dense(3, 4), 3)), "is_integer"));
      assertFalse(is(op.apply(dm(new double[][] {{3, 3}, {3, 0.7}})), "is_integer"));
      assertTrue(is(op.apply(dm(new double[][] {{3, 3}, {3, 7}})), "is_integer"));

      assertFalse(is(op.apply(new DM(Sparsity.upper(3), 7)), "is_value", 7.0));
      assertTrue(is(op.apply(new DM(Sparsity.dense(3, 4), 7)), "is_value", 7.0));
      assertFalse(is(op.apply(dm(new double[][] {{7, 7}, {7, 1}})), "is_value", 7.0));
      assertTrue(is(op.apply(new DM(Sparsity.upper(3), 0)), "is_value", 0.0));
      assertTrue(is(op.apply(new DM(Sparsity.dense(3, 4), 0)), "is_value", 0.0));
      assertTrue(is(op.apply(new DM(3, 4)), "is_value", 0.0));
      assertFalse(is(op.apply(dm(new double[][] {{0, 0}, {0, 7}})), "is_value", 0.0));

      assertFalse(is(op.apply(new DM(Sparsity.upper(3), 1)), "is_eye"));
      assertTrue(is(op.apply(new DM(Sparsity.diag(3), 1)), "is_eye"));
      assertFalse(is(op.apply(dm(new double[][] {{1, 0}, {0, 1}})), "is_eye"));

      assertFalse(is(op.apply(new DM(Sparsity.upper(3), -1)), "is_nonnegative"));
      assertTrue(is(op.apply(new DM(Sparsity.upper(3), 1)), "is_nonnegative"));
      assertTrue(is(op.apply(new DM(Sparsity.dense(3, 4), 1)), "is_nonnegative"));
      assertTrue(is(op.apply(dm(new double[][] {{1, 1}, {1, 7}})), "is_nonnegative"));
      assertFalse(is(op.apply(dm(new double[][] {{1, 1}, {1, -7}})), "is_nonnegative"));
    }
  }

  @Test
  public void test_mxnulloutput() {
    MX a = new MX(5, 0);
    MX b = MX.sym("x", 2);
    Function f = fun("f", List.of(b), List.of(a));
    MX c = (MX) f.invoke(b);
    assertEquals(5, c.size1());
    assertEquals(0, c.size2());
    c = f.call(List.of(b), true).get(0);
    assertEquals(5, c.size1());
    assertEquals(0, c.size2());

    a = new MX(0, 0);
    b = MX.sym("x", 2);
    f = fun("f", List.of(b), List.of(a));
    c = (MX) f.invoke(b);
    assertEquals(0, c.size1());
    assertEquals(0, c.size2());
    c = f.call(List.of(b), true).get(0);
    assertEquals(0, c.size1());
    assertEquals(0, c.size2());
  }

  @Test
  public void test_mxnull() {
    MX c = ca.mtimes(new MX(5, 0), new MX(0, 3));
    assertEquals(0, c.nnz());
    c = ca.mtimes(new MX(5, 3), new MX(3, 4));
    assertEquals(0, c.nnz());
  }

  @Test
  public void test_mxnullop() {
    MX c = new MX(0, 0);
    MX x = MX.sym("x", 2, 3);
    assertThrows(CasadiException.class, () -> x.plus(c));
    assertThrows(CasadiException.class, () -> x.div(c));
  }

  private static double minNz(DM c) {
    double m = Double.POSITIVE_INFINITY;
    for (double v : c.nonzeros()) m = Math.min(m, v);
    return m;
  }

  @Test
  public void test_MX_shapes() {
    Assumptions.assumeTrue(RUN_SLOW);
    message("MX unary operations");
    Sparsity sp43 = new Sparsity(4, 3, new long[] {0, 2, 2, 3}, new long[] {1, 2, 1});
    for (Sparsity sp : List.of(Sparsity.dense(0, 0), Sparsity.dense(0, 2), Sparsity.dense(2, 0), Sparsity.dense(1, 1),
        Sparsity.dense(2, 2), sp43)) {
      for (double v : new double[] {0, 1, 0.2}) {
        DM x_ = new DM(sp, v);
        MX xx = MX.sym("x", sp.size1(), sp.size2());
        MX x = xx.get(sp);
        for (Op1 op : pool) {
          if (op.nozero && v == 0) continue;
          MX r = op.cas.apply(x);
          Function f = fun("f", List.of(xx), List.of(r));
          DM rr = ev(f, v);
          System.out.println(r + " " + rr + " " + op.ref.apply(x_));
          checkarray(rr, op.ref.apply(x_));
          DM a = new DM(f.sparsity_out(0), 1);
          DM b = DM.ones(new DM(op.ref.apply(x_)).sparsity());
          DM c = b.minus(a);
          if (c.nnz() > 0) assertTrue(minNz(c) >= 0, sp + " " + v + " " + op.name);
        }
      }
    }
    List<Sparsity> sps = List.of(new Sparsity(1, 1), Sparsity.dense(1, 1), new Sparsity(3, 4), Sparsity.dense(3, 4), sp43.T());
    for (Sparsity sp : sps) {
      for (double v1 : new double[] {0, 1, 0.2, -0.2}) {
        DM x1_ = new DM(sp, v1);
        MX xx1 = MX.sym("x", sp.size1(), sp.size2());
        MX x1 = xx1.get(sp);
        SX xx1s = SX.sym("x", sp.size1(), sp.size2());
        xx1s.get(sp);
        for (Sparsity sp2 : sps) {
          for (double v2 : new double[] {0, 1, 0.2, -0.2}) {
            DM x2_ = new DM(sp2, v2);
            MX xx2 = MX.sym("x", sp2.size1(), sp2.size2());
            MX x2 = xx2.get(sp2);
            SX xx2s = SX.sym("x", sp2.size1(), sp2.size2());
            xx2s.get(sp2);
            for (Op2 op : matrixbinarypool) {
              if ((op.name.contains("mul") || op.name.contains("mtimes")) && (sp.numel() == 1 || sp2.numel() == 1)) continue;
              MX r = op.cas.apply(x1, x2);
              Function f = fun("f", List.of(xx1, xx2), List.of(r));
              List<DM> fIn = List.of(new DM(v1), new DM(v2));
              DM fOut = ev(f, v1, v2);
              Function g = fun("g", List.of(xx1, xx2), List.of(r));
              ev(g, v1, v2);
              String msg = sp + " " + sp2 + " " + v1 + " " + v2 + " " + op.name;
              checkarray(fOut, op.ref.apply(x1_, x2_), msg);
              if (!op.name.contains("mul")) {
                DM a = DM.ones(f.sparsity_out(0));
                DM b = DM.ones(g.sparsity_out(0));
                DM c = b.minus(a);
                if (c.nnz() > 0) assertTrue(minNz(c) >= 0, msg);
              }
              if (sp.nnz() > 0 && sp2.nnz() > 0 && v1 != 0 && v2 != 0) checkfunction(f, g, fIn);
            }
          }
        }
      }
    }
  }

  @Test
  public void test_MXConstant() {
    message("MX unary operations, constant");
    Sparsity sp43 = new Sparsity(4, 3, new long[] {0, 2, 2, 3}, new long[] {1, 2, 1});
    for (Sparsity sp : List.of(Sparsity.dense(0, 0), Sparsity.dense(0, 2), Sparsity.dense(2, 0), Sparsity.dense(1, 1),
        Sparsity.dense(2, 2), sp43)) {
      for (double v : new double[] {0, 1, 0.2}) {
        DM x_ = new DM(sp, v);
        MX x = new MX(sp, v);
        for (Op1 op : pool) {
          if (op.nozero && (v == 0 || !sp.is_dense())) continue;
          MX r = op.cas.apply(x);
          System.out.println(r);
          assertTrue(r.is_constant());
          checkarray(r.to_DM(), op.ref.apply(x_), x_ + " " + op.name);
          DM a = DM.ones(r.to_DM().sparsity());
          DM b = DM.ones(new DM(op.ref.apply(x_)).sparsity());
          DM c = b.minus(a);
          if (c.nnz() > 0) assertTrue(minNz(c) >= 0, sp + " " + v + " " + op.name);
        }
      }
    }
    List<Sparsity> sps = List.of(Sparsity.dense(1, 1), new Sparsity(1, 1), new Sparsity(3, 4), Sparsity.dense(3, 4), sp43.T());
    for (Sparsity sp : sps) {
      for (double v1 : new double[] {0, 1, 0.2, -0.2}) {
        DM x1_ = new DM(sp, v1);
        MX x1 = new MX(sp, v1);
        for (Sparsity sp2 : sps) {
          for (double v2 : new double[] {0, 1, 0.2, -0.2}) {
            DM x2_ = new DM(sp2, v2);
            MX x2 = new MX(sp2, v2);
            for (Op2 op : matrixbinarypool) {
              if ((op.name.contains("mul") || op.name.contains("mtimes")) && (sp.numel() == 1 || sp2.numel() == 1)) continue;
              MX r = op.cas.apply(x1, x2);
              Function f = fun("f", List.of(), List.of(r));
              List<DM> fOut = f.call(List.<DM>of());
              String msg = sp + " " + sp2 + " " + v1 + " " + v2 + " " + op.name;
              checkarray(fOut.get(0), op.ref.apply(x1_, x2_), msg);
              if (!op.name.contains("mul") && !op.name.contains("mtimes")) {
                DM a = DM.ones(f.sparsity_out(0));
                DM b = DM.ones(new DM(op.ref.apply(x1_, x2_)).sparsity());
                DM c = b.minus(a);
                if (c.nnz() > 0) assertTrue(minNz(c) >= 0, msg);
              }
            }
          }
        }
      }
    }
  }

  @Test
  public void test_graph_substitute() {
    MX x = MX.sym("X", 4, 4);
    MX y = MX.sym("Y", 4, 4);
    MX b = MX.sym("B", 4, 4);
    MX c = x.times(y);
    MX d = b.times(y);
    MX f = c.plus(d);
    MX C = MX.sym("C", 4, 4);
    f = ca.graph_substitute(f, List.of(c), List.of(C));
    Function F = fun("F", List.of(y, b, C), List.of(f));
    checkarray(ev(F, 1, 2, 3), ca.times(5.0, DM.ones(4, 4)));

    MX D = MX.sym("D", 4, 4);
    f = ca.graph_substitute(f, List.of(d), List.of(D));
    F = fun("F", List.of(D, C), List.of(f));
    checkarray(ev(F, 4, 5), ca.times(9.0, DM.ones(4, 4)));

    MX z = ca.vertcat(List.of(c, d, c.plus(d)));
    f = ca.graph_substitute(z, List.of(d), List.of(D));
    F = fun("F", List.of(x, y, b, D), List.of(f));
    Function Fref = fun("F", List.of(x, y, b, D), List.of(ca.vertcat(List.of(c, D, c.plus(D)))));
    checkfunctionLight(F, Fref, List.of(new DM(3), new DM(5), new DM(7), new DM(12)));
  }

  @Test
  public void test_graph_substitute_function() {
    MX x = MX.sym("x");
    MX y = MX.sym("y");
    Function f = new Function("forig", List.of(x, y), List.of(x.times(y)), Map.of("never_inline", true));
    Function g = new Function("fsubs", List.of(x, y), List.of(x.times(y)), Map.of("never_inline", true));
    MX sx = ca.sin(x);
    MX y2 = ca.times(2.0, y);
    MX a = (MX) f.invoke(sx, y2);
    MX x2 = x.power(2.0);
    MX z = ca.sqrt(x2.times(a)).plus(a.power(2.0));
    System.out.println(z);
    MX z2 = ca.graph_substitute(z, List.of(x2), List.of(x));
    System.out.println(z2);

    f = new Function("forig", List.of(x, y), List.of(x.times(y), x.minus(y), x.plus(y)), Map.of("never_inline", true));
    g = new Function("fsubs", List.of(x, y), List.of(x.times(y), x.minus(y), x.plus(y)), Map.of("never_inline", true));
    sx = ca.sin(x);
    y2 = ca.times(2.0, y);
    List<Object> abc = evs(f, sx, y2);
    a = (MX) abc.get(0);
    MX b = (MX) abc.get(1), c = (MX) abc.get(2);
    z = ca.sqrt(a).plus(c.power(2.0));
    System.out.println(z);
    z2 = ca.graph_substitute(z, List.of(), List.of());
    System.out.println(z2);

    z2 = ca.graph_substitute(z, List.of(a), List.of((MX) evs(g, sx, y2).get(0)));
    System.out.println(z2);
    assertTrue(z2.toString().contains("fsubs"));
    assertTrue(z2.toString().contains("forig"));
    Function Fref = fun("F_ref", List.of(x, y), List.of(z));
    Function F = fun("F", List.of(x, y), List.of(z2));
    checkfunctionLight(F, Fref, List.of(new DM(3.1), new DM(7.1)));

    z2 = ca.graph_substitute(z2, List.of(b), List.of(new MX(Double.NaN)));
    assertTrue(z2.toString().contains("fsubs"));
    assertTrue(z2.toString().contains("forig"));
    Fref = fun("F_ref", List.of(x, y), List.of(z));
    F = fun("F", List.of(x, y), List.of(z2));
    checkfunctionLight(F, Fref, List.of(new DM(3.1), new DM(7.1)));

    z2 = ca.graph_substitute(z2, List.of(c), List.of((MX) evs(g, sx, y2).get(2)));
    assertTrue(z2.toString().contains("fsubs"));
    assertFalse(z2.toString().contains("forig"));
    Fref = fun("F_ref", List.of(x, y), List.of(z));
    F = fun("F", List.of(x, y), List.of(z2));
    checkfunctionLight(F, Fref, List.of(new DM(3.1), new DM(7.1)));

    z2 = ca.graph_substitute(z, List.of(a.dep(0)), List.of(((MX) evs(g, sx, y2).get(0)).dep(0)));
    assertTrue(z2.toString().contains("fsubs"));
    assertFalse(z2.toString().contains("forig"));
    Fref = fun("F_ref", List.of(x, y), List.of(z));
    F = fun("F", List.of(x, y), List.of(z2));
    checkfunctionLight(F, Fref, List.of(new DM(3.1), new DM(7.1)));
  }

  @Test
  public void test_graph_substitute_call() {
    MX x = MX.sym("x");
    MX y = MX.sym("y");
    MX z = MX.sym("z");
    Function f = fun("f", List.of(x, y, z), List.of(x.power(2.0)));
    MX sx = ca.sin(x);
    MX x1 = x.plus(1.0);
    List<MX> res = ca.graph_substitute(List.of((MX) f.invoke(x1, y, z), ca.sin((MX) f.invoke(x1, y, z)),
        (MX) f.invoke(sx, y, z), ca.sin((MX) f.invoke(sx, y, z))), List.of(sx), List.of(ca.cos(x)));
    String s = res.toString();
    assertEquals(4, s.split("\\{0\\}", -1).length - 1);
  }

  @Test
  public void test_matrix_expand() {
    int n = 2;
    MX a = MX.sym("a", n, n);
    MX b = MX.sym("b", n, n);
    MX c = MX.sym("c", n, n);
    MX d = a.plus(b);
    MX e = d.times(c);
    assertEquals(6, ca.n_nodes(e));
    MX t0 = ca.matrix_expand(e);
    assertEquals(5, ca.n_nodes(t0));
    MX t1 = ca.matrix_expand(e, List.of(d));
    assertEquals(6, ca.n_nodes(t1));
    System.out.println(e + " " + t0 + " " + t1);
    List<DM> outs = new ArrayList<>();
    for (MX x : List.of(e, t0, t1)) {
      Function f = fun("f", List.of(a, b, c), List.of(x));
      outs.add(ev(f, 1.1, 2.2, 3.3));
      if (outs.size() > 1) checkarray(outs.get(0), outs.get(outs.size() - 1));
    }
    System.out.println(outs);
  }

  /** numpy.kron on dense arrays */
  private static double[][] npKron(double[][] a, double[][] b) {
    int ma = a.length, na = a[0].length, mb = b.length, nb = b[0].length;
    double[][] r = new double[ma * mb][na * nb];
    for (int i = 0; i < ma; i++)
      for (int j = 0; j < na; j++)
        for (int k = 0; k < mb; k++)
          for (int l = 0; l < nb; l++) r[i * mb + k][j * nb + l] = a[i][j] * b[k][l];
    return r;
  }

  @Test
  public void test_kron() {
    DM a = ca.sparsify(dm(new double[][] {{1, 0, 6}, {2, 7, 0}}));
    DM b = ca.sparsify(dm(new double[][] {{1, 0, 0}, {2, 3, 7}, {0, 0, 9}, {1, 12, 13}}));
    MX A = MX.sym("A", a.sparsity());
    MX B = MX.sym("B", b.sparsity());
    MX C = ca.kron(A, B);
    Function f = fun("f", List.of(A, B), List.of(C));
    List<DM> fIn = List.of(a, b);
    DM c_ = ev(f, a, b);
    assertEquals(a.size1() * b.size1(), c_.size1());
    assertEquals(a.size2() * b.size2(), c_.size2());
    assertEquals(a.nnz() * b.nnz(), c_.nnz());
    checkarray(c_, npKron(a.full(), b.full()));
    assertEquals(4, f.n_nodes());

    MX dA = MX.sym("dA", A.sparsity());
    MX dB = MX.sym("dB", B.sparsity());
    MX fwd = ca.jtimes(C, ca.vertcat(List.of(ca.vec(A), ca.vec(B))), ca.vertcat(List.of(ca.vec(dA), ca.vec(dB))));
    MX expected = ca.kron(dA, B).plus(ca.kron(A, dB));
    Function g = fun("g", List.of(A, B, dA, dB), List.of(fwd.minus(expected)));
    DM.rng(0);
    DM diff = ev(g, a, b, DM.rand(A.sparsity()), DM.rand(B.sparsity()));
    assertTrue(ca.norm_inf(diff).to_double() < 1e-14);

    checkCodegen(f, fIn);
    checkSerialize(f, fIn);
    checkfunction(f, f.expand(), fIn);

    SX As = SX.sym("A", a.sparsity());
    SX Bs = SX.sym("B", b.sparsity());
    SX Cs = ca.kron(As, Bs);
    Function fSX = new Function("f", List.of(As, Bs), List.of(Cs));
    checkfunction(f, fSX, fIn);
  }

  @Test
  public void test_kron_flavors() {
    DM.rng(11);
    Sparsity spD2 = Sparsity.dense(2, 3);
    Sparsity spD3 = Sparsity.dense(3, 2);
    Sparsity spS2 = ca.sparsify(dm(new double[][] {{1, 0, 1}, {0, 1, 0}})).sparsity();
    Sparsity spS3 = ca.sparsify(dm(new double[][] {{1, 0}, {1, 1}, {0, 1}})).sparsity();
    String[] labels = {"dense+dense", "dense+sparse", "sparse+dense", "sparse+sparse"};
    Sparsity[][] cases = {{spD2, spD3}, {spD2, spS3}, {spS2, spD3}, {spS2, spS3}};
    for (int i = 0; i < 4; i++) {
      String label = labels[i];
      MX A = MX.sym("A", cases[i][0]);
      MX B = MX.sym("B", cases[i][1]);
      MX K = ca.kron(A, B);
      Function f = fun("f", List.of(A, B), List.of(K));
      assertEquals(4, f.n_nodes(), label + ": graph not collapsed");
      DM aVal = DM.rand(cases[i][0]);
      DM bVal = DM.rand(cases[i][1]);
      checkarray(ev(f, aVal, bVal), new DM(npKron(aVal.full(), bVal.full())), label);
      checkCodegen(f, List.of(aVal, bVal));
      checkSerialize(f, List.of(aVal, bVal));
      checkfunction(f, f.expand(), List.of(aVal, bVal));
    }
  }

  private static double[][] refKronContract(double[][] m, double[][] x, boolean inner) {
    int mAmB = m.length, nAnB = m[0].length, xrow = x.length, xcol = x[0].length;
    if (inner) {
      int mB = xrow, nB = xcol, mA = mAmB / mB, nA = nAnB / nB;
      double[][] y = new double[mA][nA];
      for (int i = 0; i < mA; i++)
        for (int j = 0; j < nA; j++)
          for (int r = 0; r < mB; r++)
            for (int s = 0; s < nB; s++) y[i][j] += m[i * mB + r][j * nB + s] * x[r][s];
      return y;
    } else {
      int mA = xrow, nA = xcol, mB = mAmB / mA, nB = nAnB / nA;
      double[][] y = new double[mB][nB];
      for (int r = 0; r < mB; r++)
        for (int s = 0; s < nB; s++)
          for (int i = 0; i < mA; i++)
            for (int j = 0; j < nA; j++) y[r][s] += x[i][j] * m[i * mB + r][j * nB + s];
      return y;
    }
  }

  @Test
  public void test_kron_contract_flavors() {
    DM.rng(13);
    Sparsity spD2 = Sparsity.dense(2, 3);
    Sparsity spD3 = Sparsity.dense(3, 2);
    Sparsity spS2 = ca.sparsify(dm(new double[][] {{1, 0, 1}, {0, 1, 0}})).sparsity();
    Sparsity spS3 = ca.sparsify(dm(new double[][] {{1, 0}, {1, 1}, {0, 1}})).sparsity();
    String[] labels = {"dense_M+dense_X", "dense_M+sparse_X", "sparse_M+dense_X", "sparse_M+sparse_X"};
    Sparsity[][] combos = {{spD2, spD3}, {spD2, spS3}, {spS2, spD3}, {spS2, spS3}};
    for (int i = 0; i < 4; i++) {
      String label = labels[i];
      Sparsity spA = combos[i][0], spB = combos[i][1];
      Sparsity spM = ca.kron(spA, spB);
      MX Mi = MX.sym("M", spM);
      MX Xi = MX.sym("X", spB);
      MX Yi = ca.kron_contract(Mi, Xi, true);
      Function fi = fun("fi", List.of(Mi, Xi), List.of(Yi));
      DM mv = DM.rand(spM), xv = DM.rand(spB);
      double[][] ref = refKronContract(mv.full(), xv.full(), true);
      checkarray(new DM(ev(fi, mv, xv)), new DM(ref).times(DM.ones(Yi.sparsity())), label + " inner");
      checkCodegen(fi, List.of(mv, xv));
      checkSerialize(fi, List.of(mv, xv));
      checkfunction(fi, fi.expand(), List.of(mv, xv));

      MX Mo = MX.sym("M", spM);
      MX Xo = MX.sym("X", spA);
      MX Yo = ca.kron_contract(Mo, Xo, false);
      Function fo = fun("fo", List.of(Mo, Xo), List.of(Yo));
      DM mv2 = DM.rand(spM), xv2 = DM.rand(spA);
      double[][] ref2 = refKronContract(mv2.full(), xv2.full(), false);
      checkarray(new DM(ev(fo, mv2, xv2)), new DM(ref2).times(DM.ones(Yo.sparsity())), label + " outer");
      checkCodegen(fo, List.of(mv2, xv2));
      checkSerialize(fo, List.of(mv2, xv2));
      checkfunction(fo, fo.expand(), List.of(mv2, xv2));
    }
  }

  @Test
  public void test_kron_contract_empty_x() {
    Sparsity spXEmpty = new Sparsity(2, 2);
    Sparsity spM = Sparsity.dense(4, 4);
    for (boolean inner : new boolean[] {true, false}) {
      MX M = MX.sym("M", spM);
      MX X = MX.sym("X", spXEmpty);
      MX Y = ca.kron_contract(M, X, inner);
      assertEquals(2, Y.size1());
      assertEquals(2, Y.size2());
      assertEquals(0, Y.nnz());
      Function f = fun("f", List.of(M, X), List.of(Y));
      DM.rng(0);
      DM mv = DM.rand(spM), xv = new DM(spXEmpty);
      DM out = ev(f, mv, xv);
      assertEquals(2, out.size1());
      assertEquals(2, out.size2());
      assertEquals(0, out.nnz());
      checkCodegen(f, List.of(mv, xv));
      checkSerialize(f, List.of(mv, xv));
    }
  }

  @Test
  public void test_kron_contract() {
    DM.rng(7);
    DM a = ca.sparsify(dm(new double[][] {{1, 0, 6}, {2, 7, 0}}));
    DM b = ca.sparsify(dm(new double[][] {{1, 0, 0}, {2, 3, 7}, {0, 0, 9}, {1, 12, 13}}));
    MX A = MX.sym("A", a.sparsity());
    MX B = MX.sym("B", b.sparsity());
    MX K = ca.kron(A, B);
    MX M = MX.sym("M", K.sparsity());

    MX Yi = ca.kron_contract(M, B, true);
    assertEquals(A.size1(), Yi.size1());
    assertEquals(A.size2(), Yi.size2());
    MX Yo = ca.kron_contract(M, A, false);
    assertEquals(B.size1(), Yo.size1());
    assertEquals(B.size2(), Yo.size2());

    Function fInner = fun("fi", List.of(M, B), List.of(Yi));
    Function fOuter = fun("fo", List.of(M, A), List.of(Yo));

    DM mVal = DM.rand(K.sparsity());
    double[][] mNp = mVal.full(), aNp = a.full(), bNp = b.full();
    int mB1 = (int) b.size1(), nB1 = (int) b.size2(), mA1 = (int) a.size1(), nA1 = (int) a.size2();
    double[][] refInner = new double[mA1][nA1], refOuter = new double[mB1][nB1];
    for (int i = 0; i < mA1; i++)
      for (int j = 0; j < nA1; j++)
        for (int r = 0; r < mB1; r++)
          for (int s = 0; s < nB1; s++) refInner[i][j] += mNp[i * mB1 + r][j * nB1 + s] * bNp[r][s];
    for (int r = 0; r < mB1; r++)
      for (int s = 0; s < nB1; s++)
        for (int i = 0; i < mA1; i++)
          for (int j = 0; j < nA1; j++) refOuter[r][s] += aNp[i][j] * mNp[i * mB1 + r][j * nB1 + s];
    checkarray(ev(fInner, mVal, b), refInner);
    checkarray(ev(fOuter, mVal, a), refOuter);

    int[] sizes = {3, 5, 10, 20};
    List<Long> nodeCounts = new ArrayList<>();
    for (int n : sizes) {
      MX Ax = MX.sym("Ax", n, n);
      MX Bx = MX.sym("Bx", 3, 3);
      MX Kx = ca.kron(Ax, Bx);
      MX gA = ca.jacobian(ca.vec(Kx), ca.vec(Ax));
      nodeCounts.add(fun("f", List.of(Ax, Bx), List.of(gA)).n_nodes());
    }
    assertEquals(java.util.Collections.min(nodeCounts), java.util.Collections.max(nodeCounts),
        "kron jacobian graph must be scale-invariant: " + nodeCounts);

    List<Long> hessCounts = new ArrayList<>();
    for (int n : sizes) {
      MX Ax = MX.sym("Ax", n, n);
      MX Bx = MX.sym("Bx", 2, 2);
      MX Kx = ca.kron(Ax, Bx);
      MX obj = ca.times(0.5, ca.sumsqr(Kx));
      MX x = ca.vertcat(List.of(ca.vec(Ax), ca.vec(Bx)));
      MX H = (MX) ca.hessian(obj, x).get(0);
      hessCounts.add(fun("H", List.of(Ax, Bx), List.of(H)).n_nodes());
    }
    assertEquals(java.util.Collections.min(hessCounts), java.util.Collections.max(hessCounts),
        "kron hessian graph must be scale-invariant: " + hessCounts);

    MX A2 = MX.sym("A2", 3, 4);
    MX B2 = MX.sym("B2", 2, 5);
    MX M2 = MX.sym("M2", ca.kron(A2, B2).sparsity());
    List<List<MX>> cases = List.of(List.of(ca.kron_contract(M2, B2, true), M2, B2),
        List.of(ca.kron_contract(M2, A2, false), M2, A2));
    for (List<MX> cs : cases) {
      List<MX> args = cs.subList(1, 3);
      Function f = fun("f", args, List.of(cs.get(0)));
      DM.rng(11);
      List<DM> inputs = new ArrayList<>();
      for (MX arg : args) inputs.add(DM.rand(arg.sparsity()));
      checkfunctionLight(f, f.expand(), inputs);
    }

    MX A3 = MX.sym("A3", 3, 3);
    MX B3 = MX.sym("B3", 2, 2);
    MX M3 = MX.sym("M3", ca.kron(A3, B3).sparsity());
    MX Y = ca.kron_contract(M3, B3, true);
    MX dM = MX.sym("dM", M3.sparsity());
    MX dB = MX.sym("dB", B3.sparsity());
    MX fwd = ca.jtimes(Y, ca.vertcat(List.of(ca.vec(M3), ca.vec(B3))), ca.vertcat(List.of(ca.vec(dM), ca.vec(dB))));
    MX expected = ca.kron_contract(dM, B3, true).plus(ca.kron_contract(M3, dB, true));
    Function g = fun("g", List.of(M3, B3, dM, dB), List.of(fwd.minus(expected)));
    DM.rng(3);
    DM diff = ev(g, DM.rand(M3.sparsity()), DM.rand(B3.sparsity()), DM.rand(M3.sparsity()), DM.rand(B3.sparsity()));
    assertTrue(ca.norm_inf(diff).to_double() < 1e-14);

    DM Av = DM.rand(A2.sparsity()), Bv = DM.rand(B2.sparsity());
    DM Mv = DM.rand(M2.sparsity());
    Function fi = fun("fi", List.of(M2, B2), List.of(ca.kron_contract(M2, B2, true)));
    Function fo = fun("fo", List.of(M2, A2), List.of(ca.kron_contract(M2, A2, false)));
    checkCodegen(fi, List.of(Mv, Bv));
    checkCodegen(fo, List.of(Mv, Av));
    checkSerialize(fi, List.of(Mv, Bv));
    checkSerialize(fo, List.of(Mv, Av));
  }

  @Test
  public void test_project() {
    MX x = MX.sym("x", Sparsity.lower(3));
    MX y = ca.project(x, Sparsity.lower(3).T());
    Function f = fun("f", List.of(x), List.of(y));
    DM fIn = new DM(f.sparsity_in(0), dm(1, 2, 3, 4, 5, 6));
    DM fOut = ev(f, fIn);
    checkarray(fOut, dm(new double[][] {{1, 0, 0}, {0, 4, 0}, {0, 0, 6}}));
    checkarray(DM.ones(f.sparsity_out(0)), DM.ones(Sparsity.lower(3).T()));
    checkCodegen(f, List.of(dm(new double[][] {{1, 0, 0}, {0, 4, 0}, {0, 0, 6}})));
    checkSerialize(f, List.of(dm(new double[][] {{1, 0, 0}, {0, 4, 0}, {0, 0, 6}})));
  }
}

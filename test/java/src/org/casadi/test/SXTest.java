package org.casadi.test;

import static org.junit.jupiter.api.Assertions.assertArrayEquals;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;
import static org.junit.jupiter.api.Assertions.fail;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.function.BiFunction;
import java.util.function.BinaryOperator;
import java.util.function.DoubleBinaryOperator;
import java.util.function.DoubleUnaryOperator;
import java.util.function.UnaryOperator;
import java.util.regex.Pattern;

import org.casadi.CodeGenerator;
import org.casadi.DM;
import org.casadi.Function;
import org.casadi.GlobalOptions;
import org.casadi.MX;
import org.casadi.Operation;
import org.casadi.SX;
import org.casadi.Slice;
import org.casadi.Sparsity;
import org.casadi.ca;
import org.junit.jupiter.api.Assumptions;
import org.junit.jupiter.api.Disabled;
import org.junit.jupiter.api.Test;

/** Port of test/python/sx.py */
public class SXTest extends CasadiTestCase {

  private static final class UnaryOp {
    final UnaryOperator<SX> cas;
    final DoubleUnaryOperator num;
    final String name;
    final boolean nozero;

    UnaryOp(UnaryOperator<SX> cas, DoubleUnaryOperator num, String name, boolean nozero) {
      this.cas = cas;
      this.num = num;
      this.name = name;
      this.nozero = nozero;
    }
  }

  private static final class BinaryOp {
    final BinaryOperator<SX> cas;
    final BinaryOperator<DM> casDM;
    final BinaryOperator<double[][]> num;
    final String name;

    BinaryOp(BinaryOperator<SX> cas, BinaryOperator<DM> casDM, BinaryOperator<double[][]> num, String name) {
      this.cas = cas;
      this.casDM = casDM;
      this.num = num;
      this.name = name;
    }
  }

  private final List<UnaryOp> pool = new ArrayList<>();
  private final List<UnaryOp> Jpool = new ArrayList<>();
  private final List<BinaryOp> matrixbinarypool = new ArrayList<>();

  public SXTest() {
    pool.add(new UnaryOp(ca::sqrt, Math::sqrt, "sqrt", false));
    pool.add(new UnaryOp(ca::sin, Math::sin, "sin", false));
    pool.add(new UnaryOp(ca::cos, Math::cos, "cos", false));
    pool.add(new UnaryOp(ca::tan, Math::tan, "tan", false));
    pool.add(new UnaryOp(ca::abs, Math::abs, "fabs", false));
    pool.add(new UnaryOp(ca::sign, Math::signum, "sign", false));
    pool.add(new UnaryOp(ca::atan, Math::atan, "arctan", false));
    pool.add(new UnaryOp(ca::asin, Math::asin, "arcsin", false));
    pool.add(new UnaryOp(ca::acos, Math::acos, "arccos", false));
    pool.add(new UnaryOp(ca::exp, Math::exp, "exp", false));
    pool.add(new UnaryOp(ca::log, Math::log, "log", false));
    pool.add(new UnaryOp(x -> x.power(0), x -> Math.pow(x, 0), "x^0", true));
    pool.add(new UnaryOp(x -> x.power(1), x -> Math.pow(x, 1), "^1", false));
    pool.add(new UnaryOp(x -> x.power(-2), x -> Math.pow(x, -2), "^-2", true));
    pool.add(new UnaryOp(x -> x.power(0.3), x -> Math.pow(x, 0.3), "^0.3", false));
    pool.add(new UnaryOp(ca::floor, Math::floor, "floor", false));
    pool.add(new UnaryOp(ca::ceil, Math::ceil, "ceil", false));
    pool.add(new UnaryOp(ca::log1p, Math::log1p, "log1p", false));
    pool.add(new UnaryOp(ca::expm1, Math::expm1, "log1p", false));
    // num is the derivative; the reference Jacobian is its diagonal matrix
    Jpool.add(new UnaryOp(ca::sqrt, x -> 1 / (2.0 * Math.sqrt(x)), "sqrt", false));
    Jpool.add(new UnaryOp(ca::sin, Math::cos, "sin", false));
    Jpool.add(new UnaryOp(ca::abs, Math::signum, "fabs", false));
    Jpool.add(new UnaryOp(ca::sign, x -> x * 0, "fabs", false));
    Jpool.add(new UnaryOp(ca::cos, x -> -Math.sin(x), "cos", false));
    Jpool.add(new UnaryOp(ca::tan, x -> 1.0 / Math.pow(Math.cos(x), 2), "tan", false));
    Jpool.add(new UnaryOp(ca::atan, x -> 1.0 / (x * x + 1), "arctan", false));
    Jpool.add(new UnaryOp(ca::asin, x -> 1.0 / Math.sqrt(1 - x * x), "arcsin", false));
    Jpool.add(new UnaryOp(ca::acos, x -> -1.0 / Math.sqrt(1 - x * x), "arccos", false));
    Jpool.add(new UnaryOp(ca::exp, Math::exp, "exp", false));
    Jpool.add(new UnaryOp(ca::log, x -> 1.0 / x, "log", false));
    Jpool.add(new UnaryOp(x -> x.power(0), x -> 0, "x^0", false));
    Jpool.add(new UnaryOp(x -> x.power(1), x -> 1, "^1", false));
    Jpool.add(new UnaryOp(x -> x.power(-2), x -> -2.0 / Math.pow(x, 3), "^-2", false));
    Jpool.add(new UnaryOp(x -> x.power(0.3), x -> 0.3 / Math.pow(x, 0.7), "^0.3", false));
    Jpool.add(new UnaryOp(ca::log1p, x -> 1.0 / (1 + x), "log1p", false));
    Jpool.add(new UnaryOp(ca::expm1, Math::exp, "log1p", false));
    matrixbinarypool.add(new BinaryOp(SX::plus, DM::plus, (a, b) -> elementwise(a, b, Double::sum), "Matrix+Matrix"));
    matrixbinarypool.add(new BinaryOp(SX::minus, DM::minus, (a, b) -> elementwise(a, b, (u, v) -> u - v), "Matrix-Matrix"));
    matrixbinarypool.add(new BinaryOp(SX::times, DM::times, (a, b) -> elementwise(a, b, (u, v) -> u * v), "Matrix*Matrix"));
    matrixbinarypool.add(new BinaryOp(ca::fmax, ca::fmax, (a, b) -> elementwise(a, b, Math::max), "fmin"));
    matrixbinarypool.add(new BinaryOp(ca::fmin, ca::fmin, (a, b) -> elementwise(a, b, Math::min), "fmax"));
    matrixbinarypool.add(new BinaryOp(ca::hypot, ca::hypot, (a, b) -> elementwise(a, b, Math::hypot), "hypot"));
    matrixbinarypool.add(new BinaryOp((a, b) -> ca.mtimes(a, b.T()), (a, b) -> ca.mtimes(a, b.T()), SXTest::timesTransposed,
        "dot(Matrix,Matrix.T)"));
  }

  private static double[][] elementwise(double[][] a, double[][] b, DoubleBinaryOperator op) {
    double[][] r = new double[a.length][a[0].length];
    for (int i = 0; i < a.length; i++)
      for (int j = 0; j < a[0].length; j++) r[i][j] = op.applyAsDouble(a[i][j], b[i][j]);
    return r;
  }

  private static double[][] timesTransposed(double[][] a, double[][] b) {
    double[][] r = new double[a.length][b.length];
    for (int i = 0; i < a.length; i++)
      for (int j = 0; j < b.length; j++)
        for (int k = 0; k < a[0].length; k++) r[i][j] += a[i][k] * b[j][k];
    return r;
  }

  private static SX vcat(SX... a) {
    return ca.vertcat(List.of(a));
  }

  private static MX vcat(MX... a) {
    return ca.vertcat(List.of(a));
  }

  private static DM vcat(DM... a) {
    return ca.vertcat(List.of(a));
  }

  private static DM bools(int... v) {
    double[] d = new double[v.length];
    for (int i = 0; i < v.length; i++) d[i] = v[i];
    return new DM(d);
  }

  private static DM range(int n) {
    double[] d = new double[n];
    for (int i = 0; i < n; i++) d[i] = i;
    return new DM(d);
  }

  private static void assertBools(boolean[] actual, int... expected) {
    assertEquals(expected.length, actual.length);
    for (int i = 0; i < expected.length; i++) assertEquals(expected[i] != 0, actual[i], "entry " + i);
  }

  private static Function jacobianOld(Function f, int i, int j) {
    List<String> out = new ArrayList<>(List.of("jac:" + f.name_out(j) + ":" + f.name_in(i)));
    out.addAll(f.name_out());
    return f.factory(f.name() + "_jac", f.name_in(), out);
  }

  private static Function hessianOld(Function f, int i, int j) {
    List<String> out = new ArrayList<>(List.of("hess:" + f.name_out(j) + ":" + f.name_in(i) + ":" + f.name_in(i),
        "grad:" + f.name_out(j) + ":" + f.name_in(i)));
    out.addAll(f.name_out());
    return f.factory(f.name() + "_hess", f.name_in(), out);
  }

  private void evaluationCheck(List<SX> yt, DM yr, List<SX> x, List<DM> setx0, String name,
      BiFunction<Function, List<SX>, Function> fmod) {
    message(":" + name);
    Function f = new Function("f", x, yt);
    if (fmod != null) f = fmod.apply(f, x);
    List<DM> f_in = new ArrayList<>();
    for (int i = 0; i < f.n_in(); i++) f_in.add(i < setx0.size() ? setx0.get(i) : new DM(0));
    List<DM> f_out = f.call(f_in);
    checkarray(new DM(f_out.get(0).full()), yr, name);
  }

  private void numpyEvaluationCheckPool(List<UnaryOp> p, SX x, DM x0, String name,
      BiFunction<Function, List<SX>, Function> fmod, DM setx0, boolean excludeNozero, boolean jacobian) {
    for (UnaryOp op : p) {
      if (excludeNozero && op.nozero) continue;
      double[][] a = x0.full();
      int n = a.length * a[0].length;
      double[][] r = jacobian ? new double[n][n] : new double[a.length][a[0].length];
      for (int j = 0, k = 0; j < a[0].length; j++)
        for (int i = 0; i < a.length; i++, k++) {
          double v = op.num.applyAsDouble(a[i][j]);
          if (jacobian) r[k][k] = v;
          else r[i][j] = v;
        }
      evaluationCheck(List.of(op.cas.apply(x)), new DM(r), List.of(x), List.of(setx0 == null ? x0 : setx0),
          name + ":" + op.name, fmod);
    }
  }

  private void numpyEvaluationCheckBinaryPool(SX x, SX y, DM x0, DM y0, DM setx0, DM sety0, String name) {
    for (BinaryOp op : matrixbinarypool)
      evaluationCheck(List.of(op.cas.apply(x, y)), new DM(op.num.apply(x0.full(), y0.full())), List.of(x, y),
          List.of(setx0, sety0), name + ":" + op.name, null);
  }

  private static void checkSerialize(Function F, List<DM> inputs) {
    Function F2 = Function.deserialize(F.serialize(Map.of("debug", true)));
    List<DM> Fout = F.call(inputs), Fout2 = F2.call(inputs);
    for (int i = 0; i < F.n_out(); i++) checkarray(Fout.get(i), Fout2.get(i), 1e-16, "serialize");
  }

  /** As check_codegen in helpers.py: generate, compile with gcc, load as external, compare; only when slow */
  private static void checkCodegen(Function F, List<DM> inputs, Map<String, ?> opts, String std) {
    if (!RUN_SLOW) return;
    try {
      File dir = Files.createTempDirectory("casadi_java_codegen").toFile();
      String name = "codegen_" + Long.toHexString(Double.doubleToLongBits(Math.random()));
      CodeGenerator cg = new CodeGenerator(name, opts);
      cg.add(F);
      cg.generate(dir.getPath() + File.separator);
      String libdir = GlobalOptions.getCasadiPath(), inc = GlobalOptions.getCasadiIncludePath();
      String lib = new File(dir, name + ".so").getPath();
      List<String> cmd = new ArrayList<>(List.of("gcc", "-pedantic", "-std=" + std, "-fPIC", "-shared", "-Wall",
          "-Werror", "-Wextra", "-I" + inc, "-I" + inc + "/highs", "-Wno-dangling-pointer",
          "-Wno-unknown-warning-option", "-Wno-unknown-pragmas", "-Wno-long-long", "-Wno-unused-parameter", "-O3",
          new File(dir, name + ".c").getPath(), "-o", lib, "-L" + libdir, "-Wl,-rpath," + libdir));
      Process p = new ProcessBuilder(cmd).inheritIO().start();
      assertEquals(0, p.waitFor(), "compile " + cmd);
      Function F2 = ca.external(F.name(), lib);
      List<DM> Fout = F.call(inputs), Fout2 = F2.call(inputs);
      for (int i = 0; i < F.n_out(); i++) checkarray(Fout.get(i), Fout2.get(i), 1e-15, "codegen");
      deleteDir(dir);
    } catch (IOException | InterruptedException e) {
      throw new RuntimeException(e);
    }
  }

  private static void deleteDir(File dir) {
    File[] files = dir.listFiles();
    if (files != null)
      for (File f : files) f.delete();
    dir.delete();
  }

  private static void checkCodegen(Function F, List<DM> inputs, Map<String, ?> opts) {
    checkCodegen(F, inputs, opts, "c89");
  }

  private static void checkfunctionLight(Function f, Function g, List<DM> inputs) {
    List<DM> a = f.call(inputs), b = g.call(inputs);
    assertEquals(g.n_out(), f.n_out());
    for (int i = 0; i < a.size(); i++) checkarray(a.get(i), b.get(i), 1e-9, "output " + i);
  }

  private static void assertRaisesRegex(String regex, Runnable r) {
    try {
      r.run();
    } catch (RuntimeException e) {
      assertTrue(Pattern.compile(regex, Pattern.DOTALL).matcher(String.valueOf(e.getMessage())).find(),
          "expected /" + regex + "/ in: " + e.getMessage());
      return;
    }
    fail("expected an exception matching /" + regex + "/");
  }

  @Test
  public void test_equivalence() {
    SX x = SX.sym("x");
    SX y = SX.sym("y");
    SX[][] cases = {{ca.log1p(x), ca.log(x.plus(1))}, {ca.expm1(x), ca.exp(x).minus(1)},
        {ca.hypot(x, y), ca.sqrt(x.power(2).plus(y.power(2)))}};
    for (SX[] c : cases) {
      Function f = new Function("f", List.of(x, y), List.of(c[0]));
      Function equiv_f = new Function("equiv_f", List.of(x, y), List.of(c[1]));
      checkfunction(f, equiv_f, dms(new DM(1.3), new DM(1.7)));
    }
  }

  @Test
  public void test_scalarSX() {
    SX x = SX.sym("x");
    DM x0 = new DM(0.738);
    numpyEvaluationCheckPool(pool, x, x0, "scalarSX", null, null, false, false);
  }

  @Test
  public void test_gradient() {
    message("jacobian of SX**number");
    SX x = SX.sym("x");
    double x0 = 1;
    int p = 3;
    SX y = x.power(p);
    SX dx = ca.jacobian(y, x);
    double dxr = p;
    evaluationCheck(List.of(dx), new DM(dxr), List.of(x), List.of(new DM(x0)), "jacobian", null);
    dxr = 1;
    for (int i = 0; i < p; i++) {
      y = ca.jacobian(y, x);
      dxr = dxr * (p - i);
    }
    evaluationCheck(List.of(y), new DM(dxr), List.of(x), List.of(new DM(x0)), "recursive jacobian", null);
  }

  @Test
  public void test_gradient2() {
    message("jacobian of SX**SX");
    SX x = SX.sym("x");
    SX p = SX.sym("p");
    double x0 = 1;
    int p0 = 3;
    SX y = x.power(p);
    SX dx = ca.jacobian(y, x);
    double dxr = p0;
    evaluationCheck(List.of(dx), new DM(dxr), List.of(x, p), List.of(new DM(x0), new DM(p0)), "jacobian", null);
    dxr = 1;
    for (int i = 0; i < p0; i++) {
      y = ca.jacobian(y, x);
      dxr = dxr * (p0 - i);
    }
    evaluationCheck(List.of(y), new DM(dxr), List.of(x, p), List.of(new DM(x0), new DM(p0)), "jacobian", null);
  }

  @Test
  public void test_SXJacobian() {
    message("SX(1,1) unary operation, jacobian");
    SX x = SX.sym("x");
    DM x0 = dm(new double[][] {{0.738}});
    numpyEvaluationCheckPool(Jpool, x, x0, "SX unary operations, jacobian", (f, xx) -> jacobianOld(f, 0, 0), null,
        false, true);
  }

  @Test
  public void test_SXJac() {
    message("SX(1,1) unary operation, jac");
    SX x = SX.sym("x");
    DM x0 = dm(new double[][] {{0.738}});
    numpyEvaluationCheckPool(Jpool, x, x0, "SX unary operations, jac", (f, xx) -> {
      List<SX> y = f.call(xx);
      return new Function("J", xx, List.of(ca.jacobian(y.get(0), xx.get(0))));
    }, null, false, true);
  }

  @Test
  public void test_SXJacobians() {
    message("SX(3,1) unary operation, jacobian");
    SX x = SX.sym("x", 3);
    DM x0 = dm(0.738, 0.9, 0.3);
    numpyEvaluationCheckPool(Jpool, x, x0, "SX unary operations, jacobian", (f, xx) -> jacobianOld(f, 0, 0), null,
        false, true);
  }

  @Test
  public void test_SXJacobians2() {
    message("SX(1,3) unary operation, jacobian");
    SX x = SX.sym("x", 1, 3);
    DM x0 = dm(0.738, 0.9, 0.3);
    numpyEvaluationCheckPool(Jpool, x, x0, "SX unary operations, jacobian", (f, xx) -> jacobianOld(f, 0, 0), null,
        false, true);
  }

  @Test
  public void test_SX() {
    message("SX unary operations");
    SX x = SX.sym("x", 3, 2);
    DM x0 = dm(new double[][] {{0.738, 0.2}, {0.1, 0.39}, {0.99, 0.999999}});
    numpyEvaluationCheckPool(pool, x, x0, "SX", null, null, false, false);

    x = SX.sym("x", 3, 3);
    double[][] a = {{0.738, 0.2, 0.3}, {0.1, 0.39, -6}, {0.99, 0.999999, -12}};
    x0 = dm(a);
    double det = a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1]) - a[0][1] * (a[1][0] * a[2][2] - a[1][2] * a[2][0])
        + a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0]);
    double[][] inv = new double[3][3];
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++) {
        int r1 = (j + 1) % 3, r2 = (j + 2) % 3, c1 = (i + 1) % 3, c2 = (i + 2) % 3;
        inv[i][j] = (a[r1][c1] * a[r2][c2] - a[r1][c2] * a[r2][c1]) / det;
      }
    message(":det(SX)");
    evaluationCheck(List.of(new SX(ca.det(x))), new DM(det), List.of(x), List.of(x0), "det(SX)", null);
    message(":inv(SX)");
    evaluationCheck(List.of(ca.inv(x)), new DM(inv), List.of(x), List.of(x0), "inv(SX)", null);
  }

  @Test
  public void test_SXSparse() {
    message("SX unary operations, sparse");
    SX x = SX.sym("x");
    SX y = SX.sym("y");
    SX z = SX.sym("z");
    Sparsity sp = new Sparsity(4, 3, new long[] {0, 2, 2, 3}, new long[] {1, 2, 1});
    x = new SX(sp, vcat(x, y, z));
    DM x0 = new DM(sp, dm(0.738, 0.1, 0.99));
    numpyEvaluationCheckPool(pool, x, new DM(x0.full()), "SX", null, x0, true, false);
  }

  @Test
  public void test_SXbinary() {
    message("SX binary operations");
    SX x = SX.sym("x", 3, 2);
    SX y = SX.sym("x", 3, 2);
    DM x0 = dm(new double[][] {{0.738, 0.2}, {0.1, 0.39}, {0.99, 0.999999}});
    DM y0 = dm(new double[][] {{1.738, 0.6}, {0.7, 12}, {0, -6}});
    numpyEvaluationCheckBinaryPool(x, y, x0, y0, x0, y0, "SX");
    assertThrows(RuntimeException.class, () -> ca.mtimes(x, y));
  }

  @Test
  public void test_SXbinary_codegen() {
    Assumptions.assumeTrue(RUN_SLOW);
    message("SX binary operations");
    SX x = SX.sym("x", 4, 2);
    SX y = SX.sym("x", 4, 2);
    DM x0 = dm(new double[][] {{0.738, 0.2}, {0.1, 0.39}, {0.99, 0.999999}, {1, 2}});
    DM y0 = dm(new double[][] {{1.738, 0.6}, {0.7, 12}, {0, -6}, {1, 2}});
    for (BinaryOp op : matrixbinarypool) {
      Function f = new Function("f", List.of(x, y), List.of(op.cas.apply(x, y)));
      checkCodegen(f, dms(x0, y0), Map.of(), "c99");
      checkCodegen(f, dms(x0, y0), Map.of(), "c89");
    }
  }

  @Test
  public void test_SXbinary_diff() {
    message("SX binary operations");
    SX x = SX.sym("x", 4, 2);
    SX y = SX.sym("x", 4, 2);
    SX dx = SX.sym("x", 4, 2);
    SX dy = SX.sym("x", 4, 2);
    DM x0 = dm(new double[][] {{0.738, 0.2}, {0.1, 0.39}, {0.99, 0.999999}, {1, 2}});
    DM y0 = dm(new double[][] {{1.738, 0.6}, {0.7, 12}, {0, -6}, {1, 2}});
    for (BinaryOp op : matrixbinarypool) {
      String name = op.name;
      Function f = new Function("f", List.of(x, y), List.of(op.cas.apply(x, y)));
      f.disp(true);
      SX fxy = (SX) f.invoke(x, y);
      Function df = new Function("f", List.of(x, y, dx, dy), List.of(ca.jtimes(fxy, x, dx).plus(ca.jtimes(fxy, y, dy))));
      double eps = 1e-7;
      DM.rng(0);
      DM dx0 = DM.rand(x.sparsity());
      DM dy0 = name.equals("fmin") || name.equals("fmax") ? dx0 : DM.rand(y.sparsity());
      dy0 = name.equals("fmin") || name.equals("fmax") ? dx0 : DM.rand(y.sparsity());
      DM df_experiment = ((DM) f.invoke(x0.plus(dx0.times(eps)), y0.plus(dy0.times(eps)))).minus((DM) f.invoke(x0, y0))
          .div(eps);
      DM df0 = (DM) df.invoke(x0, y0, dx0, dy0);
      System.out.println(df_experiment + " " + df0);
      checkarray(df0, df_experiment, 0.5e-4, name);
    }
  }

  @Test
  public void test_DMbinary() {
    message("SX binary operations");
    DM x0 = dm(new double[][] {{0.738, 0.2}, {0.1, 0.39}, {0.99, 0.999999}});
    DM y0 = dm(new double[][] {{1.738, 0.6}, {0.7, 12}, {0, -6}});
    DM a = vcat(x0, y0);
    for (BinaryOp op : matrixbinarypool) {
      double[][] a0 = {{a.get(0).to_double()}}, a1 = {{a.get(1).to_double()}};
      checkarray(op.casDM.apply(a.get(0), a.get(1)), new DM(op.num.apply(a0, a1)), op.name);
    }
  }

  @Test
  public void test_SXbinarySparse() {
    message("SX binary operations");
    SX x = SX.sym("x");
    SX y = SX.sym("y");
    SX z = SX.sym("z");
    SX x2 = SX.sym("x2");
    SX y2 = SX.sym("y2");
    SX z2 = SX.sym("z2");
    Sparsity spx = new Sparsity(4, 3, new long[] {0, 2, 2, 3}, new long[] {1, 2, 1});
    Sparsity spy = new Sparsity(4, 3, new long[] {0, 2, 2, 3}, new long[] {0, 2, 3});
    SX xx = new SX(spx, vcat(x, y, z));
    SX yy = new SX(spy, vcat(x2, z2, y2));
    DM x0 = new DM(spx, dm(0.738, 0.1, 0.99));
    DM y0 = new DM(spy, dm(1.738, 0.7, -6));
    numpyEvaluationCheckBinaryPool(xx, yy, new DM(x0.full()), new DM(y0.full()), x0, y0, "SX");
    assertThrows(RuntimeException.class, () -> ca.mtimes(xx, yy));
  }

  private void slicingCheck(UnaryOperator<SX> ft, UnaryOperator<DM> fr, SX x, DM x0, DM setx0, String name) {
    evaluationCheck(List.of(ft.apply(x)), fr.apply(x0), List.of(x), List.of(setx0), name, null);
  }

  @Disabled("known bug in Python too: test refactoring, cf. #1436")
  @Test
  public void test_SXslicing() {
    message("SX slicing/indexing");
    SX x = SX.sym("x", 3, 2);
    DM x0 = dm(new double[][] {{0.738, 0.2}, {0.1, 0.39}, {0.99, 0.999999}});

    message(":dense");
    for (DM s0 : new DM[] {null}) {
      DM sx = s0 == null ? x0 : s0;
      slicingCheck(v -> new SX(v.get(0, 0)), v -> v.get(0, 0), x, x0, sx, "x[0,0]");
      slicingCheck(v -> new SX(v.get(1, 0)), v -> v.get(1, 0), x, x0, sx, "x[1,0]");
      slicingCheck(v -> new SX(v.get(0, 1)), v -> v.get(0, 1), x, x0, sx, "x[1,0]");
      slicingCheck(v -> new SX(v.get(0, -1)), v -> v.get(0, 1), x, x0, sx, "x[0,-1]");
      slicingCheck(v -> v.get(new Slice(), 0), v -> v.get(new Slice(), 0), x, x0, sx, "x[:,0]");
      slicingCheck(v -> v.get(new Slice(), 1), v -> v.get(new Slice(), 1), x, x0, sx, "x[:,1]");
      slicingCheck(v -> v.get(1, new Slice()), v -> v.get(1, new Slice()), x, x0, sx, "x[1,:]");
      slicingCheck(v -> v.get(0, new Slice()), v -> v.get(0, new Slice()), x, x0, sx, "x[0,:]");
      slicingCheck(v -> v.get(-1, new Slice()), v -> v.get(2, new Slice()), x, x0, sx, "x[-1,:]");
      slicingCheck(v -> v.get(new Slice(), -2), v -> v.get(new Slice(), 0), x, x0, sx, "x[:,-2]");
      slicingCheck(v -> v.get(new Slice(0, -2), new Slice(0, -1)), v -> v.get(new Slice(0, 1), new Slice(0, 1)), x, x0,
          sx, "x[0:-2,0:-1]");
      slicingCheck(v -> v.get(new Slice(0, 2), new Slice(0, 2)), v -> v.get(new Slice(0, 2), new Slice(0, 2)), x, x0,
          sx, "x[0:2,0:2]");
      slicingCheck(v -> v.get(dm(0, 1), new Slice(0, 2)), v -> v.get(new Slice(0, 2), new Slice(0, 2)), x, x0, sx,
          "x[[0,1],0:2]");
      slicingCheck(v -> v.get_nz(false, dm(0, 2, 3)), v -> dm(v.get(0, 0).to_double(), v.get(2, 0).to_double(),
          v.get(0, 1).to_double()), x, x0, sx, "x[[0,2,3]]");
      slicingCheck(v -> v.get_nz(false, new Slice(0, 2)), v -> dm(v.get(0, 0).to_double(), v.get(0, 1).to_double()), x,
          x0, sx, "x[0:2] on dense matrix");
      slicingCheck(v -> v.get_nz(false, 1), v -> v.get(0, 1), x, x0, sx, "x[1]");
      slicingCheck(v -> v.get_nz(false, -1), v -> v.get(2, 1), x, x0, sx, "x[-1]");
    }

    message(":sparse");
    Sparsity sp = new Sparsity(4, 3, new long[] {0, 2, 2, 3}, new long[] {1, 2, 1});
    x = new SX(sp, vcat(SX.sym("x"), SX.sym("y"), SX.sym("z")));
    DM sx0 = dm(0.738, 0.39, 0.99);
    x0 = new DM(new DM(sp, sx0).full());
    slicingCheck(v -> new SX(v.get(0, 0)), v -> v.get(0, 0), x, x0, sx0, "x[0,0]");
    slicingCheck(v -> new SX(v.get(0, 0)), v -> v.get(0, 0), x, x0, sx0, "x[0,0]");
    slicingCheck(v -> new SX(v.get(1, 0)), v -> v.get(1, 0), x, x0, sx0, "x[1,0]");
    slicingCheck(v -> new SX(v.get(0, 1)), v -> v.get(0, 1), x, x0, sx0, "x[1,0]");
    slicingCheck(v -> new SX(v.get(0, -1)), v -> v.get(0, 2), x, x0, sx0, "x[0,-1]");
    slicingCheck(v -> v.get(new Slice(), 0), v -> v.get(new Slice(), 0), x, x0, sx0, "x[:,0]");
    slicingCheck(v -> v.get(new Slice(), 1), v -> v.get(new Slice(), 1), x, x0, sx0, "x[:,1]");
    slicingCheck(v -> v.get(1, new Slice()), v -> v.get(1, new Slice()), x, x0, sx0, "x[1,:]");
    slicingCheck(v -> v.get(0, new Slice()), v -> v.get(0, new Slice()), x, x0, sx0, "x[0,:]");
    slicingCheck(v -> v.get(-1, new Slice()), v -> v.get(3, new Slice()), x, x0, sx0, "x[-1,:]");
    slicingCheck(v -> v.get(new Slice(), -2), v -> v.get(new Slice(), 1), x, x0, sx0, "x[:,-2]");
    slicingCheck(v -> v.get(new Slice(0, -2), new Slice(0, -1)), v -> v.get(new Slice(0, 2), new Slice(0, 2)), x, x0,
        sx0, "x[0:-2,0:-1]");
    slicingCheck(v -> v.get(new Slice(0, 2), new Slice(0, 2)), v -> v.get(new Slice(0, 2), new Slice(0, 2)), x, x0,
        sx0, "x[0:2,0:2]");
    slicingCheck(v -> v.get(dm(0, 1), new Slice(0, 2)), v -> v.get(new Slice(0, 2), new Slice(0, 2)), x, x0, sx0,
        "x[[0,1],0:2]");
    slicingCheck(v -> v.get_nz(false, dm(2, 1)), v -> dm(v.get(1, 2).to_double(), v.get(2, 0).to_double()), x, x0,
        sx0, "x[[2,1]]");
    slicingCheck(v -> v.get_nz(false, new Slice(0, 2)), v -> sx0.get(new Slice(0, 2)), x, x0, sx0,
        "x[0:2] on dense matrix");
    slicingCheck(v -> v.get_nz(false, 1), v -> sx0.get(1), x, x0, sx0, "x[1]");
    slicingCheck(v -> v.get_nz(false, -1), v -> sx0.get(2), x, x0, sx0, "x[-1]");
  }

  @Test
  public void test_SX1() {
    message("SXFunction evaluation");
    SX x = SX.sym("x");
    SX y = SX.sym("y");
    Function f = new Function("f", List.of(vcat(x, y)), List.of(vcat(x.plus(y), x.times(y), x.power(2).plus(y.power(3)))));
    DM L = dm(2, 3);
    List<DM> f_out = f.call(List.of(L));
    double[] z = f_out.get(0).nonzeros();
    double[] zr = {2 + 3, 2 * 3, 2 * 2 + 3 * 3 * 3};
    for (int i = 0; i < 3; i++) checkarray(z[i], zr[i], 1e-10);
    message("SXFunction jacobian evaluation");
    Function J = jacobianOld(f, 0, 0);
    List<DM> J_in = new ArrayList<>(List.of(L));
    for (int i = 1; i < J.n_in(); i++) J_in.add(new DM(0));
    List<DM> J_out = J.call(J_in);
    checkarray(J_out.get(0), new double[][] {{1, 1}, {3, 2}, {4, 27}});
  }

  @Test
  public void test_SX2() {
    message("SXFunction evalution 2");
    SX x = SX.sym("x");
    SX y = SX.sym("y");
    List<SX> f = List.of(ca.minus(3.0, ca.sin(x.times(x))).minus(y), ca.sqrt(y).times(x));
    if (GlobalOptions.getSimplificationOnTheFly()) {
      assertEquals("((3-sin(sq(x)))-y)", f.get(0).toString(), "SX representation is wrong " + f);
    } else {
      assertEquals("((3-sin((x*x)))-y)", f.get(0).toString(), "SX representation is wrong " + f);
    }
    assertEquals("(sqrt(y)*x)", f.get(1).toString(), "SX representation is wrong " + f);
    Function fcn = new Function("fcn", List.of(vcat(x, y)), List.of(ca.vertcat(f)));
    assertEquals("fcn:(i0[2])->(o0[2]) SXFunction", fcn.toString(), "SX representation is wrong");
    List<DM> fcn_out = fcn.call(List.of(dm(2, 3)));
    double[] res = fcn_out.get(0).nonzeros();
    checkarray(res[0], 3 - Math.sin(2 * 2) - 3, 1e-10);
    checkarray(res[1], Math.sqrt(3) * 2, 1e-10);
  }

  @Test
  public void test_const_folding_on_the_fly() {
    SX x = SX.sym("x");
    for (int a : new int[] {2, 7})
      for (int b : new int[] {2, 7}) {
        SX ref = ca.times((double) (a * b), x);
        assertEquals(ref.toString(), ca.times((double) a, ca.times((double) b, x)).toString());
        assertEquals(ref.toString(), ca.times((double) b, x).times(a).toString());
      }
  }

  @Test
  public void test_SX_func() {
    message("Function constructors");
    SX y = SX.sym("y", 2, 3);
    Function f = new Function("f", List.of(y), List.of(y));
    assertArrayEquals(new long[] {2, 3}, f.size_in(0), "Function constructors");
    assertArrayEquals(new long[] {2, 3}, f.size_out(0), "Function constructors");
    // Function("f", y, [y,y]) / Function("f", x0, [x0,x1]): python-only (non-list inputs do not compile in Java)
  }

  // test_evalfail: python-only (f.call(SX) with a non-list does not compile in Java)

  @Test
  public void test_SXconversion() {
    message("Conversions from and to SX");
    SX y = SX.sym("y");
    SX x = SX.sym("x", 3, 3);
    new SX(y);
    new SX(x);
    ca.det(x);
    DM yd = new DM(x);
    ca.det(yd);
  }

  @Test
  public void test_SXbool() {
    message("bool");
    SX x = SX.sym("x");
    SX y = SX.sym("y");
    Function f = new Function("f", List.of(vcat(x, y)), List.of(vcat(ca.logic_and(x, y), ca.logic_or(x, y), ca.logic_not(x))));
    for (int t1 : new int[] {0, 1})
      for (int t2 : new int[] {0, 1}) {
        boolean T1 = t1 != 0, T2 = t2 != 0;
        List<DM> f_out = f.call(List.of(dm(t1, t2)));
        checkarray(f_out.get(0), bools(T1 && T2 ? 1 : 0, T1 || T2 ? 1 : 0, !T1 ? 1 : 0),
            "bool(" + t1 + "," + t2 + "): " + f_out.get(0));
      }
  }

  @Test
  public void test_not_not() {
    SX x = SX.sym("x");
    SX y = SX.sym("y");
    Function f = new Function("f", List.of(x), List.of(ca.logic_not(ca.logic_not(x))));
    checkarray((DM) f.invoke(3), new DM(1));
    checkarray((DM) f.invoke(0), new DM(0));
    for (SX b : List.of(ca.lt(x, y), ca.le(x, y), ca.eq(x, y), ca.ne(x, y), ca.logic_not(x), ca.logic_and(x, y), ca.logic_or(x, y)))
      assertTrue(ca.is_equal(ca.logic_not(ca.logic_not(b)), b));
    assertFalse(ca.is_equal(ca.logic_not(ca.logic_not(x)), x));
  }

  @Test
  public void test_is_nonnegative_binary() {
    SX x = SX.sym("x");
    SX y = SX.sym("y");
    for (SX e : List.of(ca.lt(x, y), ca.le(x, y), ca.eq(x, y), ca.ne(x, y), ca.logic_and(x, y), ca.logic_or(x, y), ca.hypot(x, y),
        ca.sqrt(x)))
      assertTrue(ca.is_equal(ca.abs(e), e));
    for (SX e : List.of(x.plus(y), x.times(y), x.minus(y))) assertFalse(ca.is_equal(ca.abs(e), e));
  }

  @Test
  public void test_chained_inequality_kept() {
    SX y = SX.sym("y");
    // Python 0<=(y<=1) and (1>=y)>=0
    for (SX e : List.of(ca.le(0.0, ca.le(y, 1.0)), ca.ge(ca.ge(1.0, y), 0.0))) {
      assertFalse(e.is_constant());
      Function f = new Function("f", List.of(y), List.of(e)).transform(List.of(List.of("simplify", 0, "ref_count")));
      assertTrue(f.call(List.of(y), true, false).get(0).is_one());
    }
  }

  @Test
  public void test_and_zero() {
    SX xs = SX.sym("x");
    assertTrue(ca.logic_and(xs, 0.0).is_zero());
    assertTrue(ca.logic_and(0.0, xs).is_zero());
    MX xm = MX.sym("x");
    assertTrue(ca.logic_and(xm, 0.0).is_zero());
    assertTrue(ca.logic_and(0.0, xm).is_zero());
    assertTrue(ca.logic_and(new SX(0), new SX(0)).is_zero());
  }

  @Test
  public void test_SXineq() {
    message("SX ineq");
    SX x = SX.sym("x");
    SX y = SX.sym("y");
    Function f = new Function("f", List.of(vcat(x, y)),
        List.of(vcat(ca.lt(x, y), ca.le(x, y), ca.ge(x, y), ca.eq(x, y), ca.ne(x, y))));
    double[] vals = {-10, 0.1, 0, 1, 10};
    for (double T1 : vals)
      for (double T2 : vals) {
        List<DM> f_out = f.call(List.of(dm(T1, T2)));
        checkarray(f_out.get(0), bools(T1 < T2 ? 1 : 0, T1 <= T2 ? 1 : 0, T1 >= T2 ? 1 : 0, T1 == T2 ? 1 : 0,
            T1 != T2 ? 1 : 0), "ineq(" + T1 + "," + T2 + ")");
      }
  }

  @Test
  public void test_SX_func2() {
    message("SXmatrix typemaps constructors");
    ca.simplify(SX.sym("x"));
    message(":number");
    DM i = ca.transpose(ca.transpose(2.3));
    assertEquals(1, i.size1(), "shape mismatch");
    assertEquals(1, i.size2(), "shape mismatch");
    new SX(2.3).is_empty();
    message(":SX");
    SX s = SX.sym("x");
    SX j = ca.transpose(ca.transpose(s));
    assertEquals(1, j.size1(), "shape mismatch");
    assertEquals(1, j.size2(), "shape mismatch");
    new SX(s).is_empty();
  }

  @Test
  public void test_SX_func3() {
    message("vector(SXmatrix) typemaps constructors");
    SX y = SX.sym("y");
    SX x = SX.sym("x", 3, 1);
    vcat(x, x);
    vcat(y, y);
    vcat(x, new SX());
  }

  @Test
  public void test_eval() {
    message("Function eval");
    SX x = SX.sym("x", 2, 2);
    SX y = SX.sym("y", 2, 2);
    Function f = new Function("f", List.of(x, y), List.of(x.times(y)));
    f.invoke(x, y);
  }

  @Test
  public void test_symbolcheck() {
    message("Check if non-symbolic inputs are caught");
    assertThrows(RuntimeException.class, () -> new Function("f", List.of(new SX(0)), List.of(SX.sym("x"))));
  }

  @Test
  public void test_sparseconstr() {
    message("Check sparsity constructors");
    checkarray(new DM(DM.ones(Sparsity.lower(3)).full()), new double[][] {{1, 0, 0}, {1, 1, 0}, {1, 1, 1}});
    checkarray(new DM(DM.ones(Sparsity.diag(3)).full()), new double[][] {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}});
  }

  @Test
  public void test_subsassignment() {
    message("Check subscripted assignment");
    double[][] xn = new double[3][4];
    for (double[] r : xn)
      for (int j = 0; j < 4; j++) r[j] = random.nextDouble();
    DM x = new DM(xn);

    DM y = new DM(7, 8);
    double[][] z = new double[7][8];
    y.set(0, 0, 12);
    z[0][0] = 12;
    checkarray(y, z);
    int[] cols = {2, 4, 5, 6};
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 4; j++) z[1 + i][cols[j]] = xn[i][j];
    y.set(new Slice(1, 4), dm(2, 4, 5, 6), x);
    checkarray(y, z);

    long[] kl = {2, 4, 5, 8};
    DM kld = dm(2, 4, 5, 8);
    y.set_nz(1.0, false, kld);
    Sparsity s = y.sparsity();
    for (long k : kl) z[(int) s.row()[(int) k]][(int) s.get_col()[(int) k]] = 1.0;
    checkarray(y, z);
    y.set_nz(kld, false, kld);
    int cnt = 0;
    for (long k : kl) {
      z[(int) s.row()[(int) k]][(int) s.get_col()[(int) k]] = kl[cnt];
      cnt += 1;
    }
    checkarray(y, z);
  }

  @Test
  public void test_substitute() {
    Assumptions.assumeTrue(GlobalOptions.getSimplificationOnTheFly());
    message("Basic symbolic algebra: substitute");
    SX x = SX.sym("x");
    SX y = SX.sym("y");
    SX z = ca.cos(x).times(y);
    assertTrue(ca.depends_on(z, y));
    assertTrue(ca.depends_on(z, x));
    SX w = ca.substitute(z, x, 0.0);
    assertTrue(w.is_symbolic());
    assertTrue(ca.depends_on(w, y));
    assertFalse(ca.depends_on(w, x));
    assertTrue(ca.is_equal(w, y));
    SX r = w.minus(y);
    assertFalse(r.is_symbolic());
    assertTrue(r.is_zero());
    assertEquals(0, r.to_double());
    assertEquals(0, r.to_double());
    SX y2 = SX.sym("y", 2);
    y2 = ca.substitute(y2.plus(6), y2, 0.0);
    assertEquals(6, y2.get(0).to_int());
    assertEquals(6, y2.get(1).to_int());
  }

  @Test
  public void test_primitivefunctions() {
    message("Primitive functions");
    SX x = SX.sym("x");
    double[] nums = {-2, -1.5, -1, -0.5, -0.25, 0, 0.25, 0.5, 1, 1.5, 2};
    primitiveTest(x, ca::sign, "sign", nums, new double[] {-1, -1, -1, -1, -1, 0, 1, 1, 1, 1, 1});
    primitiveTest(x, ca::heaviside, "heaviside", nums, new double[] {0, 0, 0, 0, 0, 0.5, 1, 1, 1, 1, 1});
    primitiveTest(x, ca::ramp, "ramp", nums, new double[] {0, 0, 0, 0, 0, 0, 0.25, 0.50, 1, 1.5, 2});
    primitiveTest(x, ca::rectangle, "rectangle", nums, new double[] {0, 0, 0, 0.5, 1, 1, 1, 0.5, 0, 0, 0});
    primitiveTest(x, ca::triangle, "triangle", nums, new double[] {0, 0, 0, 0.5, 0.75, 1, 0.75, 0.5, 0, 0, 0});
  }

  private void primitiveTest(SX x, UnaryOperator<SX> fun, String comment, double[] nums, double[] reference) {
    message(":" + comment);
    Function f = new Function("f", List.of(x), List.of(fun.apply(x)));
    for (int i = 0; i < nums.length; i++) {
      List<DM> f_out = f.call(List.of(new DM(nums[i])));
      double r = f_out.get(0).get(0).to_double();
      assertTrue(reference[i] == r, reference[i] + " vs " + r);
    }
  }

  @Test
  public void test_taylor() {
    message("univariate taylor expansion");
    SX x = SX.sym("x");
    if (GlobalOptions.getSimplificationOnTheFly()) assertTrue(ca.is_equal(ca.taylor(ca.sin(x), x), x));
    double a_ = 0.13, x_ = 0.15;
    SX a = SX.sym("a");
    BiFunction<SX, Double, Void> test = (e, r) -> {
      Function f = new Function("f", List.of(x, a), List.of(e));
      List<DM> f_out = f.call(dms(new DM(x_), new DM(a_)));
      checkarray(f_out.get(0).get(0).to_double(), r, 1e-10);
      return null;
    };
    double d = x_ - a_;
    test.apply(ca.taylor(ca.sin(x), x, a, 0), Math.sin(a_));
    test.apply(ca.taylor(ca.sin(x), x, a, 1), Math.sin(a_) + Math.cos(a_) * d);
    test.apply(ca.taylor(ca.sin(x), x, a, 2), Math.sin(a_) + Math.cos(a_) * d - (Math.sin(a_) * d * d) / 2.0);
    test.apply(ca.taylor(ca.sin(x), x, a, 3),
        Math.sin(a_) + Math.cos(a_) * d - (Math.sin(a_) * d * d) / 2.0 - (Math.cos(a_) * d * d * d) / 6.0);

    SX M = ca.blockcat(List.of(List.of(a.times(ca.sin(x)), a.times(ca.cos(x))), List.of(ca.exp(a.times(x)),
        a.times(x.power(2))), List.of(ca.cos(x), new SX(0))));
    Function f = new Function("f", List.of(x, a), List.of(ca.taylor(M, x)));
    List<DM> f_out = f.call(dms(new DM(x_), new DM(a_)));
    checkarray(f_out.get(0), new double[][] {{x_ * a_, a_}, {1 + a_ * x_, 0}, {1, 0}}, 1e-10, "taylor on dense matrices");
  }

  private static void checkarray(DM a, double[][] b, double tol, String msg) {
    checkarray(a, new DM(b), tol, msg);
  }

  @Test
  public void test_null() {
    message("Function null");
    SX x = SX.sym("x");

    Function f = new Function("f", List.of(x), List.of(x.power(2), new SX()));
    List<DM> f_out = f.call(dms(new DM(0)));
    assertTrue(f_out.get(1).is_empty());

    f = new Function("f", List.of(x, new SX()), List.of(x.power(2), new SX()));
    f_out = f.call(dms(new DM(0), new DM(0)));
    assertTrue(f_out.get(1).is_empty());
    f_out = f.call(dms(new DM(0), new DM(0)));

    List<SX> r = f.call(List.of(x, new SX()));
    assertTrue(r.get(1).is_empty());
    r = f.call(List.of(x, new SX()));
    assertTrue(r.get(1).is_empty());
    r = f.call(List.of(x, new SX(0, 1)));
    assertTrue(r.get(1).is_empty());
    r = f.call(List.of(x, new SX(1, 0)));
    assertTrue(r.get(1).is_empty());
  }

  @Test
  public void test_mtaylor() {
    message("multivariate taylor expansions");
    SX x = SX.sym("x");
    SX y = SX.sym("y");
    SX a = SX.sym("a");
    SX b = SX.sym("b");
    double a_ = 0.13, x_ = 0.15, b_ = 0.73, y_ = 0.75;
    SX xy = vcat(x, y), ab = vcat(a, b);
    BiFunction<SX, Double, Void> test = (e, r) -> {
      Function f = new Function("f", List.of(xy, ab), List.of(e));
      List<DM> f_out = f.call(dms(dm(x_, y_), dm(a_, b_)));
      checkarray(f_out.get(0).get(0).to_double(), r, 1e-10);
      return null;
    };
    double s = Math.sin(a_ + b_), c = Math.cos(a_ + b_), dx = x_ - a_, dy = y_ - b_;
    SX e = ca.sin(x.plus(y));
    test.apply(ca.mtaylor(e, xy, ab, 0), s);
    test.apply(ca.mtaylor(e, xy, ab, 1), s + (c * dx + c * dy));
    double sol2 = s + (c * dx + c * dy) - (s * dx * dx + 2 * s * dy * dx + s * dy * dy) / 2;
    test.apply(ca.mtaylor(e, xy, ab, 2), sol2);
    double sol3 = sol2 - (c * Math.pow(dx, 3) + 3 * c * dy * dx * dx + 3 * c * dy * dy * dx + c * Math.pow(dy, 3)) / 6;
    test.apply(ca.mtaylor(e, xy, ab, 3), sol3);
    double sol3b = (-2 * s * dx * dy - s * dx * dx) / 2 + c * dy - (c * Math.pow(dx, 3)) / 6 + c * dx + s;
    test.apply(ca.mtaylor(e, xy, ab, 3, new long[] {1, 2}), sol3b);
    test.apply(ca.mtaylor(e, xy, vcat(new SX(0), new SX(0)), 4, new long[] {1, 2}),
        (-3 * x_ * x_ * y_ - Math.pow(x_, 3)) / 6 + y_ + x_);
  }

  @Test
  public void test_issue107() {
    message("Regression test for issue 107: +=");
    SX x = SX.sym("x");
    SX y = SX.sym("y");
    SX z = x;
    z = z.plus(y);
    assertTrue(x.is_symbolic());
    assertFalse(z.is_symbolic());

    x = SX.sym("x");
    y = SX.sym("y");
    z = x;
    z = z.plus(y);
    assertTrue(x.is_symbolic());
    assertFalse(z.is_symbolic());
  }

  @Test
  public void test_evalchecking() {
    SX x = SX.sym("x", 1, 5);
    SX y = SX.sym("y", 1, 3);
    SX z = SX.sym("z", 5, 1);
    SX q = SX.sym("z", 1, 6);
    Function f = new Function("f", List.of(x), List.of(x.power(2)));
    assertThrows(RuntimeException.class, () -> f.invoke(y));
    assertThrows(RuntimeException.class, () -> f.invoke(q));
    f.invoke(z);
  }

  @Test
  public void test_indexinglimits() {
    message("Limits of indexing");
    SX y = SX.sym("y", 3);
    assertThrows(RuntimeException.class, () -> y.get(dm(0, 5)));
    assertThrows(RuntimeException.class, () -> y.set(dm(0, 5), SX.sym("a")));
    y.get(dm(0, 2));
    y.set(dm(0, 2), SX.sym("a"));
  }

  // test_issue181: python-only (Function with nested lists containing None does not compile in Java)

  @Disabled("known bug in Python too: not implemented")
  @Test
  public void test_is_equal() {
    message("equivalent");
    SX x = SX.sym("x");
    SX a = x.times(x);
    SX b = x.times(x);
    assertTrue(ca.is_equal(a, b, 1));
  }

  @Test
  public void test_SXsimplifications() {
    Assumptions.assumeTrue(GlobalOptions.getSimplificationOnTheFly());
    message("simplifications");
    SX x = SX.sym("x");
    SX y0 = SX.sym("x");
    List<UnaryOperator<SX>> ops = new ArrayList<>();
    ops.add(v -> {
      SX y = ca.times(0.5, v);
      return y.plus(y);
    });
    ops.add(v -> {
      SX y = v.div(2);
      return y.plus(y);
    });
    ops.add(v -> {
      SX y = v.times(0.5);
      return y.plus(y);
    });
    ops.add(v -> {
      SX y = v.times(v);
      return y.unaryMinus().div(y).times(v.unaryMinus());
    });
    ops.add(v -> v.times(v).unaryMinus().div(v.times(v)).times(v.unaryMinus()));
    ops.add(v -> {
      SX y = v.times(v);
      return y.div(y.unaryMinus()).times(v.unaryMinus());
    });
    ops.add(v -> {
      SX y = v.times(v);
      return y.unaryMinus().div(y.unaryMinus()).times(v);
    });
    ops.add(v -> v.minus(v).plus(v));
    ops.add(v -> v.times(v).minus(v.times(v)).plus(v));
    ops.add(v -> ca.times(4.0, ca.times(0.25, v)));
    ops.add(v -> ca.times(4.0, v.times(0.25)));
    ops.add(v -> ca.times(4.0, ca.times(0.25, v)));
    ops.add(v -> ca.times(4.0, v.times(0.25)));
    ops.add(v -> ca.times(0.25, v).times(4));
    ops.add(v -> v.times(0.25).times(4));
    ops.add(v -> ca.times(4.0, v).div(4));
    ops.add(v -> ca.times(4.0, v.div(4)));
    ops.add(v -> v.div(4).div(0.25));
    ops.add(v -> v.times(ca.rdivide(4.0, v).times(v).div(4)));
    ops.add(v -> v.times(v.times(ca.rdivide(2.0, v)).div(2)));
    ops.add(v -> v.times(ca.times(2.0, v).div(v).div(2)));
    ops.add(v -> v.times(v.div(ca.times(2.0, v)).times(2)));
    ops.add(v -> v.plus(0));
    ops.add(v -> ca.plus(0.0, v));
    ops.add(v -> v.minus(0));
    ops.add(v -> ca.minus(0.0, v.unaryMinus()));
    ops.add(v -> v.times(1));
    ops.add(v -> ca.times(1.0, v));
    ops.add(v -> ca.times(1.0, v.times(1)));
    ops.add(v -> ca.times(1.0, v).times(1));
    ops.add(v -> ca.times(0.5, v).plus(ca.times(0.5, v)));
    ops.add(v -> v.div(2).plus(v.div(2)));
    ops.add(v -> v.times(0.5).plus(ca.times(0.5, v)));
    ops.add(v -> new SX(4).minus(new SX(4)).plus(v));
    ops.add(v -> v.plus(y0).minus(y0.plus(v)).plus(v));
    ops.add(v -> v.times(y0).minus(y0.times(v)).plus(v));
    ops.add(v -> v.unaryMinus().minus(v.unaryMinus()).plus(v));

    for (UnaryOperator<SX> op : ops) {
      SX y = op.apply(x);
      Function f = new Function("f", List.of(x), List.of(y));
      List<DM> f_out = f.call(dms(new DM(0.3)));
      checkarray(f_out.get(0), ca.evalf(op.apply(new SX(0.3))), "simplifications");
      assertEquals("x", y.toString());

      y = op.apply(x.unaryMinus());
      f = new Function("f", List.of(x), List.of(y));
      f_out = f.call(dms(new DM(0.3)));
      checkarray(f_out.get(0), ca.evalf(op.apply(new SX(-0.3))), "simplifications");
      assertEquals("(-x)", y.toString());
    }
  }

  @Test
  public void test_truth() {
    message("Truth values");
    assertThrows(RuntimeException.class, () -> SX.sym("x").nonzero());
    assertThrows(RuntimeException.class, () -> ca.gt(SX.sym("x"), 0.0).nonzero());
    assertTrue(new SX(1).nonzero());
    assertFalse(new SX(0).nonzero());
    assertTrue(new SX(0.2).nonzero());
    assertTrue(new SX(-0.2).nonzero());
    assertThrows(RuntimeException.class, () -> SX.sym("x").nonzero());
    assertThrows(RuntimeException.class, () -> ca.gt(SX.sym("x"), 0.0).nonzero());
    assertTrue(new SX(new SX(1)).nonzero());
    assertFalse(new SX(new SX(0)).nonzero());
    assertTrue(new SX(new SX(0.2)).nonzero());
    assertTrue(new SX(new SX(-0.2)).nonzero());
    assertThrows(RuntimeException.class, () -> new SX(new double[] {2.0, 3}).nonzero());
  }

  @Test
  public void test_if_else() {
    SX x = SX.sym("x");
    SX y = ca.if_else(x, 1.0, 2.0);
    Function f = new Function("f", List.of(x), List.of(y));
    List<DM> f_out = f.call(dms(new DM(1)));
    assertTrue(f_out.get(0).to_double() == 1, "if_else");
    f_out = f.call(dms(new DM(0)));
    assertTrue(f_out.get(0).to_double() == 2, "if_else");

    double x0 = 2.1;
    y = ca.if_else(ca.gt(x, 1.0), x.power(2), x.power(3));
    f = new Function("f", List.of(x), List.of(y));
    f_out = f.call(dms(new DM(x0)));
    checkarray(f_out.get(0), new DM(Math.pow(x0, 2)), "if_else sens");

    x0 = -2.1;
    f_out = f.call(dms(new DM(x0)));
    checkarray(f_out.get(0), new DM(Math.pow(x0, 3)), "if_else sens");
  }

  @Test
  public void test_is_regular() {
    SX x = SX.sym("x");
    assertTrue(new SX(0).is_regular());
    assertFalse(new SX(Double.POSITIVE_INFINITY).is_regular());
    assertThrows(RuntimeException.class, () -> assertTrue(x.get_nz(false, 0).nonzero()));
    assertTrue(new SX(dm(0, 1)).is_regular());
    assertFalse(new SX(dm(0, Double.POSITIVE_INFINITY)).is_regular());
    assertFalse(vcat(x, new SX(Double.POSITIVE_INFINITY)).is_regular());
    assertThrows(RuntimeException.class, () -> assertFalse(vcat(x, x).is_regular()));
  }

  @Test
  public void test_symvar() {
    SX a = SX.sym("a");
    SX b = SX.sym("b");
    SX c = SX.sym("c");
    SX e = ca.cos(a.times(b)).plus(c);
    List<SX> w = ca.symvar(e);
    assertEquals(3, w.size());
    if (GlobalOptions.getSimplificationOnTheFly()) {
      assertTrue(ca.is_equal(w.get(0), a));
      assertTrue(ca.is_equal(w.get(1), b));
      assertTrue(ca.is_equal(w.get(2), c));
    }
  }

  @Test
  public void test_poly_coeff() {
    SX x = SX.sym("x");
    SX a = SX.sym("a");
    SX c = SX.sym("c");
    SX p = ca.poly_coeff(ca.times(12.0, x.power(4)).plus(x.power(2)).plus(a.times(x)).plus(c), x);
    assertTrue(ca.is_equal(p.get(0), 12.0));
    assertTrue(ca.is_equal(p.get(1), 0.0));
    assertTrue(ca.is_equal(p.get(2), 1.0));
    assertTrue(ca.is_equal(p.get(3), a));
    assertTrue(ca.is_equal(p.get(4), c));

    p = ca.poly_coeff(x.minus(a).times(x.plus(a)), x);
    assertTrue(ca.is_equal(p.get(0), 1.0));
    assertTrue(ca.is_equal(p.get(1), 0.0));
  }

  @Test
  public void test_poly_roots() {
    SX p = SX.sym("[a,b]");
    SX r = ca.poly_roots(p);
    Function f = new Function("f", List.of(p), List.of(r));
    DM f_in = dm(2, 7);
    DM a_ = f_in.get(0), b_ = f_in.get(1);
    List<DM> f_out = f.call(List.of(f_in));
    checkarray(f_out.get(0), vcat(b_.unaryMinus().div(a_)));

    p = SX.sym("[a,b]");
    r = ca.poly_roots(vcat(p, new SX(0)));
    f = new Function("f", List.of(p), List.of(r));
    f_out = f.call(List.of(f_in));
    checkarray(f_out.get(0), vcat(b_.unaryMinus().div(a_), new DM(0)));

    p = SX.sym("[a,b,c]");
    r = ca.poly_roots(p);
    f = new Function("f", List.of(p), List.of(r));
    f_in = dm(1.13, 7, 3);
    a_ = f_in.get(0);
    b_ = f_in.get(1);
    DM c_ = f_in.get(2);
    DM d = b_.power(2).minus(ca.times(4.0, a_).times(c_));
    f_out = f.call(List.of(f_in));
    DM x0 = b_.unaryMinus().minus(ca.sqrt(d)).div(2).div(a_);
    DM x1 = b_.unaryMinus().plus(ca.sqrt(d)).div(2).div(a_);
    checkarray(f_out.get(0), vcat(x0, x1));

    p = SX.sym("[a,b,c,d]");
    r = ca.poly_roots(p);
    f = new Function("f", List.of(p), List.of(r));
    f_out = f.call(List.of(dm(11, 1.3, -1.7, 0.1)));
    checkarray(f_out.get(0), dm(0.298028, -0.479787, 0.0635774), 1e-5, "");

    p = SX.sym("[a,b,c,d,e]");
    r = ca.poly_roots(p);
    f = new Function("f", List.of(p), List.of(r));
    f_out = f.call(List.of(dm(3, 6, -123, -126, 1080)));
    checkarray(f_out.get(0), dm(5, 3, -4, -6), 1e-5, "");
  }

  @Test
  public void test_eig_symbolic() {
    SX x = SX.sym("x", 2, 2);
    Function f = new Function("f", List.of(x), List.of(ca.eig_symbolic(x)));
    List<DM> f_out = f.call(List.of(dm(new double[][] {{2, 0.1}, {0.3, 0.7}})));
    checkarray(f_out.get(0), dm(0.67732, 2.02268), 1e-5, "");

    x = SX.sym("x", 2);
    f = new Function("f", List.of(x), List.of(ca.eig_symbolic(ca.diag(x))));
    DM f_in = dm(3, 7);
    f_out = f.call(List.of(f_in));
    checkarray(f_out.get(0), f_in);

    x = SX.sym("x", 5);
    f = new Function("f", List.of(x), List.of(ca.eig_symbolic(ca.diag(x))));
    f_in = dm(3, 7, 2, 1, 6);
    f_out = f.call(List.of(f_in));
    checkarray(f_out.get(0), f_in);

    x = SX.sym("x", 2, 2);
    SX y = SX.sym("y", 2);
    f = new Function("f", List.of(x, y), List.of(ca.eig_symbolic(ca.diagcat(List.of(x, ca.diag(y))))));
    f_out = f.call(dms(dm(new double[][] {{2, 0.1}, {0.3, 0.7}}), dm(3, 7)));
    checkarray(f_out.get(0), dm(0.67732, 2.02268, 3, 7), 1e-5, "");

    x = SX.sym("x", 3, 3);
    x.set(2, 0, 0.0);
    x.set(1, 0, 0.0);
    x = ca.sparsify(x);
    SX e = ca.eig_symbolic(x);
    f = new Function("f", List.of(x), List.of(e));
    f_in = new DM(f.sparsity_in(0), dm(1, 2, 3, 4, 5, 6, 7));
    f_in.print_dense();
    f_out = f.call(List.of(f_in));
    checkarray(f_out.get(0), dm(1, -0.29150, 10.29150), 1e-5, "");

    x = SX.sym("x", 3, 3);
    x.set(2, 0, 0.0);
    x.set(1, 0, 0.0);
    x.set(2, 1, 0.0);
    x = ca.sparsify(x);
    e = ca.eig_symbolic(x);
    f = new Function("f", List.of(x), List.of(e));
    f_in = new DM(f.sparsity_in(0), dm(1, 2, 3, 4, 5, 6));
    f_in.print_dense();
    f_out = f.call(List.of(f_in));
    checkarray(f_out.get(0), dm(1, 3, 6), 1e-5, "");

    x = SX.sym("x", Sparsity.upper(5));
    f = new Function("f", List.of(x), List.of(ca.eig_symbolic(x)));
    DM fin = new DM(x.sparsity(), 0);
    fin.set(Sparsity.diag(5), ca.diag(range(5)));
    checkarray((DM) f.invoke(fin), range(5));
  }

  @Test
  public void test_jacobian_empty() {
    SX x = SX.sym("x", 3);
    SX s = ca.jacobian(new DM(0, 0), x);
    assertEquals(0, s.size1());
    assertEquals(3, s.size2());

    s = ca.jacobian(x, SX.sym("x", 0, 4));
    assertEquals(3, s.size1());
    assertEquals(0, s.size2());
  }

  @Test
  public void test_empty_SX() {
    SX s = new SX(new double[0]);
    assertEquals(0, s.size1());
    assertEquals(1, s.size2());
    vcat(SX.sym("x"), new SX(new double[0]));
    new SX(new Sparsity(3, 3), new DM(new double[0]));
    new SX(new Sparsity(3, 3), new DM(new double[] {5}));
    assertInException("fully sparse", () -> new SX(new Sparsity(3, 3), new DM(new double[] {5, 4})));
  }

  @Test
  public void test_mul_sparsity() {
    int N = 10;
    SX x = SX.sym("x", N, N);
    SX y = SX.sym("y", N, N);
    DM x_ = randDM(N, N);
    DM y_ = randDM(N, N);
    Sparsity filt = Sparsity.diag(N).plus(Sparsity.triplet(N, N, new long[] {1}, new long[] {3}));
    Function f = new Function("f", List.of(x, y), List.of(ca.mtimes(x, y)));
    Function g = new Function("g", List.of(x, y), List.of(ca.mac(x, y, SX.zeros(filt))));
    List<DM> f_out = f.call(dms(x_, y_));
    List<DM> g_out = g.call(dms(x_, y_));
    checkarray(DM.ones(filt), DM.ones(g.sparsity_out(0)));
    checkarray(f_out.get(0).get(filt), g_out.get(0));
  }

  @Test
  public void test_large_hessian() {
    Sparsity A = Sparsity.from_file("../data/apoa1-2.mtx");
    DM H = new DM(A, range((int) A.nnz()));
    H = H.plus(H.T());
    H = H.get(new Slice(0, 20000), new Slice(0, 20000));
    SX x = SX.sym("x", H.size1());
    Function f = new Function("f", List.of(x), List.of(ca.mtimes(List.of(x.T(), new SX(H), x))), Map.of("verbose", true));
    H = H.times(2);
    Function h = hessianOld(f, 0, 0);
    List<DM> h_out = h.call(dms(new DM(0)));
    assertTrue(h.sparsity_out(0).is_equal(H.sparsity()));
    checkarray(new DM(h_out.get(0).nonzeros()), new DM(H.nonzeros()));
  }

  @Test
  public void test_mxnulloutput() {
    SX a = new SX(5, 0);
    SX b = SX.sym("x", 2);
    MX bm = MX.sym("x", 2);

    Function f = new Function("f", List.of(b), List.of(a));
    MX c = (MX) f.invoke(bm);
    assertEquals(5, c.size1());
    assertEquals(0, c.size2());
    SX cs = (SX) f.invoke(b);
    assertEquals(5, cs.size1());
    assertEquals(0, cs.size2());

    a = new SX(0, 0);
    f = new Function("f", List.of(b), List.of(a));
    c = (MX) f.invoke(bm);
    assertEquals(0, c.size1());
    assertEquals(0, c.size2());
    cs = (SX) f.invoke(b);
    assertEquals(0, cs.size1());
    assertEquals(0, cs.size2());
  }

  @Test
  public void test_mxnull() {
    SX a = new SX(5, 0);
    SX b = new SX(0, 3);
    SX c = ca.mtimes(a, b);
    assertEquals(0, c.nnz());

    a = new SX(5, 3);
    b = new SX(3, 4);
    c = ca.mtimes(a, b);
    assertEquals(0, c.nnz());
  }

  @Test
  public void test_mxnullop() {
    SX c = new SX(0, 0);
    SX x = SX.sym("x", 2, 3);
    assertThrows(RuntimeException.class, () -> x.plus(c));
    assertThrows(RuntimeException.class, () -> x.div(c));
  }

  @Test
  @SuppressWarnings("unchecked")
  public void test_copysign_derivatives() {
    for (int t = 0; t < 2; t++) {
      Function f;
      if (t == 0) {
        SX x = SX.sym("x"), y = SX.sym("y"), z = ca.copysign(x, y);
        f = new Function("f", List.of(x, y), List.of(ca.jacobian(z, x), ca.gradient(z, x), ca.jtimes(z, x, 1.0),
            ca.jtimes(z, x, 1.0, true), ca.jacobian(z, y)));
      } else {
        MX x = MX.sym("x"), y = MX.sym("y"), z = ca.copysign(x, y);
        f = new Function("f", List.of(x, y), List.of(ca.jacobian(z, x), ca.gradient(z, x), ca.jtimes(z, x, 1.0),
            ca.jtimes(z, x, 1.0, true), ca.jacobian(z, y)));
      }
      for (double x0 : new double[] {-2, 0, 2})
        for (double y0 : new double[] {-0.5, -0.0, 0.0, 0.5}) {
          double d = Math.signum(x0) * Math.copySign(1, y0);
          List<DM> r = (List<DM>) (List<?>) f.invoke(x0, y0);
          checkarray(ca.vertcat(r), dm(d, d, d, d, 0), (t == 0 ? "SX" : "MX") + " x=" + x0 + " y=" + y0);
        }
    }
  }

  @Test
  public void test_copysign() {
    SX x = SX.sym("x");
    SX y = SX.sym("y");
    SX z = ca.copysign(x, y);
    Function f = new Function("f", List.of(x, y), List.of(z));
    double[][] in = {{2, 0.5}, {2, -0.5}, {-2, 0.5}, {-2, -0.5}, {2, 0}};
    double[] ref = {2, -2, 2, -2, 2};
    for (int i = 0; i < in.length; i++)
      checkarray(f.call(dms(new DM(in[i][0]), new DM(in[i][1]))).get(0), dm(ref[i]));

    Function J = jacobianOld(f, 0, 0);
    double[] jref = {1, -1, -1, 1, 1};
    for (int i = 0; i < in.length; i++)
      checkarray(J.call(dms(new DM(in[i][0]), new DM(in[i][1]))).get(0), dm(jref[i]));

    J = jacobianOld(f, 1, 0);
    for (int i = 0; i < in.length; i++)
      checkarray(J.call(dms(new DM(in[i][0]), new DM(in[i][1]))).get(0), dm(0));
  }

  @Test
  public void test_depends_on() {
    SX a = SX.sym("a");
    SX b = SX.sym("b");
    SX ab = vcat(a, b);
    SX zero = new SX(0);
    assertTrue(ca.depends_on(a.power(2), a));
    assertTrue(ca.depends_on(a, a));
    assertFalse(ca.depends_on(0.0, a));
    assertTrue(ca.depends_on(a.power(2), ab));
    assertTrue(ca.depends_on(a, ab));
    assertFalse(ca.depends_on(0.0, ab));
    assertTrue(ca.depends_on(b.power(2), ab));
    assertTrue(ca.depends_on(b, ab));
    assertTrue(ca.depends_on(a.power(2).plus(b.power(2)), ab));
    assertTrue(ca.depends_on(a.plus(b), ab));
    assertTrue(ca.depends_on(vcat(zero, a), a));
    assertTrue(ca.depends_on(vcat(a, zero), a));
    assertTrue(ca.depends_on(vcat(a.power(2), b.power(2)), ab));
    assertTrue(ca.depends_on(vcat(a, zero), ab));
    assertTrue(ca.depends_on(vcat(zero, b), ab));
    assertTrue(ca.depends_on(vcat(b, zero), ab));
    assertFalse(ca.depends_on(vcat(zero, zero), ab));
  }

  // test_is_smooth: python-only (DeprecationWarning of a free function is_smooth)

  @Test
  public void test_which_depends() {
    {
      SX x = SX.sym("x"), y = SX.sym("y"), p = SX.sym("p"), z = SX.sym("z");
      SX e = vcat(new SX(0), x, y, p, ca.times(2.0, p.power(3)), x.times(y), x.times(p), ca.sin(x), ca.cos(y),
          ca.sqrt(x.plus(y)), p.times(p).times(x), x.times(y).times(p));
      assertBools(ca.which_depends(e, vcat(x, y), 2, true), 0, 0, 0, 0, 0, 1, 0, 1, 1, 1, 0, 1);
      assertBools(ca.which_depends(e, vcat(x, y), 1, true), 0, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1);
      e = vcat(x.times(p), x.plus(y));
      assertBools(ca.which_depends(e, vcat(x, y, p, z), 2, false), 1, 0, 1, 0);
      assertBools(ca.which_depends(e, vcat(x, y, p, z), 1, false), 1, 1, 1, 0);
      e = vcat(x.times(p), x.plus(z.times(y)));
      assertBools(ca.which_depends(e, vcat(x, y, p), 2, false), 1, 0, 1);
      assertBools(ca.which_depends(e, vcat(x, y, p), 1, false), 1, 1, 1);
      assertBools(ca.which_depends(e, vcat(x, y, p, z), 2, false), 1, 1, 1, 1);
      assertBools(ca.which_depends(e, vcat(x, y, p, z), 1, false), 1, 1, 1, 1);
      e = vcat(ca.sin(x.plus(y)).plus(p));
      assertBools(ca.which_depends(e, vcat(x, y, p, z), 2, false), 1, 1, 0, 0);
      assertBools(ca.which_depends(e, vcat(x, y, p, z), 1, false), 1, 1, 1, 0);
      e = vcat(ca.sin(x).times(p.power(2)), y.power(2));
      assertBools(ca.which_depends(e, vcat(x, y, p), 2, true), 1, 1);
      assertBools(ca.which_depends(e, vcat(x, y, p), 2, false), 1, 1, 1);
      e = vcat(x.power(2).times(p), y);
      assertBools(ca.which_depends(e, vcat(x, y, p), 2, true), 1, 0);
      assertBools(ca.which_depends(e, vcat(x, y, p), 2, false), 1, 0, 1);
    }
    {
      MX x = MX.sym("x"), y = MX.sym("y"), p = MX.sym("p"), z = MX.sym("z");
      MX e = vcat(new MX(0), x, y, p, ca.times(2.0, p.power(3)), x.times(y), x.times(p), ca.sin(x), ca.cos(y),
          ca.sqrt(x.plus(y)), p.times(p).times(x), x.times(y).times(p));
      assertBools(ca.which_depends(e, vcat(x, y), 2, true), 0, 0, 0, 0, 0, 1, 0, 1, 1, 1, 0, 1);
      assertBools(ca.which_depends(e, vcat(x, y), 1, true), 0, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1);
      e = vcat(x.times(p), x.plus(y));
      assertBools(ca.which_depends(e, vcat(x, y, p, z), 2, false), 1, 0, 1, 0);
      assertBools(ca.which_depends(e, vcat(x, y, p, z), 1, false), 1, 1, 1, 0);
      e = vcat(x.times(p), x.plus(z.times(y)));
      assertBools(ca.which_depends(e, vcat(x, y, p), 2, false), 1, 0, 1);
      assertBools(ca.which_depends(e, vcat(x, y, p), 1, false), 1, 1, 1);
      assertBools(ca.which_depends(e, vcat(x, y, p, z), 2, false), 1, 1, 1, 1);
      assertBools(ca.which_depends(e, vcat(x, y, p, z), 1, false), 1, 1, 1, 1);
      e = vcat(ca.sin(x.plus(y)).plus(p));
      assertBools(ca.which_depends(e, vcat(x, y, p, z), 2, false), 1, 1, 0, 0);
      assertBools(ca.which_depends(e, vcat(x, y, p, z), 1, false), 1, 1, 1, 0);
      e = vcat(ca.sin(x).times(p.power(2)), y.power(2));
      assertBools(ca.which_depends(e, vcat(x, y, p), 2, true), 1, 1);
      assertBools(ca.which_depends(e, vcat(x, y, p), 2, false), 1, 1, 1);
      e = vcat(x.power(2).times(p), y);
      assertBools(ca.which_depends(e, vcat(x, y, p), 2, true), 1, 0);
      assertBools(ca.which_depends(e, vcat(x, y, p), 2, false), 1, 0, 1);
    }
  }

  @Test
  @SuppressWarnings("unchecked")
  public void test_if_else_zero_sens() {
    SX x = SX.sym("x");
    SX a = ca.plus(1.0, ca.times(3.0, x)).plus(ca.sqrt(ca.times(3.0, x)).times(x)).plus(ca.times(7.0, x));
    SX b = ca.plus(1.0, ca.times(2.0, x)).plus(ca.sin(ca.times(2.0, x)).times(x)).plus(x);
    SX z = ca.if_else(ca.gt(x, 0.0), a, b).times(x);

    Function f = new Function("f", List.of(x), List.of(z, ca.jacobian(z, x)));
    Function fa = new Function("f", List.of(x), List.of(a.times(x), ca.jacobian(a.times(x), x)));
    Function fb = new Function("f", List.of(x), List.of(b.times(x), ca.jacobian(b.times(x), x)));

    List<DM> r = (List<DM>) (List<?>) f.invoke(3), ra = (List<DM>) (List<?>) fa.invoke(3);
    for (int i = 0; i < r.size(); i++) checkarray(r.get(i), ra.get(i));
    r = (List<DM>) (List<?>) f.invoke(-3);
    List<DM> rb = (List<DM>) (List<?>) fb.invoke(-3);
    for (int i = 0; i < r.size(); i++) checkarray(r.get(i), rb.get(i));

    f = new Function("f", List.of(x), List.of(z));
    fa = new Function("f", List.of(x), List.of(a.times(x)));
    fb = new Function("f", List.of(x), List.of(b.times(x)));
    checkfunction(f, fa, dms(new DM(3)));
    checkfunction(f, fb, dms(new DM(-3)));
  }

  @Test
  public void test_pw_const() {
    SX t = SX.sym("t");
    SX e = ca.pw_const(t, dm(0, 2, 3), dm(7, 1, 3, 5));
    Function E = new Function("E", List.of(t), List.of(e));
    double[][] cases = {{-2, 7}, {-1, 7}, {0, 1}, {1, 1}, {1.9999, 1}, {2, 3}, {2.5, 3}, {3, 5}, {10, 5}};
    for (double[] c : cases) checkarray((DM) E.invoke(c[0]), c[1]);
  }

  @Test
  public void test_pw_lin() {
    SX t = SX.sym("t");
    SX e = ca.pw_lin(t, dm(0, 2, 3, 5), dm(7, 1, 3, 2));
    Function E = new Function("E", List.of(t), List.of(e));
    double[][] cases = {{-2, 13}, {-1, 10}, {0, 7}, {1, 4}, {2, 1}, {2.5, 2}, {3, 3}, {4, 2.5}, {5, 2}, {7, 1}};
    for (double[] c : cases) checkarray((DM) E.invoke(c[0]), c[1]);
  }

  // test_numpy_error: python-only (numpy dispatch on symbolic values)

  @Test
  @SuppressWarnings("unchecked")
  public void test_quadratic() {
    {
      SX x = SX.sym("x"), p = SX.sym("p"), y = SX.sym("y");
      assertFalse(ca.is_quadratic(ca.sin(x), x));
      assertFalse(ca.is_quadratic(x.power(3), x));
      assertTrue(ca.is_quadratic(x.power(2), x));
      assertTrue(ca.is_quadratic(ca.times(4.0, x), x));
      assertTrue(ca.is_quadratic(5.0, x));

      assertFalse(ca.is_quadratic(ca.sin(x).times(p.power(4)), x));
      assertFalse(ca.is_quadratic(x.power(3).times(p.power(4)), x));
      assertTrue(ca.is_quadratic(x.power(2).times(p.power(4)), x));
      assertTrue(ca.is_quadratic(x.times(p.power(4)), x));
      assertTrue(ca.is_quadratic(ca.times(5.0, p.power(4)), x));

      assertFalse(ca.is_linear(ca.sin(x), x));
      assertFalse(ca.is_linear(x.power(3), x));
      assertFalse(ca.is_linear(x.power(2), x));
      assertTrue(ca.is_linear(ca.times(3.0, x), x));
      assertTrue(ca.is_linear(5.0, x));

      assertFalse(ca.is_linear(ca.sin(x).times(p.power(4)), x));
      assertFalse(ca.is_linear(x.power(3).times(p.power(4)), x));
      assertFalse(ca.is_linear(x.power(2).times(p.power(4)), x));
      assertTrue(ca.is_linear(x.times(p.power(4)), x));
      assertTrue(ca.is_linear(ca.times(5.0, p.power(4)), x));

      SX xy = vcat(x, y);
      SX base = x.power(2).plus(ca.times(3.0, y.power(2))).plus(ca.times(0.5, x).times(y));
      SX z = base.plus(ca.times(7.0, x)).plus(ca.times(6.0, y)).plus(7);
      List<Object> ABc = ca.quadratic_coeff(z, xy);
      assertInException("non-quadratic", () -> ca.quadratic_coeff(z.plus(ca.sin(x)), xy));
      assertInException("scalar", () -> ca.quadratic_coeff(xy, x));

      SX z2 = base.minus(p.times(y)).plus(ca.times(7.0, x)).plus(ca.times(6.0, y)).plus(7);
      ABc = ca.quadratic_coeff(z2, xy);
      SX A = (SX) ABc.get(0), b = (SX) ABc.get(1), c = (SX) ABc.get(2);
      SX e = ca.times(0.5, ca.bilin(A, xy, xy)).plus(ca.dot(b, xy)).plus(c);
      Function f = new Function("f", List.of(xy, p), List.of(z2));
      Function f2 = new Function("f", List.of(xy, p), List.of(e));
      checkfunction(f, f2, dms(new DM(1.1), new DM(1.3)));

      assertInException("non-linear", () -> ca.linear_coeff(z, xy));
      assertInException("vector", () -> ca.linear_coeff(ca.blockcat(List.of(List.of(x, y), List.of(y, x))), x));

      SX z3 = vcat(ca.times(7.0, x).plus(ca.times(6.0, y)).plus(7), ca.minus(5.0, p.times(y)));
      List<Object> Ab = ca.linear_coeff(z3, xy);
      e = ca.mtimes((SX) Ab.get(0), xy).plus((SX) Ab.get(1));
      f = new Function("f", List.of(xy, p), List.of(z3));
      f2 = new Function("f", List.of(xy, p), List.of(e));
      checkfunction(f, f2, dms(new DM(1.1), new DM(1.3)));
    }
    {
      MX x = MX.sym("x"), p = MX.sym("p"), y = MX.sym("y");
      assertFalse(ca.is_quadratic(ca.sin(x), x));
      assertFalse(ca.is_quadratic(x.power(3), x));
      assertTrue(ca.is_quadratic(x.power(2), x));
      assertTrue(ca.is_quadratic(ca.times(4.0, x), x));
      assertTrue(ca.is_quadratic(5.0, x));

      assertFalse(ca.is_quadratic(ca.sin(x).times(p.power(4)), x));
      assertFalse(ca.is_quadratic(x.power(3).times(p.power(4)), x));
      assertTrue(ca.is_quadratic(x.power(2).times(p.power(4)), x));
      assertTrue(ca.is_quadratic(x.times(p.power(4)), x));
      assertTrue(ca.is_quadratic(ca.times(5.0, p.power(4)), x));

      assertFalse(ca.is_linear(ca.sin(x), x));
      assertFalse(ca.is_linear(x.power(3), x));
      assertFalse(ca.is_linear(x.power(2), x));
      assertTrue(ca.is_linear(ca.times(3.0, x), x));
      assertTrue(ca.is_linear(5.0, x));

      assertFalse(ca.is_linear(ca.sin(x).times(p.power(4)), x));
      assertFalse(ca.is_linear(x.power(3).times(p.power(4)), x));
      assertFalse(ca.is_linear(x.power(2).times(p.power(4)), x));
      assertTrue(ca.is_linear(x.times(p.power(4)), x));
      assertTrue(ca.is_linear(ca.times(5.0, p.power(4)), x));

      MX xy = vcat(x, y);
      MX base = x.power(2).plus(ca.times(3.0, y.power(2))).plus(ca.times(0.5, x).times(y));
      MX z = base.plus(ca.times(7.0, x)).plus(ca.times(6.0, y)).plus(7);
      List<Object> ABc = ca.quadratic_coeff(z, xy);
      assertInException("non-quadratic", () -> ca.quadratic_coeff(z.plus(ca.sin(x)), xy));
      assertInException("scalar", () -> ca.quadratic_coeff(xy, x));

      MX z2 = base.minus(p.times(y)).plus(ca.times(7.0, x)).plus(ca.times(6.0, y)).plus(7);
      ABc = ca.quadratic_coeff(z2, xy);
      MX A = (MX) ABc.get(0), b = (MX) ABc.get(1), c = (MX) ABc.get(2);
      MX e = ca.times(0.5, ca.bilin(A, xy, xy)).plus(ca.dot(b, xy)).plus(c);
      Function f = new Function("f", List.of(xy, p), List.of(z2));
      Function f2 = new Function("f", List.of(xy, p), List.of(e));
      checkfunction(f, f2, dms(new DM(1.1), new DM(1.3)));

      assertInException("non-linear", () -> ca.linear_coeff(z, xy));
      assertInException("vector", () -> ca.linear_coeff(ca.blockcat(List.of(List.of(x, y), List.of(y, x))), x));

      MX z3 = vcat(ca.times(7.0, x).plus(ca.times(6.0, y)).plus(7), ca.minus(5.0, p.times(y)));
      List<Object> Ab = ca.linear_coeff(z3, xy);
      e = ca.mtimes((MX) Ab.get(0), xy).plus((MX) Ab.get(1));
      f = new Function("f", List.of(xy, p), List.of(z3));
      f2 = new Function("f", List.of(xy, p), List.of(e));
      checkfunction(f, f2, dms(new DM(1.1), new DM(1.3)));
    }
  }

  @Test
  public void test_evalf() {
    SX x = SX.sym("x");
    SX y = new SX(5);
    checkarray(ca.evalf(y), 5);
    assertInException("since variables [x] are free", () -> ca.evalf(x));
  }

  @Test
  public void test_output_sx() {
    MX x = MX.sym("x");
    MX y = MX.sym("y");
    MX z = MX.sym("z");
    MX w = ca.mtimes(ca.sqrt(z), ca.sin(x.times(y)));
    Function f = new Function("f", List.of(x, y, z), List.of(x.times(y).times(z), w));
    System.out.println(f);

    SX xs = SX.sym("x");
    SX ys = SX.sym("y");
    SX zs = SX.sym("z");
    List<SX> args = List.of(xs.div(ys), ca.times(10.0, zs), ys.times(zs));
    List<SX> v = f.call(args, false, true);
    SX output_0 = v.get(0);
    SX output_1 = v.get(1);
    assertEquals(output_0.dep(0).element_hash(), output_1.dep(0).element_hash());
    SX call_node = output_0.dep(0);
    System.out.println("here");
    SX output_1b = call_node.get_output(1);
    assertEquals(output_1.element_hash(), output_1b.element_hash());
    output_1b = call_node.get_output(1);
    assertNotNull(output_1b);
  }

  @Test
  @SuppressWarnings("unchecked")
  public void test_call_fun() {
    DM A = ca.sparsify(dm(new double[][] {{1, 0, 1}, {0, 0, 6}, {0, 8, 9}}));
    MX xm = MX.sym("x", 3);
    MX ym = MX.sym("y");
    MX zm = MX.sym("z");
    MX wm = ca.mtimes(ca.times(A, ca.sqrt(zm)), ca.sin(xm.times(ym)));
    Function f = new Function("f", List.of(xm, ym, zm), List.of(xm.times(ym).times(zm).get(new Slice(0, 2)), wm));
    System.out.println(f);

    SX x = SX.sym("x", 3);
    SX y = SX.sym("y");
    SX z = SX.sym("z");
    List<SX> args = List.of(x.div(y), ca.times(10.0, z), y.times(z));
    List<SX> v = f.call(args);
    List<SX> v2 = f.call(args, false, true);
    System.out.println(v2);
    System.out.println(ca.vertcat(v2));

    assertTrue(v2.get(0).get(0).is_output());
    assertFalse(y.is_output());
    long k = 0;
    for (int i = 0; i < f.n_out(); i++)
      for (int j = 0; j < f.nnz_out(i); j++) {
        assertEquals(k, v2.get(i).get(j).which_output());
        k += 1;
      }
    SX callnode = v2.get(0).get(0).dep(0);
    assertTrue(callnode.is_call());
    assertFalse(y.is_call());
    assertTrue(callnode.has_output());
    assertFalse(y.has_output());
    assertTrue(callnode.which_function().hash() == f.hash());
    assertTrue(callnode.n_dep() == 5);
    callnode.dep(4);
    SX e = callnode.get_output(3);
    System.out.println(e + " " + e.element_hash());
    SX ee = callnode.get_output(3);
    System.out.println(ee + " " + e.element_hash());
    assertTrue(callnode.get_output(3).toString().contains("{3}"));

    System.out.println(ca.cos(ca.mtimes(v2.get(1), v2.get(0).T()).div(y)).dim());
    Function F2 = new Function("F", List.of(x, y, z), List.of(ca.sin(v2.get(0).times(y).minus(args.get(2))),
        ca.cos(ca.mtimes(v2.get(1), v2.get(0).T()).div(y))));
    Function F2f = F2;
    assertOutput("f:(i0[3],i1,i2)->(o0[2],o1[3])", () -> System.out.println(f));
    assertOutput("[[@2,@3],[@4,@5,@7]] = f([@2,@3,@4],@5,@6);", () -> F2f.disp(true));

    List<DM> inputs = dms(dm(0.1, 1.7, 2.3), new DM(1.13), new DM(0.11));
    Function F1 = new Function("F", List.of(x, y, z), List.of(ca.sin(v.get(0).times(y).minus(args.get(2))),
        ca.cos(ca.mtimes(v.get(1), v.get(0).T()).div(y))));
    for (boolean ad_weight : new boolean[] {true, false})
      for (boolean ad_weight_sp : new boolean[] {true, false}) {
        F2 = new Function("F", List.of(x, y, z), List.of(ca.sin(v2.get(0).times(y).minus(args.get(2))),
            ca.cos(ca.mtimes(v2.get(1), v2.get(0).T()).div(y))), Map.of("ad_weight_sp", ad_weight_sp, "ad_weight",
            ad_weight));
        checkfunction(F1, F2, inputs);
        checkSerialize(F2, inputs);
        for (boolean avoid_stack : new boolean[] {true, false})
          checkCodegen(F2, inputs, Map.of("avoid_stack", avoid_stack));
      }

    SX instr = F2.instructions_sx();
    boolean check = false;
    for (long kk = 0; kk < F2.n_instructions(); kk++) {
      long op = F2.instruction_id(kk);
      if (op == Operation.OP_CALL.swigValue()) {
        check = true;
        assertTrue(instr.get(kk).is_call());
        assertArrayEquals(new long[] {2, 3, 4, 5, 7}, F2.instruction_output(kk));
        assertArrayEquals(new long[] {2, 3, 4, 5, 6}, F2.instruction_input(kk));
      }
    }
    assertTrue(check);

    F2 = new Function("F", List.of(x, y, z), List.of(ca.sin(v2.get(0).times(y).minus(args.get(2))),
        ca.cos(ca.mtimes(v2.get(1), v2.get(0).T()).div(y))), Map.of("cse", true));
    checkfunction(F1, F2, inputs);
    List<Function> res = F2.find_functions();
    assertTrue(res.size() == 1);
    assertTrue(res.get(0).hash() == f.hash());

    F2 = new Function("F", List.of(x, y, z), (List<SX>) (List<?>) F2.invoke(x, y, z));
    checkfunction(F1, F2, inputs);

    F1 = new Function("F", List.of(x, y, z), (List<SX>) (List<?>) F1.invoke(ca.times(2.0, x), ca.times(3.0, y),
        ca.times(4.0, z)));
    F2 = new Function("F", List.of(x, y, z), (List<SX>) (List<?>) F2.invoke(ca.times(2.0, x), ca.times(3.0, y),
        ca.times(4.0, z)));
    checkfunction(F1, F2, inputs);

    F1 = new Function("F", List.of(x, y, z), (List<SX>) (List<?>) F1.invoke(x, ca.times(3.0, y), ca.times(4.0, z)));
    F2 = new Function("F", List.of(x, y, z), (List<SX>) (List<?>) F2.invoke(x, ca.times(3.0, y), ca.times(4.0, z)));
    checkfunction(F1, F2, inputs);

    Function F3 = new Function("F", List.of(x, y, z), List.of(ca.sin(v2.get(0).times(y).minus(args.get(2))),
        ca.cos(ca.mtimes(v2.get(1), v2.get(0).T()).get(1).div(y))));
    assertOutput("[[@2,@3],[NULL,@4,NULL]] = f([@2,@3,@4],@5,@6);", () -> F3.disp(true));

    F1 = new Function("F", List.of(x, y, z), List.of(ca.sin(v.get(0).times(y).minus(args.get(2))),
        ca.cos(ca.mtimes(v.get(1), v.get(0).T()).get(1).div(y))));
    for (boolean ad_weight : new boolean[] {true, false})
      for (boolean ad_weight_sp : new boolean[] {true, false}) {
        F2 = new Function("F", List.of(x, y, z), List.of(ca.sin(v2.get(0).times(y).minus(args.get(2))),
            ca.cos(ca.mtimes(v2.get(1), v2.get(0).T()).get(1).div(y))), Map.of("ad_weight_sp", ad_weight_sp,
            "ad_weight", ad_weight));
        checkfunction(F1, F2, inputs);
        checkSerialize(F2, inputs);
        for (boolean avoid_stack : new boolean[] {true, false})
          checkCodegen(F2, inputs, Map.of("avoid_stack", avoid_stack));
      }
  }

  @Test
  public void test_call_fun2() {
    System.out.println("test_call_fun2");
    for (int t = 0; t < 2; t++)
      for (boolean never_inline : new boolean[] {true, false}) {
        System.out.println((t == 0 ? "SX" : "MX") + " " + never_inline);
        DM A = ca.sparsify(dm(new double[][] {{1, 0, 1}, {0, 0, 6}, {0, 8, 9}}));
        Function f;
        if (t == 0) {
          SX x = SX.sym("x", 3), y = SX.sym("y"), z = SX.sym("z");
          SX w = ca.mtimes(ca.times(A, ca.sqrt(z)), ca.sin(x.times(y)));
          f = new Function("fun", List.of(x, y, z), List.of(x.times(y).times(z).get(new Slice(0, 2)), w),
              Map.of("never_inline", never_inline));
        } else {
          MX x = MX.sym("x", 3), y = MX.sym("y"), z = MX.sym("z");
          MX w = ca.mtimes(ca.times(A, ca.sqrt(z)), ca.sin(x.times(y)));
          f = new Function("fun", List.of(x, y, z), List.of(x.times(y).times(z).get(new Slice(0, 2)), w),
              Map.of("never_inline", never_inline));
        }
        SX x = SX.sym("x", 3);
        SX y = SX.sym("y");
        SX z = SX.sym("z");
        List<SX> args = List.of(x.div(y), ca.times(10.0, z), y.times(z));
        Object res = f.invoke(args.get(0), args.get(1), args.get(2));
        System.out.println(res);
        assertEquals(never_inline, res.toString().contains("fun"));
      }
  }

  @Test
  public void test_call_fun_nominal_out_deriv() {
    Assumptions.assumeTrue(hasRootfinder("newton"));
    SX x = SX.sym("x", 2);
    SX p = SX.sym("p", 2);
    Function rf = ca.rootfinder("rf", "newton", Map.of("x", x, "g", vcat(ca.sin(x.get(0)).minus(p.get(0)),
        ca.sin(x.get(0).plus(x.get(1))).minus(p.get(1).times(p.get(0)))), "p", p));
    Function f = new Function("f", List.of(p), List.of(ca.exp(rf.call(Map.of("p", p.power(2))).get("x"))));
    MX P = MX.sym("p", 2);
    Function fref = new Function("f", List.of(P), List.of(ca.exp(rf.call(Map.of("p", P.power(2))).get("x"))));
    checkfunction(f, fref, dms(dm(0.1, 0.2)));
  }

  @Test
  public void test_call_fun_copy_elision() {
    long old = GlobalOptions.getCopyElisionMinSize();
    GlobalOptions.setCopyElisionMinSize(0);
    try {
      MX mX = MX.sym("A", 3, 6);
      MX mY = MX.sym("Y");
      SX X = SX.sym("A", 3, 6);
      SX Y = SX.sym("Y");
      DM.rng(1);
      Function g = new Function("g", List.of(X, Y), List.of(X, Y));
      Function gref = new Function("g", List.of(mX, mY), List.of(mX, mY));
      List<DM> inputs = dms(DM.rand(3, 3), DM.rand(1));
      copyElisionChecks(g, gref, inputs);

      MX Xm = MX.sym("A", 3, 3);
      MX Ym = MX.sym("Y");
      Function f = new Function("f", List.of(Xm, Ym), List.of(ca.sumsqr(Xm).times(Ym)), Map.of("never_inline", true));
      Function fref = new Function("f", List.of(Xm, Ym), List.of(ca.sumsqr(Xm).times(Ym)));
      X = SX.sym("A", 3, 6);
      Y = SX.sym("Y");
      DM.rng(1);
      Slice all = new Slice();
      g = new Function("g", List.of(X, Y), List.of((SX) f.invoke(X.get(all, new Slice(3, 6)), ca.sin(Y)), X.get(1, 3),
          ca.times(3.0, X.get(2, 4))));
      gref = new Function("g", List.of(X, Y), List.of((SX) fref.invoke(X.get(all, new Slice(3, 6)), ca.sin(Y)),
          X.get(1, 3), ca.times(3.0, X.get(2, 4))));
      inputs = dms(DM.rand(3, 3), DM.rand(1));
      copyElisionChecks(g, gref, inputs);

      Xm = MX.sym("A", 4);
      Ym = MX.sym("Y");
      f = new Function("f", List.of(Xm, Ym), List.of(ca.sumsqr(Xm).times(Ym)), Map.of("never_inline", true));
      fref = new Function("f", List.of(Xm, Ym), List.of(ca.sumsqr(Xm).times(Ym)));
      X = SX.sym("A", 6);
      Y = SX.sym("Y");
      DM.rng(1);
      SX head = vcat(new SX(2), X.get(new Slice(0, 3)));
      g = new Function("g", List.of(X, Y), List.of((SX) f.invoke(head, ca.sin(Y))));
      gref = new Function("g", List.of(X, Y), List.of((SX) fref.invoke(head, ca.sin(Y))));
      inputs = dms(DM.rand(6), DM.rand(1));
      copyElisionChecks(g, gref, inputs);

      g = new Function("g", List.of(X, Y), List.of((SX) f.invoke(head, ca.sin(Y)), ca.times(3.0, X)));
      gref = new Function("g", List.of(X, Y), List.of((SX) fref.invoke(head, ca.sin(Y)), ca.times(3.0, X)));
      inputs = dms(DM.rand(6), DM.rand(1));
      copyElisionChecks(g, gref, inputs);
    } finally {
      GlobalOptions.setCopyElisionMinSize(old);
    }
  }

  private void copyElisionChecks(Function g, Function gref, List<DM> inputs) {
    checkfunction(g, gref, inputs);
    for (boolean avoid_stack : new boolean[] {true, false}) checkCodegen(g, inputs, Map.of("avoid_stack", avoid_stack));
    Function g_roundtrip = Function.deserialize(g.serialize());
    checkfunction(g, g_roundtrip, inputs);
    for (boolean avoid_stack : new boolean[] {true, false})
      checkCodegen(g_roundtrip, inputs, Map.of("avoid_stack", avoid_stack));
  }

  // test_ufunc: python-only (numpy ufunc on a symbolic value)

  @Test
  public void test_mmin() {
    SX x = SX.sym("X", 2);
    Function f0 = new Function("f", List.of(x), List.of(x.get(0).plus(x.get(1)).div(2)));
    Function f1 = new Function("f", List.of(x), List.of(ca.mmin(x)));
    Function f2 = new Function("f", List.of(x), List.of(ca.fmin(x.get(0), x.get(1))));
    checkfunction(f1, f2, dms(dm(0.2, 0.3)));
    checkfunction(f1, f2, dms(dm(2, 2)));
    checkfunction(f1, f0, dms(dm(2, 2)));
    f1 = new Function("f", List.of(x), List.of(ca.mmax(x)));
    f2 = new Function("f", List.of(x), List.of(ca.fmax(x.get(0), x.get(1))));
    checkfunction(f1, f2, dms(dm(0.2, 0.3)));
    checkfunction(f1, f2, dms(dm(2, 2)));
    checkfunction(f1, f0, dms(dm(2, 2)));

    SX xs = SX.sym("X");
    f0 = new Function("f", List.of(xs), List.of(xs));
    f1 = new Function("f", List.of(xs), List.of(ca.fmin(xs, xs)));
    checkfunction(f1, f0, dms(dm(1)));
    f1 = new Function("f", List.of(xs), List.of(ca.fmax(xs, xs)));
    checkfunction(f1, f0, dms(dm(1)));
  }

  @Test
  public void test_empty_broadcast() {
    for (int nc : new int[] {0, 2}) {
      SX res = ca.atan2(SX.sym("c", nc, 1), SX.sym("t", nc, 3));
      assertEquals(nc, res.size1());
      assertEquals(3, res.size2());

      assertInException("Dimension mismatch", () -> ca.atan2(MX.sym("c", nc, 2), MX.sym("t", nc, 3)));

      MX resm = ca.atan2(MX.sym("c", nc, 2), MX.sym("t", nc, 4));
      assertEquals(nc, resm.size1());
      assertEquals(4, resm.size2());

      resm = ca.atan2(MX.sym("c", nc, 4), MX.sym("t", nc, 2));
      assertEquals(nc, resm.size1());
      assertEquals(4, resm.size2());
    }
  }

  @Test
  public void test_logsumexp() {
    SX x = SX.sym("x", 3);
    Function f_ref = new Function("f_ref", List.of(x),
        List.of(ca.log(ca.exp(x.get(0)).plus(ca.exp(x.get(1))).plus(ca.exp(x.get(2))))));
    Function f = new Function("f", List.of(x), List.of(ca.logsumexp(x)));

    checkfunction(f, f_ref, dms(dm(1.1, 1.3, 1.7)));
    checkCodegen(f, dms(dm(1.1, 1.3, 1.7)), Map.of());
    checkfunction(f, f_ref, dms(dm(1.1, 1.3, 1.3)));
    checkfunction(f, f_ref, dms(dm(1.3, 1.3, 1.3)));

    checkarray(ca.logsumexp(dm(1.3, 1.3, 1.3)), (DM) f_ref.invoke(dm(1.3, 1.3, 1.3)));
    checkarray(ca.logsumexp(dm(1.1, 1.3, 1.3)), (DM) f_ref.invoke(dm(1.1, 1.3, 1.3)));
    checkarray(ca.logsumexp(dm(1.1, 1.3, 1.7)), (DM) f_ref.invoke(dm(1.1, 1.3, 1.7)));

    DM res = (DM) f.invoke(dm(100, 1000, 10000));
    checkarray(res, 10000);
    checkarray(ca.logsumexp(dm(100, 1000, 10000)), (DM) f.invoke(dm(100, 1000, 10000)));

    for (DM v : List.of(dm(0, 1000, 1), dm(1000, 0, 1), dm(0, 1, 1000), dm(5, 1e5, 3, 2, 1e5 - 1))) {
      double[] a = v.nonzeros();
      double max = Double.NEGATIVE_INFINITY, sum = 0;
      for (double ai : a) max = Math.max(max, ai);
      for (double ai : a) sum += Math.exp(ai - max);
      double ref = max + Math.log(sum);
      checkarray(ca.logsumexp(v), ref);
      MX xm = MX.sym("x", v.size1());
      Function fm = new Function("fm", List.of(xm), List.of(ca.logsumexp(xm)));
      checkarray((DM) fm.invoke(v), ref);
      checkCodegen(fm, dms(v), Map.of());
    }
  }

  @Test
  @SuppressWarnings("unchecked")
  public void test_extract_parametric_call_sx() throws IOException {
    MX xm = MX.sym("x");
    MX ym = MX.sym("y");
    MX zm = MX.sym("z");
    Function f = new Function("f", List.of(xm, ym, zm), List.of(xm.times(ym), ym.times(zm), xm.times(zm)),
        Map.of("never_inline", true));

    SX x = SX.sym("x");
    SX p = SX.sym("p");
    List<SX> exprs = f.call(List.of(new SX(1), ca.sin(p), ca.sqrt(p)));
    SX expr = x.times(exprs.get(0).plus(exprs.get(2)));

    List<Object> r = ca.extract_parametric(expr, p);
    SX expr_ret = (SX) r.get(0);
    List<SX> symbols = (List<SX>) r.get(1), parametric = (List<SX>) r.get(2);
    System.out.println("expr_ret " + expr_ret);
    System.out.println("symbols " + symbols);
    System.out.println("parametric " + parametric);
    assertTrue(expr_ret.toString().contains("x*e_0"));
    assertFalse(ca.depends_on(expr_ret, p));
    SX expr_recreated = ca.substitute(List.of(expr_ret), symbols, parametric).get(0);
    System.out.println(expr_recreated.minus(expr));
    assertTrue(ca.cse(expr_recreated.minus(expr)).is_zero());

    exprs = f.call(List.of(new SX(1), ca.sin(p), x.power(2)));
    expr = exprs.get(0).plus(exprs.get(2));
    r = ca.extract_parametric(expr, p);
    expr_ret = (SX) r.get(0);
    symbols = (List<SX>) r.get(1);
    parametric = (List<SX>) r.get(2);
    assertFalse(ca.depends_on(expr_ret, p));
    expr_recreated = ca.substitute(List.of(expr_ret), symbols, parametric).get(0);

    File dir = Files.createTempDirectory("casadi_java_sx").toFile();
    // f.generate('f1.c'), but into a temporary directory
    CodeGenerator cg = new CodeGenerator("f1");
    cg.add(new Function("f", List.of(x, p), List.of(expr_recreated)));
    cg.generate(dir.getPath() + File.separator);
    cg = new CodeGenerator("f2");
    cg.add(new Function("f", List.of(x, p), List.of(expr)));
    cg.generate(dir.getPath() + File.separator);
    deleteDir(dir);
    ca.cse(expr_recreated);
    assertTrue(ca.cse(expr_recreated.minus(expr)).is_zero());
  }

  @Test
  public void test_check_recursion() {
    SX x = SX.sym("x");
    Function f = new Function("f", List.of(x), List.of(x.power(2)), Map.of("always_inline", true));
    MX y = MX.sym("y");
    f.invoke(y);
  }

  @Test
  public void test_sx_eval_mx() {
    int n = 2;
    MX y = MX.sym("y", Sparsity.lower(n));
    Function g = new Function("g", List.of(y), List.of(ca.sin(y)), Map.of("never_inline", true));
    SX x = SX.sym("x", Sparsity.lower(n));
    Function f = new Function("f", List.of(x), List.of(((SX) g.invoke(x.power(2))).minus(2)));
    MX X = MX.sym("X", Sparsity.upper(n));
    Function F1 = new Function("F1", List.of(X), f.call(List.of(X), true));
    DM.rng(1);
    DM A = DM.rand(n, n);
    Function F2 = new Function("F2", List.of(X), f.call(List.of(X)));
    checkfunctionLight(F1, F2, dms(A));
  }

  @Test
  public void test_contains() {
    SX x = SX.sym("x");
    SX y = SX.sym("y");
    SX z = SX.sym("z");
    SX e = y.times(z);
    assertTrue(ca.contains(List.of(x, y, z), x));
    assertFalse(ca.contains(List.of(x, y), z));
    assertTrue(ca.contains_any(List.of(x, y), List.of(y, z)));
    assertFalse(ca.contains_all(List.of(x, y), List.of(y, z)));
    assertTrue(ca.contains_any(List.of(x, y), List.of(x, y)));
    assertTrue(ca.contains_all(List.of(x, y), List.of(x, y)));
    assertTrue(ca.contains(List.of(e, x), e));
    assertInException("Can only convert 1-by-1 matrices to scalars", () -> ca.contains(List.of(vcat(x, y)), x));
  }

  @Test
  public void test_pow() {
    SX x = SX.sym("x");
    SX y = new SX(4, 1);
    assertEquals(4, y.power(0).nnz());
    y.set(1, x);
    assertEquals(4, y.power(0).nnz());
    y.set(1, 0);
    assertEquals(4, y.power(0).nnz());
    y.set(1, 1);
    assertEquals(4, y.power(0).nnz());
  }

  @Test
  public void test_set_precision() {
    double pi_val = 3.141592653589793239;
    String expected_default = "3.14159";
    String expected_prec12 = "3.14159265359";

    DM xd = new DM(pi_val);
    try {
      assertEquals(expected_default, xd.toString());
      DM.set_precision(12);
      assertEquals(expected_prec12, xd.toString());
      DM.set_precision(4);
      DM.set_scientific(true);
      assertTrue(xd.toString().toLowerCase().contains("e"));
    } finally {
      DM.set_precision(6);
      DM.set_scientific(false);
    }
    SX xs = new SX(pi_val);
    try {
      assertEquals(expected_default, xs.toString());
      SX.set_precision(12);
      assertEquals(expected_prec12, xs.toString());
      SX.set_precision(4);
      SX.set_scientific(true);
      assertTrue(xs.toString().toLowerCase().contains("e"));
    } finally {
      SX.set_precision(6);
      SX.set_scientific(false);
    }

    MX mx = new MX(pi_val);
    try {
      assertEquals(expected_default, mx.toString());
      DM.set_precision(12);
      assertEquals(expected_prec12, mx.toString());
    } finally {
      DM.set_precision(6);
    }

    DM vd = dm(1.234567890123, 9.876543210987);
    try {
      assertEquals("[1.23457, 9.87654]", vd.toString());
      DM.set_precision(12);
      assertEquals("[1.23456789012, 9.87654321099]", vd.toString());
    } finally {
      DM.set_precision(6);
    }
    SX vs = new SX(new double[] {1.234567890123, 9.876543210987});
    try {
      assertEquals("[1.23457, 9.87654]", vs.toString());
      SX.set_precision(12);
      assertEquals("[1.23456789012, 9.87654321099]", vs.toString());
    } finally {
      SX.set_precision(6);
    }
  }

  @Test
  public void test_linearize() {
    SX x = SX.sym("x");
    SX y = ca.sin(x);
    double x0 = 0.2;
    System.out.println(ca.linearize(ca.sin(x), x, new SX(x0)));
    System.out.println(ca.plus(Math.sin(x0), ca.times(Math.cos(x0), x.minus(x0))));
    Function F = new Function("F", List.of(x), List.of(y));
    Function Ff = F.forward(1);
    DM n = (DM) F.invoke(x0);
    System.out.println(ca.plus(n, (SX) Ff.invoke(x0, n, x.minus(x0))));
    System.out.println(ca.taylor(y, x, x0));
  }

  // test_numpy_style_methods: python-only (%pythoncode method spellings x.sin(), round(x), divmod)

  @Test
  public void test_gauss_quadrature() {
    SX x = SX.sym("x");
    checkarray(ca.evalf(ca.gauss_quadrature(x.power(2), x, 0, 1)), new DM(1. / 3));
    checkarray(ca.evalf(ca.gauss_quadrature(x.power(9), x, -1, 1)), new DM(0));
  }

  @Test
  public void test_poly_roots_order() {
    assertRaisesRegex("order at most 4.*Got order 5", () -> ca.poly_roots(new SX(new double[] {1, 2, 3, 4, 5, 6})));
  }

  @Test
  public void test_mtaylor_dimension() {
    SX x = SX.sym("x");
    SX y = SX.sym("y");
    assertRaisesRegex("expected a of shape 2x1", () -> ca.mtaylor(x.times(y), vcat(x, y), 0, 1));
  }

  @Test
  public void test_large_whole_constants() {
    for (double v : new double[] {3e9, 6.02214076e23, -1e20, Math.pow(2, 31), -Math.pow(2, 31) - 1, Math.pow(2, 63),
        1e300}) {
      assertEquals(v, ca.evalf(new SX(v)).to_double());
      assertEquals(v, ca.evalf(new MX(v)).to_double());
    }
    for (double v : new double[] {Math.pow(2, 31) - 1, -Math.pow(2, 31)}) {
      assertTrue(new SX(v).is_integer());
      assertEquals(v, ca.evalf(new SX(v)).to_double());
    }
  }
}

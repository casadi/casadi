using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;
using System.Text.RegularExpressions;

namespace CasADi.Tests {

/// <summary>Port of test/java SXTest (test/python/sx.py)</summary>
public class SXTest : CasadiTestCase {

  sealed class UnaryOp {
    public readonly Func<SX, SX> cas;
    public readonly Func<double, double> num;
    public readonly string name;
    public readonly bool nozero;
    public UnaryOp(Func<SX, SX> cas, Func<double, double> num, string name, bool nozero) {
      this.cas = cas;
      this.num = num;
      this.name = name;
      this.nozero = nozero;
    }
  }

  sealed class BinaryOp {
    public readonly Func<SX, SX, SX> cas;
    public readonly Func<DM, DM, DM> casDM;
    public readonly Func<double[,], double[,], double[,]> num;
    public readonly string name;
    public BinaryOp(Func<SX, SX, SX> cas, Func<DM, DM, DM> casDM, Func<double[,], double[,], double[,]> num, string name) {
      this.cas = cas;
      this.casDM = casDM;
      this.num = num;
      this.name = name;
    }
  }

  readonly List<UnaryOp> pool = new List<UnaryOp>();
  readonly List<UnaryOp> Jpool = new List<UnaryOp>();
  readonly List<BinaryOp> matrixbinarypool = new List<BinaryOp>();

  static Dictionary<string, object> D(params (string, object)[] kv) => kv.ToDictionary(e => e.Item1, e => e.Item2);

  public SXTest() {
    pool.Add(new UnaryOp(ca.sqrt, Math.Sqrt, "sqrt", false));
    pool.Add(new UnaryOp(ca.sin, Math.Sin, "sin", false));
    pool.Add(new UnaryOp(ca.cos, Math.Cos, "cos", false));
    pool.Add(new UnaryOp(ca.tan, Math.Tan, "tan", false));
    pool.Add(new UnaryOp(ca.abs, Math.Abs, "abs", false));
    pool.Add(new UnaryOp(ca.sign, v => Math.Sign(v), "sign", false));
    pool.Add(new UnaryOp(ca.atan, Math.Atan, "arctan", false));
    pool.Add(new UnaryOp(ca.asin, Math.Asin, "arcsin", false));
    pool.Add(new UnaryOp(ca.acos, Math.Acos, "arccos", false));
    pool.Add(new UnaryOp(ca.exp, Math.Exp, "exp", false));
    pool.Add(new UnaryOp(ca.log, Math.Log, "log", false));
    pool.Add(new UnaryOp(x => ca.pow(x, 0.0), x => Math.Pow(x, 0), "x^0", true));
    pool.Add(new UnaryOp(x => ca.pow(x, 1.0), x => Math.Pow(x, 1), "^1", false));
    pool.Add(new UnaryOp(x => ca.pow(x, -2.0), x => Math.Pow(x, -2), "^-2", true));
    pool.Add(new UnaryOp(x => ca.pow(x, 0.3), x => Math.Pow(x, 0.3), "^0.3", false));
    pool.Add(new UnaryOp(ca.floor, Math.Floor, "floor", false));
    pool.Add(new UnaryOp(ca.ceil, Math.Ceiling, "ceil", false));
    pool.Add(new UnaryOp(ca.log1p, double.LogP1, "log1p", false));
    pool.Add(new UnaryOp(ca.expm1, double.ExpM1, "log1p", false));
    // num is the derivative; the reference Jacobian is its diagonal matrix
    Jpool.Add(new UnaryOp(ca.sqrt, x => 1 / (2.0 * Math.Sqrt(x)), "sqrt", false));
    Jpool.Add(new UnaryOp(ca.sin, Math.Cos, "sin", false));
    Jpool.Add(new UnaryOp(ca.abs, x => Math.Sign(x), "abs", false));
    Jpool.Add(new UnaryOp(ca.sign, x => x * 0, "abs", false));
    Jpool.Add(new UnaryOp(ca.cos, x => -Math.Sin(x), "cos", false));
    Jpool.Add(new UnaryOp(ca.tan, x => 1.0 / Math.Pow(Math.Cos(x), 2), "tan", false));
    Jpool.Add(new UnaryOp(ca.atan, x => 1.0 / (x * x + 1), "arctan", false));
    Jpool.Add(new UnaryOp(ca.asin, x => 1.0 / Math.Sqrt(1 - x * x), "arcsin", false));
    Jpool.Add(new UnaryOp(ca.acos, x => -1.0 / Math.Sqrt(1 - x * x), "arccos", false));
    Jpool.Add(new UnaryOp(ca.exp, Math.Exp, "exp", false));
    Jpool.Add(new UnaryOp(ca.log, x => 1.0 / x, "log", false));
    Jpool.Add(new UnaryOp(x => ca.pow(x, 0.0), x => 0, "x^0", false));
    Jpool.Add(new UnaryOp(x => ca.pow(x, 1.0), x => 1, "^1", false));
    Jpool.Add(new UnaryOp(x => ca.pow(x, -2.0), x => -2.0 / Math.Pow(x, 3), "^-2", false));
    Jpool.Add(new UnaryOp(x => ca.pow(x, 0.3), x => 0.3 / Math.Pow(x, 0.7), "^0.3", false));
    Jpool.Add(new UnaryOp(ca.log1p, x => 1.0 / (1 + x), "log1p", false));
    Jpool.Add(new UnaryOp(ca.expm1, Math.Exp, "log1p", false));
    matrixbinarypool.Add(new BinaryOp((a, b) => a + b, (a, b) => a + b, (a, b) => elementwise(a, b, (u, v) => u + v), "Matrix+Matrix"));
    matrixbinarypool.Add(new BinaryOp((a, b) => a - b, (a, b) => a - b, (a, b) => elementwise(a, b, (u, v) => u - v), "Matrix-Matrix"));
    matrixbinarypool.Add(new BinaryOp((a, b) => a * b, (a, b) => a * b, (a, b) => elementwise(a, b, (u, v) => u * v), "Matrix*Matrix"));
    matrixbinarypool.Add(new BinaryOp(ca.fmax, ca.fmax, (a, b) => elementwise(a, b, Math.Max), "fmin"));
    matrixbinarypool.Add(new BinaryOp(ca.fmin, ca.fmin, (a, b) => elementwise(a, b, Math.Min), "fmax"));
    matrixbinarypool.Add(new BinaryOp(ca.hypot, ca.hypot, (a, b) => elementwise(a, b, double.Hypot), "hypot"));
    matrixbinarypool.Add(new BinaryOp((a, b) => ca.mtimes(a, b.T()), (a, b) => ca.mtimes(a, b.T()), timesTransposed, "dot(Matrix,Matrix.T)"));
  }

  static double[,] elementwise(double[,] a, double[,] b, Func<double, double, double> op) {
    var r = new double[a.GetLength(0), a.GetLength(1)];
    for (int i = 0; i < a.GetLength(0); i++)
      for (int j = 0; j < a.GetLength(1); j++) r[i, j] = op(a[i, j], b[i, j]);
    return r;
  }

  static double[,] timesTransposed(double[,] a, double[,] b) {
    var r = new double[a.GetLength(0), b.GetLength(0)];
    for (int i = 0; i < a.GetLength(0); i++)
      for (int j = 0; j < b.GetLength(0); j++)
        for (int k = 0; k < a.GetLength(1); k++) r[i, j] += a[i, k] * b[j, k];
    return r;
  }

  static DM bools(params int[] v) => v.Select(e => (double)e).ToArray();

  static DM range(int n) => Enumerable.Range(0, n).Select(i => (double)i).ToArray();

  static void assertBools(List<bool> actual, params int[] expected) {
    assertEquals(expected.Length, actual.Count);
    for (int i = 0; i < expected.Length; i++) assertEquals(expected[i] != 0, actual[i], "entry " + i);
  }

  static Function jacobianOld(Function f, int i, int j) {
    var outs = new List<string> { "jac:" + f.name_out(j) + ":" + f.name_in(i) };
    outs.AddRange(f.name_out());
    return f.factory(f.name() + "_jac", f.name_in(), outs);
  }

  static Function hessianOld(Function f, int i, int j) {
    var outs = new List<string> { "hess:" + f.name_out(j) + ":" + f.name_in(i) + ":" + f.name_in(i),
                                  "grad:" + f.name_out(j) + ":" + f.name_in(i) };
    outs.AddRange(f.name_out());
    return f.factory(f.name() + "_hess", f.name_in(), outs);
  }

  void evaluationCheck(List<SX> yt, DM yr, List<SX> x, List<DM> setx0, string name, Func<Function, List<SX>, Function> fmod) {
    message(":" + name);
    var f = ca.Function("f", x, yt);
    if (fmod != null) f = fmod(f, x);
    var f_in = new List<DM>();
    for (int i = 0; i < f.n_in(); i++) f_in.Add(i < setx0.Count ? setx0[i] : ca.DM(0.0));
    List<DM> f_out = f.call(f_in);
    checkarray(ca.densify(f_out[0]), yr, name);
  }

  void numpyEvaluationCheckPool(List<UnaryOp> p, SX x, DM x0, string name, Func<Function, List<SX>, Function> fmod,
                                DM setx0, bool excludeNozero, bool jacobian) {
    foreach (UnaryOp op in p) {
      if (excludeNozero && op.nozero) continue;
      double[,] a = x0.ToArray();
      int n1 = a.GetLength(0), n2 = a.GetLength(1), n = n1 * n2;
      double[,] r = jacobian ? new double[n, n] : new double[n1, n2];
      for (int j = 0, k = 0; j < n2; j++)
        for (int i = 0; i < n1; i++, k++) {
          double v = op.num(a[i, j]);
          if (jacobian) r[k, k] = v;
          else r[i, j] = v;
        }
      evaluationCheck(new List<SX> { op.cas(x) }, r, new List<SX> { x }, new List<DM> { setx0 ?? x0 }, name + ":" + op.name, fmod);
    }
  }

  void numpyEvaluationCheckBinaryPool(SX x, SX y, DM x0, DM y0, DM setx0, DM sety0, string name) {
    foreach (BinaryOp op in matrixbinarypool)
      evaluationCheck(new List<SX> { op.cas(x, y) }, op.num(x0.ToArray(), y0.ToArray()), new List<SX> { x, y },
                      new List<DM> { setx0, sety0 }, name + ":" + op.name, null);
  }

  static void checkSerialize(Function F, List<DM> inputs) {
    Function F2 = Function.deserialize(F.serialize(D(("debug", true))));
    List<DM> Fout = F.call(inputs), Fout2 = F2.call(inputs);
    for (int i = 0; i < F.n_out(); i++) checkarray(Fout[i], Fout2[i], 1e-16, "serialize");
  }

  /// <summary>As check_codegen in helpers.py: generate, compile with gcc, load as external, compare; only when slow</summary>
  static void checkCodegen(Function F, List<DM> inputs, Dictionary<string, object> opts, string std = "c89") {
    if (!RunSlow) return;
    string dir = Path.Combine(Path.GetTempPath(), "casadi_dotnet_codegen_" + Guid.NewGuid().ToString("N"));
    Directory.CreateDirectory(dir);
    string name = "codegen_" + Guid.NewGuid().ToString("N").Substring(0, 16);
    var cg = ca.CodeGenerator(name, opts);
    cg.add(F);
    cg.generate(dir + "/");
    string libdir = GlobalOptions.getCasadiPath(), inc = GlobalOptions.getCasadiIncludePath();
    string lib = Path.Combine(dir, name + ".so");
    var cmd = new List<string> { "-pedantic", "-std=" + std, "-fPIC", "-shared", "-Wall",
      "-Werror", "-Wextra", "-I" + inc, "-I" + inc + "/highs", "-Wno-dangling-pointer",
      "-Wno-unknown-warning-option", "-Wno-unknown-pragmas", "-Wno-long-long", "-Wno-unused-parameter", "-O3",
      Path.Combine(dir, name + ".c"), "-o", lib, "-L" + libdir, "-Wl,-rpath," + libdir };
    var psi = new ProcessStartInfo("gcc");
    foreach (var a in cmd) psi.ArgumentList.Add(a);
    using (var p = Process.Start(psi)) {
      p.WaitForExit();
      assertEquals(0, p.ExitCode, "compile " + string.Join(" ", cmd));
    }
    Function F2 = ca.external(F.name(), lib);
    List<DM> Fout = F.call(inputs), Fout2 = F2.call(inputs);
    for (int i = 0; i < F.n_out(); i++) checkarray(Fout[i], Fout2[i], 1e-15, "codegen");
    Directory.Delete(dir, true);
  }

  static void checkfunctionLight(Function f, Function g, List<DM> inputs) {
    List<DM> a = f.call(inputs), b = g.call(inputs);
    assertEquals(g.n_out(), f.n_out());
    for (int i = 0; i < a.Count; i++) checkarray(a[i], b[i], 1e-9, "output " + i);
  }

  static void assertRaisesRegex(string regex, Action r) {
    try {
      r();
    } catch (Exception e) when (!(e is AssertFailed)) {
      assertTrue(Regex.IsMatch(e.Message ?? "", regex, RegexOptions.Singleline), "expected /" + regex + "/ in: " + e.Message);
      return;
    }
    throw new AssertFailed("expected an exception matching /" + regex + "/");
  }

  public void test_equivalence() {
    SX x = SX.sym("x");
    SX y = SX.sym("y");
    SX[][] cases = { new[] { ca.log1p(x), ca.log(x + 1) }, new[] { ca.expm1(x), ca.exp(x) - 1 },
                     new[] { ca.hypot(x, y), ca.sqrt(ca.pow(x, 2) + ca.pow(y, 2)) } };
    foreach (SX[] c in cases) {
      var f = ca.Function("f", new[] { x, y }, new[] { c[0] });
      var equiv_f = ca.Function("equiv_f", new[] { x, y }, new[] { c[1] });
      checkfunction(f, equiv_f, dms(1.3, 1.7));
    }
  }

  public void test_scalarSX() {
    SX x = SX.sym("x");
    DM x0 = 0.738;
    numpyEvaluationCheckPool(pool, x, x0, "scalarSX", null, null, false, false);
  }

  public void test_gradient() {
    message("jacobian of SX**number");
    SX x = SX.sym("x");
    double x0 = 1;
    int p = 3;
    SX y = ca.pow(x, p);
    SX dx = ca.jacobian(y, x);
    double dxr = p;
    evaluationCheck(new List<SX> { dx }, dxr, new List<SX> { x }, new List<DM> { x0 }, "jacobian", null);
    dxr = 1;
    for (int i = 0; i < p; i++) {
      y = ca.jacobian(y, x);
      dxr = dxr * (p - i);
    }
    evaluationCheck(new List<SX> { y }, dxr, new List<SX> { x }, new List<DM> { x0 }, "recursive jacobian", null);
  }

  public void test_gradient2() {
    message("jacobian of SX**SX");
    SX x = SX.sym("x");
    SX p = SX.sym("p");
    double x0 = 1;
    int p0 = 3;
    SX y = ca.pow(x, p);
    SX dx = ca.jacobian(y, x);
    double dxr = p0;
    evaluationCheck(new List<SX> { dx }, dxr, new List<SX> { x, p }, new List<DM> { x0, p0 }, "jacobian", null);
    dxr = 1;
    for (int i = 0; i < p0; i++) {
      y = ca.jacobian(y, x);
      dxr = dxr * (p0 - i);
    }
    evaluationCheck(new List<SX> { y }, dxr, new List<SX> { x, p }, new List<DM> { x0, p0 }, "jacobian", null);
  }

  public void test_SXJacobian() {
    message("SX(1,1) unary operation, jacobian");
    SX x = SX.sym("x");
    DM x0 = ca.DM(new double[,] { { 0.738 } });
    numpyEvaluationCheckPool(Jpool, x, x0, "SX unary operations, jacobian", (f, xx) => jacobianOld(f, 0, 0), null, false, true);
  }

  public void test_SXJac() {
    message("SX(1,1) unary operation, jac");
    SX x = SX.sym("x");
    DM x0 = ca.DM(new double[,] { { 0.738 } });
    numpyEvaluationCheckPool(Jpool, x, x0, "SX unary operations, jac", (f, xx) => {
      List<SX> y = f.call(xx);
      return ca.Function("J", xx, new[] { ca.jacobian(y[0], xx[0]) });
    }, null, false, true);
  }

  public void test_SXJacobians() {
    message("SX(3,1) unary operation, jacobian");
    SX x = SX.sym("x", 3);
    DM x0 = ca.vertcat(0.738, 0.9, 0.3);
    numpyEvaluationCheckPool(Jpool, x, x0, "SX unary operations, jacobian", (f, xx) => jacobianOld(f, 0, 0), null, false, true);
  }

  public void test_SXJacobians2() {
    message("SX(1,3) unary operation, jacobian");
    SX x = SX.sym("x", 1, 3);
    DM x0 = ca.vertcat(0.738, 0.9, 0.3);
    numpyEvaluationCheckPool(Jpool, x, x0, "SX unary operations, jacobian", (f, xx) => jacobianOld(f, 0, 0), null, false, true);
  }

  public void test_SX() {
    message("SX unary operations");
    SX x = SX.sym("x", 3, 2);
    DM x0 = ca.DM(new double[,] { { 0.738, 0.2 }, { 0.1, 0.39 }, { 0.99, 0.999999 } });
    numpyEvaluationCheckPool(pool, x, x0, "SX", null, null, false, false);

    x = SX.sym("x", 3, 3);
    double[,] a = { { 0.738, 0.2, 0.3 }, { 0.1, 0.39, -6 }, { 0.99, 0.999999, -12 } };
    x0 = ca.DM(a);
    double det_ = a[0, 0] * (a[1, 1] * a[2, 2] - a[1, 2] * a[2, 1]) - a[0, 1] * (a[1, 0] * a[2, 2] - a[1, 2] * a[2, 0])
                  + a[0, 2] * (a[1, 0] * a[2, 1] - a[1, 1] * a[2, 0]);
    var invr = new double[3, 3];
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++) {
        int r1 = (j + 1) % 3, r2 = (j + 2) % 3, c1 = (i + 1) % 3, c2 = (i + 2) % 3;
        invr[i, j] = (a[r1, c1] * a[r2, c2] - a[r1, c2] * a[r2, c1]) / det_;
      }
    message(":det(SX)");
    evaluationCheck(new List<SX> { ca.SX(ca.det(x)) }, det_, new List<SX> { x }, new List<DM> { x0 }, "det(SX)", null);
    message(":inv(SX)");
    evaluationCheck(new List<SX> { ca.inv(x) }, invr, new List<SX> { x }, new List<DM> { x0 }, "inv(SX)", null);
  }

  public void test_SXSparse() {
    message("SX unary operations, sparse");
    SX x = SX.sym("x");
    SX y = SX.sym("y");
    SX z = SX.sym("z");
    var sp = ca.Sparsity(4, 3, new long[] { 0, 2, 2, 3 }, new long[] { 1, 2, 1 });
    x = ca.SX(sp, ca.vertcat(x, y, z));
    DM x0 = ca.DM(sp, ca.vertcat(0.738, 0.1, 0.99));
    numpyEvaluationCheckPool(pool, x, ca.densify(x0), "SX", null, x0, true, false);
  }

  public void test_SXbinary() {
    message("SX binary operations");
    SX x = SX.sym("x", 3, 2);
    SX y = SX.sym("x", 3, 2);
    DM x0 = ca.DM(new double[,] { { 0.738, 0.2 }, { 0.1, 0.39 }, { 0.99, 0.999999 } });
    DM y0 = ca.DM(new double[,] { { 1.738, 0.6 }, { 0.7, 12 }, { 0, -6 } });
    numpyEvaluationCheckBinaryPool(x, y, x0, y0, x0, y0, "SX");
    Check.Throws<Exception>(() => ca.mtimes(x, y));
  }

  public void test_SXbinary_codegen() {
    Check.Assume(RunSlow, "slow");
    message("SX binary operations");
    SX x = SX.sym("x", 4, 2);
    SX y = SX.sym("x", 4, 2);
    DM x0 = ca.DM(new double[,] { { 0.738, 0.2 }, { 0.1, 0.39 }, { 0.99, 0.999999 }, { 1, 2 } });
    DM y0 = ca.DM(new double[,] { { 1.738, 0.6 }, { 0.7, 12 }, { 0, -6 }, { 1, 2 } });
    foreach (BinaryOp op in matrixbinarypool) {
      var f = ca.Function("f", new[] { x, y }, new[] { op.cas(x, y) });
      checkCodegen(f, dms(x0, y0), D(), "c99");
      checkCodegen(f, dms(x0, y0), D(), "c89");
    }
  }

  public void test_SXbinary_diff() {
    message("SX binary operations");
    SX x = SX.sym("x", 4, 2);
    SX y = SX.sym("x", 4, 2);
    SX dx = SX.sym("x", 4, 2);
    SX dy = SX.sym("x", 4, 2);
    DM x0 = ca.DM(new double[,] { { 0.738, 0.2 }, { 0.1, 0.39 }, { 0.99, 0.999999 }, { 1, 2 } });
    DM y0 = ca.DM(new double[,] { { 1.738, 0.6 }, { 0.7, 12 }, { 0, -6 }, { 1, 2 } });
    foreach (BinaryOp op in matrixbinarypool) {
      string name = op.name;
      var f = ca.Function("f", new[] { x, y }, new[] { op.cas(x, y) });
      f.disp(true);
      SX fxy = f.Call(x, y);
      var df = ca.Function("f", new[] { x, y, dx, dy }, new[] { ca.jtimes(fxy, x, dx) + ca.jtimes(fxy, y, dy) });
      double eps = 1e-7;
      DM.rng(0);
      DM dx0 = DM.rand(x.sparsity());
      DM dy0 = name == "fmin" || name == "fmax" ? dx0 : DM.rand(y.sparsity());
      dy0 = name == "fmin" || name == "fmax" ? dx0 : DM.rand(y.sparsity());
      DM df_experiment = (f.Call(x0 + dx0 * eps, y0 + dy0 * eps) - f.Call(x0, y0)) / eps;
      DM df0 = df.Call(x0, y0, dx0, dy0);
      Console.WriteLine(df_experiment + " " + df0);
      checkarray(df0, df_experiment, 0.5e-4, name);
    }
  }

  public void test_DMbinary() {
    message("SX binary operations");
    DM x0 = ca.DM(new double[,] { { 0.738, 0.2 }, { 0.1, 0.39 }, { 0.99, 0.999999 } });
    DM y0 = ca.DM(new double[,] { { 1.738, 0.6 }, { 0.7, 12 }, { 0, -6 } });
    DM a = ca.vertcat(x0, y0);
    foreach (BinaryOp op in matrixbinarypool) {
      double[,] a0 = { { (double)a[0] } }, a1 = { { (double)a[1] } };
      checkarray(op.casDM(a[0], a[1]), op.num(a0, a1), op.name);
    }
  }

  public void test_SXbinarySparse() {
    message("SX binary operations");
    SX x = SX.sym("x");
    SX y = SX.sym("y");
    SX z = SX.sym("z");
    SX x2 = SX.sym("x2");
    SX y2 = SX.sym("y2");
    SX z2 = SX.sym("z2");
    var spx = ca.Sparsity(4, 3, new long[] { 0, 2, 2, 3 }, new long[] { 1, 2, 1 });
    var spy = ca.Sparsity(4, 3, new long[] { 0, 2, 2, 3 }, new long[] { 0, 2, 3 });
    var xx = ca.SX(spx, ca.vertcat(x, y, z));
    var yy = ca.SX(spy, ca.vertcat(x2, z2, y2));
    DM x0 = ca.DM(spx, ca.vertcat(0.738, 0.1, 0.99));
    DM y0 = ca.DM(spy, ca.vertcat(1.738, 0.7, -6));
    numpyEvaluationCheckBinaryPool(xx, yy, ca.densify(x0), ca.densify(y0), x0, y0, "SX");
    Check.Throws<Exception>(() => ca.mtimes(xx, yy));
  }

  void slicingCheck(Func<SX, SX> ft, Func<DM, DM> fr, SX x, DM x0, DM setx0, string name) =>
    evaluationCheck(new List<SX> { ft(x) }, fr(x0), new List<SX> { x }, new List<DM> { setx0 }, name, null);

  public void test_SXslicing() {
    Check.Skip("known bug in Python too: test refactoring, cf. #1436");
    message("SX slicing/indexing");
    SX x = SX.sym("x", 3, 2);
    DM x0 = ca.DM(new double[,] { { 0.738, 0.2 }, { 0.1, 0.39 }, { 0.99, 0.999999 } });

    message(":dense");
    foreach (DM s0 in new DM[] { null }) {
      DM sx = s0 ?? x0;
      slicingCheck(v => v[0, 0], v => v[0, 0], x, x0, sx, "x[0,0]");
      slicingCheck(v => v[1, 0], v => v[1, 0], x, x0, sx, "x[1,0]");
      slicingCheck(v => v[0, 1], v => v[0, 1], x, x0, sx, "x[1,0]");
      slicingCheck(v => v[0, -1], v => v[0, 1], x, x0, sx, "x[0,-1]");
      slicingCheck(v => v[.., 0], v => v[.., 0], x, x0, sx, "x[:,0]");
      slicingCheck(v => v[.., 1], v => v[.., 1], x, x0, sx, "x[:,1]");
      slicingCheck(v => v[1, ..], v => v[1, ..], x, x0, sx, "x[1,:]");
      slicingCheck(v => v[0, ..], v => v[0, ..], x, x0, sx, "x[0,:]");
      slicingCheck(v => v[-1, ..], v => v[2, ..], x, x0, sx, "x[-1,:]");
      slicingCheck(v => v[.., -2], v => v[.., 0], x, x0, sx, "x[:,-2]");
      slicingCheck(v => v[ca.Slice(0, -2), ca.Slice(0, -1)], v => v[ca.Slice(0, 1), ca.Slice(0, 1)], x, x0, sx, "x[0:-2,0:-1]");
      slicingCheck(v => v[ca.Slice(0, 2), ca.Slice(0, 2)], v => v[ca.Slice(0, 2), ca.Slice(0, 2)], x, x0, sx, "x[0:2,0:2]");
      slicingCheck(v => v.get(ca.vertcat(0, 1), ca.Slice(0, 2)), v => v[ca.Slice(0, 2), ca.Slice(0, 2)], x, x0, sx, "x[[0,1],0:2]");
      slicingCheck(v => v.get_nz(false, ca.vertcat(0, 2, 3)), v => ca.vertcat((double)v[0, 0], (double)v[2, 0], (double)v[0, 1]), x, x0, sx, "x[[0,2,3]]");
      slicingCheck(v => v.get_nz(false, ca.Slice(0, 2)), v => ca.vertcat((double)v[0, 0], (double)v[0, 1]), x, x0, sx, "x[0:2] on dense matrix");
      slicingCheck(v => v.get_nz(false, (DM)1.0), v => v[0, 1], x, x0, sx, "x[1]");
      slicingCheck(v => v.get_nz(false, (DM)(-1.0)), v => v[2, 1], x, x0, sx, "x[-1]");
    }

    message(":sparse");
    var sp = ca.Sparsity(4, 3, new long[] { 0, 2, 2, 3 }, new long[] { 1, 2, 1 });
    x = ca.SX(sp, ca.vertcat(SX.sym("x"), SX.sym("y"), SX.sym("z")));
    DM sx0 = ca.vertcat(0.738, 0.39, 0.99);
    x0 = ca.densify(ca.DM(sp, sx0));
    slicingCheck(v => v[0, 0], v => v[0, 0], x, x0, sx0, "x[0,0]");
    slicingCheck(v => v[0, 0], v => v[0, 0], x, x0, sx0, "x[0,0]");
    slicingCheck(v => v[1, 0], v => v[1, 0], x, x0, sx0, "x[1,0]");
    slicingCheck(v => v[0, 1], v => v[0, 1], x, x0, sx0, "x[1,0]");
    slicingCheck(v => v[0, -1], v => v[0, 2], x, x0, sx0, "x[0,-1]");
    slicingCheck(v => v[.., 0], v => v[.., 0], x, x0, sx0, "x[:,0]");
    slicingCheck(v => v[.., 1], v => v[.., 1], x, x0, sx0, "x[:,1]");
    slicingCheck(v => v[1, ..], v => v[1, ..], x, x0, sx0, "x[1,:]");
    slicingCheck(v => v[0, ..], v => v[0, ..], x, x0, sx0, "x[0,:]");
    slicingCheck(v => v[-1, ..], v => v[3, ..], x, x0, sx0, "x[-1,:]");
    slicingCheck(v => v[.., -2], v => v[.., 1], x, x0, sx0, "x[:,-2]");
    slicingCheck(v => v[ca.Slice(0, -2), ca.Slice(0, -1)], v => v[ca.Slice(0, 2), ca.Slice(0, 2)], x, x0, sx0, "x[0:-2,0:-1]");
    slicingCheck(v => v[ca.Slice(0, 2), ca.Slice(0, 2)], v => v[ca.Slice(0, 2), ca.Slice(0, 2)], x, x0, sx0, "x[0:2,0:2]");
    slicingCheck(v => v.get(ca.vertcat(0, 1), ca.Slice(0, 2)), v => v[ca.Slice(0, 2), ca.Slice(0, 2)], x, x0, sx0, "x[[0,1],0:2]");
    slicingCheck(v => v.get_nz(false, ca.vertcat(2, 1)), v => ca.vertcat((double)v[1, 2], (double)v[2, 0]), x, x0, sx0, "x[[2,1]]");
    slicingCheck(v => v.get_nz(false, ca.Slice(0, 2)), v => sx0[ca.Slice(0, 2)], x, x0, sx0, "x[0:2] on dense matrix");
    slicingCheck(v => v.get_nz(false, (DM)1.0), v => sx0[1], x, x0, sx0, "x[1]");
    slicingCheck(v => v.get_nz(false, (DM)(-1.0)), v => sx0[2], x, x0, sx0, "x[-1]");
  }

  public void test_SX1() {
    message("SXFunction evaluation");
    SX x = SX.sym("x");
    SX y = SX.sym("y");
    var f = ca.Function("f", new[] { ca.vertcat(x, y) }, new[] { ca.vertcat(x + y, x * y, ca.pow(x, 2) + ca.pow(y, 3)) });
    DM L = ca.vertcat(2, 3);
    List<DM> f_out = f.call(new[] { L });
    List<double> z = f_out[0].nonzeros();
    double[] zr = { 2 + 3, 2 * 3, 2 * 2 + 3 * 3 * 3 };
    for (int i = 0; i < 3; i++) checkarray(z[i], zr[i], 1e-10);
    message("SXFunction jacobian evaluation");
    Function J = jacobianOld(f, 0, 0);
    var J_in = new List<DM> { L };
    for (int i = 1; i < J.n_in(); i++) J_in.Add(0.0);
    List<DM> J_out = J.call(J_in);
    checkarray(J_out[0], new double[,] { { 1, 1 }, { 3, 2 }, { 4, 27 } });
  }

  public void test_SX2() {
    message("SXFunction evalution 2");
    SX x = SX.sym("x");
    SX y = SX.sym("y");
    var f = new List<SX> { 3.0 - ca.sin(x * x) - y, ca.sqrt(y) * x };
    if (GlobalOptions.getSimplificationOnTheFly())
      assertEquals("((3-sin(sq(x)))-y)", f[0].ToString(), "SX representation is wrong " + f[0]);
    else
      assertEquals("((3-sin((x*x)))-y)", f[0].ToString(), "SX representation is wrong " + f[0]);
    assertEquals("(sqrt(y)*x)", f[1].ToString(), "SX representation is wrong " + f[1]);
    var fcn = ca.Function("fcn", new[] { ca.vertcat(x, y) }, new[] { ca.vertcat(f) });
    assertEquals("fcn:(i0[2])->(o0[2]) SXFunction", fcn.ToString(), "SX representation is wrong");
    List<DM> fcn_out = fcn.call(new[] { ca.vertcat(2, 3) });
    List<double> res = fcn_out[0].nonzeros();
    checkarray(res[0], 3 - Math.Sin(2 * 2) - 3, 1e-10);
    checkarray(res[1], Math.Sqrt(3) * 2, 1e-10);
  }

  public void test_const_folding_on_the_fly() {
    SX x = SX.sym("x");
    foreach (int a in new[] { 2, 7 })
      foreach (int b in new[] { 2, 7 }) {
        SX r = (double)(a * b) * x;
        assertEquals(r.ToString(), ((double)a * ((double)b * x)).ToString());
        assertEquals(r.ToString(), ((double)b * x * a).ToString());
      }
  }

  public void test_SX_func() {
    message("Function constructors");
    SX y = SX.sym("y", 2, 3);
    var f = ca.Function("f", new[] { y }, new[] { y });
    assertEquals((2L, 3L), f.size_in(0), "Function constructors");
    assertEquals((2L, 3L), f.size_out(0), "Function constructors");
  }

  public void test_SXconversion() {
    message("Conversions from and to SX");
    SX y = SX.sym("y");
    SX x = SX.sym("x", 3, 3);
    ca.SX(y);
    ca.SX(x);
    ca.det(x);
    var yd = ca.DM(x);
    ca.det(yd);
  }

  public void test_SXbool() {
    message("bool");
    SX x = SX.sym("x");
    SX y = SX.sym("y");
    var f = ca.Function("f", new[] { ca.vertcat(x, y) }, new[] { ca.vertcat(ca.logic_and(x, y), ca.logic_or(x, y), ca.logic_not(x)) });
    foreach (int t1 in new[] { 0, 1 })
      foreach (int t2 in new[] { 0, 1 }) {
        bool T1 = t1 != 0, T2 = t2 != 0;
        List<DM> f_out = f.call(new[] { ca.vertcat(t1, t2) });
        checkarray(f_out[0], bools(T1 && T2 ? 1 : 0, T1 || T2 ? 1 : 0, !T1 ? 1 : 0), "bool(" + t1 + "," + t2 + "): " + f_out[0]);
      }
  }

  public void test_not_not() {
    SX x = SX.sym("x");
    SX y = SX.sym("y");
    var f = ca.Function("f", new[] { x }, new[] { ca.logic_not(ca.logic_not(x)) });
    checkarray(f.Call((DM)3.0), ca.DM(1.0));
    checkarray(f.Call((DM)0.0), ca.DM(0.0));
    foreach (SX b in new[] { ca.lt(x, y), ca.le(x, y), ca.eq(x, y), ca.ne(x, y), ca.logic_not(x), ca.logic_and(x, y), ca.logic_or(x, y) })
      assertTrue(ca.is_equal(ca.logic_not(ca.logic_not(b)), b));
    assertFalse(ca.is_equal(ca.logic_not(ca.logic_not(x)), x));
  }

  public void test_is_nonnegative_binary() {
    SX x = SX.sym("x");
    SX y = SX.sym("y");
    foreach (SX e in new[] { ca.lt(x, y), ca.le(x, y), ca.eq(x, y), ca.ne(x, y), ca.logic_and(x, y), ca.logic_or(x, y), ca.hypot(x, y), ca.sqrt(x) })
      assertTrue(ca.is_equal(ca.abs(e), e));
    foreach (SX e in new[] { x + y, x * y, x - y }) assertFalse(ca.is_equal(ca.abs(e), e));
  }

  public void test_chained_inequality_kept() {
    SX y = SX.sym("y");
    // Python 0<=(y<=1) and (1>=y)>=0
    foreach (SX e in new[] { ca.le(0.0, ca.le(y, 1.0)), ca.ge(ca.ge(1.0, y), 0.0) }) {
      assertFalse(e.is_constant());
      Function f = ca.Function("f", new[] { y }, new[] { e }).transform(new[] { new object[] { "simplify", 0, "ref_count" } });
      assertTrue(f.call(new[] { y }, true, false)[0].is_one());
    }
  }

  public void test_and_zero() {
    SX xs = SX.sym("x");
    assertTrue(ca.logic_and(xs, 0.0).is_zero());
    assertTrue(ca.logic_and(0.0, xs).is_zero());
    MX xm = MX.sym("x");
    assertTrue(ca.logic_and(xm, 0.0).is_zero());
    assertTrue(ca.logic_and(0.0, xm).is_zero());
    assertTrue(ca.logic_and(ca.SX(0.0), ca.SX(0.0)).is_zero());
  }

  public void test_SXineq() {
    message("SX ineq");
    SX x = SX.sym("x");
    SX y = SX.sym("y");
    var f = ca.Function("f", new[] { ca.vertcat(x, y) }, new[] { ca.vertcat(ca.lt(x, y), ca.le(x, y), ca.ge(x, y), ca.eq(x, y), ca.ne(x, y)) });
    double[] vals = { -10, 0.1, 0, 1, 10 };
    foreach (double T1 in vals)
      foreach (double T2 in vals) {
        List<DM> f_out = f.call(new[] { ca.vertcat(T1, T2) });
        checkarray(f_out[0], bools(T1 < T2 ? 1 : 0, T1 <= T2 ? 1 : 0, T1 >= T2 ? 1 : 0, T1 == T2 ? 1 : 0, T1 != T2 ? 1 : 0),
                   "ineq(" + T1 + "," + T2 + ")");
      }
  }

  public void test_SX_func2() {
    message("SXmatrix typemaps constructors");
    ca.simplify(SX.sym("x"));
    message(":number");
    DM i = ca.transpose(ca.transpose((DM)2.3));
    assertEquals(1L, i.size1(), "shape mismatch");
    assertEquals(1L, i.size2(), "shape mismatch");
    ca.SX(2.3).is_empty();
    message(":SX");
    SX s = SX.sym("x");
    SX j = ca.transpose(ca.transpose(s));
    assertEquals(1L, j.size1(), "shape mismatch");
    assertEquals(1L, j.size2(), "shape mismatch");
    ca.SX(s).is_empty();
  }

  public void test_SX_func3() {
    message("vector(SXmatrix) typemaps constructors");
    SX y = SX.sym("y");
    SX x = SX.sym("x", 3, 1);
    ca.vertcat(x, x);
    ca.vertcat(y, y);
    ca.vertcat(x, ca.SX());
  }

  public void test_eval() {
    message("Function eval");
    SX x = SX.sym("x", 2, 2);
    SX y = SX.sym("y", 2, 2);
    var f = ca.Function("f", new[] { x, y }, new[] { x * y });
    f.Call(x, y);
  }

  public void test_symbolcheck() {
    message("Check if non-symbolic inputs are caught");
    Check.Throws<Exception>(() => ca.Function("f", new[] { ca.SX(0.0) }, new[] { SX.sym("x") }));
  }

  public void test_sparseconstr() {
    message("Check sparsity constructors");
    checkarray(ca.densify(DM.ones(Sparsity.lower(3))), new double[,] { { 1, 0, 0 }, { 1, 1, 0 }, { 1, 1, 1 } });
    checkarray(ca.densify(DM.ones(Sparsity.diag(3))), new double[,] { { 1, 0, 0 }, { 0, 1, 0 }, { 0, 0, 1 } });
  }

  public void test_subsassignment() {
    message("Check subscripted assignment");
    var xn = new double[3, 4];
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 4; j++) xn[i, j] = random.NextDouble();
    DM x = xn;

    var y = ca.DM(7, 8);
    var z = new double[7, 8];
    y[0, 0] = 12;
    z[0, 0] = 12;
    checkarray(y, z);
    int[] cols = { 2, 4, 5, 6 };
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 4; j++) z[1 + i, cols[j]] = xn[i, j];
    y.set(x, false, ca.Slice(1, 4), ca.vertcat(2, 4, 5, 6));
    checkarray(y, z);

    long[] kl = { 2, 4, 5, 8 };
    DM kld = ca.vertcat(2, 4, 5, 8);
    y.set_nz(1.0, false, kld);
    Sparsity s = y.sparsity();
    List<long> srow = s.row(), scol = s.get_col();
    foreach (long k in kl) z[srow[(int)k], scol[(int)k]] = 1.0;
    checkarray(y, z);
    y.set_nz(kld, false, kld);
    int cnt = 0;
    foreach (long k in kl) {
      z[srow[(int)k], scol[(int)k]] = kl[cnt];
      cnt += 1;
    }
    checkarray(y, z);
  }

  public void test_substitute() {
    Check.Assume(GlobalOptions.getSimplificationOnTheFly(), "no simplification on the fly");
    message("Basic symbolic algebra: substitute");
    SX x = SX.sym("x");
    SX y = SX.sym("y");
    SX z = ca.cos(x) * y;
    assertTrue(ca.depends_on(z, y));
    assertTrue(ca.depends_on(z, x));
    SX w = ca.substitute(z, x, 0.0);
    assertTrue(w.is_symbolic());
    assertTrue(ca.depends_on(w, y));
    assertFalse(ca.depends_on(w, x));
    assertTrue(ca.is_equal(w, y));
    SX r = w - y;
    assertFalse(r.is_symbolic());
    assertTrue(r.is_zero());
    assertEquals(0.0, r.to_double());
    assertEquals(0.0, r.to_double());
    SX y2 = SX.sym("y", 2);
    y2 = ca.substitute(y2 + 6, y2, 0.0);
    assertEquals(6L, y2[0].to_int());
    assertEquals(6L, y2[1].to_int());
  }

  public void test_primitivefunctions() {
    message("Primitive functions");
    SX x = SX.sym("x");
    double[] nums = { -2, -1.5, -1, -0.5, -0.25, 0, 0.25, 0.5, 1, 1.5, 2 };
    primitiveTest(x, ca.sign, "sign", nums, new double[] { -1, -1, -1, -1, -1, 0, 1, 1, 1, 1, 1 });
    primitiveTest(x, ca.heaviside, "heaviside", nums, new double[] { 0, 0, 0, 0, 0, 0.5, 1, 1, 1, 1, 1 });
    primitiveTest(x, ca.ramp, "ramp", nums, new double[] { 0, 0, 0, 0, 0, 0, 0.25, 0.50, 1, 1.5, 2 });
    primitiveTest(x, ca.rectangle, "rectangle", nums, new double[] { 0, 0, 0, 0.5, 1, 1, 1, 0.5, 0, 0, 0 });
    primitiveTest(x, ca.triangle, "triangle", nums, new double[] { 0, 0, 0, 0.5, 0.75, 1, 0.75, 0.5, 0, 0, 0 });
  }

  void primitiveTest(SX x, Func<SX, SX> fun, string comment, double[] nums, double[] reference) {
    message(":" + comment);
    var f = ca.Function("f", new[] { x }, new[] { fun(x) });
    for (int i = 0; i < nums.Length; i++) {
      List<DM> f_out = f.call(new DM[] { nums[i] });
      double r = (double)f_out[0][0];
      assertTrue(reference[i] == r, reference[i] + " vs " + r);
    }
  }

  public void test_taylor() {
    message("univariate taylor expansion");
    SX x = SX.sym("x");
    if (GlobalOptions.getSimplificationOnTheFly()) assertTrue(ca.is_equal(ca.taylor(ca.sin(x), x), x));
    double a_ = 0.13, x_ = 0.15;
    SX a = SX.sym("a");
    Action<SX, double> test = (e, r) => {
      var f = ca.Function("f", new[] { x, a }, new[] { e });
      List<DM> f_out = f.call(dms(x_, a_));
      checkarray((double)f_out[0][0], r, 1e-10);
    };
    double d = x_ - a_;
    test(ca.taylor(ca.sin(x), x, a, 0), Math.Sin(a_));
    test(ca.taylor(ca.sin(x), x, a, 1), Math.Sin(a_) + Math.Cos(a_) * d);
    test(ca.taylor(ca.sin(x), x, a, 2), Math.Sin(a_) + Math.Cos(a_) * d - (Math.Sin(a_) * d * d) / 2.0);
    test(ca.taylor(ca.sin(x), x, a, 3),
         Math.Sin(a_) + Math.Cos(a_) * d - (Math.Sin(a_) * d * d) / 2.0 - (Math.Cos(a_) * d * d * d) / 6.0);

    SX M = ca.blockcat(new[] { new[] { a * ca.sin(x), a * ca.cos(x) }, new[] { ca.exp(a * x), a * ca.pow(x, 2) }, new[] { ca.cos(x), ca.SX(0.0) } });
    var g = ca.Function("f", new[] { x, a }, new[] { ca.taylor(M, x) });
    List<DM> g_out = g.call(dms(x_, a_));
    checkarray(g_out[0], new double[,] { { x_ * a_, a_ }, { 1 + a_ * x_, 0 }, { 1, 0 } }, 1e-10, "taylor on dense matrices");
  }

  static void checkarray(DM a, double[,] b, double tol, string msg) => checkarray(a, (DM)b, tol, msg);

  public void test_null() {
    message("Function null");
    SX x = SX.sym("x");

    var f = ca.Function("f", new[] { x }, new[] { ca.pow(x, 2), ca.SX() });
    List<DM> f_out = f.call(dms(0.0));
    assertTrue(f_out[1].is_empty());

    f = ca.Function("f", new[] { x, ca.SX() }, new[] { ca.pow(x, 2), ca.SX() });
    f_out = f.call(dms(0.0, 0.0));
    assertTrue(f_out[1].is_empty());
    f_out = f.call(dms(0.0, 0.0));

    List<SX> r = f.call(new[] { x, ca.SX() });
    assertTrue(r[1].is_empty());
    r = f.call(new[] { x, ca.SX() });
    assertTrue(r[1].is_empty());
    r = f.call(new[] { x, ca.SX(0, 1) });
    assertTrue(r[1].is_empty());
    r = f.call(new[] { x, ca.SX(1, 0) });
    assertTrue(r[1].is_empty());
  }

  public void test_mtaylor() {
    message("multivariate taylor expansions");
    SX x = SX.sym("x");
    SX y = SX.sym("y");
    SX a = SX.sym("a");
    SX b = SX.sym("b");
    double a_ = 0.13, x_ = 0.15, b_ = 0.73, y_ = 0.75;
    SX xy = ca.vertcat(x, y), ab = ca.vertcat(a, b);
    Action<SX, double> test = (e, r) => {
      var f = ca.Function("f", new[] { xy, ab }, new[] { e });
      List<DM> f_out = f.call(dms(ca.vertcat(x_, y_), ca.vertcat(a_, b_)));
      checkarray((double)f_out[0][0], r, 1e-10);
    };
    double s = Math.Sin(a_ + b_), c = Math.Cos(a_ + b_), dx = x_ - a_, dy = y_ - b_;
    SX ex = ca.sin(x + y);
    test(ca.mtaylor(ex, xy, ab, 0), s);
    test(ca.mtaylor(ex, xy, ab, 1), s + (c * dx + c * dy));
    double sol2 = s + (c * dx + c * dy) - (s * dx * dx + 2 * s * dy * dx + s * dy * dy) / 2;
    test(ca.mtaylor(ex, xy, ab, 2), sol2);
    double sol3 = sol2 - (c * Math.Pow(dx, 3) + 3 * c * dy * dx * dx + 3 * c * dy * dy * dx + c * Math.Pow(dy, 3)) / 6;
    test(ca.mtaylor(ex, xy, ab, 3), sol3);
    double sol3b = (-2 * s * dx * dy - s * dx * dx) / 2 + c * dy - (c * Math.Pow(dx, 3)) / 6 + c * dx + s;
    test(ca.mtaylor(ex, xy, ab, 3, new long[] { 1, 2 }), sol3b);
    test(ca.mtaylor(ex, xy, ca.vertcat(ca.SX(0.0), ca.SX(0.0)), 4, new long[] { 1, 2 }),
         (-3 * x_ * x_ * y_ - Math.Pow(x_, 3)) / 6 + y_ + x_);
  }

  public void test_issue107() {
    message("Regression test for issue 107: +=");
    SX x = SX.sym("x");
    SX y = SX.sym("y");
    SX z = x;
    z += y;
    assertTrue(x.is_symbolic());
    assertFalse(z.is_symbolic());

    x = SX.sym("x");
    y = SX.sym("y");
    z = x;
    z = z + y;
    assertTrue(x.is_symbolic());
    assertFalse(z.is_symbolic());
  }

  public void test_evalchecking() {
    SX x = SX.sym("x", 1, 5);
    SX y = SX.sym("y", 1, 3);
    SX z = SX.sym("z", 5, 1);
    SX q = SX.sym("z", 1, 6);
    var f = ca.Function("f", new[] { x }, new[] { ca.pow(x, 2) });
    Check.Throws<Exception>(() => f.Call(y));
    Check.Throws<Exception>(() => f.Call(q));
    f.Call(z);
  }

  public void test_indexinglimits() {
    message("Limits of indexing");
    SX y = SX.sym("y", 3);
    Check.Throws<Exception>(() => y.get(ca.vertcat(0, 5)));
    Check.Throws<Exception>(() => y[new long[] { 0, 5 }] = SX.sym("a"));
    y.get(ca.vertcat(0, 2));
    y[new long[] { 0, 2 }] = SX.sym("a");
  }

  public void test_is_equal() {
    Check.Skip("known bug in Python too: not implemented");
    message("equivalent");
    SX x = SX.sym("x");
    SX a = x * x;
    SX b = x * x;
    assertTrue(ca.is_equal(a, b, 1));
  }

  public void test_SXsimplifications() {
    Check.Assume(GlobalOptions.getSimplificationOnTheFly(), "no simplification on the fly");
    message("simplifications");
    SX x = SX.sym("x");
    SX y0 = SX.sym("x");
    var ops = new List<Func<SX, SX>> {
      v => { SX y = 0.5 * v; return y + y; },
      v => { SX y = v / 2; return y + y; },
      v => { SX y = v * 0.5; return y + y; },
      v => { SX y = v * v; return -y / y * -v; },
      v => -(v * v) / (v * v) * -v,
      v => { SX y = v * v; return y / -y * -v; },
      v => { SX y = v * v; return -y / -y * v; },
      v => v - v + v,
      v => v * v - v * v + v,
      v => 4.0 * (0.25 * v),
      v => 4.0 * (v * 0.25),
      v => 4.0 * (0.25 * v),
      v => 4.0 * (v * 0.25),
      v => 0.25 * v * 4,
      v => v * 0.25 * 4,
      v => 4.0 * v / 4,
      v => 4.0 * (v / 4),
      v => v / 4 / 0.25,
      v => v * ((4.0 / v) * v / 4),
      v => v * (v * (2.0 / v) / 2),
      v => v * ((2.0 * v) / v / 2),
      v => v * (v / (2.0 * v) * 2),
      v => v + 0,
      v => 0.0 + v,
      v => v - 0,
      v => 0.0 - -v,
      v => v * 1,
      v => 1.0 * v,
      v => 1.0 * (v * 1),
      v => 1.0 * v * 1,
      v => 0.5 * v + 0.5 * v,
      v => v / 2 + v / 2,
      v => v * 0.5 + 0.5 * v,
      v => ca.SX(4.0) - ca.SX(4.0) + v,
      v => v + y0 - (y0 + v) + v,
      v => v * y0 - y0 * v + v,
      v => -v - -v + v,
    };

    foreach (Func<SX, SX> op in ops) {
      SX y = op(x);
      var f = ca.Function("f", new[] { x }, new[] { y });
      List<DM> f_out = f.call(dms(0.3));
      checkarray(f_out[0], ca.evalf(op(ca.SX(0.3))), "simplifications");
      assertEquals("x", y.ToString());

      y = op(-x);
      f = ca.Function("f", new[] { x }, new[] { y });
      f_out = f.call(dms(0.3));
      checkarray(f_out[0], ca.evalf(op(ca.SX(-0.3))), "simplifications");
      assertEquals("(-x)", y.ToString());
    }
  }

  public void test_truth() {
    message("Truth values");
    Check.Throws<Exception>(() => SX.sym("x").nonzero());
    Check.Throws<Exception>(() => ca.gt(SX.sym("x"), 0.0).nonzero());
    assertTrue(ca.SX(1.0).nonzero());
    assertFalse(ca.SX(0.0).nonzero());
    assertTrue(ca.SX(0.2).nonzero());
    assertTrue(ca.SX(-0.2).nonzero());
    Check.Throws<Exception>(() => SX.sym("x").nonzero());
    Check.Throws<Exception>(() => ca.gt(SX.sym("x"), 0.0).nonzero());
    assertTrue(ca.SX(ca.SX(1.0)).nonzero());
    assertFalse(ca.SX(ca.SX(0.0)).nonzero());
    assertTrue(ca.SX(ca.SX(0.2)).nonzero());
    assertTrue(ca.SX(ca.SX(-0.2)).nonzero());
    Check.Throws<Exception>(() => ((SX)new double[] { 2.0, 3 }).nonzero());
  }

  public void test_if_else() {
    SX x = SX.sym("x");
    SX y = ca.if_else(x, 1.0, 2.0);
    var f = ca.Function("f", new[] { x }, new[] { y });
    List<DM> f_out = f.call(dms(1.0));
    assertTrue((double)f_out[0] == 1, "if_else");
    f_out = f.call(dms(0.0));
    assertTrue((double)f_out[0] == 2, "if_else");

    double x0 = 2.1;
    y = ca.if_else(ca.gt(x, 1.0), ca.pow(x, 2), ca.pow(x, 3));
    f = ca.Function("f", new[] { x }, new[] { y });
    f_out = f.call(dms(x0));
    checkarray(f_out[0], ca.DM(Math.Pow(x0, 2)), "if_else sens");

    x0 = -2.1;
    f_out = f.call(dms(x0));
    checkarray(f_out[0], ca.DM(Math.Pow(x0, 3)), "if_else sens");
  }

  public void test_is_regular() {
    SX x = SX.sym("x");
    assertTrue(ca.SX(0.0).is_regular());
    assertFalse(ca.SX(double.PositiveInfinity).is_regular());
    Check.Throws<Exception>(() => assertTrue(x.get_nz(false, (DM)0.0).nonzero()));
    assertTrue(ca.SX(ca.vertcat(0, 1)).is_regular());
    assertFalse(ca.SX(ca.vertcat(0, double.PositiveInfinity)).is_regular());
    assertFalse(ca.vertcat(x, ca.SX(double.PositiveInfinity)).is_regular());
    Check.Throws<Exception>(() => assertFalse(ca.vertcat(x, x).is_regular()));
  }

  public void test_symvar() {
    SX a = SX.sym("a");
    SX b = SX.sym("b");
    SX c = SX.sym("c");
    SX e = ca.cos(a * b) + c;
    List<SX> w = ca.symvar(e);
    assertEquals(3, w.Count);
    if (GlobalOptions.getSimplificationOnTheFly()) {
      assertTrue(ca.is_equal(w[0], a));
      assertTrue(ca.is_equal(w[1], b));
      assertTrue(ca.is_equal(w[2], c));
    }
  }

  public void test_poly_coeff() {
    SX x = SX.sym("x");
    SX a = SX.sym("a");
    SX c = SX.sym("c");
    SX p = ca.poly_coeff(12.0 * ca.pow(x, 4) + ca.pow(x, 2) + a * x + c, x);
    assertTrue(ca.is_equal(p[0], 12.0));
    assertTrue(ca.is_equal(p[1], 0.0));
    assertTrue(ca.is_equal(p[2], 1.0));
    assertTrue(ca.is_equal(p[3], a));
    assertTrue(ca.is_equal(p[4], c));

    p = ca.poly_coeff((x - a) * (x + a), x);
    assertTrue(ca.is_equal(p[0], 1.0));
    assertTrue(ca.is_equal(p[1], 0.0));
  }

  public void test_poly_roots() {
    SX p = SX.sym("[a,b]");
    SX r = ca.poly_roots(p);
    var f = ca.Function("f", new[] { p }, new[] { r });
    DM f_in = ca.vertcat(2, 7);
    DM a_ = f_in[0], b_ = f_in[1];
    List<DM> f_out = f.call(new[] { f_in });
    checkarray(f_out[0], ca.vertcat(-b_ / a_));

    p = SX.sym("[a,b]");
    r = ca.poly_roots(ca.vertcat(p, ca.SX(0.0)));
    f = ca.Function("f", new[] { p }, new[] { r });
    f_out = f.call(new[] { f_in });
    checkarray(f_out[0], ca.vertcat(-b_ / a_, ca.DM(0.0)));

    p = SX.sym("[a,b,c]");
    r = ca.poly_roots(p);
    f = ca.Function("f", new[] { p }, new[] { r });
    f_in = ca.vertcat(1.13, 7, 3);
    a_ = f_in[0];
    b_ = f_in[1];
    DM c_ = f_in[2];
    DM d = ca.pow(b_, 2) - 4.0 * a_ * c_;
    f_out = f.call(new[] { f_in });
    DM x0 = (-b_ - ca.sqrt(d)) / 2 / a_;
    DM x1 = (-b_ + ca.sqrt(d)) / 2 / a_;
    checkarray(f_out[0], ca.vertcat(x0, x1));

    p = SX.sym("[a,b,c,d]");
    r = ca.poly_roots(p);
    f = ca.Function("f", new[] { p }, new[] { r });
    f_out = f.call(new[] { ca.vertcat(11, 1.3, -1.7, 0.1) });
    checkarray(f_out[0], ca.vertcat(0.298028, -0.479787, 0.0635774), 1e-5, "");

    p = SX.sym("[a,b,c,d,e]");
    r = ca.poly_roots(p);
    f = ca.Function("f", new[] { p }, new[] { r });
    f_out = f.call(new[] { ca.vertcat(3, 6, -123, -126, 1080) });
    checkarray(f_out[0], ca.vertcat(5, 3, -4, -6), 1e-5, "");
  }

  public void test_eig_symbolic() {
    SX x = SX.sym("x", 2, 2);
    var f = ca.Function("f", new[] { x }, new[] { ca.eig_symbolic(x) });
    List<DM> f_out = f.call(new[] { ca.DM(new double[,] { { 2, 0.1 }, { 0.3, 0.7 } }) });
    checkarray(f_out[0], ca.vertcat(0.67732, 2.02268), 1e-5, "");

    x = SX.sym("x", 2);
    f = ca.Function("f", new[] { x }, new[] { ca.eig_symbolic(ca.diag(x)) });
    DM f_in = ca.vertcat(3, 7);
    f_out = f.call(new[] { f_in });
    checkarray(f_out[0], f_in);

    x = SX.sym("x", 5);
    f = ca.Function("f", new[] { x }, new[] { ca.eig_symbolic(ca.diag(x)) });
    f_in = ca.vertcat(3, 7, 2, 1, 6);
    f_out = f.call(new[] { f_in });
    checkarray(f_out[0], f_in);

    x = SX.sym("x", 2, 2);
    SX y = SX.sym("y", 2);
    f = ca.Function("f", new[] { x, y }, new[] { ca.eig_symbolic(ca.diagcat(x, ca.diag(y))) });
    f_out = f.call(dms(ca.DM(new double[,] { { 2, 0.1 }, { 0.3, 0.7 } }), ca.vertcat(3, 7)));
    checkarray(f_out[0], ca.vertcat(0.67732, 2.02268, 3, 7), 1e-5, "");

    x = SX.sym("x", 3, 3);
    x[2, 0] = 0.0;
    x[1, 0] = 0.0;
    x = ca.sparsify(x);
    SX e = ca.eig_symbolic(x);
    f = ca.Function("f", new[] { x }, new[] { e });
    f_in = ca.DM(f.sparsity_in(0), ca.vertcat(1, 2, 3, 4, 5, 6, 7));
    f_in.print_dense();
    f_out = f.call(new[] { f_in });
    checkarray(f_out[0], ca.vertcat(1, -0.29150, 10.29150), 1e-5, "");

    x = SX.sym("x", 3, 3);
    x[2, 0] = 0.0;
    x[1, 0] = 0.0;
    x[2, 1] = 0.0;
    x = ca.sparsify(x);
    e = ca.eig_symbolic(x);
    f = ca.Function("f", new[] { x }, new[] { e });
    f_in = ca.DM(f.sparsity_in(0), ca.vertcat(1, 2, 3, 4, 5, 6));
    f_in.print_dense();
    f_out = f.call(new[] { f_in });
    checkarray(f_out[0], ca.vertcat(1, 3, 6), 1e-5, "");

    x = SX.sym("x", Sparsity.upper(5));
    f = ca.Function("f", new[] { x }, new[] { ca.eig_symbolic(x) });
    var fin = ca.DM(x.sparsity(), 0.0);
    fin.set(ca.diag(range(5)), false, Sparsity.diag(5));
    checkarray(f.Call(fin), range(5));
  }

  public void test_jacobian_empty() {
    SX x = SX.sym("x", 3);
    SX s = ca.jacobian(ca.DM(0, 0), x);
    assertEquals(0L, s.size1());
    assertEquals(3L, s.size2());

    s = ca.jacobian(x, SX.sym("x", 0, 4));
    assertEquals(3L, s.size1());
    assertEquals(0L, s.size2());
  }

  public void test_empty_SX() {
    var s = ca.SX((DM)new double[0]);
    assertEquals(0L, s.size1());
    assertEquals(1L, s.size2());
    ca.vertcat(SX.sym("x"), ca.SX((DM)new double[0]));
    ca.SX(ca.Sparsity(3, 3), (DM)new double[0]);
    ca.SX(ca.Sparsity(3, 3), (DM)new double[] { 5 });
    assertInException("fully sparse", () => ca.SX(ca.Sparsity(3, 3), (DM)new double[] { 5, 4 }));
  }

  public void test_mul_sparsity() {
    int N = 10;
    SX x = SX.sym("x", N, N);
    SX y = SX.sym("y", N, N);
    DM x_ = randDM(N, N);
    DM y_ = randDM(N, N);
    Sparsity filt = Sparsity.diag(N).plus(Sparsity.triplet(N, N, new long[] { 1 }, new long[] { 3 }));
    var f = ca.Function("f", new[] { x, y }, new[] { ca.mtimes(x, y) });
    var g = ca.Function("g", new[] { x, y }, new[] { ca.mac(x, y, SX.zeros(filt)) });
    List<DM> f_out = f.call(dms(x_, y_));
    List<DM> g_out = g.call(dms(x_, y_));
    checkarray(DM.ones(filt), DM.ones(g.sparsity_out(0)));
    checkarray(f_out[0].get(filt), g_out[0]);
  }

  /// <summary>A file of test/data: from $CASADI_TEST_DATA, else ../data as in the Java suite</summary>
  static string dataFile(string name) {
    foreach (string d in new[] { Environment.GetEnvironmentVariable("CASADI_TEST_DATA"), "../data", "../../test/data" }) {
      if (d == null) continue;
      string p = Path.Combine(d, name);
      if (File.Exists(p)) return p;
    }
    throw new SkipTest("test data " + name + " not found (set CASADI_TEST_DATA to casadi's test/data)");
  }

  public void test_large_hessian() {
    Sparsity A = Sparsity.from_file(dataFile("apoa1-2.mtx"));
    var H = ca.DM(A, range((int)A.nnz()));
    H = H + H.T();
    H = H[ca.Slice(0, 20000), ca.Slice(0, 20000)];
    SX x = SX.sym("x", H.size1());
    var f = ca.Function("f", new[] { x }, new[] { ca.mtimes(new[] { x.T(), ca.SX(H), x }) }, D(("verbose", true)));
    H = H * 2;
    Function h = hessianOld(f, 0, 0);
    List<DM> h_out = h.call(dms(0.0));
    assertTrue(h.sparsity_out(0).is_equal(H.sparsity()));
    checkarray((DM)h_out[0].nonzeros().ToArray(), (DM)H.nonzeros().ToArray());
  }

  public void test_mxnulloutput() {
    var a = ca.SX(5, 0);
    SX b = SX.sym("x", 2);
    MX bm = MX.sym("x", 2);

    var f = ca.Function("f", new[] { b }, new[] { a });
    MX c = f.Call(bm);
    assertEquals(5L, c.size1());
    assertEquals(0L, c.size2());
    SX cs = f.Call(b);
    assertEquals(5L, cs.size1());
    assertEquals(0L, cs.size2());

    a = ca.SX(0, 0);
    f = ca.Function("f", new[] { b }, new[] { a });
    c = f.Call(bm);
    assertEquals(0L, c.size1());
    assertEquals(0L, c.size2());
    cs = f.Call(b);
    assertEquals(0L, cs.size1());
    assertEquals(0L, cs.size2());
  }

  public void test_mxnull() {
    var a = ca.SX(5, 0);
    var b = ca.SX(0, 3);
    SX c = ca.mtimes(a, b);
    assertEquals(0L, c.nnz());

    a = ca.SX(5, 3);
    b = ca.SX(3, 4);
    c = ca.mtimes(a, b);
    assertEquals(0L, c.nnz());
  }

  public void test_mxnullop() {
    var c = ca.SX(0, 0);
    SX x = SX.sym("x", 2, 3);
    Check.Throws<Exception>(() => _ = x + c);
    Check.Throws<Exception>(() => _ = x / c);
  }

  public void test_copysign_derivatives() {
    for (int t = 0; t < 2; t++) {
      Function f;
      if (t == 0) {
        SX x = SX.sym("x"), y = SX.sym("y"), z = ca.copysign(x, y);
        f = ca.Function("f", new[] { x, y }, new[] { ca.jacobian(z, x), ca.gradient(z, x), ca.jtimes(z, x, 1.0), ca.jtimes(z, x, 1.0, true), ca.jacobian(z, y) });
      } else {
        MX x = MX.sym("x"), y = MX.sym("y"), z = ca.copysign(x, y);
        f = ca.Function("f", new[] { x, y }, new[] { ca.jacobian(z, x), ca.gradient(z, x), ca.jtimes(z, x, 1.0), ca.jtimes(z, x, 1.0, true), ca.jacobian(z, y) });
      }
      foreach (double x0 in new double[] { -2, 0, 2 })
        foreach (double y0 in new double[] { -0.5, -0.0, 0.0, 0.5 }) {
          double d = Math.Sign(x0) * Math.CopySign(1, y0);
          List<DM> r = f.call(new DM[] { x0, y0 });
          checkarray(ca.vertcat(r), ca.vertcat(d, d, d, d, 0), (t == 0 ? "SX" : "MX") + " x=" + x0 + " y=" + y0);
        }
    }
  }

  public void test_copysign() {
    SX x = SX.sym("x");
    SX y = SX.sym("y");
    SX z = ca.copysign(x, y);
    var f = ca.Function("f", new[] { x, y }, new[] { z });
    double[,] inp = { { 2, 0.5 }, { 2, -0.5 }, { -2, 0.5 }, { -2, -0.5 }, { 2, 0 } };
    double[] rf = { 2, -2, 2, -2, 2 };
    for (int i = 0; i < inp.GetLength(0); i++) checkarray(f.call(dms(inp[i, 0], inp[i, 1]))[0], ca.vertcat(rf[i]));

    Function J = jacobianOld(f, 0, 0);
    double[] jref = { 1, -1, -1, 1, 1 };
    for (int i = 0; i < inp.GetLength(0); i++) checkarray(J.call(dms(inp[i, 0], inp[i, 1]))[0], ca.vertcat(jref[i]));

    J = jacobianOld(f, 1, 0);
    for (int i = 0; i < inp.GetLength(0); i++) checkarray(J.call(dms(inp[i, 0], inp[i, 1]))[0], ca.vertcat(0));
  }

  public void test_depends_on() {
    SX a = SX.sym("a");
    SX b = SX.sym("b");
    SX ab = ca.vertcat(a, b);
    var zero = ca.SX(0.0);
    assertTrue(ca.depends_on(ca.pow(a, 2), a));
    assertTrue(ca.depends_on(a, a));
    assertFalse(ca.depends_on(0.0, a));
    assertTrue(ca.depends_on(ca.pow(a, 2), ab));
    assertTrue(ca.depends_on(a, ab));
    assertFalse(ca.depends_on(0.0, ab));
    assertTrue(ca.depends_on(ca.pow(b, 2), ab));
    assertTrue(ca.depends_on(b, ab));
    assertTrue(ca.depends_on(ca.pow(a, 2) + ca.pow(b, 2), ab));
    assertTrue(ca.depends_on(a + b, ab));
    assertTrue(ca.depends_on(ca.vertcat(zero, a), a));
    assertTrue(ca.depends_on(ca.vertcat(a, zero), a));
    assertTrue(ca.depends_on(ca.vertcat(ca.pow(a, 2), ca.pow(b, 2)), ab));
    assertTrue(ca.depends_on(ca.vertcat(a, zero), ab));
    assertTrue(ca.depends_on(ca.vertcat(zero, b), ab));
    assertTrue(ca.depends_on(ca.vertcat(b, zero), ab));
    assertFalse(ca.depends_on(ca.vertcat(zero, zero), ab));
  }

  public void test_which_depends() {
    {
      SX x = SX.sym("x"), y = SX.sym("y"), p = SX.sym("p"), z = SX.sym("z");
      SX e = ca.vertcat(ca.SX(0.0), x, y, p, 2.0 * ca.pow(p, 3), x * y, x * p, ca.sin(x), ca.cos(y), ca.sqrt(x + y), p * p * x, x * y * p);
      assertBools(ca.which_depends(e, ca.vertcat(x, y), 2, true), 0, 0, 0, 0, 0, 1, 0, 1, 1, 1, 0, 1);
      assertBools(ca.which_depends(e, ca.vertcat(x, y), 1, true), 0, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1);
      e = ca.vertcat(x * p, x + y);
      assertBools(ca.which_depends(e, ca.vertcat(x, y, p, z), 2, false), 1, 0, 1, 0);
      assertBools(ca.which_depends(e, ca.vertcat(x, y, p, z), 1, false), 1, 1, 1, 0);
      e = ca.vertcat(x * p, x + z * y);
      assertBools(ca.which_depends(e, ca.vertcat(x, y, p), 2, false), 1, 0, 1);
      assertBools(ca.which_depends(e, ca.vertcat(x, y, p), 1, false), 1, 1, 1);
      assertBools(ca.which_depends(e, ca.vertcat(x, y, p, z), 2, false), 1, 1, 1, 1);
      assertBools(ca.which_depends(e, ca.vertcat(x, y, p, z), 1, false), 1, 1, 1, 1);
      e = ca.vertcat(ca.sin(x + y) + p);
      assertBools(ca.which_depends(e, ca.vertcat(x, y, p, z), 2, false), 1, 1, 0, 0);
      assertBools(ca.which_depends(e, ca.vertcat(x, y, p, z), 1, false), 1, 1, 1, 0);
      e = ca.vertcat(ca.sin(x) * ca.pow(p, 2), ca.pow(y, 2));
      assertBools(ca.which_depends(e, ca.vertcat(x, y, p), 2, true), 1, 1);
      assertBools(ca.which_depends(e, ca.vertcat(x, y, p), 2, false), 1, 1, 1);
      e = ca.vertcat(ca.pow(x, 2) * p, y);
      assertBools(ca.which_depends(e, ca.vertcat(x, y, p), 2, true), 1, 0);
      assertBools(ca.which_depends(e, ca.vertcat(x, y, p), 2, false), 1, 0, 1);
    }
    {
      MX x = MX.sym("x"), y = MX.sym("y"), p = MX.sym("p"), z = MX.sym("z");
      MX e = ca.vertcat(ca.MX(0.0), x, y, p, 2.0 * ca.pow(p, 3), x * y, x * p, ca.sin(x), ca.cos(y), ca.sqrt(x + y), p * p * x, x * y * p);
      assertBools(ca.which_depends(e, ca.vertcat(x, y), 2, true), 0, 0, 0, 0, 0, 1, 0, 1, 1, 1, 0, 1);
      assertBools(ca.which_depends(e, ca.vertcat(x, y), 1, true), 0, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1);
      e = ca.vertcat(x * p, x + y);
      assertBools(ca.which_depends(e, ca.vertcat(x, y, p, z), 2, false), 1, 0, 1, 0);
      assertBools(ca.which_depends(e, ca.vertcat(x, y, p, z), 1, false), 1, 1, 1, 0);
      e = ca.vertcat(x * p, x + z * y);
      assertBools(ca.which_depends(e, ca.vertcat(x, y, p), 2, false), 1, 0, 1);
      assertBools(ca.which_depends(e, ca.vertcat(x, y, p), 1, false), 1, 1, 1);
      assertBools(ca.which_depends(e, ca.vertcat(x, y, p, z), 2, false), 1, 1, 1, 1);
      assertBools(ca.which_depends(e, ca.vertcat(x, y, p, z), 1, false), 1, 1, 1, 1);
      e = ca.vertcat(ca.sin(x + y) + p);
      assertBools(ca.which_depends(e, ca.vertcat(x, y, p, z), 2, false), 1, 1, 0, 0);
      assertBools(ca.which_depends(e, ca.vertcat(x, y, p, z), 1, false), 1, 1, 1, 0);
      e = ca.vertcat(ca.sin(x) * ca.pow(p, 2), ca.pow(y, 2));
      assertBools(ca.which_depends(e, ca.vertcat(x, y, p), 2, true), 1, 1);
      assertBools(ca.which_depends(e, ca.vertcat(x, y, p), 2, false), 1, 1, 1);
      e = ca.vertcat(ca.pow(x, 2) * p, y);
      assertBools(ca.which_depends(e, ca.vertcat(x, y, p), 2, true), 1, 0);
      assertBools(ca.which_depends(e, ca.vertcat(x, y, p), 2, false), 1, 0, 1);
    }
  }

  public void test_if_else_zero_sens() {
    SX x = SX.sym("x");
    SX a = 1.0 + 3.0 * x + ca.sqrt(3.0 * x) * x + 7.0 * x;
    SX b = 1.0 + 2.0 * x + ca.sin(2.0 * x) * x + x;
    SX z = ca.if_else(ca.gt(x, 0.0), a, b) * x;

    var f = ca.Function("f", new[] { x }, new[] { z, ca.jacobian(z, x) });
    var fa = ca.Function("f", new[] { x }, new[] { a * x, ca.jacobian(a * x, x) });
    var fb = ca.Function("f", new[] { x }, new[] { b * x, ca.jacobian(b * x, x) });

    List<DM> r = f.call(new DM[] { 3.0 }), ra = fa.call(new DM[] { 3.0 });
    for (int i = 0; i < r.Count; i++) checkarray(r[i], ra[i]);
    r = f.call(new DM[] { -3.0 });
    List<DM> rb = fb.call(new DM[] { -3.0 });
    for (int i = 0; i < r.Count; i++) checkarray(r[i], rb[i]);

    f = ca.Function("f", new[] { x }, new[] { z });
    fa = ca.Function("f", new[] { x }, new[] { a * x });
    fb = ca.Function("f", new[] { x }, new[] { b * x });
    checkfunction(f, fa, dms(3.0));
    checkfunction(f, fb, dms(-3.0));
  }

  public void test_pw_const() {
    SX t = SX.sym("t");
    SX e = ca.pw_const(t, ca.vertcat(0, 2, 3), ca.vertcat(7, 1, 3, 5));
    var E = ca.Function("E", new[] { t }, new[] { e });
    double[,] cases = { { -2, 7 }, { -1, 7 }, { 0, 1 }, { 1, 1 }, { 1.9999, 1 }, { 2, 3 }, { 2.5, 3 }, { 3, 5 }, { 10, 5 } };
    for (int i = 0; i < cases.GetLength(0); i++) checkarray(E.Call((DM)cases[i, 0]), cases[i, 1]);
  }

  public void test_pw_lin() {
    SX t = SX.sym("t");
    SX e = ca.pw_lin(t, ca.vertcat(0, 2, 3, 5), ca.vertcat(7, 1, 3, 2));
    var E = ca.Function("E", new[] { t }, new[] { e });
    double[,] cases = { { -2, 13 }, { -1, 10 }, { 0, 7 }, { 1, 4 }, { 2, 1 }, { 2.5, 2 }, { 3, 3 }, { 4, 2.5 }, { 5, 2 }, { 7, 1 } };
    for (int i = 0; i < cases.GetLength(0); i++) checkarray(E.Call((DM)cases[i, 0]), cases[i, 1]);
  }

  public void test_quadratic() {
    {
      SX x = SX.sym("x"), p = SX.sym("p"), y = SX.sym("y");
      assertFalse(ca.is_quadratic(ca.sin(x), x));
      assertFalse(ca.is_quadratic(ca.pow(x, 3), x));
      assertTrue(ca.is_quadratic(ca.pow(x, 2), x));
      assertTrue(ca.is_quadratic(4.0 * x, x));
      assertTrue(ca.is_quadratic(5.0, x));

      assertFalse(ca.is_quadratic(ca.sin(x) * ca.pow(p, 4), x));
      assertFalse(ca.is_quadratic(ca.pow(x, 3) * ca.pow(p, 4), x));
      assertTrue(ca.is_quadratic(ca.pow(x, 2) * ca.pow(p, 4), x));
      assertTrue(ca.is_quadratic(x * ca.pow(p, 4), x));
      assertTrue(ca.is_quadratic(5.0 * ca.pow(p, 4), x));

      assertFalse(ca.is_linear(ca.sin(x), x));
      assertFalse(ca.is_linear(ca.pow(x, 3), x));
      assertFalse(ca.is_linear(ca.pow(x, 2), x));
      assertTrue(ca.is_linear(3.0 * x, x));
      assertTrue(ca.is_linear(5.0, x));

      assertFalse(ca.is_linear(ca.sin(x) * ca.pow(p, 4), x));
      assertFalse(ca.is_linear(ca.pow(x, 3) * ca.pow(p, 4), x));
      assertFalse(ca.is_linear(ca.pow(x, 2) * ca.pow(p, 4), x));
      assertTrue(ca.is_linear(x * ca.pow(p, 4), x));
      assertTrue(ca.is_linear(5.0 * ca.pow(p, 4), x));

      SX xy = ca.vertcat(x, y);
      SX bse = ca.pow(x, 2) + 3.0 * ca.pow(y, 2) + 0.5 * x * y;
      SX z = bse + 7.0 * x + 6.0 * y + 7;
      var ABc = ca.quadratic_coeff(z, xy);
      assertInException("non-quadratic", () => ca.quadratic_coeff(z + ca.sin(x), xy));
      assertInException("scalar", () => ca.quadratic_coeff(xy, x));

      SX z2 = bse - p * y + 7.0 * x + 6.0 * y + 7;
      var (A, b, c) = ca.quadratic_coeff(z2, xy);
      SX e = 0.5 * ca.bilin(A, xy, xy) + ca.dot(b, xy) + c;
      var f = ca.Function("f", new[] { xy, p }, new[] { z2 });
      var f2 = ca.Function("f", new[] { xy, p }, new[] { e });
      checkfunction(f, f2, dms(1.1, 1.3));

      assertInException("non-linear", () => ca.linear_coeff(z, xy));
      assertInException("vector", () => ca.linear_coeff(ca.blockcat(new[] { new[] { x, y }, new[] { y, x } }), x));

      SX z3 = ca.vertcat(7.0 * x + 6.0 * y + 7, 5.0 - p * y);
      var (A2, b2) = ca.linear_coeff(z3, xy);
      e = ca.mtimes(A2, xy) + b2;
      f = ca.Function("f", new[] { xy, p }, new[] { z3 });
      f2 = ca.Function("f", new[] { xy, p }, new[] { e });
      checkfunction(f, f2, dms(1.1, 1.3));
    }
    {
      MX x = MX.sym("x"), p = MX.sym("p"), y = MX.sym("y");
      assertFalse(ca.is_quadratic(ca.sin(x), x));
      assertFalse(ca.is_quadratic(ca.pow(x, 3), x));
      assertTrue(ca.is_quadratic(ca.pow(x, 2), x));
      assertTrue(ca.is_quadratic(4.0 * x, x));
      assertTrue(ca.is_quadratic(5.0, x));

      assertFalse(ca.is_quadratic(ca.sin(x) * ca.pow(p, 4), x));
      assertFalse(ca.is_quadratic(ca.pow(x, 3) * ca.pow(p, 4), x));
      assertTrue(ca.is_quadratic(ca.pow(x, 2) * ca.pow(p, 4), x));
      assertTrue(ca.is_quadratic(x * ca.pow(p, 4), x));
      assertTrue(ca.is_quadratic(5.0 * ca.pow(p, 4), x));

      assertFalse(ca.is_linear(ca.sin(x), x));
      assertFalse(ca.is_linear(ca.pow(x, 3), x));
      assertFalse(ca.is_linear(ca.pow(x, 2), x));
      assertTrue(ca.is_linear(3.0 * x, x));
      assertTrue(ca.is_linear(5.0, x));

      assertFalse(ca.is_linear(ca.sin(x) * ca.pow(p, 4), x));
      assertFalse(ca.is_linear(ca.pow(x, 3) * ca.pow(p, 4), x));
      assertFalse(ca.is_linear(ca.pow(x, 2) * ca.pow(p, 4), x));
      assertTrue(ca.is_linear(x * ca.pow(p, 4), x));
      assertTrue(ca.is_linear(5.0 * ca.pow(p, 4), x));

      MX xy = ca.vertcat(x, y);
      MX bse = ca.pow(x, 2) + 3.0 * ca.pow(y, 2) + 0.5 * x * y;
      MX z = bse + 7.0 * x + 6.0 * y + 7;
      var ABc = ca.quadratic_coeff(z, xy);
      assertInException("non-quadratic", () => ca.quadratic_coeff(z + ca.sin(x), xy));
      assertInException("scalar", () => ca.quadratic_coeff(xy, x));

      MX z2 = bse - p * y + 7.0 * x + 6.0 * y + 7;
      var (A, b, c) = ca.quadratic_coeff(z2, xy);
      MX e = 0.5 * ca.bilin(A, xy, xy) + ca.dot(b, xy) + c;
      var f = ca.Function("f", new[] { xy, p }, new[] { z2 });
      var f2 = ca.Function("f", new[] { xy, p }, new[] { e });
      checkfunction(f, f2, dms(1.1, 1.3));

      assertInException("non-linear", () => ca.linear_coeff(z, xy));
      assertInException("vector", () => ca.linear_coeff(ca.blockcat(new[] { new[] { x, y }, new[] { y, x } }), x));

      MX z3 = ca.vertcat(7.0 * x + 6.0 * y + 7, 5.0 - p * y);
      var (A2, b2) = ca.linear_coeff(z3, xy);
      e = ca.mtimes(A2, xy) + b2;
      f = ca.Function("f", new[] { xy, p }, new[] { z3 });
      f2 = ca.Function("f", new[] { xy, p }, new[] { e });
      checkfunction(f, f2, dms(1.1, 1.3));
    }
  }

  public void test_evalf() {
    SX x = SX.sym("x");
    var y = ca.SX(5.0);
    checkarray(ca.evalf(y), 5.0);
    assertInException("since variables [x] are free", () => ca.evalf(x));
  }

  public void test_output_sx() {
    MX x = MX.sym("x");
    MX y = MX.sym("y");
    MX z = MX.sym("z");
    MX w = ca.mtimes(ca.sqrt(z), ca.sin(x * y));
    var f = ca.Function("f", new[] { x, y, z }, new[] { x * y * z, w });
    Console.WriteLine(f);

    SX xs = SX.sym("x");
    SX ys = SX.sym("y");
    SX zs = SX.sym("z");
    var args = new List<SX> { xs / ys, 10.0 * zs, ys * zs };
    List<SX> v = f.call(args, false, true);
    SX output_0 = v[0];
    SX output_1 = v[1];
    assertEquals(output_0.dep(0).element_hash(), output_1.dep(0).element_hash());
    SX call_node = output_0.dep(0);
    Console.WriteLine("here");
    SX output_1b = call_node.get_output(1);
    assertEquals(output_1.element_hash(), output_1b.element_hash());
    output_1b = call_node.get_output(1);
    assertTrue(output_1b is not null);
  }

  public void test_call_fun() {
    DM A = ca.sparsify(ca.DM(new double[,] { { 1, 0, 1 }, { 0, 0, 6 }, { 0, 8, 9 } }));
    MX xm = MX.sym("x", 3);
    MX ym = MX.sym("y");
    MX zm = MX.sym("z");
    MX wm = ca.mtimes(A * ca.sqrt(zm), ca.sin(xm * ym));
    var f = ca.Function("f", new[] { xm, ym, zm }, new[] { (xm * ym * zm)[ca.Slice(0, 2)], wm });
    Console.WriteLine(f);

    SX x = SX.sym("x", 3);
    SX y = SX.sym("y");
    SX z = SX.sym("z");
    var args = new List<SX> { x / y, 10.0 * z, y * z };
    List<SX> v = f.call(args);
    List<SX> v2 = f.call(args, false, true);
    Console.WriteLine(string.Join(", ", v2));
    Console.WriteLine(ca.vertcat(v2));

    assertTrue(v2[0][0].is_output());
    assertFalse(y.is_output());
    long k = 0;
    for (int i = 0; i < f.n_out(); i++)
      for (int j = 0; j < f.nnz_out(i); j++) {
        assertEquals(k, v2[i][j].which_output());
        k += 1;
      }
    SX callnode = v2[0][0].dep(0);
    assertTrue(callnode.is_call());
    assertFalse(y.is_call());
    assertTrue(callnode.has_output());
    assertFalse(y.has_output());
    assertTrue(callnode.which_function().hash() == f.hash());
    assertTrue(callnode.n_dep() == 5);
    callnode.dep(4);
    SX e = callnode.get_output(3);
    Console.WriteLine(e + " " + e.element_hash());
    SX ee = callnode.get_output(3);
    Console.WriteLine(ee + " " + e.element_hash());
    assertTrue(callnode.get_output(3).ToString().Contains("{3}"));

    Console.WriteLine(ca.cos(ca.mtimes(v2[1], v2[0].T()) / y).dim());
    var F2 = ca.Function("F", new[] { x, y, z }, new[] { ca.sin(v2[0] * y - args[2]), ca.cos(ca.mtimes(v2[1], v2[0].T()) / y) });
    Function F2f = F2;
    assertOutput("f:(i0[3],i1,i2)->(o0[2],o1[3])", () => Console.WriteLine(f));
    assertOutput("[[@2,@3],[@4,@5,@7]] = f([@2,@3,@4],@5,@6);", () => F2f.disp(true));

    List<DM> inputs = dms(ca.vertcat(0.1, 1.7, 2.3), 1.13, 0.11);
    var F1 = ca.Function("F", new[] { x, y, z }, new[] { ca.sin(v[0] * y - args[2]), ca.cos(ca.mtimes(v[1], v[0].T()) / y) });
    foreach (bool ad_weight in new[] { true, false })
      foreach (bool ad_weight_sp in new[] { true, false }) {
        F2 = ca.Function("F", new[] { x, y, z }, new[] { ca.sin(v2[0] * y - args[2]), ca.cos(ca.mtimes(v2[1], v2[0].T()) / y) },
                          D(("ad_weight_sp", ad_weight_sp), ("ad_weight", ad_weight)));
        checkfunction(F1, F2, inputs);
        checkSerialize(F2, inputs);
        foreach (bool avoid_stack in new[] { true, false }) checkCodegen(F2, inputs, D(("avoid_stack", avoid_stack)));
      }

    SX instr = F2.instructions_sx();
    bool check = false;
    for (long kk = 0; kk < F2.n_instructions(); kk++) {
      long op = F2.instruction_id(kk);
      if (op == ca.OP_CALL) {
        check = true;
        assertTrue(instr[kk].is_call());
        assertTrue(new long[] { 2, 3, 4, 5, 7 }.SequenceEqual(F2.instruction_output(kk)), "instruction_output");
        assertTrue(new long[] { 2, 3, 4, 5, 6 }.SequenceEqual(F2.instruction_input(kk)), "instruction_input");
      }
    }
    assertTrue(check);

    F2 = ca.Function("F", new[] { x, y, z }, new[] { ca.sin(v2[0] * y - args[2]), ca.cos(ca.mtimes(v2[1], v2[0].T()) / y) }, D(("cse", true)));
    checkfunction(F1, F2, inputs);
    List<Function> res = F2.find_functions();
    assertTrue(res.Count == 1);
    assertTrue(res[0].hash() == f.hash());

    F2 = ca.Function("F", new[] { x, y, z }, F2.call(new[] { x, y, z }));
    checkfunction(F1, F2, inputs);

    F1 = ca.Function("F", new[] { x, y, z }, F1.call(new[] { 2.0 * x, 3.0 * y, 4.0 * z }));
    F2 = ca.Function("F", new[] { x, y, z }, F2.call(new[] { 2.0 * x, 3.0 * y, 4.0 * z }));
    checkfunction(F1, F2, inputs);

    F1 = ca.Function("F", new[] { x, y, z }, F1.call(new[] { x, 3.0 * y, 4.0 * z }));
    F2 = ca.Function("F", new[] { x, y, z }, F2.call(new[] { x, 3.0 * y, 4.0 * z }));
    checkfunction(F1, F2, inputs);

    var F3 = ca.Function("F", new[] { x, y, z }, new[] { ca.sin(v2[0] * y - args[2]), ca.cos(ca.mtimes(v2[1], v2[0].T())[1] / y) });
    assertOutput("[[@2,@3],[NULL,@4,NULL]] = f([@2,@3,@4],@5,@6);", () => F3.disp(true));

    F1 = ca.Function("F", new[] { x, y, z }, new[] { ca.sin(v[0] * y - args[2]), ca.cos(ca.mtimes(v[1], v[0].T())[1] / y) });
    foreach (bool ad_weight in new[] { true, false })
      foreach (bool ad_weight_sp in new[] { true, false }) {
        F2 = ca.Function("F", new[] { x, y, z }, new[] { ca.sin(v2[0] * y - args[2]), ca.cos(ca.mtimes(v2[1], v2[0].T())[1] / y) },
                          D(("ad_weight_sp", ad_weight_sp), ("ad_weight", ad_weight)));
        checkfunction(F1, F2, inputs);
        checkSerialize(F2, inputs);
        foreach (bool avoid_stack in new[] { true, false }) checkCodegen(F2, inputs, D(("avoid_stack", avoid_stack)));
      }
  }

  public void test_call_fun2() {
    Console.WriteLine("test_call_fun2");
    for (int t = 0; t < 2; t++)
      foreach (bool never_inline in new[] { true, false }) {
        Console.WriteLine((t == 0 ? "SX" : "MX") + " " + never_inline);
        DM A = ca.sparsify(ca.DM(new double[,] { { 1, 0, 1 }, { 0, 0, 6 }, { 0, 8, 9 } }));
        Function f;
        if (t == 0) {
          SX x0 = SX.sym("x", 3), y0 = SX.sym("y"), z0 = SX.sym("z");
          SX w = ca.mtimes(A * ca.sqrt(z0), ca.sin(x0 * y0));
          f = ca.Function("fun", new[] { x0, y0, z0 }, new[] { (x0 * y0 * z0)[ca.Slice(0, 2)], w }, D(("never_inline", never_inline)));
        } else {
          MX x0 = MX.sym("x", 3), y0 = MX.sym("y"), z0 = MX.sym("z");
          MX w = ca.mtimes(A * ca.sqrt(z0), ca.sin(x0 * y0));
          f = ca.Function("fun", new[] { x0, y0, z0 }, new[] { (x0 * y0 * z0)[ca.Slice(0, 2)], w }, D(("never_inline", never_inline)));
        }
        SX x = SX.sym("x", 3);
        SX y = SX.sym("y");
        SX z = SX.sym("z");
        var args = new List<SX> { x / y, 10.0 * z, y * z };
        List<SX> res = f.call(args);
        string s = "[" + string.Join(", ", res) + "]";
        Console.WriteLine(s);
        assertEquals(never_inline, s.Contains("fun"));
      }
  }

  public void test_call_fun_nominal_out_deriv() {
    Check.Assume(hasRootfinder("newton"), "no rootfinder newton");
    SX x = SX.sym("x", 2);
    SX p = SX.sym("p", 2);
    Function rf = ca.rootfinder("rf", "newton", new Dictionary<string, SX> {
      ["x"] = x, ["g"] = ca.vertcat(ca.sin(x[0]) - p[0], ca.sin(x[0] + x[1]) - p[1] * p[0]), ["p"] = p });
    var f = ca.Function("f", new[] { p }, new[] { ca.exp(rf.call(new Dictionary<string, SX> { ["p"] = ca.pow(p, 2) })["x"]) });
    MX P = MX.sym("p", 2);
    var fref = ca.Function("f", new[] { P }, new[] { ca.exp(rf.call(new Dictionary<string, MX> { ["p"] = ca.pow(P, 2) })["x"]) });
    checkfunction(f, fref, dms(ca.vertcat(0.1, 0.2)));
  }

  public void test_call_fun_copy_elision() {
    long old = GlobalOptions.getCopyElisionMinSize();
    GlobalOptions.setCopyElisionMinSize(0);
    try {
      MX mX = MX.sym("A", 3, 6);
      MX mY = MX.sym("Y");
      SX X = SX.sym("A", 3, 6);
      SX Y = SX.sym("Y");
      DM.rng(1);
      var g = ca.Function("g", new[] { X, Y }, new[] { X, Y });
      var gref = ca.Function("g", new[] { mX, mY }, new[] { mX, mY });
      List<DM> inputs = dms(DM.rand(3, 3), DM.rand(1));
      copyElisionChecks(g, gref, inputs);

      MX Xm = MX.sym("A", 3, 3);
      MX Ym = MX.sym("Y");
      var f = ca.Function("f", new[] { Xm, Ym }, new[] { ca.sumsqr(Xm) * Ym }, D(("never_inline", true)));
      var fref = ca.Function("f", new[] { Xm, Ym }, new[] { ca.sumsqr(Xm) * Ym });
      X = SX.sym("A", 3, 6);
      Y = SX.sym("Y");
      DM.rng(1);
      g = ca.Function("g", new[] { X, Y }, new[] { f.Call(X[.., 3..6], ca.sin(Y)), X[1, 3], 3.0 * X[2, 4] });
      gref = ca.Function("g", new[] { X, Y }, new[] { fref.Call(X[.., 3..6], ca.sin(Y)), X[1, 3], 3.0 * X[2, 4] });
      inputs = dms(DM.rand(3, 3), DM.rand(1));
      copyElisionChecks(g, gref, inputs);

      Xm = MX.sym("A", 4);
      Ym = MX.sym("Y");
      f = ca.Function("f", new[] { Xm, Ym }, new[] { ca.sumsqr(Xm) * Ym }, D(("never_inline", true)));
      fref = ca.Function("f", new[] { Xm, Ym }, new[] { ca.sumsqr(Xm) * Ym });
      X = SX.sym("A", 6);
      Y = SX.sym("Y");
      DM.rng(1);
      SX head = ca.vertcat(ca.SX(2.0), X[0..3]);
      g = ca.Function("g", new[] { X, Y }, new[] { f.Call(head, ca.sin(Y)) });
      gref = ca.Function("g", new[] { X, Y }, new[] { fref.Call(head, ca.sin(Y)) });
      inputs = dms(DM.rand(6), DM.rand(1));
      copyElisionChecks(g, gref, inputs);

      g = ca.Function("g", new[] { X, Y }, new[] { f.Call(head, ca.sin(Y)), 3.0 * X });
      gref = ca.Function("g", new[] { X, Y }, new[] { fref.Call(head, ca.sin(Y)), 3.0 * X });
      inputs = dms(DM.rand(6), DM.rand(1));
      copyElisionChecks(g, gref, inputs);
    } finally {
      GlobalOptions.setCopyElisionMinSize(old);
    }
  }

  void copyElisionChecks(Function g, Function gref, List<DM> inputs) {
    checkfunction(g, gref, inputs);
    foreach (bool avoid_stack in new[] { true, false }) checkCodegen(g, inputs, D(("avoid_stack", avoid_stack)));
    Function g_roundtrip = Function.deserialize(g.serialize());
    checkfunction(g, g_roundtrip, inputs);
    foreach (bool avoid_stack in new[] { true, false }) checkCodegen(g_roundtrip, inputs, D(("avoid_stack", avoid_stack)));
  }

  public void test_mmin() {
    SX x = SX.sym("X", 2);
    var f0 = ca.Function("f", new[] { x }, new[] { (x[0] + x[1]) / 2 });
    var f1 = ca.Function("f", new[] { x }, new[] { ca.mmin(x) });
    var f2 = ca.Function("f", new[] { x }, new[] { ca.fmin(x[0], x[1]) });
    checkfunction(f1, f2, dms(ca.vertcat(0.2, 0.3)));
    checkfunction(f1, f2, dms(ca.vertcat(2, 2)));
    checkfunction(f1, f0, dms(ca.vertcat(2, 2)));
    f1 = ca.Function("f", new[] { x }, new[] { ca.mmax(x) });
    f2 = ca.Function("f", new[] { x }, new[] { ca.fmax(x[0], x[1]) });
    checkfunction(f1, f2, dms(ca.vertcat(0.2, 0.3)));
    checkfunction(f1, f2, dms(ca.vertcat(2, 2)));
    checkfunction(f1, f0, dms(ca.vertcat(2, 2)));

    SX xs = SX.sym("X");
    f0 = ca.Function("f", new[] { xs }, new[] { xs });
    f1 = ca.Function("f", new[] { xs }, new[] { ca.fmin(xs, xs) });
    checkfunction(f1, f0, dms(ca.vertcat(1)));
    f1 = ca.Function("f", new[] { xs }, new[] { ca.fmax(xs, xs) });
    checkfunction(f1, f0, dms(ca.vertcat(1)));
  }

  public void test_empty_broadcast() {
    foreach (int nc in new[] { 0, 2 }) {
      SX res = ca.atan2(SX.sym("c", nc, 1), SX.sym("t", nc, 3));
      assertEquals((long)nc, res.size1());
      assertEquals(3L, res.size2());

      assertInException("Dimension mismatch", () => ca.atan2(MX.sym("c", nc, 2), MX.sym("t", nc, 3)));

      MX resm = ca.atan2(MX.sym("c", nc, 2), MX.sym("t", nc, 4));
      assertEquals((long)nc, resm.size1());
      assertEquals(4L, resm.size2());

      resm = ca.atan2(MX.sym("c", nc, 4), MX.sym("t", nc, 2));
      assertEquals((long)nc, resm.size1());
      assertEquals(4L, resm.size2());
    }
  }

  public void test_logsumexp() {
    SX x = SX.sym("x", 3);
    var f_ref = ca.Function("f_ref", new[] { x }, new[] { ca.log(ca.exp(x[0]) + ca.exp(x[1]) + ca.exp(x[2])) });
    var f = ca.Function("f", new[] { x }, new[] { ca.logsumexp(x) });

    checkfunction(f, f_ref, dms(ca.vertcat(1.1, 1.3, 1.7)));
    checkCodegen(f, dms(ca.vertcat(1.1, 1.3, 1.7)), D());
    checkfunction(f, f_ref, dms(ca.vertcat(1.1, 1.3, 1.3)));
    checkfunction(f, f_ref, dms(ca.vertcat(1.3, 1.3, 1.3)));

    checkarray(ca.logsumexp(ca.vertcat(1.3, 1.3, 1.3)), f_ref.Call(ca.vertcat(1.3, 1.3, 1.3)));
    checkarray(ca.logsumexp(ca.vertcat(1.1, 1.3, 1.3)), f_ref.Call(ca.vertcat(1.1, 1.3, 1.3)));
    checkarray(ca.logsumexp(ca.vertcat(1.1, 1.3, 1.7)), f_ref.Call(ca.vertcat(1.1, 1.3, 1.7)));

    DM res = f.Call(ca.vertcat(100, 1000, 10000));
    checkarray(res, 10000.0);
    checkarray(ca.logsumexp(ca.vertcat(100, 1000, 10000)), f.Call(ca.vertcat(100, 1000, 10000)));

    foreach (DM v in new[] { ca.vertcat(0, 1000, 1), ca.vertcat(1000, 0, 1), ca.vertcat(0, 1, 1000), ca.vertcat(5, 1e5, 3, 2, 1e5 - 1) }) {
      List<double> a = v.nonzeros();
      double max = a.Max();
      double sum = a.Sum(ai => Math.Exp(ai - max));
      double rf = max + Math.Log(sum);
      checkarray(ca.logsumexp(v), rf);
      MX xm = MX.sym("x", v.size1());
      var fm = ca.Function("fm", new[] { xm }, new[] { ca.logsumexp(xm) });
      checkarray(fm.Call(v), rf);
      checkCodegen(fm, dms(v), D());
    }
  }

  public void test_extract_parametric_call_sx() {
    MX xm = MX.sym("x");
    MX ym = MX.sym("y");
    MX zm = MX.sym("z");
    var f = ca.Function("f", new[] { xm, ym, zm }, new[] { xm * ym, ym * zm, xm * zm }, D(("never_inline", true)));

    SX x = SX.sym("x");
    SX p = SX.sym("p");
    List<SX> exprs = f.call(new[] { ca.SX(1.0), ca.sin(p), ca.sqrt(p) });
    SX expr = x * (exprs[0] + exprs[2]);

    var (expr_ret, symbols, parametric) = ca.extract_parametric(expr, p);
    Console.WriteLine("expr_ret " + expr_ret);
    Console.WriteLine("symbols " + string.Join(", ", symbols));
    Console.WriteLine("parametric " + string.Join(", ", parametric));
    assertTrue(expr_ret.ToString().Contains("x*e_0"));
    assertFalse(ca.depends_on(expr_ret, p));
    SX expr_recreated = ca.substitute(new[] { expr_ret }, symbols, parametric)[0];
    Console.WriteLine(expr_recreated - expr);
    assertTrue(ca.cse(expr_recreated - expr).is_zero());

    exprs = f.call(new[] { ca.SX(1.0), ca.sin(p), ca.pow(x, 2) });
    expr = exprs[0] + exprs[2];
    (expr_ret, symbols, parametric) = ca.extract_parametric(expr, p);
    assertFalse(ca.depends_on(expr_ret, p));
    expr_recreated = ca.substitute(new[] { expr_ret }, symbols, parametric)[0];

    string dir = Path.Combine(Path.GetTempPath(), "casadi_dotnet_sx_" + Guid.NewGuid().ToString("N"));
    Directory.CreateDirectory(dir);
    // f.generate('f1.c'), but into a temporary directory
    var cg = ca.CodeGenerator("f1");
    cg.add(ca.Function("f", new[] { x, p }, new[] { expr_recreated }));
    cg.generate(dir + "/");
    cg = ca.CodeGenerator("f2");
    cg.add(ca.Function("f", new[] { x, p }, new[] { expr }));
    cg.generate(dir + "/");
    Directory.Delete(dir, true);
    ca.cse(expr_recreated);
    assertTrue(ca.cse(expr_recreated - expr).is_zero());
  }

  public void test_check_recursion() {
    SX x = SX.sym("x");
    var f = ca.Function("f", new[] { x }, new[] { ca.pow(x, 2) }, D(("always_inline", true)));
    MX y = MX.sym("y");
    f.Call(y);
  }

  public void test_sx_eval_mx() {
    int n = 2;
    MX y = MX.sym("y", Sparsity.lower(n));
    var g = ca.Function("g", new[] { y }, new[] { ca.sin(y) }, D(("never_inline", true)));
    SX x = SX.sym("x", Sparsity.lower(n));
    var f = ca.Function("f", new[] { x }, new[] { g.Call(ca.pow(x, 2)) - 2 });
    MX X = MX.sym("X", Sparsity.upper(n));
    var F1 = ca.Function("F1", new[] { X }, f.call(new[] { X }, true));
    DM.rng(1);
    DM A = DM.rand(n, n);
    var F2 = ca.Function("F2", new[] { X }, f.call(new[] { X }));
    checkfunctionLight(F1, F2, dms(A));
  }

  public void test_contains() {
    SX x = SX.sym("x");
    SX y = SX.sym("y");
    SX z = SX.sym("z");
    SX e = y * z;
    assertTrue(ca.contains(new[] { x, y, z }, x));
    assertFalse(ca.contains(new[] { x, y }, z));
    assertTrue(ca.contains_any(new[] { x, y }, new[] { y, z }));
    assertFalse(ca.contains_all(new[] { x, y }, new[] { y, z }));
    assertTrue(ca.contains_any(new[] { x, y }, new[] { x, y }));
    assertTrue(ca.contains_all(new[] { x, y }, new[] { x, y }));
    assertTrue(ca.contains(new[] { e, x }, e));
    assertInException("Can only convert 1-by-1 matrices to scalars", () => ca.contains(new[] { ca.vertcat(x, y) }, x));
  }

  public void test_pow() {
    SX x = SX.sym("x");
    var y = ca.SX(4, 1);
    assertEquals(4L, ca.pow(y, 0).nnz());
    y[1] = x;
    assertEquals(4L, ca.pow(y, 0).nnz());
    y[1] = 0;
    assertEquals(4L, ca.pow(y, 0).nnz());
    y[1] = 1;
    assertEquals(4L, ca.pow(y, 0).nnz());
  }

  public void test_set_precision() {
    double pi_val = 3.141592653589793239;
    string expected_default = "3.14159";
    string expected_prec12 = "3.14159265359";

    var xd = ca.DM(pi_val);
    try {
      assertEquals(expected_default, xd.ToString());
      DM.set_precision(12);
      assertEquals(expected_prec12, xd.ToString());
      DM.set_precision(4);
      DM.set_scientific(true);
      assertTrue(xd.ToString().ToLower().Contains("e"));
    } finally {
      DM.set_precision(6);
      DM.set_scientific(false);
    }
    var xs = ca.SX(pi_val);
    try {
      assertEquals(expected_default, xs.ToString());
      SX.set_precision(12);
      assertEquals(expected_prec12, xs.ToString());
      SX.set_precision(4);
      SX.set_scientific(true);
      assertTrue(xs.ToString().ToLower().Contains("e"));
    } finally {
      SX.set_precision(6);
      SX.set_scientific(false);
    }

    var mx = ca.MX(pi_val);
    try {
      assertEquals(expected_default, mx.ToString());
      DM.set_precision(12);
      assertEquals(expected_prec12, mx.ToString());
    } finally {
      DM.set_precision(6);
    }

    DM vd = ca.vertcat(1.234567890123, 9.876543210987);
    try {
      assertEquals("[1.23457, 9.87654]", vd.ToString());
      DM.set_precision(12);
      assertEquals("[1.23456789012, 9.87654321099]", vd.ToString());
    } finally {
      DM.set_precision(6);
    }
    SX vs = new double[] { 1.234567890123, 9.876543210987 };
    try {
      assertEquals("[1.23457, 9.87654]", vs.ToString());
      SX.set_precision(12);
      assertEquals("[1.23456789012, 9.87654321099]", vs.ToString());
    } finally {
      SX.set_precision(6);
    }
  }

  public void test_linearize() {
    SX x = SX.sym("x");
    SX y = ca.sin(x);
    double x0 = 0.2;
    Console.WriteLine(ca.linearize(ca.sin(x), x, ca.SX(x0)));
    Console.WriteLine(Math.Sin(x0) + Math.Cos(x0) * (x - x0));
    var F = ca.Function("F", new[] { x }, new[] { y });
    Function Ff = F.forward(1);
    DM n = F.Call((DM)x0);
    Console.WriteLine(n + Ff.Call((SX)x0, (SX)n, x - x0));
    Console.WriteLine(ca.taylor(y, x, x0));
  }

  public void test_gauss_quadrature() {
    SX x = SX.sym("x");
    checkarray(ca.evalf(ca.gauss_quadrature(ca.pow(x, 2), x, 0, 1)), ca.DM(1.0 / 3));
    checkarray(ca.evalf(ca.gauss_quadrature(ca.pow(x, 9), x, -1, 1)), ca.DM(0.0));
  }

  public void test_poly_roots_order() {
    assertRaisesRegex("order at most 4.*Got order 5", () => ca.poly_roots((SX)new double[] { 1, 2, 3, 4, 5, 6 }));
  }

  public void test_mtaylor_dimension() {
    SX x = SX.sym("x");
    SX y = SX.sym("y");
    assertRaisesRegex("expected a of shape 2x1", () => ca.mtaylor(x * y, ca.vertcat(x, y), 0, 1));
  }

  public void test_large_whole_constants() {
    foreach (double v in new[] { 3e9, 6.02214076e23, -1e20, Math.Pow(2, 31), -Math.Pow(2, 31) - 1, Math.Pow(2, 63), 1e300 }) {
      assertEquals(v, (double)ca.evalf(ca.SX(v)));
      assertEquals(v, (double)ca.evalf(ca.MX(v)));
    }
    foreach (double v in new[] { Math.Pow(2, 31) - 1, -Math.Pow(2, 31) }) {
      assertTrue(ca.SX(v).is_integer());
      assertEquals(v, (double)ca.evalf(ca.SX(v)));
    }
  }
}

}

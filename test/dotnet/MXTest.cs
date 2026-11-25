using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;

namespace CasADi.Tests {

/// <summary>Port of test/java MXTest (test/python/mx.py, first part, up to test_project_simplify)</summary>
public class MXTest : CasadiTestCase {

  static readonly bool KNOWN_BUGS = Environment.GetEnvironmentVariable("CASADI_JAVA_KNOWN_BUGS") != null;

  /// <summary>FunctionPool entry: casadi operation on MX, reference on DM</summary>
  sealed class Op1 {
    public readonly string name;
    public readonly bool nozero;
    public readonly Func<MX, MX> cas;
    public readonly Func<DM, DM> @ref;
    public Op1(string name, bool nozero, Func<MX, MX> cas, Func<DM, DM> @ref) {
      this.name = name;
      this.nozero = nozero;
      this.cas = cas;
      this.@ref = @ref;
    }
  }

  sealed class Op2 {
    public readonly string name;
    public readonly Func<MX, MX, MX> cas;
    public readonly Func<DM, DM, DM> @ref;
    public Op2(string name, Func<MX, MX, MX> cas, Func<DM, DM, DM> @ref) {
      this.name = name;
      this.cas = cas;
      this.@ref = @ref;
    }
  }

  List<Op1> pool, Jpool;
  readonly List<Op2> matrixbinarypool;

  public MXTest() {
    pool = new List<Op1> {
      new Op1("sqrt", false, x => ca.sqrt(x), x => ca.sqrt(x)),
      new Op1("sin", false, x => ca.sin(x), x => ca.sin(x)),
      new Op1("cos", false, x => ca.cos(x), x => ca.cos(x)),
      new Op1("tan", false, x => ca.tan(x), x => ca.tan(x)),
      new Op1("arctan", false, x => ca.atan(x), x => ca.atan(x)),
      new Op1("arcsin", false, x => ca.asin(x), x => ca.asin(x)),
      new Op1("arccos", false, x => ca.acos(x), x => ca.acos(x)),
      new Op1("exp", false, x => ca.exp(x), x => ca.exp(x)),
      new Op1("log", true, x => ca.log(x), x => ca.log(x)),
      new Op1("x^0", true, x => ca.pow(x, 0.0), x => ca.pow(x, 0.0)),
      new Op1("^1", false, x => ca.pow(x, 1.0), x => ca.pow(x, 1.0)),
      new Op1("^-2", true, x => ca.pow(x, -2.0), x => ca.pow(x, -2.0)),
      new Op1("^0.3", false, x => ca.pow(x, 0.3), x => ca.pow(x, 0.3)),
      new Op1("floor", false, x => ca.floor(x), x => ca.floor(x)),
      new Op1("ceil", false, x => ca.ceil(x), x => ca.ceil(x)),
    };
    Jpool = new List<Op1> {
      new Op1("sqrt", false, x => ca.sqrt(x), x => ca.diag(1.0 / (ca.sqrt(x) * 2.0))),
      new Op1("sin", false, x => ca.sin(x), x => ca.diag(ca.cos(x))),
      new Op1("cos", false, x => ca.cos(x), x => ca.diag(-ca.sin(x))),
      new Op1("tan", false, x => ca.tan(x), x => ca.diag(1.0 / ca.pow(ca.cos(x), 2.0))),
      new Op1("arctan", false, x => ca.atan(x), x => ca.diag(1.0 / (ca.pow(x, 2.0) + 1.0))),
      new Op1("arcsin", false, x => ca.asin(x), x => ca.diag(1.0 / ca.sqrt(1.0 - ca.pow(x, 2.0)))),
      new Op1("arccos", false, x => ca.acos(x), x => ca.diag(-1.0 / ca.sqrt(1.0 - ca.pow(x, 2.0)))),
      new Op1("exp", false, x => ca.exp(x), x => ca.diag(ca.exp(x))),
      new Op1("log", false, x => ca.log(x), x => ca.diag(1.0 / x)),
      new Op1("x^0", false, x => ca.pow(x, 0.0), x => ca.diag(DM.zeros(x.size1(), x.size2()))),
      new Op1("^1", false, x => ca.pow(x, 1.0), x => ca.diag(DM.ones(x.size1(), x.size2()))),
      new Op1("^-2", false, x => ca.pow(x, -2.0), x => ca.diag(-2.0 / ca.pow(x, 3.0))),
      new Op1("^0.3", false, x => ca.pow(x, 0.3), x => ca.diag(0.3 / ca.pow(x, 0.7))),
    };
    matrixbinarypool = new List<Op2> {
      new Op2("Matrix+Matrix", (a, b) => a + b, (a, b) => a + b),
      new Op2("Matrix-Matrix", (a, b) => a - b, (a, b) => a - b),
      new Op2("Matrix*Matrix", (a, b) => a * b, (a, b) => a * b),
      new Op2("fmin", (a, b) => ca.fmax(a, b), (a, b) => ca.fmax(a, b)),
      new Op2("fmax", (a, b) => ca.fmin(a, b), (a, b) => ca.fmin(a, b)),
      new Op2("mtimes(Matrix,Matrix.T)", (a, b) => ca.mtimes(a, b.T()), (a, b) => ca.mtimes(a, b.T())),
      new Op2("arctan2", (a, b) => ca.atan2(a, b), (a, b) => ca.atan2(a, b)),
      new Op2("mtimes(Matrix,Matrix.T)", (a, b) => ca.mtimes(a, b.T()), (a, b) => ca.mtimes(a, b.T())),
    };
  }

  // ---------------------------------------------------------------- helpers

  static DM dmsp(Sparsity sp, params double[] v) => ca.DM(sp, ca.DM(v));

  static DM idx(params double[] i) => ca.DM(i);

  static double[] range(int n) => Enumerable.Range(0, n).Select(i => (double)i).ToArray();

  static DM ev(Function f, params DM[] args) => f.Call(args);
  static List<DM> evs(Function f, params DM[] args) => f.call(args);
  static MX ev(Function f, params MX[] args) => f.Call(args);
  static List<MX> evs(Function f, params MX[] args) => f.call(args);

  static Function fun(string name, IEnumerable<MX> @in, IEnumerable<MX> @out) => ca.Function(name, @in, @out);

  static Function jacobianOld(Function f, int i, int j) {
    var @out = new List<string> { "jac:" + f.name_out(j) + ":" + f.name_in(i) };
    @out.AddRange(f.name_out());
    return f.factory(f.name() + "_jac", f.name_in(), @out);
  }

  static double[,] npT(double[,] a) {
    int n = a.GetLength(0), m = a.GetLength(1);
    var r = new double[m, n];
    for (int i = 0; i < n; i++)
      for (int j = 0; j < m; j++) r[j, i] = a[i, j];
    return r;
  }

  /// <summary>numpy.reshape, row-major</summary>
  static double[,] npReshape(double[,] a, int n, int m) {
    double[] flat = a.Cast<double>().ToArray();  // row-major enumeration
    var r = new double[n, m];
    for (int i = 0; i < n; i++)
      for (int j = 0; j < m; j++) r[i, j] = flat[i * m + j];
    return r;
  }

  static double[,] npScale(double[,] a, double s) {
    var r = (double[,])a.Clone();
    for (int i = 0; i < a.GetLength(0); i++)
      for (int j = 0; j < a.GetLength(1); j++) r[i, j] *= s;
    return r;
  }

  static double[,] npVstack(double[,] a, double[,] b) {
    int na = a.GetLength(0), nb = b.GetLength(0), m = a.GetLength(1);
    var r = new double[na + nb, m];
    for (int j = 0; j < m; j++) {
      for (int i = 0; i < na; i++) r[i, j] = a[i, j];
      for (int i = 0; i < nb; i++) r[na + i, j] = b[i, j];
    }
    return r;
  }

  static double[,] npHstack(double[,] a, double[,] b) => npT(npVstack(npT(a), npT(b)));

  static int Len(double[,] x) => x.GetLength(0) * x.GetLength(1);

  static void checkMXoperations(Func<MX, MX> ztf, Func<double[,], double[,]> zrf, string name) {
    MX x = MX.sym("x", 1, 3);
    var zs = new List<MX>();
    for (int i = 0; i < 8; i++) zs.Add(x * (i + 1.0));
    MX z = ca.vertcat(zs);
    Function f = fun("f", new[] { x }, new[] { ztf(z) });
    double[] L = { 1, 2, 3 };
    DM fOut = ev(f, (DM)new double[,] { { L[0], L[1], L[2] } });
    var zr = new double[8, 3];
    for (int i = 0; i < 8; i++)
      for (int j = 0; j < 3; j++) zr[i, j] = L[j] * (i + 1);
    checkarray(fOut, (DM)zrf(zr), name);
  }

  static void checkMXoperations2(Func<MX, MX> ztf, Func<double[,], double[,]> zrf, string name) {
    MX x = MX.sym("x", 3, 1);
    var zs = new List<MX>();
    for (int i = 0; i < 8; i++) zs.Add(x * (double)i);
    MX z = ca.horzcat(zs);
    Function f = fun("f", new[] { x }, new[] { ztf(z) });
    double[] L = { 1, 2, 3 };
    DM fOut = ev(f, (DM)L);
    var zr = new double[8, 3];
    for (int i = 0; i < 8; i++)
      for (int j = 0; j < 3; j++) zr[i, j] = L[j] * i;
    checkarray(fOut, (DM)zrf(npT(zr)), name);
  }

  static void checkMXoperations3(Func<MX, MX> ztf, Func<double[,], double[,]> zrf, string name) {
    MX x = MX.sym("x", 3, 1);
    MX p = ca.horzcat(x[0, 0], x[1, 0], x[2, 0]);
    var zs = new List<MX>();
    for (int i = 0; i < 8; i++) zs.Add(p * (double)i);
    MX z = ca.vertcat(zs);
    Function f = fun("f", new[] { x }, new[] { ztf(z) });
    double[] L = { 1, 2, 3 };
    DM fOut = ev(f, (DM)L);
    var zr = new double[8, 3];
    for (int i = 0; i < 8; i++)
      for (int j = 0; j < 3; j++) zr[i, j] = L[j] * i;
    checkarray(fOut, (DM)zrf(zr), name);
  }

  /// <summary>helpers.py evaluationCheck/numpyEvaluationCheck</summary>
  void numpyEvaluationCheck(Func<MX, MX> ft, DM @ref, List<MX> x, List<DM> setx0, string name,
                            Func<Function, Function> fmod) {
    message(":" + name);
    Function f = fun("f", x, new[] { ft(x[0]) });
    if (fmod != null) f = fmod(f);
    var fIn = new List<DM>(setx0);
    while (fIn.Count < f.n_in()) fIn.Add(ca.DM());
    DM zt = ca.densify(f.call(fIn)[0]);
    checkarray(zt, @ref, name);
  }

  void numpyEvaluationCheckPool(List<Op1> pool, MX x, DM x0, string name, Func<Function, Function> fmod,
                                DM setx0, bool excludeNozero) {
    foreach (Op1 op in pool) {
      if (excludeNozero && op.nozero) continue;
      numpyEvaluationCheck(op.cas, op.@ref(ca.densify(x0)), new List<MX> { x }, new List<DM> { setx0 ?? x0 },
                           name + ":" + op.name, fmod);
    }
  }

  void numpyEvaluationCheckPool2(List<Op2> pool, MX x, MX y, DM x0, DM y0, string name) {
    foreach (Op2 op in pool) {
      message(":" + name + ":" + op.name);
      Function f = fun("f", new[] { x, y }, new[] { op.cas(x, y) });
      DM zt = ca.densify(f.call(new[] { x0, y0 })[0]);
      checkarray(zt, op.@ref(ca.densify(x0), ca.densify(y0)), name + ":" + op.name);
    }
  }

  /// <summary>helpers.py checkfunction_light: outputs agree</summary>
  static void checkfunctionLight(Function f, Function g, List<DM> inputs) {
    List<DM> a = f.call(inputs), b = g.call(inputs);
    assertEquals(a.Count, b.Count);
    for (int i = 0; i < a.Count; i++) checkarray(a[i], b[i], 1e-9, "output " + i);
  }

  static void checkSerialize(Function f, List<DM> inputs) {
    Function f2 = Function.deserialize(f.serialize(new Dictionary<string, object> { ["debug"] = true }));
    List<DM> a = f.call(inputs), b = f2.call(inputs);
    for (int i = 0; i < f.n_out(); i++) checkarray(a[i], b[i], 1e-15, "serialize output " + i);
  }

  void checkCodegen(Function f, List<DM> inputs) => checkCodegen(f, inputs, new Dictionary<string, object>(), new List<string>());

  /// <summary>helpers.py check_codegen: generate, compile with gcc, load as external, compare (RUN_SLOW only)</summary>
  void checkCodegen(Function f, List<DM> inputs, Dictionary<string, object> opts, List<string> extralibs) {
    if (!RunSlow) return;
    string dir = tmpDir();
    string name = "codegen_" + ((ulong)random.NextInt64(long.MaxValue)).ToString("x");
    var cg = ca.CodeGenerator(name, opts);
    cg.add(f);
    cg.generate(dir + Path.DirectorySeparatorChar);
    string libdir = GlobalOptions.getCasadiPath();
    string inc = GlobalOptions.getCasadiIncludePath();
    var cmd = new List<string> { "-pedantic", "-std=c89", "-fPIC", "-shared", "-Wall",
      "-Werror", "-Wextra", "-I" + inc, "-I" + inc + "/highs", "-Wno-unknown-pragmas", "-Wno-long-long",
      "-Wno-unused-parameter", "-O3", name + ".c", "-o", name + ".so", "-L" + libdir, "-Wl,-rpath," + libdir };
    foreach (string lib in extralibs) cmd.Add(lib.StartsWith("-") || lib.Contains(".") ? lib : "-l" + lib);
    var psi = new ProcessStartInfo("gcc") { WorkingDirectory = dir, UseShellExecute = false };
    foreach (string c in cmd) psi.ArgumentList.Add(c);
    using (var p = Process.Start(psi)) {
      p.WaitForExit();
      assertEquals(0, p.ExitCode, "compile gcc " + string.Join(" ", cmd));
    }
    Function f2 = ca.external(f.name(), Path.Combine(dir, name + ".so"));
    List<DM> a = f.call(inputs), b = f2.call(inputs);
    for (int i = 0; i < f.n_out(); i++) checkarray(a[i], b[i], 1e-15, "codegen output " + i);
  }

  static string tmpDir() {
    string d = Path.Combine(Path.GetTempPath(), "casadi_mx_" + Guid.NewGuid().ToString("N"));
    Directory.CreateDirectory(d);
    return d;
  }

  static string readFile(string p) => File.ReadAllText(p).Replace("\r\n", "\n");

  // ---------------------------------------------------------------- tests

  public void test_MX1() {
    message("MX constructor");
    MX x = MX.sym("x", 2, 3);
    assertEquals(2L, x.size1(), "MX fails to indicate its size1");
    assertEquals(3L, x.size2(), "MX fails to indicate its size2");
  }

  public void test_MXvertcat() {
    message("MX vertcat");
    MX x = MX.sym("x", 1, 3);
    MX y = MX.sym("y", 1, 3);
    MX z = ca.vertcat(x, y);
    assertEquals(2L, z.size1(), "MX fails to indicate its size1");
    assertEquals(3L, z.size2(), "MX fails to indicate its size2");
  }

  public void test_MX_fun1() {
    message("MXFunction single input, single output");
    MX x = MX.sym("x");
    MX y = 2.0 * x;
    Function f = fun("f", new[] { x }, new[] { y });
    assertEquals(1L, f.n_in(), "Function fails to indicate correct number of inputs");
    assertEquals(1L, f.n_out(), "Function fails to indicate correct number of outputs");
    DM fOut = ev(f, 3);
    List<double> yt = fOut.nonzeros();
    assertEquals(1, yt.Count, "Output of Function was tuple of floats, as expected, but length is incorrect.");
    checkarray(yt[0], 2 * 3, 1e-10);
  }

  public void test_MXfunction2() {
    message("Function multi input, multi output");
    MX x = MX.sym("x");
    MX y = MX.sym("y");
    Function f = fun("f", new[] { x, y }, new[] { x + y, y * x });
    assertEquals(2L, f.n_in(), "Function fails to indicate correct number of inputs");
    assertEquals(2L, f.n_out(), "Function fails to indicate correct number of outputs");
    List<DM> fOut = evs(f, 3, 7);
    List<double> zt1 = fOut[0].nonzeros(), zt2 = fOut[1].nonzeros();
    assertEquals(1, zt1.Count);
    assertEquals(1, zt2.Count);
    checkarray(zt2[0], 21, 1e-10);
    checkarray(zt1[0], 10, 1e-10);
  }

  public void test_MXfunction3() {
    message("Function single input, multi output (1)");
    MX xy = MX.sym("xy", 2);
    Function f = fun("f", new[] { xy }, new[] { xy[0] + xy[1], xy[0] * xy[1] });
    assertEquals(1L, f.n_in(), "Function fails to indicate correct number of inputs");
    assertEquals(2L, f.n_out(), "Function fails to indicate correct number of outputs");
    List<DM> fOut = evs(f, new double[] { 3, 7 });
    List<double> zt1 = fOut[0].nonzeros(), zt2 = fOut[1].nonzeros();
    assertEquals(1, zt1.Count);
    assertEquals(1, zt2.Count);
    checkarray(zt2[0], 21, 1e-10);
    checkarray(zt1[0], 10, 1e-10);
  }

  public void test_MXfunction3b() {
    message("Function single input, multi output (2)");
    MX xy = MX.sym("xy", 1, 2);
    Function f = fun("f", new[] { xy }, new[] { xy[0, 0] + xy[0, 1], xy[0, 0] * xy[0, 1] });
    assertEquals(1L, f.n_in(), "Function fails to indicate correct number of inputs");
    assertEquals(2L, f.n_out(), "Function fails to indicate correct number of outputs");
    List<DM> fOut = evs(f, new double[,] { { 3, 7 } });
    double[,] zt1 = fOut[0].ToArray(), zt2 = fOut[1].ToArray();
    assertEquals(1, zt1.GetLength(0), "Output of Function is of wrong shape.");
    assertEquals(1, zt1.GetLength(1), "Output of Function is of wrong shape.");
    assertEquals(1, zt2.GetLength(0), "Output of Function is of wrong shape.");
    assertEquals(1, zt2.GetLength(1), "Output of Function is of wrong shape.");
    checkarray(zt2[0, 0], 21, 1e-10);
    checkarray(zt1[0, 0], 10, 1e-10);
  }

  public void test_MXfunction4() {
    message("Function single input, single output , using vertcat");
    MX xy = MX.sym("xy", 2);
    MX z = ca.vertcat(xy[0] + xy[1], xy[0] * xy[1]);
    Function f = fun("f", new[] { xy }, new[] { z });
    assertEquals(1L, f.n_in(), "Function fails to indicate correct number of inputs");
    assertEquals(1L, f.n_out(), "Function fails to indicate correct number of outputs");
    double[,] zt = ev(f, new double[] { 3, 7 }).ToArray();
    assertEquals(2, zt.GetLength(0), "Output of Function is of wrong shape.");
    assertEquals(1, zt.GetLength(1), "Output of Function is of wrong shape.");
    checkarray(zt[1, 0], 21, 1e-10);
    checkarray(zt[0, 0], 10, 1e-10);
  }

  public void test_MXfunction5() {
    message("Function single input, single output , using horzcat");
    MX xy = MX.sym("xy", 2);
    MX z = ca.horzcat(xy[0] + xy[1], xy[0] * xy[1]);
    Function f = fun("f", new[] { xy }, new[] { z });
    assertEquals(1L, f.n_in(), "Function fails to indicate correct number of inputs");
    assertEquals(1L, f.n_out(), "Function fails to indicate correct number of outputs");
    double[,] zt = ev(f, new double[] { 3, 7 }).ToArray();
    assertEquals(1, zt.GetLength(0), "Output of Function is of wrong shape.");
    assertEquals(2, zt.GetLength(1), "Output of Function is of wrong shape.");
    checkarray(zt[0, 1], 21, 1e-10);
    checkarray(zt[0, 0], 10, 1e-10);
  }

  public void test_which_depends_empty() {
    SX xs = SX.sym("x");
    MX xm = MX.sym("x");
    foreach (bool tr in new[] { true, false }) {
      for (int i = 0; i < 3; i++) {
        assertEquals(tr ? 1 : 0, ca.which_depends(xs, ca.SX(0, 1), i, tr).Count);
        assertFalse(ca.which_depends(xs, ca.SX(0, 1), i, tr).Any(b => b));
        assertEquals(tr ? 0 : 1, ca.which_depends(ca.SX(0, 1), xs, i, tr).Count);
        assertFalse(ca.which_depends(ca.SX(0, 1), xs, i, tr).Any(b => b));
        assertEquals(0, ca.which_depends(ca.SX(0, 1), ca.SX(0, 1), i, tr).Count);

        assertEquals(tr ? 1 : 0, ca.which_depends(xm, ca.MX(0, 1), i, tr).Count);
        assertFalse(ca.which_depends(xm, ca.MX(0, 1), i, tr).Any(b => b));
        assertEquals(tr ? 0 : 1, ca.which_depends(ca.MX(0, 1), xm, i, tr).Count);
        assertFalse(ca.which_depends(ca.MX(0, 1), xm, i, tr).Any(b => b));
        assertEquals(0, ca.which_depends(ca.MX(0, 1), ca.MX(0, 1), i, tr).Count);
      }
    }
  }

  public void test_issue83() {
    MX x = MX.sym("x");
    MX y = MX.sym("y");
    MX z = x + y;
    Function f = fun("f", new[] { x, y }, new[] { z });
    MX fc = f.Call(ca.MX(3), y);
    Function g = fun("g", new[] { y }, new[] { fc });
    DM gOut = ev(g, 7);
    checkarray((double)gOut, 10, 1e-10);

    fc = f.Call(x, ca.MX(7));
    g = fun("g", new[] { x }, new[] { fc });
    gOut = ev(g, 3);
    checkarray((double)gOut, 10, 1e-10);
  }

  public void test_identitySX() {
    message("identity SXFunction");
    SX x = SX.sym("x");
    Function f = ca.Function("f", new[] { x }, new[] { x });
    DM fOut = ev(f, 3);
    checkarray((double)fOut[0, 0], 3, 1e-10);
  }

  public void test_identityMX() {
    message("identity Function");
    MX x = MX.sym("x");
    Function f = fun("f", new[] { x }, new[] { x });
    DM fOut = ev(f, 3);
    checkarray((double)fOut[0, 0], 3, 1e-10);
  }

  public void test_MXorder() {
    message("Function order of non-zero elements");
    MX x = MX.sym("x", 2, 3);
    Function f = fun("f", new[] { x }, new[] { x + x });
    assertEquals(1L, f.n_in(), "Function fails to indicate correct number of inputs");
    assertEquals(1L, f.n_out(), "Function fails to indicate correct number of outputs");
    double[] L = { 1, 2, 3, 4, 5, 6 };
    DM fIn = ca.DM(f.sparsity_in(0), L);
    double[,] zt = ev(f, fIn).ToArray();
    assertEquals(2, zt.GetLength(0), "Output of Function is of wrong shape.");
    assertEquals(3, zt.GetLength(1), "Output of Function is of wrong shape.");
    for (int i = 0; i < 2; i++)
      for (int j = 0; j < 3; j++) checkarray(L[i + 2 * j] * 2, zt[i, j], 1e-10);
  }

  public void test_trans() {
    message("trans");
    MX a = ca.MX(0, 1);
    MX b = a.T();
    assertEquals(1L, b.size1());
    assertEquals(0L, b.size2());
  }

  public void test_MXtrans() {
    message("trans(MX)");
    MX x = MX.sym("x", 2, 3);
    MX z = x.T();
    assertEquals(3L, z.size1(), "Vec returns MX of wrong dimension");
    assertEquals(2L, z.size2(), "Vec returns MX of wrong dimension");
    Function f = fun("f", new[] { x }, new[] { z });
    assertEquals(1L, f.n_in(), "Function fails to indicate correct number of inputs");
    assertEquals(1L, f.n_out(), "Function fails to indicate correct number of outputs");
    double[] L = { 1, 2, 3, 4, 5, 6 };
    DM fIn = ca.DM(f.sparsity_in(0), L);
    double[,] ztr = ev(f, fIn).ToArray();
    for (int i = 0; i < 2; i++)
      for (int j = 0; j < 3; j++) checkarray(L[i + 2 * j], ztr[j, i], 1e-10);
  }

  public void test_MXvec() {
    var u0 = new double[4, 3];
    for (int j = 0; j < 4; j++)
      for (int i = 0; i < 3; i++) u0[j, i] = 10 * j + i;
    DM u = u0;
    MX U = MX.sym("u", u.size());
    Function f = fun("f", new[] { U }, new[] { ca.vec(U) });
    DM fOut = ev(f, u);
    checkarray(ca.vec(u), fOut, "vec");
  }

  public void test_MXreshape() {
    message("reshape(MX)");
    MX x = MX.sym("x", 2, 3);
    MX z = ca.reshape(x, (1, 6));
    assertEquals(1L, z.size1(), "Vec returns MX of wrong dimension");
    assertEquals(6L, z.size2(), "Vec returns MX of wrong dimension");
    Function f = fun("f", new[] { x }, new[] { z });
    assertEquals(1L, f.n_in(), "Function fails to indicate correct number of inputs");
    assertEquals(1L, f.n_out(), "Function fails to indicate correct number of outputs");
    double[] L = { 1, 2, 3, 4, 5, 6 };
    DM fIn = ca.DM(f.sparsity_in(0), L);
    double[,] zt = ev(f, fIn).ToArray();
    for (int i = 0; i < L.Length; i++) checkarray(L[i], zt[0, i], 1e-10);
  }

  public void test_sparsity_cast() {
    Sparsity spSource = ca.sparsify(ca.DM(new double[,] { { 1, 0, 1 }, { 0, 0, 1 }, { 1, 0, 0 } })).sparsity();
    MX x = MX.sym("x", spSource);
    SX xsx = SX.sym("x", spSource);
    Sparsity sp = ca.sparsify(ca.DM(new double[,] { { 1, 0, 1 }, { 1, 0, 1 } })).sparsity();
    assertInException("mismatch", () => ca.reshape(x, sp));
    assertInException("Mismatching", () => ca.sparsity_cast(x, ca.horzcat(new[] { sp, sp })));
    MX y = ca.sparsity_cast(x, sp);
    SX ysx = ca.sparsity_cast(xsx, sp);
    Function fsx = ca.Function("fsx", new[] { xsx }, new[] { ysx });
    Function f = fun("f", new[] { x }, new[] { y });
    DM inp = ca.DM(new double[,] { { 1, 0, 2 }, { 0, 0, 3 }, { 4, 0, 0 } });
    checkfunction(f, fsx, dms(inp));
    checkCodegen(f, dms(inp));
    DM y2 = ev(fsx, inp);
    checkarray(y2, ca.DM(new double[,] { { 1, 0, 2 }, { 4, 0, 3 } }));

    Sparsity sp2 = ca.sparsify(ca.DM(new double[,] { { 1, 0 }, { 0, 1 } })).sparsity();
    MX y3 = ca.sparsity_cast(MX.sym("y", 2), sp2);
    ca.vec(y3);
  }

  public void test_MXcompose() {
    message("compositions of vec, trans, reshape with vertcat");
    checkMXoperations(x => x, x => x, "vertcat");
    checkMXoperations(x => x.T(), npT, "trans(vertcat)");
    checkMXoperations(x => x.T().T(), x => x, "trans(trans(vertcat))");
    checkMXoperations(x => ca.vec(x.T()), x => npReshape(x, Len(x), 1), "vec(trans(vertcat))");
    checkMXoperations(x => ca.vec(x).T(), x => npT(npReshape(npT(x), Len(x), 1)), "vec(trans(vertcat))");
    checkMXoperations(x => ca.reshape(x.T(), 6, 4).T(), x => npReshape(x, 4, 6), "reshape(vertcat)");
    checkMXoperations(x => ca.reshape(x, 6, 4).T(), x => npReshape(npT(x), 4, 6), "reshape(trans(vertcat))");
    checkMXoperations(x => ca.reshape(x.T(), 6, 4), x => npT(npReshape(x, 4, 6)), "trans(reshape(vertcat))");
  }

  public void test_MXcompose2() {
    message("compositions of vec, trans, reshape with horzcat");
    checkMXoperations2(x => x, x => x, "horzcat");
    checkMXoperations2(x => x.T(), npT, "trans(horzcat)");
    checkMXoperations2(x => x.T().T(), x => x, "trans(trans(horzcat))");
    checkMXoperations2(x => ca.vec(x.T()), x => npReshape(x, Len(x), 1), "vec(trans(horzcat))");
    checkMXoperations2(x => ca.vec(x).T(), x => npT(npReshape(npT(x), Len(x), 1)), "vec(trans(horzcat))");
    checkMXoperations2(x => ca.reshape(x.T(), 6, 4).T(), x => npReshape(x, 4, 6), "reshape(horzcat)");
    checkMXoperations2(x => ca.reshape(x, 6, 4).T(), x => npReshape(npT(x), 4, 6), "reshape(trans(horzcat))");
    checkMXoperations2(x => ca.reshape(x.T(), 6, 4), x => npT(npReshape(x, 4, 6)), "trans(reshape(horzcat))");
  }

  public void test_MXcompose3() {
    message("compositions of vec, trans, reshape with vertcat");
    checkMXoperations3(x => x, x => x, "snippet");
    checkMXoperations3(x => x.T(), npT, "trans(snippet)");
    checkMXoperations3(x => x.T().T(), x => x, "trans(trans(snippet))");
    checkMXoperations3(x => ca.vec(x.T()), x => npReshape(x, Len(x), 1), "vec(trans(snippet))");
    checkMXoperations3(x => ca.vec(x).T(), x => npT(npReshape(npT(x), Len(x), 1)), "vec(trans(snippet))");
    checkMXoperations3(x => ca.reshape(x.T(), 6, 4).T(), x => npReshape(x, 4, 6), "reshape(snippet)");
    checkMXoperations3(x => ca.reshape(x, 6, 4).T(), x => npReshape(npT(x), 4, 6), "reshape(trans(snippet))");
    checkMXoperations3(x => ca.reshape(x.T(), 6, 4), x => npT(npReshape(x, 4, 6)), "trans(reshape(snippet))");
  }

  public void test_MXcompose4() {
    message("compositions of horzcat + vertcat");
    checkMXoperations(x => ca.vertcat(x), x => x, "vertcat(*vertcat)");
    checkMXoperations(x => ca.vertcat(x, x * 2.0), x => npVstack(x, npScale(x, 2)), "vertcat(*vertcat,vertcat)");
    checkMXoperations(x => ca.horzcat(x), x => x, "horzcat(*vertcat)");
    checkMXoperations(x => ca.horzcat(x, x * 2.0), x => npHstack(x, npScale(x, 2)), "horzcat(*vertcat,vertcat)");

    checkMXoperations2(x => ca.vertcat(x), x => x, "vertcat(*horzcat)");
    checkMXoperations2(x => ca.vertcat(x, x * 2.0), x => npVstack(x, npScale(x, 2)), "vertcat(*horzcat,horzcat)");
    checkMXoperations2(x => ca.horzcat(x), x => x, "horzcat(*horzcat)");
    checkMXoperations2(x => ca.horzcat(x, x * 2.0), x => npHstack(x, npScale(x, 2)), "horzcat(*horzcat,horzcat)");

    checkMXoperations3(x => ca.vertcat(x), x => x, "vertcat(*snippet)");
    checkMXoperations3(x => ca.vertcat(x, x * 2.0), x => npVstack(x, npScale(x, 2)), "vertcat(*snippet,snippet)");
    checkMXoperations3(x => ca.horzcat(x), x => x, "horzcat(*snippet)");
    checkMXoperations3(x => ca.horzcat(x, x * 2.0), x => npHstack(x, npScale(x, 2)), "horzcat(*snippet,snippet)");
  }

  public void test_MXslicingnew() {
    Check.Assume(KNOWN_BUGS, "known bug");
    message("MX slicing new");
    message(":dense");
    MX x = MX.sym("x", 3, 2);
    DM x0 = ca.DM(new double[,] { { 0.738, 0.2 }, { 0.1, 0.39 }, { 0.99, 0.999999 } });
    var X = new List<MX> { x };
    var X0 = new List<DM> { x0 };
    numpyEvaluationCheck(a => a[1, 0], x0[1, 0], X, X0, "x[1,0]", null);
    numpyEvaluationCheck(a => a[0, 1], x0[0, 1], X, X0, "x[0,1]", null);
    numpyEvaluationCheck(a => a[0, -1], x0[0, -1], X, X0, "x[0,-1]", null);
    numpyEvaluationCheck(a => a[.., 0], x0[.., 0], X, X0, "x[:,0]", null);
    numpyEvaluationCheck(a => a[.., 1], x0[.., 1], X, X0, "x[:,1]", null);
    numpyEvaluationCheck(a => a[1, ..], x0[1, ..], X, X0, "x[1,:]", null);
    numpyEvaluationCheck(a => a[0, ..], x0[0, ..], X, X0, "x[0,:]", null);
    numpyEvaluationCheck(a => a[-1, ..], x0[-1, ..], X, X0, "x[-1,:]", null);
    numpyEvaluationCheck(a => a[.., -2], x0[.., -2], X, X0, "x[:,-2]", null);
    numpyEvaluationCheck(a => a[0..^2, 0..^1], x0[0..^2, 0..^1], X, X0, "x[0:-2,0:-1]", null);
    numpyEvaluationCheck(a => a[0..2, 0..2], x0[0..2, 0..2], X, X0, "x[0:2,0:2]", null);
    numpyEvaluationCheck(a => a.get(idx(0, 1), ca.Slice(0, 2)), x0.get(idx(0, 1), ca.Slice(0, 2)), X, X0, "x[[0,1],0:2]", null);
    numpyEvaluationCheck(a => a.get_nz(false, idx(0, 2, 3)), ca.vertcat(0.738, 0.99, 0.2), X, X0, "x[[0,2,3]]", null);
    numpyEvaluationCheck(a => a.get_nz(false, idx(1)), ca.DM(0.1), X, X0, "x[1] on dense matrix", null);
    numpyEvaluationCheck(a => a.get_nz(false, idx(-1)), ca.DM(0.999999), X, X0, "x[-1] on dense matrix", null);
    numpyEvaluationCheck(a => a.get(idx(0, 1), ca.Slice(0, 1)), x0.get(idx(0, 1), ca.Slice(0, 1)), X, X0, "x[:,0:1]", null);
    numpyEvaluationCheck(a => a.get(ca.Slice(0, 1), idx(0, 1)), x0.get(ca.Slice(0, 1), idx(0, 1)), X, X0, "x[0:1,:]", null);

    message(":sparse");
    var sp = ca.Sparsity(4, 3, new long[] { 0, 2, 2, 3 }, new long[] { 1, 2, 1 });
    x = MX.sym("X", sp);
    X = new List<MX> { x };
    DM sx0 = ca.vertcat(0.738, 0.39, 0.99);
    x0 = ca.densify(dmsp(sp, 0.738, 0.39, 0.99));
    var S0 = new List<DM> { sx0 };
    numpyEvaluationCheck(a => a[0, 0], x0[0, 0], X, S0, "x[0,0]", null);
    numpyEvaluationCheck(a => a[1, 0], x0[1, 0], X, S0, "x[1,0]", null);
    numpyEvaluationCheck(a => a[0, 1], x0[0, 1], X, S0, "x[1,0]", null);
    numpyEvaluationCheck(a => a[0, -1], x0[0, -1], X, S0, "x[0,-1]", null);
    numpyEvaluationCheck(a => a[.., 0], x0[.., 0], X, S0, "x[:,0]", null);
    numpyEvaluationCheck(a => a[.., 1], x0[.., 1], X, S0, "x[:,1]", null);
    numpyEvaluationCheck(a => a[1, ..], x0[1, ..], X, S0, "x[1,:]", null);
    numpyEvaluationCheck(a => a[0, ..], x0[0, ..], X, S0, "x[0,:]", null);
    numpyEvaluationCheck(a => a[-1, ..], x0[-1, ..], X, S0, "x[-1,:]", null);
    numpyEvaluationCheck(a => a[.., -2], x0[.., -2], X, S0, "x[:,-2]", null);
    numpyEvaluationCheck(a => a[0..^2, 0..^1], x0[0..^2, 0..^1], X, S0, "x[0:-2,0:-1]", null);
    numpyEvaluationCheck(a => a[0..2, 0..2], x0[0..2, 0..2], X, S0, "x[0:2,0:2]", null);
    numpyEvaluationCheck(a => a.get(idx(0, 1), ca.Slice(0, 2)), x0.get(idx(0, 1), ca.Slice(0, 2)), X, S0, "x[[0,1],0:2]", null);
    numpyEvaluationCheck(a => a.get_nz(false, idx(2, 1)), ca.vertcat(x0[1, 2], x0[2, 0]), X, new List<DM> { x0 }, "x[[2,1]]", null);
    numpyEvaluationCheck(a => a.get_nz(false, ca.Slice(0, 2)), ca.vertcat(0.738, 0.39), X, new List<DM> { x0 }, "x[0:2] on dense matrix", null);
    numpyEvaluationCheck(a => a.get_nz(false, idx(1)), ca.DM(0.39), X, S0, "x[1]", null);
    numpyEvaluationCheck(a => a.get_nz(false, idx(-1)), ca.DM(0.99), X, S0, "x[-1]", null);
  }

  public void test_mx_in() {
    message("mx_out/mx_in");
    MX x = MX.sym("x", 2, 3);
    Function f = fun("f", new[] { x }, new[] { 3.0 * x });
    List<MX> xIn = f.mx_in();
    List<MX> xOut = f.call(xIn);
    Function g = fun("g", new[] { xIn[0] }, new[] { 6.0 * xOut[0] });
    double[] n = { 1, 2, 3, 4, 5, 6 };
    DM fOut = ev(f, ca.DM(f.sparsity_in(0), n));
    DM gOut = ev(g, ca.DM(g.sparsity_in(0), n));
    checkarray(gOut, 6.0 * fOut, "slicing(trans)");
  }

  public void test_scalarMX() {
    MX x = MX.sym("x");
    DM x0 = ca.DM(0.738);
    numpyEvaluationCheckPool(pool, x, x0, "scalarMX", null, null, false);
    // matrixpool is empty in mx.py
  }

  public void test_MXJacobian() {
    message("MX(1,1) unary operation, jacobian");
    Jpool = new List<Op1>();
    message("SX(1,1) unary operation, jacobian");
    MX x = MX.sym("x");
    DM x0 = new double[,] { { 0.738 } };
    numpyEvaluationCheckPool(Jpool, x, x0, "MX unary operations, jacobian", f => jacobianOld(f, 0, 0), null, false);
  }

  public void test_MXJacobians() {
    message("MX(3,1) unary operation, jacobian");
    MX x = MX.sym("x", 3, 1);
    DM x0 = ca.vertcat(0.738, 0.9, 0.3);
    numpyEvaluationCheckPool(Jpool, x, x0, "MX unary operations, jacobian", f => jacobianOld(f, 0, 0), null, false);
  }

  public void test_MXbinary() {
    message("MX binary operations");
    MX x = MX.sym("x", 3, 2);
    MX y = MX.sym("x", 3, 2);
    DM x0 = ca.DM(new double[,] { { 0.738, 0.2 }, { 0.1, 0.39 }, { 0.99, 0.999999 } });
    DM y0 = ca.DM(new double[,] { { 1.738, 0.6 }, { 0.7, 12 }, { 0, -6 } });
    numpyEvaluationCheckPool2(matrixbinarypool, x, y, x0, y0, "MX");
  }

  public void test_MXSparse() {
    message("MX unary operations, sparse");
    var sp = ca.Sparsity(4, 3, new long[] { 0, 2, 2, 3 }, new long[] { 1, 2, 1 });
    MX x = MX.sym("x", sp);
    DM x0 = dmsp(sp, 0.738, 0.1, 0.99);
    numpyEvaluationCheckPool(pool, x, x0, "MX", null, x0, true);
  }

  public void test_MXbinarySparse() {
    message("SX binary operations");
    var spx = ca.Sparsity(4, 3, new long[] { 0, 2, 2, 3 }, new long[] { 1, 2, 1 });
    var spy = ca.Sparsity(4, 3, new long[] { 0, 2, 2, 3 }, new long[] { 0, 2, 3 });
    MX xx = MX.sym("x", spx);
    MX yy = MX.sym("y", spy);
    DM x0 = dmsp(spx, 0.738, 0.1, 0.99);
    DM y0 = dmsp(spy, 1.738, 0.7, -6);
    numpyEvaluationCheckPool2(matrixbinarypool, xx, yy, x0, y0, "MX");
  }

  public void test_symbolcheck() {
    message("Check if non-symbolic inputs are caught");
    Check.Throws<SwigException>(() => fun("f", new[] { ca.MX(0) }, new[] { MX.sym("x") }));
  }

  public void test_unite() {
    message("unite operation");
    var xn = new double[3, 4];
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 4; j++) xn[i, j] = random.NextDouble();
    MX x = ca.MX(3, 4);
    MX y = MX.sym("x", 3, 4);
    MX z = ca.unite(x, y);
    Function f = fun("f", new[] { y }, new[] { z });
    checkarray(ev(f, xn), xn);

    var spx = ca.Sparsity(4, 3, new long[] { 0, 2, 2, 3 }, new long[] { 1, 2, 1 });
    var spy = ca.Sparsity(4, 3, new long[] { 0, 1, 2, 3 }, new long[] { 0, 2, 2 });
    DM nx = DM.zeros(spx);
    for (long k = 0; k < nx.nnz(); k++) nx.set_nz(random.NextDouble(), false, idx(k));
    DM ny = DM.zeros(spy);
    for (long k = 0; k < nx.nnz(); k++) ny.set_nz(random.NextDouble(), false, idx(k));
    double[,] nxn = nx.ToArray(), nyn = ny.ToArray();
    var sum = new double[4, 3];
    for (int i = 0; i < 4; i++)
      for (int j = 0; j < 3; j++) sum[i, j] = nxn[i, j] + nyn[i, j];
    MX xs = MX.sym("x", spx);
    MX ys = MX.sym("y", spy);
    z = ca.unite(xs, ys);
    f = fun("f", new[] { xs, ys }, new[] { z });
    checkarray(ev(f, nx, ny), sum);
  }

  public void test_imatrix_index() {
    message("IM indexing");
    MX X = MX.sym("x", 2, 2);
    MX Y = X.get_nz(false, (DM)new double[,] { { 0, 2 }, { 1, 1 }, { 3, 3 } });
    Function f = fun("f", new[] { X }, new[] { Y });
    DM fOut = ev(f, ca.DM(f.sparsity_in(0), ca.vertcat(1, 2, 3, 4)));
    checkarray(fOut, new double[,] { { 1, 3 }, { 2, 2 }, { 4, 4 } });

    MX Y2 = X[.., ..];
    Y2.set_nz(ca.DM(new double[,] { { 9, 8 } }), false, (DM)new double[,] { { 0, 2 } });
    f = fun("f", new[] { X }, new[] { Y2 });
    fOut = ev(f, ca.DM(f.sparsity_in(0), ca.vertcat(1, 2, 3, 4)));
    checkarray(fOut, new double[,] { { 9, 8 }, { 2, 4 } });
  }

  public void test_subsass() {
    message("Check subscripted assignment");
    MX X = MX.sym("x", 2, 2);
    X[0, 0] = ca.MX(5);
    X[0, 0] = 5.0;
    X[.., 0] = 8.0;

    MX x = MX.sym("X", 3, 4);
    var xn = new double[3, 4];
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 4; j++) xn[i, j] = random.NextDouble();
    var r = new double[7, 8];
    long[] cols = { 2, 4, 6, 7 };
    MX y = MX.zeros(7, 8);
    y.set(ca.Slice(1, 4), idx(2, 4, 6, 7), x);
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 4; j++) r[1 + i, cols[j]] = xn[i, j];
    Function fy = fun("fy", new[] { x }, new[] { y });
    checkarray(ev(fy, xn), r);

    y = ca.MX(7, 8);
    y.set(ca.Slice(1, 4), idx(2, 4, 6, 7), x);
    fy = fun("fy", new[] { x }, new[] { y });
    checkarray(ev(fy, xn), r);

    long[] kl = { 2, 4, 5, 8 };
    Sparsity s = y.sparsity();
    foreach (long k in kl) r[s.row()[(int)k], s.get_col()[(int)k]] = 1.0;
    y.set_nz(ca.MX(1), false, idx(2, 4, 5, 8));
    fy = fun("fy", new[] { x }, new[] { y });
    checkarray(ev(fy, xn), r);

    y.set_nz(x.get_nz(false, idx(0, 1, 2, 3)), false, idx(2, 4, 5, 8));
    s = y.sparsity();
    Sparsity sx = x.sparsity();
    int cnt = 0;
    foreach (long k in kl) {
      r[s.row()[(int)k], s.get_col()[(int)k]] = xn[sx.row()[cnt], sx.get_col()[cnt]];
      cnt++;
    }
    fy = fun("fy", new[] { x }, new[] { y });
    checkarray(ev(fy, xn), r);
  }

  public void test_erase() {
    message("Erase function");
    message(":dense");
    MX y = MX.sym("Y", 7, 8);
    MX z = y * 2.0;
    z.erase(new long[] { 1, 2, 3 }, new long[] { 2, 4, 6, 7 });
    Function f = fun("f", new[] { y }, new[] { z });
    double[] ones = Enumerable.Repeat(1.0, 56).ToArray();
    ca.DM(f.sparsity_in(0), ones);
    DM fOut = ev(f, 0);
    DM e = fOut;
    checkarray(fOut, e, "erase");
    message(":sparse");
  }

  static Dictionary<string, object> opts(params (string, object)[] kv) => kv.ToDictionary(p => p.Item1, p => p.Item2);

  void algebraCheck(DM A_, DM b_, DM C_, DM D_, DM e_, DM x_, MX A, MX b, MX C, MX D, MX e, MX x) {
    MX Axb = ca.mtimes(A, x) + b;
    MX Dxe = ca.mtimes(D, x) + e;
    MX a = ca.mtimes(ca.mtimes(Axb.T(), C), Dxe);
    Function f = fun("f", new[] { x, A, b, C, D, e }, new[] { a });
    DM fOut = ev(f, x_, A_, b_, C_, D_, e_);
    DM Axb_ = ca.mtimes(A_, x_) + b_;
    DM Dxe_ = ca.mtimes(D_, x_) + e_;
    DM f_ = ca.mtimes(ca.mtimes(Axb_.T(), C_), Dxe_);
    checkarray(fOut, f_, "evaluation");
    DM J_ = ca.mtimes(ca.mtimes(Dxe_.T(), C_.T()), A_) + ca.mtimes(ca.mtimes(Axb_.T(), C_), D_);
    for (int w = 0; w < 2; w++) {
      f = ca.Function("f", new[] { x, A, b, C, D, e }, new[] { a }, opts(("ad_weight", w), ("ad_weight_sp", w)));
      Function J = jacobianOld(f, 0, 0);
      List<DM> JOut = J.call(new[] { x_, A_, b_, C_, D_, e_ });
      checkarray(JOut[0], J_, "evaluation");
    }
  }

  public void test_MXalgebraDense() {
    message("Test some dense algebraic properties of matrices");
    int n = 3, m = 4;
    algebraCheck(randDM(m, n), randDM(m, 1), randDM(m, m), randDM(m, n), randDM(m, 1), randDM(n, 1),
                 MX.sym("A", m, n), MX.sym("b", m, 1), MX.sym("C", m, m), MX.sym("D", m, n), MX.sym("e", m, 1),
                 MX.sym("x", n, 1));
  }

  Sparsity randsparsity(int m, int n, bool skipMid) {
    var sp = ca.Sparsity(m, n);
    for (int k = 0; k < (n * m) / 2; k++) {
      int i = random.Next(m), j = random.Next(n);
      if (!skipMid || (i != m / 2 && (n == 1 || j != n / 2))) sp.add_nz(i, j);
    }
    return sp;
  }

  DM gentest(int m, int n, bool skipMid, out MX sym) {
    Sparsity As = randsparsity(m, n, skipMid);
    DM A_ = DM.zeros(As);
    for (long k = 0; k < As.nnz(); k++) A_.set_nz(random.NextDouble(), false, idx(k));
    sym = MX.sym("A", As);
    return A_;
  }

  void sparseAlgebra(int n, int m, bool skipMid) {
    DM A_ = gentest(m, n, skipMid, out MX A), b_ = gentest(m, 1, skipMid, out MX b), C_ = gentest(m, m, skipMid, out MX C),
       D_ = gentest(m, n, skipMid, out MX D), e_ = gentest(m, 1, skipMid, out MX e);
    var x0 = new double[n];
    for (int i = 0; i < n; i++) x0[i] = random.NextDouble();
    algebraCheck(A_, b_, C_, D_, e_, x0, A, b, C, D, e, MX.sym("x", n, 1));
  }

  public void test_MXalgebraSparse() {
    message("Test some sparse algebraic properties of matrices");
    sparseAlgebra(3, 4, false);
  }

  public void test_MXalgebraSparseSparse() {
    message("Test some sparse algebraic properties of matrices");
    sparseAlgebra(8, 10, true);
  }

  public void test_chaining() {
    message("Chaining SX and MX together");
    SX x = SX.sym("x");
    SX y = ca.pow(x, 3.0);
    Function f = ca.Function("f", new[] { x }, new[] { y });
    Function J = jacobianOld(f, 0, 0);
    MX X = MX.sym("X");
    Function F = fun("F", new[] { X }, J.call(new[] { X }));
    double x_ = 1.7;
    List<DM> FOut = evs(F, x_);
    checkarray(FOut[0], 3 * x_ * x_);
  }

  public void test_issue107() {
    message("Regression test for issue 107: +=");
    MX x = MX.sym("x");
    MX y = MX.sym("y");
    MX z = x;
    z += y;
    assertTrue(x.is_symbolic());
    assertFalse(z.is_symbolic());
  }

  static readonly DM[] none = new DM[0];

  public void test_MXd_trivial() {
    message("symbolic variables and constants jac");
    MX X = MX.sym("X", 10);
    MX V = MX.sym("V");
    MX J = ca.jacobian(X, X);
    assertEquals(10L, J.nnz());
    assertEquals(10L, J.size1());
    assertEquals(10L, J.size2());

    Function g = fun("g", new MX[0], new[] { J });
    checkarray(g.call(none)[0], DM.eye(10), "unit matrix");
    g = fun("g", new MX[0], new[] { ca.jacobian(MX.eye(3), X) });
    checkarray(g.call(none)[0], DM.zeros(9, 10), "zero matrix");
    g = fun("g", new MX[0], new[] { ca.jacobian(X, V) });
    checkarray(g.call(none)[0], DM.zeros(10, 1), "zero matrix");
    g = fun("g", new MX[0], new[] { ca.jacobian(MX.eye(3), V) });
    checkarray(g.call(none)[0], DM.zeros(9, 1), "zero matrix");
  }

  public void test_MXd_substractionl() {
    message("substraction jac");
    MX V = MX.sym("V");
    MX X = MX.sym("X");
    Function g = fun("g", new MX[0], new[] { ca.jacobian(X - V, X) });
    checkarray(g.call(none)[0], DM.ones(1, 1), "one");
    g = fun("g", new MX[0], new[] { ca.jacobian(X - V, V) });
    checkarray(g.call(none)[0], -DM.ones(1, 1), "one");
    g = fun("g", new MX[0], new[] { ca.jacobian(V - X, X) });
    checkarray(g.call(none)[0], -DM.ones(1, 1), "one");
    g = fun("g", new MX[0], new[] { ca.jacobian(V - X, V) });
    checkarray(g.call(none)[0], DM.ones(1, 1), "one");
  }

  public void test_MXd_mapping() {
    message("mapping jac");
    MX X = MX.sym("X", 3);
    MX Y = MX.sym("Y", 2);
    MX J = ca.jacobian(ca.vertcat(X, Y), X);
    DM JJ = DM.ones(J.sparsity());
    checkarray(JJ, ca.vertcat(DM.eye(3), DM.zeros(2, 3)), "diag");
    J = ca.jacobian(ca.vertcat(X, Y), Y);
    JJ = DM.ones(J.sparsity());
    checkarray(JJ, ca.vertcat(DM.zeros(3, 2), DM.eye(2)), "diag");
  }

  public void test_null() {
    message("Function null");
    MX x = MX.sym("x");
    Function f = fun("f", new[] { x }, new[] { ca.pow(x, 2.0), ca.MX() });
    assertEquals(0L, f.size1_out(1));
    assertEquals(0L, f.size2_out(1));
    evs(f, 0);

    f = fun("f", new[] { x, ca.MX() }, new[] { ca.pow(x, 2.0), ca.MX() });
    assertEquals(0L, f.size1_out(1));
    assertEquals(0L, f.size2_out(1));
    evs(f, 0, 0);

    List<MX> r = evs(f, x, ca.MX());
    assertTrue(r[1].is_empty(true));

    r = evs(f, ca.MX(), ca.MX());
    assertTrue(r[1].is_empty(true));
  }

  public void test_issue184() {
    message("Regression test issue #184");
    MX x = MX.sym("x", 3);
    MX y = x[0..0];
    assertEquals(0L, y.nnz());
  }

  public void test_indexingOutOfBounds() {
    message("Indexing out of bounds");
    DM y = DM.zeros(4, 5);
    Check.Throws<SwigException>(() => { var _ = y[12, 0]; });
    Check.Throws<SwigException>(() => { var _ = y[12, 12]; });
    Check.Throws<SwigException>(() => { var _ = y[0, 12]; });
    Check.Throws<SwigException>(() => { var _ = y[12, ..]; });
    Check.Throws<SwigException>(() => { var _ = y[12..15, 0]; });
    Check.Throws<SwigException>(() => { var _ = y[.., 12]; });
    Check.Throws<SwigException>(() => { var _ = y[0, 12..15]; });
    { var _ = y[-1, 2]; }
    Check.Throws<SwigException>(() => { var _ = y[-12, 2]; });
    { var _ = y[ca.Slice(-3, -1), 2]; }
    Check.Throws<SwigException>(() => { var _ = y[ca.Slice(-12, -9), 2]; });

    Check.Throws<SwigException>(() => y[12, 0] = 0.0);
    Check.Throws<SwigException>(() => y[12, 12] = 0.0);
    Check.Throws<SwigException>(() => y[0, 12] = 0.0);
    Check.Throws<SwigException>(() => y[12, ..] = 0.0);
    Check.Throws<SwigException>(() => y[12..15, 0] = 0.0);
    Check.Throws<SwigException>(() => y[.., 12] = 0.0);
    Check.Throws<SwigException>(() => y[0, 12..15] = 0.0);
    y[-1, 2] = 0.0;
    Check.Throws<SwigException>(() => y[-12, 2] = 0.0);
    y[ca.Slice(-3, -1), 2] = 0.0;
    Check.Throws<SwigException>(() => y[ca.Slice(-12, -9), 2] = 0.0);
  }

  public void test_indexinglimits() {
    message("Limits of indexing");
    MX y = MX.sym("y", 3);
    Check.Throws<SwigException>(() => y.get(idx(0, 5)));
    Check.Throws<SwigException>(() => y.set(idx(0, 5), MX.sym("a")));
    y.get(idx(0, 2));
    y.set(idx(0, 2), MX.sym("a"));
  }

  public void test_mtimes() {
    MX A = ca.MX(DM.ones((4, 3)));
    MX B = ca.MX(DM.ones((3, 8)));
    MX C = ca.MX(DM.ones((8, 7)));
    Check.Throws<SwigException>(() => ca.mtimes(new List<MX>()));
    MX D = ca.mtimes(new[] { A });
    assertEquals(4L, D.size1());
    assertEquals(3L, D.size2());
    D = ca.mtimes(new[] { A, B });
    assertEquals(4L, D.size1());
    assertEquals(8L, D.size2());
    D = ca.mtimes(new[] { A, B, C });
    assertEquals(4L, D.size1());
    assertEquals(7L, D.size2());
  }

  static Sparsity parseSp(string art) {
    var rows = art.Trim().Split('\n').Select(r => r.Trim()).Where(r => r.Length > 0).ToList();
    var a = new double[rows.Count, rows[0].Length];
    for (int i = 0; i < rows.Count; i++)
      for (int j = 0; j < rows[i].Length; j++) a[i, j] = rows[i][j] == 'X' ? 1 : 0;
    return ca.sparsify(a).sparsity();
  }

  static List<string> blasPlugins() {
    var plugins = new List<string> { "reference" };
    foreach (string opt in new[] { "classic", "blasfeo" }) {
      try {
        ca.load_blas(opt);
        plugins.Add(opt);
      } catch (SwigException) {
        // not available
      }
    }
    return plugins;
  }

  static List<string> classicLibs() => CasadiMeta.lapack_libraries().Split(';').Where(t => t.Length > 0).ToList();

  static bool compactibleMatch((bool, List<long>, List<long>) a, (bool, List<long>, List<long>) b,
                               (bool, List<long>, List<long>) z, Sparsity zsp) {
    return a.Item1 && b.Item1 && z.Item1 && zsp.nnz() > 0 && a.Item3.SequenceEqual(b.Item2) &&
           a.Item2.SequenceEqual(z.Item2) && b.Item3.SequenceEqual(z.Item3);
  }

  static string predictCategory(Sparsity A, Sparsity B, Sparsity Z) {
    if (A.is_dense() && B.is_dense() && Z.is_dense()) return "dense";
    if (A.is_dense() && Z.is_dense()) return "dense_sparse";
    if (compactibleMatch(A.is_compactible(), B.is_compactible(), Z.is_compactible(), Z)) return "dense";
    return "generic";
  }

  static string generateCode(Function f, string name, Dictionary<string, object> opts, string dir) {
    var cg = ca.CodeGenerator(name, opts);
    cg.add(f);
    cg.generate(dir + Path.DirectorySeparatorChar);
    return readFile(Path.Combine(dir, name + ".c"));
  }

  public void test_mtimes_blas_plugin() {
    string dir = tmpDir();
    DM.rng(1);
    int m = 3, k = 4, n = 5;
    string[] aLabels = { "dense", "compactible", "sparse" };
    Sparsity[] aPats = { parseSp("XXXX\nXXXX\nXXXX"), parseSp("X.X.\n....\nX.X."), parseSp("X.XX\n.X.X\nXX..") };
    string[] bLabels = { "dense", "compactible", "sparse" };
    Sparsity[] bPats = { parseSp("XXXXX\nXXXXX\nXXXXX\nXXXXX"), parseSp(".X.X.\n.....\n.X.X.\n....."),
                         parseSp("X.X.X\n.X..X\nXX.X.\n.X.XX") };
    DM AdenseVal = DM.rand(m, k);
    DM BdenseVal = DM.rand(k, n);
    List<string> plugins = blasPlugins();
    var denseKernel = new Dictionary<string, string> { ["reference"] = "casadi_mtimes_dense(",
      ["classic"] = "CASADI_BLAS_DGEMM(", ["blasfeo"] = "blasfeo_blas_dgemm(" };
    var pluginExtralibs = new Dictionary<string, List<string>> { ["reference"] = new List<string>(),
      ["classic"] = classicLibs(), ["blasfeo"] = new List<string> { "blasfeo" } };
    MX AdenseSym = MX.sym("A", m, k);
    MX BdenseSym = MX.sym("B", k, n);
    MX ZdenseSym = MX.sym("Z", m, n);
    DM ZdenseVal = randDM(m, n);

    foreach (string name in plugins) {
      for (int ia = 0; ia < 3; ia++) {
        for (int ib = 0; ib < 3; ib++) {
          foreach (string opLabel in new[] { "mtimes", "mac" }) {
            Sparsity aSp = aPats[ia], bSp = bPats[ib];
            string tag = opLabel + "/" + name + "/" + aLabels[ia] + "/" + bLabels[ib];
            MX A = ca.project(AdenseSym, aSp);
            MX B = ca.project(BdenseSym, bSp);
            DM Amask = ca.densify(ca.DM(aSp, 1.0));
            DM Bmask = ca.densify(ca.DM(bSp, 1.0));
            MX outF, outRef;
            MX[] symIn;
            List<DM> valIn;
            Sparsity zSp;
            if (opLabel == "mtimes") {
              outF = ca.densify(ca.mtimes(A, B, name));
              outRef = ca.mtimes(AdenseSym * Amask, BdenseSym * Bmask, "reference");
              symIn = new[] { AdenseSym, BdenseSym };
              valIn = dms(AdenseVal, BdenseVal);
              zSp = ca.mtimes(ca.DM(aSp, 1.0), ca.DM(bSp, 1.0)).sparsity();
            } else {
              outF = ca.mac(A, B, ZdenseSym, name);
              outRef = ca.mac(AdenseSym * Amask, BdenseSym * Bmask, ZdenseSym, "reference");
              symIn = new[] { AdenseSym, BdenseSym, ZdenseSym };
              valIn = dms(AdenseVal, BdenseVal, ZdenseVal);
              zSp = Sparsity.dense(m, n);
            }
            Function f = fun("f", symIn, new[] { outF });
            Function fRef = fun("f_ref", symIn, new[] { outRef });
            checkfunction(f, fRef, valIn);

            DM Adm = ca.project(AdenseVal, aSp);
            DM Bdm = ca.project(BdenseVal, bSp);
            DM dmOut = opLabel == "mtimes" ? ca.mtimes(Adm, Bdm, name) : ca.mac(Adm, Bdm, ZdenseVal, name);
            checkarray(ca.densify(dmOut), fRef.call(valIn)[0], tag + " (DM-level)");

            string cat = predictCategory(aSp, bSp, zSp);
            string expected = cat == "dense" ? denseKernel[name]
                            : cat == "dense_sparse" ? "casadi_mtimes_dense_sparse(" : "casadi_mtimes(";

            string code = generateCode(f, "f", new Dictionary<string, object>(), dir);
            assertTrue(code.Contains(expected), tag + ": expected kernel " + expected + " in generated code");

            string fpath = Path.Combine(dir, "f.casadi");
            f.save(fpath);
            Function fdeser = Function.load(fpath);
            code = generateCode(fdeser, "fdeser", new Dictionary<string, object>(), dir);
            assertTrue(code.Contains(expected), tag + ": expected kernel " + expected + " in deserialized code");

            checkfunctionLight(f, fdeser, valIn);

            List<string> extralibs = pluginExtralibs.GetValueOrDefault(name, new List<string>());
            if (name == "reference" || extralibs.Count > 0) {
              if (expected == denseKernel[name] || name == "reference")
                checkCodegen(f, valIn, new Dictionary<string, object>(), extralibs);
            }
          }
        }
      }
    }
  }

  public void test_l1_blas_plugin() {
    string dir = tmpDir();
    DM.rng(2);
    List<string> plugins = blasPlugins();
    var pluginExtralibs = new Dictionary<string, List<string>> { ["reference"] = new List<string>(),
      ["classic"] = classicLibs(), ["blasfeo"] = new List<string> { "blasfeo" } };
    var expectedKernel = new Dictionary<string, string> {
      ["reference/dot"] = "casadi_real casadi_dot(",
      ["reference/norm_1"] = "casadi_real casadi_norm_1(",
      ["reference/norm_2"] = "casadi_real casadi_norm_2(",
      ["reference/axpy"] = "void casadi_axpy(",
      ["reference/copy"] = "memcpy(y, x",
      ["classic/dot"] = "CASADI_BLAS_DDOT",
      ["classic/norm_1"] = "CASADI_BLAS_DASUM",
      ["classic/norm_2"] = "CASADI_BLAS_DNRM2",
      ["classic/axpy"] = "CASADI_BLAS_DAXPY",
      ["classic/copy"] = "memcpy(y, x",
      ["blasfeo/dot"] = "blasfeo_blas_ddot",
      ["blasfeo/norm_1"] = "casadi_real casadi_norm_1(",
      ["blasfeo/norm_2"] = "casadi_real casadi_norm_2(",
      ["blasfeo/axpy"] = "blasfeo_blas_daxpy",
      ["blasfeo/copy"] = "memcpy(y, x",
    };

    int n = 7;
    MX xSym = MX.sym("x", n);
    MX ySym = MX.sym("y", n);
    DM xVal = DM.rand(n);
    DM yVal = DM.rand(n);
    GlobalOptions.setDefaultBlas("reference");
    var opNames = new List<string> { "dot", "norm_1", "norm_2" };
    Function fRef = ca.Function("f_ref", new[] { xSym, ySym },
      new[] { ca.dot(xSym, ySym), ca.norm_1(xSym), ca.norm_2(xSym) }, new[] { "x", "y" }, opNames);
    List<DM> refOut = fRef.call(new[] { xVal, yVal });

    MX innerX = MX.sym("x", n);
    Function inner = fun("inner", new[] { innerX }, new[] { 2.0 * innerX });
    Function fAxpyRef = inner.map(4, new[] { false }, new[] { true });
    List<DM> fAxpyInputs = dms(ca.horzcat(xVal, xVal, xVal, xVal));
    List<DM> fAxpyRefOut = fAxpyRef.call(fAxpyInputs);

    try {
      foreach (string name in plugins) {
        GlobalOptions.setDefaultBlas(name);
        string tagBase = "blas=" + name;
        Function f = ca.Function("f", new[] { xSym, ySym },
          new[] { ca.dot(xSym, ySym), ca.norm_1(xSym), ca.norm_2(xSym) }, new[] { "x", "y" }, opNames);
        List<DM> got = f.call(new[] { xVal, yVal });
        for (int i = 0; i < opNames.Count; i++)
          checkarray(got[i], refOut[i], tagBase + " op=" + opNames[i] + " numerical agreement");

        Function fAxpy = inner.map(4, new[] { false }, new[] { true });
        List<DM> gotAxpy = fAxpy.call(fAxpyInputs);
        checkarray(gotAxpy[0], fAxpyRefOut[0], tagBase + " op=axpy numerical agreement");

        var perOpCarriers = new List<(string, Function)> {
          ("dot", fun("g_dot", new[] { xSym, ySym }, new[] { ca.dot(xSym, ySym) })),
          ("norm_1", fun("g_norm_1", new[] { xSym }, new[] { ca.norm_1(xSym) })),
          ("norm_2", fun("g_norm_2", new[] { xSym }, new[] { ca.norm_2(xSym) })),
          ("axpy", inner.map(4, new[] { false }, new[] { true })),
          ("copy", fun("g_copy", new[] { xSym }, new[] { xSym })),
        };
        foreach (var (key, fun_) in perOpCarriers) {
          string code = generateCode(fun_, "l1", opts(("l1_blas", true)), dir);
          string exp = expectedKernel[name + "/" + key];
          assertTrue(code.Contains(exp), tagBase + " op=" + key + ": expected " + exp + " in generated code");
        }

        string fpath = Path.Combine(dir, "l1.casadi");
        f.save(fpath);
        Function fdeser = Function.load(fpath);
        checkfunctionLight(f, fdeser, dms(xVal, yVal));

        List<string> extralibs = pluginExtralibs.GetValueOrDefault(name, new List<string>());
        if (name == "reference" || extralibs.Count > 0) {
          checkCodegen(f, dms(xVal, yVal), opts(("l1_blas", true)), extralibs);
          checkCodegen(fAxpy, fAxpyInputs, opts(("l1_blas", true)), extralibs);
        }
      }
    } finally {
      GlobalOptions.setDefaultBlas("reference");
    }
  }

  public void test_truth() {
    message("Truth values");
    Check.Throws<SwigException>(() => MX.sym("x").nonzero());
    assertTrue(ca.MX(1).nonzero());
    assertFalse(ca.MX(0).nonzero());
    assertTrue(ca.MX(0.2).nonzero());
    assertTrue(ca.MX(-0.2).nonzero());
    Check.Throws<SwigException>(() => ca.MX(ca.vertcat(2.0, 3)).nonzero());
    Check.Throws<SwigException>(() => ca.MX().nonzero());
  }

  public void test_MXbool() {
    message("bool");
    MX xy = MX.sym("x", 2);
    MX x = xy[0];
    MX y = xy[1];
    Function f = fun("f", new[] { xy }, new[] { ca.vertcat(ca.logic_and(x, y), ca.logic_or(x, y), ca.logic_not(x)) });
    for (int t1 = 0; t1 < 2; t1++) {
      for (int t2 = 0; t2 < 2; t2++) {
        bool T1 = t1 != 0, T2 = t2 != 0;
        DM fOut = ev(f, new double[] { t1, t2 });
        checkarray(fOut, ca.vertcat(T1 && T2 ? 1 : 0, T1 || T2 ? 1 : 0, !T1 ? 1 : 0), "bool(" + t1 + "," + t2 + "): " + fOut);
      }
    }
  }

  public void test_short_circuiting_codegen() {
    MX x = MX.sym("x", 5);
    MX y = ca.if_else(x > 0.0, 0.0, x);
    Function f = ca.Function("f", new[] { x }, new[] { y }, new[] { "x" }, new[] { "y" });
    checkCodegen(f, dms(ca.vertcat(-200, 200, -100, 100, 0)));

    MX x8 = MX.sym("x", 8);
    MX y8 = MX.sym("y", 8);
    MX z = ca.logic_and(x8, y8);
    MX z2 = ca.logic_or(x8, y8);
    f = fun("f", new[] { x8, y8 }, new[] { z, z2 });
    checkCodegen(f, dms(ca.vertcat(0, 1, 0, 1, 0, 1, 0, 1), ca.vertcat(0, 0, 1, 1, 0, 0, 1, 1)));
  }

  public void test_MXineq() {
    message("SX ineq");
    MX xy = MX.sym("x", 2);
    MX x = xy[0];
    MX y = xy[1];
    Function f = fun("f", new[] { xy }, new[] { ca.vertcat(x < y, x <= y, x >= y, ca.eq(x, y), ca.ne(x, y)) });
    double[] vals = { -10, 0.1, 0, 1, 10 };
    foreach (double t1 in vals) {
      foreach (double t2 in vals) {
        DM fOut = ev(f, new double[] { t1, t2 });
        checkarray(fOut, ca.vertcat(t1 < t2 ? 1 : 0, t1 <= t2 ? 1 : 0, t1 >= t2 ? 1 : 0, t1 == t2 ? 1 : 0, t1 != t2 ? 1 : 0),
                   "ineq(" + t1 + "," + t2 + ")");
      }
    }
  }

  public void test_if_else_zero() {
    MX x = MX.sym("x");
    MX y = ca.if_else(x, 5.0, 0.0);
    Function f = fun("f", new[] { x }, new[] { y });
    DM fOut = ev(f, 1);
    assertTrue((double)fOut == 5, "if_else_zero " + fOut);
    fOut = ev(f, 0);
    assertTrue((double)fOut == 0, "if_else_zero");
  }

  public void test_if_else() {
    MX x = MX.sym("x");
    MX y = ca.if_else(x, 1.0, 2.0);
    Function f = fun("f", new[] { x }, new[] { y });
    assertTrue((double)ev(f, 1) == 1, "if_else");
    assertTrue((double)ev(f, 0) == 2, "if_else");
  }

  static bool allFinite(DM r) => r.ToArray().Cast<double>().All(double.IsFinite);

  static string arr(double[] a) => "[" + string.Join(", ", a) + "]";

  public void test_if_else_call_reverse_4387() {
    double[][] masks = { new double[] { 0, 0 }, new double[] { 1, 0 }, new double[] { 0, 1 }, new double[] { 1, 1 } };
    foreach (double[] maskValue in masks) {
      MX x = MX.sym("x", 3);
      MX mask = MX.sym("mask", 2);
      Function f = ca.Function("f", new[] { x }, new[] { ca.vertcat(ca.sumsqr(x), ca.sum1(ca.sin(x))) }, opts(("never_inline", true)));
      MX z = ca.if_else(mask > 0.0, f.Call(x), 0.0);
      Function g = fun("g", new[] { x, mask }, new[] { ca.gradient(ca.sumsqr(z), x) });
      double[] value = { 1, 2, 3 };
      double vv = 0, ss = 0;
      foreach (double v in value) {
        vv += v * v;
        ss += Math.Sin(v);
      }
      var expected = new double[3];
      for (int i = 0; i < 3; i++)
        expected[i] = maskValue[0] * 4 * vv * value[i] + maskValue[1] * 2 * ss * Math.Cos(value[i]);
      checkarray(ev(g, value, maskValue), ca.DM(expected), "rectangular nonlinear " + arr(maskValue));
    }

    double[][] nonbooleanMasks = { new double[] { 0, 0 }, new double[] { 1, 0 }, new double[] { 0, 1 }, new double[] { 1, 1 },
                                   new double[] { -1, 0 }, new double[] { 0, -2 }, new double[] { 2, -3 } };
    string[] names = { "square", "scalar input", "sparse mask" };
    DM[] @as = { ca.DM(new double[,] { { 1, 1 }, { 1, -1 } }), ca.vertcat(1, 2), ca.DM(new double[,] { { 1, 2, 3 }, { 4, 5, 6 }, { 7, 8, 9 } }) };
    Sparsity[] sps = { Sparsity.dense(2, 1), Sparsity.dense(2, 1), Sparsity.triplet(3, 1, new long[] { 0, 2 }, new long[] { 0, 0 }) };
    double[][][] maskValuess = { nonbooleanMasks, nonbooleanMasks,
      new[] { new double[] { 0, 0, 0 }, new double[] { 1, 0, 0 }, new double[] { 0, 0, 1 }, new double[] { 1, 0, 1 } } };
    for (int c = 0; c < 3; c++) {
      DM a = @as[c];
      Sparsity sp = sps[c];
      foreach (double[] maskValue in maskValuess[c]) {
        for (int ndir = 1; ndir <= 2; ndir++) {
          string sub = names[c] + " mask=" + arr(maskValue) + " ndir=" + ndir;
          MX x = MX.sym("x", a.size2());
          MX mask = MX.sym("mask", sp);
          Function f = ca.Function("f", new[] { x }, new[] { ca.mtimes(a, x) }, opts(("never_inline", true)));
          MX z = ca.if_else(mask, f.Call(x), 0.0);
          var r2 = new double[a.size1()];
          for (int i = 0; i < r2.Length; i++) r2[i] = i + 2;
          List<DM> seeds = new List<DM> { DM.ones(a.size1()), ca.DM(r2) }.GetRange(0, ndir);
          var v = seeds.Select(seed => new[] { ca.MX(seed) }).ToList();
          List<List<MX>> adj = ca.reverse(new[] { z }, new[] { x }, v);
          Function g = fun("g", new[] { x, mask }, new[] { ca.horzcat(adj.Select(d => d[0])) });
          var exp = seeds.Select(seed => ca.mtimes(a.T(), ca.ne(ca.DM(maskValue), 0.0) * seed)).ToList();
          checkarray(ev(g, DM.ones(a.size2()), ca.project(ca.DM(maskValue), sp)), ca.horzcat(exp), sub);
        }
      }
    }

    var cases = new (bool, double[], double[])[] {
      (true, new double[] { 0 }, new double[] { 0, 0, 0 }), (true, new double[] { 1 }, new double[] { 0.5, 1, 1 }),
      (false, new double[] { 0, 0 }, new double[] { 0, 0, 0 }), (false, new double[] { 1, 0 }, new double[] { 0.5, 0, 0 }) };
    foreach (var (scalarMask, maskValue, expected) in cases) {
      MX x = MX.sym("x", 3);
      MX mask = MX.sym("mask", scalarMask ? 1 : 2);
      Function f = ca.Function("f", new[] { x }, new[] { ca.vertcat(ca.sqrt(x[0]), x[1] + x[2]) }, opts(("never_inline", true)));
      MX z = ca.if_else(mask, f.Call(x), 0.0);
      Function g = fun("g", new[] { x, mask }, new[] { ca.gradient(ca.sum1(z), x) });
      bool active = maskValue[0] != 0;
      DM result = ev(g, new double[] { active ? 1 : 0, 2, 3 }, maskValue);
      assertTrue(allFinite(result));
      checkarray(result, ca.DM(expected), "singular derivative " + scalarMask + " " + arr(maskValue));
    }
  }

  public void test_if_else_reverse_partial_nan() {
    Check.Assume(KNOWN_BUGS, "known bug");
    foreach (bool inline in new[] { false, true }) {
      MX x = MX.sym("x", 3);
      MX mask = MX.sym("mask", 2);
      Function f = ca.Function("f", new[] { x }, new[] { ca.vertcat(ca.sqrt(x[0]), x[1] + x[2]) },
                                opts(("always_inline", inline), ("never_inline", !inline)));
      MX z = ca.if_else(mask, f.Call(x), 0.0);
      Function g = fun("g", new[] { x, mask }, new[] { ca.gradient(ca.sum1(z), x) });
      DM result = ev(g, new double[] { 0, 2, 3 }, new double[] { 0, 1 });
      assertTrue(allFinite(result));
      checkarray(result, ca.vertcat(0, 1, 1));
    }
  }

  public void test_regression491() {
    message("regression #491");
    SX u = SX.sym("u");
    SX x = SX.sym("x");
    Function F = ca.Function("F", new[] { u, x }, new[] { u + 1.0 / x });
    MX U = MX.sym("U");
    MX X = F.Call(U, U);
    MX G = F.Call(U, X);
    for (int kk = 0; kk < 2; kk++) {
      Function gfcn;
      if (kk == 0) gfcn = fun("gfcn", new[] { U }, new[] { G }).expand("e_gfcn", opts(("ad_weight", 1)));
      else gfcn = ca.Function("gfcn", new[] { U }, new[] { G }, opts(("ad_weight", 1)));
      Function J = jacobianOld(gfcn, 0, 0);
      List<DM> JOut = J.call(dms(1));
      checkarray((double)JOut[0], 1, 1e-9);
    }
  }

  // test_ticket: python-only (Python list + MX operator)

  public void test_jacobian_tools() {
    message("jacobian");
    MX X = MX.sym("X");
    MX Y = ca.jacobian(ca.pow(X, 2.0), X);
    Function f = fun("f", new[] { X }, new[] { Y });
    checkarray((double)ev(f, 2.3), 4.6, 1e-7);
  }

  public void test_reshape() {
    message("reshape");
    MX X = MX.sym("X", 10);
    DM i = ca.DM(Sparsity.lower(3), ca.vertcat(range(6)));
    i.print_dense();
    Console.WriteLine(i.T().get_nz(false, ca.Slice()));
    MX T = X.get_nz(false, i);
    MX q = ca.pow(T.T().get_nz(false, ca.Slice()), 2.0);
    Function f = fun("f", new[] { X }, new[] { q });
    DM fOut = ev(f, range(10));
    checkarray(ca.vertcat(0, 1, 9, 4, 16, 25), fOut);

    MX Y = MX.sym("Y", 10);
    Function ff = fun("ff", new[] { Y }, f.call(new[] { Y }, true));
    checkarray(ca.vertcat(0, 1, 9, 4, 16, 25), ev(ff, range(10)));

    Function J = fun("J", new[] { X }, new[] { ca.jacobian(q, X) });
    DM JOut = ev(J, range(10));
    DM i2 = ca.horzcat(ca.diag(ca.vertcat(0, 2, 4, 6, 8, 10)), DM.zeros(6, 4));
    i2.set(idx(2, 3), ca.Slice(), i2.get(idx(3, 2), ca.Slice()));
    checkarray(i2, JOut);

    q = ca.pow(T.T().get_nz(false, ca.Slice()), 2.0);
    J = fun("J", new[] { X }, new[] { ca.jacobian(q, X) });
    JOut = ev(J, range(10));
    i2 = ca.horzcat(ca.diag(ca.vertcat(0, 2, 4, 6, 8, 10)), DM.zeros(6, 4));
    i2.set(idx(2, 3), ca.Slice(), i2.get(idx(3, 2), ca.Slice()));
    checkarray(i2, JOut);
  }

  public void test_vertcat() {
    message("vertcat");
    MX X = MX.sym("X", 10);
    MX T = ca.vertcat(X[4], X[2]);
    MX q = ca.pow(T, 2.0);
    Function f = fun("f", new[] { X }, new[] { q });
    List<DM> fOut = f.call(dms(ca.vertcat(range(10))));
    checkarray(ca.vertcat(16, 4), fOut[0]);

    MX Y = MX.sym("Y", 10);
    Function ff = fun("ff", new[] { Y }, f.call(new[] { Y }, true));
    checkarray(ca.vertcat(16, 4), ev(ff, range(10)));
    Function J = fun("J", new[] { X }, new[] { ca.jacobian(q, X) });
    DM JOut = ev(J, range(10));
    DM i = DM.zeros(2, 10);
    i[0, 4] = 8.0;
    i[1, 2] = 4.0;
    checkarray(i, JOut);
    q = ca.pow(T, 2.0);
    J = fun("J", new[] { X }, new[] { ca.jacobian(q, X) });
    JOut = ev(J, range(10));
    checkarray(i, JOut);
  }

  public void test_blockcat() {
    MX x = MX.sym("x");
    MX y = ca.blockcat(new[] { new[] { x, 2.0 * x }, new[] { 3.0 * x, 4.0 * x } });
    Function f = fun("f", new[] { x }, new[] { y });
    checkarray(ev(f, 3), ca.DM(new double[,] { { 3, 6 }, { 9, 12 } }));
  }

  public void test_veccats() {
    MX x = MX.sym("x", 2);
    assertEquals(x.hash(), ca.vec(x).hash());
  }

  public void test_constmxmtimes() {
    var _ = 0.1 * MX.ones(2);
  }

  public void test_is_regular() {
    assertTrue(ca.MX(ca.vertcat(0, 1)).is_regular());
    assertFalse(ca.MX(ca.vertcat(0, double.PositiveInfinity)).is_regular());
    Check.Throws<SwigException>(() => MX.sym("x", 2).is_regular());
  }

  public void test_diagcat() {
    MX C = ca.diagcat(ca.MX(ca.DM(new double[,] { { -1.4, -3.2 }, { -3.2, -28 } })),
                   ca.MX(ca.DM(new double[,] { { 15, -12, 2.1 }, { -12, 16, -3.8 }, { 2.1, -3.8, 15 } })), ca.MX(1.8), ca.MX(-4.0));
    DM r = ca.DM(new double[,] { { -1.4, -3.2, 0, 0, 0, 0, 0 }, { -3.2, -28, 0, 0, 0, 0, 0 }, { 0, 0, 15, -12, 2.1, 0, 0 },
                              { 0, 0, -12, 16, -3.8, 0, 0 }, { 0, 0, 2.1, -3.8, 15, 0, 0 }, { 0, 0, 0, 0, 0, 1.8, 0 }, { 0, 0, 0, 0, 0, 0, -4 } });
    r = ca.sparsify(r);
    Function f = fun("f", new MX[0], new[] { C });
    checkarray(f.call(none)[0], r);

    int cx = 3, cy = 4;
    for (int nx = 2; nx >= 0; nx--) {
      for (int ny = 2; ny >= 0; ny--) {
        Console.WriteLine(nx + " " + ny);
        DM rd = ca.diagcat(ca.DM(nx, ny), DM.ones(cx, cy), ca.DM(nx, ny));
        assertEquals((long)(cx + 2 * nx), rd.size1());
        assertEquals((long)(cy + 2 * ny), rd.size2());
        SX rs = ca.diagcat(ca.SX(nx, ny), SX.ones(cx, cy), ca.SX(nx, ny));
        assertEquals((long)(cx + 2 * nx), rs.size1());
        assertEquals((long)(cy + 2 * ny), rs.size2());
        MX rm = ca.diagcat(ca.MX(nx, ny), MX.ones(cx, cy), ca.MX(nx, ny));
        assertEquals((long)(cx + 2 * nx), rm.size1());
        assertEquals((long)(cy + 2 * ny), rm.size2());
      }
    }
  }

  public void test_tril2symm() {
    MX x = MX.sym("x", Sparsity.lower(3));
    Function f = fun("f", new[] { x }, new[] { ca.tril2symm(x) });
    DM fOut = ev(f, ca.DM(f.sparsity_in(0), ca.vertcat(range(6))));
    checkarray(fOut, ca.DM(new double[,] { { 0, 1, 2 }, { 1, 3, 4 }, { 2, 4, 5 } }));
  }

  public void test_sparsity_indexing() {
    message("sparsity");
    DM B_ = ca.DM(new double[,] { { 1, 2, 3, 4, 5 }, { 6, 7, 8, 9, 10 } });
    MX B = MX.sym("B", 2, 5);
    DM A = ca.sparsify(ca.DM(new double[,] { { 1, 1, 0, 0, 0 }, { 0, 0, 1, 0, 0 } }));
    Sparsity sp = A.sparsity();
    Func<MX, DM> meval = m => ev(fun("f", new[] { B }, new[] { m }), B_);

    checkarray(meval(B.get(sp)), ca.DM(new double[,] { { 1, 2, 0, 0, 0 }, { 0, 0, 8, 0, 0 } }), "sparsity indexing");

    MX Bmod = ca.MX(B);
    Bmod.set(sp, -4.0);
    checkarray(meval(Bmod), ca.DM(new double[,] { { -4, -4, 3, 4, 5 }, { 6, 7, -4, 9, 10 } }), "sparsity indexing assignement");

    Bmod = ca.MX(B);
    Bmod.set(sp, 2.0 * B);
    checkarray(meval(Bmod), ca.DM(new double[,] { { 2, 4, 3, 4, 5 }, { 6, 7, 16, 9, 10 } }), "Imatrix indexing assignement");

    Check.Throws<SwigException>(() => B.get(Sparsity.dense(4, 4)));
  }

  public void test_symvar() {
    MX a = MX.sym("a");
    MX b = MX.sym("b");
    MX c = MX.sym("c");
    MX e = ca.cos(a * b) + c;
    List<MX> w = ca.symvar(e);
    assertEquals(3, w.Count);
    if (GlobalOptions.getSimplificationOnTheFly()) {
      assertTrue(ca.is_equal(w[0], a));
      assertTrue(ca.is_equal(w[1], b));
      assertTrue(ca.is_equal(w[2], c));
    }
  }

  // Elementwise algebra on SX, MX or DM through dynamic dispatch, so one expression list serves all types
  delegate dynamic Gen(dynamic x, dynamic y, dynamic z, dynamic X, dynamic Y);

  sealed class Simp {
    public readonly bool refcount, onTheFly;
    public readonly Gen genA, genB;
    public Simp(bool refcount, bool onTheFly, Gen genA, Gen genB) {
      this.refcount = refcount;
      this.onTheFly = onTheFly;
      this.genA = genA;
      this.genB = genB;
    }
  }

  static Simp s(bool refcount, bool onTheFly, Gen genA, Gen genB) => new Simp(refcount, onTheFly, genA, genB);

  /// <summary>Python str() of the expected value; a string stands for a Python number literal</summary>
  static string pyStr(object o) => o.ToString();

  static List<Simp> simplificationCases(dynamic dx) {
    double inf = double.PositiveInfinity;
    // "X is ca.MX" is always False in mx.py: X is rebound to sin(A) before the list is built
    bool isMX = false;
    return new List<Simp> {
      s(false, true, (x, y, z, X, Y) => y * x + y * (1 - x), (x, y, z, X, Y) => y),
      s(false, isMX, (x, y, z, X, Y) => 2 * x - x, (x, y, z, X, Y) => x),
      s(false, true, (x, y, z, X, Y) => 4 * x - 3 * x, (x, y, z, X, Y) => x),
      s(false, true, (x, y, z, X, Y) => 0.2 * x + 0.8 * x, (x, y, z, X, Y) => x),
      s(false, true, (x, y, z, X, Y) => x + 0, (x, y, z, X, Y) => x),
      s(false, isMX, (x, y, z, X, Y) => 0 + x, (x, y, z, X, Y) => x),
      s(false, true, (x, y, z, X, Y) => x + (-y), (x, y, z, X, Y) => x - y),
      s(false, isMX, (x, y, z, X, Y) => (-x) + y, (x, y, z, X, Y) => y - x),
      s(false, true, (x, y, z, X, Y) => 0.5 * x + 0.5 * x, (x, y, z, X, Y) => x),
      s(false, true, (x, y, z, X, Y) => x / 2 + x / 2, (x, y, z, X, Y) => x),
      s(false, true, (x, y, z, X, Y) => (x - y) + y, (x, y, z, X, Y) => x),
      s(false, true, (x, y, z, X, Y) => y + (x - y), (x, y, z, X, Y) => x),
      s(false, true, (x, y, z, X, Y) => ca.pow(ca.sin(x), 2) + ca.pow(ca.cos(x), 2), (x, y, z, X, Y) => "1"),
      s(false, true, (x, y, z, X, Y) => x - 0, (x, y, z, X, Y) => x),
      s(false, isMX, (x, y, z, X, Y) => 0 - x, (x, y, z, X, Y) => -x),
      s(false, true, (x, y, z, X, Y) => x - x, (x, y, z, X, Y) => "0"),
      s(false, true, (x, y, z, X, Y) => x - (-y), (x, y, z, X, Y) => x + y),
      s(false, true, (x, y, z, X, Y) => (x + y) - y, (x, y, z, X, Y) => x),
      s(false, true, (x, y, z, X, Y) => (x + y) - x, (x, y, z, X, Y) => y),
      s(false, true, (x, y, z, X, Y) => x - (y + x), (x, y, z, X, Y) => -y),
      s(false, true, (x, y, z, X, Y) => x - (x + y), (x, y, z, X, Y) => -y),
      s(false, true, (x, y, z, X, Y) => (-x) - y, (x, y, z, X, Y) => -(x + y)),
      s(false, true, (x, y, z, X, Y) => x * x, (x, y, z, X, Y) => ca.pow(x, 2)),
      s(false, true, (x, y, z, X, Y) => x * 3, (x, y, z, X, Y) => 3 * x),
      s(false, true, (x, y, z, X, Y) => x * 0, (x, y, z, X, Y) => "0"),
      s(false, true, (x, y, z, X, Y) => 0 * x, (x, y, z, X, Y) => "0"),
      s(false, true, (x, y, z, X, Y) => x * 1, (x, y, z, X, Y) => x),
      s(false, isMX, (x, y, z, X, Y) => 1 * x, (x, y, z, X, Y) => x),
      s(false, true, (x, y, z, X, Y) => x * (-1), (x, y, z, X, Y) => -x),
      s(false, isMX, (x, y, z, X, Y) => (-1) * x, (x, y, z, X, Y) => -x),
      s(false, true, (x, y, z, X, Y) => x * (1 / y), (x, y, z, X, Y) => x / y),
      s(false, isMX, (x, y, z, X, Y) => (1 / x) * y, (x, y, z, X, Y) => y / x),
      s(false, true, (x, y, z, X, Y) => 5 * (0.2 * x), (x, y, z, X, Y) => x),
      s(false, true, (x, y, z, X, Y) => 0.5 * (2 * x), (x, y, z, X, Y) => x),
      s(false, true, (x, y, z, X, Y) => 5 * (x / 5), (x, y, z, X, Y) => x),
      s(false, true, (x, y, z, X, Y) => (2 / x) * x, (x, y, z, X, Y) => "2"),
      s(false, true, (x, y, z, X, Y) => (2 * x) / x, (x, y, z, X, Y) => "2"),
      s(false, true, (x, y, z, X, Y) => (5 * x) / x, (x, y, z, X, Y) => "5"),
      s(false, true, (x, y, z, X, Y) => x * (2 / x), (x, y, z, X, Y) => "2"),
      s(false, isMX, (x, y, z, X, Y) => (-x) * y, (x, y, z, X, Y) => -(x * y)),
      s(false, true, (x, y, z, X, Y) => x * (-y), (x, y, z, X, Y) => -(x * y)),
      s(false, isMX, (x, y, z, X, Y) => 2 * (0.5 * x), (x, y, z, X, Y) => x),
      s(false, isMX, (x, y, z, X, Y) => 2 * (x * 0.5), (x, y, z, X, Y) => x),
      s(false, isMX, (x, y, z, X, Y) => 2 * (0.5 * X), (x, y, z, X, Y) => X),
      s(false, isMX, (x, y, z, X, Y) => 2 * (X * 0.5), (x, y, z, X, Y) => X),
      s(false, true, (x, y, z, X, Y) => 0 / x, (x, y, z, X, Y) => "0"),
      s(false, true, (x, y, z, X, Y) => x / 1, (x, y, z, X, Y) => x),
      s(false, true, (x, y, z, X, Y) => x / (-1), (x, y, z, X, Y) => -x),
      s(false, true, (x, y, z, X, Y) => x / x, (x, y, z, X, Y) => "1"),
      s(false, true, (x, y, z, X, Y) => (2 * x) / 2, (x, y, z, X, Y) => x),
      s(false, true, (x, y, z, X, Y) => (x * y) / x, (x, y, z, X, Y) => y),
      s(false, isMX, (x, y, z, X, Y) => (x * y) / y, (x, y, z, X, Y) => x),
      s(false, isMX, (x, y, z, X, Y) => 1 / x, (x, y, z, X, Y) => ca.pow(x, -1)),
      s(false, true, (x, y, z, X, Y) => x / (1 / y), (x, y, z, X, Y) => x * y),
      s(false, true, (x, y, z, X, Y) => (2 * x) / (2 * y), (x, y, z, X, Y) => x / y),
      s(false, true, (x, y, z, X, Y) => (x / 5) / 0.2, (x, y, z, X, Y) => x),
      s(false, true, (x, y, z, X, Y) => x / (2 * x), (x, y, z, X, Y) => "0.5"),
      s(false, true, (x, y, z, X, Y) => (-x) / x, (x, y, z, X, Y) => "-1"),
      s(false, true, (x, y, z, X, Y) => x / (-x), (x, y, z, X, Y) => "-1"),
      s(false, true, (x, y, z, X, Y) => (-x) / (-x), (x, y, z, X, Y) => "1"),
      s(false, true, (x, y, z, X, Y) => (x / y) / x, (x, y, z, X, Y) => 1 / y),
      s(false, isMX, (x, y, z, X, Y) => (-x) / y, (x, y, z, X, Y) => -(x / y)),
      s(false, true, (x, y, z, X, Y) => x / (-y), (x, y, z, X, Y) => -(x / y)),
      s(false, true, (x, y, z, X, Y) => ca.pow(x, 0), (x, y, z, X, Y) => "1"),
      s(false, true, (x, y, z, X, Y) => ca.pow(x, 2), (x, y, z, X, Y) => x * x),
      s(false, true, (x, y, z, X, Y) => ca.pow(x, 3), (x, y, z, X, Y) => x * (x * x)),
      s(false, true, (x, y, z, X, Y) => ca.pow(x, -3), (x, y, z, X, Y) => 1 / ca.pow(x, 3)),
      s(false, true, (x, y, z, X, Y) => ca.pow(x, 4), (x, y, z, X, Y) => ca.pow(x, 2) * ca.pow(x, 2)),
      s(false, true, (x, y, z, X, Y) => ca.pow(x, 0.5), (x, y, z, X, Y) => ca.pow(x, 0.5)),
      s(false, false, (x, y, z, X, Y) => ca.pow(x, y), (x, y, z, X, Y) => ca.pow(x, y)),
      s(false, true, (x, y, z, X, Y) => ca.pow(x, 2) >= 0, (x, y, z, X, Y) => "1"),
      s(false, true, (x, y, z, X, Y) => 2 * ca.pow(x, 2) >= ca.pow(x, 2), (x, y, z, X, Y) => "1"),
      s(false, true, (x, y, z, X, Y) => ca.fmin(x, inf), (x, y, z, X, Y) => x),
      s(false, true, (x, y, z, X, Y) => ca.fmin(inf, x), (x, y, z, X, Y) => x),
      s(false, true, (x, y, z, X, Y) => ca.fmin(-inf, x), (x, y, z, X, Y) => "-inf"),
      s(false, true, (x, y, z, X, Y) => ca.fmin(x, x), (x, y, z, X, Y) => x),
      s(false, true, (x, y, z, X, Y) => ca.fmax(-inf, x), (x, y, z, X, Y) => x),
      s(false, true, (x, y, z, X, Y) => ca.fmax(x, -inf), (x, y, z, X, Y) => x),
      s(false, true, (x, y, z, X, Y) => ca.fmax(x, inf), (x, y, z, X, Y) => "inf"),
      s(false, true, (x, y, z, X, Y) => ca.fmax(x, x), (x, y, z, X, Y) => x),
      s(false, true, (x, y, z, X, Y) => ca.pow(x, 2) < 0, (x, y, z, X, Y) => "0"),
      s(false, true, (x, y, z, X, Y) => ca.eq(x, x), (x, y, z, X, Y) => "1"),
      s(false, true, (x, y, z, X, Y) => ca.ne(x, x), (x, y, z, X, Y) => "0"),
      s(false, true, (x, y, z, X, Y) => ca.pow(ca.sqrt(x), 2), (x, y, z, X, Y) => x),
      s(false, true, (x, y, z, X, Y) => ca.pow(ca.sqrt(X), 2), (x, y, z, X, Y) => X),
      s(false, true, (x, y, z, X, Y) => ca.pow(-x, 2), (x, y, z, X, Y) => ca.pow(x, 2)),
      s(false, isMX, (x, y, z, X, Y) => ca.abs(ca.abs(x)), (x, y, z, X, Y) => ca.abs(x)),
      s(false, isMX, (x, y, z, X, Y) => ca.abs(ca.sqrt(x)), (x, y, z, X, Y) => ca.sqrt(x)),
      s(false, isMX, (x, y, z, X, Y) => ca.abs(ca.exp(x)), (x, y, z, X, Y) => ca.exp(x)),
      s(false, isMX, (x, y, z, X, Y) => ca.log(ca.exp(x)), (x, y, z, X, Y) => x),
      s(false, isMX, (x, y, z, X, Y) => 1 / (1 / x), (x, y, z, X, Y) => x),
      s(false, true, (x, y, z, X, Y) => ca.sqrt(ca.pow(x, 2)), (x, y, z, X, Y) => ca.abs(x)),
      s(false, isMX, (x, y, z, X, Y) => ca.abs(-x), (x, y, z, X, Y) => ca.abs(x)),
      s(false, true, (x, y, z, X, Y) => ca.pow(ca.abs(x), 2), (x, y, z, X, Y) => ca.pow(x, 2)),
      s(false, isMX, (x, y, z, X, Y) => ca.cos(-x), (x, y, z, X, Y) => ca.cos(x)),
      s(false, isMX, (x, y, z, X, Y) => ca.cos(ca.abs(x)), (x, y, z, X, Y) => ca.cos(x)),
      s(false, false, (x, y, z, X, Y) => -(-x), (x, y, z, X, Y) => x),
      s(false, true, (x, y, z, X, Y) => ca.cosh(x * 0), (x, y, z, X, Y) => x * 0 + 1),
      s(false, true, (x, y, z, X, Y) => x / 0.5, (x, y, z, X, Y) => 2 * x),
      s(false, true, (x, y, z, X, Y) => x + x, (x, y, z, X, Y) => 2 * x),
      s(false, true, (x, y, z, X, Y) => x < x, (x, y, z, X, Y) => "0"),
      s(false, true, (x, y, z, X, Y) => 2 * x + x, (x, y, z, X, Y) => 3 * x),
      s(false, true, (x, y, z, X, Y) => x + 2 * x, (x, y, z, X, Y) => 3 * x),
      s(false, isMX, (x, y, z, X, Y) => 2 * x - x, (x, y, z, X, Y) => x),
      s(false, true, (x, y, z, X, Y) => x - 2 * x, (x, y, z, X, Y) => -x),
      s(true, false, (x, y, z, X, Y) => -(x - y), (x, y, z, X, Y) => y - x),
      s(true, false, (x, y, z, X, Y) => -(x - y) + ca.pow(x - y, 2), (x, y, z, X, Y) => -dx + ca.pow(dx, 2)),
      s(true, false, (x, y, z, X, Y) => (x + y / z) * z, (x, y, z, X, Y) => x * z + y),
      s(true, false, (x, y, z, X, Y) => (x / z + y) * z, (x, y, z, X, Y) => x + y * z),
      s(true, false, (x, y, z, X, Y) => z * (x / z + y), (x, y, z, X, Y) => x + z * y),
      s(true, false, (x, y, z, X, Y) => z * (x + y / z), (x, y, z, X, Y) => z * x + y),
      s(true, false, (x, y, z, X, Y) => (x - y / z) * z, (x, y, z, X, Y) => x * z - y),
      s(true, false, (x, y, z, X, Y) => (x / z - y) * z, (x, y, z, X, Y) => x - y * z),
      s(true, false, (x, y, z, X, Y) => z * (x / z - y), (x, y, z, X, Y) => x - z * y),
      s(true, false, (x, y, z, X, Y) => z * (x - y / z), (x, y, z, X, Y) => z * x - y),
    };
  }

  static Function mkFunction(string name, object[] args, object e) {
    if (e is SX) return ca.Function(name, args.Cast<SX>(), new[] { (SX)e });
    return ca.Function(name, args.Cast<MX>(), new[] { (MX)e });
  }

  static string callStr(Function f, object[] args) {
    if (args[0] is SX) return "[" + string.Join(", ", f.call(args.Cast<SX>(), true, false)) + "]";
    return "[" + string.Join(", ", f.call(args.Cast<MX>(), true, false)) + "]";
  }

  public void test_simplifications() {
    foreach (bool isSX in new[] { true, false }) {
      dynamic a, b, c, A, B;
      if (isSX) {
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
      object[] args = { a, b, c, A, B };
      dynamic x = ca.sin(a), y = ca.sin(b), z = ca.sin(c);
      dynamic X = ca.sin(A), Y = ca.sin(B);
      DM.rng(1);
      DM a0 = DM.rand((Sparsity)a.sparsity());
      DM b0 = DM.rand((Sparsity)b.sparsity());
      DM c0 = DM.rand((Sparsity)c.sparsity());
      DM A0 = DM.rand((Sparsity)A.sparsity());
      DM B0 = DM.rand((Sparsity)B.sparsity());
      dynamic x0 = ca.sin(a0), y0 = ca.sin(b0), z0 = ca.sin(c0), X0 = ca.sin(A0), Y0 = ca.sin(B0);
      dynamic dx = x - y;
      int count = 1477;
      dynamic Alast = null;
      foreach (Simp cs in simplificationCases(dx)) {
        Console.WriteLine(count);
        count++;
        object Ae = cs.genA(x, y, z, X, Y);
        object Be = cs.genB(x, y, z, X, Y);
        Alast = Ae;
        Function f = mkFunction("f", args, Ae);
        if (cs.refcount) {
          assertTrue(pyStr(Be) != pyStr(Ae), "expected different: " + pyStr(Be));
          Console.WriteLine(Ae);
          f = f.transform();
          Function fref = mkFunction("f", args, Be);
          assertEquals(callStr(fref, args), callStr(f, args));
        } else {
          assertEquals(pyStr(Be), pyStr(Ae), "case " + (count - 1));
        }
        if (cs.onTheFly) {
          GlobalOptions.setSimplificationOnTheFly(false);
          object Ae2;
          try {
            Ae2 = cs.genA(x, y, z, X, Y);
          } finally {
            GlobalOptions.setSimplificationOnTheFly(true);
          }
          assertTrue(pyStr(Be) != pyStr(Ae2), "case " + (count - 1) + ": expected no simplification: " + pyStr(Ae2));
        }
        DM Aval = f.call(new[] { a0, b0, c0, A0, B0 })[0];
        DM Aref = cs.genA(x0, y0, z0, X0, Y0);
        Console.WriteLine(Aval);
        Console.WriteLine(Aref);
        checkarray(Aval, Aref);
      }
      Console.WriteLine(ca.DM(Sparsity.upper(3), 0).is_nonnegative());
      Console.WriteLine(ca.DM(Sparsity.upper(3), 0.5).is_half());
      Console.WriteLine(ca.MX(ca.DM(Sparsity.upper(3), 0.5)).is_half());
      Console.WriteLine(ca.MX(ca.DM(Sparsity.dense(3, 3), 0.5)).is_half());
      var _ = 2 * (ca.DM(new double[,] { { 2, 3 }, { 6, 8 } }) * Alast);
    }
  }

  public void test_vertcat_empty() {
    Check.Assume(KNOWN_BUGS, "known bug");
    MX a = ca.MX(ca.DM(0, 2));
    MX v = ca.vertcat(a, a);
    assertEquals(0L, v.size1());
    assertEquals(2L, v.size2());
    a = ca.MX(ca.DM(2, 0));
    v = ca.vertcat(a, a);
    assertEquals(4L, v.size1());
    assertEquals(0L, v.size2());
  }

  public void test_jacobian_empty() {
    MX x = MX.sym("x", 3);
    MX s = ca.jacobian(ca.DM(0, 0), x);
    assertEquals(0L, s.size1());
    assertEquals(3L, s.size2());
    s = ca.jacobian(x, MX.sym("x", 0, 4));
    assertEquals(3L, s.size1());
    assertEquals(0L, s.size2());
  }

  public void test_mul_sparsity() {
    int N = 10;
    MX x = MX.sym("x", N, N);
    MX y = MX.sym("y", N, N);
    DM x_ = randDM(N, N);
    DM y_ = randDM(N, N);
    Sparsity filt = Sparsity.diag(N).plus(Sparsity.triplet(N, N, new long[] { 1 }, new long[] { 3 }));
    Function f = fun("f", new[] { x, y }, new[] { ca.mtimes(x, y) });
    Function g = fun("g", new[] { x, y }, new[] { ca.mac(x, y, MX.zeros(filt)) });
    DM fOut = ev(f, x_, y_);
    DM gOut = ev(g, x_, y_);
    checkarray(DM.ones(filt), DM.ones(g.sparsity_out(0)));
    checkarray(fOut.get(filt), gOut);
  }

  public void test_mul_zero_wrong() {
    Check.Throws<SwigException>(() => ca.mtimes(MX.sym("X", 4, 5), MX.zeros(3, 2)));
  }

  static DM lowerRange(Function f) => ca.DM(f.sparsity_in(0), ca.vertcat(range(15)));

  public void test_vertsplit() {
    MX a = MX.sym("X", Sparsity.lower(5));
    List<MX> v = ca.vertsplit(a, new long[] { 0, 2, 4, 5 });
    Function f = fun("f", new[] { a }, v);
    List<DM> V = f.call(dms(lowerRange(f)));
    assertEquals(3, V.Count);
    checkarray(V[0], ca.DM(new double[,] { { 0, 0, 0, 0, 0 }, { 1, 5, 0, 0, 0 } }));
    checkarray(V[1], ca.DM(new double[,] { { 2, 6, 9, 0, 0 }, { 3, 7, 10, 12, 0 } }));
    checkarray(V[2], ca.DM(new double[,] { { 4, 8, 11, 13, 14 } }));

    v = ca.vertsplit(a);
    f = fun("f", new[] { a }, v);
    V = f.call(dms(lowerRange(f)));
    assertEquals(a.size1(), (long)V.Count);
    checkarray(V[0], ca.DM(new double[,] { { 0, 0, 0, 0, 0 } }));
    checkarray(V[1], ca.DM(new double[,] { { 1, 5, 0, 0, 0 } }));
    checkarray(V[2], ca.DM(new double[,] { { 2, 6, 9, 0, 0 } }));
    checkarray(V[3], ca.DM(new double[,] { { 3, 7, 10, 12, 0 } }));
    checkarray(V[4], ca.DM(new double[,] { { 4, 8, 11, 13, 14 } }));

    v = ca.vertsplit(a, 2);
    f = fun("f", new[] { a }, v);
    V = f.call(dms(lowerRange(f)));
    assertEquals(3, V.Count);
    checkarray(V[0], ca.DM(new double[,] { { 0, 0, 0, 0, 0 }, { 1, 5, 0, 0, 0 } }));
    checkarray(V[1], ca.DM(new double[,] { { 2, 6, 9, 0, 0 }, { 3, 7, 10, 12, 0 } }));
    checkarray(V[2], ca.DM(new double[,] { { 4, 8, 11, 13, 14 } }));

    v = ca.vertsplit(a, new long[] { 0, 0, 3, a.size1() });
    f = fun("f", new[] { a }, v);
    List<DM> W = evs(f, lowerRange(f));
    assertEquals(3, v.Count);
    assertEquals(0L, v[0].size1());
    assertEquals(5L, v[0].size2());
    checkarray(W[1], ca.DM(new double[,] { { 0, 0, 0, 0, 0 }, { 1, 5, 0, 0, 0 }, { 2, 6, 9, 0, 0 } }));
    checkarray(W[2], ca.DM(new double[,] { { 3, 7, 10, 12, 0 }, { 4, 8, 11, 13, 14 } }));
  }

  public void test_horzsplit() {
    MX a = MX.sym("X", Sparsity.lower(5));
    List<MX> v = ca.horzsplit(a, new long[] { 0, 2, 4, 5 });
    Function f = fun("f", new[] { a }, v);
    List<DM> V = evs(f, lowerRange(f));
    assertEquals(3, V.Count);
    checkarray(V[0], ca.DM(new double[,] { { 0, 0 }, { 1, 5 }, { 2, 6 }, { 3, 7 }, { 4, 8 } }));
    checkarray(V[1], ca.DM(new double[,] { { 0, 0 }, { 0, 0 }, { 9, 0 }, { 10, 12 }, { 11, 13 } }));
    checkarray(V[2], ca.DM(new double[,] { { 0 }, { 0 }, { 0 }, { 0 }, { 14 } }));

    v = ca.horzsplit(a);
    f = fun("f", new[] { a }, v);
    V = evs(f, lowerRange(f));
    assertEquals(a.size1(), (long)V.Count);
    checkarray(V[0], ca.vertcat(0, 1, 2, 3, 4));
    checkarray(V[1], ca.vertcat(0, 5, 6, 7, 8));
    checkarray(V[2], ca.vertcat(0, 0, 9, 10, 11));
    checkarray(V[3], ca.vertcat(0, 0, 0, 12, 13));
    checkarray(V[4], ca.vertcat(0, 0, 0, 0, 14));

    v = ca.horzsplit(a, 2);
    f = fun("f", new[] { a }, v);
    V = evs(f, lowerRange(f));
    assertEquals(3, V.Count);
    checkarray(V[0], ca.DM(new double[,] { { 0, 0 }, { 1, 5 }, { 2, 6 }, { 3, 7 }, { 4, 8 } }));
    checkarray(V[1], ca.DM(new double[,] { { 0, 0 }, { 0, 0 }, { 9, 0 }, { 10, 12 }, { 11, 13 } }));
    checkarray(V[2], ca.DM(new double[,] { { 0 }, { 0 }, { 0 }, { 0 }, { 14 } }));

    v = ca.horzsplit(a, new long[] { 0, 0, 3, a.size2() });
    f = fun("f", new[] { a }, v);
    V = evs(f, lowerRange(f));
    assertEquals(3, v.Count);
    assertEquals(5L, v[0].size1());
    assertEquals(0L, v[0].size2());
    checkarray(V[1], ca.DM(new double[,] { { 0, 0, 0 }, { 1, 5, 0 }, { 2, 6, 9 }, { 3, 7, 10 }, { 4, 8, 11 } }));
    checkarray(V[2], ca.DM(new double[,] { { 0, 0 }, { 0, 0 }, { 0, 0 }, { 12, 0 }, { 13, 14 } }));
  }

  public void test_blocksplit() {
    MX a = MX.sym("X", Sparsity.lower(5));
    List<List<MX>> v = ca.blocksplit(a, new long[] { 0, 2, 4, 5 }, new long[] { 0, 1, 3, 5 });
    var fs = v.Select(vr => fun("fs", new[] { a }, vr)).ToList();
    var w = new List<List<DM>>();
    for (int i = 0; i < 3; i++) w.Add(fs[i].call(dms(lowerRange(fs[i]))));
    checkarray(w[0][0], ca.vertcat(0, 1));
    checkarray(w[0][1], ca.DM(new double[,] { { 0, 0 }, { 5, 0 } }));
    checkarray(w[1][0], ca.vertcat(2, 3));
    checkarray(ca.blockcat(w), lowerRange(fs[0]));
  }

  public void test_is_one_etc() {
    double inf = double.PositiveInfinity;
    var ops = new List<Func<DM, dynamic>> { d => d, d => ca.MX(d), d => ca.SX(d) };
    foreach (var op in ops) {
      assertTrue(op(ca.DM(Sparsity.upper(3), 0)).is_zero());
      assertTrue(op(ca.DM(Sparsity.dense(3, 4), 0)).is_zero());
      assertTrue(op(ca.DM(3, 4)).is_zero());
      assertFalse(op(ca.DM(new double[,] { { 0, 0 }, { 0, 7 } })).is_zero());

      assertFalse(op(ca.DM(Sparsity.upper(3), 1)).is_one());
      assertTrue(op(ca.DM(Sparsity.dense(3, 4), 1)).is_one());
      assertFalse(op(ca.DM(new double[,] { { 1, 1 }, { 1, 7 } })).is_one());

      assertFalse(op(ca.DM(Sparsity.upper(3), 0.5)).is_half());
      assertTrue(op(ca.DM(Sparsity.dense(3, 4), 0.5)).is_half());
      assertFalse(op(ca.DM(new double[,] { { 0.5, 0.5 }, { 0.5, 7 } })).is_half());

      assertFalse(op(ca.DM(Sparsity.upper(3), -1)).is_minus_one());
      assertTrue(op(ca.DM(Sparsity.dense(3, 4), -1)).is_minus_one());
      assertFalse(op(ca.DM(new double[,] { { -1, -1 }, { -1, 7 } })).is_minus_one());

      assertFalse(op(ca.DM(Sparsity.upper(3), inf)).is_inf());
      assertTrue(op(ca.DM(Sparsity.dense(3, 4), inf)).is_inf());
      assertFalse(op(ca.DM(new double[,] { { inf, inf }, { inf, 7 } })).is_inf());

      assertFalse(op(ca.DM(Sparsity.upper(3), -inf)).is_minus_inf());
      assertTrue(op(ca.DM(Sparsity.dense(3, 4), -inf)).is_minus_inf());
      assertFalse(op(ca.DM(new double[,] { { -inf, -inf }, { -inf, 7 } })).is_minus_inf());

      assertTrue(op(ca.DM(Sparsity.upper(3), 3)).is_integer());
      assertTrue(op(ca.DM(Sparsity.dense(3, 4), 3)).is_integer());
      assertFalse(op(ca.DM(new double[,] { { 3, 3 }, { 3, 0.7 } })).is_integer());
      assertTrue(op(ca.DM(new double[,] { { 3, 3 }, { 3, 7 } })).is_integer());

      assertFalse(op(ca.DM(Sparsity.upper(3), 7)).is_value(7.0));
      assertTrue(op(ca.DM(Sparsity.dense(3, 4), 7)).is_value(7.0));
      assertFalse(op(ca.DM(new double[,] { { 7, 7 }, { 7, 1 } })).is_value(7.0));
      assertTrue(op(ca.DM(Sparsity.upper(3), 0)).is_value(0.0));
      assertTrue(op(ca.DM(Sparsity.dense(3, 4), 0)).is_value(0.0));
      assertTrue(op(ca.DM(3, 4)).is_value(0.0));
      assertFalse(op(ca.DM(new double[,] { { 0, 0 }, { 0, 7 } })).is_value(0.0));

      assertFalse(op(ca.DM(Sparsity.upper(3), 1)).is_eye());
      assertTrue(op(ca.DM(Sparsity.diag(3), 1)).is_eye());
      assertFalse(op(ca.DM(new double[,] { { 1, 0 }, { 0, 1 } })).is_eye());

      assertFalse(op(ca.DM(Sparsity.upper(3), -1)).is_nonnegative());
      assertTrue(op(ca.DM(Sparsity.upper(3), 1)).is_nonnegative());
      assertTrue(op(ca.DM(Sparsity.dense(3, 4), 1)).is_nonnegative());
      assertTrue(op(ca.DM(new double[,] { { 1, 1 }, { 1, 7 } })).is_nonnegative());
      assertFalse(op(ca.DM(new double[,] { { 1, 1 }, { 1, -7 } })).is_nonnegative());
    }
  }

  public void test_mxnulloutput() {
    MX a = ca.MX(5, 0);
    MX b = MX.sym("x", 2);
    Function f = fun("f", new[] { b }, new[] { a });
    MX c = f.Call(b);
    assertEquals(5L, c.size1());
    assertEquals(0L, c.size2());
    c = f.call(new[] { b }, true)[0];
    assertEquals(5L, c.size1());
    assertEquals(0L, c.size2());

    a = ca.MX(0, 0);
    b = MX.sym("x", 2);
    f = fun("f", new[] { b }, new[] { a });
    c = f.Call(b);
    assertEquals(0L, c.size1());
    assertEquals(0L, c.size2());
    c = f.call(new[] { b }, true)[0];
    assertEquals(0L, c.size1());
  }

  public void test_mxnull() {
    MX c = ca.mtimes(ca.MX(5, 0), ca.MX(0, 3));
    assertEquals(0L, c.nnz());
    c = ca.mtimes(ca.MX(5, 3), ca.MX(3, 4));
    assertEquals(0L, c.nnz());
  }

  public void test_mxnullop() {
    MX c = ca.MX(0, 0);
    MX x = MX.sym("x", 2, 3);
    Check.Throws<SwigException>(() => { var _ = x + c; });
    Check.Throws<SwigException>(() => { var _ = x / c; });
  }

  static double minNz(DM c) => c.nonzeros().DefaultIfEmpty(double.PositiveInfinity).Min();

  public void test_MX_shapes() {
    Check.Assume(RunSlow, "slow");
    message("MX unary operations");
    var sp43 = ca.Sparsity(4, 3, new long[] { 0, 2, 2, 3 }, new long[] { 1, 2, 1 });
    foreach (Sparsity sp in new[] { Sparsity.dense(0, 0), Sparsity.dense(0, 2), Sparsity.dense(2, 0), Sparsity.dense(1, 1),
                                    Sparsity.dense(2, 2), sp43 }) {
      foreach (double v in new double[] { 0, 1, 0.2 }) {
        DM x_ = ca.DM(sp, v);
        MX xx = MX.sym("x", sp.size1(), sp.size2());
        MX x = xx.get(sp);
        foreach (Op1 op in pool) {
          if (op.nozero && v == 0) continue;
          MX r = op.cas(x);
          Function f = fun("f", new[] { xx }, new[] { r });
          DM rr = ev(f, v);
          Console.WriteLine(r + " " + rr + " " + op.@ref(x_));
          checkarray(rr, op.@ref(x_));
          DM a = ca.DM(f.sparsity_out(0), 1);
          DM b = DM.ones(ca.DM(op.@ref(x_)).sparsity());
          DM c = b - a;
          if (c.nnz() > 0) assertTrue(minNz(c) >= 0, sp + " " + v + " " + op.name);
        }
      }
    }
    var sps = new[] { ca.Sparsity(1, 1), Sparsity.dense(1, 1), ca.Sparsity(3, 4), Sparsity.dense(3, 4), sp43.T() };
    foreach (Sparsity sp in sps) {
      foreach (double v1 in new double[] { 0, 1, 0.2, -0.2 }) {
        DM x1_ = ca.DM(sp, v1);
        MX xx1 = MX.sym("x", sp.size1(), sp.size2());
        MX x1 = xx1.get(sp);
        SX xx1s = SX.sym("x", sp.size1(), sp.size2());
        xx1s.get(sp);
        foreach (Sparsity sp2 in sps) {
          foreach (double v2 in new double[] { 0, 1, 0.2, -0.2 }) {
            DM x2_ = ca.DM(sp2, v2);
            MX xx2 = MX.sym("x", sp2.size1(), sp2.size2());
            MX x2 = xx2.get(sp2);
            SX xx2s = SX.sym("x", sp2.size1(), sp2.size2());
            xx2s.get(sp2);
            foreach (Op2 op in matrixbinarypool) {
              if ((op.name.Contains("mul") || op.name.Contains("mtimes")) && (sp.numel() == 1 || sp2.numel() == 1)) continue;
              MX r = op.cas(x1, x2);
              Function f = fun("f", new[] { xx1, xx2 }, new[] { r });
              var fIn = dms(v1, v2);
              DM fOut = ev(f, v1, v2);
              Function g = fun("g", new[] { xx1, xx2 }, new[] { r });
              ev(g, v1, v2);
              string msg = sp + " " + sp2 + " " + v1 + " " + v2 + " " + op.name;
              checkarray(fOut, op.@ref(x1_, x2_), msg);
              if (!op.name.Contains("mul")) {
                DM a = DM.ones(f.sparsity_out(0));
                DM b = DM.ones(g.sparsity_out(0));
                DM c = b - a;
                if (c.nnz() > 0) assertTrue(minNz(c) >= 0, msg);
              }
              if (sp.nnz() > 0 && sp2.nnz() > 0 && v1 != 0 && v2 != 0) checkfunction(f, g, fIn);
            }
          }
        }
      }
    }
  }

  public void test_MXConstant() {
    message("MX unary operations, constant");
    var sp43 = ca.Sparsity(4, 3, new long[] { 0, 2, 2, 3 }, new long[] { 1, 2, 1 });
    foreach (Sparsity sp in new[] { Sparsity.dense(0, 0), Sparsity.dense(0, 2), Sparsity.dense(2, 0), Sparsity.dense(1, 1),
                                    Sparsity.dense(2, 2), sp43 }) {
      foreach (double v in new double[] { 0, 1, 0.2 }) {
        DM x_ = ca.DM(sp, v);
        MX x = ca.MX(sp, v);
        foreach (Op1 op in pool) {
          if (op.nozero && (v == 0 || !sp.is_dense())) continue;
          MX r = op.cas(x);
          Console.WriteLine(r);
          assertTrue(r.is_constant());
          checkarray(r.to_DM(), op.@ref(x_), x_ + " " + op.name);
          DM a = DM.ones(r.to_DM().sparsity());
          DM b = DM.ones(ca.DM(op.@ref(x_)).sparsity());
          DM c = b - a;
          if (c.nnz() > 0) assertTrue(minNz(c) >= 0, sp + " " + v + " " + op.name);
        }
      }
    }
    var sps = new[] { Sparsity.dense(1, 1), ca.Sparsity(1, 1), ca.Sparsity(3, 4), Sparsity.dense(3, 4), sp43.T() };
    foreach (Sparsity sp in sps) {
      foreach (double v1 in new double[] { 0, 1, 0.2, -0.2 }) {
        DM x1_ = ca.DM(sp, v1);
        MX x1 = ca.MX(sp, v1);
        foreach (Sparsity sp2 in sps) {
          foreach (double v2 in new double[] { 0, 1, 0.2, -0.2 }) {
            DM x2_ = ca.DM(sp2, v2);
            MX x2 = ca.MX(sp2, v2);
            foreach (Op2 op in matrixbinarypool) {
              if ((op.name.Contains("mul") || op.name.Contains("mtimes")) && (sp.numel() == 1 || sp2.numel() == 1)) continue;
              MX r = op.cas(x1, x2);
              Function f = fun("f", new MX[0], new[] { r });
              List<DM> fOut = f.call(none);
              string msg = sp + " " + sp2 + " " + v1 + " " + v2 + " " + op.name;
              checkarray(fOut[0], op.@ref(x1_, x2_), msg);
              if (!op.name.Contains("mul") && !op.name.Contains("mtimes")) {
                DM a = DM.ones(f.sparsity_out(0));
                DM b = DM.ones(ca.DM(op.@ref(x1_, x2_)).sparsity());
                DM c = b - a;
                if (c.nnz() > 0) assertTrue(minNz(c) >= 0, msg);
              }
            }
          }
        }
      }
    }
  }

  public void test_graph_substitute() {
    MX x = MX.sym("X", 4, 4);
    MX y = MX.sym("Y", 4, 4);
    MX b = MX.sym("B", 4, 4);
    MX c = x * y;
    MX d = b * y;
    MX f = c + d;
    MX C = MX.sym("C", 4, 4);
    f = ca.graph_substitute(f, new[] { c }, new[] { C });
    Function F = fun("F", new[] { y, b, C }, new[] { f });
    checkarray(ev(F, 1, 2, 3), 5.0 * DM.ones(4, 4));

    MX D = MX.sym("D", 4, 4);
    f = ca.graph_substitute(f, new[] { d }, new[] { D });
    F = fun("F", new[] { D, C }, new[] { f });
    checkarray(ev(F, 4, 5), 9.0 * DM.ones(4, 4));

    MX z = ca.vertcat(c, d, c + d);
    f = ca.graph_substitute(z, new[] { d }, new[] { D });
    F = fun("F", new[] { x, y, b, D }, new[] { f });
    Function Fref = fun("F", new[] { x, y, b, D }, new[] { ca.vertcat(c, D, c + D) });
    checkfunctionLight(F, Fref, dms(3, 5, 7, 12));
  }

  public void test_graph_substitute_function() {
    MX x = MX.sym("x");
    MX y = MX.sym("y");
    Function f = ca.Function("forig", new[] { x, y }, new[] { x * y }, opts(("never_inline", true)));
    Function g = ca.Function("fsubs", new[] { x, y }, new[] { x * y }, opts(("never_inline", true)));
    MX sx = ca.sin(x);
    MX y2 = 2.0 * y;
    MX a = f.Call(sx, y2);
    MX x2 = ca.pow(x, 2.0);
    MX z = ca.sqrt(x2 * a) + ca.pow(a, 2.0);
    Console.WriteLine(z);
    MX z2 = ca.graph_substitute(z, new[] { x2 }, new[] { x });
    Console.WriteLine(z2);

    f = ca.Function("forig", new[] { x, y }, new[] { x * y, x - y, x + y }, opts(("never_inline", true)));
    g = ca.Function("fsubs", new[] { x, y }, new[] { x * y, x - y, x + y }, opts(("never_inline", true)));
    sx = ca.sin(x);
    y2 = 2.0 * y;
    List<MX> abc = evs(f, sx, y2);
    a = abc[0];
    MX b = abc[1], c = abc[2];
    z = ca.sqrt(a) + ca.pow(c, 2.0);
    Console.WriteLine(z);
    z2 = ca.graph_substitute(z, new MX[0], new MX[0]);
    Console.WriteLine(z2);

    z2 = ca.graph_substitute(z, new[] { a }, new[] { evs(g, sx, y2)[0] });
    Console.WriteLine(z2);
    assertTrue(z2.ToString().Contains("fsubs"));
    assertTrue(z2.ToString().Contains("forig"));
    Function Fref = fun("F_ref", new[] { x, y }, new[] { z });
    Function F = fun("F", new[] { x, y }, new[] { z2 });
    checkfunctionLight(F, Fref, dms(3.1, 7.1));

    z2 = ca.graph_substitute(z2, new[] { b }, new[] { ca.MX(double.NaN) });
    assertTrue(z2.ToString().Contains("fsubs"));
    assertTrue(z2.ToString().Contains("forig"));
    Fref = fun("F_ref", new[] { x, y }, new[] { z });
    F = fun("F", new[] { x, y }, new[] { z2 });
    checkfunctionLight(F, Fref, dms(3.1, 7.1));

    z2 = ca.graph_substitute(z2, new[] { c }, new[] { evs(g, sx, y2)[2] });
    assertTrue(z2.ToString().Contains("fsubs"));
    assertFalse(z2.ToString().Contains("forig"));
    Fref = fun("F_ref", new[] { x, y }, new[] { z });
    F = fun("F", new[] { x, y }, new[] { z2 });
    checkfunctionLight(F, Fref, dms(3.1, 7.1));

    z2 = ca.graph_substitute(z, new[] { a.dep(0) }, new[] { evs(g, sx, y2)[0].dep(0) });
    assertTrue(z2.ToString().Contains("fsubs"));
    assertFalse(z2.ToString().Contains("forig"));
    Fref = fun("F_ref", new[] { x, y }, new[] { z });
    F = fun("F", new[] { x, y }, new[] { z2 });
    checkfunctionLight(F, Fref, dms(3.1, 7.1));
  }

  public void test_graph_substitute_call() {
    MX x = MX.sym("x");
    MX y = MX.sym("y");
    MX z = MX.sym("z");
    Function f = fun("f", new[] { x, y, z }, new[] { ca.pow(x, 2.0) });
    MX sx = ca.sin(x);
    MX x1 = x + 1.0;
    List<MX> res = ca.graph_substitute(new[] { f.Call(x1, y, z), ca.sin(f.Call(x1, y, z)), f.Call(sx, y, z), ca.sin(f.Call(sx, y, z)) },
                                    new[] { sx }, new[] { ca.cos(x) });
    string s = "[" + string.Join(", ", res) + "]";
    assertEquals(4, s.Split("{0}").Length - 1);
  }

  public void test_matrix_expand() {
    int n = 2;
    MX a = MX.sym("a", n, n);
    MX b = MX.sym("b", n, n);
    MX c = MX.sym("c", n, n);
    MX d = a + b;
    MX e = d * c;
    assertEquals(6L, ca.n_nodes(e));
    MX t0 = ca.matrix_expand(e);
    assertEquals(5L, ca.n_nodes(t0));
    MX t1 = ca.matrix_expand(e, new[] { d });
    assertEquals(6L, ca.n_nodes(t1));
    Console.WriteLine(e + " " + t0 + " " + t1);
    var outs = new List<DM>();
    foreach (MX x in new[] { e, t0, t1 }) {
      Function f = fun("f", new[] { a, b, c }, new[] { x });
      outs.Add(ev(f, 1.1, 2.2, 3.3));
      if (outs.Count > 1) checkarray(outs[0], outs[^1]);
    }
    Console.WriteLine("[" + string.Join(", ", outs) + "]");
  }

  /// <summary>numpy.kron on dense arrays</summary>
  static double[,] npKron(double[,] a, double[,] b) {
    int ma = a.GetLength(0), na = a.GetLength(1), mb = b.GetLength(0), nb = b.GetLength(1);
    var r = new double[ma * mb, na * nb];
    for (int i = 0; i < ma; i++)
      for (int j = 0; j < na; j++)
        for (int k = 0; k < mb; k++)
          for (int l = 0; l < nb; l++) r[i * mb + k, j * nb + l] = a[i, j] * b[k, l];
    return r;
  }

  public void test_kron() {
    DM a = ca.sparsify(ca.DM(new double[,] { { 1, 0, 6 }, { 2, 7, 0 } }));
    DM b = ca.sparsify(ca.DM(new double[,] { { 1, 0, 0 }, { 2, 3, 7 }, { 0, 0, 9 }, { 1, 12, 13 } }));
    MX A = MX.sym("A", a.sparsity());
    MX B = MX.sym("B", b.sparsity());
    MX C = ca.kron(A, B);
    Function f = fun("f", new[] { A, B }, new[] { C });
    List<DM> fIn = dms(a, b);
    DM c_ = ev(f, a, b);
    assertEquals(a.size1() * b.size1(), c_.size1());
    assertEquals(a.size2() * b.size2(), c_.size2());
    assertEquals(a.nnz() * b.nnz(), c_.nnz());
    checkarray(c_, npKron(a.ToArray(), b.ToArray()));
    assertEquals(4L, f.n_nodes());

    MX dA = MX.sym("dA", A.sparsity());
    MX dB = MX.sym("dB", B.sparsity());
    MX fwd = ca.jtimes(C, ca.vertcat(ca.vec(A), ca.vec(B)), ca.vertcat(ca.vec(dA), ca.vec(dB)));
    MX expected = ca.kron(dA, B) + ca.kron(A, dB);
    Function g = fun("g", new[] { A, B, dA, dB }, new[] { fwd - expected });
    DM.rng(0);
    DM diff = ev(g, a, b, DM.rand(A.sparsity()), DM.rand(B.sparsity()));
    assertTrue((double)ca.norm_inf(diff) < 1e-14);

    checkCodegen(f, fIn);
    checkSerialize(f, fIn);
    checkfunction(f, f.expand(), fIn);

    SX As = SX.sym("A", a.sparsity());
    SX Bs = SX.sym("B", b.sparsity());
    SX Cs = ca.kron(As, Bs);
    Function fSX = ca.Function("f", new[] { As, Bs }, new[] { Cs });
    checkfunction(f, fSX, fIn);
  }

  public void test_kron_flavors() {
    DM.rng(11);
    Sparsity spD2 = Sparsity.dense(2, 3);
    Sparsity spD3 = Sparsity.dense(3, 2);
    Sparsity spS2 = ca.sparsify(ca.DM(new double[,] { { 1, 0, 1 }, { 0, 1, 0 } })).sparsity();
    Sparsity spS3 = ca.sparsify(ca.DM(new double[,] { { 1, 0 }, { 1, 1 }, { 0, 1 } })).sparsity();
    string[] labels = { "dense+dense", "dense+sparse", "sparse+dense", "sparse+sparse" };
    Sparsity[][] cases = { new[] { spD2, spD3 }, new[] { spD2, spS3 }, new[] { spS2, spD3 }, new[] { spS2, spS3 } };
    for (int i = 0; i < 4; i++) {
      string label = labels[i];
      MX A = MX.sym("A", cases[i][0]);
      MX B = MX.sym("B", cases[i][1]);
      MX K = ca.kron(A, B);
      Function f = fun("f", new[] { A, B }, new[] { K });
      assertEquals(4L, f.n_nodes(), label + ": graph not collapsed");
      DM aVal = DM.rand(cases[i][0]);
      DM bVal = DM.rand(cases[i][1]);
      checkarray(ev(f, aVal, bVal), (DM)npKron(aVal.ToArray(), bVal.ToArray()), label);
      checkCodegen(f, dms(aVal, bVal));
      checkSerialize(f, dms(aVal, bVal));
      checkfunction(f, f.expand(), dms(aVal, bVal));
    }
  }

  static double[,] refKronContract(double[,] m, double[,] x, bool inner) {
    int mAmB = m.GetLength(0), nAnB = m.GetLength(1), xrow = x.GetLength(0), xcol = x.GetLength(1);
    if (inner) {
      int mB = xrow, nB = xcol, mA = mAmB / mB, nA = nAnB / nB;
      var y = new double[mA, nA];
      for (int i = 0; i < mA; i++)
        for (int j = 0; j < nA; j++)
          for (int r = 0; r < mB; r++)
            for (int s = 0; s < nB; s++) y[i, j] += m[i * mB + r, j * nB + s] * x[r, s];
      return y;
    } else {
      int mA = xrow, nA = xcol, mB = mAmB / mA, nB = nAnB / nA;
      var y = new double[mB, nB];
      for (int r = 0; r < mB; r++)
        for (int s = 0; s < nB; s++)
          for (int i = 0; i < mA; i++)
            for (int j = 0; j < nA; j++) y[r, s] += x[i, j] * m[i * mB + r, j * nB + s];
      return y;
    }
  }

  public void test_kron_contract_flavors() {
    DM.rng(13);
    Sparsity spD2 = Sparsity.dense(2, 3);
    Sparsity spD3 = Sparsity.dense(3, 2);
    Sparsity spS2 = ca.sparsify(ca.DM(new double[,] { { 1, 0, 1 }, { 0, 1, 0 } })).sparsity();
    Sparsity spS3 = ca.sparsify(ca.DM(new double[,] { { 1, 0 }, { 1, 1 }, { 0, 1 } })).sparsity();
    string[] labels = { "dense_M+dense_X", "dense_M+sparse_X", "sparse_M+dense_X", "sparse_M+sparse_X" };
    Sparsity[][] combos = { new[] { spD2, spD3 }, new[] { spD2, spS3 }, new[] { spS2, spD3 }, new[] { spS2, spS3 } };
    for (int i = 0; i < 4; i++) {
      string label = labels[i];
      Sparsity spA = combos[i][0], spB = combos[i][1];
      Sparsity spM = ca.kron(spA, spB);
      MX Mi = MX.sym("M", spM);
      MX Xi = MX.sym("X", spB);
      MX Yi = ca.kron_contract(Mi, Xi, true);
      Function fi = fun("fi", new[] { Mi, Xi }, new[] { Yi });
      DM mv = DM.rand(spM), xv = DM.rand(spB);
      double[,] @ref = refKronContract(mv.ToArray(), xv.ToArray(), true);
      checkarray(ca.DM(ev(fi, mv, xv)), (DM)@ref * DM.ones(Yi.sparsity()), label + " inner");
      checkCodegen(fi, dms(mv, xv));
      checkSerialize(fi, dms(mv, xv));
      checkfunction(fi, fi.expand(), dms(mv, xv));

      MX Mo = MX.sym("M", spM);
      MX Xo = MX.sym("X", spA);
      MX Yo = ca.kron_contract(Mo, Xo, false);
      Function fo = fun("fo", new[] { Mo, Xo }, new[] { Yo });
      DM mv2 = DM.rand(spM), xv2 = DM.rand(spA);
      double[,] ref2 = refKronContract(mv2.ToArray(), xv2.ToArray(), false);
      checkarray(ca.DM(ev(fo, mv2, xv2)), (DM)ref2 * DM.ones(Yo.sparsity()), label + " outer");
      checkCodegen(fo, dms(mv2, xv2));
      checkSerialize(fo, dms(mv2, xv2));
      checkfunction(fo, fo.expand(), dms(mv2, xv2));
    }
  }

  public void test_kron_contract_empty_x() {
    var spXEmpty = ca.Sparsity(2, 2);
    Sparsity spM = Sparsity.dense(4, 4);
    foreach (bool inner in new[] { true, false }) {
      MX M = MX.sym("M", spM);
      MX X = MX.sym("X", spXEmpty);
      MX Y = ca.kron_contract(M, X, inner);
      assertEquals(2L, Y.size1());
      assertEquals(2L, Y.size2());
      assertEquals(0L, Y.nnz());
      Function f = fun("f", new[] { M, X }, new[] { Y });
      DM.rng(0);
      DM mv = DM.rand(spM), xv = ca.DM(spXEmpty);
      DM @out = ev(f, mv, xv);
      assertEquals(2L, @out.size1());
      assertEquals(2L, @out.size2());
      assertEquals(0L, @out.nnz());
      checkCodegen(f, dms(mv, xv));
      checkSerialize(f, dms(mv, xv));
    }
  }

  public void test_kron_contract() {
    DM.rng(7);
    DM a = ca.sparsify(ca.DM(new double[,] { { 1, 0, 6 }, { 2, 7, 0 } }));
    DM b = ca.sparsify(ca.DM(new double[,] { { 1, 0, 0 }, { 2, 3, 7 }, { 0, 0, 9 }, { 1, 12, 13 } }));
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

    Function fInner = fun("fi", new[] { M, B }, new[] { Yi });
    Function fOuter = fun("fo", new[] { M, A }, new[] { Yo });

    DM mVal = DM.rand(K.sparsity());
    double[,] mNp = mVal.ToArray(), aNp = a.ToArray(), bNp = b.ToArray();
    int mB1 = (int)b.size1(), nB1 = (int)b.size2(), mA1 = (int)a.size1(), nA1 = (int)a.size2();
    double[,] refInner = new double[mA1, nA1], refOuter = new double[mB1, nB1];
    for (int i = 0; i < mA1; i++)
      for (int j = 0; j < nA1; j++)
        for (int r = 0; r < mB1; r++)
          for (int s = 0; s < nB1; s++) refInner[i, j] += mNp[i * mB1 + r, j * nB1 + s] * bNp[r, s];
    for (int r = 0; r < mB1; r++)
      for (int s = 0; s < nB1; s++)
        for (int i = 0; i < mA1; i++)
          for (int j = 0; j < nA1; j++) refOuter[r, s] += aNp[i, j] * mNp[i * mB1 + r, j * nB1 + s];
    checkarray(ev(fInner, mVal, b), refInner);
    checkarray(ev(fOuter, mVal, a), refOuter);

    int[] sizes = { 3, 5, 10, 20 };
    var nodeCounts = new List<long>();
    foreach (int n in sizes) {
      MX Ax = MX.sym("Ax", n, n);
      MX Bx = MX.sym("Bx", 3, 3);
      MX Kx = ca.kron(Ax, Bx);
      MX gA = ca.jacobian(ca.vec(Kx), ca.vec(Ax));
      nodeCounts.Add(fun("f", new[] { Ax, Bx }, new[] { gA }).n_nodes());
    }
    assertEquals(nodeCounts.Min(), nodeCounts.Max(),
                 "kron jacobian graph must be scale-invariant: " + string.Join(",", nodeCounts));

    var hessCounts = new List<long>();
    foreach (int n in sizes) {
      MX Ax = MX.sym("Ax", n, n);
      MX Bx = MX.sym("Bx", 2, 2);
      MX Kx = ca.kron(Ax, Bx);
      MX obj = 0.5 * ca.sumsqr(Kx);
      MX x = ca.vertcat(ca.vec(Ax), ca.vec(Bx));
      var (H, _) = ca.hessian(obj, x);
      hessCounts.Add(fun("H", new[] { Ax, Bx }, new[] { H }).n_nodes());
    }
    assertEquals(hessCounts.Min(), hessCounts.Max(),
                 "kron hessian graph must be scale-invariant: " + string.Join(",", hessCounts));

    MX A2 = MX.sym("A2", 3, 4);
    MX B2 = MX.sym("B2", 2, 5);
    MX M2 = MX.sym("M2", ca.kron(A2, B2).sparsity());
    var cases = new List<MX[]> { new[] { ca.kron_contract(M2, B2, true), M2, B2 },
                                 new[] { ca.kron_contract(M2, A2, false), M2, A2 } };
    foreach (MX[] cs in cases) {
      MX[] args = cs[1..3];
      Function f = fun("f", args, new[] { cs[0] });
      DM.rng(11);
      var inputs = args.Select(arg => DM.rand(arg.sparsity())).ToList();
      checkfunctionLight(f, f.expand(), inputs);
    }

    MX A3 = MX.sym("A3", 3, 3);
    MX B3 = MX.sym("B3", 2, 2);
    MX M3 = MX.sym("M3", ca.kron(A3, B3).sparsity());
    MX Y = ca.kron_contract(M3, B3, true);
    MX dM = MX.sym("dM", M3.sparsity());
    MX dB = MX.sym("dB", B3.sparsity());
    MX fwd = ca.jtimes(Y, ca.vertcat(ca.vec(M3), ca.vec(B3)), ca.vertcat(ca.vec(dM), ca.vec(dB)));
    MX expected = ca.kron_contract(dM, B3, true) + ca.kron_contract(M3, dB, true);
    Function g = fun("g", new[] { M3, B3, dM, dB }, new[] { fwd - expected });
    DM.rng(3);
    DM diff = ev(g, DM.rand(M3.sparsity()), DM.rand(B3.sparsity()), DM.rand(M3.sparsity()), DM.rand(B3.sparsity()));
    assertTrue((double)ca.norm_inf(diff) < 1e-14);

    DM Av = DM.rand(A2.sparsity()), Bv = DM.rand(B2.sparsity());
    DM Mv = DM.rand(M2.sparsity());
    Function fi = fun("fi", new[] { M2, B2 }, new[] { ca.kron_contract(M2, B2, true) });
    Function fo = fun("fo", new[] { M2, A2 }, new[] { ca.kron_contract(M2, A2, false) });
    checkCodegen(fi, dms(Mv, Bv));
    checkCodegen(fo, dms(Mv, Av));
    checkSerialize(fi, dms(Mv, Bv));
    checkSerialize(fo, dms(Mv, Av));
  }

  public void test_project() {
    MX x = MX.sym("x", Sparsity.lower(3));
    MX y = ca.project(x, Sparsity.lower(3).T());
    Function f = fun("f", new[] { x }, new[] { y });
    DM fIn = ca.DM(f.sparsity_in(0), ca.vertcat(1, 2, 3, 4, 5, 6));
    DM fOut = ev(f, fIn);
    checkarray(fOut, ca.DM(new double[,] { { 1, 0, 0 }, { 0, 4, 0 }, { 0, 0, 6 } }));
    checkarray(DM.ones(f.sparsity_out(0)), DM.ones(Sparsity.lower(3).T()));
    checkCodegen(f, dms(ca.DM(new double[,] { { 1, 0, 0 }, { 0, 4, 0 }, { 0, 0, 6 } })));
    checkSerialize(f, dms(ca.DM(new double[,] { { 1, 0, 0 }, { 0, 4, 0 }, { 0, 0, 6 } })));
  }
}

}

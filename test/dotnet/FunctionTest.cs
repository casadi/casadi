using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;
using System.Text.RegularExpressions;
using System.Threading;

namespace CasADi.Tests {

/// <summary>Port of test/java FunctionTest (itself a port of test/python/function.py)</summary>
public class FunctionTest : CasadiTestCase {

  // ---------- helpers ----------

  static Dictionary<string, object> opts(params object[] kv) {
    var m = new Dictionary<string, object>();
    for (int i = 0; i < kv.Length; i += 2) m[(string)kv[i]] = kv[i + 1];
    return m;
  }

  static Function jacobianOld(Function f, long i, long j) {
    var outs = new List<string> { "jac:" + f.name_out(j) + ":" + f.name_in(i) };
    outs.AddRange(f.name_out());
    return f.factory(f.name() + "_jac", f.name_in(), outs);
  }

  static Function hessianOld(Function f, long i, long j) {
    var outs = new List<string> {
      "hess:" + f.name_out(j) + ":" + f.name_in(i) + ":" + f.name_in(i),
      "grad:" + f.name_out(j) + ":" + f.name_in(i)
    };
    outs.AddRange(f.name_out());
    return f.factory(f.name() + "_hess", f.name_in(), outs);
  }

  static Function toSXFun(Function f) {
    List<SX> ins_ = f.sx_in();
    return ca.Function("f", ins_, f.call(ins_));
  }

  /// <summary>Inputs as DM: numbers, double[] (column), double[,] or DM</summary>
  static List<DM> ins(params object[] a) {
    var r = new List<DM>();
    foreach (object o in a) {
      switch (o) {
        case DM d: r.Add(d); break;
        case double[] v: r.Add(v); break;
        case double[,] m: r.Add(m); break;
        default: r.Add(ca.DM(Convert.ToDouble(o))); break;
      }
    }
    return r;
  }

  /// <summary>As np.random.random: uniform nonzeros on sp</summary>
  DM rnd(Sparsity sp) {
    double[] v = new double[sp.nnz()];
    for (int k = 0; k < v.Length; k++) v[k] = random.NextDouble();
    return ca.DM(sp, v);
  }

  DM rnd(long n, long m) => rnd(Sparsity.dense(n, m));

  static double[] linspace(double a, double b, int n) {
    double[] r = new double[n];
    for (int i = 0; i < n; i++) r[i] = n == 1 ? a : a + (b - a) * i / (n - 1);
    return r;
  }

  static List<T> cat<T>(params IEnumerable<T>[] a) => a.SelectMany(x => x).ToList();

  static string captureStderr(Action r) => captureBoth(r)[1];

  void checkfunctionLight(Function f, Function g, List<DM> inputs) => checkfun(f, g, inputs, false, false, false, 1e-9);

  /// <summary>checkfunction with the jacobian/forward/reverse parts individually selectable</summary>
  void checkfun(Function f, Function g, List<DM> inputs, bool jac, bool fwd, bool adj, double tol) {
    assertEquals(g.n_in(), f.n_in(), "n_in");
    assertEquals(g.n_out(), f.n_out(), "n_out");
    List<DM> a = f.call(inputs), b = g.call(inputs);
    for (int i = 0; i < a.Count; i++) checkarray(a[i], b[i], tol, "output " + i);
    if (jac) {
      List<DM> va = f.jacobian().call(cat(inputs, a)), vb = g.jacobian().call(cat(inputs, b));
      for (int i = 0; i < va.Count; i++) checkarray(va[i], vb[i], tol, "jacobian block " + i);
    }
    if (fwd) {
      var s = cat(inputs, a);
      for (long i = 0; i < f.n_in(); i++) s.Add(rnd(f.sparsity_in(i)));
      List<DM> fa = f.forward(1).call(s), fb = g.forward(1).call(s);
      for (int i = 0; i < fa.Count; i++) checkarray(fa[i], fb[i], tol, "forward " + i);
    }
    if (adj) {
      var s = cat(inputs, a);
      for (long i = 0; i < f.n_out(); i++) s.Add(rnd(f.sparsity_out(i)));
      List<DM> ra = f.reverse(1).call(s), rb = g.reverse(1).call(s);
      for (int i = 0; i < ra.Count; i++) checkarray(ra[i], rb[i], tol, "reverse " + i);
    }
  }

  void checkSerialize(Function f, List<DM> inputs) {
    Function f2 = Function.deserialize(f.serialize(opts("debug", true)));
    List<DM> a = f.call(inputs), b = f2.call(inputs);
    for (int i = 0; i < f.n_out(); i++) checkarray(a[i], b[i], 1e-15, "serialize output " + i);
  }

  void checkSerialize(Function f, IDictionary<string, DM> inputs) {
    Function f2 = Function.deserialize(f.serialize(opts("debug", true)));
    Dictionary<string, DM> a = f.call(inputs), b = f2.call(inputs);
    assertTrue(f.name_out().SequenceEqual(f2.name_out()), "name_out");
    foreach (string k in f.name_out()) checkarray(a[k], b[k], 1e-15, "serialize output " + k);
  }

  static int cgCounter = 0;

  static void sh(string cmd, string dir) {
    var psi = new ProcessStartInfo("sh") { WorkingDirectory = dir, UseShellExecute = false };
    psi.ArgumentList.Add("-c");
    psi.ArgumentList.Add(cmd);
    using var p = Process.Start(psi);
    p.WaitForExit();
    assertEquals(0, p.ExitCode, "command failed: " + cmd);
  }

  Function checkCodegen(Function f, List<DM> inputs) =>
    checkCodegen(f, inputs, opts(), opts(), false, "c89", false);

  /// <summary>As helpers.check_codegen: only active with RunSlow (needs gcc); returns the loaded external or null</summary>
  Function checkCodegen(Function f, List<DM> inputs, IDictionary<string, object> opts_, IDictionary<string, object> extOpts,
      bool withJacSparsity, string std, bool main) {
    if (!RunSlow) return null;
    string lib = GlobalOptions.getCasadiPath(), inc = GlobalOptions.getCasadiIncludePath();
    if (lib == "" || inc == "")
      Check.Skip("casadi include path not set (build tree)");
    string dir = Path.Combine(Path.GetTempPath(), "casadi_dotnet_cg_" + Guid.NewGuid().ToString("N"));
    Directory.CreateDirectory(dir);
    string name = "codegen_" + (cgCounter++) + "_" + random.Next();
    var o = new Dictionary<string, object>(opts_);
    if (main) o["main"] = true;
    var cg = ca.CodeGenerator(name, o);
    cg.add(f, withJacSparsity);
    cg.generate(dir + "/");
    string flags = "gcc -pedantic -std=" + std + " -fPIC -Wall -Werror -Wextra -I" + inc + " -I" + inc
        + "/highs -Wno-dangling-pointer -Wno-unknown-warning-option -Wno-unknown-pragmas -Wno-long-long"
        + " -Wno-unused-parameter -O3 " + name + ".c";
    string libs = " -L" + lib + " -Wl,-rpath," + lib + " -Wl,-rpath,.";
    sh(flags + " -shared -o " + name + ".so" + libs, dir);
    List<DM> refs = f.call(inputs);
    if (main) {
      sh(flags + " -o " + name + libs + " -lm", dir);
      string inF = Path.Combine(dir, f.name() + "_in.txt"), outF = Path.Combine(dir, f.name() + "_out.txt");
      f.generate_in(inF, inputs);
      sh("./" + name + " " + f.name() + " < '" + inF + "' > '" + outF + "'", dir);
      List<DM> outputs = f.generate_out(outF);
      for (int i = 0; i < f.n_out(); i++) checkarray(refs[i], outputs[i], 1e-14, "main output " + i);
    }
    Function f2 = ca.external(f.name(), dir + "/" + name + ".so", extOpts);
    List<DM> res = f2.call(inputs);
    for (int i = 0; i < f.n_out(); i++) checkarray(refs[i], res[i], 1e-14, "codegen output " + i);
    return f2;
  }

  void checkCodegen(Function f, IDictionary<string, DM> inputs, bool main) =>
    checkCodegen(f, f.convert_in(inputs), opts(), opts(), false, "c89", main);

  // ---------- tests ----------

  public void test_call_empty() {
    SX xs = SX.sym("x", 2);
    var fsx = ca.Function("fsx", new[] { xs, ca.SX() }, new[] { xs });
    MX x = MX.sym("x", 2);
    var fmx1 = ca.Function("fmx1", new[] { x, ca.MX() }, new[] { x });
    var fmx2 = ca.Function("fmx2", new[] { x, ca.MX() }, new[] { x });

    foreach (Function f in new[] { fsx, fmx1, fmx2 }) {
      f.Call(0, 0);
      MX X = MX.sym("X", 2);
      MX F = f.Call(X, ca.MX());
      var g = ca.Function("g", new[] { X }, new[] { F });
      g.Call(0);
    }

    xs = SX.sym("x", 2);
    fsx = ca.Function("fsx", new[] { xs }, new[] { xs, ca.SX() });
    x = MX.sym("x", 2);
    fmx1 = ca.Function("fmx1", new[] { x }, new[] { x, ca.MX() });
    fmx2 = ca.Function("fmx2", new[] { x }, new[] { x, ca.MX() });

    foreach (Function f in new[] { fsx, fmx1 }) {
      f.call(new DM[] { 0 });
      MX X = MX.sym("X", 2);
      List<MX> F = f.call(new[] { X });
      var g = ca.Function("g", new[] { X }, F);
      g.Call(0);
    }
  }

  public void test_MX_funSeed() {
    message("MX_funSeed");
    MX x1 = MX.sym("x", 2), y1 = MX.sym("y"), x2 = MX.sym("x", 2), y2 = MX.sym("y");
    var p = ca.Function("p", new[] { x1, y1, x2, y2 }, new[] { ca.sin(x1) + y1, ca.sin(x2) + y2 });
    DM n1 = ca.vertcat(4, 5), n2 = ca.vertcat(5, 7);
    double N1 = 3, N2 = 8;
    List<DM> outs = p.call(new DM[] { n1, N1, n2, N2 });
    checkarray(ca.sin(n1) + N1, outs[0], "output");
    checkarray(ca.sin(n2) + N2, outs[1], "output");
  }

  public void test_segfault() {
    var f = ca.Function();
    Check.Throws<Exception>(() => f.stats());
  }

  public void test_issue304() {
    message("regression test for #304");
    SX x = SX.sym("x");
    var f = ca.Function("f", new[] { x }, new[] { ca.pow(x, 2.0), ca.pow(x, 3.0) });
    MX X = MX.sym("X");
    List<MX> z = f.call(new[] { X });
    ca.Function("g", new[] { X }, z).expand();
  }

  public void test_jacobian() {
    SX x = SX.sym("x", 3, 1), y = SX.sym("y", 2, 1);
    var f = ca.Function("f", new[] { x, y }, new[] { ca.pow(x, 2.0), y, x * y[0] });
    Function g = jacobianOld(f, 0, 0);
    assertEquals(f.n_in(), g.n_in());
    assertEquals(f.n_out() + 1, g.n_out());
  }

  public void test_xfunction() {
    SX x = SX.sym("x", 3, 1), y = SX.sym("y", 2, 1);
    var f = ca.Function("f", new[] { x, y }, new[] { ca.pow(x, 2.0), y, x * y[0] });
    MX X = MX.sym("x", 3, 1), Y = MX.sym("y", 2, 1);
    var F = ca.Function("F", new[] { X, Y }, new[] { ca.pow(X, 2.0), Y, X * Y[0] });
    checkfun(f, F, ins(new double[] { 0.1, 0.7, 1.3 }, new double[] { 7.1, 2.9 }), true, false, false, 1e-9);
  }

  static void jacobiansTest(Sparsity sp) {
    SX x = SX.sym("x", sp.size2());
    Sparsity sp2 = ca.jacobian(ca.mtimes(DM.ones(sp), x), x).sparsity();
    assertTrue(sp.row().SequenceEqual(sp2.row()), "row");
    assertTrue(sp.colind().SequenceEqual(sp2.colind()), "colind");
  }

  public void test_jacobians() {
    SX x = SX.sym("x");
    assertEquals(0L, ca.jacobian(5.0, x).nnz());

    for (int i = 0; i < 5; i++) {
      jacobiansTest(Sparsity.lower(i));
      jacobiansTest(Sparsity.lower(i).T());
      jacobiansTest(Sparsity.dense(i, i));
      jacobiansTest(Sparsity.diag(i));
    }

    foreach (int i in new[] { 63, 64, 65, 127, 128, 129 }) {
      Sparsity d = Sparsity.diag(i);
      jacobiansTest(d);
      jacobiansTest(d.plus(Sparsity.rowcol(new long[] { 0 }, new long[] { 5 }, i, i)));
      Sparsity b = Sparsity.band(i, -1).plus(Sparsity.band(i, 1));
      jacobiansTest(b.plus(Sparsity.rowcol(new long[] { 0 }, new long[] { 5 }, i, i)));
    }

    DM m = DM.ones(Sparsity.diag(129));
    m[0..50, 0] = 1.0;
    m[60.., 0] = 1.0;
    m[6..9, 6] = 1.0;
    m[9, 9..12] = 1.0;

    Sparsity sp = m[.., 0..120].sparsity();
    jacobiansTest(sp);

    m = DM.ones(Sparsity.diag(64));
    m[0..50, 0] = 1.0;
    m[60.., 0] = 1.0;

    foreach (DM mm in new[] { m.T(), m })
      foreach (int k in new[] { 40, 20 }) {
        sp = mm[0..k, ..].sparsity();
        jacobiansTest(sp);
        jacobiansTest(sp.T());
      }

    foreach (int i in new[] { 63, 64, 65, 127, 128, 129 }) {
      jacobiansTest(Sparsity.lower(i));
      jacobiansTest(Sparsity.lower(i).T());
    }

    int[] sizes = RunSlow ? new[] { 63, 64, 65, 127, 128, 129 } : new[] { 63, 64, 65 };
    foreach (int n in sizes)
      foreach (int mm in sizes) {
        sp = Sparsity.dense(n, mm);
        jacobiansTest(sp);
        DM I = DM.ones(sp);
        for (int i = 0; i < n; i++)
          for (int j = 0; j < mm; j++)
            if (random.NextDouble() < 0.5) I[i, j] = 0.0;
        I = ca.sparsify(I);
        Sparsity spHoles = I.sparsity();
        jacobiansTest(spHoles);
        var z = ca.DM(spHoles.size1(), spHoles.size2());
        var v = new List<DM>();
        for (int r = 0; r < 5; r++) {
          var h = Enumerable.Repeat(z, 5).ToList();
          h[r] = I;
          v.Add(ca.horzcat(h));
        }
        jacobiansTest(ca.vertcat(v).sparsity());
      }
  }

  static void hessiansTest(Sparsity sp) {
    SX x = SX.sym("x", sp.size2());
    assertTrue(sp.is_equal(sp.T()));
    var f = ca.Function("f", new[] { x }, new[] { ca.mtimes(new[] { x.T(), ca.SX(DM.ones(sp)), x }) });
    Function J = hessianOld(f, 0, 0);
    Sparsity sp2 = J.sparsity_out(0);
    assertTrue(sp.row().SequenceEqual(sp2.row()), "row");
    assertTrue(sp.colind().SequenceEqual(sp2.colind()), "colind");
  }

  public void test_hessians() {
    double[][,] As = {
      new double[,] {{1, 1, 0, 0, 0, 0}, {1, 1, 1, 0, 1, 1}, {0, 1, 1, 1, 0, 0}, {0, 0, 1, 1, 0, 1}, {0, 1, 0, 0, 1, 0}, {0, 1, 0, 1, 0, 1}},
      new double[,] {{1, 0, 0, 0, 0, 0}, {0, 1, 1, 0, 1, 1}, {0, 1, 1, 1, 0, 0}, {0, 0, 1, 1, 0, 1}, {0, 1, 0, 0, 1, 0}, {0, 1, 0, 1, 0, 1}},
      new double[,] {{1, 0, 0, 0, 0, 0}, {0, 1, 0, 0, 1, 1}, {0, 0, 1, 1, 0, 0}, {0, 0, 1, 1, 0, 1}, {0, 1, 0, 0, 1, 0}, {0, 1, 0, 1, 0, 1}},
      new double[,] {{0, 0, 0, 0, 0, 0}, {0, 1, 0, 0, 1, 1}, {0, 0, 1, 1, 0, 0}, {0, 0, 1, 1, 0, 1}, {0, 1, 0, 0, 1, 0}, {0, 1, 0, 1, 0, 1}},
      new double[,] {{0, 0, 0, 0, 0, 0}, {0, 1, 0, 0, 1, 0}, {0, 0, 1, 1, 0, 0}, {0, 0, 1, 1, 0, 1}, {0, 1, 0, 0, 1, 0}, {0, 0, 0, 1, 0, 1}}};
    foreach (double[,] A in As) hessiansTest(ca.sparsify((DM)A).sparsity());

    foreach (int i in new[] { 63, 64, 65, 100, 127, 128, 129 }) {
      Sparsity d = Sparsity.diag(i);
      hessiansTest(d);
      Sparsity rc = Sparsity.rowcol(new long[] { 0 }, new long[] { 5 }, i, i).plus(Sparsity.rowcol(new long[] { 5 }, new long[] { 0 }, i, i));
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
          if (random.NextDouble() < 0.5) {
            I[ii, jj] = 0.0;
            I[jj, ii] = 0.0;
          }
      I = ca.sparsify(I);
      Sparsity spHoles = I.sparsity();
      hessiansTest(spHoles);
      var z = ca.DM(spHoles.size1(), spHoles.size2());
      var v = new List<DM>();
      for (int r = 0; r < 5; r++) {
        var h = Enumerable.Repeat(z, 5).ToList();
        h[r] = I;
        v.Add(ca.horzcat(h));
      }
      hessiansTest(ca.vertcat(v).sparsity());
    }
  }

  public void test_customIO() {
    SX x = SX.sym("x");
    var f = ca.Function("f", new[] { x }, new[] { x * x, x }, new[] { "i0" }, new[] { "foo", "bar" });
    Dictionary<string, DM> ret = f.call(new Dictionary<string, DM> { { "i0", 12 } });
    checkarray(ca.vertcat(144), ret["foo"]);
    checkarray(ca.vertcat(12), ret["bar"]);
    // f_out["baz"]: python NameError on an undefined variable, nothing to port
    Dictionary<string, SX> rets = f.call(new Dictionary<string, SX> { { "i0", ca.SX(12) } });
    checkarray(ca.DM(rets["foo"]), ca.vertcat(144));
    checkarray(ca.DM(rets["bar"]), ca.vertcat(12));
    assertFalse(rets.ContainsKey("baz"));
  }

  public void test_derivative_simplifications() {
    int n = 1;
    SX x = SX.sym("x", n);
    double[] r = new double[n];
    for (int i = 0; i < n; i++) r[i] = i;
    var M = ca.Function("M", new[] { x }, new[] { ca.mtimes(x - (DM)r, x.T()) });
    MX P = MX.sym("P", n, n), X = MX.sym("X", n);
    MX M_X = M.Call(X);
    var Pf = ca.Function("P", new[] { X, P }, new[] { ca.mtimes(M_X, P) });
    Function P_P = jacobianOld(Pf, 1, 0);
    assertFalse(P_P.ToString().Contains("derivative"));
  }

  public void test_issue1464() {
    int n = 6, N = 9;
    SX x = SX.sym("x", n), u = SX.sym("u");
    var rk4 = ca.Function("f", new[] { x, u }, new[] { x + u });
    long[] split = { 0, (N + 1) * n, (N + 1) * n + N };
    {
      List<SX> g = new List<SX>(), g2 = new List<SX>();
      SX V = SX.sym("V", (N + 1) * n + N);
      List<SX> VV = ca.vertsplit(V, split);
      List<SX> VXk = ca.vertsplit(VV[0], n), VUk = ca.vertsplit(VV[1], 1);
      for (int k = 0; k < N; k++) {
        SX xf = rk4.Call(VXk[k], VUk[k]);
        List<SX> xfp = ca.vertsplit(xf, n / 2), vp = ca.vertsplit(VXk[k + 1], n / 2);
        g.Add(xfp[0] - vp[0]);
        g.Add(xfp[1] - vp[1]);
        g2.Add(xf - VXk[k + 1]);
      }
      for (int i = 0; i < 2; i++) {
        var f = ca.Function("nlp", new[] { V }, new[] { ca.vertcat(g) }, opts("ad_weight_sp", i));
        assertEquals(162L, f.jac_sparsity(0, 0).nnz());
        var f2 = ca.Function("nlp", new[] { V }, new[] { ca.vertcat(g2) }, opts("ad_weight_sp", i));
        assertEquals(162L, f2.jac_sparsity(0, 0).nnz());
      }
    }
    {
      List<MX> g = new List<MX>(), g2 = new List<MX>();
      MX V = MX.sym("V", (N + 1) * n + N);
      List<MX> VV = ca.vertsplit(V, split);
      List<MX> VXk = ca.vertsplit(VV[0], n), VUk = ca.vertsplit(VV[1], 1);
      for (int k = 0; k < N; k++) {
        MX xf = rk4.Call(VXk[k], VUk[k]);
        List<MX> xfp = ca.vertsplit(xf, n / 2), vp = ca.vertsplit(VXk[k + 1], n / 2);
        g.Add(xfp[0] - vp[0]);
        g.Add(xfp[1] - vp[1]);
        g2.Add(xf - VXk[k + 1]);
      }
      for (int i = 0; i < 2; i++) {
        var f = ca.Function("nlp", new[] { V }, new[] { ca.vertcat(g) }, opts("ad_weight_sp", i));
        assertEquals(162L, f.jac_sparsity(0, 0).nnz());
        var f2 = ca.Function("nlp", new[] { V }, new[] { ca.vertcat(g2) }, opts("ad_weight_sp", i));
        assertEquals(162L, f2.jac_sparsity(0, 0).nnz());
      }
    }
  }

  class SquareCallback : Callback {
    public SquareCallback(string name, IDictionary<string, object> opts) { construct(name, opts); }
    public override List<DM> eval(IEnumerable<DM> argin) => new List<DM> { ca.pow(argin.First(), 2.0) };
  }

  public void test_callback() {
    var foo = new SquareCallback("my_f", opts());
    MX x = MX.sym("x");
    MX y = foo.Call(x);
    var f = ca.Function("f", new[] { x }, new[] { y });
    checkarray(f.Call(5), 25);
  }

  public void test_callback_buffer() {
    Check.Skip("BINDING GAP: Callback.eval_buffer (FunctionBuffer) has no C# equivalent");
  }

  class ThrowingCallback : Callback {
    public ThrowingCallback(string name) { construct(name, new Dictionary<string, object>()); }
    public override List<DM> eval(IEnumerable<DM> argin) => throw new Exception("foobar");
  }

  public void test_callback_errors() {
    var foo = new ThrowingCallback("my_f");
    MX x = MX.sym("x");
    MX y = foo.Call(x);
    var f = ca.Function("f", new[] { x }, new[] { y });
    try {
      f.Call(3);
    } catch (Exception e) when (!(e is AssertFailed)) {
      assertTrue(e.Message.Contains("foobar"), e.Message);
    }
  }

  class MySVD : Callback {
    readonly long n, m, k;
    readonly bool fd;

    public MySVD(string name, long n, long m, bool fd) {
      this.n = n;
      this.m = m;
      this.fd = fd;
      k = Math.Min(n, m);
      construct(name);
    }

    public override long get_n_in() => 1;
    public override long get_n_out() => 3;
    public override Sparsity get_sparsity_in(long i) => Sparsity.dense(n, m);

    public override Sparsity get_sparsity_out(long i) {
      if (i == 0) return Sparsity.dense(n, k);
      else if (i == 1) return Sparsity.dense(k);
      else if (i == 2) return Sparsity.dense(k, m);
      throw new ArgumentException("No such index");
    }
    // eval (numpy svd) is never evaluated by the test
  }

  class Convolution : Callback {
    readonly DM A;
    readonly long n, m;
    readonly bool flagTransp;
    readonly List<Convolution> data = new List<Convolution>();

    public Convolution(DM A, bool flagTransp) {
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
      var o = new Convolution(A, !flagTransp);
      data.Add(o);
      return o;
    }

    public override Sparsity get_sparsity_in(long i) => Sparsity.dense(m, 1);
    public override Sparsity get_sparsity_out(long i) => Sparsity.dense(n, 1);

    public override List<DM> eval(IEnumerable<DM> arg) {
      DM x = arg.First();
      return new List<DM> { flagTransp ? ca.mtimes(A.T(), x) : ca.mtimes(A, x) };
    }

    public override long get_n_in() => 1;
    public override long get_n_out() => 1;
    public override bool has_forward(long nfwd) => true;
    public override bool has_reverse(long nadj) => true;

    public override Function get_forward(long nfwd, string name, IEnumerable<string> inames, IEnumerable<string> onames, IDictionary<string, object> opts) {
      Function mf = map(nfwd, "serial");
      MX fwdIn = MX.sym("x", m, nfwd);
      MX fwdOut = mf.Call(fwdIn);
      MX arg = MX.sym("x", m, 1);
      MX dummy = MX.sym("x", ca.Sparsity(n, 1));
      return ca.Function(name, new[] { arg, dummy, fwdIn }, new[] { fwdOut }, inames, onames, opts);
    }

    public override Function get_reverse(long nadj, string name, IEnumerable<string> inames, IEnumerable<string> onames, IDictionary<string, object> opts) {
      Function mf = getTranspose().map(nadj, "serial");
      MX adjIn = MX.sym("x", n, nadj);
      MX adjOut = mf.Call(adjIn);
      MX arg = MX.sym("x", m, 1);
      var dummy = ca.MX(n, 1);
      return ca.Function(name, new[] { arg, dummy, adjIn }, new[] { adjOut }, inames, onames, opts);
    }

    public override bool has_jacobian() => false;
  }

  public void test_callback_convolution() {
    int n = 3, mm = 3;
    var foo = new MySVD("foo", n, mm, true);
    DM.rng(1);
    DM.rand(n, mm);
    MX X = MX.sym("x", n, mm);
    foo.call(new[] { X });
    DM A = DM.rand(5, 5);
    var conv = new Convolution(A, false);
    conv.forward(1);
    conv.jacobian();
  }

  /// <summary>The fun, X, Y, Z, V of the map tests</summary>
  class MapSetup {
    public readonly Function fun;
    public readonly List<MX> X = new List<MX>(), Y = new List<MX>(), Z = new List<MX>(), V = new List<MX>();

    public MapSetup(int n) {
      SX x = SX.sym("x"), y = SX.sym("y", 2), z = SX.sym("z", 2, 2), v = SX.sym("z", Sparsity.upper(3));
      fun = ca.Function("f", new[] { x, y, z, v }, new[] { ca.mtimes(z, y) + x, ca.sin(y * x).T(), v / x });
      for (int i = 0; i < n; i++) {
        X.Add(MX.sym("x"));
        Y.Add(MX.sym("y", 2));
        Z.Add(MX.sym("z", 2, 2));
        V.Add(MX.sym("z", Sparsity.upper(3)));
      }
    }

    public List<MX> all() => cat(X, Y, Z, V);

    public List<MX> hcats() => new List<MX> { ca.horzcat(X), ca.horzcat(Y), ca.horzcat(Z), ca.horzcat(V) };
  }

  List<DM> rndLike(List<MX> xs) => xs.Select(x => rnd(x.sparsity())).ToList();

  static List<MX> sinAll(List<MX> xs) => xs.Select(x => ca.sin(x)).ToList();

  public void test_mapdict() {
    SX x = SX.sym("x"), y = SX.sym("y", 2), z = SX.sym("z", 2, 2), v = SX.sym("z", Sparsity.upper(3));
    var dict = new Dictionary<string, SX> {
      { "x", x }, { "y", y }, { "z", z }, { "v", v },
      { "I", ca.mtimes(z, y) + x }, { "II", ca.sin(y * x).T() }, { "III", v / x }
    };
    var fun = ca.Function("f", dict, new[] { "x", "y", "z", "v" }, new[] { "I", "II", "III" });

    var s = new MapSetup(2);
    var arg = new Dictionary<string, MX> {
      { "x", ca.horzcat(s.X) }, { "y", ca.horzcat(s.Y) }, { "z", ca.horzcat(s.Z) }, { "v", ca.horzcat(s.V) }
    };
    Dictionary<string, MX> res = fun.map(2).call(arg);
    List<MX> res2 = fun.map(2).call(s.hcats());

    var F = ca.Function("F", s.all(), res2);
    var F2 = ca.Function("F", s.all(), new[] { res["I"], res["II"], res["III"] });
    checkfunctionLight(F, F2, rndLike(s.all()));
  }

  void mapNode(string[] parallelizations, bool light) {
    var s = new MapSetup(2);
    foreach (string parallelization in parallelizations) {
      Console.WriteLine(parallelization);
      List<MX> res = s.fun.map(2, parallelization).call(s.hcats());
      var F = ca.Function("F", s.all(), sinAll(res));
      var resref = new List<List<MX>>();
      for (int i = 0; i < s.fun.n_out(); i++) resref.Add(new List<MX>());
      for (int k = 0; k < 2; k++) {
        List<MX> e = sinAll(s.fun.call(new[] { s.X[k], s.Y[k], s.Z[k], s.V[k] }));
        for (int i = 0; i < e.Count; i++) resref[i].Add(e[i]);
      }
      var Fref = ca.Function("F", s.all(), resref.Select(r => ca.horzcat(r)).ToList());
      List<DM> inputs = rndLike(s.all());
      foreach (Function f in new[] { F, F.expand("expand_" + F.name()) }) {
        if (light) checkfunctionLight(f, Fref, inputs);
        else checkfunction(f, Fref, inputs);
      }
    }
  }

  public void test_map_node() {
    mapNode(RunSlow ? new[] { "serial", "openmp", "unroll", "inline", "thread" } : new[] { "serial" }, false);
  }

  public void test_map_node_light() {
    mapNode(new[] { "serial", "openmp", "unroll", "inline", "thread" }, true);
  }

  public void test_map_node_n_threads() {
    var s = new MapSetup(1);
    Function fun = s.fun;
    List<DM> X_ = new List<DM>(), Y_ = new List<DM>(), Z_ = new List<DM>(), V_ = new List<DM>();
    for (int i = 0; i < 10; i++) {
      X_.Add(rnd(fun.sparsity_in(0)));
      Y_.Add(rnd(fun.sparsity_in(1)));
      Z_.Add(rnd(fun.sparsity_in(2)));
      V_.Add(rnd(fun.sparsity_in(3)));
    }
    Console.WriteLine(fun.map(3, "thread", 2));
    long[][] cases = { new long[] { 2, 1 }, new long[] { 3, 1 }, new long[] { 3, 2 }, new long[] { 4, 2 }, new long[] { 4, 5 } };
    foreach (long[] c in cases) {
      int k = (int)c[0];
      var inputs = new List<DM> { ca.horzcat(X_.Take(k)), ca.horzcat(Y_.Take(k)), ca.horzcat(Z_.Take(k)), ca.horzcat(V_.Take(k)) };
      checkfunctionLight(fun.map(k, "thread", c[1]), fun.map(k), inputs);
    }
  }

  // test_map_fenv: python-only (switches the FPU rounding mode through ctypes/libm)

  public void test_mapsum() {
    var s = new MapSetup(2);
    Function fun = s.fun;
    foreach (List<MX> Z_alt in new[] { s.Z, new List<MX> { ca.MX(), ca.MX(), ca.MX() } }) {
      foreach (string parallelization in new[] { "serial", "openmp", "unroll", "thread" }) {
        List<MX> res = fun.mapsum(new[] { ca.horzcat(s.X), ca.horzcat(s.Y), ca.horzcat(Z_alt), ca.horzcat(s.V) }, parallelization);
        foreach (int ad_weight_sp in new[] { 0, 1 }) {
          var F = ca.Function("F", s.all(), sinAll(res), opts("ad_weight", 0, "ad_weight_sp", ad_weight_sp));
          var resref = new List<MX>();
          for (int i = 0; i < fun.n_out(); i++) resref.Add(ca.MX(0.0));
          for (int k = 0; k < 2; k++) {
            List<MX> e = fun.call(new[] { s.X[k], s.Y[k], Z_alt[k], s.V[k] });
            for (int i = 0; i < e.Count; i++) resref[i] = resref[i] + e[i];
          }
          var Fref = ca.Function("F", s.all(), sinAll(resref));
          List<DM> inputs = rndLike(s.all());
          if (parallelization != "thread") checkCodegen(F, inputs);
          foreach (Function f in new[] { F, toSXFun(F) }) checkfunction(f, Fref, inputs);
        }
      }
    }
  }

  public void test_mapsum2() {
    var s = new MapSetup(2);
    Function fun = s.fun;
    int n = 2;
    MX Z = MX.sym("z", 2, 2), V = MX.sym("z", Sparsity.upper(3));
    foreach (string parallelization in new[] { "true_map_sum", "serial", "openmp", "unroll", "thread" })
      foreach (int ad_weight_sp in new[] { 0, 1 })
        foreach (int ad_weight in new[] { 0, 1 }) {
          var o = opts("ad_weight_sp", ad_weight_sp, "ad_weight", ad_weight);
          Function F;
          if (parallelization == "true_map_sum")
            F = fun.map(n, new[] { false, false, true, true }, new[] { true, false, false }, o);
          else F = fun.map("map", parallelization, n, new long[] { 2, 3 }, new long[] { 0 }, o);

          MX acc = ca.MX(0.0);
          List<MX> bl = new List<MX>(), cl = new List<MX>();
          for (int k = 0; k < n; k++) {
            List<MX> abc = fun.call(new[] { s.X[k], s.Y[k], Z, V });
            acc = acc + abc[0];
            bl.Add(abc[1]);
            cl.Add(abc[2]);
          }
          var Fref = ca.Function("F", new[] { ca.horzcat(s.X), ca.horzcat(s.Y), Z, V }, new[] { acc, ca.horzcat(bl), ca.horzcat(cl) });
          var inputs = new List<DM> { ca.horzcat(rndLike(s.X)), ca.horzcat(rndLike(s.Y)), rnd(Z.sparsity()), rnd(V.sparsity()) };
          if (parallelization != "thread") checkCodegen(F, inputs);
          foreach (Function f in new[] { F, toSXFun(F) }) checkfunction(f, Fref, inputs);
          checkSerialize(F, inputs);
        }
  }

  public void test_issue4274() {
    MX x = MX.sym("x");
    var c = ca.Function("x", new[] { x }, new[] { ca.pow(x.printme(4), 2.0) }, new[] { "x" }, new[] { "out" });
    var localInputs = new Dictionary<string, double> { { "x", 5.0 } };
    Function F = c.map(2, "thread");
    var wideInputs = new Dictionary<string, DM>();
    foreach (var e in localInputs)
      wideInputs[e.Key] = ca.horzcat(new List<DM> { e.Value, e.Value * 1.2 });
    Console.WriteLine(string.Join(", ", wideInputs.Select(e => e.Key + "=" + e.Value)));
    F.call(wideInputs);
  }

  public void test_repmatnode() {
    MX x = MX.sym("x", 2);
    MX y = ca.sin(ca.repmat(ca.pow(x, 2.0), 1, 3));
    MX z = MX.sym("y", 2, 2);
    var F = ca.Function("f", new[] { x, z }, new[] { ca.sum2(ca.sum1(y)) });
    SX xs = SX.sym("x", 2);
    SX ys = ca.sin(ca.repmat(ca.pow(xs, 2.0), 1, 3));
    SX zs = SX.sym("y", 2, 2);
    var Fref = ca.Function("f", new[] { xs, zs }, new[] { ca.sum2(ca.sum1(ys)) });
    var inputs = new List<DM> { ca.vertcat(1, 7), ca.DM(new double[,] { { 3, 0 }, { 2, 4 } }) };
    checkCodegen(F, inputs);
    checkfunction(F, Fref, inputs);
  }

  public void test_repsumnode() {
    MX x = MX.sym("x", 2), z = MX.sym("y", 2, 2);
    var F = ca.Function("f", new[] { x, z },
        new[] { ca.sin(ca.repsum(ca.pow(x, 2.0).T(), 1, 2)), (ca.cos(ca.pow(x, 2.0)) * 2.0 * x).T() });
    SX xs = SX.sym("x", 2), zs = SX.sym("y", 2, 2);
    var Fref = ca.Function("f", new[] { xs, zs },
        new[] { ca.sin(ca.repsum(ca.pow(xs, 2.0).T(), 1, 2)), (ca.cos(ca.pow(xs, 2.0)) * 2.0 * xs).T() });
    var inputs = new List<DM> { ca.vertcat(1, 7), ca.DM(new double[,] { { 3, 0 }, { 2, 4 } }) };
    checkCodegen(F, inputs);
    checkfunction(F, Fref, inputs);
  }

  public void test_unknown_options() {
    SX x = SX.sym("x");
    // python spells the constructor SXFunction (undefined, so trivially raises); the intent is ported
    Check.Throws<Exception>(() => ca.Function("f", new[] { x }, new[] { x }, opts("fooo", false)));
    Check.Throws<Exception>(() => ca.Function("f", new[] { x }, new[] { x }, opts("ad_weight", "foo")));
  }

  public void test_unknown_options_stringvector() {
    Check.Assume(hasNlpsol("ipopt"), "no nlpsol ipopt");
    SX x = SX.sym("x");
    ca.nlpsol("mysolver", "ipopt", new Dictionary<string, SX> { { "x", x }, { "f", ca.pow(x, 2.0) } }, opts("monitor", new[] { "nlp_f" }));
    string err = captureStderr(() => ca.nlpsol("mysolver", "ipopt", new Dictionary<string, SX> { { "x", x }, { "f", ca.pow(x, 2.0) } },
        opts("monitor", new[] { "abc" })));
    assertTrue(err.Contains("Ignoring monitor 'abc'. Available functions: nlp_f"), err);
  }

  public void test_mapaccum() {
    SX x = SX.sym("x", 2), y = SX.sym("y"), z = SX.sym("z", 2, 2), v = SX.sym("v", Sparsity.upper(3));
    var fun = ca.Function("f", new[] { x, y, z, v }, new[] { ca.mtimes(z, x) + y, ca.sin(y * x).T(), v / y });
    int n = 2;
    MX X = MX.sym("x", x.sparsity());
    List<MX> Y = new List<MX>(), Z = new List<MX>(), V = new List<MX>();
    for (int i = 0; i < n; i++) {
      Y.Add(MX.sym("y", y.sparsity()));
      Z.Add(MX.sym("z", z.sparsity()));
      V.Add(MX.sym("v", v.sparsity()));
    }
    DM X_ = rnd(x.sparsity());
    List<DM> Y_ = rndLike(Y), Z_ = rndLike(Z), V_ = rndLike(V);

    for (int ad_weight = 0; ad_weight < 2; ad_weight++)
      for (int ad_weight_sp = 0; ad_weight_sp < 2; ad_weight_sp++) {
        Function F = fun.mapaccum("map", n, new long[] { 0 }, new long[] { 0 }, opts("ad_weight_sp", ad_weight_sp, "ad_weight", ad_weight));
        F.forward(2);
        MX XP = X;
        List<MX> Y0s = new List<MX>(), Y1s = new List<MX>(), Xps = new List<MX>();
        for (int k = 0; k < n; k++) {
          List<MX> r = fun.call(new[] { XP, Y[k], Z[k], V[k] });
          XP = r[0];
          Y0s.Add(r[1]);
          Y1s.Add(r[2]);
          Xps.Add(XP);
        }
        var Fref = ca.Function("f", new[] { X, ca.horzcat(Y), ca.horzcat(Z), ca.horzcat(V) },
            new[] { ca.horzcat(Xps), ca.horzcat(Y0s), ca.horzcat(Y1s) });
        var inputs = new List<DM> { X_, ca.horzcat(Y_), ca.horzcat(Z_), ca.horzcat(V_) };
        foreach (Function f in new[] { F, toSXFun(F) }) {
          checkfunction(f, Fref, inputs);
          checkCodegen(f, inputs);
        }
      }

    fun = ca.Function("f", new[] { y, x, z, v },
        new[] { ca.mtimes(z, x) + y + ca.pow(ca.trace(v), 2.0), ca.sin(y * x).T(), v / y });
    for (int ad_weight = 0; ad_weight < 2; ad_weight++)
      for (int ad_weight_sp = 0; ad_weight_sp < 2; ad_weight_sp++) {
        Function F = fun.mapaccum("map", n, new long[] { 1, 3 }, new long[] { 0, 2 }, opts("ad_weight_sp", ad_weight_sp, "ad_weight", ad_weight));
        MX XP = X, VP = V[0];
        List<MX> Y0s = new List<MX>(), Xps = new List<MX>(), Vps = new List<MX>();
        for (int k = 0; k < n; k++) {
          List<MX> r = fun.call(new[] { Y[k], XP, Z[k], VP });
          XP = r[0];
          VP = r[2];
          Y0s.Add(r[1]);
          Xps.Add(XP);
          Vps.Add(VP);
        }
        var Fref = ca.Function("f", new[] { ca.horzcat(Y), X, ca.horzcat(Z), V[0] },
            new[] { ca.horzcat(Xps), ca.horzcat(Y0s), ca.horzcat(Vps) });
        var inputs = new List<DM> { ca.horzcat(Y_), X_, ca.horzcat(Z_), V_[0] };
        foreach (Function f in new[] { F, toSXFun(F) }) {
          checkfunction(f, Fref, inputs);
          checkCodegen(f, inputs);
        }
      }
  }

  public void test_mapaccum_schemes() {
    SX x = SX.sym("x", 2), y = SX.sym("y"), z = SX.sym("z", 2, 2), v = SX.sym("v", Sparsity.upper(3));
    var fun = ca.Function("f", new[] { y, z, x, v }, new[] { ca.mtimes(z, x) + y, ca.sin(y * x).T(), v / y },
        new[] { "y", "z", "x", "v" }, new[] { "out0", "out1", "out2" });
    int n = 2;
    Function F = fun.mapaccum("map", n, new long[] { 2 }, new long[] { 0 });
    List<string> schemeInFun = fun.name_in(), schemeOutFun = fun.name_out();
    List<string> schemeInF = F.name_in(), schemeOutF = F.name_out();
    for (int rep = 0; rep < 2; rep++) {
      // python: assertTrue(len(a), len(b)) only checks the first is nonzero
      assertTrue(schemeInFun.Count > 0);
      assertTrue(schemeOutFun.Count > 0);
      for (int i = 0; i < Math.Min(schemeInFun.Count, schemeInF.Count); i++) assertEquals(schemeInFun[i], schemeInF[i]);
      for (int i = 0; i < Math.Min(schemeOutFun.Count, schemeOutF.Count); i++) assertEquals(schemeOutFun[i], schemeOutF[i]);
      fun = ca.Function("f", new[] { x, y, z, v }, new[] { ca.mtimes(z, x) + y, ca.sin(y * x).T(), v / y },
          new[] { "x", "y", "z", "v" }, new[] { "out0", "out1", "out2" });
      F = fun.mapaccum("map", n);
    }
  }

  public void test_depends_on() {
    SX x = SX.sym("x");
    SX y = ca.pow(x, 2.0);
    string s = "";
    try {
      ca.depends_on(x, y);
    } catch (SwigException e) {
      s = e.Message;
    }
    assertTrue(s.Contains("not symbolic"), s);
    try {
      ca.Function("f", new[] { y }, new[] { x });
    } catch (SwigException e) {
      s = e.Message;
    }
    assertTrue(s.Contains("not symbolic"), s);
  }

  static bool same(DM a, double b) => Math.Abs((double)a - b) < 1e-8;

  void checkCodegenSerialize(Function f, List<DM> inputs) {
    Function f2 = checkCodegen(f, inputs);
    if (f2 != null) checkSerialize(f2, inputs);
  }

  public void test_1d_interpolant() {
    double[][] grid = { new double[] { 0, 1, 1.5, 2, 3 } };
    double[] values = { 0, 1, 2, 5, 3 };
    Function F = ca.interpolant("F", "linear", grid, values);
    double[][] pairs = { new[] { 3.4, 3 - 0.4 * 2 }, new[] { 2.4, 5 - 0.4 * 2 }, new[] { 1.6, 2 + 3 * 0.1 / 0.5 },
      new[] { 1.4, 1 + 0.4 / 0.5 }, new[] { 0.4, 0.4 }, new[] { -.6, -0.6 } };
    MX X = MX.sym("x");
    ca.Function("F", new[] { X }, new[] { F.Call(X) });
    foreach (double[] p in pairs) {
      assertTrue(same(F.Call(p[0]), p[1]));
      checkCodegenSerialize(F, ins(p[0]));
      checkSerialize(F, ins(p[0]));
    }
    Console.WriteLine(string.Join(", ", F.stats().Select(e => e.Key + ": " + e.Value)));

    X = MX.sym("x");
    var J = ca.Function("F", new[] { X }, new[] { ca.jacobian(F.Call(X), X) });
    double[][] jpairs = { new[] { 3.4, -2.0 }, new[] { 2.4, -2.0 }, new[] { 1.6, 6.0 }, new[] { 1.4, 2.0 }, new[] { 0.4, 1.0 },
      new[] { -.6, 1.0 }, new[] { 1.0, 2.0 }, new[] { 0.99, 1.0 } };
    foreach (double[] p in jpairs) {
      assertTrue(same(J.Call(p[0]), p[1]));
      checkCodegen(J, ins(p[0]));
      checkSerialize(J, ins(p[0]));
    }
  }

  public void test_2d_interpolant() {
    double[][] grid = { new double[] { 0, 1, 4, 5 }, new double[] { 0, 2, 3 } };
    double[] values = { 0, 1, 8, 3, 10, -11, 12, 13, 20, 31, -42, 53 };
    Function F = ca.interpolant("F", "linear", grid, values);
    double a0 = -11 + 0.4 * (31 + 11), a1 = 12 + 0.4 * (-42 - 12);
    double[][] pairs = { new double[] { 1, 2, -11 }, new double[] { 1, 3, 31 }, new double[] { 4, 2, 12 }, new double[] { 4, 3, -42 },
      new double[] { 1, 2.4, a0 }, new double[] { 4, 2.4, a1 },
      new double[] { 3, 2, -11 + 2.0 / 3 * (12 + 11) }, new double[] { 3, 3, 31 + 2.0 / 3 * (-42 - 31) },
      new double[] { 3, 2.4, a0 + 2.0 / 3 * (a1 - a0) } };
    foreach (double[] p in pairs) {
      DM a = ca.vertcat(p[0], p[1]);
      checkarray(F.Call(a), p[2]);
      checkCodegen(F, new List<DM> { a });
    }

    MX X = MX.sym("x", 2);
    var J = ca.Function("J", new[] { X }, new[] { ca.jacobian(F.Call(X), X) });
    double jx0 = (12 + 11) / 3.0, jx1 = (-42 - 31) / 3.0, jx2 = (13 - 12), jx3 = (53 + 42);
    double jy0 = 31 + 11, jy1 = -42 - 12;
    double[][] jpairs = { new[] { 1, 2, jx0, jy0 }, new[] { 1, 3, jx1, jy0 }, new[] { 4, 2, jx2, jy1 }, new[] { 4, 3, jx3, jy1 },
      new[] { 1, 2.4, jx0 + (jx1 - jx0) * 0.4, 31 + 11 }, new[] { 4, 2.4, jx2 + (jx3 - jx2) * 0.4, -42 - 12 },
      new[] { 3, 2, jx0, jy0 + (jy1 - jy0) * 2.0 / 3 }, new[] { 3, 3, jx1, jy0 + (jy1 - jy0) * 2.0 / 3 },
      new[] { 3, 2.4, jx0 + (jx1 - jx0) * 0.4, jy0 + (jy1 - jy0) * 2.0 / 3 } };
    foreach (double[] p in jpairs) {
      DM a = ca.vertcat(p[0], p[1]);
      checkarray(J.Call(a).T(), ca.vertcat(p[2], p[3]));
      checkCodegen(J, new List<DM> { a });
      checkSerialize(J, new List<DM> { a });
    }
  }

  public void test_1d_interpolant_uniform() {
    double[][] grid = { new double[] { 0, 1, 2 } };
    double[] values = { 0, 1, 2 };
    foreach (string mode in new[] { "linear", "exact", "binary" }) {
      Function F = ca.interpolant("F", "linear", grid, values, opts("lookup_mode", new[] { mode }));
      assertTrue(same(F.Call(2.4), 2.4));
      assertTrue(same(F.Call(1.4), 1.4));
      assertTrue(same(F.Call(0), 0));
      assertTrue(same(F.Call(1), 1));
      assertTrue(same(F.Call(2), 2));
      assertTrue(same(F.Call(6), 6));
      assertTrue(same(F.Call(0.4), 0.4));
      assertTrue(same(F.Call(-.6), -.6));
      ca.interpolant("F", "linear", new[] { linspace(0, 1, 7) }, new double[] { 0, 1, 2, 3, 4, 5, 6 }, opts("lookup_mode", new[] { "exact" }));
    }
    grid = new[] { new double[] { 2, 4, 6 } };
    values = new double[] { 10, 7, 1 };
    foreach (string mode in new[] { "linear", "exact", "binary" }) {
      Function F = ca.interpolant("F", "linear", grid, values, opts("lookup_mode", new[] { mode }));
      assertTrue(same(F.Call(1), 11.5));
      assertTrue(same(F.Call(2), 10));
      assertTrue(same(F.Call(3), 8.5));
      assertTrue(same(F.Call(4), 7));
      assertTrue(same(F.Call(5), 4));
      assertTrue(same(F.Call(6), 1));
      assertTrue(same(F.Call(7), -2));
      ca.interpolant("F", "linear", new[] { linspace(0, 1, 7) }, new double[] { 0, 1, 2, 3, 4, 5, 6 }, opts("lookup_mode", new[] { "exact" }));
    }
  }

  public void test_2d_interpolant_uniform() {
    double[][] grid = { new double[] { 0, 1, 2 }, new double[] { 0, 1, 2 } };
    double[] values = { 0, 1, 2, 10, 11, 12, 20, 21, 22 };
    foreach (string mode in new[] { "linear", "exact", "binary" }) {
      Function F = ca.interpolant("F", "linear", grid, values, opts("lookup_mode", new[] { mode, mode }));
      double[][] cases = { new[] { 2.4, 0.5, 7.4 }, new[] { 1.4, 0.5, 6.4 }, new[] { 0.4, 0.5, 5.4 }, new[] { 1, 0.5, 6 },
        new double[] { 1, 0, 1 }, new double[] { 0, 0, 0 }, new[] { 0.4, 1, 10.4 },
        new[] { -.6, 0.5, 4.4 }, new[] { -.6, 1.5, 14.4 }, new[] { -.6, 2.5, 24.4 }, new[] { -.6, 3.5, 34.4 } };
      foreach (double[] c in cases) assertTrue(same(F.Call(new double[] { c[0], c[1] }), c[2]));
    }
  }

  /// <summary>Multilinear interpolation on a rectilinear grid, data in column-major (first dimension fastest)</summary>
  static double interpn(double[][] grid, double[] d, double[] x) {
    int n = grid.Length;
    int[] lo = new int[n];
    double[] w = new double[n];
    for (int i = 0; i < n; i++) {
      double[] g = grid[i];
      int k = 0;
      while (k < g.Length - 2 && x[i] >= g[k + 1]) k++;
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
        stride *= grid[i].Length;
      }
      r += wc * d[idx];
    }
    return r;
  }

  public void test_nd_linear() {
    int N = 4;
    for (int n_dim = 1; n_dim < 5; n_dim++) {
      double[][] grid = new double[n_dim][];
      double[] x0 = new double[n_dim];
      for (int i = 0; i < n_dim; i++) {
        double[] g = new double[N];
        for (int k = 0; k < N; k++) g[k] = random.NextDouble();
        Array.Sort(g);
        grid[i] = g;
        x0[i] = g.Average();
      }
      double[] d = new double[(int)Math.Pow(N, n_dim)];
      for (int k = 0; k < d.Length; k++) d[k] = random.NextDouble();
      double xref = interpn(grid, d, x0);
      Function F = ca.interpolant("F", "linear", grid, d);
      checkarray(F.Call(x0), xref);
    }
  }

  public void test_2d_bspline() {
    double[][] d_knots = { linspace(0, 1, 5), linspace(0, 1, 6) };
    double[] d_flat = rnd(5 * 6, 1).nonzeros().ToArray();
    Function LUT = ca.interpolant("name", "bspline", d_knots, d_flat);
    jacobianOld(LUT, 0, 0);
    hessianOld(LUT, 0, 0);
    checkCodegen(LUT, new List<DM> { ca.vertcat(0.2, 0.3) });
    checkSerialize(LUT, new List<DM> { ca.vertcat(0.2, 0.3) });
    // comparison against scipy.interpolate.RectBivariateSpline: python-only (scipy reference)
  }

  public void test_1d_bspline_call_fun() {
    double[][] d_knots = { linspace(0, 1, 5) };
    double[] d_flat = rnd(5, 1).nonzeros().ToArray();
    Function LUT = ca.interpolant("name", "bspline", d_knots, d_flat);
    SX xs = SX.sym("x");
    var f0 = ca.Function("f", new[] { xs }, new[] { LUT.Call(ca.sin(xs)) * 2.0 });
    MX xm = MX.sym("x");
    var f1 = ca.Function("f", new[] { xm }, new[] { LUT.Call(ca.sin(xm)) * 2.0 });
    checkfunction(f0, f1, ins(0.32));
  }

  public void test_1d_bspline() {
    double[][] d_knots = { linspace(0, 1, 5) };
    double[] d_flat = rnd(5, 1).nonzeros().ToArray();
    Function LUT = ca.interpolant("name", "bspline", d_knots, d_flat);
    checkCodegen(LUT, ins(0.2));
    jacobianOld(LUT, 0, 0);
    hessianOld(LUT, 0, 0);
    // comparison against scipy.interpolate.InterpolatedUnivariateSpline: python-only (scipy reference)
  }

  public void test_bspline_even_degree() {
    double[] x = linspace(0, 3, 12), y = new double[12];
    for (int i = 0; i < 12; i++) y[i] = Math.Sin(x[i]);
    Function LUT = ca.interpolant("LUT", "bspline", new[] { x }, y, opts("degree", new long[] { 2 }));
    // comparison against scipy make_interp_spline: python-only (scipy reference)
    checkCodegen(LUT, ins(0.2));
    foreach (long k in new long[] { 4, 6 }) {
      LUT = ca.interpolant("LUT", "bspline", new[] { x }, y, opts("degree", new long[] { k }));
      checkarray(LUT.Call(((DM)x).T()), ((DM)y).T(), 1e-10, "");
      checkarray(LUT.Call(1.234), ca.DM(Math.Sin(1.234)), 1e-4, "");
    }
    double[] xg = linspace(0, 1, 6), yg = linspace(0, 2, 7);
    double[] V = new double[6 * 7];
    for (int j = 0; j < 7; j++)
      for (int i = 0; i < 6; i++) V[i + 6 * j] = Math.Cos(xg[i]) * yg[j];
    LUT = ca.interpolant("LUT", "bspline", new[] { xg, yg }, V, opts("degree", new long[] { 2, 3 }));
    foreach (int[] ij in new[] { new[] { 0, 0 }, new[] { 2, 3 }, new[] { 5, 6 }, new[] { 4, 1 } })
      checkarray(LUT.Call(ca.vertcat(xg[ij[0]], yg[ij[1]])), ca.DM(V[ij[0] + 6 * ij[1]]), 1e-10, "");
  }

  /// <summary>sin(x+3*y)</summary>
  class SinFun : Callback {
    public override long get_n_in() => 2;
    public override long get_n_out() => 1;

    public override List<DM> eval(IEnumerable<DM> arg) {
      var a = arg.ToList();
      DM x = a[0], y = a[1];
      DM z0 = y * 3.0;
      DM z1 = x + z0;
      return new List<DM> { ca.sin(z1) };
    }
  }

  class JacFun : SinFun {
    public JacFun() { construct("Fun", new Dictionary<string, object>()); }
    public override bool has_forward(long nfwd) => false;
    public override bool has_reverse(long nadj) => false;
    public override bool has_jacobian() => true;

    public override Function get_jacobian(string name, IEnumerable<string> inames, IEnumerable<string> onames, IDictionary<string, object> opts) {
      SX x = SX.sym("x"), y = SX.sym("y");
      SX out_g = SX.sym("out_g", ca.Sparsity(1, 1));
      SX c = ca.cos(x + y * 3.0);
      return ca.Function(name, new[] { x, y, out_g }, new[] { c, c * 3.0 }, inames, onames, opts);
    }
  }

  public void test_Callback_Jacobian() {
    MX x = MX.sym("x"), y = MX.sym("y");
    List<DM> numInputs = ins(0.2, 0.7);
    var g = ca.Function("g", new[] { x, y }, new[] { ca.sin(x + y * 3.0) });
    var f = new JacFun();
    checkfunction(f, g, numInputs);
  }

  // test_Callback_errors: python-only (overrides get_sparsity_in with a wrong signature returning an int)

  class Unstringable : Exception {
    public override string Message => throw new Exception("boom in Message");
    public override string ToString() => throw new Exception("boom in ToString");
  }

  class UnstringableFun : Callback {
    public UnstringableFun() { construct("Fun", new Dictionary<string, object>()); }
    public override List<DM> eval(IEnumerable<DM> arg) => throw new Unstringable();
  }

  public void test_Callback_unstringable_exception() {
    var f = new UnstringableFun();
    Check.Throws<SwigException>(() => f.Call(0));
  }

  class PlainSinFun : SinFun {
    public PlainSinFun() { construct("Fun", new Dictionary<string, object>()); }
  }

  public void test_Callback() {
    MX x = MX.sym("x"), y = MX.sym("y");
    List<DM> numInputs = ins(0.2, 0.7);
    var g = ca.Function("g", new[] { x, y }, new[] { ca.sin(x + y * 3.0) });
    foreach (bool indirect in new[] { true, false }) {
      var cb = new PlainSinFun();
      Function f = indirect ? ca.Function("f", new[] { x, y }, new[] { cb.Call(x, y) }) : cb;
      checkfunctionLight(f, g, numInputs);
      // f.gradient(): python AttributeError (Function has no gradient method), nothing to port
      Check.Throws<Exception>(() => jacobianOld(f, 0, 0));
      Check.Throws<Exception>(() => f.forward(1));
      Check.Throws<Exception>(() => f.reverse(1));
    }
  }

  class DimFun : Callback {
    readonly List<DM> ret;

    public DimFun(List<DM> ret) {
      this.ret = ret;
      construct("Fun");
    }

    public override long get_n_in() => 2;
    public override long get_n_out() => 1;
    public override List<DM> eval(IEnumerable<DM> arg) => ret;
  }

  public void test_Callback_dimcheck() {
    var f = new DimFun(new List<DM> { ca.DM(2), ca.DM(1) });
    string s = "";
    try {
      f.Call(2);
    } catch (SwigException e) {
      s = e.Message;
    }
    assertTrue(s.Contains("Incorrect number of inputs"), s);
    f = new DimFun(new List<DM> { ca.DM(2), ca.DM(1) });
    s = "";
    try {
      f.Call(2, 3);
    } catch (SwigException e) {
      s = e.Message;
    }
    assertTrue(s.Contains("Expected 1 output"), s);
    s = "";
    f = new DimFun(new List<DM> { DM.zeros(2, 2) });
    try {
      f.Call(2, 3);
    } catch (SwigException e) {
      s = e.Message;
    }
    assertTrue(s.Contains("Shape mismatch"), s);
  }

  class ForwardFun : Callback {
    public ForwardFun(string name) { construct(name, new Dictionary<string, object> { { "verbose", true } }); }
    public override long get_n_in() => 2 + 1 + 2;
    public override long get_n_out() => 1;

    public override List<DM> eval(IEnumerable<DM> args) {
      var arg = args.ToList();
      DM x = arg[0], y = arg[1];
      DM z1 = x + y * 3.0;
      var ret = new List<DM>();
      for (int i = 3; i < arg.Count; i += 2) {
        DM dx = arg[i], dy = arg[i + 1];
        DM dz1 = dx + dy * 3.0;
        ret.Add(ca.cos(z1) * dz1);
      }
      return ret;
    }
  }

  class BackwardFun : Callback {
    public BackwardFun(string name) { construct(name, new Dictionary<string, object> { { "verbose", true } }); }
    public override long get_n_in() => 2 + 1 + 1;
    public override long get_n_out() => 2;

    public override List<DM> eval(IEnumerable<DM> args) {
      var arg = args.ToList();
      DM x = arg[0], y = arg[1];
      DM z1 = x + y * 3.0;
      var ret = new List<DM>();
      for (int i = 3; i < arg.Count; i++) {
        DM bz2 = arg[i];
        DM bz1 = bz2 * ca.cos(z1);
        DM bx = bz1, bz0 = bz1;
        DM by = bz0 * 3.0;
        ret.Add(bx);
        ret.Add(by);
      }
      return ret;
    }
  }

  class SensFun : SinFun {
    readonly bool hasFwd, hasAdj;
    ForwardFun cbFwd;
    BackwardFun cbRev;

    public SensFun(bool hasFwd, bool hasAdj, IDictionary<string, object> opts) {
      this.hasFwd = hasFwd;
      this.hasAdj = hasAdj;
      construct("Fun", opts);
    }

    public override bool has_forward(long nfwd) => hasFwd && nfwd == 1;

    public override Function get_forward(long nfwd, string name, IEnumerable<string> inames, IEnumerable<string> onames, IDictionary<string, object> opts) {
      assertEquals(1L, nfwd);
      cbFwd = new ForwardFun(name);
      return cbFwd;
    }

    public override bool has_reverse(long nadj) => hasAdj && nadj == 1;

    public override Function get_reverse(long nadj, string name, IEnumerable<string> inames, IEnumerable<string> onames, IDictionary<string, object> opts) {
      assertEquals(1L, nadj);
      cbRev = new BackwardFun(name);
      return cbRev;
    }
  }

  public void test_Callback_sens() {
    MX x = MX.sym("x"), y = MX.sym("y");
    List<DM> numInputs = ins(0.2, 0.7);
    var g = ca.Function("g", new[] { x, y }, new[] { ca.sin(x + y * 3.0) });
    var keep = new List<SensFun>();
    foreach (bool indirect in new[] { true, false }) {
      bool[][] variants = { new[] { true, true }, new[] { true, false }, new[] { false, true }, new[] { true, false }, new[] { false, true } };
      foreach (bool[] v in variants) {
        var cb = new SensFun(v[0], v[1], opts("verbose", true));
        keep.Add(cb);
        Function f = indirect ? ca.Function("f", new[] { x, y }, new[] { cb.Call(x, y) }) : cb;
        checkfunction(f, g, numInputs);
      }
    }
  }

  class JacSpFun : Callback {
    readonly List<DM> calls;
    readonly double h;
    readonly bool withJacSparsity;

    public JacSpFun(List<DM> calls, double h, bool withJacSparsity) {
      this.calls = calls;
      this.h = h;
      this.withJacSparsity = withJacSparsity;
      construct("Fun", opts("enable_fd", true, "fd_method", "forward", "fd_options", opts("h", h, "h_iter", false)));
    }

    public override long get_n_in() => 1;
    public override long get_n_out() => 1;
    public override Sparsity get_sparsity_in(long i) => Sparsity.dense(2, 1);
    public override Sparsity get_sparsity_out(long i) => Sparsity.dense(2, 1);

    public override List<DM> eval(IEnumerable<DM> arg) {
      DM x = arg.First();
      calls.Add((x - ca.vertcat(5, 7)) / h);
      return new List<DM> { ca.pow(x, 2.0) };
    }

    public override bool has_forward(long nfwd) => false;
    public override bool has_reverse(long nadj) => false;
    public override bool has_jac_sparsity(long oind, long iind) => withJacSparsity;

    public override Sparsity get_jac_sparsity(long oind, long iind, bool symmetric) {
      assertEquals(0L, oind);
      assertEquals(0L, iind);
      return Sparsity.diag(2);
    }
  }

  public void test_callback_jacobian_sparsity() {
    MX x = MX.sym("x", 2);
    double h = 1e-7;
    foreach (bool withJacSparsity in new[] { true, false }) {
      var calls = new List<DM>();
      var f = new JacSpFun(calls, h, withJacSparsity);
      MX y = ca.jacobian(f.Call(x), x);
      var F = ca.Function("F", new[] { x }, new[] { y });
      DM J = null;
      captureStdout(() => J = F.Call(new double[] { 5, 7 }));
      DM callsM = ca.horzcat(calls);
      DM J_ref = ca.DM(new double[,] { { 5 * 2, 0 }, { 0, 7 * 2 } });
      if (withJacSparsity) {
        checkarray(callsM, ca.DM(new double[,] { { 0, 0 }, { 1, 1 } }).T(), 1e-5, "");
        J_ref = ca.sparsify(J_ref);
      } else {
        checkarray(callsM, ca.DM(new double[,] { { 0, 0 }, { 1, 0 }, { 0, 1 } }).T(), 1e-5, "");
      }
      checkarray(J, J_ref, 1e-5, "");
      assertTrue(J.sparsity().is_equal(J_ref.sparsity()));
    }
  }

  public void test_common_specific_options() {
    Check.Assume(hasNlpsol("ipopt"), "no nlpsol ipopt");
    SX x = SX.sym("x");
    var nlp = new Dictionary<string, SX> { { "x", x }, { "f", ca.pow(x, 4.0) } };
    var x0 = new Dictionary<string, DM> { { "x0", 0.1 } };

    string o = captureStdout(() => ca.nlpsol("solver", "ipopt", nlp));
    assertFalse(o.Contains("nlp_f"));
    Function solver = ca.nlpsol("solver", "ipopt", nlp, opts("common_options", opts("final_options", opts("print_in", true))));
    o = captureStdout(() => solver.call(x0));
    Console.WriteLine(o);
    assertTrue(o.Contains("Function nlp_f"));
    assertTrue(o.Contains("Function nlp_grad_f"));
    Function solver2 = ca.nlpsol("solver", "ipopt", nlp, opts("specific_options", opts("nlp_f", opts("final_options", opts("print_in", true)))));
    o = captureStdout(() => solver2.call(x0));
    assertTrue(o.Contains("Function nlp_f"));
    assertFalse(o.Contains("Function nlp_grad_f"));
    Function solver3 = ca.nlpsol("solver", "ipopt", nlp, opts("common_options", opts("final_options", opts("print_in", true)),
        "specific_options", opts("nlp_f", opts("final_options", opts("print_in", false)))));
    o = captureStdout(() => solver3.call(x0));
    assertFalse(o.Contains("Function nlp_f"));
    assertTrue(o.Contains("Function nlp_grad_f"));
    Function solver4 = ca.nlpsol("solver", "ipopt", nlp, opts("common_options", opts("final_options", opts("print_in", false)),
        "specific_options", opts("nlp_f", opts("final_options", opts("print_in", true)))));
    o = captureStdout(() => solver4.call(x0));
    assertTrue(o.Contains("Function nlp_f"));
    assertFalse(o.Contains("Function nlp_grad_f"));

    string err = captureStderr(() => ca.nlpsol("solver", "ipopt", nlp));
    assertEquals(0, err.Length, err);
    err = captureStderr(() => ca.nlpsol("solver", "ipopt", nlp, opts("specific_options", opts("nlp_foo", opts("helper_options", opts("verbose", true))))));
    assertTrue(("Ignoring" + err).Length > 0);
    assertTrue(err.Contains("nlp_g"), err);
    Check.Throws<Exception>(() => ca.nlpsol("solver", "ipopt", nlp, opts("specific_options", opts("nlp_foo", 3))));
  }

  public void test_oracle_options() {
    Check.Assume(hasNlpsol("ipopt"), "no nlpsol ipopt");
    DM.set_precision(16);
    int N = 5;
    MX x = MX.sym("x", N, 1);
    var options = opts(
      "ipopt", opts("hessian_approximation", "limited-memory"),
      "oracle_options", opts("enable_fd", true, "enable_forward", false, "enable_reverse", false, "print_in", true, "fd_method", "central"));
    Function solver = ca.nlpsol("solver", "ipopt", new Dictionary<string, MX> { { "x", x }, { "f", ca.dot(x, x) } }, options);
    string o = captureStdout(() => solver.call(new Dictionary<string, DM> { { "x0", 1 } }));
    DM.set_precision(6);
    assertTrue(o.Contains("[1, 1, 1, 1, 1]"));
  }

  public void test_expm() {
    Check.Assume(ca.has_expm("slicot"), "no expm slicot");
    Check.Assume(hasIntegrator("cvodes"), "no integrator cvodes");
    int n = 4;
    DM Anum = rnd(n, n);
    MX A = MX.sym("A", n, n), t = MX.sym("t");
    var fr = ca.Function("fr", new[] { A, t }, new[] { expmIntegrator(A * t, n) });
    var f = ca.Function("f", new[] { A, t }, new[] { ca.expm(A * t) });
    checkfunction(fr, f, ins(Anum, 1.1), 1e-8);
    fr = ca.Function("fr", new[] { t }, new[] { expmIntegrator(ca.MX(Anum) * t, n) });
    f = ca.Function("f", new[] { t }, new[] { ca.expm_const(Anum, t) });
    checkfunction(fr, f, ins(1.1), 1e-8);

    MX JA = ca.jacobian(ca.expm(A * t), A), Jt = ca.jacobian(ca.expm(A * t), t);
    assertEquals((long)(n * n * n * n), JA.nnz());
    assertEquals((long)(n * n), Jt.nnz());
    JA = ca.jacobian(ca.expm_const(A, t), A);
    Jt = ca.jacobian(ca.expm_const(A, t), t);
    assertEquals(0L, JA.nnz());
    assertEquals((long)(n * n), Jt.nnz());
  }

  static MX expmIntegrator(MX Ain, int n) {
    MX x = MX.sym("x", n), As = MX.sym("A", n, n);
    var dae = new Dictionary<string, MX> { { "x", x }, { "p", ca.vec(As) }, { "ode", ca.mtimes(As, x) } };
    Function intg = ca.integrator("intg", "cvodes", dae, opts("reltol", 1e-14, "abstol", 1e-14));
    Function Intg = intg.map("identity", "serial", n, new long[] { 1 }, new long[] { });
    Dictionary<string, MX> o = Intg.call(new Dictionary<string, MX> { { "x0", ca.MX(DM.eye(n)) }, { "p", ca.vec(As) } });
    var expmF = ca.Function("expm", new[] { As }, new[] { o["xf"] });
    return expmF.Call(Ain);
  }

  public void test_conditional() {
    MX x = MX.sym("x", 2, 2);
    MX sp1 = MX.sym("y", Sparsity.lower(2)), sp2 = MX.sym("z", Sparsity.diag(2));
    var f1 = ca.Function("f", new[] { sp2, x }, new[] { ca.pow(x, 2.0), x * sp2 });
    var f2 = ca.Function("f", new[] { sp1, x }, new[] { ca.pow(x, 2.0) * 2.0, ca.sin(sp1) });
    var f3 = ca.Function("f", new[] { sp1, sp2 }, new[] { sp1 * sp2, sp1 + sp2 });
    Function F = Function.conditional("test", new[] { f1, f2 }, f3);
    Function Fsx = F.expand();
    double[,] A = rnd(2, 2).ToArray(), B = rnd(2, 2).ToArray();
    for (int i = -1; i < 3; i++) {
      checkfunction(F, Fsx, ins(i, A, B));
      checkCodegen(F, ins(i, A, B));
    }
  }

  // test_max_num_dir (first definition): shadowed in python by the second definition of the same name

  public void test_max_num_dir() {
    MX x = MX.sym("x");
    var f = ca.Function("ffff", new[] { x }, new[] { ca.sin(x) });
    Function fm = f.mapaccum("mapaccum", 100, 1, opts("base", 4));
    var c = ca.CodeGenerator("me");
    c.add(fm);
    string code = c.dump();
    assertTrue(code.Contains("ffff_acc4_acc4_acc4"));
  }

  public void test_codegen_with_jac_sparsity() {
    Check.Assume(RunSlow, "slow");
    MX x = MX.sym("x", 3), y = MX.sym("y", 3, 3);
    var f = ca.Function("f", new[] { x, y }, new[] { ca.pow(x, 3.0), ca.mtimes(y, x) });
    var c = ca.CodeGenerator("me");
    c.add(f, true);
    Function fJ = f.jacobian();
    DM.rng(0);
    DM x0 = DM.rand(x.sparsity()), y0 = DM.rand(y.sparsity());
    Function F2 = checkCodegen(f, new List<DM> { x0, y0 }, opts(), opts("enable_fd", true), true, "c89", false);
    Function FJ = F2.jacobian();
    for (int i = 0; i < fJ.n_out(); i++) assertTrue(fJ.sparsity_out(i).is_equal(FJ.sparsity_out(i)));
    F2 = checkCodegen(f, new List<DM> { x0, y0 }, opts(), opts("enable_fd", true), false, "c89", false);
    FJ = F2.jacobian();
    for (int i = 0; i < fJ.n_out(); i++) assertTrue(FJ.sparsity_out(i).is_dense());
  }

  /// <summary>column-major flat data of two random 5x6 tables, interleaved as in the python tests</summary>
  double[][] twoTables() {
    double[] d0 = rnd(30, 1).nonzeros().ToArray(), d1 = rnd(30, 1).nonzeros().ToArray(), dflat = new double[60];
    for (int i = 0; i < 30; i++) {
      dflat[2 * i] = d0[i];
      dflat[2 * i + 1] = d1[i];
    }
    return new[] { d0, d1, dflat };
  }

  void multiout(string solver, bool all) {
    double[][] d_knots = { linspace(0, 1, 5), linspace(0, 1, 6) };
    double[][] t = twoTables();
    Function LUT0 = ca.interpolant("name", solver, d_knots, t[0]);
    Function LUT1 = ca.interpolant("name", solver, d_knots, t[1]);
    double[] d_flat = t[2];
    Function LUTf = ca.interpolant("name", solver, d_knots, d_flat);
    MX x = MX.sym("x"), y = MX.sym("y");
    MX xy = ca.vertcat(x, y);
    var LUT_sep = ca.Function("f", new[] { x, y }, new[] { ca.vertcat(LUT0.Call(xy), LUT1.Call(xy)) });
    var LUT = ca.Function("f", new[] { x, y }, new[] { LUTf.Call(xy) });
    List<DM> inputs = ins(0.2, 0.333);
    checkfunction(LUT, LUT_sep, inputs);
    checkCodegen(LUT, inputs);
    checkSerialize(LUT, inputs);

    Function LUT_param = ca.interpolant("name", solver, d_knots, 2);
    var f = ca.Function("LUTp", new[] { x, y }, new[] { LUT_param.Call(xy, (DM)d_flat) });
    checkfunction(LUT, f, inputs);
    checkCodegen(f, inputs);
    checkSerialize(f, inputs);
    if (!all) return;

    DM d_knots_cat = d_knots[0].Concat(d_knots[1]).ToArray();

    LUT_param = ca.interpolant("name", solver, new long[] { 5, 6 }, 2);
    f = ca.Function("LUTp", new[] { x, y }, new[] { LUT_param.Call(xy, d_knots_cat, (DM)d_flat) });
    checkfunction(LUT, f, inputs);
    checkCodegen(f, inputs);
    checkSerialize(f, inputs);

    LUT_param = ca.interpolant("name", solver, new long[] { 5, 6 }, d_flat);
    f = ca.Function("LUTp", new[] { x, y }, new[] { LUT_param.Call(xy, d_knots_cat) });
    checkfunction(LUT, f, inputs);
    checkCodegen(f, inputs);
    checkSerialize(f, inputs);
  }

  public void test_2d_linear_multiout() {
    multiout("linear", true);
  }

  public void test_2d_bspline_multiout() {
    multiout("bspline", false);
  }

  /// <summary>B-spline sum_i c_i B_{i,p}(x), right-continuous at knots except at the right end (as scipy BSpline)</summary>
  static double bsplineRef(double[] t, double[] c, int p, double x) {
    int n = c.Length;
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

  public void test_bspline_repeated_knots() {
    var cases = new (int, double[])[] {
      (3, new double[] { 0, 0, 0, 0, 1, 1, 2, 3, 3, 3, 4, 4, 4, 4 }),
      (2, new double[] { 0, 0, 0, 1, 1, 2, 2, 3, 3, 3 }) };
    foreach (var (degree, knots) in cases) {
      int n = knots.Length - degree - 1;
      double[] c = rnd(n, 1).nonzeros().ToArray();
      MX x = MX.sym("x"), C = MX.sym("C", n);
      var f = ca.Function("f", new[] { x }, new[] { ca.bspline(x, (DM)c, new[] { knots }, new long[] { degree }, 1) });
      var F = ca.Function("F", new[] { x, C }, new[] { ca.bspline(x, C, new[] { knots }, new long[] { degree }, 1) });
      for (double xe = 0; xe < knots[^1] + 0.25; xe += 0.5) {
        double refv = bsplineRef(knots, c, degree, xe);
        checkarray(f.Call(xe), ca.DM(refv), 1e-12, "x=" + xe);
        checkarray(F.Call(xe, c), ca.DM(refv), 1e-12, "x=" + xe);
      }
      checkCodegen(f, ins(2));
      checkCodegen(F, ins(3, c));
    }
  }

  public void test_parametric_bspline() {
    double[][] knots = { new double[] { 0, 0, 0, 0, 0.2, 0.5, 0.8, 1, 1, 1, 1 }, new double[] { 0, 0, 0, 0.1, 0.5, 0.9, 1, 1, 1 } };
    MX x = MX.sym("x", 2);
    DM data = rnd(7 * 6 * 2, 1);
    MX y = ca.bspline(x, data, knots, new long[] { 3, 2 }, 2);
    MX C = MX.sym("C", data.size1(), 1);
    MX Y = ca.bspline(x, C, knots, new long[] { 3, 2 }, 2);
    var f = ca.Function("f", new[] { x }, new[] { y });
    var F = ca.Function("f", new[] { x, C }, new[] { Y });
    F = ca.Function("f", new[] { x }, new[] { F.Call(x, data) });
    var inputs = new List<DM> { ca.vertcat(0.3, 0.4) };
    checkfunction(f, F, inputs);
    checkCodegen(F, inputs);
    checkSerialize(F, inputs);

    var K = new List<MX> { MX.sym("K0", knots[0].Length), MX.sym("K1", knots[1].Length) };
    MX Y_pk = ca.bspline(x, C, K, new long[] { 3, 2 }, 2);
    var knotsDm = new List<MX> { knots[0], knots[1] };
    MX Y_sub = ca.substitute(new[] { Y_pk }, K, knotsDm)[0];
    var F_pk = ca.Function("f_pk", new[] { x, C }, new[] { Y_sub });
    F_pk = ca.Function("f_pk", new[] { x }, new[] { F_pk.Call(x, data) });
    checkfunction(f, F_pk, inputs);
    checkCodegen(F_pk, inputs, opts(), opts(), false, "c99", false);
  }

  /// <summary>np.random.seed(0) data of test_smooth_linear, which pins a reference value</summary>
  static readonly double[] SMOOTH_LINEAR_DATA = {0.5488135039273248, 0.26455561210462697, 0.4375872112626925,
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

  public void test_smooth_linear() {
    double[][] d_knots = { linspace(0, 1, 5), linspace(0, 1, 6) };
    double[] d_flat = SMOOTH_LINEAR_DATA;
    Function LUT_linear = ca.interpolant("name", "linear", d_knots, d_flat);
    Function LUT = ca.interpolant("name", "bspline", d_knots, d_flat, opts("algorithm", "smooth_linear", "smooth_linear_frac", 0.1));
    Function LUT2 = ca.interpolant("name", "bspline", d_knots, d_flat, opts("algorithm", "smooth_linear", "smooth_linear_frac", 0.05));

    double[] p1 = { 0.2, 0.333 }, p2 = { 0.25, 0.4 };
    double diff1 = (double)ca.norm_1(LUT.Call(p1) - LUT_linear.Call(p1));
    double diff2 = (double)ca.norm_2(LUT.Call(p2) - LUT_linear.Call(p2));
    double diff22 = (double)ca.norm_2(LUT2.Call(p2) - LUT_linear.Call(p2));
    assertTrue(diff1 <= 1e-10);
    assertTrue(diff2 >= 100 * diff1);
    assertTrue(diff22 <= diff2 * 0.51);
    assertTrue(diff22 >= diff2 * 0.49);
    checkarray(LUT2.Call(new double[] { 0.26, 0.39 }), ca.vertcat(0.9261362392504342, 0.9157189108791507));

    MX xy = MX.sym("xy", 2);
    Function LUT_param = ca.interpolant("name", "bspline", d_knots, 2, opts("algorithm", "smooth_linear", "smooth_linear_frac", 0.1));
    var f = ca.Function("LUTp", new[] { xy }, new[] { LUT_param.Call(xy, (DM)d_flat) });
    var inputs = new List<DM> { ca.vertcat(0.2, 0.333) };
    checkfunction(LUT, f, inputs);
    checkCodegen(f, inputs, opts(), opts(), false, "c89", true);
    checkSerialize(f, inputs);
  }

  public void test_codegen_avoid_stack() {
    SX x = SX.sym("x", 3, 3);
    var f = ca.Function("f", new[] { x }, new[] { ca.det(x) });
    checkCodegen(f, new List<DM> { rnd(3, 3) });
    checkCodegen(f, new List<DM> { rnd(3, 3) }, opts("avoid_stack", true), opts(), false, "c89", false);
  }

  public void test_serialize() {
    foreach (Dictionary<string, object> o in new[] { opts("debug", true), opts() }) {
      SX x = SX.sym("x");
      SX z = ca.sin(x + 3.0);
      var f = ca.Function("f", new[] { x }, new[] { z });
      Function fs = Function.deserialize(f.serialize(o));
      checkfunction(f, fs, ins(2));

      x = SX.sym("x");
      z = ca.sin(x + 3.0);
      f = ca.Function("f", new[] { x }, new[] { z, ca.SX(double.NaN), ca.SX(double.NegativeInfinity), ca.SX(double.PositiveInfinity) });
      fs = Function.deserialize(f.serialize(o));
      checkfunction(f, fs, ins(2));

      x = SX.sym("x");
      SX y = SX.sym("y", Sparsity.lower(3));
      z = x + y;
      SX z1 = ca.sparsify(ca.vertcat(z[0], ca.SX(0.0), z[1]));
      SX z2 = z.T();
      f = ca.Function("f", new[] { x, y }, new[] { z1, z2, ca.pow(x, 2.0) }, new[] { "x", "y" }, new[] { "a", "b", "c" });
      fs = Function.deserialize(f.serialize(o));
      assertEquals("x", fs.name_in(0));
      assertEquals("a", fs.name_out(0));
      assertEquals("f", fs.name());
      List<DM> inp = ins(3.7, new double[,] { { 1, 0, 0 }, { 2, 3, 0 }, { 4, 5, 6 } });
      checkfunction(f, fs, inp);
      fs = Function.deserialize(f.serialize(o));
      checkfunction(f, fs, inp);

      x = SX.sym("x");
      SX p = SX.sym("p");
      ca.Function("f", new[] { x }, new[] { p }, opts("allow_free", true));
      // pickle round trip of the SXFunction with free parameters: python-only (pickle)

      MX xm = MX.sym("x");
      f = ca.Function("f", new[] { xm }, new[] { ca.pow(xm, 2.0) });
      fs = Function.deserialize(f.serialize(o));
      checkfunction(f, fs, ins(3.7));

      xm = MX.sym("x");
      MX ym = MX.sym("y", 2);
      MX w = ca.if_else(xm, ca.atan2(ca.norm_fro(ym) * 3.0 * ym, xm), xm - ym, true);
      MX zm = ca.sin(xm * 2.0) * w[0] + 1.0;
      var g = ca.Function("g", new[] { xm, ym }, new[] { w - xm });
      Function gmap = g.map(2, "thread", 2);
      Function gmapsx = gmap.expand();
      MX q = gmap.Call(ca.horzcat(xm * 2.0, xm - ym[1]), ca.horzcat(zm + ym, ca.cos(zm + ym)))
          + 1.0 / gmapsx.Call(ca.horzcat(xm * 2.0, xm - ym[1]), ca.repmat(zm + ym, 1, 2));
      q = ca.solve(q, ym * 2.0, "qr");
      q = q + ca.bilin(ca.DM(new double[,] { { 1, 3 }, { 7, 8 } }), q, q * 2.0);
      f = ca.Function("f", new[] { xm, ym }, new[] { q + 1.0, ca.jacobian(q, ca.vertcat(xm, ym)) });
      fs = Function.deserialize(f.serialize(o));
      checkfunction(f, fs, new List<DM> { ca.DM(1.1), ca.vertcat(2.7, 3) });
    }
  }

  public void test_serialize_recursion_limit() {
    SX xs = SX.sym("x");
    SX ys = xs;
    for (int i = 0; i < 10000; i++) ys = ca.sin(ys);
    Function.deserialize(ca.Function("foo", new[] { xs }, new[] { ys }).serialize());
    MX xm = MX.sym("x");
    MX ym = xm;
    for (int i = 0; i < 10000; i++) ym = ca.sin(ym);
    Function.deserialize(ca.Function("foo", new[] { xm }, new[] { ym }).serialize());
  }

  public void test_string() {
    MX x = MX.sym("x"), y = MX.sym("y");
    var f = ca.Function("f", new[] { x }, new MX[0], new[] { "x" }, new string[0]);
    assertTrue(f.ToString().Contains("(x)->()"));
    f = ca.Function("f", new MX[0], new[] { x }, new string[0], new[] { "y" }, opts("allow_free", true));
    assertTrue(f.ToString().Contains("()->(y)"));
    f = ca.Function("f", new[] { x }, new[] { ca.pow(x, 2.0) }, new[] { "x" }, new[] { "y" });
    assertTrue(f.ToString().Contains("(x)->(y)"));
    f = ca.Function("f", new[] { x, y }, new[] { ca.pow(x, 2.0), x * y }, new[] { "x", "y" }, new[] { "w", "z" });
    assertTrue(f.ToString().Contains("(x,y)->(w,z)"));
  }

  public void test_fold() {
    SX x = SX.sym("x", 2, 2);
    var f = ca.Function("f", new[] { x }, new[] { ca.sin(x) });
    Function F = f.fold(10);
    SX x0 = x;
    for (int i = 0; i < 10; i++) x = f.Call(x);
    var Fref = ca.Function("f", new[] { x0 }, new[] { x });
    checkfunction(F, Fref, new List<DM> { ca.DM(new double[,] { { 1, 2 }, { 3, 7 } }) });
  }

  public void test_thread_safety() {
    MX x = MX.sym("x"), y = MX.sym("y");
    var f = ca.Function("f", new[] { x, y }, new[] { ca.pow(x, 2.0) - y });
    double[] ones = Enumerable.Repeat(1.0, 200).ToArray();
    double[] lin = linspace(0, 10, 200);
    foreach (string rf in new[] { "newton", "fast_newton", "bisection" }) {
      Function finv = ca.rootfinder("finv", rf, f);
      Function finvPar = finv.map(50, "unroll").map(4, "thread", 4);
      DM res = finvPar.Call(ones, lin);
      checkarray(ca.norm_inf(res.T() - ca.sqrt((DM)lin)), ca.DM(0), 1e-5, "");
      finvPar = finv.map(200, "thread", 4);
      res = finvPar.Call(ones, lin);
      checkarray(ca.norm_inf(res.T() - ca.sqrt((DM)lin)), ca.DM(0), 1e-5, "");
    }
  }

  /// <summary>Whether casadi was built WITH_THREAD (WITH_THREADSAFE_SYMBOLICS requires it): detected from map's fallback warning</summary>
  static bool builtWithThread() {
    MX x = MX.sym("x");
    var f = ca.Function("f", new[] { x }, new[] { x });
    string[] o = captureBoth(() => f.map(2, "thread").Call(ca.DM(new double[,] { { 1, 2 } })));
    return !(o[0] + o[1]).Contains("WITH_THREAD=ON");
  }

  public void test_threadsafe_symbolics_weakref() {
    // Not in Java: without WITH_THREADSAFE_SYMBOLICS, concurrent symbolic construction crashes casadi itself
    Check.Assume(builtWithThread(), "casadi built without WITH_THREAD/WITH_THREADSAFE_SYMBOLICS");
    int nthreads = 24, ntrials = 8000;
    Function box = null;
    var refs = new WeakRef[nthreads];
    bool stop = false;
    using var start = new Barrier(nthreads + 1);
    using var done = new Barrier(nthreads + 1);
    var threads = new List<Thread>();
    for (int t = 0; t < nthreads; t++) {
      int i = t;
      threads.Add(new Thread(() => {
        while (true) {
          start.SignalAndWait();
          if (Volatile.Read(ref stop)) return;
          refs[i] = ca.WeakRef(box);
          done.SignalAndWait();
        }
      }));
    }
    foreach (Thread t in threads) t.Start();
    int orphans = 0;
    try {
      for (int trial = 0; trial < ntrials; trial++) {
        MX x = MX.sym("x", 2);
        box = ca.Function("f", new[] { x }, new[] { ca.sin(x) });
        start.SignalAndWait();
        done.SignalAndWait();
        box.Dispose();
        box = null;
        foreach (WeakRef r in refs)
          if (r.alive()) orphans++;
        foreach (WeakRef r in refs) r.Dispose();
      }
    } finally {
      Volatile.Write(ref stop, true);
      start.SignalAndWait();
      foreach (Thread t in threads) t.Join();
    }
    assertEquals(0, orphans, "weak() double-create race: " + orphans + " orphaned weak references survived their owner's destruction");
  }

  public void test_mapped_eval() {
    SX x = SX.sym("x"), y = SX.sym("y", 2);
    var f = ca.Function("f", new[] { x, y }, new[] { ca.sin(x) * y }, new[] { "x", "y" }, new[] { "r" });
    DM r1 = f.Call(1, ca.vertcat(3, 4)), r2 = f.Call(2, ca.vertcat(3, 4)), r3 = f.Call(3, ca.vertcat(3, 4));
    DM r_all = f.Call(ca.DM(new double[,] { { 1, 2, 3 } }), ca.vertcat(3, 4));
    checkarray(r_all, ca.horzcat(r1, r2, r3), "Mapped evaluation (DM)");
    MX z = MX.sym("z", 1, 3);
    MX rz = f.Call(z, ca.vertcat(3, 4));
    var F = ca.Function("F", new[] { z }, new[] { rz });
    DM r_mx = F.Call(ca.DM(new double[,] { { 1, 2, 3 } }));
    checkarray(r_all, r_mx, "Mapped evaluation (MX)");
  }

  public void test_default_arg() {
    MX x = MX.sym("x"), y = MX.sym("y"), z = MX.sym("z");
    var f = ca.Function("f", new[] { x, y, z }, new[] { x, y, z }, new[] { "x", "y", "z" }, new[] { "a", "b", "c" },
        opts("default_in", new double[] { 1, 2, 3 }));
    checkCodegen(f, new Dictionary<string, DM> { { "x", ca.DM(5) }, { "z", ca.DM(3) } }, true);
  }

  public void test_factory_inherit_options() {
    MX x = MX.sym("x", 5);
    foreach (string op in new[] { "grad:f:x", "jac:f:x", "triu:hess:f:x:x" }) {
      var f = ca.Function("f", new[] { x }, new[] { ca.dot(x, x) }, new[] { "x" }, new[] { "f" }, opts("verbose", true));
      string o = captureStdout(() => f.factory("fgrad", new[] { "x" }, new[] { op }));
      assertTrue(o.Contains("::init"));

      var f2 = ca.Function("f", new[] { x }, new[] { ca.dot(x, x) }, new[] { "x" }, new[] { "f" },
          opts("enable_fd", true, "enable_forward", false, "enable_reverse", false, "print_in", true, "fd_method", "central",
              "fd_options", opts("h", 1e-7, "h_iter", false)));
      Function fgrad = f2.factory("fgrad", new[] { "x" }, new[] { op });
      o = captureStdout(() => fgrad.Call(0));
      assertTrue(o.Contains("-1e-07,") || o.Contains("-1e-007,"), o);
      assertTrue(o.Contains("1e-07,") || o.Contains("1e-007,"), o);
    }
  }

  static bool hasImporter(string name) => Importer.has_plugin(name);

  public void test_jit_codegen_options() {
    Check.Assume(hasImporter("shell"), "no importer shell");
    Check.Assume(RunSlow, "slow");
    MX x = MX.sym("x");
    string dir = Path.Combine(Path.GetTempPath(), "casadi_dotnet_jit_" + Guid.NewGuid().ToString("N"));
    Directory.CreateDirectory(dir);
    var o = opts("jit", true, "compiler", "shell", "jit_cleanup", false, "jit_temp_suffix", false,
        "jit_name", "jit_codegen_options", "jit_options", opts("directory", dir), "codegen_options", opts("prefix", "cgopt"));
    var f = ca.Function("f", new[] { x }, new[] { ca.pow(x, 2.0) }, o);
    foreach (Function g in new[] { f, Function.deserialize(f.serialize()) }) {
      checkarray(g.Call(3), 9);
      string src = File.ReadAllText(Path.Combine(dir, "jit_codegen_options.c")).Replace("\r\n", "\n");
      assertTrue(src.Contains("cgopt_"));
    }
  }

  static HashSet<string> compiled(string o) {
    var found = new HashSet<string>();
    foreach (Match mt in Regex.Matches(o, @"Compiling function '(\w+?)'")) found.Add(mt.Groups[1].Value);
    return found;
  }

  public void test_inherit_jit_options() {
    Check.Assume(hasNlpsol("ipopt"), "no nlpsol ipopt");
    Check.Assume(hasImporter("shell"), "no importer shell");
    Check.Assume(RunSlow, "slow");
    MX x = MX.sym("x");
    var f = ca.Function("Q", new[] { x }, new[] { ca.pow(x, 2.0) }, opts("jit", true, "compiler", "shell", "verbose", true));
    MX g = f.call(new[] { ca.sin(x) }, false, true)[0];
    MX gc = ca.cos(g);
    Function solver = null;
    string o = captureStdout(() => solver = ca.nlpsol("solver", "ipopt", new Dictionary<string, MX> { { "x", x }, { "f", ca.pow(x, 2.0) }, { "g", gc } }));
    HashSet<string> found = compiled(o);
    assertEquals(1, found.Count, string.Join(",", found));
    assertTrue(found.Contains("fwd1_Q"));

    o = captureStdout(() => checkSerialize(solver, new Dictionary<string, DM> { { "x0", ca.DM(2) } }));
    found = compiled(o);
    assertEquals(2, found.Count, string.Join(",", found));
    assertTrue(found.Contains("Q"));
    assertTrue(found.Contains("fwd1_Q"));
  }

  public void test_codegen_external_name_collision() {
    Check.Assume(hasImporter("shell"), "no importer shell");
    Check.Assume(RunSlow, "slow");
    MX x = MX.sym("x");
    var f = ca.Function("f", new[] { x }, new[] { x + 1.0 });
    f.generate("external_name_collision.c");
    var compiler = ca.Importer("external_name_collision.c", "shell");
    Function ext = ca.external("f", compiler);
    checkarray(ext.Call(2), 3);
    var wrapper = ca.Function("wrapper", new[] { x }, new[] { ext.Call(x) });
    assertInException("name conflicts with an external dependency", () => ext.generate("external_name_collision_wrapped.c"));
    foreach (Function[] functions in new[] { new[] { f, wrapper }, new[] { wrapper, f } }) {
      var cg1 = ca.CodeGenerator("external_name_collision_module.c");
      cg1.add(functions[0]);
      assertInException("name conflicts with an external dependency", () => cg1.add(functions[1]));
    }
    var cg = ca.CodeGenerator("external_name_collision_valid.c");
    cg.add(wrapper);
    cg.add(ca.Function("other", new[] { x }, new[] { x * 2.0 }));
    assertTrue(cg.dump().Contains("wrapper"));
  }

  public void test_generate_dependencies_jit() {
    Check.Assume(hasNlpsol("ipopt"), "no nlpsol ipopt");
    Check.Assume(hasImporter("shell"), "no importer shell");
    Check.Assume(RunSlow, "slow");
    MX x = MX.sym("x", 2);
    Function solver = ca.nlpsol("solver", "ipopt",
        new Dictionary<string, MX> { { "x", x }, { "f", ca.sumsqr(x) }, { "g", x[0] + x[1] - 10.0 } },
        opts("jit", true, "compiler", "shell", "print_time", false, "ipopt.print_level", 0));
    solver.generate_dependencies("generate_dependencies_jit.c");
    var compiler = ca.Importer("generate_dependencies_jit.c", "shell");
    foreach (string name in new[] { "nlp_f", "nlp_g" })
      checkfunctionLight(solver.get_function(name), ca.external(name, compiler), new List<DM> { ca.vertcat(1, 2), ca.DM(0) });
  }
}

}

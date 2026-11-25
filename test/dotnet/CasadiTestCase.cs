using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;

namespace CasADi.Tests {

/// <summary>Shared helpers, after test/java CasadiTestCase (itself after test/python/helpers.py)</summary>
public abstract class CasadiTestCase {
  protected readonly Random random = new Random(1);

  protected static bool RunSlow => Check.RunSlow;

  protected static bool hasNlpsol(string name) => ca.has_nlpsol(name);
  protected static bool hasConic(string name) => ca.has_conic(name);
  protected static bool hasIntegrator(string name) => ca.has_integrator(name);
  protected static bool hasRootfinder(string name) => ca.has_rootfinder(name);
  protected static bool hasLinsol(string name) => ca.has_linsol(name);

  /// <summary>Skips the test unless the plugin exists</summary>
  protected static void requireNlpsol(string name) => Check.Assume(ca.has_nlpsol(name), "no nlpsol " + name);
  protected static void requireConic(string name) => Check.Assume(ca.has_conic(name), "no conic " + name);
  protected static void requireIntegrator(string name) => Check.Assume(ca.has_integrator(name), "no integrator " + name);
  protected static void requireRootfinder(string name) => Check.Assume(ca.has_rootfinder(name), "no rootfinder " + name);
  protected static void requireLinsol(string name) => Check.Assume(ca.has_linsol(name), "no linsol " + name);

  protected static void message(string s) => Console.WriteLine(s);

  protected static List<DM> dms(params DM[] x) => x.ToList();

  double NextGaussian() {
    double u1 = 1.0 - random.NextDouble(), u2 = random.NextDouble();
    return Math.Sqrt(-2.0 * Math.Log(u1)) * Math.Cos(2.0 * Math.PI * u2);
  }

  /// <summary>Random dense matrix, normally distributed entries</summary>
  protected DM randDM(long n, long m) {
    double[,] a = new double[n, m];
    for (long i = 0; i < n; i++)
      for (long j = 0; j < m; j++) a[i, j] = NextGaussian();
    return n == 0 || m == 0 ? ca.DM(n, m) : (DM)a;
  }

  /// <summary>Random matrix with about a fraction density of nonzeros</summary>
  protected DM randDM(long n, long m, double density) {
    double[,] a = randDM(n, m).ToArray();
    var r = new List<long>();
    var c = new List<long>();
    var v = new List<double>();
    for (long j = 0; j < m; j++)
      for (long i = 0; i < n; i++)
        if (random.NextDouble() < density) {
          r.Add(i);
          c.Add(j);
          v.Add(a[i, j]);
        }
    Sparsity sp = Sparsity.triplet(n, m, r, c);
    return ca.DM(sp, (DM)v.ToArray());
  }

  protected static void checkarray(DM a, DM b) => checkarray(a, b, 1e-10, "");
  protected static void checkarray(DM a, DM b, string msg) => checkarray(a, b, 1e-10, msg);
  protected static void checkarray(DM a, double[,] b) => checkarray(a, (DM)b, 1e-10, "");
  protected static void checkarray(DM a, double b) => checkarray(a, ca.DM(b), 1e-10, "");

  /// <summary>Same shape, entries equal to within tol (absolute, scaled by max(1,|b|))</summary>
  protected static void checkarray(DM a, DM b, double tol, string msg) {
    string where = msg.Length == 0 ? "" : msg + ": ";
    Check.Equal(b.size1(), a.size1(), where + "size1 of " + a + " vs " + b);
    Check.Equal(b.size2(), a.size2(), where + "size2 of " + a + " vs " + b);
    double[,] x = a.ToArray(), y = b.ToArray();
    for (int i = 0; i < x.GetLength(0); i++)
      for (int j = 0; j < x.GetLength(1); j++) {
        double u = x[i, j], w = y[i, j];
        if (double.IsNaN(u) && double.IsNaN(w)) continue;
        if (u == w) continue;
        if (!(Math.Abs(u - w) <= tol * Math.Max(1, Math.Abs(w))))
          throw new AssertFailed(where + "entry (" + i + "," + j + "): " + u + " vs " + w + "\n" + a + "\n" + b);
      }
  }

  protected static void checkarray(double a, double b, double tol) =>
    Check.True(Math.Abs(a - b) <= tol * Math.Max(1, Math.Abs(b)), a + " vs " + b);

  protected static void assertEquals<T>(T expected, T actual, string msg = "") => Check.Equal(expected, actual, msg);
  protected static void assertTrue(bool c, string msg = "expected true") => Check.True(c, msg);
  protected static void assertFalse(bool c, string msg = "expected false") => Check.True(!c, msg);

  /// <summary>Asserts that r throws and its message contains s</summary>
  protected static void assertInException(string s, Action r) {
    try {
      r();
    } catch (Exception e) when (!(e is AssertFailed)) {
      Check.True((e.Message ?? "").Contains(s), "expected '" + s + "' in: " + e.Message);
      return;
    }
    throw new AssertFailed("expected an exception containing '" + s + "'");
  }

  /// <summary>What r prints to Console.Out, where casadi printout goes</summary>
  protected static string captureStdout(Action r) => captureBoth(r)[0];

  /// <summary>What r prints to Console.Out and Console.Error</summary>
  protected static string[] captureBoth(Action r) {
    TextWriter out0 = Console.Out, err0 = Console.Error;
    var bo = new StringWriter();
    var be = new StringWriter();
    Console.SetOut(bo);
    Console.SetError(be);
    try {
      r();
    } finally {
      Console.SetOut(out0);
      Console.SetError(err0);
    }
    return new[] { bo.ToString(), be.ToString() };
  }

  /// <summary>helpers.py assertInAnyOutput: s in the exception message, stdout or stderr</summary>
  protected static void assertInAnyOutput(string s, Action r) {
    string e = "";
    string[] res = captureBoth(() => {
      try {
        r();
      } catch (Exception err) when (!(err is AssertFailed)) {
        e = err.Message ?? "";
      }
    });
    Console.Write(res[0]);
    Check.True(e.Contains(s) || res[0].Contains(s) || res[1].Contains(s),
               "error: <<<" + e + ">>> stdout: <<<" + res[0] + ">>> stderr: <<<" + res[1] + ">>>");
  }

  protected static void assertOutput(string included, Action r) {
    string o = captureStdout(r);
    Check.True(o.Contains(included), "expected '" + included + "' in output:\n" + o);
  }

  /// <summary>Random inputs matching the sparsity of f's inputs</summary>
  protected List<DM> randomInputs(Function f) {
    var r = new List<DM>();
    for (long i = 0; i < f.n_in(); i++) r.Add(randFor(f.sparsity_in(i)));
    return r;
  }

  /// <summary>f and g agree on inputs, also their first-order derivatives (as checkfunction in helpers.py)</summary>
  protected void checkfunction(Function f, Function g, List<DM> inputs) => checkfunction(f, g, inputs, 1e-9);

  protected void checkfunction(Function f, Function g, List<DM> inputs, double tol) {
    Check.Equal(g.n_in(), f.n_in(), "n_in");
    Check.Equal(g.n_out(), f.n_out(), "n_out");
    List<DM> a = f.call(inputs), b = g.call(inputs);
    for (int i = 0; i < a.Count; i++) checkarray(a[i], b[i], tol, "output " + i);
    var ja = inputs.Concat(a).ToList();
    var jb = inputs.Concat(b).ToList();
    List<DM> va = f.jacobian().call(ja), vb = g.jacobian().call(jb);
    for (int i = 0; i < va.Count; i++) checkarray(va[i], vb[i], tol, "jacobian block " + i);
    // forward and reverse sweeps with the same seeds
    var fwd = inputs.Concat(a).ToList();
    var adj = inputs.Concat(a).ToList();
    for (long i = 0; i < f.n_in(); i++) fwd.Add(randFor(f.sparsity_in(i)));
    for (long i = 0; i < f.n_out(); i++) adj.Add(randFor(f.sparsity_out(i)));
    List<DM> fa = f.forward(1).call(fwd), fb = g.forward(1).call(fwd);
    for (int i = 0; i < fa.Count; i++) checkarray(fa[i], fb[i], tol, "forward " + i);
    List<DM> ra = f.reverse(1).call(adj), rb = g.reverse(1).call(adj);
    for (int i = 0; i < ra.Count; i++) checkarray(ra[i], rb[i], tol, "reverse " + i);
  }

  DM randFor(Sparsity sp) {
    double[] v = new double[sp.nnz()];
    for (int k = 0; k < v.Length; k++) v[k] = NextGaussian();
    return ca.DM(sp, v.Length == 0 ? ca.DM(0, 1) : (DM)v);
  }

  protected static SX sxsym(string name, long n, long m) => SX.sym(name, n, m);
  protected static MX mxsym(string name, long n, long m) => MX.sym(name, n, m);
}

}

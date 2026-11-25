using System;
using System.Collections;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;
using System.Runtime.InteropServices;

namespace CasADi.Tests {

/// <summary>Port of test/java NlpTest (itself a port of test/python/nlp.py)</summary>
public class NlpTest : CasadiTestCase {

  const double inf = double.PositiveInfinity;
  const double pi = Math.PI;
  const int codegen_check_digits = 15;

  /// <summary>One row of the module-level solvers table: (Solver, solver_options, aux_options)</summary>
  sealed class S {
    public readonly string name;
    public readonly Dictionary<string, object> opts;
    public readonly Dictionary<string, object> codegen;
    public readonly bool discrete;

    public S(string name, Dictionary<string, object> opts, Dictionary<string, object> codegen, bool discrete) {
      this.name = name;
      this.opts = opts;
      this.codegen = codegen;
      this.discrete = discrete;
    }

    /// <summary>Python: s in str(solver_options)</summary>
    public bool @in(string s) => Str(opts).Contains(s);

    public override string ToString() => name + " " + Str(opts);
  }

  static string Str(object o) {
    if (o is string s) return s;
    if (o is IDictionary d) {
      var parts = new List<string>();
      foreach (DictionaryEntry e in d) parts.Add(e.Key + "=" + Str(e.Value));
      return "{" + string.Join(", ", parts) + "}";
    }
    if (o is IEnumerable l) return "[" + string.Join(", ", l.Cast<object>().Select(Str)) + "]";
    return o?.ToString() ?? "null";
  }

  static bool env(string s) => Environment.GetEnvironmentVariable(s) != null;

  static readonly List<S> solvers = MakeSolvers();

  static List<S> MakeSolvers() {
    var solvers = new List<S>();
    if (!env("SKIP_WORHP_TESTS") && hasNlpsol("worhp"))
      solvers.Add(new S("worhp", d("worhp", d("TolOpti", 1e-9)), null, false));
    if (!env("SKIP_FATROP_TESTS") && hasNlpsol("fatrop"))
      solvers.Add(new S("fatrop", d("fatrop", d()), d("std", "c99", "extralibs", new[] { "fatrop", "blasfeo" },
          "extra_options", new[] { "-Wno-strict-prototypes" }), false));
    if (!env("SKIP_SLEQP_TESTS") && hasNlpsol("sleqp"))
      solvers.Add(new S("sleqp", d("print_time", false, "sleqp", d("linesearch", "Approx", "feas_tol", 1e-7, "stat_tol",
          1e-7, "slack_tol", 1e-7, "hess_eval", "Exact")), null, false));
    if (!env("SKIP_UNO_TESTS") && hasNlpsol("uno")) {
      var uno_codegen = d("std", "c99", "extralibs", new[] { "uno" }, "extra_include", new[] { "uno" });
      foreach (var u in new[] {
          d("preset", "ipopt", "primal_tolerance", 1e-8, "dual_tolerance", 1e-8),
          d("preset", "funnelsqp", "primal_tolerance", 1e-8, "dual_tolerance", 1e-8),
          d("preset", "filtersqp", "primal_tolerance", 1e-8, "dual_tolerance", 1e-8),
          d("preset", "filtersqp", "primal_tolerance", 1e-8, "dual_tolerance", 1e-8, "hessian_model", "LBFGS") })
        solvers.Add(new S("uno", d("print_time", false, "uno", u), uno_codegen, false));
    }
    if (!env("SKIP_ALPAQA_TESTS") && hasNlpsol("alpaqa"))
      solvers.Add(new S("alpaqa", d("print_time", false, "alpaqa", d("alm.tolerance", 1e-10, "alm.dual_tolerance", 1e-10,
          "alm.penalty_update_factor", 10, "alm.max_iter", 3000, "alm.print_interval", 1, "panoc.max_iter", 500,
          "panoc.print_interval", 1, "lbfgs.memory", 2)), null, false));
    if (!env("SKIP_IPOPT_TESTS") && hasNlpsol("ipopt")) {
      var codegen = d("extralibs", new[] { "ipopt" }, "std", "c99");
      solvers.Add(new S("ipopt", d("print_time", false, "ipopt", d("tol", 1e-10, "derivative_test", "second-order",
          "print_level", 0)), codegen, false));
      solvers.Add(new S("ipopt", d("print_time", false, "ipopt", d("tol", 1e-10, "derivative_test", "first-order",
          "hessian_approximation", "limited-memory", "print_level", 0)), codegen, false));
    }
    if (!env("SKIP_IPOPT_TESTS") && hasNlpsol("ipopt") && hasNlpsol("sqpmethod")) {
      var qpsol_options = d("nlpsol", "ipopt", "nlpsol_options", d("ipopt.tol", 1e-12,
          "ipopt.tiny_step_tol", 1e-20, "ipopt.fixed_variable_treatment", "make_constraint", "ipopt.print_level", 0,
          "print_time", false));
      solvers.Add(new S("sqpmethod", d("qpsol", "nlpsol", "qpsol_options", qpsol_options, "print_header", false,
          "print_iteration", true, "print_time", false), null, false));
      solvers.Add(new S("sqpmethod", d("qpsol", "nlpsol", "qpsol_options", qpsol_options, "hessian_approximation",
          "limited-memory", "tol_du", 1e-10, "tol_pr", 1e-10, "min_step_size", 1e-14, "print_header", false,
          "print_iteration", true, "print_time", false, "max_iter", 55), null, false));
    }
    if (!env("SKIP_QRQP_TESTS") && hasConic("qrqp") && hasNlpsol("sqpmethod")) {
      var codegen = d("std", "c99");
      var qpsol_options = d("print_iter", false, "print_header", false, "error_on_fail", false);
      solvers.Add(new S("sqpmethod", d("qpsol", "qrqp", "qpsol_options", qpsol_options, "print_header", false,
          "print_iteration", false, "print_time", false), codegen, false));
      solvers.Add(new S("sqpmethod", d("qpsol", "qrqp", "max_iter_ls", 0, "qpsol_options", qpsol_options,
          "print_header", false, "print_iteration", false, "print_time", false), codegen, false));
      solvers.Add(new S("sqpmethod", d("qpsol", "qrqp", "convexify_strategy", "regularize", "max_iter", 500,
          "qpsol_options", qpsol_options, "print_header", false, "print_iteration", true, "print_time", false,
          "tol_du", 1e-8, "min_step_size", 1e-12), codegen, false));
    }
    if (!env("SKIP_DAQP_TESTS") && hasConic("daqp") && hasNlpsol("sqpmethod"))
      solvers.Add(new S("sqpmethod", d("qpsol", "daqp", "hessian_approximation", "limited-memory", "tol_du", 1e-10,
          "tol_pr", 1e-10, "min_step_size", 1e-14, "max_iter", 55, "print_header", false, "print_iteration", true,
          "print_time", false, "qpsol_options", d("print_problem", true, "print_out", true)),
          d("std", "c99", "extralibs", new[] { "daqp" }), false));
    if (!env("SKIP_BLOCKSQP_TESTS") && hasNlpsol("blocksqp")) {
      try {
        ca.load_linsol("ma27");
        solvers.Add(new S("blocksqp", d(), null, false));
      } catch (SwigException) {
        // ma27 unavailable
      }
    }
    if (!env("SKIP_BONMIN_TESTS") && hasNlpsol("bonmin")) solvers.Add(new S("bonmin", d(), null, true));
    if (!env("SKIP_KNITRO_TESTS") && hasNlpsol("knitro"))
      solvers.Add(new S("knitro", d("knitro", d("feastol", 1e-9, "opttol", 1e-9, "ftol", 1e-20)), null, false));
    if (!env("SKIP_SNOPT_TESTS") && hasNlpsol("snopt"))
      solvers.Add(new S("snopt", d("snopt", d("Verify_level", 3, "Major_optimality_tolerance", 1e-12,
          "Minor_feasibility_tolerance", 1e-12, "Major_feasibility_tolerance", 1e-12)), null, false));
    if (!env("SKIP_MADNLP_TESTS") && hasNlpsol("madnlp"))
      solvers.Add(new S("madnlp", d("madnlp", d()), libmadCodegen(), false));
    if (!env("SKIP_CONOPT_TESTS") && hasNlpsol("conopt"))
      solvers.Add(new S("conopt", d("conopt", d("Tol_Feas_Min", 1e-9, "Tol_Optimality", 1e-9)), null, false));
    return solvers;
  }

  static Dictionary<string, object> libmadCodegen() {
    string libmad_dir = Environment.GetEnvironmentVariable("LIBMADDIR") ?? "/missing";
    return d("extralibs", new[] { "Mad" }, "std", "c99", "extralibdirs", new[] { libmad_dir + "/lib" },
        "extra_include", new[] { libmad_dir + "/include" });
  }

  readonly List<object> keep = new List<object>();

  // ---------------------------------------------------------------- helpers

  /// <summary>Dict literal</summary>
  static Dictionary<string, object> d(params object[] kv) {
    var r = new Dictionary<string, object>();
    for (int i = 0; i < kv.Length; i += 2) r[(string)kv[i]] = kv[i + 1];
    return r;
  }

  /// <summary>dict(base, **kv)</summary>
  static Dictionary<string, object> with(IDictionary<string, object> b, params object[] kv) {
    var r = new Dictionary<string, object>(b);
    foreach (var e in d(kv)) r[e.Key] = e.Value;
    return r;
  }

  static DM toDM(object v) {
    switch (v) {
      case DM x: return x;
      case double x: return x;
      case int x: return (double)x;
      case long x: return (double)x;
      case double[] x: return x;
    }
    throw new ArgumentException(v?.ToString() ?? "null");
  }

  static MX toMX(object v) => v is MX m ? m : (MX)toDM(v);
  static SX toSX(object v) => v is SX s ? s : (SX)toDM(v);

  static Dictionary<string, DM> @in(params object[] kv) {
    var r = new Dictionary<string, DM>();
    for (int i = 0; i < kv.Length; i += 2) r[(string)kv[i]] = toDM(kv[i + 1]);
    return r;
  }

  static Dictionary<string, MX> min(params object[] kv) {
    var r = new Dictionary<string, MX>();
    for (int i = 0; i < kv.Length; i += 2) r[(string)kv[i]] = toMX(kv[i + 1]);
    return r;
  }

  static Dictionary<string, SX> sxnlp(params object[] kv) {
    var r = new Dictionary<string, SX>();
    for (int i = 0; i < kv.Length; i += 2) r[(string)kv[i]] = toSX(kv[i + 1]);
    return r;
  }

  static Dictionary<string, MX> mxnlp(params object[] kv) => min(kv);

  static double[] rep(double v, int n) => Enumerable.Repeat(v, n).ToArray();

  static double[] range(int n) => Enumerable.Range(0, n).Select(i => (double)i).ToArray();

  /// <summary>Python x[k] on a DM, as a float</summary>
  static double el(DM x, int k) => ca.densify(x).nonzeros()[k];

  static object stat(Function f, string k) => f.stats()[k];

  static bool success(Function f) => (bool)stat(f, "success");

  static long iterCount(Function f) => Convert.ToInt64(stat(f, "iter_count"));

  /// <summary>helpers.py assertAlmostEqual: places after scaling large values</summary>
  static void almost(double first, double second, int places, string msg = "") {
    if (double.IsNaN(first) && double.IsNaN(second)) return;
    if (first == second) return;
    double n = Math.Max(Math.Abs(first), Math.Abs(second));
    n = n > 1e3 ? Math.Pow(10, Math.Floor(Math.Log10(n))) : 1.0;
    double diff = Math.Abs(first / n - second / n);
    Check.True(Math.Round(diff * Math.Pow(10, places), MidpointRounding.AwayFromZero) == 0,
        msg + " " + first + " <-> " + second + " (places=" + places + ")");
  }

  /// <summary>helpers.py checkarray with digits</summary>
  static void chk(DM a, DM b, int digits = 10, string msg = "") {
    Check.Equal(b.size1(), a.size1(), msg + ": size1 " + a + " vs " + b);
    Check.Equal(b.size2(), a.size2(), msg + ": size2 " + a + " vs " + b);
    double[,] x = a.ToArray(), y = b.ToArray();
    for (int i = 0; i < x.GetLength(0); i++)
      for (int j = 0; j < x.GetLength(1); j++)
        almost(x[i, j], y[i, j], digits, msg + " evaluation error at (" + i + "," + j + ")\n" + a + "\n" + b);
  }

  /// <summary>Python list literal as a row, as checkarray(x, [..]) does</summary>
  static DM row(params double[] v) {
    var r = new double[1, v.Length];
    for (int i = 0; i < v.Length; i++) r[0, i] = v[i];
    return r;
  }

  static void assertOutput(IEnumerable<string> included, IEnumerable<string> excluded, Action r) {
    string o = captureStdout(r);
    foreach (string s in included) Check.True(o.Contains(s), "expected '" + s + "' in output:\n" + o);
    foreach (string s in excluded) Check.True(!o.Contains(s), "unexpected '" + s + "' in output:\n" + o);
  }

  /// <summary>Positional inputs from a dict, missing ones 0, as helpers.py checkfunction</summary>
  static List<DM> positional(Function f, IDictionary<string, DM> dd) =>
    f.name_in().Select(n => dd.TryGetValue(n, out DM v) ? v : ca.DM(0.0)).ToList();

  /// <summary>helpers.py checkfunction_light: repeated evaluation agrees</summary>
  static void checkfunctionLight(Function trial, Function solution, List<DM> inputs, int digits, string failmessage = "") {
    for (int k = 0; k < 2; k++) {
      List<DM> a = trial.call(inputs), b = solution.call(inputs);
      Check.Equal(solution.n_out(), trial.n_out(), failmessage);
      Check.Equal(solution.n_in(), trial.n_in(), failmessage);
      for (int i = 0; i < trial.n_out(); i++) {
        if (trial.sparsity_out(i).is_empty() || solution.sparsity_out(i).is_empty()) continue;
        chk(a[i], b[i], digits, failmessage + ": output(" + i + ": '" + trial.name_out(i) + "')");
      }
    }
  }

  /// <summary>helpers.py check_serialize (without the python serialization_scheme validation)</summary>
  static void checkSerialize(Function f, Dictionary<string, DM> inputs) {
    string s = f.serialize(d("debug", true));
    Function f2 = Function.deserialize(s);
    Dictionary<string, DM> a = f.call(inputs), b = f2.call(inputs);
    Check.True(f.name_out().SequenceEqual(f2.name_out()), "name_out");
    foreach (string k in f.name_out()) chk(a[k], b[k], 16, k);
  }

  static int sh(string cmd, string stdin = null, string stdout = null) {
    string c = cmd + (stdin != null ? " < " + stdin : "") + (stdout != null ? " > " + stdout : " 2>&1");
    var psi = new ProcessStartInfo("sh") { UseShellExecute = false, RedirectStandardOutput = true };
    psi.ArgumentList.Add("-c");
    psi.ArgumentList.Add(c);
    using var p = Process.Start(psi);
    Console.Write(p.StandardOutput.ReadToEnd());
    p.WaitForExit();
    return p.ExitCode;
  }

  static List<string> strs(IDictionary<string, object> cg, string key) =>
    cg.TryGetValue(key, out object v) ? ((IEnumerable)v).Cast<string>().ToList() : new List<string>();

  /// <summary>helpers.py check_codegen (gcc, linux): generate, compile, load as external and compare; run_slow only</summary>
  Function checkCodegen(Function f, Dictionary<string, DM> inputs, Dictionary<string, object> cg, int digits = 15,
      bool debug = false, bool main = false, long[] main_return_code = null) {
    if (!RunSlow) return null;
    string name = "codegen_" + Guid.NewGuid().ToString("N");
    var gen = ca.CodeGenerator(name, main ? d("main", true) : d());
    gen.add(f);
    gen.generate();
    string libdir = GlobalOptions.getCasadiPath(), includedir = GlobalOptions.getCasadiIncludePath();
    var inc = new List<string> { includedir, includedir + "/highs" };
    foreach (string e in strs(cg, "extra_include")) inc.Add(includedir + "/" + e);
    var libdirs = new List<string> { libdir };
    libdirs.AddRange(strs(cg, "extralibdirs"));
    // libraries that python_install bundles next to libcasadi live on LD_LIBRARY_PATH here
    string ld = Environment.GetEnvironmentVariable("LD_LIBRARY_PATH");
    if (ld != null)
      foreach (string p in ld.Split(':'))
        if (p.Length > 0) {
          libdirs.Add(p);
          inc.Add(Path.GetDirectoryName(p) + "/include");
        }
    string flags = string.Concat(libdirs.Select(l => " -L" + l + " -Wl,-rpath," + l));
    string extra = string.Concat(strs(cg, "extralibs").Select(l => l.StartsWith("-") || l.Contains(".") ? " " + l : " -l" + l))
      + string.Concat(strs(cg, "extra_options").Select(o => " " + o));
    string incs = string.Concat(inc.Select(i => " -I" + i));
    string std = cg.TryGetValue("std", out object sv) ? (string)sv : "c89";
    string opt = debug ? "-O0 -g" : "-O3";
    string common = " -pedantic -std=" + std + " -fPIC -Wall -Werror -Wextra" + incs
        + " -Wno-dangling-pointer -Wno-unknown-warning-option -Wno-unknown-pragmas -Wno-long-long -Wno-unused-parameter " + opt;
    string lib = "./" + name + ".so";
    try {
      sh("gcc" + common + " -shared " + name + ".c -o " + name + ".so" + flags + " -Wl,-rpath,." + extra);
      Function f2 = ca.external(f.name(), lib);
      if (main) {
        sh("gcc" + common + " " + name + ".c -o " + name + flags + " -Wl,-rpath,. -lm" + extra);
        string fin = f.name() + "_in.txt", fout = f.name() + "_out.txt";
        try {
          f.generate_in(fin, f.convert_in(inputs));
          long rc = sh("./" + name + " " + f.name(), fin, fout);
          Check.True(main_return_code.Contains(rc), "return code " + rc);
          if (rc != 0) return f2;
          List<DM> refv = f.call(f.convert_in(inputs));
          List<DM> outputs = f.generate_out(fout);
          for (int i = 0; i < f.n_out(); i++) chk(refv[i], outputs[i], digits);
        } finally {
          File.Delete(name);
          File.Delete(fin);
          File.Delete(fout);
        }
      }
      Dictionary<string, DM> a = f.call(inputs), b = f2.call(inputs);
      Check.True(f.name_out().SequenceEqual(f2.name_out()), "name_out");
      foreach (string k in f.name_out()) chk(a[k], b[k], digits, k);
      return f2;
    } finally {
      File.Delete(name + ".c");
    }
  }

  /// <summary>Symmetric eigendecomposition (Jacobi): (eigenvalues, eigenvectors as columns)</summary>
  static (double[], double[,]) eigh(double[,] A0) {
    int n = A0.GetLength(0);
    var A = (double[,])A0.Clone();
    var V = new double[n, n];
    for (int i = 0; i < n; i++) V[i, i] = 1;
    for (int sweep = 0; sweep < 100; sweep++) {
      double off = 0;
      for (int i = 0; i < n; i++)
        for (int j = i + 1; j < n; j++) off += A[i, j] * A[i, j];
      if (off < 1e-30) break;
      for (int p = 0; p < n; p++)
        for (int q = p + 1; q < n; q++) {
          if (A[p, q] == 0) continue;
          double theta = (A[q, q] - A[p, p]) / (2 * A[p, q]);
          double t = Math.Sign(theta) / (Math.Abs(theta) + Math.Sqrt(theta * theta + 1));
          if (theta == 0) t = 1;
          double c = 1 / Math.Sqrt(t * t + 1), s = t * c;
          for (int k = 0; k < n; k++) {
            double akp = A[k, p], akq = A[k, q];
            A[k, p] = c * akp - s * akq;
            A[k, q] = s * akp + c * akq;
          }
          for (int k = 0; k < n; k++) {
            double apk = A[p, k], aqk = A[q, k];
            A[p, k] = c * apk - s * aqk;
            A[q, k] = s * apk + c * aqk;
          }
          for (int k = 0; k < n; k++) {
            double vkp = V[k, p], vkq = V[k, q];
            V[k, p] = c * vkp - s * vkq;
            V[k, q] = s * vkp + c * vkq;
          }
        }
    }
    var dd = new double[n];
    for (int i = 0; i < n; i++) dd[i] = A[i, i];
    return (dd, V);
  }

  /// <summary>V diag(d) V^T</summary>
  static DM recompose(double[] dd, double[,] V) {
    int n = dd.Length;
    var r = new double[n, n];
    // numpy's np.dot(V, np.dot(D, V.T)) evaluation order; the clip cases are ill-conditioned
    for (int i = 0; i < n; i++)
      for (int j = 0; j < n; j++)
        for (int k = 0; k < n; k++) r[i, j] += V[i, k] * (dd[k] * V[j, k]);
    return r;
  }

  static DM reflect(DM A, double eps) {
    var (e, V) = eigh(ca.densify(A).ToArray());
    return recompose(e.Select(Math.Abs).ToArray(), V);
  }

  static DM clip(DM A, double eps) {
    var (e, V) = eigh(ca.densify(A).ToArray());
    return recompose(e.Select(v => v < eps ? eps : v).ToArray(), V);
  }

  static SX sq(SX x) => ca.pow(x, 2);
  static MX sq(MX x) => ca.pow(x, 2);

  /// <summary>(a-x)**2+100*(y-x**2)**2</summary>
  static SX rosen(SX a, SX x, SX y) => sq(a - x) + 100 * sq(y - sq(x));

  static DM vc(params double[] x) => x;

  // ---------------------------------------------------------------- tests

  public void test_ccopt() {
    Check.Assume(hasNlpsol("ccopt"), "no nlpsol ccopt");
    SX x = SX.sym("x"), y = SX.sym("y");
    SX f = sq(x - 2) + sq(y - 2);
    var nlp = sxnlp("x", ca.vertcat(x, y), "f", f);
    var solver_opts = d("cc_pairs", new[] { new[] { 0L, 1L } }, "madnlp", d("bound_relax_factor", 0.0));
    Function solver = ca.nlpsol("solver", "ccopt", nlp, solver_opts);
    var solver_in = @in("x0", vc(2.5, 3.0), "lbx", vc(0.0, 0.0));
    var solver_out = solver.call(solver_in);
    Check.True(success(solver));
    almost((double)solver_out["f"], 4.0, 6);
    double x_sol = el(solver_out["x"], 0), y_sol = el(solver_out["x"], 1);
    almost(Math.Min(x_sol, y_sol), 0.0, 6);
    almost(Math.Max(x_sol, y_sol), 2.0, 6);
    checkCodegen(solver, solver_in, libmadCodegen());
    checkSerialize(solver, solver_in);
  }

  public void test_alpaqa() {
    Check.Assume(hasNlpsol("alpaqa"), "no nlpsol alpaqa");
    SX x = SX.sym("x"), y = SX.sym("y");
    var nlp = sxnlp("x", ca.vertcat(x, y), "f", rosen(1, x, y), "g", x + y);
    Function solver = ca.nlpsol("mysolver", "alpaqa", nlp);
    solver.call(@in("lbx", rep(-10, 2), "ubx", rep(10, 2), "lbg", vc(-10), "ubg", vc(10)));
  }

  public void test_affine_constant_terms() {
    SX x = SX.sym("x");
    foreach (S s in solvers) {
      var solver_options = s.opts;
      if (s.name == "bonmin") solver_options = with(solver_options, "calc_g", true);
      message("affine objective with a constant term: " + s.name);
      Function solver = ca.nlpsol("mysolver", s.name, sxnlp("x", x, "f", 2 * x + 7), solver_options);
      var solver_out = solver.call(@in("x0", 0, "lbx", -10, "ubx", 10));
      chk(solver_out["x"], vc(-10), 6);
      chk(solver_out["f"], vc(-13), 6);

      message("fully linear constraint with a constant term: " + s.name);
      solver = ca.nlpsol("mysolver", s.name, sxnlp("x", x, "f", sq(x - 3), "g", x + 5), solver_options);
      solver_out = solver.call(@in("x0", 0, "lbg", 0, "ubg", 0));
      chk(solver_out["x"], vc(-5), 6);
      chk(solver_out["g"], vc(0), 6);

      message("fully linear range constraint with a constant term: " + s.name);
      solver = ca.nlpsol("mysolver", s.name, sxnlp("x", x, "f", sq(x - 3), "g", x + 5), solver_options);
      solver_out = solver.call(@in("x0", 0, "lbg", 0, "ubg", 1));
      chk(solver_out["x"], vc(-4), 6);
      chk(solver_out["g"], vc(1), 6);
    }
  }

  public void test_conopt_constant_constraints() {
    Check.Assume(hasNlpsol("conopt"), "no nlpsol conopt");
    SX x = SX.sym("x"), p = SX.sym("p");
    foreach (SX other in new[] { SX.zeros(0, 1), x * x, x }) {
      SX g = ca.vertcat(other, p);
      Function solver0 = ca.nlpsol("solver", "conopt", sxnlp("x", x, "p", p, "f", sq(x - 2), "g", g),
          d("error_on_fail", false, "print_time", false));
      foreach (Function solver in new[] { solver0, Function.deserialize(solver0.serialize()) }) {
        foreach (double[] c in new[] { new double[] { 5, -10, 10 }, new double[] { 5, 0, 0 }, new double[] { 0, 0, 0 },
                                       new double[] { 5, 4, 6 }, new double[] { 5, -inf, inf } }) {
          double value = c[0], lower = c[1], upper = c[2];
          int n = (int)other.numel();
          double[] lbg = rep(-10, n + 1), ubg = rep(10, n + 1);
          lbg[n] = lower;
          ubg[n] = upper;
          var result = solver.call(@in("x0", 1, "p", value, "lbg", lbg, "ubg", ubg));
          string sub = other + " " + value + " " + lower + " " + upper;
          Check.Equal(lower <= value && value <= upper, success(solver), sub);
          chk(result["g"][^1], ca.DM(value), 10, sub);
          if (success(solver)) chk(result["x"], ca.DM(2), 6, sub);
        }
      }
    }
  }

  public void test_conopt_evaluation_errors() {
    Check.Assume(hasNlpsol("conopt"), "no nlpsol conopt");
    SX x = SX.sym("x"), p = SX.sym("p");
    var cases = new (SX, SX, double)[] {
      (ca.sqrt(x + p), x, 0.0),
      (sq(x - 1), ca.sqrt(x + p), 0.0),
      (p * x, x, inf),
      (sq(x - 1), p * x, inf),
      (sq(x - 1), p, double.NaN),
    };
    foreach (var (f, g, bad_p) in cases) {
      Function solver = ca.nlpsol("solver", "conopt", sxnlp("x", x, "p", p, "f", f, "g", g),
          d("error_on_fail", false, "print_time", false));
      foreach (double value in new[] { 1, bad_p, 1 }) {
        solver.call(@in("x0", 0, "p", value, "lbx", 0, "ubx", 2, "lbg", -10, "ubg", 10));
        string sub = f + " " + g + " " + value;
        Check.Equal(value == 1, success(solver), sub);
        Check.Equal(value == 1 ? "SOLVER_RET_SUCCESS" : "SOLVER_RET_NAN", (string)stat(solver, "unified_return_status"), sub);
      }
    }
  }

  public void test_nonregular_point() {
    SX x = SX.sym("x");
    var nlp = sxnlp("x", x, "f", sq(x + 1), "g", ca.sqrt(x));
    foreach (S s in solvers) {
      if (new[] { "madnlp", "fatrop" }.Contains(s.name)) continue;
      message("test_nonregular_point " + s);
      Function solver = ca.nlpsol("mysolver", s.name, nlp, s.opts);
      try {
        solver.call(@in("lbx", vc(-1000), "ubx", vc(1000), "lbg", vc(-1000), "ubg", vc(1000), "x0", 1e-4));
      } catch (SwigException) {
        // ignored, as in python
      }
      if (!new[] { "ipopt", "snopt", "blocksqp", "bonmin", "knitro", "sleqp", "alpaqa", "uno" }.Contains(s.name))
        Check.Equal("SOLVER_RET_NAN", (string)stat(solver, "unified_return_status"), s.ToString());
      Check.True(!success(solver), s.ToString());
    }

    nlp = sxnlp("x", x, "f", sq(x), "g", ca.sqrt(x));
    foreach (S s in solvers) {
      if (new[] { "madnlp", "fatrop" }.Contains(s.name)) continue;
      message("test_nonregular_point " + s);
      Function solver = ca.nlpsol("mysolver", s.name, nlp, s.opts);
      var solver_in = @in("lbx", vc(-10), "ubx", vc(10), "lbg", vc(-2), "ubg", vc(2), "x0", -2);
      try {
        solver.call(solver_in);
        solver.call(solver_in);
      } catch (SwigException) {
        // ignored, as in python
      }
      if (!new[] { "ipopt", "snopt", "bonmin", "knitro", "sleqp", "alpaqa" }.Contains(s.name))
        Check.Equal("SOLVER_RET_NAN", (string)stat(solver, "unified_return_status"), s.ToString());
      Check.True(!success(solver), s.ToString());
    }
  }

  class Interrupt : Callback {
    public Interrupt(string name) {
      construct(name, d("enable_fd", true));
    }

    public override List<DM> eval(IEnumerable<DM> argin) => throw new Exception("KeyboardInterrupt");
  }


  // Python raises KeyboardInterrupt in the callbacks; C# throws an Exception instead
  public void test_iteration_interrupt() {
    var interrupt = new Interrupt("interrupt");
    MX x = MX.sym("x");
    var nlp = mxnlp("x", x, "f", interrupt.Call(x), "g", x);
    foreach (S s in solvers) {
      if (s.name == "madnlp") continue;
      var solver_options = with(s.opts, "error_on_fail", true);
      Function solver = ca.nlpsol("solver", s.name, nlp, solver_options);
      assertInAnyOutput("KeyboardInterrupt", () => solver.call(@in("lbg", -5, "ubg", 5)));
      Check.Throws<Exception>(() => solver.call(@in("lbg", -5, "ubg", 5)));
    }

    foreach (S s in solvers) {
      if (new[] { "worhp", "blocksqp", "knitro", "bonmin", "snopt", "alpaqa", "fatrop", "madnlp" }.Contains(s.name)) continue;
      message("test_iteration_interrupt " + s);
      var opti = ca.Opti();
      MX xo = opti.variable();
      opti.minimize(ca.pow(xo - 1, 4));
      opti.subject_to(-5.0 <= (xo <= 5));
      var opts = new Dictionary<string, object>(s.opts);
      if (s.name == "bonmin") opts["discrete"] = new[] { true };
      opti.solver(s.name, opts);
      opti.solve();
      opti.callback(i => throw new Exception("KeyboardInterrupt"));
      Check.Throws<Exception>(() => opti.solve());

      opts["iteration_callback_ignore_errors"] = false;
      opti.solver(s.name, opts);
      Check.Throws<Exception>(() => opti.solve());

      opts["iteration_callback_ignore_errors"] = true;
      opti.solver(s.name, opts);
      opti.solve();
      // final KeyboardInterrupt-overrides-ignore_errors check: python-only (KeyboardInterrupt semantics)
    }
    keep.Add(interrupt);
  }

  // test_eval_interrupt_stops_ipopt: python-only (KeyboardInterrupt in a Callback is special-cased to abort the solve)

  public void test_discrete() {
    message("rosenbrock, limited-memory hessian approx");
    SX x = SX.sym("x"), y = SX.sym("y");
    var nlp = sxnlp("x", ca.vertcat(x, y), "f", rosen(1.4, x, y));
    foreach (S s in solvers) {
      if (!s.discrete) continue;
      Function solver = ca.nlpsol("mysolver", s.name, nlp, with(s.opts, "discrete", new[] { true, false }));
      var solver_in = @in("lbx", rep(-10, 2), "ubx", rep(10, 2));
      var solver_out = solver.call(solver_in);
      Check.True(success(solver));
      almost(el(solver_out["f"], 0), 0.16, 10, s.name);
      almost(el(solver_out["x"], 0), 1, 6, s.name);
      almost(el(solver_out["x"], 1), 1, 6, s.name);
      checkSerialize(solver, solver_in);
    }
    foreach (S s in solvers) {
      if (!s.discrete) continue;
      Function solver = ca.nlpsol("mysolver", s.name, nlp, with(s.opts, "discrete", new[] { false, true }));
      var solver_in = @in("lbx", rep(-10, 2), "ubx", rep(10, 2));
      var solver_out = solver.call(solver_in);
      Check.True(success(solver));
      almost(el(solver_out["f"], 0), 2e-4, 3, s.name);
      almost(el(solver_out["x"], 0), Math.Sqrt(2), 3, s.name);
      almost(el(solver_out["x"], 1), 2, 3, s.name);
      checkSerialize(solver, solver_in);
    }
  }

  public void test_nan() {
    SX x = SX.sym("x");
    var nlp = sxnlp("x", x, "f", -x, "g", x);
    foreach (S s in solvers) {
      message("test_nan " + s);
      Function solver = ca.nlpsol("mysolver", s.name, nlp, s.opts);
      foreach (string v in new[] { "x", "g" }) {
        string lb = "lb" + v, ub = "ub" + v;
        foreach (var data in new[] { @in(lb, 3, ub, -3), @in(lb, inf, ub, inf), @in(lb, -inf, ub, -inf),
                                     @in(lb, double.NaN), @in(ub, double.NaN) })
          assertInException("Ill-posed", () => solver.call(data));
      }
    }
  }

  public void test_wrongdims() {
    SX x = SX.sym("x", 2);
    var nlp = sxnlp("x", x, "f", -x[0], "g", ca.diag(x));
    foreach (S s in solvers) assertInException("dense vector", () => ca.nlpsol("mysolver", s.name, nlp, s.opts));
    nlp = sxnlp("x", x, "f", -x[0], "g", ca.mtimes(x, x.T()));
    foreach (S s in solvers) assertInException("dense vector", () => ca.nlpsol("mysolver", s.name, nlp, s.opts));
    nlp = sxnlp("x", x, "f", ca.SX(1, 1), "g", x);
    foreach (S s in solvers) assertInException("dense", () => ca.nlpsol("mysolver", s.name, nlp, s.opts));
    nlp = sxnlp("x", x, "f", SX.zeros(0, 0), "g", x);
    foreach (S s in solvers) ca.nlpsol("mysolver", s.name, nlp, s.opts);
    nlp = sxnlp("x", x, "g", x);
    foreach (S s in solvers) ca.nlpsol("mysolver", s.name, nlp, s.opts);
    SX xv = ca.vec(ca.diag(SX.sym("x", 2)));
    nlp = sxnlp("x", xv, "f", ca.mtimes(xv.T(), xv), "g", xv[0]);
    foreach (S s in solvers) assertInException("dense vector", () => ca.nlpsol("mysolver", s.name, nlp, s.opts));
  }

  public void test_initialcond() {
    SX x = SX.sym("x");
    var nlp = sxnlp("x", x, "f", -ca.cos(x), "g", x);
    foreach (S s in solvers) {
      if (s.name == "fatrop") continue;
      Function solver = ca.nlpsol("mysolver", s.name, nlp, s.opts);
      var solver_out = solver.call(@in("x0", vc(6 * pi + 0.01), "lbx", -inf, "ubx", inf, "lbg", -100, "ubg", 100));
      Check.True(success(solver), s.ToString());
      almost(el(solver_out["x"], 0), 6 * pi, 6, s.name);
    }
  }

  public void test_boundsviol() {
    SX x = SX.sym("x");
    var nlp = sxnlp("x", x, "f", sq(x - 1), "g", x);
    foreach (S s in solvers) {
      Function solver = ca.nlpsol("mysolver", s.name, nlp, s.opts);
      Check.Throws<Exception>(() => solver.call(@in("lbx", vc(-10), "ubx", vc(-20), "lbg", vc(-10), "ubg", vc(10))));
    }
    foreach (S s in solvers) {
      Function solver = ca.nlpsol("mysolver", s.name, nlp, s.opts);
      Check.Throws<Exception>(() => solver.call(@in("lbx", vc(-10), "ubx", vc(10), "lbg", vc(-10), "ubg", vc(-20))));
    }
  }

  public void test_IPOPT() {
    SX x = SX.sym("x");
    var nlp = sxnlp("x", x, "f", sq(x - 1), "g", x);
    foreach (S s in solvers) {
      message("test_IPOPT " + s);
      Function solver = ca.nlpsol("mysolver", s.name, nlp, s.opts);
      var solver_in = @in("lbx", vc(-10), "ubx", vc(10), "lbg", vc(-10), "ubg", vc(10));
      var solver_out = solver.call(solver_in);
      Check.True(success(solver));
      almost(el(solver_out["f"], 0), 0, 10, s.name);
      almost(el(solver_out["x"], 0), 1, 9, s.name);
      if (!s.name.Contains("bonmin")) almost(el(solver_out["g"], 0), 1, 9, s.name);
      if (!s.name.Contains("bonmin")) almost(el(solver_out["lam_x"], 0), 0, 8, s.name);
      if (!s.name.Contains("bonmin")) almost(el(solver_out["lam_g"], 0), 0, 9, s.name);
      Check.Equal("SOLVER_RET_SUCCESS", (string)stat(solver, "unified_return_status"));
      if (s.codegen != null) checkCodegen(solver, solver_in, s.codegen, 15, true, false, null);
    }
  }

  public void test_IPOPT_par() {
    SX x = SX.sym("x"), p = SX.sym("p");
    var nlp = sxnlp("x", x, "p", p, "f", sq(x - p), "g", x);
    foreach (S s in solvers) {
      Function solver = ca.nlpsol("mysolver", s.name, nlp, s.opts);
      var solver_in = @in("lbx", vc(-10), "ubx", vc(10), "lbg", vc(-10), "ubg", vc(10), "p", 1);
      var solver_out = solver.call(solver_in);
      Check.True(success(solver));
      almost(el(solver_out["f"], 0), 0, 10, s.name);
      almost(el(solver_out["x"], 0), 1, 9, s.name);
      if (!s.name.Contains("bonmin")) almost(el(solver_out["lam_x"], 0), 0, 8, s.name);
      if (!s.name.Contains("bonmin")) almost(el(solver_out["lam_g"], 0), 0, 9, s.name);
      if (s.codegen != null) checkCodegen(solver, solver_in, s.codegen);
    }
  }

  public void test_IPOPTinf() {
    message("trivial IPOPT, infinity bounds");
    SX x = SX.sym("x");
    var nlp = sxnlp("x", x, "f", sq(x - 1), "g", x);
    foreach (S s in solvers) {
      Function solver = ca.nlpsol("mysolver", s.name, nlp, s.opts);
      var solver_in = @in("lbx", vc(-inf), "ubx", vc(inf), "lbg", vc(-inf), "ubg", vc(inf));
      if (s.name == "worhp") {
        Check.Throws<Exception>(() => solver.call(solver_in));
        return;
      }
      var solver_out = solver.call(solver_in);
      Check.True(success(solver));
      almost(el(solver_out["f"], 0), 0, 10, s.name);
      almost(el(solver_out["x"], 0), 1, 7, s.name);
      if (!s.name.Contains("bonmin")) almost(el(solver_out["lam_x"], 0), 0, 9, s.name);
      if (!s.name.Contains("bonmin")) almost(el(solver_out["lam_g"], 0), 0, 9, s.name);
      if (s.codegen != null) checkCodegen(solver, solver_in, s.codegen);

      Function solver2 = ca.nlpsol("mysolver", s.name, nlp, s.opts);
      var solver_in2 = @in();
      solver_out = solver2.call(solver_in2);
      Check.True(success(solver2));
      almost(el(solver_out["f"], 0), 0, 10, s.name);
      almost(el(solver_out["x"], 0), 1, 7, s.name);
      if (!s.name.Contains("bonmin")) almost(el(solver_out["lam_x"], 0), 0, 9, s.name);
      if (!s.name.Contains("bonmin")) almost(el(solver_out["lam_g"], 0), 0, 9, s.name);
      if (s.codegen != null) checkCodegen(solver2, solver_in2, s.codegen);
    }
  }

  public void test_IPOPTrb() {
    message("rosenbrock, limited-memory hessian approx");
    SX x = SX.sym("x"), y = SX.sym("y");
    var nlp = sxnlp("x", ca.vertcat(x, y), "f", rosen(1, x, y));
    foreach (S s in solvers) {
      if (s.name == "sqpmethod" && s.@in("regularize")) continue;
      if (s.name == "snopt") continue;
      Function solver = ca.nlpsol("mysolver", s.name, nlp, s.opts);
      var solver_in = @in("lbx", rep(-10, 2), "ubx", rep(10, 2));
      var solver_out = solver.call(solver_in);
      Check.True(success(solver));
      almost(el(solver_out["f"], 0), 0, 10, s.name);
      almost(el(solver_out["x"], 0), 1, 6, s.name);
      almost(el(solver_out["x"], 1), 1, 6, s.name);
      if (!s.name.Contains("bonmin")) almost(el(solver_out["lam_x"], 0), 0, 5, s.name);
      if (!s.name.Contains("bonmin")) almost(el(solver_out["lam_x"], 1), 0, 5, s.name);
      if (s.codegen != null) checkCodegen(solver, solver_in, s.codegen);
    }
  }

  public void test_IPOPTrb2() {
    message("rosenbrock, limited-memory hessian approx");
    SX x = SX.sym("x"), y = SX.sym("y");
    var nlp = sxnlp("x", ca.vertcat(x, y), "f", rosen(1, x, y), "g", x + y);
    foreach (S s in solvers) {
      if (s.name == "sqpmethod" && s.@in("regularize")) continue;
      Function solver = ca.nlpsol("mysolver", s.name, nlp, s.opts);
      var solver_in = @in("lbx", rep(-10, 2), "ubx", rep(10, 2), "lbg", vc(-10), "ubg", vc(10));
      var solver_out = solver.call(solver_in);
      if (!s.name.Contains("knitro")) Check.True(success(solver));
      int digits = 6;
      almost(el(solver_out["f"], 0), 0, digits, s.name);
      almost(el(solver_out["x"], 0), 1, digits, s.name);
      almost(el(solver_out["x"], 1), 1, digits, s.name);
      if (!s.name.Contains("bonmin")) almost(el(solver_out["lam_x"], 0), 0, 5, s.name);
      if (!s.name.Contains("bonmin")) almost(el(solver_out["lam_x"], 1), 0, 5, s.name);
      if (!s.name.Contains("bonmin")) almost(el(solver_out["lam_g"], 0), 0, 5, s.name);
      if (s.codegen != null) checkCodegen(solver, solver_in, s.codegen);
    }
  }

  public void test_IPOPTrbf() {
    message("rosenbrock fixed, limited-memory hessian approx");
    SX x = SX.sym("x"), y = SX.sym("y");
    var nlp = sxnlp("x", ca.vertcat(x, y), "f", rosen(1, x, y), "g", x + y);
    foreach (S s in solvers) {
      Function solver = ca.nlpsol("mysolver", s.name, nlp, s.opts);
      var solver_in = @in("x0", vc(0, 1), "lbx", vc(-10, 1), "ubx", vc(10, 1), "lbg", vc(-10), "ubg", vc(10));
      if (s.@in("qrqp")) solver_in["x0"] = vc(0.6, 1);
      var solver_out = solver.call(solver_in);
      Check.True(success(solver));
      almost(el(solver_out["f"], 0), 0, 10, s.name);
      almost(el(solver_out["x"], 0), 1, 7, s.name);
      almost(el(solver_out["x"], 1), 1, 7, s.name);
      if (!s.name.Contains("stabilizedsqp")) {
        if (!s.name.Contains("bonmin")) almost(el(solver_out["lam_x"], 0), 0, 5, s.name);
        if (!s.name.Contains("bonmin")) almost(el(solver_out["lam_x"], 1), 0, 5, s.name);
        if (!s.name.Contains("bonmin")) almost(el(solver_out["lam_g"], 0), 0, 5, s.name);
      }
      if (s.codegen != null) checkCodegen(solver, solver_in, s.codegen);
    }
  }

  public void test_knitro_options_file() {
    Check.Assume(hasNlpsol("knitro"), "no nlpsol knitro");
    MX x = MX.sym("x");
    var nlp = mxnlp("x", x, "f", sq(x));
    File.WriteAllText("my_knitro.opt", "maxit            42\n");
    Function solver = ca.nlpsol("solver", "knitro", nlp, d("options_file", "my_knitro.opt"));
    assertOutput(new[] { "maxit:", "42" }, new string[0], () => solver.call(@in()));
  }

  public void test_ipopt_linear_solvers() {
    Check.Assume(hasNlpsol("ipopt"), "no nlpsol ipopt");
    var lsolvers = new List<string> { "mumps" };
    if (CasadiMeta.feature_list().Contains("spral")) lsolvers.Add("spral");
    if (hasLinsol("ma27") && !RuntimeInformation.IsOSPlatform(OSPlatform.OSX)) lsolvers.Add("ma27");
    SX x = SX.sym("x"), y = SX.sym("y");
    var nlp = sxnlp("x", ca.vertcat(x, y), "f", rosen(1, x, y), "g", sq(x) + sq(y));
    double c_r = 4.56748075136258e-02;
    double[] x_r = { 7.86415156987791e-01, 6.17698316967954e-01 };
    foreach (string ls in lsolvers) {
      Function solver = ca.nlpsol("mysolver", "ipopt", nlp, d("ipopt.linear_solver", ls));
      var solver_out = solver.call(@in("x0", vc(0.5, 0.5), "lbx", rep(-10, 2), "ubx", rep(10, 2), "lbg", vc(0), "ubg", vc(1)));
      int digits = 7;
      almost(el(solver_out["f"], 0), c_r, digits, ls);
      almost(el(solver_out["x"], 0), x_r[0], digits, ls);
      almost(el(solver_out["x"], 1), x_r[1], digits, ls);
    }
  }

  public void test_uno_linear_solvers() {
    Check.Assume(hasNlpsol("uno"), "no nlpsol uno");
    var lsolvers = new List<string> { "MUMPS" };
    if (CasadiMeta.feature_list().Contains("spral")) lsolvers.Add("SSIDS");
    if (hasLinsol("ma27") && !RuntimeInformation.IsOSPlatform(OSPlatform.OSX)) lsolvers.Add("MA27");
    SX x = SX.sym("x"), y = SX.sym("y");
    var nlp = sxnlp("x", ca.vertcat(x, y), "f", rosen(1, x, y), "g", sq(x) + sq(y));
    double c_r = 4.56748075136258e-02;
    double[] x_r = { 7.86415156987791e-01, 6.17698316967954e-01 };
    foreach (string ls in lsolvers) {
      Function mysolver = ca.nlpsol("mysolver", "uno", nlp, d("print_time", false, "uno", d("preset", "ipopt",
          "linear_solver", ls, "primal_tolerance", 1e-8, "dual_tolerance", 1e-8)));
      var solver_out = mysolver.call(@in("x0", vc(0.5, 0.5), "lbx", rep(-10, 2), "ubx", rep(10, 2), "lbg", vc(0), "ubg", vc(1)));
      int digits = 6;
      almost(el(solver_out["f"], 0), c_r, digits, ls);
      almost(el(solver_out["x"], 0), x_r[0], digits, ls);
      almost(el(solver_out["x"], 1), x_r[1], digits, ls);
    }
  }

  public void test_ipopt_stats_per_call() {
    Check.Assume(hasNlpsol("ipopt"), "no nlpsol ipopt");
    MX x = MX.sym("x");
    Function solver = ca.nlpsol("solver", "ipopt", mxnlp("x", x, "f", sq(x)));
    solver.call(@in());
    var s1 = solver.stats();
    solver.call(@in());
    var s2 = solver.stats();
    Check.Equal(s1["n_call_nlp_f"], s2["n_call_nlp_f"]);
  }

  public void test_warmstart() {
    SX x = SX.sym("x"), y = SX.sym("y");
    var nlp = sxnlp("x", ca.vertcat(x, y), "f", rosen(1, x, y), "g", sq(x) + sq(y));
    double c_r = 4.56748075136258e-02;
    double[] x_r = { 7.86415156987791e-01, 6.17698316967954e-01 };
    foreach (S s in solvers) {
      Function solver = ca.nlpsol("mysolver", s.name, nlp, s.opts);
      var solver_in = @in("x0", vc(0.5, 0.5), "lbx", rep(-10, 2), "ubx", rep(10, 2), "lbg", vc(0), "ubg", vc(1));
      var solver_out = solver.call(solver_in);
      var oldsolver_out = solver_out;
      int digits = 5;
      almost(el(solver_out["f"], 0), c_r, digits, s.name);
      almost(el(solver_out["x"], 0), x_r[0], digits, s.name);
      almost(el(solver_out["x"], 1), x_r[1], digits, s.name);
      if (!s.name.Contains("bonmin")) almost(el(solver_out["lam_x"], 0), 0, s.name == "snopt" ? 5 : 8, s.name);
      if (!s.name.Contains("bonmin")) almost(el(solver_out["lam_x"], 1), 0, s.name == "snopt" ? 5 : 8, s.name);
      if (!new[] { "bonmin", "sleqp" }.Contains(s.name)) almost(el(solver_out["lam_g"], 0), 0.12149655447670, 6, s.name);
      if (s.codegen != null) checkCodegen(solver, solver_in, s.codegen, codegen_check_digits);

      message(":warmstart");
      if (s.name.Contains("ipopt")) {
        var options = with(s.opts, "ipopt.warm_start_init_point", "yes", "ipopt.warm_start_bound_push", 1e-6,
            "ipopt.warm_start_slack_bound_push", 1e-6, "ipopt.warm_start_mult_bound_push", 1e-6, "ipopt.mu_init", 1e-6);
        solver = ca.nlpsol("mysolver", s.name, nlp, options);
        solver_in["lbx"] = vc(rep(-10, 2));
        solver_in["ubx"] = vc(rep(10, 2));
        solver_in["lbg"] = vc(0);
        solver_in["ubg"] = vc(1);
        solver_in["x0"] = oldsolver_out["x"];
        solver_in["lam_g0"] = oldsolver_out["lam_g"];
        if (!s.name.Contains("bonmin")) solver_in["lam_x0"] = oldsolver_out["lam_x"];
        solver.call(solver_in);
        if (s.codegen != null) checkCodegen(solver, solver_in, s.codegen, codegen_check_digits);
      }
    }
  }

  public void test_IPOPTrhb2_gen() {
    message("rosenbrock, exact hessian generated, constrained");
    SX x = SX.sym("x"), y = SX.sym("y");
    var nlp = sxnlp("x", ca.vertcat(x, y), "f", rosen(1, x, y), "g", sq(x) + sq(y));
    double c_r = 4.56748075136258e-02;
    double[] x_r = { 7.86415156987791e-01, 6.17698316967954e-01 };
    foreach (S s in solvers) {
      if (s.name.Contains("madnlp")) continue;
      Function solver = ca.nlpsol("mysolver", s.name, nlp, s.opts);
      var solver_in = @in("x0", vc(0.5, 0.5), "lbx", rep(-10, 2), "ubx", rep(10, 2), "lbg", vc(0), "ubg", vc(1));
      var solver_out = solver.call(solver_in);
      int digits = 5;
      almost(el(solver_out["f"], 0), c_r, digits, s.name + el(solver_out["f"], 0) + ":" + c_r);
      almost(el(solver_out["x"], 0), x_r[0], digits, s.name);
      almost(el(solver_out["x"], 1), x_r[1], digits, s.name);
      if (!s.name.Contains("bonmin")) almost(el(solver_out["lam_x"], 0), 0, s.name == "snopt" ? 5 : 8, s.name);
      if (!s.name.Contains("bonmin")) almost(el(solver_out["lam_x"], 1), 0, s.name == "snopt" ? 5 : 8, s.name);
      if (!new[] { "bonmin", "sleqp" }.Contains(s.name)) almost(el(solver_out["lam_g"], 0), 0.12149655447670, 6, s.name);
      checkSerialize(solver, solver_in);
      if (s.codegen != null) checkCodegen(solver, solver_in, s.codegen, codegen_check_digits);
    }
  }

  public void test_jacG_empty() {
    SX x = SX.sym("x"), y = SX.sym("y");
    var nlp = sxnlp("x", ca.vertcat(x, y), "f", rosen(1, x, y), "g", 1);
    foreach (S s in solvers) {
      if (new[] { "worhp", "sqpmethod", "snopt", "fatrop", "madnlp" }.Contains(s.name)) continue;
      Function solver = ca.nlpsol("mysolver", s.name, nlp, s.opts);
      var solver_in = @in("x0", vc(0.5, 0.5), "lbx", rep(-10, 2), "ubx", rep(10, 2), "lbg", vc(0), "ubg", vc(2));
      var solver_out = solver.call(solver_in);
      int digits = 5;
      chk(solver_out["f"], vc(0), digits, s.name);
      chk(solver_out["x"], vc(1, 1), digits, s.name);
      if (!s.name.Contains("bonmin")) chk(solver_out["lam_x"], vc(0, 0), digits, s.name);
      if (!s.name.Contains("bonmin")) chk(solver_out["lam_g"], vc(0), digits, s.name);
      if (s.codegen != null) {
        if (s.name == "ipopt") continue;
        checkCodegen(solver, solver_in, s.codegen);
      }
    }
  }

  public void test_IPOPTrhb2_gen_par() {
    message("rosenbrock, exact hessian generated, constrained, parametric");
    SX x = SX.sym("x"), y = SX.sym("y"), p = SX.sym("p");
    var nlp = sxnlp("x", ca.vertcat(x, y), "p", p, "f", rosen(p, x, y), "g", sq(x) + sq(y));
    double c_r = 4.56748075136258e-02;
    double[] x_r = { 7.86415156987791e-01, 6.17698316967954e-01 };
    foreach (S s in solvers) {
      Function solver = ca.nlpsol("mysolver", s.name, nlp, s.opts);
      var solver_in = @in("x0", vc(0.5, 0.5), "lbx", rep(-10, 2), "ubx", rep(10, 2), "lbg", vc(0), "ubg", vc(1), "p", vc(1));
      if (s.codegen != null) checkCodegen(solver, solver_in, s.codegen, codegen_check_digits);
      var solver_out = solver.call(solver_in);
      int digits = 5;
      almost(el(solver_out["f"], 0), c_r, digits, s.name);
      almost(el(solver_out["x"], 0), x_r[0], digits, s.name);
      almost(el(solver_out["x"], 1), x_r[1], digits, s.name);
      if (!s.name.Contains("bonmin")) almost(el(solver_out["lam_x"], 0), 0, s.name == "snopt" ? 5 : 8, s.name);
      if (!s.name.Contains("bonmin")) almost(el(solver_out["lam_x"], 1), 0, s.name == "snopt" ? 5 : 8, s.name);
      if (!new[] { "bonmin", "sleqp" }.Contains(s.name)) almost(el(solver_out["lam_g"], 0), 0.12149655447670, 6, s.name);
    }
  }

  public void test_IPOPTrhb_gen() {
    message("rosenbrock, exact hessian generated");
    SX x = SX.sym("x"), y = SX.sym("y");
    var nlp = sxnlp("x", ca.vertcat(x, y), "f", rosen(1, x, y));
    foreach (S s in solvers) {
      if (s.name == "sqpmethod" && s.@in("regularize")) continue;
      if (s.name == "snopt") continue;
      Function solver = ca.nlpsol("mysolver", s.name, nlp, s.opts);
      var solver_in = @in("lbx", rep(-10, 2), "ubx", rep(10, 2));
      var solver_out = solver.call(solver_in);
      int xdig = 6;
      if (s.name == "uno" && s.@in("LBFGS")) xdig = 7;
      almost(el(solver_out["f"], 0), 0, 10, s.name);
      almost(el(solver_out["x"], 0), 1, xdig, s.name);
      almost(el(solver_out["x"], 1), 1, xdig, s.name);
      if (!s.name.Contains("bonmin")) almost(el(solver_out["lam_x"], 0), 0, 8, s.name);
      if (!s.name.Contains("bonmin")) almost(el(solver_out["lam_x"], 1), 0, 8, s.name);
      if (s.codegen != null) checkCodegen(solver, solver_in, s.codegen);
    }
  }

  public void test_IPOPTrhb_gen_xnonfree() {
    message("rosenbrock, exact hessian generated, non-free x");
    SX x = SX.sym("x"), y = SX.sym("y");
    var nlp = sxnlp("x", ca.vertcat(x, y), "f", rosen(1, x, y));
    foreach (S s in solvers) {
      if (s.name == "sqpmethod" && s.@in("regularize")) continue;
      if (s.name == "snopt") continue;
      Function solver = ca.nlpsol("mysolver", s.name, nlp, s.opts);
      var solver_in = @in("lbx", vc(1, -10), "ubx", vc(1, 10));
      var solver_out = solver.call(solver_in);
      almost(el(solver_out["f"], 0), 0, 10, s.name);
      almost(el(solver_out["x"], 0), 1, 9, s.name);
      almost(el(solver_out["x"], 1), 1, 6, s.name);
      if (!s.name.Contains("bonmin")) almost(el(solver_out["lam_x"], 0), 0, 6, s.name);
      if (!s.name.Contains("bonmin")) almost(el(solver_out["lam_x"], 1), 0, 6, s.name);
      if (s.codegen != null) checkCodegen(solver, solver_in, s.codegen, codegen_check_digits);
    }
  }

  public void test_IPOPTrhb_gen_par() {
    message("rosenbrock, exact hessian generated, parametric");
    SX x = SX.sym("x"), y = SX.sym("y"), p = SX.sym("p");
    var nlp = sxnlp("x", ca.vertcat(x, y), "p", p, "f", rosen(p, x, y));
    foreach (S s in solvers) {
      if (s.name == "sqpmethod" && s.@in("regularize")) continue;
      if (s.name == "snopt") continue;
      Function solver = ca.nlpsol("mysolver", s.name, nlp, s.opts);
      var solver_in = @in("lbx", rep(-10, 2), "ubx", rep(10, 2), "p", 1);
      var solver_out = solver.call(solver_in);
      int xdig = 6;
      if (s.name == "uno" && s.@in("LBFGS")) xdig = 7;
      almost(el(solver_out["f"], 0), 0, 10, s.name);
      almost(el(solver_out["x"], 0), 1, xdig, s.name);
      almost(el(solver_out["x"], 1), 1, xdig, s.name);
      if (s.codegen != null) checkCodegen(solver, solver_in, s.codegen);
    }
  }

  public void test_IPOPTnorm() {
    message("IPOPT min ||x||^2_2");
    int N = 10;
    MX x = MX.sym("x", N);
    DM x0 = ca.linspace(ca.DM(0.0), ca.DM(1.0), N);
    MX X0 = x0;
    MX dd = x - X0;
    var nlp = mxnlp("x", x, "f", ca.dot(dd, dd), "g", 2 * x);
    foreach (S s in solvers) {
      Function solver = ca.nlpsol("mysolver", s.name, nlp, s.opts);
      var solver_in = @in("lbx", rep(-10, N), "ubx", rep(10, N), "lbg", rep(-10, N), "ubg", rep(10, N));
      var solver_out = solver.call(solver_in);
      almost(el(solver_out["f"], 0), 0, 10, s.name);
      chk(solver_out["x"], x0, 8, s.name);
      if (!s.name.Contains("bonmin")) chk(solver_out["lam_x"], DM.zeros(10, 1), 8, s.name);
      if (!s.name.Contains("bonmin")) almost(el(solver_out["lam_g"], 1), 0, 8, s.name);
      if (s.codegen != null) checkCodegen(solver, solver_in, s.codegen);
    }
  }

  public void test_IPOPTnoc() {
    message("trivial IPOPT, no constraints");
    SX x = SX.sym("x");
    var nlp = sxnlp("x", x, "f", sq(x - 1));
    foreach (S s in solvers) {
      if (s.name == "snopt") continue;
      Function solver = ca.nlpsol("mysolver", s.name, nlp, s.opts);
      var solver_in = @in("lbx", vc(-10), "ubx", vc(10));
      var solver_out = solver.call(solver_in);
      almost(el(solver_out["f"], 0), 0, 10, s.name);
      almost(el(solver_out["x"], 0), 1, 7, s.name);
      if (s.codegen != null) checkCodegen(solver, solver_in, s.codegen);
    }
  }

  public void test_IPOPTmx() {
    message("trivial IPOPT, using MX");
    MX x = MX.sym("x");
    var nlp = mxnlp("x", x, "f", sq(x - 1), "g", 2 * x);
    foreach (S s in solvers) {
      Function solver = ca.nlpsol("mysolver", s.name, nlp, s.opts);
      var solver_in = @in("lbx", vc(-10), "ubx", vc(10), "lbg", vc(-10), "ubg", vc(10));
      var solver_out = solver.call(solver_in);
      almost(el(solver_out["f"], 0), 0, 10, s.name);
      almost(el(solver_out["x"], 0), 1, 9, s.name);
      if (s.codegen != null) checkCodegen(solver, solver_in, s.codegen);
    }
  }

  public void test_IPOPTc() {
    message("trivial, overconstrained");
    SX x = SX.sym("x");
    var nlp = sxnlp("x", x, "f", sq(x - 1), "g", ca.vertcat(x, x, x));
    foreach (S s in solvers) {
      Function solver = ca.nlpsol("mysolver", s.name, nlp, s.opts);
      var solver_in = @in("lbx", vc(-10), "ubx", vc(10), "lbg", vc(-10, -10, -10), "ubg", vc(10, 10, 10));
      var solver_out = solver.call(solver_in);
      almost(el(solver_out["f"], 0), 0, 9, s.name);
      almost(el(solver_out["x"], 0), 1, 5, s.name);
      if (s.codegen != null) checkCodegen(solver, solver_in, s.codegen);
    }
  }

  public void test_IPOPTc2() {
    message("trivial2, overconstrained");
    SX x = SX.sym("x");
    var nlp = sxnlp("x", x, "f", sq(x - 1), "g", ca.vertcat(x, x, x + x));
    foreach (S s in solvers) {
      Function solver = ca.nlpsol("mysolver", s.name, nlp, s.opts);
      var solver_in = @in("lbx", vc(-10), "ubx", vc(10), "lbg", vc(-10, -10, -10), "ubg", vc(10, 10, 10));
      var solver_out = solver.call(solver_in);
      almost(el(solver_out["f"], 0), 0, 10, s.name);
      almost(el(solver_out["x"], 0), 1, 8, s.name);
      if (s.codegen != null) checkCodegen(solver, solver_in, s.codegen);
    }
  }

  public void test_IPOPTcmx() {
    message("trivial , overconstrained, using MX");
    MX x = MX.sym("x");
    var nlp = mxnlp("x", x, "f", sq(x - 1), "g", ca.vertcat(2 * x, 3 * x, 4 * x));
    foreach (S s in solvers) {
      Function solver = ca.nlpsol("mysolver", s.name, nlp, s.opts);
      var solver_in = @in("lbx", vc(-10), "ubx", vc(10), "lbg", vc(-10, -10, -10), "ubg", vc(10, 10, 10));
      var solver_out = solver.call(solver_in);
      almost(el(solver_out["f"], 0), 0, 9, s.name);
      almost(el(solver_out["x"], 0), 1, 8, s.name);
      if (s.codegen != null) checkCodegen(solver, solver_in, s.codegen);
    }
  }

  public void test_IPOPTdeg() {
    message("degenerate optimization IPOPT");
    SX x = SX.sym("x"), y = SX.sym("y");
    var nlp = sxnlp("x", ca.vertcat(x, y), "f", 0, "g", ca.vertcat(x - y, x));
    foreach (S s in solvers) {
      Function solver = ca.nlpsol("mysolver", s.name, nlp, s.opts);
      var solver_in = @in("lbx", vc(-10, -10), "ubx", vc(10, 10), "lbg", vc(0, 3), "ubg", vc(0, 3));
      var solver_out = solver.call(solver_in);
      almost(el(solver_out["x"], 0), el(solver_out["x"], 1), s.@in("sqic") ? 4 : 10, "IPOPT");
      if (s.codegen != null) checkCodegen(solver, solver_in, s.codegen);
    }
  }

  public void test_IPOPTdegc() {
    message("degenerate optimization IPOPT, overconstrained");
    SX x = SX.sym("x"), y = SX.sym("y");
    var nlp = sxnlp("x", ca.vertcat(x, y), "f", 0, "g", ca.vertcat(x - y, x, x + y));
    foreach (S s in solvers) {
      Function solver = ca.nlpsol("mysolver", s.name, nlp, s.opts);
      var solver_in = @in("lbx", vc(-10, -10), "ubx", vc(10, 10), "lbg", vc(0, 3, -10), "ubg", vc(0, 3, 10));
      var solver_out = solver.call(solver_in);
      almost(el(solver_out["x"], 0), el(solver_out["x"], 1), s.@in("sqic") ? 4 : 10, "IPOPT");
      if (s.codegen != null) checkCodegen(solver, solver_in, s.codegen);
    }
  }

  public void test_XfreeChange() {
    message("Change in X settings");
    SX x = SX.sym("x"), y = SX.sym("y");
    var nlp = sxnlp("x", ca.vertcat(x, y), "f", rosen(1, x, y), "g", x + y);
    foreach (S s in solvers) {
      Function solver = ca.nlpsol("mysolver", s.name, nlp, s.opts);
      var solver_in = @in("x0", vc(0, 1), "lbx", vc(-10, -10), "ubx", vc(10, 10), "lbg", vc(-10), "ubg", vc(10));
      if (s.@in("qrqp")) solver_in["x0"] = vc(0.6, 1);
      solver.call(solver_in);
      solver_in["lbx"] = vc(-10, 1);
      solver_in["ubx"] = vc(10, 1);
      solver_in["lbg"] = vc(-10);
      solver_in["ubg"] = vc(10);
      var solver_out = solver.call(solver_in);
      almost(el(solver_out["f"], 0), 0, 9, s.name);
      almost(el(solver_out["x"], 0), 1, 6, s.name);
      almost(el(solver_out["x"], 1), 1, 6, s.name);
      if (s.codegen != null) checkCodegen(solver, solver_in, s.codegen);
    }
  }

  public void test_activeLBX() {
    message("active LBX");
    SX x = SX.sym("x"), y = SX.sym("y");
    var nlp = sxnlp("x", ca.vertcat(x, y), "f", rosen(1, x, y), "g", x + y);
    foreach (S s in solvers) {
      Function solver = ca.nlpsol("mysolver", s.name, nlp, s.opts);
      var solver_in = @in("x0", vc(0, 1), "lbx", vc(-10, 1.2), "ubx", vc(10, 2), "lbg", vc(-10), "ubg", vc(10));
      if (s.@in("qrqp")) solver_in["x0"] = vc(0.5, 1);
      var solver_out = solver.call(solver_in);
      if (el(solver_out["x"], 0) < 0) {
        almost(el(solver_out["f"], 0), 4.3817250416084308, 6, s.name);
        almost(el(solver_out["x"], 0), -1.0910624688699295, 6, s.name);
        almost(el(solver_out["x"], 1), 1.2, 5, s.name);
        if (!s.name.Contains("bonmin")) almost(el(solver_out["lam_x"], 0), 0, s.name == "stabilizedsqp" ? 5 : 8, s.ToString());
        if (!s.name.Contains("bonmin")) almost(el(solver_out["lam_x"], 1), -1.9165378046901287, 4, s.name);
        if (!s.name.Contains("bonmin")) almost(el(solver_out["lam_g"], 0), 0, 8, s.name);
      } else {
        almost(el(solver_out["f"], 0), 9.0908263002590e-3, 6, s.name);
        almost(el(solver_out["x"], 0), 1.0952466252248, 6, s.name);
        almost(el(solver_out["x"], 1), 1.2, 5, s.name);
        if (!s.name.Contains("bonmin")) almost(el(solver_out["lam_x"], 0), 0, s.name == "stabilizedsqp" ? 5 : 8, s.ToString());
        if (!s.name.Contains("bonmin")) almost(el(solver_out["lam_x"], 1), -8.6963632695079e-2, 4, s.name);
        if (!s.name.Contains("bonmin")) almost(el(solver_out["lam_g"], 0), 0, 8, s.name);
      }
      if (s.codegen != null) checkCodegen(solver, solver_in, s.codegen);
    }
  }

  public void test_activeLBG() {
    message("active LBG");
    SX x = SX.sym("x"), y = SX.sym("y");
    var nlp = sxnlp("x", ca.vertcat(x, y), "f", rosen(1, x, y), "g", x + y);
    foreach (S s in solvers) {
      Function solver = ca.nlpsol("mysolver", s.name, nlp, s.opts);
      var solver_in = @in("x0", vc(0, 1), "lbx", vc(-10, -10), "ubx", vc(10, 10), "lbg", vc(2.2), "ubg", vc(10));
      if (s.@in("qrqp")) solver_in["x0"] = vc(1.5, 1);
      var solver_out = solver.call(solver_in);
      almost(el(solver_out["f"], 0), 4.252906468284e-3, 6, s.name);
      almost(el(solver_out["x"], 0), 1.065181061847138, 6, s.name);
      almost(el(solver_out["x"], 1), 1.1348189166291160, 6, s.name);
      if (!s.name.Contains("bonmin")) almost(el(solver_out["lam_x"], 0), 0, s.name == "snopt" ? 6 : 8, s.name);
      if (!s.name.Contains("bonmin")) almost(el(solver_out["lam_x"], 1), 0, 4, s.name);
      if (!s.name.Contains("bonmin")) almost(el(solver_out["lam_g"], 0), -4.1644422845712e-2, 3, s.name);
      if (s.codegen != null) checkCodegen(solver, solver_in, s.codegen);
    }
  }

  public void test_activeUBG() {
    message("active UBG");
    SX x = SX.sym("x"), y = SX.sym("y");
    var nlp = sxnlp("x", ca.vertcat(x, y), "f", rosen(1, x, y), "g", x + y);
    foreach (S s in solvers) {
      Function solver = ca.nlpsol("mysolver", s.name, nlp, s.opts);
      var solver_in = @in("x0", vc(0, 1), "lbx", vc(-10, -10), "ubx", vc(10, 10), "lbg", vc(0), "ubg", vc(1.8));
      if (s.@in("qrqp")) solver_in["x0"] = vc(1.5, 1);
      var solver_out = solver.call(solver_in);
      almost(el(solver_out["f"], 0), 4.64801220074552e-3, 6, s.name);
      almost(el(solver_out["x"], 0), 9.318651964592811e-1, 5, s.name);
      almost(el(solver_out["x"], 1), 8.68134821123689e-1, 5, s.name);
      if (!s.name.Contains("bonmin")) almost(el(solver_out["lam_x"], 0), 0, s.name == "snopt" ? 5 : 8, s.name);
      if (!s.name.Contains("bonmin")) almost(el(solver_out["lam_x"], 1), 0, 4, s.name);
      if (!s.name.Contains("bonmin")) almost(el(solver_out["lam_g"], 0), 4.75846495145007e-2, 5, s.name);
      if (s.codegen != null) checkCodegen(solver, solver_in, s.codegen);
    }
  }

  public void test_activeUBX() {
    message("active UBX");
    SX x = SX.sym("x"), y = SX.sym("y");
    var nlp = sxnlp("x", ca.vertcat(x, y), "f", rosen(1, x, y), "g", x + y);
    foreach (S s in solvers) {
      Function solver = ca.nlpsol("mysolver", s.name, nlp, s.opts);
      var solver_in = @in("x0", vc(0, 1), "lbx", vc(-10, 0), "ubx", vc(10, 0.9), "lbg", vc(-10), "ubg", vc(10));
      var solver_out = solver.call(solver_in);
      almost(el(solver_out["f"], 0), 2.626109721583e-3, 6, s.name);
      almost(el(solver_out["x"], 0), 9.4882542279172277e-01, 6, s.name);
      almost(el(solver_out["x"], 1), 0.9, 6, s.name);
      if (!s.name.Contains("bonmin")) almost(el(solver_out["lam_x"], 0), 0, s.name == "snopt" ? 6 : 8, s.name);
      if (!s.name.Contains("bonmin")) almost(el(solver_out["lam_x"], 1), 5.39346608659e-2, 4, s.name);
      if (!s.name.Contains("bonmin")) almost(el(solver_out["lam_g"], 0), 0, 8, s.name);
      if (s.codegen != null) checkCodegen(solver, solver_in, s.codegen);
    }
  }

  public void test_QP() {
    message("QP");
    int N = 50;
    SX x = SX.sym("x", N);
    DM x0 = range(N);
    DM H = ca.diag((DM)Enumerable.Range(1, N).Select(i => (double)i).ToArray());
    SX e = x - x0;
    SX obj = 0.5 * ca.mtimes(ca.mtimes(e.T(), H), e);
    var nlp = sxnlp("x", x, "f", obj);
    foreach (S s in solvers) {
      if (s.name == "snopt" || s.name == "fatrop") continue;
      if (s.name == "sqpmethod" && s.@in("limited-memory")) continue;
      Function solver = ca.nlpsol("mysolver", s.name, nlp, s.opts);
      var solver_out = solver.call(@in("lbx", -1000, "ubx", 1000));
      chk(solver_out["x"], x0, 2, s.name);
      almost(el(solver_out["f"], 0), 0, 3, s.name);
      if (!s.name.Contains("bonmin")) chk(solver_out["lam_x"], DM.zeros(N, 1), 4, s.name);
    }
  }

  public void test_QP2() {
    DM H = new double[,] { { 1, -1 }, { -1, 2 } };
    DM G = vc(-2, -6);
    DM A = new double[,] { { 1, 1 }, { -1, 2 }, { 2, 1 } };
    DM LBA = vc(rep(-inf, 3)), UBA = vc(2, 2, 3);
    DM LBX = vc(0.5, 0), UBX = vc(0.5, inf);
    SX x = SX.sym("x", 2);
    var nlp = sxnlp("x", x, "f", 0.5 * ca.mtimes(ca.mtimes(x.T(), H), x) + ca.mtimes(G.T(), x), "g", ca.mtimes(A, x));
    foreach (S s in solvers) {
      var options = new Dictionary<string, object>(s.opts);
      if (s.name.Contains("madnlp")) continue;
      if (s.name.Contains("ipopt")) options["ipopt.fixed_variable_treatment"] = "make_constraint";
      for (int rep = 0; rep < 2; rep++) {
        Function solver = ca.nlpsol("mysolver", s.name, nlp, options);
        var solver_in = @in("lbx", LBX, "ubx", UBX, "lbg", LBA, "ubg", UBA);
        if (s.@in("sqic") || s.@in("daqp")) break;
        var solver_out = solver.call(solver_in);
        almost(el(solver_out["x"], 0), 0.5, 6, s.name);
        almost(el(solver_out["x"], 1), 1.25, 6, s.name);
        if (!s.name.Contains("bonmin")) almost(el(solver_out["lam_x"], 0), 4.75, 6, s.name);
        if (!s.name.Contains("bonmin")) almost(el(solver_out["lam_x"], 1), 0, 6, s.name);
        if (!s.name.Contains("bonmin")) chk(solver_out["lam_g"], vc(0, 2, 0), 6, s.name);
        almost(el(solver_out["f"], 0), -7.4375, 6, s.name);
        if (s.codegen != null) checkCodegen(solver, solver_in, s.codegen);
      }
    }
  }

  public void test_QP2_unconvex() {
    DM H = new double[,] { { 1, -1 }, { -1, -2 } };
    DM G = vc(-2, -6);
    DM A = new double[,] { { 1, 1 }, { -1, 2 }, { 2, 1 } };
    DM LBA = vc(rep(-inf, 3)), UBA = vc(2, 2, 3);
    DM LBX = vc(0, 0), UBX = vc(inf, inf);
    SX x = SX.sym("x", 2);
    var nlp = sxnlp("x", x, "f", 0.5 * ca.mtimes(ca.mtimes(x.T(), H), x) + ca.mtimes(G.T(), x), "g", ca.mtimes(A, x));
    foreach (S s in solvers) {
      var options = new Dictionary<string, object>(s.opts);
      if (s.name.Contains("ipopt")) options["ipopt.fixed_variable_treatment"] = "make_constraint";
      var solver_in = @in("lbx", LBX, "ubx", UBX, "lbg", LBA, "ubg", UBA);
      if (s.@in("qrqp")) {
        solver_in["x0"] = vc(1, 1);
        solver_in["lam_g0"] = vc(1, 1, 0);
      }
      if (s.@in("daqp")) continue;
      for (int rep = 0; rep < 2; rep++) {
        Function solver = ca.nlpsol("mysolver", s.name, nlp, options);
        string msg = solver.ToString();
        var solver_out = solver.call(solver_in);
        almost(el(solver_out["x"], 0), 2.0 / 3, 6, msg);
        almost(el(solver_out["x"], 1), 4.0 / 3, 6, msg);
        if (!s.name.Contains("bonmin")) almost(el(solver_out["lam_x"], 0), 0, 6, msg);
        if (!s.name.Contains("bonmin")) almost(el(solver_out["lam_x"], 1), 0, 6, msg);
        if (!s.name.Contains("bonmin") && !s.name.Contains("worhp"))
          chk(solver_out["lam_g"], vc(4 + 8.0 / 9, 20.0 / 9, 0), 6, msg);
        almost(el(solver_out["f"], 0), -10 - 16.0 / 9, 6, msg);
        if (s.codegen != null) checkCodegen(solver, solver_in, s.codegen);
      }
    }
  }

  public void test_bug() {
    MX x = MX.sym("x", 3), y = MX.sym("y", 2);
    var f = ca.Function("f", new[] { x, y }, new[] { ca.MX(1.0) });
    MX aa = MX.sym("aa", 5);
    MX a = aa[0..3], b = aa[3..5];
    MX f_call = f.Call(a, b);
    var nlp = mxnlp("x", aa, "f", f_call);
    foreach (S s in solvers) {
      if (s.name.Contains("worhp") || s.name == "snopt") continue;
      ca.nlpsol("mysolver", s.name, nlp, s.opts);
    }
  }

  public void test_missing_symbols() {
    MX x = MX.sym("x"), p = MX.sym("p");
    foreach (S s in solvers)
      assertInException("[p] are free", () => ca.nlpsol("solver", s.name, mxnlp("x", x, "f", sq(x - p)), s.opts));
  }

  public void test_no_success() {
    Check.Assume(hasNlpsol("ipopt"), "no nlpsol ipopt");
    SX x = SX.sym("x"), y = SX.sym("y");
    SX f = rosen(1, x, y);
    for (int k = 0; k < solvers.Count; k++) {
      Function solver = ca.nlpsol("solver", "ipopt", sxnlp("x", ca.vertcat(x, y), "f", f, "g", ca.vertcat(x + 1, x - 2)));
      solver.call(@in("x0", 0, "lbg", 0, "ubg", 0));
      Check.True(!success(solver));
      Function solver2 = ca.nlpsol("solver", "ipopt", sxnlp("x", ca.vertcat(x, y), "f", f, "g", ca.vertcat(x + 1, x - 2)),
          d("error_on_fail", true));
      assertInException("process", () => solver2.call(@in("x0", 0, "lbg", 0, "ubg", 0)));
    }
  }

  static Function customJacF(MX x) {
    var J = ca.Function("jac_f", new[] { x, ca.MX(1, 1) }, new[] { -x }, new[] { "x", "out_r" }, new[] { "jac_r_x" });
    return ca.Function("f", new[] { x }, new[] { sq(x) }, new[] { "x" }, new[] { "r" }, d("custom_jacobian", J, "jac_penalty", 0));
  }

  public void test_postpone_expand() {
    Check.Assume(hasNlpsol("ipopt"), "no nlpsol ipopt");
    MX x = MX.sym("x");
    Function f = customJacF(x);
    Function solver = ca.nlpsol("solver", "ipopt", mxnlp("x", x, "f", f.Call(x - 3)), d("ipopt.max_iter", 5));
    solver.call(@in());
    Check.Equal("SOLVER_RET_LIMITED", (string)stat(solver, "unified_return_status"));
    Check.True(solver.get_function("nlp_jac_g").is_a("MXFunction"));

    f = customJacF(x);
    solver = ca.nlpsol("solver", "ipopt", mxnlp("x", x, "f", f.Call(x - 3)), d("ipopt.max_iter", 5, "expand", true));
    solver.call(@in());
    Check.Equal("Solve_Succeeded", (string)stat(solver, "return_status"));
    Check.True(solver.get_function("nlp_jac_g").is_a("SXFunction"));

    f = customJacF(x);
    solver = ca.nlpsol("solver", "ipopt", mxnlp("x", x, "f", f.Call(x - 3)), d("ipopt.max_iter", 5, "postpone_expand", true));
    solver.call(@in());
    Check.Equal("SOLVER_RET_LIMITED", (string)stat(solver, "unified_return_status"));
    Check.True(solver.get_function("nlp_jac_g").is_a("MXFunction"));

    f = customJacF(x);
    solver = ca.nlpsol("solver", "ipopt", mxnlp("x", x, "f", f.Call(x - 3)),
        d("ipopt.max_iter", 5, "expand", true, "postpone_expand", true));
    solver.call(@in());
    Check.Equal("SOLVER_RET_LIMITED", (string)stat(solver, "unified_return_status"));
    Check.True(solver.get_function("nlp_jac_g").is_a("SXFunction"));
  }

  class IterCallback : Callback {
    public readonly List<List<DM>> foo = new List<List<DM>>();
    readonly long nx, ng, np;
    readonly bool wrong;

    public IterCallback(long nx, long ng, long np, bool wrong) {
      this.nx = nx;
      this.ng = ng;
      this.np = np;
      this.wrong = wrong;
      construct("mycallback", d());
    }

    public override long get_n_in() => ca.nlpsol_n_out();
    public override long get_n_out() => 1;

    public override Sparsity get_sparsity_in(long i) {
      string n = ca.nlpsol_out(i);
      if (n == "f") return Sparsity.scalar();
      if (n == "x" || n == "lam_x") return Sparsity.dense(nx);
      if (n == "g" || n == "lam_g") return Sparsity.dense(wrong ? nx : ng);
      return ca.Sparsity(0, 0);
    }

    public override List<DM> eval(IEnumerable<DM> arg) {
      foo.Add(arg.ToList());
      return new List<DM> { ca.DM(0.0) };
    }
  }

  public void test_iteration_Callback() {
    Check.Assume(hasNlpsol("ipopt"), "no nlpsol ipopt");
    SX x = SX.sym("x"), y = SX.sym("y");
    SX f = rosen(1, x, y);
    var nlp = sxnlp("x", ca.vertcat(x, y), "f", f, "g", x + y);
    ca.Function("f", new[] { x, y }, new[] { f });

    var mycallback = new IterCallback(2, 1, 0, false);
    var opts = d("iteration_callback", mycallback, "ipopt.tol", 1e-8, "ipopt.max_iter", 50);
    Function solver = ca.nlpsol("solver", "ipopt", nlp, opts);
    solver.call(@in("lbx", -10, "ubx", 10, "lbg", -10, "ubg", 10));
    Check.Equal(iterCount(solver) + 1, (long)mycallback.foo.Count);

    var mycallback2 = new IterCallback(2, 1, 0, true);
    var opts2 = d("iteration_callback", mycallback2, "ipopt.tol", 1e-8, "ipopt.max_iter", 50);
    try {
      ca.nlpsol("solver", "ipopt", nlp, opts2);
    } catch (SwigException e) {
      Check.True(e.Message.Contains("Callback function input size mismatch"), e.Message);
    }
    keep.Add(mycallback);
    keep.Add(mycallback2);
  }

  public void test_pathological() {
    SX x = SX.sym("x"), y = SX.sym("y");
    var nlp = sxnlp("x", ca.vertcat(x, y), "f", sq(1 - x) + sq(y));
    foreach (S s in solvers) {
      if (s.name == "snopt" || s.name == "worhp" || s.name == "stabilizedsqp") continue;
      Function solver = ca.nlpsol("mysolver", s.name, nlp, s.opts);
      var solver_in = @in("x0", vc(1, 1), "lbx", vc(-10, -1), "ubx", vc(10, 2));
      var solver_out = solver.call(solver_in);
      chk(solver_out["f"], vc(0), 7);
      chk(solver_out["x"], vc(1, 0), 7, s.name);
      if (!s.name.Contains("bonmin")) chk(solver_out["lam_x"], vc(0, -0.0), 7, s.name);
      if (s.codegen != null) checkCodegen(solver, solver_in, s.codegen);
    }
  }

  public void test_pathological2() {
    SX x = SX.sym("x"), y = SX.sym("y");
    var nlp = sxnlp("x", ca.vertcat(x, y), "f", sq(1 - x) + y);
    foreach (S s in solvers) {
      if (s.name == "snopt") continue;
      Function solver = ca.nlpsol("mysolver", s.name, nlp, s.opts);
      var solver_in = @in("x0", vc(1, 1), "lbx", vc(-10, 0), "ubx", vc(10, 2));
      var solver_out = solver.call(solver_in);
      chk(solver_out["f"], vc(0), 7);
      chk(solver_out["x"], vc(1, 0), 7);
      if (!s.name.Contains("bonmin")) chk(solver_out["lam_x"], vc(0, -1), 7);
      if (s.codegen != null) checkCodegen(solver, solver_in, s.codegen);
    }
  }

  public void test_pathological3() {
    SX x = SX.sym("x"), y = SX.sym("y");
    var nlp = sxnlp("x", ca.vertcat(x, y), "f", sq(1 - x), "g", x + y);
    foreach (S s in solvers) {
      if (s.name == "worhp") continue;
      Function solver = ca.nlpsol("mysolver", s.name, nlp, s.opts);
      var solver_in = @in("x0", vc(1, 1), "lbx", vc(-10, 0), "ubx", vc(10, 2), "lbg", vc(2), "ubg", vc(2));
      var solver_out = solver.call(solver_in);
      chk(solver_out["f"], vc(0), 7);
      chk(solver_out["x"], vc(1, 1), 7);
      if (!s.name.Contains("bonmin")) chk(solver_out["lam_x"], vc(0, 0), 7);
      if (s.codegen != null) checkCodegen(solver, solver_in, s.codegen);
    }
  }

  public void test_pathological4() {
    SX x = SX.sym("x");
    var nlp = sxnlp("x", x, "f", x * x);
    foreach (S s in solvers) {
      if (new[] { "snopt", "worhp", "madnlp" }.Contains(s.name)) continue;
      Function solver = ca.nlpsol("mysolver", s.name, nlp, s.opts);
      var solver_in = @in("x0", vc(0), "lbx", vc(0), "ubx", vc(0));
      var solver_out = solver.call(solver_in);
      chk(solver_out["f"], vc(0), 7);
      chk(solver_out["x"], vc(0), 7);
      if (!s.name.Contains("bonmin")) chk(solver_out["lam_x"], vc(0), 7);
      if (s.codegen != null) checkCodegen(solver, solver_in, s.codegen);
    }
  }

  public void test_nlp_sensitivity() {
    MX x = MX.sym("x"), p = MX.sym("p");
    var nlp = mxnlp("x", x, "p", p, "f", sq(ca.sin(x) - sq(p)), "g", x);
    foreach (S s in solvers) {
      if (s.name.Contains("madnlp") || s.name.Contains("fatrop")) continue;
      if (s.@in("ipopt") || s.@in("snopt")) continue;
      if (s.name == "alpaqa") continue;
      var solver_options = s.opts;
      if (s.name.Contains("knitro")) {
        solver_options = new Dictionary<string, object>(s.opts);
        solver_options["knitro"] = with((Dictionary<string, object>)s.opts["knitro"], "algorithm", 4);
      }
      Function solver = ca.nlpsol("mysolver", s.name, nlp, solver_options);
      solver.call(@in("lbg", 0, "x0", 0, "p", 0.5));
      MX z = solver.call(min("p", p, "x0", x, "lbg", 0))["x"];
      MX z2 = ca.asin(sq(p));
      var f = ca.Function("f", new[] { x, p }, new[] { z, ca.jacobian(z, p) });
      var f2 = ca.Function("f", new[] { x, p }, new[] { z2, ca.jacobian(z2, p) });
      checkfunctionLight(f, f2, dms(0.0, 0.5), 6);
    }
  }

  public void test_nlp_sensitivity_equality_zero_multiplier() {
    Check.Assume(hasConic("qrqp"), "no conic qrqp");
    MX x = MX.sym("x", 3), p = MX.sym("p");
    var nlp = mxnlp("x", x, "p", p, "f", sq(x[0] - p), "g", x[1] - x[0]);
    Function solver = ca.nlpsol("solver", "sqpmethod", nlp, d("qpsol", "qrqp", "qpsol_options",
        d("print_iter", false, "print_header", false), "print_iteration", false, "print_header", false,
        "print_status", false, "sens_force_equality", true));
    var r = solver.call(min("p", p, "lbg", 0, "ubg", 0, "lbx", ca.vertcat(-inf, -inf, p), "ubx", ca.vertcat(inf, inf, p)));
    DM lam = ca.Function("lam", new[] { p }, new[] { ca.vertcat(r["lam_g"], r["lam_x"]) }).Call(ca.DM(0.0));
    Check.Equal(0.0, (double)ca.norm_inf(lam));
    MX z = r["x"];
    MX z2 = ca.vertcat(p, p, p);
    var f = ca.Function("f", new[] { p }, new[] { z, ca.jacobian(z, p) });
    var f2 = ca.Function("f", new[] { p }, new[] { z2, ca.jacobian(z2, p) });
    checkfunctionLight(f, f2, dms(0.0), 8);
  }

  public void test_regularize_sqpmethod() {
    Check.Assume(hasConic("qrqp"), "no conic qrqp");
    MX x = MX.sym("x", 2);
    MX f = 0.5 * ca.bilin(new double[,] { { 1, 0 }, { 0, -2 } }, x, x);
    var nlp = mxnlp("x", x, "f", f, "g", x[1]);

    Function solver = ca.nlpsol("mysolver", "sqpmethod", nlp, d("qpsol", "qrqp", "qpsol_options", d("print_problem", true)));
    Dictionary<string, DM> res = null;
    string result = captureStdout(() => res = solver.call(@in("lbg", 2, "ubg", 2)));
    Check.Equal(1L, iterCount(solver));
    Check.True(result.Contains("H:\n[[1, 0], \n [0, -2]]"), result);
    chk(res["x"], vc(0, 2), 6);

    Function solver2 = ca.nlpsol("mysolver", "sqpmethod", nlp, d("qpsol", "qrqp", "qpsol_options", d("print_problem", true),
        "convexify_strategy", "regularize", "convexify_margin", 0));
    result = captureStdout(() => res = solver2.call(@in("lbg", 2, "ubg", 2)));
    chk(res["x"], vc(0, 2), 6);
    Check.Equal(2L, iterCount(solver2));
    Check.True(result.Contains("H:\n[[3, 0], \n [0, 0]]"), result);

    Function solver3 = ca.nlpsol("mysolver", "sqpmethod", nlp, d("qpsol", "qrqp", "qpsol_options", d("print_problem", true),
        "convexify_strategy", "regularize", "convexify_margin", 1e-4));
    result = captureStdout(() => res = solver3.call(@in("lbg", 2, "ubg", 2)));
    chk(res["x"], vc(0, 2), 6);
    Check.Equal(2L, iterCount(solver3));
    Check.True(result.Contains("H:\n[[3.0001, 0], \n [0, 0.0001]]"), result);

    x = MX.sym("x", 2);
    f = 0.5 * ca.bilin(new double[,] { { 1, 0 }, { 0, 2 } }, x, x);
    var nlp2 = mxnlp("x", x, "f", f, "g", x[1]);

    Function solver4 = ca.nlpsol("mysolver", "sqpmethod", nlp2, d("qpsol", "qrqp", "qpsol_options", d("print_problem", true)));
    result = captureStdout(() => res = solver4.call(@in("lbg", 2, "ubg", 2)));
    Check.Equal(1L, iterCount(solver4));
    Check.True(result.Contains("H:\n[[1, 0], \n [0, 2]]"), result);
    chk(res["x"], vc(0, 2), 6);

    Function solver5 = ca.nlpsol("mysolver", "sqpmethod", nlp2, d("qpsol", "qrqp", "qpsol_options", d("print_problem", true),
        "convexify_strategy", "regularize"));
    result = captureStdout(() => res = solver5.call(@in("lbg", 2, "ubg", 2)));
    chk(res["x"], vc(0, 2), 6);
    Check.Equal(1L, iterCount(solver5));
    Check.True(result.Contains("H:\n[[1, 0], \n [0, 2]]"), result);
  }

  public void test_infeasible() {
    MX x = MX.sym("x");
    var nlp = mxnlp("f", sq(x), "x", x, "g", ca.sin(x));
    foreach (S s in solvers) {
      var myoptions = with(s.opts, "error_on_fail", true);
      if (s.name == "sqpmethod") continue;
      Function solver = ca.nlpsol("solver", s.name, nlp, myoptions);
      if (s.codegen != null && RunSlow) {
        Function F = checkCodegen(solver, @in("lbg", -10, "ubg", 10), s.codegen);
        Check.Throws<Exception>(() => F.call(@in("lbg", 5, "ubg", 10)));
        checkCodegen(solver, @in("lbg", 5, "ubg", 10), s.codegen, 15, false, true,
            new long[] { ca.SOLVER_RET_UNKNOWN, ca.SOLVER_RET_INFEASIBLE, ca.SOLVER_RET_EXCEPTION });
      }
    }
  }

  public void test_indefinite() {
    MX x = MX.sym("x", 2);
    MX f = 0.5 * ca.bilin(new double[,] { { 1, 0 }, { 0, -2 } }, x, x);
    var nlp = mxnlp("x", x, "f", f, "g", x[1]);
    foreach (S s in solvers) {
      var solver_in = @in("lbg", 2, "ubg", 2);
      Function solver = ca.nlpsol("mysolver", s.name, nlp, s.opts);
      var o = solver.call(solver_in);
      chk(o["x"], vc(0, 2), 6);
      if (!s.name.Contains("bonmin")) chk(o["lam_g"], vc(4), 6);
      if (s.codegen != null) {
        solver.generate("f.c", d("main", true));
        solver.generate_in("in.dat", solver.convert_in(solver_in));
        File.Delete("f.c");
        File.Delete("in.dat");
        checkCodegen(solver, solver_in, s.codegen);
      }
    }
  }

  public void test_gauss_newton_sqpmethod() {
    requireNlpsol("sqpmethod");
    SX x = SX.sym("x", 3);
    SX F = ca.sin(x) - vc(1, 2, 3) * ca.vertcat(x[2], 0, 0);
    SX J = ca.jacobian(F, x);
    SX f = 0.5 * ca.dot(F, F);
    SX p = SX.sym("x", 0, 1), lam_f = SX.sym("x"), lam_g = SX.sym("x", 0, 1);
    var GN = ca.Function("GN", new[] { x, p, lam_f, lam_g }, new[] { lam_f * ca.triu(ca.mtimes(J.T(), J)) });
    var nlp = sxnlp("x", x, "f", f);
    var options = d("hess_lag", GN);
    assertInException("Hessian must be symmetric", () => ca.nlpsol("solver", "sqpmethod", nlp, options));

    F = ca.sin(x - vc(0.22, 0.72, 0.2)) - vc(0.1, 0.5, 0.99);
    J = ca.jacobian(F, x);
    f = 0.5 * ca.dot(F, F);
    var Hf = ca.Function("H", new[] { x }, new[] { ca.hessian(f, x).Item1 });
    DM H = Hf.Call(ca.DM(0.0));
    var (eig, _) = eigh(H.ToArray());
    Check.True(eig.Any(e => e < 0));

    GN = ca.Function("GN", new[] { x, p, lam_f, lam_g }, new[] { lam_f * ca.mtimes(J.T(), J) });
    Function solver = ca.nlpsol("solver", "sqpmethod", sxnlp("x", x, "f", f), d("convexify_strategy", "regularize", "qpsol", "qrqp", "hess_lag", GN));
    solver.call(@in());
    Check.Equal(6L, iterCount(solver));

    solver = ca.nlpsol("solver", "sqpmethod", sxnlp("x", x, "f", f), d("convexify_strategy", "regularize", "qpsol", "qrqp"));
    solver.call(@in());
    Check.Equal(9L, iterCount(solver));
  }

  public void test_gauss_newton_ipopt() {
    Check.Assume(hasNlpsol("ipopt"), "no nlpsol ipopt");
    SX x = SX.sym("x", 3);
    SX F = ca.sin(x) - vc(1, 2, 3) * ca.vertcat(x[2], 0, 0);
    SX J = ca.jacobian(F, x);
    SX f = 0.5 * ca.dot(F, F);
    SX p = SX.sym("x", 0, 1), lam_f = SX.sym("x"), lam_g = SX.sym("x", 0, 1);
    var GN = ca.Function("GN", new[] { x, p, lam_f, lam_g }, new[] { lam_f * ca.mtimes(J.T(), J) });
    var nlp = sxnlp("x", x, "f", f);
    var options = d("hess_lag", GN);
    assertInException("Hessian must be upper triangular", () => ca.nlpsol("solver", "ipopt", nlp, options));

    F = ca.sin(x - vc(0.22, 0.72, 0.2)) - vc(0.1, 0.5, 0.99);
    J = ca.jacobian(F, x);
    f = 0.5 * ca.dot(F, F);
    var Hf = ca.Function("H", new[] { x }, new[] { ca.hessian(f, x).Item1 });
    DM H = Hf.Call(ca.DM(0.0));
    var (eig, _) = eigh(H.ToArray());
    Check.True(eig.Any(e => e < 0));

    GN = ca.Function("GN", new[] { x, p, lam_f, lam_g }, new[] { lam_f * ca.triu(ca.mtimes(J.T(), J)) });
    Function solver = ca.nlpsol("solver", "ipopt", sxnlp("x", x, "f", f), d("hess_lag", GN));
    solver.call(@in());
    Check.Equal(6L, iterCount(solver));

    solver = ca.nlpsol("solver", "ipopt", sxnlp("x", x, "f", f));
    solver.call(@in());
    Check.Equal(9L, iterCount(solver));
  }

  static DM blk(DM a, DM b) => ca.diagcat(a, b);

  public void test_cvx_sqpmethod() {
    Check.Assume(hasConic("qrqp"), "no conic qrqp");
    double eps = 1e-2;
    DM H0 = blk(new double[,] { { 2, 3 }, { 3, 0 } }, new double[,] { { 7, 1, 0 }, { 1, 0, 6 }, { 0, 6, 3 } });
    DM p = vc(0, 3, 2, 1, 4);
    DM Hp = H0.get(p, p);
    DM d3 = new double[,] { { 2, 0, 0 }, { 0, -1, 0 }, { 0, 0, 8 } };
    DM d3r = new double[,] { { 2, 0, 0 }, { 0, 1, 0 }, { 0, 0, 8 } };
    DM d3c = new double[,] { { 2, 0, 0 }, { 0, eps, 0 }, { 0, 0, 8 } };
    DM d3h = new double[,] { { 2, 0, 0 }, { 0, eps / 2, 0 }, { 0, 0, 8 } };
    DM b53 = new double[,] { { 5, 3 }, { 3, 6 } }, b58 = new double[,] { { 5, 8 }, { 8, 6 } }, b55 = new double[,] { { 5, 3 }, { 3, 5 } };
    var reg = d("convexify_strategy", "regularize", "convexify_margin", eps);
    var refl = d("convexify_strategy", "eigen-reflect", "convexify_margin", eps);
    var cl = d("convexify_strategy", "eigen-clip", "convexify_margin", eps);
    var none = new string[0];
    var cases = new (DM, DM, Dictionary<string, object>, string[], string[])[] {
      (ca.sparsify(new double[,] { { 0, 0 }, { 0, 2 } }), new double[,] { { eps, 0 }, { 0, 2 + eps } }, reg, none, none),
      (blk(d3, b53), blk(d3r, b53), refl, new[] { "Identified 2 blocks with maximum size 3" }, none),
      (blk(d3, b53), blk(d3c, b53), cl, none, none),
      (d3h, d3c, refl, none, none),
      (blk(d3, b58), blk(d3r, reflect(b58, eps)), refl, none, none),
      (blk(d3, b58), blk(d3c, clip(b58, eps)), cl, none, none),
      (d3h, d3c, refl, none, none),
      (blk(ca.sparsify(d3), b53), blk(d3r, b53), refl, new[] { "Identified 4 blocks with maximum size 2" }, none),
      (blk(ca.sparsify(d3), b55), blk(d3r, b55), refl, new[] { "Identified 4 blocks with maximum size 2" }, none),
      (Hp, reflect(Hp, eps), refl, none, none),
      (Hp, clip(Hp, eps), cl, none, none),
    };
    foreach (var (H, Hcvx, opts, includes_init, excludes_init) in cases) {
      long n = H.size1();
      MX x = MX.sym("H", n);
      var nlp = mxnlp("x", x, "f", 0.5 * ca.bilin(H, x, x));
      var options = with(d("max_iter", 1, "qpsol", "qrqp", "verbose", true), opts.SelectMany(e => new object[] { e.Key, e.Value }).ToArray());
      Function solver = null;
      assertOutput(includes_init, excludes_init, () => solver = ca.nlpsol("solver", "sqpmethod", nlp, options));
      DM x0 = DM.ones(n, 1);
      var res = solver.call(@in("x0", x0));
      chk(x0 - ca.solve(ca.densify(Hcvx), ca.mtimes(H, x0)), res["x"]);
      checkSerialize(solver, @in("x0", x0));
      checkCodegen(solver, @in("x0", x0), d("std", "c99"), codegen_check_digits);
    }
  }

  public void test_simple_bounds_detect() {
    SX x = SX.sym("x", 5), p = SX.sym("p", 1);
    var gl = new (double, SX, double)[] {
      (1.1, x[0] * x[1], 2.0),
      (-inf, x[4], 2.0),
      (-10.0, x[0], 10.0),
      (-5.0, x[0], 2.0),
      (-4.0, x[0], 4.0),
      (1.1, x[4] * x[1], 2.0),
      (0.0, x[4], inf),
      (7.0, x[2], 7.0),
      (-4.0, x[2], 40.0),
      (9.0, x[1], 9.0),
      (-4.0, x[0], 4.0),
      (-4.0, x[1], 9.0) };
    double[] lbg0 = gl.Select(t => t.Item1).ToArray(), ubg0 = gl.Select(t => t.Item3).ToArray();
    SX g = ca.vertcat(gl.Select(t => t.Item2));
    var (_, lbxs, ubxs, lam_f, lam_b) = ca.detect_simple_bounds(x, p, g, (DM)lbg0, (DM)ubg0);
    DM lbx = ca.DM(lbxs), ubx = ca.DM(ubxs);
    chk(lbx.T(), row(-4, 9, 7, -inf, 0));
    chk(ubx.T(), row(2, 9, 7, inf, 2));

    Func<DM, DM> round_trip_f = a => {
      List<DM> o = lam_f.call(dms(a, 0.0));
      return lam_b.call(dms(o[0], o[1], 0.0))[0];
    };
    Func<List<DM>, List<DM>> round_trip_b = a => {
      DM o = lam_b.call(dms(a[0], a[1], 0.0))[0];
      return lam_f.call(dms(o, 0.0));
    };

    double[,] G = {
      { 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
      { 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
      { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
      { 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0 },
      { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
      { 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0 },
      { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
      { 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0 },
      { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
      { 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1 },
      { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
      { 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1 } };
    for (int i = 0; i < 12; i++) {
      DM a = DM.zeros(12, 1);
      a[i] = 2.0;
      DM b = round_trip_f(a), c = round_trip_f(round_trip_f(a));
      chk(b, vc(rowOf(G, i)));
      chk(b, c);
    }

    G = new double[,] {
      { -2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
      { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
      { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
      { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
      { 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, -1, 0 },
      { 0, 0, 0, 0, 0, -2, 0, 0, 0, 0, 0, 0 },
      { 0, 0, 0, 0, 0, 0, -2, 0, 0, 0, 0, 0 },
      { 0, 0, 0, 0, 0, 0, 0, -2, 0, 0, 0, 0 },
      { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
      { 0, 0, 0, 0, 0, 0, 0, 0, 0, -2, 0, 0 },
      { 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, -1, 0 },
      { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 } };
    for (int i = 0; i < 12; i++) {
      DM a = DM.zeros(12, 1);
      a[i] = -2.0;
      DM b = round_trip_f(a), c = round_trip_f(round_trip_f(a));
      chk(b, vc(rowOf(G, i)));
      chk(b, c);
    }

    foreach (double sgn in new double[] { -2, 2 }) {
      var G5 = new double[5, 5];
      foreach (int i in new[] { 0, 1, 2, 4 }) G5[i, i] = sgn;
      for (int i = 0; i < 5; i++) {
        DM a = DM.zeros(5, 1);
        a[i] = sgn;
        List<DM> b = round_trip_b(dms(vc(3, 7), a));
        List<DM> c = round_trip_b(round_trip_b(dms(vc(3, 7), a)));
        chk(b[0].T(), row(3, 7));
        chk(b[1], vc(rowOf(G5, i)));
        chk(b[1], c[1]);
      }
    }
  }

  static double[] rowOf(double[,] m, int i) => Enumerable.Range(0, m.GetLength(1)).Select(j => m[i, j]).ToArray();

  static (double, MX, double)[] glMX(MX x) => new (double, MX, double)[] {
    (1.1, x[0] * x[1], 2.0),
    (-11.0, x[4], 2.0),
    (-10.0, x[0], 10.0),
    (-5.0, x[0], 2.0),
    (-4.0, x[0], 4.0),
    (1.1, x[4] * x[1], 2.0),
    (0.0, x[4], inf),
    (3.0, x[2], 3.0),
    (-4.0, x[2], 40.0),
    (2.0, x[1], 2.0),
    (-4.0, x[0], 4.0),
    (-4.0, x[1], 9.0) };

  public void test_simple_bounds_detect_solvers() {
    Check.Assume(hasNlpsol("ipopt"), "no nlpsol ipopt");
    foreach (S s in solvers) {
      if (s.name == "sleqp" || s.name == "madnlp") continue;
      MX x = MX.sym("x", 5);
      var gl = glMX(x);
      DM lbg = vc(gl.Select(t => t.Item1).ToArray()), ubg = vc(gl.Select(t => t.Item3).ToArray());
      MX g = ca.vertcat(gl.Select(t => t.Item2));
      DM I = ca.densify(DM.eye(5));
      DM lbx = -DM.inf(5), ubx = DM.inf(5);
      for (int i = 0; i < 3; i++) {
        foreach (double sign in new double[] { -1, 1 }) {
          DM x_target = I[.., i] * (sign * 20);
          var nlp = mxnlp("x", x, "g", g, "f", ca.sumsqr(x - x_target));
          Function solver_ref_ipopt = ca.nlpsol("mysolver", "ipopt", nlp);
          Function solver_ref = ca.nlpsol("mysolver", s.name, nlp, s.opts);
          Function solver = ca.nlpsol("mysolver", s.name, nlp, with(s.opts, "detect_simple_bounds", true));
          var solver_in = @in("lbg", lbg, "ubg", ubg, "lbx", lbx, "ubx", ubx);
          var sol_ref_ipopt = solver_ref_ipopt.call(solver_in);
          solver_in["x0"] = sol_ref_ipopt["x"] * 1.2;
          int digits = 6;
          if (s.@in("daqp")) continue;
          if (s.@in("worhp")) digits = 4;
          solver_ref.call(solver_in);
          if (s.name == "bonmin") {
            Dictionary<string, DM> solver_ref_out = solver_ref.call(solver_in), solver_out = solver.call(solver_in);
            foreach (string output in new[] { "x", "f" })
              chk(solver_out[output], solver_ref_out[output], digits, s.ToString());
          } else {
            checkfunctionLight(solver, solver_ref, positional(solver, solver_in), digits, s.ToString());
          }
          if (s.codegen != null) checkCodegen(solver, solver_in, s.codegen);
        }
      }
    }
  }

  static void signIssueCase(Function solver_ref, Function solver, double lbg, double ubg) =>
    checkfunctionLight(solver, solver_ref, solver.convert_in(@in("lbg", lbg, "ubg", ubg)), 9);

  public void test_simple_bounds_sign_issue() {
    Check.Assume(hasNlpsol("ipopt"), "no nlpsol ipopt");
    double[,] b = { { -2, 3 }, { 2, 2 }, { -2, 3 }, { 2, 2 }, { -2, 3 }, { -2, 3 }, { -2, 3 } };
    double[] c = { 1, 1, 0.3, 0.7, -1, -0.7, 0.3 };
    SX xs = SX.sym("x");
    foreach (SX f in new[] { xs, -xs })
      for (int k = 0; k < c.Length; k++) {
        var nlp = sxnlp("x", xs, "f", f, "g", c[k] * xs);
        signIssueCase(ca.nlpsol("solver", "ipopt", nlp), ca.nlpsol("solver", "ipopt", nlp, d("detect_simple_bounds", true)), b[k, 0], b[k, 1]);
      }
    MX xm = MX.sym("x");
    foreach (MX f in new[] { xm, -xm })
      for (int k = 0; k < c.Length; k++) {
        var nlp = mxnlp("x", xm, "f", f, "g", c[k] * xm);
        signIssueCase(ca.nlpsol("solver", "ipopt", nlp), ca.nlpsol("solver", "ipopt", nlp, d("detect_simple_bounds", true)), b[k, 0], b[k, 1]);
      }
  }

  public void test_simple_bounds_equality() {
    Check.Assume(hasNlpsol("ipopt"), "no nlpsol ipopt");
    MX x = MX.sym("x"), y = MX.sym("y");
    MX z = ca.vertcat(x, y);
    MX f = x;
    MX g = ca.vertcat(sq(x) + sq(y), x, x * y);
    DM lbg = vc(1, -2, -5), ubg = vc(1, 2, 5), x0 = vc(0.2, 0.6);
    Function solver_ref = ca.nlpsol("solver", "ipopt", mxnlp("x", z, "f", f, "g", g), d("equality", new[] { true, false, false }));
    Function solver = ca.nlpsol("solver", "ipopt", mxnlp("x", z, "f", f, "g", g),
        d("detect_simple_bounds", true, "equality", new[] { true, false, false }));
    checkfunctionLight(solver, solver_ref, solver.convert_in(@in("x0", x0, "lbg", lbg, "ubg", ubg)), 9);
  }

  public void test_issue_3407() {
    Check.Assume(hasNlpsol("ipopt"), "no nlpsol ipopt");
    var opti = ca.Opti();
    MX x = opti.variable();
    opti.minimize(sq(x));
    var f = ca.Function("f", new[] { x }, new[] { 3 * x, sq(x.attachAssert(x > 1)) }, d("never_inline", true));
    List<MX> yz = f.call(new[] { x });
    opti.subject_to(yz[0] <= 10);
    opti.subject_to(ca.eq(yz[1], 3));
    opti.set_initial(x, 3);
    opti.solver("ipopt", d("detect_simple_bounds", true));
    opti.solve();
  }

  public void test_issue_4235() {
    Check.Assume(hasNlpsol("ipopt"), "no nlpsol ipopt");
    foreach (bool expand in new[] { false, true }) {
      var opti = ca.Opti();
      MX x = opti.variable();
      opti.minimize(sq(x));
      var f = ca.Function("f", new[] { x }, new[] { 3 * x, sq(x.attachAssert(x > 1000)) });
      List<MX> yz = f.call(new[] { 2 * x });
      opti.subject_to(yz[0] <= 10);
      opti.subject_to(ca.eq(yz[1], 3));
      opti.set_initial(x, 3);
      if (!expand) {
        opti.solver("ipopt", d("detect_simple_bounds", true));
        Check.Throws<Exception>(() => opti.solve());
      } else {
        // Expand will have the side effect of removing the assertion
        opti.solver("ipopt", d("expand", true, "detect_simple_bounds", true));
        opti.solve();
      }
    }
  }

  public void test_simple_bounds_detect2() {
    Check.Assume(hasNlpsol("ipopt"), "no nlpsol ipopt");
    MX x = MX.sym("x", 5), p = MX.sym("p", 5);
    var gl = glMX(x);
    DM lbg = vc(gl.Select(t => t.Item1).ToArray()), ubg = vc(gl.Select(t => t.Item3).ToArray());
    MX g = ca.vertcat(gl.Select(t => t.Item2));
    var nlp = mxnlp("x", x, "p", p, "g", g, "f", ca.sumsqr(x - p));
    int print_level = 0;

    assertOutput(new[] { "Total number of inequality constraints...............:        2" }, new string[0],
        () => ca.nlpsol("solver", "ipopt", nlp, d("detect_simple_bounds", true)).call(@in("x0", 0, "lbg", lbg, "ubg", ubg)));
    assertOutput(new[] { "Total number of inequality constraints...............:       10" }, new string[0],
        () => ca.nlpsol("solver", "ipopt", nlp, d("detect_simple_bounds", false)).call(@in("x0", 0, "lbg", lbg, "ubg", ubg)));

    var options = d("detect_simple_bounds", true, "ipopt.print_level", print_level, "print_time", false,
        "ipopt.tol", 1e-12, "ipopt.fixed_variable_treatment", "relax_bounds", "ipopt.bound_relax_factor", 1e-12);
    Function solver = ca.nlpsol("solver", "ipopt", nlp, options);
    var options_nominal = d("detect_simple_bounds", false, "ipopt.print_level", print_level,
        "print_time", false, "ipopt.tol", 1e-12, "ipopt.bound_relax_factor", 1e-12);
    Function solver_nominal = ca.nlpsol("solver", "ipopt", nlp, options_nominal);
    object[] w_options = { "ipopt.warm_start_init_point", "yes", "ipopt.warm_start_bound_push", 1e-16,
        "ipopt.bound_push", 1e-16, "ipopt.warm_start_mult_bound_push", 1e-16, "ipopt.tol", 1e-7 };
    Function wsolver = ca.nlpsol("solver", "ipopt", nlp, with(options, w_options));
    Function wsolver_nominal = ca.nlpsol("solver", "ipopt", nlp, with(options_nominal, w_options));

    DM I = ca.densify(DM.eye(5));
    var bounds = new (DM, DM)[] {
      (-DM.inf(5), DM.inf(5)),
      (-3 * DM.ones(5), 3 * DM.ones(5)),
      (-4 * DM.ones(5), 4 * DM.ones(5)) };
    foreach (var (lbx, ubx) in bounds) {
      for (int i = 0; i < 5; i++) {
        foreach (double sign in new double[] { -1, 1 }) {
          DM x_target = I[.., i] * (sign * 20);
          var o = solver.call(@in("p", x_target, "lbg", lbg, "ubg", ubg, "lbx", lbx, "ubx", ubx));
          var out_nominal = solver_nominal.call(@in("p", x_target, "lbg", lbg, "ubg", ubg, "lbx", lbx, "ubx", ubx));
          if (el(lbx, 0) == -inf) chk(o["lam_x"], vc(rep(0, 5)));
          chk(o["x"], out_nominal["x"], 8);
          chk(o["g"], out_nominal["g"], 8);

          wsolver.call(@in("x0", o["x"], "lam_x0", o["lam_x"], "lam_g0", o["lam_g"], "p", x_target,
              "lbg", lbg, "ubg", ubg, "lbx", lbx, "ubx", ubx));
          Check.Equal(0L, iterCount(wsolver));
          wsolver_nominal.call(@in("x0", o["x"], "lam_x0", o["lam_x"], "lam_g0", o["lam_g"], "p", x_target,
              "lbg", lbg, "ubg", ubg, "lbx", lbx, "ubx", ubx));
          Check.Equal(0L, iterCount(wsolver_nominal));
          wsolver.call(@in("x0", out_nominal["x"], "lam_x0", out_nominal["lam_x"], "lam_g0", out_nominal["lam_g"],
              "p", x_target, "lbg", lbg, "ubg", ubg, "lbx", lbx, "ubx", ubx));
          Check.Equal(0L, iterCount(wsolver));
        }
      }
    }
  }

  public void test_derivative() {
    MX x = MX.sym("x", 3), p = MX.sym("p", 3);
    var nlp = mxnlp("x", x, "p", p, "f", ca.sumsqr(ca.sin(x) - p), "g", x);
    Function solver = ca.nlpsol("solver", "sqpmethod", nlp, d("qpsol", "qrqp"));
    var res = solver.call(min("p", p, "x0", 0, "lbg", -inf, "ubg", inf));
    var f = ca.Function("f", new[] { p }, new[] { res["x"] });
    var f_ref = ca.Function("f_ref", new[] { p }, new[] { ca.asin(p) });
    checkfunction(f, f_ref, dms(vc(0.1, 0.2, 0.3)));

    x = MX.sym("x", 3);
    p = MX.sym("p", 3);
    nlp = mxnlp("x", x, "p", p, "f", ca.sumsqr(ca.sin(x) - p));
    solver = ca.nlpsol("solver", "sqpmethod", nlp, d("qpsol", "qrqp"));
    res = solver.call(min("p", p, "x0", 0));
    f = ca.Function("f", new[] { p }, new[] { res["x"] });
    f_ref = ca.Function("f_ref", new[] { p }, new[] { ca.asin(p) });
    checkfunction(f, f_ref, dms(vc(0.1, 0.2, 0.3)));
  }

  // test_hock_schittkowski: returns immediately in python (body needs pandas and the hock_schittkowski modules)
  public void test_hock_schittkowski() {}

  public void test_exception_in_oraclefunction() {
    MX x = MX.sym("x");
    MX x_fail = x.attachAssert(ca.ne(x, x), "Cuckoo");
    var nlp = mxnlp("x", x, "f", sq(x - 1), "g", x_fail);
    foreach (S s in solvers) {
      if (s.name.Contains("madnlp")) continue;
      Function solver = ca.nlpsol("mysolver", s.name, nlp, s.opts);
      var solver_in = @in("lbx", vc(-10), "ubx", vc(10), "lbg", vc(-10), "ubg", vc(10));
      assertInAnyOutput("Cuckoo", () => solver.call(solver_in));
    }
  }

  public void test_unsolved_stats() {
    MX x = MX.sym("x");
    MX x_fail = x.attachAssert(ca.ne(x, x), "Cuckoo");
    var nlp = mxnlp("x", x, "f", sq(x - 1), "g", x_fail);
    foreach (S s in solvers) {
      Function solver = ca.nlpsol("mysolver", s.name, nlp, s.opts);
      assertInException("No stats available", () => solver.stats());
    }
  }

  public void test_ipopt_detect_simple_bounds() {
    Check.Assume(hasNlpsol("ipopt"), "no nlpsol ipopt");
    MX x = MX.sym("x"), y = MX.sym("y"), z = MX.sym("z"), w = MX.sym("w");
    var testcases = new List<(MX, DM, DM, DM, DM, double[])>();
    foreach (double d_scale in new[] { 1.0, 1e5 }) {
      foreach (double c_scale in new[] { 1.0, 1e5 }) {
        MX g = ca.vertcat(x, (x + y) * d_scale, sq(x) + sq(z), z, (z + y + w) * c_scale, sq(z) + sq(y) + w);
        DM lbg = vc(-5.0, -2 * d_scale, -3.0, 1.0, 2 * c_scale, 3.0), ubg = vc(3.0, 2 * d_scale, 3.0, 1.0, 2 * c_scale, 3.0);
        testcases.Add((g, lbg, ubg, vc(-inf, -inf, -inf, -inf), vc(inf, inf, inf, inf), new double[] { 0, 3, 5, 1, 4, 2 }));
        g = ca.vertcat((x + y) * d_scale, sq(x) + sq(z), (z + y + w) * c_scale, sq(z) + sq(y) + w);
        lbg = vc(-2 * d_scale, -3.0, 2 * c_scale, 3.0);
        ubg = vc(2 * d_scale, 3.0, 2 * c_scale, 3.0);
        testcases.Add((g, lbg, ubg, vc(-5, -inf, 1, -inf), vc(3, inf, 1, inf), new double[] { 3, 0, 2, 1 }));
      }
    }
    MX X = ca.vertcat(x, y, z, w);
    MX F = sq(x) + 2 * sq(y) + 3 * sq(z) + 8 * sq(w);
    foreach (var (g, lbg, ubg, lbx, ubx, perm) in testcases) {
      Function solver_ref = ca.nlpsol("solver", "ipopt", mxnlp("x", X, "g", g, "f", F));
      var refv = solver_ref.call(@in("lbg", lbg, "ubg", ubg, "lbx", lbx, "ubx", ubx));
      foreach (long[] g_perm in new[] { Enumerable.Range(0, (int)g.numel()).Select(i => (long)i).ToArray(), perm.Select(v => (long)v).ToArray() }) {
        var nlp = mxnlp("x", X, "g", g[g_perm], "f", F);
        Function solver_ref2 = ca.nlpsol("solver", "ipopt", nlp);
        var ref2 = solver_ref2.call(@in("lbg", lbg[g_perm], "ubg", ubg[g_perm], "lbx", lbx, "ubx", ubx));
        var arg = @in("x0", refv["x"] * 1.1, "lbg", lbg[g_perm], "ubg", ubg[g_perm], "lbx", lbx, "ubx", ubx);
        foreach (bool detect_simple_bounds in new[] { true, false }) {
          foreach (string fixed_variable_treatment in new[] { "make_constraint", "make_parameter", "make_parameter_nodual" }) {
            foreach (string start_with_resto in new[] { "no", "yes" }) {
              string[] keys2 = fixed_variable_treatment == "make_constraint" ? new[] { "x", "f", "g", "lam_g" } : new string[0];
              var o = d("ipopt.fixed_variable_treatment", fixed_variable_treatment,
                  "detect_simple_bounds", detect_simple_bounds, "ipopt.start_with_resto", start_with_resto);
              Function solver = ca.nlpsol("solver", "ipopt", nlp, o);
              var res = solver.call(arg);
              foreach (string k in new[] { "x", "f" }) chk(res[k], refv[k], 7);
              foreach (string k in keys2) chk(res[k], ref2[k], 7);
              Check.True(success(solver));

              var cache = d("nlp_jac_g", solver.get_function("nlp_jac_g"), "nlp_hess_l", solver.get_function("nlp_hess_l"));
              solver = ca.nlpsol("solver", "ipopt", nlp, with(o, "cache", cache));
              res = solver.call(arg);
              foreach (string k in new[] { "f", "x" }) chk(res[k], refv[k], 7);
              foreach (string k in keys2) chk(res[k], ref2[k], 7);

              cache = d("nlp_jac_g", solver_ref2.get_function("nlp_jac_g"), "nlp_hess_l", solver_ref2.get_function("nlp_hess_l"));
              solver = ca.nlpsol("solver", "ipopt", nlp, with(o, "cache", cache));
              res = solver.call(arg);
              foreach (string k in new[] { "f", "x" }) chk(res[k], refv[k], 7);
              foreach (string k in keys2) chk(res[k], ref2[k], 7);
              Check.True(success(solver));
            }
          }
        }
      }
    }
  }

  public void test_ipopt_custom_hess() {
    Check.Assume(hasNlpsol("ipopt"), "no nlpsol ipopt");
    SX x = SX.sym("x"), y = SX.sym("y"), w = SX.sym("w");
    SX p = ca.SX(0, 1);
    SX f = sq(1 - x) + 100 * sq(w);
    SX g = y - sq(x) - w;
    SX X = ca.vertcat(x, y, w);
    var nlp = sxnlp("x", X, "f", f, "g", g);
    SX lam_f = SX.sym("lam_f"), lam_g = SX.sym("lam_g", g.numel());
    Function ref_solver = ca.nlpsol("solver", "ipopt", nlp);
    var ref_sol = ref_solver.call(@in("x0", 0, "lbg", 0, "ubg", 0));

    var nlp_hess_l_custom = ca.Function("nlp_hess_l", new[] { X, p, lam_f, lam_g }, new SX[] { ca.triu(DM.zeros(3, 3)) });
    Function solver = ca.nlpsol("solver", "ipopt", nlp, d("cache", d("nlp_hess_l", nlp_hess_l_custom)));
    List<DM> hin = dms(vc(0.11, 0.3, 0.7), 0.0, 1.2, 3.7);
    Function s1 = solver;
    // python's assertInException catches the checkarray AssertionError ("evaluation error")
    var err = Check.Throws<AssertFailed>(
        () => checkfunctionLight(ref_solver.get_function("nlp_hess_l"), s1.get_function("nlp_hess_l"), hin, 9));
    Check.True(err.Message.Contains("evaluation error"), err.Message);

    SX lag = lam_f * f + ca.dot(lam_g, g);
    SX H = ca.jacobian(ca.gradient(lag, X), X, d("symmetric", true));
    nlp_hess_l_custom = ca.Function("nlp_hess_l", new[] { X, p, lam_f, lam_g }, new[] { ca.triu(H) });
    ref_solver.get_function("nlp_hess_l").disp(true);
    nlp_hess_l_custom.disp(true);
    solver = ca.nlpsol("solver", "ipopt", nlp, d("cache", d("nlp_hess_l", nlp_hess_l_custom)));
    var sol = solver.call(@in("x0", 0, "lbg", 0, "ubg", 0));
    Check.True(success(solver));
    chk(sol["x"], ref_sol["x"], 6);
    checkfunctionLight(ref_solver.get_function("nlp_hess_l"), solver.get_function("nlp_hess_l"), hin, 9);
  }

  public void test_ipopt_custom_jac() {
    Check.Assume(hasNlpsol("ipopt"), "no nlpsol ipopt");
    SX x = SX.sym("x"), y = SX.sym("y"), w = SX.sym("w");
    SX p = ca.SX(0, 1);
    SX f = sq(1 - x) + 100 * sq(w);
    SX g = y - sq(x) - w;
    SX X = ca.vertcat(x, y, w);
    var nlp = sxnlp("x", X, "f", f, "g", g);
    Function ref_solver = ca.nlpsol("solver", "ipopt", nlp);
    ref_solver.call(@in("x0", 0, "lbg", 0, "ubg", 0));
    var options = d();
    Function solver = ca.nlpsol("solver", "ipopt", nlp, options);
    var nlp_jac_g_custom = ca.Function("nlp_jac_g", new[] { X, p }, new[] { g, ca.jacobian(g, X) }, new[] { "x", "p" }, new[] { "g", "jac_g_x" });
    options["nlp_jac_g"] = d("nlp_jac_g", nlp_jac_g_custom);
    solver.call(@in("x0", 0, "lbg", 0, "ubg", 0));
    Check.True(success(solver));
    checkfunctionLight(ref_solver.get_function("nlp_jac_g"), solver.get_function("nlp_jac_g"), dms(vc(0.11, 0.3, 0.7), 0.0), 9);
  }

  public void test_nlp_signature() {
    MX x = MX.sym("x");
    assertInException("must have inputs (x, p) and outputs (f, g)",
        () => ca.nlpsol("solver", "ipopt", ca.Function("nlp", new[] { x }, new[] { sq(x) }, new[] { "a" }, new[] { "b" })));
  }

  public void test_option_propagation() {
    Check.Assume(hasNlpsol("ipopt"), "no nlpsol ipopt");
    MX x = MX.sym("x");
    MX y = sq(x.printme(0));
    var options = d("common_options", d("helper_options", d("enable_fd", true, "enable_forward", false,
        "enable_reverse", false)), "ipopt", d("resto.tol", 1e-8, "hessian_approximation", "limited-memory", "max_iter", 0));
    Function solver = ca.nlpsol("solver", "ipopt", mxnlp("x", x, "f", y), options);
    string[] result = captureBoth(() => solver.call(@in("x0", 1)));
    Check.Equal(0, result[1].Length, result[1]);
    var r = new List<double>();
    foreach (string l in result[0].Split('\n'))
      if (l.Contains("|> 0 : ")) r.Add(double.Parse(l.Split(':')[1].Trim(), System.Globalization.CultureInfo.InvariantCulture));
    double spread = r.Count == 0 ? 0 : r.Max() - r.Min();
    Check.True(spread > 0, string.Join(", ", r));
    Check.True(spread < 1e-4, string.Join(", ", r));
  }

  public void test_thread_map_reuses_solver_memory() {
    var cases = new (string, Dictionary<string, object>)[] {
      ("ipopt", d()),
      ("ipopt", d("ipopt.print_level", 0, "ipopt.sb", "yes", "print_time", false)),
      ("bonmin", d()) };
    Check.Assume(hasNlpsol("ipopt") || hasNlpsol("bonmin"), "needs ipopt or bonmin");
    foreach (var (plugin, opts) in cases) {
      if (!hasNlpsol(plugin)) continue;
      MX x = MX.sym("x", 2), p = MX.sym("p");
      Function solver = ca.nlpsol("solver", plugin, mxnlp("x", x, "p", p, "f", ca.sumsqr(x - p)), opts);
      MX P = MX.sym("P");
      var par = ca.Function("par", new[] { P }, new[] { solver.call(min("x0", 0, "p", P))["x"][0] });
      Function m = par.map(8, "thread", 4);
      DM ps = ((DM)range(8)).T();
      for (int k = 0; k < 5; k++) chk(m.Call(ps), ps, 6);
    }
  }
}

}

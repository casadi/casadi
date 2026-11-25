using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;
using System.Runtime.InteropServices;
using System.Text.RegularExpressions;

namespace CasADi.Tests {

/// <summary>Port of test/java ConicTest (itself a port of test/python/conic.py)</summary>
public class ConicTest : CasadiTestCase {

  const double inf = double.PositiveInfinity;
  static readonly bool isWindows = RuntimeInformation.IsOSPlatform(OSPlatform.Windows);
  static readonly bool isDarwin = RuntimeInformation.IsOSPlatform(OSPlatform.OSX);

  sealed class Conic {
    public readonly string name;
    public readonly Dictionary<string, object> options;
    public readonly Dictionary<string, object> aux;

    public Conic(string name, Dictionary<string, object> options, Dictionary<string, object> aux) {
      this.name = name;
      this.options = options;
      this.aux = aux;
    }

    /// <summary>Fresh mutable copy of the plugin options</summary>
    public Dictionary<string, object> opts() => new Dictionary<string, object>(options);

    public bool has(string k) => aux.TryGetValue(k, out object v) && v is bool b && b;

    public int lessDigits() {
      if (!aux.TryGetValue("less_digits", out object v)) return 0;
      if (v is bool b) return b ? 1 : 0;
      return Convert.ToInt32(v);
    }

    public Dictionary<string, object> codegen() =>
      aux.TryGetValue("codegen", out object v) ? v as Dictionary<string, object> : null;

    public override string ToString() => name + " " + str(options);
  }

  static Dictionary<string, object> map(params object[] kv) {
    var r = new Dictionary<string, object>();
    for (int i = 0; i < kv.Length; i += 2) r[(string)kv[i]] = kv[i + 1];
    return r;
  }

  static string str(object o) {
    switch (o) {
      case string s: return s;
      case IDictionary<string, object> d: return "{" + string.Join(", ", d.Select(kv => kv.Key + "=" + str(kv.Value))) + "}";
      case System.Collections.IEnumerable e: return "[" + string.Join(", ", e.Cast<object>().Select(str)) + "]";
      default: return o?.ToString() ?? "null";
    }
  }

  static bool enabled(string skip, bool available) => Environment.GetEnvironmentVariable(skip) == null && available;

  static readonly List<Conic> conics = new List<Conic>();

  static ConicTest() {
    if (enabled("SKIP_IPOPT_TESTS", hasNlpsol("ipopt"))) {
      var ipopt_options = map("fixed_variable_treatment", "relax_bounds", "jac_c_constant", "yes", "jac_d_constant", "yes",
          "hessian_constant", "yes", "tol", 1e-12, "print_level", 0);
      conics.Add(new Conic("nlpsol", map("nlpsol", "ipopt", "nlpsol_options.ipopt", ipopt_options),
          map("quadratic", true, "dual", true, "soc", false, "codegen", false, "binary", false, "discrete", false, "sos", false)));
    }
    if (enabled("SKIP_WORHP_TESTS", hasNlpsol("worhp")))
      conics.Add(new Conic("nlpsol", map("nlpsol", "worhp", "nlpsol_options.worhp", map("TolOpti", 1e-13)),
          map("less_digits", 1, "quadratic", true, "dual", false, "soc", false, "codegen", false, "binary", false, "discrete", false, "sos", false)));
    if (enabled("SKIP_OOQP_TESTS", hasConic("ooqp")))
      conics.Add(new Conic("ooqp", map(),
          map("less_digits", 1, "quadratic", true, "dual", true, "soc", false, "codegen", false, "binary", false, "discrete", false, "sos", false)));
    if (enabled("SKIP_QPOASES_TESTS", hasConic("qpoases")))
      conics.Add(new Conic("qpoases", map("printLevel", "low"),
          map("quadratic", true, "dual", true, "soc", false, "codegen", false, "binary", false, "discrete", false, "sos", false)));
    if (enabled("SKIP_CPLEX_TESTS", hasConic("cplex")))
      conics.Add(new Conic("cplex", map("cplex", map("CPX_PARAM_BARQCPEPCOMP", 1e-10, "CPX_PARAM_BAREPCOMP", 1e-10)),
          map("quadratic", true, "dual", true, "soc", true, "codegen", false, "binary", true, "discrete", true, "sos", true)));
    if (enabled("SKIP_OSQP_TESTS", hasConic("osqp"))) {
      var options = isWindows ? new List<string>() : new List<string> { "-Wno-unused-variable", "-Wno-strict-prototypes" };
      var codegen = map("extralibs", new List<string> { "osqp" }, "extra_options", options, "std", "c99");
      if (isDarwin) codegen["digits"] = 13;
      conics.Add(new Conic("osqp", map("osqp", map("alpha", 1, "eps_abs", 1e-8, "eps_rel", 1e-8)),
          map("quadratic", true, "dual", true, "codegen", codegen, "soc", false, "binary", false, "discrete", false)));
    }
    if (enabled("SKIP_SUPERSCS_TESTS", hasConic("superscs")))
      conics.Add(new Conic("superscs", map("superscs", map("eps", 1e-9, "do_super_scs", 1, "verbose", 0)),
          map("quadratic", true, "dual", true, "codegen", false, "soc", true, "binary", false, "discrete", false)));
    if (enabled("SKIP_GUROBI_TESTS", hasConic("gurobi")))
      conics.Add(new Conic("gurobi", map("gurobi", map("BarQCPConvTol", 1e-9)),
          map("quadratic", true, "less_digits", true, "dual", true, "soc", true, "codegen", false, "binary", true, "discrete", true, "sos", true)));
    if (enabled("SKIP_XPRESS_TESTS", hasConic("xpress"))) {
      string xpress_dir = Environment.GetEnvironmentVariable("XPRESSDIR") ?? "missing";
      var xpress_codegen = map("extralibs", new List<string> { "xprs" }, "std", "c99",
          "extralibdirs", new List<string> { Path.Combine(xpress_dir, "lib") }, "extra_include", new List<string> { Path.Combine(xpress_dir, "include") });
      conics.Add(new Conic("xpress", map("xpress", map("OUTPUTLOG", 0, "FEASTOL", 1e-9, "BARGAPSTOP", 1e-10, "BARDUALSTOP", 1e-10)),
          map("quadratic", true, "dual", true, "soc", true, "codegen", xpress_codegen, "binary", true, "discrete", true, "sos", true)));
    }
    if (enabled("SKIP_MOSEK_TESTS", hasConic("mosek"))) {
      string mosek_dir = Environment.GetEnvironmentVariable("MOSEKDIR") ?? "missing";
      var mosek_codegen = map("extralibs", new List<string> { "mosek64" }, "std", "c99",
          "extralibdirs", new List<string> { Path.Combine(mosek_dir, "bin") }, "extra_include", new List<string> { Path.Combine(mosek_dir, "h") });
      conics.Add(new Conic("mosek", map("mosek", map("MSK_IPAR_LOG", 0, "MSK_DPAR_INTPNT_QO_TOL_REL_GAP", 1e-10, "MSK_DPAR_INTPNT_QO_TOL_PFEAS", 1e-10,
          "MSK_DPAR_INTPNT_QO_TOL_DFEAS", 1e-10, "MSK_DPAR_INTPNT_CO_TOL_REL_GAP", 1e-10, "MSK_DPAR_INTPNT_CO_TOL_PFEAS", 1e-10,
          "MSK_DPAR_INTPNT_CO_TOL_DFEAS", 1e-10)),
          map("quadratic", true, "less_digits", 1, "dual", true, "soc", true, "codegen", mosek_codegen, "binary", true, "discrete", true, "sos", false)));
    }
    if (enabled("SKIP_CLP_TESTS", hasConic("clp")))
      conics.Add(new Conic("clp", map("verbose", true),
          map("quadratic", false, "dual", true, "soc", false, "codegen", false, "binary", false, "discrete", false, "sos", false)));
    if (enabled("SKIP_CBC_TESTS", hasConic("cbc")))
      conics.Add(new Conic("cbc", map("verbose", true),
          map("quadratic", false, "dual", true, "soc", false, "codegen", false, "binary", true, "discrete", true, "sos", true)));
    if (enabled("SKIP_QRQP_TESTS", hasConic("qrqp"))) {
      var codegen = map("std", "c99");
      if (isDarwin) codegen["digits"] = 13;
      conics.Add(new Conic("qrqp", map("max_iter", 20, "print_header", false, "print_iter", false),
          map("quadratic", true, "dual", true, "soc", false, "codegen", codegen, "binary", false, "discrete", false, "sos", false)));
    }
    if (enabled("SKIP_PROXQP_TESTS", hasConic("proxqp"))) {
      foreach (string backend in new[] { "sparse", "dense" })
        conics.Add(new Conic("proxqp", map("proxqp", map("eps_abs", 1e-11, "max_iter", 1e4, "backend", backend)),
            map("quadratic", true, "dual", true, "soc", false, "codegen", false, "binary", false, "discrete", false, "sos", false)));
    }
    if (enabled("SKIP_QPALM_TESTS", hasConic("qpalm"))) {
      double eps = 1e-8;
      conics.Add(new Conic("qpalm", map("qpalm", map("eps_abs", eps, "eps_rel", eps, "eps_abs_in", eps, "eps_prim_inf", eps)),
          map("quadratic", true, "dual", true, "soc", false, "codegen", false, "binary", false, "discrete", false, "sos", false)));
    }
    if (enabled("SKIP_HIGHS_TESTS", hasConic("highs"))) {
      var highs_options = isWindows ? new List<string>() : new List<string> { "-Wno-unused-function" };
      var codegen = map("extralibs", new List<string> { "highs" }, "extra_options", highs_options, "std", "c99");
      conics.Add(new Conic("highs", map("highs", map("primal_feasibility_tolerance", 1e-7, "solver", "choose", "output_flag", false, "ipm_iteration_limit", 50000)),
          map("quadratic", true, "dual", true, "soc", false, "codegen", codegen, "binary", false, "discrete", true, "sos", false)));
    }
    if (enabled("SKIP_DAQP_TESTS", hasConic("daqp")))
      conics.Add(new Conic("daqp", map("daqp", map()),
          map("quadratic", true, "dual", true, "soc", false, "codegen", map("extralibs", new List<string> { "daqp" }, "std", "c99"), "binary", true, "discrete", false, "sos", false)));
    if (enabled("SKIP_PIQP_TESTS", hasConic("piqp"))) {
      foreach (string kkt in new[] { "sparse_ldlt", "dense_cholesky" })
        conics.Add(new Conic("piqp", map("piqp", map("kkt_solver", kkt)),
            map("quadratic", true, "dual", false, "soc", false, "codegen", false, "binary", false, "discrete", false, "sos", false)));
    }
  }

  static readonly List<string> fatrop_flags = isWindows ? new List<string>() : new List<string> { "-Wno-strict-prototypes" };

  static double at(DM x, long i) => (double)x[i];

  /// <summary>Python assertAlmostEqual(a, b, places)</summary>
  static void assertAlmostEqual(double a, double b, int places, string msg) =>
    assertTrue(Math.Abs(a - b) <= 0.5 * Math.Pow(10, -places), msg + ": " + a + " != " + b + " within " + places + " places");

  static void checkarrayDigits(DM a, DM b, int digits, string msg) => checkarray(a, b, 0.5 * Math.Pow(10, -digits), msg);

  static Dictionary<string, DM> solverIn(DM H, DM G, DM A, DM LBX, DM UBX, DM LBA, DM UBA) {
    var r = new Dictionary<string, DM> { ["h"] = H, ["g"] = G, ["a"] = A, ["lbx"] = LBX, ["ubx"] = UBX };
    if (LBA is not null) r["lba"] = LBA;
    if (UBA is not null) r["uba"] = UBA;
    return r;
  }

  static Dictionary<string, Sparsity> hA(DM H, DM A) =>
    new Dictionary<string, Sparsity> { ["h"] = H.sparsity(), ["a"] = A.sparsity() };

  static DM fill(int n, double v) => Enumerable.Repeat(v, n).ToArray();

  static void checkSerialize(Function F, IDictionary<string, DM> inputs) {
    string serialized = F.serialize(new Dictionary<string, object> { ["debug"] = true });
    // validate_serialization(serialized): python-only (serialization_scheme module)
    Function F2 = Function.deserialize(serialized);
    Dictionary<string, DM> Fout = F.call(inputs), Fout2 = F2.call(inputs);
    assertTrue(F.name_out().SequenceEqual(F2.name_out()), "name_out");
    foreach (string k in F.name_out()) checkarray(Fout[k], Fout2[k], 1e-16, k);
  }

  static List<string> strings(IDictionary<string, object> o, string k) =>
    o.TryGetValue(k, out object v) ? ((IEnumerable<string>)v).ToList() : new List<string>();

  /// <summary>helpers.py check_codegen: only under RunSlow; returns the external Function</summary>
  Function checkCodegen(Function F, IDictionary<string, DM> inputs, IDictionary<string, object> cg,
                        IDictionary<string, object> opts = null) {
    if (!RunSlow) return null;
    IDictionary<string, object> o = cg ?? new Dictionary<string, object>();
    string std = o.TryGetValue("std", out object sv) ? (string)sv : "c89";
    int digits = o.TryGetValue("digits", out object dv) ? Convert.ToInt32(dv) : 15;
    string dir = Path.Combine(Path.GetTempPath(), "casadi_cg" + Guid.NewGuid().ToString("N"));
    Directory.CreateDirectory(dir);
    string name = "codegen_" + Guid.NewGuid().ToString("N");
    var g = ca.CodeGenerator(name, opts ?? new Dictionary<string, object>());
    g.add(F);
    g.generate(dir + Path.DirectorySeparatorChar);
    string libdir = GlobalOptions.getCasadiPath();
    string includedir = GlobalOptions.getCasadiIncludePath();
    if (libdir.Length == 0 || includedir.Length == 0) Check.Skip("casadi include path not set (build tree)");
    var cmd = new List<string> { "-pedantic", "-std=" + std, "-fPIC", "-shared", "-Wall", "-Werror", "-Wextra",
      "-I" + includedir, "-I" + includedir + "/highs" };
    foreach (string e in strings(o, "extra_include")) cmd.Add("-I" + includedir + "/" + e);
    cmd.AddRange(new[] { "-Wno-dangling-pointer", "-Wno-unknown-warning-option", "-Wno-unknown-pragmas", "-Wno-long-long",
      "-Wno-unused-parameter", "-O3", Path.Combine(dir, name + ".c"), "-o", Path.Combine(dir, name + ".so") });
    var libdirs = new List<string> { libdir };
    libdirs.AddRange(strings(o, "extralibdirs"));
    foreach (string d in libdirs) cmd.AddRange(new[] { "-L" + d, "-Wl,-rpath," + d });
    foreach (string l in strings(o, "extralibs")) {
      bool found = l.StartsWith("-") || l.Contains(".");
      foreach (string d in libdirs) found |= File.Exists(Path.Combine(d, "lib" + l + ".so")) || File.Exists(Path.Combine(d, "lib" + l + ".a"));
      if (!found) {
        message("check_codegen skipped: lib" + l + " not shipped in " + string.Join(", ", libdirs));
        return null;
      }
      cmd.Add(l.StartsWith("-") || l.Contains(".") ? l : "-l" + l);
    }
    cmd.AddRange(strings(o, "extra_options"));
    var psi = new ProcessStartInfo("gcc");
    foreach (string a in cmd) psi.ArgumentList.Add(a);
    using (var p = Process.Start(psi)) {
      p.WaitForExit();
      assertEquals(0, p.ExitCode, "compile: gcc " + string.Join(" ", cmd));
    }
    // check_symbols: omitted (nm-based symbol prefix check of helpers.py)
    Function F2 = ca.external(F.name(), Path.Combine(dir, name + ".so"));
    Dictionary<string, DM> r = F.call(inputs), r2 = F2.call(inputs);
    foreach (string k in F.name_out()) checkarrayDigits(r2[k], r[k], digits, "codegen " + k);
    return F2;
  }

  public void test_binary() {
    DM.rng(1);
    int n = 10, m = 8, ms = 5;
    MX x = MX.sym("x", n);
    DM M = DM.rand(n, n);
    DM H = ca.mtimes(M, M.T()) * 0.5;
    DM f = DM.rand(n, 1) * 100;
    DM A = DM.rand(m, n);
    DM bupper = DM.rand(m, 1) * 20;
    DM blower = DM.rand(m, 1) * -20;

    double[] lbx = new double[n], ubx = new double[n];
    bool[] discrete = new bool[n];
    for (int i = 0; i < n; i++) {
      lbx[i] = i < ms ? 0 : -10;
      ubx[i] = i < ms ? 1 : 10;
      discrete[i] = i < ms;
    }
    var solver_in = new Dictionary<string, DM> { ["lbx"] = lbx, ["ubx"] = ubx, ["lbg"] = blower, ["ubg"] = bupper };
    foreach (Conic c in conics) {
      if (!c.has("binary")) continue;
      if (c.name != "daqp") {
        Function solver = ca.qpsol("solver", c.name, new Dictionary<string, MX> { ["f"] = ca.mtimes((MX)f.T(), x), ["x"] = x, ["g"] = ca.mtimes((MX)A, x) },
          map("discrete", discrete));
        solver.call(solver_in);
        if (c.codegen() != null) checkCodegen(solver, solver_in, c.codegen());
      }
      if (!c.has("quadratic")) continue;
      if (c.name != "highs") {
        MX obj = ca.mtimes(ca.mtimes(x.T(), (MX)H), x) * 0.5 + ca.mtimes((MX)f.T(), x);
        Function solver = ca.qpsol("solver", c.name, new Dictionary<string, MX> { ["f"] = obj, ["x"] = x, ["g"] = ca.mtimes((MX)A, x) },
          map("discrete", discrete));
        solver.call(solver_in);
        int digits = c.name == "mosek" ? 4 : 15;
        if (c.codegen() != null) {
          var cg = new Dictionary<string, object>(c.codegen()) { ["digits"] = digits };
          checkCodegen(solver, solver_in, cg);
        }
      }
    }
  }

  public void test_mosek_codegen_option_precision() {
    Check.Assume(hasConic("mosek"), "no conic mosek");
    var options = map("MSK_DPAR_MIO_TOL_ABS_RELAX_INT", 1.234567890123456e-7, "MSK_DPAR_OPTIMIZER_MAX_TIME", 1234567890.1234567);
    Function solver = ca.conic("solver", "mosek", new Dictionary<string, Sparsity> { ["h"] = Sparsity.dense(1, 1) }, map("mosek", options));
    var cg = ca.CodeGenerator("mosek_precision.c");
    cg.add(solver);
    var emitted = new Dictionary<string, string>();
    foreach (Match mt in Regex.Matches(cg.dump(), "MSK_putnadouparam\\(d->task, \"([^\"]+)\", ([^)]+)\\)"))
      emitted[mt.Groups[1].Value] = mt.Groups[2].Value;
    assertTrue(new HashSet<string>(options.Keys).SetEquals(emitted.Keys), "emitted keys");
    foreach (var e in options)
      assertEquals((double)e.Value, double.Parse(emitted[e.Key], System.Globalization.CultureInfo.InvariantCulture), e.Key);
  }

  public void test_opti() {
    foreach (Conic c in conics) {
      if (c.name == "qrqp") continue;
      if (c.name == "daqp") continue;
      var opti = ca.Opti("conic");
      MX x = opti.variable();
      opti.minimize(x);
      opti.subject_to(-2.0 <= (x <= 2.0));
      opti.solver(c.name, c.opts());
      OptiSol sol = opti.solve();
      checkarrayDigits(sol.value(x), -2, 5, c.name);
    }
  }

  public void test_sos() {
    var H = ca.DM(4, 4);
    DM G = ca.vertcat(-1, -2, -3, -1);
    DM A = ca.DM(new double[,] { { -1, 1, 1, 10 }, { 1, -3, 1, 0 }, { 0, 1, 0, -3.5 } });
    DM LBA = ca.vertcat(-inf, -inf, 0);
    DM UBA = ca.vertcat(20, 30, 0);
    bool[] discrete = { false, true, true, true };
    DM LBX = ca.vertcat(0, 0, 0, 2);
    DM UBX = ca.vertcat(40, inf, inf, 3);
    var sos_groups = new List<long[]> { new long[] { 2, 3 } };
    var sos_weights = new List<double[]> { new double[] { 25.0, 18.0 } };

    foreach (Conic c in conics) {
      if (!c.has("discrete")) continue;
      Console.WriteLine("test_sos " + c);
      var options = c.opts();
      options["discrete"] = discrete;
      Function solver = ca.conic("mysolver", c.name, hA(H, A), options);
      var solver_in = solverIn(H, G, A, LBX, UBX, LBA, UBA);
      var solver_out = solver.call(solver_in);
      checkarray(solver_out["x"], ca.vertcat(29, 7, 22, 2));

      if (!c.has("sos")) continue;
      options["sos_groups"] = sos_groups;
      solver = ca.conic("mysolver", c.name, hA(H, A), options);
      solver_out = solver.call(solver_in);
      checkarray(solver_out["x"], ca.vertcat(40, 7, 0, 2));

      options["sos_weights"] = sos_weights;
      solver = ca.conic("mysolver", c.name, hA(H, A), options);
      solver_out = solver.call(solver_in);
      checkarray(solver_out["x"], ca.vertcat(40, 7, 0, 2));
      Console.WriteLine(solver_out["cost"]);
      checkSerialize(solver, solver_in);
    }
  }

  public void test_milp() {
    var H = ca.DM(2, 2);
    DM G = ca.vertcat(-1, -1);
    DM A = ca.DM(new double[,] { { -2, 2 }, { -8, 10 } });
    DM LBA = ca.vertcat(1, -inf);
    DM UBA = ca.vertcat(inf, 13);
    DM LBX = fill(2, 0);
    DM UBX = fill(2, inf);
    bool[] discrete = { true, true };

    foreach (Conic c in conics) {
      if (!c.has("discrete")) continue;
      Console.WriteLine("test_milp " + c);
      var options = c.opts();
      options["discrete"] = discrete;
      Function solver = ca.conic("mysolver", c.name, hA(H, A), options);
      var solver_in = solverIn(H, G, A, LBX, UBX, LBA, UBA);
      var solver_out = solver.call(solver_in);
      checkarray(solver_out["x"], ca.vertcat(1, 2));
      checkarray(solver_out["cost"], ca.vertcat(-3));
      checkSerialize(solver, solver_in);
    }

    H = ca.DM(3, 3);
    G = ca.vertcat(-1, -1, -5);
    A = ca.DM(new double[,] { { -2, 2, 0 }, { -8, 10, 0 }, { 0, 0, 7 } });
    LBA = ca.vertcat(1, -inf, -4.0 * 7);
    UBA = ca.vertcat(inf, 13, 4.0 * 7);
    LBX = ca.vertcat(0, 0, -5);
    UBX = ca.vertcat(inf, inf, 5);
    discrete = new[] { true, true, false };

    foreach (Conic c in conics) {
      if (!c.has("discrete")) continue;
      Console.WriteLine("test_milp " + c);
      var options = c.opts();
      options["discrete"] = discrete;
      Function solver = ca.conic("mysolver", c.name, hA(H, A), options);
      var solver_out = solver.call(solverIn(H, G, A, LBX, UBX, LBA, UBA));
      checkarray(solver_out["x"], ca.vertcat(1, 2, 4));
      checkarray(solver_out["cost"], ca.vertcat(-23));
    }
  }

  public void test_general_unconstrained() {
    DM H = ca.sparsify(ca.DM(new double[,] { { 1, 0 }, { 0, 1 } }));
    DM G = ca.vertcat(-0.7, -2.3);
    DM A = DM.zeros(0, 2);
    DM LBA = DM.zeros(0, 1);
    DM UBA = DM.zeros(0, 1);
    DM LBX = fill(2, -inf);
    DM UBX = fill(2, inf);

    foreach (Conic c in conics) {
      if (!c.has("quadratic")) continue;
      Console.WriteLine("test_general_unconstrained " + c);
      Function solver = ca.conic("mysolver", c.name, hA(H, A), c.opts());
      int less_digits = c.lessDigits();
      var solver_in = solverIn(H, G, A, LBX, UBX, LBA, UBA);
      var solver_out = solver.call(solver_in);
      assertTrue((bool)solver.stats()["success"], c.name);
      assertAlmostEqual(at(solver_out["x"], 0), 0.7, Math.Max(1, 6 - less_digits), c.name);
      assertAlmostEqual(at(solver_out["x"], 1), 2.3, Math.Max(1, 6 - less_digits), c.name);
      if (c.codegen() != null) checkCodegen(solver, solver_in, c.codegen());
      checkSerialize(solver, solver_in);
    }
  }

  public void test_general_convex_dense() {
    message("Convex dense QP with solvers: " + string.Join(" ", conics.Select(c => c.name)));
    DM H = ca.DM(new double[,] { { 1, -1 }, { -1, 2 } });
    DM G = ca.vertcat(-2, -6);
    DM A = ca.DM(new double[,] { { 1, 1 }, { -1, 2 }, { 2, 1 } });
    DM LBA = fill(3, -inf);
    DM UBA = ca.vertcat(2, 2, 3);
    DM LBX = fill(2, 0);
    DM UBX = fill(2, inf);

    foreach (Conic c in conics) {
      if (c.name == "nlpsol") continue;
      if (!c.has("quadratic")) continue;
      Console.WriteLine("test_general_convex_dense " + c);
      if (c.name == "daqp") continue;
      Function solver = ca.conic("mysolver", c.name, hA(H, A), c.opts());
      try {
        solver.stats();
      } catch (SwigException) {
        // nosegfault
      }
      int ld = c.lessDigits();
      string s = c.name;
      var solver_in = solverIn(H, G, A, LBX, UBX, LBA, UBA);
      var solver_out = solver.call(solver_in);
      assertTrue((bool)solver.stats()["success"]);
      assertEquals("SOLVER_RET_SUCCESS", (string)solver.stats()["unified_return_status"]);
      assertAlmostEqual(at(solver_out["x"], 0), 2.0 / 3, Math.Max(1, 6 - ld), s);
      assertAlmostEqual(at(solver_out["x"], 1), 4.0 / 3, Math.Max(1, 6 - ld), s);
      if (c.has("dual")) assertAlmostEqual(at(solver_out["lam_x"], 0), 0, Math.Max(1, 6 - ld), s);
      if (c.has("dual")) assertAlmostEqual(at(solver_out["lam_x"], 1), 0, Math.Max(1, 6 - ld), s);
      if (c.has("dual")) checkarrayDigits(solver_out["lam_a"], ca.vertcat(3 + 1.0 / 9, 4.0 / 9, 0), Math.Max(1, 6 - ld), s);
      assertAlmostEqual(at(solver_out["cost"], 0), -8 - 2.0 / 9, Math.Max(1, 6 - ld), s);

      solver_in["h"] = H * 4;
      solver_out = solver.call(solver_in);
      assertAlmostEqual(at(solver_out["x"], 0), 1, Math.Max(1, 3 - ld), s);
      assertAlmostEqual(at(solver_out["x"], 1), 1, Math.Max(1, 3 - ld), s);
      assertAlmostEqual((double)solver_out["cost"], -6, Math.Max(1, 6 - ld), s);
      if (c.has("dual")) assertAlmostEqual(at(solver_out["lam_x"], 0), 0, Math.Max(1, 6 - ld), s);
      if (c.has("dual")) assertAlmostEqual(at(solver_out["lam_x"], 1), 0, Math.Max(1, 6 - ld), s);
      if (c.has("dual")) checkarrayDigits(solver_out["lam_a"], ca.vertcat(2, 0, 0), Math.Max(1, 2 - ld), s);

      solver_in["h"] = 0;
      if (c.name.Contains("qcqp")) continue;
      solver_out = solver.call(solver_in);
      assertAlmostEqual(at(solver_out["x"], 0), 2.0 / 3, Math.Max(1, 6 - ld), s);
      assertAlmostEqual(at(solver_out["x"], 1), 4.0 / 3, Math.Max(1, 6 - ld), s);
      assertAlmostEqual((double)solver_out["cost"], -9 - 1.0 / 3, Math.Max(1, 6 - ld), s);
      if (c.has("dual")) assertAlmostEqual(at(solver_out["lam_x"], 0), 0, Math.Max(1, 6 - ld), s);
      if (c.has("dual")) assertAlmostEqual(at(solver_out["lam_x"], 1), 0, Math.Max(1, 6 - ld), s);
      if (c.has("dual")) checkarrayDigits(solver_out["lam_a"], ca.vertcat(10.0 / 3, 4.0 / 3, 0), Math.Max(1, 4 - ld), s);

      solver_in["lba"] = fill(3, -inf);
      solver_in["uba"] = fill(3, inf);
      solver_in["ubx"] = 5;
      if (str(c.options).Contains("worhp")) {
        // Java/Python pass the dict positionally (expected to fail); not expressible in typed C#
        continue;
      }
      solver_out = solver.call(solver_in);
      assertAlmostEqual(at(solver_out["x"], 0), 5, Math.Max(1, 6 - ld), s);
      assertAlmostEqual(at(solver_out["x"], 1), 5, Math.Max(1, 6 - ld), s);
      assertAlmostEqual((double)solver_out["cost"], -40, Math.Max(1, 5 - ld), s);
      if (c.has("dual")) assertAlmostEqual(at(solver_out["lam_x"], 0), 2, Math.Max(1, 6 - ld), s);
      if (c.has("dual")) assertAlmostEqual(at(solver_out["lam_x"], 1), 6, Math.Max(1, 6 - ld), s);
      if (c.has("dual")) checkarrayDigits(solver_out["lam_a"], ca.vertcat(0, 0, 0), Math.Max(1, 4 - ld), s);
    }
  }

  public void test_general_convex_sparse() {
    message("Convex sparse QP with solvers: " + string.Join(" ", conics.Select(c => c.name)));
    DM H = ca.diag(ca.vertcat(2, 1, 0.2, 0.7, 1.3));
    H[1, 2] = 0.1;
    H[2, 1] = 0.1;
    DM G = ca.vertcat(-2, -6, 1, 0, 0);
    DM A = ca.sparsify(ca.DM(new double[,] { { 1, 0, 0.1, 0.7, -1 }, { 0.1, 2, -0.3, 4, 0.1 } }));
    DM LBA = ca.vertcat(-inf);
    DM UBA = ca.vertcat(2, 2);
    DM LBX = fill(5, 0);
    DM UBX = fill(5, inf);

    foreach (Conic c in conics) {
      if (!c.has("quadratic")) continue;
      Console.WriteLine("test_general_convex_sparse " + c);
      Function solver = ca.conic("mysolver", c.name, hA(H, A), c.opts());
      int ld = c.lessDigits();
      var solver_in = solverIn(H, G, A, LBX, UBX, LBA, UBA);
      var solver_out = solver.call(solver_in);
      checkarrayDigits(solver_out["x"], ca.vertcat(0.873908, 0.95630465, 0, 0, 0), Math.Max(1, 6 - ld), c.name);
      if (c.has("dual")) checkarrayDigits(solver_out["lam_x"], ca.vertcat(0, 0, -0.339076, -10.0873907, -0.252185), Math.Max(1, 6 - ld), c.name);
      if (c.has("dual")) checkarrayDigits(solver_out["lam_a"], ca.vertcat(0, 2.52184767), Math.Max(1, 6 - ld), c.name);
      assertAlmostEqual(at(solver_out["cost"], 0), -6.264669320767, Math.Max(1, 6 - ld), c.name);
      if (c.codegen() != null) checkCodegen(solver, solver_in, c.codegen());
    }
  }

  public void test_general_nonconvex_dense() {
    message("Non convex dense QP with solvers: " + string.Join(" ", conics.Select(c => c.name)));
    DM H = ca.DM(new double[,] { { 1, -1 }, { -1, -2 } });
    DM G = ca.vertcat(-2, -6);
    DM A = ca.DM(new double[,] { { 1, 1 }, { -1, 2 }, { 2, 1 } });
    DM LBA = fill(3, -inf);
    DM UBA = ca.vertcat(2, 2, 3);
    DM LBX = fill(2, 0);
    DM UBX = fill(2, inf);

    foreach (Conic c in conics) {
      if (!c.has("quadratic")) continue;
      Console.WriteLine("test_general_nonconvex_dense " + c);
      Function solver = ca.conic("mysolver", c.name, hA(H, A), c.opts());
      var solver_in = solverIn(H, G, A, LBX, UBX, LBA, UBA);
      // Python passes the dict positionally (solver(solver_in)), which fails input conversion for every plugin;
      // Java mirrors that with invoke(Object[]{map}). Passing a dictionary as a positional DM is a compile error in C#.
      Check.Skip("not applicable: positional dict argument is a compile-time error in C#");
    }
  }

  public void test_equality() {
    message("Regression 452 test: equality constraints give wrong multipliers");
    DM H = ca.DM(new double[,] { { 1, -1 }, { -1, 2 } });
    DM G = ca.vertcat(-2, -6);

    foreach (Conic c in conics) {
      if (!c.has("quadratic")) continue;
      if (c.name.Contains("ooqp") || c.name.Contains("proxqp")) continue;
      Console.WriteLine("test_equality " + c);
      Function solver = ca.conic("mysolver", c.name,
        new Dictionary<string, Sparsity> { ["h"] = H.sparsity(), ["a"] = Sparsity.dense(3, 2) }, c.opts());
      int ld = c.lessDigits();
      string s = c.name;
      DM A = ca.DM(new double[,] { { 1, 1 }, { -1, 2 }, { 2, 1 } });
      DM LBA = fill(3, -inf);
      DM UBA = ca.vertcat(2, 2, 3);
      DM LBX = ca.vertcat(0.5, 0);
      DM UBX = ca.vertcat(0.5, inf);
      var solver_in = solverIn(H, G, A, LBX, UBX, LBA, UBA);
      if (str(c.options).Contains("worhp")) {
        // Java/Python pass the dict positionally (expected to fail); not expressible in typed C#
        continue;
      }
      var solver_out = solver.call(solver_in);
      assertAlmostEqual(at(solver_out["x"], 0), 0.5, Math.Max(1, 6 - ld), s);
      assertAlmostEqual(at(solver_out["x"], 1), 1.25, Math.Max(1, 6 - ld), s);
      if (c.has("dual")) assertAlmostEqual(at(solver_out["lam_x"], 0), 4.75, Math.Max(1, 6 - ld), s);
      if (c.has("dual")) assertAlmostEqual(at(solver_out["lam_x"], 1), 0, Math.Max(1, 6 - ld), s);
      if (c.has("dual")) checkarrayDigits(solver_out["lam_a"], ca.vertcat(0, 2, 0), Math.Max(1, 6 - ld), s);
      if (c.codegen() != null) checkCodegen(solver, solver_in, c.codegen());
      assertAlmostEqual(at(solver_out["cost"], 0), -7.4375, Math.Max(1, 6 - ld), s);

      solver_in["h"] = H;
      solver_in["g"] = G;
      solver_in["a"] = ca.DM(new double[,] { { 1, 1 }, { -1, 2 }, { 2, 1 } });
      solver_in["lbx"] = fill(2, -inf);
      solver_in["ubx"] = fill(2, inf);
      solver_in["lba"] = ca.vertcat(2, -inf, -inf);
      solver_in["uba"] = ca.vertcat(2, inf, inf);
      solver_out = solver.call(solver_in);
      assertAlmostEqual(at(solver_out["x"], 0), 0.4, Math.Max(1, 4 - ld), s);
      assertAlmostEqual(at(solver_out["x"], 1), 1.6, Math.Max(1, 4 - ld), s);
      if (c.has("dual")) assertAlmostEqual(at(solver_out["lam_x"], 0), 0, Math.Max(1, 5 - ld), s);
      if (c.has("dual")) assertAlmostEqual(at(solver_out["lam_x"], 1), 0, Math.Max(1, 5 - ld), s);
      if (c.has("dual")) checkarrayDigits(solver_out["lam_a"], ca.vertcat(3.2, 0, 0), Math.Max(1, 5 - ld), s);
      assertAlmostEqual(at(solver_out["cost"], 0), -8.4, Math.Max(1, 5 - ld), s);
      if (c.codegen() != null) checkCodegen(solver, solver_in, c.codegen());
    }
  }

  public void test_degenerate_hessian() {
    message("Degenerate hessian");
    DM H = ca.sparsify(ca.DM(new double[,] { { 1, -1, 0 }, { -1, 2, 0 }, { 0, 0, 0 } }));
    DM G = ca.vertcat(-2, -6, 1);
    DM A = ca.DM(new double[,] { { 1, 1, 1 } });
    DM LBA = ca.vertcat(0.5);
    DM UBA = ca.vertcat(0.5);
    DM LBX = ca.vertcat(-10);
    DM UBX = ca.vertcat(10);

    foreach (Conic c in conics) {
      if (c.name.Contains("qrqp")) continue;
      if (!c.has("quadratic")) continue;
      // Python mutates the shared table here (dump_in leaks into later tests); kept local
      var qp_options = c.opts();
      qp_options["dump_in"] = true;
      if (c.name.Contains("qcqp")) continue;
      if (c.name == "daqp") continue;
      Console.WriteLine("test_degenerate_hessian " + c.name + " " + str(qp_options));
      Function solver = ca.conic("mysolver", c.name, hA(H, A), qp_options);
      int ld = c.lessDigits();
      var solver_in = solverIn(H, G, A, LBX, UBX, LBA, UBA);
      var solver_out = solver.call(solver_in);
      checkarrayDigits(solver_out["x"], ca.vertcat(5.5, 5, -10), Math.Max(1, 4 - ld), c.name);
      if (c.has("dual")) checkarrayDigits(solver_out["lam_x"], ca.vertcat(0, 0, -2.5), Math.Max(1, 4 - ld), c.name);
      if (c.has("dual")) checkarrayDigits(solver_out["lam_a"], ca.vertcat(1.5), Math.Max(1, 4 - ld), c.name);
      assertAlmostEqual(at(solver_out["cost"], 0), -38.375, Math.Max(1, 5 - ld), c.name);
      if (c.codegen() != null) checkCodegen(solver, solver_in, c.codegen());
    }
  }

  public void test_no_inequality() {
    message("No inequalities present");
    DM H = ca.DM(new double[,] { { 1, -1 }, { -1, 2 } });
    DM G = ca.vertcat(-2, -6);
    DM A = ca.DM(new double[,] { { 1, 1 } });
    DM LBA = ca.vertcat(0.5);
    DM UBA = ca.vertcat(0.5);
    DM LBX = ca.vertcat(-10);
    DM UBX = ca.vertcat(10);

    foreach (Conic c in conics) {
      if (!c.has("quadratic")) continue;
      Console.WriteLine("test_no_inequality " + c);
      Function solver = ca.conic("mysolver", c.name, hA(H, A), c.opts());
      int ld = c.lessDigits();
      string s = c.name;
      var solver_out = solver.call(solverIn(H, G, A, LBX, UBX, LBA, UBA));
      assertAlmostEqual(at(solver_out["x"], 0), -0.5, Math.Max(1, 6 - ld), s);
      assertAlmostEqual(at(solver_out["x"], 1), 1, Math.Max(1, 6 - ld), s);
      if (c.has("dual")) assertAlmostEqual(at(solver_out["lam_x"], 0), 0, Math.Max(1, 6 - ld), s);
      if (c.has("dual")) assertAlmostEqual(at(solver_out["lam_x"], 1), 0, Math.Max(1, 6 - ld), s);
      if (c.has("dual")) checkarrayDigits(solver_out["lam_a"], ca.vertcat(3.5), Math.Max(1, 6 - ld), s);
      assertAlmostEqual(at(solver_out["cost"], 0), -3.375, Math.Max(1, 6 - ld), s);
    }
  }

  public void test_no_A() {
    message("No A present");
    DM H = ca.DM(new double[,] { { 1, -1 }, { -1, 2 } });
    DM G = ca.vertcat(-2, -6);
    var A = ca.DM(0, 2);
    var LBA = ca.DM(0, 1);
    var UBA = ca.DM(0, 1);
    DM LBX = ca.vertcat(-10);
    DM UBX = ca.vertcat(10);

    foreach (Conic c in conics) {
      if (!c.has("quadratic")) continue;
      Console.WriteLine("test_no_A " + c);
      Function solver = ca.conic("mysolver", c.name, hA(H, A), c.opts());
      int ld = c.lessDigits();
      var solver_in = solverIn(H, G, A, LBX, UBX, LBA, UBA);
      var solver_out = solver.call(solver_in);
      checkarrayDigits(solver_out["x"], ca.vertcat(10, 8), Math.Max(1, 3 - ld), c.name);
      if (c.has("dual")) checkarrayDigits(solver_out["lam_x"], ca.vertcat(0, 0), Math.Max(1, 4 - ld), c.name);
      if (c.has("dual")) checkarrayDigits(solver_out["lam_a"], ca.vertcat(), Math.Max(1, 5 - ld), c.name);
      assertAlmostEqual(at(solver_out["cost"], 0), -34, Math.Max(1, 5 - ld), c.name);
      if (c.codegen() != null) checkCodegen(solver, solver_in, c.codegen());
    }
  }

  public void test_standard_form() {
    DM H = ca.DM(new double[,] { { 1, -1 }, { -1, 2 } });
    DM G = ca.vertcat(-2, -6);
    DM A = ca.vertcat(1, 1).T();
    DM LBA = ca.vertcat(-inf);
    DM UBA = ca.vertcat(1);
    DM LBX = ca.vertcat(-10);
    DM UBX = ca.vertcat(10);

    foreach (Conic c in conics) {
      if (!c.has("quadratic")) continue;
      Console.WriteLine("test_standard_form " + c);
      Function solver = ca.conic("mysolver", c.name, hA(H, A), c.opts());
      int ld = c.lessDigits();
      var solver_in = solverIn(H, G, A, LBX, UBX, LBA, UBA);
      var solver_out = solver.call(solver_in);
      checkarrayDigits(solver_out["x"], ca.vertcat(-0.2, 1.2), Math.Max(1, 3 - ld), c.name);
      if (c.has("dual")) checkarrayDigits(solver_out["lam_x"], ca.vertcat(0, 0), Math.Max(1, 4 - ld), c.name);
      if (c.has("dual")) checkarrayDigits(solver_out["lam_a"], ca.vertcat(3.4), Math.Max(1, 5 - ld), c.name);
      assertAlmostEqual(at(solver_out["cost"], 0), -5.1, Math.Max(1, 5 - ld), c.name);
      if (c.codegen() != null) checkCodegen(solver, solver_in, c.codegen());
    }
  }

  public void test_badscaling() {
    message("Badly scaled problem");
    int N = 50;
    double[] d = new double[N], r = new double[N];
    for (int i = 0; i < N; i++) {
      d[i] = i + 1;
      r[i] = i;
    }
    DM H = ca.diag(ca.DM(d));
    DM x0 = ca.DM(r);
    DM G = ca.mtimes(H, x0) * -1.0;
    var A = ca.DM(0, N);
    DM LBX = fill(N, -1000);
    DM UBX = fill(N, 1000);

    foreach (Conic c in conics) {
      Console.WriteLine(c);
      if (!c.has("quadratic")) continue;
      if (c.name.Contains("worhp")) continue;
      if (c.name.Contains("superscs")) continue;
      Console.WriteLine("test_badscaling " + c);
      Function solver = ca.conic("mysolver", c.name, hA(H, A), c.opts());
      int ld = c.lessDigits();
      var solver_out = solver.call(solverIn(H, G, A, LBX, UBX, null, null));
      checkarrayDigits(solver_out["x"], x0, Math.Max(1, 2 - ld), c.ToString());
      assertAlmostEqual(at(solver_out["cost"], 0), -0.5 * (double)ca.mtimes(new[] { x0.T(), H, x0 }), Math.Max(1, 3 - ld), c.name);
      if (c.has("dual")) checkarrayDigits(solver_out["lam_x"], DM.zeros(N, 1), Math.Max(1, 4 - ld), c.name);
    }
  }

  public void test_redundant() {
    message("Redundant constraints");
    DM H = ca.DM(new double[,] { { 1, -1, 0 }, { -1, 2, 0 }, { 0, 0, 0 } });
    DM G = ca.vertcat(-2, -6, 1);
    DM a = ca.vertcat(1, 0, 1);
    DM a_ = ca.vertcat(0, 1, -2);

    foreach (double[] w in new[] { new double[] { 0, 2 }, new double[] { 1, 1 }, new double[] { 0.1, 0.6 } }) {
      double w0 = w[0], w1 = w[1];
      DM A = ca.vertcat(a.T(), a_.T(), (a * w0 + a_ * w1).T());
      DM LBA = ca.vertcat(0, 0, 0);
      DM UBA = ca.vertcat(0.5, 0.3, w0 * 0.5 + w1 * 0.3);
      DM LBX = ca.vertcat(-10);
      DM UBX = ca.vertcat(10);

      foreach (Conic c in conics) {
        if (!c.has("quadratic")) continue;
        if (c.name.Contains("qcqp")) continue;
        Console.WriteLine("test_redundant " + c);
        if (c.name == "daqp") continue;
        Function solver = ca.conic("qpsol", c.name, hA(H, A), c.opts());
        int ld = c.lessDigits();
        var solver_out = solver.call(solverIn(H, G, A, LBX, UBX, LBA, UBA));
        checkarrayDigits(solver_out["x"], ca.vertcat(-0.19230768069, 1.6846153915, 0.692307690769276), Math.Max(1, 6 - ld), c.name);
        assertAlmostEqual(at(solver_out["cost"], 0), -5.850384678537, Math.Max(1, 5 - ld), c.name);
        if (c.has("dual")) checkarrayDigits(solver_out["lam_x"], ca.vertcat(0, 0, 0), Math.Max(1, 6 - ld), c.name);
        if (c.has("dual"))
          checkarrayDigits(ca.mtimes(A.T(), solver_out["lam_a"]), ca.vertcat(3.876923073076, 2.4384615365384965, -1), Math.Max(1, 6 - ld), c.name);
      }
    }
  }

  public void test_linear() {
    var H = ca.DM(2, 2);
    DM A = ca.DM(new double[,] { { -1, 1 }, { 1, 1 }, { 1, -2 } });
    DM LBA = ca.vertcat(-inf, 2, -inf);
    DM UBA = ca.vertcat(1, inf, 4);
    DM LBX = ca.vertcat(-inf, 0);
    DM UBX = ca.vertcat(inf, inf);
    DM G = ca.vertcat(2, 1);

    foreach (Conic c in conics) {
      if (c.name.Contains("qcqp")) continue;
      if (c.name == "daqp") continue;
      Console.WriteLine("test_linear " + c);
      Function solver = ca.conic("mysolver", c.name, hA(H, A), c.opts());
      int ld = c.lessDigits();
      var solver_in = solverIn(H, G, A, LBX, UBX, LBA, UBA);
      if (!(c.name == "superscs" && isDarwin)) checkSerialize(solver, solver_in);
      var solver_out = solver.call(solver_in);
      checkarrayDigits(solver_out["x"], ca.vertcat(0.5, 1.5), Math.Max(1, 5 - ld), c.name);
      if (c.has("dual")) checkarrayDigits(solver_out["lam_x"], ca.vertcat(0, 0), Math.Max(1, 5 - ld), c.name);
      if (c.has("dual")) checkarrayDigits(solver_out["lam_a"], ca.vertcat(0.5, -1.5, 0), Math.Max(1, 5 - ld), c.name);
      assertAlmostEqual(at(solver_out["cost"], 0), 2.5, Math.Max(1, 5 - ld), c.name);
      if (c.codegen() != null) checkCodegen(solver, solver_in, c.codegen());
    }
  }

  public void test_linear2() {
    var H = ca.DM(2, 2);
    DM A = ca.DM(new double[,] { { -1, 1 }, { 1, 1 }, { 1, -2 } });
    DM LBA = ca.vertcat(-inf, 2, -inf);
    DM UBA = ca.vertcat(1, inf, 4);
    DM LBX = ca.vertcat(-inf, 3);
    DM UBX = ca.vertcat(inf, 3);
    DM G = ca.vertcat(2.0, 1.0);

    foreach (Conic c in conics) {
      if (c.name.Contains("qcqp")) continue;
      if (c.name.Contains("nlp")) continue;
      if (c.name == "daqp") continue;
      Console.WriteLine("test_linear2 " + c);
      Function solver = ca.conic("msyolver", c.name, hA(H, A), c.opts());
      int ld = c.lessDigits();
      var solver_out = solver.call(solverIn(H, G, A, LBX, UBX, LBA, UBA));
      checkarrayDigits(solver_out["x"], ca.vertcat(2, 3), Math.Max(1, 5 - ld), c.name);
      if (c.has("dual")) checkarrayDigits(solver_out["lam_x"], ca.vertcat(0, -3), Math.Max(1, 5 - ld), c.name);
      if (c.has("dual")) checkarrayDigits(solver_out["lam_a"], ca.vertcat(2, 0, 0), Math.Max(1, 5 - ld), c.name);
      assertAlmostEqual(at(solver_out["cost"], 0), 7, Math.Max(1, 5 - ld), c.name);
    }
  }

  public void test_overconstrained() {
    SX x = SX.sym("x");
    var qp = new Dictionary<string, SX> { ["x"] = x, ["f"] = ca.pow(x - 1, 2), ["g"] = ca.vertcat(x, x, x) };

    foreach (Conic c in conics) {
      if (!c.has("quadratic")) continue;
      Console.WriteLine("test_overconstrained " + c);
      Function solver = ca.qpsol("mysolver", c.name, qp, c.opts());
      var solver_in = new Dictionary<string, DM> {
        ["lbx"] = ca.vertcat(-10), ["ubx"] = ca.vertcat(10), ["lbg"] = ca.vertcat(-10, -10, -10), ["ubg"] = ca.vertcat(10, 10, 10), ["x0"] = ca.vertcat(0) };
      var solver_out = solver.call(solver_in);
      assertAlmostEqual(at(solver_out["f"], 0), 0, 9, c.name);
      assertAlmostEqual(at(solver_out["x"], 0), 1, 5, c.name);
      if (c.codegen() != null) checkCodegen(solver, solver_in, c.codegen());
    }
  }

  public void test_qrqp_prints() {
    Check.Assume(hasConic("qrqp"), "no conic qrqp");
    MX x = MX.sym("x");
    Function solver = ca.qpsol("solver", "qrqp", new Dictionary<string, MX> { ["x"] = x, ["f"] = ca.pow(x - 1, 2) });
    string o = captureStdout(() => solver.call(new Dictionary<string, DM>()));
    assertTrue(o.Contains("last_tau"), o);
    assertTrue(o.Contains("Converged"), o);
    if (RunSlow) {
      Function F = checkCodegen(solver, new Dictionary<string, DM>(), map("std", "c99"), map());
      string out2 = captureStdout(() => F.call(new Dictionary<string, DM>()));
      assertFalse(out2.Contains("last_tau"), out2);
      assertFalse(out2.Contains("Converged"), out2);
      checkCodegen(solver, new Dictionary<string, DM>(), map("std", "c99"), map("verbose_runtime", true));
    }
  }

  public void test_hpipm() {
    Check.Assume(hasConic("hpipm"), "no conic hpipm");
    Check.Assume(hasConic("qpoases"), "no conic qpoases");
    double inf = 100;
    int N = 4;
    MX x1 = MX.sym("x1");
    MX x2 = MX.sym("x2");
    MX x = ca.vertcat(x1, x2);
    MX u = MX.sym("u");
    MX xdot = ca.vertcat(x1 * 0.6 - x2 * 1.11 + u * 0.3 - 0.03, x1 * 0.7 + 0.01);
    MX L = ca.pow(x1, 2) + ca.pow(x2, 2) * 3 + ca.pow(u, 2) * 7 - x1 * x2 * 0.4 - x1 * u * 0.3 + u - x1 - x2 * 2;
    var F = ca.Function("F", new[] { x, u }, new[] { x + xdot, L });
    F.jacobian_old(0, 0);

    var w = new List<SX>();
    List<double> w0 = new List<double>(), lbw = new List<double>(), ubw = new List<double>(), lbg = new List<double>(), ubg = new List<double>();
    SX J = 0;
    var g = new List<SX>();
    List<SX> Xs = SX.sym("X", 2, 1, N + 1);
    List<SX> Us = SX.sym("U", 1, 1, N + 1);
    for (int k = 0; k < N; k++) {
      w.Add(Xs[k]);
      if (k == 0) {
        lbw.AddRange(new[] { -inf, 1.0 });
        ubw.AddRange(new[] { inf, 1.0 });
        w0.AddRange(new[] { 0.3, 0.7 });
      } else if (k == 2) {
        lbw.AddRange(new[] { 0.0, -inf });
        ubw.AddRange(new[] { 0.0, inf });
        w0.AddRange(new[] { 0.0, 0.3 });
      } else {
        lbw.AddRange(new[] { -inf, -inf });
        ubw.AddRange(new[] { inf, inf });
        w0.AddRange(new[] { 0.0, 1.0 });
      }
      w.Add(Us[k]);
      lbw.Add(-inf);
      ubw.Add(inf);
      w0.Add(0.0);
      List<SX> r = F.call(new[] { Xs[k], Us[k] });
      SX xplus = r[0], l = r[1];
      J = J + l;
      g.Add((xplus - Xs[k + 1]) * 3);
      lbg.AddRange(new[] { 0.0, 0.0 });
      ubg.AddRange(new[] { 0.0, 0.0 });
      g.Add(Xs[k][1] * 0.1 - Us[k] * 0.05);
      lbg.Add(-0.5 * k - 0.1);
      ubg.Add(2.0);
    }
    SX Xl = Xs[Xs.Count - 1];
    g.Add(Xl[1] * 0.1);
    lbg.Add(0.1);
    ubg.Add(2.0);
    J = J + ca.mtimes(Xl.T(), Xl);
    w.Add(Xl);
    lbw.AddRange(new[] { -inf, -inf });
    ubw.AddRange(new[] { inf, inf });
    w0.AddRange(new[] { 0.0, 1.0 });

    var prob = new Dictionary<string, SX> { ["f"] = J, ["x"] = ca.vertcat(w), ["g"] = ca.vertcat(g) };
    var Jf = ca.Function("J", new[] { prob["x"] }, new[] { ca.jacobian(prob["g"], prob["x"]) });
    Jf.Call(w0.ToArray()).print_dense();

    Function solver_ref = ca.qpsol("solver", "qpoases", prob);
    var options = map("hpipm", map("iter_max", 100, "res_g_max", 1e-10, "res_b_max", 1e-10, "res_d_max", 1e-10, "res_m_max", 1e-10), "dump_in", true);
    Function solver = ca.qpsol("solver", "hpipm", prob, options);
    var args = new Dictionary<string, DM> {
      ["x0"] = w0.ToArray(), ["lbx"] = lbw.ToArray(), ["ubx"] = ubw.ToArray(), ["lbg"] = lbg.ToArray(), ["ubg"] = ubg.ToArray() };
    var sol_ref = solver_ref.call(args);
    var sol = solver.call(args);
    checkCodegen(solver, args, map("std", "c99", "extralibs", new List<string> { "hpipm", "blasfeo" }, "extra_options", fatrop_flags, "digits", 10));
    checkarray(sol_ref["x"], sol["x"]);
    checkarrayDigits(sol_ref["lam_g"], sol["lam_g"], 8, "");
    checkarrayDigits(sol_ref["lam_x"], sol["lam_x"], 8, "");
    checkarrayDigits(sol_ref["f"], sol["f"], 8, "");

    solver = ca.nlpsol("solver", "sqpmethod", prob, map("qpsol", "hpipm", "qpsol_options", options));
    sol = solver.call(args);
    checkarray(sol_ref["x"], sol["x"]);
    checkarrayDigits(sol_ref["lam_g"], sol["lam_g"], 8, "");
    checkarrayDigits(sol_ref["lam_x"], sol["lam_x"], 8, "");
    checkarrayDigits(sol_ref["f"], sol["f"], 8, "");
  }

  static DM mat(string a) {
    var rows = new List<double[]>();
    foreach (string i in a.Split('\n')) {
      if (i.Length == 0) continue;
      string[] e = i.Split('\t');
      rows.Add(e.Select(s => s.Trim().Length == 0 ? 0 : double.Parse(s.Trim(), System.Globalization.CultureInfo.InvariantCulture)).ToArray());
    }
    var m = new double[rows.Count, rows[0].Length];
    for (int r = 0; r < rows.Count; r++)
      for (int k = 0; k < rows[r].Length; k++) m[r, k] = rows[r][k];
    return ca.sparsify(m);
  }

  static DM vec(string a) =>
    a.Split('\n').Select(s => double.Parse(s.Trim(), System.Globalization.CultureInfo.InvariantCulture)).ToArray();

  public void test_hpipm_timevarying() {
    Check.Assume(hasConic("hpipm"), "no conic hpipm");
    Check.Assume(hasConic("qpoases"), "no conic qpoases");
    DM A = mat(string.Join("\n",
        "1\t0.2\t1\t-1\t0\t0\t0\t0\t0\t0\t0\t0",
        "-0.1\t0.4\t0\t0\t-1\t0\t0\t0\t0\t0\t0\t0",
        "0.3\t0.2\t0\t0\t0\t-1\t0\t0\t0\t0\t0\t0",
        "2\t0\t0.3\t0\t0\t0\t0\t0\t0\t0\t0\t0",
        "1\t1\t0.4\t0\t0\t0\t0\t0\t0\t0\t0\t0",
        "\t0\t0\t1\t4\t2\t1\t0.3\t-1\t0\t0\t0",
        "\t0\t0\t3\t1\t0\t1\t0.2\t0\t-1\t0\t0",
        "\t0\t0\t1\t1\t1\t1\t1\t0\t0\t0\t0",
        "\t0\t0\t0\t0\t0\t0\t0\t2\t4\t0\t-1",
        "\t0\t0\t0\t0\t0\t0\t0\t2\t3\t1\t0",
        "\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t3"));
    DM H = mat(string.Join("\n",
        "7\t0\t0.2\t0\t0\t0\t0\t0\t0\t0\t0\t0",
        "\t7\t0.3\t0\t0\t0\t0\t0\t0\t0\t0\t0",
        "0.2\t0.3\t1\t0\t0\t0\t0\t0\t0\t0\t0\t0",
        "\t0\t0\t3\t0\t0\t0\t1\t0\t0\t0\t0",
        "\t0\t0\t0\t2\t0.1\t0\t0.7\t0\t0\t0\t0",
        "\t0\t0\t0\t0.1\t1\t0\t1\t0\t0\t0\t0",
        "\t0\t0\t0\t0\t0\t1\t0.1\t0\t0\t0\t0",
        "\t0\t0\t1\t0.7\t1\t0.1\t2\t0\t0\t0\t0",
        "\t0\t0\t0\t0\t0\t0\t0\t6\t0\t1\t0",
        "\t0\t0\t0\t0\t0\t0\t0\t0\t6\t0\t0",
        "\t0\t0\t0\t0\t0\t0\t0\t1\t0\t4\t0",
        "\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t9", ""));
    var options = map("hpipm", map("iter_max", 100, "res_g_max", 1e-10, "res_b_max", 1e-10, "res_d_max", 1e-10, "res_m_max", 1e-10));
    Function solver = ca.conic("solver", "hpipm", hA(H, A), options);
    Function solver_ref = ca.conic("solver", "qpoases", hA(H, A));

    DM g = vec("1\n1\n0.2\n0.4\n1\n0.5\n0.3\n1\n0.6\n1\n1\n0.7");
    DM lbg = vec("0\n0\n0\n-2\n-2\n0\n0\n-2\n0\n-2\n-2");
    DM ubg = vec("0\n0\n0\n2\n2\n0\n0\n2\n0\n2\n2");
    DM lbx = vec("0.5\n0.2\n-1\n-1\n-1\n-1\n-1\n-1\n-1\n-1\n-1\n-1");
    DM ubx = vec("0.5\n0.2\n1\n1\n1\n1\n1\n1\n1\n1\n1\n1");

    var args = new Dictionary<string, DM> { ["a"] = A, ["h"] = H, ["lba"] = lbg, ["uba"] = ubg, ["g"] = g, ["lbx"] = lbx, ["ubx"] = ubx };
    var sol = solver.call(args);
    var sol_ref = solver_ref.call(args);
    checkarrayDigits(sol_ref["x"], sol["x"], 7, "");
    checkarrayDigits(sol_ref["lam_a"], sol["lam_a"], 8, "");
    checkarrayDigits(sol_ref["lam_x"], sol["lam_x"], 8, "");
    checkCodegen(solver, args, map("std", "c99", "extralibs", new List<string> { "hpipm", "blasfeo" }, "extra_options", fatrop_flags, "digits", 9));

    solver = ca.conic("solver", "hpipm", hA(H, A), options);
    var args2 = new Dictionary<string, DM>(args) { ["x0"] = sol["x"], ["lam_a0"] = sol["lam_a"], ["lam_x0"] = sol["lam_x"] };
    sol = solver.call(args2);
    checkarrayDigits(sol_ref["x"], sol["x"], 7, "");
    checkarrayDigits(sol_ref["lam_a"], sol["lam_a"], 8, "");
    checkarrayDigits(sol_ref["lam_x"], sol["lam_x"], 8, "");

    var args3 = new Dictionary<string, DM>(args) { ["x0"] = sol["x"], ["lam_a0"] = sol["lam_a"], ["lam_x0"] = sol["lam_x"] };
    var sol2 = solver.call(args3);
    checkarray(sol["x"], sol2["x"]);
    Console.WriteLine("codegen starts here");
    checkCodegen(solver, args3, map("std", "c99", "extralibs", new List<string> { "hpipm", "blasfeo" }, "extra_options", fatrop_flags, "digits", 9));
  }

  public void test_SOCP() {
    Check.Assume(hasNlpsol("ipopt"), "no nlpsol ipopt");
    foreach (Conic c in conics) {
      // Python mutates the shared table here (verbose leaks into later tests); kept local
      var qp_options = c.opts();
      qp_options["verbose"] = true;
      MX x = MX.sym("x");
      MX y = MX.sym("y");
      MX z = MX.sym("z");
      if (!c.has("soc")) continue;
      Console.WriteLine("test_SOCP " + c.name + " " + str(qp_options));

      MX h = ca.soc(ca.vertcat(x - 5, y - 7), 4);
      Function solver = ca.qpsol("msyolver", c.name, new Dictionary<string, MX> { ["h"] = h, ["x"] = ca.vertcat(x, y), ["f"] = x * 2 + y }, qp_options);
      var res = solver.call(new Dictionary<string, DM>());
      checkarrayDigits(res["x"], ca.vertcat(5 - 8 / Math.Sqrt(5), 7 - 4 / Math.Sqrt(5)), 7, c.name);
      checkarrayDigits(res["f"], 10 - 16 / Math.Sqrt(5) + 7 - 4 / Math.Sqrt(5), 7, c.name);

      double a = 1.3;
      h = ca.soc(ca.vertcat(x, y), x * a + 4);
      MX xy = ca.vertcat(x, y);
      Function nsolver = ca.nlpsol("mysolver", "ipopt", new Dictionary<string, MX> {
        ["f"] = x * 2 + y, ["x"] = xy, ["g"] = ca.dot(xy, xy) - ca.pow(x * a + 4, 2) });
      nsolver.call(new Dictionary<string, DM> { ["ubg"] = 0 });
      solver = ca.qpsol("msyolver", c.name, new Dictionary<string, MX> { ["h"] = h, ["x"] = xy, ["f"] = x * 2 + y }, qp_options);
      res = solver.call(new Dictionary<string, DM>());
      double xs = -(8 * Math.Sqrt(5 - a * a) + 4 * Math.Pow(a, 3) - 20 * a) / (Math.Pow(a, 4) - 6 * a * a + 5);
      double ys = 4 * Math.Sqrt(5 - a * a) / (a * a - 5);
      checkarrayDigits(res["f"], 2 * xs + ys, 5, c.name);
      checkarrayDigits(res["x"], ca.vertcat(xs, ys), new[] { "cplex", "gurobi", "xpress" }.Contains(c.name) ? 4 : 5, c.name);

      h = ca.diagcat(ca.soc(ca.vertcat(x - 5, y - 7), 4), ca.soc(ca.vertcat(x / 6, y / 5), 1));
      solver = ca.qpsol("msyolver", c.name, new Dictionary<string, MX> { ["h"] = h, ["x"] = xy, ["f"] = x * 2 + y }, qp_options);
      res = solver.call(new Dictionary<string, DM>());
      checkarrayDigits(res["x"], ca.vertcat(1.655450403084473, 4.805919456574478), 5, c.name);

      h = ca.diagcat(
          ca.soc(ca.vertcat(x * -13 + y * 3 + z * 5 - 3, x * -12 + y * 12 - z * 6 - 2), x * -12 - y * 6 + z * 5 - 12),
          ca.soc(ca.vertcat(x * -3 + y * 6 + z * 2, x + y * 9 + z * 2 + 3, -x - y * 19 + z * 3 - 42), x * -3 + y * 6 - z * 10 + 27));
      solver = ca.qpsol("msyolver", c.name, new Dictionary<string, MX> { ["h"] = h, ["x"] = ca.vertcat(x, y, z), ["f"] = x * -2 + y + z * 5 }, qp_options);
      try {
        res = solver.call(new Dictionary<string, DM>());
      } catch (SwigException) {
        Console.WriteLine(str(solver.stats()));
      }
      checkarrayDigits(res["x"], ca.vertcat(-5.0147928622, -5.766930599, -8.52180472), 4, c.name);
      assertTrue((bool)solver.stats()["success"]);

      MX x4 = MX.sym("x", 4);
      DM H = ca.DM(new double[,] {
        { 2.834044009405148, 0.867080384259271, 0.396881849048015, 0.506784822363357 },
        { 0.867080384259271, 2.184677189537596, 0.725076945381028, 1.223678163433993 },
        { 0.396881849048015, 0.725076945381028, 2.838389028806589, 0.712607093594686 },
        { 0.506784822363357, 1.223678163433993, 0.712607093594686, 3.340935020356804 } });
      DM x0 = ca.vertcat(1, 2, 3, 4);
      MX f = -ca.mtimes(ca.mtimes(x4.T(), (MX)H), (MX)x0);
      var (D, Lt, _) = ca.ldl(H);
      DM Fm = ca.mtimes(ca.sqrt(ca.diag(D)), DM.eye(4) + Lt);
      h = ca.soc(ca.vertcat(ca.mtimes((MX)Fm, x4) * Math.Sqrt(2), -y + 1), y + 1);
      solver = ca.qpsol("msyolver", c.name, new Dictionary<string, MX> { ["x"] = ca.vertcat(x4, y), ["f"] = y + f, ["h"] = h }, qp_options);
      res = solver.call(new Dictionary<string, DM> { ["lbx"] = ca.vertcat(-inf, -inf, -inf, -inf, 1) });
      checkarrayDigits(res["x"][0..^1], x0, 4, c.name);
    }
  }

  public void test_no_success() {
    SX x = SX.sym("x");
    SX y = SX.sym("y");
    SX f = x - y;
    foreach (Conic c in conics) {
      Console.WriteLine("test_no_success " + c);
      var opts = c.opts();
      opts["error_on_fail"] = false;
      var qp = new Dictionary<string, SX> { ["x"] = ca.vertcat(x, y), ["f"] = f, ["g"] = ca.vertcat(x + 1, x - 2) };
      var args = new Dictionary<string, DM> { ["x0"] = 0, ["lbg"] = 0, ["ubg"] = 0, ["lbx"] = ca.vertcat(-10, -10), ["ubx"] = ca.vertcat(10, 10) };
      Function solver = ca.qpsol("solver", c.name, qp, opts);
      solver.call(args);
      assertFalse((bool)solver.stats()["success"], c.name);

      opts["error_on_fail"] = true;
      Function solver2 = ca.qpsol("solver", c.name, qp, opts);
      assertInException("process", () => solver2.call(args));
    }
  }

  public void test_expand() {
    Check.Assume(hasConic("osqp"), "no conic osqp");
    MX x = MX.sym("x");
    var qp = new Dictionary<string, MX> { ["x"] = x, ["f"] = x - 1 };
    Function solver = ca.qpsol("solver", "osqp", qp);
    Function solver_qp = solver.find_function("solver_qp");
    assertTrue(solver_qp.is_a("MXFunction"));
    solver = ca.qpsol("solver", "osqp", qp, map("expand", true));
    solver_qp = solver.find_function("solver_qp");
    assertTrue(solver_qp.is_a("SXFunction"));
  }

  public void test_ipqp_status_init() {
    Check.Assume(hasConic("ipqp"), "no conic ipqp");
    DM H = ca.DM(new double[,] { { 1, -1 }, { -1, 2 } });
    DM A = ca.DM(new double[,] { { 1, 1 }, { -1, 2 }, { 2, 1 } });
    Function solver = ca.conic("solver", "ipqp", hA(H, A), map("print_iter", false, "print_header", false, "print_info", false));
    var res = solver.call(new Dictionary<string, DM> {
      ["h"] = H, ["g"] = ca.vertcat(-2, -6), ["a"] = A, ["lbx"] = 0, ["ubx"] = inf, ["lba"] = -inf, ["uba"] = ca.vertcat(2, 2, 3) });
    assertTrue((bool)solver.stats()["success"]);
    checkarrayDigits(res["x"], ca.vertcat(2.0 / 3, 4.0 / 3), 6, "");
  }

  public void test_qrqp_iter_count() {
    Check.Assume(hasConic("qrqp"), "no conic qrqp");
    DM H = ca.DM(new double[,] { { 1, -1 }, { -1, 2 } });
    DM A = ca.DM(new double[,] { { 1, 1 }, { -1, 2 }, { 2, 1 } });
    Function solver = ca.conic("solver", "qrqp", hA(H, A), map("print_iter", false, "print_header", false, "print_info", false));
    solver.call(new Dictionary<string, DM> {
      ["h"] = H, ["g"] = ca.vertcat(-2, -6), ["a"] = A, ["lbx"] = 0, ["ubx"] = inf, ["lba"] = -inf, ["uba"] = ca.vertcat(2, 2, 3) });
    assertTrue((bool)solver.stats()["success"]);
    assertTrue((long)solver.stats()["iter_count"] > 0);
  }

  public void test_proxqp_mixed_multipliers() {
    Check.Assume(hasConic("proxqp"), "no conic proxqp");
    DM H = DM.eye(3);
    DM A = ca.DM(new double[,] { { 1, 1, 0 }, { 0, 1, -1 }, { 1, 0, 1 } });
    foreach (string backend in new[] { "sparse", "dense" }) {
      Function solver = ca.conic("solver", "proxqp", hA(H, A), map("proxqp", map("eps_abs", 1e-11, "backend", backend)));
      var res = solver.call(new Dictionary<string, DM> {
        ["h"] = H, ["g"] = ca.vertcat(0, -2, 1), ["a"] = A, ["lba"] = ca.vertcat(-inf, -0.5, 1.8), ["uba"] = ca.vertcat(1.5, -0.5, inf) });
      checkarrayDigits(res["x"], ca.vertcat(0.7, 0.6, 1.1), 6, backend);
      checkarrayDigits(res["lam_a"], ca.vertcat(0, 1.4, -0.7), 6, backend);
    }
  }

  public void test_sensitivity() {
    Check.Assume(hasConic("qrqp"), "no conic qrqp");
    MX x = MX.sym("x", 4);
    MX p = MX.sym("p", 2);
    MX p0 = p[0], p1 = p[1];
    MX H = ca.blockcat(new[] {
        new MX[] { p0 + 2, 0.5, 0, 0 },
        new MX[] { 0.5, 3, p1, 0 },
        new MX[] { 0, p1, 4, 0 },
        new MX[] { 0, 0, 0, 1 } });
    MX f = ca.bilin(H, x, x) * 0.5 + ca.dot(ca.vertcat(p0 * -5, 2, p1 - 5, 1), x);
    MX g = ca.vertcat(x[0] + p1 * x[1] + x[2], x[0] - x[2] + x[3], x[1] + p0 * x[2]);
    var qp = new Dictionary<string, MX> { ["x"] = x, ["p"] = p, ["f"] = f, ["g"] = g };
    var bounds = new Dictionary<string, MX> {
      ["lbg"] = ca.vertcat(-inf, p1 - 0.2, p0 + 0.3),
      ["ubg"] = ca.vertcat(p0 + 1, p1 - 0.2, inf),
      ["lbx"] = ca.vertcat(-inf, p1 * 1.5, -inf, -inf),
      ["ubx"] = inf };
    DM pv = ca.vertcat(0.7, 0.3);

    Function @ref = ca.nlpsol("ref", "sqpmethod", qp, map("qpsol", "qrqp", "qpsol_options", map("print_iter", false, "print_header", false),
        "print_iteration", false, "print_header", false, "print_status", false, "tol_pr", 1e-14, "tol_du", 1e-14));
    var args = new Dictionary<string, MX>(bounds) { ["p"] = p };
    var r = @ref.call(args);
    var R = ca.Function("R", new[] { p }, new[] { r["x"], r["f"], r["lam_g"], r["lam_x"] });
    DM lam_g = R.call(dms(pv))[2];
    assertTrue(at(lam_g, 0) > 1e-3 && at(lam_g, 2) < -1e-3 && Math.Abs(at(lam_g, 1)) > 1e-3);
    assertTrue(at(R.call(dms(pv))[3], 1) < -1e-3);

    foreach (Conic c in conics) {
      if (!c.has("quadratic") || !c.has("dual")) continue;
      if (c.name == "nlpsol") continue;
      Console.WriteLine("test_sensitivity " + c);
      int less_digits = c.lessDigits();
      int digits = new[] { "qrqp", "qpoases", "daqp", "ipqp" }.Contains(c.name) ? 7 : 4;
      var o = c.opts();
      o["sens_min_lam"] = 1e-6;
      Function solver = ca.qpsol("solver", c.name, qp, o);
      var rs = solver.call(args);
      var F = ca.Function("F", new[] { p }, new[] { rs["x"], rs["f"], rs["lam_g"], rs["lam_x"] });
      checkfunction(F, R, dms(pv), 0.5 * Math.Pow(10, -(digits - less_digits)));
    }
  }

  public void test_sensitivity_equality_zero_multiplier() {
    MX x = MX.sym("x", 3);
    MX p = MX.sym("p");
    var qp = new Dictionary<string, MX> { ["x"] = x, ["p"] = p, ["f"] = ca.sumsqr(x - p), ["g"] = x[1] - x[0] };
    var R = ca.Function("R", new[] { p }, new[] { ca.vertcat(p, p, p) });
    foreach (Conic c in conics) {
      if (!c.has("quadratic") || !c.has("dual")) continue;
      if (c.name == "nlpsol" || c.name == "superscs") continue;
      Console.WriteLine("test_sensitivity_equality_zero_multiplier " + c);
      int digits = new[] { "qrqp", "qpoases", "daqp", "ipqp" }.Contains(c.name) ? 7 : 4;
      var o = c.opts();
      o["sens_min_lam"] = 1e-6;
      o["sens_force_equality"] = true;
      Function solver = ca.qpsol("solver", c.name, qp, o);
      var args = new Dictionary<string, MX> {
        ["p"] = p, ["lbg"] = 0, ["ubg"] = 0, ["lbx"] = ca.vertcat(-inf, -inf, p), ["ubx"] = ca.vertcat(inf, inf, p) };
      var r = solver.call(args);
      var F = ca.Function("F", new[] { p }, new[] { r["x"] });
      checkfunction(F, R, dms(0.5), 0.5 * Math.Pow(10, -(digits - c.lessDigits())));
    }
  }

  public void test_postpone_expand() {
    Check.Assume(hasConic("osqp"), "no conic osqp");
    MX x = MX.sym("x");
    var J = ca.Function("jac_f", new[] { x, ca.MX(1, 1) }, new[] { x * 0.5 }, new[] { "x", "out_r" }, new[] { "jac_r_x" },
      map("always_inline", true));
    var f = ca.Function("f", new[] { x }, new[] { ca.pow(x, 2) }, new[] { "x" }, new[] { "r" }, map("custom_jacobian", J, "jac_penalty", 0));
    var qp = new Dictionary<string, MX> { ["x"] = x, ["f"] = f.Call(x - 1) };

    var cases = new (Dictionary<string, object>, double, string)[] {
      (map("print_time", true), 0.5, "MXFunction"),
      (map("expand", true, "print_time", true), 2.0, "SXFunction"),
      (map("print_time", true, "postpone_expand", true), 0.5, "MXFunction"),
      (map("expand", true, "print_time", true, "postpone_expand", true), 0.5, "SXFunction"),
    };
    foreach (var (o, h, cls) in cases) {
      Function solver = ca.qpsol("solver", "osqp", qp, o);
      solver.call(new Dictionary<string, DM> { ["x0"] = 2 });
      Function solver_qp = solver.find_function("solver_qp");
      checkarray(solver_qp.call(new Dictionary<string, DM> { ["x"] = 2 })["H"], h);
      assertTrue(solver_qp.is_a(cls));
    }
  }

  // XpressSpecificTests (class-level @requires_conic("xpress")), folded into this class

  public void test_mip_warmstart() {
    Check.Assume(hasConic("xpress"), "no conic xpress");
    MX x = MX.sym("x", 3);
    string logfile = Path.GetTempFileName();
    File.Delete(logfile);
    try {
      Function solver = ca.qpsol("s", "xpress", new Dictionary<string, MX> { ["x"] = x, ["f"] = -ca.sum1(x), ["g"] = x[0] + x[1] + x[2] },
          map("discrete", new[] { true, true, true }, "mip_start", true, "log_file", logfile, "xpress", map("OUTPUTLOG", 1)));
      var sol = solver.call(new Dictionary<string, DM> {
        ["x0"] = ca.vertcat(1, 1, 0), ["lbx"] = ca.vertcat(0, 0, 0), ["ubx"] = ca.vertcat(1, 1, 1), ["lbg"] = -inf, ["ubg"] = 2 });
      assertTrue((bool)solver.stats()["success"], "MIP with mip_start should succeed");
      assertAlmostEqual((double)sol["f"], -2.0, 6, "MILP optimal cost should be -2");
      string log = File.ReadAllText(logfile).Replace("\r\n", "\n");
      assertTrue(log.Contains("User solution"), "Xpress log must contain 'User solution' when mip_start=True");
      solver.Dispose();
    } finally {
      File.Delete(logfile);
    }
  }

  public void test_iis_infeasible_lp() {
    Check.Assume(hasConic("xpress"), "no conic xpress");
    Check.Assume(!isDarwin, "Xpress XPRSiisfirst corrupts the heap on macOS");
    MX x = MX.sym("x");
    Function solver = ca.qpsol("s", "xpress", new Dictionary<string, MX> { ["x"] = x, ["f"] = x, ["g"] = ca.vertcat(x, -x) },
        map("error_on_fail", false, "compute_iis", true, "xpress", map("OUTPUTLOG", 0)));
    var sol = solver.call(new Dictionary<string, DM> {
      ["x0"] = 0, ["lbx"] = -1e10, ["ubx"] = 1e10, ["lbg"] = ca.vertcat(5, 5), ["ubg"] = ca.vertcat(1e10, 1e10) });
    var stats = solver.stats();
    assertFalse((bool)stats["success"], "Solver should report failure for infeasible LP");
    assertTrue(stats.ContainsKey("iis_rows"), "get_stats() must contain 'iis_rows' for infeasible LP");
    assertTrue(stats.ContainsKey("iis_cols"), "get_stats() must contain 'iis_cols' for infeasible LP");
    assertTrue(((System.Collections.ICollection)stats["iis_rows"]).Count > 0, "IIS must identify at least one infeasible row");
    foreach (double v in sol["x"].ToArray()) assertTrue(double.IsNaN(v), "x must be NaN when no primal solution is available");
  }

  public void test_iis_infeasible_milp() {
    Check.Assume(hasConic("xpress"), "no conic xpress");
    Check.Assume(!isDarwin, "Xpress XPRSiisfirst corrupts the heap on macOS");
    MX x = MX.sym("x");
    Function solver = ca.qpsol("s", "xpress", new Dictionary<string, MX> { ["x"] = x, ["f"] = x, ["g"] = x },
        map("discrete", new[] { true }, "error_on_fail", false, "compute_iis", true, "xpress", map("OUTPUTLOG", 0)));
    solver.call(new Dictionary<string, DM> { ["x0"] = 0, ["lbx"] = ca.vertcat(0), ["ubx"] = ca.vertcat(1), ["lbg"] = 2, ["ubg"] = 1e10 });
    var stats = solver.stats();
    assertFalse((bool)stats["success"]);
    assertTrue(stats.ContainsKey("iis_rows"));
    assertTrue(stats.ContainsKey("iis_cols"));
    assertTrue(((System.Collections.ICollection)stats["iis_rows"]).Count + ((System.Collections.ICollection)stats["iis_cols"]).Count > 0,
      "IIS must contain at least one row or column");
  }
}

}

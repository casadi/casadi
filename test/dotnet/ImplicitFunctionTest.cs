using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;

namespace CasADi.Tests {

/// <summary>Port of test/java ImplicitFunctionTest (test/python/implicitfunction.py)</summary>
public class ImplicitFunctionTest : CasadiTestCase {

  static readonly bool KNOWN_BUGS = Environment.GetEnvironmentVariable("CASADI_DOTNET_KNOWN_BUGS") != null;

  sealed class Solver {
    public readonly string name;
    public readonly Dictionary<string, object> options;
    public readonly bool codegen;

    public Solver(string name, Dictionary<string, object> options, bool codegen) {
      this.name = name;
      this.options = options;
      this.codegen = codegen;
    }
  }

  static readonly List<Solver> solvers = makeSolvers();

  static List<Solver> makeSolvers() {
    var r = new List<Solver>();
    if (ca.has_rootfinder("kinsol")) r.Add(new Solver("kinsol", opts("abstol", 1e-10), false));
    if (ca.has_nlpsol("ipopt"))
      r.Add(new Solver("nlpsol", opts("nlpsol", "ipopt", "nlpsol_options",
          opts("print_time", false, "ipopt", opts("print_level", 0))), false));
    r.Add(new Solver("newton", opts(), false));
    r.Add(new Solver("fast_newton", opts(), true));
    return r;
  }

  // ---------------------------------------------------------------- helpers

  static Dictionary<string, object> opts(params object[] kv) {
    var m = new Dictionary<string, object>();
    for (int i = 0; i < kv.Length; i += 2) m[(string)kv[i]] = kv[i + 1];
    return m;
  }

  static Dictionary<string, object> with(IDictionary<string, object> b, params object[] kv) {
    var m = new Dictionary<string, object>(b);
    foreach (var e in opts(kv)) m[e.Key] = e.Value;
    return m;
  }

  static List<DM> ins(params double[] v) => v.Select(e => ca.DM(e)).ToList();

  static SX vs(params SX[] a) => ca.vertcat(a);

  static double num(object o) => Convert.ToDouble(o);

  /// <summary>Python casadiTestCase.assertAlmostEqual: places after scaling large values</summary>
  static void almost(double first, double second, int places, string msg) {
    if (double.IsNaN(first) && double.IsNaN(second)) return;
    if (first == second) return;
    double n = Math.Max(Math.Abs(first), Math.Abs(second));
    n = n > 1e3 ? Math.Pow(10, Math.Floor(Math.Log10(n))) : 1.0;
    double diff = Math.Abs(first / n - second / n);
    assertTrue(diff < 0.5 * Math.Pow(10, -places), msg + " " + first + " <-> " + second + " [places=" + places + "]");
  }

  /// <summary>Python casadiTestCase.checkarray(zr, zt, digits)</summary>
  static void chk(DM zr, DM zt, int digits, string msg = "") {
    assertEquals(zr.size1(), zt.size1(), msg + ": size1 " + zr + " <-> " + zt);
    assertEquals(zr.size2(), zt.size2(), msg + ": size2 " + zr + " <-> " + zt);
    double[,] a = zr.ToArray(), b = zt.ToArray();
    for (int i = 0; i < a.GetLength(0); i++)
      for (int j = 0; j < a.GetLength(1); j++)
        almost(b[i, j], a[i, j], digits, msg + " at (" + i + "," + j + "): " + zt + " <-> " + zr);
  }

  static void chk(DM zr, double zt, int digits) => chk(zr, ca.DM(zt), digits, "");

  static Function jacobianOld(Function f, long i, long j) {
    var o = new List<string> { "jac:" + f.name_out(j) + ":" + f.name_in(i) };
    o.AddRange(f.name_out());
    return f.factory(f.name() + "_jac", f.name_in(), o);
  }

  static Function hessianOld(Function f, long i, long j) {
    var o = new List<string> {
      "hess:" + f.name_out(j) + ":" + f.name_in(i) + ":" + f.name_in(i),
      "grad:" + f.name_out(j) + ":" + f.name_in(i) };
    o.AddRange(f.name_out());
    return f.factory(f.name() + "_hess", f.name_in(), o);
  }

  /// <summary>checkfunction of helpers.py with jacobian/hessian on; evals as given, digits_sens = digits</summary>
  void checkfun(Function trial, Function solution, List<DM> inputs, int digits, int evals, string fail) =>
    checkfun(trial, solution, inputs, digits, true, true, evals, fail);

  void checkfun(Function trial, Function solution, List<DM> inputs, int digits, bool jacobian, bool hessian, int evals,
                string fail) {
    for (int rep = 0; rep < 2; rep++) {
      List<DM> to = trial.call(inputs), so = solution.call(inputs);
      assertEquals(trial.n_out(), solution.n_out(), fail + ": n_out");
      assertEquals(trial.n_in(), solution.n_in(), fail + ": n_in");
      for (int i = 0; i < trial.n_out(); i++) {
        if (trial.sparsity_out(i).is_empty() || solution.sparsity_out(i).is_empty()) continue;
        chk(to[i], so[i], digits, fail + ": output(" + i + ": '" + trial.name_out(i) + "')");
      }
    }
    if (jacobian)
      for (long i = 0; i < trial.n_in(); i++) {
        if (trial.sparsity_in(i).is_empty() || solution.sparsity_in(i).is_empty()) continue;
        for (long j = 0; j < trial.n_out(); j++) {
          Function tj = jacobianOld(trial, i, j), sj = jacobianOld(solution, i, j);
          assertEquals(trial.n_in(), tj.n_in());
          assertEquals(trial.n_out() + 1, tj.n_out());
          assertEquals(solution.n_in(), sj.n_in());
          assertEquals(solution.n_out() + 1, sj.n_out());
          checkfun(tj, sj, inputs, digits, false, false, 0, "(" + fail + ").jacobian_old(" + i + "," + j + ")");
        }
      }
    if (hessian)
      for (long i = 0; i < trial.n_in(); i++) {
        if (trial.sparsity_in(i).is_empty() || solution.sparsity_in(i).is_empty()) continue;
        for (long j = 0; j < trial.n_out(); j++) {
          if (!(trial.sparsity_out(j).is_scalar() && solution.sparsity_out(j).is_scalar())) continue;
          Function th = hessianOld(trial, i, j), sh = hessianOld(solution, i, j);
          assertEquals(trial.n_in(), th.n_in());
          assertEquals(trial.n_out() + 2, th.n_out());
          assertEquals(solution.n_in(), sh.n_in());
          assertEquals(solution.n_out() + 2, sh.n_out());
          checkfun(th, sh, inputs, digits, false, false, 0, "(" + fail + ").hessian_old(" + i + "," + j + ")");
        }
      }
    if (evals > 0) checkEvals(trial, solution, inputs, evals, digits, fail);
  }

  void checkfun(Function trial, Function solution, List<DM> inputs, int digits) =>
    checkfun(trial, solution, inputs, digits, 2, "");

  static Sparsity spmod(int k, Sparsity sp) {
    DM ret = ca.DM(sp);
    if (k == 1 && ret.numel() > 0) ret[0, 0] = ca.DM(1, 1);
    if (k == 2 && ret.nnz() > 0) {
      List<long> row = ret.sparsity().row(), col = ret.sparsity().get_col();
      ret[row[^1], col[^1]] = ca.DM(1, 1);
    }
    return ret.sparsity();
  }

  static bool[] cat(params IEnumerable<bool>[] a) => a.SelectMany(e => e).ToArray();

  static List<DM> seeded(List<DM> inputs, Function vf, long from) {
    var r = new List<DM>(inputs);
    for (long i = from; i < vf.n_in(); i++) {
      var rnd = new Random((int)i);
      double[] v = new double[vf.nnz_in(i)];
      for (int k = 0; k < v.Length; k++) v[k] = rnd.NextDouble();
      r.Add(ca.DM(vf.sparsity_in(i), v.Length == 0 ? ca.DM(0, 1) : (DM)v));
    }
    return r;
  }

  static Function seedFunction(Function f, List<Sparsity> fsp, List<Sparsity> asp, List<Sparsity> isp) {
    var fseeds = fsp.Select(s => MX.sym("f", s)).ToList();
    var aseeds = asp.Select(s => MX.sym("a", s)).ToList();
    var inp = isp.Select(s => MX.sym("i", s)).ToList();
    List<MX> res = f.call(inp, !f.is_a("SXFunction"));
    var o = opts("helper_options", opts("is_diff_in", f.is_diff_in().ToArray(), "is_diff_out", f.is_diff_out().ToArray()));
    List<MX> fwd = ca.forward(res, inp, new[] { fseeds }, o)[0];
    List<MX> adj = ca.reverse(res, inp, new[] { aseeds }, o)[0];
    var vin = inp.Concat(fseeds).Concat(aseeds).ToList();
    var vout = res.Concat(fwd).Concat(adj).ToList();
    return ca.Function("vf", vin, vout, opts("is_diff_in", cat(f.is_diff_in(), f.is_diff_in(), f.is_diff_out()),
        "is_diff_out", cat(f.is_diff_out(), f.is_diff_out(), f.is_diff_in())));
  }

  void checkEvals(Function trial, Function solution, List<DM> inputs, int evals, int digits, string fail) {
    var storage = new Dictionary<string, List<List<DM>>>();
    var storage2 = new Dictionary<string, List<List<DM>>>();
    Function vfRef = null;
    foreach (Function f in new[] { solution, trial }) {
      for (int m1 = 0; m1 < 3; m1++)
        for (int m2 = 0; m2 < 3; m2++) {
          var fsp = new List<Sparsity>();
          var asp = new List<Sparsity>();
          var isp = new List<Sparsity>();
          for (long i = 0; i < f.n_in(); i++) {
            fsp.Add(spmod(m1, f.sparsity_in(i)));
            isp.Add(f.sparsity_in(i));
          }
          for (long i = 0; i < f.n_out(); i++) asp.Add(spmod(m2, f.sparsity_out(i)));
          Function vf = seedFunction(f, fsp, asp, isp);
          string key = m1 + "," + m2;
          if (!storage.ContainsKey(key)) storage[key] = new List<List<DM>>();
          storage[key].Add(vf.call(seeded(inputs, vf, f.n_in())));
          if (vfRef == null) vfRef = vf;
          if (evals > 1) {
            // python repeats this identically for 9 unused (spmod_2, spmod2_2) pairs
            var f2 = new List<Sparsity>();
            var a2 = new List<Sparsity>();
            var i2 = new List<Sparsity>();
            for (long i = 0; i < vf.n_in(); i++) {
              f2.Add(vfRef.sparsity_in(i));
              i2.Add(vfRef.sparsity_in(i));
            }
            for (long i = 0; i < vf.n_out(); i++) a2.Add(vfRef.sparsity_out(i));
            Function vf2 = seedFunction(vf, f2, a2, i2);
            if (!storage2.ContainsKey(key)) storage2[key] = new List<List<DM>>();
            storage2[key].Add(vf2.call(seeded(inputs, vf2, f.n_in())));
          }
        }
    }
    var stores = evals > 1 ? new[] { storage, storage2 } : new[] { storage };
    for (int s = 0; s < stores.Length; s++)
      foreach (List<List<DM>> st in stores[s].Values)
        for (int i = 0; i < st.Count - 1; i++)
          for (int k = 0; k < st[0].Count; k++) {
            DM a = st[0][k], b = st[i + 1][k];
            if (b.numel() == 0 && ca.sparsify(a).nnz() == 0) continue;
            if (a.numel() == 0 && ca.sparsify(b).nnz() == 0) continue;
            chk(a, b, digits, (s == 0 ? "first-order" : "second-order") + ", output(" + k + ")" + fail);
          }
  }

  /// <summary>check_serialize of helpers.py, without the python-side serialization scheme validation</summary>
  static void checkSerialize(Function f, List<DM> inputs) {
    Function f2 = Function.deserialize(f.serialize(opts("debug", true)));
    List<DM> a = f.call(inputs), b = f2.call(inputs);
    for (int i = 0; i < f.n_out(); i++) chk(a[i], b[i], 16);
  }

  static int cgCounter = 0;

  /// <summary>check_codegen of helpers.py: compiles with gcc only when slow tests run</summary>
  static void checkCodegen(Function f, List<DM> inputs) {
    if (!RunSlow) return;
    string dir = Path.Combine(Path.GetTempPath(), "casadi_dotnet_cg_" + Guid.NewGuid().ToString("N"));
    Directory.CreateDirectory(dir);
    string name = "codegen_" + DateTime.Now.Ticks.ToString("x") + "_" + (cgCounter++);
    var cg = ca.CodeGenerator(name);
    cg.add(f);
    cg.generate(dir + Path.DirectorySeparatorChar);
    string inc = GlobalOptions.getCasadiIncludePath(), lib = GlobalOptions.getCasadiPath();
    string so = Path.Combine(dir, name + ".so");
    var psi = new ProcessStartInfo("gcc") { UseShellExecute = false };
    foreach (string a in new[] { "-pedantic", "-std=c89", "-fPIC", "-shared", "-Wall", "-Werror", "-Wextra",
                                 "-I" + inc, "-I" + inc + "/highs", "-Wno-dangling-pointer",
                                 "-Wno-unknown-warning-option", "-Wno-unknown-pragmas", "-Wno-long-long",
                                 "-Wno-unused-parameter", "-O3", Path.Combine(dir, name + ".c"), "-o", so,
                                 "-L" + lib, "-Wl,-rpath," + lib })
      psi.ArgumentList.Add(a);
    using (var p = Process.Start(psi)) {
      p.WaitForExit();
      assertEquals(0, p.ExitCode, "gcc failed");
    }
    Function f2 = ca.external(f.name(), so);
    List<DM> ra = f.call(inputs), rb = f2.call(inputs);
    for (int i = 0; i < f.n_out(); i++) chk(ra[i], rb[i], 15);
  }

  // ---------------------------------------------------------------- tests

  public void test_linear() {
    foreach (Solver S in solvers) {
      message(S.name);
      SX x = SX.sym("x", 2);
      DM A_ = ca.DM(new double[,] { { 1, 2 }, { 3, 2.1 } }), C_ = ca.DM(new double[,] { { 1.6, 2.1 }, { 1, 1.3 } }), b_ = ca.vertcat(0.7, 0.6);
      var f = ca.Function("f", new[] { x }, new[] { ca.mtimes(A_, x) - b_, ca.mtimes(C_, x) });
      Function solver = ca.rootfinder("solver", S.name, f, S.options);
      solver.call(ins(0));
      DM sol = ca.solve(A_, b_);
      var refsol = ca.Function("refsol", new[] { x }, new[] { ca.SX(sol), ca.SX(ca.mtimes(C_, sol)) });
      checkfun(solver, refsol, ins(0), 10);
      if (S.name.Contains("newton")) checkSerialize(solver, ins(0));
      if (S.codegen) checkCodegen(solver, ins(0));

      SX A = SX.sym("A", 2, 2), b = SX.sym("b", 2);
      f = ca.Function("f", new[] { x, A, b }, new[] { ca.mtimes(A, x) - b });
      solver = ca.rootfinder("solver", S.name, f, S.options);
      var solverIn = new List<DM> { 0.0, A_, b_ };
      refsol = ca.Function("refsol", new[] { x, A, b }, new[] { ca.solve(A, b) });
      checkfun(solver, refsol, solverIn, 10);
      if (S.name.Contains("newton")) checkSerialize(solver, solverIn);
      if (S.codegen) checkCodegen(solver, solverIn);

      A = SX.sym("A", 2, 2);
      b = SX.sym("b", 2);
      f = ca.Function("f", new[] { x, A, b }, new[] { ca.mtimes(A, x) - b, ca.mtimes(C_, x) });
      foreach (int adWeightSp in new[] { 0, 1 })
        foreach (int adWeight in new[] { 0, 1 }) {
          message(adWeight + " " + adWeightSp);
          solver = ca.rootfinder("solver", S.name, f, with(S.options, "ad_weight_sp", adWeightSp, "ad_weight", adWeight));
          refsol = ca.Function("refsol", new[] { x, A, b }, new[] { ca.solve(A, b), ca.mtimes(C_, ca.solve(A, b)) });
          checkfun(solver, refsol, solverIn, 10);
          if (S.name.Contains("newton")) checkSerialize(solver, solverIn);
        }
      if (S.codegen) checkCodegen(solver, solverIn);
    }
  }

  public void test_missing_symbols() {
    foreach (Solver S in solvers) {
      message(S.name);
      SX x = SX.sym("x"), p = SX.sym("p");
      var f = ca.Function("f", new[] { x }, new[] { ca.sin(x + p) }, opts("allow_free", true));
      assertInException("[p] are free", () => ca.rootfinder("solver", S.name, f, S.options));
    }
  }

  public void test_scalar1() {
    message("Scalar implicit problem, n=0");
    foreach (Solver S in solvers) {
      message(S.name);
      SX x = SX.sym("x");
      var f = ca.Function("f", new[] { x }, new[] { ca.sin(x) });
      Function solver = ca.rootfinder("solver", S.name, f, S.options);
      solver.call(ins(0));
      var refsol = ca.Function("refsol", new[] { x }, new[] { ca.ceil(x / Math.PI - 0.5) * Math.PI });
      checkfun(solver, refsol, ins(6), 10);
    }
  }

  public void test_unsolved_stats() {
    foreach (Solver S in solvers) {
      message(S.name);
      SX x = SX.sym("x");
      var f = ca.Function("f", new[] { x }, new[] { ca.sin(x) });
      Function solver = ca.rootfinder("solver", S.name, f, S.options);
      try {
        solver.stats();
      } catch (SwigException) {
        // no segfault
      }
    }
  }

  public void test_scalar2() {
    message("Scalar implicit problem, n=1");
    foreach (Solver S in solvers) {
      message(S.name);
      SX x = SX.sym("x"), y = SX.sym("y");
      double n = 0.2;
      var f = ca.Function("f", new[] { y, x }, new[] { x - ca.asin(y) });
      Function solver = ca.rootfinder("solver", S.name, f, S.options);
      var refsol = ca.Function("refsol", new[] { y, x }, new[] { ca.sin(x) });
      checkfun(solver, refsol, ins(0, n), 6, 2, S.name);
    }
  }

  public void test_scalar2_indirect() {
    foreach (Solver S in solvers) {
      message(S.name);
      SX x = SX.sym("x"), y = SX.sym("y");
      double n = 0.2;
      var f = ca.Function("f", new[] { y, x }, new[] { x - ca.asin(y) });
      Function solver = ca.rootfinder("solver", S.name, f, S.options);
      MX X = MX.sym("X");
      MX R = solver.call(new[] { ca.MX(), X })[0];
      var trial = ca.Function("trial", new[] { X }, new[] { R });
      var refsol = ca.Function("refsol", new[] { x }, new[] { ca.sin(x) });
      checkfun(trial, refsol, ins(n), 6, 2, S.name);
    }
  }

  public void test_large() {
    foreach (Solver S in solvers) {
      if (S.name.Contains("kinsol")) continue;
      if (S.name.Contains("newton")) continue;
      int N = 5;
      Sparsity s = Sparsity.lower(N);
      SX x = SX.sym("x", s), y = SX.sym("y", s);
      DM y0 = ca.DM(Sparsity.diag(N), 0.1);
      SX r = (ca.mtimes(x + y0, (x + y0).T()) - ca.mtimes(y + y0, (y + y0).T())).get(false, s);
      var f = ca.Function("f", new[] { y.get_nz(false, ca.Slice()), x.get_nz(false, ca.Slice()) },
          new[] { r.get_nz(false, ca.Slice()) });
      long[] cons = Enumerable.Repeat(1L, (int)s.nnz()).ToArray();
      Function solver = ca.rootfinder("options2", S.name, f, with(S.options, "constraints", cons));
      MX X = MX.sym("X", x.sparsity());
      MX R = solver.call(new[] { ca.MX(), X.get_nz(false, ca.Slice()) })[0];
      var trial = ca.Function("trial", new[] { X }, new[] { R });
      double[] v = Enumerable.Range(0, (int)x.nnz()).Select(i => Math.Abs(Math.Cos(i))).ToArray();
      DM trialIn = ca.DM(trial.sparsity_in(0), v);
      DM trialOut = trial.call(new[] { trialIn })[0];
      f.call(new[] { trialOut, trialIn.get_nz(false, ca.Slice()) });
      f.call(new[] { trialIn.get_nz(false, ca.Slice()), trialIn.get_nz(false, ca.Slice()) });
      var refsol = ca.Function("refsol", new[] { X }, new[] { X.get_nz(false, ca.Slice()) });
      checkfun(trial, refsol, new List<DM> { trialIn[0] }, 6, 1, S.name);
    }
  }

  public void test_vector2() {
    Check.Assume(KNOWN_BUGS, "known bug (set CASADI_DOTNET_KNOWN_BUGS to include it)");
    message("Scalar implicit problem, n=1");
    foreach (Solver S in solvers) {
      message(S.name);
      SX x = SX.sym("x"), y = SX.sym("y", 2);
      DM y0 = ca.vertcat(0.1, 0.4);
      SX yy = y + y0;
      double n = 0.2;
      var f = ca.Function("f", new[] { y, x }, new[] { vs(x - ca.asin(yy[0]), ca.pow(yy[1], 2) - yy[0]) });
      Function solver = ca.rootfinder("solver", S.name, f, S.options);
      var refsol = ca.Function("refsol", new[] { y, x }, new[] { vs(ca.sin(x), ca.sqrt(ca.sin(x))) - y0 });
      checkfun(solver, refsol, ins(n, 0), 4, 2, S.name);
    }
  }

  public void testKINSol1c() {
    message("Scalar KINSol problem, n=0, constraint");
    SX x = SX.sym("x");
    var f = ca.Function("f", new[] { x }, new[] { ca.sin(x) });
    Function solver = ca.rootfinder("solver", "kinsol", f, opts("constraints", new long[] { -1 }));
    DM o = solver.call(ins(-6))[0];
    almost((double)o, -2 * Math.PI, 5, "");
  }

  public void test_kinsol_user_defined_scalar() {
    Check.Assume(hasRootfinder("kinsol"), "no rootfinder kinsol");
    SX x = SX.sym("x");
    var f = ca.Function("f", new[] { x }, new[] { ca.pow(x, 2) - 2 });
    foreach (string strategy in new[] { "none", "linesearch" }) {
      Function solver = ca.rootfinder("solver", "kinsol", f,
          opts("linear_solver_type", "user_defined", "strategy", strategy, "abstol", 1e-10));
      foreach (double x0 in new double[] { 1, 2 }) {
        chk(solver.Call((DM)x0), Math.Sqrt(2), 8);
        assertTrue((bool)solver.stats()["success"]);
      }
    }
  }

  public void test_kinsol_user_defined_vector() {
    Check.Assume(hasRootfinder("kinsol"), "no rootfinder kinsol");
    for (int symbolic = 0; symbolic < 2; symbolic++) {
      Function f;
      if (symbolic == 0) {
        SX x = SX.sym("x", 2), p = SX.sym("p", 2);
        SX g = vs(ca.pow(x[0], 2) + p[0] * x[1] - p[1], x[0] - x[1]);
        f = ca.Function("f", new[] { x, p }, new[] { g });
      } else {
        MX x = MX.sym("x", 2), p = MX.sym("p", 2);
        MX g = ca.vertcat(ca.pow(x[0], 2) + p[0] * x[1] - p[1], x[0] - x[1]);
        f = ca.Function("f", new[] { x, p }, new[] { g });
      }
      foreach (string strategy in new[] { "none", "linesearch" }) {
        var o = opts("strategy", strategy, "abstol", 1e-10, "u_scale", new double[] { 0.5, 2 },
            "f_scale", new double[] { 0.25, 4 });
        Function dense = ca.rootfinder("dense", "kinsol", f, with(o, "linear_solver_type", "dense"));
        Function solver = ca.rootfinder("solver", "kinsol", f, with(o, "linear_solver_type", "user_defined"));
        double[][][] cases = { new[] { new double[] { 1, 6 }, new double[] { 2, 2 } },
                               new[] { new double[] { 2, 15 }, new double[] { 3, 3 } } };
        foreach (double[][] c in cases) {
          DM parameters = ca.vertcat(c[0]);
          DM result = solver.Call(ca.vertcat(1, 1), parameters);
          chk(result, ca.vertcat(c[1]), 8);
          chk(result, dense.Call(ca.vertcat(1, 1), parameters), 8);
          chk(f.Call(result, parameters), DM.zeros(2), 8);
          assertTrue((bool)solver.stats()["success"]);
        }
      }
    }
  }

  public void test_constraints() {
    foreach (Solver S in solvers) {
      if (S.name.Contains("kinsol")) continue;
      if (S.name.Contains("newton")) continue;
      message(S.name + " " + S.options);
      SX x = SX.sym("x", 2);
      SX g = vs(ca.mtimes((x + 3).T(), x - 2), ca.mtimes((x - 4).T(), x + ca.vertcat(1, 2)));
      var f = ca.Function("f", new[] { x }, new[] { g });
      Function solver = ca.rootfinder("solver", S.name, f, with(S.options, "constraints", new long[] { -1, 0 }));
      DM o = solver.call(ins(0))[0];
      double r = Math.Sqrt(1201);
      chk(o, ca.vertcat(-3.0 / 50 * (r - 1), 2.0 / 25 * (r - 1)), 6);
      f = ca.Function("f", new[] { x }, new[] { g });
      solver = ca.rootfinder("solver", S.name, f, with(S.options, "constraints", new long[] { 1, 0 }));
      o = solver.call(ins(0))[0];
      chk(o, ca.vertcat(3.0 / 50 * (r + 1), -2.0 / 25 * (r + 1)), 6);
    }
  }

  public void test_implicitbug() {
    MX X0 = MX.sym("X0"), V = MX.sym("V");
    MX Veq = ca.vertcat(V[0] - X0);
    var vfcn = ca.Function("vfcn", new[] { V, X0 }, new[] { Veq }, opts("ad_weight", 0, "ad_weight_sp", 1));
    Function vfcnSx = vfcn.expand("vfcn_sx", opts("ad_weight", 0, "ad_weight_sp", 1));
    Function ifcn = ca.rootfinder("ifcn", "newton", vfcnSx, opts("linear_solver", "csparse"));
    MX Vs = ifcn.call(new[] { ca.MX(0.0), X0 }, true)[0];
    double f = 1;
    var F = ca.Function("F", new[] { X0 }, new[] { X0 * f + Vs }, opts("ad_weight", 0, "ad_weight_sp", 1));
    Function J = jacobianOld(F, 0, 0);
    List<DM> Jout = J.call(ins(1));
    message(Jout[0].ToString());
    message(J.ToString());
    chk(Jout[0], ca.vertcat(2), 10);
  }

  public void test_extra_outputs() {
    SX x = SX.sym("x"), a = SX.sym("a");
    var f = ca.Function("f", new[] { x, a }, new[] { ca.tan(x) - a, ca.sqrt(a) * ca.pow(x, 2) });
    foreach (Solver S in solvers) {
      message(S.name);
      Function solver = ca.rootfinder("solver", S.name, f, with(S.options, "ad_weight_sp", 1));
      List<DM> solverIn = ins(0.1, 0.3);
      var refsol = ca.Function("refsol", new[] { x, a }, new[] { ca.atan(a), ca.sqrt(a) * ca.pow(ca.atan(a), 2) });
      checkfun(solver, refsol, solverIn, 5);
      if (S.codegen) checkCodegen(solver, solverIn);
    }
  }

  public void test_no_success() {
    SX x = SX.sym("x"), y = SX.sym("y");
    foreach (Solver S in solvers) {
      var rfp = new Dictionary<string, SX> { ["x"] = vs(x, y), ["g"] = vs(ca.sin(x) - 2, ca.sin(y) - 2) };
      Function solver = ca.rootfinder("solver", S.name, rfp, with(S.options, "error_on_fail", false));
      solver.call(new Dictionary<string, DM> { ["x0"] = 0.0 });
      assertFalse((bool)solver.stats()["success"]);
      Function solver2 = ca.rootfinder("solver", S.name, rfp, S.options);
      assertInException("process", () => solver2.call(new Dictionary<string, DM> { ["x0"] = 0.0 }));
    }
  }

  public void test_loop() {
    SX x = SX.sym("x");
    foreach (Solver S in solvers) {
      var rfp = new Dictionary<string, SX> { ["x"] = x, ["g"] = ca.pow(x, 3) - x * 2 + 2 };
      Function solver = ca.rootfinder("solver", S.name, rfp, S.options);
      if (S.name == "kinsol") continue;
      if (S.name == "nlpsol") continue;
      if (S.name == "fast_newton") continue;
      message(S.name);
      message(solver.ToString());
      DM res = solver.call(new Dictionary<string, DM> { ["x0"] = 0.0 })["x"];
      chk(res, -1.7692923542386, 10);
    }
  }

  public void test_bisection() {
    SX x = SX.sym("x"), p = SX.sym("p");
    var f = ca.Function("f", new[] { x, p }, new[] { ca.pow(x, 2) - p });
    Function solver = ca.rootfinder("solver", "bisection", f);
    var refsol = ca.Function("refsol", new[] { x, p }, new[] { ca.sqrt(p) });
    checkfun(solver, refsol, ins(1, 2), 7);
    checkSerialize(solver, ins(1, 2));
    Dictionary<string, object> stats = solver.stats();
    assertTrue((bool)stats["success"]);
    assertEquals("converged_abstol", (string)stats["return_status"]);
    assertTrue(num(stats["search_iter"]) >= 1);
    assertTrue(num(stats["iter_count"]) >= 1);

    chk(solver.Call((DM)(-1), 2), -Math.Sqrt(2), 7);
    refsol = ca.Function("refsol", new[] { x, p }, new[] { -ca.sqrt(p) });
    checkfun(solver, refsol, ins(-1, 2), 7);

    solver = ca.rootfinder("solver", "bisection", f, opts("lb", -10, "ub", 0));
    chk(solver.Call((DM)1, 2), -Math.Sqrt(2), 7);
    solver = ca.rootfinder("solver", "bisection", f, opts("lb", 0, "ub", 10));
    chk(solver.Call((DM)(-1), 2), Math.Sqrt(2), 7);

    solver = ca.rootfinder("solver", "bisection", f, opts("abstol", 1e-12, "abstol_step", 1e-12));
    chk(solver.Call((DM)1, 2), Math.Sqrt(2), 10);

    var rfp = new Dictionary<string, SX> { ["x"] = vs(x, p), ["g"] = vs(x, p) };
    assertInException("scalar", () => ca.rootfinder("solver", "bisection", rfp));

    var f2 = ca.Function("f", new[] { x, p }, new[] { ca.pow(x, 2) + p });
    solver = ca.rootfinder("solver", "bisection", f2, opts("error_on_fail", false, "lb", -1, "ub", 1));
    solver.Call((DM)1, 1);
    assertFalse((bool)solver.stats()["success"]);
    assertEquals("failed_to_bracket_root", (string)solver.stats()["return_status"]);
    Function solver2 = ca.rootfinder("solver", "bisection", f2, opts("lb", -1, "ub", 1));
    assertInException("process", () => solver2.Call((DM)1, 1));
  }

  public void test_segfault_codegen() {
    MX x = MX.sym("x"), y = MX.sym("y");
    double[] u = new double[10], v = new double[10];
    for (int i = 0; i < 10; i++) {
      u[i] = 1 + 4.0 * i / 9;
      v[i] = u[i] * u[i];
    }
    Function interp = ca.interpolant("interp", "bspline", new[] { u }, v, opts());
    MX fx = interp.Call(x);
    var f = ca.Function("f", new[] { x }, new[] { ca.if_else(x < u[0], ca.pow(x, 2),
        ca.if_else(x > u[9], ca.pow(x, 2), fx, false), false) });
    var res = ca.Function("res", new[] { x, y }, new[] { y - f.Call(x) });
    Function rf = ca.rootfinder("rf", "fast_newton", res);
    var g = ca.Function("g", new[] { x }, new[] { x - rf.Call(x, x) });
    checkCodegen(g, ins(1.01));
  }

  public void test_kinsol_return_status() {
    Check.Assume(hasRootfinder("kinsol"), "no rootfinder kinsol");
    SX x = SX.sym("x");
    var rfp = new Dictionary<string, SX> { ["x"] = x, ["g"] = ca.pow(x, 2) - 2 };
    Function G = ca.rootfinder("G", "kinsol", rfp);
    chk(G.call(new Dictionary<string, DM> { ["x0"] = 1.0 })["x"], ca.DM(Math.Sqrt(2)), 6);
    assertEquals("KIN_SUCCESS", (string)G.stats()["return_status"]);
  }

  static Dictionary<string, DM> args(params object[] kv) {
    var m = new Dictionary<string, DM>();
    for (int i = 0; i < kv.Length; i += 2)
      m[(string)kv[i]] = kv[i + 1] is DM d ? d : ca.DM(Convert.ToDouble(kv[i + 1]));
    return m;
  }

  static void assertNames(string[] expected, List<string> actual) =>
    Check.True(expected.SequenceEqual(actual), "expected [" + string.Join(", ", expected) + "], got [" + string.Join(", ", actual) + "]");

  public void test_generic() {
    SX a = SX.sym("a"), b = SX.sym("b"), c = SX.sym("c");
    var f = ca.Function("f", new[] { a, b, c }, new[] { ca.pow(a, 2) - b, ca.pow(b, 2) - c, a * b * c + 1 },
        new[] { "a", "b", "c" }, new[] { "f", "g", "h" });
    double a0 = 1, b0 = 2, c0 = 3;
    Function rf = ca.rootfinder("rf", "newton", f, opts("implicit_input", 0, "implicit_output", 0));
    Dictionary<string, DM> res = rf.call(args("a0", a0, "b", b0, "c", c0));
    almost(Math.Pow((double)res["a"], 2) - b0, 0, 10, "");
    almost(b0 * b0 - c0 - (double)res["g"], 0, 10, "");
    almost((double)res["a"] * b0 * c0 + 1 - (double)res["h"], 0, 10, "");
    message(rf.ToString());
    Function jac = rf.get_function("jac_g_x");
    assertNames(new[] { "a", "b", "c" }, jac.name_in());
    assertNames(new[] { "jac_f_a", "f", "g", "h" }, jac.name_out());

    rf = ca.rootfinder("rf", "newton", f, opts("implicit_input", 1, "implicit_output", 0));
    res = rf.call(args("a", a0, "b0", b0, "c", c0));
    almost(a0 * a0 - (double)res["b"], 0, 10, "");
    almost(Math.Pow((double)res["b"], 2) - c0 - (double)res["g"], 0, 10, "");
    almost(a0 * (double)res["b"] * c0 + 1 - (double)res["h"], 0, 10, "");
    jac = rf.get_function("jac_g_x");
    assertNames(new[] { "a", "b", "c" }, jac.name_in());
    assertNames(new[] { "jac_f_b", "f", "g", "h" }, jac.name_out());

    rf = ca.rootfinder("rf", "newton", f, opts("implicit_input", 2, "implicit_output", 1));
    res = rf.call(args("a", a0, "b", b0, "c0", c0));
    almost(a0 * a0 - b0 - (double)res["f"], 0, 10, "");
    almost(b0 * b0 - (double)res["c"], 0, 10, "");
    almost(a0 * b0 * (double)res["c"] + 1 - (double)res["h"], 0, 10, "");
    jac = rf.get_function("jac_g_x");
    assertNames(new[] { "a", "b", "c" }, jac.name_in());
    assertNames(new[] { "jac_g_c", "f", "g", "h" }, jac.name_out());

    var rfp = new Dictionary<string, SX> { ["x"] = a, ["p"] = vs(b, c), ["g"] = ca.pow(a, 2) - b };
    rf = ca.rootfinder("rf", "newton", rfp);
    res = rf.call(args("x0", a0, "p", ca.vertcat(b0, c0)));
    almost(Math.Pow((double)res["x"], 2) - b0, 0, 10, "");
    jac = rf.get_function("jac_g_x");
    assertNames(new[] { "x", "p" }, jac.name_in());
    assertNames(new[] { "jac_g_x", "g" }, jac.name_out());
  }
}

}

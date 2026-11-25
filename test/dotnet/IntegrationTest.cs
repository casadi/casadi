using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;

namespace CasADi.Tests {

/// <summary>Port of test/java IntegrationTest (test/python/integration.py)</summary>
public class IntegrationTest : CasadiTestCase {

  // ---------------------------------------------------------------- helpers

  sealed class Intg {
    public readonly string name;
    public readonly List<string> features;
    public readonly Dictionary<string, object> options;

    public Intg(string name, string[] features, Dictionary<string, object> options) {
      this.name = name;
      this.features = features.ToList();
      this.options = options;
    }
  }

  sealed class EvIntg {
    public readonly string name;
    public readonly Dictionary<string, object> options;
    public readonly bool events, eventsDae, eventsAdjoint, codegen;
    public readonly int digits;

    public EvIntg(string name, Dictionary<string, object> options, bool events, bool eventsDae, bool eventsAdjoint,
                  bool codegen, int digits) {
      this.name = name;
      this.options = options;
      this.events = events;
      this.eventsDae = eventsDae;
      this.eventsAdjoint = eventsAdjoint;
      this.codegen = codegen;
      this.digits = digits;
    }
  }

  static readonly List<Intg> integrators = makeIntegrators();
  static readonly List<EvIntg> eventIntegrators = makeEventIntegrators();

  static List<Intg> makeIntegrators() {
    var r = new List<Intg>();
    if (ca.has_integrator("cvodes"))
      r.Add(new Intg("cvodes", new[] { "ode" },
          opts("abstol", 1e-15, "reltol", 1e-15, "fsens_err_con", true, "quad_err_con", false)));
    if (ca.has_integrator("idas"))
      r.Add(new Intg("idas", new[] { "dae", "ode" },
          opts("abstol", 1e-15, "reltol", 1e-15, "fsens_err_con", true, "calc_icB", true)));
    r.Add(new Intg("collocation", new[] { "dae", "ode" }, opts("rootfinder", "newton", "number_of_finite_elements", 18)));
    r.Add(new Intg("rk", new[] { "ode" }, opts("number_of_finite_elements", 1000)));
    r.Add(new Intg("collocation", new[] { "dae", "ode" },
        opts("number_of_finite_elements", 18, "simplify", true, "rootfinder", "fast_newton")));
    r.Add(new Intg("rk", new[] { "ode" }, opts("number_of_finite_elements", 1000, "simplify", true)));
    r.Add(new Intg("rosenbrock", new[] { "dae", "ode" }, opts("abstol", 1e-12, "reltol", 1e-12)));
    return r;
  }

  static List<EvIntg> makeEventIntegrators() {
    var r = new List<EvIntg>();
    foreach (string plugin in new[] { "cvodes", "idas" }) {
      if (!ca.has_integrator(plugin)) continue;
      r.Add(new EvIntg(plugin, opts("abstol", 1e-12, "reltol", 1e-12, "quad_err_con", true, "max_event_iter", 20),
          true, false, false, false, 5));
    }
    r.Add(new EvIntg("rk", opts("number_of_finite_elements", 100, "max_event_iter", 20), true, false, false, false, 7));
    r.Add(new EvIntg("collocation", opts("number_of_finite_elements", 20, "max_event_iter", 20), true, false, false, false, 7));
    var rt = new (string, Dictionary<string, object>)[] {
      ("erk", opts()), ("erk", opts("method", "rk4", "step_size", 0.002)),
      ("rosenbrock", opts()), ("rosenbrock", opts("step_size", 0.002)) };
    foreach (var (plugin, o0) in rt) {
      var o = with(o0, "abstol", 1e-11, "reltol", 1e-11, "max_tape_steps", 600);
      r.Add(new EvIntg(plugin, o, true, plugin == "rosenbrock", true, true, 8));
    }
    return r;
  }

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

  static string str<V>(IDictionary<string, V> d) =>
    "{" + string.Join(", ", d.Select(e => e.Key + "=" + e.Value)) + "}";

  static double[] linspace(double a, double b, int n) {
    // numpy.linspace rounding: i*step + start, last point exact
    double[] r = new double[n];
    double step = (b - a) / (n - 1);
    for (int i = 0; i < n; i++) r[i] = i * step + a;
    if (n > 1) r[n - 1] = b;
    return r;
  }

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

  /// <summary>Inputs by name, 0 for the others (as checkfunction with a dict in helpers.py)</summary>
  static List<DM> inputsOf(Function f, IDictionary<string, DM> d) {
    var r = Enumerable.Range(0, (int)f.n_in()).Select(_ => ca.DM(0.0)).ToList();
    List<string> ns = f.name_in();
    foreach (var e in d) r[ns.IndexOf(e.Key)] = e.Value;
    return r;
  }

  /// <summary>checkfunction of helpers.py (gradient/fwd/adj/sens_der flags have no effect there either)</summary>
  void checkfun(Function trial, Function solution, List<DM> inputs, bool jacobian = true, bool hessian = true,
                int evals = 2, int digits = 9, int? digitsSens = null, string fail = "") {
    int ds = digitsSens ?? digits;
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
          checkfun(tj, sj, inputs, jacobian: false, hessian: false, evals: 0, digits: ds,
                   fail: "(" + fail + ").jacobian_old(" + i + "," + j + ")");
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
          checkfun(th, sh, inputs, jacobian: false, hessian: false, evals: 0, digits: ds,
                   fail: "(" + fail + ").hessian_old(" + i + "," + j + ")");
        }
      }
    if (evals > 0) checkEvals(trial, solution, inputs, evals, ds, fail);
  }

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

  void checkfunctionLight(Function trial, Function solution, List<DM> inputs) =>
    checkfun(trial, solution, inputs, jacobian: false, hessian: false, evals: 0);

  /// <summary>check_serialize of helpers.py, without the python-side serialization scheme validation</summary>
  static void checkSerialize(Function f, List<DM> inputs) {
    Function f2 = Function.deserialize(f.serialize(opts("debug", true)));
    List<DM> a = f.call(inputs), b = f2.call(inputs);
    for (int i = 0; i < f.n_out(); i++) chk(a[i], b[i], 16);
  }

  static void checkPure(Function f, List<DM> inputs) {
    List<DM> a = f.call(inputs), b = f.call(inputs);
    for (int i = 0; i < f.n_out(); i++) chk(a[i], b[i], 16);
  }

  void checkThreadSafety(Function f, List<DM> inputs) =>
    checkfunctionLight(f.map(20, "thread", 2), f.map(20, "thread"), inputs);

  static int cgCounter = 0;

  /// <summary>check_codegen of helpers.py: compiles with gcc only when slow tests run</summary>
  static void checkCodegen(Function f, List<DM> inputs, int digits) {
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
    for (int i = 0; i < f.n_out(); i++) chk(ra[i], rb[i], digits);
  }

  static readonly string[] transitionIn = { "index", "t", "x", "z", "p", "u" };
  static readonly string[] transitionOut = { "post_x", "post_z" };

  static Function bounceTransition(long nz = 0) {
    SX x = SX.sym("x", 2), z = SX.sym("z", nz), p = SX.sym("p");
    SX index = SX.sym("index"), t = SX.sym("t"), u = ca.SX(0, 1);
    return ca.Function("transition", new[] { index, t, x, z, p, u }, new[] { ca.vertcat(x[0], -p * x[1]), z },
        transitionIn, transitionOut);
  }

  /// <summary>Central differences of the first output of G at v0</summary>
  static DM fdJacobian(Function G, DM v0, double h = 1e-6) {
    var cols = new List<DM>();
    foreach (DM e in ca.horzsplit(DM.eye(v0.numel()))) {
      DM a = G.call(new[] { v0 + e * h })[0], b = G.call(new[] { v0 - e * h })[0];
      cols.Add((a - b) / (2 * h));
    }
    return ca.horzcat(cols);
  }

  // ---------------------------------------------------------------- tests

  Function integrator;
  Function qe;
  MX qend, q0, par;
  Dictionary<string, SX> dae;
  const double tendN = 2.3, q0N = 7.1, pN = 2;

  public IntegrationTest() {
    SX t = SX.sym("t"), x = SX.sym("x"), p = SX.sym("p");
    dae = new Dictionary<string, SX> { ["t"] = t, ["x"] = x, ["p"] = p, ["ode"] = x / p * ca.pow(t, 2) };
    integrator = ca.integrator("integrator", "cvodes", dae, 0, 2.3,
        opts("reltol", 1e-15, "abstol", 1e-15, "fsens_err_con", true));
    q0 = MX.sym("q0");
    par = MX.sym("p");
    qend = integrator.call(new Dictionary<string, MX> { ["x0"] = q0, ["p"] = par })["xf"];
    qe = ca.Function("qe", new[] { q0, par }, new[] { qend });
  }

  public void test_full() {
    Check.Assume(RunSlow, "slow");
    SX t = SX.sym("t"), q = SX.sym("q"), p = SX.sym("p");
    var dae = new Dictionary<string, SX> { ["t"] = t, ["x"] = q, ["p"] = p, ["ode"] = q / p * ca.pow(t, 2) };
    Function integrator = ca.integrator("integrator", "cvodes", dae, 0, 2.3,
        opts("reltol", 1e-15, "abstol", 1e-15, "fsens_err_con", true));
    assertInException("No stats available", () => integrator.stats());
    double tf = 2.3;
    var solution = ca.Function("solution",
        new Dictionary<string, SX> { ["x0"] = q, ["p"] = p, ["xf"] = q * ca.exp(Math.Pow(tf, 3) / (p * 3)) },
        ca.integrator_in(), ca.integrator_out());
    checkfun(integrator, solution, inputsOf(integrator, new Dictionary<string, DM> { ["x0"] = 0.3, ["p"] = 0.7 }),
             digits: 6);
  }

  public void test_tools_trivial() {
    SX x = SX.sym("x"), p = SX.sym("p", 0), tf = SX.sym("tf");
    var f = ca.Function("f", new[] { x, p }, new[] { x });
    foreach (Function integrator in new[] { ca.simpleRK(f), ca.simpleIRK(f), ca.simpleIntegrator(f) }) {
      var solution = ca.Function("solution", new[] { x, p, tf }, new[] { x * ca.exp(tf) },
          new[] { "x0", "p", "h" }, new[] { "xf" });
      checkfun(integrator, solution, inputsOf(integrator, new Dictionary<string, DM> { ["x0"] = 1, ["h"] = 1 }),
               digits: 3);
    }
  }

  public void test_tools() {
    Check.Assume(RunSlow, "slow");
    SX t = SX.sym("t"), q = SX.sym("q"), p = SX.sym("p"), t0 = SX.sym("t0"), q0 = SX.sym("q0");
    var f = ca.Function("f", new[] { ca.vertcat(q, t), p }, new[] { ca.vertcat(q / p * ca.pow(t, 2), 1) });
    foreach (Function integrator in new[] { ca.simpleRK(f, 500), ca.simpleIRK(f, 50), ca.simpleIntegrator(f) }) {
      SX e = (ca.pow(t0 + t, 3) - ca.pow(t0, 3)) / (p * 3);
      var solution = ca.Function("solver",
          new Dictionary<string, SX> { ["x0"] = ca.vertcat(q0, t0), ["p"] = p, ["h"] = t,
                                       ["xf"] = ca.vertcat(q0 * ca.exp(e), t0 + t) },
          integrator.name_in(), integrator.name_out());
      checkfun(integrator, solution,
               inputsOf(integrator, new Dictionary<string, DM> { ["x0"] = ca.vertcat(0.3, 0), ["p"] = 0.7, ["h"] = 1 }),
               digits: 3);
    }
  }

  /// <summary>One row of the checks() generators</summary>
  sealed class Chk {
    public readonly string feature;
    public readonly Dictionary<string, SX> din, dout, solutionin, solution;
    public readonly Dictionary<string, DM> point;
    public readonly double t0, tf;

    public Chk(string feature, Dictionary<string, SX> din, Dictionary<string, SX> dout, Dictionary<string, SX> solutionin,
               Dictionary<string, SX> solution, Dictionary<string, DM> point, double t0, double tf) {
      this.feature = feature;
      this.din = din;
      this.dout = dout;
      this.solutionin = solutionin;
      this.solution = solution;
      this.point = point;
      this.t0 = t0;
      this.tf = tf;
    }
  }

  /// <summary>variations() of the python tests: an ode also as a dae</summary>
  static List<Chk> variations(Chk c, string zname) {
    var r = new List<Chk> { c };
    if (c.feature == "ode") {
      var din = new Dictionary<string, SX>(c.din);
      var dout = new Dictionary<string, SX>(c.dout);
      SX z = SX.sym(zname, c.din["x"].sparsity().size1(), c.din["x"].sparsity().size2());
      din["z"] = z;
      dout["ode"] = z;
      dout["alg"] = (c.dout["ode"] - z) * -0.8;
      r.Add(new Chk("dae", din, dout, c.solutionin, c.solution, c.point, c.t0, c.tf));
    }
    return r;
  }

  static Function solutionFunction(Chk c, SX tstart, SX tend) {
    foreach (string k in c.solution.Keys.ToList())
      c.solution[k] = ca.substitute(c.solution[k], ca.vertcat(tstart, tend), ca.vertcat(c.t0, c.tf));
    var sol = new Dictionary<string, SX>(c.solutionin);
    foreach (var e in c.solution) sol[e.Key] = e.Value;
    return ca.Function("fs", sol, ca.integrator_in(), ca.integrator_out());
  }

  static Dictionary<string, SX> merged(Dictionary<string, SX> a, Dictionary<string, SX> b) {
    var r = new Dictionary<string, SX>(a);
    foreach (var e in b) r[e.Key] = e.Value;
    return r;
  }

  public void test_jac() {
    message("Test exact jacobian #536");
    DM.set_precision(18);
    try {
      SX tstart = SX.sym("tstart"), tend = SX.sym("tend");
      var integrators = new[] {
        new Intg("idas", new[] { "dae", "ode" }, opts("abstol", 1e-9, "reltol", 1e-9, "fsens_err_con", true,
            "calc_ic", true, "calc_icB", true)),
        new Intg("cvodes", new[] { "ode" }, opts("abstol", 1e-5, "reltol", 1e-5, "fsens_err_con", false,
            "quad_err_con", false)) };
      SX x = SX.sym("x"), t = SX.sym("t");
      var check = new Chk("ode", new Dictionary<string, SX> { ["x"] = x, ["t"] = t },
          new Dictionary<string, SX> { ["ode"] = x }, new Dictionary<string, SX> { ["x0"] = x },
          new Dictionary<string, SX> { ["xf"] = x * ca.exp(tend) }, new Dictionary<string, DM> { ["x0"] = 1 }, 0, 0.9995);
      foreach (Chk c in variations(check, "x")) {
        foreach (Intg I in integrators) {
          message(I.name);
          SX xx = SX.sym("x");
          ca.integrator("dummyIntegrator", I.name, new Dictionary<string, SX> { ["x"] = xx, ["ode"] = xx });
          if (!I.features.Contains(c.feature)) continue;
          var dae = merged(c.din, c.dout);
          Function fs = solutionFunction(c, tstart, tend);
          foreach (string scheme in new[] { "direct", "gmres", "bcgstab" }) {
            message("f_options: newton_scheme " + scheme);
            Function integrator = ca.integrator("integrator", I.name, dae, c.t0, c.tf,
                with(I.options, "newton_scheme", scheme));
            var inp = c.point.Where(e => !fs.sparsity_in(e.Key).is_empty()).ToDictionary(e => e.Key, e => e.Value);
            Dictionary<string, DM> io = integrator.call(inp), fo = fs.call(inp);
            message("res=" + (io["xf"] - fo["xf"]) + " " + fo["xf"]);
            Dictionary<string, object> stats = integrator.stats();
            message(str(stats));
            assertTrue(num(stats["nstepsB"]) < 2500);
            assertTrue(num(stats["nsteps"]) < 1500);
            assertTrue(num(stats["nstepsB"]) < 2500);
            assertTrue(num(stats["nlinsetups"]) < 100);
            assertTrue(num(stats["nlinsetupsB"]) < 250);
          }
        }
      }
    } finally {
      DM.set_precision(6);
    }
  }

  static Chk nonlinCheck(SX tstart, SX tend) {
    double[] A = { 1, 0.1 };
    double p0 = 1.13;
    SX q = SX.sym("y", 2, 1), y0 = q[0], yc0 = q[1], p = SX.sym("p", 1, 1);
    SX arg = ca.atan(yc0 / ca.sqrt(p)) + ca.sqrt(p) * (tend - tstart);
    SX s1 = (y0 * 2 - ca.log(ca.pow(yc0, 2) / p + 1)) / 2 - ca.log(ca.cos(arg));
    SX s2 = ca.sqrt(p) * ca.tan(arg);
    return new Chk("ode", new Dictionary<string, SX> { ["x"] = q, ["p"] = p },
        new Dictionary<string, SX> { ["ode"] = ca.vertcat(q[1], p[0] + ca.pow(q[1], 2)) },
        new Dictionary<string, SX> { ["x0"] = q, ["p"] = p }, new Dictionary<string, SX> { ["xf"] = ca.vertcat(s1, s2) },
        new Dictionary<string, DM> { ["x0"] = A, ["p"] = p0 }, 0, 0.4);
  }

  public void test_lsolvers() {
    message("Test different linear solvers");
    SX tstart = SX.sym("tstart"), tend = SX.sym("tend");
    var integrators = new[] {
      new Intg("idas", new[] { "dae", "ode" }, opts("abstol", 1e-9, "reltol", 1e-9, "fsens_err_con", true,
          "calc_ic", true, "calc_icB", true)),
      new Intg("cvodes", new[] { "ode" }, opts("abstol", 1e-15, "reltol", 1e-15, "fsens_err_con", true,
          "quad_err_con", false)) };
    SX x = SX.sym("x"), p = SX.sym("p"), z = SX.sym("z");
    Chk daeCheck() => new Chk("dae", new Dictionary<string, SX> { ["x"] = x, ["z"] = z },
        new Dictionary<string, SX> { ["alg"] = x - z, ["ode"] = z }, new Dictionary<string, SX> { ["x0"] = x, ["p"] = p },
        new Dictionary<string, SX> { ["xf"] = x * ca.exp(tend - tstart) },
        new Dictionary<string, DM> { ["x0"] = 7.1, ["p"] = 2 }, 0.2, 2.3);
    var checks = new List<Chk> { daeCheck() };
    if (RunSlow) {
      checks.Add(daeCheck());
      checks.Add(new Chk("ode", new Dictionary<string, SX> { ["x"] = x }, new Dictionary<string, SX> { ["ode"] = x },
          new Dictionary<string, SX> { ["x0"] = x, ["p"] = p }, new Dictionary<string, SX> { ["xf"] = x * ca.exp(tend - tstart) },
          new Dictionary<string, DM> { ["x0"] = 7.1, ["p"] = 2 }, 0.2, 2.3));
      checks.Add(nonlinCheck(tstart, tend));
    }
    foreach (Chk c in checks) {
      foreach (Intg I in integrators) {
        message(I.name);
        SX xx = SX.sym("x");
        ca.integrator("dummyIntegrator", I.name, new Dictionary<string, SX> { ["x"] = xx, ["ode"] = xx });
        if (!I.features.Contains(c.feature)) continue;
        var dae = merged(c.din, c.dout);
        Function fs = solutionFunction(c, tstart, tend);
        foreach (string scheme in new[] { "direct", "gmres", "bcgstab" }) {
          string msg = "f_options: newton_scheme " + scheme;
          message(msg);
          Function integrator = ca.integrator("integrator", I.name, dae, c.t0, c.tf,
              with(I.options, "newton_scheme", scheme));
          var inp = c.point.Where(e => !integrator.sparsity_in(e.Key).is_empty()).ToDictionary(e => e.Key, e => e.Value);
          List<DM> inputs = inputsOf(integrator, inp);
          checkfun(integrator, fs, inputs, hessian: false, evals: 0, digits: 4, digitsSens: 4, fail: msg);
          checkPure(integrator, inputs);
          checkSerialize(integrator, inputs);
        }
      }
    }
  }

  public void test_X() {
    message("Extensive integrator tests");
    SX tstart = SX.sym("tstart"), tend = SX.sym("tend");
    foreach (Intg I in integrators) {
      message(I.name);
      double x0 = q0N, p_ = pN;
      SX t = SX.sym("t"), x = SX.sym("x"), p = SX.sym("p"), z = SX.sym("z");
      var checks = new List<Chk>();
      double tf = tendN;
      Dictionary<string, SX> d(params object[] kv) {
        var m = new Dictionary<string, SX>();
        for (int i = 0; i < kv.Length; i += 2) m[(string)kv[i]] = kv[i + 1] is SX s ? s : (SX)Convert.ToDouble(kv[i + 1]);
        return m;
      }
      var rows = new List<(string, Dictionary<string, SX>, Dictionary<string, SX>, Dictionary<string, SX>)> {
        ("ode", d("x", x), d("ode", 0), d("xf", x)),
        ("ode", d("x", x), d("ode", 1), d("xf", x + (tend - tstart))),
        ("ode", d("x", x), d("ode", x), d("xf", x * ca.exp(tend - tstart))),
        ("ode", d("x", x, "t", t), d("ode", t), d("xf", x + (ca.pow(tend, 2) / 2 - ca.pow(tstart, 2) / 2))),
        ("ode", d("x", x, "t", t), d("ode", x * t), d("xf", x * ca.exp(ca.pow(tend, 2) / 2 - ca.pow(tstart, 2) / 2))),
        ("ode", d("x", x, "p", p), d("ode", x / p), d("xf", x * ca.exp((tend - tstart) / p))),
      };
      if (RunSlow)
        rows.AddRange(new[] {
          ("ode", d("x", x), d("ode", x, "quad", 0), d("qf", 0)),
          ("ode", d("x", x), d("ode", x, "quad", 1), d("qf", tend - tstart)),
          ("ode", d("x", x), d("ode", 0, "quad", x), d("qf", x * (tend - tstart))),
          ("ode", d("x", x), d("ode", x, "quad", x), d("qf", x * (ca.exp(tend - tstart) - 1))),
          ("ode", d("x", x, "t", t), d("ode", x, "quad", t), d("qf", ca.pow(tend, 2) / 2 - ca.pow(tstart, 2) / 2)),
          ("ode", d("x", x, "t", t), d("ode", x, "quad", x * t),
              d("qf", x * (ca.exp(tend - tstart) * (tend - 1) - (tstart - 1)))),
          ("ode", d("x", x, "p", p), d("ode", x, "quad", x / p), d("qf", x * (ca.exp(tend - tstart) - 1) / p)),
          ("dae", d("x", x, "z", z), d("ode", z, "alg", (z - x) * -0.8, "quad", z), d("qf", x * (ca.exp(tend - tstart) - 1))),
        });
      foreach (var (feature, din, dout, sol) in rows)
        checks.Add(new Chk(feature, din, dout, d("x0", x, "p", p), sol,
            new Dictionary<string, DM> { ["x0"] = x0, ["p"] = p_ }, 0.2, tf));
      if (RunSlow) checks.Add(nonlinCheck(tstart, tend));
      foreach (Chk tt in checks) {
        foreach (Chk c in variations(tt, "z")) {
          if (!I.features.Contains(c.feature)) continue;
          string msg = I.name + ": " + str(c.din) + " => " + str(c.dout) + ", explicit (" + str(c.solution) + ") tstart = " + c.t0;
          message(msg);
          var dae = merged(c.din, c.dout);
          Function fs = solutionFunction(c, tstart, tend);
          var o = new Dictionary<string, object>(I.options);
          if (I.name == "cvodes" || I.name == "idas") {
            o["abstol"] = 1e-9;
            o["reltol"] = 1e-9;
          }
          if (I.name == "idas") {
            o["init_xdot"] = c.point["x0"].nonzeros();
            o["calc_icB"] = true;
            o["augmented_options"] = opts("init_xdot", null, "abstol", 1e-9, "reltol", 1e-9);
          }
          Function integrator = ca.integrator("integrator", I.name, dae, c.t0, c.tf, o);
          var inp = c.point.Where(e => !integrator.sparsity_in(e.Key).is_empty()).ToDictionary(e => e.Key, e => e.Value);
          List<DM> inputs = inputsOf(integrator, inp);
          checkfun(integrator, fs, inputs, hessian: false, evals: 0, digits: 4, digitsSens: 4, fail: msg);
          if (!str(o).Contains("kinsol")) {
            message(integrator.ToString());
            checkSerialize(integrator, inputs);
          }
        }
      }
    }
  }

  public void test_multipoint() {
    foreach (Intg I in integrators) {
      message(I.name);
      if (I.name == "cvodes" || I.name == "idas") continue;
      MX x = MX.sym("x"), p = MX.sym("p"), z = MX.sym("z"), u = MX.sym("u");
      var dae = I.features.Contains("dae")
          ? new Dictionary<string, MX> { ["x"] = x, ["z"] = z, ["p"] = p, ["u"] = u, ["ode"] = z + p * u,
                                         ["alg"] = z * ca.cos(z) - x, ["quad"] = ca.pow(x, 2) }
          : new Dictionary<string, MX> { ["x"] = x, ["p"] = p, ["u"] = u, ["ode"] = ca.cos(x + p * u), ["quad"] = ca.pow(x, 2) };
      Function integratorA = ca.integrator("integrator", I.name, dae, 0, 0.1, I.options);
      Function integratorB = ca.integrator("integrator", I.name, dae, 0.1, 0.2, I.options);
      Function integratorI = ca.integrator("integrator", I.name, dae, 0, new double[] { 0.1, 0.2 }, I.options);
      var argTest = new Dictionary<string, DM> { ["x0"] = 0, ["p"] = 0.15, ["u"] = ca.vertcat(1, 1.1).T() };
      Dictionary<string, MX> args = integratorI.convert_in(integratorI.mx_in());
      List<MX> uSplit = ca.horzsplit_n(args["u"], 2);
      var args1 = new Dictionary<string, MX>(args) { ["u"] = uSplit[0] };
      var args2 = new Dictionary<string, MX>(args) { ["u"] = uSplit[1] };
      Dictionary<string, MX> res1 = integratorA.call(args1);
      args2["x0"] = res1["xf"];
      Dictionary<string, MX> res2 = integratorB.call(args2);
      res2["qf"] = ca.horzcat(res1["qf"], res1["qf"] + res2["qf"]);
      res2["xf"] = ca.horzcat(res1["xf"], res2["xf"]);
      res2["zf"] = ca.horzcat(res1["zf"], res2["zf"]);
      foreach (var e in args) res2[e.Key] = e.Value;
      var If = ca.Function("I", res2, integratorI.name_in(), integratorI.name_out());
      If.call(argTest);
      integratorI.call(argTest);
      checkfun(integratorI, If,
               inputsOf(integratorI, new Dictionary<string, DM> { ["x0"] = 0, ["p"] = 0.15, ["u"] = ca.vertcat(1, 1.1).T() }));
    }
  }

  public void test_eval2() {
    message("CVodes integration: evaluation with Function indirection");
    var qe = ca.Function("qe", new[] { q0, par }, new[] { qend[0] });
    var f = ca.Function("f", new[] { q0 }, new[] { qe.Call(q0, ca.MX(pN)) });
    DM fOut = f.call(new[] { ca.DM(q0N) })[0];
    almost((double)fOut, q0N * Math.Exp(Math.Pow(tendN, 3) / (3 * pN)), 9, "Evaluation output mismatch");
  }

  public void test_issue92b() {
    message("regression check for issue 92");
    SX t = SX.sym("t"), x = SX.sym("x"), y = SX.sym("y");
    var dae = new Dictionary<string, SX> { ["t"] = t, ["x"] = ca.vertcat(x, y), ["ode"] = ca.vertcat(x, x * (1 + 1e-9)) };
    Function integrator = ca.integrator("integrator", "cvodes", dae, opts("fsens_err_con", true));
    Dictionary<string, DM> o = integrator.call(new Dictionary<string, DM> { ["x0"] = ca.vertcat(1, 0) });
    message(o["xf"].ToString());
  }

  public void test_issue92() {
    message("regression check for issue 92");
    SX t = SX.sym("t"), x = SX.sym("x");
    MX var = MX.sym("var", 2, 1);
    SX q = ca.vertcat(x, SX.sym("problem"));
    SX dq = ca.vertcat(x, x);
    Function integrator = ca.integrator("integrator", "cvodes",
        new Dictionary<string, SX> { ["t"] = t, ["x"] = q, ["ode"] = dq }, opts("fsens_err_con", true, "reltol", 1e-12));
    MX qend = integrator.call(new Dictionary<string, MX> { ["x0"] = var })["xf"];
    var f = ca.Function("f", new[] { var }, new[] { qend[0] });
    Function J = jacobianOld(f, 0, 0);
    List<DM> Jout = J.call(new[] { ca.vertcat(1, 0) });
    message("jac=" + (Jout[0].get_nz(false, ca.Slice(0)) - Math.Exp(1)));
    almost((double)Jout[0][0, 0], Math.Exp(1), 5, "Evaluation output mismatch");
  }

  public void test_eval() {
    message("CVodes integration: evaluation");
    DM o = qe.call(new[] { ca.DM(q0N), ca.DM(pN) })[0];
    almost((double)o, q0N * Math.Exp(Math.Pow(tendN, 3) / (3 * pN)), 9, "Evaluation output mismatch");
  }

  public void test_jac1() {
    message("CVodes integration: jacobian to q0");
    Function J = jacobianOld(qe, 0, 0);
    DM o = J.call(new[] { ca.DM(q0N), ca.DM(pN) })[0];
    almost((double)o, Math.Exp(Math.Pow(tendN, 3) / (3 * pN)), 9, "Evaluation output mismatch");
  }

  public void test_jac2() {
    message("CVodes integration: jacobian to p");
    Function J = jacobianOld(qe, 1, 0);
    DM o = J.call(new[] { ca.DM(q0N), ca.DM(pN) })[0];
    double te = tendN, q0 = q0N, p = pN;
    almost((double)o, -(q0 * Math.Pow(te, 3) * Math.Exp(Math.Pow(te, 3) / (3 * p))) / (3 * p * p), 9,
        "Evaluation output mismatch");
  }

  public void test_bug_repeat() {
    message("Bug that appears when rhs contains repeats");
    DM A = ca.vertcat(1, 0.1, 1);
    double p0 = 1.13, te = 0.4;
    SX t = SX.sym("t"), q = SX.sym("y", 3, 1), p = SX.sym("p");
    SX dh = p + ca.pow(q[0], 2);
    var o = opts("reltol", 1e-15, "abstol", 1e-15, "fsens_err_con", true, "steps_per_checkpoint", 10000);
    Function integrator = ca.integrator("integrator", "cvodes",
        new Dictionary<string, SX> { ["x"] = q, ["p"] = p, ["t"] = t, ["ode"] = ca.vertcat(dh, q[0], dh) }, 0, te, o);
    MX q0 = MX.sym("q0", 3, 1), par = MX.sym("p", 1, 1);
    MX qend = integrator.call(new Dictionary<string, MX> { ["x0"] = q0, ["p"] = par })["xf"];
    var qe = ca.Function("qe", new[] { q0, par }, new[] { qend });
    DM outA = jacobianOld(qe, 0, 0).call(new[] { A, ca.DM(p0) })[0];
    integrator = ca.integrator("integrator", "cvodes",
        new Dictionary<string, SX> { ["x"] = q, ["p"] = p, ["t"] = t, ["ode"] = ca.vertcat(dh, q[0], dh * (1 + 1e-9)) },
        0, te, o);
    q0 = MX.sym("q0", 3, 1);
    par = MX.sym("p", 1, 1);
    qend = integrator.call(new Dictionary<string, MX> { ["x0"] = q0, ["p"] = par })["xf"];
    qe = ca.Function("qe", new[] { q0, par }, new[] { qend });
    DM outB = jacobianOld(qe, 0, 0).call(new[] { A, ca.DM(p0) })[0];
    message((outA - outB).ToString());
  }

  double hessRef() {
    double te = tendN, q0 = q0N, p = pN;
    return (q0 * Math.Pow(te, 6) * Math.Exp(Math.Pow(te, 3) / (3 * p))) / (9 * Math.Pow(p, 4))
        + (2 * q0 * Math.Pow(te, 3) * Math.Exp(Math.Pow(te, 3) / (3 * p))) / (3 * Math.Pow(p, 3));
  }

  public void test_hess3() {
    message("CVodes integration: hessian to p: Jacobian of integrator.jacobian");
    Function J = jacobianOld(integrator, integrator.index_in("p"), integrator.index_out("xf"));
    Function H = jacobianOld(J, J.index_in("p"), 0);
    Dictionary<string, DM> o = H.call(new Dictionary<string, DM> { ["x0"] = q0N, ["p"] = pN });
    almost((double)o["jac_jac_xf_p_p"], hessRef(), 9, "Evaluation output mismatch");
  }

  public void test_hess4() {
    message("CVodes integration: hessian to p: Jacobian of integrator.jacobian indirect");
    Function J = jacobianOld(integrator, integrator.index_in("p"), integrator.index_out("xf"));
    MX q0 = MX.sym("q0"), p = MX.sym("p");
    Dictionary<string, MX> Ji = J.call(new Dictionary<string, MX> { ["x0"] = q0, ["p"] = p });
    Ji["q0"] = q0;
    Ji["p"] = p;
    var Jif = ca.Function("Ji", Ji, new[] { "q0", "p" }, J.name_out());
    Function H = jacobianOld(Jif, 1, 0);
    DM o = H.call(new[] { ca.DM(q0N), ca.DM(pN) })[0];
    almost((double)o, hessRef(), 9, "Evaluation output mismatch");
  }

  public void test_hess5() {
    message("CVodes integration: hessian to p in an MX tree");
    MX q0 = MX.sym("q0"), p = MX.sym("p");
    Dictionary<string, MX> sol = integrator.call(new Dictionary<string, MX> { ["x0"] = q0, ["p"] = p });
    sol["q0"] = q0;
    sol["p"] = p;
    var JT = ca.Function("JT", new[] { q0, p }, new[] { ca.jacobian(sol["xf"], sol["p"]).T() });
    message(string.Join(", ", JT.call(new[] { ca.DM(q0N), ca.DM(pN) })));
    Function H = jacobianOld(JT, 1, 0);
    DM o = H.call(new[] { ca.DM(q0N), ca.DM(pN) })[0];
    almost((double)o, hessRef(), 9, "Evaluation output mismatch");
  }

  public void test_hess6() {
    message("CVodes integration: hessian to p in an MX tree");
    MX q0 = MX.sym("q0"), p = MX.sym("p");
    Dictionary<string, MX> sol = integrator.call(new Dictionary<string, MX> { ["x0"] = q0, ["p"] = p });
    sol["q0"] = q0;
    sol["p"] = p;
    var qe = ca.Function("qe", sol, new[] { "q0", "p" }, ca.integrator_out());
    Function H = hessianOld(qe, 1, 0);
    DM o = H.call(new[] { ca.DM(q0N), ca.DM(pN) })[0];
    almost((double)o, hessRef(), 7, "Evaluation output mismatch");
  }

  static readonly double[,] B3 = { { 1, 2.3, 4 }, { -2, 1.3, 4.7 }, { -2, 6, 9 } };

  public void test_glibcbug() {
    message("former glibc error");
    DM A = ca.vertcat(2.3, 4.3, 7.6);
    double te = 0.7;
    SX t = SX.sym("t"), q = SX.sym("q", 3, 1), p = SX.sym("p", 9, 1);
    Function integrator = ca.integrator("integrator", "cvodes",
        new Dictionary<string, SX> { ["t"] = t, ["x"] = q, ["p"] = p, ["ode"] = ca.mtimes(ca.reshape(p, 3, 3), q) }, 0, te,
        opts("fsens_err_con", true, "steps_per_checkpoint", 1000));
    MX q0 = MX.sym("q0", 3, 1), par = MX.sym("p", 9, 1);
    integrator.call(new Dictionary<string, MX> { ["x0"] = q0, ["p"] = par });
    Function qeJ = jacobianOld(integrator, integrator.index_in("p"), integrator.index_out("xf"));
    MX qe = qeJ.call(new Dictionary<string, MX> { ["x0"] = q0, ["p"] = par })["jac_xf_p"];
    var qef = ca.Function("qef", new[] { q0, par }, new[] { qe });
    // B.ravel(): row-major
    double[] Br = new double[9];
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++) Br[3 * i + j] = B3[i, j];
    qef.call(new[] { A, ca.DM(Br) });
  }

  /// <summary>scipy.linalg.expm replacement: scaling and squaring with a Taylor series</summary>
  static double[,] expm(double[,] M) {
    int n = M.GetLength(0);
    double norm = 0;
    for (int i = 0; i < n; i++) {
      double s = 0;
      for (int j = 0; j < n; j++) s += Math.Abs(M[i, j]);
      norm = Math.Max(norm, s);
    }
    int sq = Math.Max(0, (int)Math.Ceiling(Math.Log(norm) / Math.Log(2)) + 4);
    double sc = Math.Pow(2, -sq);
    double[,] X = new double[n, n], term = new double[n, n], R = new double[n, n];
    for (int i = 0; i < n; i++)
      for (int j = 0; j < n; j++) X[i, j] = M[i, j] * sc;
    for (int i = 0; i < n; i++) {
      term[i, i] = 1;
      R[i, i] = 1;
    }
    for (int k = 1; k < 30; k++) {
      term = matmul(term, X);
      for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++) {
          term[i, j] /= k;
          R[i, j] += term[i, j];
        }
    }
    for (int s = 0; s < sq; s++) R = matmul(R, R);
    return R;
  }

  static double[,] matmul(double[,] a, double[,] b) {
    int n = a.GetLength(0), m = b.GetLength(1), l = b.GetLength(0);
    double[,] r = new double[n, m];
    for (int i = 0; i < n; i++)
      for (int j = 0; j < m; j++)
        for (int k = 0; k < l; k++) r[i, j] += a[i, k] * b[k, j];
    return r;
  }

  public void test_linear_system() {
    message("Linear ODE");
    DM A = ca.vertcat(2.3, 4.3, 7.6);
    double te = 0.7;
    double[,] Bte = new double[3, 3];
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++) Bte[i, j] = B3[i, j] * te;
    DM Be = expm(Bte);
    SX t = SX.sym("t"), q = SX.sym("q", 3, 1), p = SX.sym("p", 9, 1);
    Function integrator = ca.integrator("integrator", "cvodes",
        new Dictionary<string, SX> { ["t"] = t, ["x"] = q, ["p"] = p, ["ode"] = ca.mtimes(ca.reshape(p, 3, 3), q) }, 0, te,
        opts("fsens_err_con", true, "reltol", 1e-15, "abstol", 1e-15, "steps_per_checkpoint", 10000));
    MX q0 = MX.sym("q0", 3, 1), par = MX.sym("p", 9, 1);
    var mxin = new Dictionary<string, MX> { ["x0"] = q0, ["p"] = par };
    MX qend = integrator.call(mxin)["xf"];
    var qe = ca.Function("qe", new[] { q0, par }, new[] { qend });
    MX qendJ = jacobianOld(integrator, integrator.index_in("x0"), integrator.index_out("xf")).call(mxin)["jac_xf_x0"];
    var qeJ = ca.Function("qeJ", new[] { q0, par }, new[] { qendJ });
    MX qendJ2 = jacobianOld(integrator, integrator.index_in("x0"), integrator.index_out("xf")).call(mxin)["jac_xf_x0"];
    var qeJ2 = ca.Function("qeJ2", new[] { q0, par }, new[] { qendJ2 });
    DM B = B3;
    DM qeOut = qe.call(new[] { A, ca.vec(B) })[0];
    chk(ca.mtimes(Be, A) / 1e3, qeOut / 1e3, 10, "jacobian('x0','xf')");
    DM qeJOut = qeJ.call(new[] { A, ca.vec(B) })[0];
    chk(qeJOut / 1e3, Be / 1e3, 10, "jacobian('x0','xf')");
    qeJ2.call(new[] { A, ca.vec(B) });
  }

  public void test_mathieu_system() {
    message("Mathieu ODE");
    DM A = ca.vertcat(0.3, 1.2), B = ca.vertcat(1.3, 4.3, 2.7);
    double te = 0.7;
    SX t = SX.sym("t"), q = SX.sym("q", 2, 1), p = SX.sym("p", 3, 1);
    SX ode = ca.vertcat(q[1], (p[0] - p[1] * 2 * ca.cos(p[2] * 2)) * q[0]);
    Function integrator = ca.integrator("integrator", "cvodes",
        new Dictionary<string, SX> { ["x"] = q, ["p"] = p, ["t"] = t, ["ode"] = ode }, 0, te,
        opts("fsens_err_con", true, "reltol", 1e-15, "abstol", 1e-15, "steps_per_checkpoint", 10000));
    MX q0 = MX.sym("q0", 2, 1), par = MX.sym("p", 3, 1);
    var mxin = new Dictionary<string, MX> { ["x0"] = q0, ["p"] = par };
    MX qend = integrator.call(mxin)["xf"];
    var qe = ca.Function("qe", new[] { q0, par }, new[] { qend });
    MX qendJ = jacobianOld(integrator, integrator.index_in("x0"), integrator.index_out("xf")).call(mxin)["jac_xf_x0"];
    ca.Function("qeJ", new[] { q0, par }, new[] { qendJ });
    message(qe.call(new[] { A, B })[0].ToString());
  }

  public void test_nl_system() {
    message("Nonlinear ODE sys");
    DM A = ca.vertcat(1, 0.1);
    double p0 = 1.13, y0 = 1, yc0 = 0.1, dy0 = 0.1, te = 0.4;
    SX t = SX.sym("t"), q = SX.sym("y", 2, 1), p = SX.sym("p", 1, 1);
    Function integrator = ca.integrator("integrator", "cvodes",
        new Dictionary<string, SX> { ["x"] = q, ["p"] = p, ["t"] = t, ["ode"] = ca.vertcat(q[1], p[0] + ca.pow(q[1], 2)) }, 0, te,
        opts("reltol", 1e-15, "abstol", 1e-15, "steps_per_checkpoint", 10000, "fsens_err_con", true));
    MX q0 = MX.sym("q0", 2, 1), par = MX.sym("p", 1, 1);
    var mxin = new Dictionary<string, MX> { ["x0"] = q0, ["p"] = par };
    MX qend = integrator.call(mxin)["xf"];
    var qe = ca.Function("qe", new[] { q0, par }, new[] { qend });
    MX qendJ = jacobianOld(integrator, integrator.index_in("x0"), integrator.index_out("xf")).call(mxin)["jac_xf_x0"];
    var qeJ = ca.Function("qeJ", new[] { q0, par }, new[] { qendJ });
    DM qeOut = qe.call(new[] { A, ca.DM(p0) })[0];
    message(qeOut[0].ToString());
    message(qeOut[1].ToString());
    double sp = Math.Sqrt(p0), ta = Math.Tan(Math.Atan(yc0 / sp) + sp * te);
    almost((double)qeOut[0],
        (2 * y0 - Math.Log(yc0 * yc0 / p0 + 1)) / 2 - Math.Log(Math.Cos(Math.Atan(yc0 / sp) + sp * te)), 11, "Nonlin ODE");
    almost((double)qeOut[1], sp * ta, 11, "Nonlin ODE");
    DM qeJOut = qeJ.call(new[] { A, ca.DM(p0) })[0];
    double tb = Math.Tan(sp * te + Math.Atan(dy0 / sp));
    DM Jr = new double[,] { { 1, (sp * tb - dy0) / (dy0 * dy0 + p0) }, { 0, (p0 * tb * tb + p0) / (dy0 * dy0 + p0) } };
    chk(qeJOut, Jr, 10, "jacobian of Nonlin ODE");
    Function Jf = jacobianOld(qe, 0, 0);
    chk(Jf.call(new[] { A, ca.DM(p0) })[0], Jr, 10, "Jacobian of Nonlin ODE");
    Jf = jacobianOld(qe, 0, 0);
    chk(Jf.call(new[] { A, ca.DM(p0) })[0], Jr, 10, "Jacobian of Nonlin ODE");
    double k = te * yc0 * yc0 - yc0 + p0 * te;
    Jr = new double[,] { { (sp * k * ta + yc0 * yc0) / (2 * p0 * yc0 * yc0 + 2 * p0 * p0) },
        { (sp * (k * ta * ta + te * yc0 * yc0 - yc0 + p0 * te) + (yc0 * yc0 + p0) * ta) / (sp * (2 * yc0 * yc0 + 2 * p0)) } };
    Jf = jacobianOld(qe, 1, 0);
    chk(Jf.call(new[] { A, ca.DM(p0) })[0], Jr, 10, "Jacobian of Nonlin ODE");
    Jf = jacobianOld(qe, 1, 0);
    chk(Jf.call(new[] { A, ca.DM(p0) })[0], Jr, 10, "Jacobian of Nonlin ODE");
    qendJ = jacobianOld(integrator, integrator.index_in("p"), integrator.index_out("xf")).call(mxin)["jac_xf_p"];
    qeJ = ca.Function("qeJ", new[] { q0, par }, new[] { qendJ });
    chk(qeJ.call(new[] { A, ca.DM(p0) })[0], Jr, 10, "jacobian of Nonlin ODE");
    checkThreadSafety(qeJ, new List<DM> { A, ca.DM(p0) });
    var qeJf = ca.Function("qeJf", new[] { q0, par }, new[] { ca.vec(qeJ.Call(q0, par)) });
    Function H = jacobianOld(qeJf, 0, 0);
    message(H.call(new[] { A, ca.DM(p0) })[0].ToString());
  }

  public void test_missing_symbols() {
    SX x = SX.sym("x"), z = SX.sym("z"), p = SX.sym("p");
    var dae = new Dictionary<string, SX> { ["x"] = x, ["z"] = z, ["ode"] = z + p, ["alg"] = z * ca.cos(z) - x };
    assertInException("[p] are free", () => ca.integrator("F", "idas", dae));
  }

  public void test_hessian2D() {
    message("hessian");
    int N = 2;
    DM x0_ = ca.vertcat(1, 0.1), A_ = ca.DM(new double[,] { { 3, 1 }, { 0.74, 4 } });
    SX A = SX.sym("A", N, N), x = SX.sym("x", N);
    Function I = ca.integrator("I", "cvodes",
        new Dictionary<string, SX> { ["x"] = x, ["p"] = ca.vec(A), ["ode"] = ca.mtimes(A, x) },
        opts("fsens_err_con", true, "reltol", 1e-12));
    I.call(new Dictionary<string, DM> { ["x0"] = x0_, ["p"] = ca.vec(A_) });
    MX q0 = MX.sym("q0", N), p = MX.sym("p", N * N);
    Dictionary<string, MX> sol = I.call(new Dictionary<string, MX> { ["x0"] = q0, ["p"] = p });
    sol["q0"] = q0;
    sol["p"] = p;
    ca.Function("qe", sol, new[] { "q0", "p" }, ca.integrator_out());
    var JT = ca.Function("JT", new[] { q0, p }, new[] { ca.jacobian(sol["xf"], sol["p"]).T() });
    Function H = jacobianOld(JT, 1, 0);
    _ = H.call(new[] { x0_, ca.vec(A_) })[0];
  }

  public void test_collocationPoints() {
    message("collocation points");
    Check.Throws<Exception>(() => ca.collocation_points(0, "radau"));
    Check.Throws<Exception>(() => ca.collocation_points(10, "radau"));
    Check.Throws<Exception>(() => ca.collocation_points(0, "legendre"));
    Check.Throws<Exception>(() => ca.collocation_points(10, "legendre"));
    Check.Throws<Exception>(() => ca.collocation_points(1, "foo"));
    for (int k = 1; k < 10; k++) {
      List<double> r = ca.collocation_points(k, "radau");
      assertEquals(k + 1, r.Count + 1);
      chk(ca.DM(r[^1]), ca.vertcat(1), 10);
    }
    for (int k = 1; k < 10; k++) assertEquals(k + 1, ca.collocation_points(k, "legendre").Count + 1);
  }

  public void test_thread_safety() {
    MX x = MX.sym("x");
    DM lin = ca.vertcat(linspace(0, 10, 40));
    foreach (Intg I in integrators) {
      foreach (int threads in new[] { 4, 2 }) {
        Function intg = ca.integrator("Integrator", I.name, new Dictionary<string, MX> { ["x"] = x, ["ode"] = -x }, I.options);
        Function intgPar = intg.map(40, "thread", threads);
        DM xf = intgPar.call(new Dictionary<string, DM> { ["x0"] = lin })["xf"];
        chk(ca.norm_inf(xf.T() - lin * Math.Exp(-1)), 0, 5);
      }
    }
  }

  public void test_simplify_zdim() {
    MX x = MX.sym("x");
    Function intg = ca.integrator("intg", "rk", new Dictionary<string, MX> { ["x"] = x, ["ode"] = ca.pow(x, 2) },
        opts("simplify", true));
    assertEquals(0L, intg.nnz_out("zf"));
  }

  public void test_step_options_cvodes() {
    Check.Assume(hasIntegrator("cvodes"), "no integrator cvodes");
    SX x = SX.sym("x");
    Function I = ca.integrator("I", "cvodes", new Dictionary<string, SX> { ["x"] = x, ["ode"] = ca.sin(ca.pow(x * 10, 2)) },
        0.0, 0.5, opts("step0", 1e-4, "min_step_size", 1e-4, "max_step_size", 1.1e-4));
    I.call(new Dictionary<string, DM> { ["x0"] = 0 });
    double n = num(I.stats()["nsteps"]);
    assertTrue((int)(0.5 / 1e-4) >= n && n >= (int)(0.5 / 1.1e-4));
  }

  public void test_step_options_idas() {
    Check.Assume(hasIntegrator("idas"), "no integrator idas");
    SX x = SX.sym("x");
    Function I = ca.integrator("I", "cvodes", new Dictionary<string, SX> { ["x"] = x, ["ode"] = ca.sin(ca.pow(x * 10, 2)) },
        0.0, 0.5, opts("step0", 1e-4, "max_step_size", 1.1e-4));
    I.call(new Dictionary<string, DM> { ["x0"] = 0 });
    double n = num(I.stats()["nsteps"]);
    message("" + n);
    assertTrue(n >= (int)(0.5 / 1.1e-4));
  }

  public void test_adjoint_without_parameters() {
    message("Adjoint sensitivities without parameter or control quadratures");
    MX x = MX.sym("x", 3), x0 = MX.sym("x0", 3);
    DM rates = ca.vertcat(1, 2, 3);
    foreach (string backend in new[] { "cvodes", "idas" }) {
      if (!ca.has_integrator(backend)) continue;
      foreach (double[] grid in new[] { new double[] { 1 }, new double[] { 0.5, 1 } }) {
        Function I = ca.integrator("I", backend, new Dictionary<string, MX> { ["x"] = x, ["ode"] = -rates * x }, 0, grid,
            opts("enable_forward", false, "abstol", 1e-11, "reltol", 1e-11));
        MX xf = I.call(new Dictionary<string, MX> { ["x0"] = x0 })["xf"];
        var J = ca.Function("J", new[] { x0 }, new[] { ca.jacobian(xf[0..2, -1], x0) });
        DM reference = ca.diag(ca.exp(-rates))[0..2, ..];
        foreach (double[] value in new[] { new double[] { 1, 2, 3 }, new double[] { 3, 2, 1 } })
          chk(reference, J.call(new[] { ca.DM(value) })[0], 7);
      }
    }
  }

  public void test_constraints_idas() {
    Check.Assume(hasIntegrator("idas"), "no integrator idas");
    SX x = SX.sym("x"), z = SX.sym("z"), p = SX.sym("p");
    var dae = new Dictionary<string, SX> { ["x"] = x, ["z"] = z, ["p"] = p, ["ode"] = z + p, ["alg"] = z * ca.cos(z) - x };
    Function I = ca.integrator("I", "idas", dae, opts("constraints", new long[] { 1, 1 }));
    I.call(new Dictionary<string, DM> { ["x0"] = 0, ["p"] = 0.15 });
    Function I2 = ca.integrator("I", "idas", dae, opts("constraints", new long[] { -1, -1 }));
    assertInException("IDA_CONSTR_FAIL", () => I2.call(new Dictionary<string, DM> { ["x0"] = 0, ["p"] = 0.15 }));
  }

  public void test_constraints_idas_sensitivities() {
    Check.Assume(hasIntegrator("idas"), "no integrator idas");
    message("IDAS constraints apply to primal states, not sensitivities");
    MX x = MX.sym("x", 2), z = MX.sym("z"), p = MX.sym("p", 2);
    Function I = ca.integrator("I", "idas",
        new Dictionary<string, MX> { ["x"] = x, ["z"] = z, ["p"] = p, ["ode"] = -p * x, ["alg"] = z + ca.sum1(x) }, 0, 1,
        opts("constraints", new long[] { 2, 2, -2 }, "enable_reverse", false, "abstol", 1e-11, "reltol", 1e-11));
    Dictionary<string, MX> result = I.call(new Dictionary<string, MX> { ["x0"] = ca.vertcat(1, 2), ["z0"] = -3, ["p"] = p });
    MX y = ca.vertcat(result["xf"], result["zf"]);
    var J = ca.Function("J", new[] { p }, new[] { ca.jacobian(y, p) });
    double e1 = Math.Exp(-1), e2 = Math.Exp(-2);
    DM reference = new double[,] { { -e1, 0 }, { 0, -2 * e2 }, { e1, 2 * e2 } };
    chk(reference, J.call(new[] { ca.vertcat(1, 2) })[0], 7);
    chk(reference, Function.deserialize(J.serialize()).call(new[] { ca.vertcat(1, 2) })[0], 7);
  }

  static List<SX> permuted(List<SX> s, int[] perm) => perm.Select(e => s[e]).ToList();

  /// <summary>numpy's legacy RandomState (MT19937): np.random.seed(s); np.random.permutation(range(n))</summary>
  sealed class NumpyRandom {
    readonly uint[] mt = new uint[624];
    int mti;

    public NumpyRandom(long seed) {
      mt[0] = (uint)seed;
      for (mti = 1; mti < 624; mti++) mt[mti] = unchecked(1812433253u * (mt[mti - 1] ^ (mt[mti - 1] >> 30)) + (uint)mti);
    }

    uint next() {
      if (mti >= 624) {
        for (int k = 0; k < 624; k++) {
          uint y = (mt[k] & 0x80000000u) | (mt[(k + 1) % 624] & 0x7fffffffu);
          mt[k] = mt[(k + 397) % 624] ^ (y >> 1) ^ ((y & 1) != 0 ? 0x9908b0dfu : 0);
        }
        mti = 0;
      }
      uint z = mt[mti++];
      z ^= z >> 11;
      z ^= (z << 7) & 0x9d2c5680u;
      z ^= (z << 15) & 0xefc60000u;
      z ^= z >> 18;
      return z;
    }

    long interval(long max) {
      if (max == 0) return 0;
      long mask = max, v;
      for (int s = 1; s <= 32; s <<= 1) mask |= mask >> s;
      while ((v = next() & mask) > max) {}
      return v;
    }

    public int[] permutation(int n) {
      int[] p = Enumerable.Range(0, n).ToArray();
      for (int i = n - 1; i > 0; i--) {
        int j = (int)interval(i);
        (p[i], p[j]) = (p[j], p[i]);
      }
      return p;
    }
  }

  static long[] index(int[] perm) => perm.Select(e => (long)e).ToArray();

  sealed class IndexProblem {
    public Dictionary<string, SX> dae;
    public object index;  // expected index (long) or expected error message (string)
    public Dictionary<string, DM> strength, init;
    public double tol;
  }

  static List<IndexProblem> reduceIndexProblems() {
    double free = 0, force = -1, approx = 1, blank = 0.0;
    var r = new List<IndexProblem>();
    double L = 1, g = 9.81;
    {
      // 2D Pendulum
      SX x = SX.sym("x"), y = SX.sym("y"), dx = SX.sym("dx"), dy = SX.sym("dy"), u = SX.sym("u"), v = SX.sym("v");
      SX du = SX.sym("du"), dv = SX.sym("dv"), T = SX.sym("T");
      SX alg = ca.vertcat(dx - u, dy - v, du - T * x, dv - T * y + g, ca.pow(x, 2) + ca.pow(y, 2) - L * L);
      r.Add(new IndexProblem {
        dae = new Dictionary<string, SX> { ["x_impl"] = ca.vertcat(x, y, u, v), ["dx_impl"] = ca.vertcat(dx, dy, du, dv),
                                           ["z"] = T, ["alg"] = alg },
        index = 3L,
        strength = new Dictionary<string, DM> { ["x_impl"] = ca.vertcat(free, force, free, free),
                                                ["dx_impl"] = ca.vertcat(force, free, free, free), ["z"] = ca.DM(free) },
        init = new Dictionary<string, DM> { ["x_impl"] = ca.vertcat(-1.0, -0.5, blank, blank),
                                            ["dx_impl"] = ca.vertcat(-0.1, blank, blank, blank), ["z"] = ca.DM(blank) },
        tol = 5e-6 });
    }
    {
      // 3D Pendulum
      SX x = SX.sym("x"), y = SX.sym("y"), z = SX.sym("z"), vx = SX.sym("vx"), vy = SX.sym("vy"), vz = SX.sym("vz");
      SX dx = SX.sym("dx"), dy = SX.sym("dy"), dz = SX.sym("dz"), dvx = SX.sym("dvx"), dvy = SX.sym("dvy"), dvz = SX.sym("dvz");
      SX lambd = SX.sym("lambd");
      SX alg = ca.vertcat(dx - vx, dy - vy, dz - vz, dvx - lambd * x, dvy - lambd * y + g, dvz - lambd * z,
          ca.pow(x, 2) + ca.pow(y, 2) + ca.pow(z, 2) - L * L);
      r.Add(new IndexProblem {
        dae = new Dictionary<string, SX> { ["x_impl"] = ca.vertcat(x, y, z, vx, vy, vz),
                                           ["dx_impl"] = ca.vertcat(dx, dy, dz, dvx, dvy, dvz), ["z"] = lambd, ["alg"] = alg },
        index = 3L,
        strength = new Dictionary<string, DM> { ["x_impl"] = ca.vertcat(free, force, force, free, free, free),
                                                ["dx_impl"] = ca.vertcat(force, force, free, free, free, free), ["z"] = ca.DM(free) },
        init = new Dictionary<string, DM> { ["x_impl"] = ca.vertcat(-1.0, -0.5, -0.5, blank, blank, blank),
                                            ["dx_impl"] = ca.vertcat(-0.1, 0.2, blank, blank, blank, blank), ["z"] = ca.DM(blank) },
        tol = 1e-5 });
    }
    {
      // Flash system (CppADCodeGen dae_index_reduction/model/flash.hpp)
      SX nEthanol = SX.sym("nEthanol"), nWater = SX.sym("nWater"), T = SX.sym("T"), yWater = SX.sym("yWater");
      SX yEthanol = SX.sym("yEthanol"), FV = SX.sym("FV");
      double Q = 500, F_feed = 10, p = 1, xFEthanol = 0.5, T_feed = 50;
      SX dnEthanol = SX.sym("dnEthanol"), dnWater = SX.sym("dnWater"), dT = SX.sym("dT");
      var alg = new List<SX>();
      double FNEthanol = xFEthanol * F_feed;
      SX n_lWater = nWater * 1000.0, n_lEthanol = nEthanol * 1000.0;
      SX m_lEthanol = n_lEthanol * 0.04606844, m_lWater = n_lWater * 0.0180152833;
      SX m = m_lEthanol + m_lWater;
      SX wEthanol = m_lEthanol / m;
      SX w_Water = 1 - wEthanol;
      SX rho = 1 / (w_Water / 983.159471259596 + wEthanol / 743.278841365274);
      SX V = m / rho;
      SX cWater = n_lWater / V;
      SX dh = V / 0.502654824574367;
      double p_aux = p * 101325.0, dp = 101325.0 - p_aux;
      SX v = ca.sqrt((9.80665 * dh + dp / rho) * 2.0);
      SX F_VL = v * 0.0005;
      SX FNLWater = cWater * F_VL;
      SX n = n_lEthanol + n_lWater;
      SX xWater = n_lWater / n;
      SX F_NL = FNLWater / xWater;
      SX FNLEthanol = F_NL - FNLWater;
      SX FNVWater = yWater * FV;
      SX FNVEthanol = FV - FNVWater;
      SX dnE = FNEthanol - FNLEthanol - FNVEthanol;
      alg.Add(dnEthanol - dnE * 0.001);
      double FNFWater = F_feed - FNEthanol;
      SX dnW = FNFWater - FNLWater - FNVWater;
      alg.Add(dnWater - dnW * 0.001);
      double FMFEthanol = FNEthanol * 0.04606844, FMFWater = FNFWater * 0.0180152833, Tfeed = T_feed + 273.15;
      SX T_aux = T + 273.15;
      SX dQ_F = (FMFEthanol * 2898.42878374706 + FMFWater * 4186.92536027523) * (Tfeed - T_aux);
      SX dH_mVapEthanol = 50430.0 * ca.pow(1 - T_aux / 514.0, 0.4989) * ca.exp(0.4475 * T_aux / 514.0);
      SX T_r = T_aux / 647.0;
      SX dH_mVapWater = 52053.0 * ca.pow(1 - T_r, 0.3199 + -0.212 * T_r + 0.25795 * T_r * T_r);
      SX dH_mVap = yEthanol * dH_mVapEthanol + yWater * dH_mVapWater;
      SX dQ_vap = FV * dH_mVap;
      double Q_1 = Q * 1000.0;
      SX cp = wEthanol * 2898.42878374706 + w_Water * 4186.92536027523;
      alg.Add(dT - (dQ_F - dQ_vap + Q_1) / (m * cp));
      SX pVapWater = 100000.0 * ca.pow(10.0, 4.6543 - 1435.264 / (T_aux - 64.848));
      SX KWater = pVapWater / p_aux;
      alg.Add(yWater - xWater * KWater);
      SX xEthanol = 1 - xWater;
      SX pVapEthanol = ca.exp(74.475 + -7164.3 / T_aux + -7.327 * ca.log(T_aux) + 3.134e-06 * ca.pow(T_aux, 2.0));
      SX KEthanol = pVapEthanol / p_aux;
      alg.Add(yEthanol - xEthanol * KEthanol);
      alg.Add(yWater + yEthanol - 1);
      r.Add(new IndexProblem {
        dae = new Dictionary<string, SX> { ["x_impl"] = ca.veccat(nEthanol, nWater, T), ["dx_impl"] = ca.veccat(dnEthanol, dnWater, dT),
                                           ["z"] = ca.vertcat(yWater, yEthanol, FV), ["alg"] = ca.veccat(alg) },
        index = 2L,
        strength = new Dictionary<string, DM> { ["x_impl"] = ca.vertcat(force, force, approx), ["dx_impl"] = ca.vertcat(free, free, free),
                                                ["z"] = ca.vertcat(free, free, free) },
        init = new Dictionary<string, DM> { ["x_impl"] = ca.vertcat(2.5, 6.4, 91), ["dx_impl"] = ca.vertcat(0.1, blank, blank),
                                            ["z"] = ca.vertcat(0.5, 0.5, 0.5) },
        tol = 1e-8 });
    }
    {
      // Exothermic Reactor from pantelides paper
      SX C = SX.sym("C"), dC = SX.sym("dC"), T = SX.sym("T"), dT = SX.sym("dT"), R = SX.sym("R");
      double C0 = 1, T0 = 1, K1 = 1, K2 = 1, K3 = 1, K4 = 1, Tc = 1;
      SX alg = ca.vertcat(K1 * (C0 - C) - R - dC, K1 * (T0 - T) + R * K2 - (T - Tc) * K3 - dT, R - K3 * ca.exp(-K4 / T) * C);
      r.Add(new IndexProblem {
        dae = new Dictionary<string, SX> { ["x_impl"] = ca.veccat(C, T), ["dx_impl"] = ca.veccat(dC, dT), ["z"] = ca.vertcat(R),
                                           ["alg"] = alg },
        index = 1L,
        strength = new Dictionary<string, DM> { ["x_impl"] = ca.vertcat(free, force), ["dx_impl"] = ca.vertcat(free, free), ["z"] = ca.DM(free) },
        init = new Dictionary<string, DM> { ["x_impl"] = ca.vertcat(blank, 2), ["dx_impl"] = ca.vertcat(blank, blank), ["z"] = ca.DM(blank) },
        tol = 1e-8 });
    }
    {
      // S. Mattsson and G. Soederlind, Index reduction in differential-algebraic equations using dummy derivatives
      SX[] xs = new SX[7], dxs = new SX[7];
      for (int i = 0; i < 7; i++) {
        xs[i] = SX.sym("x" + (i + 1));
        dxs[i] = SX.sym("dx" + (i + 1));
      }
      double u1 = 1, u2 = 1, u3 = 1, u4 = 1;
      SX eqs = ca.vertcat(xs[4] - dxs[0], xs[5] - dxs[1], xs[6] - dxs[2], xs[0] + xs[1] + u1, xs[0] + xs[1] + xs[2] + u2,
          dxs[2] + xs[0] + xs[3] + u3, dxs[4] * 2 + dxs[5] + dxs[6] + dxs[3] + u4);
      r.Add(new IndexProblem {
        dae = new Dictionary<string, SX> { ["x_impl"] = ca.veccat(xs), ["dx_impl"] = ca.veccat(dxs), ["z"] = ca.SX(0, 1),
                                           ["alg"] = eqs },
        index = 3L });
    }
    {
      // Uncontrollable DAE Pantelides
      SX x = SX.sym("x"), dx = SX.sym("dx"), u1 = SX.sym("u1"), u2 = SX.sym("u2"), t = SX.sym("t");
      SX y1 = t, y2 = ca.sin(t);
      SX alg = ca.vertcat(x + u1 + u2, x + dx + y1, x + y2);
      r.Add(new IndexProblem {
        dae = new Dictionary<string, SX> { ["x_impl"] = ca.vertcat(x), ["dx_impl"] = ca.vertcat(dx), ["z"] = ca.vertcat(u1, u2),
                                           ["alg"] = alg, ["t"] = t },
        index = "Structural detection failed" });
    }
    return r;
  }

  public void test_reduce_index() {
    Check.Assume(hasIntegrator("idas"), "no integrator idas");
    Check.Assume(hasNlpsol("ipopt"), "no nlpsol ipopt");
    message("loop");
    foreach (IndexProblem pr in reduceIndexProblems()) {
      message("dae " + str(pr.dae));
      for (int seed = 0; seed < 5; seed++) {
        var dae = new Dictionary<string, SX>(pr.dae);
        var strength = pr.strength == null ? null : new Dictionary<string, DM>(pr.strength);
        var init = pr.init == null ? null : new Dictionary<string, DM>(pr.init);
        if (seed > 0) {
          var rng = new NumpyRandom(seed);
          foreach (string k in new[] { "x_impl", "z", "alg" }) {
            List<SX> s = ca.vertsplit(dae[k]);
            if (s.Count == 0) continue;
            int[] perm = rng.permutation(s.Count);
            dae[k] = ca.vertcat(permuted(s, perm));
            if (init != null && init.ContainsKey(k)) init[k] = init[k][index(perm)];
            if (strength != null && strength.ContainsKey(k)) strength[k] = strength[k][index(perm)];
            if (k == "x_impl") {
              s = ca.vertsplit(dae["dx_impl"]);
              dae["dx_impl"] = ca.vertcat(permuted(s, perm));
              if (init != null) init["dx_impl"] = init["dx_impl"][index(perm)];
              if (strength != null) strength["dx_impl"] = strength["dx_impl"][index(perm)];
            }
          }
        }
        if (pr.index is string err) {
          assertInException(err, () => ca.dae_reduce_index(dae,
              opts("baumgarte_pole", -1, "algorithm", "pantelides", "max_iter", 20)));
          continue;
        }
        var (daeReduced, meta) = ca.dae_reduce_index(dae, opts("baumgarte_pole", -1, "algorithm", "pantelides"));
        assertEquals((long)pr.index, Convert.ToInt64(meta["index"]));
        if (strength != null) {
          var (daeSe, _, phi) = ca.dae_map_semi_expl(dae, daeReduced);
          Function intg = ca.integrator("intg", "idas", daeSe, 0, linspace(0, 5, 1000));
          Function initGen = ca.dae_init_gen(dae, daeReduced, "ipopt", strength, opts("error_on_fail", true));
          Dictionary<string, DM> xz0 = initGen.call(init);
          Dictionary<string, DM> sol = intg.call(xz0);
          DM error = phi.call(new Dictionary<string, DM> { ["x"] = sol["xf"], ["z"] = sol["zf"] })["I"];
          double[,] e = error.ToArray();
          for (int i = 0; i < e.GetLength(0); i++) {
            double s = 0;
            for (int j = 0; j < e.GetLength(1); j++) s += e[i, j] * e[i, j];
            double rms = Math.Sqrt(s / e.GetLength(1));
            message("rms error " + rms);
            assertTrue(rms <= pr.tol);
          }
        }
      }
    }
  }

  public void test_IDACalcIC_u() {
    MX x = MX.sym("x"), u = MX.sym("u"), z = MX.sym("z");
    var ivp = new Dictionary<string, MX> { ["x"] = x, ["z"] = z, ["u"] = u, ["ode"] = z, ["alg"] = ca.sqrt(u) - z };
    Function intg = ca.integrator("intg", "idas", ivp, 0, 1);
    Dictionary<string, DM> res = intg.call(new Dictionary<string, DM> { ["x0"] = 0, ["u"] = 1 });
    almost((double)res["xf"], 1, 5, "test_IDACalcIC_u");
  }

  public void test_expand() {
    Check.Assume(hasIntegrator("cvodes"), "no integrator cvodes");
    MX t = MX.sym("t"), q = MX.sym("q"), p = MX.sym("p");
    var dae = new Dictionary<string, MX> { ["t"] = t, ["x"] = q, ["p"] = p, ["ode"] = q / p * ca.pow(t, 2) };
    Function integrator = ca.integrator("integrator", "cvodes", dae, 0, 2.3, opts("expand", true));
    assertTrue(integrator.get_function("jacF").is_a("SXFunction"));
    t = MX.sym("t");
    q = MX.sym("q");
    p = MX.sym("p");
    dae = new Dictionary<string, MX> { ["t"] = t, ["x"] = q, ["p"] = p, ["ode"] = q / p * ca.pow(t, 2) };
    integrator = ca.integrator("integrator", "cvodes", dae, 0, 2.3, opts("expand", false));
    assertTrue(integrator.get_function("jacF").is_a("MXFunction"));
  }

  public void test_postpone_expand() {
    Check.Assume(hasIntegrator("cvodes"), "no integrator cvodes");
    MX x = MX.sym("x");
    var J = ca.Function("jac_f", new[] { x, ca.MX(1, 1) }, new[] { -x }, new[] { "x", "out_r" },
        new[] { "jac_r_x" }, opts("always_inline", true));
    var f = ca.Function("f", new[] { x }, new[] { ca.pow(x, 2) }, new[] { "x" }, new[] { "r" },
        opts("custom_jacobian", J, "jac_penalty", 0));
    var dae = new Dictionary<string, MX> { ["x"] = x, ["ode"] = f.Call(x - 2) };
    var cases = new (Dictionary<string, object>, string, double)[] {
      (opts("expand", false), "MXFunction", 1.0), (opts("expand", true), "SXFunction", -2.0),
      (opts("expand", false, "postpone_expand", true), "MXFunction", 1.0),
      (opts("expand", true, "postpone_expand", true), "SXFunction", 1.0) };
    foreach (var (o, cls, val) in cases) {
      Function integrator = ca.integrator("integrator", "cvodes", dae, 0, 0.1, o);
      integrator.call(new Dictionary<string, DM> { ["x0"] = 1 });
      assertTrue(integrator.get_function("jacF").is_a(cls));
      chk(integrator.get_function("jacF").call(new Dictionary<string, DM> { ["x"] = 1 })["jac_ode_x"], val, 10);
    }
  }

  public void test_idas_adjoint_zf_p() {
    Check.Assume(hasIntegrator("idas"), "no integrator idas");
    MX x = MX.sym("x"), z = MX.sym("z"), p = MX.sym("p"), u = MX.sym("u");
    foreach (MX alg in new[] { z - p * x - ca.pow(u, 2), z - p * x }) {
      foreach (MX ode in new[] { -x, -x + u }) {
        double[] grid = { 1 };
        Function F = ca.integrator("F", "idas",
            new Dictionary<string, MX> { ["x"] = x, ["z"] = z, ["p"] = p, ["u"] = u, ["ode"] = ode, ["alg"] = alg }, 0, grid,
            opts("abstol", 1e-12, "reltol", 1e-12));
        MX P = MX.sym("P"), U = MX.sym("U", 1, grid.Length);
        Dictionary<string, MX> o = F.call(new Dictionary<string, MX> { ["x0"] = 2, ["z0"] = 0, ["p"] = P, ["u"] = U });
        MX e = ca.vertcat(ca.vec(o["zf"]), ca.vec(o["xf"]));
        MX w = ca.vertcat(P, ca.vec(U));
        var Jfwd = ca.Function("J_fwd", new[] { w }, new[] { ca.jacobian(e, w, opts("helper_options", opts("ad_weight", 0))) });
        var Jrev = ca.Function("J_rev", new[] { w }, new[] { ca.jacobian(e, w, opts("helper_options", opts("ad_weight", 1))) });
        DM v = ca.vertcat(new double[] { 0.5, 0.3, 0.7 }.Take(1 + grid.Length).ToArray());
        chk(Jfwd.call(new[] { v })[0], Jrev.call(new[] { v })[0], 6);
      }
    }
  }

  public void test_jit_fwd_factory() {
    Check.Assume(hasIntegrator("cvodes"), "no integrator cvodes");
    Check.Assume(RunSlow, "slow");
    MX x = MX.sym("x");
    Function F = ca.integrator("F", "cvodes", new Dictionary<string, MX> { ["x"] = x, ["ode"] = -x }, 0,
        new double[] { 0.5, 1 }, opts("abstol", 1e-10, "reltol", 1e-10, "jit", true, "compiler", "shell"));
    assertTrue(F.get_function("jacF").is_a("External"));
    Function Ifwd = F.factory("I_fwd", new[] { "x0", "fwd:x0" }, new[] { "xf", "fwd:xf" });
    List<DM> r = Ifwd.call(new[] { ca.vertcat(2.0), ca.vertcat(1.0) });
    DM refv = ca.exp(-ca.DM(new double[,] { { 0.5, 1 } }));
    chk(r[0], refv * 2, 7);
    chk(r[1], refv, 7);
  }

  void crossCheck(Function I, Dictionary<string, DM> b, List<Dictionary<string, DM>> seeds, string msg) {
    int N = seeds.Count;
    Function I1 = I.forward(1), IN = I.forward(N);
    var fwdNames = IN.name_in().Where(n => n.StartsWith("fwd_")).ToList();
    var batchedCall = new Dictionary<string, DM>(b);
    foreach (string n in fwdNames)
      batchedCall[n] = ca.horzcat(seeds.Select(s => s.TryGetValue(n, out DM v) ? v : DM.zeros(I1.sparsity_in(n))).ToList());
    Dictionary<string, DM> batched = IN.call(batchedCall);
    var perSeed = new Dictionary<string, List<DM>>();
    foreach (var s in seeds) {
      var single = new Dictionary<string, DM>(b);
      foreach (string n in fwdNames) single[n] = s.TryGetValue(n, out DM v) ? v : DM.zeros(I1.sparsity_in(n));
      foreach (var e in I1.call(single)) {
        if (!perSeed.ContainsKey(e.Key)) perSeed[e.Key] = new List<DM>();
        perSeed[e.Key].Add(e.Value);
      }
    }
    foreach (var e in batched) chk(e.Value, ca.horzcat(perSeed[e.Key]), 10, msg + " :: " + e.Key);
  }

  public void test_issue3371() {
    MX X = MX.sym("x", 2), U = MX.sym("u"), Z = MX.sym("z");
    var dae = new Dictionary<string, MX> { ["x"] = X, ["p"] = U, ["z"] = Z, ["ode"] = ca.vertcat(Z * X[0] - X[1] + U, X[0]),
                                           ["alg"] = Z - (1 - ca.pow(X[1], 2)) };
    Function I = ca.integrator("I", "collocation", dae);
    var b = new Dictionary<string, DM> { ["x0"] = ca.vertcat(0.5, 0.0), ["z0"] = ca.vertcat(1.0), ["p"] = ca.vertcat(0.1), ["u"] = ca.DM(0, 1) };
    Dictionary<string, DM> r = I.call(new Dictionary<string, DM> { ["x0"] = b["x0"], ["z0"] = b["z0"], ["p"] = b["p"] });
    b["out_xf"] = r["xf"];
    b["out_zf"] = r["zf"];
    b["out_qf"] = r["qf"];
    crossCheck(I, b, new List<Dictionary<string, DM>> {
      new Dictionary<string, DM> { ["fwd_x0"] = ca.vertcat(1.0, 0.0) }, new Dictionary<string, DM> { ["fwd_x0"] = ca.vertcat(0.0, 1.0) },
      new Dictionary<string, DM> { ["fwd_p"] = ca.vertcat(1.0) } }, "fwd-of-collocation");

    Function I2 = I.forward(2);
    var base0 = new Dictionary<string, DM> {
      ["x0"] = DM.zeros(2), ["z0"] = ca.vertcat(1.0), ["p"] = ca.vertcat(0.0), ["u"] = ca.DM(0, 1),
      ["out_xf"] = DM.zeros(2, 1), ["out_zf"] = ca.vertcat(1.0), ["out_qf"] = ca.DM(0, 1),
      ["fwd_x0"] = ca.DM(new double[,] { { 0, 1 }, { 0, 0 } }), ["fwd_z0"] = DM.zeros(I2.sparsity_in("fwd_z0")),
      ["fwd_p"] = ca.vertcat(1, 0).T(), ["fwd_u"] = DM.zeros(I2.sparsity_in("fwd_u")) };
    chk(I2.call(base0)["fwd_zf"], DM.zeros(1, 2), 10, "issue3371 fwd_zf");

    Function Iadj = ca.integrator("I", "collocation", dae, opts("nadj", 1));
    Dictionary<string, DM> refr = Iadj.call(new Dictionary<string, DM> {
      ["x0"] = ca.vertcat(0.5, 0.0), ["z0"] = ca.vertcat(1.0), ["p"] = ca.vertcat(0.1), ["adj_xf"] = ca.vertcat(0.0, 0.0), ["adj_zf"] = ca.vertcat(1.0),
      ["adj_qf"] = ca.DM(0, 1) });
    b = new Dictionary<string, DM> {
      ["x0"] = ca.vertcat(0.5, 0.0), ["z0"] = ca.vertcat(1.0), ["p"] = ca.vertcat(0.1), ["u"] = ca.DM(0, 1), ["adj_xf"] = ca.vertcat(0.0, 0.0),
      ["adj_zf"] = ca.vertcat(1.0), ["adj_qf"] = ca.DM(0, 1), ["out_xf"] = refr["xf"], ["out_zf"] = refr["zf"],
      ["out_qf"] = refr["qf"], ["out_adj_x0"] = refr["adj_x0"], ["out_adj_z0"] = refr["adj_z0"],
      ["out_adj_p"] = refr["adj_p"], ["out_adj_u"] = refr["adj_u"] };
    crossCheck(Iadj, b, new List<Dictionary<string, DM>> {
      new Dictionary<string, DM> { ["fwd_adj_zf"] = ca.vertcat(1.0) }, new Dictionary<string, DM> { ["fwd_x0"] = ca.vertcat(1.0, 0.0) },
      new Dictionary<string, DM> { ["fwd_p"] = ca.vertcat(1.0) } }, "fwd-of-reverse-of-collocation");
  }

  /// <summary>Forward and adjoint Jacobians of out(v) at v0 agree</summary>
  void checkJacobian(MX o, MX v, DM v0, int digits) {
    DM Jfwd = ca.Function("J", new[] { v }, new[] { ca.jtimes(o, v, MX.eye(v.numel())) }).call(new[] { v0 })[0];
    DM Jadj = ca.Function("J", new[] { v }, new[] { ca.jtimes(o, v, MX.eye(o.numel()), true).T() }).call(new[] { v0 })[0];
    chk(Jfwd, Jadj, digits, "fwd vs adj");
  }

  /// <summary>Hessian of obj(v) at v0: forward and adjoint mode over the forward and adjoint gradients</summary>
  DM checkHessian(MX obj, MX v, DM v0, int digits) {
    long n = v.numel();
    MX gFwd = ca.jtimes(obj, v, MX.eye(n)).T();
    MX gAdj = ca.jtimes(obj, v, MX.ones(1, 1), true);
    var H = new[] { ca.jtimes(gFwd, v, MX.eye(n)), ca.jtimes(gAdj, v, MX.eye(n)),
                    ca.jtimes(gFwd, v, MX.eye(n), true).T(), ca.jtimes(gAdj, v, MX.eye(n), true).T() }
      .Select(h => ca.Function("H", new[] { v }, new[] { h }).call(new[] { v0 })[0]).ToList();
    foreach (DM h in H.Skip(1)) chk(H[0], h, digits, "Hessian");
    return H[0];
  }

  static DM randUniform(Random rng, Sparsity sp) {
    double[] v = new double[sp.nnz()];
    for (int k = 0; k < v.Length; k++) v[k] = rng.NextDouble();
    return ca.DM(sp, v.Length == 0 ? ca.DM(0, 1) : (DM)v);
  }

  static List<DM> join(params IEnumerable<DM>[] a) => a.SelectMany(e => e).ToList();

  public void test_reverse_grid_directions() {
    message("Several adjoint directions with several output times");
    SX x = SX.sym("x", 2), p = SX.sym("p");
    var dae = new Dictionary<string, SX> { ["x"] = x, ["p"] = p, ["ode"] = ca.vertcat(x[1], -p * ca.pow(x[0], 3)),
                                           ["quad"] = ca.pow(x[0], 2) };
    var cases = new (string, Dictionary<string, object>)[] {
      ("cvodes", opts("reltol", 1e-10, "abstol", 1e-10)), ("idas", opts("reltol", 1e-10, "abstol", 1e-10)),
      ("rosenbrock", opts("step_size", 0.125)), ("erk", opts("method", "rk45", "reltol", 1e-8, "abstol", 1e-8)) };
    foreach (var (plugin, o) in cases) {
      if (!ca.has_integrator(plugin)) continue;
      Function F = ca.integrator("F", plugin, dae, 0, new double[] { 0.5, 1.0 }, o);
      var inputs = new List<DM> { ca.vertcat(0.3, -0.2), ca.DM(), ca.DM(0.7), ca.DM(1, 0) };
      for (long i = 4; i < F.n_in(); i++) inputs.Add(ca.DM(F.sparsity_in(i)));
      List<DM> outputs = F.call(inputs);
      var seeds = new List<List<DM>>();
      for (int d = 0; d < 2; d++)
        seeds.Add(Enumerable.Range(0, (int)F.n_out()).Select(i => randUniform(new Random(10 * d + i), F.sparsity_out(i))).ToList());
      var both = Enumerable.Range(0, (int)F.n_out()).Select(i => ca.horzcat(seeds[0][i], seeds[1][i])).ToList();
      List<DM> adj2 = F.reverse(2).call(join(inputs, outputs, both));
      for (int d = 0; d < 2; d++) {
        List<DM> adj1 = F.reverse(1).call(join(inputs, outputs, seeds[d]));
        for (int i = 0; i < F.n_in(); i++) {
          if (F.nnz_in(i) == 0) continue;
          chk(adj1[i], ca.horzsplit_n(adj2[i], 2)[d], 8, plugin);
        }
      }
    }
  }

  public void test_events_sparsity() {
    message("zero-crossing events: dependencies through the transition");
    double g = 9.81;
    SX x = SX.sym("x", 2), p = SX.sym("p");
    var dae = new Dictionary<string, SX> { ["x"] = x, ["p"] = p, ["ode"] = ca.vertcat(x[1], -g), ["zero"] = x[0] };
    double[] grid = linspace(0.1, 4, 40);
    var cases = new (string, Dictionary<string, object>)[] {
      ("rk", opts("number_of_finite_elements", 100)), ("cvodes", opts("abstol", 1e-12, "reltol", 1e-12)) };
    foreach (var (plugin, o) in cases) {
      if (!ca.has_integrator(plugin)) continue;
      Function F = ca.integrator("F", plugin, dae, 0, grid, with(o, "transition", bounceTransition(), "event_tol", 1e-12,
          "max_event_iter", 20));
      MX v = MX.sym("v", 3);
      MX xf = F.call(new Dictionary<string, MX> { ["x0"] = v[0..2], ["p"] = v[2] })["xf"][.., -1];
      var G = ca.Function("G", new[] { v }, new[] { xf, ca.jacobian(xf, v) });
      assertTrue(G.sparsity_out(1).is_dense());
      DM v0 = ca.vertcat(2, 0.3, 0.8);
      chk(G.call(new[] { v0 })[1], fdJacobian(G, v0), 5, plugin);
    }
  }

  public void test_events_idas() {
    Check.Assume(hasIntegrator("idas"), "no integrator idas");
    message("idas: zero-crossing events, also of a DAE");
    double g = 9.81;
    SX x = SX.sym("x", 2), z = SX.sym("z"), p = SX.sym("p");
    double[] grid = linspace(0.1, 1.5, 15);
    double tb = Math.Sqrt(4 / g), ub = 0.8 * g * tb;
    double[,] hv = new double[2, grid.Length];
    for (int i = 0; i < grid.Length; i++) {
      double tt = grid[i];
      hv[0, i] = tt < tb ? 2 - g * tt * tt / 2 : ub * (tt - tb) - g * Math.Pow(tt - tb, 2) / 2;
      hv[1, i] = tt < tb ? -g * tt : ub - g * (tt - tb);
    }
    var ode = new Dictionary<string, SX> { ["x"] = x, ["p"] = p, ["ode"] = ca.vertcat(x[1], -g), ["zero"] = x[0] };
    var dae = new Dictionary<string, SX> { ["x"] = x, ["z"] = z, ["p"] = p, ["ode"] = ca.vertcat(z, -g), ["alg"] = z - x[1],
                                           ["zero"] = x[0] };
    foreach (var (sys, transition) in new[] { (ode, bounceTransition()), (dae, bounceTransition(1)) }) {
      Function F = ca.integrator("F", "idas", sys, 0, grid, opts("transition", transition,
          "abstol", 1e-12, "reltol", 1e-12, "event_tol", 1e-12, "max_event_iter", 20));
      Dictionary<string, DM> r = F.call(new Dictionary<string, DM> { ["x0"] = ca.vertcat(2, 0), ["z0"] = 0, ["p"] = 0.8 });
      chk(r["xf"], hv, 7, "xf");
    }
  }

  public void test_reverse_over_adjoint_directions() {
    message("Derivatives of adjoint sensitivities with several directions (#3293)");
    SX x = SX.sym("x", 3), p = SX.sym("p", 2);
    var dae = new Dictionary<string, SX> { ["x"] = x, ["p"] = p,
      ["ode"] = ca.vertcat(-p[0] * x[0], -p[1] * x[1] * x[0], -p[1] * ca.pow(x[2], 2)), ["quad"] = x[0] * x[2] };
    var rng = new Random(1);
    var cases = new (string, Dictionary<string, object>, int)[] {
      ("cvodes", opts("reltol", 1e-12, "abstol", 1e-12), 7), ("idas", opts("reltol", 1e-12, "abstol", 1e-12), 7),
      ("rosenbrock", opts("step_size", 1.0 / 6), 12), ("rosenbrock", opts("reltol", 1e-6, "abstol", 1e-6), 11),
      ("erk", opts("step_size", 0.2), 12), ("erk", opts("method", "rk45", "reltol", 1e-6, "abstol", 1e-6), 11) };
    foreach (var (plugin, o, digits) in cases) {
      if (!ca.has_integrator(plugin)) continue;
      foreach (double[] grid in new[] { new[] { 1.0 }, new[] { 0.5, 1.0 } }) {
        Function F = ca.integrator("F", plugin, dae, 0, grid, o);
        Function B = F.reverse(2);
        var inputs = new List<DM>();
        var v = new List<DM>();
        for (int i = 0; i < B.n_in(); i++) inputs.Add(randUniform(rng, B.sparsity_in(i)));
        List<DM> outputs = B.call(inputs);
        for (int i = 0; i < B.n_in(); i++) v.Add(randUniform(rng, B.sparsity_in(i)));
        var w = new List<List<DM>>();
        for (int d = 0; d < 2; d++) {
          var s = new List<DM>();
          for (int i = 0; i < B.n_out(); i++) s.Add(randUniform(rng, B.sparsity_out(i)));
          w.Add(s);
        }
        List<DM> fwd = B.forward(1).call(join(inputs, outputs, v));
        var both = Enumerable.Range(0, (int)B.n_out()).Select(i => ca.horzcat(w[0][i], w[1][i])).ToList();
        List<DM> adj2 = B.reverse(2).call(join(inputs, outputs, both));
        for (int d = 0; d < 2; d++) {
          List<DM> adj = B.reverse(1).call(join(inputs, outputs, w[d]));
          double lhs = 0, rhs = 0;
          for (int i = 0; i < B.n_out(); i++) lhs += (double)ca.dot(w[d][i], fwd[i]);
          for (int i = 0; i < B.n_in(); i++) rhs += (double)ca.dot(adj[i], v[i]);
          chk(ca.DM(lhs), ca.DM(rhs), digits, plugin + " duality");
          for (int i = 0; i < B.n_in(); i++)
            if (B.nnz_in(i) > 0) chk(adj[i], ca.horzsplit_n(adj2[i], 2)[d], digits, plugin);
        }
        long[] offsets = new long[(int)B.n_in() + 1];
        for (int i = 0; i < B.n_in(); i++) offsets[i + 1] = offsets[i] + B.nnz_in(i);
        long total = offsets[^1];
        MX vin = MX.sym("v", total);
        List<MX> args = ca.vertsplit(vin, offsets);
        var resIn = Enumerable.Range(0, (int)B.n_in()).Select(i => ca.MX(B.sparsity_in(i), args[i])).ToList();
        MX vout = ca.vertcat(B.call(resIn).Select(e => ca.vec(e.get_nz(false, ca.Slice()))).ToList());
        var J = ca.Function("J", new[] { vin }, new[] { ca.jacobian(vout, vin), ca.jtimes(vout, vin, MX.eye(vin.numel())) });
        List<DM> Jo = J.call(new[] { ca.vertcat(inputs.SelectMany(e => e.nonzeros()).ToArray()) });
        chk(Jo[1], Jo[0], digits, plugin + " sparsity");
      }
    }
  }

  public void test_runtime_integrator_quadrature_error_control() {
    SX t = SX.sym("t"), x = SX.sym("x");
    var dae = new Dictionary<string, SX> { ["t"] = t, ["x"] = x, ["ode"] = 0 * x, ["quad"] = ca.exp(t * 10) };
    double exact = (Math.Exp(10) - 1) / 10;
    var o = opts("abstol", 1e-9, "reltol", 1e-9, "step0", 1.0);
    foreach (string plugin in new[] { "erk", "rosenbrock" }) {
      if (!ca.has_integrator(plugin)) continue;
      Function def = ca.integrator("default", plugin, dae, 0, 1, o);
      Function enabled = ca.integrator("enabled", plugin, dae, 0, 1, with(o, "quad_err_con", true));
      Function disabled = ca.integrator("disabled", plugin, dae, 0, 1, with(o, "quad_err_con", false));
      var x00 = new Dictionary<string, DM> { ["x0"] = 0 };
      double qd = (double)def.call(x00)["qf"];
      double qe = (double)enabled.call(x00)["qf"];
      double qn = (double)disabled.call(x00)["qf"];
      assertEquals(qd, qe);
      assertTrue(Math.Abs(qd - exact) < 1e-4);
      assertTrue(Math.Abs(qn - exact) > 1e-2);
      assertTrue(num(def.stats()["n_steps"]) > num(disabled.stats()["n_steps"]));
    }
  }

  public void test_rosenbrock_order() {
    message("rosenbrock: order 5 on fixed meshes, ODE and index-1 DAE");
    SX t = SX.sym("t"), x = SX.sym("x", 2), z = SX.sym("z");
    var dae = new Dictionary<string, SX> { ["t"] = t, ["x"] = x, ["z"] = z,
      ["ode"] = ca.vertcat(-z + x[1] * 0.5, -x[0] * z + ca.sin(t)),
      ["alg"] = ca.pow(z, 3) + z - ca.pow(x[0], 2) * (1 + ca.cos(t) * 0.5), ["quad"] = z * x[1] };
    var a = new Dictionary<string, DM> { ["x0"] = ca.vertcat(1, 0.5), ["z0"] = 0.6 };
    Dictionary<string, DM> refr = ca.integrator("ref", "idas", dae, 0, 2, opts("reltol", 1e-14, "abstol", 1e-14)).call(a);
    var err = new List<double>();
    foreach (int N in new[] { 8, 16, 32 }) {
      Function F = ca.integrator("F", "rosenbrock", dae, 0, 2, opts("step_size", 2.0 / N));
      Dictionary<string, DM> r = F.call(a);
      double e = 0;
      foreach (string k in new[] { "xf", "zf", "qf" }) e = Math.Max(e, (double)ca.norm_inf(r[k] - refr[k]));
      err.Add(e);
    }
    for (int i = 0; i + 1 < err.Count; i++) assertTrue(Math.Log(err[i] / err[i + 1]) / Math.Log(2) > 4.6);
  }

  public void test_rosenbrock_ind() {
    message("rosenbrock: derivatives are exact derivatives of the scheme");
    SX t = SX.sym("t"), x = SX.sym("x", 2), z = SX.sym("z"), p = SX.sym("p", 2), u = SX.sym("u");
    var dae = new Dictionary<string, SX> { ["t"] = t, ["x"] = x, ["z"] = z, ["p"] = p, ["u"] = u,
      ["ode"] = ca.vertcat(-z + p[0] * x[1] + u, -x[0] * z + p[1] * ca.sin(t)),
      ["alg"] = ca.pow(z, 3) + z - ca.pow(x[0], 2) * (1 + ca.cos(t) * 0.5),
      ["quad"] = z * x[1] + p[1] * ca.pow(x[0], 2) };
    foreach (var o in new[] { opts("step_size", 0.1), opts("reltol", 1e-6, "abstol", 1e-6) })
      foreach (double[] grid in new[] { new double[] { 0.5, 1, 2 }, new double[] { 0, 1, 2 } }) {
        Function F = ca.integrator("F", "rosenbrock", dae, 0, grid, o);
        MX v = MX.sym("v", 7);
        Dictionary<string, MX> r = F.call(new Dictionary<string, MX> {
          ["x0"] = v[0..2], ["z0"] = 0.5, ["p"] = v[2..4], ["u"] = v[4..7].T() });
        MX o2 = ca.vertcat(ca.vec(r["xf"]), ca.vec(r["zf"]), ca.vec(r["qf"]));
        MX obj = ca.sumsqr(o2);
        DM v0 = ca.vertcat(1, 0.5, 0.5, 0.3, 0.1, -0.2, 0.3);
        checkJacobian(o2, v, v0, 12);
        DM H = checkHessian(obj, v, v0, 12);
        chk(H, H.T(), 12, "Hessian symmetry");
      }
  }

  public void test_rosenbrock_dense_output() {
    message("rosenbrock: outputs inside steps");
    SX t = SX.sym("t"), x = SX.sym("x", 2), z = SX.sym("z"), p = SX.sym("p", 2), u = SX.sym("u");
    var dae = new Dictionary<string, SX> { ["t"] = t, ["x"] = x, ["z"] = z, ["p"] = p, ["u"] = u,
      ["ode"] = ca.vertcat(-z + p[0] * x[1] + u, -x[0] * z + p[1] * ca.sin(t)),
      ["alg"] = ca.pow(z, 3) + z - ca.pow(x[0], 2) * (1 + ca.cos(t) * 0.5) - u,
      ["quad"] = z * x[1] + p[1] * ca.pow(x[0], 2) };
    var ode = new Dictionary<string, SX> { ["t"] = t, ["x"] = x, ["p"] = p, ["u"] = u,
      ["ode"] = ca.vertcat(x[1] + p[0] * u, -p[1] * x[0] + ca.sin(t)), ["quad"] = x[0] * ca.pow(x[1], 2) };
    foreach (var (sys, refPlugin) in new[] { (ode, "cvodes"), (dae, "idas") }) {
      if (!ca.has_integrator(refPlugin)) continue;
      var args = new Dictionary<string, DM> { ["x0"] = ca.vertcat(1, 0.5), ["p"] = ca.vertcat(0.5, 0.3) };
      if (sys.ContainsKey("z")) {
        // Consistent z0 for the reference: the real root of z^3 + z - 1.3 (numpy.roots in python)
        double zr = 1;
        for (int i = 0; i < 50; i++) zr -= (zr * zr * zr + zr - 1.3) / (3 * zr * zr + 1);
        args["z0"] = ca.DM(zr);
      }
      foreach (var o in new[] { opts("step_size", 0.04), opts("reltol", 1e-6, "abstol", 1e-6) }) {
        var nSteps = new List<double>();
        foreach (double[] grid0 in new[] { new[] { 2.0 }, linspace(0.002, 2, 1000) }) {
          Function F = ca.integrator("F", "rosenbrock", sys, 0, grid0,
              with(o, "control_group_size", new long[] { grid0.Length }));
          var a = new Dictionary<string, DM>(args) { ["u"] = ca.vertcat(0.2) };
          Dictionary<string, DM> r = F.call(a);
          nSteps.Add(num(F.stats()["n_steps"]));
          Dictionary<string, DM> refr = ca.integrator("ref", refPlugin, sys, 0, grid0,
              opts("reltol", 1e-10, "abstol", 1e-10)).call(a);
          foreach (string k in new[] { "xf", "zf", "qf" }) chk(refr[k], r[k], 5, k);
        }
        assertTrue(nSteps[1] <= nSteps[0] + 2);
        double[] grid = { 0.13, 0.31, 0.37, 0.52, 0.9, 1.15, 1.2, 1.71, 2.0 };
        foreach (double[] u0 in new[] { new[] { 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1 },
                                        new[] { 0.1, 0.1, 0.1, 0.1, 0.3, 0.3, 0.3, 0.3, 0.3 } }) {
          Function F = ca.integrator("F", "rosenbrock", sys, 0, grid, o);
          MX v = MX.sym("v", 13);
          Dictionary<string, MX> r = sys.ContainsKey("z")
              ? F.call(new Dictionary<string, MX> { ["x0"] = v[0..2], ["z0"] = 0.5, ["p"] = v[2..4], ["u"] = v[4..13].T() })
              : F.call(new Dictionary<string, MX> { ["x0"] = v[0..2], ["p"] = v[2..4], ["u"] = v[4..13].T() });
          MX o2 = ca.vertcat(ca.vec(r["xf"]), ca.vec(r["zf"]), ca.vec(r["qf"]));
          double[] v0 = new double[] { 1, 0.5, 0.5, 0.3 }.Concat(u0).ToArray();
          checkJacobian(o2, v, ca.DM(v0), 12);
          checkHessian(ca.sumsqr(o2), v, ca.DM(v0), 10);
        }
      }
    }
  }

  public void test_rosenbrock_codegen() {
    message("rosenbrock: code generation, including derivatives");
    SX x = SX.sym("x", 2), z = SX.sym("z"), p = SX.sym("p");
    var dae = new Dictionary<string, SX> { ["x"] = x, ["z"] = z, ["p"] = p, ["ode"] = ca.vertcat(x[1], -p * z),
      ["alg"] = ca.pow(z, 3) + z - x[0], ["quad"] = ca.pow(x[0], 2) };
    var cases = new (Dictionary<string, object>, double[])[] {
      (opts("reltol", 1e-8, "abstol", 1e-8, "max_tape_steps", 1000), new double[] { 0.5, 1 }),
      (opts("step_size", 0.05), new double[] { 0.5, 1 }), (opts("step_size", 0.03), linspace(0.07, 1, 15)),
      (opts("reltol", 1e-6, "abstol", 1e-6, "max_tape_steps", 1000), linspace(0.07, 1, 15)) };
    foreach (var (o, grid) in cases) {
      Function F = ca.integrator("F", "rosenbrock", dae, 0, grid, o);
      MX v = MX.sym("v", 3);
      Dictionary<string, MX> r = F.call(new Dictionary<string, MX> { ["x0"] = v[0..2], ["p"] = v[2] });
      MX obj = ca.sumsqr(r["xf"]) + ca.sumsqr(r["qf"]);
      var outs = new List<MX> { r["xf"], r["zf"], r["qf"], ca.gradient(obj, v) };
      if (grid.Length > 2) outs.Add(ca.jtimes(ca.jtimes(obj, v, MX.eye(3)).T(), v, MX.eye(3)));
      var G = ca.Function("G", new[] { v }, outs);
      checkCodegen(G, new List<DM> { ca.vertcat(0.4, -0.3, 2) }, 12);
    }
  }

  public void test_erk() {
    message("erk: order 4, dense output, steps independent of the outputs");
    SX t = SX.sym("t"), x = SX.sym("x", 2), p = SX.sym("p"), u = SX.sym("u");
    var dae = new Dictionary<string, SX> { ["t"] = t, ["x"] = x, ["p"] = p, ["u"] = u,
      ["ode"] = ca.vertcat(x[1], -p * x[0] + u * ca.cos(t)), ["quad"] = ca.pow(x[0], 2) };
    double w = 1.3;
    DM exact(double tt) => ca.vertcat(Math.Cos(w * tt), -w * Math.Sin(w * tt));
    var err = new List<double[]>();
    foreach (double h in new[] { 0.1, 0.05, 0.025 }) {
      Dictionary<string, DM> r0 = ca.integrator("F", "erk", dae, 0, new[] { 0.37, 2.0 }, opts("method", "rk4", "step_size", h))
          .call(new Dictionary<string, DM> { ["x0"] = ca.vertcat(1, 0), ["p"] = w * w });
      err.Add(new[] { (double)ca.norm_inf(r0["xf"][.., 1] - exact(2.0)), (double)ca.norm_inf(r0["xf"][.., 0] - exact(0.37)) });
    }
    for (int i = 0; i + 1 < err.Count; i++) {
      assertTrue(Math.Log(err[i][0] / err[i + 1][0]) / Math.Log(2) > 3.8);
      assertTrue(Math.Log(err[i][1] / err[i + 1][1]) / Math.Log(2) > 3.4);
    }
    var args = new Dictionary<string, DM> { ["x0"] = ca.vertcat(1, 0.2), ["p"] = 2, ["u"] = ca.DM(new double[,] { { 0.3, -0.1 } }) };
    Dictionary<string, DM> r = ca.integrator("F", "erk", dae, 0, new[] { 0.5, 1.0 }, opts("method", "rk4", "step_size", 0.1))
        .call(args);
    Dictionary<string, DM> refr = ca.integrator("F", "rk", dae, 0, new[] { 0.5, 1.0 }, opts("number_of_finite_elements", 10))
        .call(args);
    foreach (string k in new[] { "xf", "qf" }) chk(refr[k], r[k], 14, k);
    foreach (int n in new[] { 10, 1000 }) {
      Function F = ca.integrator("F", "erk", dae, 0, linspace(0, 10, n),
          opts("step_size", 0.05, "control_group_size", new long[] { n }));
      F.call(new Dictionary<string, DM> { ["x0"] = ca.vertcat(1, 0), ["p"] = w * w });
      assertEquals(200L, Convert.ToInt64(F.stats()["n_steps"]));
    }
  }

  public void test_erk_derivatives() {
    message("erk: derivatives are exact derivatives of the scheme");
    SX t = SX.sym("t"), x = SX.sym("x", 2), p = SX.sym("p"), u = SX.sym("u");
    var dae = new Dictionary<string, SX> { ["t"] = t, ["x"] = x, ["p"] = p, ["u"] = u,
      ["ode"] = ca.vertcat(x[1], -p * ca.pow(x[0], 3) + u * ca.cos(t)), ["quad"] = ca.pow(x[0], 2) + u * x[1] };
    double[] grid = { 0, 0.13, 0.31, 0.37, 1.0, 1.15, 2.0 };
    foreach (double[] u0 in new[] { new[] { 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1 }, new[] { 0.1, 0.1, 0.1, 0.3, 0.3, 0.3, 0.3 },
                                    linspace(0.1, 0.4, 7) }) {
      Function F = ca.integrator("F", "erk", dae, 0, grid, opts("step_size", 0.1));
      MX v = MX.sym("v", 10);
      Dictionary<string, MX> r = F.call(new Dictionary<string, MX> { ["x0"] = v[0..2], ["p"] = v[2], ["u"] = v[3..10].T() });
      MX o = ca.vertcat(ca.vec(r["xf"]), ca.vec(r["qf"]));
      DM v0 = ca.vertcat(new[] { 1, 0.5, 1.7 }.Concat(u0).ToArray());
      checkJacobian(o, v, v0, 13);
      MX obj = ca.sumsqr(o);
      checkHessian(obj, v, v0, 12);
      DM g = ca.Function("g", new[] { v }, new[] { ca.gradient(obj, v) }).call(new[] { v0 })[0];
      chk(g, fdJacobian(ca.Function("f", new[] { v }, new[] { obj }), v0).T(), 7, "FD");
    }
  }

  public void test_erk_codegen() {
    message("erk: code generation, including derivatives");
    SX x = SX.sym("x", 2), p = SX.sym("p"), u = SX.sym("u");
    var dae = new Dictionary<string, SX> { ["x"] = x, ["p"] = p, ["u"] = u, ["ode"] = ca.vertcat(x[1], -p * ca.sin(x[0]) + u),
      ["quad"] = ca.pow(x[0], 2) };
    var cases = new (double[], Dictionary<string, object>)[] {
      (new double[] { 0.5, 1 }, opts("step_size", 0.04)), (linspace(0.07, 1, 15), opts("step_size", 0.04)),
      (linspace(0.07, 1, 15), opts("method", "rk4", "step_size", 0.1)),
      (linspace(0.07, 1, 15), opts("method", "rk45", "reltol", 1e-6, "abstol", 1e-6, "max_tape_steps", 200)) };
    foreach (var (grid, o) in cases) {
      Function F = ca.integrator("F", "erk", dae, 0, grid, o);
      int n = 3 + grid.Length;
      MX v = MX.sym("v", n);
      Dictionary<string, MX> r = F.call(new Dictionary<string, MX> { ["x0"] = v[0..2], ["p"] = v[2], ["u"] = v[3..n].T() });
      MX obj = ca.sumsqr(r["xf"]) + ca.sumsqr(r["qf"]);
      MX S = ca.MX(DM.eye(n).get(false, ca.Slice(), ca.vertcat(0, 2, 3)));
      var outs = new[] { r["xf"], r["qf"], ca.gradient(obj, v), ca.jtimes(ca.jtimes(obj, v, S).T(), v, S) };
      var G = ca.Function("G", new[] { v }, outs);
      double[] v0 = new double[n];
      v0[0] = 0.4;
      v0[1] = -0.3;
      v0[2] = 2;
      for (int i = 0; i < grid.Length; i++) v0[3 + i] = 0.1 * (i % 3);
      checkCodegen(G, new List<DM> { ca.DM(v0) }, 12);
    }
  }

  public void test_control_groups() {
    message("erk, rosenbrock: one control per control group, derivatives for each control");
    SX t = SX.sym("t"), x = SX.sym("x", 2), z = SX.sym("z"), p = SX.sym("p"), u = SX.sym("u");
    var ode = new Dictionary<string, SX> { ["t"] = t, ["x"] = x, ["p"] = p, ["u"] = u,
      ["ode"] = ca.vertcat(x[1], -p * ca.pow(x[0], 3) + u * ca.cos(t)), ["quad"] = ca.pow(x[0], 2) + u * x[1] };
    var dae = new Dictionary<string, SX> { ["t"] = t, ["x"] = x, ["z"] = z, ["p"] = p, ["u"] = u,
      ["ode"] = ca.vertcat(x[1], -p * x[0] + u * z), ["alg"] = ca.pow(z, 3) + z - u - x[0] * 0.5, ["quad"] = ca.pow(x[0], 2) + u * z };
    double[] grid = { 0.2, 0.3, 0.5, 1.0, 1.2, 2.0 };
    var cases = new (string, Dictionary<string, object>)[] {
      ("erk", opts("method", "rk4", "step_size", 0.1)), ("rosenbrock", opts("step_size", 0.1)),
      ("erk", opts("reltol", 1e-10, "abstol", 1e-10, "max_tape_steps", 500)),
      ("rosenbrock", opts("reltol", 1e-10, "abstol", 1e-10, "max_tape_steps", 500)) };
    foreach (var (plugin, o) in cases) {
      var d = plugin == "rosenbrock" ? dae : ode;
      Dictionary<string, DM> r = ca.integrator("F", plugin, d, 0, grid, o).call(new Dictionary<string, DM> {
        ["u"] = ca.DM(new double[,] { { 0.1, 0.1, 0.3, 0.2, 0.2, 0.2 } }), ["x0"] = ca.vertcat(1, 0.5), ["p"] = 1.7 });
      Function F = ca.integrator("F", plugin, d, 0, grid, with(o, "control_group_size", new long[] { 2, 1, 3 }));
      Dictionary<string, DM> rGroup = F.call(new Dictionary<string, DM> {
        ["u"] = ca.DM(new double[,] { { 0.1, 0.3, 0.2 } }), ["x0"] = ca.vertcat(1, 0.5), ["p"] = 1.7 });
      foreach (string k in new[] { "xf", "zf", "qf" }) chk(r[k], rGroup[k], 6, k);
      foreach (long[] cgs in new[] { null, new long[] { 2, 1, 3 } }) {
        long[] groups = cgs ?? new long[] { 1, 1, 1, 1, 1, 1 };
        F = ca.integrator("F", plugin, d, 0, grid, with(o, "control_group_size", groups));
        int n = 3 + (int)F.size2_in(3);
        MX v = MX.sym("v", n);
        Dictionary<string, MX> rv = F.call(new Dictionary<string, MX> { ["x0"] = v[0..2], ["p"] = v[2], ["u"] = v[3..n].T() });
        MX outv = ca.vertcat(ca.vec(rv["xf"]), ca.vec(rv["zf"]), ca.vec(rv["qf"]));
        var G = ca.Function("G", new[] { v }, new[] { outv, ca.jtimes(outv, v, MX.eye(n)),
            ca.jtimes(outv, v, MX.eye(outv.numel()), true).T() });
        double[] v0a = new double[n];
        v0a[0] = 1;
        v0a[1] = 0.5;
        v0a[2] = 1.7;
        for (int i = 3; i < n; i++) v0a[i] = 0.2;
        DM v0 = ca.DM(v0a);
        List<DM> Go = G.call(new[] { v0 });
        DM Jfd = fdJacobian(G, v0);
        chk(Jfd, Go[1], 7, "fwd vs FD");
        chk(Jfd, Go[2], 7, "adj vs FD");
        if (cgs != null) checkCodegen(G, new List<DM> { v0 }, 12);
        long[] start = new long[groups.Length + 1];
        for (int i = 0; i < groups.Length; i++) start[i + 1] = start[i] + groups[i];
        foreach (int wgt in new[] { 0, 1 }) {
          Sparsity J = ca.Function("J", new[] { v }, new[] { outv }, opts("ad_weight_sp", wgt)).jac_sparsity(0, 0);
          double[,] dep = ca.DM(J, 1.0).ToArray();
          for (int k = 0; k < 6; k++)
            for (int j = 0; j < n - 3; j++) {
              bool any = dep[2 * k, 3 + j] != 0 || dep[2 * k + 1, 3 + j] != 0;
              assertEquals(start[j] <= k, any);
            }
        }
      }
    }
    assertInException("must sum to the number of output times",
        () => ca.integrator("F", "erk", ode, 0, grid, opts("control_group_size", new long[] { 2, 1 })));
  }

  public void test_erk_rk45() {
    message("erk: pair of order 5(4), fixed and adaptive steps");
    SX t = SX.sym("t"), x = SX.sym("x", 2), p = SX.sym("p");
    var dae = new Dictionary<string, SX> { ["t"] = t, ["x"] = x, ["p"] = p, ["ode"] = ca.vertcat(x[1], -p * x[0]),
      ["quad"] = ca.pow(x[0], 2) };
    double w = 1.3;
    double[] exact(double tt) => new[] { Math.Cos(w * tt), -w * Math.Sin(w * tt) };
    var err = new List<double>();
    foreach (double h in new[] { 0.1, 0.05, 0.025 }) {
      List<double> xf = ca.integrator("F", "erk", dae, 0, 2, opts("method", "rk45", "step_size", h))
          .call(new Dictionary<string, DM> { ["x0"] = ca.vertcat(1, 0), ["p"] = w * w })["xf"].nonzeros();
      double[] ex = exact(2.0);
      err.Add(Math.Max(Math.Abs(xf[0] - ex[0]), Math.Abs(xf[1] - ex[1])));
    }
    for (int i = 0; i + 1 < err.Count; i++) assertTrue(Math.Log(err[i] / err[i + 1]) / Math.Log(2) > 4.7);
    double[] grid = linspace(0.1, 10, 50);
    foreach (double tol in new[] { 1e-6, 1e-9 }) {
      Function F = ca.integrator("F", "erk", dae, 0, grid, opts("method", "rk45", "reltol", tol, "abstol", tol));
      double[,] xf = F.call(new Dictionary<string, DM> { ["x0"] = ca.vertcat(1, 0), ["p"] = w * w })["xf"].ToArray();
      double e = 0;
      for (int j = 0; j < grid.Length; j++) {
        double[] ex = exact(grid[j]);
        for (int i = 0; i < 2; i++) e = Math.Max(e, Math.Abs(xf[i, j] - ex[i]));
      }
      assertTrue(e < 10 * tol);
      assertTrue(e > 0.01 * tol);
      Dictionary<string, object> st = F.stats();
      assertEquals(1 + 6 * (num(st["n_steps"]) + num(st["n_reject"])), num(st["n_rhs"]));
    }
    var dae2 = new Dictionary<string, SX> { ["x"] = x, ["ode"] = ca.vertcat(-50 * (x[0] - ca.cos(x[1])), 1) };
    Function F2 = ca.integrator("F", "erk", dae2, 0, 1, opts("method", "rk45", "step0", 0.5));
    Dictionary<string, DM> r = F2.call(new Dictionary<string, DM> { ["x0"] = ca.vertcat(0, 0) });
    Dictionary<string, DM> refr = ca.integrator("F", "cvodes", dae2, 0, 1, opts("reltol", 1e-12, "abstol", 1e-12))
        .call(new Dictionary<string, DM> { ["x0"] = ca.vertcat(0, 0) });
    chk(refr["xf"], r["xf"], 5, "rejections");
    assertTrue(num(F2.stats()["n_reject"]) > 0);
    assertInException("unknown method", () => ca.integrator("F", "erk", dae, 0, 1, opts("method", "rk5")));
  }

  public void test_erk_rk45_derivatives() {
    message("erk: adaptive steps, derivatives of the scheme on the nominal mesh");
    SX t = SX.sym("t"), x = SX.sym("x", 2), p = SX.sym("p"), u = SX.sym("u");
    var dae = new Dictionary<string, SX> { ["t"] = t, ["x"] = x, ["p"] = p, ["u"] = u,
      ["ode"] = ca.vertcat(x[1], -p * ca.pow(x[0], 3) + u * ca.cos(t)), ["quad"] = ca.pow(x[0], 2) + u * x[1] };
    double[] grid = { 0, 0.13, 0.31, 0.37, 1.0, 1.15, 2.0 };
    foreach (double[] u0 in new[] { new[] { 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1 }, new[] { 0.1, 0.1, 0.1, 0.3, 0.3, 0.3, 0.3 } }) {
      Function F = ca.integrator("F", "erk", dae, 0, grid, opts("method", "rk45", "reltol", 1e-5, "abstol", 1e-5));
      MX v = MX.sym("v", 10);
      Dictionary<string, MX> r = F.call(new Dictionary<string, MX> { ["x0"] = v[0..2], ["p"] = v[2], ["u"] = v[3..10].T() });
      MX o = ca.vertcat(ca.vec(r["xf"]), ca.vec(r["qf"]));
      double[] v0a = new[] { 1, 0.5, 1.7 }.Concat(u0).ToArray();
      checkJacobian(o, v, ca.DM(v0a), 13);
      checkHessian(ca.sumsqr(o), v, ca.DM(v0a), 11);
    }
  }

  const double G_ = 9.81;

  /// <summary>Height, velocity and integral of the height of the bouncing ball</summary>
  static double[] exactBounce(double h0, double p0, double tt) {
    double g = G_, tb = 0, q = 0, h = h0, v = 0;
    double? ub = null;
    while (true) {
      double tau = ub == null ? Math.Sqrt(2 * h / g) : 2 * ub.Value / g;
      if (tb + tau > tt) break;
      q += ub == null ? h * tau - g * Math.Pow(tau, 3) / 6 : ub.Value * tau * tau / 2 - g * Math.Pow(tau, 3) / 6;
      v = ub == null ? -g * tau : ub.Value - g * tau;
      tb = tb + tau;
      ub = -p0 * v;
      h = 0.0;
    }
    double tau2 = tt - tb, u0 = ub ?? 0.0;
    return new[] { h + u0 * tau2 - g * tau2 * tau2 / 2, u0 - g * tau2,
                   q + h * tau2 + u0 * tau2 * tau2 / 2 - g * Math.Pow(tau2, 3) / 6 };
  }

  public void test_events_bouncing_ball() {
    message("zero-crossing events: values, forward sensitivities, code generation");
    SX x = SX.sym("x", 2), p = SX.sym("p");
    var dae = new Dictionary<string, SX> { ["x"] = x, ["p"] = p, ["ode"] = ca.vertcat(x[1], -G_), ["quad"] = x[0],
                                           ["zero"] = x[0] };
    Function transition = bounceTransition();
    double[] grid = linspace(0.1, 4, 40);
    double h0 = 2.0, p0 = 0.8;
    double[,] xfRef = new double[2, grid.Length];
    double[] qfRef = new double[grid.Length];
    for (int i = 0; i < grid.Length; i++) {
      double[] e = exactBounce(h0, p0, grid[i]);
      xfRef[0, i] = e[0];
      xfRef[1, i] = e[1];
      qfRef[i] = e[2];
    }
    foreach (EvIntg E in eventIntegrators) {
      if (!E.events) continue;
      Function F = ca.integrator("F", E.name, dae, 0, grid, with(E.options, "transition", transition, "event_tol", 1e-12));
      Dictionary<string, DM> r = F.call(new Dictionary<string, DM> { ["x0"] = ca.vertcat(h0, 0), ["p"] = p0 });
      chk(r["xf"], xfRef, E.digits, E.name);
      chk(r["qf"], ca.DM(qfRef).T(), E.digits, E.name);
      if (F.stats().ContainsKey("n_events")) assertEquals(5L, Convert.ToInt64(F.stats()["n_events"]));
      MX v = MX.sym("v", 3);
      Dictionary<string, MX> rv = F.call(new Dictionary<string, MX> { ["x0"] = v[0..2], ["p"] = v[2] });
      MX o = ca.vertcat(rv["xf"][.., -1], rv["qf"][-1]);
      var G = ca.Function("G", new[] { v }, new[] { o, ca.jacobian(o, v) });
      DM v0 = ca.vertcat(h0, 0.3, p0);
      chk(G.call(new[] { v0 })[1], fdJacobian(G, v0), Math.Min(6, E.digits), E.name);
      if (E.codegen) checkCodegen(G, new List<DM> { v0 }, 10);
    }
  }

  public void test_events_multiple() {
    message("zero-crossing events: several functions, transition by event index");
    double g = G_, ceiling = 1.5;
    SX x = SX.sym("x", 2), p = SX.sym("p"), index = SX.sym("index"), t = SX.sym("t");
    var dae = new Dictionary<string, SX> { ["x"] = x, ["p"] = p, ["ode"] = ca.vertcat(x[1], -g),
                                           ["zero"] = ca.vertcat(x[0], ceiling - x[0]) };
    var transition = ca.Function("transition", new[] { index, t, x, ca.SX(0, 1), p, ca.SX(0, 1) },
        new[] { ca.vertcat(x[0], ca.if_else(ca.eq(index, 0), -p, -1) * x[1]), ca.SX(0, 1) }, transitionIn, transitionOut);
    DM v0 = ca.vertcat(1, 4, 0.8);
    List<double> v0a = v0.nonzeros();
    double[] grid = linspace(0.1, 3, 30);
    double[,] xr = new double[2, grid.Length];
    for (int i = 0; i < grid.Length; i++) {
      double tt = grid[i], t0 = 0, h = v0a[0], v = v0a[1];
      while (true) {
        double bestTau = double.PositiveInfinity, bestLevel = 0;
        foreach (var (level, sign) in new[] { (0.0, -1.0), (ceiling, 1.0) }) {
          double disc = v * v - 2 * g * (level - h);
          if (disc < 0) continue;
          foreach (double tau in new[] { (v - Math.Sqrt(disc)) / g, (v + Math.Sqrt(disc)) / g })
            if (tau > 1e-12 && sign * (v - g * tau) > 0 && (tau < bestTau || (tau == bestTau && level < bestLevel))) {
              bestTau = tau;
              bestLevel = level;
            }
        }
        if (t0 + bestTau > tt) break;
        v = (v - g * bestTau) * (bestLevel == 0 ? -v0a[2] : -1);
        t0 = t0 + bestTau;
        h = bestLevel;
      }
      double tau2 = tt - t0;
      xr[0, i] = h + v * tau2 - g * tau2 * tau2 / 2;
      xr[1, i] = v - g * tau2;
    }
    DM xfRef = xr;
    foreach (EvIntg E in eventIntegrators) {
      if (!E.events) continue;
      Function F = ca.integrator("F", E.name, dae, 0, grid, with(E.options, "transition", transition, "event_tol", 1e-12));
      MX v = MX.sym("v", 3);
      MX xf = F.call(new Dictionary<string, MX> { ["x0"] = v[0..2], ["p"] = v[2] })["xf"];
      var G = ca.Function("G", new[] { v }, new[] { ca.vec(xf), ca.jacobian(ca.vec(xf), v) });
      chk(G.call(new[] { v0 })[0], ca.vec(xfRef), E.digits, E.name);
      chk(G.call(new[] { v0 })[1], fdJacobian(G, v0), Math.Min(6, E.digits), E.name);
      if (E.codegen) checkCodegen(G, new List<DM> { v0 }, 10);
    }
  }

  public void test_events_dae() {
    message("zero-crossing events of a DAE, algebraic states consistent after the transition");
    double g = G_;
    SX x = SX.sym("x", 2), z = SX.sym("z"), p = SX.sym("p");
    var dae = new Dictionary<string, SX> { ["x"] = x, ["z"] = z, ["p"] = p, ["ode"] = ca.vertcat(z, -g), ["alg"] = z - x[1],
                                           ["zero"] = x[0] };
    Function transition = bounceTransition(1);
    double[] grid = linspace(0.1, 3, 30);
    Function refF = ca.integrator("ref", "erk",
        new Dictionary<string, SX> { ["x"] = x, ["p"] = p, ["ode"] = ca.vertcat(x[1], -g), ["zero"] = x[0] }, 0, grid,
        opts("transition", bounceTransition(), "event_tol", 1e-12, "abstol", 1e-11, "reltol", 1e-11));
    DM v0 = ca.vertcat(2, 0.3, 0.8);
    DM xfRef = refF.call(new Dictionary<string, DM> { ["x0"] = v0[0..2], ["p"] = v0[2] })["xf"];
    foreach (EvIntg E in eventIntegrators) {
      if (!E.eventsDae) continue;
      Function F = ca.integrator("F", E.name, dae, 0, grid, with(E.options, "transition", transition, "event_tol", 1e-12));
      MX v = MX.sym("v", 3);
      Dictionary<string, MX> r = F.call(new Dictionary<string, MX> { ["x0"] = v[0..2], ["z0"] = v[1], ["p"] = v[2] });
      MX o = ca.vertcat(ca.vec(r["xf"]), ca.vec(r["zf"]));
      var G = ca.Function("G", new[] { v }, new[] { o, ca.jacobian(o, v) });
      chk(G.call(new[] { v0 })[0], ca.vertcat(ca.vec(xfRef), xfRef[1, ..].T()), E.digits, E.name);
      chk(G.call(new[] { v0 })[1], fdJacobian(G, v0), Math.Min(6, E.digits), E.name);
      if (E.codegen) checkCodegen(G, new List<DM> { v0 }, 10);
    }
  }

  public void test_events_breaking_spring() {
    message("zero-crossing events with controls and parameters");
    SX x = SX.sym("x", 3), p = SX.sym("p", 3), u = SX.sym("u"), index = SX.sym("index"), t = SX.sym("t");
    SX f = ca.if_else(x[2] < 0.5, -p[1] * x[1] + u, 0);
    var dae = new Dictionary<string, SX> { ["x"] = x, ["p"] = p, ["u"] = u,
      ["ode"] = ca.vertcat((f - p[2] * x[0]) / p[0], x[0], 0), ["zero"] = 2 - x[1] };
    var transition = ca.Function("transition", new[] { index, t, x, ca.SX(0, 1), p, u },
        new[] { ca.vertcat(x[0], x[1], 1), ca.SX(0, 1) }, transitionIn, transitionOut);
    double[] grid = linspace(0.04, 4, 100);
    double[] d = new double[100];
    for (int i = 20; i < 30; i++) d[i] = -5;
    for (int i = 40; i < 50; i++) d[i] = 15;
    DM x0 = ca.vertcat(-5, -1, 0), p0 = ca.vertcat(1, 2, 0.1), dd = ca.DM(d);
    MX v = MX.sym("v", 3);
    Function outputs(Function F) {
      MX xf = F.call(new Dictionary<string, MX> { ["x0"] = x0, ["p"] = v, ["u"] = dd })["xf"];
      return ca.Function("G", new[] { v }, new[] { xf, ca.jacobian(xf, v) });
    }
    List<DM> refr = outputs(ca.integrator("ref", "erk", dae, 0, grid, opts("transition", transition,
        "event_tol", 1e-12, "abstol", 1e-12, "reltol", 1e-12))).call(new[] { p0 });
    var call0 = new Dictionary<string, DM> { ["x0"] = x0, ["p"] = p0, ["u"] = dd };
    foreach (EvIntg E in eventIntegrators) {
      if (!E.events) continue;
      Function F = ca.integrator("F", E.name, dae, 0, grid, with(E.options, "transition", transition, "event_tol", 1e-12));
      Function G = outputs(F);
      List<DM> r = G.call(new[] { p0 });
      chk(r[0], refr[0], E.digits - 1, E.name);
      chk(r[1], refr[1], E.digits - 2, E.name);
      F.call(call0);
      if (F.stats().ContainsKey("n_events")) assertEquals(1L, Convert.ToInt64(F.stats()["n_events"]));
      if (E.codegen) checkCodegen(G, new List<DM> { p0 }, 12);
    }
    Function F0 = ca.integrator("F", "erk", dae, 0, grid, opts("transition", transition, "max_events", 0));
    assertInException("max_events", () => F0.call(call0));
  }

  public void test_events_adjoint() {
    message("zero-crossing events: adjoint and forward-over-adjoint sensitivities");
    double g = G_;
    SX x = SX.sym("x", 2), z = SX.sym("z"), p = SX.sym("p");
    var ode = new Dictionary<string, SX> { ["x"] = x, ["p"] = p, ["ode"] = ca.vertcat(x[1], -g), ["quad"] = x[0],
                                           ["zero"] = x[0] };
    var dae = new Dictionary<string, SX> { ["x"] = x, ["z"] = z, ["p"] = p, ["ode"] = ca.vertcat(z, -g), ["alg"] = z - x[1],
                                           ["quad"] = x[0], ["zero"] = x[0] };
    double[] grid = linspace(0.1, 3, 30);
    // numpy.random.default_rng(0).uniform(0.5, 1.5, (2, 30)) in python: System.Random here
    var rng = new Random(0);
    double[,] wa = new double[2, 30];
    for (int i = 0; i < 2; i++)
      for (int j = 0; j < 30; j++) wa[i, j] = 0.5 + rng.NextDouble();
    DM wts = wa;
    foreach (EvIntg E in eventIntegrators) {
      if (!E.eventsAdjoint) continue;
      var problems = new List<(Dictionary<string, SX>, Function)> { (ode, bounceTransition()) };
      if (E.eventsDae) problems.Add((dae, bounceTransition(1)));
      foreach (var (sys, transition) in problems) {
        Function F = ca.integrator("F", E.name, sys, 0, grid, with(E.options, "transition", transition, "event_tol", 1e-12));
        MX v = MX.sym("v", 3);
        Dictionary<string, MX> r = F.call(new Dictionary<string, MX> { ["x0"] = v[0..2], ["z0"] = v[1], ["p"] = v[2] });
        MX obj = ca.sum1(ca.sum2(wts * r["xf"])) + ca.sumsqr(r["qf"]);
        var G = ca.Function("G", new[] { v }, new[] { ca.jtimes(obj, v, MX.ones(1, 1), true),
            ca.jtimes(obj, v, MX.eye(3)).T(), ca.hessian(obj, v).Item1 });
        DM v0 = ca.vertcat(2, 0.3, 0.8);
        List<DM> Go = G.call(new[] { v0 });
        chk(Go[0], Go[1], 9, E.name);
        DM fd = (fdJacobian(G, v0, 1e-4) * 4 - fdJacobian(G, v0, 2e-4)) / 3;
        chk(Go[2], fd, 6, E.name);
        if (E.codegen) checkCodegen(G, new List<DM> { v0 }, 9);
      }
    }
  }

  public void test_erk_tape() {
    message("erk: adaptive steps and adjoints, bounded tape");
    SX x = SX.sym("x", 2);
    var dae = new Dictionary<string, SX> { ["x"] = x, ["ode"] = ca.vertcat(x[1], -x[0]) };
    foreach (int nTape in new[] { 0, 3, 1000 }) {
      Function F = ca.integrator("F", "erk", dae, 0, 10, opts("method", "rk45", "max_tape_steps", nTape));
      MX x0 = MX.sym("x0", 2);
      var G = ca.Function("G", new[] { x0 }, new[] { ca.jtimes(ca.sumsqr(F.call(new Dictionary<string, MX> { ["x0"] = x0 })["xf"]),
          x0, MX.ones(1, 1), true) });
      if (nTape == 3) {
        assertInException("increase max_tape_steps", () => G.call(new[] { ca.vertcat(1, 0.5) }));
      } else {
        chk(G.call(new[] { ca.vertcat(1, 0.5) })[0], ca.vertcat(2, 1), 5, "adjoint");
      }
      if (nTape == 0)
        assertInException("need max_tape_steps for code generation", () => ca.CodeGenerator("tape_cg.c").add(G));
    }
  }
}

}

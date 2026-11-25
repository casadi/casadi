using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;
using System.Runtime.InteropServices;

namespace CasADi.Tests {

/// <summary>Port of test/java OcpTest (itself a port of test/python/ocp.py)</summary>
public class OcpTest : CasadiTestCase {
  const double inf = double.PositiveInfinity;

  static void assertAlmostEqual(double a, double b, int places, string msg = "") =>
    assertTrue(Math.Abs(a - b) < 0.5 * Math.Pow(10, -places), msg + ": " + a + " vs " + b);

  static bool[] rep(bool v, long n) => Enumerable.Repeat(v, (int)n).ToArray();

  static DM proj(DM a, IDictionary<string, Sparsity> sp, string k) => sp.ContainsKey(k) ? ca.project(a, sp[k]) : a;

  static string str(IDictionary<string, object> d) => "{" + string.Join(", ", d.Select(kv => kv.Key + "=" + kv.Value)) + "}";

  void fatrop_case(long N = 2, long nx0 = 2, long nu0 = 2, long nx1 = 2, long nu1 = 2, long nx2 = 2, long nu2 = 2,
                   long ng1 = 2, long ng2 = 2, long ng3 = 2, IDictionary<string, Sparsity> sp = null,
                   ICollection<string> eq = null) {
    sp ??= new Dictionary<string, Sparsity>();
    eq ??= new HashSet<string>();
    Console.WriteLine($"fatrop_case N={N} nx0={nx0} nu0={nu0} nx1={nx1} nu1={nu1} nx2={nx2} nu2={nu2} ng1={ng1} ng2={ng2} ng3={ng3}");
    long[] nx = { nx0, nx1, nx2 }, nu = { nu0, nu1, nu2 }, ng = { ng1, ng2, ng3 };

    DM.rng(1);

    DM A0 = proj(DM.rand(nx1, nx0), sp, "A0");
    DM B0 = proj(DM.rand(nx1, nu0), sp, "B0");
    DM C0 = proj(DM.rand(ng1, nx0), sp, "C0");
    DM D0 = DM.rand(ng1, nu0);
    DM I0 = DM.eye(nx1);

    DM A1 = proj(DM.rand(nx2, nx1), sp, "A1");
    DM B1 = proj(DM.rand(nx2, nu1), sp, "B1");
    DM C1 = proj(DM.rand(ng2, nx1), sp, "C1");
    DM D1 = proj(DM.rand(ng2, nu1), sp, "D1");
    DM I1 = DM.eye(nx2);

    DM C2 = proj(DM.rand(ng3, nx2), sp, "C2");
    DM D2 = proj(DM.rand(ng3, nu2), sp, "D2");

    DM A = ca.blockcat(new[] {
        new[] { A0, B0, I0, ca.DM(nx1, nu1 + nx2 + nu2) },
        new[] { C0, D0, ca.DM(ng1, nx1 + nu1 + nx2 + nu2) },
        new[] { ca.DM(nx2, nx0 + nu0), A1, B1, I1, ca.DM(nx2, nu2) },
        new[] { ca.DM(ng2, nx0 + nu0), C1, D1, ca.DM(ng2, nx2 + nu2) },
        new[] { ca.DM(ng3, nx0 + nu0 + nx1 + nu1), C2, D2 } });

    bool[] equality = rep(true, nx1).Concat(rep(eq.Contains("ng1"), ng1)).Concat(rep(true, nx2))
        .Concat(rep(eq.Contains("ng2"), ng2)).Concat(rep(eq.Contains("ng3"), ng3)).ToArray();

    A.sparsity().spy();
    Console.WriteLine(A);

    MX x = ca.vertcat(MX.sym("x0", nx0), MX.sym("u0", nu0), MX.sym("x1", nx1), MX.sym("u1", nu1),
        MX.sym("x2", nx2), MX.sym("u2", nu2));
    var nlp = new Dictionary<string, MX> {
      ["x"] = x,
      ["g"] = DM.zeros(A.size1(), 1) + ca.mtimes((MX)A, x),
      ["f"] = ca.sumsqr(x - DM.rand(x.numel(), 1)) };

    double a = 10;
    DM lbg = ca.vertcat(DM.zeros(nx1, 1), DM.ones(ng1, 1) * -a, DM.zeros(nx2, 1), DM.ones(ng2, 1) * -a, DM.ones(ng3, 1) * -a);
    DM ubg = ca.vertcat(DM.zeros(nx1, 1), DM.ones(ng1, 1) * a, DM.zeros(nx2, 1), DM.ones(ng2, 1) * a, DM.ones(ng3, 1) * a);
    Console.WriteLine(lbg);

    var options = new Dictionary<string, object> {
      ["structure_detection"] = "manual",
      ["N"] = N,
      ["nx"] = nx,
      ["nu"] = nu,
      ["ng"] = ng,
      ["equality"] = equality,
      ["fatrop"] = new Dictionary<string, object> { ["tol"] = 1e-7 } };
    Function solver = ca.nlpsol("solver", "fatrop", nlp, options);
    var sol = solver.call(new Dictionary<string, DM> { ["lbg"] = lbg, ["ubg"] = ubg });

    solver = ca.nlpsol("solver", "fatrop", nlp, new Dictionary<string, object> {
      ["structure_detection"] = "none", ["error_on_fail"] = true, ["equality"] = equality,
      ["fatrop"] = new Dictionary<string, object> { ["tol"] = 1e-7 } });
    var @ref = solver.call(new Dictionary<string, DM> { ["lbg"] = lbg, ["ubg"] = ubg });

    foreach (string k in sol.Keys) checkarray(sol[k], @ref[k], 1e-6, k + str(options));

    options = new Dictionary<string, object> {
      ["structure_detection"] = "auto", ["debug"] = true, ["equality"] = equality,
      ["fatrop"] = new Dictionary<string, object> { ["tol"] = 1e-7 } };
    Console.WriteLine(str(options));
    solver = ca.nlpsol("solver", "fatrop", nlp, options);
    sol = solver.call(new Dictionary<string, DM> { ["lbg"] = lbg, ["ubg"] = ubg });

    var stats = solver.stats();
    Console.WriteLine(str(stats));
    if (nx2 > 0) assertTrue(Convert.ToInt64(stats["N"]) > 1);

    foreach (string k in sol.Keys) checkarray(sol[k], @ref[k], 1e-6, k + str(options));
  }

  public void testdiscrete() {
    Check.Assume(hasNlpsol("ipopt"), "nlpsol ipopt not available");
    message("Linear-quadratic problem, discrete, using IPOPT");
    double a = 1.0, b = 1.0, q = 1.0, s = 1.0, r = 1.0;
    double x0 = 100;

    int N = 100;

    SX X = SX.sym("X", N + 1);
    SX U = SX.sym("U", N);

    SX V = ca.vertcat(X, U);

    SX cost = 0.0;
    for (int i = 0; i < N; i++) cost = cost + ca.pow(X[i], 2) * s + ca.pow(U[i], 2) * r;
    cost = cost + ca.pow(X[N], 2) * q;

    var nlp = new Dictionary<string, SX> {
      ["x"] = V, ["f"] = cost, ["g"] = ca.vertcat(X[0] - x0, X[1..(N + 1), 0..1] - (X[0..N, 0..1] * a + U * b)) };
    var opts = new Dictionary<string, object> {
      ["ipopt.tol"] = 1e-5, ["ipopt.hessian_approximation"] = "limited-memory", ["ipopt.max_iter"] = 100,
      ["ipopt.print_level"] = 0 };
    Function solver = ca.nlpsol("solver", "ipopt", nlp, opts);
    double[] lbx = Enumerable.Repeat(-1000.0, (int)V.nnz()).ToArray(), ubx = Enumerable.Repeat(1000.0, (int)V.nnz()).ToArray();
    var solver_out = solver.call(new Dictionary<string, DM> {
      ["lbx"] = lbx, ["ubx"] = ubx, ["lbg"] = new double[N + 1], ["ubg"] = new double[N + 1] });
    double ocp_sol = (double)solver_out["f"][0];
    double K = q + 0.0;
    for (int i = 0; i < N; i++) K = s + r * a * a * K / (r + b * b * K);
    double exact_sol = K * x0 * x0;
    assertAlmostEqual(ocp_sol, exact_sol, 10, "Linear-quadratic problem solution using IPOPT");
  }

  static double dfdp(double p0, double yc0, double te) =>
    (Math.Sqrt(p0) * (te * yc0 * yc0 - yc0 + p0 * te) * Math.Tan(Math.Atan(yc0 / Math.Sqrt(p0)) + Math.Sqrt(p0) * te)
        + yc0 * yc0) / (2 * p0 * yc0 * yc0 + 2 * p0 * p0);

  static double cost(double p0, double y0, double yc0, double te) =>
    (2 * y0 - Math.Log(yc0 * yc0 / p0 + 1)) / 2 - Math.Log(Math.Cos(Math.Atan(yc0 / Math.Sqrt(p0)) + Math.Sqrt(p0) * te));

  static Function shootingIntegrator(double te) {
    SX t = SX.sym("t");
    SX q = SX.sym("y", 2, 1);
    SX p = SX.sym("p", 1, 1);
    var dae = new Dictionary<string, SX> { ["x"] = q, ["p"] = p, ["t"] = t, ["ode"] = ca.vertcat(q[1], p[0] + ca.pow(q[1], 2)) };
    var opts = new Dictionary<string, object> {
      ["reltol"] = 1e-15, ["abstol"] = 1e-15, ["verbose"] = false, ["steps_per_checkpoint"] = 10000 };
    return ca.integrator("integrator", "cvodes", dae, 0, te, opts);
  }

  public void test_singleshooting() {
    Check.Assume(hasNlpsol("ipopt"), "nlpsol ipopt not available");
    message("Single shooting");
    double p0 = 0.2, y0 = 1, yc0 = 0, te = 0.4;

    Function integr = shootingIntegrator(te);

    MX var = MX.sym("var", 2, 1);
    MX par = MX.sym("par", 1, 1);
    MX parMX = par;

    MX q0 = ca.vertcat(var[0], par);
    par = var[1];
    MX qend = integr.call(new Dictionary<string, MX> { ["x0"] = q0, ["p"] = par })["xf"];

    MX parc = 0.0;

    var f = ca.Function("f", new[] { var, parMX }, new[] { qend[0] });
    var nlp = new Dictionary<string, MX> { ["x"] = var, ["f"] = -f.call(new[] { var, parc })[0] };
    var opts = new Dictionary<string, object> {
      ["ipopt.tol"] = 1e-12, ["ipopt.hessian_approximation"] = "limited-memory", ["ipopt.max_iter"] = 10,
      ["ipopt.derivative_test"] = "first-order", ["ipopt.print_level"] = 0 };
    Function solver = ca.nlpsol("solver", "ipopt", nlp, opts);
    var solver_out = solver.call(new Dictionary<string, DM> { ["lbx"] = new double[] { -1, -1 }, ["ubx"] = new double[] { 1, 0.2 } });
    DM x = solver_out["x"], lam_x = solver_out["lam_x"];
    Console.WriteLine(x);
    assertAlmostEqual((double)x[0], 1, 7, "X_opt");
    assertAlmostEqual((double)x[1], 0.2, 7, "X_opt");
    assertAlmostEqual((double)ca.fmax(lam_x, 0)[0], 1, 8, "Cost should be linear in y0");
    assertAlmostEqual((double)ca.fmax(lam_x, 0)[1], dfdp(p0, yc0, te), 8, "Cost should be linear in y0");
    assertAlmostEqual(-(double)solver_out["f"][0], cost(p0, y0, yc0, te), 7, "Cost");
    assertAlmostEqual((double)ca.fmax(-lam_x, 0)[0], 0, 8, "Constraint is supposed to be unactive");
    assertAlmostEqual((double)ca.fmax(-lam_x, 0)[1], 0, 8, "Constraint is supposed to be unactive");
  }

  public void test_singleshooting2() {
    Check.Assume(hasNlpsol("ipopt"), "nlpsol ipopt not available");
    message("Single shooting 2");
    double p0 = 0.2, y0 = 0.2, yc0 = 0.1, dy0 = 0.1, te = 0.4;

    Function integr = shootingIntegrator(te);

    MX var = MX.sym("var", 2, 1);
    MX par = MX.sym("par", 1, 1);

    MX q0 = ca.vertcat(var[0], par);
    MX parl = var[1];
    MX qend = integr.call(new Dictionary<string, MX> { ["x0"] = q0, ["p"] = parl })["xf"];

    MX parc = dy0;

    var f = ca.Function("f", new[] { var, par }, new[] { qend[0] });
    var nlp = new Dictionary<string, MX> { ["x"] = var, ["f"] = -f.call(new[] { var, parc })[0], ["g"] = var[0] - var[1] };
    var opts = new Dictionary<string, object> {
      ["ipopt.tol"] = 1e-12, ["ipopt.hessian_approximation"] = "limited-memory", ["ipopt.max_iter"] = 10,
      ["ipopt.derivative_test"] = "first-order" };
    Function solver = ca.nlpsol("solver", "ipopt", nlp, opts);
    var solver_out = solver.call(new Dictionary<string, DM> {
      ["lbx"] = new double[] { -1, -1 }, ["ubx"] = new double[] { 1, 0.2 }, ["lbg"] = new double[] { -1 }, ["ubg"] = new double[] { 0 } });
    DM x = solver_out["x"], lam_x = solver_out["lam_x"];

    assertAlmostEqual((double)x[0], 0.2, 6, "X_opt");
    assertAlmostEqual((double)x[1], 0.2, 6, "X_opt");

    assertAlmostEqual((double)ca.fmax(lam_x, 0)[0], 0, 8, "Constraint is supposed to be unactive");
    double dfdp0 = dfdp(p0, yc0, te);
    assertAlmostEqual((double)ca.fmax(lam_x, 0)[1], 1 + dfdp0, 8);
    assertAlmostEqual((double)solver_out["lam_g"][0], 1, 8);
    assertAlmostEqual(-(double)solver_out["f"][0], cost(p0, y0, yc0, te), 7, "Cost");
    assertAlmostEqual((double)ca.fmax(-lam_x, 0)[0], 0, 8, "Constraint is supposed to be unactive");
    assertAlmostEqual((double)ca.fmax(-lam_x, 0)[1], 0, 8, "Constraint is supposed to be unactive");
  }

  /// <summary>A multiple-shooting test problem: NLP, solver arguments, equality flags</summary>
  sealed class Problem {
    public Dictionary<string, MX> prob;
    public Dictionary<string, DM> args;
    public bool[] equality;
  }

  /// <summary>Accumulates a multiple-shooting NLP like the Python w/w0/lbw/... lists</summary>
  sealed class Nlp {
    readonly List<MX> w_ = new List<MX>(), g_ = new List<MX>();
    readonly List<double> w0 = new List<double>(), lbw = new List<double>(), ubw = new List<double>();
    readonly List<double> lbg = new List<double>(), ubg = new List<double>();
    readonly List<bool> equality = new List<bool>();
    public MX J = 0.0;

    public void w(MX v, double[] lb, double[] ub, double[] x0) {
      w_.Add(v);
      lbw.AddRange(lb);
      ubw.AddRange(ub);
      w0.AddRange(x0);
    }

    public void g(MX v, double[] lb, double[] ub, params bool[] eq) {
      g_.Add(v);
      lbg.AddRange(lb);
      ubg.AddRange(ub);
      equality.AddRange(eq);
    }

    public Problem yield(MX p) => new Problem {
      prob = new Dictionary<string, MX> { ["f"] = J, ["x"] = ca.vertcat(w_), ["g"] = ca.vertcat(g_), ["p"] = p },
      args = new Dictionary<string, DM> {
        ["x0"] = w0.ToArray(), ["lbx"] = lbw.ToArray(), ["ubx"] = ubw.ToArray(), ["lbg"] = lbg.ToArray(),
        ["ubg"] = ubg.ToArray(), ["p"] = 0 },
      equality = equality.ToArray() };
  }

  static double[] d(params double[] x) => x;

  static MX vdp(MX x1, MX x2, MX u, MX p) => ca.vertcat((1 - ca.pow(x2, 2)) * x1 - x2 + u + p, x1);

  static Function rkF(MX x, MX p, MX u, MX xdot) =>
    ca.integrator("F", "rk", new Dictionary<string, MX> { ["x"] = x, ["p"] = p, ["u"] = u, ["ode"] = xdot }, 0, 1,
        new Dictionary<string, object> { ["simplify"] = true, ["number_of_finite_elements"] = 1 });

  static MX xf(Function F, MX Xk, MX Uk, MX p) =>
    F.call(new Dictionary<string, MX> { ["x0"] = Xk, ["u"] = Uk, ["p"] = p })["xf"];

  static List<Problem> fatrop_test_problems() {
    var r = new List<Problem>();
    for (int i = 0; i < 2; i++) {
      int N = 10;
      MX x1 = MX.sym("x1"), x2 = MX.sym("x2");
      MX x = ca.vertcat(x1, x2);
      MX u = MX.sym("u");
      MX p = MX.sym("p");
      Function F = rkF(x, p, u, vdp(x1, x2, u, p));
      var n = new Nlp();

      MX Xk = MX.sym("X0", 2);
      n.w(Xk, d(0, 1), d(0, 1), d(0.1, 0.2));
      MX Xk_next = null;
      for (int k = 0; k < N; k++) {
        MX Uk = MX.sym("U_" + k);
        n.w(Uk, d(-1), d(1), d(0.3));
        MX Xk_end = xf(F, Xk, Uk, p);
        n.J = n.J + ca.sumsqr(Xk) + ca.sumsqr(Uk);
        Xk_next = MX.sym("X_" + (k + 1), 2);
        n.w(Xk_next, d(i == 0 ? -0.25 : -inf, -inf), d(inf, inf), d(0.1, 0.2));
        n.g(Xk_next - Xk_end, d(0, 0), d(0, 0), true, true);
        if (i >= 1) n.g(ca.sin(Xk[0]), d(-0.25), d(inf), false);
        Xk = Xk_next;
      }
      if (i >= 2) {
        Xk = MX.sym("X0", 2);
        n.w(Xk, d(-inf, -inf), d(inf, inf), d(0.1, 0.2));
        n.g(Xk_next - Xk, d(0, 0), d(0, 0), true, true);
        for (int k = 0; k < N; k++) {
          MX Uk = MX.sym("U_" + k);
          n.w(Uk, d(-0.1), d(0.1), d(0.3));
          MX Xk_end = xf(F, Xk, Uk, p);
          n.J = n.J + ca.sumsqr(Xk) * 3 + ca.sumsqr(Uk);
          Xk_next = MX.sym("X_" + (k + 1), 2);
          n.w(Xk_next, d(-inf, -inf), d(inf, inf), d(0.1, 0.2));
          n.g(Xk_next - Xk_end, d(0, 0), d(0, 0), true, true);
          Xk = Xk_next;
        }
      }
      if (i >= 3) {
        x = MX.sym("x", 3);
        u = MX.sym("u", 2);
        DM A = new double[,] { { 1, 0, 0.3 }, { 0, 1, 0.7 }, { 0.2, 0, 1 } };
        DM B = new double[,] { { 1, 0 }, { 0, 1 }, { 0.5, 0.5 } };
        DM D = new double[,] { { 0.2, 0.3 }, { 0.8, 0.7 }, { 0.1, 1 } };
        F = ca.Function("F", new[] { x, u }, new[] { ca.mtimes((MX)A, x) + ca.mtimes((MX)B, u) });
        Xk = MX.sym("X0", 3);
        n.w(Xk, d(-inf, -inf, -inf), d(inf, inf, inf), d(0.1, 0.2, 0.3));
        n.g(ca.mtimes((MX)D, Xk_next) - Xk, d(0, 0, 0), d(0, 0, 0), true, true, true);
        for (int k = 0; k < N; k++) {
          MX Uk = MX.sym("U_" + k, 2);
          n.w(Uk, d(-1, -1), d(1, 1), d(0.3, 0.3));
          MX Xk_end = F.Call(Xk, Uk);
          n.J = n.J + ca.sumsqr(Xk) + ca.sumsqr(Uk);
          Xk_next = MX.sym("X_" + (k + 1), 3);
          n.w(Xk_next, d(-inf, -inf, -inf), d(inf, inf, inf), d(0.1, 0.2, 0.3));
          n.g(Xk_next - Xk_end, d(0, 0, 0), d(0, 0, 0), true, true, true);
          Xk = Xk_next;
        }
      }
      r.Add(n.yield(p));
    }

    {
      // Multi-stage with varying number of inequalities
      int N = 10;
      MX x1 = MX.sym("x1"), x2 = MX.sym("x2");
      MX x = ca.vertcat(x1, x2);
      MX u = MX.sym("u");
      MX p = MX.sym("p");
      Function F = rkF(x, p, u, vdp(x1, x2, u, p));
      var n = new Nlp();

      MX Xk = MX.sym("X0", 2);
      n.w(Xk, d(0, 1), d(0, 1), d(0.1, 0.2));
      for (int k = 0; k < N; k++) {
        MX Uk1 = MX.sym("U1_" + k);
        n.w(Uk1, d(-1), d(1), d(0.3));
        MX Xk_end1 = xf(F, Xk, Uk1, p);
        n.J = n.J + ca.sumsqr(Xk) + ca.sumsqr(Uk1);
        MX Xk_next = MX.sym("X1_" + (k + 1), 2);
        n.w(Xk_next, d(-0.25, -inf), d(inf, inf), d(0.1, 0.2));
        n.g(Xk_next - Xk_end1, d(0, 0), d(0, 0), true, true);
        Xk = Xk_next;
      }

      MX Uk = MX.sym("U1");
      n.w(Uk, d(-1), d(1), d(0.3));
      MX Xk_end = xf(F, Xk, Uk, p);
      n.J = n.J + ca.pow(Xk[0], 2);

      Xk = MX.sym("X1", 2);
      n.w(Xk, d(-inf, -inf), d(inf, inf), d(0.1, 0.2));
      n.g(Xk - Xk_end, d(0, 0), d(0, 0), true, true);

      DM A = new double[,] { { 1, 0.1 }, { 0.2, 1.1 } };
      DM B = new double[,] { { 0.2 }, { 0.7 } };
      F = ca.Function("F", new[] { x, u }, new[] { ca.mtimes((MX)A, x) + ca.mtimes((MX)B, u) });

      for (int k = 0; k < N; k++) {
        Uk = MX.sym("U2_" + k);
        n.w(Uk, d(-inf), d(inf), d(0.3));
        Xk_end = F.Call(Xk, Uk);
        n.J = n.J + ca.sumsqr(Xk) * 3 + ca.sumsqr(Uk);
        MX Xk_next = MX.sym("X2_" + (k + 1), 2);
        n.w(Xk_next, d(-inf, -inf), d(inf, inf), d(0.1, 0.2));
        n.g(Xk_next - Xk_end, d(0, 0), d(0, 0), true, true);
        n.g(Uk * 2, d(-0.1), d(0.1), false);
        Xk = Xk_next;
      }
      r.Add(n.yield(p));
    }

    {
      // Equality constraints
      int N = 10;
      MX x1 = MX.sym("x1"), x2 = MX.sym("x2");
      MX x = ca.vertcat(x1, x2);
      MX u1 = MX.sym("u1"), u2 = MX.sym("u2");
      MX u = ca.vertcat(u1, u2);
      MX p = MX.sym("p");
      MX xdot = ca.vertcat((1 - ca.pow(x2, 2)) * x1 - x2 + u1 + p, x1 + u2);
      Function F = rkF(x, p, u, xdot);
      var n = new Nlp();

      MX Xk = MX.sym("X0", 2);
      n.w(Xk, d(0, 1), d(0, 1), d(0.1, 0.2));
      int i = 0;
      for (int k = 0; k < N; k++) {
        MX Uk = MX.sym("U_" + k, 2);
        n.w(Uk, d(-1, 0.1), d(1, 0.1), d(0.3, 0));
        MX Xk_end = xf(F, Xk, Uk, p);
        n.J = n.J + ca.sumsqr(Xk) + ca.sumsqr(Uk);
        MX Xk_next = MX.sym("X_" + (k + 1), 2);
        n.w(Xk_next, d(i == 0 ? -0.25 : -inf, -inf), d(inf, inf), d(0.1, 0.2));
        n.g(Xk_next - Xk_end, d(0, 0), d(0, 0), true, true);
        Xk = Xk_next;
      }
      r.Add(n.yield(p));
    }

    {
      int N = 10;
      MX x1 = MX.sym("x1"), x2 = MX.sym("x2");
      MX x = ca.vertcat(x1, x2);
      MX u = MX.sym("u");
      MX p = MX.sym("p");
      Function F = rkF(x, p, u, vdp(x1, x2, u, p));
      var n = new Nlp();

      MX Xk = MX.sym("X0", 2);
      n.w(Xk, d(0, 1), d(0, 1), d(0.1, 0.2));
      MX Xk_next = null;
      for (int k = 0; k < N; k++) {
        MX Uk = MX.sym("U1_" + k);
        n.w(Uk, d(-1), d(1), d(0.3));
        MX Xk_end = xf(F, Xk, Uk, p);
        n.J = n.J + ca.sumsqr(Xk) + ca.sumsqr(Uk);
        Xk_next = MX.sym("X1_" + (k + 1), 2);
        n.w(Xk_next, d(-0.25, -inf), d(inf, inf), d(0.1, 0.2));
        n.g(Xk_next - Xk_end, d(0, 0), d(0, 0), true, true);
        Xk = Xk_next;
      }

      n.J = n.J + ca.pow(Xk[0], 2);

      DM A = new double[,] { { 1, 0, 0.3 }, { 0, 1, 0.7 }, { 0.2, 0, 1 } };
      DM B = new double[,] { { 1, 0 }, { 0, 1 }, { 0.5, 0.5 } };
      DM D = new double[,] { { 0.2, 0.3 }, { 0.8, 0.7 }, { 0.1, 1 } };

      Xk = MX.sym("X1", 3);
      n.w(Xk, d(-inf, -inf, -inf), d(inf, inf, inf), d(0.7, 0.8, 0.9));
      n.g(Xk - ca.mtimes((MX)D, Xk_next), d(0, 0, 0), d(0, 0, 0), true, true, true);

      u = MX.sym("u", 2);
      x = MX.sym("x", 3);
      F = ca.Function("F", new[] { x, u }, new[] { ca.mtimes((MX)A, x) + ca.mtimes((MX)B, u) });

      for (int k = 0; k < N; k++) {
        MX Uk = MX.sym("U2_" + k, 2);
        n.w(Uk, d(-0.1, -0.1), d(0.1, 0.1), d(1.3, 1.2));
        MX Xk_end = F.Call(Xk, Uk);
        n.J = n.J + ca.sumsqr(Xk) * 3 + ca.sumsqr(Uk);
        Xk_next = MX.sym("X2_" + (k + 1), 3);
        n.w(Xk_next, d(-inf, -inf, -inf), d(inf, inf, inf), d(0.7, 0.8, 0.9));
        n.g(Xk_next - Xk_end, d(0, 0, 0), d(0, 0, 0), true, true, true);
        Xk = Xk_next;
      }
      r.Add(n.yield(p));
    }
    return r;
  }

  /// <summary>As helpers.py check_serialize: a deserialized copy gives the same outputs</summary>
  static void check_serialize(Function f, IDictionary<string, DM> args) {
    Function f2 = Function.deserialize(f.serialize());
    Dictionary<string, DM> a = f.call(args), b = f2.call(args);
    foreach (string k in a.Keys) checkarray(b[k], a[k], 1e-15, "serialize " + k);
  }

  /// <summary>As helpers.py check_codegen: compile generated C into a shared library and compare outputs</summary>
  void check_codegen(Function f, IDictionary<string, DM> args, string std, IEnumerable<string> extralibs,
                     IEnumerable<string> extra_options, int digits) {
    if (!RunSlow) return;
    string name = "codegen_" + random.Next().ToString("x");
    var cg = ca.CodeGenerator(name, new Dictionary<string, object>());
    cg.add(f);
    cg.generate();
    string libdir = GlobalOptions.getCasadiPath(), incdir = GlobalOptions.getCasadiIncludePath();
    if (libdir.Length == 0 || incdir.Length == 0) Check.Skip("casadi include path not set (build tree)");
    var psi = new ProcessStartInfo("gcc");
    foreach (string a in new[] { "-pedantic", "-std=" + std, "-fPIC", "-shared", "-O3", "-I" + incdir, "-I" + incdir + "/highs",
                                 name + ".c", "-o", name + ".so", "-L" + libdir, "-Wl,-rpath," + libdir })
      psi.ArgumentList.Add(a);
    foreach (string l in extralibs) psi.ArgumentList.Add("-l" + l);
    foreach (string o in extra_options) psi.ArgumentList.Add(o);
    using (var pr = Process.Start(psi)) {
      pr.WaitForExit();
      assertEquals(0, pr.ExitCode, "compile " + name);
    }
    Function f2 = ca.external(f.name(), "./" + name + ".so");
    Dictionary<string, DM> ra = f.call(args), rb = f2.call(args);
    foreach (string k in ra.Keys) checkarray(rb[k], ra[k], Math.Pow(10, -digits), "codegen " + k);
    File.Delete(name + ".c");
    File.Delete(name + ".so");
  }

  public void test_fatrop() {
    Check.Assume(hasNlpsol("fatrop"), "nlpsol fatrop not available");
    Check.Assume(hasNlpsol("ipopt"), "nlpsol ipopt not available");
    var flags = RuntimeInformation.IsOSPlatform(OSPlatform.Windows) ? new string[0] : new[] { "-Wno-strict-prototypes" };
    int local_codegen_check_digits = 15;
    foreach (Problem t in fatrop_test_problems()) {
      var prob = t.prob;
      var args = t.args;
      bool[] equality = t.equality;

      ca.jacobian_sparsity(prob["g"], prob["x"]).spy();

      var solutions = new Dictionary<string, Dictionary<string, DM>>();
      var stats = new Dictionary<string, Dictionary<string, object>>();
      var solvers = new (string, Dictionary<string, object>)[] {
        ("ipopt", new Dictionary<string, object>()),
        ("fatrop", new Dictionary<string, object> {
          ["structure_detection"] = "auto",
          ["fatrop"] = new Dictionary<string, object> { ["tol"] = 1e-8, ["max_iter"] = 100 },
          ["equality"] = equality }) };
      foreach (var (solver, sopts) in solvers) {
        Function f = ca.nlpsol("solver", solver, prob, sopts);
        solutions[solver] = f.call(args);
        stats[solver] = f.stats();
        if (solver != "ipopt") {
          check_codegen(f, args, "c99", new[] { "fatrop", "blasfeo" }, flags, local_codegen_check_digits);
          check_serialize(f, args);
        }
      }
      foreach (string k in solutions["ipopt"].Keys) {
        if (new[] { "x", "f", "g", "lam_g", "lam_x", "lam_p" }.Contains(k)) {
          DM v_ref = solutions["ipopt"][k];
          DM v = solutions["fatrop"][k];
          checkarray(v, v_ref, 1e-5, k);
        }
      }
      long ii = Convert.ToInt64(stats["ipopt"]["iter_count"]);
      long fi = Convert.ToInt64(stats["fatrop"]["iter_count"]);
      assertTrue(Math.Abs(ii - fi) <= 2);
    }
  }

  public void test_fatrop_sanitize() {
    Check.Assume(hasNlpsol("fatrop"), "nlpsol fatrop not available");
    int N = 10;
    MX x1 = MX.sym("x1"), x2 = MX.sym("x2");
    MX x = ca.vertcat(x1, x2);
    MX u = MX.sym("u");
    MX p = MX.sym("p");
    Function F = rkF(x, p, u, vdp(x1, x2, u, p));
    var n = new Nlp();

    MX Xk = MX.sym("X0", 2);
    n.w(Xk, d(0, 1), d(0, 1), d(0.1, 0.2));
    for (int k = 0; k < N; k++) {
      MX Uk = MX.sym("U_" + k);
      n.w(Uk, d(-1), d(1), d(0.3));
      MX Xk_end = xf(F, Xk, Uk, p);
      n.J = n.J + ca.sumsqr(Xk) + ca.sumsqr(Uk);
      MX Xk_next = MX.sym("X_" + (k + 1), 2);
      n.w(Xk_next, d(-0.25, -inf), d(inf, inf), d(0.1, 0.2));
      n.g(Xk_end - Xk_next, d(0, 0), d(0, 0), true, true);
      Xk = Xk_next;
    }
    Problem t = n.yield(p);

    ca.jacobian_sparsity(t.prob["g"], t.prob["x"]).spy();
    Function f = ca.nlpsol("solver", "fatrop", t.prob,
        new Dictionary<string, object> { ["structure_detection"] = "auto", ["verbose"] = true, ["equality"] = t.equality });
    assertInAnyOutput("gap-closing", () => f.call(t.args));
  }

  public void test_detect() {
    Check.Assume(hasNlpsol("fatrop"), "nlpsol fatrop not available");
    foreach (long nx1 in new long[] { 2 })
      foreach (long nx2 in new long[] { 2, 0 })
        foreach (long nu1 in new long[] { 2, 0 })
          foreach (long nu2 in new long[] { 2, 0 })
            foreach (long ng1 in new long[] { 2, 0 })
              foreach (long ng2 in new long[] { 2, 0 })
                foreach (long ng3 in new long[] { 2, 0 }) {
                  Console.WriteLine($"test_detect {nx1} {nx2} {nu1} {nu2} {ng1} {ng2} {ng3}");
                  fatrop_case(N: 2, nx0: 2, nu0: 2, nx1: nx1, nu1: nu1, nx2: nx2, nu2: nu2, ng1: ng1, ng2: ng2, ng3: ng3);
                }
  }

  public void test_detect_adversarial() {
    Check.Assume(hasNlpsol("fatrop"), "nlpsol fatrop not available");
    Sparsity D2 = ca.sparsify((DM)new double[,] { { 1, 0, 0 }, { 1, 1, 1 } }).sparsity();
    fatrop_case(nu2: 3, sp: new Dictionary<string, Sparsity> { ["D2"] = D2 });

    D2 = ca.sparsify((DM)new double[,] { { 1, 0, 0 }, { 1, 1, 1 } }).sparsity();
    Sparsity C2 = ca.sparsify((DM)new double[,] { { 0, 1 }, { 0, 0 } }).sparsity();
    fatrop_case(nu2: 3, sp: new Dictionary<string, Sparsity> { ["D2"] = D2, ["C2"] = C2 });
    assertInAnyOutput("gap-closing constraints must be like",
      () => fatrop_case(nu2: 3, sp: new Dictionary<string, Sparsity> { ["D2"] = D2, ["C2"] = C2 }, eq: new[] { "ng3" }));

    fatrop_case(nu0: 0, nx0: 1);

    fatrop_case(nx2: 0);

    fatrop_case(nu2: 3, sp: new Dictionary<string, Sparsity> { ["A1"] = ca.Sparsity(2, 2), ["B1"] = ca.Sparsity(2, 2) });
    fatrop_case(nu2: 3, sp: new Dictionary<string, Sparsity> { ["A1"] = ca.Sparsity(2, 2) });

    fatrop_case(nx2: 1, ng1: 0, nu1: 0);
  }

  public void test_bug() {
    Check.Assume(hasNlpsol("fatrop"), "nlpsol fatrop not available");
    MX x = MX.sym("x");

    foreach (string structure_detection in new[] { "none", "auto" }) {
      var opts = new Dictionary<string, object> {
        ["expand"] = true, ["structure_detection"] = structure_detection, ["equality"] = new[] { true } };

      Function solver = ca.nlpsol("solver", "fatrop", new Dictionary<string, MX> { ["x"] = x, ["g"] = x - 1 }, opts);
      assertAlmostEqual((double)solver.call(new Dictionary<string, DM> { ["lbg"] = 0, ["ubg"] = 0 })["x"], 1, 5);

      solver = ca.nlpsol("solver", "fatrop", new Dictionary<string, MX> { ["x"] = x, ["g"] = x }, opts);
      assertAlmostEqual((double)solver.call(new Dictionary<string, DM> { ["lbg"] = 1, ["ubg"] = 1 })["x"], 1, 5);

      solver = ca.nlpsol("solver", "fatrop", new Dictionary<string, MX> { ["x"] = x, ["g"] = x - 2 }, opts);
      assertAlmostEqual((double)solver.call(new Dictionary<string, DM> { ["lbg"] = 3, ["ubg"] = 3 })["x"], 5, 5);
    }
  }
}

}

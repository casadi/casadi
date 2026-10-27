package org.casadi.test;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.casadi.CodeGenerator;
import org.casadi.DM;
import org.casadi.Function;
import org.casadi.GenericMatrix;
import org.casadi.GlobalOptions;
import org.casadi.MX;
import org.casadi.SX;
import org.casadi.Slice;
import org.casadi.Sparsity;
import org.casadi.ca;
import org.junit.jupiter.api.Assumptions;
import org.junit.jupiter.api.Test;

/** Port of test/python/ocp.py */
public class OcpTest extends CasadiTestCase {
  private static final double inf = Double.POSITIVE_INFINITY;

  private static void assertAlmostEqual(double a, double b, int places, String msg) {
    assertTrue(Math.abs(a - b) < 0.5 * Math.pow(10, -places), msg + ": " + a + " vs " + b);
  }

  private static void assertAlmostEqual(double a, double b, int places) {
    assertAlmostEqual(a, b, places, "");
  }

  private static DM out(Map<String, GenericMatrix> r, String k) {
    return (DM) r.get(k);
  }

  private static boolean[] rep(boolean v, int n) {
    boolean[] r = new boolean[n];
    Arrays.fill(r, v);
    return r;
  }

  private static boolean[] cat(boolean[]... a) {
    int n = 0;
    for (boolean[] b : a) n += b.length;
    boolean[] r = new boolean[n];
    int k = 0;
    for (boolean[] b : a) for (boolean v : b) r[k++] = v;
    return r;
  }

  private static DM proj(DM a, Map<String, Sparsity> sp, String k) {
    return sp.containsKey(k) ? ca.project(a, sp.get(k)) : a;
  }

  private void fatrop_case(Map<String, Object> kw) {
    long N = (long) kw.getOrDefault("N", 2L);
    long nx0 = (long) kw.getOrDefault("nx0", 2L), nu0 = (long) kw.getOrDefault("nu0", 2L);
    long nx1 = (long) kw.getOrDefault("nx1", 2L), nu1 = (long) kw.getOrDefault("nu1", 2L);
    long nx2 = (long) kw.getOrDefault("nx2", 2L), nu2 = (long) kw.getOrDefault("nu2", 2L);
    long ng1 = (long) kw.getOrDefault("ng1", 2L), ng2 = (long) kw.getOrDefault("ng2", 2L);
    long ng3 = (long) kw.getOrDefault("ng3", 2L);
    @SuppressWarnings("unchecked")
    Map<String, Sparsity> sp = (Map<String, Sparsity>) kw.getOrDefault("sp", Map.of());
    @SuppressWarnings("unchecked")
    Set<String> eq = (Set<String>) kw.getOrDefault("eq", Set.of());
    System.out.println("fatrop_case " + kw);
    long[] nx = {nx0, nx1, nx2}, nu = {nu0, nu1, nu2}, ng = {ng1, ng2, ng3};

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

    DM A = ca.blockcat(List.of(
        List.of(A0, B0, I0, new DM(nx1, nu1 + nx2 + nu2)),
        List.of(C0, D0, new DM(ng1, nx1 + nu1 + nx2 + nu2)),
        List.of(new DM(nx2, nx0 + nu0), A1, B1, I1, new DM(nx2, nu2)),
        List.of(new DM(ng2, nx0 + nu0), C1, D1, new DM(ng2, nx2 + nu2)),
        List.of(new DM(ng3, nx0 + nu0 + nx1 + nu1), C2, D2)));

    boolean[] equality = cat(rep(true, (int) nx1), rep(eq.contains("ng1"), (int) ng1), rep(true, (int) nx2),
        rep(eq.contains("ng2"), (int) ng2), rep(eq.contains("ng3"), (int) ng3));

    A.sparsity().spy();
    System.out.println(A);

    MX x = ca.vertcat(List.of(MX.sym("x0", nx0), MX.sym("u0", nu0), MX.sym("x1", nx1), MX.sym("u1", nu1),
        MX.sym("x2", nx2), MX.sym("u2", nu2)));
    Map<String, MX> nlp = new HashMap<>();
    nlp.put("x", x);
    nlp.put("g", ca.plus(DM.zeros(A.size1(), 1), ca.mtimes(A, x)));
    nlp.put("f", ca.sumsqr(x.minus(DM.rand(x.numel(), 1))));

    double a = 10;
    DM lbg = ca.vertcat(List.of(DM.zeros(nx1, 1), DM.ones(ng1, 1).times(-a), DM.zeros(nx2, 1),
        DM.ones(ng2, 1).times(-a), DM.ones(ng3, 1).times(-a)));
    DM ubg = ca.vertcat(List.of(DM.zeros(nx1, 1), DM.ones(ng1, 1).times(a), DM.zeros(nx2, 1),
        DM.ones(ng2, 1).times(a), DM.ones(ng3, 1).times(a)));
    System.out.println(lbg);

    Map<String, Object> options = new LinkedHashMap<>();
    options.put("structure_detection", "manual");
    options.put("N", N);
    options.put("nx", nx);
    options.put("nu", nu);
    options.put("ng", ng);
    options.put("equality", equality);
    options.put("fatrop", Map.of("tol", 1e-7));
    Function solver = ca.nlpsol("solver", "fatrop", nlp, options);
    Map<String, GenericMatrix> sol = solver.invoke(Map.of("lbg", lbg, "ubg", ubg));

    solver = ca.nlpsol("solver", "fatrop", nlp, Map.of("structure_detection", "none", "error_on_fail", true,
        "equality", equality, "fatrop", Map.of("tol", 1e-7)));
    Map<String, GenericMatrix> ref = solver.invoke(Map.of("lbg", lbg, "ubg", ubg));

    for (String k : sol.keySet()) checkarray(out(sol, k), out(ref, k), 1e-6, k + options);

    options = Map.of("structure_detection", "auto", "debug", true, "equality", equality, "fatrop",
        Map.of("tol", 1e-7));
    System.out.println(options);
    solver = ca.nlpsol("solver", "fatrop", nlp, options);
    sol = solver.invoke(Map.of("lbg", lbg, "ubg", ubg));

    Map<String, Object> stats = solver.stats();
    System.out.println(stats);
    if (nx2 > 0) assertTrue(((Number) stats.get("N")).longValue() > 1);

    for (String k : sol.keySet()) checkarray(out(sol, k), out(ref, k), 1e-6, k + options);
  }

  @Test
  public void testdiscrete() {
    Assumptions.assumeTrue(hasNlpsol("ipopt"), "nlpsol ipopt not available");
    message("Linear-quadratic problem, discrete, using IPOPT");
    double a = 1.0, b = 1.0, q = 1.0, s = 1.0, r = 1.0;
    double x0 = 100;

    int N = 100;

    SX X = SX.sym("X", N + 1);
    SX U = SX.sym("U", N);

    SX V = ca.vertcat(List.of(X, U));

    SX cost = new SX(0.0);
    for (int i = 0; i < N; i++)
      cost = cost.plus(X.get(i).power(2).times(s)).plus(U.get(i).power(2).times(r));
    cost = cost.plus(X.get(N).power(2).times(q));

    Map<String, SX> nlp = Map.of("x", V, "f", cost, "g", ca.vertcat(List.of(X.get(0).minus(x0),
        X.get(new Slice(1, N + 1), new Slice(0, 1)).minus(X.get(new Slice(0, N), new Slice(0, 1)).times(a).plus(U.times(b))))));
    Map<String, Object> opts = new HashMap<>();
    opts.put("ipopt.tol", 1e-5);
    opts.put("ipopt.hessian_approximation", "limited-memory");
    opts.put("ipopt.max_iter", 100);
    opts.put("ipopt.print_level", 0);
    Function solver = ca.nlpsol("solver", "ipopt", nlp, opts);
    double[] lbx = new double[(int) V.nnz()], ubx = new double[(int) V.nnz()];
    Arrays.fill(lbx, -1000);
    Arrays.fill(ubx, 1000);
    Map<String, GenericMatrix> solver_out = solver.invoke(
        Map.of("lbx", lbx, "ubx", ubx, "lbg", new double[N + 1], "ubg", new double[N + 1]));
    double ocp_sol = out(solver_out, "f").get(0).to_double();
    double K = q + 0.0;
    for (int i = 0; i < N; i++) K = s + r * a * a * K / (r + b * b * K);
    double exact_sol = K * x0 * x0;
    assertAlmostEqual(ocp_sol, exact_sol, 10, "Linear-quadratic problem solution using IPOPT");
  }

  private static double dfdp(double p0, double yc0, double te) {
    return (Math.sqrt(p0) * (te * yc0 * yc0 - yc0 + p0 * te) * Math.tan(Math.atan(yc0 / Math.sqrt(p0)) + Math.sqrt(p0) * te)
        + yc0 * yc0) / (2 * p0 * yc0 * yc0 + 2 * p0 * p0);
  }

  private static double cost(double p0, double y0, double yc0, double te) {
    return (2 * y0 - Math.log(yc0 * yc0 / p0 + 1)) / 2 - Math.log(Math.cos(Math.atan(yc0 / Math.sqrt(p0)) + Math.sqrt(p0) * te));
  }

  private static Function shootingIntegrator(double te) {
    SX t = SX.sym("t");
    SX q = SX.sym("y", 2, 1);
    SX p = SX.sym("p", 1, 1);
    Map<String, SX> dae = Map.of("x", q, "p", p, "t", t, "ode", ca.vertcat(List.of(q.get(1), p.get(0).plus(q.get(1).power(2)))));
    Map<String, Object> opts = new HashMap<>();
    opts.put("reltol", 1e-15);
    opts.put("abstol", 1e-15);
    opts.put("verbose", false);
    opts.put("steps_per_checkpoint", 10000);
    return ca.integrator("integrator", "cvodes", dae, 0, te, opts);
  }

  @Test
  public void test_singleshooting() {
    Assumptions.assumeTrue(hasNlpsol("ipopt"), "nlpsol ipopt not available");
    message("Single shooting");
    double p0 = 0.2, y0 = 1, yc0 = 0, te = 0.4;

    Function integrator = shootingIntegrator(te);

    MX var = MX.sym("var", 2, 1);
    MX par = MX.sym("par", 1, 1);
    MX parMX = par;

    MX q0 = ca.vertcat(List.of(var.get(0), par));
    par = var.get(1);
    MX qend = integrator.call(Map.of("x0", q0, "p", par)).get("xf");

    MX parc = new MX(0.0);

    Function f = new Function("f", List.of(var, parMX), List.of(qend.get(0)));
    Map<String, MX> nlp = Map.of("x", var, "f", f.call(List.of(var, parc)).get(0).unaryMinus());
    Map<String, Object> opts = new HashMap<>();
    opts.put("ipopt.tol", 1e-12);
    opts.put("ipopt.hessian_approximation", "limited-memory");
    opts.put("ipopt.max_iter", 10);
    opts.put("ipopt.derivative_test", "first-order");
    opts.put("ipopt.print_level", 0);
    Function solver = ca.nlpsol("solver", "ipopt", nlp, opts);
    Map<String, GenericMatrix> solver_out = solver.invoke(Map.of("lbx", new double[] {-1, -1}, "ubx", new double[] {1, 0.2}));
    DM x = out(solver_out, "x"), lam_x = out(solver_out, "lam_x");
    System.out.println(x);
    assertAlmostEqual(x.get(0).to_double(), 1, 7, "X_opt");
    assertAlmostEqual(x.get(1).to_double(), 0.2, 7, "X_opt");
    assertAlmostEqual(ca.fmax(lam_x, 0).get(0).to_double(), 1, 8, "Cost should be linear in y0");
    assertAlmostEqual(ca.fmax(lam_x, 0).get(1).to_double(), dfdp(p0, yc0, te), 8, "Cost should be linear in y0");
    assertAlmostEqual(-out(solver_out, "f").get(0).to_double(), cost(p0, y0, yc0, te), 7, "Cost");
    assertAlmostEqual(ca.fmax(lam_x.unaryMinus(), 0).get(0).to_double(), 0, 8, "Constraint is supposed to be unactive");
    assertAlmostEqual(ca.fmax(lam_x.unaryMinus(), 0).get(1).to_double(), 0, 8, "Constraint is supposed to be unactive");
  }

  @Test
  public void test_singleshooting2() {
    Assumptions.assumeTrue(hasNlpsol("ipopt"), "nlpsol ipopt not available");
    message("Single shooting 2");
    double p0 = 0.2, y0 = 0.2, yc0 = 0.1, dy0 = 0.1, te = 0.4;

    Function integrator = shootingIntegrator(te);

    MX var = MX.sym("var", 2, 1);
    MX par = MX.sym("par", 1, 1);

    MX q0 = ca.vertcat(List.of(var.get(0), par));
    MX parl = var.get(1);
    MX qend = integrator.call(Map.of("x0", q0, "p", parl)).get("xf");

    MX parc = new MX(dy0);

    Function f = new Function("f", List.of(var, par), List.of(qend.get(0)));
    Map<String, MX> nlp = Map.of("x", var, "f", f.call(List.of(var, parc)).get(0).unaryMinus(), "g",
        var.get(0).minus(var.get(1)));
    Map<String, Object> opts = new HashMap<>();
    opts.put("ipopt.tol", 1e-12);
    opts.put("ipopt.hessian_approximation", "limited-memory");
    opts.put("ipopt.max_iter", 10);
    opts.put("ipopt.derivative_test", "first-order");
    Function solver = ca.nlpsol("solver", "ipopt", nlp, opts);
    Map<String, GenericMatrix> solver_out = solver.invoke(Map.of("lbx", new double[] {-1, -1}, "ubx",
        new double[] {1, 0.2}, "lbg", new double[] {-1}, "ubg", new double[] {0}));
    DM x = out(solver_out, "x"), lam_x = out(solver_out, "lam_x");

    assertAlmostEqual(x.get(0).to_double(), 0.2, 6, "X_opt");
    assertAlmostEqual(x.get(1).to_double(), 0.2, 6, "X_opt");

    assertAlmostEqual(ca.fmax(lam_x, 0).get(0).to_double(), 0, 8, "Constraint is supposed to be unactive");
    double dfdp0 = dfdp(p0, yc0, te);
    assertAlmostEqual(ca.fmax(lam_x, 0).get(1).to_double(), 1 + dfdp0, 8);
    assertAlmostEqual(out(solver_out, "lam_g").get(0).to_double(), 1, 8);
    assertAlmostEqual(-out(solver_out, "f").get(0).to_double(), cost(p0, y0, yc0, te), 7, "Cost");
    assertAlmostEqual(ca.fmax(lam_x.unaryMinus(), 0).get(0).to_double(), 0, 8, "Constraint is supposed to be unactive");
    assertAlmostEqual(ca.fmax(lam_x.unaryMinus(), 0).get(1).to_double(), 0, 8, "Constraint is supposed to be unactive");
  }

  /** Accumulates a multiple-shooting NLP like the Python w/w0/lbw/... lists */
  private static final class Nlp {
    final List<MX> w = new ArrayList<>(), g = new ArrayList<>();
    final List<Double> w0 = new ArrayList<>(), lbw = new ArrayList<>(), ubw = new ArrayList<>();
    final List<Double> lbg = new ArrayList<>(), ubg = new ArrayList<>();
    final List<Boolean> equality = new ArrayList<>();
    MX J = new MX(0.0);

    void w(MX v, double[] lb, double[] ub, double[] x0) {
      w.add(v);
      for (double e : lb) lbw.add(e);
      for (double e : ub) ubw.add(e);
      for (double e : x0) w0.add(e);
    }

    void g(MX v, double[] lb, double[] ub, boolean... eq) {
      g.add(v);
      for (double e : lb) lbg.add(e);
      for (double e : ub) ubg.add(e);
      for (boolean e : eq) equality.add(e);
    }

    private static double[] arr(List<Double> l) {
      return l.stream().mapToDouble(Double::doubleValue).toArray();
    }

    Object[] yield(MX p) {
      Map<String, MX> prob = Map.of("f", J, "x", ca.vertcat(w), "g", ca.vertcat(g), "p", p);
      Map<String, Object> args = Map.of("x0", arr(w0), "lbx", arr(lbw), "ubx", arr(ubw), "lbg", arr(lbg), "ubg", arr(ubg), "p", 0);
      boolean[] eq = new boolean[equality.size()];
      for (int i = 0; i < eq.length; i++) eq[i] = equality.get(i);
      return new Object[] {prob, args, eq};
    }
  }

  private static double[] d(double... x) {
    return x;
  }

  private static MX vdp(MX x1, MX x2, MX u, MX p) {
    return ca.vertcat(List.of(ca.times(ca.minus(1, x2.power(2)), x1).minus(x2).plus(u).plus(p), x1));
  }

  private static Function rkF(MX x, MX p, MX u, MX xdot) {
    return ca.integrator("F", "rk", Map.of("x", x, "p", p, "u", u, "ode", xdot), 0, 1,
        Map.of("simplify", true, "number_of_finite_elements", 1));
  }

  private static MX xf(Function F, MX Xk, MX Uk, MX p) {
    return F.call(Map.of("x0", Xk, "u", Uk, "p", p)).get("xf");
  }

  private static MX call1(Function F, MX a, MX b) {
    return F.call(List.of(a, b)).get(0);
  }

  private static List<Object[]> fatrop_test_problems() {
    List<Object[]> r = new ArrayList<>();
    for (int i = 0; i < 2; i++) {
      int N = 10;
      MX x1 = MX.sym("x1"), x2 = MX.sym("x2");
      MX x = ca.vertcat(List.of(x1, x2));
      MX u = MX.sym("u");
      MX p = MX.sym("p");
      Function F = rkF(x, p, u, vdp(x1, x2, u, p));
      Nlp n = new Nlp();

      MX Xk = MX.sym("X0", 2);
      n.w(Xk, d(0, 1), d(0, 1), d(0.1, 0.2));
      MX Xk_next = null;
      for (int k = 0; k < N; k++) {
        MX Uk = MX.sym("U_" + k);
        n.w(Uk, d(-1), d(1), d(0.3));
        MX Xk_end = xf(F, Xk, Uk, p);
        n.J = n.J.plus(ca.sumsqr(Xk)).plus(ca.sumsqr(Uk));
        Xk_next = MX.sym("X_" + (k + 1), 2);
        n.w(Xk_next, d(i == 0 ? -0.25 : -inf, -inf), d(inf, inf), d(0.1, 0.2));
        n.g(Xk_next.minus(Xk_end), d(0, 0), d(0, 0), true, true);
        if (i >= 1) n.g(ca.sin(Xk.get(0)), d(-0.25), d(inf), false);
        Xk = Xk_next;
      }
      if (i >= 2) {
        Xk = MX.sym("X0", 2);
        n.w(Xk, d(-inf, -inf), d(inf, inf), d(0.1, 0.2));
        n.g(Xk_next.minus(Xk), d(0, 0), d(0, 0), true, true);
        for (int k = 0; k < N; k++) {
          MX Uk = MX.sym("U_" + k);
          n.w(Uk, d(-0.1), d(0.1), d(0.3));
          MX Xk_end = xf(F, Xk, Uk, p);
          n.J = n.J.plus(ca.sumsqr(Xk).times(3)).plus(ca.sumsqr(Uk));
          Xk_next = MX.sym("X_" + (k + 1), 2);
          n.w(Xk_next, d(-inf, -inf), d(inf, inf), d(0.1, 0.2));
          n.g(Xk_next.minus(Xk_end), d(0, 0), d(0, 0), true, true);
          Xk = Xk_next;
        }
      }
      if (i >= 3) {
        x = MX.sym("x", 3);
        u = MX.sym("u", 2);
        DM A = new DM(new double[][] {{1, 0, 0.3}, {0, 1, 0.7}, {0.2, 0, 1}});
        DM B = new DM(new double[][] {{1, 0}, {0, 1}, {0.5, 0.5}});
        DM D = new DM(new double[][] {{0.2, 0.3}, {0.8, 0.7}, {0.1, 1}});
        F = new Function("F", List.of(x, u), List.of(ca.mtimes(A, x).plus(ca.mtimes(B, u))));
        Xk = MX.sym("X0", 3);
        n.w(Xk, d(-inf, -inf, -inf), d(inf, inf, inf), d(0.1, 0.2, 0.3));
        n.g(ca.mtimes(D, Xk_next).minus(Xk), d(0, 0, 0), d(0, 0, 0), true, true, true);
        for (int k = 0; k < N; k++) {
          MX Uk = MX.sym("U_" + k, 2);
          n.w(Uk, d(-1, -1), d(1, 1), d(0.3, 0.3));
          MX Xk_end = call1(F, Xk, Uk);
          n.J = n.J.plus(ca.sumsqr(Xk)).plus(ca.sumsqr(Uk));
          Xk_next = MX.sym("X_" + (k + 1), 3);
          n.w(Xk_next, d(-inf, -inf, -inf), d(inf, inf, inf), d(0.1, 0.2, 0.3));
          n.g(Xk_next.minus(Xk_end), d(0, 0, 0), d(0, 0, 0), true, true, true);
          Xk = Xk_next;
        }
      }
      r.add(n.yield(p));
    }

    {
      // Multi-stage with varying number of inequalities
      int N = 10;
      MX x1 = MX.sym("x1"), x2 = MX.sym("x2");
      MX x = ca.vertcat(List.of(x1, x2));
      MX u = MX.sym("u");
      MX p = MX.sym("p");
      Function F = rkF(x, p, u, vdp(x1, x2, u, p));
      Nlp n = new Nlp();

      MX Xk = MX.sym("X0", 2);
      n.w(Xk, d(0, 1), d(0, 1), d(0.1, 0.2));
      for (int k = 0; k < N; k++) {
        MX Uk = MX.sym("U1_" + k);
        n.w(Uk, d(-1), d(1), d(0.3));
        MX Xk_end = xf(F, Xk, Uk, p);
        n.J = n.J.plus(ca.sumsqr(Xk)).plus(ca.sumsqr(Uk));
        MX Xk_next = MX.sym("X1_" + (k + 1), 2);
        n.w(Xk_next, d(-0.25, -inf), d(inf, inf), d(0.1, 0.2));
        n.g(Xk_next.minus(Xk_end), d(0, 0), d(0, 0), true, true);
        Xk = Xk_next;
      }

      MX Uk = MX.sym("U1");
      n.w(Uk, d(-1), d(1), d(0.3));
      MX Xk_end = xf(F, Xk, Uk, p);
      n.J = n.J.plus(Xk.get(0).power(2));

      Xk = MX.sym("X1", 2);
      n.w(Xk, d(-inf, -inf), d(inf, inf), d(0.1, 0.2));
      n.g(Xk.minus(Xk_end), d(0, 0), d(0, 0), true, true);

      DM A = new DM(new double[][] {{1, 0.1}, {0.2, 1.1}});
      DM B = new DM(new double[][] {{0.2}, {0.7}});
      F = new Function("F", List.of(x, u), List.of(ca.mtimes(A, x).plus(ca.mtimes(B, u))));

      for (int k = 0; k < N; k++) {
        Uk = MX.sym("U2_" + k);
        n.w(Uk, d(-inf), d(inf), d(0.3));
        Xk_end = call1(F, Xk, Uk);
        n.J = n.J.plus(ca.sumsqr(Xk).times(3)).plus(ca.sumsqr(Uk));
        MX Xk_next = MX.sym("X2_" + (k + 1), 2);
        n.w(Xk_next, d(-inf, -inf), d(inf, inf), d(0.1, 0.2));
        n.g(Xk_next.minus(Xk_end), d(0, 0), d(0, 0), true, true);
        n.g(Uk.times(2), d(-0.1), d(0.1), false);
        Xk = Xk_next;
      }
      r.add(n.yield(p));
    }

    {
      // Equality constraints
      int N = 10;
      MX x1 = MX.sym("x1"), x2 = MX.sym("x2");
      MX x = ca.vertcat(List.of(x1, x2));
      MX u1 = MX.sym("u1"), u2 = MX.sym("u2");
      MX u = ca.vertcat(List.of(u1, u2));
      MX p = MX.sym("p");
      MX xdot = ca.vertcat(List.of(ca.times(ca.minus(1, x2.power(2)), x1).minus(x2).plus(u1).plus(p), x1.plus(u2)));
      Function F = rkF(x, p, u, xdot);
      Nlp n = new Nlp();

      MX Xk = MX.sym("X0", 2);
      n.w(Xk, d(0, 1), d(0, 1), d(0.1, 0.2));
      int i = 0;
      for (int k = 0; k < N; k++) {
        MX Uk = MX.sym("U_" + k, 2);
        n.w(Uk, d(-1, 0.1), d(1, 0.1), d(0.3, 0));
        MX Xk_end = xf(F, Xk, Uk, p);
        n.J = n.J.plus(ca.sumsqr(Xk)).plus(ca.sumsqr(Uk));
        MX Xk_next = MX.sym("X_" + (k + 1), 2);
        n.w(Xk_next, d(i == 0 ? -0.25 : -inf, -inf), d(inf, inf), d(0.1, 0.2));
        n.g(Xk_next.minus(Xk_end), d(0, 0), d(0, 0), true, true);
        Xk = Xk_next;
      }
      r.add(n.yield(p));
    }

    int N = 10;
    MX x1 = MX.sym("x1"), x2 = MX.sym("x2");
    MX x = ca.vertcat(List.of(x1, x2));
    MX u = MX.sym("u");
    MX p = MX.sym("p");
    Function F = rkF(x, p, u, vdp(x1, x2, u, p));
    Nlp n = new Nlp();

    MX Xk = MX.sym("X0", 2);
    n.w(Xk, d(0, 1), d(0, 1), d(0.1, 0.2));
    MX Xk_next = null;
    for (int k = 0; k < N; k++) {
      MX Uk = MX.sym("U1_" + k);
      n.w(Uk, d(-1), d(1), d(0.3));
      MX Xk_end = xf(F, Xk, Uk, p);
      n.J = n.J.plus(ca.sumsqr(Xk)).plus(ca.sumsqr(Uk));
      Xk_next = MX.sym("X1_" + (k + 1), 2);
      n.w(Xk_next, d(-0.25, -inf), d(inf, inf), d(0.1, 0.2));
      n.g(Xk_next.minus(Xk_end), d(0, 0), d(0, 0), true, true);
      Xk = Xk_next;
    }

    n.J = n.J.plus(Xk.get(0).power(2));

    DM A = new DM(new double[][] {{1, 0, 0.3}, {0, 1, 0.7}, {0.2, 0, 1}});
    DM B = new DM(new double[][] {{1, 0}, {0, 1}, {0.5, 0.5}});
    DM D = new DM(new double[][] {{0.2, 0.3}, {0.8, 0.7}, {0.1, 1}});

    Xk = MX.sym("X1", 3);
    n.w(Xk, d(-inf, -inf, -inf), d(inf, inf, inf), d(0.7, 0.8, 0.9));
    n.g(Xk.minus(ca.mtimes(D, Xk_next)), d(0, 0, 0), d(0, 0, 0), true, true, true);

    u = MX.sym("u", 2);
    x = MX.sym("x", 3);
    F = new Function("F", List.of(x, u), List.of(ca.mtimes(A, x).plus(ca.mtimes(B, u))));

    for (int k = 0; k < N; k++) {
      MX Uk = MX.sym("U2_" + k, 2);
      n.w(Uk, d(-0.1, -0.1), d(0.1, 0.1), d(1.3, 1.2));
      MX Xk_end = call1(F, Xk, Uk);
      n.J = n.J.plus(ca.sumsqr(Xk).times(3)).plus(ca.sumsqr(Uk));
      Xk_next = MX.sym("X2_" + (k + 1), 3);
      n.w(Xk_next, d(-inf, -inf, -inf), d(inf, inf, inf), d(0.7, 0.8, 0.9));
      n.g(Xk_next.minus(Xk_end), d(0, 0, 0), d(0, 0, 0), true, true, true);
      Xk = Xk_next;
    }
    r.add(n.yield(p));
    return r;
  }

  /** As helpers.py check_serialize: a deserialized copy gives the same outputs */
  private void check_serialize(Function f, Map<String, Object> args) {
    Function f2 = Function.deserialize(f.serialize());
    Map<String, GenericMatrix> a = f.invoke(args), b = f2.invoke(args);
    for (String k : a.keySet()) checkarray((DM) b.get(k), (DM) a.get(k), 1e-15, "serialize " + k);
  }

  /** As helpers.py check_codegen: compile generated C into a shared library and compare outputs */
  private void check_codegen(Function f, Map<String, Object> args, String std, List<String> extralibs,
      List<String> extra_options, int digits) {
    if (!RUN_SLOW) return;
    try {
      String name = "codegen_" + Long.toHexString(random.nextLong() & Long.MAX_VALUE);
      CodeGenerator cg = new CodeGenerator(name, Map.of());
      cg.add(f);
      cg.generate();
      String libdir = GlobalOptions.getCasadiPath(), incdir = GlobalOptions.getCasadiIncludePath();
      List<String> cmd = new ArrayList<>(List.of("gcc", "-pedantic", "-std=" + std, "-fPIC", "-shared", "-O3",
          "-I" + incdir, "-I" + incdir + "/highs", name + ".c", "-o", name + ".so", "-L" + libdir, "-Wl,-rpath," + libdir));
      for (String l : extralibs) cmd.add("-l" + l);
      cmd.addAll(extra_options);
      Process pr = new ProcessBuilder(cmd).inheritIO().start();
      assertEquals(0, pr.waitFor(), "compile " + name);
      Function f2 = ca.external(f.name(), "./" + name + ".so");
      Map<String, GenericMatrix> a = f.invoke(args), b = f2.invoke(args);
      for (String k : a.keySet()) checkarray((DM) b.get(k), (DM) a.get(k), Math.pow(10, -digits), "codegen " + k);
      new File(name + ".c").delete();
      new File(name + ".so").delete();
    } catch (java.io.IOException | InterruptedException e) {
      throw new RuntimeException(e);
    }
  }

  @Test
  @SuppressWarnings("unchecked")
  public void test_fatrop() {
    Assumptions.assumeTrue(hasNlpsol("fatrop"), "nlpsol fatrop not available");
    Assumptions.assumeTrue(hasNlpsol("ipopt"), "nlpsol ipopt not available");
    List<String> flags = System.getProperty("os.name").startsWith("Windows") ? List.of() : List.of("-Wno-strict-prototypes");
    int local_codegen_check_digits = 15;
    for (Object[] t : fatrop_test_problems()) {
      Map<String, MX> prob = (Map<String, MX>) t[0];
      Map<String, Object> args = (Map<String, Object>) t[1];
      boolean[] equality = (boolean[]) t[2];

      ca.jacobian_sparsity(prob.get("g"), prob.get("x")).spy();

      Map<String, Map<String, GenericMatrix>> solutions = new HashMap<>();
      Map<String, Map<String, Object>> stats = new HashMap<>();
      List<Object[]> solvers = List.of(new Object[] {"ipopt", Map.of()}, new Object[] {"fatrop",
          Map.of("structure_detection", "auto", "fatrop", Map.of("tol", 1e-8, "max_iter", 100), "equality", equality)});
      for (Object[] so : solvers) {
        String solver = (String) so[0];
        Function f = ca.nlpsol("solver", solver, prob, (Map<String, ?>) so[1]);
        solutions.put(solver, f.invoke(args));
        stats.put(solver, f.stats());
        if (!solver.equals("ipopt")) {
          check_codegen(f, args, "c99", List.of("fatrop", "blasfeo"), flags, local_codegen_check_digits);
          check_serialize(f, args);
        }
      }
      for (String k : solutions.get("ipopt").keySet()) {
        if (List.of("x", "f", "g", "lam_g", "lam_x", "lam_p").contains(k)) {
          DM v_ref = (DM) solutions.get("ipopt").get(k);
          DM v = (DM) solutions.get("fatrop").get(k);
          checkarray(v, v_ref, 1e-5, k);
        }
      }
      long ii = ((Number) stats.get("ipopt").get("iter_count")).longValue();
      long fi = ((Number) stats.get("fatrop").get("iter_count")).longValue();
      assertTrue(Math.abs(ii - fi) <= 2);
    }
  }

  @Test
  public void test_fatrop_sanitize() {
    Assumptions.assumeTrue(hasNlpsol("fatrop"), "nlpsol fatrop not available");
    int N = 10;
    MX x1 = MX.sym("x1"), x2 = MX.sym("x2");
    MX x = ca.vertcat(List.of(x1, x2));
    MX u = MX.sym("u");
    MX p = MX.sym("p");
    Function F = rkF(x, p, u, vdp(x1, x2, u, p));
    Nlp n = new Nlp();

    MX Xk = MX.sym("X0", 2);
    n.w(Xk, d(0, 1), d(0, 1), d(0.1, 0.2));
    for (int k = 0; k < N; k++) {
      MX Uk = MX.sym("U_" + k);
      n.w(Uk, d(-1), d(1), d(0.3));
      MX Xk_end = xf(F, Xk, Uk, p);
      n.J = n.J.plus(ca.sumsqr(Xk)).plus(ca.sumsqr(Uk));
      MX Xk_next = MX.sym("X_" + (k + 1), 2);
      n.w(Xk_next, d(-0.25, -inf), d(inf, inf), d(0.1, 0.2));
      n.g(Xk_end.minus(Xk_next), d(0, 0), d(0, 0), true, true);
      Xk = Xk_next;
    }
    Object[] t = n.yield(p);
    @SuppressWarnings("unchecked")
    Map<String, MX> prob = (Map<String, MX>) t[0];
    @SuppressWarnings("unchecked")
    Map<String, Object> args = (Map<String, Object>) t[1];
    boolean[] equality = (boolean[]) t[2];

    ca.jacobian_sparsity(prob.get("g"), prob.get("x")).spy();
    Function f = ca.nlpsol("solver", "fatrop", prob,
        Map.of("structure_detection", "auto", "verbose", true, "equality", equality));
    assertInAnyOutput("gap-closing", () -> f.invoke(args));
  }

  @Test
  public void test_detect() {
    Assumptions.assumeTrue(hasNlpsol("fatrop"), "nlpsol fatrop not available");
    for (long nx1 : new long[] {2})
      for (long nx2 : new long[] {2, 0})
        for (long nu1 : new long[] {2, 0})
          for (long nu2 : new long[] {2, 0})
            for (long ng1 : new long[] {2, 0})
              for (long ng2 : new long[] {2, 0})
                for (long ng3 : new long[] {2, 0}) {
                  System.out.println("test_detect " + nx1 + " " + nx2 + " " + nu1 + " " + nu2 + " " + ng1 + " " + ng2 + " " + ng3);
                  Map<String, Object> kw = new HashMap<>();
                  kw.put("N", 2L);
                  kw.put("nx0", 2L);
                  kw.put("nu0", 2L);
                  kw.put("nx1", nx1);
                  kw.put("nu1", nu1);
                  kw.put("nx2", nx2);
                  kw.put("nu2", nu2);
                  kw.put("ng1", ng1);
                  kw.put("ng2", ng2);
                  kw.put("ng3", ng3);
                  fatrop_case(kw);
                }
  }

  @Test
  public void test_detect_adversarial() {
    Assumptions.assumeTrue(hasNlpsol("fatrop"), "nlpsol fatrop not available");
    Sparsity D2 = ca.sparsify(new DM(new double[][] {{1, 0, 0}, {1, 1, 1}})).sparsity();
    fatrop_case(Map.of("nu2", 3L, "sp", Map.of("D2", D2)));

    D2 = ca.sparsify(new DM(new double[][] {{1, 0, 0}, {1, 1, 1}})).sparsity();
    Sparsity C2 = ca.sparsify(new DM(new double[][] {{0, 1}, {0, 0}})).sparsity();
    fatrop_case(Map.of("nu2", 3L, "sp", Map.of("D2", D2, "C2", C2)));
    Map<String, Object> kw = Map.of("nu2", 3L, "sp", Map.of("D2", D2, "C2", C2), "eq", Set.of("ng3"));
    assertInAnyOutput("gap-closing constraints must be like", () -> fatrop_case(kw));

    fatrop_case(Map.of("nu0", 0L, "nx0", 1L));

    fatrop_case(Map.of("nx2", 0L));

    fatrop_case(Map.of("nu2", 3L, "sp", Map.of("A1", new Sparsity(2, 2), "B1", new Sparsity(2, 2))));
    fatrop_case(Map.of("nu2", 3L, "sp", Map.of("A1", new Sparsity(2, 2))));

    fatrop_case(Map.of("nx2", 1L, "ng1", 0L, "nu1", 0L));
  }

  @Test
  public void test_bug() {
    Assumptions.assumeTrue(hasNlpsol("fatrop"), "nlpsol fatrop not available");
    MX x = MX.sym("x");

    for (String structure_detection : List.of("none", "auto")) {
      Map<String, Object> opts = Map.of("expand", true, "structure_detection", structure_detection, "equality",
          new boolean[] {true});

      Function solver = ca.nlpsol("solver", "fatrop", Map.of("x", x, "g", x.minus(1)), opts);
      assertAlmostEqual(out(solver.invoke(Map.of("lbg", 0, "ubg", 0)), "x").to_double(), 1, 5);

      solver = ca.nlpsol("solver", "fatrop", Map.of("x", x, "g", x), opts);
      assertAlmostEqual(out(solver.invoke(Map.of("lbg", 1, "ubg", 1)), "x").to_double(), 1, 5);

      solver = ca.nlpsol("solver", "fatrop", Map.of("x", x, "g", x.minus(2)), opts);
      assertAlmostEqual(out(solver.invoke(Map.of("lbg", 3, "ubg", 3)), "x").to_double(), 5, 5);
    }
  }
}

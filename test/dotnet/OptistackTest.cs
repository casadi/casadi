using System;
using System.Collections.Generic;
using System.Linq;

namespace CasADi.Tests {

/// <summary>Port of test/java OptistackTest (itself a port of test/python/optistack.py)</summary>
public class OptistackTest : CasadiTestCase {

  const string nlpsolver = "ipopt";
  static readonly Dictionary<string, object> nlpsolver_options = new Dictionary<string, object>();
  const double pi = Math.PI;
  const double inf = double.PositiveInfinity;

  // JUnit's @BeforeEach: a SkipTest thrown here skips every test of the class
  public OptistackTest() {
    requireNlpsol("ipopt");
  }

  static double tol(int digits) => Math.Pow(10, -digits);

  static Dictionary<string, object> d(params object[] kv) {
    var r = new Dictionary<string, object>();
    for (int i = 0; i < kv.Length; i += 2) r[(string)kv[i]] = kv[i + 1];
    return r;
  }

  static Dictionary<string, object> dumpOptions(params string[] fs) {
    var spec = new Dictionary<string, object>();
    foreach (string f in fs) {
      var fo = f == "nlp_f" ? d("dump_in", true, "dump_out", true) : d("dump_out", true);
      spec[f] = d("final_options", fo);
    }
    return d("specific_options", spec, "ipopt.nlp_scaling_method", "none", "ipopt.mumps_permuting_scaling", 0,
             "ipopt.mumps_scaling", 0);
  }

  static void checkfunctionLight(Function trial, Function solution, List<DM> inputs) {
    for (int k = 0; k < 2; k++) {
      List<DM> a = trial.call(inputs), b = solution.call(inputs);
      assertEquals(solution.n_out(), trial.n_out());
      assertEquals(solution.n_in(), trial.n_in());
      for (int i = 0; i < trial.n_out(); i++) {
        if (trial.sparsity_out(i).is_empty() || solution.sparsity_out(i).is_empty()) continue;
        checkarray(a[i], b[i], 1e-9, "output " + i);
      }
    }
  }

  static long iterCount(OptiSol sol) => Convert.ToInt64(sol.stats()["iter_count"]);


  // test_conic (qrqp): shadowed in Python by the later test_conic; ported under its own name
  public void test_conic_qrqp() {
    requireConic("qrqp");
    var opti = ca.Opti("conic");
    MX x = opti.variable(2);
    opti.minimize(ca.dot(x, x));
    opti.subject_to(x[0] >= 2);
    opti.subject_to(x[1] >= 3);
    opti.solver("qrqp");
    OptiSol sol = opti.solve();
    checkarray(sol.value(x), ca.vertcat(2, 3));
    assertInException("conic", () => ca.Opti("foo"));
  }

  public void test_lookup() {
    var opti = ca.Opti();
    MX x = opti.variable(2);
    MX y = opti.variable(2);
    opti.subject_to(x >= y);
    string str0 = opti.debug.x_lookup(0).ToString();
    string str1 = opti.debug.x_lookup(1).ToString();
    string str2 = opti.debug.x_lookup(2).ToString();
    string str3 = opti.debug.x_lookup(3).ToString();
    assertEquals(str0, str1);
    assertEquals(str2, str3);
    assertTrue(str1 != str2);
  }

  public void test_n() {
    var opti = ca.Opti();
    MX x = opti.variable();
    MX y = opti.variable();
    MX z = opti.variable();
    opti.minimize(ca.pow(x - ca.pow(y, 2), 2));
    opti.subject_to(x + y >= 1);
    opti.subject_to(z + x >= 1);
    opti.solver("ipopt");
    OptiSol sol = opti.solve();
    assertEquals(3L, sol.opti.nx);
    assertEquals(2L, sol.opti.ng);
    assertEquals(0L, sol.opti.np);
  }

  public void test_subject_to_list() {
    var opti = ca.Opti();
    MX x = opti.variable();
    MX y = opti.variable();
    MX z = opti.variable();
    opti.minimize(ca.pow(x - ca.pow(y, 2), 2));
    opti.subject_to(new[] { x + y >= 1, z + x >= 1 });
    opti.solver("ipopt");
    OptiSol sol = opti.solve();
    assertEquals(3L, sol.opti.nx);
    assertEquals(2L, sol.opti.ng);
    assertEquals(0L, sol.opti.np);
  }

  public void test_sol_opti() {
    var opti = ca.Opti();
    MX x = opti.variable();
    MX y = opti.variable();
    MX z = opti.variable();
    opti.minimize(ca.pow(x - ca.pow(y, 2), 2));
    opti.subject_to(x + y >= 1);
    opti.subject_to(z + x >= 1);
    opti.solver(nlpsolver, nlpsolver_options);
    OptiSol sol = opti.solve();
    _ = sol.opti.nx;
  }

  public void test_to_function() {
    var opti = ca.Opti();
    MX x = opti.variable();
    MX y = opti.variable();
    MX p = opti.parameter();
    opti.minimize(ca.pow(y, 2) + ca.pow(ca.sin(x - y - p), 2));
    opti.subject_to(x + y >= 1);
    opti.solver(nlpsolver, nlpsolver_options);
    Function F = opti.to_function("F", new[] { x, p, opti.lam_g }, new[] { x, y });
    F.call(dms(0, 0.1, 0));
    Console.WriteLine(string.Join(", ", F.call(dms(0.1, 0.1, 0))));
    Console.WriteLine(string.Join(", ", F.call(dms(0, 2, 0))));
    Console.WriteLine(string.Join(", ", F.call(dms(100, 1, 0))));

    opti = ca.Opti();
    x = opti.variable();
    y = opti.variable();
    MX z = opti.variable();
    p = opti.parameter();
    opti.minimize(ca.pow(ca.sin((x + y - p) * (2 * pi)), 2) + ca.pow(ca.sin(z * pi), 2));
    opti.subject_to(ca.eq(x, y * 2));
    opti.solver("ipopt");
    F = opti.to_function("F", new[] { ca.vertcat(x, y, z), p }, new[] { ca.vertcat(x, y, z) });
    checkarray(F.Call(ca.vertcat(3.9, 1.8, 0), 6), ca.vertcat(4, 2, 0));
    checkarray(F.Call(ca.vertcat(0, 1, 0), 6), ca.vertcat(2.0 / 3, 1.0 / 3, 0));
    checkarray(F.Call(ca.vertcat(3.9, 1.8, 7.8), 6), ca.vertcat(4, 2, 8));
    checkarray(F.Call(ca.vertcat(0, 1, 7.8), 6), ca.vertcat(2.0 / 3, 1.0 / 3, 8));

    F = opti.to_function("F", new[] { ca.vertcat(z, x, y), p }, new[] { ca.vertcat(x, y, z) });
    checkarray(F.Call(ca.vertcat(7.8, 3.9, 1.8), 6), ca.vertcat(4, 2, 8));
    checkarray(F.Call(ca.vertcat(7.8, 0, 1), 6), ca.vertcat(2.0 / 3, 1.0 / 3, 8));

    MX b = opti.variable();
    var opti2 = ca.Opti();
    MX x2 = opti2.variable();
    MX y2 = opti2.variable();
    MX z2 = opti2.variable();
    MX a = opti2.variable();
    MX p2 = opti2.parameter();
    opti2.minimize(ca.pow(ca.sin((x2 + y2 - p2) * (2 * pi)), 2) + ca.pow(ca.sin(z2 * pi), 2));
    opti2.subject_to(ca.eq(x2, y2 * 2));
    opti2.solver("ipopt");
    assertInException("Argument 0", () => opti2.to_function("F", new[] { ca.vertcat(x2, 2) }, new[] { ca.vertcat(x2, y2, z2) }));
    assertInException("Argument 1", () => opti2.to_function("F", new[] { y2, ca.vertcat(x2, 2) }, new[] { ca.vertcat(x2, y2, z2) }));
    assertInException("not independent", () => opti2.to_function("F", new[] { x2, x2 }, new[] { ca.vertcat(x2, y2, z2) }));
    opti2.to_function("F", new[] { a }, new[] { ca.vertcat(x2, y2, z2) });
    assertInException("belonging to a different instance", () => opti2.to_function("F", new[] { b }, new[] { ca.vertcat(x2, y2, z2) }));
  }

  // test_dual (first definition): shadowed in Python by the later test_dual; ported under its own name
  public void test_dual_ipopt() {
    var opti = ca.Opti();
    MX x = opti.variable();
    MX y = opti.variable();
    MX z = opti.variable();
    opti.minimize(ca.pow(x - ca.pow(y, 2), 2));
    opti.subject_to(x + y >= 1);
    opti.subject_to(z + x >= 1);
    opti.solver(nlpsolver, nlpsolver_options);
    opti.solve();
  }

  public void test_sparse() {
    var opti = ca.Opti();
    MX x = opti.variable(3, 1);
    MX p = opti.parameter();
    opti.minimize(ca.pow(p * x[1] - ca.pow(x[0], 2), 2));
    opti.subject_to(ca.eq(x[0] + x[1] * 3, 1));
    opti.solver(nlpsolver, nlpsolver_options);
    opti.set_value(p, 3);
    opti.solve();
    // isinstance(sol.value(DM.eye(3)), csc_matrix): python-only (scipy)
  }

  public void test_structure() {
    var opti = ca.Opti();
    MX x = opti.variable(3, 1);
    MX p = opti.parameter();
    opti.minimize(ca.pow(p * x[1] - ca.pow(x[0], 2), 2));
    opti.subject_to(ca.eq(x[0] + x[1] * 3, 1));
    opti.solver(nlpsolver, nlpsolver_options);
    opti.set_value(p, 3);
    OptiSol sol = opti.solve();
    Function solver = sol.opti.debug.solver();

    MX jac_g = ca.jacobian(opti.g, opti.x);
    checkarray(sol.value(jac_g), sol.value(solver.get_function("nlp_jac_g").call(new[] { opti.x, opti.p })[1]));
    MX grad_f = ca.gradient(opti.f, opti.x);
    checkarray(sol.value(grad_f), sol.value(solver.get_function("nlp_grad_f").call(new[] { opti.x, opti.p })[1]));
    DM hess_lag = sol.value(ca.hessian(opti.f + ca.dot(opti.lam_g, opti.g), opti.x).Item1);
    MX h = solver.get_function("nlp_hess_l").Call(opti.x, opti.p, 1, opti.lam_g);
    checkarray(sol.value(hess_lag), sol.value(ca.tril2symm(h.T())));
  }

  public void test_warmstart() {
    var opti = ca.Opti();
    MX x = opti.variable(2, 1);
    MX p = opti.parameter();
    opti.minimize(ca.pow(p * x[1] - ca.pow(x[0], 2), 2));
    opti.subject_to(ca.eq(x[0] + x[1] * 3, 1));
    opti.solver(nlpsolver, nlpsolver_options);

    opti.set_value(p, 3);
    OptiSol sol = opti.solve();
    OptiSol sol1 = sol;

    opti.set_value(p, 4);
    OptiSol solb = opti.solve();
    Check.Throws<AssertFailed>(() => checkarray(solb.value(x), sol1.value(x)));

    long it = iterCount(solb);
    assertTrue(it != 0);
    sol = opti.solve();
    assertTrue(iterCount(sol) == it);

    opti.set_initial(sol.value_variables());
    sol = opti.solve();
    assertTrue(iterCount(sol) != it);
    opti.set_value(sol1.value_parameters());
    sol = opti.solve();
    checkarray(sol.value(x), sol1.value(x), tol(6), "");
  }

  public void test_set_value_expr() {
    var opti = ca.Opti();
    MX p = opti.parameter(3, 1);
    MX v = opti.parameter(2, 1);
    MX x = opti.variable(2);

    opti.set_value(p, 0);
    opti.set_value(v, 0);
    opti.set_value(p[0], 3);
    checkarray(opti.debug.value(p), ca.vertcat(3, 0, 0));
    opti.set_value(p[new long[] { 0, 2 }], 2);
    checkarray(opti.debug.value(p), ca.vertcat(2, 0, 2));
    opti.set_value(p[new long[] { 0, 2 }], ca.vertcat(1, 2));
    checkarray(opti.debug.value(p), ca.vertcat(1, 0, 2));
    opti.set_value(p[new long[] { 2, 0 }], ca.vertcat(1, 2));
    checkarray(opti.debug.value(p), ca.vertcat(2, 0, 1));

    opti.set_value(ca.veccat(v, p), ca.vertcat(1, 2, 3, 4, 5));
    checkarray(opti.debug.value(p), ca.vertcat(3, 4, 5));
    checkarray(opti.debug.value(v), ca.vertcat(1, 2));

    opti.set_value(p, 0);
    opti.set_value(ca.veccat(p[0], p[0]), ca.vertcat(4, 4));
    checkarray(opti.debug.value(p), ca.vertcat(4, 0, 0));
    assertInException("ambiguous", () => opti.set_value(ca.veccat(p[0], p[0]), ca.vertcat(4, 5)));
    assertInException("cannot set a value for a variable", () => opti.set_value(ca.veccat(p, x, v), ca.vertcat(4, 5)));
    opti.set_value(p, 0);
    opti.set_value(p[0] * 3, 3);
    checkarray(opti.debug.value(p), ca.vertcat(1, 0, 0));
    assertInException("cannot set initial/value of an arbitrary expression", () => opti.set_value(p[0] + p[1], 3));

    opti.set_value(ca.veccat(v, 3), ca.vertcat(1, 2, 3));
    opti.set_value(ca.veccat(2, v, 5, p), ca.vertcat(2, 3, 4, 5, 6, 7, 8));
    opti.set_value(ca.veccat(2, v, 5, p, v), ca.vertcat(2, 3, 4, 5, 6, 7, 8, 3, 4));

    assertInException("inconsistent", () => opti.set_value(ca.veccat(v, 3), ca.vertcat(1, 2, 4)));
    assertInException("inconsistent", () => opti.set_value(ca.veccat(2, v, 5, p), ca.vertcat(5, 3, 4, 2, 6, 7, 8)));
    assertInException("ambiguous", () => opti.set_value(ca.veccat(2, v, 5, p, v), ca.vertcat(2, 3, 4, 5, 6, 7, 8, 4, 3)));
  }

  public void test_shapes() {
    var opti = ca.Opti();
    MX F = opti.variable(10, 1);
    MX x = opti.variable(2, 11);
    opti.subject_to(ca.eq(x[.., 0], x[.., 1]));
    opti.subject_to(opti.bounded(-1, F, 1));
    opti.solver(nlpsolver, nlpsolver_options);
    opti.minimize(ca.sum1(F));
    opti.solve();
  }

  public void test_symb_boundedshapes() {
    var opti = ca.Opti();
    MX y = opti.variable();
    MX x = opti.variable();
    opti.subject_to(opti.bounded(0, y, x));
    opti.solver(nlpsolver, nlpsolver_options);
    opti.minimize(ca.pow(x, 2) + ca.pow(y - 6, 2));
    opti.subject_to(x <= 3);
    OptiSol sol = opti.solve();
    checkarray(sol.value(y), ca.DM(3), tol(3), "");
  }

  public void test_symb_parametric() {
    var opti = ca.Opti();
    MX y = opti.parameter(3, 1);
    MX x = opti.parameter(3, 1);
    assertInException("!parametric[0]", () => opti.subject_to(x <= y));
  }

  public void test_callback() {
    var opti = ca.Opti();
    MX x = opti.variable();
    MX y = opti.variable();
    MX p = opti.parameter();
    opti.minimize(ca.pow(x - 1, 4) + ca.pow(y - p, 4));
    opti.solver(nlpsolver, nlpsolver_options);

    var ret = new Dictionary<string, DM>();
    opti.callback(i => ret.TryAdd("a", opti.debug.value(p)));
    opti.set_value(p, 3);
    ret.Clear();
    opti.solve();
    assertTrue((double)ret["a"] == 3);
    opti.set_value(p, 2);
    ret.Clear();
    opti.solve();
    assertTrue((double)ret["a"] == 2);
    opti.set_value(p, 3);
    ret.Clear();
    opti.solve();
    assertTrue((double)ret["a"] == 3);

    opti.callback();
    ret.Clear();
    opti.solve();
    assertFalse(ret.ContainsKey("a"));

    opti.callback(i => ret.TryAdd("b", opti.debug.value(p)));
    opti.solve();
    assertTrue((double)ret["b"] == 3);

    Opti opti2 = opti.copy();
    ret.Clear();
    opti2.solve();
    assertTrue(ret.Count == 0);
  }

  public void test_debug_value() {
    var opti = ca.Opti();
    MX p = opti.parameter();
    opti.set_value(p, 3);
    checkarray(opti.debug.value(ca.pow(p, 2)), 9);
    MX x = opti.variable();

    assertInException("This action is forbidden since you have not solved the Opti stack yet", () => opti.debug.value(ca.pow(x, 2)));
    assertInException("You cannot set a value for a variable", () => opti.set_value(x, 2));
    assertInException("You cannot set a value for a variable", () => opti.set_value(new[] { ca.eq(x, 2) }));
    opti.set_initial(x, 2);
    assertInException("You cannot set an initial value for a parameter", () => opti.set_initial(p, 2));
    assertInException("You cannot set an initial value for a parameter", () => opti.set_initial(new[] { ca.eq(p, 2) }));

    assertInException("This action is forbidden since you have not solved the Opti stack yet", () => opti.debug.value(ca.pow(x, 2)));
    checkarray(opti.debug.value(ca.pow(x, 2), opti.initial()), 4);

    MX y = opti.variable();
    assertInException("This action is forbidden since you have not solved the Opti stack yet", () => opti.debug.value(x * y));
  }

  public void test_print() {
    var opti = ca.Opti();
    Console.WriteLine(opti);
    assertTrue(opti.ToString().Contains("variables"));
    assertTrue(opti.debug.ToString().Contains("variables"));
    MX x = opti.variable();
    MX y = opti.variable();
    opti.parameter();
    Console.WriteLine(opti);
    opti.minimize(ca.pow(ca.pow(x, 2) - y, 2));
    Console.WriteLine(opti);
    opti.solver(nlpsolver, nlpsolver_options);
    Console.WriteLine(opti);
    OptiSol sol = opti.solve();
    Console.WriteLine(sol);
    Console.WriteLine(opti);

    opti = ca.Opti();
    Console.WriteLine(opti);
    x = opti.variable();
    y = opti.variable();
    opti.parameter();
    Console.WriteLine(opti);
    opti.subject_to(ca.pow(ca.pow(x, 2) - y, 2) <= 1);
    Console.WriteLine(opti);
    opti.solver(nlpsolver, nlpsolver_options);
    Console.WriteLine(opti);
    sol = opti.solve();
    Console.WriteLine(sol);
    Console.WriteLine(opti);
    assertTrue(sol.ToString().Contains("variables"));
  }

  public void test_flow() {
    var opti = ca.Opti();
    MX x = opti.variable();
    MX y = opti.variable();
    MX p = opti.parameter();
    MX w = MX.sym("w");

    opti.minimize(ca.pow(x, 2));
    opti.solver(nlpsolver, nlpsolver_options);
    OptiSol sol = opti.solve();

    assertInException("do not appear in the constraints and objective", () => sol.value(y));
    assertInException("decision variable", () => sol.value(y));
    assertInException("OptistackTest.cs", () => sol.value(y));
    assertInException("MX symbol 'w' of shape 1x1, declared outside of Opti.", () => sol.value(w));
    assertInException("MX symbol 'w' of shape 1x1, declared outside of Opti.", () => opti.subject_to(w <= 2));
    assertInException("MX symbol 'w' of shape 1x1, declared outside of Opti.", () => opti.minimize(w));

    opti.minimize(ca.pow(x - p, 2));
    assertInException("You have forgotten to assign a value to a parameter", () => opti.solve());
    assertInException("Opti parameter", () => opti.solve());
    assertInException("OptistackTest.cs", () => opti.solve());
    opti.set_value(p, 5);
    opti.solve();
    opti.debug.value(x);

    opti.set_value(p, 9);
    opti.debug.value(p);
    assertInException("This action is forbidden since you have not solved the Opti stack yet", () => opti.debug.value(x));

    opti.solve();
    opti.set_initial(x, 3);
    opti.debug.value(p);
    assertInException("This action is forbidden since you have not solved the Opti stack yet", () => opti.debug.value(x));

    MX con0 = x >= 0;
    opti.subject_to(con0);
    opti.solve();

    Opti opti0 = opti;
    var optib = ca.Opti();
    MX X = optib.variable();
    MX yb = optib.variable();
    optib.parameter();
    assertInException("belonging to a different instance of Opti", () => optib.subject_to(x <= 2));
    assertInException("decision variable", () => optib.subject_to(x <= 2));
    assertInException("belonging to a different instance of Opti", () => optib.subject_to(p <= 2));
    assertInException("parameter", () => optib.subject_to(p <= 2));
    assertInException("belonging to a different instance of Opti", () => optib.subject_to(opti0.dual(con0) <= 2));
    assertInException("dual variable", () => optib.subject_to(opti0.dual(con0) <= 2));

    MX con = X >= 0;
    optib.subject_to(con);
    MX dual0 = optib.dual(con);
    assertInException("Dual variables forbidden in this context", () => optib.subject_to(dual0 >= 2));
    assertInException("Dual variables forbidden in this context", () => optib.minimize(dual0));

    optib.subject_to();
    optib.minimize(ca.pow(yb, 2));
    MX dual = optib.dual(con);
    optib.solver(nlpsolver, nlpsolver_options);
    OptiSol solb = optib.solve();
    assertInException("This expression has a dual for a constraint that is not given to Opti", () => solb.value(dual));
    assertInException("OptistackTest.cs", () => solb.value(dual));
  }

  public void test_simple() {
    var opti = ca.Opti();
    double eps = 1e-5;
    MX x = opti.variable();
    MX y = opti.variable();
    MX f = ca.pow(x - 1, 2) + ca.pow(y - 2, 2);
    opti.minimize(f);
    opti.solver(nlpsolver, nlpsolver_options);
    OptiSol sol = opti.solve();
    checkarray(sol.value(x), 1);
    checkarray(sol.value(y), 2);

    var cases = new (MX, MX, double, double, double?)[] {
      (y >= 2.5, y + eps >= 2.5, 1.0, 2.5, 1.0),
      (y >= 1.5, y + eps >= 1.5, 1.0, 2.0, 0.0),
      (2.5 >= y, 2.5 + eps >= y, 1.0, 2.0, 0.0),
      (1.5 >= y, 1.5 + eps >= y, 1.0, 1.5, 1.0),
      (y <= 2.5, y <= 2.5 + eps, 1.0, 2.0, 0.0),
      (y <= 1.5, y <= 1.5 + eps, 1.0, 1.5, 1.0),
      (2.5 <= y, 2.5 <= y + eps, 1.0, 2.5, 1.0),
      (1.5 <= y, 1.5 <= y + eps, 1.0, 2.0, 0.0),
      (y >= x, y + eps >= x, 1.0, 2.0, 0.0),
      (y <= x, y <= x + eps, 1.5, 1.5, 1.0),
      (y <= 0, y <= eps, 1.0, 0.0, 4.0),
      ((3.0 <= y) <= 4, (3.0 <= y + eps) <= 4, 1.0, 3.0, null),
      (3.0 <= (y <= 4), 3.0 <= (y + eps <= 4), 1.0, 3.0, null),
      (0.0 <= (y <= 1), 0.0 <= (y <= 1 + eps), 1.0, 1.0, null),
      ((0.0 <= y) <= 1, (0.0 <= y) <= 1 + eps, 1.0, 1.0, null),
      ((4.0 >= y) >= 3, (4.0 >= y + eps) >= 3, 1.0, 3.0, null),
      ((1.0 >= y) >= 0, (1 + eps >= y) >= 0, 1.0, 1.0, null),
      (opti.bounded(3, y, 4), opti.bounded(3, y + eps, 4), 1.0, 3.0, null),
      (opti.bounded(0, y, 1), opti.bounded(0, y + eps, 1), 1.0, 1.0, null),
    };
    int count = 0;
    foreach (var (con, coneps, xs, ys, mul) in cases) {
      opti.subject_to();
      opti.subject_to(con);
      sol = opti.solve();
      checkarray(sol.value(x), ca.DM(xs), tol(7), "");
      checkarray(sol.value(y), ca.DM(ys), tol(7), "");

      DM dual_value = sol.value(opti.dual(con));
      opti.subject_to();
      opti.subject_to(coneps);
      OptiSol sol2 = opti.solve();

      DM s = -ca.sign(sol2.value(f) - sol.value(f));
      count++;
      if (mul != null) checkarray(dual_value, s * Math.Abs(mul.Value), tol(6), "case " + count);
    }

    opti.subject_to();
    opti.subject_to(ca.eq(1.5, y));
    sol = opti.solve();
    checkarray(sol.value(x), ca.DM(1), tol(7), "");
    checkarray(sol.value(y), ca.DM(1.5), tol(7), "");

    opti.subject_to();
    opti.subject_to(ca.eq(y, 1.5));
    sol = opti.solve();
    checkarray(sol.value(x), ca.DM(1), tol(7), "");
    checkarray(sol.value(y), ca.DM(1.5), tol(7), "");

    opti.subject_to();
    opti.subject_to(ca.eq(y, x));
    sol = opti.solve();

    MX x3 = opti.variable(3, 3);
    ca.trace(x3);
    opti.subject_to();
    opti.subject_to(x3 >= 0);
    Check.Throws<SwigException>(() => opti.solve());

    MX x2 = opti.variable(2, 1);
    opti.subject_to();
    opti.minimize(ca.pow(x2[0] - 1, 2) + ca.pow(x2[1] - 2, 2));
    opti.subject_to(x2[1] >= x2[0]);
    sol = opti.solve();
    checkarray(sol.value(x2), ca.vertcat(1, 2), tol(7), "");

    opti.subject_to();
    opti.subject_to(x2[1] <= x2[0]);
    sol = opti.solve();
    checkarray(sol.value(x2), ca.vertcat(1.5, 1.5), tol(7), "");

    opti.subject_to();
    opti.subject_to(x2 <= 0.5);
    sol = opti.solve();
    Console.WriteLine("test " + sol.value(x2));
    checkarray(sol.value(x2), ca.vertcat(0.5, 0.5), tol(7), "");
  }

  public void test_parametric() {
    var opti = ca.Opti();
    MX x = opti.variable();
    MX y = opti.variable();
    MX p = opti.parameter();
    opti.minimize(ca.pow(x - 1, 2) + ca.pow(y - p, 2));
    opti.solver(nlpsolver, nlpsolver_options);
    opti.set_value(p, 2);
    OptiSol sol = opti.solve();
    checkarray(sol.value(x), 1);
    checkarray(sol.value(y), 2);

    opti.set_value(p, 3);
    sol = opti.solve();
    checkarray(sol.value(x), 1);
    checkarray(sol.value(y), 3);

    opti.minimize(ca.pow(x - 1, 2) + ca.pow(y - 2, 2));

    var cases = new (MX, double, double, double, double)[] {
      (y >= p, 2.5, 2.5, 1.5, 2.0),
      (p >= y, 2.5, 2.0, 1.5, 1.5),
      (y <= p, 2.5, 2.0, 1.5, 1.5),
      (p <= y, 2.5, 2.5, 1.5, 2.0),
    };
    foreach (var (con, p1, y1, p2, y2) in cases) {
      opti.subject_to();
      opti.subject_to(con);
      opti.set_value(p, p1);
      sol = opti.solve();
      checkarray(sol.value(x), ca.DM(1), tol(7), "");
      checkarray(sol.value(y), ca.DM(y1), tol(7), "");
      opti.set_value(p, p2);
      sol = opti.solve();
      checkarray(sol.value(x), ca.DM(1), tol(7), "");
      checkarray(sol.value(y), ca.DM(y2), tol(7), "");
    }

    opti.subject_to();
    opti.subject_to(ca.eq(p, y));
    opti.set_value(p, 1.5);
    sol = opti.solve();
    checkarray(sol.value(x), ca.DM(1), tol(7), "");
    checkarray(sol.value(y), ca.DM(1.5), tol(7), "");

    opti.subject_to();
    opti.subject_to(ca.eq(y, p));
    opti.set_value(p, 1.5);
    sol = opti.solve();
    checkarray(sol.value(x), ca.DM(1), tol(7), "");
    checkarray(sol.value(y), ca.DM(1.5), tol(7), "");
  }

  public void test_symm() {
    var opti = ca.Opti();
    MX P = opti.variable(2, 2, "symmetric");
    DM R = new double[,] { { 1, 2 }, { 4, 4 } };
    MX f = ca.sum2(ca.sum1(ca.pow(P - R, 2)));
    MX x = ca.symvar(f)[0];
    opti.minimize(f);
    opti.solver(nlpsolver, nlpsolver_options);
    OptiSol sol = opti.solve();
    ca.Function("f", new[] { x }, new[] { f });
    checkarray(sol.value(P), new double[,] { { 1, 3 }, { 3, 4 } });

    P = opti.variable(2, 2);
    f = ca.sum2(ca.sum1(ca.pow(P - R, 2)));
    x = ca.symvar(f)[0];
    opti.minimize(f);
    sol = opti.solve();
    ca.Function("f", new[] { x }, new[] { f });
    checkarray(sol.value(P), R);
  }

  public void test_symm_parameter() {
    var opti = ca.Opti();
    MX Q = opti.parameter(2, 2, "symmetric");
    MX x = opti.variable(2);
    opti.minimize(ca.bilin(Q, x, x) - ca.sum1(x));
    opti.solver(nlpsolver, nlpsolver_options);
    DM Qv = new double[,] { { 3, 1 }, { 1, 2 } };
    opti.set_value(Q, Qv);
    OptiSol sol = opti.solve();
    checkarray(sol.value(Q), Qv);
    checkarray(sol.value(x), ca.solve(Qv * 2, DM.ones(2)), tol(6), "");
    assertInException("not even square", () => opti.parameter(2, 3, "symmetric"));
    assertInException("Unknown attribute", () => opti.parameter(2, 2, "foo"));
  }

  public void test_broadcast() {
    var opti = ca.Opti();
    MX x = opti.variable(3, 1);
    MX f = ca.mtimes(x.T(), x);
    opti.minimize(f);
    opti.subject_to(x[0] >= ca.vertcat(1, 2, 3));
    opti.solver(nlpsolver, nlpsolver_options);
    OptiSol sol = opti.solve();
    checkarray(sol.value(x[0]), ca.DM(3), tol(7), "");

    f = ca.mtimes(x.T(), x);
    opti.minimize(f);
    opti.subject_to();
    opti.subject_to(x[0] >= ca.vertcat(x[1] + 1, 0.5, 0));
    opti.solver(nlpsolver, nlpsolver_options);
    sol = opti.solve();
    checkarray(sol.value(x[0]), ca.DM(0.5), tol(4), "");
  }

  public void test_constraint_dim_mismatch() {
    var opti = ca.Opti();
    MX x = opti.variable(5, 1);
    MX p = opti.parameter(1, 1);
    MX q = opti.parameter(1, 1);
    opti.minimize(ca.dot(x, x));
    assertInException("Constraint must contain decision variables.", () => opti.subject_to(ca.eq(p, q)));
    assertInException("Constraint shape mismatch.", () => opti.subject_to(ca.eq(x[0], ca.vertcat(1, 2, 3))));
    assertInException("Constraint shape mismatch.", () => opti.subject_to(ca.eq(ca.vertcat(1, 2, 3), x[0])));
  }

  public void test_parametric_constraint() {
    var opti = ca.Opti();
    MX x = opti.variable();
    MX p = opti.parameter();
    opti.set_value(p, -1);
    opti.minimize(ca.pow(x - 1, 2));
    foreach (MX con in new[] { p <= 0, -inf <= (p <= 0), -1.0 <= (p <= 0) })
      assertInException("Constraint must contain decision variables.", () => opti.subject_to(con));
    opti.subject_to(x <= p + 3);
    opti.solver(nlpsolver, nlpsolver_options);
    OptiSol sol = opti.solve();
    checkarray(sol.value(x), ca.DM(1), tol(6), "");
    assertEquals(1L, opti.ng);
  }

  public void test_value() {
    var opti = ca.Opti();
    MX x = opti.variable();
    MX y = opti.variable();
    MX z = opti.variable();
    MX p = opti.parameter();
    MX q = opti.parameter();
    MX g = x - p >= 0;
    opti.subject_to(g);
    opti.minimize(ca.pow(x, 2) + ca.pow(y - 3, 2));
    MX lam = opti.dual(g);

    assertEquals(8.0, (double)opti.debug.value(ca.pow(p, 3), new[] { ca.eq(p, 2) }));
    assertEquals(8.0, (double)opti.debug.value(ca.pow(x, 3), new[] { ca.eq(x, 2) }));
    assertEquals(5.0, (double)opti.debug.value(x + p, new[] { ca.eq(x, 2), ca.eq(p, 3) }));
    assertEquals(8.0, (double)opti.debug.value(ca.pow(lam, 3), new[] { ca.eq(opti.lam_g, 2) }));
    assertInException("This expression depends on a parameter with unset value", () => opti.debug.value(ca.pow(p, 3)));
    assertInException("This action is forbidden since you have not solved", () => opti.debug.value(ca.pow(x, 3)));
    assertInException("This action is forbidden since you have not solved", () => opti.debug.value(ca.pow(lam, 3)));
    opti.solver("ipopt");

    opti.set_value(p, 2);
    opti.solve();
    assertEquals(4.0, (double)opti.debug.value(ca.pow(p, 2)));
    assertEquals(9.0, (double)opti.debug.value(ca.pow(p, 2), new[] { ca.eq(p, 3) }));
    checkarray(opti.debug.value(ca.pow(x, 2)), ca.DM(4), tol(5), "");
    checkarray(opti.debug.value(ca.pow(x, 2), new[] { ca.eq(x, 3) }), 9);
    checkarray(opti.debug.value(ca.pow(lam, 2)), ca.DM(16), tol(5), "");
    checkarray(opti.debug.value(ca.pow(lam, 2), new[] { ca.eq(opti.lam_g, 1) }), ca.DM(1), tol(5), "");
    assertInException("This expression has symbols that do not appear in the constraints and objective:", () => opti.debug.value(ca.pow(z, 3)));
    assertEquals(8.0, (double)opti.debug.value(ca.pow(z, 3), new[] { ca.eq(z, 2) }));
    assertInException("This expression depends on a parameter with unset value", () => opti.debug.value(q));
  }

  public void test_introspection() {
    var opti = ca.Opti();
    MX x = opti.variable();
    MX y = opti.variable();
    MX z = opti.variable();
    MX p = opti.parameter();
    opti.minimize(ca.pow(x - ca.pow(y, 2), 2));
    opti.subject_to(x + y + p >= 1);
    opti.subject_to(z + x >= 1);
    opti.solver("ipopt");
    opti.set_value(p, 3);
    OptiSol sol = opti.solve();

    Function solver = ca.nlpsol("solver", "ipopt", new Dictionary<string, MX> { ["x"] = opti.x, ["f"] = opti.f, ["g"] = opti.g, ["p"] = opti.p });
    var sol2 = solver.call(new Dictionary<string, DM> { ["p"] = sol.value(p), ["lbg"] = sol.value(opti.lbg), ["ubg"] = sol.value(opti.ubg) });
    checkarray(sol2["x"], sol.value(opti.x));
  }

  public void test_max_iter() {
    var opti = ca.Opti();
    MX x = opti.variable();
    opti.minimize(ca.pow(x - 1, 6));
    opti.solver("ipopt", d(), d("max_iter", 1));
    assertInException("Maximum_Iterations_Exceeded", () => opti.solve());

    opti.solve_limited();
    opti.solver("ipopt", d(), d("max_iter", 1000));
    opti.subject_to(x >= 1);
    opti.subject_to(x <= -1);
    assertInException("Infeasible", () => opti.solve());
    assertInException("Infeasible", () => opti.solve_limited());
  }

  public void test_conic() {
    requireConic("superscs");
    var options = d("eps", 1e-9, "do_super_scs", 1, "verbose", 0);
    var opti = ca.Opti("conic");
    MX x = opti.variable();
    MX y = opti.variable();
    foreach (double x_scale in new double[] { 1, 7 }) {
      foreach (double g_scale in new double[] { 1, 9 }) {
        opti.set_linear_scale(x, x_scale);
        MX h = ca.soc(ca.vertcat(x - 5, y - 7), 4);
        opti.subject_to(h > 0, g_scale);
        opti.minimize(x * 2 + y);
        opti.solver("superscs", d(), options);
        OptiSol sol = opti.solve();
        DM res = sol.value(ca.vertcat(x, y));
        checkarray(res, ca.vertcat(5 - 8 / Math.Sqrt(5), 7 - 4 / Math.Sqrt(5)), tol(7), "conic");
        checkarray(sol.value(opti.f), ca.DM(10 - 16 / Math.Sqrt(5) + 7 - 4 / Math.Sqrt(5)), tol(7), "conic");
      }
    }
  }

  public void test_discrete_linear() {
    requireConic("cbc");
    var opti = ca.Opti("conic");
    MX x = opti.variable();
    MX y = opti.variable();
    opti.minimize(x + y);
    opti.subject_to(y >= x * 0.5 + 1);
    opti.subject_to(y <= x * 2 + 3);
    opti.set_domain(x, "real");
    opti.solver("cbc");
    OptiSol sol = opti.solve();
    double f_real = (double)sol.value(opti.f);
    DM xsol_real = sol.value(x);
    Console.WriteLine(f_real + " " + xsol_real);
    opti.set_domain(x, "integer");
    sol = opti.solve();
    double f_integer = (double)sol.value(opti.f);
    DM xsol_integer = sol.value(x);
    Console.WriteLine(f_integer + " " + xsol_integer);
    checkarray(ca.DM(Math.Round((double)xsol_integer, MidpointRounding.AwayFromZero)), xsol_integer);
    assertTrue(f_integer > f_real);
  }

  public void test_discrete() {
    requireNlpsol("bonmin");
    var opti = ca.Opti();
    MX x = opti.variable();
    MX y = opti.variable();
    MX p = opti.parameter();
    opti.set_value(p, -1);
    opti.minimize(ca.pow(x, 2) + ca.pow(y, 2));
    opti.subject_to(x + y + p >= 0);
    opti.set_domain(x, "real");
    opti.solver("bonmin");
    OptiSol sol = opti.solve();
    double f_real = (double)sol.value(opti.f);
    sol.value(x);
    opti.set_domain(x, "integer");
    sol = opti.solve();
    double f_integer = (double)sol.value(opti.f);
    DM xsol_integer = sol.value(x);
    checkarray(ca.DM(Math.Round((double)xsol_integer, MidpointRounding.AwayFromZero)), xsol_integer);
    assertTrue(f_integer > f_real);
    Function f = opti.to_function("f", new[] { p }, new[] { x });
    checkarray(f.Call(ca.DM(-1)), xsol_integer);
  }

  public void test_linear_scale() {
    var opti = ca.Opti();
    MX x = opti.variable();
    MX y = opti.variable();
    MX z = opti.variable();
    opti.minimize(ca.pow(x - y, 2));
    opti.subject_to(x + y >= 1);
    opti.subject_to(z + x >= 1);
    opti.set_initial(x, 1);
    opti.set_initial(y, 2);
    opti.set_initial(z, 1);
    opti.solver("ipopt", dumpOptions("nlp_f", "nlp_grad_f", "nlp_jac_g"));
    OptiSol sol = opti.solve();
    DM x_sol_unscaled = sol.value(opti.x);
    DM x_unscaled = DM.from_file("nlp_f.000000.in.x.mtx");
    DM f_unscaled = DM.from_file("nlp_f.000000.out.f.mtx");
    DM grad_f_unscaled = DM.from_file("nlp_grad_f.000000.out.grad_f_x.mtx");
    DM jac_g_unscaled = DM.from_file("nlp_jac_g.000000.out.jac_g_x.mtx");

    opti = ca.Opti();
    x = opti.variable();
    y = opti.variable();
    z = opti.variable();
    opti.minimize(ca.pow(x - y, 2));
    opti.subject_to(x + y >= 1);
    opti.subject_to(z + x >= 1);
    opti.set_initial(x, 1);
    opti.set_initial(y, 2);
    opti.set_initial(z, 1);
    opti.set_linear_scale(x, 10);
    opti.set_linear_scale(y, 5, 1);
    opti.set_linear_scale(z, 3, 1);
    opti.solver("ipopt", dumpOptions("nlp_f", "nlp_grad_f"));
    sol = opti.solve();
    DM x_sol_scaled = sol.value(opti.x);
    DM x_scaled = DM.from_file("nlp_f.000000.in.x.mtx");
    DM f_scaled = DM.from_file("nlp_f.000000.out.f.mtx");
    DM grad_f_scaled = DM.from_file("nlp_grad_f.000000.out.grad_f_x.mtx");

    checkarray((x_unscaled - ca.vertcat(0, 1, 1)) / ca.vertcat(10, 5, 3), x_scaled);
    checkarray(f_unscaled, f_scaled);
    checkarray(x_sol_scaled, x_sol_unscaled);
    checkarray(grad_f_unscaled * ca.vertcat(10, 5, 3), grad_f_scaled);

    opti = ca.Opti();
    x = opti.variable();
    y = opti.variable();
    z = opti.variable();
    opti.minimize(ca.pow(x - y, 2));
    opti.subject_to(x + y >= 1, 10);
    opti.subject_to(z + x >= 1, 3);
    opti.solver("ipopt", dumpOptions("nlp_f", "nlp_jac_g"));
    sol = opti.solve();
    DM.from_file("nlp_f.000000.in.x.mtx");
    DM.from_file("nlp_f.000000.out.f.mtx");
    DM jac_g_scaled = DM.from_file("nlp_jac_g.000000.out.jac_g_x.mtx");
    checkarray(x_sol_scaled, x_sol_unscaled);
    checkarray(jac_g_scaled[0, ..] * 10, jac_g_unscaled[0, ..]);
    checkarray(jac_g_scaled[1, ..] * 3, jac_g_unscaled[1, ..]);
  }

  void linearScaleColumns(bool chained) {
    DM x_sol_unscaled = null, jac_g_unscaled = null;
    for (int pass = 0; pass < 2; pass++) {
      var opti = ca.Opti();
      MX x = opti.variable();
      MX y = opti.variable();
      MX z = opti.variable();
      opti.minimize(ca.pow(x - y - z, 2));
      if (chained) {
        opti.subject_to(z - 3 <= (ca.vertcat(x, y) <= z + 3));
      } else {
        opti.subject_to(z - 3 <= ca.vertcat(x, y));
        opti.subject_to(ca.vertcat(x, y) <= z + 3);
      }
      opti.set_initial(x, 1);
      opti.set_initial(y, 2);
      opti.set_initial(z, 1);
      if (pass == 0) {
        opti.solver("ipopt", dumpOptions("nlp_f", "nlp_grad_f", "nlp_jac_g"));
        OptiSol sol = opti.solve();
        x_sol_unscaled = sol.value(opti.x);
        DM.from_file("nlp_f.000000.in.x.mtx");
        DM.from_file("nlp_f.000000.out.f.mtx");
        DM.from_file("nlp_grad_f.000000.out.grad_f_x.mtx");
        jac_g_unscaled = DM.from_file("nlp_jac_g.000000.out.jac_g_x.mtx");
      } else {
        opti.set_linear_scale(x, 10);
        opti.set_linear_scale(y, 5, 1);
        opti.set_linear_scale(z, 3, 1);
        opti.solver("ipopt", dumpOptions("nlp_f", "nlp_jac_g"));
        OptiSol sol = opti.solve();
        DM x_sol_scaled = sol.value(opti.x);
        DM.from_file("nlp_f.000000.in.x.mtx");
        DM.from_file("nlp_f.000000.out.f.mtx");
        DM jac_g_scaled = DM.from_file("nlp_jac_g.000000.out.jac_g_x.mtx");
        checkarray(x_sol_scaled, x_sol_unscaled);
        checkarray(jac_g_scaled[.., 0], jac_g_unscaled[.., 0] * 10);
        checkarray(jac_g_scaled[.., 1], jac_g_unscaled[.., 1] * 5);
        checkarray(jac_g_scaled[.., 2], jac_g_unscaled[.., 2] * 3);
      }
    }
  }

  public void test_linear_scale2() {
    linearScaleColumns(false);
  }

  public void test_linear_scale3() {
    linearScaleColumns(true);
  }

  public void test_detect_simple_bounds() {
    string[] checks = { "f", "g", "lbg", "ubg", "lam_g" };
    var get = new Dictionary<string, Func<Opti, MX>> {
      ["f"] = o => o.f, ["g"] = o => o.g, ["lbg"] = o => o.lbg, ["ubg"] = o => o.ubg, ["lam_g"] = o => o.lam_g };
    var refv = new Dictionary<string, DM>();
    foreach (double? scale_x in new double?[] { null, 1.0, 9.0 }) {
      foreach (double? scale_g in new double?[] { null, 1.0, 100.0 }) {
        foreach (bool detect_simple_bounds in new[] { false, true }) {
          var opti = ca.Opti();
          MX x = opti.variable();
          MX y = opti.variable();
          MX z = opti.variable();
          if (scale_x != null) opti.set_linear_scale(x, scale_x.Value);
          opti.set_initial(x, 12);
          opti.set_initial(y, 12);
          opti.set_initial(z, 12);
          if (scale_g == null) opti.subject_to(-5.0 <= (x <= 7));
          else opti.subject_to(-5.0 <= (x <= 7), scale_g.Value);
          opti.subject_to(-3.0 <= (-y <= 2));
          opti.subject_to(-2.0 <= (-y <= 1));
          opti.subject_to(-3.0 <= (z * -0.3 <= 2));
          opti.subject_to(-1.0 <= (z * -0.3 <= 0));
          opti.subject_to(-4.0 <= (x <= 6));
          opti.minimize(ca.pow(x - 10, 2) + ca.pow(y - 10, 2) + ca.pow(z - 10, 2));
          opti.solver("ipopt", d("detect_simple_bounds", detect_simple_bounds));
          OptiSol sol = opti.solve();
          foreach (string k in checks) {
            if (scale_x == null && scale_g == null && !detect_simple_bounds) refv[k] = sol.value(get[k](opti));
            else checkarray(refv[k], sol.value(get[k](opti)), tol(6), k);
          }
        }
      }
    }
  }

  public void test_ipopt_custom_jac() {
    foreach (DM x_scale in new[] { ca.vertcat(1, 1, 1), ca.vertcat(7, 0.11, 0.13) }) {
      foreach (DM x_scale_offset in new[] { ca.vertcat(0, 0, 0), ca.vertcat(1, 2, 3) }) {
        foreach (DM g_scale in new[] { ca.vertcat(1, 1), ca.vertcat(3, 5) }) {
          foreach (double f_scale in new double[] { 1, 1.7 }) {
            foreach (bool scale_helper in new[] { true, false }) {
              var opti = ca.Opti();
              MX x = opti.variable();
              MX y = opti.variable();
              MX w = opti.variable();
              opti.set_linear_scale(x, x_scale[0], x_scale_offset[0]);
              opti.set_linear_scale(y, x_scale[1], x_scale_offset[1]);
              opti.set_linear_scale(w, x_scale[2], x_scale_offset[2]);
              opti.minimize(ca.pow(1 - x, 2) + ca.pow(w, 2) * 100, f_scale);
              opti.subject_to(ca.eq(y - ca.pow(x, 2) - w, 0), g_scale[0]);
              opti.subject_to(ca.eq(ca.pow(w, 2) - y, 7), g_scale[1]);
              opti.solver("ipopt");
              opti.solve();
              Function ref_solver = opti.debug.solver();

              opti = ca.Opti();
              x = opti.variable();
              y = opti.variable();
              w = opti.variable();
              opti.set_linear_scale(x, x_scale[0], x_scale_offset[0]);
              opti.set_linear_scale(y, x_scale[1], x_scale_offset[1]);
              opti.set_linear_scale(w, x_scale[2], x_scale_offset[2]);
              opti.minimize(ca.pow(1 - x, 2) + ca.pow(w, 2) * 100, f_scale);
              opti.subject_to(ca.eq(y - ca.pow(x, 2) - w, 0), g_scale[0]);
              opti.subject_to(ca.eq(ca.pow(w, 2) - y, 7), g_scale[1]);

              MX ox = opti.x, op = opti.p, og = opti.g;
              DM xls = opti.x_linear_scale, xlo = opti.x_linear_scale_offset, gls = opti.g_linear_scale;
              Function nlp_jac_g_custom;
              if (scale_helper) {
                nlp_jac_g_custom = opti.scale_helper(ca.Function("nlp_jac_g", new[] { ox, op }, new[] { og, ca.jacobian(og, ox) },
                    new[] { "x", "p" }, new[] { "g", "jac_g_x" }));
              } else {
                MX gs = og / gls;
                nlp_jac_g_custom = ca.Function("nlp_jac_g", new[] { ox, op },
                    ca.substitute(new[] { gs, ca.mtimes(ca.jacobian(gs, ox), ca.diag((MX)xls)) }, new[] { ox }, new[] { ox * xls + xlo }),
                    new[] { "x", "p" }, new[] { "g", "jac_g_x" });
              }
              opti.solver("ipopt", d("cache", d("nlp_jac_g", nlp_jac_g_custom)));
              try {
                opti.solve();
              } catch (SwigException) {
                // as the bare except in Python
              }
              Function solver = opti.debug.solver();

              Console.WriteLine(opti.x_linear_scale);
              Console.WriteLine(opti.x_linear_scale_offset);
              Console.WriteLine(opti.g_linear_scale);
              Console.WriteLine(opti.f_linear_scale);
              checkfunctionLight(ref_solver.get_function("nlp_jac_g"), solver.get_function("nlp_jac_g"), dms(ca.vertcat(0.11, 0.3, 0.7), 0));
              MX lam_f = MX.sym("lam_f");
              Function nlp_hess_l_custom;
              if (scale_helper) {
                MX lag = lam_f * opti.f + ca.dot(opti.lam_g, og);
                MX H = ca.jacobian(ca.gradient(lag, ox), ox, d("symmetric", true));
                nlp_hess_l_custom = opti.scale_helper(ca.Function("nlp_hess_l", new[] { ox, op, lam_f, opti.lam_g }, new[] { ca.triu(H) },
                    new[] { "x", "p", "lam_f", "lam_g" }, new[] { "triu_hess_gamma_x_x" }));
              } else {
                MX lag = lam_f * opti.f / f_scale + ca.dot(opti.lam_g, og / gls);
                MX H = ca.mtimes(ca.jacobian(ca.gradient(lag, ox), ox, d("symmetric", true)), (MX)ca.pow(ca.diag(xls), 2));
                H = ca.substitute(ca.triu(H), ox, ox * xls + xlo);
                nlp_hess_l_custom = ca.Function("nlp_hess_l", new[] { ox, op, lam_f, opti.lam_g }, new[] { H });
              }
              Console.WriteLine(ref_solver.get_function("nlp_hess_l"));
              checkfunctionLight(ref_solver.get_function("nlp_hess_l"), nlp_hess_l_custom, dms(ca.vertcat(0.11, 0.3, 0.7), 0, 3, 17));
            }
          }
        }
      }
    }
  }

  public void test_scaling() {
    DM f_ref = null, g_ref = null, x_ref = null, lam_g_ref = null, grad_f_ref = null, jac_g_ref = null, dual_ref = null, lbg_ref = null;
    foreach (double x_scale in new double[] { 1, 100 }) {
      foreach (double g_scale in new double[] { 1, 37 }) {
        foreach (double f_scale in new double[] { 1, 11 }) {
          Console.WriteLine(x_scale + " " + g_scale);
          var opti = ca.Opti();
          MX x = opti.variable();
          MX y = opti.variable();
          opti.variable();
          opti.set_linear_scale(x, x_scale);
          opti.minimize(ca.pow(1 - x, 2) + ca.pow(y - ca.pow(x, 2), 2) * 100, f_scale);
          MX my_g = ca.eq(ca.pow(x, 2) + ca.pow(y, 2), 1);
          opti.subject_to(my_g, g_scale);
          opti.solver("ipopt", d("ipopt.print_level", 0, "print_time", false));
          OptiSol sol = opti.solve();
          Console.WriteLine(sol.value(opti.lbg));
          if (x_scale == 1 && g_scale == 1 && f_scale == 1) {
            f_ref = sol.value(opti.f);
            g_ref = sol.value(opti.g);
            x_ref = sol.value(opti.x);
            lam_g_ref = sol.value(opti.lam_g);
            grad_f_ref = sol.value(ca.gradient(opti.f, opti.x));
            jac_g_ref = sol.value(ca.jacobian(opti.g, opti.x));
            dual_ref = sol.value(opti.dual(my_g));
            lbg_ref = sol.value(opti.lbg);
          } else {
            checkarray(f_ref, sol.value(opti.f));
            checkarray(g_ref, sol.value(opti.g));
            checkarray(x_ref, sol.value(opti.x));
            checkarray(lam_g_ref, sol.value(opti.lam_g), tol(6), "");
            checkarray(grad_f_ref, sol.value(ca.gradient(opti.f, opti.x)), tol(6), "");
            checkarray(jac_g_ref, sol.value(ca.jacobian(opti.g, opti.x)), tol(6), "");
            checkarray(dual_ref, sol.value(opti.dual(my_g)), tol(6), "");
            checkarray(lbg_ref, sol.value(opti.lbg));
          }
          DM lag_grad = sol.value(ca.gradient(opti.f + ca.dot(opti.lam_g, opti.g), opti.x));
          checkarray(sol.value(x), ca.DM(0.78641515), tol(6), "");
          checkarray(lag_grad, ca.vertcat(0, 0), tol(6), "");
        }
      }
    }
  }

  public void test_dual() {
    requireConic("highs");
    var refv = new Dictionary<string, DM>();
    Dictionary<string, DM> data = null;
    foreach (string t in new[] { "nlp", "conic" }) {
      var opti = ca.Opti(t);
      MX x = opti.variable();
      MX y = opti.variable();
      MX z = opti.variable();
      opti.minimize(ca.pow(x, 2) + ca.pow(y, 2) + ca.pow(z, 2));
      MX g1 = x * 2 + y * 3 + 4 <= 0;
      opti.subject_to(g1);
      MX g2 = ca.eq(x * 7 + z * 3 + 6, 0);
      opti.subject_to(g2);
      opti.solver(t == "conic" ? "highs" : "ipopt");
      OptiSol sol = opti.solve();
      data = t == "conic" ? new Dictionary<string, DM>() : refv;
      data["x"] = sol.value(opti.x);
      data["lam_g"] = sol.value(opti.lam_g);
      data["g1"] = sol.value(opti.dual(g1));
      data["g2"] = sol.value(opti.dual(g2));
    }
    Console.WriteLine(string.Join(", ", data));
    Console.WriteLine(string.Join(", ", refv));
    foreach (string k in data.Keys) checkarray(data[k], refv[k], tol(6), k);
  }

  public void test_dual_soc() {
    requireConic("superscs");
    foreach (string t in new[] { "conic" }) {
      var opti = ca.Opti(t);
      MX x = opti.variable();
      MX y = opti.variable();
      MX z = opti.variable();
      MX h = ca.soc(ca.vertcat(x - 5, y - 7), 4);
      opti.subject_to(h > 0);
      MX g1 = ca.eq(x + y * 2 + z * 3 + 3, 0);
      opti.subject_to(g1);
      opti.minimize(x * 2 + y);
      if (t == "conic") opti.solver("superscs", d("superscs", d("eps", 1e-9, "do_super_scs", 1, "verbose", 0)));
      else opti.solver("ipopt");
      OptiSol sol = opti.solve();
      var data = new Dictionary<string, DM> { ["f"] = sol.value(opti.f), ["x"] = sol.value(opti.x), ["lam_g"] = sol.value(opti.lam_g) };
      Console.WriteLine(string.Join(", ", data));
    }
  }

  public void test_value_sx() {
    var opti = ca.Opti();
    checkarray(opti.value(ca.SX(3)), ca.DM(3));
    assertInException("Opti works with MX", () => opti.value(SX.sym("z")));
  }

  public void test_subject_to_unknown() {
    var opti = ca.Opti();
    MX x = opti.variable();
    MX y = opti.variable();
    assertInException("Constraint type unknown", () => opti.subject_to(x + y));
  }

  public void test_optisol_no_debug() {
    assertTrue(typeof(OptiSol).GetMethod("debug") == null && typeof(OptiSol).GetProperty("debug") == null);
  }

  public void test_describe_range() {
    var opti = ca.Opti();
    MX x = opti.variable();
    MX y = opti.variable();
    opti.minimize(ca.pow(x, 2) + ca.pow(y, 2));
    assertInException("Index 5 out of range for 2 decision variables", () => opti.debug.x_describe(5));
    assertInException("Index 0 out of range for 0 constraints", () => opti.debug.g_describe(0));
  }


  public void test_issue4020() {
    var constr_funs = new Func<MX, MX, MX>[] {
      (x, y) => ca.eq(ca.pow(x, 2) - y, 0),
      (x, y) => ca.eq(0, ca.pow(x, 2) - y),
      (x, y) => ca.eq(ca.pow(x, 2), y),
      (x, y) => ca.eq(y, ca.pow(x, 2)),
    };
    foreach (var constr_fun in constr_funs) {
      var opti = ca.Opti();
      MX x = opti.variable();
      MX y = opti.variable();
      opti.minimize(y);
      MX constr = constr_fun(x, y);
      opti.subject_to(constr);
      opti.solver("ipopt");
      OptiSol sol = opti.solve();
      DM dual_opti = sol.value(opti.dual(constr));
      DM dual_all = sol.value(opti.lam_g);
      checkarray(dual_opti, dual_all);
    }
  }
}


/// <summary>OptistackTest.test_stacktrace, kept apart because it needs no solver (OptistackTest skips entirely without ipopt)</summary>
public class OptistackNoSolverTest : CasadiTestCase {
  static IDictionary<string, object> caller(Dictionary<string, object> user_dict) =>
    (IDictionary<string, object>)((System.Collections.IList)user_dict["stacktrace"])[0];

  public void test_stacktrace() {
    var opti = ca.Opti();
    MX x = opti.variable();
    MX p = opti.parameter(2);
    MX g = x >= 1;
    opti.subject_to(g);
    opti.subject_to(new[] { x <= 3 }, new Dictionary<string, object>());
    foreach (MX e in new[] { x, p, g }) {
      IDictionary<string, object> f = caller(opti.user_dict(e));
      assertEquals<object>("test_stacktrace", f["name"]);
      assertEquals<object>("OptistackTest.cs", f["file"]);
    }
    assertEquals(2L, opti.ng);
    opti.subject_to();
    assertEquals(0L, opti.ng);
  }
}

}

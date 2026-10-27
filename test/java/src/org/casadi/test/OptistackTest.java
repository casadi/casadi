package org.casadi.test;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertNotEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.function.BiFunction;

import org.casadi.DM;
import org.casadi.Function;
import org.casadi.MX;
import org.casadi.Opti;
import org.casadi.OptiCallback;
import org.casadi.OptiSol;
import org.casadi.SX;
import org.casadi.Slice;
import org.casadi.ca;
import org.junit.jupiter.api.Assumptions;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;

/** Port of test/python/optistack.py */
public class OptistackTest extends CasadiTestCase {

  private static final String nlpsolver = "ipopt";
  private static final Map<String, Object> nlpsolver_options = Map.of();
  private static final double pi = Math.PI;
  private static final double inf = Double.POSITIVE_INFINITY;

  @BeforeEach
  public void requireIpopt() {
    Assumptions.assumeTrue(hasNlpsol("ipopt"));
  }

  private static MX pow(MX x, double n) {
    return ca.power(x, n);
  }

  private static MX vertcat(MX... x) {
    return ca.vertcat(List.of(x));
  }

  private static MX veccat(MX... x) {
    return ca.veccat(List.of(x));
  }

  private static MX c(double x) {
    return new MX(x);
  }

  private static DM ix(double... i) {
    return new DM(i);
  }

  private static double tol(int digits) {
    return Math.pow(10, -digits);
  }

  private static Map<String, Object> dumpOptions(String... fs) {
    Map<String, Object> spec = new LinkedHashMap<>();
    for (String f : fs) {
      Map<String, Object> fo = f.equals("nlp_f") ? Map.of("dump_in", true, "dump_out", true) : Map.of("dump_out", true);
      spec.put(f, Map.of("final_options", fo));
    }
    Map<String, Object> o = new LinkedHashMap<>();
    o.put("specific_options", spec);
    o.put("ipopt.nlp_scaling_method", "none");
    o.put("ipopt.mumps_permuting_scaling", 0);
    o.put("ipopt.mumps_scaling", 0);
    return o;
  }

  private static void checkfunctionLight(Function trial, Function solution, List<DM> inputs) {
    for (int k = 0; k < 2; k++) {
      List<DM> a = trial.call(inputs), b = solution.call(inputs);
      assertEquals(solution.n_out(), trial.n_out());
      assertEquals(solution.n_in(), trial.n_in());
      for (int i = 0; i < trial.n_out(); i++) {
        if (trial.sparsity_out(i).is_empty() || solution.sparsity_out(i).is_empty()) continue;
        checkarray(a.get(i), b.get(i), 1e-9, "output " + i);
      }
    }
  }

  // test_conic (qrqp): shadowed in Python by the later test_conic; ported under its own name
  @Test
  public void test_conic_qrqp() {
    Assumptions.assumeTrue(hasConic("qrqp"));
    Opti opti = new Opti("conic");
    MX x = opti.variable(2);
    opti.minimize(ca.dot(x, x));
    opti.subject_to(ca.ge(x.get(0), 2.0));
    opti.subject_to(ca.ge(x.get(1), 3.0));
    opti.solver("qrqp");
    OptiSol sol = opti.solve();
    checkarray(sol.value(x), dm(2, 3));
    assertInException("conic", () -> new Opti("foo"));
  }

  @Test
  public void test_lookup() {
    Opti opti = new Opti();
    MX x = opti.variable(2);
    MX y = opti.variable(2);
    opti.subject_to(ca.ge(x, y));
    String str0 = opti.debug().x_lookup(0).toString();
    String str1 = opti.debug().x_lookup(1).toString();
    String str2 = opti.debug().x_lookup(2).toString();
    String str3 = opti.debug().x_lookup(3).toString();
    assertEquals(str0, str1);
    assertEquals(str2, str3);
    assertNotEquals(str1, str2);
  }

  @Test
  public void test_n() {
    Opti opti = new Opti();
    MX x = opti.variable();
    MX y = opti.variable();
    MX z = opti.variable();
    opti.minimize(pow(x.minus(pow(y, 2)), 2));
    opti.subject_to(ca.ge(x.plus(y), 1.0));
    opti.subject_to(ca.ge(z.plus(x), 1.0));
    opti.solver("ipopt");
    OptiSol sol = opti.solve();
    assertEquals(3, sol.opti().nx());
    assertEquals(2, sol.opti().ng());
    assertEquals(0, sol.opti().np());
  }

  @Test
  public void test_subject_to_list() {
    Opti opti = new Opti();
    MX x = opti.variable();
    MX y = opti.variable();
    MX z = opti.variable();
    opti.minimize(pow(x.minus(pow(y, 2)), 2));
    opti.subject_to(List.of(ca.ge(x.plus(y), 1.0), ca.ge(z.plus(x), 1.0)));
    opti.solver("ipopt");
    OptiSol sol = opti.solve();
    assertEquals(3, sol.opti().nx());
    assertEquals(2, sol.opti().ng());
    assertEquals(0, sol.opti().np());
  }

  @Test
  public void test_sol_opti() {
    Opti opti = new Opti();
    MX x = opti.variable();
    MX y = opti.variable();
    MX z = opti.variable();
    opti.minimize(pow(x.minus(pow(y, 2)), 2));
    opti.subject_to(ca.ge(x.plus(y), 1.0));
    opti.subject_to(ca.ge(z.plus(x), 1.0));
    opti.solver(nlpsolver, nlpsolver_options);
    OptiSol sol = opti.solve();
    sol.opti().nx();
  }

  @Test
  public void test_to_function() {
    Opti opti = new Opti();
    MX x = opti.variable();
    MX y = opti.variable();
    MX p = opti.parameter();
    opti.minimize(pow(y, 2).plus(pow(ca.sin(x.minus(y).minus(p)), 2)));
    opti.subject_to(ca.ge(x.plus(y), 1.0));
    opti.solver(nlpsolver, nlpsolver_options);
    Function F = opti.to_function("F", List.of(x, p, opti.lam_g()), List.of(x, y));
    F.invoke(0, 0.1, 0);
    System.out.println(F.invoke(0.1, 0.1, 0));
    System.out.println(F.invoke(0, 2, 0));
    System.out.println(F.invoke(100, 1, 0));

    opti = new Opti();
    x = opti.variable();
    y = opti.variable();
    MX z = opti.variable();
    p = opti.parameter();
    opti.minimize(pow(ca.sin(x.plus(y).minus(p).times(2 * pi)), 2).plus(pow(ca.sin(z.times(pi)), 2)));
    opti.subject_to(ca.eq(x, y.times(2)));
    opti.solver("ipopt");
    F = opti.to_function("F", List.of(vertcat(x, y, z), p), List.of(vertcat(x, y, z)));
    checkarray((DM) F.invoke(new double[] {3.9, 1.8, 0}, 6), dm(4, 2, 0));
    checkarray((DM) F.invoke(new double[] {0, 1, 0}, 6), dm(2.0 / 3, 1.0 / 3, 0));
    checkarray((DM) F.invoke(new double[] {3.9, 1.8, 7.8}, 6), dm(4, 2, 8));
    checkarray((DM) F.invoke(new double[] {0, 1, 7.8}, 6), dm(2.0 / 3, 1.0 / 3, 8));

    F = opti.to_function("F", List.of(vertcat(z, x, y), p), List.of(vertcat(x, y, z)));
    checkarray((DM) F.invoke(new double[] {7.8, 3.9, 1.8}, 6), dm(4, 2, 8));
    checkarray((DM) F.invoke(new double[] {7.8, 0, 1}, 6), dm(2.0 / 3, 1.0 / 3, 8));

    MX b = opti.variable();
    Opti opti2 = new Opti();
    MX x2 = opti2.variable();
    MX y2 = opti2.variable();
    MX z2 = opti2.variable();
    MX a = opti2.variable();
    MX p2 = opti2.parameter();
    opti2.minimize(pow(ca.sin(x2.plus(y2).minus(p2).times(2 * pi)), 2).plus(pow(ca.sin(z2.times(pi)), 2)));
    opti2.subject_to(ca.eq(x2, y2.times(2)));
    opti2.solver("ipopt");
    assertInException("Argument 0", () -> opti2.to_function("F", List.of(vertcat(x2, c(2))), List.of(vertcat(x2, y2, z2))));
    assertInException("Argument 1", () -> opti2.to_function("F", List.of(y2, vertcat(x2, c(2))), List.of(vertcat(x2, y2, z2))));
    assertInException("not independent", () -> opti2.to_function("F", List.of(x2, x2), List.of(vertcat(x2, y2, z2))));
    opti2.to_function("F", List.of(a), List.of(vertcat(x2, y2, z2)));
    assertInException("belonging to a different instance", () -> opti2.to_function("F", List.of(b), List.of(vertcat(x2, y2, z2))));
  }

  // test_dual (first definition): shadowed in Python by the later test_dual; ported under its own name
  @Test
  public void test_dual_ipopt() {
    Opti opti = new Opti();
    MX x = opti.variable();
    MX y = opti.variable();
    MX z = opti.variable();
    opti.minimize(pow(x.minus(pow(y, 2)), 2));
    opti.subject_to(ca.ge(x.plus(y), 1.0));
    opti.subject_to(ca.ge(z.plus(x), 1.0));
    opti.solver(nlpsolver, nlpsolver_options);
    opti.solve();
  }

  @Test
  public void test_sparse() {
    Opti opti = new Opti();
    MX x = opti.variable(3, 1);
    MX p = opti.parameter();
    opti.minimize(pow(p.times(x.get(1)).minus(pow(x.get(0), 2)), 2));
    opti.subject_to(ca.eq(x.get(0).plus(x.get(1).times(3)), 1.0));
    opti.solver(nlpsolver, nlpsolver_options);
    opti.set_value(p, 3);
    opti.solve();
    // isinstance(sol.value(DM.eye(3)), csc_matrix): python-only (scipy)
  }

  @Test
  public void test_structure() {
    Opti opti = new Opti();
    MX x = opti.variable(3, 1);
    MX p = opti.parameter();
    opti.minimize(pow(p.times(x.get(1)).minus(pow(x.get(0), 2)), 2));
    opti.subject_to(ca.eq(x.get(0).plus(x.get(1).times(3)), 1.0));
    opti.solver(nlpsolver, nlpsolver_options);
    opti.set_value(p, 3);
    OptiSol sol = opti.solve();
    Function solver = sol.opti().debug().solver();

    MX jac_g = ca.jacobian(opti.g(), opti.x());
    checkarray(sol.value(jac_g), sol.value((MX) ((List<?>) solver.get_function("nlp_jac_g").invoke(opti.x(), opti.p())).get(1)));
    MX grad_f = ca.gradient(opti.f(), opti.x());
    checkarray(sol.value(grad_f), sol.value((MX) ((List<?>) solver.get_function("nlp_grad_f").invoke(opti.x(), opti.p())).get(1)));
    DM hess_lag = sol.value((MX) ca.hessian(opti.f().plus(ca.dot(opti.lam_g(), opti.g())), opti.x()).get(0));
    MX h = (MX) solver.get_function("nlp_hess_l").invoke(opti.x(), opti.p(), 1, opti.lam_g());
    checkarray(sol.value(hess_lag), sol.value(ca.tril2symm(h.T())));
  }

  @Test
  public void test_warmstart() {
    Opti opti = new Opti();
    MX x = opti.variable(2, 1);
    MX p = opti.parameter();
    opti.minimize(pow(p.times(x.get(1)).minus(pow(x.get(0), 2)), 2));
    opti.subject_to(ca.eq(x.get(0).plus(x.get(1).times(3)), 1.0));
    opti.solver(nlpsolver, nlpsolver_options);

    opti.set_value(p, 3);
    OptiSol sol = opti.solve();
    OptiSol sol1 = sol;

    opti.set_value(p, 4);
    OptiSol solb = opti.solve();
    assertThrows(AssertionError.class, () -> checkarray(solb.value(x), sol1.value(x)));

    long it = (Long) solb.stats().get("iter_count");
    assertTrue(it != 0);
    sol = opti.solve();
    assertTrue((Long) sol.stats().get("iter_count") == it);

    opti.set_initial(sol.value_variables());
    sol = opti.solve();
    assertTrue((Long) sol.stats().get("iter_count") != it);
    opti.set_value(sol1.value_parameters());
    sol = opti.solve();
    checkarray(sol.value(x), sol1.value(x), tol(6), "");
  }

  @Test
  public void test_set_value_expr() {
    Opti opti = new Opti();
    MX p = opti.parameter(3, 1);
    MX v = opti.parameter(2, 1);
    MX x = opti.variable(2);

    opti.set_value(p, 0);
    opti.set_value(v, 0);
    opti.set_value(p.get(0), 3);
    checkarray(opti.debug().value(p), dm(3, 0, 0));
    opti.set_value(p.get(ix(0, 2)), 2);
    checkarray(opti.debug().value(p), dm(2, 0, 2));
    opti.set_value(p.get(ix(0, 2)), dm(1, 2));
    checkarray(opti.debug().value(p), dm(1, 0, 2));
    opti.set_value(p.get(ix(2, 0)), dm(1, 2));
    checkarray(opti.debug().value(p), dm(2, 0, 1));

    opti.set_value(veccat(v, p), dm(1, 2, 3, 4, 5));
    checkarray(opti.debug().value(p), dm(3, 4, 5));
    checkarray(opti.debug().value(v), dm(1, 2));

    opti.set_value(p, 0);
    opti.set_value(veccat(p.get(0), p.get(0)), dm(4, 4));
    checkarray(opti.debug().value(p), dm(4, 0, 0));
    assertInException("ambiguous", () -> opti.set_value(veccat(p.get(0), p.get(0)), dm(4, 5)));
    assertInException("cannot set a value for a variable", () -> opti.set_value(veccat(p, x, v), dm(4, 5)));
    opti.set_value(p, 0);
    opti.set_value(p.get(0).times(3), 3);
    checkarray(opti.debug().value(p), dm(1, 0, 0));
    assertInException("cannot set initial/value of an arbitrary expression", () -> opti.set_value(p.get(0).plus(p.get(1)), 3));

    opti.set_value(veccat(v, c(3)), dm(1, 2, 3));
    opti.set_value(veccat(c(2), v, c(5), p), dm(2, 3, 4, 5, 6, 7, 8));
    opti.set_value(veccat(c(2), v, c(5), p, v), dm(2, 3, 4, 5, 6, 7, 8, 3, 4));

    assertInException("inconsistent", () -> opti.set_value(veccat(v, c(3)), dm(1, 2, 4)));
    assertInException("inconsistent", () -> opti.set_value(veccat(c(2), v, c(5), p), dm(5, 3, 4, 2, 6, 7, 8)));
    assertInException("ambiguous", () -> opti.set_value(veccat(c(2), v, c(5), p, v), dm(2, 3, 4, 5, 6, 7, 8, 4, 3)));
  }

  @Test
  public void test_shapes() {
    Opti opti = new Opti();
    MX F = opti.variable(10, 1);
    MX x = opti.variable(2, 11);
    opti.subject_to(ca.eq(x.get(new Slice(), 0), x.get(new Slice(), 1)));
    opti.subject_to(Opti.bounded(-1, F, 1));
    opti.solver(nlpsolver, nlpsolver_options);
    opti.minimize(ca.sum1(F));
    opti.solve();
  }

  @Test
  public void test_symb_boundedshapes() {
    Opti opti = new Opti();
    MX y = opti.variable();
    MX x = opti.variable();
    opti.subject_to(Opti.bounded(0, y, x));
    opti.solver(nlpsolver, nlpsolver_options);
    opti.minimize(pow(x, 2).plus(pow(y.minus(6), 2)));
    opti.subject_to(ca.le(x, 3.0));
    OptiSol sol = opti.solve();
    checkarray(sol.value(y), new DM(3), tol(3), "");
  }

  @Test
  public void test_symb_parametric() {
    Opti opti = new Opti();
    MX y = opti.parameter(3, 1);
    MX x = opti.parameter(3, 1);
    assertInException("!parametric[0]", () -> opti.subject_to(ca.le(x, y)));
  }

  @Test
  public void test_callback() {
    Opti opti = new Opti();
    MX x = opti.variable();
    MX y = opti.variable();
    MX p = opti.parameter();
    opti.minimize(pow(x.minus(1), 4).plus(pow(y.minus(p), 4)));
    opti.solver(nlpsolver, nlpsolver_options);

    Map<String, DM> ret = new HashMap<>();
    OptiCallback cb = new OptiCallback() {
      @Override
      public void call(long i) {
        ret.putIfAbsent("a", opti.debug().value(p));
      }
    };
    opti.callback_class(cb);
    opti.set_value(p, 3);
    ret.clear();
    opti.solve();
    assertTrue(ret.get("a").to_double() == 3);
    opti.set_value(p, 2);
    ret.clear();
    opti.solve();
    assertTrue(ret.get("a").to_double() == 2);
    opti.set_value(p, 3);
    ret.clear();
    opti.solve();
    assertTrue(ret.get("a").to_double() == 3);

    opti.callback_class();
    ret.clear();
    opti.solve();
    assertFalse(ret.containsKey("a"));

    OptiCallback cb2 = new OptiCallback() {
      @Override
      public void call(long i) {
        ret.putIfAbsent("b", opti.debug().value(p));
      }
    };
    opti.callback_class(cb2);
    opti.solve();
    assertTrue(ret.get("b").to_double() == 3);

    Opti opti2 = opti.copy();
    ret.clear();
    opti2.solve();
    assertTrue(ret.size() == 0);
  }

  @Test
  public void test_debug_value() {
    Opti opti = new Opti();
    MX p = opti.parameter();
    opti.set_value(p, 3);
    checkarray(opti.debug().value(pow(p, 2)), 9);
    MX x = opti.variable();

    assertInException("This action is forbidden since you have not solved the Opti stack yet", () -> opti.debug().value(pow(x, 2)));
    assertInException("You cannot set a value for a variable", () -> opti.set_value(x, 2));
    assertInException("You cannot set a value for a variable", () -> opti.set_value(List.of(ca.eq(x, 2.0))));
    opti.set_initial(x, 2);
    assertInException("You cannot set an initial value for a parameter", () -> opti.set_initial(p, 2));
    assertInException("You cannot set an initial value for a parameter", () -> opti.set_initial(List.of(ca.eq(p, 2.0))));

    assertInException("This action is forbidden since you have not solved the Opti stack yet", () -> opti.debug().value(pow(x, 2)));
    checkarray(opti.debug().value(pow(x, 2), opti.initial()), 4);

    MX y = opti.variable();
    assertInException("This action is forbidden since you have not solved the Opti stack yet", () -> opti.debug().value(x.times(y)));
  }

  @Test
  public void test_print() {
    Opti opti = new Opti();
    System.out.println(opti);
    assertTrue(opti.toString().contains("variables"));
    assertTrue(opti.debug().toString().contains("variables"));
    MX x = opti.variable();
    MX y = opti.variable();
    opti.parameter();
    System.out.println(opti);
    opti.minimize(pow(pow(x, 2).minus(y), 2));
    System.out.println(opti);
    opti.solver(nlpsolver, nlpsolver_options);
    System.out.println(opti);
    OptiSol sol = opti.solve();
    System.out.println(sol);
    System.out.println(opti);

    opti = new Opti();
    System.out.println(opti);
    x = opti.variable();
    y = opti.variable();
    opti.parameter();
    System.out.println(opti);
    opti.subject_to(ca.le(pow(pow(x, 2).minus(y), 2), 1.0));
    System.out.println(opti);
    opti.solver(nlpsolver, nlpsolver_options);
    System.out.println(opti);
    sol = opti.solve();
    System.out.println(sol);
    System.out.println(opti);
    assertTrue(sol.toString().contains("variables"));
  }

  @Test
  public void test_flow() {
    Opti opti = new Opti();
    MX x = opti.variable();
    MX y = opti.variable();
    MX p = opti.parameter();
    MX w = MX.sym("w");

    opti.minimize(pow(x, 2));
    opti.solver(nlpsolver, nlpsolver_options);
    OptiSol sol = opti.solve();

    assertInException("do not appear in the constraints and objective", () -> sol.value(y));
    assertInException("decision variable", () -> sol.value(y));
    assertInException("OptistackTest.java", () -> sol.value(y));
    assertInException("MX symbol 'w' of shape 1x1, declared outside of Opti.", () -> sol.value(w));
    assertInException("MX symbol 'w' of shape 1x1, declared outside of Opti.", () -> opti.subject_to(ca.le(w, 2.0)));
    assertInException("MX symbol 'w' of shape 1x1, declared outside of Opti.", () -> opti.minimize(w));

    opti.minimize(pow(x.minus(p), 2));
    assertInException("You have forgotten to assign a value to a parameter", () -> opti.solve());
    assertInException("Opti parameter", () -> opti.solve());
    assertInException("OptistackTest.java", () -> opti.solve());
    opti.set_value(p, 5);
    opti.solve();
    opti.debug().value(x);

    opti.set_value(p, 9);
    opti.debug().value(p);
    assertInException("This action is forbidden since you have not solved the Opti stack yet", () -> opti.debug().value(x));

    opti.solve();
    opti.set_initial(x, 3);
    opti.debug().value(p);
    assertInException("This action is forbidden since you have not solved the Opti stack yet", () -> opti.debug().value(x));

    MX con0 = ca.ge(x, 0.0);
    opti.subject_to(con0);
    opti.solve();

    Opti opti0 = opti;
    Opti optib = new Opti();
    MX X = optib.variable();
    MX yb = optib.variable();
    optib.parameter();
    assertInException("belonging to a different instance of Opti", () -> optib.subject_to(ca.le(x, 2.0)));
    assertInException("decision variable", () -> optib.subject_to(ca.le(x, 2.0)));
    assertInException("belonging to a different instance of Opti", () -> optib.subject_to(ca.le(p, 2.0)));
    assertInException("parameter", () -> optib.subject_to(ca.le(p, 2.0)));
    assertInException("belonging to a different instance of Opti", () -> optib.subject_to(ca.le(opti0.dual(con0), 2.0)));
    assertInException("dual variable", () -> optib.subject_to(ca.le(opti0.dual(con0), 2.0)));

    MX con = ca.ge(X, 0.0);
    optib.subject_to(con);
    MX dual0 = optib.dual(con);
    assertInException("Dual variables forbidden in this context", () -> optib.subject_to(ca.ge(dual0, 2.0)));
    assertInException("Dual variables forbidden in this context", () -> optib.minimize(dual0));

    optib.subject_to();
    optib.minimize(pow(yb, 2));
    MX dual = optib.dual(con);
    optib.solver(nlpsolver, nlpsolver_options);
    OptiSol solb = optib.solve();
    assertInException("This expression has a dual for a constraint that is not given to Opti", () -> solb.value(dual));
    assertInException("OptistackTest.java", () -> solb.value(dual));
  }

  @Test
  public void test_simple() {
    Opti opti = new Opti();
    double eps = 1e-5;
    MX x = opti.variable();
    MX y = opti.variable();
    MX f = pow(x.minus(1), 2).plus(pow(y.minus(2), 2));
    opti.minimize(f);
    opti.solver(nlpsolver, nlpsolver_options);
    OptiSol sol = opti.solve();
    checkarray(sol.value(x), 1);
    checkarray(sol.value(y), 2);

    Object[][] cases = {
      {ca.ge(y, 2.5), ca.ge(y.plus(eps), 2.5), 1.0, 2.5, 1.0},
      {ca.ge(y, 1.5), ca.ge(y.plus(eps), 1.5), 1.0, 2.0, 0.0},
      {ca.ge(2.5, y), ca.ge(2.5 + eps, y), 1.0, 2.0, 0.0},
      {ca.ge(1.5, y), ca.ge(1.5 + eps, y), 1.0, 1.5, 1.0},
      {ca.le(y, 2.5), ca.le(y, 2.5 + eps), 1.0, 2.0, 0.0},
      {ca.le(y, 1.5), ca.le(y, 1.5 + eps), 1.0, 1.5, 1.0},
      {ca.le(2.5, y), ca.le(2.5, y.plus(eps)), 1.0, 2.5, 1.0},
      {ca.le(1.5, y), ca.le(1.5, y.plus(eps)), 1.0, 2.0, 0.0},
      {ca.ge(y, x), ca.ge(y.plus(eps), x), 1.0, 2.0, 0.0},
      {ca.le(y, x), ca.le(y, x.plus(eps)), 1.5, 1.5, 1.0},
      {ca.le(y, 0.0), ca.le(y, eps), 1.0, 0.0, 4.0},
      {ca.le(ca.le(3.0, y), 4.0), ca.le(ca.le(3.0, y.plus(eps)), 4.0), 1.0, 3.0, null},
      {ca.le(3.0, ca.le(y, 4.0)), ca.le(3.0, ca.le(y.plus(eps), 4.0)), 1.0, 3.0, null},
      {ca.le(0.0, ca.le(y, 1.0)), ca.le(0.0, ca.le(y, 1 + eps)), 1.0, 1.0, null},
      {ca.le(ca.le(0.0, y), 1.0), ca.le(ca.le(0.0, y), 1 + eps), 1.0, 1.0, null},
      {ca.ge(ca.ge(4.0, y), 3.0), ca.ge(ca.ge(4.0, y.plus(eps)), 3.0), 1.0, 3.0, null},
      {ca.ge(ca.ge(1.0, y), 0.0), ca.ge(ca.ge(1 + eps, y), 0.0), 1.0, 1.0, null},
      {Opti.bounded(3, y, 4), Opti.bounded(3, y.plus(eps), 4), 1.0, 3.0, null},
      {Opti.bounded(0, y, 1), Opti.bounded(0, y.plus(eps), 1), 1.0, 1.0, null},
    };
    int count = 0;
    for (Object[] cs : cases) {
      MX con = (MX) cs[0], coneps = (MX) cs[1];
      double xs = (Double) cs[2], ys = (Double) cs[3];
      Double mul = (Double) cs[4];
      opti.subject_to();
      opti.subject_to(con);
      sol = opti.solve();
      checkarray(sol.value(x), new DM(xs), tol(7), "");
      checkarray(sol.value(y), new DM(ys), tol(7), "");

      DM dual_value = sol.value(opti.dual(con));
      opti.subject_to();
      opti.subject_to(coneps);
      OptiSol sol2 = opti.solve();

      DM s = ca.sign(sol2.value(f).minus(sol.value(f))).unaryMinus();
      count++;
      if (mul != null) checkarray(dual_value, s.times(Math.abs(mul)), tol(6), "case " + count);
    }

    opti.subject_to();
    opti.subject_to(ca.eq(1.5, y));
    sol = opti.solve();
    checkarray(sol.value(x), new DM(1), tol(7), "");
    checkarray(sol.value(y), new DM(1.5), tol(7), "");

    opti.subject_to();
    opti.subject_to(ca.eq(y, 1.5));
    sol = opti.solve();
    checkarray(sol.value(x), new DM(1), tol(7), "");
    checkarray(sol.value(y), new DM(1.5), tol(7), "");

    opti.subject_to();
    opti.subject_to(ca.eq(y, x));
    sol = opti.solve();

    MX x3 = opti.variable(3, 3);
    ca.trace(x3);
    opti.subject_to();
    opti.subject_to(ca.ge(x3, 0.0));
    assertThrows(RuntimeException.class, () -> opti.solve());

    MX x2 = opti.variable(2, 1);
    opti.subject_to();
    opti.minimize(pow(x2.get(0).minus(1), 2).plus(pow(x2.get(1).minus(2), 2)));
    opti.subject_to(ca.ge(x2.get(1), x2.get(0)));
    sol = opti.solve();
    checkarray(sol.value(x2), dm(1, 2), tol(7), "");

    opti.subject_to();
    opti.subject_to(ca.le(x2.get(1), x2.get(0)));
    sol = opti.solve();
    checkarray(sol.value(x2), dm(1.5, 1.5), tol(7), "");

    opti.subject_to();
    opti.subject_to(ca.le(x2, 0.5));
    sol = opti.solve();
    System.out.println("test " + sol.value(x2));
    checkarray(sol.value(x2), dm(0.5, 0.5), tol(7), "");
  }

  @Test
  public void test_parametric() {
    Opti opti = new Opti();
    MX x = opti.variable();
    MX y = opti.variable();
    MX p = opti.parameter();
    opti.minimize(pow(x.minus(1), 2).plus(pow(y.minus(p), 2)));
    opti.solver(nlpsolver, nlpsolver_options);
    opti.set_value(p, 2);
    OptiSol sol = opti.solve();
    checkarray(sol.value(x), 1);
    checkarray(sol.value(y), 2);

    opti.set_value(p, 3);
    sol = opti.solve();
    checkarray(sol.value(x), 1);
    checkarray(sol.value(y), 3);

    opti.minimize(pow(x.minus(1), 2).plus(pow(y.minus(2), 2)));

    Object[][] cases = {
      {ca.ge(y, p), 2.5, 2.5, 1.5, 2.0},
      {ca.ge(p, y), 2.5, 2.0, 1.5, 1.5},
      {ca.le(y, p), 2.5, 2.0, 1.5, 1.5},
      {ca.le(p, y), 2.5, 2.5, 1.5, 2.0},
    };
    for (Object[] cs : cases) {
      opti.subject_to();
      opti.subject_to((MX) cs[0]);
      opti.set_value(p, (Double) cs[1]);
      sol = opti.solve();
      checkarray(sol.value(x), new DM(1), tol(7), "");
      checkarray(sol.value(y), new DM((Double) cs[2]), tol(7), "");
      opti.set_value(p, (Double) cs[3]);
      sol = opti.solve();
      checkarray(sol.value(x), new DM(1), tol(7), "");
      checkarray(sol.value(y), new DM((Double) cs[4]), tol(7), "");
    }

    opti.subject_to();
    opti.subject_to(ca.eq(p, y));
    opti.set_value(p, 1.5);
    sol = opti.solve();
    checkarray(sol.value(x), new DM(1), tol(7), "");
    checkarray(sol.value(y), new DM(1.5), tol(7), "");

    opti.subject_to();
    opti.subject_to(ca.eq(y, p));
    opti.set_value(p, 1.5);
    sol = opti.solve();
    checkarray(sol.value(x), new DM(1), tol(7), "");
    checkarray(sol.value(y), new DM(1.5), tol(7), "");
  }

  @Test
  public void test_symm() {
    Opti opti = new Opti();
    MX P = opti.variable(2, 2, "symmetric");
    DM R = dm(new double[][] {{1, 2}, {4, 4}});
    MX f = ca.sum2(ca.sum1(pow(P.minus(R), 2)));
    MX x = ca.symvar(f).get(0);
    opti.minimize(f);
    opti.solver(nlpsolver, nlpsolver_options);
    OptiSol sol = opti.solve();
    new Function("f", List.of(x), List.of(f));
    checkarray(sol.value(P), dm(new double[][] {{1, 3}, {3, 4}}));

    P = opti.variable(2, 2);
    f = ca.sum2(ca.sum1(pow(P.minus(R), 2)));
    x = ca.symvar(f).get(0);
    opti.minimize(f);
    sol = opti.solve();
    new Function("f", List.of(x), List.of(f));
    checkarray(sol.value(P), R);
  }

  @Test
  public void test_symm_parameter() {
    Opti opti = new Opti();
    MX Q = opti.parameter(2, 2, "symmetric");
    MX x = opti.variable(2);
    opti.minimize(ca.bilin(Q, x, x).minus(ca.sum1(x)));
    opti.solver(nlpsolver, nlpsolver_options);
    DM Qv = dm(new double[][] {{3, 1}, {1, 2}});
    opti.set_value(Q, Qv);
    OptiSol sol = opti.solve();
    checkarray(sol.value(Q), Qv);
    checkarray(sol.value(x), ca.solve(Qv.times(2), DM.ones(2)), tol(6), "");
    assertInException("not even square", () -> opti.parameter(2, 3, "symmetric"));
    assertInException("Unknown attribute", () -> opti.parameter(2, 2, "foo"));
  }

  @Test
  public void test_broadcast() {
    Opti opti = new Opti();
    MX x = opti.variable(3, 1);
    MX f = ca.mtimes(x.T(), x);
    opti.minimize(f);
    opti.subject_to(ca.ge(x.get(0), dm(1, 2, 3)));
    opti.solver(nlpsolver, nlpsolver_options);
    OptiSol sol = opti.solve();
    checkarray(sol.value(x.get(0)), new DM(3), tol(7), "");

    f = ca.mtimes(x.T(), x);
    opti.minimize(f);
    opti.subject_to();
    opti.subject_to(ca.ge(x.get(0), vertcat(x.get(1).plus(1), c(0.5), c(0))));
    opti.solver(nlpsolver, nlpsolver_options);
    sol = opti.solve();
    checkarray(sol.value(x.get(0)), new DM(0.5), tol(4), "");
  }

  @Test
  public void test_constraint_dim_mismatch() {
    Opti opti = new Opti();
    MX x = opti.variable(5, 1);
    MX p = opti.parameter(1, 1);
    MX q = opti.parameter(1, 1);
    opti.minimize(ca.dot(x, x));
    assertInException("Constraint must contain decision variables.", () -> opti.subject_to(ca.eq(p, q)));
    assertInException("Constraint shape mismatch.", () -> opti.subject_to(ca.eq(x.get(0), dm(1, 2, 3))));
    assertInException("Constraint shape mismatch.", () -> opti.subject_to(ca.eq(dm(1, 2, 3), x.get(0))));
  }

  @Test
  public void test_parametric_constraint() {
    Opti opti = new Opti();
    MX x = opti.variable();
    MX p = opti.parameter();
    opti.set_value(p, -1);
    opti.minimize(pow(x.minus(1), 2));
    for (MX con : List.of(ca.le(p, 0.0), ca.le(-inf, ca.le(p, 0.0)), ca.le(-1.0, ca.le(p, 0.0))))
      assertInException("Constraint must contain decision variables.", () -> opti.subject_to(con));
    opti.subject_to(ca.le(x, p.plus(3)));
    opti.solver(nlpsolver, nlpsolver_options);
    OptiSol sol = opti.solve();
    checkarray(sol.value(x), new DM(1), tol(6), "");
    assertEquals(1, opti.ng());
  }

  @Test
  public void test_value() {
    Opti opti = new Opti();
    MX x = opti.variable();
    MX y = opti.variable();
    MX z = opti.variable();
    MX p = opti.parameter();
    MX q = opti.parameter();
    MX g = ca.ge(x.minus(p), 0.0);
    opti.subject_to(g);
    opti.minimize(pow(x, 2).plus(pow(y.minus(3), 2)));
    MX lam = opti.dual(g);

    assertEquals(8.0, opti.debug().value(pow(p, 3), List.of(ca.eq(p, 2.0))).to_double());
    assertEquals(8.0, opti.debug().value(pow(x, 3), List.of(ca.eq(x, 2.0))).to_double());
    assertEquals(5.0, opti.debug().value(x.plus(p), List.of(ca.eq(x, 2.0), ca.eq(p, 3.0))).to_double());
    assertEquals(8.0, opti.debug().value(pow(lam, 3), List.of(ca.eq(opti.lam_g(), 2.0))).to_double());
    assertInException("This expression depends on a parameter with unset value", () -> opti.debug().value(pow(p, 3)));
    assertInException("This action is forbidden since you have not solved", () -> opti.debug().value(pow(x, 3)));
    assertInException("This action is forbidden since you have not solved", () -> opti.debug().value(pow(lam, 3)));
    opti.solver("ipopt");

    opti.set_value(p, 2);
    opti.solve();
    assertEquals(4.0, opti.debug().value(pow(p, 2)).to_double());
    assertEquals(9.0, opti.debug().value(pow(p, 2), List.of(ca.eq(p, 3.0))).to_double());
    checkarray(opti.debug().value(pow(x, 2)), new DM(4), tol(5), "");
    checkarray(opti.debug().value(pow(x, 2), List.of(ca.eq(x, 3.0))), 9);
    checkarray(opti.debug().value(pow(lam, 2)), new DM(16), tol(5), "");
    checkarray(opti.debug().value(pow(lam, 2), List.of(ca.eq(opti.lam_g(), 1.0))), new DM(1), tol(5), "");
    assertInException("This expression has symbols that do not appear in the constraints and objective:", () -> opti.debug().value(pow(z, 3)));
    assertEquals(8.0, opti.debug().value(pow(z, 3), List.of(ca.eq(z, 2.0))).to_double());
    assertInException("This expression depends on a parameter with unset value", () -> opti.debug().value(q));
  }

  @Test
  public void test_introspection() {
    Opti opti = new Opti();
    MX x = opti.variable();
    MX y = opti.variable();
    MX z = opti.variable();
    MX p = opti.parameter();
    opti.minimize(pow(x.minus(pow(y, 2)), 2));
    opti.subject_to(ca.ge(x.plus(y).plus(p), 1.0));
    opti.subject_to(ca.ge(z.plus(x), 1.0));
    opti.solver("ipopt");
    opti.set_value(p, 3);
    OptiSol sol = opti.solve();

    Function solver = ca.nlpsol("solver", "ipopt", Map.of("x", opti.x(), "f", opti.f(), "g", opti.g(), "p", opti.p()));
    Map<String, DM> sol2 = solver.call(Map.of("p", sol.value(p), "lbg", sol.value(opti.lbg()), "ubg", sol.value(opti.ubg())));
    checkarray(sol2.get("x"), sol.value(opti.x()));
  }

  @Test
  public void test_max_iter() {
    Opti opti = new Opti();
    MX x = opti.variable();
    opti.minimize(pow(x.minus(1), 6));
    opti.solver("ipopt", Map.of(), Map.of("max_iter", 1));
    assertInException("Maximum_Iterations_Exceeded", () -> opti.solve());

    opti.solve_limited();
    opti.solver("ipopt", Map.of(), Map.of("max_iter", 1000));
    opti.subject_to(ca.ge(x, 1.0));
    opti.subject_to(ca.le(x, -1.0));
    assertInException("Infeasible", () -> opti.solve());
    assertInException("Infeasible", () -> opti.solve_limited());
  }

  @Test
  public void test_conic() {
    Assumptions.assumeTrue(hasConic("superscs"));
    Map<String, Object> options = Map.of("eps", 1e-9, "do_super_scs", 1, "verbose", 0);
    Opti opti = new Opti("conic");
    MX x = opti.variable();
    MX y = opti.variable();
    for (double x_scale : new double[] {1, 7}) {
      for (double g_scale : new double[] {1, 9}) {
        opti.set_linear_scale(x, x_scale);
        MX h = ca.soc(vertcat(x.minus(5), y.minus(7)), c(4));
        opti.subject_to(ca.gt(h, 0.0), g_scale);
        opti.minimize(x.times(2).plus(y));
        opti.solver("superscs", Map.of(), options);
        OptiSol sol = opti.solve();
        DM res = sol.value(vertcat(x, y));
        checkarray(res, dm(5 - 8 / Math.sqrt(5), 7 - 4 / Math.sqrt(5)), tol(7), "conic");
        checkarray(sol.value(opti.f()), new DM(10 - 16 / Math.sqrt(5) + 7 - 4 / Math.sqrt(5)), tol(7), "conic");
      }
    }
  }

  @Test
  public void test_discrete_linear() {
    Assumptions.assumeTrue(hasConic("cbc"));
    Opti opti = new Opti("conic");
    MX x = opti.variable();
    MX y = opti.variable();
    opti.minimize(x.plus(y));
    opti.subject_to(ca.ge(y, x.times(0.5).plus(1)));
    opti.subject_to(ca.le(y, x.times(2).plus(3)));
    opti.set_domain(x, "real");
    opti.solver("cbc");
    OptiSol sol = opti.solve();
    double f_real = sol.value(opti.f()).to_double();
    DM xsol_real = sol.value(x);
    System.out.println(f_real + " " + xsol_real);
    opti.set_domain(x, "integer");
    sol = opti.solve();
    double f_integer = sol.value(opti.f()).to_double();
    DM xsol_integer = sol.value(x);
    System.out.println(f_integer + " " + xsol_integer);
    checkarray(new DM(Math.round(xsol_integer.to_double())), xsol_integer);
    assertTrue(f_integer > f_real);
  }

  @Test
  public void test_discrete() {
    Assumptions.assumeTrue(hasNlpsol("bonmin"));
    Opti opti = new Opti();
    MX x = opti.variable();
    MX y = opti.variable();
    MX p = opti.parameter();
    opti.set_value(p, -1);
    opti.minimize(pow(x, 2).plus(pow(y, 2)));
    opti.subject_to(ca.ge(x.plus(y).plus(p), 0.0));
    opti.set_domain(x, "real");
    opti.solver("bonmin");
    OptiSol sol = opti.solve();
    double f_real = sol.value(opti.f()).to_double();
    sol.value(x);
    opti.set_domain(x, "integer");
    sol = opti.solve();
    double f_integer = sol.value(opti.f()).to_double();
    DM xsol_integer = sol.value(x);
    checkarray(new DM(Math.round(xsol_integer.to_double())), xsol_integer);
    assertTrue(f_integer > f_real);
    Function f = opti.to_function("f", List.of(p), List.of(x));
    checkarray((DM) f.invoke(-1), xsol_integer);
  }

  @Test
  public void test_linear_scale() {
    Opti opti = new Opti();
    MX x = opti.variable();
    MX y = opti.variable();
    MX z = opti.variable();
    opti.minimize(pow(x.minus(y), 2));
    opti.subject_to(ca.ge(x.plus(y), 1.0));
    opti.subject_to(ca.ge(z.plus(x), 1.0));
    opti.set_initial(x, 1);
    opti.set_initial(y, 2);
    opti.set_initial(z, 1);
    opti.solver("ipopt", dumpOptions("nlp_f", "nlp_grad_f", "nlp_jac_g"));
    OptiSol sol = opti.solve();
    DM x_sol_unscaled = sol.value(opti.x());
    DM x_unscaled = DM.from_file("nlp_f.000000.in.x.mtx");
    DM f_unscaled = DM.from_file("nlp_f.000000.out.f.mtx");
    DM grad_f_unscaled = DM.from_file("nlp_grad_f.000000.out.grad_f_x.mtx");
    DM jac_g_unscaled = DM.from_file("nlp_jac_g.000000.out.jac_g_x.mtx");

    opti = new Opti();
    x = opti.variable();
    y = opti.variable();
    z = opti.variable();
    opti.minimize(pow(x.minus(y), 2));
    opti.subject_to(ca.ge(x.plus(y), 1.0));
    opti.subject_to(ca.ge(z.plus(x), 1.0));
    opti.set_initial(x, 1);
    opti.set_initial(y, 2);
    opti.set_initial(z, 1);
    opti.set_linear_scale(x, 10);
    opti.set_linear_scale(y, 5, 1);
    opti.set_linear_scale(z, 3, 1);
    opti.solver("ipopt", dumpOptions("nlp_f", "nlp_grad_f"));
    sol = opti.solve();
    DM x_sol_scaled = sol.value(opti.x());
    DM x_scaled = DM.from_file("nlp_f.000000.in.x.mtx");
    DM f_scaled = DM.from_file("nlp_f.000000.out.f.mtx");
    DM grad_f_scaled = DM.from_file("nlp_grad_f.000000.out.grad_f_x.mtx");

    checkarray(x_unscaled.minus(dm(0, 1, 1)).div(dm(10, 5, 3)), x_scaled);
    checkarray(f_unscaled, f_scaled);
    checkarray(x_sol_scaled, x_sol_unscaled);
    checkarray(grad_f_unscaled.times(dm(10, 5, 3)), grad_f_scaled);

    opti = new Opti();
    x = opti.variable();
    y = opti.variable();
    z = opti.variable();
    opti.minimize(pow(x.minus(y), 2));
    opti.subject_to(ca.ge(x.plus(y), 1.0), 10);
    opti.subject_to(ca.ge(z.plus(x), 1.0), 3);
    opti.solver("ipopt", dumpOptions("nlp_f", "nlp_jac_g"));
    sol = opti.solve();
    DM.from_file("nlp_f.000000.in.x.mtx");
    DM.from_file("nlp_f.000000.out.f.mtx");
    DM jac_g_scaled = DM.from_file("nlp_jac_g.000000.out.jac_g_x.mtx");
    checkarray(x_sol_scaled, x_sol_unscaled);
    checkarray(jac_g_scaled.get(0, new Slice()).times(10), jac_g_unscaled.get(0, new Slice()));
    checkarray(jac_g_scaled.get(1, new Slice()).times(3), jac_g_unscaled.get(1, new Slice()));
  }

  private void linearScaleColumns(boolean chained) {
    DM x_sol_unscaled = null, jac_g_unscaled = null;
    for (int pass = 0; pass < 2; pass++) {
      Opti opti = new Opti();
      MX x = opti.variable();
      MX y = opti.variable();
      MX z = opti.variable();
      opti.minimize(pow(x.minus(y).minus(z), 2));
      if (chained) {
        opti.subject_to(ca.le(z.minus(3), ca.le(vertcat(x, y), z.plus(3))));
      } else {
        opti.subject_to(ca.le(z.minus(3), vertcat(x, y)));
        opti.subject_to(ca.le(vertcat(x, y), z.plus(3)));
      }
      opti.set_initial(x, 1);
      opti.set_initial(y, 2);
      opti.set_initial(z, 1);
      if (pass == 0) {
        opti.solver("ipopt", dumpOptions("nlp_f", "nlp_grad_f", "nlp_jac_g"));
        OptiSol sol = opti.solve();
        x_sol_unscaled = sol.value(opti.x());
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
        DM x_sol_scaled = sol.value(opti.x());
        DM.from_file("nlp_f.000000.in.x.mtx");
        DM.from_file("nlp_f.000000.out.f.mtx");
        DM jac_g_scaled = DM.from_file("nlp_jac_g.000000.out.jac_g_x.mtx");
        checkarray(x_sol_scaled, x_sol_unscaled);
        checkarray(jac_g_scaled.get(new Slice(), 0), jac_g_unscaled.get(new Slice(), 0).times(10));
        checkarray(jac_g_scaled.get(new Slice(), 1), jac_g_unscaled.get(new Slice(), 1).times(5));
        checkarray(jac_g_scaled.get(new Slice(), 2), jac_g_unscaled.get(new Slice(), 2).times(3));
      }
    }
  }

  @Test
  public void test_linear_scale2() {
    linearScaleColumns(false);
  }

  @Test
  public void test_linear_scale3() {
    linearScaleColumns(true);
  }

  @Test
  public void test_detect_simple_bounds() {
    List<String> checks = List.of("f", "g", "lbg", "ubg", "lam_g");
    Map<String, java.util.function.Function<Opti, MX>> get = Map.of(
        "f", Opti::f, "g", Opti::g, "lbg", Opti::lbg, "ubg", Opti::ubg, "lam_g", Opti::lam_g);
    Map<String, DM> ref = new HashMap<>();
    for (Double scale_x : new Double[] {null, 1.0, 9.0}) {
      for (Double scale_g : new Double[] {null, 1.0, 100.0}) {
        for (boolean detect_simple_bounds : new boolean[] {false, true}) {
          Opti opti = new Opti();
          MX x = opti.variable();
          MX y = opti.variable();
          MX z = opti.variable();
          if (scale_x != null) opti.set_linear_scale(x, scale_x);
          opti.set_initial(x, 12);
          opti.set_initial(y, 12);
          opti.set_initial(z, 12);
          if (scale_g == null) opti.subject_to(ca.le(-5.0, ca.le(x, 7.0)));
          else opti.subject_to(ca.le(-5.0, ca.le(x, 7.0)), scale_g);
          opti.subject_to(ca.le(-3.0, ca.le(y.unaryMinus(), 2.0)));
          opti.subject_to(ca.le(-2.0, ca.le(y.unaryMinus(), 1.0)));
          opti.subject_to(ca.le(-3.0, ca.le(z.times(-0.3), 2.0)));
          opti.subject_to(ca.le(-1.0, ca.le(z.times(-0.3), 0.0)));
          opti.subject_to(ca.le(-4.0, ca.le(x, 6.0)));
          opti.minimize(pow(x.minus(10), 2).plus(pow(y.minus(10), 2)).plus(pow(z.minus(10), 2)));
          opti.solver("ipopt", Map.of("detect_simple_bounds", detect_simple_bounds));
          OptiSol sol = opti.solve();
          for (String k : checks) {
            if (scale_x == null && scale_g == null && !detect_simple_bounds) ref.put(k, sol.value(get.get(k).apply(opti)));
            else checkarray(ref.get(k), sol.value(get.get(k).apply(opti)), tol(6), k);
          }
        }
      }
    }
  }

  @Test
  public void test_ipopt_custom_jac() {
    for (DM x_scale : List.of(dm(1, 1, 1), dm(7, 0.11, 0.13))) {
      for (DM x_scale_offset : List.of(dm(0, 0, 0), dm(1, 2, 3))) {
        for (DM g_scale : List.of(dm(1, 1), dm(3, 5))) {
          for (double f_scale : new double[] {1, 1.7}) {
            for (boolean scale_helper : new boolean[] {true, false}) {
              Opti opti = new Opti();
              MX x = opti.variable();
              MX y = opti.variable();
              MX w = opti.variable();
              opti.set_linear_scale(x, x_scale.get(0), x_scale_offset.get(0));
              opti.set_linear_scale(y, x_scale.get(1), x_scale_offset.get(1));
              opti.set_linear_scale(w, x_scale.get(2), x_scale_offset.get(2));
              opti.minimize(pow(c(1).minus(x), 2).plus(pow(w, 2).times(100)), f_scale);
              opti.subject_to(ca.eq(y.minus(pow(x, 2)).minus(w), 0.0), g_scale.get(0));
              opti.subject_to(ca.eq(pow(w, 2).minus(y), 7.0), g_scale.get(1));
              opti.solver("ipopt");
              opti.solve();
              Function ref_solver = opti.debug().solver();

              opti = new Opti();
              x = opti.variable();
              y = opti.variable();
              w = opti.variable();
              opti.set_linear_scale(x, x_scale.get(0), x_scale_offset.get(0));
              opti.set_linear_scale(y, x_scale.get(1), x_scale_offset.get(1));
              opti.set_linear_scale(w, x_scale.get(2), x_scale_offset.get(2));
              opti.minimize(pow(c(1).minus(x), 2).plus(pow(w, 2).times(100)), f_scale);
              opti.subject_to(ca.eq(y.minus(pow(x, 2)).minus(w), 0.0), g_scale.get(0));
              opti.subject_to(ca.eq(pow(w, 2).minus(y), 7.0), g_scale.get(1));

              MX ox = opti.x(), op = opti.p(), og = opti.g();
              DM xls = opti.x_linear_scale(), xlo = opti.x_linear_scale_offset(), gls = opti.g_linear_scale();
              Function nlp_jac_g_custom;
              if (scale_helper) {
                nlp_jac_g_custom = opti.scale_helper(new Function("nlp_jac_g", List.of(ox, op), List.of(og, ca.jacobian(og, ox)), List.of("x", "p"), List.of("g", "jac_g_x")));
              } else {
                MX gs = og.div(gls);
                nlp_jac_g_custom = new Function("nlp_jac_g", List.of(ox, op),
                    ca.substitute(List.of(gs, ca.mtimes(ca.jacobian(gs, ox), ca.diag(new MX(xls)))), List.of(ox), List.of(ox.times(xls).plus(xlo))),
                    List.of("x", "p"), List.of("g", "jac_g_x"));
              }
              Map<String, Object> options = new HashMap<>();
              options.put("cache", Map.of("nlp_jac_g", nlp_jac_g_custom));
              opti.solver("ipopt", options);
              try {
                opti.solve();
              } catch (RuntimeException e) {
                // as the bare except in Python
              }
              Function solver = opti.debug().solver();

              System.out.println(opti.x_linear_scale());
              System.out.println(opti.x_linear_scale_offset());
              System.out.println(opti.g_linear_scale());
              System.out.println(opti.f_linear_scale());
              checkfunctionLight(ref_solver.get_function("nlp_jac_g"), solver.get_function("nlp_jac_g"), dms(dm(0.11, 0.3, 0.7), new DM(0)));
              MX lam_f = MX.sym("lam_f");
              Function nlp_hess_l_custom;
              if (scale_helper) {
                MX lag = lam_f.times(opti.f()).plus(ca.dot(opti.lam_g(), og));
                MX H = ca.jacobian(ca.gradient(lag, ox), ox, Map.of("symmetric", true));
                nlp_hess_l_custom = opti.scale_helper(new Function("nlp_hess_l", List.of(ox, op, lam_f, opti.lam_g()), List.of(ca.triu(H)), List.of("x", "p", "lam_f", "lam_g"), List.of("triu_hess_gamma_x_x")));
              } else {
                MX lag = lam_f.times(opti.f()).div(f_scale).plus(ca.dot(opti.lam_g(), og.div(gls)));
                MX H = ca.mtimes(ca.jacobian(ca.gradient(lag, ox), ox, Map.of("symmetric", true)), new MX(ca.power(ca.diag(xls), 2)));
                H = ca.substitute(ca.triu(H), ox, ox.times(xls).plus(xlo));
                nlp_hess_l_custom = new Function("nlp_hess_l", List.of(ox, op, lam_f, opti.lam_g()), List.of(H));
              }
              System.out.println(ref_solver.get_function("nlp_hess_l"));
              checkfunctionLight(ref_solver.get_function("nlp_hess_l"), nlp_hess_l_custom, dms(dm(0.11, 0.3, 0.7), new DM(0), new DM(3), new DM(17)));
            }
          }
        }
      }
    }
  }

  @Test
  public void test_scaling() {
    DM f_ref = null, g_ref = null, x_ref = null, lam_g_ref = null, grad_f_ref = null, jac_g_ref = null, dual_ref = null, lbg_ref = null;
    for (double x_scale : new double[] {1, 100}) {
      for (double g_scale : new double[] {1, 37}) {
        for (double f_scale : new double[] {1, 11}) {
          System.out.println(x_scale + " " + g_scale);
          Opti opti = new Opti();
          MX x = opti.variable();
          MX y = opti.variable();
          opti.variable();
          opti.set_linear_scale(x, x_scale);
          opti.minimize(pow(c(1).minus(x), 2).plus(pow(y.minus(pow(x, 2)), 2).times(100)), f_scale);
          MX my_g = ca.eq(pow(x, 2).plus(pow(y, 2)), 1.0);
          opti.subject_to(my_g, g_scale);
          opti.solver("ipopt", Map.of("ipopt.print_level", 0, "print_time", false));
          OptiSol sol = opti.solve();
          System.out.println(sol.value(opti.lbg()));
          if (x_scale == 1 && g_scale == 1 && f_scale == 1) {
            f_ref = sol.value(opti.f());
            g_ref = sol.value(opti.g());
            x_ref = sol.value(opti.x());
            lam_g_ref = sol.value(opti.lam_g());
            grad_f_ref = sol.value(ca.gradient(opti.f(), opti.x()));
            jac_g_ref = sol.value(ca.jacobian(opti.g(), opti.x()));
            dual_ref = sol.value(opti.dual(my_g));
            lbg_ref = sol.value(opti.lbg());
          } else {
            checkarray(f_ref, sol.value(opti.f()));
            checkarray(g_ref, sol.value(opti.g()));
            checkarray(x_ref, sol.value(opti.x()));
            checkarray(lam_g_ref, sol.value(opti.lam_g()), tol(6), "");
            checkarray(grad_f_ref, sol.value(ca.gradient(opti.f(), opti.x())), tol(6), "");
            checkarray(jac_g_ref, sol.value(ca.jacobian(opti.g(), opti.x())), tol(6), "");
            checkarray(dual_ref, sol.value(opti.dual(my_g)), tol(6), "");
            checkarray(lbg_ref, sol.value(opti.lbg()));
          }
          DM lag_grad = sol.value(ca.gradient(opti.f().plus(ca.dot(opti.lam_g(), opti.g())), opti.x()));
          checkarray(sol.value(x), new DM(0.78641515), tol(6), "");
          checkarray(lag_grad, dm(0, 0), tol(6), "");
        }
      }
    }
  }

  @Test
  public void test_dual() {
    Assumptions.assumeTrue(hasConic("highs"));
    Map<String, DM> ref = new HashMap<>();
    Map<String, DM> data = null;
    for (String t : List.of("nlp", "conic")) {
      Opti opti = new Opti(t);
      MX x = opti.variable();
      MX y = opti.variable();
      MX z = opti.variable();
      opti.minimize(pow(x, 2).plus(pow(y, 2)).plus(pow(z, 2)));
      MX g1 = ca.le(x.times(2).plus(y.times(3)).plus(4), 0.0);
      opti.subject_to(g1);
      MX g2 = ca.eq(x.times(7).plus(z.times(3)).plus(6), 0.0);
      opti.subject_to(g2);
      opti.solver(t.equals("conic") ? "highs" : "ipopt");
      OptiSol sol = opti.solve();
      data = t.equals("conic") ? new HashMap<>() : ref;
      data.put("x", sol.value(opti.x()));
      data.put("lam_g", sol.value(opti.lam_g()));
      data.put("g1", sol.value(opti.dual(g1)));
      data.put("g2", sol.value(opti.dual(g2)));
    }
    System.out.println(data);
    System.out.println(ref);
    for (String k : data.keySet()) checkarray(data.get(k), ref.get(k), tol(6), k);
  }

  @Test
  public void test_dual_soc() {
    Assumptions.assumeTrue(hasConic("superscs"));
    for (String t : List.of("conic")) {
      Opti opti = new Opti(t);
      MX x = opti.variable();
      MX y = opti.variable();
      MX z = opti.variable();
      MX h = ca.soc(vertcat(x.minus(5), y.minus(7)), c(4));
      opti.subject_to(ca.gt(h, 0.0));
      MX g1 = ca.eq(x.plus(y.times(2)).plus(z.times(3)).plus(3), 0.0);
      opti.subject_to(g1);
      opti.minimize(x.times(2).plus(y));
      if (t.equals("conic")) opti.solver("superscs", Map.of("superscs", Map.of("eps", 1e-9, "do_super_scs", 1, "verbose", 0)));
      else opti.solver("ipopt");
      OptiSol sol = opti.solve();
      Map<String, DM> data = new LinkedHashMap<>();
      data.put("f", sol.value(opti.f()));
      data.put("x", sol.value(opti.x()));
      data.put("lam_g", sol.value(opti.lam_g()));
      System.out.println(data);
    }
  }

  @Test
  public void test_value_sx() {
    Opti opti = new Opti();
    checkarray(opti.value(new SX(3)), new DM(3));
    assertInException("Opti works with MX", () -> opti.value(SX.sym("z")));
  }

  @Test
  public void test_subject_to_unknown() {
    Opti opti = new Opti();
    MX x = opti.variable();
    MX y = opti.variable();
    assertInException("Constraint type unknown", () -> opti.subject_to(x.plus(y)));
  }

  @Test
  public void test_optisol_no_debug() {
    assertThrows(NoSuchMethodException.class, () -> OptiSol.class.getMethod("debug"));
  }

  @Test
  public void test_describe_range() {
    Opti opti = new Opti();
    MX x = opti.variable();
    MX y = opti.variable();
    opti.minimize(pow(x, 2).plus(pow(y, 2)));
    assertInException("Index 5 out of range for 2 decision variables", () -> opti.debug().x_describe(5));
    assertInException("Index 0 out of range for 0 constraints", () -> opti.debug().g_describe(0));
  }

  @Test
  public void test_issue4020() {
    List<BiFunction<MX, MX, MX>> constr_funs = List.of(
        (x, y) -> ca.eq(pow(x, 2).minus(y), 0.0),
        (x, y) -> ca.eq(0.0, pow(x, 2).minus(y)),
        (x, y) -> ca.eq(pow(x, 2), y),
        (x, y) -> ca.eq(y, pow(x, 2)));
    for (BiFunction<MX, MX, MX> constr_fun : constr_funs) {
      Opti opti = new Opti();
      MX x = opti.variable();
      MX y = opti.variable();
      opti.minimize(y);
      MX constr = constr_fun.apply(x, y);
      opti.subject_to(constr);
      opti.solver("ipopt");
      OptiSol sol = opti.solve();
      DM dual_opti = sol.value(opti.dual(constr));
      DM dual_all = sol.value(opti.lam_g());
      checkarray(dual_opti, dual_all);
    }
  }
}

package org.casadi.test;

import static org.junit.jupiter.api.Assertions.assertEquals;

import java.util.LinkedHashMap;
import java.util.Map;

import org.casadi.DM;
import org.casadi.Function;
import org.casadi.GenericMatrix;
import org.casadi.SX;
import org.casadi.ca;
import org.junit.jupiter.api.Assumptions;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Disabled;
import org.junit.jupiter.api.Test;

/** Port of test/python/feasiblesqpmethod.py */
public class FeasiblesqpmethodTest extends CasadiTestCase {

  @BeforeEach
  public void requirePlugins() {
    Assumptions.assumeTrue(hasNlpsol("feasiblesqpmethod"), "nlpsol feasiblesqpmethod not available");
    Assumptions.assumeTrue(hasNlpsol("ipopt"), "nlpsol ipopt not available");
  }

  private static Map<String, Object> qpOpts() {
    Map<String, Object> opts = new LinkedHashMap<>();
    opts.put("expand", true);
    opts.put("hessian_approximation", "exact");
    opts.put("qpsol", "nlpsol");
    opts.put("qpsol_options.nlpsol", "ipopt");
    opts.put("qpsol_options.error_on_fail", false);
    opts.put("qpsol_options.nlpsol_options.ipopt.print_level", 0);
    opts.put("qpsol_options.nlpsol_options.ipopt.sb", "yes");
    opts.put("qpsol_options.nlpsol_options.print_time", 0);
    return opts;
  }

  private static long iterCount(Function solver) {
    return ((Number) solver.stats().get("iter_count")).longValue();
  }

  private static DM x(Map<String, GenericMatrix> res) {
    return (DM) res.get("x");
  }

  @Test
  public void test_feas_iter_test() {
    Map<String, Object> opts = qpOpts();
    opts.put("max_inner_iter", 30);
    opts.put("max_iter", 20);
    opts.put("feas_tol", 1e-6);
    opts.put("tr_scale_vector", new double[] {1.0, 1.0});

    SX x = SX.sym("x"), y = SX.sym("y");
    SX f = x.times(x).plus(y.minus(2.0).times(y.minus(2.0)));
    SX g = x.times(x).plus(y.times(y));
    Function solver = ca.nlpsol("solver", "feasiblesqpmethod", Map.of("x", ca.vertcat(java.util.List.of(x, y)), "f", f, "g", g), opts);

    Map<String, Object> arg = new LinkedHashMap<>();
    arg.put("x0", new double[] {2.0, 0.0});
    arg.put("lbg", 4.0);
    arg.put("ubg", 4.0);
    Map<String, GenericMatrix> res = solver.invoke(arg);
    assertEquals(9, iterCount(solver));
    checkarray(x(res), dm(0, 2), 1e-7, "");
  }

  @Test
  public void test_local_conv_test() {
    Map<String, Object> opts = qpOpts();
    SX x = SX.sym("x");
    SX f = ca.log(ca.exp(x).plus(ca.exp(x.unaryMinus())));
    Function solver = ca.nlpsol("solver", "feasiblesqpmethod", Map.of("x", x, "f", f, "g", x), opts);

    System.out.println("Here goes the non-converging initial point for local SQP");
    Map<String, Object> arg = new LinkedHashMap<>();
    arg.put("x0", 1.4);
    arg.put("lbg", -2.0);
    arg.put("ubg", 2.0);
    Map<String, GenericMatrix> res = solver.invoke(arg);
    assertEquals(4, iterCount(solver));
    checkarray(x(res), dm(0), 1e-3, "");
    arg.put("x0", 1.0);
    res = solver.invoke(arg);
    checkarray(x(res), dm(0), 1e-4, "");
    assertEquals(2, iterCount(solver));
  }

  @Test
  @Disabled("parked in Python (class not a casadiTestCase): tr_rad0=1 stops at x=-3.99996, fails digits=6 in Python too")
  public void test_qp_test() {
    double[][] cases = {{1.0, 4}, {100.0, 2}};
    for (double[] c : cases) {
      Map<String, Object> opts = qpOpts();
      opts.put("tr_rad0", c[0]);
      SX x = SX.sym("x");
      SX f = x.plus(4.0).times(x.plus(4.0)).times(2.0).plus(2.0);
      Function solver = ca.nlpsol("solver", "feasiblesqpmethod", Map.of("x", x, "f", f, "g", x), opts);

      Map<String, Object> arg = new LinkedHashMap<>();
      arg.put("x0", 1.0);
      arg.put("lbg", -5.0);
      arg.put("ubg", 5.0);
      Map<String, GenericMatrix> res = solver.invoke(arg);
      assertEquals((long) c[1], iterCount(solver));
      checkarray(x(res), dm(-4), 1e-6, "");
    }
  }
}

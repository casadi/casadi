using System;
using System.Collections.Generic;

namespace CasADi.Tests {

/// <summary>Port of test/java FeasiblesqpmethodTest (itself a port of test/python/feasiblesqpmethod.py)</summary>
public class FeasiblesqpmethodTest : CasadiTestCase {

  // JUnit's @BeforeEach: a SkipTest thrown here skips every test of the class
  public FeasiblesqpmethodTest() {
    Check.Assume(hasNlpsol("feasiblesqpmethod"), "nlpsol feasiblesqpmethod not available");
    Check.Assume(hasNlpsol("ipopt"), "nlpsol ipopt not available");
  }

  static Dictionary<string, object> qpOpts() => new Dictionary<string, object> {
    ["expand"] = true,
    ["hessian_approximation"] = "exact",
    ["qpsol"] = "nlpsol",
    ["qpsol_options.nlpsol"] = "ipopt",
    ["qpsol_options.error_on_fail"] = false,
    ["qpsol_options.nlpsol_options.ipopt.print_level"] = 0,
    ["qpsol_options.nlpsol_options.ipopt.sb"] = "yes",
    ["qpsol_options.nlpsol_options.print_time"] = 0,
  };

  static long iterCount(Function solver) => Convert.ToInt64(solver.stats()["iter_count"]);

  public void test_feas_iter_test() {
    var opts = qpOpts();
    opts["max_inner_iter"] = 30;
    opts["max_iter"] = 20;
    opts["feas_tol"] = 1e-6;
    opts["tr_scale_vector"] = new[] { 1.0, 1.0 };

    SX x = SX.sym("x"), y = SX.sym("y");
    SX f = x * x + (y - 2) * (y - 2);
    SX g = x * x + y * y;
    Function solver = ca.nlpsol("solver", "feasiblesqpmethod", new Dictionary<string, SX> { ["x"] = ca.vertcat(x, y), ["f"] = f, ["g"] = g }, opts);

    var res = solver.call(new Dictionary<string, DM> { ["x0"] = ca.vertcat(2, 0), ["lbg"] = 4, ["ubg"] = 4 });
    assertEquals(9L, iterCount(solver));
    checkarray(res["x"], ca.vertcat(0, 2), 1e-7, "");
  }

  public void test_local_conv_test() {
    var opts = qpOpts();
    SX x = SX.sym("x");
    SX f = ca.log(ca.exp(x) + ca.exp(-x));
    Function solver = ca.nlpsol("solver", "feasiblesqpmethod", new Dictionary<string, SX> { ["x"] = x, ["f"] = f, ["g"] = x }, opts);

    Console.WriteLine("Here goes the non-converging initial point for local SQP");
    var arg = new Dictionary<string, DM> { ["x0"] = 1.4, ["lbg"] = -2, ["ubg"] = 2 };
    var res = solver.call(arg);
    assertEquals(4L, iterCount(solver));
    checkarray(res["x"], ca.vertcat(0), 1e-3, "");
    arg["x0"] = 1.0;
    res = solver.call(arg);
    checkarray(res["x"], ca.vertcat(0), 1e-4, "");
    assertEquals(2L, iterCount(solver));
  }

  // @Disabled in Java: parked in Python (class not a casadiTestCase): tr_rad0=1 stops at x=-3.99996, fails digits=6 in Python too
  public void test_qp_test() {
    Check.Skip("disabled: parked in Python (tr_rad0=1 stops at x=-3.99996, fails digits=6 in Python too)");
    foreach (var (tr_rad0, iters) in new[] { (1.0, 4L), (100.0, 2L) }) {
      var opts = qpOpts();
      opts["tr_rad0"] = tr_rad0;
      SX x = SX.sym("x");
      SX f = (x + 4) * (x + 4) * 2 + 2;
      Function solver = ca.nlpsol("solver", "feasiblesqpmethod", new Dictionary<string, SX> { ["x"] = x, ["f"] = f, ["g"] = x }, opts);

      var res = solver.call(new Dictionary<string, DM> { ["x0"] = 1, ["lbg"] = -5, ["ubg"] = 5 });
      assertEquals(iters, iterCount(solver));
      checkarray(res["x"], ca.vertcat(-4), 1e-6, "");
    }
  }
}

}

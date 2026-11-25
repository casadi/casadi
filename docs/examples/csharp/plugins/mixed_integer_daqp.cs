using System;
using System.Collections.Generic;
using System.Linq;
using CasADi;

public static class mixed_integer_daqp {
  public static void Main() {
    if (Environment.GetEnvironmentVariable("SKIP_DAQP_TESTS") != null || !ca.has_conic("daqp")) {
      Console.WriteLine("SKIP: conic daqp not available");
      return;
    }

    // System.Random replaces numpy's default_rng: same distributions, different numbers
    var rng = new Random(0);
    int n = 100;
    int m = 80;
    int ms = 55;
    SX x = SX.sym("x", n);
    // --- Hessian: dense symmetric matrix ---
    DM M = standard_normal(rng, n, n);
    DM H = 0.5 * ca.mtimes(M, M.T());  // symmetrize

    // --- Linear term ---
    DM f = 100 * standard_normal(rng, n, 1);

    // --- Constraint matrix and bounds ---
    DM A = standard_normal(rng, m, n);
    DM bupper = 20 * random(rng, m, 1);
    DM blower = -20 * random(rng, m, 1);

    // Discrete pattern and bounds: [1]*ms + [0]*(n-ms) etc.
    var discrete = Enumerable.Range(0, n).Select(i => i < ms).ToList();
    double[] lbx = Enumerable.Range(0, n).Select(i => i < ms ? 0.0 : -10).ToArray();
    double[] ubx = Enumerable.Range(0, n).Select(i => i < ms ? 1.0 : 10).ToArray();

    // Using casadi qpsol
    Function solver = ca.qpsol("solver", "daqp",
        new Dictionary<string, SX> { ["f"] = 0.5 * ca.mtimes(ca.mtimes(x.T(), H), x) + ca.mtimes(f.T(), x), ["x"] = x, ["g"] = ca.mtimes(A, x) },
        new Dictionary<string, object> { ["discrete"] = discrete });
    var daqp_sol = solver.call(new Dictionary<string, DM> { ["lbx"] = lbx, ["ubx"] = ubx,
                                                            ["lbg"] = blower, ["ubg"] = bupper });
    Console.WriteLine($"Optimal solution: [{string.Join(", ", daqp_sol["x"].nonzeros())}]");
    Console.WriteLine($"Optimal objective: {(double)daqp_sol["f"]}");

    // Using casadi conic
    // Function solver = conic("solver", "daqp", ...);
    // var daqp_sol = solver.call(...);
    // Console.WriteLine($"Optimal solution: [{string.Join(", ", daqp_sol["x"].nonzeros())}]");
    // Console.WriteLine($"Optimal objective: {(double)daqp_sol["cost"]}");

    Console.WriteLine($"DAQP BnB infos: solver.stats()['bnb_itercount']={solver.stats()["bnb_itercount"]}"
        + $", solver.stats()['bnb_nodecount']={solver.stats()["bnb_nodecount"]}");

    // The remainder cross-checks DAQP against gurobi; skip it alone so that the DAQP
    // part above still runs when gurobi is unavailable.
    if (Environment.GetEnvironmentVariable("SKIP_GUROBI_TESTS") != null || !ca.has_conic("gurobi")) {
      return;
    }

    solver = ca.qpsol("solver", "gurobi",
        new Dictionary<string, SX> { ["f"] = 0.5 * ca.mtimes(ca.mtimes(x.T(), H), x) + ca.mtimes(f.T(), x), ["x"] = x, ["g"] = ca.mtimes(A, x) },
        new Dictionary<string, object> { ["discrete"] = discrete });
    var gu_sol = solver.call(new Dictionary<string, DM> { ["lbx"] = lbx, ["ubx"] = ubx,
                                                          ["lbg"] = blower, ["ubg"] = bupper });
    Console.WriteLine($"Optimal solution: [{string.Join(", ", gu_sol["x"].nonzeros())}]");
    Console.WriteLine($"Optimal objective: {(double)gu_sol["f"]}");

    if (!allclose(daqp_sol["x"].nonzeros(), gu_sol["x"].nonzeros(), 1e-6, 1e-6)) {
      throw new Exception("solutions are different");
    }
    if (!allclose(new[] { (double)gu_sol["f"] }, new[] { (double)daqp_sol["f"] }, 1e-6, 1e-6)) {
      throw new Exception("optimal objectives are different");
    }
  }

  // Box-Muller: a standard normal sample from two uniform ones
  static double gaussian(Random rng) => Math.Sqrt(-2 * Math.Log(1 - rng.NextDouble())) * Math.Cos(2 * Math.PI * rng.NextDouble());

  static DM standard_normal(Random rng, int n, int m) {
    var r = new double[n, m];
    for (int i = 0; i < n; i++) for (int j = 0; j < m; j++) r[i, j] = gaussian(rng);
    return r;
  }

  static DM random(Random rng, int n, int m) {
    var r = new double[n, m];
    for (int i = 0; i < n; i++) for (int j = 0; j < m; j++) r[i, j] = rng.NextDouble();
    return r;
  }

  // numpy.allclose: |a-b| <= atol + rtol*|b|
  static bool allclose(IList<double> a, IList<double> b, double rtol, double atol) {
    for (int i = 0; i < a.Count; i++) if (Math.Abs(a[i] - b[i]) > atol + rtol * Math.Abs(b[i])) return false;
    return true;
  }
}

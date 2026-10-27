import org.casadi.*;
import java.util.Arrays;
import java.util.Map;
import java.util.Random;

public class mixed_integer_daqp {
  public static void main(String[] args) {
    if (System.getenv("SKIP_DAQP_TESTS") != null || !ca.has_conic("daqp")) {
      System.out.println("SKIP: conic daqp not available");
      return;
    }

    // java.util.Random replaces numpy's default_rng: same distributions, different numbers
    Random rng = new Random(0);
    int n = 100;
    int m = 80;
    int ms = 55;
    SX x = SX.sym("x", n);
    // --- Hessian: dense symmetric matrix ---
    DM M = standard_normal(rng, n, n);
    DM H = ca.times(0.5, ca.mtimes(M, M.T()));  // symmetrize

    // --- Linear term ---
    DM f = ca.times(100, standard_normal(rng, n, 1));

    // --- Constraint matrix and bounds ---
    DM A = standard_normal(rng, m, n);
    DM bupper = ca.times(20, random(rng, m, 1));
    DM blower = ca.times(-20, random(rng, m, 1));

    // Discrete pattern and bounds: [1]*ms + [0]*(n-ms) etc.
    boolean[] discrete = new boolean[n];
    double[] lbx = new double[n], ubx = new double[n];
    for (int i = 0; i < n; i++) {
      discrete[i] = i < ms;
      lbx[i] = i < ms ? 0 : -10;
      ubx[i] = i < ms ? 1 : 10;
    }

    // Using casadi qpsol
    Function solver = ca.qpsol("solver", "daqp",
        Map.of("f", ca.times(0.5, ca.mtimes(ca.mtimes(x.T(), H), x)).plus(ca.mtimes(f.T(), x)), "x", x, "g", ca.mtimes(A, x)),
        Map.of("discrete", discrete));
    Map<String, DM> daqp_sol = solver.call(Map.of("lbx", new DM(lbx), "ubx", new DM(ubx),
                                                  "lbg", blower, "ubg", bupper));
    System.out.println("Optimal solution: " + Arrays.toString(daqp_sol.get("x").nonzeros()));
    System.out.println("Optimal objective: " + daqp_sol.get("f").to_double());

    // Using casadi conic
    // Function solver = ca.conic("solver", "daqp", ...);
    // Map<String, DM> daqp_sol = solver.call(...);
    // System.out.println("Optimal solution: " + Arrays.toString(daqp_sol.get("x").nonzeros()));
    // System.out.println("Optimal objective: " + daqp_sol.get("cost").to_double());

    System.out.println("DAQP BnB infos: solver.stats()['bnb_itercount']=" + solver.stats().get("bnb_itercount")
        + ", solver.stats()['bnb_nodecount']=" + solver.stats().get("bnb_nodecount"));

    // The remainder cross-checks DAQP against gurobi; skip it alone so that the DAQP
    // part above still runs when gurobi is unavailable.
    if (System.getenv("SKIP_GUROBI_TESTS") != null || !ca.has_conic("gurobi")) {
      return;
    }

    solver = ca.qpsol("solver", "gurobi",
        Map.of("f", ca.times(0.5, ca.mtimes(ca.mtimes(x.T(), H), x)).plus(ca.mtimes(f.T(), x)), "x", x, "g", ca.mtimes(A, x)),
        Map.of("discrete", discrete));
    Map<String, DM> gu_sol = solver.call(Map.of("lbx", new DM(lbx), "ubx", new DM(ubx),
                                                "lbg", blower, "ubg", bupper));
    System.out.println("Optimal solution: " + Arrays.toString(gu_sol.get("x").nonzeros()));
    System.out.println("Optimal objective: " + gu_sol.get("f").to_double());

    if (!allclose(daqp_sol.get("x").nonzeros(), gu_sol.get("x").nonzeros(), 1e-6, 1e-6)) {
      throw new RuntimeException("solutions are different");
    }
    if (!allclose(new double[] {gu_sol.get("f").to_double()}, new double[] {daqp_sol.get("f").to_double()}, 1e-6, 1e-6)) {
      throw new RuntimeException("optimal objectives are different");
    }
  }

  static DM standard_normal(Random rng, int n, int m) {
    double[][] r = new double[n][m];
    for (double[] row : r) for (int j = 0; j < m; j++) row[j] = rng.nextGaussian();
    return new DM(r);
  }

  static DM random(Random rng, int n, int m) {
    double[][] r = new double[n][m];
    for (double[] row : r) for (int j = 0; j < m; j++) row[j] = rng.nextDouble();
    return new DM(r);
  }

  // numpy.allclose: |a-b| <= atol + rtol*|b|
  static boolean allclose(double[] a, double[] b, double rtol, double atol) {
    for (int i = 0; i < a.length; i++) if (Math.abs(a[i] - b[i]) > atol + rtol * Math.abs(b[i])) return false;
    return true;
  }
}

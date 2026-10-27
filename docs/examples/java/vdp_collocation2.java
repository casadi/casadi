//
//     MIT No Attribution
//
//     Copyright (C) 2010-2023 Joel Andersson, Joris Gillis, Moritz Diehl, KU Leuven.
//
//     Permission is hereby granted, free of charge, to any person obtaining a copy of this
//     software and associated documentation files (the "Software"), to deal in the Software
//     without restriction, including without limitation the rights to use, copy, modify,
//     merge, publish, distribute, sublicense, and/or sell copies of the Software, and to
//     permit persons to whom the Software is furnished to do so.
//
//     THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
//     INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
//     PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
//     HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
//     OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
//     SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
//

import org.casadi.*;
import java.util.*;

public class vdp_collocation2 {
  public static void main(String[] args) {
    if (!ca.has_nlpsol("ipopt")) { System.out.println("SKIP: nlpsol ipopt not available"); return; }
    if (!hasMa27()) { System.out.println("SKIP: ipopt linear solver ma27 (libhsl) not available"); return; }
    double inf = Double.POSITIVE_INFINITY;

    int nk = 20;      // Control discretization
    double tf = 10.0; // End time

    // Declare variables (use scalar graph)
    SX t = SX.sym("t"); // time
    SX u = SX.sym("u"); // control

    // states = [x (vdp oscillator states, 2); L (helper state: Langrange integrand)]
    SX states = SX.sym("states", 3);

    // Create a structure for the right hand side
    SX x = states.get(new Slice(0, 2));
    SX rhs = ca.vertcat(List.of(
        ca.vertcat(List.of(x.get(1).times(x.get(1)).unaryMinus().plus(1).times(x.get(0)).minus(x.get(1)).plus(u), x.get(0))),
        x.get(0).times(x.get(0)).plus(x.get(1).times(x.get(1))).plus(u.times(u))));

    // ODE right hand side function
    Function f = new Function("f", List.of(t, states, u), List.of(rhs));

    // Objective function (meyer term)
    Function m = new Function("m", List.of(t, states, u), List.of(states.get(2)));

    // Control bounds
    double u_min = -0.75;
    double u_max = 1.0;
    double u_init = 0.0;

    // State bounds and initial guess
    double[] x_min = {-inf, -inf, -inf};
    double[] x_max = { inf,  inf,  inf};
    double[] xi_min = {0.0, 1.0, 0.0};
    double[] xi_max = {0.0, 1.0, 0.0};
    double[] xf_min = {0.0, 0.0, -inf};
    double[] xf_max = {0.0, 0.0,  inf};
    double[] x_init = {0.0, 0.0, 0.0};

    // Dimensions
    int nx = 3;
    int nu = 1;

    // Choose collocation points
    double[] cp = ca.collocation_points(3, "radau");
    double[] tau_root = new double[cp.length + 1];
    System.arraycopy(cp, 0, tau_root, 1, cp.length);

    // Degree of interpolating polynomial
    int d = tau_root.length - 1;

    // Size of the finite elements
    double h = tf / nk;

    // Coefficients of the collocation equation
    double[][] C = new double[d + 1][d + 1];

    // Coefficients of the continuity equation
    double[] D = new double[d + 1];

    // Dimensionless time inside one control interval
    SX tau = SX.sym("tau");

    // All collocation time points
    double[][] T = new double[nk][d + 1];
    for (int k = 0; k < nk; k++)
      for (int j = 0; j < d + 1; j++)
        T[k][j] = h * (k + tau_root[j]);

    // For all collocation points
    for (int j = 0; j < d + 1; j++) {
      // Construct Lagrange polynomials to get the polynomial basis at the collocation point
      SX L = new SX(1);
      for (int r = 0; r < d + 1; r++) {
        if (r != j) {
          L = L.times(tau.minus(tau_root[r]).div(tau_root[j] - tau_root[r]));
        }
      }

      // Evaluate the polynomial at the final time to get the coefficients of the continuity equation
      Function lfcn = new Function("lfcn", List.of(tau), List.of(L));
      D[j] = ((DM) lfcn.invoke(1.0)).to_double();

      // Evaluate the time derivative of the polynomial at all collocation points to get the coefficients of the continuity equation
      Function tfcn = new Function("tfcn", List.of(tau), List.of(ca.tangent(L, tau)));
      for (int r = 0; r < d + 1; r++) {
        C[j][r] = ((DM) tfcn.invoke(tau_root[r])).to_double();
      }
    }

    // Structure holding NLP variables: V = [X[0][0..d], U[0], X[1][0..d], U[1], ..., X[nk][0..d]]
    int blk = (d + 1) * nx + nu;
    int NV = nk * blk + (d + 1) * nx;
    MX V = MX.sym("V", NV);
    int[][] iX = new int[nk + 1][d + 1];
    int[] iU = new int[nk];
    for (int k = 0; k < nk + 1; k++)
      for (int j = 0; j < d + 1; j++) iX[k][j] = k * blk + j * nx;
    for (int k = 0; k < nk; k++) iU[k] = k * blk + (d + 1) * nx;

    double[] vars_lb = new double[NV];
    double[] vars_ub = new double[NV];
    double[] vars_init = new double[NV];

    // Set states and its bounds
    for (int k = 0; k < nk + 1; k++)
      for (int j = 0; j < d + 1; j++)
        for (int i = 0; i < nx; i++) {
          vars_init[iX[k][j] + i] = x_init[i];
          vars_lb[iX[k][j] + i] = x_min[i];
          vars_ub[iX[k][j] + i] = x_max[i];
        }

    // Set controls and its bounds
    for (int k = 0; k < nk; k++) {
      vars_init[iU[k]] = u_init;
      vars_lb[iU[k]] = u_min;
      vars_ub[iU[k]] = u_max;
    }

    // State at initial time
    for (int i = 0; i < nx; i++) {
      vars_lb[iX[0][0] + i] = xi_min[i];
      vars_ub[iX[0][0] + i] = xi_max[i];
    }

    // State at end time
    for (int i = 0; i < nx; i++) {
      vars_lb[iX[nk][0] + i] = xf_min[i];
      vars_ub[iX[nk][0] + i] = xf_max[i];
    }

    // Constraint function for the NLP
    List<MX> g = new ArrayList<>();
    List<Double> lbg = new ArrayList<>();
    List<Double> ubg = new ArrayList<>();

    // For all finite elements
    for (int k = 0; k < nk; k++) {

      // For all collocation points
      for (int j = 1; j < d + 1; j++) {

        // Get an expression for the state derivative at the collocation point
        MX xp_jk = new MX(0);
        for (int r = 0; r < d + 1; r++) {
          xp_jk = xp_jk.plus(VX(V, iX, nx, k, r).times(C[r][j]));
        }

        // Add collocation equations to the NLP
        MX fk = (MX) f.invoke(T[k][j], VX(V, iX, nx, k, j), V.get(iU[k]));
        g.add(fk.times(h).minus(xp_jk));
        for (int i = 0; i < nx; i++) { lbg.add(0.0); ubg.add(0.0); } // equality constraints
      }

      // Get an expression for the state at the end of the finite element
      MX xf_k = new MX(0);
      for (int r = 0; r < d + 1; r++) {
        xf_k = xf_k.plus(VX(V, iX, nx, k, r).times(D[r]));
      }

      // Add continuity equation to NLP
      g.add(VX(V, iX, nx, k + 1, 0).minus(xf_k));
      for (int i = 0; i < nx; i++) { lbg.add(0.0); ubg.add(0.0); }
    }

    // Concatenate constraints
    MX gc = ca.vertcat(g);

    // Objective function
    MX fobj = (MX) m.invoke(T[nk - 1][d], VX(V, iX, nx, nk, 0), V.get(iU[nk - 1]));

    // NLP
    Map<String, MX> nlp = Map.of("x", V, "f", fobj, "g", gc);

    // ----
    // SOLVE THE NLP
    // ----

    // Set options
    Map<String, Object> opts = new LinkedHashMap<>();
    opts.put("expand", true);
    // opts.put("ipopt.max_iter", 4);
    opts.put("ipopt.linear_solver", "ma27");

    // Allocate an NLP solver
    Function solver = ca.nlpsol("solver", "ipopt", nlp, opts);
    Map<String, Object> arg = new LinkedHashMap<>();

    // Initial condition
    arg.put("x0", vars_init);

    // Bounds on x
    arg.put("lbx", vars_lb);
    arg.put("ubx", vars_ub);

    // Bounds on g
    arg.put("lbg", lbg.stream().mapToDouble(Double::doubleValue).toArray());
    arg.put("ubg", ubg.stream().mapToDouble(Double::doubleValue).toArray());

    // Solve the problem
    Map<String, GenericMatrix> res = solver.invoke(arg);

    // Print the optimal cost
    System.out.println("optimal cost:  " + ((DM) res.get("f")).to_double());

    // Retrieve the solution
    double[] opt = ((DM) res.get("x")).nonzeros();

    // Get values at the beginning of each finite element
    double[] x0_opt = new double[nk + 1], x1_opt = new double[nk + 1], x2_opt = new double[nk + 1];
    for (int k = 0; k < nk + 1; k++) {
      x0_opt[k] = opt[iX[k][0]];
      x1_opt[k] = opt[iX[k][0] + 1];
      x2_opt[k] = opt[iX[k][0] + 2];
    }

    double[] u_opt = new double[nk];
    for (int k = 0; k < nk; k++) u_opt[k] = opt[iU[k]];

    DM tgrid = ca.linspace(new DM(0), new DM(tf), nk + 1);
    DM tgrid_u = ca.linspace(new DM(0), new DM(tf), nk);

    // plot: x0_opt, x1_opt over tgrid and u_opt (step) over tgrid_u, "Van der Pol optimization"
  }

  static boolean hasMa27() {
    SX z = SX.sym("z");
    Function s = ca.nlpsol("probe", "ipopt", Map.of("x", z, "f", z.times(z)),
        Map.of("ipopt.linear_solver", "ma27", "ipopt.print_level", 0, "print_time", false));
    s.invoke(Map.of("x0", 1));
    return (Boolean) s.stats().get("success");
  }

  // V["X",k,j]
  static MX VX(MX V, int[][] iX, int nx, int k, int j) {
    return V.get(new Slice(iX[k][j], iX[k][j] + nx));
  }
}

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
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import org.casadi.DM;
import org.casadi.Function;
import org.casadi.MX;
import org.casadi.SX;
import org.casadi.Slice;
import org.casadi.ca;

public class vdp_indirect_multiple_shooting {
  public static void main(String[] args) {
    if (!ca.has_integrator("cvodes")) { System.out.println("SKIP: integrator cvodes not available"); return; }
    if (!ca.has_rootfinder("nlpsol") || !ca.has_nlpsol("ipopt")) {
      System.out.println("SKIP: rootfinder nlpsol with ipopt not available");
      return;
    }

    // Declare variables (use simple, efficient DAG)
    SX x0 = SX.sym("x0"); SX x1 = SX.sym("x1");
    SX x = ca.vertcat(List.of(x0, x1));

    // Control
    SX u = SX.sym("u");

    // ODE right hand side
    SX xdot = ca.vertcat(List.of(ca.minus(1, x1.times(x1)).times(x0).minus(x1).plus(u), x0));

    // Lagrangian function
    SX L = x0.times(x0).plus(x1.times(x1)).plus(u.times(u));

    // Costate
    SX lam = SX.sym("lam", 2);

    // Hamiltonian function
    SX H = ca.dot(lam, xdot).plus(L);

    // Costate equations
    SX ldot = ca.gradient(H, x).unaryMinus();

    //# The control must minimize the Hamiltonian, which is:
    System.out.println("Hamiltonian:  " + H);

    // H is of a convex quadratic form in u: H = u*u + p*u + q, let's get the coefficient p
    SX p = ca.gradient(H, u);     // this gives us 2*u + p
    p = ca.substitute(p, u, new SX(0)); // replace u with zero: gives us p

    // H's unconstrained minimizer is: u = -p/2
    SX u_opt = p.unaryMinus().div(2);

    // We must constrain u to the interval [-0.75, 1.0], convexity of H ensures that the optimum is obtain at the bound when u_opt is outside the interval
    u_opt = ca.fmin(u_opt, 1.0);
    u_opt = ca.fmax(u_opt, -0.75);
    System.out.println("optimal control:  " + u_opt);

    // Augment f with lam_dot and subtitute in the value for the optimal control
    SX f = ca.vertcat(List.of(xdot, ldot));
    f = ca.substitute(f, u, u_opt);

    // Function for obtaining the optimal control from the augmented state
    Function u_fcn = new Function("ufcn", List.of(ca.vertcat(List.of(x, lam))), List.of(u_opt));

    // Formulate the DAE
    Map<String, SX> dae = Map.of("x", ca.vertcat(List.of(x, lam)), "ode", f);

    // Augmented DAE state dimension
    int nX = 4;

    // End time
    double tf = 10.0;

    // Number of shooting nodes
    int num_nodes = 20;

    // Create an integrator (CVodes)
    Map<String, Object> iopts = new LinkedHashMap<>();
    iopts.put("abstol", 1e-8); // abs. tolerance
    iopts.put("reltol", 1e-8); // rel. tolerance
    Function I = ca.integrator("I", "cvodes", dae, 0, tf / num_nodes, iopts);

    // Variables for the states at each shooting node
    MX X = MX.sym("X", nX, num_nodes + 1);

    // Formulate the root finding problem
    List<MX> G = new ArrayList<>();
    G.add(X.get(new Slice(0, 2), 0).minus(new DM(new double[] {0, 1}))); // states fixed, costates free at initial time
    for (int k = 0; k < num_nodes; k++) {
      MX XF = I.call(Map.of("x0", X.get(new Slice(), k))).get("xf");
      G.add(XF.minus(X.get(new Slice(), k + 1)));
    }
    G.add(X.get(new Slice(2, nX), num_nodes).minus(new DM(new double[] {0, 0}))); // costates fixed, states free at final time

    // Terminal constraints: lam = 0
    Map<String, MX> rfp = Map.of("x", ca.vec(X), "g", ca.vertcat(G));

    // Select a solver for the root-finding problem
    String Solver = "nlpsol";
    //String Solver = "newton";
    //String Solver = "kinsol";

    // Solver options
    Map<String, Object> opts = new LinkedHashMap<>();
    if (Solver.equals("nlpsol")) {
      opts.put("nlpsol", "ipopt");
      opts.put("nlpsol_options", Map.of("ipopt.hessian_approximation", "limited-memory"));
    } else if (Solver.equals("newton")) {
      opts.put("linear_solver", "csparse");
    } else if (Solver.equals("kinsol")) {
      opts.put("linear_solver_type", "user_defined");
      opts.put("linear_solver", "csparse");
      opts.put("max_iter", 1000);
    }

    // Allocate a solver
    Function solver = ca.rootfinder("solver", Solver, rfp, opts);

    // Solve the problem
    DM X_sol = (DM) solver.invoke(Map.of("x0", 0)).get("x");
    System.out.println(X_sol);

    // Time grid for visualization
    double[] tgrid = ca.linspace(0, tf, 100).nonzeros();

    // Simulator to get optimal state and control trajectories
    Function simulator = ca.integrator("simulator", "cvodes", dae, 0, tgrid);

    // Simulate to get the trajectories
    DM sol = (DM) simulator.invoke(Map.of("x0", X_sol.get(new Slice(0, 4)))).get("xf");

    // Calculate the optimal control
    DM u_opt_traj = (DM) u_fcn.invoke(sol);

    // plot: sol[0,:], sol[1,:] and u_opt versus tgrid, 'Van der Pol optimization - indirect multiple shooting'

    System.out.println(sol.get(0, new Slice()));
  }
}

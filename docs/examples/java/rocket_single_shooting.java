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
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;

import org.casadi.DM;
import org.casadi.Function;
import org.casadi.MX;
import org.casadi.SX;
import org.casadi.ca;

public class rocket_single_shooting {
  static boolean sundials_integrator = true;
  static boolean explicit_integrator = false;
  static boolean lifted_newton = false;

  public static void main(String[] args) {
    if (!ca.has_integrator("cvodes") || !ca.has_integrator("idas") || !ca.has_integrator("rk")
        || !ca.has_nlpsol("ipopt") || !ca.has_nlpsol("scpgen")) {
      System.out.println("SKIP: integrator cvodes/idas/rk or nlpsol ipopt/scpgen not available");
      return;
    }

    // Load integrators manually
    ca.load_integrator("cvodes");
    ca.load_integrator("idas");
    ca.load_integrator("rk");
    ca.load_nlpsol("ipopt");
    ca.load_nlpsol("scpgen");

    // Time length
    double T = 10.0;

    // Shooting length
    int nu = 20; // Number of control segments

    // Time horizon for integrator
    double t0 = 0;
    double tf = T / nu;

    // Initial position
    double[] X0 = new double[3];
    X0[0] = 0; // initial position
    X0[1] = 0; // initial speed
    X0[2] = 1; // initial mass

    // Time
    SX t = SX.sym("t");

    // Differential states
    SX s = SX.sym("s"), v = SX.sym("v"), m = SX.sym("m");
    SX x = ca.vertcat(s, v, m);

    // Control
    SX u = SX.sym("u");

    SX alpha = new SX(0.05); // friction
    SX beta = new SX(0.1); // fuel consumption rate

    // Differential equation
    SX rhs = ca.vertcat(v, u.minus(alpha.times(v).times(v)).div(m), beta.unaryMinus().times(u).times(u));

    // Initial conditions
    double[] x0 = {0, 0, 1};

    // DAE
    Map<String, SX> dae = Map.of("x", x, "p", u, "t", t, "ode", rhs);

    // Integrator options
    String plugin;
    Map<String, Object> opts = new HashMap<>();
    if (sundials_integrator) {
      if (explicit_integrator) {
        // Explicit integrator (CVODES)
        plugin = "cvodes";
        // opts.put("exact_jacobian", true);
        // opts.put("linear_multistep_method", "bdf"); // adams or bdf
        // opts.put("nonlinear_solver_iteration", "newton"); // newton or functional
      } else {
        // Implicit integrator (IDAS)
        plugin = "idas";
        opts.put("calc_ic", false);
      }
      opts.put("fsens_err_con", true);
      opts.put("quad_err_con", true);
      opts.put("abstol", 1e-6);
      opts.put("reltol", 1e-6);
      //  opts.put("fsens_all_at_once", false);
      opts.put("steps_per_checkpoint", 100); // BUG: Too low number causes segfaults
    } else {
      // An explicit Euler integrator
      plugin = "rk";
      opts.put("expand_f", true);
      opts.put("interpolation_order", 1);
      opts.put("number_of_finite_elements", 1000);
    }

    // Create integrator
    Function F = ca.integrator("integrator", plugin, dae, t0, tf, opts);

    // control for all segments
    MX U = MX.sym("U", nu);

    // Integrate over all intervals
    MX X = new MX(X0);
    for (int k = 0; k < nu; ++k) {
      // Integrate
      X = F.call(Map.of("x0", X, "p", U.get(k))).get("xf");

      // Lift X
      if (lifted_newton) {
        // BINDING: MX lift(x, x_guess) is not wrapped; it should read
        // X = ca.lift(X, X);
        System.out.println("SKIP: BINDING: ca.lift(MX, MX) is not wrapped");
        return;
      }
    }

    // Objective function
    MX J = ca.dot(U, U);

    // Terminal constraints
    MX G = ca.vertcat(X.get(0), X.get(1));

    // Create the NLP
    Map<String, MX> nlp = Map.of("x", U, "f", J, "g", G);

    // NLP solver options
    Map<String, Object> solver_opts = new HashMap<>();
    String solver_name;
    if (lifted_newton) {
      solver_name = "scpgen";
      solver_opts.put("verbose", true);
      solver_opts.put("regularize", false);
      solver_opts.put("max_iter_ls", 1);
      solver_opts.put("max_iter", 100);
      solver_opts.put("qpsol", "nlp"); // Use IPOPT as QP solver
      Map<String, Object> ipopt_options = new HashMap<>();
      ipopt_options.put("tol", 1e-12);
      ipopt_options.put("print_level", 0);
      ipopt_options.put("print_time", false);
      solver_opts.put("qpsol_options",
        Map.of("nlpsol_options", Map.of("nlpsol", "ipopt")));
    } else {
      solver_name = "ipopt";
      solver_opts.put("ipopt.tol", 1e-10);
      solver_opts.put("ipopt.hessian_approximation", "limited-memory");
    }

    // Create NLP solver
    Function solver = ca.nlpsol("nlpsol", solver_name, nlp, solver_opts);

    // Bounds on u and initial condition
    double[] umin = new double[nu], umax = new double[nu], u0 = new double[nu];
    Arrays.fill(umin, -10);
    Arrays.fill(umax, 10);
    Arrays.fill(u0, 0.4);

    // Bounds on g
    double[] gmin = {10, 0}, gmax = {10, 0};

    // Solve NLP
    double[] Usol = solver.call(Map.of("lbx", new DM(umin), "ubx", new DM(umax), "x0", new DM(u0),
                                       "lbg", new DM(gmin), "ubg", new DM(gmax))).get("x").nonzeros();
    System.out.println("optimal solution: " + str(Usol));
  }

  // C++ str(): 6 significant digits, as std::ostream
  static String str(double x) { return new DM(x).toString(); }
  static String str(double[] x) {
    StringBuilder r = new StringBuilder("[");
    for (int i = 0; i < x.length; ++i) r.append(i > 0 ? ", " : "").append(str(x[i]));
    return r.append("]").toString();
  }
}

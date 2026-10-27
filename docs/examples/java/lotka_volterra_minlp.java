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
import java.util.Collections;
import java.util.List;
import java.util.Map;

import org.casadi.DM;
import org.casadi.Function;
import org.casadi.GenericMatrix;
import org.casadi.MX;
import org.casadi.SX;
import org.casadi.ca;

// Lotka-Volerra problem from mintoc.de
public class lotka_volterra_minlp {
  public static void main(String[] args) {
    if (!ca.has_nlpsol("bonmin")) { System.out.println("SKIP: nlpsol bonmin not available"); return; }

    double T = 12.; // Time horizon
    int N = 40; // number of control intervals

    // Model parameters
    double c0 = .4;
    double c1 = .2;

    // Initial condition for x
    List<Double> x0 = List.of(0.5, 0.7);

    // Bounds on x
    List<Double> lbx = List.of(0., 0.);
    List<Double> ubx = List.of(2., 2.);

    // Bounds on u
    double lbu = 0.;
    double ubu = 1.;

    // Declare model variables
    SX x = SX.sym("x", 2);
    SX u = SX.sym("u");

    // Model equations
    SX xdot = ca.vertcat(List.of(x.get(0).minus(x.get(0).times(x.get(1))).minus(ca.times(c0, x.get(0)).times(u)),
        x.get(1).unaryMinus().plus(x.get(0).times(x.get(1))).minus(ca.times(c1, x.get(1)).times(u))));

    // Objective term
    SX L = x.get(0).minus(1).power(2).plus(x.get(1).minus(1).power(2)).plus(ca.times(1e-4, u.power(2)));

    // Formulate discrete time dynamics
    Function F;
    if (false) {
      // CVODES from the SUNDIALS suite
      Map<String, SX> dae = Map.of("x", x, "p", u, "ode", xdot, "quad", L);
      F = ca.integrator("F", "cvodes", dae, 0, T / N);
    } else {
      // Fixed step Runge-Kutta 4 integrator
      int M = 4; // RK4 steps per interval
      double DT = T / N / M;
      Function f = new Function("f", List.of(x, u), List.of(xdot, L));
      MX X0 = MX.sym("X0", 2);
      MX U = MX.sym("U");
      MX X = X0;
      MX Q = new MX(0);
      for (int j = 0; j < M; j++) {
        List<MX> k1 = f.call(List.of(X, U));
        List<MX> k2 = f.call(List.of(X.plus(k1.get(0).times(DT / 2)), U));
        List<MX> k3 = f.call(List.of(X.plus(k2.get(0).times(DT / 2)), U));
        List<MX> k4 = f.call(List.of(X.plus(k3.get(0).times(DT)), U));
        X = X.plus(k1.get(0).plus(k2.get(0).times(2)).plus(k3.get(0).times(2)).plus(k4.get(0)).times(DT / 6));
        Q = Q.plus(k1.get(1).plus(k2.get(1).times(2)).plus(k3.get(1).times(2)).plus(k4.get(1)).times(DT / 6));
      }
      F = new Function("F", List.of(X0, U), List.of(X, Q), List.of("x0", "p"), List.of("xf", "qf"));
    }

    // Initial guess for u
    List<DM> u_start = Collections.nCopies(N, new DM(0.));

    // Get a feasible trajectory as an initial guess
    DM xk = new DM(toArray(x0));
    List<DM> x_start = new ArrayList<>(List.of(xk));
    for (int k = 0; k < N; k++) {
      xk = (DM) F.invoke(Map.of("x0", xk, "p", u_start.get(k))).get("xf");
      x_start.add(xk);
    }

    // Start with an empty NLP
    List<MX> w = new ArrayList<>();
    List<DM> w0 = new ArrayList<>();
    List<Double> lbw = new ArrayList<>();
    List<Double> ubw = new ArrayList<>();
    List<Boolean> discrete = new ArrayList<>();
    MX J = new MX(0);
    List<MX> g = new ArrayList<>();
    List<Double> lbg = new ArrayList<>();
    List<Double> ubg = new ArrayList<>();

    // "Lift" initial conditions
    MX X0 = MX.sym("X0", 2);
    w.add(X0);
    lbw.addAll(x0);
    ubw.addAll(x0);
    w0.add(x_start.get(0));
    discrete.addAll(List.of(false, false));

    // Formulate the NLP
    MX Xk = X0;
    for (int k = 0; k < N; k++) {
      // New NLP variable for the control
      MX Uk = MX.sym("U_" + k);
      w.add(Uk);
      lbw.add(lbu);
      ubw.add(ubu);
      w0.add(u_start.get(k));
      discrete.add(true);

      // Integrate till the end of the interval
      Map<String, MX> Fk = F.call(Map.of("x0", Xk, "p", Uk));
      MX Xk_end = Fk.get("xf");
      J = J.plus(Fk.get("qf"));

      // New NLP variable for state at end of interval
      Xk = MX.sym("X_" + (k + 1), 2);
      w.add(Xk);
      lbw.addAll(lbx);
      ubw.addAll(ubx);
      w0.add(x_start.get(k + 1));
      discrete.addAll(List.of(false, false));

      // Add equality constraint
      g.add(Xk_end.minus(Xk));
      lbg.addAll(List.of(0.0, 0.0));
      ubg.addAll(List.of(0.0, 0.0));
    }

    // Concatenate decision variables and constraint terms
    MX w_all = ca.vertcat(w);
    MX g_all = ca.vertcat(g);

    // Create an NLP solver
    Map<String, MX> nlp_prob = Map.of("f", J, "x", w_all, "g", g_all);
    boolean[] discrete_arr = new boolean[discrete.size()];
    for (int i = 0; i < discrete_arr.length; i++) discrete_arr[i] = discrete.get(i);
    Function nlp_solver = ca.nlpsol("nlp_solver", "bonmin", nlp_prob, Map.of("discrete", discrete_arr));
    //Function nlp_solver = ca.nlpsol("nlp_solver", "knitro", nlp_prob, Map.of("discrete", discrete_arr));
    //Function nlp_solver = ca.nlpsol("nlp_solver", "ipopt", nlp_prob); // Solve relaxed problem

    // Plot the solution
    double[] tgrid = new double[N + 1];
    for (int k = 0; k <= N; k++) tgrid[k] = T / N * k;

    // Solve the NLP
    Map<String, GenericMatrix> sol = nlp_solver.invoke(Map.of("x0", ca.vertcat(w0), "lbx", toArray(lbw),
        "ubx", toArray(ubw), "lbg", toArray(lbg), "ubg", toArray(ubg)));

    System.out.println(nlp_solver.stats());

    DM w1_opt = (DM) sol.get("x");
    DM lam_w_opt = (DM) sol.get("lam_x");
    DM lam_g_opt = (DM) sol.get("lam_g");
    // plot: w1_opt[0::3], w1_opt[1::3] and (step) w1_opt[2::3] versus tgrid
  }

  static double[] toArray(List<Double> v) {
    return v.stream().mapToDouble(Double::doubleValue).toArray();
  }
}

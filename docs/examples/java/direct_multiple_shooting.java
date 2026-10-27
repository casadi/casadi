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
import java.util.List;
import java.util.Map;

import org.casadi.DM;
import org.casadi.Function;
import org.casadi.GenericMatrix;
import org.casadi.MX;
import org.casadi.ca;

public class direct_multiple_shooting {
  public static void main(String[] args) {
    if (!ca.has_nlpsol("ipopt")) { System.out.println("SKIP: nlpsol ipopt not available"); return; }
    double inf = Double.POSITIVE_INFINITY;

    double T = 10.; // Time horizon
    int N = 20; // number of control intervals

    // Declare model variables
    MX x1 = MX.sym("x1");
    MX x2 = MX.sym("x2");
    MX x = ca.vertcat(List.of(x1, x2));
    MX u = MX.sym("u");

    // Model equations
    MX xdot = ca.vertcat(List.of(ca.minus(1, x2.power(2)).times(x1).minus(x2).plus(u), x1));

    // Objective term
    MX L = x1.power(2).plus(x2.power(2)).plus(u.power(2));

    // Formulate discrete time dynamics
    Function F;
    if (false) {
      // CVODES from the SUNDIALS suite
      Map<String, MX> dae = Map.of("x", x, "p", u, "ode", xdot, "quad", L);
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

    // Evaluate at a test point
    Map<String, GenericMatrix> Fk = F.invoke(Map.of("x0", new double[] {0.2, 0.3}, "p", 0.4));
    System.out.println(Fk.get("xf"));
    System.out.println(Fk.get("qf"));

    // Start with an empty NLP
    List<MX> w = new ArrayList<>();
    List<Double> w0 = new ArrayList<>();
    List<Double> lbw = new ArrayList<>();
    List<Double> ubw = new ArrayList<>();
    MX J = new MX(0);
    List<MX> g = new ArrayList<>();
    List<Double> lbg = new ArrayList<>();
    List<Double> ubg = new ArrayList<>();

    // "Lift" initial conditions
    MX Xk = MX.sym("X0", 2);
    w.add(Xk);
    lbw.addAll(List.of(0.0, 1.0));
    ubw.addAll(List.of(0.0, 1.0));
    w0.addAll(List.of(0.0, 1.0));

    // Formulate the NLP
    for (int k = 0; k < N; k++) {
      // New NLP variable for the control
      MX Uk = MX.sym("U_" + k);
      w.add(Uk);
      lbw.add(-1.0);
      ubw.add(1.0);
      w0.add(0.0);

      // Integrate till the end of the interval
      Map<String, MX> Fk_sym = F.call(Map.of("x0", Xk, "p", Uk));
      MX Xk_end = Fk_sym.get("xf");
      J = J.plus(Fk_sym.get("qf"));

      // New NLP variable for state at end of interval
      Xk = MX.sym("X_" + (k + 1), 2);
      w.add(Xk);
      lbw.addAll(List.of(-0.25, -inf));
      ubw.addAll(List.of(inf, inf));
      w0.addAll(List.of(0.0, 0.0));

      // Add equality constraint
      g.add(Xk_end.minus(Xk));
      lbg.addAll(List.of(0.0, 0.0));
      ubg.addAll(List.of(0.0, 0.0));
    }

    // Create an NLP solver
    Map<String, MX> prob = Map.of("f", J, "x", ca.vertcat(w), "g", ca.vertcat(g));
    Function solver = ca.nlpsol("solver", "ipopt", prob);

    // Solve the NLP
    Map<String, GenericMatrix> sol = solver.invoke(Map.of("x0", toArray(w0), "lbx", toArray(lbw), "ubx", toArray(ubw),
        "lbg", toArray(lbg), "ubg", toArray(ubg)));
    double[] w_opt = ((DM) sol.get("x")).nonzeros();

    // Plot the solution
    double[] x1_opt = every(w_opt, 0, 3);
    double[] x2_opt = every(w_opt, 1, 3);
    double[] u_opt = every(w_opt, 2, 3);

    double[] tgrid = new double[N + 1];
    for (int k = 0; k <= N; k++) tgrid[k] = T / N * k;
    // plot: x1, x2 and (step) u versus tgrid
  }

  static double[] toArray(List<Double> v) {
    return v.stream().mapToDouble(Double::doubleValue).toArray();
  }

  // Python w[start::step]
  static double[] every(double[] w, int start, int step) {
    double[] r = new double[(w.length - start + step - 1) / step];
    for (int i = 0; i < r.length; i++) r[i] = w[start + i * step];
    return r;
  }
}

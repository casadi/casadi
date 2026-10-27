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
import java.util.ArrayList;
import java.util.List;
import java.util.Map;

public class direct_single_shooting {
  @SuppressWarnings("unchecked")
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
        List<MX> k2 = f.call(List.of(X.plus(ca.times(DT / 2, k1.get(0))), U));
        List<MX> k3 = f.call(List.of(X.plus(ca.times(DT / 2, k2.get(0))), U));
        List<MX> k4 = f.call(List.of(X.plus(ca.times(DT, k3.get(0))), U));
        X = X.plus(ca.times(DT / 6, k1.get(0).plus(k2.get(0).times(2)).plus(k3.get(0).times(2)).plus(k4.get(0))));
        Q = Q.plus(ca.times(DT / 6, k1.get(1).plus(k2.get(1).times(2)).plus(k3.get(1).times(2)).plus(k4.get(1))));
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

    // Formulate the NLP
    MX Xk = new MX(new DM(new double[] {0, 1}));
    for (int k = 0; k < N; k++) {
      // New NLP variable for the control
      MX Uk = MX.sym("U_" + k);
      w.add(Uk);
      lbw.add(-1.0);
      ubw.add(1.0);
      w0.add(0.0);

      // Integrate till the end of the interval
      Map<String, MX> Fkk = F.call(Map.of("x0", Xk, "p", Uk));
      Xk = Fkk.get("xf");
      J = J.plus(Fkk.get("qf"));

      // Add inequality constraint
      g.add(Xk.get(0));
      lbg.add(-.25);
      ubg.add(inf);
    }

    // Create an NLP solver
    Map<String, MX> prob = Map.of("f", J, "x", ca.vertcat(w), "g", ca.vertcat(g));
    Function solver = ca.nlpsol("solver", "ipopt", prob);

    // Solve the NLP
    Map<String, GenericMatrix> sol = solver.invoke(Map.of("x0", arr(w0), "lbx", arr(lbw), "ubx", arr(ubw),
                                                          "lbg", arr(lbg), "ubg", arr(ubg)));
    DM w_opt = (DM) sol.get("x");

    // Plot the solution
    DM u_opt = w_opt;
    List<DM> x_opt = new ArrayList<>(List.of(new DM(new double[] {0, 1})));
    for (int k = 0; k < N; k++) {
      Fk = F.invoke(Map.of("x0", x_opt.get(x_opt.size() - 1), "p", u_opt.get(k)));
      x_opt.add((DM) Fk.get("xf"));
    }
    List<DM> x1s = new ArrayList<>(), x2s = new ArrayList<>();
    for (DM r : x_opt) { x1s.add(r.get(0)); x2s.add(r.get(1)); }
    DM x1_opt = ca.vertcat(x1s);
    DM x2_opt = ca.vertcat(x2s);

    double[] tgrid = new double[N + 1];
    for (int k = 0; k <= N; k++) tgrid[k] = T / N * k;
    // plot: x1_opt (dashed), x2_opt over tgrid, and step plot of vertcat(nan, u_opt)
  }

  static double[] arr(List<Double> v) {
    return v.stream().mapToDouble(Double::doubleValue).toArray();
  }
}

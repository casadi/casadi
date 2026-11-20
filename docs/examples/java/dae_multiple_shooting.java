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
//
import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import org.casadi.DM;
import org.casadi.Function;
import org.casadi.GenericMatrix;
import org.casadi.MX;
import org.casadi.SX;
import org.casadi.Slice;
import org.casadi.ca;

/*
Solves the following optimal control problem (OCP) in differential-algebraic
equations (DAE)

minimize     integral_{t=0}^{10} x0^2 + x1^2 + u^2  dt
x0,x1,z,u

subject to   dot(x0) == z*x0-x1+u     \
             dot(x1) == x0             }  for 0 <= t <= 10
                   0 == x1^2 + z - 1  /
             x0(t=0) == 0
             x1(t=0) == 1
             x0(t=10) == 0
             x1(t=10) == 0
             -0.75 <= u <= 1  for 0 <= t <= 10

The method used is direct multiple shooting.

Joel Andersson, 2015
*/
public class dae_multiple_shooting {
  public static void main(String[] args) {
    if (!ca.has_integrator("idas")) { System.out.println("SKIP: integrator idas not available"); return; }
    if (!ca.has_nlpsol("ipopt")) { System.out.println("SKIP: nlpsol ipopt not available"); return; }
    if (!hasLibrary("hsl")) { System.out.println("SKIP: libhsl (ipopt linear_solver ma27) not available"); return; }
    double inf = Double.POSITIVE_INFINITY;

    // Declare variables
    SX x0 = SX.sym("x0");
    SX x1 = SX.sym("x1");
    SX x = ca.vertcat(List.of(x0, x1)); // Differential states
    SX z = SX.sym("z");       // Algebraic variable
    SX u = SX.sym("u");       // Control

    // Differential equation
    SX f_x = ca.vertcat(List.of(z.times(x0).minus(x1).plus(u), x0));

    // Algebraic equation
    SX f_z = x1.power(2).plus(z).minus(1);

    // Lagrange cost term (quadrature)
    SX f_q = x0.power(2).plus(x1.power(2)).plus(u.power(2));

    // Create an integrator
    Map<String, SX> dae = Map.of("x", x, "z", z, "p", u, "ode", f_x, "alg", f_z, "quad", f_q);
    // interval length 0.5s
    Function I = ca.integrator("I", "idas", dae, 0, 0.5);

    // Number of intervals
    int nk = 20;

    // Start with an empty NLP
    List<MX> w = new ArrayList<>();         // List of variables
    List<Double> lbw = new ArrayList<>();   // Lower bounds on w
    List<Double> ubw = new ArrayList<>();   // Upper bounds on w
    List<MX> G = new ArrayList<>();         // Constraints
    MX J = new MX(0);                       // Cost function

    // Initial conditions
    MX Xk = MX.sym("X0", 2);
    w.add(Xk);
    lbw.addAll(List.of(0.0, 1.0));
    ubw.addAll(List.of(0.0, 1.0));

    // Loop over all intervals
    for (int k = 0; k < nk; k++) {
      // Local control
      MX Uk = MX.sym("U" + k);
      w.add(Uk);
      lbw.add(-0.75);
      ubw.add(1.00);

      // Call integrator function
      Map<String, MX> Ik = I.call(Map.of("x0", Xk, "p", Uk));
      Xk = Ik.get("xf");
      J = J.plus(Ik.get("qf")); // Sum quadratures

      // "Lift" the variable
      MX X_prev = Xk;
      Xk = MX.sym("X" + (k + 1), 2);
      w.add(Xk);
      lbw.addAll(List.of(-inf, -inf));
      ubw.addAll(List.of(inf, inf));
      G.add(X_prev.minus(Xk));
    }

    // Allocate an NLP solver
    Map<String, MX> nlp = Map.of("x", ca.vertcat(w), "f", J, "g", ca.vertcat(G));
    Map<String, Object> opts = Map.of("ipopt.linear_solver", "ma27");
    Function solver = ca.nlpsol("solver", "ipopt", nlp, opts);

    // Pass bounds, initial guess and solve NLP
    Map<String, GenericMatrix> sol = solver.invoke(Map.of(
        "lbx", toArray(lbw), // Lower variable bound
        "ubx", toArray(ubw), // Upper variable bound
        "lbg", 0.0,          // Lower constraint bound
        "ubg", 0.0,          // Upper constraint bound
        "x0", 0.0));         // Initial guess

    // Plot the results
    DM x_opt = (DM) sol.get("x");
    DM tgrid = ca.linspace(0., 10., nk + 1);
    DM x0_opt = x_opt.get(new Slice(0, x_opt.numel(), 3));
    DM x1_opt = x_opt.get(new Slice(1, x_opt.numel(), 3));
    DM u_opt = x_opt.get(new Slice(2, x_opt.numel(), 3));
    // plot: x0, x1 (vs tgrid) and u trajectories, 'Van der Pol optimization - multiple shooting'
  }

  static boolean hasLibrary(String name) {
    try {
      System.loadLibrary(name);
      return true;
    } catch (UnsatisfiedLinkError e) {
      return false;
    }
  }

  static double[] toArray(List<Double> v) {
    return v.stream().mapToDouble(Double::doubleValue).toArray();
  }
}

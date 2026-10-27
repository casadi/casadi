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
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

public class rocket {
  public static void main(String[] args) {
    if (!ca.has_nlpsol("ipopt")) { System.out.println("SKIP: nlpsol ipopt not available"); return; }

    // Control
    MX u = MX.sym("u");

    // State
    MX x = MX.sym("x", 3);
    MX s = x.get(0); // position
    MX v = x.get(1); // speed
    MX m = x.get(2); // mass

    // ODE right hand side
    MX sdot = v;
    MX vdot = u.minus(ca.times(0.05, v).times(v)).div(m);
    MX mdot = ca.times(-0.1, u).times(u);
    MX xdot = ca.vertcat(List.of(sdot, vdot, mdot));

    // ODE right hand side function
    Function f = new Function("f", List.of(x, u), List.of(xdot));

    // Integrate with Explicit Euler over 0.2 seconds
    double dt = 0.01;  // Time step
    MX xj = x;
    for (int j = 0; j < 20; j++) {
      MX fj = (MX) f.invoke(xj, u);
      xj = xj.plus(ca.times(dt, fj));
    }

    // Discrete time dynamics function
    Function F = new Function("F", List.of(x, u), List.of(xj));

    // Number of control segments
    int nu = 50;

    // Control for all segments
    MX U = MX.sym("U", nu);

    // Initial conditions
    MX X0 = new MX(new DM(new double[] {0, 0, 1}));

    // Integrate over all intervals
    MX X = X0;
    for (int k = 0; k < nu; k++) {
      X = (MX) F.invoke(X, U.get(k));
    }

    // Objective function and constraints
    MX J = ca.mtimes(U.T(), U); // u'*u in Matlab
    MX G = X.get(new Slice(0, 2));   // x(1:2) in Matlab

    // NLP
    Map<String, MX> nlp = Map.of("x", U, "f", J, "g", G);

    // Allocate an NLP solver
    Map<String, Object> opts = Map.of("ipopt.tol", 1e-10, "expand", true);
    Function solver = ca.nlpsol("solver", "ipopt", nlp, opts);
    Map<String, Object> arg = new LinkedHashMap<>();

    // Bounds on u and initial condition
    arg.put("lbx", -0.5);
    arg.put("ubx",  0.5);
    arg.put("x0",   0.4);

    // Bounds on g
    arg.put("lbg", new double[] {10, 0});
    arg.put("ubg", new double[] {10, 0});

    // Solve the problem
    Map<String, GenericMatrix> res = solver.invoke(arg);

    // Get the solution
    // plot: res.get("x") and res.get("lam_x")
  }
}

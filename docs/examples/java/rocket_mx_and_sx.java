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
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.casadi.DM;
import org.casadi.Function;
import org.casadi.MX;
import org.casadi.SX;
import org.casadi.ca;

public class rocket_mx_and_sx {
  // CONSTRUCT THE INTEGRATOR
  static Function create_integrator(int nj, int nu) {
    SX u = SX.sym("u"); // control for one segment

    // Initial position
    SX s0 = SX.sym("s0"); // initial position
    SX v0 = SX.sym("v0"); // initial speed
    SX m0 = SX.sym("m0"); // initial mass

    SX dt = new SX(10.0 / (nj * nu)); // time step
    SX alpha = new SX(0.05); // friction
    SX beta = new SX(0.1); // fuel consumption rate

    // Integrate over the interval with Euler forward
    SX s = s0, v = v0, m = m0;

    SX dm = dt.unaryMinus().times(beta).times(u).times(u);
    for (int j = 0; j < nj; ++j) {
      s = s.plus(dt.times(v));
      v = v.plus(dt.div(m).times(u.minus(alpha.times(v).times(v))));
      m = m.plus(dm);
    }

    // State vector
    SX x = ca.vertcat(s, v, m);
    SX x0 = ca.vertcat(s0, v0, m0);

    // Integrator
    return new Function("integrator", List.of(u, x0), List.of(x));
  }

  public static void main(String[] args) {
    if (!ca.has_nlpsol("ipopt")) { System.out.println("SKIP: nlpsol ipopt not available"); return; }

    // Dimensions
    int nj = 1000; // Number of integration steps per control segment
    int nu = 1000; // Number of control segments

    // Create an integrator
    Function integrator = create_integrator(nj, nu);

    // PART 2: CONSTRUCT THE NLP
    MX U = MX.sym("U", nu); // control for all segments

    // Initial position
    double[] X0 = new double[3];
    X0[0] = 0; // initial position
    X0[1] = 0; // initial speed
    X0[2] = 1; // initial mass

    // Integrate over all intervals
    MX X = new MX(X0);
    for (int k = 0; k < nu; ++k) {
      // Assemble the input
      List<MX> input = new ArrayList<>(List.of(new MX(), new MX()));
      input.set(0, U.get(k));
      input.set(1, X);

      // Integrate
      X = integrator.call(input).get(0);
    }

    // Objective function
    MX F = ca.dot(U, U);

    // Terminal constraints
    MX G = ca.vertcat(X.get(0), X.get(1));

    // Create the NLP
    Map<String, MX> nlp = Map.of("x", U, "f", F, "g", G);

    // Allocate an NLP solver and buffers
    Map<String, Object> opts = Map.of("ipopt.tol", 1e-10,
                                      "ipopt.hessian_approximation", "limited-memory");
    Function solver = ca.nlpsol("solver", "ipopt", nlp, opts);
    Map<String, DM> arg = new HashMap<>(), res;

    // Bounds on u and initial condition
    arg.put("lbx", new DM(-10));
    arg.put("ubx", new DM(10));
    arg.put("x0", new DM(0.4));

    // Bounds on g
    double[] Gmin = new double[2], Gmax = new double[2];
    Gmin[0] = Gmax[0] = 10;
    Gmin[1] = Gmax[1] =  0;
    arg.put("lbg", new DM(Gmin));
    arg.put("ubg", new DM(Gmax));

    // Solve the problem
    res = solver.call(arg);

    // Get the solution
    System.out.println("optimal solution: " + res.get("x"));
  }
}

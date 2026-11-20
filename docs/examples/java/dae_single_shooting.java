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
import java.util.List;
import java.util.Map;

/*
This example mainly intended for CasADi presentations. It contains a compact
implementation of a direct single shooting method for DAEs using a minimal
number of CasADi concepts.

It solves the following optimal control problem (OCP) in differential-algebraic
equations (DAE):

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

Note that other methods such as direct collocation or direct multiple shooting
are usually preferably to the direct single shooting method in practise.

Joel Andersson, 2012-2015
*/

public class dae_single_shooting {
  public static void main(String[] args) {
    if (!ca.has_integrator("idas")) { System.out.println("SKIP: integrator idas not available"); return; }
    if (!ca.has_nlpsol("ipopt")) { System.out.println("SKIP: nlpsol ipopt not available"); return; }

    // Declare variables
    SX x = SX.sym("x", 2); // Differential states
    SX z = SX.sym("z");    // Algebraic variable
    SX u = SX.sym("u");    // Control

    // Differential equation
    SX f_x = ca.vertcat(List.of(z.times(x.get(0)).minus(x.get(1)).plus(u), x.get(0)));

    // Algebraic equation
    SX f_z = x.get(1).power(2).plus(z).minus(1);

    // Lagrange cost term (quadrature)
    SX f_q = x.get(0).power(2).plus(x.get(1).power(2)).plus(u.power(2));

    // Create an integrator
    Map<String, SX> dae = Map.of("x", x, "z", z, "p", u, "ode", f_x, "alg", f_z, "quad", f_q);
    // interval length 0.5s
    Function I = ca.integrator("I", "idas", dae, 0, 0.5);

    // All controls
    MX U = MX.sym("U", 20);

    // Construct graph of integrator calls
    MX X = new MX(new DM(new double[] {0, 1}));
    MX J = new MX(0);
    for (int k = 0; k < 20; k++) {
      Map<String, MX> Ik = I.call(Map.of("x0", X, "p", U.get(k)));
      X = Ik.get("xf");
      J = J.plus(Ik.get("qf"));   // Sum up quadratures
    }

    // Allocate an NLP solver
    Map<String, MX> nlp = Map.of("x", U, "f", J, "g", X);
    Map<String, Object> opts = Map.of("ipopt.linear_solver", "ma27");
    Function solver = ca.nlpsol("solver", "ipopt", nlp, opts);

    // Pass bounds, initial guess and solve NLP
    Map<String, GenericMatrix> sol = solver.invoke(Map.of(
        "lbx", -0.75, // Lower variable bound
        "ubx",  1.0,  // Upper variable bound
        "lbg",  0.0,  // Lower constraint bound
        "ubg",  0.0,  // Upper constraint bound
        "x0",   0.0)); // Initial guess
    System.out.println(sol);
  }
}

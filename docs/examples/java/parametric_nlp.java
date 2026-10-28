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
import java.util.Map;

/**
 *  Example program demonstrating parametric NLPs in CasADi
 *  Note that there is currently no support for parametric sensitivities via this feature (although it would make a lot of sense).
 *  For parametric sensitivities, see the parametric_sensitivities.cpp example which calculates sensitivitied via the sIPOPT extension
 *  to IPOPT.
 *
 *  Joel Andersson, KU Leuven 2012
 */

public class parametric_nlp {
  public static void main(String[] args) {
    if (!ca.has_nlpsol("ipopt")) { System.out.println("SKIP: nlpsol ipopt not available"); return; }

    /* Test problem (Ganesh & Biegler, A reduced Hessian strategy for sensitivity analysis of optimal flowsheets, AIChE 33, 1987, pp. 282-296)
     *
     *    min     x1^2 + x2^2 + x3^2
     *    s.t.    6*x1 + 3&x2 + 2*x3 - pi = 0
     *            p2*x1 + x2 - x3 - 1 = 0
     *            x1, x2, x3 >= 0
     *
     */

    // Optimization variables
    SX x = SX.sym("x", 3);

    // Parameters
    SX p = SX.sym("p", 2);

    // Objective
    SX f = x.get(0).times(x.get(0)).plus(x.get(1).times(x.get(1))).plus(x.get(2).times(x.get(2)));

    // Constraints
    SX g = ca.vertcat(
        ca.times(6, x.get(0)).plus(ca.times(3, x.get(1))).plus(ca.times(2, x.get(2))).minus(p.get(0)),
        p.get(1).times(x.get(0)).plus(x.get(1)).minus(x.get(2)).minus(1)
    );

    // Initial guess and bounds for the optimization variables
    double[] x0  = {0.15, 0.15, 0.00};
    double[] lbx = {0.00, 0.00, 0.00};
    double[] ubx = {ca.inf, ca.inf, ca.inf};

    // Nonlinear bounds
    double[] lbg = {0.00, 0.00};
    double[] ubg = {0.00, 0.00};

    // Original parameter values
    double[] p0  = {5.00, 1.00};

    // NLP
    Map<String, SX> nlp = Map.of("x", x, "p", p, "f", f, "g", g);

    // Create NLP solver and buffers
    Function solver = ca.nlpsol("solver", "ipopt", nlp);
    Map<String, DM> arg = new LinkedHashMap<>();
    Map<String, GenericMatrix> res;

    // Solve the NLP
    arg.put("lbx", new DM(lbx));
    arg.put("ubx", new DM(ubx));
    arg.put("lbg", new DM(lbg));
    arg.put("ubg", new DM(ubg));
    arg.put("x0", new DM(x0));
    arg.put("p", new DM(p0));
    res = solver.invoke(arg);

    // Print the solution
    System.out.println("-----");
    System.out.println("Optimal solution for p = " + arg.get("p") + ":");
    System.out.println(String.format("%30s", "Objective: ") + res.get("f"));
    System.out.println(String.format("%30s", "Primal solution: ") + res.get("x"));
    System.out.println(String.format("%30s", "Dual solution (x): ") + res.get("lam_x"));
    System.out.println(String.format("%30s", "Dual solution (g): ") + res.get("lam_g"));

    // Change the parameter and resolve
    arg.put("p", new DM(4.5));
    res = solver.invoke(arg);

    // Print the new solution
    System.out.println("-----");
    System.out.println("Optimal solution for p = " + arg.get("p") + ":");
    System.out.println(String.format("%30s", "Objective: ") + res.get("f"));
    System.out.println(String.format("%30s", "Primal solution: ") + res.get("x"));
    System.out.println(String.format("%30s", "Dual solution (x): ") + res.get("lam_x"));
    System.out.println(String.format("%30s", "Dual solution (g): ") + res.get("lam_g"));
  }
}

//
//     MIT No Attribution
//
//     Copyright (C) 2010-2026 Joel Andersson, Joris Gillis, Moritz Diehl, KU Leuven.
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
import java.util.Map;

/*
Solve the hs015 problem, formulated as the NLP:

minimize     100*(x[1] - x[0]**2)**2 + (1 - x[0])**2
subject to   x[0]*x[1] >= 1
             x[0] + x[1]**2 >= 0
             x[0] <= 1/2

David Kiessling, 2026
*/

public class hs015 {
  public static void main(String[] args) {
    if (!ca.has_nlpsol("uno")) { System.out.println("SKIP: nlpsol uno not available"); return; }
    double inf = Double.POSITIVE_INFINITY;

    // Declare variable x
    MX x = MX.sym("x", 2);

    // Formulate the NLP
    MX f = x.get(1).minus(x.get(0).power(2)).power(2).times(100).plus(ca.minus(1, x.get(0)).power(2));
    MX g = MX.zeros(2, 1);
    g.set(0, x.get(0).times(x.get(1)));
    g.set(1, x.get(0).plus(x.get(1).power(2)));

    DM ubx = new DM(new double[] {1.0 / 2, inf});
    DM lbx = new DM(new double[] {-inf, -inf});
    DM lbg = new DM(new double[] {1, 0});
    DM ubg = new DM(new double[] {inf, inf});
    Map<String, MX> nlp = Map.of("x", x, "f", f, "g", g);

    // Create an NLP solver
    Function solver = ca.nlpsol("solver", "uno", nlp, Map.of("uno", Map.of("preset", "ipopt")));

    // Solve the problem
    Map<String, GenericMatrix> res = solver.invoke(Map.of(
        "x0", new double[] {-2.0, 1.0},
        "ubg", ubg,
        "lbg", lbg,
        "ubx", ubx));

    // Print solution
    System.out.println();
    System.out.println(String.format("%50s  ", "Optimal cost:") + res.get("f"));
    System.out.println(String.format("%50s  ", "Primal solution:") + res.get("x"));
    System.out.println(String.format("%50s  ", "Dual solution (simple bounds):") + res.get("lam_x"));
    System.out.println(String.format("%50s  ", "Dual solution (nonlinear bounds):") + res.get("lam_g"));
  }
}

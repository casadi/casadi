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
Solve the Rosenbrock problem, formulated as the NLP:

minimize     x^2 + 100*z^2
subject to   z+(1-x)^2-y == 0

Joel Andersson, 2015
*/

public class rosenbrock {
  public static void main(String[] args) {
    if (!ca.has_nlpsol("ipopt")) { System.out.println("SKIP: nlpsol ipopt not available"); return; }

    // Declare variables
    SX x = SX.sym("x");
    SX y = SX.sym("y");
    SX z = SX.sym("z");

    // Formulate the NLP
    SX f = x.power(2).plus(z.power(2).times(100));
    SX g = z.plus(ca.minus(1, x).power(2)).minus(y);
    Map<String, SX> nlp = Map.of("x", ca.vertcat(List.of(x, y, z)), "f", f, "g", g);

    // Create an NLP solver
    Function solver = ca.nlpsol("solver", "ipopt", nlp);

    // Solve the Rosenbrock problem
    Map<String, GenericMatrix> res = solver.invoke(Map.of(
        "x0", new double[] {2.5, 3.0, 0.75},
        "ubg", 0,
        "lbg", 0));

    // Print solution
    System.out.println();
    System.out.println(String.format("%50s  ", "Optimal cost:") + res.get("f"));
    System.out.println(String.format("%50s  ", "Primal solution:") + res.get("x"));
    System.out.println(String.format("%50s  ", "Dual solution (simple bounds):") + res.get("lam_x"));
    System.out.println(String.format("%50s  ", "Dual solution (nonlinear bounds):") + res.get("lam_g"));
  }
}

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
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.casadi.DM;
import org.casadi.Function;
import org.casadi.SX;
import org.casadi.ca;

/*
Solve the Rosenbrock problem, formulated as the NLP:

minimize     x^2 + 100*z^2
subject to   z+(1-x)^2-y == 0

Joel Andersson, 2015-2016
*/

public class rosenbrock_madnlp {
  public static void main(String[] args) {
    if (!ca.has_nlpsol("madnlp")) { System.out.println("SKIP: nlpsol madnlp not available"); return; }

    // Declare variables
    SX x = SX.sym("x");
    SX y = SX.sym("y");
    SX z = SX.sym("z");

    // Formulate the NLP
    SX f = x.power(2).plus(z.power(2).times(100));
    SX g = z.plus(ca.minus(1, x).power(2)).minus(y);
    Map<String, SX> nlp = Map.of("x", ca.vertcat(List.of(x, y, z)),
                                 "f", f,
                                 "g", g);

    Map<String, Object> solver_opts = new HashMap<>();
    //solver_opts.put("madnlp.print_level", 1);
    solver_opts.put("madnlp.linear_solver", "MumpsSolver");
    solver_opts.put("madnlp.tol", 1e-4);
    // Create an NLP solver
    Function solver = ca.nlpsol("solver", "madnlp", nlp, solver_opts);

    // Solve the Rosenbrock problem
    Map<String, DM> arg = new HashMap<>();
    arg.put("x0", new DM(new double[] {2.5, 3.0, 0.75}));
    arg.put("lbg", new DM(0));
    arg.put("ubg", new DM(0));
    Map<String, DM> res = solver.call(arg);

    //  Print solution
    System.out.println("Optimal cost:                     " + str(res.get("f").to_double()));
    System.out.println("Primal solution:                  " + str(res.get("x").nonzeros()));
    System.out.println("Dual solution (simple bounds):    " + str(res.get("lam_x").nonzeros()));
    System.out.println("Dual solution (nonlinear bounds): " + str(res.get("lam_g").nonzeros()));
  }

  // C++ str(): 6 significant digits, as std::ostream
  static String str(double x) { return new DM(x).toString(); }
  static String str(double[] x) {
    StringBuilder r = new StringBuilder("[");
    for (int i = 0; i < x.length; ++i) r.append(i > 0 ? ", " : "").append(str(x[i]));
    return r.append("]").toString();
  }
}

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
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
Solve a basic mpcc with ccopt

minimize     (x-2)^2 + (y-2)^2
subject to   0 <= x _|_ y >= 0

*/

public class basic_mpcc {
  public static void main(String[] args) {
    if (!ca.has_nlpsol("ccopt")) { System.out.println("SKIP: nlpsol ccopt not available"); return; }

    // Declare variables
    SX x = SX.sym("x");
    SX y = SX.sym("y");

    // Formulate the NLP
    SX f = ca.power(x.minus(2), 2).plus(ca.power(y.minus(2), 2));
    Map<String, SX> nlp = Map.of("x", ca.vertcat(List.of(x, y)),
                                 "f", f);

    Map<String, Object> solver_opts = new HashMap<>();
    //solver_opts.put("madnlp.print_level", 1);
    solver_opts.put("cc_pairs", new long[][] {{0, 1}});
    solver_opts.put("madnlp.bound_relax_factor", 0.0);
    // Create an NLP solver
    Function solver = ca.nlpsol("solver", "ccopt", nlp, solver_opts);

    // Solve the Rosenbrock problem
    Map<String, Object> arg = new HashMap<>();
    arg.put("x0", new double[] {2.5, 3.0});
    arg.put("lbx", new double[] {0.0, 0.0});
    Map<String, GenericMatrix> res = solver.invoke(arg);

    //  Print solution
    System.out.println("Optimal cost:                     " + ((DM) res.get("f")).to_double());
    System.out.println("Primal solution:                  " + res.get("x"));
    System.out.println("Dual solution (simple bounds):    " + res.get("lam_x"));
  }
}

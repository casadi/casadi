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
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import org.casadi.DM;
import org.casadi.Function;
import org.casadi.GenericMatrix;
import org.casadi.SX;
import org.casadi.ca;

/** Solves the following NLP with block structured Hessian:
    minimize      x^2 - 0.5*y^2
    subject to    x - y = 0

    Joel Andersson, 2016
*/
public class test_blocksqp {
  public static void main(String[] args) {
    if (!ca.has_nlpsol("blocksqp")) { System.out.println("SKIP: nlpsol blocksqp not available"); return; }

    // Declare variables
    SX x = SX.sym("x");
    SX y = SX.sym("y");

    // Formulate the NLP
    SX f = x.power(2).minus(y.power(2).times(0.5));
    SX g = x.minus(y);
    Map<String, SX> nlp = Map.of("x", ca.vertcat(List.of(x, y)),
                                 "f", f,
                                 "g", g);

    // Create an NLP solver
    Map<String, Object> opts = new LinkedHashMap<>();
    opts.put("opttol", 1.0e-12);
    opts.put("nlinfeastol", 1.0e-12);
    opts.put("globalization", 0);
    opts.put("hess_update", 0);
    opts.put("hess_scaling", 0);
    opts.put("fallback_scaling", 0);
    opts.put("hess_lim_mem", 0);
    opts.put("max_consec_skipped_updates", 200);
    opts.put("block_hess", 0);
    opts.put("which_second_derv", 0);
    opts.put("schur", false);
    Function solver = ca.nlpsol("solver", "blocksqp", nlp, opts);

    // Solve the Rosenbrock problem
    Map<String, Object> arg = new LinkedHashMap<>();
    arg.put("x0", new double[] {10, 10});
    arg.put("lbg", 0);
    arg.put("ubg", 0);
    Map<String, GenericMatrix> res = solver.invoke(arg);

    //  Print solution
    System.out.println("Optimal cost:                     " + ((DM) res.get("f")).to_double());
    System.out.println("Primal solution:                  " + res.get("x"));
    System.out.println("Dual solution (simple bounds):    " + res.get("lam_x"));
    System.out.println("Dual solution (nonlinear bounds): " + res.get("lam_g"));
  }
}

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

public class simple_nlp {
  public static void main(String[] args) {
    // Pick an NLP solver
    String MySolver = "ipopt";
    //String MySolver = "worhp";
    //String MySolver = "sqpmethod";
    if (!ca.has_nlpsol(MySolver)) { System.out.println("SKIP: nlpsol " + MySolver + " not available"); return; }

    // Declare variables
    SX x = SX.sym("x", 2);

    // Form the NLP
    SX f = x.get(0).power(2).plus(x.get(1).power(2)); // objective
    SX g = x.get(0).plus(x.get(1)).minus(10);         // constraint
    Map<String, SX> nlp = Map.of("x", x, "f", f, "g", g);

    // Solver options
    Map<String, Object> opts = new LinkedHashMap<>();
    if (MySolver.equals("sqpmethod")) {
      opts.put("qpsol", "qpoases");
      opts.put("qpsol_options", Map.of("printLevel", "none"));
    }

    // Allocate a solver
    Function solver = ca.nlpsol("solver", MySolver, nlp, opts);

    // Solve the NLP
    Map<String, GenericMatrix> sol = solver.invoke(Map.of("lbg", 0));

    // Print solution
    System.out.println("-----");
    System.out.println("objective at solution =  " + sol.get("f"));
    System.out.println("primal solution =  " + sol.get("x"));
    System.out.println("dual solution (x) =  " + sol.get("lam_x"));
    System.out.println("dual solution (g) =  " + sol.get("lam_g"));
  }
}

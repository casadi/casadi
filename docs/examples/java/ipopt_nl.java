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

/*
This example demonstrates how NL-files, which can be generated
by AMPl or Pyomo, can be imported in CasADi and solved using
e.g. the interface to AMPL

Joel Andersson
2012
*/

public class ipopt_nl {
  public static void main(String[] args) {
    if (!ca.has_nlpsol("ipopt")) { System.out.println("SKIP: nlpsol ipopt not available"); return; }

    // Create an NLP instance
    NlpBuilder nl = new NlpBuilder();

    // Parse an NL-file
    nl.import_nl("../nl_files/hs107.nl", Map.of("verbose", false));

    // NLP solver options
    Map<String, Object> opts = new LinkedHashMap<>();
    opts.put("expand", true);
    // opts.put("verbose", true);
    // opts.put("ipopt", Map.of("max_iter", 10, "linear_solver", "ma57", "hessian_approximation", "limited-memory"));

    // Create an NLP solver
    Function nlpsol = ca.nlpsol("nlpsol", "ipopt", nl, opts);

    // Solve NLP
    Map<String, GenericMatrix> res = nlpsol.invoke(Map.of(
        "lbx", nl.getX_lb(),
        "ubx", nl.getX_ub(),
        "lbg", nl.getG_lb(),
        "ubg", nl.getG_ub(),
        "x0", nl.getX_init()));
  }
}

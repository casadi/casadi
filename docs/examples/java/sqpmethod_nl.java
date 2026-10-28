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
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.Map;
import java.util.TreeMap;

/*
 * This example demonstrates how NL-files, which can be generated
 * by AMPl or Pyomo, can be imported in CasADi and solved using
 * e.g. the interface to AMPL

 \author Joel Andersson, Vyacheslav Kungurtsev
 \date 2013
*/

public class sqpmethod_nl {
  public static void main(String[] args) {
    if (!ca.has_nlpsol("sqpmethod")) { System.out.println("SKIP: nlpsol sqpmethod not available"); return; }

    // Get the problem
    String problem = (args.length == 1) ? args[0] : "../nl_files/hs107.nl";

    // Parse an NL-file
    NlpBuilder nl = new NlpBuilder();
    nl.import_nl(problem);

    // Set options
    Map<String, Object> opts = new LinkedHashMap<>();
    opts.put("expand", true);
    // opts.put("max_iter", 10);
    // opts.put("verbose", true);
    // opts.put("linear_solver", "ma57");
    opts.put("hessian_approximation", "exact");
    // opts.put("derivative_test", "second-order");

    // Specify QP solver
    opts.put("qpsol", "nlpsol");
    opts.put("qpsol_options.nlpsol", "ipopt");
    opts.put("qpsol_options.error_on_fail", false);
    opts.put("qpsol_options.nlpsol_options.ipopt.print_level", 0);
    opts.put("qpsol_options.nlpsol_options.print_time", 0);


    // Allocate NLP solver and buffers
    Function solver = ca.nlpsol("nlpsol", "sqpmethod", nl, opts);
    Map<String, Object> arg = new LinkedHashMap<>();
    Map<String, GenericMatrix> res;

    // Solve NLP
    arg.put("lbx", nl.getX_lb());
    arg.put("ubx", nl.getX_ub());
    arg.put("lbg", nl.getG_lb());
    arg.put("ubg", nl.getG_ub());
    arg.put("x0", nl.getX_init());
    res = solver.invoke(arg);
  }
}

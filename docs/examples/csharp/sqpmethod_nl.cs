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
using System;
using System.Collections.Generic;
using CasADi;

/*
 * This example demonstrates how NL-files, which can be generated
 * by AMPl or Pyomo, can be imported in CasADi and solved using
 * e.g. the interface to AMPL

 \author Joel Andersson, Vyacheslav Kungurtsev
 \date 2013
*/

public static class sqpmethod_nl {
  public static void Main(string[] args) {
    if (!ca.has_nlpsol("sqpmethod")) { Console.WriteLine("SKIP: nlpsol sqpmethod not available"); return; }
    // The QP subproblems are solved with ipopt (not checked by the Java example)
    if (!ca.has_nlpsol("ipopt")) { Console.WriteLine("SKIP: nlpsol ipopt not available"); return; }

    // Get the problem
    string problem = args.Length == 1 ? args[0] : "../nl_files/hs107.nl";

    // Parse an NL-file
    var nl = ca.NlpBuilder();
    nl.import_nl(problem);

    // Set options
    var opts = new Dictionary<string, object>();
    opts["expand"] = true;
    // opts["max_iter"] = 10L;
    // opts["verbose"] = true;
    // opts["linear_solver"] = "ma57";
    opts["hessian_approximation"] = "exact";
    // opts["derivative_test"] = "second-order";

    // Specify QP solver
    opts["qpsol"] = "nlpsol";
    opts["qpsol_options.nlpsol"] = "ipopt";
    opts["qpsol_options.error_on_fail"] = false;
    opts["qpsol_options.nlpsol_options.ipopt.print_level"] = 0L;
    opts["qpsol_options.nlpsol_options.print_time"] = 0L;


    // Allocate NLP solver and buffers
    Function solver = ca.nlpsol("nlpsol", "sqpmethod", nl, opts);
    var arg = new Dictionary<string, DM>();

    // Solve NLP
    arg["lbx"] = nl.x_lb.ToArray();
    arg["ubx"] = nl.x_ub.ToArray();
    arg["lbg"] = nl.g_lb.ToArray();
    arg["ubg"] = nl.g_ub.ToArray();
    arg["x0"] = nl.x_init.ToArray();
    var res = solver.call(arg);
  }
}

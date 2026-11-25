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
using System.Linq;
using CasADi;

/*
 * This example demonstrates how NL-files, which can be generated
 * by AMPl or Pyomo, can be imported in CasADi and solved using
 * e.g. the interface to AMPL

 \author Joel Andersson, Vyacheslav Kungurtsev
 \date 2013
*/

public static class knitro_nl {
  public static void Main(string[] args) {
    if (Environment.GetEnvironmentVariable("SKIP_KNITRO_TESTS") != null) { Console.WriteLine("SKIP: SKIP_KNITRO_TESTS set"); return; }
    if (!ca.has_nlpsol("knitro")) { Console.WriteLine("SKIP: nlpsol knitro not available"); return; }

    // Get the problem
    string problem = (args.Length == 1) ? args[0] : "../nl_files/hs107.nl";

    // Parse an NL-file
    var nl = ca.NlpBuilder();
    nl.import_nl(problem);

    // Set options
    var opts = new Dictionary<string, object> {
      ["expand"] = true,
      // ["verbose"] = true,
    };

    // Allocate NLP solver and buffers
    Function solver = ca.nlpsol("nlpsol", "knitro", nl, opts);
    var arg = new Dictionary<string, DM>();
    Dictionary<string, DM> res;

    // Structure with bounds and initial guess
    arg["lbx"] = nl.x_lb.ToArray();
    arg["ubx"] = nl.x_ub.ToArray();
    arg["lbg"] = nl.g_lb.ToArray();
    arg["ubg"] = nl.g_ub.ToArray();
    arg["x0"] = nl.x_init.ToArray();

    // Solve the NLP
    res = solver.call(arg);
    foreach (var s in res.OrderBy(e => e.Key, StringComparer.Ordinal)) {
      Console.WriteLine($"{s.Key,10}: [{string.Join(", ", s.Value.nonzeros())}]");
    }
  }
}

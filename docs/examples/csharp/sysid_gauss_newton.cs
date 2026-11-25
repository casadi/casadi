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

public static class sysid_gauss_newton {
  static Function sysid_gauss_newton_(MX e, IDictionary<string, MX> nlp, MX V) {
    // Helper funciton for sysid.m
    MX J = ca.jacobian(e, V);
    MX H = ca.triu(ca.mtimes(J.T(), J));
    MX sigma = MX.sym("sigma");

    var io = new Dictionary<string, MX> { ["x"] = V, ["lam_f"] = sigma, ["hess_gamma_x_x"] = sigma * H };

    var opts = new Dictionary<string, object>();
    Console.WriteLine("WARNING: on Windows, JIT may require a special environment cfr https://github.com/casadi/casadi/wiki/FAQ:-how-to-perform-jit-for-function-evaluations-of-my-optimization-problem%3F");
    opts["jit"] = true;
    opts["compiler"] = "shell";
    opts["jit_options"] = new Dictionary<string, object> { ["verbose"] = true };
    var hessLag = ca.Function("nlp_hess_l", io, ["x", "p", "lam_f", "lam_g"], ["hess_gamma_x_x"], opts);
    opts["hess_lag"] = hessLag;
    Function solver = ca.nlpsol("solver", "ipopt", nlp, opts);
    return solver;
  }

  // The MATLAB source is a helper function without a script; this driver fits a line y = a + b*t with it.
  public static void Main() {
    if (!ca.has_nlpsol("ipopt")) { Console.WriteLine("SKIP: nlpsol ipopt not available"); return; }
    if (!Importer.has_plugin("shell")) { Console.WriteLine("SKIP: importer shell not available"); return; }

    DM t = new double[] { 0, 1, 2, 3 };
    DM y_data = new double[] { 1, 3, 5, 7 };
    MX V = MX.sym("V", 2);
    MX e = y_data - (V[0] + V[1] * t);
    var nlp = new Dictionary<string, MX> { ["x"] = V, ["f"] = 0.5 * ca.dot(e, e) };

    Function solver = sysid_gauss_newton_(e, nlp, V);
    var sol = solver.call(new Dictionary<string, DM> { ["x0"] = new double[] { 0, 0 } });
    Console.WriteLine(sol["x"]);
  }
}

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
import org.casadi.Importer;
import org.casadi.MX;
import org.casadi.ca;

public class sysid_gauss_newton {
  static Function sysid_gauss_newton(MX e, Map<String, MX> nlp, MX V) {
    // Helper funciton for sysid.m
    MX J = ca.jacobian(e, V);
    MX H = ca.triu(ca.mtimes(J.T(), J));
    MX sigma = MX.sym("sigma");

    Map<String, MX> io = new LinkedHashMap<>();
    io.put("x", V);
    io.put("lam_f", sigma);
    io.put("hess_gamma_x_x", sigma.times(H));

    Map<String, Object> opts = new LinkedHashMap<>();
    System.out.println("WARNING: on Windows, JIT may require a special environment cfr https://github.com/casadi/casadi/wiki/FAQ:-how-to-perform-jit-for-function-evaluations-of-my-optimization-problem%3F");
    opts.put("jit", true);
    opts.put("compiler", "shell");
    opts.put("jit_options", Map.of("verbose", true));
    Function hessLag = new Function("nlp_hess_l", io, List.of("x", "p", "lam_f", "lam_g"), List.of("hess_gamma_x_x"), opts);
    opts.put("hess_lag", hessLag);
    Function solver = ca.nlpsol("solver", "ipopt", nlp, opts);
    return solver;
  }

  // The MATLAB source is a helper function without a script; this driver fits a line y = a + b*t with it.
  public static void main(String[] args) {
    if (!ca.has_nlpsol("ipopt")) { System.out.println("SKIP: nlpsol ipopt not available"); return; }
    if (!Importer.has_plugin("shell")) { System.out.println("SKIP: importer shell not available"); return; }

    DM t = new DM(new double[] {0, 1, 2, 3});
    DM y_data = new DM(new double[] {1, 3, 5, 7});
    MX V = MX.sym("V", 2);
    MX e = ca.minus(y_data, V.get(0).plus(V.get(1).times(t)));
    Map<String, MX> nlp = Map.of("x", V, "f", ca.times(0.5, ca.dot(e, e)));

    Function solver = sysid_gauss_newton(e, nlp, V);
    Map<String, GenericMatrix> sol = solver.invoke(Map.of("x0", new double[] {0, 0}));
    System.out.println(sol.get("x"));
  }
}

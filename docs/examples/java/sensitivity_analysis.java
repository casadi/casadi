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
// -*- coding: utf-8 -*-

import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import org.casadi.DM;
import org.casadi.Function;
import org.casadi.GenericMatrix;
import org.casadi.SX;
import org.casadi.ca;

public class sensitivity_analysis {
  // Python print('%50s' % label, a, x, b, y)
  static void print(String label, String a, Object x, String b, Object y) {
    System.out.println(String.format("%50s", label) + " " + a + " " + x + " " + b + " " + y);
  }

  public static void main(String[] args) {
    for (String p : List.of("cvodes", "rk", "idas", "collocation")) {
      if (!ca.has_integrator(p)) { System.out.println("SKIP: integrator " + p + " not available"); return; }
    }
    if (!ca.has_rootfinder("kinsol")) { System.out.println("SKIP: rootfinder kinsol not available"); return; }
    System.out.println("Testing sensitivity analysis in CasADi");

    // All ODE and DAE integrators to be tested
    List<String> DAE_integrators = List.of("idas", "collocation");
    List<String> ODE_integrators = List.of("cvodes", "rk", "idas", "collocation");

    for (List<String> Integrators : List.of(ODE_integrators, DAE_integrators)) {
      Map<String, SX> dae;
      double tf;
      Object x0;
      double u0;
      if (Integrators == ODE_integrators) { // rocket example
        System.out.println("******");
        System.out.println("Testing ODE example");

        // Time
        SX t = SX.sym("t");

        // Parameter
        SX u = SX.sym("u");

        // Differential states
        SX s = SX.sym("s"); SX v = SX.sym("v"); SX m = SX.sym("m");
        SX x = ca.vertcat(List.of(s, v, m));

        // Constants
        double alpha = 0.05; // friction
        double beta = 0.1;   // fuel consumption rate

        // Differential equation
        SX ode = ca.vertcat(List.of(
          v,
          u.minus(ca.times(alpha, v).times(v)).div(m),
          ca.times(-beta, u).times(u)));

        // Quadrature
        SX quad = v.power(3).plus(ca.minus(3, ca.sin(t)).minus(u).power(2));

        // DAE callback function
        dae = Map.of("t", t, "x", x, "p", u, "ode", ode, "quad", quad);

        // Time length
        tf = 0.5;

        // Initial position
        x0 = new double[] {0., 0., 1.};

        // Parameter
        u0 = 0.4;

      } else { // Simple DAE example
        System.out.println("******");
        System.out.println("Testing DAE example");

        // Differential state
        SX x = SX.sym("x");

        // Algebraic variable
        SX z = SX.sym("z");

        // Parameter
        SX u = SX.sym("u");

        // Differential equation
        SX ode = x.unaryMinus().plus(ca.times(0.5, x).times(x)).plus(u).plus(ca.times(0.5, z));

        // Algebraic constraint
        SX alg = z.plus(ca.exp(z)).minus(1.0).plus(x);

        // Quadrature
        SX quad = x.times(x).plus(ca.times(3.0, u).times(u));

        // DAE callback function
        dae = Map.of("x", x, "z", z, "p", u, "ode", ode, "alg", alg, "quad", quad);

        // End time
        tf = 5.;

        // Initial position
        x0 = 1.;

        // Parameter
        u0 = 0.4;
      }

      // Integrator
      for (String MyIntegrator : Integrators) {
        System.out.println("========");
        System.out.println("Integrator:  " + MyIntegrator);
        System.out.println("========");

        // Integrator options
        Map<String, Object> opts = new LinkedHashMap<>();
        if (MyIntegrator.equals("collocation")) {
          opts.put("rootfinder", "kinsol");
        }

        // Integrator
        Function I = ca.integrator("I", MyIntegrator, dae, 0, tf, opts);

        // Integrate to get results
        Map<String, GenericMatrix> res = I.invoke(Map.of("x0", x0, "p", u0));
        DM xf = (DM) res.get("xf");
        DM qf = (DM) res.get("qf");
        print("Unperturbed solution:", "xf  = ", xf, ", qf  = ", qf);

        // Perturb solution to get a finite difference approximation
        double h = 0.001;
        res = I.invoke(Map.of("x0", x0, "p", u0 + h));
        DM fd_xf = ((DM) res.get("xf")).minus(xf).div(h);
        DM fd_qf = ((DM) res.get("qf")).minus(qf).div(h);
        print("Finite difference approximation:", "d(xf)/d(p) = ", fd_xf, ", d(qf)/d(p) = ", fd_qf);

        // Calculate once, forward
        Function I_fwd = I.factory("I_fwd", List.of("x0", "z0", "p", "fwd:p"), List.of("fwd:xf", "fwd:qf"));
        res = I_fwd.invoke(Map.of("x0", x0, "p", u0, "fwd_p", 1));
        DM fwd_xf = (DM) res.get("fwd_xf");
        DM fwd_qf = (DM) res.get("fwd_qf");
        print("Forward sensitivities:", "d(xf)/d(p) = ", fwd_xf, ", d(qf)/d(p) = ", fwd_qf);

        // Calculate once, adjoint
        Function I_adj = I.factory("I_adj", List.of("x0", "z0", "p", "adj:qf"), List.of("adj:x0", "adj:p"));
        res = I_adj.invoke(Map.of("x0", x0, "p", u0, "adj_qf", 1));
        DM adj_x0 = (DM) res.get("adj_x0");
        DM adj_p = (DM) res.get("adj_p");
        print("Adjoint sensitivities:", "d(qf)/d(x0) = ", adj_x0, ", d(qf)/d(p) = ", adj_p);

        // Perturb adjoint solution to get a finite difference approximation of the second order sensitivities
        res = I_adj.invoke(Map.of("x0", x0, "p", u0 + h, "adj_qf", 1));
        DM fd_adj_x0 = ((DM) res.get("adj_x0")).minus(adj_x0).div(h);
        DM fd_adj_p = ((DM) res.get("adj_p")).minus(adj_p).div(h);
        print("FD of adjoint sensitivities:", "d2(qf)/d(x0)d(p) = ", fd_adj_x0, ", d2(qf)/d(p)d(p) = ", fd_adj_p);

        // Forward over adjoint to get the second order sensitivities
        Function I_foa = I_adj.factory("I_foa", List.of("x0", "z0", "p", "adj_qf", "fwd:p"), List.of("fwd:adj_x0", "fwd:adj_p"));
        res = I_foa.invoke(Map.of("x0", x0, "p", u0, "adj_qf", 1, "fwd_p", 1));
        DM fwd_adj_x0 = (DM) res.get("fwd_adj_x0");
        DM fwd_adj_p = (DM) res.get("fwd_adj_p");
        print("Forward over adjoint sensitivities:", "d2(qf)/d(x0)d(p) = ", fwd_adj_x0, ", d2(qf)/d(p)d(p) = ", fwd_adj_p);

        // Adjoint over adjoint to get the second order sensitivities
        Function I_aoa = I_adj.factory("I_aoa", List.of("x0", "z0", "p", "adj_qf", "adj:adj_p"), List.of("adj:x0", "adj:p"));
        res = I_aoa.invoke(Map.of("x0", x0, "p", u0, "adj_qf", 1, "adj_adj_p", 1));
        adj_x0 = (DM) res.get("adj_x0");
        adj_p = (DM) res.get("adj_p");
        print("Adjoint over adjoint sensitivities:", "d2(qf)/d(x0)d(p) = ", adj_x0, ", d2(qf)/d(p)d(p) = ", adj_p);
      }
    }
  }
}

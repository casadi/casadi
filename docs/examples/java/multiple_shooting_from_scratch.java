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
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.casadi.CasadiException;
import org.casadi.DM;
import org.casadi.Function;
import org.casadi.MX;
import org.casadi.SX;
import org.casadi.Slice;
import org.casadi.ca;

/** Writing a multiple shooting code from scratch

  This example demonstrates how to write a simple, yet powerful multiple shooting code from
  scratch using CasADi. For clarity, the code below uses a simple formulation with only states
  and controls, and no path constrants. It relies on CasADi machinery to keep track of sparsity,
  formulate ODE sensitivity equations and build up the Jacobian of the NLP constraint function.

  By extending the code below, it should be possible for a user to solve ever more complex
  problems. For example, one can easily make a multi-stage formulation by simply allocating
  another integrator instance and use the two integrators instances for different shooting nodes.
  By replacing the explicit CVodes integrator with the fully-implicit IDAS integrator, one can
  solve optimal control problems in differential-algebraic equations.

  Joel Andersson, 2011-2012
*/

public class multiple_shooting_from_scratch {
  public static void main(String[] args) {
    if (!ca.has_integrator("cvodes")) { System.out.println("SKIP: integrator cvodes not available"); return; }
    if (!ca.has_nlpsol("ipopt")) { System.out.println("SKIP: nlpsol ipopt not available"); return; }
    double inf = ca.inf;

    // Declare variables
    SX u = SX.sym("u"); // control
    SX r = SX.sym("r"), s = SX.sym("s"); // states
    SX x = ca.vertcat(r, s);

    // Number of differential states
    long nx = x.size1();

    // Number of controls
    long nu = u.size1();

    // Bounds and initial guess for the control
    List<Double> u_min =  List.of( -0.75 );
    List<Double> u_max  = List.of(  1.0  );
    List<Double> u_init = List.of(  0.0  );

    // Bounds and initial guess for the state
    List<Double> x0_min = List.of(   0.0,  1.0 );
    List<Double> x0_max = List.of(   0.0,  1.0 );
    List<Double> x_min  = List.of(  -inf, -inf );
    List<Double> x_max  = List.of(   inf,  inf );
    List<Double> xf_min = List.of(   0.0,  0.0 );
    List<Double> xf_max = List.of(   0.0,  0.0 );
    List<Double> x_init = List.of(   0.0,  0.0 );

    // Final time
    double tf = 20.0;

    // Number of shooting nodes
    int ns = 50;

    // ODE right hand side and quadrature
    SX ode = ca.vertcat(ca.minus(1, s.times(s)).times(r).minus(s).plus(u), r);
    SX quad = r.times(r).plus(s.times(s)).plus(u.times(u));
    Map<String, SX> dae = Map.of("x", x, "p", u, "ode", ode, "quad", quad);

    // Create an integrator (CVodes)
    Function F = ca.integrator("integrator", "cvodes", dae, 0, tf / ns);

    // Total number of NLP variables
    long NV = nx * (ns + 1) + nu * ns;

    // Declare variable vector for the NLP
    MX V = MX.sym("V", NV);

    // NLP variable bounds and initial guess
    List<Double> v_min = new ArrayList<>(), v_max = new ArrayList<>(), v_init = new ArrayList<>();

    // Offset in V
    long offset = 0;

    // State at each shooting node and control for each shooting interval
    List<MX> X = new ArrayList<>(), U = new ArrayList<>();
    for (int k = 0; k < ns; ++k) {
      // Local state
      X.add(V.get_nz(new Slice(offset, offset + nx)));
      if (k == 0) {
        v_min.addAll(x0_min);
        v_max.addAll(x0_max);
      } else {
        v_min.addAll(x_min);
        v_max.addAll(x_max);
      }
      v_init.addAll(x_init);
      offset += nx;

      // Local control
      U.add(V.get_nz(new Slice(offset, offset + nu)));
      v_min.addAll(u_min);
      v_max.addAll(u_max);
      v_init.addAll(u_init);
      offset += nu;
    }

    // State at end
    X.add(V.get_nz(new Slice(offset, offset + nx)));
    v_min.addAll(xf_min);
    v_max.addAll(xf_max);
    v_init.addAll(x_init);
    offset += nx;

    // Make sure that the size of the variable vector is consistent with the number of variables that we have referenced
    if (offset != NV) throw new CasadiException("");

    // Objective function
    MX J = new MX(0);

    //Constraint function and bounds
    List<MX> g = new ArrayList<>();

    // Loop over shooting nodes
    for (int k = 0; k < ns; ++k) {
      // Create an evaluation node
      Map<String, MX> I_out = F.call(Map.of("x0", X.get(k), "p", U.get(k)));

      // Save continuity constraints
      g.add(I_out.get("xf").minus(X.get(k + 1)));

      // Add objective function contribution
      J = J.plus(I_out.get("qf"));
    }

    // NLP
    Map<String, MX> nlp = Map.of("x", V, "f", J, "g", ca.vertcat(g));

    // Set options
    Map<String, Object> opts = new HashMap<>();
    opts.put("ipopt.tol", 1e-5);
    opts.put("ipopt.max_iter", 100);

    // Create an NLP solver and buffers
    Function solver = ca.nlpsol("nlpsol", "ipopt", nlp, opts);
    Map<String, DM> arg = new HashMap<>(), res;

    // Bounds and initial guess
    arg.put("lbx", new DM(toArray(v_min)));
    arg.put("ubx", new DM(toArray(v_max)));
    arg.put("lbg", new DM(0));
    arg.put("ubg", new DM(0));
    arg.put("x0", new DM(toArray(v_init)));

    // Solve the problem
    res = solver.call(arg);

    // Optimal solution of the NLP
    double[] V_opt = res.get("x").nonzeros();

    // Get the optimal state trajectory
    double[] r_opt = new double[ns + 1], s_opt = new double[ns + 1];
    for (int i = 0; i <= ns; ++i) {
      r_opt[i] = V_opt[(int) (i * (nx + 1))];
      s_opt[i] = V_opt[(int) (1 + i * (nx + 1))];
    }
    System.out.println("r_opt = ");
    System.out.println(str(r_opt));
    System.out.println("s_opt = ");
    System.out.println(str(s_opt));

    // Get the optimal control
    double[] u_opt = new double[ns];
    for (int i = 0; i < ns; ++i) {
      u_opt[i] = V_opt[(int) (nx + i * (nx + 1))];
    }
    System.out.println("u_opt = ");
    System.out.println(str(u_opt));
  }

  static double[] toArray(List<Double> v) {
    return v.stream().mapToDouble(Double::doubleValue).toArray();
  }

  // C++ str(): 6 significant digits, as std::ostream
  static String str(double x) { return new DM(x).toString(); }
  static String str(double[] x) {
    StringBuilder r = new StringBuilder("[");
    for (int i = 0; i < x.length; ++i) r.append(i > 0 ? ", " : "").append(str(x[i]));
    return r.append("]").toString();
  }
}

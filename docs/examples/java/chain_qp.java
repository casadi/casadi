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
//
// -*- coding: utf-8 -*-
/*
We want to model a chain attached to two supports and hanging in between. Let us discretise
it with N mass points connected by N-1 springs. Each mass i has position (yi,zi), i=1,...,N.
The equilibrium point of the system minimises the potential energy.

The potential energy of each spring is
Vi=D_i/2 * ((y_i-y_{i+1})^2 + (z_i-z_{i+1})^2)

The gravitational potential energy of each mass is
Vg_i = m_i*g0*z_i

The total potential energy is thus given by:

Vchain(y,z) = 1/2*sum{i=1,...,N-1} D_i ((y_i-y_{i+1})^2+(z_i-z_{i+1})^2) + g0 * sum{i=1,...,N} m_i * z_i

where y=[y_1,...,y_N] and z=[z_1,...,z_N]

We wish to solve
minimize{y,z} Vchain(y, z)

Subject to the piecewise linear ground constraints:
z_i >= zin
z_i - 0.1*y_i >= 0.5
*/
import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import org.casadi.DM;
import org.casadi.Function;
import org.casadi.GenericMatrix;
import org.casadi.SX;
import org.casadi.Slice;
import org.casadi.ca;

public class chain_qp {
  public static void main(String[] args) {
    if (!ca.has_conic("qpoases")) { System.out.println("SKIP: conic qpoases not available"); return; }
    double inf = Double.POSITIVE_INFINITY;

    // Constants
    int N = 40;
    double m_i = 40.0 / N;
    double D_i = 70.0 * N;
    double g0 = 9.81;
    //double zmin = -inf; // unbounded
    double zmin = 0.5; // ground

    // Objective function
    SX Vchain = new SX(0);

    // Variables
    List<SX> x = new ArrayList<>();

    // Variable bounds
    List<Double> lbx = new ArrayList<>();
    List<Double> ubx = new ArrayList<>();

    // Constraints
    List<SX> g = new ArrayList<>();

    // Constraint bounds
    List<Double> lbg = new ArrayList<>();
    List<Double> ubg = new ArrayList<>();

    // Loop over all chain elements
    SX y_i = null, z_i = null, y_prev = null, z_prev = null;
    for (int i = 1; i <= N; i++) {
      // Previous point
      if (i > 1) {
        y_prev = y_i;
        z_prev = z_i;
      }

      // Create variables for the (y_i, z_i) coordinates
      y_i = SX.sym("y_" + i);
      z_i = SX.sym("z_" + i);

      // Add to the list of variables
      x.addAll(List.of(y_i, z_i));
      if (i == 1) {
        lbx.addAll(List.of(-2., 1.));
        ubx.addAll(List.of(-2., 1.));
      } else if (i == N) {
        lbx.addAll(List.of(2., 1.));
        ubx.addAll(List.of(2., 1.));
      } else {
        lbx.addAll(List.of(-inf, zmin));
        ubx.addAll(List.of(inf, inf));
      }

      // Spring potential
      if (i > 1) {
        Vchain = Vchain.plus(y_prev.minus(y_i).power(2).plus(z_prev.minus(z_i).power(2)).times(D_i / 2));
      }

      // Graviational potential
      Vchain = Vchain.plus(z_i.times(g0 * m_i));

      // Slanted ground constraints
      g.add(z_i.minus(y_i.times(0.1)));
      lbg.add(0.5);
      ubg.add(inf);
    }

    // Formulate QP
    Map<String, SX> qp = Map.of("x", ca.vertcat(x), "f", Vchain, "g", ca.vertcat(g));

    // Solve with IPOPT
    Function solver = ca.qpsol("solver", "qpoases", qp, Map.of("sparse", true));
    //Function solver = ca.qpsol("solver", "gurobi", qp);
    //Function solver = ca.nlpsol("solver", "ipopt", qp);

    // Get the optimal solution
    Map<String, GenericMatrix> sol = solver.invoke(Map.of("lbx", toArray(lbx), "ubx", toArray(ubx),
        "lbg", toArray(lbg), "ubg", toArray(ubg)));
    DM x_opt = (DM) sol.get("x");
    DM f_opt = (DM) sol.get("f");
    System.out.println("f_opt =  " + f_opt);

    // Retrieve the result
    DM Y0 = x_opt.get(new Slice(0, x_opt.numel(), 2));
    DM Z0 = x_opt.get(new Slice(1, x_opt.numel(), 2));

    // Plot the result
    DM ys = ca.linspace(-2., 2., 100);
    DM zs = ca.plus(0.5, ys.times(0.1));
    // plot: chain (Y0, Z0) and ground line (ys, zs), 'hanging chain QP'
  }

  static double[] toArray(List<Double> v) {
    return v.stream().mapToDouble(Double::doubleValue).toArray();
  }
}

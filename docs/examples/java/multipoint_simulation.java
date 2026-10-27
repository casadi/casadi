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
import java.util.List;
import java.util.Map;

public class multipoint_simulation {
  public static void main(String[] args) {
    if (!ca.has_integrator("cvodes")) { System.out.println("SKIP: integrator cvodes not available"); return; }

    // Consider the following nonlinear system with two states and one control input:
    MX x1 = MX.sym("x1");
    MX x2 = MX.sym("x2");
    MX x = ca.vertcat(List.of(x1, x2));
    MX u = MX.sym("u");
    MX xdot = ca.vertcat(List.of(ca.minus(1, x2.power(2)).times(x1).minus(x2).plus(u), x1));

    // In addition to the dynamics, we are also interested in the sum-of-squares distance from the origin:
    MX y = x1.power(2).plus(x2.power(2));

    // We want to simulate this system over T seconds, with the control discretized into N intervals:
    // Time horizon, discretization
    int N = 10;
    double T = 10.;
    double[] tgrid = linspace(0, T, N + 1);

    // We will integrate the problem using CVODES, while also calculating the integral of y
    Map<String, MX> dae = Map.of("x", x, "u", u, "ode", xdot, "quad", y);
    Function F = ca.integrator("F", "cvodes", dae, tgrid[0], Arrays.copyOfRange(tgrid, 1, tgrid.length));

    // We are also interested in the pointwise values of y, so let's create a function for that too
    Function yfun = new Function("yfun", List.of(x), List.of(y), List.of("x"), List.of("y"));

    // Let's define initial conditions for x and values for the piecewise constant controls
    double[] u_val = linspace(-1, 1, N);  // Python rebinds u; Java needs a new name
    double[] x0 = {0, 0};

    // We can simulate the system as follows:
    Map<String, GenericMatrix> Fk = F.invoke(Map.of("x0", x0, "u", u_val));
    DM xf = (DM) Fk.get("xf");
    DM qf = (DM) Fk.get("qf");

    // Let's add the state at the initial time to xf and qf
    xf = ca.horzcat(List.of(new DM(x0), xf));
    qf = ca.horzcat(List.of(new DM(0), qf));

    // Call yfun to get y at each grid point
    DM yf = (DM) yfun.invoke(xf);

    // Convert CasADi matrices to Java arrays
    double[][] x_sim = xf.full();
    double[][] y_sim = yf.full();
    double[][] q_sim = qf.full();

    // Plot simulation trajectories and pointwise sum of squares contribution to y, which we compare with the integral of y
    // plot: x_sim[0], x_sim[1], y_sim[0] over tgrid, and step plot of u_val
    double[] y_cumsum = new double[N + 1];
    double acc = 0;
    for (int k = 0; k <= N; k++) { acc += y_sim[0][k]; y_cumsum[k] = T / N * acc; }
    // plot: q_sim[0] and y_cumsum over tgrid (integral cost)
  }

  static double[] linspace(double a, double b, int n) {
    double[] r = new double[n];
    for (int k = 0; k < n; k++) r[k] = a + (b - a) * k / (n - 1);
    return r;
  }
}

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

// This example how an Linear Quadratic Regulator (LQR) can be designed
// and simulated for a linear-time-invariant (LTI) system

import org.casadi.*;
import java.util.*;

public class lqr_control {

  // Eigenvalues of a real 3x3 matrix as {re, im} pairs: roots of the characteristic polynomial
  // l^3 - c2 l^2 + c1 l - c0, with c2 = trace, c1 = sum of principal 2x2 minors, c0 = det
  static double[][] eig3(double[][] a) {
    double c2 = a[0][0] + a[1][1] + a[2][2];
    double c1 = a[0][0] * a[1][1] - a[0][1] * a[1][0] + a[0][0] * a[2][2] - a[0][2] * a[2][0]
        + a[1][1] * a[2][2] - a[1][2] * a[2][1];
    double c0 = a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1])
        - a[0][1] * (a[1][0] * a[2][2] - a[1][2] * a[2][0])
        + a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0]);
    // Depressed cubic t^3 + p t + q with l = t + c2/3
    double s = c2 / 3;
    double p = c1 - c2 * c2 / 3;
    double q = -2 * s * s * s + c1 * s - c0;
    double disc = q * q / 4 + p * p * p / 27;
    if (disc > 0) {
      double sq = Math.sqrt(disc);
      double u = Math.cbrt(-q / 2 + sq), v = Math.cbrt(-q / 2 - sq);
      double re = -(u + v) / 2 + s, im = Math.sqrt(3) / 2 * (u - v);
      return new double[][] {{re, im}, {re, -im}, {u + v + s, 0}};
    }
    double r = 2 * Math.sqrt(-p / 3);
    double phi = Math.acos(Math.max(-1, Math.min(1, 3 * q / (p * r))));
    double[][] l = new double[3][];
    for (int k = 0; k < 3; k++) l[k] = new double[] {r * Math.cos((phi - 2 * Math.PI * k) / 3) + s, 0};
    return l;
  }

  public static void main(String[] args) {
    // Problem statement
    // -----------------------------------

    // System dimensions
    int ns = 3;  // number of states
    int nu = 2;  // number of controls
    int ny = 2;  // number of outputs
    int N = 21;  // number of control intervals
    double te = 10; // end time [s]

    // The system equations:   x' = A.x + B.u
    DM A = new DM(new double[][] {{0.4, 0.1, -2}, {0, -0.3, 4}, {1, 0, 0}});
    DM B = new DM(new double[][] {{1, 1}, {0, 1}, {1, 0}});

    // Inspect the open-loop system
    double[][] D = eig3(A.full());
    StringJoiner sj = new StringJoiner(" ", "[", "]");
    for (double[] d : D) sj.add(String.format(Locale.ROOT, "%.8f%+.8fj", d[0], d[1]));
    System.out.println("Open-loop eigenvalues:  " + sj);

    DM tn = ca.linspace(new DM(0), new DM(te), N);

    // Check if the system is controllable
    // -----------------------------------
    List<DM> Rl = new ArrayList<>();

    DM p = B;
    for (int i = 0; i < ns; i++) {
      Rl.add(p);
      p = ca.mtimes(List.of(A, p));
    }

    DM R = ca.horzcat(Rl);

    // R must be of full rank
    // singular values of R: square roots of the eigenvalues of R R^T
    double[][] RRt = ca.mtimes(R, R.T()).full();
    double eps = 1e-9;
    int rank = 0;
    for (double[] l : eig3(RRt)) if (Math.abs(Math.sqrt(Math.abs(l[0]))) > eps) rank++;
    if (rank != ns) throw new AssertionError();

    // Simulation of the open-loop system
    // -----------------------------------

    SX y = SX.sym("y", ns);
    SX u = SX.sym("u", nu);

    DM x0 = new DM(new double[] {1, 0, 0});
    // no control
    DM u_ = ca.repmat(new DM(new double[][] {{-1, 1}, {1, -1}}), 1, (N - 1) / 2);

    System.out.println(u_);

    SX ps = SX.sym("p");

    // The remainder of the Python example is commented-out code for the removed ControlSimulator API.
  }
}

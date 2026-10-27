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
// This example uses no CasADi functionality: plain Java arrays replace numpy.
public class vdp_dynamic_programming {
  // System dynamics, can be called with matricex
  static double[][][] f(double[][] x1, double[][] x2, double u) {
    int n = x1.length, m = x1[0].length;
    double[][] x1_dot = new double[n][m], q_dot = new double[n][m];
    double[][] x2_dot = x1; // as in numpy, x2_dot is the very array x1 (the in-place updates below see this)
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < m; j++) {
        x1_dot[i][j] = (1 - x2[i][j] * x2[i][j]) * x1[i][j] - x2[i][j] + u;
        q_dot[i][j] = x1[i][j] * x1[i][j] + x2[i][j] * x2[i][j] + u * u;
      }
    }
    return new double[][][] {x1_dot, x2_dot, q_dot};
  }

  // a + c*b, elementwise
  static double[][] axpy(double[][] a, double c, double[][] b) {
    double[][] r = new double[a.length][a[0].length];
    for (int i = 0; i < a.length; i++) for (int j = 0; j < a[0].length; j++) r[i][j] = a[i][j] + c * b[i][j];
    return r;
  }

  // a += c*(k1 + 2*k2 + 2*k3 + k4), in place
  static void rk4_update(double[][] a, double c, double[][] k1, double[][] k2, double[][] k3, double[][] k4) {
    for (int i = 0; i < a.length; i++)
      for (int j = 0; j < a[0].length; j++)
        a[i][j] += c * (k1[i][j] + 2 * k2[i][j] + 2 * k3[i][j] + k4[i][j]);
  }

  static double[][] copy(double[][] a) {
    double[][] r = new double[a.length][];
    for (int i = 0; i < a.length; i++) r[i] = a[i].clone();
    return r;
  }

  static double[] linspace(double a, double b, int n) {
    double[] r = new double[n];
    double step = (b - a) / (n - 1);
    for (int i = 0; i < n; i++) r[i] = i * step + a;
    r[n - 1] = b;
    return r;
  }

  public static void main(String[] args) {
    double inf = Double.POSITIVE_INFINITY;

    // End time
    double T = 10.;

    // Number of control intervals
    int N = 20;

    // Number of Runge-Kutta 4 steps per interval and step size
    int NK = 20;
    double DT = T / (N * NK);

    // Number of discrete control values
    int NU = 101;

    // Number of discrete state values
    int NX = 101;

    // Control enumeration
    double[] U = linspace(-1, 1, NU);

    // State space enumeration
    double[] x1 = linspace(-1, 1, NX);
    double[] x2 = linspace(-1, 1, NX);
    double[][] X1 = new double[NX][NX], X2 = new double[NX][NX];
    for (int i = 0; i < NX; i++) {
      for (int j = 0; j < NX; j++) {
        X1[i][j] = x1[j];
        X2[i][j] = x2[i];
      }
    }

    // For each control action and state, precalculate next state and stage cost
    double[][][] stage_J = new double[NU][][];
    int[][][] next_x1 = new int[NU][][];
    int[][][] next_x2 = new int[NU][][];
    for (int uind = 0; uind < NU; uind++) {
      double u = U[uind];
      // Take number of integration steps
      double[][] X1_k = copy(X1);
      double[][] X2_k = copy(X2);
      double[][] Q_k = new double[NX][NX];
      for (int k = 0; k < NK; k++) {
        // RK4 integration for x1, x2 and q
        double[][][] k1 = f(X1_k, X2_k, u);
        double[][][] k2 = f(axpy(X1_k, DT / 2, k1[0]), axpy(X2_k, DT / 2, k1[1]), u);
        double[][][] k3 = f(axpy(X1_k, DT / 2, k2[0]), axpy(X2_k, DT / 2, k2[1]), u);
        double[][][] k4 = f(axpy(X1_k, DT, k3[0]), axpy(X2_k, DT, k3[1]), u);
        rk4_update(X1_k, DT / 6, k1[0], k2[0], k3[0], k4[0]);
        rk4_update(X2_k, DT / 6, k1[1], k2[1], k3[1], k4[1]);
        rk4_update(Q_k, DT / 6, k1[2], k2[2], k3[2], k4[2]);
      }

      // Find out which state comes next (index)
      int[][] I1 = new int[NX][NX], I2 = new int[NX][NX];
      for (int i = 0; i < NX; i++) {
        for (int j = 0; j < NX; j++) {
          I1[i][j] = (int) Math.rint((X1_k[i][j] + 1) / 2 * (NX - 1));
          I2[i][j] = (int) Math.rint((X2_k[i][j] + 1) / 2 * (NX - 1));

          // Infinite cost if state gets out-of-bounds
          if (I1[i][j] < 0) { Q_k[i][j] = inf; I1[i][j] = 0; }
          if (I2[i][j] < 0) { Q_k[i][j] = inf; I2[i][j] = 0; }
          if (I1[i][j] >= NX) { Q_k[i][j] = inf; I1[i][j] = 0; }
          if (I2[i][j] >= NX) { Q_k[i][j] = inf; I2[i][j] = 0; }
        }
      }

      // Save the stage cost and next state
      next_x1[uind] = I1;
      next_x2[uind] = I2;
      stage_J[uind] = Q_k;
    }

    // Calculate cost-to-go (no end cost) and optimal control
    double[][] J = new double[NX][NX];
    int[][][] U_opt = new int[N][][];
    for (int k = N - 1; k >= 0; k--) {
      // Cost to go for the previous step, optimal control action
      double[][] J_prev = new double[NX][NX];
      int[][] u_prev = new int[NX][NX];
      for (int i = 0; i < NX; i++) {
        java.util.Arrays.fill(J_prev[i], inf);
        java.util.Arrays.fill(u_prev[i], -1);
      }

      // Test all control actions
      for (int uind = 0; uind < NU; uind++) {
        for (int i = 0; i < NX; i++) {
          for (int j = 0; j < NX; j++) {
            double J_prev_test = J[next_x2[uind][i][j]][next_x1[uind][i][j]] + stage_J[uind][i][j];
            if (J_prev_test < J_prev[i][j]) {
              u_prev[i][j] = uind;
              J_prev[i][j] = J_prev_test;
            }
          }
        }
      }

      // Update cost-to-go and save optimal control (stored by stage)
      J = J_prev;
      U_opt[k] = u_prev;
    }

    // Find optimal control starting at x1=0, x2=1
    int i1 = NX / 2;
    int i2 = NX - 1;
    double[] u_opt = new double[N];
    double[] x1_opt = new double[N + 1];
    double[] x2_opt = new double[N + 1];
    x1_opt[0] = x1[i1];
    x2_opt[0] = x2[i2];
    double cost = 0;
    for (int k = 0; k < N; k++) {
      // Get the optimal control and go to next step
      int u_ind = U_opt[k][i2][i1];
      cost += stage_J[u_ind][i2][i1];
      int i1_next = next_x1[u_ind][i2][i1], i2_next = next_x2[u_ind][i2][i1];
      i1 = i1_next;
      i2 = i2_next;

      // Save the trajectories
      u_opt[k] = U[u_ind];
      x1_opt[k + 1] = x1[i1];
      x2_opt[k + 1] = x2[i2];
    }

    // Optimal cost
    System.out.println("Minimal cost:  " + cost);
    if (!(Math.abs(cost - J[NX - 1][NX / 2]) < 1e-8)) throw new AssertionError(); // Consistency check

    // plot: cost-to-go contour J over (X1, X2); x1_opt, x2_opt and (step) u_opt versus time
  }
}

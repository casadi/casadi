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
import java.util.List;
import java.util.Map;

import org.casadi.DM;
import org.casadi.Function;
import org.casadi.GenericMatrix;
import org.casadi.Importer;
import org.casadi.MX;
import org.casadi.Slice;
import org.casadi.ca;

// In this example, we fit a nonlinear model to measurements
//
// This example uses more advanced constructs than the vdp* examples:
// Since the number of control intervals is potentially very large here,
// we use memory-efficient Map and MapAccum, in combination with
// codegeneration.
//
// We will be working with a 2-norm objective:
// || y_measured - y_simulated ||_2^2
//
// This form is well-suited for the Gauss-Newton Hessian approximation.
public class sysid {
  static boolean with_jit;
  static String compiler;

  //########### Create a Gauss-Newton solver ##########
  static Function gauss_newton(MX e, Map<String, MX> nlp, MX V) {
    MX J = ca.jacobian(e, V);
    MX H = ca.triu(ca.mtimes(J.T(), J));
    MX sigma = MX.sym("sigma");
    Function hessLag = new Function("nlp_hess_l", Map.of("x", V, "lam_f", sigma, "hess_gamma_x_x", sigma.times(H)),
        List.of("x", "p", "lam_f", "lam_g"), List.of("hess_gamma_x_x"),
        Map.of("jit", with_jit, "compiler", compiler));
    return ca.nlpsol("solver", "ipopt", nlp, Map.of("hess_lag", hessLag, "jit", with_jit, "compiler", compiler));
  }

  public static void main(String[] args) {
    if (!ca.has_nlpsol("ipopt")) { System.out.println("SKIP: nlpsol ipopt not available"); return; }

    //########### SETTINGS #####################
    int N = 10000;  // Number of samples
    double fs = 610.1; // Sampling frequency [hz]

    DM param_truth = new DM(new double[] {5.625e-6, 2.3e-4, 1, 4.69});
    DM param_guess = new DM(new double[] {5, 2, 1, 5});
    DM scale = new DM(new double[] {1e-6, 1e-4, 1, 1});

    //########### MODELING #####################
    MX y = MX.sym("y");
    MX dy = MX.sym("dy");
    MX u = MX.sym("u");

    MX states = ca.vertcat(List.of(y, dy));
    MX controls = u;

    MX M = MX.sym("M");
    MX c = MX.sym("c");
    MX k = MX.sym("k");
    MX k_NL = MX.sym("k_NL");

    MX params = ca.vertcat(List.of(M, c, k, k_NL));

    MX rhs = ca.vertcat(List.of(dy, u.minus(k_NL.times(y.power(3))).minus(k.times(y)).minus(c.times(dy)).div(M)));

    // Form an ode function
    Function ode = new Function("ode", List.of(states, controls, params), List.of(rhs));

    //########### Creating a simulator ##########
    int N_steps_per_sample = 10;
    double dt = 1 / fs / N_steps_per_sample;

    // Build an integrator for this system: Runge Kutta 4 integrator
    MX k1 = ode.call(List.of(states, controls, params)).get(0);
    MX k2 = ode.call(List.of(states.plus(k1.times(dt / 2.0)), controls, params)).get(0);
    MX k3 = ode.call(List.of(states.plus(k2.times(dt / 2.0)), controls, params)).get(0);
    MX k4 = ode.call(List.of(states.plus(k3.times(dt)), controls, params)).get(0);

    MX states_final = states.plus(k1.plus(k2.times(2)).plus(k3.times(2)).plus(k4).times(dt / 6.0));

    // Create a function that simulates one step propagation in a sample
    Function one_step = new Function("one_step", List.of(states, controls, params), List.of(states_final));

    MX X = states;

    for (int i = 0; i < N_steps_per_sample; i++) {
      X = one_step.call(List.of(X, controls, params)).get(0);
    }

    // Create a function that simulates all step propagation on a sample
    Function one_sample = new Function("one_sample", List.of(states, controls, params), List.of(X));

    //########### Simulating the system ##########
    Function all_samples = one_sample.mapaccum("all_samples", N);

    // Choose an excitation signal
    MT19937 rng = new MT19937(0); // numpy.random.seed(0)
    double[] r = new double[N];
    for (int i = 0; i < N; i++) r[i] = 0.1 * rng.random();
    DM u_data = new DM(r);

    DM x0 = new DM(new double[] {0, 0});
    DM X_measured = (DM) all_samples.invoke(x0, u_data, ca.repmat(param_truth, 1, N));

    DM y_data = X_measured.get(0, new Slice()).T();

    // You may add some noise here
    //y_data = y_data.plus(...0.001*random(N)...);
    // When noise is absent, the fit will be perfect.

    // Use just-in-time compilation to speed up the evaluation
    if (Importer.has_plugin("clang")) {
      with_jit = true;
      compiler = "clang";
    } else if (Importer.has_plugin("shell")) {
      with_jit = true;
      System.out.println("WARNING: on Windows, JIT may require a special environment cfr https://github.com/casadi/casadi/wiki/FAQ:-how-to-perform-jit-for-function-evaluations-of-my-optimization-problem%3F");
      compiler = "shell";
    } else {
      System.out.println("WARNING; running without jit. This may result in very slow evaluation times");
      with_jit = false;
      compiler = "";
    }

    //########### Identifying the simulated system: single shooting strategy ##########

    // Note, it is in general a good idea to scale your decision variables such
    // that they are in the order of ~0.1..100
    MX X_symbolic = (MX) all_samples.invoke(x0, u_data, ca.repmat(params.times(scale), 1, N));

    MX e = ca.minus(y_data, X_symbolic.get(0, new Slice()).T());
    Map<String, MX> nlp = Map.of("x", params, "f", ca.times(0.5, ca.dot(e, e)));

    Function solver = gauss_newton(e, nlp, params);

    Map<String, GenericMatrix> sol = solver.invoke(Map.of("x0", param_guess));

    System.out.println(((DM) sol.get("x")).times(scale));

    if (!(ca.norm_inf(((DM) sol.get("x")).times(scale).minus(param_truth)).to_double() < 1e-8)) throw new AssertionError();

    //########### Identifying the simulated system: multiple shooting strategy ##########

    // All states become decision variables
    X = MX.sym("X", 2, N);

    MX Xn = (MX) one_sample.map(N, "openmp").invoke(X, u_data.T(), ca.repmat(params.times(scale), 1, N));

    MX gaps = Xn.get(new Slice(), new Slice(0, -1)).minus(X.get(new Slice(), new Slice(1, N)));

    e = ca.minus(y_data, Xn.get(0, new Slice()).T());

    MX V = ca.veccat(List.of(params, X));

    nlp = Map.of("x", V, "f", ca.times(0.5, ca.dot(e, e)), "g", ca.vec(gaps));

    // Multipleshooting allows for careful initialization
    DM yd = y_data.get(new Slice(1, N)).minus(y_data.get(new Slice(0, N - 1))).times(fs); // numpy.diff
    DM X_guess = ca.horzcat(List.of(y_data, ca.vertcat(List.of(yd, yd.get(-1))))).T();

    DM x0_ms = ca.veccat(List.of(param_guess, X_guess));

    solver = gauss_newton(e, nlp, V);

    sol = solver.invoke(Map.of("x0", x0_ms, "lbg", 0, "ubg", 0));

    System.out.println(((DM) sol.get("x")).get(new Slice(0, 4)).times(scale));

    if (!(ca.norm_inf(((DM) sol.get("x")).get(new Slice(0, 4)).times(scale).minus(param_truth)).to_double() < 1e-8)) throw new AssertionError();
  }

  // numpy's legacy RandomState (MT19937): seed(int) and random()
  static final class MT19937 {
    private final int[] mt = new int[624];
    private int mti;

    MT19937(long seed) {
      mt[0] = (int) seed;
      for (mti = 1; mti < 624; mti++) mt[mti] = 1812433253 * (mt[mti - 1] ^ (mt[mti - 1] >>> 30)) + mti;
    }

    int next() {
      if (mti >= 624) {
        for (int kk = 0; kk < 624; kk++) {
          int yv = (mt[kk] & 0x80000000) | (mt[(kk + 1) % 624] & 0x7fffffff);
          mt[kk] = mt[(kk + 397) % 624] ^ (yv >>> 1) ^ ((yv & 1) != 0 ? 0x9908b0df : 0);
        }
        mti = 0;
      }
      int yv = mt[mti++];
      yv ^= yv >>> 11;
      yv ^= (yv << 7) & 0x9d2c5680;
      yv ^= (yv << 15) & 0xefc60000;
      yv ^= yv >>> 18;
      return yv;
    }

    double random() {
      long a = (next() >>> 5) & 0xffffffffL, b = (next() >>> 6) & 0xffffffffL;
      return (a * 67108864.0 + b) / 9007199254740992.0;
    }
  }
}

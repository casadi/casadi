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
public static class sysid {
  static bool with_jit;
  static string compiler;

  //########### Create a Gauss-Newton solver ##########
  static Function gauss_newton(MX e, IDictionary<string, MX> nlp, MX V) {
    MX J = ca.jacobian(e, V);
    MX H = ca.triu(ca.mtimes(J.T(), J));
    MX sigma = MX.sym("sigma");
    var hessLag = ca.Function("nlp_hess_l",
        new Dictionary<string, MX> { ["x"] = V, ["lam_f"] = sigma, ["hess_gamma_x_x"] = sigma * H },
        ["x", "p", "lam_f", "lam_g"], ["hess_gamma_x_x"],
        new Dictionary<string, object> { ["jit"] = with_jit, ["compiler"] = compiler });
    return ca.nlpsol("solver", "ipopt", nlp,
        new Dictionary<string, object> { ["hess_lag"] = hessLag, ["jit"] = with_jit, ["compiler"] = compiler });
  }

  public static void Main() {
    if (!ca.has_nlpsol("ipopt")) { Console.WriteLine("SKIP: nlpsol ipopt not available"); return; }

    //########### SETTINGS #####################
    int N = 10000;  // Number of samples
    double fs = 610.1; // Sampling frequency [hz]

    DM param_truth = new double[] { 5.625e-6, 2.3e-4, 1, 4.69 };
    DM param_guess = new double[] { 5, 2, 1, 5 };
    DM scale = new double[] { 1e-6, 1e-4, 1, 1 };

    //########### MODELING #####################
    MX y = MX.sym("y");
    MX dy = MX.sym("dy");
    MX u = MX.sym("u");

    MX states = ca.vertcat(y, dy);
    MX controls = u;

    MX M = MX.sym("M");
    MX c = MX.sym("c");
    MX k = MX.sym("k");
    MX k_NL = MX.sym("k_NL");

    MX @params = ca.vertcat(M, c, k, k_NL);

    MX rhs = ca.vertcat(dy, (u - k_NL * ca.pow(y, 3) - k * y - c * dy) / M);

    // Form an ode function
    var ode = ca.Function("ode", [states, controls, @params], [rhs]);

    //########### Creating a simulator ##########
    int N_steps_per_sample = 10;
    double dt = 1 / fs / N_steps_per_sample;

    // Build an integrator for this system: Runge Kutta 4 integrator
    MX k1 = ode.Call(states, controls, @params);
    MX k2 = ode.Call(states + dt / 2.0 * k1, controls, @params);
    MX k3 = ode.Call(states + dt / 2.0 * k2, controls, @params);
    MX k4 = ode.Call(states + dt * k3, controls, @params);

    MX states_final = states + dt / 6.0 * (k1 + 2 * k2 + 2 * k3 + k4);

    // Create a function that simulates one step propagation in a sample
    var one_step = ca.Function("one_step", [states, controls, @params], [states_final]);

    MX X = states;

    for (int i = 0; i < N_steps_per_sample; i++) {
      X = one_step.Call(X, controls, @params);
    }

    // Create a function that simulates all step propagation on a sample
    var one_sample = ca.Function("one_sample", [states, controls, @params], [X]);

    //########### Simulating the system ##########
    Function all_samples = one_sample.mapaccum("all_samples", N);

    // Choose an excitation signal
    var rng = new MT19937(0); // numpy.random.seed(0)
    var r = new double[N];
    for (int i = 0; i < N; i++) r[i] = 0.1 * rng.random();
    DM u_data = r;

    DM x0 = new double[] { 0, 0 };
    DM X_measured = all_samples.Call(x0, u_data, ca.repmat(param_truth, 1, N));

    DM y_data = X_measured[0, ..].T();

    // You may add some noise here
    //y_data = y_data + 0.001*random(N);
    // When noise is absent, the fit will be perfect.

    // Use just-in-time compilation to speed up the evaluation
    if (Importer.has_plugin("clang")) {
      with_jit = true;
      compiler = "clang";
    } else if (Importer.has_plugin("shell")) {
      with_jit = true;
      Console.WriteLine("WARNING: on Windows, JIT may require a special environment cfr https://github.com/casadi/casadi/wiki/FAQ:-how-to-perform-jit-for-function-evaluations-of-my-optimization-problem%3F");
      compiler = "shell";
    } else {
      Console.WriteLine("WARNING; running without jit. This may result in very slow evaluation times");
      with_jit = false;
      compiler = "";
    }

    //########### Identifying the simulated system: single shooting strategy ##########

    // Note, it is in general a good idea to scale your decision variables such
    // that they are in the order of ~0.1..100
    MX X_symbolic = all_samples.Call(x0, u_data, ca.repmat(@params * scale, 1, N));

    MX e = y_data - X_symbolic[0, ..].T();
    var nlp = new Dictionary<string, MX> { ["x"] = @params, ["f"] = 0.5 * ca.dot(e, e) };

    Function solver = gauss_newton(e, nlp, @params);

    var sol = solver.call(new Dictionary<string, DM> { ["x0"] = param_guess });

    Console.WriteLine(sol["x"] * scale);

    if (!((double)ca.norm_inf(sol["x"] * scale - param_truth) < 1e-8)) throw new Exception("assertion failed");

    //########### Identifying the simulated system: multiple shooting strategy ##########

    // All states become decision variables
    X = MX.sym("X", 2, N);

    MX Xn = one_sample.map(N, "openmp").Call(X, u_data.T(), ca.repmat(@params * scale, 1, N));

    MX gaps = Xn[.., 0..^1] - X[.., 1..N];

    e = y_data - Xn[0, ..].T();

    MX V = ca.veccat(@params, X);

    nlp = new Dictionary<string, MX> { ["x"] = V, ["f"] = 0.5 * ca.dot(e, e), ["g"] = ca.vec(gaps) };

    // Multipleshooting allows for careful initialization
    DM yd = (y_data[1..N] - y_data[0..(N - 1)]) * fs; // numpy.diff
    DM X_guess = ca.horzcat(y_data, ca.vertcat(yd, yd[^1])).T();

    DM x0_ms = ca.veccat(param_guess, X_guess);

    solver = gauss_newton(e, nlp, V);

    sol = solver.call(new Dictionary<string, DM> { ["x0"] = x0_ms, ["lbg"] = 0, ["ubg"] = 0 });

    Console.WriteLine(sol["x"][0..4] * scale);

    if (!((double)ca.norm_inf(sol["x"][0..4] * scale - param_truth) < 1e-8)) throw new Exception("assertion failed");
  }

  // numpy's legacy RandomState (MT19937): seed(int) and random()
  sealed class MT19937 {
    private readonly uint[] mt = new uint[624];
    private int mti;

    public MT19937(uint seed) {
      mt[0] = seed;
      for (mti = 1; mti < 624; mti++) mt[mti] = 1812433253u * (mt[mti - 1] ^ (mt[mti - 1] >> 30)) + (uint)mti;
    }

    uint next() {
      if (mti >= 624) {
        for (int kk = 0; kk < 624; kk++) {
          uint yv = (mt[kk] & 0x80000000u) | (mt[(kk + 1) % 624] & 0x7fffffffu);
          mt[kk] = mt[(kk + 397) % 624] ^ (yv >> 1) ^ ((yv & 1) != 0 ? 0x9908b0dfu : 0u);
        }
        mti = 0;
      }
      uint y = mt[mti++];
      y ^= y >> 11;
      y ^= (y << 7) & 0x9d2c5680u;
      y ^= (y << 15) & 0xefc60000u;
      y ^= y >> 18;
      return y;
    }

    public double random() {
      ulong a = next() >> 5, b = next() >> 6;
      return (a * 67108864.0 + b) / 9007199254740992.0;
    }
  }
}

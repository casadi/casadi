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
// For documentation about this examples, check http://docs.casadi.org/documents/mhe_spring_damper.pdf

using System;
using System.Collections.Generic;
using CasADi;

public static class mhe_spring_damper {
  // Settings of the filter
  const int N = 10; // Horizon length
  const double dt = 0.05; // Time step

  // casadi.tools structures are emulated by vertcat-ing the entries in the same order:
  // shooting = [X_0, W_0, X_1, W_1, ..., X_{N-1}], parameters = [U_0, Y_0, ..., Y_{N-1}, vec(S), x0]
  static List<T> shootingList<T>(T[] X, T[] W) {
    var e = new List<T>();
    for (int i = 0; i < N; i++) {
      e.Add(X[i]);
      if (i < N - 1) e.Add(W[i]);
    }
    return e;
  }
  static SX shootingVec(SX[] X, SX[] W) => ca.vertcat(shootingList(X, W));
  static DM shootingVec(DM[] X, DM[] W) => ca.vertcat(shootingList(X, W));

  static List<T> parametersList<T>(T[] U, T[] Y, T vecS, T x0) {
    var e = new List<T>();
    for (int i = 0; i < N; i++) {
      if (i < N - 1) e.Add(U[i]);
      e.Add(Y[i]);
    }
    e.Add(vecS);
    e.Add(x0);
    return e;
  }
  static SX parametersVec(SX[] U, SX[] Y, SX vecS, SX x0) => ca.vertcat(parametersList(U, Y, vecS, x0));
  static DM parametersVec(DM[] U, DM[] Y, DM vecS, DM x0) => ca.vertcat(parametersList(U, Y, vecS, x0));

  // solution["X",i] and solution["W",i]
  static DM solX(DM v, int i) => v[(3 * i)..(3 * i + 2)];
  static DM solW(DM v, int i) => v[3 * i + 2];

  // Standard normal samples (Box-Muller)
  static DM randn(Random rng, int n, int m) {
    DM r = DM.zeros(n, m);
    for (int j = 0; j < m; j++)
      for (int i = 0; i < n; i++)
        r[i, j] = Math.Sqrt(-2 * Math.Log(1 - rng.NextDouble())) * Math.Cos(2 * Math.PI * rng.NextDouble());
    return r;
  }

  public static void Main() {
    if (!ca.has_nlpsol("ipopt")) { Console.WriteLine("SKIP: nlpsol ipopt not available"); return; }
    var rng = new Random(0);

    double sigma_p = 0.005; // Standard deviation of the position measurements
    double sigma_w = 0.1; // Standard deviation for the process noise
    DM R = 1 / (sigma_p * sigma_p); // resulting weighting matrix for the position measurements
    DM Q = 1 / (sigma_w * sigma_w); // resulting weighting matrix for the process noise

    int Nsimulation = 1000; // Lenght of the simulation

    // Parameters of the system
    double m = 1; // The weight of the mass
    double k = 1; // The spring constant
    double c = 0.5; // The damping of the system
    // The state
    SX x = SX.sym("x"), dx = SX.sym("dx");
    SX states = ca.vertcat(x, dx); // Full state vector of the system: position x and velocity dx
    int Nstates = (int)states.numel(); // Number of states

    // The control input
    SX F = SX.sym("F");
    SX controls = F; // Full control vector of the system: Input force F
    int Ncontrols = (int)controls.numel(); // Number of control inputs

    // Disturbances
    SX w = SX.sym("w");
    SX disturbances = w; // Process noise vector
    int Ndisturbances = (int)disturbances.numel(); // Number of disturbances

    // Measurements
    SX y = SX.sym("y");
    SX measurements = y; // Measurement vector
    int Nmeas = (int)measurements.numel(); // Number of measurements

    // Create Structure for the entire horizon

    // Structure that will be degrees of freedom for the optimizer
    SX[] Xs = new SX[N], Ws = new SX[N - 1];
    for (int i = 0; i < N; i++) Xs[i] = SX.sym($"X_{i}", Nstates);
    for (int i = 0; i < N - 1; i++) Ws[i] = SX.sym($"W_{i}", Ndisturbances);
    SX shooting = shootingVec(Xs, Ws);
    // Structure that will be fixed parameters for the optimizer
    SX[] Us = new SX[N - 1], Ys = new SX[N];
    for (int i = 0; i < N - 1; i++) Us[i] = SX.sym($"U_{i}", Ncontrols);
    for (int i = 0; i < N; i++) Ys[i] = SX.sym($"Y_{i}", Nmeas);
    SX S = SX.sym("S", Nstates, Nstates);
    SX x0s = SX.sym("x0", Nstates, 1);
    SX parameters = parametersVec(Us, Ys, ca.vec(S), x0s);
    // Define the ODE right hand side
    SX rhs = ca.vertcat(dx, (-k * x - c * dx + F) / m + w);

    var f = ca.Function("f", [states, controls, disturbances], [rhs]);

    // Build an integrator for this system: Runge Kutta 4 integrator
    SX k1 = f.Call(states, controls, disturbances);
    SX k2 = f.Call(states + dt / 2.0 * k1, controls, disturbances);
    SX k3 = f.Call(states + dt / 2.0 * k2, controls, disturbances);
    SX k4 = f.Call(states + dt * k3, controls, disturbances);

    SX states_1 = states + dt / 6.0 * (k1 + 2 * k2 + 2 * k3 + k4);
    var phi = ca.Function("phi", [states, controls, disturbances], [states_1],
        ["x", "u", "d"], ["x1"]);
    Function PHI = phi.factory("PHI", ["x", "u", "d"], ["jac:x1:x"]);
    // Define the measurement system
    var h = ca.Function("h", [states], [x], ["x"], ["y"]); // We have measurements of the position
    Function H = h.factory("H", ["x"], ["jac:y:x"]);
    // Build the objective
    SX obj = 0;
    // First the arrival cost
    obj += ca.mtimes([(Xs[0] - x0s).T(), S, Xs[0] - x0s]);
    // Next the cost for the measurement noise
    for (int i = 0; i < N; i++) {
      SX vm = h.Call(Xs[i]) - Ys[i];
      obj += ca.mtimes([vm.T(), R, vm]);
    }
    // And also the cost for the process noise
    for (int i = 0; i < N - 1; i++) {
      obj += ca.mtimes([Ws[i].T(), Q, Ws[i]]);
    }

    // Build the multiple shooting constraints
    var g = new List<SX>();
    for (int i = 0; i < N - 1; i++) {
      g.Add(Xs[i + 1] - phi.Call(Xs[i], Us[i], Ws[i]));
    }

    // Formulate the NLP
    var nlp = new Dictionary<string, SX> { ["x"] = shooting, ["p"] = parameters, ["f"] = obj, ["g"] = ca.vertcat(g) };

    // Make a simulation to create the data for the problem
    DM simulated_X = DM.zeros(Nstates, Nsimulation);
    simulated_X[.., 0] = new double[] { 1, 0 }; // Initial state
    DM t = ca.linspace(0, (Nsimulation - 1) * dt, Nsimulation); // Time grid
    DM simulated_U = ca.cos(t[0..(Nsimulation - 1)]).T(); // control input for the simulation
    simulated_U[.., (Nsimulation / 2)..(Nsimulation - 1)] = 0.0;
    DM simulated_W = sigma_w * randn(rng, Ndisturbances, Nsimulation - 1); // Process noise for the simulation
    for (int i = 0; i < Nsimulation - 1; i++) {
      simulated_X[.., i + 1] = phi.Call(simulated_X[.., i], simulated_U[.., i], simulated_W[.., i]);
    }
    // Create the measurements from these states
    DM simulated_Y = DM.zeros(Nmeas, Nsimulation); // Holder for the measurements
    for (int i = 0; i < Nsimulation; i++) {
      simulated_Y[.., i] = h.Call(simulated_X[.., i]);
    }
    // Add noise the the position measurements
    simulated_Y += sigma_p * randn(rng, (int)simulated_Y.size1(), (int)simulated_Y.size2());

    // The initial estimate and related covariance, which will be used for the arrival cost
    double sigma_x0 = 0.01;
    DM P = sigma_x0 * sigma_x0 * DM.eye(Nstates);
    DM x0 = simulated_X[.., 0] + sigma_x0 * randn(rng, Nstates, 1);
    // Create the solver
    var opts = new Dictionary<string, object> { ["ipopt.print_level"] = 0L, ["print_time"] = false, ["ipopt.max_iter"] = 100L };
    Function nlpsol = ca.nlpsol("nlpsol", "ipopt", nlp, opts);

    // Create a holder for the estimated states and disturbances
    DM estimated_X = DM.zeros(Nstates, Nsimulation);
    DM estimated_W = DM.zeros(Ndisturbances, Nsimulation - 1);
    // For the first instance we run the filter, we need to initialize it.
    DM[] cur_U = new DM[N - 1], cur_Y = new DM[N];
    for (int j = 0; j < N - 1; j++) cur_U[j] = simulated_U[.., j];
    for (int j = 0; j < N; j++) cur_Y[j] = simulated_Y[.., j];
    DM cur_S = ca.inv(P); // Arrival cost is the inverse of the initial covariance
    DM cur_x0 = x0;
    DM[] init_X = new DM[N], init_W = new DM[N - 1];
    for (int j = 0; j < N; j++) init_X[j] = simulated_X[.., j];
    for (int j = 0; j < N - 1; j++) init_W[j] = DM.zeros(Ndisturbances, 1);
    var res = nlpsol.call(new Dictionary<string, DM> {
      ["p"] = parametersVec(cur_U, cur_Y, ca.vec(cur_S), cur_x0),
      ["x0"] = shootingVec(init_X, init_W), ["lbg"] = 0, ["ubg"] = 0 });

    // Get the solution
    DM solution = res["x"];
    for (int j = 0; j < N; j++) estimated_X[.., j] = solX(solution, j);
    for (int j = 0; j < N - 1; j++) estimated_W[.., j] = solW(solution, j);

    // Now make a loop for the rest of the simulation
    for (int i = 1; i < Nsimulation - N + 1; i++) {

      // Update the arrival cost, using linearisations around the estimate of MHE at the beginning of the horizon (according to the 'Smoothed EKF Update'): first update the state and covariance with the measurement that will be deleted, and next propagate the state and covariance because of the shifting of the horizon
      Console.WriteLine($"step {i}/{Nsimulation - N} ({nlpsol.stats()["return_status"]})");
      DM H0 = H.Call(solX(solution, 0));
      DM K = ca.mtimes(new[] { P, H0.T(), ca.inv(ca.mtimes(new[] { H0, P, H0.T() }) + R) });
      P = ca.mtimes(DM.eye(Nstates) - ca.mtimes(K, H0), P);
      DM h0 = h.Call(solX(solution, 0));
      x0 += ca.mtimes(K, cur_Y[0] - h0 - ca.mtimes(H0, x0 - solX(solution, 0)));
      x0 = phi.Call(x0, cur_U[0], solW(solution, 0));
      DM Fk = PHI.Call(solX(solution, 0), cur_U[0], solW(solution, 0));
      P = ca.mtimes(new[] { Fk, P, Fk.T() }) + ca.inv(Q);
      // Get the measurements and control inputs
      for (int j = 0; j < N - 1; j++) cur_U[j] = simulated_U[.., i + j];
      for (int j = 0; j < N; j++) cur_Y[j] = simulated_Y[.., i + j];
      cur_S = ca.inv(P);
      cur_x0 = x0;
      // Initialize the system with the shifted solution
      for (int j = 0; j < N - 2; j++) init_W[j] = estimated_W[.., i + j]; // The shifted solution for the disturbances
      init_W[N - 2] = DM.zeros(Ndisturbances, 1); // The last node for the disturbances is initialized with zeros
      for (int j = 0; j < N - 1; j++) init_X[j] = estimated_X[.., i + j]; // The shifted solution for the state estimates
      // The last node for the state is initialized with a forward simulation
      DM phi0 = phi.Call(init_X[N - 1], cur_U[N - 2], init_W[N - 2]);
      init_X[N - 1] = phi0;
      // And now initialize the solver and solve the problem
      res = nlpsol.call(new Dictionary<string, DM> {
        ["p"] = parametersVec(cur_U, cur_Y, ca.vec(cur_S), cur_x0),
        ["x0"] = shootingVec(init_X, init_W), ["lbg"] = 0, ["ubg"] = 0 });
      solution = res["x"];

      // Now get the state estimate. Note that we are only interested in the last node of the horizon
      estimated_X[.., N - 1 + i] = solX(solution, N - 1);
      estimated_W[.., N - 2 + i] = solW(solution, N - 2);
    }
    // plot: estimated vs real position, and the position error, over t

    DM error = estimated_X[0, ..] - simulated_X[0, ..];
    Console.WriteLine(ca.mtimes(error, error.T()));
    if (!((double)ca.mtimes(error, error.T()) < 0.01)) throw new Exception("assertion failed");
  }
}

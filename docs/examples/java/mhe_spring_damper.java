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

import org.casadi.*;
import java.util.*;

public class mhe_spring_damper {
  // Settings of the filter
  static final int N = 10; // Horizon length
  static final double dt = 0.05; // Time step

  // casadi.tools structures are emulated by vertcat-ing the entries in the same order:
  // shooting = [X_0, W_0, X_1, W_1, ..., X_{N-1}], parameters = [U_0, Y_0, ..., Y_{N-1}, vec(S), x0]
  static <T extends GenericMatrix> T shootingVec(T[] X, T[] W) {
    List<T> e = new ArrayList<>();
    for (int i = 0; i < N; i++) {
      e.add(X[i]);
      if (i < N - 1) e.add(W[i]);
    }
    return ca.vertcat(e);
  }

  static <T extends GenericMatrix> T parametersVec(T[] U, T[] Y, T vecS, T x0) {
    List<T> e = new ArrayList<>();
    for (int i = 0; i < N; i++) {
      if (i < N - 1) e.add(U[i]);
      e.add(Y[i]);
    }
    e.add(vecS);
    e.add(x0);
    return ca.vertcat(e);
  }

  // solution["X",i] and solution["W",i]
  static DM solX(DM v, int i) { return v.get(new Slice(3 * i, 3 * i + 2)); }
  static DM solW(DM v, int i) { return v.get(3 * i + 2); }

  static DM randn(Random rng, int n, int m) {
    DM r = DM.zeros(n, m);
    for (int j = 0; j < m; j++)
      for (int i = 0; i < n; i++) r.set(rng.nextGaussian(), false, i, j);
    return r;
  }

  public static void main(String[] args) {
    if (!ca.has_nlpsol("ipopt")) { System.out.println("SKIP: nlpsol ipopt not available"); return; }
    Random rng = new Random(0);

    double sigma_p = 0.005; // Standard deviation of the position measurements
    double sigma_w = 0.1; // Standard deviation for the process noise
    DM R = new DM(1 / (sigma_p * sigma_p)); // resulting weighting matrix for the position measurements
    DM Q = new DM(1 / (sigma_w * sigma_w)); // resulting weighting matrix for the process noise

    int Nsimulation = 1000; // Lenght of the simulation

    // Parameters of the system
    double m = 1; // The weight of the mass
    double k = 1; // The spring constant
    double c = 0.5; // The damping of the system
    // The state
    SX x = SX.sym("x"), dx = SX.sym("dx");
    SX states = ca.vertcat(List.of(x, dx)); // Full state vector of the system: position x and velocity dx
    int Nstates = (int) states.numel(); // Number of states

    // The control input
    SX F = SX.sym("F");
    SX controls = F; // Full control vector of the system: Input force F
    int Ncontrols = (int) controls.numel(); // Number of control inputs

    // Disturbances
    SX w = SX.sym("w");
    SX disturbances = w; // Process noise vector
    int Ndisturbances = (int) disturbances.numel(); // Number of disturbances

    // Measurements
    SX y = SX.sym("y");
    SX measurements = y; // Measurement vector
    int Nmeas = (int) measurements.numel(); // Number of measurements

    // Create Structure for the entire horizon

    // Structure that will be degrees of freedom for the optimizer
    SX[] Xs = new SX[N], Ws = new SX[N - 1];
    for (int i = 0; i < N; i++) Xs[i] = SX.sym("X_" + i, Nstates);
    for (int i = 0; i < N - 1; i++) Ws[i] = SX.sym("W_" + i, Ndisturbances);
    SX shooting = shootingVec(Xs, Ws);
    // Structure that will be fixed parameters for the optimizer
    SX[] Us = new SX[N - 1], Ys = new SX[N];
    for (int i = 0; i < N - 1; i++) Us[i] = SX.sym("U_" + i, Ncontrols);
    for (int i = 0; i < N; i++) Ys[i] = SX.sym("Y_" + i, Nmeas);
    SX S = SX.sym("S", Nstates, Nstates);
    SX x0s = SX.sym("x0", Nstates, 1);
    SX parameters = parametersVec(Us, Ys, ca.vec(S), x0s);
    // Define the ODE right hand side
    SX rhs = ca.vertcat(List.of(dx, x.times(-k).minus(dx.times(c)).plus(F).div(m).plus(w)));

    Function f = new Function("f", List.of(states, controls, disturbances), List.of(rhs));

    // Build an integrator for this system: Runge Kutta 4 integrator
    SX k1 = (SX) f.invoke(states, controls, disturbances);
    SX k2 = (SX) f.invoke(states.plus(k1.times(dt / 2.0)), controls, disturbances);
    SX k3 = (SX) f.invoke(states.plus(k2.times(dt / 2.0)), controls, disturbances);
    SX k4 = (SX) f.invoke(states.plus(k3.times(dt)), controls, disturbances);

    SX states_1 = states.plus(k1.plus(k2.times(2)).plus(k3.times(2)).plus(k4).times(dt / 6.0));
    Function phi = new Function("phi", List.of(states, controls, disturbances), List.of(states_1),
        List.of("x", "u", "d"), List.of("x1"));
    Function PHI = phi.factory("PHI", List.of("x", "u", "d"), List.of("jac:x1:x"));
    // Define the measurement system
    Function h = new Function("h", List.of(states), List.of(x), List.of("x"), List.of("y")); // We have measurements of the position
    Function H = h.factory("H", List.of("x"), List.of("jac:y:x"));
    // Build the objective
    SX obj = new SX(0);
    // First the arrival cost
    obj = obj.plus(ca.mtimes(List.of(Xs[0].minus(x0s).T(), S, Xs[0].minus(x0s))));
    // Next the cost for the measurement noise
    for (int i = 0; i < N; i++) {
      SX vm = ((SX) h.invoke(Xs[i])).minus(Ys[i]);
      obj = obj.plus(ca.mtimes(List.of(vm.T(), new SX(R), vm)));
    }
    // And also the cost for the process noise
    for (int i = 0; i < N - 1; i++) {
      obj = obj.plus(ca.mtimes(List.of(Ws[i].T(), new SX(Q), Ws[i])));
    }

    // Build the multiple shooting constraints
    List<SX> g = new ArrayList<>();
    for (int i = 0; i < N - 1; i++) {
      g.add(Xs[i + 1].minus((SX) phi.invoke(Xs[i], Us[i], Ws[i])));
    }

    // Formulate the NLP
    Map<String, SX> nlp = Map.of("x", shooting, "p", parameters, "f", obj, "g", ca.vertcat(g));

    // Make a simulation to create the data for the problem
    DM simulated_X = DM.zeros(Nstates, Nsimulation);
    simulated_X.set(new Slice(), 0, new DM(new double[] {1, 0})); // Initial state
    DM t = ca.linspace(new DM(0), new DM((Nsimulation - 1) * dt), Nsimulation); // Time grid
    DM simulated_U = ca.cos(t.get(new Slice(0, Nsimulation - 1))).T(); // control input for the simulation
    simulated_U.set(new Slice(), new Slice(Nsimulation / 2, Nsimulation - 1), new DM(0.0));
    DM simulated_W = randn(rng, Ndisturbances, Nsimulation - 1).times(sigma_w); // Process noise for the simulation
    for (int i = 0; i < Nsimulation - 1; i++) {
      simulated_X.set(new Slice(), i + 1, (DM) phi.invoke(simulated_X.get(new Slice(), i),
          simulated_U.get(new Slice(), i), simulated_W.get(new Slice(), i)));
    }
    // Create the measurements from these states
    DM simulated_Y = DM.zeros(Nmeas, Nsimulation); // Holder for the measurements
    for (int i = 0; i < Nsimulation; i++) {
      simulated_Y.set(new Slice(), i, (DM) h.invoke(simulated_X.get(new Slice(), i)));
    }
    // Add noise the the position measurements
    simulated_Y = simulated_Y.plus(randn(rng, (int) simulated_Y.size1(), (int) simulated_Y.size2()).times(sigma_p));

    // The initial estimate and related covariance, which will be used for the arrival cost
    double sigma_x0 = 0.01;
    DM P = DM.eye(Nstates).times(sigma_x0 * sigma_x0);
    DM x0 = simulated_X.get(new Slice(), 0).plus(randn(rng, Nstates, 1).times(sigma_x0));
    // Create the solver
    Map<String, Object> opts = Map.of("ipopt.print_level", 0, "print_time", false, "ipopt.max_iter", 100);
    Function nlpsol = ca.nlpsol("nlpsol", "ipopt", nlp, opts);

    // Create a holder for the estimated states and disturbances
    DM estimated_X = DM.zeros(Nstates, Nsimulation);
    DM estimated_W = DM.zeros(Ndisturbances, Nsimulation - 1);
    // For the first instance we run the filter, we need to initialize it.
    DM[] cur_U = new DM[N - 1], cur_Y = new DM[N];
    for (int j = 0; j < N - 1; j++) cur_U[j] = simulated_U.get(new Slice(), j);
    for (int j = 0; j < N; j++) cur_Y[j] = simulated_Y.get(new Slice(), j);
    DM cur_S = ca.inv(P); // Arrival cost is the inverse of the initial covariance
    DM cur_x0 = x0;
    DM[] init_X = new DM[N], init_W = new DM[N - 1];
    for (int j = 0; j < N; j++) init_X[j] = simulated_X.get(new Slice(), j);
    for (int j = 0; j < N - 1; j++) init_W[j] = DM.zeros(Ndisturbances, 1);
    Map<String, GenericMatrix> res = nlpsol.invoke(Map.of("p", parametersVec(cur_U, cur_Y, ca.vec(cur_S), cur_x0),
        "x0", shootingVec(init_X, init_W), "lbg", 0, "ubg", 0));

    // Get the solution
    DM solution = (DM) res.get("x");
    for (int j = 0; j < N; j++) estimated_X.set(new Slice(), j, solX(solution, j));
    for (int j = 0; j < N - 1; j++) estimated_W.set(new Slice(), j, solW(solution, j));

    // Now make a loop for the rest of the simulation
    for (int i = 1; i < Nsimulation - N + 1; i++) {

      // Update the arrival cost, using linearisations around the estimate of MHE at the beginning of the horizon (according to the 'Smoothed EKF Update'): first update the state and covariance with the measurement that will be deleted, and next propagate the state and covariance because of the shifting of the horizon
      System.out.println(String.format("step %d/%d (%s)", i, Nsimulation - N, nlpsol.stats().get("return_status")));
      DM H0 = (DM) H.invoke(solX(solution, 0));
      DM K = ca.mtimes(List.of(P, H0.T(), ca.inv(ca.mtimes(List.of(H0, P, H0.T())).plus(R))));
      P = ca.mtimes(DM.eye(Nstates).minus(ca.mtimes(K, H0)), P);
      DM h0 = (DM) h.invoke(solX(solution, 0));
      x0 = x0.plus(ca.mtimes(K, cur_Y[0].minus(h0).minus(ca.mtimes(H0, x0.minus(solX(solution, 0))))));
      x0 = (DM) phi.invoke(x0, cur_U[0], solW(solution, 0));
      DM Fk = (DM) PHI.invoke(solX(solution, 0), cur_U[0], solW(solution, 0));
      P = ca.mtimes(List.of(Fk, P, Fk.T())).plus(ca.inv(Q));
      // Get the measurements and control inputs
      for (int j = 0; j < N - 1; j++) cur_U[j] = simulated_U.get(new Slice(), i + j);
      for (int j = 0; j < N; j++) cur_Y[j] = simulated_Y.get(new Slice(), i + j);
      cur_S = ca.inv(P);
      cur_x0 = x0;
      // Initialize the system with the shifted solution
      for (int j = 0; j < N - 2; j++) init_W[j] = estimated_W.get(new Slice(), i + j); // The shifted solution for the disturbances
      init_W[N - 2] = DM.zeros(Ndisturbances, 1); // The last node for the disturbances is initialized with zeros
      for (int j = 0; j < N - 1; j++) init_X[j] = estimated_X.get(new Slice(), i + j); // The shifted solution for the state estimates
      // The last node for the state is initialized with a forward simulation
      DM phi0 = (DM) phi.invoke(init_X[N - 1], cur_U[N - 2], init_W[N - 2]);
      init_X[N - 1] = phi0;
      // And now initialize the solver and solve the problem
      res = nlpsol.invoke(Map.of("p", parametersVec(cur_U, cur_Y, ca.vec(cur_S), cur_x0),
          "x0", shootingVec(init_X, init_W), "lbg", 0, "ubg", 0));
      solution = (DM) res.get("x");

      // Now get the state estimate. Note that we are only interested in the last node of the horizon
      estimated_X.set(new Slice(), N - 1 + i, solX(solution, N - 1));
      estimated_W.set(new Slice(), N - 2 + i, solW(solution, N - 2));
    }
    // plot: estimated vs real position, and the position error, over t

    DM error = estimated_X.get(0, new Slice()).minus(simulated_X.get(0, new Slice()));
    System.out.println(ca.mtimes(error, error.T()));
    if (!(ca.mtimes(error, error.T()).to_double() < 0.01)) throw new AssertionError();
  }
}

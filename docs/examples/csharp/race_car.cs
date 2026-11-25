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
using CasADi;

// Car race along a track
// ----------------------
// An optimal control problem (OCP),
// solved with direct multiple-shooting.
//
// For more information see: http://labs.casadi.org/OCP

public static class race_car {
  public static void Main() {
    if (!ca.has_nlpsol("ipopt")) { Console.WriteLine("SKIP: nlpsol ipopt not available"); return; }

    int N = 100; // number of control intervals

    var opti = ca.Opti(); // Optimization problem

    // ---- decision variables ---------
    MX X = opti.variable(2, N + 1); // state trajectory
    MX pos   = X[0, ..];
    MX speed = X[1, ..];
    MX U = opti.variable(1, N);   // control trajectory (throttle)
    MX T = opti.variable();       // final time

    // ---- objective          ---------
    opti.minimize(T); // race in minimal time

    // ---- dynamic constraints --------
    Func<MX, MX, MX> f = (x, u) => ca.vertcat(x[1], u - x[1]); // dx/dt = f(x,u)

    MX dt = T / N; // length of a control interval
    for (int k = 0; k < N; k++) { // loop over control intervals
      // Runge-Kutta 4 integration
      MX Xk = X[.., k], Uk = U[.., k];
      MX k1 = f(Xk,               Uk);
      MX k2 = f(Xk + dt / 2 * k1, Uk);
      MX k3 = f(Xk + dt / 2 * k2, Uk);
      MX k4 = f(Xk + dt * k3,     Uk);
      MX x_next = Xk + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4);
      opti.subject_to(ca.eq(X[.., k + 1], x_next)); // close the gaps
    }

    // ---- path constraints -----------
    opti.subject_to(speed <= limit(pos));   // track speed limit
    opti.subject_to(opti.bounded(0, U, 1)); // control is limited

    // ---- boundary conditions --------
    opti.subject_to(ca.eq(pos[0], 0));   // start at position 0 ...
    opti.subject_to(ca.eq(speed[0], 0)); // ... from stand-still
    opti.subject_to(ca.eq(pos[^1], 1));  // finish line at position 1

    // ---- misc. constraints  ----------
    opti.subject_to(T >= 0); // Time must be positive

    // ---- initial values for solver ---
    opti.set_initial(speed, 1);
    opti.set_initial(T, 1);

    // ---- solve NLP              ------
    opti.solver("ipopt"); // set numerical backend
    OptiSol sol = opti.solve();   // actual solve

    // ---- post-processing        ------
    DM speed_opt = sol.value(speed);
    DM pos_opt = sol.value(pos);
    DM limit_opt = limit(pos_opt);
    DM U_opt = sol.value(U);
    // plot: speed_opt, pos_opt, limit_opt (speed limit) and step plot of U_opt (throttle)

    DM jac_g = sol.value(ca.jacobian(opti.g, opti.x));
    // plot: spy(jac_g)
    DM hess_l = sol.value(ca.hessian(opti.f + ca.dot(opti.lam_g, opti.g), opti.x).Item1);
    // plot: spy(hess_l)
  }

  // track speed limit (Python uses one lambda for both MX and DM)
  static MX limit(MX pos) => 1 - ca.sin(2 * ca.pi * pos) / 2;
  static DM limit(DM pos) => 1 - ca.sin(2 * ca.pi * pos) / 2;
}

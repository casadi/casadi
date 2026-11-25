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

// casadi.tools structures (struct_symSX, struct_SX) are Python-only: plain SX/DM vectors replace them here.
public static class simulation {
  public static void Main() {
    // System states
    SX x = SX.sym("x"), y = SX.sym("y"), dx = SX.sym("dx"), dy = SX.sym("dy");
    SX states = ca.vertcat(x, y, dx, dy);

    // System controls
    SX u = SX.sym("u"), v = SX.sym("v");
    SX controls = ca.vertcat(u, v);

    // System parameters
    SX k = SX.sym("k"), c = SX.sym("c"), beta = SX.sym("beta");
    SX parameters = ca.vertcat(k, c, beta);

    // Provide some numerical values (order: k, c, beta)
    DM parameters_ = DM.zeros(3);
    parameters_[0] = 10.0; // k
    parameters_[2] = 1.0;  // beta
    parameters_[1] = 1.0;  // c

    SX vel = ca.vertcat(dx, dy);
    SX p = ca.vertcat(x, y);
    SX q = ca.vertcat(u, v);

    // System dynamics
    SX F = -k * (p - q) - beta * v * ca.sqrt(ca.sumsqr(vel) + ca.pow(c, 2));

    // System right hand side
    SX rhs = ca.vertcat(dx, dy, F[0], F[1]);

    // The remainder of the Python example is commented-out legacy code (ControlSimulator, covariance
    // propagation by integrator sensitivities, Lyapunov equations and unscented propagation).
  }
}

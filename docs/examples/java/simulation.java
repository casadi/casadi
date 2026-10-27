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

import org.casadi.DM;
import org.casadi.SX;
import org.casadi.ca;

// casadi.tools structures (struct_symSX, struct_SX) are Python-only: plain SX/DM vectors replace them here.
public class simulation {
  public static void main(String[] args) {
    // System states
    SX x = SX.sym("x"), y = SX.sym("y"), dx = SX.sym("dx"), dy = SX.sym("dy");
    SX states = ca.vertcat(List.of(x, y, dx, dy));

    // System controls
    SX u = SX.sym("u"), v = SX.sym("v");
    SX controls = ca.vertcat(List.of(u, v));

    // System parameters
    SX k = SX.sym("k"), c = SX.sym("c"), beta = SX.sym("beta");
    SX parameters = ca.vertcat(List.of(k, c, beta));

    // Provide some numerical values (order: k, c, beta)
    DM parameters_ = DM.zeros(3);
    parameters_.set(0, 10.0); // k
    parameters_.set(2, 1.0);  // beta
    parameters_.set(1, 1.0);  // c

    SX vel = ca.vertcat(List.of(dx, dy));
    SX p = ca.vertcat(List.of(x, y));
    SX q = ca.vertcat(List.of(u, v));

    // System dynamics
    SX F = k.unaryMinus().times(p.minus(q)).minus(beta.times(v).times(ca.sqrt(ca.sumsqr(vel).plus(c.power(2)))));

    // System right hand side
    SX rhs = ca.vertcat(List.of(dx, dy, F.get(0), F.get(1)));

    // The remainder of the Python example is commented-out legacy code (ControlSimulator, covariance
    // propagation by integrator sensitivities, Lyapunov equations and unscented propagation).
  }
}

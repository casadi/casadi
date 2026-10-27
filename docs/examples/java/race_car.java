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

// Car race along a track
// ----------------------
// An optimal control problem (OCP),
// solved with direct multiple-shooting.
//
// For more information see: http://labs.casadi.org/OCP
import org.casadi.*;
import java.util.List;
import java.util.function.BinaryOperator;

public class race_car {
  public static void main(String[] args) {
    if (!ca.has_nlpsol("ipopt")) { System.out.println("SKIP: nlpsol ipopt not available"); return; }

    int N = 100; // number of control intervals

    Opti opti = new Opti(); // Optimization problem

    // ---- decision variables ---------
    MX X = opti.variable(2, N + 1); // state trajectory
    MX pos   = X.get(0, new Slice());
    MX speed = X.get(1, new Slice());
    MX U = opti.variable(1, N);   // control trajectory (throttle)
    MX T = opti.variable();       // final time

    // ---- objective          ---------
    opti.minimize(T); // race in minimal time

    // ---- dynamic constraints --------
    BinaryOperator<MX> f = (x, u) -> ca.vertcat(List.of(x.get(1), u.minus(x.get(1)))); // dx/dt = f(x,u)

    MX dt = T.div(N); // length of a control interval
    for (int k = 0; k < N; k++) { // loop over control intervals
      // Runge-Kutta 4 integration
      MX Xk = X.get(new Slice(), k), Uk = U.get(new Slice(), k);
      MX k1 = f.apply(Xk,                          Uk);
      MX k2 = f.apply(Xk.plus(dt.div(2).times(k1)), Uk);
      MX k3 = f.apply(Xk.plus(dt.div(2).times(k2)), Uk);
      MX k4 = f.apply(Xk.plus(dt.times(k3)),        Uk);
      MX x_next = Xk.plus(dt.div(6).times(k1.plus(k2.times(2)).plus(k3.times(2)).plus(k4)));
      opti.subject_to(ca.eq(X.get(new Slice(), k + 1), x_next)); // close the gaps
    }

    // ---- path constraints -----------
    opti.subject_to(ca.le(speed, limit(pos)));   // track speed limit
    opti.subject_to(Opti.bounded(0, U, 1));      // control is limited

    // ---- boundary conditions --------
    opti.subject_to(ca.eq(pos.get(0), 0));   // start at position 0 ...
    opti.subject_to(ca.eq(speed.get(0), 0)); // ... from stand-still
    opti.subject_to(ca.eq(pos.get(-1), 1));  // finish line at position 1

    // ---- misc. constraints  ----------
    opti.subject_to(ca.ge(T, 0)); // Time must be positive

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

    DM jac_g = sol.value(ca.jacobian(opti.g(), opti.x()));
    // plot: spy(jac_g)
    DM hess_l = sol.value((MX) ca.hessian(opti.f().plus(ca.dot(opti.lam_g(), opti.g())), opti.x()).get(0));
    // plot: spy(hess_l)
  }

  // track speed limit (Python uses one lambda for both MX and DM)
  static MX limit(MX pos) { return ca.minus(1, ca.sin(pos.times(2 * Math.PI)).div(2)); }
  static DM limit(DM pos) { return ca.minus(1, ca.sin(pos.times(2 * Math.PI)).div(2)); }
}

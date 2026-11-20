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

use casadi as ca;
use std::f64::consts::PI;

fn main() {
    let n = 100; // number of control intervals

    let opti = ca::Opti::default(); // Optimization problem

    // ---- decision variables ---------
    let x = opti.variable_with(2, n + 1, "full"); // state trajectory
    let pos = x.at((0, ..));
    let speed = x.at((1, ..));
    let u = opti.variable_with(1, n, "full"); // control trajectory (throttle)
    let t = opti.variable(); // final time

    // ---- objective          ---------
    opti.minimize(&t); // race in minimal time

    // ---- dynamic constraints --------
    let f = |x: &ca::MX, u: &ca::MX| ca::vertcat(vec![x.at(1), u - x.at(1)]); // dx/dt = f(x,u)

    let dt = &t / n as f64; // length of a control interval
    for k in 0..n {
        // Runge-Kutta 4 integration
        let xk = x.at((.., k));
        let uk = u.at((.., k));
        let k1 = f(&xk, &uk);
        let k2 = f(&(&xk + &dt / 2.0 * &k1), &uk);
        let k3 = f(&(&xk + &dt / 2.0 * &k2), &uk);
        let k4 = f(&(&xk + &dt * &k3), &uk);
        let x_next = &xk + &dt / 6.0 * (&k1 + 2.0 * &k2 + 2.0 * &k3 + &k4);
        opti.subject_to(&ca::eq(&x.at((.., k + 1)), &x_next)); // close the gaps
    }

    // ---- path constraints -----------
    let limit = |pos: &ca::MX| 1.0 - ca::sin(&(2.0 * PI * pos)) / 2.0;
    opti.subject_to(&ca::le(&speed, &limit(&pos))); // track speed limit
    opti.subject_to(&ca::Opti::bounded(0.0, &u, 1.0)); // control is limited

    // ---- boundary conditions --------
    opti.subject_to(&ca::eq(&pos.at(0), 0.0)); // start at position 0 ...
    opti.subject_to(&ca::eq(&speed.at(0), 0.0)); // ... from stand-still
    opti.subject_to(&ca::eq(&pos.at(-1), 1.0)); // finish line at position 1

    // ---- misc. constraints  ----------
    opti.subject_to(&ca::ge(&t, 0.0)); // Time must be positive

    // ---- initial values for solver ---
    opti.set_initial(&speed, 1.0);
    opti.set_initial(&t, 1.0);

    // ---- solve NLP              ------
    opti.solver("ipopt"); // set numerical backend
    let sol = opti.solve(); // actual solve

    // ---- post-processing        ------
    println!("minimal lap time: {}", sol.value(&t));
    println!("speed: {}", sol.value(&speed));
}

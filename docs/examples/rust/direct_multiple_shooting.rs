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

use casadi as ca;
use casadi::dict;
use std::collections::HashMap;

const INF: f64 = f64::INFINITY;

fn main() {
    let t = 10.0; // Time horizon
    let n = 20; // number of control intervals

    // Declare model variables
    let x1 = ca::MX::sym("x1");
    let x2 = ca::MX::sym("x2");
    let x = ca::vertcat(vec![&x1, &x2]);
    let u = ca::MX::sym("u");

    // Model equations
    let xdot = ca::vertcat(vec![(1.0 - ca::power(&x2, 2.0)) * &x1 - &x2 + &u, x1.clone()]);

    // Objective term
    let l = ca::power(&x1, 2.0) + ca::power(&x2, 2.0) + ca::power(&u, 2.0);

    // Formulate discrete time dynamics
    let big_f = if false {
        // CVODES from the SUNDIALS suite
        let dae = dict! {"x" => &x, "p" => &u, "ode" => &xdot, "quad" => &l};
        ca::integrator_with("F", "cvodes", &dae, 0.0, t / n as f64, dict! {})
    } else {
        // Fixed step Runge-Kutta 4 integrator
        let m = 4; // RK4 steps per interval
        let dt = t / n as f64 / m as f64;
        let f = ca::Function::new("f", vec![&x, &u], vec![&xdot, &l]);
        let x0 = ca::MX::sym_with("X0", 2, 1);
        let uu = ca::MX::sym("U");
        let mut xx = x0.clone();
        let mut q = ca::MX::new(0.0);
        for _ in 0..m {
            let (k1, k1_q): (ca::MX, ca::MX) = f.eval((&xx, &uu));
            let (k2, k2_q): (ca::MX, ca::MX) = f.eval((&xx + dt / 2.0 * &k1, &uu));
            let (k3, k3_q): (ca::MX, ca::MX) = f.eval((&xx + dt / 2.0 * &k2, &uu));
            let (k4, k4_q): (ca::MX, ca::MX) = f.eval((&xx + dt * &k3, &uu));
            xx = &xx + dt / 6.0 * (&k1 + 2.0 * &k2 + 2.0 * &k3 + &k4);
            q = &q + dt / 6.0 * (&k1_q + 2.0 * &k2_q + 2.0 * &k3_q + &k4_q);
        }
        ca::Function::named("F", vec![&x0, &uu], vec![&xx, &q], vec!["x0", "p"], vec!["xf", "qf"])
    };

    // Evaluate at a test point
    let fk = big_f.call(&dict! {"x0" => vec![0.2, 0.3], "p" => 0.4}).unwrap();
    println!("{}", fk["xf"]);
    println!("{}", fk["qf"]);

    // Start with an empty NLP
    let mut w: Vec<ca::MX> = vec![];
    let mut w0: Vec<f64> = vec![];
    let mut lbw: Vec<f64> = vec![];
    let mut ubw: Vec<f64> = vec![];
    let mut j = ca::MX::new(0.0);
    let mut g: Vec<ca::MX> = vec![];
    let mut lbg: Vec<f64> = vec![];
    let mut ubg: Vec<f64> = vec![];

    // "Lift" initial conditions
    let mut xk = ca::MX::sym_with("X0", 2, 1);
    w.push(xk.clone());
    lbw.extend([0.0, 1.0]);
    ubw.extend([0.0, 1.0]);
    w0.extend([0.0, 1.0]);

    // Formulate the NLP
    for k in 0..n {
        // New NLP variable for the control
        let uk = ca::MX::sym(&format!("U_{}", k));
        w.push(uk.clone());
        lbw.push(-1.0);
        ubw.push(1.0);
        w0.push(0.0);

        // Integrate till the end of the interval
        let fk = big_f.call_symbolic(&HashMap::from([("x0", &xk), ("p", &uk)]));
        let xk_end = &fk["xf"];
        j = &j + &fk["qf"];

        // New NLP variable for state at end of interval
        xk = ca::MX::sym_with(&format!("X_{}", k + 1), 2, 1);
        w.push(xk.clone());
        lbw.extend([-0.25, -INF]);
        ubw.extend([INF, INF]);
        w0.extend([0.0, 0.0]);

        // Add equality constraint
        g.push(xk_end - &xk);
        lbg.extend([0.0, 0.0]);
        ubg.extend([0.0, 0.0]);
    }

    // Create an NLP solver
    let prob = dict! {"f" => &j, "x" => ca::vertcat(&w), "g" => ca::vertcat(&g)};
    let solver = ca::nlpsol("solver", "ipopt", &prob);

    // Solve the NLP
    let sol = solver.call(&dict! {"x0" => &w0, "lbx" => &lbw, "ubx" => &ubw, "lbg" => &lbg, "ubg" => &ubg}).unwrap();
    let w_opt = sol["x"].nonzeros();

    // Print the solution
    let x1_opt: Vec<f64> = w_opt.iter().step_by(3).copied().collect();
    let x2_opt: Vec<f64> = w_opt.iter().skip(1).step_by(3).copied().collect();
    let u_opt: Vec<f64> = w_opt.iter().skip(2).step_by(3).copied().collect();

    let tgrid: Vec<f64> = (0..=n).map(|k| t / n as f64 * k as f64).collect();
    println!("optimal cost: {}", sol["f"]);
    println!("t:  {:?}", tgrid);
    println!("x1: {:?}", x1_opt);
    println!("x2: {:?}", x2_opt);
    println!("u:  {:?}", u_opt);
}

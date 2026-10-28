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
//
// Solves the following optimal control problem (OCP) in differential-algebraic
// equations (DAE)
//
// minimize     integral_{t=0}^{10} x0^2 + x1^2 + u^2  dt
// x0,x1,z,u
//
// subject to   dot(x0) == z*x0-x1+u     \
//              dot(x1) == x0             }  for 0 <= t <= 10
//                    0 == x1^2 + z - 1  /
//              x0(t=0) == 0
//              x1(t=0) == 1
//              x0(t=10) == 0
//              x1(t=10) == 0
//              -0.75 <= u <= 1  for 0 <= t <= 10
//
// The method used is direct multiple shooting.
//
// Joel Andersson, 2015

use casadi as ca;
use casadi::dict;
use std::collections::HashMap;

const INF: f64 = f64::INFINITY;

fn main() {
    // Declare variables
    let x0 = ca::SX::sym("x0");
    let x1 = ca::SX::sym("x1");
    let x = ca::vertcat(vec![&x0, &x1]); // Differential states
    let z = ca::SX::sym("z"); // Algebraic variable
    let u = ca::SX::sym("u"); // Control

    // Differential equation
    let f_x = ca::vertcat(vec![&z * &x0 - &x1 + &u, x0.clone()]);

    // Algebraic equation
    let f_z = ca::power(&x1, 2.0) + &z - 1.0;

    // Lagrange cost term (quadrature)
    let f_q = ca::power(&x0, 2.0) + ca::power(&x1, 2.0) + ca::power(&u, 2.0);

    // Create an integrator
    let dae = dict! {"x" => &x, "z" => &z, "p" => &u, "ode" => &f_x, "alg" => &f_z, "quad" => &f_q};
    // interval length 0.5s
    let i = ca::integrator_with("I", "idas", &dae, 0.0, 0.5, dict! {});

    // Number of intervals
    let nk = 20;

    // Start with an empty NLP
    let mut w: Vec<ca::MX> = vec![]; // List of variables
    let mut lbw: Vec<f64> = vec![]; // Lower bounds on w
    let mut ubw: Vec<f64> = vec![]; // Upper bounds on w
    let mut gg: Vec<ca::MX> = vec![]; // Constraints
    let mut j = ca::MX::new(0.0); // Cost function

    // Initial conditions
    let mut xk = ca::MX::sym_with("X0", 2, 1);
    w.push(xk.clone());
    lbw.extend([0.0, 1.0]);
    ubw.extend([0.0, 1.0]);

    // Loop over all intervals
    for k in 0..nk {
        // Local control
        let uk = ca::MX::sym(&format!("U{}", k));
        w.push(uk.clone());
        lbw.push(-0.75);
        ubw.push(1.00);

        // Call integrator function
        let ik = i.call_symbolic(&HashMap::from([("x0", &xk), ("p", &uk)]));
        xk = ik["xf"].clone();
        j = j + &ik["qf"]; // Sum quadratures

        // "Lift" the variable
        let x_prev = xk;
        xk = ca::MX::sym_with(&format!("X{}", k + 1), 2, 1);
        w.push(xk.clone());
        lbw.extend([-INF, -INF]);
        ubw.extend([INF, INF]);
        gg.push(x_prev - &xk);
    }

    // Allocate an NLP solver
    let nlp = dict! {"x" => ca::vertcat(&w), "f" => &j, "g" => ca::vertcat(&gg)};
    let opts = dict! {"ipopt.linear_solver" => "ma27"};
    let solver = ca::nlpsol_with("solver", "ipopt", &nlp, &opts);

    // Pass bounds, initial guess and solve NLP
    let sol = solver
        .call(&dict! {
            "lbx" => &lbw, // Lower variable bound
            "ubx" => &ubw, // Upper variable bound
            "lbg" => 0.0,  // Lower constraint bound
            "ubg" => 0.0,  // Upper constraint bound
            "x0" => 0.0,   // Initial guess
        })
        .unwrap();

    // Print the results
    let every3 = |k: i64| ca::Slice::new_with(k, i64::MAX, 3);
    println!("optimal cost: {}", sol["f"]);
    println!("t:  {}", ca::linspace(0.0, 10.0, nk + 1));
    println!("x0 trajectory: {}", sol["x"].at(every3(0)));
    println!("x1 trajectory: {}", sol["x"].at(every3(1)));
    println!("u trajectory:  {}", sol["x"].at(every3(2)));
}

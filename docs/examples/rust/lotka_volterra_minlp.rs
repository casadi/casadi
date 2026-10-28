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

// Lotka-Volerra problem from mintoc.de

use casadi as ca;
use casadi::dict;
use std::collections::HashMap;

fn main() {
    let t = 12.0; // Time horizon
    let n = 40; // number of control intervals

    // Model parameters
    let c0 = 0.4;
    let c1 = 0.2;

    // Initial condition for x
    let x0 = vec![0.5, 0.7];

    // Bounds on x
    let lbx = vec![0.0, 0.0];
    let ubx = vec![2.0, 2.0];

    // Bounds on u
    let lbu = 0.0;
    let ubu = 1.0;

    // Declare model variables
    let x = ca::SX::sym_with("x", 2, 1);
    let u = ca::SX::sym("u");

    // Model equations
    let xdot = ca::vertcat(vec![x.at(0) - x.at(0) * x.at(1) - c0 * x.at(0) * &u,
                                -x.at(1) + x.at(0) * x.at(1) - c1 * x.at(1) * &u]);

    // Objective term
    let l = ca::power(&(x.at(0) - 1.0), 2.0) + ca::power(&(x.at(1) - 1.0), 2.0) + 1e-4 * ca::power(&u, 2.0);

    // Formulate discrete time dynamics
    let f_dyn = if false {
        // CVODES from the SUNDIALS suite
        let dae = dict! {"x" => &x, "p" => &u, "ode" => &xdot, "quad" => &l};
        ca::integrator_with("F", "cvodes", &dae, 0.0, t / n as f64, dict! {})
    } else {
        // Fixed step Runge-Kutta 4 integrator
        let m = 4; // RK4 steps per interval
        let dt = t / n as f64 / m as f64;
        let f = ca::Function::new("f", vec![&x, &u], vec![&xdot, &l]);
        let x0_sym = ca::MX::sym_with("X0", 2, 1);
        let u_sym = ca::MX::sym("U");
        let mut xs = x0_sym.clone();
        let mut q = ca::MX::new(0.0);
        for _ in 0..m {
            let (k1, k1_q): (ca::MX, ca::MX) = f.eval((&xs, &u_sym));
            let (k2, k2_q): (ca::MX, ca::MX) = f.eval((&(&xs + dt / 2.0 * &k1), &u_sym));
            let (k3, k3_q): (ca::MX, ca::MX) = f.eval((&(&xs + dt / 2.0 * &k2), &u_sym));
            let (k4, k4_q): (ca::MX, ca::MX) = f.eval((&(&xs + dt * &k3), &u_sym));
            xs = &xs + dt / 6.0 * (&k1 + 2.0 * &k2 + 2.0 * &k3 + &k4);
            q = &q + dt / 6.0 * (&k1_q + 2.0 * &k2_q + 2.0 * &k3_q + &k4_q);
        }
        ca::Function::named("F", vec![&x0_sym, &u_sym], vec![&xs, &q], vec!["x0", "p"], vec!["xf", "qf"])
    };

    // Initial guess for u
    let u_start = vec![ca::DM::new(0.0); n];

    // Get a feasible trajectory as an initial guess
    let mut xk = ca::DM::new(x0.clone());
    let mut x_start = vec![xk.clone()];
    for k in 0..n {
        let res = f_dyn.call(&dict! {"x0" => &xk, "p" => &u_start[k]}).unwrap();
        xk = res["xf"].clone();
        x_start.push(xk.clone());
    }

    // Start with an empty NLP
    let mut w: Vec<ca::MX> = vec![];
    let mut w0: Vec<ca::DM> = vec![];
    let mut lbw: Vec<f64> = vec![];
    let mut ubw: Vec<f64> = vec![];
    let mut discrete: Vec<bool> = vec![];
    let mut j = ca::MX::new(0.0);
    let mut g: Vec<ca::MX> = vec![];
    let mut lbg: Vec<f64> = vec![];
    let mut ubg: Vec<f64> = vec![];

    // "Lift" initial conditions
    let x0_sym = ca::MX::sym_with("X0", 2, 1);
    w.push(x0_sym.clone());
    lbw.extend(&x0);
    ubw.extend(&x0);
    w0.push(x_start[0].clone());
    discrete.extend([false, false]);

    // Formulate the NLP
    let mut xk = x0_sym;
    for k in 0..n {
        // New NLP variable for the control
        let uk = ca::MX::sym(format!("U_{}", k));
        w.push(uk.clone());
        lbw.push(lbu);
        ubw.push(ubu);
        w0.push(u_start[k].clone());
        discrete.push(true);

        // Integrate till the end of the interval
        let fk = f_dyn.call_symbolic(HashMap::from([("x0", xk.clone()), ("p", uk)]));
        let xk_end = &fk["xf"];
        j = &j + &fk["qf"];

        // New NLP variable for state at end of interval
        xk = ca::MX::sym_with(format!("X_{}", k + 1), 2, 1);
        w.push(xk.clone());
        lbw.extend(&lbx);
        ubw.extend(&ubx);
        w0.push(x_start[k + 1].clone());
        discrete.extend([false, false]);

        // Add equality constraint
        g.push(xk_end - &xk);
        lbg.extend([0.0, 0.0]);
        ubg.extend([0.0, 0.0]);
    }

    // Concatenate decision variables and constraint terms
    let w = ca::vertcat(w);
    let g = ca::vertcat(g);

    // Create an NLP solver
    let nlp_prob = dict! {"f" => &j, "x" => &w, "g" => &g};
    let nlp_solver = ca::nlpsol_with("nlp_solver", "bonmin", &nlp_prob, dict! {"discrete" => &discrete});
    //let nlp_solver = ca::nlpsol_with("nlp_solver", "knitro", &nlp_prob, dict! {"discrete" => &discrete});
    //let nlp_solver = ca::nlpsol("nlp_solver", "ipopt", &nlp_prob); // Solve relaxed problem

    // Print the solution (instead of plotting it)
    let tgrid: Vec<f64> = (0..=n).map(|k| t / n as f64 * k as f64).collect();
    let plot_sol = |w_opt: &ca::DM| {
        let w_opt = w_opt.nonzeros();
        let x0_opt: Vec<f64> = w_opt.iter().step_by(3).copied().collect();
        let x1_opt: Vec<f64> = w_opt.iter().skip(1).step_by(3).copied().collect();
        let u_opt: Vec<f64> = w_opt.iter().skip(2).step_by(3).copied().collect();
        println!("t:  {:?}", tgrid);
        println!("x0: {:?}", x0_opt);
        println!("x1: {:?}", x1_opt);
        println!("u:  {:?}", u_opt);
    };

    // Solve the NLP
    let sol = nlp_solver.call(&dict! {"x0" => ca::vertcat(w0), "lbx" => &lbw, "ubx" => &ubw,
                                      "lbg" => &lbg, "ubg" => &ubg}).unwrap();

    println!("{:?}", nlp_solver.stats());

    let w1_opt = &sol["x"];
    let _lam_w_opt = &sol["lam_x"];
    let _lam_g_opt = &sol["lam_g"];
    plot_sol(w1_opt);
}

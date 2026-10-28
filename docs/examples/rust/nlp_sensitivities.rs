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
// Direct-collocation Van der Pol OCP solved with an SQP method (qrqp QP
// sub-solver) for accurate multipliers, then NLP-solution sensitivities w.r.t.
// a perturbation parameter P: Hessian via factory, forward AD vs finite
// differences, and reverse AD.

use casadi as ca;
use casadi::dict;

// Polynomial helpers (stand in for numpy.poly1d/polyder/polyint), coefficients highest power first
fn polymul(a: &[f64], b: &[f64]) -> Vec<f64> {
    let mut c = vec![0.0; a.len() + b.len() - 1];
    for (i, ai) in a.iter().enumerate() {
        for (j, bj) in b.iter().enumerate() {
            c[i + j] += ai * bj;
        }
    }
    c
}
fn polyval(p: &[f64], x: f64) -> f64 {
    p.iter().fold(0.0, |acc, c| acc * x + c)
}
fn polyder(p: &[f64]) -> Vec<f64> {
    let n = p.len() - 1;
    p[..n].iter().enumerate().map(|(i, c)| c * (n - i) as f64).collect()
}
fn polyint(p: &[f64]) -> Vec<f64> {
    let n = p.len();
    let mut q: Vec<f64> = p.iter().enumerate().map(|(i, c)| c / (n - i) as f64).collect();
    q.push(0.0);
    q
}

fn main() {
    // Degree of interpolating polynomial
    let d = 3;
    // Get collocation points
    let mut tau_root = vec![0.0];
    tau_root.extend(ca::collocation_points_with(d, "legendre"));
    let d = d as usize;
    // Coefficients of the collocation equation
    let mut c = vec![vec![0.0; d + 1]; d + 1];
    // Coefficients of the continuity equation
    let mut dd = vec![0.0; d + 1];
    // Coefficients of the quadrature function
    let mut b = vec![0.0; d + 1];
    // Construct polynomial basis
    for j in 0..d + 1 {
        // Construct Lagrange polynomials to get the polynomial basis at the collocation point
        let mut p = vec![1.0];
        for r in 0..d + 1 {
            if r != j {
                let den = tau_root[j] - tau_root[r];
                p = polymul(&p, &[1.0 / den, -tau_root[r] / den]);
            }
        }
        // Evaluate the polynomial at the final time to get the coefficients of the continuity equation
        dd[j] = polyval(&p, 1.0);
        // Evaluate the time derivative of the polynomial at all collocation points to get the coefficients of the continuity equation
        let pder = polyder(&p);
        for r in 0..d + 1 {
            c[j][r] = polyval(&pder, tau_root[r]);
        }

        // Evaluate the integral of the polynomial to get the coefficients of the quadrature function
        let pint = polyint(&p);
        b[j] = polyval(&pint, 1.0);
    }

    // Time horizon
    let t = 10.0;

    // Declare model variables
    let x1 = ca::SX::sym("x1");
    let x2 = ca::SX::sym("x2");
    let x = ca::vertcat(vec![&x1, &x2]);
    let u = ca::SX::sym("u");

    // Model equations
    let xdot = ca::vertcat(vec![(1.0 - ca::power(&x2, 2.0)) * &x1 - &x2 + &u, x1.clone()]);

    // Objective term
    let l = ca::power(&x1, 2.0) + ca::power(&x2, 2.0) + ca::power(&u, 2.0);

    // Continuous time dynamics
    let f = ca::Function::named("f", vec![&x, &u], vec![&xdot, &l], vec!["x", "u"], vec!["xdot", "L"]);

    // Control discretization
    let n = 20; // number of control intervals
    let h = t / n as f64;

    // Start with an empty NLP
    let mut w: Vec<ca::MX> = vec![];
    let mut w0: Vec<f64> = vec![];
    let mut lbw: Vec<f64> = vec![];
    let mut ubw: Vec<f64> = vec![];
    let mut j_obj = ca::MX::new(0.0);
    let mut g: Vec<ca::MX> = vec![];
    let mut lbg: Vec<f64> = vec![];
    let mut ubg: Vec<f64> = vec![];

    // For plotting x and u given w
    let mut x_plot: Vec<ca::MX> = vec![];
    let mut u_plot: Vec<ca::MX> = vec![];

    // "Lift" initial conditions
    let mut xk = ca::MX::sym_with("X0", 2, 1);
    w.push(xk.clone());
    lbw.extend([0.0, 1.0]);
    ubw.extend([0.0, 1.0]);
    w0.extend([0.0, 1.0]);
    x_plot.push(xk.clone());

    // Perturb with P
    let p = ca::MX::sym_with("P", 2, 1);
    xk = &xk + &p;

    // Formulate the NLP
    for k in 0..n {
        // New NLP variable for the control
        let uk = ca::MX::sym(&format!("U_{}", k));
        w.push(uk.clone());
        lbw.push(-1.0);
        ubw.push(0.85);
        w0.push(0.0);
        u_plot.push(uk.clone());

        // State at collocation points
        let mut xc: Vec<ca::MX> = vec![];
        for j in 0..d {
            let xkj = ca::MX::sym_with(&format!("X_{}_{}", k, j), 2, 1);
            xc.push(xkj.clone());
            w.push(xkj);
            lbw.extend([-0.25, -f64::INFINITY]);
            ubw.extend([f64::INFINITY, f64::INFINITY]);
            w0.extend([0.0, 0.0]);
        }

        // Loop over collocation points
        let mut xk_end = dd[0] * &xk;
        for j in 1..d + 1 {
            // Expression for the state derivative at the collocation point
            let mut xp = c[0][j] * &xk;
            for r in 0..d {
                xp = xp + c[r + 1][j] * &xc[r];
            }

            // Append collocation equations
            let (fj, qj): (ca::MX, ca::MX) = f.eval((&xc[j - 1], &uk));
            g.push(h * &fj - &xp);
            lbg.extend([0.0, 0.0]);
            ubg.extend([0.0, 0.0]);

            // Add contribution to the end state
            xk_end = xk_end + dd[j] * &xc[j - 1];

            // Add contribution to quadrature function
            j_obj = j_obj + b[j] * &qj * h;
        }

        // New NLP variable for state at end of interval
        xk = ca::MX::sym_with(&format!("X_{}", k + 1), 2, 1);
        w.push(xk.clone());
        lbw.extend([-0.25, -f64::INFINITY]);
        ubw.extend([f64::INFINITY, f64::INFINITY]);
        w0.extend([0.0, 0.0]);
        x_plot.push(xk.clone());

        // Add equality constraint
        g.push(&xk_end - &xk);
        lbg.extend([0.0, 0.0]);
        ubg.extend([0.0, 0.0]);
    }

    // Concatenate vectors
    let w = ca::vertcat(w);
    let g = ca::vertcat(g);
    let x_plot = ca::horzcat(x_plot);
    let u_plot = ca::horzcat(u_plot);
    let (w0, lbw, ubw, lbg, ubg) = (ca::DM::new(w0), ca::DM::new(lbw), ca::DM::new(ubw), ca::DM::new(lbg), ca::DM::new(ubg));

    // Create an NLP solver
    let prob = dict! {"f" => &j_obj, "x" => &w, "g" => &g, "p" => &p};

    // Function to get x and u trajectories from w
    let trajectories = ca::Function::named("trajectories", vec![&w], vec![&x_plot, &u_plot], vec!["w"], vec!["x", "u"]);

    // Create an NLP solver, using SQP and active-set QP for accurate multipliers
    let opts = dict! {"qpsol" => "qrqp", "qpsol_options" => dict! {"print_iter" => false, "error_on_fail" => false}, "print_time" => false};
    let solver = ca::nlpsol_with("solver", "sqpmethod", &prob, &opts);

    // Solve the NLP
    let sol = solver.call(&dict! {"x0" => &w0, "lbx" => &lbw, "ubx" => &ubw, "lbg" => &lbg, "ubg" => &ubg, "p" => 0.0}).unwrap();

    // Extract trajectories
    let (x_opt, u_opt): (ca::DM, ca::DM) = trajectories.eval([&sol["x"]]);

    // Plot the result (printed instead)
    println!("x1: {}", x_opt.at((0, ..)));
    println!("x2: {}", x_opt.at((1, ..)));
    println!("u: {}", u_opt);

    // High-level approach:
    // Use factory to e.g. to calculate Hessian of optimal f w.r.t. p
    let hsolver = solver.factory("h", solver.name_in(), vec!["hess:f:p:p"]);
    println!("hsolver generated");
    let hsol = hsolver
        .call(&dict! {"x0" => &sol["x"], "lam_x0" => &sol["lam_x"], "lam_g0" => &sol["lam_g"],
        "lbx" => &lbw, "ubx" => &ubw, "lbg" => &lbg, "ubg" => &ubg, "p" => 0.0})
        .unwrap();
    println!("Hessian of f w.r.t. p:");
    println!("{}", hsol["hess_f_p_p"]);

    // Low-level approach: Calculate directional derivatives.
    // We will calculate two directions simultaneously
    let nfwd = 2;

    // Forward mode AD for the NLP solver object
    let fwd_solver = solver.forward(nfwd);
    println!("fwd_solver generated");

    // Seeds, initalized to zero
    let zeros = |s: &ca::DM| (0..nfwd).map(|_| ca::DM::zeros_sparsity(&s.sparsity())).collect::<Vec<_>>();
    let fwd_lbx = zeros(&sol["x"]);
    let fwd_ubx = zeros(&sol["x"]);
    let mut fwd_p: Vec<ca::DM> = (0..nfwd).map(|_| ca::DM::zeros_sparsity(&p.sparsity())).collect();
    let fwd_lbg = zeros(&sol["g"]);
    let fwd_ubg = zeros(&sol["g"]);

    // Let's preturb P
    fwd_p[0].set_at(0, 1.0); // first nonzero of P
    fwd_p[1].set_at(1, 1.0); // second nonzero of P

    // Calculate sensitivities using AD
    let sol_fwd = fwd_solver
        .call(&dict! {"out_x" => &sol["x"], "out_lam_g" => &sol["lam_g"], "out_lam_x" => &sol["lam_x"],
        "out_f" => &sol["f"], "out_g" => &sol["g"], "lbx" => &lbw, "ubx" => &ubw, "lbg" => &lbg, "ubg" => &ubg,
        "fwd_lbx" => ca::horzcat(fwd_lbx.clone()), "fwd_ubx" => ca::horzcat(fwd_ubx.clone()),
        "fwd_lbg" => ca::horzcat(fwd_lbg.clone()), "fwd_ubg" => ca::horzcat(fwd_ubg.clone()),
        "p" => 0.0, "fwd_p" => ca::horzcat(fwd_p.clone())})
        .unwrap();

    // Calculate the same thing using finite differences
    let h = 1e-3;
    let perturbed = |s: f64, d: usize| {
        solver
            .call(&dict! {"x0" => &sol["x"], "lam_g0" => &sol["lam_g"], "lam_x0" => &sol["lam_x"],
            "lbx" => &lbw + s * h * (&fwd_lbx[d] + &fwd_ubx[d]),
            "ubx" => &ubw + s * h * (&fwd_lbx[d] + &fwd_ubx[d]),
            "lbg" => &lbg + s * h * (&fwd_lbg[d] + &fwd_ubg[d]),
            "ubg" => &ubg + s * h * (&fwd_lbg[d] + &fwd_ubg[d]),
            "p" => 0.0 + s * h * &fwd_p[d]})
            .unwrap()
    };
    let pert: Vec<_> = (0..nfwd as usize).map(|d| perturbed(1.0, d)).collect();

    // Print the result
    for s in ["f"] {
        println!("==========");
        println!("Checking {}", s);
        println!("finite differences");
        for d in 0..nfwd as usize {
            println!("{}", (&pert[d][s] - &sol[s]) / h);
        }
        println!("AD fwd");
        let m = &sol_fwd[&format!("fwd_{}", s)];
        for d in 0..nfwd {
            println!("{}", m.at((.., d)));
        }
    }

    // Perturb again, in the opposite direction for second order derivatives
    let pert2: Vec<_> = (0..nfwd as usize).map(|d| perturbed(-1.0, d)).collect();

    // Print the result
    for s in ["f"] {
        println!("finite differences, second order: {}", s);
        for d in 0..nfwd as usize {
            println!("{}", (&pert[d][s] - 2.0 * &sol[s] + &pert2[d][s]) / (h * h));
        }
    }

    // Reverse mode AD for the NLP solver object
    let nadj = 1;
    let adj_solver = solver.reverse(nadj);
    println!("adj_solver generated");

    // Seeds, initalized to zero
    let mut adj_f: Vec<ca::DM> = (0..nadj).map(|_| ca::DM::zeros_sparsity(&sol["f"].sparsity())).collect();
    let adj_g: Vec<ca::DM> = (0..nadj).map(|_| ca::DM::zeros_sparsity(&sol["g"].sparsity())).collect();
    let adj_x: Vec<ca::DM> = (0..nadj).map(|_| ca::DM::zeros_sparsity(&sol["x"].sparsity())).collect();

    // Study which inputs influence f
    adj_f[0].set_at(0, 1.0);

    // Calculate sensitivities using AD
    let sol_adj = adj_solver
        .call(&dict! {"out_x" => &sol["x"], "out_lam_g" => &sol["lam_g"], "out_lam_x" => &sol["lam_x"],
        "out_f" => &sol["f"], "out_g" => &sol["g"], "lbx" => &lbw, "ubx" => &ubw, "lbg" => &lbg, "ubg" => &ubg,
        "adj_f" => ca::horzcat(adj_f), "adj_g" => ca::horzcat(adj_g),
        "p" => 0.0, "adj_x" => ca::horzcat(adj_x)})
        .unwrap();

    // Print the result
    for s in ["p"] {
        println!("==========");
        println!("Checking {}", s);
        println!("Reverse mode AD");
        println!("{}", sol_adj[&format!("adj_{}", s)]);
    }
}

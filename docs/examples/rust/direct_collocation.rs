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

const INF: f64 = f64::INFINITY;

// Polynomial with coefficients in decreasing powers, like numpy's poly1d
#[derive(Clone)]
struct Poly(Vec<f64>);

impl Poly {
    fn mul(&self, o: &Poly) -> Poly {
        let mut c = vec![0.0; self.0.len() + o.0.len() - 1];
        for (i, a) in self.0.iter().enumerate() {
            for (j, b) in o.0.iter().enumerate() {
                c[i + j] += a * b;
            }
        }
        Poly(c)
    }
    fn eval(&self, x: f64) -> f64 {
        self.0.iter().fold(0.0, |acc, c| acc * x + c)
    }
    fn der(&self) -> Poly {
        let n = self.0.len() - 1;
        Poly(self.0[..n].iter().enumerate().map(|(i, c)| c * (n - i) as f64).collect())
    }
    fn int(&self) -> Poly {
        let n = self.0.len();
        let mut c: Vec<f64> = self.0.iter().enumerate().map(|(i, c)| c / (n - i) as f64).collect();
        c.push(0.0);
        Poly(c)
    }
}

fn main() {
    // Degree of interpolating polynomial
    let d = 3;

    // Get collocation points
    let mut tau_root = vec![0.0];
    tau_root.extend(ca::collocation_points_with(d as i64, "legendre"));

    // Coefficients of the collocation equation
    let mut c = vec![vec![0.0; d + 1]; d + 1];

    // Coefficients of the continuity equation
    let mut dd = vec![0.0; d + 1];

    // Coefficients of the quadrature function
    let mut b = vec![0.0; d + 1];

    // Construct polynomial basis
    for j in 0..=d {
        // Construct Lagrange polynomials to get the polynomial basis at the collocation point
        let mut p = Poly(vec![1.0]);
        for r in 0..=d {
            if r != j {
                let s = tau_root[j] - tau_root[r];
                p = p.mul(&Poly(vec![1.0 / s, -tau_root[r] / s]));
            }
        }

        // Evaluate the polynomial at the final time to get the coefficients of the continuity equation
        dd[j] = p.eval(1.0);

        // Evaluate the time derivative of the polynomial at all collocation points to get the coefficients of the continuity equation
        let pder = p.der();
        for r in 0..=d {
            c[j][r] = pder.eval(tau_root[r]);
        }

        // Evaluate the integral of the polynomial to get the coefficients of the quadrature function
        let pint = p.int();
        b[j] = pint.eval(1.0);
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
    let mut jj = ca::MX::new(0.0);
    let mut g: Vec<ca::MX> = vec![];
    let mut lbg: Vec<f64> = vec![];
    let mut ubg: Vec<f64> = vec![];

    // For plotting x and u given w
    let mut x_plot: Vec<ca::MX> = vec![];
    let mut u_plot: Vec<ca::MX> = vec![];

    // "Lift" initial conditions
    let mut xk = ca::MX::sym_with("X0", 2, None);
    w.push(xk.clone());
    lbw.extend([0.0, 1.0]);
    ubw.extend([0.0, 1.0]);
    w0.extend([0.0, 1.0]);
    x_plot.push(xk.clone());

    // Formulate the NLP
    for k in 0..n {
        // New NLP variable for the control
        let uk = ca::MX::sym(&format!("U_{}", k));
        w.push(uk.clone());
        lbw.push(-1.0);
        ubw.push(1.0);
        w0.push(0.0);
        u_plot.push(uk.clone());

        // State at collocation points
        let mut xc: Vec<ca::MX> = vec![];
        for j in 0..d {
            let xkj = ca::MX::sym_with(&format!("X_{}_{}", k, j), 2, 1);
            xc.push(xkj.clone());
            w.push(xkj);
            lbw.extend([-0.25, -INF]);
            ubw.extend([INF, INF]);
            w0.extend([0.0, 0.0]);
        }

        // Loop over collocation points
        let mut xk_end = dd[0] * &xk;
        for j in 1..=d {
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
            jj = jj + b[j] * qj * h;
        }

        // New NLP variable for state at end of interval
        xk = ca::MX::sym_with(&format!("X_{}", k + 1), 2, 1);
        w.push(xk.clone());
        lbw.extend([-0.25, -INF]);
        ubw.extend([INF, INF]);
        w0.extend([0.0, 0.0]);
        x_plot.push(xk.clone());

        // Add equality constraint
        g.push(&xk_end - &xk);
        lbg.extend([0.0, 0.0]);
        ubg.extend([0.0, 0.0]);
    }

    // Concatenate vectors
    let w = ca::vertcat(&w);
    let g = ca::vertcat(&g);
    let x_plot = ca::horzcat(&x_plot);
    let u_plot = ca::horzcat(&u_plot);

    // Create an NLP solver
    let prob = dict! {"f" => &jj, "x" => &w, "g" => &g};
    let solver = ca::nlpsol("solver", "ipopt", &prob);

    // Function to get x and u trajectories from w
    let trajectories = ca::Function::named("trajectories", vec![&w], vec![&x_plot, &u_plot], vec!["w"], vec!["x", "u"]);

    // Solve the NLP
    let sol = solver.call(&dict! {"x0" => &w0, "lbx" => &lbw, "ubx" => &ubw, "lbg" => &lbg, "ubg" => &ubg}).unwrap();
    let (x_opt, u_opt): (ca::DM, ca::DM) = trajectories.eval([&sol["x"]]);

    // Print the result
    let tgrid: Vec<f64> = (0..=n).map(|k| t * k as f64 / n as f64).collect();
    println!("optimal cost: {}", sol["f"]);
    println!("t:  {:?}", tgrid);
    println!("x1: {}", x_opt.at((0, ..)));
    println!("x2: {}", x_opt.at((1, ..)));
    println!("u:  {}", u_opt);
}

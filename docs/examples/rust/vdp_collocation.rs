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

// numpy.poly1d stand-ins: coefficients, highest power first
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

    // Choose collocation points
    let tau_root: Vec<f64> = std::iter::once(0.0).chain(ca::collocation_points_with(d as i64, "radau")).collect();

    // Coefficients of the collocation equation
    let mut c = vec![vec![0.0; d + 1]; d + 1];

    // Coefficients of the continuity equation
    let mut dd = vec![0.0; d + 1];

    // Coefficients of the quadrature function
    let mut ff = vec![0.0; d + 1];

    // Construct polynomial basis
    for j in 0..=d {
        // Construct Lagrange polynomials to get the polynomial basis at the collocation point
        let mut p = vec![1.0];
        for r in 0..=d {
            if r != j {
                let den = tau_root[j] - tau_root[r];
                p = polymul(&p, &[1.0 / den, -tau_root[r] / den]);
            }
        }

        // Evaluate the polynomial at the final time to get the coefficients of the continuity equation
        dd[j] = polyval(&p, 1.0);

        // Evaluate the time derivative of the polynomial at all collocation points to get the coefficients of the continuity equation
        let pder = polyder(&p);
        for r in 0..=d {
            c[j][r] = polyval(&pder, tau_root[r]);
        }

        // Evaluate the integral of the polynomial to get the coefficients of the quadrature function
        let pint = polyint(&p);
        ff[j] = polyval(&pint, 1.0);
    }

    // Control discretization
    let nk = 20;

    // End time
    let tf = 10.0;

    // Size of the finite elements
    let h = tf / nk as f64;

    // All collocation time points
    let t_grid: Vec<Vec<f64>> = (0..nk).map(|k| (0..=d).map(|j| h * (k as f64 + tau_root[j])).collect()).collect();

    // Declare variables (use scalar graph)
    let t = ca::SX::sym("t"); // time
    let u = ca::SX::sym("u"); // control
    let x = ca::SX::sym_with("x", 2, 1); // state

    // ODE rhs function and quadratures
    let xdot = ca::vertcat(vec![(1.0 - x.at(1) * x.at(1)) * x.at(0) - x.at(1) + &u, x.at(0)]);
    let qdot = x.at(0) * x.at(0) + x.at(1) * x.at(1) + &u * &u;
    let f = ca::Function::new("f", vec![&t, &x, &u], vec![&xdot, &qdot]);

    // Control bounds
    let u_min = -0.75;
    let u_max = 1.0;
    let u_init = 0.0;

    // State bounds and initial guess
    let x_min = [-INF, -INF];
    let x_max = [INF, INF];
    let xi_min = [0.0, 1.0];
    let xi_max = [0.0, 1.0];
    let xf_min = [0.0, 0.0];
    let xf_max = [0.0, 0.0];
    let x_init = [0.0, 0.0];

    // Dimensions
    let nx = 2;
    let nu = 1;

    // Total number of variables
    let n_x = nk * (d + 1) * nx; // Collocated states
    let n_u = nk * nu; // Parametrized controls
    let n_xf = nx; // Final state
    let nv = n_x + n_u + n_xf;

    // NLP variable vector
    let v = ca::MX::sym_with("V", nv as i64, 1);

    // All variables with bounds and initial guess
    let mut vars_lb = vec![0.0; nv];
    let mut vars_ub = vec![0.0; nv];
    let mut vars_init = vec![0.0; nv];
    let mut offset = 0;

    // Get collocated states and parametrized control
    let mut xs: Vec<Vec<ca::MX>> = vec![vec![ca::MX::default(); d + 1]; nk + 1];
    let mut us: Vec<ca::MX> = vec![ca::MX::default(); nk];
    for k in 0..nk {
        // Collocated states
        for j in 0..=d {
            // Get the expression for the state vector
            xs[k][j] = v.at(offset..offset + nx);

            // Add the initial condition
            vars_init[offset..offset + nx].copy_from_slice(&x_init);

            // Add bounds
            if k == 0 && j == 0 {
                vars_lb[offset..offset + nx].copy_from_slice(&xi_min);
                vars_ub[offset..offset + nx].copy_from_slice(&xi_max);
            } else {
                vars_lb[offset..offset + nx].copy_from_slice(&x_min);
                vars_ub[offset..offset + nx].copy_from_slice(&x_max);
            }
            offset += nx;
        }

        // Parametrized controls
        us[k] = v.at(offset..offset + nu);
        vars_lb[offset..offset + nu].fill(u_min);
        vars_ub[offset..offset + nu].fill(u_max);
        vars_init[offset..offset + nu].fill(u_init);
        offset += nu;
    }

    // State at end time
    xs[nk][0] = v.at(offset..offset + nx);
    vars_lb[offset..offset + nx].copy_from_slice(&xf_min);
    vars_ub[offset..offset + nx].copy_from_slice(&xf_max);
    vars_init[offset..offset + nx].copy_from_slice(&x_init);
    offset += nx;
    assert_eq!(offset, nv);

    // Constraint function for the NLP
    let mut g: Vec<ca::MX> = vec![];
    let mut lbg: Vec<f64> = vec![];
    let mut ubg: Vec<f64> = vec![];

    // Objective function
    let mut jj = ca::MX::new(0.0);

    // For all finite elements
    for k in 0..nk {
        // For all collocation points
        for j in 1..=d {
            // Get an expression for the state derivative at the collocation point
            let mut xp_jk = ca::MX::new(0.0);
            for r in 0..=d {
                xp_jk += c[r][j] * &xs[k][r];
            }

            // Add collocation equations to the NLP
            let (fk, qk): (ca::MX, ca::MX) = f.eval((t_grid[k][j], &xs[k][j], &us[k]));
            g.push(h * fk - xp_jk);
            lbg.extend(vec![0.0; nx]); // equality constraints
            ubg.extend(vec![0.0; nx]); // equality constraints

            // Add contribution to objective
            jj += ff[j] * qk * h;
        }

        // Get an expression for the state at the end of the finite element
        let mut xf_k = ca::MX::new(0.0);
        for r in 0..=d {
            xf_k += dd[r] * &xs[k][r];
        }

        // Add continuity equation to NLP
        g.push(&xs[k + 1][0] - xf_k);
        lbg.extend(vec![0.0; nx]);
        ubg.extend(vec![0.0; nx]);
    }

    // Concatenate constraints
    let g = ca::vertcat(g);

    // NLP
    let nlp = dict! {"x" => &v, "f" => &jj, "g" => &g};

    // ----
    // SOLVE THE NLP
    // ----

    // Set options
    let opts = dict! {
        "expand" => true,
        //"ipopt.max_iter" => 4,
        "ipopt.linear_solver" => "ma27",
    };

    // Allocate an NLP solver
    let solver = ca::nlpsol_with("solver", "ipopt", &nlp, &opts);
    let arg = dict! {
        // Initial condition
        "x0" => vars_init,
        // Bounds on x
        "lbx" => vars_lb,
        "ubx" => vars_ub,
        // Bounds on g
        "lbg" => lbg,
        "ubg" => ubg,
    };

    // Solve the problem
    let res = solver.call(&arg).unwrap();

    // Print the optimal cost
    println!("optimal cost:  {}", res["f"].to_double());

    // Retrieve the solution
    let v_opt = res["x"].nonzeros();

    // Get values at the beginning of each finite element
    let stride = (d + 1) * nx + nu;
    let x0_opt: Vec<f64> = v_opt.iter().step_by(stride).copied().collect();
    let x1_opt: Vec<f64> = v_opt[1..].iter().step_by(stride).copied().collect();
    let u_opt: Vec<f64> = v_opt[(d + 1) * nx..].iter().step_by(stride).copied().collect();

    // Print the results
    println!("x[0] trajectory: {:?}", x0_opt);
    println!("x[1] trajectory: {:?}", x1_opt);
    println!("u trajectory: {:?}", u_opt);
}

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

fn main() {
    let nk = 20; // Control discretization
    let tf = 10.0; // End time

    // Declare variables (use scalar graph)
    let t = ca::SX::sym("t"); // time
    let u = ca::SX::sym("u"); // control

    // casadi.tools structures are not available: the state vector is laid out by hand
    // states = [x (2), L (1)]
    let states = ca::SX::sym_with("states", 3, 1);
    let x = states.at(0..2); //  vdp oscillator states
    let l = states.at(2); //  helper state: Langrange integrand

    // Create the right hand side
    let rhs = ca::vertcat(vec![
        (1.0 - x.at(1) * x.at(1)) * x.at(0) - x.at(1) + &u,
        x.at(0),
        x.at(0) * x.at(0) + x.at(1) * x.at(1) + &u * &u,
    ]);

    // ODE right hand side function
    let f = ca::Function::new("f", vec![&t, &states, &u], vec![&rhs]);

    // Objective function (meyer term)
    let m = ca::Function::new("m", vec![&t, &states, &u], vec![&l]);

    // Control bounds
    let u_min = -0.75;
    let u_max = 1.0;
    let u_init = 0.0;

    // State bounds and initial guess
    let x_min = [-INF, -INF, -INF];
    let x_max = [INF, INF, INF];
    let xi_min = [0.0, 1.0, 0.0];
    let xi_max = [0.0, 1.0, 0.0];
    let xf_min = [0.0, 0.0, -INF];
    let xf_max = [0.0, 0.0, INF];
    let x_init = [0.0, 0.0, 0.0];

    // Dimensions
    let nx = 3;
    let nu = 1;

    // Choose collocation points
    let tau_root: Vec<f64> = std::iter::once(0.0).chain(ca::collocation_points_with(3, "radau")).collect();

    // Degree of interpolating polynomial
    let d = tau_root.len() - 1;

    // Size of the finite elements
    let h = tf / nk as f64;

    // Coefficients of the collocation equation
    let mut c = vec![vec![0.0; d + 1]; d + 1];

    // Coefficients of the continuity equation
    let mut dd = vec![0.0; d + 1];

    // Dimensionless time inside one control interval
    let tau = ca::SX::sym("tau");

    // All collocation time points
    let t_grid: Vec<Vec<f64>> = (0..nk).map(|k| (0..=d).map(|j| h * (k as f64 + tau_root[j])).collect()).collect();

    // For all collocation points
    for j in 0..=d {
        // Construct Lagrange polynomials to get the polynomial basis at the collocation point
        let mut lp = ca::SX::new(1.0);
        for r in 0..=d {
            if r != j {
                lp *= (&tau - tau_root[r]) / (tau_root[j] - tau_root[r]);
            }
        }

        // Evaluate the polynomial at the final time to get the coefficients of the continuity equation
        let lfcn = ca::Function::new("lfcn", vec![&tau], vec![&lp]);
        let v: ca::DM = lfcn.eval([1.0]);
        dd[j] = v.to_double();

        // Evaluate the time derivative of the polynomial at all collocation points to get the coefficients of the continuity equation
        let tfcn = ca::Function::new("tfcn", vec![&tau], vec![ca::tangent(&lp, &tau)]);
        for r in 0..=d {
            let v: ca::DM = tfcn.eval([tau_root[r]]);
            c[j][r] = v.to_double();
        }
    }

    // NLP variables, in the order of the casadi.tools structure
    //   struct_symMX([(entry("X",repeat=[nk+1,d+1],struct=states), entry("U",repeat=[nk],shape=nu))])
    // which interleaves X[k,:] and U[k]
    let nv = (nk + 1) * (d + 1) * nx + nk * nu;
    let vv = ca::MX::sym_with("V", nv as i64, 1);
    let ix = |k: usize, j: usize| k * ((d + 1) * nx + nu) + j * nx;
    let iu = |k: usize| k * ((d + 1) * nx + nu) + (d + 1) * nx;
    let vx = |k: usize, j: usize| vv.at(ix(k, j)..ix(k, j) + nx);
    let vu = |k: usize| vv.at(iu(k)..iu(k) + nu);

    let mut vars_lb = vec![0.0; nv];
    let mut vars_ub = vec![0.0; nv];
    let mut vars_init = vec![0.0; nv];

    // Set states and its bounds
    for k in 0..=nk {
        for j in 0..=d {
            vars_init[ix(k, j)..ix(k, j) + nx].copy_from_slice(&x_init);
            vars_lb[ix(k, j)..ix(k, j) + nx].copy_from_slice(&x_min);
            vars_ub[ix(k, j)..ix(k, j) + nx].copy_from_slice(&x_max);
        }
    }

    // Set controls and its bounds
    for k in 0..nk {
        vars_init[iu(k)] = u_init;
        vars_lb[iu(k)] = u_min;
        vars_ub[iu(k)] = u_max;
    }

    // State at initial time
    vars_lb[ix(0, 0)..ix(0, 0) + nx].copy_from_slice(&xi_min);
    vars_ub[ix(0, 0)..ix(0, 0) + nx].copy_from_slice(&xi_max);

    // State at end time
    vars_lb[ix(nk, 0)..ix(nk, 0) + nx].copy_from_slice(&xf_min);
    vars_ub[ix(nk, 0)..ix(nk, 0) + nx].copy_from_slice(&xf_max);

    // Constraint function for the NLP
    let mut g: Vec<ca::MX> = vec![];
    let mut lbg: Vec<f64> = vec![];
    let mut ubg: Vec<f64> = vec![];

    // For all finite elements
    for k in 0..nk {
        // For all collocation points
        for j in 1..=d {
            // Get an expression for the state derivative at the collocation point
            let mut xp_jk = ca::MX::new(0.0);
            for r in 0..=d {
                xp_jk += c[r][j] * vx(k, r);
            }

            // Add collocation equations to the NLP
            let fk: ca::MX = f.eval((t_grid[k][j], vx(k, j), vu(k)));
            g.push(h * fk - xp_jk);
            lbg.extend(vec![0.0; nx]); // equality constraints
            ubg.extend(vec![0.0; nx]); // equality constraints
        }

        // Get an expression for the state at the end of the finite element
        let mut xf_k = ca::MX::new(0.0);
        for r in 0..=d {
            xf_k += dd[r] * vx(k, r);
        }

        // Add continuity equation to NLP
        g.push(vx(k + 1, 0) - xf_k);
        lbg.extend(vec![0.0; nx]);
        ubg.extend(vec![0.0; nx]);
    }

    // Concatenate constraints
    let g = ca::vertcat(g);

    // Objective function
    let obj: ca::MX = m.eval((t_grid[nk - 1][d], vx(nk, 0), vu(nk - 1)));

    // NLP
    let nlp = dict! {"x" => &vv, "f" => &obj, "g" => &g};

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
    let opt = res["x"].nonzeros();

    // Get values at the beginning of each finite element
    let x0_opt: Vec<f64> = (0..=nk).map(|k| opt[ix(k, 0)]).collect();
    let x1_opt: Vec<f64> = (0..=nk).map(|k| opt[ix(k, 0) + 1]).collect();
    let x2_opt: Vec<f64> = (0..=nk).map(|k| opt[ix(k, 0) + 2]).collect();
    let u_opt: Vec<f64> = (0..nk).map(|k| opt[iu(k)]).collect();

    // Print the results
    println!("x[0] trajectory: {:?}", x0_opt);
    println!("x[1] trajectory: {:?}", x1_opt);
    println!("L trajectory: {:?}", x2_opt);
    println!("u trajectory: {:?}", u_opt);
}

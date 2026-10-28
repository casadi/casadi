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
// Demonstration on how to construct a fixed-step implicit Runge-Kutta integrator
// @author: Joel Andersson, KU Leuven 2013

use casadi as ca;
use casadi::dict;
use std::collections::HashMap;

fn main() {
    // End time
    let tf = 10.0;

    // Dimensions
    let nx = 3;
    let np = 1;

    // Declare variables
    let x = ca::SX::sym_with("x", nx, 1); // state
    let p = ca::SX::sym_with("u", np, 1); // control

    // ODE right hand side function
    let ode = ca::vertcat(vec![
        (1.0 - x.at(1) * x.at(1)) * x.at(0) - x.at(1) + &p,
        x.at(0),
        x.at(0) * x.at(0) + x.at(1) * x.at(1) + &p * &p,
    ]);
    let dae = dict! {"x" => &x, "p" => &p, "ode" => &ode};
    let f = ca::Function::new("f", vec![&x, &p], vec![&ode]);

    // Number of finite elements
    let n = 100;

    // Size of the finite elements
    let h = tf / n as f64;

    // Degree of interpolating polynomial
    let d = 4;

    // Choose collocation points
    let mut tau_root = vec![0.0];
    tau_root.extend(ca::collocation_points_with(d as i64, "legendre"));

    // Coefficients of the collocation equation
    let mut c = vec![vec![0.0; d + 1]; d + 1];

    // Coefficients of the continuity equation
    let mut dd = vec![0.0; d + 1];

    // Dimensionless time inside one control interval
    let tau = ca::SX::sym("tau");

    // For all collocation points
    for j in 0..=d {
        // Construct Lagrange polynomials to get the polynomial basis at the collocation point
        let mut l = ca::SX::new(1.0);
        for r in 0..=d {
            if r != j {
                l *= (&tau - tau_root[r]) / (tau_root[j] - tau_root[r]);
            }
        }

        // Evaluate the polynomial at the final time to get the coefficients of the continuity equation
        let lfcn = ca::Function::new("lfcn", vec![&tau], vec![&l]);
        dd[j] = lfcn.eval::<ca::DM>([1.0]).to_double();

        // Evaluate the time derivative of the polynomial at all collocation points to get the coefficients of the continuity equation
        let tfcn = ca::Function::new("tfcn", vec![&tau], vec![ca::tangent(&l, &tau)]);
        for r in 0..=d {
            c[j][r] = tfcn.eval::<ca::DM>([tau_root[r]]).to_double();
        }
    }

    // Total number of variables for one finite element
    let x0 = ca::MX::sym_with("X0", nx, 1);
    let pp = ca::MX::sym_with("P", np, 1);
    let v = ca::MX::sym_with("V", d as i64 * nx, 1);

    // Get the state at each collocation point
    let mut xx = vec![x0.clone()];
    xx.extend(ca::vertsplit_with(&v, (0..=d as i64).map(|r| r * nx).collect::<Vec<_>>()));

    // Get the collocation quations (that define V)
    let mut v_eq: Vec<ca::MX> = vec![];
    for j in 1..=d {
        // Expression for the state derivative at the collocation point
        let mut xp_j = ca::MX::new(0.0);
        for r in 0..=d {
            xp_j += c[r][j] * &xx[r];
        }

        // Append collocation equations
        let f_j: ca::MX = f.eval((&xx[j], &pp));
        v_eq.push(h * f_j - xp_j);
    }

    // Concatenate constraints
    let v_eq = ca::vertcat(&v_eq);

    // Root-finding function, implicitly defines V as a function of X0 and P
    let vfcn = ca::Function::new("vfcn", vec![&v, &x0, &pp], vec![&v_eq]);

    // Convert to SX to decrease overhead
    let vfcn_sx = vfcn.expand();

    // Create a implicit function instance to solve the system of equations
    let ifcn = ca::rootfinder("ifcn", "newton", &vfcn_sx);
    let v: ca::MX = ifcn.eval((ca::MX::default(), &x0, &pp));
    let xx: Vec<ca::MX> =
        (0..=d as i64).map(|r| if r == 0 { x0.clone() } else { v.at((r - 1) * nx..r * nx) }).collect();

    // Get an expression for the state at the end of the finie element
    let mut xf = ca::MX::new(0.0);
    for r in 0..=d {
        xf += dd[r] * &xx[r];
    }

    // Get the discrete time dynamics
    let ff = ca::Function::new("F", vec![&x0, &pp], vec![&xf]);

    // Do this iteratively for all finite elements
    let mut xx = x0.clone();
    for _ in 0..n {
        xx = ff.eval((&xx, &pp));
    }

    // Fixed-step integrator
    let irk_integrator = ca::Function::from_dict(
        "irk_integrator",
        &HashMap::from([("x0", &x0), ("p", &pp), ("xf", &xx)]),
        ca::integrator_in(),
        ca::integrator_out(),
    );

    // Create a convensional integrator for reference
    let ref_integrator = ca::integrator_with("ref_integrator", "cvodes", &dae, 0.0, tf, dict! {});

    // Test values
    let x0_val = vec![0.0, 1.0, 0.0];
    let p_val = 0.2;

    // Make sure that both integrators give consistent results
    for ff in [&irk_integrator, &ref_integrator] {
        println!("-------");
        println!("Testing {}", ff.name());
        println!("-------");

        // Generate a new function that calculates forward and reverse directional derivatives
        let df =
            ff.factory("dF", vec!["x0", "p", "fwd:x0", "fwd:p", "adj:xf"], vec!["xf", "fwd:xf", "adj:x0", "adj:p"]);
        let arg = dict! {
            // Pass arguments
            "x0" => &x0_val,
            "p" => p_val,

            // Forward sensitivity analysis, first direction: seed p and x0[0]
            "fwd_x0" => vec![1.0, 0.0, 0.0],
            "fwd_p" => 1.0,

            // Adjoint sensitivity analysis, seed xf[2]
            "adj_xf" => vec![0.0, 0.0, 1.0],
        };

        // Integrate
        let res = df.call(&arg).unwrap();

        // Get the nondifferentiated results
        println!("{:>30} = {}", "xf", res["xf"]);

        // Get the forward sensitivities
        println!("{:>30} = {}", "d(xf)/d(p)+d(xf)/d(x0[0])", res["fwd_xf"]);

        // Get the adjoint sensitivities
        println!("{:>30} = {}", "d(xf[2])/d(x0)", res["adj_x0"]);
        println!("{:>30} = {}", "d(xf[2])/d(p)", res["adj_p"]);
    }
}

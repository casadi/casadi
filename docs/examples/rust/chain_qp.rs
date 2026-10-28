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
// We want to model a chain attached to two supports and hanging in between. Let us discretise
// it with N mass points connected by N-1 springs. Each mass i has position (yi,zi), i=1,...,N.
// The equilibrium point of the system minimises the potential energy.
//
// The potential energy of each spring is
// Vi=D_i/2 * ((y_i-y_{i+1})^2 + (z_i-z_{i+1})^2)
//
// The gravitational potential energy of each mass is
// Vg_i = m_i*g0*z_i
//
// The total potential energy is thus given by:
//
// Vchain(y,z) = 1/2*sum{i=1,...,N-1} D_i ((y_i-y_{i+1})^2+(z_i-z_{i+1})^2) + g0 * sum{i=1,...,N} m_i * z_i
//
// where y=[y_1,...,y_N] and z=[z_1,...,z_N]
//
// We wish to solve
// minimize{y,z} Vchain(y, z)
//
// Subject to the piecewise linear ground constraints:
// z_i >= zin
// z_i - 0.1*y_i >= 0.5

use casadi as ca;
use casadi::dict;

const INF: f64 = f64::INFINITY;

fn main() {
    // Constants
    let n = 40;
    let m_i = 40.0 / n as f64;
    let d_i = 70.0 * n as f64;
    let g0 = 9.81;
    //let zmin = -INF; // unbounded
    let zmin = 0.5; // ground

    // Objective function
    let mut vchain = ca::SX::zeros();

    // Variables
    let mut x: Vec<ca::SX> = vec![];

    // Variable bounds
    let mut lbx: Vec<f64> = vec![];
    let mut ubx: Vec<f64> = vec![];

    // Constraints
    let mut g: Vec<ca::SX> = vec![];

    // Constraint bounds
    let mut lbg: Vec<f64> = vec![];
    let mut ubg: Vec<f64> = vec![];

    // Loop over all chain elements
    let mut prev: Option<(ca::SX, ca::SX)> = None;
    for i in 1..=n {
        // Create variables for the (y_i, z_i) coordinates
        let y_i = ca::SX::sym(&format!("y_{}", i));
        let z_i = ca::SX::sym(&format!("z_{}", i));

        // Add to the list of variables
        x.push(y_i.clone());
        x.push(z_i.clone());
        if i == 1 {
            lbx.extend([-2., 1.]);
            ubx.extend([-2., 1.]);
        } else if i == n {
            lbx.extend([2., 1.]);
            ubx.extend([2., 1.]);
        } else {
            lbx.extend([-INF, zmin]);
            ubx.extend([INF, INF]);
        }

        // Spring potential
        if let Some((y_prev, z_prev)) = &prev {
            vchain += d_i / 2.0 * (ca::power(y_prev - &y_i, 2.0) + ca::power(z_prev - &z_i, 2.0));
        }

        // Graviational potential
        vchain += g0 * m_i * &z_i;

        // Slanted ground constraints
        g.push(&z_i - 0.1 * &y_i);
        lbg.push(0.5);
        ubg.push(INF);

        // Previous point
        prev = Some((y_i, z_i));
    }

    // Formulate QP
    let qp = dict! {"x" => ca::vertcat(&x), "f" => &vchain, "g" => ca::vertcat(&g)};

    // Solve with IPOPT
    let solver = ca::qpsol_with("solver", "qpoases", &qp, &dict! {"sparse" => true});
    //let solver = ca::qpsol("solver", "gurobi", &qp);
    //let solver = ca::nlpsol("solver", "ipopt", &qp);

    // Get the optimal solution
    let sol = solver
        .call(&dict! {"lbx" => &lbx, "ubx" => &ubx, "lbg" => &lbg, "ubg" => &ubg})
        .unwrap();
    let x_opt = &sol["x"];
    let f_opt = &sol["f"];
    println!("f_opt =  {}", f_opt);

    // Retrieve the result
    let y0 = x_opt.at(ca::Slice::new_with(0, i64::MAX, 2));
    let z0 = x_opt.at(ca::Slice::new_with(1, i64::MAX, 2));

    // Plot the result (plotting dropped: print the chain instead)
    println!("y = {}", y0);
    println!("z = {}", z0);
    let ys: ca::DM = ca::linspace(-2.0, 2.0, 100);
    let zs = 0.5 + 0.1 * &ys;
    println!("ground line z - 0.1y >= 0.5 from ({}, {}) to ({}, {})", ys.at(0), zs.at(0), ys.at(-1), zs.at(-1));
}

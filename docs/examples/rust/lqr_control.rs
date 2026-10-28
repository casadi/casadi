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
// This example how an Linear Quadratic Regulator (LQR) can be designed
// and simulated for a linear-time-invariant (LTI) system
//
// Only the live part of the Python original is ported: everything after the
// open-loop setup is commented out there (it uses the removed ControlSimulator).

use casadi as ca;

fn main() {
    // Problem statement
    // -----------------------------------

    // System dimensions
    let ns = 3; // number of states
    let nu = 2; // number of controls
    let _ny = 2; // number of outputs
    let n = 21; // number of control intervals
    let te = 10.0; // end time [s]

    // The system equations:   x' = A.x + B.u
    let a = ca::DM::new(vec![vec![0.4, 0.1, -2.0], vec![0.0, -0.3, 4.0], vec![1.0, 0.0, 0.0]]);
    let b = ca::DM::new(vec![vec![1.0, 1.0], vec![0.0, 1.0], vec![1.0, 0.0]]);

    // Inspect the open-loop system
    // (numpy's linalg.eig: casadi has no numeric non-symmetric eigensolver, so
    //  take the roots of the characteristic polynomial det(A - l*I))
    let l = ca::SX::sym("l");
    let charpoly = ca::det(&(ca::SX::new(&a) - &l * ca::SX::eye(ns)));
    let coeff = ca::evalf(&ca::poly_coeff(&charpoly, &l)).nonzeros();
    let d = cubic_roots(&coeff);
    println!("Open-loop eigenvalues:  {:?}", d);

    let _tn = ca::linspace(0.0, te, n);

    // Check if the system is controllable
    // -----------------------------------
    let mut r: Vec<ca::DM> = vec![];

    let mut p = b.clone();
    for _i in 0..ns {
        r.push(p.clone());
        p = ca::mtimes_list(vec![a.clone(), p]);
    }

    let r = ca::horzcat(r);

    // R must be of full rank
    // (no numeric SVD in casadi: count the nonzero diagonal entries of the QR factor of R')
    let (_q, rr): (ca::DM, ca::DM) = ca::qr(&r.T());
    let eps = 1e-9;
    let rank = (0..ns).filter(|&i| rr.at((i, i)).to_double().abs() > eps).count();
    assert_eq!(rank, ns as usize);

    // Simulation of the open-loop system
    // -----------------------------------

    let _y = ca::SX::sym_with("y", ns, None);
    let _u = ca::SX::sym_with("u", nu, None);

    let _x0 = ca::DM::new(vec![1.0, 0.0, 0.0]);
    // no control
    let u_ = ca::repmat_with(&ca::DM::new(vec![vec![-1.0, 1.0], vec![1.0, -1.0]]), 1, (n - 1) / 2);

    println!("{}", u_);

    let _p = ca::SX::sym("p");
}

// Roots of c[0] x^3 + c[1] x^2 + c[2] x + c[3], as (re, im) pairs
fn cubic_roots(c: &[f64]) -> Vec<(f64, f64)> {
    let (a, b, cc) = (c[1] / c[0], c[2] / c[0], c[3] / c[0]);
    // One real root by Newton's method
    let mut x = 0.0;
    for _ in 0..100 {
        x -= (((x + a) * x + b) * x + cc) / ((3.0 * x + 2.0 * a) * x + b);
    }
    // Deflate to x^2 + p x + q
    let (p, q) = (a + x, b + (a + x) * x);
    let disc = p * p / 4.0 - q;
    let mut roots = if disc >= 0.0 {
        vec![(-p / 2.0 + disc.sqrt(), 0.0), (-p / 2.0 - disc.sqrt(), 0.0)]
    } else {
        vec![(-p / 2.0, (-disc).sqrt()), (-p / 2.0, -(-disc).sqrt())]
    };
    roots.push((x, 0.0));
    roots
}

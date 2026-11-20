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

fn main() {
    // Excercise 1, chapter 10 from Larry Biegler's book
    println!("program started");

    // Test with different number of elements
    for n in 1..11 {
        println!("N =  {}", n);

        // Degree of interpolating polynomial
        let kk = 2;

        // Legrandre roots
        let tau_root = [0., 0.211325, 0.788675];

        // Radau roots (K=3)
        //tau_root = [0, 0.155051, 0.644949, 1]

        // Time
        let t = ca::SX::sym("t");

        // Differential equation
        let z = ca::SX::sym("z");
        let f = ca::Function::new("dz_dt", vec![&z], vec![&z * &z - 2.0 * &z + 1.0]);

        let z0 = -3.0;

        // Analytic solution
        let _z_analytic = ca::Function::new("z_analytic", vec![&t], vec![(4.0 * &t - 3.0) / (3.0 * &t + 1.0)]);

        // Collocation point
        let tau = ca::SX::sym("tau");

        // Step size
        let h = 1.0 / n as f64;

        // Get the coefficients of the continuity and collocation equations
        let mut d = ca::DM::zeros_with(kk + 1, 1);
        let mut c = ca::DM::zeros_with(kk + 1, kk + 1);
        for j in 0..kk + 1 {
            // Lagrange polynomial
            let mut l = ca::SX::new(1.0);
            for k in 0..kk + 1 {
                if k != j {
                    l *= (&tau - tau_root[k as usize]) / (tau_root[j as usize] - tau_root[k as usize]);
                }
            }

            // Evaluate at end for coefficients of continuity equation
            let lfcn = ca::Function::new("lfcn", vec![&tau], vec![&l]);
            let lj: ca::DM = lfcn.eval([1.0]);
            d.set_at(j, &lj);

            // Differentiate and evaluate at collocation points
            let tfcn = ca::Function::new("tfcn", vec![&tau], vec![ca::tangent(&l, &tau)]);
            for k in 0..kk + 1 {
                let cjk: ca::DM = tfcn.eval([tau_root[k as usize]]);
                c.set_at((j, k), &cjk);
            }
        }
        println!("C =  {}", c);
        println!("D =  {}", d);

        // Collocated states
        let zz = ca::SX::sym_with("Z", n, kk + 1);

        // Construct the NLP
        let x = ca::vec(&zz.T());
        let mut g: Vec<ca::SX> = vec![];
        for i in 0..n {
            for k in 1..kk + 1 {
                // Add collocation equations to NLP
                let mut rhs = ca::SX::new(0.0);
                for j in 0..kk + 1 {
                    rhs += zz.at((i, j)) * c.at((j, k));
                }
                let ff: ca::SX = f.eval([zz.at((i, k))]);
                g.push(h * ff - rhs);
            }

            // Add continuity equation to NLP
            let mut rhs = ca::SX::new(0.0);
            for j in 0..kk + 1 {
                rhs += d.at(j) * zz.at((i, j));
            }

            if i < n - 1 {
                g.push(zz.at((i + 1, 0)) - rhs);
            }
        }

        let g = ca::vertcat(g);

        println!("g =  {}", g);

        // NLP
        let nlp = dict! {"x" => &x, "f" => ca::power(&x.at(0), 2.0), "g" => &g};

        // ## ----
        // ## SOLVE THE NLP
        // ## ----

        // NLP solver options
        let opts = dict! {"ipopt.tol" => 1e-10};

        // Allocate an NLP solver and buffer
        let solver = ca::nlpsol_with("solver", "ipopt", &nlp, &opts);

        // Initial condition
        let x0 = vec![0.0; x.nnz() as usize];

        // Bounds on x
        let mut lbx = vec![-100.0; x.nnz() as usize];
        let mut ubx = vec![100.0; x.nnz() as usize];
        lbx[0] = z0;
        ubx[0] = z0;

        // Bounds on the constraints, solve the problem
        let res = solver.call(&dict! {"x0" => x0, "lbx" => lbx, "ubx" => ubx, "lbg" => 0.0, "ubg" => 0.0}).unwrap();

        // ## Print the time points
        let n = n as usize;
        let kk = kk as usize;
        let mut t_opt = vec![0.0; n * (kk + 1)];
        for i in 0..n {
            for j in 0..kk + 1 {
                t_opt[j + (kk + 1) * i] = h * (i as f64 + tau_root[j]);
            }
        }

        println!("time points:  {:?}", t_opt);

        // Print the optimal cost
        println!("optimal cost:  {:?}", res["f"].to_double());

        // Print the optimal solution
        let xopt = res["x"].nonzeros();
        println!("optimal solution:  {:?}", xopt);
    }
}

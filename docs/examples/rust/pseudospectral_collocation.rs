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
// Legendre-Gauss-Lobatto transcription of an optimal-control problem.

use casadi as ca;
use std::f64::consts::PI;

/// Computes legendre polynomial P_N(x) and its derivative P'_N(x), and P''_N(x) using numerically
/// stable three-term recurrence relations.
fn legendre_all(n: usize, x: &[f64]) -> (Vec<f64>, Vec<f64>, Vec<f64>) {
    x.iter()
        .map(|&x| {
            let (mut p0, mut p1) = (1.0, x);
            let (mut dp0, mut dp1) = (0.0, 1.0);
            let (mut d2p0, mut d2p1) = (0.0, 0.0);

            if n == 0 {
                return (p0, dp0, d2p0);
            }

            for k in 1..n {
                let k = k as f64;
                let p2 = ((2.0 * k + 1.0) * x * p1 - k * p0) / (k + 1.0);
                let dp2 = ((2.0 * k + 1.0) * (p1 + x * dp1) - k * dp0) / (k + 1.0);
                let d2p2 = ((2.0 * k + 1.0) * (2.0 * dp1 + x * d2p1) - k * d2p0) / (k + 1.0);

                (p0, p1) = (p1, p2);
                (dp0, dp1) = (dp1, dp2);
                (d2p0, d2p1) = (d2p1, d2p2);
            }
            (p1, dp1, d2p1)
        })
        .fold((vec![], vec![], vec![]), |(mut a, mut b, mut c), (p, dp, d2p)| {
            a.push(p);
            b.push(dp);
            c.push(d2p);
            (a, b, c)
        })
}

/// Given N intervals, computes N+1 Legendre-Gauss-Lobatto (LGL) points on [-1, 1]
///
/// Roots of (1-x**2)P'_N(x) are LGL points
fn lgl_nodes(n: usize) -> Vec<f64> {
    let (tol, max_iter) = (1e-15, 100);
    assert!(n >= 1, "Degree N must be at least 1.");

    // inital guess
    let mut x: Vec<f64> = (0..=n).map(|k| -(PI * k as f64 / n as f64).cos()).collect();

    // begin newton iterations
    if n > 1 {
        let mut x_int = x[1..n].to_vec();
        for _ in 0..max_iter {
            let (_, dp, d2p) = legendre_all(n, &x_int);
            let dx: Vec<f64> = dp.iter().zip(&d2p).map(|(a, b)| a / b).collect();
            for (xi, d) in x_int.iter_mut().zip(&dx) {
                *xi -= d;
            }
            if dx.iter().fold(0.0f64, |m, d| m.max(d.abs())) < tol {
                break;
            }
        }
        x[1..n].copy_from_slice(&x_int);
    }

    x[0] = -1.0;
    x[n] = 1.0;

    // use symmetry to remove round-off errors
    (0..=n).map(|i| 0.5 * (x[i] - x[n - i])).collect()
}

/// Computes LGL quadrature weights for N+1 nodes on [-1, 1].
///
/// wi=2/[(N)(N+1)(P_N(x_i))**2]
fn lgl_weights(n: usize) -> Vec<f64> {
    let x = lgl_nodes(n);
    let (p_n, _, _) = legendre_all(n, &x);
    let w: Vec<f64> = p_n.iter().map(|p| 2.0 / ((n * (n + 1)) as f64 * p * p)).collect();
    // symmetry to avoid round-off errors
    (0..=n).map(|i| 0.5 * (w[i] + w[n - i])).collect()
}

/// Computes the (N+1)x(N+1) LGL differentiation matrix using
/// exact polynomial ratios and the negative-sum trick.
///
/// Derivative of the Lagrange interpolation polynomial passing through points
/// (x1,y1),.....(x_(N+1),y_(N+1)), evaluated at x_1...x_(N+1) given by
///
/// Yd=DY, where Y=[y1....y_(N+1)].T and Yd is the derivative of Y at LGL points.
fn lgl_diff_matrix(n: usize) -> Vec<Vec<f64>> {
    let x = lgl_nodes(n);
    let (p_n, _, _) = legendre_all(n, &x);

    let mut d = vec![vec![0.0; n + 1]; n + 1];
    for i in 0..=n {
        for j in 0..=n {
            if i != j {
                d[i][j] = p_n[i] / p_n[j] / (x[i] - x[j]);
            }
        }
        d[i][i] = -d[i].iter().sum::<f64>();
    }
    let mut ds = d.clone();
    for i in 0..=n {
        for j in 0..=n {
            ds[i][j] = 0.5 * (d[i][j] - d[n - i][n - j]);
        }
    }
    ds
}

/// Convenience wrapper returning LGL nodes (x), quadrature weights (w),
/// and differentiation matrix (D).
fn lgl_setup(n: usize) -> (Vec<f64>, Vec<f64>, ca::DM) {
    (lgl_nodes(n), lgl_weights(n), ca::DM::new(lgl_diff_matrix(n)))
}

struct Data {
    var_analytical: Vec<[f64; 5]>,
    var_numerical: Vec<[f64; 5]>,
    ts: Vec<f64>,
    hamiltonian: ca::DM,
    constraint_jacobian_sparsity: ca::Sparsity,
    state_variable_count: i64,
    defect_constraint_count: i64,
}

/// template function to solve the OCP
///
/// N is number of LGL intervals or degree of polynomial used to approximate the
/// state
///
/// c scaling factor to generate a denser grid
///
/// returns solution and data needed for plots
fn solve_ocp(n: usize, c: usize) -> Data {
    // OCP definition
    // Example 1 of "A Pseudospectral Method for the Optimal Control of
    // Constrained Feedback Linearizable Systems"
    // Qi Gong , Wei Kang, and I. Michael Ross
    // number of states, controls, dynamics, lagrange term, intial state
    // initial and final time
    let nx = 2;
    let nu = 1;
    let x = ca::MX::sym_with("x", nx, 1);
    let u = ca::MX::sym_with("u", nu, 1);
    let xd = ca::Function::new("xd", vec![&x, &u], vec![ca::vertcat(vec![ca::power(&x.at(1), 3.0), u.clone()])]);
    let lag = ca::Function::new("lag", vec![&x, &u], vec![4.0 * ca::power(&u, 2.0)]);
    let x_t0 = ca::DM::new(vec![0.0, 1.0]);
    let t0 = 0.0;
    let tf = 2.0;

    // decision variables
    // N+1 LGL points [-1,1] tau mapped to physical [t0,tf]
    // t=0.5(tf-t0)tau+0.5(tf+t0)
    // state and control at these points decision variables
    let nlp = ca::Opti::default();
    let xx = nlp.variable_with(nx, n as i64 + 1, "full");
    let uu = nlp.variable_with(nu, n as i64 + 1, "full");
    let x0 = xx.at((.., 0));
    let xf = xx.at((.., -1));

    // obtain LGL nodes, quadrature weights, differentiation matrices
    let (tau, wi, d) = lgl_setup(n);
    let wi_row = ca::DM::new(vec![wi.clone()]);

    // define event, defect constraints
    let event = &x0 - &x_t0;
    // differentiation matrix
    let xd_map: ca::MX = xd.map(n as i64 + 1).eval((&xx, &uu));
    let defect: ca::MX = ca::eq(&(ca::mtimes(&xx, &d.T()) - (tf - t0) / 2.0 * xd_map), 0.0); // to extract dual vars later

    // intergral approx. with quadrature
    let mayer = 4.0 * xf.at(0) + xf.at(1);
    let lag_map: ca::MX = lag.map(n as i64 + 1).eval((&xx, &uu));
    let lagrange = 0.5 * (tf - t0) * ca::dot(&ca::MX::new(&wi_row), &lag_map);
    let objective = lagrange + mayer;

    nlp.minimize(&objective);
    nlp.subject_to(&defect);
    nlp.subject_to(&ca::eq(&event, 0.0));

    // ipopt solver and initial guess from linear interpolation
    nlp.solver("ipopt");
    let lin: ca::DM = ca::linspace(0.0, 1.0, n as i64 + 1);
    nlp.set_initial(&xx, &ca::vertcat(vec![lin.T(), lin.T()]));
    let (r, cc) = uu.size();
    nlp.set_initial(&uu, &ca::DM::ones_with(r, cc));
    let sol = nlp.solve();

    // numerical state control state derivative time costates
    let xs = sol.value(&xx);
    let us = ca::reshape_with(&sol.value(&uu), nu, -1);
    let xds = sol.value(&(2.0 / (tf - t0) * ca::mtimes(&xx, &d.T())));
    let ts: Vec<f64> = tau.iter().map(|t| (tf - t0) / 2.0 * t + 0.5 * (tf + t0)).collect();

    // covector mapping principle (-Lagrange multipliers / quadrature weights)
    let adjoint = -sol.value(&(nlp.dual(&defect) / ca::repmat_with(&wi_row, nx, 1)));

    // hamiltonian for the OCP
    let lag_num: ca::DM = lag.map(n as i64 + 1).eval((&xs, &us));
    let hamiltonian = lag_num + ca::sum1(&(&adjoint * &xds));

    // Analytical reference for the local solution computed here.
    let x1a = |t: f64| -64.0 / (5.0 * (2.0 + t).powi(5)) + 2.0 / 5.0;
    let x2a = |t: f64| 4.0 / (2.0 + t).powi(2);
    let ua = |t: f64| -8.0 / (2.0 + t).powi(3);
    let lam1a = |_t: f64| 4.0;
    let lam2a = |t: f64| 64.0 / (2.0 + t).powi(3);
    let td: Vec<f64> = lgl_nodes(n * c).iter().map(|t| (tf - t0) / 2.0 * t + 0.5 * (tf + t0)).collect();
    let var_analytical = td.iter().map(|&t| [x1a(t), x2a(t), ua(t), lam1a(t), lam2a(t)]).collect();
    let var_numerical = (0..=n as i64)
        .map(|k| {
            [
                xs.at((0, k)).to_double(),
                xs.at((1, k)).to_double(),
                us.at((0, k)).to_double(),
                adjoint.at((0, k)).to_double(),
                adjoint.at((1, k)).to_double(),
            ]
        })
        .collect();

    Data {
        var_analytical,
        var_numerical,
        ts,
        hamiltonian,
        constraint_jacobian_sparsity: ca::jacobian(&nlp.g(), &nlp.x()).sparsity(),
        state_variable_count: xx.numel(),
        defect_constraint_count: defect.numel(),
    }
}

/// primal dual solution vs analytical
fn print_solution_verify(data: &Data) {
    let labels = ["x_1", "x_2", "u", "lambda_1", "lambda_2"];
    println!("{:>8} {:>12} {:>12} {:>12} {:>12} {:>12}", "t", labels[0], labels[1], labels[2], labels[3], labels[4]);
    for (t, v) in data.ts.iter().zip(&data.var_numerical) {
        println!("{:>8.4} {:>12.6} {:>12.6} {:>12.6} {:>12.6} {:>12.6}", t, v[0], v[1], v[2], v[3], v[4]);
    }

    // Hamiltonian consistency check (a necessary optimality condition)
    // constant for time-invariant problem
    // zero for this specific problem
    println!("Hamiltonian: {}", data.hamiltonian);
}

/// Show the dense state coupling from global polynomial differentiation.
fn print_jacobian_sparsity(data: &Data) {
    let sparsity = &data.constraint_jacobian_sparsity;
    // Opti stacks X before U and the dynamics before the initial conditions.
    println!(
        "Constraint Jacobian sparsity: {}x{}, nnz {} (states: {}, dynamics constraints: {})",
        sparsity.size1(),
        sparsity.size2(),
        sparsity.nnz(),
        data.state_variable_count,
        data.defect_constraint_count
    );
}

/// solve the given OCP for increasing grid sizes. A sudden dip in
/// error with small grid sizes indicate spectral convergence
fn spectral_convergence() {
    let (nmin, nmax) = (5, 50);
    for i in (nmin..nmax).step_by(5) {
        let data = solve_ocp(i, 1);
        let mut err = [0.0f64; 5];
        for (a, b) in data.var_analytical.iter().zip(&data.var_numerical) {
            for k in 0..5 {
                err[k] = err[k].max((a[k] - b[k]).abs());
            }
        }
        println!("N = {:2}: max. absolute error {:?}", i, err);
    }
}

fn main() {
    let data = solve_ocp(25, 3);
    print_solution_verify(&data);
    print_jacobian_sparsity(&data);
    spectral_convergence();
}

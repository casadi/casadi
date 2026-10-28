//
//     MIT No Attribution
//
//     Copyright (C) 2010-2026 Joel Andersson, Joris Gillis, Moritz Diehl, KU Leuven.
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
// Solve the hs015 problem, formulated as the NLP:
//
// minimize     100*(x[1] - x[0]**2)**2 + (1 - x[0])**2
// subject to   x[0]*x[1] >= 1
//              x[0] + x[1]**2 >= 0
//              x[0] <= 1/2
//
// David Kiessling, 2026

use casadi as ca;
use casadi::dict;

const INF: f64 = f64::INFINITY;

fn main() {
    // Declare variable x
    let x = ca::MX::sym_with("x", 2, None);

    // Formulate the NLP
    let f = 100.0 * ca::power(x.at(1) - ca::power(&x.at(0), 2.0), 2.0) + ca::power(1.0 - x.at(0), 2.0);
    let mut g = ca::MX::zeros_with(2, 1);
    g.set_at(0, x.at(0) * x.at(1));
    g.set_at(1, x.at(0) + ca::power(&x.at(1), 2.0));

    let ubx = ca::DM::new(vec![1.0 / 2.0, INF]);
    let _lbx = ca::DM::new(vec![-INF, -INF]);
    let lbg = ca::DM::new(vec![1.0, 0.0]);
    let ubg = ca::DM::new(vec![INF, INF]);
    let nlp = dict! {"x" => &x, "f" => &f, "g" => &g};

    // Create an NLP solver
    let solver = ca::nlpsol_with("solver", "uno", &nlp, &dict! {"uno" => dict! {"preset" => "ipopt"}});

    // Solve the problem
    let res = solver
        .call(&dict! {"x0" => vec![-2.0, 1.0], "ubg" => &ubg, "lbg" => &lbg, "ubx" => &ubx})
        .unwrap();

    // Print solution
    println!();
    println!("{:>50}  {}", "Optimal cost:", res["f"]);
    println!("{:>50}  {}", "Primal solution:", res["x"]);
    println!("{:>50}  {}", "Dual solution (simple bounds):", res["lam_x"]);
    println!("{:>50}  {}", "Dual solution (nonlinear bounds):", res["lam_g"]);
}

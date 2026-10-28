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
// Solve the Rosenbrock problem, formulated as the NLP:
//   minimize     x^2 + 100*z^2
//   subject to   z + (1-x)^2 - y == 0

use casadi as ca;
use casadi::dict;

fn main() {
    // Declare variables
    let x = ca::SX::sym("x");
    let y = ca::SX::sym("y");
    let z = ca::SX::sym("z");

    // Formulate the NLP
    let f = ca::power(&x, 2.0) + 100.0 * ca::power(&z, 2.0);
    let g = &z + ca::power(1.0 - &x, 2.0) - &y;
    let nlp = dict! {"x" => ca::vertcat(vec![&x, &y, &z]), "f" => &f, "g" => &g};

    // Create an NLP solver
    let solver = ca::nlpsol("solver", "ipopt", &nlp);

    // Solve the Rosenbrock problem
    let res = solver.call(&dict! {"x0" => vec![2.5, 3.0, 0.75], "ubg" => 0.0, "lbg" => 0.0}).unwrap();

    // Print solution
    println!();
    println!("{:>50} {}", "Optimal cost:", res["f"]);
    println!("{:>50} {}", "Primal solution:", res["x"]);
    println!("{:>50} {}", "Dual solution (simple bounds):", res["lam_x"]);
    println!("{:>50} {}", "Dual solution (nonlinear bounds):", res["lam_g"]);
}

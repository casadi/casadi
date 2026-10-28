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
    // Control
    let u = ca::MX::sym("u");

    // State
    let x = ca::MX::sym_with("x", 3, 1);
    let _s = x.at(0); // position
    let v = x.at(1); // speed
    let m = x.at(2); // mass

    // ODE right hand side
    let sdot = v.clone();
    let vdot = (&u - 0.05 * &v * &v) / &m;
    let mdot = -0.1 * &u * &u;
    let xdot = ca::vertcat(vec![&sdot, &vdot, &mdot]);

    // ODE right hand side function
    let f = ca::Function::new("f", vec![&x, &u], vec![&xdot]);

    // Integrate with Explicit Euler over 0.2 seconds
    let dt = 0.01; // Time step
    let mut xj = x.clone();
    for _j in 0..20 {
        let fj: ca::MX = f.eval((&xj, &u));
        xj += dt * &fj;
    }

    // Discrete time dynamics function
    let ff = ca::Function::new("F", vec![&x, &u], vec![&xj]);

    // Number of control segments
    let nu = 50;

    // Control for all segments
    let uu = ca::MX::sym_with("U", nu, 1);

    // Initial conditions
    let x0 = ca::MX::new(vec![0.0, 0.0, 1.0]);

    // Integrate over all intervals
    let mut xx = x0;
    for k in 0..nu {
        xx = ff.eval((&xx, &uu.at(k)));
    }

    // Objective function and constraints
    let j = ca::mtimes(&uu.T(), &uu); // u'*u in Matlab
    let g = xx.at(0..2); // x(1:2) in Matlab

    // NLP
    let nlp = dict! {"x" => &uu, "f" => &j, "g" => &g};

    // Allocate an NLP solver
    let opts = dict! {"ipopt.tol" => 1e-10, "expand" => true};
    let solver = ca::nlpsol_with("solver", "ipopt", &nlp, &opts);

    // Bounds on u and initial condition, bounds on g
    let arg = dict! {"lbx" => -0.5, "ubx" => 0.5, "x0" => 0.4, "lbg" => vec![10.0, 0.0], "ubg" => vec![10.0, 0.0]};

    // Solve the problem
    let res = solver.call(&arg).unwrap();

    // Get the solution (Python plots x and lam_x)
    println!("f: {}", res["f"]);
    println!("x: {}", res["x"]);
    println!("lam_x: {}", res["lam_x"]);
}

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
    // Consider the following nonlinear system with two states and one control input:
    let x1 = ca::MX::sym("x1");
    let x2 = ca::MX::sym("x2");
    let x = ca::vertcat(vec![&x1, &x2]);
    let u = ca::MX::sym("u");
    let xdot = ca::vertcat(vec![(1.0 - ca::power(&x2, 2.0)) * &x1 - &x2 + &u, x1.clone()]);

    // In addition to the dynamics, we are also interested in the sum-of-squares distance from the origin:
    let y = ca::power(&x1, 2.0) + ca::power(&x2, 2.0);

    // We want to simulate this system over T seconds, with the control discretized into N intervals:
    // Time horizon, discretization
    let n = 10;
    let t_end = 10.0;
    let tgrid: ca::DM = ca::linspace(0.0, t_end, n + 1);
    let tgrid = tgrid.nonzeros();

    // We will integrate the problem using CVODES, while also calculating the integral of y
    let dae = dict! {"x" => &x, "u" => &u, "ode" => &xdot, "quad" => &y};
    let f = ca::integrator_with("F", "cvodes", &dae, tgrid[0], tgrid[1..].to_vec(), None);

    // We are also interested in the pointwise values of y, so let's create a function for that too
    let yfun = ca::Function::named("yfun", vec![&x], vec![&y], vec!["x"], vec!["y"]);

    // Let's define initial conditions for x and values for the piecewise constant controls
    let u: ca::DM = ca::linspace(-1.0, 1.0, n);
    let u = u.nonzeros();
    let x0 = vec![0.0, 0.0];

    // We can simulate the system as follows:
    let fk = f.call(&dict! {"x0" => x0.clone(), "u" => u.clone()}).unwrap();
    let xf = &fk["xf"];
    let qf = &fk["qf"];

    // Let's add the state at the initial time to xf and qf
    let xf = ca::horzcat(vec![&ca::DM::new(x0), xf]);
    let qf = ca::horzcat(vec![&ca::DM::new(0.0), qf]);

    // Call yfun to get y at each grid point
    let yf: ca::DM = yfun.eval([&xf]);

    // Print simulation trajectories and pointwise sum of squares contribution to y, which we compare with the integral of y
    let y_sim = yf.nonzeros();
    let cum: Vec<f64> = y_sim.iter().scan(0.0, |s, v| { *s += t_end / n as f64 * v; Some(*s) }).collect();
    println!("tgrid: {:?}", tgrid);
    println!("u: {:?}", u);
    println!("x_sim: {}", xf);
    println!("y_sim: {}", yf);
    println!("q_sim: {}", qf);
    println!("T/N*cumsum(y_sim): {:?}", cum);
}

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
    // Height and velocity of ball
    let h = ca::SX::sym("h");
    let v = ca::SX::sym("v");
    let x = ca::vertcat(vec![&h, &v]);

    // ODE right-hand-side
    let hdot = v.clone();
    let vdot = ca::SX::new(-9.81);
    let xdot = ca::vertcat(vec![&hdot, &vdot]);

    // Event indicator, trigger when it becomes negative
    let event_indicator = &h;

    // DAE problem structure, with zero-crossing output
    let dae = dict! {"x" => &x, "ode" => &xdot, "zero" => event_indicator};

    // Event transition function
    let post_x = ca::vertcat(vec![h.clone(), -0.8 * &v]);
    let transition = ca::Function::from_dict("transition", dict! {"x" => &x, "post_x" => &post_x},
                                             ca::event_in(), ca::event_out());

    // Create an integrator instance for integrating over 7s
    let tgrid: ca::DM = ca::linspace(0.0, 7.0, 100);
    let tgrid = tgrid.nonzeros();
    let sim = ca::integrator_with("sim", "cvodes", &dae, 0.0, &tgrid,
                                  dict! {"transition" => &transition});

    // Simulate with initial height of 5
    let x0 = vec![5.0, 0.0];
    let simres = sim.call(&dict! {"x0" => x0}).unwrap();

    // Visualize the solution: print the height instead
    let hsol = simres["xf"].at((0, ..)).nonzeros();
    for k in (0..tgrid.len()).step_by(11) {
        println!("t = {:.4}: h = {:.6}", tgrid[k], hsol[k]);
    }
}

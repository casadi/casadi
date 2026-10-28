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

// Simulating and optimizing a spring that breaks with DaeBuilder and event handling
// Joel Andersson, 2024-2025

use casadi as ca;
use casadi::dict;

fn main() {
    // Start with an empty DaeBuilder instance
    let dae = ca::DaeBuilder::new("breaking_spring");

    // Model variables
    let m = dae.add_var_with("m", "parameter", "tunable", dict! {"start" => 1, "description" => "Mass"});
    let v = dae.add_3("v", "local", dict! {"start" => -5, "description" => "Velocity"});
    let x = dae.add_3("x", "local", dict! {"start" => -1, "description" => "Displacement"});
    let k = dae.add_var_with("k", "parameter", "tunable", dict! {"start" => 2, "description" => "Spring constant"});
    let c = dae.add_var_with("c", "parameter", "tunable", dict! {"start" => 0.1, "description" => "Damping constant"});
    let d = dae.add_3("d", "input", dict! {"start" => 0, "description" => "Disturbance"});
    let f = dae.add_3("f", "local", dict! {"description" => "Spring force"});
    let b = dae.add_3("b", "local", dict! {"type" => "Boolean", "start" => false, "description" => "Is the spring broken?"});

    // Equations
    dae.eq(&dae.der_with(&x), &v);
    let spring: ca::MX = ca::if_else(&ca::not(&b), &(-&k * &x + &d), 0.0);
    dae.eq(&f, &spring);
    dae.eq(&dae.der_with(&v), &((&f - &c * &v) / &m));
    dae.when(&ca::gt(&x, 2.0), vec![dae.assign("b", 1.0)]);
    dae.disp_with(true);

    // Create an integrator in CasADi for hybrid simulation
    let tgrid: ca::DM = ca::linspace(0.0, 4.0, 101);
    let tgrid = tgrid.nonzeros();
    let simopts = dict! {"transition" => dae.transition(), "verbose" => false, "event_tol" => 1e-12};
    let sim = ca::integrator_with("sim", "cvodes", &dae.create(), 0.0, &tgrid, simopts);

    // Disturbance trajectory
    let mut d_traj = vec![0.0; tgrid.len()];
    d_traj[20..30].fill(-5.0);
    d_traj[40..50].fill(15.0);

    // Simulate the system
    let res = sim.call(&dict! {"x0" => dae.start(dae.x()), "p" => dae.start(dae.p()), "u" => &d_traj}).unwrap();

    // Visualize the solution: print the states instead
    println!("states: {:?}", dae.x());
    let xf = &res["xf"];
    for k in (0..tgrid.len()).step_by(10) {
        println!("t = {:.2}: {:?}", tgrid[k], xf.at((.., k)).nonzeros());
    }
    println!("disturbance: {}", dae.u()[0]);
}

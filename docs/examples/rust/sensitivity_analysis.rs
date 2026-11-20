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

// Integrate the DAE with each integrator and compare sensitivities
fn test_integrators(integrators: &[&str], dae: &ca::rt::Dict, tf: f64, x0: &ca::DM, u0: f64) {
    // Integrator
    for my_integrator in integrators {
        println!("========");
        println!("Integrator:  {}", my_integrator);
        println!("========");

        // Integrator options
        let opts = if *my_integrator == "collocation" { dict! {"rootfinder" => "kinsol"} } else { dict! {} };

        // Integrator
        let i = ca::integrator_with("I", my_integrator, dae, 0.0, tf, &opts);

        // Integrate to get results
        let res = i.call(&dict! {"x0" => x0, "p" => u0}).unwrap();
        let xf = &res["xf"];
        let qf = &res["qf"];
        println!("{:>50} xf  =  {} , qf  =  {}", "Unperturbed solution:", xf, qf);

        // Perturb solution to get a finite difference approximation
        let h = 0.001;
        let res = i.call(&dict! {"x0" => x0, "p" => u0 + h}).unwrap();
        let fd_xf = (&res["xf"] - xf) / h;
        let fd_qf = (&res["qf"] - qf) / h;
        println!("{:>50} d(xf)/d(p) =  {} , d(qf)/d(p) =  {}", "Finite difference approximation:", fd_xf, fd_qf);

        // Calculate once, forward
        let i_fwd = i.factory("I_fwd", vec!["x0", "z0", "p", "fwd:p"], vec!["fwd:xf", "fwd:qf"]);
        let res = i_fwd.call(&dict! {"x0" => x0, "p" => u0, "fwd_p" => 1.0}).unwrap();
        let fwd_xf = &res["fwd_xf"];
        let fwd_qf = &res["fwd_qf"];
        println!("{:>50} d(xf)/d(p) =  {} , d(qf)/d(p) =  {}", "Forward sensitivities:", fwd_xf, fwd_qf);

        // Calculate once, adjoint
        let i_adj = i.factory("I_adj", vec!["x0", "z0", "p", "adj:qf"], vec!["adj:x0", "adj:p"]);
        let res = i_adj.call(&dict! {"x0" => x0, "p" => u0, "adj_qf" => 1.0}).unwrap();
        let adj_x0 = res["adj_x0"].clone();
        let adj_p = res["adj_p"].clone();
        println!("{:>50} d(qf)/d(x0) =  {} , d(qf)/d(p) =  {}", "Adjoint sensitivities:", adj_x0, adj_p);

        // Perturb adjoint solution to get a finite difference approximation of the second order sensitivities
        let res = i_adj.call(&dict! {"x0" => x0, "p" => u0 + h, "adj_qf" => 1.0}).unwrap();
        let fd_adj_x0 = (&res["adj_x0"] - &adj_x0) / h;
        let fd_adj_p = (&res["adj_p"] - &adj_p) / h;
        println!("{:>50} d2(qf)/d(x0)d(p) =  {} , d2(qf)/d(p)d(p) =  {}", "FD of adjoint sensitivities:", fd_adj_x0, fd_adj_p);

        // Forward over adjoint to get the second order sensitivities
        let i_foa = i_adj.factory("I_foa", vec!["x0", "z0", "p", "adj_qf", "fwd:p"], vec!["fwd:adj_x0", "fwd:adj_p"]);
        let res = i_foa.call(&dict! {"x0" => x0, "p" => u0, "adj_qf" => 1.0, "fwd_p" => 1.0}).unwrap();
        let fwd_adj_x0 = &res["fwd_adj_x0"];
        let fwd_adj_p = &res["fwd_adj_p"];
        println!("{:>50} d2(qf)/d(x0)d(p) =  {} , d2(qf)/d(p)d(p) =  {}", "Forward over adjoint sensitivities:", fwd_adj_x0, fwd_adj_p);

        // Adjoint over adjoint to get the second order sensitivities
        let i_aoa = i_adj.factory("I_aoa", vec!["x0", "z0", "p", "adj_qf", "adj:adj_p"], vec!["adj:x0", "adj:p"]);
        let res = i_aoa.call(&dict! {"x0" => x0, "p" => u0, "adj_qf" => 1.0, "adj_adj_p" => 1.0}).unwrap();
        let adj_x0 = &res["adj_x0"];
        let adj_p = &res["adj_p"];
        println!("{:>50} d2(qf)/d(x0)d(p) =  {} , d2(qf)/d(p)d(p) =  {}", "Adjoint over adjoint sensitivities:", adj_x0, adj_p);
    }
}

fn main() {
    println!("Testing sensitivity analysis in CasADi");

    // All ODE and DAE integrators to be tested
    let dae_integrators = ["idas", "collocation"];
    let ode_integrators = ["cvodes", "rk", "idas", "collocation"];

    {
        // rocket example
        println!("******");
        println!("Testing ODE example");

        // Time
        let t = ca::SX::sym("t");

        // Parameter
        let u = ca::SX::sym("u");

        // Differential states
        let s = ca::SX::sym("s");
        let v = ca::SX::sym("v");
        let m = ca::SX::sym("m");
        let x = ca::vertcat(vec![&s, &v, &m]);

        // Constants
        let alpha = 0.05; // friction
        let beta = 0.1; // fuel consumption rate

        // Differential equation
        let ode = ca::vertcat(vec![v.clone(), (&u - alpha * &v * &v) / &m, -beta * &u * &u]);

        // Quadrature
        let quad = ca::power(&v, 3.0) + ca::power(&((3.0 - ca::sin(&t)) - &u), 2.0);

        // DAE callback function
        let dae = dict! {"t" => &t, "x" => &x, "p" => &u, "ode" => &ode, "quad" => &quad};

        // Time length
        let tf = 0.5;

        // Initial position
        let x0 = ca::DM::new(vec![0.0, 0.0, 1.0]);

        // Parameter
        let u0 = 0.4;

        test_integrators(&ode_integrators, &dae, tf, &x0, u0);
    }

    {
        // Simple DAE example
        println!("******");
        println!("Testing DAE example");

        // Differential state
        let x = ca::SX::sym("x");

        // Algebraic variable
        let z = ca::SX::sym("z");

        // Parameter
        let u = ca::SX::sym("u");

        // Differential equation
        let ode = -&x + 0.5 * &x * &x + &u + 0.5 * &z;

        // Algebraic constraint
        let alg = &z + ca::exp(&z) - 1.0 + &x;

        // Quadrature
        let quad = &x * &x + 3.0 * &u * &u;

        // DAE callback function
        let dae = dict! {"x" => &x, "z" => &z, "p" => &u, "ode" => &ode, "alg" => &alg, "quad" => &quad};

        // End time
        let tf = 5.0;

        // Initial position
        let x0 = ca::DM::new(1.0);

        // Parameter
        let u0 = 0.4;

        test_integrators(&dae_integrators, &dae, tf, &x0, u0);
    }
}

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
    // Declare variables (use simple, efficient DAG)
    let x0 = ca::SX::sym("x0");
    let x1 = ca::SX::sym("x1");
    let x = ca::vertcat(vec![&x0, &x1]);

    // Control
    let u = ca::SX::sym("u");

    // ODE right hand side
    let xdot = ca::vertcat(vec![(1.0 - &x1 * &x1) * &x0 - &x1 + &u, x0.clone()]);

    // Lagrangian function
    let l = &x0 * &x0 + &x1 * &x1 + &u * &u;

    // Costate
    let lam = ca::SX::sym_with("lam", 2, 1);

    // Hamiltonian function
    let h = ca::dot(&lam, &xdot) + &l;

    // Costate equations
    let ldot = -ca::gradient(&h, &x);

    // The control must minimize the Hamiltonian, which is:
    println!("Hamiltonian:  {}", h);

    // H is of a convex quadratic form in u: H = u*u + p*u + q, let's get the coefficient p
    let p = ca::gradient(&h, &u); // this gives us 2*u + p
    let p = ca::substitute(&p, &u, 0.0); // replace u with zero: gives us p

    // H's unconstrained minimizer is: u = -p/2
    let u_opt = -&p / 2.0;

    // We must constrain u to the interval [-0.75, 1.0], convexity of H ensures that the optimum is obtain at the bound when u_opt is outside the interval
    let u_opt = ca::fmin(&u_opt, 1.0);
    let u_opt = ca::fmax(&u_opt, -0.75);
    println!("optimal control:  {}", u_opt);

    // Augment f with lam_dot and substitute in the value for the optimal control
    let f = ca::vertcat(vec![&xdot, &ldot]);
    let f = ca::substitute(&f, &u, &u_opt);

    // Function for obtaining the optimal control from the augmented state
    let u_fcn = ca::Function::new("ufcn", vec![ca::vertcat(vec![&x, &lam])], vec![&u_opt]);

    // Formulate the DAE
    let dae = dict! {"x" => ca::vertcat(vec![&x, &lam]), "ode" => &f};

    // Create an integrator (CVodes)
    let opts = dict! {
        "abstol" => 1e-8, // abs. tolerance
        "reltol" => 1e-8, // rel. tolerance
    };
    let i = ca::integrator_with("I", "cvodes", &dae, 0.0, 10.0, &opts);

    // The initial state
    let x_init = vec![0.0, 1.0];

    // The initial costate
    let l_init = ca::MX::sym_with("l_init", 2, 1);

    // The initial condition for the shooting
    let xx = ca::vertcat(vec![&ca::MX::new(ca::DM::new(x_init.clone())), &l_init]);

    // Call the integrator
    let xx = i.call_symbolic(&std::collections::HashMap::from([("x0", &xx)]))["xf"].clone();

    // Costate at the final time should be zero (cf. Bryson and Ho)
    let lam_f = xx.at(2..4);
    let g = lam_f;

    // Formulate root-finding problem
    let rfp = ca::Function::new("rfp", vec![&l_init], vec![&g]);

    // Select a solver for the root-finding problem
    let solver_name = "nlpsol";
    //let solver_name = "newton";
    //let solver_name = "kinsol";

    // Allocate an implict solver
    let opts = match solver_name {
        "nlpsol" => dict! {"nlpsol" => "ipopt", "nlpsol_options" => dict! {"ipopt.hessian_approximation" => "limited-memory"}},
        "kinsol" => dict! {"linear_solver_type" => "user_defined", "max_iter" => 1000},
        _ => dict! {},
    };
    let solver = ca::rootfinder_with("solver", solver_name, &rfp, &opts);

    // Solve the problem
    let l_init_opt: ca::DM = solver.eval([0.0]);
    let l_init_opt = l_init_opt.nonzeros();

    // Time grid for visualization
    let tgrid: ca::DM = ca::linspace(0.0, 10.0, 100);
    let tgrid = tgrid.nonzeros();

    // Simulator to get optimal state and control trajectories
    let simulator = ca::integrator_with("simulator", "cvodes", &dae, 0.0, tgrid.clone(), dict! {});

    // Simulate to get the state trajectory
    let x0_sim: Vec<f64> = x_init.iter().chain(l_init_opt.iter()).copied().collect();
    let sol = simulator.call(&dict! {"x0" => x0_sim}).unwrap()["xf"].clone();

    // Calculate the optimal control
    let ufcn_all = u_fcn.map(tgrid.len() as i64);
    let u_opt: ca::DM = ufcn_all.eval([&sol]);

    // Print the results
    println!("l_init_opt: {:?}", l_init_opt);
    println!("x trajectory: {}", sol.at((0, ..)));
    println!("y trajectory: {}", sol.at((1, ..)));
    println!("u trajectory: {}", u_opt);
}

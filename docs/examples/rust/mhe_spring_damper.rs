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
// For documentation about this examples, check http://docs.casadi.org/documents/mhe_spring_damper.pdf
//
// The Python original uses casadi.tools structures (struct_symSX) for the
// horizon; here the decision variables and parameters are plain matrices,
// concatenated with vertcat/vec and split again with reshape.

use casadi as ca;
use casadi::dict;

// Normally distributed random numbers (stand in for numpy.random.randn)
struct Rng(u64);
impl Rng {
    fn uniform(&mut self) -> f64 {
        self.0 = self.0.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        ((self.0 >> 11) as f64 + 0.5) / (1u64 << 53) as f64
    }
    fn randn(&mut self, n: usize, m: usize) -> ca::DM {
        let mut rows = vec![vec![0.0; m]; n];
        for r in rows.iter_mut() {
            for v in r.iter_mut() {
                *v = (-2.0 * self.uniform().ln()).sqrt() * (2.0 * std::f64::consts::PI * self.uniform()).cos();
            }
        }
        ca::DM::new(rows)
    }
}

fn main() {
    let mut rng = Rng(0);

    // Settings of the filter
    let n = 10; // Horizon length
    let dt = 0.05; // Time step

    let sigma_p = 0.005; // Standard deviation of the position measurements
    let sigma_w = 0.1; // Standard deviation for the process noise
    let r = ca::DM::new(1.0 / (sigma_p * sigma_p)); // resulting weighting matrix for the position measurements
    let q = ca::DM::new(1.0 / (sigma_w * sigma_w)); // resulting weighting matrix for the process noise

    let nsimulation = 1000; // Lenght of the simulation

    // Parameters of the system
    let m = 1.0; // The weight of the mass
    let k = 1.0; // The spring constant
    let c = 0.5; // The damping of the system
    // The state
    let x = ca::SX::sym("x");
    let dx = ca::SX::sym("dx");
    let states = ca::vertcat(vec![&x, &dx]); // Full state vector of the system: position x and velocity dx
    let nstates = states.nnz(); // Number of states

    // The control input
    let f_in = ca::SX::sym("F"); // Full control vector of the system: Input force F
    let controls = f_in.clone();

    // Disturbances
    let w = ca::SX::sym("w"); // Process noise vector
    let disturbances = w.clone();
    let ndisturbances = disturbances.nnz(); // Number of disturbances

    // Measurements
    let nmeas = 1; // Number of measurements: y

    // Create Structure for the entire horizon

    // Degrees of freedom for the optimizer
    let xs = ca::SX::sym_with("X", nstates, n);
    let ws = ca::SX::sym_with("W", ndisturbances, n - 1);
    let shooting = ca::vertcat(vec![ca::vec(&xs), ca::vec(&ws)]);
    // Fixed parameters for the optimizer
    let us = ca::SX::sym_with("U", 1, n - 1);
    let ys = ca::SX::sym_with("Y", nmeas, n);
    let s = ca::SX::sym_with("S", nstates, nstates);
    let x0s = ca::SX::sym_with("x0", nstates, 1);
    let parameters = ca::vertcat(vec![ca::vec(&us), ca::vec(&ys), ca::vec(&s), x0s.clone()]);
    // Numeric counterparts
    let params_num = |u: &ca::DM, y: &ca::DM, s: &ca::DM, x0: &ca::DM| ca::vertcat(vec![ca::vec(u), ca::vec(y), ca::vec(s), x0.clone()]);
    let split_shooting = |v: &ca::DM| {
        let nx = (nstates * n) as i64;
        (ca::reshape(&v.at(..nx), nstates, n), ca::reshape(&v.at(nx..), ndisturbances, n - 1))
    };

    // Define the ODE right hand side
    let rhs = ca::vertcat(vec![dx.clone(), (-k * &x - c * &dx + &f_in) / m + &w]);

    let f = ca::Function::new("f", vec![&states, &controls, &disturbances], vec![&rhs]);

    // Build an integrator for this system: Runge Kutta 4 integrator
    let k1: ca::SX = f.eval((&states, &controls, &disturbances));
    let k2: ca::SX = f.eval((&states + dt / 2.0 * &k1, &controls, &disturbances));
    let k3: ca::SX = f.eval((&states + dt / 2.0 * &k2, &controls, &disturbances));
    let k4: ca::SX = f.eval((&states + dt * &k3, &controls, &disturbances));

    let states_1 = &states + dt / 6.0 * (&k1 + 2.0 * &k2 + 2.0 * &k3 + &k4);
    let phi = ca::Function::named("phi", vec![&states, &controls, &disturbances], vec![&states_1], vec!["x", "u", "d"], vec!["x1"]);
    let phi_jac = phi.factory("PHI", vec!["x", "u", "d"], vec!["jac:x1:x"]);
    // Define the measurement system
    let h = ca::Function::named("h", vec![&states], vec![&x], vec!["x"], vec!["y"]); // We have measurements of the position
    let h_jac = h.factory("H", vec!["x"], vec!["jac:y:x"]);
    // Build the objective
    let mut obj = ca::SX::new(0.0);
    // First the arrival cost
    let ex0 = xs.at((.., 0)) - &x0s;
    obj += ca::mtimes_list(vec![ex0.T(), s.clone(), ex0.clone()]);
    //Next the cost for the measurement noise
    for i in 0..n {
        let vm: ca::SX = h.eval([xs.at((.., i))]);
        let vm = vm - ys.at((.., i));
        obj += ca::mtimes(&ca::mtimes(&vm.T(), &r), &vm);
    }
    //And also the cost for the process noise
    for i in 0..n - 1 {
        obj += ca::mtimes(&ca::mtimes(&ws.at((.., i)).T(), &q), &ws.at((.., i)));
    }

    // Build the multiple shooting constraints
    let mut g: Vec<ca::SX> = vec![];
    for i in 0..n - 1 {
        let xn: ca::SX = phi.eval((xs.at((.., i)), us.at((.., i)), ws.at((.., i))));
        g.push(xs.at((.., i + 1)) - xn);
    }

    // Formulate the NLP
    let nlp = dict! {"x" => &shooting, "p" => &parameters, "f" => &obj, "g" => ca::vertcat(g)};

    // Make a simulation to create the data for the problem
    let mut simulated_x = ca::DM::zeros_with(nstates, nsimulation);
    simulated_x.set_at((.., 0), ca::DM::new(vec![1.0, 0.0])); // Initial state
    let t: Vec<f64> = (0..nsimulation).map(|i| i as f64 * dt).collect(); // Time grid
    let mut simulated_u = ca::DM::new(t[..t.len() - 1].iter().map(|ti| ti.cos()).collect::<Vec<_>>()).T(); // control input for the simulation
    simulated_u.set_at((.., nsimulation / 2..), 0.0);
    let simulated_w = sigma_w * rng.randn(ndisturbances as usize, nsimulation as usize - 1); // Process noise for the simulation
    for i in 0..nsimulation - 1 {
        let xn: ca::DM = phi.eval((simulated_x.at((.., i)), simulated_u.at((.., i)), simulated_w.at((.., i))));
        simulated_x.set_at((.., i + 1), xn);
    }
    //Create the measurements from these states
    let mut simulated_y = ca::DM::zeros_with(nmeas, nsimulation); // Holder for the measurements
    for i in 0..nsimulation {
        let yi: ca::DM = h.eval([simulated_x.at((.., i))]);
        simulated_y.set_at((.., i), yi);
    }
    // Add noise the the position measurements
    simulated_y += sigma_p * rng.randn(nmeas as usize, nsimulation as usize);

    //The initial estimate and related covariance, which will be used for the arrival cost
    let sigma_x0 = 0.01;
    let mut p = sigma_x0 * sigma_x0 * ca::DM::eye(nstates);
    let mut x0 = simulated_x.at((.., 0)) + sigma_x0 * rng.randn(nstates as usize, 1);
    // Create the solver
    let opts = dict! {"ipopt.print_level" => 0, "print_time" => false, "ipopt.max_iter" => 100};
    let nlpsol = ca::nlpsol_with("nlpsol", "ipopt", &nlp, &opts);

    // Create a holder for the estimated states and disturbances
    let mut estimated_x = ca::DM::zeros_with(nstates, nsimulation);
    let mut estimated_w = ca::DM::zeros_with(ndisturbances, nsimulation - 1);
    // For the first instance we run the filter, we need to initialize it.
    let mut cur_u = simulated_u.at((.., 0..n - 1));
    let mut cur_y = simulated_y.at((.., 0..n));
    let mut cur_s = ca::inv(&p); // Arrival cost is the inverse of the initial covariance
    let mut init_x = simulated_x.at((.., 0..n));
    let mut init_w = ca::DM::zeros_with(ndisturbances, n - 1);
    let res = nlpsol
        .call(&dict! {"p" => params_num(&cur_u, &cur_y, &cur_s, &x0), "x0" => ca::vertcat(vec![ca::vec(&init_x), ca::vec(&init_w)]), "lbg" => 0.0, "ubg" => 0.0})
        .unwrap();

    // Get the solution
    let (mut sol_x, mut sol_w) = split_shooting(&res["x"]);
    estimated_x.set_at((.., 0..n), &sol_x);
    estimated_w.set_at((.., 0..n - 1), &sol_w);

    // Now make a loop for the rest of the simulation
    for i in 1..nsimulation - n + 1 {
        // Update the arrival cost, using linearisations around the estimate of MHE at the beginning of the horizon (according to the 'Smoothed EKF Update'): first update the state and covariance with the measurement that will be deleted, and next propagate the state and covariance because of the shifting of the horizon
        let status = match &nlpsol.stats()["return_status"] {
            ca::rt::Value::Str(s) => s.clone(),
            v => format!("{:?}", v),
        };
        println!("step {}/{} ({})", i, nsimulation - n, status);
        let h0_jac: ca::DM = h_jac.eval([sol_x.at((.., 0))]);
        let kk = ca::mtimes_list(vec![p.clone(), h0_jac.T(), ca::inv(&(ca::mtimes_list(vec![h0_jac.clone(), p.clone(), h0_jac.T()]) + &r))]);
        p = ca::mtimes(&(ca::DM::eye(nstates) - ca::mtimes(&kk, &h0_jac)), &p);
        let h0: ca::DM = h.eval([sol_x.at((.., 0))]);
        x0 = &x0 + ca::mtimes(&kk, &(cur_y.at((.., 0)) - &h0 - ca::mtimes(&h0_jac, &(&x0 - sol_x.at((.., 0))))));
        x0 = phi.eval((&x0, cur_u.at((.., 0)), sol_w.at((.., 0))));
        let ff: ca::DM = phi_jac.eval((sol_x.at((.., 0)), cur_u.at((.., 0)), sol_w.at((.., 0))));
        p = ca::mtimes_list(vec![ff.clone(), p.clone(), ff.T()]) + ca::inv(&q);
        // Get the measurements and control inputs
        cur_u = simulated_u.at((.., i..i + n - 1));
        cur_y = simulated_y.at((.., i..i + n));
        cur_s = ca::inv(&p);
        // Initialize the system with the shifted solution
        init_w.set_at((.., 0..n - 2), estimated_w.at((.., i..i + n - 2))); // The shifted solution for the disturbances
        init_w.set_at((.., n - 2), ca::DM::zeros_with(ndisturbances, 1)); // The last node for the disturbances is initialized with zeros
        init_x.set_at((.., 0..n - 1), estimated_x.at((.., i..i + n - 1))); // The shifted solution for the state estimates
        // The last node for the state is initialized with a forward simulation
        let phi0: ca::DM = phi.eval((init_x.at((.., n - 1)), cur_u.at((.., -1)), init_w.at((.., -1))));
        init_x.set_at((.., n - 1), phi0);
        // And now initialize the solver and solve the problem
        let res = nlpsol
            .call(&dict! {"p" => params_num(&cur_u, &cur_y, &cur_s, &x0), "x0" => ca::vertcat(vec![ca::vec(&init_x), ca::vec(&init_w)]), "lbg" => 0.0, "ubg" => 0.0})
            .unwrap();
        (sol_x, sol_w) = split_shooting(&res["x"]);

        // Now get the state estimate. Note that we are only interested in the last node of the horizon
        estimated_x.set_at((.., n - 1 + i), sol_x.at((.., n - 1)));
        estimated_w.set_at((.., n - 2 + i), sol_w.at((.., n - 2)));
    }
    // Plot the results (printed instead)
    println!("Estimated position (last 5): {}", estimated_x.at((0, nsimulation - 5..)));
    println!("Real position (last 5):      {}", simulated_x.at((0, nsimulation - 5..)));

    let error = estimated_x.at((0, ..)) - simulated_x.at((0, ..));
    let e2 = ca::mtimes(&error, &error.T());
    println!("{}", e2);
    assert!(e2.to_double() < 0.01);
}

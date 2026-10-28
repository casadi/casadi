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
// In this example, we fit a nonlinear model to measurements
//
// This example uses more advanced constructs than the vdp* examples:
// Since the number of control intervals is potentially very large here,
// we use memory-efficient Map and MapAccum, in combination with
// codegeneration.
//
// We will be working with a 2-norm objective:
// || y_measured - y_simulated ||_2^2
//
// This form is well-suited for the Gauss-Newton Hessian approximation.

use casadi as ca;
use casadi::dict;

// Uniform random numbers in [0, 1) (stands in for numpy.random.random)
fn random(n: usize, seed: u64) -> Vec<f64> {
    let mut s = seed;
    (0..n)
        .map(|_| {
            s = s.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            (s >> 11) as f64 / (1u64 << 53) as f64
        })
        .collect()
}

fn main() {
    // ############ SETTINGS #####################
    let n = 10000; // Number of samples
    let fs = 610.1; // Sampling frequency [hz]

    let param_truth = ca::DM::new(vec![5.625e-6, 2.3e-4, 1.0, 4.69]);
    let param_guess = ca::DM::new(vec![5.0, 2.0, 1.0, 5.0]);
    let scale = ca::DM::new(vec![1e-6, 1e-4, 1.0, 1.0]);

    // ############ MODELING #####################
    let y = ca::MX::sym("y");
    let dy = ca::MX::sym("dy");
    let u = ca::MX::sym("u");

    let states = ca::vertcat(vec![&y, &dy]);
    let controls = u.clone();

    let m = ca::MX::sym("M");
    let c = ca::MX::sym("c");
    let k = ca::MX::sym("k");
    let k_nl = ca::MX::sym("k_NL");

    let params = ca::vertcat(vec![&m, &c, &k, &k_nl]);

    let rhs = ca::vertcat(vec![dy.clone(), (&u - &k_nl * ca::power(&y, 3.0) - &k * &y - &c * &dy) / &m]);

    // Form an ode function
    let ode = ca::Function::new("ode", vec![&states, &controls, &params], vec![&rhs]);

    // ############ Creating a simulator ##########
    let n_steps_per_sample = 10;
    let dt = 1.0 / fs / n_steps_per_sample as f64;

    // Build an integrator for this system: Runge Kutta 4 integrator
    let k1: ca::MX = ode.eval((&states, &controls, &params));
    let k2: ca::MX = ode.eval((&states + dt / 2.0 * &k1, &controls, &params));
    let k3: ca::MX = ode.eval((&states + dt / 2.0 * &k2, &controls, &params));
    let k4: ca::MX = ode.eval((&states + dt * &k3, &controls, &params));

    let states_final = &states + dt / 6.0 * (&k1 + 2.0 * &k2 + 2.0 * &k3 + &k4);

    // Create a function that simulates one step propagation in a sample
    let one_step = ca::Function::new("one_step", vec![&states, &controls, &params], vec![&states_final]);

    let mut x = states.clone();

    for _i in 0..n_steps_per_sample {
        x = one_step.eval((&x, &controls, &params));
    }

    // Create a function that simulates all step propagation on a sample
    let one_sample = ca::Function::new("one_sample", vec![&states, &controls, &params], vec![&x]);

    // ############ Simulating the system ##########
    let all_samples = one_sample.mapaccum_named("all_samples", n);

    // Choose an excitation signal
    let u_data = 0.1 * ca::DM::new(random(n as usize, 0));

    let x0 = ca::DM::new(vec![0.0, 0.0]);
    let x_measured: ca::DM = all_samples.eval((&x0, &u_data, ca::repmat_with(&param_truth, 1, n)));

    let y_data = x_measured.at((0, ..)).T();

    // You may add some noise here
    //y_data+= 0.001*numpy.random.random(N)
    // When noise is absent, the fit will be perfect.

    // Use just-in-time compilation to speed up the evaluation
    let (with_jit, compiler) = if ca::Importer::has_plugin("clang") {
        (true, "clang")
    } else if ca::Importer::has_plugin("shell") {
        println!("WARNING: on Windows, JIT may require a special environment cfr https://github.com/casadi/casadi/wiki/FAQ:-how-to-perform-jit-for-function-evaluations-of-my-optimization-problem%3F");
        (true, "shell")
    } else {
        println!("WARNING; running without jit. This may result in very slow evaluation times");
        (false, "")
    };

    // ############ Create a Gauss-Newton solver ##########
    let gauss_newton = |e: &ca::MX, nlp: &ca::rt::Dict, v: &ca::MX| -> ca::Function {
        let j = ca::jacobian(e, v);
        let h = ca::triu(&ca::mtimes(&j.T(), &j));
        let sigma = ca::MX::sym("sigma");
        let hess_lag = ca::Function::from_dict_with(
            "nlp_hess_l",
            dict! {"x" => v, "lam_f" => &sigma, "hess_gamma_x_x" => &sigma * &h},
            vec!["x", "p", "lam_f", "lam_g"],
            vec!["hess_gamma_x_x"],
            &dict! {"jit" => with_jit, "compiler" => compiler},
        );
        let opts = dict! {"hess_lag" => &hess_lag, "jit" => with_jit, "compiler" => compiler};
        ca::nlpsol_with("solver", "ipopt", nlp, &opts)
    };

    // ############ Identifying the simulated system: single shooting strategy ##########

    // Note, it is in general a good idea to scale your decision variables such
    // that they are in the order of ~0.1..100
    let x_symbolic: ca::MX = all_samples.eval((&x0, &u_data, ca::repmat_with(&params * &scale, 1, n)));

    let e = &y_data - x_symbolic.at((0, ..)).T();
    let nlp = dict! {"x" => &params, "f" => 0.5 * ca::dot(&e, &e)};

    let solver = gauss_newton(&e, &nlp, &params);

    let sol = solver.call(&dict! {"x0" => &param_guess}).unwrap();

    println!("{}", &sol["x"] * &scale);

    assert!(ca::norm_inf(&(&sol["x"] * &scale - &param_truth)).to_double() < 1e-8);

    // ############ Identifying the simulated system: multiple shooting strategy ##########

    // All states become decision variables
    let xs = ca::MX::sym_with("X", 2, n);

    let xn: ca::MX = one_sample.map_with(n, "openmp").eval((&xs, u_data.T(), ca::repmat_with(&params * &scale, 1, n)));

    let gaps = xn.at((.., ..-1)) - xs.at((.., 1..));

    let e = &y_data - xn.at((0, ..)).T();

    let v = ca::veccat(vec![&params, &xs]);

    let nlp = dict! {"x" => &v, "f" => 0.5 * ca::dot(&e, &e), "g" => ca::vec(&gaps)};

    // Multipleshooting allows for careful initialization
    let yd = (y_data.at(1..) - y_data.at(..-1)) * fs;
    let x_guess = ca::horzcat(vec![y_data.clone(), ca::vertcat(vec![yd.clone(), yd.at(-1)])]).T();

    let x0 = ca::veccat(vec![&param_guess, &x_guess]);

    let solver = gauss_newton(&e, &nlp, &v);

    let sol = solver.call(&dict! {"x0" => &x0, "lbg" => 0.0, "ubg" => 0.0}).unwrap();

    println!("{}", sol["x"].at(..4) * &scale);

    assert!(ca::norm_inf(&(sol["x"].at(..4) * &scale - &param_truth)).to_double() < 1e-8);
}

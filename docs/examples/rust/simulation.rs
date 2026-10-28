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
// casadi.tools structures (struct_symSX, struct_SX) are not available in Rust:
// the states, controls and parameters are plain symbolic vectors of named scalars.
// The remainder of the Python original (ControlSimulator, covariance propagation)
// is commented out there and is not ported.

use casadi as ca;

fn main() {
    // System states
    let [x, y, dx, dy] = ["x", "y", "dx", "dy"].map(ca::SX::sym);
    let states = ca::vertcat(vec![&x, &y, &dx, &dy]);
    let state_keys = ["x", "y", "dx", "dy"];

    // System controls
    let [u, v] = ["u", "v"].map(ca::SX::sym);
    let controls = ca::vertcat(vec![&u, &v]);

    // System parameters
    let [k, c, beta] = ["k", "c", "beta"].map(ca::SX::sym);
    let parameters = ca::vertcat(vec![&k, &c, &beta]);

    // Provide some numerical values
    let mut parameters_ = ca::DM::zeros_with(3, 1);
    parameters_.set_at(0, 10.0); // k
    parameters_.set_at(2, 1.0); // beta
    parameters_.set_at(1, 1.0); // c

    let vel = ca::vertcat(vec![&dx, &dy]);
    let p = ca::vertcat(vec![&x, &y]);
    let q = ca::vertcat(vec![&u, &v]);

    // System dynamics
    let f: ca::SX = -&k * (&p - &q) - &beta * &v * ca::sqrt(&(ca::sumsqr(&vel) + ca::power(&c, 2.0)));

    // System right hand side
    let rhs = ca::vertcat(vec![dx.clone(), dy.clone(), f.at(0), f.at(1)]);

    println!("states: {:?} = {}", state_keys, states);
    println!("controls: {}", controls);
    println!("parameters: {}", parameters);
    println!("parameters_: {}", parameters_);
    println!("rhs: {}", rhs);
}

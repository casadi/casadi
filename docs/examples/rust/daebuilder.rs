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

// Example on how to use the DaeBuilder class
// Joel Andersson, UW Madison 2017

use casadi as ca;
use casadi::dict;

fn main() {
    // Start with an empty DaeBuilder instance
    let dae = ca::DaeBuilder::new("rocket");

    // Add input expressions
    let a = dae.add_var("a", "parameter", "tunable");
    let b = dae.add_var("b", "parameter", "tunable");
    let u = dae.add_3("u", "input", dict! {});
    let h = dae.add("h");
    let v = dae.add("v");
    let m = dae.add("m");

    // Constants
    let g = 9.81; // gravity

    // Set ODE right-hand-side
    dae.eq(&dae.der_with(&h), &v);
    dae.eq(&dae.der_with(&v), &((&u - &a * ca::power(&v, 2.0)) / &m - g));
    dae.eq(&dae.der_with(&m), &(-&b * ca::power(&u, 2.0)));

    // Specify initial conditions
    dae.set_start("h", 0.0);
    dae.set_start("v", 0.0);
    dae.set_start("m", 1.0);

    // Add meta information
    dae.set_unit("h", "m");
    dae.set_unit("v", "m/s");
    dae.set_unit("m", "kg");

    // Print DAE
    dae.disp_with(true);
}

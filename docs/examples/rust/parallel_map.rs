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
use std::f64::consts::PI;
use std::time::Instant;

fn main() {
    println!("importing casadi...");

    // number of inputs to evaluate in parallel
    let n = 300;

    // dummy input
    let dummy_input: Vec<f64> = (0..n).map(|k| 2.0 * PI * k as f64 / (n - 1) as f64).collect();

    // make a dummy function that's moderately expensive to evaluate
    println!("creating dummy function....");
    let x = ca::SX::sym("x");
    let mut y = x.clone();
    for _ in 0..100000 {
        y = ca::sin(&y);
    }
    let f0 = ca::Function::new("f", vec![&x], vec![&y]);

    // evaluate it serially, the old-fasioned way
    let xs = ca::MX::sym_with("x", n, None);
    let ys: Vec<ca::MX> = (0..n).map(|k| f0.eval([xs.at(k)])).collect();
    let f_naive_parallel = ca::Function::new("fParallel", vec![&xs], vec![ca::vertcat(&ys)]);

    println!("evaluating naive parallel function...");
    let t0 = Instant::now();
    let _out_naive: ca::DM = f_naive_parallel.eval([&dummy_input]);
    println!("evaluated naive parallel function in {:.3} seconds", t0.elapsed().as_secs_f64());

    // evaluate it using new serial map construct
    let f_map = f0.map(n);

    println!("evaluating serial map function...");
    let t0 = Instant::now();
    let _out_map: ca::DM = f_map.eval([&dummy_input]);
    println!("evaluated serial map function in {:.3} seconds", t0.elapsed().as_secs_f64());
    // the following has different shaped outputs, so it's commented out
    //println!("{}", out_naive == out_map);

    // evaluate it using new parallel map construct
    let f_map = f0.map_with(n, "openmp");

    println!("evaluating parallel map function...");
    let t0 = Instant::now();
    let _out_map: ca::DM = f_map.eval([&dummy_input]);
    println!("evaluated parallel map function in {:.3} seconds", t0.elapsed().as_secs_f64());
    // the following has different shaped outputs, so it's commented out
    //println!("{}", out_naive == out_map);
}

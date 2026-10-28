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

use casadi as ca;
use std::process::Command;
use std::time::Instant;

fn system(cmd: &str) {
    Command::new("sh").arg("-c").arg(cmd).status().expect("failed to run the compiler");
}

fn main() {
    let mut compileme = true; // Flag wether compiling should take place or not

    if std::env::args().nth(1).as_deref() == Some("nc") {
        compileme = false;
    } else {
        println!("Info: Use 'cargo run --example c_code_generation nc' to omit compiling");
    }

    // Form an expression for the gradient of the determinant
    let x = ca::SX::sym_with("x", 7, 7);
    let gd = ca::gradient(&ca::det(&x), &x);

    // Random point to evaluate it
    let x0 = ca::DM::rand_with(7, 7);

    // Form a function and generate C code
    let name = "grad_det";
    let grad_det = ca::Function::named(name, vec![&x], vec![&gd], vec!["x"], vec!["gd"]);
    let cname = grad_det.generate();

    let oname_no_opt = format!("{}_no_opt.so", name);
    println!("Compiling without optimization:  {}", oname_no_opt);
    let t1 = Instant::now();
    if compileme {
        system(&format!("gcc -fPIC -shared {} -o {}", cname, oname_no_opt));
    }
    println!("time =  {}  ms", t1.elapsed().as_secs_f64() * 1e3);

    let oname_o3 = format!("{}_O3.so", name);
    println!("Compiling with O3 optimization:  {}", oname_o3);
    let t1 = Instant::now();
    if compileme {
        system(&format!("gcc -fPIC -shared -O3 {} -o {}", cname, oname_o3));
    }
    println!("time =  {}  ms", t1.elapsed().as_secs_f64() * 1e3);

    let oname_os = format!("{}_Os.so", name);
    println!("Compiling with Os optimization:  {}", oname_os);
    let t1 = Instant::now();
    if compileme {
        system(&format!("gcc -fPIC -shared -Os {} -o {}", cname, oname_os));
    }
    println!("time =  {}  ms", t1.elapsed().as_secs_f64() * 1e3);

    // Read function
    let grad_det_no_opt = ca::external_with(name, format!("./{}", oname_no_opt), None);
    let grad_det_o3 = ca::external_with(name, format!("./{}", oname_o3), None);
    let grad_det_os = ca::external_with(name, format!("./{}", oname_os), None);
    let f_test = [&grad_det, &grad_det_no_opt, &grad_det_o3, &grad_det_os];

    for f in f_test {
        let t1 = Instant::now();
        let nrep = 10000;
        let mut r = ca::DM::default();
        for _ in 0..nrep {
            r = f.eval([&x0]);
        }
        let dt = t1.elapsed().as_secs_f64() / nrep as f64;
        println!("result =  {:?}", r.nonzeros());
        println!("time =  {}  ms", dt * 1e3);

        let num_op = grad_det.n_nodes();
        println!("number of elementary operations:  {}", num_op);
        println!("time per elementary operations:  {}  ns", dt / num_op as f64 * 1e9);
    }
}

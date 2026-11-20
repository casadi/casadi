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
use std::process::Command;

fn main() {
    // Test problem
    //
    //    min x^2 + y^2
    //    s.t.    x + y - 10 = 0
    //

    // Optimization variables
    let x = ca::MX::sym("x");
    let y = ca::MX::sym("y");

    // Objective
    let f = &x * &x + &y * &y;

    // Constraints
    let g = &x + &y - 10.0;

    // Create an NLP problem structure
    let nlp = dict! {"x" => ca::vertcat(vec![&x, &y]), "f" => &f, "g" => &g};

    // Pick a compiler
    let compiler = "gcc"; // Linux
    // let compiler = "clang";  // OSX
    // let compiler = "cl.exe"; // Windows

    // Run this program in an environment that recognised the compiler as command.
    // On Windows, the suggested way is to run this program from a "x64 Native Tools Command Promt for VS" (Requires Visual C++ components or Build Tools for Visual Studio, available from Visual Studio installer. You also need SDK libraries in order to access stdio and math.)

    let flags = vec!["-O3"]; // Linux/OSX

    for mode in ["jit", "external"] {
        let solver = if mode == "jit" {
            // By default, the compiler will be gcc or cl.exe
            let jit_options = dict! {"flags" => &flags, "verbose" => true, "compiler" => compiler};
            let options = dict! {"jit" => true, "compiler" => "shell", "jit_options" => jit_options};

            // Create an NLP solver instance
            ca::nlpsol_with("solver", "ipopt", &nlp, &options)
        } else {
            // Create an NLP solver instance
            let solver = ca::nlpsol("solver", "ipopt", &nlp);

            // Generate C code for the NLP functions
            solver.generate_dependencies("nlp.c");

            // On Windows, use other flags
            Command::new(compiler)
                .args(["-fPIC", "-shared"])
                .args(&flags)
                .args(["nlp.c", "-o", "nlp.so"])
                .status()
                .expect("failed to run the compiler");

            // Create a new NLP solver instance from the compiled code
            ca::nlpsol("solver", "ipopt", "./nlp.so")
        };

        let arg = dict! {
            "lbx" => -ca::DM::inf(),
            "ubx" => ca::DM::inf(),
            "lbg" => 0,
            "ubg" => 0,
            "x0" => 0,
        };

        // Solve the NLP
        let res = solver.call(&arg).unwrap();

        // Print solution
        println!("-----");
        println!("objective at solution = {}", res["f"]);
        println!("primal solution = {}", res["x"]);
        println!("dual solution (x) = {}", res["lam_x"]);
        println!("dual solution (g) = {}", res["lam_g"]);
    }
}

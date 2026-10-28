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
    // Declare variables
    let x = ca::SX::sym_with("x", 2, None);

    // Form the NLP
    let f = ca::power(&x.at(0), 2.0) + ca::power(&x.at(1), 2.0); // objective
    let g = x.at(0) + x.at(1) - 10.0; // constraint
    let nlp = dict! {"x" => &x, "f" => &f, "g" => &g};

    // Pick an NLP solver
    let my_solver = "ipopt";
    //let my_solver = "worhp";
    //let my_solver = "sqpmethod";

    // Solver options
    let opts = if my_solver == "sqpmethod" {
        dict! {"qpsol" => "qpoases", "qpsol_options" => dict! {"printLevel" => "none"}}
    } else {
        dict! {}
    };

    // Allocate a solver
    let solver = ca::nlpsol_with("solver", my_solver, &nlp, &opts);

    // Solve the NLP
    let sol = solver.call(&dict! {"lbg" => 0}).unwrap();

    // Print solution
    println!("-----");
    println!("objective at solution =  {}", sol["f"]);
    println!("primal solution =  {}", sol["x"]);
    println!("dual solution (x) =  {}", sol["lam_x"]);
    println!("dual solution (g) =  {}", sol["lam_g"]);
}

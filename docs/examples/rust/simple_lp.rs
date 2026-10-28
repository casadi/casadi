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

const INF: f64 = f64::INFINITY;

// minimize    3x + 4y
// subject to  x + 2y <= 14
//            3x -  y >= 0
//             x -  y <= 2

fn main() {
    // Sparsity of the LP linear term
    let a_sp = ca::Sparsity::dense_with(3, 2);

    // Create solver
    let solver = ca::conic("solver", "qpoases", &dict! {"a" => &a_sp});
    //let solver = ca::conic("solver", "clp", &dict! {"a" => &a_sp}); // Use clp

    let g = ca::DM::new(vec![3.0, 4.0]);
    let a = ca::DM::new(vec![vec![1.0, 2.0], vec![3.0, -1.0], vec![1.0, -1.0]]);
    let lba = ca::DM::new(vec![-INF, 0.0, -INF]);
    let uba = ca::DM::new(vec![14.0, INF, 2.0]);

    let sol = solver.call(&dict! {"g" => &g, "a" => &a, "lba" => &lba, "uba" => &uba}).unwrap();
    let mut keys: Vec<_> = sol.keys().collect();
    keys.sort();
    let items: Vec<String> = keys.iter().map(|k| format!("'{}': {}", k, sol[*k])).collect();
    println!("{{{}}}", items.join(", "));
}

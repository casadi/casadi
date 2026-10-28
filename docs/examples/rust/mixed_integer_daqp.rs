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
use casadi::dict;

// numpy's default_rng(0) is not reproduced: a small seeded generator gives standard normals instead
struct Rng(u64);

impl Rng {
    fn random(&mut self) -> f64 {
        // SplitMix64
        self.0 = self.0.wrapping_add(0x9E3779B97F4A7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58476D1CE4E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D049BB133111EB);
        z ^= z >> 31;
        (z >> 11) as f64 / (1u64 << 53) as f64
    }
    fn standard_normal(&mut self) -> f64 {
        // Box-Muller
        let u1 = 1.0 - self.random();
        let u2 = self.random();
        (-2.0 * u1.ln()).sqrt() * (2.0 * std::f64::consts::PI * u2).cos()
    }
    fn matrix(&mut self, nrow: usize, ncol: usize, mut gen: impl FnMut(&mut Self) -> f64) -> ca::DM {
        ca::DM::new((0..nrow).map(|_| (0..ncol).map(|_| gen(self)).collect()).collect::<Vec<Vec<f64>>>())
    }
}

fn allclose(a: &[f64], b: &[f64], rtol: f64, atol: f64) -> bool {
    a.len() == b.len() && a.iter().zip(b).all(|(x, y)| (x - y).abs() <= atol + rtol * y.abs())
}

fn main() {
    if std::env::var_os("SKIP_DAQP_TESTS").is_some() || !ca::has_conic("daqp") {
        std::process::exit(0);
    }

    let mut rng = Rng(0);
    let n = 100;
    let m = 80;
    let ms = 55;
    let x = ca::SX::sym_with("x", n, 1);
    // --- Hessian: dense symmetric matrix ---
    let mm = rng.matrix(n, n, Rng::standard_normal);
    let h = 0.5 * ca::mtimes(&mm, &mm.T()); // symmetrize

    // --- Linear term ---
    let f = 100.0 * rng.matrix(n, 1, Rng::standard_normal);

    // --- Constraint matrix and bounds ---
    let a = rng.matrix(m, n, Rng::standard_normal);
    let bupper = 20.0 * rng.matrix(m, 1, Rng::random);
    let blower = -20.0 * rng.matrix(m, 1, Rng::random);

    let discrete: Vec<bool> = (0..n).map(|i| i < ms).collect();
    let lbx: Vec<f64> = (0..n).map(|i| if i < ms { 0.0 } else { -10.0 }).collect();
    let ubx: Vec<f64> = (0..n).map(|i| if i < ms { 1.0 } else { 10.0 }).collect();
    let qp = dict! {
        "f" => 0.5 * ca::mtimes(&ca::mtimes(&x.T(), &h), &x) + ca::mtimes(&f.T(), &x),
        "x" => &x,
        "g" => ca::mtimes(&a, &x),
    };

    // Using casadi qpsol
    let solver = ca::qpsol_with("solver", "daqp", &qp, &dict! {"discrete" => &discrete});
    let daqp_sol = solver
        .call(&dict! {"lbx" => &lbx, "ubx" => &ubx, "lbg" => &blower, "ubg" => &bupper})
        .unwrap();
    println!("Optimal solution: {:?}", daqp_sol["x"].nonzeros());
    println!("Optimal objective: {}", daqp_sol["f"].to_double());

    // Using casadi conic
    // let solver = ca::conic_with("solver", "daqp",
    //                   &dict! {"h" => h.sparsity(), "a" => a.sparsity()},
    //                   &dict! {"discrete" => &discrete});
    // let daqp_sol = solver.call(&dict! {
    //         "h" => &h,
    //         "g" => &f,
    //         "a" => &a,
    //         "lbx" => &lbx, "ubx" => &ubx,
    //         "lba" => &blower, "uba" => &bupper}).unwrap();
    // println!("Optimal solution: {:?}", daqp_sol["x"].nonzeros());
    // println!("Optimal objective: {}", daqp_sol["cost"].to_double());

    let stats = solver.stats();
    println!(
        "DAQP BnB infos: solver.stats()['bnb_itercount']={:?}, solver.stats()['bnb_nodecount']={:?}",
        stats["bnb_itercount"], stats["bnb_nodecount"]
    );

    // The remainder cross-checks DAQP against gurobi; skip it alone so that the DAQP
    // part above still runs when gurobi is unavailable.
    if std::env::var_os("SKIP_GUROBI_TESTS").is_some() || !ca::has_conic("gurobi") {
        std::process::exit(0);
    }

    let solver = ca::qpsol_with("solver", "gurobi", &qp, &dict! {"discrete" => &discrete});
    let gu_sol = solver
        .call(&dict! {"lbx" => &lbx, "ubx" => &ubx, "lbg" => &blower, "ubg" => &bupper})
        .unwrap();
    println!("Optimal solution: {:?}", gu_sol["x"].nonzeros());
    println!("Optimal objective: {}", gu_sol["f"].to_double());

    if !allclose(&daqp_sol["x"].nonzeros(), &gu_sol["x"].nonzeros(), 1e-6, 1e-6) {
        panic!("solutions are different");
    }
    if !allclose(&[gu_sol["f"].to_double()], &[daqp_sol["f"].to_double()], 1e-6, 1e-6) {
        panic!("optimal objectives are different");
    }
}

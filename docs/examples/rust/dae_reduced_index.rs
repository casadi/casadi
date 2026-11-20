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
use std::collections::{BTreeMap, HashMap};

// Largest absolute entry of each row of a matrix
fn max_abs_rows(m: &ca::DM) -> Vec<f64> {
    (0..m.size1()).map(|i| m.at((i, ..)).nonzeros().iter().fold(0.0, |a: f64, b| a.max(b.abs()))).collect()
}

fn main() {
    // Pendulum: point mass on massless rod
    // x'' = T*x, y'' = T*y-g, x**2+y**2-L**2
    //
    //
    //     .
    //  L /
    //   /
    //  /
    // o (x,y)
    //
    let x = ca::SX::sym("x"); // horizontal position of point
    let y = ca::SX::sym("y"); // vertical position of point

    let dx = ca::SX::sym("dx");
    let dy = ca::SX::sym("dy");

    let u = ca::SX::sym("u"); // helper states to bring second order system to first order
    let v = ca::SX::sym("v");

    let du = ca::SX::sym("du");
    let dv = ca::SX::sym("dv");

    // Force
    let t = ca::SX::sym("T");

    let l = 1.0;
    let g = 9.81;

    let alg = ca::vertcat(vec![&dx - &u, &dy - &v, &du - &t * &x, &dv - &t * &y + g,
                               ca::power(&x, 2.0) + ca::power(&y, 2.0) - l * l]);

    // Implicit differential states
    let x_impl = ca::vertcat(vec![&x, &y, &u, &v]);
    // Derivatives of implict differential states
    let dx_impl = ca::vertcat(vec![&dx, &dy, &du, &dv]);

    let z = t.clone();

    let dae = HashMap::from([("x_impl", x_impl.clone()), ("dx_impl", dx_impl), ("z", z), ("alg", alg)]);

    // Perform index reduction
    let (dae_reduced, stats) = ca::dae_reduce_index(&dae);
    println!("We had an index 3 system (reduced to index 1 now)");
    println!("{:?}", stats);
    println!("Here's the reduced DAE:");
    println!("{:?}", dae_reduced.iter().collect::<BTreeMap<_, _>>());
    println!("Notice how the algebraic equation now contains a second derivative of x**2+y**2-L**2");
    println!("Also note how x**2+y**2-L**2 and its first derivate are absent from the algebraic equations.");
    println!("Both are kept in 'I' though: the resulting DAE invariants.");

    // The DAE is not yet in a form that CasADi integrator can handle (semi-explicit).
    // Let's convert it, ad obtain some adaptors
    let (dae_se, state_to_orig, phi) = ca::dae_map_semi_expl(&dae, &dae_reduced);

    let grid: ca::DM = ca::linspace(0.0, 5.0, 1000);
    let grid = grid.nonzeros();
    let intg = ca::integrator_with("intg", "idas", &dae_se, 0.0, &grid, None);

    // Before we can integrate, we need a consistant initial guess
    // Let's say we start at y=-0.5, dx=-0.1;
    // consistency equations should automatically identify everything else

    // Encode the desire to only enforce y and dx
    let free = 0.0; // default
    let force = -1.0;
    let init_strength = HashMap::from([
        ("x_impl", ca::DM::new(vec![free, force, free, free])),
        ("dx_impl", ca::DM::new(vec![force, free, free, free])),
        ("z", ca::DM::new(vec![free])),
    ]);

    // Obtain a generating Function for initial guesses
    let init_gen = ca::dae_init_gen_with(&dae, &dae_reduced, "ipopt", &init_strength, None);

    let blank = 0.0; // Will not be enforced but still used as initial guess
    // suggest to initialize with the left-hanging solution by having x=-1 as initial guess
    let init = HashMap::from([
        ("x_impl", ca::DM::new(vec![-1.0, -0.5, blank, blank])),
        ("dx_impl", ca::DM::new(vec![-0.1, blank, blank, blank])),
        ("z", ca::DM::new(blank)),
    ]);

    let xz0 = init_gen.call(&init).unwrap();
    println!("We have a consistent initial guess now to feed to the integrator:");
    println!("{:?}", xz0.iter().collect::<BTreeMap<_, _>>());

    println!("Look, we found that force in the pendulum at t=0 equals:");
    println!("{}", state_to_orig.call(&dict! {"xf" => &xz0["x0"], "zf" => &xz0["z0"]}).unwrap()["z"]);

    println!("A consistent initial guess should make the invariants zero (up to solver precision):");
    let phi0 = phi.call(&dict! {"x" => &xz0["x0"], "z" => &xz0["z0"]}).unwrap();
    println!("{:?}", phi0.iter().collect::<BTreeMap<_, _>>());

    // Integrate and get resultant xf,zf on grid
    let sol = intg.call(&xz0).unwrap();

    println!("({}, {})", sol["xf"].size1(), sol["xf"].size2());

    // Solution projected into original DAE space
    let sol_orig = state_to_orig.call(&dict! {"xf" => &sol["xf"], "zf" => &sol["zf"]}).unwrap();

    // We can see the large-angle pendulum motion play out well in the u state
    let names: Vec<String> = ca::vertsplit(&x_impl).iter().map(|e| e.name()).collect();
    println!("{:?} at t = 5: {}", names, sol_orig["x_impl"].at((.., -1)));

    // A perfect integrator will perfectly preserve the values of the invariants over time
    // Integrator errors make the invariants drift in practice
    // This is not a detail; the pendulum length is in fact growing!
    let error = &phi.call(&dict! {"x" => &sol["xf"], "z" => &sol["zf"]}).unwrap()["I"];
    println!("Invariants at t = 5: {}", error.at((.., -1)));

    // There are some techniques to avoid the drifting of invariants associated with index reduction
    //  - Method of Dummy Derviatives (not implemented)
    //    hard implementation details: may need to switch between different choices online
    //     -> integration in parts triggered by zero-crossing events
    //  - Stop the integration every once in a while and project back (not implemented)
    //  - Baumgarte stabilization (implemented): build into the equations a form of continuous feedback that brings back deviations in invariants back into the origin

    // Demonstrate Baumgarte stabilization with a pole of -1.
    // Drift is absent now.
    let (dae_reduced, _stats) = ca::dae_reduce_index_with(&dae, dict! {"baumgarte_pole" => -1});
    let (dae_se, _state_to_orig, phi) = ca::dae_map_semi_expl(&dae, &dae_reduced);
    let intg = ca::integrator_with("intg", "idas", &dae_se, 0.0, &grid, None);
    let init_gen = ca::dae_init_gen_with(&dae, &dae_reduced, "ipopt", &init_strength, None);
    let xz0 = init_gen.call(&init).unwrap();
    let sol = intg.call(&xz0).unwrap();
    let error = &phi.call(&dict! {"x" => &sol["xf"], "z" => &sol["zf"]}).unwrap()["I"];
    println!("Invariants at t = 5 with Baumgarte pole=-1: {}", error.at((.., -1)));

    // Sweep across different integrator precisions and Baumgarte poles.
    // The choice of pole is not really straightforward
    let poles = [-0.1, -1.0, -10.0];
    let precisions = [0.1, 1.0, 10.0];

    let abstol_default = 1e-8;
    let reltol_default = 1e-6;

    for &pole in &poles {
        for &precision in &precisions {
            let (dae_reduced, _stats) = ca::dae_reduce_index_with(&dae, dict! {"baumgarte_pole" => pole});
            let (dae_se, _state_to_orig, phi) = ca::dae_map_semi_expl(&dae, &dae_reduced);
            let intg = ca::integrator_with("intg", "idas", &dae_se, 0.0, &grid,
                dict! {"abstol" => abstol_default * precision, "reltol" => reltol_default * precision});
            let init_gen = ca::dae_init_gen_with(&dae, &dae_reduced, "ipopt", &init_strength,
                dict! {"ipopt.print_level" => 0, "print_time" => false});
            let xz0 = init_gen.call(&init).unwrap();
            let sol = intg.call(&xz0).unwrap();
            let ca::rt::Value::Int(nfevals) = intg.stats()["nfevals"] else { panic!("nfevals") };
            let error = &phi.call(&dict! {"x" => &sol["xf"], "z" => &sol["zf"]}).unwrap()["I"];
            println!("pole={:.1}, abstol={:.0e} => nfevals={}, max |I| = {:?}",
                     pole, abstol_default * precision, nfevals, max_abs_rows(error));
        }
    }
}

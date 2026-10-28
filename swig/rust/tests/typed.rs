use casadi::*;
use std::collections::HashMap;

#[test]
fn typed_rosenbrock() {
    let x = MX::sym("x");
    let y = MX::sym("y");
    let d = minus(&y, &times(&x, &x));
    let f = plus(&times(&minus(1.0, &x), &minus(1.0, &x)), &times(100.0, &times(&d, &d)));
    let xy: MX = vertcat(vec![&x, &y]);
    let solver: Function = nlpsol_with(
        "S",
        "ipopt",
        &dict! {"x" => &xy, "f" => &f},
        &dict! {"print_time" => false, "ipopt.print_level" => 0},
    );
    let res: HashMap<String, DM> = solver.call(&dict! {"x0" => vec![-1.0, 1.0]}).unwrap();
    let xopt: Vec<f64> = res["x"].nonzeros();
    assert!((xopt[0] - 1.0).abs() < 1e-6 && (xopt[1] - 1.0).abs() < 1e-6, "{:?}", xopt);
}

#[test]
fn promotion_and_results() {
    let x = SX::sym("x");
    let a = DM::new(vec![1.0, 2.0]);
    // DM + SX -> SX, statically
    let s: SX = plus(&a, &x);
    assert_eq!(s.str(), "[(1+x), (2+x)]");
    // DM + DM -> DM
    let b: DM = times(&a, &a);
    assert_eq!(b.nonzeros(), vec![1.0, 4.0]);
    let f = Function::named("f", vec![&x], vec![&sin(&x)], vec!["x"], vec!["y"]);
    let r: Vec<DM> = f.call(vec![DM::new(0.5)]).unwrap();
    assert!((r[0].nonzeros()[0] - 0.5f64.sin()).abs() < 1e-12);
    // wrong number of inputs: an error, not a panic
    assert!(f.call(vec![DM::new(0.5), DM::new(1.0)]).is_err());
    let e: Vec<SX> = f.call_symbolic(vec![&x]);
    assert_eq!(e[0].str(), "sin(x)");
}

#[test]
fn naming_rules() {
    let x = MX::sym("x");
    let sp = Sparsity::dense_with(2, 3);
    let y = MX::sym_sparsity("y", &sp);
    assert_eq!((y.size1(), y.size2()), (2, 3));
    let z = MX::sym_with("z", 3, 1);
    assert_eq!(z.numel(), 3);
    let v: Vec<MX> = MX::sym_vec("v", &sp, 4);
    assert_eq!(v.len(), 4);
    assert_eq!(DM::zeros_sparsity(&sp).nnz(), 6);
    let f = Function::new("f", vec![&x], vec![&x]);
    let g = f.clone();
    assert_eq!(g.name(), "f");
    assert!(Function::default().is_null());
    let opti = Opti::default();
    let w = opti.variable();
    opti.subject_to(&ge(&w, 0.0));
    assert_eq!(opti.ng(), 1);
    opti.clear_constraints();
    assert_eq!(opti.ng(), 0);
}

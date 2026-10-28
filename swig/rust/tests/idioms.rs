use casadi::*;

#[test]
fn operators() {
    let x = SX::sym("x");
    let a = DM::new(vec![1.0, 2.0]);
    let e = -(&a + &x) * 2.0 - &x / 3.0;
    assert_eq!(e.str(), "@1=(x/3), [((2.*(-(1+x)))-@1), ((2.*(-(2+x)))-@1)]");
    let mut b = &a * &a;
    b += 1.0;
    assert_eq!(b.nonzeros(), vec![2.0, 5.0]);
    let m = MX::sym("m");
    let y: MX = 1.0 + &m * &m;
    assert_eq!(format!("{}", y), "(1+sq(m))");
    let p = mtimes(&DM::new(vec![vec![1.0, 2.0]]), &DM::new(vec![3.0, 4.0]));
    assert_eq!(p.nonzeros(), vec![11.0]);
}

#[test]
fn eval_positional() {
    let x = MX::sym("x");
    let y = MX::sym("y");
    let f = Function::new("f", vec![&x, &y], vec![&x + &y, &x * &y]);
    let (s, p): (DM, DM) = f.eval((2.0, 3.0));
    assert_eq!((s.nonzeros(), p.nonzeros()), (vec![5.0], vec![6.0]));
    let (s, _): (MX, MX) = f.eval([&x, &x]);
    assert_eq!(s.str(), "f(x, x){0}");
    let g = Function::new("g", vec![&x], vec![sin(&x)]);
    let r: DM = g.eval([0.0]);
    assert_eq!(r.nonzeros(), vec![0.0]);
    let z = SX::sym("z");
    let h = Function::new("h", vec![&z], vec![&z * 2.0]);
    let r: SX = h.eval((DM::new(1.0) + &z,));
    assert_eq!(r.str(), "(2.*(1+z))");
    let all: Vec<DM> = f.eval([1.0, 1.0]);
    assert_eq!(all.len(), 2);
}

#[test]
#[should_panic(expected = "Function has 2 outputs")]
fn eval_wrong_output_count() {
    let x = MX::sym("x");
    let f = Function::new("f", vec![&x], vec![&x, &x]);
    let _: DM = f.eval([1.0]);
}

#[test]
fn function_is_send_sync() {
    fn check<T: Send + Sync>() {}
    check::<Function>();
    let x = SX::sym("x");
    let f = Function::new("f", vec![&x], vec![&x * &x]);
    let t = std::thread::spawn(move || {
        let r: DM = f.eval([3.0]);
        r.nonzeros()[0]
    });
    assert_eq!(t.join().unwrap(), 9.0);
}

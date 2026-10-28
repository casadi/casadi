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

#[test]
fn indexing() {
    let mut a = DM::new(vec![vec![1.0, 2.0, 3.0], vec![4.0, 5.0, 6.0]]);
    assert_eq!(a.at((1, 2)).nonzeros(), vec![6.0]);
    assert_eq!(a.at((.., 0)).nonzeros(), vec![1.0, 4.0]);
    assert_eq!(a.at((0, 1..)).nonzeros(), vec![2.0, 3.0]);
    assert_eq!(a.at(-1).nonzeros(), vec![6.0]);
    assert_eq!(a.at(vec![0, 3]).nonzeros(), vec![1.0, 5.0]);
    assert_eq!(a.nz(1..=2).nonzeros(), vec![4.0, 2.0]);
    a.set_at((0, ..), 0.0);
    assert_eq!(a.nonzeros(), vec![0.0, 4.0, 0.0, 5.0, 0.0, 6.0]);
    let mut x = MX::sym_with("x", 3, 1);
    assert_eq!(x.at(..2).str(), "x[:2]");
    x.set_nz_at(0, 1.0);
    let s = SX::sym_with("s", 4, 1);
    assert_eq!(s.at(1).str(), "s_1");
}

#[test]
fn gap_fixes() {
    // enum values are constants
    assert!(OP_CONST > 0 && OP_ADD != OP_SUB);
    let x = MX::sym("x");
    let p = MX::sym("p");
    let f = Function::named("f", vec![&x, &p], vec![&x * &p], vec!["x", "p"], vec!["y"]);
    // keyword call, symbolic: numbers are promoted to MX
    let r: std::collections::HashMap<String, MX> = f.eval(dict! {"x" => &x, "p" => 2.0});
    assert_eq!(r["y"].str(), "f(x, 2){0}");
    // keyword call, numeric
    let r: std::collections::HashMap<String, DM> = f.eval(dict! {"x" => 3.0, "p" => 2.0});
    assert_eq!(r["y"].nonzeros(), vec![6.0]);
    // scalar result, slice arguments
    let a = vec![DM::new(3.0), DM::new(4.0)];
    let y: f64 = f.eval(&a);
    assert_eq!(y, 12.0);
    // float literals infer without f32
    let t = linspace(0.0, 1.0, 3).nonzeros();
    assert_eq!(t, vec![0.0, 0.5, 1.0]);
    assert_eq!(eq(&x, 0.0).numel(), 1);
    // strided slices, bool as a matrix, reshape(x, n, m)
    let v = DM::new(vec![0.0, 1.0, 2.0, 3.0, 4.0]);
    assert_eq!(v.at(step(0.., 2)).nonzeros(), vec![0.0, 2.0, 4.0]);
    assert_eq!(DM::new(true).nonzeros(), vec![1.0]);
    assert_eq!(reshape(&v.at(..4), 2, 2).size2(), 2);
    assert_eq!(logic_not(&DM::new(0.0)).nonzeros(), vec![1.0]);
    // Value prints like Python
    let s = Value::Dict([("a".to_string(), Value::Str("b".into())), ("n".to_string(), Value::Int(3))].into());
    assert_eq!(s.to_string(), r#"{"a": "b", "n": 3}"#);
    assert_eq!(Value::Str("ok".into()).to_string(), "ok");
    assert_eq!(Value::Int(3).as_f64(), Some(3.0));
}

#[test]
fn defaults_with_none() {
    let x = MX::sym_with("x", 3, None);
    assert_eq!((x.size1(), x.size2()), (3, 1));
    let opti = Opti::default();
    let v = opti.variable_with(2, 4, None);
    assert_eq!((v.size1(), v.size2()), (2, 4));
    let s = Sparsity::dense_with(3, None);
    assert_eq!(s.size2(), 1);
    let f = Function::new("f", vec![&x], vec![&x * 2.0]);
    let g = f.map_with(4, "serial");
    assert_eq!(g.size2_in(0), 4);
    let r = std::panic::catch_unwind(|| MX::sym_with("y", None, 2));
    assert!(r.is_err(), "a value after None must panic");
}

fn make_opts(h: &MX) -> Dict {
    // a dict built from locals can be returned: it owns its values
    let local = MX::sym("local");
    dict! {"h" => h, "l" => &local, "nested" => dict! {"a" => 1, "b" => "s"}, "v" => vec![1.0, 2.0]}
}

#[test]
fn dict_owns_values() {
    let x = MX::sym("x");
    let f = Function::named("f", vec![&x], vec![&x * 2.0], vec!["x0"], vec!["xf"]);
    // reassigning a variable that the dict reads (was E0506)
    let mut xk = MX::sym("xk");
    for _ in 0..2 {
        let r: std::collections::HashMap<String, MX> = f.eval(dict! {"x0" => &xk});
        xk = r["xf"].clone();
        xk = f.eval::<std::collections::HashMap<String, MX>>(dict! {"x0" => &xk})["xf"].clone();
    }
    assert_eq!(xk.str(), "f(f(f(f(xk){0}){0}){0}){0}");
    // a temporary stored in the dict (was E0716)
    let d = dict! {"t" => &Sparsity::dense(2)};
    assert_eq!(d.0.len(), 1);
    // a dict returned from a function (was E0597)
    let o = make_opts(&x);
    assert_eq!(o.0.len(), 4);
}

struct Square;

impl CallbackImpl for Square {
    fn eval(&self, arg: Vec<DM>) -> Vec<DM> {
        vec![&arg[0] * &arg[0]]
    }
    fn get_name_in(&self, _i: i64) -> String {
        "x".into()
    }
}

#[test]
fn callbacks() {
    // Callback: overridden methods run in Rust, the others in C++
    let cb = Callback::from_impl(Square);
    cb.construct("sq");
    let y: DM = cb.eval([3.0]);
    assert_eq!(y.nonzeros(), vec![9.0]);
    assert_eq!(cb.name_in_with(0), "x");
    // a panic in Rust becomes a casadi error
    struct Bad;
    impl CallbackImpl for Bad {
        fn eval(&self, _arg: Vec<DM>) -> Vec<DM> {
            panic!("bad input")
        }
    }
    let bad = Callback::from_impl(Bad);
    bad.construct("bad");
    let r = Function::from(&bad).call(vec![DM::new(1.0)]);
    assert!(r.unwrap_err().to_string().contains("bad input"));
    // Opti iteration callback capturing non-Send values
    let opti = Opti::default();
    let x = opti.variable();
    opti.minimize(&power(&(&x - 3.0), 2.0));
    opti.solver_with("ipopt", dict! {"print_time" => false, "ipopt.print_level" => 0, "ipopt.sb" => "yes"}, None);
    let seen = std::rc::Rc::new(std::cell::RefCell::new(Vec::new()));
    let (s, o, xc) = (seen.clone(), opti.clone(), x.clone());
    opti.callback(move |i| s.borrow_mut().push((i, o.debug().value(&xc).to_double())));
    let sol = opti.solve();
    assert!((sol.value(&x).to_double() - 3.0).abs() < 1e-8);
    assert!(!seen.borrow().is_empty());
    opti.clear_callback();
}

use casadi::*;

// One model for numeric and symbolic use, as a Python function on any casadi type
fn model<M: Matrix>(w: &M) -> M {
    let r = M::vertcat(vec![w.at(0).sin(), w.at(1) * 2.0]);
    r.clone() * &r + M::from(1.0)
}

#[test]
fn generic_model() {
    let w0 = DM::new(vec![0.5, 3.0]);
    let r: DM = model(&w0);
    let x = SX::sym_with("x", 2, 1);
    let fs = Function::new("fs", vec![&x], vec![&model(&x)]);
    let y = MX::sym_with("y", 2, 1);
    let fm = Function::new("fm", vec![&y], vec![&model(&y)]);
    let rs: DM = fs.eval([&w0]);
    let rm: DM = fm.eval([&w0]);
    assert_eq!(r.nonzeros(), vec![0.5f64.sin().powi(2) + 1.0, 37.0]);
    assert_eq!(rs.nonzeros(), r.nonzeros());
    assert_eq!(rm.nonzeros(), r.nonzeros());
}

fn rk4<M: Matrix>(f: impl Fn(&M) -> M, x: &M, h: f64) -> M {
    let k1 = f(x);
    let k2 = f(&(x.clone() + k1.clone() * (h / 2.0)));
    let k3 = f(&(x.clone() + k2.clone() * (h / 2.0)));
    let k4 = f(&(x.clone() + k3.clone() * h));
    x.clone() + (k1 + k2 * 2.0 + k3 * 2.0 + k4) * (h / 6.0)
}

#[test]
fn generic_integrator_step() {
    let ode = |x: &MX| MX::vertcat(vec![x.at(1), -x.at(0)]);
    let x = MX::sym_with("x", 2, 1);
    let f = Function::new("f", vec![&x], vec![&rk4(ode, &x, 0.1)]);
    let num: DM = rk4(|x: &DM| DM::vertcat(vec![x.at(1), -x.at(0)]), &DM::new(vec![1.0, 0.0]), 0.1);
    let sym: DM = f.eval([vec![1.0, 0.0]]);
    assert_eq!(num.nonzeros(), sym.nonzeros());
}

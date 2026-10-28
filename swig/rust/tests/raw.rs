use casadi::{raw, rt, Function, DM, MX, SX};
use rt::{take, ToRV, RV};
use std::collections::HashMap;

fn one<T: rt::FromRV>(r: Result<Vec<RV>, rt::Error>) -> T {
    take(r.unwrap().remove(0))
}

#[test]
fn symbolic_and_numeric() {
    let x: SX = one(raw::GenSX_sym(vec!["x".to_rv(), 2i64.to_rv()]));
    let y: SX = one(raw::sin(vec![x.to_rv()]));
    let z: SX = one(raw::plus(vec![y.to_rv(), 1.5f64.to_rv()]));
    let s: String = one(raw::SX_str(vec![z.to_rv()]));
    assert_eq!(s, "@1=1.5, [(sin(x_0)+@1), (sin(x_1)+@1)]");
    let f: Function = one(raw::new_Function(vec!["f".to_rv(), RV::list(vec![x.to_rv()]), RV::list(vec![z.to_rv()])]));
    let r: Vec<DM> = one(raw::Function_call(vec![f.to_rv(), RV::list(vec![vec![0.0f64, 1.0].to_rv()])]));
    let nz: Vec<f64> = one(raw::DM_nonzeros(vec![r[0].to_rv()]));
    assert!((nz[0] - 1.5).abs() < 1e-12 && (nz[1] - 1.5 - 1f64.sin()).abs() < 1e-12, "{:?}", nz);
}

#[test]
fn matrix_literals() {
    let a: DM = one(raw::new_DM(vec![vec![vec![1.0f64, 2.0], vec![3.0, 4.0]].to_rv()]));
    assert_eq!(one::<String>(raw::DM_str(vec![a.to_rv()])), "\n[[1, 2], \n [3, 4]]");
    let f: Vec<f64> = one(raw::DM_nonzeros(vec![a.to_rv()]));
    assert_eq!(f, vec![1.0, 3.0, 2.0, 4.0]);
}

#[test]
fn rosenbrock_ipopt() {
    let x: MX = one(raw::GenMX_sym(vec!["x".to_rv()]));
    let y: MX = one(raw::GenMX_sym(vec!["y".to_rv()]));
    let x2: MX = one(raw::times(vec![x.to_rv(), x.to_rv()]));
    let d: MX = one(raw::minus(vec![y.to_rv(), x2.to_rv()]));
    let d2: MX = one(raw::times(vec![d.to_rv(), d.to_rv()]));
    let t: MX = one(raw::times(vec![100.0f64.to_rv(), d2.to_rv()]));
    let om: MX = one(raw::minus(vec![1.0f64.to_rv(), x.to_rv()]));
    let om2: MX = one(raw::times(vec![om.to_rv(), om.to_rv()]));
    let f: MX = one(raw::plus(vec![om2.to_rv(), t.to_rv()]));
    let xy: MX = one(raw::vertcat(vec![RV::list(vec![x.to_rv(), y.to_rv()])]));
    let nlp = RV::dict([("x", xy.to_rv()), ("f", f.to_rv())]);
    let opts = RV::dict([("print_time", false.to_rv()), ("ipopt.print_level", 0i64.to_rv())]);
    let solver: Function = one(raw::nlpsol(vec!["S".to_rv(), "ipopt".to_rv(), nlp, opts]));
    let arg = RV::dict([("x0", vec![-1.0f64, 1.0].to_rv())]);
    let res: HashMap<String, DM> = one(raw::Function_call(vec![solver.to_rv(), arg]));
    let xopt: Vec<f64> = one(raw::DM_nonzeros(vec![res["x"].to_rv()]));
    assert!((xopt[0] - 1.0).abs() < 1e-6 && (xopt[1] - 1.0).abs() < 1e-6, "{:?}", xopt);
    let stats: HashMap<String, rt::Value> = one(raw::Function_stats(vec![solver.to_rv()]));
    let _ = stats;
}

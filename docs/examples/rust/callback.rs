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
// How to use Callback
// Joel Andersson

use casadi as ca;
use casadi::dict;
use std::sync::Mutex;

struct MyCallback {
    d: f64,
}

impl ca::CallbackImpl for MyCallback {
    // Number of inputs and outputs
    fn get_n_in(&self) -> i64 {
        1
    }
    fn get_n_out(&self) -> i64 {
        1
    }

    // Initialize the object
    fn init(&self) {
        println!("initializing object");
    }

    // Evaluate numerically
    fn eval(&self, arg: Vec<ca::DM>) -> Vec<ca::DM> {
        let x = &arg[0];
        let f = ca::sin(&(self.d * x));
        vec![f]
    }
}

fn my_callback(name: &str, d: f64, opts: ca::Dict) -> ca::Callback {
    let cb = ca::Callback::from_impl(MyCallback { d });
    cb.construct_with(name, opts);
    cb
}

// Derivates OPTION 2-4: a 4-by-3 example, with forward, reverse or Jacobian information
// Example from https://www.youtube.com/watch?v=mYOkLkS5yqc&t=4s
fn split4(x: &ca::DM) -> (ca::DM, ca::DM, ca::DM, ca::DM) {
    let v = ca::vertsplit(x);
    (v[0].clone(), v[1].clone(), v[2].clone(), v[3].clone())
}

fn example_eval(arg: &[ca::DM]) -> Vec<ca::DM> {
    let (a, b, c, d) = split4(&arg[0]);
    let ret = ca::vertcat(vec![ca::sin(&c) * &d + ca::power(&d, 2.0), 2.0 * &a + &c, ca::power(&b, 2.0) + 5.0 * &c]);
    vec![ret]
}

// Shared by the Example4To3 variants, like the Python base class
macro_rules! example_4_to_3 {
    () => {
        fn get_n_in(&self) -> i64 {
            1
        }
        fn get_n_out(&self) -> i64 {
            1
        }
        fn get_sparsity_in(&self, _i: i64) -> ca::Sparsity {
            ca::Sparsity::dense_with(4, 1)
        }
        fn get_sparsity_out(&self, _i: i64) -> ca::Sparsity {
            ca::Sparsity::dense_with(3, 1)
        }
        // Evaluate numerically
        fn eval(&self, arg: Vec<ca::DM>) -> Vec<ca::DM> {
            example_eval(&arg)
        }
    };
}

struct ForwardFun;

impl ca::CallbackImpl for ForwardFun {
    fn get_n_in(&self) -> i64 {
        3
    }
    fn get_n_out(&self) -> i64 {
        1
    }
    fn get_sparsity_in(&self, i: i64) -> ca::Sparsity {
        match i {
            0 => ca::Sparsity::dense_with(4, 1), // nominal input
            1 => ca::Sparsity::new(3, 1),        // nominal output
            _ => ca::Sparsity::dense_with(4, 1), // Forward seed
        }
    }
    fn get_sparsity_out(&self, _i: i64) -> ca::Sparsity {
        // Forward sensitivity
        ca::Sparsity::dense_with(3, 1)
    }
    // Evaluate numerically
    fn eval(&self, arg: Vec<ca::DM>) -> Vec<ca::DM> {
        let (a, b, c, d) = split4(&arg[0]);
        let (a_dot, b_dot, c_dot, d_dot) = split4(&arg[2]);
        println!("Forward sweep with {} {} {} {}", a_dot, b_dot, c_dot, d_dot);
        let w0 = ca::sin(&c);
        let w0_dot = ca::cos(&c) * &c_dot;
        let w1 = &w0 * &d;
        let w1_dot = &w0_dot * &d + &w0 * &d_dot;
        let w2 = ca::power(&d, 2.0);
        let w2_dot = 2.0 * &d_dot * &d;
        let _r0 = &w1 + &w2;
        let r0_dot = &w1_dot + &w2_dot;
        let w3 = 2.0 * &a;
        let w3_dot = 2.0 * &a_dot;
        let _r1 = &w3 + &c;
        let r1_dot = &w3_dot + &c_dot;
        let w4 = ca::power(&b, 2.0);
        let w4_dot = 2.0 * &b_dot * &b;
        let w5 = 5.0 * &w0;
        let w5_dot = 5.0 * &w0_dot;
        let _r2 = &w4 + &w5;
        let r2_dot = &w4_dot + &w5_dot;
        vec![ca::vertcat(vec![r0_dot, r1_dot, r2_dot])]
    }
}

struct Example4To3Fwd {
    // You are required to keep the returned Callback object alive
    fwd_callback: Mutex<Option<ca::Callback>>,
}

impl ca::CallbackImpl for Example4To3Fwd {
    example_4_to_3!();
    fn has_forward(&self, nfwd: i64) -> bool {
        // This example is written to work with a single forward seed vector
        // For efficiency, you may allow more seeds at once
        nfwd == 1
    }
    fn get_forward(&self, _nfwd: i64, name: String, _inames: Vec<String>, _onames: Vec<String>, opts: std::collections::HashMap<String, ca::Value>) -> ca::Function {
        let cb = ca::Callback::from_impl(ForwardFun);
        cb.construct_with(&name, opts);
        let f = ca::Function::from(&cb);
        *self.fwd_callback.lock().unwrap() = Some(cb);
        f
    }
}

struct ReverseFun;

impl ca::CallbackImpl for ReverseFun {
    fn get_n_in(&self) -> i64 {
        3
    }
    fn get_n_out(&self) -> i64 {
        1
    }
    fn get_sparsity_in(&self, i: i64) -> ca::Sparsity {
        match i {
            0 => ca::Sparsity::dense_with(4, 1), // nominal input
            1 => ca::Sparsity::new(3, 1),        // nominal output
            _ => ca::Sparsity::dense_with(3, 1), // Reverse seed
        }
    }
    fn get_sparsity_out(&self, _i: i64) -> ca::Sparsity {
        // Reverse sensitivity
        ca::Sparsity::dense_with(4, 1)
    }
    // Evaluate numerically
    fn eval(&self, arg: Vec<ca::DM>) -> Vec<ca::DM> {
        let (a, b, c, d) = split4(&arg[0]);
        let r_bar = ca::vertsplit(&arg[2]);
        let (r0_bar, r1_bar, r2_bar) = (&r_bar[0], &r_bar[1], &r_bar[2]);
        println!("Reverse sweep with {} {} {}", r0_bar, r1_bar, r2_bar);
        let w0 = ca::sin(&c);
        let _w4 = ca::power(&b, 2.0);
        let w4_bar = r2_bar;
        let w5_bar = r2_bar;
        let w0_bar = 5.0 * w5_bar;
        let b_bar = 2.0 * &b * w4_bar;
        let w3_bar = r1_bar;
        let c_bar = r1_bar.clone();
        let a_bar = 2.0 * w3_bar;
        let w1_bar = r0_bar;
        let w2_bar = r0_bar;
        let d_bar = 2.0 * &d * w2_bar;
        let w0_bar = &w0_bar + w1_bar * &d;
        let d_bar = &d_bar + &w0 * w1_bar;
        let c_bar = &c_bar + ca::cos(&c) * &w0_bar;
        let _ = &a;
        vec![ca::vertcat(vec![a_bar, b_bar, c_bar, d_bar])]
    }
}

struct Example4To3Rev {
    rev_callback: Mutex<Option<ca::Callback>>,
}

impl ca::CallbackImpl for Example4To3Rev {
    example_4_to_3!();
    fn has_reverse(&self, nadj: i64) -> bool {
        nadj == 1
    }
    fn get_reverse(&self, _nadj: i64, name: String, _inames: Vec<String>, _onames: Vec<String>, opts: std::collections::HashMap<String, ca::Value>) -> ca::Function {
        let cb = ca::Callback::from_impl(ReverseFun);
        cb.construct_with(&name, opts);
        let f = ca::Function::from(&cb);
        *self.rev_callback.lock().unwrap() = Some(cb);
        f
    }
}

struct JacFun;

impl ca::CallbackImpl for JacFun {
    fn get_n_in(&self) -> i64 {
        2
    }
    fn get_n_out(&self) -> i64 {
        1
    }
    fn get_sparsity_in(&self, i: i64) -> ca::Sparsity {
        match i {
            0 => ca::Sparsity::dense_with(4, 1), // nominal input
            _ => ca::Sparsity::new(3, 1),        // nominal output
        }
    }
    fn get_sparsity_out(&self, _i: i64) -> ca::Sparsity {
        ca::sparsify(&ca::DM::new(vec![vec![0.0, 0.0, 1.0, 1.0], vec![1.0, 0.0, 1.0, 0.0], vec![0.0, 1.0, 1.0, 0.0]])).sparsity()
    }
    // Evaluate numerically
    fn eval(&self, arg: Vec<ca::DM>) -> Vec<ca::DM> {
        let (_a, b, c, d) = split4(&arg[0]);
        let mut ret = ca::DM::zeros_with(3, 4);
        ret.set_at((0, 2), &d * ca::cos(&c));
        ret.set_at((0, 3), ca::sin(&c) + 2.0 * &d);
        ret.set_at((1, 0), 2.0);
        ret.set_at((1, 2), 1.0);
        ret.set_at((2, 1), 2.0 * &b);
        ret.set_at((2, 2), 5.0);
        vec![ret]
    }
}

struct Example4To3Jac {
    jac_callback: Mutex<Option<ca::Callback>>,
}

impl ca::CallbackImpl for Example4To3Jac {
    example_4_to_3!();
    fn has_jacobian(&self) -> bool {
        true
    }
    fn get_jacobian(&self, name: String, _inames: Vec<String>, _onames: Vec<String>, opts: std::collections::HashMap<String, ca::Value>) -> ca::Function {
        let cb = ca::Callback::from_impl(JacFun);
        cb.construct_with(&name, opts);
        let f = ca::Function::from(&cb);
        *self.jac_callback.lock().unwrap() = Some(cb);
        f
    }
}

fn jacobian_at(f: &ca::Callback) -> ca::DM {
    let x = ca::MX::sym_with("x", 4, None);
    let fx: ca::MX = f.eval([&x]);
    let j = ca::Function::new("J", vec![&x], vec![ca::jacobian(&fx, &x)]);
    j.eval([ca::DM::new(vec![1.0, 2.0, 0.0, 3.0])])
}

fn main() {
    // Use the function
    let f = my_callback("f", 0.5, dict! {});
    let res: ca::DM = f.eval([2.0]);
    println!("{}", res);

    // You may call the Callback symbolically
    let x = ca::MX::sym("x");
    let fx: ca::MX = f.eval([&x]);
    println!("{}", fx);

    // Derivates OPTION 1: finite-differences
    let eps = 1e-5;
    let f2: ca::DM = f.eval([2.0 + eps]);
    println!("{}", (f2 - &res) / eps);

    let f = my_callback("f", 0.5, dict! {"enable_fd" => true});
    let fx: ca::MX = f.eval([&x]);
    let j = ca::Function::new("J", vec![&x], vec![ca::jacobian(&fx, &x)]);
    let r: ca::DM = j.eval([2.0]);
    println!("{}", r);

    // Derivates OPTION 2: Supply forward mode
    let f = ca::Callback::from_impl(Example4To3Fwd { fwd_callback: Mutex::new(None) });
    f.construct("f");
    println!("{}", jacobian_at(&f));

    // Derivates OPTION 3: Supply reverse mode
    let f = ca::Callback::from_impl(Example4To3Rev { rev_callback: Mutex::new(None) });
    f.construct("f");
    println!("{}", jacobian_at(&f));

    // Derivates OPTION 4: Supply full Jacobian
    let f = ca::Callback::from_impl(Example4To3Jac { jac_callback: Mutex::new(None) });
    f.construct("f");
    println!("{}", jacobian_at(&f));
}

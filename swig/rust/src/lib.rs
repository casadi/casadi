//! CasADi: nonlinear optimization and algorithmic differentiation
#![allow(clippy::all, non_camel_case_types, non_snake_case)]

include!(concat!(env!("CASADI_RUST_GEN"), "/casadi.rs"));

mod catch;
mod conv;
mod eval;
mod fmt;
mod index;
mod load;
mod matrix;
mod ops;
mod opti;

pub use catch::catch;
pub use eval::{Args, Elem, EvalArgs, EvalOut, Target};
pub use rt::{Dict, Error, Value};
pub use index::{step, At, Idx, Index};
pub use matrix::Matrix;

// The generated code resolves the wrapper's symbols through this
fn swig_resolve(name: &str) -> *mut std::ffi::c_void {
    load::resolve(name)
}

// casadi is built with atomic reference counts and thread-safe evaluation
unsafe impl Send for Function {}
unsafe impl Sync for Function {}
unsafe impl Send for Callback {}
unsafe impl Sync for Callback {}

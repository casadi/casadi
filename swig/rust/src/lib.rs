//! CasADi: nonlinear optimization and algorithmic differentiation
#![allow(clippy::all, non_camel_case_types, non_snake_case)]

include!(concat!(env!("CASADI_RUST_GEN"), "/casadi.rs"));

mod eval;
mod fmt;
mod index;
mod ops;

pub use eval::{Args, Elem, EvalArgs, EvalOut, Target};
pub use rt::{Dict, Value};
pub use index::{step, At, Idx, Index};

// casadi is built with atomic reference counts and thread-safe evaluation
unsafe impl Send for Function {}
unsafe impl Sync for Function {}

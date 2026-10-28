//! CasADi: nonlinear optimization and algorithmic differentiation
#![allow(clippy::all, non_camel_case_types, non_snake_case)]

include!(concat!(env!("CASADI_RUST_GEN"), "/casadi.rs"));

mod eval;
mod fmt;
mod ops;

pub use eval::{EvalArgs, EvalOut};

// casadi is built with atomic reference counts and thread-safe evaluation
unsafe impl Send for Function {}
unsafe impl Sync for Function {}

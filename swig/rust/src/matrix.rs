//! Code over DM, SX and MX alike, as Python's duck typing: `fn f<M: Matrix>(w: &M) -> M`
use crate::{rt, At, DM, MX, SX};
use std::fmt;
use std::ops::{Add, Div, Mul, Neg, Sub};

/// What DM, SX and MX have in common: elementwise operators (also with `&M` and `f64` on the
/// right), constructors and concatenations as associated functions, math as methods
pub trait Matrix:
    Clone
    + Default
    + fmt::Display
    + rt::ToRV
    + From<f64>
    + From<Vec<f64>>
    + Add<Output = Self>
    + Sub<Output = Self>
    + Mul<Output = Self>
    + Div<Output = Self>
    + Neg<Output = Self>
    + for<'a> Add<&'a Self, Output = Self>
    + for<'a> Sub<&'a Self, Output = Self>
    + for<'a> Mul<&'a Self, Output = Self>
    + for<'a> Div<&'a Self, Output = Self>
    + Add<f64, Output = Self>
    + Sub<f64, Output = Self>
    + Mul<f64, Output = Self>
    + Div<f64, Output = Self>
{
    fn zeros(nrow: i64, ncol: i64) -> Self;
    fn ones(nrow: i64, ncol: i64) -> Self;
    fn vertcat(v: Vec<Self>) -> Self;
    fn horzcat(v: Vec<Self>) -> Self;
    fn if_else(cond: &Self, if_true: &Self, if_false: &Self) -> Self;

    fn size1(&self) -> i64;
    fn size2(&self) -> i64;
    fn numel(&self) -> i64;
    /// Python's x[idx]
    fn at(&self, idx: impl At) -> Self;
    fn T(&self) -> Self;
    fn reshape(&self, nrow: i64, ncol: i64) -> Self;
    fn mtimes(&self, y: &Self) -> Self;
    fn dot(&self, y: &Self) -> Self;
    fn power(&self, y: &Self) -> Self;
    fn fmin(&self, y: &Self) -> Self;
    fn fmax(&self, y: &Self) -> Self;
    fn atan2(&self, x: &Self) -> Self;

    fn sin(&self) -> Self;
    fn cos(&self) -> Self;
    fn tan(&self) -> Self;
    fn asin(&self) -> Self;
    fn acos(&self) -> Self;
    fn atan(&self) -> Self;
    fn sinh(&self) -> Self;
    fn cosh(&self) -> Self;
    fn tanh(&self) -> Self;
    fn exp(&self) -> Self;
    fn log(&self) -> Self;
    fn sqrt(&self) -> Self;
    fn abs(&self) -> Self;
    fn floor(&self) -> Self;
    fn ceil(&self) -> Self;
    fn sign(&self) -> Self;
    fn sum1(&self) -> Self;
    fn sum2(&self) -> Self;
    fn sumsqr(&self) -> Self;
    fn norm_1(&self) -> Self;
    fn norm_2(&self) -> Self;
    fn norm_inf(&self) -> Self;
    fn mmin(&self) -> Self;
    fn mmax(&self) -> Self;
}

macro_rules! matrix {
    ($T:ident; $($un:ident)*; $($bin:ident)*) => {
        impl Matrix for $T {
            fn zeros(nrow: i64, ncol: i64) -> Self { $T::zeros_with(nrow, ncol) }
            fn ones(nrow: i64, ncol: i64) -> Self { $T::ones_with(nrow, ncol) }
            fn vertcat(v: Vec<Self>) -> Self { crate::vertcat(v) }
            fn horzcat(v: Vec<Self>) -> Self { crate::horzcat(v) }
            fn if_else(cond: &Self, if_true: &Self, if_false: &Self) -> Self { crate::if_else(cond, if_true, if_false) }
            fn size1(&self) -> i64 { $T::size1(self) }
            fn size2(&self) -> i64 { $T::size2(self) }
            fn numel(&self) -> i64 { $T::numel(self) }
            #[track_caller]
            fn at(&self, idx: impl At) -> Self { $T::at(self, idx) }
            fn T(&self) -> Self { $T::T(self) }
            fn reshape(&self, nrow: i64, ncol: i64) -> Self { crate::reshape(self, nrow, ncol) }
            fn mtimes(&self, y: &Self) -> Self { crate::mtimes(self, y) }
            $(fn $bin(&self, y: &Self) -> Self { crate::$bin(self, y) })*
            $(fn $un(&self) -> Self { crate::$un(self) })*
        }
    };
}

macro_rules! matrices {
    ($($T:ident),*) => {$(
        matrix!($T; sin cos tan asin acos atan sinh cosh tanh exp log sqrt abs floor ceil sign
                sum1 sum2 sumsqr norm_1 norm_2 norm_inf mmin mmax; dot power fmin fmax atan2);
    )*};
}

matrices!(DM, SX, MX);

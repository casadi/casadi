//! Elementwise arithmetic operators; mtimes is the matrix product
use crate::{minus, plus, rdivide, times, DM, MX, SX};
use std::ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Neg, Sub, SubAssign};

macro_rules! binop {
    ($L:ty, $R:ty => $O:ty) => {
        binop!(@ Add add plus $L, $R => $O);
        binop!(@ Sub sub minus $L, $R => $O);
        binop!(@ Mul mul times $L, $R => $O);
        binop!(@ Div div rdivide $L, $R => $O);
    };
    (@ $Tr:ident $m:ident $f:ident $L:ty, $R:ty => $O:ty) => {
        impl $Tr<&$R> for &$L {
            type Output = $O;
            fn $m(self, r: &$R) -> $O {
                $f(self, r)
            }
        }
        impl $Tr<$R> for &$L {
            type Output = $O;
            fn $m(self, r: $R) -> $O {
                $f(self, &r)
            }
        }
        impl $Tr<&$R> for $L {
            type Output = $O;
            fn $m(self, r: &$R) -> $O {
                $f(&self, r)
            }
        }
        impl $Tr<$R> for $L {
            type Output = $O;
            fn $m(self, r: $R) -> $O {
                $f(&self, &r)
            }
        }
    };
}

macro_rules! assignop {
    ($L:ty, $R:ty) => {
        assignop!(@ AddAssign add_assign plus $L, $R);
        assignop!(@ SubAssign sub_assign minus $L, $R);
        assignop!(@ MulAssign mul_assign times $L, $R);
        assignop!(@ DivAssign div_assign rdivide $L, $R);
    };
    (@ $Tr:ident $m:ident $f:ident $L:ty, $R:ty) => {
        impl $Tr<&$R> for $L {
            fn $m(&mut self, r: &$R) {
                *self = $f(&*self, r);
            }
        }
        impl $Tr<$R> for $L {
            fn $m(&mut self, r: $R) {
                *self = $f(&*self, &r);
            }
        }
    };
}

macro_rules! neg {
    ($T:ty) => {
        impl Neg for &$T {
            type Output = $T;
            fn neg(self) -> $T {
                self.uminus()
            }
        }
        impl Neg for $T {
            type Output = $T;
            fn neg(self) -> $T {
                self.uminus()
            }
        }
    };
}

binop!(DM, DM => DM);
binop!(DM, f64 => DM);
binop!(f64, DM => DM);
binop!(SX, SX => SX);
binop!(SX, DM => SX);
binop!(DM, SX => SX);
binop!(SX, f64 => SX);
binop!(f64, SX => SX);
binop!(MX, MX => MX);
binop!(MX, DM => MX);
binop!(DM, MX => MX);
binop!(MX, f64 => MX);
binop!(f64, MX => MX);

assignop!(DM, DM);
assignop!(DM, f64);
assignop!(SX, SX);
assignop!(SX, DM);
assignop!(SX, f64);
assignop!(MX, MX);
assignop!(MX, DM);
assignop!(MX, f64);

neg!(DM);
neg!(SX);
neg!(MX);

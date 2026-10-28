//! Display and Debug through casadi's printing
use crate::{Function, Sparsity, DM, MX, SX};
use std::fmt;

macro_rules! display {
    ($($T:ty),*) => {$(
        impl fmt::Display for $T {
            fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
                f.write_str(&self.str())
            }
        }
        impl fmt::Debug for $T {
            fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
                f.write_str(&self.str())
            }
        }
    )*};
}

display!(DM, SX, MX, Sparsity, Function);

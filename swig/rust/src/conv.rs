//! Conversions to matrices, so `.into()` works where Python takes a number or a list
use crate::{DM, MX, SX};

macro_rules! from {
    ($T:ty: $($S:ty),*) => {$(
        impl From<$S> for $T {
            fn from(x: $S) -> $T {
                <$T>::new(x)
            }
        }
    )*};
}

from!(DM: f64, Vec<f64>, Vec<Vec<f64>>);
from!(SX: f64, Vec<f64>, Vec<Vec<f64>>, DM, &DM);
from!(MX: f64, Vec<f64>, Vec<Vec<f64>>, DM, &DM);

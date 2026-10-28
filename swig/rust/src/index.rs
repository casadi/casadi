//! Indexing: x.at(i), x.at((.., j)), x.at(1..3), x.set_at(i, v), x.nz(k)
use crate::{ffi, rt, Slice, DM, MX, SX};
use rt::{ToRV, RV};
use std::ops::{Range, RangeFrom, RangeFull, RangeInclusive, RangeTo, RangeToInclusive};

/// One index: an integer (negative counts from the end), a range or a list
pub enum Index {
    Int(i64),
    Slice(Slice),
    List(Vec<i64>),
}

impl ToRV for Index {
    fn to_rv(&self) -> RV {
        match self {
            Index::Int(i) => i.to_rv(),
            Index::Slice(s) => s.to_rv(),
            Index::List(l) => l.to_rv(),
        }
    }
}

fn slice(start: i64, stop: i64) -> Index {
    Index::Slice(Slice::new(start, stop))
}

/// Types usable as one index
pub trait Idx {
    fn index(&self) -> Index;
}

/// Types usable as a matrix index: one index (linear) or a (row, column) pair
pub trait At {
    fn indices(&self) -> Vec<Index>;
}

macro_rules! int_idx {
    ($($T:ty),*) => {$(
        impl Idx for $T {
            fn index(&self) -> Index {
                Index::Int(*self as i64)
            }
        }
        impl Idx for Range<$T> {
            fn index(&self) -> Index {
                slice(self.start as i64, self.end as i64)
            }
        }
        impl Idx for RangeInclusive<$T> {
            fn index(&self) -> Index {
                slice(*self.start() as i64, *self.end() as i64 + 1)
            }
        }
        impl Idx for RangeFrom<$T> {
            fn index(&self) -> Index {
                slice(self.start as i64, i64::MAX)
            }
        }
        impl Idx for RangeTo<$T> {
            fn index(&self) -> Index {
                slice(0, self.end as i64)
            }
        }
        impl Idx for RangeToInclusive<$T> {
            fn index(&self) -> Index {
                slice(0, self.end as i64 + 1)
            }
        }
        impl Idx for Vec<$T> {
            fn index(&self) -> Index {
                Index::List(self.iter().map(|&i| i as i64).collect())
            }
        }
        impl<const N: usize> Idx for [$T; N] {
            fn index(&self) -> Index {
                Index::List(self.iter().map(|&i| i as i64).collect())
            }
        }
        at_one!(Range<$T>, RangeInclusive<$T>, RangeFrom<$T>, RangeTo<$T>, RangeToInclusive<$T>, $T, Vec<$T>);
        impl<const N: usize> At for [$T; N] {
            fn indices(&self) -> Vec<Index> {
                vec![self.index()]
            }
        }
    )*};
}

macro_rules! at_one {
    ($($T:ty),*) => {$(
        impl At for $T {
            fn indices(&self) -> Vec<Index> {
                vec![self.index()]
            }
        }
    )*};
}

impl Idx for RangeFull {
    fn index(&self) -> Index {
        Index::Slice(Slice::default())
    }
}

impl Idx for Slice {
    fn index(&self) -> Index {
        Index::Slice(self.clone())
    }
}

at_one!(RangeFull, Slice);
int_idx!(i32, i64, usize);

impl<A: Idx, B: Idx> At for (A, B) {
    fn indices(&self) -> Vec<Index> {
        vec![self.0.index(), self.1.index()]
    }
}

fn args(this: &impl ToRV, ind: &[Index]) -> Vec<RV> {
    let mut a = vec![this.to_rv(), false.to_rv()];
    a.extend(ind.iter().map(|i| i.to_rv()));
    a
}

macro_rules! indexing {
    ($T:ident, $get:ident, $set:ident, $get_nz:ident, $set_nz:ident) => {
        impl $T {
            /// Python's x[idx]: a linear index or a (row, column) pair
            pub fn at(&self, idx: impl At) -> $T {
                let ind = idx.indices();
                rt::out1(rt::call_or_panic(ffi::$get, args(self, &ind), 1))
            }
            /// Python's x[idx] = v
            pub fn set_at(&mut self, idx: impl At, v: impl ToRV) {
                let ind = idx.indices();
                let mut a = vec![self.to_rv(), v.to_rv(), false.to_rv()];
                a.extend(ind.iter().map(|i| i.to_rv()));
                rt::call_or_panic(ffi::$set, a, 1);
            }
            /// Python's x.nz[k]
            pub fn nz(&self, k: impl Idx) -> $T {
                let ind = [k.index()];
                rt::out1(rt::call_or_panic(ffi::$get_nz, args(self, &ind), 1))
            }
            /// Python's x.nz[k] = v
            pub fn set_nz_at(&mut self, k: impl Idx, v: impl ToRV) {
                let ind = k.index();
                rt::call_or_panic(ffi::$set_nz, vec![self.to_rv(), v.to_rv(), false.to_rv(), ind.to_rv()], 1);
            }
        }
    };
}

indexing!(DM, swigrs_DM_get, swigrs_DM_set, swigrs_DM_get_nz, swigrs_DM_set_nz);
indexing!(SX, swigrs_SX_get, swigrs_SX_set, swigrs_SX_get_nz, swigrs_SX_set_nz);
indexing!(MX, swigrs_MX_get, swigrs_MX_set, swigrs_MX_get_nz, swigrs_MX_set_nz);

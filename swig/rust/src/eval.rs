//! Function evaluation like Python's f(a, b) and f(x0=a, p=b)
use crate::{ffi, raw, rt, Callback, Function, DM, MX, SX};
use rt::ToRV;
use std::collections::{BTreeMap, HashMap};

/// How a Function is evaluated: numerically or on SX/MX expressions
#[derive(Clone, Copy, PartialEq)]
pub enum Target {
    Num,
    SX,
    MX,
}

/// Output element types; they select the evaluation
pub trait Elem: rt::FromRV {
    const TARGET: Target;
}

impl Elem for DM {
    const TARGET: Target = Target::Num;
}
impl Elem for SX {
    const TARGET: Target = Target::SX;
}
impl Elem for MX {
    const TARGET: Target = Target::MX;
}

/// Arguments of `Function::eval`
pub enum Args {
    Positional(Vec<rt::RV>),
    Named(Vec<(String, rt::RV)>),
}

/// Positional (array, Vec, slice, tuple) or named (`dict!`, HashMap, BTreeMap) arguments
pub trait EvalArgs {
    fn to_args(&self) -> Args;
}

impl<T: ToRV, const N: usize> EvalArgs for [T; N] {
    fn to_args(&self) -> Args {
        Args::Positional(self.iter().map(|a| a.to_rv()).collect())
    }
}

impl<T: ToRV> EvalArgs for [T] {
    fn to_args(&self) -> Args {
        Args::Positional(self.iter().map(|a| a.to_rv()).collect())
    }
}

impl<T: ToRV> EvalArgs for Vec<T> {
    fn to_args(&self) -> Args {
        self.as_slice().to_args()
    }
}

impl<A: EvalArgs + ?Sized> EvalArgs for &A {
    fn to_args(&self) -> Args {
        (**self).to_args()
    }
}

impl EvalArgs for rt::Dict {
    fn to_args(&self) -> Args {
        Args::Named(self.0.iter().map(|(k, v)| (k.clone(), v.to_rv())).collect())
    }
}

impl<K: AsRef<str>, V: ToRV> EvalArgs for HashMap<K, V> {
    fn to_args(&self) -> Args {
        Args::Named(self.iter().map(|(k, v)| (k.as_ref().to_string(), v.to_rv())).collect())
    }
}

impl<K: AsRef<str>, V: ToRV> EvalArgs for BTreeMap<K, V> {
    fn to_args(&self) -> Args {
        Args::Named(self.iter().map(|(k, v)| (k.as_ref().to_string(), v.to_rv())).collect())
    }
}

/// Result of `Function::eval`: one matrix or number for one output, a tuple or Vec otherwise,
/// a HashMap for named arguments
pub trait EvalOut: Sized {
    const TARGET: Target;
    fn from_outputs(o: rt::RVRef) -> Option<Self>;
}

macro_rules! eval_out_single {
    ($($T:ty),*) => {$(
        impl EvalOut for $T {
            const TARGET: Target = <$T as Elem>::TARGET;
            fn from_outputs(o: rt::RVRef) -> Option<Self> {
                assert_eq!(o.kind(), rt::LIST, "named arguments give a HashMap");
                assert_eq!(o.len(), 1, "Function has {} outputs, a tuple or Vec is needed", o.len());
                <$T as rt::FromRV>::from_rv(o.item(0))
            }
        }
    )*};
}

eval_out_single!(DM, SX, MX);

impl EvalOut for f64 {
    const TARGET: Target = Target::Num;
    fn from_outputs(o: rt::RVRef) -> Option<Self> {
        DM::from_outputs(o).map(|d| d.to_double())
    }
}

impl<T: Elem> EvalOut for Vec<T> {
    const TARGET: Target = T::TARGET;
    fn from_outputs(o: rt::RVRef) -> Option<Self> {
        rt::FromRV::from_rv(o)
    }
}

impl<T: Elem> EvalOut for HashMap<String, T> {
    const TARGET: Target = T::TARGET;
    fn from_outputs(o: rt::RVRef) -> Option<Self> {
        assert_eq!(o.kind(), rt::DICT, "positional arguments give a matrix, tuple or Vec");
        rt::FromRV::from_rv(o)
    }
}

macro_rules! tuples {
    ($n:expr; $A0:ident $($A:ident $i:tt),*) => {
        impl<$A0: ToRV, $($A: ToRV),*> EvalArgs for ($A0, $($A,)*) {
            fn to_args(&self) -> Args {
                Args::Positional(vec![self.0.to_rv(), $(self.$i.to_rv()),*])
            }
        }
        impl<$A0: Elem, $($A: Elem),*> EvalOut for ($A0, $($A,)*) {
            const TARGET: Target = $A0::TARGET;
            fn from_outputs(o: rt::RVRef) -> Option<Self> {
                assert_eq!(o.kind(), rt::LIST, "named arguments give a HashMap");
                assert_eq!(o.len(), $n, "Function has {} outputs, not {}", o.len(), $n);
                Some(($A0::from_rv(o.item(0))?, $($A::from_rv(o.item($i))?,)*))
            }
        }
    };
}

tuples!(1; A);
tuples!(2; A B 1);
tuples!(3; A B 1, C 2);
tuples!(4; A B 1, C 2, D 3);
tuples!(5; A B 1, C 2, D 3, E 4);
tuples!(6; A B 1, C 2, D 3, E 4, F 5);
tuples!(7; A B 1, C 2, D 3, E 4, F 5, G 6);
tuples!(8; A B 1, C 2, D 3, E 4, F 5, G 6, H 7);

/// Converts an argument to the symbolic type of the evaluation
fn coerce(v: rt::RV, target: Target) -> rt::RV {
    let (ty, new): (_, fn(Vec<rt::RV>) -> Result<Vec<rt::RV>, rt::Error>) = match target {
        Target::Num => return v,
        Target::SX => (SX::type_info(), raw::new_SX),
        Target::MX => (MX::type_info(), raw::new_MX),
    };
    if v.kind() == rt::PTR && !unsafe { rt::swig_rv_cast_ptr(v.0, ty) }.is_null() {
        return v;
    }
    match new(vec![v]) {
        Ok(mut o) => o.remove(0),
        Err(e) => panic!("{}", e),
    }
}

macro_rules! eval_impl {
    ($($T:ty),*) => {$(
impl $T {
    /// Evaluate like Python's `f(a, b)` or `f(x0=a, p=b)`; the result type selects numeric, SX or MX
    /// evaluation, and arguments are converted to it
    ///
    /// Panics on errors; `call` returns a `Result` for numeric evaluation.
    pub fn eval<T: EvalOut>(&self, args: impl EvalArgs) -> T {
        // args outlives the call: the converted values borrow its objects
        let arg = match args.to_args() {
            Args::Positional(v) => rt::RV::list(v.into_iter().map(|a| coerce(a, T::TARGET)).collect()),
            Args::Named(v) => {
                let (keys, vals): (Vec<String>, Vec<rt::RV>) = v.into_iter().map(|(k, a)| (k, coerce(a, T::TARGET))).unzip();
                rt::RV::dict(keys.iter().map(|k| k.as_str()).zip(vals))
            }
        };
        let f = if T::TARGET == Target::Num { ffi::swigrs_Function_call } else { ffi::swigrs_Function_call_symbolic };
        let o = rt::call_or_panic(f, vec![self.to_rv(), arg], 1);
        drop(args);
        let o = o.into_iter().next().expect("missing output");
        let desc = o.type_desc();
        match T::from_outputs(rt::RVRef(o.0)) {
            Some(r) => r,
            None => panic!("cannot convert {} to {}", desc, std::any::type_name::<T>()),
        }
    }
}
    )*};
}

eval_impl!(Function, Callback);

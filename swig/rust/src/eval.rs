//! Function evaluation with positional arguments, like Python's f(a, b)
use crate::{ffi, rt, Function, DM};
use rt::ToRV;

/// Positional arguments of `Function::eval`: an array, a Vec or a tuple
pub trait EvalArgs {
    fn to_rvs(&self) -> Vec<rt::RV>;
}

impl<T: ToRV, const N: usize> EvalArgs for [T; N] {
    fn to_rvs(&self) -> Vec<rt::RV> {
        self.iter().map(|a| a.to_rv()).collect()
    }
}

impl<T: ToRV> EvalArgs for Vec<T> {
    fn to_rvs(&self) -> Vec<rt::RV> {
        self.iter().map(|a| a.to_rv()).collect()
    }
}

/// Result of `Function::eval`: a matrix for one output, a tuple or Vec otherwise
pub trait EvalOut: Sized {
    fn from_outputs(o: rt::RVRef) -> Option<Self>;
}

fn single<T: rt::FromRV>(o: rt::RVRef) -> Option<T> {
    assert_eq!(o.len(), 1, "Function has {} outputs, a tuple or Vec is needed", o.len());
    T::from_rv(o.item(0))
}

macro_rules! eval_out_single {
    ($($T:ty),*) => {$(
        impl EvalOut for $T {
            fn from_outputs(o: rt::RVRef) -> Option<Self> {
                single(o)
            }
        }
    )*};
}

eval_out_single!(DM, crate::SX, crate::MX);

impl<T: rt::FromRV> EvalOut for Vec<T> {
    fn from_outputs(o: rt::RVRef) -> Option<Self> {
        rt::FromRV::from_rv(o)
    }
}

macro_rules! tuples {
    ($n:expr; $($A:ident $i:tt),*) => {
        impl<$($A: ToRV),*> EvalArgs for ($($A,)*) {
            fn to_rvs(&self) -> Vec<rt::RV> {
                vec![$(self.$i.to_rv()),*]
            }
        }
        impl<$($A: rt::FromRV),*> EvalOut for ($($A,)*) {
            fn from_outputs(o: rt::RVRef) -> Option<Self> {
                assert_eq!(o.len(), $n, "Function has {} outputs, not {}", o.len(), $n);
                Some(($($A::from_rv(o.item($i))?,)*))
            }
        }
    };
}

tuples!(1; A 0);
tuples!(2; A 0, B 1);
tuples!(3; A 0, B 1, C 2);
tuples!(4; A 0, B 1, C 2, D 3);
tuples!(5; A 0, B 1, C 2, D 3, E 4);
tuples!(6; A 0, B 1, C 2, D 3, E 4, F 5);
tuples!(7; A 0, B 1, C 2, D 3, E 4, F 5, G 6);
tuples!(8; A 0, B 1, C 2, D 3, E 4, F 5, G 6, H 7);

impl Function {
    /// Evaluate with positional arguments; numeric if all are numeric, else symbolic
    ///
    /// Panics on errors; `call` returns a `Result` for numeric evaluation.
    pub fn eval<T: EvalOut>(&self, args: impl EvalArgs) -> T {
        // args outlives the call: the converted values borrow its objects
        let rvs = args.to_rvs();
        let dm = DM::type_info();
        let numeric = rvs.iter().all(|a| a.kind() != rt::PTR || unsafe { !rt::swig_rv_cast_ptr(a.0, dm).is_null() });
        let f = if numeric { ffi::swigrs_Function_call } else { ffi::swigrs_Function_call_symbolic };
        let o = rt::call_or_panic(f, vec![self.to_rv(), rt::RV::list(rvs)], 1);
        drop(args);
        let o = o.into_iter().next().expect("missing output");
        let desc = o.type_desc();
        match T::from_outputs(rt::RVRef(o.0)) {
            Some(r) => r,
            None => panic!("cannot convert {} to {}", desc, std::any::type_name::<T>()),
        }
    }
}

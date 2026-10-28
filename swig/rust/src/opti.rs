//! Opti iteration callbacks from closures, like Python's opti.callback(f)
use crate::{Opti, OptiCallback, OptiCallbackImpl};
use std::cell::RefCell;
use std::collections::HashMap;

// Opti is not Send, and solve() runs the callback on the calling thread, so the closure
// never leaves the thread that created it.
struct SameThread<F>(F);
unsafe impl<F> Send for SameThread<F> {}
unsafe impl<F> Sync for SameThread<F> {}

impl<F: Fn(i64) + 'static> OptiCallbackImpl for SameThread<F> {
    fn call(&self, i: i64) {
        (self.0)(i)
    }
}

thread_local! {
    // The callback of each Opti problem: Opti only keeps a raw pointer to it
    static CALLBACKS: RefCell<HashMap<i64, OptiCallback>> = RefCell::new(HashMap::new());
}

impl Opti {
    /// Runs f(i) after every solver iteration i, e.g. to print `opti.debug().value(&x)`
    ///
    /// The callback lives until it is replaced or `clear_callback()` is called; a closure
    /// that captures this Opti keeps the problem alive until then.
    pub fn callback(&self, f: impl Fn(i64) + 'static) {
        let cb = OptiCallback::from_impl(SameThread(f));
        self.callback_class_with(&cb);
        let old = CALLBACKS.with(|c| c.borrow_mut().insert(self.__hash__(), cb));
        drop(old);
    }

    /// Removes the iteration callback
    pub fn clear_callback(&self) {
        self.callback_class();
        let old = CALLBACKS.with(|c| c.borrow_mut().remove(&self.__hash__()));
        drop(old);
    }
}

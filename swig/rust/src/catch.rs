//! Errors as values, like Python's try/except around a call that raises
use crate::rt;
use std::cell::Cell;
use std::panic;
use std::sync::Once;

thread_local! {
    static QUIET: Cell<bool> = Cell::new(false);
}

/// Runs f and returns the error instead of panicking, e.g. `ca::catch(|| opti.solve())`
///
/// The panic message is not printed; panics on other threads are unaffected.
pub fn catch<T>(f: impl FnOnce() -> T) -> Result<T, rt::Error> {
    static HOOK: Once = Once::new();
    HOOK.call_once(|| {
        let prev = panic::take_hook();
        panic::set_hook(Box::new(move |info| {
            if !QUIET.with(|q| q.get()) {
                prev(info)
            }
        }));
    });
    let quiet = QUIET.with(|q| q.replace(true));
    let r = panic::catch_unwind(panic::AssertUnwindSafe(f));
    QUIET.with(|q| q.set(quiet));
    r.map_err(|e| {
        let message = match e.downcast::<String>() {
            Ok(s) => *s,
            Err(e) => match e.downcast::<&str>() {
                Ok(s) => s.to_string(),
                // Not an error message, e.g. a director's call of the C++ base
                Err(e) => panic::resume_unwind(e),
            },
        };
        rt::Error { code: -1, message }
    })
}

//! Loads libcasadi_rust_wrap on first use: $CASADI_RUST_LIB, else where casadi-sys put it
use libloading::Library;
use std::ffi::c_void;
use std::path::PathBuf;
use std::sync::OnceLock;

fn library() -> &'static Library {
    static LIB: OnceLock<Library> = OnceLock::new();
    LIB.get_or_init(|| {
        let path = std::env::var_os("CASADI_RUST_LIB").map(PathBuf::from).unwrap_or_else(|| {
            PathBuf::from(env!("CASADI_LIB_DIR")).join(format!("libcasadi_rust_wrap{}", std::env::consts::DLL_SUFFIX))
        });
        unsafe { open(&path) }.unwrap_or_else(|e| panic!("cannot load {}: {}", path.display(), e))
    })
}

#[cfg(unix)]
unsafe fn open(path: &PathBuf) -> Result<Library, libloading::Error> {
    Library::new(path)
}

// Dependent DLLs are found next to the wrapper, without PATH
#[cfg(windows)]
unsafe fn open(path: &PathBuf) -> Result<Library, libloading::Error> {
    use libloading::os::windows::{LOAD_LIBRARY_SEARCH_DEFAULT_DIRS, LOAD_LIBRARY_SEARCH_DLL_LOAD_DIR};
    libloading::os::windows::Library::load_with_flags(path, LOAD_LIBRARY_SEARCH_DLL_LOAD_DIR | LOAD_LIBRARY_SEARCH_DEFAULT_DIRS)
        .map(Library::from)
}

/// Address of a wrapper symbol, or null
pub(crate) fn resolve(name: &str) -> *mut c_void {
    unsafe { library().get::<*mut c_void>(name.as_bytes()).map(|s| *s).unwrap_or(std::ptr::null_mut()) }
}

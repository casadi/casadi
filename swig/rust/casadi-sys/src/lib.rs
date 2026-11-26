//! CasADi native libraries for the `casadi` crate
//!
//! The build script provides libcasadi, its plugins and libcasadi_rust_wrap, and tells the
//! `casadi` crate where they and the generated casadi.rs are. Nothing is linked: `casadi`
//! loads the wrapper library at run time. Sources, first match wins:
//! - `CASADI_RUST_GEN`: a cmake build tree with `-DWITH_RUST=ON` (`<build>/swig/rust`)
//! - an installed casadi directory around this crate, e.g. an unpacked casadi-<os>-rust.zip
//! - `CASADI_DIR`: such a casadi directory elsewhere
//! - casadi-source-v<version>.zip, for casadi.rs, plus the prebuilt casadi-<version>-<os>-rust.zip,
//!   or a cmake build of the source archive (feature `from-source`, or no prebuilt for the target).
//!   `CASADI_SOURCE_ZIP` and `CASADI_PREBUILT_ZIP` replace the downloads by a local file or URL,
//!   `CASADI_CMAKE_ARGS` adds cmake options to a source build.

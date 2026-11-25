fn main() {
    // Build tree: CASADI_RUST_GEN holds casadi.rs and the wrapper; installed: casadi.rs here, the libraries one up
    println!("cargo:rerun-if-env-changed=CASADI_RUST_GEN");
    let (gen, lib) = match std::env::var("CASADI_RUST_GEN") {
        Ok(gen) => (gen.clone(), gen),
        Err(_) => {
            let dir = std::path::PathBuf::from(std::env::var("CARGO_MANIFEST_DIR").unwrap());
            (dir.display().to_string(), dir.parent().unwrap().display().to_string())
        }
    };
    println!("cargo:rustc-env=CASADI_RUST_GEN={}", gen);
    println!("cargo:rerun-if-changed={}/casadi.rs", gen);
    println!("cargo:rustc-link-search=native={}", lib);
    println!("cargo:rustc-link-lib=dylib=casadi_rust_wrap");
    // Windows has no rpath: the DLLs are found on PATH
    if std::env::var("CARGO_CFG_TARGET_OS").unwrap() == "windows" {
        return;
    }
    println!("cargo:rustc-link-arg=-Wl,-rpath,{}", lib);
    if let Ok(lib) = std::env::var("CASADI_LIB") {
        println!("cargo:rustc-link-arg=-Wl,-rpath,{}", lib);
    }
}

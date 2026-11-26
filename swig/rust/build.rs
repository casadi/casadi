fn main() {
    // casadi-sys: where casadi.rs and libcasadi_rust_wrap are
    let gen = std::env::var("DEP_CASADI_RUST_GEN").expect("casadi-sys provided no casadi.rs (docs.rs builds offline)");
    let lib = std::env::var("DEP_CASADI_LIB").unwrap();
    println!("cargo:rustc-env=CASADI_RUST_GEN={}", gen);
    println!("cargo:rustc-env=CASADI_LIB_DIR={}", lib);
    println!("cargo:rerun-if-changed={}/casadi.rs", gen);
}

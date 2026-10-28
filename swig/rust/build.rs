fn main() {
    let gen = std::env::var("CASADI_RUST_GEN").expect("CASADI_RUST_GEN: directory with casadi.rs and the wrapper library");
    println!("cargo:rerun-if-env-changed=CASADI_RUST_GEN");
    println!("cargo:rerun-if-changed={}/casadi.rs", gen);
    println!("cargo:rustc-link-search=native={}", gen);
    println!("cargo:rustc-link-lib=dylib=casadi_rust_wrap");
    println!("cargo:rustc-link-arg=-Wl,-rpath,{}", gen);
    if let Ok(lib) = std::env::var("CASADI_LIB") {
        println!("cargo:rustc-link-arg=-Wl,-rpath,{}", lib);
    }
}

// Locates or provides the CasADi native libraries and casadi.rs; see src/lib.rs for the sources
use sha2::{Digest, Sha256};
use std::env;
use std::fs;
use std::io;
use std::path::{Path, PathBuf};

const RELEASES: &str = "https://github.com/casadi/casadi/releases/download";

// Plugin features of a source build
const PLUGINS: &[&str] = &["ipopt", "osqp", "highs", "daqp", "hpipm", "fatrop"];

fn main() {
    for v in ["CASADI_RUST_GEN", "CASADI_DIR", "CASADI_SOURCE_ZIP", "CASADI_PREBUILT_ZIP", "CASADI_CMAKE_ARGS"] {
        println!("cargo:rerun-if-env-changed={}", v);
    }
    println!("cargo:rerun-if-changed=checksums.txt");
    // docs.rs builds offline
    if env::var_os("DOCS_RS").is_some() {
        return;
    }
    let (gen, lib) = locate();
    println!("cargo:rust_gen={}", gen.display());
    println!("cargo:lib={}", lib.display());
}

/// Directory holding casadi.rs, and directory holding libcasadi_rust_wrap
fn locate() -> (PathBuf, PathBuf) {
    if let Some(gen) = env::var_os("CASADI_RUST_GEN") {
        let gen = PathBuf::from(gen);
        return (gen.clone(), gen);
    }
    // Installed: <casadi>/rust/casadi-sys
    let manifest = PathBuf::from(env::var("CARGO_MANIFEST_DIR").unwrap());
    if let Some(casadi) = manifest.parent().and_then(Path::parent) {
        if casadi.join(wrapper()).exists() {
            return (casadi.join("rust"), casadi.to_path_buf());
        }
    }
    if let Some(casadi) = env::var_os("CASADI_DIR") {
        let casadi = PathBuf::from(casadi);
        assert!(casadi.join(wrapper()).exists(), "CASADI_DIR={} holds no {}", casadi.display(), wrapper());
        return (casadi.join("rust"), casadi);
    }
    provide()
}

/// File name of the wrapper library on the target
fn wrapper() -> String {
    let ext = match env::var("CARGO_CFG_TARGET_OS").unwrap().as_str() {
        "windows" => "dll",
        "macos" => "dylib",
        _ => "so",
    };
    format!("libcasadi_rust_wrap.{}", ext)
}

/// Platform tag of the prebuilt archive for the target, if there is one
fn prebuilt_platform() -> Option<&'static str> {
    let os = env::var("CARGO_CFG_TARGET_OS").unwrap();
    let arch = env::var("CARGO_CFG_TARGET_ARCH").unwrap();
    let target_env = env::var("CARGO_CFG_TARGET_ENV").unwrap_or_default();
    match (os.as_str(), arch.as_str(), target_env.as_str()) {
        ("linux", "x86_64", "gnu") => Some("linux64"),
        ("windows", "x86_64", _) => Some("windows64"),
        ("macos", "aarch64", _) => Some("osx-arm64"),
        _ => None,
    }
}

/// Downloads the release archives into OUT_DIR; casadi.rs always comes from the source archive
fn provide() -> (PathBuf, PathBuf) {
    let out = PathBuf::from(env::var("OUT_DIR").unwrap());
    let version = env::var("CARGO_PKG_VERSION").unwrap();
    let source = fetch(&out, &format!("casadi-source-v{}.zip", version), "CASADI_SOURCE_ZIP");
    let gen = out.join("gen");
    fs::create_dir_all(&gen).unwrap();
    let mut archive = open_zip(&source);
    let mut rs = archive.by_name("swig/rust/gen/casadi.rs").expect("no swig/rust/gen/casadi.rs in the source archive");
    io::copy(&mut rs, &mut fs::File::create(gen.join("casadi.rs")).unwrap()).unwrap();
    drop(rs);
    let from_source = env::var_os("CARGO_FEATURE_FROM_SOURCE").is_some();
    let lib = match prebuilt_platform() {
        Some(platform) if !from_source => {
            let prebuilt = fetch(&out, &format!("casadi-{}-{}-rust.zip", version, platform), "CASADI_PREBUILT_ZIP");
            let casadi = out.join("casadi");
            if !casadi.join(wrapper()).exists() {
                open_zip(&prebuilt).extract(&casadi).expect("cannot unpack the prebuilt archive");
            }
            casadi
        }
        _ => build(archive, &out),
    };
    (gen, lib)
}

/// Builds the source archive with cmake into OUT_DIR
fn build(mut archive: zip::ZipArchive<fs::File>, out: &Path) -> PathBuf {
    let src = out.join("source");
    if !src.join("CMakeLists.txt").exists() {
        archive.extract(&src).expect("cannot unpack the source archive");
    }
    let mut c = cmake::Config::new(&src);
    c.profile("Release")
        .define("WITH_RUST", "ON")
        .define("SWIG_IMPORT", "ON")
        .define("WITH_SELFCONTAINED", "ON")
        .define("WITH_BUILD_REQUIRED", "ON")
        // casadi declares Function Send + Sync, as the prebuilt libraries allow
        .define("WITH_THREAD", "ON")
        .define("WITH_THREADSAFE_SYMBOLICS", "ON");
    for p in PLUGINS {
        if env::var_os(format!("CARGO_FEATURE_{}", p.to_uppercase())).is_some() {
            c.define(format!("WITH_{}", p.to_uppercase()), "ON");
            c.define(format!("WITH_BUILD_{}", p.to_uppercase()), "ON");
        }
    }
    if let Ok(args) = env::var("CASADI_CMAKE_ARGS") {
        for a in args.split_whitespace() {
            c.configure_arg(a);
        }
    }
    c.build().join("casadi")
}

/// OUT_DIR/name: the file or URL in env var `over`, else the release asset, checked against checksums.txt
fn fetch(out: &Path, name: &str, over: &str) -> PathBuf {
    let dst = out.join(name);
    if dst.exists() {
        return dst;
    }
    let tmp = out.join(format!("{}.part", name));
    match env::var(over) {
        Ok(src) if !src.starts_with("http") => {
            fs::copy(&src, &tmp).unwrap_or_else(|e| panic!("{}={}: {}", over, src, e));
        }
        Ok(url) => download(&url, &tmp),
        Err(_) => {
            let sum = checksum(name).unwrap_or_else(|| {
                panic!("checksums.txt has no entry for {}: set CASADI_DIR, or {} to a local file or URL", name, over)
            });
            download(&format!("{}/{}/{}", RELEASES, env::var("CARGO_PKG_VERSION").unwrap(), name), &tmp);
            let got = format!("{:x}", Sha256::digest(fs::read(&tmp).unwrap()));
            assert!(got == sum, "{}: SHA256 {} does not match checksums.txt ({})", name, got, sum);
        }
    }
    fs::rename(&tmp, &dst).unwrap();
    dst
}

fn download(url: &str, dst: &Path) {
    let mut resp = ureq::get(url).call().unwrap_or_else(|e| panic!("cannot download {}: {}", url, e));
    let mut f = fs::File::create(dst).unwrap();
    io::copy(&mut resp.body_mut().as_reader(), &mut f).unwrap_or_else(|e| panic!("cannot download {}: {}", url, e));
}

/// Pinned SHA256 of a release asset, from lines `<sha256>  <name>` in checksums.txt
fn checksum(name: &str) -> Option<String> {
    let sums = fs::read_to_string(Path::new(&env::var("CARGO_MANIFEST_DIR").unwrap()).join("checksums.txt")).ok()?;
    sums.lines().find_map(|l| {
        let mut it = l.split_whitespace();
        let sum = it.next()?;
        (it.next()? == name).then(|| sum.to_lowercase())
    })
}

fn open_zip(path: &Path) -> zip::ZipArchive<fs::File> {
    zip::ZipArchive::new(fs::File::open(path).unwrap()).unwrap_or_else(|e| panic!("{}: {}", path.display(), e))
}

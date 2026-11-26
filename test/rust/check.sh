#!/bin/bash
# Rust binding checks against a build tree, as CI runs them.
#
# Usage: check.sh <build dir>|<casadi dir>|<source zip> <prebuilt zip> [example ...]
#   <build dir>  a cmake build with -DWITH_RUST=ON: <build>/swig/rust holds casadi.rs and
#                libcasadi_rust_wrap, <build>/lib libcasadi and the plugins
#   <casadi dir> an unpacked casadi-<os>-rust.zip: the crate in <casadi dir>/rust, no environment needed
#   <source zip> <prebuilt zip>  casadi-source.zip and casadi-<os>-rust.zip as casadi-sys downloads
#                them: a copy of swig/rust, as published; prebuilt zip "-": casadi-sys builds the
#                source zip with cmake (feature from-source), crate tests only
# CASADI_RUST_FEATURES: cargo features of the crate
# Copies the crate tests (test/rust/*.rs) into the crate's tests/ and runs them, then every
# docs/examples/rust example in a scratch directory; modelica_fmu_import needs the OpenModelica
# compiler (omc) and is skipped without it.
set -e -o pipefail
root=$(cd "$(dirname "$0")/../.." && pwd)
features=$CASADI_RUST_FEATURES
if [ -f "$1" ]; then
  export CASADI_SOURCE_ZIP=$(realpath "$1")
  if [ "$2" = - ]; then
    features="from-source $features" examples=no
  else
    export CASADI_PREBUILT_ZIP=$(realpath "$2")
  fi
  shift 2
  # two levels below the repository, where the examples look for docs/ and external_packages/
  tmp=$(mktemp -d -p "$root")
  trap 'rm -rf "$tmp"' EXIT
  crate=$tmp/rust
  cp -rL "$root/swig/rust" "$crate"
  rm -rf "$crate/casadi.rs" "$crate/gen" "$crate/target"
  export CARGO_TARGET_DIR=$crate/target
elif build=$(cd "$1" && pwd) && shift && [ -f "$build/rust/Cargo.toml" ]; then
  export CARGO_TARGET_DIR=$build/rust/target
  crate=$build/rust
else
  export CASADI_RUST_GEN=$build/swig/rust CASADI_PLUGIN_SEARCH_PATH=$build/lib
  export CARGO_TARGET_DIR=$build/swig/rust/cargo
  # plugins built from external projects (uno -> libhighs) find their dependencies here
  export LD_LIBRARY_PATH=$build/lib:$build/external_projects/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}
  crate=$root/swig/rust
fi
mkdir -p "$crate/tests"
cp "$root"/test/rust/*.rs "$crate/tests/"
cd "$crate"
cargo test --features "$features"
[ "$examples" = no ] && exit 0
cargo build --examples --features "$features"

run=$(mktemp -d)
cp "$root/docs/examples/python/rocket.mo" "$run/"
fail=0
# macOS has no timeout(1)
to=$(command -v timeout || command -v gtimeout || true)
[ $# -eq 0 ] && set -- $(cd "$root/docs/examples/rust" && ls *.rs | sed 's/\.rs$//')
for n in "$@"; do
  if [ "$n" = modelica_fmu_import ] && ! command -v omc > /dev/null; then
    echo "skip $n (no omc)"
    continue
  fi
  if (cd "$run" && ${to:+$to 1200} "$CARGO_TARGET_DIR/debug/examples/$n" > "$n.log" 2>&1); then
    echo "ok   $n"
  else
    echo "FAIL $n"
    tail -30 "$run/$n.log"
    fail=1
  fi
done
exit $fail

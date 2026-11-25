#!/bin/bash
# Rust binding checks against a build tree, as CI runs them.
#
# Usage: check.sh <build dir>|<casadi dir> [example ...]
#   <build dir>  a cmake build with -DWITH_RUST=ON: <build>/swig/rust holds casadi.rs and
#                libcasadi_rust_wrap, <build>/lib libcasadi and the plugins
#   <casadi dir> an unpacked casadi-<os>-rust.zip: the crate in <casadi dir>/rust, no environment needed
# Copies the crate tests (test/rust/*.rs) into the crate's tests/ and runs them, then every
# docs/examples/rust example in a scratch directory; modelica_fmu_import needs the OpenModelica
# compiler (omc) and is skipped without it.
set -e -o pipefail
root=$(cd "$(dirname "$0")/../.." && pwd)
build=$(cd "$1" && pwd); shift
if [ -f "$build/rust/Cargo.toml" ]; then
  export CARGO_TARGET_DIR=$build/rust/target
  # Windows: the DLLs are found on PATH
  export PATH=$build:$PATH
  crate=$build/rust
else
  export CASADI_RUST_GEN=$build/swig/rust CASADI_LIB=$build/lib CASADI_PLUGIN_SEARCH_PATH=$build/lib
  export CARGO_TARGET_DIR=$build/swig/rust/cargo
  # plugins built from external projects (uno -> libhighs) find their dependencies here
  export LD_LIBRARY_PATH=$build/lib:$build/external_projects/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}
  crate=$root/swig/rust
fi
mkdir -p "$crate/tests"
cp "$root"/test/rust/*.rs "$crate/tests/"
cd "$crate"
cargo test
cargo build --examples

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

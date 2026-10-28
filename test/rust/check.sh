#!/bin/bash
# Rust binding checks against a build tree, as CI runs them.
#
# Usage: check.sh <build dir>
#   <build dir>  a cmake build with -DWITH_RUST=ON: <build>/swig/rust holds casadi.rs and
#                libcasadi_rust_wrap, <build>/lib libcasadi and the plugins
# Runs the crate tests (swig/rust/tests), then every docs/examples/rust example in a scratch
# directory; modelica_fmu_import needs the OpenModelica compiler (omc) and is skipped without it.
set -e -o pipefail
root=$(cd "$(dirname "$0")/../.." && pwd)
build=$(cd "$1" && pwd)
export CASADI_RUST_GEN=$build/swig/rust CASADI_LIB=$build/lib CASADI_PLUGIN_SEARCH_PATH=$build/lib
export CARGO_TARGET_DIR=$build/swig/rust/cargo
# plugins built from external projects (uno -> libhighs) find their dependencies here
export LD_LIBRARY_PATH=$build/lib:$build/external_projects/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}
cd "$root/swig/rust"
cargo test
cargo build --examples

run=$(mktemp -d)
cp "$root/docs/examples/python/rocket.mo" "$run/"
fail=0
for f in "$root"/docs/examples/rust/*.rs; do
  n=$(basename "$f" .rs)
  if [ "$n" = modelica_fmu_import ] && ! command -v omc > /dev/null; then
    echo "skip $n (no omc)"
    continue
  fi
  if (cd "$run" && timeout 1200 "$CARGO_TARGET_DIR/debug/examples/$n" > "$n.log" 2>&1); then
    echo "ok   $n"
  else
    echo "FAIL $n"
    tail -30 "$run/$n.log"
    fail=1
  fi
done
exit $fail

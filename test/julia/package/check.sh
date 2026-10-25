#!/bin/bash
# Julia binding checks against the CasADiNative package built from a build tree, as CI runs them.
#
# Usage: check.sh <module dir> <env dir> [full]
#   <module dir>  the build's CasADiNative.jl + casadi.jl + libcasadi_wrap (build/swig/julia/CasADiNative)
#   <env dir>     dev environment to (re)build (make_package.jl); afterwards
#                 `julia --project=<env dir> docs/examples/julia/rocket.jl` runs like a user's
#   JULIA         julia executable (default: julia)
#   plus the casadi runtime env (CASADI_PLUGIN_SEARCH_PATH, LD_LIBRARY_PATH).
# Fails on any "Method definition ... overwritten" warning. Always: a script-mode load of the
# layer file (as source builds may include it), the strict package precompile, the package's
# Pkg.test runner (package.jl, runtests.jl, docs/examples/julia) with the plugin path left to
# __init__ (no casadi runtime env), and Aqua; with `full` also the header doctests. Locally, point
# JULIA_DEPOT_PATH at a scratch depot: Aqua is installed from the General registry.
set -e -o pipefail
root=$(cd "$(dirname "$0")/../../.." && pwd)
J=${JULIA:-julia}
mod=$(cd "$1" && pwd)
env=$(mkdir -p "$2" && cd "$2" && pwd)

run() {  # run a command, fail on its exit code or on a method overwrite warning
  "$@" 2>&1 | tee "$env/run.log"
  if grep -q "WARNING: Method definition .* overwritten" "$env/run.log"; then
    echo "check.sh: method overwritten (see above)"; exit 1
  fi
}

echo "== load the layer file as a script"
run "$J" --warn-overwrite=yes -e "include(\"$mod/CasADiNative.jl\")"
echo "== build + precompile the package"
run "$J" --warn-overwrite=yes "$root/test/julia/package/make_package.jl" "$mod" "$env" Aqua
export JULIA_PROJECT=$env
echo "== package runtests.jl (binaries found via __init__ and RUNPATH, as from the JLLs)"
run env -u CASADI_PLUGIN_SEARCH_PATH -u LD_LIBRARY_PATH "$J" --warn-overwrite=yes "$root/test/julia/package/runtests.jl"
echo "== Aqua"
run "$J" "$root/test/julia/package/aqua.jl"
[ "$3" = full ] || exit 0
echo "== header doctests"
cd "$root" && JULIA=$J python3 docs/api/extra/doctest_headers.py --lang julia

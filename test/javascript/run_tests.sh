#!/usr/bin/env bash
# run_tests.sh -- top-level harness for the wasm+TS bindings.
#
# Two layers:
#   (1) TypeScript declaration-file typechecking via `tsc --noEmit`
#       against test_typecheck.ts + casadi.d.ts.  This is the main
#       guardrail when extending typemap macros / casadi.i type
#       annotations -- a regression in `xTsStub` plumbing surfaces
#       here as a tsc error.
#   (2) Node runtime smoke check via test_runtime.cjs.  Loads the
#       generated module, asserts module shape + a small symbolic
#       eval (SX.sym -> sin -> jacobian -> Func).
#
# Both layers depend on a populated build-wasm/swig/wasm-js
# directory; the script bails with SKIP if .d.ts / .js are missing.

set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WASM_DIR="$(cd "$HERE/../../build-wasm/swig/wasm-js" 2>/dev/null && pwd || true)"

if [[ -z "$WASM_DIR" || ! -f "$WASM_DIR/casadi.d.ts" || ! -f "$WASM_DIR/casadi.js" ]]; then
  echo "SKIP: build-wasm output missing.  Build it first:"
  echo "  (cd $(realpath "$HERE/../../build-wasm") && ./run.sh build)"
  exit 2
fi

echo "==> casadi.d.ts: $(wc -l < "$WASM_DIR/casadi.d.ts") lines"
echo "==> tsc typecheck (test_typecheck.ts)"
cd "$HERE"
tsc --noEmit -p tsconfig.json
echo "ok -- typecheck passed"
echo

echo "==> node runtime test (test_runtime.cjs)"
node "$HERE/test_runtime.cjs"
echo

echo "==> node mx.js (ports test/python/mx.py)"
node "$HERE/mx.js"
echo

echo "==> node sparsity.js (ports test/python/sparsity.py subset)"
node "$HERE/sparsity.js"
echo

echo "==> node matrix.js (ports test/python/matrix.py subset)"
node "$HERE/matrix.js"
echo

echo "==> node multiplication.js (ports test/python/multiplication.py subset)"
node "$HERE/multiplication.js"
echo

echo "==> node ad.js (ports test/python/ad.py subset)"
node "$HERE/ad.js"
echo

# Phase 0.2 audit: dispatch divergence cases.  Phase 3.1-3.4 land the
# wasm-side `swig_can_<X>` probe wrappers + fuzzy coercion + arity
# dispatch for default args; all 6 cases now pass.
echo "==> node phase0_dispatch.js"
node "$HERE/phase0_dispatch.js"
echo

# Phase 4.3: director / Callback subclass support (full end-to-end).
echo "==> node director.js (Phase 4.3 directors)"
node "$HERE/director.js"
echo

# Conic / qpsol across all loaded plugins (HiGHS, qrqp, ipqp, qpoases,
# daqp, nlpsol).  Subset port of test/python/conic.py.
echo "==> node conic.js (ports test/python/conic.py subset, multi-plugin)"
node "$HERE/conic.js"
echo

# Linsol across all loaded plugins (qr, ldl, lsqr, symbolicqr, csparse,
# lapacklu, lapackqr).  Subset port of test/python/linearsolver.py.
echo "==> node linsol.js (ports test/python/linearsolver.py subset)"
node "$HERE/linsol.js"
echo

# Rootfinder across all loaded plugins (newton, fast_newton, nlpsol,
# kinsol).  Subset port of test/python/implicitfunction.py.
echo "==> node rootfinder.js (ports test/python/implicitfunction.py subset)"
node "$HERE/rootfinder.js"
echo

# SuperSCS is not built in the wasm job (WITH_SUPERSCS defaults OFF; its
# sources still need the F77 hidden CHARACTER-length patch the vendored
# qpOASES got).  Excluded from the suite until patched -- see
# superscs_smoke.js for the intended SOCP smoke.

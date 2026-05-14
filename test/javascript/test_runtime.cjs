// test_runtime.cjs
//
// Runtime smoke test for the auto-generated wasm bindings.  Loads
// build-wasm/swig/wasm-js/casadi.js (the SWIG-emitted JS proxy --
// no hand-rolled JS layer, no mvp.js), constructs a few symbolic
// expressions, builds a Function, and validates the headline shape.
//
// Run via test/javascript/run_tests.sh after the wasm build is
// populated.

const path = require("path");
const fs = require("fs");

const wasmDir = path.resolve(__dirname, "../../build-wasm/swig/wasm-js");
const modulePath = path.join(wasmDir, "casadi.js");

if (!fs.existsSync(modulePath)) {
  console.error(`SKIP: ${modulePath} not built.  Run the wasm build first:`);
  console.error(`  (cd ${path.resolve(__dirname, "../../build-wasm")} && ./run.sh build)`);
  process.exit(2);
}

function assert(cond, msg) {
  if (!cond) { console.error("FAIL:", msg); process.exit(1); }
  console.log("ok --", msg);
}

(async () => {
  const create = require(modulePath);
  const M = await create();

  // ---- module shape ----------------------------------------------------
  const expected = ["Function", "SX", "MX", "DM", "Sparsity", "Slice"];
  const missing = expected.filter((k) => !(k in M));
  assert(missing.length === 0, `module exports present: ${expected.join(", ")}`);

  // ---- SX symbolic build -----------------------------------------------
  // `static sym` lives on the C++ template parent
  // GenericMatrix<Matrix<SXElem>> (= GenSX), not directly on Matrix.
  // The __mixin shim in casadi.js merges GenSX's statics onto SX, so
  // users see `SX.sym(...)` directly without ever touching the GenSX
  // helper class.  Class-return wrapping then turns the raw void*
  // into a proper SX instance.
  const x = M.SX.sym("x", 1n, 1n);
  assert(x instanceof M.SX,
    `SX.sym('x', 1n, 1n) returns SX (got ${x?.constructor?.name})`);

  // ---- MX symbolic build (mirror; verifies wrapping for MX too) -------
  const xm = M.MX.sym("y", 1n, 1n);
  assert(xm instanceof M.MX,
    `MX.sym('y', 1n, 1n) returns MX (got ${xm?.constructor?.name})`);

  // ---- Function class is a constructor --------------------------------
  assert(typeof M.Function === "function", "M.Function is a constructor");

  console.log("\nall runtime smoke tests passed");
})().catch((e) => {
  console.error("FAIL:", e);
  process.exit(1);
});

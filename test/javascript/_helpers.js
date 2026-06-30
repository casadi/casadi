// _helpers.js -- shared assertion / module-load helpers for the JS test suite.
// Used across the per-module ports of test/python/*.py.

const path = require("path");
const fs   = require("fs");

const wasmDir    = path.resolve(__dirname, "../../build-wasm/swig/wasm-js");
const modulePath = path.join(wasmDir, "casadi.js");

function checkModuleOrSkip() {
  if (!fs.existsSync(modulePath)) {
    console.error(`SKIP: ${modulePath} not built.`);
    process.exit(2);
  }
}

function assertEqual(actual, expected, msg) {
  const same = (actual === expected) ||
    ((typeof actual === "bigint" || typeof actual === "number") &&
     (typeof expected === "bigint" || typeof expected === "number") &&
     BigInt(actual) === BigInt(expected));
  if (!same) throw new Error(`${msg}: expected ${expected}, got ${actual}`);
}

function assertTrue(cond, msg) {
  if (!cond) throw new Error(`${msg}: expected truthy, got ${cond}`);
}

function assertAlmost(actual, expected, places, msg) {
  if (Math.abs(actual - expected) > Math.pow(10, -places)) {
    throw new Error(`${msg}: expected ${expected}, got ${actual}`);
  }
}

function assertArrayAlmostEqual(actual, expected, places, msg) {
  if (actual.length !== expected.length) {
    throw new Error(`${msg}: length expected ${expected.length}, got ${actual.length}`);
  }
  const tol = Math.pow(10, -places);
  for (let i = 0; i < actual.length; ++i) {
    if (Math.abs(actual[i] - expected[i]) > tol) {
      throw new Error(`${msg}[${i}]: expected ${expected[i]}, got ${actual[i]}`);
    }
  }
}

function checkArray(actual, expected, places, msg) {
  const exp = Array.isArray(expected) ? expected
            : (typeof expected === "number") ? [expected]
            : expected.nonzeros();
  const act = (actual && actual.nonzeros) ? actual.nonzeros() : actual;
  assertArrayAlmostEqual(act, exp, places, msg);
}

// JS port of the simplest checkfunction variants from
// test/python/helpers.py.  Used heavily by ad.py/function.py tests to
// confirm two `Function` objects (one SX-built, one MX-built, or one
// with `jit=true` vs without, ...) produce numerically identical output
// on a list of inputs.
//
// Strategy mirrors checkfunction_light: call each Function with the
// same inputs, compare outputs element-wise, throw on mismatch.  We
// stop short of porting the full forward/reverse/jacobian/hessian
// equivalence checks -- those rely on `forward(...)`/`reverse(...)`
// global builders that aren't exposed yet in the wasm-js Function
// dispatcher, and the simpler form already unlocks the bulk of the
// ad-equivalence tests.
function checkfunction_light(M, trial, solution, inputs, digits) {
  if (digits === undefined) digits = 9;
  if (!inputs) inputs = [];
  // Sanity-check arities up-front.
  assertEqual(trial.n_in(),  solution.n_in(),  "checkfunction: n_in differs");
  assertEqual(trial.n_out(), solution.n_out(), "checkfunction: n_out differs");
  // Repeated evaluation (Python helpers does 2 passes too).
  for (let pass = 0; pass < 2; ++pass) {
    const t_out = trial.call(inputs);
    const s_out = solution.call(inputs);
    if (t_out.length !== s_out.length) {
      throw new Error(`checkfunction: output count differs (${t_out.length} vs ${s_out.length})`);
    }
    for (let i = 0; i < t_out.length; ++i) {
      const a = t_out[i].nonzeros();
      const b = s_out[i].nonzeros();
      assertArrayAlmostEqual(a, b, digits,
        `checkfunction[pass=${pass}] output ${i}`);
    }
  }
}

// Build a Function with MX symbols mirroring a Function built with SX.
// Used by ad-equivalence patterns: build the same expression once in
// SX and once in MX, then checkfunction_light(M, mx_fn, sx_fn, inputs).
function checkfunction(M, trial, solution, inputs, opts) {
  // For now alias to the light variant; full checkfunction (with
  // forward/reverse/jacobian/hessian sweeps) requires
  // M.forward/M.reverse which aren't dispatcher-stable yet.
  return checkfunction_light(M, trial, solution, inputs,
                             opts && opts.digits);
}

function assertThrows(fn, msg) {
  let threw = false;
  try { fn(); } catch (e) { threw = true; }
  if (!threw) throw new Error(`${msg}: expected to throw`);
}

// Runner: takes a tests array of [name, fn] and runs them against the
// loaded module M.
async function runTests(tests) {
  checkModuleOrSkip();
  const create = require(modulePath);
  const M = await create();
  let pass = 0, fail = 0;
  const failures = [];
  for (const [name, fn] of tests) {
    try { fn(M); console.log(`ok   -- ${name}`); pass++; }
    catch (e) { console.log(`FAIL -- ${name}: ${e.message}`); fail++; failures.push([name, e.message]); }
  }
  console.log(`\n${pass}/${pass + fail} passed`);
  if (failures.length > 0) {
    console.log("\n--- failures ---");
    for (const [n, m] of failures) console.log(`  ${n}: ${m}`);
  }
  process.exit(fail === 0 ? 0 : 1);
}

module.exports = {
  modulePath,
  assertEqual,
  assertTrue,
  assertAlmost,
  assertArrayAlmostEqual,
  checkArray,
  assertThrows,
  checkfunction,
  checkfunction_light,
  runTests,
};

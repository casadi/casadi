// phase0_dispatch.js -- failing tests that expose `__can` vs `to_ptr`
// divergence (Phase 0.2 audit; expected to pass after Phase 3.1 lands).
//
// These cases all route to the wrong constructor in the current
// `constructor.name`-based dispatcher. See PHASE0_AUDIT.md for analysis.
//
// Until Phase 3.1: `node phase0_dispatch.js` should report 0/N passed and
// exit non-zero.  The harness invokes us with `EXPECT_FAIL=1` so test
// runs return 0 when these tests fail as expected.

const path = require("path");
const fs = require("fs");

const wasmDir = path.resolve(__dirname, "../../build-wasm/swig/wasm-js");
const modulePath = path.join(wasmDir, "casadi.js");
if (!fs.existsSync(modulePath)) {
  console.error(`SKIP: ${modulePath} not built.`);
  process.exit(2);
}

function assertThrows(fn, msg) {
  let threw = false;
  try { fn(); } catch (_) { threw = true; }
  if (!threw) throw new Error(`${msg}: expected throw, got success`);
}

function assertAlmostEqual(actual, expected, places, msg) {
  const tol = Math.pow(10, -places);
  if (Math.abs(actual - expected) > tol) {
    throw new Error(`${msg}: expected ${expected}, got ${actual} (|diff|>${tol})`);
  }
}

// ----- Test cases (each currently FAILS; Phase 3.1 should make them pass) -----

// 1. MX(string) — no ctor accepts string; should TypeError, currently
//    silently routes to MX(double) → garbage.
function test_mx_from_string_should_throw(M) {
  assertThrows(() => { new M.MX("hello"); },
    "MX(string) must reject");
}

// 2. MX(null) — no ctor accepts null; should TypeError, currently
//    silently routes to MX(double) → MX(0.0).
function test_mx_from_null_should_throw(M) {
  assertThrows(() => { new M.MX(null); },
    "MX(null) must reject");
}

// 3. MX(other_MX) — the copy ctor.  Should produce an MX equal to the input.
//    Currently falls through to MX(double) with the proxy object as the double.
function test_mx_copy_ctor(M) {
  const a = M.MX.sym("a");
  const b = new M.MX(a);             // expected: copy
  // The copy should have the same nonzero count and same symbolic structure.
  // (Don't compare _ptr identities — they're separate C++ instances.)
  if (b.is_constant()) throw new Error("MX(MX) copy lost symbolic-ness");
}

// 4. MX([1.0, 2.0]) — fuzzy coercion via DM.  matlab/python accept this and
//    construct an MX wrapping a DM column.  Currently the array falls through
//    to MX(double).
function test_mx_from_array_coerces(M) {
  const a = new M.MX([1.0, 2.0, 3.0]);
  if (Number(a.numel()) !== 3) {
    throw new Error(`MX([1,2,3]) numel: expected 3, got ${a.numel()}`);
  }
}

// 5. vertcat with mixed numbers / DMs — vertcat overloads include
//    `vertcat(vector<DM>)`.  Passing JS array of numbers should coerce.
function test_vertcat_numbers_coerce(M) {
  // vertcat(array) is the simplest fuzzy case: each entry must coerce to DM.
  const v = M.vcat([1.0, 2.0, 3.0]);
  if (Number(v.numel()) !== 3) {
    throw new Error(`vertcat([1,2,3]) numel: expected 3, got ${v.numel()}`);
  }
}

// 6. null first-arg dispatch — declaration-order tie-break is wrong when
//    multiple overloads have arity ≥ 1 and first parm is class-typed.
//    With null, __can returns true for ALL of them; currently picks the
//    FIRST registered.  matlab/python: rejects with "no overload matches".
function test_null_first_arg_rejects(M) {
  // vertcat(null) shouldn't silently route to vertcat<MX> or whatever
  // declared first.
  assertThrows(() => { M.vertcat(null); }, "vertcat(null) must reject");
}

// ----- Driver -----

const tests = [
  ["test_mx_from_string_should_throw",   test_mx_from_string_should_throw],
  ["test_mx_from_null_should_throw",     test_mx_from_null_should_throw],
  ["test_mx_copy_ctor",                  test_mx_copy_ctor],
  ["test_mx_from_array_coerces",         test_mx_from_array_coerces],
  ["test_vertcat_numbers_coerce",        test_vertcat_numbers_coerce],
  ["test_null_first_arg_rejects",        test_null_first_arg_rejects],
];

(async () => {
  const create = require(modulePath);
  const M = await create();
  let pass = 0, fail = 0;
  for (const [name, fn] of tests) {
    try {
      await fn(M);
      console.log(`ok   -- ${name}`);
      pass++;
    } catch (e) {
      console.log(`FAIL -- ${name}: ${e.message}`);
      fail++;
    }
  }
  console.log(`\nphase0_dispatch: ${pass}/${pass + fail} passed`);
  // EXPECT_FAIL=1 inverts the exit code (harness uses this until Phase 3.1).
  const expectFail = process.env.EXPECT_FAIL === "1";
  if (expectFail) {
    process.exit(fail > 0 ? 0 : 1);  // expect failures; success if some failed
  } else {
    process.exit(fail === 0 ? 0 : 1);
  }
})().catch((e) => { console.error("FATAL:", e); process.exit(1); });

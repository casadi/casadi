// director.js -- Phase 4.3 director / Callback subclass test.

const path = require("path");
const fs = require("fs");

const wasmDir = path.resolve(__dirname, "../../build-wasm/swig/wasm-js");
const modulePath = path.join(wasmDir, "casadi.js");
if (!fs.existsSync(modulePath)) {
  console.error(`SKIP: ${modulePath} not built.`);
  process.exit(2);
}

function assertEqual(actual, expected, msg) {
  if (actual !== expected) throw new Error(`${msg}: expected ${expected}, got ${actual}`);
}
function assertAlmost(actual, expected, places, msg) {
  if (Math.abs(actual - expected) > Math.pow(10, -places)) {
    throw new Error(`${msg}: expected ${expected}, got ${actual}`);
  }
}

// ----- tests -----

function test_director_construct(M) {
  // Subclass without overriding methods.  Constructor must complete
  // and produce a valid _ptr.
  class PlainCB extends M.Callback {
    constructor() { super(); }
  }
  const cb = new PlainCB();
  if (typeof cb._ptr !== "number") throw new Error("missing numeric _ptr");
  if (cb._ptr === 0) throw new Error("director ctor returned null _ptr");
}

function test_director_is_null(M) {
  // A freshly-constructed director (before construct() runs) must
  // report is_null=true.  This used to return false due to the
  // C++ multi-inheritance offset bug; Phase 4.3's cast-chain converter
  // in SWIG_WASMJS_ConvertPtr fixes it.
  class PlainCB extends M.Callback {
    constructor() { super(); }
  }
  const cb = new PlainCB();
  if (!cb.is_null()) throw new Error("is_null should be true before construct()");
}

function test_director_eval_override(M) {
  // JS subclass overrides eval to double the input.  C++ Function::call
  // invokes virtual eval on the SwigDirector_Callback, which dispatches
  // to the JS subclass.
  class DoubleCB extends M.Callback {
    constructor() {
      super();
      this.construct("DoubleCB", null);
    }
    get_n_in()  { return 1n; }
    get_n_out() { return 1n; }
    eval(args) {
      // C++ Callback::eval takes ONE arg: std::vector<DM>.  Our
      // director marshaling wraps that as a JS array of DM proxies
      // and pushes it onto the JS args array.  So:
      //   args[0]    = the JS array of DM proxies (the vector)
      //   args[0][0] = the first DM
      const x = args[0][0].nonzeros()[0];
      return [M.DM(2 * x)];
    }
  }
  const cb = new DoubleCB();
  const out = cb.call([M.DM(3)]);
  assertAlmost(out.nonzeros()[0], 6, 8, "eval(3) should be 6");
}

// ----- driver -----

const tests = [
  ["test_director_construct",     test_director_construct],
  ["test_director_is_null",       test_director_is_null],
  ["test_director_eval_override", test_director_eval_override],
];

(async () => {
  const create = require(modulePath);
  const M = await create();
  let pass = 0, fail = 0;
  for (const [name, fn] of tests) {
    try { fn(M); console.log(`ok -- ${name}`); pass++; }
    catch (e) { console.log(`FAIL ${name}: ${e.message}\n${e.stack}`); fail++; }
  }
  console.log(`\n${pass}/${pass + fail} passed`);
  process.exit(fail === 0 ? 0 : 1);
})().catch((e) => { console.error("FATAL:", e); process.exit(1); });

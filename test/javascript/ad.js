// ad.js -- JS port of test/python/ad.py (subset).
// Covers forward/reverse AD primitives: jacobian, gradient, hessian
// on SX/MX, evaluated via Function::call against numeric inputs.

const path = require("path");
const fs = require("fs");

const wasmDir = path.resolve(__dirname, "../../build-wasm/swig/wasm-js");
const modulePath = path.join(wasmDir, "casadi.js");
if (!fs.existsSync(modulePath)) {
  console.error(`SKIP: ${modulePath} not built.`);
  process.exit(2);
}

function assertEqual(actual, expected, msg) {
  const same = (actual === expected) ||
    ((typeof actual === 'bigint' || typeof actual === 'number') &&
     (typeof expected === 'bigint' || typeof expected === 'number') &&
     BigInt(actual) === BigInt(expected));
  if (!same) throw new Error(`${msg}: expected ${expected}, got ${actual}`);
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

// ----- tests -----

function test_sx_jacobian_scalar(M) {
  // f(x) = sin(x), df/dx = cos(x). At x=1: cos(1).
  const x = M.SX.sym("x", 1n, 1n);
  const f_sym = M.sin(x);
  const df = M.jacobian(f_sym, x);
  const f = M.Function("df", [x], [df]);
  const out = f(M.DM(1));
  assertArrayAlmostEqual(out.nonzeros(), [Math.cos(1)], 8, "d(sin)/dx at 1");
}

function test_mx_jacobian_scalar(M) {
  const x = M.MX.sym("x", 1n, 1n);
  const f_sym = M.times(x, x);    // f(x) = x^2
  const df = M.jacobian(f_sym, x);
  const f = M.Function("df", [x], [df]);
  // df/dx = 2x.  At x=3, expect 6.
  const out = f(M.DM(3));
  assertArrayAlmostEqual(out.nonzeros(), [6], 8, "d(x^2)/dx at 3");
}

function test_sx_gradient_vector(M) {
  // f(x, y) = x^2 + 3*x*y + y^2.  grad = [2x + 3y, 3x + 2y].
  const x = M.SX.sym("x");
  const y = M.SX.sym("y");
  const f_sym = M.plus(
      M.plus(M.times(x, x), M.times(M.times(M.SX(3), x), y)),
      M.times(y, y));
  const gx = M.gradient(f_sym, x);
  const gy = M.gradient(f_sym, y);
  const f = M.Function("g", [x, y], [gx, gy]);
  // At (x=2, y=5): gx = 2*2 + 3*5 = 19; gy = 3*2 + 2*5 = 16.
  const outs = f(M.DM(2), M.DM(5));
  assertArrayAlmostEqual(outs[0].nonzeros(), [19], 8, "gx at (2,5)");
  assertArrayAlmostEqual(outs[1].nonzeros(), [16], 8, "gy at (2,5)");
}

function test_mx_jacobian_vector(M) {
  // f(x) = [x0*x1, x0+x1].  J = [[x1, x0], [1, 1]].
  // Build inputs as a 2-vector MX, output as a 2-vector MX.
  const x = M.MX.sym("x", 2n, 1n);
  const x0 = x.get(false, M.Slice(0n, false));
  const x1 = x.get(false, M.Slice(1n, false));
  const f_sym = M.vcat([M.times(x0, x1), M.plus(x0, x1)]);
  const J = M.jacobian(f_sym, x);
  assertEqual(J.size1(), 2n, "J size1");
  assertEqual(J.size2(), 2n, "J size2");

  const f = M.Function("J", [x], [J]);
  // At x = [3, 5]: J = [[5, 3], [1, 1]].
  // Column-major: [5, 1, 3, 1].
  const xv = M.DM([3, 5]);
  const out = f(xv);
  assertArrayAlmostEqual(out.nonzeros(), [5, 1, 3, 1], 8, "J at (3,5)");
}

function test_sx_jacobian_of_gradient(M) {
  // 2nd-derivative via jacobian-of-gradient (avoids hessian's
  // [H, g] dual output, which currently has marshalling issues).
  // f(x) = x^4.  d2f/dx2 = 12*x^2.  At x=2: 48.
  const x = M.SX.sym("x");
  const f_sym = M.power(x, M.SX(4));
  const g = M.gradient(f_sym, x);
  const H = M.jacobian(g, x);
  const f = M.Function("H", [x], [H]);
  const out = f(M.DM(2));
  assertArrayAlmostEqual(out.nonzeros(), [48], 8, "d2(x^4)/dx2 at 2");
}

// ----- SX-vs-MX equivalence: the basic checkfunction_light pattern.
// Builds the same expression in both SX and MX, then verifies they
// produce identical numerical output on the same DM input.  This is
// the bulk of the ad.py "test_*eval_*" pattern adapted to wasm-js.

const h = require("./_helpers");

function _build_pair(M, builder, sym_in, n) {
  // Build the same Function in SX (with M.SX.sym) and MX (with M.MX.sym),
  // each of vector size `n`.  `builder(M, x)` returns the output expression.
  const xs = M.SX.sym("x", BigInt(n));
  const xm = M.MX.sym("x", BigInt(n));
  return [
    M.Function("fs", [xs], [builder(M, xs)]),
    M.Function("fm", [xm], [builder(M, xm)]),
  ];
}

function test_ad_equiv_x_squared(M) {
  // f(x) = x^2 elementwise, on a length-3 input.
  const [fs, fm] = _build_pair(M, (M, x) => M.times(x, x), null, 3);
  h.checkfunction_light(M, fs, fm, [M.DM([1.5, -2, 0.5])], 10);
}

function test_ad_equiv_trig_sum(M) {
  // f(x) = sin(x) + cos(x).
  const [fs, fm] = _build_pair(M, (M, x) => M.plus(M.sin(x), M.cos(x)), null, 4);
  h.checkfunction_light(M, fs, fm, [M.DM([0, 0.7, 1.4, 2.1])], 10);
}

function test_ad_equiv_dot_product(M) {
  // f(x) = dot(x, x).  SX vs MX should agree numerically.
  const builder = (M, x) => M.dot(x, x);
  const [fs, fm] = _build_pair(M, builder, null, 5);
  h.checkfunction_light(M, fs, fm, [M.DM([1, 2, 3, 4, 5])], 10);
}

function test_ad_equiv_jacobian(M) {
  // jacobian(sin(x), x) should agree as a Function between SX and MX.
  const xs = M.SX.sym("x", 3n);
  const fs = M.Function("Js", [xs], [M.jacobian(M.sin(xs), xs)]);
  const xm = M.MX.sym("x", 3n);
  const fm = M.Function("Jm", [xm], [M.jacobian(M.sin(xm), xm)]);
  h.checkfunction_light(M, fs, fm, [M.DM([0.1, 0.4, 0.9])], 10);
}

// ----- driver -----

const tests = [
  ["test_sx_jacobian_scalar",  test_sx_jacobian_scalar],
  ["test_mx_jacobian_scalar",  test_mx_jacobian_scalar],
  ["test_sx_gradient_vector",  test_sx_gradient_vector],
  ["test_mx_jacobian_vector",  test_mx_jacobian_vector],
  ["test_sx_jacobian_of_gradient", test_sx_jacobian_of_gradient],
  ["test_ad_equiv_x_squared",  test_ad_equiv_x_squared],
  ["test_ad_equiv_trig_sum",   test_ad_equiv_trig_sum],
  ["test_ad_equiv_dot_product", test_ad_equiv_dot_product],
  ["test_ad_equiv_jacobian",   test_ad_equiv_jacobian],
];

(async () => {
  const create = require(modulePath);
  const M = await create();
  let pass = 0, fail = 0;
  for (const [name, fn] of tests) {
    try { fn(M); console.log(`ok -- ${name}`); pass++; }
    catch (e) { console.log(`FAIL ${name}: ${e.message}`); fail++; }
  }
  console.log(`\n${pass}/${pass + fail} passed`);
  process.exit(fail === 0 ? 0 : 1);
})().catch((e) => { console.error("FATAL:", e); process.exit(1); });

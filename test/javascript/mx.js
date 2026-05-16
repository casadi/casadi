// mx.js -- JS port of test/python/mx.py.
//
// API conventions:
//   - Vectors cross the wasm boundary as plain JS arrays in both
//     directions (XVector classes stay an internal SWIG detail).
//   - Function returns auto-wrap class types as proxies; nonzeros()
//     returns a plain array of numbers.
//   - Raw JS arrays auto-promote to DM via the casadi to_ptr<DM>
//     infra (same as Python's `f([3, 7])`): `[3, 7]` becomes a 2x1
//     column DM, `[[1, 2], [3, 4]]` becomes a 2x2 DM (row-major in
//     JS, column-major in DM).  So `f([2.3])` and `f([[3, 7]])` work
//     directly -- no explicit M.DM wrap.  Uses the callable form for
//     the natural n_out=1 unwrap; use `f.call([...])` for stable list
//     return.
//
// Each test_* function is `(M) => void` that throws on failure.

const path = require("path");
const fs = require("fs");

const wasmDir = path.resolve(__dirname, "../../build-wasm/swig/wasm-js");
const modulePath = path.join(wasmDir, "casadi.js");
if (!fs.existsSync(modulePath)) {
  console.error(`SKIP: ${modulePath} not built.`);
  process.exit(2);
}

function assertEqual(actual, expected, msg) {
  // Allow cross-type number/bigint equality (size_t crosses as Number,
  // casadi_int as BigInt).
  const same = (actual === expected) ||
    ((typeof actual === 'bigint' || typeof actual === 'number') &&
     (typeof expected === 'bigint' || typeof expected === 'number') &&
     BigInt(actual) === BigInt(expected));
  if (!same) throw new Error(`${msg}: expected ${expected}, got ${actual}`);
}

function assertAlmostEqual(actual, expected, places, msg) {
  const tol = Math.pow(10, -places);
  if (Math.abs(actual - expected) > tol) {
    throw new Error(`${msg}: expected ${expected}, got ${actual} (|diff|>${tol})`);
  }
}

// Compare a DM's nonzeros() (a JS number array) to an expected array.
function dmArrayEqual(dm, expected, places, msg) {
  const nz = dm.nonzeros();
  if (nz.length !== expected.length) {
    throw new Error(`${msg}: nnz expected ${expected.length}, got ${nz.length}`);
  }
  const tol = Math.pow(10, -places);
  for (let i = 0; i < nz.length; ++i) {
    if (Math.abs(nz[i] - expected[i]) > tol) {
      throw new Error(`${msg}[${i}]: expected ${expected[i]}, got ${nz[i]}`);
    }
  }
}


// ============================================================================
//  tests
// ============================================================================

function test_MX1(M) {
  const x = M.MX.sym("x", 2n, 3n);
  assertEqual(x.size1(), 2n, "MX.size1");
  assertEqual(x.size2(), 3n, "MX.size2");
}

function test_MXvertcat(M) {
  const x = M.MX.sym("x", 1n, 3n);
  const y = M.MX.sym("y", 1n, 3n);
  const z = M.vcat([x, y]);
  assertEqual(z.size1(), 2n, "vertcat.size1");
  assertEqual(z.size2(), 3n, "vertcat.size2");
}

function test_MX_fun1(M) {
  const x = M.MX.sym("x", 1, 1);
  const y = M.times(M.MX(2), x);
  const f = M.Function("f", [x], [y]);
  assertEqual(f.n_in(),  1n, "f.n_in");
  assertEqual(f.n_out(), 1n, "f.n_out");
  const out = f(3);
  assertAlmostEqual(out.nonzeros()[0], 6, 10, "f(3) == 6");
}

function test_MXfunction2(M) {
  const x = M.MX.sym("x", 1n, 1n);
  const y = M.MX.sym("y", 1n, 1n);
  const f = M.Function("f", [x, y],
    [M.plus(x, y), M.times(y, x)]);
  assertEqual(f.n_in(),  2n, "f.n_in");
  assertEqual(f.n_out(), 2n, "f.n_out");
  const outs = f(3, 7);
  assertAlmostEqual(outs[0].nonzeros()[0], 10, 10, "x+y");
  assertAlmostEqual(outs[1].nonzeros()[0], 21, 10, "y*x");
}

function test_MXfunction3(M) {
  const xy = M.MX.sym("xy", 2n, 1n);
  const xy0 = xy.get(false, M.Slice(0n, false));
  const xy1 = xy.get(false, M.Slice(1n, false));
  const f = M.Function("f", [xy],
    [M.plus(xy0, xy1), M.times(xy0, xy1)]);
  const outs = f([3, 7]);
  assertAlmostEqual(outs[0].nonzeros()[0], 10, 10, "xy[0]+xy[1]");
  assertAlmostEqual(outs[1].nonzeros()[0], 21, 10, "xy[0]*xy[1]");
}

function test_MXfunction3b(M) {
  const xy = M.MX.sym("xy", 1n, 2n);
  const xy00 = xy.get(false, M.Slice(0n, false));
  const xy01 = xy.get(false, M.Slice(1n, false));
  const f = M.Function("f", [xy],
    [M.plus(xy00, xy01), M.times(xy00, xy01)]);
  // Nested array [[3, 7]] is row-major in JS: a 1x2 row DM.
  const outs = f([[3, 7]]);
  assertAlmostEqual(outs[0].nonzeros()[0], 10, 10, "[0,0]+[0,1]");
  assertAlmostEqual(outs[1].nonzeros()[0], 21, 10, "[0,0]*[0,1]");
}

function test_MXfunction4(M) {
  const xy = M.MX.sym("xy", 2n, 1n);
  const xy0 = xy.get(false, M.Slice(0n, false));
  const xy1 = xy.get(false, M.Slice(1n, false));
  const z = M.vcat([M.plus(xy0, xy1), M.times(xy0, xy1)]);
  const f = M.Function("f", [xy], [z]);
  const out = f([3, 7]);
  const out0 = out;
  assertEqual(out0.size1(), 2n, "out.size1");
  assertEqual(out0.size2(), 1n, "out.size2");
  const nz = out0.nonzeros();
  assertAlmostEqual(nz[0], 10, 10, "z[0]");
  assertAlmostEqual(nz[1], 21, 10, "z[1]");
}

function test_MXfunction5(M) {
  const xy = M.MX.sym("xy", 2n, 1n);
  const xy0 = xy.get(false, M.Slice(0n, false));
  const xy1 = xy.get(false, M.Slice(1n, false));
  const z = M.hcat([M.plus(xy0, xy1), M.times(xy0, xy1)]);
  const f = M.Function("f", [xy], [z]);
  const out = f([3, 7]);
  const out0 = out;
  assertEqual(out0.size1(), 1n, "out.size1");
  assertEqual(out0.size2(), 2n, "out.size2");
}

function test_identitySX(M) {
  const x = M.SX.sym("x", 1n, 1n);
  const f = M.Function("f", [x], [x]);
  const out = f(3);
  assertAlmostEqual(out.nonzeros()[0], 3, 10, "identity SX");
}

function test_identityMX(M) {
  const x = M.MX.sym("x", 1n, 1n);
  const f = M.Function("f", [x], [x]);
  const out = f(3);
  assertAlmostEqual(out.nonzeros()[0], 3, 10, "identity MX");
}

function test_MXtrans(M) {
  const x = M.MX.sym("x", 2n, 3n);
  const z = x.T;
  assertEqual(z.size1(), 3n, "z.size1");
  assertEqual(z.size2(), 2n, "z.size2");
}

function test_issue107(M) {
  const x = M.MX.sym("x", 1n, 1n);
  const y = M.MX.sym("y", 1n, 1n);
  let z = x;
  z = M.plus(z, y);
  assertEqual(x.is_symbolic(), true,  "x symbolic");
  assertEqual(z.is_symbolic(), false, "z not symbolic");
}

function test_MXd_trivial(M) {
  const X = M.MX.sym("X", 10n, 1n);
  const J = M.jacobian(X, X);
  assertEqual(J.nnz(), 10, "J.nnz");
  assertEqual(J.size1(), 10n, "J.size1");
  assertEqual(J.size2(), 10n, "J.size2");
}

function test_issue184(M) {
  const x = M.MX.sym("x", 3n, 1n);
  const y = x.get(false, M.Slice(0n, 0n, 1n));
  assertEqual(y.nnz(), 0, "empty slice nnz");
}

function test_jacobian_tools(M) {
  const X = M.MX.sym("X", 1n, 1n);
  const Y = M.jacobian(M.power(X, M.MX(2)), X, null);
  const f = M.Function("f", [X], [Y]);
  const out = f(2.3);
  assertAlmostEqual(out.nonzeros()[0], 4.6, 6, "d(X^2)/dX at 2.3");
}

function test_if_else(M) {
  const x = M.MX.sym("x", 1n, 1n);
  const y = M.if_else(x, M.MX(1), M.MX(2));
  const f = M.Function("f", [x], [y]);
  let out = f(1);
  assertAlmostEqual(out.nonzeros()[0], 1, 10, "if_else(true)");
  out = f(0);
  assertAlmostEqual(out.nonzeros()[0], 2, 10, "if_else(false)");
}

function test_if_else_zero(M) {
  const x = M.MX.sym("x", 1n, 1n);
  const y = M.if_else(x, M.MX(5), M.MX(0));
  const f = M.Function("f", [x], [y]);
  let out = f(1);
  assertAlmostEqual(out.nonzeros()[0], 5, 10, "if_else_zero(1)");
  out = f(0);
  assertAlmostEqual(out.nonzeros()[0], 0, 10, "if_else_zero(0)");
}

function test_is_regular(M) {
  const a = M.MX([0, 1]);
  assertEqual(a.is_regular(), true,  "[0,1] is_regular");
  const b = M.MX([0, Infinity]);
  assertEqual(b.is_regular(), false, "[0, inf] is_regular");
}

function test_vertcat(M) {
  const X = M.MX.sym("X", 10n, 1n);
  const T = M.vcat([
    X.get(false, M.Slice(4n, false)),
    X.get(false, M.Slice(2n, false)),
  ]);
  const q = M.times(T, T);
  const f = M.Function("f", [X], [q]);
  const out = f([0,1,2,3,4,5,6,7,8,9]);
  dmArrayEqual(out, [16, 4], 10, "f([0..9])");
}

function test_symvar(M) {
  const a = M.MX.sym("a", 1n, 1n);
  const b = M.MX.sym("b", 1n, 1n);
  const c = M.MX.sym("c", 1n, 1n);
  const e = M.plus(M.cos(M.times(a, b)), c);
  const w = M.symvar(e);
  assertEqual(w.length, 3, "symvar length");
}

function test_pow(M) {
  const y = M.MX(4n, 1n);
  const p = M.power(y, M.MX(0));
  assertEqual(p.nnz(), 4, "(empty 4x1)**0 nnz");
}

function test_contains(M) {
  const x = M.MX.sym("x", 2n, 2n);
  const y = M.MX.sym("y", 1n, 1n);
  const z = M.MX.sym("z", 3n, 3n);
  assertEqual(M.contains([x, y, z], x), true,  "contains x");
  assertEqual(M.contains([x, y],    z), false, "!contains z");
  assertEqual(M.contains_any([x, y], [y, z]), true,  "contains_any [y,z]");
  assertEqual(M.contains_all([x, y], [y, z]), false, "!contains_all [y,z]");
}

function test_null(M) {
  const x = M.MX.sym("x", 1n, 1n);
  const empty = M.MX(0n, 0n);
  const f = M.Function("f", [x],
    [M.power(x, M.MX(2)), empty]);
  assertEqual(f.size1_out(1n), 0n, "size1_out(1)");
  assertEqual(f.size2_out(1n), 0n, "size2_out(1)");
}

function test_constmxmtimes(M) {
  // No assertion -- just exercise expression building.
  const _ = M.times(M.MX(0.1), M.MX.ones(2n, 1n));
}

function test_blockcat(M) {
  const x = M.MX.sym("x", 1n, 1n);
  const y = M.blockcat(x, M.times(M.MX(2), x),
                          M.times(M.MX(3), x),
                          M.times(M.MX(4), x));
  const f = M.Function("f", [x], [y]);
  const out = f(3);
  // Column-major flatten of [[3,6],[9,12]] = [3, 9, 6, 12].
  dmArrayEqual(out, [3, 9, 6, 12], 10, "blockcat 2x2 at x=3");
}

function test_MXvec(M) {
  const u = M.DM.ones(2n, 3n);
  const U = M.MX.sym("u", 2n, 3n);
  const f = M.Function("f", [U], [M.vec(U)]);
  const out = f(u);
  const out0 = out;
  assertEqual(out0.size1(), 6n, "vec.size1");
  assertEqual(out0.size2(), 1n, "vec.size2");
}

function test_unite(M) {
  const x = M.MX(3n, 4n);
  const y = M.MX.sym("y", 3n, 4n);
  const z = M.unite(x, y);
  const f = M.Function("f", [y], [z]);
  const out = f(M.DM.ones(3n, 4n));
  assertEqual(out.nnz(), 12, "unite nnz");
}

function test_norms(M) {
  const a = M.MX(3);
  const ni = M.Function("f", [], [M.norm_inf(a)]);
  const out = ni();
  assertAlmostEqual(out.nonzeros()[0], 3, 10, "norm_inf(3)");
}

function test_MXd_substractionl(M) {
  const V = M.MX.sym("V", 1n, 1n);
  const X = M.MX.sym("X", 1n, 1n);
  const g1 = M.Function("g", [],
    [M.jacobian(M.minus(X, V), X, null)]);
  const o1 = g1();
  assertAlmostEqual(o1.nonzeros()[0], 1, 10, "d(X-V)/dX");
  const g2 = M.Function("g", [],
    [M.jacobian(M.minus(X, V), V, null)]);
  const o2 = g2();
  assertAlmostEqual(o2.nonzeros()[0], -1, 10, "d(X-V)/dV");
}

function test_reshape(M) {
  const X = M.MX.sym("X", 10n, 1n);
  const Y = M.reshape(X, 2n, 5n);
  assertEqual(Y.size1(), 2n, "reshape size1");
  assertEqual(Y.size2(), 5n, "reshape size2");
}

function test_substitute(M) {
  const x = M.MX.sym("x", 1n, 1n);
  const y = M.times(x, x);
  const r = M.substitute(y, x, M.MX(3));
  const f = M.Function("f", [], [r]);
  const o = f();
  assertAlmostEqual(o.nonzeros()[0], 9, 10, "subst x=3 in x^2");
}

function test_empty_broadcast(M) {
  for (const nc of [0n, 2n]) {
    const c = M.MX.sym("c", nc, 1n);
    const t = M.MX.sym("t", nc, 3n);
    const res = M.atan2(c, t);
    assertEqual(res.size1(), nc, `atan2 size1 nc=${nc}`);
    assertEqual(res.size2(), 3n, `atan2 size2 nc=${nc}`);
  }
}

function test_truth(M) {
  assertEqual(M.MX(1).is_one(),  true,  "MX(1).is_one");
  assertEqual(M.MX(0).is_zero(), true,  "MX(0).is_zero");
  assertEqual(M.MX(1).is_zero(), false, "MX(1).is_zero");
  assertEqual(M.MX(0).is_one(),  false, "MX(0).is_one");
}

function test_diagcat(M) {
  const a = M.MX.ones(3n, 4n);
  const b = M.MX.ones(3n, 4n);
  const c = M.dcat([a, b]);
  assertEqual(c.size1(), 6n, "diagcat rows");
  assertEqual(c.size2(), 8n, "diagcat cols");
}

function test_ticket(M) {
  const x = M.MX.sym("x", 1n, 1n);
  const _ = M.plus(x, M.MX(0n, 0n));
}

function test_pow_zero(M) {
  const z = M.power(M.MX(0), M.MX(0));
  const f = M.Function("f", [], [z]);
  const out = f();
  assertAlmostEqual(out.nonzeros()[0], 1, 10, "0**0 = 1");
}

function test_vertcat_empty(M) {
  const a1 = M.MX(0n, 2n);
  const v1 = M.vcat([a1, a1]);
  assertEqual(v1.size1(), 0n, "0x2 stack rows");
  assertEqual(v1.size2(), 2n, "0x2 stack cols");
  const a2 = M.MX(2n, 0n);
  const v2 = M.vcat([a2, a2]);
  assertEqual(v2.size1(), 4n, "2x0 stack rows");
  assertEqual(v2.size2(), 0n, "2x0 stack cols");
}

function test_jacobian_empty(M) {
  const x = M.MX.sym("x", 3n, 1n);
  const J1 = M.jacobian(M.MX(0n, 0n), x, null);
  assertEqual(J1.size1(), 0n, "J(empty,x) rows");
  assertEqual(J1.size2(), 3n, "J(empty,x) cols");
  const J2 = M.jacobian(x, M.MX.sym("x", 0n, 4n), null);
  assertEqual(J2.size1(), 3n, "J(x,empty) rows");
  assertEqual(J2.size2(), 0n, "J(x,empty) cols");
}

function test_kron(M) {
  const A = M.MX.sym("A", 3n, 2n);
  const B = M.MX.sym("B", 4n, 3n);
  const C = M.kron(A, B);
  assertEqual(C.size1(), 12n, "kron rows");
  assertEqual(C.size2(),  6n, "kron cols");
}

function test_depends_on(M) {
  const a = M.MX.sym("a", 1n, 1n);
  const b = M.MX.sym("b", 1n, 1n);
  assertEqual(M.depends_on(M.times(a, a), a), true,  "depends_on a^2,a");
  assertEqual(M.depends_on(M.times(a, a), b), false, "!depends_on a^2,b");
}

function test_vertsplit(M) {
  const a = M.MX.sym("a", 5n, 5n);
  const v = M.vertsplit(a, [0n, 2n, 4n, 5n]);
  assertEqual(v.length, 3, "vertsplit length");
}

function test_horzsplit(M) {
  const a = M.MX.sym("a", 5n, 5n);
  const v = M.horzsplit(a, [0n, 2n, 4n, 5n]);
  assertEqual(v.length, 3, "horzsplit length");
}

function test_mxnulloutput(M) {
  const a = M.MX(5n, 0n);
  const b = M.MX.sym("b", 2n, 1n);
  const f = M.Function("f", [b], [a]);
  assertEqual(f.size1_out(0n), 5n, "size1_out=5");
  assertEqual(f.size2_out(0n), 0n, "size2_out=0");
}

function test_MXorder(M) {
  // mx.py:test_MXorder -- output of f(x) = x+x must preserve element order.
  const x = M.MX.sym("x", 2n, 3n);
  const f = M.Function("f", [x], [M.plus(x, x)]);
  assertEqual(f.n_in(),  1n, "n_in");
  assertEqual(f.n_out(), 1n, "n_out");
  const in0 = M.DM(f.sparsity_in(0n), [1, 2, 3, 4, 5, 6]);
  const out = f(in0);
  dmArrayEqual(out, [2, 4, 6, 8, 10, 12], 10, "x+x");
}

function test_MXreshape(M) {
  // mx.py:test_MXreshape -- reshape MX from 2x3 to 1x6.
  const x = M.MX.sym("x", 2n, 3n);
  const z = M.reshape(x, 1n, 6n);
  assertEqual(z.size1(), 1n, "reshape rows");
  assertEqual(z.size2(), 6n, "reshape cols");
  const f = M.Function("f", [x], [z]);
  const in0 = M.DM(f.sparsity_in(0n), [1, 2, 3, 4, 5, 6]);
  const out = f(in0);
  dmArrayEqual(out, [1, 2, 3, 4, 5, 6], 10, "reshape values");
}

function test_mx_in(M) {
  // mx.py:test_mx_in -- compose two Function objects via mx_in / call.
  const x = M.MX.sym("x", 2n, 3n);
  const f = M.Function("f", [x], [M.times(M.MX(3), x)]);
  const x_in = f.mx_in();
  const x_out = f.call(x_in);
  const g = M.Function("g", [x_in[0]], [M.times(M.MX(6), x_out[0])]);
  const n = [1, 2, 3, 4, 5, 6];
  const f_in = M.DM(f.sparsity_in(0n), n);
  const g_in = M.DM(g.sparsity_in(0n), n);
  const f_out = f(f_in);
  const g_out = g(g_in);
  const sixfo = f_out.nonzeros().map((v) => 6 * v);
  dmArrayEqual(g_out, sixfo, 10, "g = 6*f");
}

function test_indexingOutOfBounds(M) {
  // mx.py:test_indexingOutOfBounds -- only some patterns. Demonstrate the
  // error-throwing direction; mirrors the DM y[12, 0] case via wasm-js
  // slice indexing.
  const y = M.DM.zeros(4n, 5n);
  let threw = false;
  try { y.get(false, M.Slice(12n, 13n), M.Slice(0n, 1n)); } catch (e) { threw = true; }
  if (!threw) throw new Error("y[12, 0] should throw");
}

// ============================================================================
//  runner
// ============================================================================

const tests = [
  ["test_MX1",                test_MX1],
  ["test_MXvertcat",          test_MXvertcat],
  ["test_MX_fun1",            test_MX_fun1],
  ["test_MXfunction2",        test_MXfunction2],
  ["test_MXfunction3",        test_MXfunction3],
  ["test_MXfunction3b",       test_MXfunction3b],
  ["test_MXfunction4",        test_MXfunction4],
  ["test_MXfunction5",        test_MXfunction5],
  ["test_identitySX",         test_identitySX],
  ["test_identityMX",         test_identityMX],
  ["test_MXtrans",            test_MXtrans],
  ["test_issue107",           test_issue107],
  ["test_MXd_trivial",        test_MXd_trivial],
  ["test_issue184",           test_issue184],
  ["test_jacobian_tools",     test_jacobian_tools],
  ["test_if_else",            test_if_else],
  ["test_if_else_zero",       test_if_else_zero],
  ["test_is_regular",         test_is_regular],
  ["test_vertcat",            test_vertcat],
  ["test_symvar",             test_symvar],
  ["test_pow",                test_pow],
  ["test_contains",           test_contains],
  ["test_null",               test_null],
  ["test_constmxmtimes",      test_constmxmtimes],
  ["test_blockcat",           test_blockcat],
  ["test_MXvec",              test_MXvec],
  ["test_unite",              test_unite],
  ["test_norms",              test_norms],
  ["test_MXd_substractionl",  test_MXd_substractionl],
  ["test_reshape",            test_reshape],
  ["test_substitute",         test_substitute],
  ["test_empty_broadcast",    test_empty_broadcast],
  ["test_truth",              test_truth],
  ["test_diagcat",            test_diagcat],
  ["test_ticket",             test_ticket],
  ["test_pow_zero",           test_pow_zero],
  ["test_vertcat_empty",      test_vertcat_empty],
  ["test_jacobian_empty",     test_jacobian_empty],
  ["test_kron",               test_kron],
  ["test_depends_on",         test_depends_on],
  ["test_vertsplit",          test_vertsplit],
  ["test_horzsplit",          test_horzsplit],
  ["test_mxnulloutput",       test_mxnulloutput],
  ["test_MXorder",            test_MXorder],
  ["test_MXreshape",          test_MXreshape],
  ["test_mx_in",              test_mx_in],
  ["test_indexingOutOfBounds", test_indexingOutOfBounds],
];

(async () => {
  const create = require(modulePath);
  const M = await create();
  let pass = 0, fail = 0;
  for (const [name, fn] of tests) {
    try {
      await fn(M);
      console.log(`ok -- ${name}`);
      pass++;
    } catch (e) {
      console.error(`FAIL ${name}: ${e.message}`);
      fail++;
    }
  }
  console.log(`\n${pass}/${pass + fail} passed`);
  process.exit(fail === 0 ? 0 : 1);
})().catch((e) => { console.error("FATAL:", e); process.exit(1); });

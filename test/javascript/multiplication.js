// multiplication.js -- JS port of test/python/multiplication.py (subset).
// Most tests require numpy + kernel_used helper; this subset covers
// the portable mtimes correctness checks.

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

function assertTrue(cond, msg) { if (!cond) throw new Error(`${msg}: expected truthy`); }

// ----- tests -----

function test_mtimes_dm_dense(M) {
  // 2x3 * 3x2 = 2x2
  const A = M.DM([[1, 2, 3], [4, 5, 6]]);
  const B = M.DM([[7, 8], [9, 10], [11, 12]]);
  const C = M.mtimes(A, B);
  assertEqual(C.size1(), 2n, "size1");
  assertEqual(C.size2(), 2n, "size2");
  // Expected (column-major): [[58, 64], [139, 154]]
  // C(0,0) = 1*7 + 2*9 + 3*11 = 58
  // C(0,1) = 1*8 + 2*10 + 3*12 = 64
  // C(1,0) = 4*7 + 5*9 + 6*11 = 139
  // C(1,1) = 4*8 + 5*10 + 6*12 = 154
  // Column-major: [58, 139, 64, 154]
  assertArrayAlmostEqual(C.nonzeros(), [58, 139, 64, 154], 10, "mtimes values");
}

function test_mtimes_mx_eval(M) {
  // Build a symbolic mtimes, eval with a Function.
  const A = M.MX.sym("A", 2n, 3n);
  const B = M.MX.sym("B", 3n, 2n);
  const C = M.mtimes(A, B);
  const f = M.Function("f", [A, B], [C]);
  // Evaluate with concrete matrices.
  const Av = M.DM([[1, 2, 3], [4, 5, 6]]);
  const Bv = M.DM([[7, 8], [9, 10], [11, 12]]);
  const out = f.call([Av, Bv]);
  assertArrayAlmostEqual(out.nonzeros(), [58, 139, 64, 154], 8, "MX mtimes eval");
}

function test_mtimes_sparse_dense(M) {
  // Sparse 3x3 lower-triangular times dense 3x2.
  const sp = M.Sparsity.lower(3n);
  // Build DM with the lower-triangular sparsity, all-ones values.
  // Using sparsify path: DM([[1,0,0],[1,1,0],[1,1,1]]) then sparsify drops zeros.
  const A_dense = M.DM([[1, 0, 0], [1, 1, 0], [1, 1, 1]]);
  const A = M.sparsify(A_dense);
  const B = M.DM([[1, 2], [3, 4], [5, 6]]);
  const C = M.mtimes(A, B);
  assertEqual(C.size1(), 3n, "sparse mtimes size1");
  assertEqual(C.size2(), 2n, "sparse mtimes size2");
  // Expected dense result: [[1,2], [4,6], [9,12]]
  // Column-major: [1, 4, 9, 2, 6, 12]
  // (After sparsify result is dense.)
  const expect_col_major = [1, 4, 9, 2, 6, 12];
  assertArrayAlmostEqual(C.nonzeros(), expect_col_major, 8, "sparse mtimes values");
}

function test_is_compactible_dense(M) {
  const sp = M.Sparsity.dense(3n, 4n);
  const out = sp.is_compactible();
  // out is a [bool, vec<int>, vec<int>] tuple wrapped as array.
  // The signature returns is_compact, row_support, col_support.
  assertTrue(Array.isArray(out) ? out[0] : !!out, "is_compactible(dense) -> true");
}

// ----- driver -----

const tests = [
  ["test_mtimes_dm_dense",      test_mtimes_dm_dense],
  ["test_mtimes_mx_eval",       test_mtimes_mx_eval],
  ["test_mtimes_sparse_dense",  test_mtimes_sparse_dense],
  ["test_is_compactible_dense", test_is_compactible_dense],
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

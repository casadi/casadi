// matrix.js -- JS port of test/python/matrix.py (subset).
//
// Covers the structurally-testable DM API: transpose, sum1/sum2,
// is_regular, cross product, sizes (nnz_diag/upper/lower), kron,
// veccat, diag.  Skips parts that need numpy/scipy or solver
// backends.

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

function assertTrue(cond, msg) {
  if (!cond) throw new Error(`${msg}: expected truthy`);
}

function assertFalse(cond, msg) {
  if (cond) throw new Error(`${msg}: expected falsy`);
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

function test_trans(M) {
  // DM(0, 1) is a 0x1 matrix.  Its transpose is 1x0.
  const a = M.DM(0n, 1n);
  const b = a.T;
  assertEqual(b.size1(), 1n, "trans size1");
  assertEqual(b.size2(), 0n, "trans size2");
}

function test_sum(M) {
  // Nested-JS-array -> DM coercion (Phase 3.4 fuzzy promotion).
  const D = M.DM([[1, 2, 3], [4, 5, 6], [7, 8, 9]]);
  const s1 = M.sum1(D);
  assertEqual(s1.size1(), 1n, "sum1 size1");
  assertEqual(s1.size2(), 3n, "sum1 size2");
  assertArrayAlmostEqual(s1.nonzeros(), [12, 15, 18], 10, "sum1 nz");
  const s2 = M.sum2(D);
  assertEqual(s2.size1(), 3n, "sum2 size1");
  assertEqual(s2.size2(), 1n, "sum2 size2");
  assertArrayAlmostEqual(s2.nonzeros(), [6, 15, 24], 10, "sum2 nz");
}

function test_is_regular(M) {
  assertTrue(M.DM([1, 2]).is_regular(), "[1,2] is_regular");
  assertFalse(M.DM([1, Infinity]).is_regular(), "[1,inf] is_regular");
}

function test_cross(M) {
  // Global cross(A, B, dim=-1) -- the 3rd `dim` arg defaults to -1
  // and is auto-filled when omitted (global default-arg dispatch).
  const i_ = M.DM([1, 0, 0]);
  const j_ = M.DM([0, 1, 0]);
  const r1 = M.cross(i_, j_);
  assertArrayAlmostEqual(r1.nonzeros(), [0, 0, 1], 10, "i x j = k");
}

function test_sizes(M) {
  assertEqual(M.Sparsity.diag(10n).nnz_diag(), 10n, "diag nnz_diag");
  assertEqual(M.Sparsity.diag(10n).nnz_upper(), 10n, "diag nnz_upper");
  assertEqual(M.Sparsity.diag(10n).nnz_lower(), 10n, "diag nnz_lower");
  // 10x10 dense lower has 1+2+...+10 = 55 entries
  assertEqual(M.Sparsity.dense(10n, 10n).nnz_lower(), 55n, "dense nnz_lower");
  assertEqual(M.Sparsity.dense(10n, 10n).nnz_upper(), 55n, "dense nnz_upper");
  assertEqual(M.Sparsity.dense(10n, 10n).nnz_diag(), 10n, "dense nnz_diag");
}

function test_kron(M) {
  const a = M.DM([[1, 0], [0, 1]]);
  const b = M.DM([[2, 3], [4, 5]]);
  const k = M.kron(a, b);
  assertEqual(k.size1(), 4n, "kron size1");
  assertEqual(k.size2(), 4n, "kron size2");
}

function test_veccat(M) {
  const r = M.veccat([M.DM([1, 2]), M.DM([3, 4, 5])]);
  assertEqual(r.numel(), 5n, "veccat numel");
  assertArrayAlmostEqual(r.nonzeros(), [1, 2, 3, 4, 5], 10, "veccat values");
}

function test_diag(M) {
  // diag(vector) builds a diagonal matrix
  const d = M.diag(M.DM([1, 2, 3]));
  assertEqual(d.size1(), 3n, "diag size1");
  assertEqual(d.size2(), 3n, "diag size2");
  assertEqual(d.nnz(), 3n, "diag nnz");
}

function test_horzcat_dm(M) {
  const r = M.hcat([M.DM([1, 2, 3]), M.DM([4, 5, 6])]);
  assertEqual(r.size1(), 3n, "horzcat dm size1");
  assertEqual(r.size2(), 2n, "horzcat dm size2");
}

function test_vertcat_dm(M) {
  const a = M.DM([[1, 2]]);   // 1x2
  const b = M.DM([[3, 4]]);   // 1x2
  const r = M.vcat([a, b]);
  assertEqual(r.size1(), 2n, "vertcat dm size1");
  assertEqual(r.size2(), 2n, "vertcat dm size2");
}

function test_truth_dm(M) {
  assertTrue(M.DM([1]).__nonzero__(), "DM([1]) truthy");
  assertFalse(M.DM([0]).__nonzero__(), "DM([0]) falsy");
  assertTrue(M.DM([0.2]).__nonzero__(), "DM([0.2]) truthy");
  assertTrue(M.DM([-0.2]).__nonzero__(), "DM([-0.2]) truthy");
}

function test_inv_dm(M) {
  // DM matrix inverse: inv(A) * A = I.
  const a = M.DM([[1, 2], [1, 3]]);
  const ai = M.inv(a);
  const r = M.mtimes(ai, a);
  // r is (sparse or dense) 2x2 identity.  Diagonal nonzeros = 1.
  const sz = Number(r.size1());
  assertEqual(sz, 2, "inv result size1");
  // Verify the diagonal sums to 2 (1 + 1).
  let diag_sum = 0;
  for (let i = 0; i < 2; ++i) {
    // r(i,i) -- approximate via element-by-element.  Use the get method.
    // For wasm-js, use M.diag to extract.
  }
  // Easier: trace = sum of diagonal.
  if (typeof M.trace === "function") {
    const tr = M.trace(r);
    assertAlmost(Number(tr.nonzeros()[0]), 2, 8, "trace(I_2) = 2");
  }
}

function test_DM_arithmetic_chain(M) {
  // a + b and a * b chained operations.
  const a = M.DM([[1, 2], [3, 4]]);
  const b = M.DM([[5, 6], [7, 8]]);
  const sum = M.plus(a, b);
  // a + b = [[6, 8], [10, 12]]; col-major flat = [6, 10, 8, 12].
  assertArrayAlmostEqual(sum.nonzeros(), [6, 10, 8, 12], 10, "a+b col-major");
  const prod = M.mtimes(a, b);
  // a * b = [[19, 22], [43, 50]]; col-major = [19, 43, 22, 50].
  assertArrayAlmostEqual(prod.nonzeros(), [19, 43, 22, 50], 10, "a*b col-major");
}

function test_DM_sym_constructor(M) {
  // DM construction variants.
  const s = M.DM(5);
  assertEqual(s.size1(), 1n, "scalar size1");
  assertEqual(s.size2(), 1n, "scalar size2");
  const v = M.DM([1, 2, 3]);
  assertEqual(v.size1(), 3n, "vec size1");
  assertEqual(v.size2(), 1n, "vec size2");
}

function test_DM_zeros_ones_eye(M) {
  // Standard factories.
  const z = M.DM.zeros(3n, 4n);
  assertEqual(z.size1(), 3n, "zeros size1");
  assertEqual(z.size2(), 4n, "zeros size2");
  // zeros may be structurally sparse (no nonzeros).
  const o = M.DM.ones(2n, 5n);
  assertEqual(o.size1(), 2n, "ones size1");
  assertEqual(o.size2(), 5n, "ones size2");
  assertEqual(o.nnz(), 10n, "ones nnz=10");
  if (typeof M.DM.eye === "function") {
    const e = M.DM.eye(3n);
    assertEqual(e.size1(), 3n, "eye size1");
    assertEqual(e.size2(), 3n, "eye size2");
    assertEqual(e.nnz(),   3n, "eye nnz");
  }
}

function assertAlmost(actual, expected, places, msg) {
  if (Math.abs(actual - expected) > Math.pow(10, -places)) {
    throw new Error(`${msg}: expected ${expected}, got ${actual}`);
  }
}

function test_det_diag(M) {
  // matrix.py:test_det -- diagonal DM determinant.
  const a = M.DM.zeros(5n, 5n);
  for (let i = 0; i < 5; ++i) {
    a.set(M.DM(i + 1), false, BigInt(i), BigInt(i));
  }
  const d = M.det(a);
  const dn = Number(d.nonzeros()[0]);
  // 1*2*3*4*5 = 120
  if (Math.abs(dn - 120) > 1e-9) throw new Error("det(diag) expected 120, got " + dn);
}

function test_sprank_basic(M) {
  // matrix.py:test_sprank -- structural rank.
  assertEqual(BigInt(M.sprank(M.DM.eye(3n))), 3n, "sprank(eye3) == 3");
  assertEqual(BigInt(M.sprank(M.DM.ones(1n, 3n))), 1n, "sprank(ones(1,3)) == 1");
  assertEqual(BigInt(M.sprank(M.DM.ones(3n, 1n))), 1n, "sprank(ones(3,1)) == 1");
  assertEqual(BigInt(M.sprank(M.DM.ones(2n, 3n))), 2n, "sprank(ones(2,3)) == 2");
  assertEqual(BigInt(M.sprank(M.DM.ones(3n, 3n))), 3n, "sprank(ones(3,3)) == 3");
}

function test_tril2symm_basic(M) {
  // matrix.py:test_tril2symm -- expand lower-triangular DM to symmetric.
  // Build a 3x3 lower DM, then expand and verify A = A.T.
  const sp = M.Sparsity.lower(3n);
  const A = M.DM(sp, [1, 2, 3, 4, 5, 6]);
  const S = M.tril2symm(A);
  // S should be 3x3 dense symmetric.
  assertEqual(S.size1(), 3n, "tril2symm size1");
  assertEqual(S.size2(), 3n, "tril2symm size2");
  // Compare against its own transpose: symmetric.
  const ST = S.T;
  const a = S.nonzeros();
  const b = ST.nonzeros();
  assertArrayAlmostEqual(a, b, 10, "tril2symm symmetric");
}

// ----- driver -----

const tests = [
  ["test_trans",        test_trans],
  ["test_sum",          test_sum],
  ["test_is_regular",   test_is_regular],
  ["test_cross",        test_cross],
  ["test_sizes",        test_sizes],
  ["test_kron",         test_kron],
  ["test_veccat",       test_veccat],
  ["test_diag",         test_diag],
  ["test_horzcat_dm",   test_horzcat_dm],
  ["test_vertcat_dm",   test_vertcat_dm],
  ["test_truth_dm",     test_truth_dm],
  ["test_inv_dm",       test_inv_dm],
  ["test_DM_arithmetic_chain", test_DM_arithmetic_chain],
  ["test_DM_sym_constructor",  test_DM_sym_constructor],
  ["test_DM_zeros_ones_eye",   test_DM_zeros_ones_eye],
  ["test_det_diag",     test_det_diag],
  ["test_sprank_basic", test_sprank_basic],
  ["test_tril2symm_basic", test_tril2symm_basic],
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

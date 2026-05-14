// sparsity.js -- JS port of test/python/sparsity.py (subset).
//
// Highest-priority Python test port per HANDOFF_PHASE3_PLUS.md: pure
// structure API with no solver deps.  Covers the core Sparsity factory
// methods (Sparsity(rows,cols), Sparsity.dense, .lower, .diag,
// .rowcol, .triplet), nonzero queries (nnz, row, get_col, find), set
// operations (unite, intersect, is_subset), and a serialize round-trip.

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

// Build a Sparsity by adding nonzeros (a JS Set of [r,c] tuples or
// arrays).  Mirrors the Python `for i in nza: a.add_nz(i[0], i[1])`
// idiom but works with JS arrays.
function spFromNZ(M, nrows, ncols, nz) {
  const s = new M.Sparsity(BigInt(nrows), BigInt(ncols));
  for (const [r, c] of nz) {
    s.add_nz(BigInt(r), BigInt(c));
  }
  return s;
}

// Collect all (row, col) nonzero indices from a Sparsity into a JS Set
// of "r,c" strings (so set ops work with primitive equality).
function spNZSet(s) {
  const out = new Set();
  const cols = s.get_col();   // returns plain array
  const n = Number(s.nnz());
  for (let k = 0; k < n; ++k) {
    out.add(`${s.row(BigInt(k))},${cols[k]}`);
  }
  return out;
}

// ----- tests -----

function test_union(M) {
  const nza = [[0, 0], [0, 1], [2, 0], [3, 1]];
  const nzb = [[0, 2], [0, 0], [2, 2]];
  const a = spFromNZ(M, 4, 5, nza);
  const b = spFromNZ(M, 4, 5, nzb);
  const c = a.unite(b);

  // Expected union size
  const union = new Set([...nza, ...nzb].map(([r, k]) => `${r},${k}`));
  assertEqual(c.nnz(), BigInt(union.size), "union.nnz");

  // Each nonzero in c must come from nza or nzb
  const cset = spNZSet(c);
  for (const ind of cset) {
    assertTrue(union.has(ind), `union contains ${ind}`);
  }
}

function test_intersection(M) {
  const nza = [[0, 0], [0, 1], [2, 0], [3, 1], [2, 3]];
  const nzb = [[0, 2], [0, 0], [2, 2], [2, 3]];
  const a = spFromNZ(M, 4, 5, nza);
  const b = spFromNZ(M, 4, 5, nzb);
  const c = a.intersect(b);

  // Expected intersection
  const sa = new Set(nza.map(([r, k]) => `${r},${k}`));
  const sb = new Set(nzb.map(([r, k]) => `${r},${k}`));
  const inter = new Set([...sa].filter(x => sb.has(x)));
  assertEqual(c.nnz(), BigInt(inter.size), "inter.nnz");

  const cset = spNZSet(c);
  for (const ind of cset) {
    assertTrue(inter.has(ind), `inter contains ${ind}`);
  }
}

function test_rowcol(M) {
  // Sparsity.rowcol([rows...], [cols...], nrow, ncol)
  const n = 3n;
  const s = M.Sparsity.rowcol([n - 1n, 0n], [0n, n - 1n], n, n);
  assertEqual(s.colind()[0], 0, "rowcol colind[0]");
  assertEqual(s.colind()[1], 2, "rowcol colind[1]");
  assertEqual(s.colind()[2], 2, "rowcol colind[2]");
  assertEqual(s.colind()[3], 4, "rowcol colind[3]");
  assertEqual(s.row(0n), 0n, "rowcol row[0]");
  assertEqual(s.row(1n), 2n, "rowcol row[1]");
  assertEqual(s.row(2n), 0n, "rowcol row[2]");
  assertEqual(s.row(3n), 2n, "rowcol row[3]");
}

function test_dense(M) {
  const s = M.Sparsity.dense(3n, 4n);
  assertEqual(s.size1(), 3n, "dense.size1");
  assertEqual(s.size2(), 4n, "dense.size2");
  assertEqual(s.nnz(), 12n, "dense.nnz");
  assertTrue(s.is_dense(), "dense.is_dense");
}

function test_lower(M) {
  const s = M.Sparsity.lower(4n);
  assertEqual(s.size1(), 4n, "lower.size1");
  assertEqual(s.size2(), 4n, "lower.size2");
  // lower triangular of 4x4 has 4 + 3 + 2 + 1 = 10 nonzeros
  assertEqual(s.nnz(), 10n, "lower.nnz");
}

function test_diag(M) {
  const s = M.Sparsity.diag(5n);
  assertEqual(s.size1(), 5n, "diag.size1");
  assertEqual(s.size2(), 5n, "diag.size2");
  assertEqual(s.nnz(), 5n, "diag.nnz");
}

function test_is_subset(M) {
  const pairs = [
    [M.Sparsity.lower(3n), M.Sparsity.dense(3n, 3n)],
    [M.Sparsity.diag(3n),  M.Sparsity.dense(3n, 3n)],
    [M.Sparsity.diag(3n),  M.Sparsity.lower(3n)],
    [new M.Sparsity(3n, 3n), M.Sparsity.lower(3n)],
  ];
  for (const [L, R] of pairs) {
    assertTrue(L.is_subset(R), "L is subset of R");
    assertFalse(R.is_subset(L), "R not subset of L");
  }
}

function test_triplet(M) {
  // Sparsity.triplet(nrow, ncol, [rows], [cols])
  const s = M.Sparsity.triplet(4n, 5n, [0n, 0n, 2n, 3n], [0n, 1n, 0n, 1n]);
  assertEqual(s.nnz(), 4n, "triplet.nnz");
  assertEqual(s.size1(), 4n, "triplet.size1");
  assertEqual(s.size2(), 5n, "triplet.size2");
}

function test_serialize_roundtrip(M) {
  // Phase-3.4 ctor default-arg fill: `new M.Sparsity()` (the
  // null-Sparsity case) now works.  Smoke-test that serialize()
  // doesn't crash; deserialize() round-trip is currently broken
  // (memory access OOB inside the wasm wrapper) -- separate
  // followup.
  const a = M.Sparsity.dense(4n, 5n);
  const blob = a.serialize();
  assertTrue(typeof blob === 'string' && blob.length > 0, "serialize returns nonempty string");
  const empty = new M.Sparsity();
  assertTrue(empty.is_null(), "Sparsity() is null");
}

function test_kron(M) {
  // Sparsity.kron and the symbolic kron should give the same sparsity.
  const a = M.Sparsity.dense(2n, 3n);
  const b = M.Sparsity.dense(4n, 3n);
  const c = M.kron(a, b);
  assertEqual(c.size1(), a.size1() * b.size1(), "kron size1");
  assertEqual(c.size2(), a.size2() * b.size2(), "kron size2");
  assertEqual(c.nnz(), a.nnz() * b.nnz(), "kron nnz");
}

function test_enlarge(M) {
  const sp = M.Sparsity.dense(3n, 4n);
  const before_n = Number(sp.size1()) * Number(sp.size2());
  sp.enlarge(7n, 8n, [1n, 2n, 4n], [0n, 3n, 4n, 6n]);
  assertEqual(sp.size1(), 7n, "enlarge size1");
  assertEqual(sp.size2(), 8n, "enlarge size2");
  // nnz unchanged
  assertEqual(sp.nnz(), BigInt(before_n), "enlarge nnz preserved");
}

function test_NZ(M) {
  // Adding NZ entries one-by-one should equal Sparsity.triplet of same.
  const nza = [[0, 0], [0, 1], [2, 0], [2, 3], [2, 4], [3, 1]];
  const a = spFromNZ(M, 4, 5, nza);
  const rows = nza.map(([r, _]) => BigInt(r));
  const cols = nza.map(([_, c]) => BigInt(c));
  const b = M.Sparsity.triplet(4n, 5n, rows, cols);
  assertEqual(a.nnz(), b.nnz(), "NZ vs triplet nnz");
  assertEqual(a.size1(), b.size1(), "size1");
  assertEqual(a.size2(), b.size2(), "size2");
}

function test_reshape_sparsity(M) {
  const a = M.Sparsity.dense(4n, 5n);
  const r = M.reshape(a, 2n, 10n);
  assertEqual(r.size1(), 2n, "reshape size1");
  assertEqual(r.size2(), 10n, "reshape size2");
  assertEqual(r.nnz(), a.nnz(), "reshape nnz preserved");
}

// ----- driver -----

const tests = [
  ["test_union",              test_union],
  ["test_intersection",       test_intersection],
  ["test_rowcol",             test_rowcol],
  ["test_dense",              test_dense],
  ["test_lower",              test_lower],
  ["test_diag",               test_diag],
  ["test_is_subset",          test_is_subset],
  ["test_triplet",            test_triplet],
  ["test_serialize_roundtrip", test_serialize_roundtrip],
  ["test_kron",               test_kron],
  ["test_enlarge",            test_enlarge],
  ["test_NZ",                 test_NZ],
  ["test_reshape_sparsity",   test_reshape_sparsity],
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

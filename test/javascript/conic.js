// conic.js -- subset port of test/python/conic.py.
//
// Exercises the conic() / qpsol() solver-construction API against the
// HiGHS plugin, which is loaded at runtime from libcasadi_conic_highs.so
// (a SIDE_MODULE built and preloaded into MEMFS by swig/wasm-js/CMakeLists.txt).
//
// Ports:
//   test_general_unconstrained  -- diag-H, no inequality A, free-x QP
//   test_general_convex_dense   -- standard 2D QP with 3 inequality rows
//   test_linear                 -- LP (H=0, 3 inequality rows)

const path = require("path");
const fs = require("fs");

const wasmDir = path.resolve(__dirname, "../../build-wasm/swig/wasm-js");
const modulePath = path.join(wasmDir, "casadi.js");
if (!fs.existsSync(modulePath)) {
  console.error(`SKIP: ${modulePath} not built.`);
  process.exit(2);
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
  for (let i = 0; i < actual.length; ++i) {
    assertAlmost(actual[i], expected[i], places, `${msg}[${i}]`);
  }
}

// ----- tests -----

function test_general_unconstrained(M) {
  // min  0.5 [x y]'*I*[x y] + [-0.7, -2.3]'*[x y]
  //  =>  optimum at x=0.7, y=2.3, cost = -0.5*(0.7^2 + 2.3^2) = -2.89
  const H = M.DM([[1, 0], [0, 1]]);
  const G = M.DM([[-0.7], [-2.3]]);
  const inf = Infinity;
  const LBX = M.DM([[-inf], [-inf]]);
  const UBX = M.DM([[inf], [inf]]);

  const qp = { h: H.sparsity() };
  const solver = M.conic("solver", "highs", qp,
    { highs: { output_flag: false } });
  const res = solver.call({ h: H, g: G, lbx: LBX, ubx: UBX });
  assertArrayAlmostEqual(res.x.nonzeros(), [0.7, 2.3], 5, "x");
  assertAlmost(res.cost.nonzeros()[0], -0.5 * (0.7 * 0.7 + 2.3 * 2.3),
               5, "cost");
}

function test_general_convex_dense(M) {
  // min  0.5 [x y]'*H*[x y] + g'*[x y]
  //  s.t. A*[x y] <= [2, 2, 3]   (and x,y >= 0)
  // From conic.py: optimum at x=2/3, y=4/3, cost = -8 - 2/9.
  const H = M.DM([[1, -1], [-1, 2]]);
  const G = M.DM([[-2], [-6]]);
  const A = M.DM([[1, 1], [-1, 2], [2, 1]]);
  const inf = Infinity;
  const LBA = M.DM([[-inf], [-inf], [-inf]]);
  const UBA = M.DM([[2], [2], [3]]);
  const LBX = M.DM([[0], [0]]);
  const UBX = M.DM([[inf], [inf]]);

  const qp = { h: H.sparsity(), a: A.sparsity() };
  const solver = M.conic("solver", "highs", qp,
    { highs: { output_flag: false } });
  const res = solver.call({
    h: H, g: G, a: A, lbx: LBX, ubx: UBX, lba: LBA, uba: UBA
  });
  assertArrayAlmostEqual(res.x.nonzeros(), [2 / 3, 4 / 3], 5, "x");
  assertAlmost(res.cost.nonzeros()[0], -8 - 2 / 9, 5, "cost");
}

function test_linear(M) {
  // Pure LP (H=0).  From conic.py: optimum at x=0.5, y=1.5, cost=2.5.
  const H = M.DM(2n, 2n);  // 2x2 of zeros (empty sparsity)
  const A = M.DM([[-1, 1], [1, 1], [1, -2]]);
  const G = M.DM([[2], [1]]);
  const inf = Infinity;
  const LBA = M.DM([[-inf], [2], [-inf]]);
  const UBA = M.DM([[1], [inf], [4]]);
  const LBX = M.DM([[-inf], [0]]);
  const UBX = M.DM([[inf], [inf]]);

  const qp = { h: H.sparsity(), a: A.sparsity() };
  const solver = M.conic("solver", "highs", qp,
    { highs: { output_flag: false } });
  const res = solver.call({
    h: H, g: G, a: A, lbx: LBX, ubx: UBX, lba: LBA, uba: UBA
  });
  assertArrayAlmostEqual(res.x.nonzeros(), [0.5, 1.5], 4, "x");
  assertAlmost(res.cost.nonzeros()[0], 2.5, 4, "cost");
}

// ----- driver -----

const tests = [
  ["test_general_unconstrained", test_general_unconstrained],
  ["test_general_convex_dense",  test_general_convex_dense],
  ["test_linear",                test_linear],
];

(async () => {
  const create = require(modulePath);
  const M = await create();
  if (!M.has_conic("highs")) {
    console.error("SKIP: HiGHS plugin not loaded.");
    process.exit(2);
  }
  let pass = 0, fail = 0;
  for (const [name, fn] of tests) {
    try { fn(M); console.log(`ok -- ${name}`); pass++; }
    catch (e) { console.log(`FAIL ${name}: ${e.message}`); fail++; }
  }
  console.log(`\n${pass}/${pass + fail} passed`);
  process.exit(fail === 0 ? 0 : 1);
})().catch((e) => { console.error("FATAL:", e); process.exit(1); });

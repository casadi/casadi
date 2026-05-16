// superscs_smoke.js -- smoke test for the SuperSCS Conic plugin.
//
// Until the SuperSCS source's Fortran-77 hidden CHARACTER-length args
// trap is patched, this test verifies the plugin only registers.  With
// the patch applied to {linAlg.c, cones.c}, the test actually solves a
// 2-D SOCP and verifies the optimum analytically.
//
// SOCP:
//   min  2x + y
//    s.t. || x-5 , y-7 ||_2  <=  4
//
// Optimum: x* = 5 - 8/sqrt(5),  y* = 7 - 4/sqrt(5)
// (test_SOCP first case from test/python/conic.py)

const { assertAlmost, assertArrayAlmostEqual, runTests } =
  require("./_helpers");

function test_superscs_register(M) {
  if (!M.has_conic("superscs")) {
    throw new Error("superscs plugin not registered");
  }
}

function test_superscs_socp_basic(M) {
  if (!M.has_conic("superscs")) {
    throw new Error("superscs plugin not registered");
  }
  // min 2x + y  s.t. ||x-5, y-7||_2 <= 4
  // Build h = soc([x-5, y-7], 4) at MX level then assemble qpsol input.
  const x = M.MX.sym("x");
  const y = M.MX.sym("y");
  const xy = M.vertcat(x, y);
  // soc(z, t) returns a "second-order cone" matrix block; M.soc exists.
  const z = M.vertcat(M.minus(x, M.MX(5)), M.minus(y, M.MX(7)));
  const h = M.soc(z, M.MX(4));
  const f = M.plus(M.times(M.MX(2), x), y);
  const qp = { h: h, x: xy, f: f };
  const opts = { superscs: { eps: 1e-9, do_super_scs: 1, verbose: 0 } };
  const solver = M.qpsol("S", "superscs", qp, opts);
  const r = solver.call({});
  const sqrt5 = Math.sqrt(5);
  assertArrayAlmostEqual(
    r.x.nonzeros(), [5 - 8 / sqrt5, 7 - 4 / sqrt5], 5,
    "superscs SOCP x"
  );
  assertAlmost(
    Number(r.f.nonzeros()[0]),
    10 - 16 / sqrt5 + 7 - 4 / sqrt5, 5,
    "superscs SOCP f"
  );
}

runTests([
  ["test_superscs_register",   test_superscs_register],
  ["test_superscs_socp_basic", test_superscs_socp_basic],
]);

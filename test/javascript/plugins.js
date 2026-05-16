// plugins.js -- exhaustive smoke tests covering every wasm-js side-loaded
// plugin.  Each test exercises a different plugin via its plugin-name
// string; the side modules ship preloaded into MEMFS by
// swig/wasm-js/CMakeLists.txt's wasmjs_side_plugin() invocations.
//
// The 21 plugins covered:
//   conic:       qrqp, ipqp, nlpsol (qp_to_nlp), highs
//   nlpsol:      sqpmethod, qrsqp, scpgen, feasiblesqpmethod, ipopt
//   integrator:  rk, collocation
//   linsol:      qr, ldl, tridiag, lsqr, symbolicqr
//   interpolant: linear, bspline
//   rootfinder:  newton, fast_newton, nlpsol (implicit_to_nlp)
//
// Most tests deliberately keep the problem trivial -- the goal is plugin
// coverage, not algorithm validation (per-algorithm tests live in the
// dedicated test files conic.js / integration.js / etc.).

const { assertEqual, assertTrue, assertAlmost, assertArrayAlmostEqual,
        runTests } = require("./_helpers");

// ---------------------------------------------------------------------
// linsol plugins:  Ax = b  with A = diag([1,2,3]) -> x = b ./ [1,2,3]
// ---------------------------------------------------------------------
function _linsol(M, plugin) {
  const A = M.DM([[1, 0, 0], [0, 2, 0], [0, 0, 3]]);
  const ls = M.Linsol("ls", plugin, A.sparsity());
  ls.sfact(A);
  ls.nfact(A);
  const b = M.DM([6, 8, 9]);
  const x = ls.solve(A, b, false);
  assertArrayAlmostEqual(x.nonzeros(), [6, 4, 3], 9, plugin + ": x = b/[1,2,3]");
}
function test_linsol_qr        (M) { _linsol(M, "qr"); }
function test_linsol_ldl       (M) {
  // ldl is symmetric-only; symmetric A = diag works.
  _linsol(M, "ldl");
}
function test_linsol_lsqr      (M) {
  // lsqr is iterative; same trivial diag input works.
  _linsol(M, "lsqr");
}
function test_linsol_symbolicqr(M) { _linsol(M, "symbolicqr"); }
function test_linsol_tridiag   (M) {
  // Tridiagonal solver: needs ACTUAL tri-banded (3 nonzeros per
  // interior column).  linsol_tridiag.cpp's solve() indexes
  // A[a_colind[i] + {0,1,2}] for sub/main/super-diag, so a diag-only
  // sparsity reads out of bounds and produces garbage.
  // Build a 3x3 symmetric tridiag pattern via triplet (col-major fill).
  // Pattern (rows, cols): (0,0),(1,0),(0,1),(1,1),(2,1),(1,2),(2,2).
  const rows = [0n, 1n, 0n, 1n, 2n, 1n, 2n];
  const cols = [0n, 0n, 1n, 1n, 1n, 2n, 2n];
  const sp = M.Sparsity.triplet(3n, 3n, rows, cols);
  // A = [[2,1,0],[1,2,1],[0,1,2]], values in CSC order.
  const A = M.DM(sp, [2, 1, 1, 2, 1, 1, 2]);
  const ls = M.Linsol("ls", "tridiag", sp);
  ls.sfact(A);
  ls.nfact(A);
  // Pick b so that x = [1, 2, 3]: A*x = [2+2, 1+4+3, 2+6] = [4, 8, 8].
  const b = M.DM([4, 8, 8]);
  const x = ls.solve(A, b, false);
  assertArrayAlmostEqual(x.nonzeros(), [1, 2, 3], 6, "tridiag x");
}

// ---------------------------------------------------------------------
// integrator plugins
// ---------------------------------------------------------------------
function _integrator(M, plugin) {
  // dx/dt = x, x(0)=1 -> x(1) = e.
  const x   = M.SX.sym("x");
  const dae = { x: x, ode: x };
  const F   = M.integrator("F", plugin, dae);
  const r   = F({ x0: M.DM(1) });
  assertAlmost(Number(r.xf.nonzeros()[0]), Math.E, 3,
    plugin + ": x(1) ~ e");
}
function test_integrator_rk         (M) { _integrator(M, "rk"); }
function test_integrator_collocation(M) {
  // Collocation defaults to implicit Radau IIA; same trivial ODE.
  // Higher-order, so accuracy is much better than rk's default 4-step.
  const x   = M.SX.sym("x");
  const dae = { x: x, ode: x };
  const F   = M.integrator("F", "collocation", dae);
  const r   = F({ x0: M.DM(1) });
  assertAlmost(Number(r.xf.nonzeros()[0]), Math.E, 4,
    "collocation: x(1) ~ e");
}

// ---------------------------------------------------------------------
// interpolant plugins
// ---------------------------------------------------------------------
function test_interpolant_linear(M) {
  // 1-D linear interpolant of f(0)=0, f(1)=1, f(2)=4, f(3)=9.
  const grid   = [[0, 1, 2, 3]];
  const values = [0, 1, 4, 9];
  const f = M.interpolant("f", "linear", grid, values);
  // Mid-point of [1, 2] -> linear blend = (1+4)/2 = 2.5.
  const r1 = f(M.DM(1.5));
  assertAlmost(Number(r1.nonzeros()[0]), 2.5, 9, "linear mid 1.5");
  // Knot point: exact.
  const r2 = f(M.DM(2));
  assertAlmost(Number(r2.nonzeros()[0]), 4, 9, "linear at knot 2");
}
function test_interpolant_bspline(M) {
  // 1-D B-spline of the same 4-point sample.  Default degree 3.
  // Don't check exact mid-point match (bspline is smooth, not piecewise
  // linear); just verify the call succeeds and recovers knot values
  // approximately.  Bspline knots use clamped boundary conditions, so
  // f(0) = 0, f(3) = 9 exactly.
  const grid   = [[0, 1, 2, 3]];
  const values = [0, 1, 4, 9];
  const f = M.interpolant("f", "bspline", grid, values);
  const r0 = f(M.DM(0));
  assertAlmost(Number(r0.nonzeros()[0]), 0, 6, "bspline at left knot");
  const r3 = f(M.DM(3));
  assertAlmost(Number(r3.nonzeros()[0]), 9, 6, "bspline at right knot");
}

// ---------------------------------------------------------------------
// rootfinder plugins:  find x s.t. f(x) = 0
// ---------------------------------------------------------------------
function _rootfinder(M, plugin) {
  // Solve x^2 - 2 = 0 starting from x0 = 1.5 -> root ~ 1.41421356.
  const x = M.SX.sym("x");
  const f = M.Function("f", [x], [M.minus(M.times(x, x), M.SX(2))]);
  const opts = {};
  // implicit_to_nlp's "nlpsol" plugin requires a sub-NLP solver;
  // route through ipopt.
  if (plugin === "nlpsol") opts["nlpsol"] = "ipopt";
  const F = M.rootfinder("F", plugin, f, opts);
  const r = F(M.DM(1.5));
  assertAlmost(Number(r.nonzeros()[0]), Math.SQRT2, 6,
    plugin + ": sqrt(2)");
}
function test_rootfinder_newton     (M) { _rootfinder(M, "newton"); }
function test_rootfinder_fast_newton(M) { _rootfinder(M, "fast_newton"); }
function test_rootfinder_implicit_to_nlp(M) {
  // implicit_to_nlp wraps an Nlpsol; test only if ipopt is present.
  if (!M.has_nlpsol("ipopt")) {
    console.log("    (skip -- ipopt not loaded)");
    return;
  }
  _rootfinder(M, "nlpsol");
}

// ---------------------------------------------------------------------
// conic plugins:  min 0.5 x'Hx + g'x  subject to lba <= Ax <= uba
// ---------------------------------------------------------------------
function _conic(M, plugin, opts) {
  // Same simple QP used in conic.js's test_general_unconstrained.
  // min 0.5 (x^2 + y^2) - 0.7 x - 2.3 y -> [0.7, 2.3].
  const H = M.DM([[1, 0], [0, 1]]);
  const G = M.DM([[-0.7], [-2.3]]);
  const inf = Infinity;
  const LBX = M.DM([[-inf], [-inf]]);
  const UBX = M.DM([[inf], [inf]]);
  const qp = { h: H.sparsity() };
  const solver = M.conic("solver", plugin, qp, opts || {});
  const res = solver.call({ h: H, g: G, lbx: LBX, ubx: UBX });
  assertArrayAlmostEqual(res.x.nonzeros(), [0.7, 2.3], 4,
    plugin + ": x");
}
function test_conic_qrqp (M) { _conic(M, "qrqp",  { print_iter: false, print_header: false }); }
function test_conic_ipqp (M) { _conic(M, "ipqp",  { print_iter: false, print_header: false }); }
function test_conic_nlpsol(M) {
  // qp_to_nlp wrapper requires an Nlpsol; route through ipopt.
  if (!M.has_nlpsol("ipopt")) {
    console.log("    (skip -- ipopt not loaded)");
    return;
  }
  _conic(M, "nlpsol", { nlpsol: "ipopt", nlpsol_options: { ipopt: { print_level: 0 }, print_time: false } });
}
function test_conic_highs(M) { _conic(M, "highs", { highs: { output_flag: false } }); }

// ---------------------------------------------------------------------
// nlpsol plugins:  min f(x) s.t. lbg <= g(x) <= ubg
// ---------------------------------------------------------------------
function _nlpsol(M, plugin, opts) {
  // Rosenbrock-like simple convex bowl: f(x,y) = (x-1)^2 + (y-2)^2.
  // Optimum at (1, 2) with f* = 0.  No constraints.
  const xy = M.SX.sym("xy", 2);
  const x  = xy.get(false, M.Slice(0n, false));
  const y  = xy.get(false, M.Slice(1n, false));
  const f  = M.plus(M.times(M.minus(x, M.SX(1)), M.minus(x, M.SX(1))),
                    M.times(M.minus(y, M.SX(2)), M.minus(y, M.SX(2))));
  const nlp = { x: xy, f: f };
  const solver = M.nlpsol("S", plugin, nlp, opts || {});
  const res = solver.call({ x0: M.DM([0, 0]) });
  assertArrayAlmostEqual(res.x.nonzeros(), [1, 2], 3, plugin + ": x");
  assertAlmost(Number(res.f.nonzeros()[0]), 0, 4, plugin + ": f");
}
function test_nlpsol_sqpmethod(M) {
  // SQP needs a QP sub-solver; route through qrqp.
  _nlpsol(M, "sqpmethod", {
    qpsol: "qrqp",
    qpsol_options: { print_iter: false, print_header: false },
    print_iteration: false, print_header: false, print_status: false, print_time: false
  });
}
function test_nlpsol_qrsqp(M) {
  _nlpsol(M, "qrsqp", { print_iteration: false, print_header: false, print_time: false });
}
function test_nlpsol_feasiblesqpmethod(M) {
  _nlpsol(M, "feasiblesqpmethod", {
    qpsol: "qrqp",
    qpsol_options: { print_iter: false, print_header: false },
    print_iteration: false, print_header: false, print_status: false, print_time: false
  });
}
function test_nlpsol_ipopt(M) {
  _nlpsol(M, "ipopt", { ipopt: { print_level: 0, sb: "yes" }, print_time: false });
}
function test_nlpsol_scpgen(M) {
  // scpgen needs MXFunction (calls Function::generate_lifted internally,
  // which is MXFunction-only).  Build the NLP from MX symbolics so the
  // emitted Function is an MXFunction.  Same Rosenbrock-style bowl
  // f = (x-1)^2 + (y-2)^2 as _nlpsol, optimum at (1,2).
  const xy = M.MX.sym("xy", 2);
  const x  = xy.get(false, M.Slice(0n, false));
  const y  = xy.get(false, M.Slice(1n, false));
  const f  = M.plus(M.times(M.minus(x, M.MX(1)), M.minus(x, M.MX(1))),
                    M.times(M.minus(y, M.MX(2)), M.minus(y, M.MX(2))));
  const nlp = { x: xy, f: f };
  const solver = M.nlpsol("S", "scpgen", nlp, {
    qpsol: "qrqp",
    qpsol_options: { print_iter: false, print_header: false },
    print_header: false, print_time: false
  });
  const res = solver.call({ x0: M.DM([0, 0]) });
  assertArrayAlmostEqual(res.x.nonzeros(), [1, 2], 3, "scpgen: x");
  assertAlmost(Number(res.f.nonzeros()[0]), 0, 4, "scpgen: f");
}

runTests([
  ["test_linsol_qr",                test_linsol_qr],
  ["test_linsol_ldl",               test_linsol_ldl],
  ["test_linsol_lsqr",              test_linsol_lsqr],
  ["test_linsol_symbolicqr",        test_linsol_symbolicqr],
  ["test_linsol_tridiag",           test_linsol_tridiag],
  ["test_integrator_rk",            test_integrator_rk],
  ["test_integrator_collocation",   test_integrator_collocation],
  ["test_interpolant_linear",       test_interpolant_linear],
  ["test_interpolant_bspline",      test_interpolant_bspline],
  ["test_rootfinder_newton",        test_rootfinder_newton],
  ["test_rootfinder_fast_newton",   test_rootfinder_fast_newton],
  ["test_rootfinder_implicit_to_nlp", test_rootfinder_implicit_to_nlp],
  ["test_conic_qrqp",               test_conic_qrqp],
  ["test_conic_ipqp",               test_conic_ipqp],
  ["test_conic_nlpsol",             test_conic_nlpsol],
  ["test_conic_highs",              test_conic_highs],
  ["test_nlpsol_sqpmethod",         test_nlpsol_sqpmethod],
  ["test_nlpsol_qrsqp",             test_nlpsol_qrsqp],
  ["test_nlpsol_feasiblesqpmethod", test_nlpsol_feasiblesqpmethod],
  ["test_nlpsol_ipopt",             test_nlpsol_ipopt],
  ["test_nlpsol_scpgen",            test_nlpsol_scpgen],
]);

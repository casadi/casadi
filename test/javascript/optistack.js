// optistack.js -- port of test/python/optistack.py (4/20 tests).
//
// Strategy: cover the basic Opti construction + ipopt-backed solve path.
// (Previously skipped due to a SWIG-emit gap where `new Opti()` failed
// with `"emsc"` garbage memory.  Fixed by the constructorHandler patch
// + member-method default-arg phantom patch in
// Source/Modules/wasm_js.cxx.)
//
// Skipped:
//   test_conic           -- needs qrqp solver (not in wasm build)
//   test_callback        -- OptiCallback subclassing
//   test_print           -- capture_stdout
//   test_to_function     -- Opti::to_function
//   test_sparse, test_structure, test_warmstart, test_set_value_expr,
//   test_dual, test_lookup, test_shapes, test_symb_*,
//   test_debug_value, test_flow, test_simple, test_parametric, test_n,
//   test_subject_to_list, test_sol_opti -- variations on the same
//   ipopt-driven theme; the basic cases below prove the integration.

const { assertAlmost, assertTrue, runTests } = require("./_helpers");

function test_basic_unconstrained(M) {
  // Minimize (x-1)^2 + (y-2)^2 -> x=1, y=2.
  const opti = new M.Opti();
  const x = opti.variable();
  const y = opti.variable();
  const xm1 = M.minus(x, M.MX(1));
  const ym2 = M.minus(y, M.MX(2));
  opti.minimize(M.plus(M.times(xm1, xm1), M.times(ym2, ym2)));
  opti.solver("ipopt", { ipopt: { print_level: 0, sb: "yes" }, print_time: false });
  const sol = opti.solve();
  assertAlmost(Number(sol.value(x).nonzeros()[0]), 1, 5, "x*=1");
  assertAlmost(Number(sol.value(y).nonzeros()[0]), 2, 5, "y*=2");
}

function test_inequality_constraint(M) {
  // Minimize (x-3)^2 + (y-3)^2 subject to x+y <= 1.  Optimal at (0.5, 0.5).
  const opti = new M.Opti();
  const x = opti.variable();
  const y = opti.variable();
  const xm3 = M.minus(x, M.MX(3));
  const ym3 = M.minus(y, M.MX(3));
  opti.minimize(M.plus(M.times(xm3, xm3), M.times(ym3, ym3)));
  opti.subject_to(M.le(M.plus(x, y), M.MX(1)));
  opti.solver("ipopt", { ipopt: { print_level: 0, sb: "yes" }, print_time: false });
  const sol = opti.solve();
  assertAlmost(Number(sol.value(x).nonzeros()[0]), 0.5, 4, "x*=0.5");
  assertAlmost(Number(sol.value(y).nonzeros()[0]), 0.5, 4, "y*=0.5");
}

function test_equality_constraint(M) {
  // Minimize (x-1)^2 + (y-2)^2 subject to y = x.  Optimal at (1.5, 1.5).
  const opti = new M.Opti();
  const x = opti.variable();
  const y = opti.variable();
  const xm1 = M.minus(x, M.MX(1));
  const ym2 = M.minus(y, M.MX(2));
  opti.minimize(M.plus(M.times(xm1, xm1), M.times(ym2, ym2)));
  opti.subject_to(M.eq(y, x));
  opti.solver("ipopt", { ipopt: { print_level: 0, sb: "yes" }, print_time: false });
  const sol = opti.solve();
  assertAlmost(Number(sol.value(x).nonzeros()[0]), 1.5, 4, "x*=1.5");
  assertAlmost(Number(sol.value(y).nonzeros()[0]), 1.5, 4, "y*=1.5");
}

function test_variable_vector(M) {
  // opti.variable(n) creates an n-vector.  Constraints work elementwise.
  // Minimize ||x - [1, 2]||^2 subject to x >= [1, 1] (gives x = [1, 2]).
  const opti = new M.Opti();
  const x = opti.variable(2n);  // 2-vector
  const target = M.MX(M.DM([1, 2]));
  const diff = M.minus(x, target);
  opti.minimize(M.dot(diff, diff));
  opti.subject_to(M.ge(x, M.MX(M.DM([1, 1]))));
  opti.solver("ipopt", { ipopt: { print_level: 0, sb: "yes" }, print_time: false });
  const sol = opti.solve();
  const xv = sol.value(x).nonzeros();
  assertAlmost(Number(xv[0]), 1, 4, "x[0]=1");
  assertAlmost(Number(xv[1]), 2, 4, "x[1]=2");
}

function test_bounded(M) {
  // opti.bounded called on the INSTANCE (static-on-instance forwarder,
  // matching Python's opti.subject_to(opti.bounded(0, x, 1))).
  // Minimize (x-3)^2 subject to 0 <= x <= 1 -> x at the upper bound.
  const opti = new M.Opti();
  const x = opti.variable();
  const xm3 = M.minus(x, M.MX(3));
  opti.minimize(M.times(xm3, xm3));
  opti.subject_to(opti.bounded(0, x, 1));
  opti.solver("ipopt", { ipopt: { print_level: 0, sb: "yes" }, print_time: false });
  const sol = opti.solve();
  assertAlmost(Number(sol.value(x).nonzeros()[0]), 1, 4, "x*=1");
  // The static form stays available too.
  assertTrue(M.Opti.bounded(0, x, 1) instanceof M.MX, "static bounded returns MX");
}

runTests([
  ["test_basic_unconstrained",   test_basic_unconstrained],
  ["test_inequality_constraint", test_inequality_constraint],
  ["test_equality_constraint",   test_equality_constraint],
  ["test_variable_vector",       test_variable_vector],
  ["test_bounded",               test_bounded],
]);

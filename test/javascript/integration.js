// integration.js -- partial port of test/python/integration.py.
//
// Strategy: only the `rk` integrator plugin is built into the wasm-js
// build (no cvodes/idas).  Test the basic ODE solve path via the 3-arg
// `integrator(name, plugin, dae)` form which uses default t0/tf.
//
// Skipped (most of integration.py's 30+ tests):
//   - Anything using cvodes/idas/collocation (not in wasm build)
//   - Sensitivity-analysis tests (use checkfunction helper)
//   - DAE/sundials-specific tests
//   - test_jacFwd/Adj, test_simulator, test_eular
//
// The 5-arg form `integrator(name, plugin, dae, t0, tf)` currently
// fails due to the global-function dispatcher's first-wins-per-arity
// trap (5-arg has multiple overloads: tf as scalar vs. tf as
// vector<double> grid).  Use the 3-arg form for now.

const { assertAlmost, assertArrayAlmostEqual, assertEqual, assertTrue, runTests } = require("./_helpers");

function test_rk_simple_ode(M) {
  // dx/dt = x, x(0)=1 -> x(1) = e ~= 2.71828.
  const x = M.SX.sym("x");
  const dae = { x: x, ode: x };
  const F = M.integrator("F1", "rk", dae);
  // 3-arg form: default t0=0, tf=1.
  const r = F({ x0: M.DM(1) });
  // r.xf is the state at tf.
  assertAlmost(Number(r.xf.nonzeros()[0]), Math.E, 3, "x(1) ≈ e");
}

function test_rk_with_param(M) {
  // dx/dt = -k*x, x(0)=1, k=0.5 -> x(1) = exp(-0.5).
  const x = M.SX.sym("x");
  const k = M.SX.sym("k");
  const dae = { x: x, p: k, ode: M.times(M.SX(-1), M.times(k, x)) };
  const F = M.integrator("F2", "rk", dae);
  const r = F({ x0: M.DM(1), p: M.DM(0.5) });
  assertAlmost(Number(r.xf.nonzeros()[0]), Math.exp(-0.5), 3, "x(1) = exp(-0.5)");
}

function test_rk_2state(M) {
  // 2-state harmonic oscillator: x1' = x2, x2' = -x1.  At t=1 from
  // x(0)=[1,0]: x1(1) = cos(1), x2(1) = -sin(1).
  const x  = M.SX.sym("x", 2);
  const x1 = x.get(false, M.Slice(0n, false));
  const x2 = x.get(false, M.Slice(1n, false));
  const ode = M.vcat([x2, M.times(M.SX(-1), x1)]);
  const dae = { x: x, ode: ode };
  const F = M.integrator("F3", "rk", dae);
  const r = F({ x0: M.DM([1, 0]) });
  const xf = r.xf.nonzeros();
  assertAlmost(Number(xf[0]),  Math.cos(1), 3, "x1(1) = cos(1)");
  assertAlmost(Number(xf[1]), -Math.sin(1), 3, "x2(1) = -sin(1)");
}

function test_rk_5arg_scalar_tf(M) {
  // 5-arg form integrator(name, plugin, dae, t0, tf): tf is scalar.
  // Was previously broken by the global-function dispatcher first-wins-
  // per-arity trap (multiple 5-arg overloads with different tf types).
  // Now fixed via the type_checks discriminator in
  // wasm_js.cxx::globalfunctionHandler.
  const x = M.SX.sym("x");
  const F = M.integrator("Fs5", "rk", { x: x, ode: x }, 0, 1);
  const r = F({ x0: M.DM(1) });
  assertAlmost(Number(r.xf.nonzeros()[0]), Math.E, 3, "x(1) ≈ e (5-arg form)");
}

function test_rk_5arg_grid_tf(M) {
  // 5-arg form with tf as vector<double> grid -- different overload
  // than scalar tf.  Exercises the type_checks discrimination.
  const x = M.SX.sym("x");
  const F = M.integrator("Fg5", "rk", { x: x, ode: x }, 0, [0.5, 1.0]);
  const r = F({ x0: M.DM(1) });
  // xf is a 1x2 matrix [exp(0.5), exp(1)].
  const vals = r.xf.nonzeros();
  assertAlmost(Number(vals[0]), Math.exp(0.5), 3, "x(0.5) ≈ e^0.5");
  assertAlmost(Number(vals[1]), Math.E,        3, "x(1) ≈ e");
}

function test_integrator_inputs_outputs(M) {
  // Integrator-as-Function exposes standard inputs/outputs.
  const x = M.SX.sym("x");
  const F = M.integrator("F4", "rk", { x: x, ode: x });
  assertTrue(Number(F.n_in()) >= 1, "n_in >= 1 (at least x0)");
  assertTrue(Number(F.n_out()) >= 1, "n_out >= 1 (at least xf)");
  // Standard input names include "x0", "p", etc.
  const names = F.name_in();
  assertTrue(names.includes("x0"), "n_in includes 'x0'");
}

runTests([
  ["test_rk_simple_ode",         test_rk_simple_ode],
  ["test_rk_with_param",         test_rk_with_param],
  ["test_rk_2state",             test_rk_2state],
  ["test_rk_5arg_scalar_tf",     test_rk_5arg_scalar_tf],
  ["test_rk_5arg_grid_tf",       test_rk_5arg_grid_tf],
  ["test_integrator_inputs_outputs", test_integrator_inputs_outputs],
]);

// rootfinder.js -- subset port of test/python/implicitfunction.py.
//
// Exercises the rootfinder() factory against every rootfinder plugin
// loaded in the wasm-js build: newton, fast_newton, nlpsol (implicit_to_nlp),
// kinsol.
//
// Coverage:
//   test_scalar1        -- sin(x) = 0, x0 = 6                    -> 2*pi
//   test_scalar2        -- x - arcsin(y) = 0 with param x       -> y = sin(x)
//   test_sqrt2          -- x^2 - 2 = 0                            -> sqrt(2)
//   test_linear_2x2     -- Ax - b = 0 for fixed A, b              -> A\b
//   test_linear_with_par -- Ax - b = 0 with A, b as params       -> A\b
//   test_cubic          -- x^3 - x - 1 = 0 starting x0=1.5        -> 1.32472...
//   test_with_constraint -- exp(x) - 1 = 0, x0 = 1                -> 0
//   test_2d_system      -- 2-D system x^2+y^2-2=0, x-y=0          -> (1, 1)
//   test_multi_eval     -- same solver called twice with different x0
//   test_indirect       -- solver embedded in a Function, eval at MX level
//
// Each test runs against each loaded plugin.

const { assertAlmost, assertArrayAlmostEqual,
        modulePath } = require("./_helpers");

const fs = require("fs");
const path = require("path");

function rootfinderOpts(plugin) {
  if (plugin === "nlpsol") return { nlpsol: "ipopt",
    nlpsol_options: { print_time: false, ipopt: { print_level: 0 } } };
  if (plugin === "newton")      return { print_iteration: false };
  if (plugin === "fast_newton") return {};
  if (plugin === "kinsol")      return { print_level: 0 };
  return {};
}

// ----- tests -----

function test_scalar1(M, plugin) {
  // sin(x) = 0; x0=6 closest root is 2*pi.
  const x = M.SX.sym("x");
  const f = M.Function("f", [x], [M.sin(x)]);
  const F = M.rootfinder("F", plugin, f, rootfinderOpts(plugin));
  const r = F(M.DM(6));
  assertAlmost(Number(r.nonzeros()[0]), 2 * Math.PI, 6,
    `${plugin} scalar1`);
}

function test_scalar2(M, plugin) {
  // x - arcsin(y) = 0  (residual in y, parameter in x):
  //   -> y* = sin(x).  For x = 0.2 -> y = sin(0.2) ~ 0.19867
  const y = M.SX.sym("y");
  const x = M.SX.sym("x");
  const f = M.Function("f", [y, x],
    [M.minus(x, M.asin(y))]);
  const F = M.rootfinder("F", plugin, f, rootfinderOpts(plugin));
  // First arg is initial guess for y, then any extra params (x).
  const r = F(M.DM(0), M.DM(0.2));
  assertAlmost(Number(r.nonzeros()[0]), Math.sin(0.2), 6,
    `${plugin} scalar2`);
}

function test_sqrt2(M, plugin) {
  // x^2 - 2 = 0 with x0 = 1.5
  const x = M.SX.sym("x");
  const f = M.Function("f", [x],
    [M.minus(M.times(x, x), M.SX(2))]);
  const F = M.rootfinder("F", plugin, f, rootfinderOpts(plugin));
  const r = F(M.DM(1.5));
  assertAlmost(Number(r.nonzeros()[0]), Math.SQRT2, 6,
    `${plugin} sqrt2`);
}

function test_linear_2x2(M, plugin) {
  // Solve A*x - b = 0 with A = [[1,2],[3,2.1]], b = [0.7,0.6]
  // Built as separate scalars to dodge DM/SX-mixing mtimes overloads.
  const x0 = M.SX.sym("x0");
  const x1 = M.SX.sym("x1");
  const x = M.vertcat(x0, x1);
  // res = A*x - b, written out scalar by scalar:
  const r0 = M.minus(
    M.plus(M.times(M.SX(1), x0), M.times(M.SX(2), x1)), M.SX(0.7));
  const r1 = M.minus(
    M.plus(M.times(M.SX(3), x0), M.times(M.SX(2.1), x1)), M.SX(0.6));
  const f = M.Function("f", [x], [M.vertcat(r0, r1)]);
  const F = M.rootfinder("F", plugin, f, rootfinderOpts(plugin));
  const r = F(M.DM([0, 0]));
  // Analytical: det = 1*2.1 - 2*3 = -3.9
  assertArrayAlmostEqual(r.nonzeros(),
    [(2.1 * 0.7 - 2 * 0.6) / -3.9,
     (1 * 0.6 - 3 * 0.7) / -3.9], 6, `${plugin} linear_2x2`);
}

function test_linear_with_par(M, plugin) {
  // Solve A*x - b = 0 with A, b as PARAMETERS.
  const x = M.SX.sym("x", 2);
  const A = M.SX.sym("A", 2, 2);
  const b = M.SX.sym("b", 2);
  const res = M.minus(M.mtimes(A, x), b);
  const f = M.Function("f", [x, A, b], [res]);
  const F = M.rootfinder("F", plugin, f, rootfinderOpts(plugin));
  const A_ = M.DM([[1, 2], [3, 2.1]]);
  const b_ = M.DM([0.7, 0.6]);
  const r = F(M.DM([0, 0]), A_, b_);
  assertArrayAlmostEqual(r.nonzeros(),
    [(2.1 * 0.7 - 2 * 0.6) / -3.9,
     (1 * 0.6 - 3 * 0.7) / -3.9], 6, `${plugin} linear_with_par`);
}

function test_cubic(M, plugin) {
  // x^3 - x - 1 = 0; real root ~ 1.32471795724475
  const x = M.SX.sym("x");
  const x3 = M.times(M.times(x, x), x);
  const f = M.Function("f", [x],
    [M.minus(M.minus(x3, x), M.SX(1))]);
  const F = M.rootfinder("F", plugin, f, rootfinderOpts(plugin));
  const r = F(M.DM(1.5));
  assertAlmost(Number(r.nonzeros()[0]), 1.32471795724474602596, 5,
    `${plugin} cubic`);
}

function test_with_constraint(M, plugin) {
  // exp(x) - 1 = 0  -> x = 0;  test convergence from x0 = 1
  // kinsol with constraints would gate x>=0, but we skip that here.
  const x = M.SX.sym("x");
  const f = M.Function("f", [x],
    [M.minus(M.exp(x), M.SX(1))]);
  const F = M.rootfinder("F", plugin, f, rootfinderOpts(plugin));
  const r = F(M.DM(0.5));
  assertAlmost(Number(r.nonzeros()[0]), 0, 5,
    `${plugin} exp(x)-1`);
}

function test_2d_system(M, plugin) {
  // f0 = x^2 + y^2 - 2 = 0
  // f1 = x - y = 0
  // -> x = y = 1
  const xv = M.SX.sym("x");
  const yv = M.SX.sym("y");
  const xy = M.vertcat(xv, yv);
  const f0 = M.minus(
    M.plus(M.times(xv, xv), M.times(yv, yv)), M.SX(2));
  const f1 = M.minus(xv, yv);
  const f = M.Function("f", [xy], [M.vertcat(f0, f1)]);
  const F = M.rootfinder("F", plugin, f, rootfinderOpts(plugin));
  const r = F(M.DM([0.5, 0.5]));
  assertArrayAlmostEqual(r.nonzeros(), [1, 1], 5, `${plugin} 2d_system`);
}

function test_multi_eval(M, plugin) {
  // Same rootfinder called twice with different starting points
  // (multimodal: f = sin(x)).  x0 = 1 -> 0, x0 = 4 -> pi.
  const x = M.SX.sym("x");
  const f = M.Function("f", [x], [M.sin(x)]);
  const F = M.rootfinder("F", plugin, f, rootfinderOpts(plugin));
  const r1 = F(M.DM(1));
  assertAlmost(Number(r1.nonzeros()[0]), 0, 5,
    `${plugin} multi_eval 1`);
  const r2 = F(M.DM(4));
  assertAlmost(Number(r2.nonzeros()[0]), Math.PI, 5,
    `${plugin} multi_eval pi`);
}

function test_indirect(M, plugin) {
  // Wrap rootfinder in another Function and evaluate at MX level.
  const y = M.SX.sym("y");
  const x = M.SX.sym("x");
  const f = M.Function("f", [y, x],
    [M.minus(x, M.asin(y))]);
  const F = M.rootfinder("F", plugin, f, rootfinderOpts(plugin));
  // Build trial = Function(X) -> R = F(0, X)
  const X = M.MX.sym("X");
  const R = F(M.MX(0), X);
  const trial = M.Function("trial", [X], [R]);
  const out = trial.call([M.DM(0.2)])[0];
  assertAlmost(Number(out.nonzeros()[0]), Math.sin(0.2), 6,
    `${plugin} indirect`);
}

// ----- driver -----

const ALL_PLUGINS = ["newton", "fast_newton", "nlpsol", "kinsol"];

const tests = [
  ["test_scalar1",          test_scalar1],
  ["test_scalar2",          test_scalar2],
  ["test_sqrt2",            test_sqrt2],
  ["test_linear_2x2",       test_linear_2x2],
  ["test_linear_with_par",  test_linear_with_par],
  ["test_cubic",            test_cubic],
  ["test_with_constraint",  test_with_constraint],
  ["test_2d_system",        test_2d_system],
  ["test_multi_eval",       test_multi_eval],
  ["test_indirect",         test_indirect],
];

(async () => {
  if (!fs.existsSync(modulePath)) {
    console.error(`SKIP: ${modulePath} not built.`);
    process.exit(2);
  }
  const create = require(modulePath);
  const M = await create();
  const plugins = ALL_PLUGINS.filter((p) => M.has_rootfinder(p));
  if (plugins.length === 0) {
    console.error("SKIP: no rootfinder plugin loaded.");
    process.exit(2);
  }
  let pass = 0, fail = 0;
  const failures = [];
  for (const plugin of plugins) {
    for (const [name, fn] of tests) {
      const label = `[${plugin}] ${name}`;
      try { fn(M, plugin); console.log(`ok   -- ${label}`); pass++; }
      catch (e) {
        console.log(`FAIL -- ${label}: ${e.message}`);
        fail++; failures.push([label, e.message]);
      }
    }
  }
  console.log(`\n${pass}/${pass + fail} passed`);
  if (failures.length > 0) {
    console.log("\n--- failures ---");
    for (const [n, m] of failures) console.log(`  ${n}: ${m}`);
  }
  process.exit(fail === 0 ? 0 : 1);
})().catch((e) => { console.error("FATAL:", e); process.exit(1); });

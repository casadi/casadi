// nlp.js -- subset port of test/python/nlp.py.
//
// Exercises nlpsol() against ipopt + sqpmethod (the two NLP solvers
// preloaded as SIDE_MODULEs in the wasm-js build).  ipopt is the
// primary target; sqpmethod is a secondary smoke (smaller solver
// surface, doesn't need the Fortran/MUMPS cascade).
//
// Ports a representative subset of nlp.py's IPOPT family:
//   test_basic           -- (x-1)^2, bounds [-10,10], linear g=x       -> x=1, f=0
//   test_basic_par       -- parametric (x-p)^2, p=1                    -> x=1
//   test_inf_bounds      -- infinity bounds                            -> x=1
//   test_rosenbrock      -- 2-D Rosenbrock, no g                       -> x=[1,1]
//   test_rosenbrock_g    -- Rosenbrock with linear g=x+y, lbg=ubg=2    -> x=[1,1]
//   test_rosenbrock_fix  -- Rosenbrock with y fixed via lbx=ubx=1      -> x=[1,1]
//   test_initial_cond    -- x0 effect on multi-modal -cos(x)
//   test_warmstart       -- two-call warmstart (basic)
//   test_x0_default      -- x0 omitted -> zero default
//   test_lbg_default     -- lbg omitted -> +/- inf default
//
// Driver runs each test for each available solver.

const path = require("path");
const fs = require("fs");

const wasmDir = path.resolve(__dirname, "../../build-wasm/swig/wasm-js");
const modulePath = path.join(wasmDir, "casadi.js");
if (!fs.existsSync(modulePath)) {
  console.error(`SKIP: ${modulePath} not built.`);
  process.exit(2);
}

function assertAlmost(actual, expected, places, msg) {
  if (!isFinite(actual) || Math.abs(actual - expected) > Math.pow(10, -places)) {
    throw new Error(`${msg}: expected ${expected}, got ${actual}`);
  }
}
function assertArrayAlmostEqual(actual, expected, places, msg) {
  if (actual.length !== expected.length) {
    throw new Error(`${msg}: length expected ${expected.length}, got ${actual.length}`);
  }
  for (let i = 0; i < actual.length; ++i) {
    assertAlmost(Number(actual[i]), expected[i], places, `${msg}[${i}]`);
  }
}

// ----- helpers -----

// Per-solver quiet-mode options.  All zero-print for clean output.
function quietOpts(plugin) {
  if (plugin === "ipopt") {
    return { ipopt: { print_level: 0, sb: "yes" }, print_time: false };
  }
  // sqpmethod / qrsqp: only `print_header` + `print_time` are common.
  return { print_header: false, print_time: false };
}

// Build a fresh nlp dict from an SX-symbolic recipe.  Keeps the test
// bodies focused on numerical assertions.
function makeNlp(M, recipe) {
  const SX = M.SX;
  return recipe(SX);
}

// ----- tests -----

function test_basic(M, plugin) {
  // (x-1)^2, lbx=-10, ubx=10, lbg=-10, ubg=10
  const nlp = makeNlp(M, (SX) => {
    const x = SX.sym("x");
    const xm1 = M.minus(x, SX(1));
    return { x: x, f: M.times(xm1, xm1), g: x };
  });
  const solver = M.nlpsol("S", plugin, nlp, quietOpts(plugin));
  const r = solver.call({ lbx: -10, ubx: 10, lbg: -10, ubg: 10 });
  assertAlmost(Number(r.f.nonzeros()[0]), 0, 8, `${plugin}: f`);
  assertAlmost(Number(r.x.nonzeros()[0]), 1, 6, `${plugin}: x`);
  assertAlmost(Number(r.g.nonzeros()[0]), 1, 6, `${plugin}: g`);
}

function test_basic_par(M, plugin) {
  // (x-p)^2 with p=1 -> x=1
  const nlp = makeNlp(M, (SX) => {
    const x = SX.sym("x");
    const p = SX.sym("p");
    const xmp = M.minus(x, p);
    return { x: x, p: p, f: M.times(xmp, xmp), g: x };
  });
  const solver = M.nlpsol("S", plugin, nlp, quietOpts(plugin));
  const r = solver.call({ lbx: -10, ubx: 10, lbg: -10, ubg: 10, p: 1 });
  assertAlmost(Number(r.f.nonzeros()[0]), 0, 8, `${plugin}: f`);
  assertAlmost(Number(r.x.nonzeros()[0]), 1, 6, `${plugin}: x`);
}

function test_inf_bounds(M, plugin) {
  // Same as basic but with +/-inf bounds.  Tests that the wasm boundary
  // marshalling preserves JS Infinity through to the C++ side.
  const nlp = makeNlp(M, (SX) => {
    const x = SX.sym("x");
    const xm1 = M.minus(x, SX(1));
    return { x: x, f: M.times(xm1, xm1), g: x };
  });
  const inf = Infinity;
  const solver = M.nlpsol("S", plugin, nlp, quietOpts(plugin));
  const r = solver.call({ lbx: -inf, ubx: inf, lbg: -inf, ubg: inf });
  assertAlmost(Number(r.f.nonzeros()[0]), 0, 8, `${plugin}: f`);
  assertAlmost(Number(r.x.nonzeros()[0]), 1, 6, `${plugin}: x`);
}

function test_rosenbrock(M, plugin) {
  // f = (1-x)^2 + 100 (y - x^2)^2; optimum (1,1), f*=0.  No g.
  const nlp = makeNlp(M, (SX) => {
    const x = SX.sym("x");
    const y = SX.sym("y");
    const xy = M.vertcat(x, y);
    const t1 = M.minus(SX(1), x);
    const t2 = M.minus(y, M.times(x, x));
    const f = M.plus(M.times(t1, t1), M.times(SX(100), M.times(t2, t2)));
    return { x: xy, f: f };
  });
  const solver = M.nlpsol("S", plugin, nlp, quietOpts(plugin));
  const r = solver.call({ lbx: [-10, -10], ubx: [10, 10] });
  assertAlmost(Number(r.f.nonzeros()[0]), 0, 6, `${plugin}: f`);
  assertArrayAlmostEqual(r.x.nonzeros(), [1, 1], 4, `${plugin}: x`);
}

function test_rosenbrock_g(M, plugin) {
  // Rosenbrock with linear g=x+y, lbg=ubg=2.  Constrained solution
  // happens to coincide with unconstrained (1+1=2 is exactly active).
  const nlp = makeNlp(M, (SX) => {
    const x = SX.sym("x");
    const y = SX.sym("y");
    const xy = M.vertcat(x, y);
    const t1 = M.minus(SX(1), x);
    const t2 = M.minus(y, M.times(x, x));
    const f = M.plus(M.times(t1, t1), M.times(SX(100), M.times(t2, t2)));
    return { x: xy, f: f, g: M.plus(x, y) };
  });
  const solver = M.nlpsol("S", plugin, nlp, quietOpts(plugin));
  const r = solver.call({
    lbx: [-10, -10], ubx: [10, 10], lbg: [2], ubg: [2]
  });
  assertAlmost(Number(r.f.nonzeros()[0]), 0, 5, `${plugin}: f`);
  assertArrayAlmostEqual(r.x.nonzeros(), [1, 1], 4, `${plugin}: x`);
}

function test_rosenbrock_fix(M, plugin) {
  // Rosenbrock with y fixed to 1 via lbx[1]=ubx[1]=1.  Optimum
  // x[0]=1, y=1, f=0.  Tests bound-equality handling.
  const nlp = makeNlp(M, (SX) => {
    const x = SX.sym("x");
    const y = SX.sym("y");
    const xy = M.vertcat(x, y);
    const t1 = M.minus(SX(1), x);
    const t2 = M.minus(y, M.times(x, x));
    const f = M.plus(M.times(t1, t1), M.times(SX(100), M.times(t2, t2)));
    return { x: xy, f: f };
  });
  const solver = M.nlpsol("S", plugin, nlp, quietOpts(plugin));
  const r = solver.call({
    x0: [0, 1], lbx: [-10, 1], ubx: [10, 1]
  });
  assertAlmost(Number(r.f.nonzeros()[0]), 0, 6, `${plugin}: f`);
  assertArrayAlmostEqual(r.x.nonzeros(), [1, 1], 4, `${plugin}: x`);
}

function test_initial_cond(M, plugin) {
  // f = -cos(x); two distant initial guesses converge to two distinct
  // local minima (2*pi*k).  Tests x0 actually controls the search.
  // Use loose constraints so both runs are feasible.
  const nlp = makeNlp(M, (SX) => {
    const x = SX.sym("x");
    return { x: x, f: M.times(SX(-1), M.cos(x)), g: x };
  });
  const solver = M.nlpsol("S", plugin, nlp, quietOpts(plugin));
  const inf = Infinity;
  // Near x=0 -> converges to 0.
  const r1 = solver.call({ x0: [0.01], lbx: -inf, ubx: inf, lbg: -100, ubg: 100 });
  // Near x=2*pi -> converges to 2*pi (sqpmethod can be off-by-2pi
  // depending on regularisation; check it's AT SOME local min of -cos
  // i.e. a multiple of 2*pi).
  const r2 = solver.call({ x0: [6 * Math.PI + 0.01], lbx: -inf, ubx: inf, lbg: -100, ubg: 100 });
  // -cos at any 2*pi*k is -1.
  assertAlmost(Number(r1.f.nonzeros()[0]), -1, 5, `${plugin} r1: f`);
  assertAlmost(Number(r2.f.nonzeros()[0]), -1, 5, `${plugin} r2: f`);
  // r1.x should be near 0; r2.x should be far from 0.
  const x1 = Number(r1.x.nonzeros()[0]);
  const x2 = Number(r2.x.nonzeros()[0]);
  if (Math.abs(x1) > 0.5) {
    throw new Error(`${plugin}: r1 x expected near 0, got ${x1}`);
  }
  if (Math.abs(x2) < 1) {
    throw new Error(`${plugin}: r2 x expected far from 0, got ${x2}`);
  }
}

function test_warmstart(M, plugin) {
  // ipopt-only: solve once, then re-solve from the previous solution
  // with warm-start options enabled.  Just smoke-checks that the
  // warmstart options round-trip through Function.call.
  if (plugin !== "ipopt") return; // sqpmethod doesn't support warmstart
  const nlp = makeNlp(M, (SX) => {
    const x = SX.sym("x");
    const y = SX.sym("y");
    const xy = M.vertcat(x, y);
    const t1 = M.minus(SX(1), x);
    const t2 = M.minus(y, M.times(x, x));
    const f = M.plus(M.times(t1, t1), M.times(SX(100), M.times(t2, t2)));
    return { x: xy, f: f, g: M.plus(M.times(x, x), M.times(y, y)) };
  });
  const solver1 = M.nlpsol("S1", plugin, nlp, quietOpts(plugin));
  const r1 = solver1.call({
    x0: [0.5, 0.5], lbx: [-10, -10], ubx: [10, 10], lbg: [0], ubg: [1]
  });
  // Use r1.x as x0 for the warm-started second solve.
  const opts2 = {
    ipopt: {
      print_level: 0, sb: "yes",
      warm_start_init_point: "yes",
      warm_start_bound_push: 1e-6,
      warm_start_slack_bound_push: 1e-6,
      warm_start_mult_bound_push: 1e-6,
      mu_init: 1e-6,
    },
    print_time: false,
  };
  const solver2 = M.nlpsol("S2", plugin, nlp, opts2);
  const r2 = solver2.call({
    x0: r1.x, lam_g0: r1.lam_g, lam_x0: r1.lam_x,
    lbx: [-10, -10], ubx: [10, 10], lbg: [0], ubg: [1]
  });
  // Both should converge to the same minimum.
  assertAlmost(Number(r2.f.nonzeros()[0]),
               Number(r1.f.nonzeros()[0]), 5,
               `${plugin}: warmstart f matches cold`);
}

function test_x0_default(M, plugin) {
  // Omitting x0 should default to zeros.  For f=(x-1)^2 starting from
  // x=0 we still converge to x=1.
  const nlp = makeNlp(M, (SX) => {
    const x = SX.sym("x");
    const xm1 = M.minus(x, SX(1));
    return { x: x, f: M.times(xm1, xm1) };
  });
  const solver = M.nlpsol("S", plugin, nlp, quietOpts(plugin));
  const r = solver.call({ lbx: -10, ubx: 10 });
  assertAlmost(Number(r.x.nonzeros()[0]), 1, 6, `${plugin}: x`);
}

function test_lbg_default(M, plugin) {
  // Omitting lbg/ubg should default to +/- inf (unconstrained g).
  // Linear g=x; solution still drives x to 1 (the f-minimiser).
  const nlp = makeNlp(M, (SX) => {
    const x = SX.sym("x");
    const xm1 = M.minus(x, SX(1));
    return { x: x, f: M.times(xm1, xm1), g: x };
  });
  const solver = M.nlpsol("S", plugin, nlp, quietOpts(plugin));
  const r = solver.call({ lbx: -10, ubx: 10 });
  assertAlmost(Number(r.x.nonzeros()[0]), 1, 6, `${plugin}: x`);
}

// ----- driver -----

const tests = [
  ["test_basic",            test_basic],
  ["test_basic_par",        test_basic_par],
  ["test_inf_bounds",       test_inf_bounds],
  ["test_rosenbrock",       test_rosenbrock],
  ["test_rosenbrock_g",     test_rosenbrock_g],
  ["test_rosenbrock_fix",   test_rosenbrock_fix],
  ["test_initial_cond",     test_initial_cond],
  ["test_warmstart",        test_warmstart],
  ["test_x0_default",       test_x0_default],
  ["test_lbg_default",      test_lbg_default],
];

(async () => {
  const create = require(modulePath);
  const M = await create();

  // Run against ipopt -- the gold-standard NLP solver in the wasm
  // build.  sqpmethod loads but needs BLAS (`dgemm_`) at solve time
  // which isn't wired as a side module, so skip it here.  Other
  // wasm-native NLP solvers (qrsqp, feasiblesqpmethod) have the
  // same BLAS dependency.
  const solvers = [];
  if (M.has_nlpsol("ipopt")) solvers.push("ipopt");
  if (solvers.length === 0) {
    console.error("SKIP: no nlpsol plugin loaded.");
    process.exit(2);
  }

  let pass = 0, fail = 0;
  for (const plugin of solvers) {
    for (const [name, fn] of tests) {
      const label = `[${plugin}] ${name}`;
      try { fn(M, plugin); console.log(`ok   -- ${label}`); pass++; }
      catch (e) {
        console.log(`FAIL -- ${label}: ${e.message}`);
        fail++;
      }
    }
  }
  console.log(`\n${pass}/${pass + fail} passed`);
  process.exit(fail === 0 ? 0 : 1);
})().catch((e) => { console.error("FATAL:", e); process.exit(1); });

// daebuilder.ts -- Type-level coverage of the DaeBuilder API used in
// test/python/daebuilder.py.  DaeBuilder is the central object for
// constructing DAE / ODE / model-exchange descriptions in CasADi,
// used both as a parametric builder (add_x / add_u / add_eq) and as
// an FMU (.fmu) loader.  Focus here is on the builder surface --
// FMU loading needs filesystem fixtures we can't ship.

import {
  DaeBuilder, Function as CFn, MX, SX, DM,
  vertcat, plus, times, rdivide,
} from "./casadi";

function expectType<T>(_v: T): void {}

// ============================================================
// Construction
// ============================================================
expectType<DaeBuilder>(DaeBuilder());
const dae = DaeBuilder("rocket");
expectType<DaeBuilder>(dae);
const dae_with_path = DaeBuilder("rocket", "./does_not_exist_skip");
expectType<DaeBuilder>(dae_with_path);

// Add a variable with full causality/variability triple
const x = dae.add("x", "input", "continuous");
expectType<MX>(x);
const u = dae.add("u", "input", "discrete");
expectType<MX>(u);
const p = dae.add("p");                                 // bare-name form
expectType<MX>(p);
const y = dae.add("y", "output");                       // 2-arg form
expectType<MX>(y);

// Add a variable bound to an expression (the 5-arg form, returns void)
dae.add("c", "output", "continuous", plus(x, p));

// ============================================================
// Introspection: enumerate variables by category
// ============================================================
expectType<string[]>(dae.x());           // diff states
expectType<string[]>(dae.u());           // controls / inputs
expectType<string[]>(dae.y());           // outputs
expectType<string[]>(dae.p());           // parameters
expectType<string[]>(dae.w());           // dependent vars
expectType<string[]>(dae.t_new());       // independent var
expectType<MX>(dae.time());

// ============================================================
// Symbolic equations: set_init for IC
// ============================================================
dae.set_init("x", MX(1.0));

// ============================================================
// Build a Function from the DAE description
// ============================================================
const f = dae.create("f", ["x", "u"], ["c"]);
expectType<CFn>(f);
const f_opts = dae.create("f2", ["x", "u"], ["c"], { jit: false });
expectType<CFn>(f_opts);
// Default-name overload (4-arg with sx flag + lifted_calls)
const f_sx = dae.create("f_sx", ["x", "u"], ["c"], true);
expectType<CFn>(f_sx);
const f_sx_lifted = dae.create("f_sxL", ["x", "u"], ["c"], true, true);
expectType<CFn>(f_sx_lifted);

// ============================================================
// eliminate / register_lc / add_fun / etc.
// ============================================================
dae.eliminate("y");

// add_lc: compose a linear combination of outputs into a new alias
dae.add_lc("g_lin", ["c"]);

// add_fun: pull an existing Function into the DAE namespace
const helper = CFn("helper", [MX.sym("h")], [MX.sym("h")]);
const added = dae.add_fun(helper);
expectType<CFn>(added);

// ============================================================
// Demonstrate full builder -> Function -> evaluate flow
// ============================================================
const dae2 = DaeBuilder("cstr");
const C_A = dae2.add("C_A", "input", "continuous");
const C_B = dae2.add("C_B", "input", "continuous");
const q_in = dae2.add("q_in", "input", "discrete");
const C_A_in = dae2.add("C_A_in", "input", "discrete");
const C_B_in = dae2.add("C_B_in", "input", "discrete");
// Add an ODE-like output expression
const k = 0.5, V = 1;
dae2.add("dC_A", "output", "continuous",
  rdivide(plus(times(q_in, minusFn(C_A_in, C_A)), times(times(-V, k), times(C_A, C_B))), MX(V)));

// Helper to bridge JS-side `_DM` coercion gaps where minus(M, scalar)
// must use the function form.
function minusFn(a: MX, b: MX): MX { return plus(a, times(b, MX(-1))); }

void [dae, dae_with_path, x, u, p, y, f, f_opts, f_sx, f_sx_lifted, helper, added,
      dae2, C_A, C_B, q_in, C_A_in, C_B_in];

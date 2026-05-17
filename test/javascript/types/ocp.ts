// ocp.ts -- Type-level coverage of optimal-control problem patterns
// drawn from test/python/ocp.py.  Single-shooting, multi-shooting,
// integrator composition into nlpsol, detect_simple_bounds tuple
// return, sparsify+blockcat sparsity patterns.
//
// Goal: every tsc error here is either a real d.ts gap (fix in SWIG)
// or a real C++ API ergonomics gap (document in fixture).  No
// @ts-ignore, no `as any`.

import {
  Function as CFn, SX, MX, DM, Sparsity,
  integrator, nlpsol,
  blockcat, sparsify, detect_simple_bounds,
  vertcat, plus, minus, times, fmax,
} from "./casadi";

function expectType<T>(_v: T): void {}

// ============================================================
// Single-shooting NLP: build integrator over [0, te], then nlpsol on
// (initial-state, parameter) -> cost.
// ============================================================
const te = 0.4;

const t = SX.sym("t");
const q = SX.sym("y", 2);
const p = SX.sym("p", 1);
// ODE: dy/dt = [y1, p + y1^2]
const q0 = q.get(false, 0);  // q[0]
const q1 = q.get(false, 1);  // q[1]
const p0 = p.get(false, 0);
const ode = vertcat([q1, plus(p0, times(q1, q1))]);
const dae = { x: q, p: p, t: t, ode: ode };
const F = integrator("F", "rk", dae, 0, te, {
  reltol: 1e-15,
  abstol: 1e-15,
  verbose: false,
  steps_per_checkpoint: 10000,
});
expectType<CFn>(F);

// Integrator output: dict-form call returns Record<string, MX>.
// (Integrator was built with SX dae, but calling it with MX inputs
// produces MX outputs since the symbolic chaining picks the MX overload.)
const var_mx = MX.sym("var", 2);
const par_mx = MX.sym("par", 1);
const var0 = var_mx.get(false, 0);
const var1 = var_mx.get(false, 1);
const x0_mx = vertcat([var0, par_mx]);
const intg_call = F.call({ x0: x0_mx, p: var1 });
expectType<Record<string, MX>>(intg_call);
const xf = intg_call["xf"];
expectType<MX>(xf);

// Wrap (var, par) -> cost in a Function, then build an NLP on it.
const parc = MX(0);
const cost_expr = xf.get(false, 0);  // first component of xf
const cost_fn = CFn("f", [var_mx, par_mx], [cost_expr]);
const nlp = { x: var_mx, f: cost_fn.call([var_mx, parc])[0] };
const solver = nlpsol("solver", "ipopt", nlp, {
  "ipopt.tol": 1e-12,
  "ipopt.hessian_approximation": "limited-memory",
  "ipopt.max_iter": 10,
  "ipopt.print_level": 0,
});
expectType<CFn>(solver);

// Call solver with numerical bounds: returns DM dict.
const sol = solver.call({
  lbx: DM([-1, -1]),
  ubx: DM([1, 0.2]),
});
expectType<Record<string, DM>>(sol);
expectType<DM>(sol["x"]);
expectType<DM>(sol["f"]);
expectType<DM>(sol["lam_x"]);

// ============================================================
// Multi-shooting NLP: integrator + constraint g(var)
// ============================================================
const nlp_g = {
  x: var_mx,
  f: cost_fn.call([var_mx, parc])[0],
  g: minus(var0, var1),
};
const solver_g = nlpsol("solver_g", "ipopt", nlp_g, {});
expectType<CFn>(solver_g);
const sol_g = solver_g.call({
  lbx: DM([-1, -1]),
  ubx: DM([1, 0.2]),
  lbg: DM([-1]),
  ubg: DM([0]),
});
expectType<Record<string, DM>>(sol_g);
expectType<DM>(sol_g["lam_g"]);

// fmax post-processing
expectType<DM>(fmax(sol["lam_x"], DM(0)));

// ============================================================
// detect_simple_bounds: 5-tuple return via argout aggregation
// (the C++ signature has FOUR &OUTPUT params + bigint[] return).
// ============================================================
const xX = SX.sym("xX", 5);
const pSX = SX.sym("p", 2);
const gSX = vertcat([xX.get(false, 0), xX.get(false, 1)]);
const lbg = SX.zeros(2, 1);
const ubg = SX.zeros(2, 1);
const det_out = detect_simple_bounds(xX, pSX, gSX, lbg, ubg);
expectType<[bigint[], SX, SX, CFn, CFn]>(det_out);
const [idx_g, lb_x, ub_x, lam_to_lam, lin_g] = det_out;
expectType<bigint[]>(idx_g);
expectType<SX>(lb_x);
expectType<SX>(ub_x);
expectType<CFn>(lam_to_lam);
expectType<CFn>(lin_g);

// ============================================================
// blockcat + sparsify (sparsity construction for fatrop-style OCP
// detection adversarial cases)
// ============================================================
const D2 = sparsify(DM([[1, 0, 0], [1, 1, 1]])).sparsity();
expectType<Sparsity>(D2);
const block = blockcat([[DM([[1]]), DM([[2]])], [DM([[3]]), DM([[4]])]]);
expectType<DM>(block);
// Note: sparsify(MX) does NOT exist (only DM/SX) -- C++ API gap.
// At runtime MX expressions carry their own sparsity already.

// ============================================================
// fatrop options: structure_detection, equality flags, expand
// ============================================================
const xfat = MX.sym("x");
// Pattern from ocp.py test_bug: structure-detection auto/none for
// fatrop with equality-constraint flags.
for (const structure_detection of ["none", "auto"] as const) {
  const opts = {
    expand: true,
    structure_detection,
    equality: [true],
  };
  // g = x-1 constraint
  const sf1 = nlpsol("solver", "fatrop", { x: xfat, g: minus(xfat, MX(1)) }, opts);
  expectType<CFn>(sf1);
  const r1 = sf1.call({ lbg: DM(0), ubg: DM(0) });
  expectType<Record<string, DM>>(r1);

  // g = x
  const sf2 = nlpsol("solver", "fatrop", { x: xfat, g: xfat }, opts);
  expectType<CFn>(sf2);
  const r2 = sf2.call({ lbg: DM(1), ubg: DM(1) });
  expectType<Record<string, DM>>(r2);
}

// ============================================================
// Adversarial sparsity construction for fatrop detection
// (from ocp.py test_detect_adversarial)
// ============================================================
const D2_adv = sparsify(DM([[1, 0, 0], [1, 1, 1]])).sparsity();
const C2_adv = sparsify(DM([[0, 1], [0, 0]])).sparsity();
expectType<Sparsity>(D2_adv);
expectType<Sparsity>(C2_adv);
// Zero-pattern sparsity construction
const A1_zero = Sparsity(2, 2);
const B1_zero = Sparsity(2, 2);
expectType<Sparsity>(A1_zero);
expectType<Sparsity>(B1_zero);

// ============================================================
// Multi-shot integrator call with parameter sweep
// ============================================================
const t0_val = 0.0;
const tf_vals: number[] = [0.25, 0.5, 0.75, 1.0];
const F_multi = integrator("Fm", "rk", dae, t0_val, tf_vals);
expectType<CFn>(F_multi);
const r_multi = F_multi.call({ x0: DM([1, 0]), p: DM(0.2) });
expectType<Record<string, DM>>(r_multi);

// Pin locals
void [t, q, p, q0, q1, p0, ode, dae, F, var_mx, par_mx, var0, var1, x0_mx, intg_call, xf,
      parc, cost_expr, cost_fn, nlp, solver, sol, nlp_g, solver_g, sol_g,
      xX, pSX, gSX, lbg, ubg, det_out, idx_g, lb_x, ub_x, lam_to_lam, lin_g,
      D2, block, xfat, D2_adv, C2_adv, A1_zero, B1_zero, t0_val, tf_vals, F_multi, r_multi];

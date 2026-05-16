// nlpsol.ts -- Type-level coverage of solver constructors: nlpsol,
// qpsol, conic, integrator, rootfinder.  All take a name + plugin +
// problem-spec dict + options; the dict carries SX/MX expressions
// describing the problem (objective, constraints, ODE, ...).

import {
  Function as CFn, SX, MX, DM, Sparsity,
  nlpsol, qpsol, conic, integrator, rootfinder, interpolant,
  plus, minus, times, sumsqr,
} from "./casadi";

function expectType<T>(_v: T): void {}

// ============================================================
// NLP: min (x-1)^2 s.t. nothing
// ============================================================
const x_sx = SX.sym("x");
const f_sx = times(minus(x_sx, new SX(1)), minus(x_sx, new SX(1)));
const nlp_sx = { x: x_sx, f: f_sx };
const solver_sx = nlpsol("solver", "ipopt", nlp_sx, { print_time: false });
expectType<CFn>(solver_sx);

// MX-side problem
const x_mx = MX.sym("x");
const f_mx = times(minus(x_mx, new MX(1)), minus(x_mx, new MX(1)));
const nlp_mx = { x: x_mx, f: f_mx };
const solver_mx = nlpsol("solver", "ipopt", nlp_mx, {});
expectType<CFn>(solver_mx);

// With constraints g and parameters p
const p = MX.sym("p");
const g = plus(x_mx, p);
const nlp_full = { x: x_mx, f: f_mx, g: g, p: p };
expectType<CFn>(nlpsol("s", "ipopt", nlp_full));

// Call the solver: keys x0/lbx/ubx/lbg/ubg/p; returns {x, f, g, lam_x, lam_g, lam_p}.
// Solvers numerically evaluate regardless of the symbolic flavor they were
// constructed with -- passing DM-shaped input yields DM-shaped output.
const res_num = solver_mx.call({ x0: new DM(0) });
expectType<Record<string, DM>>(res_num);
// Symbolic chaining: pass MX expressions to compose into a bigger MX graph.
const x_outer = MX.sym("xo");
const res_sym = solver_mx.call({ x0: x_outer });
expectType<Record<string, MX>>(res_sym);

// ============================================================
// QP: conic (sparsity-only) and qpsol (SX/MX-typed)
// ============================================================
const H = Sparsity.diag(2n);
const A = Sparsity.dense(1n, 2n);
const qp_conic = { h: H, a: A };
expectType<CFn>(conic("qp", "qrqp", qp_conic));

const x_qp = MX.sym("x", 2n);
const qp_mx = { x: x_qp, f: sumsqr(x_qp), g: plus(x_qp.get(false, 0n), x_qp.get(false, 1n)) };
expectType<CFn>(qpsol("qp", "qrqp", qp_mx));

// ============================================================
// Integrator (ODE): dx/dt = x, x(0)=1 -> x(1)=e
// ============================================================
const x_ode = SX.sym("x");
const dae = { x: x_ode, ode: x_ode };
const F = integrator("F", "rk", dae);
expectType<CFn>(F);
const F2 = integrator("F", "rk", dae, 0, 1);
expectType<CFn>(F2);
const F3 = integrator("F", "rk", dae, 0, [0.0, 0.5, 1.0]);
expectType<CFn>(F3);

// ============================================================
// Rootfinder: x s.t. g(x, p) = 0
// ============================================================
const x_rf = SX.sym("x");
const p_rf = SX.sym("p");
const rfp = { x: x_rf, p: p_rf, g: minus(times(x_rf, x_rf), p_rf) };
expectType<CFn>(rootfinder("R", "newton", rfp));

// ============================================================
// Interpolant
// ============================================================
const interp = interpolant("LUT", "linear", [[0, 1, 2, 3]], [0, 1, 4, 9]);
expectType<CFn>(interp);

void [x_sx, f_sx, nlp_sx, solver_sx, x_mx, f_mx, nlp_mx, solver_mx, p, g, nlp_full,
      res_num, res_sym, x_outer, H, A, qp_conic, x_qp, qp_mx, x_ode, dae, F, F2, F3,
      x_rf, p_rf, rfp, interp];

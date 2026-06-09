//
//     MIT No Attribution
//
//     Copyright (C) 2010-2023 Joel Andersson, Joris Gillis, Moritz Diehl, KU Leuven.
//
//     Permission is hereby granted, free of charge, to any person obtaining a copy of this
//     software and associated documentation files (the "Software"), to deal in the Software
//     without restriction, including without limitation the rights to use, copy, modify,
//     merge, publish, distribute, sublicense, and/or sell copies of the Software, and to
//     permit persons to whom the Software is furnished to do so.
//
//     THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
//     INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
//     PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
//     HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
//     OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
//     SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
//
//
// JS port of docs/examples/python/implicit_runge-kutta.py.
//
// Build a fixed-step implicit Runge-Kutta integrator by hand: a collocation
// system solved with a rootfinder, iterated over finite elements, then
// forward/adjoint directional derivatives via Function.factory.
//
// JS notes (see README.md):
//   * The Python example also builds a 'cvodes' reference integrator and runs
//     the same sensitivity checks on it.  cvodes is not in the wasm build, so
//     we keep only the hand-built IRK integrator.
//   * The rootfinder uses 'newton' (kinsol is not available in wasm).

async function example(ca, log) {
  // End time
  const tf = 10.0;

  // Dimensions
  const nx = 3;
  const np = 1;

  // Declare variables
  const x = ca.SX.sym("x", nx);  // state
  const p = ca.SX.sym("u", np);  // control
  const [xa, xb] = ca.vertsplit(x); // x[0], x[1]

  // ODE right hand side
  const ode = ca.vertcat(
    ca.plus(ca.minus(ca.times(ca.minus(1, ca.times(xb, xb)), xa), xb), p),
    xa,
    ca.plus(ca.plus(ca.times(xa, xa), ca.times(xb, xb)), ca.times(p, p)));
  const f = ca.Function("f", [x, p], [ode]);

  // Number of finite elements
  const n = 100;
  // Size of the finite elements
  const h = tf / n;
  // Degree of interpolating polynomial
  const d = 4;

  // Collocation points
  const tau_root = [0, ...ca.collocation_points(d, "legendre")];

  // Coefficients of the collocation (C) and continuity (D) equations
  const C = Array.from({ length: d + 1 }, () => new Array(d + 1).fill(0));
  const D = new Array(d + 1).fill(0);

  // Dimensionless time inside one control interval
  const tau = ca.SX.sym("tau");

  for (let j = 0; j <= d; j++) {
    // Lagrange polynomial basis at collocation point j
    let Lp = ca.SX(1);
    for (let r = 0; r <= d; r++)
      if (r !== j) Lp = ca.times(Lp, ca.times(ca.minus(tau, tau_root[r]), 1 / (tau_root[j] - tau_root[r])));
    // Continuity coefficient
    const lfcn = ca.Function("lfcn", [tau], [Lp]);
    D[j] = Number(lfcn(1.0));
    // Collocation coefficients
    const tfcn = ca.Function("tfcn", [tau], [ca.tangent(Lp, tau)]);
    for (let r = 0; r <= d; r++) C[j][r] = Number(tfcn(tau_root[r]));
  }

  // Variables for one finite element
  const X0 = ca.MX.sym("X0", nx);
  const Pm = ca.MX.sym("P", np);
  const V = ca.MX.sym("V", d * nx);

  // State at each collocation point: [X0, V split into d blocks of nx]
  const offs = [];
  for (let r = 0; r <= d; r++) offs.push(r * nx);
  let X = [X0, ...ca.vertsplit(V, offs)];

  // Collocation equations defining V
  const V_eq = [];
  for (let j = 1; j <= d; j++) {
    let xp_j = ca.times(C[0][j], X[0]);
    for (let r = 1; r <= d; r++) xp_j = ca.plus(xp_j, ca.times(C[r][j], X[r]));
    const f_j = f(X[j], Pm);
    V_eq.push(ca.minus(ca.times(h, f_j), xp_j));
  }
  const V_eq_cat = ca.vcat(V_eq);

  // Root-finding function, implicitly defines V as a function of X0 and P
  const vfcn = ca.Function("vfcn", [V, X0, Pm], [V_eq_cat]);
  const vfcn_sx = vfcn.expand();  // convert to SX to decrease overhead

  // Implicit function instance (newton, not kinsol)
  const ifcn = ca.rootfinder("ifcn", "newton", vfcn_sx);
  const Vsol = ifcn(ca.MX(), X0, Pm);

  // Recover states and form end-of-element state
  const Vparts = ca.vertsplit(Vsol, offs);
  X = [X0, ...Vparts];
  let XF = ca.times(D[0], X[0]);
  for (let r = 1; r <= d; r++) XF = ca.plus(XF, ca.times(D[r], X[r]));

  // Discrete time dynamics for one finite element
  const F = ca.Function("F", [X0, Pm], [XF]);

  // Iterate over all finite elements
  let Xacc = X0;
  for (let i = 0; i < n; i++) Xacc = F(Xacc, Pm);

  // Fixed-step integrator (as a plain Function over x0, p -> xf)
  const irk = ca.Function("irk_integrator", [X0, Pm], [Xacc], ["x0", "p"], ["xf"]);

  // Test values
  const x0_val = ca.DM([0, 1, 0]);
  const p_val = ca.DM(0.2);

  log("-------");
  log("Testing " + irk.name());
  log("-------");

  // Forward and reverse directional derivatives
  const dF = irk.factory("dF",
    ["x0", "p", "fwd:x0", "fwd:p", "adj:xf"],
    ["xf", "fwd:xf", "adj:x0", "adj:p"]);

  const res = dF.call({
    x0: x0_val, p: p_val,
    fwd_x0: ca.DM([1, 0, 0]), fwd_p: ca.DM(1),
    adj_xf: ca.DM([0, 0, 1]),
  });

  log("xf = " + res["xf"]);
  log("d(xf)/d(p)+d(xf)/d(x0[0]) = " + res["fwd_xf"]);
  log("d(xf[2])/d(x0) = " + res["adj_x0"]);
  log("d(xf[2])/d(p) = " + res["adj_p"]);
}

if (typeof require !== "undefined" && typeof module !== "undefined" && require.main === module) {
  const path = require("path");
  const casadiPath = process.env.CASADI_JS
    || path.resolve(__dirname, "../../../build-wasm/swig/wasm-js/casadi.js");
  require(casadiPath)()
    .then((ca) => example(ca, (...a) => console.log(...a)))
    .catch((e) => { console.error("FATAL:", e.message || e); process.exit(1); });
}

if (typeof module !== "undefined" && module.exports) module.exports = example;

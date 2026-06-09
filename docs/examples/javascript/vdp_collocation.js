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
// JS port of docs/examples/python/vdp_collocation.py.
//
// Van der Pol OCP transcribed by hand with Radau collocation; the NLP
// variable vector V is sliced manually (no integrator plugin) -> portable.
//
// JS notes (see README.md):
//   * numpy poly1d/polyder/polyint are replaced by small array helpers.
//   * matplotlib output is dropped; we log the optimum instead.
//   * The Python `expand=True` and `ipopt.linear_solver='ma27'` options
//     are dropped (defaults used; expand crashes the wasm GC teardown).

// Polynomial helpers (coeffs highest-degree-first, like numpy.poly1d).
function polymul(a, b) {
  const r = new Array(a.length + b.length - 1).fill(0);
  for (let i = 0; i < a.length; i++)
    for (let j = 0; j < b.length; j++) r[i + j] += a[i] * b[j];
  return r;
}
function polyval(p, x) { let r = 0; for (const c of p) r = r * x + c; return r; }
function polyder(p) {
  const n = p.length - 1, r = [];
  for (let i = 0; i < n; i++) r.push(p[i] * (n - i));
  return r.length ? r : [0];
}
function polyint(p) {
  const n = p.length, r = [];
  for (let i = 0; i < n; i++) r.push(p[i] / (n - i));
  r.push(0);
  return r;
}

async function example(ca, log) {
  const inf = Infinity;

  // Degree of interpolating polynomial
  const d = 3;

  // Choose collocation points
  const tau_root = [0, ...ca.collocation_points(d, "radau")];

  // Coefficients of the collocation (C), continuity (D), quadrature (F) eqs
  const C = Array.from({ length: d + 1 }, () => new Array(d + 1).fill(0));
  const D = new Array(d + 1).fill(0);
  const Fc = new Array(d + 1).fill(0);

  // Construct polynomial basis
  for (let j = 0; j <= d; j++) {
    let p = [1];
    for (let r = 0; r <= d; r++)
      if (r !== j) p = polymul(p, [1, -tau_root[r]]).map((c) => c / (tau_root[j] - tau_root[r]));
    D[j] = polyval(p, 1.0);
    const pder = polyder(p);
    for (let r = 0; r <= d; r++) C[j][r] = polyval(pder, tau_root[r]);
    Fc[j] = polyval(polyint(p), 1.0);
  }

  // Control discretization
  const nk = 20;
  const tf = 10.0;     // End time
  const h = tf / nk;   // Size of the finite elements

  // Declare variables (use scalar graph)
  const t = ca.SX.sym("t");
  const u = ca.SX.sym("u");
  const x = ca.SX.sym("x", 2);
  const [x0e, x1e] = ca.vertsplit(x);

  // ODE rhs function and quadratures
  const xdot = ca.vertcat(
    ca.plus(ca.minus(ca.times(ca.minus(1, ca.times(x1e, x1e)), x0e), x1e), u),
    x0e);
  const qdot = ca.plus(ca.plus(ca.times(x0e, x0e), ca.times(x1e, x1e)), ca.times(u, u));
  const f = ca.Function("f", [t, x, u], [xdot, qdot], ["t", "x", "u"], ["xdot", "qdot"]);

  // Control bounds
  const u_min = -0.75, u_max = 1.0, u_init = 0.0;

  // State bounds and initial guess
  const x_min = [-inf, -inf], x_max = [inf, inf];
  const xi_min = [0.0, 1.0], xi_max = [0.0, 1.0];
  const xf_min = [0.0, 0.0], xf_max = [0.0, 0.0];
  const x_init = [0.0, 0.0];

  const nx = 2, nu = 1;

  // Total number of variables
  const NX = nk * (d + 1) * nx;  // Collocated states
  const NU = nk * nu;            // Parametrized controls
  const NXF = nx;                // Final state
  const NV = NX + NU + NXF;

  // NLP variable vector
  const V = ca.MX.sym("V", NV);
  const Vs = ca.vertsplit(V); // scalar entries

  // All variables with bounds and initial guess
  const vars_lb = new Array(NV).fill(0);
  const vars_ub = new Array(NV).fill(0);
  const vars_init = new Array(NV).fill(0);
  let offset = 0;

  // Get collocated states and parametrized control as MX sub-vectors
  const X = Array.from({ length: nk + 1 }, () => new Array(d + 1));
  const U = new Array(nk);
  const slice = (from, n) => ca.vcat(Vs.slice(from, from + n));
  for (let k = 0; k < nk; k++) {
    for (let j = 0; j <= d; j++) {
      X[k][j] = slice(offset, nx);
      vars_init.splice(offset, nx, ...x_init);
      if (k === 0 && j === 0) {
        vars_lb.splice(offset, nx, ...xi_min); vars_ub.splice(offset, nx, ...xi_max);
      } else {
        vars_lb.splice(offset, nx, ...x_min); vars_ub.splice(offset, nx, ...x_max);
      }
      offset += nx;
    }
    U[k] = slice(offset, nu);
    vars_lb[offset] = u_min; vars_ub[offset] = u_max; vars_init[offset] = u_init;
    offset += nu;
  }
  // State at end time
  X[nk][0] = slice(offset, nx);
  vars_lb.splice(offset, nx, ...xf_min); vars_ub.splice(offset, nx, ...xf_max);
  vars_init.splice(offset, nx, ...x_init);
  offset += nx;

  // Constraints and objective
  const g = [], lbg = [], ubg = [];
  let J = ca.MX(0);

  for (let k = 0; k < nk; k++) {
    for (let j = 1; j <= d; j++) {
      // State derivative at the collocation point
      let xp_jk = ca.times(C[0][j], X[k][0]);
      for (let r = 1; r <= d; r++) xp_jk = ca.plus(xp_jk, ca.times(C[r][j], X[k][r]));

      // Collocation equations
      const Tkj = h * (k + tau_root[j]);
      const [ode_jk, ltermv] = f(ca.DM(Tkj), X[k][j], U[k]);
      g.push(ca.minus(ca.times(h, ode_jk), xp_jk));
      lbg.push(0, 0); ubg.push(0, 0);

      // Objective contribution
      J = ca.plus(J, ca.times(Fc[j], ca.times(ltermv, h)));
    }
    // State at the end of the finite element
    let xf_k = ca.times(D[0], X[k][0]);
    for (let r = 1; r <= d; r++) xf_k = ca.plus(xf_k, ca.times(D[r], X[k][r]));

    // Continuity equation
    g.push(ca.minus(X[k + 1][0], xf_k));
    lbg.push(0, 0); ubg.push(0, 0);
  }

  const G = ca.vcat(g);
  const nlp = { x: V, f: J, g: G };
  const solver = ca.nlpsol("solver", "ipopt", nlp);

  const res = solver.call({
    x0: ca.DM(vars_init), lbx: ca.DM(vars_lb), ubx: ca.DM(vars_ub),
    lbg: ca.DM(lbg), ubg: ca.DM(ubg),
  });

  log("-----");
  log("optimal cost: " + res["f"]);

  // Get values at the beginning of each finite element (Python: v_opt[0::stride])
  const v_opt = res["x"];
  const stride = (d + 1) * nx + nu;
  log("x0 trajectory = " + v_opt["0::" + stride]);
  log("x1 trajectory = " + v_opt["1::" + stride]);
  log("u  trajectory = " + v_opt[(d + 1) * nx + "::" + stride]);
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

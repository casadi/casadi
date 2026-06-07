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

async function example(M, log) {
  const inf = Infinity;

  // Degree of interpolating polynomial
  const d = 3;

  // Choose collocation points
  const tau_root = [0, ...M.collocation_points(d, "radau")];

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
  const t = M.SX.sym("t");
  const u = M.SX.sym("u");
  const x = M.SX.sym("x", 2);
  const [x0e, x1e] = M.vertsplit(x);

  // ODE rhs function and quadratures
  const xdot = M.vertcat(
    M.plus(M.minus(M.times(M.minus(1, M.times(x1e, x1e)), x0e), x1e), u),
    x0e);
  const qdot = M.plus(M.plus(M.times(x0e, x0e), M.times(x1e, x1e)), M.times(u, u));
  const f = new M.Function("f", [t, x, u], [xdot, qdot], ["t", "x", "u"], ["xdot", "qdot"]);

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
  const V = M.MX.sym("V", NV);
  const Vs = M.vertsplit(V); // scalar entries

  // All variables with bounds and initial guess
  const vars_lb = new Array(NV).fill(0);
  const vars_ub = new Array(NV).fill(0);
  const vars_init = new Array(NV).fill(0);
  let offset = 0;

  // Get collocated states and parametrized control as MX sub-vectors
  const X = Array.from({ length: nk + 1 }, () => new Array(d + 1));
  const U = new Array(nk);
  const slice = (from, n) => M.vertcat(...Vs.slice(from, from + n));
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
  let J = M.MX(0);

  for (let k = 0; k < nk; k++) {
    for (let j = 1; j <= d; j++) {
      // State derivative at the collocation point
      let xp_jk = M.times(C[0][j], X[k][0]);
      for (let r = 1; r <= d; r++) xp_jk = M.plus(xp_jk, M.times(C[r][j], X[k][r]));

      // Collocation equations
      const Tkj = h * (k + tau_root[j]);
      const out = f.call([M.DM(Tkj), X[k][j], U[k]]);
      g.push(M.minus(M.times(h, out[0]), xp_jk));
      lbg.push(0, 0); ubg.push(0, 0);

      // Objective contribution
      J = M.plus(J, M.times(Fc[j], M.times(out[1], h)));
    }
    // State at the end of the finite element
    let xf_k = M.times(D[0], X[k][0]);
    for (let r = 1; r <= d; r++) xf_k = M.plus(xf_k, M.times(D[r], X[k][r]));

    // Continuity equation
    g.push(M.minus(X[k + 1][0], xf_k));
    lbg.push(0, 0); ubg.push(0, 0);
  }

  const G = M.vcat(g);
  const nlp = { x: V, f: J, g: G };
  const solver = M.nlpsol("solver", "ipopt", nlp);

  const res = solver.call({
    x0: M.DM(vars_init), lbx: M.DM(vars_lb), ubx: M.DM(vars_ub),
    lbg: M.DM(lbg), ubg: M.DM(ubg),
  });

  log("-----");
  log("optimal cost: " + res["f"].nonzeros().join(" "));

  // Get values at the beginning of each finite element
  const v_opt = res["x"].nonzeros();
  const stride = (d + 1) * nx + nu;
  const x0_opt = [], x1_opt = [], u_opt = [];
  for (let i = 0; i < v_opt.length; i += stride) {
    x0_opt.push(v_opt[i]);
    if (i + 1 < v_opt.length) x1_opt.push(v_opt[i + 1]);
  }
  for (let i = (d + 1) * nx; i < v_opt.length; i += stride) u_opt.push(v_opt[i]);
  log("x0 trajectory = " + x0_opt.map((v) => v.toFixed(4)).join(" "));
  log("x1 trajectory = " + x1_opt.map((v) => v.toFixed(4)).join(" "));
  log("u  trajectory = " + u_opt.map((v) => v.toFixed(4)).join(" "));
}

if (typeof require !== "undefined" && typeof module !== "undefined" && require.main === module) {
  const path = require("path");
  const casadiPath = process.env.CASADI_JS
    || path.resolve(__dirname, "../../../build-wasm/swig/wasm-js/casadi.js");
  require(casadiPath)()
    .then((M) => example(M, (...a) => console.log(...a)))
    .catch((e) => { console.error("FATAL:", e.message || e); process.exit(1); });
}

if (typeof module !== "undefined" && module.exports) module.exports = example;

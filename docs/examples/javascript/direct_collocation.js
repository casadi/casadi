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
// JS port of docs/examples/python/direct_collocation.py.
//
// Direct collocation transcription of a Van der Pol OCP, building the
// collocation coefficients by hand (no integrator plugin) -> fully portable.
//
// JS notes (see README.md):
//   * numpy poly1d/polyder/polyint are replaced by small array-based helpers.
//   * ca.collocation_points(d, ...) returns a plain JS number[] here.
//   * matplotlib output is dropped; we log the optimum instead.

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
  // Degree of interpolating polynomial
  const d = 3;

  // Collocation points (prepend 0)
  const tau_root = [0, ...ca.collocation_points(d, "legendre")];

  // Coefficients of the collocation, continuity and quadrature equations
  const C = Array.from({ length: d + 1 }, () => new Array(d + 1).fill(0));
  const D = new Array(d + 1).fill(0);
  const B = new Array(d + 1).fill(0);

  // Construct polynomial basis
  for (let j = 0; j <= d; j++) {
    let p = [1];
    for (let r = 0; r <= d; r++)
      if (r !== j) p = polymul(p, [1, -tau_root[r]]).map((c) => c / (tau_root[j] - tau_root[r]));
    // Continuity coefficient
    D[j] = polyval(p, 1.0);
    // Collocation coefficients
    const pder = polyder(p);
    for (let r = 0; r <= d; r++) C[j][r] = polyval(pder, tau_root[r]);
    // Quadrature coefficient
    B[j] = polyval(polyint(p), 1.0);
  }

  // Time horizon
  const T = 10.0;

  // Declare model variables
  const x1 = ca.SX.sym("x1");
  const x2 = ca.SX.sym("x2");
  const x = ca.vertcat(x1, x2);
  const u = ca.SX.sym("u");

  // Model equations
  const xdot = ca.vertcat(
    ca.plus(ca.minus(ca.times(ca.minus(1, ca.times(x2, x2)), x1), x2), u),
    x1);

  // Objective term
  const L = ca.plus(ca.plus(ca.times(x1, x1), ca.times(x2, x2)), ca.times(u, u));

  // Continuous time dynamics
  const f = ca.Function("f", [x, u], [xdot, L], ["x", "u"], ["xdot", "L"]);

  // Control discretization
  const N = 20;       // number of control intervals
  const h = T / N;

  // Start with an empty NLP
  const w = [], w0 = [], lbw = [], ubw = [];
  let J = ca.MX(0);
  const g = [], lbg = [], ubg = [];

  // For plotting/logging x and u given w
  const x_plot = [], u_plot = [];

  // "Lift" initial conditions
  let Xk = ca.MX.sym("X0", 2);
  w.push(Xk); lbw.push(0, 1); ubw.push(0, 1); w0.push(0, 1); x_plot.push(Xk);

  // Formulate the NLP
  for (let k = 0; k < N; k++) {
    // New NLP variable for the control
    const Uk = ca.MX.sym("U_" + k);
    w.push(Uk); lbw.push(-1); ubw.push(1); w0.push(0); u_plot.push(Uk);

    // State at collocation points
    const Xc = [];
    for (let j = 0; j < d; j++) {
      const Xkj = ca.MX.sym("X_" + k + "_" + j, 2);
      Xc.push(Xkj); w.push(Xkj);
      lbw.push(-0.25, -Infinity); ubw.push(Infinity, Infinity); w0.push(0, 0);
    }

    // Loop over collocation points
    let Xk_end = ca.times(D[0], Xk);
    for (let j = 1; j <= d; j++) {
      // Expression for the state derivative at the collocation point
      let xp = ca.times(C[0][j], Xk);
      for (let r = 0; r < d; r++) xp = ca.plus(xp, ca.times(C[r + 1][j], Xc[r]));

      // Append collocation equations
      const [fj, qj] = f(Xc[j - 1], Uk);
      g.push(ca.minus(ca.times(h, fj), xp)); lbg.push(0, 0); ubg.push(0, 0);

      // Add contribution to the end state
      Xk_end = ca.plus(Xk_end, ca.times(D[j], Xc[j - 1]));

      // Add contribution to quadrature function
      J = ca.plus(J, ca.times(B[j], ca.times(qj, h)));
    }

    // New NLP variable for state at end of interval
    Xk = ca.MX.sym("X_" + (k + 1), 2);
    w.push(Xk); lbw.push(-0.25, -Infinity); ubw.push(Infinity, Infinity); w0.push(0, 0);
    x_plot.push(Xk);

    // Add equality constraint
    g.push(ca.minus(Xk_end, Xk)); lbg.push(0, 0); ubg.push(0, 0);
  }

  // Concatenate vectors
  const W = ca.vcat(w);
  const G = ca.vcat(g);
  const x_plot_m = ca.hcat(x_plot);
  const u_plot_m = ca.hcat(u_plot);

  // Create an NLP solver
  const prob = { f: J, x: W, g: G };
  const solver = ca.nlpsol("solver", "ipopt", prob);

  // Function to get x and u trajectories from w
  const trajectories = ca.Function("trajectories", [W], [x_plot_m, u_plot_m], ["w"], ["x", "u"]);

  // Solve the NLP
  const sol = solver.call({
    x0: ca.DM(w0), lbx: ca.DM(lbw), ubx: ca.DM(ubw), lbg: ca.DM(lbg), ubg: ca.DM(ubg),
  });
  const [x_opt, u_opt] = trajectories(sol["x"]);

  // Report the solution
  log("-----");
  log("objective at solution = " + sol["f"]);
  log("x1 trajectory = " + x_opt["0,:"]);
  log("x2 trajectory = " + x_opt["1,:"]);
  log("u trajectory = " + u_opt);
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

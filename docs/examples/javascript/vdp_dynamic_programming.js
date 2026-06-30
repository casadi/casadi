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
// JS port of docs/examples/python/vdp_dynamic_programming.py.
//
// Dynamic programming over a discretized state/control grid for the Van der
// Pol problem. The Python original is pure numpy (no casadi symbolics); this
// port keeps it as plain JS array math.
//
// JS notes (see README.md):
//   * numpy meshgrid/linspace become flat Float64Array grids.
//   * matplotlib output is dropped; we log the optimal cost and trajectory.

function linspace(a, b, n) {
  const r = new Array(n);
  for (let i = 0; i < n; i++) r[i] = a + (b - a) * i / (n - 1);
  return r;
}

async function example(ca, log) {
  const inf = Infinity;

  const T = 10.0;     // End time
  const N = 20;       // Number of control intervals
  const NK = 20;      // RK4 steps per interval
  const DT = T / (N * NK);
  const NU = 101;     // Number of discrete control values
  const NX = 101;     // Number of discrete state values per axis

  // System dynamics
  const f = (x1, x2, u) => [
    (1 - x2 * x2) * x1 - x2 + u, // x1_dot
    x1,                          // x2_dot
    x1 * x1 + x2 * x2 + u * u,   // q_dot
  ];

  const U = linspace(-1, 1, NU);
  const x1 = linspace(-1, 1, NX);
  const x2 = linspace(-1, 1, NX);

  const NG = NX * NX;            // total grid points (row=x2, col=x1)
  const idx = (i2, i1) => i2 * NX + i1;

  // For each control action, precompute next-state index and stage cost
  const stage_J = [], next_x1 = [], next_x2 = [];
  for (let uind = 0; uind < NU; uind++) {
    const u = U[uind];
    const nx1 = new Int32Array(NG), nx2 = new Int32Array(NG), Q = new Float64Array(NG);
    for (let i2 = 0; i2 < NX; i2++) {
      for (let i1 = 0; i1 < NX; i1++) {
        let X1 = x1[i1], X2 = x2[i2], Qk = 0;
        for (let k = 0; k < NK; k++) {
          const [a1, a2, aq] = f(X1, X2, u);
          const [b1, b2, bq] = f(X1 + DT / 2 * a1, X2 + DT / 2 * a2, u);
          const [c1, c2, cq] = f(X1 + DT / 2 * b1, X2 + DT / 2 * b2, u);
          const [d1, d2, dq] = f(X1 + DT * c1, X2 + DT * c2, u);
          X1 += DT / 6 * (a1 + 2 * b1 + 2 * c1 + d1);
          X2 += DT / 6 * (a2 + 2 * b2 + 2 * c2 + d2);
          Qk += DT / 6 * (aq + 2 * bq + 2 * cq + dq);
        }
        // Round to nearest grid index
        let r1 = Math.round((X1 + 1) / 2 * (NX - 1));
        let r2 = Math.round((X2 + 1) / 2 * (NX - 1));
        // Infinite cost if out-of-bounds
        if (r1 < 0 || r1 >= NX || r2 < 0 || r2 >= NX) { Qk = inf; r1 = 0; r2 = 0; }
        const g = idx(i2, i1);
        nx1[g] = r1; nx2[g] = r2; Q[g] = Qk;
      }
    }
    next_x1.push(nx1); next_x2.push(nx2); stage_J.push(Q);
  }

  // Cost-to-go (no end cost) and optimal control
  let J = new Float64Array(NG);
  const U_opt = [];
  for (let k = N - 1; k >= 0; k--) {
    const J_prev = new Float64Array(NG).fill(inf);
    const u_prev = new Int32Array(NG).fill(-1);
    for (let uind = 0; uind < NU; uind++) {
      const nx1 = next_x1[uind], nx2 = next_x2[uind], sj = stage_J[uind];
      for (let g = 0; g < NG; g++) {
        const test = J[idx(nx2[g], nx1[g])] + sj[g];
        if (test < J_prev[g]) { J_prev[g] = test; u_prev[g] = uind; }
      }
    }
    J = J_prev;
    U_opt.push(u_prev);
  }
  U_opt.reverse();

  // Optimal control starting at x1=0, x2=1
  let i1 = Math.floor(NX / 2);
  let i2 = NX - 1;
  const u_opt = [], x1_opt = [x1[i1]], x2_opt = [x2[i2]];
  let cost = 0;
  for (let k = 0; k < N; k++) {
    const u_ind = U_opt[k][idx(i2, i1)];
    cost += stage_J[u_ind][idx(i2, i1)];
    const ni1 = next_x1[u_ind][idx(i2, i1)], ni2 = next_x2[u_ind][idx(i2, i1)];
    i1 = ni1; i2 = ni2;
    u_opt.push(U[u_ind]); x1_opt.push(x1[i1]); x2_opt.push(x2[i2]);
  }

  log("-----");
  log("Minimal cost: " + cost);
  // Consistency check (cf. Python assert)
  if (Math.abs(cost - J[idx(NX - 1, Math.floor(NX / 2))]) >= 1e-8)
    throw new Error("consistency check failed");
  log("x1 trajectory = " + x1_opt);
  log("x2 trajectory = " + x2_opt);
  log("u  trajectory = " + u_opt);
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

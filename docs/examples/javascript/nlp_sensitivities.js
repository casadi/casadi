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
// JS port of docs/examples/python/nlp_sensitivities.py.
//
// Direct-collocation Van der Pol OCP parametrized by a perturbation P of
// the initial state, then sensitivities of the optimum w.r.t. P via three
// routes: solver.factory (Hessian), solver.forward, solver.reverse, with a
// finite-difference cross-check.
//
// JS notes (see README.md):
//   * numpy poly1d/polyder/polyint replaced by small array-based helpers.
//   * matplotlib output is dropped; we log the sensitivities instead.

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
  // Degree of interpolating polynomial
  const d = 3;
  const tau_root = [0, ...M.collocation_points(d, "legendre")];

  // Collocation / continuity / quadrature coefficients
  const C = Array.from({ length: d + 1 }, () => new Array(d + 1).fill(0));
  const D = new Array(d + 1).fill(0);
  const B = new Array(d + 1).fill(0);
  for (let j = 0; j <= d; j++) {
    let p = [1];
    for (let r = 0; r <= d; r++)
      if (r !== j) p = polymul(p, [1, -tau_root[r]]).map((c) => c / (tau_root[j] - tau_root[r]));
    D[j] = polyval(p, 1.0);
    const pder = polyder(p);
    for (let r = 0; r <= d; r++) C[j][r] = polyval(pder, tau_root[r]);
    B[j] = polyval(polyint(p), 1.0);
  }

  // Time horizon
  const T = 10.0;

  // Model variables
  const x1 = M.SX.sym("x1");
  const x2 = M.SX.sym("x2");
  const x = M.vertcat(x1, x2);
  const u = M.SX.sym("u");
  const xdot = M.vertcat(M.plus(M.minus(M.times(M.minus(1, M.times(x2, x2)), x1), x2), u), x1);
  const L = M.plus(M.plus(M.times(x1, x1), M.times(x2, x2)), M.times(u, u));
  const f = new M.Function("f", [x, u], [xdot, L], ["x", "u"], ["xdot", "L"]);

  // Control discretization
  const N = 20;
  const h = T / N;

  // Empty NLP
  const w = [], w0 = [], lbw = [], ubw = [];
  let J = M.MX(0);
  const g = [], lbg = [], ubg = [];
  const x_plot = [], u_plot = [];

  // "Lift" initial conditions
  let Xk = M.MX.sym("X0", 2);
  w.push(Xk); lbw.push(0, 1); ubw.push(0, 1); w0.push(0, 1); x_plot.push(Xk);

  // Perturb with P
  const P = M.MX.sym("P", 2);
  Xk = M.plus(Xk, P);

  // Formulate the NLP
  for (let k = 0; k < N; k++) {
    const Uk = M.MX.sym("U_" + k);
    w.push(Uk); lbw.push(-1); ubw.push(0.85); w0.push(0); u_plot.push(Uk);

    const Xc = [];
    for (let j = 0; j < d; j++) {
      const Xkj = M.MX.sym("X_" + k + "_" + j, 2);
      Xc.push(Xkj); w.push(Xkj);
      lbw.push(-0.25, -Infinity); ubw.push(Infinity, Infinity); w0.push(0, 0);
    }

    let Xk_end = M.times(D[0], Xk);
    for (let j = 1; j <= d; j++) {
      let xp = M.times(C[0][j], Xk);
      for (let r = 0; r < d; r++) xp = M.plus(xp, M.times(C[r + 1][j], Xc[r]));
      const fj = f.call([Xc[j - 1], Uk]);
      g.push(M.minus(M.times(h, fj[0]), xp)); lbg.push(0, 0); ubg.push(0, 0);
      Xk_end = M.plus(Xk_end, M.times(D[j], Xc[j - 1]));
      J = M.plus(J, M.times(B[j], M.times(fj[1], h)));
    }

    Xk = M.MX.sym("X_" + (k + 1), 2);
    w.push(Xk); lbw.push(-0.25, -Infinity); ubw.push(Infinity, Infinity); w0.push(0, 0);
    x_plot.push(Xk);
    g.push(M.minus(Xk_end, Xk)); lbg.push(0, 0); ubg.push(0, 0);
  }

  const W = M.vcat(w);
  const G = M.vcat(g);

  // NLP, using SQP + active-set QP for accurate multipliers
  const prob = { f: J, x: W, g: G, p: P };
  const opts = {
    qpsol: "qrqp",
    qpsol_options: { print_iter: false, error_on_fail: false },
    print_time: false,
  };
  const solver = M.nlpsol("solver", "sqpmethod", prob, opts);

  const DMlbw = M.DM(lbw), DMubw = M.DM(ubw), DMlbg = M.DM(lbg), DMubg = M.DM(ubg);

  // Solve the NLP
  const sol = solver.call({ x0: M.DM(w0), lbx: DMlbw, ubx: DMubw, lbg: DMlbg, ubg: DMubg, p: M.DM(0) });
  log("-----");
  log("objective at solution = " + sol["f"].nonzeros().join(" "));

  const nx = Number(W.size1());
  const ng = Number(G.size1());

  // High-level: Hessian of optimal f w.r.t. p via factory
  const hsolver = solver.factory("h", solver.name_in(), ["hess:f:p:p"]);
  log("hsolver generated");
  const hsol = hsolver.call({
    x0: sol["x"], lam_x0: sol["lam_x"], lam_g0: sol["lam_g"],
    lbx: DMlbw, ubx: DMubw, lbg: DMlbg, ubg: DMubg, p: M.DM(0),
  });
  log("Hessian of f w.r.t. p (2x2) = " + hsol["hess_f_p_p"].nonzeros().map((v) => v.toFixed(6)).join(" "));

  // Low-level forward AD: two directions at once
  const nfwd = 2;
  const zx = M.DM.zeros(nx, 1);
  const zg = M.DM.zeros(ng, 1);
  const fwd_lbx = M.DM.zeros(nx, nfwd);
  const fwd_ubx = M.DM.zeros(nx, nfwd);
  const fwd_lbg = M.DM.zeros(ng, nfwd);
  const fwd_ubg = M.DM.zeros(ng, nfwd);
  const fwd_p = M.DM.zeros(2, nfwd);   // perturb P
  fwd_p.set(1, false, 0, 0);           // fwd_p[0][0] = 1
  fwd_p.set(1, false, 1, 1);           // fwd_p[1][1] = 1

  const fwd_solver = solver.forward(nfwd);
  log("fwd_solver generated");
  const sol_fwd = fwd_solver.call({
    out_x: sol["x"], out_lam_g: sol["lam_g"], out_lam_x: sol["lam_x"],
    out_f: sol["f"], out_g: sol["g"], lbx: DMlbw, ubx: DMubw, lbg: DMlbg, ubg: DMubg,
    fwd_lbx: fwd_lbx, fwd_ubx: fwd_ubx, fwd_lbg: fwd_lbg, fwd_ubg: fwd_ubg,
    p: M.DM(0), fwd_p: fwd_p,
  });

  // Finite-difference cross-check of d f / d P
  const hfd = 1e-3;
  const fdpert = [];
  for (let dir = 0; dir < nfwd; dir++) {
    const dp = M.DM([dir === 0 ? hfd : 0, dir === 1 ? hfd : 0]);
    const r = solver.call({
      x0: sol["x"], lam_g0: sol["lam_g"], lam_x0: sol["lam_x"],
      lbx: DMlbw, ubx: DMubw, lbg: DMlbg, ubg: DMubg, p: dp,
    });
    fdpert.push(r["f"].nonzeros()[0]);
  }
  const f0 = sol["f"].nonzeros()[0];
  log("-----");
  log("d f / d P (finite differences) = " + fdpert.map((v) => ((v - f0) / hfd).toFixed(6)).join(" "));
  log("d f / d P (forward AD)         = " + sol_fwd["fwd_f"].nonzeros().map((v) => v.toFixed(6)).join(" "));

  // Reverse AD: which inputs influence f
  const nadj = 1;
  const adj_f = M.DM.zeros(1, nadj);
  adj_f.set(1, false, 0, 0);
  const adj_solver = solver.reverse(nadj);
  log("adj_solver generated");
  const sol_adj = adj_solver.call({
    out_x: sol["x"], out_lam_g: sol["lam_g"], out_lam_x: sol["lam_x"],
    out_f: sol["f"], out_g: sol["g"], lbx: DMlbw, ubx: DMubw, lbg: DMlbg, ubg: DMubg,
    adj_f: adj_f, adj_g: M.DM.zeros(ng, nadj), p: M.DM(0), adj_x: M.DM.zeros(nx, nadj),
  });
  log("d f / d P (reverse AD)         = " + sol_adj["adj_p"].nonzeros().map((v) => v.toFixed(6)).join(" "));
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

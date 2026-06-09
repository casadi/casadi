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
// JS port of docs/examples/python/biegler_10_1.py.
//
// Exercise 1, chapter 10 from Larry Biegler's book: orthogonal collocation
// applied to a scalar ODE, solved as an NLP with ipopt for N = 1..10 elements.
//
// JS notes (see README.md):
//   * Collocation coefficients C/D are computed numerically (Function eval)
//     into plain JS arrays so they can be indexed directly.
//   * The collocated states Z are kept as a JS [i][j] grid of scalar SX
//     symbols; x is their column-major vectorisation (vec(Z.T)).
//   * Plotting is dropped; the optimal cost + solution are logged instead.

async function example(ca, log) {
  for (let N = 1; N <= 10; ++N) {
    log("N = " + N);

    // Degree of interpolating polynomial
    const K = 2;

    // Legendre roots
    const tau_root = [0.0, 0.211325, 0.788675];

    // Differential equation  dz/dt = z^2 - 2z + 1
    const z = ca.SX.sym("z");
    const F = ca.Function("dz_dt", [z], [ca.plus(ca.minus(ca.times(z, z), ca.times(2, z)), 1)]);

    const z0 = -3;

    // Collocation point and step size
    const tau = ca.SX.sym("tau");
    const h = 1.0 / N;

    // Coefficients of continuity (D) and collocation (C) equations
    const D = new Array(K + 1).fill(0);
    const C = Array.from({ length: K + 1 }, () => new Array(K + 1).fill(0));
    for (let j = 0; j <= K; ++j) {
      let L = ca.SX(1);
      for (let k = 0; k <= K; ++k) {
        if (k !== j) {
          L = ca.times(L, ca.rdivide(ca.minus(tau, tau_root[k]), tau_root[j] - tau_root[k]));
        }
      }
      const lfcn = ca.Function("lfcn", [tau], [L]);
      D[j] = Number(lfcn(1.0));

      const tfcn = ca.Function("tfcn", [tau], [ca.tangent(L, tau)]);
      for (let k = 0; k <= K; ++k) {
        C[j][k] = Number(tfcn(tau_root[k]));
      }
    }

    // Collocated states as an [i][j] grid of scalar symbols
    const Z = Array.from({ length: N }, (_, i) =>
      Array.from({ length: K + 1 }, (_, j) => ca.SX.sym("Z_" + i + "_" + j)));

    // x = vec(Z.T): column-major, i.e. row i then column j
    const xparts = [];
    for (let i = 0; i < N; ++i) for (let j = 0; j <= K; ++j) xparts.push(Z[i][j]);
    const x = ca.vcat(xparts);

    // Construct the NLP constraints
    const g = [];
    for (let i = 0; i < N; ++i) {
      for (let k = 1; k <= K; ++k) {
        let rhs = ca.SX(0);
        for (let j = 0; j <= K; ++j) rhs = ca.plus(rhs, ca.times(Z[i][j], C[j][k]));
        const FF = F(Z[i][k]);
        g.push(ca.minus(ca.times(h, FF), rhs));
      }
      let rhs = ca.SX(0);
      for (let j = 0; j <= K; ++j) rhs = ca.plus(rhs, ca.times(D[j], Z[i][j]));
      if (i < N - 1) g.push(ca.minus(Z[i + 1][0], rhs));
    }
    const gv = ca.vcat(g);

    // NLP: minimize x[0]^2 s.t. collocation/continuity equations
    const nlp = { x: x, f: ca.power(xparts[0], 2), g: gv };

    const solver = ca.nlpsol("solver", "ipopt", nlp, { "ipopt.tol": 1e-10, "ipopt.print_level": 0, print_time: false });

    const n = Number(x.nnz());
    const lbx = new Array(n).fill(-100);
    const ubx = new Array(n).fill(100);
    lbx[0] = ubx[0] = z0;

    const res = solver.call({
      x0: ca.DM.zeros(n),
      lbx: ca.DM(lbx),
      ubx: ca.DM(ubx),
      lbg: ca.DM(0),
      ubg: ca.DM(0),
    });

    log("  optimal cost: " + res["f"]);
    log("  optimal solution: " + res["x"]);
  }
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

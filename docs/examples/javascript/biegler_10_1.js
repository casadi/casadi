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

async function example(M, log) {
  for (let N = 1; N <= 10; ++N) {
    log("N = " + N);

    // Degree of interpolating polynomial
    const K = 2;

    // Legendre roots
    const tau_root = [0.0, 0.211325, 0.788675];

    // Differential equation  dz/dt = z^2 - 2z + 1
    const z = M.SX.sym("z");
    const F = new M.Function("dz_dt", [z], [M.plus(M.minus(M.times(z, z), M.times(2, z)), 1)]);

    const z0 = -3;

    // Collocation point and step size
    const tau = M.SX.sym("tau");
    const h = 1.0 / N;

    // Coefficients of continuity (D) and collocation (C) equations
    const D = new Array(K + 1).fill(0);
    const C = Array.from({ length: K + 1 }, () => new Array(K + 1).fill(0));
    for (let j = 0; j <= K; ++j) {
      let L = M.SX(1);
      for (let k = 0; k <= K; ++k) {
        if (k !== j) {
          L = M.times(L, M.rdivide(M.minus(tau, tau_root[k]), tau_root[j] - tau_root[k]));
        }
      }
      const lfcn = new M.Function("lfcn", [tau], [L]);
      D[j] = lfcn.call([M.DM(1.0)])[0].nonzeros()[0];

      const tfcn = new M.Function("tfcn", [tau], [M.tangent(L, tau)]);
      for (let k = 0; k <= K; ++k) {
        C[j][k] = tfcn.call([M.DM(tau_root[k])])[0].nonzeros()[0];
      }
    }

    // Collocated states as an [i][j] grid of scalar symbols
    const Z = Array.from({ length: N }, (_, i) =>
      Array.from({ length: K + 1 }, (_, j) => M.SX.sym("Z_" + i + "_" + j)));

    // x = vec(Z.T): column-major, i.e. row i then column j
    const xparts = [];
    for (let i = 0; i < N; ++i) for (let j = 0; j <= K; ++j) xparts.push(Z[i][j]);
    const x = M.vertcat.apply(null, xparts);

    // Construct the NLP constraints
    const g = [];
    for (let i = 0; i < N; ++i) {
      for (let k = 1; k <= K; ++k) {
        let rhs = M.SX(0);
        for (let j = 0; j <= K; ++j) rhs = M.plus(rhs, M.times(Z[i][j], C[j][k]));
        const FF = F.call([Z[i][k]])[0];
        g.push(M.minus(M.times(h, FF), rhs));
      }
      let rhs = M.SX(0);
      for (let j = 0; j <= K; ++j) rhs = M.plus(rhs, M.times(D[j], Z[i][j]));
      if (i < N - 1) g.push(M.minus(Z[i + 1][0], rhs));
    }
    const gv = M.vertcat.apply(null, g);

    // NLP: minimize x[0]^2 s.t. collocation/continuity equations
    const nlp = { x: x, f: M.power(xparts[0], 2), g: gv };

    const solver = M.nlpsol("solver", "ipopt", nlp, { "ipopt.tol": 1e-10, "ipopt.print_level": 0, print_time: false });

    const n = Number(x.nnz());
    const lbx = new Array(n).fill(-100);
    const ubx = new Array(n).fill(100);
    lbx[0] = ubx[0] = z0;

    const res = solver.call({
      x0: M.DM.zeros(n),
      lbx: M.DM(lbx),
      ubx: M.DM(ubx),
      lbg: M.DM(0),
      ubg: M.DM(0),
    });

    log("  optimal cost: " + res["f"].nonzeros()[0]);
    log("  optimal solution: [" + res["x"].nonzeros().map((v) => v.toFixed(4)).join(", ") + "]");
  }
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

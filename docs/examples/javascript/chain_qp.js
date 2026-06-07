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
// JS port of docs/examples/python/chain_qp.py.
//
// Hanging chain of N masses joined by N-1 springs; the equilibrium minimises
// the total potential energy subject to piecewise-linear ground constraints.
// Solved as a QP with qpoases.
//
// JS notes (see README.md):
//   * Python's `inf` becomes JS `Infinity`.
//   * Each (y_i, z_i) is a scalar SX symbol kept in JS arrays; the variable
//     vector and constraints are assembled with vertcat.
//   * Plotting is dropped; the optimal cost and a few coordinates are logged.

async function example(M, log) {
  // Constants
  const N = 40;
  const m_i = 40.0 / N;
  const D_i = 70.0 * N;
  const g0 = 9.81;
  const zmin = 0.5; // ground

  let Vchain = M.SX(0);
  const x = [];
  const lbx = [];
  const ubx = [];
  const g = [];
  const lbg = [];
  const ubg = [];

  let y_prev = null, z_prev = null;
  for (let i = 1; i <= N; ++i) {
    const y_i = M.SX.sym("y_" + i);
    const z_i = M.SX.sym("z_" + i);
    x.push(y_i, z_i);

    if (i === 1) { lbx.push(-2.0, 1.0); ubx.push(-2.0, 1.0); }
    else if (i === N) { lbx.push(2.0, 1.0); ubx.push(2.0, 1.0); }
    else { lbx.push(-Infinity, zmin); ubx.push(Infinity, Infinity); }

    if (i > 1) {
      Vchain = M.plus(Vchain,
        M.times(D_i / 2,
          M.plus(M.power(M.minus(y_prev, y_i), 2), M.power(M.minus(z_prev, z_i), 2))));
    }
    Vchain = M.plus(Vchain, M.times(g0 * m_i, z_i));

    // Slanted ground constraint:  z_i - 0.1*y_i >= 0.5
    g.push(M.minus(z_i, M.times(0.1, y_i)));
    lbg.push(0.5);
    ubg.push(Infinity);

    y_prev = y_i; z_prev = z_i;
  }

  // Formulate and solve the QP
  const qp = { x: M.vertcat.apply(null, x), f: Vchain, g: M.vertcat.apply(null, g) };
  const solver = M.qpsol("solver", "qpoases", qp, { sparse: true, printLevel: "none" });

  const sol = solver.call({
    lbx: M.DM(lbx), ubx: M.DM(ubx), lbg: M.DM(lbg), ubg: M.DM(ubg),
  });

  const xopt = sol["x"].nonzeros();
  log("f_opt = " + sol["f"].nonzeros()[0]);

  // Retrieve coordinates (even -> y, odd -> z)
  const Y = xopt.filter((_, k) => k % 2 === 0);
  const Z = xopt.filter((_, k) => k % 2 === 1);
  log("first mass (y, z) = (" + Y[0].toFixed(4) + ", " + Z[0].toFixed(4) + ")");
  const mid = Math.floor(N / 2);
  log("middle mass (y, z) = (" + Y[mid].toFixed(4) + ", " + Z[mid].toFixed(4) + ")");
  log("last mass (y, z) = (" + Y[N - 1].toFixed(4) + ", " + Z[N - 1].toFixed(4) + ")");
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

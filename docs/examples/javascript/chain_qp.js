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

async function example(ca, log) {
  // Constants
  const N = 40;
  const m_i = 40.0 / N;
  const D_i = 70.0 * N;
  const g0 = 9.81;
  const zmin = 0.5; // ground

  let Vchain = ca.SX(0);
  const x = [];
  const lbx = [];
  const ubx = [];
  const g = [];
  const lbg = [];
  const ubg = [];

  let y_prev = null, z_prev = null;
  for (let i = 1; i <= N; ++i) {
    const y_i = ca.SX.sym("y_" + i);
    const z_i = ca.SX.sym("z_" + i);
    x.push(y_i, z_i);

    if (i === 1) { lbx.push(-2.0, 1.0); ubx.push(-2.0, 1.0); }
    else if (i === N) { lbx.push(2.0, 1.0); ubx.push(2.0, 1.0); }
    else { lbx.push(-Infinity, zmin); ubx.push(Infinity, Infinity); }

    if (i > 1) {
      Vchain = ca.plus(Vchain,
        ca.times(D_i / 2,
          ca.plus(ca.power(ca.minus(y_prev, y_i), 2), ca.power(ca.minus(z_prev, z_i), 2))));
    }
    Vchain = ca.plus(Vchain, ca.times(g0 * m_i, z_i));

    // Slanted ground constraint:  z_i - 0.1*y_i >= 0.5
    g.push(ca.minus(z_i, ca.times(0.1, y_i)));
    lbg.push(0.5);
    ubg.push(Infinity);

    y_prev = y_i; z_prev = z_i;
  }

  // Formulate and solve the QP
  const qp = { x: ca.vcat(x), f: Vchain, g: ca.vcat(g) };
  const solver = ca.qpsol("solver", "qpoases", qp, { sparse: true, printLevel: "none" });

  const sol = solver.call({
    lbx: ca.DM(lbx), ubx: ca.DM(ubx), lbg: ca.DM(lbg), ubg: ca.DM(ubg),
  });

  log("f_opt = " + sol["f"]);

  // Retrieve coordinates: column k of the 2xN reshape is (y_k, z_k)
  const YZ = ca.reshape(sol["x"], 2, N);
  const mid = Math.floor(N / 2);
  log("first mass (y, z) = " + YZ[":,0"]);
  log("middle mass (y, z) = " + YZ[":," + mid]);
  log("last mass (y, z) = " + YZ[":," + (N - 1)]);
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

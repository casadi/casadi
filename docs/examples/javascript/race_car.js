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
// JS port of docs/examples/python/race_car.py.
//
// Car race along a track: a minimal-time OCP solved with direct
// multiple-shooting via the Opti stack (fully portable).
//
// JS notes (see README.md):
//   * Column slices X[:,k] come from M.horzsplit(X); rows from M.vertsplit.
//   * matplotlib output is dropped; we log the optimum instead.

async function example(M, log) {
  const pi = Math.PI;
  const N = 100; // number of control intervals

  const opti = new M.Opti();

  // ---- decision variables ---------
  const X = opti.variable(2, N + 1);      // state trajectory
  const Xc = M.horzsplit(X);              // columns X[:,k]
  const [posRow, speedRow] = M.vertsplit(X); // rows
  const U = opti.variable(1, N);          // control trajectory (throttle)
  const Uc = M.horzsplit(U);
  const T = opti.variable();              // final time

  // ---- objective ---------
  opti.minimize(T); // race in minimal time

  // ---- dynamic constraints --------
  // dx/dt = f(x,u) = [x1, u - x1]
  const f = (x, u) => {
    const [, x1] = M.vertsplit(x);
    return M.vertcat(x1, M.minus(u, x1));
  };

  const dt = M.rdivide(T, N); // length of a control interval
  for (let k = 0; k < N; k++) {
    const xk = Xc[k], uk = Uc[k];
    // Runge-Kutta 4 integration
    const k1 = f(xk, uk);
    const k2 = f(M.plus(xk, M.times(M.rdivide(dt, 2), k1)), uk);
    const k3 = f(M.plus(xk, M.times(M.rdivide(dt, 2), k2)), uk);
    const k4 = f(M.plus(xk, M.times(dt, k3)), uk);
    const sum = M.plus(M.plus(k1, M.times(2, k2)), M.plus(M.times(2, k3), k4));
    const x_next = M.plus(xk, M.times(M.rdivide(dt, 6), sum));
    opti.subject_to(M.eq(Xc[k + 1], x_next)); // close the gaps
  }

  // ---- path constraints -----------
  const limit = (p) => M.minus(1, M.rdivide(M.sin(M.times(2 * pi, p)), 2));
  opti.subject_to(M.le(speedRow, limit(posRow)));   // track speed limit
  opti.subject_to(M.le(0, U));                      // control is limited ...
  opti.subject_to(M.le(U, 1));                      // ... to [0, 1]

  // ---- boundary conditions --------
  const posSplit = M.horzsplit(posRow); // scalars along the row
  const speedSplit = M.horzsplit(speedRow);
  opti.subject_to(M.eq(posSplit[0], 0));      // start at position 0 ...
  opti.subject_to(M.eq(speedSplit[0], 0));    // ... from stand-still
  opti.subject_to(M.eq(posSplit[N], 1));      // finish line at position 1

  // ---- misc. constraints ----------
  opti.subject_to(M.ge(T, 0)); // Time must be positive

  // ---- initial values for solver ---
  opti.set_initial(speedRow, 1);
  opti.set_initial(T, 1);

  // ---- solve NLP ------
  opti.solver("ipopt");
  const sol = opti.solve();

  // ---- post-processing ------
  log("-----");
  log("minimal time T = " + sol.value(T).nonzeros().join(" "));
  log("speed = " + sol.value(speedRow).nonzeros().map((v) => v.toFixed(4)).join(" "));
  log("pos   = " + sol.value(posRow).nonzeros().map((v) => v.toFixed(4)).join(" "));
  log("throttle U = " + sol.value(U).nonzeros().map((v) => v.toFixed(4)).join(" "));
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

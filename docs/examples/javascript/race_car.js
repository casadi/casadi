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
//   * Column slices X[:,k] come from ca.horzsplit(X); rows from ca.vertsplit.
//   * matplotlib output is dropped; we log the optimum instead.

async function example(ca, log) {
  const pi = Math.PI;
  const N = 100; // number of control intervals

  const opti = ca.Opti();

  // ---- decision variables ---------
  const X = opti.variable(2, N + 1);      // state trajectory
  const Xc = ca.horzsplit(X);              // columns X[:,k]
  const [posRow, speedRow] = ca.vertsplit(X); // rows
  const U = opti.variable(1, N);          // control trajectory (throttle)
  const Uc = ca.horzsplit(U);
  const T = opti.variable();              // final time

  // ---- objective ---------
  opti.minimize(T); // race in minimal time

  // ---- dynamic constraints --------
  // dx/dt = f(x,u) = [x1, u - x1]
  const f = (x, u) => {
    const [, x1] = ca.vertsplit(x);
    return ca.vertcat(x1, ca.minus(u, x1));
  };

  const dt = ca.rdivide(T, N); // length of a control interval
  for (let k = 0; k < N; k++) {
    const xk = Xc[k], uk = Uc[k];
    // Runge-Kutta 4 integration
    const k1 = f(xk, uk);
    const k2 = f(ca.plus(xk, ca.times(ca.rdivide(dt, 2), k1)), uk);
    const k3 = f(ca.plus(xk, ca.times(ca.rdivide(dt, 2), k2)), uk);
    const k4 = f(ca.plus(xk, ca.times(dt, k3)), uk);
    const sum = ca.plus(ca.plus(k1, ca.times(2, k2)), ca.plus(ca.times(2, k3), k4));
    const x_next = ca.plus(xk, ca.times(ca.rdivide(dt, 6), sum));
    opti.subject_to(ca.eq(Xc[k + 1], x_next)); // close the gaps
  }

  // ---- path constraints -----------
  const limit = (p) => ca.minus(1, ca.rdivide(ca.sin(ca.times(2 * pi, p)), 2));
  opti.subject_to(ca.le(speedRow, limit(posRow)));   // track speed limit
  opti.subject_to(opti.bounded(0, U, 1));            // control is limited

  // ---- boundary conditions --------
  const posSplit = ca.horzsplit(posRow); // scalars along the row
  const speedSplit = ca.horzsplit(speedRow);
  opti.subject_to(ca.eq(posSplit[0], 0));      // start at position 0 ...
  opti.subject_to(ca.eq(speedSplit[0], 0));    // ... from stand-still
  opti.subject_to(ca.eq(posSplit[N], 1));      // finish line at position 1

  // ---- misc. constraints ----------
  opti.subject_to(ca.ge(T, 0)); // Time must be positive

  // ---- initial values for solver ---
  opti.set_initial(speedRow, 1);
  opti.set_initial(T, 1);

  // ---- solve NLP ------
  opti.solver("ipopt");
  const sol = opti.solve();

  // ---- post-processing ------
  log("-----");
  log("minimal time T = " + sol.value(T));
  log("speed = " + sol.value(speedRow));
  log("pos   = " + sol.value(posRow));
  log("throttle U = " + sol.value(U));
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

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
// JS port of docs/examples/python/direct_multiple_shooting.py.
//
// Direct multiple shooting on a Van der Pol OCP. The Python version offers a
// CVODES branch (disabled by `if False`); we keep the fixed-step RK4 branch,
// which is the active one and fully portable.
//
// JS notes (see README.md): matplotlib output is dropped; key numbers logged.

async function example(ca, log) {
  const T = 10.0;  // Time horizon
  const N = 20;    // number of control intervals

  // Declare model variables
  const x1 = ca.MX.sym("x1");
  const x2 = ca.MX.sym("x2");
  const x = ca.vertcat(x1, x2);
  const u = ca.MX.sym("u");

  // Model equations
  const xdot = ca.vertcat(
    ca.plus(ca.minus(ca.times(ca.minus(1, ca.times(x2, x2)), x1), x2), u),
    x1);

  // Objective term
  const L = ca.plus(ca.plus(ca.times(x1, x1), ca.times(x2, x2)), ca.times(u, u));

  // Formulate discrete time dynamics: fixed step Runge-Kutta 4 integrator
  const Msteps = 4;            // RK4 steps per interval
  const DT = T / N / Msteps;
  const f = ca.Function("f", [x, u], [xdot, L]);
  const X0 = ca.MX.sym("X0", 2);
  const U = ca.MX.sym("U");
  let X = X0;
  let Q = ca.MX(0);
  for (let j = 0; j < Msteps; j++) {
    const [k1, k1q] = f(X, U);
    const [k2, k2q] = f(ca.plus(X, ca.times(DT / 2, k1)), U);
    const [k3, k3q] = f(ca.plus(X, ca.times(DT / 2, k2)), U);
    const [k4, k4q] = f(ca.plus(X, ca.times(DT, k3)), U);
    X = ca.plus(X, ca.times(DT / 6, ca.plus(ca.plus(k1, ca.times(2, k2)), ca.plus(ca.times(2, k3), k4))));
    Q = ca.plus(Q, ca.times(DT / 6, ca.plus(ca.plus(k1q, ca.times(2, k2q)), ca.plus(ca.times(2, k3q), k4q))));
  }
  const F = ca.Function("F", [X0, U], [X, Q], ["x0", "p"], ["xf", "qf"]);

  // Evaluate at a test point
  const Fk0 = F.call({ x0: ca.DM([0.2, 0.3]), p: ca.DM(0.4) });
  log("test xf = " + Fk0["xf"]);
  log("test qf = " + Fk0["qf"]);

  // Start with an empty NLP
  const w = [], w0 = [], lbw = [], ubw = [];
  let J = ca.MX(0);
  const g = [], lbg = [], ubg = [];

  // "Lift" initial conditions
  let Xk = ca.MX.sym("X0", 2);
  w.push(Xk); lbw.push(0, 1); ubw.push(0, 1); w0.push(0, 1);

  // Formulate the NLP
  for (let k = 0; k < N; k++) {
    // New NLP variable for the control
    const Uk = ca.MX.sym("U_" + k);
    w.push(Uk); lbw.push(-1); ubw.push(1); w0.push(0);

    // Integrate till the end of the interval
    const Fk = F.call({ x0: Xk, p: Uk });
    const Xk_end = Fk["xf"];
    J = ca.plus(J, Fk["qf"]);

    // New NLP variable for state at end of interval
    Xk = ca.MX.sym("X_" + (k + 1), 2);
    w.push(Xk); lbw.push(-0.25, -Infinity); ubw.push(Infinity, Infinity); w0.push(0, 0);

    // Add equality constraint
    g.push(ca.minus(Xk_end, Xk)); lbg.push(0, 0); ubg.push(0, 0);
  }

  // Create an NLP solver
  const prob = { f: J, x: ca.vcat(w), g: ca.vcat(g) };
  const solver = ca.nlpsol("solver", "ipopt", prob);

  // Solve the NLP
  const sol = solver.call({
    x0: ca.DM(w0), lbx: ca.DM(lbw), ubx: ca.DM(ubw), lbg: ca.DM(lbg), ubg: ca.DM(ubg),
  });
  const w_opt = sol["x"];

  // Extract the solution: layout is [x1, x2, u] repeating
  const x1_opt = w_opt["0::3"];
  const x2_opt = w_opt["1::3"];
  const u_opt = w_opt["2::3"];

  log("-----");
  log("objective at solution = " + sol["f"]);
  log("x1_opt = " + x1_opt);
  log("x2_opt = " + x2_opt);
  log("u_opt = " + u_opt);
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

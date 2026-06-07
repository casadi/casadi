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
// JS port of docs/examples/python/direct_single_shooting.py.
//
// Direct single shooting on a Van der Pol OCP. The Python version offers a
// CVODES branch (disabled by `if False`); we keep the fixed-step RK4 branch,
// which is the active one and fully portable.
//
// JS notes (see README.md): matplotlib output is dropped; key numbers logged.

async function example(M, log) {
  const T = 10.0;  // Time horizon
  const N = 20;    // number of control intervals

  // Declare model variables
  const x1 = M.MX.sym("x1");
  const x2 = M.MX.sym("x2");
  const x = M.vertcat(x1, x2);
  const u = M.MX.sym("u");

  // Model equations
  const xdot = M.vertcat(
    M.plus(M.minus(M.times(M.minus(1, M.times(x2, x2)), x1), x2), u),
    x1);

  // Objective term
  const L = M.plus(M.plus(M.times(x1, x1), M.times(x2, x2)), M.times(u, u));

  // Formulate discrete time dynamics: fixed step Runge-Kutta 4 integrator
  const Msteps = 4;            // RK4 steps per interval
  const DT = T / N / Msteps;
  const f = new M.Function("f", [x, u], [xdot, L]);
  const X0 = M.MX.sym("X0", 2);
  const U = M.MX.sym("U");
  let X = X0;
  let Q = M.MX(0);
  for (let j = 0; j < Msteps; j++) {
    const [k1, k1q] = f.call([X, U]);
    const [k2, k2q] = f.call([M.plus(X, M.times(DT / 2, k1)), U]);
    const [k3, k3q] = f.call([M.plus(X, M.times(DT / 2, k2)), U]);
    const [k4, k4q] = f.call([M.plus(X, M.times(DT, k3)), U]);
    X = M.plus(X, M.times(DT / 6, M.plus(M.plus(k1, M.times(2, k2)), M.plus(M.times(2, k3), k4))));
    Q = M.plus(Q, M.times(DT / 6, M.plus(M.plus(k1q, M.times(2, k2q)), M.plus(M.times(2, k3q), k4q))));
  }
  const F = new M.Function("F", [X0, U], [X, Q], ["x0", "p"], ["xf", "qf"]);

  // Evaluate at a test point
  const Fk0 = F.call({ x0: M.DM([0.2, 0.3]), p: M.DM(0.4) });
  log("test xf = " + Fk0["xf"].nonzeros().join(" "));
  log("test qf = " + Fk0["qf"].nonzeros().join(" "));

  // Start with an empty NLP
  const w = [], w0 = [], lbw = [], ubw = [];
  let J = M.MX(0);
  const g = [], lbg = [], ubg = [];

  // Formulate the NLP
  let Xk = M.MX([0, 1]);
  for (let k = 0; k < N; k++) {
    // New NLP variable for the control
    const Uk = M.MX.sym("U_" + k);
    w.push(Uk); lbw.push(-1); ubw.push(1); w0.push(0);

    // Integrate till the end of the interval
    const Fk = F.call({ x0: Xk, p: Uk });
    Xk = Fk["xf"];
    J = M.plus(J, Fk["qf"]);

    // Add inequality constraint
    const [Xk0] = M.vertsplit(Xk);
    g.push(Xk0); lbg.push(-0.25); ubg.push(Infinity);
  }

  // Create an NLP solver
  const prob = { f: J, x: M.vcat(w), g: M.vcat(g) };
  const solver = M.nlpsol("solver", "ipopt", prob);

  // Solve the NLP
  const sol = solver.call({
    x0: M.DM(w0), lbx: M.DM(lbw), ubx: M.DM(ubw), lbg: M.DM(lbg), ubg: M.DM(ubg),
  });
  const u_opt = sol["x"].nonzeros();

  // Re-simulate to recover the state trajectory
  let xk = [0, 1];
  const x1_opt = [xk[0]], x2_opt = [xk[1]];
  for (let k = 0; k < N; k++) {
    const Fk = F.call({ x0: M.DM(xk), p: M.DM(u_opt[k]) });
    xk = Fk["xf"].nonzeros();
    x1_opt.push(xk[0]); x2_opt.push(xk[1]);
  }

  log("-----");
  log("objective at solution = " + sol["f"].nonzeros().join(" "));
  log("u_opt = " + u_opt.map((v) => v.toFixed(4)).join(" "));
  log("x1_opt = " + x1_opt.map((v) => v.toFixed(4)).join(" "));
  log("x2_opt = " + x2_opt.map((v) => v.toFixed(4)).join(" "));
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

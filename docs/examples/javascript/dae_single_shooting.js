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
// JS port of docs/examples/python/dae_single_shooting.py.
//
// Direct single shooting on a small DAE OCP.
//
//   minimize   integral_{0}^{10} x0^2 + x1^2 + u^2 dt
//   s.t.       dot(x0) == z*x0 - x1 + u
//              dot(x1) == x0
//                    0 == x1^2 + z - 1
//              x0(0)=0, x1(0)=1, x0(10)=x1(10)=0,  -0.75 <= u <= 1
//
// JS notes (see README.md):
//   * The Python version uses the "idas" DAE integrator, which aborts the
//     wasm runtime. We substitute "collocation", which handles the implicit
//     algebraic equation and produces an equivalent solution.
//   * The Python `ma27` linear solver option is dropped (default mumps used).

async function example(ca, log) {
  // Declare variables
  const x = ca.SX.sym("x", 2);   // Differential states
  const z = ca.SX.sym("z");      // Algebraic variable
  const u = ca.SX.sym("u");      // Control
  const [x0e, x1e] = ca.vertsplit(x);

  // Differential equation
  const f_x = ca.vertcat(ca.plus(ca.minus(ca.times(z, x0e), x1e), u), x0e);

  // Algebraic equation
  const f_z = ca.minus(ca.plus(ca.times(x1e, x1e), z), 1);

  // Lagrange cost term (quadrature)
  const f_q = ca.plus(ca.plus(ca.times(x0e, x0e), ca.times(x1e, x1e)), ca.times(u, u));

  // Create an integrator (interval length 0.5s) -- "collocation" stands in for "idas"
  const dae = { x: x, z: z, p: u, ode: f_x, alg: f_z, quad: f_q };
  const I = ca.integrator("I", "collocation", dae, 0, 0.5);

  // All controls
  const U = ca.MX.sym("U", 20);
  const Us = ca.vertsplit(U);

  // Construct graph of integrator calls
  let X = ca.MX([0, 1]);
  let J = ca.MX(0);
  for (let k = 0; k < 20; k++) {
    const Ik = I.call({ x0: X, p: Us[k] });
    X = Ik["xf"];
    J = ca.plus(J, Ik["qf"]);   // Sum up quadratures
  }

  // Allocate an NLP solver
  const nlp = { x: U, f: J, g: X };
  const solver = ca.nlpsol("solver", "ipopt", nlp);

  // Pass bounds, initial guess and solve NLP
  const sol = solver.call({
    lbx: -0.75,   // Lower variable bound
    ubx: 1.0,     // Upper variable bound
    lbg: 0.0,     // Lower constraint bound
    ubg: 0.0,     // Upper constraint bound
    x0: 0.0,      // Initial guess
  });

  log("-----");
  log("objective at solution = " + sol["f"]);
  log("u_opt = " + sol["x"]);
  log("terminal state x(10) = " + sol["g"]);
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

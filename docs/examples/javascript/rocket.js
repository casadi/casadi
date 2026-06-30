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
// JS port of docs/examples/python/rocket.py.
//
// Rocket trajectory NLP: discretize the ODE with explicit Euler and let
// ipopt pick the throttle for each of 50 segments (fully portable).
//
// JS notes (see README.md):
//   * matplotlib output is dropped; we log the optimum instead.
//   * The Python `expand=True` option is dropped: under the wasm build it
//     crashes the GC teardown. The solution is identical without it.

async function example(ca, log) {
  // Control
  const u = ca.MX.sym("u");

  // State
  const x = ca.MX.sym("x", 3);
  const [s, v, m] = ca.vertsplit(x); // position, speed, mass

  // ODE right hand side
  const sdot = v;
  const vdot = ca.rdivide(ca.minus(u, ca.times(0.05, ca.times(v, v))), m);
  const mdot = ca.times(-0.1, ca.times(u, u));
  const xdot = ca.vertcat(sdot, vdot, mdot);

  // ODE right hand side function
  const f = ca.Function("f", [x, u], [xdot]);

  // Integrate with Explicit Euler over 0.2 seconds
  const dt = 0.01; // Time step
  let xj = x;
  for (let j = 0; j < 20; j++) {
    const fj = f(xj, u);
    xj = ca.plus(xj, ca.times(dt, fj));
  }

  // Discrete time dynamics function
  const F = ca.Function("F", [x, u], [xj]);

  // Number of control segments
  const nu = 50;

  // Control for all segments
  const U = ca.MX.sym("U", nu);
  const Us = ca.vertsplit(U);

  // Initial conditions
  const X0 = ca.MX([0, 0, 1]);

  // Integrate over all intervals
  let X = X0;
  for (let k = 0; k < nu; k++) X = F(X, Us[k]);

  // Objective function and constraints
  const J = ca.dot(U, U);                 // u'*u
  const G = ca.vcat(ca.vertsplit(X).slice(0, 2)); // x(0:2)

  // NLP
  const nlp = { x: U, f: J, g: G };

  // Allocate an NLP solver
  const solver = ca.nlpsol("solver", "ipopt", nlp, { "ipopt.tol": 1e-10 });

  // Bounds on u, initial condition, and bounds on g; solve
  const res = solver.call({
    lbx: -0.5, ubx: 0.5, x0: 0.4,
    lbg: ca.DM([10, 0]), ubg: ca.DM([10, 0]),
  });

  log("-----");
  log("objective at solution = " + res["f"]);
  log("u_opt = " + res["x"]);
  log("lam_x = " + res["lam_x"]);
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

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

async function example(M, log) {
  // Control
  const u = M.MX.sym("u");

  // State
  const x = M.MX.sym("x", 3);
  const [s, v, m] = M.vertsplit(x); // position, speed, mass

  // ODE right hand side
  const sdot = v;
  const vdot = M.rdivide(M.minus(u, M.times(0.05, M.times(v, v))), m);
  const mdot = M.times(-0.1, M.times(u, u));
  const xdot = M.vertcat(sdot, vdot, mdot);

  // ODE right hand side function
  const f = new M.Function("f", [x, u], [xdot]);

  // Integrate with Explicit Euler over 0.2 seconds
  const dt = 0.01; // Time step
  let xj = x;
  for (let j = 0; j < 20; j++) {
    const fj = f.call([xj, u])[0];
    xj = M.plus(xj, M.times(dt, fj));
  }

  // Discrete time dynamics function
  const F = new M.Function("F", [x, u], [xj]);

  // Number of control segments
  const nu = 50;

  // Control for all segments
  const U = M.MX.sym("U", nu);
  const Us = M.vertsplit(U);

  // Initial conditions
  const X0 = M.MX([0, 0, 1]);

  // Integrate over all intervals
  let X = X0;
  for (let k = 0; k < nu; k++) X = F.call([X, Us[k]])[0];

  // Objective function and constraints
  const J = M.dot(U, U);                 // u'*u
  const G = M.vertcat(...M.vertsplit(X).slice(0, 2)); // x(0:2)

  // NLP
  const nlp = { x: U, f: J, g: G };

  // Allocate an NLP solver
  const solver = M.nlpsol("solver", "ipopt", nlp, { "ipopt.tol": 1e-10 });

  // Bounds on u, initial condition, and bounds on g; solve
  const res = solver.call({
    lbx: -0.5, ubx: 0.5, x0: 0.4,
    lbg: M.DM([10, 0]), ubg: M.DM([10, 0]),
  });

  log("-----");
  log("objective at solution = " + res["f"].nonzeros().join(" "));
  log("u_opt = " + res["x"].nonzeros().map((v) => v.toFixed(4)).join(" "));
  log("lam_x = " + res["lam_x"].nonzeros().map((v) => v.toFixed(4)).join(" "));
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

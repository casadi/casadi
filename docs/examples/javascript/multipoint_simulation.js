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
// JS port of docs/examples/python/multipoint_simulation.py.
//
// Simulate a Van der Pol ODE over a grid, accumulating the integral of
// the sum-of-squares output via the integrator's quadrature.
//
// JS notes (see README.md):
//   * cvodes is not available in the wasm build; we substitute the "rk"
//     (fixed-step Runge-Kutta) integrator, which handles this ODE fine.
//   * numpy/matplotlib output is dropped; we log the trajectories.

async function example(M, log) {
  // Nonlinear system with two states and one control input
  const x1 = M.MX.sym("x1");
  const x2 = M.MX.sym("x2");
  const x = M.vertcat(x1, x2);
  const u = M.MX.sym("u");
  const xdot = M.vertcat(M.plus(M.minus(M.times(M.minus(1, M.times(x2, x2)), x1), x2), u), x1);

  // Sum-of-squares distance from the origin
  const y = M.plus(M.times(x1, x1), M.times(x2, x2));

  // Time horizon, discretization
  const N = 10;
  const T = 10.0;
  const tgrid = [];
  for (let k = 0; k <= N; k++) tgrid.push((T * k) / N);

  // Integrate, also calculating the integral of y (substitute rk for cvodes)
  const dae = { x: x, u: u, ode: xdot, quad: y };
  const F = M.integrator("F", "rk", dae, tgrid[0], tgrid.slice(1), { simplify: true });

  // Pointwise values of y
  const yfun = new M.Function("yfun", [x], [y], ["x"], ["y"]);

  // Initial conditions and piecewise-constant controls
  const uvals = [];
  for (let k = 0; k < N; k++) uvals.push(-1 + (2 * k) / (N - 1));
  const x0 = M.DM([0, 0]);

  // Simulate
  const Fk = F.call({ x0: x0, u: M.DM([uvals]) });
  let xf = Fk["xf"];
  let qf = Fk["qf"];

  // Prepend the state at the initial time
  xf = M.horzcat(x0, xf);
  qf = M.horzcat(M.DM(0), qf);

  // y at each grid point
  const yf = yfun.call([xf])[0];

  const x_sim = xf.nonzeros();         // column-major: [x1_0,x2_0,x1_1,x2_1,...]
  const y_sim = yf.nonzeros();
  const q_sim = qf.nonzeros();

  log("-----");
  log("x1 = " + x_sim.filter((_, i) => i % 2 === 0).map((v) => v.toFixed(4)).join(" "));
  log("x2 = " + x_sim.filter((_, i) => i % 2 === 1).map((v) => v.toFixed(4)).join(" "));
  log("y  = " + y_sim.map((v) => v.toFixed(4)).join(" "));
  log("integral of y (quadrature) = " + q_sim.map((v) => v.toFixed(4)).join(" "));

  // Compare integral with a Riemann sum of the pointwise y values
  let acc = 0;
  const cum = y_sim.map((v) => (acc += (T / N) * v));
  log("integral of y (cumsum)     = " + cum.map((v) => v.toFixed(4)).join(" "));
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

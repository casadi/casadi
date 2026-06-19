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

async function example(ca, log) {
  // Nonlinear system with two states and one control input
  const x1 = ca.MX.sym("x1");
  const x2 = ca.MX.sym("x2");
  const x = ca.vertcat(x1, x2);
  const u = ca.MX.sym("u");
  const xdot = ca.vertcat(ca.plus(ca.minus(ca.times(ca.minus(1, ca.times(x2, x2)), x1), x2), u), x1);

  // Sum-of-squares distance from the origin
  const y = ca.plus(ca.times(x1, x1), ca.times(x2, x2));

  // Time horizon, discretization
  const N = 10;
  const T = 10.0;
  const tgrid = [];
  for (let k = 0; k <= N; k++) tgrid.push((T * k) / N);

  // Integrate, also calculating the integral of y (substitute rk for cvodes)
  const dae = { x: x, u: u, ode: xdot, quad: y };
  const F = ca.integrator("F", "rk", dae, tgrid[0], tgrid.slice(1), { simplify: true });

  // Pointwise values of y
  const yfun = ca.Function("yfun", [x], [y], ["x"], ["y"]);

  // Initial conditions and piecewise-constant controls
  const uvals = [];
  for (let k = 0; k < N; k++) uvals.push(-1 + (2 * k) / (N - 1));
  const x0 = ca.DM([0, 0]);

  // Simulate
  const Fk = F.call({ x0: x0, u: ca.DM([uvals]) });
  let xf = Fk["xf"];
  let qf = Fk["qf"];

  // Prepend the state at the initial time
  xf = ca.horzcat(x0, xf);
  qf = ca.horzcat(ca.DM(0), qf);

  // y at each grid point
  const yf = yfun(xf);

  log("-----");
  log("x1 = " + xf["0,:"]);
  log("x2 = " + xf["1,:"]);
  log("y  = " + yf);
  log("integral of y (quadrature) = " + qf);

  // Compare integral with a Riemann sum of the pointwise y values
  log("integral of y (cumsum)     = " + ca.cumsum(ca.times(T / N, yf)));
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

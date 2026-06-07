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
// JS port of docs/examples/python/sensitivity_analysis.py.
//
// Forward / adjoint / second-order integrator sensitivities obtained with
// Function.factory and compared against finite differences.
//
// JS notes (see README.md):
//   * The Python script also loops over a DAE example with the "idas" and
//     "cvodes" integrators and a "kinsol" root-finder; all three abort the
//     wasm runtime, so only the ODE (rocket) example is ported here, run
//     with the "rk" and "collocation" integrators (which support sensitivity
//     factories just fine).

async function example(M, log) {
  log("Testing sensitivity analysis in CasADi");
  log("******");
  log("Testing ODE example");

  // Time
  const t = M.SX.sym("t");
  // Parameter
  const u = M.SX.sym("u");
  // Differential states
  const s = M.SX.sym("s"), v = M.SX.sym("v"), m = M.SX.sym("m");
  const x = M.vertcat(s, v, m);

  // Constants
  const alpha = 0.05; // friction
  const beta = 0.1;   // fuel consumption rate

  // Differential equation
  const ode = M.vertcat(
    v,
    M.rdivide(M.minus(u, M.times(alpha, M.times(v, v))), m),
    M.times(-beta, M.times(u, u)));

  // Quadrature
  const quad = M.plus(M.power(v, 3), M.power(M.minus(M.minus(3, M.sin(t)), u), 2));

  // DAE callback
  const dae = { t: t, x: x, p: u, ode: ode, quad: quad };

  const tf = 0.5;        // Time length
  const x0 = M.DM([0, 0, 1]); // Initial position
  const u0 = 0.4;        // Parameter
  const h = 0.001;       // FD step

  const fmt = (dm) => dm.nonzeros().map((z) => z.toPrecision(8)).join(", ");

  for (const MyIntegrator of ["rk", "collocation"]) {
    log("========");
    log("Integrator: " + MyIntegrator);
    log("========");

    // Integrator
    const I = M.integrator("I", MyIntegrator, dae, 0, tf, { number_of_finite_elements: 100 });

    // Integrate to get results
    let res = I.call({ x0: x0, p: u0 });
    const xf = res["xf"], qf = res["qf"];
    log("Unperturbed solution: xf = [" + fmt(xf) + "], qf = [" + fmt(qf) + "]");

    // Finite-difference approximation
    res = I.call({ x0: x0, p: u0 + h });
    const fd_xf = M.rdivide(M.minus(res["xf"], xf), h);
    const fd_qf = M.rdivide(M.minus(res["qf"], qf), h);
    log("Finite difference: d(xf)/d(p) = [" + fmt(fd_xf) + "], d(qf)/d(p) = [" + fmt(fd_qf) + "]");

    // Forward sensitivities
    const I_fwd = I.factory("I_fwd", ["x0", "z0", "p", "fwd:p"], ["fwd:xf", "fwd:qf"]);
    res = I_fwd.call({ x0: x0, p: u0, fwd_p: 1 });
    log("Forward: d(xf)/d(p) = [" + fmt(res["fwd_xf"]) + "], d(qf)/d(p) = [" + fmt(res["fwd_qf"]) + "]");

    // Adjoint sensitivities
    const I_adj = I.factory("I_adj", ["x0", "z0", "p", "adj:qf"], ["adj:x0", "adj:p"]);
    res = I_adj.call({ x0: x0, p: u0, adj_qf: 1 });
    let adj_x0 = res["adj_x0"], adj_p = res["adj_p"];
    log("Adjoint: d(qf)/d(x0) = [" + fmt(adj_x0) + "], d(qf)/d(p) = [" + fmt(adj_p) + "]");

    // FD of adjoint -> second order
    res = I_adj.call({ x0: x0, p: u0 + h, adj_qf: 1 });
    const fd_adj_x0 = M.rdivide(M.minus(res["adj_x0"], adj_x0), h);
    const fd_adj_p = M.rdivide(M.minus(res["adj_p"], adj_p), h);
    log("FD of adjoint: d2(qf)/d(x0)d(p) = [" + fmt(fd_adj_x0) + "], d2(qf)/d(p)d(p) = [" + fmt(fd_adj_p) + "]");

    // Forward over adjoint -> second order
    const I_foa = I_adj.factory("I_foa", ["x0", "z0", "p", "adj_qf", "fwd:p"], ["fwd:adj_x0", "fwd:adj_p"]);
    res = I_foa.call({ x0: x0, p: u0, adj_qf: 1, fwd_p: 1 });
    log("Forward over adjoint: d2(qf)/d(x0)d(p) = [" + fmt(res["fwd_adj_x0"]) + "], d2(qf)/d(p)d(p) = [" + fmt(res["fwd_adj_p"]) + "]");

    // Adjoint over adjoint -> second order
    const I_aoa = I_adj.factory("I_aoa", ["x0", "z0", "p", "adj_qf", "adj:adj_p"], ["adj:x0", "adj:p"]);
    res = I_aoa.call({ x0: x0, p: u0, adj_qf: 1, adj_adj_p: 1 });
    log("Adjoint over adjoint: d2(qf)/d(x0)d(p) = [" + fmt(res["adj_x0"]) + "], d2(qf)/d(p)d(p) = [" + fmt(res["adj_p"]) + "]");
  }
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

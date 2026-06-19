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
// JS port of docs/examples/python/vdp_indirect_multiple_shooting.py.
//
// Indirect multiple shooting for the Van der Pol OCP: derive the optimal
// control from the Hamiltonian, then root-find the node states of the
// augmented (state, costate) ODE so all shooting gaps close and the
// state/costate boundary conditions hold.
//
// JS notes (see README.md):
//   * The Python "cvodes" integrator aborts the wasm runtime; we substitute
//     "rk". Because each shooting interval is short (tf/num_nodes = 0.5s),
//     explicit RK integrates the augmented ODE stably here.
//   * The root-finding problem uses the "nlpsol"/"ipopt" plugin exactly as
//     the Python original (its globalization is needed; the plain "newton"
//     root-finder lands on spurious roots, and "kinsol" aborts the runtime).
//   * For the trajectory plot the Python code re-integrates the full 10s in
//     one shot; under explicit RK that blows up on the unstable costate ODE,
//     so we instead re-integrate per shooting interval (each stable).
//   * matplotlib output is dropped; we log the optimum instead.

async function example(ca, log) {
  // Declare variables
  const x0 = ca.SX.sym("x0"), x1 = ca.SX.sym("x1");
  const x = ca.vertcat(x0, x1);

  // Control
  const u = ca.SX.sym("u");

  // ODE right hand side
  const xdot = ca.vertcat(
    ca.plus(ca.minus(ca.times(ca.minus(1, ca.times(x1, x1)), x0), x1), u),
    x0);

  // Lagrangian
  const L = ca.plus(ca.plus(ca.times(x0, x0), ca.times(x1, x1)), ca.times(u, u));

  // Costate
  const lam = ca.SX.sym("lam", 2);

  // Hamiltonian
  const H = ca.plus(ca.dot(lam, xdot), L);

  // Costate equations
  const ldot = ca.times(-1, ca.gradient(H, x));

  log("Hamiltonian: " + H.str());

  // H is convex quadratic in u: H = u*u + p*u + q; extract p
  let p = ca.gradient(H, u);
  p = ca.substitute(p, u, ca.SX(0));

  // Unconstrained minimizer u = -p/2, clipped to [-0.75, 1.0]
  let u_opt = ca.times(-0.5, p);
  u_opt = ca.fmin(u_opt, ca.SX(1.0));
  u_opt = ca.fmax(u_opt, ca.SX(-0.75));
  log("optimal control: " + u_opt.str());

  // Augment f with ldot and substitute the optimal control
  let f = ca.vertcat(xdot, ldot);
  f = ca.substitute(f, u, u_opt);

  // Function for the optimal control given the augmented state
  const xlam = ca.vertcat(x, lam);
  const u_fcn = ca.Function("ufcn", [xlam], [u_opt]);

  // DAE for the augmented dynamics
  const dae = { x: xlam, ode: f };

  const nX = 4;              // augmented state dimension
  const tf = 10.0;          // end time
  const num_nodes = 20;     // shooting nodes

  // Integrator over one interval -- "rk" stands in for "cvodes"
  const I = ca.integrator("I", "rk", dae, 0, tf / num_nodes, { number_of_finite_elements: 50 });

  // States at each shooting node (4 x (num_nodes+1))
  const X = ca.MX.sym("X", nX, num_nodes + 1);
  const Xc = ca.horzsplit(X); // columns X[:,k]

  // Root finding problem
  const G = [];
  // states fixed, costates free at initial time
  G.push(ca.minus(ca.vcat(ca.vertsplit(Xc[0]).slice(0, 2)), ca.MX([0, 1])));
  for (let k = 0; k < num_nodes; k++) {
    const XF = I.call({ x0: Xc[k] })["xf"];
    G.push(ca.minus(XF, Xc[k + 1]));
  }
  // costates fixed, states free at final time
  G.push(ca.minus(ca.vcat(ca.vertsplit(Xc[num_nodes]).slice(2, 4)), ca.MX([0, 0])));

  // Root-finding problem (vectorize X)
  const rfp = ca.Function("rfp", [ca.vec(X)], [ca.vcat(G)]);

  // Solver: "nlpsol"/"ipopt" exactly as in the Python original
  const solver = ca.rootfinder("solver", "nlpsol", rfp, {
    nlpsol: "ipopt",
    nlpsol_options: {
      "ipopt.hessian_approximation": "limited-memory",
      "ipopt.print_level": 0,
      print_time: false,
    },
  });

  // Solve
  const X_sol = solver(ca.DM.zeros(nX * (num_nodes + 1)));
  log("node-state solution (4 x " + (num_nodes + 1) + "):");
  log("X_sol = " + X_sol);

  // Visualize per-node (each 0.5s interval is stable under explicit RK; a
  // single 10s re-integration would blow up on the unstable costate ODE).
  const per = 5; // sub-samples per shooting interval
  const tgrid_node = [];
  for (let i = 0; i <= per; i++) tgrid_node.push((tf / num_nodes) * i / per);
  const vis = ca.integrator("vis", "rk", dae, 0, tgrid_node, { number_of_finite_elements: 20 });
  const segs = [];
  for (let k = 0; k < num_nodes; k++) {
    const xk = X_sol[`${k * nX}:${(k + 1) * nX}`];        // node k costate/state
    segs.push(vis.call({ x0: xk })["xf"][":,:" + per]);   // first `per` sub-samples
  }
  const sol = ca.hcat(segs);   // 4 x (num_nodes*per)
  const u_traj = u_fcn(sol);   // columnwise evaluation, as in the Python original

  log("-----");
  log("x trajectory = " + sol["0,:"]);
  log("y trajectory = " + sol["1,:"]);
  log("u trajectory = " + u_traj);
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

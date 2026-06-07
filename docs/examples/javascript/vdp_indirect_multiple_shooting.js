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

async function example(M, log) {
  // Declare variables
  const x0 = M.SX.sym("x0"), x1 = M.SX.sym("x1");
  const x = M.vertcat(x0, x1);

  // Control
  const u = M.SX.sym("u");

  // ODE right hand side
  const xdot = M.vertcat(
    M.plus(M.minus(M.times(M.minus(1, M.times(x1, x1)), x0), x1), u),
    x0);

  // Lagrangian
  const L = M.plus(M.plus(M.times(x0, x0), M.times(x1, x1)), M.times(u, u));

  // Costate
  const lam = M.SX.sym("lam", 2);

  // Hamiltonian
  const H = M.plus(M.dot(lam, xdot), L);

  // Costate equations
  const ldot = M.times(-1, M.gradient(H, x));

  log("Hamiltonian: " + H.str());

  // H is convex quadratic in u: H = u*u + p*u + q; extract p
  let p = M.gradient(H, u);
  p = M.substitute(p, u, M.SX(0));

  // Unconstrained minimizer u = -p/2, clipped to [-0.75, 1.0]
  let u_opt = M.times(-0.5, p);
  u_opt = M.fmin(u_opt, M.SX(1.0));
  u_opt = M.fmax(u_opt, M.SX(-0.75));
  log("optimal control: " + u_opt.str());

  // Augment f with ldot and substitute the optimal control
  let f = M.vertcat(xdot, ldot);
  f = M.substitute(f, u, u_opt);

  // Function for the optimal control given the augmented state
  const xlam = M.vertcat(x, lam);
  const u_fcn = new M.Function("ufcn", [xlam], [u_opt]);

  // DAE for the augmented dynamics
  const dae = { x: xlam, ode: f };

  const nX = 4;              // augmented state dimension
  const tf = 10.0;          // end time
  const num_nodes = 20;     // shooting nodes

  // Integrator over one interval -- "rk" stands in for "cvodes"
  const I = M.integrator("I", "rk", dae, 0, tf / num_nodes, { number_of_finite_elements: 50 });

  // States at each shooting node (4 x (num_nodes+1))
  const X = M.MX.sym("X", nX, num_nodes + 1);
  const Xc = M.horzsplit(X); // columns X[:,k]

  // Root finding problem
  const G = [];
  // states fixed, costates free at initial time
  G.push(M.minus(M.vertcat(...M.vertsplit(Xc[0]).slice(0, 2)), M.MX([0, 1])));
  for (let k = 0; k < num_nodes; k++) {
    const XF = I.call({ x0: Xc[k] })["xf"];
    G.push(M.minus(XF, Xc[k + 1]));
  }
  // costates fixed, states free at final time
  G.push(M.minus(M.vertcat(...M.vertsplit(Xc[num_nodes]).slice(2, 4)), M.MX([0, 0])));

  // Root-finding problem (vectorize X)
  const rfp = new M.Function("rfp", [M.vec(X)], [M.vcat(G)]);

  // Solver: "nlpsol"/"ipopt" exactly as in the Python original
  const solver = M.rootfinder("solver", "nlpsol", rfp, {
    nlpsol: "ipopt",
    nlpsol_options: {
      "ipopt.hessian_approximation": "limited-memory",
      "ipopt.print_level": 0,
      print_time: false,
    },
  });

  // Solve
  const X_sol = solver.call([M.DM.zeros(nX * (num_nodes + 1))])[0];
  log("node-state solution (4 x " + (num_nodes + 1) + "):");
  log("X_sol = " + X_sol.nonzeros().map((v) => v.toFixed(4)).join(" "));

  // Visualize per-node (each 0.5s interval is stable under explicit RK; a
  // single 10s re-integration would blow up on the unstable costate ODE).
  const per = 5; // sub-samples per shooting interval
  const tgrid_node = [];
  for (let i = 0; i <= per; i++) tgrid_node.push((tf / num_nodes) * i / per);
  const vis = M.integrator("vis", "rk", dae, 0, tgrid_node, { number_of_finite_elements: 20 });
  const nz = X_sol.nonzeros();
  const x_traj = [], y_traj = [], u_list = [], tgrid = [];
  for (let k = 0; k < num_nodes; k++) {
    const xk = M.DM(nz.slice(k * nX, k * nX + nX));
    const seg = vis.call({ x0: xk })["xf"].nonzeros(); // 4 x per
    for (let c = 0; c < per; c++) {
      const xa = seg.slice(c * nX, c * nX + nX);
      x_traj.push(xa[0]); y_traj.push(xa[1]);
      u_list.push(u_fcn.call([M.DM(xa)])[0].nonzeros()[0]);
      tgrid.push(k * (tf / num_nodes) + tgrid_node[c]);
    }
  }

  log("-----");
  log("x trajectory = " + x_traj.map((v) => v.toFixed(4)).join(" "));
  log("y trajectory = " + y_traj.map((v) => v.toFixed(4)).join(" "));
  log("u trajectory = " + u_list.map((v) => v.toFixed(4)).join(" "));
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

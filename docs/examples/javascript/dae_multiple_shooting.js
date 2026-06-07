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
// JS port of docs/examples/python/dae_multiple_shooting.py.
//
// Direct multiple shooting on a DAE optimal control problem:
//   min  int_0^10 (x0^2 + x1^2 + u^2) dt
//   s.t. dot(x0) = z*x0 - x1 + u,  dot(x1) = x0,  0 = x1^2 + z - 1
//        x0(0)=0, x1(0)=1, x0(10)=0, x1(10)=0,  -0.75 <= u <= 1
// The Python version integrates with "idas"; that DAE integrator is not
// linked in the wasm build, so this port SUBSTITUTES the "collocation"
// integrator (which also supports algebraic states).
//
// JS notes (see README.md):
//   * Python's `inf` becomes JS `Infinity`.
//   * The integrator is called per shooting interval; xf/qf are read from
//     the returned name-keyed dict.

async function example(M, log) {
  // Declare variables
  const x0 = M.SX.sym("x0");
  const x1 = M.SX.sym("x1");
  const x = M.vertcat(x0, x1); // differential states
  const z = M.SX.sym("z");     // algebraic variable
  const u = M.SX.sym("u");     // control

  // Differential equation
  const f_x = M.vertcat(M.plus(M.minus(M.times(z, x0), x1), u), x0);
  // Algebraic equation
  const f_z = M.minus(M.plus(M.power(x1, 2), z), 1);
  // Lagrange cost term (quadrature)
  const f_q = M.plus(M.plus(M.power(x0, 2), M.power(x1, 2)), M.power(u, 2));

  // Create an integrator (collocation substituted for idas), interval 0.5s
  const dae = { x: x, z: z, p: u, ode: f_x, alg: f_z, quad: f_q };
  const I = M.integrator("I", "collocation", dae, 0, 0.5);

  // Number of intervals
  const nk = 20;

  // Build the NLP via direct multiple shooting
  const w = [];   // variables
  const lbw = [];
  const ubw = [];
  const G = [];   // constraints
  let J = M.MX(0);

  // Initial conditions
  let Xk = M.MX.sym("X0", 2);
  w.push(Xk);
  lbw.push(0, 1); ubw.push(0, 1);

  for (let k = 0; k < nk; ++k) {
    // Local control
    const Uk = M.MX.sym("U" + k);
    w.push(Uk);
    lbw.push(-0.75); ubw.push(1.0);

    // Integrate over the interval
    const Ik = I.call({ x0: Xk, p: Uk });
    const X_prev = Ik["xf"];
    J = M.plus(J, Ik["qf"]);

    // Lift the state variable
    Xk = M.MX.sym("X" + (k + 1), 2);
    w.push(Xk);
    lbw.push(-Infinity, -Infinity); ubw.push(Infinity, Infinity);
    G.push(M.minus(X_prev, Xk));
  }

  // Allocate an NLP solver
  const nlp = { x: M.vertcat.apply(null, w), f: J, g: M.vertcat.apply(null, G) };
  const solver = M.nlpsol("solver", "ipopt", nlp,
    { "ipopt.print_level": 0, print_time: false });

  const sol = solver.call({
    lbx: M.DM(lbw), ubx: M.DM(ubw), lbg: M.DM(0), ubg: M.DM(0), x0: M.DM(0),
  });

  log("optimal cost = " + sol["f"].nonzeros()[0].toFixed(6));

  // Layout per interval: [x0, x1, u]; final block has [x0, x1].
  const nz = sol["x"].nonzeros();
  log("interval samples (t, x0, x1, u):");
  for (let k = 0; k <= nk; k += 4) {
    const t = (10.0 * k) / nk;
    const xa = nz[3 * k], xb = nz[3 * k + 1];
    const uu = k < nk ? nz[3 * k + 2].toFixed(4) : "  -  ";
    log("  t=" + t.toFixed(2) + "  x0=" + xa.toFixed(4) + "  x1=" + xb.toFixed(4) + "  u=" + uu);
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

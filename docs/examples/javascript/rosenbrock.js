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
// JS port of docs/examples/python/rosenbrock.py.
//
// Solve the Rosenbrock problem, formulated as the NLP:
//
//   minimize     x^2 + 100*z^2
//   subject to   z + (1-x)^2 - y == 0
//
// CasADi-in-JS notes (see README.md for the full list):
//   * No operator overloading.  Use free-function forms:
//       x*x   -> ca.times(x, x)        1 - x -> ca.minus(ca.SX(1), x)
//       a + b -> ca.plus(a, b)         x**2  -> ca.power(x, ca.SX(2))
//     Number literals in binary ops must be wrapped explicitly
//     (ca.SX/ca.MX/ca.DM) -- there is no Python-style __radd__ coercion.
//   * Concat: ca.vertcat(x, y, z) is variadic; ca.vcat([x, y, z]) takes
//     a list.  (Bare ca.vertcat([x, y, z]) is a 1-arg list call.)
//   * Solver IO are name-keyed objects (x0/lbg/ubg... -> x/f/g/...),
//     same as Python dicts.
//   * Print casadi values via string concat: "x = " + dm (toString).
//     Pull a scalar into JS with Number(dm); element/slice access uses
//     Python-style indexing: x[0], x["1:3"], x["0::2"], x["0,:"].

// The example body is environment-agnostic: it receives the loaded
// casadi module `ca` and a `log(...)` sink, so the SAME file runs both
// under Node (see the self-run block below) and in the browser (where
// rosenbrock.html calls example()).  Keep this function pure casadi.
async function example(ca, log) {
  // Declare variables
  const x = ca.SX.sym("x");
  const y = ca.SX.sym("y");
  const z = ca.SX.sym("z");

  // Formulate the NLP
  const one = ca.SX(1), c100 = ca.SX(100);
  const f = ca.plus(ca.times(x, x), ca.times(c100, ca.times(z, z)));
  const oneMinusX = ca.minus(one, x);
  const g = ca.minus(ca.plus(z, ca.times(oneMinusX, oneMinusX)), y);
  const nlp = { x: ca.vertcat(x, y, z), f: f, g: g };

  // Create an NLP solver
  const solver = ca.nlpsol("solver", "ipopt", nlp);

  // Solve the Rosenbrock problem
  const res = solver.call({
    x0:  ca.DM([2.5, 3.0, 0.75]),
    ubg: ca.DM(0),
    lbg: ca.DM(0),
  });

  // Print solution
  const fmt = (label, val) => log(label.padStart(50) + " " + val);
  log("");
  fmt("Optimal cost:",                     res["f"]);
  fmt("Primal solution:",                  res["x"]);
  fmt("Dual solution (simple bounds):",    res["lam_x"]);
  fmt("Dual solution (nonlinear bounds):", res["lam_g"]);
}

// ---- Node entry point: `node rosenbrock.js` --------------------------
// Skipped in the browser, where `require` is undefined and rosenbrock.html
// drives example() instead.
if (typeof require !== "undefined" && typeof module !== "undefined" && require.main === module) {
  const path = require("path");
  const casadiPath = process.env.CASADI_JS
    || path.resolve(__dirname, "../../../build-wasm/swig/wasm-js/casadi.js");
  require(casadiPath)()
    .then((ca) => example(ca, (...a) => console.log(...a)))
    .catch((e) => { console.error("FATAL:", e.message || e); process.exit(1); });
}

if (typeof module !== "undefined" && module.exports) module.exports = example;

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
//       x*x   -> M.times(x, x)        1 - x -> M.minus(M.SX(1), x)
//       a + b -> M.plus(a, b)         x**2  -> M.power(x, M.SX(2))
//     Number literals in binary ops must be wrapped explicitly
//     (M.SX/M.MX/M.DM) -- there is no Python-style __radd__ coercion.
//   * Concat: M.vertcat(x, y, z) is variadic; M.vcat([x, y, z]) takes
//     a list.  (Bare M.vertcat([x, y, z]) is a 1-arg list call.)
//   * Solver IO are name-keyed objects (x0/lbg/ubg... -> x/f/g/...),
//     same as Python dicts.
//   * M.DM values read out via .nonzeros() -> number[].

// The example body is environment-agnostic: it receives the loaded
// casadi module `M` and a `log(...)` sink, so the SAME file runs both
// under Node (see the self-run block below) and in the browser (where
// rosenbrock.html calls example()).  Keep this function pure casadi.
async function example(M, log) {
  // Declare variables
  const x = M.SX.sym("x");
  const y = M.SX.sym("y");
  const z = M.SX.sym("z");

  // Formulate the NLP
  const one = M.SX(1), c100 = M.SX(100);
  const f = M.plus(M.times(x, x), M.times(c100, M.times(z, z)));
  const oneMinusX = M.minus(one, x);
  const g = M.minus(M.plus(z, M.times(oneMinusX, oneMinusX)), y);
  const nlp = { x: M.vertcat(x, y, z), f: f, g: g };

  // Create an NLP solver
  const solver = M.nlpsol("solver", "ipopt", nlp);

  // Solve the Rosenbrock problem
  const res = solver.call({
    x0:  M.DM([2.5, 3.0, 0.75]),
    ubg: M.DM(0),
    lbg: M.DM(0),
  });

  // Print solution
  const fmt = (label, val) => log(label.padStart(50) + " " + val.nonzeros().join(" "));
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
    .then((M) => example(M, (...a) => console.log(...a)))
    .catch((e) => { console.error("FATAL:", e.message || e); process.exit(1); });
}

if (typeof module !== "undefined" && module.exports) module.exports = example;

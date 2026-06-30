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
// JS port of docs/examples/python/simple_nlp.py.
//
//   minimize    x0^2 + x1^2
//   subject to  x0 + x1 - 10 >= 0
//
// JS notes (see README.md):
//   * Element access: Python's x[0]/x[1] becomes destructuring of
//     ca.vertsplit(x) -> [x0, x1] (a list of scalar sub-expressions).

async function example(ca, log) {
  // Declare variables
  const x = ca.SX.sym("x", 2);
  const [x0, x1] = ca.vertsplit(x);

  // Form the NLP
  const f = ca.plus(ca.times(x0, x0), ca.times(x1, x1));      // objective
  const g = ca.minus(ca.plus(x0, x1), ca.SX(10));             // constraint
  const nlp = { x: x, f: f, g: g };

  // Allocate an ipopt solver and solve with g >= 0
  const solver = ca.nlpsol("solver", "ipopt", nlp);
  const sol = solver.call({ lbg: ca.DM(0) });

  // Print solution
  log("-----");
  log("objective at solution = " + sol["f"]);
  log("primal solution = " + sol["x"]);
  log("dual solution (x) = " + sol["lam_x"]);
  log("dual solution (g) = " + sol["lam_g"]);
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

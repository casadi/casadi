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
// JS port of docs/examples/python/simple_lp.py.
//
//   minimize    3x + 4y
//   subject to  x + 2y <= 14
//               3x -  y >= 0
//                x -  y <= 2
//
// JS notes (see README.md):
//   * Infinity works directly inside ca.DM (renders as inf).
//   * A linear program is a Conic (QP) problem with no quadratic term:
//     pass only the sparsity of the constraint matrix `a`.

async function example(ca, log) {
  // Sparsity of the LP linear term (3 constraints x 2 variables)
  const A = ca.Sparsity.dense(3, 2);

  // Create solver (qpoases handles LPs as a degenerate QP)
  const solver = ca.conic("solver", "qpoases", { a: A });

  const sol = solver.call({
    g:   ca.DM([3, 4]),
    a:   ca.DM([[1, 2], [3, -1], [1, -1]]),
    lba: ca.DM([-Infinity, 0, -Infinity]),
    uba: ca.DM([14, Infinity, 2]),
  });

  for (const k of Object.keys(sol)) {
    log(k + " = " + sol[k]);
  }
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

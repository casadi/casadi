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
// JS port of docs/examples/python/c_code_generation.py.
//
// Generates C source for the gradient of a 7x7 determinant.  The Python
// version then compiles that C with gcc and loads it back via external();
// that part is dropped here -- the wasm runtime has no C compiler, so we
// only demonstrate CodeGenerator.dump() (the generated source text) plus a
// direct numeric evaluation of the function.
//
// JS notes (see README.md):
//   * CodeGenerator(name).add(f).dump() returns the C source as a string.
//   * n_nodes() returns bigint -> Number() for arithmetic/printing.

async function example(ca, log) {
  // Expression for the gradient of the determinant of a 7x7 matrix
  const x = ca.SX.sym("x", 7, 7);
  const gd = ca.gradient(ca.det(x), x);

  // Form a function
  const name = "grad_det";
  const grad_det = ca.Function(name, [x], [gd], ["x"], ["gd"]);

  // Generate C code (returned as text -- not compiled/loaded here)
  const cg = ca.CodeGenerator(name + ".c");
  cg.add(grad_det);
  const code = cg.dump();
  log("generated C source for '" + name + "': " + code.length + " characters");
  log("--- first lines ---");
  log(code.split("\n").slice(0, 10).join("\n"));
  log("...");

  // Number of elementary operations
  const num_op = Number(grad_det.n_nodes());
  log("number of elementary operations: " + num_op);

  // Evaluate numerically at a random point
  const x0 = ca.DM.rand(7, 7);
  const r = grad_det(x0);
  log("result (first 7 entries): " + r[":7"]);
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

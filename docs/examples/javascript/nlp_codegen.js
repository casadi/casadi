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
// JS port of docs/examples/python/nlp_codegen.py.
//
// Solve a tiny NLP and generate self-contained C code for its functions.
//
// JS notes (see README.md):
//   * The Python example JIT-compiles the NLP and reloads it as an
//     external() shared library.  The wasm build ships no C compiler, so
//     we cannot compile+load.  Instead we show the portable half: emit the
//     C source via CodeGenerator and log a preview of it.

async function example(ca, log) {
  // Test problem:  min x^2 + y^2  s.t.  x + y - 10 = 0
  const x = ca.MX.sym("x");
  const y = ca.MX.sym("y");
  const f = ca.plus(ca.times(x, x), ca.times(y, y));
  const g = ca.minus(ca.plus(x, y), 10);
  const nlp = { x: ca.vertcat(x, y), f: f, g: g };

  // Create an NLP solver instance
  const solver = ca.nlpsol("solver", "ipopt", nlp);

  // Solve the NLP
  const res = solver.call({ lbx: ca.DM(-Infinity), ubx: ca.DM(Infinity), lbg: ca.DM(0), ubg: ca.DM(0), x0: ca.DM(0) });

  log("-----");
  log("objective at solution = " + res["f"]);
  log("primal solution = " + res["x"]);
  log("dual solution (x) = " + res["lam_x"]);
  log("dual solution (g) = " + res["lam_g"]);

  // Generate C code for the NLP functions (cannot compile/load here).
  const cg = ca.CodeGenerator("nlp.c");
  cg.add(solver);
  const code = cg.dump();
  const lines = code.split("\n");
  log("-----");
  log("generated nlp.c (" + lines.length + " lines); first lines:");
  for (const line of lines.slice(0, 12)) log("  " + line);
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

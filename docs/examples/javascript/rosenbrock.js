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
// JS notes:
//   * casadi loads via an async factory: `const M = await
//     require('.../casadi.js')();`.
//   * No operator overloading in JS.  Free-fn forms:
//       x**2       -> M.times(x, x)   (or M.power(x, M.SX(2)))
//       1 - x      -> M.minus(M.SX(1), x)
//       a + b      -> M.plus(a, b)
//     Number args to binary ops must be EXPLICITLY wrapped in M.SX
//     / M.MX / M.DM at construction time -- the JS-side overload
//     dispatcher doesn't auto-coerce `2` to `M.DM(2)` the way
//     Python's __rmul__/__radd__ do.  (Tracked as a separate runtime
//     ergonomics task.)
//   * Concat: `vertcat(x, y, z)` variadic, `vcat([x, y, z])` list
//     form.  `vertcat([x, y, z])` would be a 1-arg call with a list
//     of matrices -- which Python `vertcat(*[...])` and JS `vcat`
//     handle but the bare `vertcat([...])` does NOT.
//   * Solver inputs/outputs are dicts keyed by name (x0/lbg/ubg/...
//     and x/f/g/lam_x/lam_g) -- same shape as Python.
//   * `M.DM` values render via `.nonzeros()` -> number[].

const path = require("path");
const casadiPath = process.env.CASADI_JS
  || path.resolve(__dirname, "../../../build-wasm/swig/wasm-js/casadi.js");

(async () => {
  const M = await require(casadiPath)();

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
  const fmt = (label, val) => console.log(
    label.padStart(50) + " " + val.nonzeros().join(" "));
  console.log();
  fmt("Optimal cost:",                  res["f"]);
  fmt("Primal solution:",               res["x"]);
  fmt("Dual solution (simple bounds):", res["lam_x"]);
  fmt("Dual solution (nonlinear bounds):", res["lam_g"]);
})().catch((e) => { console.error("FATAL:", e.message || e); process.exit(1); });

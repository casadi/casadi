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
// JS port of docs/examples/python/parallel_map.py.
//
// Evaluate a function over many inputs, comparing a hand-written vertcat
// of N calls against Function.map(N).
//
// JS notes (see README.md):
//   * The Python example also tries f.map(N, "openmp"); the wasm build is
//     single-threaded, so we use only the default (serial) map and note it.
//   * The expensive function uses 2000 nested sin() here (vs 100000) to keep
//     the in-browser/node run snappy.

async function example(M, log) {
  // Number of inputs to evaluate
  const N = 300;

  // Dummy input
  const dummyInput = [];
  for (let k = 0; k < N; k++) dummyInput.push((2.0 * Math.PI * k) / (N - 1));

  // A moderately expensive function: nested sines
  log("creating dummy function...");
  const x = M.SX.sym("x");
  let y = x;
  for (let k = 0; k < 2000; k++) y = M.sin(y);
  const f0 = new M.Function("f", [x], [y]);

  // Evaluate serially, the old-fashioned way (vertcat of N scalar calls)
  const X = M.MX.sym("x", N);
  const [...Xs] = M.vertsplit(X);
  const Y = M.vertcat(...Xs.map((xk) => f0.call([xk])[0]));
  const fNaiveParallel = new M.Function("fParallel", [X], [Y]);

  log("evaluating naive parallel function...");
  let t0 = Date.now();
  const outNaive = fNaiveParallel.call([M.DM(dummyInput)])[0];
  log("evaluated naive parallel function in " + ((Date.now() - t0) / 1000).toFixed(3) + " seconds");

  // Evaluate it using the serial map construct
  const fMap = f0.map(N);
  log("evaluating serial map function...");
  t0 = Date.now();
  const outMap = fMap.call([M.DM([dummyInput])])[0];
  log("evaluated serial map function in " + ((Date.now() - t0) / 1000).toFixed(3) + " seconds");

  // The two have differently-shaped outputs (column vs row); compare values.
  const a = outNaive.nonzeros();
  const b = outMap.nonzeros();
  let maxdiff = 0;
  for (let i = 0; i < a.length; i++) maxdiff = Math.max(maxdiff, Math.abs(a[i] - b[i]));
  log("-----");
  log("max |naive - map| = " + maxdiff);
  log("first 5 outputs = " + b.slice(0, 5).map((v) => v.toFixed(6)).join(" "));
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

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
// JS port of docs/examples/python/accessing_mx_algorithm.py.
//
// Demonstrates how the algorithm of an MX Function can be accessed and
// its atomic operations traversed.  Unlike the SX variant, the work
// vector holds matrices; here we keep them as DM and use casadi's
// free-function arithmetic instead of numpy.
//
// JS notes (see README.md):
//   * Opcode enums (M.Operation.OP_*) and instruction_id() are both
//     bigint -- compare directly with `===`.
//   * MX has no to_DM(); evaluate a constant node with M.evalf(mx) -> DM.
//   * Element-wise * is M.times; matrix product is M.mtimes.

async function example(M, log) {
  const OP = M.Operation;
  const dm = (v) => "[" + v.nonzeros().join(", ") + "]";

  // Create a function:  r = 3*(c @ b)*a + b
  const a = M.MX.sym("a");
  const b = M.MX.sym("b", 2);
  const c = M.MX.sym("c", 2, 2);
  const r = M.plus(M.times(M.times(M.MX(3), M.mtimes(c, b)), a), b);
  const f = new M.Function("f", [a, b, c], [r], ["a", "b", "c"], ["r"]);

  // Input values of matching dimensions
  const input_val = [
    M.DM([2.0]),
    M.DM([3.0, 4.0]),
    M.DM([[5.0, 1.0], [8.0, 4.0]]),
  ];
  const output_val = [M.DM.zeros(2)];

  // Work vector
  const work = new Array(Number(f.sz_w())).fill(null);

  const n = Number(f.n_instructions());
  for (let k = 0; k < n; ++k) {
    const op = f.instruction_id(k);
    const o = f.instruction_output(k);
    const i = f.instruction_input(k);

    if (op === OP.OP_CONST) {
      const v = M.evalf(f.instruction_MX(k));
      work[o[0]] = v;
      log(`work[${o[0]}] = ${dm(v)}`);
    } else if (op === OP.OP_INPUT) {
      work[o[0]] = input_val[i[0]];
      log(`work[${o[0]}] = input[${i[0]}]            ---> ${dm(work[o[0]])}`);
    } else if (op === OP.OP_OUTPUT) {
      output_val[o[0]] = work[i[0]];
      log(`output[${o[0]}] = work[${i[0]}]             ---> ${dm(output_val[o[0]])}`);
    } else if (op === OP.OP_ADD) {
      work[o[0]] = M.plus(work[i[0]], work[i[1]]);
      log(`work[${o[0]}] = work[${i[0]}] + work[${i[1]}]      ---> ${dm(work[o[0]])}`);
    } else if (op === OP.OP_MUL) {
      work[o[0]] = M.times(work[i[0]], work[i[1]]);
      log(`work[${o[0]}] = work[${i[0]}] * work[${i[1]}]        ---> ${dm(work[o[0]])}`);
    } else if (op === OP.OP_MTIMES) {
      work[o[0]] = M.plus(M.mtimes(work[i[1]], work[i[2]]), work[i[0]]);
      log(`work[${o[0]}] = work[${i[1]}] @ work[${i[2]}] + work[${i[0]}]        ---> ${dm(work[o[0]])}`);
    } else {
      throw new Error("Unknown operation id " + op + " at instruction " + k);
    }
  }

  log("------");
  log("Evaluated " + f.str());
  log("Expected: " + dm(f.call(input_val)[0]));
  log("Got:      " + dm(output_val[0]));
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

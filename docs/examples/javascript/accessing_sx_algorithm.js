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
// JS port of docs/examples/python/accessing_sx_algorithm.py.
//
// Demonstrates how the algorithm of an SX Function can be accessed and
// its atomic operations traversed, re-evaluating it by hand.
//
// JS notes (see README.md):
//   * Opcode enums (ca.Operation.OP_*) and instruction_id() are both
//     bigint, so compare them directly with `===`.
//   * instruction_output(k) / instruction_input(k) return number[].
//   * sz_w() / n_instructions() are bigint -> wrap with Number() for loops.

async function example(ca, log) {
  const OP = ca.Operation;

  // Create a function
  const a = ca.SX.sym("a");
  const b = ca.SX.sym("b", 2);
  const f = ca.Function("f", [a, b], [ca.plus(ca.times(ca.SX(2), a), b)], ["a", "b"], ["r"]);

  // Input / output values of matching dimensions
  const input_val = [[2.0], [3.0, 4.0]];
  const output_val = [[0.0, 0.0]];

  // Work vector
  const work = new Array(Number(f.sz_w())).fill(0);

  // Loop over the algorithm
  const n = Number(f.n_instructions());
  for (let k = 0; k < n; ++k) {
    const op = f.instruction_id(k);
    const o = f.instruction_output(k);
    const i = f.instruction_input(k);

    if (op === OP.OP_CONST) {
      work[o[0]] = f.instruction_constant(k);
      log(`work[${o[0]}] = ${work[o[0]]}`);
    } else if (op === OP.OP_INPUT) {
      work[o[0]] = input_val[i[0]][i[1]];
      log(`work[${o[0]}] = input[${i[0]}][${i[1]}]            ---> ${work[o[0]]}`);
    } else if (op === OP.OP_OUTPUT) {
      output_val[o[0]][o[1]] = work[i[0]];
      log(`output[${o[0]}][${o[1]}] = work[${i[0]}]             ---> ${output_val[o[0]][o[1]]}`);
    } else if (op === OP.OP_TWICE) {
      work[o[0]] = 2 * work[i[0]];
      log(`work[${o[0]}] = 2*work[${i[0]}]                         ---> ${work[o[0]]}`);
    } else if (op === OP.OP_ADD) {
      work[o[0]] = work[i[0]] + work[i[1]];
      log(`work[${o[0]}] = work[${i[0]}] + work[${i[1]}]        ---> ${work[o[0]]}`);
    } else if (op === OP.OP_MUL) {
      work[o[0]] = work[i[0]] * work[i[1]];
      log(`work[${o[0]}] = work[${i[0]}] * work[${i[1]}]        ---> ${work[o[0]]}`);
    } else {
      throw new Error("Unknown operation id " + op + " at instruction " + k);
    }
  }

  log("------");
  log("Evaluated " + f.str());
  log("Expected: " + f(ca.DM(input_val[0]), ca.DM(input_val[1])));
  log("Got:      " + output_val[0]);
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

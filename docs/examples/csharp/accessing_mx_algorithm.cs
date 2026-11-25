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
using System;
using System.Collections.Generic;
using System.Linq;
using CasADi;

/*
Demonstration on how the algorithm of an MX function can be accessed and its operations can be transversed.
*/

public static class accessing_mx_algorithm {
  public static void Main() {
    // Create a function
    MX a = MX.sym("a");
    MX b = MX.sym("b", 2);
    MX c = MX.sym("c", 2, 2);
    var f = ca.Function("f", [a, b, c], [3 * ca.mtimes(c, b) * a + b],
                         ["a", "b", "c"], ["r"]);

    // Input values of the same dimensions as the above
    var input_val = new List<DM> { 2.0,
                                   new double[] { 3.0, 4.0 },
                                   new double[,] { { 5.0, 1.0 }, { 8.0, 4.0 } } };

    // Output values to be calculated of the same dimensions as the above
    var output_val = new List<DM> { DM.zeros(2, 1) };

    // Work vector
    var work = new DM[f.sz_w()];

    // Loop over the algorithm
    for (long k = 0; k < f.n_instructions(); k++) {

      // Get the atomic operation
      long op = f.instruction_id(k);

      List<long> o = f.instruction_output(k);
      List<long> i = f.instruction_input(k);

      if (op == ca.OP_CONST) {
        DM v = f.instruction_MX(k).to_DM();
        if (!v.is_dense()) throw new Exception("assertion failed");
        work[o[0]] = v;
        Console.WriteLine($"work[{o[0]}] = {v}");
      } else {
        if (op == ca.OP_INPUT) {
          work[o[0]] = input_val[(int)i[0]];
          Console.WriteLine($"work[{o[0]}] = input[{i[0]}]            ---> {work[o[0]]}");
        } else if (op == ca.OP_OUTPUT) {
          output_val[(int)o[0]] = work[i[0]];
          Console.WriteLine($"output[{o[0]}] = work[{i[0]}]             ---> {output_val[(int)o[0]]}");
        } else if (op == ca.OP_ADD) {
          work[o[0]] = work[i[0]] + work[i[1]];
          Console.WriteLine($"work[{o[0]}] = work[{i[0]}] + work[{i[1]}]      ---> {work[o[0]]}");
        } else if (op == ca.OP_MUL) {
          work[o[0]] = work[i[0]] * work[i[1]];
          Console.WriteLine($"work[{o[0]}] = work[{i[0]}] * work[{i[1]}]        ---> {work[o[0]]}");
        } else if (op == ca.OP_MTIMES) {
          work[o[0]] = ca.mtimes(work[i[1]], work[i[2]]) + work[i[0]];
          Console.WriteLine($"work[{o[0]}] = work[{i[1]}] @ work[{i[2]}] + work[{i[0]}]        ---> {work[o[0]]}");
        } else {
          var disp_in = i.Select(e => $"work[{e}]").ToList();
          string debug_str = ca.print_operator(f.instruction_MX(k), disp_in);
          throw new Exception($"Unknown operation: {op} -- {debug_str}");
          // When evaluating with MX, you may use an identity transformation:
          //MX node = f.instruction_MX(k);
          //List<MX> res = node.eval_mx(i.Select(ik => work[ik]).ToList());
          //for (int j = 0; j < o.Count; j++) work[o[j]] = res[j];
        }
      }
    }

    Console.WriteLine("------");
    Console.WriteLine($"Evaluated {f}");
    Console.WriteLine($"Expected:  [{string.Join(", ", f.call(input_val))}]");
    Console.WriteLine($"Got:       [{string.Join(", ", output_val)}]");
  }
}

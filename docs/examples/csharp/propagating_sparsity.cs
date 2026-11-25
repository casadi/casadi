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
using System.Text;
using CasADi;

/** \brief Demonstration of sparsity propagation in CasADi
 */

public static class propagating_sparsity {
  // Print a 64-bit unsigned integer (bvec_t is a 64-bit ulong)
  static void print_binary(ulong v) {
    var s = new StringBuilder();
    for (int k = 0; k < 64; ++k) {
      if (k % 4 == 0) s.Append(' ');
      if (k % 16 == 0) s.Append(' ');
      s.Append((v & (1UL << k)) != 0 ? 1 : 0);
    }
    Console.WriteLine(s);
  }

  public static void Main() {
    // BINDING GAP: bit-vector (bvec_t) propagation f(arg, res) / f.rev(arg, res) is not wrapped in the C# binding (as in Java)
    Console.WriteLine("SKIP: not wrapped in the C# binding (as in Java): no forward/reverse bit-vector (bvec_t) dependency propagation on Function");
  }

  // The example proper, to be called from Main once the propagation calls below exist
  static void run() {
    // Test both SX and MX
    for (int test = 0; test < 2; ++test) {

      // Create a simple function
      Function f;
      if (test == 0) {
        Console.WriteLine("SX:");
        SX x = SX.sym("x", 3);
        SX z = x[0] * x[0] + x[2] + 3;
        f = ca.Function("f", [x], [z]);
      } else {
        Console.WriteLine("MX:");
        MX x = MX.sym("x", 3);
        MX z = x[0] * x[0] + x[2] + 3;
        f = ca.Function("f", [x], [z]);
      }

      // Memory for inputs and outputs
      var f_in = new ulong[f.nnz_in(0)];
      var f_out = new ulong[f.nnz_out(0)];

      // Propagate from input to output (forward mode)
      Console.WriteLine("forward mode");

      // Make sure that the class is able to support the dependency propagation
      if (!f.has_spfwd()) throw new Exception("Forward sparsity propagation not supported");

      // Pass seeds
      f_in[0] = 1UL << 0; // seed in direction 0
      f_in[1] = 1UL << 2; // seed in direction 2
      f_in[2] = (1UL << 4) | (1UL << 63); // seed in direction 4 and 63

      // Reset sensitivities
      f_out[0] = 0;

      // Propagate dependencies
      // BINDING GAP (not wrapped in the C# binding, as in Java): should read  f.sp_forward([f_in], [f_out]);

      // Print the result
      print_binary(f_out[0]);

      // Propagate from output to input (adjoint/reverse/backward mode)
      Console.WriteLine("backward mode");

      // Make sure that the class is able to support the dependency propagation
      if (!f.has_sprev()) throw new Exception("Backward sparsity propagation not supported");

      // Pass seeds
      f_out[0] = (1UL << 5) | (1UL << 6); // seed in direction 5 and 6

      // Reset sensitivities
      f_in[0] = 0;
      f_in[1] = 0;
      f_in[2] = 0;

      // Propagate dependencies
      // BINDING GAP (not wrapped in the C# binding, as in Java): should read  f.sp_reverse([f_in], [f_out]);

      // Print the result
      print_binary(f_in[0]);
      print_binary(f_in[1]);
      print_binary(f_in[2]);
    }
  }
}

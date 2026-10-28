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
import org.casadi.*;
import java.util.List;

/** \brief Demonstration of sparsity propagation in CasADi
 */

public class propagating_sparsity {
  // Print a 64-bit unsigned integer (bvec_t is a 64-bit Java long)
  static void print_binary(long v) {
    StringBuilder s = new StringBuilder();
    for (int k = 0; k < Long.SIZE; ++k) {
      if (k % 4 == 0) s.append(" ");
      if (k % 16 == 0) s.append(" ");
      s.append((v & (1L << k)) != 0 ? 1 : 0);
    }
    System.out.println(s);
  }

  public static void main(String[] args) {
    // No Java counterpart of the bvec_t overloads f(arg, res) and f.rev(arg, res)
    System.out.println("SKIP: BINDING: no forward/reverse bit-vector (bvec_t) dependency propagation on Function");
  }

  // The example proper, to be called from main once the propagation calls below exist
  static void run() {
    // Test both SX and MX
    for (int test = 0; test < 2; ++test) {

      // Create a simple function
      Function f;
      if (test == 0) {
        System.out.println("SX:");
        SX x = SX.sym("x", 3);
        SX z = x.get(0).times(x.get(0)).plus(x.get(2)).plus(3);
        f = new Function("f", List.of(x), List.of(z));
      } else {
        System.out.println("MX:");
        MX x = MX.sym("x", 3);
        MX z = x.get(0).times(x.get(0)).plus(x.get(2)).plus(3);
        f = new Function("f", List.of(x), List.of(z));
      }

      // Memory for inputs and outputs
      long[] f_in = new long[(int) f.nnz_in(0)];
      long[] f_out = new long[(int) f.nnz_out(0)];

      // Propagate from input to output (forward mode)
      System.out.println("forward mode");

      // Make sure that the class is able to support the dependency propagation
      if (!f.has_spfwd()) throw new CasadiException("Forward sparsity propagation not supported");

      // Pass seeds
      f_in[0] = 1L << 0; // seed in direction 0
      f_in[1] = 1L << 2; // seed in direction 2
      f_in[2] = (1L << 4) | (1L << 63); // seed in direction 4 and 63

      // Reset sensitivities
      f_out[0] = 0;

      // Propagate dependencies
      // BINDING: should read  f.sp_forward(List.of(f_in), List.of(f_out));

      // Print the result
      print_binary(f_out[0]);

      // Propagate from output to input (adjoint/reverse/backward mode)
      System.out.println("backward mode");

      // Make sure that the class is able to support the dependency propagation
      if (!f.has_sprev()) throw new CasadiException("Backward sparsity propagation not supported");

      // Pass seeds
      f_out[0] = (1L << 5) | (1L << 6); // seed in direction 5 and 6

      // Reset sensitivities
      f_in[0] = 0;
      f_in[1] = 0;
      f_in[2] = 0;

      // Propagate dependencies
      // BINDING: should read  f.sp_reverse(List.of(f_in), List.of(f_out));

      // Print the result
      print_binary(f_in[0]);
      print_binary(f_in[1]);
      print_binary(f_in[2]);
    }
  }
}

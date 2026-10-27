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
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/*
Demonstration on how the algorithm of an SX function can be accessed and its operations can be transversed.
*/

public class accessing_sx_algorithm {
  // Python-style print: arguments separated by a space
  static void print(Object... a) {
    StringBuilder s = new StringBuilder();
    for (int k = 0; k < a.length; k++) s.append(k > 0 ? " " : "").append(a[k]);
    System.out.println(s);
  }

  public static void main(String[] args) {
    // Create a function
    SX a = SX.sym("a");
    SX b = SX.sym("b", 2);
    Function f = new Function("f", List.of(a, b), List.of(ca.times(2, a).plus(b)), List.of("a", "b"), List.of("r"));

    // Input values of the same dimensions as the above
    double[][] input_val = {{2.0},
                            {3.0, 4.0}};

    // Output values to be calculated of the same dimensions as the above
    double[][] output_val = {new double[2]};

    // Work vector
    double[] work = new double[(int) f.sz_w()];

    // For debugging
    SX instr = f.instructions_sx();

    // Loop over the algorithm
    for (long k = 0; k < f.n_instructions(); k++) {

      // Get the atomic operation
      Operation op = Operation.swigToEnum((int) f.instruction_id(k));
      int[] o = Arrays.stream(f.instruction_output(k)).mapToInt(e -> (int) e).toArray();
      int[] i = Arrays.stream(f.instruction_input(k)).mapToInt(e -> (int) e).toArray();

      if (op == Operation.OP_CONST) {
        work[o[0]] = f.instruction_constant(k);
        print("work[", o[0], "] = ", f.instruction_constant(k));
      } else {
        if (op == Operation.OP_INPUT) {
          work[o[0]] = input_val[i[0]][i[1]];
          print("work[", o[0], "] = input[", i[0], "][", i[1], "]", "            ---> ", work[o[0]]);
        } else if (op == Operation.OP_OUTPUT) {
          output_val[o[0]][o[1]] = work[i[0]];
          print("output[", o[0], "][", o[1], "] = work[", i[0], "]", "             ---> ", output_val[o[0]][o[1]]);
        } else if (op == Operation.OP_TWICE) {
          work[o[0]] = 2 * work[i[0]];
          print("work[", o[0], "] = 2*work[", i[0], "]", "                         ---> ", work[o[0]]);
        } else if (op == Operation.OP_ADD) {
          work[o[0]] = work[i[0]] + work[i[1]];
          print("work[", o[0], "] = work[", i[0], "] + work[", i[1], "]", "        ---> ", work[o[0]]);
        } else if (op == Operation.OP_MUL) {
          work[o[0]] = work[i[0]] * work[i[1]];
          print("work[", o[0], "] = work[", i[0], "] * work[", i[1], "]", "        ---> ", work[o[0]]);
        } else {
          List<String> disp_in = new ArrayList<>();
          for (int e : i) disp_in.add("work[" + e + "]");
          String debug_str = ca.print_operator(instr.get(k), disp_in);
          throw new RuntimeException("Unknown operation: " + op + " -- " + debug_str);
        }
      }
    }

    print("------");
    print("Evaluated " + f);
    print("Expected: ", f.call(List.of(new DM(input_val[0]), new DM(input_val[1]))));
    print("Got:      ", Arrays.deepToString(output_val));
  }
}

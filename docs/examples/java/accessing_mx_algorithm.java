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
import java.util.List;

/*
Demonstration on how the algorithm of an MX function can be accessed and its operations can be transversed.
*/

public class accessing_mx_algorithm {
  public static void main(String[] args) {
    // Create a function
    MX a = MX.sym("a");
    MX b = MX.sym("b", 2);
    MX c = MX.sym("c", 2, 2);
    Function f = new Function("f", List.of(a, b, c), List.of(ca.times(3, ca.mtimes(c, b)).times(a).plus(b)),
                              List.of("a", "b", "c"), List.of("r"));

    // Input values of the same dimensions as the above
    List<DM> input_val = List.of(new DM(2.0),
                                 new DM(new double[] {3.0, 4.0}),
                                 new DM(new double[][] {{5.0, 1.0}, {8.0, 4.0}}));

    // Output values to be calculated of the same dimensions as the above
    List<DM> output_val = new ArrayList<>(List.of(DM.zeros(2, 1)));

    // Work vector
    DM[] work = new DM[(int) f.sz_w()];

    // Loop over the algorithm
    for (long k = 0; k < f.n_instructions(); k++) {

      // Get the atomic operation
      Operation op = Operation.swigToEnum((int) f.instruction_id(k));

      long[] o = f.instruction_output(k);
      long[] i = f.instruction_input(k);

      if (op == Operation.OP_CONST) {
        DM v = f.instruction_MX(k).to_DM();
        if (!v.is_dense()) throw new AssertionError();
        work[(int) o[0]] = v;
        System.out.println("work[" + o[0] + "] = " + v);
      } else {
        if (op == Operation.OP_INPUT) {
          work[(int) o[0]] = input_val.get((int) i[0]);
          System.out.println("work[" + o[0] + "] = input[" + i[0] + "]            ---> " + work[(int) o[0]]);
        } else if (op == Operation.OP_OUTPUT) {
          output_val.set((int) o[0], work[(int) i[0]]);
          System.out.println("output[" + o[0] + "] = work[" + i[0] + "]             ---> " + output_val.get((int) o[0]));
        } else if (op == Operation.OP_ADD) {
          work[(int) o[0]] = work[(int) i[0]].plus(work[(int) i[1]]);
          System.out.println("work[" + o[0] + "] = work[" + i[0] + "] + work[" + i[1] + "]      ---> " + work[(int) o[0]]);
        } else if (op == Operation.OP_MUL) {
          work[(int) o[0]] = work[(int) i[0]].times(work[(int) i[1]]);
          System.out.println("work[" + o[0] + "] = work[" + i[0] + "] * work[" + i[1] + "]        ---> " + work[(int) o[0]]);
        } else if (op == Operation.OP_MTIMES) {
          work[(int) o[0]] = ca.mtimes(work[(int) i[1]], work[(int) i[2]]).plus(work[(int) i[0]]);
          System.out.println("work[" + o[0] + "] = work[" + i[1] + "] @ work[" + i[2] + "] + work[" + i[0] + "]        ---> " + work[(int) o[0]]);
        } else {
          List<String> disp_in = new ArrayList<>();
          for (long e : i) disp_in.add("work[" + e + "]");
          String debug_str = ca.print_operator(f.instruction_MX(k), disp_in);
          throw new RuntimeException("Unknown operation: " + op + " -- " + debug_str);
          // When evaluating with MX, you may use an identity transformation:
          //MX node = f.instruction_MX(k);
          //List<MX> margs = new ArrayList<>(); for (long ik : i) margs.add(work[(int) ik]);
          //List<MX> res = node.eval_mx(margs);
          //for (int j = 0; j < o.length; j++) work[(int) o[j]] = res.get(j);
        }
      }
    }

    System.out.println("------");
    System.out.println("Evaluated " + f);
    System.out.println("Expected:  " + f.call(input_val));
    System.out.println("Got:       " + output_val);
  }
}

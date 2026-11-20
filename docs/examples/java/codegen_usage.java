//
//     MIT No Attribution
//
//     Copyright 2023 Joel Andersson, Joris Gillis, Moritz Diehl, KU Leuven.
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
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.DoubleBuffer;
import java.util.List;
import java.util.Map;

/**
 *  Example program demonstrating the usage of C-code generated from CasADi
 *  Generated code is encapsulated in self-contained code with the entry points
 *  defined below.
 *  Note that other usage, e.g. accessing the internal data structures in the
 *  generated files is not recommended and subject to change.
 *
 *  We show how the generated code can be used from C (or C++), without requiring
 *  any CasADi classes as well as how to use it from C++ using CasADi's external.
 *
 *  Joel Andersson, 2013-2015
 */

public class codegen_usage {

  // Usage from C
  // Java has no dlopen/dlsym (java.lang.foreign is final only from Java 22), so the C entry points
  // f_n_in, f_work, f_sparsity_in, f, f_checkout/f_release are reached through external and a FunctionBuffer
  static int usage_c() {
    System.out.println("---");
    System.out.println("Standalone usage from C/C++:");
    System.out.println();

    /* Load the dll */
    Function f;
    try {
      f = ca.external("f", "./f.so");
    } catch (CasadiException e) {
      System.out.println("Cannot open f.so, error " + e.getMessage());
      return 1;
    }

    /* Number of inputs */
    long n_in = f.n_in();

    /* Number of outputs */
    long n_out = f.n_out();

    /* Get sizes of the required work vectors */
    // Function::sz_arg/sz_res include the n_in/n_out call overhead, so they exceed f_work's 2/2
    long sz_arg = f.sz_arg(), sz_res = f.sz_res(), sz_iw = f.sz_iw(), sz_w = f.sz_w();
    System.out.println("Work vector sizes:");
    System.out.println("sz_arg = " + sz_arg + ", sz_res = " + sz_res + ", sz_iw = " + sz_iw + ", sz_w = " + sz_w);
    System.out.println();

    /* Print the sparsities of the inputs and outputs */
    for (long i = 0; i < n_in + n_out; ++i) {
      // Retrieve the sparsity pattern - CasADi uses column compressed storage (CCS)
      Sparsity sp_i;
      if (i < n_in) {
        System.out.println("Input " + i);
        sp_i = f.sparsity_in(i);
      } else {
        System.out.println("Output " + (i - n_in));
        sp_i = f.sparsity_out(i - n_in);
      }
      long nrow = sp_i.size1(); /* Number of rows */
      long ncol = sp_i.size2(); /* Number of columns */
      if (sp_i.is_dense()) {
        System.out.println("  Dimension: " + nrow + "-by-" + ncol + " (dense)");
      } else {
        long[] colind = sp_i.colind(); /* Column offsets */
        long[] row = sp_i.row(); /* Row nonzero */
        long nnz = sp_i.nnz(); /* Number of nonzeros */

        /* Print the pattern */
        System.out.println("  Dimension: " + nrow + "-by-" + ncol + " (" + nnz + " nonzeros)");
        StringBuilder s = new StringBuilder("  Nonzeros: {");
        for (int cc = 0; cc < ncol; ++cc) {                       /* loop over columns */
          for (int el = (int) colind[cc]; el < colind[cc + 1]; ++el) { /* loop over the nonzeros entries of the column */
            if (el != 0) s.append(", ");                          /* Separate the entries */
            long rr = row[el];                                    /* Get the row */
            s.append("{" + rr + "," + cc + "}");                  /* Print the nonzero */
          }
        }
        System.out.println(s + "}");
        System.out.println();
      }
    }

    /* Function input and output */
    DoubleBuffer x_val = doubles(4).put(new double[] {1, 2, 3, 4});
    DoubleBuffer y_val = doubles(1).put(new double[] {5});
    DoubleBuffer res0 = doubles(1);
    DoubleBuffer res1 = doubles(4);

    // Allocate work vectors and checkout thread-local memory
    try (FunctionBuffer buf = f.buffer()) {
      /* Evaluate the function */
      buf.set_arg(0, x_val);
      buf.set_arg(1, y_val);
      buf.set_res(0, res0);
      buf.set_res(1, res1);

      // Evaluation is thread-safe
      buf.eval();
      if (buf.ret() != 0) return 1;
    } // Release thread-local memory

    /* Print result of evaluation */
    System.out.println("result (0): " + g(res0.get(0)));
    System.out.println("result (1): [" + g(res1.get(0)) + "," + g(res1.get(1)) + ";" + g(res1.get(2)) + "," + g(res1.get(3)) + "]");

    return 0;
  }

  // CasADi from here on
  static void usage_cplusplus() {
    System.out.println("---");
    System.out.println("Usage from CasADi C++:");
    System.out.println();

    // Use CasADi's "external" to load the compiled function (compiled as f.so on every platform)
    Function f = ca.external("f", "./f.so");

    // Use like any other CasADi function
    double[] x = {1, 2, 3, 4};
    List<DM> arg = List.of(ca.reshape(new DM(x), 2, 2), new DM(5));
    List<DM> res = f.call(arg);

    System.out.println("result (0): " + res.get(0));
    System.out.println("result (1): " + res.get(1));
  }

  // Usage from C using mem.h
  // casadi/mem.h (casadi_alloc/casadi_eval/casadi_free on f_functions) is a C-only API; from Java the
  // library is again loaded with external and evaluated through a FunctionBuffer
  static int usage_c_with_mem() {
    System.out.println("---");
    System.out.println("Usage from C/C++ with casadi/mem.h:");
    System.out.println();

    /* Load the dll */
    Function f;
    try {
      f = ca.external("f", "./f_with_mem.so");
    } catch (CasadiException e) {
      System.out.println("Cannot open f_with_mem.so, error " + e.getMessage());
      return 1;
    }

    /* Function input and output */
    DoubleBuffer x_val = doubles(4).put(new double[] {1, 2, 3, 4});
    DoubleBuffer y_val = doubles(1).put(new double[] {5});
    DoubleBuffer res0 = doubles(1);
    DoubleBuffer res1 = doubles(4);

    try (FunctionBuffer mem = f.buffer()) {
      /* Evaluate the function */
      mem.set_arg(0, x_val);
      mem.set_arg(1, y_val);
      mem.set_res(0, res0);
      mem.set_res(1, res1);

      mem.eval();
    }

    /* Print result of evaluation */
    System.out.println("result (0): " + g(res0.get(0)));
    System.out.println("result (1): [" + g(res1.get(0)) + "," + g(res1.get(1)) + ";" + g(res1.get(2)) + "," + g(res1.get(3)) + "]");

    /* Success */
    return 0;
  }

  public static void main(String[] args) throws Exception {
    if (!hasCommand("gcc")) { System.out.println("SKIP: no C compiler (gcc) found"); return; }

    // Variables
    SX x = SX.sym("x", 2, 2);
    SX y = SX.sym("y");

    // Simple function
    Function f = new Function("f", List.of(x, y), List.of(ca.sqrt(y).minus(1), ca.sin(x).minus(y)));

    // Generate C-code
    f.generate("f");

    // Compile the C-code to a shared library
    String compile_command = "gcc -fPIC -shared -O3 f.c -o f.so";
    int flag = system(compile_command);
    if (flag != 0) throw new CasadiException("Compilation failed");

    // Usage from C
    flag = usage_c();
    if (flag != 0) throw new CasadiException("Example failed");

    // Usage from C++
    usage_cplusplus();

    // Generate C-code
    f.generate("f_with_mem", Map.of("with_mem", true, "force_canonical", true));

    // Compile the C-code to a shared library
    String env_inc = System.getenv("CASADI_INCLUDE_DIR");
    // sh drops backslashes: pass a Windows path with forward slashes
    String include_dir = (env_inc != null ? env_inc : GlobalOptions.getCasadiIncludePath()).replace('\\', '/');
    compile_command = "gcc -fPIC -I" + include_dir + " -shared -g f_with_mem.c -o f_with_mem.so";
    flag = system(compile_command);
    if (flag != 0) throw new CasadiException("Compilation failed");

    // Usage from c with mem.h
    usage_c_with_mem();
  }

  // A direct buffer of n doubles, as FunctionBuffer requires
  static DoubleBuffer doubles(int n) {
    return ByteBuffer.allocateDirect(8 * n).order(ByteOrder.nativeOrder()).asDoubleBuffer();
  }

  // Format like printf's %g
  static String g(double v) {
    String s = String.format("%g", v);
    String m = s.contains("e") ? s.substring(0, s.indexOf('e')) : s;
    String e = s.contains("e") ? s.substring(s.indexOf('e')) : "";
    if (m.contains(".")) m = m.replaceAll("0+$", "").replaceAll("\\.$", "");
    return m + e;
  }

  // Run a shell command, like system()
  static int system(String cmd) throws Exception {
    return new ProcessBuilder("sh", "-c", cmd).inheritIO().start().waitFor();
  }

  static boolean hasCommand(String cmd) {
    try {
      return new ProcessBuilder(cmd, "--version").redirectErrorStream(true)
          .redirectOutput(ProcessBuilder.Redirect.DISCARD).start().waitFor() == 0;
    } catch (Exception e) {
      return false;
    }
  }
}

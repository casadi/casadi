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
import java.util.Arrays;
import java.util.List;

public class c_code_generation {
  public static void main(String[] args) throws Exception {
    if (!hasCommand("gcc")) { System.out.println("SKIP: no C compiler (gcc) found"); return; }

    boolean compileme = true;  // Flag wether compiling should take place or not

    if (args.length > 0 && args[0].equals("nc")) {
      compileme = false;
    } else {
      System.out.println("Info: Use 'java -cp casadi.jar c_code_generation.java nc' to omit compiling");
    }

    // Form an expression for the gradient of the determinant
    SX x = SX.sym("x", 7, 7);
    SX gd = ca.gradient(ca.det(x), x);

    // Random point to evaluate it
    DM x0 = DM.rand(7, 7);

    // Form a function and generate C code
    String name = "grad_det";
    Function grad_det = new Function(name, List.of(x), List.of(gd), List.of("x"), List.of("gd"));
    String cname = grad_det.generate();

    String oname_no_opt = name + "_no_opt.so";
    System.out.println("Compiling without optimization:  " + oname_no_opt);
    double t1 = time();
    if (compileme) {
      system("gcc -fPIC -shared " + cname + " -o " + oname_no_opt);
    }
    double t2 = time();
    System.out.println("time =  " + (t2 - t1) * 1e3 + "  ms");

    String oname_O3 = name + "_O3.so";
    System.out.println("Compiling with O3 optimization:  " + oname_O3);
    t1 = time();
    if (compileme) {
      system("gcc -fPIC -shared -O3 " + cname + " -o " + oname_O3);
    }
    t2 = time();
    System.out.println("time =  " + (t2 - t1) * 1e3 + "  ms");

    String oname_Os = name + "_Os.so";
    System.out.println("Compiling with Os optimization:  " + oname_Os);
    t1 = time();
    if (compileme) {
      system("gcc -fPIC -shared -Os " + cname + " -o " + oname_Os);
    }
    t2 = time();
    System.out.println("time =  " + (t2 - t1) * 1e3 + "  ms");

    // Read function
    Function grad_det_no_opt = ca.external(name, "./" + oname_no_opt);
    Function grad_det_O3 = ca.external(name, "./" + oname_O3);
    Function grad_det_Os = ca.external(name, "./" + oname_Os);
    List<Function> f_test = List.of(grad_det, grad_det_no_opt, grad_det_O3, grad_det_Os);

    for (Function f : f_test) {
      t1 = time();
      int nrep = 10000;
      DM r = null;
      for (int rep = 0; rep < nrep; rep++) {
        r = (DM) f.invoke(x0);
      }
      t2 = time();
      System.out.println("result =  " + Arrays.toString(r.nonzeros()));
      double dt = (t2 - t1) / nrep;
      System.out.println("time =  " + dt * 1e3 + "  ms");

      long num_op = grad_det.n_nodes();
      System.out.println("number of elementary operations:  " + num_op);
      System.out.println("time per elementary operations:  " + dt / num_op * 1e9 + "  ns");
    }
  }

  // Wall time in seconds
  static double time() { return System.nanoTime() * 1e-9; }

  // Run a shell command, like os.system
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

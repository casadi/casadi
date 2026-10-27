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
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

public class nlp_codegen {
  public static void main(String[] args) throws Exception {
    if (!ca.has_nlpsol("ipopt")) { System.out.println("SKIP: nlpsol ipopt not available"); return; }

    // Test problem
    //
    //    min x^2 + y^2
    //    s.t.    x + y - 10 = 0
    //

    // Optimization variables
    MX x = MX.sym("x");
    MX y = MX.sym("y");

    // Objective
    MX f = x.times(x).plus(y.times(y));

    // Constraints
    MX g = x.plus(y).minus(10);

    // Create an NLP problem structure
    Map<String, MX> nlp = Map.of("x", ca.vertcat(List.of(x, y)), "f", f, "g", g);

    // Pick a compiler
    String compiler = "gcc";    // Linux
    // String compiler = "clang";  // OSX
    // String compiler = "cl.exe"; // Windows
    if (!hasCommand(compiler) || !Importer.has_plugin("shell")) {
      System.out.println("SKIP: no C compiler (" + compiler + ") or importer shell");
      return;
    }

    // Run this program in an environment that recognised the compiler as command.
    // On Windows, the suggested way is to run this program from a "x64 Native Tools Command Promt for VS" (Requires Visual C++ components or Build Tools for Visual Studio, available from Visual Studio installer. You also need SDK libraries in order to access stdio and math.)

    List<String> flags = List.of("-O3"); // Linux/OSX

    for (String mode : List.of("jit", "external")) {
      Function solver;

      if (mode.equals("jit")) {
        // By default, the compiler will be gcc or cl.exe
        Map<String, Object> jit_options = Map.of("flags", flags, "verbose", true, "compiler", compiler);
        Map<String, Object> options = Map.of("jit", true, "compiler", "shell", "jit_options", jit_options);

        // Create an NLP solver instance
        solver = ca.nlpsol("solver", "ipopt", nlp, options);

      } else {
        // Create an NLP solver instance
        solver = ca.nlpsol("solver", "ipopt", nlp);

        // Generate C code for the NLP functions
        solver.generate_dependencies("nlp.c");

        // On Windows, use other flags
        List<String> cmd_args = new java.util.ArrayList<>(List.of(compiler, "-fPIC", "-shared"));
        cmd_args.addAll(flags);
        cmd_args.addAll(List.of("nlp.c", "-o", "nlp.so"));
        new ProcessBuilder(cmd_args).inheritIO().start().waitFor();

        // Create a new NLP solver instance from the compiled code
        solver = ca.nlpsol("solver", "ipopt", "./nlp.so");
      }

      Map<String, Object> arg = new LinkedHashMap<>();

      arg.put("lbx", DM.inf().unaryMinus());
      arg.put("ubx", DM.inf());
      arg.put("lbg", 0);
      arg.put("ubg", 0);
      arg.put("x0", 0);

      // Solve the NLP
      Map<String, GenericMatrix> res = solver.invoke(arg);

      // Print solution
      System.out.println("-----");
      System.out.println("objective at solution = " + res.get("f"));
      System.out.println("primal solution = " + res.get("x"));
      System.out.println("dual solution (x) = " + res.get("lam_x"));
      System.out.println("dual solution (g) = " + res.get("lam_g"));
    }
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

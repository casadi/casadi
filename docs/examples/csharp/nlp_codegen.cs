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
using System.Diagnostics;
using CasADi;

public static class nlp_codegen {
  public static void Main() {
    if (!ca.has_nlpsol("ipopt")) { Console.WriteLine("SKIP: nlpsol ipopt not available"); return; }

    // Test problem
    //
    //    min x^2 + y^2
    //    s.t.    x + y - 10 = 0
    //

    // Optimization variables
    MX x = MX.sym("x");
    MX y = MX.sym("y");

    // Objective
    MX f = x * x + y * y;

    // Constraints
    MX g = x + y - 10;

    // Create an NLP problem structure
    var nlp = new Dictionary<string, MX> { ["x"] = ca.vertcat(x, y), ["f"] = f, ["g"] = g };

    // Pick a compiler
    string compiler = "gcc";    // Linux
    // string compiler = "clang";  // OSX
    // string compiler = "cl.exe"; // Windows
    if (!hasCommand(compiler) || !Importer.has_plugin("shell")) {
      Console.WriteLine($"SKIP: no C compiler ({compiler}) or importer shell");
      return;
    }

    // Run this program in an environment that recognised the compiler as command.
    // On Windows, the suggested way is to run this program from a "x64 Native Tools Command Promt for VS" (Requires Visual C++ components or Build Tools for Visual Studio, available from Visual Studio installer. You also need SDK libraries in order to access stdio and math.)

    var flags = new List<string> { "-O3" }; // Linux/OSX

    foreach (string mode in new[] { "jit", "external" }) {
      Function solver;

      if (mode == "jit") {
        // By default, the compiler will be gcc or cl.exe
        var jit_options = new Dictionary<string, object> { ["flags"] = flags, ["verbose"] = true, ["compiler"] = compiler };
        var options = new Dictionary<string, object> { ["jit"] = true, ["compiler"] = "shell", ["jit_options"] = jit_options };

        // Create an NLP solver instance
        solver = ca.nlpsol("solver", "ipopt", nlp, options);

      } else {
        // Create an NLP solver instance
        solver = ca.nlpsol("solver", "ipopt", nlp);

        // Generate C code for the NLP functions
        solver.generate_dependencies("nlp.c");

        // On Windows, use other flags
        var cmd_args = new List<string> { "-fPIC", "-shared" };
        cmd_args.AddRange(flags);
        cmd_args.AddRange(["nlp.c", "-o", "nlp.so"]);
        Process.Start(compiler, cmd_args).WaitForExit();

        // Create a new NLP solver instance from the compiled code
        solver = ca.nlpsol("solver", "ipopt", "./nlp.so");
      }

      var arg = new Dictionary<string, DM>();

      arg["lbx"] = -DM.inf();
      arg["ubx"] = DM.inf();
      arg["lbg"] = 0;
      arg["ubg"] = 0;
      arg["x0"] = 0;

      // Solve the NLP
      var res = solver.call(arg);

      // Print solution
      Console.WriteLine("-----");
      Console.WriteLine($"objective at solution = {res["f"]}");
      Console.WriteLine($"primal solution = {res["x"]}");
      Console.WriteLine($"dual solution (x) = {res["lam_x"]}");
      Console.WriteLine($"dual solution (g) = {res["lam_g"]}");
    }
  }

  static bool hasCommand(string cmd) {
    try {
      var psi = new ProcessStartInfo(cmd, "--version") { RedirectStandardOutput = true, RedirectStandardError = true };
      using var p = Process.Start(psi);
      p.StandardOutput.ReadToEnd();
      p.StandardError.ReadToEnd();
      p.WaitForExit();
      return p.ExitCode == 0;
    } catch (Exception) {
      return false;
    }
  }
}

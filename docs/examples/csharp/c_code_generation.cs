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

public static class c_code_generation {
  public static void Main(string[] args) {
    if (!hasCommand("gcc")) { Console.WriteLine("SKIP: no C compiler (gcc) found"); return; }

    bool compileme = true;  // Flag wether compiling should take place or not

    if (args.Length > 0 && args[0] == "nc") {
      compileme = false;
    } else {
      Console.WriteLine("Info: Use 'dotnet Examples.dll c_code_generation nc' to omit compiling");
    }

    // Form an expression for the gradient of the determinant
    SX x = SX.sym("x", 7, 7);
    SX gd = ca.gradient(ca.det(x), x);

    // Random point to evaluate it
    DM x0 = DM.rand(7, 7);

    // Form a function and generate C code
    string name = "grad_det";
    var grad_det = ca.Function(name, [x], [gd], ["x"], ["gd"]);
    string cname = grad_det.generate();

    string oname_no_opt = name + "_no_opt.so";
    Console.WriteLine($"Compiling without optimization:  {oname_no_opt}");
    double t1 = time();
    if (compileme) {
      system($"gcc -fPIC -shared {cname} -o {oname_no_opt}");
    }
    double t2 = time();
    Console.WriteLine($"time =  {(t2 - t1) * 1e3}  ms");

    string oname_O3 = name + "_O3.so";
    Console.WriteLine($"Compiling with O3 optimization:  {oname_O3}");
    t1 = time();
    if (compileme) {
      system($"gcc -fPIC -shared -O3 {cname} -o {oname_O3}");
    }
    t2 = time();
    Console.WriteLine($"time =  {(t2 - t1) * 1e3}  ms");

    string oname_Os = name + "_Os.so";
    Console.WriteLine($"Compiling with Os optimization:  {oname_Os}");
    t1 = time();
    if (compileme) {
      system($"gcc -fPIC -shared -Os {cname} -o {oname_Os}");
    }
    t2 = time();
    Console.WriteLine($"time =  {(t2 - t1) * 1e3}  ms");

    // Read function
    Function grad_det_no_opt = ca.external(name, "./" + oname_no_opt);
    Function grad_det_O3 = ca.external(name, "./" + oname_O3);
    Function grad_det_Os = ca.external(name, "./" + oname_Os);
    var f_test = new List<Function> { grad_det, grad_det_no_opt, grad_det_O3, grad_det_Os };

    foreach (Function f in f_test) {
      t1 = time();
      int nrep = 10000;
      DM r = null;
      for (int rep = 0; rep < nrep; rep++) {
        r = f.Call(x0);
      }
      t2 = time();
      Console.WriteLine($"result =  [{string.Join(", ", r.nonzeros())}]");
      double dt = (t2 - t1) / nrep;
      Console.WriteLine($"time =  {dt * 1e3}  ms");

      long num_op = grad_det.n_nodes();
      Console.WriteLine($"number of elementary operations:  {num_op}");
      Console.WriteLine($"time per elementary operations:  {dt / num_op * 1e9}  ns");
    }
  }

  // Wall time in seconds
  static double time() => Stopwatch.GetTimestamp() / (double)Stopwatch.Frequency;

  // Run a shell command, like os.system
  static int system(string cmd) {
    using var p = Process.Start("sh", ["-c", cmd]);
    p.WaitForExit();
    return p.ExitCode;
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

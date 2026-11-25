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
using System.IO;
using System.IO.Compression;
using CasADi;

public static class modelica_fmu_import {
  public static void Main() {
    if (!hasCommand("omc")) { Console.WriteLine("SKIP: OpenModelica compiler (omc) not found"); return; }

    // Consider the following simple Modelica file
    Console.WriteLine(File.ReadAllText("../python/rocket.mo"));

    // Different tools exist to compile Modelica models like these into FMUs.
    // The following code shows how we could do it using the open-source tool
    // OpenModelica (here its command line compiler omc, running a .mos script),
    // cf. the OpenModelica User's Guide.
    File.WriteAllText("rocket.mos", string.Join("\n",
        "loadFile(\"../python/rocket.mo\");",
        "getErrorString();",
        "setDebugFlags(\"-disableDirectionalDerivatives\");",
        "buildModelFMU(rocket, version=\"2.0\", fmuType=\"me\");",
        "getErrorString();", ""));
    var psi = new ProcessStartInfo("omc", "rocket.mos") { RedirectStandardOutput = true };
    string[] @out;
    using (var omc = Process.Start(psi)) {
      @out = omc.StandardOutput.ReadToEnd().Split('\n', 5);
      omc.WaitForExit();
    }
    if (@out[0] != "true") {
      throw new Exception("Modelica compilation failed: " + @out[1]);
    }
    string fmu_file = @out[3].Replace("\"", "");
    string flag = @out.Length > 4 ? @out[4] : "";
    if (!fmu_file.EndsWith(".fmu")) throw new Exception("FMU generation failed: " + flag);
    Console.WriteLine("translateModelFMU warnings:\n" + flag);

    // Regardless of how the FMU was obtained, we must unzip it before CasADi can process it
    // If there is an existing directory with the same name it is better to explicitly
    // delete it to prevent potential memory corruption on overwriting.
    string unzipped_path = Path.GetFullPath("rocket_fmu_unzipped");
    if (Directory.Exists(unzipped_path)) Directory.Delete(unzipped_path, true);
    ZipFile.ExtractToDirectory("rocket.fmu", unzipped_path);

    // Create a CasADi/DaeBuilder instance from the unzipped FMU
    var dae = ca.DaeBuilder("rocket", unzipped_path);
    dae.disp(true);

    // Get state vector, initial conditions, bounds
    List<string> x = dae.x();
    List<double> lbx = dae.min(x);
    List<double> ubx = dae.max(x);
    List<double> x0 = dae.start(x);
    Console.WriteLine($"x:  [{string.Join(", ", x)}]");
    Console.WriteLine($"lbx:  [{string.Join(", ", lbx)}]");
    Console.WriteLine($"ubx:  [{string.Join(", ", ubx)}]");
    Console.WriteLine($"x0:  [{string.Join(", ", x0)}]");

    // Get control vector, initial conditions, bounds
    List<string> u = dae.u();
    List<double> lbu = dae.min(u);
    List<double> ubu = dae.max(u);
    List<double> u0 = dae.start(u);
    Console.WriteLine($"u:  [{string.Join(", ", u)}]");
    Console.WriteLine($"lbu:  [{string.Join(", ", lbu)}]");
    Console.WriteLine($"ubu:  [{string.Join(", ", ubu)}]");
    Console.WriteLine($"u0:  [{string.Join(", ", u0)}]");

    // Let's create a CasADi function for evaluating the ODE right-hand-side.
    // We only need to specify the expressions that are varying:
    Function f = dae.create("f", ["x", "u"], ["ode"]);

    // This is a standard CasADi function that can be embedded into CasADi expressions
    Console.WriteLine(f);

    // We can evaluate it and calculate derivatives as with other CasADi functions
    DM ode0 = f.Call(x0.ToArray(), u0.ToArray());
    Console.WriteLine($"ode0:  {ode0}");

    // Analytic first derivatives with sparsities are also available
    Function jac = f.factory("jac_f", ["x", "u"], ["jac:ode:x", "jac:ode:u"]);
    List<DM> jac_res = jac.call(new DM[] { x0.ToArray(), u0.ToArray() });
    DM jac_ode_x0 = jac_res[0], jac_ode_u0 = jac_res[1];
    Console.WriteLine($"jac_ode_x0:  {jac_ode_x0}");
    Console.WriteLine($"jac_ode_u0:  {jac_ode_u0}");

    // Code generation for the functions are currently work in progress
    f.generate("fmu_codegen");
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

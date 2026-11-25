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
using System.Linq;
using CasADi;

// Simulating a bouncing ball with DaeBuilder and event handling
// Joel Andersson, 2025

public static class bouncing_ball_daebuilder {
  public static void Main() {
    if (!ca.has_integrator("cvodes")) { Console.WriteLine("SKIP: integrator cvodes not available"); return; }
    if (!hasCommand("gcc")) { Console.WriteLine("SKIP: no C compiler (gcc) found"); return; }

    // Start with an empty DaeBuilder instance
    var dae = ca.DaeBuilder("bouncing_ball");

    // Model variables
    MX t = dae.add("t", "independent");
    MX h = dae.add("h", "output", new Dictionary<string, object> { ["start"] = 5.0, ["initial"] = "exact" });
    MX v = dae.add("v", "output", new Dictionary<string, object> { ["start"] = 0.0, ["initial"] = "exact" });

    // Dynamic equations
    dae.eq(dae.der(h), v);
    dae.eq(dae.der(v), -9.81);
    dae.disp(true);

    // Event dynamics: When h < 0, reinitialize v to -0.8*v
    dae.when(h < 0, [dae.reinit("v", -0.8 * dae.pre(v))]);
    dae.disp(true);

    // Default experiment
    dae.set_start_time(0);
    dae.set_stop_time(7);

    // Simulate
    double[] tgrid = Enumerable.Range(0, 100)
        .Select(i => dae.start_time() + (dae.stop_time() - dae.start_time()) * i / 99).ToArray();
    Function sim = ca.integrator("sim", "cvodes", dae.create(), 0, tgrid,
                              new Dictionary<string, object> { ["transition"] = dae.transition() });
    var simres = sim.call(new Dictionary<string, DM> { ["x0"] = dae.start(dae.x()).ToArray() });

    // Visualize the solution
    DM h_traj = simres["xf"][0, ..].T();
    // plot: h_traj over tgrid

    // Export FMU
    var fmu_files = dae.export_fmu();
    Console.WriteLine($"Generated files: {{{string.Join(", ", fmu_files.Select(e => $"{e.Key}={e.Value}"))}}}");

    // Compile DLL
    string fmi_headers = "../../../external_packages/FMI-Standard-3.0/headers";
    string cfiles = string.Join(" ", fmu_files.Keys.Where(f => f.EndsWith(".c")));
    string sofile = dae.name() + ".so";
    Process.Start("sh", ["-c", $"gcc --shared -fPIC -I{fmi_headers} {cfiles} -o {sofile}"]).WaitForExit();
    Console.WriteLine($"Compiled {sofile}");
    fmu_files[sofile] = "binaries/x86_64-linux/" + sofile;

    // Package into an FMU
    string fmuname = dae.name() + ".fmu";
    File.Delete(fmuname);  // overwrite, as Python's zipfile.ZipFile(fmuname, "w")
    using (ZipArchive fmufile = ZipFile.Open(fmuname, ZipArchiveMode.Create)) {
      foreach (var e in fmu_files) {
        fmufile.CreateEntryFromFile(e.Key, (string)e.Value);
        File.Delete(e.Key);
      }
    }
    Console.WriteLine($"Created FMU: {fmuname}");

    // Load the FMU in FMPy (a Python package, not available from C#)
    Console.WriteLine("FMPy not installed. Skipping FMU simulation.");
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

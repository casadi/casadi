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

// Example how to export FMUs from CasADi
// Joel Andersson, joel@jaeandersson.com

public static class fmu_export {
  public static void Main() {
    if (!hasCommand("gcc")) { Console.WriteLine("SKIP: no C compiler (gcc) found"); return; }

    // Start with an empty DaeBuilder instance
    var dae = ca.DaeBuilder("vdp");

    // States
    MX t = dae.add("t", "independent");
    MX x1 = dae.add("x1", "output", new Dictionary<string, object> { ["start"] = 1.0 });
    MX x2 = dae.add("x2", "output", new Dictionary<string, object> { ["start"] = 0.0 });
    MX q = dae.add("q");

    // Input
    MX u = dae.add("u", "input");

    // Set ODE right-hand-sides
    dae.eq(dae.der(x1), (1 - x2 * x2) * x1 - x2 + u);
    dae.eq(dae.der(x2), x1);
    dae.eq(dae.der(q), ca.pow(x1, 2) + ca.pow(x2, 2) + ca.pow(u, 2));

    // Add bounds
    dae.set_min("u", -0.75);
    dae.set_max("u", 1);
    dae.set_start("u", 0.5);

    // Print DAE
    dae.disp(true);

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

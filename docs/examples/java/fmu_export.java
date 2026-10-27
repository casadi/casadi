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
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;
import java.util.zip.ZipEntry;
import java.util.zip.ZipOutputStream;

// Example how to export FMUs from CasADi
// Joel Andersson, joel@jaeandersson.com

public class fmu_export {
  public static void main(String[] args) throws Exception {
    if (!hasCommand("gcc")) { System.out.println("SKIP: no C compiler (gcc) found"); return; }

    // Start with an empty DaeBuilder instance
    DaeBuilder dae = new DaeBuilder("vdp");

    // States
    MX t = dae.add("t", "independent");
    MX x1 = dae.add("x1", "output", Map.of("start", 1));
    MX x2 = dae.add("x2", "output", Map.of("start", 0));
    MX q = dae.add("q");

    // Input
    MX u = dae.add("u", "input");

    // Set ODE right-hand-sides
    dae.eq(dae.der(x1), ca.minus(1, x2.times(x2)).times(x1).minus(x2).plus(u));
    dae.eq(dae.der(x2), x1);
    dae.eq(dae.der(q), x1.power(2).plus(x2.power(2)).plus(u.power(2)));

    // Add bounds
    dae.set_min("u", -0.75);
    dae.set_max("u", 1);
    dae.set_start("u", 0.5);

    // Print DAE
    dae.disp(true);

    // Export FMU
    Map<String, Object> fmu_files = new LinkedHashMap<>(dae.export_fmu());
    System.out.println("Generated files: " + fmu_files);

    // Compile DLL
    String fmi_headers = "../../../external_packages/FMI-Standard-3.0/headers";
    String cfiles = fmu_files.keySet().stream().filter(f -> f.endsWith(".c")).collect(Collectors.joining(" "));
    String sofile = dae.name() + ".so";
    new ProcessBuilder("sh", "-c", "gcc --shared -fPIC -I" + fmi_headers + " " + cfiles + " -o " + sofile)
        .inheritIO().start().waitFor();
    System.out.println("Compiled " + sofile);
    fmu_files.put(sofile, "binaries/x86_64-linux/" + sofile);

    // Package into an FMU
    String fmuname = dae.name() + ".fmu";
    try (ZipOutputStream fmufile = new ZipOutputStream(Files.newOutputStream(Paths.get(fmuname)))) {
      for (Map.Entry<String, Object> e : fmu_files.entrySet()) {
        fmufile.putNextEntry(new ZipEntry((String) e.getValue()));
        Files.copy(Paths.get(e.getKey()), fmufile);
        fmufile.closeEntry();
        Files.delete(Paths.get(e.getKey()));
      }
    }
    System.out.println("Created FMU: " + fmuname);

    // Load the FMU in FMPy (a Python package, not available from Java)
    System.out.println("FMPy not installed. Skipping FMU simulation.");
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

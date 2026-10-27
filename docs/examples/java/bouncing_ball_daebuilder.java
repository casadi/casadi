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

// Simulating a bouncing ball with DaeBuilder and event handling
// Joel Andersson, 2025

public class bouncing_ball_daebuilder {
  public static void main(String[] args) throws Exception {
    if (!ca.has_integrator("cvodes")) { System.out.println("SKIP: integrator cvodes not available"); return; }
    if (!hasCommand("gcc")) { System.out.println("SKIP: no C compiler (gcc) found"); return; }

    // Start with an empty DaeBuilder instance
    DaeBuilder dae = new DaeBuilder("bouncing_ball");

    // Model variables
    MX t = dae.add("t", "independent");
    MX h = dae.add("h", "output", Map.of("start", 5, "initial", "exact"));
    MX v = dae.add("v", "output", Map.of("start", 0, "initial", "exact"));

    // Dynamic equations
    dae.eq(dae.der(h), v);
    dae.eq(dae.der(v), -9.81);
    dae.disp(true);

    // Event dynamics: When h < 0, reinitialize v to -0.8*v
    dae.when(ca.lt(h, 0), List.of(dae.reinit("v", ca.times(-0.8, dae.pre(v)))));
    dae.disp(true);

    // Default experiment
    dae.set_start_time(0);
    dae.set_stop_time(7);

    // Simulate
    double[] tgrid = new double[100];
    for (int i = 0; i < 100; i++) tgrid[i] = dae.start_time() + (dae.stop_time() - dae.start_time()) * i / 99;
    Function sim = ca.integrator("sim", "cvodes", dae.create(), 0, tgrid,
                                 Map.of("transition", dae.transition()));
    Map<String, GenericMatrix> simres = sim.invoke(Map.of("x0", dae.start(dae.x())));

    // Visualize the solution
    DM h_traj = ((DM) simres.get("xf")).get(0, new Slice()).T();
    // plot: h_traj over tgrid

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

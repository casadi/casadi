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
import java.io.InputStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.Comparator;
import java.util.List;
import java.util.stream.Stream;
import java.util.zip.ZipEntry;
import java.util.zip.ZipInputStream;

public class modelica_fmu_import {
  public static void main(String[] args) throws Exception {
    if (!hasCommand("omc")) { System.out.println("SKIP: OpenModelica compiler (omc) not found"); return; }

    // Consider the following simple Modelica file
    System.out.println(Files.readString(Paths.get("../python/rocket.mo")));

    // Different tools exist to compile Modelica models like these into FMUs.
    // The following code shows how we could do it using the open-source tool
    // OpenModelica (here its command line compiler omc, running a .mos script),
    // cf. the OpenModelica User's Guide.
    Files.writeString(Paths.get("rocket.mos"), String.join("\n",
        "loadFile(\"../python/rocket.mo\");",
        "getErrorString();",
        "setDebugFlags(\"-disableDirectionalDerivatives\");",
        "buildModelFMU(rocket, version=\"2.0\", fmuType=\"me\");",
        "getErrorString();", ""));
    Process omc = new ProcessBuilder("omc", "rocket.mos").redirectErrorStream(true).start();
    List<String> out = Arrays.asList(new String(omc.getInputStream().readAllBytes()).split("\n", 5));
    omc.waitFor();
    if (!out.get(0).equals("true")) {
      throw new RuntimeException("Modelica compilation failed: " + out.get(1));
    }
    String fmu_file = out.get(3).replace("\"", "");
    String flag = out.size() > 4 ? out.get(4) : "";
    if (!fmu_file.endsWith(".fmu")) throw new RuntimeException("FMU generation failed: " + flag);
    System.out.println("translateModelFMU warnings:\n" + flag);

    // Regardless of how the FMU was obtained, we must unzip it before CasADi can process it
    // If there is an existing directory with the same name it is better to explicitly
    // delete it to prevent potential memory corruption on overwriting.
    Path unzipped_path = Paths.get("rocket_fmu_unzipped").toAbsolutePath();
    if (Files.isDirectory(unzipped_path)) rmtree(unzipped_path);
    unzip(Paths.get("rocket.fmu"), unzipped_path);

    // Create a CasADi/DaeBuilder instance from the unzipped FMU
    DaeBuilder dae = new DaeBuilder("rocket", unzipped_path.toString());
    dae.disp(true);

    // Get state vector, initial conditions, bounds
    List<String> x = dae.x();
    double[] lbx = dae.min(x);
    double[] ubx = dae.max(x);
    double[] x0 = dae.start(x);
    System.out.println("x:  " + x);
    System.out.println("lbx:  " + Arrays.toString(lbx));
    System.out.println("ubx:  " + Arrays.toString(ubx));
    System.out.println("x0:  " + Arrays.toString(x0));

    // Get control vector, initial conditions, bounds
    List<String> u = dae.u();
    double[] lbu = dae.min(u);
    double[] ubu = dae.max(u);
    double[] u0 = dae.start(u);
    System.out.println("u:  " + u);
    System.out.println("lbu:  " + Arrays.toString(lbu));
    System.out.println("ubu:  " + Arrays.toString(ubu));
    System.out.println("u0:  " + Arrays.toString(u0));

    // Let's create a CasADi function for evaluating the ODE right-hand-side.
    // We only need to specify the expressions that are varying:
    Function f = dae.create("f", List.of("x", "u"), List.of("ode"));

    // This is a standard CasADi function that can be embedded into CasADi expressions
    System.out.println(f);

    // We can evaluate it and calculate derivatives as with other CasADi functions
    Object ode0 = f.invoke(x0, u0);
    System.out.println("ode0:  " + ode0);

    // Analytic first derivatives with sparsities are also available
    Function jac = f.factory("jac_f", List.of("x", "u"), List.of("jac:ode:x", "jac:ode:u"));
    List<?> jac_res = (List<?>) jac.invoke(x0, u0);
    Object jac_ode_x0 = jac_res.get(0), jac_ode_u0 = jac_res.get(1);
    System.out.println("jac_ode_x0:  " + jac_ode_x0);
    System.out.println("jac_ode_u0:  " + jac_ode_u0);

    // Code generation for the functions are currently work in progress
    f.generate("fmu_codegen");
  }

  static boolean hasCommand(String cmd) {
    try {
      return new ProcessBuilder(cmd, "--version").redirectErrorStream(true)
          .redirectOutput(ProcessBuilder.Redirect.DISCARD).start().waitFor() == 0;
    } catch (Exception e) {
      return false;
    }
  }

  static void rmtree(Path dir) throws Exception {
    try (Stream<Path> s = Files.walk(dir)) {
      for (Path p : s.sorted(Comparator.reverseOrder()).toList()) Files.delete(p);
    }
  }

  static void unzip(Path zip, Path dest) throws Exception {
    try (InputStream in = Files.newInputStream(zip); ZipInputStream z = new ZipInputStream(in)) {
      for (ZipEntry e; (e = z.getNextEntry()) != null; ) {
        Path p = dest.resolve(e.getName()).normalize();
        if (!p.startsWith(dest)) throw new RuntimeException("Bad zip entry: " + e.getName());
        if (e.isDirectory()) {
          Files.createDirectories(p);
        } else {
          Files.createDirectories(p.getParent());
          Files.copy(z, p);
        }
      }
    }
  }
}

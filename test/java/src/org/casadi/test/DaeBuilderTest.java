package org.casadi.test;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.DirectoryStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardCopyOption;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.List;
import java.util.Map;
import java.util.stream.Stream;
import java.util.zip.ZipEntry;
import java.util.zip.ZipInputStream;

import org.casadi.CasadiMeta;
import org.casadi.DM;
import org.casadi.DaeBuilder;
import org.casadi.Function;
import org.casadi.GenericMatrix;
import org.casadi.GlobalOptions;
import org.casadi.MX;
import org.casadi.SX;
import org.casadi.Slice;
import org.casadi.ca;
import org.junit.jupiter.api.Assumptions;
import org.junit.jupiter.api.Test;

/** Port of test/python/daebuilder.py */
public class DaeBuilderTest extends CasadiTestCase {

  private static boolean hasFeature(String f) {
    return CasadiMeta.feature_list().contains(f);
  }

  private static void assumeResource(String path) {
    Assumptions.assumeTrue(Files.exists(Paths.get(path)), "resource not available: " + path);
  }

  private static SX vertcatSX(SX... x) {
    return ca.vertcat(List.of(x));
  }

  private static MX vertcatMX(MX... x) {
    return ca.vertcat(List.of(x));
  }

  /** Outputs agree (helpers.py checkfunction_light) */
  private static void checkfunctionLight(Function f, Function g, List<DM> inputs, double tol) {
    List<DM> a = f.call(inputs), b = g.call(inputs);
    assertEquals(b.size(), a.size());
    for (int i = 0; i < a.size(); i++) checkarray(a.get(i), b.get(i), tol, "output " + i);
  }

  /** Serialization roundtrip evaluates the same (helpers.py check_serialize) */
  private static void checkSerialize(Function f, List<DM> inputs) {
    Function f2 = Function.deserialize(f.serialize());
    checkfunctionLight(f, f2, inputs, 1e-12);
  }

  private static Path unzip(String zip, Path dest) throws IOException {
    if (Files.isDirectory(dest)) rmtree(dest);
    try (ZipInputStream z = new ZipInputStream(Files.newInputStream(Paths.get(zip)))) {
      for (ZipEntry e; (e = z.getNextEntry()) != null; ) {
        Path p = dest.resolve(e.getName()).normalize();
        if (e.isDirectory()) {
          Files.createDirectories(p);
        } else {
          Files.createDirectories(p.getParent());
          Files.copy(z, p, StandardCopyOption.REPLACE_EXISTING);
        }
      }
    }
    return dest;
  }

  private static void rmtree(Path p) throws IOException {
    try (Stream<Path> s = Files.walk(p)) {
      for (Path q : (Iterable<Path>) s.sorted(Comparator.reverseOrder())::iterator) Files.delete(q);
    }
  }

  private static List<Path> glob(String pattern) throws IOException {
    List<Path> r = new ArrayList<>();
    try (DirectoryStream<Path> d = Files.newDirectoryStream(Paths.get("."), pattern)) {
      for (Path p : d) r.add(p);
    }
    return r;
  }

  /** Van der Pol reference ode for x = [x0; x1] */
  private static Function vdpRef(String name) {
    SX x0 = SX.sym("x0"), x1 = SX.sym("x1");
    SX x = vertcatSX(x0, x1);
    return new Function(name, List.of(x),
        List.of(vertcatSX(x1, ca.times(1.0, ca.times(ca.minus(1.0, ca.times(x0, x0)), x1)).minus(x0))));
  }

  @Test
  public void test_reference_fmus() {
    Assumptions.assumeTrue(hasFeature("ghc-filesystem"), "ghc-filesystem feature not available");
    for (String name : List.of("VanDerPol2", "VanDerPol3")) {
      String fmu_file = "../data/" + name + ".fmu";
      assumeResource(fmu_file);
      DaeBuilder dae = new DaeBuilder("car", fmu_file);
      dae.disp(true);
      Function f = dae.create("f", List.of("x"), List.of("ode"));
      Function f_ref = vdpRef("f");
      List<DM> test_point = dms(dm(1.1, 1.3));
      checkfunction(f, f_ref, test_point, 1e-7);
      checkSerialize(f, test_point);
    }
  }

  @Test
  public void test_fmu_export_compile_pack_roundtrip() throws IOException {
    Assumptions.assumeTrue(hasFeature("ghc-filesystem") && hasFeature("libzip"), "ghc-filesystem/libzip not available");
    String fmi_headers = "../../external_packages/FMI-Standard-3.0/headers";
    assumeResource(fmi_headers + "/fmi3Functions.h");
    Assumptions.assumeTrue(RUN_SLOW, "needs a C compiler (slow)");

    DaeBuilder dae = new DaeBuilder("vdp");
    MX x0 = dae.add("x0");
    MX x1 = dae.add("x1");
    dae.eq(dae.der(x0), x1);
    dae.eq(dae.der(x1), ca.minus(ca.times(ca.minus(1.0, ca.times(x0, x0)), x1), x0));
    dae.set_start("x0", 0);
    dae.set_start("x1", 0);

    Map<String, Object> files = dae.export_fmu();
    try {
      files = dae.compile_fmu(files, Map.of("include_dirs", List.of(fmi_headers)));
    } catch (RuntimeException e) {
      for (String local : files.keySet()) Files.deleteIfExists(Paths.get(local));
      Assumptions.abort("C compilation unavailable: " + e.getMessage());
    }
    String fmu_file = dae.pack_fmu(files, Map.of("path", "vdp_roundtrip.fmu"));
    assertTrue(Files.exists(Paths.get(fmu_file)));
    assertFalse(Files.exists(Paths.get("vdp.c")));

    DaeBuilder dae2 = new DaeBuilder("vdp", fmu_file);
    Function f = dae2.create("f", List.of("x"), List.of("ode"));
    checkfunction(f, vdpRef("f_ref"), dms(dm(1.1, 1.3)), 1e-7);

    f.close();
    dae2.close();
    Files.deleteIfExists(Paths.get(fmu_file));
  }

  @Test
  public void test_cstr() throws IOException {
    String fmu_file = "../data/cstr.fmu";
    assumeResource(fmu_file);
    Path unzipped = unzip(fmu_file, Files.createTempDirectory("cstr"));
    DaeBuilder dae = new DaeBuilder("cstr", unzipped.toString());
    Function f = dae.create("f", List.of("x", "u"), List.of("ode"));
    f.stats();
    dae.close();

    SX C_A = SX.sym("C_A"), C_B = SX.sym("C_B");
    double V = 1, k = 0.5;
    SX q_in = SX.sym("q_in"), C_A_in = SX.sym("C_A_in"), C_B_in = SX.sym("C_B_in");
    SX x = vertcatSX(C_A, C_B);
    SX u = vertcatSX(C_A_in, C_B_in, q_in);
    SX r = ca.times(V * k, ca.times(C_A, C_B));
    SX ode = vertcatSX(ca.times(q_in, C_A_in.minus(C_A)).minus(r).div(V), ca.times(q_in, C_B_in.minus(C_B)).plus(r).div(V));
    Function f_ref = new Function("f", List.of(x, u), List.of(ode), List.of("x", "u"), List.of("ode"));

    List<DM> test_point = dms(dm(1.1, 1.3), dm(1.7, 1.11, 1.13));
    checkfunction(f, f_ref, test_point, 1e-4);
    System.out.println(f.stats());
    rmtree(unzipped);
  }

  @Test
  public void test_fmi2_independent_var_order() throws IOException {
    for (String time_name : List.of("time", "t")) {
      for (Integer time_position : Arrays.asList(null, 0, 2, 4)) {
        String sub = "time_name=" + time_name + " time_position=" + time_position;
        List<String> names = new ArrayList<>(List.of("u", "x", "dx", "y"));
        if (time_position != null) names.add(time_position, time_name);
        StringBuilder xml = new StringBuilder();
        xml.append("<fmiModelDescription fmiVersion=\"2.0\" modelName=\"test\" guid=\"test\" numberOfEventIndicators=\"0\">");
        xml.append("<ModelExchange modelIdentifier=\"test\" /><ModelVariables>");
        for (String name : names) {
          String causality = name.equals("u") ? "input" : name.equals("y") ? "output" : name.equals(time_name) ? "independent" : "local";
          xml.append("<ScalarVariable name=\"" + name + "\" valueReference=\"" + names.indexOf(name) + "\" causality=\""
              + causality + "\" variability=\"continuous\">");
          if (name.equals("dx")) xml.append("<Real derivative=\"" + (names.indexOf("x") + 1) + "\" />");
          else if (name.equals("u") || name.equals("x")) xml.append("<Real start=\"0\" />");
          else xml.append("<Real />");
          xml.append("</ScalarVariable>");
        }
        xml.append("</ModelVariables><ModelStructure>");
        for (String[] gn : new String[][] {{"Outputs", "y"}, {"Derivatives", "dx"}}) {
          xml.append("<" + gn[0] + "><Unknown index=\"" + (names.indexOf(gn[1]) + 1) + "\" dependencies=\""
              + (names.indexOf("x") + 1) + "\" /></" + gn[0] + ">");
        }
        xml.append("</ModelStructure></fmiModelDescription>");
        Path directory = Files.createTempDirectory("daebuilder");
        try {
          Files.write(directory.resolve("modelDescription.xml"), xml.toString().getBytes(StandardCharsets.UTF_8));
          DaeBuilder dae = new DaeBuilder("test", directory.toString());
          assertEquals(List.of("x"), dae.x(), sub);
          assertEquals(List.of("u"), dae.u(), sub);
          assertEquals(List.of("y"), dae.y(), sub);
          assertEquals(List.of("dx", "y"), dae.w(), sub);
          assertEquals(time_position == null ? List.of("time") : List.of(time_name), dae.all("t"), sub);
          assertEquals(List.of("dx"), dae.der(), sub);
          dae = new DaeBuilder("test", directory.toString(), Map.of("ignore_time", true));
          assertEquals(List.of(), dae.all("t"), sub);
        } finally {
          rmtree(directory);
        }
      }
    }
  }

  @Test
  public void test_indendent_var() throws IOException {
    String os = System.getProperty("os.name").toLowerCase();
    Assumptions.assumeTrue(os.contains("linux"), "FMU binary was produced on linux");
    String fmu_file = "../data/car_t_fmu2.fmu";
    assumeResource(fmu_file);
    Path unzipped = unzip(fmu_file, Files.createTempDirectory("car"));
    DaeBuilder dae = new DaeBuilder("car", unzipped.toString());
    dae.disp(true);

    assertEquals(List.of("u_susp", "v_road", "x_road"), dae.u());
    assertEquals(List.of("x_body_squared"), dae.y());
    assertEquals(List.of("ContState1", "ContState2", "ContState3", "ContState4"), dae.x());
    assertEquals(List.of("Fc_susp", "d_susp", "d_tire", "k_susp", "k_tire", "m_body", "m_susp"), dae.p());
    assertEquals(List.of("x_body_squared", "Deriv1", "Deriv2", "Deriv3", "Deriv4"), dae.w());
    Function f = dae.create("f", List.of("t", "u", "x"), List.of("ode", "y"));
    SX t = SX.sym("t");

    SX u_susp = SX.sym("u_susp"), v_road = SX.sym("v_road"), x_road = SX.sym("x_road");
    SX u = vertcatSX(u_susp, v_road, x_road);

    SX x_body = SX.sym("x_body"), x_susp = SX.sym("x_susp"), v_body = SX.sym("v_body"), v_susp = SX.sym("v_susp");
    SX x = vertcatSX(x_body, x_susp, v_body, v_susp);

    SX Fc_susp = SX.sym("Fc_susp"), d_susp = SX.sym("d_susp"), d_tire = SX.sym("d_tire"), k_susp = SX.sym("k_susp"),
        k_tire = SX.sym("k_tire"), m_body = SX.sym("m_body"), m_susp = SX.sym("m_susp");
    SX p = vertcatSX(Fc_susp, d_susp, d_tire, k_susp, k_tire, m_body, m_susp);
    DM p_num = dm(100, 350, 15020, 80000, 0.5 * 500000, 2500, 320);

    SX road_force = k_tire.times(x_road.minus(x_susp)).plus(d_tire.times(v_road.minus(v_susp)));
    SX x_body_der = v_body;
    SX x_susp_der = v_susp;
    SX inner = d_susp.times(v_body.minus(v_susp)).plus(Fc_susp.times(ca.tanh(ca.times(10.0, v_body.minus(v_susp)))))
        .plus(k_susp.times(x_body.minus(x_susp)));
    SX v_body_der = ca.rdivide(1.0, m_body).times(inner.unaryMinus().plus(u_susp)).plus(t);
    SX v_susp_der = ca.rdivide(1.0, m_susp).times(inner.plus(road_force).minus(u_susp));

    SX x_der = vertcatSX(x_body_der, x_susp_der, v_body_der, v_susp_der);
    SX y = ca.power(x_body, 2.0);

    Function f_ref = new Function("f", List.of(t, u, x), List.of(ca.substitute(x_der, p, new SX(p_num)), ca.substitute(y, p, new SX(p_num))));

    DM u_num = dm(0.1, 0.3, 1.7);
    DM x_num = dm(1.8, 1.5, 1.9, 1.9);

    checkfunctionLight(f, f_ref, dms(new DM(0), u_num, x_num), 1e-10);
    checkfunctionLight(f, f_ref, dms(new DM(0.5), u_num, x_num), 1e-10);

    f = dae.create("f", List.of("t", "u", "x", "p"), List.of("ode", "y"));
    f_ref = new Function("f", List.of(t, u, x, p), List.of(x_der, y));

    checkfunctionLight(f, f_ref, dms(new DM(0.5), u_num, x_num, p_num.times(1.1)), 1e-10);
    SX w = vertcatSX(SX.sym("x_body_squared"), SX.sym("Deriv1"), SX.sym("Deriv2"), SX.sym("Deriv3"), SX.sym("Deriv4"));
    f_ref = new Function("f", List.of(t, u, x, p, w), List.of(w.get(new Slice(1, 5)), w.get(0)));

    assertInException("Unsupported input: 'w'", () -> dae.create("f", List.of("t", "u", "x", "p", "w"), List.of("ode", "y")));
    rmtree(unzipped);
  }

  @Test
  public void test_rumoca() throws IOException, InterruptedException {
    Assumptions.assumeTrue(hasFeature("rumoca"), "rumoca feature not available");
    String rumoca = Paths.get(GlobalOptions.getCasadiPath(), "rumoca").toString();
    String rumoca_exe = Paths.get(GlobalOptions.getCasadiPath(), "rumoca.exe").toString();
    if (Files.exists(Paths.get(rumoca)) || Files.exists(Paths.get(rumoca_exe))) {
      new ProcessBuilder(rumoca, "-t", "../assets/casadi_daebuilder.jinja", "-m", "../assets/hello_world.mo").inheritIO().start().waitFor();
    }
  }

  @Test
  public void test_fmu_zip() throws IOException {
    String fmu_file = "../data/VanDerPol2.fmu";
    assumeResource(fmu_file);
    String unzipped_name = "VanDerPol2";
    Integer unzipped_count = null;
    for (String serialize_mode : List.of("link", "embed")) {
      DaeBuilder dae;
      if (hasFeature("ghc-filesystem")) {
        unzipped_count = glob(unzipped_name + "*.unzipped").size();
        dae = new DaeBuilder("cstr", fmu_file, Map.of("resource_serialize_mode", serialize_mode));
        assertEquals(unzipped_count + 1, glob(unzipped_name + "*.unzipped").size());
      } else {
        assertInException("passing fmu files to DaeBuilder is unsupported", () -> new DaeBuilder("cstr", fmu_file));
        continue;
      }
      dae.disp(true);
      Function f = dae.create("f", List.of("x"), List.of("ode"));
      dae.close();

      Function f_ref = vdpRef("f");
      List<DM> test_point = dms(dm(1.1, 1.3));

      checkfunctionLight(f, f_ref, test_point, 1e-7);
      f.save("f.casadi");
      f.close();
      f = Function.load("f.casadi");
      checkfunctionLight(f, f_ref, test_point, 1e-7);
      // Type decay, so test twice
      f.save("f.casadi");
      checkfunctionLight(f, f_ref, test_point, 1e-7);
      f.close();
      f = Function.load("f.casadi");
      checkfunction(f, f_ref, test_point, 1e-7);
      f.close();
      assertEquals((int) unzipped_count, glob(unzipped_name + "*.unzipped").size());
    }
    Files.deleteIfExists(Paths.get("f.casadi"));
  }

  /** Coefficients (ascending powers) of the product of two polynomials */
  private static double[] polymul(double[] a, double[] b) {
    double[] r = new double[a.length + b.length - 1];
    for (int i = 0; i < a.length; i++)
      for (int j = 0; j < b.length; j++) r[i + j] += a[i] * b[j];
    return r;
  }

  private static double polyval(double[] a, double x) {
    double r = 0;
    for (int i = a.length - 1; i >= 0; i--) r = r * x + a[i];
    return r;
  }

  private static DM dmOf(GenericMatrix x) {
    return (DM) x;
  }

  @SuppressWarnings("unchecked")
  @Test
  public void test_fmu_demo() {
    Assumptions.assumeTrue(hasRootfinder("kinsol") && hasRootfinder("newton") && hasNlpsol("ipopt") && hasIntegrator("cvodes"),
        "needs kinsol, newton, ipopt, cvodes");
    String fmu_file = "../data/vdp.fmu";
    assumeResource(fmu_file);
    DaeBuilder dae = new DaeBuilder("vdp", fmu_file);
    dae.disp(true);

    System.out.println("FMU (" + fmu_file + ") provides directional derivatives: " + dae.provides_directional_derivative());

    List<String> x = dae.x();
    double[] lbx = dae.min(x), ubx = dae.max(x), x0 = dae.start(x);
    System.out.println("x: " + x + " lbx: " + Arrays.toString(lbx) + " ubx: " + Arrays.toString(ubx) + " x0: " + Arrays.toString(x0));

    List<String> u = dae.u();
    double[] lbu = dae.min(u), ubu = dae.max(u);
    System.out.println("u: " + u + " lbu: " + Arrays.toString(lbu) + " ubu: " + Arrays.toString(ubu));

    Function f = dae.create("f", List.of("x", "u"), List.of("ode"));
    System.out.println(f);
    double u_test = 0.4;
    Object xdot_test = f.invoke(x0, u_test);
    System.out.println("xdot_test: " + xdot_test);

    dae.set("u", 0.4);
    Function f_no_u = dae.create("f_no_u", List.of("x"), List.of("ode"));
    System.out.println(f_no_u);
    System.out.println("xdot_test_no_u: " + f_no_u.invoke(x0));

    Function f_with_aux = dae.create("f_with_aux", List.of("x", "u"), List.of("ode"), Map.of("aux", List.of("x1", "x2", "u")));
    System.out.println(f_with_aux);
    System.out.println("xdot_test_with_aux: " + f_with_aux.invoke(x0, u_test));
    System.out.println(f_with_aux.stats());

    Function J = dae.create("J", List.of("x", "u"), List.of("ode", "jac_ode_x", "jac_ode_u"), Map.of("verbose", true, "parallelization", "openmp"));
    System.out.println(J);
    List<Object> r3 = (List<Object>) J.invoke(x0, u_test);
    System.out.println("xdot_test: " + r3.get(0) + " A_test: " + r3.get(1) + " B_test: " + r3.get(2));

    Function J_fd = dae.create("J_fd", List.of("x", "u"), List.of("ode", "jac_ode_x", "jac_ode_u"),
        Map.of("verbose", true, "parallelization", "openmp", "enable_ad", false));
    System.out.println(J_fd);
    r3 = (List<Object>) J_fd.invoke(x0, u_test);
    System.out.println("xdot_fd: " + r3.get(0) + " A_fd: " + r3.get(1) + " B_fd: " + r3.get(2));

    Function H = dae.create("H", List.of("x", "u", "adj_ode"), List.of("jac_adj_x_x", "jac_adj_x_u", "jac_adj_u_u"),
        Map.of("verbose", true, "parallelization", "openmp"));
    System.out.println(H);
    r3 = (List<Object>) H.invoke(x0, u_test, 1);
    System.out.println("H_xx: " + r3.get(0) + " H_xu: " + r3.get(1) + " H_uu: " + r3.get(2));

    double T = 10.;
    double[] tgrid = new double[100];
    for (int i = 0; i < 100; i++) tgrid[i] = T * i / 99;
    dae.set_min(x.get(0), Double.NEGATIVE_INFINITY);
    Function daefun = dae.create("daefun");
    Function sim = ca.integrator("sim", "cvodes", daefun, 0, tgrid);

    Map<String, GenericMatrix> r = sim.invoke(Map.of("x0", x0, "u", u_test));
    DM x_traj_sim = dmOf(r.get("xf"));

    Function fwd_sim = ca.integrator("fwd_sim", "cvodes", daefun, 0, tgrid,
        Map.of("nfwd", 1, "verbose", false, "print_stats", true, "second_order_correction", false));
    System.out.println(fwd_sim);
    DM seed_x0 = dm(1, 0);
    Map<String, GenericMatrix> fwd_test = fwd_sim.invoke(Map.of("x0", ca.horzcat(List.of(dm(x0), seed_x0)),
        "u", ca.horzcat(List.of(new DM(u_test), new DM(0)))));
    DM fxf = dmOf(fwd_test.get("xf"));
    System.out.println("dxf/dx0 * seed_x0 [AD, forward] = " + fxf.get(new Slice(), new Slice(1, fxf.size2(), 2)));
    double pert = 1e-3;
    Map<String, GenericMatrix> pert_sim = sim.invoke(Map.of("x0", dm(x0).plus(seed_x0.times(pert)), "u", u_test));
    DM x_pert_sim = dmOf(pert_sim.get("xf"));
    System.out.println("dxf/dx0 * seed_x0 [FD] = " + x_pert_sim.minus(x_traj_sim).div(pert));

    Function adj_sim = ca.integrator("adj_sim", "cvodes", daefun, 0, tgrid,
        Map.of("nadj", 1, "verbose", false, "print_stats", true, "second_order_correction", true, "print_time", true));
    System.out.println(adj_sim);
    DM seed_xf = DM.zeros(2, 100);
    seed_xf.set(new Slice(), new DM(-1), dm(1, 0));
    Map<String, GenericMatrix> adj_test = adj_sim.invoke(Map.of("x0", x0, "adj_xf", seed_xf, "u", u_test));
    System.out.println("trans(dxf/du) * seed_xf (AD, adjoint) = " + adj_test.get("adj_u"));

    Function fwd_adj_sim = ca.integrator("fwd_adj_sim", "cvodes", daefun, 0, tgrid,
        Map.of("nfwd", 1, "nadj", 1, "verbose", false, "print_stats", true, "second_order_correction", false));
    System.out.println(fwd_adj_sim);
    Map<String, GenericMatrix> fwd_adj_test = fwd_adj_sim.invoke(Map.of("x0", ca.horzcat(List.of(dm(x0), seed_x0)),
        "adj_xf", ca.horzcat(List.of(seed_xf, DM.zeros(seed_xf.size1(), seed_xf.size2()))),
        "u", ca.horzcat(List.of(new DM(u_test), new DM(0)))));
    DM fau = dmOf(fwd_adj_test.get("adj_u"));
    System.out.println("d(trans(dxf/du) * seed_xf)/dx0 * seed_x0 (AD, forward-over-adjoint) = "
        + fau.get(new Slice(), new Slice(1, fau.size2(), 2)));

    int N = 50;
    double h = T / N;

    int d = 3;
    double[] tau_root = new double[d + 1];
    System.arraycopy(ca.collocation_points(d, "legendre"), 0, tau_root, 1, d);
    double[][] C = new double[d + 1][d + 1];
    double[] D = new double[d + 1];
    double[] B = new double[d + 1];
    for (int j = 0; j < d + 1; j++) {
      double[] p = {1};
      for (int rr = 0; rr < d + 1; rr++)
        if (rr != j) {
          double den = tau_root[j] - tau_root[rr];
          p = polymul(p, new double[] {-tau_root[rr] / den, 1 / den});
        }
      D[j] = polyval(p, 1.0);
      double[] pder = new double[Math.max(1, p.length - 1)];
      for (int i = 1; i < p.length; i++) pder[i - 1] = i * p[i];
      for (int rr = 0; rr < d + 1; rr++) C[j][rr] = polyval(pder, tau_root[rr]);
      double[] pint = new double[p.length + 1];
      for (int i = 0; i < p.length; i++) pint[i + 1] = p[i] / (i + 1);
      B[j] = polyval(pint, 1.0);
    }
    System.out.println("tau_root: " + Arrays.toString(tau_root) + " C: " + Arrays.deepToString(C) + " D: "
        + Arrays.toString(D) + " B: " + Arrays.toString(B));
    DM Ddm = dm(D);
    int nx = x.size(), nu = u.size();

    MX U = MX.sym("U", nu);
    List<MX> X = new ArrayList<>();
    for (int j = 0; j < d + 1; j++) X.add(MX.sym("X_" + j, nx));
    List<MX> g = new ArrayList<>();
    for (int j = 1; j < d + 1; j++) {
      MX xdot_j = new MX(0);
      for (int rr = 0; rr < d + 1; rr++) xdot_j = xdot_j.plus(ca.times(C[rr][j], X.get(rr)));
      g.add(ca.times(h, (MX) f.invoke(X.get(j), U)).minus(xdot_j));
    }
    MX gg = ca.vertcat(g);
    MX X_unknown = ca.vertcat(X.subList(1, d + 1));
    Function rfp = new Function("rfp", List.of(X_unknown, X.get(0), U), List.of(gg), List.of("V", "X0", "U"), List.of("res"));

    Function ifcn = ca.rootfinder("ifcn", "kinsol", rfp, Map.of("print_level", 1));
    DM v0 = ca.repmat(dm(x0), d, 1);
    v0 = (DM) ifcn.invoke(v0, x0, u_test);
    DM x_all = ca.horzcat(List.of(dm(x0), ca.reshape(v0, nx, d)));
    DM xf = ca.mtimes(x_all, Ddm);
    System.out.println("x_all: " + x_all + " xf: " + xf);

    ifcn = ca.rootfinder("ifcn", "newton", rfp, Map.of("print_iteration", true));
    v0 = ca.repmat(dm(x0), d, 1);
    v0 = (DM) ifcn.invoke(v0, x0, u_test);
    x_all = ca.horzcat(List.of(dm(x0), ca.reshape(v0, nx, d)));
    xf = ca.mtimes(x_all, Ddm);
    System.out.println("x_all: " + x_all + " xf: " + xf);

    System.out.println(rfp);
    ifcn = ca.rootfinder("ifcn", "kinsol", rfp, Map.of("print_level", 1));
    System.out.println(ifcn);
    Function jac_ifcn = ifcn.jacobian();
    Map<String, GenericMatrix> sol = jac_ifcn.invoke(Map.of("V0", v0, "X0", x0, "U", u_test));
    System.out.println("jac_V_U = " + sol.get("jac_V_U") + " jac_V_X0 = " + sol.get("jac_V_X0"));

    DM v0_pert = (DM) ifcn.invoke(ca.repmat(dm(x0), d, 1), x0, u_test + pert);
    System.out.println("jac_V_U (FD) = " + v0_pert.minus(v0).div(pert));
    v0_pert = (DM) ifcn.invoke(ca.repmat(dm(x0), d, 1), dm(x0).plus(dm(pert, 0)), u_test);
    DM jac_V_X0_fd1 = v0_pert.minus(v0).div(pert);
    v0_pert = (DM) ifcn.invoke(ca.repmat(dm(x0), d, 1), dm(x0).plus(dm(0, pert)), u_test);
    DM jac_V_X0_fd2 = v0_pert.minus(v0).div(pert);
    System.out.println("jac_V_X0 (FD) = " + ca.horzcat(List.of(jac_V_X0_fd1, jac_V_X0_fd2)));

    ifcn = ca.rootfinder("ifcn", "kinsol", rfp);
    MX x0_in = MX.sym("x0", nx);
    List<MX> u_all = new ArrayList<>();
    MX xk = x0_in;
    List<MX> x_all_l = new ArrayList<>(List.of(x0_in));
    for (int k = 0; k < N; k++) {
      MX uk = MX.sym("u" + k, nu);
      u_all.add(uk);
      MX vk = (MX) ifcn.invoke(ca.repmat(xk, d, 1), xk, uk);
      vk = ca.horzcat(List.of(xk, ca.reshape(vk, nx, d)));
      xk = ca.mtimes(vk, Ddm);
      x_all_l.add(xk);
    }
    Function irksim = new Function("irksim", List.of(x0_in, ca.vertcat(u_all)), List.of(ca.horzcat(x_all_l)),
        List.of("x0", "u"), List.of("x"));
    DM x_traj = (DM) irksim.invoke(x0, u_test);

    Function fwd_irksim = irksim.forward(1);
    Map<String, GenericMatrix> fwd_t = fwd_irksim.invoke(Map.of("x0", x0, "u", u_test, "fwd_x0", seed_x0));
    System.out.println("fwd_x (AD, forward) = " + fwd_t.get("fwd_x"));
    DM x_pert = (DM) irksim.invoke(dm(x0).plus(seed_x0.times(pert)), u_test);
    System.out.println("fwd_x (FD) = " + x_pert.minus(x_traj).div(pert));

    Function adj_irksim = irksim.reverse(1);
    DM seed_x = DM.zeros(2, N + 1);
    seed_x.set(new DM(1), new DM(N), new DM(1));
    Map<String, GenericMatrix> adj_t = adj_irksim.invoke(Map.of("x0", x0, "u", u_test, "adj_x", seed_x));
    System.out.println("adj_x0 (AD, reverse) = " + adj_t.get("adj_x0") + " adj_u (AD, reverse) = " + adj_t.get("adj_u"));
    System.out.println("Forward mode: derivative of x[1, N] w.r.t. x0[0] = " + dmOf(fwd_t.get("fwd_x")).get(1, N));
    System.out.println("Reverse mode: derivative of x[1, N] w.r.t. x0[0] = " + dmOf(adj_t.get("adj_x0")).get(0));

    Function jac_irksim = irksim.jacobian();
    Map<String, GenericMatrix> jac_test = jac_irksim.invoke(Map.of("x0", x0, "u", u_test));
    System.out.println("jac_x_x0 = " + jac_test.get("jac_x_x0") + " jac_x_u = " + jac_test.get("jac_x_u"));

    List<MX> w = new ArrayList<>();
    g = new ArrayList<>();
    MX Jobj = new MX(0);
    List<MX> U_all = new ArrayList<>(), X_all = new ArrayList<>();
    xk = new MX(dm(x0));
    for (int k = 0; k < N; k++) {
      MX uk = MX.sym("u" + k, nu);
      U_all.add(uk);
      w.add(uk);
      List<MX> Xc = new ArrayList<>(List.of(xk));
      for (int j = 1; j < d + 1; j++) {
        MX xkj = MX.sym("x" + k + "_" + j, nx);
        Xc.add(xkj);
        X_all.add(xkj);
        w.add(xkj);
      }
      for (int j = 1; j < d + 1; j++) {
        MX xdot_j = new MX(0);
        for (int rr = 0; rr < d + 1; rr++) xdot_j = xdot_j.plus(ca.times(C[rr][j], Xc.get(rr)));
        g.add(ca.times(h, (MX) f.invoke(Xc.get(j), uk)).minus(xdot_j));
      }
      MX xk_end = new MX(0);
      for (int j = 0; j < d + 1; j++) xk_end = xk_end.plus(ca.times(D[j], Xc.get(j)));
      xk = MX.sym("x" + (k + 1), nx);
      X_all.add(xk);
      w.add(xk);
      g.add(xk_end.minus(xk));
      Jobj = Jobj.plus(ca.sumsqr(uk)).plus(ca.sumsqr(xk));
    }
    MX wv = ca.vertcat(w), gv = ca.vertcat(g), X_allm = ca.horzcat(X_all), U_allm = ca.horzcat(U_all);
    Function to_w = new Function("to_w", List.of(X_allm, U_allm), List.of(wv), List.of("x", "u"), List.of("w"));
    Function from_w = new Function("from_w", List.of(wv), List.of(X_allm, U_allm), List.of("w"), List.of("x", "u"));
    System.out.println(from_w);
    Map<String, MX> nlp = Map.of("x", wv, "f", Jobj, "g", gv);

    long t0 = System.nanoTime();
    MX jac_g_w = ca.jacobian(gv, wv);
    System.out.println("Symbolic expression for Jacobian of g w.r.t w formed in " + (System.nanoTime() - t0) * 1e-9 + " s");
    Function jfcn = new Function("jfcn", List.of(wv), List.of(jac_g_w), List.of("w"), List.of("jac_g_w"));
    GenericMatrix w_test = to_w.invoke(Map.of("u", u_test, "x", x0)).get("w");
    jfcn.invoke(w_test);

    Map<String, Object> opts = Map.of("ipopt", Map.of("hessian_approximation", "exact", "linear_solver", "mumps"));
    Function solver = ca.nlpsol("solver", "ipopt", nlp, opts);

    GenericMatrix w0 = to_w.invoke(Map.of("u", 0, "x", x0)).get("w");
    GenericMatrix lbw = to_w.invoke(Map.of("u", lbu, "x", lbx)).get("w");
    GenericMatrix ubw = to_w.invoke(Map.of("u", ubu, "x", ubx)).get("w");
    solver.invoke(Map.of("x0", w0, "lbx", lbw, "ubx", ubw, "lbg", 0, "ubg", 0));
  }

  @Test
  public void test_stats_available_bug() {
    String fmu_file = "../data/VanDerPol2.fmu";
    assumeResource(fmu_file);
    DaeBuilder dae = new DaeBuilder("vdp", fmu_file);
    dae.get("x1");
    dae.disp(true);
  }

  @SuppressWarnings("unchecked")
  @Test
  public void test_w_fmu() throws IOException {
    String fmu_file = "../data/VanDerPol2.fmu";
    assumeResource(fmu_file);
    DaeBuilder dae = new DaeBuilder("vdp", fmu_file);
    dae.get("x1");
    dae.disp(true);

    DM.rng(1);
    DM x0 = DM.rand(2);
    DM u0 = DM.rand();

    Function f = dae.create("f", List.of("x", "u"), List.of("ode"), Map.of("aux", dae.w()));

    DM r1 = (DM) f.invoke(x0, u0);
    Map<String, Object> a1 = (Map<String, Object>) f.stats().get("aux");
    DM s1 = ca.vertcat(List.of(auxDM(a1.get("der(x0)")), auxDM(a1.get("der(x1)"))));

    f.save("f.casadi");
    Function f2 = Function.load("f.casadi");
    DM r2 = (DM) f2.invoke(x0, u0);
    Map<String, Object> a2 = (Map<String, Object>) f2.stats().get("aux");
    DM s2 = ca.vertcat(List.of(auxDM(a2.get("der(x0)")), auxDM(a2.get("der(x1)"))));

    checkarray(r1, r2);
    checkarray(s1, s2);
    Files.deleteIfExists(Paths.get("f.casadi"));
  }

  private static DM auxDM(Object o) {
    if (o instanceof Number) return new DM(((Number) o).doubleValue());
    if (o instanceof double[]) return new DM((double[]) o);
    if (o instanceof List) return new DM(((List<?>) o).stream().mapToDouble(v -> ((Number) v).doubleValue()).toArray());
    return (DM) o;
  }

  @Test
  public void test_w() throws IOException {
    DaeBuilder dae = new DaeBuilder("test");

    dae.disp(true);

    MX u = dae.add("u", "input", "continuous");
    double u_val = 1.1;

    MX c0 = dae.add("c0", "parameter", "fixed");
    double c0_val = 0.97;
    dae.set("c0", c0_val);
    dae.add("c1", "output", "constant");
    double c1_val = 0.17;
    dae.set("c1", c1_val);
    dae.add("c2", "local", "constant");
    double c2_val = 0.19;
    dae.set("c2", c2_val);
    DM c_val = dm(c0_val, c1_val, c2_val);

    dae.get("c2");
    checkarray(new DM(((Number) dae.get("c2")).doubleValue()), c2_val);

    MX p = dae.add("p", "parameter", "tunable");
    double p_val = 1.3;
    dae.set("p", p_val);

    MX x0 = dae.add("x0", "output", "continuous");
    double x0_val = 1.7;
    MX x1 = dae.add("x1", "local", "continuous");
    double x1_val = 1.9;

    dae.eq(dae.der(x0), x0.times(x1));
    dae.eq(dae.der(x1), x1.plus(u).plus(p.times(c0)));

    MX w0 = dae.add("w0", "output", "continuous");
    MX w1 = dae.add("w1", "local", "continuous");

    dae.eq(w0, ca.sin(x1));
    dae.eq(w1, ca.power(w0, 2.0).plus(u));

    MX x2 = dae.add("x2", "local", "continuous");
    double x2_val = 1.11;
    dae.eq(dae.der(x2), ca.times(3.0, w1));

    DM x_val = dm(x0_val, x1_val, x2_val);

    MX z0 = dae.add("z0", "output", "continuous");
    MX z1 = dae.add("z1", "local", "continuous");
    double z0_val = 2.11;
    double z1_val = 3.13;
    DM z_val = dm(z0_val, z1_val);

    dae.eq(0, x1.minus(ca.times(2.0, z0)));
    dae.eq(0, x1.minus(ca.times(3.0, z1)).plus(z0));

    dae.disp(true);

    Function f = dae.create("f", List.of("p", "u", "x"), List.of("y"));

    f.disp(true);

    assertInException("Initialization failed since variables [p, c0] are free",
        () -> dae.dependent_fun("f", List.of("u", "x"), List.of("w")));
    f = dae.dependent_fun("f", List.of("c", "p", "u", "x"), List.of("w"));
    System.out.println(f);

    double s = Math.sin(x1_val);
    DM w_val = (DM) f.invoke(c_val, p_val, u_val, x_val);
    checkarray(w_val, dm(x0_val * x1_val, x1_val + u_val + p_val * c0_val, s, s * s + u_val, 3 * (s * s + u_val)));
    f = dae.create("f", List.of("u", "x"), List.of("ode"));
    double CORRECTION = 0; // Looks like a bug
    checkarray((DM) f.invoke(u_val, x_val), dm(x0_val * x1_val, x1_val + u_val + CORRECTION * p_val * c0_val, 3 * (s * s + u_val)));

    f = dae.create("f", List.of("p", "u", "x"), List.of("ode"));
    checkarray((DM) f.invoke(p_val, u_val, x_val), dm(x0_val * x1_val, x1_val + u_val + p_val * c0_val, 3 * (s * s + u_val)));

    f = dae.create("f", List.of("p", "u", "x", "w"), List.of("ode"));
    checkarray((DM) f.invoke(p_val, u_val, x_val, w_val), dm(x0_val * x1_val, x1_val + u_val + p_val * c0_val, 3 * (s * s + u_val)));
    CORRECTION = 1; // Why is this not 2?
    checkarray((DM) f.invoke(p_val, u_val, x_val, w_val.times(2.0)),
        dm(CORRECTION * x0_val * x1_val, CORRECTION * (x1_val + u_val + p_val * c0_val), 2 * 3 * (s * s + u_val)));
    f = dae.create("f", List.of("p", "u", "x", "z"), List.of("alg"));
    checkarray((DM) f.invoke(p_val, u_val, x_val, z_val), dm(-(x1_val - (2. * z0_val)), -((x1_val - (3 * z1_val)) + z0_val)));
    f = dae.create("f", List.of("p", "u", "x", "z"), List.of("y"));
    checkarray((DM) f.invoke(p_val, u_val, x_val, z_val), dm(c1_val, x0_val, s, z0_val));
    Map<String, Object> files = dae.export_fmu();
    System.out.println(files);
    // python leaves the exported sources in cwd; clean them up here
    for (String local : files.keySet()) Files.deleteIfExists(Paths.get(local));
  }

  @Test
  public void test_basemodelica() throws IOException {
    Assumptions.assumeTrue(ca.has_modelicaparser("lacemodelica"), "modelicaparser lacemodelica not available");
    DaeBuilder dae_bmo = new DaeBuilder("from_bmo", "../data/DerVariable.bmo");

    DaeBuilder dae_sym = new DaeBuilder("symbolic");
    MX x = dae_sym.add("x", "local", "continuous");
    MX y = dae_sym.add("y", "local", "continuous");
    MX dx_dt = dae_sym.add("dx_dt", "local", "continuous");
    MX dy_dt = dae_sym.add("dy_dt", "local", "continuous");
    dae_sym.eq(dae_sym.der(x), ca.times(-0.5, x).plus(y));
    dae_sym.eq(dae_sym.der(y), x.minus(ca.times(0.3, y)));
    dae_sym.eq(dx_dt, dae_sym.der(x));
    dae_sym.eq(dy_dt, dae_sym.der(y));

    Function f_bmo = dae_bmo.create("f_bmo", List.of("x"), List.of("ode"));
    Function f_sym = dae_sym.create("f_sym", List.of("x"), List.of("ode"));

    for (double[] x_test : new double[][] {{1.0, 2.0}, {0.0, 0.0}, {2.5, 1.5}, {-1.0, 3.0}})
      checkfunctionLight(f_bmo, f_sym, dms(dm(x_test)), 1e-10);

    List<Path> temp_dirs = glob("DerVariable.bmo.*.temp");
    assertEquals(1, temp_dirs.size());
    Path temp_dir = temp_dirs.get(0);

    dae_bmo.close();
    dae_sym.close();
    f_bmo.close();
    f_sym.close();

    assertFalse(Files.exists(temp_dir), "Temp dir not cleaned up: " + temp_dir);
  }

}

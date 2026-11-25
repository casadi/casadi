using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.IO.Compression;
using System.Linq;
using System.Runtime.InteropServices;

namespace CasADi.Tests {

/// <summary>Port of test/java DaeBuilderTest (itself a port of test/python/daebuilder.py)</summary>
public class DaeBuilderTest : CasadiTestCase {

  /// <summary>test/data (not copied into WSL: reached through /mnt/c)</summary>
  static string DataDir =>
    Environment.GetEnvironmentVariable("CASADI_TEST_DATA") ?? "/mnt/c/Users/Gebruiker/casadi/test/data";

  static string data(string name) => Path.Combine(DataDir, name);

  static bool hasFeature(string f) => CasadiMeta.feature_list().Contains(f);

  static void assumeResource(string path) =>
    Check.Assume(File.Exists(path) || Directory.Exists(path), "resource not available: " + path);

  /// <summary>Outputs agree (helpers.py checkfunction_light)</summary>
  static void checkfunctionLight(Function f, Function g, List<DM> inputs, double tol) {
    List<DM> a = f.call(inputs), b = g.call(inputs);
    assertEquals(b.Count, a.Count);
    for (int i = 0; i < a.Count; i++) checkarray(a[i], b[i], tol, "output " + i);
  }

  /// <summary>Serialization roundtrip evaluates the same (helpers.py check_serialize)</summary>
  static void checkSerialize(Function f, List<DM> inputs) {
    Function f2 = Function.deserialize(f.serialize());
    checkfunctionLight(f, f2, inputs, 1e-12);
  }

  static string unzip(string zip, string dest) {
    if (Directory.Exists(dest)) Directory.Delete(dest, true);
    ZipFile.ExtractToDirectory(zip, dest);
    return dest;
  }

  static string tempDir(string prefix) {
    string d = Path.Combine(Path.GetTempPath(), prefix + Guid.NewGuid().ToString("N"));
    Directory.CreateDirectory(d);
    return d;
  }

  static string[] glob(string pattern) =>
    Directory.GetFileSystemEntries(".", pattern);

  static bool listEq(IEnumerable<string> expected, IEnumerable<string> actual) => expected.SequenceEqual(actual);

  static void assertList(IEnumerable<string> expected, IEnumerable<string> actual, string msg = "") =>
    assertTrue(listEq(expected, actual),
      msg + " expected [" + string.Join(", ", expected) + "], got [" + string.Join(", ", actual) + "]");

  /// <summary>Van der Pol reference ode for x = [x0; x1]</summary>
  static Function vdpRef(string name) {
    SX x0 = SX.sym("x0"), x1 = SX.sym("x1");
    SX x = ca.vertcat(x0, x1);
    return ca.Function(name, new[] { x }, new[] { ca.vertcat(x1, 1.0 * ((1.0 - x0 * x0) * x1) - x0) });
  }

  public void test_reference_fmus() {
    Check.Assume(hasFeature("ghc-filesystem"), "ghc-filesystem feature not available");
    foreach (string name in new[] { "VanDerPol2", "VanDerPol3" }) {
      string fmu_file = data(name + ".fmu");
      assumeResource(fmu_file);
      var dae = ca.DaeBuilder("car", fmu_file);
      dae.disp(true);
      Function f = dae.create("f", new[] { "x" }, new[] { "ode" });
      Function f_ref = vdpRef("f");
      List<DM> test_point = dms(ca.vertcat(1.1, 1.3));
      checkfunction(f, f_ref, test_point, 1e-7);
      checkSerialize(f, test_point);
    }
  }

  public void test_fmu_export_compile_pack_roundtrip() {
    Check.Assume(hasFeature("ghc-filesystem") && hasFeature("libzip"), "ghc-filesystem/libzip not available");
    string fmi_headers = "/mnt/c/Users/Gebruiker/casadi/external_packages/FMI-Standard-3.0/headers";
    assumeResource(fmi_headers + "/fmi3Functions.h");
    Check.Assume(RunSlow, "needs a C compiler (slow)");

    var dae = ca.DaeBuilder("vdp");
    MX x0 = dae.add("x0");
    MX x1 = dae.add("x1");
    dae.eq(dae.der(x0), x1);
    dae.eq(dae.der(x1), (1.0 - x0 * x0) * x1 - x0);
    dae.set_start("x0", 0);
    dae.set_start("x1", 0);

    Dictionary<string, object> files = dae.export_fmu();
    try {
      files = dae.compile_fmu(files, new Dictionary<string, object> { ["include_dirs"] = new List<string> { fmi_headers } });
    } catch (SwigException e) {
      foreach (string local in files.Keys) File.Delete(local);
      Check.Skip("C compilation unavailable: " + e.Message);
    }
    string fmu_file = dae.pack_fmu(files, new Dictionary<string, object> { ["path"] = "vdp_roundtrip.fmu" });
    assertTrue(File.Exists(fmu_file));
    assertFalse(File.Exists("vdp.c"));

    var dae2 = ca.DaeBuilder("vdp", fmu_file);
    Function f = dae2.create("f", new[] { "x" }, new[] { "ode" });
    checkfunction(f, vdpRef("f_ref"), dms(ca.vertcat(1.1, 1.3)), 1e-7);

    f.Dispose();
    dae2.Dispose();
    File.Delete(fmu_file);
  }

  public void test_cstr() {
    string fmu_file = data("cstr.fmu");
    assumeResource(fmu_file);
    string unzipped = unzip(fmu_file, tempDir("cstr"));
    var dae = ca.DaeBuilder("cstr", unzipped);
    Function f = dae.create("f", new[] { "x", "u" }, new[] { "ode" });
    f.stats();
    dae.Dispose();

    SX C_A = SX.sym("C_A"), C_B = SX.sym("C_B");
    double V = 1, k = 0.5;
    SX q_in = SX.sym("q_in"), C_A_in = SX.sym("C_A_in"), C_B_in = SX.sym("C_B_in");
    SX x = ca.vertcat(C_A, C_B);
    SX u = ca.vertcat(C_A_in, C_B_in, q_in);
    SX r = V * k * (C_A * C_B);
    SX ode = ca.vertcat((q_in * (C_A_in - C_A) - r) / V, (q_in * (C_B_in - C_B) + r) / V);
    var f_ref = ca.Function("f", new[] { x, u }, new[] { ode }, new[] { "x", "u" }, new[] { "ode" });

    List<DM> test_point = dms(ca.vertcat(1.1, 1.3), ca.vertcat(1.7, 1.11, 1.13));
    checkfunction(f, f_ref, test_point, 1e-4);
    Console.WriteLine(string.Join(", ", f.stats().Select(kv => kv.Key + ": " + kv.Value)));
    Directory.Delete(unzipped, true);
  }

  public void test_fmi2_independent_var_order() {
    foreach (string time_name in new[] { "time", "t" }) {
      foreach (int? time_position in new int?[] { null, 0, 2, 4 }) {
        string sub = "time_name=" + time_name + " time_position=" + time_position;
        var names = new List<string> { "u", "x", "dx", "y" };
        if (time_position != null) names.Insert(time_position.Value, time_name);
        var xml = new System.Text.StringBuilder();
        xml.Append("<fmiModelDescription fmiVersion=\"2.0\" modelName=\"test\" guid=\"test\" numberOfEventIndicators=\"0\">");
        xml.Append("<ModelExchange modelIdentifier=\"test\" /><ModelVariables>");
        foreach (string name in names) {
          string causality = name == "u" ? "input" : name == "y" ? "output" : name == time_name ? "independent" : "local";
          xml.Append("<ScalarVariable name=\"" + name + "\" valueReference=\"" + names.IndexOf(name) + "\" causality=\""
              + causality + "\" variability=\"continuous\">");
          if (name == "dx") xml.Append("<Real derivative=\"" + (names.IndexOf("x") + 1) + "\" />");
          else if (name == "u" || name == "x") xml.Append("<Real start=\"0\" />");
          else xml.Append("<Real />");
          xml.Append("</ScalarVariable>");
        }
        xml.Append("</ModelVariables><ModelStructure>");
        foreach (string[] gn in new[] { new[] { "Outputs", "y" }, new[] { "Derivatives", "dx" } }) {
          xml.Append("<" + gn[0] + "><Unknown index=\"" + (names.IndexOf(gn[1]) + 1) + "\" dependencies=\""
              + (names.IndexOf("x") + 1) + "\" /></" + gn[0] + ">");
        }
        xml.Append("</ModelStructure></fmiModelDescription>");
        string directory = tempDir("daebuilder");
        try {
          File.WriteAllText(Path.Combine(directory, "modelDescription.xml"), xml.ToString());
          var dae = ca.DaeBuilder("test", directory);
          assertList(new[] { "x" }, dae.x(), sub);
          assertList(new[] { "u" }, dae.u(), sub);
          assertList(new[] { "y" }, dae.y(), sub);
          assertList(new[] { "dx", "y" }, dae.w(), sub);
          assertList(time_position == null ? new[] { "time" } : new[] { time_name }, dae.all("t"), sub);
          assertList(new[] { "dx" }, dae.der(), sub);
          dae = ca.DaeBuilder("test", directory, new Dictionary<string, object> { ["ignore_time"] = true });
          assertList(new string[0], dae.all("t"), sub);
        } finally {
          Directory.Delete(directory, true);
        }
      }
    }
  }

  public void test_indendent_var() {
    Check.Assume(RuntimeInformation.IsOSPlatform(OSPlatform.Linux), "FMU binary was produced on linux");
    string fmu_file = data("car_t_fmu2.fmu");
    assumeResource(fmu_file);
    string unzipped = unzip(fmu_file, tempDir("car"));
    var dae = ca.DaeBuilder("car", unzipped);
    dae.disp(true);

    assertList(new[] { "u_susp", "v_road", "x_road" }, dae.u());
    assertList(new[] { "x_body_squared" }, dae.y());
    assertList(new[] { "ContState1", "ContState2", "ContState3", "ContState4" }, dae.x());
    assertList(new[] { "Fc_susp", "d_susp", "d_tire", "k_susp", "k_tire", "m_body", "m_susp" }, dae.p());
    assertList(new[] { "x_body_squared", "Deriv1", "Deriv2", "Deriv3", "Deriv4" }, dae.w());
    Function f = dae.create("f", new[] { "t", "u", "x" }, new[] { "ode", "y" });
    SX t = SX.sym("t");

    SX u_susp = SX.sym("u_susp"), v_road = SX.sym("v_road"), x_road = SX.sym("x_road");
    SX u = ca.vertcat(u_susp, v_road, x_road);

    SX x_body = SX.sym("x_body"), x_susp = SX.sym("x_susp"), v_body = SX.sym("v_body"), v_susp = SX.sym("v_susp");
    SX x = ca.vertcat(x_body, x_susp, v_body, v_susp);

    SX Fc_susp = SX.sym("Fc_susp"), d_susp = SX.sym("d_susp"), d_tire = SX.sym("d_tire"), k_susp = SX.sym("k_susp"),
        k_tire = SX.sym("k_tire"), m_body = SX.sym("m_body"), m_susp = SX.sym("m_susp");
    SX p = ca.vertcat(Fc_susp, d_susp, d_tire, k_susp, k_tire, m_body, m_susp);
    DM p_num = ca.vertcat(100, 350, 15020, 80000, 0.5 * 500000, 2500, 320);

    SX road_force = k_tire * (x_road - x_susp) + d_tire * (v_road - v_susp);
    SX x_body_der = v_body;
    SX x_susp_der = v_susp;
    SX inner = d_susp * (v_body - v_susp) + Fc_susp * ca.tanh(10.0 * (v_body - v_susp)) + k_susp * (x_body - x_susp);
    SX v_body_der = 1.0 / m_body * (-inner + u_susp) + t;
    SX v_susp_der = 1.0 / m_susp * (inner + road_force - u_susp);

    SX x_der = ca.vertcat(x_body_der, x_susp_der, v_body_der, v_susp_der);
    SX y = ca.pow(x_body, 2.0);

    var f_ref = ca.Function("f", new[] { t, u, x }, new[] { ca.substitute(x_der, p, (SX)p_num), ca.substitute(y, p, (SX)p_num) });

    DM u_num = ca.vertcat(0.1, 0.3, 1.7);
    DM x_num = ca.vertcat(1.8, 1.5, 1.9, 1.9);

    checkfunctionLight(f, f_ref, dms(0, u_num, x_num), 1e-10);
    checkfunctionLight(f, f_ref, dms(0.5, u_num, x_num), 1e-10);

    f = dae.create("f", new[] { "t", "u", "x", "p" }, new[] { "ode", "y" });
    f_ref = ca.Function("f", new[] { t, u, x, p }, new[] { x_der, y });

    checkfunctionLight(f, f_ref, dms(0.5, u_num, x_num, p_num * 1.1), 1e-10);
    SX w = ca.vertcat(SX.sym("x_body_squared"), SX.sym("Deriv1"), SX.sym("Deriv2"), SX.sym("Deriv3"), SX.sym("Deriv4"));
    f_ref = ca.Function("f", new[] { t, u, x, p, w }, new[] { w[1..5], w[0] });

    assertInException("Unsupported input: 'w'", () => dae.create("f", new[] { "t", "u", "x", "p", "w" }, new[] { "ode", "y" }));
    Directory.Delete(unzipped, true);
  }

  public void test_rumoca() {
    Check.Assume(hasFeature("rumoca"), "rumoca feature not available");
    string rumoca = Path.Combine(GlobalOptions.getCasadiPath(), "rumoca");
    string rumoca_exe = Path.Combine(GlobalOptions.getCasadiPath(), "rumoca.exe");
    if (File.Exists(rumoca) || File.Exists(rumoca_exe)) {
      using var p = Process.Start(rumoca, new[] { "-t", "/mnt/c/Users/Gebruiker/casadi/test/assets/casadi_daebuilder.jinja",
        "-m", "/mnt/c/Users/Gebruiker/casadi/test/assets/hello_world.mo" });
      p.WaitForExit();
    }
  }

  public void test_fmu_zip() {
    string fmu_file = data("VanDerPol2.fmu");
    assumeResource(fmu_file);
    string unzipped_name = "VanDerPol2";
    int unzipped_count = 0;
    foreach (string serialize_mode in new[] { "link", "embed" }) {
      DaeBuilder dae;
      if (hasFeature("ghc-filesystem")) {
        unzipped_count = glob(unzipped_name + "*.unzipped").Length;
        dae = ca.DaeBuilder("cstr", fmu_file, new Dictionary<string, object> { ["resource_serialize_mode"] = serialize_mode });
        assertEquals(unzipped_count + 1, glob(unzipped_name + "*.unzipped").Length);
      } else {
        assertInException("passing fmu files to DaeBuilder is unsupported", () => ca.DaeBuilder("cstr", fmu_file));
        continue;
      }
      dae.disp(true);
      Function f = dae.create("f", new[] { "x" }, new[] { "ode" });
      dae.Dispose();

      Function f_ref = vdpRef("f");
      List<DM> test_point = dms(ca.vertcat(1.1, 1.3));

      checkfunctionLight(f, f_ref, test_point, 1e-7);
      f.save("f.casadi");
      f.Dispose();
      f = Function.load("f.casadi");
      checkfunctionLight(f, f_ref, test_point, 1e-7);
      // Type decay, so test twice
      f.save("f.casadi");
      checkfunctionLight(f, f_ref, test_point, 1e-7);
      f.Dispose();
      f = Function.load("f.casadi");
      checkfunction(f, f_ref, test_point, 1e-7);
      f.Dispose();
      // the unzipped directory goes away with the last reference: let finalizer-queued deletes run
      GC.Collect();
      GC.WaitForPendingFinalizers();
      SwigRT.RunQueuedDeletes();
      assertEquals(unzipped_count, glob(unzipped_name + "*.unzipped").Length);
    }
    File.Delete("f.casadi");
  }

  /// <summary>Coefficients (ascending powers) of the product of two polynomials</summary>
  static double[] polymul(double[] a, double[] b) {
    double[] r = new double[a.Length + b.Length - 1];
    for (int i = 0; i < a.Length; i++)
      for (int j = 0; j < b.Length; j++) r[i + j] += a[i] * b[j];
    return r;
  }

  static double polyval(double[] a, double x) {
    double r = 0;
    for (int i = a.Length - 1; i >= 0; i--) r = r * x + a[i];
    return r;
  }

  public void test_fmu_demo() {
    Check.Assume(hasRootfinder("kinsol") && hasRootfinder("newton") && hasNlpsol("ipopt") && hasIntegrator("cvodes"),
        "needs kinsol, newton, ipopt, cvodes");
    string fmu_file = data("vdp.fmu");
    assumeResource(fmu_file);
    var dae = ca.DaeBuilder("vdp", fmu_file);
    dae.disp(true);

    Console.WriteLine("FMU (" + fmu_file + ") provides directional derivatives: " + dae.provides_directional_derivative());

    List<string> x = dae.x();
    double[] lbx = dae.min(x).ToArray(), ubx = dae.max(x).ToArray(), x0 = dae.start(x).ToArray();
    Console.WriteLine("x: " + string.Join(",", x) + " lbx: " + string.Join(",", lbx) + " ubx: " + string.Join(",", ubx)
      + " x0: " + string.Join(",", x0));

    List<string> u = dae.u();
    double[] lbu = dae.min(u).ToArray(), ubu = dae.max(u).ToArray();
    Console.WriteLine("u: " + string.Join(",", u) + " lbu: " + string.Join(",", lbu) + " ubu: " + string.Join(",", ubu));

    Function f = dae.create("f", new[] { "x", "u" }, new[] { "ode" });
    Console.WriteLine(f);
    double u_test = 0.4;
    DM xdot_test = f.Call(x0, u_test);
    Console.WriteLine("xdot_test: " + xdot_test);

    dae.set("u", 0.4);
    Function f_no_u = dae.create("f_no_u", new[] { "x" }, new[] { "ode" });
    Console.WriteLine(f_no_u);
    Console.WriteLine("xdot_test_no_u: " + f_no_u.Call(x0));

    Function f_with_aux = dae.create("f_with_aux", new[] { "x", "u" }, new[] { "ode" },
      new Dictionary<string, object> { ["aux"] = new List<string> { "x1", "x2", "u" } });
    Console.WriteLine(f_with_aux);
    Console.WriteLine("xdot_test_with_aux: " + f_with_aux.Call(x0, u_test));
    Console.WriteLine(f_with_aux.stats());

    Function J = dae.create("J", new[] { "x", "u" }, new[] { "ode", "jac_ode_x", "jac_ode_u" },
      new Dictionary<string, object> { ["verbose"] = true, ["parallelization"] = "openmp" });
    Console.WriteLine(J);
    List<DM> r3 = J.call(new DM[] { x0, u_test });
    Console.WriteLine("xdot_test: " + r3[0] + " A_test: " + r3[1] + " B_test: " + r3[2]);

    Function J_fd = dae.create("J_fd", new[] { "x", "u" }, new[] { "ode", "jac_ode_x", "jac_ode_u" },
      new Dictionary<string, object> { ["verbose"] = true, ["parallelization"] = "openmp", ["enable_ad"] = false });
    Console.WriteLine(J_fd);
    r3 = J_fd.call(new DM[] { x0, u_test });
    Console.WriteLine("xdot_fd: " + r3[0] + " A_fd: " + r3[1] + " B_fd: " + r3[2]);

    Function H = dae.create("H", new[] { "x", "u", "adj_ode" }, new[] { "jac_adj_x_x", "jac_adj_x_u", "jac_adj_u_u" },
      new Dictionary<string, object> { ["verbose"] = true, ["parallelization"] = "openmp" });
    Console.WriteLine(H);
    r3 = H.call(new DM[] { x0, u_test, 1 });
    Console.WriteLine("H_xx: " + r3[0] + " H_xu: " + r3[1] + " H_uu: " + r3[2]);

    double T = 10.0;
    double[] tgrid = new double[100];
    for (int i = 0; i < 100; i++) tgrid[i] = T * i / 99;
    dae.set_min(x[0], double.NegativeInfinity);
    Function daefun = dae.create("daefun");
    Function sim = ca.integrator("sim", "cvodes", daefun, 0, tgrid);

    var r = sim.call(new Dictionary<string, DM> { ["x0"] = x0, ["u"] = u_test });
    DM x_traj_sim = r["xf"];

    Function fwd_sim = ca.integrator("fwd_sim", "cvodes", daefun, 0, tgrid, new Dictionary<string, object> {
      ["nfwd"] = 1, ["verbose"] = false, ["print_stats"] = true, ["second_order_correction"] = false });
    Console.WriteLine(fwd_sim);
    DM seed_x0 = ca.vertcat(1, 0);
    var fwd_test = fwd_sim.call(new Dictionary<string, DM> { ["x0"] = ca.horzcat(x0, seed_x0), ["u"] = ca.horzcat(u_test, 0) });
    DM fxf = fwd_test["xf"];
    Console.WriteLine("dxf/dx0 * seed_x0 [AD, forward] = " + fxf[.., ca.Slice(1, fxf.size2(), 2)]);
    double pert = 1e-3;
    var pert_sim = sim.call(new Dictionary<string, DM> { ["x0"] = (DM)x0 + seed_x0 * pert, ["u"] = u_test });
    DM x_pert_sim = pert_sim["xf"];
    Console.WriteLine("dxf/dx0 * seed_x0 [FD] = " + (x_pert_sim - x_traj_sim) / pert);

    Function adj_sim = ca.integrator("adj_sim", "cvodes", daefun, 0, tgrid, new Dictionary<string, object> {
      ["nadj"] = 1, ["verbose"] = false, ["print_stats"] = true, ["second_order_correction"] = true, ["print_time"] = true });
    Console.WriteLine(adj_sim);
    DM seed_xf = DM.zeros(2, 100);
    seed_xf[.., seed_xf.size2() - 1] = ca.vertcat(1, 0);
    var adj_test = adj_sim.call(new Dictionary<string, DM> { ["x0"] = x0, ["adj_xf"] = seed_xf, ["u"] = u_test });
    Console.WriteLine("trans(dxf/du) * seed_xf (AD, adjoint) = " + adj_test["adj_u"]);

    Function fwd_adj_sim = ca.integrator("fwd_adj_sim", "cvodes", daefun, 0, tgrid, new Dictionary<string, object> {
      ["nfwd"] = 1, ["nadj"] = 1, ["verbose"] = false, ["print_stats"] = true, ["second_order_correction"] = false });
    Console.WriteLine(fwd_adj_sim);
    var fwd_adj_test = fwd_adj_sim.call(new Dictionary<string, DM> {
      ["x0"] = ca.horzcat(x0, seed_x0),
      ["adj_xf"] = ca.horzcat(seed_xf, DM.zeros(seed_xf.size1(), seed_xf.size2())),
      ["u"] = ca.horzcat(u_test, 0) });
    DM fau = fwd_adj_test["adj_u"];
    Console.WriteLine("d(trans(dxf/du) * seed_xf)/dx0 * seed_x0 (AD, forward-over-adjoint) = "
        + fau[.., ca.Slice(1, fau.size2(), 2)]);

    int N = 50;
    double h = T / N;

    int d = 3;
    double[] tau_root = new double[d + 1];
    ca.collocation_points(d, "legendre").CopyTo(tau_root, 1);
    double[,] C = new double[d + 1, d + 1];
    double[] D = new double[d + 1];
    double[] B = new double[d + 1];
    for (int j = 0; j < d + 1; j++) {
      double[] p = { 1 };
      for (int rr = 0; rr < d + 1; rr++)
        if (rr != j) {
          double den = tau_root[j] - tau_root[rr];
          p = polymul(p, new[] { -tau_root[rr] / den, 1 / den });
        }
      D[j] = polyval(p, 1.0);
      double[] pder = new double[Math.Max(1, p.Length - 1)];
      for (int i = 1; i < p.Length; i++) pder[i - 1] = i * p[i];
      for (int rr = 0; rr < d + 1; rr++) C[j, rr] = polyval(pder, tau_root[rr]);
      double[] pint = new double[p.Length + 1];
      for (int i = 0; i < p.Length; i++) pint[i + 1] = p[i] / (i + 1);
      B[j] = polyval(pint, 1.0);
    }
    Console.WriteLine("tau_root: " + string.Join(",", tau_root) + " C: " + (DM)C + " D: " + string.Join(",", D)
      + " B: " + string.Join(",", B));
    DM Ddm = D;
    int nx = x.Count, nu = u.Count;

    MX U = MX.sym("U", nu);
    var X = new List<MX>();
    for (int j = 0; j < d + 1; j++) X.Add(MX.sym("X_" + j, nx));
    var g = new List<MX>();
    for (int j = 1; j < d + 1; j++) {
      MX xdot_j = 0;
      for (int rr = 0; rr < d + 1; rr++) xdot_j = xdot_j + C[rr, j] * X[rr];
      g.Add(h * f.Call(X[j], U) - xdot_j);
    }
    MX gg = ca.vertcat(g);
    MX X_unknown = ca.vertcat(X.GetRange(1, d));
    var rfp = ca.Function("rfp", new[] { X_unknown, X[0], U }, new[] { gg }, new[] { "V", "X0", "U" }, new[] { "res" });

    Function ifcn = ca.rootfinder("ifcn", "kinsol", rfp, new Dictionary<string, object> { ["print_level"] = 1 });
    DM v0 = ca.repmat((DM)x0, d, 1);
    v0 = ifcn.Call(v0, x0, u_test);
    DM x_all = ca.horzcat(x0, ca.reshape(v0, nx, d));
    DM xf = ca.mtimes(x_all, Ddm);
    Console.WriteLine("x_all: " + x_all + " xf: " + xf);

    ifcn = ca.rootfinder("ifcn", "newton", rfp, new Dictionary<string, object> { ["print_iteration"] = true });
    v0 = ca.repmat((DM)x0, d, 1);
    v0 = ifcn.Call(v0, x0, u_test);
    x_all = ca.horzcat(x0, ca.reshape(v0, nx, d));
    xf = ca.mtimes(x_all, Ddm);
    Console.WriteLine("x_all: " + x_all + " xf: " + xf);

    Console.WriteLine(rfp);
    ifcn = ca.rootfinder("ifcn", "kinsol", rfp, new Dictionary<string, object> { ["print_level"] = 1 });
    Console.WriteLine(ifcn);
    Function jac_ifcn = ifcn.jacobian();
    var sol = jac_ifcn.call(new Dictionary<string, DM> { ["V0"] = v0, ["X0"] = x0, ["U"] = u_test });
    Console.WriteLine("jac_V_U = " + sol["jac_V_U"] + " jac_V_X0 = " + sol["jac_V_X0"]);

    DM v0_pert = ifcn.Call(ca.repmat((DM)x0, d, 1), x0, u_test + pert);
    Console.WriteLine("jac_V_U (FD) = " + (v0_pert - v0) / pert);
    v0_pert = ifcn.Call(ca.repmat((DM)x0, d, 1), (DM)x0 + ca.vertcat(pert, 0), u_test);
    DM jac_V_X0_fd1 = (v0_pert - v0) / pert;
    v0_pert = ifcn.Call(ca.repmat((DM)x0, d, 1), (DM)x0 + ca.vertcat(0, pert), u_test);
    DM jac_V_X0_fd2 = (v0_pert - v0) / pert;
    Console.WriteLine("jac_V_X0 (FD) = " + ca.horzcat(jac_V_X0_fd1, jac_V_X0_fd2));

    ifcn = ca.rootfinder("ifcn", "kinsol", rfp);
    MX x0_in = MX.sym("x0", nx);
    var u_all = new List<MX>();
    MX xk = x0_in;
    var x_all_l = new List<MX> { x0_in };
    for (int k = 0; k < N; k++) {
      MX uk = MX.sym("u" + k, nu);
      u_all.Add(uk);
      MX vk = ifcn.Call(ca.repmat(xk, d, 1), xk, uk);
      vk = ca.horzcat(xk, ca.reshape(vk, nx, d));
      xk = ca.mtimes(vk, (MX)Ddm);
      x_all_l.Add(xk);
    }
    var irksim = ca.Function("irksim", new[] { x0_in, ca.vertcat(u_all) }, new[] { ca.horzcat(x_all_l) },
        new[] { "x0", "u" }, new[] { "x" });
    DM x_traj = irksim.Call(x0, u_test);

    Function fwd_irksim = irksim.forward(1);
    var fwd_t = fwd_irksim.call(new Dictionary<string, DM> { ["x0"] = x0, ["u"] = u_test, ["fwd_x0"] = seed_x0 });
    Console.WriteLine("fwd_x (AD, forward) = " + fwd_t["fwd_x"]);
    DM x_pert = irksim.Call((DM)x0 + seed_x0 * pert, u_test);
    Console.WriteLine("fwd_x (FD) = " + (x_pert - x_traj) / pert);

    Function adj_irksim = irksim.reverse(1);
    DM seed_x = DM.zeros(2, N + 1);
    seed_x[1, N] = 1;
    var adj_t = adj_irksim.call(new Dictionary<string, DM> { ["x0"] = x0, ["u"] = u_test, ["adj_x"] = seed_x });
    Console.WriteLine("adj_x0 (AD, reverse) = " + adj_t["adj_x0"] + " adj_u (AD, reverse) = " + adj_t["adj_u"]);
    Console.WriteLine("Forward mode: derivative of x[1, N] w.r.t. x0[0] = " + fwd_t["fwd_x"][1, N]);
    Console.WriteLine("Reverse mode: derivative of x[1, N] w.r.t. x0[0] = " + adj_t["adj_x0"][0]);

    Function jac_irksim = irksim.jacobian();
    var jac_test = jac_irksim.call(new Dictionary<string, DM> { ["x0"] = x0, ["u"] = u_test });
    Console.WriteLine("jac_x_x0 = " + jac_test["jac_x_x0"] + " jac_x_u = " + jac_test["jac_x_u"]);

    var w = new List<MX>();
    g = new List<MX>();
    MX Jobj = 0;
    var U_all = new List<MX>();
    var X_all = new List<MX>();
    xk = (DM)x0;
    for (int k = 0; k < N; k++) {
      MX uk = MX.sym("u" + k, nu);
      U_all.Add(uk);
      w.Add(uk);
      var Xc = new List<MX> { xk };
      for (int j = 1; j < d + 1; j++) {
        MX xkj = MX.sym("x" + k + "_" + j, nx);
        Xc.Add(xkj);
        X_all.Add(xkj);
        w.Add(xkj);
      }
      for (int j = 1; j < d + 1; j++) {
        MX xdot_j = 0;
        for (int rr = 0; rr < d + 1; rr++) xdot_j = xdot_j + C[rr, j] * Xc[rr];
        g.Add(h * f.Call(Xc[j], uk) - xdot_j);
      }
      MX xk_end = 0;
      for (int j = 0; j < d + 1; j++) xk_end = xk_end + D[j] * Xc[j];
      xk = MX.sym("x" + (k + 1), nx);
      X_all.Add(xk);
      w.Add(xk);
      g.Add(xk_end - xk);
      Jobj = Jobj + ca.sumsqr(uk) + ca.sumsqr(xk);
    }
    MX wv = ca.vertcat(w), gv = ca.vertcat(g), X_allm = ca.horzcat(X_all), U_allm = ca.horzcat(U_all);
    var to_w = ca.Function("to_w", new[] { X_allm, U_allm }, new[] { wv }, new[] { "x", "u" }, new[] { "w" });
    var from_w = ca.Function("from_w", new[] { wv }, new[] { X_allm, U_allm }, new[] { "w" }, new[] { "x", "u" });
    Console.WriteLine(from_w);
    var nlp = new Dictionary<string, MX> { ["x"] = wv, ["f"] = Jobj, ["g"] = gv };

    var sw = Stopwatch.StartNew();
    MX jac_g_w = ca.jacobian(gv, wv);
    Console.WriteLine("Symbolic expression for Jacobian of g w.r.t w formed in " + sw.Elapsed.TotalSeconds + " s");
    var jfcn = ca.Function("jfcn", new[] { wv }, new[] { jac_g_w }, new[] { "w" }, new[] { "jac_g_w" });
    DM w_test = to_w.call(new Dictionary<string, DM> { ["u"] = u_test, ["x"] = x0 })["w"];
    jfcn.Call(w_test);

    var opts = new Dictionary<string, object> {
      ["ipopt"] = new Dictionary<string, object> { ["hessian_approximation"] = "exact", ["linear_solver"] = "mumps" } };
    Function solver = ca.nlpsol("solver", "ipopt", nlp, opts);

    DM w0 = to_w.call(new Dictionary<string, DM> { ["u"] = 0, ["x"] = x0 })["w"];
    DM lbw = to_w.call(new Dictionary<string, DM> { ["u"] = lbu, ["x"] = lbx })["w"];
    DM ubw = to_w.call(new Dictionary<string, DM> { ["u"] = ubu, ["x"] = ubx })["w"];
    solver.call(new Dictionary<string, DM> { ["x0"] = w0, ["lbx"] = lbw, ["ubx"] = ubw, ["lbg"] = 0, ["ubg"] = 0 });
  }

  public void test_stats_available_bug() {
    string fmu_file = data("VanDerPol2.fmu");
    assumeResource(fmu_file);
    // without WITH_GHC_FILESYSTEM casadi cannot open .fmu archives: unzip first, as its error message suggests
    if (!hasFeature("ghc-filesystem")) fmu_file = unzip(fmu_file, tempDir("vdp"));
    var dae = ca.DaeBuilder("vdp", fmu_file);
    dae.get("x1");
    dae.disp(true);
  }

  public void test_w_fmu() {
    string fmu_file = data("VanDerPol2.fmu");
    assumeResource(fmu_file);
    if (!hasFeature("ghc-filesystem")) fmu_file = unzip(fmu_file, tempDir("vdp"));
    var dae = ca.DaeBuilder("vdp", fmu_file);
    dae.get("x1");
    dae.disp(true);

    DM.rng(1);
    DM x0 = DM.rand(2);
    DM u0 = DM.rand();

    Function f = dae.create("f", new[] { "x", "u" }, new[] { "ode" }, new Dictionary<string, object> { ["aux"] = dae.w() });

    DM r1 = f.Call(x0, u0);
    var a1 = (IDictionary<string, object>)f.stats()["aux"];
    DM s1 = ca.vertcat(auxDM(a1["der(x0)"]), auxDM(a1["der(x1)"]));

    f.save("f.casadi");
    Function f2 = Function.load("f.casadi");
    DM r2 = f2.Call(x0, u0);
    var a2 = (IDictionary<string, object>)f2.stats()["aux"];
    DM s2 = ca.vertcat(auxDM(a2["der(x0)"]), auxDM(a2["der(x1)"]));

    checkarray(r1, r2);
    checkarray(s1, s2);
    File.Delete("f.casadi");
  }

  static DM auxDM(object o) {
    switch (o) {
      case double v: return v;
      case long v: return (double)v;
      case double[] v: return v;
      case IEnumerable<double> v: return v.ToArray();
      case System.Collections.IEnumerable v when !(o is DM): return v.Cast<object>().Select(Convert.ToDouble).ToArray();
      default: return (DM)o;
    }
  }

  public void test_w() {
    var dae = ca.DaeBuilder("test");

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
    DM c_val = ca.vertcat(c0_val, c1_val, c2_val);

    dae.get("c2");
    checkarray(Convert.ToDouble(dae.get("c2")), c2_val);

    MX p = dae.add("p", "parameter", "tunable");
    double p_val = 1.3;
    dae.set("p", p_val);

    MX x0 = dae.add("x0", "output", "continuous");
    double x0_val = 1.7;
    MX x1 = dae.add("x1", "local", "continuous");
    double x1_val = 1.9;

    dae.eq(dae.der(x0), x0 * x1);
    dae.eq(dae.der(x1), x1 + u + p * c0);

    MX w0 = dae.add("w0", "output", "continuous");
    MX w1 = dae.add("w1", "local", "continuous");

    dae.eq(w0, ca.sin(x1));
    dae.eq(w1, ca.pow(w0, 2.0) + u);

    MX x2 = dae.add("x2", "local", "continuous");
    double x2_val = 1.11;
    dae.eq(dae.der(x2), 3.0 * w1);

    DM x_val = ca.vertcat(x0_val, x1_val, x2_val);

    MX z0 = dae.add("z0", "output", "continuous");
    MX z1 = dae.add("z1", "local", "continuous");
    double z0_val = 2.11;
    double z1_val = 3.13;
    DM z_val = ca.vertcat(z0_val, z1_val);

    dae.eq(0, x1 - 2.0 * z0);
    dae.eq(0, x1 - 3.0 * z1 + z0);

    dae.disp(true);

    Function f = dae.create("f", new[] { "p", "u", "x" }, new[] { "y" });

    f.disp(true);

    assertInException("Initialization failed since variables [p, c0] are free",
        () => dae.dependent_fun("f", new[] { "u", "x" }, new[] { "w" }));
    f = dae.dependent_fun("f", new[] { "c", "p", "u", "x" }, new[] { "w" });
    Console.WriteLine(f);

    double s = Math.Sin(x1_val);
    DM w_val = f.Call(c_val, p_val, u_val, x_val);
    checkarray(w_val, ca.vertcat(x0_val * x1_val, x1_val + u_val + p_val * c0_val, s, s * s + u_val, 3 * (s * s + u_val)));
    f = dae.create("f", new[] { "u", "x" }, new[] { "ode" });
    double CORRECTION = 0; // Looks like a bug
    checkarray(f.Call(u_val, x_val), ca.vertcat(x0_val * x1_val, x1_val + u_val + CORRECTION * p_val * c0_val, 3 * (s * s + u_val)));

    f = dae.create("f", new[] { "p", "u", "x" }, new[] { "ode" });
    checkarray(f.Call(p_val, u_val, x_val), ca.vertcat(x0_val * x1_val, x1_val + u_val + p_val * c0_val, 3 * (s * s + u_val)));

    f = dae.create("f", new[] { "p", "u", "x", "w" }, new[] { "ode" });
    checkarray(f.Call(p_val, u_val, x_val, w_val), ca.vertcat(x0_val * x1_val, x1_val + u_val + p_val * c0_val, 3 * (s * s + u_val)));
    CORRECTION = 1; // Why is this not 2?
    checkarray(f.Call(p_val, u_val, x_val, w_val * 2.0),
        ca.vertcat(CORRECTION * x0_val * x1_val, CORRECTION * (x1_val + u_val + p_val * c0_val), 2 * 3 * (s * s + u_val)));
    f = dae.create("f", new[] { "p", "u", "x", "z" }, new[] { "alg" });
    checkarray(f.Call(p_val, u_val, x_val, z_val), ca.vertcat(-(x1_val - (2.0 * z0_val)), -((x1_val - (3 * z1_val)) + z0_val)));
    f = dae.create("f", new[] { "p", "u", "x", "z" }, new[] { "y" });
    checkarray(f.Call(p_val, u_val, x_val, z_val), ca.vertcat(c1_val, x0_val, s, z0_val));
    Dictionary<string, object> files = dae.export_fmu();
    Console.WriteLine(string.Join(", ", files.Keys));
    // python leaves the exported sources in cwd; clean them up here
    foreach (string local in files.Keys) File.Delete(local);
  }

  public void test_basemodelica() {
    Check.Assume(ca.has_modelicaparser("lacemodelica"), "modelicaparser lacemodelica not available");
    var dae_bmo = ca.DaeBuilder("from_bmo", data("DerVariable.bmo"));

    var dae_sym = ca.DaeBuilder("symbolic");
    MX x = dae_sym.add("x", "local", "continuous");
    MX y = dae_sym.add("y", "local", "continuous");
    MX dx_dt = dae_sym.add("dx_dt", "local", "continuous");
    MX dy_dt = dae_sym.add("dy_dt", "local", "continuous");
    dae_sym.eq(dae_sym.der(x), -0.5 * x + y);
    dae_sym.eq(dae_sym.der(y), x - 0.3 * y);
    dae_sym.eq(dx_dt, dae_sym.der(x));
    dae_sym.eq(dy_dt, dae_sym.der(y));

    Function f_bmo = dae_bmo.create("f_bmo", new[] { "x" }, new[] { "ode" });
    Function f_sym = dae_sym.create("f_sym", new[] { "x" }, new[] { "ode" });

    foreach (double[] x_test in new[] { new[] { 1.0, 2.0 }, new[] { 0.0, 0.0 }, new[] { 2.5, 1.5 }, new[] { -1.0, 3.0 } })
      checkfunctionLight(f_bmo, f_sym, dms(ca.DM(x_test)), 1e-10);

    // Java globs DerVariable.bmo.*.temp in the cwd, i.e. next to the .bmo file (test/data is the cwd's ../data)
    string[] temp_dirs = Directory.GetDirectories(DataDir, "DerVariable.bmo.*.temp");
    assertEquals(1, temp_dirs.Length);
    string temp_dir = temp_dirs[0];

    dae_bmo.Dispose();
    dae_sym.Dispose();
    f_bmo.Dispose();
    f_sym.Dispose();

    assertFalse(Directory.Exists(temp_dir), "Temp dir not cleaned up: " + temp_dir);
  }
}

}

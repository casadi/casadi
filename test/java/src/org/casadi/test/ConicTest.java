package org.casadi.test;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.io.File;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.UUID;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.casadi.CodeGenerator;
import org.casadi.DM;
import org.casadi.Function;
import org.casadi.GlobalOptions;
import org.casadi.MX;
import org.casadi.Opti;
import org.casadi.OptiSol;
import org.casadi.SX;
import org.casadi.Slice;
import org.casadi.Sparsity;
import org.casadi.ca;
import org.junit.jupiter.api.Assumptions;
import org.junit.jupiter.api.Test;

/** Port of test/python/conic.py */
public class ConicTest extends CasadiTestCase {

  private static final double inf = Double.POSITIVE_INFINITY;
  private static final boolean isWindows = System.getProperty("os.name").startsWith("Windows");
  private static final boolean isDarwin = System.getProperty("os.name").startsWith("Mac");

  private static final class Conic {
    final String name;
    final Map<String, Object> options;
    final Map<String, Object> aux;

    Conic(String name, Map<String, Object> options, Map<String, Object> aux) {
      this.name = name;
      this.options = options;
      this.aux = aux;
    }

    /** Fresh mutable copy of the plugin options */
    Map<String, Object> opts() {
      return new LinkedHashMap<>(options);
    }

    boolean has(String k) {
      return Boolean.TRUE.equals(aux.get(k));
    }

    int lessDigits() {
      Object v = aux.get("less_digits");
      if (v == null) return 0;
      if (v instanceof Boolean) return (Boolean) v ? 1 : 0;
      return ((Number) v).intValue();
    }

    @SuppressWarnings("unchecked")
    Map<String, Object> codegen() {
      Object v = aux.get("codegen");
      return v instanceof Map ? (Map<String, Object>) v : null;
    }
  }

  private static Map<String, Object> map(Object... kv) {
    Map<String, Object> r = new LinkedHashMap<>();
    for (int i = 0; i < kv.length; i += 2) r.put((String) kv[i], kv[i + 1]);
    return r;
  }

  private static boolean enabled(String skip, boolean available) {
    return System.getenv(skip) == null && available;
  }

  private static final List<Conic> conics = new ArrayList<>();

  static {
    if (enabled("SKIP_IPOPT_TESTS", hasNlpsol("ipopt"))) {
      Map<String, Object> ipopt_options = map("fixed_variable_treatment", "relax_bounds", "jac_c_constant", "yes", "jac_d_constant", "yes",
          "hessian_constant", "yes", "tol", 1e-12, "print_level", 0);
      conics.add(new Conic("nlpsol", map("nlpsol", "ipopt", "nlpsol_options.ipopt", ipopt_options),
          map("quadratic", true, "dual", true, "soc", false, "codegen", false, "binary", false, "discrete", false, "sos", false)));
    }
    if (enabled("SKIP_WORHP_TESTS", hasNlpsol("worhp")))
      conics.add(new Conic("nlpsol", map("nlpsol", "worhp", "nlpsol_options.worhp", map("TolOpti", 1e-13)),
          map("less_digits", 1, "quadratic", true, "dual", false, "soc", false, "codegen", false, "binary", false, "discrete", false, "sos", false)));
    if (enabled("SKIP_OOQP_TESTS", hasConic("ooqp")))
      conics.add(new Conic("ooqp", map(),
          map("less_digits", 1, "quadratic", true, "dual", true, "soc", false, "codegen", false, "binary", false, "discrete", false, "sos", false)));
    if (enabled("SKIP_QPOASES_TESTS", hasConic("qpoases")))
      conics.add(new Conic("qpoases", map("printLevel", "low"),
          map("quadratic", true, "dual", true, "soc", false, "codegen", false, "binary", false, "discrete", false, "sos", false)));
    if (enabled("SKIP_CPLEX_TESTS", hasConic("cplex")))
      conics.add(new Conic("cplex", map("cplex", map("CPX_PARAM_BARQCPEPCOMP", 1e-10, "CPX_PARAM_BAREPCOMP", 1e-10)),
          map("quadratic", true, "dual", true, "soc", true, "codegen", false, "binary", true, "discrete", true, "sos", true)));
    if (enabled("SKIP_OSQP_TESTS", hasConic("osqp"))) {
      List<String> options = isWindows ? List.of() : List.of("-Wno-unused-variable", "-Wno-strict-prototypes");
      Map<String, Object> codegen = map("extralibs", List.of("osqp"), "extra_options", options, "std", "c99");
      if (isDarwin) codegen.put("digits", 13);
      conics.add(new Conic("osqp", map("osqp", map("alpha", 1, "eps_abs", 1e-8, "eps_rel", 1e-8)),
          map("quadratic", true, "dual", true, "codegen", codegen, "soc", false, "binary", false, "discrete", false)));
    }
    if (enabled("SKIP_SUPERSCS_TESTS", hasConic("superscs")))
      conics.add(new Conic("superscs", map("superscs", map("eps", 1e-9, "do_super_scs", 1, "verbose", 0)),
          map("quadratic", true, "dual", true, "codegen", false, "soc", true, "binary", false, "discrete", false)));
    if (enabled("SKIP_GUROBI_TESTS", hasConic("gurobi")))
      conics.add(new Conic("gurobi", map("gurobi", map("BarQCPConvTol", 1e-9)),
          map("quadratic", true, "less_digits", true, "dual", true, "soc", true, "codegen", false, "binary", true, "discrete", true, "sos", true)));
    if (enabled("SKIP_XPRESS_TESTS", hasConic("xpress"))) {
      String xpress_dir = System.getenv().getOrDefault("XPRESSDIR", "missing");
      Map<String, Object> xpress_codegen = map("extralibs", List.of("xprs"), "std", "c99",
          "extralibdirs", List.of(new File(xpress_dir, "lib").getPath()), "extra_include", List.of(new File(xpress_dir, "include").getPath()));
      conics.add(new Conic("xpress", map("xpress", map("OUTPUTLOG", 0, "FEASTOL", 1e-9, "BARGAPSTOP", 1e-10, "BARDUALSTOP", 1e-10)),
          map("quadratic", true, "dual", true, "soc", true, "codegen", xpress_codegen, "binary", true, "discrete", true, "sos", true)));
    }
    if (enabled("SKIP_MOSEK_TESTS", hasConic("mosek"))) {
      String mosek_dir = System.getenv().getOrDefault("MOSEKDIR", "missing");
      Map<String, Object> mosek_codegen = map("extralibs", List.of("mosek64"), "std", "c99",
          "extralibdirs", List.of(new File(mosek_dir, "bin").getPath()), "extra_include", List.of(new File(mosek_dir, "h").getPath()));
      conics.add(new Conic("mosek", map("mosek", map("MSK_IPAR_LOG", 0, "MSK_DPAR_INTPNT_QO_TOL_REL_GAP", 1e-10, "MSK_DPAR_INTPNT_QO_TOL_PFEAS", 1e-10,
          "MSK_DPAR_INTPNT_QO_TOL_DFEAS", 1e-10, "MSK_DPAR_INTPNT_CO_TOL_REL_GAP", 1e-10, "MSK_DPAR_INTPNT_CO_TOL_PFEAS", 1e-10,
          "MSK_DPAR_INTPNT_CO_TOL_DFEAS", 1e-10)),
          map("quadratic", true, "less_digits", 1, "dual", true, "soc", true, "codegen", mosek_codegen, "binary", true, "discrete", true, "sos", false)));
    }
    if (enabled("SKIP_CLP_TESTS", hasConic("clp")))
      conics.add(new Conic("clp", map("verbose", true),
          map("quadratic", false, "dual", true, "soc", false, "codegen", false, "binary", false, "discrete", false, "sos", false)));
    if (enabled("SKIP_CBC_TESTS", hasConic("cbc")))
      conics.add(new Conic("cbc", map("verbose", true),
          map("quadratic", false, "dual", true, "soc", false, "codegen", false, "binary", true, "discrete", true, "sos", true)));
    if (enabled("SKIP_QRQP_TESTS", hasConic("qrqp"))) {
      Map<String, Object> codegen = map("std", "c99");
      if (isDarwin) codegen.put("digits", 13);
      conics.add(new Conic("qrqp", map("max_iter", 20, "print_header", false, "print_iter", false),
          map("quadratic", true, "dual", true, "soc", false, "codegen", codegen, "binary", false, "discrete", false, "sos", false)));
    }
    if (enabled("SKIP_PROXQP_TESTS", hasConic("proxqp"))) {
      for (String backend : List.of("sparse", "dense"))
        conics.add(new Conic("proxqp", map("proxqp", map("eps_abs", 1e-11, "max_iter", 1e4, "backend", backend)),
            map("quadratic", true, "dual", true, "soc", false, "codegen", false, "binary", false, "discrete", false, "sos", false)));
    }
    if (enabled("SKIP_QPALM_TESTS", hasConic("qpalm"))) {
      double eps = 1e-8;
      conics.add(new Conic("qpalm", map("qpalm", map("eps_abs", eps, "eps_rel", eps, "eps_abs_in", eps, "eps_prim_inf", eps)),
          map("quadratic", true, "dual", true, "soc", false, "codegen", false, "binary", false, "discrete", false, "sos", false)));
    }
    if (enabled("SKIP_HIGHS_TESTS", hasConic("highs"))) {
      List<String> highs_options = isWindows ? List.of() : List.of("-Wno-unused-function");
      Map<String, Object> codegen = map("extralibs", List.of("highs"), "extra_options", highs_options, "std", "c99");
      conics.add(new Conic("highs", map("highs", map("primal_feasibility_tolerance", 1e-7, "solver", "choose", "output_flag", false, "ipm_iteration_limit", 50000)),
          map("quadratic", true, "dual", true, "soc", false, "codegen", codegen, "binary", false, "discrete", true, "sos", false)));
    }
    if (enabled("SKIP_DAQP_TESTS", hasConic("daqp")))
      conics.add(new Conic("daqp", map("daqp", map()),
          map("quadratic", true, "dual", true, "soc", false, "codegen", map("extralibs", List.of("daqp"), "std", "c99"), "binary", true, "discrete", false, "sos", false)));
    if (enabled("SKIP_PIQP_TESTS", hasConic("piqp"))) {
      for (String kkt : List.of("sparse_ldlt", "dense_cholesky"))
        conics.add(new Conic("piqp", map("piqp", map("kkt_solver", kkt)),
            map("quadratic", true, "dual", false, "soc", false, "codegen", false, "binary", false, "discrete", false, "sos", false)));
    }
  }

  private static final List<String> fatrop_flags = isWindows ? List.of() : List.of("-Wno-strict-prototypes");

  private static double at(DM x, long i) {
    return x.get(i).to_double();
  }

  /** Python assertAlmostEqual(a, b, places) */
  private static void assertAlmostEqual(double a, double b, int places, String msg) {
    assertTrue(Math.abs(a - b) <= 0.5 * Math.pow(10, -places), msg + ": " + a + " != " + b + " within " + places + " places");
  }

  private static void checkarrayDigits(DM a, DM b, int digits, String msg) {
    checkarray(a, b, 0.5 * Math.pow(10, -digits), msg);
  }

  private static Map<String, DM> solverIn(DM H, DM G, DM A, DM LBX, DM UBX, DM LBA, DM UBA) {
    Map<String, DM> r = new LinkedHashMap<>();
    r.put("h", H);
    r.put("g", G);
    r.put("a", A);
    r.put("lbx", LBX);
    r.put("ubx", UBX);
    if (LBA != null) r.put("lba", LBA);
    if (UBA != null) r.put("uba", UBA);
    return r;
  }

  private static Map<String, Sparsity> hA(DM H, DM A) {
    return Map.of("h", H.sparsity(), "a", A.sparsity());
  }

  private static DM fill(int n, double v) {
    double[] r = new double[n];
    java.util.Arrays.fill(r, v);
    return new DM(r);
  }

  private static void checkSerialize(Function F, Map<String, DM> inputs) {
    String serialized = F.serialize(Map.of("debug", true));
    // validate_serialization(serialized): python-only (serialization_scheme module)
    Function F2 = Function.deserialize(serialized);
    Map<String, DM> Fout = F.call(inputs), Fout2 = F2.call(inputs);
    assertEquals(F.name_out(), F2.name_out());
    for (String k : F.name_out()) checkarray(Fout.get(k), Fout2.get(k), 1e-16, k);
  }

  /** helpers.py check_codegen: only under RUN_SLOW; returns the external Function */
  @SuppressWarnings("unchecked")
  private Function checkCodegen(Function F, Map<String, DM> inputs, Map<String, Object> cg, Map<String, Object> opts) {
    if (!RUN_SLOW) return null;
    Map<String, Object> o = cg == null ? Map.of() : cg;
    String std = (String) o.getOrDefault("std", "c89");
    int digits = ((Number) o.getOrDefault("digits", 15)).intValue();
    try {
      Path dir = Files.createTempDirectory("casadi_cg");
      String name = "codegen_" + UUID.randomUUID().toString().replace("-", "");
      CodeGenerator g = new CodeGenerator(name, opts == null ? Map.of() : opts);
      g.add(F);
      g.generate(dir.toString() + File.separator);
      String libdir = GlobalOptions.getCasadiPath();
      String includedir = GlobalOptions.getCasadiIncludePath();
      // BINDING: the Java init sets CasadiPath but leaves CasadiIncludePath empty
      if (includedir.isEmpty()) includedir = new File(libdir).getParent() + File.separator + "include";
      List<String> cmd = new ArrayList<>(List.of("gcc", "-pedantic", "-std=" + std, "-fPIC", "-shared", "-Wall", "-Werror", "-Wextra",
          "-I" + includedir, "-I" + includedir + "/highs"));
      for (String e : (List<String>) o.getOrDefault("extra_include", List.of())) cmd.add("-I" + includedir + "/" + e);
      cmd.addAll(List.of("-Wno-dangling-pointer", "-Wno-unknown-warning-option", "-Wno-unknown-pragmas", "-Wno-long-long",
          "-Wno-unused-parameter", "-O3", dir.resolve(name + ".c").toString(), "-o", dir.resolve(name + ".so").toString()));
      List<String> libdirs = new ArrayList<>(List.of(libdir));
      libdirs.addAll((List<String>) o.getOrDefault("extralibdirs", List.of()));
      for (String d : libdirs) cmd.addAll(List.of("-L" + d, "-Wl,-rpath," + d));
      for (String l : (List<String>) o.getOrDefault("extralibs", List.of())) {
        boolean found = l.startsWith("-") || l.contains(".");
        for (String d : libdirs) found |= new File(d, "lib" + l + ".so").exists() || new File(d, "lib" + l + ".a").exists();
        if (!found) {
          message("check_codegen skipped: lib" + l + " not shipped in " + libdirs);
          return null;
        }
        cmd.add(l.startsWith("-") || l.contains(".") ? l : "-l" + l);
      }
      cmd.addAll((List<String>) o.getOrDefault("extra_options", List.of()));
      Process p = new ProcessBuilder(cmd).inheritIO().start();
      assertEquals(0, p.waitFor(), "compile: " + String.join(" ", cmd));
      // check_symbols: omitted (nm-based symbol prefix check of helpers.py)
      Function F2 = ca.external(F.name(), dir.resolve(name + ".so").toString());
      Map<String, DM> r = F.call(inputs), r2 = F2.call(inputs);
      for (String k : F.name_out()) checkarrayDigits(r2.get(k), r.get(k), digits, "codegen " + k);
      return F2;
    } catch (IOException | InterruptedException e) {
      throw new RuntimeException(e);
    }
  }

  private Function checkCodegen(Function F, Map<String, DM> inputs, Map<String, Object> cg) {
    return checkCodegen(F, inputs, cg, null);
  }

  @Test
  public void test_binary() {
    DM.rng(1);
    int n = 10, m = 8, ms = 5;
    MX x = MX.sym("x", n);
    DM M = DM.rand(n, n);
    DM H = ca.mtimes(M, M.T()).times(0.5);
    DM f = DM.rand(n, 1).times(100);
    DM A = DM.rand(m, n);
    DM bupper = DM.rand(m, 1).times(20);
    DM blower = DM.rand(m, 1).times(-20);

    double[] lbx = new double[n], ubx = new double[n];
    boolean[] discrete = new boolean[n];
    for (int i = 0; i < n; i++) {
      lbx[i] = i < ms ? 0 : -10;
      ubx[i] = i < ms ? 1 : 10;
      discrete[i] = i < ms;
    }
    Map<String, DM> solver_in = Map.of("lbx", dm(lbx), "ubx", dm(ubx), "lbg", blower, "ubg", bupper);
    for (Conic c : conics) {
      if (!c.has("binary")) continue;
      if (!c.name.equals("daqp")) {
        Function solver = ca.qpsol("solver", c.name, Map.of("f", ca.mtimes(f.T(), x), "x", x, "g", ca.mtimes(A, x)), Map.of("discrete", discrete));
        solver.call(solver_in);
        if (c.codegen() != null) checkCodegen(solver, solver_in, c.codegen());
      }
      if (!c.has("quadratic")) continue;
      if (!c.name.equals("highs")) {
        MX obj = ca.mtimes(ca.mtimes(x.T(), H), x).times(0.5).plus(ca.mtimes(f.T(), x));
        Function solver = ca.qpsol("solver", c.name, Map.of("f", obj, "x", x, "g", ca.mtimes(A, x)), Map.of("discrete", discrete));
        solver.call(solver_in);
        int digits = c.name.equals("mosek") ? 4 : 15;
        if (c.codegen() != null) {
          Map<String, Object> cg = new LinkedHashMap<>(c.codegen());
          cg.put("digits", digits);
          checkCodegen(solver, solver_in, cg);
        }
      }
    }
  }

  @Test
  public void test_mosek_codegen_option_precision() {
    Assumptions.assumeTrue(hasConic("mosek"));
    Map<String, Object> options = map("MSK_DPAR_MIO_TOL_ABS_RELAX_INT", 1.234567890123456e-7, "MSK_DPAR_OPTIMIZER_MAX_TIME", 1234567890.1234567);
    Function solver = ca.conic("solver", "mosek", Map.of("h", Sparsity.dense(1, 1)), Map.of("mosek", options));
    CodeGenerator cg = new CodeGenerator("mosek_precision.c");
    cg.add(solver);
    Map<String, String> emitted = new HashMap<>();
    Matcher mt = Pattern.compile("MSK_putnadouparam\\(d->task, \"([^\"]+)\", ([^)]+)\\)").matcher(cg.dump());
    while (mt.find()) emitted.put(mt.group(1), mt.group(2));
    assertEquals(options.keySet(), emitted.keySet());
    for (Map.Entry<String, Object> e : options.entrySet())
      assertEquals((Double) e.getValue(), Double.parseDouble(emitted.get(e.getKey())), e.getKey());
  }

  @Test
  public void test_opti() {
    for (Conic c : conics) {
      if (c.name.equals("qrqp")) continue;
      if (c.name.equals("daqp")) continue;
      Opti opti = new Opti("conic");
      MX x = opti.variable();
      opti.minimize(x);
      opti.subject_to(ca.le(-2.0, ca.le(x, 2.0)));
      opti.solver(c.name, c.opts());
      OptiSol sol = opti.solve();
      checkarrayDigits(sol.value(x), new DM(-2), 5, c.name);
    }
  }

  @Test
  public void test_sos() {
    DM H = new DM(4, 4);
    DM G = dm(-1, -2, -3, -1);
    DM A = dm(new double[][] {{-1, 1, 1, 10}, {1, -3, 1, 0}, {0, 1, 0, -3.5}});
    DM LBA = dm(-inf, -inf, 0);
    DM UBA = dm(20, 30, 0);
    boolean[] discrete = {false, true, true, true};
    DM LBX = dm(0, 0, 0, 2);
    DM UBX = dm(40, inf, inf, 3);
    List<long[]> sos_groups = List.of(new long[] {2, 3});
    List<double[]> sos_weights = List.of(new double[] {25.0, 18.0});

    for (Conic c : conics) {
      if (!c.has("discrete")) continue;
      System.out.println("test_sos " + c.name + " " + c.options);
      Map<String, Object> options = c.opts();
      options.put("discrete", discrete);
      Function solver = ca.conic("mysolver", c.name, hA(H, A), options);
      Map<String, DM> solver_in = solverIn(H, G, A, LBX, UBX, LBA, UBA);
      Map<String, DM> solver_out = solver.call(solver_in);
      checkarray(solver_out.get("x"), dm(29, 7, 22, 2));

      if (!c.has("sos")) continue;
      options.put("sos_groups", sos_groups);
      solver = ca.conic("mysolver", c.name, hA(H, A), options);
      solver_out = solver.call(solver_in);
      checkarray(solver_out.get("x"), dm(40, 7, 0, 2));

      options.put("sos_weights", sos_weights);
      solver = ca.conic("mysolver", c.name, hA(H, A), options);
      solver_out = solver.call(solver_in);
      checkarray(solver_out.get("x"), dm(40, 7, 0, 2));
      System.out.println(solver_out.get("cost"));
      checkSerialize(solver, solver_in);
    }
  }

  @Test
  public void test_milp() {
    DM H = new DM(2, 2);
    DM G = dm(-1, -1);
    DM A = dm(new double[][] {{-2, 2}, {-8, 10}});
    DM LBA = dm(1, -inf);
    DM UBA = dm(inf, 13);
    DM LBX = fill(2, 0);
    DM UBX = fill(2, inf);
    boolean[] discrete = {true, true};

    for (Conic c : conics) {
      if (!c.has("discrete")) continue;
      System.out.println("test_milp " + c.name + " " + c.options);
      Map<String, Object> options = c.opts();
      options.put("discrete", discrete);
      Function solver = ca.conic("mysolver", c.name, hA(H, A), options);
      Map<String, DM> solver_in = solverIn(H, G, A, LBX, UBX, LBA, UBA);
      Map<String, DM> solver_out = solver.call(solver_in);
      checkarray(solver_out.get("x"), dm(1, 2));
      checkarray(solver_out.get("cost"), dm(-3));
      checkSerialize(solver, solver_in);
    }

    H = new DM(3, 3);
    G = dm(-1, -1, -5);
    A = dm(new double[][] {{-2, 2, 0}, {-8, 10, 0}, {0, 0, 7}});
    LBA = dm(1, -inf, -4.0 * 7);
    UBA = dm(inf, 13, 4.0 * 7);
    LBX = dm(0, 0, -5);
    UBX = dm(inf, inf, 5);
    discrete = new boolean[] {true, true, false};

    for (Conic c : conics) {
      if (!c.has("discrete")) continue;
      System.out.println("test_milp " + c.name + " " + c.options);
      Map<String, Object> options = c.opts();
      options.put("discrete", discrete);
      Function solver = ca.conic("mysolver", c.name, hA(H, A), options);
      Map<String, DM> solver_out = solver.call(solverIn(H, G, A, LBX, UBX, LBA, UBA));
      checkarray(solver_out.get("x"), dm(1, 2, 4));
      checkarray(solver_out.get("cost"), dm(-23));
    }
  }

  @Test
  public void test_general_unconstrained() {
    DM H = ca.sparsify(dm(new double[][] {{1, 0}, {0, 1}}));
    DM G = dm(-0.7, -2.3);
    DM A = DM.zeros(0, 2);
    DM LBA = DM.zeros(0, 1);
    DM UBA = DM.zeros(0, 1);
    DM LBX = fill(2, -inf);
    DM UBX = fill(2, inf);

    for (Conic c : conics) {
      if (!c.has("quadratic")) continue;
      System.out.println("test_general_unconstrained " + c.name + " " + c.options);
      Function solver = ca.conic("mysolver", c.name, hA(H, A), c.opts());
      int less_digits = c.lessDigits();
      Map<String, DM> solver_in = solverIn(H, G, A, LBX, UBX, LBA, UBA);
      Map<String, DM> solver_out = solver.call(solver_in);
      assertTrue((Boolean) solver.stats().get("success"), c.name);
      assertAlmostEqual(at(solver_out.get("x"), 0), 0.7, Math.max(1, 6 - less_digits), c.name);
      assertAlmostEqual(at(solver_out.get("x"), 1), 2.3, Math.max(1, 6 - less_digits), c.name);
      if (c.codegen() != null) checkCodegen(solver, solver_in, c.codegen());
      checkSerialize(solver, solver_in);
    }
  }

  @Test
  public void test_general_convex_dense() {
    StringBuilder names = new StringBuilder();
    for (Conic c : conics) names.append(c.name).append(" ");
    message("Convex dense QP with solvers: " + names);
    DM H = dm(new double[][] {{1, -1}, {-1, 2}});
    DM G = dm(-2, -6);
    DM A = dm(new double[][] {{1, 1}, {-1, 2}, {2, 1}});
    DM LBA = fill(3, -inf);
    DM UBA = dm(2, 2, 3);
    DM LBX = fill(2, 0);
    DM UBX = fill(2, inf);

    for (Conic c : conics) {
      if (c.name.equals("nlpsol")) continue;
      if (!c.has("quadratic")) continue;
      System.out.println("test_general_convex_dense " + c.name + " " + c.options);
      if (c.name.equals("daqp")) continue;
      Function solver = ca.conic("mysolver", c.name, hA(H, A), c.opts());
      try {
        solver.stats();
      } catch (RuntimeException e) {
        // nosegfault
      }
      int ld = c.lessDigits();
      String s = c.name;
      Map<String, DM> solver_in = solverIn(H, G, A, LBX, UBX, LBA, UBA);
      Map<String, DM> solver_out = solver.call(solver_in);
      assertTrue((Boolean) solver.stats().get("success"));
      assertEquals("SOLVER_RET_SUCCESS", solver.stats().get("unified_return_status"));
      assertAlmostEqual(at(solver_out.get("x"), 0), 2.0 / 3, Math.max(1, 6 - ld), s);
      assertAlmostEqual(at(solver_out.get("x"), 1), 4.0 / 3, Math.max(1, 6 - ld), s);
      if (c.has("dual")) assertAlmostEqual(at(solver_out.get("lam_x"), 0), 0, Math.max(1, 6 - ld), s);
      if (c.has("dual")) assertAlmostEqual(at(solver_out.get("lam_x"), 1), 0, Math.max(1, 6 - ld), s);
      if (c.has("dual")) checkarrayDigits(solver_out.get("lam_a"), dm(3 + 1.0 / 9, 4.0 / 9, 0), Math.max(1, 6 - ld), s);
      assertAlmostEqual(at(solver_out.get("cost"), 0), -8 - 2.0 / 9, Math.max(1, 6 - ld), s);

      solver_in.put("h", H.times(4));
      solver_out = solver.call(solver_in);
      assertAlmostEqual(at(solver_out.get("x"), 0), 1, Math.max(1, 3 - ld), s);
      assertAlmostEqual(at(solver_out.get("x"), 1), 1, Math.max(1, 3 - ld), s);
      assertAlmostEqual(solver_out.get("cost").to_double(), -6, Math.max(1, 6 - ld), s);
      if (c.has("dual")) assertAlmostEqual(at(solver_out.get("lam_x"), 0), 0, Math.max(1, 6 - ld), s);
      if (c.has("dual")) assertAlmostEqual(at(solver_out.get("lam_x"), 1), 0, Math.max(1, 6 - ld), s);
      if (c.has("dual")) checkarrayDigits(solver_out.get("lam_a"), dm(2, 0, 0), Math.max(1, 2 - ld), s);

      solver_in.put("h", new DM(0));
      if (c.name.contains("qcqp")) continue;
      solver_out = solver.call(solver_in);
      assertAlmostEqual(at(solver_out.get("x"), 0), 2.0 / 3, Math.max(1, 6 - ld), s);
      assertAlmostEqual(at(solver_out.get("x"), 1), 4.0 / 3, Math.max(1, 6 - ld), s);
      assertAlmostEqual(solver_out.get("cost").to_double(), -9 - 1.0 / 3, Math.max(1, 6 - ld), s);
      if (c.has("dual")) assertAlmostEqual(at(solver_out.get("lam_x"), 0), 0, Math.max(1, 6 - ld), s);
      if (c.has("dual")) assertAlmostEqual(at(solver_out.get("lam_x"), 1), 0, Math.max(1, 6 - ld), s);
      if (c.has("dual")) checkarrayDigits(solver_out.get("lam_a"), dm(10.0 / 3, 4.0 / 3, 0), Math.max(1, 4 - ld), s);

      solver_in.put("lba", fill(3, -inf));
      solver_in.put("uba", fill(3, inf));
      solver_in.put("ubx", new DM(5));
      if (c.options.toString().contains("worhp")) {
        assertThrows(RuntimeException.class, () -> solver.invoke(new Object[] {solver_in}));
        continue;
      }
      solver_out = solver.call(solver_in);
      assertAlmostEqual(at(solver_out.get("x"), 0), 5, Math.max(1, 6 - ld), s);
      assertAlmostEqual(at(solver_out.get("x"), 1), 5, Math.max(1, 6 - ld), s);
      assertAlmostEqual(solver_out.get("cost").to_double(), -40, Math.max(1, 5 - ld), s);
      if (c.has("dual")) assertAlmostEqual(at(solver_out.get("lam_x"), 0), 2, Math.max(1, 6 - ld), s);
      if (c.has("dual")) assertAlmostEqual(at(solver_out.get("lam_x"), 1), 6, Math.max(1, 6 - ld), s);
      if (c.has("dual")) checkarrayDigits(solver_out.get("lam_a"), dm(0, 0, 0), Math.max(1, 4 - ld), s);
    }
  }

  @Test
  public void test_general_convex_sparse() {
    StringBuilder names = new StringBuilder();
    for (Conic c : conics) names.append(c.name).append(" ");
    message("Convex sparse QP with solvers: " + names);
    DM H = ca.diag(dm(2, 1, 0.2, 0.7, 1.3));
    H.set(1, 2, 0.1);
    H.set(2, 1, 0.1);
    DM G = dm(-2, -6, 1, 0, 0);
    DM A = ca.sparsify(dm(new double[][] {{1, 0, 0.1, 0.7, -1}, {0.1, 2, -0.3, 4, 0.1}}));
    DM LBA = dm(-inf);
    DM UBA = dm(2, 2);
    DM LBX = fill(5, 0);
    DM UBX = fill(5, inf);

    for (Conic c : conics) {
      if (!c.has("quadratic")) continue;
      System.out.println("test_general_convex_sparse " + c.name + " " + c.options);
      Function solver = ca.conic("mysolver", c.name, hA(H, A), c.opts());
      int ld = c.lessDigits();
      Map<String, DM> solver_in = solverIn(H, G, A, LBX, UBX, LBA, UBA);
      Map<String, DM> solver_out = solver.call(solver_in);
      checkarrayDigits(solver_out.get("x"), dm(0.873908, 0.95630465, 0, 0, 0), Math.max(1, 6 - ld), c.name);
      if (c.has("dual")) checkarrayDigits(solver_out.get("lam_x"), dm(0, 0, -0.339076, -10.0873907, -0.252185), Math.max(1, 6 - ld), c.name);
      if (c.has("dual")) checkarrayDigits(solver_out.get("lam_a"), dm(0, 2.52184767), Math.max(1, 6 - ld), c.name);
      assertAlmostEqual(at(solver_out.get("cost"), 0), -6.264669320767, Math.max(1, 6 - ld), c.name);
      if (c.codegen() != null) checkCodegen(solver, solver_in, c.codegen());
    }
  }

  @Test
  public void test_general_nonconvex_dense() {
    StringBuilder names = new StringBuilder();
    for (Conic c : conics) names.append(c.name).append(" ");
    message("Non convex dense QP with solvers: " + names);
    DM H = dm(new double[][] {{1, -1}, {-1, -2}});
    DM G = dm(-2, -6);
    DM A = dm(new double[][] {{1, 1}, {-1, 2}, {2, 1}});
    DM LBA = fill(3, -inf);
    DM UBA = dm(2, 2, 3);
    DM LBX = fill(2, 0);
    DM UBX = fill(2, inf);

    for (Conic c : conics) {
      if (!c.has("quadratic")) continue;
      System.out.println("test_general_nonconvex_dense " + c.name + " " + c.options);
      Function solver = ca.conic("mysolver", c.name, hA(H, A), c.opts());
      Map<String, DM> solver_in = solverIn(H, G, A, LBX, UBX, LBA, UBA);
      // Python passes the dict positionally (solver(solver_in)), which fails input conversion for every plugin
      assertThrows(RuntimeException.class, () -> solver.invoke(new Object[] {solver_in}), c.name);
    }
  }

  @Test
  public void test_equality() {
    message("Regression 452 test: equality constraints give wrong multipliers");
    DM H = dm(new double[][] {{1, -1}, {-1, 2}});
    DM G = dm(-2, -6);

    for (Conic c : conics) {
      if (!c.has("quadratic")) continue;
      if (c.name.contains("ooqp") || c.name.contains("proxqp")) continue;
      System.out.println("test_equality " + c.name + " " + c.options);
      Function solver = ca.conic("mysolver", c.name, Map.of("h", H.sparsity(), "a", Sparsity.dense(3, 2)), c.opts());
      int ld = c.lessDigits();
      String s = c.name;
      DM A = dm(new double[][] {{1, 1}, {-1, 2}, {2, 1}});
      DM LBA = fill(3, -inf);
      DM UBA = dm(2, 2, 3);
      DM LBX = dm(0.5, 0);
      DM UBX = dm(0.5, inf);
      Map<String, DM> solver_in = solverIn(H, G, A, LBX, UBX, LBA, UBA);
      if (c.options.toString().contains("worhp")) {
        assertThrows(RuntimeException.class, () -> solver.invoke(new Object[] {solver_in}));
        continue;
      }
      Map<String, DM> solver_out = solver.call(solver_in);
      assertAlmostEqual(at(solver_out.get("x"), 0), 0.5, Math.max(1, 6 - ld), s);
      assertAlmostEqual(at(solver_out.get("x"), 1), 1.25, Math.max(1, 6 - ld), s);
      if (c.has("dual")) assertAlmostEqual(at(solver_out.get("lam_x"), 0), 4.75, Math.max(1, 6 - ld), s);
      if (c.has("dual")) assertAlmostEqual(at(solver_out.get("lam_x"), 1), 0, Math.max(1, 6 - ld), s);
      if (c.has("dual")) checkarrayDigits(solver_out.get("lam_a"), dm(0, 2, 0), Math.max(1, 6 - ld), s);
      if (c.codegen() != null) checkCodegen(solver, solver_in, c.codegen());
      assertAlmostEqual(at(solver_out.get("cost"), 0), -7.4375, Math.max(1, 6 - ld), s);

      solver_in.put("h", H);
      solver_in.put("g", G);
      solver_in.put("a", dm(new double[][] {{1, 1}, {-1, 2}, {2, 1}}));
      solver_in.put("lbx", fill(2, -inf));
      solver_in.put("ubx", fill(2, inf));
      solver_in.put("lba", dm(2, -inf, -inf));
      solver_in.put("uba", dm(2, inf, inf));
      solver_out = solver.call(solver_in);
      assertAlmostEqual(at(solver_out.get("x"), 0), 0.4, Math.max(1, 4 - ld), s);
      assertAlmostEqual(at(solver_out.get("x"), 1), 1.6, Math.max(1, 4 - ld), s);
      if (c.has("dual")) assertAlmostEqual(at(solver_out.get("lam_x"), 0), 0, Math.max(1, 5 - ld), s);
      if (c.has("dual")) assertAlmostEqual(at(solver_out.get("lam_x"), 1), 0, Math.max(1, 5 - ld), s);
      if (c.has("dual")) checkarrayDigits(solver_out.get("lam_a"), dm(3.2, 0, 0), Math.max(1, 5 - ld), s);
      assertAlmostEqual(at(solver_out.get("cost"), 0), -8.4, Math.max(1, 5 - ld), s);
      if (c.codegen() != null) checkCodegen(solver, solver_in, c.codegen());
    }
  }

  @Test
  public void test_degenerate_hessian() {
    message("Degenerate hessian");
    DM H = ca.sparsify(dm(new double[][] {{1, -1, 0}, {-1, 2, 0}, {0, 0, 0}}));
    DM G = dm(-2, -6, 1);
    DM A = dm(new double[][] {{1, 1, 1}});
    DM LBA = dm(0.5);
    DM UBA = dm(0.5);
    DM LBX = dm(-10);
    DM UBX = dm(10);

    for (Conic c : conics) {
      if (c.name.contains("qrqp")) continue;
      if (!c.has("quadratic")) continue;
      // Python mutates the shared table here (dump_in leaks into later tests); kept local
      Map<String, Object> qp_options = c.opts();
      qp_options.put("dump_in", true);
      if (c.name.contains("qcqp")) continue;
      if (c.name.equals("daqp")) continue;
      System.out.println("test_degenerate_hessian " + c.name + " " + qp_options);
      Function solver = ca.conic("mysolver", c.name, hA(H, A), qp_options);
      int ld = c.lessDigits();
      Map<String, DM> solver_in = solverIn(H, G, A, LBX, UBX, LBA, UBA);
      Map<String, DM> solver_out = solver.call(solver_in);
      checkarrayDigits(solver_out.get("x"), dm(5.5, 5, -10), Math.max(1, 4 - ld), c.name);
      if (c.has("dual")) checkarrayDigits(solver_out.get("lam_x"), dm(0, 0, -2.5), Math.max(1, 4 - ld), c.name);
      if (c.has("dual")) checkarrayDigits(solver_out.get("lam_a"), dm(1.5), Math.max(1, 4 - ld), c.name);
      assertAlmostEqual(at(solver_out.get("cost"), 0), -38.375, Math.max(1, 5 - ld), c.name);
      if (c.codegen() != null) checkCodegen(solver, solver_in, c.codegen());
    }
  }

  @Test
  public void test_no_inequality() {
    message("No inequalities present");
    DM H = dm(new double[][] {{1, -1}, {-1, 2}});
    DM G = dm(-2, -6);
    DM A = dm(new double[][] {{1, 1}});
    DM LBA = dm(0.5);
    DM UBA = dm(0.5);
    DM LBX = dm(-10);
    DM UBX = dm(10);

    for (Conic c : conics) {
      if (!c.has("quadratic")) continue;
      System.out.println("test_no_inequality " + c.name + " " + c.options);
      Function solver = ca.conic("mysolver", c.name, hA(H, A), c.opts());
      int ld = c.lessDigits();
      String s = c.name;
      Map<String, DM> solver_out = solver.call(solverIn(H, G, A, LBX, UBX, LBA, UBA));
      assertAlmostEqual(at(solver_out.get("x"), 0), -0.5, Math.max(1, 6 - ld), s);
      assertAlmostEqual(at(solver_out.get("x"), 1), 1, Math.max(1, 6 - ld), s);
      if (c.has("dual")) assertAlmostEqual(at(solver_out.get("lam_x"), 0), 0, Math.max(1, 6 - ld), s);
      if (c.has("dual")) assertAlmostEqual(at(solver_out.get("lam_x"), 1), 0, Math.max(1, 6 - ld), s);
      if (c.has("dual")) checkarrayDigits(solver_out.get("lam_a"), dm(3.5), Math.max(1, 6 - ld), s);
      assertAlmostEqual(at(solver_out.get("cost"), 0), -3.375, Math.max(1, 6 - ld), s);
    }
  }

  @Test
  public void test_no_A() {
    message("No A present");
    DM H = dm(new double[][] {{1, -1}, {-1, 2}});
    DM G = dm(-2, -6);
    DM A = new DM(0, 2);
    DM LBA = new DM(0, 1);
    DM UBA = new DM(0, 1);
    DM LBX = dm(-10);
    DM UBX = dm(10);

    for (Conic c : conics) {
      if (!c.has("quadratic")) continue;
      System.out.println("test_no_A " + c.name + " " + c.options);
      Function solver = ca.conic("mysolver", c.name, hA(H, A), c.opts());
      int ld = c.lessDigits();
      Map<String, DM> solver_in = solverIn(H, G, A, LBX, UBX, LBA, UBA);
      Map<String, DM> solver_out = solver.call(solver_in);
      checkarrayDigits(solver_out.get("x"), dm(10, 8), Math.max(1, 3 - ld), c.name);
      if (c.has("dual")) checkarrayDigits(solver_out.get("lam_x"), dm(0, 0), Math.max(1, 4 - ld), c.name);
      if (c.has("dual")) checkarrayDigits(solver_out.get("lam_a"), dm(), Math.max(1, 5 - ld), c.name);
      assertAlmostEqual(at(solver_out.get("cost"), 0), -34, Math.max(1, 5 - ld), c.name);
      if (c.codegen() != null) checkCodegen(solver, solver_in, c.codegen());
    }
  }

  @Test
  public void test_standard_form() {
    DM H = dm(new double[][] {{1, -1}, {-1, 2}});
    DM G = dm(-2, -6);
    DM A = dm(1, 1).T();
    DM LBA = dm(-inf);
    DM UBA = dm(1);
    DM LBX = dm(-10);
    DM UBX = dm(10);

    for (Conic c : conics) {
      if (!c.has("quadratic")) continue;
      System.out.println("test_standard_form " + c.name + " " + c.options);
      Function solver = ca.conic("mysolver", c.name, hA(H, A), c.opts());
      int ld = c.lessDigits();
      Map<String, DM> solver_in = solverIn(H, G, A, LBX, UBX, LBA, UBA);
      Map<String, DM> solver_out = solver.call(solver_in);
      checkarrayDigits(solver_out.get("x"), dm(-0.2, 1.2), Math.max(1, 3 - ld), c.name);
      if (c.has("dual")) checkarrayDigits(solver_out.get("lam_x"), dm(0, 0), Math.max(1, 4 - ld), c.name);
      if (c.has("dual")) checkarrayDigits(solver_out.get("lam_a"), dm(3.4), Math.max(1, 5 - ld), c.name);
      assertAlmostEqual(at(solver_out.get("cost"), 0), -5.1, Math.max(1, 5 - ld), c.name);
      if (c.codegen() != null) checkCodegen(solver, solver_in, c.codegen());
    }
  }

  @Test
  public void test_badscaling() {
    message("Badly scaled problem");
    int N = 50;
    double[] d = new double[N], r = new double[N];
    for (int i = 0; i < N; i++) {
      d[i] = i + 1;
      r[i] = i;
    }
    DM H = ca.diag(dm(d));
    DM x0 = dm(r);
    DM G = ca.mtimes(H, x0).times(-1.0);
    DM A = new DM(0, N);
    DM LBX = fill(N, -1000);
    DM UBX = fill(N, 1000);

    for (Conic c : conics) {
      System.out.println(c.name + " " + c.options);
      if (!c.has("quadratic")) continue;
      if (c.name.contains("worhp")) continue;
      if (c.name.contains("superscs")) continue;
      System.out.println("test_badscaling " + c.name + " " + c.options);
      Function solver = ca.conic("mysolver", c.name, hA(H, A), c.opts());
      int ld = c.lessDigits();
      Map<String, DM> solver_out = solver.call(solverIn(H, G, A, LBX, UBX, null, null));
      checkarrayDigits(solver_out.get("x"), x0, Math.max(1, 2 - ld), c.name + c.options);
      assertAlmostEqual(at(solver_out.get("cost"), 0), -0.5 * ca.mtimes(List.of(x0.T(), H, x0)).to_double(), Math.max(1, 3 - ld), c.name);
      if (c.has("dual")) checkarrayDigits(solver_out.get("lam_x"), DM.zeros(N, 1), Math.max(1, 4 - ld), c.name);
    }
  }

  @Test
  public void test_redundant() {
    message("Redundant constraints");
    DM H = dm(new double[][] {{1, -1, 0}, {-1, 2, 0}, {0, 0, 0}});
    DM G = dm(-2, -6, 1);
    DM a = dm(1, 0, 1);
    DM a_ = dm(0, 1, -2);

    for (double[] w : new double[][] {{0, 2}, {1, 1}, {0.1, 0.6}}) {
      double w0 = w[0], w1 = w[1];
      DM A = ca.vertcat(List.of(a.T(), a_.T(), a.times(w0).plus(a_.times(w1)).T()));
      DM LBA = dm(0, 0, 0);
      DM UBA = dm(0.5, 0.3, w0 * 0.5 + w1 * 0.3);
      DM LBX = dm(-10);
      DM UBX = dm(10);

      for (Conic c : conics) {
        if (!c.has("quadratic")) continue;
        if (c.name.contains("qcqp")) continue;
        System.out.println("test_redundant " + c.name + " " + c.options);
        if (c.name.equals("daqp")) continue;
        Function solver = ca.conic("qpsol", c.name, hA(H, A), c.opts());
        int ld = c.lessDigits();
        Map<String, DM> solver_out = solver.call(solverIn(H, G, A, LBX, UBX, LBA, UBA));
        checkarrayDigits(solver_out.get("x"), dm(-0.19230768069, 1.6846153915, 0.692307690769276), Math.max(1, 6 - ld), c.name);
        assertAlmostEqual(at(solver_out.get("cost"), 0), -5.850384678537, Math.max(1, 5 - ld), c.name);
        if (c.has("dual")) checkarrayDigits(solver_out.get("lam_x"), dm(0, 0, 0), Math.max(1, 6 - ld), c.name);
        if (c.has("dual"))
          checkarrayDigits(ca.mtimes(A.T(), solver_out.get("lam_a")), dm(3.876923073076, 2.4384615365384965, -1), Math.max(1, 6 - ld), c.name);
      }
    }
  }

  @Test
  public void test_linear() {
    DM H = new DM(2, 2);
    DM A = dm(new double[][] {{-1, 1}, {1, 1}, {1, -2}});
    DM LBA = dm(-inf, 2, -inf);
    DM UBA = dm(1, inf, 4);
    DM LBX = dm(-inf, 0);
    DM UBX = dm(inf, inf);
    DM G = dm(2, 1);

    for (Conic c : conics) {
      if (c.name.contains("qcqp")) continue;
      if (c.name.equals("daqp")) continue;
      System.out.println("test_linear " + c.name + " " + c.options);
      Function solver = ca.conic("mysolver", c.name, hA(H, A), c.opts());
      int ld = c.lessDigits();
      Map<String, DM> solver_in = solverIn(H, G, A, LBX, UBX, LBA, UBA);
      if (!(c.name.equals("superscs") && isDarwin)) checkSerialize(solver, solver_in);
      Map<String, DM> solver_out = solver.call(solver_in);
      checkarrayDigits(solver_out.get("x"), dm(0.5, 1.5), Math.max(1, 5 - ld), c.name);
      if (c.has("dual")) checkarrayDigits(solver_out.get("lam_x"), dm(0, 0), Math.max(1, 5 - ld), c.name);
      if (c.has("dual")) checkarrayDigits(solver_out.get("lam_a"), dm(0.5, -1.5, 0), Math.max(1, 5 - ld), c.name);
      assertAlmostEqual(at(solver_out.get("cost"), 0), 2.5, Math.max(1, 5 - ld), c.name);
      if (c.codegen() != null) checkCodegen(solver, solver_in, c.codegen());
    }
  }

  @Test
  public void test_linear2() {
    DM H = new DM(2, 2);
    DM A = dm(new double[][] {{-1, 1}, {1, 1}, {1, -2}});
    DM LBA = dm(-inf, 2, -inf);
    DM UBA = dm(1, inf, 4);
    DM LBX = dm(-inf, 3);
    DM UBX = dm(inf, 3);
    DM G = dm(2.0, 1.0);

    for (Conic c : conics) {
      if (c.name.contains("qcqp")) continue;
      if (c.name.contains("nlp")) continue;
      if (c.name.equals("daqp")) continue;
      System.out.println("test_linear2 " + c.name + " " + c.options);
      Function solver = ca.conic("msyolver", c.name, hA(H, A), c.opts());
      int ld = c.lessDigits();
      Map<String, DM> solver_out = solver.call(solverIn(H, G, A, LBX, UBX, LBA, UBA));
      checkarrayDigits(solver_out.get("x"), dm(2, 3), Math.max(1, 5 - ld), c.name);
      if (c.has("dual")) checkarrayDigits(solver_out.get("lam_x"), dm(0, -3), Math.max(1, 5 - ld), c.name);
      if (c.has("dual")) checkarrayDigits(solver_out.get("lam_a"), dm(2, 0, 0), Math.max(1, 5 - ld), c.name);
      assertAlmostEqual(at(solver_out.get("cost"), 0), 7, Math.max(1, 5 - ld), c.name);
    }
  }

  @Test
  public void test_overconstrained() {
    SX x = SX.sym("x");
    Map<String, SX> qp = Map.of("x", x, "f", ca.power(x.minus(1), 2), "g", ca.vertcat(List.of(x, x, x)));

    for (Conic c : conics) {
      if (!c.has("quadratic")) continue;
      System.out.println("test_overconstrained " + c.name + " " + c.options);
      Function solver = ca.qpsol("mysolver", c.name, qp, c.opts());
      Map<String, DM> solver_in = Map.of("lbx", dm(-10), "ubx", dm(10), "lbg", dm(-10, -10, -10), "ubg", dm(10, 10, 10), "x0", dm(0));
      Map<String, DM> solver_out = solver.call(solver_in);
      assertAlmostEqual(at(solver_out.get("f"), 0), 0, 9, c.name);
      assertAlmostEqual(at(solver_out.get("x"), 0), 1, 5, c.name);
      if (c.codegen() != null) checkCodegen(solver, solver_in, c.codegen());
    }
  }

  @Test
  public void test_qrqp_prints() {
    Assumptions.assumeTrue(hasConic("qrqp"));
    MX x = MX.sym("x");
    Function solver = ca.qpsol("solver", "qrqp", Map.of("x", x, "f", ca.power(x.minus(1), 2)));
    String out = captureStdout(() -> solver.call(new HashMap<String, DM>()));
    assertTrue(out.contains("last_tau"), out);
    assertTrue(out.contains("Converged"), out);
    if (RUN_SLOW) {
      Function F = checkCodegen(solver, new HashMap<>(), map("std", "c99"), map());
      String out2 = captureStdout(() -> F.call(new HashMap<String, DM>()));
      assertFalse(out2.contains("last_tau"), out2);
      assertFalse(out2.contains("Converged"), out2);
      checkCodegen(solver, new HashMap<>(), map("std", "c99"), map("verbose_runtime", true));
    }
  }

  @Test
  public void test_hpipm() {
    Assumptions.assumeTrue(hasConic("hpipm"));
    Assumptions.assumeTrue(hasConic("qpoases"));
    double inf = 100;
    int N = 4;
    MX x1 = MX.sym("x1");
    MX x2 = MX.sym("x2");
    MX x = ca.vertcat(List.of(x1, x2));
    MX u = MX.sym("u");
    MX xdot = ca.vertcat(List.of(x1.times(0.6).minus(x2.times(1.11)).plus(u.times(0.3)).minus(0.03), x1.times(0.7).plus(0.01)));
    MX L = ca.power(x1, 2).plus(ca.power(x2, 2).times(3)).plus(ca.power(u, 2).times(7)).minus(x1.times(x2).times(0.4))
        .minus(x1.times(u).times(0.3)).plus(u).minus(x1).minus(x2.times(2));
    Function F = new Function("F", List.of(x, u), List.of(x.plus(xdot), L));
    F.jacobian_old(0, 0);

    List<SX> w = new ArrayList<>();
    List<Double> w0 = new ArrayList<>(), lbw = new ArrayList<>(), ubw = new ArrayList<>(), lbg = new ArrayList<>(), ubg = new ArrayList<>();
    SX J = new SX(0);
    List<SX> g = new ArrayList<>();
    List<SX> Xs = SX.sym("X", 2, 1, N + 1);
    List<SX> Us = SX.sym("U", 1, 1, N + 1);
    for (int k = 0; k < N; k++) {
      w.add(Xs.get(k));
      if (k == 0) {
        lbw.addAll(List.of(-inf, 1.0));
        ubw.addAll(List.of(inf, 1.0));
        w0.addAll(List.of(0.3, 0.7));
      } else if (k == 2) {
        lbw.addAll(List.of(0.0, -inf));
        ubw.addAll(List.of(0.0, inf));
        w0.addAll(List.of(0.0, 0.3));
      } else {
        lbw.addAll(List.of(-inf, -inf));
        ubw.addAll(List.of(inf, inf));
        w0.addAll(List.of(0.0, 1.0));
      }
      w.add(Us.get(k));
      lbw.add(-inf);
      ubw.add(inf);
      w0.add(0.0);
      List<SX> r = F.call(List.of(Xs.get(k), Us.get(k)));
      SX xplus = r.get(0), l = r.get(1);
      J = J.plus(l);
      g.add(xplus.minus(Xs.get(k + 1)).times(3));
      lbg.addAll(List.of(0.0, 0.0));
      ubg.addAll(List.of(0.0, 0.0));
      g.add(Xs.get(k).get(1).times(0.1).minus(Us.get(k).times(0.05)));
      lbg.add(-0.5 * k - 0.1);
      ubg.add(2.0);
    }
    SX Xl = Xs.get(Xs.size() - 1);
    g.add(Xl.get(1).times(0.1));
    lbg.add(0.1);
    ubg.add(2.0);
    J = J.plus(ca.mtimes(Xl.T(), Xl));
    w.add(Xl);
    lbw.addAll(List.of(-inf, -inf));
    ubw.addAll(List.of(inf, inf));
    w0.addAll(List.of(0.0, 1.0));

    Map<String, SX> prob = Map.of("f", J, "x", ca.vertcat(w), "g", ca.vertcat(g));
    Function Jf = new Function("J", List.of(prob.get("x")), List.of(ca.jacobian(prob.get("g"), prob.get("x"))));
    ((DM) Jf.invoke(toDM(w0))).print_dense();

    Function solver_ref = ca.qpsol("solver", "qpoases", prob);
    Map<String, Object> options = map("hpipm", map("iter_max", 100, "res_g_max", 1e-10, "res_b_max", 1e-10, "res_d_max", 1e-10, "res_m_max", 1e-10), "dump_in", true);
    Function solver = ca.qpsol("solver", "hpipm", prob, options);
    Map<String, DM> args = Map.of("x0", toDM(w0), "lbx", toDM(lbw), "ubx", toDM(ubw), "lbg", toDM(lbg), "ubg", toDM(ubg));
    Map<String, DM> sol_ref = solver_ref.call(args);
    Map<String, DM> sol = solver.call(args);
    checkCodegen(solver, args, map("std", "c99", "extralibs", List.of("hpipm", "blasfeo"), "extra_options", fatrop_flags, "digits", 10));
    checkarray(sol_ref.get("x"), sol.get("x"));
    checkarrayDigits(sol_ref.get("lam_g"), sol.get("lam_g"), 8, "");
    checkarrayDigits(sol_ref.get("lam_x"), sol.get("lam_x"), 8, "");
    checkarrayDigits(sol_ref.get("f"), sol.get("f"), 8, "");

    solver = ca.nlpsol("solver", "sqpmethod", prob, map("qpsol", "hpipm", "qpsol_options", options));
    sol = solver.call(args);
    checkarray(sol_ref.get("x"), sol.get("x"));
    checkarrayDigits(sol_ref.get("lam_g"), sol.get("lam_g"), 8, "");
    checkarrayDigits(sol_ref.get("lam_x"), sol.get("lam_x"), 8, "");
    checkarrayDigits(sol_ref.get("f"), sol.get("f"), 8, "");
  }

  private static DM toDM(List<Double> v) {
    return new DM(v.stream().mapToDouble(Double::doubleValue).toArray());
  }

  private static DM mat(String a) {
    List<double[]> rows = new ArrayList<>();
    for (String i : a.split("\n")) {
      if (i.isEmpty()) continue;
      String[] e = i.split("\t", -1);
      double[] r = new double[e.length];
      for (int k = 0; k < e.length; k++) r[k] = e[k].trim().isEmpty() ? 0 : Double.parseDouble(e[k].trim());
      rows.add(r);
    }
    return ca.sparsify(new DM(rows.toArray(new double[0][])));
  }

  private static DM vec(String a) {
    String[] e = a.split("\n");
    double[] r = new double[e.length];
    for (int k = 0; k < e.length; k++) r[k] = Double.parseDouble(e[k].trim());
    return new DM(r);
  }

  @Test
  public void test_hpipm_timevarying() {
    Assumptions.assumeTrue(hasConic("hpipm"));
    Assumptions.assumeTrue(hasConic("qpoases"));
    DM A = mat(String.join("\n",
        "1\t0.2\t1\t-1\t0\t0\t0\t0\t0\t0\t0\t0",
        "-0.1\t0.4\t0\t0\t-1\t0\t0\t0\t0\t0\t0\t0",
        "0.3\t0.2\t0\t0\t0\t-1\t0\t0\t0\t0\t0\t0",
        "2\t0\t0.3\t0\t0\t0\t0\t0\t0\t0\t0\t0",
        "1\t1\t0.4\t0\t0\t0\t0\t0\t0\t0\t0\t0",
        "\t0\t0\t1\t4\t2\t1\t0.3\t-1\t0\t0\t0",
        "\t0\t0\t3\t1\t0\t1\t0.2\t0\t-1\t0\t0",
        "\t0\t0\t1\t1\t1\t1\t1\t0\t0\t0\t0",
        "\t0\t0\t0\t0\t0\t0\t0\t2\t4\t0\t-1",
        "\t0\t0\t0\t0\t0\t0\t0\t2\t3\t1\t0",
        "\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t3"));
    DM H = mat(String.join("\n",
        "7\t0\t0.2\t0\t0\t0\t0\t0\t0\t0\t0\t0",
        "\t7\t0.3\t0\t0\t0\t0\t0\t0\t0\t0\t0",
        "0.2\t0.3\t1\t0\t0\t0\t0\t0\t0\t0\t0\t0",
        "\t0\t0\t3\t0\t0\t0\t1\t0\t0\t0\t0",
        "\t0\t0\t0\t2\t0.1\t0\t0.7\t0\t0\t0\t0",
        "\t0\t0\t0\t0.1\t1\t0\t1\t0\t0\t0\t0",
        "\t0\t0\t0\t0\t0\t1\t0.1\t0\t0\t0\t0",
        "\t0\t0\t1\t0.7\t1\t0.1\t2\t0\t0\t0\t0",
        "\t0\t0\t0\t0\t0\t0\t0\t6\t0\t1\t0",
        "\t0\t0\t0\t0\t0\t0\t0\t0\t6\t0\t0",
        "\t0\t0\t0\t0\t0\t0\t0\t1\t0\t4\t0",
        "\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t9", ""));
    Map<String, Object> options = map("hpipm", map("iter_max", 100, "res_g_max", 1e-10, "res_b_max", 1e-10, "res_d_max", 1e-10, "res_m_max", 1e-10));
    Function solver = ca.conic("solver", "hpipm", hA(H, A), options);
    Function solver_ref = ca.conic("solver", "qpoases", hA(H, A));

    DM g = vec("1\n1\n0.2\n0.4\n1\n0.5\n0.3\n1\n0.6\n1\n1\n0.7");
    DM lbg = vec("0\n0\n0\n-2\n-2\n0\n0\n-2\n0\n-2\n-2");
    DM ubg = vec("0\n0\n0\n2\n2\n0\n0\n2\n0\n2\n2");
    DM lbx = vec("0.5\n0.2\n-1\n-1\n-1\n-1\n-1\n-1\n-1\n-1\n-1\n-1");
    DM ubx = vec("0.5\n0.2\n1\n1\n1\n1\n1\n1\n1\n1\n1\n1");

    Map<String, DM> args = new LinkedHashMap<>(Map.of("a", A, "h", H, "lba", lbg, "uba", ubg, "g", g, "lbx", lbx, "ubx", ubx));
    Map<String, DM> sol = solver.call(args);
    Map<String, DM> sol_ref = solver_ref.call(args);
    checkarrayDigits(sol_ref.get("x"), sol.get("x"), 7, "");
    checkarrayDigits(sol_ref.get("lam_a"), sol.get("lam_a"), 8, "");
    checkarrayDigits(sol_ref.get("lam_x"), sol.get("lam_x"), 8, "");
    checkCodegen(solver, args, map("std", "c99", "extralibs", List.of("hpipm", "blasfeo"), "extra_options", fatrop_flags, "digits", 9));

    solver = ca.conic("solver", "hpipm", hA(H, A), options);
    Map<String, DM> args2 = new LinkedHashMap<>(args);
    args2.put("x0", sol.get("x"));
    args2.put("lam_a0", sol.get("lam_a"));
    args2.put("lam_x0", sol.get("lam_x"));
    sol = solver.call(args2);
    checkarrayDigits(sol_ref.get("x"), sol.get("x"), 7, "");
    checkarrayDigits(sol_ref.get("lam_a"), sol.get("lam_a"), 8, "");
    checkarrayDigits(sol_ref.get("lam_x"), sol.get("lam_x"), 8, "");

    Map<String, DM> args3 = new LinkedHashMap<>(args);
    args3.put("x0", sol.get("x"));
    args3.put("lam_a0", sol.get("lam_a"));
    args3.put("lam_x0", sol.get("lam_x"));
    Map<String, DM> sol2 = solver.call(args3);
    checkarray(sol.get("x"), sol2.get("x"));
    System.out.println("codegen starts here");
    checkCodegen(solver, args3, map("std", "c99", "extralibs", List.of("hpipm", "blasfeo"), "extra_options", fatrop_flags, "digits", 9));
  }

  @Test
  public void test_SOCP() {
    Assumptions.assumeTrue(hasNlpsol("ipopt"));
    for (Conic c : conics) {
      // Python mutates the shared table here (verbose leaks into later tests); kept local
      Map<String, Object> qp_options = c.opts();
      qp_options.put("verbose", true);
      MX x = MX.sym("x");
      MX y = MX.sym("y");
      MX z = MX.sym("z");
      if (!c.has("soc")) continue;
      System.out.println("test_SOCP " + c.name + " " + qp_options);

      MX h = ca.soc(ca.vertcat(List.of(x.minus(5), y.minus(7))), new MX(4));
      Function solver = ca.qpsol("msyolver", c.name, Map.of("h", h, "x", ca.vertcat(List.of(x, y)), "f", x.times(2).plus(y)), qp_options);
      Map<String, DM> res = solver.call(new HashMap<String, DM>());
      checkarrayDigits(res.get("x"), dm(5 - 8 / Math.sqrt(5), 7 - 4 / Math.sqrt(5)), 7, c.name);
      checkarrayDigits(res.get("f"), new DM(10 - 16 / Math.sqrt(5) + 7 - 4 / Math.sqrt(5)), 7, c.name);

      double a = 1.3;
      h = ca.soc(ca.vertcat(List.of(x, y)), x.times(a).plus(4));
      MX xy = ca.vertcat(List.of(x, y));
      Function nsolver = ca.nlpsol("mysolver", "ipopt", Map.of("f", x.times(2).plus(y), "x", xy, "g", ca.dot(xy, xy).minus(ca.power(x.times(a).plus(4), 2))));
      nsolver.call(Map.of("ubg", new DM(0)));
      solver = ca.qpsol("msyolver", c.name, Map.of("h", h, "x", xy, "f", x.times(2).plus(y)), qp_options);
      res = solver.call(new HashMap<String, DM>());
      double xs = -(8 * Math.sqrt(5 - a * a) + 4 * Math.pow(a, 3) - 20 * a) / (Math.pow(a, 4) - 6 * a * a + 5);
      double ys = 4 * Math.sqrt(5 - a * a) / (a * a - 5);
      checkarrayDigits(res.get("f"), new DM(2 * xs + ys), 5, c.name);
      checkarrayDigits(res.get("x"), dm(xs, ys), List.of("cplex", "gurobi", "xpress").contains(c.name) ? 4 : 5, c.name);

      h = ca.diagcat(List.of(ca.soc(ca.vertcat(List.of(x.minus(5), y.minus(7))), new MX(4)), ca.soc(ca.vertcat(List.of(x.div(6), y.div(5))), new MX(1))));
      solver = ca.qpsol("msyolver", c.name, Map.of("h", h, "x", xy, "f", x.times(2).plus(y)), qp_options);
      res = solver.call(new HashMap<String, DM>());
      checkarrayDigits(res.get("x"), dm(1.655450403084473, 4.805919456574478), 5, c.name);

      h = ca.diagcat(List.of(
          ca.soc(ca.vertcat(List.of(x.times(-13).plus(y.times(3)).plus(z.times(5)).minus(3), x.times(-12).plus(y.times(12)).minus(z.times(6)).minus(2))),
              x.times(-12).minus(y.times(6)).plus(z.times(5)).minus(12)),
          ca.soc(ca.vertcat(List.of(x.times(-3).plus(y.times(6)).plus(z.times(2)), x.plus(y.times(9)).plus(z.times(2)).plus(3), x.unaryMinus().minus(y.times(19)).plus(z.times(3)).minus(42))),
              x.times(-3).plus(y.times(6)).minus(z.times(10)).plus(27))));
      solver = ca.qpsol("msyolver", c.name, Map.of("h", h, "x", ca.vertcat(List.of(x, y, z)), "f", x.times(-2).plus(y).plus(z.times(5))), qp_options);
      try {
        res = solver.call(new HashMap<String, DM>());
      } catch (RuntimeException e) {
        System.out.println(solver.stats());
      }
      checkarrayDigits(res.get("x"), dm(-5.0147928622, -5.766930599, -8.52180472), 4, c.name);
      assertTrue((Boolean) solver.stats().get("success"));

      MX x4 = MX.sym("x", 4);
      DM H = dm(new double[][] {
        {2.834044009405148, 0.867080384259271, 0.396881849048015, 0.506784822363357},
        {0.867080384259271, 2.184677189537596, 0.725076945381028, 1.223678163433993},
        {0.396881849048015, 0.725076945381028, 2.838389028806589, 0.712607093594686},
        {0.506784822363357, 1.223678163433993, 0.712607093594686, 3.340935020356804}});
      DM x0 = dm(1, 2, 3, 4);
      MX f = ca.mtimes(ca.mtimes(x4.T(), H), x0).unaryMinus();
      List<Object> ldl = ca.ldl(H);
      DM D = (DM) ldl.get(0), Lt = (DM) ldl.get(1);
      DM F = ca.mtimes(ca.sqrt(ca.diag(D)), DM.eye(4).plus(Lt));
      h = ca.soc(ca.vertcat(List.of(ca.mtimes(F, x4).times(Math.sqrt(2)), y.unaryMinus().plus(1))), y.plus(1));
      solver = ca.qpsol("msyolver", c.name, Map.of("x", ca.vertcat(List.of(x4, y)), "f", y.plus(f), "h", h), qp_options);
      res = solver.call(Map.of("lbx", dm(-inf, -inf, -inf, -inf, 1)));
      checkarrayDigits(res.get("x").get(new Slice(0, -1)), x0, 4, c.name);
    }
  }

  @Test
  public void test_no_success() {
    SX x = SX.sym("x");
    SX y = SX.sym("y");
    SX f = x.minus(y);
    for (Conic c : conics) {
      System.out.println("test_no_success " + c.name + " " + c.options);
      Map<String, Object> opts = c.opts();
      opts.put("error_on_fail", false);
      Map<String, SX> qp = Map.of("x", ca.vertcat(List.of(x, y)), "f", f, "g", ca.vertcat(List.of(x.plus(1), x.minus(2))));
      Map<String, DM> args = Map.of("x0", new DM(0), "lbg", new DM(0), "ubg", new DM(0), "lbx", dm(-10, -10), "ubx", dm(10, 10));
      Function solver = ca.qpsol("solver", c.name, qp, opts);
      solver.call(args);
      assertFalse((Boolean) solver.stats().get("success"), c.name);

      opts.put("error_on_fail", true);
      Function solver2 = ca.qpsol("solver", c.name, qp, opts);
      assertInException("process", () -> solver2.call(args));
    }
  }

  @Test
  public void test_expand() {
    Assumptions.assumeTrue(hasConic("osqp"));
    MX x = MX.sym("x");
    Map<String, MX> qp = Map.of("x", x, "f", x.minus(1));
    Function solver = ca.qpsol("solver", "osqp", qp);
    Function solver_qp = solver.find_function("solver_qp");
    assertTrue(solver_qp.is_a("MXFunction"));
    solver = ca.qpsol("solver", "osqp", qp, Map.of("expand", true));
    solver_qp = solver.find_function("solver_qp");
    assertTrue(solver_qp.is_a("SXFunction"));
  }

  @Test
  public void test_ipqp_status_init() {
    Assumptions.assumeTrue(hasConic("ipqp"));
    DM H = dm(new double[][] {{1, -1}, {-1, 2}});
    DM A = dm(new double[][] {{1, 1}, {-1, 2}, {2, 1}});
    Function solver = ca.conic("solver", "ipqp", hA(H, A), Map.of("print_iter", false, "print_header", false, "print_info", false));
    Map<String, DM> res = solver.call(Map.of("h", H, "g", dm(-2, -6), "a", A, "lbx", new DM(0), "ubx", new DM(inf), "lba", new DM(-inf), "uba", dm(2, 2, 3)));
    assertTrue((Boolean) solver.stats().get("success"));
    checkarrayDigits(res.get("x"), dm(2. / 3, 4. / 3), 6, "");
  }

  @Test
  public void test_qrqp_iter_count() {
    Assumptions.assumeTrue(hasConic("qrqp"));
    DM H = dm(new double[][] {{1, -1}, {-1, 2}});
    DM A = dm(new double[][] {{1, 1}, {-1, 2}, {2, 1}});
    Function solver = ca.conic("solver", "qrqp", hA(H, A), Map.of("print_iter", false, "print_header", false, "print_info", false));
    solver.call(Map.of("h", H, "g", dm(-2, -6), "a", A, "lbx", new DM(0), "ubx", new DM(inf), "lba", new DM(-inf), "uba", dm(2, 2, 3)));
    assertTrue((Boolean) solver.stats().get("success"));
    assertTrue((Long) solver.stats().get("iter_count") > 0);
  }

  @Test
  public void test_proxqp_mixed_multipliers() {
    Assumptions.assumeTrue(hasConic("proxqp"));
    DM H = DM.eye(3);
    DM A = dm(new double[][] {{1, 1, 0}, {0, 1, -1}, {1, 0, 1}});
    for (String backend : List.of("sparse", "dense")) {
      Function solver = ca.conic("solver", "proxqp", hA(H, A), Map.of("proxqp", Map.of("eps_abs", 1e-11, "backend", backend)));
      Map<String, DM> res = solver.call(Map.of("h", H, "g", dm(0, -2, 1), "a", A, "lba", dm(-inf, -0.5, 1.8), "uba", dm(1.5, -0.5, inf)));
      checkarrayDigits(res.get("x"), dm(0.7, 0.6, 1.1), 6, backend);
      checkarrayDigits(res.get("lam_a"), dm(0, 1.4, -0.7), 6, backend);
    }
  }

  private static MX m(double v) {
    return new MX(v);
  }

  @Test
  public void test_sensitivity() {
    Assumptions.assumeTrue(hasConic("qrqp"));
    MX x = MX.sym("x", 4);
    MX p = MX.sym("p", 2);
    MX p0 = p.get(0), p1 = p.get(1);
    MX H = ca.blockcat(List.of(
        List.of(p0.plus(2), m(0.5), m(0), m(0)),
        List.of(m(0.5), m(3), p1, m(0)),
        List.of(m(0), p1, m(4), m(0)),
        List.of(m(0), m(0), m(0), m(1))));
    MX f = ca.bilin(H, x, x).times(0.5).plus(ca.dot(ca.vertcat(List.of(p0.times(-5), m(2), p1.minus(5), m(1))), x));
    MX g = ca.vertcat(List.of(x.get(0).plus(p1.times(x.get(1))).plus(x.get(2)), x.get(0).minus(x.get(2)).plus(x.get(3)), x.get(1).plus(p0.times(x.get(2)))));
    Map<String, MX> qp = Map.of("x", x, "p", p, "f", f, "g", g);
    Map<String, MX> bounds = new LinkedHashMap<>();
    bounds.put("lbg", ca.vertcat(List.of(m(-inf), p1.minus(0.2), p0.plus(0.3))));
    bounds.put("ubg", ca.vertcat(List.of(p0.plus(1), p1.minus(0.2), m(inf))));
    bounds.put("lbx", ca.vertcat(List.of(m(-inf), p1.times(1.5), m(-inf), m(-inf))));
    bounds.put("ubx", m(inf));
    DM pv = dm(0.7, 0.3);

    Function ref = ca.nlpsol("ref", "sqpmethod", qp, map("qpsol", "qrqp", "qpsol_options", map("print_iter", false, "print_header", false),
        "print_iteration", false, "print_header", false, "print_status", false, "tol_pr", 1e-14, "tol_du", 1e-14));
    Map<String, MX> args = new LinkedHashMap<>(bounds);
    args.put("p", p);
    Map<String, MX> r = ref.call(args);
    Function R = new Function("R", List.of(p), List.of(r.get("x"), r.get("f"), r.get("lam_g"), r.get("lam_x")));
    DM lam_g = R.call(dms(pv)).get(2);
    assertTrue(at(lam_g, 0) > 1e-3 && at(lam_g, 2) < -1e-3 && Math.abs(at(lam_g, 1)) > 1e-3);
    assertTrue(at(R.call(dms(pv)).get(3), 1) < -1e-3);

    for (Conic c : conics) {
      if (!c.has("quadratic") || !c.has("dual")) continue;
      if (c.name.equals("nlpsol")) continue;
      System.out.println("test_sensitivity " + c.name + " " + c.options);
      int less_digits = c.lessDigits();
      int digits = List.of("qrqp", "qpoases", "daqp", "ipqp").contains(c.name) ? 7 : 4;
      Map<String, Object> o = c.opts();
      o.put("sens_min_lam", 1e-6);
      Function solver = ca.qpsol("solver", c.name, qp, o);
      Map<String, MX> rs = solver.call(args);
      Function F = new Function("F", List.of(p), List.of(rs.get("x"), rs.get("f"), rs.get("lam_g"), rs.get("lam_x")));
      checkfunction(F, R, dms(pv), 0.5 * Math.pow(10, -(digits - less_digits)));
    }
  }

  @Test
  public void test_sensitivity_equality_zero_multiplier() {
    MX x = MX.sym("x", 3);
    MX p = MX.sym("p");
    Map<String, MX> qp = Map.of("x", x, "p", p, "f", ca.sumsqr(x.minus(p)), "g", x.get(1).minus(x.get(0)));
    Function R = new Function("R", List.of(p), List.of(ca.vertcat(List.of(p, p, p))));
    for (Conic c : conics) {
      if (!c.has("quadratic") || !c.has("dual")) continue;
      if (c.name.equals("nlpsol") || c.name.equals("superscs")) continue;
      System.out.println("test_sensitivity_equality_zero_multiplier " + c.name + " " + c.options);
      int digits = List.of("qrqp", "qpoases", "daqp", "ipqp").contains(c.name) ? 7 : 4;
      Map<String, Object> o = c.opts();
      o.put("sens_min_lam", 1e-6);
      o.put("sens_force_equality", true);
      Function solver = ca.qpsol("solver", c.name, qp, o);
      Map<String, MX> args = new LinkedHashMap<>();
      args.put("p", p);
      args.put("lbg", m(0));
      args.put("ubg", m(0));
      args.put("lbx", ca.vertcat(List.of(m(-inf), m(-inf), p)));
      args.put("ubx", ca.vertcat(List.of(m(inf), m(inf), p)));
      Map<String, MX> r = solver.call(args);
      Function F = new Function("F", List.of(p), List.of(r.get("x")));
      checkfunction(F, R, dms(new DM(0.5)), 0.5 * Math.pow(10, -(digits - c.lessDigits())));
    }
  }

  @Test
  public void test_postpone_expand() {
    Assumptions.assumeTrue(hasConic("osqp"));
    MX x = MX.sym("x");
    Function J = new Function("jac_f", List.of(x, new MX(1, 1)), List.of(x.times(0.5)), List.of("x", "out_r"), List.of("jac_r_x"), Map.of("always_inline", true));
    Function f = new Function("f", List.of(x), List.of(ca.power(x, 2)), List.of("x"), List.of("r"), Map.of("custom_jacobian", J, "jac_penalty", 0));
    Map<String, MX> qp = Map.of("x", x, "f", (MX) f.invoke(x.minus(1)));

    Object[][] cases = {
      {Map.of("print_time", true), 0.5, "MXFunction"},
      {Map.of("expand", true, "print_time", true), 2.0, "SXFunction"},
      {Map.of("print_time", true, "postpone_expand", true), 0.5, "MXFunction"},
      {Map.of("expand", true, "print_time", true, "postpone_expand", true), 0.5, "SXFunction"},
    };
    for (Object[] cs : cases) {
      @SuppressWarnings("unchecked")
      Map<String, Object> o = (Map<String, Object>) cs[0];
      Function solver = ca.qpsol("solver", "osqp", qp, o);
      solver.call(Map.of("x0", new DM(2)));
      Function solver_qp = solver.find_function("solver_qp");
      checkarray(solver_qp.call(Map.of("x", new DM(2))).get("H"), (Double) cs[1]);
      assertTrue(solver_qp.is_a((String) cs[2]));
    }
  }

  // XpressSpecificTests (class-level @requires_conic("xpress")), folded into this class

  @Test
  public void test_mip_warmstart() throws IOException {
    Assumptions.assumeTrue(hasConic("xpress"));
    MX x = MX.sym("x", 3);
    Path logfile = Files.createTempFile("casadi", ".log");
    Files.delete(logfile);
    try {
      Function solver = ca.qpsol("s", "xpress", Map.of("x", x, "f", ca.sum1(x).unaryMinus(), "g", x.get(0).plus(x.get(1)).plus(x.get(2))),
          map("discrete", new boolean[] {true, true, true}, "mip_start", true, "log_file", logfile.toString(), "xpress", map("OUTPUTLOG", 1)));
      Map<String, DM> sol = solver.call(Map.of("x0", dm(1, 1, 0), "lbx", dm(0, 0, 0), "ubx", dm(1, 1, 1), "lbg", new DM(-inf), "ubg", new DM(2)));
      assertTrue((Boolean) solver.stats().get("success"), "MIP with mip_start should succeed");
      assertAlmostEqual(sol.get("f").to_double(), -2.0, 6, "MILP optimal cost should be -2");
      String log = Files.readString(logfile, StandardCharsets.UTF_8).replace("\r\n", "\n");
      assertTrue(log.contains("User solution"), "Xpress log must contain 'User solution' when mip_start=True");
      solver.close();
    } finally {
      Files.deleteIfExists(logfile);
    }
  }

  @Test
  public void test_iis_infeasible_lp() {
    Assumptions.assumeTrue(hasConic("xpress"));
    Assumptions.assumeFalse(isDarwin, "Xpress XPRSiisfirst corrupts the heap on macOS");
    MX x = MX.sym("x");
    Function solver = ca.qpsol("s", "xpress", Map.of("x", x, "f", x, "g", ca.vertcat(List.of(x, x.unaryMinus()))),
        map("error_on_fail", false, "compute_iis", true, "xpress", map("OUTPUTLOG", 0)));
    Map<String, DM> sol = solver.call(Map.of("x0", new DM(0), "lbx", new DM(-1e10), "ubx", new DM(1e10), "lbg", dm(5, 5), "ubg", dm(1e10, 1e10)));
    Map<String, Object> stats = solver.stats();
    assertFalse((Boolean) stats.get("success"), "Solver should report failure for infeasible LP");
    assertTrue(stats.containsKey("iis_rows"), "get_stats() must contain 'iis_rows' for infeasible LP");
    assertTrue(stats.containsKey("iis_cols"), "get_stats() must contain 'iis_cols' for infeasible LP");
    assertTrue(((List<?>) stats.get("iis_rows")).size() > 0, "IIS must identify at least one infeasible row");
    for (double[] row : sol.get("x").full())
      for (double v : row) assertTrue(Double.isNaN(v), "x must be NaN when no primal solution is available");
  }

  @Test
  public void test_iis_infeasible_milp() {
    Assumptions.assumeTrue(hasConic("xpress"));
    Assumptions.assumeFalse(isDarwin, "Xpress XPRSiisfirst corrupts the heap on macOS");
    MX x = MX.sym("x");
    Function solver = ca.qpsol("s", "xpress", Map.of("x", x, "f", x, "g", x),
        map("discrete", new boolean[] {true}, "error_on_fail", false, "compute_iis", true, "xpress", map("OUTPUTLOG", 0)));
    solver.call(Map.of("x0", new DM(0), "lbx", dm(0), "ubx", dm(1), "lbg", new DM(2), "ubg", new DM(1e10)));
    Map<String, Object> stats = solver.stats();
    assertFalse((Boolean) stats.get("success"));
    assertTrue(stats.containsKey("iis_rows"));
    assertTrue(stats.containsKey("iis_cols"));
    assertTrue(((List<?>) stats.get("iis_rows")).size() + ((List<?>) stats.get("iis_cols")).size() > 0, "IIS must contain at least one row or column");
  }
}

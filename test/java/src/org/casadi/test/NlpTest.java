package org.casadi.test;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.UUID;

import org.casadi.CasadiMeta;
import org.casadi.Callback;
import org.casadi.CodeGenerator;
import org.casadi.DM;
import org.casadi.Function;
import org.casadi.GenericMatrix;
import org.casadi.GlobalOptions;
import org.casadi.MX;
import org.casadi.Opti;
import org.casadi.OptiCallback;
import org.casadi.SX;
import org.casadi.Slice;
import org.casadi.Sparsity;
import org.casadi.UnifiedReturnStatus;
import org.casadi.ca;
import org.junit.jupiter.api.Assumptions;
import org.junit.jupiter.api.Test;

/** Port of test/python/nlp.py */
public class NlpTest extends CasadiTestCase {

  private static final double inf = Double.POSITIVE_INFINITY;
  private static final double pi = Math.PI;
  private static final int codegen_check_digits = 15;

  /** One row of the module-level solvers table: (Solver, solver_options, aux_options) */
  private static final class S {
    final String name;
    final Map<String, Object> opts;
    final Map<String, Object> codegen;
    final boolean discrete;

    S(String name, Map<String, Object> opts, Map<String, Object> codegen, boolean discrete) {
      this.name = name;
      this.opts = opts;
      this.codegen = codegen;
      this.discrete = discrete;
    }

    /** Python: s in str(solver_options) */
    boolean in(String s) {
      return opts.toString().contains(s);
    }

    @Override
    public String toString() {
      return name + " " + opts;
    }
  }

  private static boolean env(String s) {
    return System.getenv(s) != null;
  }

  private static final List<S> solvers = new ArrayList<>();

  static {
    if (!env("SKIP_WORHP_TESTS") && hasNlpsol("worhp"))
      solvers.add(new S("worhp", d("worhp", d("TolOpti", 1e-9)), null, false));
    if (!env("SKIP_FATROP_TESTS") && hasNlpsol("fatrop"))
      solvers.add(new S("fatrop", d("fatrop", d()), d("std", "c99", "extralibs", List.of("fatrop", "blasfeo"),
          "extra_options", List.of("-Wno-strict-prototypes")), false));
    if (!env("SKIP_SLEQP_TESTS") && hasNlpsol("sleqp"))
      solvers.add(new S("sleqp", d("print_time", false, "sleqp", d("linesearch", "Approx", "feas_tol", 1e-7, "stat_tol",
          1e-7, "slack_tol", 1e-7, "hess_eval", "Exact")), null, false));
    if (!env("SKIP_UNO_TESTS") && hasNlpsol("uno")) {
      Map<String, Object> uno_codegen = d("std", "c99", "extralibs", List.of("uno"), "extra_include", List.of("uno"));
      for (Map<String, Object> u : List.of(
          d("preset", "ipopt", "primal_tolerance", 1e-8, "dual_tolerance", 1e-8),
          d("preset", "funnelsqp", "primal_tolerance", 1e-8, "dual_tolerance", 1e-8),
          d("preset", "filtersqp", "primal_tolerance", 1e-8, "dual_tolerance", 1e-8),
          d("preset", "filtersqp", "primal_tolerance", 1e-8, "dual_tolerance", 1e-8, "hessian_model", "LBFGS")))
        solvers.add(new S("uno", d("print_time", false, "uno", u), uno_codegen, false));
    }
    if (!env("SKIP_ALPAQA_TESTS") && hasNlpsol("alpaqa"))
      solvers.add(new S("alpaqa", d("print_time", false, "alpaqa", d("alm.tolerance", 1e-10, "alm.dual_tolerance", 1e-10,
          "alm.penalty_update_factor", 10, "alm.max_iter", 3000, "alm.print_interval", 1, "panoc.max_iter", 500,
          "panoc.print_interval", 1, "lbfgs.memory", 2)), null, false));
    if (!env("SKIP_IPOPT_TESTS") && hasNlpsol("ipopt")) {
      Map<String, Object> codegen = d("extralibs", List.of("ipopt"), "std", "c99");
      solvers.add(new S("ipopt", d("print_time", false, "ipopt", d("tol", 1e-10, "derivative_test", "second-order",
          "print_level", 0)), codegen, false));
      solvers.add(new S("ipopt", d("print_time", false, "ipopt", d("tol", 1e-10, "derivative_test", "first-order",
          "hessian_approximation", "limited-memory", "print_level", 0)), codegen, false));
    }
    if (!env("SKIP_IPOPT_TESTS") && hasNlpsol("ipopt") && hasNlpsol("sqpmethod")) {
      Map<String, Object> qpsol_options = d("nlpsol", "ipopt", "nlpsol_options", d("ipopt.tol", 1e-12,
          "ipopt.tiny_step_tol", 1e-20, "ipopt.fixed_variable_treatment", "make_constraint", "ipopt.print_level", 0,
          "print_time", false));
      solvers.add(new S("sqpmethod", d("qpsol", "nlpsol", "qpsol_options", qpsol_options, "print_header", false,
          "print_iteration", true, "print_time", false), null, false));
      solvers.add(new S("sqpmethod", d("qpsol", "nlpsol", "qpsol_options", qpsol_options, "hessian_approximation",
          "limited-memory", "tol_du", 1e-10, "tol_pr", 1e-10, "min_step_size", 1e-14, "print_header", false,
          "print_iteration", true, "print_time", false, "max_iter", 55), null, false));
    }
    if (!env("SKIP_QRQP_TESTS") && hasConic("qrqp") && hasNlpsol("sqpmethod")) {
      Map<String, Object> codegen = d("std", "c99");
      Map<String, Object> qpsol_options = d("print_iter", false, "print_header", false, "error_on_fail", false);
      solvers.add(new S("sqpmethod", d("qpsol", "qrqp", "qpsol_options", qpsol_options, "print_header", false,
          "print_iteration", false, "print_time", false), codegen, false));
      solvers.add(new S("sqpmethod", d("qpsol", "qrqp", "max_iter_ls", 0, "qpsol_options", qpsol_options,
          "print_header", false, "print_iteration", false, "print_time", false), codegen, false));
      solvers.add(new S("sqpmethod", d("qpsol", "qrqp", "convexify_strategy", "regularize", "max_iter", 500,
          "qpsol_options", qpsol_options, "print_header", false, "print_iteration", true, "print_time", false,
          "tol_du", 1e-8, "min_step_size", 1e-12), codegen, false));
    }
    if (!env("SKIP_DAQP_TESTS") && hasConic("daqp") && hasNlpsol("sqpmethod"))
      solvers.add(new S("sqpmethod", d("qpsol", "daqp", "hessian_approximation", "limited-memory", "tol_du", 1e-10,
          "tol_pr", 1e-10, "min_step_size", 1e-14, "max_iter", 55, "print_header", false, "print_iteration", true,
          "print_time", false, "qpsol_options", d("print_problem", true, "print_out", true)),
          d("std", "c99", "extralibs", List.of("daqp")), false));
    if (!env("SKIP_BLOCKSQP_TESTS") && hasNlpsol("blocksqp")) {
      try {
        ca.load_linsol("ma27");
        solvers.add(new S("blocksqp", d(), null, false));
      } catch (RuntimeException e) {
        // ma27 unavailable
      }
    }
    if (!env("SKIP_BONMIN_TESTS") && hasNlpsol("bonmin")) solvers.add(new S("bonmin", d(), null, true));
    if (!env("SKIP_KNITRO_TESTS") && hasNlpsol("knitro"))
      solvers.add(new S("knitro", d("knitro", d("feastol", 1e-9, "opttol", 1e-9, "ftol", 1e-20)), null, false));
    if (!env("SKIP_SNOPT_TESTS") && hasNlpsol("snopt"))
      solvers.add(new S("snopt", d("snopt", d("Verify_level", 3, "Major_optimality_tolerance", 1e-12,
          "Minor_feasibility_tolerance", 1e-12, "Major_feasibility_tolerance", 1e-12)), null, false));
    if (!env("SKIP_MADNLP_TESTS") && hasNlpsol("madnlp"))
      solvers.add(new S("madnlp", d("madnlp", d()), libmadCodegen(), false));
    if (!env("SKIP_CONOPT_TESTS") && hasNlpsol("conopt"))
      solvers.add(new S("conopt", d("conopt", d("Tol_Feas_Min", 1e-9, "Tol_Optimality", 1e-9)), null, false));
  }

  private static Map<String, Object> libmadCodegen() {
    String libmad_dir = System.getenv().getOrDefault("LIBMADDIR", "/missing");
    return d("extralibs", List.of("Mad"), "std", "c99", "extralibdirs", List.of(libmad_dir + "/lib"),
        "extra_include", List.of(libmad_dir + "/include"));
  }

  private final List<Object> keep = new ArrayList<>();

  // ---------------------------------------------------------------- helpers

  /** Dict literal */
  private static Map<String, Object> d(Object... kv) {
    Map<String, Object> r = new LinkedHashMap<>();
    for (int i = 0; i < kv.length; i += 2) r.put((String) kv[i], kv[i + 1]);
    return r;
  }

  /** dict(base, **kv) */
  private static Map<String, Object> with(Map<String, Object> base, Object... kv) {
    Map<String, Object> r = new LinkedHashMap<>(base);
    r.putAll(d(kv));
    return r;
  }

  private static DM toDM(Object v) {
    if (v instanceof DM) return (DM) v;
    if (v instanceof Number) return new DM(((Number) v).doubleValue());
    if (v instanceof double[]) return new DM((double[]) v);
    throw new IllegalArgumentException(String.valueOf(v));
  }

  private static MX toMX(Object v) {
    if (v instanceof MX) return (MX) v;
    return new MX(toDM(v));
  }

  private static SX toSX(Object v) {
    if (v instanceof SX) return (SX) v;
    return new SX(toDM(v));
  }

  private static Map<String, DM> in(Object... kv) {
    Map<String, DM> r = new LinkedHashMap<>();
    for (int i = 0; i < kv.length; i += 2) r.put((String) kv[i], toDM(kv[i + 1]));
    return r;
  }

  private static Map<String, MX> min(Object... kv) {
    Map<String, MX> r = new LinkedHashMap<>();
    for (int i = 0; i < kv.length; i += 2) r.put((String) kv[i], toMX(kv[i + 1]));
    return r;
  }

  private static Map<String, SX> sxnlp(Object... kv) {
    Map<String, SX> r = new LinkedHashMap<>();
    for (int i = 0; i < kv.length; i += 2) r.put((String) kv[i], toSX(kv[i + 1]));
    return r;
  }

  private static Map<String, MX> mxnlp(Object... kv) {
    return min(kv);
  }

  private static double[] rep(double v, int n) {
    double[] r = new double[n];
    Arrays.fill(r, v);
    return r;
  }

  private static double[] range(int n) {
    double[] r = new double[n];
    for (int i = 0; i < n; i++) r[i] = i;
    return r;
  }

  /** Python x[k] on a DM, as a float */
  private static double el(DM x, int k) {
    return ca.densify(x).nonzeros()[k];
  }

  private static Object stat(Function f, String k) {
    return f.stats().get(k);
  }

  private static boolean success(Function f) {
    return (Boolean) stat(f, "success");
  }

  private static long iterCount(Function f) {
    return ((Number) stat(f, "iter_count")).longValue();
  }

  /** helpers.py assertAlmostEqual: places after scaling large values */
  private static void almost(double first, double second, int places, String msg) {
    if (Double.isNaN(first) && Double.isNaN(second)) return;
    if (first == second) return;
    double n = Math.max(Math.abs(first), Math.abs(second));
    n = n > 1e3 ? Math.pow(10, Math.floor(Math.log10(n))) : 1.0;
    double diff = Math.abs(first / n - second / n);
    assertTrue(Math.round(diff * Math.pow(10, places)) == 0,
        msg + " " + first + " <-> " + second + " (places=" + places + ")");
  }

  private static void almost(double first, double second, int places) {
    almost(first, second, places, "");
  }

  /** helpers.py checkarray with digits */
  private static void chk(DM a, DM b, int digits, String msg) {
    assertEquals(b.size1(), a.size1(), msg + ": size1 " + a + " vs " + b);
    assertEquals(b.size2(), a.size2(), msg + ": size2 " + a + " vs " + b);
    double[][] x = a.full(), y = b.full();
    for (int i = 0; i < x.length; i++)
      for (int j = 0; j < x[i].length; j++)
        almost(x[i][j], y[i][j], digits, msg + " evaluation error at (" + i + "," + j + ")\n" + a + "\n" + b);
  }

  private static void chk(DM a, DM b, int digits) {
    chk(a, b, digits, "");
  }

  private static void chk(DM a, DM b) {
    chk(a, b, 10, "");
  }

  /** Python list literal as a row, as checkarray(x, [..]) does */
  private static DM row(double... v) {
    return new DM(new double[][] {v});
  }

  private static void assertOutput(List<String> included, List<String> excluded, Runnable r) {
    String out = captureStdout(r);
    for (String s : included) assertTrue(out.contains(s), "expected '" + s + "' in output:\n" + out);
    for (String s : excluded) assertFalse(out.contains(s), "unexpected '" + s + "' in output:\n" + out);
  }

  /** Positional inputs from a dict, missing ones 0, as helpers.py checkfunction */
  private static List<DM> positional(Function f, Map<String, DM> d) {
    List<DM> r = new ArrayList<>();
    for (String n : f.name_in()) r.add(d.containsKey(n) ? d.get(n) : new DM(0.0));
    return r;
  }

  /** helpers.py checkfunction_light: repeated evaluation agrees */
  private static void checkfunctionLight(Function trial, Function solution, List<DM> inputs, int digits, String failmessage) {
    for (int k = 0; k < 2; k++) {
      List<DM> a = trial.call(inputs), b = solution.call(inputs);
      assertEquals(solution.n_out(), trial.n_out(), failmessage);
      assertEquals(solution.n_in(), trial.n_in(), failmessage);
      for (int i = 0; i < trial.n_out(); i++) {
        if (trial.sparsity_out(i).is_empty() || solution.sparsity_out(i).is_empty()) continue;
        chk(a.get(i), b.get(i), digits, failmessage + ": output(" + i + ": '" + trial.name_out(i) + "')");
      }
    }
  }

  private static void checkfunctionLight(Function trial, Function solution, List<DM> inputs, int digits) {
    checkfunctionLight(trial, solution, inputs, digits, "");
  }

  /** helpers.py check_serialize (without the python serialization_scheme validation) */
  private static void checkSerialize(Function f, Map<String, DM> inputs) {
    String s = f.serialize(d("debug", true));
    Function f2 = Function.deserialize(s);
    Map<String, DM> a = f.call(inputs), b = f2.call(inputs);
    assertEquals(f.name_out(), f2.name_out());
    for (String k : f.name_out()) chk(a.get(k), b.get(k), 16, k);
  }

  private Function checkCodegen(Function f, Map<String, DM> inputs, Map<String, Object> cg) {
    return checkCodegen(f, inputs, cg, 15, false, false, null);
  }

  private Function checkCodegen(Function f, Map<String, DM> inputs, Map<String, Object> cg, int digits) {
    return checkCodegen(f, inputs, cg, digits, false, false, null);
  }

  private static String sh(String cmd, File stdin, File stdout) {
    try {
      ProcessBuilder pb = new ProcessBuilder("sh", "-c", cmd).redirectErrorStream(stdout == null);
      if (stdin != null) pb.redirectInput(stdin);
      if (stdout != null) pb.redirectOutput(stdout);
      else pb.redirectOutput(ProcessBuilder.Redirect.INHERIT);
      return Integer.toString(pb.start().waitFor());
    } catch (IOException | InterruptedException e) {
      throw new RuntimeException(e);
    }
  }

  /** helpers.py check_codegen (gcc, linux): generate, compile, load as external and compare; run_slow only */
  @SuppressWarnings("unchecked")
  private Function checkCodegen(Function f, Map<String, DM> inputs, Map<String, Object> cg, int digits, boolean debug,
      boolean main, long[] main_return_code) {
    if (!RUN_SLOW) return null;
    String name = "codegen_" + UUID.randomUUID().toString().replace("-", "");
    CodeGenerator gen = new CodeGenerator(name, main ? d("main", true) : d());
    gen.add(f);
    gen.generate();
    new File(name + ".c").deleteOnExit();
    String libdir = GlobalOptions.getCasadiPath(), includedir = GlobalOptions.getCasadiIncludePath();
    List<String> inc = new ArrayList<>(List.of(includedir, includedir + "/highs"));
    for (String e : (List<String>) cg.getOrDefault("extra_include", List.of())) inc.add(includedir + "/" + e);
    List<String> libdirs = new ArrayList<>(List.of(libdir));
    libdirs.addAll((List<String>) cg.getOrDefault("extralibdirs", List.of()));
    // libraries that python_install bundles next to libcasadi live on LD_LIBRARY_PATH here
    String ld = System.getenv("LD_LIBRARY_PATH");
    if (ld != null)
      for (String p : ld.split(":"))
        if (!p.isEmpty()) {
          libdirs.add(p);
          inc.add(new File(p).getParent() + "/include");
        }
    StringBuilder flags = new StringBuilder();
    for (String l : libdirs) flags.append(" -L").append(l).append(" -Wl,-rpath,").append(l);
    StringBuilder extra = new StringBuilder();
    for (String l : (List<String>) cg.getOrDefault("extralibs", List.of())) extra.append(l.startsWith("-") || l.contains(".") ? " " + l : " -l" + l);
    for (String o : (List<String>) cg.getOrDefault("extra_options", List.of())) extra.append(" ").append(o);
    StringBuilder incs = new StringBuilder();
    for (String i : inc) incs.append(" -I").append(i);
    String std = (String) cg.getOrDefault("std", "c89");
    String opt = debug ? "-O0 -g" : "-O3";
    String common = " -pedantic -std=" + std + " -fPIC -Wall -Werror -Wextra" + incs
        + " -Wno-dangling-pointer -Wno-unknown-warning-option -Wno-unknown-pragmas -Wno-long-long -Wno-unused-parameter " + opt;
    String lib = "./" + name + ".so";
    sh("gcc" + common + " -shared " + name + ".c -o " + name + ".so" + flags + " -Wl,-rpath,." + extra, null, null);
    new File(lib).deleteOnExit();
    Function f2 = ca.external(f.name(), lib);
    if (main) {
      sh("gcc" + common + " " + name + ".c -o " + name + flags + " -Wl,-rpath,. -lm" + extra, null, null);
      new File(name).deleteOnExit();
      File fin = new File(f.name() + "_in.txt"), fout = new File(f.name() + "_out.txt");
      fin.deleteOnExit();
      fout.deleteOnExit();
      f.generate_in(fin.getPath(), f.convert_in(inputs));
      long rc = Long.parseLong(sh("./" + name + " " + f.name(), fin, fout));
      boolean ok = false;
      for (long c : main_return_code) ok |= c == rc;
      assertTrue(ok, "return code " + rc);
      if (rc != 0) return f2;
      List<DM> ref = f.call(f.convert_in(inputs));
      List<DM> outputs = f.generate_out(fout.getPath());
      for (int i = 0; i < f.n_out(); i++) chk(ref.get(i), outputs.get(i), digits);
    }
    Map<String, DM> a = f.call(inputs), b = f2.call(inputs);
    assertEquals(f.name_out(), f2.name_out());
    for (String k : f.name_out()) chk(a.get(k), b.get(k), digits, k);
    return f2;
  }

  /** Symmetric eigendecomposition (Jacobi): {eigenvalues as one row, eigenvectors as columns} */
  private static double[][][] eigh(double[][] A0) {
    int n = A0.length;
    double[][] A = new double[n][], V = new double[n][n];
    for (int i = 0; i < n; i++) {
      A[i] = A0[i].clone();
      V[i][i] = 1;
    }
    for (int sweep = 0; sweep < 100; sweep++) {
      double off = 0;
      for (int i = 0; i < n; i++)
        for (int j = i + 1; j < n; j++) off += A[i][j] * A[i][j];
      if (off < 1e-30) break;
      for (int p = 0; p < n; p++)
        for (int q = p + 1; q < n; q++) {
          if (A[p][q] == 0) continue;
          double theta = (A[q][q] - A[p][p]) / (2 * A[p][q]);
          double t = Math.signum(theta) / (Math.abs(theta) + Math.sqrt(theta * theta + 1));
          if (theta == 0) t = 1;
          double c = 1 / Math.sqrt(t * t + 1), s = t * c;
          for (int k = 0; k < n; k++) {
            double akp = A[k][p], akq = A[k][q];
            A[k][p] = c * akp - s * akq;
            A[k][q] = s * akp + c * akq;
          }
          for (int k = 0; k < n; k++) {
            double apk = A[p][k], aqk = A[q][k];
            A[p][k] = c * apk - s * aqk;
            A[q][k] = s * apk + c * aqk;
          }
          for (int k = 0; k < n; k++) {
            double vkp = V[k][p], vkq = V[k][q];
            V[k][p] = c * vkp - s * vkq;
            V[k][q] = s * vkp + c * vkq;
          }
        }
    }
    double[] dd = new double[n];
    for (int i = 0; i < n; i++) dd[i] = A[i][i];
    return new double[][][] {{dd}, V};
  }

  /** V diag(d) V^T */
  private static DM recompose(double[] dd, double[][] V) {
    int n = dd.length;
    double[][] r = new double[n][n];
    // numpy's np.dot(V, np.dot(D, V.T)) evaluation order; the clip cases are ill-conditioned
    for (int i = 0; i < n; i++)
      for (int j = 0; j < n; j++)
        for (int k = 0; k < n; k++) r[i][j] += V[i][k] * (dd[k] * V[j][k]);
    return new DM(r);
  }

  private static DM reflect(DM A, double eps) {
    double[][][] e = eigh(ca.densify(A).full());
    double[] dd = e[0][0].clone();
    for (int i = 0; i < dd.length; i++) dd[i] = Math.abs(dd[i]);
    return recompose(dd, e[1]);
  }

  private static DM clip(DM A, double eps) {
    double[][][] e = eigh(ca.densify(A).full());
    double[] dd = e[0][0].clone();
    for (int i = 0; i < dd.length; i++) if (dd[i] < eps) dd[i] = eps;
    return recompose(dd, e[1]);
  }

  private static SX sq(SX x) {
    return x.power(2.0);
  }

  private static MX sq(MX x) {
    return x.power(2.0);
  }

  /** (a-x)**2+100*(y-x**2)**2 */
  private static SX rosen(SX a, SX x, SX y) {
    return sq(a.minus(x)).plus(ca.times(100.0, sq(y.minus(sq(x)))));
  }

  private static SX rosen(double a, SX x, SX y) {
    return rosen(new SX(a), x, y);
  }

  private static SX vc(SX... x) {
    return ca.vertcat(List.of(x));
  }

  private static MX vc(MX... x) {
    return ca.vertcat(List.of(x));
  }

  private static DM vc(double... x) {
    return new DM(x);
  }

  // ---------------------------------------------------------------- tests

  @Test
  public void test_ccopt() {
    Assumptions.assumeTrue(hasNlpsol("ccopt"));
    SX x = SX.sym("x"), y = SX.sym("y");
    SX f = sq(x.minus(2.0)).plus(sq(y.minus(2.0)));
    Map<String, SX> nlp = sxnlp("x", vc(x, y), "f", f);
    Map<String, Object> solver_opts = d("cc_pairs", List.of(List.of(0, 1)), "madnlp", d("bound_relax_factor", 0.0));
    Function solver = ca.nlpsol("solver", "ccopt", nlp, solver_opts);
    Map<String, DM> solver_in = in("x0", vc(2.5, 3.0), "lbx", vc(0.0, 0.0));
    Map<String, DM> solver_out = solver.call(solver_in);
    assertTrue(success(solver));
    almost(solver_out.get("f").to_double(), 4.0, 6);
    double x_sol = el(solver_out.get("x"), 0), y_sol = el(solver_out.get("x"), 1);
    almost(Math.min(x_sol, y_sol), 0.0, 6);
    almost(Math.max(x_sol, y_sol), 2.0, 6);
    checkCodegen(solver, solver_in, libmadCodegen());
    checkSerialize(solver, solver_in);
  }

  @Test
  public void test_alpaqa() {
    Assumptions.assumeTrue(hasNlpsol("alpaqa"));
    SX x = SX.sym("x"), y = SX.sym("y");
    Map<String, SX> nlp = sxnlp("x", vc(x, y), "f", rosen(1, x, y), "g", x.plus(y));
    Function solver = ca.nlpsol("mysolver", "alpaqa", nlp);
    solver.call(in("lbx", rep(-10, 2), "ubx", rep(10, 2), "lbg", vc(-10), "ubg", vc(10)));
  }

  @Test
  public void test_affine_constant_terms() {
    SX x = SX.sym("x");
    for (S s : solvers) {
      Map<String, Object> solver_options = s.opts;
      if (s.name.equals("bonmin")) solver_options = with(solver_options, "calc_g", true);
      message("affine objective with a constant term: " + s.name);
      Function solver = ca.nlpsol("mysolver", s.name, sxnlp("x", x, "f", ca.times(2.0, x).plus(7.0)), solver_options);
      Map<String, DM> solver_out = solver.call(in("x0", 0, "lbx", -10, "ubx", 10));
      chk(solver_out.get("x"), vc(-10), 6);
      chk(solver_out.get("f"), vc(-13), 6);

      message("fully linear constraint with a constant term: " + s.name);
      solver = ca.nlpsol("mysolver", s.name, sxnlp("x", x, "f", sq(x.minus(3.0)), "g", x.plus(5.0)), solver_options);
      solver_out = solver.call(in("x0", 0, "lbg", 0, "ubg", 0));
      chk(solver_out.get("x"), vc(-5), 6);
      chk(solver_out.get("g"), vc(0), 6);

      message("fully linear range constraint with a constant term: " + s.name);
      solver = ca.nlpsol("mysolver", s.name, sxnlp("x", x, "f", sq(x.minus(3.0)), "g", x.plus(5.0)), solver_options);
      solver_out = solver.call(in("x0", 0, "lbg", 0, "ubg", 1));
      chk(solver_out.get("x"), vc(-4), 6);
      chk(solver_out.get("g"), vc(1), 6);
    }
  }

  @Test
  public void test_conopt_constant_constraints() {
    Assumptions.assumeTrue(hasNlpsol("conopt"));
    SX x = SX.sym("x"), p = SX.sym("p");
    for (SX other : List.of(SX.zeros(0, 1), x.times(x), x)) {
      SX g = vc(other, p);
      Function solver0 = ca.nlpsol("solver", "conopt", sxnlp("x", x, "p", p, "f", sq(x.minus(2.0)), "g", g),
          d("error_on_fail", false, "print_time", false));
      for (Function solver : List.of(solver0, Function.deserialize(solver0.serialize()))) {
        for (double[] c : new double[][] {{5, -10, 10}, {5, 0, 0}, {0, 0, 0}, {5, 4, 6}, {5, -inf, inf}}) {
          double value = c[0], lower = c[1], upper = c[2];
          int n = (int) other.numel();
          double[] lbg = rep(-10, n + 1), ubg = rep(10, n + 1);
          lbg[n] = lower;
          ubg[n] = upper;
          Map<String, DM> result = solver.call(in("x0", 1, "p", value, "lbg", lbg, "ubg", ubg));
          String sub = other + " " + value + " " + lower + " " + upper;
          assertEquals(lower <= value && value <= upper, success(solver), sub);
          chk(result.get("g").get(-1), new DM(value), 10, sub);
          if (success(solver)) chk(result.get("x"), new DM(2), 6, sub);
        }
      }
    }
  }

  @Test
  public void test_conopt_evaluation_errors() {
    Assumptions.assumeTrue(hasNlpsol("conopt"));
    SX x = SX.sym("x"), p = SX.sym("p");
    Object[][] cases = {
      {ca.sqrt(x.plus(p)), x, 0.0},
      {sq(x.minus(1.0)), ca.sqrt(x.plus(p)), 0.0},
      {p.times(x), x, inf},
      {sq(x.minus(1.0)), p.times(x), inf},
      {sq(x.minus(1.0)), p, Double.NaN},
    };
    for (Object[] c : cases) {
      SX f = (SX) c[0], g = (SX) c[1];
      double bad_p = (Double) c[2];
      Function solver = ca.nlpsol("solver", "conopt", sxnlp("x", x, "p", p, "f", f, "g", g),
          d("error_on_fail", false, "print_time", false));
      for (double value : new double[] {1, bad_p, 1}) {
        solver.call(in("x0", 0, "p", value, "lbx", 0, "ubx", 2, "lbg", -10, "ubg", 10));
        String sub = f + " " + g + " " + value;
        assertEquals(value == 1, success(solver), sub);
        assertEquals(value == 1 ? "SOLVER_RET_SUCCESS" : "SOLVER_RET_NAN", stat(solver, "unified_return_status"), sub);
      }
    }
  }

  @Test
  public void test_nonregular_point() {
    SX x = SX.sym("x");
    Map<String, SX> nlp = sxnlp("x", x, "f", sq(x.plus(1.0)), "g", ca.sqrt(x));
    for (S s : solvers) {
      if (List.of("madnlp", "fatrop").contains(s.name)) continue;
      message("test_nonregular_point " + s);
      Function solver = ca.nlpsol("mysolver", s.name, nlp, s.opts);
      try {
        solver.call(in("lbx", vc(-1000), "ubx", vc(1000), "lbg", vc(-1000), "ubg", vc(1000), "x0", 1e-4));
      } catch (RuntimeException e) {
        // ignored, as in python
      }
      if (!List.of("ipopt", "snopt", "blocksqp", "bonmin", "knitro", "sleqp", "alpaqa", "uno").contains(s.name))
        assertEquals("SOLVER_RET_NAN", stat(solver, "unified_return_status"), s.toString());
      assertFalse(success(solver), s.toString());
    }

    nlp = sxnlp("x", x, "f", sq(x), "g", ca.sqrt(x));
    for (S s : solvers) {
      if (List.of("madnlp", "fatrop").contains(s.name)) continue;
      message("test_nonregular_point " + s);
      Function solver = ca.nlpsol("mysolver", s.name, nlp, s.opts);
      Map<String, DM> solver_in = in("lbx", vc(-10), "ubx", vc(10), "lbg", vc(-2), "ubg", vc(2), "x0", -2);
      try {
        solver.call(solver_in);
        solver.call(solver_in);
      } catch (RuntimeException e) {
        // ignored, as in python
      }
      if (!List.of("ipopt", "snopt", "bonmin", "knitro", "sleqp", "alpaqa").contains(s.name))
        assertEquals("SOLVER_RET_NAN", stat(solver, "unified_return_status"), s.toString());
      assertFalse(success(solver), s.toString());
    }
  }

  private static final class Interrupt extends Callback {
    Interrupt(String name) {
      construct(name, d("enable_fd", true));
    }

    @Override
    public List<DM> eval(List<DM> argin) {
      throw new RuntimeException("KeyboardInterrupt");
    }
  }

  private static final class Raiser extends OptiCallback {
    @Override
    public void call(long i) {
      throw new RuntimeException("KeyboardInterrupt");
    }
  }

  // Python raises KeyboardInterrupt in the callbacks; Java throws a RuntimeException instead
  @Test
  public void test_iteration_interrupt() {
    Interrupt interrupt = new Interrupt("interrupt");
    MX x = MX.sym("x");
    Map<String, MX> nlp = mxnlp("x", x, "f", interrupt.call(List.of(x)).get(0), "g", x);
    for (S s : solvers) {
      if (s.name.equals("madnlp")) continue;
      Map<String, Object> solver_options = with(s.opts, "error_on_fail", true);
      Function solver = ca.nlpsol("solver", s.name, nlp, solver_options);
      assertInAnyOutput("KeyboardInterrupt", () -> solver.call(in("lbg", -5, "ubg", 5)));
      assertThrows(RuntimeException.class, () -> solver.call(in("lbg", -5, "ubg", 5)), s.toString());
    }

    for (S s : solvers) {
      if (List.of("worhp", "blocksqp", "knitro", "bonmin", "snopt", "alpaqa", "fatrop", "madnlp").contains(s.name)) continue;
      message("test_iteration_interrupt " + s);
      Opti opti = new Opti();
      MX xo = opti.variable();
      opti.minimize(xo.minus(1.0).power(4.0));
      opti.subject_to(ca.le(-5.0, xo.le(5.0)));
      Map<String, Object> opts = new LinkedHashMap<>(s.opts);
      if (s.name.equals("bonmin")) opts.put("discrete", new boolean[] {true});
      opti.solver(s.name, opts);
      opti.solve();
      Raiser r = new Raiser();
      keep.add(r);
      opti.callback_class(r);
      assertThrows(RuntimeException.class, opti::solve, s.toString());

      opts.put("iteration_callback_ignore_errors", false);
      opti.solver(s.name, opts);
      assertThrows(RuntimeException.class, opti::solve, s.toString());

      opts.put("iteration_callback_ignore_errors", true);
      opti.solver(s.name, opts);
      opti.solve();
      // final KeyboardInterrupt-overrides-ignore_errors check: python-only (KeyboardInterrupt semantics)
    }
    keep.add(interrupt);
  }

  // test_eval_interrupt_stops_ipopt: python-only (KeyboardInterrupt in a Callback is special-cased to abort the solve)

  @Test
  public void test_discrete() {
    message("rosenbrock, limited-memory hessian approx");
    SX x = SX.sym("x"), y = SX.sym("y");
    Map<String, SX> nlp = sxnlp("x", vc(x, y), "f", rosen(1.4, x, y));
    for (S s : solvers) {
      if (!s.discrete) continue;
      Function solver = ca.nlpsol("mysolver", s.name, nlp, with(s.opts, "discrete", new boolean[] {true, false}));
      Map<String, DM> solver_in = in("lbx", rep(-10, 2), "ubx", rep(10, 2));
      Map<String, DM> solver_out = solver.call(solver_in);
      assertTrue(success(solver));
      almost(el(solver_out.get("f"), 0), 0.16, 10, s.name);
      almost(el(solver_out.get("x"), 0), 1, 6, s.name);
      almost(el(solver_out.get("x"), 1), 1, 6, s.name);
      checkSerialize(solver, solver_in);
    }
    for (S s : solvers) {
      if (!s.discrete) continue;
      Function solver = ca.nlpsol("mysolver", s.name, nlp, with(s.opts, "discrete", new boolean[] {false, true}));
      Map<String, DM> solver_in = in("lbx", rep(-10, 2), "ubx", rep(10, 2));
      Map<String, DM> solver_out = solver.call(solver_in);
      assertTrue(success(solver));
      almost(el(solver_out.get("f"), 0), 2e-4, 3, s.name);
      almost(el(solver_out.get("x"), 0), Math.sqrt(2), 3, s.name);
      almost(el(solver_out.get("x"), 1), 2, 3, s.name);
      checkSerialize(solver, solver_in);
    }
  }

  @Test
  public void test_nan() {
    SX x = SX.sym("x");
    Map<String, SX> nlp = sxnlp("x", x, "f", x.unaryMinus(), "g", x);
    for (S s : solvers) {
      message("test_nan " + s);
      Function solver = ca.nlpsol("mysolver", s.name, nlp, s.opts);
      for (String v : List.of("x", "g")) {
        String lb = "lb" + v, ub = "ub" + v;
        for (Map<String, DM> data : List.of(in(lb, 3, ub, -3), in(lb, inf, ub, inf), in(lb, -inf, ub, -inf),
            in(lb, Double.NaN), in(ub, Double.NaN))) {
          assertInException("Ill-posed", () -> solver.call(data));
        }
      }
    }
  }

  @Test
  public void test_wrongdims() {
    SX x = SX.sym("x", 2);
    Map<String, SX> nlp = sxnlp("x", x, "f", x.get(0).unaryMinus(), "g", ca.diag(x));
    for (S s : solvers) {
      Map<String, SX> n = nlp;
      assertInException("dense vector", () -> ca.nlpsol("mysolver", s.name, n, s.opts));
    }
    nlp = sxnlp("x", x, "f", x.get(0).unaryMinus(), "g", ca.mtimes(x, x.T()));
    for (S s : solvers) {
      Map<String, SX> n = nlp;
      assertInException("dense vector", () -> ca.nlpsol("mysolver", s.name, n, s.opts));
    }
    nlp = sxnlp("x", x, "f", new SX(1, 1), "g", x);
    for (S s : solvers) {
      Map<String, SX> n = nlp;
      assertInException("dense", () -> ca.nlpsol("mysolver", s.name, n, s.opts));
    }
    nlp = sxnlp("x", x, "f", SX.zeros(0, 0), "g", x);
    for (S s : solvers) ca.nlpsol("mysolver", s.name, nlp, s.opts);
    nlp = sxnlp("x", x, "g", x);
    for (S s : solvers) ca.nlpsol("mysolver", s.name, nlp, s.opts);
    SX xv = ca.vec(ca.diag(SX.sym("x", 2)));
    nlp = sxnlp("x", xv, "f", ca.mtimes(xv.T(), xv), "g", xv.get(0));
    for (S s : solvers) {
      Map<String, SX> n = nlp;
      assertInException("dense vector", () -> ca.nlpsol("mysolver", s.name, n, s.opts));
    }
  }

  @Test
  public void test_initialcond() {
    SX x = SX.sym("x");
    Map<String, SX> nlp = sxnlp("x", x, "f", ca.cos(x).unaryMinus(), "g", x);
    for (S s : solvers) {
      if (s.name.equals("fatrop")) continue;
      Function solver = ca.nlpsol("mysolver", s.name, nlp, s.opts);
      Map<String, DM> solver_out = solver.call(in("x0", vc(6 * pi + 0.01), "lbx", -inf, "ubx", inf, "lbg", -100, "ubg", 100));
      assertTrue(success(solver), s.toString());
      almost(el(solver_out.get("x"), 0), 6 * pi, 6, s.name);
    }
  }

  @Test
  public void test_boundsviol() {
    SX x = SX.sym("x");
    Map<String, SX> nlp = sxnlp("x", x, "f", sq(x.minus(1.0)), "g", x);
    for (S s : solvers) {
      Function solver = ca.nlpsol("mysolver", s.name, nlp, s.opts);
      assertThrows(RuntimeException.class, () -> solver.call(in("lbx", vc(-10), "ubx", vc(-20), "lbg", vc(-10), "ubg", vc(10))));
    }
    for (S s : solvers) {
      Function solver = ca.nlpsol("mysolver", s.name, nlp, s.opts);
      assertThrows(RuntimeException.class, () -> solver.call(in("lbx", vc(-10), "ubx", vc(10), "lbg", vc(-10), "ubg", vc(-20))));
    }
  }

  @Test
  public void test_IPOPT() {
    SX x = SX.sym("x");
    Map<String, SX> nlp = sxnlp("x", x, "f", sq(x.minus(1.0)), "g", x);
    for (S s : solvers) {
      message("test_IPOPT " + s);
      Function solver = ca.nlpsol("mysolver", s.name, nlp, s.opts);
      Map<String, DM> solver_in = in("lbx", vc(-10), "ubx", vc(10), "lbg", vc(-10), "ubg", vc(10));
      Map<String, DM> solver_out = solver.call(solver_in);
      assertTrue(success(solver));
      almost(el(solver_out.get("f"), 0), 0, 10, s.name);
      almost(el(solver_out.get("x"), 0), 1, 9, s.name);
      if (!s.name.contains("bonmin")) almost(el(solver_out.get("g"), 0), 1, 9, s.name);
      if (!s.name.contains("bonmin")) almost(el(solver_out.get("lam_x"), 0), 0, 8, s.name);
      if (!s.name.contains("bonmin")) almost(el(solver_out.get("lam_g"), 0), 0, 9, s.name);
      assertEquals("SOLVER_RET_SUCCESS", stat(solver, "unified_return_status"));
      if (s.codegen != null) checkCodegen(solver, solver_in, s.codegen, 15, true, false, null);
    }
  }

  @Test
  public void test_IPOPT_par() {
    SX x = SX.sym("x"), p = SX.sym("p");
    Map<String, SX> nlp = sxnlp("x", x, "p", p, "f", sq(x.minus(p)), "g", x);
    for (S s : solvers) {
      Function solver = ca.nlpsol("mysolver", s.name, nlp, s.opts);
      Map<String, DM> solver_in = in("lbx", vc(-10), "ubx", vc(10), "lbg", vc(-10), "ubg", vc(10), "p", 1);
      Map<String, DM> solver_out = solver.call(solver_in);
      assertTrue(success(solver));
      almost(el(solver_out.get("f"), 0), 0, 10, s.name);
      almost(el(solver_out.get("x"), 0), 1, 9, s.name);
      if (!s.name.contains("bonmin")) almost(el(solver_out.get("lam_x"), 0), 0, 8, s.name);
      if (!s.name.contains("bonmin")) almost(el(solver_out.get("lam_g"), 0), 0, 9, s.name);
      if (s.codegen != null) checkCodegen(solver, solver_in, s.codegen);
    }
  }

  @Test
  public void test_IPOPTinf() {
    message("trivial IPOPT, infinity bounds");
    SX x = SX.sym("x");
    Map<String, SX> nlp = sxnlp("x", x, "f", sq(x.minus(1.0)), "g", x);
    for (S s : solvers) {
      Function solver = ca.nlpsol("mysolver", s.name, nlp, s.opts);
      Map<String, DM> solver_in = in("lbx", vc(-inf), "ubx", vc(inf), "lbg", vc(-inf), "ubg", vc(inf));
      if (s.name.equals("worhp")) {
        assertThrows(RuntimeException.class, () -> solver.call(solver_in));
        return;
      }
      Map<String, DM> solver_out = solver.call(solver_in);
      assertTrue(success(solver));
      almost(el(solver_out.get("f"), 0), 0, 10, s.name);
      almost(el(solver_out.get("x"), 0), 1, 7, s.name);
      if (!s.name.contains("bonmin")) almost(el(solver_out.get("lam_x"), 0), 0, 9, s.name);
      if (!s.name.contains("bonmin")) almost(el(solver_out.get("lam_g"), 0), 0, 9, s.name);
      if (s.codegen != null) checkCodegen(solver, solver_in, s.codegen);

      Function solver2 = ca.nlpsol("mysolver", s.name, nlp, s.opts);
      Map<String, DM> solver_in2 = in();
      solver_out = solver2.call(solver_in2);
      assertTrue(success(solver2));
      almost(el(solver_out.get("f"), 0), 0, 10, s.name);
      almost(el(solver_out.get("x"), 0), 1, 7, s.name);
      if (!s.name.contains("bonmin")) almost(el(solver_out.get("lam_x"), 0), 0, 9, s.name);
      if (!s.name.contains("bonmin")) almost(el(solver_out.get("lam_g"), 0), 0, 9, s.name);
      if (s.codegen != null) checkCodegen(solver2, solver_in2, s.codegen);
    }
  }

  @Test
  public void test_IPOPTrb() {
    message("rosenbrock, limited-memory hessian approx");
    SX x = SX.sym("x"), y = SX.sym("y");
    Map<String, SX> nlp = sxnlp("x", vc(x, y), "f", rosen(1, x, y));
    for (S s : solvers) {
      if (s.name.equals("sqpmethod") && s.in("regularize")) continue;
      if (s.name.equals("snopt")) continue;
      Function solver = ca.nlpsol("mysolver", s.name, nlp, s.opts);
      Map<String, DM> solver_in = in("lbx", rep(-10, 2), "ubx", rep(10, 2));
      Map<String, DM> solver_out = solver.call(solver_in);
      assertTrue(success(solver));
      almost(el(solver_out.get("f"), 0), 0, 10, s.name);
      almost(el(solver_out.get("x"), 0), 1, 6, s.name);
      almost(el(solver_out.get("x"), 1), 1, 6, s.name);
      if (!s.name.contains("bonmin")) almost(el(solver_out.get("lam_x"), 0), 0, 5, s.name);
      if (!s.name.contains("bonmin")) almost(el(solver_out.get("lam_x"), 1), 0, 5, s.name);
      if (s.codegen != null) checkCodegen(solver, solver_in, s.codegen);
    }
  }

  @Test
  public void test_IPOPTrb2() {
    message("rosenbrock, limited-memory hessian approx");
    SX x = SX.sym("x"), y = SX.sym("y");
    Map<String, SX> nlp = sxnlp("x", vc(x, y), "f", rosen(1, x, y), "g", x.plus(y));
    for (S s : solvers) {
      if (s.name.equals("sqpmethod") && s.in("regularize")) continue;
      Function solver = ca.nlpsol("mysolver", s.name, nlp, s.opts);
      Map<String, DM> solver_in = in("lbx", rep(-10, 2), "ubx", rep(10, 2), "lbg", vc(-10), "ubg", vc(10));
      Map<String, DM> solver_out = solver.call(solver_in);
      if (!s.name.contains("knitro")) assertTrue(success(solver));
      int digits = 6;
      almost(el(solver_out.get("f"), 0), 0, digits, s.name);
      almost(el(solver_out.get("x"), 0), 1, digits, s.name);
      almost(el(solver_out.get("x"), 1), 1, digits, s.name);
      if (!s.name.contains("bonmin")) almost(el(solver_out.get("lam_x"), 0), 0, 5, s.name);
      if (!s.name.contains("bonmin")) almost(el(solver_out.get("lam_x"), 1), 0, 5, s.name);
      if (!s.name.contains("bonmin")) almost(el(solver_out.get("lam_g"), 0), 0, 5, s.name);
      if (s.codegen != null) checkCodegen(solver, solver_in, s.codegen);
    }
  }

  @Test
  public void test_IPOPTrbf() {
    message("rosenbrock fixed, limited-memory hessian approx");
    SX x = SX.sym("x"), y = SX.sym("y");
    Map<String, SX> nlp = sxnlp("x", vc(x, y), "f", rosen(1, x, y), "g", x.plus(y));
    for (S s : solvers) {
      Function solver = ca.nlpsol("mysolver", s.name, nlp, s.opts);
      Map<String, DM> solver_in = in("x0", vc(0, 1), "lbx", vc(-10, 1), "ubx", vc(10, 1), "lbg", vc(-10), "ubg", vc(10));
      if (s.in("qrqp")) solver_in.put("x0", vc(0.6, 1));
      Map<String, DM> solver_out = solver.call(solver_in);
      assertTrue(success(solver));
      almost(el(solver_out.get("f"), 0), 0, 10, s.name);
      almost(el(solver_out.get("x"), 0), 1, 7, s.name);
      almost(el(solver_out.get("x"), 1), 1, 7, s.name);
      if (!s.name.contains("stabilizedsqp")) {
        if (!s.name.contains("bonmin")) almost(el(solver_out.get("lam_x"), 0), 0, 5, s.name);
        if (!s.name.contains("bonmin")) almost(el(solver_out.get("lam_x"), 1), 0, 5, s.name);
        if (!s.name.contains("bonmin")) almost(el(solver_out.get("lam_g"), 0), 0, 5, s.name);
      }
      if (s.codegen != null) checkCodegen(solver, solver_in, s.codegen);
    }
  }

  @Test
  public void test_knitro_options_file() throws IOException {
    Assumptions.assumeTrue(hasNlpsol("knitro"));
    MX x = MX.sym("x");
    Map<String, MX> nlp = mxnlp("x", x, "f", sq(x));
    java.nio.file.Files.writeString(java.nio.file.Path.of("my_knitro.opt"), "maxit            42\n");
    Function solver = ca.nlpsol("solver", "knitro", nlp, d("options_file", "my_knitro.opt"));
    assertOutput(List.of("maxit:", "42"), List.of(), () -> solver.call(in()));
  }

  @Test
  public void test_ipopt_linear_solvers() {
    Assumptions.assumeTrue(hasNlpsol("ipopt"));
    List<String> lsolvers = new ArrayList<>(List.of("mumps"));
    if (CasadiMeta.feature_list().contains("spral")) lsolvers.add("spral");
    if (hasLinsol("ma27") && !System.getProperty("os.name").toLowerCase().contains("mac")) lsolvers.add("ma27");
    SX x = SX.sym("x"), y = SX.sym("y");
    Map<String, SX> nlp = sxnlp("x", vc(x, y), "f", rosen(1, x, y), "g", sq(x).plus(sq(y)));
    double c_r = 4.56748075136258e-02;
    double[] x_r = {7.86415156987791e-01, 6.17698316967954e-01};
    for (String ls : lsolvers) {
      Function solver = ca.nlpsol("mysolver", "ipopt", nlp, d("ipopt.linear_solver", ls));
      Map<String, DM> solver_out = solver.call(in("x0", vc(0.5, 0.5), "lbx", rep(-10, 2), "ubx", rep(10, 2), "lbg", vc(0), "ubg", vc(1)));
      int digits = 7;
      almost(el(solver_out.get("f"), 0), c_r, digits, ls);
      almost(el(solver_out.get("x"), 0), x_r[0], digits, ls);
      almost(el(solver_out.get("x"), 1), x_r[1], digits, ls);
    }
  }

  @Test
  public void test_uno_linear_solvers() {
    Assumptions.assumeTrue(hasNlpsol("uno"));
    List<String> lsolvers = new ArrayList<>(List.of("MUMPS"));
    if (CasadiMeta.feature_list().contains("spral")) lsolvers.add("SSIDS");
    if (hasLinsol("ma27") && !System.getProperty("os.name").toLowerCase().contains("mac")) lsolvers.add("MA27");
    SX x = SX.sym("x"), y = SX.sym("y");
    Map<String, SX> nlp = sxnlp("x", vc(x, y), "f", rosen(1, x, y), "g", sq(x).plus(sq(y)));
    double c_r = 4.56748075136258e-02;
    double[] x_r = {7.86415156987791e-01, 6.17698316967954e-01};
    for (String ls : lsolvers) {
      Function mysolver = ca.nlpsol("mysolver", "uno", nlp, d("print_time", false, "uno", d("preset", "ipopt",
          "linear_solver", ls, "primal_tolerance", 1e-8, "dual_tolerance", 1e-8)));
      Map<String, DM> solver_out = mysolver.call(in("x0", vc(0.5, 0.5), "lbx", rep(-10, 2), "ubx", rep(10, 2), "lbg", vc(0), "ubg", vc(1)));
      int digits = 6;
      almost(el(solver_out.get("f"), 0), c_r, digits, ls);
      almost(el(solver_out.get("x"), 0), x_r[0], digits, ls);
      almost(el(solver_out.get("x"), 1), x_r[1], digits, ls);
    }
  }

  @Test
  public void test_ipopt_stats_per_call() {
    Assumptions.assumeTrue(hasNlpsol("ipopt"));
    MX x = MX.sym("x");
    Function solver = ca.nlpsol("solver", "ipopt", mxnlp("x", x, "f", sq(x)));
    solver.call(in());
    Map<String, Object> s1 = solver.stats();
    solver.call(in());
    Map<String, Object> s2 = solver.stats();
    assertEquals(s1.get("n_call_nlp_f"), s2.get("n_call_nlp_f"));
  }

  @Test
  public void test_warmstart() {
    SX x = SX.sym("x"), y = SX.sym("y");
    Map<String, SX> nlp = sxnlp("x", vc(x, y), "f", rosen(1, x, y), "g", sq(x).plus(sq(y)));
    double c_r = 4.56748075136258e-02;
    double[] x_r = {7.86415156987791e-01, 6.17698316967954e-01};
    for (S s : solvers) {
      Function solver = ca.nlpsol("mysolver", s.name, nlp, s.opts);
      Map<String, DM> solver_in = in("x0", vc(0.5, 0.5), "lbx", rep(-10, 2), "ubx", rep(10, 2), "lbg", vc(0), "ubg", vc(1));
      Map<String, DM> solver_out = solver.call(solver_in);
      Map<String, DM> oldsolver_out = solver_out;
      int digits = 5;
      almost(el(solver_out.get("f"), 0), c_r, digits, s.name);
      almost(el(solver_out.get("x"), 0), x_r[0], digits, s.name);
      almost(el(solver_out.get("x"), 1), x_r[1], digits, s.name);
      if (!s.name.contains("bonmin")) almost(el(solver_out.get("lam_x"), 0), 0, s.name.equals("snopt") ? 5 : 8, s.name);
      if (!s.name.contains("bonmin")) almost(el(solver_out.get("lam_x"), 1), 0, s.name.equals("snopt") ? 5 : 8, s.name);
      if (!List.of("bonmin", "sleqp").contains(s.name)) almost(el(solver_out.get("lam_g"), 0), 0.12149655447670, 6, s.name);
      if (s.codegen != null) checkCodegen(solver, solver_in, s.codegen, codegen_check_digits);

      message(":warmstart");
      if (s.name.contains("ipopt")) {
        Map<String, Object> options = with(s.opts, "ipopt.warm_start_init_point", "yes", "ipopt.warm_start_bound_push", 1e-6,
            "ipopt.warm_start_slack_bound_push", 1e-6, "ipopt.warm_start_mult_bound_push", 1e-6, "ipopt.mu_init", 1e-6);
        solver = ca.nlpsol("mysolver", s.name, nlp, options);
        solver_in.put("lbx", vc(rep(-10, 2)));
        solver_in.put("ubx", vc(rep(10, 2)));
        solver_in.put("lbg", vc(0));
        solver_in.put("ubg", vc(1));
        solver_in.put("x0", oldsolver_out.get("x"));
        solver_in.put("lam_g0", oldsolver_out.get("lam_g"));
        if (!s.name.contains("bonmin")) solver_in.put("lam_x0", oldsolver_out.get("lam_x"));
        solver.call(solver_in);
        if (s.codegen != null) checkCodegen(solver, solver_in, s.codegen, codegen_check_digits);
      }
    }
  }

  @Test
  public void test_IPOPTrhb2_gen() {
    message("rosenbrock, exact hessian generated, constrained");
    SX x = SX.sym("x"), y = SX.sym("y");
    Map<String, SX> nlp = sxnlp("x", vc(x, y), "f", rosen(1, x, y), "g", sq(x).plus(sq(y)));
    double c_r = 4.56748075136258e-02;
    double[] x_r = {7.86415156987791e-01, 6.17698316967954e-01};
    for (S s : solvers) {
      if (s.name.contains("madnlp")) continue;
      Function solver = ca.nlpsol("mysolver", s.name, nlp, s.opts);
      Map<String, DM> solver_in = in("x0", vc(0.5, 0.5), "lbx", rep(-10, 2), "ubx", rep(10, 2), "lbg", vc(0), "ubg", vc(1));
      Map<String, DM> solver_out = solver.call(solver_in);
      int digits = 5;
      almost(el(solver_out.get("f"), 0), c_r, digits, s.name + el(solver_out.get("f"), 0) + ":" + c_r);
      almost(el(solver_out.get("x"), 0), x_r[0], digits, s.name);
      almost(el(solver_out.get("x"), 1), x_r[1], digits, s.name);
      if (!s.name.contains("bonmin")) almost(el(solver_out.get("lam_x"), 0), 0, s.name.equals("snopt") ? 5 : 8, s.name);
      if (!s.name.contains("bonmin")) almost(el(solver_out.get("lam_x"), 1), 0, s.name.equals("snopt") ? 5 : 8, s.name);
      if (!List.of("bonmin", "sleqp").contains(s.name)) almost(el(solver_out.get("lam_g"), 0), 0.12149655447670, 6, s.name);
      checkSerialize(solver, solver_in);
      if (s.codegen != null) checkCodegen(solver, solver_in, s.codegen, codegen_check_digits);
    }
  }

  @Test
  public void test_jacG_empty() {
    SX x = SX.sym("x"), y = SX.sym("y");
    Map<String, SX> nlp = sxnlp("x", vc(x, y), "f", rosen(1, x, y), "g", 1);
    for (S s : solvers) {
      if (List.of("worhp", "sqpmethod", "snopt", "fatrop", "madnlp").contains(s.name)) continue;
      Function solver = ca.nlpsol("mysolver", s.name, nlp, s.opts);
      Map<String, DM> solver_in = in("x0", vc(0.5, 0.5), "lbx", rep(-10, 2), "ubx", rep(10, 2), "lbg", vc(0), "ubg", vc(2));
      Map<String, DM> solver_out = solver.call(solver_in);
      int digits = 5;
      chk(solver_out.get("f"), vc(0), digits, s.name);
      chk(solver_out.get("x"), vc(1, 1), digits, s.name);
      if (!s.name.contains("bonmin")) chk(solver_out.get("lam_x"), vc(0, 0), digits, s.name);
      if (!s.name.contains("bonmin")) chk(solver_out.get("lam_g"), vc(0), digits, s.name);
      if (s.codegen != null) {
        if (s.name.equals("ipopt")) continue;
        checkCodegen(solver, solver_in, s.codegen);
      }
    }
  }

  @Test
  public void test_IPOPTrhb2_gen_par() {
    message("rosenbrock, exact hessian generated, constrained, parametric");
    SX x = SX.sym("x"), y = SX.sym("y"), p = SX.sym("p");
    Map<String, SX> nlp = sxnlp("x", vc(x, y), "p", p, "f", rosen(p, x, y), "g", sq(x).plus(sq(y)));
    double c_r = 4.56748075136258e-02;
    double[] x_r = {7.86415156987791e-01, 6.17698316967954e-01};
    for (S s : solvers) {
      Function solver = ca.nlpsol("mysolver", s.name, nlp, s.opts);
      Map<String, DM> solver_in = in("x0", vc(0.5, 0.5), "lbx", rep(-10, 2), "ubx", rep(10, 2), "lbg", vc(0), "ubg", vc(1), "p", vc(1));
      if (s.codegen != null) checkCodegen(solver, solver_in, s.codegen, codegen_check_digits);
      Map<String, DM> solver_out = solver.call(solver_in);
      int digits = 5;
      almost(el(solver_out.get("f"), 0), c_r, digits, s.name);
      almost(el(solver_out.get("x"), 0), x_r[0], digits, s.name);
      almost(el(solver_out.get("x"), 1), x_r[1], digits, s.name);
      if (!s.name.contains("bonmin")) almost(el(solver_out.get("lam_x"), 0), 0, s.name.equals("snopt") ? 5 : 8, s.name);
      if (!s.name.contains("bonmin")) almost(el(solver_out.get("lam_x"), 1), 0, s.name.equals("snopt") ? 5 : 8, s.name);
      if (!List.of("bonmin", "sleqp").contains(s.name)) almost(el(solver_out.get("lam_g"), 0), 0.12149655447670, 6, s.name);
    }
  }

  @Test
  public void test_IPOPTrhb_gen() {
    message("rosenbrock, exact hessian generated");
    SX x = SX.sym("x"), y = SX.sym("y");
    Map<String, SX> nlp = sxnlp("x", vc(x, y), "f", rosen(1, x, y));
    for (S s : solvers) {
      if (s.name.equals("sqpmethod") && s.in("regularize")) continue;
      if (s.name.equals("snopt")) continue;
      Function solver = ca.nlpsol("mysolver", s.name, nlp, s.opts);
      Map<String, DM> solver_in = in("lbx", rep(-10, 2), "ubx", rep(10, 2));
      Map<String, DM> solver_out = solver.call(solver_in);
      int xdig = 6;
      if (s.name.equals("uno") && s.in("LBFGS")) xdig = 7;
      almost(el(solver_out.get("f"), 0), 0, 10, s.name);
      almost(el(solver_out.get("x"), 0), 1, xdig, s.name);
      almost(el(solver_out.get("x"), 1), 1, xdig, s.name);
      if (!s.name.contains("bonmin")) almost(el(solver_out.get("lam_x"), 0), 0, 8, s.name);
      if (!s.name.contains("bonmin")) almost(el(solver_out.get("lam_x"), 1), 0, 8, s.name);
      if (s.codegen != null) checkCodegen(solver, solver_in, s.codegen);
    }
  }

  @Test
  public void test_IPOPTrhb_gen_xnonfree() {
    message("rosenbrock, exact hessian generated, non-free x");
    SX x = SX.sym("x"), y = SX.sym("y");
    Map<String, SX> nlp = sxnlp("x", vc(x, y), "f", rosen(1, x, y));
    for (S s : solvers) {
      if (s.name.equals("sqpmethod") && s.in("regularize")) continue;
      if (s.name.equals("snopt")) continue;
      Function solver = ca.nlpsol("mysolver", s.name, nlp, s.opts);
      Map<String, DM> solver_in = in("lbx", vc(1, -10), "ubx", vc(1, 10));
      Map<String, DM> solver_out = solver.call(solver_in);
      almost(el(solver_out.get("f"), 0), 0, 10, s.name);
      almost(el(solver_out.get("x"), 0), 1, 9, s.name);
      almost(el(solver_out.get("x"), 1), 1, 6, s.name);
      if (!s.name.contains("bonmin")) almost(el(solver_out.get("lam_x"), 0), 0, 6, s.name);
      if (!s.name.contains("bonmin")) almost(el(solver_out.get("lam_x"), 1), 0, 6, s.name);
      if (s.codegen != null) checkCodegen(solver, solver_in, s.codegen, codegen_check_digits);
    }
  }

  @Test
  public void test_IPOPTrhb_gen_par() {
    message("rosenbrock, exact hessian generated, parametric");
    SX x = SX.sym("x"), y = SX.sym("y"), p = SX.sym("p");
    Map<String, SX> nlp = sxnlp("x", vc(x, y), "p", p, "f", rosen(p, x, y));
    for (S s : solvers) {
      if (s.name.equals("sqpmethod") && s.in("regularize")) continue;
      if (s.name.equals("snopt")) continue;
      Function solver = ca.nlpsol("mysolver", s.name, nlp, s.opts);
      Map<String, DM> solver_in = in("lbx", rep(-10, 2), "ubx", rep(10, 2), "p", 1);
      Map<String, DM> solver_out = solver.call(solver_in);
      int xdig = 6;
      if (s.name.equals("uno") && s.in("LBFGS")) xdig = 7;
      almost(el(solver_out.get("f"), 0), 0, 10, s.name);
      almost(el(solver_out.get("x"), 0), 1, xdig, s.name);
      almost(el(solver_out.get("x"), 1), 1, xdig, s.name);
      if (s.codegen != null) checkCodegen(solver, solver_in, s.codegen);
    }
  }

  @Test
  public void test_IPOPTnorm() {
    message("IPOPT min ||x||^2_2");
    int N = 10;
    MX x = MX.sym("x", N);
    DM x0 = ca.linspace(new DM(0.0), new DM(1.0), N);
    MX X0 = new MX(x0);
    MX d = x.minus(X0);
    Map<String, MX> nlp = mxnlp("x", x, "f", ca.dot(d, d), "g", ca.times(2.0, x));
    for (S s : solvers) {
      Function solver = ca.nlpsol("mysolver", s.name, nlp, s.opts);
      Map<String, DM> solver_in = in("lbx", rep(-10, N), "ubx", rep(10, N), "lbg", rep(-10, N), "ubg", rep(10, N));
      Map<String, DM> solver_out = solver.call(solver_in);
      almost(el(solver_out.get("f"), 0), 0, 10, s.name);
      chk(solver_out.get("x"), x0, 8, s.name);
      if (!s.name.contains("bonmin")) chk(solver_out.get("lam_x"), DM.zeros(10, 1), 8, s.name);
      if (!s.name.contains("bonmin")) almost(el(solver_out.get("lam_g"), 1), 0, 8, s.name);
      if (s.codegen != null) checkCodegen(solver, solver_in, s.codegen);
    }
  }

  @Test
  public void test_IPOPTnoc() {
    message("trivial IPOPT, no constraints");
    SX x = SX.sym("x");
    Map<String, SX> nlp = sxnlp("x", x, "f", sq(x.minus(1.0)));
    for (S s : solvers) {
      if (s.name.equals("snopt")) continue;
      Function solver = ca.nlpsol("mysolver", s.name, nlp, s.opts);
      Map<String, DM> solver_in = in("lbx", vc(-10), "ubx", vc(10));
      Map<String, DM> solver_out = solver.call(solver_in);
      almost(el(solver_out.get("f"), 0), 0, 10, s.name);
      almost(el(solver_out.get("x"), 0), 1, 7, s.name);
      if (s.codegen != null) checkCodegen(solver, solver_in, s.codegen);
    }
  }

  @Test
  public void test_IPOPTmx() {
    message("trivial IPOPT, using MX");
    MX x = MX.sym("x");
    Map<String, MX> nlp = mxnlp("x", x, "f", sq(x.minus(1.0)), "g", ca.times(2.0, x));
    for (S s : solvers) {
      Function solver = ca.nlpsol("mysolver", s.name, nlp, s.opts);
      Map<String, DM> solver_in = in("lbx", vc(-10), "ubx", vc(10), "lbg", vc(-10), "ubg", vc(10));
      Map<String, DM> solver_out = solver.call(solver_in);
      almost(el(solver_out.get("f"), 0), 0, 10, s.name);
      almost(el(solver_out.get("x"), 0), 1, 9, s.name);
      if (s.codegen != null) checkCodegen(solver, solver_in, s.codegen);
    }
  }

  @Test
  public void test_IPOPTc() {
    message("trivial, overconstrained");
    SX x = SX.sym("x");
    Map<String, SX> nlp = sxnlp("x", x, "f", sq(x.minus(1.0)), "g", vc(x, x, x));
    for (S s : solvers) {
      Function solver = ca.nlpsol("mysolver", s.name, nlp, s.opts);
      Map<String, DM> solver_in = in("lbx", vc(-10), "ubx", vc(10), "lbg", vc(-10, -10, -10), "ubg", vc(10, 10, 10));
      Map<String, DM> solver_out = solver.call(solver_in);
      almost(el(solver_out.get("f"), 0), 0, 9, s.name);
      almost(el(solver_out.get("x"), 0), 1, 5, s.name);
      if (s.codegen != null) checkCodegen(solver, solver_in, s.codegen);
    }
  }

  @Test
  public void test_IPOPTc2() {
    message("trivial2, overconstrained");
    SX x = SX.sym("x");
    Map<String, SX> nlp = sxnlp("x", x, "f", sq(x.minus(1.0)), "g", vc(x, x, x.plus(x)));
    for (S s : solvers) {
      Function solver = ca.nlpsol("mysolver", s.name, nlp, s.opts);
      Map<String, DM> solver_in = in("lbx", vc(-10), "ubx", vc(10), "lbg", vc(-10, -10, -10), "ubg", vc(10, 10, 10));
      Map<String, DM> solver_out = solver.call(solver_in);
      almost(el(solver_out.get("f"), 0), 0, 10, s.name);
      almost(el(solver_out.get("x"), 0), 1, 8, s.name);
      if (s.codegen != null) checkCodegen(solver, solver_in, s.codegen);
    }
  }

  @Test
  public void test_IPOPTcmx() {
    message("trivial , overconstrained, using MX");
    MX x = MX.sym("x");
    Map<String, MX> nlp = mxnlp("x", x, "f", sq(x.minus(1.0)), "g", vc(ca.times(2.0, x), ca.times(3.0, x), ca.times(4.0, x)));
    for (S s : solvers) {
      Function solver = ca.nlpsol("mysolver", s.name, nlp, s.opts);
      Map<String, DM> solver_in = in("lbx", vc(-10), "ubx", vc(10), "lbg", vc(-10, -10, -10), "ubg", vc(10, 10, 10));
      Map<String, DM> solver_out = solver.call(solver_in);
      almost(el(solver_out.get("f"), 0), 0, 9, s.name);
      almost(el(solver_out.get("x"), 0), 1, 8, s.name);
      if (s.codegen != null) checkCodegen(solver, solver_in, s.codegen);
    }
  }

  @Test
  public void test_IPOPTdeg() {
    message("degenerate optimization IPOPT");
    SX x = SX.sym("x"), y = SX.sym("y");
    Map<String, SX> nlp = sxnlp("x", vc(x, y), "f", 0, "g", vc(x.minus(y), x));
    for (S s : solvers) {
      Function solver = ca.nlpsol("mysolver", s.name, nlp, s.opts);
      Map<String, DM> solver_in = in("lbx", vc(-10, -10), "ubx", vc(10, 10), "lbg", vc(0, 3), "ubg", vc(0, 3));
      Map<String, DM> solver_out = solver.call(solver_in);
      almost(el(solver_out.get("x"), 0), el(solver_out.get("x"), 1), s.in("sqic") ? 4 : 10, "IPOPT");
      if (s.codegen != null) checkCodegen(solver, solver_in, s.codegen);
    }
  }

  @Test
  public void test_IPOPTdegc() {
    message("degenerate optimization IPOPT, overconstrained");
    SX x = SX.sym("x"), y = SX.sym("y");
    Map<String, SX> nlp = sxnlp("x", vc(x, y), "f", 0, "g", vc(x.minus(y), x, x.plus(y)));
    for (S s : solvers) {
      Function solver = ca.nlpsol("mysolver", s.name, nlp, s.opts);
      Map<String, DM> solver_in = in("lbx", vc(-10, -10), "ubx", vc(10, 10), "lbg", vc(0, 3, -10), "ubg", vc(0, 3, 10));
      Map<String, DM> solver_out = solver.call(solver_in);
      almost(el(solver_out.get("x"), 0), el(solver_out.get("x"), 1), s.in("sqic") ? 4 : 10, "IPOPT");
      if (s.codegen != null) checkCodegen(solver, solver_in, s.codegen);
    }
  }

  @Test
  public void test_XfreeChange() {
    message("Change in X settings");
    SX x = SX.sym("x"), y = SX.sym("y");
    Map<String, SX> nlp = sxnlp("x", vc(x, y), "f", rosen(1, x, y), "g", x.plus(y));
    for (S s : solvers) {
      Function solver = ca.nlpsol("mysolver", s.name, nlp, s.opts);
      Map<String, DM> solver_in = in("x0", vc(0, 1), "lbx", vc(-10, -10), "ubx", vc(10, 10), "lbg", vc(-10), "ubg", vc(10));
      if (s.in("qrqp")) solver_in.put("x0", vc(0.6, 1));
      solver.call(solver_in);
      solver_in.put("lbx", vc(-10, 1));
      solver_in.put("ubx", vc(10, 1));
      solver_in.put("lbg", vc(-10));
      solver_in.put("ubg", vc(10));
      Map<String, DM> solver_out = solver.call(solver_in);
      almost(el(solver_out.get("f"), 0), 0, 9, s.name);
      almost(el(solver_out.get("x"), 0), 1, 6, s.name);
      almost(el(solver_out.get("x"), 1), 1, 6, s.name);
      if (s.codegen != null) checkCodegen(solver, solver_in, s.codegen);
    }
  }

  @Test
  public void test_activeLBX() {
    message("active LBX");
    SX x = SX.sym("x"), y = SX.sym("y");
    Map<String, SX> nlp = sxnlp("x", vc(x, y), "f", rosen(1, x, y), "g", x.plus(y));
    for (S s : solvers) {
      Function solver = ca.nlpsol("mysolver", s.name, nlp, s.opts);
      Map<String, DM> solver_in = in("x0", vc(0, 1), "lbx", vc(-10, 1.2), "ubx", vc(10, 2), "lbg", vc(-10), "ubg", vc(10));
      if (s.in("qrqp")) solver_in.put("x0", vc(0.5, 1));
      Map<String, DM> solver_out = solver.call(solver_in);
      if (el(solver_out.get("x"), 0) < 0) {
        almost(el(solver_out.get("f"), 0), 4.3817250416084308, 6, s.name);
        almost(el(solver_out.get("x"), 0), -1.0910624688699295, 6, s.name);
        almost(el(solver_out.get("x"), 1), 1.2, 5, s.name);
        if (!s.name.contains("bonmin")) almost(el(solver_out.get("lam_x"), 0), 0, s.name.equals("stabilizedsqp") ? 5 : 8, s.toString());
        if (!s.name.contains("bonmin")) almost(el(solver_out.get("lam_x"), 1), -1.9165378046901287, 4, s.name);
        if (!s.name.contains("bonmin")) almost(el(solver_out.get("lam_g"), 0), 0, 8, s.name);
      } else {
        almost(el(solver_out.get("f"), 0), 9.0908263002590e-3, 6, s.name);
        almost(el(solver_out.get("x"), 0), 1.0952466252248, 6, s.name);
        almost(el(solver_out.get("x"), 1), 1.2, 5, s.name);
        if (!s.name.contains("bonmin")) almost(el(solver_out.get("lam_x"), 0), 0, s.name.equals("stabilizedsqp") ? 5 : 8, s.toString());
        if (!s.name.contains("bonmin")) almost(el(solver_out.get("lam_x"), 1), -8.6963632695079e-2, 4, s.name);
        if (!s.name.contains("bonmin")) almost(el(solver_out.get("lam_g"), 0), 0, 8, s.name);
      }
      if (s.codegen != null) checkCodegen(solver, solver_in, s.codegen);
    }
  }

  @Test
  public void test_activeLBG() {
    message("active LBG");
    SX x = SX.sym("x"), y = SX.sym("y");
    Map<String, SX> nlp = sxnlp("x", vc(x, y), "f", rosen(1, x, y), "g", x.plus(y));
    for (S s : solvers) {
      Function solver = ca.nlpsol("mysolver", s.name, nlp, s.opts);
      Map<String, DM> solver_in = in("x0", vc(0, 1), "lbx", vc(-10, -10), "ubx", vc(10, 10), "lbg", vc(2.2), "ubg", vc(10));
      if (s.in("qrqp")) solver_in.put("x0", vc(1.5, 1));
      Map<String, DM> solver_out = solver.call(solver_in);
      almost(el(solver_out.get("f"), 0), 4.252906468284e-3, 6, s.name);
      almost(el(solver_out.get("x"), 0), 1.065181061847138, 6, s.name);
      almost(el(solver_out.get("x"), 1), 1.1348189166291160, 6, s.name);
      if (!s.name.contains("bonmin")) almost(el(solver_out.get("lam_x"), 0), 0, s.name.equals("snopt") ? 6 : 8, s.name);
      if (!s.name.contains("bonmin")) almost(el(solver_out.get("lam_x"), 1), 0, 4, s.name);
      if (!s.name.contains("bonmin")) almost(el(solver_out.get("lam_g"), 0), -4.1644422845712e-2, 3, s.name);
      if (s.codegen != null) checkCodegen(solver, solver_in, s.codegen);
    }
  }

  @Test
  public void test_activeUBG() {
    message("active UBG");
    SX x = SX.sym("x"), y = SX.sym("y");
    Map<String, SX> nlp = sxnlp("x", vc(x, y), "f", rosen(1, x, y), "g", x.plus(y));
    for (S s : solvers) {
      Function solver = ca.nlpsol("mysolver", s.name, nlp, s.opts);
      Map<String, DM> solver_in = in("x0", vc(0, 1), "lbx", vc(-10, -10), "ubx", vc(10, 10), "lbg", vc(0), "ubg", vc(1.8));
      if (s.in("qrqp")) solver_in.put("x0", vc(1.5, 1));
      Map<String, DM> solver_out = solver.call(solver_in);
      almost(el(solver_out.get("f"), 0), 4.64801220074552e-3, 6, s.name);
      almost(el(solver_out.get("x"), 0), 9.318651964592811e-1, 5, s.name);
      almost(el(solver_out.get("x"), 1), 8.68134821123689e-1, 5, s.name);
      if (!s.name.contains("bonmin")) almost(el(solver_out.get("lam_x"), 0), 0, s.name.equals("snopt") ? 5 : 8, s.name);
      if (!s.name.contains("bonmin")) almost(el(solver_out.get("lam_x"), 1), 0, 4, s.name);
      if (!s.name.contains("bonmin")) almost(el(solver_out.get("lam_g"), 0), 4.75846495145007e-2, 5, s.name);
      if (s.codegen != null) checkCodegen(solver, solver_in, s.codegen);
    }
  }

  @Test
  public void test_activeUBX() {
    message("active UBX");
    SX x = SX.sym("x"), y = SX.sym("y");
    Map<String, SX> nlp = sxnlp("x", vc(x, y), "f", rosen(1, x, y), "g", x.plus(y));
    for (S s : solvers) {
      Function solver = ca.nlpsol("mysolver", s.name, nlp, s.opts);
      Map<String, DM> solver_in = in("x0", vc(0, 1), "lbx", vc(-10, 0), "ubx", vc(10, 0.9), "lbg", vc(-10), "ubg", vc(10));
      Map<String, DM> solver_out = solver.call(solver_in);
      almost(el(solver_out.get("f"), 0), 2.626109721583e-3, 6, s.name);
      almost(el(solver_out.get("x"), 0), 9.4882542279172277e-01, 6, s.name);
      almost(el(solver_out.get("x"), 1), 0.9, 6, s.name);
      if (!s.name.contains("bonmin")) almost(el(solver_out.get("lam_x"), 0), 0, s.name.equals("snopt") ? 6 : 8, s.name);
      if (!s.name.contains("bonmin")) almost(el(solver_out.get("lam_x"), 1), 5.39346608659e-2, 4, s.name);
      if (!s.name.contains("bonmin")) almost(el(solver_out.get("lam_g"), 0), 0, 8, s.name);
      if (s.codegen != null) checkCodegen(solver, solver_in, s.codegen);
    }
  }

  @Test
  public void test_QP() {
    message("QP");
    int N = 50;
    SX x = SX.sym("x", N);
    DM x0 = new DM(range(N));
    double[] h = new double[N];
    for (int i = 0; i < N; i++) h[i] = i + 1;
    DM H = ca.diag(new DM(h));
    SX e = x.minus(x0);
    SX obj = ca.times(0.5, ca.mtimes(ca.mtimes(e.T(), H), e));
    Map<String, SX> nlp = sxnlp("x", x, "f", obj);
    for (S s : solvers) {
      if (s.name.equals("snopt") || s.name.equals("fatrop")) continue;
      if (s.name.equals("sqpmethod") && s.in("limited-memory")) continue;
      Function solver = ca.nlpsol("mysolver", s.name, nlp, s.opts);
      Map<String, DM> solver_out = solver.call(in("lbx", -1000, "ubx", 1000));
      chk(solver_out.get("x"), x0, 2, s.name);
      almost(el(solver_out.get("f"), 0), 0, 3, s.name);
      if (!s.name.contains("bonmin")) chk(solver_out.get("lam_x"), DM.zeros(N, 1), 4, s.name);
    }
  }

  @Test
  public void test_QP2() {
    DM H = dm(new double[][] {{1, -1}, {-1, 2}});
    DM G = vc(-2, -6);
    DM A = dm(new double[][] {{1, 1}, {-1, 2}, {2, 1}});
    DM LBA = vc(rep(-inf, 3)), UBA = vc(2, 2, 3);
    DM LBX = vc(0.5, 0), UBX = vc(0.5, inf);
    SX x = SX.sym("x", 2);
    Map<String, SX> nlp = sxnlp("x", x, "f", ca.times(0.5, ca.mtimes(ca.mtimes(x.T(), H), x)).plus(ca.mtimes(G.T(), x)),
        "g", ca.mtimes(A, x));
    for (S s : solvers) {
      Map<String, Object> options = new LinkedHashMap<>(s.opts);
      if (s.name.contains("madnlp")) continue;
      if (s.name.contains("ipopt")) options.put("ipopt.fixed_variable_treatment", "make_constraint");
      for (int rep = 0; rep < 2; rep++) {
        Function solver = ca.nlpsol("mysolver", s.name, nlp, options);
        Map<String, DM> solver_in = in("lbx", LBX, "ubx", UBX, "lbg", LBA, "ubg", UBA);
        if (s.in("sqic") || s.in("daqp")) break;
        Map<String, DM> solver_out = solver.call(solver_in);
        almost(el(solver_out.get("x"), 0), 0.5, 6, s.name);
        almost(el(solver_out.get("x"), 1), 1.25, 6, s.name);
        if (!s.name.contains("bonmin")) almost(el(solver_out.get("lam_x"), 0), 4.75, 6, s.name);
        if (!s.name.contains("bonmin")) almost(el(solver_out.get("lam_x"), 1), 0, 6, s.name);
        if (!s.name.contains("bonmin")) chk(solver_out.get("lam_g"), vc(0, 2, 0), 6, s.name);
        almost(el(solver_out.get("f"), 0), -7.4375, 6, s.name);
        if (s.codegen != null) checkCodegen(solver, solver_in, s.codegen);
      }
    }
  }

  @Test
  public void test_QP2_unconvex() {
    DM H = dm(new double[][] {{1, -1}, {-1, -2}});
    DM G = vc(-2, -6);
    DM A = dm(new double[][] {{1, 1}, {-1, 2}, {2, 1}});
    DM LBA = vc(rep(-inf, 3)), UBA = vc(2, 2, 3);
    DM LBX = vc(0, 0), UBX = vc(inf, inf);
    SX x = SX.sym("x", 2);
    Map<String, SX> nlp = sxnlp("x", x, "f", ca.times(0.5, ca.mtimes(ca.mtimes(x.T(), H), x)).plus(ca.mtimes(G.T(), x)),
        "g", ca.mtimes(A, x));
    for (S s : solvers) {
      Map<String, Object> options = new LinkedHashMap<>(s.opts);
      if (s.name.contains("ipopt")) options.put("ipopt.fixed_variable_treatment", "make_constraint");
      Map<String, DM> solver_in = in("lbx", LBX, "ubx", UBX, "lbg", LBA, "ubg", UBA);
      if (s.in("qrqp")) {
        solver_in.put("x0", vc(1, 1));
        solver_in.put("lam_g0", vc(1, 1, 0));
      }
      if (s.in("daqp")) continue;
      for (int rep = 0; rep < 2; rep++) {
        Function solver = ca.nlpsol("mysolver", s.name, nlp, options);
        String msg = solver.toString();
        Map<String, DM> solver_out = solver.call(solver_in);
        almost(el(solver_out.get("x"), 0), 2.0 / 3, 6, msg);
        almost(el(solver_out.get("x"), 1), 4.0 / 3, 6, msg);
        if (!s.name.contains("bonmin")) almost(el(solver_out.get("lam_x"), 0), 0, 6, msg);
        if (!s.name.contains("bonmin")) almost(el(solver_out.get("lam_x"), 1), 0, 6, msg);
        if (!s.name.contains("bonmin") && !s.name.contains("worhp"))
          chk(solver_out.get("lam_g"), vc(4 + 8.0 / 9, 20.0 / 9, 0), 6, msg);
        almost(el(solver_out.get("f"), 0), -10 - 16.0 / 9, 6, msg);
        if (s.codegen != null) checkCodegen(solver, solver_in, s.codegen);
      }
    }
  }

  @Test
  public void test_bug() {
    MX x = MX.sym("x", 3), y = MX.sym("y", 2);
    Function f = new Function("f", List.of(x, y), List.of(new MX(1.0)));
    MX aa = MX.sym("aa", 5);
    MX a = aa.get(new Slice(0, 3)), b = aa.get(new Slice(3, 5));
    MX f_call = f.call(List.of(a, b)).get(0);
    Map<String, MX> nlp = mxnlp("x", aa, "f", f_call);
    for (S s : solvers) {
      if (s.name.contains("worhp") || s.name.equals("snopt")) continue;
      ca.nlpsol("mysolver", s.name, nlp, s.opts);
    }
  }

  @Test
  public void test_missing_symbols() {
    MX x = MX.sym("x"), p = MX.sym("p");
    for (S s : solvers)
      assertInException("[p] are free", () -> ca.nlpsol("solver", s.name, mxnlp("x", x, "f", sq(x.minus(p))), s.opts));
  }

  @Test
  public void test_no_success() {
    Assumptions.assumeTrue(hasNlpsol("ipopt"));
    SX x = SX.sym("x"), y = SX.sym("y");
    SX f = rosen(1, x, y);
    for (int k = 0; k < solvers.size(); k++) {
      Function solver = ca.nlpsol("solver", "ipopt", sxnlp("x", vc(x, y), "f", f, "g", vc(x.plus(1.0), x.minus(2.0))));
      solver.call(in("x0", 0, "lbg", 0, "ubg", 0));
      assertFalse(success(solver));
      Function solver2 = ca.nlpsol("solver", "ipopt", sxnlp("x", vc(x, y), "f", f, "g", vc(x.plus(1.0), x.minus(2.0))),
          d("error_on_fail", true));
      assertInException("process", () -> solver2.call(in("x0", 0, "lbg", 0, "ubg", 0)));
    }
  }

  private static Function customJacF(MX x) {
    Function J = new Function("jac_f", List.of(x, new MX(1, 1)), List.of(x.unaryMinus()), List.of("x", "out_r"), List.of("jac_r_x"));
    return new Function("f", List.of(x), List.of(sq(x)), List.of("x"), List.of("r"), d("custom_jacobian", J, "jac_penalty", 0));
  }

  @Test
  public void test_postpone_expand() {
    Assumptions.assumeTrue(hasNlpsol("ipopt"));
    MX x = MX.sym("x");
    Function f = customJacF(x);
    Function solver = ca.nlpsol("solver", "ipopt", mxnlp("x", x, "f", f.call(List.of(x.minus(3.0))).get(0)), d("ipopt.max_iter", 5));
    solver.call(in());
    assertEquals("SOLVER_RET_LIMITED", stat(solver, "unified_return_status"));
    assertTrue(solver.get_function("nlp_jac_g").is_a("MXFunction"));

    f = customJacF(x);
    solver = ca.nlpsol("solver", "ipopt", mxnlp("x", x, "f", f.call(List.of(x.minus(3.0))).get(0)), d("ipopt.max_iter", 5, "expand", true));
    solver.call(in());
    assertEquals("Solve_Succeeded", stat(solver, "return_status"));
    assertTrue(solver.get_function("nlp_jac_g").is_a("SXFunction"));

    f = customJacF(x);
    solver = ca.nlpsol("solver", "ipopt", mxnlp("x", x, "f", f.call(List.of(x.minus(3.0))).get(0)), d("ipopt.max_iter", 5, "postpone_expand", true));
    solver.call(in());
    assertEquals("SOLVER_RET_LIMITED", stat(solver, "unified_return_status"));
    assertTrue(solver.get_function("nlp_jac_g").is_a("MXFunction"));

    f = customJacF(x);
    solver = ca.nlpsol("solver", "ipopt", mxnlp("x", x, "f", f.call(List.of(x.minus(3.0))).get(0)),
        d("ipopt.max_iter", 5, "expand", true, "postpone_expand", true));
    solver.call(in());
    assertEquals("SOLVER_RET_LIMITED", stat(solver, "unified_return_status"));
    assertTrue(solver.get_function("nlp_jac_g").is_a("SXFunction"));
  }

  private static final class IterCallback extends Callback {
    final List<List<DM>> foo = new ArrayList<>();
    final long nx, ng, np;
    final boolean wrong;

    IterCallback(long nx, long ng, long np, boolean wrong) {
      this.nx = nx;
      this.ng = ng;
      this.np = np;
      this.wrong = wrong;
      construct("mycallback", d());
    }

    @Override
    public long get_n_in() {
      return ca.nlpsol_n_out();
    }

    @Override
    public long get_n_out() {
      return 1;
    }

    @Override
    public Sparsity get_sparsity_in(long i) {
      String n = ca.nlpsol_out(i);
      if (n.equals("f")) return Sparsity.scalar();
      if (n.equals("x") || n.equals("lam_x")) return Sparsity.dense(nx);
      if (n.equals("g") || n.equals("lam_g")) return Sparsity.dense(wrong ? nx : ng);
      return new Sparsity(0, 0);
    }

    @Override
    public List<DM> eval(List<DM> arg) {
      foo.add(arg);
      return List.of(new DM(0.0));
    }
  }

  @Test
  public void test_iteration_Callback() {
    Assumptions.assumeTrue(hasNlpsol("ipopt"));
    SX x = SX.sym("x"), y = SX.sym("y");
    SX f = rosen(1, x, y);
    Map<String, SX> nlp = sxnlp("x", vc(x, y), "f", f, "g", x.plus(y));
    new Function("f", List.of(x, y), List.of(f));

    IterCallback mycallback = new IterCallback(2, 1, 0, false);
    Map<String, Object> opts = d("iteration_callback", mycallback, "ipopt.tol", 1e-8, "ipopt.max_iter", 50);
    Function solver = ca.nlpsol("solver", "ipopt", nlp, opts);
    solver.call(in("lbx", -10, "ubx", 10, "lbg", -10, "ubg", 10));
    assertEquals(iterCount(solver) + 1, mycallback.foo.size());

    IterCallback mycallback2 = new IterCallback(2, 1, 0, true);
    Map<String, Object> opts2 = d("iteration_callback", mycallback2, "ipopt.tol", 1e-8, "ipopt.max_iter", 50);
    try {
      ca.nlpsol("solver", "ipopt", nlp, opts2);
    } catch (RuntimeException e) {
      assertTrue(String.valueOf(e.getMessage()).contains("Callback function input size mismatch"), e.getMessage());
    }
    keep.add(mycallback);
    keep.add(mycallback2);
  }

  @Test
  public void test_pathological() {
    SX x = SX.sym("x"), y = SX.sym("y");
    Map<String, SX> nlp = sxnlp("x", vc(x, y), "f", sq(ca.minus(1.0, x)).plus(sq(y)));
    for (S s : solvers) {
      if (s.name.equals("snopt") || s.name.equals("worhp") || s.name.equals("stabilizedsqp")) continue;
      Function solver = ca.nlpsol("mysolver", s.name, nlp, s.opts);
      Map<String, DM> solver_in = in("x0", vc(1, 1), "lbx", vc(-10, -1), "ubx", vc(10, 2));
      Map<String, DM> solver_out = solver.call(solver_in);
      chk(solver_out.get("f"), vc(0), 7);
      chk(solver_out.get("x"), vc(1, 0), 7, s.name);
      if (!s.name.contains("bonmin")) chk(solver_out.get("lam_x"), vc(0, -0.0), 7, s.name);
      if (s.codegen != null) checkCodegen(solver, solver_in, s.codegen);
    }
  }

  @Test
  public void test_pathological2() {
    SX x = SX.sym("x"), y = SX.sym("y");
    Map<String, SX> nlp = sxnlp("x", vc(x, y), "f", sq(ca.minus(1.0, x)).plus(y));
    for (S s : solvers) {
      if (s.name.equals("snopt")) continue;
      Function solver = ca.nlpsol("mysolver", s.name, nlp, s.opts);
      Map<String, DM> solver_in = in("x0", vc(1, 1), "lbx", vc(-10, 0), "ubx", vc(10, 2));
      Map<String, DM> solver_out = solver.call(solver_in);
      chk(solver_out.get("f"), vc(0), 7);
      chk(solver_out.get("x"), vc(1, 0), 7);
      if (!s.name.contains("bonmin")) chk(solver_out.get("lam_x"), vc(0, -1), 7);
      if (s.codegen != null) checkCodegen(solver, solver_in, s.codegen);
    }
  }

  @Test
  public void test_pathological3() {
    SX x = SX.sym("x"), y = SX.sym("y");
    Map<String, SX> nlp = sxnlp("x", vc(x, y), "f", sq(ca.minus(1.0, x)), "g", x.plus(y));
    for (S s : solvers) {
      if (s.name.equals("worhp")) continue;
      Function solver = ca.nlpsol("mysolver", s.name, nlp, s.opts);
      Map<String, DM> solver_in = in("x0", vc(1, 1), "lbx", vc(-10, 0), "ubx", vc(10, 2), "lbg", vc(2), "ubg", vc(2));
      Map<String, DM> solver_out = solver.call(solver_in);
      chk(solver_out.get("f"), vc(0), 7);
      chk(solver_out.get("x"), vc(1, 1), 7);
      if (!s.name.contains("bonmin")) chk(solver_out.get("lam_x"), vc(0, 0), 7);
      if (s.codegen != null) checkCodegen(solver, solver_in, s.codegen);
    }
  }

  @Test
  public void test_pathological4() {
    SX x = SX.sym("x");
    Map<String, SX> nlp = sxnlp("x", x, "f", x.times(x));
    for (S s : solvers) {
      if (List.of("snopt", "worhp", "madnlp").contains(s.name)) continue;
      Function solver = ca.nlpsol("mysolver", s.name, nlp, s.opts);
      Map<String, DM> solver_in = in("x0", vc(0), "lbx", vc(0), "ubx", vc(0));
      Map<String, DM> solver_out = solver.call(solver_in);
      chk(solver_out.get("f"), vc(0), 7);
      chk(solver_out.get("x"), vc(0), 7);
      if (!s.name.contains("bonmin")) chk(solver_out.get("lam_x"), vc(0), 7);
      if (s.codegen != null) checkCodegen(solver, solver_in, s.codegen);
    }
  }

  @Test
  @SuppressWarnings("unchecked")
  public void test_nlp_sensitivity() {
    MX x = MX.sym("x"), p = MX.sym("p");
    Map<String, MX> nlp = mxnlp("x", x, "p", p, "f", sq(ca.sin(x).minus(sq(p))), "g", x);
    for (S s : solvers) {
      if (s.name.contains("madnlp") || s.name.contains("fatrop")) continue;
      if (s.in("ipopt") || s.in("snopt")) continue;
      if (s.name.equals("alpaqa")) continue;
      Map<String, Object> solver_options = s.opts;
      if (s.name.contains("knitro")) {
        solver_options = new LinkedHashMap<>(s.opts);
        solver_options.put("knitro", with((Map<String, Object>) s.opts.get("knitro"), "algorithm", 4));
      }
      Function solver = ca.nlpsol("mysolver", s.name, nlp, solver_options);
      solver.call(in("lbg", 0, "x0", 0, "p", 0.5));
      MX z = solver.call(min("p", p, "x0", x, "lbg", 0)).get("x");
      MX z2 = ca.asin(sq(p));
      Function f = new Function("f", List.of(x, p), List.of(z, ca.jacobian(z, p)));
      Function f2 = new Function("f", List.of(x, p), List.of(z2, ca.jacobian(z2, p)));
      checkfunctionLight(f, f2, List.of(new DM(0.0), new DM(0.5)), 6);
    }
  }

  @Test
  public void test_nlp_sensitivity_equality_zero_multiplier() {
    Assumptions.assumeTrue(hasConic("qrqp"));
    MX x = MX.sym("x", 3), p = MX.sym("p");
    Map<String, MX> nlp = mxnlp("x", x, "p", p, "f", sq(x.get(0).minus(p)), "g", x.get(1).minus(x.get(0)));
    Function solver = ca.nlpsol("solver", "sqpmethod", nlp, d("qpsol", "qrqp", "qpsol_options",
        d("print_iter", false, "print_header", false), "print_iteration", false, "print_header", false,
        "print_status", false, "sens_force_equality", true));
    Map<String, MX> r = solver.call(min("p", p, "lbg", 0, "ubg", 0, "lbx", vc(new MX(-inf), new MX(-inf), p),
        "ubx", vc(new MX(inf), new MX(inf), p)));
    DM lam = (DM) new Function("lam", List.of(p), List.of(vc(r.get("lam_g"), r.get("lam_x")))).invoke(0);
    assertEquals(0.0, ca.norm_inf(lam).to_double());
    MX z = r.get("x");
    MX z2 = vc(p, p, p);
    Function f = new Function("f", List.of(p), List.of(z, ca.jacobian(z, p)));
    Function f2 = new Function("f", List.of(p), List.of(z2, ca.jacobian(z2, p)));
    checkfunctionLight(f, f2, List.of(new DM(0.0)), 8);
  }

  @Test
  public void test_regularize_sqpmethod() {
    Assumptions.assumeTrue(hasConic("qrqp"));
    MX x = MX.sym("x", 2);
    MX f = ca.times(0.5, ca.bilin(dm(new double[][] {{1, 0}, {0, -2}}), x, x));
    Map<String, MX> nlp = mxnlp("x", x, "f", f, "g", x.get(1));

    Function solver = ca.nlpsol("mysolver", "sqpmethod", nlp, d("qpsol", "qrqp", "qpsol_options", d("print_problem", true)));
    Map<String, DM>[] res = resHolder();
    String result = captureStdout(() -> res[0] = solver.call(in("lbg", 2, "ubg", 2)));
    assertEquals(1, iterCount(solver));
    assertTrue(result.contains("H:\n[[1, 0], \n [0, -2]]"), result);
    chk(res[0].get("x"), vc(0, 2), 6);

    Function solver2 = ca.nlpsol("mysolver", "sqpmethod", nlp, d("qpsol", "qrqp", "qpsol_options", d("print_problem", true),
        "convexify_strategy", "regularize", "convexify_margin", 0));
    result = captureStdout(() -> res[0] = solver2.call(in("lbg", 2, "ubg", 2)));
    chk(res[0].get("x"), vc(0, 2), 6);
    assertEquals(2, iterCount(solver2));
    assertTrue(result.contains("H:\n[[3, 0], \n [0, 0]]"), result);

    Function solver3 = ca.nlpsol("mysolver", "sqpmethod", nlp, d("qpsol", "qrqp", "qpsol_options", d("print_problem", true),
        "convexify_strategy", "regularize", "convexify_margin", 1e-4));
    result = captureStdout(() -> res[0] = solver3.call(in("lbg", 2, "ubg", 2)));
    chk(res[0].get("x"), vc(0, 2), 6);
    assertEquals(2, iterCount(solver3));
    assertTrue(result.contains("H:\n[[3.0001, 0], \n [0, 0.0001]]"), result);

    x = MX.sym("x", 2);
    f = ca.times(0.5, ca.bilin(dm(new double[][] {{1, 0}, {0, 2}}), x, x));
    Map<String, MX> nlp2 = mxnlp("x", x, "f", f, "g", x.get(1));

    Function solver4 = ca.nlpsol("mysolver", "sqpmethod", nlp2, d("qpsol", "qrqp", "qpsol_options", d("print_problem", true)));
    result = captureStdout(() -> res[0] = solver4.call(in("lbg", 2, "ubg", 2)));
    assertEquals(1, iterCount(solver4));
    assertTrue(result.contains("H:\n[[1, 0], \n [0, 2]]"), result);
    chk(res[0].get("x"), vc(0, 2), 6);

    Function solver5 = ca.nlpsol("mysolver", "sqpmethod", nlp2, d("qpsol", "qrqp", "qpsol_options", d("print_problem", true),
        "convexify_strategy", "regularize"));
    result = captureStdout(() -> res[0] = solver5.call(in("lbg", 2, "ubg", 2)));
    chk(res[0].get("x"), vc(0, 2), 6);
    assertEquals(1, iterCount(solver5));
    assertTrue(result.contains("H:\n[[1, 0], \n [0, 2]]"), result);
  }

  @SuppressWarnings("unchecked")
  private static Map<String, DM>[] resHolder() {
    return (Map<String, DM>[]) new Map<?, ?>[1];
  }

  @Test
  public void test_infeasible() {
    MX x = MX.sym("x");
    Map<String, MX> nlp = mxnlp("f", sq(x), "x", x, "g", ca.sin(x));
    for (S s : solvers) {
      Map<String, Object> myoptions = with(s.opts, "error_on_fail", true);
      if (s.name.equals("sqpmethod")) continue;
      Function solver = ca.nlpsol("solver", s.name, nlp, myoptions);
      if (s.codegen != null && RUN_SLOW) {
        Function F = checkCodegen(solver, in("lbg", -10, "ubg", 10), s.codegen);
        assertThrows(RuntimeException.class, () -> F.call(in("lbg", 5, "ubg", 10)));
        checkCodegen(solver, in("lbg", 5, "ubg", 10), s.codegen, 15, false, true, new long[] {
            UnifiedReturnStatus.SOLVER_RET_UNKNOWN.swigValue(), UnifiedReturnStatus.SOLVER_RET_INFEASIBLE.swigValue(),
            UnifiedReturnStatus.SOLVER_RET_EXCEPTION.swigValue()});
      }
    }
  }

  @Test
  public void test_indefinite() {
    MX x = MX.sym("x", 2);
    MX f = ca.times(0.5, ca.bilin(dm(new double[][] {{1, 0}, {0, -2}}), x, x));
    Map<String, MX> nlp = mxnlp("x", x, "f", f, "g", x.get(1));
    for (S s : solvers) {
      Map<String, DM> solver_in = in("lbg", 2, "ubg", 2);
      Function solver = ca.nlpsol("mysolver", s.name, nlp, s.opts);
      Map<String, DM> out = solver.call(solver_in);
      chk(out.get("x"), vc(0, 2), 6);
      if (!s.name.contains("bonmin")) chk(out.get("lam_g"), vc(4), 6);
      if (s.codegen != null) {
        solver.generate("f.c", d("main", true));
        solver.generate_in("in.dat", solver.convert_in(solver_in));
        new File("f.c").deleteOnExit();
        new File("in.dat").deleteOnExit();
        checkCodegen(solver, solver_in, s.codegen);
      }
    }
  }

  @Test
  public void test_gauss_newton_sqpmethod() {
    Assumptions.assumeTrue(hasNlpsol("sqpmethod"));
    SX x = SX.sym("x", 3);
    SX F = ca.sin(x).minus(vc(1, 2, 3).times(vc(x.get(2), new SX(0), new SX(0))));
    SX J = ca.jacobian(F, x);
    SX f = ca.times(0.5, ca.dot(F, F));
    SX p = SX.sym("x", 0, 1), lam_f = SX.sym("x"), lam_g = SX.sym("x", 0, 1);
    Function GN = new Function("GN", List.of(x, p, lam_f, lam_g), List.of(lam_f.times(ca.triu(ca.mtimes(J.T(), J)))));
    Map<String, SX> nlp = sxnlp("x", x, "f", f);
    Map<String, Object> options = d("hess_lag", GN);
    assertInException("Hessian must be symmetric", () -> ca.nlpsol("solver", "sqpmethod", nlp, options));

    F = ca.sin(x.minus(vc(0.22, 0.72, 0.2))).minus(vc(0.1, 0.5, 0.99));
    J = ca.jacobian(F, x);
    f = ca.times(0.5, ca.dot(F, F));
    Function Hf = new Function("H", List.of(x), List.of((SX) ca.hessian(f, x).get(0)));
    DM H = (DM) Hf.invoke(0);
    double[] eig = eigh(H.full())[0][0];
    assertTrue(Arrays.stream(eig).anyMatch(e -> e < 0));

    GN = new Function("GN", List.of(x, p, lam_f, lam_g), List.of(lam_f.times(ca.mtimes(J.T(), J))));
    Function solver = ca.nlpsol("solver", "sqpmethod", sxnlp("x", x, "f", f), d("convexify_strategy", "regularize", "qpsol", "qrqp", "hess_lag", GN));
    solver.call(in());
    assertEquals(6, iterCount(solver));

    solver = ca.nlpsol("solver", "sqpmethod", sxnlp("x", x, "f", f), d("convexify_strategy", "regularize", "qpsol", "qrqp"));
    solver.call(in());
    assertEquals(9, iterCount(solver));
  }

  @Test
  public void test_gauss_newton_ipopt() {
    Assumptions.assumeTrue(hasNlpsol("ipopt"));
    SX x = SX.sym("x", 3);
    SX F = ca.sin(x).minus(vc(1, 2, 3).times(vc(x.get(2), new SX(0), new SX(0))));
    SX J = ca.jacobian(F, x);
    SX f = ca.times(0.5, ca.dot(F, F));
    SX p = SX.sym("x", 0, 1), lam_f = SX.sym("x"), lam_g = SX.sym("x", 0, 1);
    Function GN = new Function("GN", List.of(x, p, lam_f, lam_g), List.of(lam_f.times(ca.mtimes(J.T(), J))));
    Map<String, SX> nlp = sxnlp("x", x, "f", f);
    Map<String, Object> options = d("hess_lag", GN);
    assertInException("Hessian must be upper triangular", () -> ca.nlpsol("solver", "ipopt", nlp, options));

    F = ca.sin(x.minus(vc(0.22, 0.72, 0.2))).minus(vc(0.1, 0.5, 0.99));
    J = ca.jacobian(F, x);
    f = ca.times(0.5, ca.dot(F, F));
    Function Hf = new Function("H", List.of(x), List.of((SX) ca.hessian(f, x).get(0)));
    DM H = (DM) Hf.invoke(0);
    double[] eig = eigh(H.full())[0][0];
    assertTrue(Arrays.stream(eig).anyMatch(e -> e < 0));

    GN = new Function("GN", List.of(x, p, lam_f, lam_g), List.of(lam_f.times(ca.triu(ca.mtimes(J.T(), J)))));
    Function solver = ca.nlpsol("solver", "ipopt", sxnlp("x", x, "f", f), d("hess_lag", GN));
    solver.call(in());
    assertEquals(6, iterCount(solver));

    solver = ca.nlpsol("solver", "ipopt", sxnlp("x", x, "f", f));
    solver.call(in());
    assertEquals(9, iterCount(solver));
  }

  private static DM blk(DM a, DM b) {
    return ca.diagcat(List.of(a, b));
  }

  @Test
  public void test_cvx_sqpmethod() {
    Assumptions.assumeTrue(hasConic("qrqp"));
    double eps = 1e-2;
    DM H0 = blk(dm(new double[][] {{2, 3}, {3, 0}}), dm(new double[][] {{7, 1, 0}, {1, 0, 6}, {0, 6, 3}}));
    DM p = vc(0, 3, 2, 1, 4);
    DM Hp = H0.get(p, p);
    DM d3 = dm(new double[][] {{2, 0, 0}, {0, -1, 0}, {0, 0, 8}});
    DM d3r = dm(new double[][] {{2, 0, 0}, {0, 1, 0}, {0, 0, 8}});
    DM d3c = dm(new double[][] {{2, 0, 0}, {0, eps, 0}, {0, 0, 8}});
    DM d3h = dm(new double[][] {{2, 0, 0}, {0, eps / 2, 0}, {0, 0, 8}});
    DM b53 = dm(new double[][] {{5, 3}, {3, 6}}), b58 = dm(new double[][] {{5, 8}, {8, 6}}), b55 = dm(new double[][] {{5, 3}, {3, 5}});
    Map<String, Object> reg = d("convexify_strategy", "regularize", "convexify_margin", eps);
    Map<String, Object> refl = d("convexify_strategy", "eigen-reflect", "convexify_margin", eps);
    Map<String, Object> cl = d("convexify_strategy", "eigen-clip", "convexify_margin", eps);
    Object[][] cases = {
      {ca.sparsify(dm(new double[][] {{0, 0}, {0, 2}})), dm(new double[][] {{eps, 0}, {0, 2 + eps}}), reg, List.of(), List.of()},
      {blk(d3, b53), blk(d3r, b53), refl, List.of("Identified 2 blocks with maximum size 3"), List.of()},
      {blk(d3, b53), blk(d3c, b53), cl, List.of(), List.of()},
      {d3h, d3c, refl, List.of(), List.of()},
      {blk(d3, b58), blk(d3r, reflect(b58, eps)), refl, List.of(), List.of()},
      {blk(d3, b58), blk(d3c, clip(b58, eps)), cl, List.of(), List.of()},
      {d3h, d3c, refl, List.of(), List.of()},
      {blk(ca.sparsify(d3), b53), blk(d3r, b53), refl, List.of("Identified 4 blocks with maximum size 2"), List.of()},
      {blk(ca.sparsify(d3), b55), blk(d3r, b55), refl, List.of("Identified 4 blocks with maximum size 2"), List.of()},
      {Hp, reflect(Hp, eps), refl, List.of(), List.of()},
      {Hp, clip(Hp, eps), cl, List.of(), List.of()},
    };
    for (Object[] c : cases) {
      DM H = (DM) c[0], Hcvx = (DM) c[1];
      @SuppressWarnings("unchecked")
      Map<String, Object> opts = (Map<String, Object>) c[2];
      @SuppressWarnings("unchecked")
      List<String> includes_init = (List<String>) c[3], excludes_init = (List<String>) c[4];
      long n = H.size1();
      MX x = MX.sym("H", n);
      Map<String, MX> nlp = mxnlp("x", x, "f", ca.times(0.5, ca.bilin(H, x, x)));
      Map<String, Object> options = d("max_iter", 1, "qpsol", "qrqp", "verbose", true);
      options.putAll(opts);
      Function[] solver = new Function[1];
      assertOutput(includes_init, excludes_init, () -> solver[0] = ca.nlpsol("solver", "sqpmethod", nlp, options));
      DM x0 = DM.ones(n, 1);
      Map<String, DM> res = solver[0].call(in("x0", x0));
      chk(x0.minus(ca.solve(ca.densify(Hcvx), ca.mtimes(H, x0))), res.get("x"));
      checkSerialize(solver[0], in("x0", x0));
      checkCodegen(solver[0], in("x0", x0), d("std", "c99"), codegen_check_digits);
    }
  }

  private static DM asDM(Object o) {
    if (o instanceof DM) return (DM) o;
    if (o instanceof SX) return new DM((SX) o);
    if (o instanceof double[]) return new DM((double[]) o);
    throw new IllegalArgumentException(o == null ? "null" : o.getClass().getName());
  }

  @Test
  public void test_simple_bounds_detect() {
    SX x = SX.sym("x", 5), p = SX.sym("p", 1);
    Object[][] gl = {
      {1.1, x.get(0).times(x.get(1)), 2.0},
      {-inf, x.get(4), 2.0},
      {-10.0, x.get(0), 10.0},
      {-5.0, x.get(0), 2.0},
      {-4.0, x.get(0), 4.0},
      {1.1, x.get(4).times(x.get(1)), 2.0},
      {0.0, x.get(4), inf},
      {7.0, x.get(2), 7.0},
      {-4.0, x.get(2), 40.0},
      {9.0, x.get(1), 9.0},
      {-4.0, x.get(0), 4.0},
      {-4.0, x.get(1), 9.0}};
    double[] lbg0 = new double[gl.length], ubg0 = new double[gl.length];
    List<SX> gs = new ArrayList<>();
    for (int i = 0; i < gl.length; i++) {
      lbg0[i] = (Double) gl[i][0];
      gs.add((SX) gl[i][1]);
      ubg0[i] = (Double) gl[i][2];
    }
    SX g = ca.vertcat(gs);
    List<Object> r = ca.detect_simple_bounds(x, p, g, new SX(new DM(lbg0)), new SX(new DM(ubg0)));
    DM lbx = asDM(r.get(1)), ubx = asDM(r.get(2));
    Function lam_f = (Function) r.get(3), lam_b = (Function) r.get(4);
    chk(lbx.T(), row(-4, 9, 7, -inf, 0));
    chk(ubx.T(), row(2, 9, 7, inf, 2));

    java.util.function.Function<DM, DM> round_trip_f = a -> {
      List<DM> o = lam_f.call(List.of(a, new DM(0.0)));
      return lam_b.call(List.of(o.get(0), o.get(1), new DM(0.0))).get(0);
    };
    java.util.function.Function<List<DM>, List<DM>> round_trip_b = a -> {
      DM o = lam_b.call(List.of(a.get(0), a.get(1), new DM(0.0))).get(0);
      return lam_f.call(List.of(o, new DM(0.0)));
    };

    double[][] G = {
      {2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      {0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      {0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0},
      {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      {0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0},
      {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      {0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0},
      {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      {0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1},
      {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      {0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1}};
    for (int i = 0; i < 12; i++) {
      DM a = DM.zeros(12, 1);
      a.set(i, 2.0);
      DM b = round_trip_f.apply(a), c = round_trip_f.apply(round_trip_f.apply(a));
      chk(b, vc(G[i]));
      chk(b, c);
    }

    G = new double[][] {
      {-2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      {0, 0, 0, 0, -1, 0, 0, 0, 0, 0, -1, 0},
      {0, 0, 0, 0, 0, -2, 0, 0, 0, 0, 0, 0},
      {0, 0, 0, 0, 0, 0, -2, 0, 0, 0, 0, 0},
      {0, 0, 0, 0, 0, 0, 0, -2, 0, 0, 0, 0},
      {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      {0, 0, 0, 0, 0, 0, 0, 0, 0, -2, 0, 0},
      {0, 0, 0, 0, -1, 0, 0, 0, 0, 0, -1, 0},
      {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}};
    for (int i = 0; i < 12; i++) {
      DM a = DM.zeros(12, 1);
      a.set(i, -2.0);
      DM b = round_trip_f.apply(a), c = round_trip_f.apply(round_trip_f.apply(a));
      chk(b, vc(G[i]));
      chk(b, c);
    }

    for (double sgn : new double[] {-2, 2}) {
      double[][] G5 = new double[5][5];
      for (int i : new int[] {0, 1, 2, 4}) G5[i][i] = sgn;
      for (int i = 0; i < 5; i++) {
        DM a = DM.zeros(5, 1);
        a.set(i, sgn);
        List<DM> b = round_trip_b.apply(List.of(vc(3, 7), a));
        List<DM> c = round_trip_b.apply(round_trip_b.apply(List.of(vc(3, 7), a)));
        chk(b.get(0).T(), row(3, 7));
        chk(b.get(1), vc(G5[i]));
        chk(b.get(1), c.get(1));
      }
    }
  }

  @Test
  public void test_simple_bounds_detect_solvers() {
    Assumptions.assumeTrue(hasNlpsol("ipopt"));
    for (S s : solvers) {
      if (s.name.equals("sleqp") || s.name.equals("madnlp")) continue;
      MX x = MX.sym("x", 5);
      Object[][] gl = {
        {1.1, x.get(0).times(x.get(1)), 2.0},
        {-11.0, x.get(4), 2.0},
        {-10.0, x.get(0), 10.0},
        {-5.0, x.get(0), 2.0},
        {-4.0, x.get(0), 4.0},
        {1.1, x.get(4).times(x.get(1)), 2.0},
        {0.0, x.get(4), inf},
        {3.0, x.get(2), 3.0},
        {-4.0, x.get(2), 40.0},
        {2.0, x.get(1), 2.0},
        {-4.0, x.get(0), 4.0},
        {-4.0, x.get(1), 9.0}};
      double[] lbg0 = new double[gl.length], ubg0 = new double[gl.length];
      List<MX> gs = new ArrayList<>();
      for (int i = 0; i < gl.length; i++) {
        lbg0[i] = (Double) gl[i][0];
        gs.add((MX) gl[i][1]);
        ubg0[i] = (Double) gl[i][2];
      }
      DM lbg = vc(lbg0), ubg = vc(ubg0);
      MX g = ca.vertcat(gs);
      DM I = ca.densify(DM.eye(5));
      DM lbx = DM.inf(5).unaryMinus(), ubx = DM.inf(5);
      for (int i = 0; i < 3; i++) {
        for (double sign : new double[] {-1, 1}) {
          DM x_target = I.get(new Slice(), i).times(sign * 20);
          Map<String, MX> nlp = mxnlp("x", x, "g", g, "f", ca.sumsqr(x.minus(x_target)));
          Function solver_ref_ipopt = ca.nlpsol("mysolver", "ipopt", nlp);
          Function solver_ref = ca.nlpsol("mysolver", s.name, nlp, s.opts);
          Function solver = ca.nlpsol("mysolver", s.name, nlp, with(s.opts, "detect_simple_bounds", true));
          Map<String, DM> solver_in = in("lbg", lbg, "ubg", ubg, "lbx", lbx, "ubx", ubx);
          Map<String, DM> sol_ref_ipopt = solver_ref_ipopt.call(solver_in);
          solver_in.put("x0", sol_ref_ipopt.get("x").times(1.2));
          int digits = 6;
          if (s.in("daqp")) continue;
          if (s.in("worhp")) digits = 4;
          solver_ref.call(solver_in);
          if (s.name.equals("bonmin")) {
            Map<String, DM> solver_ref_out = solver_ref.call(solver_in), solver_out = solver.call(solver_in);
            for (String output : List.of("x", "f"))
              chk(solver_out.get(output), solver_ref_out.get(output), digits, s.toString());
          } else {
            checkfunctionLight(solver, solver_ref, positional(solver, solver_in), digits, s.toString());
          }
          if (s.codegen != null) checkCodegen(solver, solver_in, s.codegen);
        }
      }
    }
  }

  private static <T extends GenericMatrix> void signIssueCase(Map<String, T> nlp, double lbg, double ubg) {
    Function solver_ref = ca.nlpsol("solver", "ipopt", nlp);
    Function solver = ca.nlpsol("solver", "ipopt", nlp, d("detect_simple_bounds", true));
    checkfunctionLight(solver, solver_ref, solver.convert_in(in("lbg", lbg, "ubg", ubg)), 9);
  }

  @Test
  public void test_simple_bounds_sign_issue() {
    Assumptions.assumeTrue(hasNlpsol("ipopt"));
    double[][] b = {{-2, 3}, {2, 2}, {-2, 3}, {2, 2}, {-2, 3}, {-2, 3}, {-2, 3}};
    double[] c = {1, 1, 0.3, 0.7, -1, -0.7, 0.3};
    SX xs = SX.sym("x");
    for (SX f : List.of(xs, xs.unaryMinus()))
      for (int k = 0; k < c.length; k++) signIssueCase(sxnlp("x", xs, "f", f, "g", ca.times(c[k], xs)), b[k][0], b[k][1]);
    MX xm = MX.sym("x");
    for (MX f : List.of(xm, xm.unaryMinus()))
      for (int k = 0; k < c.length; k++) signIssueCase(mxnlp("x", xm, "f", f, "g", ca.times(c[k], xm)), b[k][0], b[k][1]);
  }

  @Test
  public void test_simple_bounds_equality() {
    Assumptions.assumeTrue(hasNlpsol("ipopt"));
    MX x = MX.sym("x"), y = MX.sym("y");
    MX z = vc(x, y);
    MX f = x;
    MX g = vc(sq(x).plus(sq(y)), x, x.times(y));
    DM lbg = vc(1, -2, -5), ubg = vc(1, 2, 5), x0 = vc(0.2, 0.6);
    Function solver_ref = ca.nlpsol("solver", "ipopt", mxnlp("x", z, "f", f, "g", g), d("equality", new boolean[] {true, false, false}));
    Function solver = ca.nlpsol("solver", "ipopt", mxnlp("x", z, "f", f, "g", g),
        d("detect_simple_bounds", true, "equality", new boolean[] {true, false, false}));
    checkfunctionLight(solver, solver_ref, solver.convert_in(in("x0", x0, "lbg", lbg, "ubg", ubg)), 9);
  }

  @Test
  public void test_issue_3407() {
    Assumptions.assumeTrue(hasNlpsol("ipopt"));
    Opti opti = new Opti();
    MX x = opti.variable();
    opti.minimize(sq(x));
    Function f = new Function("f", List.of(x), List.of(ca.times(3.0, x), sq(x.attachAssert(x.gt(1.0)))), d("never_inline", true));
    List<MX> yz = f.call(List.of(x));
    opti.subject_to(yz.get(0).le(10.0));
    opti.subject_to(ca.eq(yz.get(1), 3.0));
    opti.set_initial(x, 3);
    opti.solver("ipopt", d("detect_simple_bounds", true));
    opti.solve();
  }

  @Test
  public void test_issue_4235() {
    Assumptions.assumeTrue(hasNlpsol("ipopt"));
    for (boolean expand : new boolean[] {false, true}) {
      Opti opti = new Opti();
      MX x = opti.variable();
      opti.minimize(sq(x));
      Function f = new Function("f", List.of(x), List.of(ca.times(3.0, x), sq(x.attachAssert(x.gt(1000.0)))));
      List<MX> yz = f.call(List.of(ca.times(2.0, x)));
      opti.subject_to(yz.get(0).le(10.0));
      opti.subject_to(ca.eq(yz.get(1), 3.0));
      opti.set_initial(x, 3);
      if (!expand) {
        opti.solver("ipopt", d("detect_simple_bounds", true));
        assertThrows(RuntimeException.class, opti::solve);
      } else {
        // Expand will have the side effect of removing the assertion
        opti.solver("ipopt", d("expand", true, "detect_simple_bounds", true));
        opti.solve();
      }
    }
  }

  @Test
  public void test_simple_bounds_detect2() {
    Assumptions.assumeTrue(hasNlpsol("ipopt"));
    MX x = MX.sym("x", 5), p = MX.sym("p", 5);
    Object[][] gl = {
      {1.1, x.get(0).times(x.get(1)), 2.0},
      {-11.0, x.get(4), 2.0},
      {-10.0, x.get(0), 10.0},
      {-5.0, x.get(0), 2.0},
      {-4.0, x.get(0), 4.0},
      {1.1, x.get(4).times(x.get(1)), 2.0},
      {0.0, x.get(4), inf},
      {3.0, x.get(2), 3.0},
      {-4.0, x.get(2), 40.0},
      {2.0, x.get(1), 2.0},
      {-4.0, x.get(0), 4.0},
      {-4.0, x.get(1), 9.0}};
    double[] lbg0 = new double[gl.length], ubg0 = new double[gl.length];
    List<MX> gs = new ArrayList<>();
    for (int i = 0; i < gl.length; i++) {
      lbg0[i] = (Double) gl[i][0];
      gs.add((MX) gl[i][1]);
      ubg0[i] = (Double) gl[i][2];
    }
    DM lbg = vc(lbg0), ubg = vc(ubg0);
    MX g = ca.vertcat(gs);
    Map<String, MX> nlp = mxnlp("x", x, "p", p, "g", g, "f", ca.sumsqr(x.minus(p)));
    int print_level = 0;

    assertOutput(List.of("Total number of inequality constraints...............:        2"), List.of(),
        () -> ca.nlpsol("solver", "ipopt", nlp, d("detect_simple_bounds", true)).call(in("x0", 0, "lbg", lbg, "ubg", ubg)));
    assertOutput(List.of("Total number of inequality constraints...............:       10"), List.of(),
        () -> ca.nlpsol("solver", "ipopt", nlp, d("detect_simple_bounds", false)).call(in("x0", 0, "lbg", lbg, "ubg", ubg)));

    Map<String, Object> options = d("detect_simple_bounds", true, "ipopt.print_level", print_level, "print_time", false,
        "ipopt.tol", 1e-12, "ipopt.fixed_variable_treatment", "relax_bounds", "ipopt.bound_relax_factor", 1e-12);
    Function solver = ca.nlpsol("solver", "ipopt", nlp, options);
    Map<String, Object> options_nominal = d("detect_simple_bounds", false, "ipopt.print_level", print_level,
        "print_time", false, "ipopt.tol", 1e-12, "ipopt.bound_relax_factor", 1e-12);
    Function solver_nominal = ca.nlpsol("solver", "ipopt", nlp, options_nominal);
    Map<String, Object> w_options = d("ipopt.warm_start_init_point", "yes", "ipopt.warm_start_bound_push", 1e-16,
        "ipopt.bound_push", 1e-16, "ipopt.warm_start_mult_bound_push", 1e-16, "ipopt.tol", 1e-7);
    Map<String, Object> m1 = new LinkedHashMap<>(options), m2 = new LinkedHashMap<>(options_nominal);
    m1.putAll(w_options);
    m2.putAll(w_options);
    Function wsolver = ca.nlpsol("solver", "ipopt", nlp, m1);
    Function wsolver_nominal = ca.nlpsol("solver", "ipopt", nlp, m2);

    DM I = ca.densify(DM.eye(5));
    DM[][] bounds = {
      {DM.inf(5).unaryMinus(), DM.inf(5)},
      {ca.times(-3.0, DM.ones(5)), ca.times(3.0, DM.ones(5))},
      {ca.times(-4.0, DM.ones(5)), ca.times(4.0, DM.ones(5))}};
    for (DM[] bnd : bounds) {
      DM lbx = bnd[0], ubx = bnd[1];
      for (int i = 0; i < 5; i++) {
        for (double sign : new double[] {-1, 1}) {
          DM x_target = I.get(new Slice(), i).times(sign * 20);
          Map<String, DM> out = solver.call(in("p", x_target, "lbg", lbg, "ubg", ubg, "lbx", lbx, "ubx", ubx));
          Map<String, DM> out_nominal = solver_nominal.call(in("p", x_target, "lbg", lbg, "ubg", ubg, "lbx", lbx, "ubx", ubx));
          if (el(lbx, 0) == -inf) chk(out.get("lam_x"), vc(rep(0, 5)));
          chk(out.get("x"), out_nominal.get("x"), 8);
          chk(out.get("g"), out_nominal.get("g"), 8);

          wsolver.call(in("x0", out.get("x"), "lam_x0", out.get("lam_x"), "lam_g0", out.get("lam_g"), "p", x_target,
              "lbg", lbg, "ubg", ubg, "lbx", lbx, "ubx", ubx));
          assertEquals(0, iterCount(wsolver));
          wsolver_nominal.call(in("x0", out.get("x"), "lam_x0", out.get("lam_x"), "lam_g0", out.get("lam_g"), "p", x_target,
              "lbg", lbg, "ubg", ubg, "lbx", lbx, "ubx", ubx));
          assertEquals(0, iterCount(wsolver_nominal));
          wsolver.call(in("x0", out_nominal.get("x"), "lam_x0", out_nominal.get("lam_x"), "lam_g0", out_nominal.get("lam_g"),
              "p", x_target, "lbg", lbg, "ubg", ubg, "lbx", lbx, "ubx", ubx));
          assertEquals(0, iterCount(wsolver));
        }
      }
    }
  }

  @Test
  public void test_derivative() {
    MX x = MX.sym("x", 3), p = MX.sym("p", 3);
    Map<String, MX> nlp = mxnlp("x", x, "p", p, "f", ca.sumsqr(ca.sin(x).minus(p)), "g", x);
    Function solver = ca.nlpsol("solver", "sqpmethod", nlp, d("qpsol", "qrqp"));
    Map<String, MX> res = solver.call(min("p", p, "x0", 0, "lbg", -inf, "ubg", inf));
    Function f = new Function("f", List.of(p), List.of(res.get("x")));
    Function f_ref = new Function("f_ref", List.of(p), List.of(ca.asin(p)));
    checkfunction(f, f_ref, List.of(vc(0.1, 0.2, 0.3)));

    x = MX.sym("x", 3);
    p = MX.sym("p", 3);
    nlp = mxnlp("x", x, "p", p, "f", ca.sumsqr(ca.sin(x).minus(p)));
    solver = ca.nlpsol("solver", "sqpmethod", nlp, d("qpsol", "qrqp"));
    res = solver.call(min("p", p, "x0", 0));
    f = new Function("f", List.of(p), List.of(res.get("x")));
    f_ref = new Function("f_ref", List.of(p), List.of(ca.asin(p)));
    checkfunction(f, f_ref, List.of(vc(0.1, 0.2, 0.3)));
  }

  // test_hock_schittkowski: returns immediately in python (body needs pandas and the hock_schittkowski modules)
  @Test
  public void test_hock_schittkowski() {}

  @Test
  public void test_exception_in_oraclefunction() {
    MX x = MX.sym("x");
    MX x_fail = x.attachAssert(ca.ne(x, x), "Cuckoo");
    Map<String, MX> nlp = mxnlp("x", x, "f", sq(x.minus(1.0)), "g", x_fail);
    for (S s : solvers) {
      if (s.name.contains("madnlp")) continue;
      Function solver = ca.nlpsol("mysolver", s.name, nlp, s.opts);
      Map<String, DM> solver_in = in("lbx", vc(-10), "ubx", vc(10), "lbg", vc(-10), "ubg", vc(10));
      assertInAnyOutput("Cuckoo", () -> solver.call(solver_in));
    }
  }

  @Test
  public void test_unsolved_stats() {
    MX x = MX.sym("x");
    MX x_fail = x.attachAssert(ca.ne(x, x), "Cuckoo");
    Map<String, MX> nlp = mxnlp("x", x, "f", sq(x.minus(1.0)), "g", x_fail);
    for (S s : solvers) {
      Function solver = ca.nlpsol("mysolver", s.name, nlp, s.opts);
      assertInException("No stats available", solver::stats);
    }
  }

  @Test
  public void test_ipopt_detect_simple_bounds() {
    Assumptions.assumeTrue(hasNlpsol("ipopt"));
    MX x = MX.sym("x"), y = MX.sym("y"), z = MX.sym("z"), w = MX.sym("w");
    List<Object[]> testcases = new ArrayList<>();
    for (double d_scale : new double[] {1.0, 1e5}) {
      for (double c_scale : new double[] {1.0, 1e5}) {
        MX g = vc(x, x.plus(y).times(d_scale), sq(x).plus(sq(z)), z, z.plus(y).plus(w).times(c_scale), sq(z).plus(sq(y)).plus(w));
        DM lbg = vc(-5.0, -2 * d_scale, -3.0, 1.0, 2 * c_scale, 3.0), ubg = vc(3.0, 2 * d_scale, 3.0, 1.0, 2 * c_scale, 3.0);
        testcases.add(new Object[] {g, lbg, ubg, vc(-inf, -inf, -inf, -inf), vc(inf, inf, inf, inf), new double[] {0, 3, 5, 1, 4, 2}});
        g = vc(x.plus(y).times(d_scale), sq(x).plus(sq(z)), z.plus(y).plus(w).times(c_scale), sq(z).plus(sq(y)).plus(w));
        lbg = vc(-2 * d_scale, -3.0, 2 * c_scale, 3.0);
        ubg = vc(2 * d_scale, 3.0, 2 * c_scale, 3.0);
        testcases.add(new Object[] {g, lbg, ubg, vc(-5, -inf, 1, -inf), vc(3, inf, 1, inf), new double[] {3, 0, 2, 1}});
      }
    }
    MX X = vc(x, y, z, w);
    MX F = sq(x).plus(ca.times(2.0, sq(y))).plus(ca.times(3.0, sq(z))).plus(ca.times(8.0, sq(w)));
    for (Object[] tc : testcases) {
      MX g = (MX) tc[0];
      DM lbg = (DM) tc[1], ubg = (DM) tc[2], lbx = (DM) tc[3], ubx = (DM) tc[4];
      double[] perm = (double[]) tc[5];
      Function solver_ref = ca.nlpsol("solver", "ipopt", mxnlp("x", X, "g", g, "f", F));
      Map<String, DM> ref = solver_ref.call(in("lbg", lbg, "ubg", ubg, "lbx", lbx, "ubx", ubx));
      for (DM g_perm : List.of(new DM(range((int) g.numel())), new DM(perm))) {
        Map<String, MX> nlp = mxnlp("x", X, "g", g.get(g_perm), "f", F);
        Function solver_ref2 = ca.nlpsol("solver", "ipopt", nlp);
        Map<String, DM> ref2 = solver_ref2.call(in("lbg", lbg.get(g_perm), "ubg", ubg.get(g_perm), "lbx", lbx, "ubx", ubx));
        Map<String, DM> arg = in("x0", ref.get("x").times(1.1), "lbg", lbg.get(g_perm), "ubg", ubg.get(g_perm), "lbx", lbx, "ubx", ubx);
        for (boolean detect_simple_bounds : new boolean[] {true, false}) {
          for (String fixed_variable_treatment : List.of("make_constraint", "make_parameter", "make_parameter_nodual")) {
            for (String start_with_resto : List.of("no", "yes")) {
              List<String> keys2 = fixed_variable_treatment.equals("make_constraint") ? List.of("x", "f", "g", "lam_g") : List.of();
              Map<String, Object> o = d("ipopt.fixed_variable_treatment", fixed_variable_treatment,
                  "detect_simple_bounds", detect_simple_bounds, "ipopt.start_with_resto", start_with_resto);
              Function solver = ca.nlpsol("solver", "ipopt", nlp, o);
              Map<String, DM> res = solver.call(arg);
              for (String k : List.of("x", "f")) chk(res.get(k), ref.get(k), 7);
              for (String k : keys2) chk(res.get(k), ref2.get(k), 7);
              assertTrue(success(solver));

              Map<String, Object> cache = d("nlp_jac_g", solver.get_function("nlp_jac_g"), "nlp_hess_l", solver.get_function("nlp_hess_l"));
              solver = ca.nlpsol("solver", "ipopt", nlp, with(o, "cache", cache));
              res = solver.call(arg);
              for (String k : List.of("f", "x")) chk(res.get(k), ref.get(k), 7);
              for (String k : keys2) chk(res.get(k), ref2.get(k), 7);

              cache = d("nlp_jac_g", solver_ref2.get_function("nlp_jac_g"), "nlp_hess_l", solver_ref2.get_function("nlp_hess_l"));
              solver = ca.nlpsol("solver", "ipopt", nlp, with(o, "cache", cache));
              res = solver.call(arg);
              for (String k : List.of("f", "x")) chk(res.get(k), ref.get(k), 7);
              for (String k : keys2) chk(res.get(k), ref2.get(k), 7);
              assertTrue(success(solver));
            }
          }
        }
      }
    }
  }

  @Test
  public void test_ipopt_custom_hess() {
    Assumptions.assumeTrue(hasNlpsol("ipopt"));
    SX x = SX.sym("x"), y = SX.sym("y"), w = SX.sym("w");
    SX p = new SX(0, 1);
    SX f = sq(ca.minus(1.0, x)).plus(ca.times(100.0, sq(w)));
    SX g = y.minus(sq(x)).minus(w);
    SX X = vc(x, y, w);
    Map<String, SX> nlp = sxnlp("x", X, "f", f, "g", g);
    SX lam_f = SX.sym("lam_f"), lam_g = SX.sym("lam_g", g.numel());
    Function ref_solver = ca.nlpsol("solver", "ipopt", nlp);
    Map<String, DM> ref_sol = ref_solver.call(in("x0", 0, "lbg", 0, "ubg", 0));

    Function nlp_hess_l_custom = new Function("nlp_hess_l", List.of(X, p, lam_f, lam_g), List.of(new SX(ca.triu(DM.zeros(3, 3)))));
    Function solver = ca.nlpsol("solver", "ipopt", nlp, d("cache", d("nlp_hess_l", nlp_hess_l_custom)));
    List<DM> hin = List.of(vc(0.11, 0.3, 0.7), new DM(0.0), new DM(1.2), new DM(3.7));
    Function s1 = solver;
    // python's assertInException catches the checkarray AssertionError ("evaluation error")
    AssertionError err = assertThrows(AssertionError.class,
        () -> checkfunctionLight(ref_solver.get_function("nlp_hess_l"), s1.get_function("nlp_hess_l"), hin, 9));
    assertTrue(err.getMessage().contains("evaluation error"), err.getMessage());

    SX lag = lam_f.times(f).plus(ca.dot(lam_g, g));
    SX H = ca.jacobian(ca.gradient(lag, X), X, d("symmetric", true));
    nlp_hess_l_custom = new Function("nlp_hess_l", List.of(X, p, lam_f, lam_g), List.of(ca.triu(H)));
    ref_solver.get_function("nlp_hess_l").disp(true);
    nlp_hess_l_custom.disp(true);
    solver = ca.nlpsol("solver", "ipopt", nlp, d("cache", d("nlp_hess_l", nlp_hess_l_custom)));
    Map<String, DM> sol = solver.call(in("x0", 0, "lbg", 0, "ubg", 0));
    assertTrue(success(solver));
    chk(sol.get("x"), ref_sol.get("x"), 6);
    checkfunctionLight(ref_solver.get_function("nlp_hess_l"), solver.get_function("nlp_hess_l"), hin, 9);
  }

  @Test
  public void test_ipopt_custom_jac() {
    Assumptions.assumeTrue(hasNlpsol("ipopt"));
    SX x = SX.sym("x"), y = SX.sym("y"), w = SX.sym("w");
    SX p = new SX(0, 1);
    SX f = sq(ca.minus(1.0, x)).plus(ca.times(100.0, sq(w)));
    SX g = y.minus(sq(x)).minus(w);
    SX X = vc(x, y, w);
    Map<String, SX> nlp = sxnlp("x", X, "f", f, "g", g);
    Function ref_solver = ca.nlpsol("solver", "ipopt", nlp);
    ref_solver.call(in("x0", 0, "lbg", 0, "ubg", 0));
    Map<String, Object> options = d();
    Function solver = ca.nlpsol("solver", "ipopt", nlp, options);
    Function nlp_jac_g_custom = new Function("nlp_jac_g", List.of(X, p), List.of(g, ca.jacobian(g, X)), List.of("x", "p"), List.of("g", "jac_g_x"));
    options.put("nlp_jac_g", d("nlp_jac_g", nlp_jac_g_custom));
    solver.call(in("x0", 0, "lbg", 0, "ubg", 0));
    assertTrue(success(solver));
    checkfunctionLight(ref_solver.get_function("nlp_jac_g"), solver.get_function("nlp_jac_g"), List.of(vc(0.11, 0.3, 0.7), new DM(0.0)), 9);
  }

  @Test
  public void test_nlp_signature() {
    MX x = MX.sym("x");
    assertInException("must have inputs (x, p) and outputs (f, g)",
        () -> ca.nlpsol("solver", "ipopt", new Function("nlp", List.of(x), List.of(sq(x)), List.of("a"), List.of("b"))));
  }

  @Test
  public void test_option_propagation() {
    Assumptions.assumeTrue(hasNlpsol("ipopt"));
    MX x = MX.sym("x");
    MX y = sq(x.printme(0.0));
    Map<String, Object> options = d("common_options", d("helper_options", d("enable_fd", true, "enable_forward", false,
        "enable_reverse", false)), "ipopt", d("resto.tol", 1e-8, "hessian_approximation", "limited-memory", "max_iter", 0));
    Function solver = ca.nlpsol("solver", "ipopt", mxnlp("x", x, "f", y), options);
    String[] result = captureBoth(() -> solver.call(in("x0", 1)));
    assertEquals(0, result[1].length(), result[1]);
    List<Double> r = new ArrayList<>();
    for (String l : result[0].split("\n"))
      if (l.contains("|> 0 : ")) r.add(Double.parseDouble(l.split(":")[1].trim()));
    double spread = r.stream().mapToDouble(Double::doubleValue).max().orElse(0)
        - r.stream().mapToDouble(Double::doubleValue).min().orElse(0);
    assertTrue(spread > 0, r.toString());
    assertTrue(spread < 1e-4, r.toString());
  }

  @Test
  public void test_thread_map_reuses_solver_memory() {
    List<Object[]> cases = List.of(
        new Object[] {"ipopt", d()},
        new Object[] {"ipopt", d("ipopt.print_level", 0, "ipopt.sb", "yes", "print_time", false)},
        new Object[] {"bonmin", d()});
    Assumptions.assumeTrue(hasNlpsol("ipopt") || hasNlpsol("bonmin"), "needs ipopt or bonmin");
    for (Object[] c : cases) {
      String plugin = (String) c[0];
      @SuppressWarnings("unchecked")
      Map<String, Object> opts = (Map<String, Object>) c[1];
      if (!hasNlpsol(plugin)) continue;
      MX x = MX.sym("x", 2), p = MX.sym("p");
      Function solver = ca.nlpsol("solver", plugin, mxnlp("x", x, "p", p, "f", ca.sumsqr(x.minus(p))), opts);
      MX P = MX.sym("P");
      Function par = new Function("par", List.of(P), List.of(solver.call(min("x0", 0, "p", P)).get("x").get(0)));
      Function m = par.map(8, "thread", 4);
      DM ps = new DM(range(8)).T();
      for (int k = 0; k < 5; k++) chk((DM) m.invoke(ps), ps, 6);
    }
  }
}

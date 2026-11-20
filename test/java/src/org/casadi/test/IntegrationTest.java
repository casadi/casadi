package org.casadi.test;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.io.File;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;

import org.casadi.CodeGenerator;
import org.casadi.DM;
import org.casadi.Function;
import org.casadi.GlobalOptions;
import org.casadi.MX;
import org.casadi.SX;
import org.casadi.Slice;
import org.casadi.Sparsity;
import org.casadi.ca;
import org.junit.jupiter.api.Assumptions;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;

/** Port of test/python/integration.py */
public class IntegrationTest extends CasadiTestCase {

  // ---------------------------------------------------------------- helpers

  private static final class Intg {
    final String name;
    final List<String> features;
    final Map<String, Object> options;

    Intg(String name, List<String> features, Map<String, Object> options) {
      this.name = name;
      this.features = features;
      this.options = options;
    }
  }

  private static final class EvIntg {
    final String name;
    final Map<String, Object> options;
    final boolean events, eventsDae, eventsAdjoint, codegen;
    final int digits;

    EvIntg(String name, Map<String, Object> options, boolean events, boolean eventsDae, boolean eventsAdjoint,
        boolean codegen, int digits) {
      this.name = name;
      this.options = options;
      this.events = events;
      this.eventsDae = eventsDae;
      this.eventsAdjoint = eventsAdjoint;
      this.codegen = codegen;
      this.digits = digits;
    }
  }

  private static final List<Intg> integrators = new ArrayList<>();
  private static final List<EvIntg> eventIntegrators = new ArrayList<>();

  static {
    if (ca.has_integrator("cvodes"))
      integrators.add(new Intg("cvodes", List.of("ode"),
          opts("abstol", 1e-15, "reltol", 1e-15, "fsens_err_con", true, "quad_err_con", false)));
    if (ca.has_integrator("idas"))
      integrators.add(new Intg("idas", List.of("dae", "ode"),
          opts("abstol", 1e-15, "reltol", 1e-15, "fsens_err_con", true, "calc_icB", true)));
    integrators.add(new Intg("collocation", List.of("dae", "ode"),
        opts("rootfinder", "newton", "number_of_finite_elements", 18)));
    integrators.add(new Intg("rk", List.of("ode"), opts("number_of_finite_elements", 1000)));
    integrators.add(new Intg("collocation", List.of("dae", "ode"),
        opts("number_of_finite_elements", 18, "simplify", true, "rootfinder", "fast_newton")));
    integrators.add(new Intg("rk", List.of("ode"), opts("number_of_finite_elements", 1000, "simplify", true)));
    integrators.add(new Intg("rosenbrock", List.of("dae", "ode"), opts("abstol", 1e-12, "reltol", 1e-12)));

    for (String plugin : List.of("cvodes", "idas")) {
      if (!ca.has_integrator(plugin)) continue;
      eventIntegrators.add(new EvIntg(plugin,
          opts("abstol", 1e-12, "reltol", 1e-12, "quad_err_con", true, "max_event_iter", 20),
          true, false, false, false, 5));
    }
    eventIntegrators.add(new EvIntg("rk", opts("number_of_finite_elements", 100, "max_event_iter", 20),
        true, false, false, false, 7));
    eventIntegrators.add(new EvIntg("collocation", opts("number_of_finite_elements", 20, "max_event_iter", 20),
        true, false, false, false, 7));
    Object[][] rt = {{"erk", opts()}, {"erk", opts("method", "rk4", "step_size", 0.002)},
        {"rosenbrock", opts()}, {"rosenbrock", opts("step_size", 0.002)}};
    for (Object[] e : rt) {
      String plugin = (String) e[0];
      @SuppressWarnings("unchecked")
      Map<String, Object> o = with((Map<String, Object>) e[1], "abstol", 1e-11, "reltol", 1e-11,
          "max_tape_steps", 600);
      eventIntegrators.add(new EvIntg(plugin, o, true, plugin.equals("rosenbrock"), true, true, 8));
    }
  }

  private static Map<String, Object> opts(Object... kv) {
    Map<String, Object> m = new LinkedHashMap<>();
    for (int i = 0; i < kv.length; i += 2) m.put((String) kv[i], kv[i + 1]);
    return m;
  }

  private static Map<String, Object> with(Map<String, ?> base, Object... kv) {
    Map<String, Object> m = new LinkedHashMap<>(base);
    m.putAll(opts(kv));
    return m;
  }

  private static DM toDM(Object v) {
    if (v instanceof DM) return (DM) v;
    if (v instanceof Number) return new DM(((Number) v).doubleValue());
    if (v instanceof double[]) return new DM((double[]) v);
    if (v instanceof double[][]) return new DM((double[][]) v);
    throw new IllegalArgumentException(String.valueOf(v));
  }

  private static SX toSX(Object v) {
    if (v instanceof SX) return (SX) v;
    return new SX(toDM(v));
  }

  private static MX toMX(Object v) {
    if (v instanceof MX) return (MX) v;
    return new MX(toDM(v));
  }

  private static Map<String, SX> sxd(Object... kv) {
    Map<String, SX> m = new LinkedHashMap<>();
    for (int i = 0; i < kv.length; i += 2) m.put((String) kv[i], toSX(kv[i + 1]));
    return m;
  }

  private static Map<String, MX> mxd(Object... kv) {
    Map<String, MX> m = new LinkedHashMap<>();
    for (int i = 0; i < kv.length; i += 2) m.put((String) kv[i], toMX(kv[i + 1]));
    return m;
  }

  private static Map<String, DM> dmd(Object... kv) {
    Map<String, DM> m = new LinkedHashMap<>();
    for (int i = 0; i < kv.length; i += 2) m.put((String) kv[i], toDM(kv[i + 1]));
    return m;
  }

  private static SX vs(Object... a) {
    List<SX> l = new ArrayList<>();
    for (Object e : a) l.add(toSX(e));
    return ca.vertcat(l);
  }

  private static MX vm(Object... a) {
    List<MX> l = new ArrayList<>();
    for (Object e : a) l.add(toMX(e));
    return ca.vertcat(l);
  }

  private static double[] linspace(double a, double b, int n) {
    // numpy.linspace rounding: i*step + start, last point exact
    double[] r = new double[n];
    double step = (b - a) / (n - 1);
    for (int i = 0; i < n; i++) r[i] = i * step + a;
    if (n > 1) r[n - 1] = b;
    return r;
  }

  private static double num(Object o) {
    return ((Number) o).doubleValue();
  }

  /** Python casadiTestCase.assertAlmostEqual: places after scaling large values */
  private static void almost(double first, double second, int places, String msg) {
    if (Double.isNaN(first) && Double.isNaN(second)) return;
    if (first == second) return;
    double n = Math.max(Math.abs(first), Math.abs(second));
    n = n > 1e3 ? Math.pow(10, Math.floor(Math.log10(n))) : 1.0;
    double diff = Math.abs(first / n - second / n);
    assertTrue(diff < 0.5 * Math.pow(10, -places),
        msg + " " + first + " <-> " + second + " [places=" + places + "]");
  }

  /** Python casadiTestCase.checkarray(zr, zt, digits) */
  private static void chk(DM zr, DM zt, int digits, String msg) {
    assertEquals(zr.size1(), zt.size1(), msg + ": size1 " + zr + " <-> " + zt);
    assertEquals(zr.size2(), zt.size2(), msg + ": size2 " + zr + " <-> " + zt);
    double[][] a = zr.full(), b = zt.full();
    for (int i = 0; i < a.length; i++)
      for (int j = 0; j < a[i].length; j++)
        almost(b[i][j], a[i][j], digits, msg + " at (" + i + "," + j + "): " + zt + " <-> " + zr);
  }

  private static void chk(DM zr, DM zt, int digits) {
    chk(zr, zt, digits, "");
  }

  private static void chk(DM zr, double zt, int digits) {
    chk(zr, new DM(zt), digits, "");
  }

  private static Function jacobianOld(Function f, long i, long j) {
    List<String> out = new ArrayList<>();
    out.add("jac:" + f.name_out(j) + ":" + f.name_in(i));
    out.addAll(f.name_out());
    return f.factory(f.name() + "_jac", f.name_in(), out);
  }

  private static Function hessianOld(Function f, long i, long j) {
    List<String> out = new ArrayList<>();
    out.add("hess:" + f.name_out(j) + ":" + f.name_in(i) + ":" + f.name_in(i));
    out.add("grad:" + f.name_out(j) + ":" + f.name_in(i));
    out.addAll(f.name_out());
    return f.factory(f.name() + "_hess", f.name_in(), out);
  }

  /** Inputs by name, 0 for the others (as checkfunction with a dict in helpers.py) */
  private static List<DM> inputsOf(Function f, Map<String, DM> d) {
    List<DM> r = new ArrayList<>();
    for (long i = 0; i < f.n_in(); i++) r.add(new DM(0.0));
    List<String> ns = f.name_in();
    for (Map.Entry<String, DM> e : d.entrySet()) r.set(ns.indexOf(e.getKey()), e.getValue());
    return r;
  }

  /** Options of checkfunction in helpers.py */
  private static final class CF {
    boolean jacobian = true, hessian = true;
    int evals = 2, digits = 9;
    Integer digitsSens = null;
    String failmessage = "";

    CF jacobian(boolean b) { jacobian = b; return this; }
    CF hessian(boolean b) { hessian = b; return this; }
    CF evals(int e) { evals = e; return this; }
    CF digits(int d) { digits = d; return this; }
    CF digitsSens(int d) { digitsSens = d; return this; }
    CF fail(String s) { failmessage = s; return this; }
  }

  private static CF cf() {
    return new CF();
  }

  /** checkfunction of helpers.py (gradient/fwd/adj/sens_der flags have no effect there either) */
  private void checkfun(Function trial, Function solution, List<DM> inputs, CF o) {
    int digitsSens = o.digitsSens == null ? o.digits : o.digitsSens;
    for (int rep = 0; rep < 2; rep++) {
      List<DM> to = trial.call(inputs), so = solution.call(inputs);
      assertEquals(trial.n_out(), solution.n_out(), o.failmessage + ": n_out");
      assertEquals(trial.n_in(), solution.n_in(), o.failmessage + ": n_in");
      for (int i = 0; i < trial.n_out(); i++) {
        if (trial.sparsity_out(i).is_empty() || solution.sparsity_out(i).is_empty()) continue;
        chk(to.get(i), so.get(i), o.digits, o.failmessage + ": output(" + i + ": '" + trial.name_out(i) + "')");
      }
    }
    if (o.jacobian) {
      for (long i = 0; i < trial.n_in(); i++) {
        if (trial.sparsity_in(i).is_empty() || solution.sparsity_in(i).is_empty()) continue;
        for (long j = 0; j < trial.n_out(); j++) {
          Function tj = jacobianOld(trial, i, j), sj = jacobianOld(solution, i, j);
          assertEquals(trial.n_in(), tj.n_in());
          assertEquals(trial.n_out() + 1, tj.n_out());
          assertEquals(solution.n_in(), sj.n_in());
          assertEquals(solution.n_out() + 1, sj.n_out());
          checkfun(tj, sj, inputs, cf().jacobian(false).hessian(false).evals(0).digits(digitsSens)
              .fail("(" + o.failmessage + ").jacobian_old(" + i + "," + j + ")"));
        }
      }
    }
    if (o.hessian) {
      for (long i = 0; i < trial.n_in(); i++) {
        if (trial.sparsity_in(i).is_empty() || solution.sparsity_in(i).is_empty()) continue;
        for (long j = 0; j < trial.n_out(); j++) {
          if (!(trial.sparsity_out(j).is_scalar() && solution.sparsity_out(j).is_scalar())) continue;
          Function th = hessianOld(trial, i, j), sh = hessianOld(solution, i, j);
          assertEquals(trial.n_in(), th.n_in());
          assertEquals(trial.n_out() + 2, th.n_out());
          assertEquals(solution.n_in(), sh.n_in());
          assertEquals(solution.n_out() + 2, sh.n_out());
          checkfun(th, sh, inputs, cf().jacobian(false).hessian(false).evals(0).digits(digitsSens)
              .fail("(" + o.failmessage + ").hessian_old(" + i + "," + j + ")"));
        }
      }
    }
    if (o.evals > 0) checkEvals(trial, solution, inputs, o.evals, digitsSens, o.failmessage);
  }

  private static Sparsity spmod(int k, Sparsity sp) {
    DM ret = new DM(sp);
    if (k == 1 && ret.numel() > 0) ret.set(0, 0, new DM(1, 1));
    if (k == 2 && ret.nnz() > 0) {
      long[] row = ret.sparsity().row(), col = ret.sparsity().get_col();
      ret.set(row[row.length - 1], col[col.length - 1], new DM(1, 1));
    }
    return ret.sparsity();
  }

  private static boolean[] cat(boolean[]... a) {
    int n = 0;
    for (boolean[] e : a) n += e.length;
    boolean[] r = new boolean[n];
    int k = 0;
    for (boolean[] e : a) for (boolean b : e) r[k++] = b;
    return r;
  }

  private static List<DM> seeded(List<DM> inputs, Function vf, long from) {
    List<DM> r = new ArrayList<>(inputs);
    for (long i = from; i < vf.n_in(); i++) {
      Random rnd = new Random(i);
      double[] v = new double[(int) vf.nnz_in(i)];
      for (int k = 0; k < v.length; k++) v[k] = rnd.nextDouble();
      r.add(new DM(vf.sparsity_in(i), new DM(v)));
    }
    return r;
  }

  private static Function seedFunction(Function f, List<Sparsity> fsp, List<Sparsity> asp, List<Sparsity> isp) {
    List<MX> fseeds = new ArrayList<>(), aseeds = new ArrayList<>(), ins = new ArrayList<>();
    for (Sparsity s : fsp) fseeds.add(MX.sym("f", s));
    for (Sparsity s : asp) aseeds.add(MX.sym("a", s));
    for (Sparsity s : isp) ins.add(MX.sym("i", s));
    List<MX> res = f.call(ins, !f.is_a("SXFunction"));
    Map<String, Object> o = opts("helper_options", opts("is_diff_in", f.is_diff_in(), "is_diff_out", f.is_diff_out()));
    List<MX> fwd = ca.forward(res, ins, List.of(fseeds), o).get(0);
    List<MX> adj = ca.reverse(res, ins, List.of(aseeds), o).get(0);
    List<MX> vin = new ArrayList<>(ins), vout = new ArrayList<>(res);
    vin.addAll(fseeds);
    vin.addAll(aseeds);
    vout.addAll(fwd);
    vout.addAll(adj);
    return new Function("vf", vin, vout,
        opts("is_diff_in", cat(f.is_diff_in(), f.is_diff_in(), f.is_diff_out()),
            "is_diff_out", cat(f.is_diff_out(), f.is_diff_out(), f.is_diff_in())));
  }

  private void checkEvals(Function trial, Function solution, List<DM> inputs, int evals, int digits, String fail) {
    Map<String, List<List<DM>>> storage = new LinkedHashMap<>(), storage2 = new LinkedHashMap<>();
    Function vfRef = null;
    for (Function f : List.of(solution, trial)) {
      for (int m1 = 0; m1 < 3; m1++)
        for (int m2 = 0; m2 < 3; m2++) {
          List<Sparsity> fsp = new ArrayList<>(), asp = new ArrayList<>(), isp = new ArrayList<>();
          for (long i = 0; i < f.n_in(); i++) {
            fsp.add(spmod(m1, f.sparsity_in(i)));
            isp.add(f.sparsity_in(i));
          }
          for (long i = 0; i < f.n_out(); i++) asp.add(spmod(m2, f.sparsity_out(i)));
          Function vf = seedFunction(f, fsp, asp, isp);
          String key = m1 + "," + m2;
          storage.computeIfAbsent(key, k -> new ArrayList<>()).add(vf.call(seeded(inputs, vf, f.n_in())));
          if (vfRef == null) vfRef = vf;
          if (evals > 1) {
            // python repeats this identically for 9 unused (spmod_2, spmod2_2) pairs
            List<Sparsity> f2 = new ArrayList<>(), a2 = new ArrayList<>(), i2 = new ArrayList<>();
            for (long i = 0; i < vf.n_in(); i++) {
              f2.add(vfRef.sparsity_in(i));
              i2.add(vfRef.sparsity_in(i));
            }
            for (long i = 0; i < vf.n_out(); i++) a2.add(vfRef.sparsity_out(i));
            Function vf2 = seedFunction(vf, f2, a2, i2);
            storage2.computeIfAbsent(key, k -> new ArrayList<>()).add(vf2.call(seeded(inputs, vf2, f.n_in())));
          }
        }
    }
    List<Map<String, List<List<DM>>>> stores = evals > 1 ? List.of(storage, storage2) : List.of(storage);
    for (int s = 0; s < stores.size(); s++)
      for (List<List<DM>> st : stores.get(s).values())
        for (int i = 0; i < st.size() - 1; i++)
          for (int k = 0; k < st.get(0).size(); k++) {
            DM a = st.get(0).get(k), b = st.get(i + 1).get(k);
            if (b.numel() == 0 && ca.sparsify(a).nnz() == 0) continue;
            if (a.numel() == 0 && ca.sparsify(b).nnz() == 0) continue;
            chk(a, b, digits, (s == 0 ? "first-order" : "second-order") + ", output(" + k + ")" + fail);
          }
  }

  private void checkfunctionLight(Function trial, Function solution, List<DM> inputs) {
    checkfun(trial, solution, inputs, cf().jacobian(false).hessian(false).evals(0));
  }

  /** check_serialize of helpers.py, without the python-side serialization scheme validation */
  private static void checkSerialize(Function f, List<DM> inputs) {
    Function f2 = Function.deserialize(f.serialize(opts("debug", true)));
    List<DM> a = f.call(inputs), b = f2.call(inputs);
    for (int i = 0; i < f.n_out(); i++) chk(a.get(i), b.get(i), 16);
  }

  private static void checkPure(Function f, List<DM> inputs) {
    List<DM> a = f.call(inputs), b = f.call(inputs);
    for (int i = 0; i < f.n_out(); i++) chk(a.get(i), b.get(i), 16);
  }

  private void checkThreadSafety(Function f, List<DM> inputs) {
    checkfunctionLight(f.map(20, "thread", 2), f.map(20, "thread"), inputs);
  }

  private static int cgCounter = 0;

  /** check_codegen of helpers.py: compiles with gcc only when slow tests run */
  private static void checkCodegen(Function f, List<DM> inputs, int digits) {
    if (!RUN_SLOW) return;
    try {
      Path dir = Files.createTempDirectory("casadi_java_cg");
      String name = "codegen_" + Long.toHexString(System.nanoTime()) + "_" + (cgCounter++);
      CodeGenerator cg = new CodeGenerator(name);
      cg.add(f);
      cg.generate(dir.toString() + File.separator);
      String inc = GlobalOptions.getCasadiIncludePath(), lib = GlobalOptions.getCasadiPath();
      String so = dir.resolve(name + ".so").toString();
      List<String> cmd = new ArrayList<>(Arrays.asList("gcc", "-pedantic", "-std=c89", "-fPIC", "-shared", "-Wall",
          "-Werror", "-Wextra", "-I" + inc, "-I" + inc + "/highs", "-Wno-dangling-pointer",
          "-Wno-unknown-warning-option", "-Wno-unknown-pragmas", "-Wno-long-long", "-Wno-unused-parameter", "-O3",
          dir.resolve(name + ".c").toString(), "-o", so, "-L" + lib, "-Wl,-rpath," + lib));
      Process p = new ProcessBuilder(cmd).inheritIO().start();
      assertEquals(0, p.waitFor(), "gcc failed");
      Function f2 = ca.external(f.name(), so);
      List<DM> a = f.call(inputs), b = f2.call(inputs);
      for (int i = 0; i < f.n_out(); i++) chk(a.get(i), b.get(i), digits);
    } catch (java.io.IOException | InterruptedException e) {
      throw new RuntimeException(e);
    }
  }

  private static Function bounceTransition(long nz) {
    SX x = SX.sym("x", 2), z = SX.sym("z", nz), p = SX.sym("p");
    SX index = SX.sym("index"), t = SX.sym("t"), u = new SX(0, 1);
    return new Function("transition", List.of(index, t, x, z, p, u), List.of(vs(x.get(0), p.unaryMinus().times(x.get(1))), z),
        List.of("index", "t", "x", "z", "p", "u"), List.of("post_x", "post_z"));
  }

  private static Function bounceTransition() {
    return bounceTransition(0);
  }

  /** Central differences of the first output of G at v0 */
  private static DM fdJacobian(Function G, DM v0, double h) {
    List<DM> cols = new ArrayList<>();
    for (DM e : ca.horzsplit(DM.eye(v0.numel()))) {
      DM a = G.call(List.of(v0.plus(e.times(h)))).get(0), b = G.call(List.of(v0.minus(e.times(h)))).get(0);
      cols.add(a.minus(b).div(2 * h));
    }
    return ca.horzcat(cols);
  }

  private static DM fdJacobian(Function G, DM v0) {
    return fdJacobian(G, v0, 1e-6);
  }

  // ---------------------------------------------------------------- tests

  private Function integrator;
  private Function qe;
  private MX qend, q0, par;
  private Map<String, SX> dae;
  private final double tendN = 2.3, q0N = 7.1, pN = 2;

  @BeforeEach
  public void setUp() {
    SX t = SX.sym("t"), x = SX.sym("x"), p = SX.sym("p");
    dae = sxd("t", t, "x", x, "p", p, "ode", x.div(p).times(t.power(2)));
    integrator = ca.integrator("integrator", "cvodes", dae, 0, 2.3,
        opts("reltol", 1e-15, "abstol", 1e-15, "fsens_err_con", true));
    q0 = MX.sym("q0");
    par = MX.sym("p");
    qend = integrator.call(mxd("x0", q0, "p", par)).get("xf");
    qe = new Function("qe", List.of(q0, par), List.of(qend));
  }

  @Test
  public void test_full() {
    Assumptions.assumeTrue(RUN_SLOW);
    SX t = SX.sym("t"), q = SX.sym("q"), p = SX.sym("p");
    Map<String, SX> dae = sxd("t", t, "x", q, "p", p, "ode", q.div(p).times(t.power(2)));
    Function integrator = ca.integrator("integrator", "cvodes", dae, 0, 2.3,
        opts("reltol", 1e-15, "abstol", 1e-15, "fsens_err_con", true));
    assertInException("No stats available", () -> integrator.stats());
    double tf = 2.3;
    Function solution = new Function("solution",
        sxd("x0", q, "p", p, "xf", q.times(ca.exp(new SX(Math.pow(tf, 3)).div(p.times(3))))),
        ca.integrator_in(), ca.integrator_out());
    checkfun(integrator, solution, inputsOf(integrator, dmd("x0", 0.3, "p", 0.7)), cf().digits(6));
  }

  @Test
  public void test_tools_trivial() {
    SX x = SX.sym("x"), p = SX.sym("p", 0), tf = SX.sym("tf");
    Function f = new Function("f", List.of(x, p), List.of(x));
    for (Function integrator : List.of(ca.simpleRK(f), ca.simpleIRK(f), ca.simpleIntegrator(f))) {
      Function solution = new Function("solution", List.of(x, p, tf), List.of(x.times(ca.exp(tf))),
          List.of("x0", "p", "h"), List.of("xf"));
      checkfun(integrator, solution, inputsOf(integrator, dmd("x0", 1, "h", 1)), cf().digits(3));
    }
  }

  @Test
  public void test_tools() {
    Assumptions.assumeTrue(RUN_SLOW);
    SX t = SX.sym("t"), q = SX.sym("q"), p = SX.sym("p"), t0 = SX.sym("t0"), q0 = SX.sym("q0");
    Function f = new Function("f", List.of(vs(q, t), p), List.of(vs(q.div(p).times(t.power(2)), 1)));
    for (Function integrator : List.of(ca.simpleRK(f, 500), ca.simpleIRK(f, 50), ca.simpleIntegrator(f))) {
      SX e = t0.plus(t).power(3).minus(t0.power(3)).div(p.times(3));
      Function solution = new Function("solver",
          sxd("x0", vs(q0, t0), "p", p, "h", t, "xf", vs(q0.times(ca.exp(e)), t0.plus(t))),
          integrator.name_in(), integrator.name_out());
      checkfun(integrator, solution, inputsOf(integrator, dmd("x0", dm(0.3, 0), "p", 0.7, "h", 1)), cf().digits(3));
    }
  }

  /** One row of the checks() generators */
  private static final class Chk {
    final String feature;
    final Map<String, SX> din, dout, solutionin, solution;
    final Map<String, DM> point;
    final double t0, tf;

    Chk(String feature, Map<String, SX> din, Map<String, SX> dout, Map<String, SX> solutionin,
        Map<String, SX> solution, Map<String, DM> point, double t0, double tf) {
      this.feature = feature;
      this.din = din;
      this.dout = dout;
      this.solutionin = solutionin;
      this.solution = solution;
      this.point = point;
      this.t0 = t0;
      this.tf = tf;
    }
  }

  /** variations() of the python tests: an ode also as a dae */
  private static List<Chk> variations(Chk c, String zname) {
    List<Chk> r = new ArrayList<>();
    r.add(c);
    if (c.feature.equals("ode")) {
      Map<String, SX> din = new LinkedHashMap<>(c.din), dout = new LinkedHashMap<>(c.dout);
      SX z = SX.sym(zname, c.din.get("x").sparsity().size1(), c.din.get("x").sparsity().size2());
      din.put("z", z);
      dout.put("ode", z);
      dout.put("alg", c.dout.get("ode").minus(z).times(-0.8));
      r.add(new Chk("dae", din, dout, c.solutionin, c.solution, c.point, c.t0, c.tf));
    }
    return r;
  }

  private static Function solutionFunction(Chk c, SX tstart, SX tend) {
    for (String k : new ArrayList<>(c.solution.keySet()))
      c.solution.put(k, ca.substitute(c.solution.get(k), vs(tstart, tend), vs(c.t0, c.tf)));
    Map<String, SX> sol = new LinkedHashMap<>(c.solutionin);
    sol.putAll(c.solution);
    return new Function("fs", sol, ca.integrator_in(), ca.integrator_out());
  }

  @Test
  public void test_jac() {
    message("Test exact jacobian #536");
    DM.set_precision(18);
    try {
      SX tstart = SX.sym("tstart"), tend = SX.sym("tend");
      List<Intg> integrators = List.of(
          new Intg("idas", List.of("dae", "ode"), opts("abstol", 1e-9, "reltol", 1e-9, "fsens_err_con", true,
              "calc_ic", true, "calc_icB", true)),
          new Intg("cvodes", List.of("ode"), opts("abstol", 1e-5, "reltol", 1e-5, "fsens_err_con", false,
              "quad_err_con", false)));
      SX x = SX.sym("x"), t = SX.sym("t");
      Chk check = new Chk("ode", sxd("x", x, "t", t), sxd("ode", x), sxd("x0", x),
          sxd("xf", x.times(ca.exp(tend))), dmd("x0", 1), 0, 0.9995);
      for (Chk c : variations(check, "x")) {
        for (Intg I : integrators) {
          message(I.name);
          SX xx = SX.sym("x");
          ca.integrator("dummyIntegrator", I.name, sxd("x", xx, "ode", xx));
          if (!I.features.contains(c.feature)) continue;
          Map<String, SX> dae = new LinkedHashMap<>(c.din);
          dae.putAll(c.dout);
          Function fs = solutionFunction(c, tstart, tend);
          for (String scheme : List.of("direct", "gmres", "bcgstab")) {
            message("f_options: newton_scheme " + scheme);
            Function integrator = ca.integrator("integrator", I.name, dae, c.t0, c.tf,
                with(I.options, "newton_scheme", scheme));
            Map<String, DM> in = new LinkedHashMap<>();
            for (Map.Entry<String, DM> e : c.point.entrySet())
              if (!fs.sparsity_in(e.getKey()).is_empty()) in.put(e.getKey(), e.getValue());
            Map<String, DM> io = integrator.call(in), fo = fs.call(in);
            message("res=" + io.get("xf").minus(fo.get("xf")) + " " + fo.get("xf"));
            Map<String, Object> stats = integrator.stats();
            message(stats.toString());
            assertTrue(num(stats.get("nstepsB")) < 2500);
            assertTrue(num(stats.get("nsteps")) < 1500);
            assertTrue(num(stats.get("nstepsB")) < 2500);
            assertTrue(num(stats.get("nlinsetups")) < 100);
            assertTrue(num(stats.get("nlinsetupsB")) < 250);
          }
        }
      }
    } finally {
      DM.set_precision(6);
    }
  }

  private static Chk nonlinCheck(SX tstart, SX tend) {
    double[] A = {1, 0.1};
    double p0 = 1.13;
    SX q = SX.sym("y", 2, 1), y0 = q.get(0), yc0 = q.get(1), p = SX.sym("p", 1, 1);
    SX arg = ca.atan(yc0.div(ca.sqrt(p))).plus(ca.sqrt(p).times(tend.minus(tstart)));
    SX s1 = y0.times(2).minus(ca.log(yc0.power(2).div(p).plus(1))).div(2).minus(ca.log(ca.cos(arg)));
    SX s2 = ca.sqrt(p).times(ca.tan(arg));
    return new Chk("ode", sxd("x", q, "p", p), sxd("ode", vs(q.get(1), p.get(0).plus(q.get(1).power(2)))),
        sxd("x0", q, "p", p), sxd("xf", vs(s1, s2)), dmd("x0", A, "p", p0), 0, 0.4);
  }

  @Test
  public void test_lsolvers() {
    message("Test different linear solvers");
    SX tstart = SX.sym("tstart"), tend = SX.sym("tend");
    List<Intg> integrators = List.of(
        new Intg("idas", List.of("dae", "ode"), opts("abstol", 1e-9, "reltol", 1e-9, "fsens_err_con", true,
            "calc_ic", true, "calc_icB", true)),
        new Intg("cvodes", List.of("ode"), opts("abstol", 1e-15, "reltol", 1e-15, "fsens_err_con", true,
            "quad_err_con", false)));
    SX x = SX.sym("x"), p = SX.sym("p"), z = SX.sym("z");
    List<Chk> checks = new ArrayList<>();
    checks.add(new Chk("dae", sxd("x", x, "z", z), sxd("alg", x.minus(z), "ode", z), sxd("x0", x, "p", p),
        sxd("xf", x.times(ca.exp(tend.minus(tstart)))), dmd("x0", 7.1, "p", 2), 0.2, 2.3));
    if (RUN_SLOW) {
      checks.add(new Chk("dae", sxd("x", x, "z", z), sxd("alg", x.minus(z), "ode", z), sxd("x0", x, "p", p),
          sxd("xf", x.times(ca.exp(tend.minus(tstart)))), dmd("x0", 7.1, "p", 2), 0.2, 2.3));
      checks.add(new Chk("ode", sxd("x", x), sxd("ode", x), sxd("x0", x, "p", p),
          sxd("xf", x.times(ca.exp(tend.minus(tstart)))), dmd("x0", 7.1, "p", 2), 0.2, 2.3));
      checks.add(nonlinCheck(tstart, tend));
    }
    for (Chk c : checks) {
      for (Intg I : integrators) {
        message(I.name);
        SX xx = SX.sym("x");
        ca.integrator("dummyIntegrator", I.name, sxd("x", xx, "ode", xx));
        if (!I.features.contains(c.feature)) continue;
        Map<String, SX> dae = new LinkedHashMap<>(c.din);
        dae.putAll(c.dout);
        Function fs = solutionFunction(c, tstart, tend);
        for (String scheme : List.of("direct", "gmres", "bcgstab")) {
          String msg = "f_options: newton_scheme " + scheme;
          message(msg);
          Function integrator = ca.integrator("integrator", I.name, dae, c.t0, c.tf,
              with(I.options, "newton_scheme", scheme));
          Map<String, DM> in = new LinkedHashMap<>();
          for (Map.Entry<String, DM> e : c.point.entrySet())
            if (!integrator.sparsity_in(e.getKey()).is_empty()) in.put(e.getKey(), e.getValue());
          List<DM> inputs = inputsOf(integrator, in);
          checkfun(integrator, fs, inputs, cf().hessian(false).evals(0).digits(4).digitsSens(4).fail(msg));
          checkPure(integrator, inputs);
          checkSerialize(integrator, inputs);
        }
      }
    }
  }

  @Test
  public void test_X() {
    message("Extensive integrator tests");
    SX tstart = SX.sym("tstart"), tend = SX.sym("tend");
    for (Intg I : integrators) {
      message(I.name);
      double x0 = q0N, p_ = pN;
      SX t = SX.sym("t"), x = SX.sym("x"), p = SX.sym("p"), z = SX.sym("z");
      List<Chk> checks = new ArrayList<>();
      double tf = tendN;
      Object[][] rows = {
          {"ode", sxd("x", x), sxd("ode", 0), sxd("xf", x)},
          {"ode", sxd("x", x), sxd("ode", 1), sxd("xf", x.plus(tend.minus(tstart)))},
          {"ode", sxd("x", x), sxd("ode", x), sxd("xf", x.times(ca.exp(tend.minus(tstart))))},
          {"ode", sxd("x", x, "t", t), sxd("ode", t), sxd("xf", x.plus(tend.power(2).div(2).minus(tstart.power(2).div(2))))},
          {"ode", sxd("x", x, "t", t), sxd("ode", x.times(t)),
              sxd("xf", x.times(ca.exp(tend.power(2).div(2).minus(tstart.power(2).div(2)))))},
          {"ode", sxd("x", x, "p", p), sxd("ode", x.div(p)), sxd("xf", x.times(ca.exp(tend.minus(tstart).div(p))))},
      };
      Object[][] slowRows = {
          {"ode", sxd("x", x), sxd("ode", x, "quad", 0), sxd("qf", 0)},
          {"ode", sxd("x", x), sxd("ode", x, "quad", 1), sxd("qf", tend.minus(tstart))},
          {"ode", sxd("x", x), sxd("ode", 0, "quad", x), sxd("qf", x.times(tend.minus(tstart)))},
          {"ode", sxd("x", x), sxd("ode", x, "quad", x), sxd("qf", x.times(ca.exp(tend.minus(tstart)).minus(1)))},
          {"ode", sxd("x", x, "t", t), sxd("ode", x, "quad", t), sxd("qf", tend.power(2).div(2).minus(tstart.power(2).div(2)))},
          {"ode", sxd("x", x, "t", t), sxd("ode", x, "quad", x.times(t)),
              sxd("qf", x.times(ca.exp(tend.minus(tstart)).times(tend.minus(1)).minus(tstart.minus(1))))},
          {"ode", sxd("x", x, "p", p), sxd("ode", x, "quad", x.div(p)),
              sxd("qf", x.times(ca.exp(tend.minus(tstart)).minus(1)).div(p))},
          {"dae", sxd("x", x, "z", z), sxd("ode", z, "alg", z.minus(x).times(-0.8), "quad", z),
              sxd("qf", x.times(ca.exp(tend.minus(tstart)).minus(1)))},
      };
      List<Object[]> all = new ArrayList<>(Arrays.asList(rows));
      if (RUN_SLOW) all.addAll(Arrays.asList(slowRows));
      for (Object[] r : all) {
        @SuppressWarnings("unchecked")
        Chk c = new Chk((String) r[0], (Map<String, SX>) r[1], (Map<String, SX>) r[2], sxd("x0", x, "p", p),
            (Map<String, SX>) r[3], dmd("x0", x0, "p", p_), 0.2, tf);
        checks.add(c);
      }
      if (RUN_SLOW) checks.add(nonlinCheck(tstart, tend));
      for (Chk tt : checks) {
        for (Chk c : variations(tt, "z")) {
          if (!I.features.contains(c.feature)) continue;
          String msg = I.name + ": " + c.din + " => " + c.dout + ", explicit (" + c.solution + ") tstart = " + c.t0;
          message(msg);
          Map<String, SX> dae = new LinkedHashMap<>(c.din);
          dae.putAll(c.dout);
          Function fs = solutionFunction(c, tstart, tend);
          Map<String, Object> o = new LinkedHashMap<>(I.options);
          if (I.name.equals("cvodes") || I.name.equals("idas")) {
            o.put("abstol", 1e-9);
            o.put("reltol", 1e-9);
          }
          if (I.name.equals("idas")) {
            o.put("init_xdot", c.point.get("x0").nonzeros());
            o.put("calc_icB", true);
            Map<String, Object> aug = new LinkedHashMap<>();
            aug.put("init_xdot", null);
            aug.put("abstol", 1e-9);
            aug.put("reltol", 1e-9);
            o.put("augmented_options", aug);
          }
          Function integrator = ca.integrator("integrator", I.name, dae, c.t0, c.tf, o);
          Map<String, DM> in = new LinkedHashMap<>();
          for (Map.Entry<String, DM> e : c.point.entrySet())
            if (!integrator.sparsity_in(e.getKey()).is_empty()) in.put(e.getKey(), e.getValue());
          List<DM> inputs = inputsOf(integrator, in);
          checkfun(integrator, fs, inputs, cf().hessian(false).evals(0).digits(4).digitsSens(4).fail(msg));
          if (!o.toString().contains("kinsol")) {
            message(integrator.toString());
            checkSerialize(integrator, inputs);
          }
        }
      }
    }
  }

  @Test
  public void test_multipoint() {
    for (Intg I : integrators) {
      message(I.name);
      if (I.name.equals("cvodes") || I.name.equals("idas")) continue;
      MX x = MX.sym("x"), p = MX.sym("p"), z = MX.sym("z"), u = MX.sym("u");
      Map<String, MX> dae = I.features.contains("dae")
          ? mxd("x", x, "z", z, "p", p, "u", u, "ode", z.plus(p.times(u)), "alg", z.times(ca.cos(z)).minus(x), "quad", x.power(2))
          : mxd("x", x, "p", p, "u", u, "ode", ca.cos(x.plus(p.times(u))), "quad", x.power(2));
      Function integratorA = ca.integrator("integrator", I.name, dae, 0, 0.1, I.options);
      Function integratorB = ca.integrator("integrator", I.name, dae, 0.1, 0.2, I.options);
      Function integratorI = ca.integrator("integrator", I.name, dae, 0, new double[] {0.1, 0.2}, I.options);
      Map<String, DM> argTest = dmd("x0", 0, "p", 0.15, "u", dm(1, 1.1).T());
      Map<String, MX> args = integratorI.convert_in(integratorI.mx_in());
      List<MX> uSplit = ca.horzsplit_n(args.get("u"), 2);
      Map<String, MX> args1 = new LinkedHashMap<>(args), args2 = new LinkedHashMap<>(args);
      args1.put("u", uSplit.get(0));
      args2.put("u", uSplit.get(1));
      Map<String, MX> res1 = integratorA.call(args1);
      args2.put("x0", res1.get("xf"));
      Map<String, MX> res2 = integratorB.call(args2);
      res2.put("qf", ca.horzcat(List.of(res1.get("qf"), res1.get("qf").plus(res2.get("qf")))));
      res2.put("xf", ca.horzcat(List.of(res1.get("xf"), res2.get("xf"))));
      res2.put("zf", ca.horzcat(List.of(res1.get("zf"), res2.get("zf"))));
      res2.putAll(args);
      Function If = new Function("I", res2, integratorI.name_in(), integratorI.name_out());
      If.call(argTest);
      integratorI.call(argTest);
      checkfun(integratorI, If, inputsOf(integratorI, dmd("x0", 0, "p", 0.15, "u", dm(1, 1.1).T())), cf());
    }
  }

  @Test
  public void test_eval2() {
    message("CVodes integration: evaluation with Function indirection");
    Function qe = new Function("qe", List.of(q0, par), List.of(qend.get(0)));
    Function f = new Function("f", List.of(q0), List.of((MX) qe.invoke(q0, new MX(pN))));
    DM fOut = f.call(List.of(dm(q0N))).get(0);
    almost(fOut.to_double(), q0N * Math.exp(Math.pow(tendN, 3) / (3 * pN)), 9, "Evaluation output mismatch");
  }

  @Test
  public void test_issue92b() {
    message("regression check for issue 92");
    SX t = SX.sym("t"), x = SX.sym("x"), y = SX.sym("y");
    Map<String, SX> dae = sxd("t", t, "x", vs(x, y), "ode", vs(x, x.times(1 + 1e-9)));
    Function integrator = ca.integrator("integrator", "cvodes", dae, opts("fsens_err_con", true));
    Map<String, DM> out = integrator.call(dmd("x0", dm(1, 0)));
    message(out.get("xf").toString());
  }

  @Test
  public void test_issue92() {
    message("regression check for issue 92");
    SX t = SX.sym("t"), x = SX.sym("x");
    MX var = MX.sym("var", 2, 1);
    SX q = vs(x, SX.sym("problem"));
    SX dq = vs(x, x);
    Function integrator = ca.integrator("integrator", "cvodes", sxd("t", t, "x", q, "ode", dq),
        opts("fsens_err_con", true, "reltol", 1e-12));
    MX qend = integrator.call(mxd("x0", var)).get("xf");
    Function f = new Function("f", List.of(var), List.of(qend.get(0)));
    Function J = jacobianOld(f, 0, 0);
    List<DM> Jout = J.call(List.of(dm(1, 0)));
    message("jac=" + Jout.get(0).get_nz(false, 0).minus(Math.exp(1)));
    almost(Jout.get(0).get(0, 0).to_double(), Math.exp(1), 5, "Evaluation output mismatch");
  }

  @Test
  public void test_eval() {
    message("CVodes integration: evaluation");
    DM out = qe.call(List.of(dm(q0N), dm(pN))).get(0);
    almost(out.to_double(), q0N * Math.exp(Math.pow(tendN, 3) / (3 * pN)), 9, "Evaluation output mismatch");
  }

  @Test
  public void test_jac1() {
    message("CVodes integration: jacobian to q0");
    Function J = jacobianOld(qe, 0, 0);
    DM out = J.call(List.of(dm(q0N), dm(pN))).get(0);
    almost(out.to_double(), Math.exp(Math.pow(tendN, 3) / (3 * pN)), 9, "Evaluation output mismatch");
  }

  @Test
  public void test_jac2() {
    message("CVodes integration: jacobian to p");
    Function J = jacobianOld(qe, 1, 0);
    DM out = J.call(List.of(dm(q0N), dm(pN))).get(0);
    double te = tendN, q0 = q0N, p = pN;
    almost(out.to_double(), -(q0 * Math.pow(te, 3) * Math.exp(Math.pow(te, 3) / (3 * p))) / (3 * p * p), 9,
        "Evaluation output mismatch");
  }

  @Test
  public void test_bug_repeat() {
    message("Bug that appears when rhs contains repeats");
    DM A = dm(1, 0.1, 1);
    double p0 = 1.13, te = 0.4;
    SX t = SX.sym("t"), q = SX.sym("y", 3, 1), p = SX.sym("p");
    SX dh = p.plus(q.get(0).power(2));
    Map<String, Object> o = opts("reltol", 1e-15, "abstol", 1e-15, "fsens_err_con", true, "steps_per_checkpoint", 10000);
    Function integrator = ca.integrator("integrator", "cvodes", sxd("x", q, "p", p, "t", t, "ode", vs(dh, q.get(0), dh)),
        0, te, o);
    MX q0 = MX.sym("q0", 3, 1), par = MX.sym("p", 1, 1);
    MX qend = integrator.call(mxd("x0", q0, "p", par)).get("xf");
    Function qe = new Function("qe", List.of(q0, par), List.of(qend));
    DM outA = jacobianOld(qe, 0, 0).call(List.of(A, dm(p0))).get(0);
    integrator = ca.integrator("integrator", "cvodes",
        sxd("x", q, "p", p, "t", t, "ode", vs(dh, q.get(0), dh.times(1 + 1e-9))), 0, te, o);
    q0 = MX.sym("q0", 3, 1);
    par = MX.sym("p", 1, 1);
    qend = integrator.call(mxd("x0", q0, "p", par)).get("xf");
    qe = new Function("qe", List.of(q0, par), List.of(qend));
    DM outB = jacobianOld(qe, 0, 0).call(List.of(A, dm(p0))).get(0);
    message(outA.minus(outB).toString());
  }

  private double hessRef() {
    double te = tendN, q0 = q0N, p = pN;
    return (q0 * Math.pow(te, 6) * Math.exp(Math.pow(te, 3) / (3 * p))) / (9 * Math.pow(p, 4))
        + (2 * q0 * Math.pow(te, 3) * Math.exp(Math.pow(te, 3) / (3 * p))) / (3 * Math.pow(p, 3));
  }

  @Test
  public void test_hess3() {
    message("CVodes integration: hessian to p: Jacobian of integrator.jacobian");
    Function J = jacobianOld(integrator, integrator.index_in("p"), integrator.index_out("xf"));
    Function H = jacobianOld(J, J.index_in("p"), 0);
    Map<String, DM> out = H.call(dmd("x0", q0N, "p", pN));
    almost(out.get("jac_jac_xf_p_p").to_double(), hessRef(), 9, "Evaluation output mismatch");
  }

  @Test
  public void test_hess4() {
    message("CVodes integration: hessian to p: Jacobian of integrator.jacobian indirect");
    Function J = jacobianOld(integrator, integrator.index_in("p"), integrator.index_out("xf"));
    MX q0 = MX.sym("q0"), p = MX.sym("p");
    Map<String, MX> Ji = J.call(mxd("x0", q0, "p", p));
    Ji.put("q0", q0);
    Ji.put("p", p);
    Function Jif = new Function("Ji", Ji, List.of("q0", "p"), J.name_out());
    Function H = jacobianOld(Jif, 1, 0);
    DM out = H.call(List.of(dm(q0N), dm(pN))).get(0);
    almost(out.to_double(), hessRef(), 9, "Evaluation output mismatch");
  }

  @Test
  public void test_hess5() {
    message("CVodes integration: hessian to p in an MX tree");
    MX q0 = MX.sym("q0"), p = MX.sym("p");
    Map<String, MX> sol = integrator.call(mxd("x0", q0, "p", p));
    sol.put("q0", q0);
    sol.put("p", p);
    Function JT = new Function("JT", List.of(q0, p), List.of(ca.jacobian(sol.get("xf"), sol.get("p")).T()));
    message(JT.call(List.of(dm(q0N), dm(pN))).toString());
    Function H = jacobianOld(JT, 1, 0);
    DM out = H.call(List.of(dm(q0N), dm(pN))).get(0);
    almost(out.to_double(), hessRef(), 9, "Evaluation output mismatch");
  }

  @Test
  public void test_hess6() {
    message("CVodes integration: hessian to p in an MX tree");
    MX q0 = MX.sym("q0"), p = MX.sym("p");
    Map<String, MX> sol = integrator.call(mxd("x0", q0, "p", p));
    sol.put("q0", q0);
    sol.put("p", p);
    Function qe = new Function("qe", sol, List.of("q0", "p"), ca.integrator_out());
    Function H = hessianOld(qe, 1, 0);
    DM out = H.call(List.of(dm(q0N), dm(pN))).get(0);
    almost(out.to_double(), hessRef(), 7, "Evaluation output mismatch");
  }

  private static final double[][] B3 = {{1, 2.3, 4}, {-2, 1.3, 4.7}, {-2, 6, 9}};

  @Test
  public void test_glibcbug() {
    message("former glibc error");
    DM A = dm(2.3, 4.3, 7.6);
    double te = 0.7;
    SX t = SX.sym("t"), q = SX.sym("q", 3, 1), p = SX.sym("p", 9, 1);
    Function integrator = ca.integrator("integrator", "cvodes",
        sxd("t", t, "x", q, "p", p, "ode", ca.mtimes(ca.reshape(p, 3, 3), q)), 0, te,
        opts("fsens_err_con", true, "steps_per_checkpoint", 1000));
    MX q0 = MX.sym("q0", 3, 1), par = MX.sym("p", 9, 1);
    integrator.call(mxd("x0", q0, "p", par));
    Function qeJ = jacobianOld(integrator, integrator.index_in("p"), integrator.index_out("xf"));
    MX qe = qeJ.call(mxd("x0", q0, "p", par)).get("jac_xf_p");
    Function qef = new Function("qef", List.of(q0, par), List.of(qe));
    // B.ravel(): row-major
    double[] Br = new double[9];
    for (int i = 0; i < 3; i++) for (int j = 0; j < 3; j++) Br[3 * i + j] = B3[i][j];
    qef.call(List.of(A, dm(Br)));
  }

  /** scipy.linalg.expm replacement: scaling and squaring with a Taylor series */
  private static double[][] expm(double[][] M) {
    int n = M.length;
    double norm = 0;
    for (double[] r : M) {
      double s = 0;
      for (double v : r) s += Math.abs(v);
      norm = Math.max(norm, s);
    }
    int sq = Math.max(0, (int) Math.ceil(Math.log(norm) / Math.log(2)) + 4);
    double sc = Math.pow(2, -sq);
    double[][] X = new double[n][n], term = new double[n][n], R = new double[n][n];
    for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) X[i][j] = M[i][j] * sc;
    for (int i = 0; i < n; i++) { term[i][i] = 1; R[i][i] = 1; }
    for (int k = 1; k < 30; k++) {
      term = matmul(term, X);
      for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) { term[i][j] /= k; R[i][j] += term[i][j]; }
    }
    for (int s = 0; s < sq; s++) R = matmul(R, R);
    return R;
  }

  private static double[][] matmul(double[][] a, double[][] b) {
    int n = a.length, m = b[0].length, l = b.length;
    double[][] r = new double[n][m];
    for (int i = 0; i < n; i++) for (int j = 0; j < m; j++) for (int k = 0; k < l; k++) r[i][j] += a[i][k] * b[k][j];
    return r;
  }

  @Test
  public void test_linear_system() {
    message("Linear ODE");
    DM A = dm(2.3, 4.3, 7.6);
    double te = 0.7;
    double[][] Bte = new double[3][3];
    for (int i = 0; i < 3; i++) for (int j = 0; j < 3; j++) Bte[i][j] = B3[i][j] * te;
    DM Be = new DM(expm(Bte));
    SX t = SX.sym("t"), q = SX.sym("q", 3, 1), p = SX.sym("p", 9, 1);
    Function integrator = ca.integrator("integrator", "cvodes",
        sxd("t", t, "x", q, "p", p, "ode", ca.mtimes(ca.reshape(p, 3, 3), q)), 0, te,
        opts("fsens_err_con", true, "reltol", 1e-15, "abstol", 1e-15, "steps_per_checkpoint", 10000));
    MX q0 = MX.sym("q0", 3, 1), par = MX.sym("p", 9, 1);
    MX qend = integrator.call(mxd("x0", q0, "p", par)).get("xf");
    Function qe = new Function("qe", List.of(q0, par), List.of(qend));
    MX qendJ = jacobianOld(integrator, integrator.index_in("x0"), integrator.index_out("xf"))
        .call(mxd("x0", q0, "p", par)).get("jac_xf_x0");
    Function qeJ = new Function("qeJ", List.of(q0, par), List.of(qendJ));
    MX qendJ2 = jacobianOld(integrator, integrator.index_in("x0"), integrator.index_out("xf"))
        .call(mxd("x0", q0, "p", par)).get("jac_xf_x0");
    Function qeJ2 = new Function("qeJ2", List.of(q0, par), List.of(qendJ2));
    DM B = new DM(B3);
    DM qeOut = qe.call(List.of(A, ca.vec(B))).get(0);
    chk(ca.mtimes(Be, A).div(1e3), qeOut.div(1e3), 10, "jacobian('x0','xf')");
    DM qeJOut = qeJ.call(List.of(A, ca.vec(B))).get(0);
    chk(qeJOut.div(1e3), Be.div(1e3), 10, "jacobian('x0','xf')");
    qeJ2.call(List.of(A, ca.vec(B)));
  }

  @Test
  public void test_mathieu_system() {
    message("Mathieu ODE");
    DM A = dm(0.3, 1.2), B = dm(1.3, 4.3, 2.7);
    double te = 0.7;
    SX t = SX.sym("t"), q = SX.sym("q", 2, 1), p = SX.sym("p", 3, 1);
    SX ode = vs(q.get(1), p.get(0).minus(p.get(1).times(2).times(ca.cos(p.get(2).times(2)))).times(q.get(0)));
    Function integrator = ca.integrator("integrator", "cvodes", sxd("x", q, "p", p, "t", t, "ode", ode), 0, te,
        opts("fsens_err_con", true, "reltol", 1e-15, "abstol", 1e-15, "steps_per_checkpoint", 10000));
    MX q0 = MX.sym("q0", 2, 1), par = MX.sym("p", 3, 1);
    MX qend = integrator.call(mxd("x0", q0, "p", par)).get("xf");
    Function qe = new Function("qe", List.of(q0, par), List.of(qend));
    MX qendJ = jacobianOld(integrator, integrator.index_in("x0"), integrator.index_out("xf"))
        .call(mxd("x0", q0, "p", par)).get("jac_xf_x0");
    new Function("qeJ", List.of(q0, par), List.of(qendJ));
    message(qe.call(List.of(A, B)).get(0).toString());
  }

  @Test
  public void test_nl_system() {
    message("Nonlinear ODE sys");
    DM A = dm(1, 0.1);
    double p0 = 1.13, y0 = 1, yc0 = 0.1, dy0 = 0.1, te = 0.4;
    SX t = SX.sym("t"), q = SX.sym("y", 2, 1), p = SX.sym("p", 1, 1);
    Function integrator = ca.integrator("integrator", "cvodes",
        sxd("x", q, "p", p, "t", t, "ode", vs(q.get(1), p.get(0).plus(q.get(1).power(2)))), 0, te,
        opts("reltol", 1e-15, "abstol", 1e-15, "steps_per_checkpoint", 10000, "fsens_err_con", true));
    MX q0 = MX.sym("q0", 2, 1), par = MX.sym("p", 1, 1);
    MX qend = integrator.call(mxd("x0", q0, "p", par)).get("xf");
    Function qe = new Function("qe", List.of(q0, par), List.of(qend));
    MX qendJ = jacobianOld(integrator, integrator.index_in("x0"), integrator.index_out("xf"))
        .call(mxd("x0", q0, "p", par)).get("jac_xf_x0");
    Function qeJ = new Function("qeJ", List.of(q0, par), List.of(qendJ));
    DM qeOut = qe.call(List.of(A, dm(p0))).get(0);
    message(qeOut.get(0).toString());
    message(qeOut.get(1).toString());
    double sp = Math.sqrt(p0), ta = Math.tan(Math.atan(yc0 / sp) + sp * te);
    almost(qeOut.get(0).to_double(),
        (2 * y0 - Math.log(yc0 * yc0 / p0 + 1)) / 2 - Math.log(Math.cos(Math.atan(yc0 / sp) + sp * te)), 11, "Nonlin ODE");
    almost(qeOut.get(1).to_double(), sp * ta, 11, "Nonlin ODE");
    DM qeJOut = qeJ.call(List.of(A, dm(p0))).get(0);
    double tb = Math.tan(sp * te + Math.atan(dy0 / sp));
    DM Jr = dm(new double[][] {{1, (sp * tb - dy0) / (dy0 * dy0 + p0)}, {0, (p0 * tb * tb + p0) / (dy0 * dy0 + p0)}});
    chk(qeJOut, Jr, 10, "jacobian of Nonlin ODE");
    Function Jf = jacobianOld(qe, 0, 0);
    chk(Jf.call(List.of(A, dm(p0))).get(0), Jr, 10, "Jacobian of Nonlin ODE");
    Jf = jacobianOld(qe, 0, 0);
    chk(Jf.call(List.of(A, dm(p0))).get(0), Jr, 10, "Jacobian of Nonlin ODE");
    double k = te * yc0 * yc0 - yc0 + p0 * te;
    Jr = dm(new double[][] {{(sp * k * ta + yc0 * yc0) / (2 * p0 * yc0 * yc0 + 2 * p0 * p0)},
        {(sp * (k * ta * ta + te * yc0 * yc0 - yc0 + p0 * te) + (yc0 * yc0 + p0) * ta) / (sp * (2 * yc0 * yc0 + 2 * p0))}});
    Jf = jacobianOld(qe, 1, 0);
    chk(Jf.call(List.of(A, dm(p0))).get(0), Jr, 10, "Jacobian of Nonlin ODE");
    Jf = jacobianOld(qe, 1, 0);
    chk(Jf.call(List.of(A, dm(p0))).get(0), Jr, 10, "Jacobian of Nonlin ODE");
    qendJ = jacobianOld(integrator, integrator.index_in("p"), integrator.index_out("xf"))
        .call(mxd("x0", q0, "p", par)).get("jac_xf_p");
    qeJ = new Function("qeJ", List.of(q0, par), List.of(qendJ));
    chk(qeJ.call(List.of(A, dm(p0))).get(0), Jr, 10, "jacobian of Nonlin ODE");
    checkThreadSafety(qeJ, List.of(A, dm(p0)));
    Function qeJf = new Function("qeJf", List.of(q0, par), List.of(ca.vec((MX) qeJ.invoke(q0, par))));
    Function H = jacobianOld(qeJf, 0, 0);
    message(H.call(List.of(A, dm(p0))).get(0).toString());
  }

  @Test
  public void test_missing_symbols() {
    SX x = SX.sym("x"), z = SX.sym("z"), p = SX.sym("p");
    Map<String, SX> dae = sxd("x", x, "z", z, "ode", z.plus(p), "alg", z.times(ca.cos(z)).minus(x));
    assertInException("[p] are free", () -> ca.integrator("F", "idas", dae));
  }

  @Test
  public void test_hessian2D() {
    message("hessian");
    int N = 2;
    DM x0_ = dm(1, 0.1), A_ = dm(new double[][] {{3, 1}, {0.74, 4}});
    SX A = SX.sym("A", N, N), x = SX.sym("x", N);
    Function I = ca.integrator("I", "cvodes", sxd("x", x, "p", ca.vec(A), "ode", ca.mtimes(A, x)),
        opts("fsens_err_con", true, "reltol", 1e-12));
    I.call(dmd("x0", x0_, "p", ca.vec(A_)));
    MX q0 = MX.sym("q0", N), p = MX.sym("p", N * N);
    Map<String, MX> sol = I.call(mxd("x0", q0, "p", p));
    sol.put("q0", q0);
    sol.put("p", p);
    new Function("qe", sol, List.of("q0", "p"), ca.integrator_out());
    Function JT = new Function("JT", List.of(q0, p), List.of(ca.jacobian(sol.get("xf"), sol.get("p")).T()));
    Function H = jacobianOld(JT, 1, 0);
    H.call(List.of(x0_, ca.vec(A_))).get(0);
  }

  @Test
  public void test_collocationPoints() {
    message("collocation points");
    assertThrows(Exception.class, () -> ca.collocation_points(0, "radau"));
    assertThrows(Exception.class, () -> ca.collocation_points(10, "radau"));
    assertThrows(Exception.class, () -> ca.collocation_points(0, "legendre"));
    assertThrows(Exception.class, () -> ca.collocation_points(10, "legendre"));
    assertThrows(Exception.class, () -> ca.collocation_points(1, "foo"));
    for (int k = 1; k < 10; k++) {
      double[] r = ca.collocation_points(k, "radau");
      assertEquals(k + 1, r.length + 1);
      chk(new DM(r[r.length - 1]), dm(1), 10);
    }
    for (int k = 1; k < 10; k++) assertEquals(k + 1, ca.collocation_points(k, "legendre").length + 1);
  }

  @Test
  public void test_thread_safety() {
    MX x = MX.sym("x");
    DM lin = dm(linspace(0, 10, 40));
    for (Intg I : integrators) {
      for (int threads : new int[] {4, 2}) {
        Function intg = ca.integrator("Integrator", I.name, mxd("x", x, "ode", x.unaryMinus()), I.options);
        Function intgPar = intg.map(40, "thread", threads);
        DM xf = intgPar.call(dmd("x0", lin)).get("xf");
        chk(ca.norm_inf(xf.T().minus(lin.times(Math.exp(-1)))), 0, 5);
      }
    }
  }

  @Test
  public void test_simplify_zdim() {
    MX x = MX.sym("x");
    Function intg = ca.integrator("intg", "rk", mxd("x", x, "ode", x.power(2)), opts("simplify", true));
    assertEquals(0, intg.nnz_out("zf"));
  }

  @Test
  public void test_step_options_cvodes() {
    Assumptions.assumeTrue(hasIntegrator("cvodes"));
    SX x = SX.sym("x");
    Function I = ca.integrator("I", "cvodes", sxd("x", x, "ode", ca.sin(x.times(10).power(2))), 0.0, 0.5,
        opts("step0", 1e-4, "min_step_size", 1e-4, "max_step_size", 1.1e-4));
    I.call(dmd("x0", 0));
    double n = num(I.stats().get("nsteps"));
    assertTrue((int) (0.5 / 1e-4) >= n && n >= (int) (0.5 / 1.1e-4));
  }

  @Test
  public void test_step_options_idas() {
    Assumptions.assumeTrue(hasIntegrator("idas"));
    SX x = SX.sym("x");
    Function I = ca.integrator("I", "cvodes", sxd("x", x, "ode", ca.sin(x.times(10).power(2))), 0.0, 0.5,
        opts("step0", 1e-4, "max_step_size", 1.1e-4));
    I.call(dmd("x0", 0));
    double n = num(I.stats().get("nsteps"));
    message("" + n);
    assertTrue(n >= (int) (0.5 / 1.1e-4));
  }

  @Test
  public void test_adjoint_without_parameters() {
    message("Adjoint sensitivities without parameter or control quadratures");
    MX x = MX.sym("x", 3), x0 = MX.sym("x0", 3);
    DM rates = dm(1, 2, 3);
    for (String backend : List.of("cvodes", "idas")) {
      if (!ca.has_integrator(backend)) continue;
      for (double[] grid : new double[][] {{1}, {0.5, 1}}) {
        Function I = ca.integrator("I", backend, mxd("x", x, "ode", ca.times(rates.unaryMinus(), x)), 0, grid,
            opts("enable_forward", false, "abstol", 1e-11, "reltol", 1e-11));
        MX xf = I.call(mxd("x0", x0)).get("xf");
        Function J = new Function("J", List.of(x0), List.of(ca.jacobian(xf.get(new Slice(0, 2), -1), x0)));
        DM reference = ca.diag(ca.exp(rates.unaryMinus())).get(new Slice(0, 2), new Slice());
        for (double[] value : new double[][] {{1, 2, 3}, {3, 2, 1}})
          chk(reference, J.call(List.of(dm(value))).get(0), 7);
      }
    }
  }

  @Test
  public void test_constraints_idas() {
    Assumptions.assumeTrue(hasIntegrator("idas"));
    SX x = SX.sym("x"), z = SX.sym("z"), p = SX.sym("p");
    Map<String, SX> dae = sxd("x", x, "z", z, "p", p, "ode", z.plus(p), "alg", z.times(ca.cos(z)).minus(x));
    Function I = ca.integrator("I", "idas", dae, opts("constraints", new long[] {1, 1}));
    I.call(dmd("x0", 0, "p", 0.15));
    Function I2 = ca.integrator("I", "idas", dae, opts("constraints", new long[] {-1, -1}));
    assertInException("IDA_CONSTR_FAIL", () -> I2.call(dmd("x0", 0, "p", 0.15)));
  }

  @Test
  public void test_constraints_idas_sensitivities() {
    Assumptions.assumeTrue(hasIntegrator("idas"));
    message("IDAS constraints apply to primal states, not sensitivities");
    MX x = MX.sym("x", 2), z = MX.sym("z"), p = MX.sym("p", 2);
    Function I = ca.integrator("I", "idas",
        mxd("x", x, "z", z, "p", p, "ode", p.unaryMinus().times(x), "alg", z.plus(ca.sum1(x))), 0, 1,
        opts("constraints", new long[] {2, 2, -2}, "enable_reverse", false, "abstol", 1e-11, "reltol", 1e-11));
    Map<String, MX> result = I.call(mxd("x0", dm(1, 2), "z0", -3, "p", p));
    MX y = vm(result.get("xf"), result.get("zf"));
    Function J = new Function("J", List.of(p), List.of(ca.jacobian(y, p)));
    double e1 = Math.exp(-1), e2 = Math.exp(-2);
    DM reference = dm(new double[][] {{-e1, 0}, {0, -2 * e2}, {e1, 2 * e2}});
    chk(reference, J.call(List.of(dm(1, 2))).get(0), 7);
    chk(reference, Function.deserialize(J.serialize()).call(List.of(dm(1, 2))).get(0), 7);
  }

  private static List<SX> permuted(List<SX> s, int[] perm) {
    List<SX> r = new ArrayList<>();
    for (int e : perm) r.add(s.get(e));
    return r;
  }

  /** numpy's legacy RandomState (MT19937): np.random.seed(s); np.random.permutation(range(n)) */
  private static final class NumpyRandom {
    private final int[] mt = new int[624];
    private int mti;

    NumpyRandom(long seed) {
      mt[0] = (int) seed;
      for (mti = 1; mti < 624; mti++) mt[mti] = 1812433253 * (mt[mti - 1] ^ (mt[mti - 1] >>> 30)) + mti;
    }

    private int next() {
      if (mti >= 624) {
        for (int k = 0; k < 624; k++) {
          int y = (mt[k] & 0x80000000) | (mt[(k + 1) % 624] & 0x7fffffff);
          mt[k] = mt[(k + 397) % 624] ^ (y >>> 1) ^ ((y & 1) != 0 ? 0x9908b0df : 0);
        }
        mti = 0;
      }
      int y = mt[mti++];
      y ^= y >>> 11;
      y ^= (y << 7) & 0x9d2c5680;
      y ^= (y << 15) & 0xefc60000;
      y ^= y >>> 18;
      return y;
    }

    private long interval(long max) {
      if (max == 0) return 0;
      long mask = max, v;
      for (int s = 1; s <= 32; s <<= 1) mask |= mask >>> s;
      while ((v = (next() & 0xffffffffL) & mask) > max) {}
      return v;
    }

    int[] permutation(int n) {
      int[] p = new int[n];
      for (int i = 0; i < n; i++) p[i] = i;
      for (int i = n - 1; i > 0; i--) {
        int j = (int) interval(i), t = p[i];
        p[i] = p[j];
        p[j] = t;
      }
      return p;
    }
  }

  private static DM index(int[] perm) {
    double[] d = new double[perm.length];
    for (int i = 0; i < d.length; i++) d[i] = perm[i];
    return dm(d);
  }

  private static List<Object[]> reduceIndexProblems() {
    double free = 0, force = -1, approx = 1, blank = 0.0;
    List<Object[]> r = new ArrayList<>();
    double L = 1, g = 9.81;
    {
      // 2D Pendulum
      SX x = SX.sym("x"), y = SX.sym("y"), dx = SX.sym("dx"), dy = SX.sym("dy"), u = SX.sym("u"), v = SX.sym("v");
      SX du = SX.sym("du"), dv = SX.sym("dv"), T = SX.sym("T");
      SX alg = vs(dx.minus(u), dy.minus(v), du.minus(T.times(x)), dv.minus(T.times(y)).plus(g),
          x.power(2).plus(y.power(2)).minus(L * L));
      Map<String, DM> strength = dmd("x_impl", dm(free, force, free, free), "dx_impl", dm(force, free, free, free), "z", dm(free));
      Map<String, DM> init = dmd("x_impl", dm(-1.0, -0.5, blank, blank), "dx_impl", dm(-0.1, blank, blank, blank), "z", dm(blank));
      r.add(new Object[] {sxd("x_impl", vs(x, y, u, v), "dx_impl", vs(dx, dy, du, dv), "z", T, "alg", alg), 3L, strength, init, 5e-6});
    }
    {
      // 3D Pendulum
      SX x = SX.sym("x"), y = SX.sym("y"), z = SX.sym("z"), vx = SX.sym("vx"), vy = SX.sym("vy"), vz = SX.sym("vz");
      SX dx = SX.sym("dx"), dy = SX.sym("dy"), dz = SX.sym("dz"), dvx = SX.sym("dvx"), dvy = SX.sym("dvy"), dvz = SX.sym("dvz");
      SX lambd = SX.sym("lambd");
      SX alg = vs(dx.minus(vx), dy.minus(vy), dz.minus(vz), dvx.minus(lambd.times(x)), dvy.minus(lambd.times(y)).plus(g),
          dvz.minus(lambd.times(z)), x.power(2).plus(y.power(2)).plus(z.power(2)).minus(L * L));
      Map<String, DM> strength = dmd("x_impl", dm(free, force, force, free, free, free),
          "dx_impl", dm(force, force, free, free, free, free), "z", dm(free));
      Map<String, DM> init = dmd("x_impl", dm(-1.0, -0.5, -0.5, blank, blank, blank),
          "dx_impl", dm(-0.1, 0.2, blank, blank, blank, blank), "z", dm(blank));
      r.add(new Object[] {sxd("x_impl", vs(x, y, z, vx, vy, vz), "dx_impl", vs(dx, dy, dz, dvx, dvy, dvz), "z", lambd,
          "alg", alg), 3L, strength, init, 1e-5});
    }
    {
      // Flash system (CppADCodeGen dae_index_reduction/model/flash.hpp)
      SX nEthanol = SX.sym("nEthanol"), nWater = SX.sym("nWater"), T = SX.sym("T"), yWater = SX.sym("yWater");
      SX yEthanol = SX.sym("yEthanol"), FV = SX.sym("FV");
      double Q = 500, F_feed = 10, p = 1, xFEthanol = 0.5, T_feed = 50;
      SX dnEthanol = SX.sym("dnEthanol"), dnWater = SX.sym("dnWater"), dT = SX.sym("dT");
      List<SX> alg = new ArrayList<>();
      double FNEthanol = xFEthanol * F_feed;
      SX n_lWater = nWater.times(1000.), n_lEthanol = nEthanol.times(1000.);
      SX m_lEthanol = n_lEthanol.times(0.04606844), m_lWater = n_lWater.times(0.0180152833);
      SX m = m_lEthanol.plus(m_lWater);
      SX wEthanol = m_lEthanol.div(m);
      SX w_Water = ca.minus(1, wEthanol);
      SX rho = ca.rdivide(1, w_Water.div(983.159471259596).plus(wEthanol.div(743.278841365274)));
      SX V = m.div(rho);
      SX cWater = n_lWater.div(V);
      SX dh = V.div(0.502654824574367);
      double p_aux = p * 101325., dp = 101325. - p_aux;
      SX v = ca.sqrt(ca.times(9.80665, dh).plus(ca.rdivide(dp, rho)).times(2.));
      SX F_VL = v.times(0.0005);
      SX FNLWater = cWater.times(F_VL);
      SX n = n_lEthanol.plus(n_lWater);
      SX xWater = n_lWater.div(n);
      SX F_NL = FNLWater.div(xWater);
      SX FNLEthanol = F_NL.minus(FNLWater);
      SX FNVWater = yWater.times(FV);
      SX FNVEthanol = FV.minus(FNVWater);
      SX dnE = ca.minus(FNEthanol, FNLEthanol).minus(FNVEthanol);
      alg.add(dnEthanol.minus(dnE.times(0.001)));
      double FNFWater = F_feed - FNEthanol;
      SX dnW = ca.minus(FNFWater, FNLWater).minus(FNVWater);
      alg.add(dnWater.minus(dnW.times(0.001)));
      double FMFEthanol = FNEthanol * 0.04606844, FMFWater = FNFWater * 0.0180152833, Tfeed = T_feed + 273.15;
      SX T_aux = T.plus(273.15);
      SX dQ_F = ca.times(FMFEthanol * 2898.42878374706 + FMFWater * 4186.92536027523, ca.minus(Tfeed, T_aux));
      SX dH_mVapEthanol = ca.times(50430., ca.power(ca.minus(1, T_aux.div(514.)), new SX(0.4989)))
          .times(ca.exp(ca.times(0.4475, T_aux).div(514.)));
      SX T_r = T_aux.div(647.);
      SX dH_mVapWater = ca.times(52053., ca.power(ca.minus(1, T_r),
          ca.plus(0.3199, ca.times(-0.212, T_r)).plus(ca.times(0.25795, T_r).times(T_r))));
      SX dH_mVap = yEthanol.times(dH_mVapEthanol).plus(yWater.times(dH_mVapWater));
      SX dQ_vap = FV.times(dH_mVap);
      double Q_1 = Q * 1000.;
      SX cp = wEthanol.times(2898.42878374706).plus(w_Water.times(4186.92536027523));
      alg.add(dT.minus(dQ_F.minus(dQ_vap).plus(Q_1).div(m.times(cp))));
      SX pVapWater = ca.times(100000., ca.power(new SX(10.), ca.minus(4.6543, ca.rdivide(1435.264, T_aux.minus(64.848)))));
      SX KWater = pVapWater.div(p_aux);
      alg.add(yWater.minus(xWater.times(KWater)));
      SX xEthanol = ca.minus(1, xWater);
      SX pVapEthanol = ca.exp(ca.plus(74.475, ca.rdivide(-7164.3, T_aux)).plus(ca.times(-7.327, ca.log(T_aux)))
          .plus(ca.times(3.134e-06, ca.power(T_aux, new SX(2.)))));
      SX KEthanol = pVapEthanol.div(p_aux);
      alg.add(yEthanol.minus(xEthanol.times(KEthanol)));
      alg.add(yWater.plus(yEthanol).minus(1));
      Map<String, DM> strength = dmd("x_impl", dm(force, force, approx), "dx_impl", dm(free, free, free),
          "z", dm(free, free, free));
      Map<String, DM> init = dmd("x_impl", dm(2.5, 6.4, 91), "dx_impl", dm(0.1, blank, blank), "z", dm(0.5, 0.5, 0.5));
      r.add(new Object[] {sxd("x_impl", ca.veccat(List.of(nEthanol, nWater, T)), "dx_impl",
          ca.veccat(List.of(dnEthanol, dnWater, dT)), "z", vs(yWater, yEthanol, FV), "alg", ca.veccat(alg)),
          2L, strength, init, 1e-8});
    }
    {
      // Exothermic Reactor from pantelides paper
      SX C = SX.sym("C"), dC = SX.sym("dC"), T = SX.sym("T"), dT = SX.sym("dT"), R = SX.sym("R");
      double C0 = 1, T0 = 1, K1 = 1, K2 = 1, K3 = 1, K4 = 1, Tc = 1;
      SX alg = vs(ca.times(K1, ca.minus(C0, C)).minus(R).minus(dC),
          ca.times(K1, ca.minus(T0, T)).plus(R.times(K2)).minus(T.minus(Tc).times(K3)).minus(dT),
          R.minus(ca.times(K3, ca.exp(ca.rdivide(-K4, T))).times(C)));
      Map<String, DM> strength = dmd("x_impl", dm(free, force), "dx_impl", dm(free, free), "z", dm(free));
      Map<String, DM> init = dmd("x_impl", dm(blank, 2), "dx_impl", dm(blank, blank), "z", dm(blank));
      r.add(new Object[] {sxd("x_impl", ca.veccat(List.of(C, T)), "dx_impl", ca.veccat(List.of(dC, dT)), "z", vs(R),
          "alg", alg), 1L, strength, init, 1e-8});
    }
    {
      // S. Mattsson and G. Soederlind, Index reduction in differential-algebraic equations using dummy derivatives
      SX[] xs = new SX[7], dxs = new SX[7];
      for (int i = 0; i < 7; i++) {
        xs[i] = SX.sym("x" + (i + 1));
        dxs[i] = SX.sym("dx" + (i + 1));
      }
      double u1 = 1, u2 = 1, u3 = 1, u4 = 1;
      SX eqs = vs(xs[4].minus(dxs[0]), xs[5].minus(dxs[1]), xs[6].minus(dxs[2]), xs[0].plus(xs[1]).plus(u1),
          xs[0].plus(xs[1]).plus(xs[2]).plus(u2), dxs[2].plus(xs[0]).plus(xs[3]).plus(u3),
          dxs[4].times(2).plus(dxs[5]).plus(dxs[6]).plus(dxs[3]).plus(u4));
      r.add(new Object[] {sxd("x_impl", ca.veccat(Arrays.asList(xs)), "dx_impl", ca.veccat(Arrays.asList(dxs)),
          "z", new SX(0, 1), "alg", eqs), 3L, null, null, null});
    }
    {
      // Uncontrollable DAE Pantelides
      SX x = SX.sym("x"), dx = SX.sym("dx"), u1 = SX.sym("u1"), u2 = SX.sym("u2"), t = SX.sym("t");
      SX y1 = t, y2 = ca.sin(t);
      SX alg = vs(x.plus(u1).plus(u2), x.plus(dx).plus(y1), x.plus(y2));
      r.add(new Object[] {sxd("x_impl", vs(x), "dx_impl", vs(dx), "z", vs(u1, u2), "alg", alg, "t", t),
          "Structural detection failed", null, null, null});
    }
    return r;
  }

  @Test
  @SuppressWarnings("unchecked")
  public void test_reduce_index() {
    Assumptions.assumeTrue(hasIntegrator("idas"));
    Assumptions.assumeTrue(hasNlpsol("ipopt"));
    message("loop");
    for (Object[] pr : reduceIndexProblems()) {
      Map<String, SX> daeOrig = (Map<String, SX>) pr[0];
      Map<String, DM> strengthOrig = (Map<String, DM>) pr[2], initOrig = (Map<String, DM>) pr[3];
      message("dae " + daeOrig);
      for (int seed = 0; seed < 5; seed++) {
        Map<String, SX> dae = new LinkedHashMap<>(daeOrig);
        Map<String, DM> strength = strengthOrig == null ? null : new LinkedHashMap<>(strengthOrig);
        Map<String, DM> init = initOrig == null ? null : new LinkedHashMap<>(initOrig);
        if (seed > 0) {
          NumpyRandom rng = new NumpyRandom(seed);
          for (String k : List.of("x_impl", "z", "alg")) {
            List<SX> s = ca.vertsplit(dae.get(k));
            if (s.isEmpty()) continue;
            int[] perm = rng.permutation(s.size());
            dae.put(k, ca.vertcat(permuted(s, perm)));
            if (init != null && init.containsKey(k)) init.put(k, init.get(k).get(index(perm)));
            if (strength != null && strength.containsKey(k)) strength.put(k, strength.get(k).get(index(perm)));
            if (k.equals("x_impl")) {
              s = ca.vertsplit(dae.get("dx_impl"));
              dae.put("dx_impl", ca.vertcat(permuted(s, perm)));
              if (init != null) init.put("dx_impl", init.get("dx_impl").get(index(perm)));
              if (strength != null) strength.put("dx_impl", strength.get("dx_impl").get(index(perm)));
            }
          }
        }
        if (pr[1] instanceof String) {
          Map<String, SX> d = dae;
          assertInException((String) pr[1], () -> ca.dae_reduce_index(d,
              opts("baumgarte_pole", -1, "algorithm", "pantelides", "max_iter", 20)));
          continue;
        }
        List<Object> red = ca.dae_reduce_index(dae, opts("baumgarte_pole", -1, "algorithm", "pantelides"));
        Map<String, SX> daeReduced = (Map<String, SX>) red.get(0);
        Map<String, Object> meta = (Map<String, Object>) red.get(1);
        assertEquals(pr[1], meta.get("index"));
        if (strength != null) {
          List<Object> se = ca.dae_map_semi_expl(dae, daeReduced);
          Map<String, SX> daeSe = (Map<String, SX>) se.get(0);
          Function phi = (Function) se.get(2);
          Function intg = ca.integrator("intg", "idas", daeSe, 0, linspace(0, 5, 1000));
          Function initGen = ca.dae_init_gen(dae, daeReduced, "ipopt", strength, opts("error_on_fail", true));
          Map<String, DM> xz0 = initGen.call(init);
          Map<String, DM> sol = intg.call(xz0);
          DM error = phi.call(dmd("x", sol.get("xf"), "z", sol.get("zf"))).get("I");
          double[][] e = error.full();
          for (double[] row : e) {
            double s = 0;
            for (double v : row) s += v * v;
            double rms = Math.sqrt(s / row.length);
            message("rms error " + rms);
            assertTrue(rms <= (Double) pr[4]);
          }
        }
      }
    }
  }

  @Test
  public void test_IDACalcIC_u() {
    MX x = MX.sym("x"), u = MX.sym("u"), z = MX.sym("z");
    Map<String, MX> ivp = mxd("x", x, "z", z, "u", u, "ode", z, "alg", ca.sqrt(u).minus(z));
    Function intg = ca.integrator("intg", "idas", ivp, 0, 1);
    Map<String, DM> res = intg.call(dmd("x0", 0, "u", 1));
    almost(res.get("xf").to_double(), 1, 5, "test_IDACalcIC_u");
  }

  @Test
  public void test_expand() {
    Assumptions.assumeTrue(hasIntegrator("cvodes"));
    MX t = MX.sym("t"), q = MX.sym("q"), p = MX.sym("p");
    Map<String, MX> dae = mxd("t", t, "x", q, "p", p, "ode", q.div(p).times(t.power(2)));
    Function integrator = ca.integrator("integrator", "cvodes", dae, 0, 2.3, opts("expand", true));
    assertTrue(integrator.get_function("jacF").is_a("SXFunction"));
    t = MX.sym("t");
    q = MX.sym("q");
    p = MX.sym("p");
    dae = mxd("t", t, "x", q, "p", p, "ode", q.div(p).times(t.power(2)));
    integrator = ca.integrator("integrator", "cvodes", dae, 0, 2.3, opts("expand", false));
    assertTrue(integrator.get_function("jacF").is_a("MXFunction"));
  }

  @Test
  public void test_postpone_expand() {
    Assumptions.assumeTrue(hasIntegrator("cvodes"));
    MX x = MX.sym("x");
    Function J = new Function("jac_f", List.of(x, new MX(1, 1)), List.of(x.unaryMinus()), List.of("x", "out_r"),
        List.of("jac_r_x"), opts("always_inline", true));
    Function f = new Function("f", List.of(x), List.of(x.power(2)), List.of("x"), List.of("r"),
        opts("custom_jacobian", J, "jac_penalty", 0));
    Map<String, MX> dae = mxd("x", x, "ode", (MX) f.invoke(x.minus(2)));
    Object[][] cases = {{opts("expand", false), "MXFunction", 1.0}, {opts("expand", true), "SXFunction", -2.0},
        {opts("expand", false, "postpone_expand", true), "MXFunction", 1.0},
        {opts("expand", true, "postpone_expand", true), "SXFunction", 1.0}};
    for (Object[] c : cases) {
      @SuppressWarnings("unchecked")
      Function integrator = ca.integrator("integrator", "cvodes", dae, 0, 0.1, (Map<String, Object>) c[0]);
      integrator.call(dmd("x0", 1));
      assertTrue(integrator.get_function("jacF").is_a((String) c[1]));
      chk(integrator.get_function("jacF").call(dmd("x", 1)).get("jac_ode_x"), (Double) c[2], 10);
    }
  }

  @Test
  public void test_idas_adjoint_zf_p() {
    Assumptions.assumeTrue(hasIntegrator("idas"));
    MX x = MX.sym("x"), z = MX.sym("z"), p = MX.sym("p"), u = MX.sym("u");
    for (MX alg : List.of(z.minus(p.times(x)).minus(u.power(2)), z.minus(p.times(x)))) {
      for (MX ode : List.of(x.unaryMinus(), x.unaryMinus().plus(u))) {
        double[] grid = {1};
        Function F = ca.integrator("F", "idas", mxd("x", x, "z", z, "p", p, "u", u, "ode", ode, "alg", alg), 0, grid,
            opts("abstol", 1e-12, "reltol", 1e-12));
        MX P = MX.sym("P"), U = MX.sym("U", 1, grid.length);
        Map<String, MX> out = F.call(mxd("x0", 2, "z0", 0, "p", P, "u", U));
        MX e = vm(ca.vec(out.get("zf")), ca.vec(out.get("xf")));
        MX w = vm(P, ca.vec(U));
        Function Jfwd = new Function("J_fwd", List.of(w),
            List.of(ca.jacobian(e, w, opts("helper_options", opts("ad_weight", 0)))));
        Function Jrev = new Function("J_rev", List.of(w),
            List.of(ca.jacobian(e, w, opts("helper_options", opts("ad_weight", 1)))));
        DM v = dm(Arrays.copyOf(new double[] {0.5, 0.3, 0.7}, 1 + grid.length));
        chk(Jfwd.call(List.of(v)).get(0), Jrev.call(List.of(v)).get(0), 6);
      }
    }
  }

  @Test
  public void test_jit_fwd_factory() {
    Assumptions.assumeTrue(hasIntegrator("cvodes"));
    Assumptions.assumeTrue(RUN_SLOW);
    MX x = MX.sym("x");
    Function F = ca.integrator("F", "cvodes", mxd("x", x, "ode", x.unaryMinus()), 0, new double[] {0.5, 1},
        opts("abstol", 1e-10, "reltol", 1e-10, "jit", true, "compiler", "shell"));
    assertTrue(F.get_function("jacF").is_a("External"));
    Function Ifwd = F.factory("I_fwd", List.of("x0", "fwd:x0"), List.of("xf", "fwd:xf"));
    List<DM> r = Ifwd.call(List.of(dm(2.0), dm(1.0)));
    DM ref = ca.exp(dm(new double[][] {{0.5, 1}}).unaryMinus());
    chk(r.get(0), ref.times(2), 7);
    chk(r.get(1), ref, 7);
  }

  private void crossCheck(Function I, Map<String, DM> base, List<Map<String, DM>> seeds, String msg) {
    int N = seeds.size();
    Function I1 = I.forward(1), IN = I.forward(N);
    List<String> fwdNames = new ArrayList<>();
    for (String n : IN.name_in()) if (n.startsWith("fwd_")) fwdNames.add(n);
    Map<String, DM> batchedCall = new LinkedHashMap<>(base);
    for (String n : fwdNames) {
      List<DM> cols = new ArrayList<>();
      for (Map<String, DM> s : seeds) cols.add(s.getOrDefault(n, DM.zeros(I1.sparsity_in(n))));
      batchedCall.put(n, ca.horzcat(cols));
    }
    Map<String, DM> batched = IN.call(batchedCall);
    Map<String, List<DM>> perSeed = new LinkedHashMap<>();
    for (Map<String, DM> s : seeds) {
      Map<String, DM> single = new LinkedHashMap<>(base);
      for (String n : fwdNames) single.put(n, s.getOrDefault(n, DM.zeros(I1.sparsity_in(n))));
      for (Map.Entry<String, DM> e : I1.call(single).entrySet())
        perSeed.computeIfAbsent(e.getKey(), k -> new ArrayList<>()).add(e.getValue());
    }
    for (Map.Entry<String, DM> e : batched.entrySet())
      chk(e.getValue(), ca.horzcat(perSeed.get(e.getKey())), 10, msg + " :: " + e.getKey());
  }

  @Test
  public void test_issue3371() {
    MX X = MX.sym("x", 2), U = MX.sym("u"), Z = MX.sym("z");
    Map<String, MX> dae = mxd("x", X, "p", U, "z", Z, "ode", vm(Z.times(X.get(0)).minus(X.get(1)).plus(U), X.get(0)),
        "alg", Z.minus(ca.minus(1, X.get(1).power(2))));
    Function I = ca.integrator("I", "collocation", dae);
    Map<String, DM> base = dmd("x0", dm(0.5, 0.0), "z0", dm(1.0), "p", dm(0.1), "u", new DM(0, 1));
    Map<String, DM> r = I.call(dmd("x0", base.get("x0"), "z0", base.get("z0"), "p", base.get("p")));
    base.put("out_xf", r.get("xf"));
    base.put("out_zf", r.get("zf"));
    base.put("out_qf", r.get("qf"));
    crossCheck(I, base, List.of(dmd("fwd_x0", dm(1.0, 0.0)), dmd("fwd_x0", dm(0.0, 1.0)), dmd("fwd_p", dm(1.0))),
        "fwd-of-collocation");

    Function I2 = I.forward(2);
    Map<String, DM> base0 = dmd("x0", DM.zeros(2), "z0", dm(1.0), "p", dm(0.0), "u", new DM(0, 1),
        "out_xf", DM.zeros(2, 1), "out_zf", dm(1.0), "out_qf", new DM(0, 1),
        "fwd_x0", dm(new double[][] {{0, 1}, {0, 0}}), "fwd_z0", DM.zeros(I2.sparsity_in("fwd_z0")),
        "fwd_p", dm(1, 0).T(), "fwd_u", DM.zeros(I2.sparsity_in("fwd_u")));
    chk(I2.call(base0).get("fwd_zf"), DM.zeros(1, 2), 10, "issue3371 fwd_zf");

    Function Iadj = ca.integrator("I", "collocation", dae, opts("nadj", 1));
    Map<String, DM> ref = Iadj.call(dmd("x0", dm(0.5, 0.0), "z0", dm(1.0), "p", dm(0.1), "adj_xf", dm(0.0, 0.0),
        "adj_zf", dm(1.0), "adj_qf", new DM(0, 1)));
    base = dmd("x0", dm(0.5, 0.0), "z0", dm(1.0), "p", dm(0.1), "u", new DM(0, 1), "adj_xf", dm(0.0, 0.0),
        "adj_zf", dm(1.0), "adj_qf", new DM(0, 1), "out_xf", ref.get("xf"), "out_zf", ref.get("zf"),
        "out_qf", ref.get("qf"), "out_adj_x0", ref.get("adj_x0"), "out_adj_z0", ref.get("adj_z0"),
        "out_adj_p", ref.get("adj_p"), "out_adj_u", ref.get("adj_u"));
    crossCheck(Iadj, base, List.of(dmd("fwd_adj_zf", dm(1.0)), dmd("fwd_x0", dm(1.0, 0.0)), dmd("fwd_p", dm(1.0))),
        "fwd-of-reverse-of-collocation");
  }

  /** Forward and adjoint Jacobians of out(v) at v0 agree */
  private void checkJacobian(MX out, MX v, DM v0, int digits) {
    DM Jfwd = new Function("J", List.of(v), List.of(ca.jtimes(out, v, MX.eye(v.numel())))).call(List.of(v0)).get(0);
    DM Jadj = new Function("J", List.of(v), List.of(ca.jtimes(out, v, MX.eye(out.numel()), true).T()))
        .call(List.of(v0)).get(0);
    chk(Jfwd, Jadj, digits, "fwd vs adj");
  }

  /** Hessian of obj(v) at v0: forward and adjoint mode over the forward and adjoint gradients */
  private DM checkHessian(MX obj, MX v, DM v0, int digits) {
    long n = v.numel();
    MX gFwd = ca.jtimes(obj, v, MX.eye(n)).T();
    MX gAdj = ca.jtimes(obj, v, MX.ones(1, 1), true);
    List<DM> H = new ArrayList<>();
    for (MX h : List.of(ca.jtimes(gFwd, v, MX.eye(n)), ca.jtimes(gAdj, v, MX.eye(n)),
        ca.jtimes(gFwd, v, MX.eye(n), true).T(), ca.jtimes(gAdj, v, MX.eye(n), true).T()))
      H.add(new Function("H", List.of(v), List.of(h)).call(List.of(v0)).get(0));
    for (DM h : H.subList(1, H.size())) chk(H.get(0), h, digits, "Hessian");
    return H.get(0);
  }

  private static MX sl(MX v, int a, int b) {
    return v.get(new Slice(a, b));
  }

  private static DM randFor(Random rng, Sparsity sp) {
    double[] v = new double[(int) sp.nnz()];
    for (int k = 0; k < v.length; k++) v[k] = rng.nextDouble();
    return new DM(sp, new DM(v));
  }

  @SafeVarargs
  private static List<DM> join(List<DM>... a) {
    List<DM> r = new ArrayList<>();
    for (List<DM> e : a) r.addAll(e);
    return r;
  }

  @Test
  public void test_reverse_grid_directions() {
    message("Several adjoint directions with several output times");
    SX x = SX.sym("x", 2), p = SX.sym("p");
    Map<String, SX> dae = sxd("x", x, "p", p, "ode", vs(x.get(1), p.unaryMinus().times(x.get(0).power(3))),
        "quad", x.get(0).power(2));
    Object[][] cases = {{"cvodes", opts("reltol", 1e-10, "abstol", 1e-10)}, {"idas", opts("reltol", 1e-10, "abstol", 1e-10)},
        {"rosenbrock", opts("step_size", 0.125)}, {"erk", opts("method", "rk45", "reltol", 1e-8, "abstol", 1e-8)}};
    for (Object[] c : cases) {
      String plugin = (String) c[0];
      if (!ca.has_integrator(plugin)) continue;
      @SuppressWarnings("unchecked")
      Function F = ca.integrator("F", plugin, dae, 0, new double[] {0.5, 1.0}, (Map<String, Object>) c[1]);
      List<DM> inputs = new ArrayList<>(List.of(dm(0.3, -0.2), new DM(), new DM(0.7), new DM(1, 0)));
      for (long i = 4; i < F.n_in(); i++) inputs.add(new DM(F.sparsity_in(i)));
      List<DM> outputs = F.call(inputs);
      List<List<DM>> seeds = new ArrayList<>();
      for (int d = 0; d < 2; d++) {
        List<DM> s = new ArrayList<>();
        for (int i = 0; i < F.n_out(); i++) s.add(randFor(new Random(10 * d + i), F.sparsity_out(i)));
        seeds.add(s);
      }
      List<DM> both = new ArrayList<>();
      for (int i = 0; i < F.n_out(); i++) both.add(ca.horzcat(List.of(seeds.get(0).get(i), seeds.get(1).get(i))));
      List<DM> adj2 = F.reverse(2).call(join(inputs, outputs, both));
      for (int d = 0; d < 2; d++) {
        List<DM> adj1 = F.reverse(1).call(join(inputs, outputs, seeds.get(d)));
        for (int i = 0; i < F.n_in(); i++) {
          if (F.nnz_in(i) == 0) continue;
          chk(adj1.get(i), ca.horzsplit_n(adj2.get(i), 2).get(d), 8, plugin);
        }
      }
    }
  }

  @Test
  public void test_events_sparsity() {
    message("zero-crossing events: dependencies through the transition");
    double g = 9.81;
    SX x = SX.sym("x", 2), p = SX.sym("p");
    Map<String, SX> dae = sxd("x", x, "p", p, "ode", vs(x.get(1), -g), "zero", x.get(0));
    double[] grid = linspace(0.1, 4, 40);
    Object[][] cases = {{"rk", opts("number_of_finite_elements", 100)}, {"cvodes", opts("abstol", 1e-12, "reltol", 1e-12)}};
    for (Object[] c : cases) {
      String plugin = (String) c[0];
      if (!ca.has_integrator(plugin)) continue;
      @SuppressWarnings("unchecked")
      Function F = ca.integrator("F", plugin, dae, 0, grid, with((Map<String, Object>) c[1], "transition",
          bounceTransition(), "event_tol", 1e-12, "max_event_iter", 20));
      MX v = MX.sym("v", 3);
      MX xf = F.call(mxd("x0", sl(v, 0, 2), "p", v.get(2))).get("xf").get(new Slice(), -1);
      Function G = new Function("G", List.of(v), List.of(xf, ca.jacobian(xf, v)));
      assertTrue(G.sparsity_out(1).is_dense());
      DM v0 = dm(2, 0.3, 0.8);
      chk(G.call(List.of(v0)).get(1), fdJacobian(G, v0), 5, plugin);
    }
  }

  @Test
  public void test_events_idas() {
    Assumptions.assumeTrue(hasIntegrator("idas"));
    message("idas: zero-crossing events, also of a DAE");
    double g = 9.81;
    SX x = SX.sym("x", 2), z = SX.sym("z"), p = SX.sym("p");
    double[] grid = linspace(0.1, 1.5, 15);
    double tb = Math.sqrt(4 / g), ub = 0.8 * g * tb;
    double[] h = new double[grid.length], v = new double[grid.length];
    for (int i = 0; i < grid.length; i++) {
      double tt = grid[i];
      h[i] = tt < tb ? 2 - g * tt * tt / 2 : ub * (tt - tb) - g * Math.pow(tt - tb, 2) / 2;
      v[i] = tt < tb ? -g * tt : ub - g * (tt - tb);
    }
    Map<String, SX> ode = sxd("x", x, "p", p, "ode", vs(x.get(1), -g), "zero", x.get(0));
    Map<String, SX> dae = sxd("x", x, "z", z, "p", p, "ode", vs(z, -g), "alg", z.minus(x.get(1)), "zero", x.get(0));
    List<Object[]> problems = List.of(new Object[] {ode, bounceTransition()}, new Object[] {dae, bounceTransition(1)});
    for (Object[] pt : problems) {
      @SuppressWarnings("unchecked")
      Function F = ca.integrator("F", "idas", (Map<String, SX>) pt[0], 0, grid, opts("transition", pt[1],
          "abstol", 1e-12, "reltol", 1e-12, "event_tol", 1e-12, "max_event_iter", 20));
      Map<String, DM> r = F.call(dmd("x0", dm(2, 0), "z0", 0, "p", 0.8));
      chk(r.get("xf"), new DM(new double[][] {h, v}), 7, "xf");
    }
  }

  @Test
  public void test_reverse_over_adjoint_directions() {
    message("Derivatives of adjoint sensitivities with several directions (#3293)");
    SX x = SX.sym("x", 3), p = SX.sym("p", 2);
    Map<String, SX> dae = sxd("x", x, "p", p, "ode", vs(p.get(0).unaryMinus().times(x.get(0)),
        p.get(1).unaryMinus().times(x.get(1)).times(x.get(0)), p.get(1).unaryMinus().times(x.get(2).power(2))),
        "quad", x.get(0).times(x.get(2)));
    Random rng = new Random(1);
    Object[][] cases = {{"cvodes", opts("reltol", 1e-12, "abstol", 1e-12), 7},
        {"idas", opts("reltol", 1e-12, "abstol", 1e-12), 7}, {"rosenbrock", opts("step_size", 1. / 6), 12},
        {"rosenbrock", opts("reltol", 1e-6, "abstol", 1e-6), 11}, {"erk", opts("step_size", 0.2), 12},
        {"erk", opts("method", "rk45", "reltol", 1e-6, "abstol", 1e-6), 11}};
    for (Object[] c : cases) {
      String plugin = (String) c[0];
      int digits = (Integer) c[2];
      if (!ca.has_integrator(plugin)) continue;
      for (double[] grid : new double[][] {{1.0}, {0.5, 1.0}}) {
        @SuppressWarnings("unchecked")
        Function F = ca.integrator("F", plugin, dae, 0, grid, (Map<String, Object>) c[1]);
        Function B = F.reverse(2);
        List<DM> inputs = new ArrayList<>(), v = new ArrayList<>();
        for (int i = 0; i < B.n_in(); i++) inputs.add(randFor(rng, B.sparsity_in(i)));
        List<DM> outputs = B.call(inputs);
        for (int i = 0; i < B.n_in(); i++) v.add(randFor(rng, B.sparsity_in(i)));
        List<List<DM>> w = new ArrayList<>();
        for (int d = 0; d < 2; d++) {
          List<DM> s = new ArrayList<>();
          for (int i = 0; i < B.n_out(); i++) s.add(randFor(rng, B.sparsity_out(i)));
          w.add(s);
        }
        List<DM> fwd = B.forward(1).call(join(inputs, outputs, v));
        List<DM> both = new ArrayList<>();
        for (int i = 0; i < B.n_out(); i++) both.add(ca.horzcat(List.of(w.get(0).get(i), w.get(1).get(i))));
        List<DM> adj2 = B.reverse(2).call(join(inputs, outputs, both));
        for (int d = 0; d < 2; d++) {
          List<DM> adj = B.reverse(1).call(join(inputs, outputs, w.get(d)));
          double lhs = 0, rhs = 0;
          for (int i = 0; i < B.n_out(); i++) lhs += ca.dot(w.get(d).get(i), fwd.get(i)).to_double();
          for (int i = 0; i < B.n_in(); i++) rhs += ca.dot(adj.get(i), v.get(i)).to_double();
          chk(new DM(lhs), new DM(rhs), digits, plugin + " duality");
          for (int i = 0; i < B.n_in(); i++)
            if (B.nnz_in(i) > 0) chk(adj.get(i), ca.horzsplit_n(adj2.get(i), 2).get(d), digits, plugin);
        }
        long total = 0;
        long[] offsets = new long[(int) B.n_in() + 1];
        for (int i = 0; i < B.n_in(); i++) offsets[i + 1] = offsets[i] + B.nnz_in(i);
        total = offsets[offsets.length - 1];
        MX vin = MX.sym("v", total);
        List<MX> args = ca.vertsplit(vin, offsets);
        List<MX> resIn = new ArrayList<>();
        for (int i = 0; i < B.n_in(); i++) resIn.add(new MX(B.sparsity_in(i), args.get(i)));
        List<MX> vouts = new ArrayList<>();
        for (MX e : B.call(resIn)) vouts.add(ca.vec(e.get_nz(false, new Slice())));
        MX vout = ca.vertcat(vouts);
        Function J = new Function("J", List.of(vin), List.of(ca.jacobian(vout, vin), ca.jtimes(vout, vin, MX.eye(vin.numel()))));
        List<Double> nz = new ArrayList<>();
        for (DM e : inputs) for (double d : e.nonzeros()) nz.add(d);
        List<DM> Jo = J.call(List.of(dm(nz.stream().mapToDouble(Double::doubleValue).toArray())));
        chk(Jo.get(1), Jo.get(0), digits, plugin + " sparsity");
      }
    }
  }

  @Test
  public void test_runtime_integrator_quadrature_error_control() {
    SX t = SX.sym("t"), x = SX.sym("x");
    Map<String, SX> dae = sxd("t", t, "x", x, "ode", ca.times(0, x), "quad", ca.exp(t.times(10)));
    double exact = Math.expm1(10) / 10;
    Map<String, Object> o = opts("abstol", 1e-9, "reltol", 1e-9, "step0", 1.0);
    for (String plugin : List.of("erk", "rosenbrock")) {
      if (!ca.has_integrator(plugin)) continue;
      Function def = ca.integrator("default", plugin, dae, 0, 1, o);
      Function enabled = ca.integrator("enabled", plugin, dae, 0, 1, with(o, "quad_err_con", true));
      Function disabled = ca.integrator("disabled", plugin, dae, 0, 1, with(o, "quad_err_con", false));
      double qd = def.call(dmd("x0", 0)).get("qf").to_double();
      double qe = enabled.call(dmd("x0", 0)).get("qf").to_double();
      double qn = disabled.call(dmd("x0", 0)).get("qf").to_double();
      assertEquals(qd, qe);
      assertTrue(Math.abs(qd - exact) < 1e-4);
      assertTrue(Math.abs(qn - exact) > 1e-2);
      assertTrue(num(def.stats().get("n_steps")) > num(disabled.stats().get("n_steps")));
    }
  }

  @Test
  public void test_rosenbrock_order() {
    message("rosenbrock: order 5 on fixed meshes, ODE and index-1 DAE");
    SX t = SX.sym("t"), x = SX.sym("x", 2), z = SX.sym("z");
    Map<String, SX> dae = sxd("t", t, "x", x, "z", z, "ode", vs(z.unaryMinus().plus(x.get(1).times(0.5)),
        x.get(0).unaryMinus().times(z).plus(ca.sin(t))), "alg",
        z.power(3).plus(z).minus(x.get(0).power(2).times(ca.plus(1, ca.cos(t).times(0.5)))), "quad", z.times(x.get(1)));
    Map<String, DM> ref = ca.integrator("ref", "idas", dae, 0, 2, opts("reltol", 1e-14, "abstol", 1e-14))
        .call(dmd("x0", dm(1, 0.5), "z0", 0.6));
    List<Double> err = new ArrayList<>();
    for (int N : new int[] {8, 16, 32}) {
      Function F = ca.integrator("F", "rosenbrock", dae, 0, 2, opts("step_size", 2. / N));
      Map<String, DM> r = F.call(dmd("x0", dm(1, 0.5), "z0", 0.6));
      double e = 0;
      for (String k : List.of("xf", "zf", "qf")) e = Math.max(e, ca.norm_inf(r.get(k).minus(ref.get(k))).to_double());
      err.add(e);
    }
    for (int i = 0; i + 1 < err.size(); i++) assertTrue(Math.log(err.get(i) / err.get(i + 1)) / Math.log(2) > 4.6);
  }

  @Test
  public void test_rosenbrock_ind() {
    message("rosenbrock: derivatives are exact derivatives of the scheme");
    SX t = SX.sym("t"), x = SX.sym("x", 2), z = SX.sym("z"), p = SX.sym("p", 2), u = SX.sym("u");
    Map<String, SX> dae = sxd("t", t, "x", x, "z", z, "p", p, "u", u,
        "ode", vs(z.unaryMinus().plus(p.get(0).times(x.get(1))).plus(u), x.get(0).unaryMinus().times(z).plus(p.get(1).times(ca.sin(t)))),
        "alg", z.power(3).plus(z).minus(x.get(0).power(2).times(ca.plus(1, ca.cos(t).times(0.5)))),
        "quad", z.times(x.get(1)).plus(p.get(1).times(x.get(0).power(2))));
    for (Map<String, Object> o : List.of(opts("step_size", 0.1), opts("reltol", 1e-6, "abstol", 1e-6)))
      for (double[] grid : new double[][] {{0.5, 1, 2}, {0, 1, 2}}) {
        Function F = ca.integrator("F", "rosenbrock", dae, 0, grid, o);
        MX v = MX.sym("v", 7);
        Map<String, MX> r = F.call(mxd("x0", sl(v, 0, 2), "z0", 0.5, "p", sl(v, 2, 4), "u", sl(v, 4, 7).T()));
        MX out = vm(ca.vec(r.get("xf")), ca.vec(r.get("zf")), ca.vec(r.get("qf")));
        MX obj = ca.sumsqr(out);
        DM v0 = dm(1, 0.5, 0.5, 0.3, 0.1, -0.2, 0.3);
        checkJacobian(out, v, v0, 12);
        DM H = checkHessian(obj, v, v0, 12);
        chk(H, H.T(), 12, "Hessian symmetry");
      }
  }

  @Test
  public void test_rosenbrock_dense_output() {
    message("rosenbrock: outputs inside steps");
    SX t = SX.sym("t"), x = SX.sym("x", 2), z = SX.sym("z"), p = SX.sym("p", 2), u = SX.sym("u");
    Map<String, SX> dae = sxd("t", t, "x", x, "z", z, "p", p, "u", u,
        "ode", vs(z.unaryMinus().plus(p.get(0).times(x.get(1))).plus(u), x.get(0).unaryMinus().times(z).plus(p.get(1).times(ca.sin(t)))),
        "alg", z.power(3).plus(z).minus(x.get(0).power(2).times(ca.plus(1, ca.cos(t).times(0.5)))).minus(u),
        "quad", z.times(x.get(1)).plus(p.get(1).times(x.get(0).power(2))));
    Map<String, SX> ode = sxd("t", t, "x", x, "p", p, "u", u, "ode", vs(x.get(1).plus(p.get(0).times(u)),
        p.get(1).unaryMinus().times(x.get(0)).plus(ca.sin(t))), "quad", x.get(0).times(x.get(1).power(2)));
    List<Object[]> systems = List.of(new Object[] {ode, "cvodes"}, new Object[] {dae, "idas"});
    for (Object[] s : systems) {
      @SuppressWarnings("unchecked")
      Map<String, SX> sys = (Map<String, SX>) s[0];
      String refPlugin = (String) s[1];
      if (!ca.has_integrator(refPlugin)) continue;
      Map<String, DM> args = dmd("x0", dm(1, 0.5), "p", dm(0.5, 0.3));
      if (sys.containsKey("z")) {
        // Consistent z0 for the reference: the real root of z^3 + z - 1.3 (numpy.roots in python)
        double zr = 1;
        for (int i = 0; i < 50; i++) zr -= (zr * zr * zr + zr - 1.3) / (3 * zr * zr + 1);
        args.put("z0", dm(zr));
      }
      for (Map<String, Object> o : List.of(opts("step_size", 0.04), opts("reltol", 1e-6, "abstol", 1e-6))) {
        List<Double> nSteps = new ArrayList<>();
        for (double[] grid : new double[][] {{2.}, linspace(0.002, 2, 1000)}) {
          Function F = ca.integrator("F", "rosenbrock", sys, 0, grid, with(o, "control_group_size", new long[] {grid.length}));
          Map<String, DM> a = new LinkedHashMap<>(args);
          a.put("u", dm(0.2));
          Map<String, DM> r = F.call(a);
          nSteps.add(num(F.stats().get("n_steps")));
          Map<String, DM> ref = ca.integrator("ref", refPlugin, sys, 0, grid, opts("reltol", 1e-10, "abstol", 1e-10)).call(a);
          for (String k : List.of("xf", "zf", "qf")) chk(ref.get(k), r.get(k), 5, k);
        }
        assertTrue(nSteps.get(1) <= nSteps.get(0) + 2);
        double[] grid = {0.13, 0.31, 0.37, 0.52, 0.9, 1.15, 1.2, 1.71, 2.0};
        for (double[] u0 : new double[][] {{0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1}, {0.1, 0.1, 0.1, 0.1, 0.3, 0.3, 0.3, 0.3, 0.3}}) {
          Function F = ca.integrator("F", "rosenbrock", sys, 0, grid, o);
          MX v = MX.sym("v", 13);
          Map<String, MX> r = sys.containsKey("z")
              ? F.call(mxd("x0", sl(v, 0, 2), "z0", 0.5, "p", sl(v, 2, 4), "u", sl(v, 4, 13).T()))
              : F.call(mxd("x0", sl(v, 0, 2), "p", sl(v, 2, 4), "u", sl(v, 4, 13).T()));
          MX out = vm(ca.vec(r.get("xf")), ca.vec(r.get("zf")), ca.vec(r.get("qf")));
          double[] v0 = new double[13];
          v0[0] = 1;
          v0[1] = 0.5;
          v0[2] = 0.5;
          v0[3] = 0.3;
          System.arraycopy(u0, 0, v0, 4, 9);
          checkJacobian(out, v, dm(v0), 12);
          checkHessian(ca.sumsqr(out), v, dm(v0), 10);
        }
      }
    }
  }

  @Test
  public void test_rosenbrock_codegen() {
    message("rosenbrock: code generation, including derivatives");
    SX x = SX.sym("x", 2), z = SX.sym("z"), p = SX.sym("p");
    Map<String, SX> dae = sxd("x", x, "z", z, "p", p, "ode", vs(x.get(1), p.unaryMinus().times(z)),
        "alg", z.power(3).plus(z).minus(x.get(0)), "quad", x.get(0).power(2));
    Object[][] cases = {{opts("reltol", 1e-8, "abstol", 1e-8, "max_tape_steps", 1000), new double[] {0.5, 1}},
        {opts("step_size", 0.05), new double[] {0.5, 1}}, {opts("step_size", 0.03), linspace(0.07, 1, 15)},
        {opts("reltol", 1e-6, "abstol", 1e-6, "max_tape_steps", 1000), linspace(0.07, 1, 15)}};
    for (Object[] c : cases) {
      double[] grid = (double[]) c[1];
      @SuppressWarnings("unchecked")
      Function F = ca.integrator("F", "rosenbrock", dae, 0, grid, (Map<String, Object>) c[0]);
      MX v = MX.sym("v", 3);
      Map<String, MX> r = F.call(mxd("x0", sl(v, 0, 2), "p", v.get(2)));
      MX obj = ca.sumsqr(r.get("xf")).plus(ca.sumsqr(r.get("qf")));
      List<MX> outs = new ArrayList<>(List.of(r.get("xf"), r.get("zf"), r.get("qf"), ca.gradient(obj, v)));
      if (grid.length > 2) outs.add(ca.jtimes(ca.jtimes(obj, v, MX.eye(3)).T(), v, MX.eye(3)));
      Function G = new Function("G", List.of(v), outs);
      checkCodegen(G, List.of(dm(0.4, -0.3, 2)), 12);
    }
  }

  @Test
  public void test_erk() {
    message("erk: order 4, dense output, steps independent of the outputs");
    SX t = SX.sym("t"), x = SX.sym("x", 2), p = SX.sym("p"), u = SX.sym("u");
    Map<String, SX> dae = sxd("t", t, "x", x, "p", p, "u", u, "ode", vs(x.get(1),
        p.unaryMinus().times(x.get(0)).plus(u.times(ca.cos(t)))), "quad", x.get(0).power(2));
    double w = 1.3;
    java.util.function.DoubleFunction<DM> exact = tt -> dm(Math.cos(w * tt), -w * Math.sin(w * tt));
    List<double[]> err = new ArrayList<>();
    for (double h : new double[] {0.1, 0.05, 0.025}) {
      Map<String, DM> r = ca.integrator("F", "erk", dae, 0, new double[] {0.37, 2.0}, opts("method", "rk4", "step_size", h))
          .call(dmd("x0", dm(1, 0), "p", w * w));
      err.add(new double[] {ca.norm_inf(r.get("xf").get(new Slice(), 1).minus(exact.apply(2.0))).to_double(),
          ca.norm_inf(r.get("xf").get(new Slice(), 0).minus(exact.apply(0.37))).to_double()});
    }
    for (int i = 0; i + 1 < err.size(); i++) {
      assertTrue(Math.log(err.get(i)[0] / err.get(i + 1)[0]) / Math.log(2) > 3.8);
      assertTrue(Math.log(err.get(i)[1] / err.get(i + 1)[1]) / Math.log(2) > 3.4);
    }
    Map<String, DM> args = dmd("x0", dm(1, 0.2), "p", 2, "u", dm(new double[][] {{0.3, -0.1}}));
    Map<String, DM> r = ca.integrator("F", "erk", dae, 0, new double[] {0.5, 1.0}, opts("method", "rk4", "step_size", 0.1))
        .call(args);
    Map<String, DM> ref = ca.integrator("F", "rk", dae, 0, new double[] {0.5, 1.0}, opts("number_of_finite_elements", 10))
        .call(args);
    for (String k : List.of("xf", "qf")) chk(ref.get(k), r.get(k), 14, k);
    for (int n : new int[] {10, 1000}) {
      Function F = ca.integrator("F", "erk", dae, 0, linspace(0, 10, n),
          opts("step_size", 0.05, "control_group_size", new long[] {n}));
      F.call(dmd("x0", dm(1, 0), "p", w * w));
      assertEquals(200L, F.stats().get("n_steps"));
    }
  }

  @Test
  public void test_erk_derivatives() {
    message("erk: derivatives are exact derivatives of the scheme");
    SX t = SX.sym("t"), x = SX.sym("x", 2), p = SX.sym("p"), u = SX.sym("u");
    Map<String, SX> dae = sxd("t", t, "x", x, "p", p, "u", u, "ode", vs(x.get(1),
        p.unaryMinus().times(x.get(0).power(3)).plus(u.times(ca.cos(t)))), "quad", x.get(0).power(2).plus(u.times(x.get(1))));
    double[] grid = {0, 0.13, 0.31, 0.37, 1.0, 1.15, 2.0};
    for (double[] u0 : new double[][] {{0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1}, {0.1, 0.1, 0.1, 0.3, 0.3, 0.3, 0.3},
        linspace(0.1, 0.4, 7)}) {
      Function F = ca.integrator("F", "erk", dae, 0, grid, opts("step_size", 0.1));
      MX v = MX.sym("v", 10);
      Map<String, MX> r = F.call(mxd("x0", sl(v, 0, 2), "p", v.get(2), "u", sl(v, 3, 10).T()));
      MX out = vm(ca.vec(r.get("xf")), ca.vec(r.get("qf")));
      double[] v0a = new double[10];
      v0a[0] = 1;
      v0a[1] = 0.5;
      v0a[2] = 1.7;
      System.arraycopy(u0, 0, v0a, 3, 7);
      DM v0 = dm(v0a);
      checkJacobian(out, v, v0, 13);
      MX obj = ca.sumsqr(out);
      checkHessian(obj, v, v0, 12);
      DM g = new Function("g", List.of(v), List.of(ca.gradient(obj, v))).call(List.of(v0)).get(0);
      chk(g, fdJacobian(new Function("f", List.of(v), List.of(obj)), v0).T(), 7, "FD");
    }
  }

  @Test
  public void test_erk_codegen() {
    message("erk: code generation, including derivatives");
    SX x = SX.sym("x", 2), p = SX.sym("p"), u = SX.sym("u");
    Map<String, SX> dae = sxd("x", x, "p", p, "u", u, "ode", vs(x.get(1), p.unaryMinus().times(ca.sin(x.get(0))).plus(u)),
        "quad", x.get(0).power(2));
    Object[][] cases = {{new double[] {0.5, 1}, opts("step_size", 0.04)}, {linspace(0.07, 1, 15), opts("step_size", 0.04)},
        {linspace(0.07, 1, 15), opts("method", "rk4", "step_size", 0.1)},
        {linspace(0.07, 1, 15), opts("method", "rk45", "reltol", 1e-6, "abstol", 1e-6, "max_tape_steps", 200)}};
    for (Object[] c : cases) {
      double[] grid = (double[]) c[0];
      @SuppressWarnings("unchecked")
      Function F = ca.integrator("F", "erk", dae, 0, grid, (Map<String, Object>) c[1]);
      int n = 3 + grid.length;
      MX v = MX.sym("v", n);
      Map<String, MX> r = F.call(mxd("x0", sl(v, 0, 2), "p", v.get(2), "u", sl(v, 3, n).T()));
      MX obj = ca.sumsqr(r.get("xf")).plus(ca.sumsqr(r.get("qf")));
      MX S = new MX(DM.eye(n).get(new Slice(), dm(0, 2, 3)));
      List<MX> outs = List.of(r.get("xf"), r.get("qf"), ca.gradient(obj, v), ca.jtimes(ca.jtimes(obj, v, S).T(), v, S));
      Function G = new Function("G", List.of(v), outs);
      double[] v0 = new double[n];
      v0[0] = 0.4;
      v0[1] = -0.3;
      v0[2] = 2;
      for (int i = 0; i < grid.length; i++) v0[3 + i] = 0.1 * (i % 3);
      checkCodegen(G, List.of(dm(v0)), 12);
    }
  }

  @Test
  @SuppressWarnings("unchecked")
  public void test_control_groups() {
    message("erk, rosenbrock: one control per control group, derivatives for each control");
    SX t = SX.sym("t"), x = SX.sym("x", 2), z = SX.sym("z"), p = SX.sym("p"), u = SX.sym("u");
    Map<String, SX> ode = sxd("t", t, "x", x, "p", p, "u", u, "ode", vs(x.get(1),
        p.unaryMinus().times(x.get(0).power(3)).plus(u.times(ca.cos(t)))), "quad", x.get(0).power(2).plus(u.times(x.get(1))));
    Map<String, SX> dae = sxd("t", t, "x", x, "z", z, "p", p, "u", u, "ode", vs(x.get(1),
        p.unaryMinus().times(x.get(0)).plus(u.times(z))), "alg", z.power(3).plus(z).minus(u).minus(x.get(0).times(0.5)),
        "quad", x.get(0).power(2).plus(u.times(z)));
    double[] grid = {0.2, 0.3, 0.5, 1.0, 1.2, 2.0};
    Object[][] cases = {{"erk", opts("method", "rk4", "step_size", 0.1)}, {"rosenbrock", opts("step_size", 0.1)},
        {"erk", opts("reltol", 1e-10, "abstol", 1e-10, "max_tape_steps", 500)},
        {"rosenbrock", opts("reltol", 1e-10, "abstol", 1e-10, "max_tape_steps", 500)}};
    for (Object[] c : cases) {
      String plugin = (String) c[0];
      Map<String, Object> o = (Map<String, Object>) c[1];
      Map<String, SX> d = plugin.equals("rosenbrock") ? dae : ode;
      Map<String, DM> r = ca.integrator("F", plugin, d, 0, grid, o)
          .call(dmd("u", dm(new double[][] {{0.1, 0.1, 0.3, 0.2, 0.2, 0.2}}), "x0", dm(1, 0.5), "p", 1.7));
      Function F = ca.integrator("F", plugin, d, 0, grid, with(o, "control_group_size", new long[] {2, 1, 3}));
      Map<String, DM> rGroup = F.call(dmd("u", dm(new double[][] {{0.1, 0.3, 0.2}}), "x0", dm(1, 0.5), "p", 1.7));
      for (String k : List.of("xf", "zf", "qf")) chk(r.get(k), rGroup.get(k), 6, k);
      for (long[] cgs : new long[][] {null, {2, 1, 3}}) {
        long[] groups = cgs == null ? new long[] {1, 1, 1, 1, 1, 1} : cgs;
        F = ca.integrator("F", plugin, d, 0, grid, with(o, "control_group_size", groups));
        int n = 3 + (int) F.size2_in(3);
        MX v = MX.sym("v", n);
        Map<String, MX> rv = F.call(mxd("x0", sl(v, 0, 2), "p", v.get(2), "u", sl(v, 3, n).T()));
        MX out = vm(ca.vec(rv.get("xf")), ca.vec(rv.get("zf")), ca.vec(rv.get("qf")));
        Function G = new Function("G", List.of(v), List.of(out, ca.jtimes(out, v, MX.eye(n)),
            ca.jtimes(out, v, MX.eye(out.numel()), true).T()));
        double[] v0a = new double[n];
        v0a[0] = 1;
        v0a[1] = 0.5;
        v0a[2] = 1.7;
        for (int i = 3; i < n; i++) v0a[i] = 0.2;
        DM v0 = dm(v0a);
        List<DM> Go = G.call(List.of(v0));
        DM Jfd = fdJacobian(G, v0);
        chk(Jfd, Go.get(1), 7, "fwd vs FD");
        chk(Jfd, Go.get(2), 7, "adj vs FD");
        if (cgs != null) checkCodegen(G, List.of(v0), 12);
        long[] start = new long[groups.length + 1];
        for (int i = 0; i < groups.length; i++) start[i + 1] = start[i] + groups[i];
        for (int wgt : new int[] {0, 1}) {
          Sparsity J = new Function("J", List.of(v), List.of(out), opts("ad_weight_sp", wgt)).jac_sparsity(0, 0);
          double[][] dep = new DM(J, 1).full();
          for (int k = 0; k < 6; k++)
            for (int j = 0; j < n - 3; j++) {
              boolean any = dep[2 * k][3 + j] != 0 || dep[2 * k + 1][3 + j] != 0;
              assertEquals(start[j] <= k, any);
            }
        }
      }
    }
    assertInException("must sum to the number of output times",
        () -> ca.integrator("F", "erk", ode, 0, grid, opts("control_group_size", new long[] {2, 1})));
  }

  @Test
  public void test_erk_rk45() {
    message("erk: pair of order 5(4), fixed and adaptive steps");
    SX t = SX.sym("t"), x = SX.sym("x", 2), p = SX.sym("p");
    Map<String, SX> dae = sxd("t", t, "x", x, "p", p, "ode", vs(x.get(1), p.unaryMinus().times(x.get(0))),
        "quad", x.get(0).power(2));
    double w = 1.3;
    java.util.function.DoubleFunction<double[]> exact = tt -> new double[] {Math.cos(w * tt), -w * Math.sin(w * tt)};
    List<Double> err = new ArrayList<>();
    for (double h : new double[] {0.1, 0.05, 0.025}) {
      double[] xf = ca.integrator("F", "erk", dae, 0, 2, opts("method", "rk45", "step_size", h))
          .call(dmd("x0", dm(1, 0), "p", w * w)).get("xf").nonzeros();
      double[] ex = exact.apply(2.0);
      err.add(Math.max(Math.abs(xf[0] - ex[0]), Math.abs(xf[1] - ex[1])));
    }
    for (int i = 0; i + 1 < err.size(); i++) assertTrue(Math.log(err.get(i) / err.get(i + 1)) / Math.log(2) > 4.7);
    double[] grid = linspace(0.1, 10, 50);
    for (double tol : new double[] {1e-6, 1e-9}) {
      Function F = ca.integrator("F", "erk", dae, 0, grid, opts("method", "rk45", "reltol", tol, "abstol", tol));
      double[][] xf = F.call(dmd("x0", dm(1, 0), "p", w * w)).get("xf").full();
      double e = 0;
      for (int j = 0; j < grid.length; j++) {
        double[] ex = exact.apply(grid[j]);
        for (int i = 0; i < 2; i++) e = Math.max(e, Math.abs(xf[i][j] - ex[i]));
      }
      assertTrue(e < 10 * tol);
      assertTrue(e > 0.01 * tol);
      Map<String, Object> st = F.stats();
      assertEquals(1 + 6 * (num(st.get("n_steps")) + num(st.get("n_reject"))), num(st.get("n_rhs")));
    }
    Map<String, SX> dae2 = sxd("x", x, "ode", vs(ca.times(-50, x.get(0).minus(ca.cos(x.get(1)))), 1));
    Function F = ca.integrator("F", "erk", dae2, 0, 1, opts("method", "rk45", "step0", 0.5));
    Map<String, DM> r = F.call(dmd("x0", dm(0, 0)));
    Map<String, DM> ref = ca.integrator("F", "cvodes", dae2, 0, 1, opts("reltol", 1e-12, "abstol", 1e-12))
        .call(dmd("x0", dm(0, 0)));
    chk(ref.get("xf"), r.get("xf"), 5, "rejections");
    assertTrue(num(F.stats().get("n_reject")) > 0);
    assertInException("unknown method", () -> ca.integrator("F", "erk", dae, 0, 1, opts("method", "rk5")));
  }

  @Test
  public void test_erk_rk45_derivatives() {
    message("erk: adaptive steps, derivatives of the scheme on the nominal mesh");
    SX t = SX.sym("t"), x = SX.sym("x", 2), p = SX.sym("p"), u = SX.sym("u");
    Map<String, SX> dae = sxd("t", t, "x", x, "p", p, "u", u, "ode", vs(x.get(1),
        p.unaryMinus().times(x.get(0).power(3)).plus(u.times(ca.cos(t)))), "quad", x.get(0).power(2).plus(u.times(x.get(1))));
    double[] grid = {0, 0.13, 0.31, 0.37, 1.0, 1.15, 2.0};
    for (double[] u0 : new double[][] {{0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1}, {0.1, 0.1, 0.1, 0.3, 0.3, 0.3, 0.3}}) {
      Function F = ca.integrator("F", "erk", dae, 0, grid, opts("method", "rk45", "reltol", 1e-5, "abstol", 1e-5));
      MX v = MX.sym("v", 10);
      Map<String, MX> r = F.call(mxd("x0", sl(v, 0, 2), "p", v.get(2), "u", sl(v, 3, 10).T()));
      MX out = vm(ca.vec(r.get("xf")), ca.vec(r.get("qf")));
      double[] v0a = new double[10];
      v0a[0] = 1;
      v0a[1] = 0.5;
      v0a[2] = 1.7;
      System.arraycopy(u0, 0, v0a, 3, 7);
      checkJacobian(out, v, dm(v0a), 13);
      checkHessian(ca.sumsqr(out), v, dm(v0a), 11);
    }
  }

  private static final double G_ = 9.81;

  /** Height, velocity and integral of the height of the bouncing ball */
  private static double[] exactBounce(double h0, double p0, double tt) {
    double g = G_, tb = 0, q = 0, h = h0, v = 0;
    Double ub = null;
    while (true) {
      double tau = ub == null ? Math.sqrt(2 * h / g) : 2 * ub / g;
      if (tb + tau > tt) break;
      q += ub == null ? h * tau - g * Math.pow(tau, 3) / 6 : ub * tau * tau / 2 - g * Math.pow(tau, 3) / 6;
      v = ub == null ? -g * tau : ub - g * tau;
      tb = tb + tau;
      ub = -p0 * v;
      h = 0.;
    }
    double tau = tt - tb, u0 = ub == null ? 0. : ub;
    return new double[] {h + u0 * tau - g * tau * tau / 2, u0 - g * tau, q + h * tau + u0 * tau * tau / 2 - g * Math.pow(tau, 3) / 6};
  }

  @Test
  public void test_events_bouncing_ball() {
    message("zero-crossing events: values, forward sensitivities, code generation");
    SX x = SX.sym("x", 2), p = SX.sym("p");
    Map<String, SX> dae = sxd("x", x, "p", p, "ode", vs(x.get(1), -G_), "quad", x.get(0), "zero", x.get(0));
    Function transition = bounceTransition();
    double[] grid = linspace(0.1, 4, 40);
    double h0 = 2., p0 = 0.8;
    double[][] xfRef = new double[2][grid.length];
    double[] qfRef = new double[grid.length];
    for (int i = 0; i < grid.length; i++) {
      double[] e = exactBounce(h0, p0, grid[i]);
      xfRef[0][i] = e[0];
      xfRef[1][i] = e[1];
      qfRef[i] = e[2];
    }
    for (EvIntg E : eventIntegrators) {
      if (!E.events) continue;
      Function F = ca.integrator("F", E.name, dae, 0, grid, with(E.options, "transition", transition, "event_tol", 1e-12));
      Map<String, DM> r = F.call(dmd("x0", dm(h0, 0), "p", p0));
      chk(r.get("xf"), new DM(xfRef), E.digits, E.name);
      chk(r.get("qf"), dm(qfRef).T(), E.digits, E.name);
      if (F.stats().containsKey("n_events")) assertEquals(5L, F.stats().get("n_events"));
      MX v = MX.sym("v", 3);
      Map<String, MX> rv = F.call(mxd("x0", sl(v, 0, 2), "p", v.get(2)));
      MX out = vm(rv.get("xf").get(new Slice(), -1), rv.get("qf").get(-1));
      Function G = new Function("G", List.of(v), List.of(out, ca.jacobian(out, v)));
      DM v0 = dm(h0, 0.3, p0);
      chk(G.call(List.of(v0)).get(1), fdJacobian(G, v0), Math.min(6, E.digits), E.name);
      if (E.codegen) checkCodegen(G, List.of(v0), 10);
    }
  }

  @Test
  public void test_events_multiple() {
    message("zero-crossing events: several functions, transition by event index");
    double g = G_, ceiling = 1.5;
    SX x = SX.sym("x", 2), p = SX.sym("p"), index = SX.sym("index"), t = SX.sym("t");
    Map<String, SX> dae = sxd("x", x, "p", p, "ode", vs(x.get(1), -g), "zero", vs(x.get(0), ca.minus(ceiling, x.get(0))));
    Function transition = new Function("transition", List.of(index, t, x, new SX(0, 1), p, new SX(0, 1)),
        List.of(vs(x.get(0), ca.if_else(ca.eq(index, new SX(0)), p.unaryMinus(), new SX(-1)).times(x.get(1))), new SX(0, 1)),
        List.of("index", "t", "x", "z", "p", "u"), List.of("post_x", "post_z"));
    DM v0 = dm(1, 4, 0.8);
    double[] v0a = v0.nonzeros();
    double[] grid = linspace(0.1, 3, 30);
    double[][] xr = new double[2][grid.length];
    for (int i = 0; i < grid.length; i++) {
      double tt = grid[i], t0 = 0, h = v0a[0], v = v0a[1];
      while (true) {
        double bestTau = Double.POSITIVE_INFINITY, bestLevel = 0;
        for (double[] ls : new double[][] {{0., -1}, {ceiling, 1}}) {
          double level = ls[0], sign = ls[1];
          double disc = v * v - 2 * g * (level - h);
          if (disc < 0) continue;
          for (double tau : new double[] {(v - Math.sqrt(disc)) / g, (v + Math.sqrt(disc)) / g})
            if (tau > 1e-12 && sign * (v - g * tau) > 0 && (tau < bestTau || (tau == bestTau && level < bestLevel))) {
              bestTau = tau;
              bestLevel = level;
            }
        }
        if (t0 + bestTau > tt) break;
        v = (v - g * bestTau) * (bestLevel == 0 ? -v0a[2] : -1);
        t0 = t0 + bestTau;
        h = bestLevel;
      }
      double tau = tt - t0;
      xr[0][i] = h + v * tau - g * tau * tau / 2;
      xr[1][i] = v - g * tau;
    }
    DM xfRef = new DM(xr);
    for (EvIntg E : eventIntegrators) {
      if (!E.events) continue;
      Function F = ca.integrator("F", E.name, dae, 0, grid, with(E.options, "transition", transition, "event_tol", 1e-12));
      MX v = MX.sym("v", 3);
      MX xf = F.call(mxd("x0", sl(v, 0, 2), "p", v.get(2))).get("xf");
      Function G = new Function("G", List.of(v), List.of(ca.vec(xf), ca.jacobian(ca.vec(xf), v)));
      chk(G.call(List.of(v0)).get(0), ca.vec(xfRef), E.digits, E.name);
      chk(G.call(List.of(v0)).get(1), fdJacobian(G, v0), Math.min(6, E.digits), E.name);
      if (E.codegen) checkCodegen(G, List.of(v0), 10);
    }
  }

  @Test
  public void test_events_dae() {
    message("zero-crossing events of a DAE, algebraic states consistent after the transition");
    double g = G_;
    SX x = SX.sym("x", 2), z = SX.sym("z"), p = SX.sym("p");
    Map<String, SX> dae = sxd("x", x, "z", z, "p", p, "ode", vs(z, -g), "alg", z.minus(x.get(1)), "zero", x.get(0));
    Function transition = bounceTransition(1);
    double[] grid = linspace(0.1, 3, 30);
    Function ref = ca.integrator("ref", "erk", sxd("x", x, "p", p, "ode", vs(x.get(1), -g), "zero", x.get(0)), 0, grid,
        opts("transition", bounceTransition(), "event_tol", 1e-12, "abstol", 1e-11, "reltol", 1e-11));
    DM v0 = dm(2, 0.3, 0.8);
    DM xfRef = ref.call(dmd("x0", v0.get(new Slice(0, 2)), "p", v0.get(2))).get("xf");
    for (EvIntg E : eventIntegrators) {
      if (!E.eventsDae) continue;
      Function F = ca.integrator("F", E.name, dae, 0, grid, with(E.options, "transition", transition, "event_tol", 1e-12));
      MX v = MX.sym("v", 3);
      Map<String, MX> r = F.call(mxd("x0", sl(v, 0, 2), "z0", v.get(1), "p", v.get(2)));
      MX out = vm(ca.vec(r.get("xf")), ca.vec(r.get("zf")));
      Function G = new Function("G", List.of(v), List.of(out, ca.jacobian(out, v)));
      chk(G.call(List.of(v0)).get(0), ca.vertcat(List.of(ca.vec(xfRef), xfRef.get(1, new Slice()).T())), E.digits, E.name);
      chk(G.call(List.of(v0)).get(1), fdJacobian(G, v0), Math.min(6, E.digits), E.name);
      if (E.codegen) checkCodegen(G, List.of(v0), 10);
    }
  }

  @Test
  public void test_events_breaking_spring() {
    message("zero-crossing events with controls and parameters");
    SX x = SX.sym("x", 3), p = SX.sym("p", 3), u = SX.sym("u"), index = SX.sym("index"), t = SX.sym("t");
    SX f = ca.if_else(ca.lt(x.get(2), new SX(0.5)), p.get(1).unaryMinus().times(x.get(1)).plus(u), new SX(0));
    Map<String, SX> dae = sxd("x", x, "p", p, "u", u, "ode", vs(f.minus(p.get(2).times(x.get(0))).div(p.get(0)), x.get(0), 0),
        "zero", ca.minus(2, x.get(1)));
    Function transition = new Function("transition", List.of(index, t, x, new SX(0, 1), p, u),
        List.of(vs(x.get(0), x.get(1), 1), new SX(0, 1)), List.of("index", "t", "x", "z", "p", "u"),
        List.of("post_x", "post_z"));
    double[] grid = linspace(0.04, 4, 100);
    double[] d = new double[100];
    for (int i = 20; i < 30; i++) d[i] = -5;
    for (int i = 40; i < 50; i++) d[i] = 15;
    DM x0 = dm(-5, -1, 0), p0 = dm(1, 2, 0.1), dd = dm(d);
    MX v = MX.sym("v", 3);
    java.util.function.Function<Function, Function> outputs = F -> {
      MX xf = F.call(mxd("x0", x0, "p", v, "u", dd)).get("xf");
      return new Function("G", List.of(v), List.of(xf, ca.jacobian(xf, v)));
    };
    List<DM> ref = outputs.apply(ca.integrator("ref", "erk", dae, 0, grid, opts("transition", transition,
        "event_tol", 1e-12, "abstol", 1e-12, "reltol", 1e-12))).call(List.of(p0));
    for (EvIntg E : eventIntegrators) {
      if (!E.events) continue;
      Function F = ca.integrator("F", E.name, dae, 0, grid, with(E.options, "transition", transition, "event_tol", 1e-12));
      Function G = outputs.apply(F);
      List<DM> r = G.call(List.of(p0));
      chk(r.get(0), ref.get(0), E.digits - 1, E.name);
      chk(r.get(1), ref.get(1), E.digits - 2, E.name);
      F.call(dmd("x0", x0, "p", p0, "u", dd));
      if (F.stats().containsKey("n_events")) assertEquals(1L, F.stats().get("n_events"));
      if (E.codegen) checkCodegen(G, List.of(p0), 12);
    }
    Function F = ca.integrator("F", "erk", dae, 0, grid, opts("transition", transition, "max_events", 0));
    assertInException("max_events", () -> F.call(dmd("x0", x0, "p", p0, "u", dd)));
  }

  @Test
  public void test_events_adjoint() {
    message("zero-crossing events: adjoint and forward-over-adjoint sensitivities");
    double g = G_;
    SX x = SX.sym("x", 2), z = SX.sym("z"), p = SX.sym("p");
    Map<String, SX> ode = sxd("x", x, "p", p, "ode", vs(x.get(1), -g), "quad", x.get(0), "zero", x.get(0));
    Map<String, SX> dae = sxd("x", x, "z", z, "p", p, "ode", vs(z, -g), "alg", z.minus(x.get(1)), "quad", x.get(0),
        "zero", x.get(0));
    double[] grid = linspace(0.1, 3, 30);
    // numpy.random.default_rng(0).uniform(0.5, 1.5, (2, 30)) in python: java.util.Random here
    Random rng = new Random(0);
    double[][] wa = new double[2][30];
    for (int i = 0; i < 2; i++) for (int j = 0; j < 30; j++) wa[i][j] = 0.5 + rng.nextDouble();
    DM wts = new DM(wa);
    for (EvIntg E : eventIntegrators) {
      if (!E.eventsAdjoint) continue;
      List<Object[]> problems = new ArrayList<>();
      problems.add(new Object[] {ode, bounceTransition()});
      if (E.eventsDae) problems.add(new Object[] {dae, bounceTransition(1)});
      for (Object[] pt : problems) {
        @SuppressWarnings("unchecked")
        Function F = ca.integrator("F", E.name, (Map<String, SX>) pt[0], 0, grid, with(E.options, "transition", pt[1],
            "event_tol", 1e-12));
        MX v = MX.sym("v", 3);
        Map<String, MX> r = F.call(mxd("x0", sl(v, 0, 2), "z0", v.get(1), "p", v.get(2)));
        MX obj = ca.sum1(ca.sum2(ca.times(wts, r.get("xf")))).plus(ca.sumsqr(r.get("qf")));
        Function G = new Function("G", List.of(v), List.of(ca.jtimes(obj, v, MX.ones(1, 1), true),
            ca.jtimes(obj, v, MX.eye(3)).T(), (MX) ca.hessian(obj, v).get(0)));
        DM v0 = dm(2, 0.3, 0.8);
        List<DM> Go = G.call(List.of(v0));
        chk(Go.get(0), Go.get(1), 9, E.name);
        DM fd = fdJacobian(G, v0, 1e-4).times(4).minus(fdJacobian(G, v0, 2e-4)).div(3);
        chk(Go.get(2), fd, 6, E.name);
        if (E.codegen) checkCodegen(G, List.of(v0), 9);
      }
    }
  }

  @Test
  public void test_erk_tape() {
    message("erk: adaptive steps and adjoints, bounded tape");
    SX x = SX.sym("x", 2);
    Map<String, SX> dae = sxd("x", x, "ode", vs(x.get(1), x.get(0).unaryMinus()));
    for (int nTape : new int[] {0, 3, 1000}) {
      Function F = ca.integrator("F", "erk", dae, 0, 10, opts("method", "rk45", "max_tape_steps", nTape));
      MX x0 = MX.sym("x0", 2);
      Function G = new Function("G", List.of(x0), List.of(ca.jtimes(ca.sumsqr(F.call(mxd("x0", x0)).get("xf")), x0,
          MX.ones(1, 1), true)));
      if (nTape == 3) {
        assertInException("increase max_tape_steps", () -> G.call(List.of(dm(1, 0.5))));
      } else {
        chk(G.call(List.of(dm(1, 0.5))).get(0), dm(2, 1), 5, "adjoint");
      }
      if (nTape == 0)
        assertInException("need max_tape_steps for code generation", () -> new CodeGenerator("tape_cg.c").add(G));
    }
  }
}

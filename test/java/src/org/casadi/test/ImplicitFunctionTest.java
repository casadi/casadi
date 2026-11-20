package org.casadi.test;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
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
import org.junit.jupiter.api.Test;

/** Port of test/python/implicitfunction.py */
public class ImplicitFunctionTest extends CasadiTestCase {

  private static final boolean KNOWN_BUGS = System.getenv("CASADI_JAVA_KNOWN_BUGS") != null;

  private static final class Solver {
    final String name;
    final Map<String, Object> options;
    final boolean codegen;

    Solver(String name, Map<String, Object> options, boolean codegen) {
      this.name = name;
      this.options = options;
      this.codegen = codegen;
    }
  }

  private static final List<Solver> solvers = new ArrayList<>();

  static {
    if (ca.has_rootfinder("kinsol")) solvers.add(new Solver("kinsol", opts("abstol", 1e-10), false));
    if (ca.has_nlpsol("ipopt"))
      solvers.add(new Solver("nlpsol", opts("nlpsol", "ipopt", "nlpsol_options",
          opts("print_time", false, "ipopt", opts("print_level", 0))), false));
    solvers.add(new Solver("newton", opts(), false));
    solvers.add(new Solver("fast_newton", opts(), true));
  }

  // ---------------------------------------------------------------- helpers

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

  private static List<DM> ins(double... v) {
    List<DM> r = new ArrayList<>();
    for (double e : v) r.add(new DM(e));
    return r;
  }

  private static SX vs(SX... a) {
    return ca.vertcat(Arrays.asList(a));
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
    assertTrue(diff < 0.5 * Math.pow(10, -places), msg + " " + first + " <-> " + second + " [places=" + places + "]");
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

  /** checkfunction of helpers.py with jacobian/hessian on; evals as given, digits_sens = digits */
  private void checkfun(Function trial, Function solution, List<DM> inputs, int digits, int evals, String fail) {
    checkfun(trial, solution, inputs, digits, true, true, evals, fail);
  }

  private void checkfun(Function trial, Function solution, List<DM> inputs, int digits, boolean jacobian,
      boolean hessian, int evals, String fail) {
    for (int rep = 0; rep < 2; rep++) {
      List<DM> to = trial.call(inputs), so = solution.call(inputs);
      assertEquals(trial.n_out(), solution.n_out(), fail + ": n_out");
      assertEquals(trial.n_in(), solution.n_in(), fail + ": n_in");
      for (int i = 0; i < trial.n_out(); i++) {
        if (trial.sparsity_out(i).is_empty() || solution.sparsity_out(i).is_empty()) continue;
        chk(to.get(i), so.get(i), digits, fail + ": output(" + i + ": '" + trial.name_out(i) + "')");
      }
    }
    if (jacobian)
      for (long i = 0; i < trial.n_in(); i++) {
        if (trial.sparsity_in(i).is_empty() || solution.sparsity_in(i).is_empty()) continue;
        for (long j = 0; j < trial.n_out(); j++) {
          Function tj = jacobianOld(trial, i, j), sj = jacobianOld(solution, i, j);
          assertEquals(trial.n_in(), tj.n_in());
          assertEquals(trial.n_out() + 1, tj.n_out());
          assertEquals(solution.n_in(), sj.n_in());
          assertEquals(solution.n_out() + 1, sj.n_out());
          checkfun(tj, sj, inputs, digits, false, false, 0, "(" + fail + ").jacobian_old(" + i + "," + j + ")");
        }
      }
    if (hessian)
      for (long i = 0; i < trial.n_in(); i++) {
        if (trial.sparsity_in(i).is_empty() || solution.sparsity_in(i).is_empty()) continue;
        for (long j = 0; j < trial.n_out(); j++) {
          if (!(trial.sparsity_out(j).is_scalar() && solution.sparsity_out(j).is_scalar())) continue;
          Function th = hessianOld(trial, i, j), sh = hessianOld(solution, i, j);
          assertEquals(trial.n_in(), th.n_in());
          assertEquals(trial.n_out() + 2, th.n_out());
          assertEquals(solution.n_in(), sh.n_in());
          assertEquals(solution.n_out() + 2, sh.n_out());
          checkfun(th, sh, inputs, digits, false, false, 0, "(" + fail + ").hessian_old(" + i + "," + j + ")");
        }
      }
    if (evals > 0) checkEvals(trial, solution, inputs, evals, digits, fail);
  }

  private void checkfun(Function trial, Function solution, List<DM> inputs, int digits) {
    checkfun(trial, solution, inputs, digits, 2, "");
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
    List<MX> fseeds = new ArrayList<>(), aseeds = new ArrayList<>(), in = new ArrayList<>();
    for (Sparsity s : fsp) fseeds.add(MX.sym("f", s));
    for (Sparsity s : asp) aseeds.add(MX.sym("a", s));
    for (Sparsity s : isp) in.add(MX.sym("i", s));
    List<MX> res = f.call(in, !f.is_a("SXFunction"));
    Map<String, Object> o = opts("helper_options", opts("is_diff_in", f.is_diff_in(), "is_diff_out", f.is_diff_out()));
    List<MX> fwd = ca.forward(res, in, List.of(fseeds), o).get(0);
    List<MX> adj = ca.reverse(res, in, List.of(aseeds), o).get(0);
    List<MX> vin = new ArrayList<>(in), vout = new ArrayList<>(res);
    vin.addAll(fseeds);
    vin.addAll(aseeds);
    vout.addAll(fwd);
    vout.addAll(adj);
    return new Function("vf", vin, vout, opts("is_diff_in", cat(f.is_diff_in(), f.is_diff_in(), f.is_diff_out()),
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

  /** check_serialize of helpers.py, without the python-side serialization scheme validation */
  private static void checkSerialize(Function f, List<DM> inputs) {
    Function f2 = Function.deserialize(f.serialize(opts("debug", true)));
    List<DM> a = f.call(inputs), b = f2.call(inputs);
    for (int i = 0; i < f.n_out(); i++) chk(a.get(i), b.get(i), 16);
  }

  private static int cgCounter = 0;

  /** check_codegen of helpers.py: compiles with gcc only when slow tests run */
  private static void checkCodegen(Function f, List<DM> inputs) {
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
      for (int i = 0; i < f.n_out(); i++) chk(a.get(i), b.get(i), 15);
    } catch (java.io.IOException | InterruptedException e) {
      throw new RuntimeException(e);
    }
  }

  // ---------------------------------------------------------------- tests

  @Test
  public void test_linear() {
    for (Solver S : solvers) {
      message(S.name);
      SX x = SX.sym("x", 2);
      DM A_ = dm(new double[][] {{1, 2}, {3, 2.1}}), C_ = dm(new double[][] {{1.6, 2.1}, {1, 1.3}}), b_ = dm(0.7, 0.6);
      Function f = new Function("f", List.of(x), List.of(ca.mtimes(A_, x).minus(b_), ca.mtimes(C_, x)));
      Function solver = ca.rootfinder("solver", S.name, f, S.options);
      solver.call(ins(0));
      DM sol = ca.solve(A_, b_);
      Function refsol = new Function("refsol", List.of(x), List.of(new SX(sol), new SX(ca.mtimes(C_, sol))));
      checkfun(solver, refsol, ins(0), 10);
      if (S.name.contains("newton")) checkSerialize(solver, ins(0));
      if (S.codegen) checkCodegen(solver, ins(0));

      SX A = SX.sym("A", 2, 2), b = SX.sym("b", 2);
      f = new Function("f", List.of(x, A, b), List.of(ca.mtimes(A, x).minus(b)));
      solver = ca.rootfinder("solver", S.name, f, S.options);
      List<DM> solverIn = List.of(new DM(0.0), A_, b_);
      refsol = new Function("refsol", List.of(x, A, b), List.of(ca.solve(A, b)));
      checkfun(solver, refsol, solverIn, 10);
      if (S.name.contains("newton")) checkSerialize(solver, solverIn);
      if (S.codegen) checkCodegen(solver, solverIn);

      A = SX.sym("A", 2, 2);
      b = SX.sym("b", 2);
      f = new Function("f", List.of(x, A, b), List.of(ca.mtimes(A, x).minus(b), ca.mtimes(C_, x)));
      for (int adWeightSp : new int[] {0, 1})
        for (int adWeight : new int[] {0, 1}) {
          message(adWeight + " " + adWeightSp);
          solver = ca.rootfinder("solver", S.name, f, with(S.options, "ad_weight_sp", adWeightSp, "ad_weight", adWeight));
          refsol = new Function("refsol", List.of(x, A, b), List.of(ca.solve(A, b), ca.mtimes(C_, ca.solve(A, b))));
          checkfun(solver, refsol, solverIn, 10);
          if (S.name.contains("newton")) checkSerialize(solver, solverIn);
        }
      if (S.codegen) checkCodegen(solver, solverIn);
    }
  }

  @Test
  public void test_missing_symbols() {
    for (Solver S : solvers) {
      message(S.name);
      SX x = SX.sym("x"), p = SX.sym("p");
      Function f = new Function("f", List.of(x), List.of(ca.sin(x.plus(p))), opts("allow_free", true));
      assertInException("[p] are free", () -> ca.rootfinder("solver", S.name, f, S.options));
    }
  }

  @Test
  public void test_scalar1() {
    message("Scalar implicit problem, n=0");
    for (Solver S : solvers) {
      message(S.name);
      SX x = SX.sym("x");
      Function f = new Function("f", List.of(x), List.of(ca.sin(x)));
      Function solver = ca.rootfinder("solver", S.name, f, S.options);
      solver.call(ins(0));
      Function refsol = new Function("refsol", List.of(x), List.of(ca.ceil(x.div(Math.PI).minus(0.5)).times(Math.PI)));
      checkfun(solver, refsol, ins(6), 10);
    }
  }

  @Test
  public void test_unsolved_stats() {
    for (Solver S : solvers) {
      message(S.name);
      SX x = SX.sym("x");
      Function f = new Function("f", List.of(x), List.of(ca.sin(x)));
      Function solver = ca.rootfinder("solver", S.name, f, S.options);
      try {
        solver.stats();
      } catch (RuntimeException e) {
        // no segfault
      }
    }
  }

  @Test
  public void test_scalar2() {
    message("Scalar implicit problem, n=1");
    for (Solver S : solvers) {
      message(S.name);
      SX x = SX.sym("x"), y = SX.sym("y");
      double n = 0.2;
      Function f = new Function("f", List.of(y, x), List.of(x.minus(ca.asin(y))));
      Function solver = ca.rootfinder("solver", S.name, f, S.options);
      Function refsol = new Function("refsol", List.of(y, x), List.of(ca.sin(x)));
      checkfun(solver, refsol, ins(0, n), 6, 2, S.name);
    }
  }

  @Test
  public void test_scalar2_indirect() {
    for (Solver S : solvers) {
      message(S.name);
      SX x = SX.sym("x"), y = SX.sym("y");
      double n = 0.2;
      Function f = new Function("f", List.of(y, x), List.of(x.minus(ca.asin(y))));
      Function solver = ca.rootfinder("solver", S.name, f, S.options);
      MX X = MX.sym("X");
      MX R = solver.call(List.of(new MX(), X)).get(0);
      Function trial = new Function("trial", List.of(X), List.of(R));
      Function refsol = new Function("refsol", List.of(x), List.of(ca.sin(x)));
      checkfun(trial, refsol, ins(n), 6, 2, S.name);
    }
  }

  @Test
  public void test_large() {
    for (Solver S : solvers) {
      if (S.name.contains("kinsol")) continue;
      if (S.name.contains("newton")) continue;
      int N = 5;
      Sparsity s = Sparsity.lower(N);
      SX x = SX.sym("x", s), y = SX.sym("y", s);
      DM y0 = new DM(Sparsity.diag(N), 0.1);
      SX r = ca.mtimes(x.plus(y0), x.plus(y0).T()).minus(ca.mtimes(y.plus(y0), y.plus(y0).T())).get(s);
      Function f = new Function("f", List.of(y.get_nz(false, new Slice()), x.get_nz(false, new Slice())),
          List.of(r.get_nz(false, new Slice())));
      long[] cons = new long[(int) s.nnz()];
      Arrays.fill(cons, 1);
      Function solver = ca.rootfinder("options2", S.name, f, with(S.options, "constraints", cons));
      MX X = MX.sym("X", x.sparsity());
      MX R = solver.call(List.of(new MX(), X.get_nz(false, new Slice()))).get(0);
      Function trial = new Function("trial", List.of(X), List.of(R));
      double[] v = new double[(int) x.nnz()];
      for (int i = 0; i < v.length; i++) v[i] = Math.abs(Math.cos(i));
      DM trialIn = new DM(trial.sparsity_in(0), new DM(v));
      DM trialOut = trial.call(List.of(trialIn)).get(0);
      f.call(List.of(trialOut, trialIn.get_nz(false, new Slice())));
      f.call(List.of(trialIn.get_nz(false, new Slice()), trialIn.get_nz(false, new Slice())));
      Function refsol = new Function("refsol", List.of(X), List.of(X.get_nz(false, new Slice())));
      checkfun(trial, refsol, List.of(trialIn.get(0)), 6, 1, S.name);
    }
  }

  @Test
  public void test_vector2() {
    Assumptions.assumeTrue(KNOWN_BUGS, "known bug (set CASADI_JAVA_KNOWN_BUGS to include it)");
    message("Scalar implicit problem, n=1");
    for (Solver S : solvers) {
      message(S.name);
      SX x = SX.sym("x"), y = SX.sym("y", 2);
      DM y0 = dm(0.1, 0.4);
      SX yy = y.plus(y0);
      double n = 0.2;
      Function f = new Function("f", List.of(y, x),
          List.of(vs(x.minus(ca.asin(yy.get(0))), yy.get(1).power(2).minus(yy.get(0)))));
      Function solver = ca.rootfinder("solver", S.name, f, S.options);
      Function refsol = new Function("refsol", List.of(y, x), List.of(vs(ca.sin(x), ca.sqrt(ca.sin(x))).minus(y0)));
      checkfun(solver, refsol, ins(n, 0), 4, 2, S.name);
    }
  }

  @Test
  public void testKINSol1c() {
    message("Scalar KINSol problem, n=0, constraint");
    SX x = SX.sym("x");
    Function f = new Function("f", List.of(x), List.of(ca.sin(x)));
    Function solver = ca.rootfinder("solver", "kinsol", f, opts("constraints", new long[] {-1}));
    DM out = solver.call(ins(-6)).get(0);
    almost(out.to_double(), -2 * Math.PI, 5, "");
  }

  @Test
  public void test_kinsol_user_defined_scalar() {
    Assumptions.assumeTrue(hasRootfinder("kinsol"));
    SX x = SX.sym("x");
    Function f = new Function("f", List.of(x), List.of(x.power(2).minus(2)));
    for (String strategy : List.of("none", "linesearch")) {
      Function solver = ca.rootfinder("solver", "kinsol", f,
          opts("linear_solver_type", "user_defined", "strategy", strategy, "abstol", 1e-10));
      for (double x0 : new double[] {1, 2}) {
        chk((DM) solver.invoke(x0), Math.sqrt(2), 8);
        assertTrue((Boolean) solver.stats().get("success"));
      }
    }
  }

  @Test
  public void test_kinsol_user_defined_vector() {
    Assumptions.assumeTrue(hasRootfinder("kinsol"));
    for (int symbolic = 0; symbolic < 2; symbolic++) {
      Function f;
      if (symbolic == 0) {
        SX x = SX.sym("x", 2), p = SX.sym("p", 2);
        SX g = vs(x.get(0).power(2).plus(p.get(0).times(x.get(1))).minus(p.get(1)), x.get(0).minus(x.get(1)));
        f = new Function("f", List.of(x, p), List.of(g));
      } else {
        MX x = MX.sym("x", 2), p = MX.sym("p", 2);
        MX g = ca.vertcat(List.of(x.get(0).power(2).plus(p.get(0).times(x.get(1))).minus(p.get(1)), x.get(0).minus(x.get(1))));
        f = new Function("f", List.of(x, p), List.of(g));
      }
      for (String strategy : List.of("none", "linesearch")) {
        Map<String, Object> o = opts("strategy", strategy, "abstol", 1e-10, "u_scale", new double[] {0.5, 2},
            "f_scale", new double[] {0.25, 4});
        Function dense = ca.rootfinder("dense", "kinsol", f, with(o, "linear_solver_type", "dense"));
        Function solver = ca.rootfinder("solver", "kinsol", f, with(o, "linear_solver_type", "user_defined"));
        double[][][] cases = {{{1, 6}, {2, 2}}, {{2, 15}, {3, 3}}};
        for (double[][] c : cases) {
          DM parameters = dm(c[0]);
          DM result = (DM) solver.invoke(dm(1, 1), parameters);
          chk(result, dm(c[1]), 8);
          chk(result, (DM) dense.invoke(dm(1, 1), parameters), 8);
          chk((DM) f.invoke(result, parameters), DM.zeros(2), 8);
          assertTrue((Boolean) solver.stats().get("success"));
        }
      }
    }
  }

  @Test
  public void test_constraints() {
    for (Solver S : solvers) {
      if (S.name.contains("kinsol")) continue;
      if (S.name.contains("newton")) continue;
      message(S.name + " " + S.options);
      SX x = SX.sym("x", 2);
      SX g = vs(ca.mtimes(x.plus(3).T(), x.minus(2)), ca.mtimes(x.minus(4).T(), x.plus(dm(1, 2))));
      Function f = new Function("f", List.of(x), List.of(g));
      Function solver = ca.rootfinder("solver", S.name, f, with(S.options, "constraints", new long[] {-1, 0}));
      DM out = solver.call(ins(0)).get(0);
      double r = Math.sqrt(1201);
      chk(out, dm(-3.0 / 50 * (r - 1), 2.0 / 25 * (r - 1)), 6);
      f = new Function("f", List.of(x), List.of(g));
      solver = ca.rootfinder("solver", S.name, f, with(S.options, "constraints", new long[] {1, 0}));
      out = solver.call(ins(0)).get(0);
      chk(out, dm(3.0 / 50 * (r + 1), -2.0 / 25 * (r + 1)), 6);
    }
  }

  @Test
  public void test_implicitbug() {
    MX X0 = MX.sym("X0"), V = MX.sym("V");
    MX Veq = ca.vertcat(List.of(V.get(0).minus(X0)));
    Function vfcn = new Function("vfcn", List.of(V, X0), List.of(Veq), opts("ad_weight", 0, "ad_weight_sp", 1));
    Function vfcnSx = vfcn.expand("vfcn_sx", opts("ad_weight", 0, "ad_weight_sp", 1));
    Function ifcn = ca.rootfinder("ifcn", "newton", vfcnSx, opts("linear_solver", "csparse"));
    MX Vs = ifcn.call(List.of(new MX(0.0), X0), true).get(0);
    double f = 1;
    Function F = new Function("F", List.of(X0), List.of(X0.times(f).plus(Vs)), opts("ad_weight", 0, "ad_weight_sp", 1));
    Function J = jacobianOld(F, 0, 0);
    List<DM> Jout = J.call(ins(1));
    message(Jout.get(0).toString());
    message(J.toString());
    chk(Jout.get(0), dm(2), 10);
  }

  @Test
  public void test_extra_outputs() {
    SX x = SX.sym("x"), a = SX.sym("a");
    Function f = new Function("f", List.of(x, a), List.of(ca.tan(x).minus(a), ca.sqrt(a).times(x.power(2))));
    for (Solver S : solvers) {
      message(S.name);
      Function solver = ca.rootfinder("solver", S.name, f, with(S.options, "ad_weight_sp", 1));
      List<DM> solverIn = ins(0.1, 0.3);
      Function refsol = new Function("refsol", List.of(x, a),
          List.of(ca.atan(a), ca.sqrt(a).times(ca.atan(a).power(2))));
      checkfun(solver, refsol, solverIn, 5);
      if (S.codegen) checkCodegen(solver, solverIn);
    }
  }

  @Test
  public void test_no_success() {
    SX x = SX.sym("x"), y = SX.sym("y");
    for (Solver S : solvers) {
      Map<String, SX> rfp = new LinkedHashMap<>();
      rfp.put("x", vs(x, y));
      rfp.put("g", vs(ca.sin(x).minus(2), ca.sin(y).minus(2)));
      Function solver = ca.rootfinder("solver", S.name, rfp, with(S.options, "error_on_fail", false));
      solver.call(Map.of("x0", new DM(0.0)));
      assertFalse((Boolean) solver.stats().get("success"));
      Function solver2 = ca.rootfinder("solver", S.name, rfp, S.options);
      assertInException("process", () -> solver2.call(Map.of("x0", new DM(0.0))));
    }
  }

  @Test
  public void test_loop() {
    SX x = SX.sym("x");
    for (Solver S : solvers) {
      Map<String, SX> rfp = new LinkedHashMap<>();
      rfp.put("x", x);
      rfp.put("g", x.power(3).minus(x.times(2)).plus(2));
      Function solver = ca.rootfinder("solver", S.name, rfp, S.options);
      if (S.name.equals("kinsol")) continue;
      if (S.name.equals("nlpsol")) continue;
      if (S.name.equals("fast_newton")) continue;
      message(S.name);
      message(solver.toString());
      DM res = solver.call(Map.of("x0", new DM(0.0))).get("x");
      chk(res, -1.7692923542386, 10);
    }
  }

  @Test
  public void test_bisection() {
    SX x = SX.sym("x"), p = SX.sym("p");
    Function f = new Function("f", List.of(x, p), List.of(x.power(2).minus(p)));
    Function solver = ca.rootfinder("solver", "bisection", f);
    Function refsol = new Function("refsol", List.of(x, p), List.of(ca.sqrt(p)));
    checkfun(solver, refsol, ins(1, 2), 7);
    checkSerialize(solver, ins(1, 2));
    Map<String, Object> stats = solver.stats();
    assertTrue((Boolean) stats.get("success"));
    assertEquals("converged_abstol", stats.get("return_status"));
    assertTrue(num(stats.get("search_iter")) >= 1);
    assertTrue(num(stats.get("iter_count")) >= 1);

    chk((DM) solver.invoke(-1, 2), -Math.sqrt(2), 7);
    refsol = new Function("refsol", List.of(x, p), List.of(ca.sqrt(p).unaryMinus()));
    checkfun(solver, refsol, ins(-1, 2), 7);

    solver = ca.rootfinder("solver", "bisection", f, opts("lb", -10, "ub", 0));
    chk((DM) solver.invoke(1, 2), -Math.sqrt(2), 7);
    solver = ca.rootfinder("solver", "bisection", f, opts("lb", 0, "ub", 10));
    chk((DM) solver.invoke(-1, 2), Math.sqrt(2), 7);

    solver = ca.rootfinder("solver", "bisection", f, opts("abstol", 1e-12, "abstol_step", 1e-12));
    chk((DM) solver.invoke(1, 2), Math.sqrt(2), 10);

    Map<String, SX> rfp = new LinkedHashMap<>();
    rfp.put("x", vs(x, p));
    rfp.put("g", vs(x, p));
    assertInException("scalar", () -> ca.rootfinder("solver", "bisection", rfp));

    Function f2 = new Function("f", List.of(x, p), List.of(x.power(2).plus(p)));
    solver = ca.rootfinder("solver", "bisection", f2, opts("error_on_fail", false, "lb", -1, "ub", 1));
    solver.invoke(1, 1);
    assertFalse((Boolean) solver.stats().get("success"));
    assertEquals("failed_to_bracket_root", solver.stats().get("return_status"));
    Function solver2 = ca.rootfinder("solver", "bisection", f2, opts("lb", -1, "ub", 1));
    assertInException("process", () -> solver2.invoke(1, 1));
  }

  @Test
  public void test_segfault_codegen() {
    MX x = MX.sym("x"), y = MX.sym("y");
    double[] u = new double[10], v = new double[10];
    for (int i = 0; i < 10; i++) {
      u[i] = 1 + 4.0 * i / 9;
      v[i] = u[i] * u[i];
    }
    Function interp = ca.interpolant("interp", "bspline", new double[][] {u}, v, opts());
    MX fx = (MX) interp.invoke(x);
    Function f = new Function("f", List.of(x), List.of(ca.if_else(ca.lt(x, u[0]), x.power(2),
        ca.if_else(ca.gt(x, u[9]), x.power(2), fx, false), false)));
    Function res = new Function("res", List.of(x, y), List.of(y.minus((MX) f.invoke(x))));
    Function rf = ca.rootfinder("rf", "fast_newton", res);
    Function g = new Function("g", List.of(x), List.of(x.minus((MX) rf.invoke(x, x))));
    checkCodegen(g, ins(1.01));
  }

  @Test
  public void test_kinsol_return_status() {
    Assumptions.assumeTrue(hasRootfinder("kinsol"));
    SX x = SX.sym("x");
    Map<String, SX> rfp = new LinkedHashMap<>();
    rfp.put("x", x);
    rfp.put("g", x.power(2).minus(2));
    Function G = ca.rootfinder("G", "kinsol", rfp);
    chk(G.call(Map.of("x0", new DM(1.0))).get("x"), new DM(Math.sqrt(2)), 6);
    assertEquals("KIN_SUCCESS", G.stats().get("return_status"));
  }

  private static Map<String, DM> args(Object... kv) {
    Map<String, DM> m = new LinkedHashMap<>();
    for (int i = 0; i < kv.length; i += 2)
      m.put((String) kv[i], kv[i + 1] instanceof DM ? (DM) kv[i + 1] : new DM(((Number) kv[i + 1]).doubleValue()));
    return m;
  }

  @Test
  public void test_generic() {
    SX a = SX.sym("a"), b = SX.sym("b"), c = SX.sym("c");
    Function f = new Function("f", List.of(a, b, c), List.of(a.power(2).minus(b), b.power(2).minus(c), a.times(b).times(c).plus(1)),
        List.of("a", "b", "c"), List.of("f", "g", "h"));
    double a0 = 1, b0 = 2, c0 = 3;
    Function rf = ca.rootfinder("rf", "newton", f, opts("implicit_input", 0, "implicit_output", 0));
    Map<String, DM> res = rf.call(args("a0", a0, "b", b0, "c", c0));
    almost(Math.pow(res.get("a").to_double(), 2) - b0, 0, 10, "");
    almost(b0 * b0 - c0 - res.get("g").to_double(), 0, 10, "");
    almost(res.get("a").to_double() * b0 * c0 + 1 - res.get("h").to_double(), 0, 10, "");
    message(rf.toString());
    Function jac = rf.get_function("jac_g_x");
    assertEquals(List.of("a", "b", "c"), jac.name_in());
    assertEquals(List.of("jac_f_a", "f", "g", "h"), jac.name_out());

    rf = ca.rootfinder("rf", "newton", f, opts("implicit_input", 1, "implicit_output", 0));
    res = rf.call(args("a", a0, "b0", b0, "c", c0));
    almost(a0 * a0 - res.get("b").to_double(), 0, 10, "");
    almost(Math.pow(res.get("b").to_double(), 2) - c0 - res.get("g").to_double(), 0, 10, "");
    almost(a0 * res.get("b").to_double() * c0 + 1 - res.get("h").to_double(), 0, 10, "");
    jac = rf.get_function("jac_g_x");
    assertEquals(List.of("a", "b", "c"), jac.name_in());
    assertEquals(List.of("jac_f_b", "f", "g", "h"), jac.name_out());

    rf = ca.rootfinder("rf", "newton", f, opts("implicit_input", 2, "implicit_output", 1));
    res = rf.call(args("a", a0, "b", b0, "c0", c0));
    almost(a0 * a0 - b0 - res.get("f").to_double(), 0, 10, "");
    almost(b0 * b0 - res.get("c").to_double(), 0, 10, "");
    almost(a0 * b0 * res.get("c").to_double() + 1 - res.get("h").to_double(), 0, 10, "");
    jac = rf.get_function("jac_g_x");
    assertEquals(List.of("a", "b", "c"), jac.name_in());
    assertEquals(List.of("jac_g_c", "f", "g", "h"), jac.name_out());

    Map<String, SX> rfp = new LinkedHashMap<>();
    rfp.put("x", a);
    rfp.put("p", vs(b, c));
    rfp.put("g", a.power(2).minus(b));
    rf = ca.rootfinder("rf", "newton", rfp);
    res = rf.call(args("x0", a0, "p", dm(b0, c0)));
    almost(Math.pow(res.get("x").to_double(), 2) - b0, 0, 10, "");
    jac = rf.get_function("jac_g_x");
    assertEquals(List.of("x", "p"), jac.name_in());
    assertEquals(List.of("jac_g_x", "g"), jac.name_out());
  }
}

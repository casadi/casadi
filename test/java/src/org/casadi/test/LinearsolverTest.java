package org.casadi.test;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.io.File;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.casadi.CodeGenerator;
import org.casadi.DM;
import org.casadi.Function;
import org.casadi.GlobalOptions;
import org.casadi.Linsol;
import org.casadi.MX;
import org.casadi.SX;
import org.casadi.Sparsity;
import org.casadi.ca;
import org.junit.jupiter.api.Assumptions;
import org.junit.jupiter.api.Test;

/** Port of test/python/linearsolver.py */
public class LinearsolverTest extends CasadiTestCase {

  private static final class LSolver {
    final String name;
    final Map<String, Object> options;
    final Set<String> req;

    LSolver(String name, Set<String> req) {
      this.name = name;
      this.options = new HashMap<>();
      this.req = req;
    }
  }

  private static final List<LSolver> lsolvers = new ArrayList<>();

  static {
    // csparsecholesky is commented out in the Python original
    String[][] cands = {{"csparse"}, {"lapacklu"}, {"lapackqr"}, {"symbolicqr"}, {"lsqr"}, {"ma27", "symmetry"},
      {"mumps", "symmetry"}, {"qr"}, {"ldl", "posdef", "symmetry"}};
    for (String[] c : cands) {
      try {
        ca.load_linsol(c[0]);
        lsolvers.add(new LSolver(c[0], Set.of(java.util.Arrays.copyOfRange(c, 1, c.length))));
      } catch (RuntimeException e) {
        // plugin not available
      }
    }
    StringBuilder s = new StringBuilder("linear solvers");
    for (LSolver l : lsolvers) s.append(" ").append(l.name);
    System.out.println(s);
  }

  private static Function nullspacewrapper(String name, Sparsity sp, Map<String, Object> options) {
    SX a = SX.sym("a", sp);
    return new Function(name, List.of(a), List.of(ca.nullspace(a)), options);
  }

  private static Function jacobianOld(Function f, long i, long j) {
    List<String> out = new ArrayList<>();
    out.add("jac:" + f.name_out(j) + ":" + f.name_in(i));
    out.addAll(f.name_out());
    return f.factory(f.name() + "_jac", f.name_in(), out);
  }

  private DM uniform(int n, int m) {
    double[][] a = new double[n][m];
    for (int i = 0; i < n; i++)
      for (int j = 0; j < m; j++) a[i][j] = random.nextDouble();
    return new DM(a);
  }

  private static DM ev(Function f, DM... args) {
    return (DM) f.invoke((Object[]) args);
  }

  private static List<DM> inputList(Function f, Map<String, DM> inputs) {
    List<DM> r = new ArrayList<>();
    for (String n : f.name_in()) r.add(inputs.get(n));
    return r;
  }

  private static void checkSerialize(Function F, List<DM> inputs) {
    Function F2 = Function.deserialize(F.serialize(Map.of("debug", true)));
    List<DM> a = F.call(inputs), b = F2.call(inputs);
    for (int i = 0; i < F.n_out(); i++) checkarray(a.get(i), b.get(i), 1e-16, "serialize output " + i);
  }

  private static void checkSerialize(Function F, Map<String, DM> inputs) {
    Function F2 = Function.deserialize(F.serialize(Map.of("debug", true)));
    Map<String, DM> a = F.call(inputs), b = F2.call(inputs);
    assertEquals(F.name_out(), F2.name_out());
    for (String k : F.name_out()) checkarray(a.get(k), b.get(k), 1e-16, "serialize output " + k);
  }

  // as helpers.py check_codegen: compile generated C into a shared library, compare via external
  private static void checkCodegen(Function F, List<DM> inputs) {
    if (!RUN_SLOW) return;
    try {
      Path dir = Files.createTempDirectory("casadi_cg");
      String name = "codegen_" + Long.toHexString(System.nanoTime());
      CodeGenerator cg = new CodeGenerator(name + ".c");
      cg.add(F);
      cg.generate(dir.toString() + File.separator);
      String so = dir.resolve(name + ".so").toString();
      String inc = GlobalOptions.getCasadiIncludePath();
      List<String> cmd = new ArrayList<>(List.of("gcc", "-std=c89", "-fPIC", "-shared", "-O0"));
      if (!inc.isEmpty()) cmd.add("-I" + inc);
      cmd.addAll(List.of(dir.resolve(name + ".c").toString(), "-o", so, "-lm"));
      Process p = new ProcessBuilder(cmd).inheritIO().start();
      assertEquals(0, p.waitFor(), "gcc failed");
      Function G = ca.external(F.name(), so);
      List<DM> a = F.call(inputs), b = G.call(inputs);
      for (int i = 0; i < F.n_out(); i++) checkarray(a.get(i), b.get(i), 1e-13, "codegen output " + i);
      dir.toFile().deleteOnExit();
      for (File e : dir.toFile().listFiles()) e.deleteOnExit();
    } catch (java.io.IOException | InterruptedException e) {
      throw new RuntimeException(e);
    }
  }

  private static void checkCodegen(Function F, Map<String, DM> inputs) {
    checkCodegen(F, inputList(F, inputs));
  }

  @Test
  public void test_nullspace() {
    DM[] As = {
      dm(new double[][] {{1, 1.3}, {2, 5}, {1, 0.5}, {1.8, 1.7}}),
      dm(new double[][] {{1, 1.3}, {2, 5}, {1, 0.5}}),
      dm(new double[][] {{1, 1.3}, {2, 5}, {1, 0.5}, {0.2, 0.3}, {-0.3, 0.7}}),
      dm(new double[][] {{1, 0}, {0, 0}, {0, 1}, {0, 0}}),
      dm(new double[][] {{1.3, 0, 0.4, 1}, {0.2, 0.1, 11, 0}, {0, 1, 0, 0}, {0.7, 0.9, 0, 0}, {1.1, 0.99, 0, 0}})};
    Map<String, Object> options = new HashMap<>();
    for (DM A : As) {
      long n = A.size1(), m = A.size2();
      Function solver = nullspacewrapper("solver", A.T().sparsity(), options);
      List<DM> solver_out = solver.call(List.of(A.T()));

      checkarray(ca.mtimes(A.T(), solver_out.get(0)), DM.zeros(m, n - m));
      checkarray(ca.mtimes(solver_out.get(0).T(), solver_out.get(0)), DM.eye(n - m));

      options.put("ad_weight", 0);
      options.put("ad_weight_sp", 0);
      solver = nullspacewrapper("solver", A.T().sparsity(), options);
      Function Jf = jacobianOld(solver, 0, 0);

      options.put("ad_weight", 1);
      options.put("ad_weight_sp", 1);
      solver = nullspacewrapper("solver", A.T().sparsity(), options);
      Function Jb = jacobianOld(solver, 0, 0);

      List<DM> Jf_out = Jf.call(List.of(A.T()));
      List<DM> Jb_out = Jb.call(List.of(A.T()));
      checkarray(Jf_out.get(0), Jb_out.get(0));
      checkarray(Jf_out.get(1), Jb_out.get(1));

      Function d = solver.forward(1);
      DM r = uniform((int) n, (int) m);
      List<DM> d_out = d.call(List.of(A.T(), new DM(5), r.T()));
      DM exact = d_out.get(0);

      DM nom = ev(solver, A.T());
      double eps = 1e-6;
      DM pert = ev(solver, A.plus(r.times(eps)).T());
      DM fd = pert.minus(nom).div(eps);

      DM V = uniform((int) (n - m), (int) (n - m));
      V = V.plus(V.T());
      System.out.println(V);
      DM a = ca.mtimes(List.of(nom, V, fd.T())).plus(ca.mtimes(List.of(fd, V, nom.T())));
      DM b = ca.mtimes(List.of(nom, V, exact.T())).plus(ca.mtimes(List.of(exact, V, nom.T())));
      System.out.println("here: " + a.minus(b));

      V = DM.eye(n);
      a = ca.mtimes(List.of(nom.T(), V, fd)).plus(ca.mtimes(List.of(fd.T(), V, nom)));
      b = ca.mtimes(List.of(nom.T(), V, exact)).plus(ca.mtimes(List.of(exact.T(), V, nom)));
      checkarray(a, b, 1e-5, "");
    }
  }

  @Test
  public void test_simple_solve() {
    DM A_ = dm(new double[][] {{3, 7}, {1, 2}});
    DM b_ = dm(1, 0.5);
    MX A = MX.sym("A", A_.sparsity()), b = MX.sym("b", b_.sparsity());
    for (LSolver s : lsolvers) {
      MX C = ca.solve(A, b, s.name, s.options);
      DM A0 = s.req.contains("symmetry") ? A_.T().plus(A_) : A_;
      Function f = new Function("f", List.of(A, b), List.of(C));
      DM f_out = ev(f, A0, b_);
      checkarray(f_out, ca.solve(A0, b_));
      checkarray(ca.mtimes(A0, f_out), b_);
    }
  }

  @Test
  public void test_pseudo_inverse() {
    random.setSeed(0);
    DM A_ = uniform(4, 6);
    MX A = MX.sym("A", A_.sparsity());
    SX As = SX.sym("A", A_.sparsity());
    for (LSolver s : lsolvers) {
      if (s.req.contains("symmetry")) continue;
      System.out.println(s.name);
      MX B = ca.pinv(A, s.name, s.options);
      Function f = new Function("f", List.of(A), List.of(B));
      checkarray(ca.mtimes(A_, ev(f, A_)), DM.eye(4));
      f = new Function("f", List.of(As), List.of(ca.pinv(As)));
      checkarray(ca.mtimes(A_, ev(f, A_)), DM.eye(4));
      ca.solve(ca.mtimes(A, A.T()), A, s.name, s.options);
      ca.pinv(A_, s.name, s.options);
    }

    A_ = uniform(3, 5);
    A = MX.sym("A", A_.sparsity());
    As = SX.sym("A", A_.sparsity());
    for (LSolver s : lsolvers) {
      System.out.println(s.name);
      if (s.req.contains("symmetry")) continue;
      MX B = ca.pinv(A, s.name, s.options);
      Function f = new Function("f", List.of(A), List.of(B));
      checkarray(ca.mtimes(A_, ev(f, A_)), DM.eye(3));
      f = new Function("f", List.of(As), List.of(ca.pinv(As)));
      checkarray(ca.mtimes(A_, ev(f, A_)), DM.eye(3));
    }
  }

  @Test
  public void test_simple_solve_dmatrix() {
    DM A = dm(new double[][] {{3, 7}, {1, 2}});
    DM b = dm(1, 0.5);
    for (LSolver s : lsolvers) {
      System.out.println(s.name);
      DM A0 = s.req.contains("symmetry") ? A.T().plus(A) : A;
      DM C = ca.solve(A0, b, s.name, s.options);
      DM sol = ca.solve(A0, b);
      checkarray(C, sol);
      checkarray(ca.mtimes(A0, sol), b);
    }
  }

  @Test
  public void test_simple_trans() {
    DM A = dm(new double[][] {{3, 1}, {7, 2}});
    for (LSolver s : lsolvers) {
      DM A0 = s.req.contains("symmetry") ? A.T().plus(A) : A;
      Linsol solver = new Linsol("solver", s.name, A0.sparsity(), s.options);
      assertInException("No stats available since Linsol did not solve a problem yet.", solver::stats);
      DM b = dm(1, 0.5);
      DM x = solver.solve(A0.T(), b);
      DM res = ca.solve(A0.T(), b);
      System.out.println(solver.stats());
      checkarray(x, res);
    }
  }

  @Test
  public void test_simple() {
    DM A = dm(new double[][] {{3, 1}, {7, 2}});
    for (LSolver s : lsolvers) {
      System.out.println(s.name);
      DM A0 = s.req.contains("symmetry") ? A.T().plus(A) : A;
      Linsol solver = new Linsol("solver", s.name, A0.sparsity(), s.options);
      DM b = dm(1, 0.5);
      DM x = solver.solve(A0, b);
      checkarray(x, ca.solve(A0, b));
    }
  }

  @Test
  public void test_simple_function_indirect() {
    for (LSolver s : lsolvers) {
      DM A_ = dm(new double[][] {{3, 1}, {7, 2}});
      MX A = MX.sym("A", A_.sparsity());
      DM b_ = dm(1, 0.5);
      MX b = MX.sym("b", b_.sparsity());
      if (s.req.contains("symmetry")) A_ = A_.T().plus(A_);
      System.out.println(s.name);
      Linsol solver = new Linsol("solver", s.name, A.sparsity(), s.options);
      Map<String, DM> solver_in = new LinkedHashMap<>();
      solver_in.put("A", A_);
      solver_in.put("B", b_);

      Map<String, MX> sol = new LinkedHashMap<>();
      sol.put("X", solver.solve(A, b));
      sol.put("A", A);
      sol.put("B", b);
      Function relay = new Function("relay", sol, List.of("A", "B"), List.of("X"));

      MX A_0 = A.get(0, 0), A_1 = A.get(0, 1), A_2 = A.get(1, 0), A_3 = A.get(1, 1);
      MX b_0 = b.get(0), b_1 = b.get(1);
      MX det = A_0.times(A_3).minus(A_2.times(A_1));
      MX X = ca.vertcat(List.of(
        A_3.div(det).times(b_0).plus(A_1.unaryMinus().div(det).times(b_1)),
        A_2.unaryMinus().div(det).times(b_0).plus(A_0.div(det).times(b_1))));
      Map<String, MX> sd = new LinkedHashMap<>();
      sd.put("A", A);
      sd.put("B", b);
      sd.put("X", X);
      Function solution = new Function("solution", sd, List.of("A", "B"), List.of("X"));

      checkfunction(relay, solution, inputList(relay, solver_in));
      checkSerialize(relay, solver_in);

      if (s.name.equals("qr") || s.name.equals("ldl")) checkCodegen(relay, solver_in);
    }
  }

  @Test
  public void test_simple_solve_node() {
    for (LSolver s : lsolvers) {
      DM[][] cases = {
        {dm(new double[][] {{3, 1}, {7, 2}}), dm(new double[][] {{1, 0.3}, {0.5, 0.7}})},
        {ca.sparsify(dm(new double[][] {{3, 0}, {7, 2}})), dm(new double[][] {{1, 0.3}, {0.5, 0.7}})},
        {dm(new double[][] {{3, 1}, {7, 2}}), ca.sparsify(dm(new double[][] {{1, 0}, {0, 0.7}}))}};
      for (DM[] c : cases) {
        DM A_ = c[0], b_ = c[1];
        if (s.req.contains("symmetry")) A_ = A_.T().plus(A_);
        MX A = MX.sym("A", A_.sparsity()), b = MX.sym("b", b_.sparsity());
        System.out.println(s.name);
        Linsol solver = new Linsol("solver", s.name, A.sparsity(), s.options);
        for (boolean tr : new boolean[] {true, false}) {
          MX x = solver.solve(A, b, tr);
          Function f = new Function("f", List.of(A, b), List.of(x));
          DM f_out = ev(f, A_, b_);
          DM x_ = solver.solve(A_, b_, tr);
          checkarray(x_, f_out);

          MX A_0, A_1, A_2, A_3;
          if (tr) {
            A_0 = A.get(0, 0);
            A_1 = A.get(1, 0);
            A_2 = A.get(0, 1);
            A_3 = A.get(1, 1);
          } else {
            A_0 = A.get(0, 0);
            A_1 = A.get(0, 1);
            A_2 = A.get(1, 0);
            A_3 = A.get(1, 1);
          }
          MX b_0 = b.get(0, 0), b_1 = b.get(1, 0), c_0 = b.get(0, 1), c_1 = b.get(1, 1);
          MX det = A_0.times(A_3).minus(A_2.times(A_1));
          MX p = A_3.div(det), q = A_1.unaryMinus().div(det), u = A_2.unaryMinus().div(det), v = A_0.div(det);
          Function solution = new Function("solution", List.of(A, b), List.of(ca.blockcat(List.of(
            List.of(p.times(b_0).plus(q.times(b_1)), p.times(c_0).plus(q.times(c_1))),
            List.of(u.times(b_0).plus(v.times(b_1)), u.times(c_0).plus(v.times(c_1)))))));

          checkfunction(f, solution, List.of(A_, b_));

          if (!s.name.contains("SymbolicQR")) continue;
          Function solversx = f.expand("expand_" + f.name());
          checkfunction(solversx, solution, randomInputs(solversx), 1e-7);
        }
      }
    }
  }

  @Test
  public void test_large_sparse() {
    int n = 10;
    for (LSolver s : lsolvers) {
      random.setSeed(3); // numpy seed 1 analogue: ldl (no pivoting) hits a zero pivot for java.util.Random seeds 1,2,4,5
      DM A = randDM(n, n, 0.5);
      if (s.req.contains("symmetry")) A = A.T().plus(A);
      DM b = randDM(n, 3, 0.5);
      MX As = MX.sym("A", A.sparsity()), bs = MX.sym("B", b.sparsity());
      DM C = ca.solve(A, b, s.name, s.options);
      checkarray(ca.mtimes(A, C), b, s.name);
      Function f = new Function("f", List.of(As, bs), List.of(ca.solve(As, bs, s.name, s.options)));
      checkarray(ca.mtimes(A, ev(f, A, b)), b, s.name);
    }
  }

  @Test
  public void test_ma27() {
    double n = Double.NaN;
    DM A = dm(new double[][] {{2, 3, n, n, n}, {3, 0, 4, n, 6}, {n, 4, 1, 5, n}, {n, n, 5, 0, n}, {n, 6, n, n, 1}});
    A = A.get(ca.sparsify(ca.eq(A, A)).sparsity());
    DM b = dm(8, 45, 31, 15, 17);
    DM ref = ca.solve(ca.densify(A), b);
    for (LSolver s : lsolvers) {
      if (s.req.contains("posdef")) continue;
      MX.sym("A", A.sparsity());
      MX.sym("B", b.sparsity());
      DM C = ca.solve(A, b, s.name, s.options);
      checkarray(ref, C, s.name);
    }
  }

  @Test
  public void test_large_sparse2() {
    random.setSeed(1);
    int n = 10;
    for (LSolver s : lsolvers) {
      System.out.println(s.name);
      DM A = randDM(n, n, 0.5);
      DM b = randDM(n, 3);
      if (s.req.contains("posdef")) {
        A = ca.mtimes(A.T(), A);
        A = ca.densify(A);
      } else if (s.req.contains("symmetry")) {
        A = A.T().plus(A);
        A.set(Sparsity.diag(n), new DM(1e-8));
        A = ca.densify(A);
      }
      MX As = MX.sym("A", A.sparsity()), bs = MX.sym("B", b.sparsity());
      DM C = ca.solve(A, b, s.name, s.options);
      double tol = s.name.contains("ma") ? 1e-7 : 1e-10;
      checkarray(ca.mtimes(A, C), b, tol, s.name);

      List<MX> sym = List.of(As, ca.densify(As), ca.densify(As).T(), ca.densify(As.T()), As.T());
      List<DM> num = List.of(A, ca.densify(A), ca.densify(A).T(), ca.densify(A.T()), A.T());
      for (int k = 0; k < sym.size(); k++) {
        Function f = new Function("f", List.of(As, bs), List.of(ca.solve(sym.get(k), bs, s.name, s.options)));
        checkarray(ca.mtimes(num.get(k), ev(f, A, b)), b, tol, s.name + " case " + k);
      }
    }
  }

  @Test
  public void test_dimmismatch() {
    DM A = DM.eye(5);
    DM b = DM.ones(4, 1);
    MX As = MX.sym("A", A.sparsity()), bs = MX.sym("b", b.sparsity());
    for (LSolver s : lsolvers) {
      assertThrows(Exception.class, () -> ca.solve(A, b, s.name, s.options));
      assertThrows(Exception.class, () -> ca.solve(As, bs, s.name, s.options));
    }
  }

  @Test
  public void test_cache() {
    int n = 5;
    List<DM> As = new ArrayList<>();
    for (int i = 0; i < 10; i++) As.add(uniform(n, n));

    MX A = MX.sym("A", n, n);
    Function f = new Function("f", List.of(A), List.of(ca.solve(A, DM.ones(n), "qr", Map.of("cache", 2))));

    DM S1 = ev(f, As.get(0));
    checkarray(S1, ev(f, As.get(0)));

    DM S2 = ev(f, As.get(1));
    checkarray(S2, ev(f, As.get(1)));
    checkarray(S1, ev(f, As.get(0)));

    DM S3 = ev(f, As.get(2));
    checkarray(S3, ev(f, As.get(2)));
    checkarray(S1, ev(f, As.get(0)));
    checkarray(S2, ev(f, As.get(1)));

    checkCodegen(f, List.of(As.get(0)));
    checkSerialize(f, List.of(As.get(0)));
  }

  @Test
  public void test_thread_safety() {
    MX x = MX.sym("x"), y = MX.sym("y");
    for (LSolver s : lsolvers) {
      if (s.name.contains("mumps")) continue;
      System.out.println(s.name);
      Function f = new Function("f", List.of(x, y), List.of(x.times(x).minus(y)));
      Function finv = ca.rootfinder("finv", "newton", f, Map.of("linear_solver", s.name, "linear_solver_options", s.options));

      Function finv_par = finv.map(200, "thread", 4);
      DM lin = ca.linspace(0.0, 10.0, 200);
      DM res = ev(finv_par, DM.ones(1, 200), lin.T());
      checkarray(ca.norm_inf(res.T().minus(ca.sqrt(lin))), new DM(0), 1e-5, s.name);

      MX z = ca.solve(ca.blockcat(List.of(List.of(x.plus(1.0), new MX(0)), List.of(new MX(0), y.plus(1.0)))),
        ca.vertcat(List.of(x, y)), s.name, s.options);
      f = new Function("f", List.of(x, y), List.of(z));
      Function f_par = f.map(200, "thread", 4);
      f_par.invoke(ca.linspace(10.0, 0.0, 200).T(), lin.T());
    }
  }

  @Test
  public void test_issue2664() {
    DM bnum = DM.rand(3, 3);
    MX x = MX.sym("x", 3, 3);
    for (MX b : List.of(MX.sym("b", 3, 3), MX.sym("b", Sparsity.lower(3)), MX.sym("b", Sparsity.upper(3)))) {
      Function f = new Function("f", List.of(x, b), List.of(ca.solve(x, b, "symbolicqr")));
      for (DM bn : List.of(bnum, bnum.get(Sparsity.lower(3)), bnum.get(Sparsity.upper(3)))) {
        SX d = ca.diagcat(List.of(new SX(1), new SX(2), new SX(3)));
        DM A = ca.evalf(f.call(List.of(d, new SX(bn))).get(0));
        DM B = ev(f, ca.diagcat(List.of(new DM(1), new DM(2), new DM(3))), bn);
        checkarray(A, B);
      }
    }
  }

  @Test
  public void test_issue3489() {
    SX E_f = SX.sym("E_f"), E_r = SX.sym("E_r"), T = SX.sym("T"), k_f = SX.sym("k_f"), k_r = SX.sym("k_r");
    SX k0_f = SX.sym("k0_f"), k0_r = SX.sym("k0_r");
    SX x_1 = SX.sym("x_1"), x_2 = SX.sym("x_2"), x_3 = SX.sym("x_3"), x_4 = SX.sym("x_4");
    SX r_1 = SX.sym("r_1"), r_2 = SX.sym("r_2"), r_3 = SX.sym("r_3"), r_4 = SX.sym("r_4");
    SX y_1 = SX.sym("y_1"), y_2 = SX.sym("y_2"), y_3 = SX.sym("y_3"), y_4 = SX.sym("y_4");

    SX eq1 = k_f.minus(k0_f).minus(E_f).plus(T);
    SX eq2 = k_r.minus(k0_r).minus(E_r).plus(T);
    java.util.function.UnaryOperator<SX> rx = r -> r.minus(k_f).plus(x_1.unaryMinus()).plus(x_2.unaryMinus()).minus(k_r).plus(x_3).plus(x_4);
    SX eq3 = rx.apply(r_1);
    SX eq4 = rx.apply(r_2);
    SX eq5 = rx.apply(r_3);
    SX eq6 = rx.apply(r_4);
    SX eq7 = ca.plus(1 * 0.2, x_1.times(-1.0)).plus(r_1);
    SX eq8 = ca.plus(1 * 0.8, x_2.times(-1.0)).plus(r_2);
    SX eq9 = x_3.times(-1.0).plus(r_3);
    SX eq10 = x_4.times(-1.0).plus(r_4);
    SX eq11 = x_1.minus(y_1), eq12 = x_2.minus(y_2), eq13 = x_3.minus(y_3), eq14 = x_4.minus(y_4);

    List<SX> all_vars = new ArrayList<>(List.of(k_f, k_r, x_1, x_2, x_3, x_4, r_1, r_2, r_3, r_4));
    all_vars.addAll(List.of(y_1, y_2, y_3, y_4));
    List<SX> all_eq = new ArrayList<>(List.of(eq11, eq12, eq13, eq14));
    all_eq.addAll(List.of(eq3, eq4, eq5, eq6, eq7, eq8, eq9, eq10));
    all_eq.addAll(List.of(eq1, eq2));

    SX x = ca.vertcat(all_vars);
    SX p = ca.vertcat(List.of(E_f, E_r, k0_f, k0_r, T));
    SX g = ca.vertcat(all_eq);

    DM A = ca.evalf(ca.jacobian(g, x));
    DM B = ca.evalf(ca.jacobian(g, p));
    MX.sym("A", A.sparsity());
    MX.sym("B", B.sparsity());

    Map<String, SX> rf = Map.of("x", x, "p", p, "g", g);
    Function sim = ca.rootfinder("s", "fast_newton", rf, Map.of("ad_weight_sp", 0));
    Sparsity S1 = sim.jac_sparsity(0, 1);
    sim = ca.rootfinder("s", "fast_newton", rf, Map.of("ad_weight_sp", 1));
    Sparsity S2 = sim.jac_sparsity(0, 1);
    assertTrue(S1.is_equal(S2));

    random.setSeed(1);
    int m = 14, np = 13;
    for (int i = 0; i < 25; i++) {
      DM Ar;
      while (true) {
        Ar = randDM(m, m, 0.15);
        if (ca.sprank(Ar) == m) break;
      }
      DM Br = randDM(m, np, 0.1);
      Sparsity C = ca.sparsify(ca.solve(Ar, Br), 1e-7).sparsity();
      MX P = MX.sym("p", np);
      Function f0 = new Function("f", List.of(P), List.of(ca.solve(Ar, ca.mtimes(Br, P))), Map.of("ad_weight_sp", 0));
      S1 = f0.jac_sparsity(0, 0);
      Function f1 = new Function("f", List.of(P), List.of(ca.solve(Ar, ca.mtimes(Br, P))), Map.of("ad_weight_sp", 1));
      S2 = f1.jac_sparsity(0, 0);
      assertTrue(S1.is_equal(S2));
      assertTrue(C.is_subset(S1));
    }
  }

  @Test
  public void test_empty_system() {
    Linsol L = new Linsol("L", "qr", new Sparsity(0, 0));
    DM r = L.solve(new DM(0, 0), new DM(0, 2));
    assertEquals(List.of(0L, 2L), List.of(r.size1(), r.size2()));
  }

  @Test
  public void test_issue2734() {
    DM example_input = ca.blockcat(List.of(List.of(new DM(8), new DM(1), new DM(6)), List.of(new DM(3), new DM(5), new DM(7)),
      List.of(new DM(4), new DM(9), new DM(2))));

    SX a = SX.sym("a", 3, 3);
    Function f = new Function("f", List.of(a), List.of(ca.pinv(a)));
    DM res1 = ev(f, example_input);

    MX am = MX.sym("am", 3, 3);
    Function fm = new Function("f", List.of(am), List.of(ca.pinv(am, "symbolicqr")));
    DM res2 = ev(fm, example_input);

    am = MX.sym("am", 3, 3);
    fm = new Function("f", List.of(am), List.of(ca.pinv(am, "symbolicqr")));
    fm = fm.expand();
    DM res3 = ev(fm, example_input);

    checkarray(res1, res2);
    checkarray(res1, res3);
  }
}

package org.casadi.test;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.function.UnaryOperator;

import org.casadi.CodeGenerator;
import org.casadi.DM;
import org.casadi.Function;
import org.casadi.GlobalOptions;
import org.casadi.MX;
import org.casadi.SX;
import org.casadi.Slice;
import org.casadi.Sparsity;
import org.casadi.ca;
import org.junit.jupiter.api.AfterAll;
import org.junit.jupiter.api.Assumptions;
import org.junit.jupiter.api.Test;

/** Port of test/python/mx.py, second part (from test_project_simplify) */
public class MX2Test extends CasadiTestCase {

  private static final double[][] A15 = {
    {-26.9091, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    {0, -26.9091, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    {0, 0, -26.9091, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    {0, 0, 0, -14.1393, -0.348109, -0.382033, 3.6491, 9.40422, -5.05449, 3.03347, -0.00126949, -0.000414104, 0, 0.00498232, 0.00495617},
    {0, 0, 0, -0.348109, -13.6315, -0.194202, -6.41868, -0.0189287, 3.34725, 1.86987, -0.000645326, -0.000210504, 0, 0.0025327, 0.0025194},
    {0, 0, 0, -0.382033, -0.194202, -13.6677, 10.1342, -13.0101, -0.890078, -8.99621, -0.000708215, -0.000231018, 0, 0.00277951, 0.00276493},
    {0, 0, 0, 0, 0, 0, -27.0487, 0.36775, -0.277186, 0.564995, 0.424389, 0.0235157, 0.227774, 0, 0},
    {0, 0, 0, 0, 0, 0, 0.113776, -27.3241, 0.159254, -0.520988, -0.0235157, 0.424389, 0.132128, 0, 0},
    {0, 0, 0, 0, 0, 0, 0.227472, -0.0735535, -26.9135, 0.206825, -0.227774, -0.132128, 0.424389, 0, 0},
    {0, 0, 0, 0, 0, 0, 0.332188, -1.02565, -0.0471482, -28.3499, 0.132128, -0.227774, 0.0235157, 0, 0},
    {0, 0, 0, -0.00126949, -0.000645326, -0.000708215, 0.00763923, 0.0016513, -0.0044564, -0.00284703, -0.119342, -0.000941301, -0.00160961, -1.62792e-05, 0.00156863},
    {0, 0, 0, -0.000414104, -0.000210504, -0.000231018, 0.0024919, 0.000538651, -0.00145367, -0.000928698, 0.000941301, -0.119493, 0.000742542, -0.00155303, -7.50989e-06},
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2.1684e-19, -2.1684e-19, -0.171654, -0.000868032, 0.000868032},
    {0, 0, 0, 0.00181259, 0.000921407, 0.0010112, -0.0109074, -0.00235775, 0.00636292, 0.00406504, 0, -0.000567731, -0.000868032, -0.00087529, 0},
    {0, 0, 0, 0.00166876, 0.000848291, 0.000930959, -0.0100419, -0.00217066, 0.005858, 0.00374247, 0.00052268, 0, -0.000868032, 0, -0.000874062}
  };

  private static int codegenCounter = 0;
  private static final List<File> tmpDirs = new ArrayList<>();

  @AfterAll
  static void removeCodegenDirs() {
    for (File d : tmpDirs) rmtree(d);
  }

  private static void checkCodegen(Function f, List<DM> inputs) {
    checkCodegen(f, inputs, Map.of(), 1e-15, "c89", false);
  }

  private static String readFile(String p) {
    try {
      return java.nio.file.Files.readString(java.nio.file.Path.of(p)).replace("\r\n", "\n");
    } catch (IOException e) {
      throw new RuntimeException(e);
    }
  }

  /** As check_codegen in helpers.py: generate, compile with gcc, compare (only with RUN_SLOW); returns the main-mode output file */
  private static String checkCodegen(Function f, List<DM> inputs, Map<String, ?> opts, double tol, String std, boolean main) {
    return checkCodegen(f, inputs, opts, tol, std, main, false);
  }

  private static String lastLib;
  private static Function lastExternal;

  /** extras: with_jac_sparsity, with_forward and with_reverse of check_codegen */
  private static String checkCodegen(Function f, List<DM> inputs, Map<String, ?> opts, double tol, String std, boolean main, boolean extras) {
    if (!RUN_SLOW) return null;
    try {
      File tmp = java.nio.file.Files.createTempDirectory("casadi_java_codegen").toFile();
      tmpDirs.add(tmp);
      String dir = tmp.toString() + "/";
      String name = "codegen_java_" + (codegenCounter++);
      Map<String, Object> o = new LinkedHashMap<>(opts);
      if (main) o.put("main", true);
      CodeGenerator cg = new CodeGenerator(name, o);
      cg.add(f, extras);
      if (extras) {
        cg.add(f.reverse(1), true);
        cg.add(f.forward(1), true);
      }
      cg.generate(dir);
      String libdir = GlobalOptions.getCasadiPath(), inc = GlobalOptions.getCasadiIncludePath();
      for (boolean shared : main ? new boolean[] {true, false} : new boolean[] {true}) {
        List<String> cmd = new ArrayList<>(List.of("gcc", "-pedantic", "-std=" + std, "-fPIC", "-Wall", "-Werror", "-Wextra",
            "-I" + inc, "-Wno-unknown-pragmas", "-Wno-long-long", "-Wno-unused-parameter", "-O3", dir + name + ".c", "-o",
            dir + name + (shared ? ".so" : ""), "-L" + libdir, "-Wl,-rpath," + libdir));
        if (shared) cmd.add(4, "-shared");
        else cmd.add("-lm");
        assertEquals(0, new ProcessBuilder(cmd).inheritIO().start().waitFor(), "compile " + String.join(" ", cmd));
      }
      String out = null;
      if (main) {
        f.generate_in(dir + f.name() + "_in.txt", inputs);
        out = dir + f.name() + "_out.txt";
        lastMainOut = out;
        ProcessBuilder pb = new ProcessBuilder(dir + name, f.name());
        pb.redirectInput(new File(dir + f.name() + "_in.txt"));
        pb.redirectOutput(new File(out));
        assertEquals(0, pb.start().waitFor());
        List<DM> outputs = f.generate_out(out);
        List<DM> r = f.call(inputs);
        for (int i = 0; i < r.size(); i++) checkarray(r.get(i), outputs.get(i), tol, "main output " + i);
      }
      Function f2 = ca.external(f.name(), dir + name + ".so");
      lastLib = dir + name + ".so";
      lastExternal = f2;
      List<DM> r = f.call(inputs), r2 = f2.call(inputs);
      for (int i = 0; i < r.size(); i++) checkarray(r.get(i), r2.get(i), tol, "codegen output " + i);
      return out;
    } catch (IOException | InterruptedException e) {
      throw new RuntimeException(e);
    }
  }

  private static String lastMainOut;

  /** check_codegen(main=True) when generate_out is expected to fail (the program prints extra lines) */
  private static void checkCodegenMainParsingError(Function f, List<DM> inputs) {
    lastMainOut = null;
    if (!RUN_SLOW) return;
    try {
      checkCodegen(f, inputs, Map.of(), 1e-15, "c89", true);
    } catch (RuntimeException e) {
      assertTrue(String.valueOf(e.getMessage()).contains("Parsing error"), e.getMessage());
      return;
    }
    org.junit.jupiter.api.Assertions.fail("expected Parsing error");
  }

  /** As check_serialize in helpers.py */
  private static void checkSerialize(Function f, List<DM> inputs) {
    Function f2 = Function.deserialize(f.serialize(Map.of("debug", true)));
    List<DM> r = f.call(inputs), r2 = f2.call(inputs);
    for (int i = 0; i < r.size(); i++) checkarray(r.get(i), r2.get(i), 1e-16, "serialize output " + i);
  }

  private static Function hessianOld(Function f, long i, long j) {
    List<String> out = new ArrayList<>();
    out.add("hess:" + f.name_out(j) + ":" + f.name_in(i) + ":" + f.name_in(i));
    out.add("grad:" + f.name_out(j) + ":" + f.name_in(i));
    out.addAll(f.name_out());
    return f.factory(f.name() + "_hess", f.name_in(), out);
  }

  private static DM range(int a, int b) {
    double[] r = new double[b - a];
    for (int i = 0; i < r.length; i++) r[i] = a + i;
    return new DM(r);
  }

  @SafeVarargs
  private static <T> List<T> cat(List<T>... l) {
    List<T> r = new ArrayList<>();
    for (List<T> e : l) r.addAll(e);
    return r;
  }

  @Test
  public void test_project_simplify() {
    MX x = MX.sym("x", Sparsity.lower(3));
    DM.rng(1);
    DM X0 = DM.rand(x.sparsity());

    MX y = ca.project(x, Sparsity.upper(3));
    MX z = ca.project(y, ca.sparsify(dm(new double[][] {{0, 0, 0}, {1, 1, 1}, {1, 1, 1}})).sparsity());

    Function f = new Function("f", List.of(x), List.of(z));
    Function fs = f.transform();
    assertTrue(fs.n_nodes() > 3);

    checkfunction(f, fs, dms(X0));

    z = ca.project(y, ca.sparsify(dm(new double[][] {{1, 0, 0}, {0, 0, 0}, {0, 0, 1}})).sparsity());

    f = new Function("f", List.of(x), List.of(z));
    fs = f.transform();

    assertTrue(fs.n_nodes() != 0);

    fs.disp(true);

    checkfunction(f, fs, dms(X0));
  }

  @Test
  public void test_repmat() {
    DM a = dm(new double[][] {{1, 2}, {3, 4}, {5, 6}});
    checkarray(ca.repmat(a, 2, 3), ca.kron(DM.ones(2, 3), a));
  }

  @Test
  public void test_transpose_bug() {
    DM A = ca.sparsify(dm(A15));

    MX As = MX.sym("As", A.sparsity());

    Function f = new Function("f", List.of(As), List.of(ca.densify(As.T()), ca.densify(As).T(), As.T(), As, ca.densify(As)));

    List<DM> f_out = f.call(dms(A));

    checkarray(f_out.get(0), A.T());
    checkarray(f_out.get(1), A.T());
    checkarray(f_out.get(2), A.T());
    checkarray(f_out.get(3), A);
    checkarray(f_out.get(4), A);
  }

  @Test
  public void test_bizarre_bug() {
    Assumptions.assumeTrue(hasLinsol("csparse"));
    DM A = ca.sparsify(dm(A15));

    DM b = range(0, 15);

    MX Bs = MX.sym("Bs", b.sparsity());
    MX As = MX.sym("As", A.sparsity());

    MX Ast = As.T();

    Function r = new Function("r", List.of(As, Bs), List.of(ca.solve(Ast, Bs, "csparse")));
    Function R = new Function("R", List.of(As, Bs), List.of(ca.solve(ca.densify(Ast), Bs, "csparse")));

    DM r_out = r.call(dms(A, b)).get(0);
    DM R_out = R.call(dms(A, b)).get(0);

    r.sparsity_out(0).spy();
    R.sparsity_out(0).spy();

    checkarray(ca.mtimes(A.T(), R_out), b, 1e-9, "");
    checkarray(r_out, R_out);
  }

  @Test
  public void test_depends_on() {
    MX a = MX.sym("a");
    MX b = MX.sym("b");

    assertTrue(ca.depends_on(ca.power(a, 2), a));
    assertTrue(ca.depends_on(a, a));
    assertFalse(ca.depends_on(new MX(0), a));

    MX ab = ca.vertcat(List.of(a, b));
    assertTrue(ca.depends_on(ca.power(a, 2), ab));
    assertTrue(ca.depends_on(a, ab));
    assertFalse(ca.depends_on(new MX(0), ab));
    assertTrue(ca.depends_on(ca.power(b, 2), ab));
    assertTrue(ca.depends_on(b, ab));
    assertTrue(ca.depends_on(ca.power(a, 2).plus(ca.power(b, 2)), ab));
    assertTrue(ca.depends_on(a.plus(b), ab));
    MX zero = new MX(0);
    assertTrue(ca.depends_on(ca.vertcat(List.of(zero, a)), a));
    assertTrue(ca.depends_on(ca.vertcat(List.of(a, zero)), a));
    assertFalse(ca.depends_on(ca.vertcat(List.of(zero, b)), a));
    assertFalse(ca.depends_on(ca.vertcat(List.of(b, zero)), a));
    assertTrue(ca.depends_on(ca.vertcat(List.of(ca.power(a, 2), ca.power(b, 2))), ab));
    assertTrue(ca.depends_on(ca.vertcat(List.of(a, zero)), ab));
    assertTrue(ca.depends_on(ca.vertcat(List.of(zero, b)), ab));
    assertTrue(ca.depends_on(ca.vertcat(List.of(b, zero)), ab));
    assertFalse(ca.depends_on(ca.vertcat(List.of(zero, zero)), ab));
  }

  @Test
  public void test_vertcat_simp() {
    MX x = MX.sym("x", 10);
    MX y = MX.sym("y");
    MX z = MX.sym("z");
    DM x_ = range(0, 10);
    DM y_ = dm(20);
    DM z_ = dm(30);

    java.util.function.Function<List<MX>, DM> evalvertcat = a -> {
      Function f = new Function("f", List.of(x, y, z), List.of(ca.vertcat(a)));
      return f.call(dms(x_, y_, z_)).get(0);
    };
    List<MX> xs = ca.vertsplit(x);
    List<MX> ly = List.of(y), lz = List.of(z);
    List<MX> xm = xs.subList(0, 9), xp = xs.subList(1, 10);
    DM x_m = x_.get(new Slice(0, -1)), x_p = x_.get(new Slice(1, 10));

    checkarray(evalvertcat.apply(xs), x_);
    checkarray(evalvertcat.apply(cat(xs, ly)), ca.vertcat(List.of(x_, y_)));
    checkarray(evalvertcat.apply(cat(lz, xs, ly, xs, lz)), ca.vertcat(List.of(z_, x_, y_, x_, z_)));
    checkarray(evalvertcat.apply(xm), x_m);
    checkarray(evalvertcat.apply(cat(xm, ly)), ca.vertcat(List.of(x_m, y_)));
    checkarray(evalvertcat.apply(cat(lz, xm, ly, xm, lz)), ca.vertcat(List.of(z_, x_m, y_, x_m, z_)));
    checkarray(evalvertcat.apply(xp), x_p);
    checkarray(evalvertcat.apply(cat(xp, ly)), ca.vertcat(List.of(x_p, y_)));
    checkarray(evalvertcat.apply(cat(lz, xp, ly, xp, lz)), ca.vertcat(List.of(z_, x_p, y_, x_p, z_)));
    List<MX> g = cat(xs.subList(5, 10), xs.subList(0, 5));
    DM x5 = x_.get(new Slice(5, 10)), x0 = x_.get(new Slice(0, 5));
    checkarray(evalvertcat.apply(g), ca.vertcat(List.of(x5, x0)));
    checkarray(evalvertcat.apply(cat(g, ly)), ca.vertcat(List.of(x5, x0, y_)));
    checkarray(evalvertcat.apply(cat(lz, g, ly, g, lz)), ca.vertcat(List.of(z_, x5, x0, y_, x5, x0, z_)));

    List<MX> w = ca.vertsplit(x, 2);
    List<MX> r = new ArrayList<>();
    for (MX i : w) r.addAll(ca.vertsplit(i));

    checkarray(evalvertcat.apply(r), x_);

    r = new ArrayList<>();
    for (MX i : w) r.addAll(cat(ca.vertsplit(i), ly));
    System.out.println("vertcat: " + r);
    System.out.println("result: " + ca.vertcat(r));

    r = new ArrayList<>();
    for (MX i : w) r.addAll(ca.vertsplit(i));
    System.out.println("vertcat: " + r);
    System.out.println("result: " + ca.vertcat(cat(r, ly)));

    assertTrue(ca.is_equal(ca.vertcat(ca.vertsplit(x)), x));
  }

  @Test
  public void test_horzcat_simp() {
    MX x = MX.sym("x", 1, 10);
    MX y = MX.sym("y");
    MX z = MX.sym("z");
    DM x_ = range(0, 10).T();
    DM y_ = dm(20);
    DM z_ = dm(30);

    java.util.function.Function<List<MX>, DM> evalhorzcat = a -> {
      Function f = new Function("f", List.of(x, y, z), List.of(ca.horzcat(a)));
      return (DM) f.invoke(x_, y_, z_);
    };
    List<MX> xs = ca.horzsplit(x);
    List<MX> ly = List.of(y), lz = List.of(z);
    List<MX> xm = xs.subList(0, 9), xp = xs.subList(1, 10);
    DM x_m = x_.get(0, new Slice(0, -1)), x_p = x_.get(0, new Slice(1, 10));

    checkarray(evalhorzcat.apply(xs), x_);
    checkarray(evalhorzcat.apply(cat(xs, ly)), ca.horzcat(List.of(x_, y_)));
    checkarray(evalhorzcat.apply(cat(lz, xs, ly, xs, lz)), ca.horzcat(List.of(z_, x_, y_, x_, z_)));
    checkarray(evalhorzcat.apply(xm), x_m);
    checkarray(evalhorzcat.apply(cat(xm, ly)), ca.horzcat(List.of(x_m, y_)));
    checkarray(evalhorzcat.apply(cat(lz, xm, ly, xm, lz)), ca.horzcat(List.of(z_, x_m, y_, x_m, z_)));
    checkarray(evalhorzcat.apply(xp), x_p);
    checkarray(evalhorzcat.apply(cat(xp, ly)), ca.horzcat(List.of(x_p, y_)));
    checkarray(evalhorzcat.apply(cat(lz, xp, ly, xp, lz)), ca.horzcat(List.of(z_, x_p, y_, x_p, z_)));
    List<MX> g = cat(xs.subList(5, 10), xs.subList(0, 5));
    DM x5 = x_.get(0, new Slice(5, 10)), x0 = x_.get(0, new Slice(0, 5));
    checkarray(evalhorzcat.apply(g), ca.horzcat(List.of(x5, x0)));
    checkarray(evalhorzcat.apply(cat(g, ly)), ca.horzcat(List.of(x5, x0, y_)));
    checkarray(evalhorzcat.apply(cat(lz, g, ly, g, lz)), ca.horzcat(List.of(z_, x5, x0, y_, x5, x0, z_)));

    List<MX> w = ca.horzsplit(x, 2);
    List<MX> r = new ArrayList<>();
    for (MX i : w) r.addAll(ca.horzsplit(i));

    checkarray(evalhorzcat.apply(r), x_);

    r = new ArrayList<>();
    for (MX i : w) r.addAll(cat(ca.horzsplit(i), ly));
    System.out.println("vertcat: " + r);
    System.out.println("result: " + ca.horzcat(r));

    r = new ArrayList<>();
    for (MX i : w) r.addAll(ca.horzsplit(i));
    System.out.println("vertcat: " + r);
    System.out.println("result: " + ca.horzcat(cat(r, ly)));

    assertTrue(ca.is_equal(ca.horzcat(ca.horzsplit(x)), x));
  }

  @Test
  public void test_vertsplit_simp() {
    List<MX> dvars = new ArrayList<>();
    for (int i = 0; i < 5; i++) dvars.add(MX.sym("abcdefghijklm".substring(i, i + 1)));
    double[] dvars_ = {0, 1, 2, 3, 4};

    MX zz = MX.sym("zz", 2);
    DM zz_ = dm(11, 12);
    MX y = MX.sym("y");
    MX z = MX.sym("z");
    DM y_ = dm(20);
    DM z_ = dm(30);

    MX aa = MX.sym("aa", 5);
    DM aa_ = range(100, 105);

    java.util.function.BiFunction<MX, Long, List<DM>> evalvertsplit = (a, incr) -> {
      List<MX> sp = incr == null ? ca.vertsplit(a) : ca.vertsplit(a, incr);
      System.out.println(sp);
      Function f = new Function("f", cat(dvars, List.of(y, z, zz, aa)), sp);
      List<DM> f_in = new ArrayList<>();
      for (int i = 0; i < 5; i++) f_in.add(dm(dvars_[i]));
      f_in.addAll(List.of(y_, z_, zz_, aa_));
      return f.call(f_in);
    };

    List<DM> s = evalvertsplit.apply(ca.vertcat(cat(List.of(y), dvars, List.of(z))), null);
    checkarray(s.get(0), y_);
    for (int i = 0; i < 5; i++) checkarray(s.get(1 + i), dvars_[i]);
    checkarray(s.get(6), z_);

    s = evalvertsplit.apply(ca.vertcat(cat(List.of(y), dvars, List.of(z))), 2L);

    checkarray(s.get(0), ca.vertcat(List.of(y_, dm(dvars_[0]))));
    checkarray(s.get(1), dm(dvars_[1], dvars_[2]));
    checkarray(s.get(2), dm(dvars_[3], dvars_[4]));
    checkarray(s.get(3), z_);

    s = evalvertsplit.apply(ca.vertcat(List.of(y, zz, z, zz)), 2L);

    checkarray(s.get(0), ca.vertcat(List.of(y_, zz_.get(0))));
    checkarray(s.get(1), ca.vertcat(List.of(zz_.get(1), z_)));
    checkarray(s.get(2), zz_);

    s = evalvertsplit.apply(ca.vertcat(List.of(y, zz, z, zz)), 3L);

    checkarray(s.get(0), ca.vertcat(List.of(y_, zz_.get(0), zz_.get(1))));
    checkarray(s.get(1), ca.vertcat(List.of(z_, zz_.get(0), zz_.get(1))));

    s = evalvertsplit.apply(ca.vertcat(List.of(zz, zz)), 2L);
    checkarray(s.get(0), zz_);
    checkarray(s.get(1), zz_);

    s = evalvertsplit.apply(ca.vertcat(cat(List.of(zz), dvars)), null);
    checkarray(s.get(0), zz_.get(0));
    checkarray(s.get(1), zz_.get(1));

    for (int i = 0; i < 5; i++) checkarray(s.get(2 + i), dvars_[i]);

    s = evalvertsplit.apply(ca.vertcat(cat(dvars, List.of(aa))), 5L);
    checkarray(s.get(0), dm(dvars_));
    checkarray(s.get(1), aa_);

    s = evalvertsplit.apply(ca.vertcat(cat(dvars, List.of(aa))), 4L);
    checkarray(s.get(0), dm(0, 1, 2, 3));
    checkarray(s.get(1), dm(4, 100, 101, 102));
    checkarray(s.get(2), dm(103, 104));

    s = evalvertsplit.apply(ca.vertcat(cat(dvars, List.of(aa))), 6L);
    checkarray(s.get(0), dm(0, 1, 2, 3, 4, 100));
    checkarray(s.get(1), dm(101, 102, 103, 104));

    for (int i = 0; i < 5; i++) assertTrue(ca.is_equal(ca.vertsplit(ca.vertcat(dvars)).get(i), dvars.get(i)));
  }

  @Test
  public void test_horzsplit_simp() {
    List<MX> dvars = new ArrayList<>();
    for (int i = 0; i < 5; i++) dvars.add(MX.sym("abcdefghijklm".substring(i, i + 1)));
    double[] dvars_ = {0, 1, 2, 3, 4};

    MX zz = MX.sym("zz", 1, 2);
    DM zz_ = dm(11, 12).T();
    MX y = MX.sym("y");
    MX z = MX.sym("z");
    DM y_ = dm(20);
    DM z_ = dm(30);

    MX aa = MX.sym("aa", 1, 5);
    DM aa_ = range(100, 105);

    java.util.function.BiFunction<MX, Long, List<DM>> evalhorzsplit = (a, incr) -> {
      List<MX> sp = incr == null ? ca.horzsplit(a) : ca.horzsplit(a, incr);
      System.out.println(sp);
      Function f = new Function("f", cat(dvars, List.of(y, z, zz, aa)), sp);
      List<DM> f_in = new ArrayList<>();
      for (int i = 0; i < 5; i++) f_in.add(dm(dvars_[i]));
      f_in.addAll(List.of(y_, z_, zz_, aa_));
      return f.call(f_in);
    };

    List<DM> s = evalhorzsplit.apply(ca.horzcat(cat(List.of(y), dvars, List.of(z))), null);
    checkarray(s.get(0), y_);
    for (int i = 0; i < 5; i++) checkarray(s.get(1 + i), dvars_[i]);
    checkarray(s.get(6), z_);

    s = evalhorzsplit.apply(ca.horzcat(cat(List.of(y), dvars, List.of(z))), 2L);

    checkarray(s.get(0), ca.vertcat(List.of(y_, dm(dvars_[0]))).T());
    checkarray(s.get(1), dm(dvars_[1], dvars_[2]).T());
    checkarray(s.get(2), dm(dvars_[3], dvars_[4]).T());
    checkarray(s.get(3), z_.T());

    s = evalhorzsplit.apply(ca.horzcat(List.of(y, zz, z, zz)), 2L);

    checkarray(s.get(0), ca.vertcat(List.of(y_, zz_.get(0, 0))).T());
    checkarray(s.get(1), ca.vertcat(List.of(zz_.get(0, 1), z_)).T());
    checkarray(s.get(2), zz_);

    s = evalhorzsplit.apply(ca.horzcat(List.of(y, zz, z, zz)), 3L);

    checkarray(s.get(0), ca.vertcat(List.of(y_, zz_.get(0, 0), zz_.get(0, 1))).T());
    checkarray(s.get(1), ca.vertcat(List.of(z_, zz_.get(0, 0), zz_.get(0, 1))).T());

    s = evalhorzsplit.apply(ca.horzcat(List.of(zz, zz)), 2L);
    checkarray(s.get(0), zz_);
    checkarray(s.get(1), zz_);

    s = evalhorzsplit.apply(ca.horzcat(cat(List.of(zz), dvars)), null);
    checkarray(s.get(0), zz_.get(0, 0));
    checkarray(s.get(1), zz_.get(0, 1));

    for (int i = 0; i < 5; i++) checkarray(s.get(2 + i), dvars_[i]);

    s = evalhorzsplit.apply(ca.horzcat(cat(dvars, List.of(aa))), 5L);
    checkarray(s.get(0), dm(dvars_).T());
    checkarray(s.get(1), aa_.T());

    s = evalhorzsplit.apply(ca.horzcat(cat(dvars, List.of(aa))), 4L);
    checkarray(s.get(0), dm(0, 1, 2, 3).T());
    checkarray(s.get(1), dm(4, 100, 101, 102).T());
    checkarray(s.get(2), dm(103, 104).T());

    s = evalhorzsplit.apply(ca.horzcat(cat(dvars, List.of(aa))), 6L);
    checkarray(s.get(0), dm(0, 1, 2, 3, 4, 100).T());
    checkarray(s.get(1), dm(101, 102, 103, 104).T());

    for (int i = 0; i < 5; i++) assertTrue(ca.is_equal(ca.horzsplit(ca.horzcat(dvars)).get(i), dvars.get(i)));
  }

  @Test
  public void test_vertsplit_derivative() {
    MX m = MX.sym("X", 10);

    Function f = new Function("f", List.of(m), List.of(ca.vertsplit(m).get(0)));

    f.reverse(1);
  }

  @Test
  public void test_MX_const_sp() {
    MX x = MX.sym("x", 4, 1);

    Sparsity sp = Sparsity.triplet(3, 3, new long[] {0, 1, 2, 2}, new long[] {0, 0, 1, 2});

    Function f = new Function("f", List.of(x), List.of(x.get_nz(false, new DM(sp, range(0, (int) sp.nnz())))));

    Function g = new Function("g", List.of(x), List.of(new MX(sp, x)));

    List<DM> f_in = dms(range(1, 5));

    checkfunction(f, g, f_in);
    checkCodegen(f, f_in);
    checkCodegen(g, f_in);
  }

  @Test
  public void test_reshape_sp() {
    MX x = MX.sym("x", 4, 1);

    Function f = new Function("f", List.of(x), List.of(ca.reshape(x, new long[] {2, 2})));

    SX sx = SX.sym("x", 4, 1);

    Function g = new Function("g", List.of(sx), List.of(ca.reshape(sx, new long[] {2, 2})));

    List<DM> f_in = dms(range(1, 5));

    checkfunction(f, g, f_in);
    checkCodegen(f, f_in);
    checkCodegen(g, f_in);
  }

  @Test
  public void test_issue1041() {
    MX x = MX.sym("x", 2);

    MX y = ca.vertsplit(x, new long[] {0, 1, 2}).get(1);

    Function f = new Function("f", List.of(x), List.of(y));

    hessianOld(f, 0, 0);
  }

  @Test
  public void test_bug_1042() {
    MX x = MX.sym("x", 2, 1);

    Function mf = new Function("mf", List.of(x), List.of(x.times(x.get(0, 0))));

    Function mfunction = mf.expand("expand_" + mf.name());

    Function mfg = mf.reverse(1);

    Function mfunctiong = mfunction.reverse(1);

    List<DM> f_in = dms(dm(0), dm(5), dm(1, 2));

    checkfunction(mfg, mfunctiong, f_in);
  }

  @Test
  public void test_bug_1042bis() {
    MX x = MX.sym("x", 2, 1);
    MX a = MX.sym("ax", 2, 1);
    MX i1 = x.get(0, 0);
    MX z = i1.times(x);
    MX i3 = i1.times(a);
    i3 = ca.dot(x, a);
    Function d = new Function("d", List.of(x, a), List.of(z, i3));

    Function dx = d.expand("expand_" + d.name());
    List<DM> dx_in = dms(dm(1, 2), dm(3, 4));

    checkfunction(d, dx, dx_in);
  }

  @Test
  public void test_bug_1042tris() {
    MX x = MX.sym("x", 2, 1);
    MX a = MX.sym("ax", 2, 1);
    Function d = new Function("d", List.of(x, a), List.of(ca.dot(x, a)));

    Function dx = d.expand("expand_" + d.name());
    List<DM> dx_in = dms(dm(1, 2), dm(3, 4));

    checkfunction(d, dx, dx_in);
  }

  @Test
  public void test_bug_1046() {
    MX x = MX.sym("x", 1, 1);
    MX y = MX.sym("y", 1, 1);
    MX z = ca.jacobian(x, y);

    assertEquals(0, z.nnz());
  }

  @Test
  public void test_expand_free() {
    MX x = MX.sym("x");
    MX y = MX.sym("y");
    assertInException("free", () -> new Function("f", List.of(x), List.of(x.plus(y))));

    Function f = new Function("f", List.of(x), List.of(x.plus(y)), Map.of("allow_free", true));
    Function g = new Function("g", List.of(x, y), List.of((MX) f.invoke(x)));
    g.expand();
    checkfunction(g, g, dms(dm(1), dm(2)));
  }

  private interface Kind<T> {
    T zeros(long n, long m);
    List<T> vertsplit(T x, long[] o);
    List<T> horzsplit(T x, long[] o);
    long size1(T x);
    long size2(T x);
  }

  private static <T extends org.casadi.GenericMatrix> void singularcat(Kind<T> c) {
    long[] o = {0, 5, 10};
    T x0 = c.zeros(10, 0);
    List<T> x1s = c.vertsplit(x0, o);
    for (T x : x1s) shape(c, x, 5, 0);
    shape(c, ca.vertcat(x1s), 10, 0);
    shape(c, ca.vertcat(cat(List.of(c.zeros(0, 0)), x1s, List.of(c.zeros(0, 0)))), 10, 0);
    x0 = c.zeros(0, 1);
    shape(c, ca.vertcat(List.of(x0, c.zeros(0, 0), x0)), 0, 1);

    x0 = c.zeros(0, 10);
    x1s = c.horzsplit(x0, o);
    for (T x : x1s) shape(c, x, 0, 5);
    shape(c, ca.horzcat(x1s), 0, 10);
    shape(c, ca.horzcat(cat(List.of(c.zeros(0, 0)), x1s, List.of(c.zeros(0, 0)))), 0, 10);
    x0 = c.zeros(1, 0);
    shape(c, ca.horzcat(List.of(x0, c.zeros(0, 0), x0)), 1, 0);

    x1s = c.vertsplit(c.zeros(10, 0), o);
    List<T> x1st = c.horzsplit(c.zeros(0, 10), o);
    List<T> e = List.of(c.zeros(0, 0));
    shape(c, ca.diagcat(x1s), 10, 0);
    shape(c, ca.diagcat(cat(e, x1s, e)), 10, 0);
    shape(c, ca.diagcat(x1st), 0, 10);
    shape(c, ca.diagcat(cat(e, x1st, e)), 0, 10);
    shape(c, ca.diagcat(cat(x1s, x1st)), 10, 10);
    shape(c, ca.diagcat(cat(e, x1s, x1st, e)), 10, 10);
  }

  private static <T> void shape(Kind<T> c, T x, long n, long m) {
    assertEquals(n, c.size1(x));
    assertEquals(m, c.size2(x));
  }

  @Test
  public void test_singularcat() {
    singularcat(new Kind<MX>() {
      public MX zeros(long n, long m) { return MX.zeros(n, m); }
      public List<MX> vertsplit(MX x, long[] o) { return ca.vertsplit(x, o); }
      public List<MX> horzsplit(MX x, long[] o) { return ca.horzsplit(x, o); }
      public long size1(MX x) { return x.size1(); }
      public long size2(MX x) { return x.size2(); }
    });
    singularcat(new Kind<SX>() {
      public SX zeros(long n, long m) { return SX.zeros(n, m); }
      public List<SX> vertsplit(SX x, long[] o) { return ca.vertsplit(x, o); }
      public List<SX> horzsplit(SX x, long[] o) { return ca.horzsplit(x, o); }
      public long size1(SX x) { return x.size1(); }
      public long size2(SX x) { return x.size2(); }
    });
    singularcat(new Kind<DM>() {
      public DM zeros(long n, long m) { return DM.zeros(n, m); }
      public List<DM> vertsplit(DM x, long[] o) { return ca.vertsplit(x, o); }
      public List<DM> horzsplit(DM x, long[] o) { return ca.horzsplit(x, o); }
      public long size1(DM x) { return x.size1(); }
      public long size2(DM x) { return x.size2(); }
    });
  }

  @Test
  public void test_exprjac() {
    List<UnaryOperator<MX>> funmx = List.of(x -> x.get(0).times(x.get(1)), x -> x.get(0).times(ca.sin(x.get(1))));
    List<UnaryOperator<SX>> funsx = List.of(x -> x.get(0).times(x.get(1)), x -> x.get(0).times(ca.sin(x.get(1))));
    List<java.util.function.BinaryOperator<MX>> opmx = List.of(ca::gradient, ca::jacobian, (f, x) -> (MX) ca.hessian(f, x).get(0));
    List<java.util.function.BinaryOperator<SX>> opsx = List.of(ca::gradient, ca::jacobian, (f, x) -> (SX) ca.hessian(f, x).get(0));
    for (int k = 0; k < 2; k++) {
      for (int o = 0; o < 3; o++) {
        MX x = MX.sym("x", 2);
        System.out.println(funmx.get(k).apply(x));
        System.out.println(opmx.get(o).apply(funmx.get(k).apply(x), x));
        Function f = new Function("f", List.of(x), List.of(opmx.get(o).apply(funmx.get(k).apply(x), x)));

        SX xs = SX.sym("x", 2);
        Function fr = new Function("fr", List.of(xs), List.of(opsx.get(o).apply(funsx.get(k).apply(xs), xs)));

        checkfunction(f, fr, dms(dm(0)));
        checkCodegen(f, dms(dm(0)));
      }
    }
  }

  private static long prod(long[] d) {
    long r = 1;
    for (long e : d) r *= e;
    return r;
  }

  private static Map<Long, Long> dimMap(long[] dim_a, long[] dim_b, long[] ind_a, long[] ind_b) {
    Map<Long, Long> m = new LinkedHashMap<>();
    long[][] dims = {dim_a, dim_b}, inds = {ind_a, ind_b};
    for (int k = 0; k < 2; k++)
      for (int i = 0; i < inds[k].length; i++) {
        long j = inds[k][i];
        if (j < 0) {
          if (m.containsKey(j)) assertEquals((long) m.get(j), dims[k][i]);
          else m.put(j, dims[k][i]);
        }
      }
    return m;
  }

  private static long linIndex(long[] dims, long[] ind) {
    long ret = 0, cumprod = 1;
    for (int i = 0; i < dims.length; i++) {
      ret += ind[i] * cumprod;
      cumprod *= dims[i];
    }
    return ret;
  }

  private static long[] subs(long[] e, List<Long> old, long[] n) {
    long[] r = new long[e.length];
    for (int i = 0; i < e.length; i++) r[i] = old.contains(e[i]) ? n[old.indexOf(e[i])] : e[i];
    return r;
  }

  private static MX myEinstein(MX A, MX B, MX C, long[] dim_a, long[] dim_b, long[] dim_c, long[] ind_a, long[] ind_b, long[] ind_c) {
    MX R = new MX(C);
    Map<Long, Long> d = dimMap(dim_a, dim_b, ind_a, ind_b);
    List<Long> dk = new ArrayList<>(d.keySet());
    long[] val = new long[dk.size()];
    while (true) {
      long ai = linIndex(dim_a, subs(ind_a, dk, val));
      long bi = linIndex(dim_b, subs(ind_b, dk, val));
      long ci = linIndex(dim_c, subs(ind_c, dk, val));
      R.set(R.get(ci).plus(A.get(ai).times(B.get(bi))), false, ci);
      int k = val.length - 1;
      while (k >= 0 && ++val[k] == d.get(dk.get(k))) val[k--] = 0;
      if (k < 0) break;
    }
    return R;
  }

  private void einsteinTests(long[] dim_a, long[] dim_b, long[] dim_c, long[] ind_a, long[] ind_b, long[] ind_c) {
    MX A = MX.sym("A", prod(dim_a));
    MX B = MX.sym("B", prod(dim_b));
    MX C = MX.sym("C", prod(dim_c), 1);

    SX A_sx = SX.sym("A", prod(dim_a));
    SX B_sx = SX.sym("B", prod(dim_b));
    SX C_sx = SX.sym("C", prod(dim_c), 1);

    List<DM> in = dms(randDM(A.size1(), 1), randDM(B.size1(), 1), randDM(C.size1(), 1));

    MX out = ca.einstein(A, B, C, dim_a, dim_b, dim_c, ind_a, ind_b, ind_c);
    SX out_sx = ca.einstein(A_sx, B_sx, C_sx, dim_a, dim_b, dim_c, ind_a, ind_b, ind_c);
    Function f = new Function("f", List.of(A, B, C), List.of(out), Map.of("ad_weight_sp", 0));
    Function frev = new Function("f", List.of(A, B, C), List.of(out), Map.of("ad_weight_sp", 1));
    Function fr = new Function("fr", List.of(A, B, C), List.of(myEinstein(A, B, C, dim_a, dim_b, dim_c, ind_a, ind_b, ind_c)));
    Function f_sx = new Function("f", List.of(A_sx, B_sx, C_sx), List.of(out_sx));

    Function fsx = f.expand();
    checkfunction(f, fr, in);
    checkfunction(fsx, fr, in);
    checkCodegen(f, in);
    checkSerialize(f, in);
    checkfunction(fsx, f_sx, in);

    for (int i = 0; i < 3; i++) {
      assertTrue(f.jac_sparsity(0, i).is_equal(fsx.jac_sparsity(0, i)));
      assertTrue(frev.jac_sparsity(0, i).is_equal(fsx.jac_sparsity(0, i)));
    }
  }

  private static void permutations(List<Long> rest, List<Long> prefix, List<long[]> out) {
    if (rest.isEmpty()) out.add(prefix.stream().mapToLong(Long::longValue).toArray());
    for (int i = 0; i < rest.size(); i++) {
      List<Long> r = new ArrayList<>(rest), p = new ArrayList<>(prefix);
      p.add(r.remove(i));
      permutations(r, p, out);
    }
  }

  @Test
  public void test_einstein() {
    long[] dim_b = {3, 3};
    long[] dim_a = {3, 3};
    long[][] ind = {{0, -1}, {1, -1}, {0, 1}, {-1, -2}, {-2, -1}};

    for (long[] ind_a : ind) {
      for (long[] ind_b : ind) {
        java.util.Set<Long> fa = new java.util.TreeSet<>(), fb = new java.util.TreeSet<>();
        for (long i : ind_a) if (i < 0) fa.add(i);
        for (long i : ind_b) if (i < 0) fb.add(i);
        List<Long> free = new ArrayList<>();
        for (long i : fa) if (!fb.contains(i)) free.add(i);
        for (long i : fb) if (!fa.contains(i)) free.add(i);
        List<long[]> perms = new ArrayList<>();
        permutations(free, new ArrayList<>(), perms);
        for (long[] ind_c : perms) {
          Map<Long, Long> m = dimMap(dim_a, dim_b, ind_a, ind_b);
          long[] dim_c = new long[ind_c.length];
          for (int i = 0; i < ind_c.length; i++) dim_c[i] = m.get(ind_c[i]);
          einsteinTests(dim_a, dim_b, dim_c, ind_a, ind_b, ind_c);
        }
      }
    }

    einsteinTests(new long[] {2, 4, 3}, new long[] {2, 5, 3}, new long[] {5, 4}, new long[] {-1, -2, -3},
        new long[] {-1, -4, -3}, new long[] {-4, -2});
  }

  private static Object timesAny(Object a, Object b) {
    if (a instanceof DM && b instanceof DM) return ((DM) a).times((DM) b);
    MX x = a instanceof DM ? new MX((DM) a) : (MX) a, y = b instanceof DM ? new MX((DM) b) : (MX) b;
    return x.times(y);
  }

  private static long nnzAny(Object a) {
    return a instanceof DM ? ((DM) a).nnz() : ((MX) a).nnz();
  }

  @Test
  public void test_sparsity_operation() {
    List<Object> L = List.of(new MX(new Sparsity(1, 1)), new MX(new Sparsity(2, 1)), MX.sym("x", 1, 1),
        MX.sym("x", new Sparsity(1, 1)), new DM(1), new DM(new Sparsity(1, 1), 1), new DM(new Sparsity(2, 1), 1),
        new DM(Sparsity.dense(2, 1), 1));

    for (Object a : L) {
      for (Object b : L) {
        Object c = timesAny(a, b);

        if (nnzAny(a) == 0 || nnzAny(b) == 0) assertEquals(0, nnzAny(c));
        else assertTrue(nnzAny(c) > 0);
      }
    }
  }

  @Test
  public void test_solve() {
    Assumptions.assumeTrue(hasLinsol("lapackqr"));
    System.out.println(123);
    DM A = DM.rand(3, 3);
    DM B = DM.rand(3, 50);

    DM C = ca.solve(A, B, "lapackqr", Map.of("max_nrhs", 50));
    DM C1 = ca.solve(A, B, "lapackqr", Map.of("max_nrhs", 20));
    DM C2 = ca.solve(A, B, "lapackqr");

    checkarray(C, C1);
    checkarray(C1, C2);
  }

  @Test
  public void test_sparse_lt() {
    MX x = MX.sym("x", Sparsity.lower(5));
    assertEquals(5 * 6 / 2, ca.gt(x, 0).nnz());
    assertEquals(5 * 5, ca.ge(x, 0).nnz());
  }

  @Test
  public void test_inv() {
    DM Av = DM.rand(3, 3);
    {
      SX A = SX.sym("x", 3, 3);
      Function f = new Function("f", List.of(A), List.of(ca.inv(A), new SX(ca.inv(Av)), ca.mpower(A, -1), new SX(ca.mpower(Av, -1))));
      for (DM o : f.call(dms(Av))) checkarray(ca.mtimes(o, Av), DM.eye(3), 1e-9, "");
    }
    {
      MX A = MX.sym("x", 3, 3);
      Function f = new Function("f", List.of(A), List.of(ca.inv(A), new MX(ca.inv(Av)), ca.mpower(A, -1), new MX(ca.mpower(Av, -1))));
      for (DM o : f.call(dms(Av))) checkarray(ca.mtimes(o, Av), DM.eye(3), 1e-9, "");
    }
  }

  /** numpy.interp: piecewise linear, clamped to the end values */
  private static DM npInterp(double[] xq, double[] x, double[] v) {
    double[] r = new double[xq.length];
    for (int k = 0; k < xq.length; k++) {
      double q = xq[k];
      if (q <= x[0]) r[k] = v[0];
      else if (q >= x[x.length - 1]) r[k] = v[v.length - 1];
      else {
        int i = 0;
        while (q > x[i + 1]) i++;
        r[k] = v[i] + (v[i + 1] - v[i]) * (q - x[i]) / (x[i + 1] - x[i]);
      }
    }
    return new DM(r);
  }

  @Test
  public void test_interp1d() {
    double[] v = {7, 3, 4, -3};
    MX vs = MX.sym("x", 4, 1);

    double[] xq = {10, 0, 1, 2, 4, 8, 7, 5, 1.5};
    double[] x = {1, 2, 4, 8};
    Function F = new Function("f", List.of(vs), List.of(ca.interp1d(x, vs, xq)));

    checkarray((DM) F.invoke(v), npInterp(xq, x, v));

    F = new Function("f", List.of(vs), List.of(ca.interp1d(x, vs, xq, "floor")));

    checkarray((DM) F.invoke(v), dm(-3, 7, 7, 3, 4, -3, 4, 4, 7));

    F = new Function("f", List.of(vs), List.of(ca.interp1d(x, vs, xq, "ceil")));

    checkarray((DM) F.invoke(v), dm(-3, 7, 7, 3, 4, -3, -3, -3, 3));

    xq = new double[] {10, 0, 1, 2, 3, 4, 3.5, 3.25, 1.5};
    x = new double[] {1, 2, 3, 4};
    F = new Function("f", List.of(vs), List.of(ca.interp1d(x, vs, xq, "linear", true)));

    checkarray((DM) F.invoke(v), npInterp(xq, x, v));
  }

  @Test
  public void test_bilin_etc() {
    MX x = MX.sym("x", 3, 3);
    MX y = MX.sym("y", 3, 1);
    MX z = MX.sym("z", 3, 1);

    List<DM> in = dms(DM.rand(3, 3), DM.rand(3, 1), DM.rand(3, 1));

    MX[][] es = {{ca.bilin(x, y, z), ca.mtimes(ca.mtimes(y.T(), x), z)},
                 {ca.rank1(x, 0.3, y, z), x.plus(ca.mtimes(y, z.T()).times(0.3))}};
    for (MX[] e : es) {
      Function f = new Function("f", List.of(x, y, z), List.of(e[0]));
      Function fref = new Function("fref", List.of(x, y, z), List.of(e[1]));
      Function f_sx = f.expand();
      checkfunction(f, fref, in);
      checkfunction(f_sx, fref, in);
      checkCodegen(f, in);
    }
  }

  @Test
  public void test_bilin_short() {
    List<DM> in = dms(DM.rand(3, 3), DM.rand(3, 1));
    {
      SX x = SX.sym("x", 3, 3), y = SX.sym("y", 3, 1);
      Function f = new Function("f", List.of(x, y), List.of(ca.bilin(x, y)));
      Function fref = new Function("f", List.of(x, y), List.of(ca.bilin(x, y, y)));
      checkfunction(f, fref, in);
      checkCodegen(f, in);
      checkCodegen(fref, in);
    }
    {
      MX x = MX.sym("x", 3, 3), y = MX.sym("y", 3, 1);
      Function f = new Function("f", List.of(x, y), List.of(ca.bilin(x, y)));
      Function fref = new Function("f", List.of(x, y), List.of(ca.bilin(x, y, y)));
      checkfunction(f, fref, in);
      checkCodegen(f, in);
      checkCodegen(fref, in);
    }
  }

  @Test
  public void test_det_shape() {
    MX X = MX.sym("x", 2, 3);
    org.junit.jupiter.api.Assertions.assertThrows(RuntimeException.class, () -> ca.det(X));
    ca.det(MX.sym("x", 3, 3));
  }

  @Test
  public void test_det() {
    MX X = MX.sym("x", 3, 3);
    SX x = SX.sym("x", 3, 3);

    List<DM> in = dms(DM.rand(3, 3));

    Function f = new Function("f", List.of(x), List.of(ca.det(x)));
    Function F = new Function("F", List.of(X), List.of(ca.det(X)));
    checkfunction(f, F, in);
    checkCodegen(F, in);
  }

  /** Determinant by Gaussian elimination with partial pivoting (reference for numpy.linalg.det) */
  private static double detRef(double[][] a0) {
    int n = a0.length;
    double[][] a = new double[n][];
    for (int i = 0; i < n; i++) a[i] = a0[i].clone();
    double d = 1;
    for (int k = 0; k < n; k++) {
      int p = k;
      for (int i = k + 1; i < n; i++) if (Math.abs(a[i][k]) > Math.abs(a[p][k])) p = i;
      if (a[p][k] == 0) return 0;
      if (p != k) {
        double[] t = a[p];
        a[p] = a[k];
        a[k] = t;
        d = -d;
      }
      d *= a[k][k];
      for (int i = k + 1; i < n; i++) {
        double l = a[i][k] / a[k][k];
        for (int j = k; j < n; j++) a[i][j] -= l * a[k][j];
      }
    }
    return d;
  }

  private double[][] randArray(int n, int m) {
    double[][] a = new double[n][m];
    for (int i = 0; i < n; i++)
      for (int j = 0; j < m; j++) a[i][j] = random.nextDouble();
    return a;
  }

  @Test
  public void test_det_solvers() {
    List<String> solvers = new ArrayList<>(List.of("qr"));
    if (hasLinsol("csparse")) solvers.add("csparse");

    List<double[][]> dense_cases = new ArrayList<>();
    for (int n = 1; n <= 5; n++) dense_cases.add(randArray(n, n));

    double[][] S = randArray(5, 5);
    S[0][2] = S[0][4] = S[2][0] = S[3][1] = S[4][2] = 0;
    for (int i = 0; i < 5; i++) S[i][i] += 3;

    for (double[][] A : dense_cases) {
      int n = A.length;
      DM ref = new DM(detRef(A));
      MX X = MX.sym("x", n, n);
      SX x = SX.sym("x", n, n);

      checkarray(ca.det(new DM(A)), ref, 1e-9, "");
      checkarray(ca.evalf(ca.det(new SX(new DM(A)))), ref, 1e-9, "");

      SX xs = SX.sym("xs", n, n);
      checkarray(ca.evalf(ca.substitute(ca.det(xs, "symbolicqr"), xs, new SX(new DM(A)))), ref, 1e-9, "");

      Function f = new Function("f", List.of(x), List.of(ca.det(x)));
      for (String ls : solvers) {
        Function F = new Function("F", List.of(X), List.of(ca.det(X, ls)));
        checkfunction(f, F, dms(new DM(A)), 1e-8);
        checkarray(ca.det(new DM(A), ls), ref, 1e-9, "");
        if (ls.equals("qr")) checkCodegen(F, dms(new DM(A)));
      }

      checkarray(ca.det(new DM(A), "qr", Map.of("eps", 1e-13)), ref, 1e-9, "");
    }

    DM refS = new DM(detRef(S));
    checkarray(ca.det(new DM(S)), refS, 1e-9, "");
    for (String ls : solvers) checkarray(ca.det(ca.sparsify(new DM(S)), ls), refS, 1e-9, "");

    assertInException("does not provide a symbolic determinant", () -> ca.det(SX.sym("x", 2, 2), "qr"));

    if (hasLinsol("csparse"))
      assertInException("does not provide a symbolic determinant", () -> ca.det(SX.sym("x", 2, 2), "csparse"));
  }

  private static List<Object> sxdet(double[][] Anp) {
    DM Asp = ca.sparsify(new DM(Anp));
    SX x = SX.sym("x", Asp.sparsity());
    Function f = new Function("f", List.of(x), List.of(ca.det(x, "symbolicqr")));
    return List.of(f.call(dms(Asp)).get(0), f.n_nodes());
  }

  @Test
  public void test_det_symbolicqr_btf() {
    int n = 8;
    double[][] B = new double[n][n];
    double[][] b1 = randArray(4, 4), b2 = randArray(4, 4);
    for (int i = 0; i < 4; i++)
      for (int j = 0; j < 4; j++) {
        B[i][j] = b1[i][j];
        B[4 + i][4 + j] = b2[i][j];
      }
    List<Object> rb = sxdet(B);
    long nodes_b = (Long) rb.get(1);
    checkarray((DM) rb.get(0), new DM(detRef(B)), 1e-9, "");

    SX xd = SX.sym("x", n, n);
    long nodes_d = new Function("f", List.of(xd), List.of(ca.det(xd, "symbolicqr"))).n_nodes();
    assertTrue(nodes_b < nodes_d);

    double[][] L = randArray(6, 6);
    for (int i = 0; i < 6; i++) L[i][i] += 3;
    for (int i = 0; i < 3; i++)
      for (int j = 3; j < 6; j++) L[i][j] = 0;
    DM Lsp = ca.sparsify(new DM(L));
    SX xL = SX.sym("x", Lsp.sparsity());
    SX dL = ca.det(xL, "symbolicqr");
    checkarray(ca.evalf(ca.substitute(dL, xL, new SX(Lsp))), new DM(detRef(L)), 1e-9, "");
    long[] rows = Lsp.sparsity().row(), cols = Lsp.sparsity().get_col();
    List<Double> coupling = new ArrayList<>();
    for (int k = 0; k < Lsp.nnz(); k++) if (rows[k] >= 3 && cols[k] < 3) coupling.add((double) k);
    assertTrue(coupling.size() > 0);
    DM cidx = new DM(coupling.stream().mapToDouble(Double::doubleValue).toArray());
    assertEquals(0, ca.jacobian(dL, xL).get(new Slice(), cidx).nnz());

    double[][] U = randArray(5, 5);
    for (int i = 0; i < 5; i++) {
      U[i][i] += 1;
      for (int j = 0; j < i; j++) U[i][j] = 0;
    }
    List<Object> rU = sxdet(U);
    checkarray((DM) rU.get(0), new DM(detRef(U)), 1e-9, "");
    assertTrue((Long) rU.get(1) < 30);

    double[][] Sg = randArray(4, 4);
    for (int j = 0; j < 4; j++) Sg[2][j] = 0;
    checkarray((DM) sxdet(Sg).get(0), new DM(0), 1e-12, "");
  }

  @Test
  public void test_mtimes_mismatch_segfault() {
    assertInException("incompatible dimensions", () -> ca.mtimes(new DM(Sparsity.lower(5)), MX.sym("x", 100)));
  }

  @Test
  public void test_monitor() {
    MX x = MX.sym("x");
    MX y = ca.sqrt(x.monitor("hey"));

    Function f = new Function("f", List.of(x), List.of(y));
    String out = captureStdout(() -> f.invoke(3));
    System.out.println(out);

    String ref = "hey:\n3.0000000000000000e+00\n";
    assertEquals(ref, out);

    if (RUN_SLOW) {
      checkCodegenMainParsingError(f, dms(dm(3)));
      String out2 = readFile(lastMainOut);
      System.out.println(out2);
      assertTrue(out2.contains(ref));
    }

    MX x2 = MX.sym("x", 2);
    MX y2 = ca.sqrt(x2.monitor("hey"));

    Function f2 = new Function("f", List.of(x2), List.of(y2));
    out = captureStdout(() -> f2.invoke(3));
    System.out.println(out);
    ref = "hey:\n2x1: [3.0000000000000000e+00, 3.0000000000000000e+00]\n";
    assertEquals(ref, out);

    if (RUN_SLOW) {
      checkCodegenMainParsingError(f2, dms(dm(3)));
      String out2 = readFile(lastMainOut);
      System.out.println(out2);
      assertTrue(out2.contains(ref));
    }
  }

  private static void removeGlob(String dir, String prefix) {
    File[] fs = new File(dir).listFiles();
    if (fs == null) return;
    for (File f : fs) if (f.getName().startsWith(prefix) && f.getName().endsWith(".mtx")) assertTrue(f.delete());
  }

  private static void rmtree(File d) {
    File[] fs = d.listFiles();
    if (fs != null) for (File f : fs) rmtree(f);
    d.delete();
  }

  @Test
  public void test_dump() {
    removeGlob(".", "myvar.");

    MX x = MX.sym("x");
    MX y = ca.sqrt(x.dump("myvar"));

    Function f = new Function("f", List.of(x), List.of(y));

    f.invoke(9);

    assertTrue(new File("myvar.000000.mtx").exists());
    DM val = DM.from_file("myvar.000000.mtx");
    checkarray(val, 9);

    f.invoke(16);
    assertTrue(new File("myvar.000001.mtx").exists());
    val = DM.from_file("myvar.000001.mtx");
    checkarray(val, 16);

    removeGlob(".", "myvar.");

    x = MX.sym("x", Sparsity.lower(2));
    y = x.dump("mysp").times(2);
    f = new Function("f", List.of(x), List.of(y));

    DM inp = ca.sparsify(dm(new double[][] {{1, 0}, {3, 4}}));
    f.invoke(inp);

    assertTrue(new File("mysp.000000.mtx").exists());
    val = DM.from_file("mysp.000000.mtx");
    checkarray(val, inp);

    removeGlob(".", "mysp.");

    String dump_dir = "dump_test_dir";
    new File(dump_dir).mkdirs();
    x = MX.sym("x");
    y = x.dump("myvar", Map.of("dir", dump_dir));
    f = new Function("f", List.of(x), List.of(y));

    f.invoke(42);

    String dump_path = dump_dir + "/myvar.000000.mtx";
    assertTrue(new File(dump_path).exists());
    val = DM.from_file(dump_path);
    checkarray(val, 42);

    removeGlob(dump_dir, "myvar.");
    assertTrue(new File(dump_dir).delete());

    if (RUN_SLOW && org.casadi.CasadiMeta.feature_list().contains("ghc-filesystem")) {
      for (String d : List.of("mydump", "mydump_codegen")) rmtree(new File(d));

      x = MX.sym("x");
      y = ca.sqrt(x.dump("x", Map.of("dir", "mydump")));
      f = new Function("f", List.of(x), List.of(y));
      checkCodegen(f, dms(dm(9)), Map.of("dump_dir_suffix", "_codegen"), 1e-15, "c99", true);

      assertEquals(readFile("mydump/x.000000.mtx"), readFile("mydump_codegen/x.000000.mtx"));

      for (String d : List.of("mydump", "mydump_codegen")) rmtree(new File(d));
    }
  }

  @Test
  public void test_codegen_specials() {
    MX x = MX.sym("x");
    MX y = MX.sym("y");

    for (MX z : List.of(ca.power(x, 2), ca.if_else(ca.gt(y, 0), x.times(2), x.times(y)), ca.fmin(x, y), ca.fmax(x, y), ca.sign(x.times(y)))) {
      Function f = new Function("f", List.of(x, y), List.of(z));
      checkCodegen(f, dms(dm(1), dm(2)));
      checkCodegen(f, dms(dm(1), dm(-2)));
    }
  }

  private static MX getAny(MX m, Object r, Object c) {
    if (r instanceof Slice && c instanceof Slice) return m.get(false, (Slice) r, (Slice) c);
    if (r instanceof Slice) return m.get(false, (Slice) r, (DM) c);
    if (c instanceof Slice) return m.get(false, (DM) r, (Slice) c);
    return m.get(false, (DM) r, (DM) c);
  }

  private static void setAny(MX m, Object r, Object c, MX v) {
    if (r instanceof Slice && c instanceof Slice) m.set(v, false, (Slice) r, (Slice) c);
    else if (r instanceof Slice) m.set(v, false, (Slice) r, (DM) c);
    else if (c instanceof Slice) m.set(v, false, (DM) r, (Slice) c);
    else m.set(v, false, (DM) r, (DM) c);
  }

  @Test
  public void test_getsetnonzeros() {
    for (Sparsity S : List.of(Sparsity.lower(5), Sparsity.dense(5, 5))) {
      MX M = MX.sym("X", S.nnz());

      MX m = ca.sin(new MX(S, M));
      for (int i : new int[] {0, 2}) {
        Object[][] inds = {{new Slice(i), new Slice(i)}, {new Slice(1), new Slice(i)}, {new Slice(i), new Slice(1)},
          {new Slice(), new Slice(i)}, {new Slice(i), new Slice()}, {new Slice(i, i + 2), new Slice(i, i + 2)},
          {new Slice(i, i + 2), new Slice()}, {dm(i, i), dm(0, 0)}, {new DM(new double[0]), new DM(new double[0])}};
        for (Object[] ind : inds) {
          MX E = getAny(m, ind[0], ind[1]);
          MX e = ca.cos(E);
          e = ca.dot(e, e);
          Function f = new Function("f", List.of(M), List.of(e));
          checkfunction(f, f.expand(), dms(DM.rand(S.nnz(), 1)));
          checkCodegen(f, dms(DM.rand(S.nnz(), 1)), Map.of(), 1e-13, "c89", false);

          MX mc = m.plus(0);

          MX Y = MX.sym("y", E.nnz());
          MX y = new MX(E.sparsity(), Y);
          setAny(mc, ind[0], ind[1], y);
          e = ca.cos(mc);
          e = ca.dot(e, e);

          f = new Function("f", List.of(M, Y), List.of(e));
          checkfunction(f, f.expand(), dms(DM.rand(S.nnz(), 1), DM.rand(E.nnz(), 1)));
          checkSerialize(f, dms(DM.rand(S.nnz(), 1), DM.rand(E.nnz(), 1)));
          checkCodegen(f, dms(DM.rand(S.nnz(), 1), DM.rand(E.nnz(), 1)), Map.of(), 1e-13, "c89", false);
        }
      }
    }
  }

  @Test
  public void test_evalf() {
    MX x = MX.sym("x");

    MX p = MX.sym("p");
    Function f = new Function("f", List.of(x), List.of(ca.sin(x)));
    MX y0 = f.call(List.of(p), false, true).get(0);
    MX y = ca.substitute(y0, p, 3);

    assertInException("not defined", () -> y.to_DM());
    checkarray(ca.evalf(y), ca.sin(new DM(3)));
    assertInException("since variables [x] are free", () -> ca.evalf(x));
  }

  private void mminChecks(java.util.function.Function<Sparsity, Function> make) {
    Function f = make.apply(Sparsity.dense(2, 2));
    List<DM> r = f.call(dms(dm(new double[][] {{-2, 3}, {-3, 5}})));
    checkarray(r.get(0), -3);
    checkarray(r.get(1), 5);

    f = make.apply(Sparsity.lower(2));
    r = f.call(dms(ca.sparsify(dm(new double[][] {{-2, 0}, {-3, -5}}))));
    checkarray(r.get(0), -5);
    checkarray(r.get(1), 0);

    r = f.call(dms(ca.sparsify(dm(new double[][] {{2, 0}, {3, 5}}))));
    checkarray(r.get(0), 0);
    checkarray(r.get(1), 5);

    f = make.apply(new Sparsity(2, 2));
    r = f.call(dms(new DM(2, 2)));
    checkarray(r.get(0), 0);
    checkarray(r.get(1), 0);

    f = make.apply(new Sparsity(2, 0));
    r = f.call(dms(new DM(2, 0)));
    assertTrue(r.get(0).is_empty());
    assertTrue(r.get(1).is_empty());

    f = make.apply(new Sparsity(0, 0));
    r = f.call(dms(new DM(0, 0)));
    assertTrue(r.get(0).is_empty());
    assertTrue(r.get(1).is_empty());
  }

  @Test
  public void test_mmin() {
    java.util.function.Function<Sparsity, Function> sx = sp -> {
      SX x = SX.sym("x", sp);
      return new Function("f", List.of(x), List.of(ca.mmin(x), ca.mmax(x)));
    };
    java.util.function.Function<Sparsity, Function> mx = sp -> {
      MX x = MX.sym("x", sp);
      return new Function("f", List.of(x), List.of(ca.mmin(x), ca.mmax(x)));
    };
    mminChecks(sx);
    mminChecks(sp -> sx.apply(sp).expand());
    mminChecks(sp -> {
      Function f = sx.apply(sp);
      SX y = SX.sym("y", sp);
      return new Function("f", List.of(y), f.call(List.of(y), false, false));
    });
    mminChecks(mx);
    mminChecks(sp -> mx.apply(sp).expand());
    mminChecks(sp -> {
      Function f = mx.apply(sp);
      MX y = MX.sym("y", sp);
      return new Function("f", List.of(y), f.call(List.of(y), true, false));
    });

    for (int k = 0; k < 2; k++) {
      MX x = MX.sym("X", 2);
      Function f = new Function("f", List.of(x), List.of(k == 0 ? ca.mmin(x) : ca.mmax(x)));
      f.expand().disp(true);
      Function J = f.expand().jacobian();
      f.expand().jacobian().disp(true);
      System.out.println(J.call(dms(dm(2, 2), dm(0))));
      checkfunction(f, f.expand(), dms(dm(2, 2)));
      checkCodegen(f, dms(dm(2, 2)), Map.of(), 1e-15, "c99", false);
    }
  }

  // test_doc_expression_tools: python-only (runtime __doc__ docstrings)

  @Test
  public void test_densify() {
    Sparsity I = ca.sparsify(dm(new double[][] {{0, 1, 0, 1}, {1, 1, 0, 1}, {0, 1, 1, 0}})).sparsity();

    MX x = MX.sym("x", I);

    Function f = new Function("f", List.of(x), List.of(ca.densify(x)));

    DM a = dm(new double[][] {{0, 1, 0, 2}, {3, 4, 0, 5}, {0, 6, 7, 0}});

    DM y = (DM) f.invoke(ca.sparsify(a));
    assertTrue(y.sparsity().is_equal(a.sparsity()));
    checkarray(y, a);

    checkCodegen(f, dms(a));

    x = MX.sym("x", 3, 4);

    f = new Function("f", List.of(x), List.of(x.get(I)));

    a = dm(new double[][] {{0, 1, 0, 2}, {3, 4, 0, 5}, {0, 6, 7, 0}});
    DM b = ca.sparsify(a);
    a = dm(new double[][] {{9, 1, 9, 2}, {3, 4, 9, 5}, {9, 6, 7, 9}});

    y = (DM) f.invoke(a);
    assertTrue(y.sparsity().is_equal(b.sparsity()));
    checkarray(y, b);

    checkCodegen(f, dms(a));
    checkSerialize(f, dms(a));
  }

  private static void writeFile(String p, String s) {
    try {
      java.nio.file.Files.writeString(java.nio.file.Path.of(p), s);
    } catch (IOException e) {
      throw new RuntimeException(e);
    }
  }

  @Test
  public void test_constant_from_file() {
    writeFile("test.txt", "5.6 1e-5");
    assertInException("Failed to read a double from 'test.txt'. Expected 3 doubles.", () -> new MX(Sparsity.lower(2), "test.txt"));

    writeFile("test.txt", "5.6 abc 1e-5");
    assertInException("Failed to read a double from 'test.txt'. Expected 3 doubles.", () -> new MX(Sparsity.lower(2), "test.txt"));

    assertInException("Cannot open file 'testQWHAL567p.txt'.", () -> new MX(Sparsity.lower(2), "testQWHAL567p.txt"));

    writeFile("test.txt", "5.6 1e-5 -12");
    MX x = new MX(Sparsity.lower(2), "test.txt");

    MX a = MX.sym("a");

    Function f = new Function("f", List.of(a), List.of(a.times(x)));
    checkarray((DM) f.invoke(2), new DM(Sparsity.lower(2), dm(2 * 5.6, 2 * 1e-5, -2 * 12)));

    checkCodegen(f, dms(dm(2)));

    assertInException("Not defined for ConstantFile", () -> x.to_DM());
    new File("test.txt").delete();
  }

  @Test
  public void test_nonzeros_param() {
    MX x = MX.sym("x", 2);
    MX y = MX.sym("y");

    for (int ad_weight_sp : new int[] {0, 1}) {
      Map<String, Object> opts = Map.of("helper_options", Map.of("ad_weight_sp", ad_weight_sp));
      MX z = x.get_nz(false, y);
      MX J = ca.jacobian(z, x, opts);
      assertEquals(1, J.size1());
      assertEquals(2, J.size2());
      assertEquals(2, J.nnz());
      J = ca.jacobian(z, y, opts);
      assertEquals(1, J.size1());
      assertEquals(1, J.size2());
      assertEquals(0, J.nnz());

      MX q = MX.sym("q");
      z = new MX(x);
      z.set_nz(q, false, y);

      J = ca.jacobian(z, x, opts);
      assertEquals(2, J.size1());
      assertEquals(2, J.size2());
      assertEquals(4, J.nnz());
      J = ca.jacobian(z, y, opts);
      assertEquals(2, J.size1());
      assertEquals(1, J.size2());
      assertEquals(0, J.nnz());
      J = ca.jacobian(z, q, opts);
      assertEquals(2, J.size1());
      assertEquals(1, J.size2());
      assertEquals(2, J.nnz());
    }
  }

  private static List<String> lines(String s, int drop) {
    List<String> r = new ArrayList<>(List.of(s.split("\n", -1)));
    return r.subList(0, r.size() - drop);
  }

  private static List<String> printInstructionsCheck(Function f) {
    DM.rng(1);
    List<DM> inputs = new ArrayList<>();
    for (long i = 0; i < f.n_in(); i++) inputs.add(DM.rand(f.sparsity_in(i)));
    List<String> out = lines(captureStdout(() -> f.call(inputs)), 1);

    assertTrue(out.get(0).contains("inputs"));

    if (RUN_SLOW) {
      checkCodegenMainParsingError(f, inputs);
      List<String> out2 = lines(readFile(lastMainOut), 2);
      assertEquals(out, out2);
      return out2;
    }
    return null;
  }

  @Test
  public void test_print_instructions() {
    Map<String, Object> popts = Map.of("print_instructions", true, "print_canonical", true);
    {
      MX x = MX.sym("x", Sparsity.diag(5)), y = MX.sym("y", 5, 5), z = MX.sym("z", 5), w = MX.sym("w");
      Function g = new Function("g", List.of(z), List.of(z.get(0).plus(z.get(2)), z.get(1).plus(z.get(4))), Map.of("never_inline", true));
      List<MX> ab = g.call(List.of(ca.mtimes(y.plus(x), z)));
      Function f = new Function("f", List.of(x, y, z, w), List.of(ca.sin(ab.get(0).times(w.times(7)).plus(6)).times(ab.get(1))), popts);
      printInstructionsCheck(f);
    }
    {
      SX x = SX.sym("x", Sparsity.diag(5)), y = SX.sym("y", 5, 5), z = SX.sym("z", 5), w = SX.sym("w");
      Function g = new Function("g", List.of(z), List.of(z.get(0).plus(z.get(2)), z.get(1).plus(z.get(4))), Map.of("never_inline", true));
      List<SX> ab = g.call(List.of(ca.mtimes(y.plus(x), z)));
      Function f = new Function("f", List.of(x, y, z, w), List.of(ca.sin(ab.get(0).times(w.times(7)).plus(6)).times(ab.get(1))), popts);
      printInstructionsCheck(f);
    }

    Map<String, List<String>> results = new LinkedHashMap<>(), results_slow = new LinkedHashMap<>();
    for (String X : List.of("MX", "SX")) {
      Function f;
      if (X.equals("MX")) {
        MX x = MX.sym("x"), y = MX.sym("y");
        f = new Function("f", List.of(x, y), List.of(x.plus(y)), popts);
      } else {
        SX x = SX.sym("x"), y = SX.sym("y");
        f = new Function("f", List.of(x, y), List.of(x.plus(y)), popts);
      }
      f.disp(true);
      DM.rng(1);
      List<DM> inputs = new ArrayList<>();
      for (long i = 0; i < f.n_in(); i++) inputs.add(DM.rand(f.sparsity_in(i)));
      List<String> out = lines(captureStdout(() -> f.call(inputs)), 1);
      results.put(X, out);

      assertTrue(out.get(0).contains("inputs"));

      if (RUN_SLOW) {
        checkCodegenMainParsingError(f, inputs);
        List<String> out2 = lines(readFile(lastMainOut), 2);
        results_slow.put(X, out2);
        assertEquals(out, out2);
      }
    }

    assertEquals(results.get("SX"), results.get("MX"));
    if (RUN_SLOW) assertEquals(results_slow.get("SX"), results_slow.get("MX"));
  }

  @Test
  public void test_low() {
    MX v = MX.sym("v", 5);
    double[] p0 = {-1, 0, 0.2, 1, 1.3, 2, 2.5, 3, 3.5, 4, 4.5};
    MX p = MX.sym("p", p0.length);

    Function f = new Function("f", List.of(v, p), List.of(ca.low(v, p)));

    List<DM> inputs = dms(range(0, 5), dm(p0));
    checkCodegen(f, inputs);
    checkSerialize(f, inputs);

    checkarray(f.call(inputs).get(0), dm(0, 0, 0, 1, 1, 2, 2, 3, 3, 3, 3));
  }

  @Test
  public void test_linear_interpn() {
    int N = 4, nq = 3, n_out = 2;
    double[] g = new double[N];
    for (int i = 0; i < N; i++) g[i] = random.nextDouble();
    java.util.Arrays.sort(g);
    double[] x0 = new double[nq];
    for (int i = 0; i < nq; i++) x0[i] = g[0] + (g[N - 1] - g[0]) * (i + 1) / (nq + 1);

    double[] d = new double[n_out * N];
    for (int i = 0; i < d.length; i++) d[i] = random.nextDouble();

    Function F = ca.interpolant("F", "linear", new double[][] {g}, d);

    List<DM> rows = new ArrayList<>();
    for (double xi : x0) rows.add(((DM) F.invoke(xi)).T());
    DM ref = ca.vertcat(rows);

    MX x = MX.sym("x0", 1, N);
    MX v = MX.sym("v", N * n_out);
    MX xq = MX.sym("xq0", nq);

    MX e = MX.interpn_linear(List.of(x), v, List.of(xq));

    Function f = new Function("f", List.of(x, v, xq), List.of(e));
    List<DM> inputs = dms(dm(g).T(), dm(d), dm(x0));
    DM r = f.call(inputs).get(0);
    checkarray(r, ref);
    checkCodegen(f, inputs, Map.of(), 1e-15, "c99", false);
    checkSerialize(f, inputs);
  }

  private static DM sliceVals(Slice s, long len) {
    long[] a = s.all(len);
    double[] r = new double[a.length];
    for (int i = 0; i < a.length; i++) r[i] = a[i];
    return new DM(r);
  }

  private static MX get2(MX m, Object r, Object c) {
    if (r instanceof Slice && c instanceof Slice) return m.get(false, (Slice) r, (Slice) c);
    if (r instanceof Slice) return m.get(false, (Slice) r, (MX) c);
    if (c instanceof Slice) return m.get(false, (MX) r, (Slice) c);
    return m.get(false, (MX) r, (MX) c);
  }

  @Test
  public void test_paramslice() {
    int N = 4, M = 6;
    DM A = DM.rand(N, M);
    MX As = MX.sym("A", A.size1(), A.size2());
    MX ASs = MX.sym("AS", ca.tril(As.sparsity()));
    ASs.sparsity().spy();
    for (Slice r : List.of(new Slice(1), new Slice(1, 3), new Slice(0, N, 2), new Slice())) {
      MX rs = MX.sym("r", A.get(r, new Slice()).size1());
      for (Slice c : List.of(new Slice(1), new Slice(1, 3), new Slice(0, M, 2), new Slice())) {
        MX cs = MX.sym("c", A.get(new Slice(), c).size2());
        for (Object rr : List.of(r, rs)) {
          for (Object cc : List.of(c, cs)) {
            Function f = new Function("f", List.of(As, rs, cs), List.of(get2(As, rr, cc)));
            checkarray(f.call(dms(A, sliceVals(r, N), sliceVals(c, M))).get(0), A.get(r, c));
            if (rr instanceof MX || cc instanceof MX)
              assertInException("Parametric slicing only supported for dense matrices", () -> get2(ASs, rr, cc));
          }
        }
      }
    }
    for (Slice k : List.of(new Slice(1), new Slice(1, 3), new Slice(0, N, 2), new Slice(0, 2 * N, 2), new Slice())) {
      MX ks = MX.sym("r", A.get(k).size1());
      for (Object kk : List.of(k, ks)) {
        Function f = new Function("f", List.of(As, ks), List.of(kk instanceof Slice ? As.get((Slice) kk) : As.get(false, (MX) kk)));
        checkarray(f.call(dms(A, sliceVals(k, N * M))).get(0), A.get(k));
        assertInException("Parametric slicing only supported for dense matrices", () -> ASs.get(false, ks));
      }
    }
  }

  @Test
  public void test_mapping() {
    MX A = MX.sym("A", 4, 4);
    DM i = dm(new double[][] {{0, 3}, {1, 2}});
    checkarray(i, A.get(i).mapping());
  }

  /** Symmetric eigendecomposition by cyclic Jacobi (reference for numpy.linalg.eigh): {D (1 x n), V (n x n)} */
  private static double[][][] eigh(double[][] A0) {
    int n = A0.length;
    double[][] a = new double[n][];
    for (int i = 0; i < n; i++) a[i] = A0[i].clone();
    double[][] v = new double[n][n];
    for (int i = 0; i < n; i++) v[i][i] = 1;
    for (int sweep = 0; sweep < 100; sweep++) {
      double off = 0;
      for (int i = 0; i < n; i++)
        for (int j = i + 1; j < n; j++) off += a[i][j] * a[i][j];
      if (off < 1e-30) break;
      for (int p = 0; p < n; p++)
        for (int q = p + 1; q < n; q++) {
          if (a[p][q] == 0) continue;
          double theta = (a[q][q] - a[p][p]) / (2 * a[p][q]);
          double t = Math.signum(theta) / (Math.abs(theta) + Math.sqrt(theta * theta + 1));
          if (theta == 0) t = 1;
          double c = 1 / Math.sqrt(t * t + 1), s = t * c;
          for (int k = 0; k < n; k++) {
            double akp = a[k][p], akq = a[k][q];
            a[k][p] = c * akp - s * akq;
            a[k][q] = s * akp + c * akq;
          }
          for (int k = 0; k < n; k++) {
            double apk = a[p][k], aqk = a[q][k];
            a[p][k] = c * apk - s * aqk;
            a[q][k] = s * apk + c * aqk;
          }
          for (int k = 0; k < n; k++) {
            double vkp = v[k][p], vkq = v[k][q];
            v[k][p] = c * vkp - s * vkq;
            v[k][q] = s * vkp + c * vkq;
          }
        }
    }
    double[][] d = new double[1][n];
    for (int i = 0; i < n; i++) d[0][i] = a[i][i];
    return new double[][][] {d, v};
  }

  /** V diag(g(D)) V' */
  private static DM eigRecompose(double[][][] DV, java.util.function.DoubleUnaryOperator g) {
    double[] D = DV[0][0];
    double[] Dm = new double[D.length];
    for (int i = 0; i < D.length; i++) Dm[i] = g.applyAsDouble(D[i]);
    DM V = new DM(DV[1]);
    return ca.mtimes(ca.mtimes(V, ca.diag(new DM(Dm))), V.T());
  }

  @Test
  public void test_convexify() {
    DM A = ca.diagcat(List.of(new DM(1), new DM(2), new DM(-1), dm(new double[][] {{1.2, 1.3}, {1.3, 4}}),
        ca.sparsify(dm(new double[][] {{0, 1, 0}, {1, 4, 7}, {0, 7, 9}})), new DM(2, 2)));

    double margin = 1e-7;

    List<Double> perm = new ArrayList<>();
    for (int i = 0; i < A.size1(); i++) perm.add((double) i);
    java.util.Collections.shuffle(perm, random);
    DM p = new DM(perm.stream().mapToDouble(Double::doubleValue).toArray());

    double[][][] DV = eigh(A.full());
    DM Ar_ref = eigRecompose(DV, d -> Math.max(Math.abs(d), 1e-7));
    DM Ac_ref = eigRecompose(DV, d -> Math.max(d, 1e-7));
    MX As = MX.sym("As", A.sparsity());

    List<Map<String, Object>> optss = List.of(Map.of("strategy", "eigen-reflect"), Map.of("strategy", "eigen-clip"), Map.of("strategy", "regularize"));
    List<DM> refs = List.of(Ar_ref, Ac_ref, A.plus(DM.eye(A.size1()).times(4 + margin)));
    for (int k = 0; k < 3; k++) {
      Map<String, Object> opts = optss.get(k);
      DM ref = refs.get(k);
      List<UnaryOperator<DM>> ops = List.of(e -> e, e -> ca.triu(e), e -> ca.tril(e),
          e -> e.get(p, p), e -> ca.triu(e.get(p, p)), e -> ca.tril(e.get(p, p)));
      if (opts.toString().contains("regularize")) ops = List.of(e -> e, e -> e.get(p, p));
      for (UnaryOperator<DM> op : ops) {
        Function f = new Function("f", List.of(As), List.of(ca.convexify(op.apply(A), opts)));

        checkarray((DM) f.invoke(A), op.apply(ref), 1e-8, "");

        checkSerialize(f, dms(A));

        checkCodegen(f, dms(A));
      }
    }
  }

  @Test
  public void test_convexify_bugs() {
    for (double eps : new double[] {0, 1e-100, 1}) {
      double q = 4;
      double Q = 10;
      for (double s : new double[] {-4, 0.1, 0, 0.1, 4}) {
        for (UnaryOperator<DM> trans : List.<UnaryOperator<DM>>of(e -> e, e -> ca.sparsify(e))) {
          DM A = trans.apply(dm(new double[][] {{q, s, 0 + eps}, {s, Q, 0}, {0 + eps, 0, 2}}));
          System.out.println(A);
          for (double d : eigh(A.full())[0][0]) assertTrue(d > 0);
          DM Ac = ca.evalf(ca.convexify(A, Map.of("strategy", "eigen-reflect")));
          checkarray(A, Ac, 1e-8, "");
        }
      }
    }
  }

  @Test
  public void test_convexify_large() {
    for (int n : new int[] {3, 120}) {
      DM B = randDM(n, n);
      DM A = B.plus(B.T());
      double[][][] DV = eigh(A.full());
      MX As = MX.sym("As", A.sparsity());
      for (String strategy : List.of("eigen-reflect", "eigen-clip")) {
        DM ref = eigRecompose(DV, strategy.equals("eigen-reflect") ? d -> Math.max(Math.abs(d), 1e-7) : d -> Math.max(d, 1e-7));
        Function f = new Function("f", List.of(As), List.of(ca.convexify(As, Map.of("strategy", strategy, "max_iter_eig", 2000))));
        checkarray((DM) f.invoke(A), ref, 1e-8, "");
        if (n < 10) checkCodegen(f, dms(A));
      }
    }
  }

  private static DM idx(double... i) {
    return new DM(i);
  }

  private static DM revRange(long n) {
    double[] r = new double[(int) n];
    for (int i = 0; i < n; i++) r[i] = n - 1 - i;
    return new DM(r);
  }

  private static List<MX> repmatSlices(MX a, long m, long n) {
    MX Y = ca.repmat(a, m, n);
    long R = Y.size1(), C = Y.size2();
    return List.of(Y.get(new Slice(), revRange(C)), Y.get(revRange(R), new Slice()), Y.get(idx(R - 1, 0), idx(C - 1, 0)),
        Y.get_nz(false, revRange(Y.nnz())), Y.T());
  }

  private static List<SX> repmatSlices(SX a, long m, long n) {
    SX Y = ca.repmat(a, m, n);
    long R = Y.size1(), C = Y.size2();
    return List.of(Y.get(new Slice(), revRange(C)), Y.get(revRange(R), new Slice()), Y.get(idx(R - 1, 0), idx(C - 1, 0)),
        Y.get_nz(false, revRange(Y.nnz())), Y.T());
  }

  @Test
  public void test_repmat_slice() {
    MX x = MX.sym("x", 2);
    MX y = ca.repmat(x, 1, 2);
    assertTrue(y.get(new Slice(), 0).minus(y.get(new Slice(), 1)).is_zero());
    assertTrue(ca.is_equal(y.get(new Slice(), 1), x));
    SX xs = SX.sym("x", 2);
    MX X = MX.sym("X", Sparsity.lower(2));
    SX Xs = SX.sym("X", Sparsity.lower(2));
    List<MX> as = List.of(x, X);
    List<SX> bs = List.of(xs, Xs);
    for (int k = 0; k < 2; k++) {
      MX a = as.get(k);
      SX b = bs.get(k);
      Function f = new Function("f", List.of(a), List.of(ca.repmat(a, 1, 3).get(new Slice(), idx(2, 0, 1)),
          ca.repmat(a, 1, 3).get(idx(1, 0), idx(1, 2)), ca.repmat(a, 1, 3).get_nz(false, idx(5, 0, 2, 1))));
      Function f_ref = new Function("f", List.of(b), List.of(ca.repmat(b, 1, 3).get(new Slice(), idx(2, 0, 1)),
          ca.repmat(b, 1, 3).get(idx(1, 0), idx(1, 2)), ca.repmat(b, 1, 3).get_nz(false, idx(5, 0, 2, 1))));
      checkfunction(f, f_ref, dms(new DM(a.sparsity(), DM.rand(a.nnz(), 1))));
    }
    List<Sparsity> sps = List.of(Sparsity.dense(2, 3), Sparsity.triplet(3, 3, new long[] {0, 2, 1}, new long[] {0, 0, 2}),
        Sparsity.triplet(3, 2, new long[] {0, 2, 2}, new long[] {0, 0, 1}), new Sparsity(2, 3));
    for (Sparsity sp : sps) {
      for (long[] mn : new long[][] {{1, 3}, {2, 2}}) {
        long m = mn[0], n = mn[1];
        MX a = MX.sym("a", sp);
        SX b = SX.sym("a", sp);
        MX Y = ca.repmat(a, m, n);
        long w = Y.size2() / n;
        assertTrue(Y.get(new Slice(), new Slice(w * (n - 1), Y.size2())).minus(Y.get(new Slice(), new Slice(0, w))).is_zero());
        Function f = new Function("f", List.of(a), repmatSlices(a, m, n));
        Function f_ref = new Function("f", List.of(b), repmatSlices(b, m, n));
        checkfunction(f, f_ref, dms(new DM(sp, DM.rand(sp.nnz(), 1))));
      }
    }
  }

  @Test
  public void test_logsumexp() {
    MX x = MX.sym("x", 3);

    Function f_ref = new Function("f_ref", List.of(x), List.of(ca.log(ca.exp(x.get(0)).plus(ca.exp(x.get(1))).plus(ca.exp(x.get(2))))));
    Function f = new Function("f", List.of(x), List.of(ca.logsumexp(x)));

    checkfunction(f, f_ref, dms(dm(1.1, 1.3, 1.7)));
    checkCodegen(f, dms(dm(1.1, 1.3, 1.7)));
    checkfunction(f, f_ref, dms(dm(1.1, 1.3, 1.3)));
    checkfunction(f, f_ref, dms(dm(1.3, 1.3, 1.3)));

    DM res = (DM) f.invoke(dm(100, 1000, 10000));
    checkarray(res, 10000);
  }

  @Test
  public void test_empty_broadcast() {
    for (long nc : new long[] {0, 2}) {
      MX res = ca.atan2(MX.sym("c", nc, 1), MX.sym("t", nc, 3));

      assertEquals(nc, res.size1());
      assertEquals(3, res.size2());

      assertInException("Dimension mismatch", () -> ca.atan2(MX.sym("c", nc, 2), MX.sym("t", nc, 3)));

      res = ca.atan2(MX.sym("c", nc, 2), MX.sym("t", nc, 4));
      assertEquals(nc, res.size1());
      assertEquals(4, res.size2());

      res = ca.atan2(MX.sym("c", nc, 4), MX.sym("t", nc, 2));
      assertEquals(nc, res.size1());
      assertEquals(4, res.size2());
    }
  }

  @Test
  public void test_horzcat_sparsity() {
    MX x = MX.sym("x", 2, 2);
    Function f = null, f_alt = null;

    for (MX y : List.of(MX.sym("y", new Sparsity(2, 2)), MX.sym("y", Sparsity.lower(2)))) {
      MX z = ca.sin(ca.horzcat(List.of(x, ca.sin(y), x)));
      MX z_alt = ca.sin(ca.sparsity_cast(ca.vertcat(List.of(ca.vec(x), ca.vec(ca.sin(y)), ca.vec(x))), z.sparsity()));
      System.out.println(z_alt);

      f = new Function("f", List.of(x, y), List.of(z));
      f_alt = new Function("f_alt", List.of(x, y), List.of(z));

      for (Function F : List.of(f, f_alt)) {
        assertEquals(z.nnz(), F.call(List.of(x, y), true, false).get(0).nnz());
        assertEquals(z.nnz(), F.call(List.of(x, y), false, true).get(0).nnz());
        assertEquals(z.nnz(), F.call(List.of(x, ca.power(x, 2)), true, false).get(0).nnz());
        assertEquals(z.nnz(), F.call(List.of(x, ca.power(x, 2)), false, true).get(0).nnz());
      }
    }

    List<DM> in = dms(dm(new double[][] {{1, 2}, {3, 4}}), dm(new double[][] {{5, 6}, {7, 8}}));
    checkfunction(f, f.expand(), in);
    checkfunction(f, f_alt, in);
    checkCodegen(f, in);
  }

  @Test
  public void test_sparsity_cast_ad() {
    MX y = MX.sym("y", 3);

    MX xy = ca.vertcat(List.of(y.get(0), y.get(1)));
    MX w = ca.mtimes(new MX(ca.sparsify(dm(1, 1, 0)).sparsity().T(), xy), y);

    Function f = new Function("f", List.of(y), List.of(w));
    checkfunction(f, f.expand(), dms(dm(1.1, 1.3, 1.7)));
    checkCodegen(f, dms(dm(1.1, 1.3, 1.7)));
  }

  @Test
  public void test_fractional_slicing() {
    MX t = MX.sym("x");

    Function f = new Function("f", List.of(t), List.of(new MX(dm(1, 2)).get(false, t)));

    checkarray((DM) f.invoke(0), 1);
    checkarray((DM) f.invoke(1), 2);
    checkarray((DM) f.invoke(0.5), 1);
    checkarray((DM) f.invoke(0.9), 1);

    f = new Function("f", List.of(t), List.of(new MX(dm(new double[][] {{1, 2}, {3, 4}})).get(false, new Slice(), t)));

    checkarray((DM) f.invoke(0), dm(1, 3));
    checkarray((DM) f.invoke(1), dm(2, 4));
    checkarray((DM) f.invoke(0.5), dm(1, 3));
    checkarray((DM) f.invoke(0.9), dm(1, 3));
  }

  private static void checkfunctionLight(Function f, Function g, List<DM> inputs) {
    List<DM> a = f.call(inputs), b = g.call(inputs);
    assertEquals(a.size(), b.size());
    for (int i = 0; i < a.size(); i++) {
      if (f.sparsity_out(i).is_empty() || g.sparsity_out(i).is_empty()) continue;
      checkarray(a.get(i), b.get(i), 1e-9, "output " + i);
    }
  }

  @Test
  public void test_primitives() {
    MX x = MX.sym("x", 3);
    MX y = MX.sym("y", 3);
    List<MX> experiments = List.of(ca.vertcat(List.of(x, y)), ca.horzcat(List.of(x, y)), ca.diagcat(List.of(x, y)),
        ca.vertcat(List.of(x, y, new MX(0, 1))), ca.reshape(ca.vertcat(List.of(x, y)), 2, 3),
        ca.sparsity_cast(ca.vertcat(List.of(x, y)), Sparsity.diag(6)));

    for (MX z : experiments) {
      System.out.println(z);
      assertTrue(z.is_valid_input());

      System.out.println(z.primitives());
      Function f = new Function("f", List.of(z), z.primitives());
      DM.rng(1);
      DM a = DM.rand(z.sparsity());

      List<DM> ref = f.call(dms(a));
      List<DM> res = z.split_primitives(a);

      assertEquals(ref.size(), res.size());
      for (int i = 0; i < ref.size(); i++) checkarray(ref.get(i), res.get(i));

      DM a_back = z.join_primitives(res);
      checkarray(a, a_back);

      SX zsx = SX.sym("x", z.sparsity());
      Function f2 = new Function("f", List.of(zsx), z.split_primitives(zsx));

      checkfunctionLight(f, f2, dms(a));
      checkCodegen(f, dms(a));
      checkCodegen(f2, dms(a));

      Function f3 = new Function("f", List.of(zsx), List.of(z.join_primitives(z.split_primitives(zsx))));

      DM res3 = (DM) f3.invoke(a);
      checkarray(a, res3);
    }
  }

  /** Python hash(MX) */
  private static long hashOf(MX x) {
    return x.hash();
  }

  @Test
  public void test_cache_output_node() {
    MX x = MX.sym("x");
    MX y = MX.sym("y");
    Function f = new Function("f", List.of(x, y), List.of(x.plus(y), x.div(y)));

    List<MX> ab = f.call(List.of(ca.sin(x), ca.cos(y)));
    MX a = ab.get(0), b = ab.get(1);
    MX base = a.dep(0);
    assertTrue(hashOf(a) != 0, String.valueOf(hashOf(base.get_output(0))));
    assertTrue(hashOf(b) != 0, String.valueOf(hashOf(base.get_output(1))));
  }

  @Test
  public void test_naked_call_function() {
    MX x = MX.sym("x");
    MX y = MX.sym("y");
    Function f = new Function("f", List.of(x, y), List.of(x.plus(y), x.div(y)));

    List<MX> ab = f.call(List.of(ca.sin(x), ca.cos(y)));

    MX base = ab.get(0).dep(0);

    assertInException("get_output", () -> new Function("g", List.of(x, y), List.of(base)));

    assertInException("get_output", () -> base.times(3));
  }

  @Test
  public void test_substitute() {
    MX x = MX.sym("x");
    MX y = MX.sym("y");
    MX w = MX.sym("w");

    MX xy = x.times(y);

    MX z = xy.plus(w);
    MX z2 = ca.substitute(List.of(z), List.of(w), List.of(new MX(3))).get(0);
    System.out.println(xy);
    System.out.println(z2.dep(1));
    assertTrue(hashOf(xy) != hashOf(z2.dep(1)));
  }

  /** The SX twin of MX expressions over syms: evaluate them symbolically on SX symbols */
  @SuppressWarnings("unchecked")
  private static List<org.casadi.GenericMatrix> conv(boolean isMX, List<MX> syms, List<org.casadi.GenericMatrix> tsyms, List<MX> exprs) {
    if (isMX) return new ArrayList<>(exprs);
    List<SX> s = new ArrayList<>();
    for (org.casadi.GenericMatrix e : tsyms) s.add((SX) e);
    return new ArrayList<>(new Function("conv", syms, exprs).call(s));
  }

  /** ms itself for MX, else SX symbols of the same names and sparsities */
  private static List<org.casadi.GenericMatrix> twins(boolean isMX, List<MX> ms) {
    List<org.casadi.GenericMatrix> r = new ArrayList<>();
    for (MX m : ms) r.add(isMX ? m : SX.sym(m.name(), m.sparsity()));
    return r;
  }

  private static org.casadi.GenericMatrix add(org.casadi.GenericMatrix a, org.casadi.GenericMatrix b) {
    return a instanceof MX ? ((MX) a).plus((MX) b) : ((SX) a).plus((SX) b);
  }

  private static List<Object> separateLinear(org.casadi.GenericMatrix e, org.casadi.GenericMatrix x, org.casadi.GenericMatrix p) {
    return e instanceof MX ? ca.separate_linear((MX) e, (MX) x, (MX) p) : ca.separate_linear((SX) e, (SX) x, (SX) p);
  }

  private static boolean isLinear(Object e, Object x) {
    return e instanceof MX ? ca.is_linear((MX) e, (MX) x) : ca.is_linear((SX) e, (SX) x);
  }

  private static boolean dependsOn(Object e, Object x) {
    return e instanceof MX ? ca.depends_on((MX) e, (MX) x) : ca.depends_on((SX) e, (SX) x);
  }

  private static List<Object> extractParametric(Object e, Object p, Map<String, ?> opts) {
    return e instanceof MX ? ca.extract_parametric((MX) e, (MX) p, opts) : ca.extract_parametric((SX) e, (SX) p, opts);
  }

  private static String nameOf(Object e) {
    return e instanceof MX ? ((MX) e).name() : ((SX) e).name();
  }

  /** Same values at random inputs */
  private static boolean numEqual(List<org.casadi.GenericMatrix> syms, List<org.casadi.GenericMatrix> a, List<org.casadi.GenericMatrix> b) {
    Function f1 = new Function("f", syms, a);
    Function f2 = new Function("f", syms, b);
    DM.rng(1);
    List<DM> args = new ArrayList<>();
    for (long i = 0; i < f1.n_in(); i++) args.add(DM.rand(f1.sparsity_in(i)));
    List<DM> r1 = f1.call(args, false, false), r2 = f2.call(args, false, false);
    for (int i = 0; i < Math.min(r1.size(), r2.size()); i++)
      if (ca.norm_1(ca.vec(ca.densify(r1.get(i).minus(r2.get(i))))).to_double() > 1e-12) return false;
    return true;
  }

  @SuppressWarnings("unchecked")
  private static List<org.casadi.GenericMatrix> gl(Object o) {
    return o instanceof List ? new ArrayList<>((List<org.casadi.GenericMatrix>) o) : new ArrayList<>(List.of((org.casadi.GenericMatrix) o));
  }

  @Test
  public void test_extract_parametric() {
    MX x = MX.sym("x"), y = MX.sym("y"), p = MX.sym("p");
    List<MX> ms = List.of(x, y, p);
    for (boolean isMX : new boolean[] {false, true}) {
      List<org.casadi.GenericMatrix> ts = twins(isMX, ms);
      org.casadi.GenericMatrix tp = ts.get(2);

      List<org.casadi.GenericMatrix> e = conv(isMX, ms, ts, List.of(x.times(2).times(p).plus(y.times(2).times(ca.power(p, 2))), ca.power(p, 2),
          ca.vertcat(List.of(x.times(2).plus(y.times(3)).plus(p.times(8)), x.times(3).plus(y.times(2)).plus(4), new MX(7))), p.times(8)));
      List<Object> r = extractParametric(e.get(0), tp, Map.of("extract_trivial", true));
      System.out.println(r.get(0) + " " + r.get(2));
      assertTrue(numEqual(ts, gl(r.get(2)), List.of(tp, e.get(1))));

      r = extractParametric(e.get(2), tp, Map.of("extract_trivial", true));

      System.out.println(r.get(0) + " " + r.get(2));

      assertTrue(numEqual(ts, List.of(gl(r.get(2)).get(0)), List.of(e.get(3))));
    }

    for (boolean isMX : new boolean[] {false, true}) {
      List<org.casadi.GenericMatrix> ts = twins(isMX, ms);
      org.casadi.GenericMatrix tp = ts.get(2);

      MX cp = ca.cos(p);
      MX expr = x.minus(ca.sqrt(p)).times(y.plus(ca.power(cp, 2))).plus(cp.div(x.plus(cp)).div(y.minus(p)));

      for (org.casadi.GenericMatrix ex : conv(isMX, ms, ts, List.of(expr, ca.vertcat(List.of(expr, ca.power(p, 2)))))) {
        List<Object> r = extractParametric(ex, tp, Map.of("extract_trivial", true));

        System.out.println(r.get(0) + " " + r.get(1) + " " + r.get(2));

        assertFalse(dependsOn(r.get(0), tp));

        org.casadi.GenericMatrix expr_recreated = ca.substitute(gl(r.get(0)), gl(r.get(1)), gl(r.get(2))).get(0);
        assertTrue(numEqual(ts, List.of(ex), List.of(expr_recreated)));
      }
    }
  }

  private static void separateLinearCases(boolean isMX, List<MX> ms, List<org.casadi.GenericMatrix> ts, org.casadi.GenericMatrix xy,
      List<MX[]> cases) {
    org.casadi.GenericMatrix tp = ts.get(2);
    for (MX[] c : cases) {
      List<org.casadi.GenericMatrix> e = conv(isMX, ms, ts, List.of(c));
      System.out.println("ref " + e);

      List<Object> r = separateLinear(e.get(0), xy, tp);
      System.out.println("actual " + r);
      org.casadi.GenericMatrix expr_const = (org.casadi.GenericMatrix) r.get(0), expr_lin = (org.casadi.GenericMatrix) r.get(1),
          expr_nonlin = (org.casadi.GenericMatrix) r.get(2);

      assertTrue(numEqual(ts, List.of(e.get(0)), List.of(add(add(expr_const, expr_lin), expr_nonlin))));

      assertTrue(numEqual(ts, List.of(expr_const), List.of(e.get(1))));
      assertTrue(numEqual(ts, List.of(expr_lin), List.of(e.get(2))));
      assertTrue(numEqual(ts, List.of(expr_nonlin), List.of(e.get(3))));

      assertTrue(isLinear(expr_lin, xy));
      assertFalse(dependsOn(expr_const, xy));
    }
  }

  @Test
  public void test_separate_linear() {
    for (boolean isMX : new boolean[] {false, true}) {
      MX x = MX.sym("x"), y = MX.sym("y"), p = MX.sym("p");
      MX O = new MX(0);
      MX cp = ca.cos(p);
      MX xy3 = x.times(3).plus(y);

      List<MX[]> cases = List.of(
          new MX[] {x, O, x, O},
          new MX[] {x.plus(p), p, x, O},
          new MX[] {x.plus(y), O, x.plus(y), O},
          new MX[] {x.minus(p), p.unaryMinus(), x, O},
          new MX[] {p.minus(x), p, x.unaryMinus(), O},
          new MX[] {x.times(2), O, x.times(2), O},
          new MX[] {x.unaryMinus(), O, x.unaryMinus(), O},
          new MX[] {xy3, O, xy3, O},
          new MX[] {cp, cp, O, O},
          new MX[] {xy3.plus(cp), cp, xy3, O},
          new MX[] {x.times(3).plus(7), new MX(7), x.times(3), O},
          new MX[] {p.times(x.times(3).plus(7)), p.times(7), p.times(x.times(3)), O},
          new MX[] {y.times(x.times(3).plus(7)), O, y.times(7), y.times(x.times(3))},
          new MX[] {ca.cos(ca.plus(7, x.times(3))), O, O, ca.cos(ca.plus(7, x.times(3)))},
          new MX[] {ca.cos(x.times(3).plus(7)), O, O, ca.cos(x.times(3).plus(7))},
          new MX[] {ca.cos(p.times(7)), ca.cos(p.times(7)), O, O},
          new MX[] {ca.cos(p.plus(x)), O, O, ca.cos(p.plus(x))},
          new MX[] {x.times(cp).plus(7), new MX(7), x.times(cp), O},
          new MX[] {y.times(x.times(cp).plus(7)), O, y.times(7), y.times(x).times(cp)},
          new MX[] {p.times(x.plus(y).plus(x.times(y))), O, p.times(x.plus(y)), p.times(x.times(y))},
          new MX[] {ca.cos(x.plus(y).plus(x.times(y))), O, O, ca.cos(x.plus(y).plus(x.times(y)))},
          new MX[] {p.times(p.plus(cp.times(x))), ca.power(p, 2), p.times(cp).times(x), O},
          new MX[] {p.plus(x).plus(x.times(y)).times(cp.plus(x.times(7)).plus(x.times(3).times(y))), p.times(cp),
            p.times(x.times(7)).plus(x.times(cp)),
            p.times(x.times(3).times(y)).plus(x.times(y).times(cp.plus(x.times(7)).plus(x.times(3).times(y))))
              .plus(ca.power(x, 2).times(7)).plus(ca.power(x, 2).times(3).times(y))},
          new MX[] {p.plus(x).plus(x.times(y)).div(cp), p.div(cp), x.div(cp), x.times(y).div(cp)},
          new MX[] {p.plus(x).plus(x.times(y)).div(cp.plus(x.times(7))), O, O, p.plus(x).plus(x.times(y)).div(cp.plus(x.times(7)))},
          new MX[] {p.plus(x).plus(x.times(y)).div(cp.plus(x.times(7)).plus(x.times(3).times(y))), O, O,
            p.plus(x).plus(x.times(y)).div(cp.plus(x.times(7)).plus(x.times(3).times(y)))});
      List<MX> ms = List.of(x, y, p);
      List<org.casadi.GenericMatrix> ts = twins(isMX, ms);
      separateLinearCases(isMX, ms, ts, conv(isMX, ms, ts, List.of(ca.vertcat(List.of(x, y)))).get(0), cases);

      x = MX.sym("x", 2);
      y = MX.sym("y", 2);
      p = MX.sym("p", 2);
      cp = ca.cos(p);
      MX xy = x.times(y);
      cases = List.of(
          new MX[] {ca.dot(p, x), O, ca.dot(p, x), O},
          new MX[] {ca.dot(p.plus(x).plus(xy), cp.plus(x.times(7)).plus(xy.times(3))), ca.dot(p, cp),
            ca.dot(p, x.times(7)).plus(ca.dot(x, cp)),
            ca.dot(p, xy.times(3)).plus(ca.dot(x, x.times(7).plus(xy.times(3)))).plus(ca.dot(xy, cp.plus(x.times(7)).plus(xy.times(3))))});
      ms = List.of(x, y, p);
      ts = twins(isMX, ms);
      separateLinearCases(isMX, ms, ts, conv(isMX, ms, ts, List.of(ca.vertcat(List.of(x, y)))).get(0), cases);

      x = MX.sym("x", 2, 2);
      y = MX.sym("y", 2, 2);
      p = MX.sym("p", 2, 2);
      cp = ca.cos(p);
      xy = x.times(y);
      cases = List.of(
          new MX[] {p.plus(x).plus(xy).T(), p.T(), x.T(), xy.T()},
          new MX[] {ca.mtimes(p, x), O, ca.mtimes(p, x), O},
          new MX[] {ca.mtimes(p.plus(x).plus(xy), cp.plus(x.times(7)).plus(xy.times(3))), ca.mtimes(p, cp),
            ca.mtimes(p, x.times(7)).plus(ca.mtimes(x, cp)),
            ca.mtimes(p, xy.times(3)).plus(ca.mtimes(x, x.times(7).plus(xy.times(3)))).plus(ca.mtimes(xy, cp.plus(x.times(7)).plus(xy.times(3))))});
      ms = List.of(x, y, p);
      ts = twins(isMX, ms);
      separateLinearCases(isMX, ms, ts, conv(isMX, ms, ts, List.of(ca.veccat(List.of(x, y)))).get(0), cases);
    }
  }

  @Test
  public void test_separate_linear_mat() {
    MX x = MX.sym("x", 3);
    MX y = MX.sym("y", 3, 3);
    MX z = MX.sym("z", 2);
    MX p = MX.sym("p");

    DM A = dm(new double[][] {{0, 1, 3}, {2, 0, 0}, {0, 9, 0}});
    DM Asp = ca.sparsify(A);

    DM B = dm(new double[][] {{1.7, 1.1}, {6, 0}, {0, 1.2}});
    DM Bsp = ca.sparsify(B);
    long[] off = {0, 2, 7};
    MX yx = ca.mtimes(y, x);
    MX Z3 = new MX(DM.zeros(3, 1));
    MX full = ca.vertcat(List.of(new MX(3), x.times(3), yx)).times(9);
    MX fullc = ca.vertcat(List.of(new MX(3), Z3, Z3)).times(9);
    MX fulll = ca.vertcat(List.of(new MX(0), x.times(3), Z3)).times(9);
    MX fulln = ca.vertcat(List.of(new MX(0), Z3, yx)).times(9);
    MX p2 = ca.vertsplit(full, off).get(1);
    MX p2c = ca.vertsplit(fullc, off).get(1);
    MX p2l = ca.vertsplit(fulll, off).get(1);
    MX p2n = ca.vertsplit(fulln, off).get(1);
    DM nzi = idx(2, 3, 0, 4, 6, 6);
    MX E = full.get_nz(false, nzi);
    MX Ec = fullc.get_nz(false, nzi);
    MX El = fulll.get_nz(false, nzi);
    MX En = fulln.get_nz(false, nzi);

    MX Z5 = new MX(DM.zeros(5, 1)), Z3b = new MX(DM.zeros(3, 1)), Z2 = new MX(DM.zeros(2, 1));
    MX[][] cases = {
      {ca.mtimes(A, x).plus(ca.mtimes(B, z)), Z3b, ca.mtimes(A, x), ca.mtimes(B, z)},
      {ca.mtimes(Asp, x).plus(ca.mtimes(Bsp, z)), Z3b, ca.mtimes(Asp, x), ca.mtimes(Bsp, z)},
      {ca.vertcat(List.of(x.times(3).plus(yx), z.times(9))), Z5, ca.vertcat(List.of(x.times(3), Z2)), ca.vertcat(List.of(yx, z.times(9)))},
      {ca.vertcat(List.of(x.times(3).plus(yx), z.times(9))).times(7), Z5, ca.vertcat(List.of(x.times(3), Z2)).times(7),
        ca.vertcat(List.of(yx, z.times(9))).times(7)},
      {p2.times(9), p2c.times(9), p2l.times(9), p2n.times(9)},
      {E, Ec, El, En}};

    for (MX[] c : cases) {
      System.out.println("separate_linear(" + c[0] + ")");
      System.out.println("ref");
      System.out.println("  const " + c[1]);
      System.out.println("  lin " + c[2]);
      System.out.println("  nonlin " + c[3]);

      List<Object> r = ca.separate_linear(c[0], List.of(x, y), List.of(p));
      MX expr_const = (MX) r.get(0), expr_lin = (MX) r.get(1), expr_nonlin = (MX) r.get(2);
      System.out.println("actual");
      System.out.println("  const " + expr_const);
      System.out.println("  lin " + expr_lin);
      System.out.println("  nonlin " + expr_nonlin);

      java.util.function.BiConsumer<MX, MX> test_equal = (a, b) -> {
        Function f1 = new Function("f", List.of(x, y, z, p), List.of(a));
        Function f2 = new Function("f", List.of(x, y, z, p), List.of(b));
        DM.rng(1);
        List<DM> args = new ArrayList<>();
        for (long i = 0; i < f1.n_in(); i++) args.add(DM.rand(f1.sparsity_in(i)));
        checkarray(f1.call(args).get(0), f2.call(args).get(0));
      };

      test_equal.accept(c[0], expr_const.plus(expr_lin).plus(expr_nonlin));

      test_equal.accept(expr_const, c[1]);
      test_equal.accept(expr_lin, c[2]);
      test_equal.accept(expr_nonlin, c[3]);

      assertTrue(ca.is_linear(expr_lin, ca.veccat(List.of(x, y))));
      assertFalse(ca.depends_on(expr_const, ca.veccat(List.of(x, y))));
    }
  }

  @Test
  public void test_extract_parametric_opts() {
    MX x = MX.sym("x"), y = MX.sym("y"), p = MX.sym("p");
    List<MX> ms = List.of(x, y, p);
    for (boolean isMX : new boolean[] {false, true}) {
      List<org.casadi.GenericMatrix> ts = twins(isMX, ms);
      org.casadi.GenericMatrix tp = ts.get(2);
      org.casadi.GenericMatrix expr = conv(isMX, ms, ts, List.of(x.times(2).times(p).plus(y.times(2).times(ca.power(p, 2))))).get(0);

      List<Object> r = extractParametric(expr, tp, Map.of("extract_trivial", true));
      List<org.casadi.GenericMatrix> symbols = gl(r.get(1));

      assertEquals(2, symbols.size());

      assertEquals("e_0", nameOf(symbols.get(0)));
      assertEquals("e_1", nameOf(symbols.get(1)));

      r = extractParametric(expr, tp, Map.of("offset", 5, "prefix", "foo_", "suffix", "bar", "extract_trivial", true));
      symbols = gl(r.get(1));

      assertEquals("foo_5bar", nameOf(symbols.get(0)));
      assertEquals("foo_6bar", nameOf(symbols.get(1)));

      System.out.println(symbols);
      System.out.println(r.get(2));

      r = extractParametric(expr, tp, Map.of("extract_trivial", false));
      assertEquals(1, gl(r.get(1)).size());
      assertTrue(String.valueOf(r.get(2)).contains("sq(p)"));
      System.out.println(r.get(2));

      r = extractParametric(expr, tp, Map.of("extract_trivial", true));
      assertEquals(2, gl(r.get(1)).size());
      assertTrue(String.valueOf(r.get(2)).contains("sq(p)"));
      System.out.println(r.get(2));
    }
  }

  @Test
  public void test_weakref() {
    MX x = MX.sym("x");

    org.casadi.WeakRef y = new org.casadi.WeakRef(x);

    assertTrue(y.alive());

    x.close();

    assertFalse(y.alive());

    x = MX.sym("x");
    y = new org.casadi.WeakRef(x);
    org.casadi.SharedObject xx = y.shared();
    xx.close();

    assertTrue(y.alive());
  }

  @Test
  public void test_nonzeros() {
    MX A = MX.sym("A", 3, 3);

    MX x = MX.sym("x");
    MX y = MX.sym("y");
    MX z = MX.sym("z");
    MX T = MX.sym("T", Sparsity.upper(3));
    MX w = MX.sym("w", 3, 3);

    DM.rng(1);

    for (MX expr : List.of(ca.vertcat(List.of(x, y, z)), ca.veccat(List.of(x, y, z, w)), ca.veccat(List.of(x, T, z, w)), x.times(2),
        T.times(4).times(x), ca.vertcat(List.of(x, ca.vertcat(List.of(y, z)))))) {
      Function fref = new Function("fref", List.of(x, y, z, w, A, T), List.of(expr.get_nz(false, new Slice())));
      Function f = new Function("f", List.of(x, y, z, w, A, T), List.of(ca.vertcat(expr.nonzeros())));

      List<DM> in = new ArrayList<>();
      for (long i = 0; i < f.n_in(); i++) in.add(DM.rand(f.sparsity_in(i)));
      checkfunctionLight(f, fref, in);
    }

    for (MX expr : List.of(ca.vertcat(List.of(x, y, z)), ca.vertcat(List.of(x, ca.vertcat(List.of(y, z)))))) {
      String s = String.valueOf(expr.nonzeros());
      // Python checks the list repr "[MX(x), MX(y), MX(z)]"
      assertTrue(s.contains("[x, y, z]"), s);
    }
  }

  @Test
  public void test_printme_codegen() {
    {
      SX x = SX.sym("x");
      Function f = new Function("f", List.of(x), List.of(ca.power(x, 2).printme(2).plus(6)));
      checkCodegen(f, dms(dm(3)));
    }
    {
      MX x = MX.sym("x");
      Function f = new Function("f", List.of(x), List.of(ca.power(x, 2).printme(2).plus(6)));
      checkCodegen(f, dms(dm(3)));
    }
  }

  private void jtimesEmptyMX() {
    for (long nx : new long[] {2, 1, 0}) {
      for (long mx : new long[] {3, 1, 0}) {
        MX x = MX.sym("x", nx, mx);
        for (long ny : new long[] {3, 1, 0}) {
          for (long my : new long[] {4, 1, 0}) {
            for (long bm : new long[] {1, 3, 0}) {
              if (bm != 1 && my == 0) continue;
              MX y = ca.mtimes(ca.mtimes(DM.rand(ny, nx), x), DM.rand(mx, my));
              assertEquals(ny, y.size1());
              assertEquals(my, y.size2());
              DM yb = DM.rand(ny, my * bm);

              MX r = ca.jtimes(y, x, new MX(yb), true);
              List<MX> parts = new ArrayList<>();
              for (DM e : ca.horzsplit(yb, Math.max(1, my)))
                parts.add(ca.reshape(ca.mtimes(ca.jacobian(y, ca.vec(x)).T(), new MX(ca.vec(e))), x.size1(), x.size2()));
              // Python hcat([]) is DM(1,0); ca.horzcat(empty list) throws ClassCastException
              MX r_ref = parts.isEmpty() ? new MX(1, 0) : ca.horzcat(parts);
              assertEquals(nx, r.size1());
              assertEquals(mx * bm, r.size2());
              Function f = new Function("f", List.of(x), List.of(r));
              Function f_ref = new Function("f", List.of(x), List.of(r_ref));
              checkfunctionLight(f, f_ref, dms(DM.rand(nx, mx)));
            }
          }
        }
      }
    }
    for (long nx : new long[] {2, 1, 0}) {
      for (long mx : new long[] {3, 1, 0}) {
        MX x = MX.sym("x", nx, mx);
        for (long ny : new long[] {3, 1, 0}) {
          for (long my : new long[] {4, 1, 0}) {
            for (long bm : new long[] {1, 3, 0}) {
              if (bm != 1 && mx == 0) continue;
              MX y = ca.mtimes(ca.mtimes(DM.rand(ny, nx), x), DM.rand(mx, my));
              assertEquals(ny, y.size1());
              assertEquals(my, y.size2());
              DM dx = DM.rand(nx, mx * bm);

              MX r = ca.jtimes(y, x, new MX(dx));
              List<MX> parts = new ArrayList<>();
              for (DM e : ca.horzsplit(dx, Math.max(1, mx)))
                parts.add(ca.reshape(ca.mtimes(ca.jacobian(y, ca.vec(x)), new MX(ca.vec(e))), y.size1(), y.size2()));
              // Python hcat([]) is DM(1,0); ca.horzcat(empty list) throws ClassCastException
              MX r_ref = parts.isEmpty() ? new MX(1, 0) : ca.horzcat(parts);

              assertEquals(ny, r.size1());
              assertEquals(my * bm, r.size2());
              Function f = new Function("f", List.of(x), List.of(r));
              Function f_ref = new Function("f", List.of(x), List.of(r_ref));
              checkfunctionLight(f, f_ref, dms(DM.rand(nx, mx)));
            }
          }
        }
      }
    }
  }

  private void powZeroMX() {
    MX x = MX.sym("x", Sparsity.diag(4));
    assertTrue(ca.power(x, 0).is_dense());
    assertFalse(ca.power(x, 0.3).is_dense());
    assertFalse(ca.power(x, 1).is_dense());
    assertFalse(ca.power(x, 2).is_dense());
    assertFalse(ca.power(x, 2.3).is_dense());
    assertTrue(ca.power(x, -2).is_dense());
    MX p = MX.sym("p");

    assertTrue(ca.power(x, p).is_dense());

    MX s = MX.sym("x");

    checkarray(ca.evalf(ca.substitute(ca.vertcat(List.of(ca.power(s, 2), ca.power(s, 1), ca.power(s, 0))), s, new MX(0))), dm(0, 0, 1));
    checkarray(ca.evalf(ca.substitute(ca.power(s, dm(2, 1, 0)), s, new MX(0))), dm(0, 0, 1));
  }

  private void jtimesEmptySX() {
    for (long nx : new long[] {2, 1, 0}) {
      for (long mx : new long[] {3, 1, 0}) {
        SX x = SX.sym("x", nx, mx);
        for (long ny : new long[] {3, 1, 0}) {
          for (long my : new long[] {4, 1, 0}) {
            for (long bm : new long[] {1, 3, 0}) {
              if (bm != 1 && my == 0) continue;
              SX y = ca.mtimes(ca.mtimes(DM.rand(ny, nx), x), DM.rand(mx, my));
              assertEquals(ny, y.size1());
              assertEquals(my, y.size2());
              DM yb = DM.rand(ny, my * bm);

              SX r = ca.jtimes(y, x, new SX(yb), true);
              List<SX> parts = new ArrayList<>();
              for (DM e : ca.horzsplit(yb, Math.max(1, my)))
                parts.add(ca.reshape(ca.mtimes(ca.jacobian(y, ca.vec(x)).T(), new SX(ca.vec(e))), x.size1(), x.size2()));
              // Python hcat([]) is DM(1,0); ca.horzcat(empty list) throws ClassCastException
              SX r_ref = parts.isEmpty() ? new SX(1, 0) : ca.horzcat(parts);
              assertEquals(nx, r.size1());
              assertEquals(mx * bm, r.size2());
              Function f = new Function("f", List.of(x), List.of(r));
              Function f_ref = new Function("f", List.of(x), List.of(r_ref));
              checkfunctionLight(f, f_ref, dms(DM.rand(nx, mx)));
            }
          }
        }
      }
    }
    for (long nx : new long[] {2, 1, 0}) {
      for (long mx : new long[] {3, 1, 0}) {
        SX x = SX.sym("x", nx, mx);
        for (long ny : new long[] {3, 1, 0}) {
          for (long my : new long[] {4, 1, 0}) {
            for (long bm : new long[] {1, 3, 0}) {
              if (bm != 1 && mx == 0) continue;
              SX y = ca.mtimes(ca.mtimes(DM.rand(ny, nx), x), DM.rand(mx, my));
              assertEquals(ny, y.size1());
              assertEquals(my, y.size2());
              DM dx = DM.rand(nx, mx * bm);

              SX r = ca.jtimes(y, x, new SX(dx));
              List<SX> parts = new ArrayList<>();
              for (DM e : ca.horzsplit(dx, Math.max(1, mx)))
                parts.add(ca.reshape(ca.mtimes(ca.jacobian(y, ca.vec(x)), new SX(ca.vec(e))), y.size1(), y.size2()));
              // Python hcat([]) is DM(1,0); ca.horzcat(empty list) throws ClassCastException
              SX r_ref = parts.isEmpty() ? new SX(1, 0) : ca.horzcat(parts);

              assertEquals(ny, r.size1());
              assertEquals(my * bm, r.size2());
              Function f = new Function("f", List.of(x), List.of(r));
              Function f_ref = new Function("f", List.of(x), List.of(r_ref));
              checkfunctionLight(f, f_ref, dms(DM.rand(nx, mx)));
            }
          }
        }
      }
    }
  }

  private void powZeroSX() {
    SX x = SX.sym("x", Sparsity.diag(4));
    assertTrue(ca.power(x, 0).is_dense());
    assertFalse(ca.power(x, 0.3).is_dense());
    assertFalse(ca.power(x, 1).is_dense());
    assertFalse(ca.power(x, 2).is_dense());
    assertFalse(ca.power(x, 2.3).is_dense());
    assertTrue(ca.power(x, -2).is_dense());
    SX p = SX.sym("p");

    assertTrue(ca.power(x, p).is_dense());

    SX s = SX.sym("x");

    checkarray(ca.evalf(ca.substitute(ca.vertcat(List.of(ca.power(s, 2), ca.power(s, 1), ca.power(s, 0))), s, new SX(0))), dm(0, 0, 1));
    checkarray(ca.evalf(ca.substitute(ca.power(s, dm(2, 1, 0)), s, new SX(0))), dm(0, 0, 1));
  }

  @Test
  public void test_norms() {
    for (int k = 0; k < 2; k++) {
      final int kk = k;
      UnaryOperator<SX> fs = x -> kk == 0 ? ca.vertcat(List.of(x.get(0), new SX(1, 1), x.get(new Slice(1, 5)))) : ca.diag(x);
      UnaryOperator<MX> fm = x -> kk == 0 ? ca.vertcat(List.of(x.get(0), new MX(1, 1), x.get(new Slice(1, 5)))) : ca.diag(x);
      for (int n = 0; n < 3; n++) {
        System.out.println("norm " + n);
        SX xs = SX.sym("x", 5);
        SX es = fs.apply(ca.sin(xs));
        es = ca.log(n == 0 ? ca.norm_fro(es) : n == 1 ? ca.norm_1(es) : ca.norm_inf(es));
        Function fsx = new Function("f", List.of(xs), List.of(es));

        MX xm = MX.sym("x", 5);
        MX em = fm.apply(ca.sin(xm));
        em = ca.log(n == 0 ? ca.norm_fro(em) : n == 1 ? ca.norm_1(em) : ca.norm_inf(em));
        Function fmx = new Function("f", List.of(xm), List.of(em));
        for (DM args : List.of(dm(-0.3, 0.3, 0.21, 0.17, 0), dm(-1, -2, -7, -7, -7))) checkfunction(fsx, fmx, dms(args));
        checkCodegen(fmx, dms(dm(-0.3, 0.3, 0.21, 0.17, 0)));
      }
    }
  }

  @Test
  public void test_contains() {
    MX x = MX.sym("x", 2, 2);
    MX y = MX.sym("y");
    MX z = MX.sym("z", 3, 3);

    MX e = x.times(y);
    assertTrue(ca.contains(List.of(x, y, z), x));
    assertFalse(ca.contains(List.of(x, y), z));
    assertTrue(ca.contains_any(List.of(x, y), List.of(y, z)));
    assertFalse(ca.contains_all(List.of(x, y), List.of(y, z)));
    assertTrue(ca.contains_any(List.of(x, y), List.of(x, y)));
    assertTrue(ca.contains_all(List.of(x, y), List.of(x, y)));
    assertTrue(ca.contains(List.of(e, x), e));
  }

  @Test
  public void test_issue4041() {
    MX x = MX.sym("x");
    MX y = MX.sym("y");

    Function f = new Function("f", List.of(x, y), List.of(ca.power(x, 2)), Map.of("never_inline", true));

    MX expr = (MX) f.invoke(x, y);

    assertFalse(ca.depends_on(expr, y));
    assertTrue(ca.contains(ca.symvar(expr), y));

    f = new Function("f", List.of(x, y), List.of(ca.power(x, 2)), Map.of("always_inline", true));

    expr = (MX) f.invoke(x, y);

    assertFalse(ca.depends_on(expr, y));
    assertFalse(ca.contains(ca.symvar(expr), y));
  }

  @Test
  public void test_symmetric_jacobian() {
    MX x = MX.sym("x", 4);
    MX u = MX.sym("u", 2);

    MX h = x.get(1).times(5).times(x.get(2));

    long nh = 1;
    MX lam_h = MX.sym("lam_h", nh, 1);

    MX ux = ca.vertcat(List.of(u, x));
    MX adj_ux = ca.densify(ca.jtimes(h, ux, lam_h, true));
    ca.jacobian(adj_ux, ux, Map.of("symmetric", true));

    MX adj_ux2 = ca.jtimes(h, ux, lam_h, true);
    assertInException("Symmetry exploitation", () -> ca.jacobian(adj_ux2, ux, Map.of("symmetric", true)));
  }

  @Test
  public void test_pow_zero() {
    powZeroSX();
    powZeroMX();
  }

  @Test
  public void test_pow() {
    MX x = MX.sym("x");
    DM r = ca.evalf(ca.substitute(ca.power(x, dm(2, 1, 0)), x, new MX(0)));
    checkarray(r, dm(0, 0, 1));

    MX y = new MX(4, 1);
    assertEquals(4, ca.power(y, 0).nnz());
    y = new MX(4, 1);
    y.set(x, false, 1);
    assertEquals(4, ca.power(y, 0).nnz());
    y = new MX(4, 1);
    y.set(0, false, 1);
    assertEquals(4, ca.power(y, 0).nnz());
    y = new MX(4, 1);
    y.set(1, false, 1);
    assertEquals(4, ca.power(y, 0).nnz());
  }

  // test_sum: python-only (numpy np.sum dispatch on casadi types)

  @Test
  public void test_jtimes_empty() {
    assertInException("Ambiguous", () -> ca.jtimes(new SX(5), SX.sym("x", 2, 0), new SX(DM.ones(2, 5))));
    assertInException("Ambiguous", () -> ca.jtimes(new SX(DM.ones(4, 0)), SX.sym("x", 2, 3), new SX(DM.ones(4, 6)), true));

    jtimesEmptySX();
    jtimesEmptyMX();
  }

  @Test
  public void test_issue2934() {
    MX M = MX.sym("X", 5, 5);
    MX Y = MX.sym("y", 2, 2);

    MX m = new MX(M);

    m.set(Y, false, idx(0, 0), idx(0, 0));

    System.out.println(m);
    Function f = new Function("f", List.of(M, Y), List.of(m));
    System.out.println(f.call(List.of(M, Y), true, false).get(0));
    assertTrue(f.call(List.of(M.times(2), Y), true, false).get(0).toString().contains("((2.*X)[-1, -1, -1, 0] = y)"));
  }

  @Test
  public void test_issue2935() {
    MX M = new MX(DM.rand(5, 5));
    MX m = M;

    MX Y = MX.sym("y", 2, 2);
    Sparsity E = Y.sparsity();
    MX y = Y;

    m.set(y, false, idx(0, 0), idx(0, 0));

    MX e = ca.cos(m);
    e = ca.dot(e, e);

    MX eb = MX.sym("eb");
    MX ge = ca.jtimes(e, Y, eb, true);

    Function f = new Function("f", List.of(eb, Y), List.of(ge));

    MX e1 = f.call(List.of(eb, Y), true, false).get(0);
    MX e2 = f.call(List.of(eb, Y.plus(1e-300)), true, false).get(0);

    System.out.println(e1);
    System.out.println(e2);

    Function f1 = new Function("f1", List.of(eb, Y), List.of(ge));
    Function f2 = new Function("f2", List.of(eb, Y), List.of(e2));

    List<DM> args = dms(dm(1), DM.rand(E.size1(), E.size2()));

    System.out.println(f1.call(args));
    System.out.println(f2.call(args));

    assertTrue(f1.call(args).get(0).sparsity().is_equal(f2.call(args).get(0).sparsity()));
  }

  @SuppressWarnings("unchecked")
  @Test
  public void test_extract() {
    MX x = MX.sym("x", 2);
    MX y = MX.sym("y", 2);

    MX u = ca.vertcat(List.of(x, y));
    DM.rng(1);
    DM u0 = DM.rand(4, 1);
    Function f = new Function("f", List.of(x, y), List.of(ca.sin(x.plus(y))), Map.of("never_inline", true));

    Function g = new Function("g", List.of(x), List.of(ca.cos(ca.sum(x).times(3))), Map.of("never_inline", true));

    MX e = ((MX) g.invoke(((MX) f.invoke(x.div(y), x.times(2).times(y))).times(x))).div(x.plus(y));

    List<Object> res = ca.extract(List.of(e), Map.of("lift_shared", false, "lift_calls", true));

    List<MX> vexpr = (List<MX>) res.get(0), v = (List<MX>) res.get(1), vdef = (List<MX>) res.get(2);

    System.out.println(vexpr);
    System.out.println(v);
    System.out.println(vdef);

    assertEquals(5, v.size());
    assertEquals(v.size(), vdef.size());
    assertEquals(1, vexpr.size());
    Function f_ref = new Function("f", List.of(u), List.of(e));
    System.out.println("e " + e);
    List<MX> e_reconstructed = (List<MX>) ca.substitute_inplace(v, vdef, List.of(e)).get(1);
    System.out.println(e_reconstructed);
    Function fr = new Function("f", List.of(u), e_reconstructed);
    checkarray((DM) f_ref.invoke(u0), (DM) fr.invoke(u0));

    System.out.println(e);

    System.out.println(ca.substitute_inplace(v, vdef, vexpr));

    MX vexprc = ca.vertcat(vexpr);
    MX vc = ca.vertcat(v);
    MX vdefc = ca.vertcat(vdef);

    Function Jf_ref = new Function("J_ref", List.of(u), List.of(ca.jacobian(e, u)));

    MX J_v_u = ca.solve(new MX(DM.eye(vdefc.size1())).minus(ca.jacobian(vdefc, vc)), ca.jacobian(vdefc, u));
    MX J = ca.jacobian(vexprc, u).plus(ca.mtimes(ca.jacobian(vexprc, vc), J_v_u));

    vexpr = (List<MX>) res.get(0);
    v = (List<MX>) res.get(1);
    vdef = (List<MX>) res.get(2);

    J = ((List<MX>) ca.substitute_inplace(v, vdef, List.of(J)).get(1)).get(0);

    Function Jf = new Function("J2", List.of(u), List.of(J));

    checkarray((DM) Jf_ref.invoke(u0), (DM) Jf.invoke(u0));
  }

  @Test
  public void test_const_folding_on_the_fly() {
    MX x = MX.sym("x");

    for (int a : new int[] {2, 7}) {
      for (int b : new int[] {2, 7}) {
        MX ref = ca.times(a * b, x);
        if (a == 7 && b == 7) {
          assertEquals(ref.toString(), ca.times(a, ca.times(b, x)).toString());
          assertEquals(ref.toString(), ca.times(ca.times(b, x), a).toString());
        } else {
          assertFalse(ref.toString().equals(ca.times(a, ca.times(b, x)).toString()));
          assertFalse(ref.toString().equals(ca.times(ca.times(b, x), a).toString()));
        }
      }
    }
  }

  @Test
  public void test_cse_reuse() {
    for (int X = 0; X < 2; X++) {
      MX x = MX.sym("x");
      MX y = ca.sin(MX.sym("y"));
      MX z = MX.sym("z");

      MX e = x.plus(y.div(z));
      MX r = ca.vertcat(List.of(ca.cse(e), y.times(2)));
      System.out.println(r);
      assertTrue(r.toString().contains("@1=sin(y)"));
    }
  }

  @Test
  public void test_is_slice2() {
    int m = 13;
    int n = 11;

    MX x = MX.sym("x", m, n);
    DM.rng(1);
    DM A = DM.rand(m, n);

    for (int step_col : new int[] {1, 2, 3}) {
      for (int step_row : new int[] {1, 2, 3}) {
        List<Long> nz = new ArrayList<>();
        List<Double> cslice = new ArrayList<>(), rslice = new ArrayList<>();
        for (int col = 3; col < 7; col += step_col) cslice.add((double) col);
        for (int row = 1; row < 9; row += step_row) rslice.add((double) row);
        for (double col : cslice)
          for (double row : rslice) nz.add((long) (row + m * col));

        assertTrue(ca.is_slice2(toLongs(nz)));

        for (int j = 0; j < nz.size(); j++) {
          List<Long> nz_mod = new ArrayList<>(nz);
          long val = j == 0 ? nz.get(0) - 1 : (long) Math.floor((nz.get(j) + nz.get(j - 1)) / 2.0);
          nz_mod.add(j, val);
          assertFalse(ca.is_slice2(toLongs(nz_mod)));
        }
        DM ri = new DM(rslice.stream().mapToDouble(Double::doubleValue).toArray());
        DM ci = new DM(cslice.stream().mapToDouble(Double::doubleValue).toArray());
        MX r = x.get(ri, ci);
        assertTrue(r.toString().contains(";"));
        Function f = new Function("f", List.of(x), List.of(r));
        checkarray((DM) f.invoke(A), A.get(ri, ci));
      }
    }
  }

  @Test
  public void test_is_diff_sparsity_propagation() {
    for (boolean isMX : new boolean[] {true, false}) {
      boolean[] is_diff_in = {true, true, false};
      boolean[] is_diff_out = {true, false};
      Map<String, Object> o = Map.of("never_inline", true, "is_diff_in", is_diff_in, "is_diff_out", is_diff_out);

      Function F, F2;
      if (isMX) {
        MX x = MX.sym("x"), u = MX.sym("u"), q = MX.sym("q");
        F = new Function("F", List.of(x, u, q), List.of(ca.sin(u.plus(x).plus(q)), u.times(q)), o);
        F2 = new Function("F", List.of(x, u, q), List.of(ca.sin(u.plus(x).plus(q)), u.times(q)));
      } else {
        SX x = SX.sym("x"), u = SX.sym("u"), q = SX.sym("q");
        F = new Function("F", List.of(x, u, q), List.of(ca.sin(u.plus(x).plus(q)), u.times(q)), o);
        F2 = new Function("F", List.of(x, u, q), List.of(ca.sin(u.plus(x).plus(q)), u.times(q)));
      }
      List<Function> Fs = new ArrayList<>(List.of(F));
      if (RUN_SLOW) {
        List<DM> in = dms(dm(1), dm(2), dm(3));
        checkCodegen(F, in);
        Fs.add(lastExternal);

        checkCodegen(F2, in, Map.of(), 1e-15, "c89", false, true);
        Fs.add(ca.external("F", lastLib, Map.of("is_diff_in", is_diff_in, "is_diff_out", is_diff_out)));

        checkCodegen(F2, in);
        Fs.add(ca.external("F", lastLib, Map.of("is_diff_in", is_diff_in, "is_diff_out", is_diff_out)));
      }

      for (Function Fi : Fs) {
        for (boolean isMX2 : new boolean[] {true, false}) {
          for (int ad_weight_sp : new int[] {0, 1}) {
            Function H;
            if (isMX2) {
              MX X0 = MX.sym("X0"), U0 = MX.sym("U0"), Q0 = MX.sym("Q0");
              List<MX> r = Fi.call(List.of(X0.times(5), U0.times(6), Q0.times(7)));
              H = new Function("H", List.of(ca.vertcat(List.of(X0, U0, Q0))), List.of(r.get(0).times(2), r.get(1).times(3)),
                  Map.of("ad_weight_sp", ad_weight_sp));
            } else {
              SX X0 = SX.sym("X0"), U0 = SX.sym("U0"), Q0 = SX.sym("Q0");
              List<SX> r = Fi.call(List.of(X0.times(5), U0.times(6), Q0.times(7)));
              H = new Function("H", List.of(ca.vertcat(List.of(X0, U0, Q0))), List.of(r.get(0).times(2), r.get(1).times(3)),
                  Map.of("ad_weight_sp", ad_weight_sp));
            }
            assertEquals(2, H.sparsity_jac(0, 0).nnz());
            assertEquals(0, H.sparsity_jac(0, 1).nnz());
          }
        }
      }
    }
  }

  @Test
  public void test_is_nonnegative_binary() {
    MX x = MX.sym("x");
    MX y = MX.sym("y");
    for (MX e : List.of(ca.lt(x, y), ca.le(x, y), ca.eq(x, y), ca.ne(x, y), ca.logic_and(x, y), ca.logic_or(x, y), ca.hypot(x, y), ca.sqrt(x)))
      assertTrue(ca.is_equal(ca.abs(e), e));
    for (MX e : List.of(x.plus(y), x.times(y), x.minus(y))) assertFalse(ca.is_equal(ca.abs(e), e));
  }

  @Test
  public void test_chained_inequality_kept() {
    MX y = MX.sym("y");
    for (MX e : List.of(ca.le(0, ca.le(y, 1)), ca.ge(ca.ge(1, y), 0))) {
      assertFalse(e.is_constant());
      Function f = new Function("f", List.of(y), List.of(e)).transform(List.of(List.of("simplify", 0, "ref_count")));
      assertTrue(f.call(List.of(y), true, false).get(0).is_one());
    }
  }

  @Test
  public void test_linear_slice_assign() {
    MX Y0 = MX.sym("Y", 3, 3);
    MX Y = new MX(Y0);
    Y.set(1, false, new Slice());
    checkarray((DM) new Function("f", List.of(Y0), List.of(Y)).invoke(DM.zeros(3, 3)), DM.ones(3, 3));
    MX Z0 = MX.sym("Z", 3, 3);
    MX Z = new MX(Z0);
    Z.set(1, false, 4);
    Function f = new Function("f", List.of(Z0), List.of(Z));
    checkarray((DM) f.invoke(DM.zeros(3, 3)), dm(new double[][] {{0, 0, 0}, {0, 1, 0}, {0, 0, 0}}));
  }

  @Test
  public void test_bspline_checks() {
    MX x = MX.sym("x");
    double[][] knots = {{0, 0, 0, 0, 1, 1, 1, 1}};
    assertInException("Expected 4 coefficients, but got 2", () -> ca.bspline(x, dm(1, 2), knots, new long[] {3}, 1, Map.of()));
    assertInException("Unknown option 'foo'", () -> ca.bspline(x, dm(1, 2, 3, 4), knots, new long[] {3}, 1, Map.of("foo", 1)));
  }

  @Test
  public void test_user_errors() {
    MX x = MX.sym("x", 2);
    MX y = MX.sym("y");
    assertInException("Expected a column vector with 3 nonzeros", () -> new MX(Sparsity.lower(2), dm(1, 2)));
    assertInException("Expected a primitive of dimension", () -> ca.vertcat(List.of(x, y)).join_primitives(List.of(DM.ones(3, 1), new DM(1))));
    assertInException("Expected 2 arguments, but got 1", () -> x.get(0).times(y).eval_mx(List.of(y)));
    assertInException("Not a binary operation", () -> MX.binary(999, y, y));
    MX A = MX.sym("A", 4);
    assertInException("Length of B must match dim_b",
        () -> MX.einstein(A, A, new long[] {2, 2}, new long[] {3, 2}, new long[] {2}, new long[] {-1, -2}, new long[] {-2, -3}, new long[] {-1}));
    assertInException("A must be a dense vector", () -> MX.einstein(MX.sym("A", 2, 2), A, new long[] {2, 2}, new long[] {2, 2},
        new long[] {2}, new long[] {-1, -2}, new long[] {-2, -3}, new long[] {-1}));
    assertInException("Expected 2 query vectors, but got 1",
        () -> MX.interpn_linear(List.of(new MX(dm(0, 1)), new MX(dm(0, 1))), MX.sym("v", 4), List.of(y)));
    assertInException("Grid must be a dense vector", () -> ca.low(MX.sym("v", Sparsity.lower(2)), 0.5));
    assertInException("Expected 2 row indices, but got 3", () -> x.enlarge(4, 1, new long[] {0, 1, 2}, new long[] {0}));
    double[] k = {0, 0, 0, 0, 1, 1, 1, 1};
    assertInException("multiple of the number of dimensions", () -> MX.bspline_dual(new double[] {0.5, 0.5, 0.1}, new double[][] {k, k}, new long[] {3, 3}));
  }

  @Test
  public void test_unary_op_check() {
    MX y = MX.sym("y");
    for (long op : new long[] {org.casadi.Operation.OP_ADD.swigValue(), 999})
      assertInException("Not a unary operation", () -> MX.unary(op, y));
  }

  @Test
  public void test_print_operator_args() {
    MX x = MX.sym("x");
    MX y = MX.sym("y");
    assertInException("Expected 2 argument names, but got 1", () -> ca.print_operator(x.plus(y), List.of("a")));
    assertEquals("sin(a)", ca.print_operator(ca.sin(x), List.of("a", "b")));
  }

  @Test
  public void test_set_precision() {
    double pi_val = 3.141592653589793239;
    MX mx_scalar = new MX(pi_val);
    MX mx_dm = new MX(dm(1.234567890123, 9.876543210987));
    try {
      assertEquals("3.14159", mx_scalar.toString());
      DM.set_precision(12);
      assertEquals("3.14159265359", mx_scalar.toString());
      assertEquals("[1.23456789012, 9.87654321099]", mx_dm.toString());
      DM.set_precision(4);
      DM.set_scientific(true);
      assertTrue(mx_scalar.toString().toLowerCase().contains("e"));
    } finally {
      DM.set_precision(6);
      DM.set_scientific(false);
    }
  }
}

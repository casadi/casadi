package org.casadi.test;

import static org.junit.jupiter.api.Assertions.assertArrayEquals;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertNotSame;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.StandardCopyOption;
import java.util.Arrays;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import org.casadi.CasadiMeta;
import org.casadi.DM;
import org.casadi.DeserializerBase;
import org.casadi.FileDeserializer;
import org.casadi.FileSerializer;
import org.casadi.Function;
import org.casadi.GenericMatrix;
import org.casadi.Linsol;
import org.casadi.MX;
import org.casadi.SX;
import org.casadi.SerializerBase;
import org.casadi.Sparsity;
import org.casadi.StringDeserializer;
import org.casadi.StringSerializer;
import org.casadi.ca;
import org.junit.jupiter.api.Assumptions;
import org.junit.jupiter.api.Test;

/** Port of test/python/misc.py */
public class MiscTest extends CasadiTestCase {


  /** helpers.checkfunction_light: outputs agree */
  private static void checkfunctionLight(Function f, Function g, List<DM> inputs) {
    List<DM> a = f.call(inputs), b = g.call(inputs);
    assertEquals(b.size(), a.size());
    for (int i = 0; i < a.size(); i++) checkarray(a.get(i), b.get(i), 1e-9, "output " + i);
  }

  private static Map<String, DM> call0(Function f) {
    return f.call(new HashMap<String, DM>());
  }

  private static void assertOutput(List<String> included, List<String> excluded, Runnable r) {
    String out = captureStdout(r);
    for (String s : included) assertTrue(out.contains(s), "expected '" + s + "' in output:\n" + out);
    for (String s : excluded) assertFalse(out.contains(s), "unexpected '" + s + "' in output:\n" + out);
  }

  private static String captureStderr(Runnable r) {
    PrintStream old = System.err;
    ByteArrayOutputStream b = new ByteArrayOutputStream();
    System.setErr(new PrintStream(b, true, StandardCharsets.UTF_8));
    try {
      r.run();
    } finally {
      System.err.flush();
      System.setErr(old);
    }
    return b.toString(StandardCharsets.UTF_8);
  }

  @Test
  public void test_issue179B() {
    message("Regression test #179 (B)");
    SX x = SX.sym("x");
    Function f = new Function("f", List.of(x), List.of(ca.power(x, 2.0)));
    Sparsity sparsity = f.jac_sparsity(0, 0);
    sparsity.toString();
  }

  // test_issue4216: python-only (complex numbers / numpy complex arrays)

  @Test
  public void test_sanity() {
    double[] v = {0.738, 0.39, 0.99};
    new DM(new Sparsity(4, 3, new long[] {0, 2, 2, 3}, new long[] {1, 2, 1}), new DM(v));
    assertThrows(RuntimeException.class, () -> new DM(new Sparsity(4, 4, new long[] {0, 2, 2, 3}, new long[] {1, 2, 1}), new DM(v)));
    assertThrows(RuntimeException.class, () -> new DM(new Sparsity(4, 3, new long[] {0, 2, 2, 12}, new long[] {1, 2, 1}), new DM(v)));
    assertThrows(RuntimeException.class, () -> new DM(new Sparsity(4, 3, new long[] {-10, 2, 2, 3}, new long[] {1, 2, 1}), new DM(v)));
    assertThrows(RuntimeException.class, () -> new DM(new Sparsity(4, 3, new long[] {0, 2, 2, 3}, new long[] {8, 2, 1}), new DM(v)));
    assertThrows(RuntimeException.class, () -> new DM(new Sparsity(4, 3, new long[] {0, 2, 2, 3}, new long[] {-3, 2, 1}), new DM(v)));
    assertThrows(RuntimeException.class, () -> new DM(new Sparsity(4, 3, new long[] {0, 2, 2, 3}, new long[] {1, 2, 1, 2}), new DM(v)));
    assertThrows(RuntimeException.class, () -> new DM(new Sparsity(4, 3, new long[] {0, 2, 0, 3}, new long[] {1, 2, 1}), new DM(v)));
  }

  @Test
  public void test_copyconstr_norefcount() {
    message("Copy constructor for non-refcounted classes");
    DM x = DM.ones(2, 3);
    DM y = new DM(x);
    x.set(0, 0, 5);
    assertNotSame(x, y);
    assertEquals(5, x.get(0, 0).to_double());
    assertEquals(1, y.get(0, 0).to_double());
  }

  @Test
  public void test_copyconstr_refcount() {
    message("Copy constructor for refcounted classes");
    Sparsity x = Sparsity.diag(4);
    Sparsity y = new Sparsity(x);
    x.resize(2, 8);
    assertNotSame(x, y);
    assertTrue(x.numel() != 0 && y.numel() != 0);
    assertArrayEquals(new long[] {2, 8}, x.size(), "shape");
    assertArrayEquals(new long[] {4, 4}, y.size(), "shape");
  }

  // test_copy_norefcount: python-only (copy.copy protocol)
  // test_copy_refcount: python-only (copy.copy protocol)
  // test_deepcopy_norefcount: python-only (copy.deepcopy protocol)
  // test_deepcopy_refcount: python-only (copy.deepcopy protocol)
  // test_options_introspection: not ported (optionNames/optionDescription/optionTypeName exist in neither binding)
  // test_pickling: python-only (pickle protocol)
  // test_exceptions: python-only except the mixed vertcat (overload-resolution errors are compile errors in Java)

  @Test
  public void test_exceptions() {
    List<GenericMatrix> v = Arrays.asList(new DM(1), SX.sym("x"), MX.sym("x"));
    assertInException("vertcat(", () -> ca.vertcat(v));
  }

  @Test
  public void test_getscheme() {
    SX x = SX.sym("x"), p = SX.sym("p");
    Function F = new Function("F", List.of(x, p), List.of(x.plus(p), ca.power(x, 2.0)), List.of("x", "p"), List.of("f", "g"));
    Map<String, DM> fc = F.call(Map.of("x", new DM(3), "p", new DM(4)));
    DM f = fc.get("f");
    checkarray(f, dm(7));
    DM g = fc.get("g");
    checkarray(g, dm(9));
  }

  @Test
  public void test_assertions() {
    MX x = MX.sym("x");
    MX z = ca.power(x, 2.0);
    z = z.attachAssert(ca.gt(z, 3.0), "x must be larger than 3");
    MX v = ca.sin(z);
    Function f = new Function("f", List.of(x), List.of(v));
    System.out.println(f);
    f.invoke(-6.0);
    try {
      f.invoke(1.0);
    } catch (RuntimeException e) {
      System.out.println(e.getMessage());
      assertTrue(e.getMessage().contains("x must be larger than 3"));
    }
  }

  // test_doc: python-only (__doc__ introspection; docs live in Javadoc)

  @Test
  public void test_output() {
    Assumptions.assumeTrue(hasNlpsol("ipopt"));
    String result = captureStdout(() -> dm(1, 2).print_dense());
    assertTrue(result.contains("2"));
    SX x = SX.sym("x");
    Map<String, SX> f = Map.of("x", x, "f", ca.power(x, 2.0));
    Function solver = ca.nlpsol("solver", "ipopt", f);
    result = captureStdout(() -> solver.call(Map.of("x0", new DM(0))));
    assertTrue(result.contains("Number of nonzeros in equality constraint"), result);
    assertTrue(result.contains("iter    objective    inf_pr"), result);
    String err = captureStderr(() -> {
      try {
        ca.nlpsol("solver", "foo", f);
      } catch (RuntimeException e) {
        // ignored as in Python
      }
    });
    assertTrue(err.contains("casadi_nlpsol_foo"), err);
  }

  @Test
  public void test_serialize() throws IOException {
    Sparsity x = Sparsity.upper(3);
    SX y = SX.sym("y");
    StringSerializer si = new StringSerializer();
    si.pack(List.of(x));
    SX z = y;
    for (int i = 0; i < 10000; i++) z = ca.sin(z);
    SX e = ca.vertcat(List.of(y, z, z.times(2.0)));
    Function fref = new Function("f", List.of(y), List.of(e));
    si.pack(e);
    String data = si.encode();
    si.pack(List.of(x));

    StringDeserializer sd = new StringDeserializer(data);
    List<?> spx = (List<?>) sd.unpack();
    assertTrue(((Sparsity) spx.get(0)).is_equal(x));
    SX A = (SX) sd.unpack();
    Function f = new Function("f", List.of(A.get(0)), List.of(A));
    checkfunctionLight(f, fref, List.of(new DM(7)));
    assertInException("end of stream", () -> sd.unpack());
    assertInException("end of stream", () -> sd.unpack());

    FileSerializer fs = new FileSerializer("foo.dat");
    fs.pack(List.of(x));
    z = ca.sin(y);
    e = ca.vertcat(List.of(y, z, z.times(2.0)));
    fref = new Function("f", List.of(y), List.of(e));
    fs.pack(e);
    fs.close();
    FileDeserializer fd = new FileDeserializer("foo.dat");
    spx = (List<?>) fd.unpack();
    assertTrue(((Sparsity) spx.get(0)).is_equal(x));
    A = (SX) fd.unpack();
    f = new Function("f", List.of(A.get(0)), List.of(A));
    checkfunctionLight(f, fref, List.of(new DM(7)));
    assertInException("end of stream", () -> fd.unpack());
    assertInException("end of stream", () -> fd.unpack());
    fd.close();

    SX xs = SX.sym("x");
    StringSerializer s = new StringSerializer();
    s.pack(xs);
    String data1 = s.encode();
    s.pack(ca.sin(xs));
    String data2 = s.encode();
    StringDeserializer s2 = new StringDeserializer(data1);
    s2.unpack();
    assertInException("end of stream", () -> s2.unpack());
    s2.decode(data2);
    s2.unpack();
    assertInException("end of stream", () -> s2.unpack());

    fs = new FileSerializer("foo.dat");
    Function f1 = new Function("f", List.of(xs), List.of(ca.power(xs, 2.0)));
    Function g1 = new Function("g", List.of(xs), List.of(ca.power(xs, 3.0)));
    fs.pack(List.of(f1, g1));
    fs.close();
    assertInException("File is not loadable with 'load'. Use 'FileDeserializer' instead.", () -> Function.load("foo.dat"));
    FileDeserializer fd2 = new FileDeserializer("foo.dat");
    System.out.println(fd2.unpack());
    fd2.close();

    fs = new FileSerializer("foo.dat");
    fs.pack(f1);
    fs.close();
    Function r = Function.load("foo.dat");
    System.out.println(r);

    f1.save("foo.dat");
    FileDeserializer fd3 = new FileDeserializer("foo.dat");
    System.out.println(fd3.unpack());
    fd3.close();
  }

  @Test
  public void test_print_time() {
    MX x = MX.sym("x");
    MX f = ca.power(x, 2.0);
    for (boolean print_time : new boolean[] {true, false}) {
      List<String> included = print_time ? List.of("t_wall") : List.of();
      List<String> excluded = print_time ? List.of() : List.of("t_wall");
      Map<String, Object> opts = new LinkedHashMap<>();
      opts.put("print_time", print_time);
      Function ff = new Function("f", List.of(x), List.of(f), opts);
      assertOutput(included, excluded, () -> ff.invoke(3.0));

      if (hasNlpsol("ipopt")) {
        Function solver = ca.nlpsol("nlpsol", "ipopt", Map.of("x", x, "f", f), opts);
        assertOutput(included, excluded, () -> call0(solver));
      }

      if (hasConic("ipopt")) {
        Function solver = ca.qpsol("qpsol", "qpoases", Map.of("x", x, "f", f), opts);
        assertOutput(included, excluded, () -> call0(solver));
      }

      Function rf = ca.rootfinder("rootfinder", "newton", Map.of("x", x, "g", x), opts);
      assertOutput(included, excluded, () -> call0(rf));

      Function in1 = ca.integrator("integrator", "rk", Map.of("x", x, "ode", f), opts);
      assertOutput(included, excluded, () -> call0(in1));

      Map<String, Object> integr_options = new LinkedHashMap<>(opts);
      integr_options.put("simplify", true);
      Function in2 = ca.integrator("integrator", "rk", Map.of("x", x, "ode", f), integr_options);
      assertOutput(included, excluded, () -> call0(in2));

      DM A = DM.rand(3, 3);
      DM b = dm(1, 2, 3);
      if (hasLinsol("lapacklu")) assertOutput(included, excluded, () -> ca.solve(A, b, "lapacklu", opts));

      MX Amx = new MX(DM.rand(3, 3));
      if (hasLinsol("lapacklu"))
        assertOutput(included, excluded, () -> ca.evalf(ca.solve(Amx, new MX(b), "lapacklu", opts)));

      if (hasLinsol("lapacklu")) {
        Linsol solver = new Linsol("linsol", "lapacklu", A.sparsity(), opts);
        assertOutput(included, excluded, () -> solver.solve(A, b));
      }
    }
  }

  @Test
  public void test_error_formatting() {
    MX x = MX.sym("x");
    Function f = new Function("f", List.of(x), List.of(x.attachAssert(0.0, "Hey \\W %d")));
    assertInException("Hey \\W %d", () -> call0(f));
  }

  @Test
  public void test_record_time() {
    MX x = MX.sym("x");
    MX f = ca.power(x, 2.0);
    for (int i = 0; i < 10000; i++) f = ca.sin(f).times(f);
    f = f.plus(ca.power(x, 2.0));

    Map<String, Object> opts = Map.of("record_time", true);
    Function ff = new Function("f", List.of(x), List.of(f), opts);
    ff.invoke(3.0);
    assertTrue(ff.stats().containsKey("t_proc_total"));
    assertTrue((Double) ff.stats().get("t_wall_total") >= 0);

    if (hasNlpsol("ipopt")) {
      Function solver = ca.nlpsol("nlpsol", "ipopt", Map.of("x", x, "f", f), opts);
      call0(solver);
      assertTrue(solver.stats().containsKey("t_proc_total"));
      assertTrue((Double) solver.stats().get("t_wall_total") > 0);
    }

    if (hasConic("qpoases")) {
      Function solver = ca.qpsol("qpsol", "qpoases", Map.of("x", x, "f", f), opts);
      call0(solver);
      assertTrue(solver.stats().containsKey("t_proc_total"));
      assertTrue((Double) solver.stats().get("t_wall_total") > 0);
    }

    Function solver = ca.rootfinder("rootfinder", "newton", Map.of("x", x, "g", x), opts);
    call0(solver);
    assertTrue(solver.stats().containsKey("t_proc_total"));
    assertTrue((Double) solver.stats().get("t_wall_total") >= 0);

    solver = ca.integrator("integrator", "rk", Map.of("x", x, "ode", f), opts);
    call0(solver);
    assertTrue(solver.stats().containsKey("t_proc_total"));
    assertTrue((Double) solver.stats().get("t_wall_total") >= 0);

    Map<String, Object> integr_options = new LinkedHashMap<>();
    integr_options.put("simplify", true);
    integr_options.put("simplify_options", opts);
    solver = ca.integrator("integrator", "rk", Map.of("x", x, "ode", f), integr_options);
    call0(solver);
    assertTrue(solver.stats().containsKey("t_proc_total"));
    assertTrue((Double) solver.stats().get("t_wall_total") >= 0);

    DM A = DM.rand(3, 3);
    if (hasLinsol("lapacklu")) {
      Linsol ls = new Linsol("linsol", "lapacklu", A.sparsity(), opts);
      ls.solve(A, dm(1, 2, 3));
      assertTrue(ls.stats().containsKey("t_proc_total"));
      assertTrue((Double) ls.stats().get("t_wall_total") >= 0);
    }
  }

  // test_stub_no_duplicate_parameter_names: python-only (casadi.pyi stubs)
  // test_stub_no_positional_only_syntax: python-only (casadi.pyi stubs)
  // test_stub_no_list_literal_annotations: python-only (casadi.pyi stubs)

  @Test
  public void test_unicode() throws IOException {
    String path = "üni/tüst.casadi";
    String path2 = "üni/stüst.casadi";
    if (!CasadiMeta.feature_list().contains("ghc-filesystem")) {
      path = "tüst.casadi";
      path2 = "stüst.casadi";
    }
    MX x = MX.sym("x");
    Function f = new Function("f", List.of(x), List.of(ca.power(x, 2.0)));
    f.save(path);
    Files.copy(new File(path).toPath(), new File(path2).toPath(), StandardCopyOption.REPLACE_EXISTING);
    Function f2 = Function.load(path2);
    checkfunctionLight(f, f2, List.of(new DM(3)));
  }
}

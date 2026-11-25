using System;
using System.Collections;
using System.Collections.Generic;
using System.IO;
using System.Linq;

namespace CasADi.Tests {

/// <summary>Port of test/java MiscTest (itself a port of test/python/misc.py)</summary>
public class MiscTest : CasadiTestCase {

  /// <summary>helpers.checkfunction_light: outputs agree</summary>
  static void checkfunctionLight(Function f, Function g, List<DM> inputs) {
    List<DM> a = f.call(inputs), b = g.call(inputs);
    assertEquals(b.Count, a.Count);
    for (int i = 0; i < a.Count; i++) checkarray(a[i], b[i], 1e-9, "output " + i);
  }

  static Dictionary<string, DM> call0(Function f) => f.call(new Dictionary<string, DM>());

  static void assertOutput(IEnumerable<string> included, IEnumerable<string> excluded, Action r) {
    string o = captureStdout(r);
    foreach (string s in included) assertTrue(o.Contains(s), "expected '" + s + "' in output:\n" + o);
    foreach (string s in excluded) assertFalse(o.Contains(s), "unexpected '" + s + "' in output:\n" + o);
  }

  static string captureStderr(Action r) => captureBoth(r)[1];

  public void test_issue179B() {
    message("Regression test #179 (B)");
    SX x = SX.sym("x");
    var f = ca.Function("f", new[] { x }, new[] { ca.pow(x, 2.0) });
    Sparsity sparsity = f.jac_sparsity(0, 0);
    sparsity.ToString();
  }

  // test_issue4216: python-only (complex numbers / numpy complex arrays)

  public void test_sanity() {
    DM v = new double[] { 0.738, 0.39, 0.99 };
    ca.DM(ca.Sparsity(4, 3, new long[] { 0, 2, 2, 3 }, new long[] { 1, 2, 1 }), v);
    Check.Throws<Exception>(() => ca.DM(ca.Sparsity(4, 4, new long[] { 0, 2, 2, 3 }, new long[] { 1, 2, 1 }), v));
    Check.Throws<Exception>(() => ca.DM(ca.Sparsity(4, 3, new long[] { 0, 2, 2, 12 }, new long[] { 1, 2, 1 }), v));
    Check.Throws<Exception>(() => ca.DM(ca.Sparsity(4, 3, new long[] { -10, 2, 2, 3 }, new long[] { 1, 2, 1 }), v));
    Check.Throws<Exception>(() => ca.DM(ca.Sparsity(4, 3, new long[] { 0, 2, 2, 3 }, new long[] { 8, 2, 1 }), v));
    Check.Throws<Exception>(() => ca.DM(ca.Sparsity(4, 3, new long[] { 0, 2, 2, 3 }, new long[] { -3, 2, 1 }), v));
    Check.Throws<Exception>(() => ca.DM(ca.Sparsity(4, 3, new long[] { 0, 2, 2, 3 }, new long[] { 1, 2, 1, 2 }), v));
    Check.Throws<Exception>(() => ca.DM(ca.Sparsity(4, 3, new long[] { 0, 2, 0, 3 }, new long[] { 1, 2, 1 }), v));
  }

  public void test_copyconstr_norefcount() {
    message("Copy constructor for non-refcounted classes");
    DM x = DM.ones(2, 3);
    var y = ca.DM(x);
    x[0, 0] = 5;
    assertTrue(!ReferenceEquals(x, y));
    assertEquals(5.0, (double)x[0, 0]);
    assertEquals(1.0, (double)y[0, 0]);
  }

  public void test_copyconstr_refcount() {
    message("Copy constructor for refcounted classes");
    Sparsity x = Sparsity.diag(4);
    var y = ca.Sparsity(x);
    x.resize(2, 8);
    assertTrue(!ReferenceEquals(x, y));
    assertTrue(x.numel() != 0 && y.numel() != 0);
    assertEquals((2L, 8L), x.size(), "shape");
    assertEquals((4L, 4L), y.size(), "shape");
  }

  // test_copy_norefcount / test_copy_refcount / test_deepcopy_*: python-only (copy protocol)
  // test_options_introspection: not ported (as in Java)
  // test_pickling: python-only (pickle protocol)

  public void test_exceptions() {
    // Java: ca.vertcat(List<GenericMatrix>{DM, SX, MX}) fails at runtime with "vertcat(".
    // In C# there is no common matrix base type, so a mixed SX/MX vertcat does not compile.
    Check.Skip("not applicable: mixed SX/MX vertcat is a compile-time error in C#");
  }

  public void test_getscheme() {
    SX x = SX.sym("x"), p = SX.sym("p");
    var F = ca.Function("F", new[] { x, p }, new[] { x + p, ca.pow(x, 2.0) }, new[] { "x", "p" }, new[] { "f", "g" });
    var fc = F.call(new Dictionary<string, DM> { ["x"] = 3, ["p"] = 4 });
    checkarray(fc["f"], ca.vertcat(7));
    checkarray(fc["g"], ca.vertcat(9));
  }

  public void test_assertions() {
    MX x = MX.sym("x");
    MX z = ca.pow(x, 2.0);
    z = z.attachAssert(z > 3.0, "x must be larger than 3");
    MX v = ca.sin(z);
    var f = ca.Function("f", new[] { x }, new[] { v });
    Console.WriteLine(f);
    f.Call(-6.0);
    try {
      f.Call(1.0);
    } catch (SwigException e) {
      Console.WriteLine(e.Message);
      assertTrue(e.Message.Contains("x must be larger than 3"));
    }
  }

  // test_doc: python-only

  public void test_output() {
    Check.Assume(hasNlpsol("ipopt"), "no nlpsol ipopt");
    string result = captureStdout(() => ca.vertcat(1, 2).print_dense());
    assertTrue(result.Contains("2"));
    SX x = SX.sym("x");
    var f = new Dictionary<string, SX> { ["x"] = x, ["f"] = ca.pow(x, 2.0) };
    Function solver = ca.nlpsol("solver", "ipopt", f);
    result = captureStdout(() => solver.call(new Dictionary<string, DM> { ["x0"] = 0 }));
    assertTrue(result.Contains("Number of nonzeros in equality constraint"), result);
    assertTrue(result.Contains("iter    objective    inf_pr"), result);
    string err = captureStderr(() => {
      try {
        ca.nlpsol("solver", "foo", f);
      } catch (SwigException) {
        // ignored as in Python
      }
    });
    assertTrue(err.Contains("casadi_nlpsol_foo"), err);
  }

  public void test_serialize() {
    Sparsity x = Sparsity.upper(3);
    SX y = SX.sym("y");
    var si = ca.StringSerializer();
    si.pack(new[] { x });
    SX z = y;
    for (int i = 0; i < 10000; i++) z = ca.sin(z);
    SX e = ca.vertcat(y, z, z * 2.0);
    var fref = ca.Function("f", new[] { y }, new[] { e });
    si.pack(e);
    string data = si.encode();
    si.pack(new[] { x });

    var sd = ca.StringDeserializer(data);
    var spx = (IList)sd.unpack();
    assertTrue(((Sparsity)spx[0]).is_equal(x));
    SX A = (SX)sd.unpack();
    var f = ca.Function("f", new[] { A[0] }, new[] { A });
    checkfunctionLight(f, fref, new List<DM> { ca.DM(7) });
    assertInException("end of stream", () => sd.unpack());
    assertInException("end of stream", () => sd.unpack());

    var fs = ca.FileSerializer("foo.dat");
    fs.pack(new[] { x });
    z = ca.sin(y);
    e = ca.vertcat(y, z, z * 2.0);
    fref = ca.Function("f", new[] { y }, new[] { e });
    fs.pack(e);
    fs.Dispose();
    var fd = ca.FileDeserializer("foo.dat");
    spx = (IList)fd.unpack();
    assertTrue(((Sparsity)spx[0]).is_equal(x));
    A = (SX)fd.unpack();
    f = ca.Function("f", new[] { A[0] }, new[] { A });
    checkfunctionLight(f, fref, new List<DM> { ca.DM(7) });
    assertInException("end of stream", () => fd.unpack());
    assertInException("end of stream", () => fd.unpack());
    fd.Dispose();

    SX xs = SX.sym("x");
    var s = ca.StringSerializer();
    s.pack(xs);
    string data1 = s.encode();
    s.pack(ca.sin(xs));
    string data2 = s.encode();
    var s2 = ca.StringDeserializer(data1);
    s2.unpack();
    assertInException("end of stream", () => s2.unpack());
    s2.decode(data2);
    s2.unpack();
    assertInException("end of stream", () => s2.unpack());

    fs = ca.FileSerializer("foo.dat");
    var f1 = ca.Function("f", new[] { xs }, new[] { ca.pow(xs, 2.0) });
    var g1 = ca.Function("g", new[] { xs }, new[] { ca.pow(xs, 3.0) });
    fs.pack(new[] { f1, g1 });
    fs.Dispose();
    assertInException("File is not loadable with 'load'. Use 'FileDeserializer' instead.", () => Function.load("foo.dat"));
    var fd2 = ca.FileDeserializer("foo.dat");
    Console.WriteLine(fd2.unpack());
    fd2.Dispose();

    fs = ca.FileSerializer("foo.dat");
    fs.pack(f1);
    fs.Dispose();
    Function r = Function.load("foo.dat");
    Console.WriteLine(r);

    f1.save("foo.dat");
    var fd3 = ca.FileDeserializer("foo.dat");
    Console.WriteLine(fd3.unpack());
    fd3.Dispose();
  }

  public void test_print_time() {
    MX x = MX.sym("x");
    MX f = ca.pow(x, 2.0);
    foreach (bool print_time in new[] { true, false }) {
      var included = print_time ? new[] { "t_wall" } : new string[0];
      var excluded = print_time ? new string[0] : new[] { "t_wall" };
      var opts = new Dictionary<string, object> { ["print_time"] = print_time };
      var ff = ca.Function("f", new[] { x }, new[] { f }, opts);
      assertOutput(included, excluded, () => ff.Call(3.0));

      if (hasNlpsol("ipopt")) {
        Function solver = ca.nlpsol("nlpsol", "ipopt", new Dictionary<string, MX> { ["x"] = x, ["f"] = f }, opts);
        assertOutput(included, excluded, () => call0(solver));
      }

      if (hasConic("ipopt")) {
        Function solver = ca.qpsol("qpsol", "qpoases", new Dictionary<string, MX> { ["x"] = x, ["f"] = f }, opts);
        assertOutput(included, excluded, () => call0(solver));
      }

      Function rf = ca.rootfinder("rootfinder", "newton", new Dictionary<string, MX> { ["x"] = x, ["g"] = x }, opts);
      assertOutput(included, excluded, () => call0(rf));

      Function in1 = ca.integrator("integrator", "rk", new Dictionary<string, MX> { ["x"] = x, ["ode"] = f }, opts);
      assertOutput(included, excluded, () => call0(in1));

      var integr_options = new Dictionary<string, object>(opts) { ["simplify"] = true };
      Function in2 = ca.integrator("integrator", "rk", new Dictionary<string, MX> { ["x"] = x, ["ode"] = f }, integr_options);
      assertOutput(included, excluded, () => call0(in2));

      DM A = DM.rand(3, 3);
      DM b = ca.vertcat(1, 2, 3);
      if (hasLinsol("lapacklu")) assertOutput(included, excluded, () => ca.solve(A, b, "lapacklu", opts));

      var Amx = ca.MX(DM.rand(3, 3));
      if (hasLinsol("lapacklu"))
        assertOutput(included, excluded, () => ca.evalf(ca.solve(Amx, ca.MX(b), "lapacklu", opts)));

      if (hasLinsol("lapacklu")) {
        var solver = ca.Linsol("linsol", "lapacklu", A.sparsity(), opts);
        assertOutput(included, excluded, () => solver.solve(A, b));
      }
    }
  }

  public void test_error_formatting() {
    MX x = MX.sym("x");
    var f = ca.Function("f", new[] { x }, new[] { x.attachAssert(0.0, "Hey \\W %d") });
    assertInException("Hey \\W %d", () => call0(f));
  }

  public void test_record_time() {
    MX x = MX.sym("x");
    MX f = ca.pow(x, 2.0);
    for (int i = 0; i < 10000; i++) f = ca.sin(f) * f;
    f = f + ca.pow(x, 2.0);

    var opts = new Dictionary<string, object> { ["record_time"] = true };
    var ff = ca.Function("f", new[] { x }, new[] { f }, opts);
    ff.Call(3.0);
    assertTrue(ff.stats().ContainsKey("t_proc_total"));
    assertTrue((double)ff.stats()["t_wall_total"] >= 0);

    if (hasNlpsol("ipopt")) {
      Function s = ca.nlpsol("nlpsol", "ipopt", new Dictionary<string, MX> { ["x"] = x, ["f"] = f }, opts);
      call0(s);
      assertTrue(s.stats().ContainsKey("t_proc_total"));
      assertTrue((double)s.stats()["t_wall_total"] > 0);
    }

    if (hasConic("qpoases")) {
      Function s = ca.qpsol("qpsol", "qpoases", new Dictionary<string, MX> { ["x"] = x, ["f"] = f }, opts);
      call0(s);
      assertTrue(s.stats().ContainsKey("t_proc_total"));
      assertTrue((double)s.stats()["t_wall_total"] > 0);
    }

    Function solver = ca.rootfinder("rootfinder", "newton", new Dictionary<string, MX> { ["x"] = x, ["g"] = x }, opts);
    call0(solver);
    assertTrue(solver.stats().ContainsKey("t_proc_total"));
    assertTrue((double)solver.stats()["t_wall_total"] >= 0);

    solver = ca.integrator("integrator", "rk", new Dictionary<string, MX> { ["x"] = x, ["ode"] = f }, opts);
    call0(solver);
    assertTrue(solver.stats().ContainsKey("t_proc_total"));
    assertTrue((double)solver.stats()["t_wall_total"] >= 0);

    var integr_options = new Dictionary<string, object> { ["simplify"] = true, ["simplify_options"] = opts };
    solver = ca.integrator("integrator", "rk", new Dictionary<string, MX> { ["x"] = x, ["ode"] = f }, integr_options);
    call0(solver);
    assertTrue(solver.stats().ContainsKey("t_proc_total"));
    assertTrue((double)solver.stats()["t_wall_total"] >= 0);

    DM A = DM.rand(3, 3);
    if (hasLinsol("lapacklu")) {
      var ls = ca.Linsol("linsol", "lapacklu", A.sparsity(), opts);
      ls.solve(A, ca.vertcat(1, 2, 3));
      assertTrue(ls.stats().ContainsKey("t_proc_total"));
      assertTrue((double)ls.stats()["t_wall_total"] >= 0);
    }
  }

  // test_stub_*: python-only (casadi.pyi stubs)

  public void test_unicode() {
    string path = "üni/tüst.casadi";
    string path2 = "üni/stüst.casadi";
    if (!CasadiMeta.feature_list().Contains("ghc-filesystem")) {
      path = "tüst.casadi";
      path2 = "stüst.casadi";
    }
    MX x = MX.sym("x");
    var f = ca.Function("f", new[] { x }, new[] { ca.pow(x, 2.0) });
    f.save(path);
    File.Copy(path, path2, true);
    Function f2 = Function.load(path2);
    checkfunctionLight(f, f2, new List<DM> { 3 });
  }
}

}

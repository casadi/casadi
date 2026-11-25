using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;

namespace CasADi.Tests {

/// <summary>Port of test/java TypemapsTest (test/python/typemaps.py)</summary>
public class TypemapsTest : CasadiTestCase {
  /// <summary>Python bool(x)</summary>
  static bool truth(dynamic x) => x.nonzero();

  /// <summary>Python's dynamic operator dispatch: C# overload resolution at runtime, on the runtime types</summary>
  static object op(string name, params object[] args) {
    dynamic[] a = args.Select(x => (dynamic)x).ToArray();
    switch (name) {
      case "plus": return ca.plus(a[0], a[1]);
      case "minus": return ca.minus(a[0], a[1]);
      case "times": return ca.times(a[0], a[1]);
      case "rdivide": return ca.rdivide(a[0], a[1]);
      case "power": return ca.power(a[0], a[1]);
      case "fmin": return ca.fmin(a[0], a[1]);
      case "fmax": return ca.fmax(a[0], a[1]);
      case "constpow": return ca.constpow(a[0], a[1]);
      case "atan2": return ca.atan2(a[0], a[1]);
      case "copysign": return ca.copysign(a[0], a[1]);
      case "mtimes": return ca.mtimes(a[0], a[1]);
      case "floor": return ca.floor(a[0]);
      case "ceil": return ca.ceil(a[0]);
      case "ge": return ca.ge(a[0], a[1]);
      case "gt": return ca.gt(a[0], a[1]);
      case "le": return ca.le(a[0], a[1]);
      case "lt": return ca.lt(a[0], a[1]);
      case "eq": return ca.eq(a[0], a[1]);
      case "ne": return ca.ne(a[0], a[1]);
    }
    throw new AssertFailed("no operation " + name);
  }

  static object neg(object z) => -(dynamic)z;

  static Function function(string name, List<object> inputs, List<object> outputs) {
    if (inputs[0] is SX) return ca.Function(name, inputs.Cast<SX>(), outputs.Cast<SX>());
    return ca.Function(name, inputs.Cast<MX>(), outputs.Cast<MX>());
  }

  static DM dmWith(Sparsity sp, double[] values) {
    var v = new double[sp.nnz()];
    Array.Copy(values, v, Math.Min(v.Length, values.Length));
    return ca.DM(sp, v);
  }

  static bool symbolic(object s) => s is SX || s is MX;

  static DM toDM(object s) => s is DM d ? d : ca.DM((double)s);

  /// <summary>Python doit(z, s, fun) of test_autoconversionMX</summary>
  void doit(object z, object s, Func<object, object, object> fun, string what) {
    object r = fun(z, s);
    assertEquals(z.GetType(), r.GetType(), what + ": expected " + z.GetType().Name);
    double[] dummy = { 1.3, 2.7, 9.4, 1.0 }, dummy2 = { 0.3, 2.4, 1.4, 1.7 };
    if (!symbolic(s)) {
      Function f = function("f", new List<object> { z }, new List<object> { r });
      DM f_out = f.call(new[] { dmWith(f.sparsity_in(0), dummy) })[0];
      Function f_ = function("f", new List<object> { z }, new List<object> { z });
      DM f__out = f_.call(new[] { dmWith(f_.sparsity_in(0), dummy) })[0];
      checkarray((DM)fun(f__out, toDM(s)), f_out, what);
    } else {
      Function f = function("f", new List<object> { z, s }, new List<object> { r });
      var f_in = new[] { dmWith(f.sparsity_in(0), dummy), dmWith(f.sparsity_in(1), dummy2) };
      DM f_out = f.call(f_in)[0];
      Function f_ = function("f", new List<object> { z, s }, new List<object> { z, s });
      List<DM> f__out = f_.call(f_in);
      checkarray((DM)fun(f__out[0], f__out[1]), f_out, what);
    }
  }

  static readonly List<(string, Func<object, object, object>)> ARITH = new List<(string, Func<object, object, object>)> {
    ("-z", (z, s) => neg(z)),
    ("z+s", (z, s) => op("plus", z, s)),
    ("s+z", (z, s) => op("plus", s, z)),
    ("s*z", (z, s) => op("times", s, z)),
    ("z*s", (z, s) => op("times", z, s)),
    ("z-s", (z, s) => op("minus", z, s)),
    ("s-z", (z, s) => op("minus", s, z)),
    ("z/s", (z, s) => op("rdivide", z, s)),
    ("s/z", (z, s) => op("rdivide", s, z)),
    ("z**s", (z, s) => op("power", z, s)),
    ("s**z", (z, s) => op("power", s, z)),
    ("fmin(s,z)", (z, s) => op("fmin", s, z)),
    ("fmax(s,z)", (z, s) => op("fmax", s, z)),
    ("constpow(s,z)", (z, s) => op("constpow", s, z)),
    ("constpow(z,s)", (z, s) => op("constpow", z, s)),
    ("atan2(s,z)", (z, s) => op("atan2", s, z)),
    ("atan2(z,s)", (z, s) => op("atan2", z, s)),
    ("copysign(z,s)", (z, s) => op("copysign", z, s)),
    ("copysign(s,z)", (z, s) => op("copysign", s, z)),
    ("z@s", (z, s) => op("mtimes", z, s)),
    ("s@z", (z, s) => op("mtimes", s, z)),
    ("floor(z/s)", (z, s) => op("floor", op("rdivide", z, s))),
    ("floor(z)", (z, s) => op("floor", z)),
    ("ceil(z)", (z, s) => op("ceil", z)),
  };

  static Func<object, object, object> arith(string k) => ARITH.First(e => e.Item1 == k).Item2;

  static readonly List<(string, Func<object, object, object>)> COMPARE = new List<(string, Func<object, object, object>)> {
    ("s>=z", (z, s) => op("ge", s, z)),
    ("s>z", (z, s) => op("gt", s, z)),
    ("s<=z", (z, s) => op("le", s, z)),
    ("s<z", (z, s) => op("lt", s, z)),
    ("s==z", (z, s) => op("eq", s, z)),
    ("s!=z", (z, s) => op("ne", s, z)),
  };

  /// <summary>Python nums: np.array 2x2 and DM 2x2 are both DM here, as are DM(4) and np.array(4); 4 and 4.0 are double</summary>
  static List<object> nums() => new List<object> { ca.DM(new double[,] { { 1, 2 }, { 3, 4 } }), ca.DM(4), 4.0 };

  static void assertSeq<T>(IEnumerable<T> expected, IEnumerable actual, string what = "") {
    var a = actual.Cast<object>().ToList();
    var e = expected.Cast<object>().ToList();
    Check.True(e.SequenceEqual(a), what + " expected [" + string.Join(", ", e) + "], got [" + string.Join(", ", a) + "]");
  }

  /// <summary>Java's IllegalArgumentException on null: rejected cleanly, not dereferenced</summary>
  static void rejectsNull(Action r) {
    try {
      r();
    } catch (Exception e) when (!(e is AssertFailed)) {
      Check.True(!(e is NullReferenceException), "null dereferenced: " + e.GetType().Name + " " + e.StackTrace);
      return;
    }
    throw new AssertFailed("null accepted");
  }

  public void test_memleak() {
    // C# analog of the refcount check: converting the same array repeatedly leaves it untouched
    double[,] a = { { 0, 0 } };
    for (int i = 0; i < 2; i++) {
      using (DM d = a)
      using (SX s = a) {
        assertEquals(2L, d.size2());
        assertEquals(2L, s.size2());
        assertEquals(2L, d.numel());
        assertEquals(2L, s.numel());
      }
    }
    assertEquals(0.0, a[0, 0]);
    assertEquals(0.0, a[0, 1]);
    DM d2 = a;
    d2.Dispose();
    d2.Dispose();
  }

  public void test_memleak_use_after_close() {
    var d = ca.DM(3);
    d.Dispose();
    Check.Throws<ObjectDisposedException>(() => d.size1());
  }

  public void test_sanity() {
    DM a = ca.vertcat(1, 2, 3);
    assertTrue(truth(ca.abs(ca.norm_inf(a)) > 1));
    var s = ca.SX(1);
    assertTrue(truth(s > 0));
    assertTrue(truth(ca.abs(s) > 0));
    var m = ca.MX(1);
    assertTrue(truth(m > 0));
    assertTrue(truth(ca.abs(m) > 0));
  }

  public void test_0() {
    message("Typemap double[,] -> DM");
    double[][,] arrays = { new double[,] { { 1, 2, 3 }, { 4, 5, 6 } }, new double[,] { { 1, 2 }, { 3, 4 }, { 5, 6 } },
                           new double[,] { { 3.2, 4.6, 9.9 } }, new double[,] { { 7 } }, new double[,] { { 1 }, { 2 } } };
    foreach (double[,] m in arrays) {
      DM zt = ca.transpose(ca.transpose(m));
      checkarray(zt, m);
      double[,] full = zt.ToArray();
      assertEquals(m.GetLength(0), full.GetLength(0), "rows of ToArray()");
      assertSeq(m.Cast<double>(), full, "DM(double[,]).ToArray()");
      assertSeq(m.Cast<double>(), zt.full(), "DM(double[,]).full()");
    }
    // python-only: int32/int dtypes and .sparse() (scipy)
  }

  public void test_1() {
    message("DM -> DM");
    var m = ca.DM(ca.Sparsity(4, 3, new long[] { 0, 2, 2, 3 }, new long[] { 1, 2, 1 }), ca.vertcat(3, 2.3, 8));
    DM zt = ca.transpose(ca.transpose(m));
    checkarray(m, zt, "DM(DM)");
    checkarray(zt.ToArray(), m, "DM(DM).ToArray()");
  }

  // test_2: python-only (scipy csr_matrix)
  // test_sparse_array: python-only (scipy sparse arrays)
  // test_setget: returns immediately in Python (get/set by reference dropped)

  public void test_conversion() {
    message("DM conversions");
    var w = ca.DM(ca.Sparsity(4, 3, new long[] { 0, 2, 2, 3 }, new long[] { 1, 2, 1 }), ca.vertcat(3, 2.3, 8));
    assertSeq(new double[] { 3, 2.3, 8 }, w.nonzeros());
    double[,] f = w.ToArray();
    assertEquals(4, f.GetLength(0));
    assertSeq(new double[] { 0, 0, 0, 3, 0, 8, 2.3, 0, 0, 0, 0, 0 }, f);
    assertSeq(new double[] { 0, 3, 2.3, 0, 0, 0, 0, 0, 0, 8, 0, 0 }, w.elements());
    double[,] d = { { 1, 2, 3 }, { 4, 5, 6 } };
    checkarray(d, d);
    checkarray(new double[] { 1, 2 }, new double[,] { { 1 }, { 2 } });
    checkarray(2.5, 2.5);
    assertEquals(2.5, ca.DM(2.5).to_double());
    assertEquals(2.5, (double)ca.DM(2.5));
    assertEquals(7L, ca.DM(7.0).to_int());
    // python-only: np.array(w), numpy.matrix(w), w.sparse()
  }

  public void test_conversion_dict() {
    // every Dict value type the typemaps accept, round-tripped through casadi
    SX x = SX.sym("x");
    var fx = ca.Function("fx", new[] { x }, new[] { x });
    var inp = new Dictionary<string, object> {
      ["bool"] = true,
      ["int"] = 3,
      ["long"] = 3000000000L,
      ["double"] = 2.5,
      ["string"] = "säλ",
      ["doubles"] = new double[] { 1.5, 2.5 },
      ["longs"] = new long[] { 1, 2 },
      ["ints"] = new int[] { 1, 2 },
      ["bools"] = new bool[] { true, false },
      ["strings"] = new List<string> { "a", "b" },
      ["string_array"] = new string[] { "a", "b" },
      ["list_long"] = new List<long> { 1, 2 },
      ["list_double"] = new List<double> { 1.5, 2.5 },
      ["list_bool"] = new List<bool> { true, false },
      ["doubles2"] = new double[][] { new double[] { 1, 2 }, new double[] { 3 } },
      ["longs2"] = new List<List<long>> { new List<long> { 1, 2 }, new List<long> { 3 } },
      ["dict"] = new Dictionary<string, object> { ["k"] = 1, ["n"] = new Dictionary<string, object> { ["j"] = "v" } },
      ["function"] = fx,
      ["functions"] = new List<Function> { fx, fx },
      ["mixed"] = new List<object> { "a", 1 },
    };
    Dictionary<string, object> o = Options.sanitize(inp);
    assertEquals<object>(true, o["bool"]);
    assertEquals<object>(3L, o["int"]);
    assertEquals<object>(3000000000L, o["long"]);
    assertEquals<object>(2.5, o["double"]);
    assertEquals<object>("säλ", o["string"]);
    assertTrue(o["doubles"] is double[], "doubles: " + o["doubles"].GetType());
    assertSeq(new double[] { 1.5, 2.5 }, (IEnumerable)o["doubles"]);
    assertTrue(o["longs"] is long[], "longs: " + o["longs"].GetType());
    assertSeq(new long[] { 1, 2 }, (IEnumerable)o["longs"]);
    assertTrue(o["ints"] is long[], "ints: " + o["ints"].GetType());
    assertSeq(new long[] { 1, 2 }, (IEnumerable)o["ints"]);
    assertTrue(o["bools"] is bool[], "bools: " + o["bools"].GetType());
    assertSeq(new bool[] { true, false }, (IEnumerable)o["bools"]);
    assertSeq(new[] { "a", "b" }, (IEnumerable)o["strings"]);
    assertSeq(new[] { "a", "b" }, (IEnumerable)o["string_array"]);
    assertTrue(o["list_long"] is long[], "list_long: " + o["list_long"].GetType());
    assertSeq(new long[] { 1, 2 }, (IEnumerable)o["list_long"]);
    assertTrue(o["list_double"] is double[], "list_double: " + o["list_double"].GetType());
    assertSeq(new double[] { 1.5, 2.5 }, (IEnumerable)o["list_double"]);
    assertTrue(o["list_bool"] is bool[], "list_bool: " + o["list_bool"].GetType());
    assertSeq(new bool[] { true, false }, (IEnumerable)o["list_bool"]);
    var doubles2 = ((IEnumerable)o["doubles2"]).Cast<object>().ToList();
    assertEquals(2, doubles2.Count);
    assertSeq(new double[] { 1, 2 }, (IEnumerable)doubles2[0], "doubles2[0]");
    assertSeq(new double[] { 3 }, (IEnumerable)doubles2[1], "doubles2[1]");
    var longs2 = ((IEnumerable)o["longs2"]).Cast<object>().ToList();
    assertEquals(2, longs2.Count);
    assertSeq(new long[] { 1, 2 }, (IEnumerable)longs2[0], "longs2[0]");
    assertSeq(new long[] { 3 }, (IEnumerable)longs2[1], "longs2[1]");
    var dict = (IDictionary<string, object>)o["dict"];
    assertEquals<object>(1L, dict["k"]);
    var n = (IDictionary<string, object>)dict["n"];
    assertEquals(1, n.Count);
    assertEquals<object>("v", n["j"]);
    assertEquals("fx", ((Function)o["function"]).name());
    assertEquals(2, ((IList)o["functions"]).Count);
    assertSeq(new object[] { "a", 1L }, (IEnumerable)o["mixed"]);
    Check.Throws<ArgumentException>(() => Options.sanitize(new Dictionary<string, object> { ["o"] = new object() }));
    Check.Throws<Exception>(() => Options.sanitize(new Dictionary<string, object> { ["o"] = ca.DM(3) }));
    // options with typed values reach casadi
    SX y = SX.sym("y");
    var f = ca.Function("f", new[] { x, y }, new[] { x + y },
                         new Dictionary<string, object> { ["default_in"] = new double[] { 3, 4 },
                                                          ["is_diff_in"] = new bool[] { true, false } });
    assertEquals(4.0, f.default_in(1));
    assertSeq(new bool[] { true, false }, f.is_diff_in());
    f = ca.Function("f", new[] { x, y }, new[] { x + y },
                     new Dictionary<string, object> { ["default_in"] = new List<int> { 3, 4 } });
    assertEquals(4.0, f.default_in(1));
    f = ca.Function("f", new[] { x, y }, new[] { x + y },
                     new Dictionary<string, object> { ["default_in"] = new long[] { 3, 4 } });
    assertEquals(4.0, f.default_in(1));
    assertInException("OT_STRING cannot be cast to OT_BOOL",
                      () => ca.Function("f", new[] { x }, new[] { x },
                                         new Dictionary<string, object> { ["verbose"] = "true" }));
  }

  public void test_conversion_dict_empty_arrays() {
    object sanitized(object v) => Options.sanitize(new Dictionary<string, object> { ["o"] = v })["o"];
    object r = sanitized(new double[0]);
    assertTrue(r is double[], "double[0] -> " + r?.GetType());
    r = sanitized(new long[0]);
    assertTrue(r is long[], "long[0] -> " + r?.GetType());
    r = sanitized(new string[0]);
    assertTrue(r is IEnumerable && !((IEnumerable)r).Cast<object>().Any(), "string[0] -> " + r?.GetType());
    r = sanitized(new double[0][]);
    assertTrue(r is IEnumerable && !((IEnumerable)r).Cast<object>().Any(), "double[0][] -> " + r?.GetType());
    ca.Function("f", new SX[0], new[] { ca.SX(1) }, new Dictionary<string, object> { ["default_in"] = new double[0] });
  }

  public void test_autoconversion() {
    message("Auto conversion DM");
    double x = 2.3;
    DM s = ca.DM(new double[,] { { 1, 2 }, { 3, 4 } });
    double[,] n = s.ToArray();
    double[,] times_ = new double[2, 2], div = new double[2, 2], minus_ = new double[2, 2], plus_ = new double[2, 2];
    for (int i = 0; i < 2; i++)
      for (int j = 0; j < 2; j++) {
        times_[i, j] = x * n[i, j];
        div[i, j] = x / n[i, j];
        minus_[i, j] = x - n[i, j];
        plus_[i, j] = x + n[i, j];
      }
    checkarray(x * s, s * x);
    checkarray(x * s, times_);
    checkarray(x / s, 1 / (s / x));
    checkarray(x / s, div);
    checkarray(x - s, -(s - x));
    checkarray(x - s, minus_);
    checkarray(x + s, s + x);
    checkarray(x + s, plus_);
    // Python in-place w += s on a float
    DM w = 2.3;
    w += s;
    checkarray(w, plus_);
    w = 2.3;
    w -= s;
    checkarray(w, minus_);
    w = 2.3;
    w *= s;
    checkarray(w, times_);
    w = 2.3;
    w /= s;
    checkarray(w, div);
  }

  public void test_matmul() {
    DM A = ca.DM(new double[,] { { 1, 3 }, { 4, 5 } });
    DM B = ca.DM(new double[,] { { 7, 2 }, { 0, 9 } });
    DM AB = ca.mtimes(A, B);
    checkarray(AB, new double[,] { { 7, 29 }, { 28, 53 } });
    SX As = SX.sym("x", 2, 2), Bs = SX.sym("x", 2, 2);
    foreach (SX y in new SX[] { ca.mtimes(As, Bs), ca.mtimes(As, B), ca.mtimes(A, Bs), ca.SX(ca.mtimes(A, B)) }) {
      var yf = ca.Function("y", new[] { As, Bs }, new[] { y });
      checkarray(yf.Call(A, B), AB);
    }
    MX Am = MX.sym("x", 2, 2), Bm = MX.sym("x", 2, 2);
    foreach (MX y in new MX[] { ca.mtimes(Am, Bm), ca.mtimes(Am, B), ca.mtimes(A, Bm), ca.MX(ca.mtimes(A, B)) }) {
      var yf = ca.Function("y", new[] { Am, Bm }, new[] { y });
      checkarray(yf.Call(A, B), AB);
    }
  }

  public void test_autoconversionMX() {
    message("Auto conversion MX");
    // python-only: in-place dunders, divmod, pow(z, s, mod), round, math.trunc
    foreach (object s in nums())
      foreach (object z in new object[] { SX.sym("x"), SX.sym("x", 2, 2) })
        foreach (var (k, fun) in ARITH) doit(z, s, fun, k);
    foreach (object s in nums())
      foreach (var (k, fun) in ARITH) doit(MX.sym("x", 2, 2), s, fun, k);
    foreach (object s in new object[] { SX.sym("x"), SX.sym("x"), SX.sym("x", 2, 2) })
      foreach (object z in new object[] { SX.sym("x"), SX.sym("x"), SX.sym("x", 2, 2) })
        foreach (var (k, fun) in ARITH) doit(z, s, fun, k);
    foreach (object s in new object[] { MX.sym("x"), MX.sym("x", 2, 2) })
      foreach (object z in new object[] { MX.sym("x"), MX.sym("x", 2, 2) })
        foreach (var (k, fun) in ARITH) doit(z, s, fun, k);
    DM sd = ca.DM(new double[,] { { 1, 2 }, { 3, 4 } });
    foreach (object z in new object[] { SX.sym("x", 2, 2), MX.sym("x", 2, 2) }) {
      doit(z, sd, (a, b) => neg(a), "-z");
      _ = -sd;
      foreach (string k in new[] { "z+s", "z+s", "s*z", "s*z", "z-s", "z-s", "z/s", "z/s" }) doit(z, sd, arith(k), k);
    }
    // statically typed int/double widening in both operand positions
    SX x = SX.sym("x");
    MX y = MX.sym("y");
    assertEquals(typeof(SX), ca.plus(4, x).GetType());
    assertEquals(typeof(SX), ca.plus(x, 4).GetType());
    assertEquals(typeof(SX), (x + 4).GetType());
    assertEquals(typeof(MX), ca.times(4, y).GetType());
    assertEquals(typeof(MX), (y * 4).GetType());
    assertEquals(typeof(SX), ca.plus(ca.DM(4), x).GetType());
    assertEquals(typeof(MX), ca.plus(ca.DM(4), y).GetType());
    assertEquals(typeof(SX), (ca.DM(4) + x).GetType());
    assertEquals(typeof(MX), (ca.DM(4) + y).GetType());
  }

  // test_numpy_preserve_type, test_numpy_mode_behaviour, test_numpy_warning_shown_once,
  // test_numpy_no_global_filter, test_numpy_unary_ufuncs_warn: python-only (numpy)

  public void test_operator_array_priority_deferral() {
    // python-only: __array_priority__ / ArrayInterface; the interop checks are ported
    DM D = ca.DM(new double[,] { { 1.0, 2.0 } });
    MX M = MX.sym("m", 1, 2);
    assertEquals(typeof(MX), ca.plus(D, M).GetType());
    assertEquals(typeof(MX), ca.plus(M, D).GetType());
    assertEquals(typeof(MX), (D + M).GetType());
    assertEquals(typeof(MX), (M + D).GetType());
    assertEquals(typeof(DM), ca.plus(D, ca.DM(new double[,] { { 3.0, 4.0 } })).GetType());
    assertEquals(typeof(DM), ca.plus(D, 2).GetType());
    assertEquals(typeof(DM), (D + 2).GetType());
    assertEquals(typeof(DM), ca.times(2, D).GetType());
  }

  // test_array_alias, test_squeeze, test_flat_flatten_ravel, test_experimental_numpy_array,
  // test_experimental_ndarray: python-only (ArrayInterface)
  // test_issue4268: python-only (reflected operators of foreign classes)

  public void test_conversion_operators() {
    message("COnversion operations");
    foreach (object s in nums())
      foreach (object z in new object[] { SX.sym("x"), SX.sym("x"), SX.sym("x", 2, 2) })
        foreach (var (k, fun) in COMPARE) assertEquals(typeof(SX), fun(z, s).GetType(), k);
    foreach (object s in nums())
      foreach (var (k, fun) in COMPARE) assertEquals(typeof(MX), fun(MX.sym("x", 2, 2), s).GetType(), k);
    foreach (object s in new object[] { SX.sym("x"), SX.sym("x"), SX.sym("x", 2, 2) })
      foreach (object z in new object[] { SX.sym("x"), SX.sym("x"), SX.sym("x", 2, 2) })
        foreach (var (k, fun) in COMPARE) assertEquals(typeof(SX), fun(z, s).GetType(), k);
    foreach (object s in new object[] { MX.sym("x"), MX.sym("x", 2, 2) })
      foreach (object z in new object[] { MX.sym("x"), MX.sym("x", 2, 2) })
        foreach (var (k, fun) in COMPARE) assertEquals(typeof(MX), fun(z, s).GetType(), k);
  }

  public void test_set() {
    message("DM set on dense matrices");
    double[,] goal = { { 1 }, { 2 }, { 3 } };
    DM w = goal;
    checkarray(w, goal);
    foreach (DM value in new[] { ca.vertcat(1, 2, 3), ca.DM(new double[,] { { 1, 2, 3 } }) }) {
      w.set_nz(.., value);
      checkarray(w, goal);
    }
    w = goal;
    checkarray(w, goal);
    w[.., ..] = goal;
    checkarray(w, goal);
  }

  public void test_DMSXcast() {
    message("Casting DM to SX");
    var W = ca.SX(ca.DM(new double[,] { { 1, 2, 3 }, { 4, 5, 6 } }));
    assertEquals(2L, W.size1());
    assertEquals(3L, W.size2());
  }

  public void test_DMMXcast() {
    message("Casting DM to MX");
    var W = ca.MX(ca.DM(new double[,] { { 1, 2, 3 }, { 4, 5, 6 } }));
    assertEquals(2L, W.size1());
    assertEquals(3L, W.size2());
  }

  public void test_DMSX() {
    message("Casting DM to SX");
    DM w = ca.DM(new double[,] { { 1, 2, 3 }, { 4, 5, 6 } });
    SX x = SX.sym("x");
    var f = ca.Function("f", new[] { x }, new SX[] { w });
    SX W = f.call(f.sx_in())[0];
    assertEquals(2L, W.size1());
    assertEquals(3L, W.size2());
  }

  public void test_DMMX() {
    message("Casting DM to MX");
    DM w = ca.DM(new double[,] { { 1, 2, 3 }, { 4, 5, 6 } });
    MX x = MX.sym("x");
    var f = ca.Function("f", new[] { x }, new MX[] { w });
    MX W = f.call(f.mx_in())[0];
    assertEquals(2L, W.size1());
    assertEquals(3L, W.size2());
  }

  public void test_setgetslice() {
    message("set/get on DM using slices");
    DM w = ca.DM(new double[,] { { 0, 0 } });
    double[,] A = { { 1.0, 2 }, { 3, 4 } };
    double[,] B = { { 4.0, 5 }, { 6, 7 } };
    w[.., ..] = new double[,] { { A[0, 0], A[0, 1] } };
    checkarray(w, new double[,] { { A[0, 0], A[0, 1] } });
    double[,] wf = w.ToArray();
    B[0, 0] = wf[0, 0];
    B[0, 1] = wf[0, 1];
    assertEquals(A[0, 0], B[0, 0], "get");
    assertEquals(A[0, 1], B[0, 1], "get");
    w = ca.DM(new double[,] { { 0 }, { 0 } });
    w[.., ..] = new double[] { A[0, 0], A[1, 0] };
    checkarray(w, new double[,] { { A[0, 0] }, { A[1, 0] } });
    wf = w.ToArray();
    for (int i = 0; i < 2; i++) B[i, 0] = wf[i, 0];
    for (int i = 0; i < 2; i++) assertEquals(A[i, 0], B[i, 0], "get");
    // python-only: numpy fancy slice assignment A[2:7:3,:7:4] = w
  }

  public void test_vertcatprecedence() {
    message("Argument precedence DM");
    DM a = ca.vertcat(1, 2);
    assertEquals(typeof(DM), ca.vertcat(a, a).GetType());
    assertEquals(typeof(DM), ca.vertcat(a, new double[] { 1, 2, 3 }).GetType());
    var b = ca.MX(ca.vertcat(1, 2));
    // Python vertcat(MX, [1,2,3]): the array converts to MX
    Console.WriteLine(ca.vertcat(b, new double[] { 1, 2, 3 }));
    assertEquals(typeof(MX), ca.vertcat(b, new double[] { 1, 2, 3 }).GetType());
    assertEquals(typeof(MX), ca.vertcat(b, ca.vertcat(1, 2, 3)).GetType());
  }

  public void test_issue190() {
    message("regression test issue #190");
    SX x = SX.sym("x");
    _ = x * 1;
    _ = x * 1.2;
    _ = SX.sym("x") * 1.0;
    _ = MX.sym("x") * 1.0;
  }

  public void test_array_cat() {
    ca.horzcat(SX.sym("x", 4, 3), ca.SX(DM.ones(4, 3)));
    assertEquals(typeof(SX), ca.horzcat(SX.sym("x", 4, 3), DM.ones(4, 3)).GetType());
  }

  public void test_issue() {
    message("std::vector<double> typemap.");
    long[] a = { 0, 2, 2, 3 };
    double[] b = { 0.738, 0.39, 0.99 };
    ca.DM(ca.Sparsity(4, 3, new long[] { 0, 2, 2, 3 }, new long[] { 1, 2, 1 }), new double[] { 0.738, 0.39, 0.99 });
    ca.DM(ca.Sparsity(4, 3, a, new long[] { 1, 2, 1 }), b);
    ca.DM(ca.Sparsity(4, 3, a, new long[] { 1, 2, 1 }), new double[,] { { b[0], b[1], b[2] } });
  }

  public void test_issue314() {
    message("regression test for #314: SX sparsity constructor");
    ca.SX(Sparsity.diag(3), ca.vertcat(1, 2, 3));
  }

  public void test_setAll_365() {
    message("ticket #365: DMAtrix.setAll does not work for 1x1 Matrices as input");
    foreach (DM i in new DM[] { ca.DM(4), new double[] { 4 }, new double[,] { { 4 } } }) {
      DM m = DM.ones(5, 5);
      m[.., ..] = i;
      checkarray(m, DM.ones(5, 5) * 4);
    }
    DM m2 = DM.ones(5, 5);
    m2[.., ..] = 4;
    checkarray(m2, DM.ones(5, 5) * 4);
    m2 = DM.ones(5, 5);
    m2[.., ..] = 4.0;
    checkarray(m2, DM.ones(5, 5) * 4);
  }

  // test_issue570: skipped on Python 3

  public void test_casting_DM() {
    message("casting DM");
    SX x = SX.sym("x");
    var f = ca.Function("f", new[] { x }, new[] { x });
    // python-only: the __DM__ protocol; C# analog: what a call converts, and what it rejects
    assertEquals(4.0, (double)f.Call(ca.DM(4)));
    assertEquals(4.0, (double)f.Call(4));
    assertEquals(4.0, (double)f.Call(4L));
    assertEquals(4.0, (double)f.Call(4.0f));
    assertEquals(4.0, (double)f.Call(new double[] { 4 }));
    assertEquals(4.0, (double)f.Call(new double[,] { { 4 } }));
    assertEquals(typeof(SX), f.Call(SX.sym("y")).GetType());
    assertEquals(typeof(MX), f.Call(MX.sym("y")).GetType());
    // new Object(), "4" and true do not compile
    rejectsNull(() => f.Call((DM)null));
    assertEquals(4.0, (double)f.call(new Dictionary<string, DM> { ["i0"] = 4 })["o0"]);
  }

  // test_casting_SX, test_casting_MX: python-only (__SX__/__MX__ protocols)
  // test_shape_attr_no_refleak: python-only (refcounts of a duck-typed .shape)

  public void test_bytes_typemaps() {
    message("bytes is neither a string nor a sequence of numbers");
    // python-only: bytes; C# analog: strings round-trip exactly, null is rejected
    SX x = SX.sym("x");
    assertEquals("y", SX.sym("y").name());
    assertEquals("y", MX.sym("y").name());
    ca.Function("f", new[] { x }, new[] { x }, new Dictionary<string, object> { ["jit"] = false });
    string u = "äλ€😀";
    assertEquals(u, SX.sym(u).name());
    assertEquals(u, MX.sym(u).name());
    assertEquals(u, SX.sym(u).ToString());
    assertEquals<object>(u, Options.sanitize(new Dictionary<string, object> { ["o"] = u })["o"]);
    assertTrue(Options.sanitize(new Dictionary<string, object> { [u] = 1 }).ContainsKey(u));
    assertEquals("a\0b", SX.sym("a\0b").name());
    rejectsNull(() => SX.sym(null));
    rejectsNull(() => MX.sym(null));
    assertEquals("f", ca.Function("f", new[] { x }, new[] { x }).name());
  }

  public void test_bytes_typemaps_exception_message() {
    SX x = SX.sym("x");
    string u = "äλ€😀";
    assertInException("Got '" + u + "'", () => ca.Function(u, new[] { x }, new[] { x }));
  }

  public void test_OUTPUT() {
    message("OUTPUT typemap");
    SX a = SX.sym("A", 3, 3);
    var (q, r) = ca.qr(a);
    assertEquals(typeof(SX), q.GetType());
    assertEquals(typeof(SX), r.GetType());
  }

  // test_cvar: python-only (SWIG cvar)

  public void test_ufuncsum() {
    message("ufunc.add");
    checkarray(ca.vertcat(1, 2, 3).nonzeros().Sum(), ca.DM(6));
  }

  static DM val(SX a) {
    var f = ca.Function("f", new SX[0], new[] { a });
    return f.call(new DM[0])[0];
  }

  public void test_sxmatrix() {
    double[,] a = { { 1, 2 }, { 3, 4 } }, aT = { { 1, 3 }, { 2, 4 } };
    for (int i = 0; i < 3; i++) {
      Console.WriteLine(val(a));
      Console.WriteLine(val(aT));
      checkarray(val(a), ca.DM(a));
      checkarray(val(((SX)aT).T()), ca.DM(a));
    }
  }

  public void test_issue1158() {
    var a = ca.DM(0, 2);
    assertEquals((0L, 2L), a.size());
    assertEquals((2L, 0L), ((DM)new double[2, 0]).size());
    assertEquals((1L, 0L), ((DM)new double[1, 0]).size());
    assertEquals((0L, 2L), ((DM)new double[0, 2]).size());
    assertEquals((0L, 1L), ((DM)new double[0]).size());
  }

  public void test_matrices() {
    foreach (double[,] D in new[] { new double[,] { { 1, 2 }, { 3, 4 } }, new double[,] { { 1, 2 }, { 3, 4.0 } } }) {
      DM d = DM.ones(2, 2);
      SX x = SX.sym("x", d.sparsity());
      ca.Function("f", new[] { x }, new[] { x });
      var fin = ca.DM(d.sparsity(), 0);
      fin[.., ..] = D;
      checkarray(fin, ca.DM(new double[,] { { 1, 2 }, { 3, 4 } }));
      d[.., ..] = D;
      checkarray(d, ca.DM(new double[,] { { 1, 2 }, { 3, 4 } }));
    }
    // python-only: numpy.matrix and scipy csc_matrix inputs
  }

  // test_issue1217: python-only (numpy object array of SX)
  // test_issue1373: python-only (np.array(DM))

  public void test_None() {
    // python-only: np.atleast_2d(None); C# analog: null matrices are rejected, not dereferenced
    rejectsNull(() => ca.repmat((DM)null, 1, 1));
    rejectsNull(() => ca.repmat((SX)null, 1, 1));
    rejectsNull(() => ca.repmat((MX)null, 1, 1));
    rejectsNull(() => ca.plus((DM)null, ca.DM(1)));
    rejectsNull(() => ca.DM((Sparsity)null));
    rejectsNull(() => ca.DM((double[])null));
    Check.Throws<SwigException>(() => ca.repmat((DM)(double[])null, 1, 1));
    rejectsNull(() => ca.DM((double[,])null));
    rejectsNull(() => ca.DM(new double[][] { new double[] { 1, 2 }, null }));
    rejectsNull(() => Sparsity.triplet(2, 2, null, null));
    rejectsNull(() => ca.vertcat(ca.DM(1), null));
    rejectsNull(() => ca.DM(3).get((Slice)null));
    SX x = SX.sym("x");
    var f = ca.Function("f", new[] { x }, new[] { x });
    rejectsNull(() => f.call(new DM[] { null }));
    var m = new Dictionary<string, DM> { ["i0"] = null };
    rejectsNull(() => f.call(m));
    rejectsNull(() => ca.Function(null, new[] { x }, new[] { x }));
  }

  public void testGenericTypeBoolean() {
    Check.Assume(hasNlpsol("ipopt"), "no nlpsol ipopt");
    SX x = SX.sym("x");
    Check.Throws<Exception>(() => ca.nlpsol("mysolver", "ipopt",
        new Dictionary<string, SX> { ["x"] = x, ["f"] = x * x },
        new Dictionary<string, object> { ["ipopt"] = new Dictionary<string, object> { ["acceptable_tol"] = SX.sym("x") } }));
    ca.nlpsol("mysolver", "ipopt", new Dictionary<string, SX> { ["x"] = x, ["f"] = x * x },
           new Dictionary<string, object> { ["ipopt"] = new Dictionary<string, object> { ["acceptable_tol"] = 1 } });
  }

  // test_to_longlong: commented out in Python

  public void test_issue_3214() {
    Check.Throws<Exception>(() => Console.WriteLine(ca.vertcat(-1.0, 0.0) - ca.DM(new double[,] { { 0.1, 0.0 } })));
    DM r = ca.vertcat(-1.0, 0.0) - ca.vertcat(0.1, 0.0);
    checkarray(r, ca.vertcat(-1.1, 0));
  }

  public void test_buglonglong() {
    SX x = SX.sym("x");
    ca.jacobian(x / 1.458151064450277e-12, x);
  }

  public void test_issue273() {
    // Python builds a numpy object array of SX; C# analog: blockcat of SX entries and numbers
    SX state = SX.sym("H_body_world", 6, 6);
    SX s0 = state[0, 0], s1 = state[1, 0], s2 = state[2, 0];
    var H = new SX[][] {
      new SX[] { ca.cos(s2) * ca.cos(s1), -ca.sin(s2) * ca.cos(s1), ca.sin(s1), state[3, 0] },
      new SX[] { ca.sin(s0) * ca.sin(s1) * ca.cos(s2) + ca.sin(s2) * ca.cos(s0), -ca.sin(s0) * ca.sin(s2) * ca.sin(s1) + ca.cos(s0) * ca.cos(s2),
                 -ca.sin(s0) * ca.cos(s1), state[4, 0] },
      new SX[] { ca.sin(s0) * ca.sin(s2) - ca.sin(s1) * ca.cos(s0) * ca.cos(s2), ca.sin(s0) * ca.cos(s2) + ca.sin(s2) * ca.sin(s1) * ca.cos(s0),
                 ca.cos(s0) * ca.cos(s1), state[5, 0] },
      new SX[] { 0, 0, 0, 1 },
    };
    SX M = ca.blockcat(H);
    Console.WriteLine(M);
    assertEquals((4L, 4L), M.size());
  }

  // test_issue_2625: python-only (numpy mode / ArrayInterface)
}

}

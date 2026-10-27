package org.casadi.test;

import static org.junit.jupiter.api.Assertions.assertArrayEquals;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertInstanceOf;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;
import static org.junit.jupiter.api.Assertions.fail;

import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;
import java.lang.reflect.Modifier;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.function.BinaryOperator;

import org.casadi.DM;
import org.casadi.Function;
import org.casadi.GenericMatrix;
import org.casadi.MX;
import org.casadi.Options;
import org.casadi.SX;
import org.casadi.Slice;
import org.casadi.Sparsity;
import org.casadi.ca;
import org.junit.jupiter.api.Assumptions;
import org.junit.jupiter.api.Test;

/** Port of test/python/typemaps.py */
public class TypemapsTest extends CasadiTestCase {

  /** Python bool(x) */
  private static boolean truth(Object x) {
    if (x instanceof org.casadi.DM) return ((org.casadi.DM) x).nonzero();
    if (x instanceof org.casadi.SX) return ((org.casadi.SX) x).nonzero();
    return ((org.casadi.MX) x).nonzero();
  }

  /** Python's dynamic operator dispatch: the static ca.name overload whose parameters match the runtime types exactly */
  private static Object op(String name, Object... args) {
    for (Method m : ca.class.getMethods()) {
      if (!m.getName().equals(name) || !Modifier.isStatic(m.getModifiers())) continue;
      Class<?>[] p = m.getParameterTypes();
      if (p.length != args.length) continue;
      boolean ok = true;
      for (int i = 0; i < p.length; i++)
        ok &= p[i] == args[i].getClass() || (p[i] == double.class && args[i] instanceof Double);
      if (!ok) continue;
      try {
        return m.invoke(null, args);
      } catch (InvocationTargetException e) {
        throw (RuntimeException) e.getCause();
      } catch (IllegalAccessException e) {
        throw new RuntimeException(e);
      }
    }
    StringBuilder t = new StringBuilder();
    for (Object a : args) t.append(a.getClass().getSimpleName()).append(",");
    fail("no overload ca." + name + "(" + t + ")");
    return null;
  }

  private static Object neg(Object z) {
    if (z instanceof DM) return ((DM) z).unaryMinus();
    if (z instanceof SX) return ((SX) z).unaryMinus();
    return ((MX) z).unaryMinus();
  }

  @SuppressWarnings("unchecked")
  private static Function function(String name, List<Object> in, List<Object> out) {
    if (in.get(0) instanceof SX) return new Function(name, (List<SX>) (List<?>) in, (List<SX>) (List<?>) out);
    return new Function(name, (List<MX>) (List<?>) in, (List<MX>) (List<?>) out);
  }

  private static DM dmWith(Sparsity sp, double[] values) {
    return new DM(sp, new DM(Arrays.copyOf(values, (int) sp.nnz())));
  }

  private static boolean symbolic(Object s) {
    return s instanceof SX || s instanceof MX;
  }

  private static DM toDM(Object s) {
    return s instanceof DM ? (DM) s : new DM((Double) s);
  }

  /** Python doit(z, s, fun) of test_autoconversionMX */
  private void doit(Object z, Object s, BinaryOperator<Object> fun, String what) {
    Object r = fun.apply(z, s);
    assertEquals(z.getClass(), r.getClass(), what + ": expected " + z.getClass().getSimpleName());
    double[] dummy = {1.3, 2.7, 9.4, 1.0}, dummy2 = {0.3, 2.4, 1.4, 1.7};
    if (!symbolic(s)) {
      Function f = function("f", List.of(z), List.of(r));
      DM f_out = f.call(List.of(dmWith(f.sparsity_in(0), dummy))).get(0);
      Function f_ = function("f", List.of(z), List.of(z));
      DM f__out = f_.call(List.of(dmWith(f_.sparsity_in(0), dummy))).get(0);
      checkarray((DM) fun.apply(f__out, toDM(s)), f_out, what);
    } else {
      Function f = function("f", List.of(z, s), List.of(r));
      List<DM> f_in = List.of(dmWith(f.sparsity_in(0), dummy), dmWith(f.sparsity_in(1), dummy2));
      DM f_out = f.call(f_in).get(0);
      Function f_ = function("f", List.of(z, s), List.of(z, s));
      List<DM> f__out = f_.call(f_in);
      checkarray((DM) fun.apply(f__out.get(0), f__out.get(1)), f_out, what);
    }
  }

  private static final Map<String, BinaryOperator<Object>> ARITH = new LinkedHashMap<>();
  static {
    ARITH.put("-z", (z, s) -> neg(z));
    ARITH.put("z+s", (z, s) -> op("plus", z, s));
    ARITH.put("s+z", (z, s) -> op("plus", s, z));
    ARITH.put("s*z", (z, s) -> op("times", s, z));
    ARITH.put("z*s", (z, s) -> op("times", z, s));
    ARITH.put("z-s", (z, s) -> op("minus", z, s));
    ARITH.put("s-z", (z, s) -> op("minus", s, z));
    ARITH.put("z/s", (z, s) -> op("rdivide", z, s));
    ARITH.put("s/z", (z, s) -> op("rdivide", s, z));
    ARITH.put("z**s", (z, s) -> op("power", z, s));
    ARITH.put("s**z", (z, s) -> op("power", s, z));
    ARITH.put("fmin(s,z)", (z, s) -> op("fmin", s, z));
    ARITH.put("fmax(s,z)", (z, s) -> op("fmax", s, z));
    ARITH.put("constpow(s,z)", (z, s) -> op("constpow", s, z));
    ARITH.put("constpow(z,s)", (z, s) -> op("constpow", z, s));
    ARITH.put("atan2(s,z)", (z, s) -> op("atan2", s, z));
    ARITH.put("atan2(z,s)", (z, s) -> op("atan2", z, s));
    ARITH.put("copysign(z,s)", (z, s) -> op("copysign", z, s));
    ARITH.put("copysign(s,z)", (z, s) -> op("copysign", s, z));
    ARITH.put("z@s", (z, s) -> op("mtimes", z, s));
    ARITH.put("s@z", (z, s) -> op("mtimes", s, z));
    ARITH.put("floor(z/s)", (z, s) -> op("floor", op("rdivide", z, s)));
    ARITH.put("floor(z)", (z, s) -> op("floor", z));
    ARITH.put("ceil(z)", (z, s) -> op("ceil", z));
  }

  private static final Map<String, BinaryOperator<Object>> COMPARE = new LinkedHashMap<>();
  static {
    COMPARE.put("s>=z", (z, s) -> op("ge", s, z));
    COMPARE.put("s>z", (z, s) -> op("gt", s, z));
    COMPARE.put("s<=z", (z, s) -> op("le", s, z));
    COMPARE.put("s<z", (z, s) -> op("lt", s, z));
    COMPARE.put("s==z", (z, s) -> op("eq", s, z));
    COMPARE.put("s!=z", (z, s) -> op("ne", s, z));
  }

  /** Python nums: np.array 2x2 and DM 2x2 are both DM here, as are DM(4) and np.array(4); 4 and 4.0 are double */
  private static List<Object> nums() {
    return List.of(dm(new double[][] {{1, 2}, {3, 4}}), new DM(4), 4.0);
  }

  @Test
  public void test_memleak() {
    // Java analog of the refcount check: converting the same array repeatedly leaves it untouched
    double[][] a = {{0, 0}};
    for (int i = 0; i < 2; i++) {
      try (DM d = new DM(a); SX s = new SX(a)) {
        assertEquals(2, d.size2());
        assertEquals(2, s.size2());
      }
    }
    assertArrayEquals(new double[] {0, 0}, a[0]);
    DM d = new DM(a);
    d.close();
    d.close();
  }

  @Test
  public void test_memleak_use_after_close() {
    DM d = new DM(3);
    d.close();
    assertThrows(RuntimeException.class, () -> d.size1());
  }

  @Test
  public void test_sanity() {
    DM a = dm(1, 2, 3);
    assertTrue(truth(ca.gt(ca.abs(ca.norm_inf(a)), 1)));
    SX s = new SX(1);
    assertTrue(truth(ca.gt(s, 0)));
    assertTrue(truth(ca.gt(ca.abs(s), 0)));
    MX m = new MX(1);
    assertTrue(truth(ca.gt(m, 0)));
    assertTrue(truth(ca.gt(ca.abs(m), 0)));
  }

  @Test
  public void test_0() {
    message("Typemap double[][] -> DM");
    double[][][] arrays = {{{1, 2, 3}, {4, 5, 6}}, {{1, 2}, {3, 4}, {5, 6}}, {{3.2, 4.6, 9.9}}, {{7}}, {{1}, {2}}};
    for (double[][] m : arrays) {
      DM zt = ca.transpose(ca.transpose(new DM(m)));
      assertInstanceOf(DM.class, zt, "DM expected");
      checkarray(zt, m);
      assertArrayEquals(m, zt.full(), "DM(double[][]).full()");
    }
    // python-only: int32/int dtypes and .sparse() (scipy); Java has no int[][] constructor
  }

  @Test
  public void test_1() {
    message("DM -> DM");
    DM m = new DM(new Sparsity(4, 3, new long[] {0, 2, 2, 3}, new long[] {1, 2, 1}), dm(3, 2.3, 8));
    DM zt = ca.transpose(ca.transpose(m));
    assertInstanceOf(DM.class, zt, "DM expected");
    checkarray(m, zt, "DM(DM)");
    checkarray(new DM(zt.full()), m, "DM(DM).full()");
  }

  // test_2: python-only (scipy csr_matrix)

  // test_sparse_array: python-only (scipy sparse arrays)

  // test_setget: returns immediately in Python (get/set by reference dropped)

  @Test
  public void test_conversion() {
    message("DM conversions");
    DM w = new DM(new Sparsity(4, 3, new long[] {0, 2, 2, 3}, new long[] {1, 2, 1}), dm(3, 2.3, 8));
    assertArrayEquals(new double[] {3, 2.3, 8}, w.nonzeros());
    double[][] f = w.full();
    assertEquals(4, f.length);
    assertArrayEquals(new double[] {0, 0, 0}, f[0]);
    assertArrayEquals(new double[] {3, 0, 8}, f[1]);
    assertArrayEquals(new double[] {2.3, 0, 0}, f[2]);
    assertArrayEquals(new double[] {0, 3, 2.3, 0, 0, 0, 0, 0, 0, 8, 0, 0}, w.elements());
    double[][] d = {{1, 2, 3}, {4, 5, 6}};
    checkarray(new DM(d), d);
    checkarray(new DM(new double[] {1, 2}), new double[][] {{1}, {2}});
    checkarray(new DM(2.5), 2.5);
    assertEquals(2.5, new DM(2.5).to_double());
    assertEquals(7, new DM(7.0).to_int());
    // python-only: np.array(w), numpy.matrix(w), w.sparse()
  }

  @Test
  public void test_conversion_dict() {
    // Java-specific: every Dict value type the typemaps accept, round-tripped through casadi
    SX x = SX.sym("x");
    Function fx = new Function("fx", List.of(x), List.of(x));
    Map<String, Object> in = new LinkedHashMap<>();
    in.put("bool", true);
    in.put("int", 3);
    in.put("long", 3000000000L);
    in.put("double", 2.5);
    in.put("string", "säλ");
    in.put("doubles", new double[] {1.5, 2.5});
    in.put("longs", new long[] {1, 2});
    in.put("ints", new int[] {1, 2});
    in.put("bools", new boolean[] {true, false});
    in.put("strings", List.of("a", "b"));
    in.put("string_array", new String[] {"a", "b"});
    in.put("list_long", List.of(1L, 2L));
    in.put("list_double", List.of(1.5, 2.5));
    in.put("list_bool", List.of(true, false));
    in.put("doubles2", new double[][] {{1, 2}, {3}});
    in.put("longs2", List.of(List.of(1L, 2L), List.of(3L)));
    in.put("dict", Map.of("k", 1, "n", Map.of("j", "v")));
    in.put("function", fx);
    in.put("functions", List.of(fx, fx));
    in.put("mixed", List.of("a", 1));
    Map<String, Object> out = Options.sanitize(in);
    assertEquals(Boolean.TRUE, out.get("bool"));
    assertEquals(3L, out.get("int"));
    assertEquals(3000000000L, out.get("long"));
    assertEquals(2.5, out.get("double"));
    assertEquals("säλ", out.get("string"));
    assertArrayEquals(new double[] {1.5, 2.5}, (double[]) out.get("doubles"));
    assertArrayEquals(new long[] {1, 2}, (long[]) out.get("longs"));
    assertArrayEquals(new long[] {1, 2}, (long[]) out.get("ints"));
    assertArrayEquals(new boolean[] {true, false}, (boolean[]) out.get("bools"));
    assertEquals(List.of("a", "b"), out.get("strings"));
    assertEquals(List.of("a", "b"), out.get("string_array"));
    assertArrayEquals(new long[] {1, 2}, (long[]) out.get("list_long"));
    assertArrayEquals(new double[] {1.5, 2.5}, (double[]) out.get("list_double"));
    assertArrayEquals(new boolean[] {true, false}, (boolean[]) out.get("list_bool"));
    assertArrayEquals(new double[][] {{1, 2}, {3}}, (double[][]) out.get("doubles2"));
    assertArrayEquals(new long[][] {{1, 2}, {3}}, (long[][]) out.get("longs2"));
    Map<?, ?> dict = (Map<?, ?>) out.get("dict");
    assertEquals(1L, dict.get("k"));
    assertEquals(Map.of("j", "v"), dict.get("n"));
    assertEquals("fx", ((Function) out.get("function")).name());
    assertEquals(2, ((List<?>) out.get("functions")).size());
    assertEquals(List.of("a", 1L), out.get("mixed"));
    assertThrows(IllegalArgumentException.class, () -> Options.sanitize(Map.of("o", new Object())));
    assertThrows(IllegalArgumentException.class, () -> Options.sanitize(Map.of("o", new DM(3))));
    // options with typed values reach casadi
    SX y = SX.sym("y");
    Function f = new Function("f", List.of(x, y), List.of(x.plus(y)),
        Map.of("default_in", new double[] {3, 4}, "is_diff_in", new boolean[] {true, false}));
    assertEquals(4.0, f.default_in(1));
    assertArrayEquals(new boolean[] {true, false}, f.is_diff_in());
    f = new Function("f", List.of(x, y), List.of(x.plus(y)), Map.of("default_in", List.of(3, 4)));
    assertEquals(4.0, f.default_in(1));
    f = new Function("f", List.of(x, y), List.of(x.plus(y)), Map.of("default_in", new long[] {3, 4}));
    assertEquals(4.0, f.default_in(1));
    assertInException("OT_STRING cannot be cast to OT_BOOL",
        () -> new Function("f", List.of(x), List.of(x), Map.of("verbose", "true")));
  }

  @Test
  public void test_conversion_dict_empty_arrays() {
    assertInstanceOf(double[].class, Options.sanitize(Map.of("o", new double[0])).get("o"));
    assertInstanceOf(long[].class, Options.sanitize(Map.of("o", new long[0])).get("o"));
    assertInstanceOf(List.class, Options.sanitize(Map.of("o", new String[0])).get("o"));
    assertInstanceOf(double[][].class, Options.sanitize(Map.of("o", new double[0][])).get("o"));
    new Function("f", List.<SX>of(), List.of(new SX(1)), Map.of("default_in", new double[0]));
  }

  @Test
  public void test_autoconversion() {
    message("Auto conversion DM");
    double x = 2.3;
    DM s = dm(new double[][] {{1, 2}, {3, 4}});
    double[][] n = s.full();
    double[][] times = new double[2][2], div = new double[2][2], minus = new double[2][2], plus = new double[2][2];
    for (int i = 0; i < 2; i++)
      for (int j = 0; j < 2; j++) {
        times[i][j] = x * n[i][j];
        div[i][j] = x / n[i][j];
        minus[i][j] = x - n[i][j];
        plus[i][j] = x + n[i][j];
      }
    checkarray(ca.times(x, s), s.times(x));
    checkarray(ca.times(x, s), times);
    checkarray(ca.rdivide(x, s), ca.rdivide(1, s.div(x)));
    checkarray(ca.rdivide(x, s), div);
    checkarray(ca.minus(x, s), s.minus(x).unaryMinus());
    checkarray(ca.minus(x, s), minus);
    checkarray(ca.plus(x, s), s.plus(x));
    checkarray(ca.plus(x, s), plus);
    // Python in-place w += s on a float: same operator overloads
    checkarray(ca.plus(2.3, s), plus);
    checkarray(ca.minus(2.3, s), minus);
    checkarray(ca.times(2.3, s), times);
    checkarray(ca.rdivide(2.3, s), div);
  }

  @Test
  public void test_matmul() {
    DM A = dm(new double[][] {{1, 3}, {4, 5}});
    DM B = dm(new double[][] {{7, 2}, {0, 9}});
    DM AB = ca.mtimes(A, B);
    checkarray(AB, new double[][] {{7, 29}, {28, 53}});
    checkarray(A.mtimes(B), AB);
    SX As = SX.sym("x", 2, 2), Bs = SX.sym("x", 2, 2);
    for (SX y : List.of(ca.mtimes(As, Bs), ca.mtimes(As, B), ca.mtimes(A, Bs), new SX(ca.mtimes(A, B)), As.mtimes(Bs),
        As.mtimes(B), A.mtimes(Bs))) {
      Function yf = new Function("y", List.of(As, Bs), List.of(y));
      checkarray((DM) yf.invoke(A, B), AB);
    }
    MX Am = MX.sym("x", 2, 2), Bm = MX.sym("x", 2, 2);
    for (MX y : List.of(ca.mtimes(Am, Bm), ca.mtimes(Am, B), ca.mtimes(A, Bm), new MX(ca.mtimes(A, B)), Am.mtimes(Bm),
        Am.mtimes(B), A.mtimes(Bm))) {
      Function yf = new Function("y", List.of(Am, Bm), List.of(y));
      checkarray((DM) yf.invoke(A, B), AB);
    }
  }

  @Test
  public void test_autoconversionMX() {
    message("Auto conversion MX");
    // python-only: in-place dunders, divmod, pow(z, s, mod), round, math.trunc
    for (Object s : nums())
      for (Object z : List.of(SX.sym("x"), SX.sym("x", 2, 2)))
        for (Map.Entry<String, BinaryOperator<Object>> e : ARITH.entrySet()) doit(z, s, e.getValue(), e.getKey());
    for (Object s : nums())
      for (Map.Entry<String, BinaryOperator<Object>> e : ARITH.entrySet())
        doit(MX.sym("x", 2, 2), s, e.getValue(), e.getKey());
    for (Object s : List.of(SX.sym("x"), SX.sym("x"), SX.sym("x", 2, 2)))
      for (Object z : List.of(SX.sym("x"), SX.sym("x"), SX.sym("x", 2, 2)))
        for (Map.Entry<String, BinaryOperator<Object>> e : ARITH.entrySet()) doit(z, s, e.getValue(), e.getKey());
    for (Object s : List.of(MX.sym("x"), MX.sym("x", 2, 2)))
      for (Object z : List.of(MX.sym("x"), MX.sym("x", 2, 2)))
        for (Map.Entry<String, BinaryOperator<Object>> e : ARITH.entrySet()) doit(z, s, e.getValue(), e.getKey());
    DM s = dm(new double[][] {{1, 2}, {3, 4}});
    for (Object z : List.of(SX.sym("x", 2, 2), MX.sym("x", 2, 2))) {
      doit(z, s, (a, b) -> neg(a), "-z");
      s.unaryMinus();
      for (String k : new String[] {"z+s", "z+s", "s*z", "s*z", "z-s", "z-s", "z/s", "z/s"}) doit(z, s, ARITH.get(k), k);
    }
    // statically typed int/double widening in both operand positions
    SX x = SX.sym("x");
    MX y = MX.sym("y");
    assertInstanceOf(SX.class, ca.plus(4, x));
    assertInstanceOf(SX.class, ca.plus(x, 4));
    assertInstanceOf(SX.class, x.plus(4));
    assertInstanceOf(MX.class, ca.times(4, y));
    assertInstanceOf(MX.class, y.times(4));
    assertInstanceOf(SX.class, ca.plus(new DM(4), x));
    assertInstanceOf(MX.class, ca.plus(new DM(4), y));
    assertInstanceOf(SX.class, new DM(4).plus(x));
    assertInstanceOf(MX.class, new DM(4).plus(y));
  }

  // test_numpy_preserve_type: python-only (numpy mode / ArrayInterface)
  // test_numpy_mode_behaviour: python-only (numpy mode)
  // test_numpy_warning_shown_once: python-only (Python warnings)
  // test_numpy_no_global_filter: python-only (Python warnings)
  // test_numpy_unary_ufuncs_warn: python-only (numpy ufuncs)

  @Test
  public void test_operator_array_priority_deferral() {
    // python-only: __array_priority__ / ArrayInterface; the interop checks are ported
    DM D = dm(new double[][] {{1.0, 2.0}});
    MX M = MX.sym("m", 1, 2);
    assertInstanceOf(MX.class, ca.plus(D, M));
    assertInstanceOf(MX.class, ca.plus(M, D));
    assertInstanceOf(MX.class, D.plus(M));
    assertInstanceOf(MX.class, M.plus(D));
    assertInstanceOf(DM.class, ca.plus(D, dm(new double[][] {{3.0, 4.0}})));
    assertInstanceOf(DM.class, ca.plus(D, 2));
    assertInstanceOf(DM.class, D.plus(2));
    assertInstanceOf(DM.class, ca.times(2, D));
  }

  // test_array_alias: python-only (casadi.array / star import)
  // test_squeeze: python-only (ArrayInterface)
  // test_flat_flatten_ravel: python-only (ArrayInterface)
  // test_experimental_numpy_array: python-only (ArrayInterface)
  // test_experimental_ndarray: python-only (ArrayInterface)
  // test_issue4268: python-only (reflected operators of foreign classes; SX == MX is a compile error in Java)

  @Test
  public void test_conversion_operators() {
    message("COnversion operations");
    for (Object s : nums())
      for (Object z : List.of(SX.sym("x"), SX.sym("x"), SX.sym("x", 2, 2)))
        for (Map.Entry<String, BinaryOperator<Object>> e : COMPARE.entrySet())
          assertEquals(SX.class, e.getValue().apply(z, s).getClass(), e.getKey());
    for (Object s : nums())
      for (Map.Entry<String, BinaryOperator<Object>> e : COMPARE.entrySet())
        assertEquals(MX.class, e.getValue().apply(MX.sym("x", 2, 2), s).getClass(), e.getKey());
    for (Object s : List.of(SX.sym("x"), SX.sym("x"), SX.sym("x", 2, 2)))
      for (Object z : List.of(SX.sym("x"), SX.sym("x"), SX.sym("x", 2, 2)))
        for (Map.Entry<String, BinaryOperator<Object>> e : COMPARE.entrySet())
          assertEquals(SX.class, e.getValue().apply(z, s).getClass(), e.getKey());
    for (Object s : List.of(MX.sym("x"), MX.sym("x", 2, 2)))
      for (Object z : List.of(MX.sym("x"), MX.sym("x", 2, 2)))
        for (Map.Entry<String, BinaryOperator<Object>> e : COMPARE.entrySet())
          assertEquals(MX.class, e.getValue().apply(z, s).getClass(), e.getKey());
  }

  @Test
  public void test_set() {
    message("DM set on dense matrices");
    double[][] goal = {{1}, {2}, {3}};
    DM w = new DM(goal);
    checkarray(w, goal);
    for (DM value : List.of(dm(1, 2, 3), dm(new double[][] {{1, 2, 3}}))) {
      w.set_nz(value, false, new Slice());
      checkarray(w, goal);
    }
    w = new DM(goal);
    checkarray(w, goal);
    w.set(new Slice(), new Slice(), new DM(goal));
    checkarray(w, goal);
  }

  @Test
  public void test_DMSXcast() {
    message("Casting DM to SX");
    SX W = new SX(dm(new double[][] {{1, 2, 3}, {4, 5, 6}}));
    assertEquals(2, W.size1());
    assertEquals(3, W.size2());
  }

  @Test
  public void test_DMMXcast() {
    message("Casting DM to MX");
    MX W = new MX(dm(new double[][] {{1, 2, 3}, {4, 5, 6}}));
    assertEquals(2, W.size1());
    assertEquals(3, W.size2());
  }

  @Test
  public void test_DMSX() {
    message("Casting DM to SX");
    DM w = dm(new double[][] {{1, 2, 3}, {4, 5, 6}});
    SX x = SX.sym("x");
    Function f = new Function("f", List.of(x), List.of(new SX(w)));
    SX W = f.call(f.sx_in()).get(0);
    assertEquals(2, W.size1());
    assertEquals(3, W.size2());
  }

  @Test
  public void test_DMMX() {
    message("Casting DM to MX");
    DM w = dm(new double[][] {{1, 2, 3}, {4, 5, 6}});
    MX x = MX.sym("x");
    Function f = new Function("f", List.of(x), List.of(new MX(w)));
    MX W = f.call(f.mx_in()).get(0);
    assertEquals(2, W.size1());
    assertEquals(3, W.size2());
  }

  @Test
  public void test_setgetslice() {
    message("set/get on DM using slices");
    DM w = dm(new double[][] {{0, 0}});
    double[][] A = {{1.0, 2}, {3, 4}};
    double[][] B = {{4.0, 5}, {6, 7}};
    w.set(new Slice(), new Slice(), new DM(new double[][] {A[0]}));
    checkarray(w, new double[][] {A[0]});
    B[0] = w.full()[0];
    assertArrayEquals(A[0], B[0], "get");
    w = dm(new double[][] {{0}, {0}});
    w.set(new Slice(), new Slice(), new DM(new double[] {A[0][0], A[1][0]}));
    checkarray(w, new double[][] {{A[0][0]}, {A[1][0]}});
    double[][] wf = w.full();
    for (int i = 0; i < 2; i++) B[i][0] = wf[i][0];
    for (int i = 0; i < 2; i++) assertEquals(A[i][0], B[i][0], "get");
    // python-only: numpy fancy slice assignment A[2:7:3,:7:4] = w
  }

  @Test
  public void test_vertcatprecedence() {
    message("Argument precedence DM");
    DM a = dm(1, 2);
    assertInstanceOf(DM.class, ca.vertcat(List.of(a, a)));
    assertInstanceOf(DM.class, ca.vertcat(List.of(a, new DM(new double[] {1, 2, 3}))));
    MX b = new MX(dm(1, 2));
    // Python vertcat(MX, [1,2,3]): Java lists are homogeneous, the numbers are wrapped explicitly
    System.out.println(ca.vertcat(List.of(b, new MX(new DM(new double[] {1, 2, 3})))));
    assertInstanceOf(MX.class, ca.vertcat(List.of(b, new MX(new DM(new double[] {1, 2, 3})))));
    List<GenericMatrix> mixed = List.of(b, new DM(new double[] {1, 2, 3}));
    assertInstanceOf(MX.class, ca.vertcat(mixed));
  }

  @Test
  public void test_issue190() {
    message("regression test issue #190");
    SX x = SX.sym("x");
    x.times(1);
    x.times(1.2);
    SX.sym("x").times(1.0);
    MX.sym("x").times(1.0);
  }

  @Test
  public void test_array_cat() {
    ca.horzcat(List.of(SX.sym("x", 4, 3), new SX(DM.ones(4, 3))));
    List<GenericMatrix> mixed = List.of(SX.sym("x", 4, 3), DM.ones(4, 3));
    assertInstanceOf(SX.class, ca.horzcat(mixed));
  }

  @Test
  public void test_issue() {
    message("std::vector<double> typemap.");
    long[] a = {0, 2, 2, 3};
    double[] b = {0.738, 0.39, 0.99};
    new DM(new Sparsity(4, 3, new long[] {0, 2, 2, 3}, new long[] {1, 2, 1}), new DM(new double[] {0.738, 0.39, 0.99}));
    new DM(new Sparsity(4, 3, a, new long[] {1, 2, 1}), new DM(b));
    new DM(new Sparsity(4, 3, a, new long[] {1, 2, 1}), new DM(new double[][] {b}));
  }

  @Test
  public void test_issue314() {
    message("regression test for #314: SX sparsity constructor");
    new SX(Sparsity.diag(3), dm(1, 2, 3));
  }

  @Test
  public void test_setAll_365() {
    message("ticket #365: DMAtrix.setAll does not work for 1x1 Matrices as input");
    for (DM i : List.of(new DM(4), new DM(new double[] {4}), new DM(new double[][] {{4}}))) {
      DM m = DM.ones(5, 5);
      m.set(new Slice(), new Slice(), i);
      checkarray(m, DM.ones(5, 5).times(4));
    }
    DM m = DM.ones(5, 5);
    m.set(new Slice(), new Slice(), 4);
    checkarray(m, DM.ones(5, 5).times(4));
    m = DM.ones(5, 5);
    m.set(new Slice(), new Slice(), 4.0);
    checkarray(m, DM.ones(5, 5).times(4));
  }

  // test_issue570: skipped on Python 3

  @Test
  public void test_casting_DM() {
    message("casting DM");
    SX x = SX.sym("x");
    Function f = new Function("f", List.of(x), List.of(x));
    // python-only: the __DM__ protocol; Java analog: what invoke converts, and what it rejects
    assertEquals(4.0, ((DM) f.invoke(new DM(4))).to_double());
    assertEquals(4.0, ((DM) f.invoke(4)).to_double());
    assertEquals(4.0, ((DM) f.invoke(4L)).to_double());
    assertEquals(4.0, ((DM) f.invoke(4.0f)).to_double());
    assertEquals(4.0, ((DM) f.invoke(new double[] {4})).to_double());
    assertEquals(4.0, ((DM) f.invoke((Object) new double[][] {{4}})).to_double());
    assertInstanceOf(SX.class, f.invoke(SX.sym("y")));
    assertInstanceOf(MX.class, f.invoke(MX.sym("y")));
    assertThrows(IllegalArgumentException.class, () -> f.invoke(new Object()));
    assertThrows(IllegalArgumentException.class, () -> f.invoke("4"));
    assertThrows(IllegalArgumentException.class, () -> f.invoke(true));
    assertThrows(IllegalArgumentException.class, () -> f.invoke((Object) null));
    assertEquals(4.0, ((DM) f.invoke(Map.of("i0", 4)).get("o0")).to_double());
    assertThrows(IllegalArgumentException.class, () -> f.invoke(Map.of("i0", "4")));
  }

  // test_casting_SX: python-only (__SX__ protocol; Function outputs are statically typed in Java)
  // test_casting_MX: python-only (__MX__ protocol; Function outputs are statically typed in Java)
  // test_shape_attr_no_refleak: python-only (refcounts of a duck-typed .shape)

  @Test
  public void test_bytes_typemaps() {
    message("bytes is neither a string nor a sequence of numbers");
    // python-only: bytes; Java analog: strings round-trip exactly, null is rejected
    SX x = SX.sym("x");
    assertEquals("y", SX.sym("y").name());
    assertEquals("y", MX.sym("y").name());
    new Function("f", List.of(x), List.of(x), Map.of("jit", false));
    String u = "äλ€😀";
    assertEquals(u, SX.sym(u).name());
    assertEquals(u, MX.sym(u).name());
    assertEquals(u, SX.sym(u).toString());
    assertEquals(u, Options.sanitize(Map.of("o", u)).get("o"));
    assertTrue(Options.sanitize(Map.of(u, 1)).containsKey(u));
    assertEquals("a\u0000b", SX.sym("a\u0000b").name());
    assertThrows(IllegalArgumentException.class, () -> SX.sym(null));
    assertThrows(IllegalArgumentException.class, () -> MX.sym(null));
    assertEquals("f", new Function("f", List.of(x), List.of(x)).name());
  }

  @Test
  public void test_bytes_typemaps_exception_message() {
    SX x = SX.sym("x");
    String u = "\u00e4\u03bb\u20ac\ud83d\ude00";
    assertInException("Got '" + u + "'", () -> new Function(u, List.of(x), List.of(x)));
  }

  @Test
  public void test_OUTPUT() {
    message("OUTPUT typemap");
    SX a = SX.sym("A", 3, 3);
    List<Object> r = ca.qr(a);
    assertEquals(2, r.size());
    assertInstanceOf(SX.class, r.get(0));
    assertInstanceOf(SX.class, r.get(1));
  }

  // test_cvar: python-only (SWIG cvar)

  @Test
  public void test_ufuncsum() {
    message("ufunc.add");
    checkarray(new DM(Arrays.stream(dm(1, 2, 3).nonzeros()).sum()), new DM(6));
  }

  private static DM val(SX a) {
    Function f = new Function("f", new ArrayList<SX>(), List.of(a));
    return f.call(new ArrayList<DM>()).get(0);
  }

  @Test
  public void test_sxmatrix() {
    double[][] a = {{1, 2}, {3, 4}}, aT = {{1, 3}, {2, 4}};
    for (int i = 0; i < 3; i++) {
      System.out.println(val(new SX(a)));
      System.out.println(val(new SX(aT)));
      checkarray(val(new SX(a)), dm(a));
      checkarray(val(new SX(aT).T()), dm(a));
    }
  }

  @Test
  public void test_issue1158() {
    // Java cannot express numpy.zeros((0,2)) as double[][] (no rows, so no column count)
    DM a = new DM(0, 2);
    assertEquals(0, a.size()[0]);
    assertEquals(2, a.size()[1]);
    assertArrayEquals(new long[] {2, 0}, new DM(new double[2][0]).size());
    assertArrayEquals(new long[] {1, 0}, new DM(new double[][] {{}}).size());
    assertArrayEquals(new long[] {0, 1}, new DM(new double[0]).size());
  }

  @Test
  public void test_matrices() {
    for (double[][] D : new double[][][] {{{1, 2}, {3, 4}}, {{1, 2}, {3, 4.0}}}) {
      DM d = DM.ones(2, 2);
      SX x = SX.sym("x", d.sparsity());
      new Function("f", List.of(x), List.of(x));
      DM fin = new DM(d.sparsity(), 0);
      fin.set(new Slice(), new Slice(), new DM(D));
      checkarray(fin, dm(new double[][] {{1, 2}, {3, 4}}));
      d.set(new Slice(), new Slice(), new DM(D));
      checkarray(d, dm(new double[][] {{1, 2}, {3, 4}}));
    }
    // python-only: numpy.matrix and scipy csc_matrix inputs
  }

  // test_issue1217: python-only (numpy object array of SX)
  // test_issue1373: python-only (np.array(DM))

  @Test
  public void test_None() {
    // python-only: np.atleast_2d(None); Java analog: null matrices are rejected, not dereferenced
    assertThrows(IllegalArgumentException.class, () -> ca.repmat((DM) null, 1, 1));
    assertThrows(IllegalArgumentException.class, () -> ca.repmat((SX) null, 1, 1));
    assertThrows(IllegalArgumentException.class, () -> ca.repmat((MX) null, 1, 1));
    assertThrows(IllegalArgumentException.class, () -> ca.plus((DM) null, new DM(1)));
    assertThrows(IllegalArgumentException.class, () -> new DM((Sparsity) null));
    assertThrows(IllegalArgumentException.class, () -> new DM((double[]) null));
    assertThrows(IllegalArgumentException.class, () -> new DM((double[][]) null));
    assertThrows(IllegalArgumentException.class, () -> new DM(new double[][] {{1, 2}, null}));
    assertThrows(IllegalArgumentException.class, () -> Sparsity.triplet(2, 2, null, null));
    assertThrows(IllegalArgumentException.class, () -> ca.vertcat(Arrays.asList(new DM(1), null)));
    assertThrows(IllegalArgumentException.class, () -> new DM(3).get((Slice) null));
    SX x = SX.sym("x");
    Function f = new Function("f", List.of(x), List.of(x));
    assertThrows(IllegalArgumentException.class, () -> f.call(Arrays.asList((DM) null)));
    Map<String, DM> m = new LinkedHashMap<>();
    m.put("i0", null);
    assertThrows(IllegalArgumentException.class, () -> f.call(m));
    assertThrows(IllegalArgumentException.class, () -> new Function(null, List.of(x), List.of(x)));
  }

  @Test
  public void testGenericTypeBoolean() {
    Assumptions.assumeTrue(hasNlpsol("ipopt"));
    SX x = SX.sym("x");
    assertThrows(RuntimeException.class, () -> ca.nlpsol("mysolver", "ipopt", Map.of("x", x, "f", x.times(x)),
        Map.of("ipopt", Map.of("acceptable_tol", SX.sym("x")))));
    ca.nlpsol("mysolver", "ipopt", Map.of("x", x, "f", x.times(x)), Map.of("ipopt", Map.of("acceptable_tol", 1)));
  }

  // test_to_longlong: commented out in Python

  @Test
  public void test_issue_3214() {
    assertThrows(RuntimeException.class, () -> System.out.println(ca.minus(dm(-1.0, 0.0), dm(new double[][] {{0.1, 0.0}}))));
    DM r = ca.minus(dm(-1.0, 0.0), dm(0.1, 0.0));
    checkarray(r, dm(-1.1, 0));
  }

  @Test
  public void test_buglonglong() {
    SX x = SX.sym("x");
    ca.jacobian(x.div(1.458151064450277e-12), x);
  }

  @Test
  public void test_issue273() {
    // Python builds a numpy object array of SX; Java analog: blockcat of SX entries and numbers
    SX state = SX.sym("H_body_world", 6, 6);
    SX s0 = state.get(0, 0), s1 = state.get(1, 0), s2 = state.get(2, 0);
    List<List<SX>> H = List.of(
        List.of(ca.cos(s2).times(ca.cos(s1)), ca.sin(s2).unaryMinus().times(ca.cos(s1)), ca.sin(s1), state.get(3, 0)),
        List.of(ca.sin(s0).times(ca.sin(s1)).times(ca.cos(s2)).plus(ca.sin(s2).times(ca.cos(s0))),
            ca.sin(s0).unaryMinus().times(ca.sin(s2)).times(ca.sin(s1)).plus(ca.cos(s0).times(ca.cos(s2))),
            ca.sin(s0).unaryMinus().times(ca.cos(s1)), state.get(4, 0)),
        List.of(ca.sin(s0).times(ca.sin(s2)).minus(ca.sin(s1).times(ca.cos(s0)).times(ca.cos(s2))),
            ca.sin(s0).times(ca.cos(s2)).plus(ca.sin(s2).times(ca.sin(s1)).times(ca.cos(s0))),
            ca.cos(s0).times(ca.cos(s1)), state.get(5, 0)),
        List.of(new SX(0), new SX(0), new SX(0), new SX(1)));
    SX M = ca.blockcat(H);
    System.out.println(M);
    assertArrayEquals(new long[] {4, 4}, M.size());
  }

  // test_issue_2625: python-only (numpy mode / ArrayInterface)
}

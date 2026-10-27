package org.casadi.test;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;
import static org.junit.jupiter.api.Assertions.fail;

import java.io.ByteArrayOutputStream;
import java.io.PrintStream;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import org.casadi.DM;
import org.casadi.Function;
import org.casadi.MX;
import org.casadi.SX;
import org.casadi.Sparsity;
import org.casadi.ca;

/** Shared helpers, after test/python/helpers.py casadiTestCase */
public abstract class CasadiTestCase {
  protected final Random random = new Random(1);

  protected static final boolean RUN_SLOW = System.getenv("CASADI_JAVA_RUN_SLOW") != null;

  protected static boolean hasNlpsol(String name) {
    return ca.has_nlpsol(name);
  }

  protected static boolean hasConic(String name) {
    return ca.has_conic(name);
  }

  protected static boolean hasIntegrator(String name) {
    return ca.has_integrator(name);
  }

  protected static boolean hasRootfinder(String name) {
    return ca.has_rootfinder(name);
  }

  protected static boolean hasLinsol(String name) {
    return ca.has_linsol(name);
  }

  protected static void message(String s) {
    System.out.println(s);
  }

  /** Column vector */
  protected static DM dm(double... x) {
    return new DM(x);
  }

  /** Dense matrix from rows */
  protected static DM dm(double[][] rows) {
    return new DM(rows);
  }

  protected static List<DM> dms(DM... x) {
    return new ArrayList<>(List.of(x));
  }

  /** Random dense matrix, normally distributed entries */
  protected DM randDM(long n, long m) {
    double[][] a = new double[(int) n][(int) m];
    for (int i = 0; i < n; i++)
      for (int j = 0; j < m; j++) a[i][j] = random.nextGaussian();
    return new DM(a);
  }

  /** Random matrix with about a fraction density of nonzeros */
  protected DM randDM(long n, long m, double density) {
    DM d = randDM(n, m);
    double[][] a = d.full();
    List<Long> r = new ArrayList<>(), c = new ArrayList<>();
    List<Double> v = new ArrayList<>();
    for (int j = 0; j < m; j++)
      for (int i = 0; i < n; i++)
        if (random.nextDouble() < density) {
          r.add((long) i);
          c.add((long) j);
          v.add(a[i][j]);
        }
    Sparsity sp = Sparsity.triplet(n, m, toLongs(r), toLongs(c));
    return new DM(sp, new DM(v.stream().mapToDouble(Double::doubleValue).toArray()));
  }

  protected static long[] toLongs(List<Long> x) {
    return x.stream().mapToLong(Long::longValue).toArray();
  }

  protected static void checkarray(DM a, DM b) {
    checkarray(a, b, 1e-10, "");
  }

  protected static void checkarray(DM a, DM b, String msg) {
    checkarray(a, b, 1e-10, msg);
  }

  protected static void checkarray(DM a, double[][] b) {
    checkarray(a, new DM(b), 1e-10, "");
  }

  protected static void checkarray(DM a, double b) {
    checkarray(a, new DM(b), 1e-10, "");
  }

  /** Same shape, entries equal to within tol (absolute, scaled by max(1,|b|)) */
  protected static void checkarray(DM a, DM b, double tol, String msg) {
    String where = msg.isEmpty() ? "" : msg + ": ";
    assertEquals(b.size1(), a.size1(), where + "size1 of " + a + " vs " + b);
    assertEquals(b.size2(), a.size2(), where + "size2 of " + a + " vs " + b);
    double[][] x = a.full(), y = b.full();
    for (int i = 0; i < x.length; i++)
      for (int j = 0; j < x[i].length; j++) {
        double u = x[i][j], v = y[i][j];
        if (Double.isNaN(u) && Double.isNaN(v)) continue;
        if (u == v) continue;
        if (!(Math.abs(u - v) <= tol * Math.max(1, Math.abs(v))))
          fail(where + "entry (" + i + "," + j + "): " + u + " vs " + v + "\n" + a + "\n" + b);
      }
  }

  protected static void checkarray(double a, double b, double tol) {
    assertTrue(Math.abs(a - b) <= tol * Math.max(1, Math.abs(b)), a + " vs " + b);
  }

  /** Asserts that r throws and its message contains s */
  protected static void assertInException(String s, Runnable r) {
    try {
      r.run();
    } catch (RuntimeException e) {
      assertTrue(String.valueOf(e.getMessage()).contains(s), "expected '" + s + "' in: " + e.getMessage());
      return;
    }
    fail("expected an exception containing '" + s + "'");
  }

  /** What r prints to System.out, where casadi printout goes */
  protected static String captureStdout(Runnable r) {
    PrintStream old = System.out;
    ByteArrayOutputStream b = new ByteArrayOutputStream();
    System.setOut(new PrintStream(b, true, StandardCharsets.UTF_8));
    try {
      r.run();
    } finally {
      System.out.flush();
      System.setOut(old);
    }
    return b.toString(StandardCharsets.UTF_8);
  }

  /** What r prints to System.out and System.err */
  protected static String[] captureBoth(Runnable r) {
    PrintStream out = System.out, err = System.err;
    ByteArrayOutputStream bo = new ByteArrayOutputStream(), be = new ByteArrayOutputStream();
    System.setOut(new PrintStream(bo, true, StandardCharsets.UTF_8));
    System.setErr(new PrintStream(be, true, StandardCharsets.UTF_8));
    try {
      r.run();
    } finally {
      System.out.flush();
      System.err.flush();
      System.setOut(out);
      System.setErr(err);
    }
    return new String[] {bo.toString(StandardCharsets.UTF_8), be.toString(StandardCharsets.UTF_8)};
  }

  /** helpers.py assertInAnyOutput: s in the exception message, stdout or stderr */
  protected static void assertInAnyOutput(String s, Runnable r) {
    String[] e = {""};
    String[] res = captureBoth(() -> {
      try {
        r.run();
      } catch (RuntimeException err) {
        e[0] = String.valueOf(err.getMessage());
      }
    });
    System.out.print(res[0]);
    assertTrue(e[0].contains(s) || res[0].contains(s) || res[1].contains(s),
        "error: <<<" + e[0] + ">>> stdout: <<<" + res[0] + ">>> stderr: <<<" + res[1] + ">>>");
  }

  protected static void assertOutput(String included, Runnable r) {
    String out = captureStdout(r);
    assertTrue(out.contains(included), "expected '" + included + "' in output:\n" + out);
  }

  /** Random inputs matching the sparsity of f's inputs */
  protected List<DM> randomInputs(Function f) {
    List<DM> r = new ArrayList<>();
    for (long i = 0; i < f.n_in(); i++) {
      Sparsity sp = f.sparsity_in(i);
      double[] v = new double[(int) sp.nnz()];
      for (int k = 0; k < v.length; k++) v[k] = random.nextGaussian();
      r.add(new DM(sp, new DM(v)));
    }
    return r;
  }

  /** f and g agree on inputs, also their first-order derivatives (as checkfunction in helpers.py) */
  protected void checkfunction(Function f, Function g, List<DM> inputs) {
    checkfunction(f, g, inputs, 1e-9);
  }

  protected void checkfunction(Function f, Function g, List<DM> inputs, double tol) {
    assertEquals(g.n_in(), f.n_in(), "n_in");
    assertEquals(g.n_out(), f.n_out(), "n_out");
    List<DM> a = f.call(inputs), b = g.call(inputs);
    for (int i = 0; i < a.size(); i++) checkarray(a.get(i), b.get(i), tol, "output " + i);
    List<DM> ja = jacobianInputs(f, inputs, a), jb = jacobianInputs(g, inputs, b);
    List<DM> va = f.jacobian().call(ja), vb = g.jacobian().call(jb);
    for (int i = 0; i < va.size(); i++) checkarray(va.get(i), vb.get(i), tol, "jacobian block " + i);
    // forward and reverse sweeps with the same seeds
    List<DM> fwd = new ArrayList<>(inputs), adj = new ArrayList<>(inputs);
    fwd.addAll(a);
    adj.addAll(a);
    for (long i = 0; i < f.n_in(); i++) fwd.add(randFor(f.sparsity_in(i)));
    for (long i = 0; i < f.n_out(); i++) adj.add(randFor(f.sparsity_out(i)));
    List<DM> fa = f.forward(1).call(fwd), fb = g.forward(1).call(fwd);
    for (int i = 0; i < fa.size(); i++) checkarray(fa.get(i), fb.get(i), tol, "forward " + i);
    List<DM> ra = f.reverse(1).call(adj), rb = g.reverse(1).call(adj);
    for (int i = 0; i < ra.size(); i++) checkarray(ra.get(i), rb.get(i), tol, "reverse " + i);
  }

  private DM randFor(Sparsity sp) {
    double[] v = new double[(int) sp.nnz()];
    for (int k = 0; k < v.length; k++) v[k] = random.nextGaussian();
    return new DM(sp, new DM(v));
  }

  private static List<DM> jacobianInputs(Function f, List<DM> inputs, List<DM> outputs) {
    List<DM> r = new ArrayList<>(inputs);
    r.addAll(outputs);
    return r;
  }

  protected static SX sxsym(String name, long n, long m) {
    return SX.sym(name, n, m);
  }

  protected static MX mxsym(String name, long n, long m) {
    return MX.sym(name, n, m);
  }
}

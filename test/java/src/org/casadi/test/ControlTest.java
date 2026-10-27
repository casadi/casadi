package org.casadi.test;

import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import org.casadi.DM;
import org.casadi.Function;
import org.casadi.MX;
import org.casadi.Slice;
import org.casadi.Sparsity;
import org.casadi.ca;
import org.junit.jupiter.api.Assumptions;
import org.junit.jupiter.api.Test;

/** Port of test/python/control.py */
public class ControlTest extends CasadiTestCase {

  private final List<Object[]> dplesolvers = new ArrayList<>();

  public ControlTest() {
    if (ca.has_dple("slicot")) dplesolvers.add(new Object[] {"slicot", Map.of("linear_solver", "csparse")});
  }

  private DM uniform(int n, int m) {
    double[][] a = new double[n][m];
    for (int i = 0; i < n; i++)
      for (int j = 0; j < m; j++) a[i][j] = random.nextDouble();
    return new DM(a);
  }

  private DM randstable(int n) {
    return randstable(n, 0.8, 0);
  }

  // scipy.linalg.orth(rand) -> Q factor of a QR decomposition
  private DM randstable(int n, double margin, double minimal) {
    double r = margin;
    DM A_ = ca.tril(uniform(n, n));
    for (int i = 0; i < n; i++) {
      double eig = 0;
      while (Math.abs(eig) <= minimal) eig = random.nextDouble() * r * 2 - r;
      A_.set(i, i, new DM(eig));
    }
    DM Q = (DM) ca.qr(uniform(n, n)).get(0);
    return ca.mtimes(List.of(Q, A_, Q.T()));
  }

  private static <T> List<T> isigma(List<T> a) {
    List<T> r = new ArrayList<>();
    r.add(a.get(a.size() - 1));
    r.addAll(a.subList(0, a.size() - 1));
    return r;
  }

  private static List<DM> inputList(Function f, Map<String, DM> inputs) {
    List<DM> r = new ArrayList<>();
    for (String n : f.name_in()) r.add(inputs.containsKey(n) ? inputs.get(n) : DM.zeros(f.sparsity_in(n)));
    return r;
  }

  private List<DM> symm(int n, int K) {
    List<DM> r = new ArrayList<>();
    for (int i = 0; i < K; i++) {
      DM v = uniform(n, n);
      r.add(ca.mtimes(v, v.T()));
    }
    return r;
  }

  private static MX stackShifted(MX AA, int n, int K) {
    long nn = (long) n * n, tot = nn * K;
    return ca.vertcat(List.of(AA.get(new Slice(tot - nn, tot), new Slice()), AA.get(new Slice(0, tot - nn), new Slice())));
  }

  @Test
  @SuppressWarnings("unchecked")
  public void test_dple_small() {
    Assumptions.assumeFalse(dplesolvers.isEmpty(), "no dple plugin (slicot) available");
    for (Object[] so : dplesolvers) {
      String Solver = (String) so[0];
      Map<String, Object> options = (Map<String, Object>) so[1];
      for (int K : RUN_SLOW ? new int[] {3, 4} : new int[] {2, 3}) {
        for (int n : new int[] {2, 3, 4}) {
          System.out.println(Solver + " " + options);
          random.setSeed(1);
          System.out.println(n + " " + K);
          List<DM> A_ = new ArrayList<>();
          for (int i = 0; i < K; i++) A_.add(randstable(n));
          List<DM> V_ = symm(n, K), V2_ = symm(n, K);
          Sparsity S = ca.kron(Sparsity.diag(K), Sparsity.dense(n, n));
          Function solver = ca.dplesol("solver", Solver, Map.of("a", S, "v", ca.repmat(S, 1, 2)), options);

          Map<String, DM> inputs = Map.of("a", ca.diagcat(A_), "v", ca.horzcat(List.of(ca.diagcat(V_), ca.diagcat(V2_))));

          MX As = MX.sym("A", S), Vs = MX.sym("V", S);
          List<MX> Vl = new ArrayList<>();
          for (MX i : isigma(ca.diagsplit(Vs, n))) Vl.add(i.plus(i.T()).div(2.0));
          MX Vss = ca.horzcat(Vl);
          List<MX> kr = new ArrayList<>();
          for (MX i : ca.diagsplit(As, n)) kr.add(ca.kron(i, i));
          MX AA = ca.diagcat(kr);
          MX A_total = ca.minus(DM.eye((long) n * n * K), stackShifted(AA, n, K));
          MX Pf = ca.solve(A_total, ca.vec(Vss), "csparse");
          MX P = ca.reshape(Pf, n, (long) K * n);
          P = ca.diagcat(ca.horzsplit(P, n));

          Map<String, MX> d = new LinkedHashMap<>();
          d.put("a", As);
          d.put("v", Vs);
          d.put("p", P);
          Function refsol = new Function("refsol", d, ca.dple_in(), ca.dple_out()).map("map", "serial", 2, List.of("a"), List.<String>of());

          checkfunction(solver, refsol, inputList(solver, inputs));
        }
      }
    }
  }

  @Test
  @SuppressWarnings("unchecked")
  public void test_dple_alt_small() {
    Assumptions.assumeFalse(dplesolvers.isEmpty(), "no dple plugin (slicot) available");
    for (Object[] so : dplesolvers) {
      String Solver = (String) so[0];
      Map<String, Object> options = (Map<String, Object>) so[1];
      for (int K : RUN_SLOW ? new int[] {3, 4} : new int[] {2, 3}) {
        for (int n : new int[] {2, 3, 4}) {
          System.out.println(Solver + " " + options);
          random.setSeed(1);
          System.out.println(n + " " + K);
          List<DM> A_ = new ArrayList<>();
          for (int i = 0; i < K; i++) A_.add(randstable(n));
          List<DM> V_ = symm(n, K);

          Map<String, DM> inputs = Map.of("a", ca.horzcat(A_), "v", ca.horzcat(V_));

          MX As = MX.sym("A", n, (long) n * K), Vs = MX.sym("V", n, (long) n * K);
          List<MX> Vl = new ArrayList<>();
          for (MX i : isigma(ca.horzsplit(Vs, n))) Vl.add(i.plus(i.T()).div(2.0));
          MX Vss = ca.horzcat(Vl);
          List<MX> kr = new ArrayList<>();
          for (MX i : ca.horzsplit(As, n)) kr.add(ca.kron(i, i));
          MX AA = ca.diagcat(kr);
          MX A_total = ca.minus(DM.eye((long) n * n * K), stackShifted(AA, n, K));
          MX Pf = ca.solve(A_total, ca.vec(Vss), "csparse");
          MX P = ca.reshape(Pf, n, (long) K * n);

          Map<String, MX> ds = new LinkedHashMap<>();
          ds.put("a", As);
          ds.put("v", Vs);
          ds.put("p", ca.horzcat(ca.dplesol(ca.horzsplit(As, n), ca.horzsplit(Vs, n), Solver, options)));
          Function solver = new Function("solver", ds, ca.dple_in(), ca.dple_out());
          Map<String, MX> dr = new LinkedHashMap<>();
          dr.put("a", As);
          dr.put("v", Vs);
          dr.put("p", P);
          Function refsol = new Function("refsol", dr, ca.dple_in(), ca.dple_out());

          checkfunction(solver, refsol, inputList(solver, inputs));
        }
      }
    }
  }
}

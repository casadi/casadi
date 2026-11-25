using System;
using System.Collections.Generic;
using System.Linq;

namespace CasADi.Tests {

/// <summary>Port of test/java ControlTest (test/python/control.py)</summary>
public class ControlTest : CasadiTestCase {
  readonly List<(string, Dictionary<string, object>)> dplesolvers = new List<(string, Dictionary<string, object>)>();

  // Reseedable generator (Java's random.setSeed(1))
  Random rng = new Random(1);

  public ControlTest() {
    if (ca.has_dple("slicot"))
      dplesolvers.Add(("slicot", new Dictionary<string, object> { ["linear_solver"] = "csparse" }));
  }

  DM uniform(int n, int m) {
    double[,] a = new double[n, m];
    for (int i = 0; i < n; i++)
      for (int j = 0; j < m; j++) a[i, j] = rng.NextDouble();
    return a;
  }

  DM randstable(int n) => randstable(n, 0.8, 0);

  // scipy.linalg.orth(rand) -> Q factor of a QR decomposition
  DM randstable(int n, double margin, double minimal) {
    double r = margin;
    DM A_ = ca.tril(uniform(n, n));
    for (int i = 0; i < n; i++) {
      double eig = 0;
      while (Math.Abs(eig) <= minimal) eig = rng.NextDouble() * r * 2 - r;
      A_[i, i] = eig;
    }
    var (Q, _) = ca.qr(uniform(n, n));
    return ca.mtimes(new[] { Q, A_, Q.T() });
  }

  static List<T> isigma<T>(List<T> a) {
    var r = new List<T> { a[a.Count - 1] };
    r.AddRange(a.Take(a.Count - 1));
    return r;
  }

  static List<DM> inputList(Function f, Dictionary<string, DM> inputs) =>
    f.name_in().Select(n => inputs.TryGetValue(n, out DM v) ? v : DM.zeros(f.sparsity_in(n))).ToList();

  List<DM> symm(int n, int K) {
    var r = new List<DM>();
    for (int i = 0; i < K; i++) {
      DM v = uniform(n, n);
      r.Add(ca.mtimes(v, v.T()));
    }
    return r;
  }

  static MX stackShifted(MX AA, int n, int K) {
    long nn = (long)n * n, tot = nn * K;
    return ca.vertcat(AA[ca.Slice(tot - nn, tot), ca.Slice()], AA[ca.Slice(0, tot - nn), ca.Slice()]);
  }

  public void test_dple_small() {
    Check.Assume(dplesolvers.Count > 0, "no dple plugin (slicot) available");
    foreach (var (Solver, options) in dplesolvers) {
      foreach (int K in RunSlow ? new[] { 3, 4 } : new[] { 2, 3 }) {
        foreach (int n in new[] { 2, 3, 4 }) {
          Console.WriteLine(Solver + " " + options);
          rng = new Random(1);
          Console.WriteLine(n + " " + K);
          var A_ = Enumerable.Range(0, K).Select(_ => randstable(n)).ToList();
          List<DM> V_ = symm(n, K), V2_ = symm(n, K);
          Sparsity S = ca.kron(Sparsity.diag(K), Sparsity.dense(n, n));
          Function solver = ca.dplesol("solver", Solver,
            new Dictionary<string, Sparsity> { ["a"] = S, ["v"] = ca.repmat(S, 1, 2) }, options);

          var inputs = new Dictionary<string, DM> {
            ["a"] = ca.diagcat(A_), ["v"] = ca.horzcat(ca.diagcat(V_), ca.diagcat(V2_)) };

          MX As = MX.sym("A", S), Vs = MX.sym("V", S);
          MX Vss = ca.horzcat(isigma(ca.diagsplit(Vs, n)).Select(i => (i + i.T()) / 2).ToList());
          MX AA = ca.diagcat(ca.diagsplit(As, n).Select(i => ca.kron(i, i)).ToList());
          MX A_total = DM.eye((long)n * n * K) - stackShifted(AA, n, K);
          MX Pf = ca.solve(A_total, ca.vec(Vss), "csparse");
          MX P = ca.reshape(Pf, n, (long)K * n);
          P = ca.diagcat(ca.horzsplit(P, n));

          var d = new Dictionary<string, MX> { ["a"] = As, ["v"] = Vs, ["p"] = P };
          Function refsol = ca.Function("refsol", d, ca.dple_in(), ca.dple_out())
            .map("map", "serial", 2, new[] { "a" }, new string[0]);

          checkfunction(solver, refsol, inputList(solver, inputs));
        }
      }
    }
  }

  public void test_dple_alt_small() {
    Check.Assume(dplesolvers.Count > 0, "no dple plugin (slicot) available");
    foreach (var (Solver, options) in dplesolvers) {
      foreach (int K in RunSlow ? new[] { 3, 4 } : new[] { 2, 3 }) {
        foreach (int n in new[] { 2, 3, 4 }) {
          Console.WriteLine(Solver + " " + options);
          rng = new Random(1);
          Console.WriteLine(n + " " + K);
          var A_ = Enumerable.Range(0, K).Select(_ => randstable(n)).ToList();
          List<DM> V_ = symm(n, K);

          var inputs = new Dictionary<string, DM> { ["a"] = ca.horzcat(A_), ["v"] = ca.horzcat(V_) };

          MX As = MX.sym("A", n, (long)n * K), Vs = MX.sym("V", n, (long)n * K);
          MX Vss = ca.horzcat(isigma(ca.horzsplit(Vs, n)).Select(i => (i + i.T()) / 2).ToList());
          MX AA = ca.diagcat(ca.horzsplit(As, n).Select(i => ca.kron(i, i)).ToList());
          MX A_total = DM.eye((long)n * n * K) - stackShifted(AA, n, K);
          MX Pf = ca.solve(A_total, ca.vec(Vss), "csparse");
          MX P = ca.reshape(Pf, n, (long)K * n);

          var ds = new Dictionary<string, MX> {
            ["a"] = As, ["v"] = Vs,
            ["p"] = ca.horzcat(ca.dplesol(ca.horzsplit(As, n), ca.horzsplit(Vs, n), Solver, options)) };
          var solver = ca.Function("solver", ds, ca.dple_in(), ca.dple_out());
          var dr = new Dictionary<string, MX> { ["a"] = As, ["v"] = Vs, ["p"] = P };
          var refsol = ca.Function("refsol", dr, ca.dple_in(), ca.dple_out());

          checkfunction(solver, refsol, inputList(solver, inputs));
        }
      }
    }
  }
}

}

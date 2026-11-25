using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;

namespace CasADi.Tests {

/// <summary>Port of test/java LinearsolverTest (itself a port of test/python/linearsolver.py)</summary>
public class LinearsolverTest : CasadiTestCase {

  sealed class LSolver {
    public readonly string name;
    public readonly Dictionary<string, object> options = new Dictionary<string, object>();
    public readonly HashSet<string> req;
    public LSolver(string name, IEnumerable<string> req) {
      this.name = name;
      this.req = new HashSet<string>(req);
    }
  }

  static readonly List<LSolver> lsolvers = new List<LSolver>();

  static LinearsolverTest() {
    // csparsecholesky is commented out in the Python original
    string[][] cands = { new[] { "csparse" }, new[] { "lapacklu" }, new[] { "lapackqr" }, new[] { "symbolicqr" },
      new[] { "lsqr" }, new[] { "ma27", "symmetry" }, new[] { "mumps", "symmetry" }, new[] { "qr" },
      new[] { "ldl", "posdef", "symmetry" } };
    foreach (string[] c in cands) {
      try {
        ca.load_linsol(c[0]);
        lsolvers.Add(new LSolver(c[0], c.Skip(1)));
      } catch (SwigException) {
        // plugin not available
      }
    }
    Console.WriteLine("linear solvers " + string.Join(" ", lsolvers.Select(l => l.name)));
  }

  static Function nullspacewrapper(string name, Sparsity sp, IDictionary<string, object> options) {
    SX a = SX.sym("a", sp);
    return ca.Function(name, new[] { a }, new[] { ca.nullspace(a) }, options);
  }

  static Function jacobianOld(Function f, long i, long j) {
    var o = new List<string> { "jac:" + f.name_out(j) + ":" + f.name_in(i) };
    o.AddRange(f.name_out());
    return f.factory(f.name() + "_jac", f.name_in(), o);
  }

  DM uniform(int n, int m) {
    double[,] a = new double[n, m];
    for (int i = 0; i < n; i++)
      for (int j = 0; j < m; j++) a[i, j] = random.NextDouble();
    return a;
  }

  static DM ev(Function f, params DM[] args) => f.Call(args);

  static List<DM> inputList(Function f, IDictionary<string, DM> inputs) => f.name_in().Select(n => inputs[n]).ToList();

  static void checkSerialize(Function F, List<DM> inputs) {
    Function F2 = Function.deserialize(F.serialize(new Dictionary<string, object> { ["debug"] = true }));
    List<DM> a = F.call(inputs), b = F2.call(inputs);
    for (int i = 0; i < F.n_out(); i++) checkarray(a[i], b[i], 1e-16, "serialize output " + i);
  }

  static void checkSerialize(Function F, IDictionary<string, DM> inputs) {
    Function F2 = Function.deserialize(F.serialize(new Dictionary<string, object> { ["debug"] = true }));
    Dictionary<string, DM> a = F.call(inputs), b = F2.call(inputs);
    assertTrue(F.name_out().SequenceEqual(F2.name_out()), "name_out");
    foreach (string k in F.name_out()) checkarray(a[k], b[k], 1e-16, "serialize output " + k);
  }

  // as helpers.py check_codegen: compile generated C into a shared library, compare via external
  static void checkCodegen(Function F, List<DM> inputs) {
    if (!RunSlow) return;
    string dir = Path.Combine(Path.GetTempPath(), "casadi_cg" + Guid.NewGuid().ToString("N"));
    Directory.CreateDirectory(dir);
    string name = "codegen_" + DateTime.Now.Ticks.ToString("x");
    var cg = ca.CodeGenerator(name + ".c");
    cg.add(F);
    cg.generate(dir + Path.DirectorySeparatorChar);
    string so = Path.Combine(dir, name + ".so");
    string inc = GlobalOptions.getCasadiIncludePath();
    var args = new List<string> { "-std=c89", "-fPIC", "-shared", "-O0" };
    if (inc.Length > 0) args.Add("-I" + inc);
    args.AddRange(new[] { Path.Combine(dir, name + ".c"), "-o", so, "-lm" });
    var psi = new ProcessStartInfo("gcc");
    foreach (string a in args) psi.ArgumentList.Add(a);
    using (var p = Process.Start(psi)) {
      p.WaitForExit();
      assertEquals(0, p.ExitCode, "gcc failed");
    }
    Function G = ca.external(F.name(), so);
    List<DM> ra = F.call(inputs), rb = G.call(inputs);
    for (int i = 0; i < F.n_out(); i++) checkarray(ra[i], rb[i], 1e-13, "codegen output " + i);
  }

  static void checkCodegen(Function F, IDictionary<string, DM> inputs) => checkCodegen(F, inputList(F, inputs));

  public void test_nullspace() {
    DM[] As = {
      ca.DM(new double[,] { { 1, 1.3 }, { 2, 5 }, { 1, 0.5 }, { 1.8, 1.7 } }),
      ca.DM(new double[,] { { 1, 1.3 }, { 2, 5 }, { 1, 0.5 } }),
      ca.DM(new double[,] { { 1, 1.3 }, { 2, 5 }, { 1, 0.5 }, { 0.2, 0.3 }, { -0.3, 0.7 } }),
      ca.DM(new double[,] { { 1, 0 }, { 0, 0 }, { 0, 1 }, { 0, 0 } }),
      ca.DM(new double[,] { { 1.3, 0, 0.4, 1 }, { 0.2, 0.1, 11, 0 }, { 0, 1, 0, 0 }, { 0.7, 0.9, 0, 0 }, { 1.1, 0.99, 0, 0 } }) };
    var options = new Dictionary<string, object>();
    foreach (DM A in As) {
      long n = A.size1(), m = A.size2();
      Function solver = nullspacewrapper("solver", A.T().sparsity(), options);
      List<DM> solver_out = solver.call(new[] { A.T() });

      checkarray(ca.mtimes(A.T(), solver_out[0]), DM.zeros(m, n - m));
      checkarray(ca.mtimes(solver_out[0].T(), solver_out[0]), DM.eye(n - m));

      options["ad_weight"] = 0;
      options["ad_weight_sp"] = 0;
      solver = nullspacewrapper("solver", A.T().sparsity(), options);
      Function Jf = jacobianOld(solver, 0, 0);

      options["ad_weight"] = 1;
      options["ad_weight_sp"] = 1;
      solver = nullspacewrapper("solver", A.T().sparsity(), options);
      Function Jb = jacobianOld(solver, 0, 0);

      List<DM> Jf_out = Jf.call(new[] { A.T() });
      List<DM> Jb_out = Jb.call(new[] { A.T() });
      checkarray(Jf_out[0], Jb_out[0]);
      checkarray(Jf_out[1], Jb_out[1]);

      Function d = solver.forward(1);
      DM r = uniform((int)n, (int)m);
      List<DM> d_out = d.call(new[] { A.T(), 5, r.T() });
      DM exact = d_out[0];

      DM nom = ev(solver, A.T());
      double eps = 1e-6;
      DM pert = ev(solver, (A + r * eps).T());
      DM fd = (pert - nom) / eps;

      DM V = uniform((int)(n - m), (int)(n - m));
      V = V + V.T();
      Console.WriteLine(V);
      DM a = ca.mtimes(new[] { nom, V, fd.T() }) + ca.mtimes(new[] { fd, V, nom.T() });
      DM b = ca.mtimes(new[] { nom, V, exact.T() }) + ca.mtimes(new[] { exact, V, nom.T() });
      Console.WriteLine("here: " + (a - b));

      V = DM.eye(n);
      a = ca.mtimes(new[] { nom.T(), V, fd }) + ca.mtimes(new[] { fd.T(), V, nom });
      b = ca.mtimes(new[] { nom.T(), V, exact }) + ca.mtimes(new[] { exact.T(), V, nom });
      checkarray(a, b, 1e-5, "");
    }
  }

  public void test_simple_solve() {
    DM A_ = ca.DM(new double[,] { { 3, 7 }, { 1, 2 } });
    DM b_ = ca.vertcat(1, 0.5);
    MX A = MX.sym("A", A_.sparsity()), b = MX.sym("b", b_.sparsity());
    foreach (LSolver s in lsolvers) {
      MX C = ca.solve(A, b, s.name, s.options);
      DM A0 = s.req.Contains("symmetry") ? A_.T() + A_ : A_;
      var f = ca.Function("f", new[] { A, b }, new[] { C });
      DM f_out = ev(f, A0, b_);
      checkarray(f_out, ca.solve(A0, b_));
      checkarray(ca.mtimes(A0, f_out), b_);
    }
  }

  public void test_pseudo_inverse() {
    // random.setSeed(0) in Java: the C# Random of CasadiTestCase cannot be reseeded
    DM A_ = uniform(4, 6);
    MX A = MX.sym("A", A_.sparsity());
    SX As = SX.sym("A", A_.sparsity());
    foreach (LSolver s in lsolvers) {
      if (s.req.Contains("symmetry")) continue;
      Console.WriteLine(s.name);
      MX B = ca.pinv(A, s.name, s.options);
      var f = ca.Function("f", new[] { A }, new[] { B });
      checkarray(ca.mtimes(A_, ev(f, A_)), DM.eye(4));
      f = ca.Function("f", new[] { As }, new[] { ca.pinv(As) });
      checkarray(ca.mtimes(A_, ev(f, A_)), DM.eye(4));
      ca.solve(ca.mtimes(A, A.T()), A, s.name, s.options);
      ca.pinv(A_, s.name, s.options);
    }

    A_ = uniform(3, 5);
    A = MX.sym("A", A_.sparsity());
    As = SX.sym("A", A_.sparsity());
    foreach (LSolver s in lsolvers) {
      Console.WriteLine(s.name);
      if (s.req.Contains("symmetry")) continue;
      MX B = ca.pinv(A, s.name, s.options);
      var f = ca.Function("f", new[] { A }, new[] { B });
      checkarray(ca.mtimes(A_, ev(f, A_)), DM.eye(3));
      f = ca.Function("f", new[] { As }, new[] { ca.pinv(As) });
      checkarray(ca.mtimes(A_, ev(f, A_)), DM.eye(3));
    }
  }

  public void test_simple_solve_dmatrix() {
    DM A = ca.DM(new double[,] { { 3, 7 }, { 1, 2 } });
    DM b = ca.vertcat(1, 0.5);
    foreach (LSolver s in lsolvers) {
      Console.WriteLine(s.name);
      DM A0 = s.req.Contains("symmetry") ? A.T() + A : A;
      DM C = ca.solve(A0, b, s.name, s.options);
      DM sol = ca.solve(A0, b);
      checkarray(C, sol);
      checkarray(ca.mtimes(A0, sol), b);
    }
  }

  public void test_simple_trans() {
    DM A = ca.DM(new double[,] { { 3, 1 }, { 7, 2 } });
    foreach (LSolver s in lsolvers) {
      DM A0 = s.req.Contains("symmetry") ? A.T() + A : A;
      var solver = ca.Linsol("solver", s.name, A0.sparsity(), s.options);
      assertInException("No stats available since Linsol did not solve a problem yet.", () => solver.stats());
      DM b = ca.vertcat(1, 0.5);
      DM x = solver.solve(A0.T(), b);
      DM res = ca.solve(A0.T(), b);
      Console.WriteLine(string.Join(", ", solver.stats().Select(kv => kv.Key + ": " + kv.Value)));
      checkarray(x, res);
    }
  }

  public void test_simple() {
    DM A = ca.DM(new double[,] { { 3, 1 }, { 7, 2 } });
    foreach (LSolver s in lsolvers) {
      Console.WriteLine(s.name);
      DM A0 = s.req.Contains("symmetry") ? A.T() + A : A;
      var solver = ca.Linsol("solver", s.name, A0.sparsity(), s.options);
      DM b = ca.vertcat(1, 0.5);
      DM x = solver.solve(A0, b);
      checkarray(x, ca.solve(A0, b));
    }
  }

  public void test_simple_function_indirect() {
    foreach (LSolver s in lsolvers) {
      DM A_ = ca.DM(new double[,] { { 3, 1 }, { 7, 2 } });
      MX A = MX.sym("A", A_.sparsity());
      DM b_ = ca.vertcat(1, 0.5);
      MX b = MX.sym("b", b_.sparsity());
      if (s.req.Contains("symmetry")) A_ = A_.T() + A_;
      Console.WriteLine(s.name);
      var solver = ca.Linsol("solver", s.name, A.sparsity(), s.options);
      var solver_in = new Dictionary<string, DM> { ["A"] = A_, ["B"] = b_ };

      var sol = new Dictionary<string, MX> { ["X"] = solver.solve(A, b), ["A"] = A, ["B"] = b };
      var relay = ca.Function("relay", sol, new[] { "A", "B" }, new[] { "X" });

      MX A_0 = A[0, 0], A_1 = A[0, 1], A_2 = A[1, 0], A_3 = A[1, 1];
      MX b_0 = b[0], b_1 = b[1];
      MX det = A_0 * A_3 - A_2 * A_1;
      MX X = ca.vertcat(A_3 / det * b_0 + -A_1 / det * b_1, -A_2 / det * b_0 + A_0 / det * b_1);
      var sd = new Dictionary<string, MX> { ["A"] = A, ["B"] = b, ["X"] = X };
      var solution = ca.Function("solution", sd, new[] { "A", "B" }, new[] { "X" });

      checkfunction(relay, solution, inputList(relay, solver_in));
      checkSerialize(relay, solver_in);

      if (s.name == "qr" || s.name == "ldl") checkCodegen(relay, solver_in);
    }
  }

  public void test_simple_solve_node() {
    foreach (LSolver s in lsolvers) {
      DM[][] cases = {
        new[] { ca.DM(new double[,] { { 3, 1 }, { 7, 2 } }), ca.DM(new double[,] { { 1, 0.3 }, { 0.5, 0.7 } }) },
        new[] { ca.sparsify(ca.DM(new double[,] { { 3, 0 }, { 7, 2 } })), ca.DM(new double[,] { { 1, 0.3 }, { 0.5, 0.7 } }) },
        new[] { ca.DM(new double[,] { { 3, 1 }, { 7, 2 } }), ca.sparsify(ca.DM(new double[,] { { 1, 0 }, { 0, 0.7 } })) } };
      foreach (DM[] c in cases) {
        DM A_ = c[0], b_ = c[1];
        if (s.req.Contains("symmetry")) A_ = A_.T() + A_;
        MX A = MX.sym("A", A_.sparsity()), b = MX.sym("b", b_.sparsity());
        Console.WriteLine(s.name);
        var solver = ca.Linsol("solver", s.name, A.sparsity(), s.options);
        foreach (bool tr in new[] { true, false }) {
          MX x = solver.solve(A, b, tr);
          var f = ca.Function("f", new[] { A, b }, new[] { x });
          DM f_out = ev(f, A_, b_);
          DM x_ = solver.solve(A_, b_, tr);
          checkarray(x_, f_out);

          MX A_0, A_1, A_2, A_3;
          if (tr) {
            A_0 = A[0, 0];
            A_1 = A[1, 0];
            A_2 = A[0, 1];
            A_3 = A[1, 1];
          } else {
            A_0 = A[0, 0];
            A_1 = A[0, 1];
            A_2 = A[1, 0];
            A_3 = A[1, 1];
          }
          MX b_0 = b[0, 0], b_1 = b[1, 0], c_0 = b[0, 1], c_1 = b[1, 1];
          MX det = A_0 * A_3 - A_2 * A_1;
          MX p = A_3 / det, q = -A_1 / det, u = -A_2 / det, v = A_0 / det;
          var solution = ca.Function("solution", new[] { A, b }, new[] { ca.blockcat(new[] {
            new[] { p * b_0 + q * b_1, p * c_0 + q * c_1 },
            new[] { u * b_0 + v * b_1, u * c_0 + v * c_1 } }) });

          checkfunction(f, solution, new List<DM> { A_, b_ });

          if (!s.name.Contains("SymbolicQR")) continue;
          Function solversx = f.expand("expand_" + f.name());
          checkfunction(solversx, solution, randomInputs(solversx), 1e-7);
        }
      }
    }
  }

  public void test_large_sparse() {
    int n = 10;
    // Java reseeds (setSeed(3)) per solver so all solvers see the same matrix, one for which ldl (no
    // pivoting) does not hit a zero pivot. The C# Random cannot be reseeded: draw until such a matrix comes up.
    DM A0, b0;
    while (true) {
      A0 = randDM(n, n, 0.5);
      b0 = randDM(n, 3, 0.5);
      if (!hasLinsol("ldl")) break;
      DM t = ca.solve(A0.T() + A0, b0, "ldl");
      if (t.nonzeros().All(v => !double.IsNaN(v) && !double.IsInfinity(v))) break;
    }
    foreach (LSolver s in lsolvers) {
      DM A = A0;
      if (s.req.Contains("symmetry")) A = A.T() + A;
      DM b = b0;
      MX As = MX.sym("A", A.sparsity()), bs = MX.sym("B", b.sparsity());
      DM C = ca.solve(A, b, s.name, s.options);
      checkarray(ca.mtimes(A, C), b, s.name);
      var f = ca.Function("f", new[] { As, bs }, new[] { ca.solve(As, bs, s.name, s.options) });
      checkarray(ca.mtimes(A, ev(f, A, b)), b, s.name);
    }
  }

  public void test_ma27() {
    double n = double.NaN;
    DM A = ca.DM(new double[,] { { 2, 3, n, n, n }, { 3, 0, 4, n, 6 }, { n, 4, 1, 5, n }, { n, n, 5, 0, n }, { n, 6, n, n, 1 } });
    A = A.get(ca.sparsify(ca.eq(A, A)).sparsity());
    DM b = ca.vertcat(8, 45, 31, 15, 17);
    DM @ref = ca.solve(ca.densify(A), b);
    foreach (LSolver s in lsolvers) {
      if (s.req.Contains("posdef")) continue;
      MX.sym("A", A.sparsity());
      MX.sym("B", b.sparsity());
      DM C = ca.solve(A, b, s.name, s.options);
      checkarray(@ref, C, s.name);
    }
  }

  public void test_large_sparse2() {
    int n = 10;
    foreach (LSolver s in lsolvers) {
      Console.WriteLine(s.name);
      DM A = randDM(n, n, 0.5);
      DM b = randDM(n, 3);
      if (s.req.Contains("posdef")) {
        A = ca.mtimes(A.T(), A);
        A = ca.densify(A);
      } else if (s.req.Contains("symmetry")) {
        A = A.T() + A;
        A.set(Sparsity.diag(n), 1e-8);
        A = ca.densify(A);
      }
      MX As = MX.sym("A", A.sparsity()), bs = MX.sym("B", b.sparsity());
      DM C = ca.solve(A, b, s.name, s.options);
      double tol = s.name.Contains("ma") ? 1e-7 : 1e-10;
      checkarray(ca.mtimes(A, C), b, tol, s.name);

      var sym = new List<MX> { As, ca.densify(As), ca.densify(As).T(), ca.densify(As.T()), As.T() };
      var num = new List<DM> { A, ca.densify(A), ca.densify(A).T(), ca.densify(A.T()), A.T() };
      for (int k = 0; k < sym.Count; k++) {
        var f = ca.Function("f", new[] { As, bs }, new[] { ca.solve(sym[k], bs, s.name, s.options) });
        checkarray(ca.mtimes(num[k], ev(f, A, b)), b, tol, s.name + " case " + k);
      }
    }
  }

  public void test_dimmismatch() {
    DM A = DM.eye(5);
    DM b = DM.ones(4, 1);
    MX As = MX.sym("A", A.sparsity()), bs = MX.sym("b", b.sparsity());
    foreach (LSolver s in lsolvers) {
      Check.Throws<Exception>(() => ca.solve(A, b, s.name, s.options));
      Check.Throws<Exception>(() => ca.solve(As, bs, s.name, s.options));
    }
  }

  public void test_cache() {
    int n = 5;
    var As = new List<DM>();
    for (int i = 0; i < 10; i++) As.Add(uniform(n, n));

    MX A = MX.sym("A", n, n);
    var f = ca.Function("f", new[] { A }, new[] { ca.solve(A, (MX)DM.ones(n), "qr", new Dictionary<string, object> { ["cache"] = 2 }) });

    DM S1 = ev(f, As[0]);
    checkarray(S1, ev(f, As[0]));

    DM S2 = ev(f, As[1]);
    checkarray(S2, ev(f, As[1]));
    checkarray(S1, ev(f, As[0]));

    DM S3 = ev(f, As[2]);
    checkarray(S3, ev(f, As[2]));
    checkarray(S1, ev(f, As[0]));
    checkarray(S2, ev(f, As[1]));

    checkCodegen(f, new List<DM> { As[0] });
    checkSerialize(f, new List<DM> { As[0] });
  }

  public void test_thread_safety() {
    MX x = MX.sym("x"), y = MX.sym("y");
    foreach (LSolver s in lsolvers) {
      if (s.name.Contains("mumps")) continue;
      Console.WriteLine(s.name);
      var f = ca.Function("f", new[] { x, y }, new[] { x * x - y });
      Function finv = ca.rootfinder("finv", "newton", f,
        new Dictionary<string, object> { ["linear_solver"] = s.name, ["linear_solver_options"] = s.options });

      Function finv_par = finv.map(200, "thread", 4);
      DM lin = ca.linspace(0.0, 10.0, 200);
      DM res = ev(finv_par, DM.ones(1, 200), lin.T());
      checkarray(ca.norm_inf(res.T() - ca.sqrt(lin)), 0, 1e-5, s.name);

      MX z = ca.solve(ca.blockcat(new[] { new[] { x + 1.0, (MX)0 }, new[] { (MX)0, y + 1.0 } }), ca.vertcat(x, y), s.name, s.options);
      f = ca.Function("f", new[] { x, y }, new[] { z });
      Function f_par = f.map(200, "thread", 4);
      f_par.Call(ca.linspace(10.0, 0.0, 200).T(), lin.T());
    }
  }

  public void test_issue2664() {
    DM bnum = DM.rand(3, 3);
    MX x = MX.sym("x", 3, 3);
    foreach (MX b in new[] { MX.sym("b", 3, 3), MX.sym("b", Sparsity.lower(3)), MX.sym("b", Sparsity.upper(3)) }) {
      var f = ca.Function("f", new[] { x, b }, new[] { ca.solve(x, b, "symbolicqr") });
      foreach (DM bn in new[] { bnum, bnum.get(Sparsity.lower(3)), bnum.get(Sparsity.upper(3)) }) {
        SX d = ca.diagcat((SX)1, (SX)2, (SX)3);
        DM A = ca.evalf(f.call(new[] { d, (SX)bn })[0]);
        DM B = ev(f, ca.diagcat((DM)1, (DM)2, (DM)3), bn);
        checkarray(A, B);
      }
    }
  }

  public void test_issue3489() {
    SX E_f = SX.sym("E_f"), E_r = SX.sym("E_r"), T = SX.sym("T"), k_f = SX.sym("k_f"), k_r = SX.sym("k_r");
    SX k0_f = SX.sym("k0_f"), k0_r = SX.sym("k0_r");
    SX x_1 = SX.sym("x_1"), x_2 = SX.sym("x_2"), x_3 = SX.sym("x_3"), x_4 = SX.sym("x_4");
    SX r_1 = SX.sym("r_1"), r_2 = SX.sym("r_2"), r_3 = SX.sym("r_3"), r_4 = SX.sym("r_4");
    SX y_1 = SX.sym("y_1"), y_2 = SX.sym("y_2"), y_3 = SX.sym("y_3"), y_4 = SX.sym("y_4");

    SX eq1 = k_f - k0_f - E_f + T;
    SX eq2 = k_r - k0_r - E_r + T;
    Func<SX, SX> rx = r => r - k_f + -x_1 + -x_2 - k_r + x_3 + x_4;
    SX eq3 = rx(r_1);
    SX eq4 = rx(r_2);
    SX eq5 = rx(r_3);
    SX eq6 = rx(r_4);
    SX eq7 = 1 * 0.2 + x_1 * -1.0 + r_1;
    SX eq8 = 1 * 0.8 + x_2 * -1.0 + r_2;
    SX eq9 = x_3 * -1.0 + r_3;
    SX eq10 = x_4 * -1.0 + r_4;
    SX eq11 = x_1 - y_1, eq12 = x_2 - y_2, eq13 = x_3 - y_3, eq14 = x_4 - y_4;

    var all_vars = new List<SX> { k_f, k_r, x_1, x_2, x_3, x_4, r_1, r_2, r_3, r_4, y_1, y_2, y_3, y_4 };
    var all_eq = new List<SX> { eq11, eq12, eq13, eq14, eq3, eq4, eq5, eq6, eq7, eq8, eq9, eq10, eq1, eq2 };

    SX x = ca.vertcat(all_vars);
    SX p = ca.vertcat(E_f, E_r, k0_f, k0_r, T);
    SX g = ca.vertcat(all_eq);

    DM A = ca.evalf(ca.jacobian(g, x));
    DM B = ca.evalf(ca.jacobian(g, p));
    MX.sym("A", A.sparsity());
    MX.sym("B", B.sparsity());

    var rf = new Dictionary<string, SX> { ["x"] = x, ["p"] = p, ["g"] = g };
    Function sim = ca.rootfinder("s", "fast_newton", rf, new Dictionary<string, object> { ["ad_weight_sp"] = 0 });
    Sparsity S1 = sim.jac_sparsity(0, 1);
    sim = ca.rootfinder("s", "fast_newton", rf, new Dictionary<string, object> { ["ad_weight_sp"] = 1 });
    Sparsity S2 = sim.jac_sparsity(0, 1);
    assertTrue(S1.is_equal(S2));

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
      var f0 = ca.Function("f", new[] { P }, new[] { ca.solve((MX)Ar, ca.mtimes(Br, P)) }, new Dictionary<string, object> { ["ad_weight_sp"] = 0 });
      S1 = f0.jac_sparsity(0, 0);
      var f1 = ca.Function("f", new[] { P }, new[] { ca.solve((MX)Ar, ca.mtimes(Br, P)) }, new Dictionary<string, object> { ["ad_weight_sp"] = 1 });
      S2 = f1.jac_sparsity(0, 0);
      assertTrue(S1.is_equal(S2));
      assertTrue(C.is_subset(S1));
    }
  }

  public void test_empty_system() {
    var L = ca.Linsol("L", "qr", ca.Sparsity(0, 0));
    DM r = L.solve(ca.DM(0, 0), ca.DM(0, 2));
    assertEquals((0L, 2L), (r.size1(), r.size2()));
  }

  public void test_issue2734() {
    DM example_input = ca.blockcat(new[] { new DM[] { 8, 1, 6 }, new DM[] { 3, 5, 7 }, new DM[] { 4, 9, 2 } });

    SX a = SX.sym("a", 3, 3);
    var f = ca.Function("f", new[] { a }, new[] { ca.pinv(a) });
    DM res1 = ev(f, example_input);

    MX am = MX.sym("am", 3, 3);
    var fm = ca.Function("f", new[] { am }, new[] { ca.pinv(am, "symbolicqr") });
    DM res2 = ev(fm, example_input);

    am = MX.sym("am", 3, 3);
    fm = ca.Function("f", new[] { am }, new[] { ca.pinv(am, "symbolicqr") });
    fm = fm.expand();
    DM res3 = ev(fm, example_input);

    checkarray(res1, res2);
    checkarray(res1, res3);
  }
}

}

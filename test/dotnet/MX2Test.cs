using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;

namespace CasADi.Tests {

/// <summary>Port of test/java MX2Test (test/python/mx.py, second part, from test_project_simplify)</summary>
public class MX2Test : CasadiTestCase {

  static readonly double[,] A15 = {
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

  static int codegenCounter = 0;

  static Dictionary<string, object> opts(params (string, object)[] kv) => kv.ToDictionary(p => p.Item1, p => p.Item2);

  static void checkCodegen(Function f, List<DM> inputs) => checkCodegen(f, inputs, opts(), 1e-15, "c89", false);

  static string readFile(string p) => File.ReadAllText(p).Replace("\r\n", "\n");

  static string checkCodegen(Function f, List<DM> inputs, Dictionary<string, object> o, double tol, string std, bool main) =>
    checkCodegen(f, inputs, o, tol, std, main, false);

  static string lastLib;
  static Function lastExternal;
  static string lastMainOut;

  static int run(string exe, IEnumerable<string> args, string stdin = null, string stdout = null) {
    var psi = new ProcessStartInfo(exe) { UseShellExecute = false, RedirectStandardInput = stdin != null,
                                          RedirectStandardOutput = stdout != null };
    foreach (string a in args) psi.ArgumentList.Add(a);
    using (var p = Process.Start(psi)) {
      if (stdin != null) {
        p.StandardInput.Write(File.ReadAllText(stdin));
        p.StandardInput.Close();
      }
      if (stdout != null) File.WriteAllText(stdout, p.StandardOutput.ReadToEnd());
      p.WaitForExit();
      return p.ExitCode;
    }
  }

  /// <summary>As check_codegen in helpers.py (only with RUN_SLOW); extras: with_jac_sparsity, with_forward and with_reverse</summary>
  static string checkCodegen(Function f, List<DM> inputs, Dictionary<string, object> opts0, double tol, string std, bool main, bool extras) {
    if (!RunSlow) return null;
    string tmp = Path.Combine(Path.GetTempPath(), "casadi_dotnet_codegen_" + Guid.NewGuid().ToString("N"));
    Directory.CreateDirectory(tmp);
    string dir = tmp + "/";
    string name = "codegen_dotnet_" + (codegenCounter++);
    var o = new Dictionary<string, object>(opts0);
    if (main) o["main"] = true;
    var cg = ca.CodeGenerator(name, o);
    cg.add(f, extras);
    if (extras) {
      cg.add(f.reverse(1), true);
      cg.add(f.forward(1), true);
    }
    cg.generate(dir);
    string libdir = GlobalOptions.getCasadiPath(), inc = GlobalOptions.getCasadiIncludePath();
    foreach (bool shared in main ? new[] { true, false } : new[] { true }) {
      var cmd = new List<string> { "-pedantic", "-std=" + std, "-fPIC", "-Wall", "-Werror", "-Wextra",
        "-I" + inc, "-Wno-unknown-pragmas", "-Wno-long-long", "-Wno-unused-parameter", "-O3", dir + name + ".c", "-o",
        dir + name + (shared ? ".so" : ""), "-L" + libdir, "-Wl,-rpath," + libdir };
      if (shared) cmd.Insert(3, "-shared");
      else cmd.Add("-lm");
      assertEquals(0, run("gcc", cmd), "compile gcc " + string.Join(" ", cmd));
    }
    string @out = null;
    if (main) {
      f.generate_in(dir + f.name() + "_in.txt", inputs);
      @out = dir + f.name() + "_out.txt";
      lastMainOut = @out;
      assertEquals(0, run(dir + name, new[] { f.name() }, dir + f.name() + "_in.txt", @out));
      List<DM> outputs = f.generate_out(@out);
      List<DM> r0 = f.call(inputs);
      for (int i = 0; i < r0.Count; i++) checkarray(r0[i], outputs[i], tol, "main output " + i);
    }
    Function f2 = ca.external(f.name(), dir + name + ".so");
    lastLib = dir + name + ".so";
    lastExternal = f2;
    List<DM> r = f.call(inputs), r2 = f2.call(inputs);
    for (int i = 0; i < r.Count; i++) checkarray(r[i], r2[i], tol, "codegen output " + i);
    return @out;
  }

  /// <summary>check_codegen(main=True) when generate_out is expected to fail (the program prints extra lines)</summary>
  static void checkCodegenMainParsingError(Function f, List<DM> inputs) {
    lastMainOut = null;
    if (!RunSlow) return;
    try {
      checkCodegen(f, inputs, opts(), 1e-15, "c89", true);
    } catch (SwigException e) {
      assertTrue((e.Message ?? "").Contains("Parsing error"), e.Message);
      return;
    }
    throw new AssertFailed("expected Parsing error");
  }

  /// <summary>As check_serialize in helpers.py</summary>
  static void checkSerialize(Function f, List<DM> inputs) {
    Function f2 = Function.deserialize(f.serialize(opts(("debug", true))));
    List<DM> r = f.call(inputs), r2 = f2.call(inputs);
    for (int i = 0; i < r.Count; i++) checkarray(r[i], r2[i], 1e-16, "serialize output " + i);
  }

  static Function hessianOld(Function f, long i, long j) {
    var @out = new List<string> {
      "hess:" + f.name_out(j) + ":" + f.name_in(i) + ":" + f.name_in(i),
      "grad:" + f.name_out(j) + ":" + f.name_in(i),
    };
    @out.AddRange(f.name_out());
    return f.factory(f.name() + "_hess", f.name_in(), @out);
  }

  static DM range(int a, int b) => Enumerable.Range(a, b - a).Select(i => (double)i).ToArray();

  static List<T> cat<T>(params IEnumerable<T>[] l) => l.SelectMany(e => e).ToList();

  static List<MX> L(params MX[] x) => x.ToList();

  public void test_project_simplify() {
    MX x = MX.sym("x", Sparsity.lower(3));
    DM.rng(1);
    DM X0 = DM.rand(x.sparsity());

    MX y = ca.project(x, Sparsity.upper(3));
    MX z = ca.project(y, ca.sparsify(ca.DM(new double[,] { { 0, 0, 0 }, { 1, 1, 1 }, { 1, 1, 1 } })).sparsity());

    Function f = ca.Function("f", new[] { x }, new[] { z });
    Function fs = f.transform();
    assertTrue(fs.n_nodes() > 3);

    checkfunction(f, fs, dms(X0));

    z = ca.project(y, ca.sparsify(ca.DM(new double[,] { { 1, 0, 0 }, { 0, 0, 0 }, { 0, 0, 1 } })).sparsity());

    f = ca.Function("f", new[] { x }, new[] { z });
    fs = f.transform();

    assertTrue(fs.n_nodes() != 0);

    fs.disp(true);

    checkfunction(f, fs, dms(X0));
  }

  public void test_repmat() {
    DM a = ca.DM(new double[,] { { 1, 2 }, { 3, 4 }, { 5, 6 } });
    checkarray(ca.repmat(a, 2, 3), ca.kron(DM.ones(2, 3), a));
  }

  public void test_transpose_bug() {
    DM A = ca.sparsify(ca.DM(A15));

    MX As = MX.sym("As", A.sparsity());

    Function f = ca.Function("f", new[] { As }, new[] { ca.densify(As.T()), ca.densify(As).T(), As.T(), As, ca.densify(As) });

    List<DM> f_out = f.call(dms(A));

    checkarray(f_out[0], A.T());
    checkarray(f_out[1], A.T());
    checkarray(f_out[2], A.T());
    checkarray(f_out[3], A);
    checkarray(f_out[4], A);
  }

  public void test_bizarre_bug() {
    Check.Assume(hasLinsol("csparse"), "no linsol csparse");
    DM A = ca.sparsify(ca.DM(A15));

    DM b = range(0, 15);

    MX Bs = MX.sym("Bs", b.sparsity());
    MX As = MX.sym("As", A.sparsity());

    MX Ast = As.T();

    Function r = ca.Function("r", new[] { As, Bs }, new[] { ca.solve(Ast, Bs, "csparse") });
    Function R = ca.Function("R", new[] { As, Bs }, new[] { ca.solve(ca.densify(Ast), Bs, "csparse") });

    DM r_out = r.call(dms(A, b))[0];
    DM R_out = R.call(dms(A, b))[0];

    r.sparsity_out(0).spy();
    R.sparsity_out(0).spy();

    checkarray(ca.mtimes(A.T(), R_out), b, 1e-9, "");
    checkarray(r_out, R_out);
  }

  public void test_depends_on() {
    MX a = MX.sym("a");
    MX b = MX.sym("b");

    assertTrue(ca.depends_on(ca.pow(a, 2), a));
    assertTrue(ca.depends_on(a, a));
    assertFalse(ca.depends_on(ca.MX(0), a));

    MX ab = ca.vertcat(a, b);
    assertTrue(ca.depends_on(ca.pow(a, 2), ab));
    assertTrue(ca.depends_on(a, ab));
    assertFalse(ca.depends_on(ca.MX(0), ab));
    assertTrue(ca.depends_on(ca.pow(b, 2), ab));
    assertTrue(ca.depends_on(b, ab));
    assertTrue(ca.depends_on(ca.pow(a, 2) + ca.pow(b, 2), ab));
    assertTrue(ca.depends_on(a + b, ab));
    MX zero = ca.MX(0);
    assertTrue(ca.depends_on(ca.vertcat(zero, a), a));
    assertTrue(ca.depends_on(ca.vertcat(a, zero), a));
    assertFalse(ca.depends_on(ca.vertcat(zero, b), a));
    assertFalse(ca.depends_on(ca.vertcat(b, zero), a));
    assertTrue(ca.depends_on(ca.vertcat(ca.pow(a, 2), ca.pow(b, 2)), ab));
    assertTrue(ca.depends_on(ca.vertcat(a, zero), ab));
    assertTrue(ca.depends_on(ca.vertcat(zero, b), ab));
    assertTrue(ca.depends_on(ca.vertcat(b, zero), ab));
    assertFalse(ca.depends_on(ca.vertcat(zero, zero), ab));
  }

  public void test_vertcat_simp() {
    MX x = MX.sym("x", 10);
    MX y = MX.sym("y");
    MX z = MX.sym("z");
    DM x_ = range(0, 10);
    DM y_ = ca.vertcat(20);
    DM z_ = ca.vertcat(30);

    Func<List<MX>, DM> evalvertcat = a => {
      Function f = ca.Function("f", new[] { x, y, z }, new[] { ca.vertcat(a) });
      return f.call(dms(x_, y_, z_))[0];
    };
    List<MX> xs = ca.vertsplit(x);
    List<MX> ly = L(y), lz = L(z);
    List<MX> xm = xs.GetRange(0, 9), xp = xs.GetRange(1, 9);
    DM x_m = x_[ca.Slice(0, -1)], x_p = x_[1..10];

    checkarray(evalvertcat(xs), x_);
    checkarray(evalvertcat(cat(xs, ly)), ca.vertcat(x_, y_));
    checkarray(evalvertcat(cat(lz, xs, ly, xs, lz)), ca.vertcat(z_, x_, y_, x_, z_));
    checkarray(evalvertcat(xm), x_m);
    checkarray(evalvertcat(cat(xm, ly)), ca.vertcat(x_m, y_));
    checkarray(evalvertcat(cat(lz, xm, ly, xm, lz)), ca.vertcat(z_, x_m, y_, x_m, z_));
    checkarray(evalvertcat(xp), x_p);
    checkarray(evalvertcat(cat(xp, ly)), ca.vertcat(x_p, y_));
    checkarray(evalvertcat(cat(lz, xp, ly, xp, lz)), ca.vertcat(z_, x_p, y_, x_p, z_));
    List<MX> g = cat(xs.GetRange(5, 5), xs.GetRange(0, 5));
    DM x5 = x_[5..10], x0 = x_[0..5];
    checkarray(evalvertcat(g), ca.vertcat(x5, x0));
    checkarray(evalvertcat(cat(g, ly)), ca.vertcat(x5, x0, y_));
    checkarray(evalvertcat(cat(lz, g, ly, g, lz)), ca.vertcat(z_, x5, x0, y_, x5, x0, z_));

    List<MX> w = ca.vertsplit(x, 2);
    var r = new List<MX>();
    foreach (MX i in w) r.AddRange(ca.vertsplit(i));

    checkarray(evalvertcat(r), x_);

    r = new List<MX>();
    foreach (MX i in w) r.AddRange(cat(ca.vertsplit(i), ly));
    Console.WriteLine("vertcat: [" + string.Join(", ", r) + "]");
    Console.WriteLine("result: " + ca.vertcat(r));

    r = new List<MX>();
    foreach (MX i in w) r.AddRange(ca.vertsplit(i));
    Console.WriteLine("vertcat: [" + string.Join(", ", r) + "]");
    Console.WriteLine("result: " + ca.vertcat(cat(r, ly)));

    assertTrue(ca.is_equal(ca.vertcat(ca.vertsplit(x)), x));
  }

  public void test_horzcat_simp() {
    MX x = MX.sym("x", 1, 10);
    MX y = MX.sym("y");
    MX z = MX.sym("z");
    DM x_ = range(0, 10).T();
    DM y_ = ca.vertcat(20);
    DM z_ = ca.vertcat(30);

    Func<List<MX>, DM> evalhorzcat = a => {
      Function f = ca.Function("f", new[] { x, y, z }, new[] { ca.horzcat(a) });
      return f.Call(x_, y_, z_);
    };
    List<MX> xs = ca.horzsplit(x);
    List<MX> ly = L(y), lz = L(z);
    List<MX> xm = xs.GetRange(0, 9), xp = xs.GetRange(1, 9);
    DM x_m = x_[0, ca.Slice(0, -1)], x_p = x_[0, 1..10];

    checkarray(evalhorzcat(xs), x_);
    checkarray(evalhorzcat(cat(xs, ly)), ca.horzcat(x_, y_));
    checkarray(evalhorzcat(cat(lz, xs, ly, xs, lz)), ca.horzcat(z_, x_, y_, x_, z_));
    checkarray(evalhorzcat(xm), x_m);
    checkarray(evalhorzcat(cat(xm, ly)), ca.horzcat(x_m, y_));
    checkarray(evalhorzcat(cat(lz, xm, ly, xm, lz)), ca.horzcat(z_, x_m, y_, x_m, z_));
    checkarray(evalhorzcat(xp), x_p);
    checkarray(evalhorzcat(cat(xp, ly)), ca.horzcat(x_p, y_));
    checkarray(evalhorzcat(cat(lz, xp, ly, xp, lz)), ca.horzcat(z_, x_p, y_, x_p, z_));
    List<MX> g = cat(xs.GetRange(5, 5), xs.GetRange(0, 5));
    DM x5 = x_[0, 5..10], x0 = x_[0, 0..5];
    checkarray(evalhorzcat(g), ca.horzcat(x5, x0));
    checkarray(evalhorzcat(cat(g, ly)), ca.horzcat(x5, x0, y_));
    checkarray(evalhorzcat(cat(lz, g, ly, g, lz)), ca.horzcat(z_, x5, x0, y_, x5, x0, z_));

    List<MX> w = ca.horzsplit(x, 2);
    var r = new List<MX>();
    foreach (MX i in w) r.AddRange(ca.horzsplit(i));

    checkarray(evalhorzcat(r), x_);

    r = new List<MX>();
    foreach (MX i in w) r.AddRange(cat(ca.horzsplit(i), ly));
    Console.WriteLine("vertcat: [" + string.Join(", ", r) + "]");
    Console.WriteLine("result: " + ca.horzcat(r));

    r = new List<MX>();
    foreach (MX i in w) r.AddRange(ca.horzsplit(i));
    Console.WriteLine("vertcat: [" + string.Join(", ", r) + "]");
    Console.WriteLine("result: " + ca.horzcat(cat(r, ly)));

    assertTrue(ca.is_equal(ca.horzcat(ca.horzsplit(x)), x));
  }

  public void test_vertsplit_simp() {
    var dvars = new List<MX>();
    for (int i = 0; i < 5; i++) dvars.Add(MX.sym("abcdefghijklm".Substring(i, 1)));
    double[] dvars_ = { 0, 1, 2, 3, 4 };

    MX zz = MX.sym("zz", 2);
    DM zz_ = ca.vertcat(11, 12);
    MX y = MX.sym("y");
    MX z = MX.sym("z");
    DM y_ = ca.vertcat(20);
    DM z_ = ca.vertcat(30);

    MX aa = MX.sym("aa", 5);
    DM aa_ = range(100, 105);

    Func<MX, long?, List<DM>> evalvertsplit = (a, incr) => {
      List<MX> sp = incr == null ? ca.vertsplit(a) : ca.vertsplit(a, incr.Value);
      Console.WriteLine("[" + string.Join(", ", sp) + "]");
      Function f = ca.Function("f", cat(dvars, L(y, z, zz, aa)), sp);
      var f_in = new List<DM>();
      for (int i = 0; i < 5; i++) f_in.Add(ca.vertcat(dvars_[i]));
      f_in.AddRange(new[] { y_, z_, zz_, aa_ });
      return f.call(f_in);
    };

    List<DM> s = evalvertsplit(ca.vertcat(cat(L(y), dvars, L(z))), null);
    checkarray(s[0], y_);
    for (int i = 0; i < 5; i++) checkarray(s[1 + i], dvars_[i]);
    checkarray(s[6], z_);

    s = evalvertsplit(ca.vertcat(cat(L(y), dvars, L(z))), 2L);

    checkarray(s[0], ca.vertcat(y_, ca.vertcat(dvars_[0])));
    checkarray(s[1], ca.vertcat(dvars_[1], dvars_[2]));
    checkarray(s[2], ca.vertcat(dvars_[3], dvars_[4]));
    checkarray(s[3], z_);

    s = evalvertsplit(ca.vertcat(y, zz, z, zz), 2L);

    checkarray(s[0], ca.vertcat(y_, zz_[0]));
    checkarray(s[1], ca.vertcat(zz_[1], z_));
    checkarray(s[2], zz_);

    s = evalvertsplit(ca.vertcat(y, zz, z, zz), 3L);

    checkarray(s[0], ca.vertcat(y_, zz_[0], zz_[1]));
    checkarray(s[1], ca.vertcat(z_, zz_[0], zz_[1]));

    s = evalvertsplit(ca.vertcat(zz, zz), 2L);
    checkarray(s[0], zz_);
    checkarray(s[1], zz_);

    s = evalvertsplit(ca.vertcat(cat(L(zz), dvars)), null);
    checkarray(s[0], zz_[0]);
    checkarray(s[1], zz_[1]);

    for (int i = 0; i < 5; i++) checkarray(s[2 + i], dvars_[i]);

    s = evalvertsplit(ca.vertcat(cat(dvars, L(aa))), 5L);
    checkarray(s[0], ca.DM(dvars_));
    checkarray(s[1], aa_);

    s = evalvertsplit(ca.vertcat(cat(dvars, L(aa))), 4L);
    checkarray(s[0], ca.vertcat(0, 1, 2, 3));
    checkarray(s[1], ca.vertcat(4, 100, 101, 102));
    checkarray(s[2], ca.vertcat(103, 104));

    s = evalvertsplit(ca.vertcat(cat(dvars, L(aa))), 6L);
    checkarray(s[0], ca.vertcat(0, 1, 2, 3, 4, 100));
    checkarray(s[1], ca.vertcat(101, 102, 103, 104));

    for (int i = 0; i < 5; i++) assertTrue(ca.is_equal(ca.vertsplit(ca.vertcat(dvars))[i], dvars[i]));
  }

  public void test_horzsplit_simp() {
    var dvars = new List<MX>();
    for (int i = 0; i < 5; i++) dvars.Add(MX.sym("abcdefghijklm".Substring(i, 1)));
    double[] dvars_ = { 0, 1, 2, 3, 4 };

    MX zz = MX.sym("zz", 1, 2);
    DM zz_ = ca.vertcat(11, 12).T();
    MX y = MX.sym("y");
    MX z = MX.sym("z");
    DM y_ = ca.vertcat(20);
    DM z_ = ca.vertcat(30);

    MX aa = MX.sym("aa", 1, 5);
    DM aa_ = range(100, 105);

    Func<MX, long?, List<DM>> evalhorzsplit = (a, incr) => {
      List<MX> sp = incr == null ? ca.horzsplit(a) : ca.horzsplit(a, incr.Value);
      Console.WriteLine("[" + string.Join(", ", sp) + "]");
      Function f = ca.Function("f", cat(dvars, L(y, z, zz, aa)), sp);
      var f_in = new List<DM>();
      for (int i = 0; i < 5; i++) f_in.Add(ca.vertcat(dvars_[i]));
      f_in.AddRange(new[] { y_, z_, zz_, aa_ });
      return f.call(f_in);
    };

    List<DM> s = evalhorzsplit(ca.horzcat(cat(L(y), dvars, L(z))), null);
    checkarray(s[0], y_);
    for (int i = 0; i < 5; i++) checkarray(s[1 + i], dvars_[i]);
    checkarray(s[6], z_);

    s = evalhorzsplit(ca.horzcat(cat(L(y), dvars, L(z))), 2L);

    checkarray(s[0], ca.vertcat(y_, ca.vertcat(dvars_[0])).T());
    checkarray(s[1], ca.vertcat(dvars_[1], dvars_[2]).T());
    checkarray(s[2], ca.vertcat(dvars_[3], dvars_[4]).T());
    checkarray(s[3], z_.T());

    s = evalhorzsplit(ca.horzcat(y, zz, z, zz), 2L);

    checkarray(s[0], ca.vertcat(y_, zz_[0, 0]).T());
    checkarray(s[1], ca.vertcat(zz_[0, 1], z_).T());
    checkarray(s[2], zz_);

    s = evalhorzsplit(ca.horzcat(y, zz, z, zz), 3L);

    checkarray(s[0], ca.vertcat(y_, zz_[0, 0], zz_[0, 1]).T());
    checkarray(s[1], ca.vertcat(z_, zz_[0, 0], zz_[0, 1]).T());

    s = evalhorzsplit(ca.horzcat(zz, zz), 2L);
    checkarray(s[0], zz_);
    checkarray(s[1], zz_);

    s = evalhorzsplit(ca.horzcat(cat(L(zz), dvars)), null);
    checkarray(s[0], zz_[0, 0]);
    checkarray(s[1], zz_[0, 1]);

    for (int i = 0; i < 5; i++) checkarray(s[2 + i], dvars_[i]);

    s = evalhorzsplit(ca.horzcat(cat(dvars, L(aa))), 5L);
    checkarray(s[0], ca.DM(dvars_).T());
    checkarray(s[1], aa_.T());

    s = evalhorzsplit(ca.horzcat(cat(dvars, L(aa))), 4L);
    checkarray(s[0], ca.vertcat(0, 1, 2, 3).T());
    checkarray(s[1], ca.vertcat(4, 100, 101, 102).T());
    checkarray(s[2], ca.vertcat(103, 104).T());

    s = evalhorzsplit(ca.horzcat(cat(dvars, L(aa))), 6L);
    checkarray(s[0], ca.vertcat(0, 1, 2, 3, 4, 100).T());
    checkarray(s[1], ca.vertcat(101, 102, 103, 104).T());

    for (int i = 0; i < 5; i++) assertTrue(ca.is_equal(ca.horzsplit(ca.horzcat(dvars))[i], dvars[i]));
  }

  public void test_vertsplit_derivative() {
    MX m = MX.sym("X", 10);

    Function f = ca.Function("f", new[] { m }, new[] { ca.vertsplit(m)[0] });

    f.reverse(1);
  }

  public void test_MX_const_sp() {
    MX x = MX.sym("x", 4, 1);

    Sparsity sp = Sparsity.triplet(3, 3, new long[] { 0, 1, 2, 2 }, new long[] { 0, 0, 1, 2 });

    Function f = ca.Function("f", new[] { x }, new[] { x.get_nz(false, ca.DM(sp, range(0, (int)sp.nnz()))) });

    Function g = ca.Function("g", new[] { x }, new[] { ca.MX(sp, x) });

    List<DM> f_in = dms(range(1, 5));

    checkfunction(f, g, f_in);
    checkCodegen(f, f_in);
    checkCodegen(g, f_in);
  }

  public void test_reshape_sp() {
    MX x = MX.sym("x", 4, 1);

    Function f = ca.Function("f", new[] { x }, new[] { ca.reshape(x, (2, 2)) });

    SX sx = SX.sym("x", 4, 1);

    Function g = ca.Function("g", new[] { sx }, new[] { ca.reshape(sx, (2, 2)) });

    List<DM> f_in = dms(range(1, 5));

    checkfunction(f, g, f_in);
    checkCodegen(f, f_in);
    checkCodegen(g, f_in);
  }

  public void test_issue1041() {
    MX x = MX.sym("x", 2);

    MX y = ca.vertsplit(x, new long[] { 0, 1, 2 })[1];

    Function f = ca.Function("f", new[] { x }, new[] { y });

    hessianOld(f, 0, 0);
  }

  public void test_bug_1042() {
    MX x = MX.sym("x", 2, 1);

    Function mf = ca.Function("mf", new[] { x }, new[] { x * x[0, 0] });

    Function mfunction = mf.expand("expand_" + mf.name());

    Function mfg = mf.reverse(1);

    Function mfunctiong = mfunction.reverse(1);

    List<DM> f_in = dms(ca.vertcat(0), ca.vertcat(5), ca.vertcat(1, 2));

    checkfunction(mfg, mfunctiong, f_in);
  }

  public void test_bug_1042bis() {
    MX x = MX.sym("x", 2, 1);
    MX a = MX.sym("ax", 2, 1);
    MX i1 = x[0, 0];
    MX z = i1 * x;
    MX i3 = i1 * a;
    i3 = ca.dot(x, a);
    Function d = ca.Function("d", new[] { x, a }, new[] { z, i3 });

    Function dx = d.expand("expand_" + d.name());
    List<DM> dx_in = dms(ca.vertcat(1, 2), ca.vertcat(3, 4));

    checkfunction(d, dx, dx_in);
  }

  public void test_bug_1042tris() {
    MX x = MX.sym("x", 2, 1);
    MX a = MX.sym("ax", 2, 1);
    Function d = ca.Function("d", new[] { x, a }, new[] { ca.dot(x, a) });

    Function dx = d.expand("expand_" + d.name());
    List<DM> dx_in = dms(ca.vertcat(1, 2), ca.vertcat(3, 4));

    checkfunction(d, dx, dx_in);
  }

  public void test_bug_1046() {
    MX x = MX.sym("x", 1, 1);
    MX y = MX.sym("y", 1, 1);
    MX z = ca.jacobian(x, y);

    assertEquals(0L, z.nnz());
  }

  public void test_expand_free() {
    MX x = MX.sym("x");
    MX y = MX.sym("y");
    assertInException("free", () => ca.Function("f", new[] { x }, new[] { x + y }));

    Function f = ca.Function("f", new[] { x }, new[] { x + y }, opts(("allow_free", true)));
    Function g = ca.Function("g", new[] { x, y }, new[] { f.Call(x) });
    g.expand();
    checkfunction(g, g, dms(ca.vertcat(1), ca.vertcat(2)));
  }

  /// <summary>vertsplit/horzsplit/concatenations through dynamic dispatch, one body for MX, SX and DM</summary>
  static void singularcat(Func<long, long, dynamic> zeros) {
    long[] o = { 0, 5, 10 };
    dynamic x0 = zeros(10, 0);
    dynamic x1s = ca.vertsplit(x0, o);
    foreach (dynamic x in x1s) shape(x, 5, 0);
    shape(ca.vertcat(x1s), 10, 0);
    shape(ca.vertcat(catd(L1(zeros(0, 0)), x1s, L1(zeros(0, 0)))), 10, 0);
    x0 = zeros(0, 1);
    shape(ca.vertcat(catd(L1(x0), L1(zeros(0, 0)), L1(x0))), 0, 1);

    x0 = zeros(0, 10);
    x1s = ca.horzsplit(x0, o);
    foreach (dynamic x in x1s) shape(x, 0, 5);
    shape(ca.horzcat(x1s), 0, 10);
    shape(ca.horzcat(catd(L1(zeros(0, 0)), x1s, L1(zeros(0, 0)))), 0, 10);
    x0 = zeros(1, 0);
    shape(ca.horzcat(catd(L1(x0), L1(zeros(0, 0)), L1(x0))), 1, 0);

    x1s = ca.vertsplit(zeros(10, 0), o);
    dynamic x1st = ca.horzsplit(zeros(0, 10), o);
    dynamic e = L1(zeros(0, 0));
    shape(ca.diagcat(x1s), 10, 0);
    shape(ca.diagcat(catd(e, x1s, e)), 10, 0);
    shape(ca.diagcat(x1st), 0, 10);
    shape(ca.diagcat(catd(e, x1st, e)), 0, 10);
    shape(ca.diagcat(catd(x1s, x1st)), 10, 10);
    shape(ca.diagcat(catd(e, x1s, x1st, e)), 10, 10);
  }

  /// <summary>One-element list of the runtime type of x</summary>
  static object L1(object x) {
    if (x is MX m) return new List<MX> { m };
    if (x is SX s) return new List<SX> { s };
    return new List<DM> { (DM)x };
  }

  /// <summary>Concatenation of lists of the same element type, typed after the first element</summary>
  static dynamic catd(params dynamic[] lists) {
    var items = lists.SelectMany(l => (IEnumerable<object>)l).ToList();
    if (items.Count > 0 && items[0] is MX) return items.Cast<MX>().ToList();
    if (items.Count > 0 && items[0] is SX) return items.Cast<SX>().ToList();
    return items.Cast<DM>().ToList();
  }

  static void shape(dynamic x, long n, long m) {
    assertEquals(n, (long)x.size1());
    assertEquals(m, (long)x.size2());
  }

  public void test_singularcat() {
    singularcat((n, m) => MX.zeros(n, m));
    singularcat((n, m) => SX.zeros(n, m));
    singularcat((n, m) => DM.zeros(n, m));
  }

  public void test_exprjac() {
    var funmx = new List<Func<MX, MX>> { x => x[0] * x[1], x => x[0] * ca.sin(x[1]) };
    var funsx = new List<Func<SX, SX>> { x => x[0] * x[1], x => x[0] * ca.sin(x[1]) };
    var opmx = new List<Func<MX, MX, MX>> { (f, x) => ca.gradient(f, x), (f, x) => ca.jacobian(f, x), (f, x) => ca.hessian(f, x).Item1 };
    var opsx = new List<Func<SX, SX, SX>> { (f, x) => ca.gradient(f, x), (f, x) => ca.jacobian(f, x), (f, x) => ca.hessian(f, x).Item1 };
    for (int k = 0; k < 2; k++) {
      for (int o = 0; o < 3; o++) {
        MX x = MX.sym("x", 2);
        Console.WriteLine(funmx[k](x));
        Console.WriteLine(opmx[o](funmx[k](x), x));
        Function f = ca.Function("f", new[] { x }, new[] { opmx[o](funmx[k](x), x) });

        SX xs = SX.sym("x", 2);
        Function fr = ca.Function("fr", new[] { xs }, new[] { opsx[o](funsx[k](xs), xs) });

        checkfunction(f, fr, dms(ca.vertcat(0)));
        checkCodegen(f, dms(ca.vertcat(0)));
      }
    }
  }

  static long prod(long[] d) => d.Aggregate(1L, (a, b) => a * b);

  /// <summary>Dimensions of the summed (negative) indices, in order of first appearance</summary>
  static (List<long> keys, Dictionary<long, long> map) dimMap(long[] dim_a, long[] dim_b, long[] ind_a, long[] ind_b) {
    var keys = new List<long>();
    var m = new Dictionary<long, long>();
    long[][] dims = { dim_a, dim_b }, inds = { ind_a, ind_b };
    for (int k = 0; k < 2; k++)
      for (int i = 0; i < inds[k].Length; i++) {
        long j = inds[k][i];
        if (j < 0) {
          if (m.ContainsKey(j)) assertEquals(m[j], dims[k][i]);
          else {
            m[j] = dims[k][i];
            keys.Add(j);
          }
        }
      }
    return (keys, m);
  }

  static long linIndex(long[] dims, long[] ind) {
    long ret = 0, cumprod = 1;
    for (int i = 0; i < dims.Length; i++) {
      ret += ind[i] * cumprod;
      cumprod *= dims[i];
    }
    return ret;
  }

  static long[] subs(long[] e, List<long> old, long[] n) => e.Select(v => old.Contains(v) ? n[old.IndexOf(v)] : v).ToArray();

  static MX myEinstein(MX A, MX B, MX C, long[] dim_a, long[] dim_b, long[] dim_c, long[] ind_a, long[] ind_b, long[] ind_c) {
    MX R = ca.MX(C);
    var (dk, d) = dimMap(dim_a, dim_b, ind_a, ind_b);
    var val = new long[dk.Count];
    while (true) {
      long ai = linIndex(dim_a, subs(ind_a, dk, val));
      long bi = linIndex(dim_b, subs(ind_b, dk, val));
      long ci = linIndex(dim_c, subs(ind_c, dk, val));
      R[ci] = R[ci] + A[ai] * B[bi];
      int k = val.Length - 1;
      while (k >= 0 && ++val[k] == d[dk[k]]) val[k--] = 0;
      if (k < 0) break;
    }
    return R;
  }

  void einsteinTests(long[] dim_a, long[] dim_b, long[] dim_c, long[] ind_a, long[] ind_b, long[] ind_c) {
    MX A = MX.sym("A", prod(dim_a));
    MX B = MX.sym("B", prod(dim_b));
    MX C = MX.sym("C", prod(dim_c), 1);

    SX A_sx = SX.sym("A", prod(dim_a));
    SX B_sx = SX.sym("B", prod(dim_b));
    SX C_sx = SX.sym("C", prod(dim_c), 1);

    List<DM> @in = dms(randDM(A.size1(), 1), randDM(B.size1(), 1), randDM(C.size1(), 1));

    MX @out = ca.einstein(A, B, C, dim_a, dim_b, dim_c, ind_a, ind_b, ind_c);
    SX out_sx = ca.einstein(A_sx, B_sx, C_sx, dim_a, dim_b, dim_c, ind_a, ind_b, ind_c);
    Function f = ca.Function("f", new[] { A, B, C }, new[] { @out }, opts(("ad_weight_sp", 0)));
    Function frev = ca.Function("f", new[] { A, B, C }, new[] { @out }, opts(("ad_weight_sp", 1)));
    Function fr = ca.Function("fr", new[] { A, B, C }, new[] { myEinstein(A, B, C, dim_a, dim_b, dim_c, ind_a, ind_b, ind_c) });
    Function f_sx = ca.Function("f", new[] { A_sx, B_sx, C_sx }, new[] { out_sx });

    Function fsx = f.expand();
    checkfunction(f, fr, @in);
    checkfunction(fsx, fr, @in);
    checkCodegen(f, @in);
    checkSerialize(f, @in);
    checkfunction(fsx, f_sx, @in);

    for (int i = 0; i < 3; i++) {
      assertTrue(f.jac_sparsity(0, i).is_equal(fsx.jac_sparsity(0, i)));
      assertTrue(frev.jac_sparsity(0, i).is_equal(fsx.jac_sparsity(0, i)));
    }
  }

  static void permutations(List<long> rest, List<long> prefix, List<long[]> @out) {
    if (rest.Count == 0) @out.Add(prefix.ToArray());
    for (int i = 0; i < rest.Count; i++) {
      var r = new List<long>(rest);
      var p = new List<long>(prefix) { r[i] };
      r.RemoveAt(i);
      permutations(r, p, @out);
    }
  }

  public void test_einstein() {
    long[] dim_b = { 3, 3 };
    long[] dim_a = { 3, 3 };
    long[][] ind = { new long[] { 0, -1 }, new long[] { 1, -1 }, new long[] { 0, 1 }, new long[] { -1, -2 }, new long[] { -2, -1 } };

    foreach (long[] ind_a in ind) {
      foreach (long[] ind_b in ind) {
        var fa = new SortedSet<long>(ind_a.Where(i => i < 0));
        var fb = new SortedSet<long>(ind_b.Where(i => i < 0));
        var free = fa.Where(i => !fb.Contains(i)).Concat(fb.Where(i => !fa.Contains(i))).ToList();
        var perms = new List<long[]>();
        permutations(free, new List<long>(), perms);
        foreach (long[] ind_c in perms) {
          var (_, m) = dimMap(dim_a, dim_b, ind_a, ind_b);
          long[] dim_c = ind_c.Select(i => m[i]).ToArray();
          einsteinTests(dim_a, dim_b, dim_c, ind_a, ind_b, ind_c);
        }
      }
    }

    einsteinTests(new long[] { 2, 4, 3 }, new long[] { 2, 5, 3 }, new long[] { 5, 4 }, new long[] { -1, -2, -3 },
                  new long[] { -1, -4, -3 }, new long[] { -4, -2 });
  }

  static object timesAny(object a, object b) {
    if (a is DM da && b is DM db) return da * db;
    MX x = a is DM xa ? ca.MX(xa) : (MX)a, y = b is DM yb ? ca.MX(yb) : (MX)b;
    return x * y;
  }

  static long nnzAny(object a) => a is DM d ? d.nnz() : ((MX)a).nnz();

  public void test_sparsity_operation() {
    var L = new List<object> { ca.MX(ca.Sparsity(1, 1)), ca.MX(ca.Sparsity(2, 1)), MX.sym("x", 1, 1),
      MX.sym("x", ca.Sparsity(1, 1)), ca.DM(1), ca.DM(ca.Sparsity(1, 1), 1), ca.DM(ca.Sparsity(2, 1), 1),
      ca.DM(Sparsity.dense(2, 1), 1) };

    foreach (object a in L) {
      foreach (object b in L) {
        object c = timesAny(a, b);

        if (nnzAny(a) == 0 || nnzAny(b) == 0) assertEquals(0L, nnzAny(c));
        else assertTrue(nnzAny(c) > 0);
      }
    }
  }

  public void test_solve() {
    Check.Assume(hasLinsol("lapackqr"), "no linsol lapackqr");
    Console.WriteLine(123);
    DM A = DM.rand(3, 3);
    DM B = DM.rand(3, 50);

    DM C = ca.solve(A, B, "lapackqr", opts(("max_nrhs", 50)));
    DM C1 = ca.solve(A, B, "lapackqr", opts(("max_nrhs", 20)));
    DM C2 = ca.solve(A, B, "lapackqr");

    checkarray(C, C1);
    checkarray(C1, C2);
  }

  public void test_sparse_lt() {
    MX x = MX.sym("x", Sparsity.lower(5));
    assertEquals(5L * 6 / 2, (x > 0).nnz());
    assertEquals(5L * 5, (x >= 0).nnz());
  }

  public void test_inv() {
    DM Av = DM.rand(3, 3);
    {
      SX A = SX.sym("x", 3, 3);
      Function f = ca.Function("f", new[] { A }, new[] { ca.inv(A), ca.SX(ca.inv(Av)), ca.mpower(A, -1), ca.SX(ca.mpower(Av, -1)) });
      foreach (DM o in f.call(dms(Av))) checkarray(ca.mtimes(o, Av), DM.eye(3), 1e-9, "");
    }
    {
      MX A = MX.sym("x", 3, 3);
      Function f = ca.Function("f", new[] { A }, new[] { ca.inv(A), ca.MX(ca.inv(Av)), ca.mpower(A, -1), ca.MX(ca.mpower(Av, -1)) });
      foreach (DM o in f.call(dms(Av))) checkarray(ca.mtimes(o, Av), DM.eye(3), 1e-9, "");
    }
  }

  /// <summary>numpy.interp: piecewise linear, clamped to the end values</summary>
  static DM npInterp(double[] xq, double[] x, double[] v) {
    var r = new double[xq.Length];
    for (int k = 0; k < xq.Length; k++) {
      double q = xq[k];
      if (q <= x[0]) r[k] = v[0];
      else if (q >= x[^1]) r[k] = v[^1];
      else {
        int i = 0;
        while (q > x[i + 1]) i++;
        r[k] = v[i] + (v[i + 1] - v[i]) * (q - x[i]) / (x[i + 1] - x[i]);
      }
    }
    return r;
  }

  public void test_interp1d() {
    double[] v = { 7, 3, 4, -3 };
    MX vs = MX.sym("x", 4, 1);

    double[] xq = { 10, 0, 1, 2, 4, 8, 7, 5, 1.5 };
    double[] x = { 1, 2, 4, 8 };
    Function F = ca.Function("f", new[] { vs }, new[] { ca.interp1d(x, vs, xq) });

    checkarray(F.Call(v), npInterp(xq, x, v));

    F = ca.Function("f", new[] { vs }, new[] { ca.interp1d(x, vs, xq, "floor") });

    checkarray(F.Call(v), ca.vertcat(-3, 7, 7, 3, 4, -3, 4, 4, 7));

    F = ca.Function("f", new[] { vs }, new[] { ca.interp1d(x, vs, xq, "ceil") });

    checkarray(F.Call(v), ca.vertcat(-3, 7, 7, 3, 4, -3, -3, -3, 3));

    xq = new double[] { 10, 0, 1, 2, 3, 4, 3.5, 3.25, 1.5 };
    x = new double[] { 1, 2, 3, 4 };
    F = ca.Function("f", new[] { vs }, new[] { ca.interp1d(x, vs, xq, "linear", true) });

    checkarray(F.Call(v), npInterp(xq, x, v));
  }

  public void test_bilin_etc() {
    MX x = MX.sym("x", 3, 3);
    MX y = MX.sym("y", 3, 1);
    MX z = MX.sym("z", 3, 1);

    List<DM> @in = dms(DM.rand(3, 3), DM.rand(3, 1), DM.rand(3, 1));

    MX[][] es = { new[] { ca.bilin(x, y, z), ca.mtimes(ca.mtimes(y.T(), x), z) },
                  new[] { ca.rank1(x, 0.3, y, z), x + ca.mtimes(y, z.T()) * 0.3 } };
    foreach (MX[] e in es) {
      Function f = ca.Function("f", new[] { x, y, z }, new[] { e[0] });
      Function fref = ca.Function("fref", new[] { x, y, z }, new[] { e[1] });
      Function f_sx = f.expand();
      checkfunction(f, fref, @in);
      checkfunction(f_sx, fref, @in);
      checkCodegen(f, @in);
    }
  }

  public void test_bilin_short() {
    List<DM> @in = dms(DM.rand(3, 3), DM.rand(3, 1));
    {
      SX x = SX.sym("x", 3, 3), y = SX.sym("y", 3, 1);
      Function f = ca.Function("f", new[] { x, y }, new[] { ca.bilin(x, y) });
      Function fref = ca.Function("f", new[] { x, y }, new[] { ca.bilin(x, y, y) });
      checkfunction(f, fref, @in);
      checkCodegen(f, @in);
      checkCodegen(fref, @in);
    }
    {
      MX x = MX.sym("x", 3, 3), y = MX.sym("y", 3, 1);
      Function f = ca.Function("f", new[] { x, y }, new[] { ca.bilin(x, y) });
      Function fref = ca.Function("f", new[] { x, y }, new[] { ca.bilin(x, y, y) });
      checkfunction(f, fref, @in);
      checkCodegen(f, @in);
      checkCodegen(fref, @in);
    }
  }

  public void test_det_shape() {
    MX X = MX.sym("x", 2, 3);
    Check.Throws<SwigException>(() => ca.det(X));
    ca.det(MX.sym("x", 3, 3));
  }

  public void test_det() {
    MX X = MX.sym("x", 3, 3);
    SX x = SX.sym("x", 3, 3);

    List<DM> @in = dms(DM.rand(3, 3));

    Function f = ca.Function("f", new[] { x }, new[] { ca.det(x) });
    Function F = ca.Function("F", new[] { X }, new[] { ca.det(X) });
    checkfunction(f, F, @in);
    checkCodegen(F, @in);
  }

  /// <summary>Determinant by Gaussian elimination with partial pivoting (reference for numpy.linalg.det)</summary>
  static double detRef(double[,] a0) {
    int n = a0.GetLength(0);
    var a = (double[,])a0.Clone();
    double d = 1;
    for (int k = 0; k < n; k++) {
      int p = k;
      for (int i = k + 1; i < n; i++) if (Math.Abs(a[i, k]) > Math.Abs(a[p, k])) p = i;
      if (a[p, k] == 0) return 0;
      if (p != k) {
        for (int j = 0; j < n; j++) (a[p, j], a[k, j]) = (a[k, j], a[p, j]);
        d = -d;
      }
      d *= a[k, k];
      for (int i = k + 1; i < n; i++) {
        double l = a[i, k] / a[k, k];
        for (int j = k; j < n; j++) a[i, j] -= l * a[k, j];
      }
    }
    return d;
  }

  double[,] randArray(int n, int m) {
    var a = new double[n, m];
    for (int i = 0; i < n; i++)
      for (int j = 0; j < m; j++) a[i, j] = random.NextDouble();
    return a;
  }

  public void test_det_solvers() {
    var solvers = new List<string> { "qr" };
    if (hasLinsol("csparse")) solvers.Add("csparse");

    var dense_cases = new List<double[,]>();
    for (int n = 1; n <= 5; n++) dense_cases.Add(randArray(n, n));

    double[,] S = randArray(5, 5);
    S[0, 2] = S[0, 4] = S[2, 0] = S[3, 1] = S[4, 2] = 0;
    for (int i = 0; i < 5; i++) S[i, i] += 3;

    foreach (double[,] A in dense_cases) {
      int n = A.GetLength(0);
      DM @ref = detRef(A);
      MX X = MX.sym("x", n, n);
      SX x = SX.sym("x", n, n);

      checkarray(ca.det((DM)A), @ref, 1e-9, "");
      checkarray(ca.evalf(ca.det(ca.SX((DM)A))), @ref, 1e-9, "");

      SX xs = SX.sym("xs", n, n);
      checkarray(ca.evalf(ca.substitute(ca.det(xs, "symbolicqr"), xs, ca.SX((DM)A))), @ref, 1e-9, "");

      Function f = ca.Function("f", new[] { x }, new[] { ca.det(x) });
      foreach (string ls in solvers) {
        Function F = ca.Function("F", new[] { X }, new[] { ca.det(X, ls) });
        checkfunction(f, F, dms(A), 1e-8);
        checkarray(ca.det((DM)A, ls), @ref, 1e-9, "");
        if (ls == "qr") checkCodegen(F, dms(A));
      }

      checkarray(ca.det((DM)A, "qr", opts(("eps", 1e-13))), @ref, 1e-9, "");
    }

    DM refS = detRef(S);
    checkarray(ca.det((DM)S), refS, 1e-9, "");
    foreach (string ls in solvers) checkarray(ca.det(ca.sparsify((DM)S), ls), refS, 1e-9, "");

    assertInException("does not provide a symbolic determinant", () => ca.det(SX.sym("x", 2, 2), "qr"));

    if (hasLinsol("csparse"))
      assertInException("does not provide a symbolic determinant", () => ca.det(SX.sym("x", 2, 2), "csparse"));
  }

  static (DM, long) sxdet(double[,] Anp) {
    DM Asp = ca.sparsify(Anp);
    SX x = SX.sym("x", Asp.sparsity());
    Function f = ca.Function("f", new[] { x }, new[] { ca.det(x, "symbolicqr") });
    return (f.call(dms(Asp))[0], f.n_nodes());
  }

  public void test_det_symbolicqr_btf() {
    int n = 8;
    var B = new double[n, n];
    double[,] b1 = randArray(4, 4), b2 = randArray(4, 4);
    for (int i = 0; i < 4; i++)
      for (int j = 0; j < 4; j++) {
        B[i, j] = b1[i, j];
        B[4 + i, 4 + j] = b2[i, j];
      }
    var (rb, nodes_b) = sxdet(B);
    checkarray(rb, detRef(B), 1e-9, "");

    SX xd = SX.sym("x", n, n);
    long nodes_d = ca.Function("f", new[] { xd }, new[] { ca.det(xd, "symbolicqr") }).n_nodes();
    assertTrue(nodes_b < nodes_d);

    double[,] L = randArray(6, 6);
    for (int i = 0; i < 6; i++) L[i, i] += 3;
    for (int i = 0; i < 3; i++)
      for (int j = 3; j < 6; j++) L[i, j] = 0;
    DM Lsp = ca.sparsify(L);
    SX xL = SX.sym("x", Lsp.sparsity());
    SX dL = ca.det(xL, "symbolicqr");
    checkarray(ca.evalf(ca.substitute(dL, xL, ca.SX(Lsp))), detRef(L), 1e-9, "");
    List<long> rows = Lsp.sparsity().row(), cols = Lsp.sparsity().get_col();
    var coupling = new List<double>();
    for (int k = 0; k < Lsp.nnz(); k++) if (rows[k] >= 3 && cols[k] < 3) coupling.Add(k);
    assertTrue(coupling.Count > 0);
    DM cidx = coupling.ToArray();
    assertEquals(0L, ca.jacobian(dL, xL).get(ca.Slice(), cidx).nnz());

    double[,] U = randArray(5, 5);
    for (int i = 0; i < 5; i++) {
      U[i, i] += 1;
      for (int j = 0; j < i; j++) U[i, j] = 0;
    }
    var (rU, nU) = sxdet(U);
    checkarray(rU, detRef(U), 1e-9, "");
    assertTrue(nU < 30);

    double[,] Sg = randArray(4, 4);
    for (int j = 0; j < 4; j++) Sg[2, j] = 0;
    checkarray(sxdet(Sg).Item1, ca.DM(0), 1e-12, "");
  }

  public void test_mtimes_mismatch_segfault() {
    assertInException("incompatible dimensions", () => ca.mtimes(ca.DM(Sparsity.lower(5)), MX.sym("x", 100)));
  }

  public void test_monitor() {
    MX x = MX.sym("x");
    MX y = ca.sqrt(x.monitor("hey"));

    Function f = ca.Function("f", new[] { x }, new[] { y });
    string @out = captureStdout(() => f.Call(3));
    Console.WriteLine(@out);

    string @ref = "hey:\n3.0000000000000000e+00\n";
    assertEquals(@ref, @out);

    if (RunSlow) {
      checkCodegenMainParsingError(f, dms(ca.vertcat(3)));
      string out2 = readFile(lastMainOut);
      Console.WriteLine(out2);
      assertTrue(out2.Contains(@ref));
    }

    MX x2 = MX.sym("x", 2);
    MX y2 = ca.sqrt(x2.monitor("hey"));

    Function f2 = ca.Function("f", new[] { x2 }, new[] { y2 });
    @out = captureStdout(() => f2.Call(3));
    Console.WriteLine(@out);
    @ref = "hey:\n2x1: [3.0000000000000000e+00, 3.0000000000000000e+00]\n";
    assertEquals(@ref, @out);

    if (RunSlow) {
      checkCodegenMainParsingError(f2, dms(ca.vertcat(3)));
      string out2 = readFile(lastMainOut);
      Console.WriteLine(out2);
      assertTrue(out2.Contains(@ref));
    }
  }

  static void removeGlob(string dir, string prefix) {
    if (!Directory.Exists(dir)) return;
    foreach (string f in Directory.GetFiles(dir)) {
      string n = Path.GetFileName(f);
      if (n.StartsWith(prefix) && n.EndsWith(".mtx")) File.Delete(f);
    }
  }

  static void rmtree(string d) {
    if (Directory.Exists(d)) Directory.Delete(d, true);
  }

  public void test_dump() {
    removeGlob(".", "myvar.");

    MX x = MX.sym("x");
    MX y = ca.sqrt(x.dump("myvar"));

    Function f = ca.Function("f", new[] { x }, new[] { y });

    f.Call(9);

    assertTrue(File.Exists("myvar.000000.mtx"));
    DM val = DM.from_file("myvar.000000.mtx");
    checkarray(val, 9);

    f.Call(16);
    assertTrue(File.Exists("myvar.000001.mtx"));
    val = DM.from_file("myvar.000001.mtx");
    checkarray(val, 16);

    removeGlob(".", "myvar.");

    x = MX.sym("x", Sparsity.lower(2));
    y = x.dump("mysp") * 2;
    f = ca.Function("f", new[] { x }, new[] { y });

    DM inp = ca.sparsify(ca.DM(new double[,] { { 1, 0 }, { 3, 4 } }));
    f.Call(inp);

    assertTrue(File.Exists("mysp.000000.mtx"));
    val = DM.from_file("mysp.000000.mtx");
    checkarray(val, inp);

    removeGlob(".", "mysp.");

    string dump_dir = "dump_test_dir";
    Directory.CreateDirectory(dump_dir);
    x = MX.sym("x");
    y = x.dump("myvar", opts(("dir", dump_dir)));
    f = ca.Function("f", new[] { x }, new[] { y });

    f.Call(42);

    string dump_path = dump_dir + "/myvar.000000.mtx";
    assertTrue(File.Exists(dump_path));
    val = DM.from_file(dump_path);
    checkarray(val, 42);

    removeGlob(dump_dir, "myvar.");
    Directory.Delete(dump_dir);

    if (RunSlow && CasadiMeta.feature_list().Contains("ghc-filesystem")) {
      foreach (string d in new[] { "mydump", "mydump_codegen" }) rmtree(d);

      x = MX.sym("x");
      y = ca.sqrt(x.dump("x", opts(("dir", "mydump"))));
      f = ca.Function("f", new[] { x }, new[] { y });
      checkCodegen(f, dms(ca.vertcat(9)), opts(("dump_dir_suffix", "_codegen")), 1e-15, "c99", true);

      assertEquals(readFile("mydump/x.000000.mtx"), readFile("mydump_codegen/x.000000.mtx"));

      foreach (string d in new[] { "mydump", "mydump_codegen" }) rmtree(d);
    }
  }

  public void test_codegen_specials() {
    MX x = MX.sym("x");
    MX y = MX.sym("y");

    foreach (MX z in new[] { ca.pow(x, 2), ca.if_else(y > 0, x * 2, x * y), ca.fmin(x, y), ca.fmax(x, y), ca.sign(x * y) }) {
      Function f = ca.Function("f", new[] { x, y }, new[] { z });
      checkCodegen(f, dms(ca.vertcat(1), ca.vertcat(2)));
      checkCodegen(f, dms(ca.vertcat(1), ca.vertcat(-2)));
    }
  }

  static MX getAny(MX m, object r, object c) {
    if (r is Slice rs && c is Slice cs) return m.get(false, rs, cs);
    if (r is Slice rs2) return m.get(false, rs2, (DM)c);
    if (c is Slice cs2) return m.get(false, (DM)r, cs2);
    return m.get(false, (DM)r, (DM)c);
  }

  static void setAny(MX m, object r, object c, MX v) {
    if (r is Slice rs && c is Slice cs) m.set(v, false, rs, cs);
    else if (r is Slice rs2) m.set(v, false, rs2, (DM)c);
    else if (c is Slice cs2) m.set(v, false, (DM)r, cs2);
    else m.set(v, false, (DM)r, (DM)c);
  }

  public void test_getsetnonzeros() {
    foreach (Sparsity S in new[] { Sparsity.lower(5), Sparsity.dense(5, 5) }) {
      MX M = MX.sym("X", S.nnz());

      MX m = ca.sin(ca.MX(S, M));
      foreach (int i in new[] { 0, 2 }) {
        object[][] inds = { new object[] { ca.Slice(i), ca.Slice(i) }, new object[] { ca.Slice(1), ca.Slice(i) },
          new object[] { ca.Slice(i), ca.Slice(1) }, new object[] { ca.Slice(), ca.Slice(i) },
          new object[] { ca.Slice(i), ca.Slice() }, new object[] { ca.Slice(i, i + 2), ca.Slice(i, i + 2) },
          new object[] { ca.Slice(i, i + 2), ca.Slice() }, new object[] { ca.vertcat(i, i), ca.vertcat(0, 0) },
          new object[] { (DM)new double[0], (DM)new double[0] } };
        foreach (object[] ind in inds) {
          MX E = getAny(m, ind[0], ind[1]);
          MX e = ca.cos(E);
          e = ca.dot(e, e);
          Function f = ca.Function("f", new[] { M }, new[] { e });
          checkfunction(f, f.expand(), dms(DM.rand(S.nnz(), 1)));
          checkCodegen(f, dms(DM.rand(S.nnz(), 1)), opts(), 1e-13, "c89", false);

          MX mc = m + 0;

          MX Y = MX.sym("y", E.nnz());
          MX y = ca.MX(E.sparsity(), Y);
          setAny(mc, ind[0], ind[1], y);
          e = ca.cos(mc);
          e = ca.dot(e, e);

          f = ca.Function("f", new[] { M, Y }, new[] { e });
          checkfunction(f, f.expand(), dms(DM.rand(S.nnz(), 1), DM.rand(E.nnz(), 1)));
          checkSerialize(f, dms(DM.rand(S.nnz(), 1), DM.rand(E.nnz(), 1)));
          checkCodegen(f, dms(DM.rand(S.nnz(), 1), DM.rand(E.nnz(), 1)), opts(), 1e-13, "c89", false);
        }
      }
    }
  }

  public void test_evalf() {
    MX x = MX.sym("x");

    MX p = MX.sym("p");
    Function f = ca.Function("f", new[] { x }, new[] { ca.sin(x) });
    MX y0 = f.call(new[] { p }, false, true)[0];
    MX y = ca.substitute(y0, p, 3);

    assertInException("not defined", () => y.to_DM());
    checkarray(ca.evalf(y), ca.sin(ca.DM(3)));
    assertInException("since variables [x] are free", () => ca.evalf(x));
  }

  void mminChecks(Func<Sparsity, Function> make) {
    Function f = make(Sparsity.dense(2, 2));
    List<DM> r = f.call(dms(ca.DM(new double[,] { { -2, 3 }, { -3, 5 } })));
    checkarray(r[0], -3);
    checkarray(r[1], 5);

    f = make(Sparsity.lower(2));
    r = f.call(dms(ca.sparsify(ca.DM(new double[,] { { -2, 0 }, { -3, -5 } }))));
    checkarray(r[0], -5);
    checkarray(r[1], 0);

    r = f.call(dms(ca.sparsify(ca.DM(new double[,] { { 2, 0 }, { 3, 5 } }))));
    checkarray(r[0], 0);
    checkarray(r[1], 5);

    f = make(ca.Sparsity(2, 2));
    r = f.call(dms(ca.DM(2, 2)));
    checkarray(r[0], 0);
    checkarray(r[1], 0);

    f = make(ca.Sparsity(2, 0));
    r = f.call(dms(ca.DM(2, 0)));
    assertTrue(r[0].is_empty());
    assertTrue(r[1].is_empty());

    f = make(ca.Sparsity(0, 0));
    r = f.call(dms(ca.DM(0, 0)));
    assertTrue(r[0].is_empty());
    assertTrue(r[1].is_empty());
  }

  public void test_mmin() {
    Func<Sparsity, Function> sx = sp => {
      SX x = SX.sym("x", sp);
      return ca.Function("f", new[] { x }, new[] { ca.mmin(x), ca.mmax(x) });
    };
    Func<Sparsity, Function> mx = sp => {
      MX x = MX.sym("x", sp);
      return ca.Function("f", new[] { x }, new[] { ca.mmin(x), ca.mmax(x) });
    };
    mminChecks(sx);
    mminChecks(sp => sx(sp).expand());
    mminChecks(sp => {
      Function f = sx(sp);
      SX y = SX.sym("y", sp);
      return ca.Function("f", new[] { y }, f.call(new[] { y }, false, false));
    });
    mminChecks(mx);
    mminChecks(sp => mx(sp).expand());
    mminChecks(sp => {
      Function f = mx(sp);
      MX y = MX.sym("y", sp);
      return ca.Function("f", new[] { y }, f.call(new[] { y }, true, false));
    });

    for (int k = 0; k < 2; k++) {
      MX x = MX.sym("X", 2);
      Function f = ca.Function("f", new[] { x }, new[] { k == 0 ? ca.mmin(x) : ca.mmax(x) });
      f.expand().disp(true);
      Function J = f.expand().jacobian();
      f.expand().jacobian().disp(true);
      Console.WriteLine("[" + string.Join(", ", J.call(dms(ca.vertcat(2, 2), ca.vertcat(0)))) + "]");
      checkfunction(f, f.expand(), dms(ca.vertcat(2, 2)));
      checkCodegen(f, dms(ca.vertcat(2, 2)), opts(), 1e-15, "c99", false);
    }
  }

  // test_doc_expression_tools: python-only (runtime __doc__ docstrings)

  public void test_densify() {
    Sparsity I = ca.sparsify(ca.DM(new double[,] { { 0, 1, 0, 1 }, { 1, 1, 0, 1 }, { 0, 1, 1, 0 } })).sparsity();

    MX x = MX.sym("x", I);

    Function f = ca.Function("f", new[] { x }, new[] { ca.densify(x) });

    DM a = ca.DM(new double[,] { { 0, 1, 0, 2 }, { 3, 4, 0, 5 }, { 0, 6, 7, 0 } });

    DM y = f.Call(ca.sparsify(a));
    assertTrue(y.sparsity().is_equal(a.sparsity()));
    checkarray(y, a);

    checkCodegen(f, dms(a));

    x = MX.sym("x", 3, 4);

    f = ca.Function("f", new[] { x }, new[] { x.get(I) });

    a = ca.DM(new double[,] { { 0, 1, 0, 2 }, { 3, 4, 0, 5 }, { 0, 6, 7, 0 } });
    DM b = ca.sparsify(a);
    a = ca.DM(new double[,] { { 9, 1, 9, 2 }, { 3, 4, 9, 5 }, { 9, 6, 7, 9 } });

    y = f.Call(a);
    assertTrue(y.sparsity().is_equal(b.sparsity()));
    checkarray(y, b);

    checkCodegen(f, dms(a));
    checkSerialize(f, dms(a));
  }

  public void test_constant_from_file() {
    File.WriteAllText("test.txt", "5.6 1e-5");
    assertInException("Failed to read a double from 'test.txt'. Expected 3 doubles.", () => ca.MX(Sparsity.lower(2), "test.txt"));

    File.WriteAllText("test.txt", "5.6 abc 1e-5");
    assertInException("Failed to read a double from 'test.txt'. Expected 3 doubles.", () => ca.MX(Sparsity.lower(2), "test.txt"));

    assertInException("Cannot open file 'testQWHAL567p.txt'.", () => ca.MX(Sparsity.lower(2), "testQWHAL567p.txt"));

    File.WriteAllText("test.txt", "5.6 1e-5 -12");
    MX x = ca.MX(Sparsity.lower(2), "test.txt");

    MX a = MX.sym("a");

    Function f = ca.Function("f", new[] { a }, new[] { a * x });
    checkarray(f.Call(2), ca.DM(Sparsity.lower(2), ca.vertcat(2 * 5.6, 2 * 1e-5, -2 * 12)));

    checkCodegen(f, dms(ca.vertcat(2)));

    assertInException("Not defined for ConstantFile", () => x.to_DM());
    File.Delete("test.txt");
  }

  public void test_nonzeros_param() {
    MX x = MX.sym("x", 2);
    MX y = MX.sym("y");

    foreach (int ad_weight_sp in new[] { 0, 1 }) {
      var o = opts(("helper_options", opts(("ad_weight_sp", ad_weight_sp))));
      MX z = x.get_nz(false, y);
      MX J = ca.jacobian(z, x, o);
      assertEquals(1L, J.size1());
      assertEquals(2L, J.size2());
      assertEquals(2L, J.nnz());
      J = ca.jacobian(z, y, o);
      assertEquals(1L, J.size1());
      assertEquals(1L, J.size2());
      assertEquals(0L, J.nnz());

      MX q = MX.sym("q");
      z = ca.MX(x);
      z.set_nz(q, false, y);

      J = ca.jacobian(z, x, o);
      assertEquals(2L, J.size1());
      assertEquals(2L, J.size2());
      assertEquals(4L, J.nnz());
      J = ca.jacobian(z, y, o);
      assertEquals(2L, J.size1());
      assertEquals(1L, J.size2());
      assertEquals(0L, J.nnz());
      J = ca.jacobian(z, q, o);
      assertEquals(2L, J.size1());
      assertEquals(1L, J.size2());
      assertEquals(2L, J.nnz());
    }
  }

  static List<string> lines(string s, int drop) {
    var r = s.Split('\n').ToList();
    return r.GetRange(0, r.Count - drop);
  }

  static List<string> printInstructionsCheck(Function f) {
    DM.rng(1);
    var inputs = new List<DM>();
    for (long i = 0; i < f.n_in(); i++) inputs.Add(DM.rand(f.sparsity_in(i)));
    List<string> @out = lines(captureStdout(() => f.call(inputs)), 1);

    assertTrue(@out[0].Contains("inputs"));

    if (RunSlow) {
      checkCodegenMainParsingError(f, inputs);
      List<string> out2 = lines(readFile(lastMainOut), 2);
      assertTrue(@out.SequenceEqual(out2));
      return out2;
    }
    return null;
  }

  public void test_print_instructions() {
    var popts = opts(("print_instructions", true), ("print_canonical", true));
    {
      MX x = MX.sym("x", Sparsity.diag(5)), y = MX.sym("y", 5, 5), z = MX.sym("z", 5), w = MX.sym("w");
      Function g = ca.Function("g", new[] { z }, new[] { z[0] + z[2], z[1] + z[4] }, opts(("never_inline", true)));
      List<MX> ab = g.call(new[] { ca.mtimes(y + x, z) });
      Function f = ca.Function("f", new[] { x, y, z, w }, new[] { ca.sin(ab[0] * (w * 7) + 6) * ab[1] }, popts);
      printInstructionsCheck(f);
    }
    {
      SX x = SX.sym("x", Sparsity.diag(5)), y = SX.sym("y", 5, 5), z = SX.sym("z", 5), w = SX.sym("w");
      Function g = ca.Function("g", new[] { z }, new[] { z[0] + z[2], z[1] + z[4] }, opts(("never_inline", true)));
      List<SX> ab = g.call(new[] { ca.mtimes(y + x, z) });
      Function f = ca.Function("f", new[] { x, y, z, w }, new[] { ca.sin(ab[0] * (w * 7) + 6) * ab[1] }, popts);
      printInstructionsCheck(f);
    }

    var results = new Dictionary<string, List<string>>();
    var results_slow = new Dictionary<string, List<string>>();
    foreach (string X in new[] { "MX", "SX" }) {
      Function f;
      if (X == "MX") {
        MX x = MX.sym("x"), y = MX.sym("y");
        f = ca.Function("f", new[] { x, y }, new[] { x + y }, popts);
      } else {
        SX x = SX.sym("x"), y = SX.sym("y");
        f = ca.Function("f", new[] { x, y }, new[] { x + y }, popts);
      }
      f.disp(true);
      DM.rng(1);
      var inputs = new List<DM>();
      for (long i = 0; i < f.n_in(); i++) inputs.Add(DM.rand(f.sparsity_in(i)));
      List<string> @out = lines(captureStdout(() => f.call(inputs)), 1);
      results[X] = @out;

      assertTrue(@out[0].Contains("inputs"));

      if (RunSlow) {
        checkCodegenMainParsingError(f, inputs);
        List<string> out2 = lines(readFile(lastMainOut), 2);
        results_slow[X] = out2;
        assertTrue(@out.SequenceEqual(out2));
      }
    }

    assertTrue(results["SX"].SequenceEqual(results["MX"]), string.Join("\n", results["SX"]) + "\nvs\n" + string.Join("\n", results["MX"]));
    if (RunSlow) assertTrue(results_slow["SX"].SequenceEqual(results_slow["MX"]));
  }

  public void test_low() {
    MX v = MX.sym("v", 5);
    double[] p0 = { -1, 0, 0.2, 1, 1.3, 2, 2.5, 3, 3.5, 4, 4.5 };
    MX p = MX.sym("p", p0.Length);

    Function f = ca.Function("f", new[] { v, p }, new[] { ca.low(v, p) });

    List<DM> inputs = dms(range(0, 5), ca.DM(p0));
    checkCodegen(f, inputs);
    checkSerialize(f, inputs);

    checkarray(f.call(inputs)[0], ca.vertcat(0, 0, 0, 1, 1, 2, 2, 3, 3, 3, 3));
  }

  public void test_linear_interpn() {
    int N = 4, nq = 3, n_out = 2;
    var g = new double[N];
    for (int i = 0; i < N; i++) g[i] = random.NextDouble();
    Array.Sort(g);
    var x0 = new double[nq];
    for (int i = 0; i < nq; i++) x0[i] = g[0] + (g[N - 1] - g[0]) * (i + 1) / (nq + 1);

    var d = new double[n_out * N];
    for (int i = 0; i < d.Length; i++) d[i] = random.NextDouble();

    Function F = ca.interpolant("F", "linear", new[] { g }, d);

    var rows = new List<DM>();
    foreach (double xi in x0) rows.Add(F.Call(xi).T());
    DM @ref = ca.vertcat(rows);

    MX x = MX.sym("x0", 1, N);
    MX v = MX.sym("v", N * n_out);
    MX xq = MX.sym("xq0", nq);

    MX e = MX.interpn_linear(new[] { x }, v, new[] { xq });

    Function f = ca.Function("f", new[] { x, v, xq }, new[] { e });
    List<DM> inputs = dms(ca.DM(g).T(), ca.DM(d), ca.DM(x0));
    DM r = f.call(inputs)[0];
    checkarray(r, @ref);
    checkCodegen(f, inputs, opts(), 1e-15, "c99", false);
    checkSerialize(f, inputs);
  }

  static DM sliceVals(Slice s, long len) => s.all(len).Select(v => (double)v).ToArray();

  static MX get2(MX m, object r, object c) {
    if (r is Slice rs && c is Slice cs) return m.get(false, rs, cs);
    if (r is Slice rs2) return m.get(false, rs2, (MX)c);
    if (c is Slice cs2) return m.get(false, (MX)r, cs2);
    return m.get(false, (MX)r, (MX)c);
  }

  public void test_paramslice() {
    int N = 4, M = 6;
    DM A = DM.rand(N, M);
    MX As = MX.sym("A", A.size1(), A.size2());
    MX ASs = MX.sym("AS", ca.tril(As.sparsity()));
    ASs.sparsity().spy();
    foreach (Slice r in new[] { ca.Slice(1), ca.Slice(1, 3), ca.Slice(0, N, 2), ca.Slice() }) {
      MX rs = MX.sym("r", A.get(r, ca.Slice()).size1());
      foreach (Slice c in new[] { ca.Slice(1), ca.Slice(1, 3), ca.Slice(0, M, 2), ca.Slice() }) {
        MX cs = MX.sym("c", A.get(ca.Slice(), c).size2());
        foreach (object rr in new object[] { r, rs }) {
          foreach (object cc in new object[] { c, cs }) {
            Function f = ca.Function("f", new[] { As, rs, cs }, new[] { get2(As, rr, cc) });
            checkarray(f.call(dms(A, sliceVals(r, N), sliceVals(c, M)))[0], A.get(r, c));
            if (rr is MX || cc is MX)
              assertInException("Parametric slicing only supported for dense matrices", () => get2(ASs, rr, cc));
          }
        }
      }
    }
    foreach (Slice k in new[] { ca.Slice(1), ca.Slice(1, 3), ca.Slice(0, N, 2), ca.Slice(0, 2 * N, 2), ca.Slice() }) {
      MX ks = MX.sym("r", A.get(k).size1());
      foreach (object kk in new object[] { k, ks }) {
        Function f = ca.Function("f", new[] { As, ks }, new[] { kk is Slice ks0 ? As.get(ks0) : As.get(false, (MX)kk) });
        checkarray(f.call(dms(A, sliceVals(k, N * M)))[0], A.get(k));
        assertInException("Parametric slicing only supported for dense matrices", () => ASs.get(false, ks));
      }
    }
  }

  public void test_mapping() {
    MX A = MX.sym("A", 4, 4);
    DM i = ca.DM(new double[,] { { 0, 3 }, { 1, 2 } });
    checkarray(i, A.get(i).mapping());
  }

  /// <summary>Symmetric eigendecomposition by cyclic Jacobi (reference for numpy.linalg.eigh): (D, V)</summary>
  static (double[] D, double[,] V) eigh(double[,] A0) {
    int n = A0.GetLength(0);
    var a = (double[,])A0.Clone();
    var v = new double[n, n];
    for (int i = 0; i < n; i++) v[i, i] = 1;
    for (int sweep = 0; sweep < 100; sweep++) {
      double off = 0;
      for (int i = 0; i < n; i++)
        for (int j = i + 1; j < n; j++) off += a[i, j] * a[i, j];
      if (off < 1e-30) break;
      for (int p = 0; p < n; p++)
        for (int q = p + 1; q < n; q++) {
          if (a[p, q] == 0) continue;
          double theta = (a[q, q] - a[p, p]) / (2 * a[p, q]);
          double t = Math.Sign(theta) / (Math.Abs(theta) + Math.Sqrt(theta * theta + 1));
          if (theta == 0) t = 1;
          double c = 1 / Math.Sqrt(t * t + 1), s = t * c;
          for (int k = 0; k < n; k++) {
            double akp = a[k, p], akq = a[k, q];
            a[k, p] = c * akp - s * akq;
            a[k, q] = s * akp + c * akq;
          }
          for (int k = 0; k < n; k++) {
            double apk = a[p, k], aqk = a[q, k];
            a[p, k] = c * apk - s * aqk;
            a[q, k] = s * apk + c * aqk;
          }
          for (int k = 0; k < n; k++) {
            double vkp = v[k, p], vkq = v[k, q];
            v[k, p] = c * vkp - s * vkq;
            v[k, q] = s * vkp + c * vkq;
          }
        }
    }
    var d = new double[n];
    for (int i = 0; i < n; i++) d[i] = a[i, i];
    return (d, v);
  }

  /// <summary>V diag(g(D)) V'</summary>
  static DM eigRecompose((double[] D, double[,] V) DV, Func<double, double> g) {
    DM V = DV.V;
    return ca.mtimes(ca.mtimes(V, ca.diag((DM)DV.D.Select(g).ToArray())), V.T());
  }

  public void test_convexify() {
    DM A = ca.diagcat(ca.DM(1), ca.DM(2), ca.DM(-1), ca.DM(new double[,] { { 1.2, 1.3 }, { 1.3, 4 } }),
                   ca.sparsify(ca.DM(new double[,] { { 0, 1, 0 }, { 1, 4, 7 }, { 0, 7, 9 } })), ca.DM(2, 2));

    double margin = 1e-7;

    var perm = Enumerable.Range(0, (int)A.size1()).Select(i => (double)i).ToArray();
    for (int i = perm.Length - 1; i > 0; i--) {
      int j = random.Next(i + 1);
      (perm[i], perm[j]) = (perm[j], perm[i]);
    }
    DM p = perm;

    var DV = eigh(A.ToArray());
    DM Ar_ref = eigRecompose(DV, d => Math.Max(Math.Abs(d), 1e-7));
    DM Ac_ref = eigRecompose(DV, d => Math.Max(d, 1e-7));
    MX As = MX.sym("As", A.sparsity());

    var optss = new List<Dictionary<string, object>> { opts(("strategy", "eigen-reflect")), opts(("strategy", "eigen-clip")),
                                                       opts(("strategy", "regularize")) };
    var refs = new List<DM> { Ar_ref, Ac_ref, A + DM.eye(A.size1()) * (4 + margin) };
    for (int k = 0; k < 3; k++) {
      var o = optss[k];
      DM @ref = refs[k];
      var ops = new List<Func<DM, DM>> { e => e, e => ca.triu(e), e => ca.tril(e), e => e.get(p, p), e => ca.triu(e.get(p, p)),
                                         e => ca.tril(e.get(p, p)) };
      if ((string)o["strategy"] == "regularize") ops = new List<Func<DM, DM>> { e => e, e => e.get(p, p) };
      foreach (var op in ops) {
        Function f = ca.Function("f", new[] { As }, new MX[] { ca.convexify(op(A), o) });

        checkarray(f.Call(A), op(@ref), 1e-8, "");

        checkSerialize(f, dms(A));

        checkCodegen(f, dms(A));
      }
    }
  }

  public void test_convexify_bugs() {
    foreach (double eps in new double[] { 0, 1e-100, 1 }) {
      double q = 4;
      double Q = 10;
      foreach (double s in new double[] { -4, 0.1, 0, 0.1, 4 }) {
        foreach (var trans in new List<Func<DM, DM>> { e => e, e => ca.sparsify(e) }) {
          DM A = trans(ca.DM(new double[,] { { q, s, 0 + eps }, { s, Q, 0 }, { 0 + eps, 0, 2 } }));
          Console.WriteLine(A);
          foreach (double d in eigh(A.ToArray()).D) assertTrue(d > 0);
          DM Ac = ca.evalf(ca.convexify(A, opts(("strategy", "eigen-reflect"))));
          checkarray(A, Ac, 1e-8, "");
        }
      }
    }
  }

  public void test_convexify_large() {
    foreach (int n in new[] { 3, 120 }) {
      DM B = randDM(n, n);
      DM A = B + B.T();
      var DV = eigh(A.ToArray());
      MX As = MX.sym("As", A.sparsity());
      foreach (string strategy in new[] { "eigen-reflect", "eigen-clip" }) {
        DM @ref = eigRecompose(DV, strategy == "eigen-reflect" ? d => Math.Max(Math.Abs(d), 1e-7) : d => Math.Max(d, 1e-7));
        Function f = ca.Function("f", new[] { As }, new[] { ca.convexify(As, opts(("strategy", strategy), ("max_iter_eig", 2000))) });
        checkarray(f.Call(A), @ref, 1e-8, "");
        if (n < 10) checkCodegen(f, dms(A));
      }
    }
  }

  static DM idx(params double[] i) => i;

  static DM revRange(long n) => Enumerable.Range(0, (int)n).Select(i => (double)(n - 1 - i)).ToArray();

  static List<MX> repmatSlices(MX a, long m, long n) {
    MX Y = ca.repmat(a, m, n);
    long R = Y.size1(), C = Y.size2();
    return new List<MX> { Y.get(ca.Slice(), revRange(C)), Y.get(revRange(R), ca.Slice()), Y.get(idx(R - 1, 0), idx(C - 1, 0)),
                          Y.get_nz(false, revRange(Y.nnz())), Y.T() };
  }

  static List<SX> repmatSlices(SX a, long m, long n) {
    SX Y = ca.repmat(a, m, n);
    long R = Y.size1(), C = Y.size2();
    return new List<SX> { Y.get(ca.Slice(), revRange(C)), Y.get(revRange(R), ca.Slice()), Y.get(idx(R - 1, 0), idx(C - 1, 0)),
                          Y.get_nz(false, revRange(Y.nnz())), Y.T() };
  }

  public void test_repmat_slice() {
    MX x = MX.sym("x", 2);
    MX y = ca.repmat(x, 1, 2);
    assertTrue((y[.., 0] - y[.., 1]).is_zero());
    assertTrue(ca.is_equal(y[.., 1], x));
    SX xs = SX.sym("x", 2);
    MX X = MX.sym("X", Sparsity.lower(2));
    SX Xs = SX.sym("X", Sparsity.lower(2));
    MX[] @as = { x, X };
    SX[] bs = { xs, Xs };
    for (int k = 0; k < 2; k++) {
      MX a = @as[k];
      SX b = bs[k];
      Function f = ca.Function("f", new[] { a }, new[] { ca.repmat(a, 1, 3).get(ca.Slice(), idx(2, 0, 1)),
        ca.repmat(a, 1, 3).get(idx(1, 0), idx(1, 2)), ca.repmat(a, 1, 3).get_nz(false, idx(5, 0, 2, 1)) });
      Function f_ref = ca.Function("f", new[] { b }, new[] { ca.repmat(b, 1, 3).get(ca.Slice(), idx(2, 0, 1)),
        ca.repmat(b, 1, 3).get(idx(1, 0), idx(1, 2)), ca.repmat(b, 1, 3).get_nz(false, idx(5, 0, 2, 1)) });
      checkfunction(f, f_ref, dms(ca.DM(a.sparsity(), DM.rand(a.nnz(), 1))));
    }
    var sps = new[] { Sparsity.dense(2, 3), Sparsity.triplet(3, 3, new long[] { 0, 2, 1 }, new long[] { 0, 0, 2 }),
                      Sparsity.triplet(3, 2, new long[] { 0, 2, 2 }, new long[] { 0, 0, 1 }), ca.Sparsity(2, 3) };
    foreach (Sparsity sp in sps) {
      foreach (long[] mn in new[] { new long[] { 1, 3 }, new long[] { 2, 2 } }) {
        long m = mn[0], n = mn[1];
        MX a = MX.sym("a", sp);
        SX b = SX.sym("a", sp);
        MX Y = ca.repmat(a, m, n);
        long w = Y.size2() / n;
        assertTrue((Y[.., ca.Slice(w * (n - 1), Y.size2())] - Y[.., ca.Slice(0, w)]).is_zero());
        Function f = ca.Function("f", new[] { a }, repmatSlices(a, m, n));
        Function f_ref = ca.Function("f", new[] { b }, repmatSlices(b, m, n));
        checkfunction(f, f_ref, dms(ca.DM(sp, DM.rand(sp.nnz(), 1))));
      }
    }
  }

  public void test_logsumexp() {
    MX x = MX.sym("x", 3);

    Function f_ref = ca.Function("f_ref", new[] { x }, new[] { ca.log(ca.exp(x[0]) + ca.exp(x[1]) + ca.exp(x[2])) });
    Function f = ca.Function("f", new[] { x }, new[] { ca.logsumexp(x) });

    checkfunction(f, f_ref, dms(ca.vertcat(1.1, 1.3, 1.7)));
    checkCodegen(f, dms(ca.vertcat(1.1, 1.3, 1.7)));
    checkfunction(f, f_ref, dms(ca.vertcat(1.1, 1.3, 1.3)));
    checkfunction(f, f_ref, dms(ca.vertcat(1.3, 1.3, 1.3)));

    DM res = f.Call(ca.vertcat(100, 1000, 10000));
    checkarray(res, 10000);
  }

  public void test_empty_broadcast() {
    foreach (long nc in new long[] { 0, 2 }) {
      MX res = ca.atan2(MX.sym("c", nc, 1), MX.sym("t", nc, 3));

      assertEquals(nc, res.size1());
      assertEquals(3L, res.size2());

      assertInException("Dimension mismatch", () => ca.atan2(MX.sym("c", nc, 2), MX.sym("t", nc, 3)));

      res = ca.atan2(MX.sym("c", nc, 2), MX.sym("t", nc, 4));
      assertEquals(nc, res.size1());
      assertEquals(4L, res.size2());

      res = ca.atan2(MX.sym("c", nc, 4), MX.sym("t", nc, 2));
      assertEquals(nc, res.size1());
      assertEquals(4L, res.size2());
    }
  }

  public void test_horzcat_sparsity() {
    MX x = MX.sym("x", 2, 2);
    Function f = null, f_alt = null;

    foreach (MX y in new[] { MX.sym("y", ca.Sparsity(2, 2)), MX.sym("y", Sparsity.lower(2)) }) {
      MX z = ca.sin(ca.horzcat(x, ca.sin(y), x));
      MX z_alt = ca.sin(ca.sparsity_cast(ca.vertcat(ca.vec(x), ca.vec(ca.sin(y)), ca.vec(x)), z.sparsity()));
      Console.WriteLine(z_alt);

      f = ca.Function("f", new[] { x, y }, new[] { z });
      f_alt = ca.Function("f_alt", new[] { x, y }, new[] { z });

      foreach (Function F in new[] { f, f_alt }) {
        assertEquals(z.nnz(), F.call(new[] { x, y }, true, false)[0].nnz());
        assertEquals(z.nnz(), F.call(new[] { x, y }, false, true)[0].nnz());
        assertEquals(z.nnz(), F.call(new[] { x, ca.pow(x, 2) }, true, false)[0].nnz());
        assertEquals(z.nnz(), F.call(new[] { x, ca.pow(x, 2) }, false, true)[0].nnz());
      }
    }

    List<DM> @in = dms(ca.DM(new double[,] { { 1, 2 }, { 3, 4 } }), ca.DM(new double[,] { { 5, 6 }, { 7, 8 } }));
    checkfunction(f, f.expand(), @in);
    checkfunction(f, f_alt, @in);
    checkCodegen(f, @in);
  }

  public void test_sparsity_cast_ad() {
    MX y = MX.sym("y", 3);

    MX xy = ca.vertcat(y[0], y[1]);
    MX w = ca.mtimes(ca.MX(ca.sparsify(ca.vertcat(1, 1, 0)).sparsity().T(), xy), y);

    Function f = ca.Function("f", new[] { y }, new[] { w });
    checkfunction(f, f.expand(), dms(ca.vertcat(1.1, 1.3, 1.7)));
    checkCodegen(f, dms(ca.vertcat(1.1, 1.3, 1.7)));
  }

  public void test_fractional_slicing() {
    MX t = MX.sym("x");

    Function f = ca.Function("f", new[] { t }, new[] { ca.MX(ca.vertcat(1, 2)).get(false, t) });

    checkarray(f.Call(0), 1);
    checkarray(f.Call(1), 2);
    checkarray(f.Call(0.5), 1);
    checkarray(f.Call(0.9), 1);

    f = ca.Function("f", new[] { t }, new[] { ca.MX(ca.DM(new double[,] { { 1, 2 }, { 3, 4 } })).get(false, ca.Slice(), t) });

    checkarray(f.Call(0), ca.vertcat(1, 3));
    checkarray(f.Call(1), ca.vertcat(2, 4));
    checkarray(f.Call(0.5), ca.vertcat(1, 3));
    checkarray(f.Call(0.9), ca.vertcat(1, 3));
  }

  static void checkfunctionLight(Function f, Function g, List<DM> inputs) {
    List<DM> a = f.call(inputs), b = g.call(inputs);
    assertEquals(a.Count, b.Count);
    for (int i = 0; i < a.Count; i++) {
      if (f.sparsity_out(i).is_empty() || g.sparsity_out(i).is_empty()) continue;
      checkarray(a[i], b[i], 1e-9, "output " + i);
    }
  }

  public void test_primitives() {
    MX x = MX.sym("x", 3);
    MX y = MX.sym("y", 3);
    var experiments = new List<MX> { ca.vertcat(x, y), ca.horzcat(x, y), ca.diagcat(x, y), ca.vertcat(x, y, ca.MX(0, 1)),
                                     ca.reshape(ca.vertcat(x, y), 2, 3), ca.sparsity_cast(ca.vertcat(x, y), Sparsity.diag(6)) };

    foreach (MX z in experiments) {
      Console.WriteLine(z);
      assertTrue(z.is_valid_input());

      Console.WriteLine("[" + string.Join(", ", z.primitives()) + "]");
      Function f = ca.Function("f", new[] { z }, z.primitives());
      DM.rng(1);
      DM a = DM.rand(z.sparsity());

      List<DM> @ref = f.call(dms(a));
      List<DM> res = z.split_primitives(a);

      assertEquals(@ref.Count, res.Count);
      for (int i = 0; i < @ref.Count; i++) checkarray(@ref[i], res[i]);

      DM a_back = z.join_primitives(res);
      checkarray(a, a_back);

      SX zsx = SX.sym("x", z.sparsity());
      Function f2 = ca.Function("f", new[] { zsx }, z.split_primitives(zsx));

      checkfunctionLight(f, f2, dms(a));
      checkCodegen(f, dms(a));
      checkCodegen(f2, dms(a));

      Function f3 = ca.Function("f", new[] { zsx }, new[] { z.join_primitives(z.split_primitives(zsx)) });

      DM res3 = f3.Call(a);
      checkarray(a, res3);
    }
  }

  /// <summary>Python hash(MX)</summary>
  static long hashOf(MX x) => x.hash();

  public void test_cache_output_node() {
    MX x = MX.sym("x");
    MX y = MX.sym("y");
    Function f = ca.Function("f", new[] { x, y }, new[] { x + y, x / y });

    List<MX> ab = f.call(new[] { ca.sin(x), ca.cos(y) });
    MX a = ab[0], b = ab[1];
    MX @base = a.dep(0);
    assertTrue(hashOf(a) != 0, hashOf(@base.get_output(0)).ToString());
    assertTrue(hashOf(b) != 0, hashOf(@base.get_output(1)).ToString());
  }

  public void test_naked_call_function() {
    MX x = MX.sym("x");
    MX y = MX.sym("y");
    Function f = ca.Function("f", new[] { x, y }, new[] { x + y, x / y });

    List<MX> ab = f.call(new[] { ca.sin(x), ca.cos(y) });

    MX @base = ab[0].dep(0);

    assertInException("get_output", () => ca.Function("g", new[] { x, y }, new[] { @base }));

    assertInException("get_output", () => { var _ = @base * 3; });
  }

  public void test_substitute() {
    MX x = MX.sym("x");
    MX y = MX.sym("y");
    MX w = MX.sym("w");

    MX xy = x * y;

    MX z = xy + w;
    MX z2 = ca.substitute(new[] { z }, new[] { w }, new[] { ca.MX(3) })[0];
    Console.WriteLine(xy);
    Console.WriteLine(z2.dep(1));
    assertTrue(hashOf(xy) != hashOf(z2.dep(1)));
  }

  // ---- MX expressions and their SX twins, handled as object (SX or MX) ----

  static Function mkF(string name, IEnumerable<object> ins, IEnumerable<object> outs) {
    var i = ins.ToList();
    if (i[0] is SX) return ca.Function(name, i.Cast<SX>(), outs.Cast<SX>());
    return ca.Function(name, i.Cast<MX>(), outs.Cast<MX>());
  }

  /// <summary>The SX twin of MX expressions over syms: evaluate them symbolically on SX symbols</summary>
  static List<object> conv(bool isMX, List<MX> syms, List<object> tsyms, IEnumerable<MX> exprs) {
    if (isMX) return exprs.Cast<object>().ToList();
    return ca.Function("conv", syms, exprs).call(tsyms.Cast<SX>()).Cast<object>().ToList();
  }

  /// <summary>ms itself for MX, else SX symbols of the same names and sparsities</summary>
  static List<object> twins(bool isMX, List<MX> ms) => ms.Select(m => isMX ? (object)m : SX.sym(m.name(), m.sparsity())).ToList();

  static object add(object a, object b) => a is MX ma ? ma + (MX)b : (object)((SX)a + (SX)b);

  static (object, object, object) separateLinear(object e, object x, object p) {
    if (e is MX m) {
      var (a, b, c) = ca.separate_linear(m, (MX)x, (MX)p);
      return (a, b, c);
    } else {
      var (a, b, c) = ca.separate_linear((SX)e, (SX)x, (SX)p);
      return (a, b, c);
    }
  }

  static bool isLinear(object e, object x) => e is MX m ? ca.is_linear(m, (MX)x) : ca.is_linear((SX)e, (SX)x);

  static bool dependsOn(object e, object x) => e is MX m ? ca.depends_on(m, (MX)x) : ca.depends_on((SX)e, (SX)x);

  static (object, List<object>, List<object>) extractParametric(object e, object p, Dictionary<string, object> o) {
    if (e is MX m) {
      var (a, b, c) = ca.extract_parametric(m, (MX)p, o);
      return (a, b.Cast<object>().ToList(), c.Cast<object>().ToList());
    } else {
      var (a, b, c) = ca.extract_parametric((SX)e, (SX)p, o);
      return (a, b.Cast<object>().ToList(), c.Cast<object>().ToList());
    }
  }

  static List<object> substituteAny(List<object> ex, List<object> v, List<object> vdef) {
    if (ex[0] is MX) return ca.substitute(ex.Cast<MX>(), v.Cast<MX>(), vdef.Cast<MX>()).Cast<object>().ToList();
    return ca.substitute(ex.Cast<SX>(), v.Cast<SX>(), vdef.Cast<SX>()).Cast<object>().ToList();
  }

  static string nameOf(object e) => e is MX m ? m.name() : ((SX)e).name();

  static string str(object o) => o is System.Collections.IEnumerable l && !(o is string)
    ? "[" + string.Join(", ", l.Cast<object>()) + "]" : o.ToString();

  /// <summary>Same values at random inputs</summary>
  static bool numEqual(List<object> syms, List<object> a, List<object> b) {
    Function f1 = mkF("f", syms, a);
    Function f2 = mkF("f", syms, b);
    DM.rng(1);
    var args = new List<DM>();
    for (long i = 0; i < f1.n_in(); i++) args.Add(DM.rand(f1.sparsity_in(i)));
    List<DM> r1 = f1.call(args, false, false), r2 = f2.call(args, false, false);
    for (int i = 0; i < Math.Min(r1.Count, r2.Count); i++)
      if ((double)ca.norm_1(ca.vec(ca.densify(r1[i] - r2[i]))) > 1e-12) return false;
    return true;
  }

  static List<object> O(params object[] x) => x.ToList();

  public void test_extract_parametric() {
    MX x = MX.sym("x"), y = MX.sym("y"), p = MX.sym("p");
    List<MX> ms = L(x, y, p);
    foreach (bool isMX in new[] { false, true }) {
      List<object> ts = twins(isMX, ms);
      object tp = ts[2];

      List<object> e = conv(isMX, ms, ts, L(x * 2 * p + y * 2 * ca.pow(p, 2), ca.pow(p, 2),
        ca.vertcat(x * 2 + y * 3 + p * 8, x * 3 + y * 2 + 4, ca.MX(7)), p * 8));
      var r = extractParametric(e[0], tp, opts(("extract_trivial", true)));
      Console.WriteLine(r.Item1 + " " + str(r.Item3));
      assertTrue(numEqual(ts, r.Item3, O(tp, e[1])));

      r = extractParametric(e[2], tp, opts(("extract_trivial", true)));

      Console.WriteLine(r.Item1 + " " + str(r.Item3));

      assertTrue(numEqual(ts, O(r.Item3[0]), O(e[3])));
    }

    foreach (bool isMX in new[] { false, true }) {
      List<object> ts = twins(isMX, ms);
      object tp = ts[2];

      MX cp = ca.cos(p);
      MX expr = (x - ca.sqrt(p)) * (y + ca.pow(cp, 2)) + cp / (x + cp) / (y - p);

      foreach (object ex in conv(isMX, ms, ts, L(expr, ca.vertcat(expr, ca.pow(p, 2))))) {
        var r = extractParametric(ex, tp, opts(("extract_trivial", true)));

        Console.WriteLine(r.Item1 + " " + str(r.Item2) + " " + str(r.Item3));

        assertFalse(dependsOn(r.Item1, tp));

        object expr_recreated = substituteAny(O(r.Item1), r.Item2, r.Item3)[0];
        assertTrue(numEqual(ts, O(ex), O(expr_recreated)));
      }
    }
  }

  static void separateLinearCases(bool isMX, List<MX> ms, List<object> ts, object xy, List<MX[]> cases) {
    object tp = ts[2];
    foreach (MX[] c in cases) {
      List<object> e = conv(isMX, ms, ts, c);
      Console.WriteLine("ref " + str(e));

      var (expr_const, expr_lin, expr_nonlin) = separateLinear(e[0], xy, tp);
      Console.WriteLine("actual (" + expr_const + ", " + expr_lin + ", " + expr_nonlin + ")");

      assertTrue(numEqual(ts, O(e[0]), O(add(add(expr_const, expr_lin), expr_nonlin))));

      assertTrue(numEqual(ts, O(expr_const), O(e[1])));
      assertTrue(numEqual(ts, O(expr_lin), O(e[2])));
      assertTrue(numEqual(ts, O(expr_nonlin), O(e[3])));

      assertTrue(isLinear(expr_lin, xy));
      assertFalse(dependsOn(expr_const, xy));
    }
  }

  public void test_separate_linear() {
    foreach (bool isMX in new[] { false, true }) {
      MX x = MX.sym("x"), y = MX.sym("y"), p = MX.sym("p");
      MX O = ca.MX(0);
      MX cp = ca.cos(p);
      MX xy3 = x * 3 + y;

      var cases = new List<MX[]> {
        new[] { x, O, x, O },
        new[] { x + p, p, x, O },
        new[] { x + y, O, x + y, O },
        new[] { x - p, -p, x, O },
        new[] { p - x, p, -x, O },
        new[] { x * 2, O, x * 2, O },
        new[] { -x, O, -x, O },
        new[] { xy3, O, xy3, O },
        new[] { cp, cp, O, O },
        new[] { xy3 + cp, cp, xy3, O },
        new[] { x * 3 + 7, ca.MX(7), x * 3, O },
        new[] { p * (x * 3 + 7), p * 7, p * (x * 3), O },
        new[] { y * (x * 3 + 7), O, y * 7, y * (x * 3) },
        new[] { ca.cos(7 + x * 3), O, O, ca.cos(7 + x * 3) },
        new[] { ca.cos(x * 3 + 7), O, O, ca.cos(x * 3 + 7) },
        new[] { ca.cos(p * 7), ca.cos(p * 7), O, O },
        new[] { ca.cos(p + x), O, O, ca.cos(p + x) },
        new[] { x * cp + 7, ca.MX(7), x * cp, O },
        new[] { y * (x * cp + 7), O, y * 7, y * x * cp },
        new[] { p * (x + y + x * y), O, p * (x + y), p * (x * y) },
        new[] { ca.cos(x + y + x * y), O, O, ca.cos(x + y + x * y) },
        new[] { p * (p + cp * x), ca.pow(p, 2), p * cp * x, O },
        new[] { (p + x + x * y) * (cp + x * 7 + x * 3 * y), p * cp, p * (x * 7) + x * cp,
                p * (x * 3 * y) + x * y * (cp + x * 7 + x * 3 * y) + ca.pow(x, 2) * 7 + ca.pow(x, 2) * 3 * y },
        new[] { (p + x + x * y) / cp, p / cp, x / cp, x * y / cp },
        new[] { (p + x + x * y) / (cp + x * 7), O, O, (p + x + x * y) / (cp + x * 7) },
        new[] { (p + x + x * y) / (cp + x * 7 + x * 3 * y), O, O, (p + x + x * y) / (cp + x * 7 + x * 3 * y) },
      };
      List<MX> ms = L(x, y, p);
      List<object> ts = twins(isMX, ms);
      separateLinearCases(isMX, ms, ts, conv(isMX, ms, ts, L(ca.vertcat(x, y)))[0], cases);

      x = MX.sym("x", 2);
      y = MX.sym("y", 2);
      p = MX.sym("p", 2);
      cp = ca.cos(p);
      MX xy = x * y;
      cases = new List<MX[]> {
        new[] { ca.dot(p, x), O, ca.dot(p, x), O },
        new[] { ca.dot(p + x + xy, cp + x * 7 + xy * 3), ca.dot(p, cp), ca.dot(p, x * 7) + ca.dot(x, cp),
                ca.dot(p, xy * 3) + ca.dot(x, x * 7 + xy * 3) + ca.dot(xy, cp + x * 7 + xy * 3) },
      };
      ms = L(x, y, p);
      ts = twins(isMX, ms);
      separateLinearCases(isMX, ms, ts, conv(isMX, ms, ts, L(ca.vertcat(x, y)))[0], cases);

      x = MX.sym("x", 2, 2);
      y = MX.sym("y", 2, 2);
      p = MX.sym("p", 2, 2);
      cp = ca.cos(p);
      xy = x * y;
      cases = new List<MX[]> {
        new[] { (p + x + xy).T(), p.T(), x.T(), xy.T() },
        new[] { ca.mtimes(p, x), O, ca.mtimes(p, x), O },
        new[] { ca.mtimes(p + x + xy, cp + x * 7 + xy * 3), ca.mtimes(p, cp), ca.mtimes(p, x * 7) + ca.mtimes(x, cp),
                ca.mtimes(p, xy * 3) + ca.mtimes(x, x * 7 + xy * 3) + ca.mtimes(xy, cp + x * 7 + xy * 3) },
      };
      ms = L(x, y, p);
      ts = twins(isMX, ms);
      separateLinearCases(isMX, ms, ts, conv(isMX, ms, ts, L(ca.veccat(x, y)))[0], cases);
    }
  }

  public void test_separate_linear_mat() {
    MX x = MX.sym("x", 3);
    MX y = MX.sym("y", 3, 3);
    MX z = MX.sym("z", 2);
    MX p = MX.sym("p");

    DM A = ca.DM(new double[,] { { 0, 1, 3 }, { 2, 0, 0 }, { 0, 9, 0 } });
    DM Asp = ca.sparsify(A);

    DM B = ca.DM(new double[,] { { 1.7, 1.1 }, { 6, 0 }, { 0, 1.2 } });
    DM Bsp = ca.sparsify(B);
    long[] off = { 0, 2, 7 };
    MX yx = ca.mtimes(y, x);
    MX Z3 = ca.MX(DM.zeros(3, 1));
    MX full = ca.vertcat(ca.MX(3), x * 3, yx) * 9;
    MX fullc = ca.vertcat(ca.MX(3), Z3, Z3) * 9;
    MX fulll = ca.vertcat(ca.MX(0), x * 3, Z3) * 9;
    MX fulln = ca.vertcat(ca.MX(0), Z3, yx) * 9;
    MX p2 = ca.vertsplit(full, off)[1];
    MX p2c = ca.vertsplit(fullc, off)[1];
    MX p2l = ca.vertsplit(fulll, off)[1];
    MX p2n = ca.vertsplit(fulln, off)[1];
    DM nzi = idx(2, 3, 0, 4, 6, 6);
    MX E = full.get_nz(false, nzi);
    MX Ec = fullc.get_nz(false, nzi);
    MX El = fulll.get_nz(false, nzi);
    MX En = fulln.get_nz(false, nzi);

    MX Z5 = ca.MX(DM.zeros(5, 1)), Z3b = ca.MX(DM.zeros(3, 1)), Z2 = ca.MX(DM.zeros(2, 1));
    MX[][] cases = {
      new[] { ca.mtimes(A, x) + ca.mtimes(B, z), Z3b, ca.mtimes(A, x), ca.mtimes(B, z) },
      new[] { ca.mtimes(Asp, x) + ca.mtimes(Bsp, z), Z3b, ca.mtimes(Asp, x), ca.mtimes(Bsp, z) },
      new[] { ca.vertcat(x * 3 + yx, z * 9), Z5, ca.vertcat(x * 3, Z2), ca.vertcat(yx, z * 9) },
      new[] { ca.vertcat(x * 3 + yx, z * 9) * 7, Z5, ca.vertcat(x * 3, Z2) * 7, ca.vertcat(yx, z * 9) * 7 },
      new[] { p2 * 9, p2c * 9, p2l * 9, p2n * 9 },
      new[] { E, Ec, El, En } };

    foreach (MX[] c in cases) {
      Console.WriteLine("separate_linear(" + c[0] + ")");
      Console.WriteLine("ref");
      Console.WriteLine("  const " + c[1]);
      Console.WriteLine("  lin " + c[2]);
      Console.WriteLine("  nonlin " + c[3]);

      var (expr_const, expr_lin, expr_nonlin) = ca.separate_linear(c[0], L(x, y), L(p));
      Console.WriteLine("actual");
      Console.WriteLine("  const " + expr_const);
      Console.WriteLine("  lin " + expr_lin);
      Console.WriteLine("  nonlin " + expr_nonlin);

      Action<MX, MX> test_equal = (a, b) => {
        Function f1 = ca.Function("f", new[] { x, y, z, p }, new[] { a });
        Function f2 = ca.Function("f", new[] { x, y, z, p }, new[] { b });
        DM.rng(1);
        var args = new List<DM>();
        for (long i = 0; i < f1.n_in(); i++) args.Add(DM.rand(f1.sparsity_in(i)));
        checkarray(f1.call(args)[0], f2.call(args)[0]);
      };

      test_equal(c[0], expr_const + expr_lin + expr_nonlin);

      test_equal(expr_const, c[1]);
      test_equal(expr_lin, c[2]);
      test_equal(expr_nonlin, c[3]);

      assertTrue(ca.is_linear(expr_lin, ca.veccat(x, y)));
      assertFalse(ca.depends_on(expr_const, ca.veccat(x, y)));
    }
  }

  public void test_extract_parametric_opts() {
    MX x = MX.sym("x"), y = MX.sym("y"), p = MX.sym("p");
    List<MX> ms = L(x, y, p);
    foreach (bool isMX in new[] { false, true }) {
      List<object> ts = twins(isMX, ms);
      object tp = ts[2];
      object expr = conv(isMX, ms, ts, L(x * 2 * p + y * 2 * ca.pow(p, 2)))[0];

      var r = extractParametric(expr, tp, opts(("extract_trivial", true)));
      List<object> symbols = r.Item2;

      assertEquals(2, symbols.Count);

      assertEquals("e_0", nameOf(symbols[0]));
      assertEquals("e_1", nameOf(symbols[1]));

      r = extractParametric(expr, tp, opts(("offset", 5), ("prefix", "foo_"), ("suffix", "bar"), ("extract_trivial", true)));
      symbols = r.Item2;

      assertEquals("foo_5bar", nameOf(symbols[0]));
      assertEquals("foo_6bar", nameOf(symbols[1]));

      Console.WriteLine(str(symbols));
      Console.WriteLine(str(r.Item3));

      r = extractParametric(expr, tp, opts(("extract_trivial", false)));
      assertEquals(1, r.Item2.Count);
      assertTrue(str(r.Item3).Contains("sq(p)"));
      Console.WriteLine(str(r.Item3));

      r = extractParametric(expr, tp, opts(("extract_trivial", true)));
      assertEquals(2, r.Item2.Count);
      assertTrue(str(r.Item3).Contains("sq(p)"));
      Console.WriteLine(str(r.Item3));
    }
  }

  public void test_weakref() {
    MX x = MX.sym("x");

    WeakRef y = ca.WeakRef(x);

    assertTrue(y.alive());

    x.Dispose();

    assertFalse(y.alive());

    x = MX.sym("x");
    y = ca.WeakRef(x);
    SharedObject xx = y.shared();
    xx.Dispose();

    assertTrue(y.alive());
  }

  public void test_nonzeros() {
    MX A = MX.sym("A", 3, 3);

    MX x = MX.sym("x");
    MX y = MX.sym("y");
    MX z = MX.sym("z");
    MX T = MX.sym("T", Sparsity.upper(3));
    MX w = MX.sym("w", 3, 3);

    DM.rng(1);

    foreach (MX expr in new[] { ca.vertcat(x, y, z), ca.veccat(x, y, z, w), ca.veccat(x, T, z, w), x * 2, T * 4 * x, ca.vertcat(x, ca.vertcat(y, z)) }) {
      Function fref = ca.Function("fref", new[] { x, y, z, w, A, T }, new[] { expr.get_nz(false, ca.Slice()) });
      Function f = ca.Function("f", new[] { x, y, z, w, A, T }, new[] { ca.vertcat(expr.nonzeros()) });

      var @in = new List<DM>();
      for (long i = 0; i < f.n_in(); i++) @in.Add(DM.rand(f.sparsity_in(i)));
      checkfunctionLight(f, fref, @in);
    }

    foreach (MX expr in new[] { ca.vertcat(x, y, z), ca.vertcat(x, ca.vertcat(y, z)) }) {
      string s = str(expr.nonzeros());
      // Python checks the list repr "[MX(x), MX(y), MX(z)]"
      assertTrue(s.Contains("[x, y, z]"), s);
    }
  }

  public void test_printme_codegen() {
    {
      SX x = SX.sym("x");
      Function f = ca.Function("f", new[] { x }, new[] { ca.pow(x, 2).printme(2) + 6 });
      checkCodegen(f, dms(ca.vertcat(3)));
    }
    {
      MX x = MX.sym("x");
      Function f = ca.Function("f", new[] { x }, new[] { ca.pow(x, 2).printme(2) + 6 });
      checkCodegen(f, dms(ca.vertcat(3)));
    }
  }

  void jtimesEmptyMX() {
    foreach (long nx in new long[] { 2, 1, 0 }) {
      foreach (long mx in new long[] { 3, 1, 0 }) {
        MX x = MX.sym("x", nx, mx);
        foreach (long ny in new long[] { 3, 1, 0 }) {
          foreach (long my in new long[] { 4, 1, 0 }) {
            foreach (long bm in new long[] { 1, 3, 0 }) {
              if (bm != 1 && my == 0) continue;
              MX y = ca.mtimes(ca.mtimes(DM.rand(ny, nx), x), DM.rand(mx, my));
              assertEquals(ny, y.size1());
              assertEquals(my, y.size2());
              DM yb = DM.rand(ny, my * bm);

              MX r = ca.jtimes(y, x, ca.MX(yb), true);
              var parts = new List<MX>();
              foreach (DM e in ca.horzsplit(yb, Math.Max(1, my)))
                parts.Add(ca.reshape(ca.mtimes(ca.jacobian(y, ca.vec(x)).T(), ca.MX(ca.vec(e))), x.size1(), x.size2()));
              // Python hcat([]) is DM(1,0)
              MX r_ref = parts.Count == 0 ? ca.MX(1, 0) : ca.horzcat(parts);
              assertEquals(nx, r.size1());
              assertEquals(mx * bm, r.size2());
              Function f = ca.Function("f", new[] { x }, new[] { r });
              Function f_ref = ca.Function("f", new[] { x }, new[] { r_ref });
              checkfunctionLight(f, f_ref, dms(DM.rand(nx, mx)));
            }
          }
        }
      }
    }
    foreach (long nx in new long[] { 2, 1, 0 }) {
      foreach (long mx in new long[] { 3, 1, 0 }) {
        MX x = MX.sym("x", nx, mx);
        foreach (long ny in new long[] { 3, 1, 0 }) {
          foreach (long my in new long[] { 4, 1, 0 }) {
            foreach (long bm in new long[] { 1, 3, 0 }) {
              if (bm != 1 && mx == 0) continue;
              MX y = ca.mtimes(ca.mtimes(DM.rand(ny, nx), x), DM.rand(mx, my));
              assertEquals(ny, y.size1());
              assertEquals(my, y.size2());
              DM dx = DM.rand(nx, mx * bm);

              MX r = ca.jtimes(y, x, ca.MX(dx));
              var parts = new List<MX>();
              foreach (DM e in ca.horzsplit(dx, Math.Max(1, mx)))
                parts.Add(ca.reshape(ca.mtimes(ca.jacobian(y, ca.vec(x)), ca.MX(ca.vec(e))), y.size1(), y.size2()));
              MX r_ref = parts.Count == 0 ? ca.MX(1, 0) : ca.horzcat(parts);

              assertEquals(ny, r.size1());
              assertEquals(my * bm, r.size2());
              Function f = ca.Function("f", new[] { x }, new[] { r });
              Function f_ref = ca.Function("f", new[] { x }, new[] { r_ref });
              checkfunctionLight(f, f_ref, dms(DM.rand(nx, mx)));
            }
          }
        }
      }
    }
  }

  void powZeroMX() {
    MX x = MX.sym("x", Sparsity.diag(4));
    assertTrue(ca.pow(x, 0).is_dense());
    assertFalse(ca.pow(x, 0.3).is_dense());
    assertFalse(ca.pow(x, 1).is_dense());
    assertFalse(ca.pow(x, 2).is_dense());
    assertFalse(ca.pow(x, 2.3).is_dense());
    assertTrue(ca.pow(x, -2).is_dense());
    MX p = MX.sym("p");

    assertTrue(ca.pow(x, p).is_dense());

    MX s = MX.sym("x");

    checkarray(ca.evalf(ca.substitute(ca.vertcat(ca.pow(s, 2), ca.pow(s, 1), ca.pow(s, 0)), s, ca.MX(0))), ca.vertcat(0, 0, 1));
    checkarray(ca.evalf(ca.substitute(ca.pow(s, ca.vertcat(2, 1, 0)), s, ca.MX(0))), ca.vertcat(0, 0, 1));
  }

  void jtimesEmptySX() {
    foreach (long nx in new long[] { 2, 1, 0 }) {
      foreach (long mx in new long[] { 3, 1, 0 }) {
        SX x = SX.sym("x", nx, mx);
        foreach (long ny in new long[] { 3, 1, 0 }) {
          foreach (long my in new long[] { 4, 1, 0 }) {
            foreach (long bm in new long[] { 1, 3, 0 }) {
              if (bm != 1 && my == 0) continue;
              SX y = ca.mtimes(ca.mtimes(DM.rand(ny, nx), x), DM.rand(mx, my));
              assertEquals(ny, y.size1());
              assertEquals(my, y.size2());
              DM yb = DM.rand(ny, my * bm);

              SX r = ca.jtimes(y, x, ca.SX(yb), true);
              var parts = new List<SX>();
              foreach (DM e in ca.horzsplit(yb, Math.Max(1, my)))
                parts.Add(ca.reshape(ca.mtimes(ca.jacobian(y, ca.vec(x)).T(), ca.SX(ca.vec(e))), x.size1(), x.size2()));
              SX r_ref = parts.Count == 0 ? ca.SX(1, 0) : ca.horzcat(parts);
              assertEquals(nx, r.size1());
              assertEquals(mx * bm, r.size2());
              Function f = ca.Function("f", new[] { x }, new[] { r });
              Function f_ref = ca.Function("f", new[] { x }, new[] { r_ref });
              checkfunctionLight(f, f_ref, dms(DM.rand(nx, mx)));
            }
          }
        }
      }
    }
    foreach (long nx in new long[] { 2, 1, 0 }) {
      foreach (long mx in new long[] { 3, 1, 0 }) {
        SX x = SX.sym("x", nx, mx);
        foreach (long ny in new long[] { 3, 1, 0 }) {
          foreach (long my in new long[] { 4, 1, 0 }) {
            foreach (long bm in new long[] { 1, 3, 0 }) {
              if (bm != 1 && mx == 0) continue;
              SX y = ca.mtimes(ca.mtimes(DM.rand(ny, nx), x), DM.rand(mx, my));
              assertEquals(ny, y.size1());
              assertEquals(my, y.size2());
              DM dx = DM.rand(nx, mx * bm);

              SX r = ca.jtimes(y, x, ca.SX(dx));
              var parts = new List<SX>();
              foreach (DM e in ca.horzsplit(dx, Math.Max(1, mx)))
                parts.Add(ca.reshape(ca.mtimes(ca.jacobian(y, ca.vec(x)), ca.SX(ca.vec(e))), y.size1(), y.size2()));
              SX r_ref = parts.Count == 0 ? ca.SX(1, 0) : ca.horzcat(parts);

              assertEquals(ny, r.size1());
              assertEquals(my * bm, r.size2());
              Function f = ca.Function("f", new[] { x }, new[] { r });
              Function f_ref = ca.Function("f", new[] { x }, new[] { r_ref });
              checkfunctionLight(f, f_ref, dms(DM.rand(nx, mx)));
            }
          }
        }
      }
    }
  }

  void powZeroSX() {
    SX x = SX.sym("x", Sparsity.diag(4));
    assertTrue(ca.pow(x, 0).is_dense());
    assertFalse(ca.pow(x, 0.3).is_dense());
    assertFalse(ca.pow(x, 1).is_dense());
    assertFalse(ca.pow(x, 2).is_dense());
    assertFalse(ca.pow(x, 2.3).is_dense());
    assertTrue(ca.pow(x, -2).is_dense());
    SX p = SX.sym("p");

    assertTrue(ca.pow(x, p).is_dense());

    SX s = SX.sym("x");

    checkarray(ca.evalf(ca.substitute(ca.vertcat(ca.pow(s, 2), ca.pow(s, 1), ca.pow(s, 0)), s, ca.SX(0))), ca.vertcat(0, 0, 1));
    checkarray(ca.evalf(ca.substitute(ca.pow(s, ca.vertcat(2, 1, 0)), s, ca.SX(0))), ca.vertcat(0, 0, 1));
  }

  public void test_norms() {
    for (int k = 0; k < 2; k++) {
      int kk = k;
      Func<SX, SX> fs = x => kk == 0 ? ca.vertcat(x[0], ca.SX(1, 1), x[1..5]) : ca.diag(x);
      Func<MX, MX> fm = x => kk == 0 ? ca.vertcat(x[0], ca.MX(1, 1), x[1..5]) : ca.diag(x);
      for (int n = 0; n < 3; n++) {
        Console.WriteLine("norm " + n);
        SX xs = SX.sym("x", 5);
        SX es = fs(ca.sin(xs));
        es = ca.log(n == 0 ? ca.norm_fro(es) : n == 1 ? ca.norm_1(es) : ca.norm_inf(es));
        Function fsx = ca.Function("f", new[] { xs }, new[] { es });

        MX xm = MX.sym("x", 5);
        MX em = fm(ca.sin(xm));
        em = ca.log(n == 0 ? ca.norm_fro(em) : n == 1 ? ca.norm_1(em) : ca.norm_inf(em));
        Function fmx = ca.Function("f", new[] { xm }, new[] { em });
        foreach (DM args in new[] { ca.vertcat(-0.3, 0.3, 0.21, 0.17, 0), ca.vertcat(-1, -2, -7, -7, -7) }) checkfunction(fsx, fmx, dms(args));
        checkCodegen(fmx, dms(ca.vertcat(-0.3, 0.3, 0.21, 0.17, 0)));
      }
    }
  }

  public void test_contains() {
    MX x = MX.sym("x", 2, 2);
    MX y = MX.sym("y");
    MX z = MX.sym("z", 3, 3);

    MX e = x * y;
    assertTrue(ca.contains(L(x, y, z), x));
    assertFalse(ca.contains(L(x, y), z));
    assertTrue(ca.contains_any(L(x, y), L(y, z)));
    assertFalse(ca.contains_all(L(x, y), L(y, z)));
    assertTrue(ca.contains_any(L(x, y), L(x, y)));
    assertTrue(ca.contains_all(L(x, y), L(x, y)));
    assertTrue(ca.contains(L(e, x), e));
  }

  public void test_issue4041() {
    MX x = MX.sym("x");
    MX y = MX.sym("y");

    Function f = ca.Function("f", new[] { x, y }, new[] { ca.pow(x, 2) }, opts(("never_inline", true)));

    MX expr = f.Call(x, y);

    assertFalse(ca.depends_on(expr, y));
    assertTrue(ca.contains(ca.symvar(expr), y));

    f = ca.Function("f", new[] { x, y }, new[] { ca.pow(x, 2) }, opts(("always_inline", true)));

    expr = f.Call(x, y);

    assertFalse(ca.depends_on(expr, y));
    assertFalse(ca.contains(ca.symvar(expr), y));
  }

  public void test_symmetric_jacobian() {
    MX x = MX.sym("x", 4);
    MX u = MX.sym("u", 2);

    MX h = x[1] * 5 * x[2];

    long nh = 1;
    MX lam_h = MX.sym("lam_h", nh, 1);

    MX ux = ca.vertcat(u, x);
    MX adj_ux = ca.densify(ca.jtimes(h, ux, lam_h, true));
    ca.jacobian(adj_ux, ux, opts(("symmetric", true)));

    MX adj_ux2 = ca.jtimes(h, ux, lam_h, true);
    assertInException("Symmetry exploitation", () => ca.jacobian(adj_ux2, ux, opts(("symmetric", true))));
  }

  public void test_pow_zero() {
    powZeroSX();
    powZeroMX();
  }

  public void test_pow() {
    MX x = MX.sym("x");
    DM r = ca.evalf(ca.substitute(ca.pow(x, ca.vertcat(2, 1, 0)), x, ca.MX(0)));
    checkarray(r, ca.vertcat(0, 0, 1));

    MX y = ca.MX(4, 1);
    assertEquals(4L, ca.pow(y, 0).nnz());
    y = ca.MX(4, 1);
    y[1] = x;
    assertEquals(4L, ca.pow(y, 0).nnz());
    y = ca.MX(4, 1);
    y[1] = 0;
    assertEquals(4L, ca.pow(y, 0).nnz());
    y = ca.MX(4, 1);
    y[1] = 1;
    assertEquals(4L, ca.pow(y, 0).nnz());
  }

  // test_sum: python-only (numpy np.sum dispatch on casadi types)

  public void test_jtimes_empty() {
    assertInException("Ambiguous", () => ca.jtimes(ca.SX(5), SX.sym("x", 2, 0), ca.SX(DM.ones(2, 5))));
    assertInException("Ambiguous", () => ca.jtimes(ca.SX(DM.ones(4, 0)), SX.sym("x", 2, 3), ca.SX(DM.ones(4, 6)), true));

    jtimesEmptySX();
    jtimesEmptyMX();
  }

  public void test_issue2934() {
    MX M = MX.sym("X", 5, 5);
    MX Y = MX.sym("y", 2, 2);

    MX m = ca.MX(M);

    m.set(Y, false, idx(0, 0), idx(0, 0));

    Console.WriteLine(m);
    Function f = ca.Function("f", new[] { M, Y }, new[] { m });
    Console.WriteLine(f.call(new[] { M, Y }, true, false)[0]);
    assertTrue(f.call(new[] { M * 2, Y }, true, false)[0].ToString().Contains("((2.*X)[-1, -1, -1, 0] = y)"));
  }

  public void test_issue2935() {
    MX M = ca.MX(DM.rand(5, 5));
    MX m = M;

    MX Y = MX.sym("y", 2, 2);
    Sparsity E = Y.sparsity();
    MX y = Y;

    m.set(y, false, idx(0, 0), idx(0, 0));

    MX e = ca.cos(m);
    e = ca.dot(e, e);

    MX eb = MX.sym("eb");
    MX ge = ca.jtimes(e, Y, eb, true);

    Function f = ca.Function("f", new[] { eb, Y }, new[] { ge });

    MX e1 = f.call(new[] { eb, Y }, true, false)[0];
    MX e2 = f.call(new[] { eb, Y + 1e-300 }, true, false)[0];

    Console.WriteLine(e1);
    Console.WriteLine(e2);

    Function f1 = ca.Function("f1", new[] { eb, Y }, new[] { ge });
    Function f2 = ca.Function("f2", new[] { eb, Y }, new[] { e2 });

    List<DM> args = dms(ca.vertcat(1), DM.rand(E.size1(), E.size2()));

    Console.WriteLine(str(f1.call(args)));
    Console.WriteLine(str(f2.call(args)));

    assertTrue(f1.call(args)[0].sparsity().is_equal(f2.call(args)[0].sparsity()));
  }

  public void test_extract() {
    MX x = MX.sym("x", 2);
    MX y = MX.sym("y", 2);

    MX u = ca.vertcat(x, y);
    DM.rng(1);
    DM u0 = DM.rand(4, 1);
    Function f = ca.Function("f", new[] { x, y }, new[] { ca.sin(x + y) }, opts(("never_inline", true)));

    Function g = ca.Function("g", new[] { x }, new[] { ca.cos(ca.sum(x) * 3) }, opts(("never_inline", true)));

    MX e = g.Call(f.Call(x / y, x * 2 * y) * x) / (x + y);

    var (vexpr, v, vdef) = ca.extract(L(e), opts(("lift_shared", false), ("lift_calls", true)));

    Console.WriteLine(str(vexpr));
    Console.WriteLine(str(v));
    Console.WriteLine(str(vdef));

    assertEquals(5, v.Count);
    assertEquals(v.Count, vdef.Count);
    assertEquals(1, vexpr.Count);
    Function f_ref = ca.Function("f", new[] { u }, new[] { e });
    Console.WriteLine("e " + e);
    List<MX> e_reconstructed = ca.substitute_inplace(v, vdef, L(e)).Item2;
    Console.WriteLine(str(e_reconstructed));
    Function fr = ca.Function("f", new[] { u }, e_reconstructed);
    checkarray(f_ref.Call(u0), fr.Call(u0));

    Console.WriteLine(e);

    var si = ca.substitute_inplace(v, vdef, vexpr);
    Console.WriteLine("(" + str(si.Item1) + ", " + str(si.Item2) + ")");

    MX vexprc = ca.vertcat(vexpr);
    MX vc = ca.vertcat(v);
    MX vdefc = ca.vertcat(vdef);

    Function Jf_ref = ca.Function("J_ref", new[] { u }, new[] { ca.jacobian(e, u) });

    MX J_v_u = ca.solve(ca.MX(DM.eye(vdefc.size1())) - ca.jacobian(vdefc, vc), ca.jacobian(vdefc, u));
    MX J = ca.jacobian(vexprc, u) + ca.mtimes(ca.jacobian(vexprc, vc), J_v_u);

    J = ca.substitute_inplace(v, vdef, L(J)).Item2[0];

    Function Jf = ca.Function("J2", new[] { u }, new[] { J });

    checkarray(Jf_ref.Call(u0), Jf.Call(u0));
  }

  public void test_const_folding_on_the_fly() {
    MX x = MX.sym("x");

    foreach (int a in new[] { 2, 7 }) {
      foreach (int b in new[] { 2, 7 }) {
        MX @ref = (a * b) * x;
        if (a == 7 && b == 7) {
          assertEquals(@ref.ToString(), (a * (b * x)).ToString());
          assertEquals(@ref.ToString(), ((b * x) * a).ToString());
        } else {
          assertFalse(@ref.ToString() == (a * (b * x)).ToString());
          assertFalse(@ref.ToString() == ((b * x) * a).ToString());
        }
      }
    }
  }

  public void test_cse_reuse() {
    for (int X = 0; X < 2; X++) {
      MX x = MX.sym("x");
      MX y = ca.sin(MX.sym("y"));
      MX z = MX.sym("z");

      MX e = x + y / z;
      MX r = ca.vertcat(ca.cse(e), y * 2);
      Console.WriteLine(r);
      assertTrue(r.ToString().Contains("@1=sin(y)"));
    }
  }

  public void test_is_slice2() {
    int m = 13;
    int n = 11;

    MX x = MX.sym("x", m, n);
    DM.rng(1);
    DM A = DM.rand(m, n);

    foreach (int step_col in new[] { 1, 2, 3 }) {
      foreach (int step_row in new[] { 1, 2, 3 }) {
        var nz = new List<long>();
        var cslice = new List<double>();
        var rslice = new List<double>();
        for (int col = 3; col < 7; col += step_col) cslice.Add(col);
        for (int row = 1; row < 9; row += step_row) rslice.Add(row);
        foreach (double col in cslice)
          foreach (double row in rslice) nz.Add((long)(row + m * col));

        assertTrue(ca.is_slice2(nz));

        for (int j = 0; j < nz.Count; j++) {
          var nz_mod = new List<long>(nz);
          long val = j == 0 ? nz[0] - 1 : (long)Math.Floor((nz[j] + nz[j - 1]) / 2.0);
          nz_mod.Insert(j, val);
          assertFalse(ca.is_slice2(nz_mod));
        }
        DM ri = rslice.ToArray();
        DM ci = cslice.ToArray();
        MX r = x.get(ri, ci);
        assertTrue(r.ToString().Contains(";"));
        Function f = ca.Function("f", new[] { x }, new[] { r });
        checkarray(f.Call(A), A.get(ri, ci));
      }
    }
  }

  public void test_is_diff_sparsity_propagation() {
    foreach (bool isMX in new[] { true, false }) {
      bool[] is_diff_in = { true, true, false };
      bool[] is_diff_out = { true, false };
      var o = opts(("never_inline", true), ("is_diff_in", is_diff_in), ("is_diff_out", is_diff_out));

      Function F, F2;
      if (isMX) {
        MX x = MX.sym("x"), u = MX.sym("u"), q = MX.sym("q");
        F = ca.Function("F", new[] { x, u, q }, new[] { ca.sin(u + x + q), u * q }, o);
        F2 = ca.Function("F", new[] { x, u, q }, new[] { ca.sin(u + x + q), u * q });
      } else {
        SX x = SX.sym("x"), u = SX.sym("u"), q = SX.sym("q");
        F = ca.Function("F", new[] { x, u, q }, new[] { ca.sin(u + x + q), u * q }, o);
        F2 = ca.Function("F", new[] { x, u, q }, new[] { ca.sin(u + x + q), u * q });
      }
      var Fs = new List<Function> { F };
      if (RunSlow) {
        List<DM> @in = dms(ca.vertcat(1), ca.vertcat(2), ca.vertcat(3));
        checkCodegen(F, @in);
        Fs.Add(lastExternal);

        checkCodegen(F2, @in, opts(), 1e-15, "c89", false, true);
        Fs.Add(ca.external("F", lastLib, opts(("is_diff_in", is_diff_in), ("is_diff_out", is_diff_out))));

        checkCodegen(F2, @in);
        Fs.Add(ca.external("F", lastLib, opts(("is_diff_in", is_diff_in), ("is_diff_out", is_diff_out))));
      }

      foreach (Function Fi in Fs) {
        foreach (bool isMX2 in new[] { true, false }) {
          foreach (int ad_weight_sp in new[] { 0, 1 }) {
            Function H;
            if (isMX2) {
              MX X0 = MX.sym("X0"), U0 = MX.sym("U0"), Q0 = MX.sym("Q0");
              List<MX> r = Fi.call(new[] { X0 * 5, U0 * 6, Q0 * 7 });
              H = ca.Function("H", new[] { ca.vertcat(X0, U0, Q0) }, new[] { r[0] * 2, r[1] * 3 }, opts(("ad_weight_sp", ad_weight_sp)));
            } else {
              SX X0 = SX.sym("X0"), U0 = SX.sym("U0"), Q0 = SX.sym("Q0");
              List<SX> r = Fi.call(new[] { X0 * 5, U0 * 6, Q0 * 7 });
              H = ca.Function("H", new[] { ca.vertcat(X0, U0, Q0) }, new[] { r[0] * 2, r[1] * 3 }, opts(("ad_weight_sp", ad_weight_sp)));
            }
            assertEquals(2L, H.sparsity_jac(0, 0).nnz());
            assertEquals(0L, H.sparsity_jac(0, 1).nnz());
          }
        }
      }
    }
  }

  public void test_is_nonnegative_binary() {
    MX x = MX.sym("x");
    MX y = MX.sym("y");
    foreach (MX e in new[] { x < y, x <= y, ca.eq(x, y), ca.ne(x, y), ca.logic_and(x, y), ca.logic_or(x, y), ca.hypot(x, y), ca.sqrt(x) })
      assertTrue(ca.is_equal(ca.abs(e), e));
    foreach (MX e in new[] { x + y, x * y, x - y }) assertFalse(ca.is_equal(ca.abs(e), e));
  }

  public void test_chained_inequality_kept() {
    MX y = MX.sym("y");
    foreach (MX e in new[] { ca.le(0, y <= 1), ca.ge(1 >= y, 0) }) {
      assertFalse(e.is_constant());
      Function f = ca.Function("f", new[] { y }, new[] { e }).transform(new[] { new object[] { "simplify", 0, "ref_count" } });
      assertTrue(f.call(new[] { y }, true, false)[0].is_one());
    }
  }

  public void test_linear_slice_assign() {
    MX Y0 = MX.sym("Y", 3, 3);
    MX Y = ca.MX(Y0);
    Y[..] = 1;
    checkarray(ca.Function("f", new[] { Y0 }, new[] { Y }).Call(DM.zeros(3, 3)), DM.ones(3, 3));
    MX Z0 = MX.sym("Z", 3, 3);
    MX Z = ca.MX(Z0);
    Z[4] = 1;
    Function f = ca.Function("f", new[] { Z0 }, new[] { Z });
    checkarray(f.Call(DM.zeros(3, 3)), ca.DM(new double[,] { { 0, 0, 0 }, { 0, 1, 0 }, { 0, 0, 0 } }));
  }

  public void test_bspline_checks() {
    MX x = MX.sym("x");
    double[][] knots = { new double[] { 0, 0, 0, 0, 1, 1, 1, 1 } };
    assertInException("Expected 4 coefficients, but got 2", () => ca.bspline(x, ca.vertcat(1, 2), knots, new long[] { 3 }, 1, opts()));
    assertInException("Unknown option 'foo'", () => ca.bspline(x, ca.vertcat(1, 2, 3, 4), knots, new long[] { 3 }, 1, opts(("foo", 1))));
  }

  public void test_user_errors() {
    MX x = MX.sym("x", 2);
    MX y = MX.sym("y");
    assertInException("Expected a column vector with 3 nonzeros", () => ca.MX(Sparsity.lower(2), ca.vertcat(1, 2)));
    assertInException("Expected a primitive of dimension", () => ca.vertcat(x, y).join_primitives(dms(DM.ones(3, 1), ca.DM(1))));
    assertInException("Expected 2 arguments, but got 1", () => (x[0] * y).eval_mx(L(y)));
    assertInException("Not a binary operation", () => MX.binary(999, y, y));
    MX A = MX.sym("A", 4);
    assertInException("Length of B must match dim_b",
      () => MX.einstein(A, A, new long[] { 2, 2 }, new long[] { 3, 2 }, new long[] { 2 }, new long[] { -1, -2 }, new long[] { -2, -3 }, new long[] { -1 }));
    assertInException("A must be a dense vector", () => MX.einstein(MX.sym("A", 2, 2), A, new long[] { 2, 2 }, new long[] { 2, 2 },
      new long[] { 2 }, new long[] { -1, -2 }, new long[] { -2, -3 }, new long[] { -1 }));
    assertInException("Expected 2 query vectors, but got 1",
      () => MX.interpn_linear(L(ca.MX(ca.vertcat(0, 1)), ca.MX(ca.vertcat(0, 1))), MX.sym("v", 4), L(y)));
    assertInException("Grid must be a dense vector", () => ca.low(MX.sym("v", Sparsity.lower(2)), 0.5));
    assertInException("Expected 2 row indices, but got 3", () => x.enlarge(4, 1, new long[] { 0, 1, 2 }, new long[] { 0 }));
    double[] k = { 0, 0, 0, 0, 1, 1, 1, 1 };
    assertInException("multiple of the number of dimensions",
      () => MX.bspline_dual(new double[] { 0.5, 0.5, 0.1 }, new[] { k, k }, new long[] { 3, 3 }));
  }

  public void test_unary_op_check() {
    MX y = MX.sym("y");
    foreach (long op in new long[] { ca.OP_ADD, 999 })
      assertInException("Not a unary operation", () => MX.unary(op, y));
  }

  public void test_print_operator_args() {
    MX x = MX.sym("x");
    MX y = MX.sym("y");
    assertInException("Expected 2 argument names, but got 1", () => ca.print_operator(x + y, new[] { "a" }));
    assertEquals("sin(a)", ca.print_operator(ca.sin(x), new[] { "a", "b" }));
  }

  public void test_set_precision() {
    double pi_val = 3.141592653589793239;
    MX mx_scalar = ca.MX(pi_val);
    MX mx_dm = ca.MX(ca.vertcat(1.234567890123, 9.876543210987));
    try {
      assertEquals("3.14159", mx_scalar.ToString());
      DM.set_precision(12);
      assertEquals("3.14159265359", mx_scalar.ToString());
      assertEquals("[1.23456789012, 9.87654321099]", mx_dm.ToString());
      DM.set_precision(4);
      DM.set_scientific(true);
      assertTrue(mx_scalar.ToString().ToLower().Contains("e"));
    } finally {
      DM.set_precision(6);
      DM.set_scientific(false);
    }
  }
}

}

using System;
using System.Collections;
using System.Collections.Generic;
using System.Diagnostics;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Text.Json;
using System.Text.RegularExpressions;

namespace CasADi.Tests {

/// <summary>Port of test/java Function2Test (test/python/function.py, second part, from test_jit_directory)</summary>
public class Function2Test : CasadiTestCase {

  static readonly string[] XS = { "SX", "MX" };

  // ---------- helpers ----------

  static Dictionary<string, object> opts(params object[] kv) {
    var m = new Dictionary<string, object>();
    for (int i = 0; i < kv.Length; i += 2) m[(string)kv[i]] = kv[i + 1];
    return m;
  }

  /// <summary>A json-like list</summary>
  static List<object> L(params object[] x) => x.ToList();

  static List<object> jl(params double[] x) => x.Select(e => (object)e).ToList();

  static double tol(int digits) => Math.Pow(10, -digits);

  static string read(string path) => File.ReadAllText(path).Replace("\r\n", "\n");

  static void write(string path, string s) => File.WriteAllText(path, s);

  static string tempDir() {
    string d = Path.Combine(Path.GetTempPath(), "casadi_cs_test_" + Guid.NewGuid().ToString("N"));
    Directory.CreateDirectory(d);
    return d;
  }

  static void rmtree(string p) {
    if (Directory.Exists(p)) Directory.Delete(p, true);
  }

  static int countPrefix(string dir, string prefix) =>
    Directory.Exists(dir) ? Directory.GetFileSystemEntries(dir).Count(e => Path.GetFileName(e).StartsWith(prefix)) : 0;

  static void removePrefix(string dir, string prefix) {
    if (!Directory.Exists(dir)) return;
    foreach (string e in Directory.GetFiles(dir))
      if (Path.GetFileName(e).StartsWith(prefix)) File.Delete(e);
  }

  /// <summary>Python json.loads: Dictionary/List/double/string/bool/null trees</summary>
  static object json(string s) {
    using var doc = JsonDocument.Parse(s);
    return conv(doc.RootElement);
  }

  static object conv(JsonElement e) {
    switch (e.ValueKind) {
      case JsonValueKind.Object: {
        var m = new Dictionary<string, object>();
        foreach (var p in e.EnumerateObject()) m[p.Name] = conv(p.Value);
        return m;
      }
      case JsonValueKind.Array: return e.EnumerateArray().Select(conv).ToList();
      case JsonValueKind.String: return e.GetString();
      case JsonValueKind.Number: return e.GetDouble();
      case JsonValueKind.True: return true;
      case JsonValueKind.False: return false;
      default: return null;
    }
  }

  /// <summary>Deep equality of json-like trees (numbers compared as double)</summary>
  static bool jeq(object a, object b) {
    if (a == null || b == null) return a == null && b == null;
    if (a is string || b is string) return Equals(a, b);
    if (a is bool || b is bool) return Equals(a, b);
    if (a is IDictionary da && b is IDictionary db) {
      if (da.Count != db.Count) return false;
      foreach (DictionaryEntry en in da)
        if (!db.Contains(en.Key) || !jeq(en.Value, db[en.Key])) return false;
      return true;
    }
    if (a is IList la && b is IList lb) {
      if (la.Count != lb.Count) return false;
      for (int i = 0; i < la.Count; i++)
        if (!jeq(la[i], lb[i])) return false;
      return true;
    }
    if (a is IConvertible && b is IConvertible && !(a is IList) && !(b is IList))
      return Convert.ToDouble(a, CultureInfo.InvariantCulture) == Convert.ToDouble(b, CultureInfo.InvariantCulture);
    return Equals(a, b);
  }

  static string jstr(object o) => o == null ? "null" : o is IDictionary d
    ? "{" + string.Join(", ", d.Keys.Cast<object>().Select(k => k + ": " + jstr(d[k]))) + "}"
    : o is IList l && !(o is string) ? "[" + string.Join(", ", l.Cast<object>().Select(jstr)) + "]" : o.ToString();

  static void assertJ(object expected, object actual) =>
    Check.True(jeq(expected, actual), "expected " + jstr(expected) + ", got " + jstr(actual));

  /// <summary>Navigate parsed json: string keys for objects, int indices (negative from the end) for lists</summary>
  static object at(object o, params object[] keys) {
    foreach (object k in keys) {
      if (k is string s) {
        o = ((IDictionary<string, object>)o)[s];
      } else {
        var l = (IList)o;
        int i = (int)k;
        o = l[i < 0 ? l.Count + i : i];
      }
    }
    return o;
  }

  static List<object> jsonLines(string path) =>
    read(path).Split('\n').Where(line => line.Length > 0).Select(json).ToList();

  static DM dmOf(object jsonList) => ((IList)jsonList).Cast<object>().Select(x => (double)x).ToArray();

  static IDictionary<string, object> graphBundle(string path) {
    string text = read(path);
    if (path.EndsWith(".html")) {
      Match m = Regex.Match(text, "<script id=\"graph-data\" type=\"application/json\">(.*?)</script>", RegexOptions.Singleline);
      assertTrue(m.Success);
      text = m.Groups[1].Value;
    }
    return (IDictionary<string, object>)json(text);
  }

  static Function unpack(object source) => ca.StringDeserializer((string)source).unpack_function();

  static void assertOutputs(IEnumerable<string> included, IEnumerable<string> excluded, Action r) {
    string o = captureStdout(r);
    foreach (string s in included) assertTrue(o.Contains(s), "expected '" + s + "' in output:\n" + o);
    foreach (string s in excluded) assertFalse(o.Contains(s), "did not expect '" + s + "' in output:\n" + o);
  }

  static List<DM> randInputs(Function f) {
    var r = new List<DM>();
    for (long i = 0; i < f.n_in(); i++) r.Add(DM.rand(f.sparsity_in(i)));
    return r;
  }

  static void checkfunctionLight(Function trial, Function solution, List<DM> inputs) => checkfunctionLight(trial, solution, inputs, 9);

  static void checkfunctionLight(Function trial, Function solution, List<DM> inputs, int digits) {
    assertEquals(solution.n_out(), trial.n_out());
    assertEquals(solution.n_in(), trial.n_in());
    for (int rep = 0; rep < 2; rep++) {
      List<DM> a = trial.call(inputs), b = solution.call(inputs);
      for (int i = 0; i < trial.n_out(); i++) {
        if (trial.sparsity_out(i).is_empty() || solution.sparsity_out(i).is_empty()) continue;
        checkarray(a[i], b[i], tol(digits), "output " + i);
      }
    }
  }

  static void checkfunctionLight(Function trial, Function solution, IDictionary<string, DM> inputs) {
    var inp = new List<DM>();
    for (long i = 0; i < trial.n_in(); i++) inp.Add(inputs.TryGetValue(trial.name_in(i), out DM v) ? v : ca.DM(0.0));
    checkfunctionLight(trial, solution, inp);
  }

  DM uniformFor(Sparsity sp) {
    double[] v = new double[sp.nnz()];
    for (int k = 0; k < v.Length; k++) v[k] = random.NextDouble();
    return ca.DM(sp, v.Length == 0 ? ca.DM(0, 1) : (DM)v);
  }

  /// <summary>checkfunction with outputs compared to tol, derivatives to tolSens, and uniform [0,1) seeds as in helpers.py</summary>
  void checkfunction(Function f, Function g, List<DM> inputs, double tol, double tolSens) {
    checkfunctionLight(f, g, inputs, (int)Math.Round(-Math.Log10(tol)));
    List<DM> a = f.call(inputs), b = g.call(inputs);
    var ja = inputs.Concat(a).ToList();
    var jb = inputs.Concat(b).ToList();
    List<DM> va = f.jacobian().call(ja), vb = g.jacobian().call(jb);
    for (int i = 0; i < va.Count; i++) checkarray(va[i], vb[i], tolSens, "jacobian block " + i);
    var fwd = new List<DM>(ja);
    var adj = new List<DM>(ja);
    for (long i = 0; i < f.n_in(); i++) fwd.Add(uniformFor(f.sparsity_in(i)));
    for (long i = 0; i < f.n_out(); i++) adj.Add(uniformFor(f.sparsity_out(i)));
    List<DM> fa = f.forward(1).call(fwd), fb = g.forward(1).call(fwd);
    for (int i = 0; i < fa.Count; i++) checkarray(fa[i], fb[i], tolSens, "forward " + i);
    List<DM> ra = f.reverse(1).call(adj), rb = g.reverse(1).call(adj);
    for (int i = 0; i < ra.Count; i++) checkarray(ra[i], rb[i], tolSens, "reverse " + i);
  }

  static void checkfunctionIdentical(Function trial, Function solution) {
    trial.generate("f.c");
    string[] a = read("f.c").Split('\n');
    solution.generate("f.c");
    string[] b = read("f.c").Split('\n');
    for (int i = 0; i < Math.Min(a.Length, b.Length); i++) assertEquals(b[i], a[i]);
    assertEquals(b.Length, a.Length);
  }

  /// <summary>Inputs from DMs and numbers</summary>
  static List<DM> ins(params object[] x) =>
    x.Select(e => e is DM d ? d : ca.DM(Convert.ToDouble(e, CultureInfo.InvariantCulture))).ToList();

  static DM sparseDM(double[,] rows) => ca.sparsify((DM)rows);

  /// <summary>Options for checkCodegen, after helpers.py check_codegen</summary>
  sealed class CG {
    public Dictionary<string, object> opts = new Dictionary<string, object>();
    public string std = "c89";
    public bool main, mainOutputCheck = true, withExternal = true, withJacSparsity, withForward, debug;
    public int digits = 15;
    public List<string> extralibs = new List<string>(), extraOptions = new List<string>(), definitions = new List<string>();
    public Dictionary<string, object> externalOpts = new Dictionary<string, object>();

    public CG Std(string s) { std = s; return this; }
    public CG Main() { main = true; return this; }
    public CG NoMainCheck() { mainOutputCheck = false; return this; }
    public CG NoExternal() { withExternal = false; return this; }
    public CG Opt(string k, object v) { opts[k] = v; return this; }
    public CG Digits(int d) { digits = d; return this; }
    public CG Lib(string l) { extralibs.Add(l); return this; }
    public CG Flag(params string[] f) { extraOptions.AddRange(f); return this; }
    public CG Def(string d) { definitions.Add(d); return this; }
    public CG JacSparsity() { withJacSparsity = true; return this; }
    public CG Forward() { withForward = true; return this; }
    public CG Debug() { debug = true; return this; }
    public CG ExternalOpts(Dictionary<string, object> o) { externalOpts = o; return this; }
  }

  static string includeDir() {
    string inc = GlobalOptions.getCasadiIncludePath();
    if (GlobalOptions.getCasadiPath().Length == 0 || inc.Length == 0) Check.Skip("casadi include path not set (build tree)");
    return inc;
  }

  static string q(string s) => "'" + s.Replace("'", "'\\''") + "'";

  /// <summary>Runs cmd through sh, optionally redirecting stdin/stdout from/to files; returns the exit code</summary>
  static int sh(string cmd, string inFile, string outFile) {
    if (inFile != null) cmd += " < " + q(inFile);
    if (outFile != null) cmd += " > " + q(outFile);
    var psi = new ProcessStartInfo("sh") { UseShellExecute = false };
    psi.ArgumentList.Add("-c");
    psi.ArgumentList.Add(cmd);
    string ld = Environment.GetEnvironmentVariable("LD_LIBRARY_PATH") ?? "";
    psi.Environment["LD_LIBRARY_PATH"] = GlobalOptions.getCasadiPath() + ":" + ld;
    using var p = Process.Start(psi);
    p.WaitForExit();
    return p.ExitCode;
  }

  static string[] run(params string[] cmd) {
    var psi = new ProcessStartInfo(cmd[0]) { UseShellExecute = false, RedirectStandardOutput = true, RedirectStandardError = true };
    foreach (string a in cmd.Skip(1)) psi.ArgumentList.Add(a);
    using var p = Process.Start(psi);
    var errTask = p.StandardError.ReadToEndAsync();
    string o = p.StandardOutput.ReadToEnd();
    string e = errTask.Result;
    p.WaitForExit();
    assertEquals(0, p.ExitCode);
    return new[] { o, e };
  }

  static bool which(string exe) => sh("command -v " + exe + " > /dev/null", null, null) == 0;

  string compileCommand(string name, string source, CG o, bool shared) {
    string libdir = GlobalOptions.getCasadiPath(), inc = includeDir();
    string libs = string.Concat(o.extralibs.Select(l => l.StartsWith("-") || l.Contains(".") ? " " + l : " -l" + l));
    string defs = string.Concat(o.definitions.Select(d => " -D" + d));
    return "gcc -pedantic -std=" + o.std + " -fPIC " + (shared ? "-shared" : "") + " -Wall -Werror -Wextra -I" + inc + " -I" + inc
        + "/highs -Wno-dangling-pointer -Wno-unknown-warning-option -Wno-unknown-pragmas -Wno-long-long -Wno-unused-parameter "
        + (o.debug ? "-O0 -g" : "-O3") + defs + " " + source + " -o " + name + (shared ? ".so" : "") + " -L" + libdir + " -Wl,-rpath,"
        + libdir + " -Wl,-rpath,." + (shared ? "" : " -lm") + libs + " " + string.Join(" ", o.extraOptions);
  }

  /// <summary>helpers.py check_codegen: generate C, compile, compare against an external; a no-op without RunSlow</summary>
  Dictionary<string, object> checkCodegen(Function F, List<DM> inputs, CG o) {
    var ret = new Dictionary<string, object>();
    if (!RunSlow) return ret;
    string name = "codegen_" + DateTime.Now.Ticks.ToString("x") + random.Next(1 << 30).ToString("x");
    var cgOpts = new Dictionary<string, object>(o.opts);
    if (o.main) cgOpts["main"] = true;
    var cg = ca.CodeGenerator(name, cgOpts);
    cg.add(F, o.withJacSparsity);
    if (o.withForward) cg.add(F.forward(1), o.withJacSparsity);
    cg.generate();
    string libname = "./" + name + ".so";
    assertEquals(0, sh(compileCommand(name, name + ".c", o, true), null, null), "compilation failed");
    Function F2 = o.withExternal ? ca.external(F.name(), libname, o.externalOpts) : null;
    if (o.main) {
      assertEquals(0, sh(compileCommand(name, name + ".c", o, false), null, null), "compilation failed");
      F.generate_in(F.name() + "_in.txt", inputs);
      int rc = sh("./" + name + " " + F.name(), F.name() + "_in.txt", F.name() + "_out.txt");
      assertEquals(0, rc, "main return code");
    }
    List<DM> outp = F.call(inputs);
    if (o.main && o.mainOutputCheck) {
      List<DM> outputs = F.generate_out(F.name() + "_out.txt");
      for (int i = 0; i < F.n_out(); i++) checkarray(outp[i], outputs[i], tol(o.digits), "main output " + i);
    }
    ret["libname"] = libname;
    if (F2 == null) return ret;
    List<DM> out2 = F2.call(inputs);
    for (int i = 0; i < F.n_out(); i++) checkarray(outp[i], out2[i], tol(o.digits), "external output " + i);
    ret["F"] = F2;
    return ret;
  }

  string externalLibname;

  Function compileExternal(string name, string source, Dictionary<string, object> opts, CG o) {
    assertEquals(0, sh(compileCommand(name, source, o, true), null, null), "compilation failed");
    externalLibname = "./" + name + ".so";
    return ca.external(name, externalLibname, opts);
  }

  static void checkSerialize(Function F, List<DM> inputs) {
    Function F2 = Function.deserialize(F.serialize(opts("debug", true)));
    List<DM> a = F.call(inputs), b = F2.call(inputs);
    for (int i = 0; i < F.n_out(); i++) checkarray(a[i], b[i], 1e-16, "serialize " + i);
  }

  static Regex vizConfig = new Regex("<script id=\"viz-config\" type=\"application/json\">(.*?)</script>", RegexOptions.Singleline);

  // ---------- tests ----------

  public void test_jit_directory() {
    Check.Assume(RunSlow, "slow");
    message(CasadiMeta.feature_list());
    bool haveMkstemps = CasadiMeta.compiler_flags().Contains("HAVE_MKSTEMPS");
    bool haveSimpleMkstemps = CasadiMeta.compiler_flags().Contains("HAVE_SIMPLE_MKSTEMPS");
    bool createDirs = !CasadiMeta.feature_list().Contains("ghc-filesystem");
    string absTempDir = Path.GetTempPath();
    string relTempDir = "jit_dir_test";
    rmtree(relTempDir);
    Directory.CreateDirectory("park");
    if (createDirs) {
      Directory.CreateDirectory(Path.Combine(relTempDir, "foo"));
      Directory.CreateDirectory("foo");
      Directory.CreateDirectory("/tmp/foo");
    }
    foreach (string tempDir in new[] { "./", relTempDir, absTempDir }) {
      GlobalOptions.setTempWorkDir(tempDir);
      foreach (string directory in new[] { "", "foo", absTempDir }) {
        message("directory " + directory);
        MX x = MX.sym("x");
        var o = opts("jit", true, "compiler", "shell", "jit_options", opts("verbose", true, "directory", directory));
        if (!haveMkstemps && !haveSimpleMkstemps) {
          assertInException("custom temporary directory", () => ca.Function("f", new[] { x }, new[] { ca.pow(x, 2) }, o));
          continue;
        }
        string dir = Path.IsPathRooted(directory) ? directory : Path.Combine(tempDir, directory);
        removePrefix(dir, "tmp_casadi");
        var f = ca.Function("f", new[] { x }, new[] { ca.pow(x, 2) }, o);
        assertTrue(countPrefix(dir, "tmp_casadi") > 1);
        // cannot chdir("park") for the ghc-filesystem variant
        f.Dispose();
        GC.Collect();
        assertEquals(0, countPrefix(dir, "tmp_casadi"));
      }
      GlobalOptions.setTempWorkDir("./");
    }
  }

  public void test_custom_jacobian() {
    MX x = MX.sym("x");
    MX p = MX.sym("p");
    var J = ca.Function("jac_Q", new[] { x, p, ca.MX(1, 1), ca.MX(1, 1) }, new[] { x * Math.PI, ca.MX(1.0), ca.MX(1.0), ca.MX(1.0) },
        new[] { "x", "p", "out_r", "out_s" }, new[] { "jac_r_x", "jac_r_p", "jac_s_x", "jac_s_p" });
    var f = ca.Function("Q", new[] { x, p }, new[] { ca.pow(x, 2), x * p * 2 }, new[] { "x", "p" }, new[] { "r", "s" },
        opts("custom_jacobian", J, "jac_penalty", 0));
    J = null;
    List<MX> g = f.call(new[] { x, p }, false, true);
    var JJ = ca.Function("JJ", new[] { x, p }, new[] { ca.jacobian(g[0] + x * 3, x) });
    checkarray(JJ.Call(5, 2), 5 * Math.PI + 3);
  }

  public void test_dump() {
    MX x = MX.sym("x", Sparsity.lower(3));
    MX y = MX.sym("y", 0, 0);
    MX z = MX.sym("z", 2, 2);
    foreach (string fmt in new[] { "txt", "mtx" }) {
      var f = ca.Function("f", new[] { x, y, z }, new[] { x * 2, z * 2 }, new[] { "x", "y", "z" }, new[] { "a", "c" },
          opts("dump", true, "dump_in", true, "dump_out", true, "dump_format", fmt));
      var inp = dms(sparseDM(new double[,] { { 1, 0, 0 }, { 2, 4, 0 }, { 7, 8, 9 } }), ca.DM(), ca.DM(new double[,] { { 1, 3 }, { 4, 5 } }));
      List<DM> outp = f.call(inp);
      Function F = Function.load("f.casadi");
      DM xNum = DM.from_file("f.000000.in.x." + fmt);
      DM yNum = DM.from_file("f.000000.in.y." + fmt);
      DM zNum = DM.from_file("f.000000.in.z." + fmt);
      checkarray(xNum, inp[0]);
      checkarray(yNum, inp[1]);
      checkarray(zNum, inp[2]);
      DM aNum = DM.from_file("f.000000.out.a." + fmt);
      DM cNum = DM.from_file("f.000000.out.c." + fmt);
      checkarray(aNum, outp[0]);
      checkarray(cNum, outp[1]);
      F.generate_in("test_in.txt", dms(xNum, yNum, zNum));
      F.generate_out("test_out.txt", F.call(dms(xNum, yNum, zNum)));
      DM X = DM.from_file("f.000000.in.txt");
      DM A = DM.from_file("f.000000.out.txt");
      DM Xr = DM.from_file("test_in.txt");
      DM Ar = DM.from_file("test_out.txt");
      checkarray(Xr, X);
      checkarray(Ar, A);
    }
    if (RunSlow && CasadiMeta.feature_list().Contains("ghc-filesystem")) {
      foreach (string d in new[] { "dump_fn", "dump_fn_codegen" }) rmtree(d);
      MX x2 = MX.sym("x", Sparsity.lower(3));
      MX z2 = MX.sym("z", 2, 2);
      var f = ca.Function("f", new[] { x2, z2 }, new[] { x2 * 2, z2 * 2 }, new[] { "x", "z" }, new[] { "a", "c" },
          opts("dump_in", true, "dump_out", true, "dump_dir", "dump_fn"));
      var inp = dms(sparseDM(new double[,] { { 1, 0, 0 }, { 2, 4, 0 }, { 7, 8, 9 } }), ca.DM(new double[,] { { 1, 3 }, { 4, 5 } }));
      f.call(inp);
      checkCodegen(f, inp, new CG().Std("c99").Main().Opt("dump_dir_suffix", "_codegen"));
      foreach (string name in new[] { "f.000000.in.txt", "f.000000.in.x.mtx", "f.000000.in.z.mtx", "f.000000.out.a.mtx", "f.000000.out.c.mtx", "f.000000.out.txt" })
        assertEquals(read(Path.Combine("dump_fn", name)), read(Path.Combine("dump_fn_codegen", name)));
      foreach (string d in new[] { "dump_fn", "dump_fn_codegen" }) rmtree(d);
    }
  }

  public void test_dump_trace() {
    message("Instruction dumps preserve intermediate values and share the dump counter");
    foreach (string X in XS) {
      string dir = tempDir();
      try {
        var o = opts("dump_trace", true, "dump_in", true, "dump_out", true, "dump_dir", dir);
        Function f;
        if (X == "SX") {
          SX x = SX.sym("x", 2), p = SX.sym("p");
          f = ca.Function("traced", new[] { x, p }, new[] { x * p + ca.pow(x, 2), ca.sum1(x) }, o);
        } else {
          MX x = MX.sym("x", 2), p = MX.sym("p");
          f = ca.Function("traced", new[] { x, p }, new[] { x * p + ca.pow(x, 2), ca.sum1(x) }, o);
        }
        double[][] xvals = { new double[] { 2, -3 }, new double[] { 1.1, 5 } };
        double[] pvals = { 4, -2 };
        for (int count = 0; count < 2; count++) {
          double[] xval = xvals[count];
          double pval = pvals[count];
          var result = new List<DM>();
          string captured = captureStdout(() => result.AddRange(f.call(dms(xval, pval))));
          assertEquals("", captured);
          string prefix = Path.Combine(dir, string.Format("traced.{0:D6}", count));
          List<object> records = jsonLines(prefix + ".trace.jsonl");
          assertJ("casadi_trace", at(records, 0, "format"));
          assertJ(1.0, at(records, 0, "version"));
          assertJ((double)count, at(records, 0, "dump_id"));
          assertJ(opts("event", "inputs", "values", L(jl(xval), jl(pval))), records[1]);
          assertJ(opts("event", "end", "status", 0.0), at(records, -1));
          var last = (IList)at(records, -2, "values");
          for (int i = 0; i < Math.Min(result.Count, last.Count); i++) checkarray(result[i].nonzeros().ToArray(), dmOf(last[i]));
          checkarray(ca.vertcat(result), DM.from_file(prefix + ".out.txt"));
          checkarray(ca.vertcat(xval, pval), DM.from_file(prefix + ".in.txt"));
          List<object> steps = records.GetRange(2, records.Count - 4);
          assertEquals(2 * f.n_instructions(), (long)steps.Count);
          var checkedOps = new HashSet<long>();
          long add = ca.OP_ADD, mul = ca.OP_MUL, sq = ca.OP_SQ;
          for (int k = 0; k < f.n_instructions(); k++) {
            object before = steps[2 * k], after = steps[2 * k + 1];
            long op = f.instruction_id(k);
            assertJ((double)k, at(before, "instruction"));
            assertJ((double)k, at(after, "instruction"));
            assertJ((double)op, at(before, "op"));
            assertJ((double)op, at(after, "op"));
            assertJ("inputs", at(before, "phase"));
            assertJ("outputs", at(after, "phase"));
            if (op == add || op == mul || op == sq) {
              DM a = dmOf(at(before, "values", 0));
              DM expected;
              if (op == sq) {
                expected = a * a;
              } else {
                DM b = dmOf(at(before, "values", 1));
                expected = op == add ? a + b : a * b;
              }
              checkarray(expected, dmOf(at(after, "values", 0)));
              checkedOps.Add(op);
            }
          }
          assertTrue(checkedOps.SetEquals(new[] { add, mul, sq }));
        }
        Function restored = Function.deserialize(f.serialize());
        restored.reset_dump_count();
        restored.call(dms(ca.vertcat(2, -3), 4));
        assertJ(0.0, at(jsonLines(Path.Combine(dir, "traced.000000.trace.jsonl")), 0, "dump_id"));
        f.change_option("dump_trace", false);
        f.call(dms(ca.vertcat(2, -3), 4));
        assertFalse(File.Exists(Path.Combine(dir, "traced.000002.trace.jsonl")));
        f.change_option("dump_trace", true);
        f.call(dms(ca.vertcat(2, -3), 4));
        assertTrue(File.Exists(Path.Combine(dir, "traced.000003.trace.jsonl")));
        var generator = ca.CodeGenerator("traced.c");
        generator.add(f);
        string tracedCode = generator.dump();
        f.change_option("dump_trace", false);
        generator = ca.CodeGenerator("traced.c");
        generator.add(f);
        assertEquals(tracedCode, generator.dump());
      } finally {
        rmtree(dir);
      }
    }
  }

  public void test_dump_trace_options() {
    message("Reconstructed XFunctions preserve instruction dumping");
    foreach (string X in XS) {
      string dir = tempDir();
      try {
        var o = opts("dump_trace", true, "dump_dir", dir);
        Function f;
        Action jitTrace;
        if (X == "SX") {
          SX x = SX.sym("x");
          f = ca.Function("original", new[] { x }, new[] { x * x }, o);
          jitTrace = () => ca.Function("jit_trace", new[] { x }, new[] { x * x }, opts("jit", true, "dump_trace", true));
        } else {
          MX x = MX.sym("x");
          f = ca.Function("original", new[] { x }, new[] { x * x }, o);
          jitTrace = () => ca.Function("jit_trace", new[] { x }, new[] { x * x }, opts("jit", true, "dump_trace", true));
        }
        Function expanded = f.expand("expanded");
        Function derivative = f.forward(1);
        expanded.Call(3);
        derivative.Call(3, 9, 1);
        foreach (var (fn, v) in new[] { (expanded, 9.0), (derivative, 6.0) }) {
          List<object> records = jsonLines(Path.Combine(dir, fn.name() + ".000000.trace.jsonl"));
          assertJ(L(L(v)), at(records, -2, "values"));
        }
        assertInException("dump_trace", jitTrace);
      } finally {
        rmtree(dir);
      }
    }
  }

  public void test_dump_trace_sparse_call() {
    message("Instruction dumps handle sparse matrices and unused call outputs");
    foreach (string X in XS) {
      string dir = tempDir();
      try {
        var innerOpts = opts("never_inline", true, "dump_trace", true, "dump_dir", dir);
        var o = opts("dump_trace", true, "dump_dir", dir);
        Function f;
        if (X == "SX") {
          SX v = SX.sym("v", 3);
          var inner = ca.Function("inner", new[] { v }, new[] { v + 1, v * 3 }, innerOpts);
          SX x = SX.sym("x", Sparsity.diag(3));
          SX a = inner.call(new[] { ca.vertcat(x[0, 0], x[1, 1], x[2, 2]) })[0];
          f = ca.Function("outer", new[] { x }, new[] { x * 2, a, ca.SX(0, 0), x }, o);
        } else {
          MX v = MX.sym("v", 3);
          var inner = ca.Function("inner", new[] { v }, new[] { v + 1, v * 3 }, innerOpts);
          MX x = MX.sym("x", Sparsity.diag(3));
          MX a = inner.call(new[] { ca.vertcat(x[0, 0], x[1, 1], x[2, 2]) })[0];
          f = ca.Function("outer", new[] { x }, new[] { x * 2, a, ca.MX(0, 0), x }, o);
        }
        f.Call(ca.diag(ca.vertcat(1, 2, 3)));
        List<object> records = jsonLines(Path.Combine(dir, "outer.000000.trace.jsonl"));
        assertJ(L(jl(1, 2, 3)), at(records, 1, "values"));
        var vals = (List<object>)at(records, -2, "values");
        assertJ(L(jl(2, 4, 6), jl(2, 3, 4)), vals.GetRange(0, 2));
        assertTrue(vals[2] == null || jeq(L(), vals[2]));
        assertJ(jl(1, 2, 3), vals[3]);
        double call = ca.OP_CALL;
        var calls = records.Where(r => ((IDictionary<string, object>)r).TryGetValue("op", out object op) && jeq(call, op)).ToList();
        assertEquals(2, calls.Count, "calls");
        assertTrue(((IList)at(calls, 1, "values")).Contains(null));
        List<object> innerRecords = jsonLines(Path.Combine(dir, "inner.000000.trace.jsonl"));
        assertJ(opts("event", "end", "status", 0.0), at(innerRecords, -1));
      } finally {
        rmtree(dir);
      }
    }
  }

  public void test_dump_trace_error() {
    message("Failed evaluations leave a readable partial trace");
    string dir = tempDir();
    try {
      MX x = MX.sym("x");
      MX y = x.attachAssert(ca.gt(x, 0.0), "positive input required");
      var f = ca.Function("checked", new[] { x }, new[] { y + 1 }, opts("dump_trace", true, "dump_dir", dir));
      assertInException("positive input required", () => f.Call(-1));
      List<object> records = jsonLines(Path.Combine(dir, "checked.000000.trace.jsonl"));
      assertJ(opts("event", "error"), at(records, -1));
      checkarray(f.Call(2), 3);
      records = jsonLines(Path.Combine(dir, "checked.000001.trace.jsonl"));
      assertJ(opts("event", "end", "status", 0.0), at(records, -1));
    } finally {
      rmtree(dir);
    }
  }

  public void test_dump_trace_nonfinite() {
    message("Instruction dumps encode non-finite values as valid JSON");
    foreach (string X in XS) {
      string dir = tempDir();
      try {
        var o = opts("dump_trace", true, "dump_dir", dir, "print_instructions", true, "print_canonical", true);
        Function f;
        if (X == "SX") {
          SX x = SX.sym("x");
          f = ca.Function("nonfinite", new[] { x }, new[] { 1.0 / x, -1.0 / x, ca.sqrt(x - 1) }, o);
        } else {
          MX x = MX.sym("x");
          f = ca.Function("nonfinite", new[] { x }, new[] { 1.0 / x, -1.0 / x, ca.sqrt(x - 1) }, o);
        }
        string captured = captureStdout(() => f.Call(0));
        assertTrue(captured.Contains("nonfinite:"), captured);
        List<object> records = jsonLines(Path.Combine(dir, "nonfinite.000000.trace.jsonl"));
        assertJ(L(L("inf"), L("-inf"), L("nan")), at(records, -2, "values"));
      } finally {
        rmtree(dir);
      }
    }
  }

  public void test_export_graph() {
    message("HTML graph export embeds the serialized Function, consistent with dump_trace");
    foreach (string X in XS) {
      string dir = tempDir();
      try {
        var o = opts("dump_trace", true, "dump_dir", dir);
        Function f;
        if (X == "SX") {
          SX x = SX.sym("x");
          f = ca.Function("graph", new[] { x }, new[] { x * x + x }, new[] { "state" }, new[] { "result" }, o);
        } else {
          MX x = MX.sym("x");
          f = ca.Function("graph", new[] { x }, new[] { x * x + x }, new[] { "state" }, new[] { "result" }, o);
        }
        string path = Path.Combine(dir, "graph.html");
        f.export_graph(path);
        var bundle = graphBundle(path);
        assertJ("casadi_viz", bundle["format"]);
        assertJ(1.0, bundle["version"]);
        assertJ("function", bundle["view"]);
        assertJ(CasadiMeta.version(), bundle["casadi_version"]);
        assertFalse(bundle.ContainsKey("nodes"));
        Function restored = unpack(bundle["source"]);
        assertEquals("graph", restored.name());
        assertEquals(X + "Function", restored.class_name());
        assertTrue(restored.name_in().SequenceEqual(new[] { "state" }));
        assertTrue(restored.name_out().SequenceEqual(new[] { "result" }));
        var ra = new List<long>();
        var fa = new List<long>();
        for (long k = 0; k < restored.n_instructions(); k++) ra.Add(restored.instruction_id(k));
        for (long k = 0; k < f.n_instructions(); k++) fa.Add(f.instruction_id(k));
        assertTrue(fa.SequenceEqual(ra));
        f.Call(3);
        List<object> records = jsonLines(Path.Combine(dir, "graph.000000.trace.jsonl"));
        assertJ(restored.name(), at(records, 0, "function"));
        assertJ(restored.class_name(), at(records, 0, "type"));
        foreach (object record in records) {
          var r = (IDictionary<string, object>)record;
          if (!r.ContainsKey("instruction")) continue;
          assertJ((double)restored.instruction_id((long)(double)r["instruction"]), r["op"]);
        }
        assertJ(L(jl(12)), at(records, -2, "values"));
        string html = read(path);
        Match m = vizConfig.Match(html);
        assertTrue(m.Success);
        var config = (IDictionary<string, object>)json(m.Groups[1].Value);
        string[] v = CasadiMeta.version().Split('.');
        assertJ("https://unpkg.com/@casadi/casadi-viz@" + v[0] + "." + v[1] + "/dist/index.js", config["viewer_url"]);
        assertTrue(config.ContainsKey("runtime"));
        assertTrue(config["runtime"] == null);
        Function ff = f;
        assertInException("Invalid export_graph direction", () => ff.export_graph(path, opts("direction", "sideways")));
        assertInException("Unknown export_graph option", () => ff.export_graph(path, opts("unknown", true)));
      } finally {
        rmtree(dir);
      }
    }
  }

  public void test_export_graph_calls() {
    message("Graph exports keep called Functions, by identity, and sparse metadata");
    foreach (string X in XS) {
      string dir = tempDir();
      try {
        string path = Path.Combine(dir, "calls.html");
        Function f, nested, sparse;
        if (X == "SX") {
          SX x = SX.sym("x");
          var pair = ca.Function("pair", new[] { x }, new[] { x + 1, x - 1 }, opts("never_inline", true));
          List<SX> ab = pair.call(new[] { x });
          f = ca.Function("calls", new[] { x }, new[] { ab[0] * ab[1] });
          var first = ca.Function("inner", new[] { x }, new[] { x + 1 }, opts("never_inline", true));
          var second = ca.Function("inner", new[] { x }, new[] { x - 2 }, opts("never_inline", true));
          var middle = ca.Function("middle", new[] { x }, new[] { first.Call(x) + second.Call(x) }, opts("never_inline", true));
          nested = ca.Function("nested", new[] { x }, new[] { middle.Call(x) + first.Call(x + 3) });
          SX xs = SX.sym("x", Sparsity.lower(3));
          sparse = ca.Function("sparse_graph", new[] { xs }, new[] { xs * 2 });
        } else {
          MX x = MX.sym("x");
          var pair = ca.Function("pair", new[] { x }, new[] { x + 1, x - 1 }, opts("never_inline", true));
          List<MX> ab = pair.call(new[] { x });
          f = ca.Function("calls", new[] { x }, new[] { ab[0] * ab[1] });
          var first = ca.Function("inner", new[] { x }, new[] { x + 1 }, opts("never_inline", true));
          var second = ca.Function("inner", new[] { x }, new[] { x - 2 }, opts("never_inline", true));
          var middle = ca.Function("middle", new[] { x }, new[] { first.Call(x) + second.Call(x) }, opts("never_inline", true));
          nested = ca.Function("nested", new[] { x }, new[] { middle.Call(x) + first.Call(x + 3) });
          MX xs = MX.sym("x", Sparsity.lower(3));
          sparse = ca.Function("sparse_graph", new[] { xs }, new[] { xs * 2 });
        }
        f.export_graph(path, opts("direction", "TB"));
        var bundle = graphBundle(path);
        assertJ("TB", bundle["direction"]);
        Function restored = unpack(bundle["source"]);
        assertTrue(restored.get_function().SequenceEqual(new[] { "pair" }));
        checkarray(restored.Call(3), f.Call(3));
        restored = unpack(((IDictionary<string, object>)json(nested.export_graph()))["source"]);
        var names = restored.find_functions().Select(g => g.name()).OrderBy(n => n, StringComparer.Ordinal).ToList();
        assertTrue(names.SequenceEqual(new[] { "inner", "inner", "middle" }), string.Join(",", names));
        checkarray(restored.Call(1.5), nested.Call(1.5));
        sparse.export_graph(path);
        restored = unpack(graphBundle(path)["source"]);
        assertTrue(restored.sparsity_in(0).is_equal(Sparsity.lower(3)));
      } finally {
        rmtree(dir);
      }
    }
  }

  public void test_export_graph_helpers() {
    message("Graph exports accept any serializable Function, with its helper Functions");
    string dir = tempDir();
    try {
      MX x = MX.sym("x", 2), p = MX.sym("p");
      var nlp = new Dictionary<string, MX> {
        ["x"] = x, ["p"] = p, ["f"] = ca.pow(x[0] - p, 2) + ca.pow(x[1], 2), ["g"] = x[0] + x[1]
      };
      Function solver = ca.nlpsol("solver", "sqpmethod", nlp, opts("qpsol", "qrqp", "print_header", false, "print_iteration", false, "print_time", false,
          "qpsol_options", opts("print_iter", false, "print_header", false)));
      SX s = SX.sym("s");
      Function mapped = ca.Function("f", new[] { s }, new[] { ca.sin(s) * s }).map(4);
      foreach (Function f in new[] { solver, mapped }) {
        foreach (string extension in new[] { ".html", ".casadi_viz" }) {
          string path = Path.Combine(dir, "helpers" + extension);
          f.export_graph(path);
          var bundle = graphBundle(path);
          assertJ("function", bundle["view"]);
          Function restored = unpack(bundle["source"]);
          assertEquals(f.class_name(), restored.class_name());
          assertTrue(f.get_function().SequenceEqual(restored.get_function()));
        }
        assertInException("DOT export requires an SXFunction or MXFunction", () => f.export_graph(Path.Combine(dir, "helpers.dot")));
      }
      assertTrue(solver.get_function().Contains("nlp_hess_l"));
    } finally {
      rmtree(dir);
    }
  }

  public void test_export_graph_expressions() {
    message("Expression and Function exports default to expression and function views");
    string directory = tempDir();
    try {
      string path = Path.Combine(directory, "expressions.html");
      foreach (string X in XS) {
        Function f;
        Action bad;
        if (X == "SX") {
          SX x = SX.sym("actual_x", 2), y = SX.sym("actual_y", 2);
          ca.export_graph(x + y, path);
          assertJ("expression", graphBundle(path)["view"]);
          foreach (var expressions in new[] { new[] { x + y, x - y }, new[] { x[0], y[0] }, new[] { ca.SX(7.0), SX.zeros(0, 2) } }) {
            ca.export_graph(expressions, path);
            assertJ("expression", graphBundle(path)["view"]);
          }
          ca.export_graph(new[] { x + y, x - y }, path, opts("view", "function"));
          assertJ("function", graphBundle(path)["view"]);
          f = ca.Function("renamed", new[] { x, y }, new[] { x + y, x - y }, new[] { "left", "right" }, new[] { "sum", "difference" });
          bad = () => ca.export_graph(x, path, opts("view", "other"));
        } else {
          MX x = MX.sym("actual_x", 2), y = MX.sym("actual_y", 2);
          ca.export_graph(x + y, path);
          assertJ("expression", graphBundle(path)["view"]);
          foreach (var expressions in new[] { new[] { x + y, x - y }, new[] { x[0], y[0] }, new[] { ca.MX(7.0), MX.zeros(0, 2) } }) {
            ca.export_graph(expressions, path);
            assertJ("expression", graphBundle(path)["view"]);
          }
          ca.export_graph(new[] { x + y, x - y }, path, opts("view", "function"));
          assertJ("function", graphBundle(path)["view"]);
          f = ca.Function("renamed", new[] { x, y }, new[] { x + y, x - y }, new[] { "left", "right" }, new[] { "sum", "difference" });
          bad = () => ca.export_graph(x, path, opts("view", "other"));
        }
        f.export_graph(path);
        var bundle = graphBundle(path);
        assertJ("function", bundle["view"]);
        Function restored = unpack(bundle["source"]);
        assertTrue(restored.name_in().SequenceEqual(new[] { "left", "right" }));
        assertTrue(restored.name_out().SequenceEqual(new[] { "sum", "difference" }));
        HashSet<string> syms;
        if (X == "SX") {
          syms = ca.symvar(ca.veccat(restored.sx_in())).Select(v => v.ToString()).ToHashSet();
          assertTrue(syms.SetEquals(new[] { "actual_x_0", "actual_x_1", "actual_y_0", "actual_y_1" }), string.Join(",", syms));
        } else {
          syms = ca.symvar(ca.veccat(restored.mx_in())).Select(v => v.ToString()).ToHashSet();
          assertTrue(syms.SetEquals(new[] { "actual_x", "actual_y" }), string.Join(",", syms));
        }
        f.export_graph(path, opts("view", "expression"));
        assertJ("expression", graphBundle(path)["view"]);
        assertInException("Invalid export_graph view", bad);
      }
      ca.export_graph(new SX[0], path);
    } finally {
      rmtree(directory);
    }
  }

  public void test_export_graph_string() {
    message("Filename-free graph exports return the same bundle as file exports");
    foreach (string X in XS) {
      string directory = tempDir();
      try {
        string path = Path.Combine(directory, "graph.casadi_viz");
        Function f;
        var toString = new List<Func<Dictionary<string, object>, string>>();
        var toFile = new List<Action<string>>();
        if (X == "SX") {
          SX x = SX.sym("x", 2);
          var child = ca.Function("child", new[] { x }, new[] { x + 1 }, opts("never_inline", true));
          f = ca.Function("root", new[] { x }, new[] { child.Call(x) / x });
          toString.Add(o => o == null ? ca.export_graph(x) : ca.export_graph(x, o));
          toFile.Add(pa => ca.export_graph(x, pa));
          foreach (var e in new List<SX>[] { new List<SX> { x }, new List<SX> { x[0] + x[1], x[0] - x[1] }, new List<SX> { ca.SX(2.0) }, new List<SX>() }) {
            toString.Add(o => o == null ? ca.export_graph(e) : ca.export_graph(e, o));
            toFile.Add(pa => ca.export_graph(e, pa));
          }
        } else {
          MX x = MX.sym("x", 2);
          var child = ca.Function("child", new[] { x }, new[] { x + 1 }, opts("never_inline", true));
          f = ca.Function("root", new[] { x }, new[] { child.Call(x) / x });
          toString.Add(o => o == null ? ca.export_graph(x) : ca.export_graph(x, o));
          toFile.Add(pa => ca.export_graph(x, pa));
          foreach (var e in new List<MX>[] { new List<MX> { x }, new List<MX> { x[0] + x[1], x[0] - x[1] }, new List<MX> { ca.MX(2.0) }, new List<MX>() }) {
            toString.Add(o => o == null ? ca.export_graph(e) : ca.export_graph(e, o));
            toFile.Add(pa => ca.export_graph(e, pa));
          }
        }
        string def = f.export_graph();
        assertEquals(def, f.export_graph(opts()));
        assertJ("function", ((IDictionary<string, object>)json(def))["view"]);
        foreach (var o in new[] { opts(), opts("include_functions", false, "view", "expression", "direction", "LR") }) {
          string text = f.export_graph(o);
          f.export_graph(path, o);
          assertEquals(text + "\n", read(path));
          var bundle = (IDictionary<string, object>)json(text);
          assertJ(o.GetValueOrDefault("include_functions", true), bundle["include_functions"]);
          assertJ(o.GetValueOrDefault("direction", "TB"), bundle["direction"]);
          assertJ(((IDictionary<string, object>)json(def))["source"], bundle["source"]);
        }
        assertEquals(def, f.export_graph(opts("viz_js", "not-needed.js")));
        for (int i = 0; i < toString.Count; i++) {
          string text = toString[i](null);
          assertJ("expression", ((IDictionary<string, object>)json(text))["view"]);
          toFile[i](path);
          assertEquals(text + "\n", read(path));
          assertJ("function", ((IDictionary<string, object>)json(toString[i](opts("view", "function"))))["view"]);
        }
        Function ff = f;
        assertInException("Invalid export_graph view", () => ff.export_graph(opts("view", "invalid")));
        assertInException("Unknown export_graph option", () => ff.export_graph(opts("unknown", true)));
      } finally {
        rmtree(directory);
      }
    }
  }

  public void test_export_graph_formats() {
    message("HTML embeds the .casadi_viz bundle; DOT stays a standalone instruction graph");
    foreach (string X in XS) {
      string directory = tempDir();
      try {
        string bas = Path.Combine(directory, "graph");
        Function f, free;
        Action exprExport;
        if (X == "SX") {
          SX x = SX.sym("x", 2);
          var child = ca.Function("child", new[] { x }, new[] { x + 1 }, opts("never_inline", true));
          f = ca.Function("root", new[] { x }, new[] { child.Call(x) / x }, new[] { "a" }, new[] { "result" });
          exprExport = () => ca.export_graph(new[] { x[0] + x[1], x[0] - x[1] }, bas + ".casadi_viz");
          SX fr = SX.sym("free");
          free = ca.Function("free_graph", new[] { x }, new[] { x + fr }, opts("allow_free", true));
        } else {
          MX x = MX.sym("x", 2);
          var child = ca.Function("child", new[] { x }, new[] { x + 1 }, opts("never_inline", true));
          f = ca.Function("root", new[] { x }, new[] { child.Call(x) / x }, new[] { "a" }, new[] { "result" });
          exprExport = () => ca.export_graph(new[] { x[0] + x[1], x[0] - x[1] }, bas + ".casadi_viz");
          MX fr = MX.sym("free");
          free = ca.Function("free_graph", new[] { x }, new[] { x + fr }, opts("allow_free", true));
        }
        f.export_graph(bas + ".html");
        var embedded = graphBundle(bas + ".html");
        f.export_graph(bas + ".casadi_viz", opts("viz_js", "not-needed.js"));
        var bundle = graphBundle(bas + ".casadi_viz");
        assertJ(embedded, bundle);
        assertJ(CasadiMeta.version(), bundle["casadi_version"]);
        foreach (string extension in new[] { ".html", ".casadi_viz" }) {
          f.export_graph(bas + extension, opts("include_functions", false));
          var reduced = graphBundle(bas + extension);
          assertJ(false, reduced["include_functions"]);
          assertJ(bundle["source"], reduced["source"]);
        }
        foreach (string view in new[] { "function", "expression" }) {
          f.export_graph(bas + ".dot", opts("view", view, "viz_js", "not-needed.js"));
          string dot = read(bas + ".dot");
          assertTrue(dot.StartsWith("digraph G {"));
          assertFalse(dot.Contains("<script"));
          assertFalse(dot.Contains("arg0"));
          assertTrue(dot.Contains(":w"));
          assertTrue(dot.Contains(":e"));
          assertTrue(dot.Contains("2-by-1"));
          assertEquals(view == "function", dot.Contains("result"));
          if (which("dot")) {
            string[] rendered = run("dot", "-Tsvg", bas + ".dot");
            assertTrue(rendered[0].Contains("<svg"));
            assertEquals("", rendered[1]);
          }
        }
        foreach (string extension in new[] { "", ".json", ".casadi_graph" }) {
          Function ff = f;
          assertInException("Supported extensions: .html, .dot, .casadi_viz", () => ff.export_graph(bas + extension));
        }
        exprExport();
        assertJ("expression", graphBundle(bas + ".casadi_viz")["view"]);
        free.export_graph(bas + ".casadi_viz");
        Function restored = unpack(graphBundle(bas + ".casadi_viz")["source"]);
        List<string> names = X == "SX" ? restored.free_sx().Select(v => v.ToString()).ToList() : restored.free_mx().Select(v => v.ToString()).ToList();
        assertTrue(names.SequenceEqual(new[] { "free" }), string.Join(",", names));
      } finally {
        rmtree(directory);
      }
    }
  }

  public void test_export_graph_escaping() {
    message("Serialized graphs and optional local Viz.js stay safely embedded in HTML");
    string directory = tempDir();
    try {
      MX x = MX.sym("x");
      string name = "</script><script>alert(\"input\")</script>";
      var f = ca.Function("escaped", new[] { x }, new[] { x }, new[] { name }, new[] { "result" });
      string source = "/* </script> & \"quoted\" \\ source */";
      string runtime = Path.Combine(directory, "viz.js");
      write(runtime, source);
      string path = Path.Combine(directory, "escaped.html");
      f.export_graph(path, opts("viz_js", runtime, "viewer_url", "https://example.com/viewer.js?a=</script>"));
      string html = read(path);
      assertFalse(html.Contains(name));
      Function restored = unpack(graphBundle(path)["source"]);
      Match m = vizConfig.Match(html);
      assertTrue(m.Success);
      assertTrue(json(m.Groups[1].Value) != null);
      assertEquals(name, restored.name_in(0));
    } finally {
      rmtree(directory);
    }
  }

  public void test_print() {
    MX x = MX.sym("x", Sparsity.lower(3));
    MX z = MX.sym("z", 2, 2);
    var f = ca.Function("f", new[] { x, z }, new[] { x * 2, z * 2 }, new[] { "x", "z" }, new[] { "a", "c" },
        opts("print_in", true, "print_out", true, "print_canonical", true));
    var inp = dms(sparseDM(new double[,] { { 1, 0, 0 }, { 2, 4, 0 }, { 7, 8, 9 } }), ca.DM(new double[,] { { 1, 3 }, { 4, 5 } }));
    var result = new List<DM>();
    string vmOutput = captureStdout(() => result.AddRange(f.call(inp)));
    checkarray(result[0], inp[0] * 2);
    checkarray(result[1], inp[1] * 2);
    var vmLines = vmOutput.Trim().Split('\n').Select(line => Regex.Replace(line, @" \(0x[0-9a-f]+\)", "")).ToList();
    assertTrue(vmLines.Count >= 6);
    checkCodegen(f, inp, new CG().Std("c99").Main().NoMainCheck());
    if (RunSlow) {
      var cgLines = read(f.name() + "_out.txt").Trim().Split('\n')
          .Where(line => line.StartsWith("Function") || line.StartsWith("Input") || line.StartsWith("Output")).ToList();
      assertTrue(vmLines.SequenceEqual(cgLines), string.Join("\n", vmLines) + "\nvs\n" + string.Join("\n", cgLines));
    }
  }

  public void test_eval_shapes() {
    MX x = MX.sym("x", Sparsity.lower(3));
    MX y = MX.sym("y", 3, 1);
    MX z = MX.sym("z", 2, 2);
    var f = ca.Function("f", new[] { x, y, z }, new[] { x * 2, y * 2, z * 2 }, opts("default_in", new double[] { 1, 2, 3 }));
    double[,] Lm = { { 1, 0, 0 }, { 2, 4, 0 }, { 7, 8, 9 } };
    double[,] Z = { { 1, 3 }, { 4, 5 } };
    DM l = sparseDM(Lm);
    DM lt2 = l * 2;
    DM y987 = ca.vertcat(9, 8, 7);
    var cases = new List<(List<DM> inp, List<DM> refs, double[] refFlat)>();
    // normal
    cases.Add((ins(l, y987, ca.DM(Z)), dms(lt2, y987 * 2, ca.DM(Z) * 2), new double[] { 1, 2, 7, 4, 8, 9, 9, 8, 7, 1, 4, 3, 5 }));
    // Project
    cases.Add((ins(sparseDM(new double[,] { { 1, 0, 0 }, { 2, 0, 0 }, { 7, 8, 9 } }), y987, ca.DM(Z)),
        dms(sparseDM(new double[,] { { 1, 0, 0 }, { 2, 1e-30, 0 }, { 7, 8, 9 } }) * 2, y987 * 2, ca.DM(Z) * 2),
        new double[] { 1, 2, 7, 0, 8, 9, 9, 8, 7, 1, 4, 3, 5 }));
    // Project away
    cases.Add((ins(sparseDM(new double[,] { { 1, 0, 1 }, { 2, 4, 0 }, { 7, 8, 9 } }), y987, ca.DM(Z)), dms(lt2, y987 * 2, ca.DM(Z) * 2),
        new double[] { 1, 2, 7, 4, 8, 9, 9, 8, 7, 1, 4, 3, 5 }));
    // Scalar expansion
    cases.Add((ins(ca.DM(4), 5, 7),
        dms(sparseDM(new double[,] { { 4, 0, 0 }, { 4, 4, 0 }, { 4, 4, 4 } }) * 2, ca.vertcat(5, 5, 5) * 2, ca.DM(new double[,] { { 7, 7 }, { 7, 7 } }) * 2),
        new double[] { 4, 4, 4, 4, 4, 4, 5, 5, 5, 7, 7, 7, 7 }));
    // Repmat
    cases.Add((ins(l, y987, ca.vertcat(1, 3)), dms(lt2, y987 * 2, ca.DM(new double[,] { { 1, 1 }, { 3, 3 } }) * 2),
        new double[] { 1, 2, 7, 4, 8, 9, 9, 8, 7, 1, 3, 1, 3 }));
    // npar
    cases.Add((ins(l, ca.DM(new double[,] { { 9, 2 }, { 8, 4 }, { 7, 4 } }), ca.DM(Z)),
        dms(sparseDM(new double[,] { { 1, 0, 0, 1, 0, 0 }, { 2, 4, 0, 2, 4, 0 }, { 7, 8, 9, 7, 8, 9 } }) * 2, ca.DM(new double[,] { { 9, 2 }, { 8, 4 }, { 7, 4 } }) * 2,
            ca.DM(new double[,] { { 1, 3, 1, 3 }, { 4, 5, 4, 5 } }) * 2),
        null));
    // growing npar
    DM z8 = ca.DM(new double[,] { { 1, 3, 4, 5, 1, 1, 1, 1 }, { 4, 5, 7, 8, 3, 3, 3, 3 } });
    cases.Add((ins(l, ca.DM(new double[,] { { 9, 2 }, { 8, 4 }, { 7, 4 } }), z8),
        dms(ca.repmat(l, 1, 4) * 2, ca.repmat(ca.DM(new double[,] { { 9, 2 }, { 8, 4 }, { 7, 4 } }), 1, 2) * 2, z8 * 2), null));
    // growing npar
    cases.Add((ins(l, ca.DM(new double[,] { { 9, 2, 3, 4 }, { 8, 4, 7, 8 }, { 7, 4, 9, 10 } }), ca.DM(new double[,] { { 1, 3, 4, 5 }, { 4, 5, 7, 8 } })), null, null));
    // Transpose
    cases.Add((ins(l, y987.T(), ca.DM(Z)), dms(lt2, y987 * 2, ca.DM(Z) * 2), new double[] { 1, 2, 7, 4, 8, 9, 9, 8, 7, 1, 4, 3, 5 }));
    // Null: note default_in applies only to dict style; the last two should be forbidden
    foreach (DM e in new[] { ca.DM(), ca.DM(3, 0), ca.DM(0, 3), ca.DM(0, 5), ca.DM(8, 0) })
      cases.Add((ins(l, e, ca.vertcat(1, 3)), dms(lt2, ca.vertcat(0, 0, 0) * 2, ca.DM(new double[,] { { 1, 1 }, { 3, 3 } }) * 2),
          new double[] { 1, 2, 7, 4, 8, 9, 0, 0, 0, 1, 3, 1, 3 }));
    // Wrong
    cases.Add((ins(l, y987.T(), ca.vertcat(1, 3).T()), null, null));
    cases.Add((ins(sparseDM(new double[,] { { 1, 0, 0 }, { 2, 4, 0 }, { 7, 8, 9 }, { 4, 5, 6 } }), y987, ca.vertcat(1, 3)), null, null));
    cases.Add((ins(l, ca.DM(new double[,] { { 9, 2, 1 }, { 8, 4, 1 }, { 7, 4, 1 } }), ca.DM(new double[,] { { 1, 1, 1, 3 }, { 3, 4, 5, 6 } })), null, null));
    cases.Add((ins(l, y987, ca.DM(new double[,] { { 1, 1, 1 }, { 3, 4, 5 } })), null, null));
    cases.Add((ins(l, ca.vertcat(9, 8, 7, 6), ca.DM(Z)), null, null));

    foreach (var (inp, refs, refFlat) in cases) {
      if (refs == null) {
        assertInException("mismatching shape", () => f.call(inp));
      } else {
        List<DM> res = f.call(inp);
        for (int i = 0; i < f.n_in(); i++) {
          assertTrue(res[i].sparsity().is_equal(refs[i].sparsity()));
          checkarray(res[i], refs[i]);
        }
        if (refFlat == null) {
          assertInException("mismatching shape", () => f.nz_from_out(res));
        } else {
          double[] resFlat = f.nz_from_out(res).ToArray();
          checkarray(resFlat, refFlat.Select(v => 2 * v).ToArray());
          List<DM> resRe = f.nz_to_out(resFlat);
          for (int i = 0; i < f.n_out(); i++) {
            assertTrue(resRe[i].sparsity().is_equal(refs[i].sparsity()));
            checkarray(resRe[i], refs[i]);
          }
          f.generate_out("test.txt", res);
          List<DM> res2 = f.generate_in("test.txt");
          for (int i = 0; i < f.n_out(); i++) {
            assertTrue(res2[i].sparsity().is_equal(refs[i].sparsity()));
            checkarray(res2[i], refs[i]);
          }
        }
      }
      if (refFlat == null) {
        assertInException("mismatching shape", () => f.nz_from_in(inp));
        assertInException("mismatching shape", () => f.generate_in("test.txt", inp));
      } else {
        double[] res = f.nz_from_in(inp).ToArray();
        checkarray(res, refFlat);
        List<DM> inRe = f.nz_to_in(res);
        for (int i = 0; i < f.n_in(); i++) {
          assertTrue(inRe[i].sparsity().is_equal(refs[i].sparsity()));
          checkarray(inRe[i], refs[i] / 2);
        }
        f.generate_in("test.txt", inp);
        List<DM> ins2 = f.generate_in("test.txt");
        for (int i = 0; i < f.n_in(); i++) {
          assertTrue(ins2[i].sparsity().is_equal(refs[i].sparsity()));
          checkarray(ins2[i], refs[i] / 2);
        }
      }
    }
    // Null dicts
    var ind = new Dictionary<string, DM> { ["i0"] = l, ["i2"] = ca.vertcat(1, 3) };
    Dictionary<string, DM> resd = f.call(ind);
    checkarray(resd["o0"], lt2);
    checkarray(resd["o1"], ca.vertcat(2, 2, 2) * 2);
    checkarray(resd["o2"], ca.DM(new double[,] { { 1, 1 }, { 3, 3 } }) * 2);
    checkarray(f.nz_from_in(f.convert_in(ind)).ToArray(), ca.vertcat(1, 2, 7, 4, 8, 9, 2, 2, 2, 1, 3, 1, 3));
  }

  public void test_issue3977() {
    foreach (string X in XS) {
      Function f, g;
      if (X == "SX") {
        SX x = SX.sym("x", 0, 1);
        f = ca.Function("f", new[] { x }, new[] { ca.pow(x, 2) });
        SX x2 = SX.sym("x", 0, 2);
        g = ca.Function("f", new[] { x2 }, new[] { ca.pow(x2, 2) });
      } else {
        MX x = MX.sym("x", 0, 1);
        f = ca.Function("f", new[] { x }, new[] { ca.pow(x, 2) });
        MX x2 = MX.sym("x", 0, 2);
        g = ca.Function("f", new[] { x2 }, new[] { ca.pow(x2, 2) });
      }
      assertEquals((0L, 3L), f.Call(ca.DM(0, 3)).size());
      assertEquals((0L, 1L), f.Call(ca.DM(0, 1)).size());
      assertEquals((0L, 6L), g.Call(ca.DM(0, 6)).size());
    }
    SX xs = SX.sym("x", 3, 1);
    var fs = ca.Function("f", new[] { xs }, new[] { ca.pow(xs, 2) });
    assertEquals((3L, 1L), fs.Call(ca.DM(0, 0)).size());
    assertEquals((3L, 1L), fs.Call(ca.DM(0, 1)).size());
  }

  public void test_convert_in() {
    MX x = MX.sym("x");
    MX y = MX.sym("y");
    MX z = MX.sym("z");
    var f = ca.Function("f", new[] { x, y, z }, new[] { x * 2, y * 2, z * 2, x + y }, new[] { "x", "y", "z" }, new[] { "a", "b", "c", "d" },
        opts("default_in", new double[] { 1, 2, 3 }));
    List<DM> res = f.convert_in(new Dictionary<string, DM> { ["x"] = 4 });
    checkarray(res[0], 4);
    checkarray(res[1], 2);
    checkarray(res[2], 3);
    assertInException("xy", () => f.convert_in(new Dictionary<string, DM> { ["xy"] = 4 }));
    Dictionary<string, DM> resd = f.convert_in(ins(7, 3, 4));
    checkarray(resd["x"], 7);
    checkarray(resd["y"], 3);
    checkarray(resd["z"], 4);
    assertInException("Incorrect", () => f.convert_in(ins(7, 3)));
    res = f.convert_out(new Dictionary<string, DM> { ["b"] = 4 });
    checkarray(res[0], double.NaN);
    checkarray(res[1], 4);
    checkarray(res[2], double.NaN);
    checkarray(res[3], double.NaN);
    assertInException("xy", () => f.convert_out(new Dictionary<string, DM> { ["xy"] = 4 }));
    resd = f.convert_out(ins(1, 2, 3, 4));
    checkarray(resd["a"], 1);
    checkarray(resd["b"], 2);
    checkarray(resd["c"], 3);
    checkarray(resd["d"], 4);
    assertInException("Incorrect", () => f.convert_out(ins(1, 2, 3)));
  }

  public void test_ad_weight_sp() {
    Check.Assume(hasNlpsol("ipopt"), "no nlpsol ipopt");
    MX x = MX.sym("x", 4);
    var f = ca.Function("f", new[] { x }, new[] { x });
    assertOutputs(new[] { "1 forward sweeps needed" }, new[] { "1 reverse sweeps needed" },
        () => ca.jacobian(f.Call(ca.sin(x)), x, opts("helper_options", opts("verbose", true, "ad_weight_sp", 0))));
    assertOutputs(new[] { "1 reverse sweeps needed" }, new[] { "1 forward sweeps needed" },
        () => ca.jacobian(f.Call(ca.sin(x)), x, opts("helper_options", opts("verbose", true, "ad_weight_sp", 1))));
    MX x2 = MX.sym("x", 4);
    var nlp = new Dictionary<string, MX> { ["x"] = x2, ["f"] = ca.dot(x2, x2), ["g"] = ca.sin(x2) };
    assertOutputs(new[] { "1 forward sweeps needed" }, new[] { "1 reverse sweeps needed" },
        () => ca.nlpsol("solver", "ipopt", nlp, opts("common_options", opts("helper_options", opts("verbose", true, "ad_weight_sp", 0)))));
    assertOutputs(new[] { "1 reverse sweeps needed" }, new[] { "1 forward sweeps needed" },
        () => ca.nlpsol("solver", "ipopt", nlp, opts("common_options", opts("helper_options", opts("verbose", true, "ad_weight_sp", 1)))));
    var f2 = ca.Function("f", new[] { x2 }, new[] { x2 }, opts("verbose", true, "ad_weight_sp", 1));
    assertOutputs(new[] { "1 reverse sweeps needed" }, new[] { "1 forward sweeps needed" }, () => ca.jacobian(f2.Call(ca.sin(x2)), x2));
    MX x3 = MX.sym("x", 100);
    var f3 = ca.Function("f", new[] { x3 }, new[] { x3[0] }, opts("verbose", true, "ad_weight_sp", 0));
    assertOutputs(new[] { "2 forward sweeps needed" }, new[] { "1 reverse sweeps needed" }, () => ca.jacobian(f3.Call(ca.sin(x3)), x3));
  }

  public void test_map_exception() {
    MX x = MX.sym("x", 4);
    MX y = MX.sym("y", 4);
    var f = ca.Function("f", new[] { x }, new[] { x + y }, opts("allow_free", true));
    string msg = CasadiMeta.compiler_flags().Contains("CASADI_WITH_THREAD") ? "Evaluation failed" : "since variables [y] are free";
    assertInException(msg, () => {
      Function F = f.map(4, "thread", 2);
      F.Call(3);
    });
  }

  public void test_DM_arg() {
    // C# has no Function(string, IEnumerable<DM>, ...) constructor; DM converts implicitly to SX and MX
    var f = ca.Function("f", new SX[] { ca.DM(0, 1) }, new SX[0]);
    message(f.ToString());
    f = ca.Function("f", new SX[] { ca.DM(0, 1), SX.sym("x") }, new SX[0]);
    message(f.ToString());
    f = ca.Function("f", new MX[] { ca.DM(0, 1), MX.sym("x") }, new MX[0]);
    message(f.ToString());
    object v = ca.vertcat();
    assertTrue(v is DM);
    assertEquals((0L, 1L), ((DM)v).size());
    object h = ca.horzcat();
    assertTrue(h is DM);
    assertEquals((1L, 0L), ((DM)h).size());
    object vc = ca.veccat();
    assertTrue(vc is DM);
    assertEquals((0L, 1L), ((DM)vc).size());
    object d = ca.diagcat();
    assertTrue(d is DM);
    assertEquals((0L, 0L), ((DM)d).size());
    // vcat/hcat/vvcat/dcat: python-only list aliases of vertcat/horzcat/veccat/diagcat
  }

  public void test_nondiff() {
    foreach (string X in XS) {
      DM xn = ca.DM(new double[,] { { 1, 3 }, { 0.1, 2 } });
      DM yn = ca.DM(new double[,] { { 4, 5 }, { 6, 7 } });
      var options = opts("is_diff_in", new[] { true, false }, "is_diff_out", new[] { true, false });
      Function f, F, G;
      if (X == "SX") {
        SX x = SX.sym("x", 2, 2), y = SX.sym("y", 2, 2);
        f = ca.Function("f", new[] { x, y }, new[] { ca.sin(x + y), x * y }, options);
        F = ca.Function("F", new[] { x, y }, f.call(new[] { ca.cos(x), ca.pow(x * y, 2) }), new[] { "x", "y" }, new[] { "z", "zz" }, options);
        var ye = ca.SX(2, 2);
        G = ca.Function("G", new[] { x, ye }, new[] { f.call(new[] { ca.cos(x), ca.pow(x * yn, 2) })[0], ca.SX(f.call(dms(ca.cos(xn), ca.pow(xn * yn, 2)))[1]) },
            new[] { "x", "y" }, new[] { "z", "zz" });
      } else {
        MX x = MX.sym("x", 2, 2), y = MX.sym("y", 2, 2);
        f = ca.Function("f", new[] { x, y }, new[] { ca.sin(x + y), x * y }, options);
        Function fexpanded = f.expand();
        assertEquals(f.ToString().Split(' ')[0], fexpanded.ToString().Split(' ')[0]);
        F = ca.Function("F", new[] { x, y }, f.call(new[] { ca.cos(x), ca.pow(x * y, 2) }), new[] { "x", "y" }, new[] { "z", "zz" }, options);
        var ye = ca.MX(2, 2);
        G = ca.Function("G", new[] { x, ye }, new[] { f.call(new[] { ca.cos(x), ca.pow(x * yn, 2) })[0], ca.MX(f.call(dms(ca.cos(xn), ca.pow(xn * yn, 2)))[1]) },
            new[] { "x", "y" }, new[] { "z", "zz" });
      }
      foreach (Function ff in new[] { F.forward(1).forward(1), F.forward(1).reverse(1), F.reverse(1).forward(1), F.reverse(1).reverse(1) }) {
        string s = ff.ToString().Replace("_y[2x2,0nz]", "foo");
        assertFalse(s.Substring(Math.Min(s.Length, "fwd1_adj1_F:(x[2x2],y[2x2]".Length)).Contains("y"), s);
      }
      checkfunction(F, G, dms(xn, yn));
      var ds = new Func<Function, Function>[] { g => g.forward(1), g => g.reverse(1) };
      foreach (var f1 in ds) {
        Function Gf = f1(G), Ff = f1(F);
        var arg = dms(xn, yn);
        for (long i = 2; i < Gf.n_in(); i++) arg.Add(DM.rand(Gf.sparsity_in(i)));
        List<DM> a = Gf.call(arg), b = Ff.call(arg);
        for (int i = 0; i < a.Count; i++) checkarray(a[i], b[i]);
      }
      foreach (var f1 in ds)
        foreach (var f2 in ds) {
          Function Gf = f1(f2(G)), Ff = f1(f2(F));
          var arg = dms(xn, yn);
          for (long i = 2; i < Gf.n_in(); i++) arg.Add(DM.rand(Gf.sparsity_in(i)));
          List<DM> a = Gf.call(arg), b = Ff.call(arg);
          for (int i = 0; i < a.Count; i++) checkarray(a[i], b[i]);
        }
    }
  }

  public void test_non_diff_call() {
    MX x = MX.sym("x", 5);
    MX y = MX.sym("x", 5, 5);
    var f = ca.Function("f", new[] { x, y }, new[] { ca.mtimes(y, x), ca.mtimes(y.T(), x) },
        opts("never_inline", true, "is_diff_in", new[] { false, true }, "is_diff_out", new[] { false, true }));
    SX xs = SX.sym("x", 5);
    SX ys = SX.sym("x", 5, 5);
    var g = ca.Function("g", new[] { xs, ys }, f.call(new[] { xs, ys }));
    var G = ca.Function("g", new[] { xs, ys }, f.call(new[] { xs, ys }));
    DM.rng(1);
    DM x0 = DM.rand(5);
    DM y0 = DM.rand(5, 5);
    checkfunction(g, G, dms(x0, y0));
  }

  public void test_inline_linear_interpolant() {
    int N = 3, M = 4;
    double[][] dKnots = { new double[N], new double[M] };
    for (int i = 0; i < N; i++) dKnots[0][i] = i / (N - 1.0);
    for (int i = 0; i < M; i++) dKnots[1][i] = i / (M - 1.0);
    double[] dFlat = new double[2 * N * M];
    for (int i = 0; i < dFlat.Length; i++) dFlat[i] = random.NextDouble();
    Function LUTparam = ca.interpolant("name", "linear", new long[] { N, M }, 2, opts("lookup_mode", new[] { "exact", "linear" }));
    var onlyFd = opts("enable_fd", true, "enable_forward", false, "enable_reverse", false);
    var options = opts("ad_weight_sp", -1);
    foreach (var kv in onlyFd) options[kv.Key] = kv.Value;
    options["forward_options"] = onlyFd;
    options["reverse_options"] = onlyFd;
    message(jstr(options));
    Function LUTparamRef = LUTparam.wrap_as_needed(options);
    Function Jref = LUTparamRef.jacobian();
    var inputs = dms(ca.vertcat(0.2, 0.333), dKnots[0].Concat(dKnots[1]).ToArray(), dFlat);
    LUTparam = ca.interpolant("name", "linear", new long[] { N, M }, 2, opts("inline", true, "lookup_mode", new[] { "exact", "linear" }));
    Function J = LUTparam.jacobian();
    checkarray(LUTparam.call(inputs)[0], LUTparamRef.call(inputs)[0]);
    inputs.Add(0.0);
    message("J_ref " + Jref + " " + string.Join(", ", Jref.call(inputs)));
    checkfunction(J, Jref, inputs, 1e-8, 0.1);
  }

  public void test_issue3079() {
    ca.interp1d(new double[] { 1, 2, 3 }, ca.MX(ca.vertcat(4, 5, 6)), new double[] { 1.1 }, "linear");
  }

  // test_functionbuffer: python-only (memoryview-backed FunctionBuffer)

  // test_jit_buffer_eval: python-only (FunctionBuffer evaluation into a numpy memoryview)

  public void test_codegen_inf_nan() {
    Check.Assume(RunSlow, "slow");
    MX x = MX.sym("x");
    var f = ca.Function("F", new[] { x }, new[] { x + double.PositiveInfinity });
    checkCodegen(f, ins(1), new CG().Std("c99").Main());
    x = MX.sym("x");
    f = ca.Function("F", new[] { x }, new[] { x + double.NaN });
    checkCodegen(f, ins(1), new CG().Std("c99").Main());
    x = MX.sym("x");
    f = ca.Function("F", new[] { x }, new[] { x + ca.vertcat(double.PositiveInfinity, double.NaN, double.NegativeInfinity) });
    checkCodegen(f, ins(1), new CG().Std("c99").Main());
  }

  public void test_codegen_with_mem() {
    Check.Assume(RunSlow, "slow");
    MX x = MX.sym("x");
    var f = ca.Function("F", new[] { x }, new[] { x * 3 });
    checkCodegen(f, ins(1), new CG().Main().Opt("with_mem", true).Opt("force_canonical", true).Def("inline=\"\""));
    checkCodegen(f, ins(1), new CG().Main().Opt("with_mem", true).Opt("with_header", true).Opt("force_canonical", true).Def("inline=\"\""));
  }

  public void test_codegen_scalars_bug() {
    Check.Assume(RunSlow, "slow");
    MX x = MX.sym("x");
    MX z = x * 3 / ca.sin(x);
    var f = ca.Function("F", new[] { x }, new[] { z, z / x }, opts("live_variables", false));
    checkCodegen(f, ins(1), new CG().Opt("codegen_scalars", true));
  }

  public void test_bug_codegen_logical() {
    Check.Assume(RunSlow, "slow");
    var a = ca.MX(ca.vertcat(1, 0, 0));
    var b = ca.MX(ca.vertcat(1, 1, 0));
    MX c = ca.logic_or(a, b);
    var f = ca.Function("f", new MX[0], new[] { c });
    checkCodegen(f, ins(), new CG());
  }

  static List<Func<Function>> jitSerializeCases(string jitSerialize) {
    var o = opts("jit", true, "compiler", "shell", "jit_options", opts("verbose", true), "verbose", false, "jit_serialize", jitSerialize);
    var r = new List<Func<Function>>();
    MX x = MX.sym("x");
    r.Add(() => ca.Function("f", new[] { x }, new[] { ca.pow(x - 3, 2) }, o));
    MX x2 = MX.sym("x", 2);
    MX f = x2[0] * x2[0] + x2[1] * x2[1];
    if (ca.has_nlpsol("ipopt")) r.Add(() => ca.nlpsol("solver", "ipopt", new Dictionary<string, MX> { ["x"] = x2, ["f"] = f }, o));
    return r;
  }

  public void test_jit_serialize() {
    Check.Assume(RunSlow, "slow");
    Check.Assume(!System.Runtime.InteropServices.RuntimeInformation.IsOSPlatform(System.Runtime.InteropServices.OSPlatform.OSX), "not on mac");
    var fg = new Function[2];
    var none = new Dictionary<string, DM>();
    var no = new string[0];
    var jt = new[] { "jit_tmp" };
    foreach (var c in jitSerializeCases("source")) {
      assertOutputs(jt, no, () => fg[0] = c());
      fg[0].save("f.casadi");
      assertOutputs(jt, no, () => fg[1] = Function.load("f.casadi"));
      checkfunctionLight(fg[0], fg[1], none);
      fg[1].save("f.casadi");
      assertOutputs(jt, no, () => fg[1] = Function.load("f.casadi"));
      checkfunctionLight(fg[0], fg[1], none);
    }
    foreach (var c in jitSerializeCases("link")) {
      assertOutputs(jt, no, () => fg[0] = c());
      fg[0].save("f.casadi");
      assertOutputs(no, jt, () => fg[1] = Function.load("f.casadi"));
      checkfunctionLight(fg[0], fg[1], none);
      fg[1].save("f.casadi");
      fg[1].Dispose();
      assertOutputs(no, jt, () => fg[1] = Function.load("f.casadi"));
      checkfunctionLight(fg[0], fg[1], none);
      fg[0].Dispose();
      fg[1].Dispose();
      GC.Collect();
      assertInException("No such file", () => Function.load("f.casadi"));
    }
    foreach (var c in jitSerializeCases("embed")) {
      assertOutputs(jt, no, () => fg[0] = c());
      fg[0].save("f.casadi");
      assertOutputs(no, jt, () => fg[1] = Function.load("f.casadi"));
      fg[1].call(none);
      fg[0].Dispose();
      assertOutputs(no, jt, () => fg[1] = Function.load("f.casadi"));
      fg[1].call(none);
      fg[1].save("f.casadi");
      assertOutputs(no, jt, () => fg[1] = Function.load("f.casadi"));
      fg[1].call(none);
      fg[1].Dispose();
      assertOutputs(no, jt, () => fg[1] = Function.load("f.casadi"));
    }
  }

  public void test_map_get_function() {
    MX x = MX.sym("x");
    var g = ca.Function("g", new[] { x }, new[] { ca.pow(x, 2) });
    Function ff = g.map(5);
    assertEquals(1, ff.get_function().Count);
    Function f2 = ff.get_function("f");
    checkfunctionLight(g, f2, ins(3));
  }

  static void assertNames(IList<Function> refs, IList<Function> res) {
    assertEquals(refs.Count, res.Count, "number of functions");
    for (int i = 0; i < refs.Count; i++) assertEquals(refs[i].name(), res[i].name());
  }

  public void test_find_function() {
    foreach (string X in new[] { "MX", "SX" }) {
      Function g, h, i, j, k, intg;
      var ni = opts("never_inline", true);
      if (X == "MX") {
        MX x = MX.sym("x");
        g = ca.Function("g", new[] { x }, new[] { ca.pow(x, 2) }, ni);
        h = ca.Function("h", new[] { x }, new[] { g.Call(ca.pow(x, 2)) }, ni);
        i = ca.Function("i", new[] { x }, new[] { h.Call(ca.pow(x, 2)) + g.Call(ca.pow(x, 2)) }, ni);
        j = ca.Function("j", new[] { x }, new[] { i.Call(ca.pow(x, 2)) + g.Call(ca.pow(x, 2)) }, ni);
        k = ca.Function("k", new[] { x }, new[] { j.Call(ca.pow(x, 2)) }, ni);
        intg = ca.integrator("intg", "rk", new Dictionary<string, MX> { ["x"] = x, ["ode"] = k.Call(x) });
      } else {
        SX x = SX.sym("x");
        g = ca.Function("g", new[] { x }, new[] { ca.pow(x, 2) }, ni);
        h = ca.Function("h", new[] { x }, new[] { g.Call(ca.pow(x, 2)) }, ni);
        i = ca.Function("i", new[] { x }, new[] { h.Call(ca.pow(x, 2)) + g.Call(ca.pow(x, 2)) }, ni);
        j = ca.Function("j", new[] { x }, new[] { i.Call(ca.pow(x, 2)) + g.Call(ca.pow(x, 2)) }, ni);
        k = ca.Function("k", new[] { x }, new[] { j.Call(ca.pow(x, 2)) }, ni);
        intg = ca.integrator("intg", "rk", new Dictionary<string, SX> { ["x"] = x, ["ode"] = k.Call(x) });
      }
      k.generate("f.c");
      assertNames(new[] { j }, k.find_functions(0));
      assertNames(new[] { j, i, g }, k.find_functions(1));
      assertNames(new[] { j, i, h, g }, k.find_functions(2));
      assertNames(new[] { j, i, h, g }, k.find_functions());
      List<Function> res = intg.find_functions();
      message(string.Join(", ", res.Select(e => e.name())));
      assertNames(new[] { intg.get_function("dae"), k, j, i, h, g, intg.get_function("daeF"), intg.get_function("step") }, res);
      Function K = k.map(5);
      assertNames(new[] { k, j, i, h, g }, K.find_functions());
    }
  }

  public void test_post_expand() {
    MX x = MX.sym("x");
    var f = ca.Function("f", new[] { x }, new[] { ca.pow(x, 2) });
    message(f.ToString());
    f = ca.Function("f", new[] { x }, new[] { ca.pow(x, 2) }, opts("post_expand", true));
    message(f.ToString());
    assertTrue(f.is_a("SXFunction"));
    var g = ca.Function("f", new[] { x }, new[] { ca.pow(x, 2) }, opts("post_expand", true, "post_expand_options", opts("print_in", true)));
    message(g.ToString());
    assertTrue(g.is_a("SXFunction"));
    assertOutputs(new[] { "Input 0 (i0): 3" }, new string[0], () => g.Call(3));
  }

  static Dictionary<string, DM> qpInputs() => new Dictionary<string, DM> {
    ["h"] = ca.DM(new double[,] { { 1, 0.2 }, { 0.2, 1 } }), ["g"] = ca.vertcat(1, 2), ["a"] = ca.DM(new double[,] { { 1, 1 } }), ["lba"] = -1, ["uba"] = 1
  };

  static Dictionary<string, Sparsity> qpStructure() => new Dictionary<string, Sparsity> { ["a"] = Sparsity.dense(1, 2), ["h"] = Sparsity.dense(2, 2) };

  static Dictionary<string, MX> withG(Dictionary<string, DM> inp, string key, MX sym) {
    var r = inp.ToDictionary(e => e.Key, e => ca.MX(e.Value));
    r[key] = sym;
    return r;
  }

  static Dictionary<string, SX> withG(Dictionary<string, DM> inp, string key, SX sym) {
    var r = inp.ToDictionary(e => e.Key, e => ca.SX(e.Value));
    r[key] = sym;
    return r;
  }

  static Function wrapper(Function c, Dictionary<string, DM> inputs, string X, string symbol, bool names) {
    if (X == "MX") {
      MX g = MX.sym(symbol, 2);
      MX r = c.call(withG(inputs, symbol, g))["x"];
      return names ? ca.Function("wrapper", new[] { g }, new[] { r }, new[] { symbol }, new[] { "res" }) : ca.Function("wrapper", new[] { g }, new[] { r });
    }
    SX gs = SX.sym(symbol, 2);
    SX rs = c.call(withG(inputs, symbol, gs))["x"];
    return names ? ca.Function("wrapper", new[] { gs }, new[] { rs }, new[] { symbol }, new[] { "res" }) : ca.Function("wrapper", new[] { gs }, new[] { rs });
  }

  static int countRefCounters(string file) {
    var p = new Regex(@"static int casadi_(\w+)_ref_counter");
    return read(file).Split('\n').Count(line => p.IsMatch(line));
  }

  static readonly string[] extraWarn = { "-Wno-endif-labels", "-Wno-unused-variable", "-Wno-strict-prototypes" };

  public void test_memful_main() {
    Check.Assume(hasConic("osqp"), "no conic osqp");
    Function c = ca.conic("conic", "osqp", qpStructure(), opts("print_problem", true, "osqp.verbose", false));
    var inputs = qpInputs();
    checkCodegen(c, c.convert_in(inputs), new CG().Main().Std("c99").Flag(extraWarn).Lib("osqp"));
    c.generate("f.c");
    assertEquals(1, countRefCounters("f.c"));
    foreach (string X in new[] { "MX", "SX" }) {
      Function f = wrapper(c, inputs, X, "g", false);
      message(f.ToString());
      checkCodegen(f, dms(ca.vertcat(1, 2)), new CG().Main().Std("c99").Flag(extraWarn).Lib("osqp"));
      f.generate("f.c");
      assertEquals(1, countRefCounters("f.c"));
    }
  }

  public void test_codegen_thread_safe() {
    Check.Assume(RunSlow, "slow");
    foreach (int variant in new[] { 0, 1 }) {
      Function c;
      var inputs = new Dictionary<string, DM>();
      List<string> extralibs;
      string symbol;
      if (variant == 0 && hasConic("osqp")) {
        c = ca.conic("conic", "osqp", qpStructure(), opts("osqp.verbose", false, "osqp.alpha", 1, "osqp.eps_abs", 1e-8, "osqp.eps_rel", 1e-8, "print_time", false));
        inputs = qpInputs();
        extralibs = new List<string> { "osqp" };
        symbol = "g";
      } else if (hasNlpsol("fatrop")) {
        MX x = MX.sym("x", 2);
        c = ca.nlpsol("solver", "fatrop", new Dictionary<string, MX> { ["x"] = x, ["f"] = ca.sumsqr(x - 10), ["g"] = x[0] + x[1] }, opts("fatrop.print_level", 0, "print_time", false));
        inputs["x0"] = ca.vertcat(1, 2);
        inputs["lbx"] = ca.vertcat(-5, -4);
        inputs["ubx"] = ca.vertcat(5, 4);
        inputs["lbg"] = -2;
        inputs["ubg"] = 2;
        extralibs = new List<string> { "fatrop", "blasfeo", "m" };
        symbol = "ubx";
      } else {
        continue;
      }
      var functions = new List<(Function, Dictionary<string, DM>)> { (c, inputs) };
      foreach (string X in new[] { "MX", "SX" }) functions.Add((wrapper(c, inputs, X, symbol, true), new Dictionary<string, DM> { [symbol] = inputs[symbol] }));
      foreach (var (f, local) in functions) {
        foreach (string mapType in new[] { "serial", "thread", "openmp" }) {
          Function F = f.map(2, mapType);
          var wide = local.ToDictionary(e => e.Key, e => ca.horzcat(e.Value, e.Value * 1.2));
          CG bas = new CG().Main().Std("c99").Flag(extraWarn).NoExternal().Digits(12);
          foreach (string l in extralibs) bas.Lib(l);
          if (mapType == "openmp") {
            CG o = bas.Flag("-fopenmp").Opt("thread_safe", true).Def("CASADI_MAX_NUM_THREADS=2").Def("CASADI_THREAD_TYPE=CASADI_THREAD_TYPE_OMP");
            checkCodegen(F, F.convert_in(wide), o);
          } else if (mapType == "thread") {
            foreach (string staticInit in new[] { "0", "1" }) {
              CG o1 = new CG().Main().Std("c99").Flag(extraWarn).Flag("-pthread").NoExternal().Digits(12).Opt("thread_safe", true).Def("CASADI_MAX_NUM_THREADS=2")
                  .Def("CASADI_THREAD_TYPE=CASADI_THREAD_TYPE_POSIX").Def("CASADI_MUTEX_USE_STATIC_INIT=" + staticInit);
              foreach (string l in extralibs) o1.Lib(l);
              checkCodegen(F, F.convert_in(wide), o1);
            }
            CG o = new CG().Main().Std("c11").Flag(extraWarn).NoExternal().Digits(12).Opt("thread_safe", true).Def("CASADI_MAX_NUM_THREADS=2").Def("CASADI_THREAD_TYPE=CASADI_THREAD_TYPE_C11");
            foreach (string l in extralibs) o.Lib(l);
            checkCodegen(F, F.convert_in(wide), o);
          } else {
            checkCodegen(f, f.convert_in(local), bas);
            var scaled = local.ToDictionary(e => e.Key, e => e.Value * 1.2);
            checkCodegen(f, f.convert_in(scaled), bas);
            checkCodegen(F, f.convert_in(wide), bas);
          }
        }
      }
    }
  }

  public void test_memful_external() {
    Check.Assume(hasConic("osqp"), "no conic osqp");
    Check.Assume(RunSlow, "slow");
    Function c = ca.conic("conic", "osqp", qpStructure(), opts("print_problem", true, "osqp.verbose", false));
    var inputs = qpInputs();
    var cg = checkCodegen(c, c.convert_in(inputs), new CG().Std("c99").Flag(extraWarn).Lib("osqp").Debug());
    Function F = ((Function)cg["F"]).wrap();
    message("memful_external");
    checkCodegen(F, c.convert_in(inputs), new CG().Main().Lib((string)cg["libname"]).Debug());
    foreach (string X in new[] { "MX", "SX" }) {
      Function f = wrapper(c, inputs, X, "g", false);
      message(f.ToString());
      cg = checkCodegen(f, dms(ca.vertcat(1, 2)), new CG().Std("c99").Flag(extraWarn).Lib("osqp").Debug());
      F = ((Function)cg["F"]).wrap();
      message("memful_external");
      checkCodegen(F, dms(ca.vertcat(1, 2)), new CG().Main().Lib((string)cg["libname"]).Debug());
    }
  }

  public void test_memful_jit() {
    Check.Assume(hasConic("osqp"), "no conic osqp");
    Check.Assume(RunSlow, "slow");
    Function c = ca.conic("conic", "osqp", qpStructure(), opts("print_problem", true, "osqp.verbose", false));
    var inputs = qpInputs();
    var extra = new List<string>(extraWarn);
    foreach (string X in new[] { "MX", "SX" }) {
      string libdir = GlobalOptions.getCasadiPath(), includedir = includeDir();
      List<string> compilerFlags = extra;
      var linkerFlags = new List<string>();
      compilerFlags.Add("-I" + includedir);
      compilerFlags.Add("-L" + libdir);
      compilerFlags.Add("-losqp");
      linkerFlags.Add("-I" + includedir);
      linkerFlags.Add("-L" + libdir);
      linkerFlags.Add("-losqp");
      var jitOpts = opts("jit", true, "compiler", "shell", "jit_options", opts("compiler_flags", compilerFlags, "linker_flags", linkerFlags, "verbose", true));
      Function f = wrapper(c, inputs, X, "g", false);
      Function fj;
      if (X == "MX") {
        MX g = MX.sym("g", 2);
        fj = ca.Function("wrapper", new[] { g }, new[] { c.call(withG(inputs, "g", g))["x"] }, jitOpts);
      } else {
        SX g = SX.sym("g", 2);
        fj = ca.Function("wrapper", new[] { g }, new[] { c.call(withG(inputs, "g", g))["x"] }, jitOpts);
      }
      checkfunctionLight(f, fj, dms(ca.vertcat(1, 2)));
      CG o = new CG().Std("c99").ExternalOpts(jitOpts).Lib("osqp").Debug();
      o.extraOptions.AddRange(extra);
      checkCodegen(f, dms(ca.vertcat(1, 2)), o);
    }
  }

  public void test_cse() {
    foreach (string X in XS) {
      bool wZero, resZero;
      Function f1, f2;
      if (X == "SX") {
        SX x = SX.sym("x", 2), y = SX.sym("y", 2), p = SX.sym("p", 2);
        SX z = p * ca.cos(x + y * 5);
        SX z2 = p * ca.cos(x + y * 5);
        SX w = z - z2;
        wZero = w.is_zero();
        resZero = ca.cse(w).is_zero();
        f1 = ca.Function("f", new[] { x, y, p }, new[] { w }, opts("cse", false));
        f2 = ca.Function("f", new[] { x, y, p }, new[] { w }, opts("cse", true));
      } else {
        MX x = MX.sym("x", 2), y = MX.sym("y", 2), p = MX.sym("p", 2);
        MX z = p * ca.cos(x + y * 5);
        MX z2 = p * ca.cos(x + y * 5);
        MX w = z - z2;
        wZero = w.is_zero();
        resZero = ca.cse(w).is_zero();
        f1 = ca.Function("f", new[] { x, y, p }, new[] { w }, opts("cse", false));
        f2 = ca.Function("f", new[] { x, y, p }, new[] { w }, opts("cse", true));
      }
      assertFalse(wZero);
      assertTrue(resZero);
      assertTrue(f1.n_instructions() > 3);
      assertTrue(f2.n_instructions() <= 3);
    }
  }

  public void test_cse_call() {
    foreach (string X in XS) {
      string y;
      var ni = opts("never_inline", true);
      if (X == "SX") {
        SX x = SX.sym("x");
        var f = ca.Function("f", new[] { x }, new[] { ca.pow(x, 2) }, new[] { "x" }, new[] { "y" }, ni);
        var fcopy = ca.Function("f", new[] { x }, new[] { ca.pow(x, 2) }, new[] { "x" }, new[] { "y" }, ni);
        y = ca.cse(ca.vertcat(f.Call(ca.sin(x)), fcopy.Call(ca.sin(x)), f.Call(ca.sin(x)))).ToString();
      } else {
        MX x = MX.sym("x");
        var f = ca.Function("f", new[] { x }, new[] { ca.pow(x, 2) }, new[] { "x" }, new[] { "y" }, ni);
        var fcopy = ca.Function("f", new[] { x }, new[] { ca.pow(x, 2) }, new[] { "x" }, new[] { "y" }, ni);
        y = ca.cse(ca.vertcat(f.Call(ca.sin(x)), fcopy.Call(ca.sin(x)), f.Call(ca.sin(x)))).ToString();
      }
      message(y);
      assertTrue(y.Contains(X == "MX" ? "vertcat(@1, @1, @1)" : "[@1, @1, @1]"), y);
    }
  }

  static List<string> prefixed(string prefix, IEnumerable<string> names) => names.Select(n => prefix + n).ToList();

  static void checkZeroPattern(Function f, Function fRef, List<DM> inputs) {
    List<DM> res = f.call(inputs, false, false), resRef = fRef.call(inputs, false, false);
    for (int k = 0; k < res.Count; k++) {
      if (resRef[k].is_zero()) assertTrue(res[k].is_zero());
      else assertFalse(res[k].is_zero());
    }
  }

  /// <summary>Factory labels for one more order of derivatives, as in test_stop_diff</summary>
  static (List<string>, List<string>) derivLabels(Function f, bool withGrad) {
    var inLabels = new List<string>(f.name_in());
    inLabels.AddRange(prefixed("fwd:", f.name_in()));
    inLabels.AddRange(prefixed("adj:", f.name_out()));
    var outLabels = prefixed("fwd:", f.name_out());
    outLabels.AddRange(prefixed("adj:", f.name_in()));
    foreach (string eIn in f.name_in())
      foreach (string eOut in f.name_out()) {
        outLabels.Add("jac:" + eOut + ":" + eIn);
        if (withGrad) outLabels.Add("grad:" + eOut + ":" + eIn);
      }
    return (inLabels, outLabels);
  }

  public void test_stop_diff() {
    MX x = MX.sym("x");
    MX y = ca.sin(x);
    int maxOrder = 2;
    for (int order = 1; order <= maxOrder; order++) {
      MX testZ = ca.pow(x, order - 1) * 3;
      var fa = ca.Function("f", new[] { x }, new[] { ca.stop_diff(testZ, order) }, new[] { "x" }, new[] { "z" });
      var fb = ca.Function("f", new[] { x }, new[] { testZ }, new[] { "x" }, new[] { "z" });
      checkfunction(fa, fb, dms(DM.rand(1, 1)));
      MX z = ca.stop_diff(y, order);
      DM.rng(1);
      DM v0 = DM.rand(1);
      DM expr0 = DM.rand(1);
      var ops = new Func<MX, MX, MX>[] { (e, v) => ca.gradient(e, v), (e, v) => ca.jacobian(e, v), (e, v) => ca.jtimes(e, v, v0), (e, v) => ca.jtimes(e, v, expr0, true) };
      foreach (var op in ops) {
        MX yd = y, zd = z;
        for (int i = 0; i < maxOrder + 2; i++) {
          var fi = ca.Function("f", new[] { x }, new[] { yd, zd });
          DM x0 = DM.rand(fi.sparsity_in(0));
          List<DM> res = fi.call(dms(x0));
          if (i >= order) assertTrue(res[1].is_zero());
          else assertEquals(res[0].to_double(), res[1].to_double());
          if (i < order) {
            yd = op(yd, x);
            zd = op(zd, x);
          }
        }
      }
      var f = ca.Function("f", new[] { x }, new[] { z }, new[] { "x" }, new[] { "z" });
      var fNormal = ca.Function("f", new[] { x }, new[] { y }, new[] { "x" }, new[] { "z" });
      var fRef = ca.Function("f", new[] { x }, new[] { testZ }, new[] { "x" }, new[] { "z" });
      for (int i = 0; i < 3; i++) {
        List<DM> inputs = randInputs(f);
        if (i >= order) checkZeroPattern(f, fRef, inputs);
        else checkfunctionLight(f, fNormal, inputs);
        if (i < order) {
          var (inLabels, outLabels) = derivLabels(f, true);
          f = f.factory("f", inLabels, outLabels);
          fNormal = fNormal.factory("f", inLabels, outLabels);
          fRef = fRef.factory("f", inLabels, outLabels);
        }
      }
    }
  }

  public void test_stop_diff_cross() {
    MX x = MX.sym("x", 2);
    MX y = ca.sin(x[0]) + ca.exp(x[1]) + ca.cos(ca.pow(x[0], 2) * ca.pow(x[1], 2));
    int maxOrder = 2;
    for (int order = 1; order <= maxOrder; order++) {
      MX testZ = ca.pow(x, order - 1) * 3;
      var fa = ca.Function("f", new[] { x }, new[] { ca.stop_diff(testZ, order) }, new[] { "x" }, new[] { "z" });
      var fb = ca.Function("f", new[] { x }, new[] { testZ }, new[] { "x" }, new[] { "z" });
      checkfunction(fa, fb, dms(DM.rand(1, 1)));
      MX z = ca.stop_diff(y, order);
      DM.rng(1);
      var ops = new Func<MX, MX, MX>[] { (e, v) => ca.jacobian(e, v), (e, v) => ca.jtimes(e, v, DM.rand(v.sparsity())),
          (e, v) => ca.jtimes(e, v, DM.rand(e.sparsity()), true) };
      foreach (var op in ops) {
        MX yd = y, zd = z;
        for (int i = 0; i < maxOrder + 1; i++) {
          var fi = ca.Function("f", new[] { x }, new[] { yd, zd });
          DM x0 = DM.rand(fi.sparsity_in(0));
          List<DM> res = fi.call(dms(x0));
          if (i >= order) assertTrue(res[1].is_zero());
          else checkarray(res[0], res[1]);
          if (i < order) {
            DM.rng(5 + i);
            yd = op(yd, x);
            DM.rng(5 + i);
            zd = op(zd, x);
          }
        }
      }
      var f = ca.Function("f", new[] { x }, new[] { z }, new[] { "x" }, new[] { "z" });
      var fNormal = ca.Function("f", new[] { x }, new[] { y }, new[] { "x" }, new[] { "z" });
      var fRef = ca.Function("f", new[] { x }, new[] { testZ }, new[] { "x" }, new[] { "z" });
      if (order <= 1) {
        for (int i = 0; i < 3; i++) {
          List<DM> inputs = randInputs(f);
          if (i >= order) checkZeroPattern(f, fRef, inputs);
          else checkfunctionLight(f, fNormal, inputs);
          if (i < order) {
            var (inLabels, outLabels) = derivLabels(f, false);
            f = f.factory("f", inLabels, outLabels);
            fNormal = fNormal.factory("f", inLabels, outLabels);
            fRef = fRef.factory("f", inLabels, outLabels);
          }
        }
      }
    }
  }

  public void test_stop_diff_p() {
    MX x = MX.sym("x");
    MX y = MX.sym("y");
    MX yy = ca.sin(y * x);
    int maxOrder = 2;
    for (int order = 1; order <= maxOrder; order++) {
      DM.rng(1);
      message("order " + order);
      MX testZ = ca.sin(y) * ca.pow(x, order - 1);
      MX z = ca.stop_diff(yy, x, order);
      DM v0 = DM.rand(1);
      DM expr0 = DM.rand(1);
      var ops = new Func<MX, MX, MX>[] { (e, v) => ca.gradient(e, v), (e, v) => ca.jacobian(e, v), (e, v) => ca.jtimes(e, v, v0), (e, v) => ca.jtimes(e, v, expr0, true) };
      foreach (var op in ops) {
        MX yd = yy, zd = z;
        for (int i = 0; i < maxOrder + 2; i++) {
          var fi = ca.Function("f", new[] { x, y }, new[] { yd, zd });
          DM x0 = DM.rand(fi.sparsity_in(0));
          DM p0 = DM.rand(fi.sparsity_in(1));
          List<DM> res = fi.call(dms(x0, p0));
          if (i >= order) assertTrue(res[1].is_zero());
          else assertEquals(res[0].to_double(), res[1].to_double());
          if (i < order) {
            yd = op(yd, x);
            zd = op(zd, x);
          }
        }
      }
      var f = ca.Function("f", new[] { x, y }, new[] { z }, new[] { "x", "y" }, new[] { "z" });
      var fNormal = ca.Function("f", new[] { x, y }, new[] { yy }, new[] { "x", "y" }, new[] { "z" });
      var fRef = ca.Function("f", new[] { x, y }, new[] { testZ }, new[] { "x", "y" }, new[] { "z" });
      for (int i = 0; i < 3; i++) {
        List<DM> inputs = randInputs(f);
        if (i >= order) {
          fRef.disp(true);
          f.generate("f.c");
          checkZeroPattern(f, fRef, inputs);
        } else {
          checkfunctionLight(f, fNormal, inputs);
        }
        if (i < order) {
          var inLabels = new List<string>(f.name_in());
          inLabels.AddRange(prefixed("adj:", f.name_out()));
          var outLabels = new List<string>(f.name_out());
          outLabels.AddRange(prefixed("adj:", f.name_in()));
          f = f.factory("f", inLabels, outLabels);
          fNormal = fNormal.factory("f", inLabels, outLabels);
          fRef = fRef.factory("f", inLabels, outLabels);
        }
      }
    }
  }

  public void test_no_hess2() {
    Check.Assume(hasConic("qpoases"), "no conic qpoases");
    MX y = MX.sym("y", 2);
    var f = ca.Function("foo", new[] { y }, new[] { ca.vertcat(ca.sin(y[0] * y[1]), ca.cos(y[0] * y[1])) }, opts("never_inline", true, "print_out", true));
    var opti = ca.Opti();
    MX x = opti.variable(2);
    opti.subject_to(ca.ge(ca.no_hess(f.Call(ca.pow(x, 2))), 0.0));
    opti.minimize(ca.sumsqr(x - 2));
    opti.solver("sqpmethod");
    Function F = opti.to_function("F", new[] { x }, new[] { x });
    F.generate("F.c");
    string code = read("F.c");
    foreach (Match m in Regex.Matches(code, @"\w+FS")) assertTrue(m.Value.Split('_').Length <= 2, m.Value);
  }

  public void test_issue_1522() {
    MX V2 = MX.sym("X", 8);
    List<MX> xt = ca.vertsplit(V2, new long[] { 0, 6, 8 });
    List<MX> xs = ca.vertsplit(xt[0], new long[] { 0, 3, 6 });
    List<MX> travels = ca.vertsplit(xt[1], new long[] { 0, 1, 2 });
    MX dist = 0.0;
    for (int j = 0; j < 2; j++) dist = dist + ca.sum1(ca.pow(xs[0] - (xs[j] + travels[j]), 2));
    var p = ca.MX(0, 1);
    var f = ca.Function("f", new[] { V2, p }, new[] { -dist });
    foreach (bool compact in new[] { true, false })
      foreach (Function F in new[] { f, f.expand() }) {
        Sparsity HF = F.reverse(1).jac_sparsity(0, 0, compact);
        assertTrue(HF.is_symmetric());
      }
  }

  public void test_issue_3074() {
    MX x = MX.sym("x", 4);
    var f = ca.Function("f", new[] { x }, new MX[] { ca.dot(1.0, 1.0) });
    message(f.ToString());
    Function fr = f.reverse(1);
    message(fr.ToString());
    x = MX.sym("x", ca.sparsify(ca.vertcat(0, 1, 1, 0)).sparsity());
    MX y = ca.dot(x, x);
    f = ca.Function("f", new[] { x }, new[] { ca.dot(y, y) });
    message(f.ToString());
    fr = f.reverse(1);
    message(fr.ToString());
    message(string.Join(", ", fr.jac_sparsity()));
    f = ca.StringDeserializer(ISSUE_3074_DATA).unpack_function();
    message(f.reverse(1).ToString());
    f.reverse(1).jac_sparsity();
  }

  const string ISSUE_3074_DATA = "jhpnnagiieahaaaadaaaaaaaaaaaaaaaaafaegaakaaaaaaaneifgefhogdgehjgpgogbaaaaaaapaaaaaaaegfgmgbgjhpfbgchhgfhngfgogehdhaaaaaacaaaaaaahaaaaaaaaaaaaaaabababababababaaaaaaaaaaaaaaaaahaaaaaaaaaaaaaaaegfaaaaaaaaaaaaaaabaaaaaaaaaaaaaaabaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaabaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaeggaaaaaaaaaaaaaaacaaaaaaaaaaaaaaabaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaacaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaabaaaaaaaaaaaaaaachbaaaaaaaaaaaaaaaegkaaaaaaaaaaaaaaagaaaaaaaaaaaaaaabaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaagaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaabaaaaaaaaaaaaaaacaaaaaaaaaaaaaaadaaaaaaaaaaaaaaaeaaaaaaaaaaaaaaafaaaaaaaaaaaaaaaegeaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaabaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaachdaaaaaaaaaaaaaaaegiaaaaaaaaaaaaaaaeaaaaaaaaaaaaaaabaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaeaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaabaaaaaaaaaaaaaaacaaaaaaaaaaaaaaadaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaahaaaaaaaaaaaaaaacaaaaaaajgadcaaaaaaajgbdcaaaaaaajgcdcaaaaaaajgddcaaaaaaajgedcaaaaaaajgfdcaaaaaaajggdaaaaaaaaaaaaaaaaaabagaaaaaaadhpgfhchdgfgbahaaaaaaakgjgehpfehngahaaaaaaaaaaaaaaaafaaaaaaadgmgbgoghgaaegbaaaaaaaaaaaaaaaaebababaaabababaaapbfilobfilobfnpdmfpicmfpicmfpnpdaaaaaeaaaaaaaaaaaaaaaabakdmiadcooijhfeodaaaaaaaaaaaaaaaabaaaaaaaocdaaaaaaangehihaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaachfaaaaaaaaaaaaaaahaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaabaaaaaaahaaaaaaaaaaaaaaaegpcaaaaaaaaaaaaaaaaaaaaaachaaaaaaaaaaaaaaaaeaaaaaaaehjgngfgegndaaaaaacaaaaaaaaaaaaaaaegpcaaaaaaaaaaaaaaaaaaaaaachaaaaaaaaaaaaaaaagaaaaaaacgpgegjhocihegpcaaaaaaaaaaaaaaaaaaaaaachaaaaaaaaaaaaaaaaiaaaaaaacgpgegjhocghpfihchbaaaaaaaaaaaaaaaegndaaaaaacaaaaaaaaaaaaaaaegpcaaaaaaaaaaaaaaaaaaaaaachaaaaaaaaaaaaaaaalaaaaaaaegfgchiccgpgegjhocihjcegpcaaaaaaaaaaaaaaaaaaaaaachaaaaaaaaaaaaaaaanaaaaaaaegfgchiccgpgegjhocghpfihjcchbaaaaaaaaaaaaaaaegndaaaaaagaaaaaaaaaaaaaaaegpcaaaaaaaaaaaaaaaaaaaaaachaaaaaaaaaaaaaaaaiaaaaaaacgpgegjhocggpfihegpcaaaaaaaaaaaaaaaaaaaaaachaaaaaaaaaaaaaaaaiaaaaaaacgpgegjhocbgpfihegpcaaaaaaaaaaaaaaaaaaaaaachaaaaaaaaaaaaaaaajaaaaaaabgdgdgfgmgocbgpfihegpcaaaaaaaaaaaaaaaaaaaaaachaaaaaaaaaaaaaaaakaaaaaaabgdgdgfgmgocngbgpfihegpcaaaaaaaaaaaaaaaaaaaaaachaaaaaaaaaaaaaaaalaaaaaaabgdgdgfgmgoccgpfihocfhegpcaaaaaaaaaaaaaaaaaaaaaachaaaaaaaaaaaaaaaalaaaaaaabgdgdgfgmgoccgpfihocjhchcaaaaaaaaaaaaaaaegmcaaaaaaadaaaaaaaaaaaaaaaachdaaaaaaaaaaaaaaaegmcaaaaaaadaaaaaaaaaaaaaaaachdaaaaaaaaaaaaaaaegndaaaaaaeaaaaaaaaaaaaaaaegpcaaaaaaaaaaaaaaaaaaaaaachaaaaaaaaaaaaaaaagaaaaaaacgpgegjhochgegpcaaaaaaaaaaaaaaaaaaaaaachaaaaaaaaaaaaaaaagaaaaaaacgpgegjhocdgegpcaaaaaaaaaaaaaaaaaaaaaachaaaaaaaaaaaaaaaagaaaaaaacgpgegjhocngegpcaaaaaaaaaaaaaaaaaaaaaachaaaaaaaaaaaaaaaalaaaaaaabgdgdgfgmgoccgpfihoccgcheaaaaaaaaaaaaaaabaaaaaaaaaaaaaaaaaaaaaaabaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaahaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaabaaaaaaaaaaaaaaaaa";

  public void test_issue_4345() {
    foreach (string X in XS) {
      var fs = new List<Function>();
      if (X == "SX") {
        SX x = SX.sym("x", Sparsity.upper(3)), y = SX.sym("y", Sparsity.upper(3));
        foreach (SX z in new[] { ca.pow(x, 2), ca.pow(x, 2) + y[0, 0], x * y }) fs.Add(ca.Function("f", new[] { x, y }, new[] { z }));
      } else {
        MX x = MX.sym("x", Sparsity.upper(3)), y = MX.sym("y", Sparsity.upper(3));
        foreach (MX z in new[] { ca.pow(x, 2), ca.pow(x, 2) + y[0, 0], x * y }) fs.Add(ca.Function("f", new[] { x, y }, new[] { z }));
      }
      long[] adjYnnz = { 0, 6, 6 };
      for (int k = 0; k < 3; k++) {
        Function adj = fs[k].reverse(1);
        assertEquals(6L, adj.nnz_out(0));
        assertEquals(adjYnnz[k], adj.nnz_out(1));
      }
    }
  }

  public void test_issue_3134() {
    Check.Assume(hasNlpsol("ipopt"), "no nlpsol ipopt");
    MX p = MX.sym("p");
    MX v = MX.sym("v");
    MX u = MX.sym("u");
    MX x0 = ca.vertcat(p, v);
    MX rhs = ca.vertcat(v, u);
    MX t0 = MX.sym("t0");
    MX DT = MX.sym("DT");
    var F = ca.Function("F", new[] { x0, u, t0, DT }, new[] { x0 + DT * rhs }, new[] { "x0", "u", "t0", "DT" }, new[] { "xf" });
    var x = new List<MX>();
    var g = new List<MX>();
    MX Xk = MX.sym("Xk", 2);
    x.Add(Xk);
    MX T = MX.sym("T");
    x.Add(T);
    g.Add(T);
    MX t0b = MX.sym("t0");
    x.Add(t0b);
    var DTs = new[] { T / 2, t0b + T - (t0b + T / 2) };
    for (int i = 0; i < 2; i++) {
      MX Uk = MX.sym("Uk");
      x.Add(Uk);
      MX XkNext = MX.sym("Xk_next", 2);
      x.Add(XkNext);
      MX Xf = F.call(new Dictionary<string, MX> { ["x0"] = Xk, ["u"] = Uk, ["t0"] = t0b, ["DT"] = DTs[i] })["xf"];
      g.Add(XkNext - Xf);
      Xk = XkNext;
    }
    var nlp = new Dictionary<string, MX> { ["x"] = ca.veccat(x), ["f"] = T, ["g"] = ca.veccat(g) };
    ca.nlpsol("solver", "ipopt", nlp);
  }

  static void assertBools(bool[] expected, List<bool> actual) =>
    assertTrue(expected.SequenceEqual(actual), "expected [" + string.Join(",", expected) + "], got [" + string.Join(",", actual) + "]");

  static void whichDependsChecks(Function f, bool exempt) {
    string[] z = { "z" }, z2 = { "z2" }, z3 = { "z3" }, z4 = { "z4" }, z5 = { "z5" }, z7 = { "z7" }, z8 = { "z8" };
    const bool T = true, F = false;
    assertBools(new[] { T, T }, f.which_depends("x", z));
    assertBools(new[] { F, F }, f.which_depends("x", z, 2));
    if (!exempt) assertBools(new[] { F, F }, f.which_depends("x", z2, 2));
    assertBools(new[] { T, T }, f.which_depends("x", z3, 2));
    assertBools(new[] { T, T }, f.which_depends("x", z4));
    assertBools(new[] { T, T }, f.which_depends("x", z4, 2));
    assertBools(new[] { F, F }, f.which_depends("x", z5));
    assertBools(new[] { F, F }, f.which_depends("x", z5, 2));
    assertBools(new[] { F, F }, f.which_depends("x", new[] { "z5", "z6" }));
    if (!exempt) {
      assertBools(new[] { T, T, T, T, T, T, T, F, F, F }, f.which_depends("x", new[] { "z", "z2", "z3", "z4", "z5", "z6" }, 1, true));
      assertBools(new[] { F, F, F, F, T, T, T, F, F, F }, f.which_depends("x", new[] { "z", "z2", "z3", "z4", "z5", "z6" }, 2, true));
    }
    assertBools(new[] { T, T, T, T, T, F, F, F }, f.which_depends("x", new[] { "z", "z3", "z4", "z5", "z6" }, 1, true));
    assertBools(new[] { F, F, T, T, T, F, F, F }, f.which_depends("x", new[] { "z", "z3", "z4", "z5", "z6" }, 2, true));
    assertBools(new[] { T, F }, f.which_depends("x", z7));
    assertBools(new[] { F, T }, f.which_depends("y", z7));
    assertBools(new[] { F, F }, f.which_depends("x", z7, 2));
    assertBools(new[] { F, F }, f.which_depends("y", z7, 2));
    assertBools(new[] { T, F }, f.which_depends("x", z8, 2));
    assertBools(new[] { F, T }, f.which_depends("y", z8, 2));
    assertBools(new[] { T, F, T, F }, f.which_depends("x", new[] { "z7", "z8" }, 1, true));
    assertBools(new[] { F, F, T, F }, f.which_depends("x", new[] { "z7", "z8" }, 2, true));
    assertBools(new[] { F, T, F, T }, f.which_depends("y", new[] { "z7", "z8" }, 1, true));
    assertBools(new[] { F, F, F, T }, f.which_depends("y", new[] { "z7", "z8" }, 2, true));
  }

  public void test_which_depends() {
    SX x = SX.sym("x", 2);
    SX y = SX.sym("y", 2);
    SX z = x + y;
    SX z2 = x * y;
    SX z3 = ca.sin(ca.dot(x, y));
    SX z4 = ca.cos(x);
    SX z5 = ca.sin(y);
    SX z6 = 12.0;
    SX z7 = ca.vertcat(x[0], y[1]);
    SX z8 = ca.vertcat(ca.pow(x[0], 2), ca.pow(y[1], 2));
    var f = ca.Function("f", new[] { x, y }, new[] { z, z2, z3, z4, z5, z6, z7, z8 }, new[] { "x", "y" }, new[] { "z", "z2", "z3", "z4", "z5", "z6", "z7", "z8" });
    whichDependsChecks(f, false);
    if (!RunSlow) return;
    var cg = checkCodegen(f, randInputs(f), new CG().JacSparsity().Forward());
    whichDependsChecks((Function)cg["F"], true);
  }

  public void test_cat_input() {
    MX x = MX.sym("x", 2);
    var fmx = ca.Function("f", new[] { ca.vertcat(ca.MX(2, 1), x) }, new[] { ca.sin(x) }, opts("always_inline", true));
    message(fmx.ToString());
    MX X = MX.sym("X", 4);
    message(fmx.Call(X).ToString());
    var fEval = ca.Function("f_eval", new[] { X }, new[] { fmx.Call(X) });
    message(fEval.Call(ca.vertcat(double.NaN, double.NaN, 3, double.NaN)).ToString());
  }

  public void test_reshape_input() {
    MX x = MX.sym("x", 6, 2);
    SX xsx = SX.sym("x", 6, 2);
    var fmx = ca.Function("f", new[] { ca.reshape(x, 3, 4) }, new[] { ca.sin(x) }, opts("always_inline", true));
    var fsx = ca.Function("f", new[] { ca.reshape(xsx, 3, 4) }, new[] { ca.sin(xsx) });
    DM.rng(1);
    foreach (DM input in new[] { DM.rand(3, 4) }) checkfunctionLight(fmx, fsx, dms(input));
  }

  public void test_sparsity_cast_input() {
    MX X0 = MX.sym("X", Sparsity.lower(3));
    ca.Function("f", new[] { X0 }, new[] { ca.pow(X0, 2) }, opts("always_inline", true));
    Sparsity sp = Sparsity.lower(3);
    long N = sp.nnz();
    MX x = MX.sym("x", N);
    SX xsx = SX.sym("x", N);
    var fmx = ca.Function("f", new[] { ca.sparsity_cast(x, sp) }, new[] { ca.sin(x) }, opts("always_inline", true));
    var fsx = ca.Function("f", new[] { ca.sparsity_cast(xsx, sp) }, new[] { ca.sin(xsx) });
    foreach (DM input in new[] { sparseDM(new double[,] { { 1, 0, 0 }, { 3, 7, 0 }, { 1, 8, 9 } }), sparseDM(new double[,] { { 1, 0, 0 }, { 3, 0, 0 }, { 0, 8, 9 } }),
        ca.DM(new double[,] { { 1, 1.1, 8 }, { 3, 1.2, 3 }, { 1, 8, 9 } }) })
      checkfunctionLight(fmx, fsx, dms(input));
    foreach (MX X in new[] { MX.sym("X", Sparsity.lower(3)), MX.sym("X", Sparsity.diag(3)), MX.sym("X", 3, 3) }) {
      var fcmx = ca.Function("f", new[] { X }, new[] { fmx.Call(X) });
      var fcsx = ca.Function("f", new[] { X }, new[] { fsx.Call(X) });
      checkfunctionLight(fcmx, fcsx, dms(DM.rand(X.sparsity())));
    }
  }

  public void test_external() {
    Check.Assume(RunSlow, "slow");
    string src = "../python/assets/externalfun1.c";
    Check.Assume(File.Exists(src), "run from test/dotnet: " + src + " not found");
    assertInException("config failed", () => compileExternal("F", src, opts(), new CG()));
    foreach (int nArgs in new[] { 3, 7 }) {
      var cfg = Enumerable.Repeat("foo", nArgs).ToList();
      Function externalfun1 = compileExternal("F", src, opts("config_args", cfg), new CG().Debug());
      string libname = externalLibname;
      message(libname);
      List<DM> r = externalfun1.call(ins(0));
      message(jstr(externalfun1.stats()));
      assertEquals(nArgs + 1, (int)r[0].to_double());
      assertEquals(0, (int)r[1].to_double());
      MX x = MX.sym("x");
      var g = ca.Function("g", new[] { x }, new[] { externalfun1.Call(ca.sin(x)) + externalfun1.Call(ca.cos(x)) });
      assertEquals(2.0 * (nArgs + 1), g.Call(0.3).to_double());
      checkCodegen(g, ins(0.3), new CG().Lib(libname));
      checkSerialize(externalfun1, ins(0.3));
      checkSerialize(g, ins(0.3));
      externalfun1.Dispose();
      g.Dispose();
      GC.Collect();
    }
  }

  public void test_cache() {
    MX x = MX.sym("x");
    var f = ca.Function("f", new[] { x }, new[] { ca.pow(x, 2) });
    Function ff = f.forward(1);
    Function ff2 = f.forward(1);
    long h = ff.hash();
    assertEquals(ff.hash(), ff2.hash());
    ff.Dispose();
    ff2.Dispose();
    GC.Collect();
    Function f2a = f.forward(2);
    Function f2b = f.forward(2);
    ff = f.forward(1);
    message(h + " " + ff.hash());
    assertEquals(f2a.hash(), f2b.hash());
  }

  static string gradientCode(Function B, string fname) {
    int N = 4;
    MX z = MX.sym("z", 2);
    MX s = z;
    for (int i = 0; i < N; i++) s = B.Call(s);
    var A = ca.Function("A", new[] { z }, new[] { ca.dot(s, s) });
    Function g = A.factory("g", new[] { "i0" }, new[] { "o0", "grad:o0:i0" });
    var cg = ca.CodeGenerator(fname);
    cg.add(g);
    cg.generate();
    cg.Dispose();
    g.Dispose();
    A.Dispose();
    return read(fname);
  }

  static int count(string s, string sub) {
    int n = 0;
    for (int i = s.IndexOf(sub, StringComparison.Ordinal); i >= 0; i = s.IndexOf(sub, i + sub.Length, StringComparison.Ordinal)) n++;
    return n;
  }

  public void test_cache_stale_entry() {
    MX x = MX.sym("x", 2);
    var B = ca.Function("B", new[] { x }, new[] { ca.vertcat(ca.sin(x[0]) * x[1], ca.cos(x[1]) + x[0]) }, new[] { "x" }, new[] { "y" });
    assertEquals(1, count(gradientCode(B, "adj_cache1.c"), "/* adj1_B:"));
    GC.Collect();
    assertEquals(1, count(gradientCode(B, "adj_cache2.c"), "/* adj1_B:"));
  }

  const string COPY_ELISION_1 =
      "casadi_int i;\n" +
      "  casadi_real *rr, w0, *w1=w+1, w2, *w3=w+6, *w4=w+16, *w5=w+26, *w6=w+31, *w7=w+36, *w8=w+41, *w9=w+66;\n" +
      "  const casadi_real *cs, *wr3, *wr4, *wr6, *wr9;\n" +
      "  /* #0: @0 = 0 */\n" +
      "  w0 = 0.;\n" +
      "  /* #1: @1 = ones(1x5) */\n" +
      "  casadi_fill(w1, 5, 1.);\n" +
      "  /* #2: @2 = 6 */\n" +
      "  w2 = 6.;\n" +
      "  /* #3: @3 = input[0][0] */\n" +
      "  wr3 = arg[0] ? arg[0] : casadi_zeros;\n" +
      "  /* #4: @4 = (@2*@3) */\n" +
      "  for (i=0, rr=w4, cs=wr3; i<10; ++i) (*rr++)  = (w2*(*cs++));\n" +
      "  /* #5: @5 = @4[2:7] */\n" +
      "  for (rr=w5, cs=w4+2; cs!=w4+7; cs+=1) *rr++ = *cs;\n" +
      "  /* #6: @6 = input[1][0] */\n" +
      "  wr6 = arg[1] ? arg[1] : casadi_zeros;\n" +
      "  /* #7: @5 = (@5+@6) */\n" +
      "  for (i=0, rr=w5, cs=wr6; i<5; ++i) (*rr++) += (*cs++);\n" +
      "  /* #8: @2 = 12 */\n" +
      "  w2 = 12.;\n" +
      "  /* #9: @4 = input[2][0] */\n" +
      "  wr4 = arg[2] ? arg[2] : casadi_zeros;\n" +
      "  /* #10: @4 = (@2*@4) */\n" +
      "  for (i=0, rr=w4, cs=wr4; i<10; ++i) (*rr++)  = (w2*(*cs++));\n" +
      "  /* #11: @4 = monitor(@4, b) */\n" +
      "  CASADI_PRINTF(\"b:\\n\");\n" +
      "  casadi_print_canonical(casadi_s0, w4);\n" +
      "  CASADI_PRINTF(\"\\n\");\n" +
      "  /* #12: @7 = @4[2:7] */\n" +
      "  for (rr=w7, cs=w4+2; cs!=w4+7; cs+=1) *rr++ = *cs;\n" +
      "  /* #13: @5 = (@5+@7) */\n" +
      "  for (i=0, rr=w5, cs=w7; i<5; ++i) (*rr++) += (*cs++);\n" +
      "  /* #14: @0 = mac(@1,@5,@0) */\n" +
      "  casadi_mtimes_dense(w1, 1, 5, w5, 1, (&w0), 0);\n" +
      "  /* #15: output[0][0] = @0 */\n" +
      "  if (res[0]) res[0][0] = w0;\n" +
      "  /* #16: @0 = 9 */\n" +
      "  w0 = 9.;\n" +
      "  /* #17: @8 = vertcat(@3, @6, @4) */\n" +
      "  rr=w8;\n" +
      "  for (i=0, cs=wr3; i<10; ++i) *rr++ = *cs++;\n" +
      "  for (i=0, cs=wr6; i<5; ++i) *rr++ = *cs++;\n" +
      "  for (i=0, cs=w4; i<10; ++i) *rr++ = *cs++;\n" +
      "  /* #18: @8 = (@0*@8) */\n" +
      "  for (i=0, rr=w8, cs=w8; i<25; ++i) (*rr++)  = (w0*(*cs++));\n" +
      "  /* #19: output[1][0] = @8 */\n" +
      "  casadi_copy(w8, 25, res[1]);\n" +
      "  /* #20: output[2][0] = @3 */\n" +
      "  casadi_copy(wr3, 10, res[2]);\n" +
      "  /* #21: {NULL, @9} = vertsplit(@3) */\n" +
      "  wr9 = wr3+7;\n" +
      "  /* #22: @9 = sqrt(@9) */\n" +
      "  for (i=0, rr=w9, cs=wr9; i<3; ++i) *rr++ = sqrt( *cs++ );\n" +
      "  /* #23: output[3][0] = @9 */\n" +
      "  casadi_copy(w9, 3, res[3]);\n" +
      "  /* #24: @3 = sin(@3) */\n" +
      "  for (i=0, rr=w3, cs=wr3; i<10; ++i) *rr++ = sin( *cs++ );\n" +
      "  /* #25: {@6, @1} = vertsplit(@3) */\n" +
      "  casadi_copy(w3, 5, w6);\n" +
      "  casadi_copy(w3+5, 5, w1);\n" +
      "  /* #26: @6 = (@6+@1) */\n" +
      "  for (i=0, rr=w6, cs=w1; i<5; ++i) (*rr++) += (*cs++);\n" +
      "  /* #27: output[4][0] = @6 */\n" +
      "  casadi_copy(w6, 5, res[4]);\n" +
      "  return 0;";

  const string COPY_ELISION_2 =
      "casadi_int i;\n" +
      "  casadi_real *rr, w1, w2, w3, *w4=w+0;\n" +
      "  const casadi_real *cs, *wr0, *wr4;\n" +
      "  /* #0: @0 = input[1][0] */\n" +
      "  wr0 = arg[1] ? arg[1] : casadi_zeros;\n" +
      "  /* #1: @1 = 100 */\n" +
      "  w1 = 100.;\n" +
      "  /* #2: @2 = input[0][0] */\n" +
      "  w2 = arg[0] ? arg[0][0] : 0;\n" +
      "  /* #3: @3 = floor(@2) */\n" +
      "  w3 = floor( w2 );\n" +
      "  /* #4: @1 = (@1*@3) */\n" +
      "  w1 *= w3;\n" +
      "  /* #5: @4 = @0[(4:90;@1)] */\n" +
      "  wr4 = wr0 + 4 + (int)w1;\n" +
      "  /* #6: @4 = sin(@4) */\n" +
      "  for (i=0, rr=w4, cs=wr4; i<86; ++i) *rr++ = sin( *cs++ );\n" +
      "  /* #7: @1 = 100 */\n" +
      "  w1 = 100.;\n" +
      "  /* #8: @2 = floor(@2) */\n" +
      "  w2 = floor( w2 );\n" +
      "  /* #9: @1 = (@1*@2) */\n" +
      "  w1 *= w2;\n" +
      "  /* #10: @2 = @0[(7;@1)] */\n" +
      "  i = 7 + (int) w1;\n" +
      "  w2 = i>=0 && i<1000 ? wr0[i] : casadi_nan;\n" +
      "  /* #11: @4 = (@4+@2) */\n" +
      "  for (i=0, rr=w4; i<86; ++i) (*rr++) += w2;\n" +
      "  /* #12: output[0][0] = @4 */\n" +
      "  casadi_copy(w4, 86, res[0]);\n" +
      "  return 0;";

  public void test_copy_elision() {
    long backup = GlobalOptions.getCopyElisionMinSize();
    try {
      foreach (bool enabled in new[] { false, true }) {
        GlobalOptions.setCopyElisionMinSize(enabled ? 2 : -1);
        MX x = MX.sym("x", 10);
        MX a = MX.sym("a", 5);
        MX bs = MX.sym("b", 10);
        MX b = (bs * 12).monitor("b");
        MX y = x * 6;
        List<MX> s01 = ca.vertsplit(x, new long[] { 0, 7, 10 });
        List<MX> ss = ca.vertsplit(ca.sin(x), new long[] { 0, 5, 10 });
        MX z = y[2..7] + a + b[2..7];
        z = ca.sum1(z);
        MX z2 = ca.vertcat(x, a, b) * 9;
        var f = ca.Function("f", new[] { x, a, bs }, new[] { z, z2, x, ca.sqrt(s01[1]), ss[0] + ss[1] });
        f.disp(true);
        f.generate("f.c");
        string code = read("f.c");
        assertEquals(enabled, code.Contains("wr3"));
        if (enabled) {
          message(code);
          assertTrue(code.Contains(COPY_ELISION_1));
        }
        checkCodegen(f, randInputs(f), new CG());
        MX X = MX.sym("x", 100, 10);
        MX i = MX.sym("i");
        MX yy = X.get(false, ca.Slice(4, 90), i);
        f = ca.Function("f", new[] { i, X }, new[] { ca.sin(yy) + X.get(false, ca.Slice(7), i) });
        f.disp(true);
        f.generate("f.c");
        code = read("f.c");
        assertEquals(enabled, code.Contains("wr"));
        Match m = Regex.Match(code, @"\*sz_w = (\d+);");
        assertTrue(m.Success);
        assertEquals(enabled ? 86 : 1086, int.Parse(m.Groups[1].Value));
        if (enabled) {
          message(code);
          assertTrue(code.Contains(COPY_ELISION_2));
        }
        checkCodegen(f, ins(3, DM.rand(f.sparsity_in(1))), new CG().Std("c99"));
      }
    } finally {
      GlobalOptions.setCopyElisionMinSize(backup);
    }
  }

  static List<List<double>> testPoints(List<double[]> knotsOrig) {
    var selector = new Func<double[], double>[] { e => e[0] - 0.1, e => e[0], e => (e[0] + e[1]) / 2, e => e[e.Length - 1], e => e[e.Length - 1] + 0.1 };
    var r = new List<List<double>> { new List<double>() };
    foreach (double[] k in knotsOrig) {
      var next = new List<List<double>>();
      foreach (var prefix in r)
        foreach (var sel in selector) next.Add(new List<double>(prefix) { sel(k) });
      r = next;
    }
    return r;
  }

  static double[] knotsExpand(double[] k) {
    double[] r = new double[k.Length + 6];
    for (int i = 0; i < 3; i++) {
      r[i] = k[0];
      r[k.Length + 3 + i] = k[k.Length - 1];
    }
    Array.Copy(k, 0, r, 3, k.Length);
    return r;
  }

  static DM vcat(List<double> a) => a.ToArray();

  /// <summary>Random sample of n elements (Java's shuffle with seed 1 is not reproducible here; any fixed sample serves)</summary>
  static List<List<double>> sample(List<List<double>> all, int n) {
    var r = new List<List<double>>(all);
    var rnd = new Random(1);
    for (int i = r.Count - 1; i > 0; i--) {
      int j = rnd.Next(i + 1);
      (r[i], r[j]) = (r[j], r[i]);
    }
    return r.GetRange(0, n);
  }

  static double[] linspace(double a, double b, int n) {
    double[] r = new double[n];
    for (int i = 0; i < n; i++) r[i] = a + (b - a) * i / (n - 1);
    return r;
  }

  static long[] fill(int n, long v) => Enumerable.Repeat(v, n).ToArray();

  static Dictionary<string, object> blazingJitOptions() => opts("flags", new[] { "-I" + includeDir(), "-g", "-ffast-math", "-march=native" });

  public void test_blazing_spline() {
    Check.Assume(CasadiMeta.feature_list().Contains("simde"), "no simde");
    Check.Assume(Importer.has_plugin("shell"), "no importer shell");
    Check.Assume(RunSlow, "slow");
    var allKnots = new List<double[]> { new double[] { 0, 0.2, 0.5, 0.8, 1 }, new double[] { 0, 0.1, 0.5, 0.9, 1 }, new double[] { 0, 0.4, 1 }, new double[] { 0, 0.3, 1 },
        new double[] { 0, 0.2, 1 } };
    var allKnotsAlt = new List<double[]> { new double[] { 0, 0.3, 0.6, 0.7, 1 }, new double[] { 0, 0.2, 0.4, 0.8, 1 }, new double[] { 0, 0.5, 1 }, new double[] { 0, 0.4, 1 },
        new double[] { 0, 0.3, 1 } };
    double[] x0all = { 0.4, 0.5, 0.6, 0.7, 0.8 };
    for (int N = 1; N <= 5; N++) {
      foreach (bool precomputeCoeff in new[] { false, true }) {
        foreach (bool precomputeGrid in new[] { false, true }) {
          if (precomputeCoeff && N > 3) continue;
          foreach (bool parametricGrid in new[] { false, true }) {
            if (parametricGrid && precomputeCoeff && N > 3) continue;
            DM x0 = x0all.Take(N).ToArray();
            List<double[]> knotsOrig = allKnots.GetRange(0, N);
            double[][] knots = knotsOrig.Select(knotsExpand).ToArray();
            MX x = MX.sym("x", N);
            int nc = 1;
            foreach (double[] k in knots) nc *= k.Length - 4;
            DM.rng(1);
            DM data = DM.rand(nc);
            MX C = MX.sym("C", nc, 1);
            MX Y = ca.bspline(x, C, knots, fill(N, 3), 1);
            var Fref = ca.Function("f", new[] { x, C }, new[] { Y });
            var o = opts("precompute_coeff", precomputeCoeff, "precompute_grid", precomputeGrid, "pedantic_mode_order", "ignore", "pedantic_mode_size", "ignore",
                "jit", true, "jit_options", blazingJitOptions());
            Function F, Ftest, Fref2 = null;
            MX K = null;
            DM knotsData = null;
            if (parametricGrid) {
              long[] knotDims = knots.Select(k => (long)k.Length).ToArray();
              var stacked = knots.SelectMany(k => k).ToList();
              F = ca.blazing_spline("F", knotDims, o);
              knotsData = vcat(stacked);
              K = MX.sym("K", stacked.Count);
              List<MX> Fcall = F.call(new[] { x, C, K });
              List<MX> Fsub = ca.substitute(new[] { Fcall[0] }, new[] { K }, new[] { ca.MX(knotsData) });
              Ftest = ca.Function("F_test", new[] { x, C }, Fsub);
            } else {
              Fref2 = ca.Function("f", new[] { x, C }, new[] { Y, ca.jacobian(Y, x) });
              Fref2.call(dms(DM.zeros(N), data));
              F = ca.blazing_spline("F", knots, o);
              Ftest = ca.Function("F_test", new[] { x, C }, new[] { F.Call(x, C) });
            }
            List<List<double>> allPoints = testPoints(knotsOrig);
            if (N > 3) allPoints = sample(allPoints, 125);
            foreach (var a in allPoints) checkfunctionLight(Ftest, Fref, dms(vcat(a), data));
            if (!parametricGrid) {
              Function F2 = Function.deserialize(F.serialize(opts("debug", true)));
              foreach (var a in allPoints) checkfunctionLight(F, F2, dms(vcat(a), data));
            }
            List<double> last = allPoints[allPoints.Count - 1];
            checkCodegen(Ftest, dms(vcat(last), data), new CG().Std("c99"));
            MX ySin = parametricGrid ? ca.substitute(ca.sin(F.Call(x, C, K)), K, ca.MX(knotsData)) : ca.sin(F.Call(x, C));
            var f = ca.Function("f", new[] { x, C }, new[] { ySin, ca.jacobian(ySin, x), ca.gradient(ySin, x) });
            var fcse = ca.Function("fcse", new[] { x, C }, new[] { ySin, ca.jacobian(ySin, x), ca.gradient(ySin, x) }, opts("cse", true));
            checkfunctionLight(f, fcse, dms(x0, data));
            var (h0, h1) = ca.hessian(ySin, x);
            var outs = new List<MX> { ySin, ca.jacobian(ySin, x), ca.gradient(ySin, x), h0, h1 };
            var f2 = ca.Function("f", new[] { x, C }, outs);
            var fcse2 = ca.Function("fcse", new[] { x, C }, outs, opts("cse", true));
            checkfunctionLight(f2, fcse2, dms(x0, data));
            string result = captureStdout(() => fcse2.disp(true));
            assertEquals(1, count(result, " F_der_der("));
            assertEquals(0, count(result, " F_der("));
            assertEquals(0, count(result, " F("));
            result = captureStdout(() => fcse.disp(true));
            assertEquals(0, count(result, " F_der_der("));
            assertEquals(1, count(result, " F_der("));
            assertEquals(0, count(result, " F("));
            if (parametricGrid) {
              MX yJ = ca.jacobian(F.Call(x, C, K), x);
              List<MX> pe = ca.extract_parametric(yJ, ca.vertcat(C, K)).Item3;
              if (precomputeCoeff || precomputeGrid) {
                assertTrue(pe.Count > 0);
                assertTrue(pe.Any(e => !e.is_zero()));
              } else {
                assertEquals(0, pe.Count);
              }
            } else {
              MX yJ = ca.jacobian(F.Call(x, C), x);
              List<MX> pe = ca.extract_parametric(yJ, C).Item3;
              if (precomputeCoeff) {
                assertTrue(pe.Count > 0);
                assertTrue(pe.Any(e => !e.is_zero()));
              } else {
                assertEquals(0, pe.Count);
              }
            }
            Function FJref = Fref.jacobian();
            Function FJtest, FJ = null;
            if (parametricGrid) {
              Function FJparam = F.jacobian();
              MX adjF = MX.sym("adj_f", 1, 1);
              List<MX> fjOut = FJparam.call(new[] { x, C, K, adjF });
              MX fjX = ca.substitute(fjOut[0], K, ca.MX(knotsData));
              MX fjC = ca.substitute(fjOut[1], K, ca.MX(knotsData));
              FJtest = ca.Function("FJ_test", new[] { x, C, adjF }, new[] { fjX, fjC });
            } else {
              FJ = F.jacobian();
              FJ.generate("FJ.c", opts("main", true));
              FJ.generate_in("FJ_in.txt", dms(vcat(last), data, 0.0));
              FJtest = FJ;
            }
            foreach (var a in allPoints) checkfunctionLight(FJtest, FJref, dms(vcat(a), data, 0.0));
            if (!parametricGrid) {
              message("ref");
              Fref2.call(dms(DM.zeros(N), data));
              FJ = FJ.jacobian();
              FJ.generate("FH.c", opts("main", true));
              FJ.generate_in("FH_in.txt", ins(x0, data, 0, 0, 0));
              message(FJ.ToString());
              FJref = FJref.jacobian();
              var sw = Stopwatch.StartNew();
              FJref.call(ins(vcat(last), data, 0, 0, 0));
              message("FJ_ref " + sw.Elapsed.TotalSeconds);
              sw.Restart();
              FJ.call(ins(vcat(last), data, 0, 0, 0));
              message("FJ " + sw.Elapsed.TotalSeconds);
              foreach (var a in allPoints) checkfunctionLight(FJ, FJref, ins(vcat(a), data, 0, 0, 0));
            }
            if (parametricGrid) {
              var altStacked = allKnotsAlt.Take(N).SelectMany(k => knotsExpand(k)).ToList();
              MX yAlt = ca.substitute(F.Call(x, C, K), K, ca.MX(vcat(altStacked)));
              var Falt = ca.Function("F_alt", new[] { x, C }, new[] { yAlt });
              double res1 = Ftest.Call(x0, data).to_double();
              double res2 = Falt.Call(x0, data).to_double();
              assertFalse(Math.Abs(res1 - res2) <= 1e-8 + 1e-5 * Math.Abs(res2), "Parametric knots should produce different results with different knot vectors");
            }
          }
          MX xLu = MX.sym("x", N);
          List<double[]> knotsUniformOrig = new List<double[]> { linspace(0, 1, 5), linspace(0, 1, 5), linspace(0, 1, 3), linspace(0, 1, 3), linspace(0, 1, 3) }.GetRange(0, N);
          double[][] knotsUniform = knotsUniformOrig.Select(knotsExpand).ToArray();
          int ncU = 1;
          foreach (double[] k in knotsUniform) ncU *= k.Length - 4;
          DM.rng(1);
          DM dataU = DM.rand(ncU);
          MX Cu = MX.sym("C", ncU, 1);
          MX Yu = ca.bspline(xLu, Cu, knotsUniform, fill(N, 3), 1);
          var FrefU = ca.Function("f", new[] { xLu, Cu }, new[] { Yu });
          List<List<double>> pointsU = testPoints(knotsUniformOrig);
          if (N > 3) pointsU = sample(pointsU, 125);
          foreach (string lookupMode in new[] { "linear", "exact", "binary" }) {
            var modes = Enumerable.Repeat(lookupMode, N).ToList();
            var optsLu = opts("precompute_coeff", precomputeCoeff, "precompute_grid", precomputeGrid, "lookup_mode", modes, "pedantic_mode_order", "ignore",
                "pedantic_mode_size", "ignore", "jit", true, "jit_options", blazingJitOptions());
            Function Flu = ca.blazing_spline("F", knotsUniform, optsLu);
            var FluTest = ca.Function("F_lu_test", new[] { xLu, Cu }, new[] { Flu.Call(xLu, Cu) });
            foreach (var a in pointsU) checkfunctionLight(FluTest, FrefU, dms(vcat(a), dataU));
            Function Flu2 = Function.deserialize(Flu.serialize(opts("debug", true)));
            foreach (var a in pointsU) checkfunctionLight(Flu, Flu2, dms(vcat(a), dataU));
          }
        }
      }
    }
  }

  static bool inEither(string s, string[] r) => r[0].Contains(s) || r[1].Contains(s);

  static double[][] range(params int[] n) => n.Select(k => Enumerable.Range(0, k).Select(i => (double)i).ToArray()).ToArray();

  public void test_blazing_spline_pedantic() {
    Check.Assume(CasadiMeta.feature_list().Contains("simde"), "no simde");
    double[][] sizeBad1d = range(12);
    double[][] orderBad2d = range(17, 13);
    double[][] cumulBad2d = range(8, 8);
    double[][] clean2d = range(13, 17);
    double[][] clean1d = range(13);
    assertInException("powers of 2", () => ca.blazing_spline("F", sizeBad1d, opts("pedantic_mode_size", "error")));
    assertInException("not increasing", () => ca.blazing_spline("F", orderBad2d, opts("pedantic_mode_order", "error", "pedantic_mode_size", "ignore")));
    assertInException("dim 0 (zero-based)", () => ca.blazing_spline("F", sizeBad1d, opts("pedantic_mode_size", "error")));
    assertInException("'pedantic_mode_size' must be one of", () => ca.blazing_spline("F", sizeBad1d, opts("pedantic_mode_size", "bogus")));
    assertInException("prefix product", () => ca.blazing_spline("F", cumulBad2d, opts("pedantic_mode_size", "error")));
    assertTrue(inEither("powers of 2", captureBoth(() => ca.blazing_spline("F", sizeBad1d, opts("pedantic_mode_size", "warn")))));
    assertTrue(inEither("not increasing", captureBoth(() => ca.blazing_spline("F", orderBad2d, opts("pedantic_mode_size", "ignore")))));
    assertFalse(inEither("powers of 2", captureBoth(() => ca.blazing_spline("F", sizeBad1d, opts("pedantic_mode_size", "ignore", "pedantic_mode_order", "ignore")))));
    ca.blazing_spline("F", clean1d, opts("pedantic_mode_order", "error", "pedantic_mode_size", "error"));
    ca.blazing_spline("F", clean2d, opts("pedantic_mode_order", "error", "pedantic_mode_size", "error"));
    ca.blazing_spline("F", clean1d, opts("pedantic_mode_order", "error"));
    double[][] n5bad1d = range(13);
    assertInException("(n_knots - 5)", () => ca.blazing_spline("F", n5bad1d).jacobian());
    assertInException("diff order 1", () => ca.blazing_spline("F", n5bad1d).jacobian());
    assertTrue(inEither("(n_knots - 5)", captureBoth(() => ca.blazing_spline("F", n5bad1d, opts("pedantic_mode_size", "warn")).jacobian())));
    assertFalse(inEither("(n_knots - 5)", captureBoth(() => ca.blazing_spline("F", n5bad1d, opts("pedantic_mode_size", "ignore")).jacobian())));
    assertFalse(inEither("(n_knots - 5)", captureBoth(() => ca.blazing_spline("F", n5bad1d, opts("jacobian_options", opts("pedantic_mode_size", "ignore"))).jacobian())));
    double[][] jacClean1d = range(15);
    ca.blazing_spline("F", jacClean1d, opts("pedantic_mode_size", "error")).jacobian();
    double[][] n6bad1d = range(14);
    assertInException("(n_knots - 6)", () => ca.blazing_spline("F", n6bad1d, opts("precompute_coeff", false)).jacobian().jacobian());
    assertInException("diff order 2", () => ca.blazing_spline("F", n6bad1d, opts("precompute_coeff", false)).jacobian().jacobian());
  }

  public void test_noncanonical_sparsity() {
    MX x = MX.sym("x", 4, 4);
    MX y = MX.sym("y");
    var f = ca.Function("f", new[] { x, y }, new[] { ca.MX(24.0) });
    checkCodegen(f, ins(DM.rand(4, 4), 1), new CG().Opt("force_canonical", false));
    checkCodegen(f, ins(DM.rand(4, 4), 1), new CG().Opt("force_canonical", true));
  }

  public void test_options_sanitize() {
    assertJ(opts("foo", opts("bar", opts("baz", opts("goot", 7L)))), Options.sanitize(opts("foo", opts("bar", opts("baz.goot", 7)))));
    assertJ(opts("foo", opts("bar", 9L, "baz", 1L)), Options.sanitize(opts("foo.baz", 1, "foo", opts("bar", 9))));
    assertJ(opts("foo", opts("baz", opts("bar", 1L, "go", 8L))), Options.sanitize(opts("foo.baz.bar", 1, "foo", opts("baz", opts("go", 8)))));
    assertJ(opts("foo", opts("bar", 9L, "baz", 1L)), Options.sanitize(opts("foo", opts("bar", 9), "foo.baz", 1)));
    assertJ(opts("bar", opts("abc", opts("def", 7L)), "foo", opts("bar", opts("baz", 1L, "mat", 3L, "got", 9L, "qint", opts("pad", 5L)), "lib", 5L)),
        Options.sanitize(opts("foo.bar.baz", 1, "foo.bar.mat", 3, "foo.lib", 5, "foo.bar.qint.pad", 5, "foo", opts("bar", opts("got", 9)), "bar", opts("abc.def", 7))));
    assertInException("update_dict error", () => Options.sanitize(opts("foo", 9, "foo.baz", 1)));
    assertInException("update_dict error", () => Options.sanitize(opts("foo", opts("baz", opts("w", 5)), "foo.baz", 1)));
    Dictionary<string, object> res = Options.sanitize(opts("foo", 1, "bar", null));
    message(jstr(res));
    assertJ(opts("foo", 1L), res);
    res = Options.sanitize(opts("foo", 1, "bar", opts("r", null)));
    message(jstr(res));
    assertJ(opts("foo", 1L, "bar", opts("r", null)), res);
  }

  static List<List<object>> passes(params object[][] p) => p.Select(e => e.ToList()).ToList();

  public void test_simplify() {
    MX x = MX.sym("x", 5);
    MX y = MX.sym("y", 5);
    var f = ca.Function("test_func", new[] { x, y }, new[] { x + y });
    Function fSimplified = f.transform();
    assertEquals("test_func", fSimplified.name());
    Function fSameName = f.transform("test_func");
    assertEquals("test_func", fSameName.name());
    Function fRenamed = f.transform("new_name");
    assertEquals("new_name", fRenamed.name());
    checkfunctionLight(f, fSimplified, ins(2.0, 3.0));
    x = MX.sym("x", 5);
    y = MX.sym("y", 5);
    f = ca.Function("test_func", new[] { x, y }, new[] { x });
    fSimplified = f.transform(passes(new object[] { "simplify", "empty_inputs" }));
    assertEquals(5L, f.nnz_in(0));
    assertEquals(5L, fSimplified.nnz_in(0));
    assertEquals(0L, fSimplified.nnz_in(1));
    checkfunctionLight(f, fSimplified, ins(2.0, 3.0));
  }

  public void test_simplify_ref_count() {
    var tests = new List<(Function, Function)>();
    foreach (string X in XS) {
      if (X == "SX") {
        SX a = SX.sym("a"), b = SX.sym("b");
        tests.Add((ca.Function("f", new[] { a, b }, new[] { -(a - b) }), ca.Function("f", new[] { a, b }, new[] { b - a })));
        SX e = a - b;
        SX z = -e;
        tests.Add((ca.Function("f", new[] { a, b }, new[] { e * 3, -e }), ca.Function("f", new[] { a, b }, new[] { e * 3, z })));
      } else {
        MX a = MX.sym("a"), b = MX.sym("b");
        tests.Add((ca.Function("f", new[] { a, b }, new[] { -(a - b) }), ca.Function("f", new[] { a, b }, new[] { b - a })));
        MX e = a - b;
        MX z = -e;
        tests.Add((ca.Function("f", new[] { a, b }, new[] { e * 3, -e }), ca.Function("f", new[] { a, b }, new[] { e * 3, z })));
      }
    }
    MX x = MX.sym("x", Sparsity.lower(3));
    MX y = ca.project(x, Sparsity.upper(3));
    tests.Add((ca.Function("f", new[] { x }, new[] { y.get_nz(false, 5) }), ca.Function("f", new[] { x }, new[] { x.get_nz(false, 5) })));
    x = MX.sym("x", Sparsity.lower(3));
    y = ca.project(x, Sparsity.upper(3));
    tests.Add((ca.Function("f", new[] { x }, new[] { y[.., -1] }), ca.Function("f", new[] { x }, new[] { x.nzref(Sparsity.dense(3, 1), new long[] { -1, -1, 5 }) })));
    foreach (var (f, refF) in tests) {
      refF.disp(true);
      Function fs = f.transform();
      fs.disp(true);
      fs.generate("f.c");
      refF.generate("ref.c");
      checkfunctionIdentical(fs, refF);
      DM.rng(1);
      checkfunctionLight(f, fs, randInputs(f));
    }
  }

  public void test_simplify_const_folding() {
    DM.rng(1);
    var A = ca.MX(DM.rand(2, 2));
    var B = ca.MX(DM.rand(2, 2));
    MX C = ca.mtimes(A, B);
    C = ca.mtimes(C, C);
    MX x = MX.sym("x");
    var f = ca.Function("f", new[] { x }, new[] { (x * (C * 2) + C * 2) * C[0] });
    Function fs = f.transform();
    checkfunctionLight(f, fs, randInputs(f));
    fs.disp(true);
    assertEquals(7L, fs.n_nodes());
  }

  public void test_simplify_combine_terms() {
    SX x = SX.sym("x");
    SX y = SX.sym("y");
    Action<SX, SX> numeq = (expr, refE) => {
      var fe = ca.Function("fe", new[] { x, y }, new[] { expr });
      var fr = ca.Function("fr", new[] { x, y }, new[] { refE });
      foreach (var (v0, v1) in new[] { (0.3, 1.7), (-2.0, 4.5) }) checkarray(fe.Call(v0, v1), fr.Call(v0, v1));
    };
    Func<SX, long> nnodes = expr => ca.Function("f", new[] { x, y }, new[] { expr }).n_nodes();
    SX e1 = y * 8 - y * 3;
    SX s1 = ca.simplify(e1);
    numeq(s1, y * 5);
    assertTrue(nnodes(s1) < nnodes(e1));
    SX e2 = x * 5 * 8 - (x + 3 + (x * 12 - x) - x * 17) * 3;
    SX s2 = ca.simplify(e2);
    numeq(s2, x * 55 - 9);
    assertTrue(nnodes(s2) < nnodes(e2));
    var f = ca.Function("f", new[] { x }, new[] { x * 8 - x * 3 });
    Function fc = f.transform(passes(new object[] { "simplify", "combine_terms" }));
    checkfunctionLight(f, fc, ins(2.0));
    assertTrue(fc.n_nodes() < f.n_nodes());
    MX a = MX.sym("a");
    MX b = MX.sym("b");
    MX e = -(a - b);
    assertTrue(ca.is_equal(ca.simplify(e), e, 2));
  }

  public void test_transform() {
    SX x = SX.sym("x");
    SX y = SX.sym("y");
    SX s = ca.transform(y * 8 - y * 3, passes(new object[] { "simplify", "combine_terms" }));
    var fe = ca.Function("fe", new[] { x, y }, new[] { s });
    var fr = ca.Function("fr", new[] { x, y }, new[] { y * 5 });
    var vals = new[] { (0.3, 1.7), (-2.0, 4.5) };
    foreach (var (v0, v1) in vals) checkarray(fe.Call(v0, v1), fr.Call(v0, v1));
    List<SX> sv = ca.transform(new[] { y * 8 - y * 3, x + x }, passes(new object[] { "simplify", "combine_terms" }));
    assertEquals(2, sv.Count);
    var fev = ca.Function("fev", new[] { x, y }, sv);
    var frv = ca.Function("frv", new[] { x, y }, new[] { y * 5, x * 2 });
    foreach (var (v0, v1) in vals) {
      checkarray(fev.call(dms(v0, v1))[0], frv.call(dms(v0, v1))[0]);
      checkarray(fev.call(dms(v0, v1))[1], frv.call(dms(v0, v1))[1]);
    }
    assertEquals(2, ca.transform(new[] { y * 8 - y * 3, x + x }, opts("combine_terms", true)).Count);
    MX a = MX.sym("a");
    MX b = MX.sym("b");
    var f = ca.Function("f", new[] { a, b }, new[] { -(a - b) });
    Function ft = f.transform(passes(new object[] { "simplify", 0, "ref_count", "const_folding" }));
    checkfunctionIdentical(ft, ca.Function("f", new[] { a, b }, new[] { b - a }));
    assertTrue(ft.n_nodes() < f.n_nodes());
    Function ft2 = f.transform(passes(new object[] { "simplify", 2, "ref_count" }));
    checkfunctionIdentical(ft2, ca.Function("f", new[] { a, b }, new[] { b - a }));
    assertEquals("f", f.transform().name());
    assertEquals("f", f.transform(passes(new object[] { "simplify", "cse" })).name());
    assertEquals("g", f.transform("g").name());
    assertEquals("g", f.transform("g", passes(new object[] { "simplify", "cse" })).name());
    assertEquals("g", f.transform("g", opts("cse", true)).name());
    Function frr = f.transform("g", passes(new object[] { "simplify", 0, "ref_count", "const_folding" }));
    checkfunctionIdentical(frr, ca.Function("g", new[] { a, b }, new[] { b - a }));
    Function fdNoarg = f.transform();
    Function fdEmpty = f.transform(opts());
    Function fdBool = f.transform(opts("combine_terms", false));
    foreach (Function fd in new[] { fdNoarg, fdEmpty, fdBool }) checkfunctionIdentical(fd, ca.Function("f", new[] { a, b }, new[] { b - a }));
    Function fdNone = f.transform(passes());
    checkfunctionIdentical(fdNone, f);
    SX sx = SX.sym("s");
    var g = ca.Function("g", new[] { sx }, new[] { sx * 8 - sx * 3 });
    Function gc = g.transform(opts("combine_terms", true));
    checkfunctionLight(g, gc, ins(2.0));
    assertTrue(gc.n_nodes() < g.n_nodes());
    Check.Throws<Exception>(() => g.transform(opts("bogus", true)));
    Check.Throws<Exception>(() => g.transform(passes(new object[] { "simplify", "cse" }), opts("cse", true)));
    Function fx = f.transform(passes(new object[] { "expand" }));
    assertEquals("SXFunction", fx.class_name());
    checkfunctionLight(f, fx, ins(1.0, 2.0));
    Check.Throws<Exception>(() => f.transform(passes(new object[] { "bogus" })));
    Check.Throws<Exception>(() => f.transform(passes(new object[] { "simplify", "nope" })));
    Check.Throws<Exception>(() => f.transform(passes(new object[] { "simplify", -1, "ref_count" })));
    if (!System.Runtime.InteropServices.RuntimeInformation.IsOSPlatform(System.Runtime.InteropServices.OSPlatform.OSX)) {
      string libcasadi = CasadiMeta.shared_library_prefix() + "casadi";
      var g2 = ca.Function("g", new[] { x }, new[] { ca.pow(x, 2) });
      Function gt = null;
      string[] r0 = captureBoth(() => gt = g2.transform(passes(new object[] { "external", libcasadi, "external_transform_test_success", opts("foo", "bar") })));
      foreach (string e in new[] { "Doing", "bar" }) assertTrue(r0[0].Contains(e), r0[0] + "<->" + e);
      assertTrue(r0[1].Contains("Warning"), r0[1] + "<->Warning");
      checkfunction(g2, gt, ins(3));
    }
    string[] r = captureBoth(() => f.transform(passes(new object[] { "simplify", 0, "ref_count" }), opts("verbose", true)));
    foreach (string e in new[] { "transform(", "n_nodes" }) assertTrue(r[0].Contains(e), r[0] + "<->" + e);
    r = captureBoth(() => f.transform(opts("verbose", true)));
    foreach (string e in new[] { "transform(", "n_nodes" }) assertTrue(r[0].Contains(e), r[0] + "<->" + e);
  }

  public void test_duplicate_check() {
    var x = ca.MX();
    var f = ca.Function("f", new[] { x }, new[] { x * 2 });
    assertEquals(0L, f.numel_in(0));
    f = ca.Function("f", new[] { MX.zeros(0, 0) }, new[] { x * 2 });
    assertEquals(0L, f.numel_in(0));
    f = ca.Function("f", new[] { MX.ones(0, 0) }, new[] { x * 2 });
    assertEquals(0L, f.numel_in(0));
    MX p = MX.sym("p");
    assertInException("input arguments must be purely symbolic", () => {
      var g = ca.Function("f", new[] { x.get_nz(false, p) }, new[] { x * 2 });
      assertEquals(0L, g.numel_in(0));
    });
    x.set_nz(3.7, false, p);
    assertInException("input arguments must be purely symbolic", () => {
      var g = ca.Function("f", new[] { x }, new[] { x * 3 });
      assertEquals(0L, g.numel_in(0));
    });
  }

  public void test_finite_diff() {
    MX x = MX.sym("x");
    foreach (string method in new[] { "forward", "backward", "central", "smoothing" }) {
      var f = ca.Function("f", new[] { x }, new[] { ca.sin(x) }, opts("enable_fd", true, "enable_forward", false, "enable_reverse", false, "fd_method", method));
      var fref = ca.Function("f", new[] { x }, new[] { ca.sin(x) });
      Function J = f.jacobian();
      Function Jref = fref.jacobian();
      List<DM> inputs = ins(1.1, 0);
      checkfunctionLight(J, Jref, inputs, 5);
      if (RunSlow) checkCodegen(J, inputs, new CG().Std("c99").Digits(15));
    }
  }

  public void test_finite_diff_serialization() {
    MX x = MX.sym("x", 2);
    MX p = MX.sym("p");
    MX y = p * ca.dot(x, x);
    foreach (string method in new[] { "forward", "backward", "central", "smoothing" }) {
      foreach (bool adaptive in new[] { false, true }) {
        double h = method == "backward" ? -1e-3 : 1e-3;
        var options = opts("h", h, "h_iter", (adaptive && (method == "central" || method == "smoothing")) ? 1 : 0, "h_min", 1e-8, "h_max", .1,
            "u_aim", 42.0, "reltol", 1e-11, "abstol", 1e-12, "smoothing", 5e-8);
        var f = ca.Function("f", new[] { x, p }, new[] { y, x + p }, opts("enable_fd", true, "enable_forward", false, "enable_reverse", false, "is_diff_in",
            new[] { true, false }, "is_diff_out", new[] { true, false }, "fd_method", method, "fd_options", options));
        var values = dms(ca.vertcat(.4, -.7), 2);
        foreach (int cnt in new[] { 1, 3 }) {
          Function df = f.forward(cnt);
          Function restored0 = Function.deserialize(df.serialize());
          assertEquals(df.class_name(), restored0.class_name());
          assertTrue(df.is_diff_in().SequenceEqual(restored0.is_diff_in()));
          assertTrue(df.is_diff_out().SequenceEqual(restored0.is_diff_out()));
          var inputs = new List<DM>(values);
          inputs.AddRange(f.call(values));
          inputs.Add(DM.ones(2, cnt));
          inputs.Add(DM.zeros(1, cnt));
          List<DM> got0 = restored0.call(inputs), expected0 = df.call(inputs);
          for (int i = 0; i < got0.Count; i++) checkarray(got0[i], expected0[i], 1e-12, "");
          for (long i = 0; i < df.n_in(); i++) assertTrue(restored0.sparsity_in(i).is_equal(df.sparsity_in(i)));
        }
        MX result = f.call(new[] { x, p })[0];
        var D = ca.Function("D", new[] { x, p }, new[] { ca.jacobian(result, x), ca.hessian(result, x).Item1, ca.jacobian(result, p) });
        Function restored = Function.deserialize(D.serialize());
        List<DM> got = restored.call(values), expected = D.call(values);
        for (int i = 0; i < got.Count; i++) checkarray(got[i], expected[i], 1e-12, "");
        if (!adaptive) checkCodegen(restored, values, new CG().Std("c99").Digits(10));
      }
    }
  }

  public void test_is_diff_fd() {
    DM.set_precision(16);
    try {
      foreach (int n in new[] { 1, 2 }) {
        foreach (int m in new[] { 1, 3 }) {
          MX x = MX.sym("x", n, m);
          MX y = MX.sym("y", n, m);
          var expr = new[] { ca.sumsqr(x) * 2 * ca.sumsqr(y), ca.pow(y, 2) * ca.sumsqr(x) };
          var fref = ca.Function("f", new[] { x, y }, expr, new[] { "x", "y" }, new[] { "a", "b" }, opts("is_diff_in", new[] { true, false }));
          var f = ca.Function("f", new[] { x, y }, expr, new[] { "x", "y" }, new[] { "a", "b" }, opts("is_diff_in", new[] { true, false }, "enable_fd", true,
              "print_in", true, "print_out", false, "enable_forward", false, "enable_reverse", false, "fd_method", "forward"));
          Function Jref = fref.jacobian();
          Function J = f.jacobian();
          assertTrue(Jref.is_diff_in(0));
          assertFalse(Jref.is_diff_in(1));
          assertEquals(0L, Jref.nnz_in(2));
          assertEquals(0L, J.nnz_in(2));
          Function Ffref = fref.forward(2);
          Function Ff = f.forward(2);
          assertEquals(0L, Ffref.nnz_in(2));
          assertEquals(1L, Ff.nnz_in(2));
          foreach (var (F, Fref) in new[] { (J, Jref), (Ff, Ffref) }) {
            message("F = " + F);
            message("Fref = " + Fref);
            assertEquals(Fref.n_in(), F.n_in());
            assertEquals(Fref.n_out(), F.n_out());
            for (long i = 0; i < F.n_in(); i++) {
              assertEquals(Fref.name_in(i), F.name_in(i));
              if (!F.name_in(i).StartsWith("out_")) assertTrue(F.sparsity_in(i).is_equal(Fref.sparsity_in(i)));
              assertEquals(Fref.is_diff_in(i), F.is_diff_in(i));
            }
            DM.rng(1);
            var inputs = new List<DM>();
            for (long i = 0; i < f.n_in(); i++) inputs.Add(DM.rand(F.sparsity_in(i)));
            inputs.AddRange(f.call(inputs));
            if (F.name().Contains("fwd")) for (long i = 0; i < f.n_in(); i++) inputs.Add(DM.rand(F.sparsity_in(i)));
            if (F.name().Contains("adj")) for (long i = 0; i < f.n_out(); i++) inputs.Add(DM.rand(F.sparsity_out(i)));
            checkfunctionLight(F, Fref, inputs, 5);
            if (RunSlow) checkCodegen(F, inputs, new CG().Std("c99").Digits(System.Runtime.InteropServices.RuntimeInformation.IsOSPlatform(System.Runtime.InteropServices.OSPlatform.OSX) ? 5 : 15));
            List<DM> in2 = ins(2, 3, 0, 0);
            if (F.name().Contains("fwd")) for (long i = 0; i < f.n_in(); i++) in2.Add(DM.rand(F.sparsity_in(i)));
            if (F.name().Contains("adj")) for (long i = 0; i < f.n_out(); i++) in2.Add(DM.rand(F.sparsity_out(i)));
            string result = captureStdout(() => F.call(in2));
            assertTrue(result.Contains("2.00"), result);
            assertFalse(result.Contains("3.00"), result);
          }
        }
      }
    } finally {
      DM.set_precision(6);
    }
  }

  public void test_forward_reverse_zero() {
    SX x = SX.sym("x");
    SX y = SX.sym("y", 2);
    var f = ca.Function("f", new[] { x, y }, new[] { x * y, ca.sin(x) });
    foreach (Function F in new[] { f, f.wrap() }) {
      Function fwd = F.forward(0);
      assertEquals((2L, 0L), fwd.size_in(5));
      assertEquals((2L, 0L), fwd.size_out(0));
      Function adj = F.reverse(0);
      assertEquals((2L, 0L), adj.size_in(4));
      assertEquals((2L, 0L), adj.size_out(1));
    }
  }

  public void test_jit_function_fd() {
    Check.Assume(Importer.has_plugin("shell"), "no importer shell");
    Check.Assume(RunSlow, "slow");
    Function f = Function.jit("f", "b[0] = a[0]*a[0];", new[] { "a" }, new[] { "b" }, opts("compiler", "shell"));
    checkarray(f.jacobian().Call(3, 0), ca.DM(6), 1e-5, "");
  }

  public void test_has_function() {
    MX x = MX.sym("x");
    var f = ca.Function("f", new[] { x }, new[] { ca.pow(x, 2) }, opts("never_inline", true));
    var g = ca.Function("g", new[] { x }, new[] { f.Call(x) + 1 });
    assertTrue(g.get_function().SequenceEqual(new[] { "f" }));
    assertTrue(g.has_function("f"));
    assertFalse(g.has_function("h"));
  }

  public void test_string_ctor_hidden() {
    assertTrue(typeof(Function).GetConstructor(new[] { typeof(string) }) == null);
  }

  static bool[] empActivity(Function f, bool[] mask) {
    bool[] emp = new bool[f.nnz_out()];
    for (int t = 1; t < 6; t++) {
      var inp = new List<DM>();
      int off = 0;
      for (long i = 0; i < f.n_in(); i++) {
        int ni = (int)f.nnz_in(i);
        double[] vals = new double[ni];
        for (int j = 0; j < ni; j++) vals[j] = mask[off + j] ? (double)t * (off + j + 1) : 0.0;
        inp.Add(ca.DM(f.sparsity_in(i), ni == 0 ? ca.DM(0, 1) : (DM)vals));
        off += ni;
      }
      var flat = f.call(inp).SelectMany(o => o.nonzeros()).ToList();
      for (int k = 0; k < emp.Length; k++) emp[k] = emp[k] || flat[k] != 0;
    }
    return emp;
  }

  static List<bool[]> allMasks(int n) {
    var r = new List<bool[]>();
    for (int b = 0; b < (1 << n); b++) r.Add(Enumerable.Range(0, n).Select(i => (b & (1 << i)) != 0).ToArray());
    return r;
  }

  static List<bool[]> prependTrue(int k, List<bool[]> masks) => masks.Select(m => Enumerable.Repeat(true, k).Concat(m).ToArray()).ToList();

  static void checkActivity(Function f, List<bool[]> masks, bool minimal) {
    foreach (bool[] m in masks) {
      bool[] emp = empActivity(f, m);
      List<bool> got = f.activity(m);
      string ms = "[" + string.Join(",", m) + "]";
      for (int k = 0; k < f.nnz_out(); k++)
        assertTrue(got[k] || !emp[k], "UNSOUND: " + f.name() + " mask " + ms + " output " + k + ": activity says inactive but empirically nonzero");
      if (minimal) assertTrue(emp.SequenceEqual(got), "not minimal: " + f.name() + " mask " + ms);
    }
  }

  public void test_activity() {
    foreach (string X in new[] { "MX", "SX" }) {
      if (X == "MX") {
        MX x = MX.sym("x"), y = MX.sym("y");
        checkActivity(ca.Function("mul", new[] { x, y }, new[] { x * y }), allMasks(2), true);
        checkActivity(ca.Function("add", new[] { x, y }, new[] { x + y }), allMasks(2), true);
        checkActivity(ca.Function("cs", new[] { x }, new[] { ca.cos(x) }), allMasks(1), true);
        checkActivity(ca.Function("sn", new[] { x }, new[] { ca.sin(x) }), allMasks(1), true);
        checkActivity(ca.Function("z", new[] { x }, new[] { x * 0.0 }), allMasks(1), true);
      } else {
        SX x = SX.sym("x"), y = SX.sym("y");
        checkActivity(ca.Function("mul", new[] { x, y }, new[] { x * y }), allMasks(2), true);
        checkActivity(ca.Function("add", new[] { x, y }, new[] { x + y }), allMasks(2), true);
        checkActivity(ca.Function("cs", new[] { x }, new[] { ca.cos(x) }), allMasks(1), true);
        checkActivity(ca.Function("sn", new[] { x }, new[] { ca.sin(x) }), allMasks(1), true);
        checkActivity(ca.Function("z", new[] { x }, new[] { x * 0.0 }), allMasks(1), true);
      }
    }
    MX x2 = MX.sym("x2", 2);
    checkActivity(ca.Function("pc", new[] { x2 }, new[] { x2 * ca.densify(ca.MX(ca.vertcat(0, 5))) }), allMasks(2), true);
    MX xr = MX.sym("xr", 1, 2), yc = MX.sym("yc", 2, 1);
    checkActivity(ca.Function("mt", new[] { xr, yc }, new[] { ca.mtimes(xr, yc) }), allMasks(4), true);
    MX A = MX.sym("A", 2, 2), b = MX.sym("b", 2, 1);
    checkActivity(ca.Function("slv", new[] { A, b }, new[] { ca.solve(A, b) }), prependTrue(4, allMasks(2)), true);
    MX ad = MX.sym("ad", 2), bd = MX.sym("bd", 2, 1);
    checkActivity(ca.Function("sd", new[] { ad, bd }, new[] { ca.solve(ca.diag(ad), bd) }), prependTrue(2, allMasks(2)), false);
    MX yv = MX.sym("y", 2);
    var foo = ca.Function("foo", new[] { yv }, new[] { ca.vertcat(ca.sin(yv[0] * yv[1]), ca.cos(yv[0] * yv[1])) });
    Function fadj = foo.reverse(1).forward(1);
    checkActivity(fadj, allMasks((int)fadj.nnz_in()), true);
    var args = new List<MX>();
    for (long i = 0; i < fadj.n_in(); i++) {
      if (fadj.name_in(i).StartsWith("fwd_")) args.Add(ca.MX(fadj.size1_in(i), fadj.size2_in(i)));
      else args.Add(MX.sym("a" + i, fadj.sparsity_in(i)));
    }
    assertTrue(fadj.call(args)[0].is_zero());
    var argsNz = new List<MX>();
    for (long i = 0; i < fadj.n_in(); i++) argsNz.Add(MX.sym("b" + i, fadj.sparsity_in(i)));
    assertFalse(fadj.call(argsNz)[0].is_zero());
    MX xx = MX.sym("x", 2);
    MX lam = MX.sym("lam", 2), sig = MX.sym("sig");
    MX L = sig * ca.sumsqr(xx - 2) + ca.dot(lam, ca.no_hess(foo.Call(ca.pow(xx, 2))));
    MX H = ca.hessian(L, xx).Item1;
    var Hf = ca.Function("H", new[] { xx, lam, sig }, new[] { H });
    assertEquals(0, Hf.get_function().Count);
    assertEquals(2L, H.nnz());
  }
}

}

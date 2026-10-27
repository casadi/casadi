package org.casadi.test;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.junit.jupiter.api.Assertions.assertNull;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.io.File;
import java.io.IOException;
import java.io.UncheckedIOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.casadi.CasadiMeta;
import org.casadi.CodeGenerator;
import org.casadi.DM;
import org.casadi.Function;
import org.casadi.GenericMatrix;
import org.casadi.GlobalOptions;
import org.casadi.Importer;
import org.casadi.MX;
import org.casadi.Operation;
import org.casadi.Opti;
import org.casadi.Options;
import org.casadi.SX;
import org.casadi.Slice;
import org.casadi.Sparsity;
import org.casadi.StringDeserializer;
import org.casadi.ca;
import org.junit.jupiter.api.Assumptions;
import org.junit.jupiter.api.Test;

/** Port of test/python/function.py (second part, from test_jit_directory) */
public class Function2Test extends CasadiTestCase {

  private static final List<String> XS = List.of("SX", "MX");

  // ---------- helpers ----------

  private static Map<String, Object> opts(Object... kv) {
    Map<String, Object> m = new LinkedHashMap<>();
    for (int i = 0; i < kv.length; i += 2) m.put((String) kv[i], kv[i + 1]);
    return m;
  }

  private static List<Double> jl(double... x) {
    List<Double> r = new ArrayList<>();
    for (double e : x) r.add(e);
    return r;
  }

  private static double tol(int digits) {
    return Math.pow(10, -digits);
  }

  private static String read(String path) {
    try {
      return new String(Files.readAllBytes(Paths.get(path)), StandardCharsets.UTF_8).replace("\r\n", "\n");
    } catch (IOException e) {
      throw new UncheckedIOException(e);
    }
  }

  private static void write(String path, String s) {
    try {
      Files.write(Paths.get(path), s.getBytes(StandardCharsets.UTF_8));
    } catch (IOException e) {
      throw new UncheckedIOException(e);
    }
  }

  private static Path tempDir() {
    try {
      return Files.createTempDirectory("casadi_java_test");
    } catch (IOException e) {
      throw new UncheckedIOException(e);
    }
  }

  private static void rmtree(Path p) {
    File f = p.toFile();
    if (!f.exists()) return;
    File[] c = f.listFiles();
    if (c != null) for (File e : c) rmtree(e.toPath());
    f.delete();
  }

  private static int countPrefix(String dir, String prefix) {
    File[] c = new File(dir).listFiles();
    int n = 0;
    if (c != null) for (File e : c) if (e.getName().startsWith(prefix)) n++;
    return n;
  }

  private static void removePrefix(String dir, String prefix) {
    File[] c = new File(dir).listFiles();
    if (c != null) for (File e : c) if (e.getName().startsWith(prefix)) e.delete();
  }

  /** Python json.loads; numbers become Double */
  private static Object json(String s) {
    int[] p = {0};
    Object v = jsonValue(s, p);
    jsonWs(s, p);
    if (p[0] != s.length()) throw new IllegalArgumentException("trailing json at " + p[0]);
    return v;
  }

  private static void jsonWs(String s, int[] p) {
    while (p[0] < s.length() && Character.isWhitespace(s.charAt(p[0]))) p[0]++;
  }

  private static Object jsonValue(String s, int[] p) {
    jsonWs(s, p);
    char c = s.charAt(p[0]);
    if (c == '{') {
      Map<String, Object> m = new LinkedHashMap<>();
      p[0]++;
      jsonWs(s, p);
      if (s.charAt(p[0]) == '}') {
        p[0]++;
        return m;
      }
      while (true) {
        jsonWs(s, p);
        String k = (String) jsonValue(s, p);
        jsonWs(s, p);
        p[0]++; // ':'
        m.put(k, jsonValue(s, p));
        jsonWs(s, p);
        if (s.charAt(p[0]++) == '}') return m;
      }
    }
    if (c == '[') {
      List<Object> l = new ArrayList<>();
      p[0]++;
      jsonWs(s, p);
      if (s.charAt(p[0]) == ']') {
        p[0]++;
        return l;
      }
      while (true) {
        l.add(jsonValue(s, p));
        jsonWs(s, p);
        if (s.charAt(p[0]++) == ']') return l;
      }
    }
    if (c == '"') {
      StringBuilder b = new StringBuilder();
      p[0]++;
      while (true) {
        char d = s.charAt(p[0]++);
        if (d == '"') return b.toString();
        if (d == '\\') {
          char e = s.charAt(p[0]++);
          switch (e) {
            case 'n': b.append('\n'); break;
            case 't': b.append('\t'); break;
            case 'r': b.append('\r'); break;
            case 'b': b.append('\b'); break;
            case 'f': b.append('\f'); break;
            case 'u':
              b.append((char) Integer.parseInt(s.substring(p[0], p[0] + 4), 16));
              p[0] += 4;
              break;
            default: b.append(e);
          }
        } else {
          b.append(d);
        }
      }
    }
    if (s.startsWith("true", p[0])) {
      p[0] += 4;
      return Boolean.TRUE;
    }
    if (s.startsWith("false", p[0])) {
      p[0] += 5;
      return Boolean.FALSE;
    }
    if (s.startsWith("null", p[0])) {
      p[0] += 4;
      return null;
    }
    int st = p[0];
    while (p[0] < s.length() && "+-0123456789.eE".indexOf(s.charAt(p[0])) >= 0) p[0]++;
    return Double.parseDouble(s.substring(st, p[0]));
  }

  /** Navigate parsed json: String keys for objects, Integer indices (negative from the end) for lists */
  private static Object at(Object o, Object... keys) {
    for (Object k : keys) {
      if (k instanceof String) {
        o = ((Map<?, ?>) o).get(k);
      } else {
        List<?> l = (List<?>) o;
        int i = (Integer) k;
        o = l.get(i < 0 ? l.size() + i : i);
      }
    }
    return o;
  }

  private static List<Object> jsonLines(String path) {
    List<Object> r = new ArrayList<>();
    for (String line : read(path).split("\n")) if (!line.isEmpty()) r.add(json(line));
    return r;
  }

  private static DM dmOf(Object jsonList) {
    List<?> l = (List<?>) jsonList;
    double[] v = new double[l.size()];
    for (int i = 0; i < v.length; i++) v[i] = (Double) l.get(i);
    return new DM(v);
  }

  private static Map<?, ?> graphBundle(String path) {
    String text = read(path);
    if (path.endsWith(".html")) {
      Matcher m = Pattern.compile("<script id=\"graph-data\" type=\"application/json\">(.*?)</script>", Pattern.DOTALL).matcher(text);
      assertTrue(m.find());
      text = m.group(1);
    }
    return (Map<?, ?>) json(text);
  }

  private static Function unpack(Object source) {
    return new StringDeserializer((String) source).unpack_function();
  }

  private static void assertOutputs(List<String> included, List<String> excluded, Runnable r) {
    String out = captureStdout(r);
    for (String s : included) assertTrue(out.contains(s), "expected '" + s + "' in output:\n" + out);
    for (String s : excluded) assertFalse(out.contains(s), "did not expect '" + s + "' in output:\n" + out);
  }

  private static DM dmCall(Function f, Object... args) {
    return (DM) f.invoke(args);
  }

  @SuppressWarnings("unchecked")
  private static List<DM> dmCallN(Function f, Object... args) {
    return (List<DM>) f.invoke(args);
  }

  private static List<DM> randInputs(Function f) {
    List<DM> r = new ArrayList<>();
    for (long i = 0; i < f.n_in(); i++) r.add(DM.rand(f.sparsity_in(i)));
    return r;
  }

  private static void checkfunctionLight(Function trial, Function solution, List<DM> inputs) {
    checkfunctionLight(trial, solution, inputs, 9);
  }

  private static void checkfunctionLight(Function trial, Function solution, List<DM> inputs, int digits) {
    assertEquals(solution.n_out(), trial.n_out());
    assertEquals(solution.n_in(), trial.n_in());
    for (int rep = 0; rep < 2; rep++) {
      List<DM> a = trial.call(inputs), b = solution.call(inputs);
      for (int i = 0; i < trial.n_out(); i++) {
        if (trial.sparsity_out(i).is_empty() || solution.sparsity_out(i).is_empty()) continue;
        checkarray(a.get(i), b.get(i), tol(digits), "output " + i);
      }
    }
  }

  private static void checkfunctionLight(Function trial, Function solution, Map<String, DM> inputs) {
    List<DM> in = new ArrayList<>();
    for (long i = 0; i < trial.n_in(); i++) in.add(inputs.containsKey(trial.name_in(i)) ? inputs.get(trial.name_in(i)) : new DM(0.0));
    checkfunctionLight(trial, solution, in);
  }

  private DM uniformFor(Sparsity sp) {
    double[] v = new double[(int) sp.nnz()];
    for (int k = 0; k < v.length; k++) v[k] = random.nextDouble();
    return new DM(sp, new DM(v));
  }

  /** checkfunction with outputs compared to tol, derivatives to tolSens, and uniform [0,1) seeds as in helpers.py */
  private void checkfunction(Function f, Function g, List<DM> inputs, double tol, double tolSens) {
    checkfunctionLight(f, g, inputs, (int) Math.round(-Math.log10(tol)));
    List<DM> a = f.call(inputs), b = g.call(inputs);
    List<DM> ja = new ArrayList<>(inputs), jb = new ArrayList<>(inputs);
    ja.addAll(a);
    jb.addAll(b);
    List<DM> va = f.jacobian().call(ja), vb = g.jacobian().call(jb);
    for (int i = 0; i < va.size(); i++) checkarray(va.get(i), vb.get(i), tolSens, "jacobian block " + i);
    List<DM> fwd = new ArrayList<>(ja), adj = new ArrayList<>(ja);
    for (long i = 0; i < f.n_in(); i++) fwd.add(uniformFor(f.sparsity_in(i)));
    for (long i = 0; i < f.n_out(); i++) adj.add(uniformFor(f.sparsity_out(i)));
    List<DM> fa = f.forward(1).call(fwd), fb = g.forward(1).call(fwd);
    for (int i = 0; i < fa.size(); i++) checkarray(fa.get(i), fb.get(i), tolSens, "forward " + i);
    List<DM> ra = f.reverse(1).call(adj), rb = g.reverse(1).call(adj);
    for (int i = 0; i < ra.size(); i++) checkarray(ra.get(i), rb.get(i), tolSens, "reverse " + i);
  }

  private static void checkfunctionIdentical(Function trial, Function solution) {
    trial.generate("f.c");
    String[] a = read("f.c").split("\n", -1);
    solution.generate("f.c");
    String[] b = read("f.c").split("\n", -1);
    for (int i = 0; i < Math.min(a.length, b.length); i++) assertEquals(b[i], a[i]);
    assertEquals(b.length, a.length);
  }

  private static List<DM> ins(Object... x) {
    List<DM> r = new ArrayList<>();
    for (Object e : x) r.add(e instanceof DM ? (DM) e : new DM(((Number) e).doubleValue()));
    return r;
  }

  private static DM sparseDM(double[][] rows) {
    return ca.sparsify(new DM(rows));
  }

  /** Options for checkCodegen, after helpers.py check_codegen */
  private static final class CG {
    Map<String, Object> opts = new LinkedHashMap<>();
    String std = "c89";
    boolean main;
    boolean mainOutputCheck = true;
    boolean withExternal = true;
    boolean withJacSparsity;
    boolean withForward;
    boolean debug;
    int digits = 15;
    List<String> extralibs = new ArrayList<>();
    List<String> extraOptions = new ArrayList<>();
    List<String> definitions = new ArrayList<>();
    Map<String, Object> externalOpts = new LinkedHashMap<>();

    CG std(String s) { std = s; return this; }
    CG main() { main = true; return this; }
    CG noMainCheck() { mainOutputCheck = false; return this; }
    CG noExternal() { withExternal = false; return this; }
    CG opt(String k, Object v) { opts.put(k, v); return this; }
    CG digits(int d) { digits = d; return this; }
    CG lib(String l) { extralibs.add(l); return this; }
    CG flag(String... f) { extraOptions.addAll(List.of(f)); return this; }
    CG def(String d) { definitions.add(d); return this; }
    CG jacSparsity() { withJacSparsity = true; return this; }
    CG forward() { withForward = true; return this; }
    CG debug() { debug = true; return this; }
    CG externalOpts(Map<String, Object> o) { externalOpts = o; return this; }
  }

  private static String includeDir() {
    String inc = GlobalOptions.getCasadiIncludePath();
    if (inc.isEmpty()) inc = Paths.get(GlobalOptions.getCasadiPath()).getParent().resolve("include").toString();
    return inc;
  }

  private static int sh(String cmd, File in, File out) {
    try {
      ProcessBuilder pb = new ProcessBuilder("sh", "-c", cmd);
      if (in != null) pb.redirectInput(in);
      if (out != null) pb.redirectOutput(out);
      else pb.redirectOutput(ProcessBuilder.Redirect.INHERIT);
      pb.redirectError(ProcessBuilder.Redirect.INHERIT);
      pb.environment().put("LD_LIBRARY_PATH", GlobalOptions.getCasadiPath() + ":" + System.getenv().getOrDefault("LD_LIBRARY_PATH", ""));
      return pb.start().waitFor();
    } catch (IOException | InterruptedException e) {
      throw new RuntimeException(e);
    }
  }

  private String compileCommand(String name, String source, CG o, boolean shared) {
    String libdir = GlobalOptions.getCasadiPath(), inc = includeDir();
    StringBuilder libs = new StringBuilder();
    for (String l : o.extralibs) {
      if (l.startsWith("-") || l.contains(".")) libs.append(" ").append(l);
      else libs.append(" -l").append(l);
    }
    StringBuilder defs = new StringBuilder();
    for (String d : o.definitions) defs.append(" -D").append(d);
    return "gcc -pedantic -std=" + o.std + " -fPIC " + (shared ? "-shared" : "") + " -Wall -Werror -Wextra -I" + inc + " -I" + inc
        + "/highs -Wno-dangling-pointer -Wno-unknown-warning-option -Wno-unknown-pragmas -Wno-long-long -Wno-unused-parameter "
        + (o.debug ? "-O0 -g" : "-O3") + defs + " " + source + " -o " + name + (shared ? ".so" : "") + " -L" + libdir + " -Wl,-rpath,"
        + libdir + " -Wl,-rpath,." + (shared ? "" : " -lm") + libs + " " + String.join(" ", o.extraOptions);
  }

  /** helpers.py check_codegen: generate C, compile, compare against an external; a no-op without RUN_SLOW */
  private Map<String, Object> checkCodegen(Function F, List<DM> inputs, CG o) {
    Map<String, Object> ret = new HashMap<>();
    if (!RUN_SLOW) return ret;
    String name = "codegen_" + Long.toHexString(System.nanoTime()) + Integer.toHexString(random.nextInt(1 << 30));
    Map<String, Object> cgOpts = new LinkedHashMap<>(o.opts);
    if (o.main) cgOpts.put("main", true);
    CodeGenerator cg = new CodeGenerator(name, cgOpts);
    cg.add(F, o.withJacSparsity);
    if (o.withForward) cg.add(F.forward(1), o.withJacSparsity);
    cg.generate();
    String libname = "./" + name + ".so";
    assertEquals(0, sh(compileCommand(name, name + ".c", o, true), null, null), "compilation failed");
    Function F2 = o.withExternal ? ca.external(F.name(), libname, o.externalOpts) : null;
    if (o.main) {
      assertEquals(0, sh(compileCommand(name, name + ".c", o, false), null, null), "compilation failed");
      F.generate_in(F.name() + "_in.txt", inputs);
      int rc = sh("./" + name + " " + F.name(), new File(F.name() + "_in.txt"), new File(F.name() + "_out.txt"));
      assertEquals(0, rc, "main return code");
    }
    List<DM> out = F.call(inputs);
    if (o.main && o.mainOutputCheck) {
      List<DM> outputs = F.generate_out(F.name() + "_out.txt");
      for (int i = 0; i < F.n_out(); i++) checkarray(out.get(i), outputs.get(i), tol(o.digits), "main output " + i);
    }
    ret.put("libname", libname);
    if (F2 == null) return ret;
    List<DM> out2 = F2.call(inputs);
    for (int i = 0; i < F.n_out(); i++) checkarray(out.get(i), out2.get(i), tol(o.digits), "external output " + i);
    ret.put("F", F2);
    return ret;
  }

  private Function compileExternal(String name, String source, Map<String, Object> opts, CG o) {
    assertEquals(0, sh(compileCommand(name, source, o, true), null, null), "compilation failed");
    externalLibname = "./" + name + ".so";
    return ca.external(name, externalLibname, opts);
  }

  private String externalLibname;

  private static void checkSerialize(Function F, List<DM> inputs) {
    Function F2 = Function.deserialize(F.serialize(opts("debug", true)));
    List<DM> a = F.call(inputs), b = F2.call(inputs);
    for (int i = 0; i < F.n_out(); i++) checkarray(a.get(i), b.get(i), 1e-16, "serialize " + i);
  }

  // ---------- tests ----------

  @Test
  public void test_jit_directory() throws IOException {
    Assumptions.assumeTrue(RUN_SLOW);
    message(CasadiMeta.feature_list());
    boolean haveMkstemps = CasadiMeta.compiler_flags().contains("HAVE_MKSTEMPS");
    boolean haveSimpleMkstemps = CasadiMeta.compiler_flags().contains("HAVE_SIMPLE_MKSTEMPS");
    boolean createDirs = !CasadiMeta.feature_list().contains("ghc-filesystem");
    String absTempDir = System.getProperty("java.io.tmpdir");
    String relTempDir = "jit_dir_test";
    rmtree(Paths.get(relTempDir));
    Files.createDirectories(Paths.get("park"));
    if (createDirs) {
      Files.createDirectories(Paths.get(relTempDir, "foo"));
      Files.createDirectories(Paths.get("foo"));
      Files.createDirectories(Paths.get("/tmp/foo"));
    }
    for (String tempDir : List.of("./", relTempDir, absTempDir)) {
      GlobalOptions.setTempWorkDir(tempDir);
      for (String directory : List.of("", "foo", absTempDir)) {
        message("directory " + directory);
        MX x = MX.sym("x");
        Map<String, Object> o = opts("jit", true, "compiler", "shell", "jit_options", opts("verbose", true, "directory", directory));
        if (!haveMkstemps && !haveSimpleMkstemps) {
          assertInException("custom temporary directory", () -> new Function("f", List.of(x), List.of(x.power(2.0)), o));
          continue;
        }
        String dir = Paths.get(directory).isAbsolute() ? directory : Paths.get(tempDir, directory).toString();
        removePrefix(dir, "tmp_casadi");
        Function f = new Function("f", List.of(x), List.of(x.power(2.0)), o);
        assertTrue(countPrefix(dir, "tmp_casadi") > 1);
        // Java cannot chdir("park") for the ghc-filesystem variant
        f.delete();
        System.gc();
        assertEquals(0, countPrefix(dir, "tmp_casadi"));
      }
      GlobalOptions.setTempWorkDir("./");
    }
  }

  @Test
  public void test_custom_jacobian() {
    MX x = MX.sym("x");
    MX p = MX.sym("p");
    Function J = new Function("jac_Q", List.of(x, p, new MX(1, 1), new MX(1, 1)), List.of(x.times(Math.PI), new MX(1.0), new MX(1.0), new MX(1.0)),
        List.of("x", "p", "out_r", "out_s"), List.of("jac_r_x", "jac_r_p", "jac_s_x", "jac_s_p"));
    Function f = new Function("Q", List.of(x, p), List.of(x.power(2.0), x.times(p).times(2.0)), List.of("x", "p"), List.of("r", "s"),
        opts("custom_jacobian", J, "jac_penalty", 0));
    J = null;
    List<MX> g = f.call(List.of(x, p), false, true);
    Function JJ = new Function("JJ", List.of(x, p), List.of(ca.jacobian(g.get(0).plus(x.times(3.0)), x)));
    checkarray(dmCall(JJ, 5, 2), 5 * Math.PI + 3);
  }

  @Test
  public void test_dump() {
    MX x = MX.sym("x", Sparsity.lower(3));
    MX y = MX.sym("y", 0, 0);
    MX z = MX.sym("z", 2, 2);
    for (String fmt : List.of("txt", "mtx")) {
      Function f = new Function("f", List.of(x, y, z), List.of(x.times(2.0), z.times(2.0)), List.of("x", "y", "z"), List.of("a", "c"),
          opts("dump", true, "dump_in", true, "dump_out", true, "dump_format", fmt));
      List<DM> in = List.of(sparseDM(new double[][] {{1, 0, 0}, {2, 4, 0}, {7, 8, 9}}), new DM(), new DM(new double[][] {{1, 3}, {4, 5}}));
      List<DM> out = f.call(in);
      Function F = Function.load("f.casadi");
      DM xNum = DM.from_file("f.000000.in.x." + fmt);
      DM yNum = DM.from_file("f.000000.in.y." + fmt);
      DM zNum = DM.from_file("f.000000.in.z." + fmt);
      checkarray(xNum, in.get(0));
      checkarray(yNum, in.get(1));
      checkarray(zNum, in.get(2));
      DM aNum = DM.from_file("f.000000.out.a." + fmt);
      DM cNum = DM.from_file("f.000000.out.c." + fmt);
      checkarray(aNum, out.get(0));
      checkarray(cNum, out.get(1));
      F.generate_in("test_in.txt", List.of(xNum, yNum, zNum));
      F.generate_out("test_out.txt", F.call(List.of(xNum, yNum, zNum)));
      DM X = DM.from_file("f.000000.in.txt");
      DM A = DM.from_file("f.000000.out.txt");
      DM Xr = DM.from_file("test_in.txt");
      DM Ar = DM.from_file("test_out.txt");
      checkarray(Xr, X);
      checkarray(Ar, A);
    }
    if (RUN_SLOW && CasadiMeta.feature_list().contains("ghc-filesystem")) {
      for (String d : List.of("dump_fn", "dump_fn_codegen")) rmtree(Paths.get(d));
      MX x2 = MX.sym("x", Sparsity.lower(3));
      MX z2 = MX.sym("z", 2, 2);
      Function f = new Function("f", List.of(x2, z2), List.of(x2.times(2.0), z2.times(2.0)), List.of("x", "z"), List.of("a", "c"),
          opts("dump_in", true, "dump_out", true, "dump_dir", "dump_fn"));
      List<DM> in = List.of(sparseDM(new double[][] {{1, 0, 0}, {2, 4, 0}, {7, 8, 9}}), new DM(new double[][] {{1, 3}, {4, 5}}));
      f.call(in);
      checkCodegen(f, in, new CG().std("c99").main().opt("dump_dir_suffix", "_codegen"));
      for (String name : List.of("f.000000.in.txt", "f.000000.in.x.mtx", "f.000000.in.z.mtx", "f.000000.out.a.mtx", "f.000000.out.c.mtx", "f.000000.out.txt"))
        assertEquals(read(Paths.get("dump_fn", name).toString()), read(Paths.get("dump_fn_codegen", name).toString()));
      for (String d : List.of("dump_fn", "dump_fn_codegen")) rmtree(Paths.get(d));
    }
  }

  @Test
  public void test_dump_trace() {
    message("Instruction dumps preserve intermediate values and share the dump counter");
    for (String X : XS) {
      Path directory = tempDir();
      try {
        String dir = directory.toString();
        Map<String, Object> o = opts("dump_trace", true, "dump_in", true, "dump_out", true, "dump_dir", dir);
        Function f;
        if (X.equals("SX")) {
          SX x = SX.sym("x", 2), p = SX.sym("p");
          f = new Function("traced", List.of(x, p), List.of(x.times(p).plus(x.power(2.0)), ca.sum1(x)), o);
        } else {
          MX x = MX.sym("x", 2), p = MX.sym("p");
          f = new Function("traced", List.of(x, p), List.of(x.times(p).plus(x.power(2.0)), ca.sum1(x)), o);
        }
        double[][] xvals = {{2, -3}, {1.1, 5}};
        double[] pvals = {4, -2};
        for (int count = 0; count < 2; count++) {
          double[] xval = xvals[count];
          double pval = pvals[count];
          List<DM> result = new ArrayList<>();
          String captured = captureStdout(() -> result.addAll(f.call(List.of(new DM(xval), new DM(pval)))));
          assertEquals("", captured);
          String prefix = Paths.get(dir, String.format("traced.%06d", count)).toString();
          List<Object> records = jsonLines(prefix + ".trace.jsonl");
          assertEquals("casadi_trace", at(records, 0, "format"));
          assertEquals(1.0, at(records, 0, "version"));
          assertEquals((double) count, at(records, 0, "dump_id"));
          assertEquals(Map.of("event", "inputs", "values", List.of(jl(xval), jl(pval))), records.get(1));
          assertEquals(Map.of("event", "end", "status", 0.0), at(records, -1));
          List<?> last = (List<?>) at(records, -2, "values");
          for (int i = 0; i < Math.min(result.size(), last.size()); i++) checkarray(new DM(result.get(i).nonzeros()), dmOf(last.get(i)));
          checkarray(ca.vertcat(result), DM.from_file(prefix + ".out.txt"));
          checkarray(ca.vertcat(List.of(new DM(xval), new DM(pval))), DM.from_file(prefix + ".in.txt"));
          List<Object> steps = records.subList(2, records.size() - 2);
          assertEquals(2 * f.n_instructions(), steps.size());
          Set<Long> checked = new HashSet<>();
          long add = Operation.OP_ADD.swigValue(), mul = Operation.OP_MUL.swigValue(), sq = Operation.OP_SQ.swigValue();
          for (int k = 0; k < f.n_instructions(); k++) {
            Object before = steps.get(2 * k), after = steps.get(2 * k + 1);
            long op = f.instruction_id(k);
            assertEquals((double) k, at(before, "instruction"));
            assertEquals((double) k, at(after, "instruction"));
            assertEquals((double) op, at(before, "op"));
            assertEquals((double) op, at(after, "op"));
            assertEquals("inputs", at(before, "phase"));
            assertEquals("outputs", at(after, "phase"));
            if (op == add || op == mul || op == sq) {
              DM a = dmOf(at(before, "values", 0));
              DM expected;
              if (op == sq) {
                expected = a.times(a);
              } else {
                DM b = dmOf(at(before, "values", 1));
                expected = op == add ? a.plus(b) : a.times(b);
              }
              checkarray(expected, dmOf(at(after, "values", 0)));
              checked.add(op);
            }
          }
          assertEquals(Set.of(add, mul, sq), checked);
        }
        Function restored = Function.deserialize(f.serialize());
        restored.reset_dump_count();
        restored.call(List.of(dm(2, -3), new DM(4)));
        assertEquals(0.0, at(jsonLines(Paths.get(dir, "traced.000000.trace.jsonl").toString()), 0, "dump_id"));
        f.change_option("dump_trace", false);
        f.call(List.of(dm(2, -3), new DM(4)));
        assertFalse(Files.exists(Paths.get(dir, "traced.000002.trace.jsonl")));
        f.change_option("dump_trace", true);
        f.call(List.of(dm(2, -3), new DM(4)));
        assertTrue(Files.exists(Paths.get(dir, "traced.000003.trace.jsonl")));
        CodeGenerator generator = new CodeGenerator("traced.c");
        generator.add(f);
        String tracedCode = generator.dump();
        f.change_option("dump_trace", false);
        generator = new CodeGenerator("traced.c");
        generator.add(f);
        assertEquals(tracedCode, generator.dump());
      } finally {
        rmtree(directory);
      }
    }
  }

  @Test
  public void test_dump_trace_options() {
    message("Reconstructed XFunctions preserve instruction dumping");
    for (String X : XS) {
      Path directory = tempDir();
      try {
        String dir = directory.toString();
        Map<String, Object> o = opts("dump_trace", true, "dump_dir", dir);
        Function f;
        Runnable jitTrace;
        if (X.equals("SX")) {
          SX x = SX.sym("x");
          f = new Function("original", List.of(x), List.of(x.times(x)), o);
          jitTrace = () -> new Function("jit_trace", List.of(x), List.of(x.times(x)), opts("jit", true, "dump_trace", true));
        } else {
          MX x = MX.sym("x");
          f = new Function("original", List.of(x), List.of(x.times(x)), o);
          jitTrace = () -> new Function("jit_trace", List.of(x), List.of(x.times(x)), opts("jit", true, "dump_trace", true));
        }
        Function expanded = f.expand("expanded");
        Function derivative = f.forward(1);
        expanded.invoke(3);
        derivative.invoke(3, 9, 1);
        Object[][] cases = {{expanded, 9.0}, {derivative, 6.0}};
        for (Object[] c : cases) {
          List<Object> records = jsonLines(Paths.get(dir, ((Function) c[0]).name() + ".000000.trace.jsonl").toString());
          assertEquals(List.of(List.of(c[1])), at(records, -2, "values"));
        }
        assertInException("dump_trace", jitTrace);
      } finally {
        rmtree(directory);
      }
    }
  }

  @Test
  public void test_dump_trace_sparse_call() {
    message("Instruction dumps handle sparse matrices and unused call outputs");
    for (String X : XS) {
      Path directory = tempDir();
      try {
        String dir = directory.toString();
        Map<String, Object> innerOpts = opts("never_inline", true, "dump_trace", true, "dump_dir", dir);
        Map<String, Object> o = opts("dump_trace", true, "dump_dir", dir);
        Function f;
        if (X.equals("SX")) {
          SX v = SX.sym("v", 3);
          Function inner = new Function("inner", List.of(v), List.of(v.plus(1.0), v.times(3.0)), innerOpts);
          SX x = SX.sym("x", Sparsity.diag(3));
          SX a = inner.call(List.of(ca.vertcat(List.of(x.get(0, 0), x.get(1, 1), x.get(2, 2))))).get(0);
          f = new Function("outer", List.of(x), List.of(x.times(2.0), a, new SX(0, 0), x), o);
        } else {
          MX v = MX.sym("v", 3);
          Function inner = new Function("inner", List.of(v), List.of(v.plus(1.0), v.times(3.0)), innerOpts);
          MX x = MX.sym("x", Sparsity.diag(3));
          MX a = inner.call(List.of(ca.vertcat(List.of(x.get(0, 0), x.get(1, 1), x.get(2, 2))))).get(0);
          f = new Function("outer", List.of(x), List.of(x.times(2.0), a, new MX(0, 0), x), o);
        }
        f.invoke(ca.diag(dm(1, 2, 3)));
        List<Object> records = jsonLines(Paths.get(dir, "outer.000000.trace.jsonl").toString());
        assertEquals(List.of(jl(1, 2, 3)), at(records, 1, "values"));
        List<?> vals = (List<?>) at(records, -2, "values");
        assertEquals(List.of(jl(2, 4, 6), jl(2, 3, 4)), vals.subList(0, 2));
        assertTrue(vals.get(2) == null || List.of().equals(vals.get(2)));
        assertEquals(jl(1, 2, 3), vals.get(3));
        double call = Operation.OP_CALL.swigValue();
        List<Object> calls = new ArrayList<>();
        for (Object r : records) if (Double.valueOf(call).equals(((Map<?, ?>) r).get("op"))) calls.add(r);
        assertEquals(2, calls.size());
        assertTrue(((List<?>) at(calls, 1, "values")).contains(null));
        List<Object> innerRecords = jsonLines(Paths.get(dir, "inner.000000.trace.jsonl").toString());
        assertEquals(Map.of("event", "end", "status", 0.0), at(innerRecords, -1));
      } finally {
        rmtree(directory);
      }
    }
  }

  @Test
  public void test_dump_trace_error() {
    message("Failed evaluations leave a readable partial trace");
    Path directory = tempDir();
    try {
      String dir = directory.toString();
      MX x = MX.sym("x");
      MX y = x.attachAssert(ca.gt(x, 0.0), "positive input required");
      Function f = new Function("checked", List.of(x), List.of(y.plus(1.0)), opts("dump_trace", true, "dump_dir", dir));
      assertInException("positive input required", () -> f.invoke(-1));
      List<Object> records = jsonLines(Paths.get(dir, "checked.000000.trace.jsonl").toString());
      assertEquals(Map.of("event", "error"), at(records, -1));
      checkarray(dmCall(f, 2), 3);
      records = jsonLines(Paths.get(dir, "checked.000001.trace.jsonl").toString());
      assertEquals(Map.of("event", "end", "status", 0.0), at(records, -1));
    } finally {
      rmtree(directory);
    }
  }

  @Test
  public void test_dump_trace_nonfinite() {
    message("Instruction dumps encode non-finite values as valid JSON");
    for (String X : XS) {
      Path directory = tempDir();
      try {
        String dir = directory.toString();
        Map<String, Object> o = opts("dump_trace", true, "dump_dir", dir, "print_instructions", true, "print_canonical", true);
        Function f;
        if (X.equals("SX")) {
          SX x = SX.sym("x");
          f = new Function("nonfinite", List.of(x), List.of(ca.rdivide(1.0, x), ca.rdivide(-1.0, x), ca.sqrt(x.minus(1.0))), o);
        } else {
          MX x = MX.sym("x");
          f = new Function("nonfinite", List.of(x), List.of(ca.rdivide(1.0, x), ca.rdivide(-1.0, x), ca.sqrt(x.minus(1.0))), o);
        }
        String captured = captureStdout(() -> f.invoke(0));
        assertTrue(captured.contains("nonfinite:"), captured);
        List<Object> records = jsonLines(Paths.get(dir, "nonfinite.000000.trace.jsonl").toString());
        assertEquals(List.of(List.of("inf"), List.of("-inf"), List.of("nan")), at(records, -2, "values"));
      } finally {
        rmtree(directory);
      }
    }
  }

  private static String[] run(String... cmd) {
    try {
      Process pr = new ProcessBuilder(cmd).start();
      String out = new String(pr.getInputStream().readAllBytes(), StandardCharsets.UTF_8);
      String err = new String(pr.getErrorStream().readAllBytes(), StandardCharsets.UTF_8);
      assertEquals(0, pr.waitFor());
      return new String[] {out, err};
    } catch (IOException | InterruptedException e) {
      throw new RuntimeException(e);
    }
  }

  private static boolean which(String exe) {
    return sh("command -v " + exe + " > /dev/null", null, null) == 0;
  }

  @Test
  public void test_export_graph() {
    message("HTML graph export embeds the serialized Function, consistent with dump_trace");
    for (String X : XS) {
      Path directory = tempDir();
      try {
        String dir = directory.toString();
        Map<String, Object> o = opts("dump_trace", true, "dump_dir", dir);
        Function f;
        if (X.equals("SX")) {
          SX x = SX.sym("x");
          f = new Function("graph", List.of(x), List.of(x.times(x).plus(x)), List.of("state"), List.of("result"), o);
        } else {
          MX x = MX.sym("x");
          f = new Function("graph", List.of(x), List.of(x.times(x).plus(x)), List.of("state"), List.of("result"), o);
        }
        String path = Paths.get(dir, "graph.html").toString();
        f.export_graph(path);
        Map<?, ?> bundle = graphBundle(path);
        assertEquals("casadi_viz", bundle.get("format"));
        assertEquals(1.0, bundle.get("version"));
        assertEquals("function", bundle.get("view"));
        assertEquals(CasadiMeta.version(), bundle.get("casadi_version"));
        assertFalse(bundle.containsKey("nodes"));
        Function restored = unpack(bundle.get("source"));
        assertEquals("graph", restored.name());
        assertEquals(X + "Function", restored.class_name());
        assertEquals(List.of("state"), restored.name_in());
        assertEquals(List.of("result"), restored.name_out());
        List<Long> ra = new ArrayList<>(), fa = new ArrayList<>();
        for (long k = 0; k < restored.n_instructions(); k++) ra.add(restored.instruction_id(k));
        for (long k = 0; k < f.n_instructions(); k++) fa.add(f.instruction_id(k));
        assertEquals(fa, ra);
        f.invoke(3);
        List<Object> records = jsonLines(Paths.get(dir, "graph.000000.trace.jsonl").toString());
        assertEquals(restored.name(), at(records, 0, "function"));
        assertEquals(restored.class_name(), at(records, 0, "type"));
        for (Object record : records) {
          Map<?, ?> r = (Map<?, ?>) record;
          if (!r.containsKey("instruction")) continue;
          assertEquals((double) restored.instruction_id(((Double) r.get("instruction")).longValue()), r.get("op"));
        }
        assertEquals(List.of(jl(12)), at(records, -2, "values"));
        String html = read(path);
        Matcher m = Pattern.compile("<script id=\"viz-config\" type=\"application/json\">(.*?)</script>", Pattern.DOTALL).matcher(html);
        assertTrue(m.find());
        Map<?, ?> config = (Map<?, ?>) json(m.group(1));
        String[] v = CasadiMeta.version().split("\\.");
        assertEquals("https://unpkg.com/@casadi/casadi-viz@" + v[0] + "." + v[1] + "/dist/index.js", config.get("viewer_url"));
        assertTrue(config.containsKey("runtime"));
        assertNull(config.get("runtime"));
        Function ff = f;
        assertInException("Invalid export_graph direction", () -> ff.export_graph(path, opts("direction", "sideways")));
        assertInException("Unknown export_graph option", () -> ff.export_graph(path, opts("unknown", true)));
      } finally {
        rmtree(directory);
      }
    }
  }

  @Test
  public void test_export_graph_calls() {
    message("Graph exports keep called Functions, by identity, and sparse metadata");
    for (String X : XS) {
      Path directory = tempDir();
      try {
        String dir = directory.toString();
        String path = Paths.get(dir, "calls.html").toString();
        Function f, nested, sparse;
        if (X.equals("SX")) {
          SX x = SX.sym("x");
          Function pair = new Function("pair", List.of(x), List.of(x.plus(1.0), x.minus(1.0)), opts("never_inline", true));
          List<SX> ab = pair.call(List.of(x));
          f = new Function("calls", List.of(x), List.of(ab.get(0).times(ab.get(1))));
          Function first = new Function("inner", List.of(x), List.of(x.plus(1.0)), opts("never_inline", true));
          Function second = new Function("inner", List.of(x), List.of(x.minus(2.0)), opts("never_inline", true));
          Function middle = new Function("middle", List.of(x), List.of(first.call(List.of(x)).get(0).plus(second.call(List.of(x)).get(0))), opts("never_inline", true));
          nested = new Function("nested", List.of(x), List.of(middle.call(List.of(x)).get(0).plus(first.call(List.of(x.plus(3.0))).get(0))));
          SX xs = SX.sym("x", Sparsity.lower(3));
          sparse = new Function("sparse_graph", List.of(xs), List.of(xs.times(2.0)));
        } else {
          MX x = MX.sym("x");
          Function pair = new Function("pair", List.of(x), List.of(x.plus(1.0), x.minus(1.0)), opts("never_inline", true));
          List<MX> ab = pair.call(List.of(x));
          f = new Function("calls", List.of(x), List.of(ab.get(0).times(ab.get(1))));
          Function first = new Function("inner", List.of(x), List.of(x.plus(1.0)), opts("never_inline", true));
          Function second = new Function("inner", List.of(x), List.of(x.minus(2.0)), opts("never_inline", true));
          Function middle = new Function("middle", List.of(x), List.of(first.call(List.of(x)).get(0).plus(second.call(List.of(x)).get(0))), opts("never_inline", true));
          nested = new Function("nested", List.of(x), List.of(middle.call(List.of(x)).get(0).plus(first.call(List.of(x.plus(3.0))).get(0))));
          MX xs = MX.sym("x", Sparsity.lower(3));
          sparse = new Function("sparse_graph", List.of(xs), List.of(xs.times(2.0)));
        }
        f.export_graph(path, opts("direction", "TB"));
        Map<?, ?> bundle = graphBundle(path);
        assertEquals("TB", bundle.get("direction"));
        Function restored = unpack(bundle.get("source"));
        assertEquals(List.of("pair"), restored.get_function());
        checkarray(dmCall(restored, 3), dmCall(f, 3));
        restored = unpack(((Map<?, ?>) json(nested.export_graph())).get("source"));
        List<String> names = new ArrayList<>();
        for (Function g : restored.find_functions()) names.add(g.name());
        names.sort(null);
        assertEquals(List.of("inner", "inner", "middle"), names);
        checkarray(dmCall(restored, 1.5), dmCall(nested, 1.5));
        sparse.export_graph(path);
        restored = unpack(graphBundle(path).get("source"));
        assertTrue(restored.sparsity_in(0).is_equal(Sparsity.lower(3)));
      } finally {
        rmtree(directory);
      }
    }
  }

  @Test
  public void test_export_graph_helpers() {
    message("Graph exports accept any serializable Function, with its helper Functions");
    Path directory = tempDir();
    try {
      String dir = directory.toString();
      MX x = MX.sym("x", 2), p = MX.sym("p");
      Map<String, MX> nlp = new LinkedHashMap<>();
      nlp.put("x", x);
      nlp.put("p", p);
      nlp.put("f", x.get(0).minus(p).power(2.0).plus(x.get(1).power(2.0)));
      nlp.put("g", x.get(0).plus(x.get(1)));
      Function solver = ca.nlpsol("solver", "sqpmethod", nlp, opts("qpsol", "qrqp", "print_header", false, "print_iteration", false, "print_time", false,
          "qpsol_options", opts("print_iter", false, "print_header", false)));
      SX s = SX.sym("s");
      Function mapped = new Function("f", List.of(s), List.of(ca.sin(s).times(s))).map(4);
      for (Function f : List.of(solver, mapped)) {
        for (String extension : List.of(".html", ".casadi_viz")) {
          String path = Paths.get(dir, "helpers" + extension).toString();
          f.export_graph(path);
          Map<?, ?> bundle = graphBundle(path);
          assertEquals("function", bundle.get("view"));
          Function restored = unpack(bundle.get("source"));
          assertEquals(f.class_name(), restored.class_name());
          assertEquals(f.get_function(), restored.get_function());
        }
        assertInException("DOT export requires an SXFunction or MXFunction", () -> f.export_graph(Paths.get(dir, "helpers.dot").toString()));
      }
      assertTrue(solver.get_function().contains("nlp_hess_l"));
    } finally {
      rmtree(directory);
    }
  }

  @Test
  public void test_export_graph_expressions() {
    message("Expression and Function exports default to expression and function views");
    Path directory = tempDir();
    try {
      String path = Paths.get(directory.toString(), "expressions.html").toString();
      for (String X : XS) {
        Function f;
        Runnable bad;
        if (X.equals("SX")) {
          SX x = SX.sym("actual_x", 2), y = SX.sym("actual_y", 2);
          ca.export_graph(x.plus(y), path);
          assertEquals("expression", graphBundle(path).get("view"));
          for (List<SX> expressions : List.of(List.of(x.plus(y), x.minus(y)), List.of(x.get(0), y.get(0)), List.of(new SX(7.0), SX.zeros(0, 2)))) {
            ca.export_graph(expressions, path);
            assertEquals("expression", graphBundle(path).get("view"));
          }
          ca.export_graph(List.of(x.plus(y), x.minus(y)), path, opts("view", "function"));
          assertEquals("function", graphBundle(path).get("view"));
          f = new Function("renamed", List.of(x, y), List.of(x.plus(y), x.minus(y)), List.of("left", "right"), List.of("sum", "difference"));
          bad = () -> ca.export_graph(x, path, opts("view", "other"));
        } else {
          MX x = MX.sym("actual_x", 2), y = MX.sym("actual_y", 2);
          ca.export_graph(x.plus(y), path);
          assertEquals("expression", graphBundle(path).get("view"));
          for (List<MX> expressions : List.of(List.of(x.plus(y), x.minus(y)), List.of(x.get(0), y.get(0)), List.of(new MX(7.0), MX.zeros(0, 2)))) {
            ca.export_graph(expressions, path);
            assertEquals("expression", graphBundle(path).get("view"));
          }
          ca.export_graph(List.of(x.plus(y), x.minus(y)), path, opts("view", "function"));
          assertEquals("function", graphBundle(path).get("view"));
          f = new Function("renamed", List.of(x, y), List.of(x.plus(y), x.minus(y)), List.of("left", "right"), List.of("sum", "difference"));
          bad = () -> ca.export_graph(x, path, opts("view", "other"));
        }
        f.export_graph(path);
        Map<?, ?> bundle = graphBundle(path);
        assertEquals("function", bundle.get("view"));
        Function restored = unpack(bundle.get("source"));
        assertEquals(List.of("left", "right"), restored.name_in());
        assertEquals(List.of("sum", "difference"), restored.name_out());
        Set<String> syms = new HashSet<>();
        if (X.equals("SX")) {
          for (SX v : ca.symvar(ca.veccat(restored.sx_in()))) syms.add(v.toString());
          assertEquals(Set.of("actual_x_0", "actual_x_1", "actual_y_0", "actual_y_1"), syms);
        } else {
          for (MX v : ca.symvar(ca.veccat(restored.mx_in()))) syms.add(v.toString());
          assertEquals(Set.of("actual_x", "actual_y"), syms);
        }
        f.export_graph(path, opts("view", "expression"));
        assertEquals("expression", graphBundle(path).get("view"));
        assertInException("Invalid export_graph view", bad);
      }
      ca.export_graph(List.of(), path);
    } finally {
      rmtree(directory);
    }
  }

  @Test
  public void test_export_graph_string() {
    message("Filename-free graph exports return the same bundle as file exports");
    for (String X : XS) {
      Path directory = tempDir();
      try {
        String path = Paths.get(directory.toString(), "graph.casadi_viz").toString();
        Function f;
        List<java.util.function.Function<Map<String, Object>, String>> toString = new ArrayList<>();
        List<java.util.function.BiConsumer<String, Map<String, Object>>> toFile = new ArrayList<>();
        if (X.equals("SX")) {
          SX x = SX.sym("x", 2);
          Function child = new Function("child", List.of(x), List.of(x.plus(1.0)), opts("never_inline", true));
          f = new Function("root", List.of(x), List.of(child.call(List.of(x)).get(0).div(x)));
          toString.add(o -> o == null ? ca.export_graph(x) : ca.export_graph(x, o));
          toFile.add((pa, o) -> ca.export_graph(x, pa));
          for (List<SX> e : List.of(List.of(x), List.of(x.get(0).plus(x.get(1)), x.get(0).minus(x.get(1))), List.of(new SX(2.0)), List.<SX>of())) {
            toString.add(o -> o == null ? ca.export_graph(e) : ca.export_graph(e, o));
            toFile.add((pa, o) -> ca.export_graph(e, pa));
          }
        } else {
          MX x = MX.sym("x", 2);
          Function child = new Function("child", List.of(x), List.of(x.plus(1.0)), opts("never_inline", true));
          f = new Function("root", List.of(x), List.of(child.call(List.of(x)).get(0).div(x)));
          toString.add(o -> o == null ? ca.export_graph(x) : ca.export_graph(x, o));
          toFile.add((pa, o) -> ca.export_graph(x, pa));
          for (List<MX> e : List.of(List.of(x), List.of(x.get(0).plus(x.get(1)), x.get(0).minus(x.get(1))), List.of(new MX(2.0)), List.<MX>of())) {
            toString.add(o -> o == null ? ca.export_graph(e) : ca.export_graph(e, o));
            toFile.add((pa, o) -> ca.export_graph(e, pa));
          }
        }
        String def = f.export_graph();
        assertEquals(def, f.export_graph(opts()));
        assertEquals("function", ((Map<?, ?>) json(def)).get("view"));
        for (Map<String, Object> o : List.of(opts(), opts("include_functions", false, "view", "expression", "direction", "LR"))) {
          String text = f.export_graph(o);
          f.export_graph(path, o);
          assertEquals(text + "\n", read(path));
          Map<?, ?> bundle = (Map<?, ?>) json(text);
          assertEquals(o.getOrDefault("include_functions", true), bundle.get("include_functions"));
          assertEquals(o.getOrDefault("direction", "TB"), bundle.get("direction"));
          assertEquals(((Map<?, ?>) json(def)).get("source"), bundle.get("source"));
        }
        assertEquals(def, f.export_graph(opts("viz_js", "not-needed.js")));
        for (int i = 0; i < toString.size(); i++) {
          String text = toString.get(i).apply(null);
          assertEquals("expression", ((Map<?, ?>) json(text)).get("view"));
          toFile.get(i).accept(path, null);
          assertEquals(text + "\n", read(path));
          assertEquals("function", ((Map<?, ?>) json(toString.get(i).apply(opts("view", "function")))).get("view"));
        }
        Function ff = f;
        assertInException("Invalid export_graph view", () -> ff.export_graph(opts("view", "invalid")));
        assertInException("Unknown export_graph option", () -> ff.export_graph(opts("unknown", true)));
      } finally {
        rmtree(directory);
      }
    }
  }

  @Test
  public void test_export_graph_formats() {
    message("HTML embeds the .casadi_viz bundle; DOT stays a standalone instruction graph");
    for (String X : XS) {
      Path directory = tempDir();
      try {
        String base = Paths.get(directory.toString(), "graph").toString();
        Function f, free;
        Runnable exprExport;
        if (X.equals("SX")) {
          SX x = SX.sym("x", 2);
          Function child = new Function("child", List.of(x), List.of(x.plus(1.0)), opts("never_inline", true));
          f = new Function("root", List.of(x), List.of(child.call(List.of(x)).get(0).div(x)), List.of("a"), List.of("result"));
          exprExport = () -> ca.export_graph(List.of(x.get(0).plus(x.get(1)), x.get(0).minus(x.get(1))), base + ".casadi_viz");
          SX fr = SX.sym("free");
          free = new Function("free_graph", List.of(x), List.of(x.plus(fr)), opts("allow_free", true));
        } else {
          MX x = MX.sym("x", 2);
          Function child = new Function("child", List.of(x), List.of(x.plus(1.0)), opts("never_inline", true));
          f = new Function("root", List.of(x), List.of(child.call(List.of(x)).get(0).div(x)), List.of("a"), List.of("result"));
          exprExport = () -> ca.export_graph(List.of(x.get(0).plus(x.get(1)), x.get(0).minus(x.get(1))), base + ".casadi_viz");
          MX fr = MX.sym("free");
          free = new Function("free_graph", List.of(x), List.of(x.plus(fr)), opts("allow_free", true));
        }
        f.export_graph(base + ".html");
        Map<?, ?> embedded = graphBundle(base + ".html");
        f.export_graph(base + ".casadi_viz", opts("viz_js", "not-needed.js"));
        Map<?, ?> bundle = graphBundle(base + ".casadi_viz");
        assertEquals(embedded, bundle);
        assertEquals(CasadiMeta.version(), bundle.get("casadi_version"));
        for (String extension : List.of(".html", ".casadi_viz")) {
          f.export_graph(base + extension, opts("include_functions", false));
          Map<?, ?> reduced = graphBundle(base + extension);
          assertEquals(Boolean.FALSE, reduced.get("include_functions"));
          assertEquals(bundle.get("source"), reduced.get("source"));
        }
        for (String view : List.of("function", "expression")) {
          f.export_graph(base + ".dot", opts("view", view, "viz_js", "not-needed.js"));
          String dot = read(base + ".dot");
          assertTrue(dot.startsWith("digraph G {"));
          assertFalse(dot.contains("<script"));
          assertFalse(dot.contains("arg0"));
          assertTrue(dot.contains(":w"));
          assertTrue(dot.contains(":e"));
          assertTrue(dot.contains("2-by-1"));
          assertEquals(view.equals("function"), dot.contains("result"));
          if (which("dot")) {
            String[] rendered = run("dot", "-Tsvg", base + ".dot");
            assertTrue(rendered[0].contains("<svg"));
            assertEquals("", rendered[1]);
          }
        }
        for (String extension : List.of("", ".json", ".casadi_graph")) {
          Function ff = f;
          assertInException("Supported extensions: .html, .dot, .casadi_viz", () -> ff.export_graph(base + extension));
        }
        exprExport.run();
        assertEquals("expression", graphBundle(base + ".casadi_viz").get("view"));
        free.export_graph(base + ".casadi_viz");
        Function restored = unpack(graphBundle(base + ".casadi_viz").get("source"));
        List<String> names = new ArrayList<>();
        if (X.equals("SX")) for (SX v : restored.free_sx()) names.add(v.toString());
        else for (MX v : restored.free_mx()) names.add(v.toString());
        assertEquals(List.of("free"), names);
      } finally {
        rmtree(directory);
      }
    }
  }

  @Test
  public void test_export_graph_escaping() {
    message("Serialized graphs and optional local Viz.js stay safely embedded in HTML");
    Path directory = tempDir();
    try {
      MX x = MX.sym("x");
      String name = "</script><script>alert(\"input\")</script>";
      Function f = new Function("escaped", List.of(x), List.of(x), List.of(name), List.of("result"));
      String source = "/* </script> & \"quoted\" \\ source */";
      String runtime = Paths.get(directory.toString(), "viz.js").toString();
      write(runtime, source);
      String path = Paths.get(directory.toString(), "escaped.html").toString();
      f.export_graph(path, opts("viz_js", runtime, "viewer_url", "https://example.com/viewer.js?a=</script>"));
      String html = read(path);
      assertFalse(html.contains(name));
      Function restored = unpack(graphBundle(path).get("source"));
      Matcher m = Pattern.compile("<script id=\"viz-config\" type=\"application/json\">(.*?)</script>", Pattern.DOTALL).matcher(html);
      assertTrue(m.find());
      assertNotNull(json(m.group(1)));
      assertEquals(name, restored.name_in(0));
    } finally {
      rmtree(directory);
    }
  }

  @Test
  public void test_print() {
    MX x = MX.sym("x", Sparsity.lower(3));
    MX z = MX.sym("z", 2, 2);
    Function f = new Function("f", List.of(x, z), List.of(x.times(2.0), z.times(2.0)), List.of("x", "z"), List.of("a", "c"),
        opts("print_in", true, "print_out", true, "print_canonical", true));
    List<DM> in = List.of(sparseDM(new double[][] {{1, 0, 0}, {2, 4, 0}, {7, 8, 9}}), new DM(new double[][] {{1, 3}, {4, 5}}));
    List<DM> result = new ArrayList<>();
    String vmOutput = captureStdout(() -> result.addAll(f.call(in)));
    checkarray(result.get(0), in.get(0).times(2.0));
    checkarray(result.get(1), in.get(1).times(2.0));
    List<String> vmLines = new ArrayList<>();
    for (String line : vmOutput.strip().split("\n")) vmLines.add(line.replaceAll(" \\(0x[0-9a-f]+\\)", ""));
    assertTrue(vmLines.size() >= 6);
    checkCodegen(f, in, new CG().std("c99").main().noMainCheck());
    if (RUN_SLOW) {
      List<String> cgLines = new ArrayList<>();
      for (String line : read(f.name() + "_out.txt").strip().split("\n"))
        if (line.startsWith("Function") || line.startsWith("Input") || line.startsWith("Output")) cgLines.add(line);
      assertEquals(vmLines, cgLines);
    }
  }

  @Test
  public void test_eval_shapes() {
    MX x = MX.sym("x", Sparsity.lower(3));
    MX y = MX.sym("y", 3, 1);
    MX z = MX.sym("z", 2, 2);
    Function f = new Function("f", List.of(x, y, z), List.of(x.times(2.0), y.times(2.0), z.times(2.0)), opts("default_in", new double[] {1, 2, 3}));
    double[][] L = {{1, 0, 0}, {2, 4, 0}, {7, 8, 9}};
    double[][] Z = {{1, 3}, {4, 5}};
    DM l = sparseDM(L);
    DM lt2 = l.times(2.0);
    DM y987 = dm(9, 8, 7);
    List<Object[]> cases = new ArrayList<>();
    // normal
    cases.add(new Object[] {ins(l, y987, dm(Z)), List.of(lt2, y987.times(2.0), dm(Z).times(2.0)), new double[] {1, 2, 7, 4, 8, 9, 9, 8, 7, 1, 4, 3, 5}});
    // Project
    cases.add(new Object[] {ins(sparseDM(new double[][] {{1, 0, 0}, {2, 0, 0}, {7, 8, 9}}), y987, dm(Z)),
        List.of(sparseDM(new double[][] {{1, 0, 0}, {2, 1e-30, 0}, {7, 8, 9}}).times(2.0), y987.times(2.0), dm(Z).times(2.0)),
        new double[] {1, 2, 7, 0, 8, 9, 9, 8, 7, 1, 4, 3, 5}});
    // Project away
    cases.add(new Object[] {ins(sparseDM(new double[][] {{1, 0, 1}, {2, 4, 0}, {7, 8, 9}}), y987, dm(Z)), List.of(lt2, y987.times(2.0), dm(Z).times(2.0)),
        new double[] {1, 2, 7, 4, 8, 9, 9, 8, 7, 1, 4, 3, 5}});
    // Scalar expansion
    cases.add(new Object[] {ins(new DM(4), 5, 7),
        List.of(sparseDM(new double[][] {{4, 0, 0}, {4, 4, 0}, {4, 4, 4}}).times(2.0), dm(5, 5, 5).times(2.0), dm(new double[][] {{7, 7}, {7, 7}}).times(2.0)),
        new double[] {4, 4, 4, 4, 4, 4, 5, 5, 5, 7, 7, 7, 7}});
    // Repmat
    cases.add(new Object[] {ins(l, y987, dm(1, 3)), List.of(lt2, y987.times(2.0), dm(new double[][] {{1, 1}, {3, 3}}).times(2.0)),
        new double[] {1, 2, 7, 4, 8, 9, 9, 8, 7, 1, 3, 1, 3}});
    // npar
    cases.add(new Object[] {ins(l, dm(new double[][] {{9, 2}, {8, 4}, {7, 4}}), dm(Z)),
        List.of(sparseDM(new double[][] {{1, 0, 0, 1, 0, 0}, {2, 4, 0, 2, 4, 0}, {7, 8, 9, 7, 8, 9}}).times(2.0), dm(new double[][] {{9, 2}, {8, 4}, {7, 4}}).times(2.0),
            dm(new double[][] {{1, 3, 1, 3}, {4, 5, 4, 5}}).times(2.0)),
        null});
    // growing npar
    DM z8 = dm(new double[][] {{1, 3, 4, 5, 1, 1, 1, 1}, {4, 5, 7, 8, 3, 3, 3, 3}});
    cases.add(new Object[] {ins(l, dm(new double[][] {{9, 2}, {8, 4}, {7, 4}}), z8),
        List.of(ca.repmat(l, 1, 4).times(2.0), ca.repmat(dm(new double[][] {{9, 2}, {8, 4}, {7, 4}}), 1, 2).times(2.0), z8.times(2.0)), null});
    // growing npar
    cases.add(new Object[] {ins(l, dm(new double[][] {{9, 2, 3, 4}, {8, 4, 7, 8}, {7, 4, 9, 10}}), dm(new double[][] {{1, 3, 4, 5}, {4, 5, 7, 8}})), null, null});
    // Transpose
    cases.add(new Object[] {ins(l, y987.T(), dm(Z)), List.of(lt2, y987.times(2.0), dm(Z).times(2.0)), new double[] {1, 2, 7, 4, 8, 9, 9, 8, 7, 1, 4, 3, 5}});
    // Null: note default_in applies only to dict style; the last two should be forbidden
    for (DM e : List.of(new DM(), new DM(3, 0), new DM(0, 3), new DM(0, 5), new DM(8, 0)))
      cases.add(new Object[] {ins(l, e, dm(1, 3)), List.of(lt2, dm(0, 0, 0).times(2.0), dm(new double[][] {{1, 1}, {3, 3}}).times(2.0)),
          new double[] {1, 2, 7, 4, 8, 9, 0, 0, 0, 1, 3, 1, 3}});
    // Wrong
    cases.add(new Object[] {ins(l, y987.T(), dm(1, 3).T()), null, null});
    cases.add(new Object[] {ins(sparseDM(new double[][] {{1, 0, 0}, {2, 4, 0}, {7, 8, 9}, {4, 5, 6}}), y987, dm(1, 3)), null, null});
    cases.add(new Object[] {ins(l, dm(new double[][] {{9, 2, 1}, {8, 4, 1}, {7, 4, 1}}), dm(new double[][] {{1, 1, 1, 3}, {3, 4, 5, 6}})), null, null});
    cases.add(new Object[] {ins(l, y987, dm(new double[][] {{1, 1, 1}, {3, 4, 5}})), null, null});
    cases.add(new Object[] {ins(l, dm(9, 8, 7, 6), dm(Z)), null, null});

    for (Object[] c : cases) {
      @SuppressWarnings("unchecked")
      List<DM> in = (List<DM>) c[0];
      @SuppressWarnings("unchecked")
      List<DM> ref = (List<DM>) c[1];
      double[] refFlat = (double[]) c[2];
      if (ref == null) {
        assertInException("mismatching shape", () -> f.call(in));
      } else {
        List<DM> res = f.call(in);
        for (int i = 0; i < f.n_in(); i++) {
          assertTrue(res.get(i).sparsity().is_equal(ref.get(i).sparsity()));
          checkarray(res.get(i), ref.get(i));
        }
        if (refFlat == null) {
          assertInException("mismatching shape", () -> f.nz_from_out(res));
        } else {
          double[] resFlat = f.nz_from_out(res);
          double[] twice = new double[refFlat.length];
          for (int i = 0; i < twice.length; i++) twice[i] = 2 * refFlat[i];
          checkarray(new DM(resFlat), new DM(twice));
          List<DM> resRe = f.nz_to_out(resFlat);
          for (int i = 0; i < f.n_out(); i++) {
            assertTrue(resRe.get(i).sparsity().is_equal(ref.get(i).sparsity()));
            checkarray(resRe.get(i), ref.get(i));
          }
          f.generate_out("test.txt", res);
          List<DM> res2 = f.generate_in("test.txt");
          for (int i = 0; i < f.n_out(); i++) {
            assertTrue(res2.get(i).sparsity().is_equal(ref.get(i).sparsity()));
            checkarray(res2.get(i), ref.get(i));
          }
        }
      }
      if (refFlat == null) {
        assertInException("mismatching shape", () -> f.nz_from_in(in));
        assertInException("mismatching shape", () -> f.generate_in("test.txt", in));
      } else {
        double[] res = f.nz_from_in(in);
        checkarray(new DM(res), new DM(refFlat));
        List<DM> inRe = f.nz_to_in(res);
        for (int i = 0; i < f.n_in(); i++) {
          assertTrue(inRe.get(i).sparsity().is_equal(ref.get(i).sparsity()));
          checkarray(inRe.get(i), ref.get(i).div(2.0));
        }
        f.generate_in("test.txt", in);
        List<DM> ins2 = f.generate_in("test.txt");
        for (int i = 0; i < f.n_in(); i++) {
          assertTrue(ins2.get(i).sparsity().is_equal(ref.get(i).sparsity()));
          checkarray(ins2.get(i), ref.get(i).div(2.0));
        }
      }
    }
    // Null dicts
    Map<String, DM> in = new LinkedHashMap<>();
    in.put("i0", l);
    in.put("i2", dm(1, 3));
    Map<String, DM> res = f.call(in);
    checkarray(res.get("o0"), lt2);
    checkarray(res.get("o1"), dm(2, 2, 2).times(2.0));
    checkarray(res.get("o2"), dm(new double[][] {{1, 1}, {3, 3}}).times(2.0));
    checkarray(new DM(f.nz_from_in(f.convert_in(in))), dm(1, 2, 7, 4, 8, 9, 2, 2, 2, 1, 3, 1, 3));
  }

  @Test
  public void test_issue3977() {
    for (String X : XS) {
      Function f, g;
      if (X.equals("SX")) {
        SX x = SX.sym("x", 0, 1);
        f = new Function("f", List.of(x), List.of(x.power(2.0)));
        SX x2 = SX.sym("x", 0, 2);
        g = new Function("f", List.of(x2), List.of(x2.power(2.0)));
      } else {
        MX x = MX.sym("x", 0, 1);
        f = new Function("f", List.of(x), List.of(x.power(2.0)));
        MX x2 = MX.sym("x", 0, 2);
        g = new Function("f", List.of(x2), List.of(x2.power(2.0)));
      }
      assertTrue(Arrays.equals(new long[] {0, 3}, dmCall(f, new DM(0, 3)).size()));
      assertTrue(Arrays.equals(new long[] {0, 1}, dmCall(f, new DM(0, 1)).size()));
      assertTrue(Arrays.equals(new long[] {0, 6}, dmCall(g, new DM(0, 6)).size()));
    }
    SX x = SX.sym("x", 3, 1);
    Function f = new Function("f", List.of(x), List.of(x.power(2.0)));
    assertTrue(Arrays.equals(new long[] {3, 1}, dmCall(f, new DM(0, 0)).size()));
    assertTrue(Arrays.equals(new long[] {3, 1}, dmCall(f, new DM(0, 1)).size()));
  }

  @Test
  public void test_convert_in() {
    MX x = MX.sym("x");
    MX y = MX.sym("y");
    MX z = MX.sym("z");
    Function f = new Function("f", List.of(x, y, z), List.of(x.times(2.0), y.times(2.0), z.times(2.0), x.plus(y)), List.of("x", "y", "z"), List.of("a", "b", "c", "d"),
        opts("default_in", new double[] {1, 2, 3}));
    List<DM> res = f.convert_in(Map.of("x", new DM(4)));
    checkarray(res.get(0), 4);
    checkarray(res.get(1), 2);
    checkarray(res.get(2), 3);
    assertInException("xy", () -> f.convert_in(Map.of("xy", new DM(4))));
    Map<String, DM> resd = f.convert_in(ins(7, 3, 4));
    checkarray(resd.get("x"), 7);
    checkarray(resd.get("y"), 3);
    checkarray(resd.get("z"), 4);
    assertInException("Incorrect", () -> f.convert_in(ins(7, 3)));
    res = f.convert_out(Map.of("b", new DM(4)));
    checkarray(res.get(0), Double.NaN);
    checkarray(res.get(1), 4);
    checkarray(res.get(2), Double.NaN);
    checkarray(res.get(3), Double.NaN);
    assertInException("xy", () -> f.convert_out(Map.of("xy", new DM(4))));
    resd = f.convert_out(ins(1, 2, 3, 4));
    checkarray(resd.get("a"), 1);
    checkarray(resd.get("b"), 2);
    checkarray(resd.get("c"), 3);
    checkarray(resd.get("d"), 4);
    assertInException("Incorrect", () -> f.convert_out(ins(1, 2, 3)));
  }

  @Test
  public void test_ad_weight_sp() {
    Assumptions.assumeTrue(hasNlpsol("ipopt"));
    MX x = MX.sym("x", 4);
    Function f = new Function("f", List.of(x), List.of(x));
    assertOutputs(List.of("1 forward sweeps needed"), List.of("1 reverse sweeps needed"),
        () -> ca.jacobian(f.call(List.of(ca.sin(x))).get(0), x, opts("helper_options", opts("verbose", true, "ad_weight_sp", 0))));
    assertOutputs(List.of("1 reverse sweeps needed"), List.of("1 forward sweeps needed"),
        () -> ca.jacobian(f.call(List.of(ca.sin(x))).get(0), x, opts("helper_options", opts("verbose", true, "ad_weight_sp", 1))));
    MX x2 = MX.sym("x", 4);
    Map<String, MX> nlp = new LinkedHashMap<>();
    nlp.put("x", x2);
    nlp.put("f", ca.dot(x2, x2));
    nlp.put("g", ca.sin(x2));
    assertOutputs(List.of("1 forward sweeps needed"), List.of("1 reverse sweeps needed"),
        () -> ca.nlpsol("solver", "ipopt", nlp, opts("common_options", opts("helper_options", opts("verbose", true, "ad_weight_sp", 0)))));
    assertOutputs(List.of("1 reverse sweeps needed"), List.of("1 forward sweeps needed"),
        () -> ca.nlpsol("solver", "ipopt", nlp, opts("common_options", opts("helper_options", opts("verbose", true, "ad_weight_sp", 1)))));
    Function f2 = new Function("f", List.of(x2), List.of(x2), opts("verbose", true, "ad_weight_sp", 1));
    assertOutputs(List.of("1 reverse sweeps needed"), List.of("1 forward sweeps needed"), () -> ca.jacobian(f2.call(List.of(ca.sin(x2))).get(0), x2));
    MX x3 = MX.sym("x", 100);
    Function f3 = new Function("f", List.of(x3), List.of(x3.get(0)), opts("verbose", true, "ad_weight_sp", 0));
    assertOutputs(List.of("2 forward sweeps needed"), List.of("1 reverse sweeps needed"), () -> ca.jacobian(f3.call(List.of(ca.sin(x3))).get(0), x3));
  }

  @Test
  public void test_map_exception() {
    MX x = MX.sym("x", 4);
    MX y = MX.sym("y", 4);
    Function f = new Function("f", List.of(x), List.of(x.plus(y)), opts("allow_free", true));
    String msg = CasadiMeta.compiler_flags().contains("CASADI_WITH_THREAD") ? "Evaluation failed" : "since variables [y] are free";
    assertInException(msg, () -> {
      Function F = f.map(4, "thread", 2);
      F.invoke(3);
    });
  }

  @Test
  public void test_DM_arg() {
    Function f = new Function("f", List.of(new DM(0, 1)), List.<DM>of());
    message(f.toString());
    f = new Function("f", List.of(new DM(0, 1), SX.sym("x")), List.of());
    message(f.toString());
    f = new Function("f", List.of(new DM(0, 1), MX.sym("x")), List.of());
    message(f.toString());
    Object v = ca.vertcat(List.of());
    assertTrue(v instanceof DM);
    assertTrue(Arrays.equals(new long[] {0, 1}, ((DM) v).size()));
    Object h = ca.horzcat(List.of());
    assertTrue(h instanceof DM);
    assertTrue(Arrays.equals(new long[] {1, 0}, ((DM) h).size()));
    Object vc = ca.veccat(List.of());
    assertTrue(vc instanceof DM);
    assertTrue(Arrays.equals(new long[] {0, 1}, ((DM) vc).size()));
    Object d = ca.diagcat(List.of());
    assertTrue(d instanceof DM);
    assertTrue(Arrays.equals(new long[] {0, 0}, ((DM) d).size()));
    // vcat/hcat/vvcat/dcat: python-only list aliases of vertcat/horzcat/veccat/diagcat
  }

  @Test
  public void test_nondiff() {
    for (String X : XS) {
      DM xn = dm(new double[][] {{1, 3}, {0.1, 2}});
      DM yn = dm(new double[][] {{4, 5}, {6, 7}});
      Map<String, Object> options = opts("is_diff_in", new boolean[] {true, false}, "is_diff_out", new boolean[] {true, false});
      Function f, F, G;
      if (X.equals("SX")) {
        SX x = SX.sym("x", 2, 2), y = SX.sym("y", 2, 2);
        f = new Function("f", List.of(x, y), List.of(ca.sin(x.plus(y)), x.times(y)), options);
        F = new Function("F", List.of(x, y), f.call(List.of(ca.cos(x), x.times(y).power(2.0))), List.of("x", "y"), List.of("z", "zz"), options);
        SX ye = new SX(2, 2);
        G = new Function("G", List.of(x, ye), List.of(f.call(List.of(ca.cos(x), x.times(yn).power(2.0))).get(0), new SX(f.call(List.of(ca.cos(xn), xn.times(yn).power(2.0))).get(1))),
            List.of("x", "y"), List.of("z", "zz"));
      } else {
        MX x = MX.sym("x", 2, 2), y = MX.sym("y", 2, 2);
        f = new Function("f", List.of(x, y), List.of(ca.sin(x.plus(y)), x.times(y)), options);
        Function fexpanded = f.expand();
        assertEquals(f.toString().split(" ")[0], fexpanded.toString().split(" ")[0]);
        F = new Function("F", List.of(x, y), f.call(List.of(ca.cos(x), x.times(y).power(2.0))), List.of("x", "y"), List.of("z", "zz"), options);
        MX ye = new MX(2, 2);
        G = new Function("G", List.of(x, ye), List.of(f.call(List.of(ca.cos(x), x.times(yn).power(2.0))).get(0), new MX(f.call(List.of(ca.cos(xn), xn.times(yn).power(2.0))).get(1))),
            List.of("x", "y"), List.of("z", "zz"));
      }
      for (Function ff : List.of(F.forward(1).forward(1), F.forward(1).reverse(1), F.reverse(1).forward(1), F.reverse(1).reverse(1))) {
        String s = ff.toString().replace("_y[2x2,0nz]", "foo");
        assertFalse(s.substring(Math.min(s.length(), "fwd1_adj1_F:(x[2x2],y[2x2]".length())).contains("y"), s);
      }
      checkfunction(F, G, List.of(xn, yn));
      List<java.util.function.UnaryOperator<Function>> ds = List.of(g -> g.forward(1), g -> g.reverse(1));
      for (java.util.function.UnaryOperator<Function> f1 : ds) {
        Function Gf = f1.apply(G), Ff = f1.apply(F);
        List<DM> arg = new ArrayList<>(List.of(xn, yn));
        for (long i = 2; i < Gf.n_in(); i++) arg.add(DM.rand(Gf.sparsity_in(i)));
        List<DM> a = Gf.call(arg), b = Ff.call(arg);
        for (int i = 0; i < a.size(); i++) checkarray(a.get(i), b.get(i));
      }
      for (java.util.function.UnaryOperator<Function> f1 : ds)
        for (java.util.function.UnaryOperator<Function> f2 : ds) {
          Function Gf = f1.apply(f2.apply(G)), Ff = f1.apply(f2.apply(F));
          List<DM> arg = new ArrayList<>(List.of(xn, yn));
          for (long i = 2; i < Gf.n_in(); i++) arg.add(DM.rand(Gf.sparsity_in(i)));
          List<DM> a = Gf.call(arg), b = Ff.call(arg);
          for (int i = 0; i < a.size(); i++) checkarray(a.get(i), b.get(i));
        }
    }
  }

  @Test
  public void test_non_diff_call() {
    MX x = MX.sym("x", 5);
    MX y = MX.sym("x", 5, 5);
    Function f = new Function("f", List.of(x, y), List.of(ca.mtimes(y, x), ca.mtimes(y.T(), x)),
        opts("never_inline", true, "is_diff_in", new boolean[] {false, true}, "is_diff_out", new boolean[] {false, true}));
    SX xs = SX.sym("x", 5);
    SX ys = SX.sym("x", 5, 5);
    Function g = new Function("g", List.of(xs, ys), f.call(List.of(xs, ys)));
    Function G = new Function("g", List.of(xs, ys), f.call(List.of(xs, ys)));
    DM.rng(1);
    DM x0 = DM.rand(5);
    DM y0 = DM.rand(5, 5);
    checkfunction(g, G, List.of(x0, y0));
  }

  @Test
  public void test_inline_linear_interpolant() {
    int N = 3, M = 4;
    double[][] dKnots = {new double[N], new double[M]};
    for (int i = 0; i < N; i++) dKnots[0][i] = i / (N - 1.0);
    for (int i = 0; i < M; i++) dKnots[1][i] = i / (M - 1.0);
    double[] dFlat = new double[2 * N * M];
    for (int i = 0; i < dFlat.length; i++) dFlat[i] = random.nextDouble();
    Function LUTparam = ca.interpolant("name", "linear", new long[] {N, M}, 2, opts("lookup_mode", List.of("exact", "linear")));
    Map<String, Object> onlyFd = opts("enable_fd", true, "enable_forward", false, "enable_reverse", false);
    Map<String, Object> options = opts("ad_weight_sp", -1);
    options.putAll(onlyFd);
    options.put("forward_options", onlyFd);
    options.put("reverse_options", onlyFd);
    message(options.toString());
    Function LUTparamRef = LUTparam.wrap_as_needed(options);
    Function Jref = LUTparamRef.jacobian();
    double[] cat = new double[N + M];
    System.arraycopy(dKnots[0], 0, cat, 0, N);
    System.arraycopy(dKnots[1], 0, cat, N, M);
    List<DM> inputs = new ArrayList<>(List.of(dm(0.2, 0.333), new DM(cat), new DM(dFlat)));
    LUTparam = ca.interpolant("name", "linear", new long[] {N, M}, 2, opts("inline", true, "lookup_mode", List.of("exact", "linear")));
    Function J = LUTparam.jacobian();
    checkarray(LUTparam.call(inputs).get(0), LUTparamRef.call(inputs).get(0));
    inputs.add(new DM(0.0));
    message("J_ref " + Jref + " " + Jref.call(inputs));
    checkfunction(J, Jref, inputs, 1e-8, 0.1);
  }

  @Test
  public void test_issue3079() {
    ca.interp1d(new double[] {1, 2, 3}, new MX(dm(4, 5, 6)), new double[] {1.1}, "linear");
  }

  // test_functionbuffer: python-only (memoryview-backed FunctionBuffer)

  // test_jit_buffer_eval: python-only (FunctionBuffer evaluation into a numpy memoryview)

  @Test
  public void test_codegen_inf_nan() {
    Assumptions.assumeTrue(RUN_SLOW);
    MX x = MX.sym("x");
    Function f = new Function("F", List.of(x), List.of(x.plus(Double.POSITIVE_INFINITY)));
    checkCodegen(f, ins(1), new CG().std("c99").main());
    x = MX.sym("x");
    f = new Function("F", List.of(x), List.of(x.plus(Double.NaN)));
    checkCodegen(f, ins(1), new CG().std("c99").main());
    x = MX.sym("x");
    f = new Function("F", List.of(x), List.of(x.plus(dm(Double.POSITIVE_INFINITY, Double.NaN, Double.NEGATIVE_INFINITY))));
    checkCodegen(f, ins(1), new CG().std("c99").main());
  }

  @Test
  public void test_codegen_with_mem() {
    Assumptions.assumeTrue(RUN_SLOW);
    MX x = MX.sym("x");
    Function f = new Function("F", List.of(x), List.of(x.times(3.0)));
    checkCodegen(f, ins(1), new CG().main().opt("with_mem", true).opt("force_canonical", true).def("inline=\"\""));
    checkCodegen(f, ins(1), new CG().main().opt("with_mem", true).opt("with_header", true).opt("force_canonical", true).def("inline=\"\""));
  }

  @Test
  public void test_codegen_scalars_bug() {
    Assumptions.assumeTrue(RUN_SLOW);
    MX x = MX.sym("x");
    MX z = x.times(3.0).div(ca.sin(x));
    Function f = new Function("F", List.of(x), List.of(z, z.div(x)), opts("live_variables", false));
    checkCodegen(f, ins(1), new CG().opt("codegen_scalars", true));
  }

  @Test
  public void test_bug_codegen_logical() {
    Assumptions.assumeTrue(RUN_SLOW);
    MX a = new MX(dm(1, 0, 0));
    MX b = new MX(dm(1, 1, 0));
    MX c = ca.logic_or(a, b);
    Function f = new Function("f", List.<MX>of(), List.of(c));
    checkCodegen(f, ins(), new CG());
  }

  private static List<java.util.function.Supplier<Function>> jitSerializeCases(String jitSerialize) {
    Map<String, Object> o = opts("jit", true, "compiler", "shell", "jit_options", opts("verbose", true), "verbose", false, "jit_serialize", jitSerialize);
    List<java.util.function.Supplier<Function>> r = new ArrayList<>();
    MX x = MX.sym("x");
    r.add(() -> new Function("f", List.of(x), List.of(x.minus(3.0).power(2.0)), o));
    MX x2 = MX.sym("x", 2);
    MX f = x2.get(0).times(x2.get(0)).plus(x2.get(1).times(x2.get(1)));
    if (ca.has_nlpsol("ipopt")) r.add(() -> ca.nlpsol("solver", "ipopt", Map.of("x", x2, "f", f), o));
    return r;
  }

  @Test
  public void test_jit_serialize() {
    Assumptions.assumeTrue(RUN_SLOW);
    Assumptions.assumeFalse(System.getProperty("os.name").toLowerCase().contains("mac"));
    Function[] fg = new Function[2];
    for (java.util.function.Supplier<Function> c : jitSerializeCases("source")) {
      assertOutputs(List.of("jit_tmp"), List.of(), () -> fg[0] = c.get());
      fg[0].save("f.casadi");
      assertOutputs(List.of("jit_tmp"), List.of(), () -> fg[1] = Function.load("f.casadi"));
      checkfunctionLight(fg[0], fg[1], Map.of());
      fg[1].save("f.casadi");
      assertOutputs(List.of("jit_tmp"), List.of(), () -> fg[1] = Function.load("f.casadi"));
      checkfunctionLight(fg[0], fg[1], Map.of());
    }
    for (java.util.function.Supplier<Function> c : jitSerializeCases("link")) {
      assertOutputs(List.of("jit_tmp"), List.of(), () -> fg[0] = c.get());
      fg[0].save("f.casadi");
      assertOutputs(List.of(), List.of("jit_tmp"), () -> fg[1] = Function.load("f.casadi"));
      checkfunctionLight(fg[0], fg[1], Map.of());
      fg[1].save("f.casadi");
      fg[1].delete();
      assertOutputs(List.of(), List.of("jit_tmp"), () -> fg[1] = Function.load("f.casadi"));
      checkfunctionLight(fg[0], fg[1], Map.of());
      fg[0].delete();
      fg[1].delete();
      System.gc();
      assertInException("No such file", () -> Function.load("f.casadi"));
    }
    for (java.util.function.Supplier<Function> c : jitSerializeCases("embed")) {
      assertOutputs(List.of("jit_tmp"), List.of(), () -> fg[0] = c.get());
      fg[0].save("f.casadi");
      assertOutputs(List.of(), List.of("jit_tmp"), () -> fg[1] = Function.load("f.casadi"));
      fg[1].invoke(Map.of());
      fg[0].delete();
      assertOutputs(List.of(), List.of("jit_tmp"), () -> fg[1] = Function.load("f.casadi"));
      fg[1].invoke(Map.of());
      fg[1].save("f.casadi");
      assertOutputs(List.of(), List.of("jit_tmp"), () -> fg[1] = Function.load("f.casadi"));
      fg[1].invoke(Map.of());
      fg[1].delete();
      assertOutputs(List.of(), List.of("jit_tmp"), () -> fg[1] = Function.load("f.casadi"));
    }
  }

  @Test
  public void test_map_get_function() {
    MX x = MX.sym("x");
    Function g = new Function("g", List.of(x), List.of(x.power(2.0)));
    Function ff = g.map(5);
    assertEquals(1, ff.get_function().size());
    Function f2 = ff.get_function("f");
    checkfunctionLight(g, f2, ins(3));
  }

  private static void assertNames(List<Function> ref, List<Function> res) {
    assertEquals(ref.size(), res.size());
    for (int i = 0; i < ref.size(); i++) assertEquals(ref.get(i).name(), res.get(i).name());
  }

  @Test
  public void test_find_function() {
    for (String X : List.of("MX", "SX")) {
      Function g, h, i, j, k, intg;
      Map<String, Object> ni = opts("never_inline", true);
      if (X.equals("MX")) {
        MX x = MX.sym("x");
        g = new Function("g", List.of(x), List.of(x.power(2.0)), ni);
        h = new Function("h", List.of(x), List.of(g.call(List.of(x.power(2.0))).get(0)), ni);
        i = new Function("i", List.of(x), List.of(h.call(List.of(x.power(2.0))).get(0).plus(g.call(List.of(x.power(2.0))).get(0))), ni);
        j = new Function("j", List.of(x), List.of(i.call(List.of(x.power(2.0))).get(0).plus(g.call(List.of(x.power(2.0))).get(0))), ni);
        k = new Function("k", List.of(x), List.of(j.call(List.of(x.power(2.0))).get(0)), ni);
        intg = ca.integrator("intg", "rk", Map.of("x", x, "ode", k.call(List.of(x)).get(0)));
      } else {
        SX x = SX.sym("x");
        g = new Function("g", List.of(x), List.of(x.power(2.0)), ni);
        h = new Function("h", List.of(x), List.of(g.call(List.of(x.power(2.0))).get(0)), ni);
        i = new Function("i", List.of(x), List.of(h.call(List.of(x.power(2.0))).get(0).plus(g.call(List.of(x.power(2.0))).get(0))), ni);
        j = new Function("j", List.of(x), List.of(i.call(List.of(x.power(2.0))).get(0).plus(g.call(List.of(x.power(2.0))).get(0))), ni);
        k = new Function("k", List.of(x), List.of(j.call(List.of(x.power(2.0))).get(0)), ni);
        intg = ca.integrator("intg", "rk", Map.of("x", x, "ode", k.call(List.of(x)).get(0)));
      }
      k.generate("f.c");
      assertNames(List.of(j), k.find_functions(0));
      assertNames(List.of(j, i, g), k.find_functions(1));
      assertNames(List.of(j, i, h, g), k.find_functions(2));
      assertNames(List.of(j, i, h, g), k.find_functions());
      List<Function> res = intg.find_functions();
      message(res.toString());
      assertNames(List.of(intg.get_function("dae"), k, j, i, h, g, intg.get_function("daeF"), intg.get_function("step")), res);
      Function K = k.map(5);
      assertNames(List.of(k, j, i, h, g), K.find_functions());
    }
  }

  @Test
  public void test_post_expand() {
    MX x = MX.sym("x");
    Function f = new Function("f", List.of(x), List.of(x.power(2.0)));
    message(f.toString());
    f = new Function("f", List.of(x), List.of(x.power(2.0)), opts("post_expand", true));
    message(f.toString());
    assertTrue(f.is_a("SXFunction"));
    Function g = new Function("f", List.of(x), List.of(x.power(2.0)), opts("post_expand", true, "post_expand_options", opts("print_in", true)));
    message(g.toString());
    assertTrue(g.is_a("SXFunction"));
    assertOutputs(List.of("Input 0 (i0): 3"), List.of(), () -> g.invoke(3));
  }

  private static Map<String, DM> qpInputs() {
    Map<String, DM> m = new LinkedHashMap<>();
    m.put("h", dm(new double[][] {{1, 0.2}, {0.2, 1}}));
    m.put("g", dm(1, 2));
    m.put("a", dm(new double[][] {{1, 1}}));
    m.put("lba", new DM(-1));
    m.put("uba", new DM(1));
    return m;
  }

  private static Map<String, Sparsity> qpStructure() {
    return Map.of("a", Sparsity.dense(1, 2), "h", Sparsity.dense(2, 2));
  }

  private static <T extends GenericMatrix> Map<String, T> withG(Map<String, DM> in, java.util.function.Function<DM, T> conv, String key, T sym) {
    Map<String, T> r = new LinkedHashMap<>();
    for (Map.Entry<String, DM> e : in.entrySet()) r.put(e.getKey(), conv.apply(e.getValue()));
    r.put(key, sym);
    return r;
  }

  private static Function wrapper(Function c, Map<String, DM> inputs, String X, String symbol, boolean names) {
    if (X.equals("MX")) {
      MX g = MX.sym(symbol, 2);
      MX r = c.call(withG(inputs, MX::new, symbol, g)).get("x");
      return names ? new Function("wrapper", List.of(g), List.of(r), List.of(symbol), List.of("res")) : new Function("wrapper", List.of(g), List.of(r));
    }
    SX g = SX.sym(symbol, 2);
    SX r = c.call(withG(inputs, SX::new, symbol, g)).get("x");
    return names ? new Function("wrapper", List.of(g), List.of(r), List.of(symbol), List.of("res")) : new Function("wrapper", List.of(g), List.of(r));
  }

  private static int countRefCounters(String file) {
    int n = 0;
    Pattern p = Pattern.compile("static int casadi_(\\w+)_ref_counter");
    for (String line : read(file).split("\n")) if (p.matcher(line).find()) n++;
    return n;
  }

  @Test
  public void test_memful_main() {
    Assumptions.assumeTrue(hasConic("osqp"));
    Function c = ca.conic("conic", "osqp", qpStructure(), opts("print_problem", true, "osqp.verbose", false));
    Map<String, DM> inputs = qpInputs();
    String[] extra = {"-Wno-endif-labels", "-Wno-unused-variable", "-Wno-strict-prototypes"};
    checkCodegen(c, c.convert_in(inputs), new CG().main().std("c99").flag(extra).lib("osqp"));
    c.generate("f.c");
    assertEquals(1, countRefCounters("f.c"));
    for (String X : List.of("MX", "SX")) {
      Function f = wrapper(c, inputs, X, "g", false);
      message(f.toString());
      checkCodegen(f, List.of(dm(1, 2)), new CG().main().std("c99").flag(extra).lib("osqp"));
      f.generate("f.c");
      assertEquals(1, countRefCounters("f.c"));
    }
  }

  @Test
  public void test_codegen_thread_safe() {
    Assumptions.assumeTrue(RUN_SLOW);
    for (int variant : new int[] {0, 1}) {
      Function c;
      Map<String, DM> inputs = new LinkedHashMap<>();
      List<String> extralibs;
      String symbol;
      if (variant == 0 && hasConic("osqp")) {
        c = ca.conic("conic", "osqp", qpStructure(), opts("osqp.verbose", false, "osqp.alpha", 1, "osqp.eps_abs", 1e-8, "osqp.eps_rel", 1e-8, "print_time", false));
        inputs = qpInputs();
        extralibs = List.of("osqp");
        symbol = "g";
      } else if (hasNlpsol("fatrop")) {
        MX x = MX.sym("x", 2);
        c = ca.nlpsol("solver", "fatrop", Map.of("x", x, "f", ca.sumsqr(x.minus(10.0)), "g", x.get(0).plus(x.get(1))), opts("fatrop.print_level", 0, "print_time", false));
        inputs.put("x0", dm(1, 2));
        inputs.put("lbx", dm(-5, -4));
        inputs.put("ubx", dm(5, 4));
        inputs.put("lbg", new DM(-2));
        inputs.put("ubg", new DM(2));
        extralibs = List.of("fatrop", "blasfeo", "m");
        symbol = "ubx";
      } else {
        continue;
      }
      String[] extra = {"-Wno-endif-labels", "-Wno-unused-variable", "-Wno-strict-prototypes"};
      List<Object[]> functions = new ArrayList<>();
      functions.add(new Object[] {c, inputs});
      for (String X : List.of("MX", "SX")) functions.add(new Object[] {wrapper(c, inputs, X, symbol, true), Map.of(symbol, inputs.get(symbol))});
      for (Object[] fl : functions) {
        Function f = (Function) fl[0];
        @SuppressWarnings("unchecked")
        Map<String, DM> local = (Map<String, DM>) fl[1];
        for (String mapType : List.of("serial", "thread", "openmp")) {
          Function F = f.map(2, mapType);
          Map<String, DM> wide = new LinkedHashMap<>();
          for (Map.Entry<String, DM> e : local.entrySet()) wide.put(e.getKey(), ca.horzcat(List.of(e.getValue(), e.getValue().times(1.2))));
          CG base = new CG().main().std("c99").flag(extra).noExternal().digits(12);
          for (String l : extralibs) base.lib(l);
          if (mapType.equals("openmp")) {
            CG o = base.flag("-fopenmp").opt("thread_safe", true).def("CASADI_MAX_NUM_THREADS=2").def("CASADI_THREAD_TYPE=CASADI_THREAD_TYPE_OMP");
            checkCodegen(F, F.convert_in(wide), o);
          } else if (mapType.equals("thread")) {
            for (String staticInit : List.of("0", "1")) {
              CG o = new CG().main().std("c99").flag(extra).flag("-pthread").noExternal().digits(12).opt("thread_safe", true).def("CASADI_MAX_NUM_THREADS=2")
                  .def("CASADI_THREAD_TYPE=CASADI_THREAD_TYPE_POSIX").def("CASADI_MUTEX_USE_STATIC_INIT=" + staticInit);
              for (String l : extralibs) o.lib(l);
              checkCodegen(F, F.convert_in(wide), o);
            }
            CG o = new CG().main().std("c11").flag(extra).noExternal().digits(12).opt("thread_safe", true).def("CASADI_MAX_NUM_THREADS=2").def("CASADI_THREAD_TYPE=CASADI_THREAD_TYPE_C11");
            for (String l : extralibs) o.lib(l);
            checkCodegen(F, F.convert_in(wide), o);
          } else {
            checkCodegen(f, f.convert_in(local), base);
            Map<String, DM> scaled = new LinkedHashMap<>();
            for (Map.Entry<String, DM> e : local.entrySet()) scaled.put(e.getKey(), e.getValue().times(1.2));
            checkCodegen(f, f.convert_in(scaled), base);
            checkCodegen(F, f.convert_in(wide), base);
          }
        }
      }
    }
  }

  @Test
  public void test_memful_external() {
    Assumptions.assumeTrue(hasConic("osqp"));
    Assumptions.assumeTrue(RUN_SLOW);
    Function c = ca.conic("conic", "osqp", qpStructure(), opts("print_problem", true, "osqp.verbose", false));
    Map<String, DM> inputs = qpInputs();
    String[] extra = {"-Wno-endif-labels", "-Wno-unused-variable", "-Wno-strict-prototypes"};
    Map<String, Object> cg = checkCodegen(c, c.convert_in(inputs), new CG().std("c99").flag(extra).lib("osqp").debug());
    Function F = ((Function) cg.get("F")).wrap();
    message("memful_external");
    checkCodegen(F, c.convert_in(inputs), new CG().main().lib((String) cg.get("libname")).debug());
    for (String X : List.of("MX", "SX")) {
      Function f = wrapper(c, inputs, X, "g", false);
      message(f.toString());
      cg = checkCodegen(f, List.of(dm(1, 2)), new CG().std("c99").flag(extra).lib("osqp").debug());
      F = ((Function) cg.get("F")).wrap();
      message("memful_external");
      checkCodegen(F, List.of(dm(1, 2)), new CG().main().lib((String) cg.get("libname")).debug());
    }
  }

  @Test
  public void test_memful_jit() {
    Assumptions.assumeTrue(hasConic("osqp"));
    Assumptions.assumeTrue(RUN_SLOW);
    Function c = ca.conic("conic", "osqp", qpStructure(), opts("print_problem", true, "osqp.verbose", false));
    Map<String, DM> inputs = qpInputs();
    List<String> extra = new ArrayList<>(List.of("-Wno-endif-labels", "-Wno-unused-variable", "-Wno-strict-prototypes"));
    for (String X : List.of("MX", "SX")) {
      String libdir = GlobalOptions.getCasadiPath(), includedir = includeDir();
      List<String> compilerFlags = extra;
      List<String> linkerFlags = new ArrayList<>();
      compilerFlags.add("-I" + includedir);
      compilerFlags.add("-L" + libdir);
      compilerFlags.add("-losqp");
      linkerFlags.add("-I" + includedir);
      linkerFlags.add("-L" + libdir);
      linkerFlags.add("-losqp");
      Map<String, Object> jitOpts = opts("jit", true, "compiler", "shell", "jit_options", opts("compiler_flags", compilerFlags, "linker_flags", linkerFlags, "verbose", true));
      Function f = wrapper(c, inputs, X, "g", false);
      Function fj;
      if (X.equals("MX")) {
        MX g = MX.sym("g", 2);
        fj = new Function("wrapper", List.of(g), List.of(c.call(withG(inputs, MX::new, "g", g)).get("x")), jitOpts);
      } else {
        SX g = SX.sym("g", 2);
        fj = new Function("wrapper", List.of(g), List.of(c.call(withG(inputs, SX::new, "g", g)).get("x")), jitOpts);
      }
      checkfunctionLight(f, fj, List.of(dm(1, 2)));
      CG o = new CG().std("c99").externalOpts(jitOpts).lib("osqp").debug();
      o.extraOptions.addAll(extra);
      checkCodegen(f, List.of(dm(1, 2)), o);
    }
  }

  @Test
  public void test_cse() {
    for (String X : XS) {
      boolean wZero, resZero;
      Function f1, f2;
      if (X.equals("SX")) {
        SX x = SX.sym("x", 2), y = SX.sym("y", 2), p = SX.sym("p", 2);
        SX z = p.times(ca.cos(x.plus(y.times(5.0))));
        SX z2 = p.times(ca.cos(x.plus(y.times(5.0))));
        SX w = z.minus(z2);
        wZero = w.is_zero();
        resZero = ca.cse(w).is_zero();
        f1 = new Function("f", List.of(x, y, p), List.of(w), opts("cse", false));
        f2 = new Function("f", List.of(x, y, p), List.of(w), opts("cse", true));
      } else {
        MX x = MX.sym("x", 2), y = MX.sym("y", 2), p = MX.sym("p", 2);
        MX z = p.times(ca.cos(x.plus(y.times(5.0))));
        MX z2 = p.times(ca.cos(x.plus(y.times(5.0))));
        MX w = z.minus(z2);
        wZero = w.is_zero();
        resZero = ca.cse(w).is_zero();
        f1 = new Function("f", List.of(x, y, p), List.of(w), opts("cse", false));
        f2 = new Function("f", List.of(x, y, p), List.of(w), opts("cse", true));
      }
      assertFalse(wZero);
      assertTrue(resZero);
      assertTrue(f1.n_instructions() > 3);
      assertTrue(f2.n_instructions() <= 3);
    }
  }

  @Test
  public void test_cse_call() {
    for (String X : XS) {
      String y;
      Map<String, Object> ni = opts("never_inline", true);
      if (X.equals("SX")) {
        SX x = SX.sym("x");
        Function f = new Function("f", List.of(x), List.of(x.power(2.0)), List.of("x"), List.of("y"), ni);
        Function fcopy = new Function("f", List.of(x), List.of(x.power(2.0)), List.of("x"), List.of("y"), ni);
        y = ca.cse(ca.vertcat(List.of(f.call(List.of(ca.sin(x))).get(0), fcopy.call(List.of(ca.sin(x))).get(0), f.call(List.of(ca.sin(x))).get(0)))).toString();
      } else {
        MX x = MX.sym("x");
        Function f = new Function("f", List.of(x), List.of(x.power(2.0)), List.of("x"), List.of("y"), ni);
        Function fcopy = new Function("f", List.of(x), List.of(x.power(2.0)), List.of("x"), List.of("y"), ni);
        y = ca.cse(ca.vertcat(List.of(f.call(List.of(ca.sin(x))).get(0), fcopy.call(List.of(ca.sin(x))).get(0), f.call(List.of(ca.sin(x))).get(0)))).toString();
      }
      message(y);
      assertTrue(y.contains(X.equals("MX") ? "vertcat(@1, @1, @1)" : "[@1, @1, @1]"), y);
    }
  }

  private static List<String> prefixed(String prefix, List<String> names) {
    List<String> r = new ArrayList<>();
    for (String n : names) r.add(prefix + n);
    return r;
  }

  private static void checkZeroPattern(Function f, Function fRef, List<DM> inputs) {
    List<DM> res = f.call(inputs, false, false), resRef = fRef.call(inputs, false, false);
    for (int k = 0; k < res.size(); k++) {
      if (resRef.get(k).is_zero()) assertTrue(res.get(k).is_zero());
      else assertFalse(res.get(k).is_zero());
    }
  }

  @Test
  public void test_stop_diff() {
    MX x = MX.sym("x");
    MX y = ca.sin(x);
    int maxOrder = 2;
    for (int order = 1; order <= maxOrder; order++) {
      MX testZ = x.power((double) (order - 1)).times(3.0);
      Function fa = new Function("f", List.of(x), List.of(ca.stop_diff(testZ, order)), List.of("x"), List.of("z"));
      Function fb = new Function("f", List.of(x), List.of(testZ), List.of("x"), List.of("z"));
      checkfunction(fa, fb, List.of(DM.rand(1, 1)));
      MX z = ca.stop_diff(y, order);
      DM.rng(1);
      DM v0 = DM.rand(1);
      DM expr0 = DM.rand(1);
      List<java.util.function.BinaryOperator<MX>> ops = List.of(ca::gradient, ca::jacobian, (e, v) -> ca.jtimes(e, v, v0), (e, v) -> ca.jtimes(e, v, expr0, true));
      for (java.util.function.BinaryOperator<MX> op : ops) {
        MX yd = y, zd = z;
        for (int i = 0; i < maxOrder + 2; i++) {
          Function f = new Function("f", List.of(x), List.of(yd, zd));
          DM x0 = DM.rand(f.sparsity_in(0));
          List<DM> res = f.call(List.of(x0));
          if (i >= order) assertTrue(res.get(1).is_zero());
          else assertEquals(res.get(0).to_double(), res.get(1).to_double());
          if (i < order) {
            yd = op.apply(yd, x);
            zd = op.apply(zd, x);
          }
        }
      }
      Function f = new Function("f", List.of(x), List.of(z), List.of("x"), List.of("z"));
      Function fNormal = new Function("f", List.of(x), List.of(y), List.of("x"), List.of("z"));
      Function fRef = new Function("f", List.of(x), List.of(testZ), List.of("x"), List.of("z"));
      for (int i = 0; i < 3; i++) {
        List<DM> inputs = randInputs(f);
        if (i >= order) checkZeroPattern(f, fRef, inputs);
        else checkfunctionLight(f, fNormal, inputs);
        if (i < order) {
          List<String> inLabels = new ArrayList<>(f.name_in());
          inLabels.addAll(prefixed("fwd:", f.name_in()));
          inLabels.addAll(prefixed("adj:", f.name_out()));
          List<String> outLabels = new ArrayList<>(prefixed("fwd:", f.name_out()));
          outLabels.addAll(prefixed("adj:", f.name_in()));
          for (String eIn : f.name_in())
            for (String eOut : f.name_out()) {
              outLabels.add("jac:" + eOut + ":" + eIn);
              outLabels.add("grad:" + eOut + ":" + eIn);
            }
          f = f.factory("f", inLabels, outLabels);
          fNormal = fNormal.factory("f", inLabels, outLabels);
          fRef = fRef.factory("f", inLabels, outLabels);
        }
      }
    }
  }

  @Test
  public void test_stop_diff_cross() {
    MX x = MX.sym("x", 2);
    MX y = ca.sin(x.get(0)).plus(ca.exp(x.get(1))).plus(ca.cos(x.get(0).power(2.0).times(x.get(1).power(2.0))));
    int maxOrder = 2;
    for (int order = 1; order <= maxOrder; order++) {
      MX testZ = x.power((double) (order - 1)).times(3.0);
      Function fa = new Function("f", List.of(x), List.of(ca.stop_diff(testZ, order)), List.of("x"), List.of("z"));
      Function fb = new Function("f", List.of(x), List.of(testZ), List.of("x"), List.of("z"));
      checkfunction(fa, fb, List.of(DM.rand(1, 1)));
      MX z = ca.stop_diff(y, order);
      DM.rng(1);
      List<java.util.function.BinaryOperator<MX>> ops = List.of(ca::jacobian, (e, v) -> ca.jtimes(e, v, DM.rand(v.sparsity())),
          (e, v) -> ca.jtimes(e, v, DM.rand(e.sparsity()), true));
      for (java.util.function.BinaryOperator<MX> op : ops) {
        MX yd = y, zd = z;
        for (int i = 0; i < maxOrder + 1; i++) {
          Function f = new Function("f", List.of(x), List.of(yd, zd));
          DM x0 = DM.rand(f.sparsity_in(0));
          List<DM> res = f.call(List.of(x0));
          if (i >= order) assertTrue(res.get(1).is_zero());
          else checkarray(res.get(0), res.get(1));
          if (i < order) {
            DM.rng(5 + i);
            yd = op.apply(yd, x);
            DM.rng(5 + i);
            zd = op.apply(zd, x);
          }
        }
      }
      Function f = new Function("f", List.of(x), List.of(z), List.of("x"), List.of("z"));
      Function fNormal = new Function("f", List.of(x), List.of(y), List.of("x"), List.of("z"));
      Function fRef = new Function("f", List.of(x), List.of(testZ), List.of("x"), List.of("z"));
      if (order <= 1) {
        for (int i = 0; i < 3; i++) {
          List<DM> inputs = randInputs(f);
          if (i >= order) checkZeroPattern(f, fRef, inputs);
          else checkfunctionLight(f, fNormal, inputs);
          if (i < order) {
            List<String> inLabels = new ArrayList<>(f.name_in());
            inLabels.addAll(prefixed("fwd:", f.name_in()));
            inLabels.addAll(prefixed("adj:", f.name_out()));
            List<String> outLabels = new ArrayList<>(prefixed("fwd:", f.name_out()));
            outLabels.addAll(prefixed("adj:", f.name_in()));
            for (String eIn : f.name_in())
              for (String eOut : f.name_out()) outLabels.add("jac:" + eOut + ":" + eIn);
            f = f.factory("f", inLabels, outLabels);
            fNormal = fNormal.factory("f", inLabels, outLabels);
            fRef = fRef.factory("f", inLabels, outLabels);
          }
        }
      }
    }
  }

  @Test
  public void test_stop_diff_p() {
    MX x = MX.sym("x");
    MX y = MX.sym("y");
    MX yy = ca.sin(y.times(x));
    int maxOrder = 2;
    for (int order = 1; order <= maxOrder; order++) {
      DM.rng(1);
      message("order " + order);
      MX testZ = ca.sin(y).times(x.power((double) (order - 1)));
      MX z = ca.stop_diff(yy, x, order);
      DM v0 = DM.rand(1);
      DM expr0 = DM.rand(1);
      List<java.util.function.BinaryOperator<MX>> ops = List.of(ca::gradient, ca::jacobian, (e, v) -> ca.jtimes(e, v, v0), (e, v) -> ca.jtimes(e, v, expr0, true));
      for (java.util.function.BinaryOperator<MX> op : ops) {
        MX yd = yy, zd = z;
        for (int i = 0; i < maxOrder + 2; i++) {
          Function f = new Function("f", List.of(x, y), List.of(yd, zd));
          DM x0 = DM.rand(f.sparsity_in(0));
          DM p0 = DM.rand(f.sparsity_in(1));
          List<DM> res = f.call(List.of(x0, p0));
          if (i >= order) assertTrue(res.get(1).is_zero());
          else assertEquals(res.get(0).to_double(), res.get(1).to_double());
          if (i < order) {
            yd = op.apply(yd, x);
            zd = op.apply(zd, x);
          }
        }
      }
      Function f = new Function("f", List.of(x, y), List.of(z), List.of("x", "y"), List.of("z"));
      Function fNormal = new Function("f", List.of(x, y), List.of(yy), List.of("x", "y"), List.of("z"));
      Function fRef = new Function("f", List.of(x, y), List.of(testZ), List.of("x", "y"), List.of("z"));
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
          List<String> inLabels = new ArrayList<>(f.name_in());
          inLabels.addAll(prefixed("adj:", f.name_out()));
          List<String> outLabels = new ArrayList<>(f.name_out());
          outLabels.addAll(prefixed("adj:", f.name_in()));
          f = f.factory("f", inLabels, outLabels);
          fNormal = fNormal.factory("f", inLabels, outLabels);
          fRef = fRef.factory("f", inLabels, outLabels);
        }
      }
    }
  }

  @Test
  public void test_no_hess2() {
    Assumptions.assumeTrue(hasConic("qpoases"));
    MX y = MX.sym("y", 2);
    Function f = new Function("foo", List.of(y), List.of(ca.vertcat(List.of(ca.sin(y.get(0).times(y.get(1))), ca.cos(y.get(0).times(y.get(1)))))),
        opts("never_inline", true, "print_out", true));
    Opti opti = new Opti();
    MX x = opti.variable(2);
    opti.subject_to(ca.ge(ca.no_hess(f.call(List.of(x.power(2.0))).get(0)), 0.0));
    opti.minimize(ca.sumsqr(x.minus(2.0)));
    opti.solver("sqpmethod");
    Function F = opti.to_function("F", List.of(x), List.of(x));
    F.generate("F.c");
    String code = read("F.c");
    Matcher m = Pattern.compile("\\w+FS").matcher(code);
    while (m.find()) assertTrue(m.group().split("_").length <= 2, m.group());
  }

  @Test
  public void test_issue_1522() {
    MX V2 = MX.sym("X", 8);
    List<MX> xt = ca.vertsplit(V2, new long[] {0, 6, 8});
    List<MX> xs = ca.vertsplit(xt.get(0), new long[] {0, 3, 6});
    List<MX> travels = ca.vertsplit(xt.get(1), new long[] {0, 1, 2});
    MX dist = new MX(0.0);
    for (int j = 0; j < 2; j++) dist = dist.plus(ca.sum1(xs.get(0).minus(xs.get(j).plus(travels.get(j))).power(2.0)));
    MX p = new MX(0, 1);
    Function f = new Function("f", List.of(V2, p), List.of(dist.unaryMinus()));
    for (boolean compact : new boolean[] {true, false})
      for (Function F : List.of(f, f.expand())) {
        Sparsity HF = F.reverse(1).jac_sparsity(0, 0, compact);
        assertTrue(HF.is_symmetric());
      }
  }

  @Test
  public void test_issue_3074() {
    MX x = MX.sym("x", 4);
    Function f = new Function("f", List.<GenericMatrix>of(x), List.<GenericMatrix>of(ca.dot(1.0, 1.0)));
    message(f.toString());
    Function fr = f.reverse(1);
    message(fr.toString());
    x = MX.sym("x", ca.sparsify(dm(0, 1, 1, 0)).sparsity());
    MX y = ca.dot(x, x);
    f = new Function("f", List.of(x), List.of(ca.dot(y, y)));
    message(f.toString());
    fr = f.reverse(1);
    message(fr.toString());
    message(fr.jac_sparsity().toString());
    f = new StringDeserializer(ISSUE_3074_DATA).unpack_function();
    message(f.reverse(1).toString());
    f.reverse(1).jac_sparsity();
  }

  private static final String ISSUE_3074_DATA = "jhpnnagiieahaaaadaaaaaaaaaaaaaaaaafaegaakaaaaaaaneifgefhogdgehjgpgogbaaaaaaapaaaaaaaegfgmgbgjhpfbgchhgfhngfgogehdhaaaaaacaaaaaaahaaaaaaaaaaaaaaabababababababaaaaaaaaaaaaaaaaahaaaaaaaaaaaaaaaegfaaaaaaaaaaaaaaabaaaaaaaaaaaaaaabaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaabaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaeggaaaaaaaaaaaaaaacaaaaaaaaaaaaaaabaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaacaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaabaaaaaaaaaaaaaaachbaaaaaaaaaaaaaaaegkaaaaaaaaaaaaaaagaaaaaaaaaaaaaaabaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaagaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaabaaaaaaaaaaaaaaacaaaaaaaaaaaaaaadaaaaaaaaaaaaaaaeaaaaaaaaaaaaaaafaaaaaaaaaaaaaaaegeaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaabaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaachdaaaaaaaaaaaaaaaegiaaaaaaaaaaaaaaaeaaaaaaaaaaaaaaabaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaeaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaabaaaaaaaaaaaaaaacaaaaaaaaaaaaaaadaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaahaaaaaaaaaaaaaaacaaaaaaajgadcaaaaaaajgbdcaaaaaaajgcdcaaaaaaajgddcaaaaaaajgedcaaaaaaajgfdcaaaaaaajggdaaaaaaaaaaaaaaaaaabagaaaaaaadhpgfhchdgfgbahaaaaaaakgjgehpfehngahaaaaaaaaaaaaaaaafaaaaaaadgmgbgoghgaaegbaaaaaaaaaaaaaaaaebababaaabababaaapbfilobfilobfnpdmfpicmfpicmfpnpdaaaaaeaaaaaaaaaaaaaaaabakdmiadcooijhfeodaaaaaaaaaaaaaaaabaaaaaaaocdaaaaaaangehihaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaachfaaaaaaaaaaaaaaahaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaabaaaaaaahaaaaaaaaaaaaaaaegpcaaaaaaaaaaaaaaaaaaaaaachaaaaaaaaaaaaaaaaeaaaaaaaehjgngfgegndaaaaaacaaaaaaaaaaaaaaaegpcaaaaaaaaaaaaaaaaaaaaaachaaaaaaaaaaaaaaaagaaaaaaacgpgegjhocihegpcaaaaaaaaaaaaaaaaaaaaaachaaaaaaaaaaaaaaaaiaaaaaaacgpgegjhocghpfihchbaaaaaaaaaaaaaaaegndaaaaaacaaaaaaaaaaaaaaaegpcaaaaaaaaaaaaaaaaaaaaaachaaaaaaaaaaaaaaaalaaaaaaaegfgchiccgpgegjhocihjcegpcaaaaaaaaaaaaaaaaaaaaaachaaaaaaaaaaaaaaaanaaaaaaaegfgchiccgpgegjhocghpfihjcchbaaaaaaaaaaaaaaaegndaaaaaagaaaaaaaaaaaaaaaegpcaaaaaaaaaaaaaaaaaaaaaachaaaaaaaaaaaaaaaaiaaaaaaacgpgegjhocggpfihegpcaaaaaaaaaaaaaaaaaaaaaachaaaaaaaaaaaaaaaaiaaaaaaacgpgegjhocbgpfihegpcaaaaaaaaaaaaaaaaaaaaaachaaaaaaaaaaaaaaaajaaaaaaabgdgdgfgmgocbgpfihegpcaaaaaaaaaaaaaaaaaaaaaachaaaaaaaaaaaaaaaakaaaaaaabgdgdgfgmgocngbgpfihegpcaaaaaaaaaaaaaaaaaaaaaachaaaaaaaaaaaaaaaalaaaaaaabgdgdgfgmgoccgpfihocfhegpcaaaaaaaaaaaaaaaaaaaaaachaaaaaaaaaaaaaaaalaaaaaaabgdgdgfgmgoccgpfihocjhchcaaaaaaaaaaaaaaaegmcaaaaaaadaaaaaaaaaaaaaaaachdaaaaaaaaaaaaaaaegmcaaaaaaadaaaaaaaaaaaaaaaachdaaaaaaaaaaaaaaaegndaaaaaaeaaaaaaaaaaaaaaaegpcaaaaaaaaaaaaaaaaaaaaaachaaaaaaaaaaaaaaaagaaaaaaacgpgegjhochgegpcaaaaaaaaaaaaaaaaaaaaaachaaaaaaaaaaaaaaaagaaaaaaacgpgegjhocdgegpcaaaaaaaaaaaaaaaaaaaaaachaaaaaaaaaaaaaaaagaaaaaaacgpgegjhocngegpcaaaaaaaaaaaaaaaaaaaaaachaaaaaaaaaaaaaaaalaaaaaaabgdgdgfgmgoccgpfihoccgcheaaaaaaaaaaaaaaabaaaaaaaaaaaaaaaaaaaaaaabaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaahaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaabaaaaaaaaaaaaaaaaa";

  @Test
  public void test_issue_4345() {
    for (String X : XS) {
      List<Function> fs = new ArrayList<>();
      if (X.equals("SX")) {
        SX x = SX.sym("x", Sparsity.upper(3)), y = SX.sym("y", Sparsity.upper(3));
        for (SX z : List.of(x.power(2.0), x.power(2.0).plus(y.get(0, 0)), x.times(y))) fs.add(new Function("f", List.of(x, y), List.of(z)));
      } else {
        MX x = MX.sym("x", Sparsity.upper(3)), y = MX.sym("y", Sparsity.upper(3));
        for (MX z : List.of(x.power(2.0), x.power(2.0).plus(y.get(0, 0)), x.times(y))) fs.add(new Function("f", List.of(x, y), List.of(z)));
      }
      long[] adjYnnz = {0, 6, 6};
      for (int k = 0; k < 3; k++) {
        Function adj = fs.get(k).reverse(1);
        assertEquals(6, adj.nnz_out(0));
        assertEquals(adjYnnz[k], adj.nnz_out(1));
      }
    }
  }

  @Test
  public void test_issue_3134() {
    Assumptions.assumeTrue(hasNlpsol("ipopt"));
    MX p = MX.sym("p");
    MX v = MX.sym("v");
    MX u = MX.sym("u");
    MX x0 = ca.vertcat(List.of(p, v));
    MX rhs = ca.vertcat(List.of(v, u));
    MX t0 = MX.sym("t0");
    MX DT = MX.sym("DT");
    Function F = new Function("F", List.of(x0, u, t0, DT), List.of(x0.plus(DT.times(rhs))), List.of("x0", "u", "t0", "DT"), List.of("xf"));
    List<MX> x = new ArrayList<>(), g = new ArrayList<>();
    MX Xk = MX.sym("Xk", 2);
    x.add(Xk);
    MX T = MX.sym("T");
    x.add(T);
    g.add(T);
    MX t0b = MX.sym("t0");
    x.add(t0b);
    List<MX> DTs = List.of(T.div(2.0), t0b.plus(T).minus(t0b.plus(T.div(2.0))));
    for (int i = 0; i < 2; i++) {
      MX Uk = MX.sym("Uk");
      x.add(Uk);
      MX XkNext = MX.sym("Xk_next", 2);
      x.add(XkNext);
      MX Xf = F.call(Map.of("x0", Xk, "u", Uk, "t0", t0b, "DT", DTs.get(i))).get("xf");
      g.add(XkNext.minus(Xf));
      Xk = XkNext;
    }
    Map<String, MX> nlp = Map.of("x", ca.veccat(x), "f", T, "g", ca.veccat(g));
    ca.nlpsol("solver", "ipopt", nlp);
  }

  private static void whichDependsChecks(Function f, boolean exempt) {
    List<String> z = List.of("z"), z2 = List.of("z2"), z3 = List.of("z3"), z4 = List.of("z4"), z5 = List.of("z5"), z7 = List.of("z7"), z8 = List.of("z8");
    boolean T = true, F = false;
    assertTrue(Arrays.equals(new boolean[] {T, T}, f.which_depends("x", z)));
    assertTrue(Arrays.equals(new boolean[] {F, F}, f.which_depends("x", z, 2)));
    if (!exempt) assertTrue(Arrays.equals(new boolean[] {F, F}, f.which_depends("x", z2, 2)));
    assertTrue(Arrays.equals(new boolean[] {T, T}, f.which_depends("x", z3, 2)));
    assertTrue(Arrays.equals(new boolean[] {T, T}, f.which_depends("x", z4)));
    assertTrue(Arrays.equals(new boolean[] {T, T}, f.which_depends("x", z4, 2)));
    assertTrue(Arrays.equals(new boolean[] {F, F}, f.which_depends("x", z5)));
    assertTrue(Arrays.equals(new boolean[] {F, F}, f.which_depends("x", z5, 2)));
    assertTrue(Arrays.equals(new boolean[] {F, F}, f.which_depends("x", List.of("z5", "z6"))));
    if (!exempt) {
      assertTrue(Arrays.equals(new boolean[] {T, T, T, T, T, T, T, F, F, F}, f.which_depends("x", List.of("z", "z2", "z3", "z4", "z5", "z6"), 1, true)));
      assertTrue(Arrays.equals(new boolean[] {F, F, F, F, T, T, T, F, F, F}, f.which_depends("x", List.of("z", "z2", "z3", "z4", "z5", "z6"), 2, true)));
    }
    assertTrue(Arrays.equals(new boolean[] {T, T, T, T, T, F, F, F}, f.which_depends("x", List.of("z", "z3", "z4", "z5", "z6"), 1, true)));
    assertTrue(Arrays.equals(new boolean[] {F, F, T, T, T, F, F, F}, f.which_depends("x", List.of("z", "z3", "z4", "z5", "z6"), 2, true)));
    assertTrue(Arrays.equals(new boolean[] {T, F}, f.which_depends("x", z7)));
    assertTrue(Arrays.equals(new boolean[] {F, T}, f.which_depends("y", z7)));
    assertTrue(Arrays.equals(new boolean[] {F, F}, f.which_depends("x", z7, 2)));
    assertTrue(Arrays.equals(new boolean[] {F, F}, f.which_depends("y", z7, 2)));
    assertTrue(Arrays.equals(new boolean[] {T, F}, f.which_depends("x", z8, 2)));
    assertTrue(Arrays.equals(new boolean[] {F, T}, f.which_depends("y", z8, 2)));
    assertTrue(Arrays.equals(new boolean[] {T, F, T, F}, f.which_depends("x", List.of("z7", "z8"), 1, true)));
    assertTrue(Arrays.equals(new boolean[] {F, F, T, F}, f.which_depends("x", List.of("z7", "z8"), 2, true)));
    assertTrue(Arrays.equals(new boolean[] {F, T, F, T}, f.which_depends("y", List.of("z7", "z8"), 1, true)));
    assertTrue(Arrays.equals(new boolean[] {F, F, F, T}, f.which_depends("y", List.of("z7", "z8"), 2, true)));
  }

  @Test
  public void test_which_depends() {
    SX x = SX.sym("x", 2);
    SX y = SX.sym("y", 2);
    SX z = x.plus(y);
    SX z2 = x.times(y);
    SX z3 = ca.sin(ca.dot(x, y));
    SX z4 = ca.cos(x);
    SX z5 = ca.sin(y);
    SX z6 = new SX(12.0);
    SX z7 = ca.vertcat(List.of(x.get(0), y.get(1)));
    SX z8 = ca.vertcat(List.of(x.get(0).power(2.0), y.get(1).power(2.0)));
    Function f = new Function("f", List.of(x, y), List.of(z, z2, z3, z4, z5, z6, z7, z8), List.of("x", "y"), List.of("z", "z2", "z3", "z4", "z5", "z6", "z7", "z8"));
    whichDependsChecks(f, false);
    if (!RUN_SLOW) return;
    Map<String, Object> cg = checkCodegen(f, randInputs(f), new CG().jacSparsity().forward());
    whichDependsChecks((Function) cg.get("F"), true);
  }

  @Test
  public void test_cat_input() {
    MX x = MX.sym("x", 2);
    Function fmx = new Function("f", List.of(ca.vertcat(List.of(new MX(2, 1), x))), List.of(ca.sin(x)), opts("always_inline", true));
    message(fmx.toString());
    MX X = MX.sym("X", 4);
    message(fmx.call(List.of(X)).get(0).toString());
    Function fEval = new Function("f_eval", List.of(X), List.of(fmx.call(List.of(X)).get(0)));
    message(fEval.invoke(dm(Double.NaN, Double.NaN, 3, Double.NaN)).toString());
  }

  @Test
  public void test_reshape_input() {
    MX x = MX.sym("x", 6, 2);
    SX xsx = SX.sym("x", 6, 2);
    Function fmx = new Function("f", List.of(ca.reshape(x, 3, 4)), List.of(ca.sin(x)), opts("always_inline", true));
    Function fsx = new Function("f", List.of(ca.reshape(xsx, 3, 4)), List.of(ca.sin(xsx)));
    DM.rng(1);
    for (DM input : List.of(DM.rand(3, 4))) checkfunctionLight(fmx, fsx, List.of(input));
  }

  @Test
  public void test_sparsity_cast_input() {
    MX X0 = MX.sym("X", Sparsity.lower(3));
    new Function("f", List.of(X0), List.of(X0.power(2.0)), opts("always_inline", true));
    Sparsity sp = Sparsity.lower(3);
    long N = sp.nnz();
    MX x = MX.sym("x", N);
    SX xsx = SX.sym("x", N);
    Function fmx = new Function("f", List.of(ca.sparsity_cast(x, sp)), List.of(ca.sin(x)), opts("always_inline", true));
    Function fsx = new Function("f", List.of(ca.sparsity_cast(xsx, sp)), List.of(ca.sin(xsx)));
    for (DM input : List.of(sparseDM(new double[][] {{1, 0, 0}, {3, 7, 0}, {1, 8, 9}}), sparseDM(new double[][] {{1, 0, 0}, {3, 0, 0}, {0, 8, 9}}),
        dm(new double[][] {{1, 1.1, 8}, {3, 1.2, 3}, {1, 8, 9}})))
      checkfunctionLight(fmx, fsx, List.of(input));
    for (MX X : List.of(MX.sym("X", Sparsity.lower(3)), MX.sym("X", Sparsity.diag(3)), MX.sym("X", 3, 3))) {
      Function fcmx = new Function("f", List.of(X), List.of(fmx.call(List.of(X)).get(0)));
      Function fcsx = new Function("f", List.of(X), List.of(fsx.call(List.of(X)).get(0)));
      checkfunctionLight(fcmx, fcsx, List.of(DM.rand(X.sparsity())));
    }
  }

  @Test
  public void test_external() {
    Assumptions.assumeTrue(RUN_SLOW);
    String src = "../python/assets/externalfun1.c";
    assertInException("config failed", () -> compileExternal("F", src, opts(), new CG()));
    for (int nArgs : new int[] {3, 7}) {
      List<String> cfg = new ArrayList<>();
      for (int i = 0; i < nArgs; i++) cfg.add("foo");
      Function externalfun1 = compileExternal("F", src, opts("config_args", cfg), new CG().debug());
      String libname = externalLibname;
      message(libname);
      List<DM> r = externalfun1.call(ins(0));
      message(externalfun1.stats().toString());
      assertEquals(nArgs + 1, (int) r.get(0).to_double());
      assertEquals(0, (int) r.get(1).to_double());
      MX x = MX.sym("x");
      Function g = new Function("g", List.of(x), List.of(externalfun1.call(List.of(ca.sin(x))).get(0).plus(externalfun1.call(List.of(ca.cos(x))).get(0))));
      assertEquals(2.0 * (nArgs + 1), dmCall(g, 0.3).to_double());
      checkCodegen(g, ins(0.3), new CG().lib(libname));
      checkSerialize(externalfun1, ins(0.3));
      checkSerialize(g, ins(0.3));
      externalfun1.delete();
      g.delete();
      System.gc();
    }
  }

  @Test
  public void test_cache() {
    MX x = MX.sym("x");
    Function f = new Function("f", List.of(x), List.of(x.power(2.0)));
    Function ff = f.forward(1);
    Function ff2 = f.forward(1);
    long h = ff.hash();
    assertEquals(ff.hash(), ff2.hash());
    ff.delete();
    ff2.delete();
    System.gc();
    Function f2a = f.forward(2);
    Function f2b = f.forward(2);
    ff = f.forward(1);
    message(h + " " + ff.hash());
    assertEquals(f2a.hash(), f2b.hash());
  }

  private static String gradientCode(Function B, String fname) {
    int N = 4;
    MX z = MX.sym("z", 2);
    MX s = z;
    for (int i = 0; i < N; i++) s = B.call(List.of(s)).get(0);
    Function A = new Function("A", List.of(z), List.of(ca.dot(s, s)));
    Function g = A.factory("g", List.of("i0"), List.of("o0", "grad:o0:i0"));
    CodeGenerator cg = new CodeGenerator(fname);
    cg.add(g);
    cg.generate();
    cg.delete();
    g.delete();
    A.delete();
    return read(fname);
  }

  private static int count(String s, String sub) {
    int n = 0;
    for (int i = s.indexOf(sub); i >= 0; i = s.indexOf(sub, i + sub.length())) n++;
    return n;
  }

  @Test
  public void test_cache_stale_entry() {
    MX x = MX.sym("x", 2);
    Function B = new Function("B", List.of(x), List.of(ca.vertcat(List.of(ca.sin(x.get(0)).times(x.get(1)), ca.cos(x.get(1)).plus(x.get(0))))), List.of("x"), List.of("y"));
    assertEquals(1, count(gradientCode(B, "adj_cache1.c"), "/* adj1_B:"));
    System.gc();
    assertEquals(1, count(gradientCode(B, "adj_cache2.c"), "/* adj1_B:"));
  }

  private static final String COPY_ELISION_1 =
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

  private static final String COPY_ELISION_2 =
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

  @Test
  public void test_copy_elision() {
    long backup = GlobalOptions.getCopyElisionMinSize();
    try {
      for (boolean enabled : new boolean[] {false, true}) {
        GlobalOptions.setCopyElisionMinSize(enabled ? 2 : -1);
        MX x = MX.sym("x", 10);
        MX a = MX.sym("a", 5);
        MX bs = MX.sym("b", 10);
        MX b = bs.times(12.0).monitor("b");
        MX y = x.times(6.0);
        List<MX> s01 = ca.vertsplit(x, new long[] {0, 7, 10});
        List<MX> ss = ca.vertsplit(ca.sin(x), new long[] {0, 5, 10});
        MX z = y.get(new Slice(2, 7)).plus(a).plus(b.get(new Slice(2, 7)));
        z = ca.sum1(z);
        MX z2 = ca.vertcat(List.of(x, a, b)).times(9.0);
        Function f = new Function("f", List.of(x, a, bs), List.of(z, z2, x, ca.sqrt(s01.get(1)), ss.get(0).plus(ss.get(1))));
        f.disp(true);
        f.generate("f.c");
        String code = read("f.c");
        assertEquals(enabled, code.contains("wr3"));
        if (enabled) {
          message(code);
          assertTrue(code.contains(COPY_ELISION_1));
        }
        checkCodegen(f, randInputs(f), new CG());
        MX X = MX.sym("x", 100, 10);
        MX i = MX.sym("i");
        MX yy = X.get(false, new Slice(4, 90), i);
        f = new Function("f", List.of(i, X), List.of(ca.sin(yy).plus(X.get(false, 7, i))));
        f.disp(true);
        f.generate("f.c");
        code = read("f.c");
        assertEquals(enabled, code.contains("wr"));
        Matcher m = Pattern.compile("\\*sz_w = (\\d+);").matcher(code);
        assertTrue(m.find());
        assertEquals(enabled ? 86 : 1086, Integer.parseInt(m.group(1)));
        if (enabled) {
          message(code);
          assertTrue(code.contains(COPY_ELISION_2));
        }
        checkCodegen(f, ins(3, DM.rand(f.sparsity_in(1))), new CG().std("c99"));
      }
    } finally {
      GlobalOptions.setCopyElisionMinSize(backup);
    }
  }

  private static List<List<Double>> testPoints(List<double[]> knotsOrig) {
    List<java.util.function.Function<double[], Double>> selector = List.of(e -> e[0] - 0.1, e -> e[0], e -> (e[0] + e[1]) / 2, e -> e[e.length - 1], e -> e[e.length - 1] + 0.1);
    List<List<Double>> r = new ArrayList<>();
    r.add(new ArrayList<>());
    for (double[] k : knotsOrig) {
      List<List<Double>> next = new ArrayList<>();
      for (List<Double> prefix : r)
        for (java.util.function.Function<double[], Double> sel : selector) {
          List<Double> p = new ArrayList<>(prefix);
          p.add(sel.apply(k));
          next.add(p);
        }
      r = next;
    }
    return r;
  }

  private static double[] knotsExpand(double[] k) {
    double[] r = new double[k.length + 6];
    for (int i = 0; i < 3; i++) {
      r[i] = k[0];
      r[k.length + 3 + i] = k[k.length - 1];
    }
    System.arraycopy(k, 0, r, 3, k.length);
    return r;
  }

  private static DM vcat(List<Double> a) {
    return new DM(a.stream().mapToDouble(Double::doubleValue).toArray());
  }

  private static List<List<Double>> sample(List<List<Double>> all, int n) {
    List<List<Double>> r = new ArrayList<>(all);
    java.util.Collections.shuffle(r, new java.util.Random(1));
    return r.subList(0, n);
  }

  private static double[] linspace(double a, double b, int n) {
    double[] r = new double[n];
    for (int i = 0; i < n; i++) r[i] = a + (b - a) * i / (n - 1);
    return r;
  }

  private static long[] fill(int n, long v) {
    long[] r = new long[n];
    Arrays.fill(r, v);
    return r;
  }

  private static Map<String, Object> blazingJitOptions() {
    return opts("flags", List.of("-I" + includeDir(), "-g", "-ffast-math", "-march=native"));
  }

  @SuppressWarnings("unchecked")
  private static List<MX> parametricExprs(List<Object> r) {
    return (List<MX>) r.get(2);
  }

  @Test
  public void test_blazing_spline() {
    Assumptions.assumeTrue(CasadiMeta.feature_list().contains("simde"));
    Assumptions.assumeTrue(Importer.has_plugin("shell"));
    Assumptions.assumeTrue(RUN_SLOW);
    List<double[]> allKnots = List.of(new double[] {0, 0.2, 0.5, 0.8, 1}, new double[] {0, 0.1, 0.5, 0.9, 1}, new double[] {0, 0.4, 1}, new double[] {0, 0.3, 1},
        new double[] {0, 0.2, 1});
    List<double[]> allKnotsAlt = List.of(new double[] {0, 0.3, 0.6, 0.7, 1}, new double[] {0, 0.2, 0.4, 0.8, 1}, new double[] {0, 0.5, 1}, new double[] {0, 0.4, 1},
        new double[] {0, 0.3, 1});
    double[] x0all = {0.4, 0.5, 0.6, 0.7, 0.8};
    for (int N = 1; N <= 5; N++) {
      for (boolean precomputeCoeff : new boolean[] {false, true}) {
        for (boolean precomputeGrid : new boolean[] {false, true}) {
          if (precomputeCoeff && N > 3) continue;
          for (boolean parametricGrid : new boolean[] {false, true}) {
            if (parametricGrid && precomputeCoeff && N > 3) continue;
            DM x0 = new DM(Arrays.copyOf(x0all, N));
            List<double[]> knotsOrig = allKnots.subList(0, N);
            double[][] knots = new double[N][];
            for (int k = 0; k < N; k++) knots[k] = knotsExpand(knotsOrig.get(k));
            MX x = MX.sym("x", N);
            int nc = 1;
            for (double[] k : knots) nc *= k.length - 4;
            DM.rng(1);
            DM data = DM.rand(nc);
            MX C = MX.sym("C", nc, 1);
            MX Y = ca.bspline(x, C, knots, fill(N, 3), 1);
            Function Fref = new Function("f", List.of(x, C), List.of(Y));
            Map<String, Object> o = opts("precompute_coeff", precomputeCoeff, "precompute_grid", precomputeGrid, "pedantic_mode_order", "ignore", "pedantic_mode_size", "ignore",
                "jit", true, "jit_options", blazingJitOptions());
            Function F, Ftest, Fref2 = null;
            MX K = null;
            DM knotsData = null;
            if (parametricGrid) {
              long[] knotDims = new long[N];
              List<Double> stacked = new ArrayList<>();
              for (int k = 0; k < N; k++) {
                knotDims[k] = knots[k].length;
                for (double e : knots[k]) stacked.add(e);
              }
              F = ca.blazing_spline("F", knotDims, o);
              knotsData = vcat(stacked);
              K = MX.sym("K", stacked.size());
              List<MX> Fcall = F.call(List.of(x, C, K));
              List<MX> Fsub = ca.substitute(List.of(Fcall.get(0)), List.of(K), List.of(new MX(knotsData)));
              Ftest = new Function("F_test", List.of(x, C), Fsub);
            } else {
              Fref2 = new Function("f", List.of(x, C), List.of(Y, ca.jacobian(Y, x)));
              Fref2.call(List.of(DM.zeros(N), data));
              F = ca.blazing_spline("F", knots, o);
              Ftest = new Function("F_test", List.of(x, C), List.of(F.call(List.of(x, C)).get(0)));
            }
            List<List<Double>> allPoints = testPoints(knotsOrig);
            if (N > 3) allPoints = sample(allPoints, 125);
            for (List<Double> a : allPoints) checkfunctionLight(Ftest, Fref, List.of(vcat(a), data));
            if (!parametricGrid) {
              Function F2 = Function.deserialize(F.serialize(opts("debug", true)));
              for (List<Double> a : allPoints) checkfunctionLight(F, F2, List.of(vcat(a), data));
            }
            List<Double> last = allPoints.get(allPoints.size() - 1);
            checkCodegen(Ftest, List.of(vcat(last), data), new CG().std("c99"));
            MX ySin = parametricGrid ? ca.substitute(ca.sin(F.call(List.of(x, C, K)).get(0)), K, new MX(knotsData)) : ca.sin(F.call(List.of(x, C)).get(0));
            Function f = new Function("f", List.of(x, C), List.of(ySin, ca.jacobian(ySin, x), ca.gradient(ySin, x)));
            Function fcse = new Function("fcse", List.of(x, C), List.of(ySin, ca.jacobian(ySin, x), ca.gradient(ySin, x)), opts("cse", true));
            checkfunctionLight(f, fcse, List.of(x0, data));
            List<Object> hs = ca.hessian(ySin, x);
            List<MX> outs = new ArrayList<>(List.of(ySin, ca.jacobian(ySin, x), ca.gradient(ySin, x), (MX) hs.get(0), (MX) hs.get(1)));
            Function f2 = new Function("f", List.of(x, C), outs);
            Function fcse2 = new Function("fcse", List.of(x, C), outs, opts("cse", true));
            checkfunctionLight(f2, fcse2, List.of(x0, data));
            String result = captureStdout(() -> fcse2.disp(true));
            assertEquals(1, count(result, " F_der_der("));
            assertEquals(0, count(result, " F_der("));
            assertEquals(0, count(result, " F("));
            result = captureStdout(() -> fcse.disp(true));
            assertEquals(0, count(result, " F_der_der("));
            assertEquals(1, count(result, " F_der("));
            assertEquals(0, count(result, " F("));
            if (parametricGrid) {
              MX yJ = ca.jacobian(F.call(List.of(x, C, K)).get(0), x);
              List<MX> pe = parametricExprs(ca.extract_parametric(yJ, ca.vertcat(List.of(C, K))));
              if (precomputeCoeff || precomputeGrid) {
                assertTrue(pe.size() > 0);
                assertTrue(pe.stream().anyMatch(e -> !e.is_zero()));
              } else {
                assertEquals(0, pe.size());
              }
            } else {
              MX yJ = ca.jacobian(F.call(List.of(x, C)).get(0), x);
              List<MX> pe = parametricExprs(ca.extract_parametric(yJ, C));
              if (precomputeCoeff) {
                assertTrue(pe.size() > 0);
                assertTrue(pe.stream().anyMatch(e -> !e.is_zero()));
              } else {
                assertEquals(0, pe.size());
              }
            }
            Function FJref = Fref.jacobian();
            Function FJtest, FJ = null;
            if (parametricGrid) {
              Function FJparam = F.jacobian();
              MX adjF = MX.sym("adj_f", 1, 1);
              List<MX> fjOut = FJparam.call(List.of(x, C, K, adjF));
              MX fjX = ca.substitute(fjOut.get(0), K, new MX(knotsData));
              MX fjC = ca.substitute(fjOut.get(1), K, new MX(knotsData));
              FJtest = new Function("FJ_test", List.of(x, C, adjF), List.of(fjX, fjC));
            } else {
              FJ = F.jacobian();
              FJ.generate("FJ.c", opts("main", true));
              FJ.generate_in("FJ_in.txt", List.of(vcat(last), data, new DM(0.0)));
              FJtest = FJ;
            }
            for (List<Double> a : allPoints) checkfunctionLight(FJtest, FJref, List.of(vcat(a), data, new DM(0.0)));
            if (!parametricGrid) {
              message("ref");
              Fref2.call(List.of(DM.zeros(N), data));
              FJ = FJ.jacobian();
              FJ.generate("FH.c", opts("main", true));
              FJ.generate_in("FH_in.txt", ins(x0, data, 0, 0, 0));
              message(FJ.toString());
              FJref = FJref.jacobian();
              long t0 = System.nanoTime();
              FJref.call(ins(vcat(last), data, 0, 0, 0));
              message("FJ_ref " + (System.nanoTime() - t0) * 1e-9);
              t0 = System.nanoTime();
              FJ.call(ins(vcat(last), data, 0, 0, 0));
              message("FJ " + (System.nanoTime() - t0) * 1e-9);
              for (List<Double> a : allPoints) checkfunctionLight(FJ, FJref, ins(vcat(a), data, 0, 0, 0));
            }
            if (parametricGrid) {
              List<Double> altStacked = new ArrayList<>();
              for (int k = 0; k < N; k++) for (double e : knotsExpand(allKnotsAlt.get(k))) altStacked.add(e);
              MX yAlt = ca.substitute(F.call(List.of(x, C, K)).get(0), K, new MX(vcat(altStacked)));
              Function Falt = new Function("F_alt", List.of(x, C), List.of(yAlt));
              double res1 = dmCall(Ftest, x0, data).to_double();
              double res2 = dmCall(Falt, x0, data).to_double();
              assertFalse(Math.abs(res1 - res2) <= 1e-8 + 1e-5 * Math.abs(res2), "Parametric knots should produce different results with different knot vectors");
            }
          }
          MX xLu = MX.sym("x", N);
          List<double[]> knotsUniformOrig = List.of(linspace(0, 1, 5), linspace(0, 1, 5), linspace(0, 1, 3), linspace(0, 1, 3), linspace(0, 1, 3)).subList(0, N);
          double[][] knotsUniform = new double[N][];
          int ncU = 1;
          for (int k = 0; k < N; k++) {
            knotsUniform[k] = knotsExpand(knotsUniformOrig.get(k));
            ncU *= knotsUniform[k].length - 4;
          }
          DM.rng(1);
          DM dataU = DM.rand(ncU);
          MX Cu = MX.sym("C", ncU, 1);
          MX Yu = ca.bspline(xLu, Cu, knotsUniform, fill(N, 3), 1);
          Function FrefU = new Function("f", List.of(xLu, Cu), List.of(Yu));
          List<List<Double>> pointsU = testPoints(knotsUniformOrig);
          if (N > 3) pointsU = sample(pointsU, 125);
          for (String lookupMode : List.of("linear", "exact", "binary")) {
            List<String> modes = new ArrayList<>();
            for (int k = 0; k < N; k++) modes.add(lookupMode);
            Map<String, Object> optsLu = opts("precompute_coeff", precomputeCoeff, "precompute_grid", precomputeGrid, "lookup_mode", modes, "pedantic_mode_order", "ignore",
                "pedantic_mode_size", "ignore", "jit", true, "jit_options", blazingJitOptions());
            Function Flu = ca.blazing_spline("F", knotsUniform, optsLu);
            Function FluTest = new Function("F_lu_test", List.of(xLu, Cu), List.of(Flu.call(List.of(xLu, Cu)).get(0)));
            for (List<Double> a : pointsU) checkfunctionLight(FluTest, FrefU, List.of(vcat(a), dataU));
            Function Flu2 = Function.deserialize(Flu.serialize(opts("debug", true)));
            for (List<Double> a : pointsU) checkfunctionLight(Flu, Flu2, List.of(vcat(a), dataU));
          }
        }
      }
    }
  }

  /** stdout and stderr of r */
  private static String[] captureOutErr(Runnable r) {
    java.io.PrintStream oldErr = System.err;
    java.io.ByteArrayOutputStream b = new java.io.ByteArrayOutputStream();
    System.setErr(new java.io.PrintStream(b, true, StandardCharsets.UTF_8));
    String out;
    try {
      out = captureStdout(r);
    } finally {
      System.err.flush();
      System.setErr(oldErr);
    }
    return new String[] {out, b.toString(StandardCharsets.UTF_8)};
  }

  private static boolean inEither(String s, String[] r) {
    return r[0].contains(s) || r[1].contains(s);
  }

  private static double[][] range(int... n) {
    double[][] r = new double[n.length][];
    for (int k = 0; k < n.length; k++) {
      r[k] = new double[n[k]];
      for (int i = 0; i < n[k]; i++) r[k][i] = i;
    }
    return r;
  }

  @Test
  public void test_blazing_spline_pedantic() {
    Assumptions.assumeTrue(CasadiMeta.feature_list().contains("simde"));
    double[][] sizeBad1d = range(12);
    double[][] orderBad2d = range(17, 13);
    double[][] cumulBad2d = range(8, 8);
    double[][] clean2d = range(13, 17);
    double[][] clean1d = range(13);
    assertInException("powers of 2", () -> ca.blazing_spline("F", sizeBad1d, opts("pedantic_mode_size", "error")));
    assertInException("not increasing", () -> ca.blazing_spline("F", orderBad2d, opts("pedantic_mode_order", "error", "pedantic_mode_size", "ignore")));
    assertInException("dim 0 (zero-based)", () -> ca.blazing_spline("F", sizeBad1d, opts("pedantic_mode_size", "error")));
    assertInException("'pedantic_mode_size' must be one of", () -> ca.blazing_spline("F", sizeBad1d, opts("pedantic_mode_size", "bogus")));
    assertInException("prefix product", () -> ca.blazing_spline("F", cumulBad2d, opts("pedantic_mode_size", "error")));
    assertTrue(inEither("powers of 2", captureOutErr(() -> ca.blazing_spline("F", sizeBad1d, opts("pedantic_mode_size", "warn")))));
    assertTrue(inEither("not increasing", captureOutErr(() -> ca.blazing_spline("F", orderBad2d, opts("pedantic_mode_size", "ignore")))));
    assertFalse(inEither("powers of 2", captureOutErr(() -> ca.blazing_spline("F", sizeBad1d, opts("pedantic_mode_size", "ignore", "pedantic_mode_order", "ignore")))));
    ca.blazing_spline("F", clean1d, opts("pedantic_mode_order", "error", "pedantic_mode_size", "error"));
    ca.blazing_spline("F", clean2d, opts("pedantic_mode_order", "error", "pedantic_mode_size", "error"));
    ca.blazing_spline("F", clean1d, opts("pedantic_mode_order", "error"));
    double[][] n5bad1d = range(13);
    assertInException("(n_knots - 5)", () -> ca.blazing_spline("F", n5bad1d).jacobian());
    assertInException("diff order 1", () -> ca.blazing_spline("F", n5bad1d).jacobian());
    assertTrue(inEither("(n_knots - 5)", captureOutErr(() -> ca.blazing_spline("F", n5bad1d, opts("pedantic_mode_size", "warn")).jacobian())));
    assertFalse(inEither("(n_knots - 5)", captureOutErr(() -> ca.blazing_spline("F", n5bad1d, opts("pedantic_mode_size", "ignore")).jacobian())));
    assertFalse(inEither("(n_knots - 5)", captureOutErr(() -> ca.blazing_spline("F", n5bad1d, opts("jacobian_options", opts("pedantic_mode_size", "ignore"))).jacobian())));
    double[][] jacClean1d = range(15);
    ca.blazing_spline("F", jacClean1d, opts("pedantic_mode_size", "error")).jacobian();
    double[][] n6bad1d = range(14);
    assertInException("(n_knots - 6)", () -> ca.blazing_spline("F", n6bad1d, opts("precompute_coeff", false)).jacobian().jacobian());
    assertInException("diff order 2", () -> ca.blazing_spline("F", n6bad1d, opts("precompute_coeff", false)).jacobian().jacobian());
  }

  @Test
  public void test_noncanonical_sparsity() {
    MX x = MX.sym("x", 4, 4);
    MX y = MX.sym("y");
    Function f = new Function("f", List.of(x, y), List.of(new MX(24.0)));
    checkCodegen(f, ins(DM.rand(4, 4), 1), new CG().opt("force_canonical", false));
    checkCodegen(f, ins(DM.rand(4, 4), 1), new CG().opt("force_canonical", true));
  }

  @Test
  public void test_options_sanitize() {
    assertEquals(opts("foo", opts("bar", opts("baz", opts("goot", 7L)))), Options.sanitize(opts("foo", opts("bar", opts("baz.goot", 7)))));
    assertEquals(opts("foo", opts("bar", 9L, "baz", 1L)), Options.sanitize(opts("foo.baz", 1, "foo", opts("bar", 9))));
    assertEquals(opts("foo", opts("baz", opts("bar", 1L, "go", 8L))), Options.sanitize(opts("foo.baz.bar", 1, "foo", opts("baz", opts("go", 8)))));
    assertEquals(opts("foo", opts("bar", 9L, "baz", 1L)), Options.sanitize(opts("foo", opts("bar", 9), "foo.baz", 1)));
    assertEquals(opts("bar", opts("abc", opts("def", 7L)), "foo", opts("bar", opts("baz", 1L, "mat", 3L, "got", 9L, "qint", opts("pad", 5L)), "lib", 5L)),
        Options.sanitize(opts("foo.bar.baz", 1, "foo.bar.mat", 3, "foo.lib", 5, "foo.bar.qint.pad", 5, "foo", opts("bar", opts("got", 9)), "bar", opts("abc.def", 7))));
    assertInException("update_dict error", () -> Options.sanitize(opts("foo", 9, "foo.baz", 1)));
    assertInException("update_dict error", () -> Options.sanitize(opts("foo", opts("baz", opts("w", 5)), "foo.baz", 1)));
    Map<String, Object> res = Options.sanitize(opts("foo", 1, "bar", null));
    message(res.toString());
    assertEquals(opts("foo", 1L), res);
    res = Options.sanitize(opts("foo", 1, "bar", opts("r", null)));
    message(res.toString());
    assertEquals(opts("foo", 1L, "bar", opts("r", null)), res);
  }

  private static List<List<Object>> passes(Object[]... p) {
    List<List<Object>> r = new ArrayList<>();
    for (Object[] e : p) r.add(new ArrayList<>(Arrays.asList(e)));
    return r;
  }

  @Test
  public void test_simplify() {
    MX x = MX.sym("x", 5);
    MX y = MX.sym("y", 5);
    Function f = new Function("test_func", List.of(x, y), List.of(x.plus(y)));
    Function fSimplified = f.transform();
    assertEquals("test_func", fSimplified.name());
    Function fSameName = f.transform("test_func");
    assertEquals("test_func", fSameName.name());
    Function fRenamed = f.transform("new_name");
    assertEquals("new_name", fRenamed.name());
    checkfunctionLight(f, fSimplified, ins(2.0, 3.0));
    x = MX.sym("x", 5);
    y = MX.sym("y", 5);
    f = new Function("test_func", List.of(x, y), List.of(x));
    fSimplified = f.transform(passes(new Object[] {"simplify", "empty_inputs"}));
    assertEquals(5, f.nnz_in(0));
    assertEquals(5, fSimplified.nnz_in(0));
    assertEquals(0, fSimplified.nnz_in(1));
    checkfunctionLight(f, fSimplified, ins(2.0, 3.0));
  }

  @Test
  public void test_simplify_ref_count() {
    List<Function[]> tests = new ArrayList<>();
    for (String X : XS) {
      if (X.equals("SX")) {
        SX a = SX.sym("a"), b = SX.sym("b");
        tests.add(new Function[] {new Function("f", List.of(a, b), List.of(a.minus(b).unaryMinus())), new Function("f", List.of(a, b), List.of(b.minus(a)))});
        SX e = a.minus(b);
        SX z = e.unaryMinus();
        tests.add(new Function[] {new Function("f", List.of(a, b), List.of(e.times(3.0), e.unaryMinus())), new Function("f", List.of(a, b), List.of(e.times(3.0), z))});
      } else {
        MX a = MX.sym("a"), b = MX.sym("b");
        tests.add(new Function[] {new Function("f", List.of(a, b), List.of(a.minus(b).unaryMinus())), new Function("f", List.of(a, b), List.of(b.minus(a)))});
        MX e = a.minus(b);
        MX z = e.unaryMinus();
        tests.add(new Function[] {new Function("f", List.of(a, b), List.of(e.times(3.0), e.unaryMinus())), new Function("f", List.of(a, b), List.of(e.times(3.0), z))});
      }
    }
    MX x = MX.sym("x", Sparsity.lower(3));
    MX y = ca.project(x, Sparsity.upper(3));
    tests.add(new Function[] {new Function("f", List.of(x), List.of(y.get_nz(false, 5))), new Function("f", List.of(x), List.of(x.get_nz(false, 5)))});
    x = MX.sym("x", Sparsity.lower(3));
    y = ca.project(x, Sparsity.upper(3));
    tests.add(new Function[] {new Function("f", List.of(x), List.of(y.get(new Slice(), -1))), new Function("f", List.of(x), List.of(x.nzref(Sparsity.dense(3, 1), new long[] {-1, -1, 5})))});
    for (Function[] t : tests) {
      Function f = t[0], ref = t[1];
      ref.disp(true);
      Function fs = f.transform();
      fs.disp(true);
      fs.generate("f.c");
      ref.generate("ref.c");
      checkfunctionIdentical(fs, ref);
      DM.rng(1);
      checkfunctionLight(f, fs, randInputs(f));
    }
  }

  @Test
  public void test_simplify_const_folding() {
    DM.rng(1);
    MX A = new MX(DM.rand(2, 2));
    MX B = new MX(DM.rand(2, 2));
    MX C = ca.mtimes(A, B);
    C = ca.mtimes(C, C);
    MX x = MX.sym("x");
    Function f = new Function("f", List.of(x), List.of(x.times(C.times(2.0)).plus(C.times(2.0)).times(C.get(0))));
    Function fs = f.transform();
    checkfunctionLight(f, fs, randInputs(f));
    fs.disp(true);
    assertEquals(7, fs.n_nodes());
  }

  @Test
  public void test_simplify_combine_terms() {
    SX x = SX.sym("x");
    SX y = SX.sym("y");
    java.util.function.BiConsumer<SX, SX> numeq = (expr, ref) -> {
      Function fe = new Function("fe", List.of(x, y), List.of(expr));
      Function fr = new Function("fr", List.of(x, y), List.of(ref));
      for (double[] v : new double[][] {{0.3, 1.7}, {-2.0, 4.5}}) checkarray(dmCall(fe, v[0], v[1]), dmCall(fr, v[0], v[1]));
    };
    java.util.function.ToLongFunction<SX> nnodes = expr -> new Function("f", List.of(x, y), List.of(expr)).n_nodes();
    SX e1 = y.times(8.0).minus(y.times(3.0));
    SX s1 = ca.simplify(e1);
    numeq.accept(s1, y.times(5.0));
    assertTrue(nnodes.applyAsLong(s1) < nnodes.applyAsLong(e1));
    SX e2 = x.times(5.0).times(8.0).minus(x.plus(3.0).plus(x.times(12.0).minus(x)).minus(x.times(17.0)).times(3.0));
    SX s2 = ca.simplify(e2);
    numeq.accept(s2, x.times(55.0).minus(9.0));
    assertTrue(nnodes.applyAsLong(s2) < nnodes.applyAsLong(e2));
    Function f = new Function("f", List.of(x), List.of(x.times(8.0).minus(x.times(3.0))));
    Function fc = f.transform(passes(new Object[] {"simplify", "combine_terms"}));
    checkfunctionLight(f, fc, ins(2.0));
    assertTrue(fc.n_nodes() < f.n_nodes());
    MX a = MX.sym("a");
    MX b = MX.sym("b");
    MX e = a.minus(b).unaryMinus();
    assertTrue(ca.is_equal(ca.simplify(e), e, 2));
  }

  @Test
  public void test_transform() {
    SX x = SX.sym("x");
    SX y = SX.sym("y");
    SX s = ca.transform(y.times(8.0).minus(y.times(3.0)), passes(new Object[] {"simplify", "combine_terms"}));
    Function fe = new Function("fe", List.of(x, y), List.of(s));
    Function fr = new Function("fr", List.of(x, y), List.of(y.times(5.0)));
    double[][] vals = {{0.3, 1.7}, {-2.0, 4.5}};
    for (double[] v : vals) checkarray(dmCall(fe, v[0], v[1]), dmCall(fr, v[0], v[1]));
    List<SX> sv = ca.transform(List.of(y.times(8.0).minus(y.times(3.0)), x.plus(x)), passes(new Object[] {"simplify", "combine_terms"}));
    assertEquals(2, sv.size());
    Function fev = new Function("fev", List.of(x, y), sv);
    Function frv = new Function("frv", List.of(x, y), List.of(y.times(5.0), x.times(2.0)));
    for (double[] v : vals) {
      checkarray(dmCallN(fev, v[0], v[1]).get(0), dmCallN(frv, v[0], v[1]).get(0));
      checkarray(dmCallN(fev, v[0], v[1]).get(1), dmCallN(frv, v[0], v[1]).get(1));
    }
    assertEquals(2, ca.transform(List.of(y.times(8.0).minus(y.times(3.0)), x.plus(x)), opts("combine_terms", true)).size());
    MX a = MX.sym("a");
    MX b = MX.sym("b");
    Function f = new Function("f", List.of(a, b), List.of(a.minus(b).unaryMinus()));
    Function ft = f.transform(passes(new Object[] {"simplify", 0, "ref_count", "const_folding"}));
    checkfunctionIdentical(ft, new Function("f", List.of(a, b), List.of(b.minus(a))));
    assertTrue(ft.n_nodes() < f.n_nodes());
    Function ft2 = f.transform(passes(new Object[] {"simplify", 2, "ref_count"}));
    checkfunctionIdentical(ft2, new Function("f", List.of(a, b), List.of(b.minus(a))));
    assertEquals("f", f.transform().name());
    assertEquals("f", f.transform(passes(new Object[] {"simplify", "cse"})).name());
    assertEquals("g", f.transform("g").name());
    assertEquals("g", f.transform("g", passes(new Object[] {"simplify", "cse"})).name());
    assertEquals("g", f.transform("g", opts("cse", true)).name());
    Function frr = f.transform("g", passes(new Object[] {"simplify", 0, "ref_count", "const_folding"}));
    checkfunctionIdentical(frr, new Function("g", List.of(a, b), List.of(b.minus(a))));
    Function fdNoarg = f.transform();
    Function fdEmpty = f.transform(opts());
    Function fdBool = f.transform(opts("combine_terms", false));
    for (Function fd : List.of(fdNoarg, fdEmpty, fdBool)) checkfunctionIdentical(fd, new Function("f", List.of(a, b), List.of(b.minus(a))));
    Function fdNone = f.transform(passes());
    checkfunctionIdentical(fdNone, f);
    SX sx = SX.sym("s");
    Function g = new Function("g", List.of(sx), List.of(sx.times(8.0).minus(sx.times(3.0))));
    Function gc = g.transform(opts("combine_terms", true));
    checkfunctionLight(g, gc, ins(2.0));
    assertTrue(gc.n_nodes() < g.n_nodes());
    assertThrows(RuntimeException.class, () -> g.transform(opts("bogus", true)));
    assertThrows(RuntimeException.class, () -> g.transform(passes(new Object[] {"simplify", "cse"}), opts("cse", true)));
    Function fx = f.transform(passes(new Object[] {"expand"}));
    assertEquals("SXFunction", fx.class_name());
    checkfunctionLight(f, fx, ins(1.0, 2.0));
    assertThrows(RuntimeException.class, () -> f.transform(passes(new Object[] {"bogus"})));
    assertThrows(RuntimeException.class, () -> f.transform(passes(new Object[] {"simplify", "nope"})));
    assertThrows(RuntimeException.class, () -> f.transform(passes(new Object[] {"simplify", -1, "ref_count"})));
    if (!System.getProperty("os.name").toLowerCase().contains("mac")) {
      String libcasadi = CasadiMeta.shared_library_prefix() + "casadi";
      Function g2 = new Function("g", List.of(x), List.of(x.power(2.0)));
      Function[] gt = new Function[1];
      String[] r = captureOutErr(() -> gt[0] = g2.transform(passes(new Object[] {"external", libcasadi, "external_transform_test_success", opts("foo", "bar")})));
      for (String e : List.of("Doing", "bar")) assertTrue(r[0].contains(e), r[0] + "<->" + e);
      assertTrue(r[1].contains("Warning"), r[1] + "<->Warning");
      checkfunction(g2, gt[0], ins(3));
    }
    String[] r = captureOutErr(() -> f.transform(passes(new Object[] {"simplify", 0, "ref_count"}), opts("verbose", true)));
    for (String e : List.of("transform(", "n_nodes")) assertTrue(r[0].contains(e), r[0] + "<->" + e);
    r = captureOutErr(() -> f.transform(opts("verbose", true)));
    for (String e : List.of("transform(", "n_nodes")) assertTrue(r[0].contains(e), r[0] + "<->" + e);
  }

  @Test
  public void test_duplicate_check() {
    MX x = new MX();
    Function f = new Function("f", List.of(x), List.of(x.times(2.0)));
    assertEquals(0, f.numel_in(0));
    f = new Function("f", List.of(MX.zeros(0, 0)), List.of(x.times(2.0)));
    assertEquals(0, f.numel_in(0));
    f = new Function("f", List.of(MX.ones(0, 0)), List.of(x.times(2.0)));
    assertEquals(0, f.numel_in(0));
    MX p = MX.sym("p");
    assertInException("input arguments must be purely symbolic", () -> {
      Function g = new Function("f", List.of(x.get_nz(false, p)), List.of(x.times(2.0)));
      assertEquals(0, g.numel_in(0));
    });
    x.set_nz(3.7, false, p);
    assertInException("input arguments must be purely symbolic", () -> {
      Function g = new Function("f", List.of(x), List.of(x.times(3.0)));
      assertEquals(0, g.numel_in(0));
    });
  }

  @Test
  public void test_finite_diff() {
    MX x = MX.sym("x");
    for (String method : List.of("forward", "backward", "central", "smoothing")) {
      Function f = new Function("f", List.of(x), List.of(ca.sin(x)), opts("enable_fd", true, "enable_forward", false, "enable_reverse", false, "fd_method", method));
      Function fref = new Function("f", List.of(x), List.of(ca.sin(x)));
      Function J = f.jacobian();
      Function Jref = fref.jacobian();
      List<DM> inputs = ins(1.1, 0);
      checkfunctionLight(J, Jref, inputs, 5);
      if (RUN_SLOW) checkCodegen(J, inputs, new CG().std("c99").digits(15));
    }
  }

  @Test
  public void test_finite_diff_serialization() {
    MX x = MX.sym("x", 2);
    MX p = MX.sym("p");
    MX y = p.times(ca.dot(x, x));
    for (String method : List.of("forward", "backward", "central", "smoothing")) {
      for (boolean adaptive : new boolean[] {false, true}) {
        double h = method.equals("backward") ? -1e-3 : 1e-3;
        Map<String, Object> options = opts("h", h, "h_iter", (adaptive && (method.equals("central") || method.equals("smoothing"))) ? 1 : 0, "h_min", 1e-8, "h_max", .1,
            "u_aim", 42., "reltol", 1e-11, "abstol", 1e-12, "smoothing", 5e-8);
        Function f = new Function("f", List.of(x, p), List.of(y, x.plus(p)), opts("enable_fd", true, "enable_forward", false, "enable_reverse", false, "is_diff_in",
            new boolean[] {true, false}, "is_diff_out", new boolean[] {true, false}, "fd_method", method, "fd_options", options));
        List<DM> values = List.of(dm(.4, -.7), new DM(2));
        for (int cnt : new int[] {1, 3}) {
          Function df = f.forward(cnt);
          Function restored = Function.deserialize(df.serialize());
          assertEquals(df.class_name(), restored.class_name());
          assertTrue(Arrays.equals(df.is_diff_in(), restored.is_diff_in()));
          assertTrue(Arrays.equals(df.is_diff_out(), restored.is_diff_out()));
          List<DM> inputs = new ArrayList<>(values);
          inputs.addAll(f.call(values));
          inputs.add(DM.ones(2, cnt));
          inputs.add(DM.zeros(1, cnt));
          List<DM> got = restored.call(inputs), expected = df.call(inputs);
          for (int i = 0; i < got.size(); i++) checkarray(got.get(i), expected.get(i), 1e-12, "");
          for (long i = 0; i < df.n_in(); i++) assertTrue(restored.sparsity_in(i).is_equal(df.sparsity_in(i)));
        }
        MX result = f.call(List.of(x, p)).get(0);
        Function D = new Function("D", List.of(x, p), List.of(ca.jacobian(result, x), (MX) ca.hessian(result, x).get(0), ca.jacobian(result, p)));
        Function restored = Function.deserialize(D.serialize());
        List<DM> got = restored.call(values), expected = D.call(values);
        for (int i = 0; i < got.size(); i++) checkarray(got.get(i), expected.get(i), 1e-12, "");
        if (!adaptive) checkCodegen(restored, values, new CG().std("c99").digits(10));
      }
    }
  }

  @Test
  public void test_is_diff_fd() {
    DM.set_precision(16);
    try {
      for (int n : new int[] {1, 2}) {
        for (int m : new int[] {1, 3}) {
          MX x = MX.sym("x", n, m);
          MX y = MX.sym("y", n, m);
          List<MX> expr = List.of(ca.sumsqr(x).times(2.0).times(ca.sumsqr(y)), y.power(2.0).times(ca.sumsqr(x)));
          Function fref = new Function("f", List.of(x, y), expr, List.of("x", "y"), List.of("a", "b"), opts("is_diff_in", new boolean[] {true, false}));
          Function f = new Function("f", List.of(x, y), expr, List.of("x", "y"), List.of("a", "b"), opts("is_diff_in", new boolean[] {true, false}, "enable_fd", true,
              "print_in", true, "print_out", false, "enable_forward", false, "enable_reverse", false, "fd_method", "forward"));
          Function Jref = fref.jacobian();
          Function J = f.jacobian();
          assertTrue(Jref.is_diff_in(0));
          assertFalse(Jref.is_diff_in(1));
          assertEquals(0, Jref.nnz_in(2));
          assertEquals(0, J.nnz_in(2));
          Function Ffref = fref.forward(2);
          Function Ff = f.forward(2);
          assertEquals(0, Ffref.nnz_in(2));
          assertEquals(1, Ff.nnz_in(2));
          for (Function[] pair : new Function[][] {{J, Jref}, {Ff, Ffref}}) {
            Function F = pair[0], Fref = pair[1];
            message("F = " + F);
            message("Fref = " + Fref);
            assertEquals(Fref.n_in(), F.n_in());
            assertEquals(Fref.n_out(), F.n_out());
            for (long i = 0; i < F.n_in(); i++) {
              assertEquals(Fref.name_in(i), F.name_in(i));
              if (!F.name_in(i).startsWith("out_")) assertTrue(F.sparsity_in(i).is_equal(Fref.sparsity_in(i)));
              assertEquals(Fref.is_diff_in(i), F.is_diff_in(i));
            }
            DM.rng(1);
            List<DM> inputs = new ArrayList<>();
            for (long i = 0; i < f.n_in(); i++) inputs.add(DM.rand(F.sparsity_in(i)));
            inputs.addAll(f.call(inputs));
            if (F.name().contains("fwd")) for (long i = 0; i < f.n_in(); i++) inputs.add(DM.rand(F.sparsity_in(i)));
            if (F.name().contains("adj")) for (long i = 0; i < f.n_out(); i++) inputs.add(DM.rand(F.sparsity_out(i)));
            checkfunctionLight(F, Fref, inputs, 5);
            if (RUN_SLOW) checkCodegen(F, inputs, new CG().std("c99").digits(System.getProperty("os.name").toLowerCase().contains("mac") ? 5 : 15));
            List<DM> in2 = ins(2, 3, 0, 0);
            if (F.name().contains("fwd")) for (long i = 0; i < f.n_in(); i++) in2.add(DM.rand(F.sparsity_in(i)));
            if (F.name().contains("adj")) for (long i = 0; i < f.n_out(); i++) in2.add(DM.rand(F.sparsity_out(i)));
            String result = captureStdout(() -> F.call(in2));
            assertTrue(result.contains("2.00"), result);
            assertFalse(result.contains("3.00"), result);
          }
        }
      }
    } finally {
      DM.set_precision(6);
    }
  }

  @Test
  public void test_forward_reverse_zero() {
    SX x = SX.sym("x");
    SX y = SX.sym("y", 2);
    Function f = new Function("f", List.of(x, y), List.of(x.times(y), ca.sin(x)));
    for (Function F : List.of(f, f.wrap())) {
      Function fwd = F.forward(0);
      assertTrue(Arrays.equals(new long[] {2, 0}, fwd.size_in(5)));
      assertTrue(Arrays.equals(new long[] {2, 0}, fwd.size_out(0)));
      Function adj = F.reverse(0);
      assertTrue(Arrays.equals(new long[] {2, 0}, adj.size_in(4)));
      assertTrue(Arrays.equals(new long[] {2, 0}, adj.size_out(1)));
    }
  }

  @Test
  public void test_jit_function_fd() {
    Assumptions.assumeTrue(Importer.has_plugin("shell"));
    Assumptions.assumeTrue(RUN_SLOW);
    Function f = Function.jit("f", "b[0] = a[0]*a[0];", List.of("a"), List.of("b"), opts("compiler", "shell"));
    checkarray(dmCall(f.jacobian(), 3, 0), new DM(6), 1e-5, "");
  }

  @Test
  public void test_has_function() {
    MX x = MX.sym("x");
    Function f = new Function("f", List.of(x), List.of(x.power(2.0)), opts("never_inline", true));
    Function g = new Function("g", List.of(x), List.of(f.call(List.of(x)).get(0).plus(1.0)));
    assertEquals(List.of("f"), g.get_function());
    assertTrue(g.has_function("f"));
    assertFalse(g.has_function("h"));
  }

  @Test
  public void test_string_ctor_hidden() {
    assertThrows(NoSuchMethodException.class, () -> Function.class.getConstructor(String.class));
  }

  private static boolean[] empActivity(Function f, boolean[] mask) {
    boolean[] emp = new boolean[(int) f.nnz_out()];
    for (int t = 1; t < 6; t++) {
      List<DM> in = new ArrayList<>();
      int off = 0;
      for (long i = 0; i < f.n_in(); i++) {
        int ni = (int) f.nnz_in(i);
        double[] vals = new double[ni];
        for (int j = 0; j < ni; j++) vals[j] = mask[off + j] ? (double) t * (off + j + 1) : 0.0;
        in.add(new DM(f.sparsity_in(i), new DM(vals)));
        off += ni;
      }
      List<Double> flat = new ArrayList<>();
      for (DM o : f.call(in)) for (double v : o.nonzeros()) flat.add(v);
      for (int k = 0; k < emp.length; k++) emp[k] = emp[k] || flat.get(k) != 0;
    }
    return emp;
  }

  private static List<boolean[]> allMasks(int n) {
    List<boolean[]> r = new ArrayList<>();
    for (int b = 0; b < (1 << n); b++) {
      boolean[] m = new boolean[n];
      for (int i = 0; i < n; i++) m[i] = (b & (1 << i)) != 0;
      r.add(m);
    }
    return r;
  }

  private static List<boolean[]> prependTrue(int k, List<boolean[]> masks) {
    List<boolean[]> r = new ArrayList<>();
    for (boolean[] m : masks) {
      boolean[] e = new boolean[k + m.length];
      Arrays.fill(e, 0, k, true);
      System.arraycopy(m, 0, e, k, m.length);
      r.add(e);
    }
    return r;
  }

  private static void checkActivity(Function f, List<boolean[]> masks, boolean minimal) {
    for (boolean[] m : masks) {
      boolean[] emp = empActivity(f, m), got = f.activity(m);
      for (int k = 0; k < f.nnz_out(); k++)
        assertTrue(got[k] || !emp[k], "UNSOUND: " + f.name() + " mask " + Arrays.toString(m) + " output " + k + ": activity says inactive but empirically nonzero");
      if (minimal) assertTrue(Arrays.equals(emp, got), "not minimal: " + f.name() + " mask " + Arrays.toString(m));
    }
  }

  @Test
  public void test_activity() {
    for (String X : List.of("MX", "SX")) {
      if (X.equals("MX")) {
        MX x = MX.sym("x"), y = MX.sym("y");
        checkActivity(new Function("mul", List.of(x, y), List.of(x.times(y))), allMasks(2), true);
        checkActivity(new Function("add", List.of(x, y), List.of(x.plus(y))), allMasks(2), true);
        checkActivity(new Function("cs", List.of(x), List.of(ca.cos(x))), allMasks(1), true);
        checkActivity(new Function("sn", List.of(x), List.of(ca.sin(x))), allMasks(1), true);
        checkActivity(new Function("z", List.of(x), List.of(x.times(0.0))), allMasks(1), true);
      } else {
        SX x = SX.sym("x"), y = SX.sym("y");
        checkActivity(new Function("mul", List.of(x, y), List.of(x.times(y))), allMasks(2), true);
        checkActivity(new Function("add", List.of(x, y), List.of(x.plus(y))), allMasks(2), true);
        checkActivity(new Function("cs", List.of(x), List.of(ca.cos(x))), allMasks(1), true);
        checkActivity(new Function("sn", List.of(x), List.of(ca.sin(x))), allMasks(1), true);
        checkActivity(new Function("z", List.of(x), List.of(x.times(0.0))), allMasks(1), true);
      }
    }
    MX x2 = MX.sym("x2", 2);
    checkActivity(new Function("pc", List.of(x2), List.of(x2.times(ca.densify(new MX(dm(0, 5)))))), allMasks(2), true);
    MX xr = MX.sym("xr", 1, 2), yc = MX.sym("yc", 2, 1);
    checkActivity(new Function("mt", List.of(xr, yc), List.of(ca.mtimes(xr, yc))), allMasks(4), true);
    MX A = MX.sym("A", 2, 2), b = MX.sym("b", 2, 1);
    checkActivity(new Function("slv", List.of(A, b), List.of(ca.solve(A, b))), prependTrue(4, allMasks(2)), true);
    MX ad = MX.sym("ad", 2), bd = MX.sym("bd", 2, 1);
    checkActivity(new Function("sd", List.of(ad, bd), List.of(ca.solve(ca.diag(ad), bd))), prependTrue(2, allMasks(2)), false);
    MX y = MX.sym("y", 2);
    Function foo = new Function("foo", List.of(y), List.of(ca.vertcat(List.of(ca.sin(y.get(0).times(y.get(1))), ca.cos(y.get(0).times(y.get(1)))))));
    Function fadj = foo.reverse(1).forward(1);
    checkActivity(fadj, allMasks((int) fadj.nnz_in()), true);
    List<MX> args = new ArrayList<>();
    for (long i = 0; i < fadj.n_in(); i++) {
      if (fadj.name_in(i).startsWith("fwd_")) args.add(new MX(fadj.size1_in(i), fadj.size2_in(i)));
      else args.add(MX.sym("a" + i, fadj.sparsity_in(i)));
    }
    assertTrue(fadj.call(args).get(0).is_zero());
    List<MX> argsNz = new ArrayList<>();
    for (long i = 0; i < fadj.n_in(); i++) argsNz.add(MX.sym("b" + i, fadj.sparsity_in(i)));
    assertFalse(fadj.call(argsNz).get(0).is_zero());
    MX x = MX.sym("x", 2);
    MX lam = MX.sym("lam", 2), sig = MX.sym("sig");
    MX L = sig.times(ca.sumsqr(x.minus(2.0))).plus(ca.dot(lam, ca.no_hess(foo.call(List.of(x.power(2.0))).get(0))));
    MX H = (MX) ca.hessian(L, x).get(0);
    Function Hf = new Function("H", List.of(x, lam, sig), List.of(H));
    assertEquals(0, Hf.get_function().size());
    assertEquals(2, H.nnz());
  }
}

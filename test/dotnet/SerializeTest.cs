using System;
using System.Collections;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Runtime.InteropServices;

namespace CasADi.Tests {

/// <summary>Port of test/java SerializeTest (itself a port of test/python/serialize.py)</summary>
public class SerializeTest : CasadiTestCase {

  static string PythonTestDir =>
    Environment.GetEnvironmentVariable("CASADI_TEST_PYTHON") ?? "/mnt/c/Users/Gebruiker/casadi/test/python";

  public void test_compat() {
    string dir = Path.Combine(PythonTestDir, "serialize_3.5.5");
    Check.Assume(Directory.Exists(dir), "no serialize_3.5.5 corpus");
    var inputs = new SortedDictionary<string, List<string>>();
    var outputs = new SortedDictionary<string, List<string>>();
    var functions = new List<string>();
    foreach (string path in Directory.GetFiles(dir)) {
      string n = Path.GetFileName(path);
      if (n.EndsWith(".casadi")) functions.Add(n);
      var m = n.Contains("_in") ? inputs : outputs;
      string k = n.Split('_')[0];
      if (!m.ContainsKey(k)) m[k] = new List<string>();
      m[k].Add(n);
    }
    foreach (string fun in functions) {
      Function f;
      try {
        f = Function.load(Path.Combine(dir, fun));
      } catch (SwigException e) {
        string msg = e.Message ?? "";
        if (msg.Contains("DllLibrary::init_handle")) continue;
        if (msg.Contains("DeSerialization of Integrator failed")) continue;
        if (msg.Contains("DeSerialization of Rootfinder failed")) continue;
        bool skip = false;
        foreach (string token in msg.Split((char[])null, StringSplitOptions.RemoveEmptyEntries)) {
          if (token.StartsWith("libcasadi_")) {
            string[] parts = token.Trim('\'', '"').Split('_', 3);
            if (parts.Length == 3 && !hasPlugin(parts[1], parts[2])) {
              skip = true;
              break;
            }
          }
        }
        if (skip) continue;
        throw;
      }
      string key = fun.Split('.')[0];
      var inl = new List<string>(inputs.GetValueOrDefault(key) ?? new List<string>());
      var outl = new List<string>(outputs.GetValueOrDefault(key) ?? new List<string>());
      inl.Sort(StringComparer.Ordinal);
      outl.Sort(StringComparer.Ordinal);
      for (int k = 0; k < Math.Min(inl.Count, outl.Count); k++) {
        List<DM> inp = f.generate_in(Path.Combine(dir, inl[k]));
        List<DM> outp_ref = f.generate_out(Path.Combine(dir, outl[k]));
        List<DM> outp = f.call(inp);
        assertEquals(outp_ref.Count, outp.Count);
        for (int i = 0; i < outp.Count; i++) {
          double tol = 1e-15;
          string name = f.name_out(i);
          if (f.ToString().Contains("SuperscsInterface")) {
            if (name.StartsWith("lam_")) continue;
            tol = 1e-7;
          }
          if (!RuntimeInformation.IsOSPlatform(OSPlatform.Linux)) tol = 1e-7;
          checkarray(outp[i], outp_ref[i], tol, fun + "/" + f);
        }
      }
    }
  }

  static bool hasPlugin(string kind, string plugin) {
    switch (kind) {
      case "nlpsol": return hasNlpsol(plugin);
      case "conic": return hasConic(plugin);
      case "integrator": return hasIntegrator(plugin);
      case "rootfinder": return hasRootfinder(plugin);
      case "linsol": return hasLinsol(plugin);
      case "interpolant": return ca.has_interpolant(plugin);
      default: return true;
    }
  }

  public void test_identity() {
    var obj = new List<object> {
      "foo",
      new Dictionary<string, object> { ["foo"] = "bar" },
      new List<object> { new Dictionary<string, object> { ["a"] = 3L }, new Dictionary<string, object> { ["b"] = 9L } },
      new List<object> { "a", 5L },
      new Dictionary<string, object> { ["foo"] = new List<object> { "a", 5L } },
      new Dictionary<string, object> { ["foo"] = new List<object> { new List<object> { "a", 5L }, new List<object> { "b", 2L } } },
      new List<object> { new List<object> { "a", 5L }, new List<object> { "b", 2L } },
      new List<object> { new List<object> { new List<object> { "a", 5L }, new List<object> { "b", 2L } } },
    };
    foreach (object e in obj) {
      var ss = ca.StringSerializer();
      if (e is string str) ss.pack(str);
      else ss.pack(e);
      var ds = ca.StringDeserializer(ss.encode());
      object r = ds.unpack();
      Console.WriteLine(e + " " + r);
      checkEqual(e, r);
    }
  }

  static void checkEqual(object a, object b) {
    if (a is IDictionary x) {
      var y = (IDictionary)b;
      assertEquals(string.Join(",", x.Keys.Cast<string>().OrderBy(k => k, StringComparer.Ordinal)),
                   string.Join(",", y.Keys.Cast<string>().OrderBy(k => k, StringComparer.Ordinal)));
      foreach (object k in x.Keys) checkEqual(x[k], y[k]);
    } else if (a is IList xl) {
      var yl = (IList)b;
      assertEquals(xl.Count, yl.Count);
      for (int i = 0; i < xl.Count; i++) checkEqual(xl[i], yl[i]);
    } else {
      assertEquals(a, b);
    }
  }

  // test_pickling_context: python-only (pickle protocol, global_pickle_context)
}

}

// Minimal test runner: every public void method of a class named *Test, without arguments, is a test.
//   Tests [substring...]   runs the tests whose Class.Method contains one of the substrings
using System;
using System.Collections.Generic;
using System.Linq;
using System.Reflection;

namespace CasADi.Tests {

public class AssertFailed : Exception {
  public AssertFailed(string msg) : base(msg) {}
}

/// <summary>Thrown by Check.Skip: the test does not apply here (missing plugin, slow test, ...)</summary>
public class SkipTest : Exception {
  public SkipTest(string msg) : base(msg) {}
}

public static class Check {
  public static void True(bool c, string msg = "expected true") { if (!c) throw new AssertFailed(msg); }
  public static void Equal<T>(T expected, T actual, string what = "") {
    if (!EqualityComparer<T>.Default.Equals(expected, actual))
      throw new AssertFailed($"{what} expected {expected}, got {actual}");
  }
  public static void Near(double expected, double actual, double tol = 1e-9, string what = "") {
    if (!(Math.Abs(expected - actual) <= tol * Math.Max(1, Math.Abs(expected))))
      throw new AssertFailed($"{what} expected {expected}, got {actual}");
  }
  public static void Near(double[,] expected, DM actual, double tol = 1e-9) {
    double[,] a = actual.ToArray();
    Equal(expected.GetLength(0), a.GetLength(0), "rows");
    Equal(expected.GetLength(1), a.GetLength(1), "columns");
    for (int i = 0; i < a.GetLength(0); ++i)
      for (int j = 0; j < a.GetLength(1); ++j) Near(expected[i, j], a[i, j], tol, $"[{i},{j}]");
  }
  public static void Near(double[] expected, DM actual, double tol = 1e-9) {
    double[] a = ca.densify(actual).nonzeros().ToArray();
    Equal(expected.Length, a.Length, "length");
    for (int i = 0; i < a.Length; ++i) Near(expected[i], a[i], tol, $"[{i}]");
  }
  /// <summary>Skips the test unless cond holds (JUnit's assumeTrue)</summary>
  public static void Assume(bool cond, string why) { if (!cond) throw new SkipTest(why); }
  public static void Skip(string why) { throw new SkipTest(why); }
  public static bool RunSlow => Environment.GetEnvironmentVariable("CASADI_DOTNET_RUN_SLOW") != null;

  public static E Throws<E>(Action a) where E : Exception {
    try {
      a();
    } catch (E e) {
      return e;
    }
    throw new AssertFailed("expected " + typeof(E).Name);
  }
}

public static class Runner {
  public static int Main(string[] args) {
    var tests = Assembly.GetExecutingAssembly().GetTypes()
      .Where(t => t.Name.EndsWith("Test") && t.IsClass && !t.IsAbstract)
      .OrderBy(t => t.Name)
      .SelectMany(t => t.GetMethods(BindingFlags.Public | BindingFlags.Instance | BindingFlags.DeclaredOnly)
        .Where(m => m.ReturnType == typeof(void) && m.GetParameters().Length == 0)
        .Select(m => (t, m)))
      .Where(x => args.Length == 0 || args.Any(a => (x.t.Name + "." + x.m.Name).Contains(a)))
      .ToList();
    int failed = 0, skipped = 0;
    // CASADI_DOTNET_TEST_TRACE=1: name each test before it runs, so a crash in native code names its test
    bool trace = Environment.GetEnvironmentVariable("CASADI_DOTNET_TEST_TRACE") == "1";
    foreach (var (t, m) in tests) {
      string name = t.Name + "." + m.Name;
      if (trace) {
        Console.Out.Flush();
        Console.Error.WriteLine("run  " + name);
        Console.Error.Flush();
      }
      try {
        m.Invoke(Activator.CreateInstance(t), null);
        Console.WriteLine("ok   " + name);
      } catch (TargetInvocationException e) {
        Exception inner = e.InnerException ?? e;
        if (inner is SkipTest) {
          skipped++;
          Console.WriteLine("skip " + name + ": " + inner.Message);
          continue;
        }
        failed++;
        Console.WriteLine("FAIL " + name + ": " + inner.GetType().Name + ": " + inner.Message);
        if (!(inner is AssertFailed)) Console.WriteLine(inner.StackTrace);
      }
    }
    // Finalizer-queued deletes must run cleanly too
    GC.Collect();
    GC.WaitForPendingFinalizers();
    SwigRT.RunQueuedDeletes();
    Console.WriteLine($"{tests.Count - failed - skipped} passed, {skipped} skipped, {failed} failed of {tests.Count}");
    return failed == 0 ? 0 : 1;
  }
}

}

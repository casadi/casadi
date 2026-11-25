// Runs one example by name: ExampleRunner rosenbrock
// (each example is a class with its own Main; this entry point lets run_all.sh build them all once)
using System;
using System.Linq;
using System.Reflection;

public static class ExampleRunner {
  public static int Main(string[] args) {
    if (args.Length == 0) {
      Console.WriteLine("usage: Examples <example name>");
      return 2;
    }
    Type t = Assembly.GetExecutingAssembly().GetTypes().FirstOrDefault(x => x.Name == args[0] && x != typeof(ExampleRunner));
    MethodInfo main = t?.GetMethod("Main", BindingFlags.Public | BindingFlags.NonPublic | BindingFlags.Static);
    if (main == null) {
      Console.WriteLine("no example " + args[0]);
      return 2;
    }
    try {
      main.Invoke(null, main.GetParameters().Length == 0 ? null : new object[] { args.Skip(1).ToArray() });
      return 0;
    } catch (TargetInvocationException e) {
      Console.WriteLine((e.InnerException ?? e).ToString());
      return 1;
    }
  }
}

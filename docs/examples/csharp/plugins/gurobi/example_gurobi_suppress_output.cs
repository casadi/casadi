// This simple MIP comes from https://docs.gurobi.com/projects/examples/en/current/examples/python/mip1.html#subsubsectionmip1-py
// Here it is re-written to match the Casadi syntax

using System;
using System.Collections.Generic;
using CasADi;

public static class example_gurobi_suppress_output {
  public static void Main() {
    if (Environment.GetEnvironmentVariable("SKIP_GUROBI_TESTS") != null || !ca.has_conic("gurobi")) {
      Console.WriteLine("SKIP: conic gurobi not available");
      return;
    }

    // Binary decision variables
    SX x = SX.sym("x");
    SX y = SX.sym("y");
    SX z = SX.sym("z");

    // Constraints
    var g = new List<SX>();
    g.Add(x + 2 * y + 3 * z);
    var lbg = new List<double> { -ca.inf };
    var ubg = new List<double> { 4.0 };
    g.Add(x + y);
    lbg.Add(1.0);
    ubg.Add(ca.inf);

    // Objective: minimize -(x + y + 2 z)
    SX f = -(x + y + 2 * z);

    var qp = new Dictionary<string, SX> { ["f"] = f, ["g"] = ca.vertcat(g), ["x"] = ca.vertcat(x, y, z) };
    var bounds = new Dictionary<string, DM> { ["lbx"] = 0, ["ubx"] = 1, ["lbg"] = lbg.ToArray(), ["ubg"] = ubg.ToArray() };

    // Standard behavior of gurobi -- printing to console is allowed
    // The option suppress_all_output prevents gurobi to print to console.
    Console.WriteLine("\n\n#################### STANDARD GUROBI BEHAVIOR ####################");
    Function solver = ca.qpsol("solver", "gurobi", qp,
        new Dictionary<string, object> { ["discrete"] = new[] { true, true, true } });

    // Solve
    var sol = solver.call(bounds);
    Console.WriteLine($"Optimal solution: [{string.Join(", ", sol["x"].nonzeros())}]");

    sol = solver.call(bounds);
    Console.WriteLine($"Optimal solution: [{string.Join(", ", sol["x"].nonzeros())}]");

    // Even if you set gurobi options OutputFlag or LogToConsole some information like the license banner will still be printed, making the log cluttered in case of repetitive optimization like in model predictive control
    Console.WriteLine("\n\n#################### GUROBI WITH OutputFlag=0 ####################");
    solver = ca.qpsol("solver", "gurobi", qp,
        new Dictionary<string, object> { ["discrete"] = new[] { true, true, true },
                                         ["gurobi.OutputFlag"] = 0L });

    // Solve
    sol = solver.call(bounds);
    Console.WriteLine($"Optimal solution: [{string.Join(", ", sol["x"].nonzeros())}]");

    sol = solver.call(bounds);
    Console.WriteLine($"Optimal solution: [{string.Join(", ", sol["x"].nonzeros())}]");

    // Setting both OutputFlag=0 and LogToConsole=0 prevents gurobi to print to console.
    // It follows the guidelines of setting Gurobi OutputFlag before creating the gurobi environment like reported here:
    // https://support.gurobi.com/hc/en-us/articles/360044784552-How-do-I-suppress-all-console-output-from-Gurobi
    Console.WriteLine("\n\n#################### GUROBI WITH both OutputFlag=0 and LogToConsole=0 ####################");
    solver = ca.qpsol("solver", "gurobi", qp,
        new Dictionary<string, object> { ["discrete"] = new[] { true, true, true },
                                         ["gurobi.OutputFlag"] = 0L,
                                         ["gurobi.LogToConsole"] = 0L });

    // Solve
    sol = solver.call(bounds);
    Console.WriteLine($"Optimal solution: [{string.Join(", ", sol["x"].nonzeros())}]");

    sol = solver.call(bounds);
    Console.WriteLine($"Optimal solution: [{string.Join(", ", sol["x"].nonzeros())}]");
  }
}

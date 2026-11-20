// This simple MIP comes from https://docs.gurobi.com/projects/examples/en/current/examples/python/mip1.html#subsubsectionmip1-py
// Here it is re-written to match the Casadi syntax

import org.casadi.*;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Map;

public class example_gurobi_suppress_output {
  public static void main(String[] args) {
    if (System.getenv("SKIP_GUROBI_TESTS") != null || !ca.has_conic("gurobi")) {
      System.out.println("SKIP: conic gurobi not available");
      return;
    }
    double inf = Double.POSITIVE_INFINITY;

    // Binary decision variables
    SX x = SX.sym("x");
    SX y = SX.sym("y");
    SX z = SX.sym("z");

    // Constraints
    List<SX> g = new ArrayList<>();
    g.add(x.plus(y.times(2)).plus(z.times(3)));
    List<Double> lbg = new ArrayList<>(List.of(-inf));
    List<Double> ubg = new ArrayList<>(List.of(4.0));
    g.add(x.plus(y));
    lbg.add(1.0);
    ubg.add(inf);

    // Objective: minimize -(x + y + 2 z)
    SX f = x.plus(y).plus(z.times(2)).unaryMinus();

    // Standard behavior of gurobi -- printing to console is allowed
    // The option suppress_all_output prevents gurobi to print to console.
    System.out.println("\n\n#################### STANDARD GUROBI BEHAVIOR ####################");
    Function solver = ca.qpsol("solver", "gurobi",
        Map.of("f", f, "g", ca.vertcat(g), "x", ca.vertcat(List.of(x, y, z))),
        Map.of("discrete", new boolean[] {true, true, true}));

    // Solve
    Map<String, GenericMatrix> sol = solver.invoke(Map.of("lbx", 0, "ubx", 1, "lbg", arr(lbg), "ubg", arr(ubg)));
    System.out.println("Optimal solution: " + Arrays.toString(((DM) sol.get("x")).nonzeros()));

    sol = solver.invoke(Map.of("lbx", 0, "ubx", 1, "lbg", arr(lbg), "ubg", arr(ubg)));
    System.out.println("Optimal solution: " + Arrays.toString(((DM) sol.get("x")).nonzeros()));

    // Even if you set gurobi options OutputFlag or LogToConsole some information like the license banner will still be printed, making the log cluttered in case of repetitive optimization like in model predictive control
    System.out.println("\n\n#################### GUROBI WITH OutputFlag=0 ####################");
    solver = ca.qpsol("solver", "gurobi",
        Map.of("f", f, "g", ca.vertcat(g), "x", ca.vertcat(List.of(x, y, z))),
        Map.of("discrete", new boolean[] {true, true, true},
               "gurobi.OutputFlag", 0));

    // Solve
    sol = solver.invoke(Map.of("lbx", 0, "ubx", 1, "lbg", arr(lbg), "ubg", arr(ubg)));
    System.out.println("Optimal solution: " + Arrays.toString(((DM) sol.get("x")).nonzeros()));

    sol = solver.invoke(Map.of("lbx", 0, "ubx", 1, "lbg", arr(lbg), "ubg", arr(ubg)));
    System.out.println("Optimal solution: " + Arrays.toString(((DM) sol.get("x")).nonzeros()));

    // Setting both OutputFlag=0 and LogToConsole=0 prevents gurobi to print to console.
    // It follows the guidelines of setting Gurobi OutputFlag before creating the gurobi environment like reported here:
    // https://support.gurobi.com/hc/en-us/articles/360044784552-How-do-I-suppress-all-console-output-from-Gurobi
    System.out.println("\n\n#################### GUROBI WITH both OutputFlag=0 and LogToConsole=0 ####################");
    solver = ca.qpsol("solver", "gurobi",
        Map.of("f", f, "g", ca.vertcat(g), "x", ca.vertcat(List.of(x, y, z))),
        Map.of("discrete", new boolean[] {true, true, true},
               "gurobi.OutputFlag", 0,
               "gurobi.LogToConsole", 0));

    // Solve
    sol = solver.invoke(Map.of("lbx", 0, "ubx", 1, "lbg", arr(lbg), "ubg", arr(ubg)));
    System.out.println("Optimal solution: " + Arrays.toString(((DM) sol.get("x")).nonzeros()));

    sol = solver.invoke(Map.of("lbx", 0, "ubx", 1, "lbg", arr(lbg), "ubg", arr(ubg)));
    System.out.println("Optimal solution: " + Arrays.toString(((DM) sol.get("x")).nonzeros()));
  }

  static double[] arr(List<Double> v) {
    return v.stream().mapToDouble(Double::doubleValue).toArray();
  }
}

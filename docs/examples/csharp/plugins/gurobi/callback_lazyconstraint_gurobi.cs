// This example shows how to implement a Gurobi callback to add LazyConstraint at MIPSOL using the Casadi/Gurobi interface
// The results obtained by this program are identical to the one in the official Gurobi documentation available at:
// https://colab.research.google.com/github/Gurobi/modeling-examples/blob/master/traveling_salesman/tsp.ipynb
// At the bottom of the Python version of this file there is an implementation of the Colab notebook

using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Net.Http;
using System.Text.Json;
using CasADi;

public static class callback_lazyconstraint_gurobi {
  // A pair of cities (Python uses tuples as dictionary keys)
  record Pair(string i, string j);

  static List<string> capitals = new List<string>();
  static Dictionary<string, double[]> coordinates = new Dictionary<string, double[]>();

  public static void Main() {
    if (Environment.GetEnvironmentVariable("SKIP_GUROBI_TESTS") != null || !ca.has_conic("gurobi")) {
      Console.WriteLine("SKIP: conic gurobi not available");
      return;
    }

    // Read capital data (same as your original example)
    string capitals_json;
    try {
      capitals_json = File.ReadAllText("../../../python/plugins/gurobi/capitals.json");
    } catch (IOException) {
      string url = "https://raw.githubusercontent.com/Gurobi/modeling-examples/master/traveling_salesman/capitals.json";
      using var http = new HttpClient();
      capitals_json = http.GetStringAsync(url).Result;
    }

    // A JSON value that may be a number or a string holding a number
    static double number(JsonElement v) => v.ValueKind == JsonValueKind.String ? double.Parse(v.GetString(), System.Globalization.CultureInfo.InvariantCulture) : v.GetDouble();
    using (JsonDocument doc = JsonDocument.Parse(capitals_json)) {
      foreach (JsonProperty e in doc.RootElement.EnumerateObject()) {
        string state = e.Name;
        if (!new[] { "AK", "HI" }.Contains(state)) {
          string capital = e.Value.GetProperty("capital").GetString();
          capitals.Add(capital);
          coordinates[capital] = new[] { number(e.Value.GetProperty("lat")), number(e.Value.GetProperty("long")) };
        }
      }
    }

    // Compute pairwise distance matrix
    var dist = new Dictionary<Pair, double>();
    var pairs = new List<Pair>();  // the keys of dist, in insertion order
    for (int a = 0; a < capitals.Count; a++)
      for (int b = a + 1; b < capitals.Count; b++) {
        var p = new Pair(capitals[a], capitals[b]);
        pairs.Add(p);
        dist[p] = distance(capitals[a], capitals[b]);
      }

    // Create the callback
    var callback = new SubtourEliminationCallback("subtour_callback", capitals, pairs);

    // Create CasADi optimization problem

    // Variables: is city 'i' adjacent to city 'j' on the tour?
    // Create binary variables for all pairs
    var vars = new Dictionary<Pair, SX>();
    var x = new List<SX>();
    foreach (Pair p in pairs) {
      x.Add(SX.sym($"x_{p.i}_{p.j}"));
      vars[p] = x[^1];
      vars[new Pair(p.j, p.i)] = vars[p];  // Make symmetric variables point to same symbol
    }

    // Add objective
    SX objective = 0;
    foreach (Pair p in pairs) {
      objective += vars[p] * dist[p];
    }

    // Constraints: two edges incident to each city
    var g = new List<SX>();
    var lbg = new List<double>();
    var ubg = new List<double>();
    foreach (string c in capitals) {
      SX constraint_sum = 0;
      foreach (Pair p in pairs) {
        if (p.i == c || p.j == c) {
          constraint_sum += vars[p];
        }
      }
      g.Add(constraint_sum - 2);  // constraint_sum == 2, so constraint_sum - 2 == 0
      lbg.Add(0.0);
      ubg.Add(0.0);
    }

    // Set up the solver with callback
    var discrete = Enumerable.Repeat(true, x.Count).ToList();
    Function solver = ca.qpsol(
        "myqp",
        "gurobi",
        new Dictionary<string, SX> { ["f"] = objective, ["x"] = ca.vertcat(x), ["g"] = ca.vertcat(g) },
        new Dictionary<string, object> {
          ["discrete"] = discrete,
          ["lazy_constraints_callback"] = callback });

    // Solve
    Dictionary<string, DM> sol;
    try {
      sol = solver.call(new Dictionary<string, DM> { ["lbx"] = 0, ["ubx"] = 1, ["lbg"] = lbg.ToArray(), ["ubg"] = ubg.ToArray() });
      Console.WriteLine("Optimization successful!");
      Console.WriteLine($"Objective value: {(double)sol["f"]}");
    } catch (SwigException e) {
      Console.WriteLine($"Optimization failed: {e.Message}");
      throw;
    }

    List<double> x_sol = sol["x"].nonzeros();

    // Find which edges are selected in the current solution
    var selected_edges_1 = new List<Pair>();
    var selected_edges_2 = new List<Pair>();
    for (int idx = 0; idx < pairs.Count; idx++) {
      Pair p = pairs[idx];
      if (x_sol[idx] > 0.5) {
        selected_edges_1.Add(p);
        selected_edges_2.Add(new Pair(p.j, p.i));
      }
    }
    List<Pair> selected_edges = selected_edges_1.Concat(selected_edges_2).ToList();
    List<string> tour = callback.find_subtour(selected_edges);

    if (tour.Count != capitals.Count) throw new Exception("assertion failed");

    // Map the solution: python-only (folium)

    // Solving using Gurobipy: python-only (Gurobi's own Python API), see the Python version of this example
  }

  static double distance(string city1, string city2) {
    double[] c1 = coordinates[city1];
    double[] c2 = coordinates[city2];
    double[] diff = { c1[0] - c2[0], c1[1] - c2[1] };
    return Math.Sqrt(diff[0] * diff[0] + diff[1] * diff[1]);
  }

  class SubtourEliminationCallback : Callback {
    readonly List<string> capitals;
    readonly List<Pair> var_names;
    readonly int n_vars;

    public SubtourEliminationCallback(string name, List<string> capitals_list, List<Pair> distances, IDictionary<string, object> opts) {
      capitals = capitals_list;
      var_names = new List<Pair>(distances);
      n_vars = distances.Count;
      construct(name, opts);
    }

    public SubtourEliminationCallback(string name, List<string> capitals_list, List<Pair> distances)
      : this(name, capitals_list, distances, new Dictionary<string, object>()) { }

    public override long get_n_in() {
      return 5;  // x_solution, obj_val, obj_best, obj_bound, sol_count
    }

    public override long get_n_out() {
      return 3;  // add_lazy_flag, A_lazy, b_lazy
    }

    public override Sparsity get_sparsity_in(long i) {
      if (i == 0) {  // x_solution - will be a dense vector of all variables
        return Sparsity.dense(n_vars, 1);
      } else if (i == 1) {  // obj_val
        return Sparsity.dense(1, 1);
      } else if (i == 2) {  // obj_best
        return Sparsity.dense(1, 1);
      } else if (i == 3) {  // obj_bound
        return Sparsity.dense(1, 1);
      } else if (i == 4) {  // sol_count
        return Sparsity.dense(1, 1);
      } else {
        return ca.Sparsity(0, 0);
      }
    }

    public override Sparsity get_sparsity_out(long i) {
      if (i == 0) {  // add_lazy_flag
        return Sparsity.dense(1, 1);
      } else if (i == 1) {  // A_lazy - sparse matrix for constraint coefficients
        return Sparsity.dense(1, n_vars);
      } else if (i == 2) {  // b_lazy
        return Sparsity.dense(1, 1);
      } else {
        return ca.Sparsity(0, 0);
      }
    }

    public override List<DM> eval(IEnumerable<DM> arg) {
      var args = arg.ToList();
      // Extract solution data
      List<double> x_solution = args[0].nonzeros();  // Vector of all variable values
      double obj_val = (double)args[1];
      double obj_best = (double)args[2];
      double obj_bound = (double)args[3];
      int sol_count = (int)(double)args[4];

      // Find which edges are selected in the current solution
      var selected_edges_1 = new List<Pair>();
      var selected_edges_2 = new List<Pair>();
      for (int idx = 0; idx < var_names.Count; idx++) {
        Pair p = var_names[idx];
        if (x_solution[idx] > 0.5) {
          selected_edges_1.Add(p);
          selected_edges_2.Add(new Pair(p.j, p.i));  // need to handle the symmetry in the variable vector
        }
      }
      List<Pair> selected_edges = selected_edges_1.Concat(selected_edges_2).ToList();
      // Find subtour
      List<string> tour = find_subtour(selected_edges);
      var cities_in_tour = new HashSet<Pair>();
      for (int a = 0; a < tour.Count; a++)
        for (int b = a + 1; b < tour.Count; b++) cities_in_tour.Add(new Pair(tour[a], tour[b]));
      var result = new double[n_vars];
      var result_symm = new double[n_vars];  // need to handle the symmetry in the variable vector
      for (int k = 0; k < n_vars; k++) {
        Pair p = var_names[k];
        result[k] = cities_in_tour.Contains(p) ? 1 : 0;
        result_symm[k] = cities_in_tour.Contains(new Pair(p.j, p.i)) ? 1 : 0;
      }
      DM resultDM = (DM)result + (DM)result_symm;
      // Check if we found a valid subtour that needs constraint
      if (tour.Count < capitals.Count) {
        // Create lazy constraint to eliminate this subtour
        // Constraint: sum of edges in subtour <= len(tour) - 1
        double add_lazy_flag = 1.0;
        // Create coefficient matrix A and right-hand side b
        DM A_lazy = resultDM.T();
        DM b_lazy = new double[] { tour.Count - 1 };
        return new List<DM> { add_lazy_flag, A_lazy, b_lazy };
      } else {
        // No constraint to add
        double add_lazy_flag = 0.0;
        DM A_lazy = ca.DM(1, n_vars);
        DM b_lazy = ca.DM(1, 1);
        return new List<DM> { add_lazy_flag, A_lazy, b_lazy };
      }
    }

    /// <summary>Find shortest subtour in the selected edges</summary>
    public List<string> find_subtour(List<Pair> edges) {
      var unvisited = new List<string>(capitals);
      var cycle = new List<string>(capitals);  // Dummy

      while (unvisited.Count > 0) {
        var thiscycle = new List<string>();
        var neighbors = new List<string>(unvisited);

        while (neighbors.Count > 0) {
          string current = neighbors[0];
          thiscycle.Add(current);
          unvisited.Remove(current);
          neighbors = edges.Where(e => e.i == current && unvisited.Contains(e.j)).Select(e => e.j).ToList();
        }

        if (thiscycle.Count <= cycle.Count) {
          cycle = thiscycle;
        }
      }

      return cycle;
    }
  }
}

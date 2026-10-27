// This example shows how to implement a Gurobi callback to add LazyConstraint at MIPSOL using the Casadi/Gurobi interface
// The results obtained by this program are identical to the one in the official Gurobi documentation available at:
// https://colab.research.google.com/github/Gurobi/modeling-examples/blob/master/traveling_salesman/tsp.ipynb
// At the bottom of the Python version of this file there is an implementation of the Colab notebook

import org.casadi.*;
import java.net.URI;
import java.net.http.HttpClient;
import java.net.http.HttpRequest;
import java.net.http.HttpResponse;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class callback_lazyconstraint_gurobi {
  // A pair of cities (Python uses tuples as dictionary keys)
  record Pair(String i, String j) {}

  static List<String> capitals = new ArrayList<>();
  static Map<String, double[]> coordinates = new LinkedHashMap<>();

  public static void main(String[] args) throws Exception {
    if (System.getenv("SKIP_GUROBI_TESTS") != null || !ca.has_conic("gurobi")) {
      System.out.println("SKIP: conic gurobi not available");
      return;
    }

    // Read capital data (same as your original example)
    String capitals_json;
    try {
      capitals_json = Files.readString(Paths.get("../../../python/plugins/gurobi/capitals.json"));
    } catch (java.io.IOException e) {
      String url = "https://raw.githubusercontent.com/Gurobi/modeling-examples/master/traveling_salesman/capitals.json";
      capitals_json = HttpClient.newHttpClient().send(HttpRequest.newBuilder(URI.create(url)).build(),
                                                      HttpResponse.BodyHandlers.ofString()).body();
    }

    for (Map.Entry<String, Map<String, String>> e : parseCapitals(capitals_json).entrySet()) {
      String state = e.getKey();
      if (!List.of("AK", "HI").contains(state)) {
        String capital = e.getValue().get("capital");
        capitals.add(capital);
        coordinates.put(capital, new double[] {Double.parseDouble(e.getValue().get("lat")),
                                               Double.parseDouble(e.getValue().get("long"))});
      }
    }

    // Compute pairwise distance matrix
    Map<Pair, Double> dist = new LinkedHashMap<>();
    for (int a = 0; a < capitals.size(); a++)
      for (int b = a + 1; b < capitals.size(); b++)
        dist.put(new Pair(capitals.get(a), capitals.get(b)), distance(capitals.get(a), capitals.get(b)));

    // Create the callback
    SubtourEliminationCallback callback = new SubtourEliminationCallback("subtour_callback", capitals, dist);

    // Create CasADi optimization problem

    // Variables: is city 'i' adjacent to city 'j' on the tour?
    // Create binary variables for all pairs
    Map<Pair, SX> vars = new LinkedHashMap<>();
    List<SX> x = new ArrayList<>();
    for (Pair p : dist.keySet()) {
      x.add(SX.sym("x_" + p.i() + "_" + p.j()));
      vars.put(p, x.get(x.size() - 1));
      vars.put(new Pair(p.j(), p.i()), vars.get(p));  // Make symmetric variables point to same symbol
    }

    // Add objective
    SX objective = new SX(0);
    for (Pair p : dist.keySet()) {
      objective = objective.plus(vars.get(p).times(dist.get(p)));
    }

    // Constraints: two edges incident to each city
    List<SX> g = new ArrayList<>();
    List<Double> lbg = new ArrayList<>();
    List<Double> ubg = new ArrayList<>();
    for (String c : capitals) {
      SX constraint_sum = new SX(0);
      for (Pair p : dist.keySet()) {
        if (p.i().equals(c) || p.j().equals(c)) {
          constraint_sum = constraint_sum.plus(vars.get(p));
        }
      }
      g.add(constraint_sum.minus(2));  // constraint_sum == 2, so constraint_sum - 2 == 0
      lbg.add(0.0);
      ubg.add(0.0);
    }

    // Set up the solver with callback
    boolean[] discrete = new boolean[x.size()];
    java.util.Arrays.fill(discrete, true);
    Function solver = ca.qpsol(
        "myqp",
        "gurobi",
        Map.of("f", objective, "x", ca.vertcat(x), "g", ca.vertcat(g)),
        Map.of(
          "discrete", discrete,
          "lazy_constraints_callback", callback));

    // Solve
    Map<String, DM> sol;
    try {
      sol = solver.call(Map.of("lbx", new DM(0), "ubx", new DM(1), "lbg", new DM(arr(lbg)), "ubg", new DM(arr(ubg))));
      System.out.println("Optimization successful!");
      System.out.println("Objective value: " + sol.get("f").to_double());
    } catch (CasadiException e) {
      System.out.println("Optimization failed: " + e.getMessage());
      throw e;
    }

    double[] x_sol = sol.get("x").nonzeros();

    // Find which edges are selected in the current solution
    List<Pair> selected_edges_1 = new ArrayList<>();
    List<Pair> selected_edges_2 = new ArrayList<>();
    int idx = 0;
    for (Pair p : dist.keySet()) {
      if (x_sol[idx++] > 0.5) {
        selected_edges_1.add(p);
        selected_edges_2.add(new Pair(p.j(), p.i()));
      }
    }
    List<Pair> selected_edges = new ArrayList<>(selected_edges_1);
    selected_edges.addAll(selected_edges_2);
    List<String> tour = callback.find_subtour(selected_edges);

    if (tour.size() != capitals.size()) throw new AssertionError();

    // Map the solution: python-only (folium)

    // Solving using Gurobipy: python-only (Gurobi's own Python API), see the Python version of this example
  }

  static double distance(String city1, String city2) {
    double[] c1 = coordinates.get(city1);
    double[] c2 = coordinates.get(city2);
    double[] diff = {c1[0] - c2[0], c1[1] - c2[1]};
    return Math.sqrt(diff[0] * diff[0] + diff[1] * diff[1]);
  }

  static class SubtourEliminationCallback extends Callback {
    List<String> capitals;
    List<Pair> var_names;
    int n_vars;

    SubtourEliminationCallback(String name, List<String> capitals_list, Map<Pair, Double> distances, Map<String, ?> opts) {
      this.capitals = capitals_list;
      this.var_names = new ArrayList<>(distances.keySet());
      this.n_vars = distances.size();
      construct(name, opts);
    }

    SubtourEliminationCallback(String name, List<String> capitals_list, Map<Pair, Double> distances) {
      this(name, capitals_list, distances, Map.of());
    }

    @Override
    public long get_n_in() {
      return 5;  // x_solution, obj_val, obj_best, obj_bound, sol_count
    }

    @Override
    public long get_n_out() {
      return 3;  // add_lazy_flag, A_lazy, b_lazy
    }

    @Override
    public Sparsity get_sparsity_in(long i) {
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
        return new Sparsity(0, 0);
      }
    }

    @Override
    public Sparsity get_sparsity_out(long i) {
      if (i == 0) {  // add_lazy_flag
        return Sparsity.dense(1, 1);
      } else if (i == 1) {  // A_lazy - sparse matrix for constraint coefficients
        return Sparsity.dense(1, n_vars);
      } else if (i == 2) {  // b_lazy
        return Sparsity.dense(1, 1);
      } else {
        return new Sparsity(0, 0);
      }
    }

    @Override
    public List<DM> eval(List<DM> arg) {
      // Extract solution data
      double[] x_solution = arg.get(0).nonzeros();  // Vector of all variable values
      double obj_val = arg.get(1).to_double();
      double obj_best = arg.get(2).to_double();
      double obj_bound = arg.get(3).to_double();
      int sol_count = (int) arg.get(4).to_double();

      // Find which edges are selected in the current solution
      List<Pair> selected_edges_1 = new ArrayList<>();
      List<Pair> selected_edges_2 = new ArrayList<>();
      for (int idx = 0; idx < var_names.size(); idx++) {
        Pair p = var_names.get(idx);
        if (x_solution[idx] > 0.5) {
          selected_edges_1.add(p);
          selected_edges_2.add(new Pair(p.j(), p.i()));  // need to handle the symmetry in the variable vector
        }
      }
      List<Pair> selected_edges = new ArrayList<>(selected_edges_1);
      selected_edges.addAll(selected_edges_2);
      // Find subtour
      List<String> tour = find_subtour(selected_edges);
      Set<Pair> cities_in_tour = new HashSet<>();
      for (int a = 0; a < tour.size(); a++)
        for (int b = a + 1; b < tour.size(); b++) cities_in_tour.add(new Pair(tour.get(a), tour.get(b)));
      double[] result = new double[n_vars];
      double[] result_symm = new double[n_vars];  // need to handle the symmetry in the variable vector
      for (int k = 0; k < n_vars; k++) {
        Pair p = var_names.get(k);
        result[k] = cities_in_tour.contains(p) ? 1 : 0;
        result_symm[k] = cities_in_tour.contains(new Pair(p.j(), p.i())) ? 1 : 0;
      }
      DM resultDM = new DM(result).plus(new DM(result_symm));
      // Check if we found a valid subtour that needs constraint
      if (tour.size() < capitals.size()) {
        // Create lazy constraint to eliminate this subtour
        // Constraint: sum of edges in subtour <= len(tour) - 1
        double add_lazy_flag = 1.0;
        // Create coefficient matrix A and right-hand side b
        DM A_lazy = resultDM.T();
        DM b_lazy = new DM(new double[] {tour.size() - 1});
        return List.of(new DM(add_lazy_flag), A_lazy, b_lazy);
      } else {
        // No constraint to add
        double add_lazy_flag = 0.0;
        DM A_lazy = new DM(1, n_vars);
        DM b_lazy = new DM(1, 1);
        return List.of(new DM(add_lazy_flag), A_lazy, b_lazy);
      }
    }

    /** Find shortest subtour in the selected edges */
    List<String> find_subtour(List<Pair> edges) {
      List<String> unvisited = new ArrayList<>(capitals);
      List<String> cycle = new ArrayList<>(capitals);  // Dummy

      while (!unvisited.isEmpty()) {
        List<String> thiscycle = new ArrayList<>();
        List<String> neighbors = new ArrayList<>(unvisited);

        while (!neighbors.isEmpty()) {
          String current = neighbors.get(0);
          thiscycle.add(current);
          unvisited.remove(current);
          neighbors = new ArrayList<>();
          for (Pair e : edges) if (e.i().equals(current) && unvisited.contains(e.j())) neighbors.add(e.j());
        }

        if (thiscycle.size() <= cycle.size()) {
          cycle = thiscycle;
        }
      }

      return cycle;
    }
  }

  // Minimal reader for capitals.json: {"AL": {"capital": "...", "lat": "...", "long": "..."}, ...}
  static Map<String, Map<String, String>> parseCapitals(String json) {
    Map<String, Map<String, String>> r = new LinkedHashMap<>();
    Matcher m = Pattern.compile("\"([^\"]+)\"\\s*:\\s*\\{([^}]*)\\}").matcher(json);
    while (m.find()) {
      Map<String, String> fields = new LinkedHashMap<>();
      Matcher f = Pattern.compile("\"([^\"]+)\"\\s*:\\s*\"?([^\",]*)\"?").matcher(m.group(2));
      while (f.find()) fields.put(f.group(1), f.group(2).trim());
      r.put(m.group(1), fields);
    }
    return r;
  }

  static double[] arr(List<Double> v) {
    return v.stream().mapToDouble(Double::doubleValue).toArray();
  }
}

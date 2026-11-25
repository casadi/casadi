using System;
using System.Collections.Generic;
using System.Linq;

namespace CasADi.Tests {

/// <summary>Port of test/java SimulatorTest (test/python/simulator.py)</summary>
public class SimulatorTest : CasadiTestCase {
  Dictionary<string, SX> dae;
  Function integrator, qe;
  MX qend, q0, par;
  const double numTend = 2.3, numQ0 = 7.1, numP = 2;

  static double[] linspace(double a, double b, int n) {
    double[] r = new double[n];
    for (int i = 0; i < n; i++) r[i] = a + (b - a) * i / (n - 1);
    return r;
  }

  static Dictionary<string, object> opts() => new Dictionary<string, object> {
    ["reltol"] = 1e-15, ["abstol"] = 1e-15, ["fsens_err_con"] = true,
  };

  /// <summary>helpers.checkfunction with adj=False, jacobian=False, sens_der=False, evals=False: outputs and forward</summary>
  void checkOutputsForward(Function trial, Function solution, List<DM> inputs, double tol) {
    assertEquals(solution.n_in(), trial.n_in());
    assertEquals(solution.n_out(), trial.n_out());
    List<DM> a = trial.call(inputs), b = solution.call(inputs);
    for (int i = 0; i < a.Count; i++) {
      if (trial.sparsity_out(i).is_empty() || solution.sparsity_out(i).is_empty()) continue;
      checkarray(a[i], b[i], tol, "output " + trial.name_out(i));
    }
    var fwd = inputs.Concat(a).Concat(randomInputs(trial)).ToList();
    List<DM> fa = trial.forward(1).call(fwd), fb = solution.forward(1).call(fwd);
    for (int i = 0; i < fa.Count; i++) {
      if (fa[i].is_empty() || fb[i].is_empty()) continue;
      checkarray(fa[i], fb[i], tol, "forward " + i);
    }
  }

  void setUp() {
    Check.Assume(hasIntegrator("cvodes"), "no integrator cvodes");
    SX t = SX.sym("t"), q = SX.sym("q"), p = SX.sym("p");
    var f = new Dictionary<string, SX> { ["t"] = t, ["x"] = q, ["p"] = p, ["ode"] = q / p * ca.pow(t, 2) };
    integrator = ca.integrator("integrator", "cvodes", f, 0, 2.3, opts());
    q0 = MX.sym("q0");
    par = MX.sym("p");
    qend = integrator.call(new Dictionary<string, MX> { ["x0"] = q0, ["p"] = par })["xf"];
    qe = ca.Function("qe", new[] { q0, par }, new[] { qend });
    dae = f;
  }

  public void test_sim_full() {
    setUp();
    message("Simulator inputs");
    int N = 4;
    double[] tc = linspace(0, numTend, N);
    SX t = SX.sym("t"), q = SX.sym("q"), p = SX.sym("p");
    var f = new Dictionary<string, SX> { ["t"] = t, ["x"] = q, ["p"] = p, ["ode"] = q / p * ca.pow(t, 2) };
    Function integrator = ca.integrator("integrator", "cvodes", f, tc[0], tc, opts());
    var cols = tc.Select(ti => q * ca.exp(Math.Pow(ti, 3) / (p * 3))).ToList();
    var sol = new Dictionary<string, SX> { ["x0"] = q, ["p"] = p, ["xf"] = ca.horzcat(cols) };
    var solution = ca.Function("solution", sol, ca.integrator_in(), ca.integrator_out());
    var fin = new Dictionary<string, DM> { ["x0"] = 0.3, ["p"] = 0.7 };
    var inputs = integrator.name_in()
      .Select(n => fin.TryGetValue(n, out DM v) ? v : DM.zeros(integrator.sparsity_in(n))).ToList();
    checkOutputsForward(integrator, solution, inputs, 1e-6);
  }

  public void test_simulator_time_offset() {
    setUp();
    message("CVodes integration: simulator time offset");
    double[] t = linspace(0.7, numTend, 100);
    Function integrator = ca.integrator("integrator", "cvodes", dae, t[0], t, opts());
    int x0_ind = ca.integrator_in().IndexOf("x0");
    var integrator_in_ = Enumerable.Range(0, (int)integrator.n_in()).Select(_ => ca.DM(0)).ToList();
    integrator_in_[x0_ind] = ca.DM(numQ0);
    int p_ind = ca.integrator_in().IndexOf("p");
    integrator_in_[p_ind] = ca.DM(numP);
    List<DM> integrator_out_ = integrator.call(integrator_in_);
    double tend = numTend, q0 = numQ0, p = numP;
    DM xf = integrator_out_[0];
    double expected = q0 * Math.Exp((Math.Pow(tend, 3) - Math.Pow(0.7, 3)) / (3 * p)), actual = (double)xf[0, xf.size2() - 1];
    Check.True(Math.Abs(expected - actual) <= 0.5e-9, $"Evaluation output mismatch: {expected} vs {actual}");
  }
}

}
